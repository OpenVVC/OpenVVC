/**
 *
 *   OpenVVC is open-source real time software decoder compliant with the 
 *   ITU-T H.266- MPEG-I - Part 3 VVC standard. OpenVVC is developed from 
 *   scratch in C as a library that provides consumers with real time and
 *   energy-aware decoding capabilities under different OS including MAC OS,
 *   Windows, Linux and Android targeting low energy real-time decoding of
 *   4K VVC videos on Intel x86 and ARM platforms.
 * 
 *   Copyright (C) 2020-2022  IETR-INSA Rennes :
 *   
 *   Pierre-Loup CABARAT
 *   Wassim HAMIDOUCHE
 *   Guillaume GAUTIER
 *   Thomas AMESTOY
 *
 *   This library is free software; you can redistribute it and/or
 *   modify it under the terms of the GNU Lesser General Public
 *   License as published by the Free Software Foundation; either
 *   version 2.1 of the License, or (at your option) any later version.
 *
 *   This library is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *   Lesser General Public License for more details.
 *
 *   You should have received a copy of the GNU Lesser General Public
 *   License along with this library; if not, write to the Free Software
 *   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301
 *   USA
 * 
 **/

#include <stdint.h>

#include "ovutils.h"
#include "rcn.h"
#include "ctudec.h"
#include "drv_utils.h"

static const int8_t g_angle2mask[GEO_NUM_ANGLES] =
{
    0, -1,  1,  2,  3,  4, -1, -1,
    5, -1, -1,  4,  3,  2,  1, -1,

    0, -1,  1,  2,  3,  4, -1, -1,
    5, -1, -1,  4,  3,  2,  1, -1
};

int16_t g_globalGeoWeights [GEO_NUM_PRESTORED_MASK][GEO_WEIGHT_MASK_SIZE * GEO_WEIGHT_MASK_SIZE];
int16_t g_weightOffset     [GEO_NUM_PARTITION_MODE][GEO_NUM_CU_SIZE][GEO_NUM_CU_SIZE][2];

const int16_t g_GeoParams[GEO_NUM_PARTITION_MODE][2] =
{
    {0, 1},
    {0, 3},

    {2, 0},
    {2, 1},
    {2, 2},
    {2, 3},

    {3, 0},
    {3, 1},
    {3, 2},
    {3, 3},

    {4, 0},
    {4, 1},
    {4, 2},
    {4, 3},

    {5, 0},
    {5, 1},
    {5, 2},
    {5, 3},

    {8, 1},
    {8, 3},

    {11, 0},
    {11, 1},
    {11, 2},
    {11, 3},

    {12, 0},
    {12, 1},
    {12, 2},
    {12, 3},

    {13, 0},
    {13, 1},
    {13, 2},
    {13, 3},

    {14, 0},
    {14, 1},
    {14, 2},
    {14, 3},

    {16, 1},
    {16, 3},

    {18, 1},
    {18, 2},
    {18, 3},

    {19, 1},
    {19, 2},
    {19, 3},

    {20, 1},
    {20, 2},
    {20, 3},

    {21, 1},
    {21, 2},
    {21, 3},

    {24, 1},
    {24, 3},

    {27, 1},
    {27, 2},
    {27, 3},

    {28, 1},
    {28, 2},
    {28, 3},

    {29, 1},
    {29, 2},
    {29, 3},

    {30, 1},
    {30, 2},
    {30, 3}
};

const int8_t g_Dis[GEO_NUM_ANGLES] =
{
    8,  8,  8,  8,  4,  4,  2,  1,
    0, -1, -2, -4, -4, -8, -8, -8,

   -8, -8, -8, -8, -4, -4, -2, -1,
    0,  1,  2,  4,  4,  8,  8,  8
};

void
rcn_init_gpm_params()
{
  int angle_idx;
  for (angle_idx = 0; angle_idx < (GEO_NUM_ANGLES >> 2) + 1; angle_idx++){
    int x, y;

    if (g_angle2mask[angle_idx] == -1){
      continue;
    }

    int dist_x = angle_idx;
    int dist_y = (dist_x + (GEO_NUM_ANGLES >> 2)) % GEO_NUM_ANGLES;
    int16_t rho = (int32_t)((uint32_t)g_Dis[dist_x] << (GEO_MAX_CU_LOG2 + 1)) +
                  (int32_t)((uint32_t)g_Dis[dist_y] << (GEO_MAX_CU_LOG2 + 1));

    static const int16_t offset_msk = (2 * GEO_MAX_CU_SIZE - GEO_WEIGHT_MASK_SIZE) >> 1;
    int idx = 0;
    for (y = 0; y < GEO_WEIGHT_MASK_SIZE; y++) {

      int16_t lookUpY = (((y + offset_msk) << 1) + 1) * g_Dis[dist_y];

      for (x = 0; x < GEO_WEIGHT_MASK_SIZE; x++, idx++) {
        int16_t sx_i = ((x + offset_msk) << 1) + 1;
        int16_t weightIdx = sx_i * g_Dis[dist_x] + lookUpY - rho;
        int weightLinearIdx = 32 + weightIdx;

        g_globalGeoWeights[g_angle2mask[angle_idx]][idx] = ov_clip((weightLinearIdx + 4) >> 3, 0, 8);
      }
    }
  }

  for( int hIdx = 0; hIdx < GEO_NUM_CU_SIZE; hIdx++ )
  {
    int16_t height = 1 << ( hIdx + GEO_MIN_CU_LOG2);
    for( int wIdx = 0; wIdx < GEO_NUM_CU_SIZE; wIdx++ ){
      int16_t width = 1 << (wIdx + GEO_MIN_CU_LOG2);
      for( int splitDir = 0; splitDir < GEO_NUM_PARTITION_MODE; splitDir++ ){
        int16_t angle         = g_GeoParams[splitDir][0];
        int16_t distance      = g_GeoParams[splitDir][1];
        int16_t offsetX       = (GEO_WEIGHT_MASK_SIZE - width) >> 1;
        int16_t offsetY       = (GEO_WEIGHT_MASK_SIZE - height) >> 1;
        if( distance > 0 ){
          if( angle % 16 == 8 || (angle % 16 != 0 && height >= width) ){
            offsetY += angle < 16 ? ((distance * (int32_t)height) >> 3) : -((distance * (int32_t)height) >> 3);
          }
          else{
            offsetX += angle < 16 ? ((distance * (int32_t)width) >> 3) : -((distance * (int32_t)width) >> 3);
          }
        }
        g_weightOffset[splitDir][hIdx][wIdx][0] = offsetX;
        g_weightOffset[splitDir][hIdx][wIdx][1] = offsetY;
      }
    }
  }
}
