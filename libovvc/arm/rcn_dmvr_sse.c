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

#include "rcn_neon.h"

#include "simde/x86/sse4.2.h"
#include "simde/x86/mmx.h"
#include "simde/x86/sse.h"
#include "simde/x86/sse2.h"
#include "simde/x86/sse3.h"
#include "simde/x86/ssse3.h"
#include "simde/x86/sse4.2.h"
#include "simde/x86/avx.h"
#include "simde/x86/avx2.h"
#include "simde/x86/avx512.h"

#include "rcn_structures.h"

#define DMVR_NUM_ITERATION 2
#define DMVR_SAD_STRIDE  ((2 * DMVR_NUM_ITERATION) + 1)
#define DMVR_NB_IDX (DMVR_SAD_STRIDE * DMVR_SAD_STRIDE)

static const int8_t dmvr_mv_x[25 + 25] = {
    -2, -1, 0, 1, 2,
    -2, -1, 0, 1, 2,
    -2, -1, 0, 1, 2,
    -2, -1, 0, 1, 2,
    -2, -1, 0, 1, 2,
    2, 1, 0, -1, -2,
    2, 1, 0, -1, -2,
    2, 1, 0, -1, -2,
    2, 1, 0, -1, -2,
    2, 1, 0, -1, -2,
};

static const int8_t dmvr_mv_y[25 + 25] = {
    -2, -2, -2, -2, -2,
    -1, -1, -1, -1, -1,
     0,  0,  0,  0,  0,
     1,  1,  1,  1,  1,
     2,  2,  2,  2,  2,
     2,  2,  2,  2,  2,
     1,  1,  1,  1,  1,
     0,  0,  0,  0,  0,
    -1, -1, -1, -1, -1,
    -2, -2, -2, -2, -2,
};

static uint64_t
rcn_dmvr_sad_16(const uint16_t *ref0, const uint16_t *ref1,
             int16_t dmvr_stride, int16_t pb_w, int16_t pb_h)
{
  uint64_t sum = 0;
  int i;
  __m128i zero = _mm_setzero_si128();
  __m128i sumV = _mm_setzero_si128();
  for (i = 0; i < (pb_h >> 1); ++i) {
      __m128i x1 = _mm_loadu_si128((__m128i *)&ref0[0]);
      __m128i x2 = _mm_loadu_si128((__m128i *)&ref0[8]);
      __m128i y1 = _mm_loadu_si128((__m128i *)&ref1[0]);
      __m128i y2 = _mm_loadu_si128((__m128i *)&ref1[8]);
      x1 = _mm_abs_epi16(_mm_sub_epi16(x1, y1));
      x2 = _mm_abs_epi16(_mm_sub_epi16(x2, y2));
      x1 = _mm_add_epi16(x1, x2);
      sumV = _mm_add_epi32(sumV, x1);

      ref0 += dmvr_stride << 1;
      ref1 += dmvr_stride << 1;
  }
  sumV = _mm_add_epi32(_mm_unpacklo_epi16(sumV,zero), _mm_unpackhi_epi16(sumV,zero));
  sumV = _mm_add_epi64(_mm_unpacklo_epi32(sumV,zero), _mm_unpackhi_epi32(sumV,zero));
  sumV = _mm_add_epi64(_mm_unpacklo_epi64(sumV,zero), _mm_unpackhi_epi64(sumV,zero));
  _mm_storel_epi64((__m128i *) &sum, sumV);
  return sum;
}

static uint64_t
rcn_dmvr_sad_8(const uint16_t *ref0, const uint16_t *ref1,
             int16_t dmvr_stride, int16_t pb_w, int16_t pb_h)
{
  uint64_t sum = 0;
  int i;
  __m128i zero = _mm_setzero_si128();
  __m128i sumV = _mm_setzero_si128();
  for (i = 0; i < (pb_h >> 1); ++i) {
      __m128i x1 = _mm_loadu_si128((__m128i *)&ref0[0]);
      __m128i y1 = _mm_loadu_si128((__m128i *)&ref1[0]);
      x1 = _mm_abs_epi16(_mm_sub_epi16(x1, y1));
      sumV = _mm_add_epi32(sumV, x1);

      ref0 += dmvr_stride << 1;
      ref1 += dmvr_stride << 1;
  }
  sumV = _mm_add_epi32(_mm_unpacklo_epi16(sumV,zero), _mm_unpackhi_epi16(sumV,zero));
  sumV = _mm_add_epi64(_mm_unpacklo_epi32(sumV,zero), _mm_unpackhi_epi32(sumV,zero));
  sumV = _mm_add_epi64(_mm_unpacklo_epi64(sumV,zero), _mm_unpackhi_epi64(sumV,zero));
  _mm_storel_epi64((__m128i *) &sum, sumV);
  return sum;
}

/*FIXME return min_dmvr_idx; */
static uint8_t
dmvr_compute_sads_16(const uint16_t *ref0, const uint16_t *ref1,
                  uint64_t *sad_array, int sb_w, int sb_h)
{
    const int32_t stride_l0 = 128 + 4;
    const int32_t stride_l1 = 128 + 4;

    const int16_t *const ref0_start = ref0;
    const int16_t *const ref1_start = ref1;
    uint64_t min_cost = (uint64_t) -1;

    uint8_t idx;
    uint8_t dmvr_idx = 12;

    for (idx = 0; idx < 12; ++idx) {
        ref0 = ref0_start + (int16_t)dmvr_mv_x[idx]
                          + (int16_t)dmvr_mv_y[idx] * stride_l0;

        ref1 = ref1_start - (int16_t)dmvr_mv_x[idx]
                          - (int16_t)dmvr_mv_y[idx] * stride_l1;

        sad_array[idx] = rcn_dmvr_sad_16(ref0, ref1, stride_l1,
                                     sb_w, sb_h);
    }

    for (idx = 13; idx < DMVR_NB_IDX; ++idx) {
        ref0 = ref0_start + (int16_t)dmvr_mv_x[idx]
                          + (int16_t)dmvr_mv_y[idx] * stride_l0;

        ref1 = ref1_start - (int16_t)dmvr_mv_x[idx]
                          - (int16_t)dmvr_mv_y[idx] * stride_l1;

        sad_array[idx] = rcn_dmvr_sad_16(ref0, ref1, stride_l1,
                                     sb_w, sb_h);
    }
    for (idx = 0; idx < DMVR_NB_IDX; ++idx) {
        if (sad_array[idx] < min_cost || (idx == 12 && sad_array[idx] <= min_cost)) {
            min_cost = sad_array[idx];
            dmvr_idx = idx;
        }
    }

    return dmvr_idx;
}

/*FIXME return min_dmvr_idx; */
static uint8_t
dmvr_compute_sads_8(const uint16_t *ref0, const uint16_t *ref1,
                  uint64_t *sad_array, int sb_w, int sb_h)
{
    const int32_t stride_l0 = 128 + 4;
    const int32_t stride_l1 = 128 + 4;

    const int16_t *const ref0_start = ref0;
    const int16_t *const ref1_start = ref1;
    uint64_t min_cost = (uint64_t) -1;

    uint8_t idx;
    uint8_t dmvr_idx = 12;

    for (idx = 0; idx < 12; ++idx) {
        ref0 = ref0_start + (int16_t)dmvr_mv_x[idx]
                          + (int16_t)dmvr_mv_y[idx] * stride_l0;

        ref1 = ref1_start - (int16_t)dmvr_mv_x[idx]
                          - (int16_t)dmvr_mv_y[idx] * stride_l1;

        sad_array[idx] = rcn_dmvr_sad_8(ref0, ref1, stride_l1,
                                     sb_w, sb_h);
    }

    for (idx = 13; idx < DMVR_NB_IDX; ++idx) {
        ref0 = ref0_start + (int16_t)dmvr_mv_x[idx]
                          + (int16_t)dmvr_mv_y[idx] * stride_l0;

        ref1 = ref1_start - (int16_t)dmvr_mv_x[idx]
                          - (int16_t)dmvr_mv_y[idx] * stride_l1;

        sad_array[idx] = rcn_dmvr_sad_8(ref0, ref1, stride_l1,
                                     sb_w, sb_h);
    }
    for (idx = 0; idx < DMVR_NB_IDX; ++idx) {
        if (sad_array[idx] < min_cost || (idx == 12 && sad_array[idx] <= min_cost)) {
            min_cost = sad_array[idx];
            dmvr_idx = idx;
        }
    }

    return dmvr_idx;
}

void
rcn_init_dmvr_functions_sse(struct RCNFunctions *const rcn_funcs)
{
    rcn_funcs->dmvr.sad[0] = &rcn_dmvr_sad_8;
    rcn_funcs->dmvr.sad[1] = &rcn_dmvr_sad_16;

    rcn_funcs->dmvr.computeSB[0] = &dmvr_compute_sads_8;
    rcn_funcs->dmvr.computeSB[1] = &dmvr_compute_sads_16;
}
