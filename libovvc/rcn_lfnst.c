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

#include <string.h>

#include "ovutils.h"

#include "rcn_structures.h"

static void
compute_lfnst_4x4(const int16_t* const src, int16_t* const dst,
                  const int8_t* const lfnst_matrix, int log2_tb_w,
                  int log2_tb_h)
{
    int i, j;
    uint64_t scan_map = 0xfbe7ad369c258140;

    int16_t tmp[16];

    for (int i = 0; i < 16; ++i) {
        tmp[i] = src[scan_map & 0xF];
        scan_map >>= 4;
    }

    uint8_t is_4x4 = (log2_tb_w == log2_tb_h);
    uint8_t log2_tr_s = 3 + !is_4x4;

    for (i = 0; i < 16; ++i) {
        int sum = 0;
        for (j = 0; j < 1 << log2_tr_s; ++j) {
            sum += tmp[j] * lfnst_matrix[i + j * 16];
        }
        dst[(i & 3) + ((i >> 2) << log2_tb_w)] =
            ov_clip((sum + 64) >> 7, -(1 << 15), (1 << 15));
    }
}

static void
compute_lfnst_8x8(const int16_t* const src, int16_t* const dst,
                  const int8_t* const lfnst_matrix, int log2_tb_w,
                  int log2_tb_h)
{
    int i, j;
    uint64_t scan_map = 0xfbe7ad369c258140;

    int16_t tmp[16];

    for (int i = 0; i < 16; ++i) {
        tmp[i] = src[scan_map & 0xF];
        scan_map >>= 4;
    }

    for (i = 0; i < 32; ++i) {
        int sum = 0;
        for (j = 0; j < 16; ++j) {
            sum += tmp[j] * lfnst_matrix[i + j * 48];
        }
        dst[(i & 7) + ((i >> 3) << log2_tb_w)] =
            ov_clip((sum + 64) >> 7, -(1 << 15), (1 << 15));
    }

    for (; i < 48; ++i) {
        int sum = 0;
        for (j = 0; j < 16; ++j) {
            sum += tmp[j] * lfnst_matrix[i + j * 48];
        }
        dst[(i & 3) + ((4 + ((i - 32) >> 2)) << log2_tb_w)] =
            ov_clip((sum + 64) >> 7, -(1 << 15), (1 << 15));
    }
}

static void
compute_lfnst_4x4_tr(const int16_t* const src, int16_t* const dst,
                     const int8_t* const lfnst_matrix, int log2_tb_w,
                     int log2_tb_h)
{
    int i, j;
    uint8_t is_4x4 = (log2_tb_w == log2_tb_h);
    uint8_t log2_tr_s = 3 + !is_4x4;

    uint64_t scan_map = 0xfbe7ad369c258140;

    int16_t tmp[16];

    for (int i = 0; i < 16; ++i) {
        tmp[i] = src[scan_map & 0xF];
        scan_map >>= 4;
    }


    for (i = 0; i < 16; ++i) {
        int sum = 0;
        for (j = 0; j < 1 << log2_tr_s; ++j) {
            sum += tmp[j] * lfnst_matrix[i + j * 16];
        }
        dst[((i & 3) << log2_tb_w) + (i >> 2)] =
            ov_clip((sum + 64) >> 7, -(1 << 15), (1 << 15));
    }
}

static void
compute_lfnst_8x8_tr(const int16_t* const src, int16_t* const dst,
                     const int8_t* const lfnst_matrix, int log2_tb_w,
                     int log2_tb_h)
{
    int i, j;
    uint64_t scan_map = 0xfbe7ad369c258140;

    int16_t tmp[16];

    for (int i = 0; i < 16; ++i) {
        tmp[i] = src[scan_map & 0xF];
        scan_map >>= 4;
    }

    for (i = 0; i < 32; ++i) {
        int sum = 0;
        for (j = 0; j < 16; ++j) {
            sum += tmp[j] * lfnst_matrix[i + j * 48];
        }
        dst[((i & 7) << log2_tb_w) + (i >> 3)] =
            ov_clip((sum + 64) >> 7, -(1 << 15), (1 << 15));
    }
    for (; i < 48; ++i) {
        int sum = 0;
        for (j = 0; j < 16; ++j) {
            sum += tmp[j] * lfnst_matrix[i + j * 48];
        }
        dst[((i & 3) << log2_tb_w) + (4 + ((i - 32) >> 2))] =
            ov_clip((sum + 64) >> 7, -(1 << 15), (1 << 15));
    }
}

void
rcn_init_lfnst_functions(struct RCNFunctions *rcn_func)
{
   rcn_func->lfnst.func[0][0] = &compute_lfnst_4x4;
   rcn_func->lfnst.func[0][1] = &compute_lfnst_8x8;
   rcn_func->lfnst.func[1][0] = &compute_lfnst_4x4_tr;
   rcn_func->lfnst.func[1][1] = &compute_lfnst_8x8_tr;
}
