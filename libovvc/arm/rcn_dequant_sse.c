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
#include <string.h>
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

#include "ovutils.h"
#include "rcn_structures.h"

#define MAX_LOG2_TR_RANGE 15
static int
derive_nb_cols(uint64_t sig_sb_map)
{
    uint8_t num_z = ov_clz64(sig_sb_map | 1);
    uint8_t z_div8 = (num_z >> 3);

    return (8 - z_div8) << 2;
}

static void
dequant_tb_4x4_neg(int16_t *dst, const int16_t *src, int scale, int shift,
                   uint8_t log2_tb_w, uint8_t log2_tb_h, uint64_t sig_sb_map)
{
    int nb_rows = derive_nb_cols(sig_sb_map);
    int nb_cols = 1 << OVMIN(5, log2_tb_w);//derive_nb_rows(sig_sb_map);
    uint8_t src_stride = 1 << (OVMIN(5, log2_tb_w));
    uint8_t dst_stride = 1 << (OVMIN(5, log2_tb_w));

    int padd_size = (1 << OVMIN(5, log2_tb_h)) - nb_rows;

    __m128i add   = _mm_set1_epi16((1 << shift) >> 1);
    __m128i scale_v = _mm_unpacklo_epi16(_mm_set1_epi16(scale), _mm_setzero_si128());
    __m128i z = _mm_setzero_si128();

    /* Force sig_sb_map to one in case of DC coefficient */
    sig_sb_map |= !sig_sb_map;

    if (dst_stride < 8) {
        for (int i = 0; i < nb_rows/4 ; i++) {
            uint8_t sig_sb_row = sig_sb_map >> (i << 3);
            for (int j = 0; j < nb_cols/4 ; j++) {
                int16_t *_dst = dst + (j << 2);
                const int16_t *_src = src + (j << 4);
                if (sig_sb_row & 0x1) {
                    __m128i sb_lo = _mm_load_si128((__m128i *)_src);
                    __m128i sb_hi = _mm_load_si128((__m128i *)_src + 1);

                    __m128i x_lo0 = _mm_unpacklo_epi16(sb_lo, add);
                    __m128i x_lo1 = _mm_unpackhi_epi16(sb_lo, add);
                    __m128i x_hi0 = _mm_unpacklo_epi16(sb_hi, add);
                    __m128i x_hi1 = _mm_unpackhi_epi16(sb_hi, add);

                    x_lo0 = _mm_madd_epi16(x_lo0, scale_v);
                    x_lo1 = _mm_madd_epi16(x_lo1, scale_v);
                    x_hi0 = _mm_madd_epi16(x_hi0, scale_v);
                    x_hi1 = _mm_madd_epi16(x_hi1, scale_v);

                    x_lo0 = _mm_slli_epi32(x_lo0, shift);
                    x_lo1 = _mm_slli_epi32(x_lo1, shift);
                    x_hi0 = _mm_slli_epi32(x_hi0, shift);
                    x_hi1 = _mm_slli_epi32(x_hi1, shift);

                    sb_lo = _mm_packs_epi32(x_lo0, x_lo1);
                    sb_hi = _mm_packs_epi32(x_hi0, x_hi1);

                    _mm_store_si128((__m128i *)_dst, sb_lo);
                    _dst += dst_stride << 1;
                    _mm_store_si128((__m128i *)_dst, sb_hi);
                } else {
                    _mm_store_si128((__m128i *)_dst, z);
                    _dst += dst_stride << 1;
                    _mm_store_si128((__m128i *)_dst, z);
                }
                sig_sb_row >>= 1;
            }
            src += src_stride << 2;
            dst += dst_stride << 2;
        }
    } else {
        __m128i lo_msk = _mm_unpackhi_epi64(_mm_set1_epi16(-1), z);
        for (int i = 0; i < nb_rows/4 ; i++) {
            uint8_t sig_sb_row = sig_sb_map >> (i << 3);
            for (int j = 0; j < nb_cols/8 ; j++) {
                int16_t *_dst = dst + (j << 3);
                const int16_t *_src = src + (j << 5);
                if ((sig_sb_row & 0x3) == 3) {
                    __m128i dst_l0, dst_l1, dst_l2, dst_l3;
                    __m128i sb0_lo = _mm_load_si128((__m128i *)_src);
                    __m128i sb0_hi = _mm_load_si128((__m128i *)_src + 1);

                    __m128i sb1_lo = _mm_load_si128((__m128i *)_src + 2);
                    __m128i sb1_hi = _mm_load_si128((__m128i *)_src + 3);

                    __m128i x_lo0 = _mm_unpacklo_epi16(sb0_lo, add);
                    __m128i x_lo1 = _mm_unpackhi_epi16(sb0_lo, add);
                    __m128i x_hi0 = _mm_unpacklo_epi16(sb0_hi, add);
                    __m128i x_hi1 = _mm_unpackhi_epi16(sb0_hi, add);

                    x_lo0 = _mm_madd_epi16(x_lo0, scale_v);
                    x_lo1 = _mm_madd_epi16(x_lo1, scale_v);
                    x_hi0 = _mm_madd_epi16(x_hi0, scale_v);
                    x_hi1 = _mm_madd_epi16(x_hi1, scale_v);

                    x_lo0 = _mm_slli_epi32(x_lo0, shift);
                    x_lo1 = _mm_slli_epi32(x_lo1, shift);
                    x_hi0 = _mm_slli_epi32(x_hi0, shift);
                    x_hi1 = _mm_slli_epi32(x_hi1, shift);

                    sb0_lo = _mm_packs_epi32(x_lo0, x_lo1);
                    sb0_hi = _mm_packs_epi32(x_hi0, x_hi1);

                    x_lo0 = _mm_unpacklo_epi16(sb1_lo, add);
                    x_lo1 = _mm_unpackhi_epi16(sb1_lo, add);
                    x_hi0 = _mm_unpacklo_epi16(sb1_hi, add);
                    x_hi1 = _mm_unpackhi_epi16(sb1_hi, add);

                    x_lo0 = _mm_madd_epi16(x_lo0, scale_v);
                    x_lo1 = _mm_madd_epi16(x_lo1, scale_v);
                    x_hi0 = _mm_madd_epi16(x_hi0, scale_v);
                    x_hi1 = _mm_madd_epi16(x_hi1, scale_v);

                    x_lo0 = _mm_slli_epi32(x_lo0, shift);
                    x_lo1 = _mm_slli_epi32(x_lo1, shift);
                    x_hi0 = _mm_slli_epi32(x_hi0, shift);
                    x_hi1 = _mm_slli_epi32(x_hi1, shift);

                    sb1_lo = _mm_packs_epi32(x_lo0, x_lo1);
                    sb1_hi = _mm_packs_epi32(x_hi0, x_hi1);

                    dst_l0 = _mm_unpacklo_epi64(sb0_lo, sb1_lo);
                    dst_l1 = _mm_unpackhi_epi64(sb0_lo, sb1_lo);
                    dst_l2 = _mm_unpacklo_epi64(sb0_hi, sb1_hi);
                    dst_l3 = _mm_unpackhi_epi64(sb0_hi, sb1_hi);

                    _mm_store_si128((__m128i *)_dst, dst_l0);
                    _dst += dst_stride;
                    _mm_store_si128((__m128i *)_dst, dst_l1);
                    _dst += dst_stride;
                    _mm_store_si128((__m128i *)_dst, dst_l2);
                    _dst += dst_stride;
                    _mm_store_si128((__m128i *)_dst, dst_l3);
                } else if (sig_sb_row & 0x3) {
                    __m128i dst_l0, dst_l1, dst_l2, dst_l3;
                    int8_t not_first_sb = !(sig_sb_row & 0x1);

                    __m128i sb_lo = _mm_load_si128((__m128i *)_src + (not_first_sb << 1));
                    __m128i sb_hi = _mm_load_si128((__m128i *)_src + (not_first_sb << 1) + 1);
                    __m128i store_msk = _mm_xor_si128(_mm_set1_epi16(-not_first_sb), lo_msk);

                    __m128i x_lo0 = _mm_unpacklo_epi16(sb_lo, add);
                    __m128i x_lo1 = _mm_unpackhi_epi16(sb_lo, add);
                    __m128i x_hi0 = _mm_unpacklo_epi16(sb_hi, add);
                    __m128i x_hi1 = _mm_unpackhi_epi16(sb_hi, add);

                    x_lo0 = _mm_madd_epi16(x_lo0, scale_v);
                    x_lo1 = _mm_madd_epi16(x_lo1, scale_v);
                    x_hi0 = _mm_madd_epi16(x_hi0, scale_v);
                    x_hi1 = _mm_madd_epi16(x_hi1, scale_v);

                    x_lo0 = _mm_slli_epi32(x_lo0, shift);
                    x_lo1 = _mm_slli_epi32(x_lo1, shift);
                    x_hi0 = _mm_slli_epi32(x_hi0, shift);
                    x_hi1 = _mm_slli_epi32(x_hi1, shift);

                    sb_lo = _mm_packs_epi32(x_lo0, x_lo1);
                    sb_hi = _mm_packs_epi32(x_hi0, x_hi1);

                    dst_l0 = _mm_unpacklo_epi64(sb_lo, sb_lo);
                    dst_l1 = _mm_unpackhi_epi64(sb_lo, sb_lo);
                    dst_l2 = _mm_unpacklo_epi64(sb_hi, sb_hi);
                    dst_l3 = _mm_unpackhi_epi64(sb_hi, sb_hi);

                    dst_l0 = _mm_and_si128(dst_l0, store_msk);
                    dst_l1 = _mm_and_si128(dst_l1, store_msk);
                    dst_l2 = _mm_and_si128(dst_l2, store_msk);
                    dst_l3 = _mm_and_si128(dst_l3, store_msk);

                    _mm_store_si128((__m128i *)_dst, dst_l0);
                    _dst += dst_stride;
                    _mm_store_si128((__m128i *)_dst, dst_l1);
                    _dst += dst_stride;
                    _mm_store_si128((__m128i *)_dst, dst_l2);
                    _dst += dst_stride;
                    _mm_store_si128((__m128i *)_dst, dst_l3);
                } else {
                    _mm_store_si128((__m128i *)_dst, z);
                    _dst += dst_stride;
                    _mm_store_si128((__m128i *)_dst, z);
                    _dst += dst_stride;
                    _mm_store_si128((__m128i *)_dst, z);
                    _dst += dst_stride;
                    _mm_store_si128((__m128i *)_dst, z);
                }
                sig_sb_row >>= 2;
            }
            src += src_stride << 2;
            dst += dst_stride << 2;
        }
    }
    memset(dst, 0, sizeof(*dst) * src_stride * padd_size);
}

static void
dequant_tb_4x4(int16_t *dst, const int16_t *src, int scale, int shift,
               uint8_t log2_tb_w, uint8_t log2_tb_h, uint64_t sig_sb_map)
{
    int nb_rows = derive_nb_cols(sig_sb_map);
    int nb_cols = 1 << OVMIN(5, log2_tb_w);//derive_nb_rows(sig_sb_map);
    uint8_t src_stride = 1 << (OVMIN(5, log2_tb_w));
    uint8_t dst_stride = 1 << (OVMIN(5, log2_tb_w));

    int padd_size = (1 << OVMIN(5, log2_tb_h)) - nb_rows;

    /* Force sig_sb_map to one in case of DC coefficient */
    sig_sb_map |= !sig_sb_map;
    __m128i add   = _mm_set1_epi16((1 << shift) >> 1);
    __m128i scale_v = _mm_unpacklo_epi16(_mm_set1_epi16(scale), _mm_set1_epi16(1));
    __m128i z = _mm_setzero_si128();

    if (src_stride > 4) {
        __m128i lo_msk = _mm_unpackhi_epi64(_mm_set1_epi16(-1), z);
        for (int i = 0; i < nb_rows/4 ; i++) {
            uint8_t sig_sb_row = sig_sb_map >> (i << 3);
            for (int j = 0; j < nb_cols/8 ; j++) {
                int16_t *_dst = dst + (j << 3);
                const int16_t *_src = src + (j << 5);
                if ((sig_sb_row & 0x3) == 3) {
                    __m128i dst_l0, dst_l1, dst_l2, dst_l3;
                    __m128i sb0_lo = _mm_load_si128((__m128i *)_src);
                    __m128i sb0_hi = _mm_load_si128((__m128i *)_src + 1);

                    __m128i sb1_lo = _mm_load_si128((__m128i *)_src + 2);
                    __m128i sb1_hi = _mm_load_si128((__m128i *)_src + 3);

                    __m128i x_lo0 = _mm_unpacklo_epi16(sb0_lo, add);
                    __m128i x_lo1 = _mm_unpackhi_epi16(sb0_lo, add);
                    __m128i x_hi0 = _mm_unpacklo_epi16(sb0_hi, add);
                    __m128i x_hi1 = _mm_unpackhi_epi16(sb0_hi, add);

                    x_lo0 = _mm_madd_epi16(x_lo0, scale_v);
                    x_lo1 = _mm_madd_epi16(x_lo1, scale_v);
                    x_hi0 = _mm_madd_epi16(x_hi0, scale_v);
                    x_hi1 = _mm_madd_epi16(x_hi1, scale_v);

                    x_lo0 = _mm_srai_epi32(x_lo0, shift);
                    x_lo1 = _mm_srai_epi32(x_lo1, shift);
                    x_hi0 = _mm_srai_epi32(x_hi0, shift);
                    x_hi1 = _mm_srai_epi32(x_hi1, shift);

                    sb0_lo = _mm_packs_epi32(x_lo0, x_lo1);
                    sb0_hi = _mm_packs_epi32(x_hi0, x_hi1);

                    x_lo0 = _mm_unpacklo_epi16(sb1_lo, add);
                    x_lo1 = _mm_unpackhi_epi16(sb1_lo, add);
                    x_hi0 = _mm_unpacklo_epi16(sb1_hi, add);
                    x_hi1 = _mm_unpackhi_epi16(sb1_hi, add);

                    x_lo0 = _mm_madd_epi16(x_lo0, scale_v);
                    x_lo1 = _mm_madd_epi16(x_lo1, scale_v);
                    x_hi0 = _mm_madd_epi16(x_hi0, scale_v);
                    x_hi1 = _mm_madd_epi16(x_hi1, scale_v);

                    x_lo0 = _mm_srai_epi32(x_lo0, shift);
                    x_lo1 = _mm_srai_epi32(x_lo1, shift);
                    x_hi0 = _mm_srai_epi32(x_hi0, shift);
                    x_hi1 = _mm_srai_epi32(x_hi1, shift);

                    sb1_lo = _mm_packs_epi32(x_lo0, x_lo1);
                    sb1_hi = _mm_packs_epi32(x_hi0, x_hi1);

                    dst_l0 = _mm_unpacklo_epi64(sb0_lo, sb1_lo);
                    dst_l1 = _mm_unpackhi_epi64(sb0_lo, sb1_lo);
                    dst_l2 = _mm_unpacklo_epi64(sb0_hi, sb1_hi);
                    dst_l3 = _mm_unpackhi_epi64(sb0_hi, sb1_hi);

                    _mm_store_si128((__m128i *)_dst, dst_l0);
                    _dst += dst_stride;
                    _mm_store_si128((__m128i *)_dst, dst_l1);
                    _dst += dst_stride;
                    _mm_store_si128((__m128i *)_dst, dst_l2);
                    _dst += dst_stride;
                    _mm_store_si128((__m128i *)_dst, dst_l3);
                } else if (sig_sb_row & 0x3) {
                    __m128i dst_l0, dst_l1, dst_l2, dst_l3;
                    int8_t not_first_sb = !(sig_sb_row & 0x1);

                    __m128i sb_lo = _mm_load_si128((__m128i *)_src + (not_first_sb << 1));
                    __m128i sb_hi = _mm_load_si128((__m128i *)_src + (not_first_sb << 1) + 1);
                    __m128i store_msk = _mm_xor_si128(_mm_set1_epi16(-not_first_sb), lo_msk);

                    __m128i x_lo0 = _mm_unpacklo_epi16(sb_lo, add);
                    __m128i x_lo1 = _mm_unpackhi_epi16(sb_lo, add);
                    __m128i x_hi0 = _mm_unpacklo_epi16(sb_hi, add);
                    __m128i x_hi1 = _mm_unpackhi_epi16(sb_hi, add);

                    x_lo0 = _mm_madd_epi16(x_lo0, scale_v);
                    x_lo1 = _mm_madd_epi16(x_lo1, scale_v);
                    x_hi0 = _mm_madd_epi16(x_hi0, scale_v);
                    x_hi1 = _mm_madd_epi16(x_hi1, scale_v);

                    x_lo0 = _mm_srai_epi32(x_lo0, shift);
                    x_lo1 = _mm_srai_epi32(x_lo1, shift);
                    x_hi0 = _mm_srai_epi32(x_hi0, shift);
                    x_hi1 = _mm_srai_epi32(x_hi1, shift);

                    sb_lo = _mm_packs_epi32(x_lo0, x_lo1);
                    sb_hi = _mm_packs_epi32(x_hi0, x_hi1);

                    dst_l0 = _mm_unpacklo_epi64(sb_lo, sb_lo);
                    dst_l1 = _mm_unpackhi_epi64(sb_lo, sb_lo);
                    dst_l2 = _mm_unpacklo_epi64(sb_hi, sb_hi);
                    dst_l3 = _mm_unpackhi_epi64(sb_hi, sb_hi);

                    dst_l0 = _mm_and_si128(dst_l0, store_msk);
                    dst_l1 = _mm_and_si128(dst_l1, store_msk);
                    dst_l2 = _mm_and_si128(dst_l2, store_msk);
                    dst_l3 = _mm_and_si128(dst_l3, store_msk);

                    _mm_store_si128((__m128i *)_dst, dst_l0);
                    _dst += dst_stride;
                    _mm_store_si128((__m128i *)_dst, dst_l1);
                    _dst += dst_stride;
                    _mm_store_si128((__m128i *)_dst, dst_l2);
                    _dst += dst_stride;
                    _mm_store_si128((__m128i *)_dst, dst_l3);
                } else {
                    _mm_store_si128((__m128i *)_dst, z);
                    _dst += dst_stride;
                    _mm_store_si128((__m128i *)_dst, z);
                    _dst += dst_stride;
                    _mm_store_si128((__m128i *)_dst, z);
                    _dst += dst_stride;
                    _mm_store_si128((__m128i *)_dst, z);
                }
                sig_sb_row >>= 2;
            }
            src += src_stride << 2;
            dst += dst_stride << 2;
        }

    } else {
        for (int i = 0; i < nb_rows/4 ; i++) {
            uint8_t sig_sb_row = sig_sb_map >> (i << 3);
            for (int j = 0; j < nb_cols/4 ; j++) {
                int16_t *_dst = dst + (j << 2);
                const int16_t *_src = src + (j << 4);
                if (sig_sb_row & 0x1) {
                    __m128i sb_lo = _mm_load_si128((__m128i *)_src);
                    __m128i sb_hi = _mm_load_si128((__m128i *)_src + 1);

                    __m128i x_lo0 = _mm_unpacklo_epi16(sb_lo, add);
                    __m128i x_lo1 = _mm_unpackhi_epi16(sb_lo, add);
                    __m128i x_hi0 = _mm_unpacklo_epi16(sb_hi, add);
                    __m128i x_hi1 = _mm_unpackhi_epi16(sb_hi, add);

                    x_lo0 = _mm_madd_epi16(x_lo0, scale_v);
                    x_lo1 = _mm_madd_epi16(x_lo1, scale_v);
                    x_hi0 = _mm_madd_epi16(x_hi0, scale_v);
                    x_hi1 = _mm_madd_epi16(x_hi1, scale_v);

                    x_lo0 = _mm_srai_epi32(x_lo0, shift);
                    x_lo1 = _mm_srai_epi32(x_lo1, shift);
                    x_hi0 = _mm_srai_epi32(x_hi0, shift);
                    x_hi1 = _mm_srai_epi32(x_hi1, shift);

                    sb_lo = _mm_packs_epi32(x_lo0, x_lo1);
                    sb_hi = _mm_packs_epi32(x_hi0, x_hi1);

                    _mm_store_si128((__m128i *)_dst, sb_lo);
                    _dst += dst_stride << 1;
                    _mm_store_si128((__m128i *)_dst, sb_hi);
                } else {
                    _mm_store_si128((__m128i *)_dst, z);
                    _dst += dst_stride << 1;
                    _mm_store_si128((__m128i *)_dst, z);
                }
                sig_sb_row >>= 1;
            }
            src += src_stride << 2;
            dst += dst_stride << 2;
        }
    }

    memset(dst, 0, sizeof(*dst) * src_stride * padd_size);
}

void
rcn_init_dequant_sse(struct RCNFunctions *rcn_funcs)
{
     rcn_funcs->tmp.dequant_tb_4x4 = &dequant_tb_4x4;
     rcn_funcs->tmp.dequant_tb_4x4_neg = &dequant_tb_4x4_neg;
}
