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
#include <immintrin.h>

#include "rcn_structures.h"
#include "ovutils.h"

#define MIP_SHIFT 6

static inline void
mip_matmult_8_8(const int16_t *bndy_line, uint16_t *dst,
                const uint8_t *matrix_mip, int16_t input_offset,
                int rnd_mip, uint8_t log2_red_h)
{
    uint8_t y;
    __m256i d[2], tmp[4], m[4], b256, a[4];
    __m128i b;

    __m256i rnd = _mm256_set1_epi32(rnd_mip);
    __m256i add = _mm256_set1_epi32(input_offset);

    b = _mm_load_si128((__m128i *)bndy_line);
    b256 = _mm256_inserti128_si256(_mm256_castsi128_si256(b), b, 1);

    for (y = 0; y < (1 << log2_red_h); y++) {
        d[0] = _mm256_load_si256((__m256i *)&matrix_mip[0 * 8]);
        d[1] = _mm256_load_si256((__m256i *)&matrix_mip[4 * 8]);

        tmp[0] = _mm256_unpacklo_epi8(d[0], _mm256_setzero_si256());
        tmp[1] = _mm256_unpackhi_epi8(d[0], _mm256_setzero_si256());
        tmp[2] = _mm256_unpacklo_epi8(d[1], _mm256_setzero_si256());
        tmp[3] = _mm256_unpackhi_epi8(d[1], _mm256_setzero_si256());

        m[0] = _mm256_madd_epi16(b256, tmp[0]);
        m[1] = _mm256_madd_epi16(b256, tmp[1]);
        m[2] = _mm256_madd_epi16(b256, tmp[2]);
        m[3] = _mm256_madd_epi16(b256, tmp[3]);

        a[0] = _mm256_unpacklo_epi32(m[0], m[1]);
        a[1] = _mm256_unpacklo_epi32(m[2], m[3]);

        a[2] = _mm256_unpackhi_epi32(m[0], m[1]);
        a[3] = _mm256_unpackhi_epi32(m[2], m[3]);

        m[0] = _mm256_add_epi32(a[0], a[2]);
        m[1] = _mm256_add_epi32(a[1], a[3]);

        tmp[0] = _mm256_permute2x128_si256(m[0], m[1], 0x20);
        tmp[1] = _mm256_permute2x128_si256(m[0], m[1], 0x31);

        a[0] = _mm256_unpacklo_epi64(tmp[0], tmp[1]);
        a[1] = _mm256_unpackhi_epi64(tmp[0], tmp[1]);

        m[0] = _mm256_add_epi32(a[0], a[1]);
        
        m[0] = _mm256_add_epi32(m[0], rnd);

        m[0] = _mm256_srai_epi32(m[0], MIP_SHIFT);

        m[0] = _mm256_add_epi32(m[0], add);

        b = _mm_packs_epi32(_mm256_extracti128_si256(m[0],0), _mm256_extracti128_si256(m[0],1));

        b = _mm_min_epi16(b, _mm_set1_epi16(1023));
        b = _mm_max_epi16(b, _mm_setzero_si128());

        _mm_store_si128((__m128i *)dst, b);
        dst        += 8;
        matrix_mip += 8 * 8;
    }
}

static inline void
mip_matmult_8_4(const int16_t *bndy_line, uint16_t *dst,
                const uint8_t *matrix_mip, int16_t input_offset,
                int rnd_mip, uint8_t log2_red_h)
{
    __m128i d[2], a[4], b, m[4], tmp[4];
    uint8_t y;

    __m128i rnd = _mm_set1_epi32(rnd_mip);
    __m128i add = _mm_set1_epi32(input_offset);

    b = _mm_loadu_si128((__m128i *)bndy_line);
    __m256i b256 = _mm256_inserti128_si256(_mm256_castsi128_si256(b), b, 1);

    for (y = 0; y < (1 << log2_red_h); y++) {
        __m256i d256 = _mm256_load_si256((__m256i *)&matrix_mip[0 * 8]);

        __m256i tmpl = _mm256_unpacklo_epi8(d256, _mm256_setzero_si256());
        __m256i tmph = _mm256_unpackhi_epi8(d256, _mm256_setzero_si256());

        __m256i ml = _mm256_madd_epi16(b256, tmpl);
        __m256i mh = _mm256_madd_epi16(b256, tmph);

        __m256i al = _mm256_unpacklo_epi32(ml, mh);
        __m256i ah = _mm256_unpackhi_epi32(ml, mh);

        ml = _mm256_add_epi32(al, ah);

        m[0] = _mm256_extracti128_si256(ml, 0);
        m[1] = _mm256_extracti128_si256(ml, 1);

        a[0] = _mm_unpacklo_epi64(m[0], m[1]);
        a[1] = _mm_unpackhi_epi64(m[0], m[1]);

        m[0] = _mm_add_epi32(a[0], a[1]);

        m[0] = _mm_add_epi32(m[0], rnd);

        m[0] = _mm_srai_epi32(m[0], MIP_SHIFT);

        m[0] = _mm_add_epi32(m[0], add);

        m[0] = _mm_packs_epi32(m[0], m[1]);

        m[0] = _mm_min_epi16(m[0], _mm_set1_epi16(1023));
        m[0] = _mm_max_epi16(m[0], _mm_setzero_si128());

        _mm_storeu_si128((__m128i *)dst, m[0]);

        dst        += 4;
        matrix_mip += 4 * 8;
    }
}

static inline void
mip_matmult_4_4(const int16_t *bndy_line, uint16_t *dst,
                const uint8_t *matrix_mip, int16_t input_offset,
                int rnd_mip, uint8_t log2_red_h)
{
    __m128i d, b, a, m[2], tmp[2];
    uint8_t y;

    __m128i rnd = _mm_set1_epi32(rnd_mip);
    __m128i add = _mm_set1_epi32(input_offset);

    b = _mm_loadu_si128((__m128i *)bndy_line);
    b = _mm_unpacklo_epi64(b,b);

    for (y = 0; y < (1 << log2_red_h); y++) {
        d = _mm_loadu_si128((__m128i *)&matrix_mip[0 * 4]);

        tmp[0] = _mm_unpacklo_epi8(d, _mm_setzero_si128());
        tmp[1] = _mm_unpackhi_epi8(d, _mm_setzero_si128());

        m[0] = _mm_madd_epi16(b, tmp[0]);
        m[1] = _mm_madd_epi16(b, tmp[1]);

        a = _mm_hadd_epi32(m[0], m[1]);

        m[0] = _mm_add_epi32(a, rnd);

        m[0] = _mm_srai_epi32(m[0], MIP_SHIFT);

        m[0] = _mm_add_epi32(m[0], add);

        m[0] = _mm_packs_epi32(m[0], _mm_setzero_si128());

        m[0] = _mm_min_epi16(m[0], _mm_set1_epi16(1023));
        m[0] = _mm_max_epi16(m[0], _mm_setzero_si128());

        _mm_storel_epi64((__m128i *)dst, m[0]);
        dst        += 4;
        matrix_mip += 4 * 4;
    }
}

static void
mip_matmult(const int16_t *bndy_line, uint16_t *dst,
            const uint8_t *matrix_mip, int16_t input_offset,
            int rnd_mip, uint8_t log2_bndy,
            uint8_t log2_red_w, uint8_t log2_red_h)
{

    if (log2_bndy - 1) {
        if (log2_red_w == 3) {
            mip_matmult_8_8(bndy_line, dst, matrix_mip, input_offset, rnd_mip, log2_red_h);
        } else if (log2_red_w == 2) {
            mip_matmult_8_4(bndy_line, dst, matrix_mip, input_offset, rnd_mip, log2_red_h);
        }
    } else {
        mip_matmult_4_4(bndy_line, dst, matrix_mip, input_offset, rnd_mip, log2_red_h);
    }
}

void
rcn_init_mip_functions_avx2(struct RCNFunctions *const rcn_funcs)
{
    rcn_funcs->mip.matmult= &mip_matmult;
}
