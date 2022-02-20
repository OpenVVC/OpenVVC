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

#include <immintrin.h>

#include "rcn_intra_angular.h"
#include "rcn_structures.h"
#include "bitdepth.h"
#include "ovutils.h"
#include "data_rcn_angular.h"

#define LOAD_GAUSS_FILTER_SSE() \
        __m128i filter01 = _mm_set1_epi32((16 - (delta_frac >> 1))&0xFFFF | ((int32_t)(32 - (delta_frac >> 1))<<16)); \
        __m128i filter23 = _mm_set1_epi32((16 + (delta_frac >> 1))&0xFFFF | ((int32_t)(delta_frac >> 1)<<16));

#define LOAD_CUBIC_FILTER_SSE() \
        __m128i filter01 = _mm_set1_epi32(filter[0]&0xFFFF | ((int32_t)filter[1]<<16)); \
        __m128i filter23 = _mm_set1_epi32(filter[2]&0xFFFF | ((int32_t)filter[3]<<16));

#define LOAD_GAUSS_FILTER_AVX2() \
        __m256i filter01 = _mm256_set1_epi32((16 - (delta_frac >> 1))&0xFFFF | ((int32_t)(32 - (delta_frac >> 1))<<16)); \
        __m256i filter23 = _mm256_set1_epi32((16 + (delta_frac >> 1))&0xFFFF | ((int32_t)(delta_frac >> 1)<<16));

#define LOAD_CUBIC_FILTER_AVX2() \
        __m256i filter01 = _mm256_set1_epi32(filter[0]&0xFFFF | ((int32_t)filter[1]<<16)); \
        __m256i filter23 = _mm256_set1_epi32(filter[2]&0xFFFF | ((int32_t)filter[3]<<16));

#define FILTER_16_SAMPLES() \
        __m256i ref0 = _mm256_loadu_si256((__m256i *)&ref[0]);          \
        __m256i ref1 = _mm256_loadu_si256((__m256i *)&ref[1]);          \
        __m256i ref2 = _mm256_loadu_si256((__m256i *)&ref[2]);          \
        __m256i ref3 = _mm256_loadu_si256((__m256i *)&ref[3]);          \
                                                                     \
        __m256i ref01L = _mm256_unpacklo_epi16(ref0, ref1);             \
        __m256i ref01H = _mm256_unpackhi_epi16(ref0, ref1);             \
        __m256i ref23L = _mm256_unpacklo_epi16(ref2, ref3);             \
        __m256i ref23H = _mm256_unpackhi_epi16(ref2, ref3);             \
                                                                     \
        ref01L = _mm256_madd_epi16(ref01L, filter01);                   \
        ref01H = _mm256_madd_epi16(ref01H, filter01);                   \
        ref23L = _mm256_madd_epi16(ref23L, filter23);                   \
        ref23H = _mm256_madd_epi16(ref23H, filter23);                   \
                                                                     \
        ref01L = _mm256_add_epi32(ref01L, ref23L);                      \
        ref01H = _mm256_add_epi32(ref01H, ref23H);                      \
                                                                     \
        ref01L = _mm256_add_epi32(ref01L, offset);                      \
        ref01H = _mm256_add_epi32(ref01H, offset);                      \
                                                                     \
        ref01L = _mm256_srai_epi32(ref01L, 6);                          \
        ref01H = _mm256_srai_epi32(ref01H, 6);                          \
                                                                     \
        ref01L = _mm256_min_epi32(ref01L, _mm256_set1_epi32(1023));        \
        ref01H = _mm256_min_epi32(ref01H, _mm256_set1_epi32(1023));        \
                                                                     \
        ref01L = _mm256_max_epi32(ref01L, _mm256_setzero_si256());         \
        ref01H = _mm256_max_epi32(ref01H, _mm256_setzero_si256());         \
                                                                     \
        ref01L = _mm256_packs_epi32(ref01L, ref01H);                    \
                                                                     \
        _mm256_storeu_si256((__m256i *)&(_dst[x]), ref01L);

#define FILTER_8_SAMPLES() \
        __m128i ref0 = _mm_loadu_si128((__m128i *)&ref[0]);          \
        __m128i ref1 = _mm_loadu_si128((__m128i *)&ref[1]);          \
        __m128i ref2 = _mm_loadu_si128((__m128i *)&ref[2]);          \
        __m128i ref3 = _mm_loadu_si128((__m128i *)&ref[3]);          \
                                                                     \
        __m128i ref01L = _mm_unpacklo_epi16(ref0, ref1);             \
        __m128i ref01H = _mm_unpackhi_epi16(ref0, ref1);             \
        __m128i ref23L = _mm_unpacklo_epi16(ref2, ref3);             \
        __m128i ref23H = _mm_unpackhi_epi16(ref2, ref3);             \
                                                                     \
        ref01L = _mm_madd_epi16(ref01L, filter01);                   \
        ref01H = _mm_madd_epi16(ref01H, filter01);                   \
        ref23L = _mm_madd_epi16(ref23L, filter23);                   \
        ref23H = _mm_madd_epi16(ref23H, filter23);                   \
                                                                     \
        ref01L = _mm_add_epi32(ref01L, ref23L);                      \
        ref01H = _mm_add_epi32(ref01H, ref23H);                      \
                                                                     \
        ref01L = _mm_add_epi32(ref01L, offset);                      \
        ref01H = _mm_add_epi32(ref01H, offset);                      \
                                                                     \
        ref01L = _mm_srai_epi32(ref01L, 6);                          \
        ref01H = _mm_srai_epi32(ref01H, 6);                          \
                                                                     \
        ref01L = _mm_min_epi32(ref01L, _mm_set1_epi32(1023));        \
        ref01H = _mm_min_epi32(ref01H, _mm_set1_epi32(1023));        \
                                                                     \
        ref01L = _mm_max_epi32(ref01L, _mm_setzero_si128());         \
        ref01H = _mm_max_epi32(ref01H, _mm_setzero_si128());         \
                                                                     \
        ref01L = _mm_packs_epi32(ref01L, ref01H);                    \
                                                                     \
        _mm_storeu_si128((__m128i *)&(_dst[x]), ref01L);

#define FILTER_4_SAMPLES() \
        __m128i ref0 = _mm_loadl_epi64((__m128i *)&ref[0]);        \
        __m128i ref1 = _mm_loadl_epi64((__m128i *)&ref[1]);        \
        __m128i ref2 = _mm_loadl_epi64((__m128i *)&ref[2]);        \
        __m128i ref3 = _mm_loadl_epi64((__m128i *)&ref[3]);        \
                                                                   \
        __m128i ref01L = _mm_unpacklo_epi16(ref0, ref1);           \
        __m128i ref23L = _mm_unpacklo_epi16(ref2, ref3);           \
                                                                   \
        ref01L = _mm_madd_epi16(ref01L, filter01);                 \
        ref23L = _mm_madd_epi16(ref23L, filter23);                 \
                                                                   \
        ref01L = _mm_add_epi32(ref01L, ref23L);                    \
                                                                   \
        ref01L = _mm_add_epi32(ref01L, offset);                    \
                                                                   \
        ref01L = _mm_srai_epi32(ref01L, 6);                        \
                                                                   \
        ref01L = _mm_min_epi32(ref01L, _mm_set1_epi32(1023));      \
                                                                   \
        ref01L = _mm_max_epi32(ref01L, _mm_setzero_si128());       \
                                                                   \
        ref01L = _mm_packs_epi32(ref01L, _mm_setzero_si128());     \
                                                                   \
        _mm_storel_epi64((__m128i *)&(_dst[0]), ref01L);

static const int8_t chroma_filter[4 * 32] =
{
     0, 64,  0,  0,
    -1, 63,  2,  0,
    -2, 62,  4,  0,
    -2, 60,  7, -1,
    -2, 58, 10, -2,
    -3, 57, 12, -2,
    -4, 56, 14, -2,
    -4, 55, 15, -2,
    -4, 54, 16, -2,
    -5, 53, 18, -2,
    -6, 52, 20, -2,
    -6, 49, 24, -3,
    -6, 46, 28, -4,
    -5, 44, 29, -4,
    -4, 42, 30, -4,
    -4, 39, 33, -4,
    -4, 36, 36, -4,
    -4, 33, 39, -4,
    -4, 30, 42, -4,
    -4, 29, 44, -5,
    -4, 28, 46, -6,
    -3, 24, 49, -6,
    -2, 20, 52, -6,
    -2, 18, 53, -5,
    -2, 16, 54, -4,
    -2, 15, 55, -4,
    -2, 14, 56, -4,
    -2, 12, 57, -3,
    -2, 10, 58, -2,
    -1,  7, 60, -2,
     0,  4, 62, -2,
     0,  2, 63, -1
};

static void
intra_angular_v_gauss_avx2_16(const OVSample* ref_abv, OVSample* dst,
                      ptrdiff_t dst_stride, int8_t log2_pb_w,
                      int8_t log2_pb_h, int angle_val)
{
    int delta_pos = angle_val;
    int pb_w = 1 << log2_pb_w;
    int pb_h = 1 << log2_pb_h;
    OVSample* _dst = dst;
    __m256i offset = _mm256_set1_epi32(32);

    for (int y = 0; y < pb_h; y++) {
        const int delta_int = delta_pos >> 5;
        const int delta_frac = delta_pos & 0x1F;

        const OVSample* ref = (OVSample*)ref_abv + delta_int;

        LOAD_GAUSS_FILTER_AVX2();
        for (int x = 0; x < pb_w; x+=16) {
            FILTER_16_SAMPLES();
            ref+=16;
        }

        delta_pos += angle_val;
        _dst += dst_stride;
    }
}

static void
intra_angular_v_gauss_avx2_8(const OVSample* ref_abv, OVSample* dst,
                      ptrdiff_t dst_stride, int8_t log2_pb_w,
                      int8_t log2_pb_h, int angle_val)
{
    int delta_pos = angle_val;
    int pb_w = 1 << log2_pb_w;
    int pb_h = 1 << log2_pb_h;
    OVSample* _dst = dst;
    __m128i offset = _mm_set1_epi32(32);

    for (int y = 0; y < pb_h; y++) {
        const int delta_int = delta_pos >> 5;
        const int delta_frac = delta_pos & 0x1F;

        const OVSample* ref = (OVSample*)ref_abv + delta_int;

        LOAD_GAUSS_FILTER_SSE();
        for (int x = 0; x < pb_w; x+=8) {
            FILTER_8_SAMPLES();
            ref+=8;
        }

        delta_pos += angle_val;
        _dst += dst_stride;
    }
}

static void
intra_angular_v_gauss_avx2_4(const OVSample* ref_abv, OVSample* dst,
                      ptrdiff_t dst_stride, int8_t log2_pb_w,
                      int8_t log2_pb_h, int angle_val)
{
    int delta_pos = angle_val;
    int pb_w = 1 << log2_pb_w;
    int pb_h = 1 << log2_pb_h;
    OVSample* _dst = dst;
    __m128i offset = _mm_set1_epi32(32);

    for (int y = 0; y < pb_h; y++) {
        const int delta_int = delta_pos >> 5;
        const int delta_frac = delta_pos & 0x1F;

        const OVSample* ref = (OVSample*)ref_abv + delta_int;
        LOAD_GAUSS_FILTER_SSE();

        FILTER_4_SAMPLES();

        delta_pos += angle_val;
        _dst += dst_stride;
    }
}

static void
intra_angular_v_gauss_avx2(const OVSample* ref_abv, OVSample* dst,
                      ptrdiff_t dst_stride, int8_t log2_pb_w,
                      int8_t log2_pb_h, int angle_val)
{
    if (log2_pb_w >=4){
        intra_angular_v_gauss_avx2_16(ref_abv, dst, dst_stride, log2_pb_w, log2_pb_h,angle_val);
    }
    else if (log2_pb_w ==3){
        intra_angular_v_gauss_avx2_8(ref_abv, dst, dst_stride, log2_pb_w, log2_pb_h,angle_val);
    }
    else{
        intra_angular_v_gauss_avx2_4(ref_abv, dst, dst_stride, log2_pb_w, log2_pb_h,angle_val);
    }
    
}

static void
intra_angular_h_gauss_avx2(const OVSample* ref_lft, OVSample* dst,
                      ptrdiff_t dst_stride, int8_t log2_pb_w,
                      int8_t log2_pb_h, int angle_val)
{
    OVSample tmp_dst[128 * 128];
    const int tmp_stride = 128;
    OVSample* _tmp = tmp_dst;
    OVSample* _dst = dst;
    int pb_w = 1 << log2_pb_w;
    int pb_h = 1 << log2_pb_h;

    if (log2_pb_h >=4){
        intra_angular_v_gauss_avx2_16(ref_lft, _tmp, tmp_stride, log2_pb_h, log2_pb_w,angle_val);
    }
    else if (log2_pb_h ==3){
        intra_angular_v_gauss_avx2_8(ref_lft, _tmp, tmp_stride, log2_pb_h, log2_pb_w,angle_val);
    }
    else{
        intra_angular_v_gauss_avx2_4(ref_lft, _tmp, tmp_stride, log2_pb_h, log2_pb_w,angle_val);
    }

    for (int y = 0; y < pb_w; y++) {
        _dst = &dst[y];
        for (int x = 0; x < pb_h; x++) {
            _dst[0] = _tmp[x];
            _dst += dst_stride;
        }
        _tmp += tmp_stride;
    }
}

static void
intra_angular_v_cubic_avx2_16(const OVSample* ref_abv, OVSample* dst,
                      ptrdiff_t dst_stride, int8_t log2_pb_w,
                      int8_t log2_pb_h, int angle_val)
{
    int delta_pos = angle_val;
    int pb_w = 1 << log2_pb_w;
    int pb_h = 1 << log2_pb_h;
    OVSample* _dst = dst;
    __m256i offset = _mm256_set1_epi32(32);

    for (int y = 0; y < pb_h; y++) {
        const int delta_int = delta_pos >> 5;
        const int delta_frac = delta_pos & 0x1F;

        const OVSample* ref = (OVSample*)ref_abv + delta_int;
        const int8_t* filter = &chroma_filter[delta_frac << 2];
        LOAD_CUBIC_FILTER_AVX2();
        for (int x = 0; x < pb_w; x+=16) {
            FILTER_16_SAMPLES();
            ref+=16;
        }

        delta_pos += angle_val;
        _dst += dst_stride;
    }
}

static void
intra_angular_v_cubic_avx2_8(const OVSample* ref_abv, OVSample* dst,
                      ptrdiff_t dst_stride, int8_t log2_pb_w,
                      int8_t log2_pb_h, int angle_val)
{
    int delta_pos = angle_val;
    int pb_w = 1 << log2_pb_w;
    int pb_h = 1 << log2_pb_h;
    OVSample* _dst = dst;
    __m128i offset = _mm_set1_epi32(32);

    for (int y = 0; y < pb_h; y++) {
        const int delta_int = delta_pos >> 5;
        const int delta_frac = delta_pos & 0x1F;

        const OVSample* ref = (OVSample*)ref_abv + delta_int;
        const int8_t* filter = &chroma_filter[delta_frac << 2];
        LOAD_CUBIC_FILTER_SSE();
        for (int x = 0; x < pb_w; x+=8) {
            FILTER_8_SAMPLES();
            ref+=8;
        }

        delta_pos += angle_val;
        _dst += dst_stride;
    }
}

static void
intra_angular_v_cubic_avx2_4(const OVSample* ref_abv, OVSample* dst,
                      ptrdiff_t dst_stride, int8_t log2_pb_w,
                      int8_t log2_pb_h, int angle_val)
{
    int delta_pos = angle_val;
    int pb_w = 1 << log2_pb_w;
    int pb_h = 1 << log2_pb_h;
    OVSample* _dst = dst;
    __m128i offset = _mm_set1_epi32(32);

    for (int y = 0; y < pb_h; y++) {
        const int delta_int = delta_pos >> 5;
        const int delta_frac = delta_pos & 0x1F;

        const OVSample* ref = (OVSample*)ref_abv + delta_int;
        const int8_t* filter = &chroma_filter[delta_frac << 2];
        LOAD_CUBIC_FILTER_SSE();

        FILTER_4_SAMPLES();

        delta_pos += angle_val;
        _dst += dst_stride;
    }
}

static void
intra_angular_v_cubic_avx2(const OVSample* ref_abv, OVSample* dst,
                      ptrdiff_t dst_stride, int8_t log2_pb_w,
                      int8_t log2_pb_h, int angle_val)
{
    if (log2_pb_w >= 4){
        intra_angular_v_cubic_avx2_16(ref_abv, dst, dst_stride, log2_pb_w, log2_pb_h,angle_val);
    }
    else if (log2_pb_w == 3){
        intra_angular_v_cubic_avx2_8(ref_abv, dst, dst_stride, log2_pb_w, log2_pb_h,angle_val);
    }
    else{
        intra_angular_v_cubic_avx2_4(ref_abv, dst, dst_stride, log2_pb_w, log2_pb_h,angle_val);
    }
    
}

static void
intra_angular_h_cubic_avx2(const OVSample* ref_lft, OVSample* dst,
                      ptrdiff_t dst_stride, int8_t log2_pb_w,
                      int8_t log2_pb_h, int angle_val)
{
    OVSample tmp_dst[128 * 128];
    const int tmp_stride = 128;
    OVSample* _tmp = tmp_dst;
    OVSample* _dst = dst;
    int pb_w = 1 << log2_pb_w;
    int pb_h = 1 << log2_pb_h;

    if (log2_pb_h >= 4){
        intra_angular_v_cubic_avx2_16(ref_lft, _tmp, tmp_stride, log2_pb_h, log2_pb_w,angle_val);
    }
    else if (log2_pb_h == 3){
        intra_angular_v_cubic_avx2_8(ref_lft, _tmp, tmp_stride, log2_pb_h, log2_pb_w,angle_val);
    }
    else{
        intra_angular_v_cubic_avx2_4(ref_lft, _tmp, tmp_stride, log2_pb_h, log2_pb_w,angle_val);
    }

    for (int y = 0; y < pb_w; y++) {
        _dst = &dst[y];
        for (int x = 0; x < pb_h; x++) {
            _dst[0] = _tmp[x];
            _dst += dst_stride;
        }
        _tmp += tmp_stride;
    }
}

static void
intra_angular_v_c_avx2_16(const OVSample* ref_abv, OVSample* dst,
                  ptrdiff_t dst_stride, int8_t log2_pb_w,
                  int8_t log2_pb_h, int angle_val)
{
    OVSample* _dst = dst;
    int pb_w = 1 << log2_pb_w;
    int pb_h = 1 << log2_pb_h;
    __m256i offset = _mm256_set1_epi32(16);

    for (int y = 0, delta_pos = angle_val; y < pb_h; y++) {
        const int delta_int = delta_pos >> 5;
        __m256i delta_frac = _mm256_set1_epi16(delta_pos & 0x1F);
        const OVSample* ref = ref_abv + delta_int + 1;
        for (int x = 0; x < pb_w; ref+=16, x+=16) {
            __m256i last_ref_val = _mm256_loadu_si256((__m256i *)&ref[0]);
            __m256i curr_ref_val = _mm256_loadu_si256((__m256i *)&ref[1]);
            __m256i val = _mm256_sub_epi16(curr_ref_val, last_ref_val);
            
            __m256i valL = _mm256_mullo_epi16(delta_frac, val);
            __m256i valH = _mm256_mulhi_epi16(delta_frac, val);

            __m256i val03 = _mm256_unpacklo_epi16(valL, valH);
            __m256i val47 = _mm256_unpackhi_epi16(valL, valH);

            val03 = _mm256_add_epi32(val03, offset);
            val47 = _mm256_add_epi32(val47, offset);

            val03 = _mm256_srai_epi32(val03, 5);
            val47 = _mm256_srai_epi32(val47, 5);

            // FIXME PACKS before clipping might cause errors
            val = _mm256_packs_epi32(val03, val47);

            val = _mm256_add_epi16(val, last_ref_val);

            val = _mm256_min_epi16(val, _mm256_set1_epi16(1023));
            val = _mm256_max_epi16(val, _mm256_setzero_si256());
            _mm256_storeu_si256((__m256i *) &_dst[x], val);
        }
        delta_pos += angle_val;
        _dst += dst_stride;
    }
}

static void
intra_angular_v_c_avx2_8(const OVSample* ref_abv, OVSample* dst,
                  ptrdiff_t dst_stride, int8_t log2_pb_w,
                  int8_t log2_pb_h, int angle_val)
{
    OVSample* _dst = dst;
    int pb_w = 1 << log2_pb_w;
    int pb_h = 1 << log2_pb_h;
    __m128i offset = _mm_set1_epi32(16);

    for (int y = 0, delta_pos = angle_val; y < pb_h; y++) {
        const int delta_int = delta_pos >> 5;
        __m128i delta_frac = _mm_set1_epi16(delta_pos & 0x1F);
        const OVSample* ref = ref_abv + delta_int + 1;
        for (int x = 0; x < pb_w; ref+=8, x+=8) {
            __m128i last_ref_val = _mm_loadu_si128((__m128i *)&ref[0]);
            __m128i curr_ref_val = _mm_loadu_si128((__m128i *)&ref[1]);
            __m128i val = _mm_sub_epi16(curr_ref_val, last_ref_val);
            
            __m128i valL = _mm_mullo_epi16(delta_frac, val);
            __m128i valH = _mm_mulhi_epi16(delta_frac, val);

            __m128i val03 = _mm_unpacklo_epi16(valL, valH);
            __m128i val47 = _mm_unpackhi_epi16(valL, valH);

            val03 = _mm_add_epi32(val03, offset);
            val47 = _mm_add_epi32(val47, offset);

            val03 = _mm_srai_epi32(val03, 5);
            val47 = _mm_srai_epi32(val47, 5);

            // FIXME PACKS before clipping might cause errors
            val = _mm_packs_epi32(val03, val47);

            val = _mm_add_epi16(val, last_ref_val);

            val = _mm_min_epi16(val, _mm_set1_epi16(1023));
            val = _mm_max_epi16(val, _mm_setzero_si128());
            _mm_storeu_si128((__m128i *) &_dst[x], val);
        }
        delta_pos += angle_val;
        _dst += dst_stride;
    }
}

static void
intra_angular_v_c_avx2_4(const OVSample* ref_abv, OVSample* dst,
                  ptrdiff_t dst_stride, int8_t log2_pb_w,
                  int8_t log2_pb_h, int angle_val)
{
    OVSample* _dst = dst;
    int pb_w = 1 << log2_pb_w;
    int pb_h = 1 << log2_pb_h;
    __m128i offset = _mm_set1_epi32(16);

    for (int y = 0, delta_pos = angle_val; y < pb_h; y++) {
        const int delta_int = delta_pos >> 5;
        __m128i delta_frac = _mm_set1_epi16(delta_pos & 0x1F);
        const OVSample* ref = ref_abv + delta_int + 1;
        __m128i last_ref_val = _mm_loadl_epi64((__m128i *)&ref[0]);
        __m128i curr_ref_val = _mm_loadl_epi64((__m128i *)&ref[1]);
        __m128i val = _mm_sub_epi16(curr_ref_val, last_ref_val);
        
        __m128i valL = _mm_mullo_epi16(delta_frac, val);
        __m128i valH = _mm_mulhi_epi16(delta_frac, val);

        __m128i val03 = _mm_unpacklo_epi16(valL, valH);

        val03 = _mm_add_epi32(val03, offset);

        val03 = _mm_srai_epi32(val03, 5);

        // FIXME PACKS before clipping might cause errors
        val = _mm_packs_epi32(val03, _mm_setzero_si128());

        val = _mm_add_epi16(val, last_ref_val);

        val = _mm_min_epi16(val, _mm_set1_epi16(1023));
        val = _mm_max_epi16(val, _mm_setzero_si128());
        _mm_storel_epi64((__m128i *) &_dst[0], val);
        delta_pos += angle_val;
        _dst += dst_stride;
    }
}

static void
intra_angular_v_c_avx2(const OVSample* ref_abv, OVSample* dst,
                  ptrdiff_t dst_stride, int8_t log2_pb_w,
                  int8_t log2_pb_h, int angle_val)
{
    OVSample* _dst = dst;
    int pb_w = 1 << log2_pb_w;
    int pb_h = 1 << log2_pb_h;

    if (pb_w >=16){
        intra_angular_v_c_avx2_16(ref_abv, dst, dst_stride, log2_pb_w, log2_pb_h, angle_val);
    } else if (pb_w ==8){
        intra_angular_v_c_avx2_8(ref_abv, dst, dst_stride, log2_pb_w, log2_pb_h, angle_val);
    } else {
        intra_angular_v_c_avx2_4(ref_abv, dst, dst_stride, log2_pb_w, log2_pb_h, angle_val);
    }
}

static void
intra_angular_h_c_avx2(const OVSample* ref_lft, OVSample* dst,
                  ptrdiff_t dst_stride, int8_t log2_pb_w,
                  int8_t log2_pb_h, int angle_val)
{
    OVSample tmp_dst[128 * 128];
    OVSample* _tmp = tmp_dst;
    OVSample* _dst = dst;
    int pb_w = 1 << log2_pb_w;
    int pb_h = 1 << log2_pb_h;
    const int tmp_stride = 128;

    if (pb_h >=16){
        intra_angular_v_c_avx2_16(ref_lft, tmp_dst, tmp_stride, log2_pb_h, log2_pb_w, angle_val);
    } else if (pb_h ==8){
        intra_angular_v_c_avx2_8(ref_lft, tmp_dst, tmp_stride, log2_pb_h, log2_pb_w, angle_val);
    } else {
        intra_angular_v_c_avx2_4(ref_lft, tmp_dst, tmp_stride, log2_pb_h, log2_pb_w, angle_val);
    }

    for (int y = 0; y < pb_w; y++) {
        _dst = &dst[y];
        for (int x = 0; x < pb_h; x++) {
            _dst[0] = _tmp[x];
            _dst += dst_stride;
        }
        _tmp += tmp_stride;
    }

}

static void
intra_angular_v_gauss_pdpc_avx2_16(const OVSample* ref_abv, const OVSample* ref_lft,
                           OVSample* const dst, ptrdiff_t dst_stride,
                           int8_t log2_pb_w, int8_t log2_pb_h, int mode_idx)
{
    OVSample* _dst = dst;
    int angle_val = angle_table[mode_idx];
    int inv_angle = inverse_angle_table[mode_idx];
    int pb_w = 1 << log2_pb_w;
    int pb_h = 1 << log2_pb_h;
    int delta_pos = angle_val;
    int scale = OVMIN(2, log2_pb_h - (floor_log2(3 * inv_angle - 2) - 8));
    __m256i offset = _mm256_set1_epi32(32);
    for (int y = 0; y < pb_h; y++) {
        const int delta_int  = delta_pos >> 5;
        const int delta_frac = delta_pos & 0x1F;
        int inv_angle_sum = 256 + inv_angle;
        const OVSample* ref = (OVSample*)ref_abv + delta_int;
        LOAD_GAUSS_FILTER_AVX2();
        for (int x = 0; x < pb_w; x+=16) {
            FILTER_16_SAMPLES();
            ref+=16;
        }
        for (int x = 0; x < OVMIN(3 << scale, pb_w); x++) {
            int wL = 32 >> ((x << 1) >> scale);
            const OVSample* p = ref_lft + y + (inv_angle_sum >> 9) + 1;

            int16_t left = p[0];
            _dst[x] =
                ov_bdclip(_dst[x] + ((wL * (left - _dst[x]) + 32) >> 6));
            inv_angle_sum += inv_angle;
        }
        delta_pos += angle_val;
        _dst += dst_stride;
    }
}

static void
intra_angular_v_gauss_pdpc_avx2_8(const OVSample* ref_abv, const OVSample* ref_lft,
                           OVSample* const dst, ptrdiff_t dst_stride,
                           int8_t log2_pb_w, int8_t log2_pb_h, int mode_idx)
{
    OVSample* _dst = dst;
    int angle_val = angle_table[mode_idx];
    int inv_angle = inverse_angle_table[mode_idx];
    int pb_w = 1 << log2_pb_w;
    int pb_h = 1 << log2_pb_h;
    int delta_pos = angle_val;
    int scale = OVMIN(2, log2_pb_h - (floor_log2(3 * inv_angle - 2) - 8));
    __m128i offset = _mm_set1_epi32(32);
    for (int y = 0; y < pb_h; y++) {
        const int delta_int  = delta_pos >> 5;
        const int delta_frac = delta_pos & 0x1F;
        int inv_angle_sum = 256 + inv_angle;
        const OVSample* ref = (OVSample*)ref_abv + delta_int;
        LOAD_GAUSS_FILTER_SSE();
        for (int x = 0; x < pb_w; x+=8) {
            FILTER_8_SAMPLES();
            ref+=8;
        }
        for (int x = 0; x < OVMIN(3 << scale, pb_w); x++) {
            int wL = 32 >> ((x << 1) >> scale);
            const OVSample* p = ref_lft + y + (inv_angle_sum >> 9) + 1;

            int16_t left = p[0];
            _dst[x] =
                ov_bdclip(_dst[x] + ((wL * (left - _dst[x]) + 32) >> 6));
            inv_angle_sum += inv_angle;
        }
        delta_pos += angle_val;
        _dst += dst_stride;
    }
}

static void
intra_angular_v_gauss_pdpc_avx2_4(const OVSample* ref_abv, const OVSample* ref_lft,
                           OVSample* const dst, ptrdiff_t dst_stride,
                           int8_t log2_pb_w, int8_t log2_pb_h, int mode_idx)
{
    OVSample* _dst = dst;
    int angle_val = angle_table[mode_idx];
    int inv_angle = inverse_angle_table[mode_idx];
    int pb_w = 1 << log2_pb_w;
    int pb_h = 1 << log2_pb_h;
    int delta_pos = angle_val;
    int scale = OVMIN(2, log2_pb_h - (floor_log2(3 * inv_angle - 2) - 8));
    __m128i offset = _mm_set1_epi32(32);
    for (int y = 0; y < pb_h; y++) {
        const int delta_int  = delta_pos >> 5;
        const int delta_frac = delta_pos & 0x1F;
        int inv_angle_sum = 256 + inv_angle;
        const OVSample* ref = (OVSample*)ref_abv + delta_int;
        LOAD_GAUSS_FILTER_SSE();
        FILTER_4_SAMPLES();
        for (int x = 0; x < OVMIN(3 << scale, pb_w); x++) {
            int wL = 32 >> ((x << 1) >> scale);
            const OVSample* p = ref_lft + y + (inv_angle_sum >> 9) + 1;

            int16_t left = p[0];
            _dst[x] =
                ov_bdclip(_dst[x] + ((wL * (left - _dst[x]) + 32) >> 6));
            inv_angle_sum += inv_angle;
        }
        delta_pos += angle_val;
        _dst += dst_stride;
    }
}


static void
intra_angular_v_gauss_pdpc_avx2(const OVSample* ref_abv, const OVSample* ref_lft,
                           OVSample* const dst, ptrdiff_t dst_stride,
                           int8_t log2_pb_w, int8_t log2_pb_h, int mode_idx)
{
    if (log2_pb_w >=4){
        intra_angular_v_gauss_pdpc_avx2_16(ref_abv, ref_lft, dst, dst_stride, log2_pb_w, log2_pb_h, mode_idx);
    }
    else if (log2_pb_w ==3){
        intra_angular_v_gauss_pdpc_avx2_8(ref_abv, ref_lft, dst, dst_stride, log2_pb_w, log2_pb_h, mode_idx);
    }
    else{
        intra_angular_v_gauss_pdpc_avx2_4(ref_abv, ref_lft, dst, dst_stride, log2_pb_w, log2_pb_h, mode_idx);
    }
    
}

static void
intra_angular_h_gauss_pdpc_avx2(const OVSample* ref_abv, const OVSample* ref_lft,
                           OVSample* const dst, ptrdiff_t dst_stride,
                           int8_t log2_pb_w, int8_t log2_pb_h, int mode_idx)
{
    OVSample tmp_dst[128 * 128];
    const int tmp_stride = 128;
    OVSample* _tmp = tmp_dst;
    OVSample* _dst = dst;
    int pb_w = 1 << log2_pb_w;
    int pb_h = 1 << log2_pb_h;

    if (log2_pb_h >=4){
        intra_angular_v_gauss_pdpc_avx2_16(ref_lft, ref_abv, _tmp, tmp_stride, log2_pb_h, log2_pb_w, mode_idx);
    }
    else if (log2_pb_h ==3){
        intra_angular_v_gauss_pdpc_avx2_8(ref_lft, ref_abv, _tmp, tmp_stride, log2_pb_h, log2_pb_w, mode_idx);
    }
    else{
        intra_angular_v_gauss_pdpc_avx2_4(ref_lft, ref_abv, _tmp, tmp_stride, log2_pb_h, log2_pb_w, mode_idx);
    }

    for (int y = 0; y < pb_w; y++) {
        _dst = &dst[y];
        for (int x = 0; x < pb_h; x++) {
            _dst[0] = _tmp[x];
            _dst += dst_stride;
        }
        _tmp += tmp_stride;
    }
}

static void
intra_angular_v_cubic_pdpc_avx2_16(const OVSample* ref_abv, const OVSample* ref_lft,
                           OVSample* const dst, ptrdiff_t dst_stride,
                           int8_t log2_pb_w, int8_t log2_pb_h, int mode_idx)
{
    OVSample* _dst = dst;
    int angle_val = angle_table[mode_idx];
    int inv_angle = inverse_angle_table[mode_idx];
    int pb_w = 1 << log2_pb_w;
    int pb_h = 1 << log2_pb_h;
    int delta_pos = angle_val;
    int scale = OVMIN(2, log2_pb_h - (floor_log2(3 * inv_angle - 2) - 8));
    __m256i offset = _mm256_set1_epi32(32);
    for (int y = 0; y < pb_h; y++) {
        const int delta_int  = delta_pos >> 5;
        const int delta_frac = delta_pos & 0x1F;
        int inv_angle_sum = 256 + inv_angle;
        const OVSample* ref = (OVSample*)ref_abv + delta_int;
        const int8_t* filter = &chroma_filter[delta_frac << 2];
        LOAD_CUBIC_FILTER_AVX2();
        for (int x = 0; x < pb_w; x+=16) {
            FILTER_16_SAMPLES();
            ref+=16;
        }
        for (int x = 0; x < OVMIN(3 << scale, pb_w); x++) {
            int wL = 32 >> ((x << 1) >> scale);
            const OVSample* p = ref_lft + y + (inv_angle_sum >> 9) + 1;

            int16_t left = p[0];
            _dst[x] =
                ov_bdclip(_dst[x] + ((wL * (left - _dst[x]) + 32) >> 6));
            inv_angle_sum += inv_angle;
        }
        delta_pos += angle_val;
        _dst += dst_stride;
    }
}

static void
intra_angular_v_cubic_pdpc_avx2_8(const OVSample* ref_abv, const OVSample* ref_lft,
                           OVSample* const dst, ptrdiff_t dst_stride,
                           int8_t log2_pb_w, int8_t log2_pb_h, int mode_idx)
{
    OVSample* _dst = dst;
    int angle_val = angle_table[mode_idx];
    int inv_angle = inverse_angle_table[mode_idx];
    int pb_w = 1 << log2_pb_w;
    int pb_h = 1 << log2_pb_h;
    int delta_pos = angle_val;
    int scale = OVMIN(2, log2_pb_h - (floor_log2(3 * inv_angle - 2) - 8));
    __m128i offset = _mm_set1_epi32(32);
    for (int y = 0; y < pb_h; y++) {
        const int delta_int  = delta_pos >> 5;
        const int delta_frac = delta_pos & 0x1F;
        int inv_angle_sum = 256 + inv_angle;
        const OVSample* ref = (OVSample*)ref_abv + delta_int;
        const int8_t* filter = &chroma_filter[delta_frac << 2];
        LOAD_CUBIC_FILTER_SSE();
        for (int x = 0; x < pb_w; x+=8) {
            FILTER_8_SAMPLES();
            ref+=8;
        }
        for (int x = 0; x < OVMIN(3 << scale, pb_w); x++) {
            int wL = 32 >> ((x << 1) >> scale);
            const OVSample* p = ref_lft + y + (inv_angle_sum >> 9) + 1;

            int16_t left = p[0];
            _dst[x] =
                ov_bdclip(_dst[x] + ((wL * (left - _dst[x]) + 32) >> 6));
            inv_angle_sum += inv_angle;
        }
        delta_pos += angle_val;
        _dst += dst_stride;
    }
}

static void
intra_angular_v_cubic_pdpc_avx2_4(const OVSample* ref_abv, const OVSample* ref_lft,
                           OVSample* const dst, ptrdiff_t dst_stride,
                           int8_t log2_pb_w, int8_t log2_pb_h, int mode_idx)
{
    OVSample* _dst = dst;
    int angle_val = angle_table[mode_idx];
    int inv_angle = inverse_angle_table[mode_idx];
    int pb_w = 1 << log2_pb_w;
    int pb_h = 1 << log2_pb_h;
    int delta_pos = angle_val;
    int scale = OVMIN(2, log2_pb_h - (floor_log2(3 * inv_angle - 2) - 8));
    __m128i offset = _mm_set1_epi32(32);
    for (int y = 0; y < pb_h; y++) {
        const int delta_int  = delta_pos >> 5;
        const int delta_frac = delta_pos & 0x1F;
        int inv_angle_sum = 256 + inv_angle;
        const OVSample* ref = (OVSample*)ref_abv + delta_int;
        const int8_t* filter = &chroma_filter[delta_frac << 2];
        LOAD_CUBIC_FILTER_SSE();
        FILTER_4_SAMPLES();
        for (int x = 0; x < OVMIN(3 << scale, pb_w); x++) {
            int wL = 32 >> ((x << 1) >> scale);
            const OVSample* p = ref_lft + y + (inv_angle_sum >> 9) + 1;

            int16_t left = p[0];
            _dst[x] =
                ov_bdclip(_dst[x] + ((wL * (left - _dst[x]) + 32) >> 6));
            inv_angle_sum += inv_angle;
        }
        delta_pos += angle_val;
        _dst += dst_stride;
    }
}


static void
intra_angular_v_cubic_pdpc_avx2(const OVSample* ref_abv, const OVSample* ref_lft,
                           OVSample* const dst, ptrdiff_t dst_stride,
                           int8_t log2_pb_w, int8_t log2_pb_h, int mode_idx)
{
    if (log2_pb_w >=4){
        intra_angular_v_cubic_pdpc_avx2_16(ref_abv, ref_lft, dst, dst_stride, log2_pb_w, log2_pb_h, mode_idx);
    }
    else if (log2_pb_w ==3){
        intra_angular_v_cubic_pdpc_avx2_8(ref_abv, ref_lft, dst, dst_stride, log2_pb_w, log2_pb_h, mode_idx);
    }
    else{
        intra_angular_v_cubic_pdpc_avx2_4(ref_abv, ref_lft, dst, dst_stride, log2_pb_w, log2_pb_h, mode_idx);
    }
    
}

static void
intra_angular_h_cubic_pdpc_avx2(const OVSample* ref_abv, const OVSample* ref_lft,
                           OVSample* const dst, ptrdiff_t dst_stride,
                           int8_t log2_pb_w, int8_t log2_pb_h, int mode_idx)
{
    OVSample tmp_dst[128 * 128];
    const int tmp_stride = 128;
    OVSample* _tmp = tmp_dst;
    OVSample* _dst = dst;
    int pb_w = 1 << log2_pb_w;
    int pb_h = 1 << log2_pb_h;

    if (log2_pb_h >=4){
        intra_angular_v_cubic_pdpc_avx2_16(ref_lft, ref_abv, _tmp, tmp_stride, log2_pb_h, log2_pb_w, mode_idx);
    }
    else if (log2_pb_h ==3){
        intra_angular_v_cubic_pdpc_avx2_8(ref_lft, ref_abv, _tmp, tmp_stride, log2_pb_h, log2_pb_w, mode_idx);
    }
    else{
        intra_angular_v_cubic_pdpc_avx2_4(ref_lft, ref_abv, _tmp, tmp_stride, log2_pb_h, log2_pb_w, mode_idx);
    }

    for (int y = 0; y < pb_w; y++) {
        _dst = &dst[y];
        for (int x = 0; x < pb_h; x++) {
            _dst[0] = _tmp[x];
            _dst += dst_stride;
        }
        _tmp += tmp_stride;
    }
}

static void
intra_angular_v_c_pdpc_avx2_16(const OVSample* const ref_abv,
                       const OVSample* const ref_lft, OVSample* const dst,
                       ptrdiff_t dst_stride, int8_t log2_pb_w,
                       int8_t log2_pb_h, int mode_idx)
{
    OVSample* _dst = dst;
    int angle_val = angle_table[mode_idx];
    int inv_angle = inverse_angle_table[mode_idx];
    int delta_pos = angle_val;
    int pb_w = 1 << log2_pb_w;
    int pb_h = 1 << log2_pb_h;
    int scale = OVMIN(2, log2_pb_h - (floor_log2(3 * inv_angle - 2) - 8));
    __m256i offset = _mm256_set1_epi32(16);

    for (int y = 0, delta_pos = angle_val; y < pb_h; y++) {
        const int delta_int = delta_pos >> 5;
        int inv_angle_sum = 256 + inv_angle;
        __m256i delta_frac = _mm256_set1_epi16(delta_pos & 0x1F);
        const OVSample* ref = ref_abv + delta_int + 1;
        for (int x = 0; x < pb_w; ref+=16, x+=16) {
            __m256i last_ref_val = _mm256_loadu_si256((__m256i *)&ref[0]);
            __m256i curr_ref_val = _mm256_loadu_si256((__m256i *)&ref[1]);
            __m256i val = _mm256_sub_epi16(curr_ref_val, last_ref_val);
            
            __m256i valL = _mm256_mullo_epi16(delta_frac, val);
            __m256i valH = _mm256_mulhi_epi16(delta_frac, val);

            __m256i val03 = _mm256_unpacklo_epi16(valL, valH);
            __m256i val47 = _mm256_unpackhi_epi16(valL, valH);

            val03 = _mm256_add_epi32(val03, offset);
            val47 = _mm256_add_epi32(val47, offset);

            val03 = _mm256_srai_epi32(val03, 5);
            val47 = _mm256_srai_epi32(val47, 5);

            // FIXME PACKS before clipping might cause errors
            val = _mm256_packs_epi32(val03, val47);

            val = _mm256_add_epi16(val, last_ref_val);

            _mm256_storeu_si256((__m256i *) &_dst[x], val);
        }
        for (int x = 0; x < OVMIN(3 << scale, pb_w); x++) {
            int wL = 32 >> ((x << 1) >> scale);
            const OVSample* p = ref_lft + y + (inv_angle_sum >> 9) + 1;

            int32_t left = p[0];
            _dst[x] = ov_bdclip(_dst[x] + ((wL * (left - _dst[x]) + 32) >> 6));
            inv_angle_sum += inv_angle;
        }
        delta_pos += angle_val;
        _dst += dst_stride;
    }
}

static void
intra_angular_v_c_pdpc_avx2_8(const OVSample* const ref_abv,
                       const OVSample* const ref_lft, OVSample* const dst,
                       ptrdiff_t dst_stride, int8_t log2_pb_w,
                       int8_t log2_pb_h, int mode_idx)
{
    OVSample* _dst = dst;
    int angle_val = angle_table[mode_idx];
    int inv_angle = inverse_angle_table[mode_idx];
    int delta_pos = angle_val;
    int pb_w = 1 << log2_pb_w;
    int pb_h = 1 << log2_pb_h;
    int scale = OVMIN(2, log2_pb_h - (floor_log2(3 * inv_angle - 2) - 8));
    __m128i offset = _mm_set1_epi32(16);

    for (int y = 0, delta_pos = angle_val; y < pb_h; y++) {
        const int delta_int = delta_pos >> 5;
        int inv_angle_sum = 256 + inv_angle;
        __m128i delta_frac = _mm_set1_epi16(delta_pos & 0x1F);
        const OVSample* ref = ref_abv + delta_int + 1;
        for (int x = 0; x < pb_w; ref+=8, x+=8) {
            __m128i last_ref_val = _mm_loadu_si128((__m128i *)&ref[0]);
            __m128i curr_ref_val = _mm_loadu_si128((__m128i *)&ref[1]);
            __m128i val = _mm_sub_epi16(curr_ref_val, last_ref_val);
            
            __m128i valL = _mm_mullo_epi16(delta_frac, val);
            __m128i valH = _mm_mulhi_epi16(delta_frac, val);

            __m128i val03 = _mm_unpacklo_epi16(valL, valH);
            __m128i val47 = _mm_unpackhi_epi16(valL, valH);

            val03 = _mm_add_epi32(val03, offset);
            val47 = _mm_add_epi32(val47, offset);

            val03 = _mm_srai_epi32(val03, 5);
            val47 = _mm_srai_epi32(val47, 5);

            // FIXME PACKS before clipping might cause errors
            val = _mm_packs_epi32(val03, val47);

            val = _mm_add_epi16(val, last_ref_val);

            _mm_storeu_si128((__m128i *) &_dst[x], val);
        }
        for (int x = 0; x < OVMIN(3 << scale, pb_w); x++) {
            int wL = 32 >> ((x << 1) >> scale);
            const OVSample* p = ref_lft + y + (inv_angle_sum >> 9) + 1;

            int32_t left = p[0];
            _dst[x] = ov_bdclip(_dst[x] + ((wL * (left - _dst[x]) + 32) >> 6));
            inv_angle_sum += inv_angle;
        }
        delta_pos += angle_val;
        _dst += dst_stride;
    }
}

static void
intra_angular_v_c_pdpc_avx2_4(const OVSample* const ref_abv,
                       const OVSample* const ref_lft, OVSample* const dst,
                       ptrdiff_t dst_stride, int8_t log2_pb_w,
                       int8_t log2_pb_h, int mode_idx)
{
    OVSample* _dst = dst;
    int angle_val = angle_table[mode_idx];
    int inv_angle = inverse_angle_table[mode_idx];
    int delta_pos = angle_val;
    int pb_w = 1 << log2_pb_w;
    int pb_h = 1 << log2_pb_h;
    int scale = OVMIN(2, log2_pb_h - (floor_log2(3 * inv_angle - 2) - 8));
    __m128i offset = _mm_set1_epi32(16);

    for (int y = 0, delta_pos = angle_val; y < pb_h; y++) {
        const int delta_int = delta_pos >> 5;
        int inv_angle_sum = 256 + inv_angle;
        __m128i delta_frac = _mm_set1_epi16(delta_pos & 0x1F);
        const OVSample* ref = ref_abv + delta_int + 1;
        __m128i last_ref_val = _mm_loadl_epi64((__m128i *)&ref[0]);
        __m128i curr_ref_val = _mm_loadl_epi64((__m128i *)&ref[1]);
        __m128i val = _mm_sub_epi16(curr_ref_val, last_ref_val);
        
        __m128i valL = _mm_mullo_epi16(delta_frac, val);
        __m128i valH = _mm_mulhi_epi16(delta_frac, val);

        __m128i val03 = _mm_unpacklo_epi16(valL, valH);

        val03 = _mm_add_epi32(val03, offset);

        val03 = _mm_srai_epi32(val03, 5);

        // FIXME PACKS before clipping might cause errors
        val = _mm_packs_epi32(val03, _mm_setzero_si128());

        val = _mm_add_epi16(val, last_ref_val);

        _mm_storel_epi64((__m128i *) &_dst[0], val);

        for (int x = 0; x < OVMIN(3 << scale, pb_w); x++) {
            int wL = 32 >> ((x << 1) >> scale);
            const OVSample* p = ref_lft + y + (inv_angle_sum >> 9) + 1;

            int32_t left = p[0];
            _dst[x] = ov_bdclip(_dst[x] + ((wL * (left - _dst[x]) + 32) >> 6));
            inv_angle_sum += inv_angle;
        }
        delta_pos += angle_val;
        _dst += dst_stride;
    }
}

static void
intra_angular_v_c_pdpc_avx2(const OVSample* const ref_abv,
                       const OVSample* const ref_lft, OVSample* const dst,
                       ptrdiff_t dst_stride, int8_t log2_pb_w,
                       int8_t log2_pb_h, int mode_idx)
{
    OVSample* _dst = dst;
    int pb_w = 1 << log2_pb_w;
    int pb_h = 1 << log2_pb_h;

    if (pb_w >=16){
        intra_angular_v_c_pdpc_avx2_16(ref_abv, ref_lft, dst, dst_stride, log2_pb_w, log2_pb_h, mode_idx);
    } else if (pb_w ==8){
        intra_angular_v_c_pdpc_avx2_8(ref_abv, ref_lft, dst, dst_stride, log2_pb_w, log2_pb_h, mode_idx);
    } else {
        intra_angular_v_c_pdpc_avx2_4(ref_abv, ref_lft, dst, dst_stride, log2_pb_w, log2_pb_h, mode_idx);
    }
}

static void
intra_angular_h_c_pdpc_avx2(const OVSample* const ref_abv,
                       const OVSample* const ref_lft, OVSample* const dst,
                       ptrdiff_t dst_stride, int8_t log2_pb_w,
                       int8_t log2_pb_h, int mode_idx)
{
    OVSample tmp_dst[128 * 128];
    OVSample* _tmp = tmp_dst;
    OVSample* _dst = dst;
    int pb_w = 1 << log2_pb_w;
    int pb_h = 1 << log2_pb_h;
    const int tmp_stride = 128;

    if (pb_h >=16){
        intra_angular_v_c_pdpc_avx2_16(ref_lft, ref_abv, tmp_dst, tmp_stride, log2_pb_h, log2_pb_w, mode_idx);
    } else if (pb_h ==8){
        intra_angular_v_c_pdpc_avx2_8(ref_lft, ref_abv, tmp_dst, tmp_stride, log2_pb_h, log2_pb_w, mode_idx);
    } else {
        intra_angular_v_c_pdpc_avx2_4(ref_lft, ref_abv, tmp_dst, tmp_stride, log2_pb_h, log2_pb_w, mode_idx);
    }

    for (int y = 0; y < pb_w; y++) {
        _dst = &dst[y];
        for (int x = 0; x < pb_h; x++) {
            _dst[0] = _tmp[x];
            _dst += dst_stride;
        }
        _tmp += tmp_stride;
    }

}

static void
intra_angular_v_cubic_mref_avx2_16(const OVSample* const ref_abv, OVSample* const dst,
                           ptrdiff_t dst_stride, int8_t log2_pb_w,
                           int8_t log2_pb_h, int angle_val,
                           uint8_t mrl_idx)
{
    int delta_pos = angle_val * (mrl_idx + 1);
    OVSample* _dst = dst;
    int pb_w = 1 << log2_pb_w;
    int pb_h = 1 << log2_pb_h;
    __m256i offset = _mm256_set1_epi32(32);
    for (int y = 0; y < pb_h; y++) {
        const int delta_int  = delta_pos >> 5;
        const int delta_frac = delta_pos & 0x1F;
        const OVSample* ref = (OVSample *)ref_abv + delta_int;
        const int8_t* filter = &chroma_filter[delta_frac << 2];
        LOAD_CUBIC_FILTER_AVX2();
        
        for (int x = 0; x < pb_w; x+=16) {
            FILTER_16_SAMPLES();
            ref+=16;
        }
        delta_pos += angle_val;
        _dst += dst_stride;
    }
}

static void
intra_angular_v_cubic_mref_avx2_8(const OVSample* const ref_abv, OVSample* const dst,
                           ptrdiff_t dst_stride, int8_t log2_pb_w,
                           int8_t log2_pb_h, int angle_val,
                           uint8_t mrl_idx)
{
    int delta_pos = angle_val * (mrl_idx + 1);
    OVSample* _dst = dst;
    int pb_w = 1 << log2_pb_w;
    int pb_h = 1 << log2_pb_h;
    __m128i offset = _mm_set1_epi32(32);
    for (int y = 0; y < pb_h; y++) {
        const int delta_int  = delta_pos >> 5;
        const int delta_frac = delta_pos & 0x1F;
        const OVSample* ref = (OVSample *)ref_abv + delta_int;
        const int8_t* filter = &chroma_filter[delta_frac << 2];
        LOAD_CUBIC_FILTER_SSE();
        
        for (int x = 0; x < pb_w; x+=8) {
            FILTER_8_SAMPLES();
            ref+=8;
        }
        delta_pos += angle_val;
        _dst += dst_stride;
    }
}

static void
intra_angular_v_cubic_mref_avx2_4(const OVSample* const ref_abv, OVSample* const dst,
                           ptrdiff_t dst_stride, int8_t log2_pb_w,
                           int8_t log2_pb_h, int angle_val,
                           uint8_t mrl_idx)
{
    int delta_pos = angle_val * (mrl_idx + 1);
    OVSample* _dst = dst;
    int pb_w = 1 << log2_pb_w;
    int pb_h = 1 << log2_pb_h;
    __m128i offset = _mm_set1_epi32(32);
    for (int y = 0; y < pb_h; y++) {
        const int delta_int  = delta_pos >> 5;
        const int delta_frac = delta_pos & 0x1F;
        const OVSample* ref = (OVSample *)ref_abv + delta_int;
        const int8_t* filter = &chroma_filter[delta_frac << 2];
        LOAD_CUBIC_FILTER_SSE();

        FILTER_4_SAMPLES();

        delta_pos += angle_val;
        _dst += dst_stride;
    }
}

static void
intra_angular_h_cubic_mref_avx2(const OVSample* const ref_lft, OVSample* const dst,
                           ptrdiff_t dst_stride,
                           int8_t log2_pb_w, int8_t log2_pb_h,
                           int angle_val, uint8_t mrl_idx)
{
    OVSample tmp_dst[128 * 128];
    const int tmp_stride = 128;
    OVSample* _tmp = tmp_dst;
    OVSample* _dst = dst;
    int pb_w = 1 << log2_pb_w;
    int pb_h = 1 << log2_pb_h;

    if (log2_pb_h >=4){
        intra_angular_v_cubic_mref_avx2_16(ref_lft, _tmp, tmp_stride, log2_pb_h, log2_pb_w,angle_val, mrl_idx);
    }
    else if (log2_pb_h ==3){
        intra_angular_v_cubic_mref_avx2_8(ref_lft, _tmp, tmp_stride, log2_pb_h, log2_pb_w,angle_val, mrl_idx);
    }
    else{
        intra_angular_v_cubic_mref_avx2_4(ref_lft, _tmp, tmp_stride, log2_pb_h, log2_pb_w,angle_val, mrl_idx);
    }

    for (int y = 0; y < pb_w; y++) {
        _dst = &dst[y];
        for (int x = 0; x < pb_h; x++) {
            _dst[0] = _tmp[x];
            _dst += dst_stride;
        }
        _tmp += tmp_stride;
    }
}

static void
intra_angular_v_cubic_mref_avx2(const OVSample* const ref_abv, OVSample* const dst,
                           ptrdiff_t dst_stride, int8_t log2_pb_w,
                           int8_t log2_pb_h, int angle_val,
                           uint8_t mrl_idx)
{
    if (log2_pb_w >=4){
        intra_angular_v_cubic_mref_avx2_16(ref_abv, dst, dst_stride, log2_pb_w, log2_pb_h,angle_val, mrl_idx);
    }
    else if (log2_pb_w ==3){
        intra_angular_v_cubic_mref_avx2_8(ref_abv, dst, dst_stride, log2_pb_w, log2_pb_h,angle_val, mrl_idx);
    }
    else{
        intra_angular_v_cubic_mref_avx2_4(ref_abv, dst, dst_stride, log2_pb_w, log2_pb_h,angle_val, mrl_idx);
    }
    
}

struct IntraAngularFunctions angular_gauss_h_10_avx2;
struct IntraAngularFunctions angular_gauss_v_10_avx2;
struct IntraAngularFunctions angular_cubic_v_10_avx2;
struct IntraAngularFunctions angular_cubic_h_10_avx2;
struct IntraAngularFunctions angular_c_h_10_avx2;
struct IntraAngularFunctions angular_c_v_10_avx2;
struct IntraAngularFunctions angular_nofrac_v_10_avx2;
struct IntraAngularFunctions angular_nofrac_h_10_avx2;
struct IntraMRLFunctions mrl_func_10_avx2;

void
rcn_init_intra_angular_functions_10_avx2(struct RCNFunctions *rcn_func)
{
    angular_gauss_h_10_avx2.pure = rcn_func->intra_angular_gauss_h->pure;
    angular_gauss_h_10_avx2.diagonal = rcn_func->intra_angular_gauss_h->diagonal;
    angular_gauss_h_10_avx2.angular = intra_angular_h_gauss_avx2;

    angular_gauss_h_10_avx2.pure_pdpc = rcn_func->intra_angular_gauss_h->pure_pdpc;
    angular_gauss_h_10_avx2.diagonal_pdpc = rcn_func->intra_angular_gauss_h->diagonal_pdpc;
    angular_gauss_h_10_avx2.angular_pdpc = intra_angular_h_gauss_pdpc_avx2;


    angular_gauss_v_10_avx2.pure = rcn_func->intra_angular_gauss_v->pure;
    angular_gauss_v_10_avx2.diagonal = rcn_func->intra_angular_gauss_v->diagonal;
    angular_gauss_v_10_avx2.angular = intra_angular_v_gauss_avx2;

    angular_gauss_v_10_avx2.pure_pdpc = rcn_func->intra_angular_gauss_v->pure_pdpc;
    angular_gauss_v_10_avx2.diagonal_pdpc = rcn_func->intra_angular_gauss_v->diagonal_pdpc;
    angular_gauss_v_10_avx2.angular_pdpc = intra_angular_v_gauss_pdpc_avx2;


    angular_cubic_v_10_avx2.pure = rcn_func->intra_angular_cubic_v->pure;
    angular_cubic_v_10_avx2.diagonal = rcn_func->intra_angular_cubic_v->diagonal;
    angular_cubic_v_10_avx2.angular = intra_angular_v_cubic_avx2;

    angular_cubic_v_10_avx2.pure_pdpc = rcn_func->intra_angular_cubic_v->pure_pdpc;
    angular_cubic_v_10_avx2.diagonal_pdpc = rcn_func->intra_angular_cubic_v->diagonal_pdpc;
    angular_cubic_v_10_avx2.angular_pdpc = intra_angular_v_cubic_pdpc_avx2;


    angular_cubic_h_10_avx2.pure = rcn_func->intra_angular_cubic_h->pure;
    angular_cubic_h_10_avx2.diagonal = rcn_func->intra_angular_cubic_h->diagonal;
    angular_cubic_h_10_avx2.angular = intra_angular_h_cubic_avx2;

    angular_cubic_h_10_avx2.pure_pdpc = rcn_func->intra_angular_cubic_h->pure_pdpc;
    angular_cubic_h_10_avx2.diagonal_pdpc = rcn_func->intra_angular_cubic_h->diagonal_pdpc;
    angular_cubic_h_10_avx2.angular_pdpc = intra_angular_h_cubic_pdpc_avx2;


    angular_c_h_10_avx2.pure = rcn_func->intra_angular_c_h->pure;
    angular_c_h_10_avx2.diagonal = rcn_func->intra_angular_c_h->diagonal;
    angular_c_h_10_avx2.angular = intra_angular_h_c_avx2;

    angular_c_h_10_avx2.pure_pdpc = rcn_func->intra_angular_c_h->pure_pdpc;
    angular_c_h_10_avx2.diagonal_pdpc = rcn_func->intra_angular_c_h->diagonal_pdpc;
    angular_c_h_10_avx2.angular_pdpc = intra_angular_h_c_pdpc_avx2;


    angular_c_v_10_avx2.pure = rcn_func->intra_angular_c_v->pure;
    angular_c_v_10_avx2.diagonal = rcn_func->intra_angular_c_v->diagonal;
    angular_c_v_10_avx2.angular = intra_angular_v_c_avx2;

    angular_c_v_10_avx2.pure_pdpc = rcn_func->intra_angular_c_v->pure_pdpc;
    angular_c_v_10_avx2.diagonal_pdpc = rcn_func->intra_angular_c_v->diagonal_pdpc;
    angular_c_v_10_avx2.angular_pdpc = intra_angular_v_c_pdpc_avx2;


    angular_nofrac_v_10_avx2.pure = rcn_func->intra_angular_nofrac_v->pure;
    angular_nofrac_v_10_avx2.diagonal = rcn_func->intra_angular_nofrac_v->diagonal;
    angular_nofrac_v_10_avx2.angular = rcn_func->intra_angular_nofrac_v->angular;

    angular_nofrac_v_10_avx2.pure_pdpc = rcn_func->intra_angular_nofrac_v->pure_pdpc;
    angular_nofrac_v_10_avx2.diagonal_pdpc = rcn_func->intra_angular_nofrac_v->diagonal_pdpc;
    angular_nofrac_v_10_avx2.angular_pdpc = rcn_func->intra_angular_nofrac_v->angular_pdpc;


    angular_nofrac_h_10_avx2.pure = rcn_func->intra_angular_nofrac_h->pure;
    angular_nofrac_h_10_avx2.diagonal = rcn_func->intra_angular_nofrac_h->diagonal;
    angular_nofrac_h_10_avx2.angular = rcn_func->intra_angular_nofrac_h->angular;

    angular_nofrac_h_10_avx2.pure_pdpc = rcn_func->intra_angular_nofrac_h->pure_pdpc;
    angular_nofrac_h_10_avx2.diagonal_pdpc = rcn_func->intra_angular_nofrac_h->diagonal_pdpc;
    angular_nofrac_h_10_avx2.angular_pdpc = rcn_func->intra_angular_nofrac_h->angular_pdpc;


    mrl_func_10_avx2.angular_h = intra_angular_h_cubic_mref_avx2;
    mrl_func_10_avx2.angular_v = intra_angular_v_cubic_mref_avx2;


    rcn_func->intra_angular_gauss_h = &angular_gauss_h_10_avx2;
    rcn_func->intra_angular_gauss_v = &angular_gauss_v_10_avx2;
    rcn_func->intra_angular_cubic_h = &angular_cubic_h_10_avx2;
    rcn_func->intra_angular_cubic_v = &angular_cubic_v_10_avx2;
    rcn_func->intra_angular_c_h     = &angular_c_h_10_avx2;
    rcn_func->intra_angular_c_v     = &angular_c_v_10_avx2;
    rcn_func->intra_angular_nofrac_v = &angular_nofrac_v_10_avx2;
    rcn_func->intra_angular_nofrac_h = &angular_nofrac_h_10_avx2;
    rcn_func->intra_mrl = &mrl_func_10_avx2;
}