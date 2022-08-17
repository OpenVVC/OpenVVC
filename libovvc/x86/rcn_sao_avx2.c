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

#include <stddef.h>
#include <stdint.h>

#include <immintrin.h>

#include "dec_structures.h"
#include "rcn_structures.h"


static void
sao_band_filter_0_10_avx2(OVSample* _dst,
                         OVSample* _src,
                         ptrdiff_t _stride_dst,
                         ptrdiff_t _stride_src,
                         struct SAOParamsCtu* sao,
                         int width,
                         int height,
                         int c_idx)
{
  int y, x;
  int shift = 10 - 5;
  int8_t* sao_offset_val = sao->offset_val[c_idx];
  uint8_t sao_left_class = sao->band_position[c_idx];
  __m256i r0, r1, r2, r3, sao1, sao2, sao3, sao4;
  uint16_t* dst = (uint16_t*)_dst;
  uint16_t* src = (uint16_t*)_src;
  ptrdiff_t stride_dst = _stride_dst;
  ptrdiff_t stride_src = _stride_src;
  r0 = _mm256_set1_epi16((sao_left_class)&31);
  r1 = _mm256_set1_epi16((sao_left_class + 1) & 31);
  r2 = _mm256_set1_epi16((sao_left_class + 2) & 31);
  r3 = _mm256_set1_epi16((sao_left_class + 3) & 31);
  sao1 = _mm256_set1_epi16(sao_offset_val[0]);
  sao2 = _mm256_set1_epi16(sao_offset_val[1]);
  sao3 = _mm256_set1_epi16(sao_offset_val[2]);
  sao4 = _mm256_set1_epi16(sao_offset_val[3]);
  for (y = 0; y < height; y++) {
    for (x = 0; x < width-width%16; x += 16) {
      __m256i src0 = _mm256_loadu_si256((__m256i*)&src[x]);
      __m256i src2 = _mm256_srai_epi16(src0, shift);
      __m256i x0 = _mm256_cmpeq_epi16(src2, r0);
      __m256i x1 = _mm256_cmpeq_epi16(src2, r1);
      __m256i x2 = _mm256_cmpeq_epi16(src2, r2);
      __m256i x3 = _mm256_cmpeq_epi16(src2, r3);
      x0 = _mm256_and_si256(x0, sao1);
      x1 = _mm256_and_si256(x1, sao2);
      x2 = _mm256_and_si256(x2, sao3);
      x3 = _mm256_and_si256(x3, sao4);
      x0 = _mm256_or_si256(x0, x1);
      x2 = _mm256_or_si256(x2, x3);
      x0 = _mm256_or_si256(x0, x2);
      src0 = _mm256_add_epi16(src0, x0);
      src0 = _mm256_max_epi16(src0, _mm256_setzero_si256());
      src0 = _mm256_min_epi16(src0, _mm256_set1_epi16(0x03FF));
      _mm256_storeu_si256((__m256i*)&dst[x], src0);
    }
    for ( ; x < width; x += 8) {
      __m128i src0 = _mm_loadu_si128((__m128i*)&src[x]);
      __m128i src2 = _mm_srai_epi16(src0, shift);
      __m128i x0 = _mm_cmpeq_epi16(src2, _mm256_castsi256_si128(r0));
      __m128i x1 = _mm_cmpeq_epi16(src2, _mm256_castsi256_si128(r1));
      __m128i x2 = _mm_cmpeq_epi16(src2, _mm256_castsi256_si128(r2));
      __m128i x3 = _mm_cmpeq_epi16(src2, _mm256_castsi256_si128(r3));
      x0 = _mm_and_si128(x0, _mm256_castsi256_si128(sao1));
      x1 = _mm_and_si128(x1, _mm256_castsi256_si128(sao2));
      x2 = _mm_and_si128(x2, _mm256_castsi256_si128(sao3));
      x3 = _mm_and_si128(x3, _mm256_castsi256_si128(sao4));
      x0 = _mm_or_si128(x0, x1);
      x2 = _mm_or_si128(x2, x3);
      x0 = _mm_or_si128(x0, x2);
      src0 = _mm_add_epi16(src0, x0);
      src0 = _mm_max_epi16(src0, _mm_setzero_si128());
      src0 = _mm_min_epi16(src0, _mm_set1_epi16(0x03FF));
      _mm_storeu_si128((__m128i*)&dst[x], src0);
    }
    dst += stride_dst;
    src += stride_src;
  }
}

static void
sao_edge_filter_10_avx2(OVSample* _dst,
                       OVSample* _src,
                       ptrdiff_t _stride_dst,
                       ptrdiff_t _stride_src,
                       SAOParamsCtu* sao,
                       int width,
                       int height,
                       int c_idx)
{
  int x, y;
  int8_t* sao_offset_val = sao->offset_val[c_idx];
  int eo = sao->eo_class[c_idx];
  static const int8_t pos[4][2][2] = {
    { { -1, 0 }, { 1, 0 } },
    { { 0, -1 }, { 0, 1 } },
    { { -1, -1 }, { 1, 1 } },
    { { 1, -1 }, { -1, 1 } },
  };
  __m256i offset0, offset1, offset2, offset3, offset4;
  uint16_t* dst = (uint16_t*)_dst;
  uint16_t* src = (uint16_t*)_src;
  ptrdiff_t stride_dst = _stride_dst;
  ptrdiff_t stride_src = _stride_src;
  {
    int a_stride = pos[eo][0][0] + pos[eo][0][1] * stride_src;
    int b_stride = pos[eo][1][0] + pos[eo][1][1] * stride_src;
    offset0 = _mm256_set1_epi16(sao_offset_val[0]);
    offset1 = _mm256_set1_epi16(sao_offset_val[1]);
    offset2 = _mm256_set1_epi16(0);
    offset3 = _mm256_set1_epi16(sao_offset_val[2]);
    offset4 = _mm256_set1_epi16(sao_offset_val[3]);
    for (y = 0; y < height; y++) {
      for (x = 0; x < width-width%16; x += 16) {
        __m256i x0, x1, x2, x3;
        __m256i cmp0, cmp1, r0, r1, r2, r3, r4;
        x0 = _mm256_loadu_si256((__m256i*)(src + x));
        cmp0 = _mm256_loadu_si256((__m256i*)(src + x + a_stride));
        cmp1 = _mm256_loadu_si256((__m256i*)(src + x + b_stride));
        r2 = _mm256_min_epu16(x0, cmp0);
        x1 = _mm256_cmpeq_epi16(cmp0, r2);
        x2 = _mm256_cmpeq_epi16(x0, r2);
        x1 = _mm256_sub_epi16(x2, x1);
        r2 = _mm256_min_epu16(x0, cmp1);
        x3 = _mm256_cmpeq_epi16(cmp1, r2);
        x2 = _mm256_cmpeq_epi16(x0, r2);
        x3 = _mm256_sub_epi16(x2, x3);
        x1 = _mm256_add_epi16(x1, x3);
        r0 = _mm256_cmpeq_epi16(x1, _mm256_set1_epi16(-2));
        r1 = _mm256_cmpeq_epi16(x1, _mm256_set1_epi16(-1));
        r2 = _mm256_cmpeq_epi16(x1, _mm256_set1_epi16(0));
        r3 = _mm256_cmpeq_epi16(x1, _mm256_set1_epi16(1));
        r4 = _mm256_cmpeq_epi16(x1, _mm256_set1_epi16(2));
        r0 = _mm256_and_si256(r0, offset0);
        r1 = _mm256_and_si256(r1, offset1);
        r2 = _mm256_and_si256(r2, offset2);
        r3 = _mm256_and_si256(r3, offset3);
        r4 = _mm256_and_si256(r4, offset4);
        r0 = _mm256_add_epi16(r0, r1);
        r2 = _mm256_add_epi16(r2, r3);
        r0 = _mm256_add_epi16(r0, r4);
        r0 = _mm256_add_epi16(r0, r2);
        r0 = _mm256_add_epi16(r0, x0);
        r1 = _mm256_set1_epi16(0x03FF);
        r0 = _mm256_max_epi16(r0, _mm256_setzero_si256());
        r0 = _mm256_min_epi16(r0, r1);
        _mm256_storeu_si256((__m256i*)(dst + x), r0);
      }
      for (; x < width; x += 8) {
        __m128i x0, x1, x2, x3;
        __m128i cmp0, cmp1, r0, r1, r2, r3, r4;
        x0 = _mm_loadu_si128((__m128i*)(src + x));
        cmp0 = _mm_loadu_si128((__m128i*)(src + x + a_stride));
        cmp1 = _mm_loadu_si128((__m128i*)(src + x + b_stride));
        r2 = _mm_min_epu16(x0, cmp0);
        x1 = _mm_cmpeq_epi16(cmp0, r2);
        x2 = _mm_cmpeq_epi16(x0, r2);
        x1 = _mm_sub_epi16(x2, x1);
        r2 = _mm_min_epu16(x0, cmp1);
        x3 = _mm_cmpeq_epi16(cmp1, r2);
        x2 = _mm_cmpeq_epi16(x0, r2);
        x3 = _mm_sub_epi16(x2, x3);
        x1 = _mm_add_epi16(x1, x3);
        r0 = _mm_cmpeq_epi16(x1, _mm_set1_epi16(-2));
        r1 = _mm_cmpeq_epi16(x1, _mm_set1_epi16(-1));
        r2 = _mm_cmpeq_epi16(x1, _mm_set1_epi16(0));
        r3 = _mm_cmpeq_epi16(x1, _mm_set1_epi16(1));
        r4 = _mm_cmpeq_epi16(x1, _mm_set1_epi16(2));
        r0 = _mm_and_si128(r0, _mm256_castsi256_si128(offset0));
        r1 = _mm_and_si128(r1, _mm256_castsi256_si128(offset1));
        r2 = _mm_and_si128(r2, _mm256_castsi256_si128(offset2));
        r3 = _mm_and_si128(r3, _mm256_castsi256_si128(offset3));
        r4 = _mm_and_si128(r4, _mm256_castsi256_si128(offset4));
        r0 = _mm_add_epi16(r0, r1);
        r2 = _mm_add_epi16(r2, r3);
        r0 = _mm_add_epi16(r0, r4);
        r0 = _mm_add_epi16(r0, r2);
        r0 = _mm_add_epi16(r0, x0);
        r1 = _mm_set1_epi16(0x03FF);
        r0 = _mm_max_epi16(r0, _mm_setzero_si128());
        r0 = _mm_min_epi16(r0, r1);
        _mm_storeu_si128((__m128i*)(dst + x), r0);
      }
      src += stride_src;
      dst += stride_dst;
    }
  }
}

static void
sao_edge_filter_7_10_avx2(OVSample* _dst,
                         OVSample* _src,
                         ptrdiff_t _stride_dst,
                         ptrdiff_t _stride_src,
                         SAOParamsCtu* sao,
                         int width,
                         int height,
                         int c_idx)
{
  int x, y;
  int8_t* sao_offset_val = sao->offset_val[c_idx];
  int eo = sao->eo_class[c_idx];
  const int8_t pos[4][2][2] = {
    { { -1, 0 }, { 1, 0 } },
    { { 0, -1 }, { 0, 1 } },
    { { -1, -1 }, { 1, 1 } },
    { { 1, -1 }, { -1, 1 } },
  };
  __m256i offset0, offset1, offset2, offset3, offset4;
  uint16_t* dst = (uint16_t*)_dst;
  uint16_t* src = (uint16_t*)_src;
  ptrdiff_t stride_dst = _stride_dst;
  ptrdiff_t stride_src = _stride_src;
  {
    int a_stride = pos[eo][0][0] + pos[eo][0][1] * stride_src;
    int b_stride = pos[eo][1][0] + pos[eo][1][1] * stride_src;
    offset0 = _mm256_set1_epi16(sao_offset_val[0]);
    offset1 = _mm256_set1_epi16(sao_offset_val[1]);
    offset2 = _mm256_set1_epi16(0);
    offset3 = _mm256_set1_epi16(sao_offset_val[2]);
    offset4 = _mm256_set1_epi16(sao_offset_val[3]);
    for (y = 0; y < height; y++) {
      for (x = 0; x < width-width%16; x += 16) {
        __m256i x0, x1, x2, x3;
        __m256i cmp0, cmp1, r0, r1, r2, r3, r4;
        x0 = _mm256_loadu_si256((__m256i*)(src + x));
        cmp0 = _mm256_loadu_si256((__m256i*)(src + x + a_stride));
        cmp1 = _mm256_loadu_si256((__m256i*)(src + x + b_stride));
        r2 = _mm256_min_epu16(x0, cmp0);
        x1 = _mm256_cmpeq_epi16(cmp0, r2);
        x2 = _mm256_cmpeq_epi16(x0, r2);
        x1 = _mm256_sub_epi16(x2, x1);
        r2 = _mm256_min_epu16(x0, cmp1);
        x3 = _mm256_cmpeq_epi16(cmp1, r2);
        x2 = _mm256_cmpeq_epi16(x0, r2);
        x3 = _mm256_sub_epi16(x2, x3);
        x1 = _mm256_add_epi16(x1, x3);
        r0 = _mm256_cmpeq_epi16(x1, _mm256_set1_epi16(-2));
        r1 = _mm256_cmpeq_epi16(x1, _mm256_set1_epi16(-1));
        r3 = _mm256_cmpeq_epi16(x1, _mm256_set1_epi16(1));
        r4 = _mm256_cmpeq_epi16(x1, _mm256_set1_epi16(2));
        r0 = _mm256_and_si256(r0, offset0);
        r1 = _mm256_and_si256(r1, offset1);
        r3 = _mm256_and_si256(r3, offset3);
        r4 = _mm256_and_si256(r4, offset4);
        r0 = _mm256_add_epi16(r0, r1);
        r0 = _mm256_add_epi16(r0, r4);
        r0 = _mm256_add_epi16(r0, r3);
        r0 = _mm256_add_epi16(r0, x0);
        r1 = _mm256_set1_epi16(0x03FF);
        r0 = _mm256_max_epi16(r0, _mm256_setzero_si256());
        r0 = _mm256_min_epi16(r0, r1);
        _mm256_storeu_si256((__m256i*)(dst + x), r0);
      }
      for (; x < width - width%8; x += 8) {
        __m128i x0, x1, x2, x3;
        __m128i cmp0, cmp1, r0, r1, r2, r3, r4;
        x0 = _mm_loadu_si128((__m128i*)(src + x));
        cmp0 = _mm_loadu_si128((__m128i*)(src + x + a_stride));
        cmp1 = _mm_loadu_si128((__m128i*)(src + x + b_stride));
        r2 = _mm_min_epu16(x0, cmp0);
        x1 = _mm_cmpeq_epi16(cmp0, r2);
        x2 = _mm_cmpeq_epi16(x0, r2);
        x1 = _mm_sub_epi16(x2, x1);
        r2 = _mm_min_epu16(x0, cmp1);
        x3 = _mm_cmpeq_epi16(cmp1, r2);
        x2 = _mm_cmpeq_epi16(x0, r2);
        x3 = _mm_sub_epi16(x2, x3);
        x1 = _mm_add_epi16(x1, x3);
        r0 = _mm_cmpeq_epi16(x1, _mm_set1_epi16(-2));
        r1 = _mm_cmpeq_epi16(x1, _mm_set1_epi16(-1));
        r2 = _mm_cmpeq_epi16(x1, _mm_set1_epi16(0));
        r3 = _mm_cmpeq_epi16(x1, _mm_set1_epi16(1));
        r4 = _mm_cmpeq_epi16(x1, _mm_set1_epi16(2));
        r0 = _mm_and_si128(r0, _mm256_castsi256_si128(offset0));
        r1 = _mm_and_si128(r1, _mm256_castsi256_si128(offset1));
        r3 = _mm_and_si128(r3, _mm256_castsi256_si128(offset3));
        r4 = _mm_and_si128(r4, _mm256_castsi256_si128(offset4));
        r0 = _mm_add_epi16(r0, r1);
        r0 = _mm_add_epi16(r0, r4);
        r0 = _mm_add_epi16(r0, r3);
        r0 = _mm_add_epi16(r0, x0);
        r1 = _mm_set1_epi16(0x03FF);
        r0 = _mm_max_epi16(r0, _mm_setzero_si128());
        r0 = _mm_min_epi16(r0, r1);
        _mm_storeu_si128((__m128i*)(dst + x), r0);
      }
      {
        __m128i x0, x1, x2, x3;
        __m128i cmp0, cmp1, r0, r1, r2, r3, r4;
        x0 = _mm_loadu_si128((__m128i*)(src + x));
        cmp0 = _mm_loadu_si128((__m128i*)(src + x + a_stride));
        cmp1 = _mm_loadu_si128((__m128i*)(src + x + b_stride));
        r2 = _mm_min_epu16(x0, cmp0);
        x1 = _mm_cmpeq_epi16(cmp0, r2);
        x2 = _mm_cmpeq_epi16(x0, r2);
        x1 = _mm_sub_epi16(x2, x1);
        r2 = _mm_min_epu16(x0, cmp1);
        x3 = _mm_cmpeq_epi16(cmp1, r2);
        x2 = _mm_cmpeq_epi16(x0, r2);
        x3 = _mm_sub_epi16(x2, x3);
        x1 = _mm_add_epi16(x1, x3);
        r0 = _mm_cmpeq_epi16(x1, _mm_set1_epi16(-2));
        r1 = _mm_cmpeq_epi16(x1, _mm_set1_epi16(-1));
        r2 = _mm_cmpeq_epi16(x1, _mm_set1_epi16(0));
        r3 = _mm_cmpeq_epi16(x1, _mm_set1_epi16(1));
        r4 = _mm_cmpeq_epi16(x1, _mm_set1_epi16(2));
        r0 = _mm_and_si128(r0, _mm256_castsi256_si128(offset0));
        r1 = _mm_and_si128(r1, _mm256_castsi256_si128(offset1));
        r2 = _mm_and_si128(r2, _mm256_castsi256_si128(offset2));
        r3 = _mm_and_si128(r3, _mm256_castsi256_si128(offset3));
        r4 = _mm_and_si128(r4, _mm256_castsi256_si128(offset4));
        r0 = _mm_add_epi16(r0, r1);
        r2 = _mm_add_epi16(r2, r3);
        r0 = _mm_add_epi16(r0, r4);
        r0 = _mm_add_epi16(r0, r2);

        //mask to remove processing on last element
        r2 = _mm_set1_epi16((int16_t)0xFFFF);
        r2 = _mm_bsrli_si128(r2,2);
        r0 = _mm_and_si128(r0, r2);

        r0 = _mm_add_epi16(r0, x0);
        r1 = _mm_set1_epi16(0x03FF);
        r0 = _mm_max_epi16(r0, _mm_setzero_si128());
        r0 = _mm_min_epi16(r0, r1);
        _mm_storeu_si128((__m128i*)(dst + x), r0);
      }
      src += stride_src;
      dst += stride_dst;
    }
  }
}

void rcn_init_sao_functions_avx2(struct RCNFunctions *const rcn_funcs){
    rcn_funcs->sao.band= &sao_band_filter_0_10_avx2;
    rcn_funcs->sao.edge[0]= &sao_edge_filter_7_10_avx2;
    rcn_funcs->sao.edge[1]= &sao_edge_filter_10_avx2;
}
