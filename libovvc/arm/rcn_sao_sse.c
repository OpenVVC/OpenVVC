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


#include "dec_structures.h"
#include "rcn_structures.h"

static void
sao_band_filter_0_10_sse(OVSample* _dst,
                         OVSample* _src,
                         ptrdiff_t _stride_dst,
                         ptrdiff_t _stride_src,
                         int width,
                         int height,
                         int8_t offset_val[],
                         uint8_t band_pos)
{
  int y, x;
  int shift = 10 - 5;
  __m128i r0, r1, r2, r3, x0, x1, x2, x3, sao1, sao2, sao3, sao4;
  __m128i src0, src2;
  uint16_t* dst = (uint16_t*)_dst;
  uint16_t* src = (uint16_t*)_src;
  ptrdiff_t stride_dst = _stride_dst;
  ptrdiff_t stride_src = _stride_src;
  r0 = _mm_set1_epi16((band_pos    ) & 31);
  r1 = _mm_set1_epi16((band_pos + 1) & 31);
  r2 = _mm_set1_epi16((band_pos + 2) & 31);
  r3 = _mm_set1_epi16((band_pos + 3) & 31);
  sao1 = _mm_set1_epi16(offset_val[0]);
  sao2 = _mm_set1_epi16(offset_val[1]);
  sao3 = _mm_set1_epi16(offset_val[2]);
  sao4 = _mm_set1_epi16(offset_val[3]);
  for (y = 0; y < height; y++) {
    for (x = 0; x < width; x += 8) {
      src0 = _mm_loadu_si128((__m128i*)&src[x]);
      src2 = _mm_srai_epi16(src0, shift);
      x0 = _mm_cmpeq_epi16(src2, r0);
      x1 = _mm_cmpeq_epi16(src2, r1);
      x2 = _mm_cmpeq_epi16(src2, r2);
      x3 = _mm_cmpeq_epi16(src2, r3);
      x0 = _mm_and_si128(x0, sao1);
      x1 = _mm_and_si128(x1, sao2);
      x2 = _mm_and_si128(x2, sao3);
      x3 = _mm_and_si128(x3, sao4);
      x0 = _mm_or_si128(x0, x1);
      x2 = _mm_or_si128(x2, x3);
      x0 = _mm_or_si128(x0, x2);
      src0 = _mm_add_epi16(src0, x0);
      src0 = _mm_max_epi16(src0, _mm_setzero_si128());
      src0 = _mm_min_epi16(src0, _mm_set1_epi16(0x03FF));
      _mm_store_si128((__m128i*)&dst[x], src0);
    }
    dst += stride_dst;
    src += stride_src;
  }
}

static void
sao_edge_filter_10_sse(OVSample* _dst,
                       OVSample* _src,
                       ptrdiff_t _stride_dst,
                       ptrdiff_t _stride_src,
                       int width,
                       int height,
                       int8_t offset_val[],
                       uint8_t eo_dir)
{
  int x, y;
  static const int8_t pos[4][2][2] = {
    { { -1, 0 }, { 1, 0 } },
    { { 0, -1 }, { 0, 1 } },
    { { -1, -1 }, { 1, 1 } },
    { { 1, -1 }, { -1, 1 } },
  };
  __m128i x0, x1, x2, x3, offset0, offset1, offset2, offset3, offset4;
  __m128i cmp0, cmp1, r0, r1, r2, r3, r4;
  uint16_t* dst = (uint16_t*)_dst;
  uint16_t* src = (uint16_t*)_src;
  ptrdiff_t stride_dst = _stride_dst;
  ptrdiff_t stride_src = _stride_src;
  {
    int a_stride = pos[eo_dir][0][0] + pos[eo_dir][0][1] * stride_src;
    int b_stride = pos[eo_dir][1][0] + pos[eo_dir][1][1] * stride_src;
    offset0 = _mm_set1_epi16(offset_val[0]);
    offset1 = _mm_set1_epi16(offset_val[1]);
    offset3 = _mm_set1_epi16(offset_val[2]);
    offset4 = _mm_set1_epi16(offset_val[3]);
    for (y = 0; y < height; y++) {
      for (x = 0; x < width; x += 8) {
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
        r3 = _mm_cmpeq_epi16(x1, _mm_set1_epi16(1));
        r4 = _mm_cmpeq_epi16(x1, _mm_set1_epi16(2));
        r0 = _mm_and_si128(r0, offset0);
        r1 = _mm_and_si128(r1, offset1);
        r3 = _mm_and_si128(r3, offset3);
        r4 = _mm_and_si128(r4, offset4);
        r0 = _mm_add_epi16(r0, r1);
        r0 = _mm_add_epi16(r0, r4);
        r0 = _mm_add_epi16(r0, r3);
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
sao_edge_filter_7_10_sse(OVSample* _dst,
                         OVSample* _src,
                         ptrdiff_t _stride_dst,
                         ptrdiff_t _stride_src,
                         int width,
                         int height,
                         int8_t offset_val[],
                         uint8_t eo_dir)
{
  int x, y;
  static const int8_t pos[4][2][2] = {
    { { -1, 0 }, { 1, 0 } },
    { { 0, -1 }, { 0, 1 } },
    { { -1, -1 }, { 1, 1 } },
    { { 1, -1 }, { -1, 1 } },
  };
  __m128i x0, x1, x2, x3, offset0, offset1, offset2, offset3, offset4;
  __m128i cmp0, cmp1, r0, r1, r2, r3, r4;
  uint16_t* dst = (uint16_t*)_dst;
  uint16_t* src = (uint16_t*)_src;
  ptrdiff_t stride_dst = _stride_dst;
  ptrdiff_t stride_src = _stride_src;
  {
    int a_stride = pos[eo_dir][0][0] + pos[eo_dir][0][1] * stride_src;
    int b_stride = pos[eo_dir][1][0] + pos[eo_dir][1][1] * stride_src;
    offset0 = _mm_set1_epi16(offset_val[0]);
    offset1 = _mm_set1_epi16(offset_val[1]);
    offset3 = _mm_set1_epi16(offset_val[2]);
    offset4 = _mm_set1_epi16(offset_val[3]);
    for (y = 0; y < height; y++) {
      for (x = 0; x < width - width%8; x += 8) {
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
        r3 = _mm_cmpeq_epi16(x1, _mm_set1_epi16(1));
        r4 = _mm_cmpeq_epi16(x1, _mm_set1_epi16(2));
        r0 = _mm_and_si128(r0, offset0);
        r1 = _mm_and_si128(r1, offset1);
        r3 = _mm_and_si128(r3, offset3);
        r4 = _mm_and_si128(r4, offset4);
        r0 = _mm_add_epi16(r0, r1);
        r0 = _mm_add_epi16(r0, r4);
        r0 = _mm_add_epi16(r0, r3);
        r0 = _mm_add_epi16(r0, x0);
        r1 = _mm_set1_epi16(0x03FF);
        r0 = _mm_max_epi16(r0, _mm_setzero_si128());
        r0 = _mm_min_epi16(r0, r1);
        _mm_storeu_si128((__m128i*)(dst + x), r0);
      }
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
      r3 = _mm_cmpeq_epi16(x1, _mm_set1_epi16(1));
      r4 = _mm_cmpeq_epi16(x1, _mm_set1_epi16(2));
      r0 = _mm_and_si128(r0, offset0);
      r1 = _mm_and_si128(r1, offset1);
      r3 = _mm_and_si128(r3, offset3);
      r4 = _mm_and_si128(r4, offset4);
      r0 = _mm_add_epi16(r0, r1);
      r0 = _mm_add_epi16(r0, r4);
      r0 = _mm_add_epi16(r0, r3);

      //mask to remove processing on last element
      r2 = _mm_set1_epi16((int16_t)0xFFFF);
      r2 = _mm_bsrli_si128(r2,2);
      r0 = _mm_and_si128(r0, r2);

      r0 = _mm_add_epi16(r0, x0);
      r1 = _mm_set1_epi16(0x03FF);
      r0 = _mm_max_epi16(r0, _mm_setzero_si128());
      r0 = _mm_min_epi16(r0, r1);
      _mm_storeu_si128((__m128i*)(dst + x), r0);
      src += stride_src;
      dst += stride_dst;
    }
  }
}

static void
sao_edge_filter_v_sse(OVSample *dst, OVSample *src_row, OVSample *src_col,
                      ptrdiff_t stride_dst, ptrdiff_t stride_src,
                      int width, int height,
                      int8_t offset_val[],
                      uint8_t eo_dir)
{
    int x, y;
    __m128i x0, x1, x2, x3;
    __m128i r0, r1, r2, r3, r4;

    const __m128i offset0 = _mm_set1_epi16(offset_val[0]);
    const __m128i offset1 = _mm_set1_epi16(offset_val[1]);
    const __m128i offset3 = _mm_set1_epi16(offset_val[2]);
    const __m128i offset4 = _mm_set1_epi16(offset_val[3]);

    for (x = 0; x < width; x += 8) {
        OVSample *src = dst;
        __m128i a = _mm_loadu_si128((__m128i*)(src_row + x));
        __m128i c = _mm_load_si128((__m128i*)(src));

        for (y = 0; y < height; y++) {

            __m128i b = _mm_load_si128((__m128i*)(src + stride_dst));

            r2 = _mm_min_epu16(c, a);
            x1 = _mm_cmpeq_epi16(a, r2);
            x2 = _mm_cmpeq_epi16(c, r2);
            x1 = _mm_sub_epi16(x2, x1);

            r2 = _mm_min_epu16(c, b);
            x3 = _mm_cmpeq_epi16(b, r2);
            x2 = _mm_cmpeq_epi16(c, r2);
            x3 = _mm_sub_epi16(x2, x3);

            x1 = _mm_add_epi16(x1, x3);

            r0 = _mm_cmpeq_epi16(x1, _mm_set1_epi16(-2));
            r1 = _mm_cmpeq_epi16(x1, _mm_set1_epi16(-1));
            r3 = _mm_cmpeq_epi16(x1, _mm_set1_epi16(1));
            r4 = _mm_cmpeq_epi16(x1, _mm_set1_epi16(2));

            r0 = _mm_and_si128(r0, offset0);
            r1 = _mm_and_si128(r1, offset1);
            r3 = _mm_and_si128(r3, offset3);
            r4 = _mm_and_si128(r4, offset4);

            r0 = _mm_add_epi16(r0, r1);
            r0 = _mm_add_epi16(r0, r4);
            r0 = _mm_add_epi16(r0, r3);
            r0 = _mm_add_epi16(r0, c);

            r1 = _mm_set1_epi16(0x03FF);

            r0 = _mm_max_epi16(r0, _mm_setzero_si128());
            r0 = _mm_min_epi16(r0, r1);

            _mm_store_si128((__m128i*)src, r0);

            a = c;
            c = b;

            src += stride_dst;
        }
        dst += 8;
    }
}

static void
sao_edge_filter_d_sse(OVSample *dst, OVSample *src_row, OVSample *src_col,
                      ptrdiff_t stride_dst, ptrdiff_t stride_src,
                      int width, int height,
                      int8_t offset_val[],
                      uint8_t eo_dir)
{
    int x, y;
    __m128i x0, x1, x2, x3;
    __m128i cmp0, cmp1, r0, r1, r2, r3, r4;

    const __m128i offset0 = _mm_set1_epi16(offset_val[0]);
    const __m128i offset1 = _mm_set1_epi16(offset_val[1]);
    const __m128i offset3 = _mm_set1_epi16(offset_val[2]);
    const __m128i offset4 = _mm_set1_epi16(offset_val[3]);

    for (x = 0; x < width; x += 8) {
        OVSample *src = dst;
        __m128i _a = _mm_loadu_si128((__m128i*)(src_row + x));
        __m128i c = _mm_load_si128((__m128i*)(src + x));
        for (y = 0; y < height; y++) {
            __m128i _b = _mm_load_si128((__m128i*)(src + stride_dst + x));
            __m128i b_ = _mm_load_si128((__m128i*)(src + stride_dst + x + 8));
            __m128i a0 = _mm_set1_epi16(src_col[y - 1]);
            src_col[y - 1] = _mm_extract_epi16(_a, 7);
            //__m128i c0 = _mm_set1_epi16(src_col[y]);

            __m128i a = _mm_alignr_epi8(_a, a0, 14);
            __m128i b = _mm_alignr_epi8(b_, _b, 2);

            r2 = _mm_min_epu16(c, a);
            x1 = _mm_cmpeq_epi16(a, r2);
            x2 = _mm_cmpeq_epi16(c, r2);
            x1 = _mm_sub_epi16(x2, x1);

            r2 = _mm_min_epu16(c, b);
            x3 = _mm_cmpeq_epi16(b, r2);
            x2 = _mm_cmpeq_epi16(c, r2);
            x3 = _mm_sub_epi16(x2, x3);

            x1 = _mm_add_epi16(x1, x3);

            r0 = _mm_cmpeq_epi16(x1, _mm_set1_epi16(-2));
            r1 = _mm_cmpeq_epi16(x1, _mm_set1_epi16(-1));
            r3 = _mm_cmpeq_epi16(x1, _mm_set1_epi16(1));
            r4 = _mm_cmpeq_epi16(x1, _mm_set1_epi16(2));

            r0 = _mm_and_si128(r0, offset0);
            r1 = _mm_and_si128(r1, offset1);
            r3 = _mm_and_si128(r3, offset3);
            r4 = _mm_and_si128(r4, offset4);

            r0 = _mm_add_epi16(r0, r1);
            r0 = _mm_add_epi16(r0, r4);
            r0 = _mm_add_epi16(r0, r3);
            r0 = _mm_add_epi16(r0, c);

            r1 = _mm_set1_epi16(0x03FF);

            r0 = _mm_max_epi16(r0, _mm_setzero_si128());
            r0 = _mm_min_epi16(r0, r1);

            _mm_store_si128((__m128i*)(src + x), r0);

            _a = c;
            c = _b;
            src += stride_dst;
        }
    }
}

static void
sao_edge_filter_b_sse(OVSample *dst, OVSample *src_row, OVSample *src_col,
                      ptrdiff_t stride_dst, ptrdiff_t stride_src,
                      int width, int height,
                      int8_t offset_val[],
                      uint8_t eo_dir)
{
    int x, y;
    __m128i x0, x1, x2, x3;
    __m128i cmp0, cmp1, r0, r1, r2, r3, r4;

    const __m128i offset0 = _mm_set1_epi16(offset_val[0]);
    const __m128i offset1 = _mm_set1_epi16(offset_val[1]);
    const __m128i offset3 = _mm_set1_epi16(offset_val[2]);
    const __m128i offset4 = _mm_set1_epi16(offset_val[3]);

    for (x = 0; x < width; x += 8) {
        OVSample *src = dst;
        __m128i _a = _mm_loadu_si128((__m128i*)(src_row + x));
        __m128i a_ = _mm_loadu_si128((__m128i*)(src_row + x + 8));
        __m128i c = _mm_load_si128((__m128i*)(src + x));
        for (y = 0; y < height; y++) {
            __m128i c_ = _mm_load_si128((__m128i*)(src + x + 8));
            __m128i _b = _mm_load_si128((__m128i*)(src + stride_dst + x));
            __m128i b0 = _mm_set1_epi16(src_col[y + 1]);

            __m128i b = _mm_alignr_epi8(_b, b0, 14);
            __m128i a = _mm_alignr_epi8(a_, _a, 2);

            src_col[y + 1] = _mm_extract_epi16(_b, 7);

            r2 = _mm_min_epu16(c, a);
            x1 = _mm_cmpeq_epi16(a, r2);
            x2 = _mm_cmpeq_epi16(c, r2);
            x1 = _mm_sub_epi16(x2, x1);

            r2 = _mm_min_epu16(c, b);
            x3 = _mm_cmpeq_epi16(b, r2);
            x2 = _mm_cmpeq_epi16(c, r2);
            x3 = _mm_sub_epi16(x2, x3);

            x1 = _mm_add_epi16(x1, x3);

            r0 = _mm_cmpeq_epi16(x1, _mm_set1_epi16(-2));
            r1 = _mm_cmpeq_epi16(x1, _mm_set1_epi16(-1));
            r3 = _mm_cmpeq_epi16(x1, _mm_set1_epi16(1));
            r4 = _mm_cmpeq_epi16(x1, _mm_set1_epi16(2));

            r0 = _mm_and_si128(r0, offset0);
            r1 = _mm_and_si128(r1, offset1);
            r3 = _mm_and_si128(r3, offset3);
            r4 = _mm_and_si128(r4, offset4);

            r0 = _mm_add_epi16(r0, r1);
            r0 = _mm_add_epi16(r0, r4);
            r0 = _mm_add_epi16(r0, r3);
            r0 = _mm_add_epi16(r0, c);

            r1 = _mm_set1_epi16(0x03FF);

            r0 = _mm_max_epi16(r0, _mm_setzero_si128());
            r0 = _mm_min_epi16(r0, r1);

            _mm_store_si128((__m128i*)(src + x), r0);

            _a = c;
            a_ = c_;
            c = _b;
            src += stride_dst;
        }
    }
}


static void
sao_edge_filter_h_sse(OVSample *dst, OVSample *src_row, OVSample *src_col,
                      ptrdiff_t stride_dst, ptrdiff_t stride_src,
                      int width, int height,
                      int8_t offset_val[],
                      uint8_t eo_dir)
{
    int x, y;
    __m128i x0, x1, x2, x3;
    __m128i cmp0, cmp1, r0, r1, r2, r3, r4;
    OVSample *src = dst;

    const __m128i offset0 = _mm_set1_epi16(offset_val[0]);
    const __m128i offset1 = _mm_set1_epi16(offset_val[1]);
    const __m128i offset3 = _mm_set1_epi16(offset_val[2]);
    const __m128i offset4 = _mm_set1_epi16(offset_val[3]);

    for (y = 0; y < height; y++) {
        __m128i prev = _mm_set1_epi16(src_col[y]);
        __m128i c = _mm_load_si128((__m128i*)(src));
        for (x = 0; x < width; x += 8) {

            __m128i next = _mm_load_si128((__m128i*)(src + x + 8));

            __m128i a = _mm_alignr_epi8(c, prev, 14);
            __m128i b = _mm_alignr_epi8(next, c, 2);

            r2 = _mm_min_epu16(c, a);
            x1 = _mm_cmpeq_epi16(a, r2);
            x2 = _mm_cmpeq_epi16(c, r2);
            x1 = _mm_sub_epi16(x2, x1);

            r2 = _mm_min_epu16(c, b);
            x3 = _mm_cmpeq_epi16(b, r2);
            x2 = _mm_cmpeq_epi16(c, r2);
            x3 = _mm_sub_epi16(x2, x3);

            x1 = _mm_add_epi16(x1, x3);

            r0 = _mm_cmpeq_epi16(x1, _mm_set1_epi16(-2));
            r1 = _mm_cmpeq_epi16(x1, _mm_set1_epi16(-1));
            r3 = _mm_cmpeq_epi16(x1, _mm_set1_epi16(1));
            r4 = _mm_cmpeq_epi16(x1, _mm_set1_epi16(2));

            r0 = _mm_and_si128(r0, offset0);
            r1 = _mm_and_si128(r1, offset1);
            r3 = _mm_and_si128(r3, offset3);
            r4 = _mm_and_si128(r4, offset4);

            r0 = _mm_add_epi16(r0, r1);
            r0 = _mm_add_epi16(r0, r4);
            r0 = _mm_add_epi16(r0, r3);
            r0 = _mm_add_epi16(r0, c);

            r1 = _mm_set1_epi16(0x03FF);

            r0 = _mm_max_epi16(r0, _mm_setzero_si128());
            r0 = _mm_min_epi16(r0, r1);

            _mm_store_si128((__m128i*)(dst + x), r0);

            prev = c;
            c = next;
        }
        src += stride_dst;
        dst += stride_dst;
    }
}

void rcn_init_sao_functions_sse(struct RCNFunctions *const rcn_funcs){
    rcn_funcs->sao.band= &sao_band_filter_0_10_sse;
    rcn_funcs->sao.edge[0]= &sao_edge_filter_7_10_sse;
    rcn_funcs->sao.edge[1]= &sao_edge_filter_10_sse;
    rcn_funcs->sao.edge2[0]= &sao_edge_filter_h_sse;
    rcn_funcs->sao.edge2[1]= &sao_edge_filter_v_sse;
    rcn_funcs->sao.edge2[2]= &sao_edge_filter_d_sse;
    rcn_funcs->sao.edge2[3]= &sao_edge_filter_b_sse;
}
