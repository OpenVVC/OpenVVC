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

#include "ovmem.h"
#include "rcn_structures.h"
#include <immintrin.h>
#include <stdint.h>
#include <string.h>

#define SIZE_BLOCK_2 0
#define SIZE_BLOCK_4 1
#define SIZE_BLOCK_8 2
#define SIZE_BLOCK_16 3
#define SIZE_BLOCK_32 4
#define SIZE_BLOCK_64 5
#define SIZE_BLOCK_128 6

#define MAX_PB_SIZE 128

#define EPEL_EXTRA_BEFORE 1
#define EPEL_EXTRA_AFTER 2
#define EPEL_EXTRA EPEL_EXTRA_BEFORE + EPEL_EXTRA_AFTER

#define QPEL_EXTRA_BEFORE 3
#define QPEL_EXTRA_AFTER 4
#define QPEL_EXTRA QPEL_EXTRA_BEFORE + QPEL_EXTRA_AFTER

#define BITDEPTH 10

static const uint32_t ov_bilinear_filters_4[15] =
{
  /*0x00000010,*/
    0x0001000f,
    0x0002000e,
    0x0003000d,
    0x0004000c,
    0x0005000b,
    0x0006000a,
    0x00070009,
    0x00080008,
    0x00090007,
    0x000a0006,
    0x000b0005,
    0x000c0004,
    0x000d0003,
    0x000e0002,
    0x000f0001,
};

DECLARE_ALIGNED(16, static const uint32_t, ov_mcp_filters_l[16][4]) = {
    {0x00010000, 0x003ffffd, 0xfffe0004, 0x00000001},
    {0x0002ffff, 0x003efffb, 0xfffd0008, 0x00000001},
    {0x0003ffff, 0x003cfff8, 0xfffc000d, 0x00000001},
    {0x0004ffff, 0x003afff6, 0xfffb0011, 0x00000001},
    {0x0004ffff, 0x0034fff5, 0xfff8001a, 0xffff0003},
    {0x0003ffff, 0x002ffff7, 0xfff6001f, 0xffff0004},
    {0x0004ffff, 0x002dfff5, 0xfff60022, 0xffff0004},
    {0x0004ffff, 0x0028fff5, 0xfff50028, 0xffff0004},
    {0x0004ffff, 0x0022fff6, 0xfff5002d, 0xffff0004},
    {0x0004ffff, 0x001ffff6, 0xfff7002f, 0xffff0003},
    {0x0003ffff, 0x001afff8, 0xfff50034, 0xffff0004},
    {0x00010000, 0x0011fffb, 0xfff6003a, 0xffff0004},
    {0x00010000, 0x000dfffc, 0xfff8003c, 0xffff0003},
    {0x00010000, 0x0008fffd, 0xfffb003e, 0xffff0002},
    {0x00010000, 0x0004fffe, 0xfffd003f, 0x00000001},
    {0x00030000, 0x00140009, 0x00090014, 0x00000003}
};

DECLARE_ALIGNED(16, static const uint32_t, ov_mcp_filters_c[31][2]) = {
  /*{0x00400000, 0x00000000},*/
    {0x003fffff, 0x00000002},
    {0x003efffe, 0x00000004},
    {0x003cfffe, 0xffff0007},
    {0x003afffe, 0xfffe000a},
    {0x0039fffd, 0xfffe000c},
    {0x0038fffc, 0xfffe000e},
    {0x0037fffc, 0xfffe000f},
    {0x0036fffc, 0xfffe0010},
    {0x0035fffb, 0xfffe0012},
    {0x0034fffa, 0xfffe0014},
    {0x0031fffa, 0xfffd0018},
    {0x002efffa, 0xfffc001c},
    {0x002cfffb, 0xfffc001d},
    {0x002afffc, 0xfffc001e},
    {0x0027fffc, 0xfffc0021},
    {0x0024fffc, 0xfffc0024},
    {0x0021fffc, 0xfffc0027},
    {0x001efffc, 0xfffc002a},
    {0x001dfffc, 0xfffb002c},
    {0x001cfffc, 0xfffa002e},
    {0x0018fffd, 0xfffa0031},
    {0x0014fffe, 0xfffa0034},
    {0x0012fffe, 0xfffb0035},
    {0x0010fffe, 0xfffc0036},
    {0x000ffffe, 0xfffc0037},
    {0x000efffe, 0xfffc0038},
    {0x000cfffe, 0xfffd0039},
    {0x000afffe, 0xfffe003a},
    {0x0007ffff, 0xfffe003c},
    {0x00040000, 0xfffe003e},
    {0x00020000, 0xffff003f}
};

static void
oh_hevc_put_hevc_bi0_pel_pixels16_10_avx2(int16_t* dst,
                                          const uint16_t* _src,
                                          ptrdiff_t _srcstride,
                                          int height,
                                          intptr_t mx,
                                          intptr_t my,
                                          int width)
{
    int x, y;
    __m256i x1;
    const uint16_t* src = _src;
    const int srcstride = _srcstride;
    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x += 16) {
            x1 = _mm256_loadu_si256((__m256i*)&src[x]);
            x1 = _mm256_slli_epi16(x1, 14 - BITDEPTH);
            _mm256_storeu_si256((__m256i*)&dst[x], x1);
        }
        src += srcstride;
        dst += MAX_PB_SIZE;
    }
}

static void
oh_hevc_put_hevc_bi1_pel_pixels16_10_avx2(uint16_t* _dst,
                                          ptrdiff_t _dststride,
                                          const uint16_t* _src,
                                          ptrdiff_t _srcstride,
                                          const int16_t* src2,
                                          int height,
                                          intptr_t mx,
                                          intptr_t my,
                                          int width)
{
    int x, y;
    __m256i x1;
    const uint16_t* src = _src;
    const int srcstride = _srcstride;
    const __m256i offset = _mm256_set1_epi16(1 << BITDEPTH);
    uint16_t* dst = (uint16_t*)_dst;
    const int dststride = _dststride;
    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x += 16) {
            __m256i r5;
            x1 = _mm256_loadu_si256((__m256i*)&src[x]);
            x1 = _mm256_slli_epi16(x1, 14 - BITDEPTH);
            r5 = _mm256_load_si256((__m256i*)&src2[x]);
            x1 = _mm256_adds_epi16(x1, r5);
            x1 = _mm256_mulhrs_epi16(x1, offset);
            x1 = _mm256_max_epi16(x1, _mm256_setzero_si256());
            x1 = _mm256_min_epi16(x1, _mm256_set1_epi16(0x03FF));
            _mm256_storeu_si256((__m256i*)&dst[x], x1);
        }
        src += srcstride;
        src2 += MAX_PB_SIZE;
        dst += dststride;
    }
}

static void
oh_hevc_put_hevc_bi0_epel_h16_10_avx2(int16_t* dst,
                                      const uint16_t* _src,
                                      ptrdiff_t _srcstride,
                                      int height,
                                      intptr_t mx,
                                      intptr_t my,
                                      int width)
{
    int x, y;
    uint16_t* src = ((uint16_t*)_src) - 1;
    ptrdiff_t srcstride = _srcstride;
    __m256i f1 = _mm256_set1_epi32(ov_mcp_filters_c[mx - 1][0]);
    __m256i f2 = _mm256_set1_epi32(ov_mcp_filters_c[mx - 1][1]);
    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x += 16) {
            __m256i x1, x2, x3, x4, t1, t2;
            x1 = _mm256_loadu_si256((__m256i*)&src[x]);
            x2 = _mm256_loadu_si256((__m256i*)&src[x + 1]);
            x3 = _mm256_loadu_si256((__m256i*)&src[x + 2 * 1]);
            x4 = _mm256_loadu_si256((__m256i*)&src[x + 3 * 1]);
            t1 = _mm256_unpackhi_epi16(x1, x2);
            x1 = _mm256_unpacklo_epi16(x1, x2);
            t2 = _mm256_unpackhi_epi16(x3, x4);
            x2 = _mm256_unpacklo_epi16(x3, x4);
            x1 = _mm256_madd_epi16(x1, f1);
            t1 = _mm256_madd_epi16(t1, f1);
            x2 = _mm256_madd_epi16(x2, f2);
            t2 = _mm256_madd_epi16(t2, f2);
            x1 = _mm256_add_epi32(x1, x2);
            t1 = _mm256_add_epi32(t1, t2);
            t1 = _mm256_srai_epi32(t1, 2);
            x1 = _mm256_srai_epi32(x1, 2);
            x1 = _mm256_packs_epi32(x1, t1);
            _mm256_storeu_si256((__m256i*)&dst[x], x1);
        }
        src += srcstride;
        dst += MAX_PB_SIZE;
    }
}

static void
oh_hevc_put_hevc_bi1_epel_h16_10_avx2(uint16_t* _dst,
                                      ptrdiff_t _dststride,
                                      const uint16_t* _src,
                                      ptrdiff_t _srcstride,
                                      const int16_t* src2,
                                      int height,
                                      intptr_t mx,
                                      intptr_t my,
                                      int width)
{
    int x, y;
    uint16_t* src = ((uint16_t*)_src) - 1;
    ptrdiff_t srcstride = _srcstride;
    const __m256i offset = _mm256_set1_epi16(1 << BITDEPTH);
    uint16_t* dst = (uint16_t*)_dst;
    const int dststride = _dststride;
    __m256i f1 = _mm256_set1_epi32(ov_mcp_filters_c[mx - 1][0]);
    __m256i f2 = _mm256_set1_epi32(ov_mcp_filters_c[mx - 1][1]);
    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x += 16) {
            __m256i x1, x2, x3, x4, t1, t2;
            __m256i r5;
            x1 = _mm256_loadu_si256((__m256i*)&src[x]);
            x2 = _mm256_loadu_si256((__m256i*)&src[x + 1]);
            x3 = _mm256_loadu_si256((__m256i*)&src[x + 2 * 1]);
            x4 = _mm256_loadu_si256((__m256i*)&src[x + 3 * 1]);
            t1 = _mm256_unpackhi_epi16(x1, x2);
            x1 = _mm256_unpacklo_epi16(x1, x2);
            t2 = _mm256_unpackhi_epi16(x3, x4);
            x2 = _mm256_unpacklo_epi16(x3, x4);
            x1 = _mm256_madd_epi16(x1, f1);
            t1 = _mm256_madd_epi16(t1, f1);
            x2 = _mm256_madd_epi16(x2, f2);
            t2 = _mm256_madd_epi16(t2, f2);
            x1 = _mm256_add_epi32(x1, x2);
            t1 = _mm256_add_epi32(t1, t2);
            t1 = _mm256_srai_epi32(t1, 2);
            x1 = _mm256_srai_epi32(x1, 2);
            x1 = _mm256_packs_epi32(x1, t1);
            r5 = _mm256_load_si256((__m256i*)&src2[x]);
            x1 = _mm256_adds_epi16(x1, r5);
            x1 = _mm256_mulhrs_epi16(x1, offset);
            x1 = _mm256_max_epi16(x1, _mm256_setzero_si256());
            x1 = _mm256_min_epi16(x1, _mm256_set1_epi16(0x03FF));
            _mm256_storeu_si256((__m256i*)&dst[x], x1);
        }
        src += srcstride;
        src2 += MAX_PB_SIZE;
        dst += dststride;
    }
}

static void
oh_hevc_put_hevc_uni_epel_h16_10_avx2(uint16_t* _dst,
                                      ptrdiff_t _dststride,
                                      const uint16_t* _src,
                                      ptrdiff_t _srcstride,
                                      int height,
                                      intptr_t mx,
                                      intptr_t my,
                                      int width)
{
    int x, y;
    uint16_t* src = ((uint16_t*)_src) - 1;
    ptrdiff_t srcstride = _srcstride;
    const __m256i offset = _mm256_set1_epi16(1 << (BITDEPTH + 1));
    uint16_t* dst = (uint16_t*)_dst;
    const int dststride = _dststride;
    __m256i f1 = _mm256_set1_epi32(ov_mcp_filters_c[mx - 1][0]);
    __m256i f2 = _mm256_set1_epi32(ov_mcp_filters_c[mx - 1][1]);
    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x += 16) {
            __m256i x1, x2, x3, x4, t1, t2;
            x1 = _mm256_loadu_si256((__m256i*)&src[x]);
            x2 = _mm256_loadu_si256((__m256i*)&src[x + 1]);
            x3 = _mm256_loadu_si256((__m256i*)&src[x + 2 * 1]);
            x4 = _mm256_loadu_si256((__m256i*)&src[x + 3 * 1]);
            t1 = _mm256_unpackhi_epi16(x1, x2);
            x1 = _mm256_unpacklo_epi16(x1, x2);
            t2 = _mm256_unpackhi_epi16(x3, x4);
            x2 = _mm256_unpacklo_epi16(x3, x4);
            x1 = _mm256_madd_epi16(x1, f1);
            t1 = _mm256_madd_epi16(t1, f1);
            x2 = _mm256_madd_epi16(x2, f2);
            t2 = _mm256_madd_epi16(t2, f2);
            x1 = _mm256_add_epi32(x1, x2);
            t1 = _mm256_add_epi32(t1, t2);
            t1 = _mm256_srai_epi32(t1, 2);
            x1 = _mm256_srai_epi32(x1, 2);
            x1 = _mm256_packs_epi32(x1, t1);
            x1 = _mm256_mulhrs_epi16(x1, offset);
            x1 = _mm256_max_epi16(x1, _mm256_setzero_si256());
            x1 = _mm256_min_epi16(x1, _mm256_set1_epi16(0x03FF));
            _mm256_storeu_si256((__m256i*)&dst[x], x1);
        }
        src += srcstride;
        dst += dststride;
    }
}

static void
oh_hevc_put_hevc_bi0_epel_v16_10_avx2(int16_t* dst,
                                      const uint16_t* _src,
                                      ptrdiff_t _srcstride,
                                      int height,
                                      intptr_t mx,
                                      intptr_t my,
                                      int width)
{
    int x, y;
    ptrdiff_t srcstride = _srcstride;
    uint16_t* src = ((uint16_t*)_src) - srcstride;
    __m256i f1 = _mm256_set1_epi32(ov_mcp_filters_c[my - 1][0]);
    __m256i f2 = _mm256_set1_epi32(ov_mcp_filters_c[my - 1][1]);
    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x += 16) {
            __m256i x1, x2, x3, x4, t1, t2;
            x1 = _mm256_loadu_si256((__m256i*)&src[x]);
            x2 = _mm256_loadu_si256((__m256i*)&src[x + srcstride]);
            x3 = _mm256_loadu_si256((__m256i*)&src[x + 2 * srcstride]);
            x4 = _mm256_loadu_si256((__m256i*)&src[x + 3 * srcstride]);
            t1 = _mm256_unpackhi_epi16(x1, x2);
            x1 = _mm256_unpacklo_epi16(x1, x2);
            t2 = _mm256_unpackhi_epi16(x3, x4);
            x2 = _mm256_unpacklo_epi16(x3, x4);
            x1 = _mm256_madd_epi16(x1, f1);
            t1 = _mm256_madd_epi16(t1, f1);
            x2 = _mm256_madd_epi16(x2, f2);
            t2 = _mm256_madd_epi16(t2, f2);
            x1 = _mm256_add_epi32(x1, x2);
            t1 = _mm256_add_epi32(t1, t2);
            t1 = _mm256_srai_epi32(t1, 2);
            x1 = _mm256_srai_epi32(x1, 2);
            x1 = _mm256_packs_epi32(x1, t1);
            _mm256_storeu_si256((__m256i*)&dst[x], x1);
        }
        src += srcstride;
        dst += MAX_PB_SIZE;
    }
}

static void
oh_hevc_put_hevc_bi1_epel_v16_10_avx2(uint16_t* _dst,
                                      ptrdiff_t _dststride,
                                      const uint16_t* _src,
                                      ptrdiff_t _srcstride,
                                      const int16_t* src2,
                                      int height,
                                      intptr_t mx,
                                      intptr_t my,
                                      int width)
{
    int x, y;
    ptrdiff_t srcstride = _srcstride;
    uint16_t* src = ((uint16_t*)_src) - srcstride;
    const __m256i offset = _mm256_set1_epi16(1 << BITDEPTH);
    uint16_t* dst = (uint16_t*)_dst;
    const int dststride = _dststride;
    __m256i f1 = _mm256_set1_epi32(ov_mcp_filters_c[my - 1][0]);
    __m256i f2 = _mm256_set1_epi32(ov_mcp_filters_c[my - 1][1]);
    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x += 16) {
            __m256i x1, x2, x3, x4, t1, t2;
            __m256i r5;
            x1 = _mm256_loadu_si256((__m256i*)&src[x]);
            x2 = _mm256_loadu_si256((__m256i*)&src[x + srcstride]);
            x3 = _mm256_loadu_si256((__m256i*)&src[x + 2 * srcstride]);
            x4 = _mm256_loadu_si256((__m256i*)&src[x + 3 * srcstride]);
            t1 = _mm256_unpackhi_epi16(x1, x2);
            x1 = _mm256_unpacklo_epi16(x1, x2);
            t2 = _mm256_unpackhi_epi16(x3, x4);
            x2 = _mm256_unpacklo_epi16(x3, x4);
            x1 = _mm256_madd_epi16(x1, f1);
            t1 = _mm256_madd_epi16(t1, f1);
            x2 = _mm256_madd_epi16(x2, f2);
            t2 = _mm256_madd_epi16(t2, f2);
            x1 = _mm256_add_epi32(x1, x2);
            t1 = _mm256_add_epi32(t1, t2);
            t1 = _mm256_srai_epi32(t1, 2);
            x1 = _mm256_srai_epi32(x1, 2);
            x1 = _mm256_packs_epi32(x1, t1);
            r5 = _mm256_load_si256((__m256i*)&src2[x]);
            x1 = _mm256_adds_epi16(x1, r5);
            x1 = _mm256_mulhrs_epi16(x1, offset);
            x1 = _mm256_max_epi16(x1, _mm256_setzero_si256());
            x1 = _mm256_min_epi16(x1, _mm256_set1_epi16(0x03FF));
            _mm256_storeu_si256((__m256i*)&dst[x], x1);
        }
        src += srcstride;
        src2 += MAX_PB_SIZE;
        dst += dststride;
    }
}

static void
oh_hevc_put_hevc_uni_epel_v16_10_avx2(uint16_t* _dst,
                                      ptrdiff_t _dststride,
                                      const uint16_t* _src,
                                      ptrdiff_t _srcstride,
                                      int height,
                                      intptr_t mx,
                                      intptr_t my,
                                      int width)
{
    int x, y;
    ptrdiff_t srcstride = _srcstride;
    uint16_t* src = ((uint16_t*)_src) - srcstride;
    const __m256i offset = _mm256_set1_epi16(1 << (BITDEPTH + 1));
    uint16_t* dst = (uint16_t*)_dst;
    const int dststride = _dststride;
    __m256i f1 = _mm256_set1_epi32(ov_mcp_filters_c[my - 1][0]);
    __m256i f2 = _mm256_set1_epi32(ov_mcp_filters_c[my - 1][1]);
    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x += 16) {
            __m256i x1, x2, x3, x4, t1, t2;
            x1 = _mm256_loadu_si256((__m256i*)&src[x]);
            x2 = _mm256_loadu_si256((__m256i*)&src[x + srcstride]);
            x3 = _mm256_loadu_si256((__m256i*)&src[x + 2 * srcstride]);
            x4 = _mm256_loadu_si256((__m256i*)&src[x + 3 * srcstride]);
            t1 = _mm256_unpackhi_epi16(x1, x2);
            x1 = _mm256_unpacklo_epi16(x1, x2);
            t2 = _mm256_unpackhi_epi16(x3, x4);
            x2 = _mm256_unpacklo_epi16(x3, x4);
            x1 = _mm256_madd_epi16(x1, f1);
            t1 = _mm256_madd_epi16(t1, f1);
            x2 = _mm256_madd_epi16(x2, f2);
            t2 = _mm256_madd_epi16(t2, f2);
            x1 = _mm256_add_epi32(x1, x2);
            t1 = _mm256_add_epi32(t1, t2);
            t1 = _mm256_srai_epi32(t1, 2);
            x1 = _mm256_srai_epi32(x1, 2);
            x1 = _mm256_packs_epi32(x1, t1);
            x1 = _mm256_mulhrs_epi16(x1, offset);
            x1 = _mm256_max_epi16(x1, _mm256_setzero_si256());
            x1 = _mm256_min_epi16(x1, _mm256_set1_epi16(0x03FF));
            _mm256_storeu_si256((__m256i*)&dst[x], x1);
        }
        src += srcstride;
        dst += dststride;
    }
}

static void
oh_hevc_put_hevc_bi0_epel_hv16_10_avx2(int16_t* dst,
                                       const uint16_t* _src,
                                       ptrdiff_t _srcstride,
                                       int height,
                                       intptr_t mx,
                                       intptr_t my,
                                       int width)
{
    int x, y;
    int16_t* dst_bis = dst;
    uint16_t *src_bis, *src = ((uint16_t*)_src) - 1;
    ptrdiff_t srcstride = _srcstride;
    src -= 1 * srcstride;
    src_bis = src;
    __m256i f1 = _mm256_set1_epi32(ov_mcp_filters_c[mx - 1][0]);
    __m256i f2 = _mm256_set1_epi32(ov_mcp_filters_c[mx - 1][1]);
    __m256i f3 = _mm256_set1_epi32(ov_mcp_filters_c[my - 1][0]);
    __m256i f4 = _mm256_set1_epi32(ov_mcp_filters_c[my - 1][1]);
    for (x = 0; x < width; x += 16) {
        __m256i x1, x2, x3, x4, t1, t2, r1, r2, r3, r4;
        x1 = _mm256_loadu_si256((__m256i*)&src[x]);
        x2 = _mm256_loadu_si256((__m256i*)&src[x + 1]);
        x3 = _mm256_loadu_si256((__m256i*)&src[x + 2 * 1]);
        x4 = _mm256_loadu_si256((__m256i*)&src[x + 3 * 1]);
        t1 = _mm256_unpackhi_epi16(x1, x2);
        x1 = _mm256_unpacklo_epi16(x1, x2);
        t2 = _mm256_unpackhi_epi16(x3, x4);
        x2 = _mm256_unpacklo_epi16(x3, x4);
        x1 = _mm256_madd_epi16(x1, f1);
        t1 = _mm256_madd_epi16(t1, f1);
        x2 = _mm256_madd_epi16(x2, f2);
        t2 = _mm256_madd_epi16(t2, f2);
        x1 = _mm256_add_epi32(x1, x2);
        t1 = _mm256_add_epi32(t1, t2);
        t1 = _mm256_srai_epi32(t1, 2);
        x1 = _mm256_srai_epi32(x1, 2);
        r1 = _mm256_packs_epi32(x1, t1);
        src += srcstride;
        x1 = _mm256_loadu_si256((__m256i*)&src[x]);
        x2 = _mm256_loadu_si256((__m256i*)&src[x + 1]);
        x3 = _mm256_loadu_si256((__m256i*)&src[x + 2 * 1]);
        x4 = _mm256_loadu_si256((__m256i*)&src[x + 3 * 1]);
        t1 = _mm256_unpackhi_epi16(x1, x2);
        x1 = _mm256_unpacklo_epi16(x1, x2);
        t2 = _mm256_unpackhi_epi16(x3, x4);
        x2 = _mm256_unpacklo_epi16(x3, x4);
        x1 = _mm256_madd_epi16(x1, f1);
        t1 = _mm256_madd_epi16(t1, f1);
        x2 = _mm256_madd_epi16(x2, f2);
        t2 = _mm256_madd_epi16(t2, f2);
        x1 = _mm256_add_epi32(x1, x2);
        t1 = _mm256_add_epi32(t1, t2);
        t1 = _mm256_srai_epi32(t1, 2);
        x1 = _mm256_srai_epi32(x1, 2);
        r2 = _mm256_packs_epi32(x1, t1);
        src += srcstride;
        x1 = _mm256_loadu_si256((__m256i*)&src[x]);
        x2 = _mm256_loadu_si256((__m256i*)&src[x + 1]);
        x3 = _mm256_loadu_si256((__m256i*)&src[x + 2 * 1]);
        x4 = _mm256_loadu_si256((__m256i*)&src[x + 3 * 1]);
        t1 = _mm256_unpackhi_epi16(x1, x2);
        x1 = _mm256_unpacklo_epi16(x1, x2);
        t2 = _mm256_unpackhi_epi16(x3, x4);
        x2 = _mm256_unpacklo_epi16(x3, x4);
        x1 = _mm256_madd_epi16(x1, f1);
        t1 = _mm256_madd_epi16(t1, f1);
        x2 = _mm256_madd_epi16(x2, f2);
        t2 = _mm256_madd_epi16(t2, f2);
        x1 = _mm256_add_epi32(x1, x2);
        t1 = _mm256_add_epi32(t1, t2);
        t1 = _mm256_srai_epi32(t1, 2);
        x1 = _mm256_srai_epi32(x1, 2);
        r3 = _mm256_packs_epi32(x1, t1);
        src += srcstride;
        for (y = 0; y < height; y++) {
            x1 = _mm256_loadu_si256((__m256i*)&src[x]);
            x2 = _mm256_loadu_si256((__m256i*)&src[x + 1]);
            x3 = _mm256_loadu_si256((__m256i*)&src[x + 2 * 1]);
            x4 = _mm256_loadu_si256((__m256i*)&src[x + 3 * 1]);
            t1 = _mm256_unpackhi_epi16(x1, x2);
            x1 = _mm256_unpacklo_epi16(x1, x2);
            t2 = _mm256_unpackhi_epi16(x3, x4);
            x2 = _mm256_unpacklo_epi16(x3, x4);
            x1 = _mm256_madd_epi16(x1, f1);
            t1 = _mm256_madd_epi16(t1, f1);
            x2 = _mm256_madd_epi16(x2, f2);
            t2 = _mm256_madd_epi16(t2, f2);
            x1 = _mm256_add_epi32(x1, x2);
            t1 = _mm256_add_epi32(t1, t2);
            t1 = _mm256_srai_epi32(t1, 2);
            x1 = _mm256_srai_epi32(x1, 2);
            r4 = _mm256_packs_epi32(x1, t1);
            src += srcstride;
            t1 = _mm256_unpackhi_epi16(r1, r2);
            x1 = _mm256_unpacklo_epi16(r1, r2);
            t2 = _mm256_unpackhi_epi16(r3, r4);
            x2 = _mm256_unpacklo_epi16(r3, r4);
            x1 = _mm256_madd_epi16(x1, f3);
            t1 = _mm256_madd_epi16(t1, f3);
            x2 = _mm256_madd_epi16(x2, f4);
            t2 = _mm256_madd_epi16(t2, f4);
            x1 = _mm256_add_epi32(x1, x2);
            t1 = _mm256_add_epi32(t1, t2);
            t1 = _mm256_srai_epi32(t1, 6);
            x1 = _mm256_srai_epi32(x1, 6);
            x1 = _mm256_packs_epi32(x1, t1);
            _mm256_storeu_si256((__m256i*)&dst[x], x1);
            dst += MAX_PB_SIZE;
            r1 = r2;
            r2 = r3;
            r3 = r4;
        }
        src = src_bis;
        dst = dst_bis;
    }
}

static void
oh_hevc_put_hevc_bi1_epel_hv16_10_avx2(uint16_t* _dst,
                                       ptrdiff_t _dststride,
                                       const uint16_t* _src,
                                       ptrdiff_t _srcstride,
                                       const int16_t* src2,
                                       int height,
                                       intptr_t mx,
                                       intptr_t my,
                                       int width)
{
    int x, y;
    const __m256i offset = _mm256_set1_epi16(1 << BITDEPTH);
    uint16_t* dst = (uint16_t*)_dst;
    const int dststride = _dststride;
    const int16_t* src2_bis = src2;
    uint16_t* dst_bis = dst;
    uint16_t *src_bis, *src = ((uint16_t*)_src) - 1;
    ptrdiff_t srcstride = _srcstride;
    src -= 1 * srcstride;
    src_bis = src;
    __m256i f1 = _mm256_set1_epi32(ov_mcp_filters_c[mx - 1][0]);
    __m256i f2 = _mm256_set1_epi32(ov_mcp_filters_c[mx - 1][1]);
    __m256i f3 = _mm256_set1_epi32(ov_mcp_filters_c[my - 1][0]);
    __m256i f4 = _mm256_set1_epi32(ov_mcp_filters_c[my - 1][1]);
    for (x = 0; x < width; x += 16) {
        __m256i x1, x2, x3, x4, t1, t2, r1, r2, r3, r4;
        __m256i r5;
        x1 = _mm256_loadu_si256((__m256i*)&src[x]);
        x2 = _mm256_loadu_si256((__m256i*)&src[x + 1]);
        x3 = _mm256_loadu_si256((__m256i*)&src[x + 2 * 1]);
        x4 = _mm256_loadu_si256((__m256i*)&src[x + 3 * 1]);
        t1 = _mm256_unpackhi_epi16(x1, x2);
        x1 = _mm256_unpacklo_epi16(x1, x2);
        t2 = _mm256_unpackhi_epi16(x3, x4);
        x2 = _mm256_unpacklo_epi16(x3, x4);
        x1 = _mm256_madd_epi16(x1, f1);
        t1 = _mm256_madd_epi16(t1, f1);
        x2 = _mm256_madd_epi16(x2, f2);
        t2 = _mm256_madd_epi16(t2, f2);
        x1 = _mm256_add_epi32(x1, x2);
        t1 = _mm256_add_epi32(t1, t2);
        t1 = _mm256_srai_epi32(t1, 2);
        x1 = _mm256_srai_epi32(x1, 2);
        r1 = _mm256_packs_epi32(x1, t1);
        src += srcstride;
        x1 = _mm256_loadu_si256((__m256i*)&src[x]);
        x2 = _mm256_loadu_si256((__m256i*)&src[x + 1]);
        x3 = _mm256_loadu_si256((__m256i*)&src[x + 2 * 1]);
        x4 = _mm256_loadu_si256((__m256i*)&src[x + 3 * 1]);
        t1 = _mm256_unpackhi_epi16(x1, x2);
        x1 = _mm256_unpacklo_epi16(x1, x2);
        t2 = _mm256_unpackhi_epi16(x3, x4);
        x2 = _mm256_unpacklo_epi16(x3, x4);
        x1 = _mm256_madd_epi16(x1, f1);
        t1 = _mm256_madd_epi16(t1, f1);
        x2 = _mm256_madd_epi16(x2, f2);
        t2 = _mm256_madd_epi16(t2, f2);
        x1 = _mm256_add_epi32(x1, x2);
        t1 = _mm256_add_epi32(t1, t2);
        t1 = _mm256_srai_epi32(t1, 2);
        x1 = _mm256_srai_epi32(x1, 2);
        r2 = _mm256_packs_epi32(x1, t1);
        src += srcstride;
        x1 = _mm256_loadu_si256((__m256i*)&src[x]);
        x2 = _mm256_loadu_si256((__m256i*)&src[x + 1]);
        x3 = _mm256_loadu_si256((__m256i*)&src[x + 2 * 1]);
        x4 = _mm256_loadu_si256((__m256i*)&src[x + 3 * 1]);
        t1 = _mm256_unpackhi_epi16(x1, x2);
        x1 = _mm256_unpacklo_epi16(x1, x2);
        t2 = _mm256_unpackhi_epi16(x3, x4);
        x2 = _mm256_unpacklo_epi16(x3, x4);
        x1 = _mm256_madd_epi16(x1, f1);
        t1 = _mm256_madd_epi16(t1, f1);
        x2 = _mm256_madd_epi16(x2, f2);
        t2 = _mm256_madd_epi16(t2, f2);
        x1 = _mm256_add_epi32(x1, x2);
        t1 = _mm256_add_epi32(t1, t2);
        t1 = _mm256_srai_epi32(t1, 2);
        x1 = _mm256_srai_epi32(x1, 2);
        r3 = _mm256_packs_epi32(x1, t1);
        src += srcstride;
        for (y = 0; y < height; y++) {
            x1 = _mm256_loadu_si256((__m256i*)&src[x]);
            x2 = _mm256_loadu_si256((__m256i*)&src[x + 1]);
            x3 = _mm256_loadu_si256((__m256i*)&src[x + 2 * 1]);
            x4 = _mm256_loadu_si256((__m256i*)&src[x + 3 * 1]);
            t1 = _mm256_unpackhi_epi16(x1, x2);
            x1 = _mm256_unpacklo_epi16(x1, x2);
            t2 = _mm256_unpackhi_epi16(x3, x4);
            x2 = _mm256_unpacklo_epi16(x3, x4);
            x1 = _mm256_madd_epi16(x1, f1);
            t1 = _mm256_madd_epi16(t1, f1);
            x2 = _mm256_madd_epi16(x2, f2);
            t2 = _mm256_madd_epi16(t2, f2);
            x1 = _mm256_add_epi32(x1, x2);
            t1 = _mm256_add_epi32(t1, t2);
            t1 = _mm256_srai_epi32(t1, 2);
            x1 = _mm256_srai_epi32(x1, 2);
            r4 = _mm256_packs_epi32(x1, t1);
            src += srcstride;
            t1 = _mm256_unpackhi_epi16(r1, r2);
            x1 = _mm256_unpacklo_epi16(r1, r2);
            t2 = _mm256_unpackhi_epi16(r3, r4);
            x2 = _mm256_unpacklo_epi16(r3, r4);
            x1 = _mm256_madd_epi16(x1, f3);
            t1 = _mm256_madd_epi16(t1, f3);
            x2 = _mm256_madd_epi16(x2, f4);
            t2 = _mm256_madd_epi16(t2, f4);
            x1 = _mm256_add_epi32(x1, x2);
            t1 = _mm256_add_epi32(t1, t2);
            t1 = _mm256_srai_epi32(t1, 6);
            x1 = _mm256_srai_epi32(x1, 6);
            x1 = _mm256_packs_epi32(x1, t1);
            r5 = _mm256_load_si256((__m256i*)&src2[x]);
            x1 = _mm256_adds_epi16(x1, r5);
            x1 = _mm256_mulhrs_epi16(x1, offset);
            x1 = _mm256_max_epi16(x1, _mm256_setzero_si256());
            x1 = _mm256_min_epi16(x1, _mm256_set1_epi16(0x03FF));
            _mm256_storeu_si256((__m256i*)&dst[x], x1);
            src2 += MAX_PB_SIZE;
            dst += dststride;
            r1 = r2;
            r2 = r3;
            r3 = r4;
        }
        src = src_bis;
        src2 = src2_bis;
        dst = dst_bis;
    }
}

static void
oh_hevc_put_hevc_uni_epel_hv16_10_avx2(uint16_t* _dst,
                                       ptrdiff_t _dststride,
                                       const uint16_t* _src,
                                       ptrdiff_t _srcstride,
                                       int height,
                                       intptr_t mx,
                                       intptr_t my,
                                       int width)
{
    int x, y;
    const __m256i offset = _mm256_set1_epi16(1 << (BITDEPTH + 1));
    uint16_t* dst = (uint16_t*)_dst;
    const int dststride = _dststride;
    uint16_t* dst_bis = dst;
    uint16_t *src_bis, *src = ((uint16_t*)_src) - 1;
    ptrdiff_t srcstride = _srcstride;
    src -= 1 * srcstride;
    src_bis = src;
    __m256i f1 = _mm256_set1_epi32(ov_mcp_filters_c[mx - 1][0]);
    __m256i f2 = _mm256_set1_epi32(ov_mcp_filters_c[mx - 1][1]);
    __m256i f3 = _mm256_set1_epi32(ov_mcp_filters_c[my - 1][0]);
    __m256i f4 = _mm256_set1_epi32(ov_mcp_filters_c[my - 1][1]);
    for (x = 0; x < width; x += 16) {
        __m256i x1, x2, x3, x4, t1, t2, r1, r2, r3, r4;
        x1 = _mm256_loadu_si256((__m256i*)&src[x]);
        x2 = _mm256_loadu_si256((__m256i*)&src[x + 1]);
        x3 = _mm256_loadu_si256((__m256i*)&src[x + 2 * 1]);
        x4 = _mm256_loadu_si256((__m256i*)&src[x + 3 * 1]);
        t1 = _mm256_unpackhi_epi16(x1, x2);
        x1 = _mm256_unpacklo_epi16(x1, x2);
        t2 = _mm256_unpackhi_epi16(x3, x4);
        x2 = _mm256_unpacklo_epi16(x3, x4);
        x1 = _mm256_madd_epi16(x1, f1);
        t1 = _mm256_madd_epi16(t1, f1);
        x2 = _mm256_madd_epi16(x2, f2);
        t2 = _mm256_madd_epi16(t2, f2);
        x1 = _mm256_add_epi32(x1, x2);
        t1 = _mm256_add_epi32(t1, t2);
        t1 = _mm256_srai_epi32(t1, 2);
        x1 = _mm256_srai_epi32(x1, 2);
        r1 = _mm256_packs_epi32(x1, t1);
        src += srcstride;
        x1 = _mm256_loadu_si256((__m256i*)&src[x]);
        x2 = _mm256_loadu_si256((__m256i*)&src[x + 1]);
        x3 = _mm256_loadu_si256((__m256i*)&src[x + 2 * 1]);
        x4 = _mm256_loadu_si256((__m256i*)&src[x + 3 * 1]);
        t1 = _mm256_unpackhi_epi16(x1, x2);
        x1 = _mm256_unpacklo_epi16(x1, x2);
        t2 = _mm256_unpackhi_epi16(x3, x4);
        x2 = _mm256_unpacklo_epi16(x3, x4);
        x1 = _mm256_madd_epi16(x1, f1);
        t1 = _mm256_madd_epi16(t1, f1);
        x2 = _mm256_madd_epi16(x2, f2);
        t2 = _mm256_madd_epi16(t2, f2);
        x1 = _mm256_add_epi32(x1, x2);
        t1 = _mm256_add_epi32(t1, t2);
        t1 = _mm256_srai_epi32(t1, 2);
        x1 = _mm256_srai_epi32(x1, 2);
        r2 = _mm256_packs_epi32(x1, t1);
        src += srcstride;
        x1 = _mm256_loadu_si256((__m256i*)&src[x]);
        x2 = _mm256_loadu_si256((__m256i*)&src[x + 1]);
        x3 = _mm256_loadu_si256((__m256i*)&src[x + 2 * 1]);
        x4 = _mm256_loadu_si256((__m256i*)&src[x + 3 * 1]);
        t1 = _mm256_unpackhi_epi16(x1, x2);
        x1 = _mm256_unpacklo_epi16(x1, x2);
        t2 = _mm256_unpackhi_epi16(x3, x4);
        x2 = _mm256_unpacklo_epi16(x3, x4);
        x1 = _mm256_madd_epi16(x1, f1);
        t1 = _mm256_madd_epi16(t1, f1);
        x2 = _mm256_madd_epi16(x2, f2);
        t2 = _mm256_madd_epi16(t2, f2);
        x1 = _mm256_add_epi32(x1, x2);
        t1 = _mm256_add_epi32(t1, t2);
        t1 = _mm256_srai_epi32(t1, 2);
        x1 = _mm256_srai_epi32(x1, 2);
        r3 = _mm256_packs_epi32(x1, t1);
        src += srcstride;
        for (y = 0; y < height; y++) {
            x1 = _mm256_loadu_si256((__m256i*)&src[x]);
            x2 = _mm256_loadu_si256((__m256i*)&src[x + 1]);
            x3 = _mm256_loadu_si256((__m256i*)&src[x + 2 * 1]);
            x4 = _mm256_loadu_si256((__m256i*)&src[x + 3 * 1]);
            t1 = _mm256_unpackhi_epi16(x1, x2);
            x1 = _mm256_unpacklo_epi16(x1, x2);
            t2 = _mm256_unpackhi_epi16(x3, x4);
            x2 = _mm256_unpacklo_epi16(x3, x4);
            x1 = _mm256_madd_epi16(x1, f1);
            t1 = _mm256_madd_epi16(t1, f1);
            x2 = _mm256_madd_epi16(x2, f2);
            t2 = _mm256_madd_epi16(t2, f2);
            x1 = _mm256_add_epi32(x1, x2);
            t1 = _mm256_add_epi32(t1, t2);
            t1 = _mm256_srai_epi32(t1, 2);
            x1 = _mm256_srai_epi32(x1, 2);
            r4 = _mm256_packs_epi32(x1, t1);
            src += srcstride;
            t1 = _mm256_unpackhi_epi16(r1, r2);
            x1 = _mm256_unpacklo_epi16(r1, r2);
            t2 = _mm256_unpackhi_epi16(r3, r4);
            x2 = _mm256_unpacklo_epi16(r3, r4);
            x1 = _mm256_madd_epi16(x1, f3);
            t1 = _mm256_madd_epi16(t1, f3);
            x2 = _mm256_madd_epi16(x2, f4);
            t2 = _mm256_madd_epi16(t2, f4);
            x1 = _mm256_add_epi32(x1, x2);
            t1 = _mm256_add_epi32(t1, t2);
            t1 = _mm256_srai_epi32(t1, 6);
            x1 = _mm256_srai_epi32(x1, 6);
            x1 = _mm256_packs_epi32(x1, t1);
            x1 = _mm256_mulhrs_epi16(x1, offset);
            x1 = _mm256_max_epi16(x1, _mm256_setzero_si256());
            x1 = _mm256_min_epi16(x1, _mm256_set1_epi16(0x03FF));
            _mm256_storeu_si256((__m256i*)&dst[x], x1);
            dst += dststride;
            r1 = r2;
            r2 = r3;
            r3 = r4;
        }
        src = src_bis;
        dst = dst_bis;
    }
}

static void
oh_hevc_put_hevc_bi0_qpel_h16_10_avx2(int16_t* dst,
                                      const uint16_t* _src,
                                      ptrdiff_t _srcstride,
                                      int height,
                                      intptr_t mx,
                                      intptr_t my,
                                      int width)
{
    int x, y;
    int shift = BITDEPTH - 8;
    const uint16_t* src = _src;
    const int srcstride = _srcstride;

    __m256i c1 = _mm256_set1_epi32(ov_mcp_filters_l[mx - 1][0]);
    __m256i c2 = _mm256_set1_epi32(ov_mcp_filters_l[mx - 1][1]);
    __m256i c3 = _mm256_set1_epi32(ov_mcp_filters_l[mx - 1][2]);
    __m256i c4 = _mm256_set1_epi32(ov_mcp_filters_l[mx - 1][3]);

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x += 16) {
            __m256i x1, x2, x3, x4, r1, r2, r3, r4, t1, t2;
            x1 = _mm256_loadu_si256((__m256i*)&src[x - 3 * 1]);
            x2 = _mm256_loadu_si256((__m256i*)&src[x - 2 * 1]);
            x3 = _mm256_loadu_si256((__m256i*)&src[x - 1 * 1]);
            x4 = _mm256_loadu_si256((__m256i*)&src[x]);
            t1 = _mm256_unpackhi_epi16(x1, x2);
            x1 = _mm256_unpacklo_epi16(x1, x2);
            t2 = _mm256_unpackhi_epi16(x3, x4);
            x2 = _mm256_unpacklo_epi16(x3, x4);
            t1 = _mm256_madd_epi16(t1, c1);
            t2 = _mm256_madd_epi16(t2, c2);
            r2 = _mm256_add_epi32(t1, t2);
            x1 = _mm256_madd_epi16(x1, c1);
            x2 = _mm256_madd_epi16(x2, c2);
            r1 = _mm256_add_epi32(x1, x2);
            x1 = _mm256_loadu_si256((__m256i*)&src[x + 1]);
            x2 = _mm256_loadu_si256((__m256i*)&src[x + 2 * 1]);
            x3 = _mm256_loadu_si256((__m256i*)&src[x + 3 * 1]);
            x4 = _mm256_loadu_si256((__m256i*)&src[x + 4 * 1]);
            t1 = _mm256_unpackhi_epi16(x1, x2);
            x1 = _mm256_unpacklo_epi16(x1, x2);
            t2 = _mm256_unpackhi_epi16(x3, x4);
            x2 = _mm256_unpacklo_epi16(x3, x4);
            t1 = _mm256_madd_epi16(t1, c3);
            t2 = _mm256_madd_epi16(t2, c4);
            r4 = _mm256_add_epi32(t1, t2);
            x1 = _mm256_madd_epi16(x1, c3);
            x2 = _mm256_madd_epi16(x2, c4);
            r3 = _mm256_add_epi32(x1, x2);
            x1 = _mm256_add_epi32(r1, r3);
            x2 = _mm256_add_epi32(r2, r4);
            x1 = _mm256_srai_epi32(x1, shift);
            x2 = _mm256_srai_epi32(x2, shift);
            x1 = _mm256_packs_epi32(x1, x2);
            _mm256_storeu_si256((__m256i*)&dst[x], x1);
        }
        src += srcstride;
        dst += MAX_PB_SIZE;
    }
}

static void
oh_hevc_put_hevc_bi1_qpel_h16_10_avx2(uint16_t* _dst,
                                      ptrdiff_t _dststride,
                                      const uint16_t* _src,
                                      ptrdiff_t _srcstride,
                                      const int16_t* src2,
                                      int height,
                                      intptr_t mx,
                                      intptr_t my,
                                      int width)
{
    int x, y;
    int shift = BITDEPTH - 8;
    const uint16_t* src = _src;
    const int srcstride = _srcstride;
    const __m256i offset = _mm256_set1_epi16(1 << BITDEPTH);
    uint16_t* dst = (uint16_t*)_dst;
    const int dststride = _dststride;

    __m256i c1 = _mm256_set1_epi32(ov_mcp_filters_l[mx - 1][0]);
    __m256i c2 = _mm256_set1_epi32(ov_mcp_filters_l[mx - 1][1]);
    __m256i c3 = _mm256_set1_epi32(ov_mcp_filters_l[mx - 1][2]);
    __m256i c4 = _mm256_set1_epi32(ov_mcp_filters_l[mx - 1][3]);

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x += 16) {
            __m256i x1, x2, x3, x4, r1, r2, r3, r4, t1, t2;
            __m256i r5;
            x1 = _mm256_loadu_si256((__m256i*)&src[x - 3 * 1]);
            x2 = _mm256_loadu_si256((__m256i*)&src[x - 2 * 1]);
            x3 = _mm256_loadu_si256((__m256i*)&src[x - 1 * 1]);
            x4 = _mm256_loadu_si256((__m256i*)&src[x]);
            t1 = _mm256_unpackhi_epi16(x1, x2);
            x1 = _mm256_unpacklo_epi16(x1, x2);
            t2 = _mm256_unpackhi_epi16(x3, x4);
            x2 = _mm256_unpacklo_epi16(x3, x4);
            t1 = _mm256_madd_epi16(t1, c1);
            t2 = _mm256_madd_epi16(t2, c2);
            r2 = _mm256_add_epi32(t1, t2);
            x1 = _mm256_madd_epi16(x1, c1);
            x2 = _mm256_madd_epi16(x2, c2);
            r1 = _mm256_add_epi32(x1, x2);
            x1 = _mm256_loadu_si256((__m256i*)&src[x + 1]);
            x2 = _mm256_loadu_si256((__m256i*)&src[x + 2 * 1]);
            x3 = _mm256_loadu_si256((__m256i*)&src[x + 3 * 1]);
            x4 = _mm256_loadu_si256((__m256i*)&src[x + 4 * 1]);
            t1 = _mm256_unpackhi_epi16(x1, x2);
            x1 = _mm256_unpacklo_epi16(x1, x2);
            t2 = _mm256_unpackhi_epi16(x3, x4);
            x2 = _mm256_unpacklo_epi16(x3, x4);
            t1 = _mm256_madd_epi16(t1, c3);
            t2 = _mm256_madd_epi16(t2, c4);
            r4 = _mm256_add_epi32(t1, t2);
            x1 = _mm256_madd_epi16(x1, c3);
            x2 = _mm256_madd_epi16(x2, c4);
            r3 = _mm256_add_epi32(x1, x2);
            x1 = _mm256_add_epi32(r1, r3);
            x2 = _mm256_add_epi32(r2, r4);
            x1 = _mm256_srai_epi32(x1, shift);
            x2 = _mm256_srai_epi32(x2, shift);
            x1 = _mm256_packs_epi32(x1, x2);
            r5 = _mm256_load_si256((__m256i*)&src2[x]);
            x1 = _mm256_adds_epi16(x1, r5);
            x1 = _mm256_mulhrs_epi16(x1, offset);
            x1 = _mm256_max_epi16(x1, _mm256_setzero_si256());
            x1 = _mm256_min_epi16(x1, _mm256_set1_epi16(0x03FF));
            _mm256_storeu_si256((__m256i*)&dst[x], x1);
        }
        src += srcstride;
        src2 += MAX_PB_SIZE;
        dst += dststride;
    }
}

static void
oh_hevc_put_hevc_uni_qpel_h16_10_avx2(uint16_t* _dst,
                                      ptrdiff_t _dststride,
                                      const uint16_t* _src,
                                      ptrdiff_t _srcstride,
                                      int height,
                                      intptr_t mx,
                                      intptr_t my,
                                      int width)
{
    int x, y;
    int shift = BITDEPTH - 8;
    const uint16_t* src = _src;
    const int srcstride = _srcstride;
    const __m256i offset = _mm256_set1_epi16(1 << (BITDEPTH + 1));
    uint16_t* dst = (uint16_t*)_dst;
    const int dststride = _dststride;

    __m256i c1 = _mm256_set1_epi32(ov_mcp_filters_l[mx - 1][0]);
    __m256i c2 = _mm256_set1_epi32(ov_mcp_filters_l[mx - 1][1]);
    __m256i c3 = _mm256_set1_epi32(ov_mcp_filters_l[mx - 1][2]);
    __m256i c4 = _mm256_set1_epi32(ov_mcp_filters_l[mx - 1][3]);

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x += 16) {
            __m256i x1, x2, x3, x4, r1, r2, r3, r4, t1, t2;
            x1 = _mm256_loadu_si256((__m256i*)&src[x - 3 * 1]);
            x2 = _mm256_loadu_si256((__m256i*)&src[x - 2 * 1]);
            x3 = _mm256_loadu_si256((__m256i*)&src[x - 1 * 1]);
            x4 = _mm256_loadu_si256((__m256i*)&src[x]);
            t1 = _mm256_unpackhi_epi16(x1, x2);
            x1 = _mm256_unpacklo_epi16(x1, x2);
            t2 = _mm256_unpackhi_epi16(x3, x4);
            x2 = _mm256_unpacklo_epi16(x3, x4);
            t1 = _mm256_madd_epi16(t1, c1);
            t2 = _mm256_madd_epi16(t2, c2);
            r2 = _mm256_add_epi32(t1, t2);
            x1 = _mm256_madd_epi16(x1, c1);
            x2 = _mm256_madd_epi16(x2, c2);
            r1 = _mm256_add_epi32(x1, x2);
            x1 = _mm256_loadu_si256((__m256i*)&src[x + 1]);
            x2 = _mm256_loadu_si256((__m256i*)&src[x + 2 * 1]);
            x3 = _mm256_loadu_si256((__m256i*)&src[x + 3 * 1]);
            x4 = _mm256_loadu_si256((__m256i*)&src[x + 4 * 1]);
            t1 = _mm256_unpackhi_epi16(x1, x2);
            x1 = _mm256_unpacklo_epi16(x1, x2);
            t2 = _mm256_unpackhi_epi16(x3, x4);
            x2 = _mm256_unpacklo_epi16(x3, x4);
            t1 = _mm256_madd_epi16(t1, c3);
            t2 = _mm256_madd_epi16(t2, c4);
            r4 = _mm256_add_epi32(t1, t2);
            x1 = _mm256_madd_epi16(x1, c3);
            x2 = _mm256_madd_epi16(x2, c4);
            r3 = _mm256_add_epi32(x1, x2);
            x1 = _mm256_add_epi32(r1, r3);
            x2 = _mm256_add_epi32(r2, r4);
            x1 = _mm256_srai_epi32(x1, shift);
            x2 = _mm256_srai_epi32(x2, shift);
            x1 = _mm256_packs_epi32(x1, x2);
            x1 = _mm256_mulhrs_epi16(x1, offset);
            x1 = _mm256_max_epi16(x1, _mm256_setzero_si256());
            x1 = _mm256_min_epi16(x1, _mm256_set1_epi16(0x03FF));
            _mm256_storeu_si256((__m256i*)&dst[x], x1);
        }
        src += srcstride;
        dst += dststride;
    }
}

static void
oh_hevc_put_hevc_bi0_qpel_v16_10_avx2(int16_t* dst,
                                      const uint16_t* _src,
                                      ptrdiff_t _srcstride,
                                      int height,
                                      intptr_t mx,
                                      intptr_t my,
                                      int width)
{
    int x, y;
    const uint16_t* src = _src;
    const int srcstride = _srcstride;

    __m256i c1 = _mm256_set1_epi32(ov_mcp_filters_l[my - 1][0]);
    __m256i c2 = _mm256_set1_epi32(ov_mcp_filters_l[my - 1][1]);
    __m256i c3 = _mm256_set1_epi32(ov_mcp_filters_l[my - 1][2]);
    __m256i c4 = _mm256_set1_epi32(ov_mcp_filters_l[my - 1][3]);

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x += 16) {
            __m256i x1, x2, x3, x4, x5, x6, x7, x8, x9;
            x1 = _mm256_loadu_si256((__m256i*)&src[x - 3 * srcstride]);
            x2 = _mm256_loadu_si256((__m256i*)&src[x - 2 * srcstride]);
            x3 = _mm256_loadu_si256((__m256i*)&src[x - 1 * srcstride]);
            x4 = _mm256_loadu_si256((__m256i*)&src[x]);
            x5 = _mm256_loadu_si256((__m256i*)&src[x + srcstride]);
            x6 = _mm256_loadu_si256((__m256i*)&src[x + 2 * srcstride]);
            x7 = _mm256_loadu_si256((__m256i*)&src[x + 3 * srcstride]);
            x8 = _mm256_loadu_si256((__m256i*)&src[x + 4 * srcstride]);
            x9 = x1;
            x1 = _mm256_unpacklo_epi16(x9, x2);
            x2 = _mm256_unpackhi_epi16(x9, x2);
            x9 = x3;
            x3 = _mm256_unpacklo_epi16(x9, x4);
            x4 = _mm256_unpackhi_epi16(x9, x4);
            x9 = x5;
            x5 = _mm256_unpacklo_epi16(x9, x6);
            x6 = _mm256_unpackhi_epi16(x9, x6);
            x9 = x7;
            x7 = _mm256_unpacklo_epi16(x9, x8);
            x8 = _mm256_unpackhi_epi16(x9, x8);
            x1 = _mm256_madd_epi16(x1, c1);
            x3 = _mm256_madd_epi16(x3, c2);
            x5 = _mm256_madd_epi16(x5, c3);
            x7 = _mm256_madd_epi16(x7, c4);
            x2 = _mm256_madd_epi16(x2, c1);
            x4 = _mm256_madd_epi16(x4, c2);
            x6 = _mm256_madd_epi16(x6, c3);
            x8 = _mm256_madd_epi16(x8, c4);
            x1 = _mm256_add_epi32(x1, x3);
            x3 = _mm256_add_epi32(x5, x7);
            x2 = _mm256_add_epi32(x2, x4);
            x4 = _mm256_add_epi32(x6, x8);
            x1 = _mm256_add_epi32(x1, x3);
            x2 = _mm256_add_epi32(x2, x4);
            x1 = _mm256_srai_epi32(x1, 2);
            x2 = _mm256_srai_epi32(x2, 2);
            x1 = _mm256_packs_epi32(x1, x2);
            _mm256_storeu_si256((__m256i*)&dst[x], x1);
        }
        src += srcstride;
        dst += MAX_PB_SIZE;
    }
}

static void
oh_hevc_put_hevc_bi1_qpel_v16_10_avx2(uint16_t* _dst,
                                      ptrdiff_t _dststride,
                                      const uint16_t* _src,
                                      ptrdiff_t _srcstride,
                                      const int16_t* src2,
                                      int height,
                                      intptr_t mx,
                                      intptr_t my,
                                      int width)
{
    int x, y;
    const uint16_t* src = _src;
    const int srcstride = _srcstride;
    const __m256i offset = _mm256_set1_epi16(1 << BITDEPTH);
    uint16_t* dst = (uint16_t*)_dst;
    const int dststride = _dststride;

    __m256i c1 = _mm256_set1_epi32(ov_mcp_filters_l[my - 1][0]);
    __m256i c2 = _mm256_set1_epi32(ov_mcp_filters_l[my - 1][1]);
    __m256i c3 = _mm256_set1_epi32(ov_mcp_filters_l[my - 1][2]);
    __m256i c4 = _mm256_set1_epi32(ov_mcp_filters_l[my - 1][3]);

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x += 16) {
            __m256i x1, x2, x3, x4, x5, x6, x7, x8, x9;
            __m256i r5;
            x1 = _mm256_loadu_si256((__m256i*)&src[x - 3 * srcstride]);
            x2 = _mm256_loadu_si256((__m256i*)&src[x - 2 * srcstride]);
            x3 = _mm256_loadu_si256((__m256i*)&src[x - 1 * srcstride]);
            x4 = _mm256_loadu_si256((__m256i*)&src[x]);
            x5 = _mm256_loadu_si256((__m256i*)&src[x + srcstride]);
            x6 = _mm256_loadu_si256((__m256i*)&src[x + 2 * srcstride]);
            x7 = _mm256_loadu_si256((__m256i*)&src[x + 3 * srcstride]);
            x8 = _mm256_loadu_si256((__m256i*)&src[x + 4 * srcstride]);
            x9 = x1;
            x1 = _mm256_unpacklo_epi16(x9, x2);
            x2 = _mm256_unpackhi_epi16(x9, x2);
            x9 = x3;
            x3 = _mm256_unpacklo_epi16(x9, x4);
            x4 = _mm256_unpackhi_epi16(x9, x4);
            x9 = x5;
            x5 = _mm256_unpacklo_epi16(x9, x6);
            x6 = _mm256_unpackhi_epi16(x9, x6);
            x9 = x7;
            x7 = _mm256_unpacklo_epi16(x9, x8);
            x8 = _mm256_unpackhi_epi16(x9, x8);
            x1 = _mm256_madd_epi16(x1, c1);
            x3 = _mm256_madd_epi16(x3, c2);
            x5 = _mm256_madd_epi16(x5, c3);
            x7 = _mm256_madd_epi16(x7, c4);
            x2 = _mm256_madd_epi16(x2, c1);
            x4 = _mm256_madd_epi16(x4, c2);
            x6 = _mm256_madd_epi16(x6, c3);
            x8 = _mm256_madd_epi16(x8, c4);
            x1 = _mm256_add_epi32(x1, x3);
            x3 = _mm256_add_epi32(x5, x7);
            x2 = _mm256_add_epi32(x2, x4);
            x4 = _mm256_add_epi32(x6, x8);
            x1 = _mm256_add_epi32(x1, x3);
            x2 = _mm256_add_epi32(x2, x4);
            x1 = _mm256_srai_epi32(x1, 2);
            x2 = _mm256_srai_epi32(x2, 2);
            x1 = _mm256_packs_epi32(x1, x2);
            r5 = _mm256_load_si256((__m256i*)&src2[x]);
            x1 = _mm256_adds_epi16(x1, r5);
            x1 = _mm256_mulhrs_epi16(x1, offset);
            x1 = _mm256_max_epi16(x1, _mm256_setzero_si256());
            x1 = _mm256_min_epi16(x1, _mm256_set1_epi16(0x03FF));
            _mm256_storeu_si256((__m256i*)&dst[x], x1);
        }
        src += srcstride;
        src2 += MAX_PB_SIZE;
        dst += dststride;
    }
}

static void
oh_hevc_put_hevc_uni_qpel_v16_10_avx2(uint16_t* _dst,
                                      ptrdiff_t _dststride,
                                      const uint16_t* _src,
                                      ptrdiff_t _srcstride,
                                      int height,
                                      intptr_t mx,
                                      intptr_t my,
                                      int width)
{
    int x, y;
    const uint16_t* src = _src;
    const int srcstride = _srcstride;
    const __m256i offset = _mm256_set1_epi16(1 << (BITDEPTH + 1));
    uint16_t* dst = (uint16_t*)_dst;
    const int dststride = _dststride;

    __m256i c1 = _mm256_set1_epi32(ov_mcp_filters_l[my - 1][0]);
    __m256i c2 = _mm256_set1_epi32(ov_mcp_filters_l[my - 1][1]);
    __m256i c3 = _mm256_set1_epi32(ov_mcp_filters_l[my - 1][2]);
    __m256i c4 = _mm256_set1_epi32(ov_mcp_filters_l[my - 1][3]);

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x += 16) {
            __m256i x1, x2, x3, x4, x5, x6, x7, x8, x9;
            x1 = _mm256_loadu_si256((__m256i*)&src[x - 3 * srcstride]);
            x2 = _mm256_loadu_si256((__m256i*)&src[x - 2 * srcstride]);
            x3 = _mm256_loadu_si256((__m256i*)&src[x - 1 * srcstride]);
            x4 = _mm256_loadu_si256((__m256i*)&src[x]);
            x5 = _mm256_loadu_si256((__m256i*)&src[x + srcstride]);
            x6 = _mm256_loadu_si256((__m256i*)&src[x + 2 * srcstride]);
            x7 = _mm256_loadu_si256((__m256i*)&src[x + 3 * srcstride]);
            x8 = _mm256_loadu_si256((__m256i*)&src[x + 4 * srcstride]);
            x9 = x1;
            x1 = _mm256_unpacklo_epi16(x9, x2);
            x2 = _mm256_unpackhi_epi16(x9, x2);
            x9 = x3;
            x3 = _mm256_unpacklo_epi16(x9, x4);
            x4 = _mm256_unpackhi_epi16(x9, x4);
            x9 = x5;
            x5 = _mm256_unpacklo_epi16(x9, x6);
            x6 = _mm256_unpackhi_epi16(x9, x6);
            x9 = x7;
            x7 = _mm256_unpacklo_epi16(x9, x8);
            x8 = _mm256_unpackhi_epi16(x9, x8);
            x1 = _mm256_madd_epi16(x1, c1);
            x3 = _mm256_madd_epi16(x3, c2);
            x5 = _mm256_madd_epi16(x5, c3);
            x7 = _mm256_madd_epi16(x7, c4);
            x2 = _mm256_madd_epi16(x2, c1);
            x4 = _mm256_madd_epi16(x4, c2);
            x6 = _mm256_madd_epi16(x6, c3);
            x8 = _mm256_madd_epi16(x8, c4);
            x1 = _mm256_add_epi32(x1, x3);
            x3 = _mm256_add_epi32(x5, x7);
            x2 = _mm256_add_epi32(x2, x4);
            x4 = _mm256_add_epi32(x6, x8);
            x1 = _mm256_add_epi32(x1, x3);
            x2 = _mm256_add_epi32(x2, x4);
            x1 = _mm256_srai_epi32(x1, 2);
            x2 = _mm256_srai_epi32(x2, 2);
            x1 = _mm256_packs_epi32(x1, x2);
            x1 = _mm256_mulhrs_epi16(x1, offset);
            x1 = _mm256_max_epi16(x1, _mm256_setzero_si256());
            x1 = _mm256_min_epi16(x1, _mm256_set1_epi16(0x03FF));
            _mm256_storeu_si256((__m256i*)&dst[x], x1);
        }
        src += srcstride;
        dst += dststride;
    }
}

static void
oh_hevc_put_hevc_bi0_qpel_v16_14_avx2(int16_t* dst,
                                      const uint16_t* _src,
                                      ptrdiff_t _srcstride,
                                      int height,
                                      intptr_t mx,
                                      intptr_t my,
                                      int width)
{
    int x, y;
    const uint16_t* src = _src;
    const int srcstride = _srcstride;

    __m256i c1 = _mm256_set1_epi32(ov_mcp_filters_l[my - 1][0]);
    __m256i c2 = _mm256_set1_epi32(ov_mcp_filters_l[my - 1][1]);
    __m256i c3 = _mm256_set1_epi32(ov_mcp_filters_l[my - 1][2]);
    __m256i c4 = _mm256_set1_epi32(ov_mcp_filters_l[my - 1][3]);

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x += 16) {
            __m256i x1, x2, x3, x4, x5, x6, x7, x8, x9;
            x1 = _mm256_loadu_si256((__m256i*)&src[x - 3 * srcstride]);
            x2 = _mm256_loadu_si256((__m256i*)&src[x - 2 * srcstride]);
            x3 = _mm256_loadu_si256((__m256i*)&src[x - 1 * srcstride]);
            x4 = _mm256_loadu_si256((__m256i*)&src[x]);
            x5 = _mm256_loadu_si256((__m256i*)&src[x + srcstride]);
            x6 = _mm256_loadu_si256((__m256i*)&src[x + 2 * srcstride]);
            x7 = _mm256_loadu_si256((__m256i*)&src[x + 3 * srcstride]);
            x8 = _mm256_loadu_si256((__m256i*)&src[x + 4 * srcstride]);
            x9 = x1;
            x1 = _mm256_unpacklo_epi16(x9, x2);
            x2 = _mm256_unpackhi_epi16(x9, x2);
            x9 = x3;
            x3 = _mm256_unpacklo_epi16(x9, x4);
            x4 = _mm256_unpackhi_epi16(x9, x4);
            x9 = x5;
            x5 = _mm256_unpacklo_epi16(x9, x6);
            x6 = _mm256_unpackhi_epi16(x9, x6);
            x9 = x7;
            x7 = _mm256_unpacklo_epi16(x9, x8);
            x8 = _mm256_unpackhi_epi16(x9, x8);
            x1 = _mm256_madd_epi16(x1, c1);
            x3 = _mm256_madd_epi16(x3, c2);
            x5 = _mm256_madd_epi16(x5, c3);
            x7 = _mm256_madd_epi16(x7, c4);
            x2 = _mm256_madd_epi16(x2, c1);
            x4 = _mm256_madd_epi16(x4, c2);
            x6 = _mm256_madd_epi16(x6, c3);
            x8 = _mm256_madd_epi16(x8, c4);
            x1 = _mm256_add_epi32(x1, x3);
            x3 = _mm256_add_epi32(x5, x7);
            x2 = _mm256_add_epi32(x2, x4);
            x4 = _mm256_add_epi32(x6, x8);
            x1 = _mm256_add_epi32(x1, x3);
            x2 = _mm256_add_epi32(x2, x4);
            x1 = _mm256_srai_epi32(x1, 6);
            x2 = _mm256_srai_epi32(x2, 6);
            x1 = _mm256_packs_epi32(x1, x2);
            _mm256_storeu_si256((__m256i*)&dst[x], x1);
        }
        src += srcstride;
        dst += MAX_PB_SIZE;
    }
}

static void
oh_hevc_put_hevc_bi1_qpel_v16_14_10_avx2(uint16_t* _dst,
                                         ptrdiff_t _dststride,
                                         const int16_t* _src,
                                         ptrdiff_t _srcstride,
                                         const int16_t* src2,
                                         int height,
                                         intptr_t mx,
                                         intptr_t my,
                                         int width)
{
    int x, y;
    const int16_t* src = _src;
    const int srcstride = _srcstride;
    const __m256i offset = _mm256_set1_epi16(1 << BITDEPTH);
    uint16_t* dst = (uint16_t*)_dst;
    const int dststride = _dststride;

    __m256i c1 = _mm256_set1_epi32(ov_mcp_filters_l[my - 1][0]);
    __m256i c2 = _mm256_set1_epi32(ov_mcp_filters_l[my - 1][1]);
    __m256i c3 = _mm256_set1_epi32(ov_mcp_filters_l[my - 1][2]);
    __m256i c4 = _mm256_set1_epi32(ov_mcp_filters_l[my - 1][3]);

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x += 16) {
            __m256i x1, x2, x3, x4, x5, x6, x7, x8, x9;
            __m256i r5;
            x1 = _mm256_loadu_si256((__m256i*)&src[x - 3 * srcstride]);
            x2 = _mm256_loadu_si256((__m256i*)&src[x - 2 * srcstride]);
            x3 = _mm256_loadu_si256((__m256i*)&src[x - 1 * srcstride]);
            x4 = _mm256_loadu_si256((__m256i*)&src[x]);
            x5 = _mm256_loadu_si256((__m256i*)&src[x + srcstride]);
            x6 = _mm256_loadu_si256((__m256i*)&src[x + 2 * srcstride]);
            x7 = _mm256_loadu_si256((__m256i*)&src[x + 3 * srcstride]);
            x8 = _mm256_loadu_si256((__m256i*)&src[x + 4 * srcstride]);
            x9 = x1;
            x1 = _mm256_unpacklo_epi16(x9, x2);
            x2 = _mm256_unpackhi_epi16(x9, x2);
            x9 = x3;
            x3 = _mm256_unpacklo_epi16(x9, x4);
            x4 = _mm256_unpackhi_epi16(x9, x4);
            x9 = x5;
            x5 = _mm256_unpacklo_epi16(x9, x6);
            x6 = _mm256_unpackhi_epi16(x9, x6);
            x9 = x7;
            x7 = _mm256_unpacklo_epi16(x9, x8);
            x8 = _mm256_unpackhi_epi16(x9, x8);
            x1 = _mm256_madd_epi16(x1, c1);
            x3 = _mm256_madd_epi16(x3, c2);
            x5 = _mm256_madd_epi16(x5, c3);
            x7 = _mm256_madd_epi16(x7, c4);
            x2 = _mm256_madd_epi16(x2, c1);
            x4 = _mm256_madd_epi16(x4, c2);
            x6 = _mm256_madd_epi16(x6, c3);
            x8 = _mm256_madd_epi16(x8, c4);
            x1 = _mm256_add_epi32(x1, x3);
            x3 = _mm256_add_epi32(x5, x7);
            x2 = _mm256_add_epi32(x2, x4);
            x4 = _mm256_add_epi32(x6, x8);
            x1 = _mm256_add_epi32(x1, x3);
            x2 = _mm256_add_epi32(x2, x4);
            x1 = _mm256_srai_epi32(x1, 6);
            x2 = _mm256_srai_epi32(x2, 6);
            x1 = _mm256_packs_epi32(x1, x2);
            r5 = _mm256_load_si256((__m256i*)&src2[x]);
            x1 = _mm256_adds_epi16(x1, r5);
            x1 = _mm256_mulhrs_epi16(x1, offset);
            x1 = _mm256_max_epi16(x1, _mm256_setzero_si256());
            x1 = _mm256_min_epi16(x1, _mm256_set1_epi16(0x03FF));
            _mm256_storeu_si256((__m256i*)&dst[x], x1);
        }
        src += srcstride;
        src2 += MAX_PB_SIZE;
        dst += dststride;
    }
}

static void
oh_hevc_put_hevc_uni_qpel_v16_14_10_avx2(uint16_t* _dst,
                                         ptrdiff_t _dststride,
                                         const int16_t* _src,
                                         ptrdiff_t _srcstride,
                                         int height,
                                         intptr_t mx,
                                         intptr_t my,
                                         int width)
{
    int x, y;
    const int16_t* src = _src;
    const int srcstride = _srcstride;
    const __m256i offset = _mm256_set1_epi16(1 << (BITDEPTH + 1));
    uint16_t* dst = (uint16_t*)_dst;
    const int dststride = _dststride;

    __m256i c1 = _mm256_set1_epi32(ov_mcp_filters_l[my - 1][0]);
    __m256i c2 = _mm256_set1_epi32(ov_mcp_filters_l[my - 1][1]);
    __m256i c3 = _mm256_set1_epi32(ov_mcp_filters_l[my - 1][2]);
    __m256i c4 = _mm256_set1_epi32(ov_mcp_filters_l[my - 1][3]);

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x += 16) {
            __m256i x1, x2, x3, x4, x5, x6, x7, x8, x9;
            x1 = _mm256_loadu_si256((__m256i*)&src[x - 3 * srcstride]);
            x2 = _mm256_loadu_si256((__m256i*)&src[x - 2 * srcstride]);
            x3 = _mm256_loadu_si256((__m256i*)&src[x - 1 * srcstride]);
            x4 = _mm256_loadu_si256((__m256i*)&src[x]);
            x5 = _mm256_loadu_si256((__m256i*)&src[x + srcstride]);
            x6 = _mm256_loadu_si256((__m256i*)&src[x + 2 * srcstride]);
            x7 = _mm256_loadu_si256((__m256i*)&src[x + 3 * srcstride]);
            x8 = _mm256_loadu_si256((__m256i*)&src[x + 4 * srcstride]);
            x9 = x1;
            x1 = _mm256_unpacklo_epi16(x9, x2);
            x2 = _mm256_unpackhi_epi16(x9, x2);
            x9 = x3;
            x3 = _mm256_unpacklo_epi16(x9, x4);
            x4 = _mm256_unpackhi_epi16(x9, x4);
            x9 = x5;
            x5 = _mm256_unpacklo_epi16(x9, x6);
            x6 = _mm256_unpackhi_epi16(x9, x6);
            x9 = x7;
            x7 = _mm256_unpacklo_epi16(x9, x8);
            x8 = _mm256_unpackhi_epi16(x9, x8);
            x1 = _mm256_madd_epi16(x1, c1);
            x3 = _mm256_madd_epi16(x3, c2);
            x5 = _mm256_madd_epi16(x5, c3);
            x7 = _mm256_madd_epi16(x7, c4);
            x2 = _mm256_madd_epi16(x2, c1);
            x4 = _mm256_madd_epi16(x4, c2);
            x6 = _mm256_madd_epi16(x6, c3);
            x8 = _mm256_madd_epi16(x8, c4);
            x1 = _mm256_add_epi32(x1, x3);
            x3 = _mm256_add_epi32(x5, x7);
            x2 = _mm256_add_epi32(x2, x4);
            x4 = _mm256_add_epi32(x6, x8);
            x1 = _mm256_add_epi32(x1, x3);
            x2 = _mm256_add_epi32(x2, x4);
            x1 = _mm256_srai_epi32(x1, 6);
            x2 = _mm256_srai_epi32(x2, 6);
            x1 = _mm256_packs_epi32(x1, x2);
            x1 = _mm256_mulhrs_epi16(x1, offset);
            x1 = _mm256_max_epi16(x1, _mm256_setzero_si256());
            x1 = _mm256_min_epi16(x1, _mm256_set1_epi16(0x03FF));
            _mm256_storeu_si256((__m256i*)&dst[x], x1);
        }
        src += srcstride;
        dst += dststride;
    }
}

static void
oh_hevc_put_hevc_bi0_qpel_hv16_10_avx2(int16_t* dst,
                                       const uint16_t* _src,
                                       ptrdiff_t _srcstride,
                                       int height,
                                       intptr_t mx,
                                       intptr_t my,
                                       int width)
{
    int16_t tmp_array[(MAX_PB_SIZE + QPEL_EXTRA) * MAX_PB_SIZE];
    int16_t* tmp = tmp_array;
    const uint16_t* src = _src;
    const int srcstride = _srcstride;
    src -= QPEL_EXTRA_BEFORE * srcstride;
    oh_hevc_put_hevc_bi0_qpel_h16_10_avx2(tmp, src, _srcstride, height + QPEL_EXTRA, mx, my, width);
    tmp = tmp_array + QPEL_EXTRA_BEFORE * MAX_PB_SIZE;
    oh_hevc_put_hevc_bi0_qpel_v16_14_avx2(dst, (const uint16_t*)tmp, MAX_PB_SIZE, height, mx, my, width);
}

static void
oh_hevc_put_hevc_bi1_qpel_hv16_10_avx2(uint16_t* dst,
                                       ptrdiff_t dststride,
                                       const uint16_t* _src,
                                       ptrdiff_t _srcstride,
                                       const int16_t* src2,
                                       int height,
                                       intptr_t mx,
                                       intptr_t my,
                                       int width)
{
    int16_t tmp_array[(MAX_PB_SIZE + QPEL_EXTRA) * MAX_PB_SIZE];
    int16_t* tmp = tmp_array;
    const uint16_t* src = _src;
    const int srcstride = _srcstride;
    src -= QPEL_EXTRA_BEFORE * srcstride;
    oh_hevc_put_hevc_bi0_qpel_h16_10_avx2(tmp, src, _srcstride, height + QPEL_EXTRA, mx, my, width);
    tmp = tmp_array + QPEL_EXTRA_BEFORE * MAX_PB_SIZE;
    oh_hevc_put_hevc_bi1_qpel_v16_14_10_avx2(dst, dststride, tmp, MAX_PB_SIZE, src2, height, mx, my, width);
}

static void
oh_hevc_put_hevc_uni_qpel_hv16_10_avx2(uint16_t* dst,
                                       ptrdiff_t dststride,
                                       const uint16_t* _src,
                                       ptrdiff_t _srcstride,
                                       int height,
                                       intptr_t mx,
                                       intptr_t my,
                                       int width)
{
    int16_t tmp_array[(MAX_PB_SIZE + QPEL_EXTRA) * MAX_PB_SIZE];
    int16_t* tmp = tmp_array;
    const uint16_t* src = _src;
    const int srcstride = _srcstride;
    src -= QPEL_EXTRA_BEFORE * srcstride;
    oh_hevc_put_hevc_bi0_qpel_h16_10_avx2(tmp, src, _srcstride, height + QPEL_EXTRA, mx, my, width);
    tmp = tmp_array + QPEL_EXTRA_BEFORE * MAX_PB_SIZE;
    oh_hevc_put_hevc_uni_qpel_v16_14_10_avx2(dst, dststride, tmp, MAX_PB_SIZE, height, mx, my, width);
}

static void
oh_hevc_put_hevc_bi0_pel_pixels32_10_avx2(int16_t* dst,
                                          const uint16_t* _src,
                                          ptrdiff_t _srcstride,
                                          int height,
                                          intptr_t mx,
                                          intptr_t my,
                                          int width)
{
    oh_hevc_put_hevc_bi0_pel_pixels16_10_avx2(dst, _src, _srcstride, height, mx, my, width);
}

static void
oh_hevc_put_hevc_bi1_pel_pixels32_10_avx2(uint16_t* dst,
                                          ptrdiff_t dststride,
                                          const uint16_t* _src,
                                          ptrdiff_t _srcstride,
                                          const int16_t* src2,
                                          int height,
                                          intptr_t mx,
                                          intptr_t my,
                                          int width)
{
    oh_hevc_put_hevc_bi1_pel_pixels16_10_avx2(dst, dststride, _src, _srcstride, src2, height, mx, my, width);
}

static void
oh_hevc_put_hevc_bi0_pel_pixels64_10_avx2(int16_t* dst,
                                          const uint16_t* _src,
                                          ptrdiff_t _srcstride,
                                          int height,
                                          intptr_t mx,
                                          intptr_t my,
                                          int width)
{
    oh_hevc_put_hevc_bi0_pel_pixels16_10_avx2(dst, _src, _srcstride, height, mx, my, width);
}

static void
oh_hevc_put_hevc_bi1_pel_pixels64_10_avx2(uint16_t* dst,
                                          ptrdiff_t dststride,
                                          const uint16_t* _src,
                                          ptrdiff_t _srcstride,
                                          const int16_t* src2,
                                          int height,
                                          intptr_t mx,
                                          intptr_t my,
                                          int width)
{
    oh_hevc_put_hevc_bi1_pel_pixels16_10_avx2(dst, dststride, _src, _srcstride, src2, height, mx, my, width);
}

static void
oh_hevc_put_hevc_bi0_qpel_h32_10_avx2(int16_t* dst,
                                      const uint16_t* _src,
                                      ptrdiff_t _srcstride,
                                      int height,
                                      intptr_t mx,
                                      intptr_t my,
                                      int width)
{
    oh_hevc_put_hevc_bi0_qpel_h16_10_avx2(dst, _src, _srcstride, height, mx, my, width);
}

static void
oh_hevc_put_hevc_bi1_qpel_h32_10_avx2(uint16_t* dst,
                                      ptrdiff_t dststride,
                                      const uint16_t* _src,
                                      ptrdiff_t _srcstride,
                                      const int16_t* src2,
                                      int height,
                                      intptr_t mx,
                                      intptr_t my,
                                      int width)
{
    oh_hevc_put_hevc_bi1_qpel_h16_10_avx2(dst, dststride, _src, _srcstride, src2, height, mx, my, width);
}

static void
oh_hevc_put_hevc_uni_qpel_h32_10_avx2(uint16_t* dst,
                                      ptrdiff_t dststride,
                                      const uint16_t* _src,
                                      ptrdiff_t _srcstride,
                                      int height,
                                      intptr_t mx,
                                      intptr_t my,
                                      int width)
{
    oh_hevc_put_hevc_uni_qpel_h16_10_avx2(dst, dststride, _src, _srcstride, height, mx, my, width);
}

static void
oh_hevc_put_hevc_bi0_qpel_h64_10_avx2(int16_t* dst,
                                      const uint16_t* _src,
                                      ptrdiff_t _srcstride,
                                      int height,
                                      intptr_t mx,
                                      intptr_t my,
                                      int width)
{
    oh_hevc_put_hevc_bi0_qpel_h16_10_avx2(dst, _src, _srcstride, height, mx, my, width);
}

static void
oh_hevc_put_hevc_bi1_qpel_h64_10_avx2(uint16_t* dst,
                                      ptrdiff_t dststride,
                                      const uint16_t* _src,
                                      ptrdiff_t _srcstride,
                                      const int16_t* src2,
                                      int height,
                                      intptr_t mx,
                                      intptr_t my,
                                      int width)
{
    oh_hevc_put_hevc_bi1_qpel_h16_10_avx2(dst, dststride, _src, _srcstride, src2, height, mx, my, width);
}

static void
oh_hevc_put_hevc_uni_qpel_h64_10_avx2(uint16_t* dst,
                                      ptrdiff_t dststride,
                                      const uint16_t* _src,
                                      ptrdiff_t _srcstride,
                                      int height,
                                      intptr_t mx,
                                      intptr_t my,
                                      int width)
{
    oh_hevc_put_hevc_uni_qpel_h16_10_avx2(dst, dststride, _src, _srcstride, height, mx, my, width);
}

static void
oh_hevc_put_hevc_bi0_qpel_v32_10_avx2(int16_t* dst,
                                      const uint16_t* _src,
                                      ptrdiff_t _srcstride,
                                      int height,
                                      intptr_t mx,
                                      intptr_t my,
                                      int width)
{
    oh_hevc_put_hevc_bi0_qpel_v16_10_avx2(dst, _src, _srcstride, height, mx, my, width);
}

static void
oh_hevc_put_hevc_bi1_qpel_v32_10_avx2(uint16_t* dst,
                                      ptrdiff_t dststride,
                                      const uint16_t* _src,
                                      ptrdiff_t _srcstride,
                                      const int16_t* src2,
                                      int height,
                                      intptr_t mx,
                                      intptr_t my,
                                      int width)
{
    oh_hevc_put_hevc_bi1_qpel_v16_10_avx2(dst, dststride, _src, _srcstride, src2, height, mx, my, width);
}

static void
oh_hevc_put_hevc_uni_qpel_v32_10_avx2(uint16_t* dst,
                                      ptrdiff_t dststride,
                                      const uint16_t* _src,
                                      ptrdiff_t _srcstride,
                                      int height,
                                      intptr_t mx,
                                      intptr_t my,
                                      int width)
{
    oh_hevc_put_hevc_uni_qpel_v16_10_avx2(dst, dststride, _src, _srcstride, height, mx, my, width);
}

static void
oh_hevc_put_hevc_bi0_qpel_v64_10_avx2(int16_t* dst,
                                      const uint16_t* _src,
                                      ptrdiff_t _srcstride,
                                      int height,
                                      intptr_t mx,
                                      intptr_t my,
                                      int width)
{
    oh_hevc_put_hevc_bi0_qpel_v16_10_avx2(dst, _src, _srcstride, height, mx, my, width);
}

static void
oh_hevc_put_hevc_bi1_qpel_v64_10_avx2(uint16_t* dst,
                                      ptrdiff_t dststride,
                                      const uint16_t* _src,
                                      ptrdiff_t _srcstride,
                                      const int16_t* src2,
                                      int height,
                                      intptr_t mx,
                                      intptr_t my,
                                      int width)
{
    oh_hevc_put_hevc_bi1_qpel_v16_10_avx2(dst, dststride, _src, _srcstride, src2, height, mx, my, width);
}

static void
oh_hevc_put_hevc_uni_qpel_v64_10_avx2(uint16_t* dst,
                                      ptrdiff_t dststride,
                                      const uint16_t* _src,
                                      ptrdiff_t _srcstride,
                                      int height,
                                      intptr_t mx,
                                      intptr_t my,
                                      int width)
{
    oh_hevc_put_hevc_uni_qpel_v16_10_avx2(dst, dststride, _src, _srcstride, height, mx, my, width);
}

static void
oh_hevc_put_hevc_bi0_qpel_hv32_10_avx2(int16_t* dst,
                                       const uint16_t* _src,
                                       ptrdiff_t _srcstride,
                                       int height,
                                       intptr_t mx,
                                       intptr_t my,
                                       int width)
{
    oh_hevc_put_hevc_bi0_qpel_hv16_10_avx2(dst, _src, _srcstride, height, mx, my, width);
}

static void
oh_hevc_put_hevc_bi1_qpel_hv32_10_avx2(uint16_t* dst,
                                       ptrdiff_t dststride,
                                       const uint16_t* _src,
                                       ptrdiff_t _srcstride,
                                       const int16_t* src2,
                                       int height,
                                       intptr_t mx,
                                       intptr_t my,
                                       int width)
{
    oh_hevc_put_hevc_bi1_qpel_hv16_10_avx2(dst, dststride, _src, _srcstride, src2, height, mx, my, width);
}

static void
oh_hevc_put_hevc_uni_qpel_hv32_10_avx2(uint16_t* dst,
                                       ptrdiff_t dststride,
                                       const uint16_t* _src,
                                       ptrdiff_t _srcstride,
                                       int height,
                                       intptr_t mx,
                                       intptr_t my,
                                       int width)
{
    oh_hevc_put_hevc_uni_qpel_hv16_10_avx2(dst, dststride, _src, _srcstride, height, mx, my, width);
}

static void
oh_hevc_put_hevc_bi0_qpel_hv64_10_avx2(int16_t* dst,
                                       const uint16_t* _src,
                                       ptrdiff_t _srcstride,
                                       int height,
                                       intptr_t mx,
                                       intptr_t my,
                                       int width)
{
    oh_hevc_put_hevc_bi0_qpel_hv16_10_avx2(dst, _src, _srcstride, height, mx, my, width);
}

static void
oh_hevc_put_hevc_bi1_qpel_hv64_10_avx2(uint16_t* dst,
                                       ptrdiff_t dststride,
                                       const uint16_t* _src,
                                       ptrdiff_t _srcstride,
                                       const int16_t* src2,
                                       int height,
                                       intptr_t mx,
                                       intptr_t my,
                                       int width)
{
    oh_hevc_put_hevc_bi1_qpel_hv16_10_avx2(dst, dststride, _src, _srcstride, src2, height, mx, my, width);
}

static void
oh_hevc_put_hevc_uni_qpel_hv64_10_avx2(uint16_t* dst,
                                       ptrdiff_t dststride,
                                       const uint16_t* _src,
                                       ptrdiff_t _srcstride,
                                       int height,
                                       intptr_t mx,
                                       intptr_t my,
                                       int width)
{
    oh_hevc_put_hevc_uni_qpel_hv16_10_avx2(dst, dststride, _src, _srcstride, height, mx, my, width);
}

static void
oh_hevc_put_hevc_bi0_epel_h32_10_avx2(int16_t* dst,
                                      const uint16_t* _src,
                                      ptrdiff_t _srcstride,
                                      int height,
                                      intptr_t mx,
                                      intptr_t my,
                                      int width)
{
    oh_hevc_put_hevc_bi0_epel_h16_10_avx2(dst, _src, _srcstride, height, mx, my, width);
}

static void
oh_hevc_put_hevc_bi1_epel_h32_10_avx2(uint16_t* dst,
                                      ptrdiff_t dststride,
                                      const uint16_t* _src,
                                      ptrdiff_t _srcstride,
                                      const int16_t* src2,
                                      int height,
                                      intptr_t mx,
                                      intptr_t my,
                                      int width)
{
    oh_hevc_put_hevc_bi1_epel_h16_10_avx2(dst, dststride, _src, _srcstride, src2, height, mx, my, width);
}

static void
oh_hevc_put_hevc_uni_epel_h32_10_avx2(uint16_t* dst,
                                      ptrdiff_t dststride,
                                      const uint16_t* _src,
                                      ptrdiff_t _srcstride,
                                      int height,
                                      intptr_t mx,
                                      intptr_t my,
                                      int width)
{
    oh_hevc_put_hevc_uni_epel_h16_10_avx2(dst, dststride, _src, _srcstride, height, mx, my, width);
}

static void
oh_hevc_put_hevc_bi0_epel_h64_10_avx2(int16_t* dst,
                                      const uint16_t* _src,
                                      ptrdiff_t _srcstride,
                                      int height,
                                      intptr_t mx,
                                      intptr_t my,
                                      int width)
{
    oh_hevc_put_hevc_bi0_epel_h16_10_avx2(dst, _src, _srcstride, height, mx, my, width);
}

static void
oh_hevc_put_hevc_bi1_epel_h64_10_avx2(uint16_t* dst,
                                      ptrdiff_t dststride,
                                      const uint16_t* _src,
                                      ptrdiff_t _srcstride,
                                      const int16_t* src2,
                                      int height,
                                      intptr_t mx,
                                      intptr_t my,
                                      int width)
{
    oh_hevc_put_hevc_bi1_epel_h16_10_avx2(dst, dststride, _src, _srcstride, src2, height, mx, my, width);
}

static void
oh_hevc_put_hevc_uni_epel_h64_10_avx2(uint16_t* dst,
                                      ptrdiff_t dststride,
                                      const uint16_t* _src,
                                      ptrdiff_t _srcstride,
                                      int height,
                                      intptr_t mx,
                                      intptr_t my,
                                      int width)
{
    oh_hevc_put_hevc_uni_epel_h16_10_avx2(dst, dststride, _src, _srcstride, height, mx, my, width);
}

static void
oh_hevc_put_hevc_bi0_epel_v32_10_avx2(int16_t* dst,
                                      const uint16_t* _src,
                                      ptrdiff_t _srcstride,
                                      int height,
                                      intptr_t mx,
                                      intptr_t my,
                                      int width)
{
    oh_hevc_put_hevc_bi0_epel_v16_10_avx2(dst, _src, _srcstride, height, mx, my, width);
}

static void
oh_hevc_put_hevc_bi1_epel_v32_10_avx2(uint16_t* dst,
                                      ptrdiff_t dststride,
                                      const uint16_t* _src,
                                      ptrdiff_t _srcstride,
                                      const int16_t* src2,
                                      int height,
                                      intptr_t mx,
                                      intptr_t my,
                                      int width)
{
    oh_hevc_put_hevc_bi1_epel_v16_10_avx2(dst, dststride, _src, _srcstride, src2, height, mx, my, width);
}

static void
oh_hevc_put_hevc_uni_epel_v32_10_avx2(uint16_t* dst,
                                      ptrdiff_t dststride,
                                      const uint16_t* _src,
                                      ptrdiff_t _srcstride,
                                      int height,
                                      intptr_t mx,
                                      intptr_t my,
                                      int width)
{
    oh_hevc_put_hevc_uni_epel_v16_10_avx2(dst, dststride, _src, _srcstride, height, mx, my, width);
}

static void
oh_hevc_put_hevc_bi0_epel_v64_10_avx2(int16_t* dst,
                                      const uint16_t* _src,
                                      ptrdiff_t _srcstride,
                                      int height,
                                      intptr_t mx,
                                      intptr_t my,
                                      int width)
{
    oh_hevc_put_hevc_bi0_epel_v16_10_avx2(dst, _src, _srcstride, height, mx, my, width);
}

static void
oh_hevc_put_hevc_bi1_epel_v64_10_avx2(uint16_t* dst,
                                      ptrdiff_t dststride,
                                      const uint16_t* _src,
                                      ptrdiff_t _srcstride,
                                      const int16_t* src2,
                                      int height,
                                      intptr_t mx,
                                      intptr_t my,
                                      int width)
{
    oh_hevc_put_hevc_bi1_epel_v16_10_avx2(dst, dststride, _src, _srcstride, src2, height, mx, my, width);
}

static void
oh_hevc_put_hevc_uni_epel_v64_10_avx2(uint16_t* dst,
                                      ptrdiff_t dststride,
                                      const uint16_t* _src,
                                      ptrdiff_t _srcstride,
                                      int height,
                                      intptr_t mx,
                                      intptr_t my,
                                      int width)
{
    oh_hevc_put_hevc_uni_epel_v16_10_avx2(dst, dststride, _src, _srcstride, height, mx, my, width);
}

static void
oh_hevc_put_hevc_bi0_epel_hv32_10_avx2(int16_t* dst,
                                       const uint16_t* _src,
                                       ptrdiff_t _srcstride,
                                       int height,
                                       intptr_t mx,
                                       intptr_t my,
                                       int width)
{
    oh_hevc_put_hevc_bi0_epel_hv16_10_avx2(dst, _src, _srcstride, height, mx, my, width);
}

static void
oh_hevc_put_hevc_bi1_epel_hv32_10_avx2(uint16_t* dst,
                                       ptrdiff_t dststride,
                                       const uint16_t* _src,
                                       ptrdiff_t _srcstride,
                                       const int16_t* src2,
                                       int height,
                                       intptr_t mx,
                                       intptr_t my,
                                       int width)
{
    oh_hevc_put_hevc_bi1_epel_hv16_10_avx2(dst, dststride, _src, _srcstride, src2, height, mx, my, width);
}

static void
oh_hevc_put_hevc_uni_epel_hv32_10_avx2(uint16_t* dst,
                                       ptrdiff_t dststride,
                                       const uint16_t* _src,
                                       ptrdiff_t _srcstride,
                                       int height,
                                       intptr_t mx,
                                       intptr_t my,
                                       int width)
{
    oh_hevc_put_hevc_uni_epel_hv16_10_avx2(dst, dststride, _src, _srcstride, height, mx, my, width);
}

static void
oh_hevc_put_hevc_bi0_epel_hv64_10_avx2(int16_t* dst,
                                       const uint16_t* _src,
                                       ptrdiff_t _srcstride,
                                       int height,
                                       intptr_t mx,
                                       intptr_t my,
                                       int width)
{
    oh_hevc_put_hevc_bi0_epel_hv16_10_avx2(dst, _src, _srcstride, height, mx, my, width);
}

static void
oh_hevc_put_hevc_bi1_epel_hv64_10_avx2(uint16_t* dst,
                                       ptrdiff_t dststride,
                                       const uint16_t* _src,
                                       ptrdiff_t _srcstride,
                                       const int16_t* src2,
                                       int height,
                                       intptr_t mx,
                                       intptr_t my,
                                       int width)
{
    oh_hevc_put_hevc_bi1_epel_hv16_10_avx2(dst, dststride, _src, _srcstride, src2, height, mx, my, width);
}

static void
oh_hevc_put_hevc_uni_epel_hv64_10_avx2(uint16_t* dst,
                                       ptrdiff_t dststride,
                                       const uint16_t* _src,
                                       ptrdiff_t _srcstride,
                                       int height,
                                       intptr_t mx,
                                       intptr_t my,
                                       int width)
{
    oh_hevc_put_hevc_uni_epel_hv16_10_avx2(dst, dststride, _src, _srcstride, height, mx, my, width);
}

static void
put_vvc_qpel_h16_10_avx2(int16_t* dst,
                         ptrdiff_t dststride,
                         OVSample* _src,
                         ptrdiff_t _srcstride,
                         int height,
                         intptr_t mx,
                         intptr_t my,
                         int width)
{
    int x, y;
    int shift = BITDEPTH - 8;
    uint16_t* src = (uint16_t*)_src;
    const int srcstride = _srcstride;

    __m256i c1 = _mm256_set1_epi32(ov_mcp_filters_l[mx - 1][0]);
    __m256i c2 = _mm256_set1_epi32(ov_mcp_filters_l[mx - 1][1]);
    __m256i c3 = _mm256_set1_epi32(ov_mcp_filters_l[mx - 1][2]);
    __m256i c4 = _mm256_set1_epi32(ov_mcp_filters_l[mx - 1][3]);

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x += 16) {
            __m256i x1, x2, x3, x4, r1, r2, r3, r4, t1, t2;
            x1 = _mm256_loadu_si256((__m256i*)&src[x - 3 * 1]);
            x2 = _mm256_loadu_si256((__m256i*)&src[x - 2 * 1]);
            x3 = _mm256_loadu_si256((__m256i*)&src[x - 1 * 1]);
            x4 = _mm256_loadu_si256((__m256i*)&src[x]);
            t1 = _mm256_unpackhi_epi16(x1, x2);
            x1 = _mm256_unpacklo_epi16(x1, x2);
            t2 = _mm256_unpackhi_epi16(x3, x4);
            x2 = _mm256_unpacklo_epi16(x3, x4);
            t1 = _mm256_madd_epi16(t1, c1);
            t2 = _mm256_madd_epi16(t2, c2);
            r2 = _mm256_add_epi32(t1, t2);
            x1 = _mm256_madd_epi16(x1, c1);
            x2 = _mm256_madd_epi16(x2, c2);
            r1 = _mm256_add_epi32(x1, x2);
            x1 = _mm256_loadu_si256((__m256i*)&src[x + 1]);
            x2 = _mm256_loadu_si256((__m256i*)&src[x + 2 * 1]);
            x3 = _mm256_loadu_si256((__m256i*)&src[x + 3 * 1]);
            x4 = _mm256_loadu_si256((__m256i*)&src[x + 4 * 1]);
            t1 = _mm256_unpackhi_epi16(x1, x2);
            x1 = _mm256_unpacklo_epi16(x1, x2);
            t2 = _mm256_unpackhi_epi16(x3, x4);
            x2 = _mm256_unpacklo_epi16(x3, x4);
            t1 = _mm256_madd_epi16(t1, c3);
            t2 = _mm256_madd_epi16(t2, c4);
            r4 = _mm256_add_epi32(t1, t2);
            x1 = _mm256_madd_epi16(x1, c3);
            x2 = _mm256_madd_epi16(x2, c4);
            r3 = _mm256_add_epi32(x1, x2);
            x1 = _mm256_add_epi32(r1, r3);
            x2 = _mm256_add_epi32(r2, r4);
            x1 = _mm256_srai_epi32(x1, shift);
            x2 = _mm256_srai_epi32(x2, shift);
            x1 = _mm256_packs_epi32(x1, x2);
            _mm256_storeu_si256((__m256i*)&dst[x], x1);
        }
        src += srcstride;
        dst += dststride;
    }
}

static void
put_vvc_bi_w_qpel_v16_14_10_avx2(OVSample* _dst,
                                 ptrdiff_t _dststride,
                                 const OVSample* _src,
                                 ptrdiff_t _srcstride,
                                 const int16_t* src2,
                                 ptrdiff_t src2stride,
                                 int height,
                                 int denom,
                                 int _wx0,
                                 int _wx1,
                                 int offset0,
                                 int offset1,

                                 intptr_t mx,
                                 intptr_t my,
                                 int width)
{
    int x, y;
    uint16_t* src = (uint16_t*)_src;
    const int srcstride = _srcstride;
    const int log2Wd = denom + 14 - BITDEPTH;
    const int shift2 = log2Wd + 1;

    const __m256i wx0 = _mm256_set1_epi16(_wx0);
    const __m256i wx1 = _mm256_set1_epi16(_wx1);
    const __m256i offset = _mm256_set1_epi32((offset0 + offset1 + 1) << log2Wd);
    uint16_t* dst = (uint16_t*)_dst;
    const int dststride = _dststride;

    __m256i c1 = _mm256_set1_epi32(ov_mcp_filters_l[my - 1][0]);
    __m256i c2 = _mm256_set1_epi32(ov_mcp_filters_l[my - 1][1]);
    __m256i c3 = _mm256_set1_epi32(ov_mcp_filters_l[my - 1][2]);
    __m256i c4 = _mm256_set1_epi32(ov_mcp_filters_l[my - 1][3]);

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x += 16) {
            __m256i x1, x2, x3, x4, x5, x6, x7, x8, x9;
            __m256i r5;
            x1 = _mm256_loadu_si256((__m256i*)&src[x - 3 * srcstride]);
            x2 = _mm256_loadu_si256((__m256i*)&src[x - 2 * srcstride]);
            x3 = _mm256_loadu_si256((__m256i*)&src[x - 1 * srcstride]);
            x4 = _mm256_loadu_si256((__m256i*)&src[x]);
            x5 = _mm256_loadu_si256((__m256i*)&src[x + srcstride]);
            x6 = _mm256_loadu_si256((__m256i*)&src[x + 2 * srcstride]);
            x7 = _mm256_loadu_si256((__m256i*)&src[x + 3 * srcstride]);
            x8 = _mm256_loadu_si256((__m256i*)&src[x + 4 * srcstride]);
            x9 = x1;
            x1 = _mm256_unpacklo_epi16(x9, x2);
            x2 = _mm256_unpackhi_epi16(x9, x2);
            x9 = x3;
            x3 = _mm256_unpacklo_epi16(x9, x4);
            x4 = _mm256_unpackhi_epi16(x9, x4);
            x9 = x5;
            x5 = _mm256_unpacklo_epi16(x9, x6);
            x6 = _mm256_unpackhi_epi16(x9, x6);
            x9 = x7;
            x7 = _mm256_unpacklo_epi16(x9, x8);
            x8 = _mm256_unpackhi_epi16(x9, x8);
            x1 = _mm256_madd_epi16(x1, c1);
            x3 = _mm256_madd_epi16(x3, c2);
            x5 = _mm256_madd_epi16(x5, c3);
            x7 = _mm256_madd_epi16(x7, c4);
            x2 = _mm256_madd_epi16(x2, c1);
            x4 = _mm256_madd_epi16(x4, c2);
            x6 = _mm256_madd_epi16(x6, c3);
            x8 = _mm256_madd_epi16(x8, c4);
            x1 = _mm256_add_epi32(x1, x3);
            x3 = _mm256_add_epi32(x5, x7);
            x2 = _mm256_add_epi32(x2, x4);
            x4 = _mm256_add_epi32(x6, x8);
            x1 = _mm256_add_epi32(x1, x3);
            x2 = _mm256_add_epi32(x2, x4);
            x1 = _mm256_srai_epi32(x1, 6);
            x2 = _mm256_srai_epi32(x2, 6);
            x1 = _mm256_packs_epi32(x1, x2);
            r5 = _mm256_load_si256((__m256i*)&src2[x]);
            {
                __m256i x3, x4, r7, r8;
                x3 = _mm256_mulhi_epi16(x1, wx1);
                x1 = _mm256_mullo_epi16(x1, wx1);
                r7 = _mm256_mulhi_epi16(r5, wx0);
                r5 = _mm256_mullo_epi16(r5, wx0);
                x4 = _mm256_unpackhi_epi16(x1, x3);
                x1 = _mm256_unpacklo_epi16(x1, x3);
                r8 = _mm256_unpackhi_epi16(r5, r7);
                r5 = _mm256_unpacklo_epi16(r5, r7);
                x4 = _mm256_add_epi32(x4, r8);
                x1 = _mm256_add_epi32(x1, r5);
                x4 = _mm256_srai_epi32(_mm256_add_epi32(x4, offset), shift2);
                x1 = _mm256_srai_epi32(_mm256_add_epi32(x1, offset), shift2);
                x1 = _mm256_packus_epi32(x1, x4);
            };
            x1 = _mm256_max_epi16(x1, _mm256_setzero_si256());
            x1 = _mm256_min_epi16(x1, _mm256_set1_epi16(0x03FF));
            _mm256_storeu_si256((__m256i*)&dst[x], x1);
        }
        src += srcstride;
        src2 += src2stride;
        dst += dststride;
    }
}

static void
put_vvc_uni_w_qpel_v16_14_10_avx2(OVSample* _dst,
                                  ptrdiff_t _dststride,
                                  const OVSample* _src,
                                  ptrdiff_t _srcstride,
                                  int height,
                                  int denom,
                                  int _wx,
                                  int _ox,
                                  intptr_t mx,
                                  intptr_t my,
                                  int width)
{
    int x, y;
    uint16_t* src = (uint16_t*)_src;
    const int srcstride = _srcstride;
    const int shift2 = denom + 14 - BITDEPTH;
    const __m256i ox = _mm256_set1_epi32(_ox);
    const __m256i wx = _mm256_set1_epi16(_wx);
    const __m256i offset = _mm256_set1_epi32(1 << (shift2 - 1));
    uint16_t* dst = (uint16_t*)_dst;
    const int dststride = _dststride;

    __m256i c1 = _mm256_set1_epi32(ov_mcp_filters_l[my - 1][0]);
    __m256i c2 = _mm256_set1_epi32(ov_mcp_filters_l[my - 1][1]);
    __m256i c3 = _mm256_set1_epi32(ov_mcp_filters_l[my - 1][2]);
    __m256i c4 = _mm256_set1_epi32(ov_mcp_filters_l[my - 1][3]);

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x += 16) {
            __m256i x1, x2, x3, x4, x5, x6, x7, x8, x9;
            x1 = _mm256_loadu_si256((__m256i*)&src[x - 3 * srcstride]);
            x2 = _mm256_loadu_si256((__m256i*)&src[x - 2 * srcstride]);
            x3 = _mm256_loadu_si256((__m256i*)&src[x - 1 * srcstride]);
            x4 = _mm256_loadu_si256((__m256i*)&src[x]);
            x5 = _mm256_loadu_si256((__m256i*)&src[x + srcstride]);
            x6 = _mm256_loadu_si256((__m256i*)&src[x + 2 * srcstride]);
            x7 = _mm256_loadu_si256((__m256i*)&src[x + 3 * srcstride]);
            x8 = _mm256_loadu_si256((__m256i*)&src[x + 4 * srcstride]);
            x9 = x1;
            x1 = _mm256_unpacklo_epi16(x9, x2);
            x2 = _mm256_unpackhi_epi16(x9, x2);
            x9 = x3;
            x3 = _mm256_unpacklo_epi16(x9, x4);
            x4 = _mm256_unpackhi_epi16(x9, x4);
            x9 = x5;
            x5 = _mm256_unpacklo_epi16(x9, x6);
            x6 = _mm256_unpackhi_epi16(x9, x6);
            x9 = x7;
            x7 = _mm256_unpacklo_epi16(x9, x8);
            x8 = _mm256_unpackhi_epi16(x9, x8);
            x1 = _mm256_madd_epi16(x1, c1);
            x3 = _mm256_madd_epi16(x3, c2);
            x5 = _mm256_madd_epi16(x5, c3);
            x7 = _mm256_madd_epi16(x7, c4);
            x2 = _mm256_madd_epi16(x2, c1);
            x4 = _mm256_madd_epi16(x4, c2);
            x6 = _mm256_madd_epi16(x6, c3);
            x8 = _mm256_madd_epi16(x8, c4);
            x1 = _mm256_add_epi32(x1, x3);
            x3 = _mm256_add_epi32(x5, x7);
            x2 = _mm256_add_epi32(x2, x4);
            x4 = _mm256_add_epi32(x6, x8);
            x1 = _mm256_add_epi32(x1, x3);
            x2 = _mm256_add_epi32(x2, x4);
            x1 = _mm256_srai_epi32(x1, 6);
            x2 = _mm256_srai_epi32(x2, 6);
            x1 = _mm256_packs_epi32(x1, x2);
            {
                __m256i x3, x4;
                x3 = _mm256_mulhi_epi16(x1, wx);
                x1 = _mm256_mullo_epi16(x1, wx);
                x4 = _mm256_unpackhi_epi16(x1, x3);
                x1 = _mm256_unpacklo_epi16(x1, x3);
                x3 = _mm256_srai_epi32(_mm256_add_epi32(x4, offset), shift2);
                x1 = _mm256_srai_epi32(_mm256_add_epi32(x1, offset), shift2);
                x3 = _mm256_add_epi32(x3, ox);
                x1 = _mm256_add_epi32(x1, ox);
                x1 = _mm256_packus_epi32(x1, x3);
            };
            x1 = _mm256_max_epi16(x1, _mm256_setzero_si256());
            x1 = _mm256_min_epi16(x1, _mm256_set1_epi16(0x03FF));
            _mm256_storeu_si256((__m256i*)&dst[x], x1);
        }
        src += srcstride;
        dst += dststride;
    }
}

static void
put_vvc_uni_w_pel_pixels8_10_avx2(OVSample* _dst,
                                  ptrdiff_t _dststride,
                                  const OVSample* _src,
                                  ptrdiff_t _srcstride,
                                  int height,
                                  int denom,
                                  int _wx,
                                  int _ox,
                                  intptr_t mx,
                                  intptr_t my,
                                  int width)
{
    int x, y;
    uint16_t* src = (uint16_t*)_src;
    const int srcstride = _srcstride;
    const int shift2 = denom + 14 - BITDEPTH;
    const __m256i ox = _mm256_set1_epi32(_ox);
    const __m256i wx = _mm256_set1_epi16(_wx);
    const __m256i offset = _mm256_set1_epi32(1 << (shift2 - 1));
    uint16_t* dst = (uint16_t*)_dst;
    const int dststride = _dststride;
    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x += 16) {
            __m256i x1;
            x1 = _mm256_loadu_si256((__m256i*)&src[x]);
            x1 = _mm256_slli_epi16(x1, 14 - BITDEPTH);
            {
                __m256i x3, x4;
                x3 = _mm256_mulhi_epi16(x1, wx);
                x1 = _mm256_mullo_epi16(x1, wx);
                x4 = _mm256_unpackhi_epi16(x1, x3);
                x1 = _mm256_unpacklo_epi16(x1, x3);
                x3 = _mm256_srai_epi32(_mm256_add_epi32(x4, offset), shift2);
                x1 = _mm256_srai_epi32(_mm256_add_epi32(x1, offset), shift2);
                x3 = _mm256_add_epi32(x3, ox);
                x1 = _mm256_add_epi32(x1, ox);
                x1 = _mm256_packus_epi32(x1, x3);
            };
            x1 = _mm256_max_epi16(x1, _mm256_setzero_si256());
            x1 = _mm256_min_epi16(x1, _mm256_set1_epi16(0x03FF));
            _mm256_storeu_si256((__m256i*)&dst[x], x1);
        }
        src += srcstride;
        dst += dststride;
    }
}

static void
put_vvc_bi_w_pel_pixels8_10_avx2(OVSample* _dst,
                                 ptrdiff_t _dststride,
                                 const OVSample* _src,
                                 ptrdiff_t _srcstride,
                                 const int16_t* src2,
                                 ptrdiff_t src2stride,
                                 int height,
                                 int denom,
                                 int _wx0,
                                 int _wx1,
                                 int offset0,
                                 int offset1,
                                 intptr_t mx,
                                 intptr_t my,
                                 int width)
{
    int x, y;
    uint16_t* src = (uint16_t*)_src;
    const int srcstride = _srcstride;
    const int log2Wd = denom + 14 - BITDEPTH;
    const int shift2 = log2Wd + 1;
    const __m256i wx0 = _mm256_set1_epi16(_wx0);
    const __m256i wx1 = _mm256_set1_epi16(_wx1);
    const __m256i offset = _mm256_set1_epi32((offset0 + offset1 + 1) << log2Wd);
    uint16_t* dst = (uint16_t*)_dst;
    const int dststride = _dststride;
    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x += 16) {
            __m256i x1;
            __m256i r5;
            x1 = _mm256_loadu_si256((__m256i*)&src[x]);
            x1 = _mm256_slli_epi16(x1, 14 - BITDEPTH);
            r5 = _mm256_load_si256((__m256i*)&src2[x]);
            {
                __m256i x3, x4, r7, r8;
                x3 = _mm256_mulhi_epi16(x1, wx1);
                x1 = _mm256_mullo_epi16(x1, wx1);
                r7 = _mm256_mulhi_epi16(r5, wx0);
                r5 = _mm256_mullo_epi16(r5, wx0);
                x4 = _mm256_unpackhi_epi16(x1, x3);
                x1 = _mm256_unpacklo_epi16(x1, x3);
                r8 = _mm256_unpackhi_epi16(r5, r7);
                r5 = _mm256_unpacklo_epi16(r5, r7);
                x4 = _mm256_add_epi32(x4, r8);
                x1 = _mm256_add_epi32(x1, r5);
                x4 = _mm256_srai_epi32(_mm256_add_epi32(x4, offset), shift2);
                x1 = _mm256_srai_epi32(_mm256_add_epi32(x1, offset), shift2);
                x1 = _mm256_packus_epi32(x1, x4);
            };
            x1 = _mm256_max_epi16(x1, _mm256_setzero_si256());
            x1 = _mm256_min_epi16(x1, _mm256_set1_epi16(0x03FF));
            _mm256_storeu_si256((__m256i*)&dst[x], x1);
        }
        src += srcstride;
        src2 += src2stride;
        dst += dststride;
    }
}

static void
put_vvc_uni_w_epel_h16_10_avx2(OVSample* _dst,
                               ptrdiff_t _dststride,
                               const OVSample* _src,
                               ptrdiff_t _srcstride,
                               int height,
                               int denom,
                               int _wx,
                               int _ox,
                               intptr_t mx,
                               intptr_t my,
                               int width)
{
    int x, y;
    uint16_t* src = ((uint16_t*)_src) - 1;
    ptrdiff_t srcstride = _srcstride;
    const int shift2 = denom + 14 - BITDEPTH;
    const __m256i ox = _mm256_set1_epi32(_ox);
    const __m256i wx = _mm256_set1_epi16(_wx);
    const __m256i offset = _mm256_set1_epi32(1 << (shift2 - 1));
    uint16_t* dst = (uint16_t*)_dst;
    const int dststride = _dststride;
    __m256i f1 = _mm256_set1_epi32(ov_mcp_filters_c[mx - 1][0]);
    __m256i f2 = _mm256_set1_epi32(ov_mcp_filters_c[mx - 1][1]);
    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x += 16) {
            __m256i x1, x2, x3, x4, t1, t2;
            x1 = _mm256_loadu_si256((__m256i*)&src[x]);
            x2 = _mm256_loadu_si256((__m256i*)&src[x + 1]);
            x3 = _mm256_loadu_si256((__m256i*)&src[x + 2 * 1]);
            x4 = _mm256_loadu_si256((__m256i*)&src[x + 3 * 1]);
            t1 = _mm256_unpackhi_epi16(x1, x2);
            x1 = _mm256_unpacklo_epi16(x1, x2);
            t2 = _mm256_unpackhi_epi16(x3, x4);
            x2 = _mm256_unpacklo_epi16(x3, x4);
            x1 = _mm256_madd_epi16(x1, f1);
            t1 = _mm256_madd_epi16(t1, f1);
            x2 = _mm256_madd_epi16(x2, f2);
            t2 = _mm256_madd_epi16(t2, f2);
            x1 = _mm256_add_epi32(x1, x2);
            t1 = _mm256_add_epi32(t1, t2);
            t1 = _mm256_srai_epi32(t1, 2);
            x1 = _mm256_srai_epi32(x1, 2);
            x1 = _mm256_packs_epi32(x1, t1);
            {
                __m256i x3, x4;
                x3 = _mm256_mulhi_epi16(x1, wx);
                x1 = _mm256_mullo_epi16(x1, wx);
                x4 = _mm256_unpackhi_epi16(x1, x3);
                x1 = _mm256_unpacklo_epi16(x1, x3);
                x3 = _mm256_srai_epi32(_mm256_add_epi32(x4, offset), shift2);
                x1 = _mm256_srai_epi32(_mm256_add_epi32(x1, offset), shift2);
                x3 = _mm256_add_epi32(x3, ox);
                x1 = _mm256_add_epi32(x1, ox);
                x1 = _mm256_packus_epi32(x1, x3);
            };
            x1 = _mm256_max_epi16(x1, _mm256_setzero_si256());
            x1 = _mm256_min_epi16(x1, _mm256_set1_epi16(0x03FF));
            _mm256_storeu_si256((__m256i*)&dst[x], x1);
        }
        src += srcstride;
        dst += dststride;
    }
}

static void
put_vvc_bi_w_epel_h16_10_avx2(OVSample* _dst,
                              ptrdiff_t _dststride,
                              const OVSample* _src,
                              ptrdiff_t _srcstride,
                              const int16_t* src2,
                              ptrdiff_t src2stride,
                              int height,
                              int denom,
                              int _wx0,
                              int _wx1,
                              int offset0,
                              int offset1,
                              intptr_t mx,
                              intptr_t my,
                              int width)
{
    int x, y;
    uint16_t* src = ((uint16_t*)_src) - 1;
    ptrdiff_t srcstride = _srcstride;
    const int log2Wd = denom + 14 - BITDEPTH;
    const int shift2 = log2Wd + 1;

    const __m256i wx0 = _mm256_set1_epi16(_wx0);
    const __m256i wx1 = _mm256_set1_epi16(_wx1);
    const __m256i offset = _mm256_set1_epi32((offset0 + offset1 + 1) << log2Wd);
    uint16_t* dst = (uint16_t*)_dst;
    const int dststride = _dststride;
    __m256i f1 = _mm256_set1_epi32(ov_mcp_filters_c[mx - 1][0]);
    __m256i f2 = _mm256_set1_epi32(ov_mcp_filters_c[mx - 1][1]);
    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x += 16) {
            __m256i x1, x2, x3, x4, t1, t2;
            __m256i r5;
            x1 = _mm256_loadu_si256((__m256i*)&src[x]);
            x2 = _mm256_loadu_si256((__m256i*)&src[x + 1]);
            x3 = _mm256_loadu_si256((__m256i*)&src[x + 2 * 1]);
            x4 = _mm256_loadu_si256((__m256i*)&src[x + 3 * 1]);
            t1 = _mm256_unpackhi_epi16(x1, x2);
            x1 = _mm256_unpacklo_epi16(x1, x2);
            t2 = _mm256_unpackhi_epi16(x3, x4);
            x2 = _mm256_unpacklo_epi16(x3, x4);
            x1 = _mm256_madd_epi16(x1, f1);
            t1 = _mm256_madd_epi16(t1, f1);
            x2 = _mm256_madd_epi16(x2, f2);
            t2 = _mm256_madd_epi16(t2, f2);
            x1 = _mm256_add_epi32(x1, x2);
            t1 = _mm256_add_epi32(t1, t2);
            t1 = _mm256_srai_epi32(t1, 2);
            x1 = _mm256_srai_epi32(x1, 2);
            x1 = _mm256_packs_epi32(x1, t1);
            r5 = _mm256_load_si256((__m256i*)&src2[x]);
            {
                __m256i x3, x4, r7, r8;
                x3 = _mm256_mulhi_epi16(x1, wx1);
                x1 = _mm256_mullo_epi16(x1, wx1);
                r7 = _mm256_mulhi_epi16(r5, wx0);
                r5 = _mm256_mullo_epi16(r5, wx0);
                x4 = _mm256_unpackhi_epi16(x1, x3);
                x1 = _mm256_unpacklo_epi16(x1, x3);
                r8 = _mm256_unpackhi_epi16(r5, r7);
                r5 = _mm256_unpacklo_epi16(r5, r7);
                x4 = _mm256_add_epi32(x4, r8);
                x1 = _mm256_add_epi32(x1, r5);
                x4 = _mm256_srai_epi32(_mm256_add_epi32(x4, offset), shift2);
                x1 = _mm256_srai_epi32(_mm256_add_epi32(x1, offset), shift2);
                x1 = _mm256_packus_epi32(x1, x4);
            };
            x1 = _mm256_max_epi16(x1, _mm256_setzero_si256());
            x1 = _mm256_min_epi16(x1, _mm256_set1_epi16(0x03FF));
            _mm256_storeu_si256((__m256i*)&dst[x], x1);
        }
        src += srcstride;
        src2 += src2stride;
        dst += dststride;
    }
}

static void
put_vvc_uni_w_epel_v16_10_avx2(OVSample* _dst,
                               ptrdiff_t _dststride,
                               const OVSample* _src,
                               ptrdiff_t _srcstride,
                               int height,
                               int denom,
                               int _wx,
                               int _ox,
                               intptr_t mx,
                               intptr_t my,
                               int width)
{
    int x, y;
    ptrdiff_t srcstride = _srcstride;
    uint16_t* src = ((uint16_t*)_src) - srcstride;
    const int shift2 = denom + 14 - BITDEPTH;
    const __m256i ox = _mm256_set1_epi32(_ox);
    const __m256i wx = _mm256_set1_epi16(_wx);
    const __m256i offset = _mm256_set1_epi32(1 << (shift2 - 1));
    uint16_t* dst = (uint16_t*)_dst;
    const int dststride = _dststride;
    __m256i f1 = _mm256_set1_epi32(ov_mcp_filters_c[my - 1][0]);
    __m256i f2 = _mm256_set1_epi32(ov_mcp_filters_c[my - 1][1]);
    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x += 16) {
            __m256i x1, x2, x3, x4, t1, t2;
            x1 = _mm256_loadu_si256((__m256i*)&src[x]);
            x2 = _mm256_loadu_si256((__m256i*)&src[x + srcstride]);
            x3 = _mm256_loadu_si256((__m256i*)&src[x + 2 * srcstride]);
            x4 = _mm256_loadu_si256((__m256i*)&src[x + 3 * srcstride]);
            t1 = _mm256_unpackhi_epi16(x1, x2);
            x1 = _mm256_unpacklo_epi16(x1, x2);
            t2 = _mm256_unpackhi_epi16(x3, x4);
            x2 = _mm256_unpacklo_epi16(x3, x4);
            x1 = _mm256_madd_epi16(x1, f1);
            t1 = _mm256_madd_epi16(t1, f1);
            x2 = _mm256_madd_epi16(x2, f2);
            t2 = _mm256_madd_epi16(t2, f2);
            x1 = _mm256_add_epi32(x1, x2);
            t1 = _mm256_add_epi32(t1, t2);
            t1 = _mm256_srai_epi32(t1, 2);
            x1 = _mm256_srai_epi32(x1, 2);
            x1 = _mm256_packs_epi32(x1, t1);
            {
                __m256i x3, x4;
                x3 = _mm256_mulhi_epi16(x1, wx);
                x1 = _mm256_mullo_epi16(x1, wx);
                x4 = _mm256_unpackhi_epi16(x1, x3);
                x1 = _mm256_unpacklo_epi16(x1, x3);
                x3 = _mm256_srai_epi32(_mm256_add_epi32(x4, offset), shift2);
                x1 = _mm256_srai_epi32(_mm256_add_epi32(x1, offset), shift2);
                x3 = _mm256_add_epi32(x3, ox);
                x1 = _mm256_add_epi32(x1, ox);
                x1 = _mm256_packus_epi32(x1, x3);
            };
            x1 = _mm256_max_epi16(x1, _mm256_setzero_si256());
            x1 = _mm256_min_epi16(x1, _mm256_set1_epi16(0x03FF));
            _mm256_storeu_si256((__m256i*)&dst[x], x1);
        }
        src += srcstride;
        dst += dststride;
    }
}

static void
put_vvc_bi_w_epel_v16_10_avx2(OVSample* _dst,
                              ptrdiff_t _dststride,
                              const OVSample* _src,
                              ptrdiff_t _srcstride,
                              const int16_t* src2,
                              ptrdiff_t src2stride,
                              int height,
                              int denom,
                              int _wx0,
                              int _wx1,
                              int offset0,
                              int offset1,
                              intptr_t mx,
                              intptr_t my,
                              int width)
{
    int x, y;
    ptrdiff_t srcstride = _srcstride;
    uint16_t* src = ((uint16_t*)_src) - srcstride;
    const int log2Wd = denom + 14 - BITDEPTH;
    const int shift2 = log2Wd + 1;

    const __m256i wx0 = _mm256_set1_epi16(_wx0);
    const __m256i wx1 = _mm256_set1_epi16(_wx1);
    const __m256i offset = _mm256_set1_epi32((offset0 + offset1 + 1) << log2Wd);
    uint16_t* dst = (uint16_t*)_dst;
    const int dststride = _dststride;
    __m256i f1 = _mm256_set1_epi32(ov_mcp_filters_c[my - 1][0]);
    __m256i f2 = _mm256_set1_epi32(ov_mcp_filters_c[my - 1][1]);
    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x += 16) {
            __m256i x1, x2, x3, x4, t1, t2;
            __m256i r5;
            x1 = _mm256_loadu_si256((__m256i*)&src[x]);
            x2 = _mm256_loadu_si256((__m256i*)&src[x + srcstride]);
            x3 = _mm256_loadu_si256((__m256i*)&src[x + 2 * srcstride]);
            x4 = _mm256_loadu_si256((__m256i*)&src[x + 3 * srcstride]);
            t1 = _mm256_unpackhi_epi16(x1, x2);
            x1 = _mm256_unpacklo_epi16(x1, x2);
            t2 = _mm256_unpackhi_epi16(x3, x4);
            x2 = _mm256_unpacklo_epi16(x3, x4);
            x1 = _mm256_madd_epi16(x1, f1);
            t1 = _mm256_madd_epi16(t1, f1);
            x2 = _mm256_madd_epi16(x2, f2);
            t2 = _mm256_madd_epi16(t2, f2);
            x1 = _mm256_add_epi32(x1, x2);
            t1 = _mm256_add_epi32(t1, t2);
            t1 = _mm256_srai_epi32(t1, 2);
            x1 = _mm256_srai_epi32(x1, 2);
            x1 = _mm256_packs_epi32(x1, t1);
            r5 = _mm256_load_si256((__m256i*)&src2[x]);
            {
                __m256i x3, x4, r7, r8;
                x3 = _mm256_mulhi_epi16(x1, wx1);
                x1 = _mm256_mullo_epi16(x1, wx1);
                r7 = _mm256_mulhi_epi16(r5, wx0);
                r5 = _mm256_mullo_epi16(r5, wx0);
                x4 = _mm256_unpackhi_epi16(x1, x3);
                x1 = _mm256_unpacklo_epi16(x1, x3);
                r8 = _mm256_unpackhi_epi16(r5, r7);
                r5 = _mm256_unpacklo_epi16(r5, r7);
                x4 = _mm256_add_epi32(x4, r8);
                x1 = _mm256_add_epi32(x1, r5);
                x4 = _mm256_srai_epi32(_mm256_add_epi32(x4, offset), shift2);
                x1 = _mm256_srai_epi32(_mm256_add_epi32(x1, offset), shift2);
                x1 = _mm256_packus_epi32(x1, x4);
            };
            x1 = _mm256_max_epi16(x1, _mm256_setzero_si256());
            x1 = _mm256_min_epi16(x1, _mm256_set1_epi16(0x03FF));
            _mm256_storeu_si256((__m256i*)&dst[x], x1);
        }
        src += srcstride;
        src2 += src2stride;
        dst += dststride;
    }
}

static void
put_vvc_uni_w_epel_hv16_10_avx2(OVSample* _dst,
                                ptrdiff_t _dststride,
                                const OVSample* _src,
                                ptrdiff_t _srcstride,
                                int height,
                                int denom,
                                int _wx,
                                int _ox,
                                intptr_t mx,
                                intptr_t my,
                                int width)
{
    int x, y;
    const int shift2 = denom + 14 - BITDEPTH;
    const __m256i ox = _mm256_set1_epi32(_ox);
    const __m256i wx = _mm256_set1_epi16(_wx);
    const __m256i offset = _mm256_set1_epi32(1 << (shift2 - 1));
    uint16_t* dst = (uint16_t*)_dst;
    const int dststride = _dststride;
    uint16_t* dst_bis = dst;
    uint16_t *src_bis, *src = ((uint16_t*)_src) - 1;
    ptrdiff_t srcstride = _srcstride;
    src -= EPEL_EXTRA_BEFORE * srcstride;
    src_bis = src;
    __m256i f1 = _mm256_set1_epi32(ov_mcp_filters_c[mx - 1][0]);
    __m256i f2 = _mm256_set1_epi32(ov_mcp_filters_c[mx - 1][1]);
    __m256i f3 = _mm256_set1_epi32(ov_mcp_filters_c[my - 1][0]);
    __m256i f4 = _mm256_set1_epi32(ov_mcp_filters_c[my - 1][1]);
    for (x = 0; x < width; x += 16) {
        __m256i x1, x2, x3, x4, t1, t2, r1, r2, r3, r4;
        x1 = _mm256_loadu_si256((__m256i*)&src[x]);
        x2 = _mm256_loadu_si256((__m256i*)&src[x + 1]);
        x3 = _mm256_loadu_si256((__m256i*)&src[x + 2 * 1]);
        x4 = _mm256_loadu_si256((__m256i*)&src[x + 3 * 1]);
        t1 = _mm256_unpackhi_epi16(x1, x2);
        x1 = _mm256_unpacklo_epi16(x1, x2);
        t2 = _mm256_unpackhi_epi16(x3, x4);
        x2 = _mm256_unpacklo_epi16(x3, x4);
        x1 = _mm256_madd_epi16(x1, f1);
        t1 = _mm256_madd_epi16(t1, f1);
        x2 = _mm256_madd_epi16(x2, f2);
        t2 = _mm256_madd_epi16(t2, f2);
        x1 = _mm256_add_epi32(x1, x2);
        t1 = _mm256_add_epi32(t1, t2);
        t1 = _mm256_srai_epi32(t1, 2);
        x1 = _mm256_srai_epi32(x1, 2);
        r1 = _mm256_packs_epi32(x1, t1);
        src += srcstride;
        x1 = _mm256_loadu_si256((__m256i*)&src[x]);
        x2 = _mm256_loadu_si256((__m256i*)&src[x + 1]);
        x3 = _mm256_loadu_si256((__m256i*)&src[x + 2 * 1]);
        x4 = _mm256_loadu_si256((__m256i*)&src[x + 3 * 1]);
        t1 = _mm256_unpackhi_epi16(x1, x2);
        x1 = _mm256_unpacklo_epi16(x1, x2);
        t2 = _mm256_unpackhi_epi16(x3, x4);
        x2 = _mm256_unpacklo_epi16(x3, x4);
        x1 = _mm256_madd_epi16(x1, f1);
        t1 = _mm256_madd_epi16(t1, f1);
        x2 = _mm256_madd_epi16(x2, f2);
        t2 = _mm256_madd_epi16(t2, f2);
        x1 = _mm256_add_epi32(x1, x2);
        t1 = _mm256_add_epi32(t1, t2);
        t1 = _mm256_srai_epi32(t1, 2);
        x1 = _mm256_srai_epi32(x1, 2);
        r2 = _mm256_packs_epi32(x1, t1);
        src += srcstride;
        x1 = _mm256_loadu_si256((__m256i*)&src[x]);
        x2 = _mm256_loadu_si256((__m256i*)&src[x + 1]);
        x3 = _mm256_loadu_si256((__m256i*)&src[x + 2 * 1]);
        x4 = _mm256_loadu_si256((__m256i*)&src[x + 3 * 1]);
        t1 = _mm256_unpackhi_epi16(x1, x2);
        x1 = _mm256_unpacklo_epi16(x1, x2);
        t2 = _mm256_unpackhi_epi16(x3, x4);
        x2 = _mm256_unpacklo_epi16(x3, x4);
        x1 = _mm256_madd_epi16(x1, f1);
        t1 = _mm256_madd_epi16(t1, f1);
        x2 = _mm256_madd_epi16(x2, f2);
        t2 = _mm256_madd_epi16(t2, f2);
        x1 = _mm256_add_epi32(x1, x2);
        t1 = _mm256_add_epi32(t1, t2);
        t1 = _mm256_srai_epi32(t1, 2);
        x1 = _mm256_srai_epi32(x1, 2);
        r3 = _mm256_packs_epi32(x1, t1);
        src += srcstride;
        for (y = 0; y < height; y++) {
            x1 = _mm256_loadu_si256((__m256i*)&src[x]);
            x2 = _mm256_loadu_si256((__m256i*)&src[x + 1]);
            x3 = _mm256_loadu_si256((__m256i*)&src[x + 2 * 1]);
            x4 = _mm256_loadu_si256((__m256i*)&src[x + 3 * 1]);
            t1 = _mm256_unpackhi_epi16(x1, x2);
            x1 = _mm256_unpacklo_epi16(x1, x2);
            t2 = _mm256_unpackhi_epi16(x3, x4);
            x2 = _mm256_unpacklo_epi16(x3, x4);
            x1 = _mm256_madd_epi16(x1, f1);
            t1 = _mm256_madd_epi16(t1, f1);
            x2 = _mm256_madd_epi16(x2, f2);
            t2 = _mm256_madd_epi16(t2, f2);
            x1 = _mm256_add_epi32(x1, x2);
            t1 = _mm256_add_epi32(t1, t2);
            t1 = _mm256_srai_epi32(t1, 2);
            x1 = _mm256_srai_epi32(x1, 2);
            r4 = _mm256_packs_epi32(x1, t1);
            src += srcstride;
            t1 = _mm256_unpackhi_epi16(r1, r2);
            x1 = _mm256_unpacklo_epi16(r1, r2);
            t2 = _mm256_unpackhi_epi16(r3, r4);
            x2 = _mm256_unpacklo_epi16(r3, r4);
            x1 = _mm256_madd_epi16(x1, f3);
            t1 = _mm256_madd_epi16(t1, f3);
            x2 = _mm256_madd_epi16(x2, f4);
            t2 = _mm256_madd_epi16(t2, f4);
            x1 = _mm256_add_epi32(x1, x2);
            t1 = _mm256_add_epi32(t1, t2);
            t1 = _mm256_srai_epi32(t1, 6);
            x1 = _mm256_srai_epi32(x1, 6);
            x1 = _mm256_packs_epi32(x1, t1);
            {
                __m256i x3, x4;
                x3 = _mm256_mulhi_epi16(x1, wx);
                x1 = _mm256_mullo_epi16(x1, wx);
                x4 = _mm256_unpackhi_epi16(x1, x3);
                x1 = _mm256_unpacklo_epi16(x1, x3);
                x3 = _mm256_srai_epi32(_mm256_add_epi32(x4, offset), shift2);
                x1 = _mm256_srai_epi32(_mm256_add_epi32(x1, offset), shift2);
                x3 = _mm256_add_epi32(x3, ox);
                x1 = _mm256_add_epi32(x1, ox);
                x1 = _mm256_packus_epi32(x1, x3);
            };
            x1 = _mm256_max_epi16(x1, _mm256_setzero_si256());
            x1 = _mm256_min_epi16(x1, _mm256_set1_epi16(0x03FF));
            _mm256_storeu_si256((__m256i*)&dst[x], x1);
            dst += dststride;
            r1 = r2;
            r2 = r3;
            r3 = r4;
        }
        src = src_bis;
        dst = dst_bis;
    }
}

static void
put_vvc_bi_w_epel_hv16_10_avx2(OVSample* _dst,
                               ptrdiff_t _dststride,
                               const OVSample* _src,
                               ptrdiff_t _srcstride,
                               const int16_t* src2,
                               ptrdiff_t src2stride,
                               int height,
                               int denom,
                               int _wx0,
                               int _wx1,
                               int offset0,
                               int offset1,

                               intptr_t mx,
                               intptr_t my,
                               int width)
{
    int x, y;
    const int log2Wd = denom + 14 - BITDEPTH;
    const int shift2 = log2Wd + 1;

    const __m256i wx0 = _mm256_set1_epi16(_wx0);
    const __m256i wx1 = _mm256_set1_epi16(_wx1);
    const __m256i offset = _mm256_set1_epi32((offset0 + offset1 + 1) << log2Wd);
    uint16_t* dst = (uint16_t*)_dst;
    const int dststride = _dststride;
    const int16_t* src2_bis = src2;
    uint16_t* dst_bis = dst;
    const uint16_t *src_bis, *src = ((uint16_t*)_src) - 1;
    ptrdiff_t srcstride = _srcstride;
    src -= EPEL_EXTRA_BEFORE * srcstride;
    src_bis = src;
    __m256i f1 = _mm256_set1_epi32(ov_mcp_filters_c[mx - 1][0]);
    __m256i f2 = _mm256_set1_epi32(ov_mcp_filters_c[mx - 1][1]);
    __m256i f3 = _mm256_set1_epi32(ov_mcp_filters_c[my - 1][0]);
    __m256i f4 = _mm256_set1_epi32(ov_mcp_filters_c[my - 1][1]);
    for (x = 0; x < width; x += 16) {
        __m256i x1, x2, x3, x4, t1, t2, r1, r2, r3, r4;
        __m256i r5;
        x1 = _mm256_loadu_si256((__m256i*)&src[x]);
        x2 = _mm256_loadu_si256((__m256i*)&src[x + 1]);
        x3 = _mm256_loadu_si256((__m256i*)&src[x + 2 * 1]);
        x4 = _mm256_loadu_si256((__m256i*)&src[x + 3 * 1]);
        t1 = _mm256_unpackhi_epi16(x1, x2);
        x1 = _mm256_unpacklo_epi16(x1, x2);
        t2 = _mm256_unpackhi_epi16(x3, x4);
        x2 = _mm256_unpacklo_epi16(x3, x4);
        x1 = _mm256_madd_epi16(x1, f1);
        t1 = _mm256_madd_epi16(t1, f1);
        x2 = _mm256_madd_epi16(x2, f2);
        t2 = _mm256_madd_epi16(t2, f2);
        x1 = _mm256_add_epi32(x1, x2);
        t1 = _mm256_add_epi32(t1, t2);
        t1 = _mm256_srai_epi32(t1, 2);
        x1 = _mm256_srai_epi32(x1, 2);
        r1 = _mm256_packs_epi32(x1, t1);
        src += srcstride;
        x1 = _mm256_loadu_si256((__m256i*)&src[x]);
        x2 = _mm256_loadu_si256((__m256i*)&src[x + 1]);
        x3 = _mm256_loadu_si256((__m256i*)&src[x + 2 * 1]);
        x4 = _mm256_loadu_si256((__m256i*)&src[x + 3 * 1]);
        t1 = _mm256_unpackhi_epi16(x1, x2);
        x1 = _mm256_unpacklo_epi16(x1, x2);
        t2 = _mm256_unpackhi_epi16(x3, x4);
        x2 = _mm256_unpacklo_epi16(x3, x4);
        x1 = _mm256_madd_epi16(x1, f1);
        t1 = _mm256_madd_epi16(t1, f1);
        x2 = _mm256_madd_epi16(x2, f2);
        t2 = _mm256_madd_epi16(t2, f2);
        x1 = _mm256_add_epi32(x1, x2);
        t1 = _mm256_add_epi32(t1, t2);
        t1 = _mm256_srai_epi32(t1, 2);
        x1 = _mm256_srai_epi32(x1, 2);
        r2 = _mm256_packs_epi32(x1, t1);
        src += srcstride;
        x1 = _mm256_loadu_si256((__m256i*)&src[x]);
        x2 = _mm256_loadu_si256((__m256i*)&src[x + 1]);
        x3 = _mm256_loadu_si256((__m256i*)&src[x + 2 * 1]);
        x4 = _mm256_loadu_si256((__m256i*)&src[x + 3 * 1]);
        t1 = _mm256_unpackhi_epi16(x1, x2);
        x1 = _mm256_unpacklo_epi16(x1, x2);
        t2 = _mm256_unpackhi_epi16(x3, x4);
        x2 = _mm256_unpacklo_epi16(x3, x4);
        x1 = _mm256_madd_epi16(x1, f1);
        t1 = _mm256_madd_epi16(t1, f1);
        x2 = _mm256_madd_epi16(x2, f2);
        t2 = _mm256_madd_epi16(t2, f2);
        x1 = _mm256_add_epi32(x1, x2);
        t1 = _mm256_add_epi32(t1, t2);
        t1 = _mm256_srai_epi32(t1, 2);
        x1 = _mm256_srai_epi32(x1, 2);
        r3 = _mm256_packs_epi32(x1, t1);
        src += srcstride;
        for (y = 0; y < height; y++) {
            x1 = _mm256_loadu_si256((__m256i*)&src[x]);
            x2 = _mm256_loadu_si256((__m256i*)&src[x + 1]);
            x3 = _mm256_loadu_si256((__m256i*)&src[x + 2 * 1]);
            x4 = _mm256_loadu_si256((__m256i*)&src[x + 3 * 1]);
            t1 = _mm256_unpackhi_epi16(x1, x2);
            x1 = _mm256_unpacklo_epi16(x1, x2);
            t2 = _mm256_unpackhi_epi16(x3, x4);
            x2 = _mm256_unpacklo_epi16(x3, x4);
            x1 = _mm256_madd_epi16(x1, f1);
            t1 = _mm256_madd_epi16(t1, f1);
            x2 = _mm256_madd_epi16(x2, f2);
            t2 = _mm256_madd_epi16(t2, f2);
            x1 = _mm256_add_epi32(x1, x2);
            t1 = _mm256_add_epi32(t1, t2);
            t1 = _mm256_srai_epi32(t1, 2);
            x1 = _mm256_srai_epi32(x1, 2);
            r4 = _mm256_packs_epi32(x1, t1);
            src += srcstride;
            t1 = _mm256_unpackhi_epi16(r1, r2);
            x1 = _mm256_unpacklo_epi16(r1, r2);
            t2 = _mm256_unpackhi_epi16(r3, r4);
            x2 = _mm256_unpacklo_epi16(r3, r4);
            x1 = _mm256_madd_epi16(x1, f3);
            t1 = _mm256_madd_epi16(t1, f3);
            x2 = _mm256_madd_epi16(x2, f4);
            t2 = _mm256_madd_epi16(t2, f4);
            x1 = _mm256_add_epi32(x1, x2);
            t1 = _mm256_add_epi32(t1, t2);
            t1 = _mm256_srai_epi32(t1, 6);
            x1 = _mm256_srai_epi32(x1, 6);
            x1 = _mm256_packs_epi32(x1, t1);
            r5 = _mm256_load_si256((__m256i*)&src2[x]);
            {
                __m256i x3, x4, r7, r8;
                x3 = _mm256_mulhi_epi16(x1, wx1);
                x1 = _mm256_mullo_epi16(x1, wx1);
                r7 = _mm256_mulhi_epi16(r5, wx0);
                r5 = _mm256_mullo_epi16(r5, wx0);
                x4 = _mm256_unpackhi_epi16(x1, x3);
                x1 = _mm256_unpacklo_epi16(x1, x3);
                r8 = _mm256_unpackhi_epi16(r5, r7);
                r5 = _mm256_unpacklo_epi16(r5, r7);
                x4 = _mm256_add_epi32(x4, r8);
                x1 = _mm256_add_epi32(x1, r5);
                x4 = _mm256_srai_epi32(_mm256_add_epi32(x4, offset), shift2);
                x1 = _mm256_srai_epi32(_mm256_add_epi32(x1, offset), shift2);
                x1 = _mm256_packus_epi32(x1, x4);
            };
            x1 = _mm256_max_epi16(x1, _mm256_setzero_si256());
            x1 = _mm256_min_epi16(x1, _mm256_set1_epi16(0x03FF));
            _mm256_storeu_si256((__m256i*)&dst[x], x1);
            src2 += src2stride;
            dst += dststride;
            r1 = r2;
            r2 = r3;
            r3 = r4;
        }
        src = src_bis;
        src2 = src2_bis;
        dst = dst_bis;
    }
}

static void
put_vvc_uni_w_qpel_h16_10_avx2(OVSample* _dst,
                               ptrdiff_t _dststride,
                               const OVSample* _src,
                               ptrdiff_t _srcstride,
                               int height,
                               int denom,
                               int _wx,
                               int _ox,
                               intptr_t mx,
                               intptr_t my,
                               int width)
{
    int x, y;
    int shift = BITDEPTH - 8;
    uint16_t* src = (uint16_t*)_src;
    const int srcstride = _srcstride;
    const int shift2 = denom + 14 - BITDEPTH;
    const __m256i ox = _mm256_set1_epi32(_ox);
    const __m256i wx = _mm256_set1_epi16(_wx);
    const __m256i offset = _mm256_set1_epi32(1 << (shift2 - 1));
    uint16_t* dst = (uint16_t*)_dst;
    const int dststride = _dststride;

    __m256i c1 = _mm256_set1_epi32(ov_mcp_filters_l[mx - 1][0]);
    __m256i c2 = _mm256_set1_epi32(ov_mcp_filters_l[mx - 1][1]);
    __m256i c3 = _mm256_set1_epi32(ov_mcp_filters_l[mx - 1][2]);
    __m256i c4 = _mm256_set1_epi32(ov_mcp_filters_l[mx - 1][3]);

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x += 16) {
            __m256i x1, x2, x3, x4, r1, r2, r3, r4, t1, t2;
            x1 = _mm256_loadu_si256((__m256i*)&src[x - 3 * 1]);
            x2 = _mm256_loadu_si256((__m256i*)&src[x - 2 * 1]);
            x3 = _mm256_loadu_si256((__m256i*)&src[x - 1 * 1]);
            x4 = _mm256_loadu_si256((__m256i*)&src[x]);
            t1 = _mm256_unpackhi_epi16(x1, x2);
            x1 = _mm256_unpacklo_epi16(x1, x2);
            t2 = _mm256_unpackhi_epi16(x3, x4);
            x2 = _mm256_unpacklo_epi16(x3, x4);
            t1 = _mm256_madd_epi16(t1, c1);
            t2 = _mm256_madd_epi16(t2, c2);
            r2 = _mm256_add_epi32(t1, t2);
            x1 = _mm256_madd_epi16(x1, c1);
            x2 = _mm256_madd_epi16(x2, c2);
            r1 = _mm256_add_epi32(x1, x2);
            x1 = _mm256_loadu_si256((__m256i*)&src[x + 1]);
            x2 = _mm256_loadu_si256((__m256i*)&src[x + 2 * 1]);
            x3 = _mm256_loadu_si256((__m256i*)&src[x + 3 * 1]);
            x4 = _mm256_loadu_si256((__m256i*)&src[x + 4 * 1]);
            t1 = _mm256_unpackhi_epi16(x1, x2);
            x1 = _mm256_unpacklo_epi16(x1, x2);
            t2 = _mm256_unpackhi_epi16(x3, x4);
            x2 = _mm256_unpacklo_epi16(x3, x4);
            t1 = _mm256_madd_epi16(t1, c3);
            t2 = _mm256_madd_epi16(t2, c4);
            r4 = _mm256_add_epi32(t1, t2);
            x1 = _mm256_madd_epi16(x1, c3);
            x2 = _mm256_madd_epi16(x2, c4);
            r3 = _mm256_add_epi32(x1, x2);
            x1 = _mm256_add_epi32(r1, r3);
            x2 = _mm256_add_epi32(r2, r4);
            x1 = _mm256_srai_epi32(x1, shift);
            x2 = _mm256_srai_epi32(x2, shift);
            x1 = _mm256_packs_epi32(x1, x2);
            {
                __m256i x3, x4;
                x3 = _mm256_mulhi_epi16(x1, wx);
                x1 = _mm256_mullo_epi16(x1, wx);
                x4 = _mm256_unpackhi_epi16(x1, x3);
                x1 = _mm256_unpacklo_epi16(x1, x3);
                x3 = _mm256_srai_epi32(_mm256_add_epi32(x4, offset), shift2);
                x1 = _mm256_srai_epi32(_mm256_add_epi32(x1, offset), shift2);
                x3 = _mm256_add_epi32(x3, ox);
                x1 = _mm256_add_epi32(x1, ox);
                x1 = _mm256_packus_epi32(x1, x3);
            };
            x1 = _mm256_max_epi16(x1, _mm256_setzero_si256());
            x1 = _mm256_min_epi16(x1, _mm256_set1_epi16(0x03FF));
            _mm256_storeu_si256((__m256i*)&dst[x], x1);
        }
        src += srcstride;
        dst += dststride;
    }
}

static void
put_vvc_bi_w_qpel_h16_10_avx2(OVSample* _dst,
                              ptrdiff_t _dststride,
                              const OVSample* _src,
                              ptrdiff_t _srcstride,
                              const int16_t* src2,
                              ptrdiff_t src2stride,
                              int height,
                              int denom,
                              int _wx0,
                              int _wx1,
                              int offset0,
                              int offset1,
                              intptr_t mx,
                              intptr_t my,
                              int width)
{
    int x, y;
    int shift = BITDEPTH - 8;
    uint16_t* src = (uint16_t*)_src;
    const int srcstride = _srcstride;
    const int log2Wd = denom + 14 - BITDEPTH;
    const int shift2 = log2Wd + 1;

    const __m256i wx0 = _mm256_set1_epi16(_wx0);
    const __m256i wx1 = _mm256_set1_epi16(_wx1);
    const __m256i offset = _mm256_set1_epi32((offset0 + offset1 + 1) << log2Wd);
    uint16_t* dst = (uint16_t*)_dst;
    const int dststride = _dststride;

    __m256i c1 = _mm256_set1_epi32(ov_mcp_filters_l[mx - 1][0]);
    __m256i c2 = _mm256_set1_epi32(ov_mcp_filters_l[mx - 1][1]);
    __m256i c3 = _mm256_set1_epi32(ov_mcp_filters_l[mx - 1][2]);
    __m256i c4 = _mm256_set1_epi32(ov_mcp_filters_l[mx - 1][3]);

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x += 16) {
            __m256i x1, x2, x3, x4, r1, r2, r3, r4, t1, t2;
            __m256i r5;
            x1 = _mm256_loadu_si256((__m256i*)&src[x - 3 * 1]);
            x2 = _mm256_loadu_si256((__m256i*)&src[x - 2 * 1]);
            x3 = _mm256_loadu_si256((__m256i*)&src[x - 1 * 1]);
            x4 = _mm256_loadu_si256((__m256i*)&src[x]);
            t1 = _mm256_unpackhi_epi16(x1, x2);
            x1 = _mm256_unpacklo_epi16(x1, x2);
            t2 = _mm256_unpackhi_epi16(x3, x4);
            x2 = _mm256_unpacklo_epi16(x3, x4);
            t1 = _mm256_madd_epi16(t1, c1);
            t2 = _mm256_madd_epi16(t2, c2);
            r2 = _mm256_add_epi32(t1, t2);
            x1 = _mm256_madd_epi16(x1, c1);
            x2 = _mm256_madd_epi16(x2, c2);
            r1 = _mm256_add_epi32(x1, x2);
            x1 = _mm256_loadu_si256((__m256i*)&src[x + 1]);
            x2 = _mm256_loadu_si256((__m256i*)&src[x + 2 * 1]);
            x3 = _mm256_loadu_si256((__m256i*)&src[x + 3 * 1]);
            x4 = _mm256_loadu_si256((__m256i*)&src[x + 4 * 1]);
            t1 = _mm256_unpackhi_epi16(x1, x2);
            x1 = _mm256_unpacklo_epi16(x1, x2);
            t2 = _mm256_unpackhi_epi16(x3, x4);
            x2 = _mm256_unpacklo_epi16(x3, x4);
            t1 = _mm256_madd_epi16(t1, c3);
            t2 = _mm256_madd_epi16(t2, c4);
            r4 = _mm256_add_epi32(t1, t2);
            x1 = _mm256_madd_epi16(x1, c3);
            x2 = _mm256_madd_epi16(x2, c4);
            r3 = _mm256_add_epi32(x1, x2);
            x1 = _mm256_add_epi32(r1, r3);
            x2 = _mm256_add_epi32(r2, r4);
            x1 = _mm256_srai_epi32(x1, shift);
            x2 = _mm256_srai_epi32(x2, shift);
            x1 = _mm256_packs_epi32(x1, x2);
            r5 = _mm256_load_si256((__m256i*)&src2[x]);
            {
                __m256i x3, x4, r7, r8;
                x3 = _mm256_mulhi_epi16(x1, wx1);
                x1 = _mm256_mullo_epi16(x1, wx1);
                r7 = _mm256_mulhi_epi16(r5, wx0);
                r5 = _mm256_mullo_epi16(r5, wx0);
                x4 = _mm256_unpackhi_epi16(x1, x3);
                x1 = _mm256_unpacklo_epi16(x1, x3);
                r8 = _mm256_unpackhi_epi16(r5, r7);
                r5 = _mm256_unpacklo_epi16(r5, r7);
                x4 = _mm256_add_epi32(x4, r8);
                x1 = _mm256_add_epi32(x1, r5);
                x4 = _mm256_srai_epi32(_mm256_add_epi32(x4, offset), shift2);
                x1 = _mm256_srai_epi32(_mm256_add_epi32(x1, offset), shift2);
                x1 = _mm256_packus_epi32(x1, x4);
            };
            x1 = _mm256_max_epi16(x1, _mm256_setzero_si256());
            x1 = _mm256_min_epi16(x1, _mm256_set1_epi16(0x03FF));
            _mm256_storeu_si256((__m256i*)&dst[x], x1);
        }
        src += srcstride;
        src2 += src2stride;
        dst += dststride;
    }
}

static void
put_vvc_uni_w_qpel_v16_10_avx2(OVSample* _dst,
                               ptrdiff_t _dststride,
                               const OVSample* _src,
                               ptrdiff_t _srcstride,
                               int height,
                               int denom,
                               int _wx,
                               int _ox,
                               intptr_t mx,
                               intptr_t my,
                               int width)
{
    int x, y;
    uint16_t* src = (uint16_t*)_src;
    const int srcstride = _srcstride;
    const int shift2 = denom + 14 - BITDEPTH;
    const __m256i ox = _mm256_set1_epi32(_ox);
    const __m256i wx = _mm256_set1_epi16(_wx);
    const __m256i offset = _mm256_set1_epi32(1 << (shift2 - 1));
    uint16_t* dst = (uint16_t*)_dst;
    const int dststride = _dststride;

    __m256i c1 = _mm256_set1_epi32(ov_mcp_filters_l[my - 1][0]);
    __m256i c2 = _mm256_set1_epi32(ov_mcp_filters_l[my - 1][1]);
    __m256i c3 = _mm256_set1_epi32(ov_mcp_filters_l[my - 1][2]);
    __m256i c4 = _mm256_set1_epi32(ov_mcp_filters_l[my - 1][3]);

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x += 16) {
            __m256i x1, x2, x3, x4, x5, x6, x7, x8, x9;
            x1 = _mm256_loadu_si256((__m256i*)&src[x - 3 * srcstride]);
            x2 = _mm256_loadu_si256((__m256i*)&src[x - 2 * srcstride]);
            x3 = _mm256_loadu_si256((__m256i*)&src[x - 1 * srcstride]);
            x4 = _mm256_loadu_si256((__m256i*)&src[x]);
            x5 = _mm256_loadu_si256((__m256i*)&src[x + srcstride]);
            x6 = _mm256_loadu_si256((__m256i*)&src[x + 2 * srcstride]);
            x7 = _mm256_loadu_si256((__m256i*)&src[x + 3 * srcstride]);
            x8 = _mm256_loadu_si256((__m256i*)&src[x + 4 * srcstride]);
            x9 = x1;
            x1 = _mm256_unpacklo_epi16(x9, x2);
            x2 = _mm256_unpackhi_epi16(x9, x2);
            x9 = x3;
            x3 = _mm256_unpacklo_epi16(x9, x4);
            x4 = _mm256_unpackhi_epi16(x9, x4);
            x9 = x5;
            x5 = _mm256_unpacklo_epi16(x9, x6);
            x6 = _mm256_unpackhi_epi16(x9, x6);
            x9 = x7;
            x7 = _mm256_unpacklo_epi16(x9, x8);
            x8 = _mm256_unpackhi_epi16(x9, x8);
            x1 = _mm256_madd_epi16(x1, c1);
            x3 = _mm256_madd_epi16(x3, c2);
            x5 = _mm256_madd_epi16(x5, c3);
            x7 = _mm256_madd_epi16(x7, c4);
            x2 = _mm256_madd_epi16(x2, c1);
            x4 = _mm256_madd_epi16(x4, c2);
            x6 = _mm256_madd_epi16(x6, c3);
            x8 = _mm256_madd_epi16(x8, c4);
            x1 = _mm256_add_epi32(x1, x3);
            x3 = _mm256_add_epi32(x5, x7);
            x2 = _mm256_add_epi32(x2, x4);
            x4 = _mm256_add_epi32(x6, x8);
            x1 = _mm256_add_epi32(x1, x3);
            x2 = _mm256_add_epi32(x2, x4);
            x1 = _mm256_srai_epi32(x1, 2);
            x2 = _mm256_srai_epi32(x2, 2);
            x1 = _mm256_packs_epi32(x1, x2);
            {
                __m256i x3, x4;
                x3 = _mm256_mulhi_epi16(x1, wx);
                x1 = _mm256_mullo_epi16(x1, wx);
                x4 = _mm256_unpackhi_epi16(x1, x3);
                x1 = _mm256_unpacklo_epi16(x1, x3);
                x3 = _mm256_srai_epi32(_mm256_add_epi32(x4, offset), shift2);
                x1 = _mm256_srai_epi32(_mm256_add_epi32(x1, offset), shift2);
                x3 = _mm256_add_epi32(x3, ox);
                x1 = _mm256_add_epi32(x1, ox);
                x1 = _mm256_packus_epi32(x1, x3);
            };
            x1 = _mm256_max_epi16(x1, _mm256_setzero_si256());
            x1 = _mm256_min_epi16(x1, _mm256_set1_epi16(0x03FF));
            _mm256_storeu_si256((__m256i*)&dst[x], x1);
        }
        src += srcstride;
        dst += dststride;
    }
}

static void
put_vvc_bi_w_qpel_v16_10_avx2(OVSample* _dst,
                              ptrdiff_t _dststride,
                              const OVSample* _src,
                              ptrdiff_t _srcstride,
                              const int16_t* src2,
                              ptrdiff_t src2stride,
                              int height,
                              int denom,
                              int _wx0,
                              int _wx1,
                              int offset0,
                              int offset1,
                              intptr_t mx,
                              intptr_t my,
                              int width)
{
    int x, y;
    uint16_t* src = (uint16_t*)_src;
    const int srcstride = _srcstride;
    const int log2Wd = denom + 14 - BITDEPTH;
    const int shift2 = log2Wd + 1;
    const __m256i wx0 = _mm256_set1_epi16(_wx0);
    const __m256i wx1 = _mm256_set1_epi16(_wx1);
    const __m256i offset = _mm256_set1_epi32((offset0 + offset1 + 1) << log2Wd);
    uint16_t* dst = (uint16_t*)_dst;
    const int dststride = _dststride;

    __m256i c1 = _mm256_set1_epi32(ov_mcp_filters_l[my - 1][0]);
    __m256i c2 = _mm256_set1_epi32(ov_mcp_filters_l[my - 1][1]);
    __m256i c3 = _mm256_set1_epi32(ov_mcp_filters_l[my - 1][2]);
    __m256i c4 = _mm256_set1_epi32(ov_mcp_filters_l[my - 1][3]);

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x += 16) {
            __m256i x1, x2, x3, x4, x5, x6, x7, x8, x9;
            __m256i r5;
            x1 = _mm256_loadu_si256((__m256i*)&src[x - 3 * srcstride]);
            x2 = _mm256_loadu_si256((__m256i*)&src[x - 2 * srcstride]);
            x3 = _mm256_loadu_si256((__m256i*)&src[x - 1 * srcstride]);
            x4 = _mm256_loadu_si256((__m256i*)&src[x]);
            x5 = _mm256_loadu_si256((__m256i*)&src[x + srcstride]);
            x6 = _mm256_loadu_si256((__m256i*)&src[x + 2 * srcstride]);
            x7 = _mm256_loadu_si256((__m256i*)&src[x + 3 * srcstride]);
            x8 = _mm256_loadu_si256((__m256i*)&src[x + 4 * srcstride]);
            x9 = x1;
            x1 = _mm256_unpacklo_epi16(x9, x2);
            x2 = _mm256_unpackhi_epi16(x9, x2);
            x9 = x3;
            x3 = _mm256_unpacklo_epi16(x9, x4);
            x4 = _mm256_unpackhi_epi16(x9, x4);
            x9 = x5;
            x5 = _mm256_unpacklo_epi16(x9, x6);
            x6 = _mm256_unpackhi_epi16(x9, x6);
            x9 = x7;
            x7 = _mm256_unpacklo_epi16(x9, x8);
            x8 = _mm256_unpackhi_epi16(x9, x8);
            x1 = _mm256_madd_epi16(x1, c1);
            x3 = _mm256_madd_epi16(x3, c2);
            x5 = _mm256_madd_epi16(x5, c3);
            x7 = _mm256_madd_epi16(x7, c4);
            x2 = _mm256_madd_epi16(x2, c1);
            x4 = _mm256_madd_epi16(x4, c2);
            x6 = _mm256_madd_epi16(x6, c3);
            x8 = _mm256_madd_epi16(x8, c4);
            x1 = _mm256_add_epi32(x1, x3);
            x3 = _mm256_add_epi32(x5, x7);
            x2 = _mm256_add_epi32(x2, x4);
            x4 = _mm256_add_epi32(x6, x8);
            x1 = _mm256_add_epi32(x1, x3);
            x2 = _mm256_add_epi32(x2, x4);
            x1 = _mm256_srai_epi32(x1, 2);
            x2 = _mm256_srai_epi32(x2, 2);
            x1 = _mm256_packs_epi32(x1, x2);
            r5 = _mm256_load_si256((__m256i*)&src2[x]);
            {
                __m256i x3, x4, r7, r8;
                x3 = _mm256_mulhi_epi16(x1, wx1);
                x1 = _mm256_mullo_epi16(x1, wx1);
                r7 = _mm256_mulhi_epi16(r5, wx0);
                r5 = _mm256_mullo_epi16(r5, wx0);
                x4 = _mm256_unpackhi_epi16(x1, x3);
                x1 = _mm256_unpacklo_epi16(x1, x3);
                r8 = _mm256_unpackhi_epi16(r5, r7);
                r5 = _mm256_unpacklo_epi16(r5, r7);
                x4 = _mm256_add_epi32(x4, r8);
                x1 = _mm256_add_epi32(x1, r5);
                x4 = _mm256_srai_epi32(_mm256_add_epi32(x4, offset), shift2);
                x1 = _mm256_srai_epi32(_mm256_add_epi32(x1, offset), shift2);
                x1 = _mm256_packus_epi32(x1, x4);
            };
            x1 = _mm256_max_epi16(x1, _mm256_setzero_si256());
            x1 = _mm256_min_epi16(x1, _mm256_set1_epi16(0x03FF));
            _mm256_storeu_si256((__m256i*)&dst[x], x1);
        }
        src += srcstride;
        src2 += src2stride;
        dst += dststride;
    }
}

static void
put_vvc_uni_w_qpel_hv16_10_avx2(OVSample* dst,
                                ptrdiff_t dststride,
                                const OVSample* _src,
                                ptrdiff_t _srcstride,
                                int height,
                                int denom,
                                int wx,
                                int ox,
                                intptr_t mx,
                                intptr_t my,
                                int width)
{
    int16_t tmp_array[(MAX_PB_SIZE + QPEL_EXTRA) * MAX_PB_SIZE];
    int16_t* tmp = tmp_array;
    uint16_t* src = (uint16_t*)_src;
    const int srcstride = _srcstride;
    src -= QPEL_EXTRA_BEFORE * srcstride;
    put_vvc_qpel_h16_10_avx2(tmp,
                             MAX_PB_SIZE,
                             src,
                             _srcstride,
                             height + QPEL_EXTRA,
                             mx,
                             my,
                             width);
    tmp = tmp_array + QPEL_EXTRA_BEFORE * MAX_PB_SIZE;
    put_vvc_uni_w_qpel_v16_14_10_avx2(dst,
                                      dststride,
                                      tmp,
                                      MAX_PB_SIZE,
                                      height,
                                      denom,
                                      wx,
                                      ox,
                                      mx,
                                      my,
                                      width);
}

static void
put_vvc_bi_w_qpel_hv16_10_avx2(OVSample* dst,
                               ptrdiff_t dststride,
                               const OVSample* _src,
                               ptrdiff_t _srcstride,
                               const int16_t* src2,
                               ptrdiff_t src2stride,
                               int height,
                               int denom,
                               int wx0,
                               int wx1,
                               int offset0,
                               int offset1,
                               intptr_t mx,
                               intptr_t my,
                               int width)
{
    int16_t tmp_array[(MAX_PB_SIZE + QPEL_EXTRA) * MAX_PB_SIZE];
    int16_t* tmp = tmp_array;
    uint16_t* src = (uint16_t*)_src;
    const int srcstride = _srcstride;
    src -= QPEL_EXTRA_BEFORE * srcstride;
    put_vvc_qpel_h16_10_avx2(tmp,
                             MAX_PB_SIZE,
                             src,
                             _srcstride,
                             height + QPEL_EXTRA,
                             mx,
                             my,
                             width);
    tmp = tmp_array + QPEL_EXTRA_BEFORE * MAX_PB_SIZE;
    put_vvc_bi_w_qpel_v16_14_10_avx2(dst,
                                     dststride,
                                     tmp,
                                     MAX_PB_SIZE,
                                     src2,
                                     src2stride,
                                     height,
                                     denom,
                                     wx0,
                                     wx1,
                                     offset0,
                                     offset1,
                                     mx,
                                     my,
                                     width);
}

static void
put_vvc_uni_w_pel_pixels16_10_avx2(OVSample* dst,
                                   ptrdiff_t dststride,
                                   const OVSample* _src,
                                   ptrdiff_t _srcstride,
                                   int height,
                                   int denom,
                                   int wx,
                                   int ox,
                                   intptr_t mx,
                                   intptr_t my,
                                   int width)
{
    put_vvc_uni_w_pel_pixels8_10_avx2(dst, dststride, _src, _srcstride, height, denom, wx, ox, mx, my, width);
}

static void
put_vvc_bi_w_pel_pixels16_10_avx2(OVSample* dst,
                                  ptrdiff_t dststride,
                                  const OVSample* _src,
                                  ptrdiff_t _srcstride,
                                  const int16_t* src2,
                                  ptrdiff_t src2stride,
                                  int height,
                                  int denom,
                                  int wx0,
                                  int wx1,
                                  int offset0,
                                  int offset1,
                                  intptr_t mx,
                                  intptr_t my,
                                  int width)
{
    put_vvc_bi_w_pel_pixels8_10_avx2(dst,
                                     dststride,
                                     _src,
                                     _srcstride,
                                     src2,
                                     src2stride,
                                     height,
                                     denom,
                                     wx0,
                                     wx1,
                                     offset0,
                                     offset1,
                                     mx,
                                     my,
                                     width);
}

static void
put_vvc_uni_w_pel_pixels32_10_avx2(OVSample* dst,
                                   ptrdiff_t dststride,
                                   const OVSample* _src,
                                   ptrdiff_t _srcstride,
                                   int height,
                                   int denom,
                                   int wx,
                                   int ox,
                                   intptr_t mx,
                                   intptr_t my,
                                   int width)
{
    put_vvc_uni_w_pel_pixels8_10_avx2(dst, dststride, _src, _srcstride, height, denom, wx, ox, mx, my, width);
}

static void
put_vvc_bi_w_pel_pixels32_10_avx2(OVSample* dst,
                                  ptrdiff_t dststride,
                                  const OVSample* _src,
                                  ptrdiff_t _srcstride,
                                  const int16_t* src2,
                                  ptrdiff_t src2stride,
                                  int height,
                                  int denom,
                                  int wx0,
                                  int wx1,
                                  int offset0,
                                  int offset1,
                                  intptr_t mx,
                                  intptr_t my,
                                  int width)
{
    put_vvc_bi_w_pel_pixels8_10_avx2(dst,
                                     dststride,
                                     _src,
                                     _srcstride,
                                     src2,
                                     src2stride,
                                     height,
                                     denom,
                                     wx0,
                                     wx1,
                                     offset0,
                                     offset1,
                                     mx,
                                     my,
                                     width);
}

static void
put_vvc_uni_w_pel_pixels64_10_avx2(OVSample* dst,
                                   ptrdiff_t dststride,
                                   const OVSample* _src,
                                   ptrdiff_t _srcstride,
                                   int height,
                                   int denom,
                                   int wx,
                                   int ox,
                                   intptr_t mx,
                                   intptr_t my,
                                   int width)
{
    put_vvc_uni_w_pel_pixels8_10_avx2(dst, dststride, _src, _srcstride, height, denom, wx, ox, mx, my, width);
}

static void
put_vvc_bi_w_pel_pixels64_10_avx2(OVSample* dst,
                                  ptrdiff_t dststride,
                                  const OVSample* _src,
                                  ptrdiff_t _srcstride,
                                  const int16_t* src2,
                                  ptrdiff_t src2stride,
                                  int height,
                                  int denom,
                                  int wx0,
                                  int wx1,
                                  int offset0,
                                  int offset1,
                                  intptr_t mx,
                                  intptr_t my,
                                  int width)
{
    put_vvc_bi_w_pel_pixels8_10_avx2(dst,
                                     dststride,
                                     _src,
                                     _srcstride,
                                     src2,
                                     src2stride,
                                     height,
                                     denom,
                                     wx0,
                                     wx1,
                                     offset0,
                                     offset1,
                                     mx,
                                     my,
                                     width);
}

static void
put_vvc_uni_w_qpel_h32_10_avx2(OVSample* dst,
                               ptrdiff_t dststride,
                               const OVSample* _src,
                               ptrdiff_t _srcstride,
                               int height,
                               int denom,
                               int wx,
                               int ox,
                               intptr_t mx,
                               intptr_t my,
                               int width)
{
    put_vvc_uni_w_qpel_h16_10_avx2(dst, dststride, _src, _srcstride, height, denom, wx, ox, mx, my, width);
}

static void
put_vvc_bi_w_qpel_h32_10_avx2(OVSample* dst,
                              ptrdiff_t dststride,
                              const OVSample* _src,
                              ptrdiff_t _srcstride,
                              const int16_t* src2,
                              ptrdiff_t src2stride,
                              int height,
                              int denom,
                              int wx0,
                              int wx1,
                              int offset0,
                              int offset1,
                              intptr_t mx,
                              intptr_t my,
                              int width)
{
    put_vvc_bi_w_qpel_h16_10_avx2(dst,
                                  dststride,
                                  _src,
                                  _srcstride,
                                  src2,
                                  src2stride,
                                  height,
                                  denom,
                                  wx0,
                                  wx1,
                                  offset0,
                                  offset1,
                                  mx,
                                  my,
                                  width);
}

static void
put_vvc_uni_w_qpel_h64_10_avx2(OVSample* dst,
                               ptrdiff_t dststride,
                               const OVSample* _src,
                               ptrdiff_t _srcstride,
                               int height,
                               int denom,
                               int wx,
                               int ox,
                               intptr_t mx,
                               intptr_t my,
                               int width)
{
    put_vvc_uni_w_qpel_h16_10_avx2(dst, dststride, _src, _srcstride, height, denom, wx, ox, mx, my, width);
}

static void
put_vvc_bi_w_qpel_h64_10_avx2(OVSample* dst,
                              ptrdiff_t dststride,
                              const OVSample* _src,
                              ptrdiff_t _srcstride,
                              const int16_t* src2,
                              ptrdiff_t src2stride,
                              int height,
                              int denom,
                              int wx0,
                              int wx1,
                              int offset0,
                              int offset1,
                              intptr_t mx,
                              intptr_t my,
                              int width)
{
    put_vvc_bi_w_qpel_h16_10_avx2(dst,
                                  dststride,
                                  _src,
                                  _srcstride,
                                  src2,
                                  src2stride,
                                  height,
                                  denom,
                                  wx0,
                                  wx1,
                                  offset0,
                                  offset1,
                                  mx,
                                  my,
                                  width);
}

static void
put_vvc_uni_w_qpel_v32_10_avx2(OVSample* dst,
                               ptrdiff_t dststride,
                               const OVSample* _src,
                               ptrdiff_t _srcstride,
                               int height,
                               int denom,
                               int wx,
                               int ox,
                               intptr_t mx,
                               intptr_t my,
                               int width)
{
    put_vvc_uni_w_qpel_v16_10_avx2(dst, dststride, _src, _srcstride, height, denom, wx, ox, mx, my, width);
}

static void
put_vvc_bi_w_qpel_v32_10_avx2(OVSample* dst,
                              ptrdiff_t dststride,
                              const OVSample* _src,
                              ptrdiff_t _srcstride,
                              const int16_t* src2,
                              ptrdiff_t src2stride,
                              int height,
                              int denom,
                              int wx0,
                              int wx1,
                              int offset0,
                              int offset1,
                              intptr_t mx,
                              intptr_t my,
                              int width)
{
    put_vvc_bi_w_qpel_v16_10_avx2(dst,
                                  dststride,
                                  _src,
                                  _srcstride,
                                  src2,
                                  src2stride,
                                  height,
                                  denom,
                                  wx0,
                                  wx1,
                                  offset0,
                                  offset1,
                                  mx,
                                  my,
                                  width);
}

static void
put_vvc_uni_w_qpel_v64_10_avx2(OVSample* dst,
                               ptrdiff_t dststride,
                               const OVSample* _src,
                               ptrdiff_t _srcstride,
                               int height,
                               int denom,
                               int wx,
                               int ox,
                               intptr_t mx,
                               intptr_t my,
                               int width)
{
    put_vvc_uni_w_qpel_v16_10_avx2(dst, dststride, _src, _srcstride, height, denom, wx, ox, mx, my, width);
}

static void
put_vvc_bi_w_qpel_v64_10_avx2(OVSample* dst,
                              ptrdiff_t dststride,
                              const OVSample* _src,
                              ptrdiff_t _srcstride,
                              const int16_t* src2,
                              ptrdiff_t src2stride,
                              int height,
                              int denom,
                              int wx0,
                              int wx1,
                              int offset0,
                              int offset1,
                              intptr_t mx,
                              intptr_t my,
                              int width)
{
    put_vvc_bi_w_qpel_v16_10_avx2(dst,
                                  dststride,
                                  _src,
                                  _srcstride,
                                  src2,
                                  src2stride,
                                  height,
                                  denom,
                                  wx0,
                                  wx1,
                                  offset0,
                                  offset1,
                                  mx,
                                  my,
                                  width);
}

static void
put_vvc_uni_w_qpel_hv32_10_avx2(OVSample* dst,
                                ptrdiff_t dststride,
                                const OVSample* _src,
                                ptrdiff_t _srcstride,
                                int height,
                                int denom,
                                int wx,
                                int ox,
                                intptr_t mx,
                                intptr_t my,
                                int width)
{
    put_vvc_uni_w_qpel_hv16_10_avx2(dst, dststride, _src, _srcstride, height, denom, wx, ox, mx, my, width);
}

static void
put_vvc_bi_w_qpel_hv32_10_avx2(OVSample* dst,
                               ptrdiff_t dststride,
                               const OVSample* _src,
                               ptrdiff_t _srcstride,
                               const int16_t* src2,
                               ptrdiff_t src2stride,
                               int height,
                               int denom,
                               int wx0,
                               int wx1,
                               int offset0,
                               int offset1,
                               intptr_t mx,
                               intptr_t my,
                               int width)
{
    put_vvc_bi_w_qpel_hv16_10_avx2(dst,
                                   dststride,
                                   _src,
                                   _srcstride,
                                   src2,
                                   src2stride,
                                   height,
                                   denom,
                                   wx0,
                                   wx1,
                                   offset0,
                                   offset1,
                                   mx,
                                   my,
                                   width);
}

static void
put_vvc_uni_w_qpel_hv64_10_avx2(OVSample* dst,
                                ptrdiff_t dststride,
                                const OVSample* _src,
                                ptrdiff_t _srcstride,
                                int height,
                                int denom,
                                int wx,
                                int ox,
                                intptr_t mx,
                                intptr_t my,
                                int width)
{
    put_vvc_uni_w_qpel_hv16_10_avx2(dst, dststride, _src, _srcstride, height, denom, wx, ox, mx, my, width);
}

static void
put_vvc_bi_w_qpel_hv64_10_avx2(OVSample* dst,
                               ptrdiff_t dststride,
                               const OVSample* _src,
                               ptrdiff_t _srcstride,
                               const int16_t* src2,
                               ptrdiff_t src2stride,
                               int height,
                               int denom,
                               int wx0,
                               int wx1,
                               int offset0,
                               int offset1,
                               intptr_t mx,
                               intptr_t my,
                               int width)
{
    put_vvc_bi_w_qpel_hv16_10_avx2(dst,
                                   dststride,
                                   _src,
                                   _srcstride,
                                   src2,
                                   src2stride,
                                   height,
                                   denom,
                                   wx0,
                                   wx1,
                                   offset0,
                                   offset1,
                                   mx,
                                   my,
                                   width);
}

static void
put_vvc_uni_w_epel_h32_10_avx2(OVSample* dst,
                               ptrdiff_t dststride,
                               const OVSample* _src,
                               ptrdiff_t _srcstride,
                               int height,
                               int denom,
                               int wx,
                               int ox,
                               intptr_t mx,
                               intptr_t my,
                               int width)
{
    put_vvc_uni_w_epel_h16_10_avx2(dst, dststride, _src, _srcstride, height, denom, wx, ox, mx, my, width);
}

static void
put_vvc_bi_w_epel_h32_10_avx2(OVSample* dst,
                              ptrdiff_t dststride,
                              const OVSample* _src,
                              ptrdiff_t _srcstride,
                              const int16_t* src2,
                              ptrdiff_t src2stride,
                              int height,
                              int denom,
                              int wx0,
                              int wx1,
                              int offset0,
                              int offset1,
                              intptr_t mx,
                              intptr_t my,
                              int width)
{
    put_vvc_bi_w_epel_h16_10_avx2(dst,
                                  dststride,
                                  _src,
                                  _srcstride,
                                  src2,
                                  src2stride,
                                  height,
                                  denom,
                                  wx0,
                                  wx1,
                                  offset0,
                                  offset1,
                                  mx,
                                  my,
                                  width);
}

static void
put_vvc_uni_w_epel_h64_10_avx2(OVSample* dst,
                               ptrdiff_t dststride,
                               const OVSample* _src,
                               ptrdiff_t _srcstride,
                               int height,
                               int denom,
                               int wx,
                               int ox,
                               intptr_t mx,
                               intptr_t my,
                               int width)
{
    put_vvc_uni_w_epel_h16_10_avx2(dst, dststride, _src, _srcstride, height, denom, wx, ox, mx, my, width);
}

static void
put_vvc_bi_w_epel_h64_10_avx2(OVSample* dst,
                              ptrdiff_t dststride,
                              const OVSample* _src,
                              ptrdiff_t _srcstride,
                              const int16_t* src2,
                              ptrdiff_t src2stride,
                              int height,
                              int denom,
                              int wx0,
                              int wx1,
                              int offset0,
                              int offset1,
                              intptr_t mx,
                              intptr_t my,
                              int width)
{
    put_vvc_bi_w_epel_h16_10_avx2(dst,
                                  dststride,
                                  _src,
                                  _srcstride,
                                  src2,
                                  src2stride,
                                  height,
                                  denom,
                                  wx0,
                                  wx1,
                                  offset0,
                                  offset1,
                                  mx,
                                  my,
                                  width);
}

static void
put_vvc_uni_w_epel_v32_10_avx2(OVSample* dst,
                               ptrdiff_t dststride,
                               const OVSample* _src,
                               ptrdiff_t _srcstride,
                               int height,
                               int denom,
                               int wx,
                               int ox,
                               intptr_t mx,
                               intptr_t my,
                               int width)
{
    put_vvc_uni_w_epel_v16_10_avx2(dst, dststride, _src, _srcstride, height, denom, wx, ox, mx, my, width);
}

static void
put_vvc_bi_w_epel_v32_10_avx2(OVSample* dst,
                              ptrdiff_t dststride,
                              const OVSample* _src,
                              ptrdiff_t _srcstride,
                              const int16_t* src2,
                              ptrdiff_t src2stride,
                              int height,
                              int denom,
                              int wx0,
                              int wx1,
                              int offset0,
                              int offset1,
                              intptr_t mx,
                              intptr_t my,
                              int width)
{
    put_vvc_bi_w_epel_v16_10_avx2(dst,
                                  dststride,
                                  _src,
                                  _srcstride,
                                  src2,
                                  src2stride,
                                  height,
                                  denom,
                                  wx0,
                                  wx1,
                                  offset0,
                                  offset1,
                                  mx,
                                  my,
                                  width);
}

static void
put_vvc_uni_w_epel_v64_10_avx2(OVSample* dst,
                               ptrdiff_t dststride,
                               const OVSample* _src,
                               ptrdiff_t _srcstride,
                               int height,
                               int denom,
                               int wx,
                               int ox,
                               intptr_t mx,
                               intptr_t my,
                               int width)
{
    put_vvc_uni_w_epel_v16_10_avx2(dst, dststride, _src, _srcstride, height, denom, wx, ox, mx, my, width);
}

static void
put_vvc_bi_w_epel_v64_10_avx2(OVSample* dst,
                              ptrdiff_t dststride,
                              const OVSample* _src,
                              ptrdiff_t _srcstride,
                              const int16_t* src2,
                              ptrdiff_t src2stride,
                              int height,
                              int denom,
                              int wx0,
                              int wx1,
                              int offset0,
                              int offset1,
                              intptr_t mx,
                              intptr_t my,
                              int width)
{
    put_vvc_bi_w_epel_v16_10_avx2(dst,
                                  dststride,
                                  _src,
                                  _srcstride,
                                  src2,
                                  src2stride,
                                  height,
                                  denom,
                                  wx0,
                                  wx1,
                                  offset0,
                                  offset1,
                                  mx,
                                  my,
                                  width);
}

static void
put_vvc_uni_w_epel_hv32_10_avx2(OVSample* dst,
                                ptrdiff_t dststride,
                                const OVSample* _src,
                                ptrdiff_t _srcstride,
                                int height,
                                int denom,
                                int wx,
                                int ox,
                                intptr_t mx,
                                intptr_t my,
                                int width)
{
    put_vvc_uni_w_epel_hv16_10_avx2(dst, dststride, _src, _srcstride, height, denom, wx, ox, mx, my, width);
}

static void
put_vvc_bi_w_epel_hv32_10_avx2(OVSample* dst,
                               ptrdiff_t dststride,
                               const OVSample* _src,
                               ptrdiff_t _srcstride,
                               const int16_t* src2,
                               ptrdiff_t src2stride,
                               int height,
                               int denom,
                               int wx0,
                               int wx1,
                               int offset0,
                               int offset1,
                               intptr_t mx,
                               intptr_t my,
                               int width)
{
    put_vvc_bi_w_epel_hv16_10_avx2(dst,
                                   dststride,
                                   _src,
                                   _srcstride,
                                   src2,
                                   src2stride,
                                   height,
                                   denom,
                                   wx0,
                                   wx1,
                                   offset0,
                                   offset1,
                                   mx,
                                   my,
                                   width);
}

static void
put_vvc_uni_w_epel_hv64_10_avx2(OVSample* dst,
                                ptrdiff_t dststride,
                                const OVSample* _src,
                                ptrdiff_t _srcstride,
                                int height,
                                int denom,
                                int wx,
                                int ox,
                                intptr_t mx,
                                intptr_t my,
                                int width)
{
    put_vvc_uni_w_epel_hv16_10_avx2(dst, dststride, _src, _srcstride, height, denom, wx, ox, mx, my, width);
}

static void
put_vvc_bi_w_epel_hv64_10_avx2(OVSample* dst,
                               ptrdiff_t dststride,
                               const OVSample* _src,
                               ptrdiff_t _srcstride,
                               const int16_t* src2,
                               ptrdiff_t src2stride,
                               int height,
                               int denom,
                               int wx0,
                               int wx1,
                               int offset0,
                               int offset1,
                               intptr_t mx,
                               intptr_t my,
                               int width)
{
    put_vvc_bi_w_epel_hv16_10_avx2(dst,
                                   dststride,
                                   _src,
                                   _srcstride,
                                   src2,
                                   src2stride,
                                   height,
                                   denom,
                                   wx0,
                                   wx1,
                                   offset0,
                                   offset1,
                                   mx,
                                   my,
                                   width);
}

static void
put_vvc_qpel_bilinear_h_avx2(uint16_t* _dst, ptrdiff_t _dststride, const uint16_t* _src,
                             ptrdiff_t _srcstride, int height, intptr_t mx, intptr_t my,
                             int width)
{
    int x, y;
    const uint16_t* src = ((uint16_t*)_src);
    ptrdiff_t srcstride = _srcstride;
    uint16_t* dst = (uint16_t*)_dst;
    ptrdiff_t dststride = _dststride;
    int shift = 14 - BITDEPTH;
    __m256i c = _mm256_set1_epi32(ov_bilinear_filters_4[mx - 1]);
    __m256i offset = _mm256_set1_epi32(1 << (shift - 1));

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x+=16) {
            __m256i x1, x2, t1, t2;
            x1 = _mm256_loadu_si256((__m256i*)&src[x]);
            x2 = _mm256_loadu_si256((__m256i*)&src[x + 1]);

            t1 = _mm256_unpacklo_epi16(x1, x2);
            t2 = _mm256_unpackhi_epi16(x1, x2);

            t1 = _mm256_madd_epi16(t1, c);
            t2 = _mm256_madd_epi16(t2, c);

            x1 = _mm256_add_epi32(t1, offset);
            x2 = _mm256_add_epi32(t2, offset);

            x1 = _mm256_srai_epi32(x1, shift);
            x2 = _mm256_srai_epi32(x2, shift);

            x1 = _mm256_max_epi32(x1, _mm256_setzero_si256());
            x2 = _mm256_max_epi32(x2, _mm256_setzero_si256());

            x1 = _mm256_min_epi32(x1, _mm256_set1_epi32(1023));
            x2 = _mm256_min_epi32(x2, _mm256_set1_epi32(1023));

            x1 = _mm256_packs_epi32(x1,x2);

            _mm256_storeu_si256((__m256i*)&dst[x], x1);
        }
        src += srcstride;
        dst += dststride;
    }
}

static void
put_vvc_qpel_bilinear_v_avx2(uint16_t* _dst, ptrdiff_t _dststride, const uint16_t* _src,
                             ptrdiff_t _srcstride, int height, intptr_t mx, intptr_t my,
                             int width)
{
    {
        int x, y;
        const uint16_t* src = ((uint16_t*)_src);
        ptrdiff_t srcstride = _srcstride;
        uint16_t* dst = (uint16_t*)_dst;
        ptrdiff_t dststride = _dststride;
        __m256i c = _mm256_set1_epi32(ov_bilinear_filters_4[my - 1]);
        int shift = 14 - BITDEPTH;
        __m256i offset = _mm256_set1_epi32(1 << (shift - 1));

        for (y = 0; y < height; y++) {
            for (x = 0; x < width; x+=16) {
                __m256i x1, x2, t1, t2;
                x1 = _mm256_loadu_si256((__m256i*)&src[x]);
                x2 = _mm256_loadu_si256((__m256i*)&src[x + srcstride]);

                t1 = _mm256_unpacklo_epi16(x1, x2);
                t2 = _mm256_unpackhi_epi16(x1, x2);

                t1 = _mm256_madd_epi16(t1, c);
                t2 = _mm256_madd_epi16(t2, c);

                x1 = _mm256_add_epi32(t1, offset);
                x2 = _mm256_add_epi32(t2, offset);

                x1 = _mm256_srai_epi32(x1, shift);
                x2 = _mm256_srai_epi32(x2, shift);

                x1 = _mm256_max_epi32(x1, _mm256_setzero_si256());
                x2 = _mm256_max_epi32(x2, _mm256_setzero_si256());

                x1 = _mm256_min_epi32(x1, _mm256_set1_epi32(1023));
                x2 = _mm256_min_epi32(x2, _mm256_set1_epi32(1023));

                x1 = _mm256_packs_epi32(x1,x2);

                _mm256_storeu_si256((__m256i*)&dst[x], x1);
            }
            src += srcstride;
            dst += dststride;
        }
    }
}

static void
put_vvc_qpel_bilinear_hv_avx2(uint16_t* _dst, ptrdiff_t _dststride, const uint16_t* _src,
                              ptrdiff_t _srcstride, int height, intptr_t mx, intptr_t my,
                              int width)
{
    uint16_t tmp_array[(MAX_PB_SIZE + QPEL_EXTRA) * MAX_PB_SIZE];
    uint16_t* tmp = tmp_array;
    put_vvc_qpel_bilinear_h_avx2(tmp, MAX_PB_SIZE, _src-1, _srcstride, height+2, mx, my, width+2);
    tmp = tmp_array + 1;
    put_vvc_qpel_bilinear_v_avx2(_dst, _dststride, tmp, MAX_PB_SIZE, height, mx, my, width);
}

static void
put_weighted_ciip_pixels_avx2(uint16_t* dst, int dststride,
                              const uint16_t* src_intra, const uint16_t* src_inter, int stride_intra, int stride_inter,
                              int width, int height, int wt)
{
    int x, y;
    int shift  = 2;
    __m256i c1 = _mm256_set1_epi16(wt);
    __m256i c2 = _mm256_set1_epi16(4-wt);
    __m256i c = _mm256_unpacklo_epi16(c1, c2);
    __m256i offset = _mm256_set1_epi32(1 << (shift - 1));

    __m256i x1l, x2l, t1l, t2l;
    __m256i c1l = _mm256_set1_epi16(wt);
    __m256i c2l = _mm256_set1_epi16(4-wt);
    __m256i cl = _mm256_unpacklo_epi16(c1l, c2l);
    __m256i offsetl = _mm256_set1_epi32(1 << (shift - 1));
    for (y = 0; y < height; y++) {
        for (x = 0; x < width /*- width%16*/; x+=16) {
            __m256i x1, x2, t1, t2;
            x1l = _mm256_loadu_si256((__m256i*)&src_intra[x]);
            x2l = _mm256_loadu_si256((__m256i*)&src_inter[x]);

            t1l = _mm256_unpacklo_epi16(x1l, x2l);
            t2l = _mm256_unpackhi_epi16(x1l, x2l);

            t1l = _mm256_madd_epi16(t1l, cl);
            t2l = _mm256_madd_epi16(t2l, cl);

            x1l = _mm256_add_epi32(t1l, offsetl);
            x2l = _mm256_add_epi32(t2l, offsetl);

            x1l = _mm256_srai_epi32(x1l, shift);
            x2l = _mm256_srai_epi32(x2l, shift);

            x1l = _mm256_max_epi32(x1l, _mm256_setzero_si256());
            x2l = _mm256_max_epi32(x2l, _mm256_setzero_si256());

            x1l = _mm256_min_epi32(x1l, _mm256_set1_epi32(1023));
            x2l = _mm256_min_epi32(x2l, _mm256_set1_epi32(1023));

            x1l = _mm256_packs_epi32(x1l,x2l);

            _mm256_storeu_si256((__m256i*)&dst[x], x1l);
        }
        #if 0
        for (; x < width; x+=16) {
            __m256i x1, x2, t1, t2;
            x1 = _mm256_loadu_si256((__m256i*)&src_intra[x]);
            x2 = _mm256_loadu_si256((__m256i*)&src_inter[x]);

            t1 = _mm256_unpacklo_epi16(x1, x2);
            t2 = _mm256_unpackhi_epi16(x1, x2);

            t1 = _mm256_madd_epi16(t1, c);
            t2 = _mm256_madd_epi16(t2, c);

            x1 = _mm256_add_epi32(t1, offset);
            x2 = _mm256_add_epi32(t2, offset);

            x1 = _mm256_srai_epi32(x1, shift);
            x2 = _mm256_srai_epi32(x2, shift);

            x1 = _mm256_max_epi32(x1, _mm256_setzero_si256());
            x2 = _mm256_max_epi32(x2, _mm256_setzero_si256());

            x1 = _mm256_min_epi32(x1, _mm256_set1_epi32(1023));
            x2 = _mm256_min_epi32(x2, _mm256_set1_epi32(1023));

            x1 = _mm256_packs_epi32(x1,x2);

            _mm256_storeu_si256((__m256i*)&dst[x], x1);
        }
        #endif
        src_intra += stride_intra;
        src_inter += stride_inter;
        dst += dststride;
    }
}

void
rcn_init_mc_functions_avx2(struct RCNFunctions* const rcn_funcs)
{
    struct MCFunctions* const mc_l = &rcn_funcs->mc_l;
    struct MCFunctions* const mc_c = &rcn_funcs->mc_c;

    /* Luma functions */
    mc_l->bidir0[0][SIZE_BLOCK_16] = &oh_hevc_put_hevc_bi0_pel_pixels16_10_avx2;
    mc_l->bidir1[0][SIZE_BLOCK_16] = &oh_hevc_put_hevc_bi1_pel_pixels16_10_avx2;
    mc_l->unidir_w[0][SIZE_BLOCK_16] = &put_vvc_uni_w_pel_pixels16_10_avx2;
    mc_l->bidir_w[0][SIZE_BLOCK_16] = &put_vvc_bi_w_pel_pixels16_10_avx2;

    mc_l->unidir[1][SIZE_BLOCK_16] = &oh_hevc_put_hevc_uni_qpel_h16_10_avx2;
    mc_l->bidir0[1][SIZE_BLOCK_16] = &oh_hevc_put_hevc_bi0_qpel_h16_10_avx2;
    mc_l->bidir1[1][SIZE_BLOCK_16] = &oh_hevc_put_hevc_bi1_qpel_h16_10_avx2;
    mc_l->bilinear[1][SIZE_BLOCK_16] = &put_vvc_qpel_bilinear_h_avx2;
    mc_l->unidir_w[1][SIZE_BLOCK_16] = &put_vvc_uni_w_qpel_h16_10_avx2;
    mc_l->bidir_w[1][SIZE_BLOCK_16] = &put_vvc_bi_w_qpel_h16_10_avx2;

    mc_l->unidir[2][SIZE_BLOCK_16] = &oh_hevc_put_hevc_uni_qpel_v16_10_avx2;
    mc_l->bidir0[2][SIZE_BLOCK_16] = &oh_hevc_put_hevc_bi0_qpel_v16_10_avx2;
    mc_l->bidir1[2][SIZE_BLOCK_16] = &oh_hevc_put_hevc_bi1_qpel_v16_10_avx2;
    mc_l->bilinear[2][SIZE_BLOCK_16] = &put_vvc_qpel_bilinear_v_avx2;
    mc_l->unidir_w[2][SIZE_BLOCK_16] = &put_vvc_uni_w_qpel_v16_10_avx2;
    mc_l->bidir_w[2][SIZE_BLOCK_16] = &put_vvc_bi_w_qpel_v16_10_avx2;

    mc_l->unidir[3][SIZE_BLOCK_16] = &oh_hevc_put_hevc_uni_qpel_hv16_10_avx2;
    mc_l->bidir0[3][SIZE_BLOCK_16] = &oh_hevc_put_hevc_bi0_qpel_hv16_10_avx2;
    mc_l->bidir1[3][SIZE_BLOCK_16] = &oh_hevc_put_hevc_bi1_qpel_hv16_10_avx2;
    mc_l->bilinear[3][SIZE_BLOCK_16] = &put_vvc_qpel_bilinear_hv_avx2;
    mc_l->unidir_w[3][SIZE_BLOCK_16] = &put_vvc_uni_w_qpel_hv16_10_avx2;
    mc_l->bidir_w[3][SIZE_BLOCK_16] = &put_vvc_bi_w_qpel_hv16_10_avx2;

    mc_l->bidir0[0][SIZE_BLOCK_32] = &oh_hevc_put_hevc_bi0_pel_pixels32_10_avx2;
    mc_l->bidir1[0][SIZE_BLOCK_32] = &oh_hevc_put_hevc_bi1_pel_pixels32_10_avx2;
    mc_l->unidir_w[0][SIZE_BLOCK_32] = &put_vvc_uni_w_pel_pixels32_10_avx2;
    mc_l->bidir_w[0][SIZE_BLOCK_32] = &put_vvc_bi_w_pel_pixels32_10_avx2;

    mc_l->unidir[1][SIZE_BLOCK_32] = &oh_hevc_put_hevc_uni_qpel_h32_10_avx2;
    mc_l->bidir0[1][SIZE_BLOCK_32] = &oh_hevc_put_hevc_bi0_qpel_h32_10_avx2;
    mc_l->bidir1[1][SIZE_BLOCK_32] = &oh_hevc_put_hevc_bi1_qpel_h32_10_avx2;
    mc_l->bilinear[1][SIZE_BLOCK_32] = &put_vvc_qpel_bilinear_h_avx2;
    mc_l->unidir_w[1][SIZE_BLOCK_32] = &put_vvc_uni_w_qpel_h32_10_avx2;
    mc_l->bidir_w[1][SIZE_BLOCK_32] = &put_vvc_bi_w_qpel_h32_10_avx2;

    mc_l->unidir[2][SIZE_BLOCK_32] = &oh_hevc_put_hevc_uni_qpel_v32_10_avx2;
    mc_l->bidir0[2][SIZE_BLOCK_32] = &oh_hevc_put_hevc_bi0_qpel_v32_10_avx2;
    mc_l->bidir1[2][SIZE_BLOCK_32] = &oh_hevc_put_hevc_bi1_qpel_v32_10_avx2;
    mc_l->bilinear[2][SIZE_BLOCK_32] = &put_vvc_qpel_bilinear_v_avx2;
    mc_l->unidir_w[2][SIZE_BLOCK_32] = &put_vvc_uni_w_qpel_v32_10_avx2;
    mc_l->bidir_w[2][SIZE_BLOCK_32] = &put_vvc_bi_w_qpel_v32_10_avx2;

    mc_l->unidir[3][SIZE_BLOCK_32] = &oh_hevc_put_hevc_uni_qpel_hv32_10_avx2;
    mc_l->bidir0[3][SIZE_BLOCK_32] = &oh_hevc_put_hevc_bi0_qpel_hv32_10_avx2;
    mc_l->bidir1[3][SIZE_BLOCK_32] = &oh_hevc_put_hevc_bi1_qpel_hv32_10_avx2;
    mc_l->bilinear[3][SIZE_BLOCK_32] = &put_vvc_qpel_bilinear_hv_avx2;
    mc_l->unidir_w[3][SIZE_BLOCK_32] = &put_vvc_uni_w_qpel_hv32_10_avx2;
    mc_l->bidir_w[3][SIZE_BLOCK_32] = &put_vvc_bi_w_qpel_hv32_10_avx2;

    mc_l->bidir0[0][SIZE_BLOCK_64] = &oh_hevc_put_hevc_bi0_pel_pixels64_10_avx2;
    mc_l->bidir1[0][SIZE_BLOCK_64] = &oh_hevc_put_hevc_bi1_pel_pixels64_10_avx2;
    mc_l->unidir_w[0][SIZE_BLOCK_64] = &put_vvc_uni_w_pel_pixels64_10_avx2;
    mc_l->bidir_w[0][SIZE_BLOCK_64] = &put_vvc_bi_w_pel_pixels64_10_avx2;

    mc_l->unidir[1][SIZE_BLOCK_64] = &oh_hevc_put_hevc_uni_qpel_h64_10_avx2;
    mc_l->bidir0[1][SIZE_BLOCK_64] = &oh_hevc_put_hevc_bi0_qpel_h64_10_avx2;
    mc_l->bidir1[1][SIZE_BLOCK_64] = &oh_hevc_put_hevc_bi1_qpel_h64_10_avx2;
    mc_l->bilinear[1][SIZE_BLOCK_64] = &put_vvc_qpel_bilinear_h_avx2;
    mc_l->unidir_w[1][SIZE_BLOCK_64] = &put_vvc_uni_w_qpel_h64_10_avx2;
    mc_l->bidir_w[1][SIZE_BLOCK_64] = &put_vvc_bi_w_qpel_h64_10_avx2;

    mc_l->unidir[2][SIZE_BLOCK_64] = &oh_hevc_put_hevc_uni_qpel_v64_10_avx2;
    mc_l->bidir0[2][SIZE_BLOCK_64] = &oh_hevc_put_hevc_bi0_qpel_v64_10_avx2;
    mc_l->bidir1[2][SIZE_BLOCK_64] = &oh_hevc_put_hevc_bi1_qpel_v64_10_avx2;
    mc_l->bilinear[2][SIZE_BLOCK_64] = &put_vvc_qpel_bilinear_v_avx2;
    mc_l->unidir_w[2][SIZE_BLOCK_64] = &put_vvc_uni_w_qpel_v64_10_avx2;
    mc_l->bidir_w[2][SIZE_BLOCK_64] = &put_vvc_bi_w_qpel_v64_10_avx2;

    mc_l->unidir[3][SIZE_BLOCK_64] = &oh_hevc_put_hevc_uni_qpel_hv64_10_avx2;
    mc_l->bidir0[3][SIZE_BLOCK_64] = &oh_hevc_put_hevc_bi0_qpel_hv64_10_avx2;
    mc_l->bidir1[3][SIZE_BLOCK_64] = &oh_hevc_put_hevc_bi1_qpel_hv64_10_avx2;
    mc_l->bilinear[3][SIZE_BLOCK_64] = &put_vvc_qpel_bilinear_hv_avx2;
    mc_l->unidir_w[3][SIZE_BLOCK_64] = &put_vvc_uni_w_qpel_hv64_10_avx2;
    mc_l->bidir_w[3][SIZE_BLOCK_64] = &put_vvc_bi_w_qpel_hv64_10_avx2;

    mc_l->bidir0[0][SIZE_BLOCK_128] = &oh_hevc_put_hevc_bi0_pel_pixels64_10_avx2;
    mc_l->bidir1[0][SIZE_BLOCK_128] = &oh_hevc_put_hevc_bi1_pel_pixels64_10_avx2;
    mc_l->unidir_w[0][SIZE_BLOCK_128] = &put_vvc_uni_w_pel_pixels64_10_avx2;
    mc_l->bidir_w[0][SIZE_BLOCK_128] = &put_vvc_bi_w_pel_pixels64_10_avx2;

    mc_l->unidir[1][SIZE_BLOCK_128] = &oh_hevc_put_hevc_uni_qpel_h64_10_avx2;
    mc_l->bidir0[1][SIZE_BLOCK_128] = &oh_hevc_put_hevc_bi0_qpel_h64_10_avx2;
    mc_l->bidir1[1][SIZE_BLOCK_128] = &oh_hevc_put_hevc_bi1_qpel_h64_10_avx2;
    mc_l->bilinear[1][SIZE_BLOCK_128] = &put_vvc_qpel_bilinear_h_avx2;
    mc_l->unidir_w[1][SIZE_BLOCK_128] = &put_vvc_uni_w_qpel_h64_10_avx2;
    mc_l->bidir_w[1][SIZE_BLOCK_128] = &put_vvc_bi_w_qpel_h64_10_avx2;

    mc_l->unidir[2][SIZE_BLOCK_128] = &oh_hevc_put_hevc_uni_qpel_v64_10_avx2;
    mc_l->bidir0[2][SIZE_BLOCK_128] = &oh_hevc_put_hevc_bi0_qpel_v64_10_avx2;
    mc_l->bidir1[2][SIZE_BLOCK_128] = &oh_hevc_put_hevc_bi1_qpel_v64_10_avx2;
    mc_l->bilinear[2][SIZE_BLOCK_128] = &put_vvc_qpel_bilinear_v_avx2;
    mc_l->unidir_w[2][SIZE_BLOCK_128] = &put_vvc_uni_w_qpel_v64_10_avx2;
    mc_l->bidir_w[2][SIZE_BLOCK_128] = &put_vvc_bi_w_qpel_v64_10_avx2;

    mc_l->unidir[3][SIZE_BLOCK_128] = &oh_hevc_put_hevc_uni_qpel_hv64_10_avx2;
    mc_l->bidir0[3][SIZE_BLOCK_128] = &oh_hevc_put_hevc_bi0_qpel_hv64_10_avx2;
    mc_l->bidir1[3][SIZE_BLOCK_128] = &oh_hevc_put_hevc_bi1_qpel_hv64_10_avx2;
    mc_l->bilinear[3][SIZE_BLOCK_128] = &put_vvc_qpel_bilinear_hv_avx2;
    mc_l->unidir_w[3][SIZE_BLOCK_128] = &put_vvc_uni_w_qpel_hv64_10_avx2;
    mc_l->bidir_w[3][SIZE_BLOCK_128] = &put_vvc_bi_w_qpel_hv64_10_avx2;

    /* Chroma functions */
    mc_c->bidir0[0][SIZE_BLOCK_16] = &oh_hevc_put_hevc_bi0_pel_pixels16_10_avx2;
    mc_c->bidir1[0][SIZE_BLOCK_16] = &oh_hevc_put_hevc_bi1_pel_pixels16_10_avx2;
    mc_c->unidir_w[0][SIZE_BLOCK_16] = &put_vvc_uni_w_pel_pixels16_10_avx2;
    mc_c->bidir_w[0][SIZE_BLOCK_16] = &put_vvc_bi_w_pel_pixels16_10_avx2;

    mc_c->unidir[1][SIZE_BLOCK_16] = &oh_hevc_put_hevc_uni_epel_h16_10_avx2;
    mc_c->bidir0[1][SIZE_BLOCK_16] = &oh_hevc_put_hevc_bi0_epel_h16_10_avx2;
    mc_c->bidir1[1][SIZE_BLOCK_16] = &oh_hevc_put_hevc_bi1_epel_h16_10_avx2;
    mc_c->unidir_w[1][SIZE_BLOCK_16] = &put_vvc_uni_w_epel_h16_10_avx2;
    mc_c->bidir_w[1][SIZE_BLOCK_16] = &put_vvc_bi_w_epel_h16_10_avx2;

    mc_c->unidir[2][SIZE_BLOCK_16] = &oh_hevc_put_hevc_uni_epel_v16_10_avx2;
    mc_c->bidir0[2][SIZE_BLOCK_16] = &oh_hevc_put_hevc_bi0_epel_v16_10_avx2;
    mc_c->bidir1[2][SIZE_BLOCK_16] = &oh_hevc_put_hevc_bi1_epel_v16_10_avx2;
    mc_c->unidir_w[2][SIZE_BLOCK_16] = &put_vvc_uni_w_epel_v16_10_avx2;
    mc_c->bidir_w[2][SIZE_BLOCK_16] = &put_vvc_bi_w_epel_v16_10_avx2;

    mc_c->unidir[3][SIZE_BLOCK_16] = &oh_hevc_put_hevc_uni_epel_hv16_10_avx2;
    mc_c->bidir0[3][SIZE_BLOCK_16] = &oh_hevc_put_hevc_bi0_epel_hv16_10_avx2;
    mc_c->bidir1[3][SIZE_BLOCK_16] = &oh_hevc_put_hevc_bi1_epel_hv16_10_avx2;
    mc_c->unidir_w[3][SIZE_BLOCK_16] = &put_vvc_uni_w_epel_hv16_10_avx2;
    mc_c->bidir_w[3][SIZE_BLOCK_16] = &put_vvc_bi_w_epel_hv16_10_avx2;

    mc_c->bidir0[0][SIZE_BLOCK_32] = &oh_hevc_put_hevc_bi0_pel_pixels32_10_avx2;
    mc_c->bidir1[0][SIZE_BLOCK_32] = &oh_hevc_put_hevc_bi1_pel_pixels32_10_avx2;
    mc_c->unidir_w[0][SIZE_BLOCK_32] = &put_vvc_uni_w_pel_pixels32_10_avx2;
    mc_c->bidir_w[0][SIZE_BLOCK_32] = &put_vvc_bi_w_pel_pixels32_10_avx2;

    mc_c->unidir[1][SIZE_BLOCK_32] = &oh_hevc_put_hevc_uni_epel_h32_10_avx2;
    mc_c->bidir0[1][SIZE_BLOCK_32] = &oh_hevc_put_hevc_bi0_epel_h32_10_avx2;
    mc_c->bidir1[1][SIZE_BLOCK_32] = &oh_hevc_put_hevc_bi1_epel_h32_10_avx2;
    mc_c->unidir_w[1][SIZE_BLOCK_32] = &put_vvc_uni_w_epel_h32_10_avx2;
    mc_c->bidir_w[1][SIZE_BLOCK_32] = &put_vvc_bi_w_epel_h32_10_avx2;

    mc_c->unidir[2][SIZE_BLOCK_32] = &oh_hevc_put_hevc_uni_epel_v32_10_avx2;
    mc_c->bidir0[2][SIZE_BLOCK_32] = &oh_hevc_put_hevc_bi0_epel_v32_10_avx2;
    mc_c->bidir1[2][SIZE_BLOCK_32] = &oh_hevc_put_hevc_bi1_epel_v32_10_avx2;
    mc_c->unidir_w[2][SIZE_BLOCK_32] = &put_vvc_uni_w_epel_v32_10_avx2;
    mc_c->bidir_w[2][SIZE_BLOCK_32] = &put_vvc_bi_w_epel_v32_10_avx2;

    mc_c->unidir[3][SIZE_BLOCK_32] = &oh_hevc_put_hevc_uni_epel_hv32_10_avx2;
    mc_c->bidir0[3][SIZE_BLOCK_32] = &oh_hevc_put_hevc_bi0_epel_hv32_10_avx2;
    mc_c->bidir1[3][SIZE_BLOCK_32] = &oh_hevc_put_hevc_bi1_epel_hv32_10_avx2;
    mc_c->unidir_w[3][SIZE_BLOCK_32] = &put_vvc_uni_w_epel_hv32_10_avx2;
    mc_c->bidir_w[3][SIZE_BLOCK_32] = &put_vvc_bi_w_epel_hv32_10_avx2;

    mc_c->bidir0[0][SIZE_BLOCK_64] = &oh_hevc_put_hevc_bi0_pel_pixels64_10_avx2;
    mc_c->bidir1[0][SIZE_BLOCK_64] = &oh_hevc_put_hevc_bi1_pel_pixels64_10_avx2;
    mc_c->unidir_w[0][SIZE_BLOCK_64] = &put_vvc_uni_w_pel_pixels64_10_avx2;
    mc_c->bidir_w[0][SIZE_BLOCK_64] = &put_vvc_bi_w_pel_pixels64_10_avx2;

    mc_c->unidir[1][SIZE_BLOCK_64] = &oh_hevc_put_hevc_uni_epel_h64_10_avx2;
    mc_c->bidir0[1][SIZE_BLOCK_64] = &oh_hevc_put_hevc_bi0_epel_h64_10_avx2;
    mc_c->bidir1[1][SIZE_BLOCK_64] = &oh_hevc_put_hevc_bi1_epel_h64_10_avx2;
    mc_c->unidir_w[1][SIZE_BLOCK_64] = &put_vvc_uni_w_epel_h64_10_avx2;
    mc_c->bidir_w[1][SIZE_BLOCK_64] = &put_vvc_bi_w_epel_h64_10_avx2;

    mc_c->unidir[2][SIZE_BLOCK_64] = &oh_hevc_put_hevc_uni_epel_v64_10_avx2;
    mc_c->bidir0[2][SIZE_BLOCK_64] = &oh_hevc_put_hevc_bi0_epel_v64_10_avx2;
    mc_c->bidir1[2][SIZE_BLOCK_64] = &oh_hevc_put_hevc_bi1_epel_v64_10_avx2;
    mc_c->unidir_w[2][SIZE_BLOCK_64] = &put_vvc_uni_w_epel_v64_10_avx2;
    mc_c->bidir_w[2][SIZE_BLOCK_64] = &put_vvc_bi_w_epel_v64_10_avx2;

    mc_c->unidir[3][SIZE_BLOCK_64] = &oh_hevc_put_hevc_uni_epel_hv64_10_avx2;
    mc_c->bidir0[3][SIZE_BLOCK_64] = &oh_hevc_put_hevc_bi0_epel_hv64_10_avx2;
    mc_c->bidir1[3][SIZE_BLOCK_64] = &oh_hevc_put_hevc_bi1_epel_hv64_10_avx2;
    mc_c->unidir_w[3][SIZE_BLOCK_64] = &put_vvc_uni_w_epel_hv64_10_avx2;
    mc_c->bidir_w[3][SIZE_BLOCK_64] = &put_vvc_bi_w_epel_hv64_10_avx2;

    mc_c->bidir0[0][SIZE_BLOCK_128] = &oh_hevc_put_hevc_bi0_pel_pixels64_10_avx2;
    mc_c->bidir1[0][SIZE_BLOCK_128] = &oh_hevc_put_hevc_bi1_pel_pixels64_10_avx2;
    mc_c->unidir_w[0][SIZE_BLOCK_128] = &put_vvc_uni_w_pel_pixels64_10_avx2;
    mc_c->bidir_w[0][SIZE_BLOCK_128] = &put_vvc_bi_w_pel_pixels64_10_avx2;

    mc_c->unidir[1][SIZE_BLOCK_128] = &oh_hevc_put_hevc_uni_epel_h64_10_avx2;
    mc_c->bidir0[1][SIZE_BLOCK_128] = &oh_hevc_put_hevc_bi0_epel_h64_10_avx2;
    mc_c->bidir1[1][SIZE_BLOCK_128] = &oh_hevc_put_hevc_bi1_epel_h64_10_avx2;
    mc_c->unidir_w[1][SIZE_BLOCK_128] = &put_vvc_uni_w_epel_h64_10_avx2;
    mc_c->bidir_w[1][SIZE_BLOCK_128] = &put_vvc_bi_w_epel_h64_10_avx2;

    mc_c->unidir[2][SIZE_BLOCK_128] = &oh_hevc_put_hevc_uni_epel_v64_10_avx2;
    mc_c->bidir0[2][SIZE_BLOCK_128] = &oh_hevc_put_hevc_bi0_epel_v64_10_avx2;
    mc_c->bidir1[2][SIZE_BLOCK_128] = &oh_hevc_put_hevc_bi1_epel_v64_10_avx2;
    mc_c->unidir_w[2][SIZE_BLOCK_128] = &put_vvc_uni_w_epel_v64_10_avx2;
    mc_c->bidir_w[2][SIZE_BLOCK_128] = &put_vvc_bi_w_epel_v64_10_avx2;

    mc_c->unidir[3][SIZE_BLOCK_128] = &oh_hevc_put_hevc_uni_epel_hv64_10_avx2;
    mc_c->bidir0[3][SIZE_BLOCK_128] = &oh_hevc_put_hevc_bi0_epel_hv64_10_avx2;
    mc_c->bidir1[3][SIZE_BLOCK_128] = &oh_hevc_put_hevc_bi1_epel_hv64_10_avx2;
    mc_c->unidir_w[3][SIZE_BLOCK_128] = &put_vvc_uni_w_epel_hv64_10_avx2;
    mc_c->bidir_w[3][SIZE_BLOCK_128] = &put_vvc_bi_w_epel_hv64_10_avx2;
}

void
rcn_init_ciip_functions_avx2(struct RCNFunctions *const rcn_funcs)
{
    struct CIIPFunctions *const ciip = &rcn_funcs->ciip;
    ciip->weighted = &put_weighted_ciip_pixels_avx2;
}
