#include "ovmem.h"
#include "rcn_structures.h"
#include <emmintrin.h>
#include <stdint.h>
#include <string.h>
#include <tmmintrin.h>

#define SIZE_BLOCK_4 2
#define SIZE_BLOCK_8 3
#define SIZE_BLOCK_16 4
#define SIZE_BLOCK_32 5
#define SIZE_BLOCK_64 6

DECLARE_ALIGNED(16, const int16_t, oh_hevc_epel_filters_sse[31][2][8]) = {
  //{  0, 64,  0,  0 },
  { { -1, 63, -1, 63, -1, 63, -1, 63 }, { 2, 0, 2, 0, 2, 0, 2, 0 } },
  { { -2, 62, -2, 62, -2, 62, -2, 62 }, { 4, 0, 4, 0, 4, 0, 4, 0 } },
  { { -2, 60, -2, 60, -2, 60, -2, 60 }, { 7, -1, 7, -1, 7, -1, 7, -1 } },
  { { -2, 58, -2, 58, -2, 58, -2, 58 }, { 10, -2, 10, -2, 10, -2, 10, -2 } },
  { { -3, 57, -3, 57, -3, 57, -3, 57 }, { 12, -2, 12, -2, 12, -2, 12, -2 } },
  { { -4, 56, -4, 56, -4, 56, -4, 56 }, { 14, -2, 14, -2, 14, -2, 14, -2 } },
  { { -4, 55, -4, 55, -4, 55, -4, 55 }, { 15, -2, 15, -2, 15, -2, 15, -2 } },
  { { -4, 54, -4, 54, -4, 54, -4, 54 }, { 16, -2, 16, -2, 16, -2, 16, -2 } },
  { { -5, 53, -5, 53, -5, 53, -5, 53 }, { 18, -2, 18, -2, 18, -2, 18, -2 } },
  { { -6, 52, -6, 52, -6, 52, -6, 52 }, { 20, -2, 20, -2, 20, -2, 20, -2 } },
  { { -6, 49, -6, 49, -6, 49, -6, 49 }, { 24, -3, 24, -3, 24, -3, 24, -3 } },
  { { -6, 46, -6, 46, -6, 46, -6, 46 }, { 28, -4, 28, -4, 28, -4, 28, -4 } },
  { { -5, 44, -5, 44, -5, 44, -5, 44 }, { 29, -4, 29, -4, 29, -4, 29, -4 } },
  { { -4, 42, -4, 42, -4, 42, -4, 42 }, { 30, -4, 30, -4, 30, -4, 30, -4 } },
  { { -4, 39, -4, 39, -4, 39, -4, 39 }, { 33, -4, 33, -4, 33, -4, 33, -4 } },
  { { -4, 36, -4, 36, -4, 36, -4, 36 }, { 36, -4, 36, -4, 36, -4, 36, -4 } },
  { { -4, 33, -4, 33, -4, 33, -4, 33 }, { 39, -4, 39, -4, 39, -4, 39, -4 } },
  { { -4, 30, -4, 30, -4, 30, -4, 30 }, { 42, -4, 42, -4, 42, -4, 42, -4 } },
  { { -4, 29, -4, 29, -4, 29, -4, 29 }, { 44, -5, 44, -5, 44, -5, 44, -5 } },
  { { -4, 28, -4, 28, -4, 28, -4, 28 }, { 46, -6, 46, -6, 46, -6, 46, -6 } },
  { { -3, 24, -3, 24, -3, 24, -3, 24 }, { 49, -6, 49, -6, 49, -6, 49, -6 } },
  { { -2, 20, -2, 20, -2, 20, -2, 20 }, { 52, -6, 52, -6, 52, -6, 52, -6 } },
  { { -2, 18, -2, 18, -2, 18, -2, 18 }, { 53, -5, 53, -5, 53, -5, 53, -5 } },
  { { -2, 16, -2, 16, -2, 16, -2, 16 }, { 54, -4, 54, -4, 54, -4, 54, -4 } },
  { { -2, 15, -2, 15, -2, 15, -2, 15 }, { 55, -4, 55, -4, 55, -4, 55, -4 } },
  { { -2, 14, -2, 14, -2, 14, -2, 14 }, { 56, -4, 56, -4, 56, -4, 56, -4 } },
  { { -2, 12, -2, 12, -2, 12, -2, 12 }, { 57, -3, 57, -3, 57, -3, 57, -3 } },
  { { -2, 10, -2, 10, -2, 10, -2, 10 }, { 58, -2, 58, -2, 58, -2, 58, -2 } },
  { { -1, 7, -1, 7, -1, 7, -1, 7 }, { 60, -2, 60, -2, 60, -2, 60, -2 } },
  { { 0, 4, 0, 4, 0, 4, 0, 4 }, { 62, -2, 62, -2, 62, -2, 62, -2 } },
  { { 0, 2, 0, 2, 0, 2, 0, 2 }, { 63, -1, 63, -1, 63, -1, 63, -1 } },
};

DECLARE_ALIGNED(16, const int16_t, oh_hevc_qpel_filters_sse[15][4][8]) = {
#if 0
  { { 0, 0, 0, 0, 0, 0, 0, 0},
  {0, 64, 0, 64, 0, 64, 0, 64},
  {0,   0, 0,   0, 0,   0, 0,   0},
  {0,  0 , 0,  0 , 0,  0 , 0,  0 }},
#endif
  { { 0, 1, 0, 1, 0, 1, 0, 1 },
    { -3, 63, -3, 63, -3, 63, -3, 63 },
    { 4, -2, 4, -2, 4, -2, 4, -2 },
    { 1, 0, 1, 0, 1, 0, 1, 0 } },
  { { -1, 2, -1, 2, -1, 2, -1, 2 },
    { -5, 62, -5, 62, -5, 62, -5, 62 },
    { 8, -3, 8, -3, 8, -3, 8, -3 },
    { 1, 0, 1, 0, 1, 0, 1, 0 } },
  { { -1, 3, -1, 3, -1, 3, -1, 3 },
    { -8, 60, -8, 60, -8, 60, -8, 60 },
    { 13, -4, 13, -4, 13, -4, 13, -4 },
    { 1, 0, 1, 0, 1, 0, 1, 0 } },
  { { -1, 4, -1, 4, -1, 4, -1, 4 },
    { -10, 58, -10, 58, -10, 58, -10, 58 },
    { 17, -5, 17, -5, 17, -5, 17, -5 },
    { 1, 0, 1, 0, 1, 0, 1, 0 } },
  { { -1, 4, -1, 4, -1, 4, -1, 4 },
    { -11, 52, -11, 52, -11, 52, -11, 52 },
    { 26, -8, 26, -8, 26, -8, 26, -8 },
    { 3, -1, 3, -1, 3, -1, 3, -1 } },
  { { -1, 3, -1, 3, -1, 3, -1, 3 },
    { -9, 47, -9, 47, -9, 47, -9, 47 },
    { 31, -10, 31, -10, 31, -10, 31, -10 },
    { 4, -1, 4, -1, 4, -1, 4, -1 } },
  { { -1, 4, -1, 4, -1, 4, -1, 4 },
    { -11, 45, -11, 45, -11, 45, -11, 45 },
    { 34, -10, 34, -10, 34, -10, 34, -10 },
    { 4, -1, 4, -1, 4, -1, 4, -1 } },
  { { -1, 4, -1, 4, -1, 4, -1, 4 },
    { -11, 40, -11, 40, -11, 40, -11, 40 },
    { 40, -11, 40, -11, 40, -11, 40, -11 },
    { 4, -1, 4, -1, 4, -1, 4, -1 } },
  { { -1, 4, -1, 4, -1, 4, -1, 4 },
    { -10, 34, -10, 34, -10, 34, -10, 34 },
    { 45, -11, 45, -11, 45, -11, 45, -11 },
    { 4, -1, 4, -1, 4, -1, 4, -1 } },
  { { -1, 4, -1, 4, -1, 4, -1, 4 },
    { -10, 31, -10, 31, -10, 31, -10, 31 },
    { 47, -9, 47, -9, 47, -9, 47, -9 },
    { 3, -1, 3, -1, 3, -1, 3, -1 } },
  { { -1, 3, -1, 3, -1, 3, -1, 3 },
    { -8, 26, -8, 26, -8, 26, -8, 26 },
    { 52, -11, 52, -11, 52, -11, 52, -11 },
    { 4, -1, 4, -1, 4, -1, 4, -1 } },
  { { 0, 1, 0, 1, 0, 1, 0, 1 },
    { -5, 17, -5, 17, -5, 17, -5, 17 },
    { 58, -10, 58, -10, 58, -10, 58, -10 },
    { 4, -1, 4, -1, 4, -1, 4, -1 } },
  { { 0, 1, 0, 1, 0, 1, 0, 1 },
    { -4, 13, -4, 13, -4, 13, -4, 13 },
    { 60, -8, 60, -8, 60, -8, 60, -8 },
    { 3, -1, 3, -1, 3, -1, 3, -1 } },
  { { 0, 1, 0, 1, 0, 1, 0, 1 },
    { -3, 8, -3, 8, -3, 8, -3, 8 },
    { 62, -5, 62, -5, 62, -5, 62, -5 },
    { 2, -1, 2, -1, 2, -1, 2, -1 } },
  { { 0, 1, 0, 1, 0, 1, 0, 1 },
    { -2, 4, -2, 4, -2, 4, -2, 4 },
    { 63, -3, 63, -3, 63, -3, 63, -3 },
    { 1, 0, 1, 0, 1, 0, 1, 0 } }
};

static void
oh_hevc_put_hevc_bi0_pel_pixels4_10_sse(int16_t* dst,
                                        const uint16_t* _src,
                                        ptrdiff_t _srcstride,
                                        int height,
                                        intptr_t mx,
                                        intptr_t my,
                                        int width)
{
  int x, y;
  __m128i x1;
  const uint16_t* src = _src;
  const int srcstride = _srcstride;
  for (y = 0; y < height; y++) {
    for (x = 0; x < width; x += 4) {
      x1 = _mm_loadu_si128((__m128i*)&src[x]);
      x1 = _mm_slli_epi16(x1, 14 - 10);
      _mm_storel_epi64((__m128i*)&dst[x], x1);
    }
    src += srcstride;
    dst += 128;
  }
}

static void
oh_hevc_put_hevc_bi1_pel_pixels4_10_sse(uint16_t* _dst,
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
  __m128i x1;
  const uint16_t* src = _src;
  const int srcstride = _srcstride;
  const __m128i offset = _mm_set1_epi16(1 << 10);
  uint16_t* dst = (uint16_t*)_dst;
  const int dststride = _dststride;
  for (y = 0; y < height; y++) {
    for (x = 0; x < width; x += 4) {
      __m128i r5;
      x1 = _mm_loadu_si128((__m128i*)&src[x]);
      x1 = _mm_slli_epi16(x1, 14 - 10);
      r5 = _mm_loadl_epi64((__m128i*)&src2[x]);
      x1 = _mm_adds_epi16(x1, r5);
      x1 = _mm_mulhrs_epi16(x1, offset);
      x1 = _mm_max_epi16(x1, _mm_setzero_si128());
      x1 = _mm_min_epi16(x1, _mm_set1_epi16(0x03FF));
      _mm_storel_epi64((__m128i*)&dst[x], x1);
    }
    src += srcstride;
    src2 += 128;
    dst += dststride;
  }
}

static void
oh_hevc_put_hevc_uni_pel_pixels4_10_sse(uint16_t* _dst,
                                        ptrdiff_t _dststride,
                                        const uint16_t* _src,
                                        ptrdiff_t _srcstride,
                                        int height,
                                        intptr_t mx,
                                        intptr_t my,
                                        int width)
{
  int y;
  const uint16_t* src = _src;
  const int srcstride = _srcstride;
  uint16_t* dst = (uint16_t*)_dst;
  const int dststride = _dststride;
  for (y = 0; y < height; y++) {
    memcpy(dst, src, width * ((10 + 7) / 8));
    src += srcstride;
    dst += dststride;
  }
}

static void
oh_hevc_put_hevc_bi0_pel_pixels8_10_sse(int16_t* dst,
                                        const uint16_t* _src,
                                        ptrdiff_t _srcstride,
                                        int height,
                                        intptr_t mx,
                                        intptr_t my,
                                        int width)
{
  int x, y;
  __m128i x1;
  const uint16_t* src = _src;
  const int srcstride = _srcstride;
  for (y = 0; y < height; y++) {
    for (x = 0; x < width; x += 8) {
      x1 = _mm_loadu_si128((__m128i*)&src[x]);
      x1 = _mm_slli_epi16(x1, 14 - 10);
      _mm_storeu_si128((__m128i*)&dst[x], x1);
    }
    src += srcstride;
    dst += 128;
  }
}

static void
oh_hevc_put_hevc_bi1_pel_pixels8_10_sse(uint16_t* _dst,
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
  __m128i x1;
  const uint16_t* src = _src;
  const int srcstride = _srcstride;
  const __m128i offset = _mm_set1_epi16(1 << 10);
  uint16_t* dst = (uint16_t*)_dst;
  const int dststride = _dststride;
  for (y = 0; y < height; y++) {
    for (x = 0; x < width; x += 8) {
      __m128i r5;
      x1 = _mm_loadu_si128((__m128i*)&src[x]);
      x1 = _mm_slli_epi16(x1, 14 - 10);
      r5 = _mm_load_si128((__m128i*)&src2[x]);
      x1 = _mm_adds_epi16(x1, r5);
      x1 = _mm_mulhrs_epi16(x1, offset);
      x1 = _mm_max_epi16(x1, _mm_setzero_si128());
      x1 = _mm_min_epi16(x1, _mm_set1_epi16(0x03FF));
      _mm_storeu_si128((__m128i*)&dst[x], x1);
    }
    src += srcstride;
    src2 += 128;
    dst += dststride;
  }
}

static void
oh_hevc_put_hevc_uni_pel_pixels8_10_sse(uint16_t* _dst,
                                        ptrdiff_t _dststride,
                                        const uint16_t* _src,
                                        ptrdiff_t _srcstride,
                                        int height,
                                        intptr_t mx,
                                        intptr_t my,
                                        int width)
{
  int y;
  const uint16_t* src = _src;
  const int srcstride = _srcstride;
  uint16_t* dst = (uint16_t*)_dst;
  const int dststride = _dststride;
  for (y = 0; y < height; y++) {
    memcpy(dst, src, width * ((10 + 7) / 8));
    src += srcstride;
    dst += dststride;
  }
}

static void
oh_hevc_put_hevc_bi0_epel_h4_10_sse(int16_t* dst,
                                    const uint16_t* _src,
                                    ptrdiff_t _srcstride,
                                    int height,
                                    intptr_t mx,
                                    intptr_t my,
                                    int width)
{
  int x, y;
  __m128i x1, x2, x3, x4, t1, t2, f1, f2;
  uint16_t* src = ((uint16_t*)_src) - 1;
  ptrdiff_t srcstride = _srcstride;
  f1 = _mm_load_si128((__m128i*)oh_hevc_epel_filters_sse[mx - 1][0]);
  f2 = _mm_load_si128((__m128i*)oh_hevc_epel_filters_sse[mx - 1][1]);
  for (y = 0; y < height; y++) {
    for (x = 0; x < width; x += 4) {
      x1 = _mm_loadu_si128((__m128i*)&src[x]);
      x2 = _mm_loadu_si128((__m128i*)&src[x + 1]);
      x3 = _mm_loadu_si128((__m128i*)&src[x + 2 * 1]);
      x4 = _mm_loadu_si128((__m128i*)&src[x + 3 * 1]);
      t1 = _mm_unpackhi_epi16(x1, x2);
      x1 = _mm_unpacklo_epi16(x1, x2);
      t2 = _mm_unpackhi_epi16(x3, x4);
      x2 = _mm_unpacklo_epi16(x3, x4);
      x1 = _mm_madd_epi16(x1, f1);
      t1 = _mm_madd_epi16(t1, f1);
      x2 = _mm_madd_epi16(x2, f2);
      t2 = _mm_madd_epi16(t2, f2);
      x1 = _mm_add_epi32(x1, x2);
      t1 = _mm_add_epi32(t1, t2);
      t1 = _mm_srai_epi32(t1, 2);
      x1 = _mm_srai_epi32(x1, 2);
      x1 = _mm_packs_epi32(x1, t1);
      _mm_storel_epi64((__m128i*)&dst[x], x1);
    }
    src += srcstride;
    dst += 128;
  }
}

static void
oh_hevc_put_hevc_bi1_epel_h4_10_sse(uint16_t* _dst,
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
  __m128i x1, x2, x3, x4, t1, t2, f1, f2;
  uint16_t* src = ((uint16_t*)_src) - 1;
  ptrdiff_t srcstride = _srcstride;
  const __m128i offset = _mm_set1_epi16(1 << 10);
  uint16_t* dst = (uint16_t*)_dst;
  const int dststride = _dststride;
  f1 = _mm_load_si128((__m128i*)oh_hevc_epel_filters_sse[mx - 1][0]);
  f2 = _mm_load_si128((__m128i*)oh_hevc_epel_filters_sse[mx - 1][1]);
  for (y = 0; y < height; y++) {
    for (x = 0; x < width; x += 4) {
      __m128i r5;
      x1 = _mm_loadu_si128((__m128i*)&src[x]);
      x2 = _mm_loadu_si128((__m128i*)&src[x + 1]);
      x3 = _mm_loadu_si128((__m128i*)&src[x + 2 * 1]);
      x4 = _mm_loadu_si128((__m128i*)&src[x + 3 * 1]);
      t1 = _mm_unpackhi_epi16(x1, x2);
      x1 = _mm_unpacklo_epi16(x1, x2);
      t2 = _mm_unpackhi_epi16(x3, x4);
      x2 = _mm_unpacklo_epi16(x3, x4);
      x1 = _mm_madd_epi16(x1, f1);
      t1 = _mm_madd_epi16(t1, f1);
      x2 = _mm_madd_epi16(x2, f2);
      t2 = _mm_madd_epi16(t2, f2);
      x1 = _mm_add_epi32(x1, x2);
      t1 = _mm_add_epi32(t1, t2);
      t1 = _mm_srai_epi32(t1, 2);
      x1 = _mm_srai_epi32(x1, 2);
      x1 = _mm_packs_epi32(x1, t1);
      r5 = _mm_loadl_epi64((__m128i*)&src2[x]);
      x1 = _mm_adds_epi16(x1, r5);
      x1 = _mm_mulhrs_epi16(x1, offset);
      x1 = _mm_max_epi16(x1, _mm_setzero_si128());
      x1 = _mm_min_epi16(x1, _mm_set1_epi16(0x03FF));
      _mm_storel_epi64((__m128i*)&dst[x], x1);
    }
    src += srcstride;
    src2 += 128;
    dst += dststride;
  }
}

static void
oh_hevc_put_hevc_uni_epel_h4_10_sse(uint16_t* _dst,
                                    ptrdiff_t _dststride,
                                    const uint16_t* _src,
                                    ptrdiff_t _srcstride,
                                    int height,
                                    intptr_t mx,
                                    intptr_t my,
                                    int width)
{
  int x, y;
  __m128i x1, x2, x3, x4, t1, t2, f1, f2;
  uint16_t* src = ((uint16_t*)_src) - 1;
  ptrdiff_t srcstride = _srcstride;
  const __m128i offset = _mm_set1_epi16(1 << (10 + 1));
  uint16_t* dst = (uint16_t*)_dst;
  const int dststride = _dststride;
  f1 = _mm_load_si128((__m128i*)oh_hevc_epel_filters_sse[mx - 1][0]);
  f2 = _mm_load_si128((__m128i*)oh_hevc_epel_filters_sse[mx - 1][1]);
  for (y = 0; y < height; y++) {
    for (x = 0; x < width; x += 4) {
      x1 = _mm_loadu_si128((__m128i*)&src[x]);
      x2 = _mm_loadu_si128((__m128i*)&src[x + 1]);
      x3 = _mm_loadu_si128((__m128i*)&src[x + 2 * 1]);
      x4 = _mm_loadu_si128((__m128i*)&src[x + 3 * 1]);
      t1 = _mm_unpackhi_epi16(x1, x2);
      x1 = _mm_unpacklo_epi16(x1, x2);
      t2 = _mm_unpackhi_epi16(x3, x4);
      x2 = _mm_unpacklo_epi16(x3, x4);
      x1 = _mm_madd_epi16(x1, f1);
      t1 = _mm_madd_epi16(t1, f1);
      x2 = _mm_madd_epi16(x2, f2);
      t2 = _mm_madd_epi16(t2, f2);
      x1 = _mm_add_epi32(x1, x2);
      t1 = _mm_add_epi32(t1, t2);
      t1 = _mm_srai_epi32(t1, 2);
      x1 = _mm_srai_epi32(x1, 2);
      x1 = _mm_packs_epi32(x1, t1);
      x1 = _mm_mulhrs_epi16(x1, offset);
      x1 = _mm_max_epi16(x1, _mm_setzero_si128());
      x1 = _mm_min_epi16(x1, _mm_set1_epi16(0x03FF));
      _mm_storel_epi64((__m128i*)&dst[x], x1);
    }
    src += srcstride;
    dst += dststride;
  }
}

static void
oh_hevc_put_hevc_bi0_epel_h8_10_sse(int16_t* dst,
                                    const uint16_t* _src,
                                    ptrdiff_t _srcstride,
                                    int height,
                                    intptr_t mx,
                                    intptr_t my,
                                    int width)
{
  int x, y;
  __m128i x1, x2, x3, x4, t1, t2, f1, f2;
  uint16_t* src = ((uint16_t*)_src) - 1;
  ptrdiff_t srcstride = _srcstride;
  f1 = _mm_load_si128((__m128i*)oh_hevc_epel_filters_sse[mx - 1][0]);
  f2 = _mm_load_si128((__m128i*)oh_hevc_epel_filters_sse[mx - 1][1]);
  for (y = 0; y < height; y++) {
    for (x = 0; x < width; x += 8) {
      x1 = _mm_loadu_si128((__m128i*)&src[x]);
      x2 = _mm_loadu_si128((__m128i*)&src[x + 1]);
      x3 = _mm_loadu_si128((__m128i*)&src[x + 2 * 1]);
      x4 = _mm_loadu_si128((__m128i*)&src[x + 3 * 1]);
      t1 = _mm_unpackhi_epi16(x1, x2);
      x1 = _mm_unpacklo_epi16(x1, x2);
      t2 = _mm_unpackhi_epi16(x3, x4);
      x2 = _mm_unpacklo_epi16(x3, x4);
      x1 = _mm_madd_epi16(x1, f1);
      t1 = _mm_madd_epi16(t1, f1);
      x2 = _mm_madd_epi16(x2, f2);
      t2 = _mm_madd_epi16(t2, f2);
      x1 = _mm_add_epi32(x1, x2);
      t1 = _mm_add_epi32(t1, t2);
      t1 = _mm_srai_epi32(t1, 2);
      x1 = _mm_srai_epi32(x1, 2);
      x1 = _mm_packs_epi32(x1, t1);
      _mm_storeu_si128((__m128i*)&dst[x], x1);
    }
    src += srcstride;
    dst += 128;
  }
}

static void
oh_hevc_put_hevc_bi1_epel_h8_10_sse(uint16_t* _dst,
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
  __m128i x1, x2, x3, x4, t1, t2, f1, f2;
  uint16_t* src = ((uint16_t*)_src) - 1;
  ptrdiff_t srcstride = _srcstride;
  const __m128i offset = _mm_set1_epi16(1 << 10);
  uint16_t* dst = (uint16_t*)_dst;
  const int dststride = _dststride;
  f1 = _mm_load_si128((__m128i*)oh_hevc_epel_filters_sse[mx - 1][0]);
  f2 = _mm_load_si128((__m128i*)oh_hevc_epel_filters_sse[mx - 1][1]);
  for (y = 0; y < height; y++) {
    for (x = 0; x < width; x += 8) {
      __m128i r5;
      x1 = _mm_loadu_si128((__m128i*)&src[x]);
      x2 = _mm_loadu_si128((__m128i*)&src[x + 1]);
      x3 = _mm_loadu_si128((__m128i*)&src[x + 2 * 1]);
      x4 = _mm_loadu_si128((__m128i*)&src[x + 3 * 1]);
      t1 = _mm_unpackhi_epi16(x1, x2);
      x1 = _mm_unpacklo_epi16(x1, x2);
      t2 = _mm_unpackhi_epi16(x3, x4);
      x2 = _mm_unpacklo_epi16(x3, x4);
      x1 = _mm_madd_epi16(x1, f1);
      t1 = _mm_madd_epi16(t1, f1);
      x2 = _mm_madd_epi16(x2, f2);
      t2 = _mm_madd_epi16(t2, f2);
      x1 = _mm_add_epi32(x1, x2);
      t1 = _mm_add_epi32(t1, t2);
      t1 = _mm_srai_epi32(t1, 2);
      x1 = _mm_srai_epi32(x1, 2);
      x1 = _mm_packs_epi32(x1, t1);
      r5 = _mm_load_si128((__m128i*)&src2[x]);
      x1 = _mm_adds_epi16(x1, r5);
      x1 = _mm_mulhrs_epi16(x1, offset);
      x1 = _mm_max_epi16(x1, _mm_setzero_si128());
      x1 = _mm_min_epi16(x1, _mm_set1_epi16(0x03FF));
      _mm_storeu_si128((__m128i*)&dst[x], x1);
    }
    src += srcstride;
    src2 += 128;
    dst += dststride;
  }
}

static void
oh_hevc_put_hevc_uni_epel_h8_10_sse(uint16_t* _dst,
                                    ptrdiff_t _dststride,
                                    const uint16_t* _src,
                                    ptrdiff_t _srcstride,
                                    int height,
                                    intptr_t mx,
                                    intptr_t my,
                                    int width)
{
  int x, y;
  __m128i x1, x2, x3, x4, t1, t2, f1, f2;
  uint16_t* src = ((uint16_t*)_src) - 1;
  ptrdiff_t srcstride = _srcstride;
  const __m128i offset = _mm_set1_epi16(1 << (10 + 1));
  uint16_t* dst = (uint16_t*)_dst;
  const int dststride = _dststride;
  f1 = _mm_load_si128((__m128i*)oh_hevc_epel_filters_sse[mx - 1][0]);
  f2 = _mm_load_si128((__m128i*)oh_hevc_epel_filters_sse[mx - 1][1]);
  for (y = 0; y < height; y++) {
    for (x = 0; x < width; x += 8) {
      x1 = _mm_loadu_si128((__m128i*)&src[x]);
      x2 = _mm_loadu_si128((__m128i*)&src[x + 1]);
      x3 = _mm_loadu_si128((__m128i*)&src[x + 2 * 1]);
      x4 = _mm_loadu_si128((__m128i*)&src[x + 3 * 1]);
      t1 = _mm_unpackhi_epi16(x1, x2);
      x1 = _mm_unpacklo_epi16(x1, x2);
      t2 = _mm_unpackhi_epi16(x3, x4);
      x2 = _mm_unpacklo_epi16(x3, x4);
      x1 = _mm_madd_epi16(x1, f1);
      t1 = _mm_madd_epi16(t1, f1);
      x2 = _mm_madd_epi16(x2, f2);
      t2 = _mm_madd_epi16(t2, f2);
      x1 = _mm_add_epi32(x1, x2);
      t1 = _mm_add_epi32(t1, t2);
      t1 = _mm_srai_epi32(t1, 2);
      x1 = _mm_srai_epi32(x1, 2);
      x1 = _mm_packs_epi32(x1, t1);
      x1 = _mm_mulhrs_epi16(x1, offset);
      x1 = _mm_max_epi16(x1, _mm_setzero_si128());
      x1 = _mm_min_epi16(x1, _mm_set1_epi16(0x03FF));
      _mm_storeu_si128((__m128i*)&dst[x], x1);
    }
    src += srcstride;
    dst += dststride;
  }
}

static void
oh_hevc_put_hevc_bi0_epel_v4_10_sse(int16_t* dst,
                                    const uint16_t* _src,
                                    ptrdiff_t _srcstride,
                                    int height,
                                    intptr_t mx,
                                    intptr_t my,
                                    int width)
{
  int x, y;
  __m128i x1, x2, x3, x4, t1, t2, f1, f2;
  ptrdiff_t srcstride = _srcstride;
  uint16_t* src = ((uint16_t*)_src) - srcstride;
  f1 = _mm_load_si128((__m128i*)oh_hevc_epel_filters_sse[my - 1][0]);
  f2 = _mm_load_si128((__m128i*)oh_hevc_epel_filters_sse[my - 1][1]);
  for (y = 0; y < height; y++) {
    for (x = 0; x < width; x += 4) {
      x1 = _mm_loadu_si128((__m128i*)&src[x]);
      x2 = _mm_loadu_si128((__m128i*)&src[x + srcstride]);
      x3 = _mm_loadu_si128((__m128i*)&src[x + 2 * srcstride]);
      x4 = _mm_loadu_si128((__m128i*)&src[x + 3 * srcstride]);
      t1 = _mm_unpackhi_epi16(x1, x2);
      x1 = _mm_unpacklo_epi16(x1, x2);
      t2 = _mm_unpackhi_epi16(x3, x4);
      x2 = _mm_unpacklo_epi16(x3, x4);
      x1 = _mm_madd_epi16(x1, f1);
      t1 = _mm_madd_epi16(t1, f1);
      x2 = _mm_madd_epi16(x2, f2);
      t2 = _mm_madd_epi16(t2, f2);
      x1 = _mm_add_epi32(x1, x2);
      t1 = _mm_add_epi32(t1, t2);
      t1 = _mm_srai_epi32(t1, 2);
      x1 = _mm_srai_epi32(x1, 2);
      x1 = _mm_packs_epi32(x1, t1);
      _mm_storel_epi64((__m128i*)&dst[x], x1);
    }
    src += srcstride;
    dst += 128;
  }
}

static void
oh_hevc_put_hevc_bi1_epel_v4_10_sse(uint16_t* _dst,
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
  __m128i x1, x2, x3, x4, t1, t2, f1, f2;
  ptrdiff_t srcstride = _srcstride;
  uint16_t* src = ((uint16_t*)_src) - srcstride;
  const __m128i offset = _mm_set1_epi16(1 << 10);
  uint16_t* dst = (uint16_t*)_dst;
  const int dststride = _dststride;
  f1 = _mm_load_si128((__m128i*)oh_hevc_epel_filters_sse[my - 1][0]);
  f2 = _mm_load_si128((__m128i*)oh_hevc_epel_filters_sse[my - 1][1]);
  for (y = 0; y < height; y++) {
    for (x = 0; x < width; x += 4) {
      __m128i r5;
      x1 = _mm_loadu_si128((__m128i*)&src[x]);
      x2 = _mm_loadu_si128((__m128i*)&src[x + srcstride]);
      x3 = _mm_loadu_si128((__m128i*)&src[x + 2 * srcstride]);
      x4 = _mm_loadu_si128((__m128i*)&src[x + 3 * srcstride]);
      t1 = _mm_unpackhi_epi16(x1, x2);
      x1 = _mm_unpacklo_epi16(x1, x2);
      t2 = _mm_unpackhi_epi16(x3, x4);
      x2 = _mm_unpacklo_epi16(x3, x4);
      x1 = _mm_madd_epi16(x1, f1);
      t1 = _mm_madd_epi16(t1, f1);
      x2 = _mm_madd_epi16(x2, f2);
      t2 = _mm_madd_epi16(t2, f2);
      x1 = _mm_add_epi32(x1, x2);
      t1 = _mm_add_epi32(t1, t2);
      t1 = _mm_srai_epi32(t1, 2);
      x1 = _mm_srai_epi32(x1, 2);
      x1 = _mm_packs_epi32(x1, t1);
      r5 = _mm_loadl_epi64((__m128i*)&src2[x]);
      x1 = _mm_adds_epi16(x1, r5);
      x1 = _mm_mulhrs_epi16(x1, offset);
      x1 = _mm_max_epi16(x1, _mm_setzero_si128());
      x1 = _mm_min_epi16(x1, _mm_set1_epi16(0x03FF));
      _mm_storel_epi64((__m128i*)&dst[x], x1);
    }
    src += srcstride;
    src2 += 128;
    dst += dststride;
  }
}

static void
oh_hevc_put_hevc_uni_epel_v4_10_sse(uint16_t* _dst,
                                    ptrdiff_t _dststride,
                                    const uint16_t* _src,
                                    ptrdiff_t _srcstride,
                                    int height,
                                    intptr_t mx,
                                    intptr_t my,
                                    int width)
{
  int x, y;
  __m128i x1, x2, x3, x4, t1, t2, f1, f2;
  ptrdiff_t srcstride = _srcstride;
  uint16_t* src = ((uint16_t*)_src) - srcstride;
  const __m128i offset = _mm_set1_epi16(1 << (10 + 1));
  uint16_t* dst = (uint16_t*)_dst;
  const int dststride = _dststride;
  f1 = _mm_load_si128((__m128i*)oh_hevc_epel_filters_sse[my - 1][0]);
  f2 = _mm_load_si128((__m128i*)oh_hevc_epel_filters_sse[my - 1][1]);
  for (y = 0; y < height; y++) {
    for (x = 0; x < width; x += 4) {
      x1 = _mm_loadu_si128((__m128i*)&src[x]);
      x2 = _mm_loadu_si128((__m128i*)&src[x + srcstride]);
      x3 = _mm_loadu_si128((__m128i*)&src[x + 2 * srcstride]);
      x4 = _mm_loadu_si128((__m128i*)&src[x + 3 * srcstride]);
      t1 = _mm_unpackhi_epi16(x1, x2);
      x1 = _mm_unpacklo_epi16(x1, x2);
      t2 = _mm_unpackhi_epi16(x3, x4);
      x2 = _mm_unpacklo_epi16(x3, x4);
      x1 = _mm_madd_epi16(x1, f1);
      t1 = _mm_madd_epi16(t1, f1);
      x2 = _mm_madd_epi16(x2, f2);
      t2 = _mm_madd_epi16(t2, f2);
      x1 = _mm_add_epi32(x1, x2);
      t1 = _mm_add_epi32(t1, t2);
      t1 = _mm_srai_epi32(t1, 2);
      x1 = _mm_srai_epi32(x1, 2);
      x1 = _mm_packs_epi32(x1, t1);
      x1 = _mm_mulhrs_epi16(x1, offset);
      x1 = _mm_max_epi16(x1, _mm_setzero_si128());
      x1 = _mm_min_epi16(x1, _mm_set1_epi16(0x03FF));
      _mm_storel_epi64((__m128i*)&dst[x], x1);
    }
    src += srcstride;
    dst += dststride;
  }
}

static void
oh_hevc_put_hevc_bi0_epel_v8_10_sse(int16_t* dst,
                                    const uint16_t* _src,
                                    ptrdiff_t _srcstride,
                                    int height,
                                    intptr_t mx,
                                    intptr_t my,
                                    int width)
{
  int x, y;
  __m128i x1, x2, x3, x4, t1, t2, f1, f2;
  ptrdiff_t srcstride = _srcstride;
  uint16_t* src = ((uint16_t*)_src) - srcstride;
  f1 = _mm_load_si128((__m128i*)oh_hevc_epel_filters_sse[my - 1][0]);
  f2 = _mm_load_si128((__m128i*)oh_hevc_epel_filters_sse[my - 1][1]);
  for (y = 0; y < height; y++) {
    for (x = 0; x < width; x += 8) {
      x1 = _mm_loadu_si128((__m128i*)&src[x]);
      x2 = _mm_loadu_si128((__m128i*)&src[x + srcstride]);
      x3 = _mm_loadu_si128((__m128i*)&src[x + 2 * srcstride]);
      x4 = _mm_loadu_si128((__m128i*)&src[x + 3 * srcstride]);
      t1 = _mm_unpackhi_epi16(x1, x2);
      x1 = _mm_unpacklo_epi16(x1, x2);
      t2 = _mm_unpackhi_epi16(x3, x4);
      x2 = _mm_unpacklo_epi16(x3, x4);
      x1 = _mm_madd_epi16(x1, f1);
      t1 = _mm_madd_epi16(t1, f1);
      x2 = _mm_madd_epi16(x2, f2);
      t2 = _mm_madd_epi16(t2, f2);
      x1 = _mm_add_epi32(x1, x2);
      t1 = _mm_add_epi32(t1, t2);
      t1 = _mm_srai_epi32(t1, 2);
      x1 = _mm_srai_epi32(x1, 2);
      x1 = _mm_packs_epi32(x1, t1);
      _mm_storeu_si128((__m128i*)&dst[x], x1);
    }
    src += srcstride;
    dst += 128;
  }
}

static void
oh_hevc_put_hevc_bi1_epel_v8_10_sse(uint16_t* _dst,
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
  __m128i x1, x2, x3, x4, t1, t2, f1, f2;
  ptrdiff_t srcstride = _srcstride;
  uint16_t* src = ((uint16_t*)_src) - srcstride;
  const __m128i offset = _mm_set1_epi16(1 << 10);
  uint16_t* dst = (uint16_t*)_dst;
  const int dststride = _dststride;
  f1 = _mm_load_si128((__m128i*)oh_hevc_epel_filters_sse[my - 1][0]);
  f2 = _mm_load_si128((__m128i*)oh_hevc_epel_filters_sse[my - 1][1]);
  for (y = 0; y < height; y++) {
    for (x = 0; x < width; x += 8) {
      __m128i r5;
      x1 = _mm_loadu_si128((__m128i*)&src[x]);
      x2 = _mm_loadu_si128((__m128i*)&src[x + srcstride]);
      x3 = _mm_loadu_si128((__m128i*)&src[x + 2 * srcstride]);
      x4 = _mm_loadu_si128((__m128i*)&src[x + 3 * srcstride]);
      t1 = _mm_unpackhi_epi16(x1, x2);
      x1 = _mm_unpacklo_epi16(x1, x2);
      t2 = _mm_unpackhi_epi16(x3, x4);
      x2 = _mm_unpacklo_epi16(x3, x4);
      x1 = _mm_madd_epi16(x1, f1);
      t1 = _mm_madd_epi16(t1, f1);
      x2 = _mm_madd_epi16(x2, f2);
      t2 = _mm_madd_epi16(t2, f2);
      x1 = _mm_add_epi32(x1, x2);
      t1 = _mm_add_epi32(t1, t2);
      t1 = _mm_srai_epi32(t1, 2);
      x1 = _mm_srai_epi32(x1, 2);
      x1 = _mm_packs_epi32(x1, t1);
      r5 = _mm_load_si128((__m128i*)&src2[x]);
      x1 = _mm_adds_epi16(x1, r5);
      x1 = _mm_mulhrs_epi16(x1, offset);
      x1 = _mm_max_epi16(x1, _mm_setzero_si128());
      x1 = _mm_min_epi16(x1, _mm_set1_epi16(0x03FF));
      _mm_storeu_si128((__m128i*)&dst[x], x1);
    }
    src += srcstride;
    src2 += 128;
    dst += dststride;
  }
}

static void
oh_hevc_put_hevc_uni_epel_v8_10_sse(uint16_t* _dst,
                                    ptrdiff_t _dststride,
                                    const uint16_t* _src,
                                    ptrdiff_t _srcstride,
                                    int height,
                                    intptr_t mx,
                                    intptr_t my,
                                    int width)
{
  int x, y;
  __m128i x1, x2, x3, x4, t1, t2, f1, f2;
  ptrdiff_t srcstride = _srcstride;
  uint16_t* src = ((uint16_t*)_src) - srcstride;
  const __m128i offset = _mm_set1_epi16(1 << (10 + 1));
  uint16_t* dst = (uint16_t*)_dst;
  const int dststride = _dststride;
  f1 = _mm_load_si128((__m128i*)oh_hevc_epel_filters_sse[my - 1][0]);
  f2 = _mm_load_si128((__m128i*)oh_hevc_epel_filters_sse[my - 1][1]);
  for (y = 0; y < height; y++) {
    for (x = 0; x < width; x += 8) {
      x1 = _mm_loadu_si128((__m128i*)&src[x]);
      x2 = _mm_loadu_si128((__m128i*)&src[x + srcstride]);
      x3 = _mm_loadu_si128((__m128i*)&src[x + 2 * srcstride]);
      x4 = _mm_loadu_si128((__m128i*)&src[x + 3 * srcstride]);
      t1 = _mm_unpackhi_epi16(x1, x2);
      x1 = _mm_unpacklo_epi16(x1, x2);
      t2 = _mm_unpackhi_epi16(x3, x4);
      x2 = _mm_unpacklo_epi16(x3, x4);
      x1 = _mm_madd_epi16(x1, f1);
      t1 = _mm_madd_epi16(t1, f1);
      x2 = _mm_madd_epi16(x2, f2);
      t2 = _mm_madd_epi16(t2, f2);
      x1 = _mm_add_epi32(x1, x2);
      t1 = _mm_add_epi32(t1, t2);
      t1 = _mm_srai_epi32(t1, 2);
      x1 = _mm_srai_epi32(x1, 2);
      x1 = _mm_packs_epi32(x1, t1);
      x1 = _mm_mulhrs_epi16(x1, offset);
      x1 = _mm_max_epi16(x1, _mm_setzero_si128());
      x1 = _mm_min_epi16(x1, _mm_set1_epi16(0x03FF));
      _mm_storeu_si128((__m128i*)&dst[x], x1);
    }
    src += srcstride;
    dst += dststride;
  }
}

static void
oh_hevc_put_hevc_bi0_epel_hv4_10_sse(int16_t* dst,
                                     const uint16_t* _src,
                                     ptrdiff_t _srcstride,
                                     int height,
                                     intptr_t mx,
                                     intptr_t my,
                                     int width)
{
  int x, y;
  __m128i x1, x2, x3, x4, t1, t2, f1, f2, f3, f4, r1, r2, r3, r4;
  int16_t* dst_bis = dst;
  uint16_t *src_bis, *src = ((uint16_t*)_src) - 1;
  ptrdiff_t srcstride = _srcstride;
  src -= 1 * srcstride;
  src_bis = src;
  f1 = _mm_load_si128((__m128i*)oh_hevc_epel_filters_sse[mx - 1][0]);
  f2 = _mm_load_si128((__m128i*)oh_hevc_epel_filters_sse[mx - 1][1]);
  f3 = _mm_load_si128((__m128i*)oh_hevc_epel_filters_sse[my - 1][0]);
  f4 = _mm_load_si128((__m128i*)oh_hevc_epel_filters_sse[my - 1][1]);
  for (x = 0; x < width; x += 4) {
    x1 = _mm_loadu_si128((__m128i*)&src[x]);
    x2 = _mm_loadu_si128((__m128i*)&src[x + 1]);
    x3 = _mm_loadu_si128((__m128i*)&src[x + 2 * 1]);
    x4 = _mm_loadu_si128((__m128i*)&src[x + 3 * 1]);
    t1 = _mm_unpackhi_epi16(x1, x2);
    x1 = _mm_unpacklo_epi16(x1, x2);
    t2 = _mm_unpackhi_epi16(x3, x4);
    x2 = _mm_unpacklo_epi16(x3, x4);
    x1 = _mm_madd_epi16(x1, f1);
    t1 = _mm_madd_epi16(t1, f1);
    x2 = _mm_madd_epi16(x2, f2);
    t2 = _mm_madd_epi16(t2, f2);
    x1 = _mm_add_epi32(x1, x2);
    t1 = _mm_add_epi32(t1, t2);
    t1 = _mm_srai_epi32(t1, 2);
    x1 = _mm_srai_epi32(x1, 2);
    r1 = _mm_packs_epi32(x1, t1);
    src += srcstride;
    x1 = _mm_loadu_si128((__m128i*)&src[x]);
    x2 = _mm_loadu_si128((__m128i*)&src[x + 1]);
    x3 = _mm_loadu_si128((__m128i*)&src[x + 2 * 1]);
    x4 = _mm_loadu_si128((__m128i*)&src[x + 3 * 1]);
    t1 = _mm_unpackhi_epi16(x1, x2);
    x1 = _mm_unpacklo_epi16(x1, x2);
    t2 = _mm_unpackhi_epi16(x3, x4);
    x2 = _mm_unpacklo_epi16(x3, x4);
    x1 = _mm_madd_epi16(x1, f1);
    t1 = _mm_madd_epi16(t1, f1);
    x2 = _mm_madd_epi16(x2, f2);
    t2 = _mm_madd_epi16(t2, f2);
    x1 = _mm_add_epi32(x1, x2);
    t1 = _mm_add_epi32(t1, t2);
    t1 = _mm_srai_epi32(t1, 2);
    x1 = _mm_srai_epi32(x1, 2);
    r2 = _mm_packs_epi32(x1, t1);
    src += srcstride;
    x1 = _mm_loadu_si128((__m128i*)&src[x]);
    x2 = _mm_loadu_si128((__m128i*)&src[x + 1]);
    x3 = _mm_loadu_si128((__m128i*)&src[x + 2 * 1]);
    x4 = _mm_loadu_si128((__m128i*)&src[x + 3 * 1]);
    t1 = _mm_unpackhi_epi16(x1, x2);
    x1 = _mm_unpacklo_epi16(x1, x2);
    t2 = _mm_unpackhi_epi16(x3, x4);
    x2 = _mm_unpacklo_epi16(x3, x4);
    x1 = _mm_madd_epi16(x1, f1);
    t1 = _mm_madd_epi16(t1, f1);
    x2 = _mm_madd_epi16(x2, f2);
    t2 = _mm_madd_epi16(t2, f2);
    x1 = _mm_add_epi32(x1, x2);
    t1 = _mm_add_epi32(t1, t2);
    t1 = _mm_srai_epi32(t1, 2);
    x1 = _mm_srai_epi32(x1, 2);
    r3 = _mm_packs_epi32(x1, t1);
    src += srcstride;
    for (y = 0; y < height; y++) {
      x1 = _mm_loadu_si128((__m128i*)&src[x]);
      x2 = _mm_loadu_si128((__m128i*)&src[x + 1]);
      x3 = _mm_loadu_si128((__m128i*)&src[x + 2 * 1]);
      x4 = _mm_loadu_si128((__m128i*)&src[x + 3 * 1]);
      t1 = _mm_unpackhi_epi16(x1, x2);
      x1 = _mm_unpacklo_epi16(x1, x2);
      t2 = _mm_unpackhi_epi16(x3, x4);
      x2 = _mm_unpacklo_epi16(x3, x4);
      x1 = _mm_madd_epi16(x1, f1);
      t1 = _mm_madd_epi16(t1, f1);
      x2 = _mm_madd_epi16(x2, f2);
      t2 = _mm_madd_epi16(t2, f2);
      x1 = _mm_add_epi32(x1, x2);
      t1 = _mm_add_epi32(t1, t2);
      t1 = _mm_srai_epi32(t1, 2);
      x1 = _mm_srai_epi32(x1, 2);
      r4 = _mm_packs_epi32(x1, t1);
      src += srcstride;
      t1 = _mm_unpackhi_epi16(r1, r2);
      x1 = _mm_unpacklo_epi16(r1, r2);
      t2 = _mm_unpackhi_epi16(r3, r4);
      x2 = _mm_unpacklo_epi16(r3, r4);
      x1 = _mm_madd_epi16(x1, f3);
      t1 = _mm_madd_epi16(t1, f3);
      x2 = _mm_madd_epi16(x2, f4);
      t2 = _mm_madd_epi16(t2, f4);
      x1 = _mm_add_epi32(x1, x2);
      t1 = _mm_add_epi32(t1, t2);
      t1 = _mm_srai_epi32(t1, 6);
      x1 = _mm_srai_epi32(x1, 6);
      x1 = _mm_packs_epi32(x1, t1);
      _mm_storel_epi64((__m128i*)&dst[x], x1);
      dst += 128;
      r1 = r2;
      r2 = r3;
      r3 = r4;
    }
    src = src_bis;
    dst = dst_bis;
  }
}

static void
oh_hevc_put_hevc_bi1_epel_hv4_10_sse(uint16_t* _dst,
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
  __m128i x1, x2, x3, x4, t1, t2, f1, f2, f3, f4, r1, r2, r3, r4;
  const __m128i offset = _mm_set1_epi16(1 << 10);
  uint16_t* dst = (uint16_t*)_dst;
  const int dststride = _dststride;
  const int16_t* src2_bis = src2;
  uint16_t* dst_bis = dst;
  uint16_t *src_bis, *src = ((uint16_t*)_src) - 1;
  ptrdiff_t srcstride = _srcstride;
  src -= 1 * srcstride;
  src_bis = src;
  f1 = _mm_load_si128((__m128i*)oh_hevc_epel_filters_sse[mx - 1][0]);
  f2 = _mm_load_si128((__m128i*)oh_hevc_epel_filters_sse[mx - 1][1]);
  f3 = _mm_load_si128((__m128i*)oh_hevc_epel_filters_sse[my - 1][0]);
  f4 = _mm_load_si128((__m128i*)oh_hevc_epel_filters_sse[my - 1][1]);
  for (x = 0; x < width; x += 4) {
    __m128i r5;
    x1 = _mm_loadu_si128((__m128i*)&src[x]);
    x2 = _mm_loadu_si128((__m128i*)&src[x + 1]);
    x3 = _mm_loadu_si128((__m128i*)&src[x + 2 * 1]);
    x4 = _mm_loadu_si128((__m128i*)&src[x + 3 * 1]);
    t1 = _mm_unpackhi_epi16(x1, x2);
    x1 = _mm_unpacklo_epi16(x1, x2);
    t2 = _mm_unpackhi_epi16(x3, x4);
    x2 = _mm_unpacklo_epi16(x3, x4);
    x1 = _mm_madd_epi16(x1, f1);
    t1 = _mm_madd_epi16(t1, f1);
    x2 = _mm_madd_epi16(x2, f2);
    t2 = _mm_madd_epi16(t2, f2);
    x1 = _mm_add_epi32(x1, x2);
    t1 = _mm_add_epi32(t1, t2);
    t1 = _mm_srai_epi32(t1, 2);
    x1 = _mm_srai_epi32(x1, 2);
    r1 = _mm_packs_epi32(x1, t1);
    src += srcstride;
    x1 = _mm_loadu_si128((__m128i*)&src[x]);
    x2 = _mm_loadu_si128((__m128i*)&src[x + 1]);
    x3 = _mm_loadu_si128((__m128i*)&src[x + 2 * 1]);
    x4 = _mm_loadu_si128((__m128i*)&src[x + 3 * 1]);
    t1 = _mm_unpackhi_epi16(x1, x2);
    x1 = _mm_unpacklo_epi16(x1, x2);
    t2 = _mm_unpackhi_epi16(x3, x4);
    x2 = _mm_unpacklo_epi16(x3, x4);
    x1 = _mm_madd_epi16(x1, f1);
    t1 = _mm_madd_epi16(t1, f1);
    x2 = _mm_madd_epi16(x2, f2);
    t2 = _mm_madd_epi16(t2, f2);
    x1 = _mm_add_epi32(x1, x2);
    t1 = _mm_add_epi32(t1, t2);
    t1 = _mm_srai_epi32(t1, 2);
    x1 = _mm_srai_epi32(x1, 2);
    r2 = _mm_packs_epi32(x1, t1);
    src += srcstride;
    x1 = _mm_loadu_si128((__m128i*)&src[x]);
    x2 = _mm_loadu_si128((__m128i*)&src[x + 1]);
    x3 = _mm_loadu_si128((__m128i*)&src[x + 2 * 1]);
    x4 = _mm_loadu_si128((__m128i*)&src[x + 3 * 1]);
    t1 = _mm_unpackhi_epi16(x1, x2);
    x1 = _mm_unpacklo_epi16(x1, x2);
    t2 = _mm_unpackhi_epi16(x3, x4);
    x2 = _mm_unpacklo_epi16(x3, x4);
    x1 = _mm_madd_epi16(x1, f1);
    t1 = _mm_madd_epi16(t1, f1);
    x2 = _mm_madd_epi16(x2, f2);
    t2 = _mm_madd_epi16(t2, f2);
    x1 = _mm_add_epi32(x1, x2);
    t1 = _mm_add_epi32(t1, t2);
    t1 = _mm_srai_epi32(t1, 2);
    x1 = _mm_srai_epi32(x1, 2);
    r3 = _mm_packs_epi32(x1, t1);
    src += srcstride;
    for (y = 0; y < height; y++) {
      x1 = _mm_loadu_si128((__m128i*)&src[x]);
      x2 = _mm_loadu_si128((__m128i*)&src[x + 1]);
      x3 = _mm_loadu_si128((__m128i*)&src[x + 2 * 1]);
      x4 = _mm_loadu_si128((__m128i*)&src[x + 3 * 1]);
      t1 = _mm_unpackhi_epi16(x1, x2);
      x1 = _mm_unpacklo_epi16(x1, x2);
      t2 = _mm_unpackhi_epi16(x3, x4);
      x2 = _mm_unpacklo_epi16(x3, x4);
      x1 = _mm_madd_epi16(x1, f1);
      t1 = _mm_madd_epi16(t1, f1);
      x2 = _mm_madd_epi16(x2, f2);
      t2 = _mm_madd_epi16(t2, f2);
      x1 = _mm_add_epi32(x1, x2);
      t1 = _mm_add_epi32(t1, t2);
      t1 = _mm_srai_epi32(t1, 2);
      x1 = _mm_srai_epi32(x1, 2);
      r4 = _mm_packs_epi32(x1, t1);
      src += srcstride;
      t1 = _mm_unpackhi_epi16(r1, r2);
      x1 = _mm_unpacklo_epi16(r1, r2);
      t2 = _mm_unpackhi_epi16(r3, r4);
      x2 = _mm_unpacklo_epi16(r3, r4);
      x1 = _mm_madd_epi16(x1, f3);
      t1 = _mm_madd_epi16(t1, f3);
      x2 = _mm_madd_epi16(x2, f4);
      t2 = _mm_madd_epi16(t2, f4);
      x1 = _mm_add_epi32(x1, x2);
      t1 = _mm_add_epi32(t1, t2);
      t1 = _mm_srai_epi32(t1, 6);
      x1 = _mm_srai_epi32(x1, 6);
      x1 = _mm_packs_epi32(x1, t1);
      r5 = _mm_loadl_epi64((__m128i*)&src2[x]);
      x1 = _mm_adds_epi16(x1, r5);
      x1 = _mm_mulhrs_epi16(x1, offset);
      x1 = _mm_max_epi16(x1, _mm_setzero_si128());
      x1 = _mm_min_epi16(x1, _mm_set1_epi16(0x03FF));
      _mm_storel_epi64((__m128i*)&dst[x], x1);
      src2 += 128;
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
oh_hevc_put_hevc_uni_epel_hv4_10_sse(uint16_t* _dst,
                                     ptrdiff_t _dststride,
                                     const uint16_t* _src,
                                     ptrdiff_t _srcstride,
                                     int height,
                                     intptr_t mx,
                                     intptr_t my,
                                     int width)
{
  int x, y;
  __m128i x1, x2, x3, x4, t1, t2, f1, f2, f3, f4, r1, r2, r3, r4;
  const __m128i offset = _mm_set1_epi16(1 << (10 + 1));
  uint16_t* dst = (uint16_t*)_dst;
  const int dststride = _dststride;
  uint16_t* dst_bis = dst;
  uint16_t *src_bis, *src = ((uint16_t*)_src) - 1;
  ptrdiff_t srcstride = _srcstride;
  src -= 1 * srcstride;
  src_bis = src;
  f1 = _mm_load_si128((__m128i*)oh_hevc_epel_filters_sse[mx - 1][0]);
  f2 = _mm_load_si128((__m128i*)oh_hevc_epel_filters_sse[mx - 1][1]);
  f3 = _mm_load_si128((__m128i*)oh_hevc_epel_filters_sse[my - 1][0]);
  f4 = _mm_load_si128((__m128i*)oh_hevc_epel_filters_sse[my - 1][1]);
  for (x = 0; x < width; x += 4) {
    x1 = _mm_loadu_si128((__m128i*)&src[x]);
    x2 = _mm_loadu_si128((__m128i*)&src[x + 1]);
    x3 = _mm_loadu_si128((__m128i*)&src[x + 2 * 1]);
    x4 = _mm_loadu_si128((__m128i*)&src[x + 3 * 1]);
    t1 = _mm_unpackhi_epi16(x1, x2);
    x1 = _mm_unpacklo_epi16(x1, x2);
    t2 = _mm_unpackhi_epi16(x3, x4);
    x2 = _mm_unpacklo_epi16(x3, x4);
    x1 = _mm_madd_epi16(x1, f1);
    t1 = _mm_madd_epi16(t1, f1);
    x2 = _mm_madd_epi16(x2, f2);
    t2 = _mm_madd_epi16(t2, f2);
    x1 = _mm_add_epi32(x1, x2);
    t1 = _mm_add_epi32(t1, t2);
    t1 = _mm_srai_epi32(t1, 2);
    x1 = _mm_srai_epi32(x1, 2);
    r1 = _mm_packs_epi32(x1, t1);
    src += srcstride;
    x1 = _mm_loadu_si128((__m128i*)&src[x]);
    x2 = _mm_loadu_si128((__m128i*)&src[x + 1]);
    x3 = _mm_loadu_si128((__m128i*)&src[x + 2 * 1]);
    x4 = _mm_loadu_si128((__m128i*)&src[x + 3 * 1]);
    t1 = _mm_unpackhi_epi16(x1, x2);
    x1 = _mm_unpacklo_epi16(x1, x2);
    t2 = _mm_unpackhi_epi16(x3, x4);
    x2 = _mm_unpacklo_epi16(x3, x4);
    x1 = _mm_madd_epi16(x1, f1);
    t1 = _mm_madd_epi16(t1, f1);
    x2 = _mm_madd_epi16(x2, f2);
    t2 = _mm_madd_epi16(t2, f2);
    x1 = _mm_add_epi32(x1, x2);
    t1 = _mm_add_epi32(t1, t2);
    t1 = _mm_srai_epi32(t1, 2);
    x1 = _mm_srai_epi32(x1, 2);
    r2 = _mm_packs_epi32(x1, t1);
    src += srcstride;
    x1 = _mm_loadu_si128((__m128i*)&src[x]);
    x2 = _mm_loadu_si128((__m128i*)&src[x + 1]);
    x3 = _mm_loadu_si128((__m128i*)&src[x + 2 * 1]);
    x4 = _mm_loadu_si128((__m128i*)&src[x + 3 * 1]);
    t1 = _mm_unpackhi_epi16(x1, x2);
    x1 = _mm_unpacklo_epi16(x1, x2);
    t2 = _mm_unpackhi_epi16(x3, x4);
    x2 = _mm_unpacklo_epi16(x3, x4);
    x1 = _mm_madd_epi16(x1, f1);
    t1 = _mm_madd_epi16(t1, f1);
    x2 = _mm_madd_epi16(x2, f2);
    t2 = _mm_madd_epi16(t2, f2);
    x1 = _mm_add_epi32(x1, x2);
    t1 = _mm_add_epi32(t1, t2);
    t1 = _mm_srai_epi32(t1, 2);
    x1 = _mm_srai_epi32(x1, 2);
    r3 = _mm_packs_epi32(x1, t1);
    src += srcstride;
    for (y = 0; y < height; y++) {
      x1 = _mm_loadu_si128((__m128i*)&src[x]);
      x2 = _mm_loadu_si128((__m128i*)&src[x + 1]);
      x3 = _mm_loadu_si128((__m128i*)&src[x + 2 * 1]);
      x4 = _mm_loadu_si128((__m128i*)&src[x + 3 * 1]);
      t1 = _mm_unpackhi_epi16(x1, x2);
      x1 = _mm_unpacklo_epi16(x1, x2);
      t2 = _mm_unpackhi_epi16(x3, x4);
      x2 = _mm_unpacklo_epi16(x3, x4);
      x1 = _mm_madd_epi16(x1, f1);
      t1 = _mm_madd_epi16(t1, f1);
      x2 = _mm_madd_epi16(x2, f2);
      t2 = _mm_madd_epi16(t2, f2);
      x1 = _mm_add_epi32(x1, x2);
      t1 = _mm_add_epi32(t1, t2);
      t1 = _mm_srai_epi32(t1, 2);
      x1 = _mm_srai_epi32(x1, 2);
      r4 = _mm_packs_epi32(x1, t1);
      src += srcstride;
      t1 = _mm_unpackhi_epi16(r1, r2);
      x1 = _mm_unpacklo_epi16(r1, r2);
      t2 = _mm_unpackhi_epi16(r3, r4);
      x2 = _mm_unpacklo_epi16(r3, r4);
      x1 = _mm_madd_epi16(x1, f3);
      t1 = _mm_madd_epi16(t1, f3);
      x2 = _mm_madd_epi16(x2, f4);
      t2 = _mm_madd_epi16(t2, f4);
      x1 = _mm_add_epi32(x1, x2);
      t1 = _mm_add_epi32(t1, t2);
      t1 = _mm_srai_epi32(t1, 6);
      x1 = _mm_srai_epi32(x1, 6);
      x1 = _mm_packs_epi32(x1, t1);
      x1 = _mm_mulhrs_epi16(x1, offset);
      x1 = _mm_max_epi16(x1, _mm_setzero_si128());
      x1 = _mm_min_epi16(x1, _mm_set1_epi16(0x03FF));
      _mm_storel_epi64((__m128i*)&dst[x], x1);
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
oh_hevc_put_hevc_bi0_epel_hv8_10_sse(int16_t* dst,
                                     const uint16_t* _src,
                                     ptrdiff_t _srcstride,
                                     int height,
                                     intptr_t mx,
                                     intptr_t my,
                                     int width)
{
  int x, y;
  __m128i x1, x2, x3, x4, t1, t2, f1, f2, f3, f4, r1, r2, r3, r4;
  int16_t* dst_bis = dst;
  uint16_t *src_bis, *src = ((uint16_t*)_src) - 1;
  ptrdiff_t srcstride = _srcstride;
  src -= 1 * srcstride;
  src_bis = src;
  f1 = _mm_load_si128((__m128i*)oh_hevc_epel_filters_sse[mx - 1][0]);
  f2 = _mm_load_si128((__m128i*)oh_hevc_epel_filters_sse[mx - 1][1]);
  f3 = _mm_load_si128((__m128i*)oh_hevc_epel_filters_sse[my - 1][0]);
  f4 = _mm_load_si128((__m128i*)oh_hevc_epel_filters_sse[my - 1][1]);
  for (x = 0; x < width; x += 8) {
    x1 = _mm_loadu_si128((__m128i*)&src[x]);
    x2 = _mm_loadu_si128((__m128i*)&src[x + 1]);
    x3 = _mm_loadu_si128((__m128i*)&src[x + 2 * 1]);
    x4 = _mm_loadu_si128((__m128i*)&src[x + 3 * 1]);
    t1 = _mm_unpackhi_epi16(x1, x2);
    x1 = _mm_unpacklo_epi16(x1, x2);
    t2 = _mm_unpackhi_epi16(x3, x4);
    x2 = _mm_unpacklo_epi16(x3, x4);
    x1 = _mm_madd_epi16(x1, f1);
    t1 = _mm_madd_epi16(t1, f1);
    x2 = _mm_madd_epi16(x2, f2);
    t2 = _mm_madd_epi16(t2, f2);
    x1 = _mm_add_epi32(x1, x2);
    t1 = _mm_add_epi32(t1, t2);
    t1 = _mm_srai_epi32(t1, 2);
    x1 = _mm_srai_epi32(x1, 2);
    r1 = _mm_packs_epi32(x1, t1);
    src += srcstride;
    x1 = _mm_loadu_si128((__m128i*)&src[x]);
    x2 = _mm_loadu_si128((__m128i*)&src[x + 1]);
    x3 = _mm_loadu_si128((__m128i*)&src[x + 2 * 1]);
    x4 = _mm_loadu_si128((__m128i*)&src[x + 3 * 1]);
    t1 = _mm_unpackhi_epi16(x1, x2);
    x1 = _mm_unpacklo_epi16(x1, x2);
    t2 = _mm_unpackhi_epi16(x3, x4);
    x2 = _mm_unpacklo_epi16(x3, x4);
    x1 = _mm_madd_epi16(x1, f1);
    t1 = _mm_madd_epi16(t1, f1);
    x2 = _mm_madd_epi16(x2, f2);
    t2 = _mm_madd_epi16(t2, f2);
    x1 = _mm_add_epi32(x1, x2);
    t1 = _mm_add_epi32(t1, t2);
    t1 = _mm_srai_epi32(t1, 2);
    x1 = _mm_srai_epi32(x1, 2);
    r2 = _mm_packs_epi32(x1, t1);
    src += srcstride;
    x1 = _mm_loadu_si128((__m128i*)&src[x]);
    x2 = _mm_loadu_si128((__m128i*)&src[x + 1]);
    x3 = _mm_loadu_si128((__m128i*)&src[x + 2 * 1]);
    x4 = _mm_loadu_si128((__m128i*)&src[x + 3 * 1]);
    t1 = _mm_unpackhi_epi16(x1, x2);
    x1 = _mm_unpacklo_epi16(x1, x2);
    t2 = _mm_unpackhi_epi16(x3, x4);
    x2 = _mm_unpacklo_epi16(x3, x4);
    x1 = _mm_madd_epi16(x1, f1);
    t1 = _mm_madd_epi16(t1, f1);
    x2 = _mm_madd_epi16(x2, f2);
    t2 = _mm_madd_epi16(t2, f2);
    x1 = _mm_add_epi32(x1, x2);
    t1 = _mm_add_epi32(t1, t2);
    t1 = _mm_srai_epi32(t1, 2);
    x1 = _mm_srai_epi32(x1, 2);
    r3 = _mm_packs_epi32(x1, t1);
    src += srcstride;
    for (y = 0; y < height; y++) {
      x1 = _mm_loadu_si128((__m128i*)&src[x]);
      x2 = _mm_loadu_si128((__m128i*)&src[x + 1]);
      x3 = _mm_loadu_si128((__m128i*)&src[x + 2 * 1]);
      x4 = _mm_loadu_si128((__m128i*)&src[x + 3 * 1]);
      t1 = _mm_unpackhi_epi16(x1, x2);
      x1 = _mm_unpacklo_epi16(x1, x2);
      t2 = _mm_unpackhi_epi16(x3, x4);
      x2 = _mm_unpacklo_epi16(x3, x4);
      x1 = _mm_madd_epi16(x1, f1);
      t1 = _mm_madd_epi16(t1, f1);
      x2 = _mm_madd_epi16(x2, f2);
      t2 = _mm_madd_epi16(t2, f2);
      x1 = _mm_add_epi32(x1, x2);
      t1 = _mm_add_epi32(t1, t2);
      t1 = _mm_srai_epi32(t1, 2);
      x1 = _mm_srai_epi32(x1, 2);
      r4 = _mm_packs_epi32(x1, t1);
      src += srcstride;
      t1 = _mm_unpackhi_epi16(r1, r2);
      x1 = _mm_unpacklo_epi16(r1, r2);
      t2 = _mm_unpackhi_epi16(r3, r4);
      x2 = _mm_unpacklo_epi16(r3, r4);
      x1 = _mm_madd_epi16(x1, f3);
      t1 = _mm_madd_epi16(t1, f3);
      x2 = _mm_madd_epi16(x2, f4);
      t2 = _mm_madd_epi16(t2, f4);
      x1 = _mm_add_epi32(x1, x2);
      t1 = _mm_add_epi32(t1, t2);
      t1 = _mm_srai_epi32(t1, 6);
      x1 = _mm_srai_epi32(x1, 6);
      x1 = _mm_packs_epi32(x1, t1);
      _mm_storeu_si128((__m128i*)&dst[x], x1);
      dst += 128;
      r1 = r2;
      r2 = r3;
      r3 = r4;
    }
    src = src_bis;
    dst = dst_bis;
  }
}

static void
oh_hevc_put_hevc_bi1_epel_hv8_10_sse(uint16_t* _dst,
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
  __m128i x1, x2, x3, x4, t1, t2, f1, f2, f3, f4, r1, r2, r3, r4;
  const __m128i offset = _mm_set1_epi16(1 << 10);
  uint16_t* dst = (uint16_t*)_dst;
  const int dststride = _dststride;
  const int16_t* src2_bis = src2;
  uint16_t* dst_bis = dst;
  uint16_t *src_bis, *src = ((uint16_t*)_src) - 1;
  ptrdiff_t srcstride = _srcstride;
  src -= 1 * srcstride;
  src_bis = src;
  f1 = _mm_load_si128((__m128i*)oh_hevc_epel_filters_sse[mx - 1][0]);
  f2 = _mm_load_si128((__m128i*)oh_hevc_epel_filters_sse[mx - 1][1]);
  f3 = _mm_load_si128((__m128i*)oh_hevc_epel_filters_sse[my - 1][0]);
  f4 = _mm_load_si128((__m128i*)oh_hevc_epel_filters_sse[my - 1][1]);
  for (x = 0; x < width; x += 8) {
    __m128i r5;
    x1 = _mm_loadu_si128((__m128i*)&src[x]);
    x2 = _mm_loadu_si128((__m128i*)&src[x + 1]);
    x3 = _mm_loadu_si128((__m128i*)&src[x + 2 * 1]);
    x4 = _mm_loadu_si128((__m128i*)&src[x + 3 * 1]);
    t1 = _mm_unpackhi_epi16(x1, x2);
    x1 = _mm_unpacklo_epi16(x1, x2);
    t2 = _mm_unpackhi_epi16(x3, x4);
    x2 = _mm_unpacklo_epi16(x3, x4);
    x1 = _mm_madd_epi16(x1, f1);
    t1 = _mm_madd_epi16(t1, f1);
    x2 = _mm_madd_epi16(x2, f2);
    t2 = _mm_madd_epi16(t2, f2);
    x1 = _mm_add_epi32(x1, x2);
    t1 = _mm_add_epi32(t1, t2);
    t1 = _mm_srai_epi32(t1, 2);
    x1 = _mm_srai_epi32(x1, 2);
    r1 = _mm_packs_epi32(x1, t1);
    src += srcstride;
    x1 = _mm_loadu_si128((__m128i*)&src[x]);
    x2 = _mm_loadu_si128((__m128i*)&src[x + 1]);
    x3 = _mm_loadu_si128((__m128i*)&src[x + 2 * 1]);
    x4 = _mm_loadu_si128((__m128i*)&src[x + 3 * 1]);
    t1 = _mm_unpackhi_epi16(x1, x2);
    x1 = _mm_unpacklo_epi16(x1, x2);
    t2 = _mm_unpackhi_epi16(x3, x4);
    x2 = _mm_unpacklo_epi16(x3, x4);
    x1 = _mm_madd_epi16(x1, f1);
    t1 = _mm_madd_epi16(t1, f1);
    x2 = _mm_madd_epi16(x2, f2);
    t2 = _mm_madd_epi16(t2, f2);
    x1 = _mm_add_epi32(x1, x2);
    t1 = _mm_add_epi32(t1, t2);
    t1 = _mm_srai_epi32(t1, 2);
    x1 = _mm_srai_epi32(x1, 2);
    r2 = _mm_packs_epi32(x1, t1);
    src += srcstride;
    x1 = _mm_loadu_si128((__m128i*)&src[x]);
    x2 = _mm_loadu_si128((__m128i*)&src[x + 1]);
    x3 = _mm_loadu_si128((__m128i*)&src[x + 2 * 1]);
    x4 = _mm_loadu_si128((__m128i*)&src[x + 3 * 1]);
    t1 = _mm_unpackhi_epi16(x1, x2);
    x1 = _mm_unpacklo_epi16(x1, x2);
    t2 = _mm_unpackhi_epi16(x3, x4);
    x2 = _mm_unpacklo_epi16(x3, x4);
    x1 = _mm_madd_epi16(x1, f1);
    t1 = _mm_madd_epi16(t1, f1);
    x2 = _mm_madd_epi16(x2, f2);
    t2 = _mm_madd_epi16(t2, f2);
    x1 = _mm_add_epi32(x1, x2);
    t1 = _mm_add_epi32(t1, t2);
    t1 = _mm_srai_epi32(t1, 2);
    x1 = _mm_srai_epi32(x1, 2);
    r3 = _mm_packs_epi32(x1, t1);
    src += srcstride;
    for (y = 0; y < height; y++) {
      x1 = _mm_loadu_si128((__m128i*)&src[x]);
      x2 = _mm_loadu_si128((__m128i*)&src[x + 1]);
      x3 = _mm_loadu_si128((__m128i*)&src[x + 2 * 1]);
      x4 = _mm_loadu_si128((__m128i*)&src[x + 3 * 1]);
      t1 = _mm_unpackhi_epi16(x1, x2);
      x1 = _mm_unpacklo_epi16(x1, x2);
      t2 = _mm_unpackhi_epi16(x3, x4);
      x2 = _mm_unpacklo_epi16(x3, x4);
      x1 = _mm_madd_epi16(x1, f1);
      t1 = _mm_madd_epi16(t1, f1);
      x2 = _mm_madd_epi16(x2, f2);
      t2 = _mm_madd_epi16(t2, f2);
      x1 = _mm_add_epi32(x1, x2);
      t1 = _mm_add_epi32(t1, t2);
      t1 = _mm_srai_epi32(t1, 2);
      x1 = _mm_srai_epi32(x1, 2);
      r4 = _mm_packs_epi32(x1, t1);
      src += srcstride;
      t1 = _mm_unpackhi_epi16(r1, r2);
      x1 = _mm_unpacklo_epi16(r1, r2);
      t2 = _mm_unpackhi_epi16(r3, r4);
      x2 = _mm_unpacklo_epi16(r3, r4);
      x1 = _mm_madd_epi16(x1, f3);
      t1 = _mm_madd_epi16(t1, f3);
      x2 = _mm_madd_epi16(x2, f4);
      t2 = _mm_madd_epi16(t2, f4);
      x1 = _mm_add_epi32(x1, x2);
      t1 = _mm_add_epi32(t1, t2);
      t1 = _mm_srai_epi32(t1, 6);
      x1 = _mm_srai_epi32(x1, 6);
      x1 = _mm_packs_epi32(x1, t1);
      r5 = _mm_load_si128((__m128i*)&src2[x]);
      x1 = _mm_adds_epi16(x1, r5);
      x1 = _mm_mulhrs_epi16(x1, offset);
      x1 = _mm_max_epi16(x1, _mm_setzero_si128());
      x1 = _mm_min_epi16(x1, _mm_set1_epi16(0x03FF));
      _mm_storeu_si128((__m128i*)&dst[x], x1);
      src2 += 128;
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
oh_hevc_put_hevc_uni_epel_hv8_10_sse(uint16_t* _dst,
                                     ptrdiff_t _dststride,
                                     const uint16_t* _src,
                                     ptrdiff_t _srcstride,
                                     int height,
                                     intptr_t mx,
                                     intptr_t my,
                                     int width)
{
  int x, y;
  __m128i x1, x2, x3, x4, t1, t2, f1, f2, f3, f4, r1, r2, r3, r4;
  const __m128i offset = _mm_set1_epi16(1 << (10 + 1));
  uint16_t* dst = (uint16_t*)_dst;
  const int dststride = _dststride;
  uint16_t* dst_bis = dst;
  uint16_t *src_bis, *src = ((uint16_t*)_src) - 1;
  ptrdiff_t srcstride = _srcstride;
  src -= 1 * srcstride;
  src_bis = src;
  f1 = _mm_load_si128((__m128i*)oh_hevc_epel_filters_sse[mx - 1][0]);
  f2 = _mm_load_si128((__m128i*)oh_hevc_epel_filters_sse[mx - 1][1]);
  f3 = _mm_load_si128((__m128i*)oh_hevc_epel_filters_sse[my - 1][0]);
  f4 = _mm_load_si128((__m128i*)oh_hevc_epel_filters_sse[my - 1][1]);
  for (x = 0; x < width; x += 8) {
    x1 = _mm_loadu_si128((__m128i*)&src[x]);
    x2 = _mm_loadu_si128((__m128i*)&src[x + 1]);
    x3 = _mm_loadu_si128((__m128i*)&src[x + 2 * 1]);
    x4 = _mm_loadu_si128((__m128i*)&src[x + 3 * 1]);
    t1 = _mm_unpackhi_epi16(x1, x2);
    x1 = _mm_unpacklo_epi16(x1, x2);
    t2 = _mm_unpackhi_epi16(x3, x4);
    x2 = _mm_unpacklo_epi16(x3, x4);
    x1 = _mm_madd_epi16(x1, f1);
    t1 = _mm_madd_epi16(t1, f1);
    x2 = _mm_madd_epi16(x2, f2);
    t2 = _mm_madd_epi16(t2, f2);
    x1 = _mm_add_epi32(x1, x2);
    t1 = _mm_add_epi32(t1, t2);
    t1 = _mm_srai_epi32(t1, 2);
    x1 = _mm_srai_epi32(x1, 2);
    r1 = _mm_packs_epi32(x1, t1);
    src += srcstride;
    x1 = _mm_loadu_si128((__m128i*)&src[x]);
    x2 = _mm_loadu_si128((__m128i*)&src[x + 1]);
    x3 = _mm_loadu_si128((__m128i*)&src[x + 2 * 1]);
    x4 = _mm_loadu_si128((__m128i*)&src[x + 3 * 1]);
    t1 = _mm_unpackhi_epi16(x1, x2);
    x1 = _mm_unpacklo_epi16(x1, x2);
    t2 = _mm_unpackhi_epi16(x3, x4);
    x2 = _mm_unpacklo_epi16(x3, x4);
    x1 = _mm_madd_epi16(x1, f1);
    t1 = _mm_madd_epi16(t1, f1);
    x2 = _mm_madd_epi16(x2, f2);
    t2 = _mm_madd_epi16(t2, f2);
    x1 = _mm_add_epi32(x1, x2);
    t1 = _mm_add_epi32(t1, t2);
    t1 = _mm_srai_epi32(t1, 2);
    x1 = _mm_srai_epi32(x1, 2);
    r2 = _mm_packs_epi32(x1, t1);
    src += srcstride;
    x1 = _mm_loadu_si128((__m128i*)&src[x]);
    x2 = _mm_loadu_si128((__m128i*)&src[x + 1]);
    x3 = _mm_loadu_si128((__m128i*)&src[x + 2 * 1]);
    x4 = _mm_loadu_si128((__m128i*)&src[x + 3 * 1]);
    t1 = _mm_unpackhi_epi16(x1, x2);
    x1 = _mm_unpacklo_epi16(x1, x2);
    t2 = _mm_unpackhi_epi16(x3, x4);
    x2 = _mm_unpacklo_epi16(x3, x4);
    x1 = _mm_madd_epi16(x1, f1);
    t1 = _mm_madd_epi16(t1, f1);
    x2 = _mm_madd_epi16(x2, f2);
    t2 = _mm_madd_epi16(t2, f2);
    x1 = _mm_add_epi32(x1, x2);
    t1 = _mm_add_epi32(t1, t2);
    t1 = _mm_srai_epi32(t1, 2);
    x1 = _mm_srai_epi32(x1, 2);
    r3 = _mm_packs_epi32(x1, t1);
    src += srcstride;
    for (y = 0; y < height; y++) {
      x1 = _mm_loadu_si128((__m128i*)&src[x]);
      x2 = _mm_loadu_si128((__m128i*)&src[x + 1]);
      x3 = _mm_loadu_si128((__m128i*)&src[x + 2 * 1]);
      x4 = _mm_loadu_si128((__m128i*)&src[x + 3 * 1]);
      t1 = _mm_unpackhi_epi16(x1, x2);
      x1 = _mm_unpacklo_epi16(x1, x2);
      t2 = _mm_unpackhi_epi16(x3, x4);
      x2 = _mm_unpacklo_epi16(x3, x4);
      x1 = _mm_madd_epi16(x1, f1);
      t1 = _mm_madd_epi16(t1, f1);
      x2 = _mm_madd_epi16(x2, f2);
      t2 = _mm_madd_epi16(t2, f2);
      x1 = _mm_add_epi32(x1, x2);
      t1 = _mm_add_epi32(t1, t2);
      t1 = _mm_srai_epi32(t1, 2);
      x1 = _mm_srai_epi32(x1, 2);
      r4 = _mm_packs_epi32(x1, t1);
      src += srcstride;
      t1 = _mm_unpackhi_epi16(r1, r2);
      x1 = _mm_unpacklo_epi16(r1, r2);
      t2 = _mm_unpackhi_epi16(r3, r4);
      x2 = _mm_unpacklo_epi16(r3, r4);
      x1 = _mm_madd_epi16(x1, f3);
      t1 = _mm_madd_epi16(t1, f3);
      x2 = _mm_madd_epi16(x2, f4);
      t2 = _mm_madd_epi16(t2, f4);
      x1 = _mm_add_epi32(x1, x2);
      t1 = _mm_add_epi32(t1, t2);
      t1 = _mm_srai_epi32(t1, 6);
      x1 = _mm_srai_epi32(x1, 6);
      x1 = _mm_packs_epi32(x1, t1);
      x1 = _mm_mulhrs_epi16(x1, offset);
      x1 = _mm_max_epi16(x1, _mm_setzero_si128());
      x1 = _mm_min_epi16(x1, _mm_set1_epi16(0x03FF));
      _mm_storeu_si128((__m128i*)&dst[x], x1);
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
oh_hevc_put_hevc_bi0_qpel_h4_10_sse(int16_t* dst,
                                    const uint16_t* _src,
                                    ptrdiff_t _srcstride,
                                    int height,
                                    intptr_t mx,
                                    intptr_t my,
                                    int width)
{
  int x, y;
  int shift = 10 - 8;
  const __m128i c0 = _mm_setzero_si128();
  __m128i x1, x2, x3, x4, r1, r3, c1, c2, c3, c4;
  const uint16_t* src = _src;
  const int srcstride = _srcstride;
  c1 = _mm_load_si128((__m128i*)oh_hevc_qpel_filters_sse[mx - 1][0]);
  c2 = _mm_load_si128((__m128i*)oh_hevc_qpel_filters_sse[mx - 1][1]);
  c3 = _mm_load_si128((__m128i*)oh_hevc_qpel_filters_sse[mx - 1][2]);
  c4 = _mm_load_si128((__m128i*)oh_hevc_qpel_filters_sse[mx - 1][3]);
  for (y = 0; y < height; y++) {
    for (x = 0; x < width; x += 4) {
      x1 = _mm_loadl_epi64((__m128i*)&src[x - 3 * 1]);
      x2 = _mm_loadl_epi64((__m128i*)&src[x - 2 * 1]);
      x3 = _mm_loadl_epi64((__m128i*)&src[x - 1 * 1]);
      x4 = _mm_loadl_epi64((__m128i*)&src[x]);
      x1 = _mm_unpacklo_epi16(x1, x2);
      x2 = _mm_unpacklo_epi16(x3, x4);
      x1 = _mm_madd_epi16(x1, c1);
      x2 = _mm_madd_epi16(x2, c2);
      r1 = _mm_add_epi32(x1, x2);
      x1 = _mm_loadl_epi64((__m128i*)&src[x + 1]);
      x2 = _mm_loadl_epi64((__m128i*)&src[x + 2 * 1]);
      x3 = _mm_loadl_epi64((__m128i*)&src[x + 3 * 1]);
      x4 = _mm_loadl_epi64((__m128i*)&src[x + 4 * 1]);
      x1 = _mm_unpacklo_epi16(x1, x2);
      x2 = _mm_unpacklo_epi16(x3, x4);
      x1 = _mm_madd_epi16(x1, c3);
      x2 = _mm_madd_epi16(x2, c4);
      r3 = _mm_add_epi32(x1, x2);
      x1 = _mm_add_epi32(r1, r3);
      x1 = _mm_srai_epi32(x1, shift);
      x1 = _mm_packs_epi32(x1, c0);
      _mm_storel_epi64((__m128i*)&dst[x], x1);
    }
    src += srcstride;
    dst += 128;
  }
}

static void
oh_hevc_put_hevc_bi1_qpel_h4_10_sse(uint16_t* _dst,
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
  int shift = 10 - 8;
  const __m128i c0 = _mm_setzero_si128();
  __m128i x1, x2, x3, x4, r1, r3, c1, c2, c3, c4;
  const uint16_t* src = _src;
  const int srcstride = _srcstride;
  const __m128i offset = _mm_set1_epi16(1 << 10);
  uint16_t* dst = (uint16_t*)_dst;
  const int dststride = _dststride;
  c1 = _mm_load_si128((__m128i*)oh_hevc_qpel_filters_sse[mx - 1][0]);
  c2 = _mm_load_si128((__m128i*)oh_hevc_qpel_filters_sse[mx - 1][1]);
  c3 = _mm_load_si128((__m128i*)oh_hevc_qpel_filters_sse[mx - 1][2]);
  c4 = _mm_load_si128((__m128i*)oh_hevc_qpel_filters_sse[mx - 1][3]);
  for (y = 0; y < height; y++) {
    for (x = 0; x < width; x += 4) {
      __m128i r5;
      x1 = _mm_loadl_epi64((__m128i*)&src[x - 3 * 1]);
      x2 = _mm_loadl_epi64((__m128i*)&src[x - 2 * 1]);
      x3 = _mm_loadl_epi64((__m128i*)&src[x - 1 * 1]);
      x4 = _mm_loadl_epi64((__m128i*)&src[x]);
      x1 = _mm_unpacklo_epi16(x1, x2);
      x2 = _mm_unpacklo_epi16(x3, x4);
      x1 = _mm_madd_epi16(x1, c1);
      x2 = _mm_madd_epi16(x2, c2);
      r1 = _mm_add_epi32(x1, x2);
      x1 = _mm_loadl_epi64((__m128i*)&src[x + 1]);
      x2 = _mm_loadl_epi64((__m128i*)&src[x + 2 * 1]);
      x3 = _mm_loadl_epi64((__m128i*)&src[x + 3 * 1]);
      x4 = _mm_loadl_epi64((__m128i*)&src[x + 4 * 1]);
      x1 = _mm_unpacklo_epi16(x1, x2);
      x2 = _mm_unpacklo_epi16(x3, x4);
      x1 = _mm_madd_epi16(x1, c3);
      x2 = _mm_madd_epi16(x2, c4);
      r3 = _mm_add_epi32(x1, x2);
      x1 = _mm_add_epi32(r1, r3);
      x1 = _mm_srai_epi32(x1, shift);
      x1 = _mm_packs_epi32(x1, c0);
      r5 = _mm_loadl_epi64((__m128i*)&src2[x]);
      x1 = _mm_adds_epi16(x1, r5);
      x1 = _mm_mulhrs_epi16(x1, offset);
      x1 = _mm_max_epi16(x1, _mm_setzero_si128());
      x1 = _mm_min_epi16(x1, _mm_set1_epi16(0x03FF));
      _mm_storel_epi64((__m128i*)&dst[x], x1);
    }
    src += srcstride;
    src2 += 128;
    dst += dststride;
  }
}

static void
oh_hevc_put_hevc_uni_qpel_h4_10_sse(uint16_t* _dst,
                                    ptrdiff_t _dststride,
                                    const uint16_t* _src,
                                    ptrdiff_t _srcstride,
                                    int height,
                                    intptr_t mx,
                                    intptr_t my,
                                    int width)
{
  int x, y;
  int shift = 10 - 8;
  const __m128i c0 = _mm_setzero_si128();
  __m128i x1, x2, x3, x4, r1, r3, c1, c2, c3, c4;
  const uint16_t* src = _src;
  const int srcstride = _srcstride;
  const __m128i offset = _mm_set1_epi16(1 << (10 + 1));
  uint16_t* dst = (uint16_t*)_dst;
  const int dststride = _dststride;
  c1 = _mm_load_si128((__m128i*)oh_hevc_qpel_filters_sse[mx - 1][0]);
  c2 = _mm_load_si128((__m128i*)oh_hevc_qpel_filters_sse[mx - 1][1]);
  c3 = _mm_load_si128((__m128i*)oh_hevc_qpel_filters_sse[mx - 1][2]);
  c4 = _mm_load_si128((__m128i*)oh_hevc_qpel_filters_sse[mx - 1][3]);
  for (y = 0; y < height; y++) {
    for (x = 0; x < width; x += 4) {
      x1 = _mm_loadl_epi64((__m128i*)&src[x - 3 * 1]);
      x2 = _mm_loadl_epi64((__m128i*)&src[x - 2 * 1]);
      x3 = _mm_loadl_epi64((__m128i*)&src[x - 1 * 1]);
      x4 = _mm_loadl_epi64((__m128i*)&src[x]);
      x1 = _mm_unpacklo_epi16(x1, x2);
      x2 = _mm_unpacklo_epi16(x3, x4);
      x1 = _mm_madd_epi16(x1, c1);
      x2 = _mm_madd_epi16(x2, c2);
      r1 = _mm_add_epi32(x1, x2);
      x1 = _mm_loadl_epi64((__m128i*)&src[x + 1]);
      x2 = _mm_loadl_epi64((__m128i*)&src[x + 2 * 1]);
      x3 = _mm_loadl_epi64((__m128i*)&src[x + 3 * 1]);
      x4 = _mm_loadl_epi64((__m128i*)&src[x + 4 * 1]);
      x1 = _mm_unpacklo_epi16(x1, x2);
      x2 = _mm_unpacklo_epi16(x3, x4);
      x1 = _mm_madd_epi16(x1, c3);
      x2 = _mm_madd_epi16(x2, c4);
      r3 = _mm_add_epi32(x1, x2);
      x1 = _mm_add_epi32(r1, r3);
      x1 = _mm_srai_epi32(x1, shift);
      x1 = _mm_packs_epi32(x1, c0);
      x1 = _mm_mulhrs_epi16(x1, offset);
      x1 = _mm_max_epi16(x1, _mm_setzero_si128());
      x1 = _mm_min_epi16(x1, _mm_set1_epi16(0x03FF));
      _mm_storel_epi64((__m128i*)&dst[x], x1);
    }
    src += srcstride;
    dst += dststride;
  }
}

static void
oh_hevc_put_hevc_bi0_qpel_h8_10_sse(int16_t* dst,
                                    const uint16_t* _src,
                                    ptrdiff_t _srcstride,
                                    int height,
                                    intptr_t mx,
                                    intptr_t my,
                                    int width)
{
  int x, y;
  int shift = 10 - 8;
  __m128i x1, x2, x3, x4, r1, r2, r3, r4, c1, c2, c3, c4, t1, t2;
  const uint16_t* src = _src;
  const int srcstride = _srcstride;
  c1 = _mm_load_si128((__m128i*)oh_hevc_qpel_filters_sse[mx - 1][0]);
  c2 = _mm_load_si128((__m128i*)oh_hevc_qpel_filters_sse[mx - 1][1]);
  c3 = _mm_load_si128((__m128i*)oh_hevc_qpel_filters_sse[mx - 1][2]);
  c4 = _mm_load_si128((__m128i*)oh_hevc_qpel_filters_sse[mx - 1][3]);
  for (y = 0; y < height; y++) {
    for (x = 0; x < width; x += 8) {
      x1 = _mm_loadu_si128((__m128i*)&src[x - 3 * 1]);
      x2 = _mm_loadu_si128((__m128i*)&src[x - 2 * 1]);
      x3 = _mm_loadu_si128((__m128i*)&src[x - 1 * 1]);
      x4 = _mm_loadu_si128((__m128i*)&src[x]);
      t1 = _mm_unpackhi_epi16(x1, x2);
      x1 = _mm_unpacklo_epi16(x1, x2);
      t2 = _mm_unpackhi_epi16(x3, x4);
      x2 = _mm_unpacklo_epi16(x3, x4);
      t1 = _mm_madd_epi16(t1, c1);
      t2 = _mm_madd_epi16(t2, c2);
      r2 = _mm_add_epi32(t1, t2);
      x1 = _mm_madd_epi16(x1, c1);
      x2 = _mm_madd_epi16(x2, c2);
      r1 = _mm_add_epi32(x1, x2);
      x1 = _mm_loadu_si128((__m128i*)&src[x + 1]);
      x2 = _mm_loadu_si128((__m128i*)&src[x + 2 * 1]);
      x3 = _mm_loadu_si128((__m128i*)&src[x + 3 * 1]);
      x4 = _mm_loadu_si128((__m128i*)&src[x + 4 * 1]);
      t1 = _mm_unpackhi_epi16(x1, x2);
      x1 = _mm_unpacklo_epi16(x1, x2);
      t2 = _mm_unpackhi_epi16(x3, x4);
      x2 = _mm_unpacklo_epi16(x3, x4);
      t1 = _mm_madd_epi16(t1, c3);
      t2 = _mm_madd_epi16(t2, c4);
      r4 = _mm_add_epi32(t1, t2);
      x1 = _mm_madd_epi16(x1, c3);
      x2 = _mm_madd_epi16(x2, c4);
      r3 = _mm_add_epi32(x1, x2);
      x1 = _mm_add_epi32(r1, r3);
      x2 = _mm_add_epi32(r2, r4);
      x1 = _mm_srai_epi32(x1, shift);
      x2 = _mm_srai_epi32(x2, shift);
      x1 = _mm_packs_epi32(x1, x2);
      _mm_storeu_si128((__m128i*)&dst[x], x1);
    }
    src += srcstride;
    dst += 128;
  }
}

static void
oh_hevc_put_hevc_bi1_qpel_h8_10_sse(uint16_t* _dst,
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
  int shift = 10 - 8;
  __m128i x1, x2, x3, x4, r1, r2, r3, r4, c1, c2, c3, c4, t1, t2;
  const uint16_t* src = _src;
  const int srcstride = _srcstride;
  const __m128i offset = _mm_set1_epi16(1 << 10);
  uint16_t* dst = (uint16_t*)_dst;
  const int dststride = _dststride;
  c1 = _mm_load_si128((__m128i*)oh_hevc_qpel_filters_sse[mx - 1][0]);
  c2 = _mm_load_si128((__m128i*)oh_hevc_qpel_filters_sse[mx - 1][1]);
  c3 = _mm_load_si128((__m128i*)oh_hevc_qpel_filters_sse[mx - 1][2]);
  c4 = _mm_load_si128((__m128i*)oh_hevc_qpel_filters_sse[mx - 1][3]);
  for (y = 0; y < height; y++) {
    for (x = 0; x < width; x += 8) {
      __m128i r5;
      x1 = _mm_loadu_si128((__m128i*)&src[x - 3 * 1]);
      x2 = _mm_loadu_si128((__m128i*)&src[x - 2 * 1]);
      x3 = _mm_loadu_si128((__m128i*)&src[x - 1 * 1]);
      x4 = _mm_loadu_si128((__m128i*)&src[x]);
      t1 = _mm_unpackhi_epi16(x1, x2);
      x1 = _mm_unpacklo_epi16(x1, x2);
      t2 = _mm_unpackhi_epi16(x3, x4);
      x2 = _mm_unpacklo_epi16(x3, x4);
      t1 = _mm_madd_epi16(t1, c1);
      t2 = _mm_madd_epi16(t2, c2);
      r2 = _mm_add_epi32(t1, t2);
      x1 = _mm_madd_epi16(x1, c1);
      x2 = _mm_madd_epi16(x2, c2);
      r1 = _mm_add_epi32(x1, x2);
      x1 = _mm_loadu_si128((__m128i*)&src[x + 1]);
      x2 = _mm_loadu_si128((__m128i*)&src[x + 2 * 1]);
      x3 = _mm_loadu_si128((__m128i*)&src[x + 3 * 1]);
      x4 = _mm_loadu_si128((__m128i*)&src[x + 4 * 1]);
      t1 = _mm_unpackhi_epi16(x1, x2);
      x1 = _mm_unpacklo_epi16(x1, x2);
      t2 = _mm_unpackhi_epi16(x3, x4);
      x2 = _mm_unpacklo_epi16(x3, x4);
      t1 = _mm_madd_epi16(t1, c3);
      t2 = _mm_madd_epi16(t2, c4);
      r4 = _mm_add_epi32(t1, t2);
      x1 = _mm_madd_epi16(x1, c3);
      x2 = _mm_madd_epi16(x2, c4);
      r3 = _mm_add_epi32(x1, x2);
      x1 = _mm_add_epi32(r1, r3);
      x2 = _mm_add_epi32(r2, r4);
      x1 = _mm_srai_epi32(x1, shift);
      x2 = _mm_srai_epi32(x2, shift);
      x1 = _mm_packs_epi32(x1, x2);
      r5 = _mm_load_si128((__m128i*)&src2[x]);
      x1 = _mm_adds_epi16(x1, r5);
      x1 = _mm_mulhrs_epi16(x1, offset);
      x1 = _mm_max_epi16(x1, _mm_setzero_si128());
      x1 = _mm_min_epi16(x1, _mm_set1_epi16(0x03FF));
      _mm_storeu_si128((__m128i*)&dst[x], x1);
    }
    src += srcstride;
    src2 += 128;
    dst += dststride;
  }
}

static void
oh_hevc_put_hevc_uni_qpel_h8_10_sse(uint16_t* _dst,
                                    ptrdiff_t _dststride,
                                    const uint16_t* _src,
                                    ptrdiff_t _srcstride,
                                    int height,
                                    intptr_t mx,
                                    intptr_t my,
                                    int width)
{
  int x, y;
  int shift = 10 - 8;
  __m128i x1, x2, x3, x4, r1, r2, r3, r4, c1, c2, c3, c4, t1, t2;
  const uint16_t* src = _src;
  const int srcstride = _srcstride;
  const __m128i offset = _mm_set1_epi16(1 << (10 + 1));
  uint16_t* dst = (uint16_t*)_dst;
  const int dststride = _dststride;
  c1 = _mm_load_si128((__m128i*)oh_hevc_qpel_filters_sse[mx - 1][0]);
  c2 = _mm_load_si128((__m128i*)oh_hevc_qpel_filters_sse[mx - 1][1]);
  c3 = _mm_load_si128((__m128i*)oh_hevc_qpel_filters_sse[mx - 1][2]);
  c4 = _mm_load_si128((__m128i*)oh_hevc_qpel_filters_sse[mx - 1][3]);
  for (y = 0; y < height; y++) {
    for (x = 0; x < width; x += 8) {
      x1 = _mm_loadu_si128((__m128i*)&src[x - 3 * 1]);
      x2 = _mm_loadu_si128((__m128i*)&src[x - 2 * 1]);
      x3 = _mm_loadu_si128((__m128i*)&src[x - 1 * 1]);
      x4 = _mm_loadu_si128((__m128i*)&src[x]);
      t1 = _mm_unpackhi_epi16(x1, x2);
      x1 = _mm_unpacklo_epi16(x1, x2);
      t2 = _mm_unpackhi_epi16(x3, x4);
      x2 = _mm_unpacklo_epi16(x3, x4);
      t1 = _mm_madd_epi16(t1, c1);
      t2 = _mm_madd_epi16(t2, c2);
      r2 = _mm_add_epi32(t1, t2);
      x1 = _mm_madd_epi16(x1, c1);
      x2 = _mm_madd_epi16(x2, c2);
      r1 = _mm_add_epi32(x1, x2);
      x1 = _mm_loadu_si128((__m128i*)&src[x + 1]);
      x2 = _mm_loadu_si128((__m128i*)&src[x + 2 * 1]);
      x3 = _mm_loadu_si128((__m128i*)&src[x + 3 * 1]);
      x4 = _mm_loadu_si128((__m128i*)&src[x + 4 * 1]);
      t1 = _mm_unpackhi_epi16(x1, x2);
      x1 = _mm_unpacklo_epi16(x1, x2);
      t2 = _mm_unpackhi_epi16(x3, x4);
      x2 = _mm_unpacklo_epi16(x3, x4);
      t1 = _mm_madd_epi16(t1, c3);
      t2 = _mm_madd_epi16(t2, c4);
      r4 = _mm_add_epi32(t1, t2);
      x1 = _mm_madd_epi16(x1, c3);
      x2 = _mm_madd_epi16(x2, c4);
      r3 = _mm_add_epi32(x1, x2);
      x1 = _mm_add_epi32(r1, r3);
      x2 = _mm_add_epi32(r2, r4);
      x1 = _mm_srai_epi32(x1, shift);
      x2 = _mm_srai_epi32(x2, shift);
      x1 = _mm_packs_epi32(x1, x2);
      x1 = _mm_mulhrs_epi16(x1, offset);
      x1 = _mm_max_epi16(x1, _mm_setzero_si128());
      x1 = _mm_min_epi16(x1, _mm_set1_epi16(0x03FF));
      _mm_storeu_si128((__m128i*)&dst[x], x1);
    }
    src += srcstride;
    dst += dststride;
  }
}

static void
oh_hevc_put_hevc_bi0_qpel_v4_10_sse(int16_t* dst,
                                    const uint16_t* _src,
                                    ptrdiff_t _srcstride,
                                    int height,
                                    intptr_t mx,
                                    intptr_t my,
                                    int width)
{
  int x, y;
  const __m128i c0 = _mm_setzero_si128();
  __m128i x1, x2, x3, x4, x5, x6, x7, x8, c1, c2, c3, c4;
  const uint16_t* src = _src;
  const int srcstride = _srcstride;
  c1 = _mm_load_si128((__m128i*)oh_hevc_qpel_filters_sse[my - 1][0]);
  c2 = _mm_load_si128((__m128i*)oh_hevc_qpel_filters_sse[my - 1][1]);
  c3 = _mm_load_si128((__m128i*)oh_hevc_qpel_filters_sse[my - 1][2]);
  c4 = _mm_load_si128((__m128i*)oh_hevc_qpel_filters_sse[my - 1][3]);
  for (y = 0; y < height; y++) {
    for (x = 0; x < width; x += 4) {
      x1 = _mm_loadu_si128((__m128i*)&src[x - 3 * srcstride]);
      x2 = _mm_loadu_si128((__m128i*)&src[x - 2 * srcstride]);
      x3 = _mm_loadu_si128((__m128i*)&src[x - 1 * srcstride]);
      x4 = _mm_loadu_si128((__m128i*)&src[x]);
      x5 = _mm_loadu_si128((__m128i*)&src[x + srcstride]);
      x6 = _mm_loadu_si128((__m128i*)&src[x + 2 * srcstride]);
      x7 = _mm_loadu_si128((__m128i*)&src[x + 3 * srcstride]);
      x8 = _mm_loadu_si128((__m128i*)&src[x + 4 * srcstride]);
      x1 = _mm_unpacklo_epi16(x1, x2);
      x2 = _mm_unpacklo_epi16(x3, x4);
      x3 = _mm_unpacklo_epi16(x5, x6);
      x4 = _mm_unpacklo_epi16(x7, x8);
      x2 = _mm_madd_epi16(x2, c2);
      x3 = _mm_madd_epi16(x3, c3);
      x1 = _mm_madd_epi16(x1, c1);
      x4 = _mm_madd_epi16(x4, c4);
      x1 = _mm_add_epi32(x1, x2);
      x2 = _mm_add_epi32(x3, x4);
      x1 = _mm_add_epi32(x1, x2);
      x1 = _mm_srai_epi32(x1, 2);
      x1 = _mm_packs_epi32(x1, c0);
      _mm_storel_epi64((__m128i*)&dst[x], x1);
    }
    src += srcstride;
    dst += 128;
  }
}

static void
oh_hevc_put_hevc_bi1_qpel_v4_10_sse(uint16_t* _dst,
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
  const __m128i c0 = _mm_setzero_si128();
  __m128i x1, x2, x3, x4, x5, x6, x7, x8, c1, c2, c3, c4;
  const uint16_t* src = _src;
  const int srcstride = _srcstride;
  const __m128i offset = _mm_set1_epi16(1 << 10);
  uint16_t* dst = (uint16_t*)_dst;
  const int dststride = _dststride;
  c1 = _mm_load_si128((__m128i*)oh_hevc_qpel_filters_sse[my - 1][0]);
  c2 = _mm_load_si128((__m128i*)oh_hevc_qpel_filters_sse[my - 1][1]);
  c3 = _mm_load_si128((__m128i*)oh_hevc_qpel_filters_sse[my - 1][2]);
  c4 = _mm_load_si128((__m128i*)oh_hevc_qpel_filters_sse[my - 1][3]);
  for (y = 0; y < height; y++) {
    for (x = 0; x < width; x += 4) {
      __m128i r5;
      x1 = _mm_loadu_si128((__m128i*)&src[x - 3 * srcstride]);
      x2 = _mm_loadu_si128((__m128i*)&src[x - 2 * srcstride]);
      x3 = _mm_loadu_si128((__m128i*)&src[x - 1 * srcstride]);
      x4 = _mm_loadu_si128((__m128i*)&src[x]);
      x5 = _mm_loadu_si128((__m128i*)&src[x + srcstride]);
      x6 = _mm_loadu_si128((__m128i*)&src[x + 2 * srcstride]);
      x7 = _mm_loadu_si128((__m128i*)&src[x + 3 * srcstride]);
      x8 = _mm_loadu_si128((__m128i*)&src[x + 4 * srcstride]);
      x1 = _mm_unpacklo_epi16(x1, x2);
      x2 = _mm_unpacklo_epi16(x3, x4);
      x3 = _mm_unpacklo_epi16(x5, x6);
      x4 = _mm_unpacklo_epi16(x7, x8);
      x2 = _mm_madd_epi16(x2, c2);
      x3 = _mm_madd_epi16(x3, c3);
      x1 = _mm_madd_epi16(x1, c1);
      x4 = _mm_madd_epi16(x4, c4);
      x1 = _mm_add_epi32(x1, x2);
      x2 = _mm_add_epi32(x3, x4);
      x1 = _mm_add_epi32(x1, x2);
      x1 = _mm_srai_epi32(x1, 2);
      x1 = _mm_packs_epi32(x1, c0);
      r5 = _mm_loadl_epi64((__m128i*)&src2[x]);
      x1 = _mm_adds_epi16(x1, r5);
      x1 = _mm_mulhrs_epi16(x1, offset);
      x1 = _mm_max_epi16(x1, _mm_setzero_si128());
      x1 = _mm_min_epi16(x1, _mm_set1_epi16(0x03FF));
      _mm_storel_epi64((__m128i*)&dst[x], x1);
    }
    src += srcstride;
    src2 += 128;
    dst += dststride;
  }
}

static void
oh_hevc_put_hevc_uni_qpel_v4_10_sse(uint16_t* _dst,
                                    ptrdiff_t _dststride,
                                    const uint16_t* _src,
                                    ptrdiff_t _srcstride,
                                    int height,
                                    intptr_t mx,
                                    intptr_t my,
                                    int width)
{
  int x, y;
  const __m128i c0 = _mm_setzero_si128();
  __m128i x1, x2, x3, x4, x5, x6, x7, x8, c1, c2, c3, c4;
  const uint16_t* src = _src;
  const int srcstride = _srcstride;
  const __m128i offset = _mm_set1_epi16(1 << (10 + 1));
  uint16_t* dst = (uint16_t*)_dst;
  const int dststride = _dststride;
  c1 = _mm_load_si128((__m128i*)oh_hevc_qpel_filters_sse[my - 1][0]);
  c2 = _mm_load_si128((__m128i*)oh_hevc_qpel_filters_sse[my - 1][1]);
  c3 = _mm_load_si128((__m128i*)oh_hevc_qpel_filters_sse[my - 1][2]);
  c4 = _mm_load_si128((__m128i*)oh_hevc_qpel_filters_sse[my - 1][3]);
  for (y = 0; y < height; y++) {
    for (x = 0; x < width; x += 4) {
      x1 = _mm_loadu_si128((__m128i*)&src[x - 3 * srcstride]);
      x2 = _mm_loadu_si128((__m128i*)&src[x - 2 * srcstride]);
      x3 = _mm_loadu_si128((__m128i*)&src[x - 1 * srcstride]);
      x4 = _mm_loadu_si128((__m128i*)&src[x]);
      x5 = _mm_loadu_si128((__m128i*)&src[x + srcstride]);
      x6 = _mm_loadu_si128((__m128i*)&src[x + 2 * srcstride]);
      x7 = _mm_loadu_si128((__m128i*)&src[x + 3 * srcstride]);
      x8 = _mm_loadu_si128((__m128i*)&src[x + 4 * srcstride]);
      x1 = _mm_unpacklo_epi16(x1, x2);
      x2 = _mm_unpacklo_epi16(x3, x4);
      x3 = _mm_unpacklo_epi16(x5, x6);
      x4 = _mm_unpacklo_epi16(x7, x8);
      x2 = _mm_madd_epi16(x2, c2);
      x3 = _mm_madd_epi16(x3, c3);
      x1 = _mm_madd_epi16(x1, c1);
      x4 = _mm_madd_epi16(x4, c4);
      x1 = _mm_add_epi32(x1, x2);
      x2 = _mm_add_epi32(x3, x4);
      x1 = _mm_add_epi32(x1, x2);
      x1 = _mm_srai_epi32(x1, 2);
      x1 = _mm_packs_epi32(x1, c0);
      x1 = _mm_mulhrs_epi16(x1, offset);
      x1 = _mm_max_epi16(x1, _mm_setzero_si128());
      x1 = _mm_min_epi16(x1, _mm_set1_epi16(0x03FF));
      _mm_storel_epi64((__m128i*)&dst[x], x1);
    }
    src += srcstride;
    dst += dststride;
  }
}

static void
oh_hevc_put_hevc_bi0_qpel_v8_10_sse(int16_t* dst,
                                    const uint16_t* _src,
                                    ptrdiff_t _srcstride,
                                    int height,
                                    intptr_t mx,
                                    intptr_t my,
                                    int width)
{
  int x, y;
  __m128i x1, x2, x3, x4, x5, x6, x7, x8, x9, c1, c2, c3, c4;
  const uint16_t* src = _src;
  const int srcstride = _srcstride;
  c1 = _mm_load_si128((__m128i*)oh_hevc_qpel_filters_sse[my - 1][0]);
  c2 = _mm_load_si128((__m128i*)oh_hevc_qpel_filters_sse[my - 1][1]);
  c3 = _mm_load_si128((__m128i*)oh_hevc_qpel_filters_sse[my - 1][2]);
  c4 = _mm_load_si128((__m128i*)oh_hevc_qpel_filters_sse[my - 1][3]);
  for (y = 0; y < height; y++) {
    for (x = 0; x < width; x += 8) {
      x1 = _mm_loadu_si128((__m128i*)&src[x - 3 * srcstride]);
      x2 = _mm_loadu_si128((__m128i*)&src[x - 2 * srcstride]);
      x3 = _mm_loadu_si128((__m128i*)&src[x - 1 * srcstride]);
      x4 = _mm_loadu_si128((__m128i*)&src[x]);
      x5 = _mm_loadu_si128((__m128i*)&src[x + srcstride]);
      x6 = _mm_loadu_si128((__m128i*)&src[x + 2 * srcstride]);
      x7 = _mm_loadu_si128((__m128i*)&src[x + 3 * srcstride]);
      x8 = _mm_loadu_si128((__m128i*)&src[x + 4 * srcstride]);
      x9 = x1;
      x1 = _mm_unpacklo_epi16(x9, x2);
      x2 = _mm_unpackhi_epi16(x9, x2);
      x9 = x3;
      x3 = _mm_unpacklo_epi16(x9, x4);
      x4 = _mm_unpackhi_epi16(x9, x4);
      x9 = x5;
      x5 = _mm_unpacklo_epi16(x9, x6);
      x6 = _mm_unpackhi_epi16(x9, x6);
      x9 = x7;
      x7 = _mm_unpacklo_epi16(x9, x8);
      x8 = _mm_unpackhi_epi16(x9, x8);
      x1 = _mm_madd_epi16(x1, c1);
      x3 = _mm_madd_epi16(x3, c2);
      x5 = _mm_madd_epi16(x5, c3);
      x7 = _mm_madd_epi16(x7, c4);
      x2 = _mm_madd_epi16(x2, c1);
      x4 = _mm_madd_epi16(x4, c2);
      x6 = _mm_madd_epi16(x6, c3);
      x8 = _mm_madd_epi16(x8, c4);
      x1 = _mm_add_epi32(x1, x3);
      x3 = _mm_add_epi32(x5, x7);
      x2 = _mm_add_epi32(x2, x4);
      x4 = _mm_add_epi32(x6, x8);
      x1 = _mm_add_epi32(x1, x3);
      x2 = _mm_add_epi32(x2, x4);
      x1 = _mm_srai_epi32(x1, 2);
      x2 = _mm_srai_epi32(x2, 2);
      x1 = _mm_packs_epi32(x1, x2);
      _mm_storeu_si128((__m128i*)&dst[x], x1);
    }
    src += srcstride;
    dst += 128;
  }
}

static void
oh_hevc_put_hevc_bi1_qpel_v8_10_sse(uint16_t* _dst,
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
  __m128i x1, x2, x3, x4, x5, x6, x7, x8, x9, c1, c2, c3, c4;
  const uint16_t* src = _src;
  const int srcstride = _srcstride;
  const __m128i offset = _mm_set1_epi16(1 << 10);
  uint16_t* dst = (uint16_t*)_dst;
  const int dststride = _dststride;
  c1 = _mm_load_si128((__m128i*)oh_hevc_qpel_filters_sse[my - 1][0]);
  c2 = _mm_load_si128((__m128i*)oh_hevc_qpel_filters_sse[my - 1][1]);
  c3 = _mm_load_si128((__m128i*)oh_hevc_qpel_filters_sse[my - 1][2]);
  c4 = _mm_load_si128((__m128i*)oh_hevc_qpel_filters_sse[my - 1][3]);
  for (y = 0; y < height; y++) {
    for (x = 0; x < width; x += 8) {
      __m128i r5;
      x1 = _mm_loadu_si128((__m128i*)&src[x - 3 * srcstride]);
      x2 = _mm_loadu_si128((__m128i*)&src[x - 2 * srcstride]);
      x3 = _mm_loadu_si128((__m128i*)&src[x - 1 * srcstride]);
      x4 = _mm_loadu_si128((__m128i*)&src[x]);
      x5 = _mm_loadu_si128((__m128i*)&src[x + srcstride]);
      x6 = _mm_loadu_si128((__m128i*)&src[x + 2 * srcstride]);
      x7 = _mm_loadu_si128((__m128i*)&src[x + 3 * srcstride]);
      x8 = _mm_loadu_si128((__m128i*)&src[x + 4 * srcstride]);
      x9 = x1;
      x1 = _mm_unpacklo_epi16(x9, x2);
      x2 = _mm_unpackhi_epi16(x9, x2);
      x9 = x3;
      x3 = _mm_unpacklo_epi16(x9, x4);
      x4 = _mm_unpackhi_epi16(x9, x4);
      x9 = x5;
      x5 = _mm_unpacklo_epi16(x9, x6);
      x6 = _mm_unpackhi_epi16(x9, x6);
      x9 = x7;
      x7 = _mm_unpacklo_epi16(x9, x8);
      x8 = _mm_unpackhi_epi16(x9, x8);
      x1 = _mm_madd_epi16(x1, c1);
      x3 = _mm_madd_epi16(x3, c2);
      x5 = _mm_madd_epi16(x5, c3);
      x7 = _mm_madd_epi16(x7, c4);
      x2 = _mm_madd_epi16(x2, c1);
      x4 = _mm_madd_epi16(x4, c2);
      x6 = _mm_madd_epi16(x6, c3);
      x8 = _mm_madd_epi16(x8, c4);
      x1 = _mm_add_epi32(x1, x3);
      x3 = _mm_add_epi32(x5, x7);
      x2 = _mm_add_epi32(x2, x4);
      x4 = _mm_add_epi32(x6, x8);
      x1 = _mm_add_epi32(x1, x3);
      x2 = _mm_add_epi32(x2, x4);
      x1 = _mm_srai_epi32(x1, 2);
      x2 = _mm_srai_epi32(x2, 2);
      x1 = _mm_packs_epi32(x1, x2);
      r5 = _mm_load_si128((__m128i*)&src2[x]);
      x1 = _mm_adds_epi16(x1, r5);
      x1 = _mm_mulhrs_epi16(x1, offset);
      x1 = _mm_max_epi16(x1, _mm_setzero_si128());
      x1 = _mm_min_epi16(x1, _mm_set1_epi16(0x03FF));
      _mm_storeu_si128((__m128i*)&dst[x], x1);
    }
    src += srcstride;
    src2 += 128;
    dst += dststride;
  }
}

static void
oh_hevc_put_hevc_uni_qpel_v8_10_sse(uint16_t* _dst,
                                    ptrdiff_t _dststride,
                                    const uint16_t* _src,
                                    ptrdiff_t _srcstride,
                                    int height,
                                    intptr_t mx,
                                    intptr_t my,
                                    int width)
{
  int x, y;
  __m128i x1, x2, x3, x4, x5, x6, x7, x8, x9, c1, c2, c3, c4;
  const uint16_t* src = _src;
  const int srcstride = _srcstride;
  const __m128i offset = _mm_set1_epi16(1 << (10 + 1));
  uint16_t* dst = (uint16_t*)_dst;
  const int dststride = _dststride;
  c1 = _mm_load_si128((__m128i*)oh_hevc_qpel_filters_sse[my - 1][0]);
  c2 = _mm_load_si128((__m128i*)oh_hevc_qpel_filters_sse[my - 1][1]);
  c3 = _mm_load_si128((__m128i*)oh_hevc_qpel_filters_sse[my - 1][2]);
  c4 = _mm_load_si128((__m128i*)oh_hevc_qpel_filters_sse[my - 1][3]);
  for (y = 0; y < height; y++) {
    for (x = 0; x < width; x += 8) {
      x1 = _mm_loadu_si128((__m128i*)&src[x - 3 * srcstride]);
      x2 = _mm_loadu_si128((__m128i*)&src[x - 2 * srcstride]);
      x3 = _mm_loadu_si128((__m128i*)&src[x - 1 * srcstride]);
      x4 = _mm_loadu_si128((__m128i*)&src[x]);
      x5 = _mm_loadu_si128((__m128i*)&src[x + srcstride]);
      x6 = _mm_loadu_si128((__m128i*)&src[x + 2 * srcstride]);
      x7 = _mm_loadu_si128((__m128i*)&src[x + 3 * srcstride]);
      x8 = _mm_loadu_si128((__m128i*)&src[x + 4 * srcstride]);
      x9 = x1;
      x1 = _mm_unpacklo_epi16(x9, x2);
      x2 = _mm_unpackhi_epi16(x9, x2);
      x9 = x3;
      x3 = _mm_unpacklo_epi16(x9, x4);
      x4 = _mm_unpackhi_epi16(x9, x4);
      x9 = x5;
      x5 = _mm_unpacklo_epi16(x9, x6);
      x6 = _mm_unpackhi_epi16(x9, x6);
      x9 = x7;
      x7 = _mm_unpacklo_epi16(x9, x8);
      x8 = _mm_unpackhi_epi16(x9, x8);
      x1 = _mm_madd_epi16(x1, c1);
      x3 = _mm_madd_epi16(x3, c2);
      x5 = _mm_madd_epi16(x5, c3);
      x7 = _mm_madd_epi16(x7, c4);
      x2 = _mm_madd_epi16(x2, c1);
      x4 = _mm_madd_epi16(x4, c2);
      x6 = _mm_madd_epi16(x6, c3);
      x8 = _mm_madd_epi16(x8, c4);
      x1 = _mm_add_epi32(x1, x3);
      x3 = _mm_add_epi32(x5, x7);
      x2 = _mm_add_epi32(x2, x4);
      x4 = _mm_add_epi32(x6, x8);
      x1 = _mm_add_epi32(x1, x3);
      x2 = _mm_add_epi32(x2, x4);
      x1 = _mm_srai_epi32(x1, 2);
      x2 = _mm_srai_epi32(x2, 2);
      x1 = _mm_packs_epi32(x1, x2);
      x1 = _mm_mulhrs_epi16(x1, offset);
      x1 = _mm_max_epi16(x1, _mm_setzero_si128());
      x1 = _mm_min_epi16(x1, _mm_set1_epi16(0x03FF));
      _mm_storeu_si128((__m128i*)&dst[x], x1);
    }
    src += srcstride;
    dst += dststride;
  }
}

static void
oh_hevc_put_hevc_bi0_qpel_v4_14_sse(int16_t* dst,
                                    const uint16_t* _src,
                                    ptrdiff_t _srcstride,
                                    int height,
                                    intptr_t mx,
                                    intptr_t my,
                                    int width)
{
  int x, y;
  const __m128i c0 = _mm_setzero_si128();
  __m128i x1, x2, x3, x4, x5, x6, x7, x8, c1, c2, c3, c4;
  const uint16_t* src = _src;
  const int srcstride = _srcstride;
  c1 = _mm_load_si128((__m128i*)oh_hevc_qpel_filters_sse[my - 1][0]);
  c2 = _mm_load_si128((__m128i*)oh_hevc_qpel_filters_sse[my - 1][1]);
  c3 = _mm_load_si128((__m128i*)oh_hevc_qpel_filters_sse[my - 1][2]);
  c4 = _mm_load_si128((__m128i*)oh_hevc_qpel_filters_sse[my - 1][3]);
  for (y = 0; y < height; y++) {
    for (x = 0; x < width; x += 4) {
      x1 = _mm_loadu_si128((__m128i*)&src[x - 3 * srcstride]);
      x2 = _mm_loadu_si128((__m128i*)&src[x - 2 * srcstride]);
      x3 = _mm_loadu_si128((__m128i*)&src[x - 1 * srcstride]);
      x4 = _mm_loadu_si128((__m128i*)&src[x]);
      x5 = _mm_loadu_si128((__m128i*)&src[x + srcstride]);
      x6 = _mm_loadu_si128((__m128i*)&src[x + 2 * srcstride]);
      x7 = _mm_loadu_si128((__m128i*)&src[x + 3 * srcstride]);
      x8 = _mm_loadu_si128((__m128i*)&src[x + 4 * srcstride]);
      x1 = _mm_unpacklo_epi16(x1, x2);
      x2 = _mm_unpacklo_epi16(x3, x4);
      x3 = _mm_unpacklo_epi16(x5, x6);
      x4 = _mm_unpacklo_epi16(x7, x8);
      x2 = _mm_madd_epi16(x2, c2);
      x3 = _mm_madd_epi16(x3, c3);
      x1 = _mm_madd_epi16(x1, c1);
      x4 = _mm_madd_epi16(x4, c4);
      x1 = _mm_add_epi32(x1, x2);
      x2 = _mm_add_epi32(x3, x4);
      x1 = _mm_add_epi32(x1, x2);
      x1 = _mm_srai_epi32(x1, 6);
      x1 = _mm_packs_epi32(x1, c0);
      _mm_storel_epi64((__m128i*)&dst[x], x1);
    }
    src += srcstride;
    dst += 128;
  }
}

static void
oh_hevc_put_hevc_bi1_qpel_v4_14_10_sse(uint16_t* _dst,
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
  const __m128i c0 = _mm_setzero_si128();
  __m128i x1, x2, x3, x4, x5, x6, x7, x8, c1, c2, c3, c4;
  const int16_t* src = _src;
  const int srcstride = _srcstride;
  const __m128i offset = _mm_set1_epi16(1 << 10);
  uint16_t* dst = (uint16_t*)_dst;
  const int dststride = _dststride;
  c1 = _mm_load_si128((__m128i*)oh_hevc_qpel_filters_sse[my - 1][0]);
  c2 = _mm_load_si128((__m128i*)oh_hevc_qpel_filters_sse[my - 1][1]);
  c3 = _mm_load_si128((__m128i*)oh_hevc_qpel_filters_sse[my - 1][2]);
  c4 = _mm_load_si128((__m128i*)oh_hevc_qpel_filters_sse[my - 1][3]);
  for (y = 0; y < height; y++) {
    for (x = 0; x < width; x += 4) {
      __m128i r5;
      x1 = _mm_loadu_si128((__m128i*)&src[x - 3 * srcstride]);
      x2 = _mm_loadu_si128((__m128i*)&src[x - 2 * srcstride]);
      x3 = _mm_loadu_si128((__m128i*)&src[x - 1 * srcstride]);
      x4 = _mm_loadu_si128((__m128i*)&src[x]);
      x5 = _mm_loadu_si128((__m128i*)&src[x + srcstride]);
      x6 = _mm_loadu_si128((__m128i*)&src[x + 2 * srcstride]);
      x7 = _mm_loadu_si128((__m128i*)&src[x + 3 * srcstride]);
      x8 = _mm_loadu_si128((__m128i*)&src[x + 4 * srcstride]);
      x1 = _mm_unpacklo_epi16(x1, x2);
      x2 = _mm_unpacklo_epi16(x3, x4);
      x3 = _mm_unpacklo_epi16(x5, x6);
      x4 = _mm_unpacklo_epi16(x7, x8);
      x2 = _mm_madd_epi16(x2, c2);
      x3 = _mm_madd_epi16(x3, c3);
      x1 = _mm_madd_epi16(x1, c1);
      x4 = _mm_madd_epi16(x4, c4);
      x1 = _mm_add_epi32(x1, x2);
      x2 = _mm_add_epi32(x3, x4);
      x1 = _mm_add_epi32(x1, x2);
      x1 = _mm_srai_epi32(x1, 6);
      x1 = _mm_packs_epi32(x1, c0);
      r5 = _mm_loadl_epi64((__m128i*)&src2[x]);
      x1 = _mm_adds_epi16(x1, r5);
      x1 = _mm_mulhrs_epi16(x1, offset);
      x1 = _mm_max_epi16(x1, _mm_setzero_si128());
      x1 = _mm_min_epi16(x1, _mm_set1_epi16(0x03FF));
      _mm_storel_epi64((__m128i*)&dst[x], x1);
    }
    src += srcstride;
    src2 += 128;
    dst += dststride;
  }
}

static void
oh_hevc_put_hevc_uni_qpel_v4_14_10_sse(uint16_t* _dst,
                                       ptrdiff_t _dststride,
                                       const int16_t* _src,
                                       ptrdiff_t _srcstride,
                                       int height,
                                       intptr_t mx,
                                       intptr_t my,
                                       int width)
{
  int x, y;
  const __m128i c0 = _mm_setzero_si128();
  __m128i x1, x2, x3, x4, x5, x6, x7, x8, c1, c2, c3, c4;
  const int16_t* src = _src;
  const int srcstride = _srcstride;
  const __m128i offset = _mm_set1_epi16(1 << (10 + 1));
  uint16_t* dst = (uint16_t*)_dst;
  const int dststride = _dststride;
  c1 = _mm_load_si128((__m128i*)oh_hevc_qpel_filters_sse[my - 1][0]);
  c2 = _mm_load_si128((__m128i*)oh_hevc_qpel_filters_sse[my - 1][1]);
  c3 = _mm_load_si128((__m128i*)oh_hevc_qpel_filters_sse[my - 1][2]);
  c4 = _mm_load_si128((__m128i*)oh_hevc_qpel_filters_sse[my - 1][3]);
  for (y = 0; y < height; y++) {
    for (x = 0; x < width; x += 4) {
      x1 = _mm_loadu_si128((__m128i*)&src[x - 3 * srcstride]);
      x2 = _mm_loadu_si128((__m128i*)&src[x - 2 * srcstride]);
      x3 = _mm_loadu_si128((__m128i*)&src[x - 1 * srcstride]);
      x4 = _mm_loadu_si128((__m128i*)&src[x]);
      x5 = _mm_loadu_si128((__m128i*)&src[x + srcstride]);
      x6 = _mm_loadu_si128((__m128i*)&src[x + 2 * srcstride]);
      x7 = _mm_loadu_si128((__m128i*)&src[x + 3 * srcstride]);
      x8 = _mm_loadu_si128((__m128i*)&src[x + 4 * srcstride]);
      x1 = _mm_unpacklo_epi16(x1, x2);
      x2 = _mm_unpacklo_epi16(x3, x4);
      x3 = _mm_unpacklo_epi16(x5, x6);
      x4 = _mm_unpacklo_epi16(x7, x8);
      x2 = _mm_madd_epi16(x2, c2);
      x3 = _mm_madd_epi16(x3, c3);
      x1 = _mm_madd_epi16(x1, c1);
      x4 = _mm_madd_epi16(x4, c4);
      x1 = _mm_add_epi32(x1, x2);
      x2 = _mm_add_epi32(x3, x4);
      x1 = _mm_add_epi32(x1, x2);
      x1 = _mm_srai_epi32(x1, 6);
      x1 = _mm_packs_epi32(x1, c0);
      x1 = _mm_mulhrs_epi16(x1, offset);
      x1 = _mm_max_epi16(x1, _mm_setzero_si128());
      x1 = _mm_min_epi16(x1, _mm_set1_epi16(0x03FF));
      _mm_storel_epi64((__m128i*)&dst[x], x1);
    }
    src += srcstride;
    dst += dststride;
  }
}

static void
oh_hevc_put_hevc_bi0_qpel_v8_14_sse(int16_t* dst,
                                    const uint16_t* _src,
                                    ptrdiff_t _srcstride,
                                    int height,
                                    intptr_t mx,
                                    intptr_t my,
                                    int width)
{
  int x, y;
  __m128i x1, x2, x3, x4, x5, x6, x7, x8, x9, c1, c2, c3, c4;
  const uint16_t* src = _src;
  const int srcstride = _srcstride;
  c1 = _mm_load_si128((__m128i*)oh_hevc_qpel_filters_sse[my - 1][0]);
  c2 = _mm_load_si128((__m128i*)oh_hevc_qpel_filters_sse[my - 1][1]);
  c3 = _mm_load_si128((__m128i*)oh_hevc_qpel_filters_sse[my - 1][2]);
  c4 = _mm_load_si128((__m128i*)oh_hevc_qpel_filters_sse[my - 1][3]);
  for (y = 0; y < height; y++) {
    for (x = 0; x < width; x += 8) {
      x1 = _mm_loadu_si128((__m128i*)&src[x - 3 * srcstride]);
      x2 = _mm_loadu_si128((__m128i*)&src[x - 2 * srcstride]);
      x3 = _mm_loadu_si128((__m128i*)&src[x - 1 * srcstride]);
      x4 = _mm_loadu_si128((__m128i*)&src[x]);
      x5 = _mm_loadu_si128((__m128i*)&src[x + srcstride]);
      x6 = _mm_loadu_si128((__m128i*)&src[x + 2 * srcstride]);
      x7 = _mm_loadu_si128((__m128i*)&src[x + 3 * srcstride]);
      x8 = _mm_loadu_si128((__m128i*)&src[x + 4 * srcstride]);
      x9 = x1;
      x1 = _mm_unpacklo_epi16(x9, x2);
      x2 = _mm_unpackhi_epi16(x9, x2);
      x9 = x3;
      x3 = _mm_unpacklo_epi16(x9, x4);
      x4 = _mm_unpackhi_epi16(x9, x4);
      x9 = x5;
      x5 = _mm_unpacklo_epi16(x9, x6);
      x6 = _mm_unpackhi_epi16(x9, x6);
      x9 = x7;
      x7 = _mm_unpacklo_epi16(x9, x8);
      x8 = _mm_unpackhi_epi16(x9, x8);
      x1 = _mm_madd_epi16(x1, c1);
      x3 = _mm_madd_epi16(x3, c2);
      x5 = _mm_madd_epi16(x5, c3);
      x7 = _mm_madd_epi16(x7, c4);
      x2 = _mm_madd_epi16(x2, c1);
      x4 = _mm_madd_epi16(x4, c2);
      x6 = _mm_madd_epi16(x6, c3);
      x8 = _mm_madd_epi16(x8, c4);
      x1 = _mm_add_epi32(x1, x3);
      x3 = _mm_add_epi32(x5, x7);
      x2 = _mm_add_epi32(x2, x4);
      x4 = _mm_add_epi32(x6, x8);
      x1 = _mm_add_epi32(x1, x3);
      x2 = _mm_add_epi32(x2, x4);
      x1 = _mm_srai_epi32(x1, 6);
      x2 = _mm_srai_epi32(x2, 6);
      x1 = _mm_packs_epi32(x1, x2);
      _mm_storeu_si128((__m128i*)&dst[x], x1);
    }
    src += srcstride;
    dst += 128;
  }
}

static void
oh_hevc_put_hevc_bi1_qpel_v8_14_10_sse(uint16_t* _dst,
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
  __m128i x1, x2, x3, x4, x5, x6, x7, x8, x9, c1, c2, c3, c4;
  const int16_t* src = _src;
  const int srcstride = _srcstride;
  const __m128i offset = _mm_set1_epi16(1 << 10);
  uint16_t* dst = (uint16_t*)_dst;
  const int dststride = _dststride;
  c1 = _mm_load_si128((__m128i*)oh_hevc_qpel_filters_sse[my - 1][0]);
  c2 = _mm_load_si128((__m128i*)oh_hevc_qpel_filters_sse[my - 1][1]);
  c3 = _mm_load_si128((__m128i*)oh_hevc_qpel_filters_sse[my - 1][2]);
  c4 = _mm_load_si128((__m128i*)oh_hevc_qpel_filters_sse[my - 1][3]);
  for (y = 0; y < height; y++) {
    for (x = 0; x < width; x += 8) {
      __m128i r5;
      x1 = _mm_loadu_si128((__m128i*)&src[x - 3 * srcstride]);
      x2 = _mm_loadu_si128((__m128i*)&src[x - 2 * srcstride]);
      x3 = _mm_loadu_si128((__m128i*)&src[x - 1 * srcstride]);
      x4 = _mm_loadu_si128((__m128i*)&src[x]);
      x5 = _mm_loadu_si128((__m128i*)&src[x + srcstride]);
      x6 = _mm_loadu_si128((__m128i*)&src[x + 2 * srcstride]);
      x7 = _mm_loadu_si128((__m128i*)&src[x + 3 * srcstride]);
      x8 = _mm_loadu_si128((__m128i*)&src[x + 4 * srcstride]);
      x9 = x1;
      x1 = _mm_unpacklo_epi16(x9, x2);
      x2 = _mm_unpackhi_epi16(x9, x2);
      x9 = x3;
      x3 = _mm_unpacklo_epi16(x9, x4);
      x4 = _mm_unpackhi_epi16(x9, x4);
      x9 = x5;
      x5 = _mm_unpacklo_epi16(x9, x6);
      x6 = _mm_unpackhi_epi16(x9, x6);
      x9 = x7;
      x7 = _mm_unpacklo_epi16(x9, x8);
      x8 = _mm_unpackhi_epi16(x9, x8);
      x1 = _mm_madd_epi16(x1, c1);
      x3 = _mm_madd_epi16(x3, c2);
      x5 = _mm_madd_epi16(x5, c3);
      x7 = _mm_madd_epi16(x7, c4);
      x2 = _mm_madd_epi16(x2, c1);
      x4 = _mm_madd_epi16(x4, c2);
      x6 = _mm_madd_epi16(x6, c3);
      x8 = _mm_madd_epi16(x8, c4);
      x1 = _mm_add_epi32(x1, x3);
      x3 = _mm_add_epi32(x5, x7);
      x2 = _mm_add_epi32(x2, x4);
      x4 = _mm_add_epi32(x6, x8);
      x1 = _mm_add_epi32(x1, x3);
      x2 = _mm_add_epi32(x2, x4);
      x1 = _mm_srai_epi32(x1, 6);
      x2 = _mm_srai_epi32(x2, 6);
      x1 = _mm_packs_epi32(x1, x2);
      r5 = _mm_load_si128((__m128i*)&src2[x]);
      x1 = _mm_adds_epi16(x1, r5);
      x1 = _mm_mulhrs_epi16(x1, offset);
      x1 = _mm_max_epi16(x1, _mm_setzero_si128());
      x1 = _mm_min_epi16(x1, _mm_set1_epi16(0x03FF));
      _mm_storeu_si128((__m128i*)&dst[x], x1);
    }
    src += srcstride;
    src2 += 128;
    dst += dststride;
  }
}

static void
oh_hevc_put_hevc_uni_qpel_v8_14_10_sse(uint16_t* _dst,
                                       ptrdiff_t _dststride,
                                       const int16_t* _src,
                                       ptrdiff_t _srcstride,
                                       int height,
                                       intptr_t mx,
                                       intptr_t my,
                                       int width)
{
  int x, y;
  __m128i x1, x2, x3, x4, x5, x6, x7, x8, x9, c1, c2, c3, c4;
  const int16_t* src = _src;
  const int srcstride = _srcstride;
  const __m128i offset = _mm_set1_epi16(1 << (10 + 1));
  uint16_t* dst = (uint16_t*)_dst;
  const int dststride = _dststride;
  c1 = _mm_load_si128((__m128i*)oh_hevc_qpel_filters_sse[my - 1][0]);
  c2 = _mm_load_si128((__m128i*)oh_hevc_qpel_filters_sse[my - 1][1]);
  c3 = _mm_load_si128((__m128i*)oh_hevc_qpel_filters_sse[my - 1][2]);
  c4 = _mm_load_si128((__m128i*)oh_hevc_qpel_filters_sse[my - 1][3]);
  for (y = 0; y < height; y++) {
    for (x = 0; x < width; x += 8) {
      x1 = _mm_loadu_si128((__m128i*)&src[x - 3 * srcstride]);
      x2 = _mm_loadu_si128((__m128i*)&src[x - 2 * srcstride]);
      x3 = _mm_loadu_si128((__m128i*)&src[x - 1 * srcstride]);
      x4 = _mm_loadu_si128((__m128i*)&src[x]);
      x5 = _mm_loadu_si128((__m128i*)&src[x + srcstride]);
      x6 = _mm_loadu_si128((__m128i*)&src[x + 2 * srcstride]);
      x7 = _mm_loadu_si128((__m128i*)&src[x + 3 * srcstride]);
      x8 = _mm_loadu_si128((__m128i*)&src[x + 4 * srcstride]);
      x9 = x1;
      x1 = _mm_unpacklo_epi16(x9, x2);
      x2 = _mm_unpackhi_epi16(x9, x2);
      x9 = x3;
      x3 = _mm_unpacklo_epi16(x9, x4);
      x4 = _mm_unpackhi_epi16(x9, x4);
      x9 = x5;
      x5 = _mm_unpacklo_epi16(x9, x6);
      x6 = _mm_unpackhi_epi16(x9, x6);
      x9 = x7;
      x7 = _mm_unpacklo_epi16(x9, x8);
      x8 = _mm_unpackhi_epi16(x9, x8);
      x1 = _mm_madd_epi16(x1, c1);
      x3 = _mm_madd_epi16(x3, c2);
      x5 = _mm_madd_epi16(x5, c3);
      x7 = _mm_madd_epi16(x7, c4);
      x2 = _mm_madd_epi16(x2, c1);
      x4 = _mm_madd_epi16(x4, c2);
      x6 = _mm_madd_epi16(x6, c3);
      x8 = _mm_madd_epi16(x8, c4);
      x1 = _mm_add_epi32(x1, x3);
      x3 = _mm_add_epi32(x5, x7);
      x2 = _mm_add_epi32(x2, x4);
      x4 = _mm_add_epi32(x6, x8);
      x1 = _mm_add_epi32(x1, x3);
      x2 = _mm_add_epi32(x2, x4);
      x1 = _mm_srai_epi32(x1, 6);
      x2 = _mm_srai_epi32(x2, 6);
      x1 = _mm_packs_epi32(x1, x2);
      x1 = _mm_mulhrs_epi16(x1, offset);
      x1 = _mm_max_epi16(x1, _mm_setzero_si128());
      x1 = _mm_min_epi16(x1, _mm_set1_epi16(0x03FF));
      _mm_storeu_si128((__m128i*)&dst[x], x1);
    }
    src += srcstride;
    dst += dststride;
  }
}

static void
oh_hevc_put_hevc_bi0_qpel_hv4_10_sse(int16_t* dst,
                                     const uint16_t* _src,
                                     ptrdiff_t _srcstride,
                                     int height,
                                     intptr_t mx,
                                     intptr_t my,
                                     int width)
{
  int16_t tmp_array[(128 + 7) * 128];
  int16_t* tmp = tmp_array;
  const uint16_t* src = _src;
  const int srcstride = _srcstride;
  src -= 3 * srcstride;
  oh_hevc_put_hevc_bi0_qpel_h4_10_sse(
    tmp, src, _srcstride, height + 7, mx, my, width);
  tmp = tmp_array + 3 * 128;
  oh_hevc_put_hevc_bi0_qpel_v4_14_sse(
    dst, (const uint16_t*)tmp, 128, height, mx, my, width);
}

static void
oh_hevc_put_hevc_bi1_qpel_hv4_10_sse(uint16_t* dst,
                                     ptrdiff_t dststride,
                                     const uint16_t* _src,
                                     ptrdiff_t _srcstride,
                                     const int16_t* src2,
                                     int height,
                                     intptr_t mx,
                                     intptr_t my,
                                     int width)
{
  int16_t tmp_array[(128 + 7) * 128];
  int16_t* tmp = tmp_array;
  const uint16_t* src = _src;
  const int srcstride = _srcstride;
  src -= 3 * srcstride;
  oh_hevc_put_hevc_bi0_qpel_h4_10_sse(
    tmp, src, _srcstride, height + 7, mx, my, width);
  tmp = tmp_array + 3 * 128;
  oh_hevc_put_hevc_bi1_qpel_v4_14_10_sse(
    dst, dststride, tmp, 128, src2, height, mx, my, width);
}

static void
oh_hevc_put_hevc_uni_qpel_hv4_10_sse(uint16_t* dst,
                                     ptrdiff_t dststride,
                                     const uint16_t* _src,
                                     ptrdiff_t _srcstride,
                                     int height,
                                     intptr_t mx,
                                     intptr_t my,
                                     int width)
{
  int16_t tmp_array[(128 + 7) * 128];
  int16_t* tmp = tmp_array;
  const uint16_t* src = _src;
  const int srcstride = _srcstride;
  src -= 3 * srcstride;
  oh_hevc_put_hevc_bi0_qpel_h4_10_sse(
    tmp, src, _srcstride, height + 7, mx, my, width);
  tmp = tmp_array + 3 * 128;
  oh_hevc_put_hevc_uni_qpel_v4_14_10_sse(
    dst, dststride, tmp, 128, height, mx, my, width);
}

static void
oh_hevc_put_hevc_bi0_qpel_hv8_10_sse(int16_t* dst,
                                     const uint16_t* _src,
                                     ptrdiff_t _srcstride,
                                     int height,
                                     intptr_t mx,
                                     intptr_t my,
                                     int width)
{
  int16_t tmp_array[(128 + 7) * 128];
  int16_t* tmp = tmp_array;
  const uint16_t* src = _src;
  const int srcstride = _srcstride;
  src -= 3 * srcstride;
  oh_hevc_put_hevc_bi0_qpel_h8_10_sse(
    tmp, src, _srcstride, height + 7, mx, my, width);
  tmp = tmp_array + 3 * 128;
  oh_hevc_put_hevc_bi0_qpel_v8_14_sse(
    dst, (const uint16_t*)tmp, 128, height, mx, my, width);
}

static void
oh_hevc_put_hevc_bi1_qpel_hv8_10_sse(uint16_t* dst,
                                     ptrdiff_t dststride,
                                     const uint16_t* _src,
                                     ptrdiff_t _srcstride,
                                     const int16_t* src2,
                                     int height,
                                     intptr_t mx,
                                     intptr_t my,
                                     int width)
{
  int16_t tmp_array[(128 + 7) * 128];
  int16_t* tmp = tmp_array;
  const uint16_t* src = _src;
  const int srcstride = _srcstride;
  src -= 3 * srcstride;
  oh_hevc_put_hevc_bi0_qpel_h8_10_sse(
    tmp, src, _srcstride, height + 7, mx, my, width);
  tmp = tmp_array + 3 * 128;
  oh_hevc_put_hevc_bi1_qpel_v8_14_10_sse(
    dst, dststride, tmp, 128, src2, height, mx, my, width);
}

static void
oh_hevc_put_hevc_uni_qpel_hv8_10_sse(uint16_t* dst,
                                     ptrdiff_t dststride,
                                     const uint16_t* _src,
                                     ptrdiff_t _srcstride,
                                     int height,
                                     intptr_t mx,
                                     intptr_t my,
                                     int width)
{
  int16_t tmp_array[(128 + 7) * 128];
  int16_t* tmp = tmp_array;
  const uint16_t* src = _src;
  const int srcstride = _srcstride;
  src -= 3 * srcstride;
  oh_hevc_put_hevc_bi0_qpel_h8_10_sse(
    tmp, src, _srcstride, height + 7, mx, my, width);
  tmp = tmp_array + 3 * 128;
  oh_hevc_put_hevc_uni_qpel_v8_14_10_sse(
    dst, dststride, tmp, 128, height, mx, my, width);
}

static void
oh_hevc_put_hevc_bi0_pel_pixels16_10_sse(int16_t* dst,
                                         const uint16_t* _src,
                                         ptrdiff_t _srcstride,
                                         int height,
                                         intptr_t mx,
                                         intptr_t my,
                                         int width)
{
  oh_hevc_put_hevc_bi0_pel_pixels8_10_sse(
    dst, _src, _srcstride, height, mx, my, width);
}

static void
oh_hevc_put_hevc_bi1_pel_pixels16_10_sse(uint16_t* dst,
                                         ptrdiff_t dststride,
                                         const uint16_t* _src,
                                         ptrdiff_t _srcstride,
                                         const int16_t* src2,
                                         int height,
                                         intptr_t mx,
                                         intptr_t my,
                                         int width)
{
  oh_hevc_put_hevc_bi1_pel_pixels8_10_sse(
    dst, dststride, _src, _srcstride, src2, height, mx, my, width);
}

static void
oh_hevc_put_hevc_uni_pel_pixels16_10_sse(uint16_t* dst,
                                         ptrdiff_t dststride,
                                         const uint16_t* _src,
                                         ptrdiff_t _srcstride,
                                         int height,
                                         intptr_t mx,
                                         intptr_t my,
                                         int width)
{
  oh_hevc_put_hevc_uni_pel_pixels8_10_sse(
    dst, dststride, _src, _srcstride, height, mx, my, width);
}

static void
oh_hevc_put_hevc_bi0_pel_pixels32_10_sse(int16_t* dst,
                                         const uint16_t* _src,
                                         ptrdiff_t _srcstride,
                                         int height,
                                         intptr_t mx,
                                         intptr_t my,
                                         int width)
{
  oh_hevc_put_hevc_bi0_pel_pixels8_10_sse(
    dst, _src, _srcstride, height, mx, my, width);
}

static void
oh_hevc_put_hevc_bi1_pel_pixels32_10_sse(uint16_t* dst,
                                         ptrdiff_t dststride,
                                         const uint16_t* _src,
                                         ptrdiff_t _srcstride,
                                         const int16_t* src2,
                                         int height,
                                         intptr_t mx,
                                         intptr_t my,
                                         int width)
{
  oh_hevc_put_hevc_bi1_pel_pixels8_10_sse(
    dst, dststride, _src, _srcstride, src2, height, mx, my, width);
}

static void
oh_hevc_put_hevc_uni_pel_pixels32_10_sse(uint16_t* dst,
                                         ptrdiff_t dststride,
                                         const uint16_t* _src,
                                         ptrdiff_t _srcstride,
                                         int height,
                                         intptr_t mx,
                                         intptr_t my,
                                         int width)
{
  oh_hevc_put_hevc_uni_pel_pixels8_10_sse(
    dst, dststride, _src, _srcstride, height, mx, my, width);
}

static void
oh_hevc_put_hevc_bi0_pel_pixels64_10_sse(int16_t* dst,
                                         const uint16_t* _src,
                                         ptrdiff_t _srcstride,
                                         int height,
                                         intptr_t mx,
                                         intptr_t my,
                                         int width)
{
  oh_hevc_put_hevc_bi0_pel_pixels8_10_sse(
    dst, _src, _srcstride, height, mx, my, width);
}

static void
oh_hevc_put_hevc_bi1_pel_pixels64_10_sse(uint16_t* dst,
                                         ptrdiff_t dststride,
                                         const uint16_t* _src,
                                         ptrdiff_t _srcstride,
                                         const int16_t* src2,
                                         int height,
                                         intptr_t mx,
                                         intptr_t my,
                                         int width)
{
  oh_hevc_put_hevc_bi1_pel_pixels8_10_sse(
    dst, dststride, _src, _srcstride, src2, height, mx, my, width);
}

static void
oh_hevc_put_hevc_uni_pel_pixels64_10_sse(uint16_t* dst,
                                         ptrdiff_t dststride,
                                         const uint16_t* _src,
                                         ptrdiff_t _srcstride,
                                         int height,
                                         intptr_t mx,
                                         intptr_t my,
                                         int width)
{
  oh_hevc_put_hevc_uni_pel_pixels8_10_sse(
    dst, dststride, _src, _srcstride, height, mx, my, width);
}

static void
oh_hevc_put_hevc_bi0_qpel_h16_10_sse(int16_t* dst,
                                     const uint16_t* _src,
                                     ptrdiff_t _srcstride,
                                     int height,
                                     intptr_t mx,
                                     intptr_t my,
                                     int width)
{
  oh_hevc_put_hevc_bi0_qpel_h8_10_sse(
    dst, _src, _srcstride, height, mx, my, width);
}

static void
oh_hevc_put_hevc_bi1_qpel_h16_10_sse(uint16_t* dst,
                                     ptrdiff_t dststride,
                                     const uint16_t* _src,
                                     ptrdiff_t _srcstride,
                                     const int16_t* src2,
                                     int height,
                                     intptr_t mx,
                                     intptr_t my,
                                     int width)
{
  oh_hevc_put_hevc_bi1_qpel_h8_10_sse(
    dst, dststride, _src, _srcstride, src2, height, mx, my, width);
}

static void
oh_hevc_put_hevc_uni_qpel_h16_10_sse(uint16_t* dst,
                                     ptrdiff_t dststride,
                                     const uint16_t* _src,
                                     ptrdiff_t _srcstride,
                                     int height,
                                     intptr_t mx,
                                     intptr_t my,
                                     int width)
{
  oh_hevc_put_hevc_uni_qpel_h8_10_sse(
    dst, dststride, _src, _srcstride, height, mx, my, width);
}

static void
oh_hevc_put_hevc_bi0_qpel_h32_10_sse(int16_t* dst,
                                     const uint16_t* _src,
                                     ptrdiff_t _srcstride,
                                     int height,
                                     intptr_t mx,
                                     intptr_t my,
                                     int width)
{
  oh_hevc_put_hevc_bi0_qpel_h8_10_sse(
    dst, _src, _srcstride, height, mx, my, width);
}

static void
oh_hevc_put_hevc_bi1_qpel_h32_10_sse(uint16_t* dst,
                                     ptrdiff_t dststride,
                                     const uint16_t* _src,
                                     ptrdiff_t _srcstride,
                                     const int16_t* src2,
                                     int height,
                                     intptr_t mx,
                                     intptr_t my,
                                     int width)
{
  oh_hevc_put_hevc_bi1_qpel_h8_10_sse(
    dst, dststride, _src, _srcstride, src2, height, mx, my, width);
}

static void
oh_hevc_put_hevc_uni_qpel_h32_10_sse(uint16_t* dst,
                                     ptrdiff_t dststride,
                                     const uint16_t* _src,
                                     ptrdiff_t _srcstride,
                                     int height,
                                     intptr_t mx,
                                     intptr_t my,
                                     int width)
{
  oh_hevc_put_hevc_uni_qpel_h8_10_sse(
    dst, dststride, _src, _srcstride, height, mx, my, width);
}

static void
oh_hevc_put_hevc_bi0_qpel_h64_10_sse(int16_t* dst,
                                     const uint16_t* _src,
                                     ptrdiff_t _srcstride,
                                     int height,
                                     intptr_t mx,
                                     intptr_t my,
                                     int width)
{
  oh_hevc_put_hevc_bi0_qpel_h8_10_sse(
    dst, _src, _srcstride, height, mx, my, width);
}

static void
oh_hevc_put_hevc_bi1_qpel_h64_10_sse(uint16_t* dst,
                                     ptrdiff_t dststride,
                                     const uint16_t* _src,
                                     ptrdiff_t _srcstride,
                                     const int16_t* src2,
                                     int height,
                                     intptr_t mx,
                                     intptr_t my,
                                     int width)
{
  oh_hevc_put_hevc_bi1_qpel_h8_10_sse(
    dst, dststride, _src, _srcstride, src2, height, mx, my, width);
}

static void
oh_hevc_put_hevc_uni_qpel_h64_10_sse(uint16_t* dst,
                                     ptrdiff_t dststride,
                                     const uint16_t* _src,
                                     ptrdiff_t _srcstride,
                                     int height,
                                     intptr_t mx,
                                     intptr_t my,
                                     int width)
{
  oh_hevc_put_hevc_uni_qpel_h8_10_sse(
    dst, dststride, _src, _srcstride, height, mx, my, width);
}

static void
oh_hevc_put_hevc_bi0_qpel_v16_10_sse(int16_t* dst,
                                     const uint16_t* _src,
                                     ptrdiff_t _srcstride,
                                     int height,
                                     intptr_t mx,
                                     intptr_t my,
                                     int width)
{
  oh_hevc_put_hevc_bi0_qpel_v8_10_sse(
    dst, _src, _srcstride, height, mx, my, width);
}

static void
oh_hevc_put_hevc_bi1_qpel_v16_10_sse(uint16_t* dst,
                                     ptrdiff_t dststride,
                                     const uint16_t* _src,
                                     ptrdiff_t _srcstride,
                                     const int16_t* src2,
                                     int height,
                                     intptr_t mx,
                                     intptr_t my,
                                     int width)
{
  oh_hevc_put_hevc_bi1_qpel_v8_10_sse(
    dst, dststride, _src, _srcstride, src2, height, mx, my, width);
}

static void
oh_hevc_put_hevc_uni_qpel_v16_10_sse(uint16_t* dst,
                                     ptrdiff_t dststride,
                                     const uint16_t* _src,
                                     ptrdiff_t _srcstride,
                                     int height,
                                     intptr_t mx,
                                     intptr_t my,
                                     int width)
{
  oh_hevc_put_hevc_uni_qpel_v8_10_sse(
    dst, dststride, _src, _srcstride, height, mx, my, width);
}

static void
oh_hevc_put_hevc_bi0_qpel_v32_10_sse(int16_t* dst,
                                     const uint16_t* _src,
                                     ptrdiff_t _srcstride,
                                     int height,
                                     intptr_t mx,
                                     intptr_t my,
                                     int width)
{
  oh_hevc_put_hevc_bi0_qpel_v8_10_sse(
    dst, _src, _srcstride, height, mx, my, width);
}

static void
oh_hevc_put_hevc_bi1_qpel_v32_10_sse(uint16_t* dst,
                                     ptrdiff_t dststride,
                                     const uint16_t* _src,
                                     ptrdiff_t _srcstride,
                                     const int16_t* src2,
                                     int height,
                                     intptr_t mx,
                                     intptr_t my,
                                     int width)
{
  oh_hevc_put_hevc_bi1_qpel_v8_10_sse(
    dst, dststride, _src, _srcstride, src2, height, mx, my, width);
}

static void
oh_hevc_put_hevc_uni_qpel_v32_10_sse(uint16_t* dst,
                                     ptrdiff_t dststride,
                                     const uint16_t* _src,
                                     ptrdiff_t _srcstride,
                                     int height,
                                     intptr_t mx,
                                     intptr_t my,
                                     int width)
{
  oh_hevc_put_hevc_uni_qpel_v8_10_sse(
    dst, dststride, _src, _srcstride, height, mx, my, width);
}

static void
oh_hevc_put_hevc_bi0_qpel_v64_10_sse(int16_t* dst,
                                     const uint16_t* _src,
                                     ptrdiff_t _srcstride,
                                     int height,
                                     intptr_t mx,
                                     intptr_t my,
                                     int width)
{
  oh_hevc_put_hevc_bi0_qpel_v8_10_sse(
    dst, _src, _srcstride, height, mx, my, width);
}

static void
oh_hevc_put_hevc_bi1_qpel_v64_10_sse(uint16_t* dst,
                                     ptrdiff_t dststride,
                                     const uint16_t* _src,
                                     ptrdiff_t _srcstride,
                                     const int16_t* src2,
                                     int height,
                                     intptr_t mx,
                                     intptr_t my,
                                     int width)
{
  oh_hevc_put_hevc_bi1_qpel_v8_10_sse(
    dst, dststride, _src, _srcstride, src2, height, mx, my, width);
}

static void
oh_hevc_put_hevc_uni_qpel_v64_10_sse(uint16_t* dst,
                                     ptrdiff_t dststride,
                                     const uint16_t* _src,
                                     ptrdiff_t _srcstride,
                                     int height,
                                     intptr_t mx,
                                     intptr_t my,
                                     int width)
{
  oh_hevc_put_hevc_uni_qpel_v8_10_sse(
    dst, dststride, _src, _srcstride, height, mx, my, width);
}

static void
oh_hevc_put_hevc_bi0_qpel_hv16_10_sse(int16_t* dst,
                                      const uint16_t* _src,
                                      ptrdiff_t _srcstride,
                                      int height,
                                      intptr_t mx,
                                      intptr_t my,
                                      int width)
{
  oh_hevc_put_hevc_bi0_qpel_hv8_10_sse(
    dst, _src, _srcstride, height, mx, my, width);
}

static void
oh_hevc_put_hevc_bi1_qpel_hv16_10_sse(uint16_t* dst,
                                      ptrdiff_t dststride,
                                      const uint16_t* _src,
                                      ptrdiff_t _srcstride,
                                      const int16_t* src2,
                                      int height,
                                      intptr_t mx,
                                      intptr_t my,
                                      int width)
{
  oh_hevc_put_hevc_bi1_qpel_hv8_10_sse(
    dst, dststride, _src, _srcstride, src2, height, mx, my, width);
}

static void
oh_hevc_put_hevc_uni_qpel_hv16_10_sse(uint16_t* dst,
                                      ptrdiff_t dststride,
                                      const uint16_t* _src,
                                      ptrdiff_t _srcstride,
                                      int height,
                                      intptr_t mx,
                                      intptr_t my,
                                      int width)
{
  oh_hevc_put_hevc_uni_qpel_hv8_10_sse(
    dst, dststride, _src, _srcstride, height, mx, my, width);
}

static void
oh_hevc_put_hevc_bi0_qpel_hv32_10_sse(int16_t* dst,
                                      const uint16_t* _src,
                                      ptrdiff_t _srcstride,
                                      int height,
                                      intptr_t mx,
                                      intptr_t my,
                                      int width)
{
  oh_hevc_put_hevc_bi0_qpel_hv8_10_sse(
    dst, _src, _srcstride, height, mx, my, width);
}

static void
oh_hevc_put_hevc_bi1_qpel_hv32_10_sse(uint16_t* dst,
                                      ptrdiff_t dststride,
                                      const uint16_t* _src,
                                      ptrdiff_t _srcstride,
                                      const int16_t* src2,
                                      int height,
                                      intptr_t mx,
                                      intptr_t my,
                                      int width)
{
  oh_hevc_put_hevc_bi1_qpel_hv8_10_sse(
    dst, dststride, _src, _srcstride, src2, height, mx, my, width);
}

static void
oh_hevc_put_hevc_uni_qpel_hv32_10_sse(uint16_t* dst,
                                      ptrdiff_t dststride,
                                      const uint16_t* _src,
                                      ptrdiff_t _srcstride,
                                      int height,
                                      intptr_t mx,
                                      intptr_t my,
                                      int width)
{
  oh_hevc_put_hevc_uni_qpel_hv8_10_sse(
    dst, dststride, _src, _srcstride, height, mx, my, width);
}

static void
oh_hevc_put_hevc_bi0_qpel_hv64_10_sse(int16_t* dst,
                                      const uint16_t* _src,
                                      ptrdiff_t _srcstride,
                                      int height,
                                      intptr_t mx,
                                      intptr_t my,
                                      int width)
{
  oh_hevc_put_hevc_bi0_qpel_hv8_10_sse(
    dst, _src, _srcstride, height, mx, my, width);
}

static void
oh_hevc_put_hevc_bi1_qpel_hv64_10_sse(uint16_t* dst,
                                      ptrdiff_t dststride,
                                      const uint16_t* _src,
                                      ptrdiff_t _srcstride,
                                      const int16_t* src2,
                                      int height,
                                      intptr_t mx,
                                      intptr_t my,
                                      int width)
{
  oh_hevc_put_hevc_bi1_qpel_hv8_10_sse(
    dst, dststride, _src, _srcstride, src2, height, mx, my, width);
}

static void
oh_hevc_put_hevc_uni_qpel_hv64_10_sse(uint16_t* dst,
                                      ptrdiff_t dststride,
                                      const uint16_t* _src,
                                      ptrdiff_t _srcstride,
                                      int height,
                                      intptr_t mx,
                                      intptr_t my,
                                      int width)
{
  oh_hevc_put_hevc_uni_qpel_hv8_10_sse(
    dst, dststride, _src, _srcstride, height, mx, my, width);
}

static void
oh_hevc_put_hevc_bi0_epel_h16_10_sse(int16_t* dst,
                                     const uint16_t* _src,
                                     ptrdiff_t _srcstride,
                                     int height,
                                     intptr_t mx,
                                     intptr_t my,
                                     int width)
{
  oh_hevc_put_hevc_bi0_epel_h8_10_sse(
    dst, _src, _srcstride, height, mx, my, width);
}

static void
oh_hevc_put_hevc_bi1_epel_h16_10_sse(uint16_t* dst,
                                     ptrdiff_t dststride,
                                     const uint16_t* _src,
                                     ptrdiff_t _srcstride,
                                     const int16_t* src2,
                                     int height,
                                     intptr_t mx,
                                     intptr_t my,
                                     int width)
{
  oh_hevc_put_hevc_bi1_epel_h8_10_sse(
    dst, dststride, _src, _srcstride, src2, height, mx, my, width);
}

static void
oh_hevc_put_hevc_uni_epel_h16_10_sse(uint16_t* dst,
                                     ptrdiff_t dststride,
                                     const uint16_t* _src,
                                     ptrdiff_t _srcstride,
                                     int height,
                                     intptr_t mx,
                                     intptr_t my,
                                     int width)
{
  oh_hevc_put_hevc_uni_epel_h8_10_sse(
    dst, dststride, _src, _srcstride, height, mx, my, width);
}

static void
oh_hevc_put_hevc_bi0_epel_h32_10_sse(int16_t* dst,
                                     const uint16_t* _src,
                                     ptrdiff_t _srcstride,
                                     int height,
                                     intptr_t mx,
                                     intptr_t my,
                                     int width)
{
  oh_hevc_put_hevc_bi0_epel_h8_10_sse(
    dst, _src, _srcstride, height, mx, my, width);
}

static void
oh_hevc_put_hevc_bi1_epel_h32_10_sse(uint16_t* dst,
                                     ptrdiff_t dststride,
                                     const uint16_t* _src,
                                     ptrdiff_t _srcstride,
                                     const int16_t* src2,
                                     int height,
                                     intptr_t mx,
                                     intptr_t my,
                                     int width)
{
  oh_hevc_put_hevc_bi1_epel_h8_10_sse(
    dst, dststride, _src, _srcstride, src2, height, mx, my, width);
}

static void
oh_hevc_put_hevc_uni_epel_h32_10_sse(uint16_t* dst,
                                     ptrdiff_t dststride,
                                     const uint16_t* _src,
                                     ptrdiff_t _srcstride,
                                     int height,
                                     intptr_t mx,
                                     intptr_t my,
                                     int width)
{
  oh_hevc_put_hevc_uni_epel_h8_10_sse(
    dst, dststride, _src, _srcstride, height, mx, my, width);
}

static void
oh_hevc_put_hevc_bi0_epel_h64_10_sse(int16_t* dst,
                                     const uint16_t* _src,
                                     ptrdiff_t _srcstride,
                                     int height,
                                     intptr_t mx,
                                     intptr_t my,
                                     int width)
{
  oh_hevc_put_hevc_bi0_epel_h8_10_sse(
    dst, _src, _srcstride, height, mx, my, width);
}

static void
oh_hevc_put_hevc_bi1_epel_h64_10_sse(uint16_t* dst,
                                     ptrdiff_t dststride,
                                     const uint16_t* _src,
                                     ptrdiff_t _srcstride,
                                     const int16_t* src2,
                                     int height,
                                     intptr_t mx,
                                     intptr_t my,
                                     int width)
{
  oh_hevc_put_hevc_bi1_epel_h8_10_sse(
    dst, dststride, _src, _srcstride, src2, height, mx, my, width);
}

static void
oh_hevc_put_hevc_uni_epel_h64_10_sse(uint16_t* dst,
                                     ptrdiff_t dststride,
                                     const uint16_t* _src,
                                     ptrdiff_t _srcstride,
                                     int height,
                                     intptr_t mx,
                                     intptr_t my,
                                     int width)
{
  oh_hevc_put_hevc_uni_epel_h8_10_sse(
    dst, dststride, _src, _srcstride, height, mx, my, width);
}

static void
oh_hevc_put_hevc_bi0_epel_v16_10_sse(int16_t* dst,
                                     const uint16_t* _src,
                                     ptrdiff_t _srcstride,
                                     int height,
                                     intptr_t mx,
                                     intptr_t my,
                                     int width)
{
  oh_hevc_put_hevc_bi0_epel_v8_10_sse(
    dst, _src, _srcstride, height, mx, my, width);
}

static void
oh_hevc_put_hevc_bi1_epel_v16_10_sse(uint16_t* dst,
                                     ptrdiff_t dststride,
                                     const uint16_t* _src,
                                     ptrdiff_t _srcstride,
                                     const int16_t* src2,
                                     int height,
                                     intptr_t mx,
                                     intptr_t my,
                                     int width)
{
  oh_hevc_put_hevc_bi1_epel_v8_10_sse(
    dst, dststride, _src, _srcstride, src2, height, mx, my, width);
}

static void
oh_hevc_put_hevc_uni_epel_v16_10_sse(uint16_t* dst,
                                     ptrdiff_t dststride,
                                     const uint16_t* _src,
                                     ptrdiff_t _srcstride,
                                     int height,
                                     intptr_t mx,
                                     intptr_t my,
                                     int width)
{
  oh_hevc_put_hevc_uni_epel_v8_10_sse(
    dst, dststride, _src, _srcstride, height, mx, my, width);
}

static void
oh_hevc_put_hevc_bi0_epel_v32_10_sse(int16_t* dst,
                                     const uint16_t* _src,
                                     ptrdiff_t _srcstride,
                                     int height,
                                     intptr_t mx,
                                     intptr_t my,
                                     int width)
{
  oh_hevc_put_hevc_bi0_epel_v8_10_sse(
    dst, _src, _srcstride, height, mx, my, width);
}

static void
oh_hevc_put_hevc_bi1_epel_v32_10_sse(uint16_t* dst,
                                     ptrdiff_t dststride,
                                     const uint16_t* _src,
                                     ptrdiff_t _srcstride,
                                     const int16_t* src2,
                                     int height,
                                     intptr_t mx,
                                     intptr_t my,
                                     int width)
{
  oh_hevc_put_hevc_bi1_epel_v8_10_sse(
    dst, dststride, _src, _srcstride, src2, height, mx, my, width);
}

static void
oh_hevc_put_hevc_uni_epel_v32_10_sse(uint16_t* dst,
                                     ptrdiff_t dststride,
                                     const uint16_t* _src,
                                     ptrdiff_t _srcstride,
                                     int height,
                                     intptr_t mx,
                                     intptr_t my,
                                     int width)
{
  oh_hevc_put_hevc_uni_epel_v8_10_sse(
    dst, dststride, _src, _srcstride, height, mx, my, width);
}

static void
oh_hevc_put_hevc_bi0_epel_v64_10_sse(int16_t* dst,
                                     const uint16_t* _src,
                                     ptrdiff_t _srcstride,
                                     int height,
                                     intptr_t mx,
                                     intptr_t my,
                                     int width)
{
  oh_hevc_put_hevc_bi0_epel_v8_10_sse(
    dst, _src, _srcstride, height, mx, my, width);
}

static void
oh_hevc_put_hevc_bi1_epel_v64_10_sse(uint16_t* dst,
                                     ptrdiff_t dststride,
                                     const uint16_t* _src,
                                     ptrdiff_t _srcstride,
                                     const int16_t* src2,
                                     int height,
                                     intptr_t mx,
                                     intptr_t my,
                                     int width)
{
  oh_hevc_put_hevc_bi1_epel_v8_10_sse(
    dst, dststride, _src, _srcstride, src2, height, mx, my, width);
}

static void
oh_hevc_put_hevc_uni_epel_v64_10_sse(uint16_t* dst,
                                     ptrdiff_t dststride,
                                     const uint16_t* _src,
                                     ptrdiff_t _srcstride,
                                     int height,
                                     intptr_t mx,
                                     intptr_t my,
                                     int width)
{
  oh_hevc_put_hevc_uni_epel_v8_10_sse(
    dst, dststride, _src, _srcstride, height, mx, my, width);
}

static void
oh_hevc_put_hevc_bi0_epel_hv16_10_sse(int16_t* dst,
                                      const uint16_t* _src,
                                      ptrdiff_t _srcstride,
                                      int height,
                                      intptr_t mx,
                                      intptr_t my,
                                      int width)
{
  oh_hevc_put_hevc_bi0_epel_hv8_10_sse(
    dst, _src, _srcstride, height, mx, my, width);
}

static void
oh_hevc_put_hevc_bi1_epel_hv16_10_sse(uint16_t* dst,
                                      ptrdiff_t dststride,
                                      const uint16_t* _src,
                                      ptrdiff_t _srcstride,
                                      const int16_t* src2,
                                      int height,
                                      intptr_t mx,
                                      intptr_t my,
                                      int width)
{
  oh_hevc_put_hevc_bi1_epel_hv8_10_sse(
    dst, dststride, _src, _srcstride, src2, height, mx, my, width);
}

static void
oh_hevc_put_hevc_uni_epel_hv16_10_sse(uint16_t* dst,
                                      ptrdiff_t dststride,
                                      const uint16_t* _src,
                                      ptrdiff_t _srcstride,
                                      int height,
                                      intptr_t mx,
                                      intptr_t my,
                                      int width)
{
  oh_hevc_put_hevc_uni_epel_hv8_10_sse(
    dst, dststride, _src, _srcstride, height, mx, my, width);
}

static void
oh_hevc_put_hevc_bi0_epel_hv32_10_sse(int16_t* dst,
                                      const uint16_t* _src,
                                      ptrdiff_t _srcstride,
                                      int height,
                                      intptr_t mx,
                                      intptr_t my,
                                      int width)
{
  oh_hevc_put_hevc_bi0_epel_hv8_10_sse(
    dst, _src, _srcstride, height, mx, my, width);
}

static void
oh_hevc_put_hevc_bi1_epel_hv32_10_sse(uint16_t* dst,
                                      ptrdiff_t dststride,
                                      const uint16_t* _src,
                                      ptrdiff_t _srcstride,
                                      const int16_t* src2,
                                      int height,
                                      intptr_t mx,
                                      intptr_t my,
                                      int width)
{
  oh_hevc_put_hevc_bi1_epel_hv8_10_sse(
    dst, dststride, _src, _srcstride, src2, height, mx, my, width);
}

static void
oh_hevc_put_hevc_uni_epel_hv32_10_sse(uint16_t* dst,
                                      ptrdiff_t dststride,
                                      const uint16_t* _src,
                                      ptrdiff_t _srcstride,
                                      int height,
                                      intptr_t mx,
                                      intptr_t my,
                                      int width)
{
  oh_hevc_put_hevc_uni_epel_hv8_10_sse(
    dst, dststride, _src, _srcstride, height, mx, my, width);
}

static void
oh_hevc_put_hevc_bi0_epel_hv64_10_sse(int16_t* dst,
                                      const uint16_t* _src,
                                      ptrdiff_t _srcstride,
                                      int height,
                                      intptr_t mx,
                                      intptr_t my,
                                      int width)
{
  oh_hevc_put_hevc_bi0_epel_hv8_10_sse(
    dst, _src, _srcstride, height, mx, my, width);
}

static void
oh_hevc_put_hevc_bi1_epel_hv64_10_sse(uint16_t* dst,
                                      ptrdiff_t dststride,
                                      const uint16_t* _src,
                                      ptrdiff_t _srcstride,
                                      const int16_t* src2,
                                      int height,
                                      intptr_t mx,
                                      intptr_t my,
                                      int width)
{
  oh_hevc_put_hevc_bi1_epel_hv8_10_sse(
    dst, dststride, _src, _srcstride, src2, height, mx, my, width);
}

static void
oh_hevc_put_hevc_uni_epel_hv64_10_sse(uint16_t* dst,
                                      ptrdiff_t dststride,
                                      const uint16_t* _src,
                                      ptrdiff_t _srcstride,
                                      int height,
                                      intptr_t mx,
                                      intptr_t my,
                                      int width)
{
  oh_hevc_put_hevc_uni_epel_hv8_10_sse(
    dst, dststride, _src, _srcstride, height, mx, my, width);
}

void
rcn_init_mc_functions_sse(struct RCNFunctions* const rcn_funcs)
{
  struct MCFunctions* const mc_l = &rcn_funcs->mc_l;
  struct MCFunctions* const mc_c = &rcn_funcs->mc_c;

  /* Luma functions */
  mc_l->unidir[0][SIZE_BLOCK_4] = &oh_hevc_put_hevc_uni_pel_pixels4_10_sse;
  mc_l->bidir0[0][SIZE_BLOCK_4] = &oh_hevc_put_hevc_bi0_pel_pixels4_10_sse;
  mc_l->bidir1[0][SIZE_BLOCK_4] = &oh_hevc_put_hevc_bi1_pel_pixels4_10_sse;

  mc_l->unidir[1][SIZE_BLOCK_4] = &oh_hevc_put_hevc_uni_qpel_h4_10_sse;
  mc_l->bidir0[1][SIZE_BLOCK_4] = &oh_hevc_put_hevc_bi0_qpel_h4_10_sse;
  mc_l->bidir1[1][SIZE_BLOCK_4] = &oh_hevc_put_hevc_bi1_qpel_h4_10_sse;

  mc_l->unidir[2][SIZE_BLOCK_4] = &oh_hevc_put_hevc_uni_qpel_v4_10_sse;
  mc_l->bidir0[2][SIZE_BLOCK_4] = &oh_hevc_put_hevc_bi0_qpel_v4_10_sse;
  mc_l->bidir1[2][SIZE_BLOCK_4] = &oh_hevc_put_hevc_bi1_qpel_v4_10_sse;

  mc_l->unidir[3][SIZE_BLOCK_4] = &oh_hevc_put_hevc_uni_qpel_hv4_10_sse;
  mc_l->bidir0[3][SIZE_BLOCK_4] = &oh_hevc_put_hevc_bi0_qpel_hv4_10_sse;
  mc_l->bidir1[3][SIZE_BLOCK_4] = &oh_hevc_put_hevc_bi1_qpel_hv4_10_sse;

  mc_l->unidir[0][SIZE_BLOCK_8] = &oh_hevc_put_hevc_uni_pel_pixels8_10_sse;
  mc_l->bidir0[0][SIZE_BLOCK_8] = &oh_hevc_put_hevc_bi0_pel_pixels8_10_sse;
  mc_l->bidir1[0][SIZE_BLOCK_8] = &oh_hevc_put_hevc_bi1_pel_pixels8_10_sse;

  mc_l->unidir[1][SIZE_BLOCK_8] = &oh_hevc_put_hevc_uni_qpel_h8_10_sse;
  mc_l->bidir0[1][SIZE_BLOCK_8] = &oh_hevc_put_hevc_bi0_qpel_h8_10_sse;
  mc_l->bidir1[1][SIZE_BLOCK_8] = &oh_hevc_put_hevc_bi1_qpel_h8_10_sse;

  mc_l->unidir[2][SIZE_BLOCK_8] = &oh_hevc_put_hevc_uni_qpel_v8_10_sse;
  mc_l->bidir0[2][SIZE_BLOCK_8] = &oh_hevc_put_hevc_bi0_qpel_v8_10_sse;
  mc_l->bidir1[2][SIZE_BLOCK_8] = &oh_hevc_put_hevc_bi1_qpel_v8_10_sse;

  mc_l->unidir[3][SIZE_BLOCK_8] = &oh_hevc_put_hevc_uni_qpel_hv8_10_sse;
  mc_l->bidir0[3][SIZE_BLOCK_8] = &oh_hevc_put_hevc_bi0_qpel_hv8_10_sse;
  mc_l->bidir1[3][SIZE_BLOCK_8] = &oh_hevc_put_hevc_bi1_qpel_hv8_10_sse;

  mc_l->unidir[0][SIZE_BLOCK_16] = &oh_hevc_put_hevc_uni_pel_pixels16_10_sse;
  mc_l->bidir0[0][SIZE_BLOCK_16] = &oh_hevc_put_hevc_bi0_pel_pixels16_10_sse;
  mc_l->bidir1[0][SIZE_BLOCK_16] = &oh_hevc_put_hevc_bi1_pel_pixels16_10_sse;

  mc_l->unidir[1][SIZE_BLOCK_16] = &oh_hevc_put_hevc_uni_qpel_h16_10_sse;
  mc_l->bidir0[1][SIZE_BLOCK_16] = &oh_hevc_put_hevc_bi0_qpel_h16_10_sse;
  mc_l->bidir1[1][SIZE_BLOCK_16] = &oh_hevc_put_hevc_bi1_qpel_h16_10_sse;

  mc_l->unidir[2][SIZE_BLOCK_16] = &oh_hevc_put_hevc_uni_qpel_v16_10_sse;
  mc_l->bidir0[2][SIZE_BLOCK_16] = &oh_hevc_put_hevc_bi0_qpel_v16_10_sse;
  mc_l->bidir1[2][SIZE_BLOCK_16] = &oh_hevc_put_hevc_bi1_qpel_v16_10_sse;

  mc_l->unidir[3][SIZE_BLOCK_16] = &oh_hevc_put_hevc_uni_qpel_hv16_10_sse;
  mc_l->bidir0[3][SIZE_BLOCK_16] = &oh_hevc_put_hevc_bi0_qpel_hv16_10_sse;
  mc_l->bidir1[3][SIZE_BLOCK_16] = &oh_hevc_put_hevc_bi1_qpel_hv16_10_sse;

  mc_l->unidir[0][SIZE_BLOCK_32] = &oh_hevc_put_hevc_uni_pel_pixels32_10_sse;
  mc_l->bidir0[0][SIZE_BLOCK_32] = &oh_hevc_put_hevc_bi0_pel_pixels32_10_sse;
  mc_l->bidir1[0][SIZE_BLOCK_32] = &oh_hevc_put_hevc_bi1_pel_pixels32_10_sse;

  mc_l->unidir[1][SIZE_BLOCK_32] = &oh_hevc_put_hevc_uni_qpel_h32_10_sse;
  mc_l->bidir0[1][SIZE_BLOCK_32] = &oh_hevc_put_hevc_bi0_qpel_h32_10_sse;
  mc_l->bidir1[1][SIZE_BLOCK_32] = &oh_hevc_put_hevc_bi1_qpel_h32_10_sse;

  mc_l->unidir[2][SIZE_BLOCK_32] = &oh_hevc_put_hevc_uni_qpel_v32_10_sse;
  mc_l->bidir0[2][SIZE_BLOCK_32] = &oh_hevc_put_hevc_bi0_qpel_v32_10_sse;
  mc_l->bidir1[2][SIZE_BLOCK_32] = &oh_hevc_put_hevc_bi1_qpel_v32_10_sse;

  mc_l->unidir[3][SIZE_BLOCK_32] = &oh_hevc_put_hevc_uni_qpel_hv32_10_sse;
  mc_l->bidir0[3][SIZE_BLOCK_32] = &oh_hevc_put_hevc_bi0_qpel_hv32_10_sse;
  mc_l->bidir1[3][SIZE_BLOCK_32] = &oh_hevc_put_hevc_bi1_qpel_hv32_10_sse;

  mc_l->unidir[0][SIZE_BLOCK_64] = &oh_hevc_put_hevc_uni_pel_pixels64_10_sse;
  mc_l->bidir0[0][SIZE_BLOCK_64] = &oh_hevc_put_hevc_bi0_pel_pixels64_10_sse;
  mc_l->bidir1[0][SIZE_BLOCK_64] = &oh_hevc_put_hevc_bi1_pel_pixels64_10_sse;

  mc_l->unidir[1][SIZE_BLOCK_64] = &oh_hevc_put_hevc_uni_qpel_h64_10_sse;
  mc_l->bidir0[1][SIZE_BLOCK_64] = &oh_hevc_put_hevc_bi0_qpel_h64_10_sse;
  mc_l->bidir1[1][SIZE_BLOCK_64] = &oh_hevc_put_hevc_bi1_qpel_h64_10_sse;

  mc_l->unidir[2][SIZE_BLOCK_64] = &oh_hevc_put_hevc_uni_qpel_v64_10_sse;
  mc_l->bidir0[2][SIZE_BLOCK_64] = &oh_hevc_put_hevc_bi0_qpel_v64_10_sse;
  mc_l->bidir1[2][SIZE_BLOCK_64] = &oh_hevc_put_hevc_bi1_qpel_v64_10_sse;

  mc_l->unidir[3][SIZE_BLOCK_64] = &oh_hevc_put_hevc_uni_qpel_hv64_10_sse;
  mc_l->bidir0[3][SIZE_BLOCK_64] = &oh_hevc_put_hevc_bi0_qpel_hv64_10_sse;
  mc_l->bidir1[3][SIZE_BLOCK_64] = &oh_hevc_put_hevc_bi1_qpel_hv64_10_sse;

  // /* Chroma functions */
  mc_c->unidir[0][SIZE_BLOCK_4] = &oh_hevc_put_hevc_uni_pel_pixels4_10_sse;
  mc_c->bidir0[0][SIZE_BLOCK_4] = &oh_hevc_put_hevc_bi0_pel_pixels4_10_sse;
  mc_c->bidir1[0][SIZE_BLOCK_4] = &oh_hevc_put_hevc_bi1_pel_pixels4_10_sse;

  mc_c->unidir[1][SIZE_BLOCK_4] = &oh_hevc_put_hevc_uni_epel_h4_10_sse;
  mc_c->bidir0[1][SIZE_BLOCK_4] = &oh_hevc_put_hevc_bi0_epel_h4_10_sse;
  mc_c->bidir1[1][SIZE_BLOCK_4] = &oh_hevc_put_hevc_bi1_epel_h4_10_sse;

  mc_c->unidir[2][SIZE_BLOCK_4] = &oh_hevc_put_hevc_uni_epel_v4_10_sse;
  mc_c->bidir0[2][SIZE_BLOCK_4] = &oh_hevc_put_hevc_bi0_epel_v4_10_sse;
  mc_c->bidir1[2][SIZE_BLOCK_4] = &oh_hevc_put_hevc_bi1_epel_v4_10_sse;

  mc_c->unidir[3][SIZE_BLOCK_4] = &oh_hevc_put_hevc_uni_epel_hv4_10_sse;
  mc_c->bidir0[3][SIZE_BLOCK_4] = &oh_hevc_put_hevc_bi0_epel_hv4_10_sse;
  mc_c->bidir1[3][SIZE_BLOCK_4] = &oh_hevc_put_hevc_bi1_epel_hv4_10_sse;

  mc_c->unidir[0][SIZE_BLOCK_8] = &oh_hevc_put_hevc_uni_pel_pixels8_10_sse;
  mc_c->bidir0[0][SIZE_BLOCK_8] = &oh_hevc_put_hevc_bi0_pel_pixels8_10_sse;
  mc_c->bidir1[0][SIZE_BLOCK_8] = &oh_hevc_put_hevc_bi1_pel_pixels8_10_sse;

  mc_c->unidir[1][SIZE_BLOCK_8] = &oh_hevc_put_hevc_uni_epel_h8_10_sse;
  mc_c->bidir0[1][SIZE_BLOCK_8] = &oh_hevc_put_hevc_bi0_epel_h8_10_sse;
  mc_c->bidir1[1][SIZE_BLOCK_8] = &oh_hevc_put_hevc_bi1_epel_h8_10_sse;

  mc_c->unidir[2][SIZE_BLOCK_8] = &oh_hevc_put_hevc_uni_epel_v8_10_sse;
  mc_c->bidir0[2][SIZE_BLOCK_8] = &oh_hevc_put_hevc_bi0_epel_v8_10_sse;
  mc_c->bidir1[2][SIZE_BLOCK_8] = &oh_hevc_put_hevc_bi1_epel_v8_10_sse;

  mc_c->unidir[3][SIZE_BLOCK_8] = &oh_hevc_put_hevc_uni_epel_hv8_10_sse;
  mc_c->bidir0[3][SIZE_BLOCK_8] = &oh_hevc_put_hevc_bi0_epel_hv8_10_sse;
  mc_c->bidir1[3][SIZE_BLOCK_8] = &oh_hevc_put_hevc_bi1_epel_hv8_10_sse;

  mc_c->unidir[0][SIZE_BLOCK_16] = &oh_hevc_put_hevc_uni_pel_pixels16_10_sse;
  mc_c->bidir0[0][SIZE_BLOCK_16] = &oh_hevc_put_hevc_bi0_pel_pixels16_10_sse;
  mc_c->bidir1[0][SIZE_BLOCK_16] = &oh_hevc_put_hevc_bi1_pel_pixels16_10_sse;

  mc_c->unidir[1][SIZE_BLOCK_16] = &oh_hevc_put_hevc_uni_epel_h16_10_sse;
  mc_c->bidir0[1][SIZE_BLOCK_16] = &oh_hevc_put_hevc_bi0_epel_h16_10_sse;
  mc_c->bidir1[1][SIZE_BLOCK_16] = &oh_hevc_put_hevc_bi1_epel_h16_10_sse;

  mc_c->unidir[2][SIZE_BLOCK_16] = &oh_hevc_put_hevc_uni_epel_v16_10_sse;
  mc_c->bidir0[2][SIZE_BLOCK_16] = &oh_hevc_put_hevc_bi0_epel_v16_10_sse;
  mc_c->bidir1[2][SIZE_BLOCK_16] = &oh_hevc_put_hevc_bi1_epel_v16_10_sse;

  mc_c->unidir[3][SIZE_BLOCK_16] = &oh_hevc_put_hevc_uni_epel_hv16_10_sse;
  mc_c->bidir0[3][SIZE_BLOCK_16] = &oh_hevc_put_hevc_bi0_epel_hv16_10_sse;
  mc_c->bidir1[3][SIZE_BLOCK_16] = &oh_hevc_put_hevc_bi1_epel_hv16_10_sse;

  mc_c->unidir[0][SIZE_BLOCK_32] = &oh_hevc_put_hevc_uni_pel_pixels32_10_sse;
  mc_c->bidir0[0][SIZE_BLOCK_32] = &oh_hevc_put_hevc_bi0_pel_pixels32_10_sse;
  mc_c->bidir1[0][SIZE_BLOCK_32] = &oh_hevc_put_hevc_bi1_pel_pixels32_10_sse;

  mc_c->unidir[1][SIZE_BLOCK_32] = &oh_hevc_put_hevc_uni_epel_h32_10_sse;
  mc_c->bidir0[1][SIZE_BLOCK_32] = &oh_hevc_put_hevc_bi0_epel_h32_10_sse;
  mc_c->bidir1[1][SIZE_BLOCK_32] = &oh_hevc_put_hevc_bi1_epel_h32_10_sse;

  mc_c->unidir[2][SIZE_BLOCK_32] = &oh_hevc_put_hevc_uni_epel_v32_10_sse;
  mc_c->bidir0[2][SIZE_BLOCK_32] = &oh_hevc_put_hevc_bi0_epel_v32_10_sse;
  mc_c->bidir1[2][SIZE_BLOCK_32] = &oh_hevc_put_hevc_bi1_epel_v32_10_sse;

  mc_c->unidir[3][SIZE_BLOCK_32] = &oh_hevc_put_hevc_uni_epel_hv32_10_sse;
  mc_c->bidir0[3][SIZE_BLOCK_32] = &oh_hevc_put_hevc_bi0_epel_hv32_10_sse;
  mc_c->bidir1[3][SIZE_BLOCK_32] = &oh_hevc_put_hevc_bi1_epel_hv32_10_sse;

  mc_c->unidir[0][SIZE_BLOCK_64] = &oh_hevc_put_hevc_uni_pel_pixels64_10_sse;
  mc_c->bidir0[0][SIZE_BLOCK_64] = &oh_hevc_put_hevc_bi0_pel_pixels64_10_sse;
  mc_c->bidir1[0][SIZE_BLOCK_64] = &oh_hevc_put_hevc_bi1_pel_pixels64_10_sse;

  mc_c->unidir[1][SIZE_BLOCK_64] = &oh_hevc_put_hevc_uni_epel_h64_10_sse;
  mc_c->bidir0[1][SIZE_BLOCK_64] = &oh_hevc_put_hevc_bi0_epel_h64_10_sse;
  mc_c->bidir1[1][SIZE_BLOCK_64] = &oh_hevc_put_hevc_bi1_epel_h64_10_sse;

  mc_c->unidir[2][SIZE_BLOCK_64] = &oh_hevc_put_hevc_uni_epel_v64_10_sse;
  mc_c->bidir0[2][SIZE_BLOCK_64] = &oh_hevc_put_hevc_bi0_epel_v64_10_sse;
  mc_c->bidir1[2][SIZE_BLOCK_64] = &oh_hevc_put_hevc_bi1_epel_v64_10_sse;

  mc_c->unidir[3][SIZE_BLOCK_64] = &oh_hevc_put_hevc_uni_epel_hv64_10_sse;
  mc_c->bidir0[3][SIZE_BLOCK_64] = &oh_hevc_put_hevc_bi0_epel_hv64_10_sse;
  mc_c->bidir1[3][SIZE_BLOCK_64] = &oh_hevc_put_hevc_bi1_epel_hv64_10_sse;
}
