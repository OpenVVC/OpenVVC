#include <stddef.h>
#include <stdint.h>
#include <string.h>

#include "ovutils.h"

#include "rcn_mc.h"
#include "rcn_structures.h"

#include "bitdepth.h"

#define MAX_PB_SIZE 128

#define EPEL_EXTRA_BEFORE 1
#define EPEL_EXTRA_AFTER 2
#define EPEL_EXTRA EPEL_EXTRA_BEFORE + EPEL_EXTRA_AFTER

#define QPEL_EXTRA_BEFORE 3
#define QPEL_EXTRA_AFTER 4
#define QPEL_EXTRA QPEL_EXTRA_BEFORE + QPEL_EXTRA_AFTER

#define MCP_FILTER_L(src, stride, filter)                                      \
        (filter[0] * src[x - 3 * stride] + filter[1] * src[x - 2 * stride] +   \
         filter[2] * src[x - stride] + filter[3] * src[x] +                    \
         filter[4] * src[x + stride] + filter[5] * src[x + 2 * stride] +       \
         filter[6] * src[x + 3 * stride] + filter[7] * src[x + 4 * stride])

#define MCP_FILTER_C(src, stride, filter)                                      \
        (filter[0] * src[x - stride] + filter[1] * src[x] +                    \
         filter[2] * src[x + stride] + filter[3] * src[x + 2 * stride])

#define BLN_FILTER_L(src, stride, filter)                                      \
        (filter[0] * src[x - 1 * stride] + filter[1] * src[x - 0 * stride])

static const int8_t ov_mcp_filters_c[31][4] =
{
    //{  0, 64,  0,  0 },
    { -1, 63,  2,  0 },
    { -2, 62,  4,  0 },
    { -2, 60,  7, -1 },
    { -2, 58, 10, -2 },
    { -3, 57, 12, -2 },
    { -4, 56, 14, -2 },
    { -4, 55, 15, -2 },
    { -4, 54, 16, -2 },
    { -5, 53, 18, -2 },
    { -6, 52, 20, -2 },
    { -6, 49, 24, -3 },
    { -6, 46, 28, -4 },
    { -5, 44, 29, -4 },
    { -4, 42, 30, -4 },
    { -4, 39, 33, -4 },
    { -4, 36, 36, -4 },
    { -4, 33, 39, -4 },
    { -4, 30, 42, -4 },
    { -4, 29, 44, -5 },
    { -4, 28, 46, -6 },
    { -3, 24, 49, -6 },
    { -2, 20, 52, -6 },
    { -2, 18, 53, -5 },
    { -2, 16, 54, -4 },
    { -2, 15, 55, -4 },
    { -2, 14, 56, -4 },
    { -2, 12, 57, -3 },
    { -2, 10, 58, -2 },
    { -1,  7, 60, -2 },
    {  0,  4, 62, -2 },
    {  0,  2, 63, -1 },
};

static const int8_t ov_mc_filters[16][8] =
{
    {   0, 1,  -3, 63,  4,  -2,  1,  0 },
    {  -1, 2,  -5, 62,  8,  -3,  1,  0 },
    {  -1, 3,  -8, 60, 13,  -4,  1,  0 },
    {  -1, 4, -10, 58, 17,  -5,  1,  0 },
    {  -1, 4, -11, 52, 26,  -8,  3, -1 },
    {  -1, 3,  -9, 47, 31, -10,  4, -1 },
    {  -1, 4, -11, 45, 34, -10,  4, -1 },
    {  -1, 4, -11, 40, 40, -11,  4, -1 },
    {  -1, 4, -10, 34, 45, -11,  4, -1 },
    {  -1, 4, -10, 31, 47,  -9,  3, -1 },
    {  -1, 3,  -8, 26, 52, -11,  4, -1 },
    {   0, 1,  -5, 17, 58, -10,  4, -1 },
    {   0, 1,  -4, 13, 60,  -8,  3, -1 },
    {   0, 1,  -3,  8, 62,  -5,  2, -1 },
    {   0, 1,  -2,  4, 63,  -3,  1,  0 },

    //Hpel for amvr
    {  0, 3, 9, 20, 20, 9, 3, 0 }
};

static const int8_t ov_mc_filters_4[15][8] =
{
    {  0, 1,  -3, 63,  4,  -2,  1,  0 },
    {  0, 1,  -5, 62,  8,  -3,  1,  0 },
    {  0, 2,  -8, 60, 13,  -4,  1,  0 },
    {  0, 3, -10, 58, 17,  -5,  1,  0 }, //1/4
    {  0, 3, -11, 52, 26,  -8,  2,  0 },
    {  0, 2,  -9, 47, 31, -10,  3,  0 },
    {  0, 3, -11, 45, 34, -10,  3,  0 },
    {  0, 3, -11, 40, 40, -11,  3,  0 }, //1/2
    {  0, 3, -10, 34, 45, -11,  3,  0 },
    {  0, 3, -10, 31, 47,  -9,  2,  0 },
    {  0, 2,  -8, 26, 52, -11,  3,  0 },
    {  0, 1,  -5, 17, 58, -10,  3,  0 }, //3/4
    {  0, 1,  -4, 13, 60,  -8,  2,  0 },
    {  0, 1,  -3,  8, 62,  -5,  1,  0 },
    {  0, 1,  -2,  4, 63,  -3,  1,  0 }
};

static const int8_t ov_bilinear_filters_4[15][2] =
{
  /*{ 16,  0, },*/
  { 15,  1, },
  { 14,  2, },
  { 13, 3, },
  { 12, 4, },
  { 11, 5, },
  { 10, 6, },
  { 9, 7, },
  { 8, 8, },
  { 7, 9, },
  { 6, 10, },
  { 5, 11, },
  { 4, 12, },
  { 3, 13, },
  { 2, 14, },
  { 1, 15, }
};

static void
put_vvc_pel_uni_pixels(uint16_t* _dst, ptrdiff_t _dststride,
                       const uint16_t* _src, ptrdiff_t _srcstride, int height,
                       intptr_t mx, intptr_t my, int width)
{
    int y;
    const uint16_t* src = _src;
    ptrdiff_t srcstride = _srcstride;
    uint16_t* dst = (uint16_t*)_dst;
    ptrdiff_t dststride = _dststride;

    for (y = 0; y < height; y++) {
        memcpy(dst, src, width * sizeof(uint16_t));
        src += srcstride;
        dst += dststride;
    }
}

/*FIXME It might actually be better to use one function for bi pred instead of
 * 2*/
static void
put_vvc_pel_pixels(int16_t* _dst, const uint16_t* _src, ptrdiff_t _srcstride,
                   int height, intptr_t mx, intptr_t my, int width)
{
    int x, y;
    const uint16_t* src = _src;

    int16_t* dst = (int16_t*)_dst;

    ptrdiff_t srcstride = _srcstride;

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; ++x) {
            dst[x] = src[x] << (14 - BITDEPTH);
        }
        src += srcstride;
        dst += MAX_PB_SIZE;
    }
}

static void
put_vvc_pel_bi_pixels(uint16_t* _dst, ptrdiff_t _dststride,
                      const uint16_t* _src0, ptrdiff_t _srcstride,
                      const int16_t* _src1, int height, intptr_t mx,
                      intptr_t my, int width)
{
    int x, y;
    const uint16_t* src0 = _src0;
    const int16_t* src1 = _src1;
    ptrdiff_t srcstride = _srcstride;
    uint16_t* dst = (uint16_t*)_dst;
    ptrdiff_t dststride = _dststride;
    int shift = 14 - BITDEPTH + 1;
    int offset = 1 << (shift - 1);

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; ++x) {
            dst[x] = ov_bdclip(((src0[x] << (14 - BITDEPTH)) + src1[x] + offset)
                               >> shift);
        }
        src0 += srcstride;
        src1 += MAX_PB_SIZE;
        dst += dststride;
    }
}

static void
put_vvc_qpel_uni_h(uint16_t* _dst, ptrdiff_t _dststride, const uint16_t* _src,
                   ptrdiff_t _srcstride, int height, intptr_t mx, intptr_t my,
                   int width)
{
    int x, y;
    const uint16_t* src = _src;
    ptrdiff_t srcstride = _srcstride;
    uint16_t* dst = (uint16_t*)_dst;
    ptrdiff_t dststride = _dststride;
    const int8_t* filter = width == 4 && height == 4 ? ov_mc_filters_4[mx - 1] : ov_mc_filters[mx - 1];
    int shift = 14 - BITDEPTH;
    int offset = 1 << (shift - 1);

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = ov_bdclip(((MCP_FILTER_L(src, 1, filter) >> (BITDEPTH - 8))
                                + offset) >> shift);
        }
        src += srcstride;
        dst += dststride;
    }
}

static void
put_vvc_qpel_uni_v(uint16_t* _dst, ptrdiff_t _dststride, const uint16_t* _src,
                   ptrdiff_t _srcstride, int height, intptr_t mx, intptr_t my,
                   int width)
{
    int x, y;
    const uint16_t* src = _src;
    ptrdiff_t srcstride = _srcstride;
    uint16_t* dst = (uint16_t*)_dst;
    ptrdiff_t dststride = _dststride;
    const int8_t* filter = width == 4 && height == 4 ? ov_mc_filters_4[my - 1] : ov_mc_filters[my - 1];
    int shift = 14 - BITDEPTH;
    int offset = 1 << (shift - 1);

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = ov_bdclip(((MCP_FILTER_L(src, srcstride, filter) >> (BITDEPTH - 8))
                                + offset) >> shift);
        }
        src += srcstride;
        dst += dststride;
    }
}

static void
put_vvc_qpel_uni_hv(uint16_t* _dst, ptrdiff_t _dststride, const uint16_t* _src,
                    ptrdiff_t _srcstride, int height, intptr_t mx, intptr_t my,
                    int width)
{
    int x, y;
    const int8_t* filter;
    const uint16_t* src = _src;
    ptrdiff_t srcstride = _srcstride;
    uint16_t* dst = (uint16_t*)_dst;
    ptrdiff_t dststride = _dststride;
    int16_t tmp_array[(MAX_PB_SIZE + QPEL_EXTRA) * MAX_PB_SIZE];
    int16_t* tmp = tmp_array;
    int shift = 14 - BITDEPTH;
    int offset = 1 << (shift - 1);

    src -= QPEL_EXTRA_BEFORE * srcstride;
    filter = width == 4 && height == 4 ? ov_mc_filters_4[mx - 1] : ov_mc_filters[mx - 1];
    for (y = 0; y < height + QPEL_EXTRA; y++) {
        for (x = 0; x < width; x++) {
            tmp[x] = MCP_FILTER_L(src, 1, filter) >> (BITDEPTH - 8);
        }
        src += srcstride;
        tmp += MAX_PB_SIZE;
    }

    tmp = tmp_array + QPEL_EXTRA_BEFORE * MAX_PB_SIZE;
    filter = width == 4 && height == 4 ? ov_mc_filters_4[my - 1] : ov_mc_filters[my - 1];

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = ov_bdclip(((MCP_FILTER_L(tmp, MAX_PB_SIZE, filter) >> 6)
                                + offset) >> shift);
        }
        tmp += MAX_PB_SIZE;
        dst += dststride;
    }
}

static void
put_vvc_pel_bilinear_pixels(uint16_t* _dst, ptrdiff_t _dststride,
                            const uint16_t* _src, ptrdiff_t _srcstride, int height,
                            intptr_t mx, intptr_t my, int width)
{
    int y;
    const uint16_t* src = _src;
    ptrdiff_t srcstride = _srcstride;
    uint16_t* dst = (uint16_t*)_dst;
    ptrdiff_t dststride = _dststride;

    for (y = 0; y < height; y++) {
        int x;
        for (x = 0; x < width; ++x) {
            #if (BITDEPTH - 10) > 0
            dst [x] = (src[x] + (1 << (BITDEPTH - 10)) >> (BITDEPTH - 10);
            #else
            dst [x] = src[x] << (10 - BITDEPTH);
            #endif
        }
        src += srcstride;
        dst += dststride;
    }
}

static void
put_vvc_qpel_bilinear_h(uint16_t* _dst, ptrdiff_t _dststride, const uint16_t* _src,
                        ptrdiff_t _srcstride, int height, intptr_t mx, intptr_t my,
                        int width)
{
    int x, y;
    const uint16_t* src = _src;
    ptrdiff_t srcstride = _srcstride;
    uint16_t* dst = (uint16_t*)_dst;
    ptrdiff_t dststride = _dststride;
    const int8_t* filter = ov_bilinear_filters_4[mx - 1];
    int shift = 4 - (10 - BITDEPTH);
    int offset = 1 << (shift - 1);

    src += 1;
    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = (BLN_FILTER_L(src, 1, filter) + offset) >> shift;
        }
        src += srcstride;
        dst += dststride;
    }
}

static void
put_vvc_qpel_bilinear_v(uint16_t* _dst, ptrdiff_t _dststride, const uint16_t* _src,
                        ptrdiff_t _srcstride, int height, intptr_t mx, intptr_t my,
                        int width)
{
    int x, y;
    const uint16_t* src = _src;
    ptrdiff_t srcstride = _srcstride;
    uint16_t* dst = (uint16_t*)_dst;
    ptrdiff_t dststride = _dststride;
    const int8_t* filter = ov_bilinear_filters_4[my - 1];
    int shift = 4 - (10 - BITDEPTH);
    int offset = 1 << (shift - 1);

    src += srcstride;
    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = (BLN_FILTER_L(src, srcstride, filter) + offset) >> shift;
        }
        src += srcstride;
        dst += dststride;
    }
}

static void
put_vvc_qpel_bilinear_hv(uint16_t* _dst, ptrdiff_t _dststride, const uint16_t* _src,
                         ptrdiff_t _srcstride, int height, intptr_t mx, intptr_t my,
                         int width)
{
    int x, y;
    const int8_t* filter;
    const uint16_t* src = _src;
    ptrdiff_t srcstride = _srcstride;
    uint16_t* dst = (uint16_t*)_dst;
    ptrdiff_t dststride = _dststride;
    int16_t tmp_array[(MAX_PB_SIZE + QPEL_EXTRA) * MAX_PB_SIZE];
    int16_t* tmp = tmp_array;

    int shift = 4 - (10 - BITDEPTH);
    int offset = 1 << (shift - 1);

    filter = ov_bilinear_filters_4[mx - 1];

    for (y = 0; y < height + 2; y++) {
        for (x = 0; x < width + 2; x++) {
            tmp[x] = (BLN_FILTER_L(src, 1, filter) + offset) >> shift;
        }
        src += srcstride;
        tmp += MAX_PB_SIZE;
    }

    tmp = tmp_array + MAX_PB_SIZE + 1;
    filter = ov_bilinear_filters_4[my - 1];

    shift = 4;
    offset = 1 << (shift - 1);
    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = (BLN_FILTER_L(tmp, MAX_PB_SIZE, filter) + offset) >> shift;
        }
        tmp += MAX_PB_SIZE;
        dst += dststride;
    }
}

static void
put_vvc_qpel_h(int16_t* _dst, const uint16_t* _src, ptrdiff_t _srcstride,
               int height, intptr_t mx, intptr_t my, int width)
{
    int x, y;
    const uint16_t* src = _src;

    int16_t* dst = (int16_t*)_dst;

    ptrdiff_t srcstride = _srcstride;

    const int8_t* filter;
    filter = width == 4 && height == 4 ? ov_mc_filters_4[mx - 1] : ov_mc_filters[mx - 1];

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = MCP_FILTER_L(src, 1, filter) >> (BITDEPTH - 8);
        }
        src += srcstride;
        dst += MAX_PB_SIZE;
    }
}

static void
put_vvc_qpel_v(int16_t* _dst, const uint16_t* _src, ptrdiff_t _srcstride,
               int height, intptr_t mx, intptr_t my, int width)
{
    int x, y;
    const uint16_t* src = (uint16_t*)_src;

    int16_t* dst = (int16_t*)_dst;

    ptrdiff_t srcstride = _srcstride;

    const int8_t* filter = ov_mc_filters[my - 1];
    filter = width == 4 && height == 4 ? ov_mc_filters_4[my - 1] : ov_mc_filters[my - 1];

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = MCP_FILTER_L(src, srcstride, filter) >>
                (BITDEPTH - 8);
        }
        src += srcstride;
        dst += MAX_PB_SIZE;
    }
}

static void
put_vvc_qpel_hv(int16_t* _dst, const uint16_t* _src, ptrdiff_t _srcstride,
                int height, intptr_t mx, intptr_t my, int width)
{
    int x, y;

    const uint16_t* src = (uint16_t*)_src;

    int16_t* dst = (int16_t*)_dst;

    ptrdiff_t srcstride = _srcstride;

    int16_t tmp_array[(MAX_PB_SIZE + QPEL_EXTRA) * MAX_PB_SIZE];
    int16_t* tmp = tmp_array;

    const int8_t* filter = ov_mc_filters[mx - 1];
    filter = width == 4 && height == 4 ? ov_mc_filters_4[mx - 1] : ov_mc_filters[mx - 1];

    src -= QPEL_EXTRA_BEFORE * srcstride;

    for (y = 0; y < height + QPEL_EXTRA; y++) {
        for (x = 0; x < width; x++) {
            tmp[x] = MCP_FILTER_L(src, 1, filter) >> (BITDEPTH - 8);
        }
        src += srcstride;
        tmp += MAX_PB_SIZE;
    }

    tmp = tmp_array + QPEL_EXTRA_BEFORE * MAX_PB_SIZE;
    filter = ov_mc_filters[my - 1];
    filter = width == 4 && height == 4 ? ov_mc_filters_4[my - 1] : ov_mc_filters[my - 1];

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = MCP_FILTER_L(tmp, MAX_PB_SIZE, filter) >> 6;
        }
        tmp += MAX_PB_SIZE;
        dst += MAX_PB_SIZE;
    }
}

static void
put_vvc_qpel_bi_h(uint16_t* _dst, ptrdiff_t _dststride, const uint16_t* _src0,
                  ptrdiff_t _srcstride, const int16_t* _src1, int height,
                  intptr_t mx, intptr_t my, int width)
{
    int x, y;
    const uint16_t* src0 = _src0;
    const int16_t* src1 = _src1;

    uint16_t* dst = (uint16_t*)_dst;

    ptrdiff_t srcstride = _srcstride;
    ptrdiff_t dststride = _dststride;

    const int8_t* filter = ov_mc_filters[mx - 1];
    filter = width == 4 && height == 4 ? ov_mc_filters_4[mx - 1] : ov_mc_filters[mx - 1];

    int shift = 14 + 1 - BITDEPTH;
    int offset = 1 << (shift - 1);

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = ov_bdclip(((MCP_FILTER_L(src0, 1, filter) >> (BITDEPTH - 8))
                                + src1[x] + offset) >> shift);
        }
        src0 += srcstride;
        src1 += MAX_PB_SIZE;
        dst += dststride;
    }
}

static void
put_vvc_qpel_bi_v(uint16_t* _dst, ptrdiff_t _dststride, const uint16_t* _src0,
                  ptrdiff_t _srcstride, const int16_t* _src1, int height,
                  intptr_t mx, intptr_t my, int width)
{
    int x, y;
    const uint16_t* src0 = _src0;
    const int16_t* src1 = _src1;

    uint16_t* dst = (uint16_t*)_dst;

    ptrdiff_t srcstride = _srcstride;
    ptrdiff_t dststride = _dststride;

    const int8_t* filter = ov_mc_filters[my - 1];
    filter = width == 4 && height == 4 ? ov_mc_filters_4[my - 1] : ov_mc_filters[my - 1];

    int shift = 14 + 1 - BITDEPTH;
    int offset = 1 << (shift - 1);

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = ov_bdclip(((MCP_FILTER_L(src0, srcstride, filter) >> (BITDEPTH - 8))
                                + src1[x] + offset) >> shift);
        }
        src0 += srcstride;
        src1 += MAX_PB_SIZE;
        dst += dststride;
    }
}

static void
put_vvc_qpel_bi_hv(uint16_t* _dst, ptrdiff_t _dststride, const uint16_t* _src0,
                   ptrdiff_t _srcstride, const int16_t* _src1, int height,
                   intptr_t mx, intptr_t my, int width)
{
    int x, y;

    const uint16_t* src0 = _src0;
    const int16_t* src1 = _src1;
    uint16_t* dst = (uint16_t*)_dst;

    ptrdiff_t srcstride = _srcstride;
    ptrdiff_t dststride = _dststride;

    int16_t tmp_array[(MAX_PB_SIZE + QPEL_EXTRA) * MAX_PB_SIZE];
    int16_t* tmp = tmp_array;

    const int8_t* filter = ov_mc_filters[mx - 1];
    filter = width == 4 && height == 4 ? ov_mc_filters_4[mx - 1] : ov_mc_filters[mx - 1];

    int shift = 14 + 1 - BITDEPTH;
    int offset = 1 << (shift - 1);

    src0 -= QPEL_EXTRA_BEFORE * srcstride;

    for (y = 0; y < height + QPEL_EXTRA; y++) {
        for (x = 0; x < width; x++) {
            tmp[x] = MCP_FILTER_L(src0, 1, filter) >> (BITDEPTH - 8);
        }
        src0 += srcstride;
        tmp += MAX_PB_SIZE;
    }

    tmp = tmp_array + QPEL_EXTRA_BEFORE * MAX_PB_SIZE;
    filter = ov_mc_filters[my - 1];
    filter = width == 4 && height == 4 ? ov_mc_filters_4[my - 1] : ov_mc_filters[my - 1];

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = ov_bdclip(((MCP_FILTER_L(tmp, MAX_PB_SIZE, filter) >> 6)
                                + src1[x] + offset) >> shift);
        }
        tmp += MAX_PB_SIZE;
        src1 += MAX_PB_SIZE;
        dst += dststride;
    }
}

static void
put_vvc_epel_uni_h(uint16_t* _dst, ptrdiff_t _dststride, const uint16_t* _src,
                   ptrdiff_t _srcstride, int height, intptr_t mx, intptr_t my,
                   int width)
{
    int x, y;
    const uint16_t* src = _src;
    ptrdiff_t srcstride = _srcstride;
    uint16_t* dst = (uint16_t*)_dst;
    ptrdiff_t dststride = _dststride;
    const int8_t* filter = ov_mcp_filters_c[mx - 1];
    int shift = 14 - BITDEPTH;
    int offset = 1 << (shift - 1);

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = ov_bdclip(((MCP_FILTER_C(src, 1, filter) >> (BITDEPTH - 8))
                                + offset) >> shift);
        }
        src += srcstride;
        dst += dststride;
    }
}

static void
put_vvc_epel_uni_v(uint16_t* _dst, ptrdiff_t _dststride, const uint16_t* _src,
                   ptrdiff_t _srcstride, int height, intptr_t mx, intptr_t my,
                   int width)
{
    int x, y;
    const uint16_t* src = _src;
    ptrdiff_t srcstride = _srcstride;
    uint16_t* dst = (uint16_t*)_dst;
    ptrdiff_t dststride = _dststride;
    const int8_t* filter = ov_mcp_filters_c[my - 1];
    int shift = 14 - BITDEPTH;
    int offset = 1 << (shift - 1);

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = ov_bdclip(((MCP_FILTER_C(src, srcstride, filter) >> (BITDEPTH - 8))
                                + offset) >> shift);
        }
        src += srcstride;
        dst += dststride;
    }
}

static void
put_vvc_epel_uni_hv(uint16_t* _dst, ptrdiff_t _dststride, const uint16_t* _src,
                    ptrdiff_t _srcstride, int height, intptr_t mx, intptr_t my,
                    int width)
{
    int x, y;
    const uint16_t* src = _src;
    ptrdiff_t srcstride = _srcstride;
    uint16_t* dst = (uint16_t*)_dst;
    ptrdiff_t dststride = _dststride;
    const int8_t* filter = ov_mcp_filters_c[mx - 1];
    int16_t tmp_array[(MAX_PB_SIZE + EPEL_EXTRA) * MAX_PB_SIZE];
    int16_t* tmp = tmp_array;
    int shift = 14 - BITDEPTH;
    int offset = 1 << (shift - 1);

    src -= EPEL_EXTRA_BEFORE * srcstride;

    for (y = 0; y < height + EPEL_EXTRA; y++) {
        for (x = 0; x < width; x++) {
            tmp[x] = MCP_FILTER_C(src, 1, filter) >> (BITDEPTH - 8);
        }
        src += srcstride;
        tmp += MAX_PB_SIZE;
    }

    tmp = tmp_array + EPEL_EXTRA_BEFORE * MAX_PB_SIZE;
    filter = ov_mcp_filters_c[my - 1];

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = ov_bdclip(((MCP_FILTER_C(tmp, MAX_PB_SIZE, filter) >> 6)
                                + offset) >> shift);
        }
        tmp += MAX_PB_SIZE;
        dst += dststride;
    }
}

static void
put_vvc_epel_h(int16_t* _dst, const uint16_t* _src, ptrdiff_t _srcstride,
               int height, intptr_t mx, intptr_t my, int width)
{
    int x, y;

    const uint16_t* src = _src;

    int16_t* dst = (int16_t*)_dst;

    ptrdiff_t srcstride = _srcstride;

    const int8_t* filter = ov_mcp_filters_c[mx - 1];

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = MCP_FILTER_C(src, 1, filter) >> (BITDEPTH - 8);
        }
        src += srcstride;
        dst += MAX_PB_SIZE;
    }
}

static void
put_vvc_epel_v(int16_t* _dst, const uint16_t* _src, ptrdiff_t _srcstride,
               int height, intptr_t mx, intptr_t my, int width)
{
    int x, y;

    const uint16_t* src = _src;

    int16_t* dst = (int16_t*)_dst;

    ptrdiff_t srcstride = _srcstride;

    const int8_t* filter = ov_mcp_filters_c[my - 1];

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = MCP_FILTER_C(src, srcstride, filter) >> (BITDEPTH - 8);
        }
        src += srcstride;
        dst += MAX_PB_SIZE;
    }
}

static void
put_vvc_epel_hv(int16_t* _dst, const uint16_t* _src, ptrdiff_t _srcstride,
                int height, intptr_t mx, intptr_t my, int width)
{
    int x, y;
    const uint16_t* src = _src;

    int16_t* dst = (int16_t*)_dst;

    ptrdiff_t srcstride = _srcstride;

    const int8_t* filter = ov_mcp_filters_c[mx - 1];

    int16_t tmp_array[(MAX_PB_SIZE + EPEL_EXTRA) * MAX_PB_SIZE];
    int16_t* tmp = tmp_array;

    src -= EPEL_EXTRA_BEFORE * srcstride;

    for (y = 0; y < height + EPEL_EXTRA; y++) {
        for (x = 0; x < width; x++) {
            tmp[x] = MCP_FILTER_C(src, 1, filter) >> (BITDEPTH - 8);
        }
        src += srcstride;
        tmp += MAX_PB_SIZE;
    }

    tmp = tmp_array + EPEL_EXTRA_BEFORE * MAX_PB_SIZE;
    filter = ov_mcp_filters_c[my - 1];

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = MCP_FILTER_C(tmp, MAX_PB_SIZE, filter) >> 6;
        }
        tmp += MAX_PB_SIZE;
        dst += MAX_PB_SIZE;
    }
}

static void
put_vvc_epel_bi_h(uint16_t* _dst, ptrdiff_t _dststride, const uint16_t* _src0,
                  ptrdiff_t _srcstride, const int16_t* _src1, int height,
                  intptr_t mx, intptr_t my, int width)
{
    int x, y;
    const uint16_t* src0 = _src0;
    const int16_t* src1 = (int16_t*)_src1;

    uint16_t* dst = (uint16_t*)_dst;

    ptrdiff_t srcstride = _srcstride;
    ptrdiff_t dststride = _dststride;

    const int8_t* filter = ov_mcp_filters_c[mx - 1];

    int shift = 14 + 1 - BITDEPTH;
    int offset = 1 << (shift - 1);

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = ov_bdclip(((MCP_FILTER_C(src0, 1, filter) >> (BITDEPTH - 8))
                                + src1[x] + offset) >> shift);
        }
        src0 += srcstride;
        src1 += MAX_PB_SIZE;
        dst += dststride;
    }
}


static void
put_vvc_epel_bi_v(uint16_t* _dst, ptrdiff_t _dststride, const uint16_t* _src0,
                  ptrdiff_t _srcstride, const int16_t* _src1, int height,
                  intptr_t mx, intptr_t my, int width)
{
    int x, y;
    const uint16_t* src0 = _src0;
    const int16_t* src1 = _src1;

    uint16_t* dst = (uint16_t*)_dst;

    ptrdiff_t srcstride = _srcstride;
    ptrdiff_t dststride = _dststride;

    const int8_t* filter = ov_mcp_filters_c[my - 1];

    int shift = 14 + 1 - BITDEPTH;
    int offset = 1 << (shift - 1);

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = ov_bdclip(((MCP_FILTER_C(src0, srcstride, filter) >> (BITDEPTH - 8))
                                + src1[x] + offset) >> shift);
        }
        src0 += srcstride;
        src1 += MAX_PB_SIZE;
        dst += dststride;
    }
}

static void
put_vvc_epel_bi_hv(uint16_t* _dst, ptrdiff_t _dststride, const uint16_t* _src0,
                   ptrdiff_t _srcstride, const int16_t* _src1, int height,
                   intptr_t mx, intptr_t my, int width)
{
    int x, y;
    const uint16_t* src0 = _src0;
    const int16_t* src1 = _src1;

    uint16_t* dst = (uint16_t*)_dst;

    ptrdiff_t srcstride = _srcstride;
    ptrdiff_t dststride = _dststride;

    int16_t tmp_array[(MAX_PB_SIZE + EPEL_EXTRA) * MAX_PB_SIZE];
    int16_t* tmp = tmp_array;

    const int8_t* filter = ov_mcp_filters_c[mx - 1];

    int shift = 14 + 1 - BITDEPTH;
    int offset = 1 << (shift - 1);

    src0 -= EPEL_EXTRA_BEFORE * srcstride;

    for (y = 0; y < height + EPEL_EXTRA; y++) {
        for (x = 0; x < width; x++) {
            tmp[x] = MCP_FILTER_C(src0, 1, filter) >> (BITDEPTH - 8);
        }
        src0 += srcstride;
        tmp += MAX_PB_SIZE;
    }

    tmp = tmp_array + EPEL_EXTRA_BEFORE * MAX_PB_SIZE;
    filter = ov_mcp_filters_c[my - 1];

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = ov_bdclip(((MCP_FILTER_C(tmp, MAX_PB_SIZE, filter) >> 6)
                                + src1[x] + offset) >> shift);
        }
        tmp += MAX_PB_SIZE;
        src1 += MAX_PB_SIZE;
        dst += dststride;
    }
}


static void
put_weighted_epel_bi_h(uint8_t* _dst, ptrdiff_t _dststride, uint8_t* _src, ptrdiff_t _srcstride,
                  int16_t* src2, ptrdiff_t src2stride, int height, int denom,
                  int wx0, int wx1, intptr_t mx, intptr_t my, int width)
{
    int x, y;
    const uint16_t* src = (uint16_t*) _src;
    uint16_t* dst = (uint16_t*)_dst;

    ptrdiff_t srcstride = _srcstride >> 1;
    ptrdiff_t dststride = _dststride >> 1;

    const int8_t* filter = ov_mcp_filters_c[mx - 1];

    denom = floor_log2(wx0 + wx1);
    int shift = 14 + denom -BITDEPTH;
    int offset = 1 << (shift - 1);

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = ov_bdclip(((MCP_FILTER_C(src, 1, filter) >> (BITDEPTH - 8)) * wx1
                                + src2[x] * wx0 + offset) >> shift);
        }
        src += srcstride;
        src2 += MAX_PB_SIZE;
        dst += dststride;
    }
}

static void
put_weighted_epel_bi_v(uint8_t* _dst, ptrdiff_t _dststride, uint8_t* _src, ptrdiff_t _srcstride,
                  int16_t* src2, ptrdiff_t src2stride, int height, int denom,
                  int wx0, int wx1, intptr_t mx, intptr_t my, int width)
{
    int x, y;
    const uint16_t* src = (uint16_t*) _src;


    uint16_t* dst = (uint16_t*)_dst;

    ptrdiff_t srcstride = _srcstride >> 1;
    ptrdiff_t dststride = _dststride >> 1;

    const int8_t* filter = ov_mcp_filters_c[my - 1];

    denom = floor_log2(wx0 + wx1);
    int shift = 14 + denom -BITDEPTH;
    int offset = 1 << (shift - 1);

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = ov_bdclip(((MCP_FILTER_C(src, srcstride, filter) >> (BITDEPTH - 8)) * wx1
                                + src2[x] * wx0 + offset) >> shift);
        }
        src += srcstride;
        src2 += MAX_PB_SIZE;
        dst += dststride;
    }
}

static void
put_weighted_epel_bi_hv(uint8_t* _dst, ptrdiff_t _dststride, uint8_t* _src, ptrdiff_t _srcstride,
                   int16_t* src2, ptrdiff_t src2stride, int height, int denom,
                   int wx0, int wx1, intptr_t mx, intptr_t my, int width)
{
    int x, y;
    const uint16_t* src = (uint16_t*) _src;


    uint16_t* dst = (uint16_t*)_dst;

    ptrdiff_t srcstride = _srcstride >> 1;
    ptrdiff_t dststride = _dststride >> 1;

    int16_t tmp_array[(MAX_PB_SIZE + EPEL_EXTRA) * MAX_PB_SIZE];
    int16_t* tmp = tmp_array;

    const int8_t* filter = ov_mcp_filters_c[mx - 1];

    denom = floor_log2(wx0 + wx1);
    int shift = 14 + denom -BITDEPTH;
    int offset = 1 << (shift - 1);

    src -= EPEL_EXTRA_BEFORE * srcstride;

    for (y = 0; y < height + EPEL_EXTRA; y++) {
        for (x = 0; x < width; x++) {
            tmp[x] = MCP_FILTER_C(src, 1, filter) >> (BITDEPTH - 8);
        }
        src += srcstride;
        tmp += MAX_PB_SIZE;
    }

    tmp = tmp_array + EPEL_EXTRA_BEFORE * MAX_PB_SIZE;
    filter = ov_mcp_filters_c[my - 1];

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = ov_bdclip(((MCP_FILTER_C(tmp, MAX_PB_SIZE, filter) >> 6) * wx1
                                + src2[x] * wx0 + offset) >> shift);
        }
        tmp += MAX_PB_SIZE;
        src2 += MAX_PB_SIZE;
        dst += dststride;
    }
}

static void
put_weighted_pel_bi_pixels(uint8_t* _dst, ptrdiff_t _dststride, uint8_t* _src, ptrdiff_t _srcstride,
                  int16_t* src2, ptrdiff_t src2stride, int height, int denom,
                  int wx0, int wx1, intptr_t mx, intptr_t my, int width)
{
    int x, y;
    const uint16_t* src = (uint16_t*) _src;

    ptrdiff_t srcstride = _srcstride >> 1;
    uint16_t* dst = (uint16_t*)_dst;
    ptrdiff_t dststride = _dststride >> 1;
    int shift = 14 - BITDEPTH + 3;
    int offset = (1 << (shift - 1)) ;
    for (y = 0; y < height; y++) {
        for (x = 0; x < width; ++x) {
            dst[x] = ov_bdclip(((src2[x] * wx0 + ((src[x] * wx1) << (14 - BITDEPTH))
                                 + offset)  >> shift));
        }
        src2 += MAX_PB_SIZE;
        src += srcstride;
        dst += dststride;
    }
}

static void
put_weighted_qpel_bi_h(uint8_t* _dst, ptrdiff_t _dststride, uint8_t* _src, ptrdiff_t _srcstride,
                  int16_t* src2, ptrdiff_t src2stride, int height, int denom,
                  int wx0, int wx1, intptr_t mx, intptr_t my, int width)
{
    int x, y;
    const uint16_t* src = (uint16_t*) _src;


    uint16_t* dst = (uint16_t*)_dst;

    ptrdiff_t srcstride = _srcstride >> 1;
    ptrdiff_t dststride = _dststride >> 1;

    const int8_t* filter;
    filter = width == 4 && height == 4 ? ov_mc_filters_4[mx - 1] : ov_mc_filters[mx - 1];

    denom = floor_log2(wx0 + wx1);
    int shift = 14 + denom -BITDEPTH;
    int offset = 1 << (shift - 1);

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = ov_bdclip(((MCP_FILTER_L(src, 1, filter) >> (BITDEPTH - 8)) * wx1
                                + src2[x] * wx0 + offset) >> shift);
        }
        src += srcstride;
        src2 += MAX_PB_SIZE;
        dst += dststride;
    }
}


static void
put_weighted_qpel_bi_v(uint8_t* _dst, ptrdiff_t _dststride, uint8_t* _src, ptrdiff_t _srcstride,
                  int16_t* src2, ptrdiff_t src2stride, int height, int denom,
                  int wx0, int wx1, intptr_t mx, intptr_t my, int width)
{
    int x, y;
    const uint16_t* src = (uint16_t*) _src;

    uint16_t* dst = (uint16_t*)_dst;

    ptrdiff_t srcstride = _srcstride >> 1;
    ptrdiff_t dststride = _dststride >> 1;

    const int8_t* filter;
    filter = width == 4 && height == 4 ? ov_mc_filters_4[my - 1] : ov_mc_filters[my - 1];

    denom = floor_log2(wx0 + wx1);
    int shift = 14 + denom -BITDEPTH;
    int offset = 1 << (shift - 1);

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = ov_bdclip(((MCP_FILTER_L(src, srcstride, filter) >> (BITDEPTH - 8)) * wx1
                                + src2[x] * wx0 + offset) >> shift);
        }
        src += srcstride;
        src2 += MAX_PB_SIZE;
        dst += dststride;
    }
}

static void
put_weighted_qpel_bi_hv(uint8_t* _dst, ptrdiff_t _dststride, uint8_t* _src, ptrdiff_t _srcstride,
                   int16_t* src2, ptrdiff_t src2stride, int height, int denom,
                   int wx0, int wx1, intptr_t mx, intptr_t my, int width)
{
    int x, y;

    const uint16_t* src = (uint16_t*) _src;

    uint16_t* dst = (uint16_t*)_dst;

    ptrdiff_t srcstride = _srcstride >> 1;
    ptrdiff_t dststride = _dststride >> 1;

    int16_t tmp_array[(MAX_PB_SIZE + QPEL_EXTRA) * MAX_PB_SIZE];
    int16_t* tmp = tmp_array;

    const int8_t* filter;
    filter = width == 4 && height == 4 ? ov_mc_filters_4[mx - 1] : ov_mc_filters[mx - 1];

    denom = floor_log2(wx0 + wx1);
    int shift = 14 + denom -BITDEPTH;
    int offset = 1 << (shift - 1);

    src -= QPEL_EXTRA_BEFORE * srcstride;

    for (y = 0; y < height + QPEL_EXTRA; y++) {
        for (x = 0; x < width; x++) {
            tmp[x] = MCP_FILTER_L(src, 1, filter) >> (BITDEPTH - 8);
        }
        src += srcstride;
        tmp += MAX_PB_SIZE;
    }

    tmp = tmp_array + QPEL_EXTRA_BEFORE * MAX_PB_SIZE;
    filter = width == 4 && height == 4 ? ov_mc_filters_4[my - 1] : ov_mc_filters[my - 1];

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = ov_bdclip(((MCP_FILTER_L(tmp, MAX_PB_SIZE, filter) >> 6) * wx1
                                + src2[x] * wx0 + offset) >> shift);
        }
        tmp += MAX_PB_SIZE;
        src2 += MAX_PB_SIZE;
        dst += dststride;
    }
}


static void
put_weighted_ciip_pixels(uint16_t* dst, int dststride,
                      const uint16_t* src_intra, const uint16_t* src_inter, int srcstride,
                      int width, int height, int wt)
{
    int x, y;
    int shift  = 2;
    int offset = (1 << (shift - 1));
    for (y = 0; y < height; y++) {
        for (x = 0; x < width; ++x) {
            dst[x] = ov_bdclip(((src_intra[x] * wt + src_inter[x] * (4 - wt)) + offset) >> shift);
        }
        src_intra += srcstride;
        src_inter += srcstride;
        dst += dststride;
    }
}


static void
put_weighted_gpm_bi_pixels(uint16_t* _dst, int _dststride, const int16_t* _src0,
                           int _srcstride, const int16_t* _src1, int height,
                           intptr_t mx, intptr_t my, int width,
                           int step_x, int step_y, int16_t* weight)
{
    int x, y;
    const int16_t* src0 = _src0;
    const int16_t* src1 = _src1;
    ptrdiff_t srcstride = _srcstride;
    uint16_t* dst = (uint16_t*)_dst;
    ptrdiff_t dststride = _dststride;
    int shift = 14 - BITDEPTH + 3;
    int offset = (1 << (shift - 1)) ;
    for (y = 0; y < height; y++) {
        for (x = 0; x < width; ++x) {
            int w1 = weight[0];
            int w0 = 8 - w1;
            dst[x] = ov_bdclip(((src1[x] * w1 + src0[x] * w0 + offset) >> shift));
            weight += step_x;
        }
        src1 += srcstride;
        src0 += srcstride;
        dst += dststride;
        weight += step_y;
    }
}

void
BD_DECL(rcn_init_mc_functions)(struct RCNFunctions *const rcn_funcs)
{
    struct MCFunctions *const mc_l = &rcn_funcs->mc_l;
    struct MCFunctions *const mc_c = &rcn_funcs->mc_c;

    int i;

    for (i = 0; i < 8; ++i) {

        /* Luma functions */
        mc_l->unidir[0][i] = &put_vvc_pel_uni_pixels;
        mc_l->unidir[1][i] = &put_vvc_qpel_uni_h;
        mc_l->unidir[2][i] = &put_vvc_qpel_uni_v;
        mc_l->unidir[3][i] = &put_vvc_qpel_uni_hv;

        mc_l->bidir0[0][i] = &put_vvc_pel_pixels;
        mc_l->bidir0[1][i] = &put_vvc_qpel_h;
        mc_l->bidir0[2][i] = &put_vvc_qpel_v;
        mc_l->bidir0[3][i] = &put_vvc_qpel_hv;

        mc_l->bidir1[0][i] = &put_vvc_pel_bi_pixels;
        mc_l->bidir1[1][i] = &put_vvc_qpel_bi_h;
        mc_l->bidir1[2][i] = &put_vvc_qpel_bi_v;
        mc_l->bidir1[3][i] = &put_vvc_qpel_bi_hv;

        mc_l->bidir_w[0][i] = &put_weighted_pel_bi_pixels;
        mc_l->bidir_w[1][i] = &put_weighted_qpel_bi_h;
        mc_l->bidir_w[2][i] = &put_weighted_qpel_bi_v;
        mc_l->bidir_w[3][i] = &put_weighted_qpel_bi_hv;

        mc_l->bilinear[0][i] = &put_vvc_pel_bilinear_pixels;
        mc_l->bilinear[1][i] = &put_vvc_qpel_bilinear_h;
        mc_l->bilinear[2][i] = &put_vvc_qpel_bilinear_v;
        mc_l->bilinear[3][i] = &put_vvc_qpel_bilinear_hv;

        /* Chroma functions */
        mc_c->unidir[0][i] = &put_vvc_pel_uni_pixels;
        mc_c->unidir[1][i] = &put_vvc_epel_uni_h;
        mc_c->unidir[2][i] = &put_vvc_epel_uni_v;
        mc_c->unidir[3][i] = &put_vvc_epel_uni_hv;

        mc_c->bidir0[0][i] = &put_vvc_pel_pixels;
        mc_c->bidir0[1][i] = &put_vvc_epel_h;
        mc_c->bidir0[2][i] = &put_vvc_epel_v;
        mc_c->bidir0[3][i] = &put_vvc_epel_hv;

        mc_c->bidir1[0][i] = &put_vvc_pel_bi_pixels;
        mc_c->bidir1[1][i] = &put_vvc_epel_bi_h;
        mc_c->bidir1[2][i] = &put_vvc_epel_bi_v;
        mc_c->bidir1[3][i] = &put_vvc_epel_bi_hv;

        mc_c->bidir_w[0][i] = &put_weighted_pel_bi_pixels;
        mc_c->bidir_w[1][i] = &put_weighted_epel_bi_h;
        mc_c->bidir_w[2][i] = &put_weighted_epel_bi_v;
        mc_c->bidir_w[3][i] = &put_weighted_epel_bi_hv;
    }
    mc_l->gpm_weighted = &put_weighted_gpm_bi_pixels;
    mc_c->gpm_weighted = &put_weighted_gpm_bi_pixels;
}

void
BD_DECL(rcn_init_ciip_functions)(struct RCNFunctions *const rcn_funcs)
{
    struct CIIPFunctions *const ciip = &rcn_funcs->ciip;
    ciip->weighted = &put_weighted_ciip_pixels;
}
