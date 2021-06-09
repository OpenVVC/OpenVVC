#include <stddef.h>
#include <stdint.h>
#include <string.h>

#include "ovutils.h"

#include "rcn_mc.h"
#include "rcn_structures.h"

#define ov_clip_pixel(a) ov_clip_uintp2(a, BIT_DEPTH)
#define MAX_PB_SIZE 128
#define BIT_DEPTH 10

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

/* FIXME reduce duplicated table */
static const int8_t ov_mc_filters[15][8] =
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
    {   0, 1,  -2,  4, 63,  -3,  1,  0 }
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
            dst[x] = src[x] << (14 - BIT_DEPTH);
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
    int shift = 14 - BIT_DEPTH + 1;
    int offset = 1 << (shift - 1);

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; ++x) {
            dst[x] = ov_clip_pixel(
                                   ((src0[x] << (14 - BIT_DEPTH)) + src1[x] + offset) >>
                                   shift);
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
    int shift = 14 - BIT_DEPTH;
    int offset = 1 << (shift - 1);

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = ov_clip_pixel(
                                   ((MCP_FILTER_L(src, 1, filter) >> (BIT_DEPTH - 8)) +
                                    offset) >> shift);
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
    int shift = 14 - BIT_DEPTH;
    int offset = 1 << (shift - 1);

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] =
                ov_clip_pixel(((MCP_FILTER_L(src, srcstride, filter) >>
                                (BIT_DEPTH - 8)) + offset) >> shift);
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
    int shift = 14 - BIT_DEPTH;
    int offset = 1 << (shift - 1);

    src -= QPEL_EXTRA_BEFORE * srcstride;
    filter = width == 4 && height == 4 ? ov_mc_filters_4[mx - 1] : ov_mc_filters[mx - 1];
    for (y = 0; y < height + QPEL_EXTRA; y++) {
        for (x = 0; x < width; x++) {
            tmp[x] = MCP_FILTER_L(src, 1, filter) >> (BIT_DEPTH - 8);
        }
        src += srcstride;
        tmp += MAX_PB_SIZE;
    }

    tmp = tmp_array + QPEL_EXTRA_BEFORE * MAX_PB_SIZE;
    filter = width == 4 && height == 4 ? ov_mc_filters_4[my - 1] : ov_mc_filters[my - 1];

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = ov_clip_pixel(
                                   ((MCP_FILTER_L(tmp, MAX_PB_SIZE, filter) >> 6) +
                                    offset) >> shift);
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

    const int8_t* filter = ov_mc_filters[mx - 1];
    filter = width == 4 && height == 4 ? ov_mc_filters_4[mx - 1] : ov_mc_filters[mx - 1];

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = MCP_FILTER_L(src, 1, filter) >> (BIT_DEPTH - 8);
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
                (BIT_DEPTH - 8);
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
            tmp[x] = MCP_FILTER_L(src, 1, filter) >> (BIT_DEPTH - 8);
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

    int shift = 14 + 1 - BIT_DEPTH;
    int offset = 1 << (shift - 1);

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = ov_clip_pixel(
                                   ((MCP_FILTER_L(src0, 1, filter) >> (BIT_DEPTH - 8)) +
                                    src1[x] + offset) >> shift);
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

    int shift = 14 + 1 - BIT_DEPTH;
    int offset = 1 << (shift - 1);

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = ov_clip_pixel(
                                   ((MCP_FILTER_L(src0, srcstride, filter) >>
                                     (BIT_DEPTH - 8)) +
                                    src1[x] + offset) >> shift);
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

    int shift = 14 + 1 - BIT_DEPTH;
    int offset = 1 << (shift - 1);

    src0 -= QPEL_EXTRA_BEFORE * srcstride;

    for (y = 0; y < height + QPEL_EXTRA; y++) {
        for (x = 0; x < width; x++) {
            tmp[x] =
                MCP_FILTER_L(src0, 1, filter) >> (BIT_DEPTH - 8);
        }
        src0 += srcstride;
        tmp += MAX_PB_SIZE;
    }

    tmp = tmp_array + QPEL_EXTRA_BEFORE * MAX_PB_SIZE;
    filter = ov_mc_filters[my - 1];
    filter = width == 4 && height == 4 ? ov_mc_filters_4[my - 1] : ov_mc_filters[my - 1];

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = ov_clip_pixel(
                                   ((MCP_FILTER_L(tmp, MAX_PB_SIZE, filter) >> 6) +
                                    src1[x] + offset) >> shift);
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
    int shift = 14 - BIT_DEPTH;
    int offset = 1 << (shift - 1);

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = ov_clip_pixel(
                                   ((MCP_FILTER_C(src, 1, filter) >> (BIT_DEPTH - 8)) +
                                    offset) >> shift);
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
    int shift = 14 - BIT_DEPTH;
    int offset = 1 << (shift - 1);

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] =
                ov_clip_pixel(((MCP_FILTER_C(src, srcstride, filter) >>
                                (BIT_DEPTH - 8)) + offset) >> shift);
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
    int shift = 14 - BIT_DEPTH;
    int offset = 1 << (shift - 1);

    src -= EPEL_EXTRA_BEFORE * srcstride;

    for (y = 0; y < height + EPEL_EXTRA; y++) {
        for (x = 0; x < width; x++) {
            tmp[x] = MCP_FILTER_C(src, 1, filter) >> (BIT_DEPTH - 8);
        }
        src += srcstride;
        tmp += MAX_PB_SIZE;
    }

    tmp = tmp_array + EPEL_EXTRA_BEFORE * MAX_PB_SIZE;
    filter = ov_mcp_filters_c[my - 1];

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = ov_clip_pixel(
                                   ((MCP_FILTER_C(tmp, MAX_PB_SIZE, filter) >> 6) +
                                    offset) >> shift);
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
            dst[x] = MCP_FILTER_C(src, 1, filter) >> (BIT_DEPTH - 8);
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
            dst[x] = MCP_FILTER_C(src, srcstride, filter) >>
                (BIT_DEPTH - 8);
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
            tmp[x] = MCP_FILTER_C(src, 1, filter) >> (BIT_DEPTH - 8);
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

    int shift = 14 + 1 - BIT_DEPTH;
    int offset = 1 << (shift - 1);

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = ov_clip_pixel(
                                   ((MCP_FILTER_C(src0, 1, filter) >> (BIT_DEPTH - 8)) +
                                    src1[x] + offset) >> shift);
        }
        src0 += srcstride;
        src1 += MAX_PB_SIZE;
        dst += dststride;
    }
}


static void
put_weighted_epel_bi_h(uint16_t* _dst, ptrdiff_t _dststride, const uint16_t* _src0,
                  ptrdiff_t _srcstride, const int16_t* _src1, int height,
                  intptr_t mx, intptr_t my, int width, int wt1, int wt0)
{
    int x, y;
    const uint16_t* src0 = _src0;
    const int16_t* src1 = (int16_t*)_src1;

    uint16_t* dst = (uint16_t*)_dst;

    ptrdiff_t srcstride = _srcstride;
    ptrdiff_t dststride = _dststride;

    const int8_t* filter = ov_mcp_filters_c[mx - 1];

    int shift = 14 + 3 - BIT_DEPTH;
    int offset = 1 << (shift - 1);

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = ov_clip_pixel(
                                   ((MCP_FILTER_C(src0, 1, filter) >> (BIT_DEPTH - 8)) * wt0 +
                                    src1[x] * wt1 + offset) >> shift);
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

    int shift = 14 + 1 - BIT_DEPTH;
    int offset = 1 << (shift - 1);

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = ov_clip_pixel(
                                   ((MCP_FILTER_C(src0, srcstride, filter) >>
                                     (BIT_DEPTH - 8)) +
                                    src1[x] + offset) >> shift);
        }
        src0 += srcstride;
        src1 += MAX_PB_SIZE;
        dst += dststride;
    }
}

static void
put_weighted_epel_bi_v(uint16_t* _dst, ptrdiff_t _dststride, const uint16_t* _src0,
                  ptrdiff_t _srcstride, const int16_t* _src1, int height,
                  intptr_t mx, intptr_t my, int width, int wt1, int wt0)
{
    int x, y;
    const uint16_t* src0 = _src0;
    const int16_t* src1 = _src1;

    uint16_t* dst = (uint16_t*)_dst;

    ptrdiff_t srcstride = _srcstride;
    ptrdiff_t dststride = _dststride;

    const int8_t* filter = ov_mcp_filters_c[my - 1];

    int shift = 14 + 3 - BIT_DEPTH;
    int offset = 1 << (shift - 1);

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = ov_clip_pixel(
                                   ((MCP_FILTER_C(src0, srcstride, filter) >>
                                     (BIT_DEPTH - 8)) * wt0 +
                                    src1[x] * wt1 + offset) >> shift);
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

    int shift = 14 + 1 - BIT_DEPTH;
    int offset = 1 << (shift - 1);

    src0 -= EPEL_EXTRA_BEFORE * srcstride;

    for (y = 0; y < height + EPEL_EXTRA; y++) {
        for (x = 0; x < width; x++) {
            tmp[x] =
                MCP_FILTER_C(src0, 1, filter) >> (BIT_DEPTH - 8);
        }
        src0 += srcstride;
        tmp += MAX_PB_SIZE;
    }

    tmp = tmp_array + EPEL_EXTRA_BEFORE * MAX_PB_SIZE;
    filter = ov_mcp_filters_c[my - 1];

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = ov_clip_pixel(
                                   ((MCP_FILTER_C(tmp, MAX_PB_SIZE, filter) >> 6) +
                                    src1[x] + offset) >> shift);
        }
        tmp += MAX_PB_SIZE;
        src1 += MAX_PB_SIZE;
        dst += dststride;
    }
}

static void
put_weighted_epel_bi_hv(uint16_t* _dst, ptrdiff_t _dststride, const uint16_t* _src0,
                   ptrdiff_t _srcstride, const int16_t* _src1, int height,
                   intptr_t mx, intptr_t my, int width, int wt1, int wt0)
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

    int shift = 14 + 3 - BIT_DEPTH;
    int offset = 1 << (shift - 1);

    src0 -= EPEL_EXTRA_BEFORE * srcstride;

    for (y = 0; y < height + EPEL_EXTRA; y++) {
        for (x = 0; x < width; x++) {
            tmp[x] =
                MCP_FILTER_C(src0, 1, filter) >> (BIT_DEPTH - 8);
        }
        src0 += srcstride;
        tmp += MAX_PB_SIZE;
    }

    tmp = tmp_array + EPEL_EXTRA_BEFORE * MAX_PB_SIZE;
    filter = ov_mcp_filters_c[my - 1];

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = ov_clip_pixel(
                                   ((MCP_FILTER_C(tmp, MAX_PB_SIZE, filter) >> 6) * wt0 +
                                    src1[x] * wt1 + offset) >> shift);
        }
        tmp += MAX_PB_SIZE;
        src1 += MAX_PB_SIZE;
        dst += dststride;
    }
}

void
put_weighted_pel_bi_pixels(uint16_t* _dst, ptrdiff_t _dststride, const uint16_t* _src0,
                  ptrdiff_t _srcstride, const int16_t* _src1, int height,
                  intptr_t mx, intptr_t my, int width, int wt1, int wt0)
{   
    int x, y;
    const uint16_t* src0 = _src0;
    const int16_t* src1 = _src1;
    ptrdiff_t srcstride = _srcstride;
    uint16_t* dst = (uint16_t*)_dst;
    ptrdiff_t dststride = _dststride;
    int shift = 14 - BIT_DEPTH + 3;
    int offset = (1 << (shift - 1)) ;
    for (y = 0; y < height; y++) {
        for (x = 0; x < width; ++x) {
            dst[x] = ov_clip_pixel( ( (src1[x]  * wt1
                                        + ((src0[x] * wt0) << (14 - BIT_DEPTH)) + offset )  >> shift ) );
        }
        src1 += MAX_PB_SIZE;
        src0 += srcstride;
        dst += dststride;
    }
}

static void
put_weighted_qpel_bi_h(uint16_t* _dst, ptrdiff_t _dststride, const uint16_t* _src0,
                  ptrdiff_t _srcstride, const int16_t* _src1, int height,
                  intptr_t mx, intptr_t my, int width, int wt1, int wt0)
{
    int x, y;
    const uint16_t* src0 = _src0;
    const int16_t* src1 = _src1;

    uint16_t* dst = (uint16_t*)_dst;

    ptrdiff_t srcstride = _srcstride;
    ptrdiff_t dststride = _dststride;

    const int8_t* filter = ov_mc_filters[mx - 1];

    int shift = 14 + 3 - BIT_DEPTH;
    int offset = 1 << (shift - 1);

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = ov_clip_pixel(
                                   ((MCP_FILTER_L(src0, 1, filter) >> (BIT_DEPTH - 8)) * wt0 +
                                    src1[x] * wt1 + offset) >> shift);
        }
        src0 += srcstride;
        src1 += MAX_PB_SIZE;
        dst += dststride;
    }
}


static void
put_weighted_qpel_bi_v(uint16_t* _dst, ptrdiff_t _dststride, const uint16_t* _src0,
                  ptrdiff_t _srcstride, const int16_t* _src1, int height,
                  intptr_t mx, intptr_t my, int width, int wt1, int wt0)
{
    int x, y;
    const uint16_t* src0 = _src0;
    const int16_t* src1 = _src1;

    uint16_t* dst = (uint16_t*)_dst;

    ptrdiff_t srcstride = _srcstride;
    ptrdiff_t dststride = _dststride;

    const int8_t* filter = ov_mc_filters[my - 1];

    int shift = 14 + 3 - BIT_DEPTH;
    int offset = 1 << (shift - 1);

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = ov_clip_pixel(
                                   ((MCP_FILTER_L(src0, srcstride, filter) >>
                                     (BIT_DEPTH - 8)) * wt0 +
                                    src1[x] * wt1 + offset) >> shift);
        }
        src0 += srcstride;
        src1 += MAX_PB_SIZE;
        dst += dststride;
    }
}


static void
put_weighted_qpel_bi_hv(uint16_t* _dst, ptrdiff_t _dststride, const uint16_t* _src0,
                   ptrdiff_t _srcstride, const int16_t* _src1, int height,
                   intptr_t mx, intptr_t my, int width, int wt1, int wt0)
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

    //TODObcw: change 3 with a parameter log2_sum_weights
    int shift = 14 + 3 - BIT_DEPTH;
    int offset = 1 << (shift - 1);

    src0 -= QPEL_EXTRA_BEFORE * srcstride;

    for (y = 0; y < height + QPEL_EXTRA; y++) {
        for (x = 0; x < width; x++) {
            tmp[x] =
                MCP_FILTER_L(src0, 1, filter) >> (BIT_DEPTH - 8);
        }
        src0 += srcstride;
        tmp += MAX_PB_SIZE;
    }

    tmp = tmp_array + QPEL_EXTRA_BEFORE * MAX_PB_SIZE;
    filter = ov_mc_filters[my - 1];

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = ov_clip_pixel(
                                   ((MCP_FILTER_L(tmp, MAX_PB_SIZE, filter) >> 6) * wt0 +
                                    src1[x] * wt1 + offset) >> shift);
        }
        tmp += MAX_PB_SIZE;
        src1 += MAX_PB_SIZE;
        dst += dststride;
    }
}

void
put_weighted_ciip_pixels(uint16_t* dst, int dststride,
                      const uint16_t* src_intra, const uint16_t* src_inter, int srcstride,
                      int width, int height, int wt)
{   
    int x, y;
    int shift  = 2;
    int offset = 2;
    for (y = 0; y < height; y++) {
        for (x = 0; x < width; ++x) {
            dst[x] = ov_clip_pixel( ( ( src_intra[x] * wt + src_inter[x] * (4 - wt)) + offset ) >> shift );
        }
        src_intra += srcstride;
        src_inter += srcstride;
        dst += dststride;
    }
}


void
put_gpm_pel_bi_pixels(uint16_t* _dst, int _dststride, const int16_t* _src0,
                  int _srcstride, const int16_t* _src1, int height,
                  intptr_t mx, intptr_t my, int width, int step_x, int step_y, int16_t* weight)
{   
    int x, y;
    const int16_t* src0 = _src0;
    const int16_t* src1 = _src1;
    ptrdiff_t srcstride = _srcstride;
    uint16_t* dst = (uint16_t*)_dst;
    ptrdiff_t dststride = _dststride;
    int shift = 14 - BIT_DEPTH + 3;
    int offset = (1 << (shift - 1)) ;
    for (y = 0; y < height; y++) {
        for (x = 0; x < width; ++x) {
            int w1 = weight[0];
            int w0 = 8 - w1;
            dst[x] = ov_clip_pixel( ( (src1[x]  * w1
                                        + src0[x] * w0 + offset )  >> shift ) );
            weight += step_x;
        }
        src1 += srcstride;
        src0 += srcstride;
        dst += dststride;
        weight += step_y;
    }
}

void rcn_init_mc_functions(struct RCNFunctions *const rcn_funcs)
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
}
