#include <libovvcutils/ovvcutils.h>
#include <stddef.h>
#include <stdint.h>

#include "rcn_mc.h"

#define ov_clip_pixel(a) ov_clip_uintp2(a, BIT_DEPTH)
#define MAX_PB_SIZE 128
#define BIT_DEPTH 10

#define EPEL_EXTRA_BEFORE 1
#define EPEL_EXTRA_AFTER 2
#define EPEL_EXTRA EPEL_EXTRA_BEFORE + EPEL_EXTRA_AFTER

#define QPEL_EXTRA_BEFORE 3
#define QPEL_EXTRA_AFTER 4
#define QPEL_EXTRA QPEL_EXTRA_BEFORE + QPEL_EXTRA_AFTER

#define POS_IN_CTB(x, y) (VVC_CTB_OFFSET + (x) + (y)*VVC_CTB_STRIDE)
#define POS_IN_CTB_C(x, y)                                                     \
        (VVC_CTB_OFFSET_CHROMA + (x) + (y)*VVC_CTB_STRIDE_CHROMA)

static const int8_t ff_vvc_epel_filters[7][4] = {
        { -2, 58, 10, -2 }, { -4, 54, 16, -2 }, { -6, 46, 28, -4 },
        { -4, 36, 36, -4 }, { -4, 28, 46, -6 }, { -2, 16, 54, -4 },
        { -2, 10, 58, -2 },
};

static const int8_t vvc_qpel_filters[3][16] = {
        { -1, 4, -10, 58, 17, -5, 1, 0, -1, 4, -10, 58, 17, -5, 1, 0 },
        { -1, 4, -11, 40, 40, -11, 4, -1, -1, 4, -11, 40, 40, -11, 4, -1 },
        { 0, 1, -5, 17, 58, -10, 4, -1, 0, 1, -5, 17, 58, -10, 4, -1 }
};

void
put_vvc_pel_uni_pixels(uint16_t* _dst, ptrdiff_t _dststride,
                       const uint16_t* _src, ptrdiff_t _srcstride, int height,
                       intptr_t mx, intptr_t my, int width)
{
        int y;
        uint16_t* src = (uint16_t*)_src;
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
void
put_vvc_pel_pixels(int16_t* _dst, const uint16_t* _src, ptrdiff_t _srcstride,
                   int height, intptr_t mx, intptr_t my, int width)
{
        int x, y;
        const uint16_t* src = (uint16_t*)_src;

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

void
put_vvc_pel_bi_pixels(uint16_t* _dst, ptrdiff_t _dststride,
                      const uint16_t* _src0, ptrdiff_t _srcstride,
                      const int16_t* _src1, int height, intptr_t mx,
                      intptr_t my, int width)
{
        int x, y;
        uint16_t* src0 = (uint16_t*)_src0;
        int16_t* src1 = (int16_t*)_src1;
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

#define QPEL_FILTER(src, stride, filter)                                       \
        (filter[0] * src[x - 3 * stride] + filter[1] * src[x - 2 * stride] +   \
         filter[2] * src[x - stride] + filter[3] * src[x] +                    \
         filter[4] * src[x + stride] + filter[5] * src[x + 2 * stride] +       \
         filter[6] * src[x + 3 * stride] + filter[7] * src[x + 4 * stride])

void
put_vvc_qpel_uni_h(uint16_t* _dst, ptrdiff_t _dststride, const uint16_t* _src,
                   ptrdiff_t _srcstride, int height, intptr_t mx, intptr_t my,
                   int width)
{
        int x, y;
        uint16_t* src = (uint16_t*)_src;
        ptrdiff_t srcstride = _srcstride;
        uint16_t* dst = (uint16_t*)_dst;
        ptrdiff_t dststride = _dststride;
        const int8_t* filter = vvc_qpel_filters[mx - 1];
        int shift = 14 - BIT_DEPTH;
        int offset = 1 << (shift - 1);

        for (y = 0; y < height; y++) {
                for (x = 0; x < width; x++) {
                        dst[x] = ov_clip_pixel(
                          ((QPEL_FILTER(src, 1, filter) >> (BIT_DEPTH - 8)) +
                           offset) >>
                          shift);
                }
                src += srcstride;
                dst += dststride;
        }
}

void
put_vvc_qpel_uni_v(uint16_t* _dst, ptrdiff_t _dststride, const uint16_t* _src,
                   ptrdiff_t _srcstride, int height, intptr_t mx, intptr_t my,
                   int width)
{
        int x, y;
        uint16_t* src = (uint16_t*)_src;
        ptrdiff_t srcstride = _srcstride;
        uint16_t* dst = (uint16_t*)_dst;
        ptrdiff_t dststride = _dststride;
        const int8_t* filter = vvc_qpel_filters[my - 1];
        int shift = 14 - BIT_DEPTH;
        int offset = 1 << (shift - 1);

        for (y = 0; y < height; y++) {
                for (x = 0; x < width; x++) {
                        dst[x] =
                          ov_clip_pixel(((QPEL_FILTER(src, srcstride, filter) >>
                                          (BIT_DEPTH - 8)) +
                                         offset) >>
                                        shift);
                }
                src += srcstride;
                dst += dststride;
        }
}

void
put_vvc_qpel_uni_hv(uint16_t* _dst, ptrdiff_t _dststride, const uint16_t* _src,
                    ptrdiff_t _srcstride, int height, intptr_t mx, intptr_t my,
                    int width)
{
        int x, y;
        const int8_t* filter;
        uint16_t* src = (uint16_t*)_src;
        ptrdiff_t srcstride = _srcstride;
        uint16_t* dst = (uint16_t*)_dst;
        ptrdiff_t dststride = _dststride;
        int16_t tmp_array[(MAX_PB_SIZE + QPEL_EXTRA) * MAX_PB_SIZE];
        int16_t* tmp = tmp_array;
        int shift = 14 - BIT_DEPTH;
        int offset = 1 << (shift - 1);

        src -= QPEL_EXTRA_BEFORE * srcstride;
        filter = vvc_qpel_filters[mx - 1];
        for (y = 0; y < height + QPEL_EXTRA; y++) {
                for (x = 0; x < width; x++) {
                        tmp[x] = QPEL_FILTER(src, 1, filter) >> (BIT_DEPTH - 8);
                }
                src += srcstride;
                tmp += MAX_PB_SIZE;
        }

        tmp = tmp_array + QPEL_EXTRA_BEFORE * MAX_PB_SIZE;
        filter = vvc_qpel_filters[my - 1];

        for (y = 0; y < height; y++) {
                for (x = 0; x < width; x++) {
                        dst[x] = ov_clip_pixel(
                          ((QPEL_FILTER(tmp, MAX_PB_SIZE, filter) >> 6) +
                           offset) >>
                          shift);
                }
                tmp += MAX_PB_SIZE;
                dst += dststride;
        }
}

void
put_vvc_qpel_h(int16_t* _dst, const uint16_t* _src, ptrdiff_t _srcstride,
               int height, intptr_t mx, intptr_t my, int width)
{
        int x, y;
        const uint16_t* src = (uint16_t*)_src;

        int16_t* dst = (int16_t*)_dst;

        ptrdiff_t srcstride = _srcstride;

        const int8_t* filter = vvc_qpel_filters[mx - 1];

        for (y = 0; y < height; y++) {
                for (x = 0; x < width; x++) {
                        dst[x] = QPEL_FILTER(src, 1, filter) >> (BIT_DEPTH - 8);
                }
                src += srcstride;
                dst += MAX_PB_SIZE;
        }
}

void
put_vvc_qpel_v(int16_t* _dst, const uint16_t* _src, ptrdiff_t _srcstride,
               int height, intptr_t mx, intptr_t my, int width)
{
        int x, y;
        const uint16_t* src = (uint16_t*)_src;

        int16_t* dst = (int16_t*)_dst;

        ptrdiff_t srcstride = _srcstride;

        const int8_t* filter = vvc_qpel_filters[my - 1];

        for (y = 0; y < height; y++) {
                for (x = 0; x < width; x++) {
                        dst[x] = QPEL_FILTER(src, srcstride, filter) >>
                                 (BIT_DEPTH - 8);
                }
                src += srcstride;
                dst += MAX_PB_SIZE;
        }
}

void
put_vvc_qpel_hv(int16_t* _dst, const uint16_t* _src, ptrdiff_t _srcstride,
                int height, intptr_t mx, intptr_t my, int width)
{
        int x, y;

        const uint16_t* src = (uint16_t*)_src;

        int16_t* dst = (int16_t*)_dst;

        ptrdiff_t srcstride = _srcstride;

        int16_t tmp_array[(MAX_PB_SIZE + QPEL_EXTRA) * MAX_PB_SIZE];
        int16_t* tmp = tmp_array;

        const int8_t* filter = vvc_qpel_filters[mx - 1];

        src -= QPEL_EXTRA_BEFORE * srcstride;

        for (y = 0; y < height + QPEL_EXTRA; y++) {
                for (x = 0; x < width; x++) {
                        tmp[x] = QPEL_FILTER(src, 1, filter) >> (BIT_DEPTH - 8);
                }
                src += srcstride;
                tmp += MAX_PB_SIZE;
        }

        tmp = tmp_array + QPEL_EXTRA_BEFORE * MAX_PB_SIZE;
        filter = vvc_qpel_filters[my - 1];

        for (y = 0; y < height; y++) {
                for (x = 0; x < width; x++) {
                        dst[x] = QPEL_FILTER(tmp, MAX_PB_SIZE, filter) >> 6;
                }
                tmp += MAX_PB_SIZE;
                dst += MAX_PB_SIZE;
        }
}

void
put_vvc_qpel_bi_h(uint16_t* _dst, ptrdiff_t _dststride, const uint16_t* _src0,
                  ptrdiff_t _srcstride, const int16_t* _src1, int height,
                  intptr_t mx, intptr_t my, int width)
{
        int x, y;
        const uint16_t* src0 = (uint16_t*)_src0;
        const int16_t* src1 = (int16_t*)_src1;

        uint16_t* dst = (uint16_t*)_dst;

        ptrdiff_t srcstride = _srcstride;
        ptrdiff_t dststride = _dststride;

        const int8_t* filter = vvc_qpel_filters[mx - 1];

        int shift = 14 + 1 - BIT_DEPTH;
        int offset = 1 << (shift - 1);

        for (y = 0; y < height; y++) {
                for (x = 0; x < width; x++) {
                        dst[x] = ov_clip_pixel(
                          ((QPEL_FILTER(src0, 1, filter) >> (BIT_DEPTH - 8)) +
                           src1[x] + offset) >>
                          shift);
                }
                src0 += srcstride;
                src1 += MAX_PB_SIZE;
                dst += dststride;
        }
}

void
put_vvc_qpel_bi_v(uint16_t* _dst, ptrdiff_t _dststride, const uint16_t* _src0,
                  ptrdiff_t _srcstride, const int16_t* _src1, int height,
                  intptr_t mx, intptr_t my, int width)
{
        int x, y;
        const uint16_t* src0 = (uint16_t*)_src0;
        const int16_t* src1 = (int16_t*)_src1;

        uint16_t* dst = (uint16_t*)_dst;

        ptrdiff_t srcstride = _srcstride;
        ptrdiff_t dststride = _dststride;

        const int8_t* filter = vvc_qpel_filters[my - 1];

        int shift = 14 + 1 - BIT_DEPTH;
        int offset = 1 << (shift - 1);

        for (y = 0; y < height; y++) {
                for (x = 0; x < width; x++) {
                        dst[x] = ov_clip_pixel(
                          ((QPEL_FILTER(src0, srcstride, filter) >>
                            (BIT_DEPTH - 8)) +
                           src1[x] + offset) >>
                          shift);
                }
                src0 += srcstride;
                src1 += MAX_PB_SIZE;
                dst += dststride;
        }
}

void
put_vvc_qpel_bi_hv(uint16_t* _dst, ptrdiff_t _dststride, const uint16_t* _src0,
                   ptrdiff_t _srcstride, const int16_t* _src1, int height,
                   intptr_t mx, intptr_t my, int width)
{
        int x, y;

        uint16_t* src0 = (uint16_t*)_src0;
        int16_t* src1 = (int16_t*)_src1;
        uint16_t* dst = (uint16_t*)_dst;

        ptrdiff_t srcstride = _srcstride;
        ptrdiff_t dststride = _dststride;

        int16_t tmp_array[(MAX_PB_SIZE + QPEL_EXTRA) * MAX_PB_SIZE];
        int16_t* tmp = tmp_array;

        const int8_t* filter = vvc_qpel_filters[mx - 1];

        int shift = 14 + 1 - BIT_DEPTH;
        int offset = 1 << (shift - 1);

        src0 -= QPEL_EXTRA_BEFORE * srcstride;

        for (y = 0; y < height + QPEL_EXTRA; y++) {
                for (x = 0; x < width; x++) {
                        tmp[x] =
                          QPEL_FILTER(src0, 1, filter) >> (BIT_DEPTH - 8);
                }
                src0 += srcstride;
                tmp += MAX_PB_SIZE;
        }

        tmp = tmp_array + QPEL_EXTRA_BEFORE * MAX_PB_SIZE;
        filter = vvc_qpel_filters[my - 1];

        for (y = 0; y < height; y++) {
                for (x = 0; x < width; x++) {
                        dst[x] = ov_clip_pixel(
                          ((QPEL_FILTER(tmp, MAX_PB_SIZE, filter) >> 6) +
                           src1[x] + offset) >>
                          shift);
                }
                tmp += MAX_PB_SIZE;
                src1 += MAX_PB_SIZE;
                dst += dststride;
        }
}

#define EPEL_FILTER(src, stride, filter)                                       \
        (filter[0] * src[x - stride] + filter[1] * src[x] +                    \
         filter[2] * src[x + stride] + filter[3] * src[x + 2 * stride])

void
put_vvc_epel_uni_h(uint16_t* _dst, ptrdiff_t _dststride, const uint16_t* _src,
                   ptrdiff_t _srcstride, int height, intptr_t mx, intptr_t my,
                   int width)
{
        int x, y;
        uint16_t* src = (uint16_t*)_src;
        ptrdiff_t srcstride = _srcstride;
        uint16_t* dst = (uint16_t*)_dst;
        ptrdiff_t dststride = _dststride;
        const int8_t* filter = ff_vvc_epel_filters[mx - 1];
        int shift = 14 - BIT_DEPTH;
        int offset = 1 << (shift - 1);

        for (y = 0; y < height; y++) {
                for (x = 0; x < width; x++) {
                        dst[x] = ov_clip_pixel(
                          ((EPEL_FILTER(src, 1, filter) >> (BIT_DEPTH - 8)) +
                           offset) >>
                          shift);
                }
                src += srcstride;
                dst += dststride;
        }
}

void
put_vvc_epel_uni_v(uint16_t* _dst, ptrdiff_t _dststride, const uint16_t* _src,
                   ptrdiff_t _srcstride, int height, intptr_t mx, intptr_t my,
                   int width)
{
        int x, y;
        uint16_t* src = (uint16_t*)_src;
        ptrdiff_t srcstride = _srcstride;
        uint16_t* dst = (uint16_t*)_dst;
        ptrdiff_t dststride = _dststride;
        const int8_t* filter = ff_vvc_epel_filters[my - 1];
        int shift = 14 - BIT_DEPTH;
        int offset = 1 << (shift - 1);

        for (y = 0; y < height; y++) {
                for (x = 0; x < width; x++) {
                        dst[x] =
                          ov_clip_pixel(((EPEL_FILTER(src, srcstride, filter) >>
                                          (BIT_DEPTH - 8)) +
                                         offset) >>
                                        shift);
                }
                src += srcstride;
                dst += dststride;
        }
}

void
put_vvc_epel_uni_hv(uint16_t* _dst, ptrdiff_t _dststride, const uint16_t* _src,
                    ptrdiff_t _srcstride, int height, intptr_t mx, intptr_t my,
                    int width)
{
        int x, y;
        uint16_t* src = (uint16_t*)_src;
        ptrdiff_t srcstride = _srcstride;
        uint16_t* dst = (uint16_t*)_dst;
        ptrdiff_t dststride = _dststride;
        const int8_t* filter = ff_vvc_epel_filters[mx - 1];
        int16_t tmp_array[(MAX_PB_SIZE + EPEL_EXTRA) * MAX_PB_SIZE];
        int16_t* tmp = tmp_array;
        int shift = 14 - BIT_DEPTH;
        int offset = 1 << (shift - 1);

        src -= EPEL_EXTRA_BEFORE * srcstride;

        for (y = 0; y < height + EPEL_EXTRA; y++) {
                for (x = 0; x < width; x++) {
                        tmp[x] = EPEL_FILTER(src, 1, filter) >> (BIT_DEPTH - 8);
                }
                src += srcstride;
                tmp += MAX_PB_SIZE;
        }

        tmp = tmp_array + EPEL_EXTRA_BEFORE * MAX_PB_SIZE;
        filter = ff_vvc_epel_filters[my - 1];

        for (y = 0; y < height; y++) {
                for (x = 0; x < width; x++) {
                        dst[x] = ov_clip_pixel(
                          ((EPEL_FILTER(tmp, MAX_PB_SIZE, filter) >> 6) +
                           offset) >>
                          shift);
                }
                tmp += MAX_PB_SIZE;
                dst += dststride;
        }
}

void
put_vvc_epel_h(int16_t* _dst, const uint16_t* _src, ptrdiff_t _srcstride,
               int height, intptr_t mx, intptr_t my, int width)
{
        int x, y;

        const uint16_t* src = (uint16_t*)_src;

        int16_t* dst = (int16_t*)_dst;

        ptrdiff_t srcstride = _srcstride;

        const int8_t* filter = ff_vvc_epel_filters[mx - 1];

        for (y = 0; y < height; y++) {
                for (x = 0; x < width; x++) {
                        dst[x] = EPEL_FILTER(src, 1, filter) >> (BIT_DEPTH - 8);
                }
                src += srcstride;
                dst += MAX_PB_SIZE;
        }
}

void
put_vvc_epel_v(int16_t* _dst, const uint16_t* _src, ptrdiff_t _srcstride,
               int height, intptr_t mx, intptr_t my, int width)
{
        int x, y;

        const uint16_t* src = (uint16_t*)_src;

        int16_t* dst = (int16_t*)_dst;

        ptrdiff_t srcstride = _srcstride;

        const int8_t* filter = ff_vvc_epel_filters[my - 1];

        for (y = 0; y < height; y++) {
                for (x = 0; x < width; x++) {
                        dst[x] = EPEL_FILTER(src, srcstride, filter) >>
                                 (BIT_DEPTH - 8);
                }
                src += srcstride;
                dst += MAX_PB_SIZE;
        }
}

void
put_vvc_epel_hv(int16_t* _dst, const uint16_t* _src, ptrdiff_t _srcstride,
                int height, intptr_t mx, intptr_t my, int width)
{
        int x, y;
        const uint16_t* src = (uint16_t*)_src;

        int16_t* dst = (int16_t*)_dst;

        ptrdiff_t srcstride = _srcstride;

        const int8_t* filter = ff_vvc_epel_filters[mx - 1];

        int16_t tmp_array[(MAX_PB_SIZE + EPEL_EXTRA) * MAX_PB_SIZE];
        int16_t* tmp = tmp_array;

        src -= EPEL_EXTRA_BEFORE * srcstride;

        for (y = 0; y < height + EPEL_EXTRA; y++) {
                for (x = 0; x < width; x++) {
                        tmp[x] = EPEL_FILTER(src, 1, filter) >> (BIT_DEPTH - 8);
                }
                src += srcstride;
                tmp += MAX_PB_SIZE;
        }

        tmp = tmp_array + EPEL_EXTRA_BEFORE * MAX_PB_SIZE;
        filter = ff_vvc_epel_filters[my - 1];

        for (y = 0; y < height; y++) {
                for (x = 0; x < width; x++) {
                        dst[x] = EPEL_FILTER(tmp, MAX_PB_SIZE, filter) >> 6;
                }
                tmp += MAX_PB_SIZE;
                dst += MAX_PB_SIZE;
        }
}

void
put_vvc_epel_bi_h(uint16_t* _dst, ptrdiff_t _dststride, const uint16_t* _src0,
                  ptrdiff_t _srcstride, const int16_t* _src1, int height,
                  intptr_t mx, intptr_t my, int width)
{
        int x, y;
        const uint16_t* src0 = (uint16_t*)_src0;
        const int16_t* src1 = (int16_t*)_src1;

        uint16_t* dst = (uint16_t*)_dst;

        ptrdiff_t srcstride = _srcstride;
        ptrdiff_t dststride = _dststride;

        const int8_t* filter = ff_vvc_epel_filters[mx - 1];

        int shift = 14 + 1 - BIT_DEPTH;
        int offset = 1 << (shift - 1);

        for (y = 0; y < height; y++) {
                for (x = 0; x < width; x++) {
                        dst[x] = ov_clip_pixel(
                          ((EPEL_FILTER(src0, 1, filter) >> (BIT_DEPTH - 8)) +
                           src1[x] + offset) >>
                          shift);
                }
                src0 += srcstride;
                src1 += MAX_PB_SIZE;
                dst += dststride;
        }
}

void
put_vvc_epel_bi_v(uint16_t* _dst, ptrdiff_t _dststride, const uint16_t* _src0,
                  ptrdiff_t _srcstride, const int16_t* _src1, int height,
                  intptr_t mx, intptr_t my, int width)
{
        int x, y;
        const uint16_t* src0 = (uint16_t*)_src0;
        const int16_t* src1 = (int16_t*)_src1;

        uint16_t* dst = (uint16_t*)_dst;

        ptrdiff_t srcstride = _srcstride;
        ptrdiff_t dststride = _dststride;

        const int8_t* filter = ff_vvc_epel_filters[my - 1];

        int shift = 14 + 1 - BIT_DEPTH;
        int offset = 1 << (shift - 1);

        for (y = 0; y < height; y++) {
                for (x = 0; x < width; x++) {
                        dst[x] = ov_clip_pixel(
                          ((EPEL_FILTER(src0, srcstride, filter) >>
                            (BIT_DEPTH - 8)) +
                           src1[x] + offset) >>
                          shift);
                }
                src0 += srcstride;
                src1 += MAX_PB_SIZE;
                dst += dststride;
        }
}

void
put_vvc_epel_bi_hv(uint16_t* _dst, ptrdiff_t _dststride, const uint16_t* _src0,
                   ptrdiff_t _srcstride, const int16_t* _src1, int height,
                   intptr_t mx, intptr_t my, int width)
{
        int x, y;
        const uint16_t* src0 = (uint16_t*)_src0;
        const int16_t* src1 = (int16_t*)_src1;

        uint16_t* dst = (uint16_t*)_dst;

        ptrdiff_t srcstride = _srcstride;
        ptrdiff_t dststride = _dststride;

        int16_t tmp_array[(MAX_PB_SIZE + EPEL_EXTRA) * MAX_PB_SIZE];
        int16_t* tmp = tmp_array;

        const int8_t* filter = ff_vvc_epel_filters[mx - 1];

        int shift = 14 + 1 - BIT_DEPTH;
        int offset = 1 << (shift - 1);

        src0 -= EPEL_EXTRA_BEFORE * srcstride;

        for (y = 0; y < height + EPEL_EXTRA; y++) {
                for (x = 0; x < width; x++) {
                        tmp[x] =
                          EPEL_FILTER(src0, 1, filter) >> (BIT_DEPTH - 8);
                }
                src0 += srcstride;
                tmp += MAX_PB_SIZE;
        }

        tmp = tmp_array + EPEL_EXTRA_BEFORE * MAX_PB_SIZE;
        filter = ff_vvc_epel_filters[my - 1];

        for (y = 0; y < height; y++) {
                for (x = 0; x < width; x++) {
                        dst[x] = ov_clip_pixel(
                          ((EPEL_FILTER(tmp, MAX_PB_SIZE, filter) >> 6) +
                           src1[x] + offset) >>
                          shift);
                }
                tmp += MAX_PB_SIZE;
                src1 += MAX_PB_SIZE;
                dst += dststride;
        }
}
