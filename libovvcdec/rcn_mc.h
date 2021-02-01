#ifndef RCN_MC_H
#define RCN_MC_H

#include <stddef.h>
#include <stdint.h>

void
put_vvc_pel_uni_pixels(uint16_t* _dst, ptrdiff_t _dststride,
                       const uint16_t* _src, ptrdiff_t _srcstride, int height,
                       intptr_t mx, intptr_t my, int width);

void
put_vvc_pel_pixels(int16_t* _dst, const uint16_t* _src, ptrdiff_t _srcstride,
                   int height, intptr_t mx, intptr_t my, int width);

void
put_vvc_pel_bi_pixels(uint16_t* _dst, ptrdiff_t _dststride,
                      const uint16_t* _src0, ptrdiff_t _srcstride,
                      const int16_t* _src1, int height, intptr_t mx,
                      intptr_t my, int width);

void
put_vvc_qpel_uni_h(uint16_t* _dst, ptrdiff_t _dststride, const uint16_t* _src,
                   ptrdiff_t _srcstride, int height, intptr_t mx, intptr_t my,
                   int width);

void
put_vvc_qpel_uni_v(uint16_t* _dst, ptrdiff_t _dststride, const uint16_t* _src,
                   ptrdiff_t _srcstride, int height, intptr_t mx, intptr_t my,
                   int width);

void
put_vvc_qpel_uni_hv(uint16_t* _dst, ptrdiff_t _dststride, const uint16_t* _src,
                    ptrdiff_t _srcstride, int height, intptr_t mx, intptr_t my,
                    int width);

void
put_vvc_qpel_h(int16_t* _dst, const uint16_t* _src, ptrdiff_t _srcstride,
               int height, intptr_t mx, intptr_t my, int width);

void
put_vvc_qpel_v(int16_t* _dst, const uint16_t* _src, ptrdiff_t _srcstride,
               int height, intptr_t mx, intptr_t my, int width);

void
put_vvc_qpel_hv(int16_t* _dst, const uint16_t* _src, ptrdiff_t _srcstride,
                int height, intptr_t mx, intptr_t my, int width);

void
put_vvc_qpel_bi_h(uint16_t* _dst, ptrdiff_t _dststride, const uint16_t* _src0,
                  ptrdiff_t _srcstride, const int16_t* _src1, int height,
                  intptr_t mx, intptr_t my, int width);

void
put_vvc_qpel_bi_v(uint16_t* _dst, ptrdiff_t _dststride, const uint16_t* _src0,
                  ptrdiff_t _srcstride, const int16_t* _src1, int height,
                  intptr_t mx, intptr_t my, int width);

void
put_vvc_qpel_bi_hv(uint16_t* _dst, ptrdiff_t _dststride, const uint16_t* _src0,
                   ptrdiff_t _srcstride, const int16_t* _src1, int height,
                   intptr_t mx, intptr_t my, int width);

void
put_vvc_epel_uni_h(uint16_t* _dst, ptrdiff_t _dststride, const uint16_t* _src,
                   ptrdiff_t _srcstride, int height, intptr_t mx, intptr_t my,
                   int width);

void
put_vvc_epel_uni_v(uint16_t* _dst, ptrdiff_t _dststride, const uint16_t* _src,
                   ptrdiff_t _srcstride, int height, intptr_t mx, intptr_t my,
                   int width);

void
put_vvc_epel_uni_hv(uint16_t* _dst, ptrdiff_t _dststride, const uint16_t* _src,
                    ptrdiff_t _srcstride, int height, intptr_t mx, intptr_t my,
                    int width);

void
put_vvc_epel_h(int16_t* _dst, const uint16_t* _src, ptrdiff_t _srcstride,
               int height, intptr_t mx, intptr_t my, int width);

void
put_vvc_epel_v(int16_t* _dst, const uint16_t* _src, ptrdiff_t _srcstride,
               int height, intptr_t mx, intptr_t my, int width);

void
put_vvc_epel_hv(int16_t* _dst, const uint16_t* _src, ptrdiff_t _srcstride,
                int height, intptr_t mx, intptr_t my, int width);

void
put_vvc_epel_bi_h(uint16_t* _dst, ptrdiff_t _dststride, const uint16_t* _src0,
                  ptrdiff_t _srcstride, const int16_t* _src1, int height,
                  intptr_t mx, intptr_t my, int width);

void
put_vvc_epel_bi_v(uint16_t* _dst, ptrdiff_t _dststride, const uint16_t* _src0,
                  ptrdiff_t _srcstride, const int16_t* _src1, int height,
                  intptr_t mx, intptr_t my, int width);

void
put_vvc_epel_bi_hv(uint16_t* _dst, ptrdiff_t _dststride, const uint16_t* _src0,
                   ptrdiff_t _srcstride, const int16_t* _src1, int height,
                   intptr_t mx, intptr_t my, int width);

#endif // RCN_MC_H
