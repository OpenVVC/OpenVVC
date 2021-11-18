#ifndef RCN_PROF_BDOF_H
#define RCN_PROF_BDOF_H
#include <stdint.h>
#include <stddef.h>

#define PROF_BUFF_STRIDE 128

void tmp_prof_mrg(uint16_t* _dst, ptrdiff_t _dststride,
                  const uint16_t* _src0, ptrdiff_t _srcstride,
                  const int16_t* _src1, int height, intptr_t mx,
                  intptr_t my, int width);

void tmp_prof_mrg_w(uint16_t* _dst, ptrdiff_t _dststride,
                    const uint16_t* _src0, ptrdiff_t _srcstride,
                    const int16_t* _src1, int height, intptr_t mx,
                    intptr_t my, int width, int wt0, int wt1);

void extend_prof_buff(const uint16_t *const src, uint16_t *dst_prof, int16_t ref_stride,
                      uint8_t ext_x, uint8_t ext_y);

void extend_bdof_buff(const uint16_t *const src, uint16_t *dst_prof,
                      int16_t ref_stride, int16_t pb_w, int16_t pb_h,
                      uint8_t ext_x, uint8_t ext_y);
struct BDOFFunctions;
void rcn_bdof(struct BDOFFunctions *const bdof, int16_t *dst, int dst_stride,
              const int16_t *ref_bdof0, const int16_t *ref_bdof1, int ref_stride,
              const int16_t *grad_x0, const int16_t *grad_y0,
              const int16_t *grad_x1, const int16_t *grad_y1,
              int grad_stride, uint8_t pb_w, uint8_t pb_h);

#endif

