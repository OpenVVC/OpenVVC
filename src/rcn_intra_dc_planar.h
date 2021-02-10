#ifndef RCN_INTRA_DC_PLANAR_H
#define RCN_INTRA_DC_PLANAR_H

#include <stdint.h>

void
vvc_intra_dc(const uint16_t* const src_above, const uint16_t* const src_left,
             uint16_t* const dst, ptrdiff_t dst_stride, int log2_pb_width,
             int log2_pb_height);

void
vvc_intra_planar(const uint16_t* const src_above,
                 const uint16_t* const src_left, uint16_t* const dst,
                 ptrdiff_t dst_stride, int log2_pb_width, int log2_pb_height);

void
vvc_intra_dc_pdpc(const uint16_t* const src_above,
                  const uint16_t* const src_left, uint16_t* const dst,
                  ptrdiff_t dst_stride, int log2_pb_w, int log2_pb_h);

void
vvc_intra_planar_pdpc(const uint16_t* const src_above,
                      const uint16_t* const src_left, uint16_t* const dst,
                      ptrdiff_t dst_stride, int log2_pb_w, int log2_pb_h);

#endif // RCN_INTRA_DC_PLANAR_H
