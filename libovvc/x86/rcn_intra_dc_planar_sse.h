#ifndef VVC_INTRA_PRED_SSE_H
#define VVC_INTRA_PRED_SSE_H

#include "stdint.h"
#include "stddef.h"

void vvc_intra_dc_pdpc_sse(const uint16_t *const src_above,
                           const uint16_t *const src_left,
                           uint16_t *const dst, ptrdiff_t dst_stride,
                           int log2_width, int log2_height);

void vvc_intra_planar_pdpc_sse(const uint16_t *const src_above,
                               const uint16_t *const src_left,
                               uint16_t *const dst, ptrdiff_t dst_stride,
                               int log2_width, int log2_height);

void vvc_intra_planar_pdpc_2_sse(const uint16_t *const src_above,
                                 const uint16_t *const src_left,
                                 uint16_t *const dst, ptrdiff_t dst_stride,
                                 int log2_width, int log2_height);

#endif
