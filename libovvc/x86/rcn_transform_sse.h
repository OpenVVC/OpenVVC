#ifndef VVC_TRANSFORM_SSE_H
#define VVC_TRANSFORM_SSE_H

#include <stdint.h>
#include <stddef.h>
#include "rcn_transform.h"
#include "rcn_structures.h"

void rcn_init_tr_functions_sse(struct RCNFunctions *const rcn_funcs);

void vvc_add_residual_sse(const int16_t *const src, uint16_t *const dst,
                          int log2_tb_w, int log2_tb_h, int scale);

// static void (*vvc_add_chroma_sse[3]) (const int16_t *const src, uint16_t *const dst,
//                                       int log2_tb_w, int log2_tb_h,
//                                       int scale)=
// {
//     &vvc_add_residual_sse,
//     &vvc_add_residual_sse,
//     &vvc_add_half_residual
// };
//
// static void (*vvc_scale_add_chroma_sse[3]) (const int16_t *const src, uint16_t *const dst,
//                                             int log2_tb_w, int log2_tb_h,
//                                             int scale)=
// {
//     &vvc_lmcs_add_residual,
//     &vvc_lmcs_add_residual,
//     &vvc_lmcs_add_half_residual
// };
//
// static void (*vvc_sub_chroma_sse[3]) (const int16_t *const src, uint16_t *const dst,
//                                       int log2_tb_w, int log2_tb_h,
//                                       int scale)=
// {
//     &vvc_add_residual_sse,
//     &vvc_sub_residual,
//     &vvc_sub_half_residual
// };
//
// static void (*vvc_scale_sub_chroma_sse[3]) (const int16_t *const src, uint16_t *const dst,
//                                             int log2_tb_w, int log2_tb_h,
//                                             int scale)=
// {
//     &vvc_lmcs_add_residual,
//     &vvc_lmcs_sub_residual,
//     &vvc_lmcs_sub_half_residual
// };
#endif//VVC_TRANSFORM_SSE_H
