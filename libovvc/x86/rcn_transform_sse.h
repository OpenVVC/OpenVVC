#ifndef VVC_TRANSFORM_SSE_H
#define VVC_TRANSFORM_SSE_H

#include <stdint.h>
#include <stddef.h>
#include "rcn_transform.h"

void vvc_inverse_dct_ii_2_sse(const int16_t *src, int16_t *dst, ptrdiff_t src_stride,
                              int num_lines, int num_columns, int shift);

void vvc_inverse_dct_ii_4_sse(const int16_t *src, int16_t *dst, ptrdiff_t src_stride,
                              int num_lines, int num_columns, int shift);

void vvc_inverse_dct_ii_8_sse(const int16_t *src, int16_t *dst, ptrdiff_t src_stride,
                              int num_lines, int num_columns, int shift);

void vvc_inverse_dct_ii_16_sse(const int16_t *src, int16_t *dst, ptrdiff_t src_stride,
                               int num_lines, int num_columns, int shift);

void vvc_inverse_dct_ii_32_sse(const int16_t *src, int16_t *dst, ptrdiff_t src_stride,
                               int num_lines, int num_columns, int shift);

void vvc_inverse_dst_vii_4_sse(const int16_t *src, int16_t *dst, ptrdiff_t src_stride,
                              int num_lines, int num_columns, int shift);

void vvc_inverse_dst_vii_8_sse(const int16_t *src, int16_t *dst, ptrdiff_t src_stride,
                              int num_lines, int num_columns, int shift);

void vvc_inverse_dst_vii_16_sse(const int16_t *src, int16_t *dst, ptrdiff_t src_stride,
                               int num_lines, int num_columns, int shift);

void vvc_inverse_dst_vii_32_sse(const int16_t *src, int16_t *dst, ptrdiff_t src_stride,
                               int num_lines, int num_columns, int shift);

void vvc_inverse_dct_viii_4_sse(const int16_t *src, int16_t *dst, ptrdiff_t src_stride,
                              int num_lines, int num_columns, int shift);

void vvc_inverse_dct_viii_8_sse(const int16_t *src, int16_t *dst, ptrdiff_t src_stride,
                              int num_lines, int num_columns, int shift);

void vvc_inverse_dct_viii_16_sse(const int16_t *src, int16_t *dst, ptrdiff_t src_stride,
                               int num_lines, int num_columns, int shift);

void vvc_inverse_dct_viii_32_sse(const int16_t *src, int16_t *dst, ptrdiff_t src_stride,
                               int num_lines, int num_columns, int shift);

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

void vvc_inverse_dct_ii_dc_sse(uint16_t *const dst, int log2_tb_w, int log2_tb_h,
                               int dc_val);
#endif//VVC_TRANSFORM_SSE_H
