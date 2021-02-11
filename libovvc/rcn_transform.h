#ifndef RCN_TRANSFORM_DATA_H
#define RCN_TRANSFORM_DATA_H

#include <stddef.h>
#include <stdint.h>

static void
matrix_multiplication(const int16_t* src, const int16_t* const tr_matrix,
                      int16_t* dst, ptrdiff_t src_stride, int tr_size,
                      int num_lines, int num_columns, int shift);

void
vvc_inverse_dct_ii_2(const int16_t* src, int16_t* dst, ptrdiff_t src_stride,
                     int num_lines, int num_columns, int shift);

void
vvc_inverse_dct_ii_4(const int16_t* src, int16_t* dst, ptrdiff_t src_stride,
                     int num_lines, int num_columns, int shift);

void
vvc_inverse_dct_ii_8(const int16_t* src, int16_t* dst, ptrdiff_t src_stride,
                     int num_lines, int num_columns, int shift);

void
vvc_inverse_dct_ii_16(const int16_t* src, int16_t* dst, ptrdiff_t src_stride,
                      int num_lines, int num_columns, int shift);

void
vvc_inverse_dct_ii_32(const int16_t* src, int16_t* dst, ptrdiff_t src_stride,
                      int num_lines, int num_columns, int shift);

void
vvc_inverse_dct_ii_64(const int16_t* src, int16_t* dst, ptrdiff_t src_stride,
                      int num_lines, int num_columns, int shift);

void
vvc_inverse_dct_viii_4(const int16_t* src, int16_t* dst, ptrdiff_t src_stride,
                       int num_lines, int num_columns, int shift);

void
vvc_inverse_dct_viii_8(const int16_t* src, int16_t* dst, ptrdiff_t src_stride,
                       int num_lines, int num_columns, int shift);

void
vvc_inverse_dct_viii_16(const int16_t* src, int16_t* dst, ptrdiff_t src_stride,
                        int num_lines, int num_columns, int shift);

void
vvc_inverse_dct_viii_32(const int16_t* src, int16_t* dst, ptrdiff_t src_stride,
                        int num_lines, int num_columns, int shift);

void
vvc_inverse_dst_vii_4(const int16_t* src, int16_t* dst, ptrdiff_t src_stride,
                      int num_lines, int num_columns, int shift);

void
vvc_inverse_dst_vii_8(const int16_t* src, int16_t* dst, ptrdiff_t src_stride,
                      int num_lines, int num_columns, int shift);

void
vvc_inverse_dst_vii_16(const int16_t* src, int16_t* dst, ptrdiff_t src_stride,
                       int num_lines, int num_columns, int shift);

void
vvc_inverse_dst_vii_32(const int16_t* src, int16_t* dst, ptrdiff_t src_stride,
                       int num_lines, int num_columns, int shift);

void
vvc_inverse_dct_ii_dc(int16_t* const dst, int log2_tb_w, int log2_tb_h,
                      int dc_val);

#endif // RCN_TRANSFORM_DATA_H
