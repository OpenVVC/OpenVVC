#ifndef RCN_TRANSFORM_H
#define RCN_TRANSFORM_H

#include <stddef.h>
#include <stdint.h>
#include "rcn_structures.h"


void rcn_init_tr_functions(struct RCNFunctions *const rcn_funcs);

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
vvc_inverse_dct_ii_dc(int16_t* const dst, int log2_tb_w, int log2_tb_h,
                      int dc_val);

#endif // RCN_TRANSFORM_H
