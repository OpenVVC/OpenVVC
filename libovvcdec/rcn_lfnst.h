#ifndef RCN_LFNST_H
#define RCN_LFNST_H

void
compute_lfnst_4x4(const int16_t* const src, uint16_t* const dst,
                  const int8_t* const lfnst_matrix, int log2_tb_w,
                  int log2_tb_h);
void
compute_lfnst_8x8(const int16_t* const src, uint16_t* const dst,
                  const int8_t* const lfnst_matrix, int log2_tb_w,
                  int log2_tb_h);
void
compute_lfnst_4x4_tr(const int16_t* const src, uint16_t* const dst,
                     const int8_t* const lfnst_matrix, int log2_tb_w,
                     int log2_tb_h);
void
compute_lfnst_8x8_tr(const int16_t* const src, uint16_t* const dst,
                     const int8_t* const lfnst_matrix, int log2_tb_w,
                     int log2_tb_h);

#endif // RCN_LFNST_H
