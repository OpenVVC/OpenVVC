#include <string.h>

#include "ovutils.h"

#include "data_rcn_transform.h"
#include "rcn_lfnst.h"

void
compute_lfnst_4x4(const int16_t* const src, uint16_t* const dst,
                  const int8_t* const lfnst_matrix, int log2_tb_w,
                  int log2_tb_h)
{
        int i, j;
        for (i = 0; i < 16; ++i) {
                int sum = 0;
                for (j = 0; j < 8 << !(log2_tb_w == log2_tb_h); ++j) {
                        sum += src[j] * lfnst_matrix[i + j * 16];
                }
                dst[(i & 3) + ((i >> 2) << log2_tb_w)] =
                  ov_clip((sum + 64) >> 7, -(1 << 15), (1 << 15));
        }
}

void
compute_lfnst_8x8(const int16_t* const src, uint16_t* const dst,
                  const int8_t* const lfnst_matrix, int log2_tb_w,
                  int log2_tb_h)
{
        int i, j;
        for (i = 0; i < 32; ++i) {
                int sum = 0;
                for (j = 0;
                     j < 8 << !(log2_tb_w == log2_tb_h && log2_tb_w == 8);
                     ++j) {
                        sum += src[j] * lfnst_matrix[i + j * 48];
                }
                dst[(i & 7) + ((i >> 3) << log2_tb_w)] =
                  ov_clip((sum + 64) >> 7, -(1 << 15), (1 << 15));
        }
        for (; i < 48; ++i) {
                int sum = 0;
                for (j = 0;
                     j < 8 << !(log2_tb_w == log2_tb_h && log2_tb_w == 8);
                     ++j) {
                        sum += src[j] * lfnst_matrix[i + j * 48];
                }
                dst[(i & 3) + (4 + ((i - 32) >> 2) << log2_tb_w)] =
                  ov_clip((sum + 64) >> 7, -(1 << 15), (1 << 15));
        }
}

void
compute_lfnst_4x4_tr(const int16_t* const src, uint16_t* const dst,
                     const int8_t* const lfnst_matrix, int log2_tb_w,
                     int log2_tb_h)
{
        int i, j;
        for (i = 0; i < 16; ++i) {
                int sum = 0;
                for (j = 0; j < 8 << !(log2_tb_w == log2_tb_h); ++j) {
                        sum += src[j] * lfnst_matrix[i + j * 16];
                }
                dst[((i & 3) << log2_tb_w) + (i >> 2)] =
                  ov_clip((sum + 64) >> 7, -(1 << 15), (1 << 15));
        }
}

void
compute_lfnst_8x8_tr(const int16_t* const src, uint16_t* const dst,
                     const int8_t* const lfnst_matrix, int log2_tb_w,
                     int log2_tb_h)
{
        int i, j;
        for (i = 0; i < 32; ++i) {
                int sum = 0;
                for (j = 0;
                     j < 8 << !(log2_tb_w == log2_tb_h && log2_tb_w == 8);
                     ++j) {
                        sum += src[j] * lfnst_matrix[i + j * 48];
                }
                dst[((i & 7) << log2_tb_w) + (i >> 3)] =
                  ov_clip((sum + 64) >> 7, -(1 << 15), (1 << 15));
        }
        for (; i < 48; ++i) {
                int sum = 0;
                for (j = 0;
                     j < 8 << !(log2_tb_w == log2_tb_h && log2_tb_w == 8);
                     ++j) {
                        sum += src[j] * lfnst_matrix[i + j * 48];
                }
                dst[((i & 3) << log2_tb_w) + (4 + ((i - 32) >> 2))] =
                  ov_clip((sum + 64) >> 7, -(1 << 15), (1 << 15));
        }
}
