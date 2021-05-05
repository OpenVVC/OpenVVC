#include <string.h>

#include "ovutils.h"

#include "rcn_structures.h"

static void
compute_lfnst_4x4(const int16_t* const src, int16_t* const dst,
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

static void
compute_lfnst_8x8(const int16_t* const src, int16_t* const dst,
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
                dst[(i & 3) + ((4 + ((i - 32) >> 2)) << log2_tb_w)] =
                  ov_clip((sum + 64) >> 7, -(1 << 15), (1 << 15));
        }
}

static void
compute_lfnst_4x4_tr(const int16_t* const src, int16_t* const dst,
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

static void
compute_lfnst_8x8_tr(const int16_t* const src, int16_t* const dst,
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

void
rcn_init_lfnst_functions(struct RCNFunctions *rcn_func)
{
   rcn_func->lfnst.func[0][0] = &compute_lfnst_4x4;
   rcn_func->lfnst.func[0][1] = &compute_lfnst_8x8;
   rcn_func->lfnst.func[1][0] = &compute_lfnst_4x4_tr;
   rcn_func->lfnst.func[1][1] = &compute_lfnst_8x8_tr;
}
