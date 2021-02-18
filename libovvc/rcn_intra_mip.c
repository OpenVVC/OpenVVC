#include <stdint.h>

#include "data_rcn_mip.h"
#include "rcn_intra_mip.h"

const struct MIPCtx
derive_mip_ctx(int log2_pb_w, int log2_pb_h, uint8_t mip_mode)
{
    struct MIPCtx mip_ctx;
    const uint8_t *mip_matrix_selector = mip_weight_16x16;
    int matrix_size = 64*8;
    if (log2_pb_h == log2_pb_w && log2_pb_h == 2) { //4x4
        mip_matrix_selector = mip_weight_4x4;
        matrix_size = 16*4;
    } else if (log2_pb_h == 2 || log2_pb_w == 2 || (log2_pb_h <= 3 && log2_pb_w <= 3)) { //8x8
        mip_matrix_selector = mip_weight_8x8;
        matrix_size = 16*8;
    }

    mip_ctx.mip_matrix = &mip_matrix_selector[(int)mip_mode * matrix_size];

    return mip_ctx;
}

void
up_sample(uint16_t *const dst, const int16_t *const src,
          const uint16_t *ref,
          int log2_upsampled_size_src, int log2_opposite_size,
          int src_step, int src_stride,
          int dst_step, int dst_stride,
          int ref_step, int log2_scale)
{
    const int up_rounding_offset = 1 << (log2_scale - 1);
    const int16_t *src_line   = src;
    const uint16_t *bndy_line = ref + ref_step;
    uint16_t *dst_line  = dst;

    for (int i = 0; i < (1 << log2_opposite_size); ++i) {
        const int16_t *before   = bndy_line;
        const int16_t *after    = src_line;
        int16_t *curr_dst = dst_line;
        for (int j = 0; j < (1 << log2_upsampled_size_src); ++j) {
            int pos = 1;
            int32_t before_value = (*before) << log2_scale;
            int32_t after_value = 0;
            while (pos <= (1 << log2_scale)) {
                before_value -= *before;
                after_value  += *after;
                *curr_dst = (before_value + after_value + up_rounding_offset) >> log2_scale;
                curr_dst += dst_step;
                pos++;
            }
            before = after;
            after += src_step;
        }
        src_line  += src_stride;
        dst_line  += dst_stride;
        bndy_line += ref_step;
    }
}
