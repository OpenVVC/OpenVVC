#include <stdint.h>

#include "data_rcn_mip.h"
#include "rcn_intra_mip.h"
#include "rcn_fill_ref.h"
#include "ovutils.h"
#include "ctudec.h"

#if 0
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
#else

#define MIP_SHIFT 6
#define MIP_OFFSET (1 << (MIP_SHIFT - 1))

struct MIPCtx{
    const uint8_t *mip_matrix;
};

static const struct MIPCtx
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

static inline void
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
        const int16_t *before   = (int16_t*)bndy_line;
        const int16_t *after    = src_line;
        int16_t *curr_dst = (int16_t*)dst_line;
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


void
vvc_intra_pred_mip(const struct OVRCNCtx *const rcn_ctx,
                   uint16_t *const dst,
                   int x0, int y0, int log2_pu_w, int log2_pu_h,
                   uint8_t mip_mode)
{
    const uint16_t *src = &rcn_ctx->ctu_buff.y[0];

    int32_t bndy_line[8];//buffer used to store averaged boundaries use int
    int16_t mip_pred[64];//buffer used to store reduced matrix vector results

    /* FIXME determine max size of those buffers */
    uint16_t ref_abv[(128<<1) + 128];
    uint16_t ref_lft[(128<<1) + 128];

    int dst_stride = RCN_CTB_STRIDE;

    //compute reduced boundaries
    uint8_t log2_bndy = 1 << ((log2_pu_w > 2) || (log2_pu_h > 2));
    uint8_t log2_bnd_x = log2_pu_w - log2_bndy;
    uint8_t log2_bnd_y = log2_pu_h - log2_bndy;
    int i, j;

    int rnd = (1 << log2_bnd_x) >> 1;

    fill_ref_left_0(src, dst_stride, ref_lft,
                    rcn_ctx->progress_field.vfield[x0 >> 2],
                    rcn_ctx->progress_field.hfield[y0 >> 2],
                    x0, y0, log2_pu_w, log2_pu_h, 0);
    fill_ref_above_0(src, dst_stride, ref_abv,
                     rcn_ctx->progress_field.hfield[y0 >> 2],
                     rcn_ctx->progress_field.vfield[x0 >> 2],
                     x0, y0, log2_pu_w, log2_pu_h, 0);

    for (j = 0; j < (1 << log2_bndy); ++j) {
        int sum = 0;
        for (i = 0; i < (1 << log2_bnd_x); ++i) {
            sum += ref_abv[i + 1 + (j << log2_bnd_x)];
        }
        bndy_line[j] = (sum + rnd) >> log2_bnd_x;
    }

    rnd = (1 << log2_bnd_y) >> 1;
    for (j = 0; j < (1 << log2_bndy); ++j) {
        int sum = 0;
        for (i = 0; i < (1 << log2_bnd_y); ++i) {
            sum += ref_lft[i + 1 + (j << log2_bnd_y)];
        }
        bndy_line[(1 << log2_bndy) + j] = (sum + rnd) >> log2_bnd_y;
    }

    int16_t input_offset = bndy_line[0];

    uint8_t red_size = log2_pu_h == 2 || log2_pu_w == 2 || (log2_pu_h <= 3 && log2_pu_w <= 3);

    if (red_size) {
        bndy_line[0] = (1 << (10 - 1));
    }

    int sum = 0;
    for (i = 0; i < (2 << log2_bndy); ++i) {
        bndy_line[i] -= input_offset;
        sum += bndy_line[i];
    }

    //compute matrix multiplication
    const int rnd_mip = MIP_OFFSET - MIP_OFFSET * sum;

    uint8_t log2_red_w;
    uint8_t log2_red_h;

    if (red_size) { // 8x8 => 4x4
        log2_red_w = 2;
        log2_red_h = 2;
    } else { //saturate to 8
        log2_red_w = OVMIN(3, log2_pu_w);
        log2_red_h = OVMIN(3, log2_pu_h);
    }

    // if 4x16 bndy_size = 8 but need to skip some lines since 16 is reduced
    //(output on 4 * 8 =>32 instead of 8 * 8
    const int stride_x = 2 << log2_bndy;
    const struct MIPCtx mip_ctx = derive_mip_ctx(log2_pu_w, log2_pu_h, mip_mode);

    const uint8_t *matrix_mip = mip_ctx.mip_matrix;

    int x, y;
    int pos = 0;

    for (y = 0; y < (1 << log2_red_h); y++) {
        for (x = 0; x < (1 << log2_red_w); x++) {
            int val;
            int tmp0 = bndy_line[0] * matrix_mip[0];
            int tmp1 = bndy_line[1] * matrix_mip[1];
            int tmp2 = bndy_line[2] * matrix_mip[2];
            int tmp3 = bndy_line[3] * matrix_mip[3];
            for (i = 4; i < (2 << log2_bndy); i += 4) {
                tmp0 += bndy_line[i    ] * matrix_mip[i    ];
                tmp1 += bndy_line[i + 1] * matrix_mip[i + 1];
                tmp2 += bndy_line[i + 2] * matrix_mip[i + 2];
                tmp3 += bndy_line[i + 3] * matrix_mip[i + 3];
            }
            val = (tmp0 + tmp1) + (tmp2 + tmp3);
            mip_pred[pos++] = ov_clip(((val + rnd_mip) >> MIP_SHIFT) + input_offset, 0, 1023);
            matrix_mip += stride_x;
        }
    }

    // compute up_sampling
    uint8_t log2_scale_x = log2_pu_w - log2_red_w;
    uint8_t log2_scale_y = log2_pu_h - log2_red_h;

    if (log2_scale_x || log2_scale_y) {
        int src_stride;
        int src_step;
        //width then height
        const uint16_t *src;
        if (log2_scale_x) {
            uint16_t *_dst = dst + ((1 << log2_scale_y) - 1) * RCN_CTB_STRIDE;
            up_sample(_dst, mip_pred, ref_lft, log2_red_w, log2_red_h,
                       1, (1 << log2_red_w),
                       1, (1 << log2_scale_y) * RCN_CTB_STRIDE,
                       (1 << log2_scale_y), log2_scale_x);
            src        = _dst;
            src_step   = (1 << log2_scale_y) * RCN_CTB_STRIDE;
            src_stride = 1;
        } else { //TODO use mip_pred directly in next ste
            src        = mip_pred;
            src_step   = (1 << log2_pu_w);
            src_stride = 1;
        }
        up_sample(dst, src, ref_abv, log2_red_h, log2_pu_w,
                   src_step, src_stride,
                   RCN_CTB_STRIDE, 1,
                   1, log2_scale_y);

    } else {//write to dst
        for (i = 0; i < (1 << log2_red_h); ++i) {
            for (j = 0; j < (1 << log2_red_w); ++j) {
                dst [j + i * RCN_CTB_STRIDE] = mip_pred[(i << log2_red_w) + j];
            }
        }
    }
}

void
vvc_intra_pred_mip_tr(const struct OVRCNCtx *const rcn_ctx,
                      uint16_t *const dst,
                      int x0, int y0, int log2_pu_w, int log2_pu_h,
                      uint8_t mip_mode)
{
    uint8_t log2_bndy = 1 << ((log2_pu_w > 2) || (log2_pu_h > 2));

    uint8_t log2_red_w;
    uint8_t log2_red_h;

    int rnd;
    int i, j;
    const uint16_t *src = &rcn_ctx->ctu_buff.y[0];

    int32_t bndy_line[8]; //buffer used to store averaged boundaries use int
    int16_t mip_pred[64];//buffer used to store reduced matrix vector results

    uint16_t ref_abv[(128<<1) + 128];
    uint16_t ref_lft[(128<<1) + 128];

    int dst_stride = RCN_CTB_STRIDE;

    fill_ref_left_0(src, dst_stride, ref_lft,
                    rcn_ctx->progress_field.vfield[x0 >> 2],
                    rcn_ctx->progress_field.hfield[y0 >> 2],
                    x0, y0, log2_pu_w, log2_pu_h, 0);

    fill_ref_above_0(src, dst_stride, ref_abv,
                     rcn_ctx->progress_field.hfield[y0 >> 2],
                     rcn_ctx->progress_field.vfield[x0 >> 2],
                     x0, y0, log2_pu_w, log2_pu_h, 0);

    uint8_t log2_bnd_x = log2_pu_w - log2_bndy;
    uint8_t log2_bnd_y = log2_pu_h - log2_bndy;
    //compute reduced boundaries
    rnd = (1 << log2_bnd_x) >> 1;
    for (j = 0; j < (1 << log2_bndy); ++j) {
        int sum = 0;
        for (i = 0; i < (1 << log2_bnd_x); ++i) {
            sum += ref_abv[1 + i + (j << log2_bnd_x)];
        }
        bndy_line[(1 << log2_bndy) + j] = (sum + rnd) >> log2_bnd_x;
    }

    rnd = (1 << log2_bnd_y) >> 1;
    for (j = 0; j < (1 << log2_bndy); ++j) {
        int sum = 0;
        for (i = 0; i < (1 << log2_bnd_y); ++i) {
            sum += ref_lft[i + 1 + (j << log2_bnd_y)];
        }
        bndy_line[j] = (sum + rnd) >> log2_bnd_y;
    }

    int16_t input_offset = bndy_line[0];

    uint8_t red_size = log2_pu_h == 2 || log2_pu_w == 2 || (log2_pu_h <= 3 && log2_pu_w <= 3);

    if (red_size) {
        bndy_line[0] = (1 << (10 - 1));
    }

    int sum = 0;
    for (i = 0; i < (2 << log2_bndy); ++i) {
        bndy_line[i] -= input_offset;
        sum += bndy_line[i];
    }

    const int rnd_mip = MIP_OFFSET - MIP_OFFSET * sum;

    if (red_size) {
        log2_red_w = 2;
        log2_red_h = 2;
    } else {
        log2_red_w = OVMIN(3, log2_pu_w);
        log2_red_h = OVMIN(3, log2_pu_h);
    }

    // if 4x16 bndy_size = 8 but need to skip some lines since 16 is reduced
    //(output on 4 * 8 =>32 instead of 8 * 8
    const int stride_x = 2 << log2_bndy;
    const struct MIPCtx mip_ctx = derive_mip_ctx(log2_pu_w, log2_pu_h, mip_mode);

    const uint8_t *matrix_mip = mip_ctx.mip_matrix;
    int x, y;
    int pos = 0;

    for (y = 0; y < (1 << log2_red_w); y++) {
        for (x = 0; x < (1 << log2_red_h); x++) {
            int val;
            int tmp0 = bndy_line[0] * matrix_mip[0];
            int tmp1 = bndy_line[1] * matrix_mip[1];
            int tmp2 = bndy_line[2] * matrix_mip[2];
            int tmp3 = bndy_line[3] * matrix_mip[3];
            for (i = 4; i < (2 << log2_bndy); i += 4) {
                tmp0 += bndy_line[i    ] * matrix_mip[i    ];
                tmp1 += bndy_line[i + 1] * matrix_mip[i + 1];
                tmp2 += bndy_line[i + 2] * matrix_mip[i + 2];
                tmp3 += bndy_line[i + 3] * matrix_mip[i + 3];
            }
            val = (tmp0 + tmp1) + (tmp2 + tmp3);
            mip_pred[pos++] = ov_clip(((val + rnd_mip) >> MIP_SHIFT) + input_offset, 0, 1023);
            matrix_mip += stride_x;
        }
    }

    // compute up_sampling
    uint8_t log2_scale_x = log2_pu_w - log2_red_w;
    uint8_t log2_scale_y = log2_pu_h - log2_red_h;

    if (log2_scale_x || log2_scale_y) {
            int src_stride;
            int src_step;
            const uint16_t *src;
            if (log2_scale_x) {
                uint16_t *_dst = dst + ((1 << log2_scale_y) - 1) * RCN_CTB_STRIDE;
                up_sample(_dst, mip_pred, ref_lft,
                           log2_red_w, log2_red_h,
                           (1 << log2_red_h), 1,
                           1, (1 << log2_scale_y) * RCN_CTB_STRIDE,
                           (1 << log2_scale_y), log2_scale_x);
                src        = _dst;
                src_step   = (1 << log2_scale_y) * RCN_CTB_STRIDE;
                src_stride = 1;
            } else { //TODO use mip_pred directly in next ste
                src        = mip_pred;
                src_step   = 1;
                src_stride = 1 << log2_red_h;
            }
            up_sample(dst, src, ref_abv, log2_red_h, log2_pu_w,
                       src_step, src_stride,
                       RCN_CTB_STRIDE, 1,
                       1, log2_scale_y);

    } else {
        for (i = 0; i < (1 << log2_red_h); ++i) {
            for (j = 0; j < (1 << log2_red_w); ++j) {
                dst [j + i * RCN_CTB_STRIDE] = mip_pred[(j << log2_red_h) + i];
            }
        }
    }
}
#endif
