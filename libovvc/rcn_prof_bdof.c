#include <stdlib.h>

#include "rcn_prof_bdof.h"
#include "rcn_structures.h"
#include "ovutils.h"

#include "bitdepth.h"

#define MAX_PB_SIZE 128

#define PROF_SMP_SHIFT (14 - BITDEPTH)
#define PROF_PREC_RND (1 << (14 - 1))
#define PROF_DELTA_LIMIT (1 << 13)

#define BDOF_WGT_LIMIT ((1 << 4) - 1)
#define BDOF_SHIFT   (14 + 1 - BITDEPTH)
#define BDOF_OFFSET  ((1 << (BDOF_SHIFT - 1)))

#define GRAD_SHIFT 6


#define SB_H 4
#define SB_W 4

static void
rcn_apply_bdof_subblock(const int16_t* src0, int src0_stride,
                        const int16_t* src1, int src1_stride,
                        int16_t *dst, int dst_stride,
                        const int16_t *gradX0, const int16_t *gradX1,
                        const int16_t *gradY0, const int16_t *gradY1, int grad_stride,
                        int wgt_x, int wgt_y)
{
    int i;

    for (i = 0; i < SB_H; i++) {
        int32_t b0, b1, b2, b3;
        int16_t val0, val1, val2, val3;

        b0 = wgt_x * (gradX0[0] - gradX1[0]) + wgt_y * (gradY0[0] - gradY1[0]);
        b1 = wgt_x * (gradX0[1] - gradX1[1]) + wgt_y * (gradY0[1] - gradY1[1]);
        b2 = wgt_x * (gradX0[2] - gradX1[2]) + wgt_y * (gradY0[2] - gradY1[2]);
        b3 = wgt_x * (gradX0[3] - gradX1[3]) + wgt_y * (gradY0[3] - gradY1[3]);

        val0 = (int16_t)((src0[0] + src1[0] + b0 + BDOF_OFFSET) >> BDOF_SHIFT);
        val1 = (int16_t)((src0[1] + src1[1] + b1 + BDOF_OFFSET) >> BDOF_SHIFT);
        val2 = (int16_t)((src0[2] + src1[2] + b2 + BDOF_OFFSET) >> BDOF_SHIFT);
        val3 = (int16_t)((src0[3] + src1[3] + b3 + BDOF_OFFSET) >> BDOF_SHIFT);

        dst[0] = ov_bdclip(val0);
        dst[1] = ov_bdclip(val1);
        dst[2] = ov_bdclip(val2);
        dst[3] = ov_bdclip(val3);

        dst += dst_stride;

        src0 += src0_stride;
        src1 += src1_stride;

        gradX0 += grad_stride;
        gradX1 += grad_stride;

        gradY0 += grad_stride;
        gradY1 += grad_stride;
    }
}

void
tmp_prof_mrg(uint16_t* _dst, ptrdiff_t _dststride,
             const uint16_t* _src0, ptrdiff_t _srcstride,
             const int16_t* _src1, int height, intptr_t mx,
             intptr_t my, int width)
{
    int x, y;
    const int16_t* src0 = (int16_t *)_src0;
    const int16_t* src1 = (int16_t *)_src1;
    ptrdiff_t srcstride = _srcstride;
    uint16_t* dst = (uint16_t*)_dst;
    ptrdiff_t dststride = _dststride;
    int shift = 14 - BITDEPTH + 1;
    int offset = 2*((1 << (13 - BITDEPTH)));

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; ++x) {
            dst[x] = ov_bdclip((src0[x] + src1[x] + offset) >> shift);
        }
        src0 += srcstride;
        src1 += MAX_PB_SIZE;
        dst += dststride;
    }
}

void
tmp_prof_mrg_w(uint16_t* _dst, ptrdiff_t _dststride,
               const uint16_t* _src0, ptrdiff_t _srcstride,
               const int16_t* _src1, int height, intptr_t mx,
               intptr_t my, int width, int wt0, int wt1)
{
    int x, y;
    const int16_t* src0 = (int16_t *)_src0;
    const int16_t* src1 = (int16_t *)_src1;
    ptrdiff_t srcstride = _srcstride;
    uint16_t* dst = (uint16_t*)_dst;
    ptrdiff_t dststride = _dststride;
    int log_weights = floor_log2(wt0 + wt1);
    int shift = 14 - BITDEPTH + log_weights;
    int offset = 2*((1 << (13 - BITDEPTH))) << (log_weights - 1) ;

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; ++x) {
            dst[x] = ov_bdclip((src0[x] * wt0 + src1[x] * wt1 + offset) >> shift);
        }
        src0 += srcstride;
        src1 += MAX_PB_SIZE;
        dst += dststride;
    }
}

static void
compute_prof_grad(const uint16_t* src, int src_stride, int sb_w, int sb_h,
                  int grad_stride, int16_t* grad_x, int16_t* grad_y)
{
    int y, x;
    const int nb_smp_h = sb_h;
    const int nb_smp_w = sb_w;

    src += src_stride + 1;

    for (y = 0; y < nb_smp_h; ++y) {
        for (x = 0; x < nb_smp_w; ++x) {
            grad_y[x]  = (((int16_t)src[x + src_stride] - PROF_PREC_RND) >> GRAD_SHIFT);
            grad_y[x] -= (((int16_t)src[x - src_stride] - PROF_PREC_RND) >> GRAD_SHIFT);
            grad_x[x]  = (((int16_t)src[x + 1]          - PROF_PREC_RND) >> GRAD_SHIFT);
            grad_x[x] -= (((int16_t)src[x - 1]          - PROF_PREC_RND) >> GRAD_SHIFT);
        }
        grad_x += grad_stride;
        grad_y += grad_stride;
        src += src_stride;
    }
}

void
extend_prof_buff(const uint16_t *const src, uint16_t *dst_prof, int16_t ref_stride, uint8_t ext_x, uint8_t ext_y)
{
    const uint16_t *ref = src  - ref_stride  - 1;
    uint16_t       *dst = dst_prof;
    uint16_t *dst_lst = dst_prof + (SB_H + 1) * PROF_BUFF_STRIDE;
    int i, j;

    /* Position ref according to precision */
    if (ext_x) {
        ref += 1;
    }

    if (ext_y) {
        ref += ref_stride;
    }

    const uint16_t *ref_lst = ref + (SB_H + 1) * ref_stride;

    /* Copy or extend upper and lower ref_line */
    for (i = 0; i < SB_W + 2; ++i) {
        dst[i]     = (ref[i]     << PROF_SMP_SHIFT);
        dst_lst[i] = (ref_lst[i] << PROF_SMP_SHIFT);
    }

    dst += PROF_BUFF_STRIDE;
    dst_lst = dst + SB_W + 1;

    ref = src - 1;

    if (ext_x) {
        ref += 1;
    }

    if (ext_y) {
        ref += ref_stride;
    }

    ref_lst = ref + SB_W + 1;

    /* Copy or extend left and right column*/
    for (j = 0; j < SB_H; ++j) {
        dst[0]     = (ref[0]     << PROF_SMP_SHIFT);
        dst_lst[0] = (ref_lst[0] << PROF_SMP_SHIFT);

        ref     += ref_stride;
        ref_lst += ref_stride;
        dst     += PROF_BUFF_STRIDE;
        dst_lst += PROF_BUFF_STRIDE;
    }
}

static void
rcn_prof(uint16_t* dst, int dst_stride, const uint16_t* src, int src_stride,
         const int16_t* grad_x, const int16_t* grad_y, int grad_stride,
         const int32_t* dmv_scale_h, const int32_t* dmv_scale_v,
         uint8_t bidir)
{
    int idx = 0;
    int x, y;

    for (y = 0; y < SB_H; ++y) {
        for (x = 0; x < SB_W; ++x) {
            int32_t add = dmv_scale_h[idx] * grad_x[x] + dmv_scale_v[idx] * grad_y[x];
            int16_t val;

            add = ov_clip(add, -PROF_DELTA_LIMIT, PROF_DELTA_LIMIT - 1);

            val = (int16_t)src[x] + add;

            /* Clipping if not bi directional */
            if (!bidir) {
                val = (val + (1 << (13 - BITDEPTH))) >> PROF_SMP_SHIFT;
                dst[x] = ov_bdclip(val);
            } else {
                dst[x] = val;
            }

            idx++;
        }

        grad_x += grad_stride;
        grad_y += grad_stride;

        dst += dst_stride;
        src += src_stride;
    }
}

void
BD_DECL(rcn_init_prof_functions)(struct RCNFunctions *const rcn_funcs)
{
    rcn_funcs->prof.grad = &compute_prof_grad;
    rcn_funcs->prof.rcn  = &rcn_prof;
}

void
extend_bdof_buff(const uint16_t *const src, uint16_t *dst_prof,
                 int16_t ref_stride, int16_t pb_w, int16_t pb_h,
                 uint8_t ext_x, uint8_t ext_y)
{
    const uint16_t *ref = src  - ref_stride  - 1;

    uint16_t     *dst = dst_prof;
    uint16_t *dst_lst = dst_prof + (pb_h + 1) * PROF_BUFF_STRIDE;
    int i, j;

    /* Position ref according to precision */
    if (ext_x) {
        ref += 1;
    }

    if (ext_y) {
        ref += ref_stride;
    }

    const uint16_t *ref_lst = ref + (pb_h + 1) * ref_stride;

    /* Copy or extend upper and lower ref_line */
    for (i = 0; i < pb_w + 2; ++i) {
        dst[i]     = (ref[i]     << PROF_SMP_SHIFT);
        dst_lst[i] = (ref_lst[i] << PROF_SMP_SHIFT);
    }

    dst += PROF_BUFF_STRIDE;
    dst_lst = dst + pb_w + 1;

    ref = src - 1;

    if (ext_x) {
        ref += 1;
    }

    if (ext_y) {
        ref += ref_stride;
    }

    ref_lst = ref + pb_w + 1;

    /* Copy or extend left and right column*/
    for (j = 0; j < pb_h; ++j) {
        dst[0]     = (ref[0]     << PROF_SMP_SHIFT);
        dst_lst[0] = (ref_lst[0] << PROF_SMP_SHIFT);

        ref     += ref_stride;
        ref_lst += ref_stride;
        dst     += PROF_BUFF_STRIDE;
        dst_lst += PROF_BUFF_STRIDE;
    }
}


static void
derive_bdof_weights(const int16_t* ref0, const int16_t* ref1,
                    const int16_t* grad_x0, const int16_t* grad_x1,
                    const int16_t* grad_y0, const int16_t* grad_y1,
                    const int src0_stride, const int src1_stride,
                    const int grad_stride,
                    int *weight_x, int *weight_y)
{
    int sum_avg_x = 0;
    int sum_avg_y = 0;

    int sum_delta_x = 0;
    int sum_delta_y = 0;
    int wgt_x = 0;
    int wgt_y = 0;

    /* FIXME understand this part */
    /* sum / substract avg_grad_x based on avg_grad_y signs*/
    int sum_avg_x_y_signs = 0;

    int i, j;

    for (i = 0; i < SB_H + 2; i++) {
        for (j = 0; j < SB_W + 2; j++) {
            int32_t avg_grad_x = (grad_x0[j] + grad_x1[j]) >> 1;
            int32_t avg_grad_y = (grad_y0[j] + grad_y1[j]) >> 1;

            int32_t delta_ref = ((ref1[j] - PROF_PREC_RND) >> 4) - ((ref0[j] - PROF_PREC_RND) >> 4);

            sum_avg_x += abs(avg_grad_x);
            sum_avg_y += abs(avg_grad_y);

            sum_avg_x_y_signs += (avg_grad_y < 0 ? -avg_grad_x : (avg_grad_y == 0 ? 0 : avg_grad_x));

            sum_delta_x += (avg_grad_x < 0 ? -delta_ref : (avg_grad_x == 0 ? 0 : delta_ref));
            sum_delta_y += (avg_grad_y < 0 ? -delta_ref : (avg_grad_y == 0 ? 0 : delta_ref));
        }

        ref1 += src1_stride;
        ref0 += src0_stride;

        grad_x0 += grad_stride;
        grad_x1 += grad_stride;

        grad_y0 += grad_stride;
        grad_y1 += grad_stride;
    }

    if (sum_avg_x) {
        int log2_renorm_x = floor_log2(sum_avg_x);

        wgt_x = (sum_delta_x << 2) >> log2_renorm_x;
        wgt_x = ov_clip(wgt_x, -BDOF_WGT_LIMIT, BDOF_WGT_LIMIT);
        *weight_x = wgt_x;
    }

    if (sum_avg_y) {
        int log2_renorm_y = floor_log2(sum_avg_y);
        int x_offset = 0;

        if (wgt_x) {
            /* FIXME understand this part */
            int high = sum_avg_x_y_signs >> 12;
            int low  = sum_avg_x_y_signs & ((1 << 12) - 1);
            x_offset = (((wgt_x * high) << 12) + (wgt_x * low)) >> 1;
        }

        wgt_y = ((sum_delta_y << 2) - x_offset) >> log2_renorm_y;
        wgt_y = ov_clip(wgt_y, -BDOF_WGT_LIMIT, BDOF_WGT_LIMIT);
        *weight_y = wgt_y;
    }
}

void
rcn_bdof(struct BDOFFunctions *const bdof, int16_t *dst, int dst_stride,
         const int16_t *ref_bdof0, const int16_t *ref_bdof1, int ref_stride,
         const int16_t *grad_x0, const int16_t *grad_y0,
         const int16_t *grad_x1, const int16_t *grad_y1,
         int grad_stride, uint8_t pb_w, uint8_t pb_h)
{
    int nb_sb_w = (pb_w >> 2);
    int nb_sb_h = (pb_h >> 2);

    const int16_t *grad_x0_ln = grad_x0;
    const int16_t *grad_y0_ln = grad_y0;
    const int16_t *grad_x1_ln = grad_x1;
    const int16_t *grad_y1_ln = grad_y1;

    const int16_t *ref0_ln = ref_bdof0 - 128 - 1;
    const int16_t *ref1_ln = ref_bdof1 - 128 - 1;
    int16_t *dst_ln = dst;

    int i, j;

    for (i = 0; i < nb_sb_h; i++) {
        const int16_t *ref0_tmp = ref0_ln;
        const int16_t *ref1_tmp = ref1_ln;

        grad_x0 = grad_x0_ln;
        grad_y0 = grad_y0_ln;
        grad_x1 = grad_x1_ln;
        grad_y1 = grad_y1_ln;

        dst = dst_ln;

        for (j = 0; j < nb_sb_w; j++) {
            int wgt_x = 0;
            int wgt_y = 0;

            derive_bdof_weights(ref0_tmp, ref1_tmp,
                                grad_x0, grad_x1, grad_y0, grad_y1,
                                ref_stride, ref_stride,
                                grad_stride,
                                &wgt_x, &wgt_y);

            bdof->subblock(ref0_tmp + ref_stride + 1, ref_stride,
                           ref1_tmp + ref_stride + 1, ref_stride,
                           dst, dst_stride,
                           grad_x0 + grad_stride + 1, grad_x1 + grad_stride + 1,
                           grad_y0 + grad_stride + 1, grad_y1 + grad_stride + 1,
                           grad_stride,
                           wgt_x, wgt_y);

            grad_x0 += 1 << 2;
            grad_x1 += 1 << 2;
            grad_y0 += 1 << 2;
            grad_y1 += 1 << 2;
            ref0_tmp += 1 << 2;
            ref1_tmp += 1 << 2;
            dst    += 1 << 2;
        }

        grad_x0_ln += grad_stride << 2;
        grad_y0_ln += grad_stride << 2;
        grad_x1_ln += grad_stride << 2;
        grad_y1_ln += grad_stride << 2;

        ref0_ln += ref_stride << 2;
        ref1_ln += ref_stride << 2;

        dst_ln += dst_stride << 2;
    }
}

void
BD_DECL(rcn_init_bdof_functions)(struct RCNFunctions *const rcn_funcs)
{
    rcn_funcs->bdof.grad     = &compute_prof_grad;
    rcn_funcs->bdof.subblock = &rcn_apply_bdof_subblock;
}
