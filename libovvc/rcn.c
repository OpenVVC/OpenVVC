#include <stdint.h>
#include <string.h>
#include "nvcl_utils.h"
#include "ovutils.h"
#include "rcn_fill_ref.h"
#include "rcn_intra_angular.h"
#include "rcn_intra_dc_planar.h"
#include "data_rcn_angular.h"
#include "ctudec.h"
#include "rcn_intra_mip.h"
#include "rcn_structures.h"
#include "rcn_mc.h"
#include "rcn.h"
#include "ovmem.h"
#include "ovconfig.h"
#include "drv.h"

#if SSE_ENABLED
#include "x86/rcn_sse.h"
#endif

static int
derive_wide_angular_mode(int log2_pb_w, int log2_pb_h, int pred_mode)
{
    static const uint8_t mode_shift_tab[6] = {0, 6, 10, 12, 14, 15};
    int mode_shift = mode_shift_tab[OVABS(log2_pb_w - log2_pb_h)];

    if (log2_pb_w > log2_pb_h && pred_mode < 2 + mode_shift) {
        pred_mode += (OVINTRA_VDIA - 1);
    } else if (log2_pb_h > log2_pb_w && pred_mode > OVINTRA_VDIA - mode_shift) {
        pred_mode -= (OVINTRA_VDIA - 1);
    }
    return pred_mode;
}

void
vvc_intra_pred_isp(const OVCTUDec *const ctudec,
                   uint16_t *const src,
                   ptrdiff_t dst_stride,
                   uint8_t intra_mode,
                   int x0, int y0,
                   int log2_pb_w, int log2_pb_h,
                   int log2_cb_w, int log2_cb_h,
                   int offset_x, int offset_y)
{

    uint16_t ref_abv[(128<<1) + 128];
    uint16_t ref_lft[(128<<1) + 128];
    uint16_t *dst = &src[x0 + (y0 * dst_stride)];
    uint16_t *ref1 = ref_abv + (1 << log2_pb_h);
    uint16_t *ref2 = ref_lft + (1 << log2_pb_w);
    const struct OVRCNCtx *const rcn_ctx = &ctudec->rcn_ctx;
    const struct DCFunctions *dc = &rcn_ctx->rcn_funcs.dc;
    const struct PlanarFunctions *planar = &rcn_ctx->rcn_funcs.planar;

    fill_ref_left_0(src, dst_stride, ref2,
                    ctudec->rcn_ctx.progress_field.vfield[(x0 >> 2) + !!(offset_x % 4)],
                    ctudec->rcn_ctx.progress_field.hfield[((y0 ) >> 2) + !!(y0 % 4)],
                    x0, y0 - offset_y, log2_cb_w, log2_cb_h, offset_y);

    fill_ref_above_0(src, dst_stride, ref1,
                     ctudec->rcn_ctx.progress_field.hfield[(y0 >> 2) + !!(offset_y % 4)],
                     ctudec->rcn_ctx.progress_field.vfield[((x0 ) >> 2) + !!(x0 % 4)],
                     x0 - offset_x, y0, log2_cb_w, log2_cb_h, offset_x);

    ref1 += offset_x;
    ref2 += offset_y;

    for (int i = 0; i < 4; ++i){
        ref2[(1 << log2_cb_h) + (1 << log2_pb_h) + 1 + i] = ref2[(1 << log2_cb_h) + (1 << log2_pb_h) + i];
    }

    for (int i = 0; i < 4 ; ++i){
        ref1[(1 << log2_cb_w) + (1 << log2_pb_w) + 1 + i] = ref1[(1 << log2_cb_w) + (1 << log2_pb_w) + i];
    }

    switch (intra_mode) {
    case OVINTRA_PLANAR:
    {
        if (log2_pb_h > 1) {
            planar->pdpc[log2_pb_w > 5 || log2_pb_h > 5](ref1, ref2, dst, dst_stride,
                                                         log2_pb_w, log2_pb_h);
        } else {
            planar->func(ref1, ref2, dst, dst_stride, log2_pb_w, log2_pb_h);
        }
        break;
    }
    case OVINTRA_DC:
    {
        if (log2_pb_h > 1) {
            dc->pdpc(ref1, ref2, dst, dst_stride, log2_pb_w, log2_pb_h);
        } else {
            dc->func(ref1, ref2, dst, dst_stride, log2_pb_w, log2_pb_h);
        }
        break;
    }
    default:
    {
        int pred_mode = derive_wide_angular_mode(log2_cb_w, log2_cb_h, intra_mode);

        int is_vertical = pred_mode >= OVINTRA_DIA ? 1 : 0;

        if(is_vertical) {
            int mode_idx = pred_mode - OVINTRA_VER;
            switch (mode_idx) {
            case 0:
                if (log2_pb_h > 1) {
                    vvc_intra_ver_pdpc(ref1, ref2, dst, dst_stride, log2_pb_w,
                                       log2_pb_h);
                } else {
                    vvc_intra_ver(ref1, ref2, dst, dst_stride, log2_pb_w,
                                  log2_pb_h);
                }
                break;
            case 16:
                vvc_intra_angular_vdia(ref1, ref2, dst, dst_stride,
                                       log2_pb_w, log2_pb_h);
                break;
            default:
                if (mode_idx < 0){
                    int abs_angle_val = angle_table[-mode_idx];
                    uint8_t req_frac = !!(abs_angle_val & 0x1F);
                    int pu_height = 1 << log2_pb_h;
                    int inv_angle = inverse_angle_table[-mode_idx];
                    int inv_angle_sum    = 256;

                    for ( int k = -1; k >= -pu_height; k-- ){
                        inv_angle_sum += inv_angle;
                        ref1[k] = ref2[OVMIN(inv_angle_sum >> 9,pu_height)];
                    }

                    if (!req_frac){
                        intra_angular_v_nofrac(ref1, dst, dst_stride,
                                               log2_pb_w, log2_pb_h,
                                               -abs_angle_val);
                    } else {
                        intra_angular_v_cubic(ref1, dst, dst_stride,
                                              log2_pb_w, log2_pb_h,
                                              -abs_angle_val);
                    }

                } else if ((mode_idx < 8 &&  OVMIN(2, log2_pb_h - (floor_log2(3*inverse_angle_table[mode_idx] - 2) - 8)) < 0) || log2_pb_h < 2){//FIXME check this
                    int abs_angle_val = angle_table[mode_idx];
                    uint8_t req_frac = !!(abs_angle_val & 0x1F);
                    if (!req_frac){
                        intra_angular_v_nofrac(ref1, dst, dst_stride,
                                                log2_pb_w, log2_pb_h,
                                                abs_angle_val);
                    } else {
                        intra_angular_v_cubic(ref1, dst, dst_stride,
                                              log2_pb_w, log2_pb_h,
                                              abs_angle_val);
                    }
                } else {
                    uint8_t req_frac = !!(angle_table[mode_idx] & 0x1F);
                    if (!req_frac){
                        intra_angular_v_nofrac_pdpc(ref1, ref2, dst, dst_stride,
                                                    log2_pb_w, log2_pb_h,
                                                    mode_idx);
                    } else {
                        intra_angular_v_cubic_pdpc(ref1, ref2, dst, dst_stride,
                                              log2_pb_w, log2_pb_h,
                                              mode_idx);
                    }
                }
                break;
            }
        } else {
            int mode_idx = -(pred_mode - OVINTRA_HOR);
            switch (mode_idx) {
            case 0:
                if (log2_pb_h > 1){
                    vvc_intra_hor_pdpc(ref1, ref2, dst, dst_stride,
                                       log2_pb_w, log2_pb_h);
                } else {
                    vvc_intra_hor(ref1, ref2, dst, dst_stride,
                                  log2_pb_w, log2_pb_h);
                }
                break;
            case 16:
                vvc_intra_angular_hdia(ref1, ref2, dst, dst_stride,
                                       log2_pb_w, log2_pb_h);
                break;
            default:
            {
                if (mode_idx < 0){
                    int abs_angle_val = angle_table[-mode_idx];
                    uint8_t req_frac = !!(abs_angle_val & 0x1F);
                    int pu_width = 1 << log2_pb_w;
                    int inv_angle = inverse_angle_table[-mode_idx];
                    int inv_angle_sum    = 256;

                    for ( int k = -1; k >= -pu_width; k-- ){
                        inv_angle_sum += inv_angle;
                        ref2[k] = ref1[OVMIN(inv_angle_sum >> 9, pu_width)];
                    }

                    if (!req_frac){
                        intra_angular_h_nofrac(ref2, dst, dst_stride,
                                               log2_pb_w, log2_pb_h,
                                               -abs_angle_val);
                    } else {
                        intra_angular_h_cubic(ref2, dst, dst_stride,
                                              log2_pb_w, log2_pb_h,
                                              -abs_angle_val);

                    }
                } else if ((mode_idx < 8 &&  OVMIN(2, log2_pb_w - (floor_log2(3*inverse_angle_table[mode_idx] - 2) - 8)) < 0) || log2_pb_h < 2){//FIXME check this
                    int abs_angle_val = angle_table[mode_idx];
                    uint8_t req_frac = !!(abs_angle_val & 0x1F);
                    if (!req_frac){
                        intra_angular_h_nofrac(ref2, dst, dst_stride,
                                                log2_pb_w, log2_pb_h,
                                                abs_angle_val);
                    } else {
                        intra_angular_h_cubic(ref2, dst, dst_stride,
                                              log2_pb_w, log2_pb_h,
                                              abs_angle_val);
                    }
                } else {
                    uint8_t req_frac = !!(angle_table[mode_idx] & 0x1F);
                    if (!req_frac){
                        intra_angular_h_nofrac_pdpc(ref1, ref2, dst, dst_stride,
                                                    log2_pb_w, log2_pb_h,
                                                    mode_idx);
                    } else {
                        intra_angular_h_cubic_pdpc(ref1, ref2, dst, dst_stride,
                                              log2_pb_w, log2_pb_h,
                                              mode_idx);
                    }
                }
            }
                break;
            }
        }
        break;
    }
    }
}

void
vvc_intra_pred_multi_ref(const OVCTUDec *const ctudec,
                         uint16_t *const src,
                         ptrdiff_t dst_stride,
                         uint8_t intra_mode, int x0, int y0,
                         int log2_pb_w, int log2_pb_h,
                         int mrl_idx)
{
    uint16_t ref_abv[(128<<1) + 128];
    uint16_t ref_lft [(128<<1) + 128];
    uint16_t *dst = &src[x0 + (y0 * dst_stride)];
    uint16_t *ref1 = ref_abv + (1 << log2_pb_h);
    uint16_t *ref2 = ref_lft  + (1 << log2_pb_w);
    const struct OVRCNCtx *const rcn_ctx = &ctudec->rcn_ctx;
    const struct DCFunctions *dc = &rcn_ctx->rcn_funcs.dc;
    const struct PlanarFunctions *planar = &rcn_ctx->rcn_funcs.planar;

    fill_ref_left_0_mref(src, dst_stride, ref2,
                         ctudec->rcn_ctx.progress_field.vfield[x0 >> 2],
                         ctudec->rcn_ctx.progress_field.hfield[y0 >> 2],
                         mrl_idx, x0, y0,
                         log2_pb_w, log2_pb_h);

    fill_ref_above_0_mref(src, dst_stride, ref1,
                          ctudec->rcn_ctx.progress_field.hfield[y0 >> 2],
                          ctudec->rcn_ctx.progress_field.vfield[x0 >> 2],
                          mrl_idx, x0 , y0,
                          log2_pb_w, log2_pb_h);

    ref1 += mrl_idx;
    ref2 += mrl_idx;

    switch (intra_mode) {
    case OVINTRA_PLANAR:
        planar->func(ref1, ref2, dst, dst_stride, log2_pb_w, log2_pb_h);
        break;
    case OVINTRA_DC:
        dc->func(ref1, ref2, dst, dst_stride, log2_pb_w, log2_pb_h);
        break;
    default:
    {
        int pred_mode = derive_wide_angular_mode(log2_pb_w, log2_pb_h, intra_mode);

        int is_vertical = pred_mode >= OVINTRA_DIA ? 1 : 0;

        if(is_vertical){
            int mode_idx = pred_mode - OVINTRA_VER ;
            switch (mode_idx) {
            case 0:
                vvc_intra_ver(ref1, ref2, dst, dst_stride,
                              log2_pb_w, log2_pb_h);
                break;
            case 16:
                intra_angular_vdia_mref(ref1, ref2, dst, dst_stride,
                                        log2_pb_w, log2_pb_h,
                                        mrl_idx);
                break;
            default:
            {
                int angle_val;
                if (mode_idx < 0){
                    int inv_angle = inverse_angle_table[-mode_idx];
                    int pb_h  = 1 << log2_pb_h;
                    int inv_angle_sum = 256;

                    uint16_t *dst_ref = ref1 - mrl_idx;
                    uint16_t *tmp_lft = ref2 - mrl_idx;

                    angle_val = -angle_table[-mode_idx];

                    /* FIXME two stage fill and broadcast last value */
                    for (int k = -1; k >= -pb_h; k--) {
                        inv_angle_sum += inv_angle;
                        dst_ref[k] = tmp_lft[OVMIN(inv_angle_sum >> 9, pb_h)];
                    }

                } else {
                    angle_val = angle_table[mode_idx];
                }

                intra_angular_v_cubic_mref(ref1, dst, dst_stride,
                                           log2_pb_w, log2_pb_h,
                                           angle_val, mrl_idx);

                break;
            }
            }
        } else {
            int mode_idx = -(pred_mode - OVINTRA_HOR);
            switch (mode_idx) {
            case 0:
                vvc_intra_hor(ref1, ref2, dst, dst_stride,
                              log2_pb_w, log2_pb_h);
                break;
            case 16:
                intra_angular_hdia_mref(ref1, ref2, dst, dst_stride,
                                        log2_pb_w, log2_pb_h,
                                        mrl_idx);
                break;
            default:
            {
                int angle_val;
                if (mode_idx < 0) {
                    int inv_angle = inverse_angle_table[-mode_idx];
                    int inv_angle_sum = 256;

                    uint16_t *dst_ref = ref2 - mrl_idx;
                    uint16_t *tmp_abv = ref1 - mrl_idx;

                    int pb_w = 1 << log2_pb_w;

                    angle_val = -angle_table[-mode_idx];

                    /* FIXME two stage fill and broadcast last value */
                    for (int k = -1; k >= -pb_w; k--) {
                        inv_angle_sum += inv_angle;
                        dst_ref[k] = tmp_abv[OVMIN((inv_angle_sum >> 9), pb_w)];
                    }

                } else {
                    angle_val = angle_table[mode_idx];
                }

                intra_angular_h_cubic_mref(ref2, dst, dst_stride,
                                           log2_pb_w, log2_pb_h,
                                           angle_val, mrl_idx);
                break;
            }
            }
        }
        break;
    }
    }
}

void
rcn_residual(OVCTUDec *const ctudec,
             int16_t *const dst, int16_t *src,
             uint8_t x0, uint8_t y0,
             unsigned int log2_tb_w, unsigned int log2_tb_h,
             unsigned int lim_cg_w,
             uint8_t cu_mts_flag, uint8_t cu_mts_idx,
             uint8_t is_dc, uint8_t lfnst_flag, uint8_t is_mip, uint8_t lfnst_idx)
{
    struct TRFunctions *TRFunc = &ctudec->rcn_ctx.rcn_funcs.tr;
    int shift_v = 6 + 1;
    int shift_h = (6 + 15 - 1) - 10;
    DECLARE_ALIGNED(32, int16_t, tmp)[64*64];
    // int16_t tmp[64*64];
    int tb_w = 1 << log2_tb_w;
    int tb_h = 1 << log2_tb_h;
    #if 0
    struct OVRCNCtx *const rcn = &ctudec->rcn_ctx;
    #endif
    /* TODO switch */

    memset(tmp, 0, sizeof(int16_t) << (log2_tb_w + log2_tb_h));

    if (lfnst_flag) {
        lim_cg_w = 8;
        is_dc = 0;
        /* FIXME separate lfnst mode derivation from lfnst reconstruction */
#if 1
        process_lfnst_luma(ctudec, src, ctudec->lfnst_subblock, log2_tb_w, log2_tb_h, x0, y0,
                           lfnst_idx);
#endif
    }

    if (!is_mip && ctudec->mts_implicit && (log2_tb_w <= 4 || log2_tb_h <= 4) && !lfnst_flag) {
        /*FIXME condition on size in the if could be removed ?*/
        enum DCTType tr_h_idx = log2_tb_w <= 4 ? DST_VII : DCT_II;
        enum DCTType tr_v_idx = log2_tb_h <= 4 ? DST_VII : DCT_II;

#if 0
        call(vvc_inverse_mts_func,
             (tr_v_idx, log2_tb_h),
             (src, tmp, tb_w, tb_w, tb_h, shift_v));
        call(vvc_inverse_mts_func,
             (tr_h_idx, log2_tb_w),
             (tmp, dst, tb_h, tb_h, tb_w, shift_h));
#endif
        TRFunc->func[tr_v_idx][log2_tb_h](src, tmp, tb_w, tb_w, tb_h, shift_v);
        TRFunc->func[tr_h_idx][log2_tb_w](tmp, dst, tb_h, tb_h, tb_w, shift_h);

    } else if (!cu_mts_flag) {

#if 1
        if (is_dc) {

            TRFunc->dc(dst, log2_tb_w, log2_tb_h, src[0]);

        } else {
#endif
            int nb_row = OVMIN(lim_cg_w, 1 << log2_tb_w);
            int nb_col = OVMIN(lim_cg_w, 1 << log2_tb_h);

            TRFunc->func[DCT_II][log2_tb_h](src, tmp, tb_w, nb_row, nb_col, shift_v);
            TRFunc->func[DCT_II][log2_tb_w](tmp, dst, tb_h, tb_h, nb_row, shift_h);
#if 1
        }
#endif
    } else {
        enum DCTType tr_h_idx = cu_mts_idx  & 1;
        enum DCTType tr_v_idx = cu_mts_idx >> 1;

        TRFunc->func[tr_v_idx][log2_tb_h](src, tmp, tb_w, tb_w, tb_h, shift_v);
        TRFunc->func[tr_h_idx][log2_tb_w](tmp, dst, tb_h, tb_h, tb_w, shift_h);

    }
}

void
rcn_residual_c(OVCTUDec *const ctudec,
               int16_t *const dst, int16_t *src,
               int16_t *const lfnst_sb,
               uint8_t x0, uint8_t y0,
               unsigned int log2_tb_w, unsigned int log2_tb_h,
               unsigned int lim_cg_w,
               uint8_t cu_mts_flag, uint8_t cu_mts_idx,
               uint8_t is_dc, uint8_t lfnst_flag, uint8_t is_mip, uint8_t lfnst_idx)
{
    struct TRFunctions *TRFunc = &ctudec->rcn_ctx.rcn_funcs.tr;
    const int shift_v = 6 + 1;
    const int shift_h = (6 + 15 - 1) - 10;
    DECLARE_ALIGNED(32, int16_t, tmp)[32*32];
    // int16_t tmp[32*32];
    int tb_w = 1 << log2_tb_w;
    int tb_h = 1 << log2_tb_h;

    memset(tmp, 0, sizeof(int16_t) << (log2_tb_w + log2_tb_h));

    if (lfnst_flag) {
        /* Update lim_cg_w since lfnst part of coeff are now non zero */
        lim_cg_w = 8;
        is_dc = 0;
        /* FIXME separate lfnst mode derivation from lfnst reconstruction */
#if 1
        process_lfnst(ctudec, src, lfnst_sb, log2_tb_w, log2_tb_h,
                      x0, y0, lfnst_idx);
#endif
    }

#if 1
    if (is_dc && !lfnst_flag) {

            TRFunc->dc(dst, log2_tb_w, log2_tb_h, src[0]);

    } else {
#endif
        int nb_row =  OVMIN(lim_cg_w, 1 << log2_tb_w);
        int nb_col =  OVMIN(lim_cg_w, 1 << log2_tb_h);
        /*FIXME might be transform SKIP */

        #if 0
        TRFunc->func[DCT_II][log2_tb_h](src, tmp, tb_w, tb_w, tb_h, shift_v);
        TRFunc->func[DCT_II][log2_tb_w](tmp, src, tb_h, tb_h, tb_w, shift_h);
        #endif
        TRFunc->func[DCT_II][log2_tb_h](src, tmp, tb_w, nb_row, nb_col, shift_v);
        TRFunc->func[DCT_II][log2_tb_w](tmp, dst, tb_h, tb_h, nb_row, shift_h);
#if 1
    }
#endif
}

void rcn_init_functions(struct RCNFunctions *rcn_func, uint8_t ict_type){
  rcn_init_mc_functions(rcn_func);
  rcn_init_tr_functions(rcn_func);
  rcn_init_dc_planar_functions(rcn_func);
  rcn_init_ict_functions(rcn_func, ict_type);
  rcn_init_alf_functions(rcn_func);
  #if ARCH_X86
    #if SSE_ENABLED
      rcn_init_mc_functions_sse(rcn_func);
      rcn_init_tr_functions_sse(rcn_func);
      rcn_init_dc_planar_functions_sse(rcn_func);
      rcn_init_ict_functions_sse(rcn_func, ict_type);
      rcn_init_alf_functions_sse(rcn_func);
    #elif AVX_ENABLED
      //Link AVX optims
    #else
      //Failover x86
    #endif
  #elif ARCH_ARM
    #if NEON_ENABLED
      //Link NEON optims
    #else
      //Failover ARM
    #endif
  #else
    //Failover other arch
  #endif
}
