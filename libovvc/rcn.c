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

#if 1
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

void vvc_intra_pred_isp(const OVCTUDec *const ctudec,
                        uint16_t *const src,
                        ptrdiff_t dst_stride,
                        uint8_t intra_mode,
                        int x0, int y0,
                        int log2_pb_width, int log2_pb_height,
                        int log2_cu_width,int log2_cu_height,
                        int offset_x, int offset_y){

    uint16_t ref_above[(128<<1) + 128]/*={512}*/;
    uint16_t ref_left [(128<<1) + 128]/*={512}*/;
    uint16_t *dst = &src[x0 + (y0*dst_stride)];
    uint16_t *ref1 = ref_above + (1 << log2_pb_height);
    uint16_t *ref2 = ref_left + (1 << log2_pb_width);
    const struct OVRCNCtx *const rcn_ctx = &ctudec->rcn_ctx;
    const struct DCFunctions *dc = &rcn_ctx->rcn_funcs.dc;
    const struct PlanarFunctions *planar = &rcn_ctx->rcn_funcs.planar;

    // FIXED? :fill_ref_left_0(src, dst_stride, ref2,
    //                 ctudec->progress_map.cols[(x0 >> 2) + !!(offset_x % 4)],
    //                 ctudec->progress_map.rows[((y0 ) >> 2) + !!(y0 % 4)],
    //                 x0, y0 - offset_y, log2_cu_width, log2_cu_height, offset_y);
    fill_ref_left_0(src, dst_stride, ref2,
                    ctudec->rcn_ctx.progress_field.vfield[(x0 >> 2) + !!(offset_x % 4)],
                    ctudec->rcn_ctx.progress_field.hfield[((y0 ) >> 2) + !!(y0 % 4)],
                    x0, y0 - offset_y, log2_cu_width, log2_cu_height, offset_y);
    //  FIXED? :fill_ref_above_0(src, dst_stride, ref1,
    //                  ctudec->progress_map.rows[(y0 >> 2) + !!(offset_y % 4)],
    //                  ctudec->progress_map.cols[((x0 ) >> 2) + !!(x0 % 4)],
    //                  x0 - offset_x, y0, log2_cu_width, log2_cu_height, offset_x);
   fill_ref_above_0(src, dst_stride, ref1,
                    ctudec->rcn_ctx.progress_field.hfield[(y0 >> 2) + !!(offset_y % 4)],
                    ctudec->rcn_ctx.progress_field.vfield[((x0 ) >> 2) + !!(x0 % 4)],
                    x0 - offset_x, y0, log2_cu_width, log2_cu_height, offset_x);
    ref1+= offset_x;
    ref2+= offset_y;
    for (int i = 0; i < 4; ++i){
        ref2[(1 << log2_cu_height) + (1 << log2_pb_height) + 1 + i] = ref2[(1 << log2_cu_height) + (1 << log2_pb_height) + i];
    }
    for (int i = 0; i < 4 ; ++i){
        ref1[(1 << log2_cu_width) + (1 << log2_pb_width) + 1 + i] = ref1[(1 << log2_cu_width) + (1 << log2_pb_width) + i];
    }

    switch (intra_mode) {
    case OVINTRA_PLANAR://PLANAR
    {
        if (log2_pb_height > 1)
            planar->pdpc[log2_pb_width > 5 || log2_pb_height > 5](ref1, ref2, dst, dst_stride,
                                                                  log2_pb_width, log2_pb_height);
        else
            planar->func(ref1, ref2, dst, dst_stride, log2_pb_width, log2_pb_height);
        break;
    }
    case OVINTRA_DC://DC
    {
        if (log2_pb_height > 1)
            dc->pdpc(ref1, ref2, dst, dst_stride, log2_pb_width, log2_pb_height);
        else
            dc->func(ref1, ref2, dst, dst_stride, log2_pb_width, log2_pb_height);
        break;
    }
    default://angular
    {
        int pred_mode = derive_wide_angular_mode(log2_cu_width, log2_cu_height,
                                                 intra_mode);

        int is_vertical = pred_mode >= OVINTRA_DIA ? 1 : 0;

        if(is_vertical){
            int mode_idx = pred_mode - OVINTRA_VER;
            switch (mode_idx) {
            case 0:
                //pure vertical
                if (log2_pb_height > 1){
                    vvc_intra_ver_pdpc(ref1, ref2, dst, dst_stride, log2_pb_width,
                                       log2_pb_height);
                } else {
                    vvc_intra_ver(ref1, ref2, dst, dst_stride, log2_pb_width,
                                  log2_pb_height);
                }
                break;
            case (16)://Pure diagonal
                vvc_intra_angular_vdia(ref1, ref2, dst, dst_stride,
                                       log2_pb_width, log2_pb_height);
                break;
            default:
                if (mode_idx < 0){
                    int abs_angle_val = angle_table[-mode_idx];
                    uint8_t req_frac = !!(abs_angle_val & 0x1F);
                    int pu_height = 1 << log2_pb_height;
                    int inv_angle = inverse_angle_table[-mode_idx];
                    int inv_angle_sum    = 256;

                    for ( int k = -1; k >= -pu_height; k-- ){
                        inv_angle_sum += inv_angle;
                        ref1[k] = ref2[OVMIN(inv_angle_sum >> 9,pu_height)];
                    }

                    if (!req_frac){
                        intra_angular_v_nofrac(ref1, dst, dst_stride,
                                                log2_pb_width, log2_pb_height,
                                                -abs_angle_val);
                    } else {
                        intra_angular_v_cubic(ref1, dst, dst_stride,
                                              log2_pb_width, log2_pb_height,
                                              -abs_angle_val);
                    }

                } else if ((mode_idx < 8 &&  OVMIN(2, log2_pb_height - (floor_log2(3*inverse_angle_table[mode_idx] - 2) - 8)) < 0) || log2_pb_height < 2){//FIXME check this
                    int abs_angle_val = angle_table[mode_idx];
                    uint8_t req_frac = !!(abs_angle_val & 0x1F);
                    if (!req_frac){
                        intra_angular_v_nofrac(ref1, dst, dst_stride,
                                                log2_pb_width, log2_pb_height,
                                                abs_angle_val);
                    } else {
                        intra_angular_v_cubic(ref1, dst, dst_stride,
                                              log2_pb_width, log2_pb_height,
                                              abs_angle_val);
                    }
                } else {
                    uint8_t req_frac = !!(angle_table[mode_idx] & 0x1F);
                    if (!req_frac){
                        intra_angular_v_nofrac_pdpc(ref1, ref2, dst, dst_stride,
                                                    log2_pb_width, log2_pb_height,
                                                    mode_idx);
                    } else {
                        intra_angular_v_cubic_pdpc(ref1, ref2, dst, dst_stride,
                                              log2_pb_width, log2_pb_height,
                                              mode_idx);
                    }
                }
                break;
            }
        } else {
            int mode_idx = -(pred_mode - OVINTRA_HOR);
            switch (mode_idx) {
            case 0:
                if (log2_pb_height > 1){
                    vvc_intra_hor_pdpc(ref1, ref2, dst, dst_stride,
                                       log2_pb_width, log2_pb_height);
                } else {
                    vvc_intra_hor(ref1, ref2, dst, dst_stride,
                                  log2_pb_width, log2_pb_height);
                }
                break;
            case (16):
                vvc_intra_angular_hdia(ref1, ref2, dst, dst_stride,
                                       log2_pb_width, log2_pb_height);
                break;
            default:
            {
                if (mode_idx < 0){
                    int abs_angle_val = angle_table[-mode_idx];
                    uint8_t req_frac = !!(abs_angle_val & 0x1F);
                    int pu_width = 1 << log2_pb_width;
                    int inv_angle = inverse_angle_table[-mode_idx];
                    int inv_angle_sum    = 256;

                    for ( int k = -1; k >= -pu_width; k-- ){
                        inv_angle_sum += inv_angle;
                        ref2[k] = ref1[OVMIN(inv_angle_sum >> 9, pu_width)];
                    }

                    if (!req_frac){
                        intra_angular_h_nofrac(ref2, dst, dst_stride,
                                               log2_pb_width, log2_pb_height,
                                               -abs_angle_val);
                    } else {
                        intra_angular_h_cubic(ref2, dst, dst_stride,
                                              log2_pb_width, log2_pb_height,
                                              -abs_angle_val);

                    }
                } else if ((mode_idx < 8 &&  OVMIN(2, log2_pb_width - (floor_log2(3*inverse_angle_table[mode_idx] - 2) - 8)) < 0) || log2_pb_height < 2){//FIXME check this
                    int abs_angle_val = angle_table[mode_idx];
                    uint8_t req_frac = !!(abs_angle_val & 0x1F);
                    if (!req_frac){
                        intra_angular_h_nofrac(ref2, dst, dst_stride,
                                                log2_pb_width, log2_pb_height,
                                                abs_angle_val);
                    } else {
                        intra_angular_h_cubic(ref2, dst, dst_stride,
                                              log2_pb_width, log2_pb_height,
                                              abs_angle_val);
                    }
                } else {
                    uint8_t req_frac = !!(angle_table[mode_idx] & 0x1F);
                    if (!req_frac){
                        intra_angular_h_nofrac_pdpc(ref1, ref2, dst, dst_stride,
                                                    log2_pb_width, log2_pb_height,
                                                    mode_idx);
                    } else {
                        intra_angular_h_cubic_pdpc(ref1, ref2, dst, dst_stride,
                                              log2_pb_width, log2_pb_height,
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

void vvc_intra_pred_multi_ref( const OVCTUDec *const ctudec,
                               uint16_t *const src,
                               ptrdiff_t dst_stride,
                               uint8_t intra_mode, int x0, int y0,
                               int log2_pb_width, int log2_pb_height,
                               int multi_ref_idx){

    uint16_t ref_above[(128<<1) + 128]/*={0}*/;
    uint16_t ref_left [(128<<1) + 128]/*={0}*/;
    uint16_t *dst = &src[x0 + (y0*dst_stride)];
    uint16_t *ref1 = ref_above + (1 << log2_pb_height);
    uint16_t *ref2 = ref_left  + (1 << log2_pb_width);
    const struct OVRCNCtx *const rcn_ctx = &ctudec->rcn_ctx;
    const struct DCFunctions *dc = &rcn_ctx->rcn_funcs.dc;
    const struct PlanarFunctions *planar = &rcn_ctx->rcn_funcs.planar;
// FIXED? fill_ref_left_0_mref(src, dst_stride, ref2,
//                      ctudec->progress_map.cols[x0 >> 2],
//         ctudec->progress_map.rows[y0 >> 2],
//         multi_ref_idx, x0, y0,
//         log2_pb_width, log2_pb_height);
fill_ref_left_0_mref(src, dst_stride, ref2,
                   ctudec->rcn_ctx.progress_field.vfield[x0 >> 2],
                    ctudec->rcn_ctx.progress_field.hfield[y0 >> 2],
                    multi_ref_idx, x0, y0,
                    log2_pb_width, log2_pb_height);
                // FIXED? :fill_ref_above_0_mref(src, dst_stride, ref1,
                //                       ctudec->progress_map.rows[y0 >> 2],
                //         ctudec->progress_map.cols[x0 >> 2],
                //         multi_ref_idx, x0 , y0,
                //         log2_pb_width, log2_pb_height);
          fill_ref_above_0_mref(src, dst_stride, ref1,
                  ctudec->rcn_ctx.progress_field.hfield[y0 >> 2],
                  ctudec->rcn_ctx.progress_field.vfield[x0 >> 2],
                  multi_ref_idx, x0 , y0,
                  log2_pb_width, log2_pb_height);

    ref1 += multi_ref_idx;
    ref2 += multi_ref_idx;

    switch (intra_mode) {
    case OVINTRA_PLANAR://PLANAR
    {
        planar->func(ref1, ref2, dst, dst_stride, log2_pb_width, log2_pb_height);
        break;
    }
    case OVINTRA_DC://DC
    {
        dc->func(ref1, ref2, dst, dst_stride, log2_pb_width, log2_pb_height);
        break;
    }
    default://angular
    {
        int pred_mode = derive_wide_angular_mode(log2_pb_width, log2_pb_height,
                                                 intra_mode);

        int is_vertical = pred_mode >= OVINTRA_DIA ? 1 : 0;

        if(is_vertical){
            int mode_idx = pred_mode - OVINTRA_VER ;
            switch (mode_idx) {
            case 0:
                //pure vertical
                vvc_intra_ver(ref1, ref2,
                              dst, dst_stride, log2_pb_width, log2_pb_height);
                break;
            case (16)://Pure diagonal
                intra_angular_vdia_mref(ref1, ref2,
                                        dst, dst_stride,
                                        log2_pb_width, log2_pb_height,
                                        multi_ref_idx);
                break;
            default:
                if (mode_idx < 0){
                    ref1 = ref_above + multi_ref_idx;
                    ref2 = ref_left  + multi_ref_idx;
                    intra_vneg_cubic_mref(ref1, ref2, dst, dst_stride,
                                          log2_pb_width, log2_pb_height,
                                          -mode_idx, multi_ref_idx);
                } else {
                    intra_angular_v_cubic_mref(ref1, dst, dst_stride,
                                          log2_pb_width, log2_pb_height,
                                          mode_idx, multi_ref_idx);

                }
                break;
            }
        } else {
            int mode_idx = -(pred_mode - OVINTRA_HOR);
            switch (mode_idx) {
            case 0:
                //pure horizontal
                vvc_intra_hor(ref1, ref2,
                              dst, dst_stride, log2_pb_width, log2_pb_height);
                break;

            case (16)://Pure diagonal
                intra_angular_hdia_mref(ref1, ref2,
                                        dst, dst_stride,
                                        log2_pb_width, log2_pb_height,
                                        multi_ref_idx);
                break;
            default:
            {
                if (mode_idx < 0){ //TODO move copy inside neg function
                    ref1 = ref_above + multi_ref_idx;
                    ref2 = ref_left  + multi_ref_idx;
                    intra_hneg_cubic_mref(ref1, ref2,
                                          dst, dst_stride,
                                          log2_pb_width, log2_pb_height,
                                          -mode_idx, multi_ref_idx);
                } else {
                    //from 0 to ref_lengths +1, 0 being top_left sample
                    intra_angular_h_cubic_mref(ref2, dst, dst_stride,mode_idx,
                                          log2_pb_width, log2_pb_height,
                                          multi_ref_idx);
                }
            }
                break;
            }
        }
        break;
    }
    }
}
#endif

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

  #if ARCH_X86
    #if SSE_ENABLED
      rcn_init_mc_functions_sse(rcn_func);
      rcn_init_tr_functions_sse(rcn_func);
      rcn_init_dc_planar_functions_sse(rcn_func);
      rcn_init_ict_functions_sse(rcn_func, ict_type);
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
