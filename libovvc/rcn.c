#include <stdint.h>
#include <string.h>
#include "rcn.h"
#include "nvcl_utils.h"
#include "ovutils.h"
#include "rcn_fill_ref.h"
#include "rcn_intra_angular.h"
#include "rcn_intra_dc_planar.h"
#include "data_rcn_angular.h"
#include "ctudec.h"
#include "rcn_intra_mip.h"
#include "rcn_struct.h"
#include "rcn.h"

enum OVIntraMode
{
    OVINTRA_PLANAR = 0,
    OVINTRA_DC = 1,
    OVINTRA_HOR = 18,
    OVINTRA_DIA = 34,
    OVINTRA_VER = 50,
    OVINTRA_VDIA = 66,
    OVINTRA_LM_CHROMA = 67,
    OVINTRA_NUM_INTRA_MODES = 67,
    OVINTRA_MDLM_LEFT = 68,
    OVINTRA_MDLM_TOP  = 69,
    OVINTRA_DM_CHROMA = 70,
    OVINTRA_MIP_MODE = 75,
    OVINTRA_NOT_AVAILABLE=128,
};


/* FIXME 
 * Remove planar and dc cases and use masks on 64 instead of % 
 * Factorize + clean redundancies in derive mpm functions
 */
static const uint8_t mode_shift_tab[6] = {0, 6, 10, 12, 14, 15};

/*
 * Modify angular mode for non square according to width / height
 * ratio.
 * WARNING: do not call if DC or PLANAR, or LM
 *          return value is not the mode to be used 
 *          for derivation but for reconstruction.
 * FIXME clean return unsigned and smaller sizes
 */
static int
derive_wide_angular_mode(int log2_pb_w, int log2_pb_h, int pred_mode)
{
    int mode_shift = mode_shift_tab[OVABS(log2_pb_w - log2_pb_h)];

    if (log2_pb_w > log2_pb_h && pred_mode < 2 + mode_shift) {
        pred_mode += (OVINTRA_VDIA - 1);
    } else if (log2_pb_h > log2_pb_w && pred_mode > OVINTRA_VDIA - mode_shift){
        pred_mode -= (OVINTRA_VDIA - 1);
    }
    return pred_mode;
}


#if 1
static void
vvc_intra_chroma_angular(const uint16_t *const src, uint16_t *const dst,
                         uint16_t *ref_left, uint16_t *ref_above,
                         uint64_t left_col_map, uint64_t top_row_map,
                         int8_t log2_pb_width, int8_t log2_pb_height,
                         int8_t x0, int8_t y0,
                         int8_t intra_mode);

void vvc_intra_pred_chroma(const OVCTUDec *const ctudec,
                           uint16_t *const dst_NUU, uint16_t *const dst_NUk,
                           ptrdiff_t dst_stride,
                           uint8_t intra_mode, int x0, int y0,
                           int log2_pb_w, int log2_pb_h){

    const struct OVRCNCtx *rcn_ctx = &ctudec->rcn_ctx;
    uint16_t *const dst_cb = &rcn_ctx->ctu_buff.cb[(x0) + (y0 * RCN_CTB_STRIDE)];
    uint16_t *const dst_cr = &rcn_ctx->ctu_buff.cr[(x0) + (y0 * RCN_CTB_STRIDE)];
    dst_stride = rcn_ctx->ctu_buff.stride_c;
    /*TODO load ref_sample for cb and cr in same function*/
    uint16_t ref_above[(128<<1) + 128]/*={0}*/;
    uint16_t ref_left [(128<<1) + 128]/*={0}*/;
    uint16_t *ref1 = ref_above;
    uint16_t *ref2 = ref_left;
    uint8_t neighbour= ctudec->ctu_ngh_flags;
    uint8_t got_left_ctu = neighbour & VVC_CTU_LEFT_FLAG;
    uint8_t got_top_ctu  = neighbour & VVC_CTU_UP_FLAG;
    const uint16_t *const src_cb = &ctudec->rcn_ctx.ctu_buff.cb[0];
    const uint16_t *const src_cr = &ctudec->rcn_ctx.ctu_buff.cr[0];

    uint64_t left_col_map = ctudec->rcn_ctx.progress_field_c.vfield[x0 >> 1];
    uint64_t top_row_map  = ctudec->rcn_ctx.progress_field_c.hfield[y0 >> 1];


    switch (intra_mode) {
    case OVINTRA_PLANAR://PLANAR
    {
        fill_ref_left_0_chroma(src_cb, dst_stride, ref_left,
                               left_col_map, top_row_map,
                               x0, y0, log2_pb_w, log2_pb_h);

        fill_ref_above_0_chroma(src_cb, dst_stride, ref_above,
                               top_row_map, left_col_map,
                               x0, y0, log2_pb_w, log2_pb_h);

        if (log2_pb_h > 1 && log2_pb_w > 1)
        // FIXME! vvc_intra_dsp_ctx.planar_pdpc[0](ref1, ref2, dst_cb, dst_stride,
                         // log2_pb_w, log2_pb_h);
                         vvc_intra_planar_pdpc(ref1, ref2, dst_cb, dst_stride, log2_pb_w,
                                          log2_pb_h);
        else

        fill_ref_left_0_chroma(src_cr, dst_stride, ref_left,
                               left_col_map, top_row_map,
                               x0, y0, log2_pb_w, log2_pb_h);

        fill_ref_above_0_chroma(src_cr, dst_stride, ref_above,
                               top_row_map, left_col_map,
                               x0, y0, log2_pb_w, log2_pb_h);

        if (log2_pb_h > 1 && log2_pb_w > 1)
        // FIXED? :vvc_intra_dsp_ctx.planar_pdpc[0](ref1, ref2, dst_cr, dst_stride,
                         // log2_pb_w, log2_pb_h);
                         vvc_intra_planar_pdpc(ref1, ref2, dst_cr, dst_stride,
                                          log2_pb_w, log2_pb_h);
        else
            vvc_intra_planar(ref1, ref2, dst_cr, dst_stride, log2_pb_w,
                             log2_pb_h);
        break;
    }
    case OVINTRA_DC://DC
    {
        fill_ref_left_0_chroma(src_cb, dst_stride, ref_left,
                               left_col_map, top_row_map,
                               x0, y0, log2_pb_w, log2_pb_h);

        fill_ref_above_0_chroma(src_cb, dst_stride, ref_above,
                                top_row_map, left_col_map,
                                x0, y0, log2_pb_w, log2_pb_h);

        /* PDPC disable for 4xX and Xx4 blocks */
        if (log2_pb_h > 1 && log2_pb_w > 1) {
            vvc_intra_dc_pdpc(ref1, ref2, dst_cb, dst_stride, log2_pb_w,
                              log2_pb_h);
        } else {
            vvc_intra_dc(ref1, ref2, dst_cb, dst_stride, log2_pb_w,
                         log2_pb_h);
        }

        fill_ref_left_0_chroma(src_cr, dst_stride, ref_left,
                               left_col_map, top_row_map,
                               x0, y0, log2_pb_w, log2_pb_h);

        fill_ref_above_0_chroma(src_cr, dst_stride, ref_above,
                                top_row_map, left_col_map,
                                x0, y0, log2_pb_w, log2_pb_h);

        /* PDPC disable for 4xX and Xx4 blocks */
        if (log2_pb_h > 1 && log2_pb_w > 1) {
            vvc_intra_dc_pdpc(ref1, ref2, dst_cr, dst_stride, log2_pb_w,
                              log2_pb_h);
        } else {
            vvc_intra_dc(ref1, ref2, dst_cr, dst_stride, log2_pb_w, log2_pb_h);
        }

        break;
    }
    case OVINTRA_LM_CHROMA:
    {
        const uint16_t  *const src_luma = &ctudec->rcn_ctx.ctu_buff.y[(x0<<1)+((y0<<1)*VVC_CTB_STRIDE)];
        // ctudec->cclm_func(src_luma, dst_cb, dst_cr, log2_pb_w, log2_pb_h,
        //                   y0, got_top_ctu || y0, got_left_ctu || x0);

        break;
    }
    case OVINTRA_MDLM_LEFT:
    {
        const uint16_t  *const src_luma = &ctudec->rcn_ctx.ctu_buff.y[(x0<<1)+((y0<<1)*VVC_CTB_STRIDE)];

        // FIXED? : ctudec->left_mdlm(src_luma, dst_cb, dst_cr,
        //                   left_col_map, log2_pb_w, log2_pb_h,
        //                   x0, y0, x0 || got_left_ctu, y0 || got_top_ctu);

        break;
    }
    case OVINTRA_MDLM_TOP:
    {
        const uint16_t  *const src_luma = &ctudec->rcn_ctx.ctu_buff.y[(x0<<1)+((y0<<1)*VVC_CTB_STRIDE)];

        // FIXED? : ctudec->top_mdlm(src_luma, dst_cb, dst_cr,
        //                  top_row_map, log2_pb_w, log2_pb_h,
        //                  x0, y0, x0 || got_left_ctu, y0 || got_top_ctu);

        break;
    }
    default://angular
    {
        vvc_intra_chroma_angular(src_cb, dst_cb, ref_left, ref_above, left_col_map,
                                 top_row_map, log2_pb_w, log2_pb_h,
                                 x0, y0, intra_mode);

        vvc_intra_chroma_angular(src_cr, dst_cr, ref_left, ref_above, left_col_map,
                                 top_row_map, log2_pb_w, log2_pb_h,
                                 x0, y0, intra_mode);
        break;

    }
    }
}
#endif


static void
vvc_intra_chroma_angular(const uint16_t *const src, uint16_t *const dst,
                         uint16_t *ref_left, uint16_t *ref_above,
                         uint64_t left_col_map, uint64_t top_row_map,
                         int8_t log2_pb_width, int8_t log2_pb_height,
                         int8_t x0, int8_t y0,
                         int8_t intra_mode)
{
    int pred_mode = derive_wide_angular_mode(log2_pb_width, log2_pb_height,
                                             intra_mode);

    int dst_stride = VVC_CTB_STRIDE_CHROMA;
    uint16_t *ref1 = ref_above + (1 << log2_pb_height);
    uint16_t *ref2 = ref_left + (1 << log2_pb_width);
    int is_vertical = pred_mode >= OVINTRA_DIA ? 1 : 0;
    //TODO check when this is useful
    //FIXME src and dst are not the same
    fill_ref_left_0_chroma(src, dst_stride, ref2,
                           left_col_map, top_row_map,
                           x0, y0, log2_pb_width, log2_pb_height);

    fill_ref_above_0_chroma(src, dst_stride, ref1,
                            top_row_map, left_col_map,
                            x0, y0, log2_pb_width, log2_pb_height);

    if(is_vertical){
        int mode_idx = pred_mode - (int)OVINTRA_VER;
        switch (mode_idx) {
            case 0:
                //pure vertical
                if (log2_pb_height > 1 && log2_pb_width > 1)
                vvc_intra_ver_pdpc(ref1, ref2, dst, dst_stride,
                                   log2_pb_width, log2_pb_height);
                else
                vvc_intra_ver(ref1, ref2, dst, dst_stride,
                              log2_pb_width, log2_pb_height);

                break;
            case (16)://Pure diagonal
            {
                    int abs_angle = angle_table[mode_idx];
                if (log2_pb_height > 1 && log2_pb_width > 1)
                vvc_intra_angular_vdia(ref1, ref2, dst, dst_stride,
                                       log2_pb_width, log2_pb_height);
                else
                    vvc_intra_angular_v_c(ref1, dst, dst_stride,
                                           log2_pb_width, log2_pb_height,
                                           abs_angle);
                    }
                break;
            default:
                if(mode_idx < 0){
                    int inv_angle = inverse_angle_table[-mode_idx];
                    int abs_angle = angle_table[-mode_idx];
                    int pb_height = 1 << log2_pb_height;
                    int inv_angle_sum    = 256;
                    uint8_t req_frac = !!(abs_angle& 0x1F);

                    for( int k = -1; k >= -pb_height; k-- ){
                        inv_angle_sum += inv_angle;
                        ref1[k] = ref2[OVMIN(inv_angle_sum >> 9, pb_height)];
                    }
                    if (!req_frac){
                        intra_angular_v_nofrac(ref1, dst, dst_stride,
                                log2_pb_width, log2_pb_height,
                                -abs_angle);
                    } else {
                        vvc_intra_angular_v_c(ref1, dst, dst_stride,
                                              log2_pb_width, log2_pb_height,
                                              -abs_angle);
                    }
                } else if (mode_idx < 8&&  OVMIN(2, log2_pb_height - (floor_log2(3*inverse_angle_table[mode_idx] - 2) - 8)) < 0){//FIXME check this
                    int abs_angle = angle_table[mode_idx];
                    uint8_t req_frac = !!(abs_angle& 0x1F);

                    if (!req_frac){
                        intra_angular_v_nofrac(ref1, dst, dst_stride,
                                log2_pb_width, log2_pb_height,
                                abs_angle);
                    } else {
                        vvc_intra_angular_v_c(ref1, dst, dst_stride,
                                              log2_pb_width, log2_pb_height,
                                              abs_angle);
                    }
                } else {//wide angular
                    int abs_angle = angle_table[mode_idx];
                    uint8_t req_frac = !!(abs_angle& 0x1F);
                    if (!req_frac){
                        if (log2_pb_height > 1 && log2_pb_width > 1)
                            intra_angular_v_nofrac_pdpc(ref1, ref2, dst, dst_stride,
                                                        log2_pb_width, log2_pb_height,
                                                        mode_idx);
                        else
                            intra_angular_v_nofrac(ref1, dst, dst_stride,
                                                   log2_pb_width, log2_pb_height,
                                                   abs_angle);
                    } else {
                        if (log2_pb_height > 1 && log2_pb_width > 1)
                        vvc_intra_angular_vpos_wide(ref1, ref2, dst, dst_stride,
                                                    log2_pb_width, log2_pb_height,
                                                    mode_idx);
                        else
                            vvc_intra_angular_v_c(ref1, dst, dst_stride,
                                                   log2_pb_width, log2_pb_height,
                                                   abs_angle);
                    }
                }
                break;
        }
    } else {
        int mode_idx = -(pred_mode - (int)OVINTRA_HOR);
        switch (mode_idx) {
            case 0:
                //pure horizontal
                if (log2_pb_height > 1 && log2_pb_width > 1)
                    vvc_intra_hor_pdpc(ref1, ref2, dst, dst_stride,
                                   log2_pb_width, log2_pb_height);
                else
                    vvc_intra_hor(ref1, ref2, dst, dst_stride,
                                   log2_pb_width, log2_pb_height);
                break;
            case (16)://Pure diagonal
                if (log2_pb_height > 1 && log2_pb_width > 1)
                vvc_intra_angular_hdia(ref1, ref2, dst, dst_stride,
                                       log2_pb_width, log2_pb_height);
                                       else
                            vvc_intra_angular_h_c(ref2, dst, dst_stride,
                                                   log2_pb_width, log2_pb_height,
                                                   32);
                break;
            default:
                {
                    if(mode_idx < 0){
                        int inv_angle = inverse_angle_table[-mode_idx];
                        int abs_angle = angle_table[-mode_idx];
                        int pb_width = 1 << log2_pb_width;
                        int inv_angle_sum    = 256;
                        uint8_t req_frac = !!(abs_angle& 0x1F);

                        for( int k = -1; k >= -pb_width; k-- ){
                            inv_angle_sum += inv_angle;
                            ref2[k] = ref1[OVMIN(inv_angle_sum >> 9, pb_width)];
                        }

                        if (!req_frac){
                            intra_angular_h_nofrac(ref2, dst, dst_stride,
                                                   log2_pb_width, log2_pb_height,
                                                   -abs_angle);
                        } else {
                            vvc_intra_angular_h_c(ref2, dst, dst_stride,
                                                  log2_pb_width, log2_pb_height,
                                                  -abs_angle);
                        }

                    } else if (mode_idx < 8 &&  OVMIN(2, log2_pb_width - (floor_log2(3*inverse_angle_table[mode_idx] - 2) - 8)) < 0 ){//FIXME check this
                        int abs_angle = angle_table[mode_idx];
                        uint8_t req_frac = !!(abs_angle& 0x1F);
                        if (!req_frac){
                            intra_angular_h_nofrac(ref2, dst, dst_stride,
                                                   log2_pb_width, log2_pb_height,
                                                   abs_angle);
                        } else {
                            vvc_intra_angular_h_c(ref2, dst, dst_stride,
                                                  log2_pb_width, log2_pb_height,
                                                  abs_angle);
                        }

                    } else {//wide angular
                        int abs_angle = angle_table[mode_idx];
                        uint8_t req_frac = !!(abs_angle& 0x1F);
                        if (!req_frac){
                            if (log2_pb_height > 1 && log2_pb_width > 1)
                            intra_angular_h_nofrac_pdpc(ref1, ref2, dst, dst_stride,
                                                        log2_pb_width, log2_pb_height,
                                                        mode_idx);
                            else
                            intra_angular_h_nofrac(ref2, dst, dst_stride,
                                                   log2_pb_width, log2_pb_height,
                                                   abs_angle);
                        } else {
                            if (log2_pb_height > 1 && log2_pb_width > 1)
                            vvc_intra_angular_hpos_wide(ref1, ref2, dst, dst_stride,
                                                        log2_pb_width, log2_pb_height,
                                                        mode_idx);
                            else
                            vvc_intra_angular_h_c(ref2, dst, dst_stride,
                                                   log2_pb_width, log2_pb_height,
                                                   abs_angle);

                        }
                    }
                }
                break;
        }
    }
}
#if 0


void
vvc_intra_pred_chroma(const OVCTUDec *const ctudec,
                      uint16_t *const dst_cb, uint16_t *const dst_cr,
                      ptrdiff_t dst_stride,
                      uint8_t intra_mode, int x0, int y0,
                      int log2_pb_width, int log2_pb_height)
{

    /*TODO load ref_sample for cb and cr in same function*/
    uint16_t ref_above[(128<<1) + 128]/*={0}*/;
    uint16_t ref_left [(128<<1) + 128]/*={0}*/;
    uint16_t *ref1 = ref_above;
    uint16_t *ref2 = ref_left;
    uint8_t neighbour= ctudec->ctu_ngh_flags;
    uint8_t got_left_ctu = neighbour & VVC_CTU_LEFT_FLAG;
    uint8_t got_top_ctu  = neighbour & VVC_CTU_UP_FLAG;
    const uint16_t *const src_cb = &ctudec->rcn_ctx.ctu_buff.data_cb[0];
    const uint16_t *const src_cr = &ctudec->rcn_ctx.ctu_buff.data_cr[0];
    // FIXED? : uint64_t left_col_map = ctudec->progress_map_c.cols[x0 >> 1];
    // FIXED? : uint64_t top_row_map  = ctudec->progress_map_c.rows[y0 >> 1];
    uint64_t left_col_map = ctudec->rcn_ctx.progress_field.vfield[x0 >> 1];
    uint64_t top_row_map  = ctudec->rcn_ctx.progress_field.hfield[y0 >> 1];


    switch (intra_mode) {
    case OVINTRA_PLANAR://PLANAR
    {
        fill_ref_left_0_chroma(src_cb, dst_stride, ref_left,
                               left_col_map, top_row_map,
                               x0, y0, log2_pb_width, log2_pb_height);

        fill_ref_above_0_chroma(src_cb, dst_stride, ref_above,
                               top_row_map, left_col_map,
                               x0, y0, log2_pb_width, log2_pb_height);

        if (log2_pb_height > 1 && log2_pb_width > 1)
        // FIXME! vvc_intra_dsp_ctx.planar_pdpc[0](ref1, ref2, dst_cb, dst_stride,
                         // log2_pb_width, log2_pb_height);
                         vvc_intra_planar_pdpc(ref1, ref2, dst_cb, dst_stride, log2_pb_width,
                                          log2_pb_height);
        else

        fill_ref_left_0_chroma(src_cr, dst_stride, ref_left,
                               left_col_map, top_row_map,
                               x0, y0, log2_pb_width, log2_pb_height);

        fill_ref_above_0_chroma(src_cr, dst_stride, ref_above,
                               top_row_map, left_col_map,
                               x0, y0, log2_pb_width, log2_pb_height);

        if (log2_pb_height > 1 && log2_pb_width > 1)
        // FIXED? :vvc_intra_dsp_ctx.planar_pdpc[0](ref1, ref2, dst_cr, dst_stride,
                         // log2_pb_width, log2_pb_height);
                         vvc_intra_planar_pdpc(ref1, ref2, dst_cr, dst_stride,
                                          log2_pb_width, log2_pb_height);
        else
            vvc_intra_planar(ref1, ref2, dst_cr, dst_stride, log2_pb_width,
                             log2_pb_height);
        break;
    }
    case OVINTRA_DC://DC
    {
        fill_ref_left_0_chroma(src_cb, dst_stride, ref_left,
                               left_col_map, top_row_map,
                               x0, y0, log2_pb_width, log2_pb_height);

        fill_ref_above_0_chroma(src_cb, dst_stride, ref_above,
                               top_row_map, left_col_map,
                               x0, y0, log2_pb_width, log2_pb_height);

        if (log2_pb_height > 1 && log2_pb_width > 1)
            // FIXED? :vvc_intra_dsp_ctx.dc_pdpc(ref1, ref2, dst_cb, dst_stride, log2_pb_width,
            //               log2_pb_height);
                          vvc_intra_dc_pdpc(ref1, ref2, dst_cb, dst_stride, log2_pb_width,
                                        log2_pb_height);
        else
            vvc_intra_dc(ref1, ref2, dst_cb, dst_stride, log2_pb_width,
                         log2_pb_height);

        fill_ref_left_0_chroma(src_cr, dst_stride, ref_left,
                               left_col_map, top_row_map,
                               x0, y0, log2_pb_width, log2_pb_height);

        fill_ref_above_0_chroma(src_cr, dst_stride, ref_above,
                               top_row_map, left_col_map,
                               x0, y0, log2_pb_width, log2_pb_height);

        if (log2_pb_height > 1 && log2_pb_width > 1)
            // FIXED? :vvc_intra_dsp_ctx.dc_pdpc(ref1, ref2, dst_cr, dst_stride, log2_pb_width,
                          // log2_pb_height);
                          vvc_intra_dc_pdpc(ref1, ref2, dst_cr, dst_stride, log2_pb_width,
                                        log2_pb_height);
        else
            vvc_intra_dc(ref1, ref2, dst_cr, dst_stride, log2_pb_width,
                         log2_pb_height);


        break;
    }
    case OVINTRA_LM_CHROMA:
    {
        const uint16_t  *const src_luma = &ctudec->rcn_ctx.ctu_buff.data_y[(x0<<1)+((y0<<1)*VVC_CTB_STRIDE)];

        // ctudec->cclm_func(src_luma, dst_cb, dst_cr, log2_pb_width, log2_pb_height,
        //                   y0, got_top_ctu || y0, got_left_ctu || x0);

        break;
    }
    case VVC_MDLM_LEFT:
    {
      const uint16_t  *const src_luma = &ctudec->rcn_ctx.ctu_buff.data_y[(x0<<1)+((y0<<1)*VVC_CTB_STRIDE)];

        // FIXED? : ctudec->left_mdlm(src_luma, dst_cb, dst_cr,
        //                   left_col_map, log2_pb_width, log2_pb_height,
        //                   x0, y0, x0 || got_left_ctu, y0 || got_top_ctu);

        break;
    }
    case VVC_MDLM_TOP:
    {
        const uint16_t  *const src_luma = &ctudec->rcn_ctx.ctu_buff.data_y[(x0<<1)+((y0<<1)*VVC_CTB_STRIDE)];

        // FIXED? : ctudec->top_mdlm(src_luma, dst_cb, dst_cr,
        //                  top_row_map, log2_pb_width, log2_pb_height,
        //                  x0, y0, x0 || got_left_ctu, y0 || got_top_ctu);

        break;
    }
    default://angular
    {
        vvc_intra_chroma_angular(src_cb, dst_cb, ref_left, ref_above, left_col_map,
                                 top_row_map, log2_pb_width, log2_pb_height,
                                 x0, y0, intra_mode);

        vvc_intra_chroma_angular(src_cr, dst_cr, ref_left, ref_above, left_col_map,
                                 top_row_map, log2_pb_width, log2_pb_height,
                                 x0, y0, intra_mode);
        break;

    }
    }
}
#endif

void
vvc_intra_pred(const OVCTUDec *const ctudec,
               uint16_t *const src,
               ptrdiff_t dst_stride,
               uint8_t intra_mode, int x0, int y0,
               int log2_pb_width, int log2_pb_height)
{

    uint16_t ref_above[(128<<1) + 128]/*={0}*/;
    uint16_t ref_left [(128<<1) + 128]/*={0}*/;
    uint16_t ref_above_filtered[(128<<1) + 128]/*={0}*/;
    uint16_t ref_left_filtered [(128<<1) + 128]/*={0}*/;
    uint16_t *dst = &src[x0 + (y0*dst_stride)];
    uint16_t *ref1 = ref_above + (1 << log2_pb_height);
    uint16_t *ref2 = ref_left + (1 << log2_pb_width);
    // FIXED?: fill_ref_left_0(src,dst_stride,ref2,
    //                 ctudec->progress_map.cols[x0 >> 2],
    //         ctudec->progress_map.rows[y0 >> 2],
    //         x0, y0, log2_pb_width, log2_pb_height,0);

            fill_ref_left_0(src,dst_stride,ref2,
                            ctudec->rcn_ctx.progress_field.vfield[x0 >> 2],
                    ctudec->rcn_ctx.progress_field.hfield[y0 >> 2],
                    x0, y0, log2_pb_width, log2_pb_height,0);

    // FIXED?: fill_ref_above_0(src, dst_stride, ref1,
            //          ctudec->progress_map.rows[y0 >> 2],
            // ctudec->progress_map.cols[x0 >> 2],
            // x0, y0, log2_pb_width, log2_pb_height,0);

            fill_ref_above_0(src, dst_stride, ref1,
                             ctudec->rcn_ctx.progress_field.hfield[y0 >> 2],
                    ctudec->rcn_ctx.progress_field.vfield[x0 >> 2],
                    x0, y0, log2_pb_width, log2_pb_height,0);
    switch (intra_mode) {
    case OVINTRA_PLANAR://PLANAR
    {
        if((log2_pb_height + log2_pb_width) > 5){
            filter_ref_samples(ref1, ref_above_filtered, ref2,
                               (1 << log2_pb_width)+4);
            filter_ref_samples(ref2, ref_left_filtered, ref1,
                               (1 << log2_pb_height)+4);
            ref1 = ref_above_filtered;
            ref2 = ref_left_filtered;
        }
        // FIXED? :vvc_intra_dsp_ctx.planar_pdpc[log2_pb_width > 5 || log2_pb_height > 5](ref1, ref2, dst, dst_stride,
                              // log2_pb_width, log2_pb_height);
        vvc_intra_planar_pdpc(ref1, ref2, dst, dst_stride,
                              log2_pb_width, log2_pb_height);
        break;
    }
    case OVINTRA_DC://DC
    {
      // FIXED? :vvc_intra_dsp_ctx.dc_pdpc(ref1, ref2, dst, dst_stride, log2_pb_width, log2_pb_height);
      vvc_intra_dc_pdpc(ref1, ref2, dst, dst_stride, log2_pb_width, log2_pb_height);

        break;
    }
    default://angular
    {
        int pred_mode = derive_wide_angular_mode(log2_pb_width, log2_pb_height,
                                                 intra_mode);

        int is_vertical = pred_mode >= OVINTRA_DIA ? 1 : 0;

        if(is_vertical){
            int mode_idx = pred_mode - (int)OVINTRA_VER;

            //FIXME recheck filter derivation
            int use_gauss_filter = log2_pb_width + log2_pb_height > 5 &&
                        (OVABS(mode_idx) > intra_filter[((log2_pb_width +
                                                          log2_pb_height) >> 1)]);

            switch (mode_idx) {
            case 0: //Pure vertical
                vvc_intra_ver_pdpc(ref1, ref2, dst, dst_stride, log2_pb_width,
                                   log2_pb_height);
                break;
            case (16)://Pure diagonal
                if(use_gauss_filter){
                    int top_ref_length  = 1 << (log2_pb_width  + 1);
                    int left_ref_length = 1 << (log2_pb_height + 1);
                    filter_ref_samples(ref1, ref_above_filtered, ref2,
                                       top_ref_length);
                    filter_ref_samples(ref2, ref_left_filtered, ref1,
                                       left_ref_length);
                    ref1 = ref_above_filtered;
                    ref2 = ref_left_filtered;
                }
                vvc_intra_angular_vdia(ref1, ref2, dst, dst_stride,
                                       log2_pb_width, log2_pb_height);
                break;
            default:{
                if(mode_idx < 0){
                    int top_ref_length  = 1 << (log2_pb_width  + 1);
                    int left_ref_length = 1 << (log2_pb_height + 1);
                    int pu_width = 1 << log2_pb_width;
                    int pu_height = 1 << log2_pb_height;
                    int abs_angle_val = angle_table[-mode_idx];
                    uint8_t req_frac = !!(abs_angle_val & 0x1F);
                    if (!req_frac){
                        if (use_gauss_filter){
                            filter_ref_samples(ref1, ref_above_filtered + pu_height,
                                               ref2, top_ref_length);
                            filter_ref_samples(ref2, ref_left_filtered + pu_width,
                                               ref1, left_ref_length);
                            ref1 = ref_above_filtered + pu_height;
                            ref2 = ref_left_filtered + pu_width;
                        }

                        int inv_angle = inverse_angle_table[-mode_idx];
                        int inv_angle_sum    = 256;
                        for ( int k = -1; k >= -pu_height; k-- ){
                            inv_angle_sum += inv_angle;
                            ref1[k] = ref2[OVMIN(inv_angle_sum >> 9,pu_height)];
                        }

                        intra_angular_v_nofrac(ref1, dst, dst_stride,
                                               log2_pb_width, log2_pb_height,
                                               -abs_angle_val);
                    } else {
                        int inv_angle = inverse_angle_table[-mode_idx];
                        int inv_angle_sum    = 256;
                        for ( int k = -1; k >= -pu_height; k-- ){
                            inv_angle_sum += inv_angle;
                            ref1[k] = ref2[OVMIN(inv_angle_sum >> 9,pu_height)];
                        }

                        if (use_gauss_filter){
                            intra_angular_v_gauss(ref1, dst, dst_stride,
                                                  log2_pb_width, log2_pb_height,
                                                  -abs_angle_val);
                        } else {
                            intra_angular_v_cubic(ref1, dst, dst_stride,
                                                  log2_pb_width, log2_pb_height,
                                                  -abs_angle_val);
                        }
                    }
                } else if (OVMIN(2, log2_pb_height - (floor_log2(3*inverse_angle_table[mode_idx] - 2) - 8)) < 0 ){
                    //FIXME check this
                    int abs_angle_val = angle_table[mode_idx];
                    uint8_t req_frac = !!(abs_angle_val & 0x1F);
                    if (!req_frac){
                        if (use_gauss_filter){
                            int top_ref_length  = 1 << (log2_pb_width  + 1);
                            filter_ref_samples(ref1, ref_above_filtered, ref2,
                                               top_ref_length);
                            ref1 = ref_above_filtered;
                        }
                        intra_angular_v_nofrac(ref1, dst, dst_stride,
                                               log2_pb_width, log2_pb_height,
                                               abs_angle_val);
                    } else {
                        if (use_gauss_filter){
                            intra_angular_v_gauss(ref1, dst, dst_stride,
                                                  log2_pb_width, log2_pb_height,
                                                  abs_angle_val);
                        } else {
                            intra_angular_v_cubic(ref1, dst, dst_stride,
                                                  log2_pb_width, log2_pb_height,
                                                  abs_angle_val);
                        }
                    }
                } else {
                    uint8_t req_frac = !!(angle_table[mode_idx] & 0x1F);
                    if (!req_frac){
                        if (use_gauss_filter){
                            int top_ref_length  = 1 << (log2_pb_width  + 1);
                            int left_ref_length = 1 << (log2_pb_height + 1);
                            filter_ref_samples(ref1, ref_above_filtered, ref2,
                                               top_ref_length);
                            filter_ref_samples(ref2, ref_left_filtered, ref1,
                                               left_ref_length);
                            ref1 = ref_above_filtered;
                            ref2 = ref_left_filtered;
                        }
                        intra_angular_v_nofrac_pdpc(ref1, ref2, dst, dst_stride,
                                                    log2_pb_width, log2_pb_height,
                                                    mode_idx);
                    } else {
                        if (use_gauss_filter){
                        intra_angular_v_gauss_pdpc(ref1, ref2, dst, dst_stride,
                                             log2_pb_width, log2_pb_height,
                                             mode_idx);
                        } else {
                            intra_angular_v_cubic_pdpc(ref1, ref2, dst, dst_stride,
                                                  log2_pb_width, log2_pb_height,
                                                  mode_idx);
                        }
                    }
                }
                break;
            }
            }
        } else {
            int mode_idx = -(pred_mode - (int)OVINTRA_HOR);

            //FIXME recheck filter derivation
            int use_gauss_filter = log2_pb_width + log2_pb_height > 5 &&
                        (OVABS(mode_idx) > intra_filter[((log2_pb_width +
                                                          log2_pb_height) >> 1)]) ? 1: 0;

            switch (mode_idx) {
            case 0: //Pure horizontal
                vvc_intra_hor_pdpc(ref1, ref2, dst, dst_stride,
                                   log2_pb_width, log2_pb_height);
                break;

            case (16)://Pure diagonal
                if (use_gauss_filter){
                    int top_ref_length  = 1 << (log2_pb_width  + 1);
                    int left_ref_length = 1 << (log2_pb_height + 1);
                    filter_ref_samples(ref1, ref_above_filtered, ref2,
                                       top_ref_length);
                    filter_ref_samples(ref2, ref_left_filtered, ref1,
                                       left_ref_length);
                    ref1 = ref_above_filtered;
                    ref2 = ref_left_filtered;
                }
                vvc_intra_angular_hdia(ref1, ref2, dst, dst_stride,
                                       log2_pb_width, log2_pb_height);
                break;
            default:
            {
                if (mode_idx < 0){
                    int top_ref_length  = 1 << (log2_pb_width  + 1);
                    int left_ref_length = 1 << (log2_pb_height + 1);
                    int pu_width  = 1 << log2_pb_width;
                    int pu_height = 1 << log2_pb_height;
                    int abs_angle_val = angle_table[-mode_idx];

                    uint8_t req_frac = !!(abs_angle_val & 0x1F);
                    if (!req_frac){
                        if (use_gauss_filter){
                            filter_ref_samples(ref1, ref_above_filtered + pu_height,
                                               ref2, top_ref_length);
                            filter_ref_samples(ref2, ref_left_filtered + pu_width,
                                               ref1, left_ref_length);
                            ref1 = ref_above_filtered + pu_height;
                            ref2 = ref_left_filtered + pu_width;
                        }

                        int inv_angle = inverse_angle_table[-mode_idx];
                        int inv_angle_sum    = 256;
                        for ( int k = -1; k >= -pu_width; k-- ){
                            inv_angle_sum += inv_angle;
                            ref2[k] = ref1[OVMIN(inv_angle_sum >> 9,pu_width)];
                        }
                        intra_angular_h_nofrac(ref2, dst, dst_stride,
                                                log2_pb_width, log2_pb_height,
                                                -abs_angle_val);
                    } else {
                        int inv_angle = inverse_angle_table[-mode_idx];
                        int inv_angle_sum    = 256;
                        for( int k = -1; k >= -pu_width; k-- ){
                            inv_angle_sum += inv_angle;
                            ref2[k] = ref1[OVMIN(inv_angle_sum >> 9,pu_width)];
                        }
                        if (use_gauss_filter){
                            intra_angular_h_gauss(ref2, dst, dst_stride,
                                                  log2_pb_width, log2_pb_height,
                                                  -abs_angle_val);
                        } else {
                            intra_angular_h_cubic(ref2, dst, dst_stride,
                                                  log2_pb_width, log2_pb_height,
                                                  -abs_angle_val);
                        }
                    }

                } else if (OVMIN(2, log2_pb_width - (floor_log2(3*inverse_angle_table[mode_idx] - 2) - 8)) < 0 ){//FIXME check this
                    //from 0 to ref_lengths +1, 0 being top_left sample
                    int abs_angle_val = angle_table[mode_idx];
                    uint8_t req_frac = !!(abs_angle_val & 0x1F);
                    if (!req_frac){
                        if (use_gauss_filter){
                            int left_ref_length = 1 << (log2_pb_height + 1);
                            filter_ref_samples(ref2, ref_left_filtered,
                                               ref1, left_ref_length);
                            ref2 = ref_left_filtered;
                        }
                        intra_angular_h_nofrac(ref2, dst, dst_stride,
                                                log2_pb_width, log2_pb_height,
                                                abs_angle_val);
                    } else {
                        if (use_gauss_filter){
                        intra_angular_h_gauss(ref2, dst, dst_stride,
                                         log2_pb_width, log2_pb_height,
                                         abs_angle_val);
                        } else {
                            intra_angular_h_cubic(ref2, dst, dst_stride,
                                                  log2_pb_width, log2_pb_height,
                                                  abs_angle_val);
                        }
                    }
                } else {
                    uint8_t req_frac = !!(angle_table[mode_idx] & 0x1F);
                    if (!req_frac){
                        if (use_gauss_filter){
                            int top_ref_length  = 1 << (log2_pb_width  + 1);
                            int left_ref_length = 1 << (log2_pb_height + 1);
                            filter_ref_samples(ref1, ref_above_filtered,
                                    ref2, top_ref_length);
                            filter_ref_samples(ref2, ref_left_filtered,
                                    ref1, left_ref_length);
                            ref1 = ref_above_filtered;
                            ref2 = ref_left_filtered;
                        }
                        intra_angular_h_nofrac_pdpc(ref1, ref2, dst, dst_stride,
                                                    log2_pb_width, log2_pb_height,
                                                    mode_idx);
                    } else {
                        if (use_gauss_filter){
                        intra_angular_h_gauss_pdpc(ref1, ref2, dst, dst_stride,
                                             log2_pb_width, log2_pb_height,
                                             mode_idx);
                        } else {
                            intra_angular_h_cubic_pdpc(ref1, ref2, dst, dst_stride,
                                                  log2_pb_width, log2_pb_height,
                                                  mode_idx);
                        }
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

#if 0
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
    case VVC_PLANAR://PLANAR
    {
        if (log2_pb_height > 1)
        // FIXED? :vvc_intra_dsp_ctx.planar_pdpc[log2_pb_width > 5 || log2_pb_height > 5](ref1, ref2, dst, dst_stride, log2_pb_width,
        //                  log2_pb_height);
                         vvc_intra_planar_pdpc(ref1, ref2, dst, dst_stride, log2_pb_width,
                                          log2_pb_height);
        else
            vvc_intra_planar(ref1, ref2, dst, dst_stride, log2_pb_width,
                             log2_pb_height);
        break;
    }
    case VVC_DC://DC
    {
        if (log2_pb_height > 1)
        // FIXED? :vvc_intra_dsp_ctx.dc_pdpc(ref1, ref2, dst, dst_stride, log2_pb_width,
        //              log2_pb_height);
                     vvc_intra_dc_pdpc(ref1, ref2, dst, dst_stride, log2_pb_width,
                                  log2_pb_height);
        else
            vvc_intra_dc(ref1, ref2, dst, dst_stride, log2_pb_width,
                         log2_pb_height);
        break;
    }
    default://angular
    {
        int pred_mode = derive_wide_angular_mode(log2_cu_width, log2_cu_height,
                                                 intra_mode);

        int is_vertical = pred_mode >= VVC_DIA ? 1 : 0;

        if(is_vertical){
            int mode_idx = pred_mode - VVC_VER;
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

                } else if (mode_idx < 8 &&  OVMIN(2, log2_pb_height - (floor_log2(3*inverse_angle_table[mode_idx] - 2) - 8)) < 0 || log2_pb_height < 2){//FIXME check this
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
            int mode_idx = -(pred_mode - VVC_HOR);
            switch (mode_idx) {
            case 0:
                //pure horizontal
                if (log2_pb_height > 1){
                    vvc_intra_hor_pdpc(ref1, ref2, dst, dst_stride,
                                       log2_pb_width, log2_pb_height);
                } else {
                    vvc_intra_hor(ref1, ref2, dst, dst_stride,
                                  log2_pb_width, log2_pb_height);
                }
                break;
            case (16)://Pure diagonal
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
                } else if (mode_idx < 8 &&  OVMIN(2, log2_pb_width - (floor_log2(3*inverse_angle_table[mode_idx] - 2) - 8)) < 0 || log2_pb_height < 2){//FIXME check this
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
        vvc_intra_planar(ref1, ref2, dst, dst_stride,
                         log2_pb_width, log2_pb_height);
        break;
    }
    case OVINTRA_DC://DC
    {
        vvc_intra_dc(ref1, ref2, dst, dst_stride, log2_pb_width, log2_pb_height);

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



void
vvc_intra_pred_mip(const OVCTUDec *const ctudec,
                   uint16_t *const dst,
                   int x0, int y0, int log2_pu_w, int log2_pu_h,
                   uint8_t mip_mode)
{
    /* FIXME used ?*/
    const uint16_t *src = &ctudec->rcn_ctx.ctu_buff.data_y[0];

    int32_t bndy_line[8];//buffer used to store averaged boundaries use int
    int16_t mip_pred[64];//buffer used to store reduced matrix vector results

    /* FIXME determine max size of those buffers */
    uint16_t ref_abv[(128<<1) + 128];
    uint16_t ref_lft[(128<<1) + 128];

    int dst_stride = RCN_CTB_STRIDE;
// FIXED? :fill_ref_left_0(src, dst_stride, ref_lft,
//                 ctudec->progress_map.cols[x0 >> 2],
//                 ctudec->progress_map.rows[y0 >> 2],
//                 x0, y0, log2_pu_w, log2_pu_h, 0);
// FIXED? :fill_ref_above_0(src, dst_stride, ref_abv,
//                  ctudec->progress_map.rows[y0 >> 2],
//                  ctudec->progress_map.cols[x0 >> 2],
//                  x0, y0, log2_pu_w, log2_pu_h, 0);
   fill_ref_left_0(src, dst_stride, ref_lft,
                   ctudec->rcn_ctx.progress_field.vfield[x0 >> 2],
                   ctudec->rcn_ctx.progress_field.hfield[y0 >> 2],
                   x0, y0, log2_pu_w, log2_pu_h, 0);
   fill_ref_above_0(src, dst_stride, ref_abv,
                    ctudec->rcn_ctx.progress_field.hfield[y0 >> 2],
                    ctudec->rcn_ctx.progress_field.vfield[x0 >> 2],
                    x0, y0, log2_pu_w, log2_pu_h, 0);

    //compute reduced boundaries
    uint8_t log2_bndy = 1 << ((log2_pu_w > 2) || (log2_pu_h > 2));
    uint8_t log2_bnd_x = log2_pu_w - log2_bndy;
    uint8_t log2_bnd_y = log2_pu_h - log2_bndy;
    int i, j;

    int rnd = (1 << log2_bnd_x) >> 1;
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
            uint16_t *_dst = dst + ((1 << log2_scale_y) - 1) * VVC_CTB_STRIDE;
            up_sample(_dst, mip_pred, ref_lft, log2_red_w, log2_red_h,
                       1, (1 << log2_red_w),
                       1, (1 << log2_scale_y) * VVC_CTB_STRIDE,
                       (1 << log2_scale_y), log2_scale_x);
            src        = _dst;
            src_step   = (1 << log2_scale_y) * VVC_CTB_STRIDE;
            src_stride = 1;
        } else { //TODO use mip_pred directly in next ste
            src        = mip_pred;
            src_step   = (1 << log2_pu_w);
            src_stride = 1;
        }
        up_sample(dst, src, ref_abv, log2_red_h, log2_pu_w,
                   src_step, src_stride,
                   VVC_CTB_STRIDE, 1,
                   1, log2_scale_y);

    } else {//write to dst
        for (i = 0; i < (1 << log2_red_h); ++i) {
            for (j = 0; j < (1 << log2_red_w); ++j) {
                dst [j + i * VVC_CTB_STRIDE] = mip_pred[(i << log2_red_w) + j];
            }
        }
    }
}

void
vvc_intra_pred_mip_tr(const OVCTUDec *const ctudec,
                      uint16_t *const dst,
                      int x0, int y0, int log2_pu_w, int log2_pu_h,
                      uint8_t mip_mode)
{
    uint8_t log2_bndy = 1 << ((log2_pu_w > 2) || (log2_pu_h > 2));

    uint8_t log2_red_w;
    uint8_t log2_red_h;

    int rnd;
    int i, j;
    const uint16_t *src = &ctudec->rcn_ctx.ctu_buff.data_y[0];

    int32_t bndy_line[8]; //buffer used to store averaged boundaries use int
    int16_t mip_pred[64];//buffer used to store reduced matrix vector results

    uint16_t ref_abv[(128<<1) + 128];
    uint16_t ref_lft[(128<<1) + 128];

    int dst_stride = VVC_CTB_STRIDE;

    // FIXED? :fill_ref_left_0(src, dst_stride, ref_lft,
    //                 ctudec->progress_map.cols[x0 >> 2],
    //                 ctudec->progress_map.rows[y0 >> 2],
    //                 x0, y0, log2_pu_w, log2_pu_h, 0);
    //
    // FIXED? :fill_ref_above_0(src, dst_stride, ref_abv,
    //                  ctudec->progress_map.rows[y0 >> 2],
    //                  ctudec->progress_map.cols[x0 >> 2],
    //                  x0, y0, log2_pu_w, log2_pu_h, 0);
   fill_ref_left_0(src, dst_stride, ref_lft,
                   ctudec->rcn_ctx.progress_field.vfield[x0 >> 2],
                   ctudec->rcn_ctx.progress_field.hfield[y0 >> 2],
                   x0, y0, log2_pu_w, log2_pu_h, 0);

   fill_ref_above_0(src, dst_stride, ref_abv,
                    ctudec->rcn_ctx.progress_field.hfield[y0 >> 2],
                    ctudec->rcn_ctx.progress_field.vfield[x0 >> 2],
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
                uint16_t *_dst = dst + ((1 << log2_scale_y) - 1) * VVC_CTB_STRIDE;
                up_sample(_dst, mip_pred, ref_lft,
                           log2_red_w, log2_red_h,
                           (1 << log2_red_h), 1,
                           1, (1 << log2_scale_y) * VVC_CTB_STRIDE,
                           (1 << log2_scale_y), log2_scale_x);
                src        = _dst;
                src_step   = (1 << log2_scale_y) * VVC_CTB_STRIDE;
                src_stride = 1;
            } else { //TODO use mip_pred directly in next ste
                src        = mip_pred;
                src_step   = 1;
                src_stride = 1 << log2_red_h;
            }
            up_sample(dst, src, ref_abv, log2_red_h, log2_pu_w,
                       src_step, src_stride,
                       VVC_CTB_STRIDE, 1,
                       1, log2_scale_y);

    } else {
        for (i = 0; i < (1 << log2_red_h); ++i) {
            for (j = 0; j < (1 << log2_red_w); ++j) {
                dst [j + i * VVC_CTB_STRIDE] = mip_pred[(j << log2_red_h) + i];
            }
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
    int shift_v = 6 + 1;
    int shift_h = (6 + 15 - 1) - 10;
    int16_t tmp[64*64];
    int tb_w = 1 << log2_tb_w;
    int tb_h = 1 << log2_tb_h;
    struct OVRCNCtx *const rcn = &ctudec->rcn_ctx;
    /* TODO switch */

    memset(tmp, 0, sizeof(int16_t) << (log2_tb_w + log2_tb_h));

    if (lfnst_flag) {
        lim_cg_w = 8;
        /* FIXME separate lfnst mode derivation from lfnst reconstruction */
#if 0
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
        tr_templates[tr_v_idx][log2_tb_w].transform(src, tmp, tb_w, tb_w, tb_h, shift_v);
        tr_templates[tr_h_idx][log2_tb_w].transform(tmp, dst, tb_h, tb_h, tb_w, shift_h);

    } else if (!cu_mts_flag) {

#if 0
        if (is_dc) {

            vvc_dsp_context.vvc_inverse_dc(dst, log2_tb_w, log2_tb_h,
                                           src[0]);

        } else {
#endif
            int nb_row = tb_w;//OVMIN(lim_cg_w, 1 << log2_tb_w);
            int nb_col = tb_h;//OVMIN(lim_cg_w, 1 << log2_tb_h);

            tr_templates[DCT_II][log2_tb_h].transform(src, tmp, tb_w, nb_row, nb_col, shift_v);
            tr_templates[DCT_II][log2_tb_w].transform(tmp, dst, tb_h, tb_h, nb_row, shift_h);
#if 0
        }
#endif
    } else {
        enum DCTType tr_h_idx = cu_mts_idx  & 1;
        enum DCTType tr_v_idx = cu_mts_idx >> 1;

        tr_templates[tr_v_idx][log2_tb_h].transform(src, tmp, tb_w, tb_w, tb_h, shift_v);
        tr_templates[tr_h_idx][log2_tb_w].transform(tmp, src, tb_h, tb_h, tb_w, shift_h);

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
    const int shift_v = 6 + 1;
    const int shift_h = (6 + 15 - 1) - 10;
    int16_t tmp[32*32];
    int tb_w = 1 << log2_tb_w;
    int tb_h = 1 << log2_tb_h;

    memset(tmp, 0, sizeof(int16_t) << (log2_tb_w + log2_tb_h));

    if (lfnst_flag) {
        /* Update lim_cg_w since lfnst part of coeff are now non zero */
        lim_cg_w = 8;
        /* FIXME separate lfnst mode derivation from lfnst reconstruction */
#if 0
        process_lfnst(ctudec, src, lfnst_sb, log2_tb_w, log2_tb_h,
                      x0, y0, lfnst_idx);
#endif
    }

#if 0
    if (is_dc && !lfnst_flag) {

        vvc_dsp_context.vvc_inverse_dc(dst, log2_tb_w, log2_tb_h,
                                       src[0]);

    } else {
#endif
        int nb_row = OVMIN(lim_cg_w, 1 << log2_tb_w);
        int nb_col = OVMIN(lim_cg_w, 1 << log2_tb_h);
        /*FIXME might be transform SKIP */

        tr_templates[DCT_II][log2_tb_h].transform(src, tmp, tb_w, tb_w, tb_h, shift_v);
        tr_templates[DCT_II][log2_tb_w].transform(tmp, src, tb_h, tb_h, tb_w, shift_v);
#if 0
    }
#endif
}
