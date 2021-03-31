/* Derivations functions for intra Most Probable Modes (MPM) and conversion
 * to "usable" mode for reconstruction 
 */

#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include "ovutils.h"
#include "ctudec.h"
#include "dec_structures.h"
#include "drv.h"
#include "vcl.h"
#include "drv_utils.h"
#include "rcn.h"
#include "rcn_fill_ref.h"
#include "rcn_intra_angular.h"
#include "rcn_intra_dc_planar.h"
#include "data_rcn_angular.h"

/* Modify angular mode for non square according to width / height
 * ratio.
 * WARNING: do not call if DC or PLANAR, or LM
 *          return value is not the mode to be used
 *          for derivation but for reconstruction.
 * FIXME clean return unsigned and smaller sizes
 * FIXME remove the + 2 if specialized angular modes
 */
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


/* FIXME
 * try to remove + 2  and use mask instead of %
 * Factorize + clean redundancies in derive mpm functions
 */
static uint8_t
vvc_derive_mpm_list(uint8_t mpm_idx, uint8_t lft_mode, uint8_t abv_mode)
{
    const int offset = 67 - 6;
    const int mod    = offset + 3;
    int8_t mpm_list[6];

    /* FIXME Is there a diff with unsorted version ?*/
    mpm_list[0] = OVINTRA_PLANAR;
    mpm_list[1] = OVINTRA_DC;
    mpm_list[2] = OVINTRA_VER;
    mpm_list[3] = OVINTRA_HOR;
    mpm_list[4] = OVINTRA_VER - 4;
    mpm_list[5] = OVINTRA_VER + 4;

    if (lft_mode == abv_mode){
        if (lft_mode > OVINTRA_DC){
            mpm_list[0] = OVINTRA_PLANAR;
            mpm_list[1] = lft_mode;
            mpm_list[2] = ((lft_mode + offset    ) % mod) + 2;
            mpm_list[3] = ((lft_mode          - 1) % mod) + 2;
            mpm_list[4] = ((lft_mode + offset - 1) % mod) + 2;
            mpm_list[5] = ((lft_mode             ) % mod) + 2;
        }
    } else { //L!=A
        if ((lft_mode > OVINTRA_DC) && (abv_mode > OVINTRA_DC)){
            mpm_list[0] = OVINTRA_PLANAR;
            mpm_list[1] = lft_mode;
            mpm_list[2] = abv_mode;
            if (lft_mode > abv_mode){
                if (lft_mode - abv_mode == 1){
                    mpm_list[3] = ((abv_mode + offset    ) % mod) + 2;
                    mpm_list[4] = ((lft_mode           - 1) % mod) + 2;
                    mpm_list[5] = ((abv_mode + offset - 1) % mod) + 2;
                } else if (lft_mode - abv_mode >= 62){
                    mpm_list[3] = ((abv_mode          - 1) % mod) + 2;
                    mpm_list[4] = ((lft_mode + offset     ) % mod) + 2;
                    mpm_list[5] = ( abv_mode               % mod) + 2;
                } else if (lft_mode - abv_mode == 2) {
                    mpm_list[3] = ((abv_mode          - 1) % mod) + 2;
                    mpm_list[4] = ((abv_mode + offset    ) % mod) + 2;
                    mpm_list[5] = ((lft_mode           - 1) % mod) + 2;
                } else {
                    mpm_list[3] = ((abv_mode + offset    ) % mod) + 2;
                    mpm_list[4] = ((abv_mode          - 1) % mod) + 2;
                    mpm_list[5] = ((lft_mode + offset     ) % mod) + 2;
                }
            } else {
                if (abv_mode - lft_mode == 1){
                    mpm_list[3] = ((lft_mode + offset    ) % mod) + 2;
                    mpm_list[4] = ((abv_mode         - 1) % mod) + 2;
                    mpm_list[5] = ((lft_mode + offset - 1) % mod) + 2;
                } else if (abv_mode - lft_mode >= 62){
                    mpm_list[3] = ((lft_mode          - 1) % mod) + 2;
                    mpm_list[4] = ((abv_mode + offset   ) % mod) + 2;
                    mpm_list[5] = ( lft_mode               % mod) + 2;
                } else if (abv_mode - lft_mode == 2) {
                    mpm_list[3] = ((lft_mode          - 1) % mod) + 2;
                    mpm_list[4] = ((lft_mode + offset    ) % mod) + 2;
                    mpm_list[5] = ((abv_mode         - 1) % mod) + 2;
                } else {
                    mpm_list[3] = ((lft_mode + offset     ) % mod) + 2;
                    mpm_list[4] = ((lft_mode           - 1) % mod) + 2;
                    mpm_list[5] = ((abv_mode + offset    ) % mod) + 2;
                }
            }
        } else if (lft_mode + abv_mode >= 2){
            if (lft_mode > abv_mode){
                mpm_list[0] = OVINTRA_PLANAR;
                mpm_list[1] = lft_mode;
                mpm_list[2] = ((lft_mode + offset    ) % mod) + 2;
                mpm_list[3] = ((lft_mode          - 1) % mod) + 2;
                mpm_list[4] = ((lft_mode + offset - 1) % mod) + 2;
                mpm_list[5] = ( lft_mode               % mod) + 2;
            } else {
                mpm_list[0] = OVINTRA_PLANAR;
                mpm_list[1] = abv_mode;
                mpm_list[2] = ((abv_mode + offset    ) % mod) + 2;
                mpm_list[3] = ((abv_mode          - 1) % mod) + 2;
                mpm_list[4] = ((abv_mode + offset - 1) % mod) + 2;
                mpm_list[5] = ( abv_mode               % mod) + 2;
            }
        }
    }

    return mpm_list[mpm_idx];
}

/* FIXME
 * write our own sorting function for faster sorting
 */
static int
cmpfunc (const void * a, const void * b)
{
   return ( *(uint8_t*)a - *(uint8_t*)b );
}

static uint8_t
vvc_derive_mpm_list_sorted(uint8_t lft_mode, uint8_t abv_mode, uint8_t intra_luma_comp)
{
    const int offset = 67 - 6;
    const int mod    = offset + 3;
    int8_t mpm_list[6];

    mpm_list[0] = OVINTRA_PLANAR;
    mpm_list[1] = OVINTRA_DC;
    mpm_list[2] = OVINTRA_VER;
    mpm_list[3] = OVINTRA_HOR;
    mpm_list[4] = OVINTRA_VER - 4;
    mpm_list[5] = OVINTRA_VER + 4;

    if (lft_mode == abv_mode){
        if (lft_mode > OVINTRA_DC){
            mpm_list[0] = OVINTRA_PLANAR;
            mpm_list[1] = lft_mode;
            mpm_list[2] = ((lft_mode + offset    ) % mod) + 2;
            mpm_list[3] = ((lft_mode          - 1) % mod) + 2;
            mpm_list[4] = ((lft_mode + offset - 1) % mod) + 2;
            mpm_list[5] = ((lft_mode             ) % mod) + 2;
        }
    } else { //L!=A
        if ((lft_mode > OVINTRA_DC) && (abv_mode > OVINTRA_DC)){
            mpm_list[0] = OVINTRA_PLANAR;
            mpm_list[1] = lft_mode;
            mpm_list[2] = abv_mode;
            if (lft_mode > abv_mode){
                if (lft_mode - abv_mode == 1){
                    mpm_list[3] = ((abv_mode + offset    ) % mod) + 2;
                    mpm_list[4] = ((lft_mode          - 1) % mod) + 2;
                    mpm_list[5] = ((abv_mode + offset - 1) % mod) + 2;
                } else if (lft_mode - abv_mode >= 62){
                    mpm_list[3] = ((abv_mode          - 1) % mod) + 2;
                    mpm_list[4] = ((lft_mode + offset    ) % mod) + 2;
                    mpm_list[5] = ( abv_mode               % mod) + 2;
                } else if (lft_mode - abv_mode == 2) {
                    mpm_list[3] = ((abv_mode          - 1) % mod) + 2;
                    mpm_list[4] = ((abv_mode + offset    ) % mod) + 2;
                    mpm_list[5] = ((lft_mode          - 1) % mod) + 2;
                } else {
                    mpm_list[3] = ((abv_mode + offset    ) % mod) + 2;
                    mpm_list[4] = ((abv_mode          - 1) % mod) + 2;
                    mpm_list[5] = ((lft_mode + offset    ) % mod) + 2;
                }
            } else {
                if (abv_mode - lft_mode == 1){
                    mpm_list[3] = ((lft_mode + offset    ) % mod) + 2;
                    mpm_list[4] = ((abv_mode          - 1) % mod) + 2;
                    mpm_list[5] = ((lft_mode + offset - 1) % mod) + 2;
                } else if (abv_mode - lft_mode >= 62){
                    mpm_list[3] = ((lft_mode         - 1) % mod) + 2;
                    mpm_list[4] = ((abv_mode + offset   ) % mod) + 2;
                    mpm_list[5] = ( lft_mode              % mod) + 2;
                } else if (abv_mode - lft_mode == 2) {
                    mpm_list[3] = ((lft_mode          - 1) % mod) + 2;
                    mpm_list[4] = ((lft_mode + offset    ) % mod) + 2;
                    mpm_list[5] = ((abv_mode          - 1) % mod) + 2;
                } else {
                    mpm_list[3] = ((lft_mode + offset    ) % mod) + 2;
                    mpm_list[4] = ((lft_mode          - 1) % mod) + 2;
                    mpm_list[5] = ((abv_mode + offset    ) % mod) + 2;
                }
            }
        } else if (lft_mode + abv_mode >= 2){
            if (lft_mode > abv_mode){
                mpm_list[0] = OVINTRA_PLANAR;
                mpm_list[1] = lft_mode;
                mpm_list[2] = ((lft_mode + offset    ) % mod) + 2;
                mpm_list[3] = ((lft_mode          - 1) % mod) + 2;
                mpm_list[4] = ((lft_mode + offset - 1) % mod) + 2;
                mpm_list[5] = ( lft_mode               % mod) + 2;
            } else {
                mpm_list[0] = OVINTRA_PLANAR;
                mpm_list[1] = abv_mode;
                mpm_list[2] = ((abv_mode + offset    ) % mod) + 2;
                mpm_list[3] = ((abv_mode          - 1) % mod) + 2;
                mpm_list[4] = ((abv_mode + offset - 1) % mod) + 2;
                mpm_list[5] = ( abv_mode               % mod) + 2;
            }
        }
    }

    qsort(mpm_list, 6, sizeof(uint8_t), cmpfunc);

    intra_luma_comp += intra_luma_comp >= mpm_list[0];
    intra_luma_comp += intra_luma_comp >= mpm_list[1];
    intra_luma_comp += intra_luma_comp >= mpm_list[2];

    intra_luma_comp += intra_luma_comp >= mpm_list[3];
    intra_luma_comp += intra_luma_comp >= mpm_list[4];
    intra_luma_comp += intra_luma_comp >= mpm_list[5];

    return intra_luma_comp;
}

uint8_t
derive_intra_mode_c(uint8_t cclm_flag, uint8_t mpm_flag, uint8_t mpm_idx,
                    uint8_t luma_mode, uint8_t cclm_idx)
{
    /* FIXME
     *     Modify lm_mode_idx so we fit in mode_list1
     *     Use a switch to clarify ?
     */
    if (cclm_flag) {
        static const uint8_t lm_list[3] = {
            OVINTRA_LM_CHROMA,
            OVINTRA_MDLM_LEFT,
            OVINTRA_MDLM_TOP
        };

        return lm_list[cclm_idx];

    } else if (mpm_flag) {
        static const uint8_t mode_list[4] = {
            OVINTRA_PLANAR,
            OVINTRA_VER,
            OVINTRA_HOR,
            OVINTRA_DC,
        };

        /* FIXME MIP should not be stored or stored as PLANAR for mode derivation
         */
        if (mode_list[mpm_idx] == luma_mode) {
            return OVINTRA_VDIA;
        }

        return mode_list[mpm_idx];
    }

    /* Note defult mode force use of luma */
    return luma_mode;
}


static uint8_t
derive_intra_angular_mode(struct IntraDRVInfo *i_info,
                          uint8_t cu_flags, uint8_t cu_mode_info,
                          int x_pu, int y_pu, int nb_pb_w, int nb_pb_h)
{
    uint8_t mpm_flag = cu_flags & flg_mpm_flag;

    int16_t tr_pos = (x_pu + nb_pb_w - 1) + ((y_pu - 1) << 5);
    int16_t bl_pos = (x_pu - 1) + ((y_pu + nb_pb_h - 1) << 5);

    uint8_t intra_mode;

#if 0
    /*FIXME do not use cclm table here */
    uint8_t mode_trght = y_pu ? pred_ctx->cclm_intra_mode[tr_pos]
                              : OVINTRA_NOT_AVAILABLE;

    uint8_t mode_blft  = x_pu ? pred_ctx->cclm_intra_mode[bl_pos]
                              : pred_ctx->intra_modes_y_luma[y_pu + nb_pb_h - 1];
#else
    /* Zero if not in ctu so there is no actual need of storage 
     * but a reset between CTUs? */
    uint8_t mode_trght = y_pu ? i_info->luma_modes[tr_pos]
                              : OVINTRA_PLANAR;

    /*  Use value from left if not inside
     *  FIXME
     *  Check if we could always use mode line instead of table
     */
    uint8_t mode_blft  = x_pu ? i_info->luma_modes[bl_pos]
                              : i_info->luma_mode_y[y_pu + nb_pb_h - 1];
 
#endif

    if (mpm_flag) {
        uint8_t mpm_idx = cu_mode_info;

        intra_mode = vvc_derive_mpm_list(mpm_idx, mode_blft, mode_trght);

    } else {
        uint8_t mpm_rem = cu_mode_info;

        intra_mode = vvc_derive_mpm_list_sorted(mode_blft, mode_trght, mpm_rem);
    }

    /* TODO store intra mode decision (do not forget to use planar when inter)
     */

    return intra_mode;
}

VVCCU
drv_intra_cu(OVCTUDec *const ctudec, const OVPartInfo *const part_ctx,
             uint8_t x0, uint8_t y0, uint8_t log2_cb_w, uint8_t log2_cb_h,
             VVCCU cu)
{
    struct IntraDRVInfo *const i_info = &ctudec->drv_ctx.intra_info;
    uint8_t x_pu = x0 >> part_ctx->log2_min_cb_s;
    uint8_t y_pu = y0 >> part_ctx->log2_min_cb_s;
    uint8_t nb_pb_w = (1 << log2_cb_w) >> part_ctx->log2_min_cb_s;
    uint8_t nb_pb_h = (1 << log2_cb_h) >> part_ctx->log2_min_cb_s;

    uint8_t pu_shift = part_ctx->log2_min_cb_s - 2;

    uint8_t mip_flag = cu.cu_flags & flg_mip_flag;

    /* Note PLANAR is the default for other modes derivaion */
    uint8_t intra_mode = OVINTRA_PLANAR;
    

    if (mip_flag){
        #if 1
        if (!(cu.cu_opaque & (1 << 7))) {
            uint8_t mip_mode_idx = cu.cu_opaque & 0x3F;
            const struct OVRCNCtx *rcn = &ctudec->rcn_ctx;
            uint16_t *dst = &rcn->ctu_buff.y[x0 + y0 * RCN_CTB_STRIDE];
            vvc_intra_pred_mip(rcn, dst, x0, y0, log2_cb_w, log2_cb_h, mip_mode_idx);
        } else {
            uint8_t mip_mode_idx = cu.cu_opaque & 0x3F;
            const struct OVRCNCtx *rcn = &ctudec->rcn_ctx;
            uint16_t *dst = &rcn->ctu_buff.y[x0 + y0 * RCN_CTB_STRIDE];
            vvc_intra_pred_mip_tr(rcn, dst, x0, y0, log2_cb_w, log2_cb_h, mip_mode_idx);
        }
        #endif

    } else {
        #if 0
        uint8_t mpm_flag = !!(cu.cu_flags & flg_mpm_flag);
        #endif
        uint8_t isp_flag = !!(cu.cu_flags & flg_isp_flag);

        intra_mode = derive_intra_angular_mode(i_info, cu.cu_flags, cu.cu_mode_info,
                                               x_pu, y_pu, nb_pb_w, nb_pb_h);

        cu.cu_mode_idx = intra_mode;

        /* FIXME
         * we do not reconstruct ISP here this will be removed when we 
         * will cal CU reconstruction after TUs of current CU are read
         */

        if (!isp_flag) {
            uint8_t mrl_flag = !!(cu.cu_flags & flg_mrl_flag);
            if (!mrl_flag){
                vvc_intra_pred(&ctudec->rcn_ctx, intra_mode, x0, y0,
                               log2_cb_w, log2_cb_h);
            }

        #if 1
            if (mrl_flag){
                uint8_t mrl_idx = cu.cu_opaque;
                vvc_intra_pred_multi_ref(ctudec, &ctudec->rcn_ctx.ctu_buff.y[0],
                                         RCN_CTB_STRIDE, intra_mode, x0, y0,
                                         log2_cb_w, log2_cb_h,
                                         mrl_idx);
            }
        #endif
        }
    }
    memset(&i_info->luma_mode_x[x_pu], intra_mode, sizeof(uint8_t) * nb_pb_w);
    memset(&i_info->luma_mode_y[y_pu], intra_mode, sizeof(uint8_t) * nb_pb_h);

    for (int i = 0; i < nb_pb_h; i++) {
        memset(&i_info->luma_modes[x_pu + (i << 5) + (y_pu << 5)], intra_mode,
               sizeof(uint8_t) * nb_pb_w);
    }

    /* Store derived actual intra mode */

    /* FIXME  can we update it befor ref is filled ? */
    ctu_field_set_rect_bitfield(&ctudec->rcn_ctx.progress_field, x_pu<<pu_shift,
                                y_pu<<pu_shift, nb_pb_w<<pu_shift, nb_pb_h<<pu_shift);

    if (!ctudec->share && ctudec->coding_tree != &dual_tree
        && ctudec->coding_tree_implicit != &dual_tree_implicit) {
        ctu_field_set_rect_bitfield(&ctudec->rcn_ctx.progress_field_c, x_pu << pu_shift,
                                    y_pu << pu_shift, nb_pb_w << pu_shift, nb_pb_h << pu_shift);
    }

    #if 0
    fill_bs_map(&ctudec->dbf_info.bs2_map, x0, y0, log2_cb_w, log2_cb_h);
    #endif

    return cu;
}

void
vvc_intra_pred(const struct OVRCNCtx *const rcn_ctx,
               uint8_t intra_mode, int x0, int y0,
               int log2_pb_width, int log2_pb_height)
{
    const struct OVBuffInfo *ctu_buff = &rcn_ctx->ctu_buff;
    const struct DCFunctions *dc = &rcn_ctx->rcn_funcs.dc;
    const struct PlanarFunctions *planar = &rcn_ctx->rcn_funcs.planar;

    uint16_t ref_above[(128<<1) + 128]/*={0}*/;
    uint16_t ref_left [(128<<1) + 128]/*={0}*/;

    uint16_t ref_above_filtered[(128<<1) + 128]/*={0}*/;
    uint16_t ref_left_filtered [(128<<1) + 128]/*={0}*/;

    ptrdiff_t dst_stride = ctu_buff->stride;

    const uint16_t *src = &ctu_buff->y[0];
    uint16_t *dst = &ctu_buff->y[x0 + (y0 * dst_stride)];

    uint16_t *ref1 = ref_above + (1 << log2_pb_height);
    uint16_t *ref2 = ref_left + (1 << log2_pb_width);

    fill_ref_left_0(src,dst_stride,ref2,
                    rcn_ctx->progress_field.vfield[x0 >> 2],
                    rcn_ctx->progress_field.hfield[y0 >> 2],
                    x0, y0, log2_pb_width, log2_pb_height, 0);

    fill_ref_above_0(src, dst_stride, ref1,
                     rcn_ctx->progress_field.hfield[y0 >> 2],
                     rcn_ctx->progress_field.vfield[x0 >> 2],
                     x0, y0, log2_pb_width, log2_pb_height, 0);

    switch (intra_mode) {
    case OVINTRA_PLANAR:
    {
        if((log2_pb_height + log2_pb_width) > 5){
            filter_ref_samples(ref1, ref_above_filtered, ref2,
                               (1 << log2_pb_width)+4);
            filter_ref_samples(ref2, ref_left_filtered, ref1,
                               (1 << log2_pb_height)+4);
            ref1 = ref_above_filtered;
            ref2 = ref_left_filtered;
        }

        planar->pdpc[log2_pb_width > 5 || log2_pb_height > 5](ref1, ref2, dst, dst_stride,
                                                              log2_pb_width, log2_pb_height);
        break;
    }
    case OVINTRA_DC:
    {
        dc->pdpc(ref1, ref2, dst, dst_stride, log2_pb_width, log2_pb_height);

        break;
    }
    default:
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
                        int inv_angle = inverse_angle_table[-mode_idx];
                        int inv_angle_sum    = 256;

                        if (use_gauss_filter){
                            filter_ref_samples(ref1, ref_above_filtered + pu_height,
                                               ref2, top_ref_length);
                            filter_ref_samples(ref2, ref_left_filtered + pu_width,
                                               ref1, left_ref_length);
                            ref1 = ref_above_filtered + pu_height;
                            ref2 = ref_left_filtered + pu_width;
                        }

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
                        int inv_angle = inverse_angle_table[-mode_idx];
                        int inv_angle_sum    = 256;

                        if (use_gauss_filter) {
                            filter_ref_samples(ref1, ref_above_filtered + pu_height,
                                               ref2, top_ref_length);
                            filter_ref_samples(ref2, ref_left_filtered + pu_width,
                                               ref1, left_ref_length);
                            ref1 = ref_above_filtered + pu_height;
                            ref2 = ref_left_filtered + pu_width;
                        }

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


static void
vvc_intra_chroma_angular(const uint16_t *const src, uint16_t *const dst,
                         uint16_t *ref_left, uint16_t *ref_above,
                         uint64_t left_col_map, uint64_t top_row_map,
                         int8_t log2_pb_width, int8_t log2_pb_height,
                         int8_t x0, int8_t y0,
                         int8_t intra_mode);

void
vvc_intra_pred_chroma(const struct OVRCNCtx *const rcn_ctx,
                      uint8_t intra_mode, int x0, int y0,
                      int log2_pb_w, int log2_pb_h){

    const struct OVBuffInfo *ctu_buff = &rcn_ctx->ctu_buff;
    const struct RCNFunctions *rcn_func = &rcn_ctx->rcn_funcs;
    const struct DCFunctions *dc = &rcn_ctx->rcn_funcs.dc;
    const struct PlanarFunctions *planar = &rcn_ctx->rcn_funcs.planar;

    uint16_t *const dst_cb = &ctu_buff->cb[(x0) + (y0 * RCN_CTB_STRIDE)];
    uint16_t *const dst_cr = &ctu_buff->cr[(x0) + (y0 * RCN_CTB_STRIDE)];

    const uint16_t *const src_cb = &ctu_buff->cb[0];
    const uint16_t *const src_cr = &ctu_buff->cr[0];

    ptrdiff_t dst_stride = ctu_buff->stride_c;

    /*TODO load ref_sample for cb and cr in same function*/
    uint16_t ref_above[(128<<1) + 128];
    uint16_t ref_left [(128<<1) + 128];
    uint16_t *ref1 = ref_above;
    uint16_t *ref2 = ref_left;

    uint64_t left_col_map = rcn_ctx->progress_field_c.vfield[x0 >> 1];
    uint64_t top_row_map  = rcn_ctx->progress_field_c.hfield[y0 >> 1];

    switch (intra_mode) {
    case OVINTRA_PLANAR://PLANAR
    {
        fill_ref_left_0_chroma(src_cb, dst_stride, ref_left,
                               left_col_map, top_row_map,
                               x0, y0, log2_pb_w, log2_pb_h);

        fill_ref_above_0_chroma(src_cb, dst_stride, ref_above,
                               top_row_map, left_col_map,
                               x0, y0, log2_pb_w, log2_pb_h);

        if (log2_pb_h > 1 && log2_pb_w > 1) {
        // FIXME! vvc_intra_dsp_ctx.planar_pdpc[0](ref1, ref2, dst_cb, dst_stride,
                         // log2_pb_w, log2_pb_h);
               planar->pdpc[0](ref1, ref2, dst_cb, dst_stride, log2_pb_w,
                                          log2_pb_h);
        } else {
                planar->func(ref1, ref2, dst_cb, dst_stride, log2_pb_w,
                                 log2_pb_h);

        }
            fill_ref_left_0_chroma(src_cr, dst_stride, ref_left,
                                   left_col_map, top_row_map,
                                   x0, y0, log2_pb_w, log2_pb_h);

            fill_ref_above_0_chroma(src_cr, dst_stride, ref_above,
                                    top_row_map, left_col_map,
                                    x0, y0, log2_pb_w, log2_pb_h);

            if (log2_pb_h > 1 && log2_pb_w > 1) {
                // FIXED? :vvc_intra_dsp_ctx.planar_pdpc[0](ref1, ref2, dst_cr, dst_stride,
                // log2_pb_w, log2_pb_h);
                planar->pdpc[0](ref1, ref2, dst_cr, dst_stride,
                                      log2_pb_w, log2_pb_h);
            } else {
                planar->func(ref1, ref2, dst_cr, dst_stride, log2_pb_w,
                                 log2_pb_h);
            }
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
            dc->pdpc(ref1, ref2, dst_cb, dst_stride, log2_pb_w,
                              log2_pb_h);
        } else {
            dc->func(ref1, ref2, dst_cb, dst_stride, log2_pb_w,
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
            dc->pdpc(ref1, ref2, dst_cr, dst_stride, log2_pb_w,
                              log2_pb_h);
        } else {
            dc->func(ref1, ref2, dst_cr, dst_stride, log2_pb_w, log2_pb_h);
        }

        break;
    }
    case OVINTRA_LM_CHROMA:
    {
        const uint16_t  *const src_luma = &rcn_ctx->ctu_buff.y[(x0<<1)+((y0<<1)*RCN_CTB_STRIDE)];
        /* FIXME to be replaced by progress fields */
        #if 0
        uint8_t got_left_ctu = left_col_map;// & ((uint64_t)1 << (y0 >> 1));//neighbour & CTU_LFT_FLG;
        uint8_t got_top_ctu  = top_row_map;// & ((uint64_t)1 << (x0 >> 1));//neighbour & CTU_UP_FLG;
        #endif
        uint8_t neighbour = rcn_ctx->ctudec->ctu_ngh_flags;
        uint8_t got_left_ctu = neighbour & CTU_LFT_FLG;
        uint8_t got_top_ctu  = neighbour & CTU_UP_FLG;

        rcn_func->cclm.cclm(src_luma, dst_cb, dst_cr, log2_pb_w, log2_pb_h,
                            y0, got_top_ctu || y0, got_left_ctu || x0);
        break;
    }
    case OVINTRA_MDLM_LEFT:
    {
        uint8_t neighbour = rcn_ctx->ctudec->ctu_ngh_flags;
        uint8_t got_left_ctu = neighbour & CTU_LFT_FLG;
        uint8_t got_top_ctu  = neighbour & CTU_UP_FLG;
        const uint16_t  *const src_luma = &rcn_ctx->ctu_buff.y[(x0<<1)+((y0<<1)*RCN_CTB_STRIDE)];
        #if 0
        /* FIXME to be replaced by progress fields */
        uint8_t got_left_ctu = left_col_map;// & ((uint64_t)1 << (y0 >> 1));//neighbour & CTU_LFT_FLG;
        uint8_t got_top_ctu  = top_row_map;// & ((uint64_t)1 << (x0 >> 1));//neighbour & CTU_UP_FLG;
        #endif

        rcn_func->cclm.mdlm_left(src_luma, dst_cb, dst_cr,
                                 left_col_map, log2_pb_w, log2_pb_h,
                                 x0, y0, x0 || got_left_ctu, y0 || got_top_ctu);
        break;
    }
    case OVINTRA_MDLM_TOP:
    {
        uint8_t neighbour = rcn_ctx->ctudec->ctu_ngh_flags;
        uint8_t got_left_ctu = neighbour & CTU_LFT_FLG;
        uint8_t got_top_ctu  = neighbour & CTU_UP_FLG;
        const uint16_t  *const src_luma = &rcn_ctx->ctu_buff.y[(x0<<1)+((y0<<1)*RCN_CTB_STRIDE)];
        #if 0
        /* FIXME to be replaced by progress fields */
        uint8_t got_left_ctu = left_col_map;// & ((uint64_t)1 << (y0 >> 1));//neighbour & CTU_LFT_FLG;
        uint8_t got_top_ctu  = top_row_map;// & ((uint64_t)1 << (x0 >> 1));//neighbour & CTU_UP_FLG;
        #endif

        rcn_func->cclm.mdlm_top(src_luma, dst_cb, dst_cr,
                                top_row_map, log2_pb_w, log2_pb_h,
                                x0, y0, x0 || got_left_ctu, y0 || got_top_ctu);
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

    int dst_stride = RCN_CTB_STRIDE;
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
