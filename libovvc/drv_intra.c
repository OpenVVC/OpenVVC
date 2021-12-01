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
#include "dbf_utils.h"

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

    int16_t tr_pos = (x_pu + nb_pb_w - 1) + (y_pu << 5) - 32;
    int16_t bl_pos = (x_pu - 1) + ((y_pu + nb_pb_h) << 5) - 32;

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

    uint8_t intra_mode;

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

uint8_t
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
        const struct OVRCNCtx *rcn = &ctudec->rcn_ctx;
        ctudec->rcn_funcs.mip.rcn_intra_mip(rcn, x0, y0, log2_cb_w, log2_cb_h, cu.cu_opaque);
    } else {
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
                ctudec->rcn_funcs.intra_pred(&ctudec->rcn_ctx, &ctudec->rcn_ctx.ctu_buff, intra_mode, x0, y0,
                               log2_cb_w, log2_cb_h);
            }

            if (mrl_flag){
                uint8_t mrl_idx = cu.cu_opaque;
                ctudec->rcn_funcs.intra_pred_mrl(ctudec, ctudec->rcn_ctx.ctu_buff.y,
                                         RCN_CTB_STRIDE, intra_mode, x0, y0,
                                         log2_cb_w, log2_cb_h,
                                         mrl_idx);
            }
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

    if (ctudec->share != 1 && ctudec->coding_tree != &dual_tree
        && ctudec->coding_tree_implicit != &dual_tree_implicit) {
        ctu_field_set_rect_bitfield(&ctudec->rcn_ctx.progress_field_c, x_pu << pu_shift,
                                    y_pu << pu_shift, nb_pb_w << pu_shift, nb_pb_h << pu_shift);
    }

    fill_bs_map(&ctudec->dbf_info.bs2_map, x0, y0, log2_cb_w, log2_cb_h);

    return intra_mode;
}
