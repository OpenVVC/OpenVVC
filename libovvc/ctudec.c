#include "ctudec.h"
#include "ovmem.h"
#include "overror.h"

static void
attach_rcn_ctu_buff(OVCTUDec *const ctudec)
{
     struct OVRCNCtx *rcn_ctx = &ctudec->rcn_ctx;
     struct CTURCNData *rcn_data = &ctu_dec->rcn_ctx.data;
     struct OVBuffInfo *ctu_binfo = &rcn_ctx->ctu_buff;

     ctu_info->data_y  = &rcn_data->y_buff [RCN_CTB_PADDING];
     ctu_info->data_cb = &rcn_data->cb_buff[RCN_CTB_PADDING];
     ctu_info->data_cr = &rcn_data->cr_buff[RCN_CTB_PADDING];

     ctu_info->stride   = RCN_CTB_STRIDE.
     ctu_info->stride_c = RCN_CTB_STRIDE.
}

int
ovdec_decode_ctu(OVVCDec *dec, OVCTUDec *ctu_dec)
{
    return 0;
}

#if 0
static int
ovdec_decode_ctu(const VVCContext *const vvc_ctx,
                 OVCTUDec *const ctudec, int ctb_x, int ctb_y)
{
    uint8_t log2_ctb_s = ctudec->part_ctx->log2_ctu_s;

    const int last_ctu_w = vvc_ctx->frame_width  & ((1 << log2_ctb_s) - 1);
    const int last_ctu_h = vvc_ctx->frame_height & ((1 << log2_ctb_s) - 1);

    int val, ret;

    int is_border_right = ctudec->ctb_x + 1 >= vvc_ctx->nb_ctu_w;
    int is_bottom       = ctudec->ctb_y + 1 >= vvc_ctx->nb_ctu_h;

    if (drv) {
        /* Update derivation maps */
        if (rcn) {
            /* Fill ctb buffers border for intra reconstruction  */
        }
    }

    if (ctudec->filter) {

        /* decode ALF and SAO info */

    }

    if ((!is_border_right || !last_ctu_w) && (!is_bottom || !last_ctu_h)) {

        ret = ctudec->coding_tree(ctudec, ctudec->part_ctx, 0, 0, log2_ctb_s, 0);

        if (rcn) {
            /* Write complete ctu _to frame */
            write_ctu_to_frame(ctudec, vvc_ctx->frame);
        }

    } else {
        //TODO overwrite intra_pred_modes
        const int remaining_w = !is_border_right || !last_ctu_w ? (1 << log2_ctb_s)
                                                                : last_ctu_w;
        const int remaining_h = !is_bottom       || !last_ctu_h ? (1 << log2_ctb_s)
                                                                : last_ctu_h;

        fill_ctu_borders(ctudec, vvc_ctx->frame, remaining_w, remaining_h);

        ret = ctudec->coding_tree_implicit(ctudec, ctudec->part_ctx, 0, 0, log2_ctb_s,
                                           0, remaining_w, remaining_h);

        ctudec->dbf_info.edge_map_ver[remaining_w >> 2] &= -!is_border_right;
        ctudec->dbf_info.edge_map_hor[remaining_h >> 2] &= -!is_bottom;

        if (rcn) {
            /* Write partial ctu _to frame */
            write_ctu_to_frame_border(ctudec, vvc_ctx->frame, remaining_w, remaining_h);
        }

    }

    return 0;
}

static int
ovdec_decode_ctu_border(const VVCContext *const vvc_ctx,
                        OVCTUDec *const ctudec, int ctb_x, int ctb_y)
{
    uint8_t log2_ctb_s = ctudec->part_ctx->log2_ctu_s;

    const int last_ctu_w = vvc_ctx->frame_width  & ((1 << log2_ctb_s) - 1);
    const int last_ctu_h = vvc_ctx->frame_height & ((1 << log2_ctb_s) - 1);

    int val, ret;

    int is_border_right = ctudec->ctb_x + 1 >= vvc_ctx->nb_ctu_w;
    int is_bottom       = ctudec->ctb_y + 1 >= vvc_ctx->nb_ctu_h;

    if (drv) {
        /* Update derivation maps */
        if (rcn) {
            /* Fill ctb buffers border for intra reconstruction  */
        }
    }

    if (ctudec->filter) {

        /* decode ALF and SAO info */

    }

    if ((!is_border_right || !last_ctu_w) && (!is_bottom || !last_ctu_h)) {

        ret = ctudec->coding_tree(ctudec, ctudec->part_ctx, 0, 0, log2_ctb_s, 0);

        if (rcn) {
            /* Write complete ctu _to frame */
            write_ctu_to_frame(ctudec, vvc_ctx->frame);
        }

    } else {
        //TODO overwrite intra_pred_modes
        const int remaining_w = !is_border_right || !last_ctu_w ? (1 << log2_ctb_s)
                                                                : last_ctu_w;
        const int remaining_h = !is_bottom       || !last_ctu_h ? (1 << log2_ctb_s)
                                                                : last_ctu_h;

        fill_ctu_borders(ctudec, vvc_ctx->frame, remaining_w, remaining_h);

        ret = ctudec->coding_tree_implicit(ctudec, ctudec->part_ctx, 0, 0, log2_ctb_s,
                                           0, remaining_w, remaining_h);

        ctudec->dbf_info.edge_map_ver[remaining_w >> 2] &= -!is_border_right;
        ctudec->dbf_info.edge_map_hor[remaining_h >> 2] &= -!is_bottom;

        if (rcn) {
            /* Write partial ctu _to frame */
            write_ctu_to_frame_border(ctudec, vvc_ctx->frame, remaining_w, remaining_h);
        }

    }

    return 0;
}
#endif

int
ctudec_init(OVCTUDec **ctudec_p)
{
     OVCTUDec *ctudec;

     ctudec = ov_mallocz(sizeof(*ctudec));

     if (!ctudec) {
         return OVVC_ENOMEM;
     }

     *ctudec_p = ctudec;

     attach_rcn_ctu_buff(ctudec);

     return 0;
}

int
ctudec_uninit(OVCTUDec *ctudec)
{
     ov_free(ctudec);
     return 0;
}
