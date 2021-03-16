#include <string.h>
#include "ctudec.h"
#include "ovmem.h"
#include "overror.h"

static void
attach_rcn_ctu_buff(OVCTUDec *const ctudec)
{
     struct OVRCNCtx *rcn_ctx = &ctudec->rcn_ctx;
     struct CTURCNData *rcn_data = &ctudec->rcn_ctx.data;
     struct OVBuffInfo *ctu_binfo = &rcn_ctx->ctu_buff;

     ctu_binfo->y  = &rcn_data->y_buff [RCN_CTB_PADDING];
     ctu_binfo->cb = &rcn_data->cb_buff[RCN_CTB_PADDING];
     ctu_binfo->cr = &rcn_data->cr_buff[RCN_CTB_PADDING];

     ctu_binfo->stride   = RCN_CTB_STRIDE;
     ctu_binfo->stride_c = RCN_CTB_STRIDE;
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



void ctudec_extend_filter_region(OVCTUDec *const ctudec)
{   
    struct OVFilterBuffers fb = ctudec->filter_buffers;
    int16_t** saved_rows = fb.saved_rows;
    int16_t** saved_cols = fb.saved_cols;
    int16_t** filter_region = fb.filter_region;

    for(int comp = 0; comp < 3; comp++)
    {
        int ratio_luma_chroma = ctudec->rcn_ctx.frame_buff.stride / ctudec->rcn_ctx.frame_buff.stride_c;
        int ratio = comp==0 ? 1 : ratio_luma_chroma;        

        int stride_pic = fb.pic_frame->linesize[comp]/2;
        int16_t* frame_comp = (int16_t*) fb.pic_frame->data[comp];

        //*******************************************************/
        //Copy of entire frame in filter buffer
        for(int ii=0; ii < ctudec->pic_h / ratio; ii++)
        {
            memcpy(&filter_region[comp][ii*fb.filter_region_stride[comp] + fb.filter_region_offset[comp]], &frame_comp[ii * stride_pic], sizeof(int16_t)* stride_pic);
        }
    }
}



void ctudec_create_filter_buffers(OVCTUDec *const ctudec, struct Frame *pic_frame, int nb_ctu_w, int margin)
{   
    const OVPartInfo *const pinfo = ctudec->part_ctx;
    uint8_t log2_ctb_size = pinfo->log2_ctu_s;
    int max_cu_width_l = 1 << log2_ctb_size;

    struct OVFilterBuffers* fb = &ctudec->filter_buffers;
    int16_t** saved_rows = fb->saved_rows;
    int16_t** saved_cols = fb->saved_cols;
    int16_t** filter_region = fb->filter_region;
    fb->margin    = margin;
    fb->pic_frame = pic_frame;

    for(int comp = 0; comp < 3; comp++)
    {
        int ratio_luma_chroma = ctudec->rcn_ctx.frame_buff.stride / ctudec->rcn_ctx.frame_buff.stride_c;
        int ratio = comp==0 ? 1 : ratio_luma_chroma;

        // fb->saved_rows_h[comp]        = margin ;
        // fb->saved_rows_stride[comp]   = nb_ctu_w * max_cu_width_l /ratio + 2 * margin ;
        // if(!saved_rows[comp]){
        //     saved_rows[comp] = ov_malloc(fb->saved_rows_h[comp] * fb->saved_rows_stride[comp] * sizeof(int16_t));
        // } else {
        //     memset(saved_rows[comp],0, fb->saved_rows_h[comp] * fb->saved_rows_stride[comp] * sizeof(int16_t));
        // }

        // fb->saved_cols_h[comp]        = max_cu_width_l /ratio ;
        // fb->saved_cols_stride[comp]   = margin ;
        // if(!saved_cols[comp]){
        //     saved_cols[comp] = ov_malloc(fb->saved_cols_h[comp] * fb->saved_cols_stride[comp] * sizeof(int16_t));
        // } else {
        //     memset(saved_cols[comp],0, fb->saved_cols_h[comp] * fb->saved_cols_stride[comp] * sizeof(int16_t));
        // }
        
        // fb->filter_region_h[comp]        = max_cu_width_l /ratio + 2*margin ;
        // fb->filter_region_stride[comp]   = max_cu_width_l /ratio + 2*margin ;
        // if(!filter_region[comp]){
        //     filter_region[comp] = ov_malloc( fb->filter_region_stride[comp] * fb->filter_region_h[comp] * sizeof(int16_t));
        // } else {
        //     memset(filter_region[comp],0, fb->filter_region_stride[comp] * fb->filter_region_h[comp] * sizeof(int16_t));
        // }

        fb->filter_region_h[comp]        = ctudec->pic_h /ratio + 2 * margin ;
        fb->filter_region_stride[comp]   = ctudec->pic_w /ratio + 2 * margin ;
        fb->filter_region_offset[comp]   = margin * fb->filter_region_stride[comp] + margin;
        if(!filter_region[comp]){
            filter_region[comp] = ov_malloc( fb->filter_region_stride[comp] * fb->filter_region_h[comp] * sizeof(int16_t));
        } else {
            memset(filter_region[comp],0, fb->filter_region_stride[comp] * fb->filter_region_h[comp] * sizeof(int16_t));
        }
    }
}

void free_filter_buffers(OVCTUDec *const ctudec)
{
    int16_t** saved_rows    = ctudec->filter_buffers.saved_rows;
    int16_t** saved_cols    = ctudec->filter_buffers.saved_cols;
    int16_t** filter_region = ctudec->filter_buffers.filter_region;

    for(int comp = 0; comp < 3; comp++)
    {
        if(saved_rows[comp])    ov_free(saved_rows[comp]);
        if(saved_cols[comp])    ov_free(saved_cols[comp]);
        if(filter_region[comp]) ov_free(filter_region[comp]);
    }
}


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
    free_filter_buffers(ctudec);
    ov_free(ctudec);
    return 0;
}

