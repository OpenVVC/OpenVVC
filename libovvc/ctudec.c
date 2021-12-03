#include <string.h>
#include "ovdefs.h"
#include "ctudec.h"
#include "rcn_lmcs.h"
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

void
ctudec_alloc_intra_line_buff(OVCTUDec *const ctudec, int nb_ctu_w)
{
    struct OVRCNCtx *rcn_ctx = &ctudec->rcn_ctx;
    struct OVBuffInfo *intra_line_b = &rcn_ctx->intra_line_buff;

    const OVPartInfo *const pinfo = ctudec->part_ctx;
    uint8_t log2_ctb_size = pinfo->log2_ctu_s;
    int max_cu_width_l = 1 << log2_ctb_size;
    intra_line_b->stride    = nb_ctu_w*max_cu_width_l ;
    intra_line_b->stride_c  = nb_ctu_w*max_cu_width_l / 2 ;

    //Free and re-alloc when new ctu width for rectangular entry. 
    ctudec_free_intra_line_buff(ctudec);
    intra_line_b->y = ov_malloc(intra_line_b->stride * sizeof(uint16_t));
    intra_line_b->cb = ov_malloc(intra_line_b->stride_c * sizeof(uint16_t));
    intra_line_b->cr = ov_malloc(intra_line_b->stride_c * sizeof(uint16_t));
}

void
ctudec_free_intra_line_buff(OVCTUDec *const ctudec)
{
    struct OVRCNCtx *rcn_ctx = &ctudec->rcn_ctx;
    struct OVBuffInfo *intra_line_b = &rcn_ctx->intra_line_buff;

   if(intra_line_b->y){
        ov_freep(&intra_line_b->y);
        ov_freep(&intra_line_b->cb);
        ov_freep(&intra_line_b->cr);
    }
}

int
ovdec_decode_ctu(OVVCDec *dec, OVCTUDec *ctu_dec)
{
    return 0;
}

void
ctudec_save_last_cols(OVCTUDec *const ctudec, int x_pic_l, int y_pic_l, uint8_t is_border_rect)
{
    if (is_border_rect & OV_BOUNDARY_RIGHT_RECT)
        return;
    
    struct OVFilterBuffers* fb = &ctudec->filter_buffers;
    const int width_l = ( x_pic_l + fb->filter_region_w[0] > ctudec->pic_w ) ? ( ctudec->pic_w - x_pic_l ) : fb->filter_region_w[0];
    const int height_l = ( y_pic_l + fb->filter_region_h[0] > ctudec->pic_h ) ? ( ctudec->pic_h - y_pic_l ) : fb->filter_region_h[0];
    const int margin = fb->margin;

    for(int comp = 0; comp < 3; comp++)
    {
        int16_t* saved_cols = fb->saved_cols[comp];
        int16_t* filter_region = fb->filter_region[comp];
        int stride_filter = fb->filter_region_stride[comp];

        int ratio_luma_chroma = 2;
        int ratio = comp==0 ? 1 : ratio_luma_chroma;        
        const int width = width_l/ratio;
        const int height = height_l/ratio;

        for(int ii=0; ii < height; ii++)
        {
            for(int jj=0; jj < margin; jj++)
            {
                saved_cols[ii*margin + jj] = filter_region[(ii+margin)*stride_filter + width + jj];
            }
        }
    }
    
}

void ctudec_save_last_rows(OVCTUDec *const ctudec, int16_t** saved_rows, int x_l, int x_pic_l, int y_pic_l, uint8_t is_border_rect)
{
    struct OVFilterBuffers* fb = &ctudec->filter_buffers;
    const int width_l = ( x_pic_l + fb->filter_region_w[0] > ctudec->pic_w ) ? ( ctudec->pic_w - x_pic_l ) : fb->filter_region_w[0];
    const int height_l = ( y_pic_l + fb->filter_region_h[0] > ctudec->pic_h ) ? ( ctudec->pic_h - y_pic_l ) : fb->filter_region_h[0];
    const int margin = fb->margin;

    for(int comp = 0; comp < 3; comp++)
    {
        int16_t* saved_rows_comp = saved_rows[comp];
        int16_t* filter_region = fb->filter_region[comp];
        int stride_filter = fb->filter_region_stride[comp];

        int ratio_luma_chroma = 2;
        int ratio = comp==0 ? 1 : ratio_luma_chroma;        
        const int width = width_l/ratio;
        const int height = height_l/ratio;
        const int x = x_l/ratio;

        int stride_rows = fb->saved_rows_stride[comp];
        //save pixels in top left corner of ctu filter
        for(int ii=0; ii < margin; ii++)
        {
            for(int jj=0; jj < margin; jj++)
            {
                // if ( is_border_rect & VVC_BOUNDARY_RIGHT_TILE)
                if ( 0 )
                    filter_region[ii*stride_filter + jj] = saved_rows_comp[ii*stride_rows];
                else
                    filter_region[ii*stride_filter + jj] = saved_rows_comp[ii*stride_rows + x + width - margin + jj];
            }
        }

        if ( is_border_rect & OV_BOUNDARY_BOTTOM_RECT)
            continue;

        for(int ii=0 ; ii < margin; ii++)
        {
            memcpy(&saved_rows_comp[ii*stride_rows + x], &filter_region[(height+ii)*stride_filter + margin], width * sizeof(int16_t));
        }
    } 
}


void
ctudec_extend_filter_region(OVCTUDec *const ctudec, int16_t** saved_rows, int x_l,
                            int x_pic_l, int y_pic_l, uint8_t bnd_msk)
{   

    struct OVFilterBuffers* fb = &ctudec->filter_buffers;

    const int width_l = (x_pic_l + fb->filter_region_w[0] > ctudec->pic_w) ? (ctudec->pic_w - x_pic_l)
                                                                           : fb->filter_region_w[0];
    const int height_l = (y_pic_l + fb->filter_region_h[0] > ctudec->pic_h) ? (ctudec->pic_h - y_pic_l)
                                                                            : fb->filter_region_h[0];
    const int margin = fb->margin;

    for (int comp = 0; comp < 3; comp++) {
        int ratio_luma_chroma = 2;
        int ratio = comp==0 ? 1 : ratio_luma_chroma;        

        const int width = width_l/ratio;
        const int height = height_l/ratio;
        const int x = x_l/ratio;
        const int x_pic = x_pic_l/ratio;
        const int y_pic = y_pic_l/ratio;

        int16_t* saved_rows_comp = saved_rows[comp];
        int16_t* saved_cols = fb->saved_cols[comp];
        int16_t* filter_region = fb->filter_region[comp];
        int stride_filter = fb->filter_region_stride[comp];

        int stride_pic = fb->pic_frame->linesize[comp]/2;
        int16_t* frame = (int16_t*) fb->pic_frame->data[comp] + y_pic*stride_pic + x_pic;
        int stride_rows = fb->saved_rows_stride[comp];
        int filter_region_offset = fb->filter_region_offset[comp];

        //*******************************************************/
        //Copy of entire CTU from frame, before border extension
        if (1) {
            int16_t *dst = &filter_region[filter_region_offset];
            const int16_t *src = &frame[0];
            int cpy_s = sizeof(int16_t) * width;
            for (int ii=0; ii < height; ii++) {
                memcpy(dst, src, cpy_s);
                dst += stride_filter;
                src += stride_pic;
            }
        }

        // //*******************************************************/
        //Left margins
        if (bnd_msk & OV_BOUNDARY_LEFT_RECT) {
            int16_t *dst = &filter_region[margin * stride_filter];
            for (int ii=0; ii < height; ii++) {
                const int16_t pad = dst[margin];
                for (int jj=0; jj < margin; jj++) {
                    dst[jj] = pad;
                }
                dst += stride_filter;
            }
        } else {
            int cpy_s = sizeof(int16_t) * margin;
            int16_t *dst = &filter_region[margin * stride_filter];
            const int16_t *src = saved_cols;
            for (int ii=0; ii < height; ii++) {
                for (int jj=0; jj < margin; jj++) {
                    memcpy(dst, src, cpy_s);
                }
                src += margin;
                dst += stride_filter;
            }
        }

        //Right margins
        int h = margin + height;
        int w = margin + width;
        if (bnd_msk & OV_BOUNDARY_RIGHT_RECT) {
            int16_t *dst = &filter_region[margin * stride_filter + w];
            for (int ii=0; ii < height; ii++) {
                const int16_t pad = dst[-1];
                for (int jj=0; jj < margin; jj++) {
                    dst[jj] = pad;
                }
                dst += stride_filter;
            }
        } else {
            int cpy_s = sizeof(int16_t) * margin;
            const int16_t *src = &frame[width];
                  int16_t *dst = &filter_region[margin * stride_filter + w];
            for (int ii=0; ii < height; ii++) {
                for (int jj=0; jj < margin; jj++) {
                    memcpy(dst, src, cpy_s);
                }
                dst += stride_filter;
                src += stride_pic;
            }
        }

        //*******************************************************/
        //Upper margins
        int x_offset_end = 0;

        if (!(bnd_msk & OV_BOUNDARY_RIGHT_RECT))
            x_offset_end = margin;

        if (bnd_msk & OV_BOUNDARY_UPPER_RECT) {
            int cpy_s = sizeof(int16_t) * (width + x_offset_end);
            int16_t *dst = &filter_region[margin];
            const int16_t *src = dst + (margin * stride_filter);
            for (int ii=0; ii < margin; ii++) {
                memcpy(dst, src, cpy_s);
                dst += stride_filter;
            }
        } else {
            int cpy_s = sizeof(int16_t) * (width + x_offset_end);
            int16_t *dst = &filter_region[margin];
            const int16_t *src = &saved_rows_comp[x];
            for (int ii=0; ii < margin; ii++) {
                memcpy(dst, src, cpy_s);
                dst += stride_filter;
                src += stride_rows;
            }
        }

        //Bottom margins
        if (bnd_msk & OV_BOUNDARY_BOTTOM_RECT) {
            int cpy_s = sizeof(int16_t) * (width + 2 * margin);
                  int16_t *dst = &filter_region[h * stride_filter];
            const int16_t *src = dst - stride_filter;
            for (int ii=0; ii < margin; ii++) {
                memcpy(dst, src, cpy_s);
                dst += stride_filter;
            }
        } else {
            int cpy_s = sizeof(int16_t) * (width + 2 * margin);
                  int16_t *dst = &filter_region[h * stride_filter];
            const int16_t *src = &frame[height * stride_pic - margin];
            for (int ii=0; ii < margin; ii++) {
                memcpy(dst, src, cpy_s);
                dst += stride_filter;
                src += stride_pic;
            }
        }

        //*******************************************************/
        //Fill all corners on boudaries
        if (bnd_msk & OV_BOUNDARY_UPPER_RECT) {
            int cpy_s = sizeof(int16_t) * margin;
                  int16_t *dst = &filter_region[0];
            const int16_t *src = &filter_region[margin * stride_filter];
            for (int ii = 0; ii < margin; ii++) {
                memcpy(dst    , src    , cpy_s);
                memcpy(dst + w, src + w, cpy_s);
                dst += stride_filter;
            }
        }

        if (bnd_msk & OV_BOUNDARY_BOTTOM_RECT) {
            int cpy_s = sizeof(int16_t) * margin;
                  int16_t *dst = &filter_region[h * stride_filter];
            const int16_t *src = dst - stride_filter;
            for (int ii = 0; ii < margin; ii++) {
                memcpy(dst    , src    , cpy_s);
                memcpy(dst + w, src + w, cpy_s);
                dst += stride_filter;
            }
        }

        if (bnd_msk & OV_BOUNDARY_LEFT_RECT) {
            int16_t *src  = &filter_region[0];
            int16_t *src2 = &filter_region[h * stride_filter];
            for (int ii = 0; ii < margin; ii++) {
                const int16_t pad  = src [margin];
                const int16_t pad2 = src2[margin];
                for (int jj = 0; jj < margin; jj++) {
                    src [jj] = pad;
                    src2[jj] = pad2;
                }
                src  += stride_filter;
                src2 += stride_filter;
            }
        }

        if (bnd_msk & OV_BOUNDARY_RIGHT_RECT) {
            int16_t *src  = &filter_region[w];
            int16_t *src2 = &filter_region[w + (h * stride_filter)];
            for (int ii=0; ii < margin; ii++) {
                const int16_t pad  = src [-1];
                const int16_t pad2 = src2[-1];
                for (int jj=0; jj < margin; jj++) {
                    src [jj] = pad;
                    src2[jj] = pad2;
                }
                src  += stride_filter;
                src2 += stride_filter;
            }
        }
    }
}

void ctudec_alloc_filter_buffers(OVCTUDec *const ctudec, int nb_ctu_w, int margin)
{   
    const OVPartInfo *const pinfo = ctudec->part_ctx;
    uint8_t log2_ctb_size = pinfo->log2_ctu_s;
    int ctu_s = 1 << log2_ctb_size;

    struct OVFilterBuffers* fb = &ctudec->filter_buffers;
    int16_t** saved_rows_sao = fb->saved_rows_sao;
    int16_t** saved_rows_alf = fb->saved_rows_alf;
    int16_t** saved_cols    = fb->saved_cols;
    int16_t** filter_region = fb->filter_region;
    fb->margin = margin;

    for(int comp = 0; comp < 3; comp++) {   
        int ratio_luma_chroma = 2;
        int ratio = comp==0 ? 1 : ratio_luma_chroma;

        fb->filter_region_w[comp]        = ctu_s / ratio ;
        fb->filter_region_h[comp]        = ctu_s / ratio ;
        fb->filter_region_stride[comp]   = ctu_s / ratio + 2 * margin ;
        fb->filter_region_offset[comp]   = margin * fb->filter_region_stride[comp] + margin;

        if (!filter_region[comp]) {
            int ext_size = fb->filter_region_stride[comp] * (fb->filter_region_h[comp] + 2 * margin + 1);
            filter_region[comp] = ov_malloc(ext_size * sizeof(int16_t));
            saved_cols[comp]    = ov_malloc(fb->filter_region_h[comp] * margin * sizeof(int16_t));
        } 

        fb->saved_rows_stride[comp] = nb_ctu_w * ctu_s / ratio; ;
        if (saved_rows_sao[comp]) {
            ov_freep(&saved_rows_sao[comp]);
        }
        saved_rows_sao[comp] = ov_mallocz(margin * fb->saved_rows_stride[comp] * sizeof(int16_t));

        if(saved_rows_alf[comp]){
            ov_freep(&saved_rows_alf[comp]);
        }
        saved_rows_alf[comp] = ov_mallocz(margin * fb->saved_rows_stride[comp] * sizeof(int16_t));
    }

}

void ctudec_free_filter_buffers(OVCTUDec *const ctudec)
{
    int16_t** saved_rows_sao    = ctudec->filter_buffers.saved_rows_sao;
    int16_t** saved_rows_alf    = ctudec->filter_buffers.saved_rows_alf;
    int16_t** saved_cols        = ctudec->filter_buffers.saved_cols;
    int16_t** filter_region     = ctudec->filter_buffers.filter_region;

    for(int comp = 0; comp < 3; comp++)
    {
        if(filter_region[comp])     ov_freep(&filter_region[comp]);
        if(saved_rows_sao[comp])    ov_freep(&saved_rows_sao[comp]);
        if(saved_rows_alf[comp])    ov_freep(&saved_rows_alf[comp]);
        if(saved_cols[comp])        ov_freep(&saved_cols[comp]);
    }
}


int
ctudec_init_in_loop_filters(OVCTUDec *const ctudec, const OVPS *const prms)
{
    const OVSPS *const sps = prms->sps;
    const OVSH *const sh = prms->sh;
    const OVPH *const ph = prms->ph;

    uint16_t pic_w = sps->sps_pic_width_max_in_luma_samples;
    uint16_t pic_h = sps->sps_pic_height_max_in_luma_samples;
    uint8_t log2_ctb_s = sps->sps_log2_ctu_size_minus5 + 5;
    int nb_ctb_pic_w = (pic_w + ((1 << log2_ctb_s) - 1)) >> log2_ctb_s;
    int nb_ctb_pic_h = (pic_h + ((1 << log2_ctb_s) - 1)) >> log2_ctb_s;

    //Init SAO info and ctu params
    struct SAOInfo* sao_info  = &ctudec->sao_info;
    sao_info->sao_luma_flag   =  sh->sh_sao_luma_used_flag;
    sao_info->sao_chroma_flag =  sh->sh_sao_chroma_used_flag;
    sao_info->chroma_format_idc = sps->sps_chroma_format_idc;
    if(sao_info->sao_luma_flag || sao_info->sao_chroma_flag){
        if(!sao_info->sao_params){
            sao_info->sao_params = ov_mallocz(sizeof(SAOParamsCtu) * nb_ctb_pic_w * nb_ctb_pic_h);
        } else {
            memset(sao_info->sao_params,0,sizeof(SAOParamsCtu) * nb_ctb_pic_w * nb_ctb_pic_h);
        }
    }

    //Init ALF info and ctu params
    struct ALFInfo* alf_info  = &ctudec->alf_info;
    alf_info->alf_luma_enabled_flag = sh->sh_alf_enabled_flag;
    alf_info->alf_cb_enabled_flag = sh->sh_alf_cb_enabled_flag;
    alf_info->alf_cr_enabled_flag = sh->sh_alf_cr_enabled_flag;

    if(alf_info->alf_luma_enabled_flag || alf_info->alf_cb_enabled_flag || alf_info->alf_cr_enabled_flag){
        alf_info->num_alf_aps_ids_luma  = sh->sh_num_alf_aps_ids_luma;
        for (int i = 0; i < alf_info->num_alf_aps_ids_luma; i++)
        {
            alf_info->aps_alf_data[i] = &prms->aps_alf[i]->aps_alf_data;
        }
        alf_info->aps_alf_data_c = &prms->aps_alf_c->aps_alf_data;
        if(!alf_info->ctb_alf_params){
            alf_info->ctb_alf_params = ov_malloc(sizeof(ALFParamsCtu) * nb_ctb_pic_w * nb_ctb_pic_h);
        } else {
            memset(alf_info->ctb_alf_params, 0, sizeof(ALFParamsCtu) * nb_ctb_pic_w * nb_ctb_pic_h);
        }

        //Initialization of ALF reconstruction structures
        RCNALF* alf = &alf_info->rcn_alf;
        uint8_t luma_flag = alf_info->alf_luma_enabled_flag;
        uint8_t chroma_flag = alf_info->alf_cb_enabled_flag || alf_info->alf_cr_enabled_flag;
        ctudec->rcn_ctx.rcn_funcs.alf.rcn_alf_reconstruct_coeff_APS(alf, ctudec, luma_flag, chroma_flag);
    }

    //Init CC ALF ctu params
    alf_info->cc_alf_cb_enabled_flag = sh->sh_alf_cc_cb_enabled_flag;
    alf_info->cc_alf_cr_enabled_flag = sh->sh_alf_cc_cr_enabled_flag;
    if(alf_info->cc_alf_cb_enabled_flag || alf_info->cc_alf_cr_enabled_flag){
        alf_info->aps_cc_alf_data_cb   = &prms->aps_cc_alf_cb->aps_alf_data;
        alf_info->aps_cc_alf_data_cr = &prms->aps_cc_alf_cr->aps_alf_data;
        if(!alf_info->ctb_cc_alf_filter_idx[0]){
            alf_info->ctb_cc_alf_filter_idx[0] = ov_malloc(sizeof(uint8_t) * nb_ctb_pic_w * nb_ctb_pic_h);
            alf_info->ctb_cc_alf_filter_idx[1] = ov_malloc(sizeof(uint8_t) * nb_ctb_pic_w * nb_ctb_pic_h);
        } else {
            memset(alf_info->ctb_cc_alf_filter_idx[0], 0, sizeof(uint8_t) * nb_ctb_pic_w * nb_ctb_pic_h);
            memset(alf_info->ctb_cc_alf_filter_idx[1], 0, sizeof(uint8_t) * nb_ctb_pic_w * nb_ctb_pic_h);
        }
    }

    //Init LMCS info and output pivots
    struct LMCSInfo* lmcs_info   = &ctudec->lmcs_info;
    lmcs_info->lmcs_enabled_flag = ph->ph_lmcs_enabled_flag;
    lmcs_info->scale_c_flag      = ph->ph_chroma_residual_scale_flag;
    if(sh->sh_lmcs_used_flag){
        const struct OVLMCSData* lmcs_data = &prms->aps_lmcs->aps_lmcs_data;
        ctudec->rcn_ctx.rcn_funcs.rcn_init_lmcs(lmcs_info, lmcs_data);
    }

    return 0;
}

void
ctudec_uninit_in_loop_filters(OVCTUDec *const ctudec)
{
    //Uninit SAO info and ctu params
    struct SAOInfo* sao_info  = &ctudec->sao_info;
    if(sao_info->sao_params){
        ov_free(sao_info->sao_params);
    }

    //Uninit ALF info and ctu params
    struct ALFInfo* alf_info  = &ctudec->alf_info;
    if(alf_info->ctb_alf_params){
        ov_free(alf_info->ctb_alf_params);
        // rcn_alf_destroy(&alf_info->rcn_alf, ctb_size);
    }

    //Uninit CC ALF ctu params
    if(alf_info->ctb_cc_alf_filter_idx[0]){
        ov_free(alf_info->ctb_cc_alf_filter_idx[0]);
        ov_free(alf_info->ctb_cc_alf_filter_idx[1]);
    }

    //Uninit LMCS info and output pivots
    struct LMCSInfo* lmcs_info  = &ctudec->lmcs_info;
    if (lmcs_info->luts) {
        ov_free(lmcs_info->luts);
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
     
     ctudec->prev_nb_ctu_w_rect_entry = 0;

     return 0;
}

int
ctudec_uninit(OVCTUDec *ctudec)
{
    ctudec_uninit_in_loop_filters(ctudec);

    ctudec_free_filter_buffers(ctudec);    
    ctudec_free_intra_line_buff(ctudec);

    ov_freep(&ctudec);

    return 0;
}
