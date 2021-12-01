#include <string.h>
#include "ovdefs.h"
#include "ctudec.h"
#include "rcn_lmcs.h"
#include "ovmem.h"
#include "overror.h"

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
    struct SAOInfo* sao_info  = &ctudec->sao_info;
    struct ALFInfo* alf_info  = &ctudec->alf_info;
    struct LMCSInfo* lmcs_info  = &ctudec->lmcs_info;

    if (sao_info->sao_params) {
        ov_free(sao_info->sao_params);
    }

    if (alf_info->ctb_alf_params) {
        ov_free(alf_info->ctb_alf_params);
    }

    if (alf_info->ctb_cc_alf_filter_idx[0]) {
        ov_free(alf_info->ctb_cc_alf_filter_idx[0]);
        ov_free(alf_info->ctb_cc_alf_filter_idx[1]);
    }

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

     ctudec->prev_nb_ctu_w_rect_entry = 0;

     return 0;
}

static void
ctudec_free_intra_line_buff(struct OVRCNCtx *const rcn_ctx)
{
    struct OVBuffInfo *intra_line_b = &rcn_ctx->intra_line_buff;

    if (intra_line_b->y) {
        ov_freep(&intra_line_b->y);
        ov_freep(&intra_line_b->cb);
        ov_freep(&intra_line_b->cr);
    }
}

static void
ctudec_free_filter_buffers(struct OVRCNCtx *const rcn_ctx)
{
    OVSample** saved_rows_sao    = rcn_ctx->filter_buffers.saved_rows_sao;
    OVSample** saved_rows_alf    = rcn_ctx->filter_buffers.saved_rows_alf;
    OVSample** saved_cols        = rcn_ctx->filter_buffers.saved_cols;
    OVSample** filter_region     = rcn_ctx->filter_buffers.filter_region;

    for(int comp = 0; comp < 3; comp++)
    {
        if(filter_region[comp])     ov_freep(&filter_region[comp]);
        if(saved_rows_sao[comp])    ov_freep(&saved_rows_sao[comp]);
        if(saved_rows_alf[comp])    ov_freep(&saved_rows_alf[comp]);
        if(saved_cols[comp])        ov_freep(&saved_cols[comp]);
    }
}

int
ctudec_uninit(OVCTUDec *ctudec)
{
    ctudec_uninit_in_loop_filters(ctudec);

    ctudec_free_filter_buffers(&ctudec->rcn_ctx);
    ctudec_free_intra_line_buff(&ctudec->rcn_ctx);

    ov_freep(&ctudec);

    return 0;
}
