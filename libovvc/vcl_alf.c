#include "cabac_internal.h"
#include "ctudec.h"
#include "dec_structures.h"
#include "nvcl_structures.h"
#include "nvcl_utils.h"
#include "nvcl.h"


static uint8_t
ovcabac_read_ae_alf_idx(OVCABACCtx *const cabac_ctx, uint64_t *const cabac_state,
                   unsigned int tile_group_num_aps)
{
    uint8_t filter_idx = 0;
    if (tile_group_num_aps) {
        // uint8_t use_latest_filter = vvc_get_cabac_inline(cabac_ctx, &cabac_state[ALF_USE_TEMPORAL_FILT_CTX_OFFSET]);
        uint8_t use_latest_filter = ovcabac_ae_read(cabac_ctx, &cabac_state[ALF_USE_TEMPORAL_FILT_CTX_OFFSET]);

        if (use_latest_filter){
            if (tile_group_num_aps > 1) {
                filter_idx = vvc_get_cabac_truncated(cabac_ctx, tile_group_num_aps );
            }
            filter_idx += NUM_FIXED_FILTER_SETS;
        } 
        else {
            filter_idx = vvc_get_cabac_truncated(cabac_ctx, NUM_FIXED_FILTER_SETS );
        }
    }
    else{
        filter_idx = vvc_get_cabac_truncated(cabac_ctx, NUM_FIXED_FILTER_SETS );
    }
    return filter_idx;
}


void 
ovcabac_read_ae_alf_ctu(OVCTUDec *const ctudec, const OVPS *const prms, uint16_t ctb_rs, uint16_t nb_ctu_w)
{
    uint8_t ret_luma = 0;
    uint8_t ret_cb = 0;
    uint8_t ret_cr = 0;
    uint8_t ret;

    uint8_t ctx;
    uint8_t alf_idx;

    OVCABACCtx *const cabac_ctx = ctudec->cabac_ctx;
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    const uint8_t ctu_neighbour_flags = ctudec->ctu_ngh_flags;
    uint8_t alf_luma_flag   =  prms->sh->sh_alf_enabled_flag;
    uint8_t alf_cb_flag     =  prms->sh->sh_alf_cb_enabled_flag;
    uint8_t alf_cr_flag     =  prms->sh->sh_alf_cr_enabled_flag;

    const uint8_t left_ctb_alf_flag = ctudec->left_ctb_alf_flag;
    int           ctb_col           = ctb_rs % nb_ctu_w;
    const uint8_t up_ctb_alf_flag   = (ctb_rs - nb_ctu_w) >= 0 ? ctudec->ctb_alf_flag_line[ctb_col] : 0;
   
    //TODO: verifier nom et signification de tile_group_num_aps.
    uint8_t tile_group_num_aps  = prms->sh->sh_num_alf_aps_ids_luma;
    uint8_t num_alf_chroma_alternative = prms->aps->aps_alf_data.alf_chroma_num_alt_filters_minus1 + 1;

    if(alf_luma_flag){
        ctx  = ctu_neighbour_flags & CTU_LFT_FLG  ? ((left_ctb_alf_flag & 4) >> 2) : 0;
        ctx += ctu_neighbour_flags & CTU_UP_FLG   ? ((up_ctb_alf_flag   & 4) >> 2) : 0;
        ret_luma = ovcabac_ae_read(cabac_ctx,&cabac_state[CTB_ALF_FLAG_CTX_OFFSET + ctx]);
        if(ret_luma){
            alf_idx = ovcabac_read_ae_alf_idx(cabac_ctx, cabac_state, tile_group_num_aps);
        }
    }

    uint8_t cb_alternative, cr_alternative = 0;
    if(alf_cb_flag){
        int decoded = 0;
        ctx  = ctu_neighbour_flags & CTU_LFT_FLG ? ((left_ctb_alf_flag & 2) >> 1) : 0;
        ctx += ctu_neighbour_flags & CTU_UP_FLG   ? ((up_ctb_alf_flag   & 2) >> 1) : 0;
        ret_cb = ovcabac_ae_read(cabac_ctx,&cabac_state[CTB_ALF_FLAG_CTX_OFFSET + 3 + ctx]);
        while (ret_cb && decoded < num_alf_chroma_alternative - 1 && ovcabac_ae_read(cabac_ctx,
                                                                         &cabac_state[CTB_ALF_ALTERNATIVE_CTX_OFFSET])){
            ++decoded;
        }
        cb_alternative = decoded;
    }
    if(alf_cr_flag){
        int decoded = 0;
        ctx  = ctu_neighbour_flags & CTU_LFT_FLG ? (left_ctb_alf_flag & 1) : 0;
        ctx += ctu_neighbour_flags & CTU_UP_FLG   ? (up_ctb_alf_flag   & 1) : 0;
        ret_cr = ovcabac_ae_read(cabac_ctx,&cabac_state[CTB_ALF_FLAG_CTX_OFFSET + 6 + ctx]);
        while (ret_cr && decoded < num_alf_chroma_alternative - 1 && ovcabac_ae_read(cabac_ctx,
                                                                         &cabac_state[CTB_ALF_ALTERNATIVE_CTX_OFFSET + 1])){
            ++decoded;
        }
        cr_alternative = decoded;
    }
    ret = (ret_luma << 2) | (ret_cb << 1) | ret_cr;
    ctudec->left_ctb_alf_flag           = ret;
    ctudec->ctb_alf_flag_line[ctb_col]  = ret; 
    // AlfCtuParam alf_ctu_param;
    // alf_ctu_param.flags_y_cb_cr = ret;
    // alf_ctu_param.alf_idx = alf_idx;
    // alf_ctu_param.cb_alternative = cb_alternative;
    // alf_ctu_param.cr_alternative = cr_alternative;
    // return alf_ctu_param;

    //TODO: decider dans quelle stucture de donnees mettre les parametres alf
    // prms->ctb_alf_flag[lc_ctx->ctb_x] = alf_ctu_param.flags_y_cb_cr;
    // prms->ctb_alf_idx[lc_ctx->ctb_x] = alf_ctu_param.alf_idx;
    // prms->ctb_alf_flag_pic[ctb_rs] = alf_ctu_param.flags_y_cb_cr;
    // prms->ctb_alf_idx_pic[ctb_rs] = alf_ctu_param.alf_idx;
    // prms->cb_alternative_pic[ctb_rs] = alf_ctu_param.cb_alternative;
    // prms->cr_alternative_pic[ctb_rs] = alf_ctu_param.cr_alternative;
}

void
ovcabac_read_ae_cc_alf_ctu(OVCTUDec *const ctudec, const OVPS *const prms, uint16_t ctb_rs, uint16_t nb_ctu_w)
{
    OVALFData alf_data         = prms->aps->aps_alf_data;
    uint8_t cc_alf_cb_flag     = prms->sh->sh_alf_cc_cb_enabled_flag;
    uint8_t cc_alf_cr_flag     = prms->sh->sh_alf_cc_cr_enabled_flag;

    OVCABACCtx *const cabac_ctx = ctudec->cabac_ctx;
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    const uint8_t ctu_neighbour_flags = ctudec->ctu_ngh_flags;

    for ( int comp_id = 0; comp_id < 2; comp_id++){

        if ((comp_id==0 && cc_alf_cb_flag) || (comp_id==1 && cc_alf_cr_flag)){

            const uint8_t left_ctb_cc_alf_flag = ctudec->left_ctb_cc_alf_flag[comp_id];
            int           ctb_col              = ctb_rs % nb_ctu_w;
            const uint8_t up_ctb_cc_alf_flag   = (ctb_rs - nb_ctu_w) >= 0 ? ctudec->ctb_cc_alf_flag_line[comp_id][ctb_col] : 0;

            const int filters_signalled = (comp_id == 0) ? alf_data.alf_cc_cb_filters_signalled_minus1 + 1
                                                            : alf_data.alf_cc_cr_filters_signalled_minus1 + 1;
            uint8_t  ctx = 0;
            if (ctu_neighbour_flags & CTU_LFT_FLG){
                ctx += left_ctb_cc_alf_flag ? 1 : 0;
            }
            if (ctu_neighbour_flags & CTU_UP_FLG){
                ctx += up_ctb_cc_alf_flag ? 1 : 0;
            }
            ctx += ( comp_id == 1 ) ? 3 : 0;

            int ret_cc_alf = ovcabac_ae_read(cabac_ctx,&cabac_state[CC_ALF_FILTER_CONTROL_FLAG_CTX_OFFSET + ctx]);
            if ( ret_cc_alf ){
                while ( ( ret_cc_alf != filters_signalled ) && ovcabac_bypass_read(cabac_ctx)){
                  ret_cc_alf++;
                }
            }
            ctudec->left_ctb_cc_alf_flag[comp_id]           = ret_cc_alf;
            ctudec->ctb_cc_alf_flag_line[comp_id][ctb_col]  = ret_cc_alf;
            //TODO: voir ou stocker le filter control id (si besoin de stocker)
            // cc_alf_filter_control_id[ctb_rs] = ret_cc_alf;
        }

    }
}
