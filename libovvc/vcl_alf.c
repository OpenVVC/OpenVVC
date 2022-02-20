/**
 *
 *   OpenVVC is open-source real time software decoder compliant with the 
 *   ITU-T H.266- MPEG-I - Part 3 VVC standard. OpenVVC is developed from 
 *   scratch in C as a library that provides consumers with real time and
 *   energy-aware decoding capabilities under different OS including MAC OS,
 *   Windows, Linux and Android targeting low energy real-time decoding of
 *   4K VVC videos on Intel x86 and ARM platforms.
 * 
 *   Copyright (C) 2020-2022  IETR-INSA Rennes :
 *   
 *   Pierre-Loup CABARAT
 *   Wassim HAMIDOUCHE
 *   Guillaume GAUTIER
 *   Thomas AMESTOY
 *
 *   This library is free software; you can redistribute it and/or
 *   modify it under the terms of the GNU Lesser General Public
 *   License as published by the Free Software Foundation; either
 *   version 2.1 of the License, or (at your option) any later version.
 *
 *   This library is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *   Lesser General Public License for more details.
 *
 *   You should have received a copy of the GNU Lesser General Public
 *   License along with this library; if not, write to the Free Software
 *   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301
 *   USA
 * 
 **/

#include "cabac_internal.h"
#include "ctudec.h"
#include "dec_structures.h"
#include "nvcl_structures.h"
#include "nvcl.h"
#include "rcn_alf.h"


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
ovcabac_read_ae_alf_ctu(OVCTUDec *const ctudec, uint16_t ctb_rs, uint16_t nb_ctu_w)
{
    uint8_t ret_luma = 0;
    uint8_t ret_cb = 0;
    uint8_t ret_cr = 0;
    uint8_t ret;

    uint8_t ctx;
    uint8_t alf_idx = 0;

    struct ALFInfo* alf_info  = &ctudec->alf_info;
    OVCABACCtx *const cabac_ctx = ctudec->cabac_ctx;
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    const uint8_t ctu_neighbour_flags = ctudec->ctu_ngh_flags;
    uint8_t alf_luma_flag   =  alf_info->alf_luma_enabled_flag;
    uint8_t alf_cb_flag     =  alf_info->alf_cb_enabled_flag;
    uint8_t alf_cr_flag     =  alf_info->alf_cr_enabled_flag;

    if(!(alf_luma_flag || alf_cb_flag || alf_cr_flag))
        return;

    const uint8_t left_ctb_alf_flag = alf_info->left_ctb_alf_flag;
    int           ctb_col           = ctb_rs % nb_ctu_w;
    const uint8_t up_ctb_alf_flag   = (ctb_rs - nb_ctu_w) >= 0 ? alf_info->ctb_alf_flag_line[ctb_col] : 0;
   
    uint8_t tile_group_num_aps  = alf_info->num_alf_aps_ids_luma;
    if(alf_luma_flag){
        ctx  = ctu_neighbour_flags & CTU_LFT_FLG  ? ((left_ctb_alf_flag & 4) >> 2) : 0;
        ctx += ctu_neighbour_flags & CTU_UP_FLG   ? ((up_ctb_alf_flag   & 4) >> 2) : 0;
        ret_luma = ovcabac_ae_read(cabac_ctx,&cabac_state[CTB_ALF_FLAG_CTX_OFFSET + ctx]);
        if(ret_luma){
            alf_idx = ovcabac_read_ae_alf_idx(cabac_ctx, cabac_state, tile_group_num_aps);
        }
    }

    uint8_t cb_alternative = 0;
    uint8_t cr_alternative = 0;
    uint8_t num_alf_chroma_alternative;
    if(alf_cb_flag){
        num_alf_chroma_alternative = alf_info->aps_alf_data_c->alf_chroma_num_alt_filters_minus1 + 1;
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
        num_alf_chroma_alternative = alf_info->aps_alf_data_c->alf_chroma_num_alt_filters_minus1 + 1;
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
    alf_info->left_ctb_alf_flag           = ret;
    alf_info->ctb_alf_flag_line[ctb_col]  = ret; 

    ALFParamsCtu* alf_params_ctu = &alf_info->ctb_alf_params[ctb_rs];
    alf_params_ctu->ctb_alf_flag = ret;
    alf_params_ctu->ctb_alf_idx = alf_idx;
    alf_params_ctu->cb_alternative = cb_alternative;
    alf_params_ctu->cr_alternative = cr_alternative;
}

void
ovcabac_read_ae_cc_alf_ctu(OVCTUDec *const ctudec, uint16_t ctb_rs, uint16_t nb_ctu_w)
{   
    struct ALFInfo* alf_info   = &ctudec->alf_info;
    uint8_t cc_alf_cb_flag     = alf_info->cc_alf_cb_enabled_flag;
    uint8_t cc_alf_cr_flag     = alf_info->cc_alf_cr_enabled_flag;

    OVCABACCtx *const cabac_ctx = ctudec->cabac_ctx;
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    const uint8_t ctu_neighbour_flags = ctudec->ctu_ngh_flags;

    for ( int comp_id = 0; comp_id < 2; comp_id++){

        if ((comp_id==0 && cc_alf_cb_flag) || (comp_id==1 && cc_alf_cr_flag)){
            const OVALFData* alf_data = (comp_id==0) ? alf_info->aps_cc_alf_data_cb : alf_info->aps_cc_alf_data_cr;

            const uint8_t left_ctb_cc_alf_flag = alf_info->left_ctb_cc_alf_flag[comp_id];
            int           ctb_col              = ctb_rs % nb_ctu_w;
            const uint8_t up_ctb_cc_alf_flag   = (ctb_rs - nb_ctu_w) >= 0 ? alf_info->ctb_cc_alf_flag_line[comp_id][ctb_col] : 0;

            const int filters_signalled = (comp_id == 0) ? alf_data->alf_cc_cb_filters_signalled_minus1 + 1
                                                            : alf_data->alf_cc_cr_filters_signalled_minus1 + 1;
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
            alf_info->left_ctb_cc_alf_flag[comp_id]              = ret_cc_alf;
            alf_info->ctb_cc_alf_flag_line[comp_id][ctb_col]     = ret_cc_alf;
            alf_info->ctb_cc_alf_filter_idx[comp_id][ctb_rs]  = ret_cc_alf;
        }
    }
}
