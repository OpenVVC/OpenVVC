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
#include "vcl.h"


static uint8_t
ovcabac_read_ae_alf_idx(OVCABACCtx *const cabac_ctx, uint64_t *const cabac_state,
                        unsigned int tile_group_num_aps)
{
    uint8_t filter_idx = 0;

    if (tile_group_num_aps) {
        uint8_t use_latest_filter = ovcabac_ae_read(cabac_ctx, &cabac_state[ALF_USE_TEMPORAL_FILT_CTX_OFFSET]);

        if (use_latest_filter) {
            if (tile_group_num_aps > 1) {
                filter_idx = vvc_get_cabac_truncated(cabac_ctx, tile_group_num_aps);
            }
            filter_idx += NUM_FIXED_FILTER_SETS;
        } else {
            filter_idx = vvc_get_cabac_truncated(cabac_ctx, NUM_FIXED_FILTER_SETS);
        }

    } else {
        filter_idx = vvc_get_cabac_truncated(cabac_ctx, NUM_FIXED_FILTER_SETS);
    }

    return filter_idx;
}


void
ovcabac_read_ae_alf_ctu(OVCTUDec *const ctudec, uint16_t ctb_rs, uint16_t nb_ctu_w)
{
    uint8_t alf_flags = 0;

    OVCABACCtx *const cabac_ctx = ctudec->cabac_ctx;
    uint64_t *const cabac_state = cabac_ctx->ctx_table;

    struct ALFInfo* alf_info  = &ctudec->alf_info;
    ALFParamsCtu* alf_params_ctu = &alf_info->ctb_alf_params[ctb_rs];

    const uint8_t ctu_ngh_ctx = ctudec->ctu_ngh_flags;

    uint8_t alf_luma_enabled = alf_info->alf_luma_enabled_flag;
    uint8_t alf_cb_enabled   = alf_info->alf_cb_enabled_flag;
    uint8_t alf_cr_enabled   = alf_info->alf_cr_enabled_flag;

    if(!(alf_luma_enabled || alf_cb_enabled || alf_cr_enabled))
        return;

    int ctb_col = ctb_rs % nb_ctu_w;

    const uint8_t ctb_alf_abv = (ctb_rs - nb_ctu_w) >= 0 ? alf_info->ctb_alf_flag_line[ctb_col] : 0;
    const uint8_t ctb_alf_lft = alf_info->left_ctb_alf_flag;

    uint8_t tile_group_num_aps  = alf_info->num_alf_aps_ids_luma;

    if (alf_luma_enabled) {
        uint8_t ctx  = ctu_ngh_ctx & CTU_LFT_FLG && (ctb_alf_lft & 4);
                ctx += ctu_ngh_ctx & CTU_UP_FLG  && (ctb_alf_abv & 4);

        uint8_t alf_luma_flag = ovcabac_ae_read(cabac_ctx, &cabac_state[CTB_ALF_FLAG_CTX_OFFSET + ctx]);

        if (alf_luma_flag) {
            alf_flags |= alf_luma_flag << 2;
            uint8_t alf_idx = ovcabac_read_ae_alf_idx(cabac_ctx, cabac_state, tile_group_num_aps);

            alf_params_ctu->ctb_alf_idx = alf_idx;
        }
    }


    if (alf_cb_enabled) {
        uint8_t nb_alt_c = alf_info->aps_alf_data_c->alf_chroma_num_alt_filters_minus1;
        uint8_t ctx  = ctu_ngh_ctx & CTU_LFT_FLG && (ctb_alf_lft & 2);
                ctx += ctu_ngh_ctx & CTU_UP_FLG  && (ctb_alf_abv & 2);

        uint8_t alf_cb_flag = ovcabac_ae_read(cabac_ctx, &cabac_state[CTB_ALF_FLAG_CTX_OFFSET + 3 + ctx]);

        if (alf_cb_flag) {

            uint8_t cb_alt_idx = 0;
            alf_flags |= alf_cb_flag << 1;

            while (cb_alt_idx < nb_alt_c &&
                   ovcabac_ae_read(cabac_ctx, &cabac_state[CTB_ALF_ALTERNATIVE_CTX_OFFSET])) {
                ++cb_alt_idx;
            }

            alf_params_ctu->cb_alternative = cb_alt_idx;

        }
    }

    if (alf_cr_enabled) {
        uint8_t nb_alt_c = alf_info->aps_alf_data_c->alf_chroma_num_alt_filters_minus1;
        uint8_t ctx  = ctu_ngh_ctx & CTU_LFT_FLG && (ctb_alf_lft & 1);
                ctx += ctu_ngh_ctx & CTU_UP_FLG  && (ctb_alf_abv & 1);

        uint8_t alf_cr_flag = ovcabac_ae_read(cabac_ctx, &cabac_state[CTB_ALF_FLAG_CTX_OFFSET + 6 + ctx]);

        if (alf_cr_flag) {

            uint8_t cr_alt_idx = 0;
            alf_flags |= alf_cr_flag;

            while (cr_alt_idx < nb_alt_c &&
                   ovcabac_ae_read(cabac_ctx, &cabac_state[CTB_ALF_ALTERNATIVE_CTX_OFFSET + 1])) {
                ++cr_alt_idx;
            }

            alf_params_ctu->cr_alternative = cr_alt_idx;

        }
    }

    alf_params_ctu->ctb_alf_flag = alf_flags;

    ovcabac_read_ae_cc_alf_ctu(ctudec, ctb_rs, nb_ctu_w);

    alf_info->left_ctb_alf_flag          = alf_params_ctu->ctb_alf_flag;
    alf_info->ctb_alf_flag_line[ctb_col] = alf_params_ctu->ctb_alf_flag;
}

#define CCALF_CB 0
#define CCALF_CR 1
void
ovcabac_read_ae_cc_alf_ctu(OVCTUDec *const ctudec, uint16_t ctb_rs, uint16_t nb_ctu_w)
{
    struct ALFInfo* alf_info = &ctudec->alf_info;
    uint8_t ccalf_enabled_cb = alf_info->cc_alf_cb_enabled_flag;
    uint8_t ccalf_enabled_cr = alf_info->cc_alf_cr_enabled_flag;

    OVCABACCtx *const cabac_ctx = ctudec->cabac_ctx;
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    const uint8_t ctu_ngh_ctx = ctudec->ctu_ngh_flags;
    ALFParamsCtu* alf_params_ctu = &alf_info->ctb_alf_params[ctb_rs];

    if (ccalf_enabled_cb) {
        const OVALFData* alf_data = alf_info->aps_cc_alf_data_cb;

        int ctb_col = ctb_rs % nb_ctu_w;
        const uint8_t ctb_alf_abv = (ctb_rs - nb_ctu_w) >= 0 ? alf_info->ctb_alf_flag_line[ctb_col] : 0;
        const uint8_t ctb_alf_lft = alf_info->left_ctb_alf_flag;

        const int nb_ccalf_alt = alf_data->alf_cc_cb_filters_signalled_minus1 + 1;

        uint8_t ctx  = (ctu_ngh_ctx & CTU_LFT_FLG) && (ctb_alf_lft & 0x8);
                ctx += (ctu_ngh_ctx & CTU_UP_FLG ) && (ctb_alf_abv & 0x8);

        uint8_t ccalf_idx = ovcabac_ae_read(cabac_ctx, &cabac_state[CC_ALF_FILTER_CONTROL_FLAG_CTX_OFFSET + ctx]);

        if (ccalf_idx) {
            while ((ccalf_idx != nb_ccalf_alt) && ovcabac_bypass_read(cabac_ctx)) {
                ccalf_idx++;
            }
            alf_params_ctu->ctb_alf_flag |= 1 << 3;
            alf_params_ctu->cb_ccalf = ccalf_idx - 1;
        }

    }

    if (ccalf_enabled_cr) {
        const OVALFData* alf_data = alf_info->aps_cc_alf_data_cr;

        int ctb_col = ctb_rs % nb_ctu_w;
        const uint8_t ctb_alf_abv = (ctb_rs - nb_ctu_w) >= 0 ? alf_info->ctb_alf_flag_line[ctb_col] : 0;
        const uint8_t ctb_alf_lft = alf_info->left_ctb_alf_flag;

        const int nb_ccalf_alt = alf_data->alf_cc_cr_filters_signalled_minus1 + 1;

        uint8_t ctx  = (ctu_ngh_ctx & CTU_LFT_FLG) && (ctb_alf_lft & 0x10);
                ctx += (ctu_ngh_ctx & CTU_UP_FLG ) && (ctb_alf_abv & 0x10);
                ctx += 3;

        uint8_t ccalf_idx = ovcabac_ae_read(cabac_ctx, &cabac_state[CC_ALF_FILTER_CONTROL_FLAG_CTX_OFFSET + ctx]);
        if (ccalf_idx) {
            while ((ccalf_idx != nb_ccalf_alt) && ovcabac_bypass_read(cabac_ctx)) {
                ccalf_idx++;
            }
            alf_params_ctu->ctb_alf_flag |= 1 << 4;
            alf_params_ctu->cr_ccalf = ccalf_idx - 1;
        }
    }
}

