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

#include <stdlib.h>

#include "vcl.h"
#include "cabac_internal.h"
#include "ctudec.h"
#include "nvcl_structures.h"

/* FIXME Separate CABAC reading from offset/parmaeters derivation */

static uint8_t
ovcabac_read_ae_sao_merge_type(OVCABACCtx *const cabac_ctx, uint64_t *const cabac_state,
                               uint8_t neighbour_flags)
{
    uint8_t sao_merge_type = 0;

    if (neighbour_flags & CTU_LFT_FLG) {
        sao_merge_type = ovcabac_ae_read(cabac_ctx, &cabac_state[SAO_MERGE_FLAG_CTX_OFFSET]);
    }

    if (!sao_merge_type && neighbour_flags & CTU_UP_FLG) {
        sao_merge_type = ovcabac_ae_read(cabac_ctx, &cabac_state[SAO_MERGE_FLAG_CTX_OFFSET]);
        sao_merge_type = sao_merge_type << 1;
    }

    return sao_merge_type;
}

static void
ovcabac_read_ae_sao_type_idx(OVCABACCtx *const cabac_ctx, uint64_t *const cabac_state,
                             SAOParamsCtu *sao_ctu,
                             uint8_t sao_enabled_l,
                             uint8_t sao_enabled_c, uint8_t bitdepth_minus8)
{

    int i, k;
    uint8_t nb_bits = (0x1F >> ((bitdepth_minus8 <= 1) + !bitdepth_minus8)) | 0x7;

    if (sao_enabled_l) {
        memset(sao_ctu, 0 ,sizeof(SAOParamsCtu));
        uint8_t ctu_sao_luma_flag = ovcabac_ae_read(cabac_ctx, &cabac_state[SAO_TYPE_IDX_CTX_OFFSET]);

        if (ctu_sao_luma_flag) {
            uint8_t offset_abs[4];
            sao_ctu->type_idx[0] = ovcabac_bypass_read(cabac_ctx) ? SAO_EDGE : SAO_BAND;

            for (i = 0; i < 4; i++) {
                for (k = 0; k < nb_bits; k++) {
                    if (!ovcabac_bypass_read(cabac_ctx)) {
                        break;
                    }
                }
                offset_abs[i] = k;
            }

            if (sao_ctu->type_idx[0] & SAO_BAND) {
                for (k = 0; k < 4; k++) {
                    uint8_t sao_offset_sign = offset_abs[k] && ovcabac_bypass_read(cabac_ctx);
                    if (sao_offset_sign) {
                        sao_ctu->offset_val[0][k] = -offset_abs[k];
                    } else {
                        sao_ctu->offset_val[0][k] = offset_abs[k];
                    }

                }

                sao_ctu->band_position[0] = 0;

                for (i = 1; i < 6; i++) {
                    sao_ctu->band_position[0] |= ovcabac_bypass_read(cabac_ctx) << (5 - i);
                }

            } else {
                sao_ctu->eo_class[0] = ovcabac_bypass_read(cabac_ctx)<<1;
                sao_ctu->eo_class[0] |= ovcabac_bypass_read(cabac_ctx);

                sao_ctu->offset_val[0][0] =  offset_abs[0];
                sao_ctu->offset_val[0][1] =  offset_abs[1];
                sao_ctu->offset_val[0][2] =  0;
                sao_ctu->offset_val[0][3] = -offset_abs[2];
                sao_ctu->offset_val[0][4] = -offset_abs[3];
            }
        }
    }

    if (sao_enabled_c) {
        uint8_t ctu_sao_c_flag = ovcabac_ae_read(cabac_ctx, &cabac_state[SAO_TYPE_IDX_CTX_OFFSET]);
        if (ctu_sao_c_flag) {
            uint8_t offset_abs[5];
            sao_ctu->type_idx[2] = sao_ctu->type_idx[1] = ovcabac_bypass_read(cabac_ctx) ? SAO_EDGE : SAO_BAND;
            for (i = 0; i < 4; i++) {
                for (k = 0; k < nb_bits; k++) {
                    if (!ovcabac_bypass_read(cabac_ctx)) {
                        break;
                    }
                }
                offset_abs[i] = k;
            }

            if (sao_ctu->type_idx[1] & SAO_BAND) {

                for (k = 0; k < 4; k++) {
                    if (offset_abs[k] && ovcabac_bypass_read(cabac_ctx)) {
                        sao_ctu->offset_val[1][k]     = -offset_abs[k];
                    } else {
                        sao_ctu->offset_val[1][k]     = offset_abs[k];
                    }
                }

                sao_ctu->band_position[1] = 0;

                for (i=1; i < 6; i++) {
                    sao_ctu->band_position[1] |= ovcabac_bypass_read(cabac_ctx)<<(5-i);
                }

            } else {//edge
                uint8_t sao_eo_class  = ovcabac_bypass_read(cabac_ctx) << 1;
                        sao_eo_class |= ovcabac_bypass_read(cabac_ctx);

                sao_ctu->eo_class[1] = sao_eo_class;

                sao_ctu->offset_val[1][0] =  offset_abs[0];
                sao_ctu->offset_val[1][1] =  offset_abs[1];
                sao_ctu->offset_val[1][2] =  0;
                sao_ctu->offset_val[1][3] = -offset_abs[2];
                sao_ctu->offset_val[1][4] = -offset_abs[3];
            }

            for (i = 0; i < 4; i++) {
                for (k = 0; k < nb_bits; k++) {
                    if (!ovcabac_bypass_read(cabac_ctx)) {
                        break;
                    }
                }
                offset_abs[i] = k;
            }

            if (sao_ctu->type_idx[2] & SAO_BAND) {

                for (k = 0; k < 4; k++) {
                    if (offset_abs[k] && ovcabac_bypass_read(cabac_ctx)) {
                        sao_ctu->offset_val[2][k]     = -offset_abs[k];
                    } else {
                        sao_ctu->offset_val[2][k]     = offset_abs[k];
                    }
                }

                sao_ctu->band_position[2] = 0;

                for (i = 1; i < 6; i++) {
                    sao_ctu->band_position[2] |= ovcabac_bypass_read(cabac_ctx)<<(5-i);
                }

            } else {
                sao_ctu->eo_class[2] = sao_ctu->eo_class[1];
                sao_ctu->offset_val[2][0] =  offset_abs[0];
                sao_ctu->offset_val[2][1] =  offset_abs[1];
                sao_ctu->offset_val[2][2] =  0;
                sao_ctu->offset_val[2][3] = -offset_abs[2];
                sao_ctu->offset_val[2][4] = -offset_abs[3];
            }
        }
    }
}


void
ovcabac_read_ae_sao_ctu(OVCTUDec *const ctudec, int ctb_rs, uint16_t nb_ctu_w)
{   
    SAOParamsCtu* sao_ctu = &ctudec->sao_info.sao_params[ctb_rs];
    uint8_t sao_enabled_l = ctudec->sao_info.sao_luma_flag;
    uint8_t sao_enabled_c = ctudec->sao_info.sao_chroma_flag;

    if (sao_enabled_l | sao_enabled_c) {
        OVCABACCtx *const cabac_ctx = ctudec->cabac_ctx;
        uint64_t *const cabac_state = cabac_ctx->ctx_table;
        const uint8_t ctu_neighbour_flags = ctudec->ctu_ngh_flags;
        uint8_t merge_type = ovcabac_read_ae_sao_merge_type(cabac_ctx, cabac_state, ctu_neighbour_flags);

        if (!merge_type) {
            uint8_t bitdepth_minus8 = ctudec->bitdepth_minus8;
            ovcabac_read_ae_sao_type_idx(cabac_ctx, cabac_state, sao_ctu,
                                         sao_enabled_l, sao_enabled_c, bitdepth_minus8);
        } else {

            if (merge_type == SAO_MERGE_LEFT) {
                int ctb_left = ctb_rs - 1;
                *sao_ctu = ctudec->sao_info.sao_params[ctb_left];
            }

            if (merge_type == SAO_MERGE_ABOVE) {
                int ctb_above = ctb_rs - nb_ctu_w; 
                *sao_ctu = ctudec->sao_info.sao_params[ctb_above];;
            }
        }
    }
}

