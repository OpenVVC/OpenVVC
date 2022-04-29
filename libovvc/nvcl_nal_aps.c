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

#include "ovmem.h"
#include "ovutils.h"
#include "overror.h"

#include "nvcl.h"
#include "nvcl_utils.h"
#include "nvcl_structures.h"

typedef struct OVScalingList
{
    uint8_t scaling_list_copy_mode_flag[28];
    uint8_t scaling_list_pred_mode_flag[28];
    int16_t scaling_list_pred_id_delta[28];
    int16_t scaling_list_dc_coef[28 - 14];
    int16_t scaling_list_delta_coef[28][64];
} OVScalingList;

enum APSType {
   APS_ALF          = 0,
   APS_LMCS         = 1,
   APS_SCALING_LIST = 2
};

static int
validate_aps(OVNVCLReader *rdr, OVAPS *const aps)
{
    /* TODO various check on limitation and max sizes */
    return 1;
}


void 
nvcl_read_alf_data(OVNVCLReader *const rdr, struct OVALFData* alf_data,
                   uint8_t aps_chroma_present_flag)
{
    uint8_t sign = 0;

    alf_data->alf_luma_filter_signal_flag = nvcl_read_flag(rdr);

    if (aps_chroma_present_flag) {
        alf_data->alf_chroma_filter_signal_flag = nvcl_read_flag(rdr);
        alf_data->alf_cc_cb_filter_signal_flag = nvcl_read_flag(rdr);
        alf_data->alf_cc_cr_filter_signal_flag = nvcl_read_flag(rdr);
    }

    if (alf_data->alf_luma_filter_signal_flag) {
        alf_data->alf_luma_clip_flag = nvcl_read_flag(rdr);
        alf_data->alf_luma_num_filters_signalled_minus1 = nvcl_read_u_expgolomb(rdr);
        if (alf_data->alf_luma_num_filters_signalled_minus1 > 0) {
            for (int filtIdx = 0; filtIdx < MAX_NUM_ALF_CLASSES; filtIdx++) {
                alf_data->alf_luma_coeff_delta_idx[filtIdx] = 
                nvcl_read_bits(rdr, 31 - __builtin_clz(alf_data->alf_luma_num_filters_signalled_minus1) + 1);
            }
        }

        for (int sfIdx = 0; sfIdx <= alf_data->alf_luma_num_filters_signalled_minus1; sfIdx++) {
            for (int j = 0; j < 12; j++) {
                alf_data->alf_luma_coeff[sfIdx][j] = nvcl_read_u_expgolomb(rdr);
                if (alf_data->alf_luma_coeff[sfIdx][j]) {
                    sign = nvcl_read_bits(rdr, 1);
                    if (sign) alf_data->alf_luma_coeff[sfIdx][j] *= -1;
                }
            }
        }

        if (alf_data->alf_luma_clip_flag) {
            for (int sfIdx = 0; sfIdx <= alf_data->alf_luma_num_filters_signalled_minus1; sfIdx++) {
                for (int j = 0; j < 12; j++) {
                    alf_data->alf_luma_clip_idx[sfIdx][j] = nvcl_read_bits(rdr, 2);
                }
            }
        }
    }

    if (alf_data->alf_chroma_filter_signal_flag) {
        alf_data->alf_chroma_clip_flag = nvcl_read_flag(rdr);
        alf_data->alf_chroma_num_alt_filters_minus1 = nvcl_read_u_expgolomb(rdr);
        for (int altIdx = 0; altIdx <= alf_data->alf_chroma_num_alt_filters_minus1; altIdx++) {
            for (int j = 0; j < 6; j++) {
                alf_data->alf_chroma_coeff[altIdx][j] =  nvcl_read_u_expgolomb(rdr);
                if (alf_data->alf_chroma_coeff[altIdx][j] > 0) {
                    sign = nvcl_read_bits(rdr, 1);
                    if (sign) alf_data->alf_chroma_coeff[altIdx][j] *= -1;
                }
            }

            if (alf_data->alf_chroma_clip_flag) {
                for (int j = 0; j < 6; j++) {
                    alf_data->alf_chroma_clip_idx[altIdx][j] = nvcl_read_bits(rdr, 2);
                }
            }
        }
    }

    if (alf_data->alf_cc_cb_filter_signal_flag) {
        alf_data->alf_cc_cb_filters_signalled_minus1 = nvcl_read_u_expgolomb(rdr);
        for (int k = 0; k < alf_data->alf_cc_cb_filters_signalled_minus1 + 1; k++) {
            for (int j = 0; j < 7; j++) {
                int code = nvcl_read_bits(rdr, 3);
                if (code) {
                    code = 1 << (code - 1);
                    sign = nvcl_read_bits(rdr, 1);
                    alf_data->alf_cc_mapped_coeff[0][k][j] = sign ? -1*code : code;
                }
                else{
                    alf_data->alf_cc_mapped_coeff[0][k][j] = 0;
                }
            }
        }
    }

    if (alf_data->alf_cc_cr_filter_signal_flag) {
        alf_data->alf_cc_cr_filters_signalled_minus1 = nvcl_read_u_expgolomb(rdr);
        for (int k = 0; k < alf_data->alf_cc_cr_filters_signalled_minus1 + 1; k++) {
            for (int j = 0; j < 7; j++) {
                int code = nvcl_read_bits(rdr, 3);
                if (code) {
                    code = 1 << (code - 1);
                    sign = nvcl_read_bits(rdr, 1);
                    alf_data->alf_cc_mapped_coeff[1][k][j] = sign ? -1*code : code;
                }
                else{
                    alf_data->alf_cc_mapped_coeff[1][k][j] = 0;
                }
            }
        }
    }
}

static int 
nvcl_read_lmcs_data(OVNVCLReader *const rdr, struct OVLMCSData* lmcs,
                    uint8_t aps_chroma_present_flag)
{
    int i;

    lmcs->lmcs_min_bin_idx          = nvcl_read_u_expgolomb(rdr);
    lmcs->lmcs_delta_max_bin_idx    = nvcl_read_u_expgolomb(rdr);
    lmcs->lmcs_delta_cw_prec_minus1 = nvcl_read_u_expgolomb(rdr);

    for (i = lmcs->lmcs_min_bin_idx; i < PIC_CODE_CW_BINS - lmcs->lmcs_delta_max_bin_idx; i++) {
        lmcs->lmcs_delta_abs_cw[i] = nvcl_read_bits(rdr, lmcs->lmcs_delta_cw_prec_minus1 + 1);
        if (lmcs->lmcs_delta_abs_cw[i]) {
            lmcs->lmcs_delta_sign_cw_flag[i] = nvcl_read_flag(rdr);
        }
    }

    if (aps_chroma_present_flag) {
        lmcs->lmcs_delta_abs_crs = nvcl_read_bits(rdr, 3);
        if (lmcs->lmcs_delta_abs_crs) {
            lmcs->lmcs_delta_sign_crs_flag = nvcl_read_flag(rdr);
        }
    }

    return 0;
}

static void
nvcl_read_scaling_list_data(OVNVCLReader *const rdr, struct OVScalingList* sl,
                             uint8_t aps_chroma_present_flag)
{
    int id;
    if (aps_chroma_present_flag) {
        for (id = 0; id < 2; id++) {
            uint8_t is_pred_or_cpy;
            sl->scaling_list_copy_mode_flag[id] = nvcl_read_flag(rdr);
            if (!sl->scaling_list_copy_mode_flag[id]) {
                sl->scaling_list_pred_mode_flag[id] = nvcl_read_flag(rdr);
            }

            is_pred_or_cpy  = sl->scaling_list_pred_mode_flag[id];
            is_pred_or_cpy |= sl->scaling_list_copy_mode_flag[id];

            if (is_pred_or_cpy && id != 0) {
                sl->scaling_list_pred_id_delta[id] = nvcl_read_u_expgolomb(rdr);
            }

            if (!sl->scaling_list_copy_mode_flag[id]) {
                int i;
                for (i = 0; i < 4; i++) {
                    sl->scaling_list_delta_coef[id][i] = nvcl_read_s_expgolomb(rdr);
                }
            }
        }

        for (; id < 8; id++) {
            uint8_t is_pred_or_cpy;
            sl->scaling_list_copy_mode_flag[id] = nvcl_read_flag(rdr);
            if (!sl->scaling_list_copy_mode_flag[id]) {
                sl->scaling_list_pred_mode_flag[id] = nvcl_read_flag(rdr);
            }

            is_pred_or_cpy  = sl->scaling_list_pred_mode_flag[id];
            is_pred_or_cpy |= sl->scaling_list_copy_mode_flag[id];

            if (is_pred_or_cpy && id != 2) {
                sl->scaling_list_pred_id_delta[id] = nvcl_read_u_expgolomb(rdr);
            }

            if (!sl->scaling_list_copy_mode_flag[id]) {
                int i;
                for (i = 0; i < 16; i++) {
                    sl->scaling_list_delta_coef[id][i] = nvcl_read_s_expgolomb(rdr);
                }
            }
        }

        for (; id < 9; id++) {
            uint8_t is_pred_or_cpy;
            sl->scaling_list_copy_mode_flag[id] = nvcl_read_flag(rdr);
            if (!sl->scaling_list_copy_mode_flag[id]) {
                sl->scaling_list_pred_mode_flag[id] = nvcl_read_flag(rdr);
            }

            is_pred_or_cpy  = sl->scaling_list_pred_mode_flag[id];
            is_pred_or_cpy |= sl->scaling_list_copy_mode_flag[id];

            if (!sl->scaling_list_copy_mode_flag[id]) {

                int i;
                for (i = 0; i < 64; i++) {
                    sl->scaling_list_delta_coef[id][i] = nvcl_read_s_expgolomb(rdr);
                }
            }
        }

        for (; id < 14; id++) {
            uint8_t is_pred_or_cpy;
            sl->scaling_list_copy_mode_flag[id] = nvcl_read_flag(rdr);
            if (!sl->scaling_list_copy_mode_flag[id]) {
                sl->scaling_list_pred_mode_flag[id] = nvcl_read_flag(rdr);
            }

            is_pred_or_cpy  = sl->scaling_list_pred_mode_flag[id];
            is_pred_or_cpy |= sl->scaling_list_copy_mode_flag[id];

            if (is_pred_or_cpy) {
                sl->scaling_list_pred_id_delta[id] = nvcl_read_u_expgolomb(rdr);
            }

            if (!sl->scaling_list_copy_mode_flag[id]) {

                int i;
                for (i = 0; i < 64; i++) {
                    sl->scaling_list_delta_coef[id][i] = nvcl_read_s_expgolomb(rdr);
                }
            }
        }

        for (; id < 26; id++) {
            uint8_t is_pred_or_cpy;
            sl->scaling_list_copy_mode_flag[id] = nvcl_read_flag(rdr);
            if (!sl->scaling_list_copy_mode_flag[id]) {
                sl->scaling_list_pred_mode_flag[id] = nvcl_read_flag(rdr);
            }

            is_pred_or_cpy  = sl->scaling_list_pred_mode_flag[id];
            is_pred_or_cpy |= sl->scaling_list_copy_mode_flag[id];

            if (is_pred_or_cpy) {
                sl->scaling_list_pred_id_delta[id] = nvcl_read_u_expgolomb(rdr);
            }

            if (!sl->scaling_list_copy_mode_flag[id]) {

                int i;
                sl->scaling_list_dc_coef[id - 14] = nvcl_read_s_expgolomb(rdr);

                for (i = 0; i < 64; i++) {
                    sl->scaling_list_delta_coef[id][i] = nvcl_read_s_expgolomb(rdr);
                }
            }
        }

        for (; id < 27; id++) {
            uint8_t is_pred_or_cpy;
            sl->scaling_list_copy_mode_flag[id] = nvcl_read_flag(rdr);
            if (!sl->scaling_list_copy_mode_flag[id]) {
                sl->scaling_list_pred_mode_flag[id] = nvcl_read_flag(rdr);
            }

            is_pred_or_cpy  = sl->scaling_list_pred_mode_flag[id];
            is_pred_or_cpy |= sl->scaling_list_copy_mode_flag[id];

            if (is_pred_or_cpy) {
                sl->scaling_list_pred_id_delta[id] = nvcl_read_u_expgolomb(rdr);
            }

            if (!sl->scaling_list_copy_mode_flag[id]) {
                int i;
                sl->scaling_list_dc_coef[id - 14] = nvcl_read_s_expgolomb(rdr);

                for (i = 0; i < 48; i++) {
                    sl->scaling_list_delta_coef[id][i] = nvcl_read_s_expgolomb(rdr);
                }
            }
        }

        for (; id < 28; id++) {
            uint8_t is_pred_or_cpy;
            sl->scaling_list_copy_mode_flag[id] = nvcl_read_flag(rdr);
            if (!sl->scaling_list_copy_mode_flag[id]) {
                sl->scaling_list_pred_mode_flag[id] = nvcl_read_flag(rdr);
            }

            is_pred_or_cpy  = sl->scaling_list_pred_mode_flag[id];
            is_pred_or_cpy |= sl->scaling_list_copy_mode_flag[id];

            if (is_pred_or_cpy) {
                sl->scaling_list_pred_id_delta[id] = nvcl_read_u_expgolomb(rdr);
            }

            if (!sl->scaling_list_copy_mode_flag[id]) {

                int i;
                sl->scaling_list_dc_coef[id - 14] = nvcl_read_s_expgolomb(rdr);

                for (i = 0; i < 48; i++) {
                    sl->scaling_list_delta_coef[id][i] = nvcl_read_s_expgolomb(rdr);
                }
            }
        }
    } else {

        for (id = 2; id < 5; id += 3) {

            sl->scaling_list_copy_mode_flag[id] = nvcl_read_flag(rdr);
            if (!sl->scaling_list_copy_mode_flag[id]) {
                sl->scaling_list_pred_mode_flag[id] = nvcl_read_flag(rdr);
            }

            if (!sl->scaling_list_copy_mode_flag[id]) {
                int i;
                for (i = 0; i < 16; i++) {
                    sl->scaling_list_delta_coef[id][i] = nvcl_read_s_expgolomb(rdr);
                }
            }
        }

        for (id = 5; id < 8; id += 3) {
            uint8_t is_pred_or_cpy;
            sl->scaling_list_copy_mode_flag[id] = nvcl_read_flag(rdr);
            if (!sl->scaling_list_copy_mode_flag[id]) {
                sl->scaling_list_pred_mode_flag[id] = nvcl_read_flag(rdr);
            }

            is_pred_or_cpy  = sl->scaling_list_pred_mode_flag[id];
            is_pred_or_cpy |= sl->scaling_list_copy_mode_flag[id];

            if (is_pred_or_cpy) {
                sl->scaling_list_pred_id_delta[id] = nvcl_read_u_expgolomb(rdr);
            }

            if (!sl->scaling_list_copy_mode_flag[id]) {
                int i;
                for (i = 0; i < 16; i++) {
                    sl->scaling_list_delta_coef[id][i] = nvcl_read_s_expgolomb(rdr);
                }
            }
        }

        for (id = 8; id < 14; id += 3) {
            sl->scaling_list_copy_mode_flag[id] = nvcl_read_flag(rdr);
            if (!sl->scaling_list_copy_mode_flag[id]) {
                sl->scaling_list_pred_mode_flag[id] = nvcl_read_flag(rdr);
            }

            if (!sl->scaling_list_copy_mode_flag[id]) {

                int i;
                for (i = 0; i < 64; i++) {
                    sl->scaling_list_delta_coef[id][i] = nvcl_read_s_expgolomb(rdr);
                }
            }
        }

        for (id = 11; id < 14; id += 3) {
            uint8_t is_pred_or_cpy;
            sl->scaling_list_copy_mode_flag[id] = nvcl_read_flag(rdr);
            if (!sl->scaling_list_copy_mode_flag[id]) {
                sl->scaling_list_pred_mode_flag[id] = nvcl_read_flag(rdr);
            }

            is_pred_or_cpy  = sl->scaling_list_pred_mode_flag[id];
            is_pred_or_cpy |= sl->scaling_list_copy_mode_flag[id];

            if (is_pred_or_cpy) {
                sl->scaling_list_pred_id_delta[id] = nvcl_read_u_expgolomb(rdr);
            }

            if (!sl->scaling_list_copy_mode_flag[id]) {

                int i;
                for (i = 0; i < 64; i++) {
                    sl->scaling_list_delta_coef[id][i] = nvcl_read_s_expgolomb(rdr);
                }
            }
        }

        for (id = 14; id < 26; id += 3) {
            uint8_t is_pred_or_cpy;

            sl->scaling_list_copy_mode_flag[id] = nvcl_read_flag(rdr);

            if (!sl->scaling_list_copy_mode_flag[id]) {
                sl->scaling_list_pred_mode_flag[id] = nvcl_read_flag(rdr);
            }

            is_pred_or_cpy  = sl->scaling_list_pred_mode_flag[id];
            is_pred_or_cpy |= sl->scaling_list_copy_mode_flag[id];

            if (is_pred_or_cpy) {
                sl->scaling_list_pred_id_delta[id] = nvcl_read_u_expgolomb(rdr);
            }

            if (!sl->scaling_list_copy_mode_flag[id]) {
                sl->scaling_list_dc_coef[id - 14] = nvcl_read_s_expgolomb(rdr);

                int i;
                for (i = 0; i < 64; i++) {
                    sl->scaling_list_delta_coef[id][i] = nvcl_read_s_expgolomb(rdr);
                }
            }
        }

        for (id = 26; id < 27; id += 3) {
            uint8_t is_pred_or_cpy;

            sl->scaling_list_copy_mode_flag[id] = nvcl_read_flag(rdr);

            if (!sl->scaling_list_copy_mode_flag[id]) {
                sl->scaling_list_pred_mode_flag[id] = nvcl_read_flag(rdr);
            }

            is_pred_or_cpy  = sl->scaling_list_pred_mode_flag[id];
            is_pred_or_cpy |= sl->scaling_list_copy_mode_flag[id];

            if (is_pred_or_cpy) {
                sl->scaling_list_pred_id_delta[id] = nvcl_read_u_expgolomb(rdr);
            }

            if (!sl->scaling_list_copy_mode_flag[id]) {

                int i;
                sl->scaling_list_dc_coef[id - 14] = nvcl_read_s_expgolomb(rdr);

                for (i = 0; i < 48; i++) {
                    sl->scaling_list_delta_coef[id][i] = nvcl_read_s_expgolomb(rdr);
                }
            }
        }

        for (id = 27; id < 28; id++) {
            uint8_t is_pred_or_cpy;
            sl->scaling_list_copy_mode_flag[id] = nvcl_read_flag(rdr);
            if (!sl->scaling_list_copy_mode_flag[id]) {
                sl->scaling_list_pred_mode_flag[id] = nvcl_read_flag(rdr);
            }

            is_pred_or_cpy  = sl->scaling_list_pred_mode_flag[id];
            is_pred_or_cpy |= sl->scaling_list_copy_mode_flag[id];

            if (is_pred_or_cpy) {
                sl->scaling_list_pred_id_delta[id] = nvcl_read_u_expgolomb(rdr);
            }

            if (!sl->scaling_list_copy_mode_flag[id]) {

                int i;
                sl->scaling_list_dc_coef[id - 14] = nvcl_read_s_expgolomb(rdr);

                for (i = 0; i < 48; i++) {
                    sl->scaling_list_delta_coef[id][i] = nvcl_read_s_expgolomb(rdr);
                }
            }
        }
    }
}

int
nvcl_aps_read(OVNVCLReader *const rdr, OVAPS *const aps,
              OVNVCLCtx *const nvcl_ctx)
{
    aps->aps_params_type                    = nvcl_read_bits(rdr, 3);
    aps->aps_adaptation_parameter_set_id    = nvcl_read_bits(rdr, 5);
    aps->aps_chroma_present_flag            = nvcl_read_flag(rdr);

    if (aps->aps_params_type == APS_ALF) {
        nvcl_read_alf_data(rdr, &aps->aps_alf_data, aps->aps_chroma_present_flag);
    } else if (aps->aps_params_type == APS_LMCS) {
        nvcl_read_lmcs_data(rdr, &aps->aps_lmcs_data, aps->aps_chroma_present_flag);
    } else if (aps->aps_params_type == APS_SCALING_LIST) {
        ov_log(NULL, OVLOG_WARNING, "Ignored unsupported scaling list APS.\n");
        OVScalingList sl = {0};
        nvcl_read_scaling_list_data(rdr, &sl, aps->aps_chroma_present_flag);

    }

    aps->aps_extension_flag = nvcl_read_flag(rdr);
    #if 0
    if(aps->aps_extension_flag) {
        while(more_rbsp_data()) {
            aps->aps_extension_data_flag = nvcl_read_flag(rdr);
        }
    }
    rbsp_trailing_bits()
    #endif
    return 0;
}

int
nvcl_decode_nalu_aps(OVNVCLCtx *const nvcl_ctx, OVNVCLReader *const rdr, uint8_t nalu_type)
{
    int ret;
    /* TODO compare RBSP data to avoid new read */

    OVAPS *aps = ov_mallocz(sizeof(*aps));
    if (!aps) {
        return OVVC_ENOMEM;
    }

    ret = nvcl_aps_read(rdr, aps, nvcl_ctx);
    if (ret < 0) {
        goto cleanup;
    }

    ret = validate_aps(rdr, aps);
    if (ret < 0) {
        goto cleanup;
    }

    uint8_t aps_id = aps->aps_adaptation_parameter_set_id;
    if (aps->aps_params_type == 0) {
        ov_free(nvcl_ctx->alf_aps_list[aps_id]);
        nvcl_ctx->alf_aps_list[aps_id] = aps;
    } else if (aps->aps_params_type == 1) {
        ov_free(nvcl_ctx->lmcs_aps_list[aps_id]);
        nvcl_ctx->lmcs_aps_list[aps_id] = aps;
    } else {
        ov_free(aps);
    }
    return aps_id;

cleanup:
    ov_free(aps);
    return ret;
}

