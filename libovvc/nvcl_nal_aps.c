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
    uint8_t scaling_list_pred_id_delta[28];
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


static void
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

struct ScalingList2x2
{
    int16_t coeff[4];
};

struct ScalingList4x4
{
    int16_t coeff[16];
};

struct ScalingList8x8
{
    int16_t coeff[64];
};

struct ScalingLists {
    struct ScalingList2x2 chroma_2x2[2];
    struct ScalingList4x4 ycbcbr_4x4[8 - 2];
    struct ScalingList8x8 ycbcbr_8x8[28 - 8];
};

static void
derive_scaling_matrix_2x2(struct ScalingLists *const sl_ctx, const OVScalingList *sl, int id)
{
    struct ScalingList2x2 *sl_dst = &sl_ctx->chroma_2x2[id];
    const struct ScalingList2x2 *sl_ref = sl_dst;
    int i;
    uint8_t init_val = 0;

    uint8_t is_pred_or_cpy  = sl->scaling_list_pred_mode_flag[id];
            is_pred_or_cpy |= sl->scaling_list_copy_mode_flag[id] << 1;


    if (!is_pred_or_cpy) {
        init_val = 8;
    } else {
        sl_ref -= sl->scaling_list_pred_id_delta[id];
        if (sl_ref == sl_dst) {
            init_val = 16;
        }
    }

    uint8_t coeff = init_val;
    for (i = 0; i < 4; i++) {
        coeff += sl->scaling_list_delta_coef[id][i];
        sl_dst->coeff[i] = sl_ref->coeff[i] + coeff;
    }
}

static void
derive_scaling_matrix_4x4(struct ScalingLists *const sl_ctx, const OVScalingList *sl, int id)
{
    struct ScalingList4x4 *sl_dst = &sl_ctx->ycbcbr_4x4[id - 2];
    const struct ScalingList4x4 *sl_ref = sl_dst;
    int i;
    uint8_t init_val = 0;

    uint8_t is_pred_or_cpy  = sl->scaling_list_pred_mode_flag[id];
    is_pred_or_cpy |= sl->scaling_list_copy_mode_flag[id] << 1;


    if (!is_pred_or_cpy) {
        init_val = 8;
    } else {
        sl_ref -= sl->scaling_list_pred_id_delta[id];
        if (sl_ref == sl_dst) {
            init_val = 16;
        }
    }

    uint8_t coeff = init_val;
    for (i = 0; i < 16; i++) {
        coeff += sl->scaling_list_delta_coef[id][i];
        sl_dst->coeff[i] = sl_ref->coeff[i] + coeff;
    }

}

static void
derive_scaling_matrix_8x8(struct ScalingLists *const sl_ctx, const OVScalingList *sl, int id)
{
    struct ScalingList8x8 *sl_dst = &sl_ctx->ycbcbr_8x8[id - 8];
    const struct ScalingList8x8 *sl_ref = sl_dst;
    int i;
    uint8_t init_val = 0;

    uint8_t is_pred_or_cpy  = sl->scaling_list_pred_mode_flag[id];
    is_pred_or_cpy |= sl->scaling_list_copy_mode_flag[id] << 1;


    if (!is_pred_or_cpy) {
        init_val = 8;
    } else {
        sl_ref -= sl->scaling_list_pred_id_delta[id];
        if (sl_ref == sl_dst) {
            init_val = 16;
        }
    }

    uint8_t coeff = init_val;

    if (id > 13) {
        coeff += sl->scaling_list_dc_coef[id - 14];
    }

    for (i = 0; i < 64; i++) {
        coeff += sl->scaling_list_delta_coef[id][i];
        sl_dst->coeff[i] = sl_ref->coeff[i] + coeff;
    }

}

static uint8_t map[64] =
{
    1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 1, 0,
        1, 1, 1, 1, 1, 0, 0, 1,
        1, 1, 0, 0, 0, 1, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
};

static void
derive_scaling_matrix_8x8_2(struct ScalingLists *const sl_ctx, const OVScalingList *sl, int id)
{
    struct ScalingList8x8 *sl_dst = &sl_ctx->ycbcbr_8x8[id - 8];
    const struct ScalingList8x8 *sl_ref = sl_dst;
    int i;
    uint8_t init_val = 0;

    uint8_t is_pred_or_cpy  = sl->scaling_list_pred_mode_flag[id];
    is_pred_or_cpy |= sl->scaling_list_copy_mode_flag[id] << 1;


    if (!is_pred_or_cpy) {
        init_val = 8;
    } else {
        sl_ref -= sl->scaling_list_pred_id_delta[id];
        if (sl_ref == sl_dst) {
            init_val = 16;
        }
    }

    int16_t coeff = init_val;

    if (id > 13) {
        coeff += sl->scaling_list_dc_coef[id - 14];
    }

    uint8_t skip = 0;
    for (i = 0; i < 64; i++) {
        if (map[i]) {
            coeff += sl->scaling_list_delta_coef[id][i - skip];
        } else {
            skip++;
        }
        sl_dst->coeff[i] = sl_ref->coeff[i] + coeff;
    }

}

#define LOG2_TB_W(id) (((id) >> 3) & 7)
#define LOG2_TB_H(id) ((id) & 7)

#define TB_EXIST(id) (LOG2_TB_H(id) != 7 && LOG2_TB_W(id) != 7 && !(LOG2_TB_H(id) == 0 && LOG2_TB_W(id) == 0))

#define TB_SIZE(id) (TB_EXIST(id) << (OVMIN(LOG2_TB_W(id), 5) + OVMIN(LOG2_TB_H(id), 5)))

static const uint16_t tb_lut_size_tb[64] = {
   TB_SIZE( 0),
   TB_SIZE( 1),
   TB_SIZE( 2),
   TB_SIZE( 3),
   TB_SIZE( 4),
   TB_SIZE( 5),
   TB_SIZE( 6),
   TB_SIZE( 7),
   TB_SIZE( 8),
   TB_SIZE( 9),
   TB_SIZE(10),
   TB_SIZE(11),
   TB_SIZE(12),
   TB_SIZE(13),
   TB_SIZE(14),
   TB_SIZE(15),
   TB_SIZE(16),
   TB_SIZE(17),
   TB_SIZE(18),
   TB_SIZE(19),
   TB_SIZE(20),
   TB_SIZE(21),
   TB_SIZE(22),
   TB_SIZE(23),
   TB_SIZE(24),
   TB_SIZE(25),
   TB_SIZE(26),
   TB_SIZE(27),
   TB_SIZE(28),
   TB_SIZE(29),
   TB_SIZE(30),
   TB_SIZE(31),
   TB_SIZE(32),
   TB_SIZE(33),
   TB_SIZE(34),
   TB_SIZE(35),
   TB_SIZE(36),
   TB_SIZE(37),
   TB_SIZE(38),
   TB_SIZE(39),
   TB_SIZE(40),
   TB_SIZE(41),
   TB_SIZE(42),
   TB_SIZE(43),
   TB_SIZE(44),
   TB_SIZE(45),
   TB_SIZE(46),
   TB_SIZE(47),
   TB_SIZE(48),
   TB_SIZE(49),
   TB_SIZE(50),
   TB_SIZE(51),
   TB_SIZE(52),
   TB_SIZE(53),
   TB_SIZE(54),
   TB_SIZE(55),
   TB_SIZE(56),
   TB_SIZE(57),
   TB_SIZE(58),
   TB_SIZE(59),
   TB_SIZE(60),
   TB_SIZE(61),
   TB_SIZE(62),
   TB_SIZE(63),
};

static const uint16_t tb_lut_offset[64] =
{      0,    0,    2,    6,   14,   30,   62,   94,
      94,   96,  100,  108,  124,  156,  220,  284,
     284,  288,  296,  312,  344,  408,  536,  664,
     664,  672,  688,  720,  784,  912, 1168, 1424,
    1424, 1440, 1472, 1536, 1664, 1920, 2432, 2944,
    2944, 2976, 3040, 3168, 3424, 3936, 4960, 5984,
    5984, 6016, 6080, 6208, 6464, 6976, 8000, 9024,
    9024, 9024, 9024, 9024, 9024, 9024, 9024, 9024
};

const uint8_t diag_scan_2x2[4] = {
    0,  2,
    1,  3
};

const uint8_t diag_scan_4x4[16] = {
    0,  2,  5,  9,
    1,  4,  8, 12,
    3,  7, 11, 14,
    6, 10, 13, 15
};

const uint8_t diag_scan_8x8[64] = {
     0,  2,  5,  9, 14, 20, 27, 35,
     1,  4,  8, 13, 19, 26, 34, 42,
     3,  7, 12, 18, 25, 33, 41, 48,
     6, 11, 17, 24, 32, 40, 47, 53,
    10, 16, 23, 31, 39, 46, 52, 57,
    15, 22, 30, 38, 45, 51, 56, 60,
    21, 29, 37, 44, 50, 55, 59, 62,
    28, 36, 43, 49, 54, 58, 61, 63
};

static void
derive_scaling_tb(const int16_t *mat, int16_t *dst,  uint8_t log2_tb_h, uint8_t log2_tb_w)
{
    const uint8_t log2_tb_s = log2_tb_w + log2_tb_h;
    uint8_t log2_max_wh = OVMAX(log2_tb_w, log2_tb_h);

    uint8_t log2_matrix_s = OVMIN(log2_max_wh, 3);

    uint8_t log2_red_w = OVMIN(log2_tb_w, 5);
    uint8_t log2_red_h = OVMIN(log2_tb_h, 5);

    uint8_t mat_w = 1 << log2_red_w;
    uint8_t mat_h = 1 << log2_red_h;

    int16_t m[32*32] = {0};

    uint16_t cnt = 0;
    const uint8_t *diag_scan = log2_matrix_s < 2 ? diag_scan_2x2 : log2_matrix_s < 3 ? diag_scan_4x4 : diag_scan_8x8;

    for (int i = 0; i < mat_h; ++i) {
        int y = (i << log2_matrix_s) >> log2_tb_h;
        for (int j = 0; j < mat_w; ++j) {
            int x = (j << log2_matrix_s) >> log2_tb_w;
            int16_t val = mat[diag_scan_8x8[x + (y << log2_matrix_s)]];
            m[cnt++] = val;
        }
    }

    if (log2_tb_w >= 2 && log2_tb_h >= 2) {
        /* Reorder per coeff subblock */
        //int16_t m2[32*32] = {0};
        uint16_t cnt = 0;
        for (int i = 0; i < mat_h; i += 4) {
            for (int j = 0; j < mat_w; j += 4) {
                dst[ 0 + cnt * 16] = m[(j + 0) + ((i + 0) << log2_red_w)];
                dst[ 1 + cnt * 16] = m[(j + 1) + ((i + 0) << log2_red_w)];
                dst[ 2 + cnt * 16] = m[(j + 2) + ((i + 0) << log2_red_w)];
                dst[ 3 + cnt * 16] = m[(j + 3) + ((i + 0) << log2_red_w)];

                dst[ 4 + cnt * 16] = m[(j + 0) + ((i + 1) << log2_red_w)];
                dst[ 5 + cnt * 16] = m[(j + 1) + ((i + 1) << log2_red_w)];
                dst[ 6 + cnt * 16] = m[(j + 2) + ((i + 1) << log2_red_w)];
                dst[ 7 + cnt * 16] = m[(j + 3) + ((i + 1) << log2_red_w)];

                dst[ 8 + cnt * 16] = m[(j + 0) + ((i + 2) << log2_red_w)];
                dst[ 9 + cnt * 16] = m[(j + 1) + ((i + 2) << log2_red_w)];
                dst[10 + cnt * 16] = m[(j + 2) + ((i + 2) << log2_red_w)];
                dst[11 + cnt * 16] = m[(j + 3) + ((i + 2) << log2_red_w)];

                dst[12 + cnt * 16] = m[(j + 0) + ((i + 3) << log2_red_w)];
                dst[13 + cnt * 16] = m[(j + 1) + ((i + 3) << log2_red_w)];
                dst[14 + cnt * 16] = m[(j + 2) + ((i + 3) << log2_red_w)];
                dst[15 + cnt * 16] = m[(j + 3) + ((i + 3) << log2_red_w)];
                cnt++;
            }
        }

    } else {
        uint16_t cnt = 0;
        for (int i = 0; i < mat_h; ++i) {
            for (int j = 0; j < mat_w; ++j) {
                dst[cnt] = m[cnt];
                cnt++;
            }
        }
    }

    return;
}


static void
derive_tbs_luma(struct ScalingLists *sl_ctx, int16_t *intra_luts, int16_t *inter_luts)
{
    for (int log2_tb_w = 0; log2_tb_w < 7; ++log2_tb_w) {
        for (int log2_tb_h = 0; log2_tb_h < 7; ++log2_tb_h) {
            uint8_t log2_tb_s = log2_tb_w + log2_tb_h;
            if (log2_tb_s >= 2 && (log2_tb_h > 1 || log2_tb_w > 1)) {
                uint8_t log2_max_wh = OVMAX(log2_tb_w, log2_tb_h);

                uint8_t intra_list_id = (log2_max_wh - 1) * 6 - 4;
                uint8_t inter_list_id = OVMIN(27, intra_list_id + 3);
                uint8_t lut_id = (log2_tb_w << 3) | log2_tb_h;

                if (intra_list_id < 8) {
                    struct ScalingList4x4 *sl_dst = &sl_ctx->ycbcbr_4x4[intra_list_id - 2];
                    derive_scaling_tb(sl_dst->coeff, &intra_luts[tb_lut_offset[lut_id]], log2_tb_w, log2_tb_h);
                    sl_dst = &sl_ctx->ycbcbr_4x4[inter_list_id - 2];
                    derive_scaling_tb(sl_dst->coeff, &inter_luts[tb_lut_offset[lut_id]], log2_tb_w, log2_tb_h);
                } else {
                    struct ScalingList8x8 *sl_dst = &sl_ctx->ycbcbr_8x8[intra_list_id - 8];
                    derive_scaling_tb(sl_dst->coeff, &intra_luts[tb_lut_offset[lut_id]], log2_tb_w, log2_tb_h);
                    sl_dst = &sl_ctx->ycbcbr_8x8[inter_list_id - 8];
                    derive_scaling_tb(sl_dst->coeff, &inter_luts[tb_lut_offset[lut_id]], log2_tb_w, log2_tb_h);

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

        int16_t intra_luts[9024];
        int16_t inter_luts[9024];
        struct ScalingLists sls = {0};
        for (int i = 0; i < 28; ++i) {
            if (i < 2) {
                derive_scaling_matrix_2x2(&sls, &sl, i);
            } else if (i < 8) {
                derive_scaling_matrix_4x4(&sls, &sl, i);
            } else if (i < 26) {
                derive_scaling_matrix_8x8(&sls, &sl, i);
            } else {
                derive_scaling_matrix_8x8_2(&sls, &sl, i);
            }
        }
        derive_tbs_luma(&sls, intra_luts, inter_luts);
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

