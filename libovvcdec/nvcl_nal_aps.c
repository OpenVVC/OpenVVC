#include "nvcl.h"
#include "nvcl_utils.h"

typedef struct OVScalingList
{
    uint8_t scaling_list_copy_mode_flag[id];
    uint8_t scaling_list_pred_mode_flag[id];
    uint8_t scaling_list_pred_id_delta[id];
    uint8_t scaling_list_dc_coef[id − 14];
    uint8_t scaling_list_delta_coef[id][i];
} OVScalingList;

typedef struct OVALFData
{
    uint8_t alf_luma_filter_signal_flag;
    uint8_t alf_chroma_filter_signal_flag;
    uint8_t alf_cc_cb_filter_signal_flag;
    uint8_t alf_cc_cr_filter_signal_flag;
    uint8_t alf_luma_clip_flag;
    uint8_t alf_luma_num_filters_signalled_minus1;
    uint8_t alf_luma_coeff_delta_idx[filtIdx];
    uint8_t alf_luma_coeff_abs[sfIdx][j];
    uint8_t alf_luma_coeff_sign[sfIdx][j];
    uint8_t alf_luma_clip_idx[sfIdx][j];
    uint8_t alf_chroma_clip_flag;
    uint8_t alf_chroma_num_alt_filters_minus1;
    uint8_t alf_chroma_coeff_abs[altIdx][j];
    uint8_t alf_chroma_coeff_sign[altIdx][j];
    uint8_t alf_chroma_clip_idx[altIdx][j];
    uint8_t alf_cc_cb_filters_signalled_minus1;
    uint8_t alf_cc_cb_mapped_coeff_abs[k][j];
    uint8_t alf_cc_cb_coeff_sign[k][j];
    uint8_t alf_cc_cr_filters_signalled_minus1;
    uint8_t alf_cc_cr_mapped_coeff_abs[k][j];
    uint8_t alf_cc_cr_coeff_sign[k][j];
} OVALFData;

typedef struct OVLMCS
{
    uint8_t lmcs_min_bin_idx;
    uint8_t lmcs_delta_max_bin_idx;
    uint8_t lmcs_delta_cw_prec_minus1;
    uint8_t lmcs_delta_abs_cw[i];
    uint8_t lmcs_delta_sign_cw_flag[i];
    uint8_t lmcs_delta_abs_crs;
    uint8_t lmcs_delta_sign_crs_flag;
} OVLMCS;

typedef struct OVAPS
{
    uint8_t aps_params_type;
    uint8_t aps_adaptation_parameter_set_id;
    uint8_t aps_chroma_present_flag;
    /* Note unused */
    uint8_t aps_extension_flag;
    uint8_t aps_extension_data_flag;
} OVAPS;

int
nvcl_aps_read(OVNVCLReader *const rdr, OVAPS *const aps,
              OVNVCLCtx *const nvcl_ctx)

#if 0
adaptation_parameter_set_rbsp()
{
    aps->aps_params_type = nvcl_read_bits(rdr, 3);
    aps->aps_adaptation_parameter_set_id = nvcl_read_bits(rdr, 5);
    aps->aps_chroma_present_flag = nvcl_read_flag(rdr);

    if(aps->aps_params_type == ALF_APS) {
        alf_data();
    } else if(aps->aps_params_type == LMCS_APS) {
        lmcs_data();
    } else if(aps->aps_params_type == SCALING_APS) {
        scaling_list_data();
    }

    aps->aps_extension_flag = nvcl_read_flag(rdr);
    if(aps->aps_extension_flag) {
        while(more_rbsp_data()) {
            aps->aps_extension_data_flag = nvcl_read_flag(rdr);
        }
    }
    rbsp_trailing_bits()
}

alf_data()
{
    alf->alf_luma_filter_signal_flag = nvcl_read_flag(rdr);

    if (aps_chroma_present_flag) {
        alf->alf_chroma_filter_signal_flag = nvcl_read_flag(rdr);
        alf->alf_cc_cb_filter_signal_flag = nvcl_read_flag(rdr);
        alf->alf_cc_cr_filter_signal_flag = nvcl_read_flag(rdr);
    }

    if (alf_luma_filter_signal_flag) {
        alf->alf_luma_clip_flag = nvcl_read_flag(rdr);
        alf->alf_luma_num_filters_signalled_minus1 = nvcl_read_u_expgolomb(rdr);
        if (alf_luma_num_filters_signalled_minus1 > 0) {
            for (filtIdx = 0; filtIdx < NumAlfFilters; filtIdx++) {
                alf->alf_luma_coeff_delta_idx[filtIdx];
            }
        }

        for (sfIdx = 0; sfIdx <= alf_luma_num_filters_signalled_minus1; sfIdx++) {
            for (j = 0; j < 12; j++) {
                alf->alf_luma_coeff_abs[sfIdx][j] = nvcl_read_u_expgolomb(rdr);
                if (alf_luma_coeff_abs[sfIdx][j]) {
                    alf->alf_luma_coeff_sign[sfIdx][j] = nvcl_read_bits(rdr, 1);
                }
            }
        }

        if (alf_luma_clip_flag) {
            for (sfIdx = 0; sfIdx <= alf_luma_num_filters_signalled_minus1; sfIdx++) {
                for (j = 0; j < 12; j++) {
                    alf->alf_luma_clip_idx[sfIdx][j] = nvcl_read_bits(rdr, 2);
                }
            }
        }
    }

    if (alf_chroma_filter_signal_flag) {
        alf->alf_chroma_clip_flag = nvcl_read_flag(rdr);
        alf->alf_chroma_num_alt_filters_minus1 = nvcl_read_u_expgolomb(rdr);
        for (altIdx = 0; altIdx <= alf_chroma_num_alt_filters_minus1; altIdx++) {
            for (j = 0; j < 6; j++) {
                alf->alf_chroma_coeff_abs[altIdx][j] = nvcl_read_u_expgolomb(rdr);
                if (alf_chroma_coeff_abs[altIdx][j] > 0) {
                    alf->alf_chroma_coeff_sign[altIdx][j] = nvcl_read_bits(rdr, 1);
                }
            }

            if (alf_chroma_clip_flag) {
                for (j = 0; j < 6; j++) {
                    alf->alf_chroma_clip_idx[altIdx][j] = nvcl_read_bits(rdr, 2);
                }
            }
        }
    }

    if (alf_cc_cb_filter_signal_flag) {
        alf->alf_cc_cb_filters_signalled_minus1 = nvcl_read_u_expgolomb(rdr);
        for (k = 0; k < alf_cc_cb_filters_signalled_minus1 + 1; k++) {
            for (j = 0; j < 7; j++) {
                alf->alf_cc_cb_mapped_coeff_abs[k][j] = nvcl_read_bits(rdr, 3);
                if (alf_cc_cb_mapped_coeff_abs[k][j]) {
                    alf->alf_cc_cb_coeff_sign[k][j] = nvcl_read_bits(rdr, 1);
                }
            }
        }
    }

    if (alf_cc_cr_filter_signal_flag) {
        alf->alf_cc_cr_filters_signalled_minus1 = nvcl_read_u_expgolomb(rdr);
        for (k = 0; k < alf_cc_cr_filters_signalled_minus1 + 1; k++) {
            for (j = 0; j < 7; j++) {
                alf->alf_cc_cr_mapped_coeff_abs[k][j] = nvcl_read_bits(rdr, 3);
                if (alf_cc_cr_mapped_coeff_abs[k][j]) {
                    alf->alf_cc_cr_coeff_sign[k][j] = nvcl_read_bits(rdr, 1);
                }
            }
        }
    }
}

lmcs_data()
{
    lmcs->lmcs_min_bin_idx = nvcl_read_u_expgolomb(rdr);
    lmcs->lmcs_delta_max_bin_idx = nvcl_read_u_expgolomb(rdr);
    lmcs->lmcs_delta_cw_prec_minus1 = nvcl_read_u_expgolomb(rdr);
    for (i = lmcs->lmcs_min_bin_idx; i <= LmcsMaxBinIdx; i++) {
        lmcs->lmcs_delta_abs_cw[i] = nvcl_read_bits(rdr, lmcs->lmcs_delta_cw_prec_minus1 + 1);
        if (lmcs->lmcs_delta_abs_cw[i] > 0) {
            lmcs->lmcs_delta_sign_cw_flag[i] = nvcl_read_flag(rdr);
        }
    }

    if (aps_chroma_present_flag) {
        lmcs->lmcs_delta_abs_crs nvcl_read_bits(rdr, 3);
        if (lmcs->lmcs_delta_abs_crs > 0) {
            lmcs->lmcs_delta_sign_crs_flag = nvcl_read_flag(rdr);
        }
    }
}

#if 0
scaling_list_data()
{
    for (id = 0; id < 28; id ++) {
        int matrixSize = id < 2 ? 2 : ( id < 8 ? 4 : 8);
        if (aps_chroma_present_flag || id % 3 == 2 || id == 27) {
            sl->scaling_list_copy_mode_flag[id] = nvcl_read_flag(rdr);
            if (!sl->scaling_list_copy_mode_flag[id]) {
                sl->scaling_list_pred_mode_flag[id] = nvcl_read_flag(rdr);
            }

            if ((scaling_list_copy_mode_flag[id] || scaling_list_pred_mode_flag[id]) && id != 0 && id != 2 && id != 8) {
                sl->scaling_list_pred_id_delta[id] = nvcl_read_u_expgolomb(rdr);
            }

            if (!scaling_list_copy_mode_flag[id]) {
                int nextCoef = 0;
                if (id > 13) {
                    sl->scaling_list_dc_coef[id − 14] = nvcl_read_s_expgolomb(rdr);
                    nextCoef += sl->scaling_list_dc_coef[id − 14]
                }

                for (i = 0; i < matrixSize * matrixSize; i++) {
                    int x = DiagScanOrder[3][3][i][0];
                    int y = DiagScanOrder[3][3][i][1];
                    if (!( id > 25 && x >= 4 && y >= 4)) {
                        sl->scaling_list_delta_coef[id][i] = nvcl_read_s_expgolomb(rdr);
                        nextCoef += scaling_list_delta_coef[id][i];
                           ScalingList[id][i] = nextCoef;
                    }
                }
            }
        }
    }
}
#endif

scaling_list_data2()
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

            if (!scaling_list_copy_mode_flag[id]) {
                int nextCoef = 0;
                for (i = 0; i < 4; i++) {
                    sl->scaling_list_delta_coef[id][i] = nvcl_read_s_expgolomb(rdr);
                    nextCoef += sl->scaling_list_delta_coef[id][i];
                    ScalingList[id][i] = nextCoef;
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
                int nextCoef = 0;
                for (i = 0; i < 16; i++) {
                    sl->scaling_list_delta_coef[id][i] = nvcl_read_s_expgolomb(rdr);
                    nextCoef += scaling_list_delta_coef[id][i];
                    ScalingList[id][i] = nextCoef;
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

            if (is_pred_or_cpy && id != 8) {
                sl->scaling_list_pred_id_delta[id] = nvcl_read_u_expgolomb(rdr);
            }

            if (!sl->scaling_list_copy_mode_flag[id]) {
                int nextCoef = 0;

                for (i = 0; i < 32; i++) {
                    /* FIXME check position */
                    sl->scaling_list_delta_coef[id][i] = nvcl_read_s_expgolomb(rdr);
                    nextCoef += sl->scaling_list_delta_coef[id][i];
                    ScalingList[id][i] = nextCoef;
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

            if (is_pred_or_cpy && id != 8) {
                sl->scaling_list_pred_id_delta[id] = nvcl_read_u_expgolomb(rdr);
            }

            if (!sl->scaling_list_copy_mode_flag[id]) {
                int nextCoef = 0;

                sl->scaling_list_dc_coef[id − 14] = nvcl_read_s_expgolomb(rdr);
                nextCoef += sl->scaling_list_dc_coef[id − 14]

                for (i = 0; i < 32; i++) {
                    /* FIXME check position */
                    sl->scaling_list_delta_coef[id][i] = nvcl_read_s_expgolomb(rdr);
                    nextCoef += sl->scaling_list_delta_coef[id][i];
                    ScalingList[id][i] = nextCoef;
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

            if (is_pred_or_cpy && id != 8) {
                sl->scaling_list_pred_id_delta[id] = nvcl_read_u_expgolomb(rdr);
            }

            if (!sl->scaling_list_copy_mode_flag[id]) {
                int nextCoef = 0;
                sl->scaling_list_dc_coef[id − 14] = nvcl_read_s_expgolomb(rdr);
                nextCoef += sl->scaling_list_dc_coef[id − 14]

                for (i = 0; i < 32; i++) {
                    /* FIXME check position */
                    int x = DiagScanOrder[3][3][i][0];
                    int y = DiagScanOrder[3][3][i][1];
                    if (!(x >= 4 && y >= 4)) {
                        sl->scaling_list_delta_coef[id][i] = nvcl_read_s_expgolomb(rdr);
                        nextCoef += sl->scaling_list_delta_coef[id][i];
                        ScalingList[id][i] = nextCoef;
                    }
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
                int nextCoef = 0;

                sl->scaling_list_dc_coef[id − 14] = nvcl_read_s_expgolomb(rdr);
                nextCoef += sl->scaling_list_dc_coef[id − 14]

                for (i = 0; i < 32; i++) {
                    /* FIXME check position */
                    int x = DiagScanOrder[3][3][i][0];
                    int y = DiagScanOrder[3][3][i][1];
                    if (!(x >= 4 && y >= 4)) {
                        sl->scaling_list_delta_coef[id][i] = nvcl_read_s_expgolomb(rdr);
                        nextCoef += sl->scaling_list_delta_coef[id][i];
                        ScalingList[id][i] = nextCoef;
                    }
                }
            }
        }
    } else {

        for (id = 2; id < 8; id += 3) {
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
                int nextCoef = 0;
                for (i = 0; i < 16; i++) {
                    sl->scaling_list_delta_coef[id][i] = nvcl_read_s_expgolomb(rdr);
                    nextCoef += scaling_list_delta_coef[id][i];
                    ScalingList[id][i] = nextCoef;
                }
            }
        }

        for (id = 8; id < 14; id += 3) {
            uint8_t is_pred_or_cpy;
            sl->scaling_list_copy_mode_flag[id] = nvcl_read_flag(rdr);
            if (!sl->scaling_list_copy_mode_flag[id]) {
                sl->scaling_list_pred_mode_flag[id] = nvcl_read_flag(rdr);
            }

            is_pred_or_cpy  = sl->scaling_list_pred_mode_flag[id];
            is_pred_or_cpy |= sl->scaling_list_copy_mode_flag[id];

            if (is_pred_or_cpy && id != 8) {
                sl->scaling_list_pred_id_delta[id] = nvcl_read_u_expgolomb(rdr);
            }

            if (!sl->scaling_list_copy_mode_flag[id]) {
                int nextCoef = 0;

                for (i = 0; i < 32; i++) {
                    /* FIXME check position */
                    sl->scaling_list_delta_coef[id][i] = nvcl_read_s_expgolomb(rdr);
                    nextCoef += sl->scaling_list_delta_coef[id][i];
                    ScalingList[id][i] = nextCoef;
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
                int nextCoef = 0;
                sl->scaling_list_dc_coef[id − 14] = nvcl_read_s_expgolomb(rdr);
                nextCoef += sl->scaling_list_dc_coef[id − 14];

                for (i = 0; i < 32; i++) {
                    /* FIXME check position */
                    sl->scaling_list_delta_coef[id][i] = nvcl_read_s_expgolomb(rdr);
                    nextCoef += sl->scaling_list_delta_coef[id][i];
                    ScalingList[id][i] = nextCoef;
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
                int nextCoef = 0;
                sl->scaling_list_dc_coef[id − 14] = nvcl_read_s_expgolomb(rdr);
                nextCoef += sl->scaling_list_dc_coef[id − 14];

                for (i = 0; i < 32; i++) {
                    /* FIXME check position */
                    int x = DiagScanOrder[3][3][i][0];
                    int y = DiagScanOrder[3][3][i][1];
                    if (!(x >= 4 && y >= 4)) {
                        sl->scaling_list_delta_coef[id][i] = nvcl_read_s_expgolomb(rdr);
                        nextCoef += sl->scaling_list_delta_coef[id][i];
                        ScalingList[id][i] = nextCoef;
                    }
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
                int nextCoef = 0;

                sl->scaling_list_dc_coef[id − 14] = nvcl_read_s_expgolomb(rdr);
                nextCoef += sl->scaling_list_dc_coef[id − 14]

                for (i = 0; i < 32; i++) {
                    /* FIXME check position */
                    int x = DiagScanOrder[3][3][i][0];
                    int y = DiagScanOrder[3][3][i][1];
                    if (!(x >= 4 && y >= 4)) {
                        sl->scaling_list_delta_coef[id][i] = nvcl_read_s_expgolomb(rdr);
                        nextCoef += sl->scaling_list_delta_coef[id][i];
                        ScalingList[id][i] = nextCoef;
                    }
                }
            }
        }
    }
}
#endif
