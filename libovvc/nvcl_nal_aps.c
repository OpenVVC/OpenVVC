#include "ovmem.h"

#include "nvcl.h"
#include "nvcl_utils.h"
#include "nvcl_structures.h"

// typedef struct OVScalingList
// {
//     uint8_t scaling_list_copy_mode_flag[id];
//     uint8_t scaling_list_pred_mode_flag[id];
//     uint8_t scaling_list_pred_id_delta[id];
//     uint8_t scaling_list_dc_coef[id − 14];
//     uint8_t scaling_list_delta_coef[id][i];
// } OVScalingList;



static int
validate_aps(OVNVCLReader *rdr, OVAPS *const aps)
{
    /* TODO various check on limitation and max sizes */
    return 1;
}


void 
nvcl_read_alf_data(OVNVCLReader *const rdr, struct OVALFData* alf_data, uint8_t aps_chroma_present_flag)
{
    alf_data->alf_luma_filter_signal_flag = nvcl_read_flag(rdr);
    uint8_t sign = 0;

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
                // alf_data->alf_luma_coeff_delta_idx[filtIdx];
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
nvcl_read_lmcs_data(OVNVCLReader *const rdr, struct OVLMCSData* lmcs, uint8_t aps_chroma_present_flag)
{
    lmcs->lmcs_min_bin_idx = nvcl_read_u_expgolomb(rdr);
    lmcs->lmcs_delta_max_bin_idx = nvcl_read_u_expgolomb(rdr);
    lmcs->lmcs_delta_cw_prec_minus1 = nvcl_read_u_expgolomb(rdr);
    for (int i = lmcs->lmcs_min_bin_idx; i <= PIC_CODE_CW_BINS-(lmcs->lmcs_delta_max_bin_idx + 1); i++) {
        lmcs->lmcs_delta_abs_cw[i] = nvcl_read_bits(rdr, lmcs->lmcs_delta_cw_prec_minus1 + 1);
        if (lmcs->lmcs_delta_abs_cw[i] > 0) {
            lmcs->lmcs_delta_sign_cw_flag[i] = nvcl_read_flag(rdr);
        }
    }

    if (aps_chroma_present_flag) {
        lmcs->lmcs_delta_abs_crs = nvcl_read_bits(rdr, 3);
        if (lmcs->lmcs_delta_abs_crs > 0) {
            lmcs->lmcs_delta_sign_crs_flag = nvcl_read_flag(rdr);
        }
    }
    return 1;
    //ATTENTION: Utilite ?
    //int signCW = code;
    //info.chrResScalingOffset = (1 - 2 * signCW) * absCW;
    //aps->setReshaperAPSInfo(info);
}


void
nvcl_aps_read(OVNVCLReader *const rdr, OVAPS *const aps,
              OVNVCLCtx *const nvcl_ctx)
{
    aps->aps_params_type                    = nvcl_read_bits(rdr, 3);
    aps->aps_adaptation_parameter_set_id    = nvcl_read_bits(rdr, 5);
    aps->aps_chroma_present_flag            = nvcl_read_flag(rdr);

    //TODO: definir les ALF_APS etc
    // if(aps->aps_params_type == ALF_APS) {
    if(aps->aps_params_type == 0) {
        nvcl_read_alf_data(rdr, &aps->aps_alf_data, aps->aps_chroma_present_flag);
    // } else if(aps->aps_params_type == LMCS_APS) {
    } else if(aps->aps_params_type == 1) {
        nvcl_read_lmcs_data(rdr, &aps->aps_lmcs_data, aps->aps_chroma_present_flag);
    // } else if(aps->aps_params_type == SCALING_APS) {
    } else if(aps->aps_params_type == 2) {
        // scaling_list_data();
    }

    aps->aps_extension_flag = nvcl_read_flag(rdr);
    //TODO: prendre en compte les extensions aps.
    #if 0
    if(aps->aps_extension_flag) {
        while(more_rbsp_data()) {
            aps->aps_extension_data_flag = nvcl_read_flag(rdr);
        }
    }
    // rbsp_trailing_bits()
    #endif
}

int
nvcl_decode_nalu_aps(OVNVCLReader *const rdr, OVNVCLCtx *const nvcl_ctx)
{
    int ret;
    /* TODO compare RBSP data to avoid new read */

    OVAPS *aps = ov_mallocz(sizeof(*aps));
    if (!aps) {
        return OV_ENOMEM;
    }

    //TODO: mettre un retour d'erreur.
    // ret = 
    nvcl_aps_read(rdr, aps, nvcl_ctx);
    // if (ret < 0) {
    //     goto cleanup;
    // }

    ret = validate_aps(rdr, aps);
    if (ret < 0) {
        goto cleanup;
    }

    uint8_t aps_id = aps->aps_adaptation_parameter_set_id;
    ov_free(nvcl_ctx->aps_list[aps_id]);
    nvcl_ctx->aps_list[aps_id] = aps;
    // av_log(avctx, AV_LOG_DEBUG, "Success decoding APS:%d \n",aps_id);
    return aps_id;

cleanup:
    ov_free(aps);
    return ret;
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
#if 0
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
