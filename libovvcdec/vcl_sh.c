#include "nvcl.h"
#include "nvcl_utils.h"

typedef struct OVSH
{
    uint8_t sh_picture_header_in_slice_header_flag;
    uint8_t sh_subpic_id;
    uint8_t sh_slice_address;
    uint8_t sh_extra_bit[i];
    uint8_t sh_num_tiles_in_slice_minus1;
    uint8_t sh_slice_type;
    uint8_t sh_no_output_of_prior_pics_flag;
    uint8_t sh_alf_enabled_flag;
    uint8_t sh_num_alf_aps_ids_luma;
    uint8_t sh_alf_aps_id_luma[i];
    uint8_t sh_alf_cb_enabled_flag;
    uint8_t sh_alf_cr_enabled_flag;
    uint8_t sh_alf_aps_id_chroma;
    uint8_t sh_alf_cc_cb_enabled_flag;
    uint8_t sh_alf_cc_cb_aps_id;
    uint8_t sh_alf_cc_cr_enabled_flag;
    uint8_t sh_alf_cc_cr_aps_id;
    uint8_t sh_lmcs_used_flag;
    uint8_t sh_explicit_scaling_list_used_flag;
    uint8_t sh_num_ref_idx_active_override_flag;
    uint8_t sh_num_ref_idx_active_minus1[i];
    uint8_t sh_cabac_init_flag;
    uint8_t sh_collocated_from_l0_flag;
    uint8_t sh_collocated_ref_idx;
    uint8_t sh_qp_delta;
    uint8_t sh_cb_qp_offset;
    uint8_t sh_cr_qp_offset;
    uint8_t sh_joint_cbcr_qp_offset;
    uint8_t sh_cu_chroma_qp_offset_enabled_flag;
    uint8_t sh_sao_luma_used_flag;
    uint8_t sh_sao_chroma_used_flag;
    uint8_t sh_deblocking_params_present_flag;
    uint8_t sh_deblocking_filter_disabled_flag;
    uint8_t sh_luma_beta_offset_div2;
    uint8_t sh_luma_tc_offset_div2;
    uint8_t sh_cb_beta_offset_div2;
    uint8_t sh_cb_tc_offset_div2;
    uint8_t sh_cr_beta_offset_div2;
    uint8_t sh_cr_tc_offset_div2;
    uint8_t sh_dep_quant_used_flag;
    uint8_t sh_sign_data_hiding_used_flag;
    uint8_t sh_ts_residual_coding_disabled_flag;
    uint8_t sh_slice_header_extension_length;
    uint8_t sh_slice_header_extension_data_byte[i];
    uint8_t sh_entry_offset_len_minus1;
    uint8_t sh_entry_point_offset_minus1[i];
} OVSH;

int
nvcl_sh_read(OVNVCLReader *const rdr, OVSH *const sh,
             OVNVCLCtx *const nvcl_ctx);
{
    int i;
    OVPH *ph;
    const OVPPS *const pps;
    const OVSPS *const sps;

    sh->sh_picture_header_in_slice_header_flag = nvcl_read_flag(rdr);
    if (sh->sh_picture_header_in_slice_header_flag) {
        #if 0
        picture_header_structure();
        #endif
    }

    if (sps->sps_subpic_info_present_flag) {
        /*FIXME check subpic_id_len_minus1 overrides */
        sh->sh_subpic_id = nvcl_read_bits(rdr, sps->sps_subpic_id_len_minus1 + 1);
    }

    if ((pps->pps_rect_slice_flag && NumSlicesInSubpic[CurrSubpicIdx] > 1) ||
        (!pps->pps_rect_slice_flag && NumTilesInPic > 1)) {
        /*TODO ceil log2_num_tiles_in_pic / num_slices_in_sub_pc*/
        int nb_bits_in_slice_address = 0;
        sh->sh_slice_address = nvcl_read_bits(rdr, nb_bits_in_slice_address);
    }

    for (i = 0; i < NumExtraShBits; i++) {
        sh->sh_extra_bit[i] = nvcl_read_bits(rdr, 1);
    }

    if (!pps->pps_rect_slice_flag && NumTilesInPic − sh->sh_slice_address > 1) {
        sh->sh_num_tiles_in_slice_minus1 = nvcl_read_u_expgolomb(rdr);
    }

    if (ph->ph_inter_slice_allowed_flag) {
        sh->sh_slice_type = nvcl_read_u_expgolomb(rdr);
    }

    if (nal_unit_type == IDR_W_RADL ||
        nal_unit_type == IDR_N_LP   ||
        nal_unit_type == CRA_NUT    ||
        nal_unit_type == GDR_NUT) {
        sh->sh_no_output_of_prior_pics_flag = nvcl_read_flag(rdr);
    }

    if (sps->sps_alf_enabled_flag && !pps->pps_alf_info_in_ph_flag) {

        sh->sh_alf_enabled_flag = nvcl_read_flag(rdr);
        if (sh->sh_alf_enabled_flag) {

            sh->sh_num_alf_aps_ids_luma;
            for (i = 0; i < sh->sh_num_alf_aps_ids_luma; i++) {
                sh->sh_alf_aps_id_luma[i] = nvcl_read_bits(rdr, 3);
            }

            if (sps->sps_chroma_format_idc != 0) {
                sh->sh_alf_cb_enabled_flag = nvcl_read_flag(rdr);
                sh->sh_alf_cr_enabled_flag = nvcl_read_flag(rdr);
            }

            if (sh->sh_alf_cb_enabled_flag || sh->sh_alf_cr_enabled_flag) {
                sh->sh_alf_aps_id_chroma = nvcl_read_bits(rdr, 3);
            }

            if (sps->sps_ccalf_enabled_flag) {

                sh->sh_alf_cc_cb_enabled_flag = nvcl_read_flag(rdr);
                if (sh->sh_alf_cc_cb_enabled_flag) {
                    sh->sh_alf_cc_cb_aps_id = nvcl_read_bits(rdr, 3);
                }

                sh->sh_alf_cc_cr_enabled_flag = nvcl_read_flag(rdr);
                if (sh->sh_alf_cc_cr_enabled_flag) {
                    sh->sh_alf_cc_cr_aps_id = nvcl_read_bits(rdr, 3);
                }
            }
        }
    }

    if (ph->ph_lmcs_enabled_flag && !sh->sh_picture_header_in_slice_header_flag) {
        sh->sh_lmcs_used_flag = nvcl_read_flag(rdr);
    }

    if (ph->ph_explicit_scaling_list_enabled_flag && !sh->sh_picture_header_in_slice_header_flag) {
        sh->sh_explicit_scaling_list_used_flag = nvcl_read_flag(rdr);
    }

    if (!pps->pps_rpl_info_in_ph_flag && ((nal_unit_type != IDR_W_RADL && nal_unit_type != IDR_N_LP) || sps->sps_idr_rpl_present_flag)) {
        #if 0
        ref_pic_lists();
        #endif
    }

    if ((sh->sh_slice_type != I && num_ref_entries[0][RplsIdx[0]] > 1) || (sh->sh_slice_type == B && num_ref_entries[1][RplsIdx[1]] > 1)) {
        sh->sh_num_ref_idx_active_override_flag = nvcl_read_flag(rdr);
        if (sh->sh_num_ref_idx_active_override_flag) {
            for (i = 0; i < (sh->sh_slice_type == B ? 2: 1); i++) {
                if (num_ref_entries[i][RplsIdx[i]] > 1) {
                    sh->sh_num_ref_idx_active_minus1[i] = nvcl_read_u_expgolomb(rdr);
                }
            }
        }
    }

    if (sh->sh_slice_type != I) {
        if (pps->pps_cabac_init_present_flag) {
            sh->sh_cabac_init_flag = nvcl_read_flag(rdr);
        }

        if (ph->ph_temporal_mvp_enabled_flag && !pps->pps_rpl_info_in_ph_flag) {
            if (sh->sh_slice_type == B) {
                sh->sh_collocated_from_l0_flag = nvcl_read_flag(rdr);
            }

            if ((sh->sh_collocated_from_l0_flag && NumRefIdxActive[0] > 1) || (! sh->sh_collocated_from_l0_flag && NumRefIdxActive[1] > 1)) {
                sh->sh_collocated_ref_idx = nvcl_read_u_expgolomb(rdr);
            }
        }

        if (!pps->pps_wp_info_in_ph_flag && ((pps->pps_weighted_pred_flag && sh->sh_slice_type == P) || (pps->pps_weighted_bipred_flag && sh->sh_slice_type == B))) {
            pred_weight_table();
        }
    }

    if (!pps->pps_qp_delta_info_in_ph_flag) {
        sh->sh_qp_delta = nvcl_read_s_expgolomb(rdr);
    }

    if (pps->pps_slice_chroma_qp_offsets_present_flag) {
        sh->sh_cb_qp_offset = nvcl_read_s_expgolomb(rdr);
        sh->sh_cr_qp_offset = nvcl_read_s_expgolomb(rdr);
        if (sps->sps_joint_cbcr_enabled_flag) {
            sh->sh_joint_cbcr_qp_offset = nvcl_read_s_expgolomb(rdr);
        }
    }

    if (pps->pps_cu_chroma_qp_offset_list_enabled_flag) {
        sh->sh_cu_chroma_qp_offset_enabled_flag = nvcl_read_flag(rdr);
    }

    if (sps->sps_sao_enabled_flag && !pps->pps_sao_info_in_ph_flag) {
        sh->sh_sao_luma_used_flag = nvcl_read_flag(rdr);
        if (sps->sps_chroma_format_idc != 0) {
            sh->sh_sao_chroma_used_flag = nvcl_read_flag(rdr);
        }
    }

    if (pps->pps_deblocking_filter_override_enabled_flag && !pps->pps_dbf_info_in_ph_flag) {
        sh->sh_deblocking_params_present_flag = nvcl_read_flag(rdr);
    }

    if (sh->sh_deblocking_params_present_flag) {
        if (!pps->pps_deblocking_filter_disabled_flag) {
            sh->sh_deblocking_filter_disabled_flag = nvcl_read_flag(rdr);
        }

        if (!sh->sh_deblocking_filter_disabled_flag) {
            sh->sh_luma_beta_offset_div2 = nvcl_read_s_expgolomb(rdr);
            sh->sh_luma_tc_offset_div2   = nvcl_read_s_expgolomb(rdr);
            if (pps->pps_chroma_tool_offsets_present_flag) {
                sh->sh_cb_beta_offset_div2 = nvcl_read_s_expgolomb(rdr);
                sh->sh_cb_tc_offset_div2   = nvcl_read_s_expgolomb(rdr);
                sh->sh_cr_beta_offset_div2 = nvcl_read_s_expgolomb(rdr);
                sh->sh_cr_tc_offset_div2   = nvcl_read_s_expgolomb(rdr);
            }
        }
    }


    if (sps->sps_dep_quant_enabled_flag) {
        sh->sh_dep_quant_used_flag = nvcl_read_flag(rdr);
    }

    if (sps->sps_sign_data_hiding_enabled_flag && !sh->sh_dep_quant_used_flag) {
        sh->sh_sign_data_hiding_used_flag = nvcl_read_flag(rdr);
    }

    if (sps->sps_transform_skip_enabled_flag && !sh->sh_dep_quant_used_flag && !sh->sh_sign_data_hiding_used_flag) {
        sh->sh_ts_residual_coding_disabled_flag = nvcl_read_flag(rdr);
    }

    if (pps->pps_slice_header_extension_present_flag) {
        sh->sh_slice_header_extension_length = nvcl_read_u_expgolomb(rdr);
        for (i = 0; i < sh->sh_slice_header_extension_length; i++) {
            sh->sh_slice_header_extension_data_byte[i] = nvcl_read_bits(rdr, 8);
        }
    }

    /*FIXME derive nb entry points */
    if (NumEntryPoints > 0) {
        sh->sh_entry_offset_len_minus1 = nvcl_read_u_expgolomb(rdr);
        for (i = 0; i < NumEntryPoints; i++) {
            sh->sh_entry_point_offset_minus1[i] = nvcl_read_bits(rdr, sh->sh_entry_offset_len_minus1 + 1);
        }
    }

    byte_alignment()
}

