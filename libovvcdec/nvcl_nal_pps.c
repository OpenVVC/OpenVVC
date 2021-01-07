#include "nvcl.h"
#include "nvcl_utils.h"

typedef struct OVPPS
{
    uint8_t pps_pic_parameter_set_id;
    uint8_t pps_seq_parameter_set_id;
    uint8_t pps_mixed_nalu_types_in_pic_flag;
    uint8_t pps_pic_width_in_luma_samples;
    uint8_t pps_pic_height_in_luma_samples;
    uint8_t pps_conformance_window_flag;
    uint8_t pps_conf_win_left_offset;
    uint8_t pps_conf_win_right_offset;
    uint8_t pps_conf_win_top_offset;
    uint8_t pps_conf_win_bottom_offset;

    uint8_t pps_scaling_window_explicit_signalling_flag;
    uint8_t pps_scaling_win_left_offset;
    uint8_t pps_scaling_win_right_offset;
    uint8_t pps_scaling_win_top_offset;
    uint8_t pps_scaling_win_bottom_offset;

    uint8_t pps_output_flag_present_flag;
    uint8_t pps_no_pic_partition_flag;
    uint8_t pps_subpic_id_mapping_present_flag;
    uint8_t pps_num_subpics_minus1;

    uint8_t pps_subpic_id_len_minus1;
    uint8_t pps_subpic_id[i];
    uint8_t pps_log2_ctu_size_minus5;
    uint8_t pps_num_exp_tile_columns_minus1;
    uint8_t pps_num_exp_tile_rows_minus1;
    uint8_t pps_tile_column_width_minus1[i];

    uint8_t pps_tile_row_height_minus1[i];

    uint8_t pps_loop_filter_across_tiles_enabled_flag;
    uint8_t pps_rect_slice_flag;

    uint8_t pps_single_slice_per_subpic_flag;

    uint8_t pps_num_slices_in_pic_minus1;
    uint8_t pps_tile_idx_delta_present_flag;

    uint8_t pps_slice_width_in_tiles_minus1[i];

    uint8_t pps_slice_height_in_tiles_minus1[i];

    uint8_t pps_num_exp_slices_in_tile[i];
    uint8_t pps_exp_slice_height_in_ctus_minus1[i][j];

    uint8_t pps_tile_idx_delta_val[i];

    uint8_t pps_loop_filter_across_slices_enabled_flag;

    uint8_t pps_cabac_init_present_flag;
    uint8_t pps_num_ref_idx_default_active_minus1[i];

    uint8_t pps_rpl1_idx_present_flag;
    uint8_t pps_weighted_pred_flag;
    uint8_t pps_weighted_bipred_flag;
    uint8_t pps_ref_wraparound_enabled_flag;
    uint8_t pps_pic_width_minus_wraparound_offset;

    uint8_t pps_init_qp_minus26;
    uint8_t pps_cu_qp_delta_enabled_flag;
    uint8_t pps_chroma_tool_offsets_present_flag;
    uint8_t pps_cb_qp_offset;
    uint8_t pps_cr_qp_offset;
    uint8_t pps_joint_cbcr_qp_offset_present_flag;
    uint8_t pps_joint_cbcr_qp_offset_value;

    uint8_t pps_slice_chroma_qp_offsets_present_flag;
    uint8_t pps_cu_chroma_qp_offset_list_enabled_flag;
    uint8_t pps_chroma_qp_offset_list_len_minus1;
    uint8_t pps_cb_qp_offset_list[i];
    uint8_t pps_cr_qp_offset_list[i];
    uint8_t pps_joint_cbcr_qp_offset_list[i];

    uint8_t pps_deblocking_filter_control_present_flag;
    uint8_t pps_deblocking_filter_override_enabled_flag;
    uint8_t pps_deblocking_filter_disabled_flag;
    uint8_t pps_dbf_info_in_ph_flag;

    uint8_t pps_luma_beta_offset_div2;
    uint8_t pps_luma_tc_offset_div2;
    uint8_t pps_cb_beta_offset_div2;
    uint8_t pps_cb_tc_offset_div2;
    uint8_t pps_cr_beta_offset_div2;
    uint8_t pps_cr_tc_offset_div2;

    uint8_t pps_rpl_info_in_ph_flag;
    uint8_t pps_sao_info_in_ph_flag;
    uint8_t pps_alf_info_in_ph_flag;
    uint8_t pps_wp_info_in_ph_flag;

    uint8_t pps_qp_delta_info_in_ph_flag;

    uint8_t pps_picture_header_extension_present_flag;
    uint8_t pps_slice_header_extension_present_flag;
    uint8_t pps_extension_flag;

    uint8_t pps_extension_data_flag;

} OVPPS;


int
nvcl_pps_read(OVNVCLReader *const rdr, OVPPS *const pps,
              OVNVCLCtx *const nvcl_ctx)
{
    int i, j;

    pps->pps_pic_parameter_set_id = nvcl_read_bits(rdr, 6);
    pps->pps_seq_parameter_set_id = nvcl_read_bits(rdr, 4);

    pps->pps_mixed_nalu_types_in_pic_flag = nvcl_read_flag(rdr);

    pps->pps_pic_width_in_luma_samples  = nvcl_read_u_expgolomb(rdr);
    pps->pps_pic_height_in_luma_samples = nvcl_read_u_expgolomb(rdr);

    pps->pps_conformance_window_flag = nvcl_read_flag(rdr);
    if (pps->pps_conformance_window_flag) {
        pps->pps_conf_win_left_offset   = nvcl_read_u_expgolomb(rdr);
        pps->pps_conf_win_right_offset  = nvcl_read_u_expgolomb(rdr);
        pps->pps_conf_win_top_offset    = nvcl_read_u_expgolomb(rdr);
        pps->pps_conf_win_bottom_offset = nvcl_read_u_expgolomb(rdr);
    }

    pps->pps_scaling_window_explicit_signalling_flag = nvcl_read_flag(rdr);
    if (pps->pps_scaling_window_explicit_signalling_flag) {
        pps->pps_scaling_win_left_offset   = nvcl_read_s_expgolomb(rdr);
        pps->pps_scaling_win_right_offset  = nvcl_read_s_expgolomb(rdr);
        pps->pps_scaling_win_top_offset    = nvcl_read_s_expgolomb(rdr);
        pps->pps_scaling_win_bottom_offset = nvcl_read_s_expgolomb(rdr);
    }

    pps->pps_output_flag_present_flag = nvcl_read_flag(rdr);

    pps->pps_no_pic_partition_flag = nvcl_read_flag(rdr);

    pps->pps_subpic_id_mapping_present_flag = nvcl_read_flag(rdr);
    if (pps->pps_subpic_id_mapping_present_flag) {
        if (!pps->pps_no_pic_partition_flag) {
            pps->pps_num_subpics_minus1 = nvcl_read_u_expgolomb(rdr);
        }

        pps->pps_subpic_id_len_minus1 = nvcl_read_u_expgolomb(rdr);
        for(i = 0; i <= pps->pps_num_subpics_minus1; i++) {
            pps->pps_subpic_id[i] = nvcl_read_bits(rdr, pps->pps_subpic_id_len_minus1 + 1);
        }
    }

    if (!pps->pps_no_pic_partition_flag) {

        pps->pps_log2_ctu_size_minus5 = nvcl_read_bits(rdr, 2);

        pps->pps_num_exp_tile_columns_minus1 = nvcl_read_u_expgolomb(rdr);
        pps->pps_num_exp_tile_rows_minus1    = nvcl_read_u_expgolomb(rdr);

        for(i = 0; i <= pps->pps_num_exp_tile_columns_minus1; i++) {
            pps->pps_tile_column_width_minus1[i] = nvcl_read_u_expgolomb(rdr);
        }

        for(i = 0; i <= pps->pps_num_exp_tile_rows_minus1; i++) {
            pps->pps_tile_row_height_minus1[i] = nvcl_read_u_expgolomb(rdr);
        }

        if (NumTilesInPic > 1) {
            pps->pps_loop_filter_across_tiles_enabled_flag = nvcl_read_flag(rdr);
            pps->pps_rect_slice_flag = nvcl_read_flag(rdr);
        }

        if (pps->pps_rect_slice_flag) {
            pps->pps_single_slice_per_subpic_flag = nvcl_read_flag(rdr);
        }

        if (pps->pps_rect_slice_flag && !pps->pps_single_slice_per_subpic_flag) {
            pps->pps_num_slices_in_pic_minus1 = nvcl_read_u_expgolomb(rdr);
            if (pps->pps_num_slices_in_pic_minus1 > 1) {
                pps->pps_tile_idx_delta_present_flag = nvcl_read_flag(rdr);
            }

            for(i = 0; i < pps->pps_num_slices_in_pic_minus1; i++) {
                if (SliceTopLeftTileIdx[i] % NumTileColumns != NumTileColumns − 1) {
                    pps->pps_slice_width_in_tiles_minus1[i] = nvcl_read_u_expgolomb(rdr);
                }

                if (SliceTopLeftTileIdx[i] / NumTileColumns != NumTileRows − 1 && (pps->pps_tile_idx_delta_present_flag || SliceTopLeftTileIdx[i] % NumTileColumns == 0)) {
                    pps->pps_slice_height_in_tiles_minus1[i] = nvcl_read_u_expgolomb(rdr);
                }

                if (pps->pps_slice_width_in_tiles_minus1[i] == 0 && pps->pps_slice_height_in_tiles_minus1[i] == 0 && RowHeightVal[SliceTopLeftTileIdx[i] / NumTileColumns] > 1) {
                    pps->pps_num_exp_slices_in_tile[i] = nvcl_read_u_expgolomb(rdr);
                    for(j = 0; j < pps->pps_num_exp_slices_in_tile[i]; j++) {
                        pps->pps_exp_slice_height_in_ctus_minus1[i][j] = nvcl_read_u_expgolomb(rdr);
                    }
                    i += NumSlicesInTile[i] − 1;
                }

                if (pps->pps_tile_idx_delta_present_flag && i < pps->pps_num_slices_in_pic_minus1) {
                    pps->pps_tile_idx_delta_val[i] = nvcl_read_s_expgolomb(rdr);
                }
            }
        }

        if (!pps->pps_rect_slice_flag || pps->pps_single_slice_per_subpic_flag || pps->pps_num_slices_in_pic_minus1 > 0) {
            pps->pps_loop_filter_across_slices_enabled_flag = nvcl_read_flag(rdr);
        }
    }

    pps->pps_cabac_init_present_flag = nvcl_read_flag(rdr);
    for(i = 0; i < 2; i++) {
        pps->pps_num_ref_idx_default_active_minus1[i] = nvcl_read_u_expgolomb(rdr);
    }

    pps->pps_rpl1_idx_present_flag = nvcl_read_flag(rdr);

    pps->pps_weighted_pred_flag   = nvcl_read_flag(rdr);
    pps->pps_weighted_bipred_flag = nvcl_read_flag(rdr);

    pps->pps_ref_wraparound_enabled_flag = nvcl_read_flag(rdr);
    if (pps->pps_ref_wraparound_enabled_flag) {
        pps->pps_pic_width_minus_wraparound_offset = nvcl_read_u_expgolomb(rdr);
    }

    pps->pps_init_qp_minus26 = nvcl_read_s_expgolomb(rdr);

    pps->pps_cu_qp_delta_enabled_flag = nvcl_read_flag(rdr);

    pps->pps_chroma_tool_offsets_present_flag = nvcl_read_flag(rdr);
    if (pps->pps_chroma_tool_offsets_present_flag) {
        pps->pps_cb_qp_offset = nvcl_read_s_expgolomb(rdr);
        pps->pps_cr_qp_offset = nvcl_read_s_expgolomb(rdr);
        pps->pps_joint_cbcr_qp_offset_present_flag = nvcl_read_flag(rdr);
        if (pps->pps_joint_cbcr_qp_offset_present_flag) {
            pps->pps_joint_cbcr_qp_offset_value = nvcl_read_s_expgolomb(rdr);
        }

        pps->pps_slice_chroma_qp_offsets_present_flag = nvcl_read_flag(rdr);
        pps->pps_cu_chroma_qp_offset_list_enabled_flag = nvcl_read_flag(rdr);
        if (pps->pps_cu_chroma_qp_offset_list_enabled_flag) {
            pps->pps_chroma_qp_offset_list_len_minus1 = nvcl_read_u_expgolomb(rdr);
            for (i = 0; i <= pps->pps_chroma_qp_offset_list_len_minus1; i++) {
                pps->pps_cb_qp_offset_list[i] = nvcl_read_s_expgolomb(rdr);
                pps->pps_cr_qp_offset_list[i] = nvcl_read_s_expgolomb(rdr);
                if (pps->pps_joint_cbcr_qp_offset_present_flag) {
                    pps->pps_joint_cbcr_qp_offset_list[i] = nvcl_read_s_expgolomb(rdr);
                }
            }
        }
    }

    pps->pps_deblocking_filter_control_present_flag = nvcl_read_flag(rdr);
    if (pps->pps_deblocking_filter_control_present_flag) {
        pps->pps_deblocking_filter_override_enabled_flag = nvcl_read_flag(rdr);
        pps->pps_deblocking_filter_disabled_flag = nvcl_read_flag(rdr);
        if (!pps->pps_no_pic_partition_flag && pps->pps_deblocking_filter_override_enabled_flag) {
            pps->pps_dbf_info_in_ph_flag = nvcl_read_flag(rdr);
        }

        if (!pps->pps_deblocking_filter_disabled_flag) {
            pps->pps_luma_beta_offset_div2 = nvcl_read_s_expgolomb(rdr);
            pps->pps_luma_tc_offset_div2   = nvcl_read_s_expgolomb(rdr);
            if (pps->pps_chroma_tool_offsets_present_flag) {
                pps->pps_cb_beta_offset_div2 = nvcl_read_s_expgolomb(rdr);
                pps->pps_cb_tc_offset_div2   = nvcl_read_s_expgolomb(rdr);
                pps->pps_cr_beta_offset_div2 = nvcl_read_s_expgolomb(rdr);
                pps->pps_cr_tc_offset_div2   = nvcl_read_s_expgolomb(rdr);
            }
        }
    }

    if (!pps->pps_no_pic_partition_flag) {
        pps->pps_rpl_info_in_ph_flag = nvcl_read_flag(rdr);
        pps->pps_sao_info_in_ph_flag = nvcl_read_flag(rdr);
        pps->pps_alf_info_in_ph_flag = nvcl_read_flag(rdr);
        if ((pps->pps_weighted_pred_flag || pps->pps_weighted_bipred_flag) && pps->pps_rpl_info_in_ph_flag) {
            pps->pps_wp_info_in_ph_flag = nvcl_read_flag(rdr);
        }

        pps->pps_qp_delta_info_in_ph_flag = nvcl_read_flag(rdr);
    }

    pps->pps_picture_header_extension_present_flag = nvcl_read_flag(rdr);
    pps->pps_slice_header_extension_present_flag = nvcl_read_flag(rdr);
    pps->pps_extension_flag = nvcl_read_flag(rdr);

    if (pps->pps_extension_flag) {
        while (more_rbsp_data()) {
            pps->pps_extension_data_flag = nvcl_read_flag(rdr);
        }
    }

    rbsp_trailing_bits()
}
