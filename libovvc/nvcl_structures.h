#ifndef NVCL_STRUCTURES_H
#define NVCL_STRUCTURES_H

#include <stdint.h>
#include "nvcl_utils.h"
#include "rcn_alf.h"

#define OV_MAX_NB_RP 16
#define PIC_CODE_CW_BINS 16

/* FIXME :
 *     -use union to shorten size and better reflect
 *    how Ref Picture are read ?
 */
/* Syntax Elements in Reference Picture
 */
struct RefPic
{
    /* Flags to determine Reference Picture Type */
    uint8_t inter_layer_ref_pic_flag;
    uint8_t st_ref_pic_flag;

    /* Short Term References Pictures */
    uint8_t abs_delta_poc_st;
    uint8_t strp_entry_sign_flag;

    /* Long Term References Pictures */
    uint16_t rpls_poc_lsb_lt;/*fixme*/

    /* Inter Layer References Pictures */
    uint8_t ilrp_idx;
};

struct OVRPL
{
    uint8_t num_ref_entries;
    uint8_t num_ref_active_entries;
    uint8_t ltrp_in_header_flag;
    struct RefPic rp_list[OV_MAX_NB_RP];
};

/* Syntax Elements specific to Headers (PH or SH)
 */
struct RPLHeader
{
    /* Flag to tell if Ref Picture List is derived from
     * SPS
     */
    uint8_t rpl_sps_flag;

    /* SPS Ref Pic List Index if rpl_sps_flag is ON
     */
    uint8_t rpl_idx;

    /* Ref Pic List to be read in Header RBSP
     */
    struct OVRPL rpl_data;

    struct LTInfo {
       /* Note this is only read in headers */
       uint8_t delta_poc_msb_cycle_present_flag;
       uint16_t delta_poc_msb_cycle_lt;
       uint16_t poc_lsb_lt;
    } lt_info[OV_MAX_NB_RP];
};

struct OVHRPL
{
    struct RPLHeader rpl_h0;
    struct RPLHeader rpl_h1;
    const struct OVRPL *rpl0;
    const struct OVRPL *rpl1;
};

/* FIXME find a better place for this structure */
struct OVDPBParams
{
    uint8_t dpb_max_dec_pic_buffering_minus1;
    uint8_t dpb_max_num_reorder_pics;
    uint8_t dpb_max_latency_increase_plus1;
};


struct OVSEIFGrain
{
    uint8_t fg_characteristics_cancel_flag;
    uint8_t fg_model_id;
    uint8_t fg_separate_colour_description_present_flag;
    uint8_t fg_bit_depth_luma_minus8;
    uint8_t fg_bit_depth_chroma_minus8;
    uint8_t fg_full_range_flag;
    uint8_t fg_colour_primaries;
    uint8_t fg_transfer_characteristics;
    uint8_t fg_matrix_coeffs;
    uint8_t fg_blending_mode_id;
    uint8_t fg_log2_scale_factor;

    uint8_t fg_comp_model_present_flag[3];
    uint8_t fg_num_intensity_intervals_minus1[3];
    uint8_t fg_num_model_values_minus1[3];

    //TODO: maybe more than 8 intensity intervals?
    uint8_t fg_intensity_interval_lower_bound[3][8];
    uint8_t fg_intensity_interval_upper_bound[3][8];
    int16_t fg_comp_model_value[3][8][3];

    uint8_t fg_characteristics_persistence_flag;
    int16_t fg_idr_pic;
};

struct OVSEISLHDR
{
    void* slhdr_context;
    uint8_t payload_array[255];
};


/*Structure containing all SEI structures*/
struct OVSEI
{
    struct OVSEIFGrain* sei_fg;
    struct OVSEISLHDR* sei_slhdr;
};

struct OVSPS
{
    uint8_t sps_seq_parameter_set_id;
    uint8_t sps_video_parameter_set_id;
    uint8_t sps_max_sublayers_minus1;
    uint8_t sps_chroma_format_idc;
    uint8_t sps_log2_ctu_size_minus5;

    uint8_t sps_ptl_dpb_hrd_params_present_flag;

    struct OVPTL *profile_tier_level;

    uint8_t sps_gdr_enabled_flag;

    uint8_t sps_ref_pic_resampling_enabled_flag;
    uint8_t sps_res_change_in_clvs_allowed_flag;

    uint16_t sps_pic_width_max_in_luma_samples;
    uint16_t sps_pic_height_max_in_luma_samples;

    uint8_t sps_conformance_window_flag;
    uint8_t sps_conf_win_left_offset;
    uint8_t sps_conf_win_right_offset;
    uint8_t sps_conf_win_top_offset;
    uint8_t sps_conf_win_bottom_offset;

    uint8_t sps_subpic_info_present_flag;
    uint8_t sps_num_subpics_minus1;
    uint8_t sps_independent_subpics_flag;
    uint8_t sps_subpic_same_size_flag;
    uint8_t sps_subpic_ctu_top_left_x[16]; /* max num_sub_pic_h ?*/
    uint8_t sps_subpic_ctu_top_left_y[16]; /* max num_sub_pic_v? */
    uint8_t sps_subpic_width_minus1[16]; /* max num_sub_pic */
    uint8_t sps_subpic_height_minus1[16]; /* max num_sub_pic */
    uint8_t sps_subpic_treated_as_pic_flag[16]; /* max num_sub_pic */
    uint8_t sps_loop_filter_across_subpic_enabled_flag[16]; /* max num_sub_pic */
    uint8_t sps_subpic_id_len_minus1;
    uint8_t sps_subpic_id_mapping_explicitly_signalled_flag;
    uint8_t sps_subpic_id_mapping_present_flag;
    uint8_t sps_subpic_id[16]; /* max num_sub_pic */

    uint8_t sps_bitdepth_minus8;
    uint8_t sps_entropy_coding_sync_enabled_flag;
    uint8_t sps_entry_point_offsets_present_flag;
    uint8_t sps_log2_max_pic_order_cnt_lsb_minus4;

    uint8_t sps_poc_msb_cycle_flag;
    uint8_t sps_poc_msb_cycle_len_minus1;

    uint8_t sps_num_extra_ph_bytes;
    uint8_t sps_extra_ph_bit_present_flag[24];
    uint8_t sps_num_extra_sh_bytes;
    uint8_t sps_extra_sh_bit_present_flag[24]; /*max num_extra_ph_bytes?*/

    uint8_t sps_sublayer_dpb_params_flag;
    struct OVDPBParams dpb_parameters[64];

    uint8_t sps_log2_min_luma_coding_block_size_minus2;
    uint8_t sps_partition_constraints_override_enabled_flag;
    uint8_t sps_log2_diff_min_qt_min_cb_intra_slice_luma;

    uint8_t sps_max_mtt_hierarchy_depth_intra_slice_luma;
    uint8_t sps_log2_diff_max_bt_min_qt_intra_slice_luma;
    uint8_t sps_log2_diff_max_tt_min_qt_intra_slice_luma;

    uint8_t sps_qtbtt_dual_tree_intra_flag;

    uint8_t sps_log2_diff_min_qt_min_cb_intra_slice_chroma;

    uint8_t sps_max_mtt_hierarchy_depth_intra_slice_chroma;
    uint8_t sps_log2_diff_max_bt_min_qt_intra_slice_chroma;
    uint8_t sps_log2_diff_max_tt_min_qt_intra_slice_chroma;

    uint8_t sps_log2_diff_min_qt_min_cb_inter_slice;

    uint8_t sps_max_mtt_hierarchy_depth_inter_slice;
    uint8_t sps_log2_diff_max_bt_min_qt_inter_slice;
    uint8_t sps_log2_diff_max_tt_min_qt_inter_slice;

    uint8_t sps_max_luma_transform_size_64_flag;

    uint8_t sps_transform_skip_enabled_flag;
    uint8_t sps_log2_transform_skip_max_size_minus2;
    uint8_t sps_bdpcm_enabled_flag;

    uint8_t sps_mts_enabled_flag;
    uint8_t sps_explicit_mts_intra_enabled_flag;
    uint8_t sps_explicit_mts_inter_enabled_flag;

    uint8_t sps_lfnst_enabled_flag;

    uint8_t sps_joint_cbcr_enabled_flag;
    uint8_t sps_same_qp_table_for_chroma_flag;
    int8_t sps_qp_table_start_minus26[3]; /* i == 1 for same_qp_table_for_chroma, 2 + joint_cbrcr_flag otherwise*/
    uint8_t sps_num_points_in_qp_table_minus1[3];
    uint8_t sps_delta_qp_in_val_minus1[3][64]; /*j = max num_points_in_qp_table 64?*/
    uint8_t sps_delta_qp_diff_val[3][64]; /*j = max num_points_in_qp_table 64?*/

    uint8_t sps_sao_enabled_flag;

    uint8_t sps_alf_enabled_flag;
    uint8_t sps_ccalf_enabled_flag;

    uint8_t sps_lmcs_enabled_flag;
    uint8_t sps_weighted_pred_flag;
    uint8_t sps_weighted_bipred_flag;
    uint8_t sps_long_term_ref_pics_flag;

    uint8_t sps_inter_layer_prediction_enabled_flag;

    uint8_t sps_idr_rpl_present_flag;

    uint8_t sps_rpl1_same_as_rpl0_flag;
    uint8_t sps_num_ref_pic_lists0;
    uint8_t sps_num_ref_pic_lists1;

    struct OVRPL rpl_s0[64]; /* max num_ref_pic_list 64*/
    struct OVRPL rpl_s1[64]; /* max num_ref_pic_list 64*/

    uint8_t sps_ref_wraparound_enabled_flag;

    uint8_t sps_temporal_mvp_enabled_flag;
    uint8_t sps_sbtmvp_enabled_flag;

    uint8_t sps_amvr_enabled_flag;

    uint8_t sps_bdof_enabled_flag;
    uint8_t sps_bdof_control_present_in_ph_flag;

    uint8_t sps_smvd_enabled_flag;

    uint8_t sps_dmvr_enabled_flag;
    uint8_t sps_dmvr_control_present_in_ph_flag;

    uint8_t sps_mmvd_enabled_flag;
    uint8_t sps_mmvd_fullpel_only_enabled_flag;

    uint8_t sps_six_minus_max_num_merge_cand;
    uint8_t sps_sbt_enabled_flag;

    uint8_t sps_affine_enabled_flag;
    uint8_t sps_five_minus_max_num_subblock_merge_cand;

    uint8_t sps_6param_affine_enabled_flag;
    uint8_t sps_affine_amvr_enabled_flag;

    uint8_t sps_affine_prof_enabled_flag;
    uint8_t sps_prof_control_present_in_ph_flag;

    uint8_t sps_bcw_enabled_flag;
    uint8_t sps_ciip_enabled_flag;

    uint8_t sps_gpm_enabled_flag;
    uint8_t sps_max_num_merge_cand_minus_max_num_gpm_cand;

    uint8_t sps_log2_parallel_merge_level_minus2;
    uint8_t sps_isp_enabled_flag;
    uint8_t sps_mrl_enabled_flag;
    uint8_t sps_mip_enabled_flag;

    uint8_t sps_cclm_enabled_flag;

    uint8_t sps_chroma_horizontal_collocated_flag;
    uint8_t sps_chroma_vertical_collocated_flag;

    uint8_t sps_palette_enabled_flag;

    uint8_t sps_act_enabled_flag;

    uint8_t sps_min_qp_prime_ts;

    uint8_t sps_ibc_enabled_flag;
    uint8_t sps_six_minus_max_num_ibc_merge_cand;

    uint8_t sps_ladf_enabled_flag;
    uint8_t sps_num_ladf_intervals_minus2;
    int8_t sps_ladf_lowest_interval_qp_offset;
    int8_t sps_ladf_qp_offset[64]; /* max_num_ladf_intervals */
    uint8_t sps_ladf_delta_threshold_minus1[64];

    uint8_t sps_explicit_scaling_list_enabled_flag;
    uint8_t sps_scaling_matrix_for_lfnst_disabled_flag;
    uint8_t sps_scaling_matrix_for_alternative_colour_space_disabled_flag;
    uint8_t sps_scaling_matrix_designated_colour_space_flag;

    uint8_t sps_dep_quant_enabled_flag;
    uint8_t sps_sign_data_hiding_enabled_flag;

    uint8_t sps_virtual_boundaries_enabled_flag;
    uint8_t sps_virtual_boundaries_present_flag;
    uint8_t sps_num_ver_virtual_boundaries;
    uint8_t sps_virtual_boundary_pos_x_minus1[16]; /* max_num_virtual_boundaries_ver */
    uint8_t sps_num_hor_virtual_boundaries;
    uint8_t sps_virtual_boundary_pos_y_minus1[16]; /* max_num_virtual_boundaries_hor */

    uint8_t sps_timing_hrd_params_present_flag;
    struct OVGHRDTiming * general_timing_hrd_parameters;
    uint8_t sps_sublayer_cpb_params_present_flag;
    struct OVOLSHRDTiming * ols_timing_hrd_parameters;

    uint8_t sps_field_seq_flag;

    uint8_t sps_vui_parameters_present_flag;
    uint8_t sps_vui_payload_size_minus1;
    uint8_t sps_vui_alignment_zero_bit;

    struct OVVUI *vui_payload;

    uint8_t sps_extension_flag;
    uint8_t sps_extension_data_flag;
};

/*FIXME check values over flow */
struct OVPPS
{
    uint8_t pps_pic_parameter_set_id;
    uint8_t pps_seq_parameter_set_id;
    uint8_t pps_mixed_nalu_types_in_pic_flag;
    uint16_t pps_pic_width_in_luma_samples;
    uint16_t pps_pic_height_in_luma_samples;

    uint8_t pps_conformance_window_flag;
    uint8_t pps_conf_win_left_offset;
    uint8_t pps_conf_win_right_offset;
    uint8_t pps_conf_win_top_offset;
    uint8_t pps_conf_win_bottom_offset;

    uint8_t pps_scaling_window_explicit_signalling_flag;
    int8_t pps_scaling_win_left_offset;
    int8_t pps_scaling_win_right_offset;
    int8_t pps_scaling_win_top_offset;
    int8_t pps_scaling_win_bottom_offset;

    uint8_t pps_output_flag_present_flag;
    uint8_t pps_no_pic_partition_flag;
    uint8_t pps_subpic_id_mapping_present_flag;
    uint8_t pps_num_subpics_minus1;

    uint8_t pps_subpic_id_len_minus1;
    uint8_t pps_subpic_id[16];
    uint8_t pps_log2_ctu_size_minus5;
    uint8_t pps_num_exp_tile_columns_minus1;
    uint8_t pps_num_exp_tile_rows_minus1;
    uint8_t pps_num_tile_columns_minus1;
    uint8_t pps_num_tile_rows_minus1;
    uint8_t pps_tile_column_width_minus1[16];

    uint8_t pps_tile_row_height_minus1[16];

    uint8_t pps_loop_filter_across_tiles_enabled_flag;
    uint8_t pps_rect_slice_flag;

    uint8_t pps_single_slice_per_subpic_flag;

    uint8_t pps_num_slices_in_pic_minus1;
    uint8_t pps_tile_idx_delta_present_flag;

    uint8_t pps_slice_width_in_tiles_minus1[16];

    uint8_t pps_slice_height_in_tiles_minus1[16];

    uint8_t pps_num_exp_slices_in_tile[16];
    uint8_t pps_exp_slice_height_in_ctus_minus1[16*16];

    int8_t pps_tile_idx_delta_val[16*16];

    uint8_t pps_loop_filter_across_slices_enabled_flag;

    uint8_t pps_cabac_init_present_flag;
    uint8_t pps_num_ref_idx_default_active_minus1[2];

    uint8_t pps_rpl1_idx_present_flag;
    uint8_t pps_weighted_pred_flag;
    uint8_t pps_weighted_bipred_flag;
    uint8_t pps_ref_wraparound_enabled_flag;
    uint8_t pps_pic_width_minus_wraparound_offset;

    int8_t pps_init_qp_minus26;
    uint8_t pps_cu_qp_delta_enabled_flag;
    uint8_t pps_chroma_tool_offsets_present_flag;
    int8_t pps_cb_qp_offset;
    int8_t pps_cr_qp_offset;
    uint8_t pps_joint_cbcr_qp_offset_present_flag;
    int8_t pps_joint_cbcr_qp_offset_value;

    uint8_t pps_slice_chroma_qp_offsets_present_flag;
    uint8_t pps_cu_chroma_qp_offset_list_enabled_flag;
    uint8_t pps_chroma_qp_offset_list_len_minus1;
    int8_t pps_cb_qp_offset_list[16];
    int8_t pps_cr_qp_offset_list[16];
    int8_t pps_joint_cbcr_qp_offset_list[16];

    uint8_t pps_deblocking_filter_control_present_flag;
    uint8_t pps_deblocking_filter_override_enabled_flag;
    uint8_t pps_deblocking_filter_disabled_flag;
    uint8_t pps_dbf_info_in_ph_flag;

    int8_t pps_luma_beta_offset_div2;
    int8_t pps_luma_tc_offset_div2;
    int8_t pps_cb_beta_offset_div2;
    int8_t pps_cb_tc_offset_div2;
    int8_t pps_cr_beta_offset_div2;
    int8_t pps_cr_tc_offset_div2;

    uint8_t pps_rpl_info_in_ph_flag;
    uint8_t pps_sao_info_in_ph_flag;
    uint8_t pps_alf_info_in_ph_flag;
    uint8_t pps_wp_info_in_ph_flag;

    uint8_t pps_qp_delta_info_in_ph_flag;

    uint8_t pps_picture_header_extension_present_flag;
    uint8_t pps_slice_header_extension_present_flag;
    uint8_t pps_extension_flag;

    uint8_t pps_extension_data_flag;
};

struct OVPH
{
    uint8_t ph_gdr_or_irap_pic_flag;
    uint8_t ph_non_ref_pic_flag;
    uint8_t ph_gdr_pic_flag;
    uint8_t ph_inter_slice_allowed_flag;
    uint8_t ph_intra_slice_allowed_flag;
    uint8_t ph_pic_parameter_set_id;
    uint8_t ph_pic_order_cnt_lsb;
    uint8_t ph_recovery_poc_cnt;
    /* FIXME ignored */
    uint8_t ph_extra_bit[24];
    uint8_t ph_poc_msb_cycle_present_flag;
    uint8_t ph_poc_msb_cycle_val;
    uint8_t ph_alf_enabled_flag;
    uint8_t ph_num_alf_aps_ids_luma;
    uint8_t ph_alf_aps_id_luma[8];
    uint8_t ph_alf_cb_enabled_flag;
    uint8_t ph_alf_cr_enabled_flag;
    uint8_t ph_alf_aps_id_chroma;
    uint8_t ph_alf_cc_cb_enabled_flag;
    uint8_t ph_alf_cc_cb_aps_id;
    uint8_t ph_alf_cc_cr_enabled_flag;
    uint8_t ph_alf_cc_cr_aps_id;
    uint8_t ph_lmcs_enabled_flag;
    uint8_t ph_lmcs_aps_id;
    uint8_t ph_chroma_residual_scale_flag;
    uint8_t ph_explicit_scaling_list_enabled_flag;
    uint8_t ph_scaling_list_aps_id;
    uint8_t ph_virtual_boundaries_present_flag;
    uint8_t ph_num_ver_virtual_boundaries;
    uint8_t ph_virtual_boundary_pos_x_minus1[32];
    uint8_t ph_num_hor_virtual_boundaries;
    uint8_t ph_virtual_boundary_pos_y_minus1[32];
    uint8_t ph_pic_output_flag;
    uint8_t ph_partition_constraints_override_flag;
    uint8_t ph_log2_diff_min_qt_min_cb_intra_slice_luma;
    uint8_t ph_max_mtt_hierarchy_depth_intra_slice_luma;
    uint8_t ph_log2_diff_max_bt_min_qt_intra_slice_luma;
    uint8_t ph_log2_diff_max_tt_min_qt_intra_slice_luma;
    uint8_t ph_log2_diff_min_qt_min_cb_intra_slice_chroma;
    uint8_t ph_max_mtt_hierarchy_depth_intra_slice_chroma;
    uint8_t ph_log2_diff_max_bt_min_qt_intra_slice_chroma;
    uint8_t ph_log2_diff_max_tt_min_qt_intra_slice_chroma;
    uint8_t ph_cu_qp_delta_subdiv_intra_slice;
    uint8_t ph_cu_chroma_qp_offset_subdiv_intra_slice;
    uint8_t ph_log2_diff_min_qt_min_cb_inter_slice;
    uint8_t ph_max_mtt_hierarchy_depth_inter_slice;
    uint8_t ph_log2_diff_max_bt_min_qt_inter_slice;
    uint8_t ph_log2_diff_max_tt_min_qt_inter_slice;
    uint8_t ph_cu_qp_delta_subdiv_inter_slice;
    uint8_t ph_cu_chroma_qp_offset_subdiv_inter_slice;
    uint8_t ph_temporal_mvp_enabled_flag;
    uint8_t ph_collocated_ref_idx;
    uint8_t ph_mmvd_fullpel_only_flag;
    uint8_t ph_mvd_l1_zero_flag;
    uint8_t ph_bdof_disabled_flag;
    uint8_t ph_dmvr_disabled_flag;
    uint8_t ph_prof_disabled_flag;
    int8_t ph_qp_delta;
    uint8_t ph_joint_cbcr_sign_flag;
    uint8_t ph_sao_luma_enabled_flag;
    uint8_t ph_sao_chroma_enabled_flag;
    uint8_t ph_deblocking_params_present_flag;
    uint8_t ph_deblocking_filter_disabled_flag;
    int8_t ph_luma_beta_offset_div2;
    int8_t ph_luma_tc_offset_div2;
    int8_t ph_cb_beta_offset_div2;
    int8_t ph_cb_tc_offset_div2;
    int8_t ph_cr_beta_offset_div2;
    int8_t ph_cr_tc_offset_div2;
    uint8_t ph_collocated_from_l0_flag;
    uint8_t ph_extension_length;
    /* FIXME IGNORED */
    uint8_t ph_extension_data_byte[8];
    /* Ref Pic List Info */
    struct OVHRPL hrpl;
};

struct OVSH
{
    uint8_t sh_picture_header_in_slice_header_flag;
    uint16_t sh_subpic_id;
    uint16_t sh_slice_address;
    /* Could be reduced unused */
    uint8_t sh_extra_bit[64];
    uint16_t sh_num_tiles_in_slice_minus1;
    uint8_t sh_slice_type;
    uint8_t sh_no_output_of_prior_pics_flag;
    uint8_t sh_alf_enabled_flag;
    uint8_t sh_num_alf_aps_ids_luma;
    uint8_t sh_alf_aps_id_luma[8];
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
    uint8_t sh_num_ref_idx_active_l0_minus1;
    uint8_t sh_num_ref_idx_active_l1_minus1;
    uint8_t sh_cabac_init_flag;
    uint8_t sh_collocated_from_l0_flag;
    uint8_t sh_collocated_ref_idx;
    int8_t sh_qp_delta;
    int8_t sh_cb_qp_offset;
    int8_t sh_cr_qp_offset;
    int8_t sh_joint_cbcr_qp_offset;
    uint8_t sh_cu_chroma_qp_offset_enabled_flag;
    uint8_t sh_sao_luma_used_flag;
    uint8_t sh_sao_chroma_used_flag;
    uint8_t sh_deblocking_params_present_flag;
    uint8_t sh_deblocking_filter_disabled_flag;
    int8_t sh_luma_beta_offset_div2;
    int8_t sh_luma_tc_offset_div2;
    int8_t sh_cb_beta_offset_div2;
    int8_t sh_cb_tc_offset_div2;
    int8_t sh_cr_beta_offset_div2;
    int8_t sh_cr_tc_offset_div2;
    uint8_t sh_dep_quant_used_flag;
    uint8_t sh_sign_data_hiding_used_flag;
    uint8_t sh_ts_residual_coding_disabled_flag;
    uint8_t sh_slice_header_extension_length;
    /* unused */
    uint8_t sh_slice_header_extension_data_byte[8];
    uint8_t sh_entry_offset_len_minus1;
    /* FIXME find bounds on number of entry points */
    uint32_t sh_entry_point_offset_minus1[64];
    /* Ref pic list info */
    struct OVHRPL hrpl;
};


typedef struct OVALFData
{
    uint8_t alf_luma_filter_signal_flag;
    uint8_t alf_chroma_filter_signal_flag;
    uint8_t alf_luma_clip_flag;
    uint8_t alf_luma_num_filters_signalled_minus1;
    uint8_t alf_cc_cb_filter_signal_flag;
    uint8_t alf_cc_cr_filter_signal_flag;

    uint8_t alf_luma_coeff_delta_idx[MAX_NUM_ALF_CLASSES];
    int16_t alf_luma_coeff[MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_LUMA_COEFF];
    int16_t alf_luma_clip_idx[MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_LUMA_COEFF];
    uint8_t alf_chroma_clip_flag;
    uint8_t alf_chroma_num_alt_filters_minus1;
    int16_t alf_chroma_coeff[MAX_NUM_ALF_ALTERNATIVES_CHROMA][MAX_NUM_ALF_CHROMA_COEFF];
    int16_t alf_chroma_clip_idx[MAX_NUM_ALF_ALTERNATIVES_CHROMA][MAX_NUM_ALF_CHROMA_COEFF];

    uint8_t alf_cc_cb_filters_signalled_minus1;
    uint8_t alf_cc_cr_filters_signalled_minus1;
    int16_t alf_cc_mapped_coeff[2][MAX_NUM_CC_ALF_FILTERS][MAX_NUM_CC_ALF_CHROMA_COEFF];
} OVALFData;

typedef struct OVLMCSData
{
    uint8_t lmcs_min_bin_idx;
    uint8_t lmcs_delta_max_bin_idx;
    uint8_t lmcs_delta_cw_prec_minus1;
    uint8_t lmcs_delta_abs_cw[PIC_CODE_CW_BINS];
    uint8_t lmcs_delta_sign_cw_flag[PIC_CODE_CW_BINS];
    uint8_t lmcs_delta_abs_crs;
    uint8_t lmcs_delta_sign_crs_flag;
} OVLMCSData;

typedef struct OVAPS
{
    uint8_t aps_params_type;
    uint8_t aps_adaptation_parameter_set_id;
    uint8_t aps_chroma_present_flag;
    /* Note unused */
    uint8_t aps_extension_flag;
    uint8_t aps_extension_data_flag;

    struct OVALFData  aps_alf_data;
    struct OVLMCSData aps_lmcs_data;
} OVAPS;


#endif
