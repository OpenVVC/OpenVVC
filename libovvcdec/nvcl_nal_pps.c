#include "libovvcutils/ovmem.h"

#include "nvcl.h"
#include "nvcl_utils.h"
#include "nvcl_structures.h"

enum DecReturn {
    OV_INVALID_DATA = -1,
    OV_ENOMEM = -2,
};


static uint8_t
probe_pps_id(OVNVCLReader *const rdr)
{
    uint8_t pps_id = fetch_bits(rdr, 6);
    return pps_id;
}


static int
validate_pps(OVNVCLReader *rdr, OVPPS *const pps)
{
    /* TODO various check on limitation and max sizes */
    return 1;
}

static void
free_pps(OVPPS *const pps)
{
    /* TODO unref and/or free dynamic structure */
    ov_free(pps);
}

static void
replace_pps(OVPPS *pps_list[], OVPPS *const pps, uint8_t pps_id)
{
    /* TODO unref and/or free dynamic structure */
    OVPPS *to_free = pps_list[pps_id];

    free_pps(to_free);

    pps_list[pps_id] = pps;
}

int
nvcl_decode_nalu_pps(OVNVCLReader *const rdr, OVNVCLCtx *const nvcl_ctx)
{
    int ret;
    uint8_t pps_id = probe_pps_id(rdr);
    OVPPS **pps_list = nvcl_ctx->pps_list;
    if (pps_list[pps_id]) {
        /* TODO compare RBSP data to avoid new read */
        uint8_t identical_rbsp = 0;
        if (identical_rbsp) {
            goto duplicated;
        }
    }

    OVPPS *pps = ov_mallocz(sizeof(*pps));
    if (!pps) {
        return OV_ENOMEM;
    }

    ret = nvcl_pps_read(rdr, pps, nvcl_ctx);
    if (ret < 0) {
        goto cleanup;
    }

    ret = validate_pps(rdr, pps);
    if (ret < 0) {
        goto cleanup;
    }

    /*FIXME unref instead of free */
    replace_pps(pps_list, pps, pps_id);

    return 0;

cleanup:
    ov_free(pps);
    return ret;

duplicated:
    #if 0
    ov_log(NULL, 3, "Ignored Duplicated PPS");
    #endif
    return 0;
}
static void
pps_read_slices_in_subpic(OVNVCLReader *const rdr, OVPPS *const pps)
{
    int i;
    int tile_id = 0;

    /* FIXME not tested */
    pps->pps_num_slices_in_pic_minus1 = nvcl_read_u_expgolomb(rdr);
    if (pps->pps_num_slices_in_pic_minus1 > 1){
        pps->pps_tile_idx_delta_present_flag = nvcl_read_flag(rdr);
    }

    for (i = 0; i < pps->pps_num_slices_in_pic_minus1; i++){
        int tile_pos_x = tile_id % (pps->pps_num_exp_tile_columns_minus1 + 1);
        int tile_pos_y = tile_id / (pps->pps_num_exp_tile_columns_minus1 + 1);

        /* Each new tile column read slice width exept for implicit last column */
        if (tile_pos_x != pps->pps_num_exp_tile_columns_minus1){
            pps->pps_slice_width_in_tiles_minus1[i] = nvcl_read_u_expgolomb(rdr);
        }

        /* Each new tile row read slice height except for implicit last row */
        if (tile_pos_y !=  pps->pps_num_exp_tile_rows_minus1){
            if(pps->pps_tile_idx_delta_present_flag || tile_pos_x == 0){
                pps->pps_slice_height_in_tiles_minus1[i] = nvcl_read_u_expgolomb(rdr);
            }
        }

        /* Multiple slices in tiles */
        if (!pps->pps_slice_width_in_tiles_minus1[i] && !pps->pps_slice_height_in_tiles_minus1[i]){
            if (pps->pps_tile_row_height_minus1[tile_pos_y] > 1){
                pps->pps_num_exp_slices_in_tile[i] = nvcl_read_u_expgolomb(rdr);
                if (pps->pps_num_exp_slices_in_tile[i]){
                    int j;
                    for (j = 0; j < pps->pps_num_exp_slices_in_tile[i]; j++){
                        pps->pps_exp_slice_height_in_ctus_minus1[i + j] = nvcl_read_u_expgolomb(rdr);
                    }
                    i += (j - 1);
                }
            }
        }

        if (pps->pps_tile_idx_delta_present_flag){
            pps->pps_tile_idx_delta_val[i] = nvcl_read_s_expgolomb(rdr);
            tile_id += pps->pps_tile_idx_delta_val[i];
        } else {
            int offset_y;
            tile_id += pps->pps_slice_width_in_tiles_minus1[i] + 1;
            offset_y  = pps->pps_slice_height_in_tiles_minus1[i];
            offset_y *= pps->pps_num_exp_tile_columns_minus1 + 1;
            if (tile_id % (pps->pps_num_exp_tile_columns_minus1 + 1) == 0){
                tile_id += offset_y;
            }
        }
    }
}

static void
pps_read_pic_partition(OVNVCLReader *const rdr, OVPPS *const pps)
{
    int i;
    int row_sum = 0;
    int col_sum = 0;
    #if 0
    const int log2_ctu_s = pps->pps_log2_ctu_size_minus5 + 5;
    const int pic_w = pps->pps_pic_width_in_luma_samples;
    const int pic_h = pps->pps_pic_height_in_luma_samples;

    const int nb_ctu_w = (pic_w >> log2_ctu_s) + !!(pic_w & ((1 << log2_ctu_s) - 1));
    const int nb_ctu_h = (pic_h >> log2_ctu_s) + !!(pic_h & ((1 << log2_ctu_s) - 1));
    #endif


    pps->pps_num_exp_tile_columns_minus1 = nvcl_read_u_expgolomb(rdr);
    pps->pps_num_exp_tile_rows_minus1    = nvcl_read_u_expgolomb(rdr);

    for (i = 0; i <=  pps->pps_num_exp_tile_columns_minus1; i++){
        pps->pps_tile_column_width_minus1[i] = nvcl_read_u_expgolomb(rdr);
        col_sum += pps->pps_tile_column_width_minus1[i] + 1;
    }

    for (i = 0; i <=  pps->pps_num_exp_tile_rows_minus1; i++){
        pps->pps_tile_row_height_minus1[i] = nvcl_read_u_expgolomb(rdr);
        row_sum += pps->pps_tile_row_height_minus1[i] + 1;
    }

    /* Default initialisation values
     * TODO confirm this is correct
     */
    pps->pps_loop_filter_across_tiles_enabled_flag = 1;
    pps->pps_rect_slice_flag                       = 1;

    /* if more than one tile tiles are enabled
     * FIXME confirm check of implicit tiles on col_sum and row_sum
     * is required / correct
     */
    if (pps->pps_num_exp_tile_columns_minus1 || pps->pps_num_exp_tile_rows_minus1
                  /*|| col_sum != nb_ctu_w || row_sum != nb_ctu_h*/) {
        pps->pps_loop_filter_across_tiles_enabled_flag = nvcl_read_flag(rdr);
        pps->pps_rect_slice_flag                       = nvcl_read_flag(rdr);
    }

    if (pps->pps_rect_slice_flag){
        pps->pps_single_slice_per_subpic_flag = nvcl_read_flag(rdr);
    }

    if (pps->pps_rect_slice_flag && !pps->pps_single_slice_per_subpic_flag){
        pps_read_slices_in_subpic(rdr, pps);
    }

    if (!pps->pps_rect_slice_flag || pps->pps_single_slice_per_subpic_flag ||
            pps->pps_num_slices_in_pic_minus1){
        pps->pps_loop_filter_across_slices_enabled_flag = nvcl_read_flag(rdr);
    }
}


int
nvcl_pps_read(OVNVCLReader *const rdr, OVPPS *const pps,
              OVNVCLCtx *const nvcl_ctx)
{
    int i;

    #if 0
    OVPPS stck_pps = {0};
    OVPPS *const pps = &stck_pps;
    #endif
    nvcl_skip_bits(rdr, 16);

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

        /* TODO confirm sub functions*/
        #if 0
        int nb_tiles_pic;
        pps->pps_log2_ctu_size_minus5 = nvcl_read_bits(rdr, 2);

        pps->pps_num_exp_tile_columns_minus1 = nvcl_read_u_expgolomb(rdr);
        pps->pps_num_exp_tile_rows_minus1    = nvcl_read_u_expgolomb(rdr);

        nb_tiles_pic = (pps->pps_num_exp_tile_columns_minus1 + 1)
                       * (pps->pps_num_exp_tile_rows_minus1 + 1);
        for(i = 0; i <= pps->pps_num_exp_tile_columns_minus1; i++) {
            pps->pps_tile_column_width_minus1[i] = nvcl_read_u_expgolomb(rdr);
        }

        for(i = 0; i <= pps->pps_num_exp_tile_rows_minus1; i++) {
            pps->pps_tile_row_height_minus1[i] = nvcl_read_u_expgolomb(rdr);
        }

        if (nb_tiles_pic > 1) {
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
    #else
    pps->pps_log2_ctu_size_minus5 = nvcl_read_bits(rdr, 2);

    pps_read_pic_partition(rdr, pps);
    #endif
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
        #if 0
        while (more_rbsp_data()) {
            pps->pps_extension_data_flag = nvcl_read_flag(rdr);
        }
        #endif
    }

    #if 0
    rbsp_trailing_bits()
    #endif
    return 0;
}
