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

#include <string.h>

#include "ovmem.h"
#include "overror.h"

#include "hls_structures.h"
#include "nvcl.h"
#include "nvcl_utils.h"
#include "nvcl_structures.h"
#include "nvcl_private.h"

//TODOgpm: do not init here
#include "rcn.h"

static uint8_t
probe_sps_id(OVNVCLReader *const rdr)
{
    uint8_t sps_id = fetch_bits(rdr, 4);
    return sps_id;
}

static struct HLSDataRef **
storage_in_nvcl_ctx(OVNVCLReader *const rdr, OVNVCLCtx *const nvcl_ctx)
{
    uint8_t id = probe_sps_id(rdr);
    struct HLSDataRef **storage = &nvcl_ctx->sps_list[id];

    return storage;
}

static int
validate_sps(OVNVCLReader *rdr, const union HLSData *const data)
{
    /* TODO various check on limitation and max sizes */
    const OVSPS *const sps =  (const OVSPS *)data;
    uint32_t nb_bits_read = nvcl_nb_bits_read(rdr) + 1;
    uint32_t stop_bit_pos = nvcl_find_rbsp_stop_bit(rdr);
    if (stop_bit_pos != nb_bits_read) {

        ov_log(NULL, OVLOG_ERROR, "rbsp_stop_bit mismatch: cursor at %d,  expected %d\n", nb_bits_read, stop_bit_pos);
        return OVVC_EINDATA;
    }


    if (sps->sps_weighted_pred_flag || sps->sps_weighted_bipred_flag) {
        ov_log(NULL, OVLOG_WARNING, "Unsupported weighted pred\n");
    }

    if (sps->sps_subpic_info_present_flag) {
        ov_log(NULL, OVLOG_ERROR, "Unsupported subpicture\n");
        return OVVC_EINDATA;
    }

    if (sps->sps_entropy_coding_sync_enabled_flag) {
        ov_log(NULL, OVLOG_ERROR, "Unsupported WPP\n");
        return OVVC_EINDATA;
    }

    if (sps->sps_explicit_scaling_list_enabled_flag) {
        ov_log(NULL, OVLOG_WARNING, "Unsupported scaling lists\n");
    }

    if (sps->sps_long_term_ref_pics_flag) {
        ov_log(NULL, OVLOG_ERROR, "Unsupported long term references\n");
        return OVVC_EINDATA;
    }

    if (sps->sps_pic_width_max_in_luma_samples & 0xE003 || sps->sps_pic_height_max_in_luma_samples & 0xE003) {
        ov_log(NULL, OVLOG_ERROR, "Invalid maximum picture dimension %dx%d in SPS.\n",
               sps->sps_pic_width_max_in_luma_samples,
               sps->sps_pic_height_max_in_luma_samples);
        return OVVC_EINDATA;
    }

    if ((sps->sps_num_points_in_qp_table_minus1[0] |
         sps->sps_num_points_in_qp_table_minus1[1] |
         sps->sps_num_points_in_qp_table_minus1[2]) & ~0x3F) {

        ov_log(NULL, OVLOG_ERROR, "Invalid sps_num_points_in_qp_table_minus1 %d %d %d in SPS.\n",
               sps->sps_num_points_in_qp_table_minus1[0],
               sps->sps_num_points_in_qp_table_minus1[1],
               sps->sps_num_points_in_qp_table_minus1[2]);
        return OVVC_EINDATA;
    }

    if (sps->sps_six_minus_max_num_merge_cand > 6) {
        ov_log(NULL, OVLOG_ERROR, "Invalid sps_six_minus_max_num_merge_cand %d in SPS.\n",
               sps->sps_six_minus_max_num_merge_cand);
        return OVVC_EINDATA;
    }

    return 1;
}

static void
free_sps(const union HLSData *const data)
{
    /* TODO unref and/or free dynamic structure */
    const OVSPS *const sps = (const OVSPS *)data;
    ov_free((void *)sps);
}

#if 0
static int
replace_sps(const struct HLSReader *const manager,
            struct HLSDataRef **storage,
            const OVHLSData *const hls_data)
{
    /* TODO unref and/or free dynamic structure */
    const union HLSData *to_free = *storage;
    union HLSData *new = ov_malloc(manager->data_size);

    if (!new) {
        return OVVC_ENOMEM;
    }

    memcpy(new, hls_data, manager->data_size);

    *storage = new;

    if (to_free) {
        manager->free(to_free);
    }

    return 0;
}
#endif

static void
subpic_info(OVNVCLReader *const rdr, OVSPS *const sps)
{
    sps->sps_num_subpics_minus1 = nvcl_read_u_expgolomb(rdr);
    if (sps->sps_num_subpics_minus1 > 0) {
        int i;
        int pic_w = sps->sps_pic_width_max_in_luma_samples;
        int pic_h = sps->sps_pic_height_max_in_luma_samples;
        int log2_ctb_s = sps->sps_log2_ctu_size_minus5 + 5;
        int nb_ctu_w = (pic_w + ((1 << log2_ctb_s) - 1)) >> log2_ctb_s;
        int nb_ctu_h = (pic_h + ((1 << log2_ctb_s) - 1)) >> log2_ctb_s;
        int ctb_s = 1 << log2_ctb_s;
        sps->sps_independent_subpics_flag = nvcl_read_flag(rdr);
        sps->sps_subpic_same_size_flag = nvcl_read_flag(rdr);
        for (i = 0; i < 1; i++) {

            /*TODO check ctb_w and ctb_h conditions once outside of loop*/
            /*FIXME ceil_log2 */
            if (sps->sps_pic_width_max_in_luma_samples > ctb_s) {
                int v = ov_ceil_log2(nb_ctu_w);
                sps->sps_subpic_width_minus1[i] = nvcl_read_bits(rdr,v);
            }

            if (sps->sps_pic_height_max_in_luma_samples > ctb_s) {
                int v = ov_ceil_log2(nb_ctu_h);
                sps->sps_subpic_height_minus1[i] = nvcl_read_bits(rdr,v);
            }

            if (!sps->sps_independent_subpics_flag) {
                sps->sps_subpic_treated_as_pic_flag[i] = nvcl_read_flag(rdr);
                sps->sps_loop_filter_across_subpic_enabled_flag[i] = nvcl_read_flag(rdr);
            }
        }

        if (!sps->sps_subpic_same_size_flag) {
            /*TODO check ctb_w and ctb_h conditions once outside of loop*/
            /*FIXME ensure sps_subpic_ctu_top_left_*[0] are init at 0*/
            for (i = 1; i <= (sps->sps_num_subpics_minus1 & 0xF); i++) {

                if (sps->sps_pic_width_max_in_luma_samples > ctb_s) {
                    int v = ov_ceil_log2(nb_ctu_w);
                    sps->sps_subpic_ctu_top_left_x[i] = nvcl_read_bits(rdr, v);
                }

                if (sps->sps_pic_height_max_in_luma_samples > ctb_s) {
                    int v = ov_ceil_log2(nb_ctu_h);
                    sps->sps_subpic_ctu_top_left_y[i] = nvcl_read_bits(rdr, v);
                }

                if (i < sps->sps_num_subpics_minus1 && sps->sps_pic_width_max_in_luma_samples > ctb_s) {
                    int v = ov_ceil_log2(nb_ctu_w);
                    sps->sps_subpic_width_minus1[i] = nvcl_read_bits(rdr, v);
                }

                if (i < sps->sps_num_subpics_minus1 && sps->sps_pic_height_max_in_luma_samples > ctb_s) {
                    int v = ov_ceil_log2(nb_ctu_h);
                    sps->sps_subpic_height_minus1[i] = nvcl_read_bits(rdr, v);
                }

                if (!sps->sps_independent_subpics_flag) {
                    sps->sps_subpic_treated_as_pic_flag[i] = nvcl_read_flag(rdr);
                    sps->sps_loop_filter_across_subpic_enabled_flag[i] = nvcl_read_flag(rdr);
                }
            }
        }
    }

    sps->sps_subpic_id_len_minus1 = nvcl_read_u_expgolomb(rdr);
    sps->sps_subpic_id_mapping_explicitly_signalled_flag = nvcl_read_flag(rdr);

    if (sps->sps_subpic_id_mapping_explicitly_signalled_flag) {
        sps->sps_subpic_id_mapping_present_flag = nvcl_read_flag(rdr);
        if (sps->sps_subpic_id_mapping_present_flag) {
            int i;
            for (i = 0; i <= (sps->sps_num_subpics_minus1 & 0xF); i++) {
                sps->sps_subpic_id[i] = nvcl_read_bits(rdr, sps->sps_subpic_id_len_minus1 + 1);
            }
        }
    }
}


static void
vui_payload(OVNVCLReader *const rdr, struct OVVUI *vui)
{
    vui->vui_progressive_source_flag        = nvcl_read_flag(rdr);
    vui->vui_interlaced_source_flag         = nvcl_read_flag(rdr);
    vui->vui_non_packed_constraint_flag     = nvcl_read_flag(rdr);
    vui->vui_non_projected_constraint_flag  = nvcl_read_flag(rdr);

    vui->vui_aspect_ratio_info_present_flag = nvcl_read_flag(rdr);
    if (vui->vui_aspect_ratio_info_present_flag) {
        vui->vui_aspect_ratio_constant_flag = nvcl_read_flag(rdr);
        vui->vui_aspect_ratio_idc = nvcl_read_bits(rdr, 8);
        if (vui->vui_aspect_ratio_idc == 255) {
            vui->vui_sar_width  = nvcl_read_bits(rdr, 16);
            vui->vui_sar_height = nvcl_read_bits(rdr, 16);
        }
    }

    vui->vui_overscan_info_present_flag = nvcl_read_flag(rdr);
    if (vui->vui_overscan_info_present_flag) {
        vui->vui_overscan_appropriate_flag = nvcl_read_flag(rdr);
    }

    vui->vui_colour_description_present_flag = nvcl_read_flag(rdr);
    if (vui->vui_colour_description_present_flag) {
        vui->vui_colour_primaries         = nvcl_read_bits(rdr, 8);
        vui->vui_transfer_characteristics = nvcl_read_bits(rdr, 8);
        vui->vui_matrix_coeffs            = nvcl_read_bits(rdr, 8);
        vui->vui_full_range_flag = nvcl_read_flag(rdr);
    }

    vui->vui_chroma_loc_info_present_flag = nvcl_read_flag(rdr);
    if (vui->vui_chroma_loc_info_present_flag) {
        if (vui->vui_progressive_source_flag && !vui->vui_interlaced_source_flag) {
            vui->vui_chroma_sample_loc_type_frame = nvcl_read_u_expgolomb(rdr);
        } else {
            vui->vui_chroma_sample_loc_type_top_field = nvcl_read_u_expgolomb(rdr);
            vui->vui_chroma_sample_loc_type_bottom_field = nvcl_read_u_expgolomb(rdr);
        }
    }
}

struct HRDTiming
{
    uint32_t num_units_in_tick;
    uint32_t time_scale;
    uint8_t general_nal_hrd_params_present_flag;
    uint8_t general_vcl_hrd_params_present_flag;
    uint8_t general_same_pic_timing_in_all_ols_flag;
    uint8_t general_du_hrd_params_present_flag;
    uint8_t tick_divisor_minus2;
    uint8_t bit_rate_scale;
    uint8_t cpb_size_scale;
    uint8_t cpb_size_du_scale;
    uint16_t hrd_cpb_cnt_minus1;
};

static void
general_timing_hrd_parameters(OVNVCLReader *const rdr, struct HRDTiming *hrd)
{
    hrd->num_units_in_tick = nvcl_read_bits(rdr, 32);
    hrd->time_scale = nvcl_read_bits(rdr, 32);
    hrd->general_nal_hrd_params_present_flag = nvcl_read_flag(rdr);
    hrd->general_vcl_hrd_params_present_flag = nvcl_read_flag(rdr);
    if (hrd->general_nal_hrd_params_present_flag || hrd->general_vcl_hrd_params_present_flag) {
        hrd->general_same_pic_timing_in_all_ols_flag = nvcl_read_flag(rdr);
        hrd->general_du_hrd_params_present_flag = nvcl_read_flag(rdr);
        if (hrd->general_du_hrd_params_present_flag) {
            hrd->tick_divisor_minus2 = nvcl_read_bits(rdr, 8);
        }

        hrd->bit_rate_scale = nvcl_read_bits(rdr, 4);
        hrd->cpb_size_scale = nvcl_read_bits(rdr, 4);

        if (hrd->general_du_hrd_params_present_flag) {
            hrd->cpb_size_du_scale  = nvcl_read_bits(rdr,4);
        }
        hrd->hrd_cpb_cnt_minus1 = nvcl_read_u_expgolomb(rdr);
    }
}

static void
sublayer_hrd_parameters(OVNVCLReader *const rdr, const struct HRDTiming *const hrd, uint8_t i) {
    int j;
    for (j = 0; j <= hrd->hrd_cpb_cnt_minus1; ++j) {
        uint16_t bit_rate_value_minus1 = nvcl_read_u_expgolomb(rdr);
        uint16_t cpb_size_value_minus1 = nvcl_read_u_expgolomb(rdr);
        if (hrd->general_du_hrd_params_present_flag) {
            uint16_t cpb_size_du_value_minus1 = nvcl_read_u_expgolomb(rdr);
            uint16_t bit_rate_du_value_minus1 = nvcl_read_u_expgolomb(rdr);
        }
        uint8_t cbr_flag = nvcl_read_flag(rdr);
    }
}

static void
ols_timing_hrd_parameters(OVNVCLReader *const rdr, const struct HRDTiming *const hrd, uint8_t first_sublayer, uint8_t max_sublayer_min1)
{
    int i;
    for (i = first_sublayer; i <= max_sublayer_min1; ++i) {
        uint8_t fixed_pic_rate_general_flag = nvcl_read_flag(rdr);

        uint8_t fixed_pic_rate_within_cvs_flag = fixed_pic_rate_general_flag;
        if (!fixed_pic_rate_general_flag) {
            fixed_pic_rate_within_cvs_flag = nvcl_read_flag(rdr);
        }

        if (fixed_pic_rate_within_cvs_flag) {
            uint32_t elemental_duration_in_tc_minus1 = nvcl_read_s_expgolomb(rdr);
        } else if ((hrd->general_nal_hrd_params_present_flag || hrd->general_vcl_hrd_params_present_flag)
                   && !hrd->hrd_cpb_cnt_minus1) {
            uint8_t low_delay_hrd_flag = nvcl_read_flag(rdr);
        }

        if (hrd->general_nal_hrd_params_present_flag) {
            sublayer_hrd_parameters(rdr, hrd, i);
        }

        if (hrd->general_vcl_hrd_params_present_flag) {
            sublayer_hrd_parameters(rdr, hrd, i);
        }
    }
}

int
nvcl_sps_read(OVNVCLReader *const rdr, OVHLSData *const hls_data,
              const OVNVCLCtx *const nvcl_ctx, uint8_t nalu_type)
{
    int i, j;
    OVSPS *const sps = &hls_data->sps;

    sps->sps_seq_parameter_set_id   = nvcl_read_bits(rdr, 4);
    sps->sps_video_parameter_set_id = nvcl_read_bits(rdr, 4);
    sps->sps_max_sublayers_minus1   = nvcl_read_bits(rdr, 3);
    sps->sps_chroma_format_idc      = nvcl_read_bits(rdr, 2);
    sps->sps_log2_ctu_size_minus5   = nvcl_read_bits(rdr, 2);

    sps->sps_ptl_dpb_hrd_params_present_flag = nvcl_read_flag(rdr);
    if (sps->sps_ptl_dpb_hrd_params_present_flag) {
        profile_tier_level_sps(rdr, sps->sps_max_sublayers_minus1);
    }

    sps->sps_gdr_enabled_flag = nvcl_read_flag(rdr);

    sps->sps_ref_pic_resampling_enabled_flag = nvcl_read_flag(rdr);
    if (sps->sps_ref_pic_resampling_enabled_flag) {
        sps->sps_res_change_in_clvs_allowed_flag = nvcl_read_flag(rdr);
    }

    sps->sps_pic_width_max_in_luma_samples  = nvcl_read_u_expgolomb(rdr);
    sps->sps_pic_height_max_in_luma_samples = nvcl_read_u_expgolomb(rdr);

    sps->sps_conformance_window_flag = nvcl_read_flag(rdr);
    if (sps->sps_conformance_window_flag) {
        sps->sps_conf_win_left_offset   = nvcl_read_u_expgolomb(rdr);
        sps->sps_conf_win_right_offset  = nvcl_read_u_expgolomb(rdr);
        sps->sps_conf_win_top_offset    = nvcl_read_u_expgolomb(rdr);
        sps->sps_conf_win_bottom_offset = nvcl_read_u_expgolomb(rdr);
    }

    sps->sps_subpic_info_present_flag = nvcl_read_flag(rdr);
    if (sps->sps_subpic_info_present_flag) {
    /* FIXME move to specialized function */
       subpic_info(rdr, sps);
    }

    sps->sps_bitdepth_minus8 = nvcl_read_u_expgolomb(rdr);
    sps->sps_entropy_coding_sync_enabled_flag = nvcl_read_flag(rdr);
    sps->sps_entry_point_offsets_present_flag = nvcl_read_flag(rdr);
    sps->sps_log2_max_pic_order_cnt_lsb_minus4 = nvcl_read_bits(rdr, 4);

    sps->sps_poc_msb_cycle_flag = nvcl_read_flag(rdr);
    if (sps->sps_poc_msb_cycle_flag) {
        sps->sps_poc_msb_cycle_len_minus1 = nvcl_read_u_expgolomb(rdr);
    }

    sps->sps_num_extra_ph_bytes = nvcl_read_bits(rdr, 2);
    for (i = 0; i < (sps->sps_num_extra_ph_bytes * 8); i++) {
        sps->sps_extra_ph_bit_present_flag[i] = nvcl_read_flag(rdr);
    }

    sps->sps_num_extra_sh_bytes = nvcl_read_bits(rdr, 2);
    for (i = 0; i < (sps->sps_num_extra_sh_bytes * 8); i++) {
        sps->sps_extra_sh_bit_present_flag[i] = nvcl_read_flag(rdr);
    }

    if (sps->sps_ptl_dpb_hrd_params_present_flag) {

        if (sps->sps_max_sublayers_minus1 > 0) {
            sps->sps_sublayer_dpb_params_flag = nvcl_read_flag(rdr);
        }

        dpb_parameters(rdr, sps->dpb_parameters, sps->sps_max_sublayers_minus1, sps->sps_sublayer_dpb_params_flag);
    }

    sps->sps_log2_min_luma_coding_block_size_minus2 = nvcl_read_u_expgolomb(rdr);
    sps->sps_partition_constraints_override_enabled_flag = nvcl_read_flag(rdr);
    sps->sps_log2_diff_min_qt_min_cb_intra_slice_luma = nvcl_read_u_expgolomb(rdr);

    sps->sps_max_mtt_hierarchy_depth_intra_slice_luma = nvcl_read_u_expgolomb(rdr);
    if (sps->sps_max_mtt_hierarchy_depth_intra_slice_luma != 0) {
        sps->sps_log2_diff_max_bt_min_qt_intra_slice_luma = nvcl_read_u_expgolomb(rdr);
        sps->sps_log2_diff_max_tt_min_qt_intra_slice_luma = nvcl_read_u_expgolomb(rdr);
    }

    if (sps->sps_chroma_format_idc) {
        sps->sps_qtbtt_dual_tree_intra_flag = nvcl_read_flag(rdr);
    }

    if (sps->sps_qtbtt_dual_tree_intra_flag) {
        sps->sps_log2_diff_min_qt_min_cb_intra_slice_chroma = nvcl_read_u_expgolomb(rdr);

        sps->sps_max_mtt_hierarchy_depth_intra_slice_chroma = nvcl_read_u_expgolomb(rdr);
        if (sps->sps_max_mtt_hierarchy_depth_intra_slice_chroma) {
            sps->sps_log2_diff_max_bt_min_qt_intra_slice_chroma = nvcl_read_u_expgolomb(rdr);
            sps->sps_log2_diff_max_tt_min_qt_intra_slice_chroma = nvcl_read_u_expgolomb(rdr);
        }
    }

    sps->sps_log2_diff_min_qt_min_cb_inter_slice = nvcl_read_u_expgolomb(rdr);

    sps->sps_max_mtt_hierarchy_depth_inter_slice = nvcl_read_u_expgolomb(rdr);
    if (sps->sps_max_mtt_hierarchy_depth_inter_slice) {
        sps->sps_log2_diff_max_bt_min_qt_inter_slice = nvcl_read_u_expgolomb(rdr);
        sps->sps_log2_diff_max_tt_min_qt_inter_slice = nvcl_read_u_expgolomb(rdr);
    }

    if (sps->sps_log2_ctu_size_minus5) {
        sps->sps_max_luma_transform_size_64_flag = nvcl_read_flag(rdr);
    }

    sps->sps_transform_skip_enabled_flag = nvcl_read_flag(rdr);
    if (sps->sps_transform_skip_enabled_flag) {
        sps->sps_log2_transform_skip_max_size_minus2 = nvcl_read_u_expgolomb(rdr);
        sps->sps_bdpcm_enabled_flag = nvcl_read_flag(rdr);
    }

    sps->sps_mts_enabled_flag = nvcl_read_flag(rdr);
    if (sps->sps_mts_enabled_flag) {
        sps->sps_explicit_mts_intra_enabled_flag = nvcl_read_flag(rdr);
        sps->sps_explicit_mts_inter_enabled_flag = nvcl_read_flag(rdr);
    }

    sps->sps_lfnst_enabled_flag = nvcl_read_flag(rdr);

    if (sps->sps_chroma_format_idc) {

        sps->sps_joint_cbcr_enabled_flag       = nvcl_read_flag(rdr);
        sps->sps_same_qp_table_for_chroma_flag = nvcl_read_flag(rdr);

        sps->sps_qp_table_start_minus26[0]        = nvcl_read_s_expgolomb(rdr);
        sps->sps_num_points_in_qp_table_minus1[0] = nvcl_read_u_expgolomb(rdr);
        for (j = 0; j <= (sps->sps_num_points_in_qp_table_minus1[0] & 0x3F); j++) {
            sps->sps_delta_qp_in_val_minus1[0][j] = nvcl_read_u_expgolomb(rdr);
            sps->sps_delta_qp_diff_val[0][j] = nvcl_read_u_expgolomb(rdr);
        }

        if (!sps->sps_same_qp_table_for_chroma_flag) {
            for (i = 1; i < 2 + sps->sps_joint_cbcr_enabled_flag; i++) {
                sps->sps_qp_table_start_minus26[i] = nvcl_read_s_expgolomb(rdr);
                sps->sps_num_points_in_qp_table_minus1[i] = nvcl_read_u_expgolomb(rdr);
                for (j = 0; j <= (sps->sps_num_points_in_qp_table_minus1[i] & 0x3F); j++) {
                    sps->sps_delta_qp_in_val_minus1[i][j] = nvcl_read_u_expgolomb(rdr);
                    sps->sps_delta_qp_diff_val[i][j] = nvcl_read_u_expgolomb(rdr);
                }
            }
        }
    }

    sps->sps_sao_enabled_flag = nvcl_read_flag(rdr);

    sps->sps_alf_enabled_flag = nvcl_read_flag(rdr);
    if (sps->sps_alf_enabled_flag && sps->sps_chroma_format_idc) {
        sps->sps_ccalf_enabled_flag = nvcl_read_flag(rdr);
    }

    sps->sps_lmcs_enabled_flag = nvcl_read_flag(rdr);
    sps->sps_weighted_pred_flag = nvcl_read_flag(rdr);
    sps->sps_weighted_bipred_flag = nvcl_read_flag(rdr);
    sps->sps_long_term_ref_pics_flag = nvcl_read_flag(rdr);

    if (sps->sps_video_parameter_set_id) {
        sps->sps_inter_layer_prediction_enabled_flag = nvcl_read_flag(rdr);
    }

    sps->sps_idr_rpl_present_flag = nvcl_read_flag(rdr);

    sps->sps_rpl1_same_as_rpl0_flag = nvcl_read_flag(rdr);

    sps->sps_num_ref_pic_lists0 = nvcl_read_u_expgolomb(rdr);
    if (sps->sps_num_ref_pic_lists0) {
        uint8_t nb_rpl_min1 = (sps->sps_num_ref_pic_lists0 - 1) & 0x3F;
        for (j = 0; j <= nb_rpl_min1; j++) {
            OVRPL *rpl = &sps->rpl_s0[j];
            nvcl_read_sps_ref_pic_list(rdr, sps, rpl);
        }
    }

    if (!sps->sps_rpl1_same_as_rpl0_flag) {
        sps->sps_num_ref_pic_lists1 = nvcl_read_u_expgolomb(rdr);
        if (sps->sps_num_ref_pic_lists1) {
            uint8_t nb_rpl_min1 = (sps->sps_num_ref_pic_lists1 - 1) & 0x3F;
            for (j = 0; j <= nb_rpl_min1; j++) {
                OVRPL *rpl = &sps->rpl_s1[j];
                nvcl_read_sps_ref_pic_list(rdr, sps, rpl);
            }
        }
    } else {
       /* TODO copy or ref to rpl0 ?*/
    }

    sps->sps_ref_wraparound_enabled_flag = nvcl_read_flag(rdr);

    sps->sps_temporal_mvp_enabled_flag = nvcl_read_flag(rdr);
    if (sps->sps_temporal_mvp_enabled_flag) {
        sps->sps_sbtmvp_enabled_flag = nvcl_read_flag(rdr);
    }

    sps->sps_amvr_enabled_flag = nvcl_read_flag(rdr);

    sps->sps_bdof_enabled_flag = nvcl_read_flag(rdr);
    if (sps->sps_bdof_enabled_flag) {
        sps->sps_bdof_control_present_in_ph_flag = nvcl_read_flag(rdr);
    }

    sps->sps_smvd_enabled_flag = nvcl_read_flag(rdr);

    sps->sps_dmvr_enabled_flag = nvcl_read_flag(rdr);
    if (sps->sps_dmvr_enabled_flag) {
        sps->sps_dmvr_control_present_in_ph_flag = nvcl_read_flag(rdr);
    }

    sps->sps_mmvd_enabled_flag = nvcl_read_flag(rdr);
    if (sps->sps_mmvd_enabled_flag) {
        sps->sps_mmvd_fullpel_only_enabled_flag = nvcl_read_flag(rdr);
    }

    sps->sps_six_minus_max_num_merge_cand = nvcl_read_u_expgolomb(rdr);
    sps->sps_sbt_enabled_flag = nvcl_read_flag(rdr);

    sps->sps_affine_enabled_flag = nvcl_read_flag(rdr);
    if (sps->sps_affine_enabled_flag) {
        sps->sps_five_minus_max_num_subblock_merge_cand = nvcl_read_u_expgolomb(rdr);

        sps->sps_6param_affine_enabled_flag = nvcl_read_flag(rdr);
        if (sps->sps_amvr_enabled_flag) {
            sps->sps_affine_amvr_enabled_flag = nvcl_read_flag(rdr);
        }

        sps->sps_affine_prof_enabled_flag = nvcl_read_flag(rdr);
        if (sps->sps_affine_prof_enabled_flag) {
            sps->sps_prof_control_present_in_ph_flag = nvcl_read_flag(rdr);
        }
    }

    sps->sps_bcw_enabled_flag = nvcl_read_flag(rdr);
    sps->sps_ciip_enabled_flag = nvcl_read_flag(rdr);

    if (sps->sps_six_minus_max_num_merge_cand <= 4) {
        sps->sps_gpm_enabled_flag = nvcl_read_flag(rdr);
        if (sps->sps_gpm_enabled_flag){
            /* FIXME remove dirty static initialization of GPM LUTs */
            rcn_init_gpm_params();
            if (sps->sps_six_minus_max_num_merge_cand <= 3) {
                sps->sps_max_num_merge_cand_minus_max_num_gpm_cand = nvcl_read_u_expgolomb(rdr);
            }
        }
    }

    sps->sps_log2_parallel_merge_level_minus2 = nvcl_read_u_expgolomb(rdr);

    sps->sps_isp_enabled_flag = nvcl_read_flag(rdr);
    sps->sps_mrl_enabled_flag = nvcl_read_flag(rdr);
    sps->sps_mip_enabled_flag = nvcl_read_flag(rdr);

    if (sps->sps_chroma_format_idc) {
        sps->sps_cclm_enabled_flag = nvcl_read_flag(rdr);
    }

    if (sps->sps_chroma_format_idc == 1) {
        sps->sps_chroma_horizontal_collocated_flag = nvcl_read_flag(rdr);
        sps->sps_chroma_vertical_collocated_flag   = nvcl_read_flag(rdr);
    }

    sps->sps_palette_enabled_flag = nvcl_read_flag(rdr);

    if (sps->sps_chroma_format_idc == 3 && !sps->sps_max_luma_transform_size_64_flag) {
        sps->sps_act_enabled_flag = nvcl_read_flag(rdr);
    }

    if (sps->sps_transform_skip_enabled_flag || sps->sps_palette_enabled_flag) {
        sps->sps_min_qp_prime_ts = nvcl_read_u_expgolomb(rdr);
    }

    sps->sps_ibc_enabled_flag = nvcl_read_flag(rdr);
    if (sps->sps_ibc_enabled_flag) {
        sps->sps_six_minus_max_num_ibc_merge_cand = nvcl_read_u_expgolomb(rdr);
    }

    sps->sps_ladf_enabled_flag = nvcl_read_flag(rdr);
    if (sps->sps_ladf_enabled_flag) {
        sps->sps_num_ladf_intervals_minus2 = nvcl_read_bits(rdr, 2);
        sps->sps_ladf_lowest_interval_qp_offset = nvcl_read_s_expgolomb(rdr);
        for (i = 0; i < ((sps->sps_num_ladf_intervals_minus2 + 1) & 0x3F); i++) {
            sps->sps_ladf_qp_offset[i] = nvcl_read_s_expgolomb(rdr);
            sps->sps_ladf_delta_threshold_minus1[i] = nvcl_read_u_expgolomb(rdr);
        }
    }

    sps->sps_explicit_scaling_list_enabled_flag = nvcl_read_flag(rdr);
    if (sps->sps_explicit_scaling_list_enabled_flag) {
        if (sps->sps_lfnst_enabled_flag) {
            sps->sps_scaling_matrix_for_lfnst_disabled_flag = nvcl_read_flag(rdr);
        }

        if (sps->sps_act_enabled_flag) {
            sps->sps_scaling_matrix_for_alternative_colour_space_disabled_flag = nvcl_read_flag(rdr);
            if (sps->sps_scaling_matrix_for_alternative_colour_space_disabled_flag) {
                sps->sps_scaling_matrix_designated_colour_space_flag = nvcl_read_flag(rdr);
            }
        }
    }


    sps->sps_dep_quant_enabled_flag = nvcl_read_flag(rdr);
    sps->sps_sign_data_hiding_enabled_flag = nvcl_read_flag(rdr);

    sps->sps_virtual_boundaries_enabled_flag = nvcl_read_flag(rdr);
    if (sps->sps_virtual_boundaries_enabled_flag) {
        sps->sps_virtual_boundaries_present_flag = nvcl_read_flag(rdr);
        if (sps->sps_virtual_boundaries_present_flag) {
            sps->sps_num_ver_virtual_boundaries = nvcl_read_u_expgolomb(rdr);
            for (i = 0; i < ((sps->sps_num_ver_virtual_boundaries - !!sps->sps_num_ver_virtual_boundaries) & 0x3FF); i++) {
                sps->sps_virtual_boundary_pos_x_minus1[i] = nvcl_read_u_expgolomb(rdr);
            }

            sps->sps_num_hor_virtual_boundaries = nvcl_read_u_expgolomb(rdr);
            for (i = 0; i < ((sps->sps_num_hor_virtual_boundaries - !!sps->sps_num_hor_virtual_boundaries) & 0x3FF); i++) {
                sps->sps_virtual_boundary_pos_y_minus1[i] = nvcl_read_u_expgolomb(rdr);
            }
        }
    }

    if (sps->sps_ptl_dpb_hrd_params_present_flag) {
        sps->sps_timing_hrd_params_present_flag = nvcl_read_flag(rdr);
        if (sps->sps_timing_hrd_params_present_flag) {
            struct HRDTiming hrd_timing = {0};
            general_timing_hrd_parameters(rdr, &hrd_timing);
            if (sps->sps_max_sublayers_minus1) {
                sps->sps_sublayer_cpb_params_present_flag = nvcl_read_flag(rdr);
            }

            int first_sublayer = sps->sps_sublayer_cpb_params_present_flag ? 0 : sps->sps_max_sublayers_minus1;
            ols_timing_hrd_parameters(rdr, &hrd_timing, first_sublayer, sps->sps_max_sublayers_minus1);
        }
    }

    sps->sps_field_seq_flag = nvcl_read_flag(rdr);

    sps->sps_vui_parameters_present_flag = nvcl_read_flag(rdr);
    if (sps->sps_vui_parameters_present_flag) {
        struct OVNVCLReader vui_rdr;
        sps->sps_vui_payload_size_minus1 = nvcl_read_u_expgolomb(rdr);

        nvcl_align(rdr);

        vui_rdr = *rdr;

        vui_payload(&vui_rdr, &sps->vui);

        nvcl_skip_bits(rdr, (sps->sps_vui_payload_size_minus1 + 1) << 3);
    }

    sps->sps_extension_flag = nvcl_read_flag(rdr);
    if (sps->sps_extension_flag) {
        int32_t nb_bits_read = nvcl_nb_bits_read(rdr) + 1;
        int32_t stop_bit_pos = nvcl_find_rbsp_stop_bit(rdr);
        int32_t nb_bits_remaining = stop_bit_pos - nb_bits_read;

        if (nb_bits_remaining < 0) {
            ov_log(NULL, OVLOG_ERROR, "Overread SPS %d", nb_bits_read, stop_bit_pos);
            return OVVC_EINDATA;
        }

        /* Ignore extension */
        nvcl_skip_bits(rdr, nb_bits_remaining);
    }

    return 0;
}

const struct HLSReader sps_manager =
{
    .name = "SPS",
    .data_size    = sizeof(struct OVSPS),
    .probe_id     = &probe_sps_id,
    .find_storage = &storage_in_nvcl_ctx,
    .read         = &nvcl_sps_read,
    .validate     = &validate_sps,
    //.replace      = &replace_sps,
    .free         = &free_sps
};
