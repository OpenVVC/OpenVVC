#include <stddef.h>

#include "libovvcutils/ovvcutils.h"
#include "libovvcutils/ovmem.h"

#include "nvcl.h"
#include "nvcl_utils.h"
#include "nvcl_structures.h"
#include "nvcl_private.h"

enum DecReturn {
    OV_INVALID_DATA = -1,
    OV_ENOMEM = -2,
};

static int
validate_ph(OVNVCLReader *rdr, OVPH *const ph)
{
    /* TODO various check on limitation and max sizes */
    return 1;
}

static void
free_ph(OVPH *const ph)
{
    /* TODO unref and/or free dynamic structure */
    ov_free(ph);
}

static void
replace_ph(OVNVCLCtx *const nvcl_ctx, OVPH *const ph)
{
    /* TODO unref and/or free dynamic structure */
    OVPH *to_free = nvcl_ctx->ph;

    free_ph(to_free);

    nvcl_ctx->ph = ph;
}

int
nvcl_decode_nalu_ph(OVNVCLReader *const rdr, OVNVCLCtx *const nvcl_ctx)
{
    int ret;
    /* TODO compare RBSP data to avoid new read */

    OVPH *ph = ov_mallocz(sizeof(*ph));
    if (!ph) {
        return OV_ENOMEM;
    }

    ret = nvcl_ph_read(rdr, ph, nvcl_ctx);
    if (ret < 0) {
        goto cleanup;
    }

    ret = validate_ph(rdr, ph);
    if (ret < 0) {
        goto cleanup;
    }

    /*FIXME unref instead of free */
    replace_ph(nvcl_ctx, ph);

    return 0;

cleanup:
    ov_free(ph);
    return ret;
}

int
nvcl_ph_read(OVNVCLReader *const rdr, OVPH *const ph,
             OVNVCLCtx *const nvcl_ctx)
{
    const OVPPS *pps = NULL;
    const OVSPS *sps = NULL;
    int i;

    ph->ph_gdr_or_irap_pic_flag = nvcl_read_flag(rdr);
    ph->ph_non_ref_pic_flag     = nvcl_read_flag(rdr);
    if (ph->ph_gdr_or_irap_pic_flag) {
        ph->ph_gdr_pic_flag = nvcl_read_flag(rdr);
    }

    ph->ph_inter_slice_allowed_flag = nvcl_read_flag(rdr);
    if (ph->ph_inter_slice_allowed_flag) {
        ph->ph_intra_slice_allowed_flag = nvcl_read_flag(rdr);
    }

    ph->ph_pic_parameter_set_id = nvcl_read_u_expgolomb(rdr);

    if (ph->ph_pic_parameter_set_id < OV_MAX_NUM_PPS) {
        uint8_t pps_id = ph->ph_pic_parameter_set_id & 0xF;
        pps = nvcl_ctx->pps_list[pps_id];
        if (pps) {
            /* We suppose sps_id already checked by sps reader */
            uint8_t sps_id = pps->pps_seq_parameter_set_id & 0xF;
            sps = nvcl_ctx->sps_list[sps_id];
        }
        if (!pps || !sps) {
            ov_log(NULL, 3, "SPS or PPS missing when trying to decode PH\n");
            return OV_INVALID_DATA;
        }
    }

    ph->ph_pic_order_cnt_lsb = nvcl_read_bits(rdr, sps->sps_log2_max_pic_order_cnt_lsb_minus4 + 4);

    if (ph->ph_gdr_pic_flag) {
        ph->ph_recovery_poc_cnt = nvcl_read_u_expgolomb(rdr);
    }

    /* FIXME this suppose to count 1 in sps_ph_extra_bits table */
    for (i = 0; i < sps->sps_num_extra_ph_bytes; i++) {
        ph->ph_extra_bit[i] = nvcl_read_bits(rdr, 1);
    }

    if (sps->sps_poc_msb_cycle_flag) {
        ph->ph_poc_msb_cycle_present_flag = nvcl_read_flag(rdr);
        if (ph->ph_poc_msb_cycle_present_flag) {
            ph->ph_poc_msb_cycle_val = nvcl_read_bits(rdr, sps->sps_poc_msb_cycle_len_minus1 + 1);
        }
    }

    if (sps->sps_alf_enabled_flag && pps->pps_alf_info_in_ph_flag) {
        ph->ph_alf_enabled_flag = nvcl_read_flag(rdr);
        if (ph->ph_alf_enabled_flag) {
            ph->ph_num_alf_aps_ids_luma = nvcl_read_bits(rdr, 3);
            for (i = 0; i < ph->ph_num_alf_aps_ids_luma; i++) {
                ph->ph_alf_aps_id_luma[i] = nvcl_read_bits(rdr, 3);
            }

            if (sps->sps_chroma_format_idc != 0) {
                ph->ph_alf_cb_enabled_flag = nvcl_read_flag(rdr);
                ph->ph_alf_cr_enabled_flag = nvcl_read_flag(rdr);
            }

            if (ph->ph_alf_cb_enabled_flag || ph->ph_alf_cr_enabled_flag) {
                ph->ph_alf_aps_id_chroma = nvcl_read_bits(rdr, 3);
            }

            if (sps->sps_ccalf_enabled_flag) {
                ph->ph_alf_cc_cb_enabled_flag = nvcl_read_flag(rdr);
                if (ph->ph_alf_cc_cb_enabled_flag) {
                    ph->ph_alf_cc_cb_aps_id = nvcl_read_bits(rdr, 3);
                }

                ph->ph_alf_cc_cr_enabled_flag = nvcl_read_flag(rdr);
                if (ph->ph_alf_cc_cr_enabled_flag) {
                    ph->ph_alf_cc_cr_aps_id = nvcl_read_bits(rdr, 3);
                }
            }
        }
    }

    if (sps->sps_lmcs_enabled_flag) {
        ph->ph_lmcs_enabled_flag = nvcl_read_flag(rdr);
        if (ph->ph_lmcs_enabled_flag) {
            ph->ph_lmcs_aps_id = nvcl_read_bits(rdr, 2);
            if (sps->sps_chroma_format_idc != 0) {
                ph->ph_chroma_residual_scale_flag = nvcl_read_flag(rdr);
            }
        }
    }

    if (sps->sps_explicit_scaling_list_enabled_flag) {
        ph->ph_explicit_scaling_list_enabled_flag = nvcl_read_flag(rdr);
        if (ph->ph_explicit_scaling_list_enabled_flag) {
            ph->ph_scaling_list_aps_id = nvcl_read_bits(rdr, 3);
        }
    }

    if (sps->sps_virtual_boundaries_enabled_flag && !sps->sps_virtual_boundaries_present_flag) {
        ph->ph_virtual_boundaries_present_flag = nvcl_read_flag(rdr);
        if (ph->ph_virtual_boundaries_present_flag) {
            ph->ph_num_ver_virtual_boundaries = nvcl_read_u_expgolomb(rdr);
            for (i = 0; i < ph->ph_num_ver_virtual_boundaries; i++) {
                ph->ph_virtual_boundary_pos_x_minus1[i] = nvcl_read_u_expgolomb(rdr);
            }

            ph->ph_num_hor_virtual_boundaries = nvcl_read_u_expgolomb(rdr);
            for (i = 0; i < ph->ph_num_hor_virtual_boundaries; i++) {
                ph->ph_virtual_boundary_pos_y_minus1[i] = nvcl_read_u_expgolomb(rdr);
            }
        }
    }

    if (pps->pps_output_flag_present_flag && !ph->ph_non_ref_pic_flag) {
        ph->ph_pic_output_flag = nvcl_read_flag(rdr);
    }

    if (pps->pps_rpl_info_in_ph_flag) {
         struct OVHRPL *const hrpl = &ph->hrpl;
         nvcl_read_header_ref_pic_lists(rdr, hrpl, sps, pps);
        #if 0
        ref_pic_lists();
        #endif
    }

    if (sps->sps_partition_constraints_override_enabled_flag) {
        ph->ph_partition_constraints_override_flag = nvcl_read_flag(rdr);
    }

    if (ph->ph_intra_slice_allowed_flag) {
        if (ph->ph_partition_constraints_override_flag) {
            ph->ph_log2_diff_min_qt_min_cb_intra_slice_luma = nvcl_read_u_expgolomb(rdr);
            ph->ph_max_mtt_hierarchy_depth_intra_slice_luma = nvcl_read_u_expgolomb(rdr);
            if (ph->ph_max_mtt_hierarchy_depth_intra_slice_luma != 0) {
                ph->ph_log2_diff_max_bt_min_qt_intra_slice_luma = nvcl_read_u_expgolomb(rdr);
                ph->ph_log2_diff_max_tt_min_qt_intra_slice_luma = nvcl_read_u_expgolomb(rdr);
            }

            if (sps->sps_qtbtt_dual_tree_intra_flag) {
                ph->ph_log2_diff_min_qt_min_cb_intra_slice_chroma = nvcl_read_u_expgolomb(rdr);
                ph->ph_max_mtt_hierarchy_depth_intra_slice_chroma = nvcl_read_u_expgolomb(rdr);
                if (ph->ph_max_mtt_hierarchy_depth_intra_slice_chroma != 0) {
                    ph->ph_log2_diff_max_bt_min_qt_intra_slice_chroma = nvcl_read_u_expgolomb(rdr);
                    ph->ph_log2_diff_max_tt_min_qt_intra_slice_chroma = nvcl_read_u_expgolomb(rdr);
                }
            }
        }

        if (pps->pps_cu_qp_delta_enabled_flag) {
            ph->ph_cu_qp_delta_subdiv_intra_slice = nvcl_read_u_expgolomb(rdr);
        }

        if (pps->pps_cu_chroma_qp_offset_list_enabled_flag) {
            ph->ph_cu_chroma_qp_offset_subdiv_intra_slice = nvcl_read_u_expgolomb(rdr);
        }
    }

    if (ph->ph_inter_slice_allowed_flag) {
        if (ph->ph_partition_constraints_override_flag) {
            ph->ph_log2_diff_min_qt_min_cb_inter_slice = nvcl_read_u_expgolomb(rdr);
            ph->ph_max_mtt_hierarchy_depth_inter_slice = nvcl_read_u_expgolomb(rdr);
            if (ph->ph_max_mtt_hierarchy_depth_inter_slice != 0) {
                ph->ph_log2_diff_max_bt_min_qt_inter_slice = nvcl_read_u_expgolomb(rdr);
                ph->ph_log2_diff_max_tt_min_qt_inter_slice = nvcl_read_u_expgolomb(rdr);
            }
        }

        if (pps->pps_cu_qp_delta_enabled_flag) {
            ph->ph_cu_qp_delta_subdiv_inter_slice = nvcl_read_u_expgolomb(rdr);
        }

        if (pps->pps_cu_chroma_qp_offset_list_enabled_flag) {
            ph->ph_cu_chroma_qp_offset_subdiv_inter_slice = nvcl_read_u_expgolomb(rdr);
        }

        if (sps->sps_temporal_mvp_enabled_flag) {
            ph->ph_temporal_mvp_enabled_flag = nvcl_read_flag(rdr);
            if (ph->ph_temporal_mvp_enabled_flag && pps->pps_rpl_info_in_ph_flag) {
                /* FIXME compute num_ref_entries */
                int num_ref_entries1 = 0; /* num_ref_entries[1][RplsIdx[1]] */
                if (num_ref_entries1 > 0) {
                    int num_ref_entries0 = 0; /* num_ref_entries[0][RplsIdx[0]] */
                    ph->ph_collocated_from_l0_flag = nvcl_read_flag(rdr);
                    if ((ph->ph_collocated_from_l0_flag && num_ref_entries0 > 1) || (!ph->ph_collocated_from_l0_flag && num_ref_entries1 > 1)) {
                        ph->ph_collocated_ref_idx = nvcl_read_u_expgolomb(rdr);
                    }
                }
            }
        }

        if (sps->sps_mmvd_fullpel_only_enabled_flag) {
            ph->ph_mmvd_fullpel_only_flag = nvcl_read_flag(rdr);
        }

        /* FIXME Num ref entries can be undefined if rpl info are in PH */
        int num_ref_entries1 = 0;
        if (!pps->pps_rpl_info_in_ph_flag) {

            ph->ph_mvd_l1_zero_flag = nvcl_read_flag(rdr);
            if (sps->sps_bdof_control_present_in_ph_flag) {
                ph->ph_bdof_disabled_flag = nvcl_read_flag(rdr);
            }

            if (sps->sps_dmvr_control_present_in_ph_flag) {
                ph->ph_dmvr_disabled_flag = nvcl_read_flag(rdr);
            }

        } else if (num_ref_entries1 > 0){

            /* FIXME compute num_ref_entries */
            if (num_ref_entries1 > 0) {
                ph->ph_mvd_l1_zero_flag = nvcl_read_flag(rdr);
                if (sps->sps_bdof_control_present_in_ph_flag) {
                    ph->ph_bdof_disabled_flag = nvcl_read_flag(rdr);
                }

                if (sps->sps_dmvr_control_present_in_ph_flag) {
                    ph->ph_dmvr_disabled_flag = nvcl_read_flag(rdr);
                }
            }
        }

        if (sps->sps_prof_control_present_in_ph_flag) {
            ph->ph_prof_disabled_flag = nvcl_read_flag(rdr);
        }

        if ((pps->pps_weighted_pred_flag || pps->pps_weighted_bipred_flag) && pps->pps_wp_info_in_ph_flag) {
            #if 0
            pred_weight_table()
            #endif
        }
    }

    if (pps->pps_qp_delta_info_in_ph_flag) {
        ph->ph_qp_delta = nvcl_read_s_expgolomb(rdr);
    }

    if (sps->sps_joint_cbcr_enabled_flag) {
        ph->ph_joint_cbcr_sign_flag = nvcl_read_flag(rdr);
    }

    if (sps->sps_sao_enabled_flag && pps->pps_sao_info_in_ph_flag) {
        ph->ph_sao_luma_enabled_flag = nvcl_read_flag(rdr);
        if (sps->sps_chroma_format_idc != 0) {
            ph->ph_sao_chroma_enabled_flag = nvcl_read_flag(rdr);
        }
    }

    if (pps->pps_dbf_info_in_ph_flag) {
        ph->ph_deblocking_params_present_flag = nvcl_read_flag(rdr);
        if (ph->ph_deblocking_params_present_flag) {
            if (!pps->pps_deblocking_filter_disabled_flag) {
                ph->ph_deblocking_filter_disabled_flag = nvcl_read_flag(rdr);
            }
            if (!ph->ph_deblocking_filter_disabled_flag) {
                ph->ph_luma_beta_offset_div2 = nvcl_read_s_expgolomb(rdr);
                ph->ph_luma_tc_offset_div2 = nvcl_read_s_expgolomb(rdr);
                if (pps->pps_chroma_tool_offsets_present_flag) {
                    ph->ph_cb_beta_offset_div2 = nvcl_read_s_expgolomb(rdr);
                    ph->ph_cb_tc_offset_div2 = nvcl_read_s_expgolomb(rdr);
                    ph->ph_cr_beta_offset_div2 = nvcl_read_s_expgolomb(rdr);
                    ph->ph_cr_tc_offset_div2 = nvcl_read_s_expgolomb(rdr);
                }
            }
        }
    }

    if (pps->pps_picture_header_extension_present_flag) {
        ph->ph_extension_length = nvcl_read_u_expgolomb(rdr);
        for (i = 0; i < ph->ph_extension_length; i++) {
            ph->ph_extension_data_byte[i] = nvcl_read_bits(rdr, 8);
        }
    }
    /*TODO decide on return checks and values */
    return 0;
}
