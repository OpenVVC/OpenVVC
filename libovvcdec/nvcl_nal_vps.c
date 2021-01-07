#include "nvcl.h"
#include "nvcl_utils.h"

struct OVVPS
{
    uint8_t vps_video_parameter_set_id;
    uint8_t vps_max_layers_minus1;
    uint8_t vps_max_sublayers_minus1;
    uint8_t vps_default_ptl_dpb_hrd_max_tid_flag;
    uint8_t vps_all_independent_layers_flag;
    uint8_t vps_layer_id[64];
    uint8_t vps_independent_layer_flag[64];
    uint8_t vps_max_tid_ref_present_flag[i];
    uint8_t vps_direct_ref_layer_flag[i][j];
    uint8_t vps_max_tid_il_ref_pics_plus1[i][j];
    uint8_t vps_each_layer_is_an_ols_flag;
    uint8_t vps_ols_mode_idc;
    uint8_t vps_num_output_layer_sets_minus2;
    uint8_t vps_ols_output_layer_flag[i][j];
    uint8_t vps_num_ptls_minus1;
    uint8_t vps_pt_present_flag[i];
    uint8_t vps_ptl_max_tid[i];
    uint8_t vps_ptl_alignment_zero_bit;
    OVPTL * profile_tier_level;
    uint8_t vps_ols_ptl_idx[i];
    uint8_t vps_num_dpb_params_minus1;
    uint8_t vps_sublayer_dpb_params_present_flag;
    uint8_t vps_dpb_max_tid[i];
    OVDPBParams *dpb_params
    uint8_t vps_ols_dpb_pic_width[i];
    uint8_t vps_ols_dpb_pic_height[i];
    uint8_t vps_ols_dpb_chroma_format[i];
    uint8_t vps_ols_dpb_bitdepth_minus8[i];
    uint8_t vps_ols_dpb_params_idx[i];
    uint8_t vps_timing_hrd_params_present_flag;
    OVGHRDTiming * ghrd;
    uint8_t vps_sublayer_cpb_params_present_flag;
    uint8_t vps_num_ols_timing_hrd_params_minus1;
    uint8_t vps_hrd_max_tid[i];
    OVOLSHRDTiming *ols_hrd;
    uint8_t vps_ols_timing_hrd_idx[i];
    /* unused */
    uint8_t vps_extension_flag;
    uint8_t vps_extension_data_flag;
} OVVPS;

int
nvcl_vps_read(OVNVCLReader *const rdr, OVVPS *const vps,
                  OVNVCLCtx *const nvcl_ctx)
{
    vps->vps_video_parameter_set_id = nvcl_read_bits(rdr, 4);
    vps->vps_max_layers_minus1      = nvcl_read_bits(rdr, 6);
    vps->vps_max_sublayers_minus1   = nvcl_read_bits(rdr, 3);

    if (vps->vps_max_layers_minus1 > 0 && vps->vps_max_sublayers_minus1 > 0) {
        vps->vps_default_ptl_dpb_hrd_max_tid_flag = nvcl_read_flag(rdr);
    }

    if (vps->vps_max_layers_minus1 > 0) {
        vps->vps_all_independent_layers_flag = nvcl_read_flag(rdr);
    }

    for (i = 0; i <= vps->vps_max_layers_minus1; i++) {
        vps->vps_layer_id[i] = nvcl_read_bits(rdr, 6);
        if (i > 0 && !vps->vps_all_independent_layers_flag) {
            vps->vps_independent_layer_flag[i] = nvcl_read_flag(rdr);
            if (!vps->vps_independent_layer_flag[i]) {
                vps->vps_max_tid_ref_present_flag[i] = nvcl_read_flag(rdr);
                for (j = 0; j < i; j++) {
                    vps->vps_direct_ref_layer_flag[i][j] = nvcl_read_flag(rdr);
                    if (vps->vps_max_tid_ref_present_flag[i] && vps->vps_direct_ref_layer_flag[i][j]) {
                        vps->vps_max_tid_il_ref_pics_plus1[i][j] = nvcl_read_bits(rdr, 3);
                    }
                }
            }
        }
    }

    if (vps->vps_max_layers_minus1 > 0) {
        if (vps->vps_all_independent_layers_flag) {
            vps->vps_each_layer_is_an_ols_flag = nvcl_read_flag(rdr);
        }

        if (!vps->vps_each_layer_is_an_ols_flag) {
            if (!vps->vps_all_independent_layers_flag) {
                vps->vps_ols_mode_idc = nvcl_read_bits(rdr, 2);
            }

            if (vps->vps_ols_mode_idc == 2) {
                vps->vps_num_output_layer_sets_minus2 = nvcl_read_bits(rdr, 8);
                for (i = 1; i <= vps->vps_num_output_layer_sets_minus2 + 1; i ++) {
                    for (j = 0; j <= vps->vps_max_layers_minus1; j++) {
                        vps->vps_ols_output_layer_flag[i][j] = nvcl_read_flag(rdr);
                    }
                }
            }
        }
        vps->vps_num_ptls_minus1 = nvcl_read_bits(rdr, 8);
    }

    for (i = 0; i <= vps->vps_num_ptls_minus1; i++) {
        if (i > 0) {
            vps->vps_pt_present_flag[i] = nvcl_read_flag(rdr);
        }

        if (!vps->vps_default_ptl_dpb_hrd_max_tid_flag) {
            vps->vps_ptl_max_tid[i] = nvcl_read_bits(rdr, 3);
        }
    }

    /* TODO align function */
    while(!byte_aligned()) {
        vps->vps_ptl_alignment_zero_bit = nvcl_read_bits(rdr, 1);
    }

    for (i = 0; i <= vps->vps_num_ptls_minus1; i++) {
        profile_tier_level(vps->vps_pt_present_flag[i], vps->vps_ptl_max_tid[i]);
    }

    if (vps->vps_num_ptls_minus1 > 0 && vps->vps_num_ptls_minus1 + 1 != TotalNumOlss) {
        for (i = 0; i < TotalNumOlss; i++) {
            vps->vps_ols_ptl_idx[i] = nvcl_read_bits(rdr, 8);
        }
    }

    if (!vps->vps_each_layer_is_an_ols_flag) {
        vps->vps_num_dpb_params_minus1 = nvcl_read_u_expgolomb(rdr);
        if (vps->vps_max_sublayers_minus1 > 0) {
            vps->vps_sublayer_dpb_params_present_flag = nvcl_read_flag(rdr);
        }

        for (i = 0; i < VpsNumDpbParams; i++) {
            if (!vps->vps_default_ptl_dpb_hrd_max_tid_flag) {
                vps->vps_dpb_max_tid[i] = nvcl_read_bits(rdr, 3);
            }
            dpb_parameters(vps->vps_dpb_max_tid[i], vps->vps_sublayer_dpb_params_present_flag);
        }

        for (i = 0; i < NumMultiLayerOlss; i++) {
            vps->vps_ols_dpb_pic_width[i]       = nvcl_read_u_expgolomb(rdr);
            vps->vps_ols_dpb_pic_height[i]      = nvcl_read_u_expgolomb(rdr);
            vps->vps_ols_dpb_chroma_format[i]   = nvcl_read_bits(rdr, 2);
            vps->vps_ols_dpb_bitdepth_minus8[i] = nvcl_read_u_expgolomb(rdr);

            if (VpsNumDpbParams > 1 && VpsNumDpbParams != NumMultiLayerOlss) {
                vps->vps_ols_dpb_params_idx[i] = nvcl_read_u_expgolomb(rdr);
            }
        }

        vps->vps_timing_hrd_params_present_flag = nvcl_read_flag(rdr);
        if (vps->vps_timing_hrd_params_present_flag) {
            general_timing_hrd_parameters();
            if (vps->vps_max_sublayers_minus1 > 0) {
                    vps->vps_sublayer_cpb_params_present_flag = nvcl_read_flag(rdr);
            }

            vps->vps_num_ols_timing_hrd_params_minus1 = nvcl_read_u_expgolomb(rdr);
            for (i = 0; i <= vps->vps_num_ols_timing_hrd_params_minus1; i++) {
                if (!vps->vps_default_ptl_dpb_hrd_max_tid_flag) {
                    vps->vps_hrd_max_tid[i] = nvcl_read_flag(rdr, 3);
                }
                int firstSubLayer = vps->vps_sublayer_cpb_params_present_flag ? 0 : vps->vps_hrd_max_tid[i];
                ols_timing_hrd_parameters(firstSubLayer, vps->vps_hrd_max_tid[i]);
            }

            if (vps->vps_num_ols_timing_hrd_params_minus1 > 0 && vps->vps_num_ols_timing_hrd_params_minus1 + 1 != NumMultiLayerOlss) {
                for (i = 0; i < NumMultiLayerOlss; i++) {
                    vps->vps_ols_timing_hrd_idx[i] = nvcl_read_u_expgolomb(rdr);
                }
            }
        }
    }

    vps->vps_extension_flag = nvcl_read_flag(rdr);
    if (vps->vps_extension_flag) {
        while(more_rbsp_data()) {
            vps->vps_extension_data_flag = nvcl_read_flag(rdr);
        }
    }

    rbsp_trailing_bits()
}
