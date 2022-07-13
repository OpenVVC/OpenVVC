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
#include "overror.h"

#include "hls_structures.h"
#include "nvcl.h"
#include "nvcl_utils.h"
#include "nvcl_structures.h"
#include "nvcl_private.h"

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

static uint8_t
probe_vps_id(OVNVCLReader *const rdr)
{
    uint8_t vps_id = fetch_bits(rdr, 4);
    return vps_id;
}

static struct HLSDataRef **
storage_in_nvcl_ctx(OVNVCLReader *const rdr, OVNVCLCtx *const nvcl_ctx)
{
    uint8_t id = probe_vps_id(rdr);
    struct HLSDataRef **storage = &nvcl_ctx->vps_list[id];

    return storage;
}

static int
validate_vps(OVNVCLReader *rdr, const union HLSData *const data)
{
    /* TODO various check on limitation and max sizes */
    const OVVPS *const vps =  (const OVVPS *)data;
    uint32_t nb_bits_read = nvcl_nb_bits_read(rdr) + 1;
    uint32_t stop_bit_pos = nvcl_find_rbsp_stop_bit(rdr);
    if (stop_bit_pos != nb_bits_read) {

        ov_log(NULL, OVLOG_ERROR, "rbsp_stop_bit mismatch: cursor at %d,  expected %d\n", nb_bits_read, stop_bit_pos);
        return OVVC_EINDATA;
    }


    return 1;
}

static void
free_vps(const union HLSData *const data)
{
    /* TODO unref and/or free dynamic structure */
    const OVVPS *const vps = (const OVVPS *)data;
    ov_free((void *)vps);
}

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
            hrd->bit_rate_scale = nvcl_read_bits(rdr, 4);
            hrd->cpb_size_scale = nvcl_read_bits(rdr, 4);
        }

        if (hrd->general_du_hrd_params_present_flag) {
            hrd->cpb_size_du_scale  = nvcl_read_bits(rdr,4);
            hrd->hrd_cpb_cnt_minus1 = nvcl_read_u_expgolomb(rdr);
        }
    }
}

static void
sublayer_hrd_parameters(OVNVCLReader *const rdr, const struct HRDTiming *const hrd, uint8_t i)
{
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
ols_timing_hrd_parameters(OVNVCLReader *const rdr, const struct HRDTiming *const hrd, uint8_t first_sublayer, uint8_t max_sublayer)
{
    int i;
    for (i = first_sublayer; i <= max_sublayer; ++i) {
        uint8_t fixed_pic_rate_general_flag = nvcl_read_flag(rdr);
        uint8_t fixed_pic_rate_within_cvs_flag = 0;;
        if (fixed_pic_rate_general_flag) {
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
#if 0
{
    num_layer_in_ols[0] = 1;
    layer_id_in_ols[0][0] = vps_layer_id[0];
    num_multi_layer_olss = 0;

    if (vps_each_layer_is_an_ols_flag) {
        /* Note this is equivalent to ols_mode_idc == 4 */
        total_num_olss = vps_max_layers_minus1 + 1;
        for (i = 1; i < total_num_olss; ++i) {
            num_layer_in_ols[i] = 1;
            layer_id_in_ols[i][0] = vps_layer_id[i];
        }
    } else if (vps_ols_mode_idc == 0 || vps_ols_mode_idc == 1) {
        total_num_olss = vps_max_layers_minus1 + 1;
        for (i = 1; i < total_num_olss; ++i) {
            num_layer_in_ols[i] = i + 1;
            for (j = 0; j < num_layer_in_ols[i]; ++j) {
                layer_id_in_ols[i][j] = vps_layer_id[j];
            }
            multi_layer_ols_idx[i] = i - 1;
        }
        num_multi_layer_olss = total_num_olss - 1;
    } else if (vps_ols_mode_idc == 2) {
        /* Multi layer */
        total_num_olss = vps_num_output_layer_sets_minus2 + 2;
        for (i = 1; i < total_num_olss; ++i) {
            for (k = 0, j = 0; k <= vps_max_layers_minus1; k++) {
                if (vps_ols_output_layer_flag[i][k]) {
                    layer_id_in_ols[i][j] = vps_layer_id[k];
                    j++;
                }
            }
            num_layer_in_ols[i] = j;
            if (j > 1) {
                multi_layer_ols_idx[i] = num_multi_layer_olss;
                num_multi_layer_olss++;
            }
        }
    }
}

{
    num_output_layer_in_ols[0] = 1;
    output_layer_id_in_ols[0][0] = vps_layer_id[0];
    num_sublayers_in_layer_in_ols[0][0] = vps_ptl_max_tid[vps_ols_ptl_idx[0]] + 1;
    layer_used_as_output_layer_flag[0] = 1;

    if (ols_mode_idc == 4 || ols_mode_idc == 0) {
        total_num_olss = vps_max_layers_minus1 + 1;
        for (i = 1; i <= vps_max_layers_minus1; i++) {
            layer_used_as_output_layer_flag[i] = 1;
        }
        for (i = 1; i < total_num_olss; i++) {
            num_output_layer_in_ols[i] = 1;
            output_layer_id_in_ols[i][0] = vps_layer_id[i];
            if (vps_each_layer_is_an_ols_flag) {
                num_sublayers_in_layer_in_ols[i][0] = vps_ptl_max_tid[vps_ols_ptl_idx[i]] + 1;
            } else {
                num_sublayers_in_layer_in_ols[i][i] = vps_ptl_max_tid[vps_ols_ptl_idx[i]] + 1;
                for (k = i - 1; k >= 0; k--) {
                    num_sublayers_in_layer_in_ols[i][k] = 0;
                    for (m = k + 1; m <= i; m++) {
                        max_sub_layer_needed = Min(num_sublayers_in_layer_in_ols[i][m],
                                                 vps_max_tid_il_ref_pics_plus1[m][k]);
                        if (vps_direct_ref_layer_flag[m][k] && num_sublayers_in_layer_in_ols[i][k] < max_sub_layer_needed) {
                            num_sublayers_in_layer_in_ols[i][k] = max_sub_layer_needed;
                        }
                    }
                }
            }
        }
    } else if (vps_ols_mode_idc == 1) {
        /* All layers are used as output */
        total_num_olss = vps_max_layers_minus1 + 1;
        for (i = 1; i <= vps_max_layers_minus1; i++) {
            layer_used_as_output_layer_flag[i] = 1;
        }
        for (i = 1; i < total_num_olss; i++) {
            num_output_layer_in_ols[i] = i + 1;
            for (j = 0; j < num_output_layer_in_ols[i]; j++) {
                output_layer_id_in_ols[i][j] = vps_layer_id[j];
                num_sublayers_in_layer_in_ols[i][j] = vps_ptl_max_tid[vps_ols_ptl_idx[i]] + 1;
            }
        }
    } else if (vps_ols_mode_idc == 2) {
        /* Multi layers */
        total_num_olss = vps_num_output_layer_sets_minus2 + 2;
        for (i = 1; i < total_num_olss; i++) {
            highest_included_layer = 0;
            for (k = 0, j = 0; k <= vps_max_layers_minus1; k++) {
                if (vps_ols_output_layer_flag[i][k]) {
                    highest_included_layer = k;
                    output_layer_idx[i][j] = k;
                    output_layer_id_in_ols[i][j] = vps_layer_id[k];
                    num_sublayers_in_layer_in_ols[i][k] = vps_ptl_max_tid[vps_ols_ptl_idx[i]] + 1;
                    j++;
                } else {
                    num_sublayers_in_layer_in_ols[i][j] = 0;
                }
                layer_included_in_ols_flag[i][k]   = vps_ols_output_layer_flag[i][k];
                layer_used_as_output_layer_flag[k] = vps_ols_output_layer_flag[i][k];
            }
            num_output_layer_in_ols[i] = j;

            for (j = 0; j < num_output_layer_in_ols[i]; j++) {
                idx = output_layer_idx[i][j];
                for (k = 0; k < num_ref_layers[idx]; k++) {
                    if (!layer_included_in_ols_flag[i][reference_layer_idx[idx][k]]) {
                        layer_included_in_ols_flag[i][reference_layer_idx[idx][k]] = 1;
                    }
                }
            }

            for (k = highest_included_layer - 1; k >= 0; k--) {
                if (layer_included_in_ols_flag[i][k] && !vps_ols_output_layer_flag[i][k]) {
                    for (m = k + 1; m <= highest_included_layer; m++) {
                        max_sub_layer_needed = Min(num_sublayers_in_layer_in_ols[i][m], vps_max_tid_il_ref_pics_plus1[m][k]);
                        if (vps_direct_ref_layer_flag[m][k] && layer_included_in_ols_flag[i][m] && num_sublayers_in_layer_in_ols[i][k] < max_sub_layer_needed) {
                            num_sublayers_in_layer_in_ols[i][k] = max_sub_layer_needed;
                        }
                    }
                }
            }
        }
    }
}

{
    /* set inter layer dependencies flags */
    for (i = 0; i <= vps_max_layers_minus1; i++) {
        for (j = 0; j <= vps_max_layers_minus1; j++) {
            dependency_flag[i][j] = vps_direct_ref_layer_flag[i][j];
            for (k = 0; k < i; k++) {
                if (vps_direct_ref_layer_flag[i][k] && dependency_flag[k][j]) {
                    dependency_flag[i][j] = 1;
                }
            }
        }
    }

    for (i = 0; i <= vps_max_layers_minus1; i++) {
        for (j = 0, d = 0, r = 0; j <= vps_max_layers_minus1; j++) {
            if (vps_direct_ref_layer_flag[i][j]) {
                direct_ref_layer_idx[i][d] = j;
                layer_used_as_ref_flag[j] = 1;
                d++;
            } else {
                layer_used_as_ref_flag[j] = 0;
            }

            if (dependency_flag[i][j]){
                reference_layer_idx[i][r] = j;
                r++;
            }
        }
        num_direct_ref_layers[i] = d;
        num_ref_layers[i] = r;
    }
}
#endif

int
nvcl_vps_read(OVNVCLReader *const rdr, OVHLSData *const hls_data,
              const OVNVCLCtx *const nvcl_ctx, uint8_t nalu_type)
{

    OVVPS *const vps = &hls_data->vps;
    vps->vps_video_parameter_set_id = nvcl_read_bits(rdr, 4);
    vps->vps_max_layers_minus1      = nvcl_read_bits(rdr, 6);
    vps->vps_max_sublayers_minus1   = nvcl_read_bits(rdr, 3);

    vps->vps_default_ptl_dpb_hrd_max_tid_flag = 1;
    if (vps->vps_max_layers_minus1 > 0 && vps->vps_max_sublayers_minus1 > 0) {
        vps->vps_default_ptl_dpb_hrd_max_tid_flag = nvcl_read_flag(rdr);
    }

    vps->vps_all_independent_layers_flag = 1;
    if (vps->vps_max_layers_minus1 > 0) {
        vps->vps_all_independent_layers_flag = nvcl_read_flag(rdr);
    }

    for (int i = 0; i <= vps->vps_max_layers_minus1; i++) {
        vps->vps_layer_id[i] = nvcl_read_bits(rdr, 6);
        if (i > 0 && !vps->vps_all_independent_layers_flag) {
            vps->vps_independent_layer_flag[i] = nvcl_read_flag(rdr);
            if (!vps->vps_independent_layer_flag[i]) {
                vps->vps_max_tid_ref_present_flag[i] = nvcl_read_flag(rdr);
                if (vps->vps_max_tid_ref_present_flag[i]) {
                    for (int j = 0; j < i; j++) {
                        vps->vps_direct_ref_layer_flag[i][j] = nvcl_read_flag(rdr);
                        if (vps->vps_direct_ref_layer_flag[i][j]) {
                            vps->vps_max_tid_il_ref_pics_plus1[i][j] = nvcl_read_bits(rdr, 3);
                        }
                    }
                } else {
                    for (int j = 0; j < i; j++) {
                        vps->vps_direct_ref_layer_flag[i][j] = nvcl_read_flag(rdr);
                    }
                }
            }
        }
    }

    if (vps->vps_max_layers_minus1) {
        if (vps->vps_all_independent_layers_flag) {
            vps->vps_each_layer_is_an_ols_flag = nvcl_read_flag(rdr);
        }

        if (!vps->vps_each_layer_is_an_ols_flag) {
            if (!vps->vps_all_independent_layers_flag) {
                vps->vps_ols_mode_idc = nvcl_read_bits(rdr, 2);
            }

            if (vps->vps_ols_mode_idc == 2) {
                /* unsupported multilayer */
                return -1;
                vps->vps_num_output_layer_sets_minus2 = nvcl_read_bits(rdr, 8);
                for (int i = 1; i <= vps->vps_num_output_layer_sets_minus2 + 1; i ++) {
                    for (int j = 0; j <= vps->vps_max_layers_minus1; j++) {
                        vps->vps_ols_output_layer_flag[i][j] = nvcl_read_flag(rdr);
                    }
                }
            }
        }
        vps->vps_num_ptls_minus1 = nvcl_read_bits(rdr, 8);
    } else {
        vps->vps_each_layer_is_an_ols_flag = 1;
    }

    if (!vps->vps_default_ptl_dpb_hrd_max_tid_flag) {
        vps->vps_pt_present_flag[0] = 1;
        vps->vps_ptl_max_tid[0] = nvcl_read_bits(rdr, 3);
        for (int i = 1; i <= vps->vps_num_ptls_minus1; i++) {
            vps->vps_pt_present_flag[i] = nvcl_read_flag(rdr);
            vps->vps_ptl_max_tid[i] = nvcl_read_bits(rdr, 3);
        }
    } else {
        vps->vps_pt_present_flag[0] = 1;
        vps->vps_ptl_max_tid[0] = vps->vps_max_sublayers_minus1;
        for (int i = 1; i <= vps->vps_num_ptls_minus1; i++) {
            vps->vps_pt_present_flag[i] = nvcl_read_flag(rdr);
            vps->vps_ptl_max_tid[i] = vps->vps_max_sublayers_minus1;
        }
    }

    nvcl_align(rdr);
    //while(!byte_aligned()) {
    //    vps->vps_ptl_alignment_zero_bit = nvcl_read_bits(rdr, 1);
    //}

    for (int i = 0; i <= vps->vps_num_ptls_minus1; i++) {
        profile_tier_level_vps(rdr, vps->vps_pt_present_flag[i], vps->vps_ptl_max_tid[i]);
    }

    /* We assume multi layer is  disabled */
    /* Note idc_mode == 2 should be equivalent in this case */
    uint16_t total_num_olss = vps->vps_max_layers_minus1 + 1;
    if (vps->vps_num_ptls_minus1 > 0 && vps->vps_num_ptls_minus1 + 1 != total_num_olss) {
        for (int i = 0; i < total_num_olss; i++) {
            vps->vps_ols_ptl_idx[i] = nvcl_read_bits(rdr, 8);
        }
    }

    if (!vps->vps_each_layer_is_an_ols_flag) {
        vps->vps_num_dpb_params_minus1 = nvcl_read_u_expgolomb(rdr);
        if (vps->vps_max_sublayers_minus1 > 0) {
            vps->vps_sublayer_dpb_params_present_flag = nvcl_read_flag(rdr);
        }

        if (!vps->vps_default_ptl_dpb_hrd_max_tid_flag) {
            for (int i = 0; i < vps->vps_num_dpb_params_minus1; i++) {
                vps->vps_dpb_max_tid[i] = nvcl_read_bits(rdr, 3);
                dpb_parameters(rdr, vps->dpb_parameters, vps->vps_dpb_max_tid[i], vps->vps_sublayer_dpb_params_present_flag);
            }
        } else {
            for (int i = 0; i < vps->vps_num_dpb_params_minus1; i++) {
                dpb_parameters(rdr, vps->dpb_parameters, vps->vps_dpb_max_tid[0], vps->vps_sublayer_dpb_params_present_flag);
            }
        }

        /* We assume multi layer_is _disabled */
        uint16_t num_multi_layer_olss = total_num_olss - 1;
        for (int i = 0; i < num_multi_layer_olss; i++) {
            vps->vps_ols_dpb_pic_width[i]       = nvcl_read_u_expgolomb(rdr);
            vps->vps_ols_dpb_pic_height[i]      = nvcl_read_u_expgolomb(rdr);
            vps->vps_ols_dpb_chroma_format[i]   = nvcl_read_bits(rdr, 2);
            vps->vps_ols_dpb_bitdepth_minus8[i] = nvcl_read_u_expgolomb(rdr);

            if (vps->vps_num_dpb_params_minus1 && vps->vps_num_dpb_params_minus1 - 1 != num_multi_layer_olss) {
                vps->vps_ols_dpb_params_idx[i] = nvcl_read_u_expgolomb(rdr);
            }
        }

        vps->vps_timing_hrd_params_present_flag = nvcl_read_flag(rdr);
        if (vps->vps_timing_hrd_params_present_flag) {
            struct HRDTiming hrd_timing = {0};
            general_timing_hrd_parameters(rdr, &hrd_timing);
            if (vps->vps_max_sublayers_minus1) {
                vps->vps_sublayer_cpb_params_present_flag = nvcl_read_flag(rdr);
            }

            vps->vps_num_ols_timing_hrd_params_minus1 = nvcl_read_u_expgolomb(rdr);
            if (!vps->vps_default_ptl_dpb_hrd_max_tid_flag) {
                for (int i = 0; i <= vps->vps_num_ols_timing_hrd_params_minus1; i++) {
                    vps->vps_hrd_max_tid[i] = nvcl_read_bits(rdr, 3);
                    int first_sublayer = vps->vps_sublayer_cpb_params_present_flag ? 0 : vps->vps_hrd_max_tid[i];
                    ols_timing_hrd_parameters(rdr, &hrd_timing, first_sublayer, vps->vps_hrd_max_tid[i]);
                }
            } else {
                for (int i = 0; i <= vps->vps_num_ols_timing_hrd_params_minus1; i++) {
                    int first_sublayer = vps->vps_sublayer_cpb_params_present_flag ? 0 : vps->vps_hrd_max_tid[i];
                    ols_timing_hrd_parameters(rdr, &hrd_timing, first_sublayer, vps->vps_hrd_max_tid[0]);
                }
            }

            if (vps->vps_num_ols_timing_hrd_params_minus1 && vps->vps_num_ols_timing_hrd_params_minus1 + 1 != num_multi_layer_olss) {
                for (int i = 0; i < num_multi_layer_olss; i++) {
                    vps->vps_ols_timing_hrd_idx[i] = nvcl_read_u_expgolomb(rdr);
                }
            }
        }
    }

    vps->vps_extension_flag = nvcl_read_flag(rdr);

    if (vps->vps_extension_flag) {
        int32_t nb_bits_read = nvcl_nb_bits_read(rdr) + 1;
        int32_t stop_bit_pos = nvcl_find_rbsp_stop_bit(rdr);
        int32_t nb_bits_remaining = stop_bit_pos - nb_bits_read;

        if (nb_bits_remaining < 0) {
            ov_log(NULL, OVLOG_ERROR, "Overread PPS %d", nb_bits_read, stop_bit_pos);
            return OVVC_EINDATA;
        }

        /* Ignore extension */
        nvcl_skip_bits(rdr, nb_bits_remaining);
    }

    return 0;
}

const struct HLSReader vps_manager =
{
    .name = "VPS",
    .data_size    = sizeof(struct OVVPS),
    .probe_id     = &probe_vps_id,
    .find_storage = &storage_in_nvcl_ctx,
    .read         = &nvcl_vps_read,
    .validate     = &validate_vps,
    //.replace      = &replace_sps,
    .free         = &free_vps
};
