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
#include <stddef.h>
#include "overror.h"

#include "ovutils.h"
#include "ovmem.h"
#include "ovunits.h"

#include "nvcl.h"
#include "nvcl_utils.h"
#include "nvcl_structures.h"
#include "nvcl_private.h"
#include "hls_structures.h"

enum SliceType {
     B = 0,
     P = 1,
     I = 2
};

static int
validate_sh(OVNVCLReader *rdr, const union HLSData *const data)
{
    const OVSH *const sh = (const OVSH *)&data->sh;
    if (sh->sh_slice_type > 3) return OVVC_EINDATA;

    return 0;
}

static void
free_sh(const union HLSData *const data)
{
    const OVSH *const sh = (const OVSH *)&data->sh;
    ov_free((void *)sh);
}

static struct HLSDataRef **
storage_in_nvcl_ctx(OVNVCLReader *const rdr, OVNVCLCtx *const nvcl_ctx)
{
    struct HLSDataRef **storage = &nvcl_ctx->sh;

    return storage;
}

static void
hls_replace_ref(const struct HLSReader *const hls_hdl, struct HLSDataRef **storage, const union HLSData *const data)
{
    union HLSData *tmp = ov_malloc(hls_hdl->data_size);
    if (!tmp) {
        return;
    }
    memcpy(tmp, data, hls_hdl->data_size);

    hlsdata_unref(storage);
    *storage = hlsdataref_create(tmp, NULL, NULL);

    if (!*storage) {
        return;
    }
}

const struct HLSReader sh_manager =
{
    .name = "Slice Header",
    .data_size    = sizeof(struct OVSH),
    .probe_id     = NULL,
    .find_storage = &storage_in_nvcl_ctx,
    .read         = &nvcl_sh_read,
    .validate     = &validate_sh,
    .replace      = NULL,
    .free         = &free_sh
};

/* Slice header is only presnt in VCL NAL Units but we use NVCL
 * denomination since it can be read by same reader
 */
int
nvcl_decode_nalu_sh(OVNVCLReader *const rdr, OVNVCLCtx *const nvcl_ctx, uint8_t nalu_type)
{
    #if 0
    int ret;

    OVSH *sh = ov_mallocz(sizeof(*sh));
    if (!sh) {
        return OVVC_ENOMEM;
    }

    ret = nvcl_sh_read(rdr, (OVHLSData *const)sh, nvcl_ctx, nalu_type);
    if (ret < 0) {
        goto cleanup;
    }

    #if 0
    ret = validate_sps(rdr, sh);
    if (ret < 0) {
        goto cleanup;
    }
    #endif

    /* FIXME this is temporary to ensure the
     * sh is freed between each slice
     */
    if (nvcl_ctx->sh) {
        ov_freep(&nvcl_ctx->sh);
    }

    nvcl_ctx->sh = hlsdataref_create((OVHLSData *const)sh, NULL, NULL);

    return 0;

cleanup:
    ov_free(sh);
    return ret;
    #endif
    struct HLSDataRef **storage = sh_manager.find_storage(rdr, nvcl_ctx);
    union HLSData data;
    int ret;

    if (*storage) {
        /* TODO compare RBSP data to avoid new read */
        uint8_t identical_rbsp = 0;
        if (identical_rbsp) goto duplicated;
    }

    memset(&data, 0, sh_manager.data_size);
    ov_log(NULL, OVLOG_TRACE, "Reading new %s\n", sh_manager.name);

    ret = sh_manager.read(rdr, &data, nvcl_ctx, nalu_type);
    if (ret < 0)  goto failread;

    #if 0
    ov_log(NULL, OVLOG_TRACE, "Checking %s\n", sh_manager.name);
    ret = sh_manager->validate(rdr, &data);
    if (ret < 0)  goto invalid;

    ov_log(NULL, OVLOG_TRACE, "Replacing new %s\n", sh_manager.name);
    ret = hls_hdl->replace(hls_hdl, storage, &data);
    #else
    ov_log(NULL, OVLOG_TRACE, "Replacing new %s\n", sh_manager.name);
    hls_replace_ref(&sh_manager, storage, &data);
    #endif

    return ret;

invalid:
    ov_log(NULL, OVLOG_ERROR, "Invalid %s\n", sh_manager.name);
    return ret;

failread:
    ov_log(NULL, OVLOG_ERROR, "Error while reading %s\n", sh_manager.name);
    return ret;

duplicated:
    ov_log(NULL, OVLOG_TRACE, "Ignored duplicated %s\n", sh_manager.name);
    return 0;

}

extern const struct HLSReader ph_manager;

static int
nvcl_decode_ph(OVNVCLReader *const rdr, OVNVCLCtx *const nvcl_ctx,
               const struct HLSReader *const hls_hdl, uint8_t nalu_type)
{
    struct HLSDataRef **storage = hls_hdl->find_storage(rdr, nvcl_ctx);
    union HLSData data;
    int ret;

    /* FIXME PH cannot be identical since it cntains POC info */
    if (*storage) {
        /* TODO compare RBSP data to avoid new read */
        uint8_t identical_rbsp = 0;
        if (identical_rbsp) goto duplicated;
    }

    memset(&data, 0, hls_hdl->data_size);

    ov_log(NULL, OVLOG_TRACE, "Reading new %s\n", hls_hdl->name);

    ret = hls_hdl->read(rdr, &data, nvcl_ctx, nalu_type);
    if (ret < 0)  goto failread;

    ret = hls_hdl->validate(rdr, &data);
    if (ret < 0)  goto invalid;

    #if 0
    ret = hls_hdl->replace(hls_hdl, storage, &data);
    #else
    hls_replace_ref(hls_hdl, storage, &data);
    #endif

    return ret;

invalid:
    ov_log(NULL, OVLOG_ERROR, "Invalid %s\n", hls_hdl->name);
    return ret;

failread:
    ov_log(NULL, OVLOG_ERROR, "Error while reading %s\n", hls_hdl->name);
    return ret;

duplicated:
    ov_log(NULL, OVLOG_TRACE, "Ignored duplicated %s\n", hls_hdl->name);
    return 0;
}

int
nvcl_sh_read(OVNVCLReader *const rdr, OVHLSData *const hls_data,
             const OVNVCLCtx *const nvcl_ctx, uint8_t nalu_type)
{
    int i;
    OVSH *sh = &hls_data->sh;
    OVPH *ph = NULL;
    const OVPPS *pps = NULL;
    const OVSPS *sps = NULL;

    sh->sh_picture_header_in_slice_header_flag = nvcl_read_flag(rdr);
    if (sh->sh_picture_header_in_slice_header_flag) {
        /* FIXME do not skip 16 first bits in this case */
        int ret = nvcl_decode_ph(rdr, nvcl_ctx, &ph_manager, 0);
        if (ret < 0) {
            ov_log(NULL, 3, "Failed reading PH from SH\n");
            return OVVC_EINDATA;
        }
    }

    /* TODO create proper structure to hold activated parameters sets */
    if (nvcl_ctx->ph) {
        ph = (OVPH *)nvcl_ctx->ph->data;
        pps = (OVPPS *)nvcl_ctx->pps_list[ph->ph_pic_parameter_set_id]->data;
        sps = (OVSPS *)nvcl_ctx->sps_list[pps->pps_seq_parameter_set_id]->data;
    }

    if (!ph || !pps || !sps) {
        ov_log(NULL, 3, "Missing parameter sets while reading SH\n");
        return OVVC_EINDATA;
    }
    /* TODO Once other parameter sets are properly activated we can derive
     * convenience variables from them including :
     *    - Sub Pictures / Tiles contexts
     *    - Numbers of extra sh bits
     */

    if (sps->sps_subpic_info_present_flag) {
        /*FIXME check subpic_id_len_minus1 overrides */
        sh->sh_subpic_id = nvcl_read_bits(rdr, sps->sps_subpic_id_len_minus1 + 1);
    }

    /* FIXME subpic/tile  support */
    int nb_slices_subpic = 1; /*NumSlicesInSubpic[CurrSubpicIdx] */
    int nb_tiles_pic = (pps->pps_num_tile_columns_minus1 + 1) * (pps->pps_num_tile_rows_minus1 + 1);

    /* FIXME would be better to distinguish in two different branches */
    if ((pps->pps_rect_slice_flag && nb_slices_subpic > 1) ||
        (!pps->pps_rect_slice_flag && nb_tiles_pic > 1)) {
        /*TODO ceil log2_num_tiles_in_pic / num_slices_in_sub_pc*/
        int nb_bits_in_slice_address = 0;
        sh->sh_slice_address = nvcl_read_bits(rdr, nb_bits_in_slice_address);
        if (sh->sh_slice_address) {
           ov_log(NULL, OVLOG_ERROR, "Invalid slice address\n");
           return -1;
        }
    }

    /* FIXME to be derived from sps */
    int nb_extra_sh_bits = 0;
    for (i = 0; i < nb_extra_sh_bits; i++) {
        sh->sh_extra_bit[i] = nvcl_read_bits(rdr, 1);
    }

    if (!pps->pps_rect_slice_flag && nb_tiles_pic - sh->sh_slice_address > 1) {
        sh->sh_num_tiles_in_slice_minus1 = nvcl_read_u_expgolomb(rdr);
    }

    sh->sh_slice_type = I;
    if (ph->ph_inter_slice_allowed_flag) {
        sh->sh_slice_type = nvcl_read_u_expgolomb(rdr);
        if (sh->sh_slice_type > 3) return OVVC_EINDATA;
    }

    if (nalu_type == OVNALU_IDR_W_RADL ||
        nalu_type == OVNALU_IDR_N_LP   ||
        nalu_type == OVNALU_CRA        ||
        nalu_type == OVNALU_GDR) {
        sh->sh_no_output_of_prior_pics_flag = nvcl_read_flag(rdr);
    }

    if (sps->sps_alf_enabled_flag && !pps->pps_alf_info_in_ph_flag) {

        sh->sh_alf_enabled_flag = nvcl_read_flag(rdr);
        if (sh->sh_alf_enabled_flag) {

            sh->sh_num_alf_aps_ids_luma = nvcl_read_bits(rdr, 3);
            for (i = 0; i < sh->sh_num_alf_aps_ids_luma; i++) {
                sh->sh_alf_aps_id_luma[i] = nvcl_read_bits(rdr, 3);
            }

            if (sps->sps_chroma_format_idc) {
                sh->sh_alf_cb_enabled_flag = nvcl_read_flag(rdr);
                sh->sh_alf_cr_enabled_flag = nvcl_read_flag(rdr);

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
    }

    sh->sh_lmcs_used_flag = ph->ph_lmcs_enabled_flag && sh->sh_picture_header_in_slice_header_flag;
    if (ph->ph_lmcs_enabled_flag && !sh->sh_picture_header_in_slice_header_flag) {
        sh->sh_lmcs_used_flag = nvcl_read_flag(rdr);
    }

    if (ph->ph_explicit_scaling_list_enabled_flag && !sh->sh_picture_header_in_slice_header_flag) {
        sh->sh_explicit_scaling_list_used_flag = nvcl_read_flag(rdr);
    }

    if (!pps->pps_rpl_info_in_ph_flag &&
        ((nalu_type != OVNALU_IDR_W_RADL && nalu_type != OVNALU_IDR_N_LP) ||
         sps->sps_idr_rpl_present_flag)) {
         OVHRPL *const hrpl = &sh->hrpl;
         nvcl_read_header_ref_pic_lists(rdr, hrpl, sps, pps);
    }

    if (sh->sh_slice_type != I) {
        const OVHRPL *const hrpl = pps->pps_rpl_info_in_ph_flag ? &ph->hrpl : &sh->hrpl;

        int nb_ref_entries0 = hrpl->rpl0 ? hrpl->rpl_h0.rpl_data.num_ref_entries : 0;
        int nb_ref_entries1 = hrpl->rpl1 ? hrpl->rpl_h1.rpl_data.num_ref_entries : 0;

        if ((nb_ref_entries0 > 1) ||
            (sh->sh_slice_type == B && nb_ref_entries1 > 1)) {
            sh->sh_num_ref_idx_active_override_flag = nvcl_read_flag(rdr);
            if (sh->sh_num_ref_idx_active_override_flag) {
                if (nb_ref_entries0 > 1) {
                    sh->sh_num_ref_idx_active_l0_minus1 = nvcl_read_u_expgolomb(rdr);
                    nb_ref_entries0 = sh->sh_num_ref_idx_active_l0_minus1 + 1;
                }

                if (sh->sh_slice_type == B && nb_ref_entries1 > 1) {
                    sh->sh_num_ref_idx_active_l1_minus1 = nvcl_read_u_expgolomb(rdr);
                    nb_ref_entries1 = sh->sh_num_ref_idx_active_l1_minus1 + 1;
                }
            } else if (nb_ref_entries0 > pps->pps_num_ref_idx_default_active_minus1[0] ||
                       nb_ref_entries1 > pps->pps_num_ref_idx_default_active_minus1[1]) {

                if (nb_ref_entries0 > pps->pps_num_ref_idx_default_active_minus1[0]) {
                    nb_ref_entries0 = pps->pps_num_ref_idx_default_active_minus1[0] + 1;
                }

                if (sh->sh_slice_type == B &&
                    nb_ref_entries1 > pps->pps_num_ref_idx_default_active_minus1[1]) {
                    nb_ref_entries1 = pps->pps_num_ref_idx_default_active_minus1[1] + 1;
                }
            }
        }

        if (pps->pps_cabac_init_present_flag) {
            sh->sh_cabac_init_flag = nvcl_read_flag(rdr);
        }

        if (ph->ph_temporal_mvp_enabled_flag && !pps->pps_rpl_info_in_ph_flag) {
            sh->sh_collocated_from_l0_flag = 1;
            if (sh->sh_slice_type == B) {
                sh->sh_collocated_from_l0_flag = nvcl_read_flag(rdr);
            }

            if ((sh->sh_collocated_from_l0_flag && nb_ref_entries0 > 1) ||
                (sh->sh_slice_type == B && (!sh->sh_collocated_from_l0_flag && nb_ref_entries1 > 1))) {
                sh->sh_collocated_ref_idx = nvcl_read_u_expgolomb(rdr);
            }
        }

        if (!pps->pps_wp_info_in_ph_flag &&
            ((pps->pps_weighted_pred_flag && sh->sh_slice_type == P) ||
             (pps->pps_weighted_bipred_flag && sh->sh_slice_type == B))) {
            #if 0
            pred_weight_table();
            #endif
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
        if (sps->sps_chroma_format_idc) {
            sh->sh_sao_chroma_used_flag = nvcl_read_flag(rdr);
        }
    }

    if (pps->pps_deblocking_filter_override_enabled_flag && !pps->pps_dbf_info_in_ph_flag) {
        sh->sh_deblocking_params_present_flag = nvcl_read_flag(rdr);

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
    int nb_entry_points_minus1 = (pps->pps_num_tile_columns_minus1 + 1) * (pps->pps_num_tile_rows_minus1 + 1) - 1;
    if (nb_entry_points_minus1) {
        if (nb_entry_points_minus1 > 64) return -1;
        sh->sh_entry_offset_len_minus1 = nvcl_read_u_expgolomb(rdr);
        for (i = 0; i < nb_entry_points_minus1; i++) {
            sh->sh_entry_point_offset_minus1[i] = nvcl_read_bits(rdr, sh->sh_entry_offset_len_minus1 + 1);
        }
    }

    nvcl_read_flag(rdr);
    nvcl_align(rdr);

    return 0;
}
