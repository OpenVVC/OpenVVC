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

/* This file is intended to translate Parameter Sets information into
 * more convenient structures to be used by the decoder contexts
 * It contains translations from nvcl structures to decoder structures
 * and Parameters sets activation functions.
 */
#include "stddef.h"
#include "ovutils.h"
#include "overror.h"
#include "ovunits.h"

#include "decinit.h"
#include "nvcl_structures.h"
#include "dec_structures.h"
#include "rcn_lmcs.h"
#include "hls_structures.h"


static int
sps_init_partition_constraint_info_intra(OVPartInfo *const pinfo, const OVSPS *const sps)
{
    uint8_t log2_ctu_s    = sps->sps_log2_ctu_size_minus5 + 5;
    uint8_t log2_min_cb_s = sps->sps_log2_min_luma_coding_block_size_minus2 + 2;
    uint8_t log2_min_qt_s = sps->sps_log2_diff_min_qt_min_cb_intra_slice_luma + log2_min_cb_s;

    pinfo->log2_ctu_s    = log2_ctu_s;
    pinfo->log2_min_cb_s = log2_min_cb_s;

    pinfo->log2_min_qt_s = log2_min_qt_s;

    pinfo->log2_max_bt_s =  log2_min_qt_s + sps->sps_log2_diff_max_bt_min_qt_intra_slice_luma;
    pinfo->log2_max_tt_s =  log2_min_qt_s + sps->sps_log2_diff_max_tt_min_qt_intra_slice_luma;

    pinfo->max_mtt_depth = sps->sps_max_mtt_hierarchy_depth_intra_slice_luma;

    pinfo->log2_max_tb_s = 5 + sps->sps_max_luma_transform_size_64_flag;

    return 0;
}

static int
sps_init_partition_constraint_info_inter(OVPartInfo *const pinfo, const OVSPS *const sps)
{
    uint8_t log2_ctu_s = sps->sps_log2_ctu_size_minus5 + 5;
    uint8_t log2_min_cb_s = sps->sps_log2_min_luma_coding_block_size_minus2 + 2;
    uint8_t log2_min_qt_s = sps->sps_log2_diff_min_qt_min_cb_inter_slice + log2_min_cb_s;

    pinfo->log2_ctu_s    = log2_ctu_s;
    pinfo->log2_min_cb_s = log2_min_cb_s;

    pinfo->log2_min_qt_s = log2_min_qt_s;

    pinfo->log2_max_bt_s =  log2_min_qt_s + sps->sps_log2_diff_max_bt_min_qt_inter_slice;
    pinfo->log2_max_tt_s =  log2_min_qt_s + sps->sps_log2_diff_max_tt_min_qt_inter_slice;

    pinfo->max_mtt_depth = sps->sps_max_mtt_hierarchy_depth_inter_slice;

    pinfo->log2_max_tb_s = 5 + sps->sps_max_luma_transform_size_64_flag;

    return 0;
}

static int
sps_init_partition_constraint_info_inter_chroma(OVPartInfo *const pinfo, const OVSPS *const sps)
{
    uint8_t log2_ctu_s = sps->sps_log2_ctu_size_minus5 + 5;
    uint8_t log2_min_cb_s = sps->sps_log2_min_luma_coding_block_size_minus2 + 2 - 1;
    uint8_t log2_min_qt_s = sps->sps_log2_diff_min_qt_min_cb_inter_slice + log2_min_cb_s;

    pinfo->log2_ctu_s    = log2_ctu_s;
    pinfo->log2_min_cb_s = log2_min_cb_s;

    pinfo->log2_min_qt_s = log2_min_qt_s;

    pinfo->log2_max_bt_s =  log2_min_qt_s + sps->sps_log2_diff_max_bt_min_qt_inter_slice;
    pinfo->log2_max_tt_s =  log2_min_qt_s + sps->sps_log2_diff_max_tt_min_qt_inter_slice;

    pinfo->max_mtt_depth = sps->sps_max_mtt_hierarchy_depth_inter_slice;

    pinfo->log2_max_tb_s = 5 + sps->sps_max_luma_transform_size_64_flag - 1;

    return 0;
}

static int
ph_override_partition_constraint_info_intra(OVPartInfo *const pinfo, const OVPH *const ph, const OVSPS *const sps)
{
    uint8_t log2_ctu_s    = sps->sps_log2_ctu_size_minus5 + 5;
    uint8_t log2_min_cb_s = sps->sps_log2_min_luma_coding_block_size_minus2 + 2;
    uint8_t log2_min_qt_s = ph->ph_log2_diff_min_qt_min_cb_intra_slice_luma + log2_min_cb_s;

    pinfo->log2_ctu_s    = log2_ctu_s;
    pinfo->log2_min_cb_s = log2_min_cb_s;

    pinfo->log2_min_qt_s = log2_min_qt_s;

    pinfo->log2_max_bt_s =  log2_min_qt_s + ph->ph_log2_diff_max_bt_min_qt_intra_slice_luma;
    pinfo->log2_max_tt_s =  log2_min_qt_s + ph->ph_log2_diff_max_tt_min_qt_intra_slice_luma;

    pinfo->max_mtt_depth = ph->ph_max_mtt_hierarchy_depth_intra_slice_luma;

    pinfo->log2_max_tb_s = 5 + sps->sps_max_luma_transform_size_64_flag;

    return 0;
}

static int
ph_override_partition_constraint_info_inter(OVPartInfo *const pinfo, const OVPH *const ph, const OVSPS *const sps)
{
    uint8_t log2_ctu_s = sps->sps_log2_ctu_size_minus5 + 5;
    uint8_t log2_min_cb_s = sps->sps_log2_min_luma_coding_block_size_minus2 + 2;
    uint8_t log2_min_qt_s = ph->ph_log2_diff_min_qt_min_cb_inter_slice + log2_min_cb_s;

    pinfo->log2_ctu_s    = log2_ctu_s;
    pinfo->log2_min_cb_s = log2_min_cb_s;

    pinfo->log2_min_qt_s = log2_min_qt_s;

    pinfo->log2_max_bt_s =  log2_min_qt_s + ph->ph_log2_diff_max_bt_min_qt_inter_slice;
    pinfo->log2_max_tt_s =  log2_min_qt_s + ph->ph_log2_diff_max_tt_min_qt_inter_slice;

    pinfo->max_mtt_depth = ph->ph_max_mtt_hierarchy_depth_inter_slice;

    pinfo->log2_max_tb_s = 5 + sps->sps_max_luma_transform_size_64_flag;

    return 0;
}

static int
ph_override_partition_constraint_info_inter_chroma(OVPartInfo *const pinfo, const OVPH *const ph, const OVSPS *const sps)
{
    uint8_t log2_ctu_s = sps->sps_log2_ctu_size_minus5 + 5;
    uint8_t log2_min_cb_s = sps->sps_log2_min_luma_coding_block_size_minus2 + 2 - 1;
    uint8_t log2_min_qt_s = ph->ph_log2_diff_min_qt_min_cb_inter_slice + log2_min_cb_s;

    pinfo->log2_ctu_s    = log2_ctu_s;
    pinfo->log2_min_cb_s = log2_min_cb_s;

    pinfo->log2_min_qt_s = log2_min_qt_s;

    pinfo->log2_max_bt_s =  log2_min_qt_s + ph->ph_log2_diff_max_bt_min_qt_inter_slice;
    pinfo->log2_max_tt_s =  log2_min_qt_s + ph->ph_log2_diff_max_tt_min_qt_inter_slice;

    pinfo->max_mtt_depth = ph->ph_max_mtt_hierarchy_depth_inter_slice;

    pinfo->log2_max_tb_s = 5 + sps->sps_max_luma_transform_size_64_flag - 1;

    return 0;
}

static int
ph_override_partition_constraint_info_chroma(OVPartInfo *const pinfo, const OVPH *const ph, const OVSPS *const sps)
{

    pinfo->log2_ctu_s    = sps->sps_log2_ctu_size_minus5 + 5;
    pinfo->log2_min_cb_s = sps->sps_log2_min_luma_coding_block_size_minus2 + 2 - 1;


    pinfo->log2_min_qt_s = ph->ph_log2_diff_min_qt_min_cb_intra_slice_chroma +
            pinfo->log2_min_cb_s;

    pinfo->log2_max_bt_s = pinfo->log2_min_qt_s +
            ph->ph_log2_diff_max_bt_min_qt_intra_slice_chroma;
    pinfo->log2_max_tt_s = pinfo->log2_min_qt_s +
            ph->ph_log2_diff_max_tt_min_qt_intra_slice_chroma;

    pinfo->max_mtt_depth = ph->ph_max_mtt_hierarchy_depth_intra_slice_chroma;

    pinfo->log2_max_tb_s = 5 + sps->sps_max_luma_transform_size_64_flag - 1;

    return 0;
}

static int
sps_init_partition_constraint_info_chroma(OVPartInfo *const pinfo, const OVSPS *const sps)
{
    /* FIXME we could handle chroma format from here */
    #if 0
    uint8_t log2_ctu_s = sps->sps_log2_ctu_size_minus5 + 5;
    uint8_t log2_min_cb_s = sps->sps_log2_min_luma_coding_block_size_minus2 + 2 - 1;
    /* FIXME check this  factorize if correct*/
    uint8_t log2_min_qt_s = sps->sps_log2_diff_min_qt_min_cb_intra_slice_chroma + log2_min_cb_s;
    #endif

    pinfo->log2_ctu_s    = sps->sps_log2_ctu_size_minus5 + 5;
    pinfo->log2_min_cb_s = sps->sps_log2_min_luma_coding_block_size_minus2 + 2 - 1;


    pinfo->log2_min_qt_s = sps->sps_log2_diff_min_qt_min_cb_intra_slice_chroma +
            pinfo->log2_min_cb_s;

    pinfo->log2_max_bt_s = pinfo->log2_min_qt_s +
            sps->sps_log2_diff_max_bt_min_qt_intra_slice_chroma;
    pinfo->log2_max_tt_s = pinfo->log2_min_qt_s +
            sps->sps_log2_diff_max_tt_min_qt_intra_slice_chroma;

    pinfo->max_mtt_depth = sps->sps_max_mtt_hierarchy_depth_intra_slice_chroma;

    pinfo->log2_max_tb_s = 5 + sps->sps_max_luma_transform_size_64_flag - 1;

    return 0;
}


static void
sps_fill_qp_table(uint8_t dst_qp_tab[],
                  const uint8_t dqp_input_val_min1[],
                  const uint8_t dqp_diff_val[],
                  int8_t start_qp, uint8_t nb_qp_min1)
{
    uint8_t nb_qp_points = nb_qp_min1 + 1;
    int idx = start_qp;
    int next_idx;
    int j;

    /*FIXME:
     *    -actual qp can be clipped to -qp_bitdepth_offset
     *    -fill first part of qp_map + check corner cases
     */

    dst_qp_tab[idx] = idx;
    for (j = idx - 1; j >= 0; --j){
        dst_qp_tab[j] = dst_qp_tab[j + 1] - 1;
    }

    for (j = 0; j <= nb_qp_points; ++j) {
        const int nb_step = dqp_input_val_min1[j] + 1;
        const int nb_qp_w = dqp_input_val_min1[j] ^ dqp_diff_val[j];
        const int first_qp = dst_qp_tab[idx];
        const int round = nb_step >> 1;
        int qp_sum = nb_qp_w;
        int k;

        next_idx = idx + nb_step;

        for (k = idx + 1;  k <= next_idx; ++k) {
            dst_qp_tab[k]  = first_qp + (qp_sum + round) / nb_step;
            qp_sum += nb_qp_w;
        }
        idx = next_idx;
    }

    for (j = next_idx;  j < 64; ++j) {
        dst_qp_tab[j] = ov_clip(dst_qp_tab[j - 1] + 1, 0, 63);
    }
}

static void
sps_init_chroma_qp_tables(struct SPSInfo *sps_info, const OVSPS *const sps)
{
    int i = 0;
    const int nb_qp_tab = sps->sps_same_qp_table_for_chroma_flag ? 1 : 2 + sps->sps_joint_cbcr_enabled_flag;
    for (i = 0; i < nb_qp_tab; ++i) {
        const uint8_t nb_qp_min1 = sps->sps_num_points_in_qp_table_minus1[i];
        const uint8_t *dqp_input_val_min1 = sps->sps_delta_qp_in_val_minus1[i];
        const uint8_t *dqp_diff_val = sps->sps_delta_qp_diff_val[i];
        int8_t start_qp = sps->sps_qp_table_start_minus26[i] + 26;
        uint8_t *dst_qp_tab = (uint8_t *)sps_info->qp_tables_c[i].qp;

        sps_fill_qp_table(dst_qp_tab, dqp_input_val_min1, dqp_diff_val, start_qp,
                          nb_qp_min1);
    }

    /* Copy first table into others in case of same qp for all chroma tables */
    for (; i < 3 - !sps->sps_joint_cbcr_enabled_flag; ++i) {
        sps_info->qp_tables_c[i] = sps_info->qp_tables_c[i-1];
    }
}

static int
update_sps_info(struct SPSInfo *const sps_info, const OVSPS *const sps)
{
    sps_init_partition_constraint_info_intra(&sps_info->part_info[0], sps);
    sps_init_partition_constraint_info_inter(&sps_info->part_info[1], sps);
    sps_init_partition_constraint_info_chroma(&sps_info->part_info[2], sps);
    sps_init_partition_constraint_info_inter_chroma(&sps_info->part_info[3], sps);

    sps_init_chroma_qp_tables(sps_info, sps);

    if (sps->sps_vui_parameters_present_flag) {
        sps_info->color_desc.colour_primaries         = sps->vui.vui_colour_primaries;
        sps_info->color_desc.transfer_characteristics = sps->vui.vui_transfer_characteristics;
        sps_info->color_desc.matrix_coeffs            = sps->vui.vui_matrix_coeffs;
        sps_info->color_desc.full_range               = sps->vui.vui_full_range_flag;
    }

    return 0;
}

int
decinit_set_entry_points(OVPS *const prms, const OVNALUnit *nal, uint32_t nb_sh_bytes)
{
    int i, j;
    struct SHInfo *const sh_info = &prms->sh_info;
    const OVSH *const sh = prms->sh;
    struct TileInfo *const tinfo = &prms->pps_info.tile_info;
    /* FIXME we consider nb_entries is nb_tiles */
    /* TODO compute and keep track of nb_tiles from pps */
    int nb_entries = tinfo->nb_tile_cols * tinfo->nb_tile_rows;
    uint32_t rbsp_offset[256];
    const int nb_rbsp_epb = nal->nb_epb;
    const uint32_t *rbsp_epb_pos = nal->epb_pos;
    const uint8_t *const rbsp_end = nal->rbsp_data + nal->rbsp_size;
    int nb_sh_epb = 0;

    rbsp_offset[0] = 0;

    for (j = 0; j < nb_rbsp_epb; ++j) {
        nb_sh_epb += rbsp_epb_pos[j] <= nb_sh_bytes;
    }

    for (i = 0; i < nb_entries; ++i) {
        uint32_t entry_offset = sh->sh_entry_point_offset_minus1[i] + 1;
        rbsp_offset[i + 1] = rbsp_offset[i] + entry_offset;
    }

    for (i = 0; i < nb_entries; ++i) {
        for (j = nb_sh_epb; j < nb_rbsp_epb; ++j) {
            uint32_t entry_offset = rbsp_offset[i + 1];
            entry_offset -= (entry_offset > (rbsp_epb_pos[j] - nb_sh_bytes));
            rbsp_offset[i + 1] = entry_offset;
        }
    }

    /* FIXME avoid using an offset tab and compute directly on entries */
    sh_info->rbsp_entry[0] = nal->rbsp_data + nb_sh_bytes;
    for (i = 1; i < nb_entries; ++i) {
        sh_info->rbsp_entry[i] = OVMIN(nal->rbsp_data + rbsp_offset[i] + nb_sh_bytes, rbsp_end);
    }

    /* Note this is so we can retrieve entry end by using rbsp_entry [i + 1] */
    sh_info->rbsp_entry[i] = rbsp_end;

    /*FIXME check entries do not exceed rpbs size */
    return 0;
}

static void
init_tile_ctx(struct TileInfo *const tinfo, const OVPPS *const pps)
{
    int nb_cols = pps->pps_num_tile_columns_minus1 + 1;
    int nb_rows = pps->pps_num_tile_rows_minus1 + 1;

    const int log2_ctu_s = pps->pps_log2_ctu_size_minus5 + 5;

    const int pic_w = pps->pps_pic_width_in_luma_samples;
    const int pic_h = pps->pps_pic_height_in_luma_samples;

    /* FIXME harmonize this with sps
    */
    const int nb_ctu_w = (pic_w >> log2_ctu_s) + !!(pic_w & ((1 << log2_ctu_s) - 1));
    const int nb_ctu_h = (pic_h >> log2_ctu_s) + !!(pic_h & ((1 << log2_ctu_s) - 1));

    int rem_ctu_w = nb_ctu_w;
    int rem_ctu_h = nb_ctu_h;

    int i;
    for (i = 1; i <=  nb_rows; ++i) {
        const int tile_nb_ctu_h = pps->pps_tile_row_height_minus1[i - 1] + 1;
        tinfo->ctu_y[i - 1] = nb_ctu_h - rem_ctu_h;
        tinfo->nb_ctu_h[i - 1] = tile_nb_ctu_h;
        rem_ctu_h -= tile_nb_ctu_h;
    }

    for (i = 1; i <= nb_cols; ++i) {
        const int tile_nb_ctu_w = pps->pps_tile_column_width_minus1[i - 1] + 1;
        tinfo->ctu_x[i - 1] = nb_ctu_w - rem_ctu_w;
        tinfo->nb_ctu_w[i - 1] = tile_nb_ctu_w;
        rem_ctu_w -= tile_nb_ctu_w;
    }

    tinfo->nb_tile_cols = nb_cols ;
    tinfo->nb_tile_rows = nb_rows ;
}

static int
check_pps_dimension(const OVPPS *const pps, const OVSPS *const sps)
{
     uint8_t pps_w = pps->pps_pic_width_in_luma_samples > sps->sps_pic_width_max_in_luma_samples;
     uint8_t pps_h = pps->pps_pic_height_in_luma_samples > sps->sps_pic_height_max_in_luma_samples;

     return -(pps_w | pps_h);
}

static int
update_pps_info(struct PPSInfo *const pps_info, const OVPPS *const pps,
                const OVSPS *const sps)
{
    struct TileInfo *const tinfo = &pps_info->tile_info;
    uint8_t have_tile = pps->pps_num_tile_columns_minus1 + pps->pps_num_tile_rows_minus1;

    if (have_tile) {
        init_tile_ctx(tinfo, pps);
    } else {
        /* Initialize tile context to 1 tile on whole picture based on
         * SPS informations
         */
        const int log2_ctu_s = sps->sps_log2_ctu_size_minus5 + 5;

        const int pic_w = pps->pps_pic_width_in_luma_samples;
        const int pic_h = pps->pps_pic_height_in_luma_samples;

        const int nb_ctu_w = (pic_w + ((1 << log2_ctu_s) - 1)) >> log2_ctu_s;
        const int nb_ctu_h = (pic_h + ((1 << log2_ctu_s) - 1)) >> log2_ctu_s;

        tinfo->ctu_x[0] = 0;
        tinfo->ctu_y[0] = 0;
        /* FIXME log2_ctu_s from pps is invalid not present
         * Use SPS instead
         */
        tinfo->nb_ctu_w[0] = nb_ctu_w;
        tinfo->nb_ctu_h[0] = nb_ctu_h;

        tinfo->nb_tile_cols = 1;
        tinfo->nb_tile_rows = 1;
    }

    if (check_pps_dimension(pps, sps)) {
        ov_log(NULL, OVLOG_ERROR,
               "PPS picture dimension %dx%d bigger than max picture dimension in SPS %dx%d.\n",
               pps->pps_pic_width_in_luma_samples, pps->pps_pic_height_in_luma_samples,
               sps->sps_pic_width_max_in_luma_samples, sps->sps_pic_height_max_in_luma_samples);

        return OVVC_EINDATA;
    }

    return 0;
}

static int
update_ph_info(struct SPSInfo *const sps_info, const OVPH *const ph, const OVSPS *const sps)
{
    if (ph->ph_partition_constraints_override_flag) {
        /*FIXME test for dual tree etc. */
        ph_override_partition_constraint_info_intra(&sps_info->part_info[0], ph, sps);
        ph_override_partition_constraint_info_inter(&sps_info->part_info[1], ph, sps);
        ph_override_partition_constraint_info_chroma(&sps_info->part_info[2], ph, sps);
        ph_override_partition_constraint_info_inter_chroma(&sps_info->part_info[3], ph, sps);
    } else {
        sps_init_partition_constraint_info_intra(&sps_info->part_info[0], sps);
        sps_init_partition_constraint_info_inter(&sps_info->part_info[1], sps);
        sps_init_partition_constraint_info_chroma(&sps_info->part_info[2], sps);
        sps_init_partition_constraint_info_inter_chroma(&sps_info->part_info[3], sps);
    }

    return 0;
}

static int
update_sh_info(struct SHInfo *const sh_info, const OVSH *const sh)
{
    return 0;
}


static OVAPS *
retrieve_aps_alf(const OVNVCLCtx *const nvcl_ctx, uint8_t aps_id)
{
    OVAPS *aps = NULL;
    if (aps_id < 16) {
        aps = nvcl_ctx->alf_aps_list[aps_id];
    } else {
        ov_log(NULL, 3, "Invalid APS ID  %d\n", aps_id);
    }
    return aps;
}


static OVAPS *
retrieve_aps_lmcs(const OVNVCLCtx *const nvcl_ctx, const OVPH *const ph)
{
    uint8_t aps_id = ph->ph_lmcs_aps_id;
    OVAPS *aps = NULL;
    if (aps_id < 16) {
        aps = nvcl_ctx->lmcs_aps_list[aps_id];
    } else {
        ov_log(NULL, 3, "Invalid APS ID  %d\n", aps_id);
    }
    return aps;
}

static OVAPS *
retrieve_aps_scaling_list(const OVNVCLCtx *const nvcl_ctx, const OVPH *const ph)
{
    uint8_t aps_id = ph->ph_scaling_list_aps_id;
    OVAPS *aps = NULL;
    if (aps_id < 16) {
        aps = nvcl_ctx->scaling_list_aps_list[aps_id];
    } else {
        ov_log(NULL, 3, "Invalid APS ID  %d\n", aps_id);
    }
    return aps;
}

static OVSPS *
retrieve_sps(const OVNVCLCtx *const nvcl_ctx, const OVPPS *const pps)
{
    uint8_t sps_id = pps->pps_seq_parameter_set_id;
    OVSPS *sps = NULL;
    if (sps_id < 16) {
        sps = (OVSPS *)nvcl_ctx->sps_list[sps_id]->data;
    } else {
        ov_log(NULL, 3, "Invalid SPS ID  %d\n", sps_id);
    }

    return sps;
}

static OVPPS *
retrieve_pps(const OVNVCLCtx *const nvcl_ctx, const OVPH *const ph)
{
    uint8_t pps_id = ph->ph_pic_parameter_set_id;
    OVPPS *pps = NULL;
    if (pps_id < 16) {
        pps = (OVPPS *)nvcl_ctx->pps_list[pps_id]->data;
    } else {
        ov_log(NULL, 3, "Invalid PPS ID  %d\n", pps_id);
    }

    return pps;
}

static void
set_pic_part_info(struct PicPartInfo *pic_info, const OVSPS *const sps, const OVPPS *const pps)
{
     /* Masks are to ensure log2_size does not exceed standard requirements */
     uint8_t log2_ctb_s    = (sps->sps_log2_ctu_size_minus5 + 5) & 0x7;
     uint8_t log2_min_cb_s = (sps->sps_log2_min_luma_coding_block_size_minus2 + 2) & 0x7;
     /* FIXME assert log2_min < log2_ctb */

     uint16_t pic_w = pps->pps_pic_width_in_luma_samples;
     uint16_t pic_h = pps->pps_pic_height_in_luma_samples;

     uint16_t nb_ctb_pic_w = (pic_w + ((1 << log2_ctb_s) - 1)) >> log2_ctb_s;
     uint16_t nb_ctb_pic_h = (pic_h + ((1 << log2_ctb_s) - 1)) >> log2_ctb_s;

     uint16_t nb_pb_pic_w = (nb_ctb_pic_w << log2_ctb_s) >> log2_min_cb_s;
     uint16_t nb_pb_pic_h = (nb_ctb_pic_h << log2_ctb_s) >> log2_min_cb_s;

     pic_info->log2_ctu_s = log2_ctb_s;
     pic_info->log2_min_cb_s = log2_min_cb_s;

     pic_info->pic_w = pic_w;
     pic_info->pic_h = pic_h;

     pic_info->nb_ctb_w = nb_ctb_pic_w;
     pic_info->nb_ctb_h = nb_ctb_pic_h;

     pic_info->nb_pb_w = nb_pb_pic_w;
     pic_info->nb_pb_h = nb_pb_pic_h;
}

static void
set_max_pic_part_info(struct PicPartInfo *pic_info, const OVSPS *const sps, const OVPPS *const pps)
{
     /* Masks are to ensure log2_size does not exceed standard requirements */
     uint8_t log2_ctb_s    = (sps->sps_log2_ctu_size_minus5 + 5) & 0x7;
     uint8_t log2_min_cb_s = (sps->sps_log2_min_luma_coding_block_size_minus2 + 2) & 0x7;
     /* FIXME assert log2_min < log2_ctb */

     uint16_t pic_w = sps->sps_pic_width_max_in_luma_samples;
     uint16_t pic_h = sps->sps_pic_height_max_in_luma_samples;

     uint16_t nb_ctb_pic_w = (pic_w + ((1 << log2_ctb_s) - 1)) >> log2_ctb_s;
     uint16_t nb_ctb_pic_h = (pic_h + ((1 << log2_ctb_s) - 1)) >> log2_ctb_s;

     uint16_t nb_pb_pic_w = (nb_ctb_pic_w << log2_ctb_s) >> log2_min_cb_s;
     uint16_t nb_pb_pic_h = (nb_ctb_pic_h << log2_ctb_s) >> log2_min_cb_s;

     pic_info->log2_ctu_s = log2_ctb_s;
     pic_info->log2_min_cb_s = log2_min_cb_s;

     pic_info->pic_w = pic_w;
     pic_info->pic_h = pic_h;

     pic_info->nb_ctb_w = nb_ctb_pic_w;
     pic_info->nb_ctb_h = nb_ctb_pic_h;

     pic_info->nb_pb_w = nb_pb_pic_w;
     pic_info->nb_pb_h = nb_pb_pic_h;

}

void
decinit_unref_params(struct OVPS *const ps)
{
    hlsdata_unref(&ps->sps_ref);
    hlsdata_unref(&ps->pps_ref);
    hlsdata_unref(&ps->ph_ref);
    hlsdata_unref(&ps->sh_ref);
    ps->sps = NULL;
    ps->pps = NULL;
    ps->ph = NULL;
    ps->sh = NULL;
}

uint8_t
sps_check_dimension_change(const OVSPS *const old, const OVSPS *const new)
{
    if (old) {
        uint8_t diff_w = old->sps_pic_width_max_in_luma_samples != new->sps_pic_width_max_in_luma_samples;
        uint8_t diff_h = old->sps_pic_height_max_in_luma_samples != new->sps_pic_height_max_in_luma_samples;
        return diff_w | diff_h;
    }
    return 1;
}

int
decinit_update_params(struct OVPS *const ps, const OVNVCLCtx *const nvcl_ctx)
{
    /* FIXME assert nvcl_ctx params sets are not NULL*/
    int ret;
    OVSH * sh = (OVSH *)nvcl_ctx->sh->data;
    OVPH * ph = (OVPH *)nvcl_ctx->ph->data;
    OVPPS * pps = retrieve_pps(nvcl_ctx, ph);
    OVSPS * sps = retrieve_sps(nvcl_ctx, pps);

    ps->sei = nvcl_ctx->sei;

    if (!sh || !ph || !pps || !sps) {
        ov_log(NULL, 3, "Missing Parameter sets for dec initialisation\n");
        return OVVC_EINDATA;
    }

    uint8_t pps_id = ph->ph_pic_parameter_set_id;
    uint8_t sps_id = pps->pps_seq_parameter_set_id;

    if (ps->sps != sps) {
        ret = update_sps_info(&ps->sps_info, sps);
        if (ret < 0) {
            goto failsps;
        }
        if (ps->sps) {
            uint8_t req_dpb_realloc = sps_check_dimension_change(ps->sps, sps);
            ps->sps_info.req_dpb_realloc = req_dpb_realloc;
        }
        hlsdata_unref(&ps->sps_ref);
        hlsdata_newref(&ps->sps_ref, nvcl_ctx->sps_list[sps_id]);
        ps->sps = sps;
    }

    if (ps->pps != pps) {
        /* FIXME pps could trigger use a new sps */
        ret = update_pps_info(&ps->pps_info, pps, ps->sps);
        if (ret < 0) {
            goto failpps;
        }
        hlsdata_unref(&ps->pps_ref);
        hlsdata_newref(&ps->pps_ref, nvcl_ctx->pps_list[pps_id]);
        ps->pps = pps;
    }

    if (ps->ph != ph) {
        /* FIXME use ph_info instead of sps_info */
        ret = update_ph_info(&ps->sps_info, ph, sps);
        if (ret < 0) {
            goto failph;
        }
        hlsdata_unref(&ps->ph_ref);
        hlsdata_newref(&ps->ph_ref,  nvcl_ctx->ph);
        ps->ph = ph;
    }

    /* We only check sh is present in ctx since we should be called
     * directely after the Slice Header is read so we are
     * reading a new slice
     */
    if (ps->sh != sh) {
        ret = update_sh_info(&ps->sh_info, (OVSH *)nvcl_ctx->sh->data);
        if (ret < 0) {
            goto failsh;
        }
        hlsdata_unref(&ps->sh_ref);
        hlsdata_newref(&ps->sh_ref,  nvcl_ctx->sh);
        ps->sh = (OVSH *)nvcl_ctx->sh->data;
    }

    for(int i = 0; i < sh->sh_num_alf_aps_ids_luma; i++){
        uint8_t aps_id = sh->sh_alf_aps_id_luma[i];
        OVAPS * aps_alf = retrieve_aps_alf(nvcl_ctx, aps_id);
        if (ps->aps_alf[i] != aps_alf) {
            ps->aps_alf[i] = aps_alf;
        }
    }

    for (int i=sh->sh_num_alf_aps_ids_luma; i < 8; i++){
        ps->aps_alf[i] = NULL;
    }

    OVAPS * aps_alf_c = retrieve_aps_alf(nvcl_ctx, sh->sh_alf_aps_id_chroma);
    if (ps->aps_alf_c != aps_alf_c) {
        ps->aps_alf_c = aps_alf_c;
    }

    OVAPS * aps_cc_alf_cb = retrieve_aps_alf(nvcl_ctx, sh->sh_alf_cc_cb_aps_id);
    if (ps->aps_cc_alf_cb != aps_cc_alf_cb) {
        ps->aps_cc_alf_cb = aps_cc_alf_cb;
    }
    OVAPS * aps_cc_alf_cr = retrieve_aps_alf(nvcl_ctx, sh->sh_alf_cc_cr_aps_id);
    if (ps->aps_cc_alf_cr != aps_cc_alf_cr) {
        ps->aps_cc_alf_cr = aps_cc_alf_cr;
    }

    OVAPS * aps_lmcs = retrieve_aps_lmcs(nvcl_ctx, ph);
    if (ps->aps_lmcs != aps_lmcs) {
        ps->aps_lmcs = aps_lmcs;
    }

    OVAPS * aps_scaling_list = retrieve_aps_scaling_list(nvcl_ctx, ph);
    if (ps->aps_scaling_list != aps_scaling_list) {
        ps->aps_scaling_list = aps_scaling_list;
    }

    set_pic_part_info(&ps->pic_info, ps->sps, ps->pps);
    set_max_pic_part_info(&ps->pic_info_max, ps->sps, ps->pps);

    return 0;

/* TODO if some alloc are done from update function free it here*/
failsps:
failpps:
failph:
failsh:
    /* On failure  we reset all active ps so next activation will try again */
    ps->sps = NULL;
    ps->pps = NULL;
    ps->ph = NULL;
    ps->sh = NULL;

    return ret;
}
