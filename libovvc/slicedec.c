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
#include "ovdefs.h"
#include "overror.h"
#include "ovmem.h"
#include "ovutils.h"
#include "slicedec.h"
#include "ctudec.h"
#include "nvcl_structures.h"
#include "vcl_cabac.h"
#include "vcl.h"
#include "drv_utils.h"
#include "rcn.h"
#include "ovdpb.h"
#include "drv_lines.h"
#include "ovthreads.h"
#include "rcn_lmcs.h"
#include "rcn_dequant.h"
#include "hls_structures.h"


/* TODO define in a header */
enum SliceType {
     SLICE_B = 0,
     SLICE_P = 1,
     SLICE_I = 2
};

static int
slicedec_decode_rect_entry(OVSliceDec *sldec, OVCTUDec *const ctudec, const OVPS *const prms,
                           uint16_t entry_idx);

static void derive_ctu_neighborhood(OVCTUDec *const ctudec,
                                    int ctb_address, int nb_ctu_w);
static void
init_slice_tree_ctx(OVCTUDec *const ctudec, const struct OVPS *prms)
{
    /* TODO use a specific structure for handling trees */
    const OVSPS *sps = prms->sps;
    const OVSH *sh = prms->sh;

    if (sh->sh_slice_type == SLICE_I) {
        if (sps->sps_qtbtt_dual_tree_intra_flag) {
            ctudec->coding_tree          = &dual_tree;
            ctudec->coding_tree_implicit = &dual_tree_implicit;
            ctudec->coding_unit          = coding_unit_intra;
            ctudec->transform_unit       = transform_unit_l;
        } else {
            ctudec->coding_tree          = coding_quadtree;
            ctudec->coding_tree_implicit = coding_quadtree_implicit;
            ctudec->coding_unit          = coding_unit_intra_st;
            ctudec->transform_unit       = transform_unit_st;
        }
    } else {
        ctudec->coding_tree          = coding_quadtree;
        ctudec->coding_tree_implicit = coding_quadtree_implicit;
        ctudec->coding_unit          = coding_unit_inter_st;
        ctudec->prediction_unit      = sh->sh_slice_type == SLICE_B ? prediction_unit_inter_b : &prediction_unit_inter_p;
        ctudec->transform_unit = transform_unit_st;
    }
    /* FIXME  move default active part map selection to somewhere else*/
    ctudec->active_part_map = &ctudec->part_map;
}


#if 0
static int
init_filter_info()
{
    uint8_t flag = 0;
    if (sh->slice_sao_luma_flag)
        flag |= VVC_SAO_LUMA_SLICE_FLAG;
    if (sh->slice_sao_chroma_flag)
        flag |= VVC_SAO_CHROMA_SLICE_FLAG;
    if (sh->slice_alf_enabled_flag) {
        flag |= VVC_ALF_LUMA_SLICE_FLAG;
        flag |= (sh->slice_alf_cb_enabled_flag & 1) ? VVC_ALF_CB_SLICE_FLAG : 0;
        flag |= (sh->slice_alf_cr_enabled_flag & 2) ? VVC_ALF_CR_SLICE_FLAG : 0;
    }
    ctudec->alf_num_chroma_alt = vvc_ctx->alf_num_alt_chroma;

    ctudec->loop_filter_state_flags = flag;
    /* ALF */
    uint8_t sh_alf_enabled_flag;
    uint8_t sh_alf_cb_enabled_flag;
    uint8_t sh_alf_cr_enabled_flag;
    uint8_t sh_alf_cc_cb_enabled_flag;
    uint8_t sh_alf_cc_cr_enabled_flag;

    /* SAO from PH */
    uint8_t ph_sao_luma_enabled_flag;
    uint8_t ph_sao_chroma_enabled_flag;

    /* SAO */
    uint8_t sh_sao_luma_used_flag;
    uint8_t sh_sao_chroma_used_flag;

    /* ALF info */
    uint8_t sh_num_alf_aps_ids_luma;
    uint8_t sh_alf_aps_id_luma[16];
    uint8_t sh_alf_aps_id_chroma;
    uint8_t sh_alf_cc_cb_aps_id;
    uint8_t sh_alf_cc_cr_aps_id;

    /* LMCS */
    uint8_t sh_lmcs_used_flag;

    /* DBF */
    uint8_t sh_deblocking_params_present_flag;
    uint8_t sh_deblocking_filter_disabled_flag;
    int8_t sh_luma_beta_offset_div2;
    int8_t sh_luma_tc_offset_div2;
    int8_t sh_cb_beta_offset_div2;
    int8_t sh_cb_tc_offset_div2;
    int8_t sh_cr_beta_offset_div2;
    int8_t sh_cr_tc_offset_div2;

    /* DBF from PH */
    uint8_t ph_deblocking_params_present_flag;
    uint8_t ph_deblocking_filter_disabled_flag;
    int8_t ph_luma_beta_offset_div2;
    int8_t ph_luma_tc_offset_div2;
    int8_t ph_cb_beta_offset_div2;
    int8_t ph_cb_tc_offset_div2;
    int8_t ph_cr_beta_offset_div2;
    int8_t ph_cr_tc_offset_div2;

}

static int
init_qp_ctx()
{
    uint8_t ph_qp_delta;
    uint8_t ph_joint_cbcr_sign_flag;
    int8_t sh_qp_delta;
    int8_t sh_cb_qp_offset;
    int8_t sh_cr_qp_offset;
    int8_t sh_joint_cbcr_qp_offset;
    uint8_t sh_cu_chroma_qp_offset_enabled_flag;
    uint8_t ph_qp_delta;
}

static int
dequant_tools()
{
}
#endif

static void
init_coding_coeff_coding_ctx(OVCTUDec *ctudec, const OVPS *prms)
{
    const OVSPS *const sps = prms->sps;
    const OVSH *const sh = prms->sh;

    #if 0
    uint8_t ph_joint_cbcr_sign_flag;

    uint8_t sh_ts_residual_coding_disabled_flag;

    uint8_t ph_chroma_residual_scale_flag;
    uint8_t ph_explicit_scaling_list_enabled_flag;
    uint8_t ph_scaling_list_aps_id;
    #endif

    /* FIXME replace this with a status on MTS */
    uint8_t mts_enabled = sps->sps_mts_enabled_flag;
    ctudec->mts_enabled = mts_enabled;
    ctudec->mts_explicit_intra = mts_enabled && sps->sps_explicit_mts_intra_enabled_flag;
    ctudec->mts_explicit_inter = mts_enabled && sps->sps_explicit_mts_inter_enabled_flag;
    ctudec->mts_implicit = sps->sps_mts_enabled_flag && !sps->sps_explicit_mts_intra_enabled_flag;

    if (sh->sh_dep_quant_used_flag) {
        ctudec->residual_coding_l = residual_coding_dpq;
        ctudec->residual_coding_c = residual_coding_chroma_dpq;
        ctudec->residual_coding_isp_h = residual_coding_isp_h_dpq;
        ctudec->residual_coding_isp_v = residual_coding_isp_v_dpq;
    } else {
        ctudec->residual_coding_isp_h = residual_coding_isp_h_sdh;
        ctudec->residual_coding_isp_v = residual_coding_isp_v_sdh;
        ctudec->residual_coding_c = residual_coding_chroma_sdh;
        ctudec->residual_coding_l = residual_coding_sdh;
        ctudec->enable_sdh = sh->sh_sign_data_hiding_used_flag;
    }
}

#if 0
static int
init_inter_ctx()
{
    struct OVHRPL hrpl;
    uint8_t ph_temporal_mvp_enabled_flag;
    uint8_t ph_collocated_ref_idx;
    uint8_t ph_mmvd_fullpel_only_flag;
    uint8_t ph_mvd_l1_zero_flag;
    uint8_t ph_bdof_disabled_flag;
    uint8_t ph_dmvr_disabled_flag;
    uint8_t ph_prof_disabled_flag;
}

static int
init_partition_info()
{
    /* Separate Picture size info from CTU sie info */
    uint8_t ph_virtual_boundaries_present_flag;
    uint8_t ph_num_ver_virtual_boundaries;
    uint8_t ph_virtual_boundary_pos_x_minus1[32];
    uint8_t ph_num_hor_virtual_boundaries;
    uint8_t ph_virtual_boundary_pos_y_minus1[32];

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
}

#endif

int
slicedec_alloc_cabac_lines(OVSliceDec *const sldec, const struct OVPS *const prms)
{
   #if 0
   int nb_pu_w;


   ov_mallocz(sizeof(
   #endif
   return 0;

}

static void
slice_init_qp_ctx(OVCTUDec *const ctudec, const struct OVPS *const prms)
{
    const OVSPS *const sps = prms->sps;
    const OVPPS *const pps = prms->pps;
    const OVSH *const sh = prms->sh;
    const OVPH *const ph = prms->ph;
    VVCQPCTX *const qp_ctx = &ctudec->qp_ctx;

    /*FIXME check if not done in dec init */
    const int8_t qp_bd_offset = 6 * sps->sps_bitdepth_minus8;
    int8_t pic_base_qp = pps->pps_init_qp_minus26 + 26;
    int8_t pic_qp = pic_base_qp + ph->ph_qp_delta;
    int8_t slice_qp = pic_qp + sh->sh_qp_delta;
    int8_t cb_qp_offset = sh->sh_cb_qp_offset + pps->pps_cb_qp_offset;
    int8_t cr_qp_offset = sh->sh_cr_qp_offset + pps->pps_cr_qp_offset;
    int8_t jcbcr_qp_offset = sh->sh_joint_cbcr_qp_offset + pps->pps_joint_cbcr_qp_offset_value;
    uint8_t cu_qp_delta_subdiv = sh->sh_slice_type == 2 ? ph->ph_cu_qp_delta_subdiv_intra_slice : ph->ph_cu_qp_delta_subdiv_inter_slice;
    uint8_t cu_qp_chroma_offset_subdiv = sh->sh_slice_type == 2 ? ph->ph_cu_chroma_qp_offset_subdiv_intra_slice : ph->ph_cu_chroma_qp_offset_subdiv_inter_slice;
    ctudec->chroma_qp_offset_enabled = sh->sh_cu_chroma_qp_offset_enabled_flag;
    ctudec->chroma_qp_offset_len = pps->pps_chroma_qp_offset_list_len_minus1 + 1;

    ctudec->slice_qp = pic_qp + sh->sh_qp_delta;

    ctudec->cu_qp_delta_subdiv = cu_qp_delta_subdiv;
    ctudec->cu_qp_chroma_offset_subdiv = cu_qp_chroma_offset_subdiv;
    /* FIXME
     * check tables are valid when same qp table for all
     */
    qp_ctx->chroma_qp_map_cb    = prms->sps_info.qp_tables_c[0].qp;
    qp_ctx->chroma_qp_map_cr    = prms->sps_info.qp_tables_c[1].qp;
    qp_ctx->chroma_qp_map_jcbcr = prms->sps_info.qp_tables_c[2].qp;

    memcpy(qp_ctx->cb_qp_offset_list, pps->pps_cb_qp_offset_list, sizeof(int8_t)*16);
    memcpy(qp_ctx->cr_qp_offset_list, pps->pps_cr_qp_offset_list, sizeof(int8_t)*16);
    memcpy(qp_ctx->joint_cbcr_qp_offset_list, pps->pps_joint_cbcr_qp_offset_list, sizeof(int8_t)*16);
    qp_ctx->current_qp = slice_qp;
    qp_ctx->cb_offset = cb_qp_offset;
    qp_ctx->cr_offset = cr_qp_offset;
    qp_ctx->jcbcr_offset = jcbcr_qp_offset;
    qp_ctx->qp_bd_offset = qp_bd_offset;
    qp_ctx->min_qp_prime_ts = 4 + 6 * sps->sps_min_qp_prime_ts;

    ctudec->dequant_luma.qp = slice_qp + qp_bd_offset;
    ctudec->dequant_cb.qp = qp_ctx->chroma_qp_map_cb[ov_clip(slice_qp, 0, 63)] + cb_qp_offset + qp_bd_offset;
    ctudec->dequant_cr.qp = qp_ctx->chroma_qp_map_cr[ov_clip(slice_qp, 0, 63)] + cr_qp_offset + qp_bd_offset;
    ctudec->dequant_joint_cb_cr.qp = qp_ctx->chroma_qp_map_jcbcr[ov_clip(slice_qp, 0, 63)] + qp_bd_offset + jcbcr_qp_offset;
    derive_dequant_ctx(ctudec, qp_ctx, 0);
}

static void
init_part_info(OVCTUDec *const ctudec, const struct OVPS *const prms)
{
    const OVSH *const sh = prms->sh;
    uint8_t slice_type = sh->sh_slice_type;
    if (slice_type == SLICE_I) {
        ctudec->part_ctx = &prms->sps_info.part_info[0];
        ctudec->part_ctx_c = &prms->sps_info.part_info[2];
    } else {
        ctudec->part_ctx = &prms->sps_info.part_info[1];
        ctudec->part_ctx_c = &prms->sps_info.part_info[3];
    }
}

static void
cabac_lines_uninit(OVSliceDec *sldec)
{
     struct CCLines *const lns = &sldec->cabac_lines[0];
     struct CCLines *const lns_c = &sldec->cabac_lines[1];

     /* FIXME avoid double free ensure it is allocated before
      * freeing it (in case of allocation failure.
      */
     ov_freep(&lns->qt_depth_map_x);
     ov_freep(&lns->log2_cu_w_map_x);
     ov_freep(&lns->cu_mode_x);

     ov_freep(&lns_c->qt_depth_map_x);
     ov_freep(&lns_c->log2_cu_w_map_x);
     ov_freep(&lns_c->cu_mode_x);
}


int
init_cabac_lines(OVSliceDec *sldec, const OVPS *const prms)
{
    uint8_t slice_type = sldec->slice_type;
    const OVPartInfo *const pinfo = slice_type == SLICE_I ? &prms->sps_info.part_info[0]
                                                          : &prms->sps_info.part_info[1];
     const struct TileInfo *const tinfo = &prms->pps_info.tile_info;

     struct CCLines *const lns   = &sldec->cabac_lines[0];
     struct CCLines *const lns_c = &sldec->cabac_lines[1];

     uint8_t log2_ctb_s = pinfo->log2_ctu_s;
     uint8_t log2_min_cb_s = pinfo->log2_min_cb_s;

     uint16_t pic_w = prms->sps->sps_pic_width_max_in_luma_samples;

     uint16_t nb_ctb_pic_w = (pic_w + ((1 << log2_ctb_s) - 1)) >> log2_ctb_s;
     uint16_t nb_pb_pic_w = nb_ctb_pic_w << (log2_ctb_s - log2_min_cb_s);

     lns->qt_depth_map_x  = ov_mallocz(sizeof(*lns->qt_depth_map_x)  * nb_pb_pic_w * tinfo->nb_tile_rows);
     lns->log2_cu_w_map_x = ov_mallocz(sizeof(*lns->log2_cu_w_map_x) * nb_pb_pic_w * tinfo->nb_tile_rows);
     lns->cu_mode_x       = ov_mallocz(sizeof(*lns->cu_mode_x)       * nb_pb_pic_w * tinfo->nb_tile_rows);

     lns_c->qt_depth_map_x  = ov_mallocz(sizeof(*lns_c->qt_depth_map_x)  * nb_pb_pic_w * tinfo->nb_tile_rows);
     lns_c->log2_cu_w_map_x = ov_mallocz(sizeof(*lns_c->log2_cu_w_map_x) * nb_pb_pic_w * tinfo->nb_tile_rows);
     lns_c->cu_mode_x       = ov_mallocz(sizeof(*lns_c->cu_mode_x)       * nb_pb_pic_w * tinfo->nb_tile_rows);

     lns->nb_pb_w   = nb_pb_pic_w;
     lns_c->nb_pb_w = nb_pb_pic_w;

     if (!lns->qt_depth_map_x || !lns->log2_cu_w_map_x || !lns->cu_mode_x ||
         !lns_c->qt_depth_map_x || !lns_c->log2_cu_w_map_x || !lns_c->cu_mode_x) {
         cabac_lines_uninit(sldec);
         return OVVC_ENOMEM;
     }

     return 0;
}


/* FIXME
 * reset according to entry info instead of whole line
 */
void
clear_cabac_lines(const OVSliceDec *sldec, const OVPS *const prms)
{
     uint8_t slice_type = sldec->slice_type;
     const OVPartInfo *const pinfo = slice_type == SLICE_I ? &prms->sps_info.part_info[0]
                                                           : &prms->sps_info.part_info[1];
     const struct TileInfo *const tinfo = &prms->pps_info.tile_info;

     const struct CCLines *const lns   = &sldec->cabac_lines[0];
     const struct CCLines *const lns_c = &sldec->cabac_lines[1];

     uint8_t log2_ctb_s    = pinfo->log2_ctu_s;
     uint8_t log2_min_cb_s = pinfo->log2_min_cb_s;

     /* TODO use active parameters such as generic pic info
      * see init_cabac_lines
      */
     const OVPPS *const pps = prms->pps;
     uint16_t pic_w = pps->pps_pic_width_in_luma_samples;
     uint16_t nb_ctb_pic_w = (pic_w + ((1 << log2_ctb_s) - 1)) >> log2_ctb_s;
     uint16_t nb_pb_pic_w = nb_ctb_pic_w << (log2_ctb_s - log2_min_cb_s);

     memset(lns->qt_depth_map_x,     0,  sizeof(*lns->qt_depth_map_x)  * nb_pb_pic_w * tinfo->nb_tile_rows);
     memset(lns->log2_cu_w_map_x, 0xFF,  sizeof(*lns->log2_cu_w_map_x) * nb_pb_pic_w * tinfo->nb_tile_rows);
     memset(lns->cu_mode_x,       0xFF,  sizeof(*lns->cu_mode_x)       * nb_pb_pic_w * tinfo->nb_tile_rows);

     memset(lns_c->qt_depth_map_x,     0,  sizeof(*lns_c->qt_depth_map_x)  * nb_pb_pic_w * tinfo->nb_tile_rows);
     memset(lns_c->log2_cu_w_map_x, 0xFF,  sizeof(*lns_c->log2_cu_w_map_x) * nb_pb_pic_w * tinfo->nb_tile_rows);
     memset(lns_c->cu_mode_x,       0xFF,  sizeof(*lns_c->cu_mode_x)       * nb_pb_pic_w * tinfo->nb_tile_rows);
}

static void
offset_cabac_lines(struct CCLines *const cc_lns, uint16_t ctb_x, uint8_t log2_ctb_s, uint8_t log2_min_cb_s)
{
     struct CCLines *const lns   = &cc_lns[0];
     struct CCLines *const lns_c = &cc_lns[1];

     int offset = ((uint32_t)ctb_x << log2_ctb_s) >> log2_min_cb_s;

     lns->qt_depth_map_x  += offset;
     lns->log2_cu_w_map_x += offset;
     lns->cu_mode_x       += offset;

     lns_c->qt_depth_map_x  += offset;
     lns_c->log2_cu_w_map_x += offset;
     lns_c->cu_mode_x       += offset;
}

/* FIXME
 * reset according to entry info instead of whole line
 */
static int
init_pic_border_info(struct RectEntryInfo *einfo, const OVPS *const prms, int entry_idx)
{
    /* Various info on pic_border */
    /* TODO check entry is border pic */
    /*TODO use derived value instead */
    const OVPPS *const pps = prms->pps;
    uint16_t pic_w = pps->pps_pic_width_in_luma_samples;
    uint16_t pic_h = pps->pps_pic_height_in_luma_samples;

    const OVSPS *const sps = prms->sps;
    uint8_t log2_ctb_s = sps->sps_log2_ctu_size_minus5 + 5;

    const int last_ctu_w = pic_w & ((1 << log2_ctb_s) - 1);
    const int last_ctu_h = pic_h & ((1 << log2_ctb_s) - 1);

    int nb_ctb_pic_w = (pic_w + ((1 << log2_ctb_s) - 1)) >> log2_ctb_s;
    int nb_ctb_pic_h = (pic_h + ((1 << log2_ctb_s) - 1)) >> log2_ctb_s;

    /* FIXME report an error when > nb_ctb_pic earlier */
    int pic_last_ctb_x = (einfo->ctb_x + einfo->nb_ctu_w) == nb_ctb_pic_w;
    int pic_last_ctb_y = (einfo->ctb_y + einfo->nb_ctu_h) == nb_ctb_pic_h;

    uint8_t full_ctb_w = (!pic_last_ctb_x || !last_ctu_w);
    uint8_t full_ctb_h = (!pic_last_ctb_y || !last_ctu_h);

    einfo->implicit_w = !full_ctb_w;
    einfo->implicit_h = !full_ctb_h;
    einfo->last_ctu_w = full_ctb_w ? (1 << log2_ctb_s) : last_ctu_w;
    einfo->last_ctu_h = full_ctb_h ? (1 << log2_ctb_s) : last_ctu_h;
    einfo->nb_ctb_pic_w = nb_ctb_pic_w;

    return 0;
}

void
slicedec_init_rect_entry(struct RectEntryInfo *einfo, const OVPS *const prms, int entry_idx)
{
    /* FIXME retrieve values from parameter sets */
    const struct SHInfo *const sh_info     = &prms->sh_info;
    const struct PPSInfo *const pps_info   = &prms->pps_info;
    const struct TileInfo *const tile_info = &pps_info->tile_info;

    int tile_x = entry_idx % tile_info->nb_tile_cols;
    int tile_y = entry_idx / tile_info->nb_tile_cols;

    einfo->tile_x = tile_x;
    einfo->tile_y = tile_y;

    einfo->nb_ctu_w = tile_info->nb_ctu_w[tile_x];
    einfo->nb_ctu_h = tile_info->nb_ctu_h[tile_y];

    einfo->nb_ctu_rect = einfo->nb_ctu_w * einfo->nb_ctu_h;
    einfo->entry_start = sh_info->rbsp_entry[entry_idx];
    einfo->entry_end   = sh_info->rbsp_entry[entry_idx + 1];

    einfo->ctb_x = tile_info->ctu_x[tile_x];
    einfo->ctb_y = tile_info->ctu_y[tile_y];

    /* FIXME test if need for init */
    einfo->ngh_flag = 0;
    einfo->implicit_h = 0;
    einfo->implicit_w = 0;

    init_pic_border_info(einfo, prms, entry_idx);
}

//TODOpar: temporary function, change with refs and ref_counts when functional
void
slicedec_copy_params(OVSliceDec *sldec, struct OVPS* dec_params)
{   
    struct OVPS* slice_params = &sldec->active_params;

    if (!slice_params->sps) {

        for (int i = 0; i < 8; i++) {
            slice_params->aps_alf[i] = ov_mallocz(sizeof(struct OVAPS));
        }

        slice_params->aps_alf_c = ov_mallocz(sizeof(struct OVAPS));
        slice_params->aps_cc_alf_cb = ov_mallocz(sizeof(struct OVAPS));
        slice_params->aps_cc_alf_cr = ov_mallocz(sizeof(struct OVAPS));
        slice_params->aps_lmcs = ov_mallocz(sizeof(struct OVAPS));
        slice_params->aps_scaling_list = ov_mallocz(sizeof(struct OVAPS));
    }

    #if 1
    hlsdata_unref(&slice_params->sps_ref);
    hlsdata_unref(&slice_params->pps_ref);
    hlsdata_unref(&slice_params->ph_ref);
    hlsdata_unref(&slice_params->sh_ref);
    #endif

    hlsdata_newref(&slice_params->sps_ref, dec_params->sps_ref);
    hlsdata_newref(&slice_params->pps_ref, dec_params->pps_ref);
    hlsdata_newref(&slice_params->ph_ref, dec_params->ph_ref);
    hlsdata_newref(&slice_params->sh_ref, dec_params->sh_ref);

    slice_params->sps = (OVSPS *)slice_params->sps_ref->data;
    slice_params->pps = (OVPPS *)slice_params->pps_ref->data;
    slice_params->ph = (OVPH *)slice_params->ph_ref->data;
    slice_params->sh = (OVSH *)slice_params->sh_ref->data;


    for (int i = 0; i < 8; i++) {
        if (dec_params->aps_alf[i]) {
            *(slice_params->aps_alf[i]) = *(dec_params->aps_alf[i]);
        }
    }

    if (dec_params->aps_alf_c) {
        *(slice_params->aps_alf_c) = *(dec_params->aps_alf_c);
    }

    if (dec_params->aps_cc_alf_cb) {
        *(slice_params->aps_cc_alf_cb) = *(dec_params->aps_cc_alf_cb);
    }

    if (dec_params->aps_cc_alf_cr) {
        *(slice_params->aps_cc_alf_cr) = *(dec_params->aps_cc_alf_cr);
    }

    if (dec_params->aps_lmcs) {
        *(slice_params->aps_lmcs) = *(dec_params->aps_lmcs);
    }

    if (dec_params->aps_scaling_list) {
        *(slice_params->aps_scaling_list) = *(dec_params->aps_scaling_list);
    }

    slice_params->sps_info = dec_params->sps_info;
    slice_params->pps_info = dec_params->pps_info;
    slice_params->ph_info = dec_params->ph_info;
    slice_params->sh_info = dec_params->sh_info;
    slice_params->pic_info = dec_params->pic_info;
}

//TODOpar: temporary function, change with refs and ref_counts when functional
void
slicedec_free_params(OVSliceDec *sldec)
{   
    struct OVPS* slice_params = &sldec->active_params;

    if(slice_params->sps){
    }

    if (slice_params->aps_alf_c) {
        for (int i = 0; i < 8; i++) {
            ov_freep(&(slice_params->aps_alf[i]));
        }
        ov_freep(&slice_params->aps_alf_c);
    }

    if (slice_params->aps_cc_alf_cb) {
        ov_freep(&slice_params->aps_cc_alf_cb);
        ov_freep(&slice_params->aps_cc_alf_cr);
    }

    if (slice_params->aps_lmcs) {
        ov_freep(&slice_params->aps_lmcs);
    }

    if (slice_params->aps_scaling_list) {
        ov_freep(&slice_params->aps_scaling_list);
    }
}

void
slicedec_finish_decoding(OVSliceDec *sldec)
{
    struct SliceSynchro *slice_sync = &sldec->slice_sync;

    /* There might be no NAL Unit attached to slicedec if
     * we failed before attaching NALU
     */
    if (slice_sync->slice_nalu) {
        ov_nalu_unref(&slice_sync->slice_nalu);
    }

    hlsdata_unref(&sldec->active_params.sps_ref);
    hlsdata_unref(&sldec->active_params.pps_ref);
    hlsdata_unref(&sldec->active_params.ph_ref);
    hlsdata_unref(&sldec->active_params.sh_ref);

    if (sldec->pic) {
        ov_log(NULL, OVLOG_DEBUG, "Decoder with POC %d, finished frame \n", sldec->pic->poc);

        ovdpb_report_decoded_frame( sldec->pic );
    }

    pthread_mutex_lock(&slice_sync->gnrl_mtx);
    slice_sync->active_state = DECODING_FINISHED;
    pthread_mutex_unlock(&slice_sync->gnrl_mtx);

    //Signal main thread that a slice thread is available
    struct MainThread* t_main = slice_sync->main_thread;
    if (t_main) {
        pthread_mutex_lock(&t_main->io_mtx);
        // ov_log(NULL, OVLOG_DEBUG,"Slice sign main\n");
        pthread_cond_signal(&t_main->io_cnd);
        pthread_mutex_unlock(&t_main->io_mtx);
    }
}

int
slicedec_decode_rect_entries(OVSliceDec *sldec, const OVPS *const prms, struct EntryThread* entry_th)
{
    /* FIXME do not recompute everywhere */
    int nb_entries = prms->pps_info.tile_info.nb_tile_cols *
                     prms->pps_info.tile_info.nb_tile_rows;

    int ret = 0;
    #if USE_THREADS
    ovthread_slice_add_entry_jobs(&sldec->slice_sync, slicedec_decode_rect_entry, nb_entries);
    #else
    int i;
    for (i = 0; i < nb_entries; ++i) {
        slicedec_update_entry_decoder(sldec, entry_th->ctudec);
        ret = slicedec_decode_rect_entry(sldec, entry_th->ctudec, prms, i);
    }
    slicedec_finish_decoding(sldec);
    ret = 0;
    #endif
    return ret;
}

static void
cabac_line_next_line(OVCTUDec *const ctudec, const struct CCLines *cclines)
{
    struct PartMap *const pmap_l = &ctudec->part_map;
    struct PartMap *const pmap_c = &ctudec->part_map_c;
    const struct CCLines *const lns_l = &cclines[0];
    const struct CCLines *const lns_c = &cclines[1];

    /* CTU Decoder part maps points to cabac_lines start */
    pmap_l->log2_cu_w_map_x = lns_l->log2_cu_w_map_x;
    pmap_l->qt_depth_map_x  = lns_l->qt_depth_map_x;
    pmap_l->cu_mode_x       = lns_l->cu_mode_x;

    pmap_c->log2_cu_w_map_x = lns_c->log2_cu_w_map_x;
    pmap_c->qt_depth_map_x  = lns_c->qt_depth_map_x;
    pmap_c->cu_mode_x       = lns_c->cu_mode_x;

    /* CTU Decoder part maps points to cabac_lines start */

    /* FIXME done twice on new entry see (reset lines function) */
    /* FIXME use nb_pb_ctb instead of size */
    memset(pmap_l->qt_depth_map_y,     0, sizeof(pmap_l->qt_depth_map_y));
    memset(pmap_l->cu_mode_y,       0xFF, sizeof(pmap_l->cu_mode_y));
    memset(pmap_l->log2_cu_h_map_y, 0xFF, sizeof(pmap_l->log2_cu_h_map_y));

    memset(pmap_c->qt_depth_map_y,     0, sizeof(pmap_c->qt_depth_map_y));
    memset(pmap_c->log2_cu_h_map_y, 0xFF, sizeof(pmap_c->log2_cu_h_map_y));
}

static void
cabac_line_next_ctu(OVCTUDec *const ctudec, uint8_t nb_pb_ctb)
{
    struct PartMap *const pmap_l = &ctudec->part_map;
    struct PartMap *const pmap_c = &ctudec->part_map_c;

    /* After each CTU move pointers to data corresponding
     * to next CTU
     */
    pmap_l->log2_cu_w_map_x += nb_pb_ctb;
    pmap_l->qt_depth_map_x  += nb_pb_ctb;
    pmap_l->cu_mode_x       += nb_pb_ctb;

    /*FIXME no diff between chroma / luma nb_pb */
    pmap_c->log2_cu_w_map_x += nb_pb_ctb;
    pmap_c->qt_depth_map_x  += nb_pb_ctb;
    pmap_c->cu_mode_x       += nb_pb_ctb;
}


/* Wrapper function around decode CTU calls so we can easily modify
 * what is to be done before and after each CTU
 * without adding many thing in each lin decoder
 */
static int
decode_ctu(OVCTUDec *const ctudec, const struct RectEntryInfo *const einfo,
           uint16_t ctb_addr_rs)
{
    uint8_t log2_ctb_s = ctudec->part_ctx->log2_ctu_s;
    int nb_ctu_w = einfo->nb_ctu_w;
    struct OVRCNCtx *rcn_ctx = &ctudec->rcn_ctx;
    int ret;

    ctudec->drv_ctx.inter_ctx.tmvp_avail = 0;
    /* FIXME pic border detection in neighbour flags ? CTU Neighbours
     * could be set according to upper level context
     */
    derive_ctu_neighborhood(ctudec, ctb_addr_rs, nb_ctu_w);

    ovcabac_read_ae_sao_ctu(ctudec, ctb_addr_rs, nb_ctu_w);
    ovcabac_read_ae_alf_ctu(ctudec, ctb_addr_rs, nb_ctu_w);
    ovcabac_read_ae_cc_alf_ctu(ctudec, ctb_addr_rs, nb_ctu_w);

    init_ctu_bitfield(rcn_ctx, ctudec->ctu_ngh_flags, log2_ctb_s);

    if ((ctb_addr_rs + 2) % nb_ctu_w == 0 && einfo->implicit_w) {
         uint64_t mask = ((uint64_t)1 << (((einfo->last_ctu_w + (1 << log2_ctb_s)) >> 2) + 1)) - 1;
         rcn_ctx->progress_field_c.hfield[0] &= mask;
         rcn_ctx->progress_field.hfield[0] &= mask;
    }

    ret = ctudec->coding_tree(ctudec, ctudec->part_ctx, 0, 0, log2_ctb_s, 0);

    ctudec->rcn_funcs.rcn_ctu_to_intra_line(rcn_ctx, ctb_addr_rs % nb_ctu_w << log2_ctb_s, log2_ctb_s);

    ctudec->rcn_funcs.rcn_write_ctu_to_frame(rcn_ctx, log2_ctb_s);

    /*FIXME call LMCS in write function instead */
    const struct OVBuffInfo *const fbuff = &rcn_ctx->frame_buff;
    ptrdiff_t stride_out_pic = fbuff->stride;
    OVSample *out_pic = fbuff->y;
    ctudec->rcn_funcs.lmcs_reshape_backward(out_pic, stride_out_pic, ctudec->lmcs_info.luts,
                                                    1 << log2_ctb_s, 1 << log2_ctb_s);

    if (!ctudec->dbf_disable) {
        uint8_t is_last_x = (ctb_addr_rs + 1) % nb_ctu_w == 0;
        uint8_t is_last_y = einfo->nb_ctu_h == (ctb_addr_rs / nb_ctu_w) + 1;
        #if 1

        ctudec->rcn_funcs.df.rcn_dbf_ctu(rcn_ctx, &ctudec->dbf_info, log2_ctb_s,
                    is_last_x, is_last_y);
                    #endif
    }

    return ret;
}

static int
decode_truncated_ctu(OVCTUDec *const ctudec, const struct RectEntryInfo *const einfo,
                     uint16_t ctb_addr_rs, int ctu_w, int ctu_h)
{
    uint8_t log2_ctb_s = ctudec->part_ctx->log2_ctu_s;
    struct OVRCNCtx *rcn_ctx = &ctudec->rcn_ctx;
    int nb_ctu_w = einfo->nb_ctu_w;
    int ret;

    ctudec->drv_ctx.inter_ctx.tmvp_avail = 0;
    /* FIXME pic border detection in neighbour flags ?*/
    derive_ctu_neighborhood(ctudec, ctb_addr_rs, nb_ctu_w);

    ovcabac_read_ae_sao_ctu(ctudec, ctb_addr_rs, nb_ctu_w);
    ovcabac_read_ae_alf_ctu(ctudec, ctb_addr_rs, nb_ctu_w);
    ovcabac_read_ae_cc_alf_ctu(ctudec, ctb_addr_rs, nb_ctu_w);

    /* FIXME pic border detection in neighbour flags ?*/
    init_ctu_bitfield_border(rcn_ctx, ctudec->ctu_ngh_flags, log2_ctb_s,
                             ctu_w, ctu_h);

    if ((ctb_addr_rs + 2) % nb_ctu_w == 0 && einfo->implicit_w) {
         uint64_t mask = ((uint64_t)1 << (((einfo->last_ctu_w + (1 << log2_ctb_s)) >> 2) + 1)) - 1;
         rcn_ctx->progress_field_c.hfield[0] &= mask;
         rcn_ctx->progress_field.hfield[0] &= mask;
    }

    ret = ctudec->coding_tree_implicit(ctudec, ctudec->part_ctx, 0, 0, log2_ctb_s,
                                       0, ctu_w, ctu_h);

    ctudec->rcn_funcs.rcn_ctu_to_intra_line(rcn_ctx, ctb_addr_rs % nb_ctu_w << log2_ctb_s, log2_ctb_s);

    ctudec->rcn_funcs.rcn_write_ctu_to_frame_border(rcn_ctx, ctu_w, ctu_h);

    const struct OVBuffInfo *const fbuff = &rcn_ctx->frame_buff;
    ptrdiff_t stride_out_pic = fbuff->stride;
    OVSample *out_pic = fbuff->y;
    ctudec->rcn_funcs.lmcs_reshape_backward(out_pic, stride_out_pic, ctudec->lmcs_info.luts,
                                                    ctu_w, ctu_h);

    if (!ctudec->dbf_disable) {
        uint8_t is_last_x = (ctb_addr_rs + 1) % nb_ctu_w == 0;
        uint8_t is_last_y = einfo->nb_ctu_h == (ctb_addr_rs / nb_ctu_w) + 1;
        ctudec->rcn_funcs.df.rcn_dbf_truncated_ctu(rcn_ctx, &ctudec->dbf_info, log2_ctb_s,
                              is_last_x, is_last_y, ctu_w, ctu_h);
    }

    return ret;
}

static int
decode_ctu_line(OVCTUDec *const ctudec, const OVSliceDec *const sldec,
                const struct DRVLines *const drv_lines,
                const struct RectEntryInfo *const einfo,
                uint16_t ctb_addr_rs)
{
    int nb_ctu_w = einfo->nb_ctu_w;
    struct OVRCNCtx *rcn_ctx = &ctudec->rcn_ctx;
    uint8_t log2_ctb_s = ctudec->part_ctx->log2_ctu_s;
    uint8_t log2_min_cb_s = ctudec->part_ctx->log2_min_cb_s;
    uint16_t nb_pb_ctb = (1 << log2_ctb_s) >> log2_min_cb_s;
    const uint8_t slice_type = sldec->slice_type;
    int ctb_x = 0;
    int ret;
    /* FIXME not really required ?*/
    uint8_t backup_qp = ctudec->drv_ctx.qp_map_x[0];

    ctudec->drv_ctx.inter_ctx.tmvp_ctx.ctu_w = 1 << log2_ctb_s;
    ctudec->drv_ctx.inter_ctx.tmvp_ctx.ctu_h = 1 << log2_ctb_s;

    ctudec->rcn_funcs.rcn_attach_ctu_buff(rcn_ctx, log2_ctb_s, 0);
    /* Do not copy on first line */
    if (ctb_addr_rs >= nb_ctu_w) {
        ctudec->rcn_funcs.rcn_intra_line_to_ctu(rcn_ctx, 0, log2_ctb_s);
    }

    while (ctb_x < nb_ctu_w - 1) {
        /*FIXME try to remove ctb_x computation */
        ctudec->ctb_x = einfo->ctb_x + ctb_x;

        ret = decode_ctu(ctudec, einfo, ctb_addr_rs);

        cabac_line_next_ctu(ctudec, nb_pb_ctb);

        /* Hackish way of keeping track of CTU last line
         * first QP to initialise delta qp for next ctu line
         * needs to be called before updating drv lines
         */
        if (ctb_x == 0) {
            backup_qp = ctudec->drv_ctx.qp_map_x[0];
        }

        ctudec->rcn_funcs.rcn_attach_ctu_buff(rcn_ctx, log2_ctb_s, ctb_x + 1);

        uint8_t ctb_addr_msk = (256 * 128 >> (2 * log2_ctb_s)) - 1;
        if (!((ctb_x + 1) & ctb_addr_msk)) {
            ctudec->rcn_funcs.rcn_update_ctu_border(rcn_ctx, log2_ctb_s);
        }
        ctudec->rcn_funcs.rcn_update_frame_buff(rcn_ctx, log2_ctb_s);

        if (slice_type != SLICE_I) {
            store_inter_maps(drv_lines, ctudec, ctb_x, 0);
        }
        if (ctudec->ibc_enabled) {
            store_ibc_maps(drv_lines, ctudec, ctb_x, 0);
        }

        if (!ctudec->dbf_disable) {
            const struct DBFLines *const dbf_lns = &drv_lines->dbf_lines;
            struct DBFInfo *const dbf_info = &ctudec->dbf_info;
            dbf_store_info(dbf_info, dbf_lns, log2_ctb_s, ctb_x);
            dbf_load_info(dbf_info, dbf_lns, log2_ctb_s, (ctb_x + 1) % nb_ctu_w);
        }

        ctb_addr_rs++;
        ctb_x++;

        /* FIXME
       * Move this somewhere else to avoid first line check
       */
        if (ctb_addr_rs >= nb_ctu_w) {
            // rcn_frame_line_to_ctu(&ctudec->rcn_ctx, log2_ctb_s);
            ctudec->rcn_funcs.rcn_intra_line_to_ctu(rcn_ctx, ctb_x << log2_ctb_s, log2_ctb_s);
        }
    }

    ctudec->ctb_x = einfo->ctb_x + ctb_x;

    /* last CTU require check on picture border for implicit splits*/
    if (!einfo->implicit_w) {

        ret = decode_ctu(ctudec, einfo, ctb_addr_rs);

        if (!ctudec->dbf_disable) {
            const struct DBFLines *const dbf_lns = &drv_lines->dbf_lines;
            struct DBFInfo *const dbf_info = &ctudec->dbf_info;
            dbf_store_info(dbf_info, dbf_lns, log2_ctb_s, ctb_x);
        }

    } else {
        int ctu_w = einfo->last_ctu_w;

        /* No horizontal implicit split since we would be
         * in decode_ctu_last_line otherwise
         */
        int ctu_h = 1 << log2_ctb_s;
        ctudec->drv_ctx.inter_ctx.tmvp_ctx.ctu_w = einfo->last_ctu_w;
        ctudec->drv_ctx.inter_ctx.tmvp_ctx.ctu_h = ctu_h;

        ret = decode_truncated_ctu(ctudec, einfo, ctb_addr_rs,
                                   ctu_w, ctu_h);

        if (!ctudec->dbf_disable) {
            const struct DBFLines *const dbf_lns = &drv_lines->dbf_lines;
            struct DBFInfo *const dbf_info = &ctudec->dbf_info;
            dbf_store_info(dbf_info, dbf_lns, log2_ctb_s, ctb_x);
        }
    }

    /* FIXME if inter only */
    if (slice_type != SLICE_I) {
        store_inter_maps(drv_lines, ctudec, ctb_x, 1);
    }

    if (ctudec->ibc_enabled) {
        store_ibc_maps(drv_lines, ctudec, ctb_x, 1);
    }
    
    //Apply in-loop filters on the available pixels of CTU line
    int ctb_y = ctudec->ctb_y - einfo->ctb_y;
    if(ctb_y == 0){
        ctudec->rcn_funcs.sao.rcn_sao_first_pix_rows(ctudec, einfo, ctb_y);
        if(einfo->nb_ctu_h == 1){
            ctudec->rcn_funcs.sao.rcn_sao_filter_line(ctudec, einfo, ctb_y);
            ctudec->rcn_funcs.alf.rcn_alf_filter_line(ctudec, einfo, ctb_y);
            ovdpb_report_decoded_ctu_line(sldec->pic, ctudec->ctb_y, einfo->ctb_x, einfo->ctb_x + nb_ctu_w - 1);
        } 
    }    
    else if(ctb_y == einfo->nb_ctu_h - 1){
        ctudec->rcn_funcs.sao.rcn_sao_filter_line(ctudec, einfo, ctb_y-1);
        ctudec->rcn_funcs.sao.rcn_sao_filter_line(ctudec, einfo, ctb_y);

        ctudec->rcn_funcs.alf.rcn_alf_filter_line(ctudec, einfo, ctb_y-1);
        ovdpb_report_decoded_ctu_line(sldec->pic, ctudec->ctb_y-1, einfo->ctb_x, einfo->ctb_x + einfo->nb_ctu_w - 1);
        ctudec->rcn_funcs.alf.rcn_alf_filter_line(ctudec, einfo, ctb_y);
        ovdpb_report_decoded_ctu_line(sldec->pic, ctudec->ctb_y, einfo->ctb_x, einfo->ctb_x + einfo->nb_ctu_w - 1);
    }
    else{
        ctudec->rcn_funcs.sao.rcn_sao_filter_line(ctudec, einfo, ctb_y-1);
        ctudec->rcn_funcs.alf.rcn_alf_filter_line(ctudec, einfo, ctb_y-1);
        ovdpb_report_decoded_ctu_line(sldec->pic, ctudec->ctb_y-1, einfo->ctb_x, einfo->ctb_x + nb_ctu_w - 1);
    }

    if (ctb_x == 0) {
        backup_qp = ctudec->drv_ctx.qp_map_x[0];
    }

    /* Next line will use the qp of the first pu as a start value
     * for qp_prediction
     */
    ctudec->qp_ctx.current_qp = backup_qp;

    ret = 0;

    /* FIXME Temporary error report on CABAC end of stream */
    if (ctudec->cabac_ctx->bytestream_end - ctudec->cabac_ctx->bytestream < -2) {
        ov_log(NULL, OVLOG_ERROR, "CABAC error diff end %d \n", ctb_addr_rs / nb_ctu_w);
    }

    return ret;
}

static int
decode_ctu_last_line(OVCTUDec *const ctudec, const OVSliceDec *const sldec,
                     const struct DRVLines *const drv_lines,
                     const struct RectEntryInfo *const einfo,
                     uint16_t ctb_addr_rs)
{
    int ret;

    struct OVRCNCtx *rcn_ctx = &ctudec->rcn_ctx;
    const int ctu_h = einfo->last_ctu_h;
    uint8_t log2_ctb_s = ctudec->part_ctx->log2_ctu_s;
    uint8_t log2_min_cb_s = ctudec->part_ctx->log2_min_cb_s;
    uint16_t nb_pb_ctb = (1 << log2_ctb_s) >> log2_min_cb_s;
    int nb_ctu_w = einfo->nb_ctu_w;
    uint8_t slice_type = sldec->slice_type;
    int ctb_x = 0;

    ctudec->rcn_funcs.rcn_attach_ctu_buff(rcn_ctx, log2_ctb_s, 0);

    ctudec->rcn_funcs.rcn_intra_line_to_ctu(rcn_ctx, 0, log2_ctb_s);

    ctudec->drv_ctx.inter_ctx.tmvp_ctx.ctu_w = 1 << log2_ctb_s;
    ctudec->drv_ctx.inter_ctx.tmvp_ctx.ctu_h = ctu_h;

    while (ctb_x < nb_ctu_w - 1) {
        uint8_t ctu_w = 1 << log2_ctb_s;

        ctudec->ctb_x = einfo->ctb_x + ctb_x;

        ret = decode_truncated_ctu(ctudec, einfo, ctb_addr_rs,
                                   ctu_w, ctu_h);

        cabac_line_next_ctu(ctudec, nb_pb_ctb);

        if (slice_type != SLICE_I) {
            store_inter_maps(drv_lines, ctudec, ctb_x, 0);
        }

        if (ctudec->ibc_enabled) {
            store_ibc_maps(drv_lines, ctudec, ctb_x, 0);
        }

        ctudec->rcn_funcs.rcn_attach_ctu_buff(rcn_ctx, log2_ctb_s, ctb_x + 1);

        uint8_t ctb_addr_msk = (256 * 128 >> (2 * log2_ctb_s)) - 1;
        if (!((ctb_x + 1) & ctb_addr_msk)) {
            ctudec->rcn_funcs.rcn_update_ctu_border(rcn_ctx, log2_ctb_s);
        }

        ctudec->rcn_funcs.rcn_update_frame_buff(rcn_ctx, log2_ctb_s);

        if (!ctudec->dbf_disable) {
            const struct DBFLines *const dbf_lns = &drv_lines->dbf_lines;
            struct DBFInfo *const dbf_info = &ctudec->dbf_info;
            dbf_store_info(dbf_info, dbf_lns, log2_ctb_s, ctb_x);
            dbf_load_info(dbf_info, dbf_lns, log2_ctb_s, (ctb_x + 1) % nb_ctu_w);
        }

        ctb_addr_rs++;
        ctb_x++;

        /* FIXME is first line check if only one line? */
        ctudec->rcn_funcs.rcn_intra_line_to_ctu(rcn_ctx, ctb_x << log2_ctb_s, log2_ctb_s);
    }

    ctudec->ctb_x = einfo->ctb_x + ctb_x;

    ctudec->drv_ctx.inter_ctx.tmvp_ctx.ctu_w = einfo->last_ctu_w;
    ctudec->drv_ctx.inter_ctx.tmvp_ctx.ctu_h = ctu_h;

    ret = decode_truncated_ctu(ctudec, einfo, ctb_addr_rs,
                               einfo->last_ctu_w, ctu_h);

    if (slice_type != SLICE_I) {
        store_inter_maps(drv_lines, ctudec, ctb_x, 1);
    }

    if (ctudec->ibc_enabled) {
        store_ibc_maps(drv_lines, ctudec, ctb_x, 1);
    }

    int ctb_y = ctudec->ctb_y - einfo->ctb_y;
    if(ctb_y == 0){
        ctudec->rcn_funcs.sao.rcn_sao_first_pix_rows(ctudec, einfo, ctb_y);
        ctudec->rcn_funcs.sao.rcn_sao_filter_line(ctudec, einfo, ctb_y);
        ctudec->rcn_funcs.alf.rcn_alf_filter_line(ctudec, einfo, ctb_y);
        ovdpb_report_decoded_ctu_line(sldec->pic, ctudec->ctb_y, einfo->ctb_x, einfo->ctb_x + nb_ctu_w - 1);
    }    
    else{
        ctudec->rcn_funcs.sao.rcn_sao_filter_line(ctudec, einfo, ctb_y-1);
        ctudec->rcn_funcs.sao.rcn_sao_filter_line(ctudec, einfo, ctb_y);

        ctudec->rcn_funcs.alf.rcn_alf_filter_line(ctudec, einfo, ctb_y-1);
        ovdpb_report_decoded_ctu_line(sldec->pic, ctudec->ctb_y-1, einfo->ctb_x, einfo->ctb_x + einfo->nb_ctu_w - 1);
        ctudec->rcn_funcs.alf.rcn_alf_filter_line(ctudec, einfo, ctb_y);
        ovdpb_report_decoded_ctu_line(sldec->pic, ctudec->ctb_y, einfo->ctb_x, einfo->ctb_x + einfo->nb_ctu_w - 1);
    }

    ret = 0;
    /* FIXME Temporary error report on CABAC end of stream */
    if (ctudec->cabac_ctx->bytestream_end - ctudec->cabac_ctx->bytestream < -2) {
        ov_log(NULL, OVLOG_ERROR, "CABAC error diff end %d \n", ctb_addr_rs / nb_ctu_w);
    }

    return ret;
}

static OVPicture *const
find_tmvp_collocated_ref(const OVSliceDec *const sldec, const OVPS *const ps)
{
    OVPicture *col_pic = NULL;
    const OVPPS *pps = ps->pps;
    const OVPH *ph = ps->ph;
    const OVSH *sh = ps->sh;

    uint8_t slice_type = sh->sh_slice_type;
    if(ph->ph_temporal_mvp_enabled_flag) {
        if (ph->ph_collocated_from_l0_flag || sh->sh_collocated_from_l0_flag || sh->sh_slice_type == SLICE_P) {
            int ref_idx = pps->pps_rpl_info_in_ph_flag ? ph->ph_collocated_ref_idx : sh->sh_collocated_ref_idx;
            col_pic = sldec->rpl0[ref_idx];
        } else if (sh->sh_slice_type != SLICE_I) {
            int ref_idx = pps->pps_rpl_info_in_ph_flag ? ph->ph_collocated_ref_idx : sh->sh_collocated_ref_idx;
            col_pic = sldec->rpl1[ref_idx];
        }
    }
    return col_pic;

}

static void
tmvp_entry_init(OVCTUDec *ctudec, const OVSliceDec *const sldec, const OVPS *const ps)
{
    /* FIXME try to remove ctu decoder reference from inter context */
    struct VVCTMVP *tmvp_ctx = &ctudec->drv_ctx.inter_ctx.tmvp_ctx;
    struct InterDRVCtx *inter_ctx = &ctudec->drv_ctx.inter_ctx;

    const OVPicture *const active_pic = sldec->pic;
    const OVPicture *const collocated_ref = find_tmvp_collocated_ref(sldec, ps);
    tmvp_ctx->col_ref = collocated_ref;

    ctudec->rcn_ctx.ctudec = ctudec;
    tmvp_ctx->ctudec = ctudec;

    tmvp_ctx->plane0 = &active_pic->mv_plane0;
    tmvp_ctx->plane1 = &active_pic->mv_plane1;

    tmvp_ctx->col_plane0 = collocated_ref ? &collocated_ref->mv_plane0 : NULL;
    tmvp_ctx->col_plane1 = collocated_ref ? &collocated_ref->mv_plane1 : NULL;

    /* FIXME used by other tools */
    memcpy(inter_ctx->dist_ref_0, sldec->dist_ref_0, sizeof(inter_ctx->dist_ref_0));
    memcpy(inter_ctx->dist_ref_1, sldec->dist_ref_1, sizeof(inter_ctx->dist_ref_1));

    memset(tmvp_ctx->dir_map_v0, 0, 34 * sizeof(uint64_t));
    memset(tmvp_ctx->dir_map_v1, 0, 34 * sizeof(uint64_t));
}

static void
init_lines(OVCTUDec *ctudec, const OVSliceDec *sldec, const struct RectEntryInfo *const einfo, const OVPS *const prms, const OVPartInfo*const part_ctx, struct DRVLines *drv_lines, struct CCLines *cc_lines)
{
    uint8_t log2_ctb_s    = part_ctx->log2_ctu_s;
    uint8_t log2_min_cb_s = part_ctx->log2_min_cb_s;

    uint32_t nb_ctb_pic_w = einfo->nb_ctb_pic_w;
    const struct TileInfo *tinfo = &prms->pps_info.tile_info;

    uint16_t ctb_offset = nb_ctb_pic_w * einfo->tile_y + einfo->ctb_x;

    offset_cabac_lines(cc_lines, ctb_offset, log2_ctb_s, log2_min_cb_s);

    *drv_lines = sldec->drv_lines;

    offset_drv_lines(drv_lines, einfo->tile_x, einfo->tile_y, einfo->ctb_x, log2_ctb_s,
                     log2_min_cb_s, tinfo->nb_tile_cols, nb_ctb_pic_w);

    cabac_line_next_line(ctudec, cc_lines);

    drv_line_next_line(ctudec, drv_lines);

}

static void
slicedec_smvd_params(OVCTUDec *const ctudec, const OVPS *const prms, int cur_poc)
{
    const struct InterDRVCtx *const inter_ctx = &ctudec->drv_ctx.inter_ctx;
    ctudec->smvd_enabled = 0;

    if (prms->sps->sps_smvd_enabled_flag && !inter_ctx->tmvp_ctx.ldc
        && !ctudec->mvd1_zero_enabled) {
        const int nb_active_ref0 = inter_ctx->nb_active_ref0;
        const int nb_active_ref1 = inter_ctx->nb_active_ref1;
        int ref = 0;
        int ref_idx0 = -1;
        int ref_idx1 = -1;
        int forw_poc = cur_poc;
        int back_poc = cur_poc;

        // search nearest forward POC in List 0
        for (ref = 0; ref < nb_active_ref0; ref++) {
            int ref_poc = inter_ctx->rpl0[ref]->poc;
            int ref_type = ST_REF;
            uint8_t is_lterm = (ref_type == LT_REF);
            if(ref_poc < cur_poc && (ref_poc > forw_poc || ref_idx0 == -1)  && !is_lterm) {
                forw_poc = ref_poc;
                ref_idx0 = ref;
            }
        }

        // search nearest backward POC in List 1
        for (ref = 0; ref < nb_active_ref1; ref++) {
            int ref_poc = inter_ctx->rpl1[ref]->poc;
            int ref_type = ST_REF;
            uint8_t is_lterm = (ref_type == LT_REF);
            if(ref_poc > cur_poc && (ref_poc < back_poc || ref_idx1 == -1)  && !is_lterm) {
                back_poc = ref_poc;
                ref_idx1 = ref;
            }
        }

        if (!(forw_poc < cur_poc && back_poc > cur_poc)) {
            forw_poc = cur_poc;
            back_poc = cur_poc;
            ref_idx0 = -1;
            ref_idx1 = -1;

            // search nearest backward POC in List 0
            for (ref = 0; ref < nb_active_ref0; ref++) {
                int ref_poc = inter_ctx->rpl0[ref]->poc;
                int ref_type = ST_REF;
                uint8_t is_lterm = (ref_type == LT_REF);
                if(ref_poc > cur_poc && (ref_poc < back_poc || ref_idx0 == -1)  && !is_lterm) {
                    back_poc = ref_poc;
                    ref_idx0 = ref;
                }
            }

            // search nearest forward POC in List 1
            for (ref = 0; ref < nb_active_ref1; ref++) {
                int ref_poc = inter_ctx->rpl1[ref]->poc;
                int ref_type = ST_REF;
                uint8_t is_lterm = (ref_type == LT_REF);
                if(ref_poc < cur_poc && (ref_poc > forw_poc || ref_idx1 == -1)  && !is_lterm) {
                    forw_poc = ref_poc;
                    ref_idx1 = ref;
                }
            }
        }

        if (forw_poc < cur_poc && back_poc > cur_poc){
            ctudec->smvd_enabled = 1;
            ctudec->drv_ctx.inter_ctx.ref_smvd_idx0 = ref_idx0;
            ctudec->drv_ctx.inter_ctx.ref_smvd_idx1 = ref_idx1;
        }
    }
}

static void
set_ref_weights_l0(struct WPInfo *wp_info,
                   const struct RPLWeightInfo *const wgt_info)
{
    int i;
    int16_t denom_l = 1 << wgt_info->luma_log2_weight_denom;
    int16_t denom_c = 1 << (wgt_info->luma_log2_weight_denom + wgt_info->delta_chroma_log2_weight_denom);

    for (i = 0; i < wgt_info->num_l0_weights; i++) {

        wp_info[i].flag = wgt_info->luma_weight_l0_flag[i];
        wp_info[i].weight_y = denom_l + wgt_info->delta_luma_weight_l0[i];
        wp_info[i].offset_y = wgt_info->luma_offset_l0[i] * (1 << 2);

        wp_info[i].flag_c = wgt_info->chroma_weight_l0_flag[i];
        wp_info[i].weight_cb = denom_c + wgt_info->delta_chroma_weight_l0[0][i];
        wp_info[i].weight_cr = denom_c + wgt_info->delta_chroma_weight_l0[1][i];

        int range = 128;
        int pred_cb  = (range - ((range * wp_info[i].weight_cb) >> (wgt_info->luma_log2_weight_denom + wgt_info->delta_chroma_log2_weight_denom)));
        int pred_cr  = (range - ((range * wp_info[i].weight_cr) >> (wgt_info->luma_log2_weight_denom + wgt_info->delta_chroma_log2_weight_denom)));
        wp_info[i].offset_cb = (pred_cb + wgt_info->delta_chroma_offset_l0[0][i]) * (1 << 2);
        wp_info[i].offset_cr = (pred_cr + wgt_info->delta_chroma_offset_l0[1][i]) * (1 << 2);
    }
}

static void
set_ref_weights_l1(struct WPInfo *wp_info,
                   const struct RPLWeightInfo *const wgt_info)
{
    int i;
    int16_t denom_l = 1 << wgt_info->luma_log2_weight_denom;
    int16_t denom_c = 1 << (wgt_info->luma_log2_weight_denom + wgt_info->delta_chroma_log2_weight_denom);

    for (i = 0; i < wgt_info->num_l1_weights; i++) {

        wp_info[i].flag = wgt_info->luma_weight_l1_flag[i];
        wp_info[i].weight_y = denom_l + wgt_info->delta_luma_weight_l1[i];
        wp_info[i].offset_y = wgt_info->luma_offset_l1[i] * (1 << 2);

        wp_info[i].flag_c = wgt_info->chroma_weight_l1_flag[i];
        wp_info[i].weight_cb = denom_c + wgt_info->delta_chroma_weight_l1[0][i];
        wp_info[i].weight_cr = denom_c + wgt_info->delta_chroma_weight_l1[1][i];

        int range = 128;
        int pred_cb  = (range - ((range * wp_info[i].weight_cb) >> (wgt_info->luma_log2_weight_denom + wgt_info->delta_chroma_log2_weight_denom)));
        int pred_cr  = (range - ((range * wp_info[i].weight_cr) >> (wgt_info->luma_log2_weight_denom + wgt_info->delta_chroma_log2_weight_denom)));
        wp_info[i].offset_cb = (pred_cb + wgt_info->delta_chroma_offset_l1[0][i]) * (1 << 2);
        wp_info[i].offset_cr = (pred_cr + wgt_info->delta_chroma_offset_l1[1][i]) * (1 << 2);
    }
}

static void
derive_ref_weights(OVCTUDec *const ctudec, const OVPS *const prms)
{
    struct InterDRVCtx *const inter_ctx = &ctudec->drv_ctx.inter_ctx;
    inter_ctx->weighted_pred_status = (prms->pps->pps_weighted_bipred_flag << 1) | prms->pps->pps_weighted_pred_flag;
    if (inter_ctx->weighted_pred_status) {
        const struct RPLWeightInfo *const wgh_info = prms->pps->pps_wp_info_in_ph_flag ? &prms->ph->wgt_info
                                                                                       : &prms->sh->wgt_info;

        set_ref_weights_l0(inter_ctx->wp_info0, wgh_info);
        set_ref_weights_l1(inter_ctx->wp_info1, wgh_info);
        inter_ctx->weighted_denom = wgh_info->luma_log2_weight_denom;
        inter_ctx->weighted_denom_c = wgh_info->luma_log2_weight_denom + wgh_info->delta_chroma_log2_weight_denom;
    }
}

#define RPR_NO_SCALE (1 << RPR_SCALE_BITS)
static void
ctudec_compute_refs_scaling(OVCTUDec *const ctudec, OVPicture *pic)
{
    struct Frame *frame = pic->frame;
    /*FIXME only 420*/
    uint16_t add_w = (pic->scale_info.scaling_win_left + pic->scale_info.scaling_win_right) << 1;
    uint16_t add_h = (pic->scale_info.scaling_win_top + pic->scale_info.scaling_win_bottom) << 1;
    uint16_t pic_w = frame->width  - add_w;
    uint16_t pic_h = frame->height - add_h;

    int scale_h, scale_v;
    int ref_pic_w, ref_pic_h;
    struct InterDRVCtx *const inter_ctx = &ctudec->drv_ctx.inter_ctx;
    inter_ctx->rpr_scale_msk0 = 0;
    inter_ctx->rpr_scale_msk1 = 0;
    for (int i = 0;  i < inter_ctx->nb_active_ref0; ++i){
        OVPicture *ref_pic  = inter_ctx->rpl0[i];
        frame = ref_pic->frame;
        if (!frame)
            continue;
        /*FIXME only 420*/
        add_w = (ref_pic->scale_info.scaling_win_left + ref_pic->scale_info.scaling_win_right) << 1;
        add_h = (ref_pic->scale_info.scaling_win_top + ref_pic->scale_info.scaling_win_bottom) << 1;
        ref_pic_w = frame->width  - add_w;
        ref_pic_h = frame->height - add_h;

        scale_h = ((ref_pic_w << RPR_SCALE_BITS) + (pic_w >> 1)) / pic_w;
        scale_v = ((ref_pic_h << RPR_SCALE_BITS) + (pic_h >> 1)) / pic_h;

        inter_ctx->scale_fact_rpl0[i][0] = scale_h;
        inter_ctx->scale_fact_rpl0[i][1] = scale_v;

        inter_ctx->rpr_scale_msk0 |= (scale_h != RPR_NO_SCALE || scale_v != RPR_NO_SCALE) << i;
    }
    for (int i = 0;  i < inter_ctx->nb_active_ref1; ++i){
        OVPicture *ref_pic = inter_ctx->rpl1[i];
        frame = ref_pic->frame;
        if (!frame)
            continue;
        /*FIXME only 420*/
        add_w = (ref_pic->scale_info.scaling_win_left + ref_pic->scale_info.scaling_win_right) << 1;
        add_h = (ref_pic->scale_info.scaling_win_top  + ref_pic->scale_info.scaling_win_bottom) << 1;
        ref_pic_w = frame->width  - add_w;
        ref_pic_h = frame->height - add_h;

        scale_h = ((ref_pic_w << RPR_SCALE_BITS) + (pic_w >> 1)) / pic_w;
        scale_v = ((ref_pic_h << RPR_SCALE_BITS) + (pic_h >> 1)) / pic_h;
        inter_ctx->scale_fact_rpl1[i][0] = scale_h;
        inter_ctx->scale_fact_rpl1[i][1] = scale_v;

        inter_ctx->rpr_scale_msk1 |= (scale_h != RPR_NO_SCALE || scale_v != RPR_NO_SCALE) << i;
    }
}

static int
slicedec_decode_rect_entry(OVSliceDec *sldec, OVCTUDec *const ctudec, const OVPS *const prms,
                           uint16_t entry_idx)
{
    int ctb_addr_rs = 0;
    int ctb_y = 0;
    int ret;

    struct RectEntryInfo einfo;
    struct OVRCNCtx *rcn_ctx = &ctudec->rcn_ctx;

    /*FIXME handle cabac alloc or keep it on the stack ? */
    OVCABACCtx cabac_ctx;
    slicedec_init_rect_entry(&einfo, prms, entry_idx);

    struct DRVLines drv_lines;
    struct CCLines cc_lines[2] = {sldec->cabac_lines[0], sldec->cabac_lines[1]};
    uint8_t log2_ctb_s = ctudec->part_ctx->log2_ctu_s;

    const int nb_ctu_w = einfo.nb_ctu_w;
    const int nb_ctu_h = einfo.nb_ctu_h;
    
    ctudec->cabac_ctx = &cabac_ctx;

    ctudec->qp_ctx.current_qp = ctudec->slice_qp;

    derive_dequant_ctx(ctudec, &ctudec->qp_ctx, 0);

    /*FIXME quick tmvp import */
    ctudec->nb_ctb_pic_w = einfo.nb_ctb_pic_w;

    tmvp_entry_init(ctudec, sldec, prms);
    ctudec->drv_ctx.inter_ctx.poc = sldec->pic->poc;
    /* FIXME tmp Reset DBF */
    memcpy(ctudec->drv_ctx.inter_ctx.rpl0, sldec->rpl0, sizeof(sldec->rpl0));
    memcpy(ctudec->drv_ctx.inter_ctx.rpl1, sldec->rpl1, sizeof(sldec->rpl1));

    ctudec->drv_ctx.inter_ctx.nb_active_ref0 = sldec->nb_active_refs0;
    ctudec->drv_ctx.inter_ctx.nb_active_ref1 = sldec->nb_active_refs1;

    ctudec_compute_refs_scaling(ctudec, sldec->pic);

    derive_ref_weights(ctudec, &sldec->active_params);

    ctudec->drv_ctx.inter_ctx.tmvp_ctx.ldc = 1;
    for (int i = 0; i < ctudec->drv_ctx.inter_ctx.nb_active_ref0; ++i) {
        uint8_t opp_ref_idx0 = 0xFF;

        if (sldec->dist_ref_0[i] < 0) {
            ctudec->drv_ctx.inter_ctx.tmvp_ctx.ldc = 0;
        }

        for (int j = 0; j < ctudec->drv_ctx.inter_ctx.nb_active_ref1; j ++) {
            if (ctudec->drv_ctx.inter_ctx.rpl0[i] == ctudec->drv_ctx.inter_ctx.rpl1[j]) {
                opp_ref_idx0 = j;
                break;
            }
        }
        ctudec->drv_ctx.inter_ctx.rpl0_opp[i] = opp_ref_idx0;
    }

    for (int i = 0; i < ctudec->drv_ctx.inter_ctx.nb_active_ref1; ++i) {
        uint8_t opp_ref_idx1 = 0xFF;

        if (sldec->dist_ref_1[i] < 0) {
            ctudec->drv_ctx.inter_ctx.tmvp_ctx.ldc = 0;
        }

        for (int j = 0; j < ctudec->drv_ctx.inter_ctx.nb_active_ref0; j ++) {
            if (ctudec->drv_ctx.inter_ctx.rpl1[i] == ctudec->drv_ctx.inter_ctx.rpl0[j]) {
                opp_ref_idx1 = j;
                break;
            }
        }
        ctudec->drv_ctx.inter_ctx.rpl1_opp[i] = opp_ref_idx1;
    }

    /* FIXME Bidir only */
    slicedec_smvd_params(ctudec, prms, sldec->pic->poc);

    /* FIXME entry might be check before attaching entry to CABAC so there
     * is no need for this check
     */
    ret = ovcabac_attach_entry(ctudec->cabac_ctx, einfo.entry_start, einfo.entry_end);
    if (ret < 0) {
        return OVVC_EINDATA;
    }

    /* FIXME Note cabac context tables could be initialised earlier
     * so we could only init once and recopy context tables to others
     * entries CABAC readers
     */
    ovcabac_init_slice_context_table(cabac_ctx.ctx_table, prms->sh->sh_slice_type ^ prms->sh->sh_cabac_init_flag,
                                     ctudec->slice_qp);

    init_lines(ctudec, sldec, &einfo, prms, ctudec->part_ctx,
               &drv_lines, cc_lines);

    ctudec->rcn_funcs.rcn_attach_ctu_buff(rcn_ctx, log2_ctb_s, 0);

    ctudec->rcn_funcs.rcn_attach_frame_buff(rcn_ctx, sldec->pic->frame, &einfo, log2_ctb_s);

    if (nb_ctu_w > ctudec->prev_nb_ctu_w_rect_entry) {
        int margin = 3;
        ctudec->rcn_funcs.rcn_alloc_filter_buffers(rcn_ctx, einfo.nb_ctu_w, margin, log2_ctb_s);
        ctudec->rcn_funcs.rcn_alloc_intra_line_buff(rcn_ctx, einfo.nb_ctu_w + 2, log2_ctb_s);
        ctudec->prev_nb_ctu_w_rect_entry = nb_ctu_w;
    }

    while (ctb_y < nb_ctu_h - 1) {

        ctudec->ctb_y = einfo.ctb_y + ctb_y;

        /* New ctu line */
        ret = decode_ctu_line(ctudec, sldec, &drv_lines, &einfo, ctb_addr_rs);

        cabac_line_next_line(ctudec, cc_lines);

        drv_line_next_line(ctudec, &drv_lines);

        ctudec->rcn_funcs.rcn_next_buff_line(rcn_ctx, log2_ctb_s);

        ctb_addr_rs += nb_ctu_w;
        ctb_y++;
    }

    ctudec->ctb_y = einfo.ctb_y + ctb_y;
    /* Last line */
    if (!einfo.implicit_h) {
        ret = decode_ctu_line(ctudec, sldec, &drv_lines, &einfo, ctb_addr_rs);
    } else {
        ret = decode_ctu_last_line(ctudec, sldec, &drv_lines, &einfo, ctb_addr_rs);
    }

    /*FIXME decide return value */
    return ctb_addr_rs;
}

static uint8_t ict_type(const OVPH *const ph)
{
    uint8_t type = (ph->ph_joint_cbcr_sign_flag << 1);
    type |= ph->ph_chroma_residual_scale_flag;
    return type;
}

static void
init_affine_status(OVCTUDec *const ctudec, const OVSPS *const sps,
                   const OVPH *const ph)
{
    ctudec->affine_nb_merge_cand = 5 - sps->sps_five_minus_max_num_subblock_merge_cand;
    ctudec->affine_status  = sps->sps_affine_amvr_enabled_flag;
    ctudec->affine_status |= sps->sps_6param_affine_enabled_flag << 1;
    ctudec->affine_status |= sps->sps_affine_prof_enabled_flag << 2;
    ctudec->drv_ctx.inter_ctx.affine_6params_enabled = sps->sps_6param_affine_enabled_flag;
}

static void
init_ladf(struct LADFParams *ladf, const OVSPS *sps)
{
    ladf->nb_intervals = sps->sps_num_ladf_intervals_minus2 + 2;
    ladf->qp_offset[0] = sps->sps_ladf_lowest_interval_qp_offset;
    ladf->threshold[0] = 0;
    for (int i = 1; i < ladf->nb_intervals; ++i) {
        ladf->qp_offset[i] = sps->sps_ladf_qp_offset[i - 1];
        ladf->threshold[i] = ladf->threshold[i - 1] + sps->sps_ladf_delta_threshold_minus1[i - 1] + 1;
    }
}

/* FIXME clean this init */
static int
slicedec_init_slice_tools(OVCTUDec *const ctudec, const OVPS *const prms)
{
    const OVSPS *const sps = prms->sps;
    const OVPPS *const pps = prms->pps;
    const OVSH *const sh = prms->sh;
    const OVPH *const ph = prms->ph;

    ctudec->max_log2_transform_skip_size = sps->sps_log2_transform_skip_max_size_minus2 + 2;

    /* FIXME dissociate SPS and SH/PH specific overrides to avoid  always resetting*/
    ctudec->enabled_mip   = sps->sps_mip_enabled_flag;
    ctudec->jcbcr_enabled = sps->sps_joint_cbcr_enabled_flag;
    ctudec->enable_lfnst  = sps->sps_lfnst_enabled_flag;
    ctudec->isp_enabled   = sps->sps_isp_enabled_flag;
    ctudec->enable_mrl    = sps->sps_mrl_enabled_flag;
    ctudec->bitdepth_minus8  = sps->sps_bitdepth_minus8;
    ctudec->bdpcm_enabled  = sps->sps_bdpcm_enabled_flag;

    ctudec->transform_skip_enabled   = sps->sps_transform_skip_enabled_flag;
    ctudec->sh_ts_disabled   = sh->sh_ts_residual_coding_disabled_flag;
    ctudec->max_num_merge_candidates = 6 - sps->sps_six_minus_max_num_merge_cand;
    ctudec->drv_ctx.inter_ctx.log2_parallel_merge_level = sps->sps_log2_parallel_merge_level_minus2 + 2;

    ctudec->delta_qp_enabled = pps->pps_cu_qp_delta_enabled_flag;
    ctudec->sbt_enabled      = sps->sps_sbt_enabled_flag;
    ctudec->affine_enabled   = sps->sps_affine_enabled_flag;
    ctudec->drv_ctx.inter_ctx.sbtmvp_enabled = sps->sps_sbtmvp_enabled_flag  && ph->ph_temporal_mvp_enabled_flag;
    ctudec->sbtmvp_enabled = sps->sps_sbtmvp_enabled_flag  && ph->ph_temporal_mvp_enabled_flag;
    ctudec->ibc_enabled = sps->sps_ibc_enabled_flag;
    ctudec->nb_ibc_cand_min1 = 6 - sps->sps_six_minus_max_num_ibc_merge_cand;
    ctudec->drv_ctx.ibc_ctx.bs1_map = &ctudec->dbf_info.bs1_map;
    ctudec->dbf_info.ibc_ctx = &ctudec->drv_ctx.ibc_ctx;

    if (ctudec->affine_enabled || sps->sps_sbtmvp_enabled_flag) {
        init_affine_status(ctudec, sps, ph);
    }

#if 1
    ctudec->dbf_disable = sh->sh_deblocking_filter_disabled_flag |
                          ph->ph_deblocking_filter_disabled_flag |
                          pps->pps_deblocking_filter_disabled_flag;
                          #else
                          ctudec->dbf_disable = 1;
#endif
    ctudec->dbf_info.beta_offset = pps->pps_luma_beta_offset_div2 * 2;
    ctudec->dbf_info.tc_offset   = pps->pps_luma_tc_offset_div2 * 2;
    if (pps->pps_chroma_tool_offsets_present_flag) {
        ctudec->dbf_info.beta_offset_cb = pps->pps_cb_beta_offset_div2 * 2;
        ctudec->dbf_info.tc_offset_cb   = pps->pps_cb_tc_offset_div2 * 2;
        ctudec->dbf_info.beta_offset_cr = pps->pps_cr_beta_offset_div2 * 2;
        ctudec->dbf_info.tc_offset_cr   = pps->pps_cr_tc_offset_div2 * 2;
    } else {
        ctudec->dbf_info.beta_offset_cb = ctudec->dbf_info.beta_offset;
        ctudec->dbf_info.tc_offset_cb   = ctudec->dbf_info.tc_offset;
        ctudec->dbf_info.beta_offset_cr = ctudec->dbf_info.beta_offset;
        ctudec->dbf_info.tc_offset_cr   = ctudec->dbf_info.tc_offset;
    }
    if (ph->ph_deblocking_params_present_flag) {
        ctudec->dbf_info.beta_offset = ph->ph_luma_beta_offset_div2 * 2;
        ctudec->dbf_info.tc_offset   = ph->ph_luma_tc_offset_div2 * 2;
        if (pps->pps_chroma_tool_offsets_present_flag) {
            ctudec->dbf_info.beta_offset_cb = ph->ph_cb_beta_offset_div2 * 2;
            ctudec->dbf_info.tc_offset_cb   = ph->ph_cb_tc_offset_div2 * 2;
            ctudec->dbf_info.beta_offset_cr = ph->ph_cr_beta_offset_div2 * 2;
            ctudec->dbf_info.tc_offset_cr   = ph->ph_cr_tc_offset_div2 * 2;
        } else {
            ctudec->dbf_info.beta_offset_cb = ctudec->dbf_info.beta_offset;
            ctudec->dbf_info.tc_offset_cb   = ctudec->dbf_info.tc_offset;
            ctudec->dbf_info.beta_offset_cr = ctudec->dbf_info.beta_offset;
            ctudec->dbf_info.tc_offset_cr   = ctudec->dbf_info.tc_offset;
        }
    }
    if (sh->sh_deblocking_params_present_flag) {
        ctudec->dbf_info.beta_offset = sh->sh_luma_beta_offset_div2 * 2;
        ctudec->dbf_info.tc_offset   = sh->sh_luma_tc_offset_div2 * 2;
        if (pps->pps_chroma_tool_offsets_present_flag) {
            ctudec->dbf_info.beta_offset_cb = sh->sh_cb_beta_offset_div2 * 2;
            ctudec->dbf_info.tc_offset_cb   = sh->sh_cb_tc_offset_div2 * 2;
            ctudec->dbf_info.beta_offset_cr = sh->sh_cr_beta_offset_div2 * 2;
            ctudec->dbf_info.tc_offset_cr   = sh->sh_cr_tc_offset_div2 * 2;
        } else {
            ctudec->dbf_info.beta_offset_cb = ctudec->dbf_info.beta_offset;
            ctudec->dbf_info.tc_offset_cb   = ctudec->dbf_info.tc_offset;
            ctudec->dbf_info.beta_offset_cr = ctudec->dbf_info.beta_offset;
            ctudec->dbf_info.tc_offset_cr   = ctudec->dbf_info.tc_offset;
        }
    }

    ctudec->lm_chroma_enabled = sps->sps_cclm_enabled_flag;

    ctudec->drv_ctx.inter_ctx.prof_enabled = sps->sps_affine_prof_enabled_flag  && !ph->ph_prof_disabled_flag;
    ctudec->drv_ctx.inter_ctx.bdof_enabled = sps->sps_bdof_enabled_flag && (!ph->ph_bdof_disabled_flag);
    ctudec->drv_ctx.inter_ctx.bdof_enabled &= sh->sh_slice_type == SLICE_B;
    ctudec->drv_ctx.inter_ctx.dmvr_enabled = sps->sps_dmvr_enabled_flag && (!ph->ph_dmvr_disabled_flag);
    ctudec->drv_ctx.inter_ctx.dmvr_enabled &= sh->sh_slice_type == SLICE_B;

    slice_init_qp_ctx(ctudec, prms);

    derive_dequant_ctx(ctudec, &ctudec->qp_ctx, 0);

    init_coding_coeff_coding_ctx(ctudec, prms);

    init_part_info(ctudec, prms);

    init_slice_tree_ctx(ctudec, prms);

    ctudec->drv_ctx.inter_ctx.mmvd_shift = ph->ph_mmvd_fullpel_only_flag << 1;
    ctudec->drv_ctx.inter_ctx.tmvp_enabled = ph->ph_temporal_mvp_enabled_flag;
    ctudec->mvd1_zero_enabled = ph->ph_mvd_l1_zero_flag;
    ctudec->drv_ctx.inter_ctx.tmvp_ctx.col_ref_l0 = ph->ph_collocated_from_l0_flag ||
                                                    sh->sh_collocated_from_l0_flag ||
                                                    sh->sh_slice_type == SLICE_P;

    ctudec->ciip_enabled = sps->sps_ciip_enabled_flag;
    ctudec->mmvd_enabled = sps->sps_mmvd_enabled_flag;
    ctudec->gpm_enabled  = sps->sps_gpm_enabled_flag;

    if (sps->sps_gpm_enabled_flag) {
        if (ctudec->max_num_merge_candidates >= 3) {
            ctudec->max_gpm_cand = ctudec->max_num_merge_candidates
                - sps->sps_max_num_merge_cand_minus_max_num_gpm_cand;
        } else if (ctudec->max_num_merge_candidates == 2) {
            ctudec->max_gpm_cand = 2;
        } else {
            ctudec->max_gpm_cand = 0;
        }
    }

    ctudec->bcw_enabled = !sps->sps_weighted_pred_flag && sps->sps_bcw_enabled_flag;
    ctudec->amvr_enabled = sps->sps_amvr_enabled_flag;
    ctudec->affine_amvr_enabled = sps->sps_affine_amvr_enabled_flag;

    rcn_init_functions(&ctudec->rcn_funcs, ict_type(ph), ctudec->lm_chroma_enabled,
                        sps->sps_chroma_vertical_collocated_flag, ph->ph_lmcs_enabled_flag,
                        sps->sps_bitdepth_minus8 + 8, sh->sh_dep_quant_used_flag);

    //In loop filter information for CTU reconstruction
    ctudec_init_in_loop_filters(ctudec, prms);
    ctudec->tmp_slice_type = sh->sh_slice_type;
    ctudec->scaling_list_enabled = ph->ph_explicit_scaling_list_enabled_flag || sh->sh_explicit_scaling_list_used_flag;
    ctudec->lfnst_scaling_list_enabled = ctudec->scaling_list_enabled && !sps->sps_scaling_matrix_for_lfnst_disabled_flag;

    if (ctudec->scaling_list_enabled) {
        uint8_t aps_id = ph->ph_scaling_list_aps_id;

        const OVAPS *aps = prms->aps_scaling_list;
        derive_scaling_tb_luts(&ctudec->tb_scaling_luts, aps);
    }

    if (sps->sps_ladf_enabled_flag) {
        init_ladf(&ctudec->dbf_info.ladf_prms, sps);
    } else {
        ctudec->dbf_info.ladf_prms.nb_intervals = 0; 
    }

    return 0;
}

int
slicedec_update_entry_decoder(OVSliceDec *sldec, OVCTUDec *ctudec)
{
    const OVPS *const prms = &sldec->active_params;
    ctudec->pic_w = prms->pps->pps_pic_width_in_luma_samples;
    ctudec->pic_h = prms->pps->pps_pic_height_in_luma_samples;
    slicedec_init_slice_tools(ctudec, prms);

    return 0;
}

int
slicedec_init_lines(OVSliceDec *const sldec, const OVPS *const prms)
{
    const OVSH *sh = prms->sh;
    sldec->slice_type = sh->sh_slice_type;

    if (!sldec->cabac_lines[0].qt_depth_map_x) {
        int ret;
        ret = init_cabac_lines(sldec, prms);
        if (ret < 0) {
            ov_log(NULL, 3, "FAILED init cabac lines\n");
            return ret;
        }
    } else {
        uint8_t slice_type = sldec->slice_type;
        const OVPartInfo *const pinfo = slice_type == SLICE_I ? &prms->sps_info.part_info[0]
            : &prms->sps_info.part_info[1];
        const struct TileInfo *const tinfo = &prms->pps_info.tile_info;

        uint8_t log2_ctb_s = pinfo->log2_ctu_s;
        uint8_t log2_min_cb_s = pinfo->log2_min_cb_s;
        uint16_t pic_w = prms->sps->sps_pic_width_max_in_luma_samples;
        uint16_t nb_ctb_pic_w = (pic_w + ((1 << log2_ctb_s) - 1)) >> log2_ctb_s;
        uint16_t nb_pb_pic_w = nb_ctb_pic_w << (log2_ctb_s - log2_min_cb_s);

        if (nb_pb_pic_w != sldec->cabac_lines[0].nb_pb_w) {
            cabac_lines_uninit(sldec);
            int ret = init_cabac_lines(sldec, prms);
            if (ret < 0) {
                ov_log(NULL, 3, "FAILED init cabac lines\n");
                return ret;
            }
        }
    }
    clear_cabac_lines(sldec, prms);

    if (!sldec->drv_lines.intra_luma_x) {
        int ret;
        ret = init_drv_lines(sldec, prms);
        if (ret < 0) {
            ov_log(NULL, 3, "FAILED init DRV lines\n");
            return ret;
        }
    } else {
        reset_drv_lines(sldec, prms);
    }

    return 0;
}

int
slicedec_init(OVSliceDec *sldec)
{
    int ret;

    sldec->slice_sync.owner = sldec;
    ret = ovthread_slice_sync_init(&sldec->slice_sync);
    if (ret < 0) {
        goto failthreads;
    }

    return 0;

failthreads:
    ov_log(NULL, OVLOG_ERROR, "Failed slice decoder initialisation\n");
    return OVVC_ENOMEM;
}

/*FIXME check is alloc in case of init failure*/
void
slicedec_uninit(OVSliceDec **sldec_p)
{
    OVSliceDec *sldec = *sldec_p;

    ovthread_slice_sync_uninit(&sldec->slice_sync);

    /*FIXME is init test */
    if (sldec->cabac_lines[0].log2_cu_w_map_x) {
        cabac_lines_uninit(sldec);
        drv_lines_uninit(sldec);
    }

    slicedec_free_params(sldec);

    ov_freep(sldec_p);

}

/* FIXME cleanup */
static void
derive_ctu_neighborhood(OVCTUDec *const ctudec,
                       int ctb_address, int nb_ctu_w)
{
    int is_left_border = ((ctb_address) % nb_ctu_w)  ? 0 : 1;
    int is_up_border   = (ctb_address < nb_ctu_w) ? 1 : 0;
    int is_right_border= ((ctb_address + 1) % nb_ctu_w);

    uint8_t ctb_flags = 0;

    /* FIXME clean */
    if (!is_left_border) {
        ctb_flags |= CTU_LFT_FLG;
    }
    if (!is_up_border) {
        ctb_flags |= CTU_UP_FLG;
    }
    if (!is_left_border && !is_up_border) {
        ctb_flags |= CTU_UPLFT_FLG;
    }
    if (!is_up_border && is_right_border) {
        ctb_flags |= CTU_UPRGT_FLG;
    }

    ctudec->ctu_ngh_flags = ctb_flags;
}

