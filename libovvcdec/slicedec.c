#include <string.h>
#include "ovdefs.h"
#include "libovvcutils/ovvcerror.h"
#include "libovvcutils/ovmem.h"
#include "libovvcutils/ovvcutils.h"
#include "slicedec.h"
#include "ctudec.h"
#include "nvcl_structures.h"
#include "dec_structures.h"
#include "vcl_cabac.h"
#include "vcl.h"

/* TODO define in a header */
enum SliceType {
     SLICE_B = 0,
     SLICE_P = 1,
     SLICE_I = 2
};

enum CTUNGHFlags
{
     CTU_UP_FLG    = 1 << 0,
     CTU_LFT_FLG   = 1 << 1,
     CTU_UPLFT_FLG = 1 << 2,
     CTU_UPRGT_FLG = 1 << 3,
};

struct RectEntryInfo {
    int nb_ctu_w;
    int nb_ctu_h;
    int nb_ctu_rect;
    const uint8_t *entry_start;
    const uint8_t *entry_end;
};

static int
slicedec_decode_rect_entry(OVSliceDec *sldec, const OVPS *const prms,
                           const struct RectEntryInfo *const einfo);

static void derive_dequant_ctx(OVCTUDec *const ctudec, const VVCQPCTX *const qp_ctx,
                               int cu_qp_delta);

static void derive_ctu_neighborhood(const OVSliceDec *const sldec,
                                    OVCTUDec *const ctudec,
                                    int ctb_address, int nb_ctu_w, int nb_ctu_h);
static void
init_slice_tree_ctx(OVCTUDec *const ctudec, const struct OVPS *prms)
{
    /* TODO use a specific structure for handling trees */
    const struct SPSInfo *sps_info = &prms->sps_info;
    const OVSPS *sps = prms->sps;
    const OVSH *sh = prms->sh;
    
    if (sh->sh_slice_type == SLICE_I) {
        if (sps->sps_qtbtt_dual_tree_intra_flag) {
            ctudec->coding_tree          = &dual_tree;
            ctudec->coding_tree_implicit = &dual_tree_implicit;
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
    ctudec->mts_enabled  = sps->sps_mts_enabled_flag && sps->sps_explicit_mts_intra_enabled_flag;
    ctudec->mts_implicit = sps->sps_mts_enabled_flag && !sps->sps_explicit_mts_intra_enabled_flag;

    if (sh->sh_dep_quant_used_flag) {
        ctudec->residual_coding = residual_coding_dpq;
        ctudec->residual_coding_chroma = residual_coding_chroma_dpq;
        ctudec->residual_coding_isp_h = residual_coding_isp_h_dpq;
        ctudec->residual_coding_isp_v = residual_coding_isp_v_dpq;
    } else {
        ctudec->residual_coding_isp_h = residual_coding_isp_h_sdh;
        ctudec->residual_coding_isp_v = residual_coding_isp_v_sdh;
        ctudec->residual_coding_chroma = residual_coding_chroma_sdh;
        ctudec->residual_coding = residual_coding_sdh;
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
    #if 0
    const OVSPS *const sps = prms->sps;
    #endif
    const OVPPS *const pps = prms->pps;
    const OVSH *const sh = prms->sh;
    const OVPH *const ph = prms->ph;
    VVCQPCTX *const qp_ctx = &ctudec->qp_ctx;

    /*FIXME check if not done in dec init */
    const uint8_t qp_bd_offset = 12;
    uint8_t pic_base_qp = pps->pps_init_qp_minus26 + 26;
    uint8_t pic_qp = pic_base_qp + ph->ph_qp_delta;
    int8_t slice_qp = pic_qp + sh->sh_qp_delta;
    int8_t cb_qp_offset = sh->sh_cb_qp_offset + pps->pps_cb_qp_offset;
    int8_t cr_qp_offset = sh->sh_cr_qp_offset + pps->pps_cr_qp_offset;
    int8_t jcbcr_qp_offset = sh->sh_joint_cbcr_qp_offset + pps->pps_joint_cbcr_qp_offset_value;

    ctudec->slice_qp = pic_qp + sh->sh_qp_delta;

    /* FIXME
     * check tables are valid when same qp table for all
     */
    qp_ctx->chroma_qp_map_cb    = prms->sps_info.qp_tables_c[0].qp;
    qp_ctx->chroma_qp_map_cr    = prms->sps_info.qp_tables_c[1].qp;
    qp_ctx->chroma_qp_map_jcbcr = prms->sps_info.qp_tables_c[2].qp;

    qp_ctx->current_qp = slice_qp;
    qp_ctx->cb_offset = cb_qp_offset;
    qp_ctx->cr_offset = cr_qp_offset;
    qp_ctx->jcbcr_offset = jcbcr_qp_offset;

    ctudec->dequant_luma.qp = slice_qp + qp_bd_offset;
    ctudec->dequant_cb.qp = qp_ctx->chroma_qp_map_cb[slice_qp + cb_qp_offset] + qp_bd_offset;
    ctudec->dequant_cr.qp = qp_ctx->chroma_qp_map_cr[slice_qp + cr_qp_offset] + qp_bd_offset;
    ctudec->dequant_joint_cb_cr.qp = qp_ctx->chroma_qp_map_jcbcr[slice_qp + jcbcr_qp_offset] + qp_bd_offset;
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
        ctudec->part_ctx_c = &prms->sps_info.part_info[2];
    }
}

static void
cabac_lines_uninit(OVSliceDec *sldec)
{
     struct CCLines *const lns = &sldec->cabac_lines[0];
     struct CCLines *const lns_c = &sldec->cabac_lines[1];

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
     const OVPartInfo *const pinfo = sldec->ctudec_list->part_ctx;
     const OVSPS *const sps = prms->sps;

     struct CCLines *const lns   = &sldec->cabac_lines[0];
     struct CCLines *const lns_c = &sldec->cabac_lines[1];

     uint8_t log2_ctb_s = pinfo->log2_ctu_s;
     uint8_t log2_min_cb_s = pinfo->log2_min_cb_s;
     /* TODO use active parameters such as generic pic info 
      * or something instead of this since this could be
      * overridden by PPS in case of sub_pic etc.
      * we could compute those values once for all earlier
      * in the decoding process
      */
     uint16_t pic_w = sps->sps_pic_width_max_in_luma_samples;

     uint16_t nb_ctb_pic_w = (pic_w + ((1 << log2_ctb_s) - 1)) >> log2_ctb_s;
     uint16_t nb_pb_pic_w = nb_ctb_pic_w << (log2_ctb_s - log2_min_cb_s);
     #if 0
     uint16_t pic_h = sps->sps_pic_height_max_in_luma_samples;
     uint16_t nb_ctb_pic_h = (pic_h + ((1 << log2_ctb_s) - 1)) >> log2_ctb_s;
     uint16_t nb_pb_pic_h = nb_ctb_pic_h << log2_min_cb_s;
     #endif


     lns->qt_depth_map_x  = ov_mallocz(sizeof(*lns->qt_depth_map_x) * nb_pb_pic_w);
     lns->log2_cu_w_map_x = ov_mallocz(sizeof(*lns->log2_cu_w_map_x) * nb_pb_pic_w);

     /*FIXME check if zero init value */
     lns->cu_mode_x       = ov_mallocz(sizeof(*lns->cu_mode_x) * nb_pb_pic_w);

     lns_c->qt_depth_map_x  = ov_mallocz(sizeof(*lns_c->qt_depth_map_x) * nb_pb_pic_w);
     lns_c->log2_cu_w_map_x = ov_mallocz(sizeof(*lns_c->log2_cu_w_map_x) * nb_pb_pic_w);

     /*FIXME check if zero init value */
     lns_c->cu_mode_x       = ov_mallocz(sizeof(*lns_c->cu_mode_x) * nb_pb_pic_w);

     if (!lns->qt_depth_map_x || !lns->log2_cu_w_map_x || !lns->cu_mode_x ||
         !lns_c->qt_depth_map_x || !lns_c->log2_cu_w_map_x || !lns_c->cu_mode_x) {
         cabac_lines_uninit(sldec);
         return OVVC_ENOMEM;
     }

     return 0;

}

void
reset_cabac_lines(OVSliceDec *sldec, const OVPS *const prms)
{
     const OVPartInfo *pinfo = sldec->ctudec_list->part_ctx;
     const OVSPS *const sps = prms->sps;
     struct CCLines *const lns = &sldec->cabac_lines[0];
     struct CCLines *const lns_c = &sldec->cabac_lines[1];
     uint8_t log2_ctb_s = pinfo->log2_ctu_s;
     uint8_t log2_min_cb_s = pinfo->log2_min_cb_s;
     /* TODO use active parameters such as generic pic info 
      * see init_cabac_lines
      */
     uint16_t pic_w = sps->sps_pic_width_max_in_luma_samples;
     uint16_t nb_ctb_pic_w = (pic_w + ((1 << log2_ctb_s) - 1)) >> log2_ctb_s;
     uint16_t nb_pb_pic_w = nb_ctb_pic_w << (log2_ctb_s - log2_min_cb_s);
     #if 0
     uint16_t nb_ctb_pic_h = (pic_h + ((1 << log2_ctb_s) - 1)) >> log2_ctb_s;
     uint16_t nb_pb_pic_h = nb_ctb_pic_h << log2_min_cb_s;
     #endif

     memset(lns->qt_depth_map_x,  0,  sizeof(*lns->qt_depth_map_x)  * nb_pb_pic_w);
     memset(lns->log2_cu_w_map_x, 0xFF,  sizeof(*lns->log2_cu_w_map_x) * nb_pb_pic_w);
     memset(lns->cu_mode_x,       0xFF,  sizeof(*lns->cu_mode_x)       * nb_pb_pic_w);

     memset(lns_c->qt_depth_map_x,  0,  sizeof(*lns_c->qt_depth_map_x)  * nb_pb_pic_w);
     memset(lns_c->log2_cu_w_map_x, 0xFF,  sizeof(*lns_c->log2_cu_w_map_x) * nb_pb_pic_w);
     memset(lns_c->cu_mode_x,       0xFF,  sizeof(*lns_c->cu_mode_x)       * nb_pb_pic_w);
}

int
slicedec_init_rect_entry(struct RectEntryInfo *entry, const OVPS *const prms, int entry_idx)
{
    
    /* FIXME retrieve values from parameter sets */
    const struct SHInfo *const sh_info= &prms->sh_info;
    entry->nb_ctu_w = 30;
    entry->nb_ctu_h = 17;
    entry->nb_ctu_rect = entry->nb_ctu_w * entry->nb_ctu_h;
    entry->entry_start = sh_info->rbsp_entry[entry_idx];
    entry->entry_end   = sh_info->rbsp_entry[entry_idx + 1];
    /* FIXME better way for last entry */

    #if 0
    struct LineCtxs line_ctx = sldec->line_ctx;
    #endif

    #if 0
    int i = 1;
    int tile_offset = 0;
    while (i <= tile_idx) {
        tile_offset += tile_ctx->nb_ctu_w[(i - 1) % num_tiles_cols] + 1;
        i++;
    }

    offset_line_map(&line_map, &line_ctx, tile_offset, log2_ctb_s);
    #endif

#if 0
    ctudec->last_ctb_alf_flag = 0;

    /*FIXME find a better way to handle tmvp reset */
    if (vvc_ctx->ps.sps_data->sps_temporal_mvp_enabled_flag) {
        ctudec->inter_ctx.tmvp_ctx.data_cur = vvc_ctx->active_pic->pic_mv_info->data;
        if (vvc_ctx->active_pic->collocated_ref)
            ctudec->inter_ctx.tmvp_ctx.data_ref = vvc_ctx->active_pic->collocated_ref->pic_mv_info->data;

        memset(&ctudec->inter_ctx.tmvp_ctx.tmvp_mv.mv_ctx0.map, 0, sizeof(ctudec->inter_ctx.tmvp_ctx.tmvp_mv.mv_ctx0.map));
        memset(&ctudec->inter_ctx.tmvp_ctx.tmvp_mv.mv_ctx1.map, 0, sizeof(ctudec->inter_ctx.tmvp_ctx.tmvp_mv.mv_ctx1.map));
    }
    memset(&ctudec->inter_ctx.mv_ctx0.map, 0, sizeof(ctudec->inter_ctx.mv_ctx0.map));
    memset(&ctudec->inter_ctx.mv_ctx1.map, 0, sizeof(ctudec->inter_ctx.mv_ctx1.map));

    memset(line_map.dbf_edge_ver, 0, sizeof(*line_map.dbf_edge_ver) * nb_ctu_w);
    memset(line_map.dbf_edge_hor, 0, sizeof(*line_map.dbf_edge_hor) * nb_ctu_w);

    init_frame_data(vvc_ctx->active_pic->frame, &ctudec->frame_data,
                    tile_ctx->ctu_x[tile_x], tile_ctx->ctu_y[tile_y], log2_ctb_s);
#endif
}

int
slicedec_decode_rect_entries(OVSliceDec *sldec, const OVPS *const prms)
{
    int nb_entries = 16;
    int i;
    for (i = 0; i < nb_entries; ++i) {
        struct RectEntryInfo entry;
        slicedec_init_rect_entry(&entry, prms, i);
        slicedec_decode_rect_entry(sldec, prms, &entry);
    }
}

static void
attach_cabac_lines(OVCTUDec *const ctudec, const OVSliceDec *const sldec)
{
    struct PartMap *const pmap_l = &ctudec->part_map;
    struct PartMap *const pmap_c = &ctudec->part_map_c;
    const struct CCLines *const lns_l = &sldec->cabac_lines[0];
    const struct CCLines *const lns_c = &sldec->cabac_lines[1];
    
    pmap_l->log2_cu_w_map_x = lns_l->log2_cu_w_map_x;
    pmap_l->qt_depth_map_x  = lns_l->qt_depth_map_x;
    pmap_l->cu_mode_x       = lns_l->cu_mode_x;

    pmap_c->log2_cu_w_map_x = lns_c->log2_cu_w_map_x;
    pmap_c->qt_depth_map_x  = lns_c->qt_depth_map_x;
    pmap_c->cu_mode_x       = lns_c->cu_mode_x;
}

static int
update_cabac_lines(OVCTUDec *const ctudec, const OVPS *const prms)
{
    const OVPartInfo *const pinfo = ctudec->part_ctx;

    uint8_t log2_ctb_s    = pinfo->log2_ctu_s;
    uint8_t log2_min_cb_s = pinfo->log2_min_cb_s;

    uint16_t nb_pb_ctb_w = (1 << log2_ctb_s) >> log2_min_cb_s;

    struct PartMap *const pmap_l = &ctudec->part_map;
    struct PartMap *const pmap_c = &ctudec->part_map_c;

    pmap_l->log2_cu_w_map_x += nb_pb_ctb_w;
    pmap_l->qt_depth_map_x  += nb_pb_ctb_w;
    pmap_l->cu_mode_x       += nb_pb_ctb_w;

    /*FIXME no diff between chroma / luma nb_pb */
    pmap_c->log2_cu_w_map_x += nb_pb_ctb_w;
    pmap_c->qt_depth_map_x  += nb_pb_ctb_w;
    pmap_c->cu_mode_x       += nb_pb_ctb_w;
}


static int
slicedec_decode_rect_entry(OVSliceDec *sldec, const OVPS *const prms,
                           const struct RectEntryInfo *const einfo)
{
    const struct SHInfo *const sh_info = &prms->sh_info;
    int ctb_addr_rs = 0;
    int ctb_y = 0;
    int ret;

    #if 0
    const struct VVCTileCTx *const tile_ctx = &vvc_ctx->tile_ctx;
    int num_tiles_cols = tile_ctx->nb_tile_cols;
    const int tile_x = tile_idx % num_tiles_cols;
    const int tile_y = tile_idx / num_tiles_cols;
    const int nb_ctu_w = tile_ctx->nb_ctu_w[tile_x];
    const int nb_ctu_h = tile_ctx->nb_ctu_h[tile_y];
    #endif

    /* FIXME handle more than one ctu dec */
    OVCTUDec *const ctudec = sldec->ctudec_list;
    /*FIXME handle cabac alloc or keep it on the stack ? */
    OVCABACCtx cabac_ctx;
    uint8_t log2_ctb_s = ctudec->part_ctx->log2_ctu_s;

    const int nb_ctu_w = einfo->nb_ctu_w;
    const int nb_ctu_h = einfo->nb_ctu_h;
    const int nb_ctu_rect = einfo->nb_ctu_rect;

    struct CCLines *line_map = sldec->cabac_lines;
    ctudec->cabac_ctx = &cabac_ctx;

    ctudec->qp_ctx.current_qp = ctudec->slice_qp;

    ret = ovcabac_attach_entry(ctudec->cabac_ctx, einfo->entry_start, einfo->entry_end);
    if (ret < 0) {
        return OVVC_EINDATA;
    }

    ovcabac_init_slice_context_table(cabac_ctx.ctx_table, prms->sh->sh_slice_type, ctudec->slice_qp);
    reset_cabac_lines(sldec, prms);
    #if 0
    attach_cabac_lines(ctudec, sldec);
    #if 1
    memset(ctudec->part_map.log2_cu_h_map_y, 0xFF, sizeof(ctudec->part_map.log2_cu_h_map_y));
    memset(ctudec->part_map_c.log2_cu_h_map_y, 0xFF, sizeof(ctudec->part_map_c.log2_cu_h_map_y));
    #endif
    #endif

    /* FIXME add a check for cabac end ?*/
    #if 0
    while (ctb_addr_rs < nb_ctu_rect) {
        #if 0
        int ctb_x = ctb_addr_rs % nb_ctu_w;
        int ctb_y = ctb_addr_rs / nb_ctu_w;
        #endif

        #if 0
        ctudec->ctb_x = tile_ctx->ctu_x[tile_x] + ctb_x;
        ctudec->ctb_y = tile_ctx->ctu_y[tile_y] + ctb_y;
        #endif

#if 0
        if (!ctb_x && ctb_y)
            ctudec->qp_ctx.current_qp = line_ctx.qp_x_map[0];
#endif

        /* FIXME pic border detection in neighbour flags ?*/
        derive_ctu_neighborhood(sldec, ctudec, ctb_addr_rs, nb_ctu_w, nb_ctu_h);

        #if 0
        load_ctu_ctx(ctudec, &line_map, &line_ctx, ctb_x, sldec->sh.slice_type);
        #endif

        if (1) {
            /*FIXME call ctudec */
            ret = ctudec->coding_tree(ctudec, ctudec->part_ctx, 0, 0, log2_ctb_s, 0);
        } else {
            #if 0
            ret = ctudec->coding_tree_implicit(ctudec, ctudec->part_ctx, 0, 0, log2_ctb_s,
                                               0, remaining_w, remaining_h);
                                               #endif
        }


#if 0
        if (!ctudec->dbf_disable)
            apply_dbf_ctu(ctudec, ctb_x, ctb_y, nb_ctu_w, nb_ctu_h);
#endif

        #if 0
        store_ctu_ctx(ctudec, &line_map, &line_ctx, ctb_x, sldec->sh.slice_type);
        #endif

        #if 0
        if (sldec->ps.sps_data->sps_temporal_mvp_enabled_flag) {
            store_tmvp(sldec, ctudec, &ctudec->inter_ctx.tmvp_ctx);
        }
        #endif

        ctb_addr_rs++;

        /*
        update_fd(&ctudec->frame_data, log2_ctb_s);
        if (!(ctb_addr_rs % nb_ctu_w)) {
            if (tile_y == num_tiles_cols - 1)
                if (vvc_ctx->threads_type & FF_THREAD_FRAME) {
                    ff_thread_report_progress(&vvc_ctx->active_pic->tf, ctb_y, 0);
                }
            update_fd_line(&ctudec->frame_data, nb_ctu_w, log2_ctb_s);
        }
        */
    }
    ret = ovcabac_end_of_slice(ctudec->cabac_ctx);
    #endif

    while (ctb_y < nb_ctu_h) {
        /* Start entry line */
        int ctb_x = 0;
        attach_cabac_lines(ctudec, sldec);
        #if 0
        reset_cabac_lines(sldec, prms);
        #endif
#if 1
        memset(ctudec->part_map.log2_cu_h_map_y, 0xFF, sizeof(ctudec->part_map.log2_cu_h_map_y));
        memset(ctudec->part_map_c.log2_cu_h_map_y, 0xFF, sizeof(ctudec->part_map_c.log2_cu_h_map_y));
#endif
        while (ctb_x < nb_ctu_w) {


#if 0
            if (!ctb_x && ctb_y)
                ctudec->qp_ctx.current_qp = line_ctx.qp_x_map[0];
#endif

            /* FIXME pic border detection in neighbour flags ?*/
            derive_ctu_neighborhood(sldec, ctudec, ctb_addr_rs, nb_ctu_w, nb_ctu_h);

#if 0
            load_ctu_ctx(ctudec, &line_map, &line_ctx, ctb_x, sldec->sh.slice_type);
#endif

            if (1) {
                /*FIXME call ctudec */
                ret = ctudec->coding_tree(ctudec, ctudec->part_ctx, 0, 0, log2_ctb_s, 0);
            } else {
#if 0
                ret = ctudec->coding_tree_implicit(ctudec, ctudec->part_ctx, 0, 0, log2_ctb_s,
                                                   0, remaining_w, remaining_h);
#endif
            }


#if 0
            if (!ctudec->dbf_disable)
                apply_dbf_ctu(ctudec, ctb_x, ctb_y, nb_ctu_w, nb_ctu_h);
#endif

#if 0
            store_ctu_ctx(ctudec, &line_map, &line_ctx, ctb_x, sldec->sh.slice_type);
#endif

#if 0
            if (sldec->ps.sps_data->sps_temporal_mvp_enabled_flag) {
                store_tmvp(sldec, ctudec, &ctudec->inter_ctx.tmvp_ctx);
            }
#endif

            ctb_addr_rs++;
            ctb_x++;
            update_cabac_lines(ctudec, prms);
        }
    }
    ctb_y++;

    /*FIXME decide return value */
    return ctb_addr_rs;
}

int
slicedec_init_slice_tools(OVSliceDec *const sldec, const OVPS *const prms)
{

    /* FIXME separate allocation outside of the scope of this function */
    OVCTUDec *const ctudec = sldec->ctudec_list;
    const OVSPS *const sps = prms->sps;
    const OVPPS *const pps = prms->pps;
    const OVSH *const sh = prms->sh;
    const OVPH *const ph = prms->ph;

    ctudec->max_log2_transform_skip_size = sps->sps_log2_transform_skip_max_size_minus2 + 2;

    #if 0
    ctudec->alf_num_chroma_alt = vvc_ctx->alf_num_alt_chroma;
    #endif
    /* FIXME dissociate SPS and SH/PH specific overrides to avoid  always resetting*/
    ctudec->lm_chroma_enabled = sps->sps_cclm_enabled_flag;
    ctudec->enabled_mip = sps->sps_mip_enabled_flag;

    ctudec->jcbcr_enabled = sps->sps_joint_cbcr_enabled_flag;
    ctudec->enable_lfnst  = sps->sps_lfnst_enabled_flag;
    ctudec->isp_enabled   = sps->sps_isp_enabled_flag;
    ctudec->enable_mrl    = sps->sps_mrl_enabled_flag;

    ctudec->transform_skip_enabled = sps->sps_transform_skip_enabled_flag;
    ctudec->max_num_merge_candidates = 6 - sps->sps_six_minus_max_num_merge_cand;

    ctudec->delta_qp_enabled = pps->pps_cu_qp_delta_enabled_flag;

    ctudec->dbf_disable = sh->sh_deblocking_filter_disabled_flag |
                          ph->ph_deblocking_filter_disabled_flag |
                          pps->pps_deblocking_filter_disabled_flag;

    slice_init_qp_ctx(ctudec, prms);

    derive_dequant_ctx(ctudec, &ctudec->qp_ctx, 0);

    init_coding_coeff_coding_ctx(ctudec, prms);

    init_part_info(ctudec, prms);

    init_slice_tree_ctx(ctudec, prms);

    /* Note it is important here that part info has already been set before calling
     * this function since it will be used to set line sizes*/

    /* FIXME
     * move this somewhere else so we can handle lines at a higher level
     */
    if (!sldec->cabac_lines[0].qt_depth_map_x) {
        int ret;
        ret = init_cabac_lines(sldec, prms);
        if (ret < 0) {
            ov_log(NULL, 3, "FAILED init cabac lines\n");
            return ret;
        }
    } else {
        reset_cabac_lines(sldec, prms);
    }

    return 0;
}

int
slicedec_init(OVSliceDec **dec_p)
{
    OVSliceDec *sldec;
    sldec = ov_mallocz(sizeof(OVSliceDec));
    if (!sldec) {
        return OVVC_ENOMEM;
    }

    *dec_p = sldec;


    #if 1
    ctudec_init(&sldec->ctudec_list);
    #endif

    return 0;
}

static void
uninit_ctudec_list(OVSliceDec *const sldec) {
     #if 0
     int nb_ctudec = sldec->nb_ctudec;
     int i;
     for (i = 0; i < nb_ctudec; ++i) {
         OVCTUDec *ctudec = &sldec->ctudec_list[i];
         ctudec_uninit(ctudec);
     }
     #endif
     ctudec_uninit(sldec->ctudec_list);
}

void
slicedec_uninit(OVSliceDec **sldec_p)
{
    OVSliceDec *sldec = *sldec_p;

    if (sldec->ctudec_list) {
        uninit_ctudec_list(sldec);
    }

    /*FIXME is init test */
    if (sldec->cabac_lines[0].log2_cu_w_map_x) {
        cabac_lines_uninit(sldec);
    }

    ov_freep(sldec_p);

}

#if 0
static int
init_slice_decoder(OVVCDec *const vvcdec, const OVNVCLCtx *const nvcl_ctx,
                   const OVSH *const sh)
{
    uint8_t first_slice_in_pic = sh->sh_slice_address == 0;
    if (first_slice_in_pic) {
        /* Activate ph parameters */
        activate_ph();

        if (drv) {
            #if 0
            init_pic_drv_ctx();

            init_pic_rcn_ctx();
            #endif
        }
    }

    init_slice_qp();

    override_slice_parameters();

    return 0;
}
#endif

static void
derive_ctu_neighborhood(const OVSliceDec *const sldec,
                       OVCTUDec *const ctudec,
                       int ctb_address, int nb_ctu_w, int nb_ctu_h)
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

#if 0
static int
ovdec_decode_slice()
{
    /* Init slice specific parameters */
    init_slice_parameters();

    /* init cabac */

    while (ctb_addr_rs < nb_ctb_slice){
        int ctb_x = ctb_addr_rs % nb_ctu_w;
        int ctb_y = ctb_addr_rs / nb_ctu_w;
        uint8_t new_line_flag = ctb_x == 0;

        if (!ctb_x && ctb_y)
            ctudec->qp_ctx.current_qp = line_ctx.qp_x_map[0];

        /* Derive border ctx */
        load_ctu_ctx(ctudec, &line_map, &line_ctx, ctb_x, vvc_ctx->sh.slice_type);

        ovdec_decode_ctu(vvc_ctx, ctudec, ctb_x, ctb_y, 0);

        /* Post ctu reconstruction */
        if (!ctudec->dbf_disable)
            apply_dbf_ctu(ctudec, ctb_x, ctb_y, nb_ctu_w, nb_ctu_h);

        store_ctu_ctx(ctudec, &line_map, &line_ctx, ctb_x, vvc_ctx->sh.slice_type);

        if (vvc_ctx->ps.sps_data->sps_temporal_mvp_enabled_flag) {
            store_tmvp(vvc_ctx, ctudec, &ctudec->inter_ctx.tmvp_ctx);
        }
        ctb_addr_rs++;

        update_fd(&ctudec->frame_data, ctudec->part_ctx->log2_ctu_s);
        if (!(ctb_addr_rs % nb_ctu_w)) {
            update_fd_line(&ctudec->frame_data, nb_ctu_w, ctudec->part_ctx->log2_ctu_s);
        }
    }
}
#endif

static void
derive_dequant_ctx(OVCTUDec *const ctudec, const VVCQPCTX *const qp_ctx,
                  int cu_qp_delta)
{
    /*FIXME avoid negative values especiallly in chroma_map derivation*/
    int base_qp = (qp_ctx->current_qp + cu_qp_delta + 64) & 63;
    ctudec->dequant_luma.qp = ((base_qp + 12) & 63);

    /*FIXME update transform skip ctx to VTM-10.0 */
    ctudec->dequant_luma_skip.qp = OVMAX(ctudec->dequant_luma.qp, qp_ctx->min_qp_prime_ts);
    ctudec->dequant_cb.qp = qp_ctx->chroma_qp_map_cb[(base_qp + qp_ctx->cb_offset + 64) & 63] + 12;
    ctudec->dequant_cr.qp = qp_ctx->chroma_qp_map_cr[(base_qp + qp_ctx->cr_offset + 64) & 63] + 12;
    ctudec->dequant_joint_cb_cr.qp = qp_ctx->chroma_qp_map_jcbcr[(base_qp + qp_ctx->jcbcr_offset + 64) & 63] + 12;
    ctudec->qp_ctx.current_qp = base_qp;
}

#if 0
static void
load_depth_maps(const struct LineCtxs *const l, OVCTUDec *const ctudec, int ctb_x)
{

    const uint8_t flags = ctudec->ctu_neighbour_flags;
    struct PartMap *part_map   = &ctudec->part_map;
    struct PartMap *part_map_c = &ctudec->part_map_c;
    int log2_nb_pb_ctu = 5 - (7 - ctudec->part_ctx->log2_ctu_s);
    int x_cb = ctb_x << log2_nb_pb_ctu;

    part_map->qt_depth_map_x  = &l->qt_depth_map_x_l [x_cb];
    part_map->log2_cu_w_map_x = &l->log2_cu_w_map_x_l[x_cb];

    part_map_c->qt_depth_map_x  = &l->qt_depth_map_x_c [x_cb];
    part_map_c->log2_cu_w_map_x = &l->log2_cu_w_map_x_c[x_cb];


    /* Reset on new line */
    if (!(flags&VVC_CTU_LEFT_FLAG)) {
        memset(part_map->qt_depth_map_y    , VVC_NA_DEPTH, sizeof(uint8_t) << log2_nb_pb_ctu);
        memset(part_map->log2_cu_h_map_y   , VVC_NA_DEPTH, sizeof(uint8_t) << log2_nb_pb_ctu);

        memset(part_map_c->qt_depth_map_y  , VVC_NA_DEPTH, sizeof(uint8_t) << log2_nb_pb_ctu);
        memset(part_map_c->log2_cu_h_map_y , VVC_NA_DEPTH, sizeof(uint8_t) << log2_nb_pb_ctu);
    }

    /* Reset on new tile / pic line */
    if (!(flags&VVC_CTU_UP_FLAG)) {
        memset(part_map->qt_depth_map_x    , VVC_NA_DEPTH, sizeof(uint8_t) << log2_nb_pb_ctu);
        memset(part_map->log2_cu_w_map_x   , VVC_NA_DEPTH, sizeof(uint8_t) << log2_nb_pb_ctu);

        memset(part_map_c->qt_depth_map_x  , VVC_NA_DEPTH, sizeof(uint8_t) << log2_nb_pb_ctu);
        memset(part_map_c->log2_cu_w_map_x , VVC_NA_DEPTH, sizeof(uint8_t) << log2_nb_pb_ctu);
    }
}

static void
offset_line_map(struct LineMaps *const line_map,
                struct LineCtxs *const line_ctx,
                int offset, int log2_ctu_s)
{
    int log2_nb_pb_ctu = 5 - (7 - log2_ctu_s);
    line_ctx->qt_depth_map_x_l  += offset << log2_nb_pb_ctu;
    line_ctx->qt_depth_map_x_c  += offset << log2_nb_pb_ctu;
    line_ctx->log2_cu_w_map_x_l += offset << log2_nb_pb_ctu;
    line_ctx->log2_cu_w_map_x_c += offset << log2_nb_pb_ctu;
}

static void
init_frame_data(const OVFrame *const frame, struct VVCFrameData *const fd,
                int ctb_x, int ctb_y, int log2_ctb_s)
{
    fd->data_y  = &frame->data[0] [(ctb_x << (log2_ctb_s + 1)) + frame->linesize[0] * (ctb_y << log2_ctb_s)];
    fd->data_cb = &frame->data[1] [(ctb_x << log2_ctb_s) + frame->linesize[1] * (ctb_y << (log2_ctb_s - 1))];
    fd->data_cr = &frame->data[2] [(ctb_x << log2_ctb_s) + frame->linesize[2] * (ctb_y << (log2_ctb_s - 1))];

    fd->ctb_stride   = frame->linesize[0] << log2_ctb_s;
    fd->ctb_stride_c = frame->linesize[1] << (log2_ctb_s - 1);
}

static void
update_fd_line(struct VVCFrameData *const fd,
               int nb_ctb_w, int log2_ctb_s)
{
    fd->data_y  += fd->ctb_stride - (nb_ctb_w << (log2_ctb_s + 1));
    fd->data_cb += fd->ctb_stride_c - (nb_ctb_w << log2_ctb_s);
    fd->data_cr += fd->ctb_stride_c - (nb_ctb_w << log2_ctb_s);
}

static void
update_fd(struct VVCFrameData *const fd, int log2_ctb_s)
{
    fd->data_y  += 1 << (log2_ctb_s + 1);
    fd->data_cb += 1 << log2_ctb_s;
    fd->data_cr += 1 << log2_ctb_s;
}
#endif

#if 0
static void
load_ctu_ctx(OVCTUDec *const ctudec,
             struct LineMaps *const line_map,
             struct LineCtxs *const line_ctx,
             int ctb_x, int slice_type)
{
    if (slice_type != SLICE_I) {
       /* specific Inter slice updates */
    }

    if (slice_type == SLICE_B) {
       /* specific B slice updates */
    }

    if (slice_type != SLICE_I) {
        #if 0
        load_inter_ngh_ctx(ctudec, line_ctx, ctb_x);
        #endif
    }

    load_depth_maps(line_ctx, ctudec, ctb_x);

    #if 0
    load_intra_ngh_ctx(ctudec, line_ctx, ctb_x);
    #endif
}

static void
store_ctu_ctx(OVCTUDec *const ctudec,
              struct LineMaps *const line_map,
              struct LineCtxs *const line_ctx,
              int ctb_x, int slice_type)
{
    store_availability_maps(line_ctx, line_map, ctudec, ctb_x);

    if (slice_type != VVC_SLICE_I) {
        store_inter_maps(line_ctx, line_map, ctudec, ctb_x);
    }

    if (!ctudec->dbf_disable) {
        store_edge_map(&ctudec->dbf_info, line_map, ctudec->part_ctx, ctb_x);
        store_bs_map(&ctudec->dbf_info, line_map, ctudec->part_ctx, ctb_x);
    }
}

static void
vvc_free_ctx_lines(VVCContext *const s)
{
    struct LineCtxs *const l = &s->line_ctx;
    struct LineMaps *const m = &s->line_map;
    ov_freep(&l->qt_depth_map_x_l);
    ov_freep(&l->log2_cu_w_map_x_l);

    ov_freep(&l->qt_depth_map_x_c);
    ov_freep(&l->log2_cu_w_map_x_c);
}

#endif
