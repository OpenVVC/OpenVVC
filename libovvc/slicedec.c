#include <string.h>
#include "ovdefs.h"
#include "overror.h"
#include "ovmem.h"
#include "ovutils.h"
#include "slicedec.h"
#include "ctudec.h"
#include "nvcl_structures.h"
#include "dec_structures.h"
#include "vcl_cabac.h"
#include "vcl.h"
#include "drv_utils.h"
#include "rcn.h"
#include "ovdpb.h"
#include "drv_lines.h"

/* TODO define in a header */
enum SliceType {
     SLICE_B = 0,
     SLICE_P = 1,
     SLICE_I = 2
};


struct RectEntryInfo {
    int ctb_x;
    int ctb_y;
    int nb_ctu_w;
    int nb_ctu_h;
    int nb_ctu_rect;
    const uint8_t *entry_start;
    const uint8_t *entry_end;
    uint8_t ngh_flag;
    uint8_t implicit_h;
    uint8_t implicit_w;
    int last_ctu_w;
    int last_ctu_h;
    int nb_ctb_pic_w;
};

static int
slicedec_decode_rect_entry(OVSliceDec *sldec, const OVPS *const prms,
                           const struct RectEntryInfo *const einfo);

static void derive_dequant_ctx(OVCTUDec *const ctudec, const VVCQPCTX *const qp_ctx,
                               int cu_qp_delta);

static void derive_ctu_neighborhood(const OVSliceDec *const sldec,
                                    OVCTUDec *const ctudec,
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


/* FIXME
 * reset according to entry info instead of whole line
 */
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

     memset(lns->qt_depth_map_x,     0,  sizeof(*lns->qt_depth_map_x)  * nb_pb_pic_w);
     memset(lns->log2_cu_w_map_x, 0xFF,  sizeof(*lns->log2_cu_w_map_x) * nb_pb_pic_w);
     memset(lns->cu_mode_x,       0xFF,  sizeof(*lns->cu_mode_x)       * nb_pb_pic_w);

     memset(lns_c->qt_depth_map_x,     0,  sizeof(*lns_c->qt_depth_map_x)  * nb_pb_pic_w);
     memset(lns_c->log2_cu_w_map_x, 0xFF,  sizeof(*lns_c->log2_cu_w_map_x) * nb_pb_pic_w);
     memset(lns_c->cu_mode_x,       0xFF,  sizeof(*lns_c->cu_mode_x)       * nb_pb_pic_w);
}

/* FIXME
 * reset according to entry info instead of whole line
 */
static int
init_pic_border_info(struct RectEntryInfo *einfo, const OVPS *const prms, int entry_idx)
{
    /* Various info on pic_border */
    /* TODO check entry is border pic */
    const OVSPS *const sps = prms->sps;
    /*TODO use derived value instead */
    uint16_t pic_w = sps->sps_pic_width_max_in_luma_samples;
    uint16_t pic_h = sps->sps_pic_height_max_in_luma_samples;
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

    #if 0
    int tile_offset = 0;
    int i = 1;

    while (i <= entry_idx) {
        tile_offset += pps_info->tile_info->nb_ctu_w[(i - 1) % num_tiles_cols] + 1;
        i++;
    }
    #endif

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

int
slicedec_decode_rect_entries(OVSliceDec *sldec, const OVPS *const prms)
{
    /* FIXME do not recompute everywhere */
    int nb_entries = prms->pps_info.tile_info.nb_tile_cols *
                     prms->pps_info.tile_info.nb_tile_rows;
    int ret = 0;
    int i;
    for (i = 0; i < nb_entries; ++i) {
        struct RectEntryInfo entry;
        slicedec_init_rect_entry(&entry, prms, i);
        ret = slicedec_decode_rect_entry(sldec, prms, &entry);
    }
    ret = 0;
    return ret;
}

static void
cabac_line_next_line(OVCTUDec *const ctudec, const OVSliceDec *const sldec)
{
    struct PartMap *const pmap_l = &ctudec->part_map;
    struct PartMap *const pmap_c = &ctudec->part_map_c;
    const struct CCLines *const lns_l = &sldec->cabac_lines[0];
    const struct CCLines *const lns_c = &sldec->cabac_lines[1];

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
    memset(pmap_l->cu_mode_y,       0xFF, sizeof(pmap_l->cu_mode_y));
    memset(pmap_l->qt_depth_map_y,  0xFF, sizeof(pmap_l->qt_depth_map_y));
    memset(pmap_l->log2_cu_h_map_y, 0xFF, sizeof(pmap_l->log2_cu_h_map_y));

    memset(pmap_c->qt_depth_map_y,  0xFF, sizeof(pmap_c->qt_depth_map_y));
    memset(pmap_c->log2_cu_h_map_y, 0xFF, sizeof(pmap_c->log2_cu_h_map_y));
}

static void
cabac_line_next_ctu(OVCTUDec *const ctudec, const OVPS *const prms)
{
    const OVPartInfo *const pinfo = ctudec->part_ctx;

    uint8_t log2_ctb_s    = pinfo->log2_ctu_s;
    uint8_t log2_min_cb_s = pinfo->log2_min_cb_s;

    uint16_t nb_pb_ctb_w = (1 << log2_ctb_s) >> log2_min_cb_s;

    struct PartMap *const pmap_l = &ctudec->part_map;
    struct PartMap *const pmap_c = &ctudec->part_map_c;

    /* After each CTU move pointers to data corresponding
     * to next CTU
     */
    pmap_l->log2_cu_w_map_x += nb_pb_ctb_w;
    pmap_l->qt_depth_map_x  += nb_pb_ctb_w;
    pmap_l->cu_mode_x       += nb_pb_ctb_w;

    /*FIXME no diff between chroma / luma nb_pb */
    pmap_c->log2_cu_w_map_x += nb_pb_ctb_w;
    pmap_c->qt_depth_map_x  += nb_pb_ctb_w;
    pmap_c->cu_mode_x       += nb_pb_ctb_w;
}

static void
tmvp_store_mv(OVCTUDec *ctudec)
{
    uint16_t ctb_x = ctudec->ctb_x;
    uint16_t ctb_y = ctudec->ctb_y;
    uint8_t log2_ctb_s = ctudec->part_ctx->log2_ctu_s;
    uint8_t log2_min_cb_s = ctudec->part_ctx->log2_min_cb_s;
    int nb_pb_ctb_w = (1 << log2_ctb_s) >> log2_min_cb_s;
    struct InterDRVCtx *const inter_ctx = &ctudec->drv_ctx.inter_ctx;

    struct MVPlane *plane0 = inter_ctx->tmvp_ctx.plane0;
    struct MVPlane *plane1 = inter_ctx->tmvp_ctx.plane1;

    int nb_ctb_w = ctudec->nb_ctb_pic_w;
    uint16_t ctb_addr_rs = ctb_x + ctb_y * nb_ctb_w;

    if (plane0->dirs) {
        uint64_t *dst_dirs = plane0->dirs + ctb_addr_rs * nb_pb_ctb_w;

        OVMV *dst_mv = plane0->mvs + ctb_x * nb_pb_ctb_w + (ctb_y * nb_pb_ctb_w* nb_pb_ctb_w) * nb_ctb_w;
        struct OVMVCtx *mv_ctx = &inter_ctx->mv_ctx0;
        int i;

        memcpy(dst_dirs, &mv_ctx->map.vfield[1], sizeof(uint64_t) * nb_pb_ctb_w);
        for (i = 0; i < nb_pb_ctb_w; ++i) {
            memcpy(dst_mv, &mv_ctx->mvs[1 + 34 * (i + 1)], sizeof(*dst_mv) * nb_pb_ctb_w);
            dst_mv += nb_pb_ctb_w * nb_ctb_w;
        }
    }

    if (plane1->dirs) {
        struct OVMVCtx *mv_ctx = &inter_ctx->mv_ctx1;
        uint64_t *dst_dirs = plane1->dirs + ctb_addr_rs * nb_pb_ctb_w;
        int i;

        OVMV *dst_mv = plane1->mvs + ctb_x * nb_pb_ctb_w + (ctb_y * nb_pb_ctb_w* nb_pb_ctb_w) * nb_ctb_w;

        /*FIXME memory could be spared with smaller map size when possible */
        memcpy(dst_dirs, &mv_ctx->map.vfield[1], sizeof(uint64_t) * nb_pb_ctb_w);
        for (i = 0; i < nb_pb_ctb_w; ++i) {
            memcpy(dst_mv, &mv_ctx->mvs[1 + 34 * (i + 1)], sizeof(*dst_mv) * nb_pb_ctb_w);
            dst_mv += nb_pb_ctb_w * nb_ctb_w;
        }

    }

}

/* Wrapper function around decode CTU calls so we can easily modify
 * what is to be done before and after each CTU
 * without adding many thing in each lin decoder
 */
static int
decode_ctu(OVCTUDec *const ctudec, const OVSliceDec *const sldec,
           const OVPS *const prms, const struct RectEntryInfo *const einfo,
           uint16_t ctb_addr_rs)
{
    uint8_t log2_ctb_s = ctudec->part_ctx->log2_ctu_s;
    int nb_ctu_w = einfo->nb_ctu_w;
    int ret;

    ctudec->drv_ctx.inter_ctx.tmvp_avail = 0;
    /* FIXME pic border detection in neighbour flags ? CTU Neighbours
     * could be set according to upper level context
     */
    derive_ctu_neighborhood(sldec, ctudec, ctb_addr_rs, nb_ctu_w);

    init_ctu_bitfield(&ctudec->rcn_ctx, ctudec->ctu_ngh_flags, log2_ctb_s);

    ret = ctudec->coding_tree(ctudec, ctudec->part_ctx, 0, 0, log2_ctb_s, 0);

    rcn_write_ctu_to_frame(&ctudec->rcn_ctx, log2_ctb_s);

    if (!ctudec->dbf_disable) {
        uint8_t is_last_x = (ctb_addr_rs + 1) % nb_ctu_w == 0;
        uint8_t is_last_y = einfo->nb_ctu_h == (ctb_addr_rs / nb_ctu_w) + 1;
        rcn_dbf_ctu(&ctudec->rcn_ctx, &ctudec->dbf_info, log2_ctb_s,
                    is_last_x, is_last_y);
    }

    tmvp_store_mv(ctudec);

    return ret;
}

static int
decode_ctu_implicit(OVCTUDec *const ctudec, const OVSliceDec *const sldec,
                    const OVPS *const prms, const struct RectEntryInfo *const einfo,
                    uint16_t ctb_addr_rs, int ctu_w, int ctu_h)
{
    uint8_t log2_ctb_s = ctudec->part_ctx->log2_ctu_s;
    int nb_ctu_w = einfo->nb_ctu_w;
    int ret;

    ctudec->drv_ctx.inter_ctx.tmvp_avail = 0;
    /* FIXME pic border detection in neighbour flags ?*/
    derive_ctu_neighborhood(sldec, ctudec, ctb_addr_rs, nb_ctu_w);

    /* FIXME pic border detection in neighbour flags ?*/
    init_ctu_bitfield_border(&ctudec->rcn_ctx, ctudec->ctu_ngh_flags, log2_ctb_s,
                             ctu_w, ctu_h);

    ret = ctudec->coding_tree_implicit(ctudec, ctudec->part_ctx, 0, 0, log2_ctb_s,
                                       0, ctu_w, ctu_h);

    rcn_write_ctu_to_frame_border(&ctudec->rcn_ctx,
                                  ctu_w, ctu_h);

    if (!ctudec->dbf_disable) {
        uint8_t is_last_x = (ctb_addr_rs + 1) % nb_ctu_w == 0;
        uint8_t is_last_y = einfo->nb_ctu_h == (ctb_addr_rs / nb_ctu_w) + 1;
        rcn_dbf_truncated_ctu(&ctudec->rcn_ctx, &ctudec->dbf_info, log2_ctb_s,
                              is_last_x, is_last_y, ctu_w, ctu_h);
    }

    tmvp_store_mv(ctudec);

    return ret;
}

static int
decode_ctu_line(OVCTUDec *const ctudec, const OVSliceDec *const sldec,
                const OVPS *const prms, const struct RectEntryInfo *const einfo,
                uint16_t ctb_addr_rs)
{
    int nb_ctu_w = einfo->nb_ctu_w;
    int ctb_x = 0;
    int ret;
    uint8_t backup_qp;

    drv_line_next_ctu(ctudec, sldec, &sldec->drv_lines, prms, 0);

    /* Do not copy on first line */
    if (ctb_addr_rs >= nb_ctu_w) {
        rcn_frame_line_to_ctu(&ctudec->rcn_ctx, ctudec->part_ctx->log2_ctu_s);
    }

    while (ctb_x < nb_ctu_w - 1) {
        /*FIXME try to remove ctb_x computation */
        ctudec->ctb_x = einfo->ctb_x + ctb_x;

        ret = decode_ctu(ctudec, sldec, prms, einfo, ctb_addr_rs);

        cabac_line_next_ctu(ctudec, prms);

        /* Hackish way of keeping track of CTU last line
         * first QP to initialise delta qp for next ctu line
         * needs to be called before updating drv lines
         */
        if (ctb_x == 0) {
            backup_qp = ctudec->drv_ctx.qp_map_x[0];
        }

        drv_line_next_ctu(ctudec, sldec, &sldec->drv_lines, prms, ctb_x);

        rcn_update_ctu_border(&ctudec->rcn_ctx, ctudec->part_ctx->log2_ctu_s);

        store_inter_maps(&sldec->drv_lines, ctudec, ctb_x);
        uint8_t log2_ctb_s = ctudec->part_ctx->log2_ctu_s;
        dbf_store_info(&ctudec->dbf_info, &sldec->drv_lines.dbf_lines, log2_ctb_s, ctb_x);

        dbf_load_info(&ctudec->dbf_info, &sldec->drv_lines.dbf_lines, log2_ctb_s, (ctb_x + 1) % nb_ctu_w);
        /* FIXME
         * Move this somewhere else to avoid first line check
         */
        if (ctb_addr_rs >= nb_ctu_w) {
            rcn_frame_line_to_ctu(&ctudec->rcn_ctx, ctudec->part_ctx->log2_ctu_s);
        }

        ctb_addr_rs++;
        ctb_x++;
    }

    ctudec->ctb_x = einfo->ctb_x + ctb_x;

    /* last CTU require check on picture border for implicit splits*/
    if (!einfo->implicit_w) {

        uint8_t log2_ctb_s = ctudec->part_ctx->log2_ctu_s;
        ret = decode_ctu(ctudec, sldec, prms, einfo, ctb_addr_rs);

        dbf_store_info(&ctudec->dbf_info, &sldec->drv_lines.dbf_lines, log2_ctb_s, ctb_x);

    } else {
        uint8_t log2_ctb_s = ctudec->part_ctx->log2_ctu_s;
        int ctu_w = einfo->last_ctu_w;

        /* No horizontal implicit split since we would be
         * in decode_ctu_last_line otherwise
         */
        int ctu_h = 1 << log2_ctb_s;

        ret = decode_ctu_implicit(ctudec, sldec, prms, einfo, ctb_addr_rs,
                                  ctu_w, ctu_h);

        dbf_store_info(&ctudec->dbf_info, &sldec->drv_lines.dbf_lines, log2_ctb_s, ctb_x);
    }

    /* Next line will use the qp of the first pu as a start value
     * for qp_prediction
     */
    /* FIXME if inter only */
    store_inter_maps(&sldec->drv_lines, ctudec, ctb_x);

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
                     const OVPS *const prms, const struct RectEntryInfo *const einfo,
                     uint16_t ctb_addr_rs)
{
    int ret;

    const int ctu_h = einfo->last_ctu_h;
    uint8_t log2_ctb_s = ctudec->part_ctx->log2_ctu_s;
    int nb_ctu_w = einfo->nb_ctu_w;
    int ctb_x = 0;

    drv_line_next_ctu(ctudec, sldec, &sldec->drv_lines, prms, 0);

    rcn_frame_line_to_ctu(&ctudec->rcn_ctx, log2_ctb_s);

    while (ctb_x < nb_ctu_w - 1) {
        uint8_t ctu_w = 1 << log2_ctb_s;

        ctudec->ctb_x = einfo->ctb_x + ctb_x;

        ret = decode_ctu_implicit(ctudec, sldec, prms, einfo, ctb_addr_rs,
                                  ctu_w, ctu_h);

        cabac_line_next_ctu(ctudec, prms);

        drv_line_next_ctu(ctudec, sldec, &sldec->drv_lines, prms, ctb_x);

        store_inter_maps(&sldec->drv_lines, ctudec, ctb_x);

        rcn_update_ctu_border(&ctudec->rcn_ctx, log2_ctb_s);

        /* FIXME is first line check if only one line? */
        rcn_frame_line_to_ctu(&ctudec->rcn_ctx, log2_ctb_s);

        dbf_store_info(&ctudec->dbf_info, &sldec->drv_lines.dbf_lines, log2_ctb_s, ctb_x);

        dbf_load_info(&ctudec->dbf_info, &sldec->drv_lines.dbf_lines, log2_ctb_s, (ctb_x + 1) % nb_ctu_w);

        ctb_addr_rs++;
        ctb_x++;
    }

    ctudec->ctb_x = einfo->ctb_x + ctb_x;

    ret = decode_ctu_implicit(ctudec, sldec, prms, einfo, ctb_addr_rs,
                              einfo->last_ctu_w, ctu_h);

    ret = 0;
    /* FIXME Temporary error report on CABAC end of stream */
    if (ctudec->cabac_ctx->bytestream_end - ctudec->cabac_ctx->bytestream < -2) {
        ov_log(NULL, OVLOG_ERROR, "CABAC error diff end %d \n", ctb_addr_rs / nb_ctu_w);
    }

    return ret;
}

static void
slicedec_attach_frame_buff(OVCTUDec *const ctudec, OVSliceDec *sldec,
                           const struct RectEntryInfo *const einfo)
{
    OVFrame *f = sldec->pic->frame;
    uint8_t log2_ctb_s = ctudec->part_ctx->log2_ctu_s;
    struct OVBuffInfo *const fbuff = &ctudec->rcn_ctx.frame_buff;
    uint32_t entry_start_offset = ((uint32_t)einfo->ctb_x << (log2_ctb_s + 1));
    uint32_t entry_start_offset_c = ((uint32_t)einfo->ctb_x << (log2_ctb_s));

    entry_start_offset   += ((uint32_t)einfo->ctb_y << log2_ctb_s) * (f->linesize[0]);
    entry_start_offset_c += ((uint32_t)einfo->ctb_y << (log2_ctb_s - 1)) * (f->linesize[1]);

    /*FIXME clean offset */
    fbuff->y  = (uint16_t *)&f->data[0][entry_start_offset];
    fbuff->cb = (uint16_t *)&f->data[1][entry_start_offset_c];
    fbuff->cr = (uint16_t *)&f->data[2][entry_start_offset_c];

    fbuff->stride   = f->linesize[0] >> 1;
    fbuff->stride_c = f->linesize[1] >> 1;
}

static void
fbuff_new_line(struct OVBuffInfo *fbuff, uint8_t log2_ctb_s)
{
    fbuff->y  += fbuff->stride << log2_ctb_s;
    fbuff->cb += fbuff->stride_c << (log2_ctb_s - 1);
    fbuff->cr += fbuff->stride_c << (log2_ctb_s - 1);
}

static int
slicedec_decode_rect_entry(OVSliceDec *sldec, const OVPS *const prms,
                           const struct RectEntryInfo *const einfo)
{
    int ctb_addr_rs = 0;
    int ctb_y = 0;
    int ret;

    /* FIXME handle more than one ctu dec */
    OVCTUDec *const ctudec = sldec->ctudec_list;
    /*FIXME handle cabac alloc or keep it on the stack ? */
    OVCABACCtx cabac_ctx;

    const int nb_ctu_w = einfo->nb_ctu_w;
    const int nb_ctu_h = einfo->nb_ctu_h;

    struct OVBuffInfo tmp_fbuff;
    ctudec->cabac_ctx = &cabac_ctx;

    ctudec->qp_ctx.current_qp = ctudec->slice_qp;
    derive_dequant_ctx(ctudec, &ctudec->qp_ctx, 0);

    /*FIXME quick tmvp import */
    ctudec->nb_ctb_pic_w = einfo->nb_ctb_pic_w;
    ctudec->drv_ctx.inter_ctx.tmvp_ctx.plane0 = &sldec->pic->mv_plane0;
    ctudec->drv_ctx.inter_ctx.tmvp_ctx.plane1 = &sldec->pic->mv_plane1;

    #if 0
    ctudec->drv_ctx.inter_ctx.tmvp_ctx.col_ref = sldec->pic->tmvp.collocated_ref;
    #endif
    #if 0
    if (sldec->pic->tmvp.collocated_ref) {
    #endif
    ctudec->drv_ctx.inter_ctx.tmvp_ctx.col_plane0 = &sldec->pic->tmvp.collocated_ref->mv_plane0;
    ctudec->drv_ctx.inter_ctx.tmvp_ctx.col_plane1 = &sldec->pic->tmvp.collocated_ref->mv_plane1;
    #if 0
    } else {
    ctudec->drv_ctx.inter_ctx.tmvp_ctx.col_plane0 = NULL;
    ctudec->drv_ctx.inter_ctx.tmvp_ctx.col_plane1 = NULL;
    }
    #endif

    ctudec->drv_ctx.inter_ctx.tmvp_ctx.ctudec = ctudec;
    ctudec->rcn_ctx.ctudec = ctudec;
    ctudec->drv_ctx.inter_ctx.tmvp_ctx.scale00 = sldec->pic->tmvp.scale00;
    ctudec->drv_ctx.inter_ctx.tmvp_ctx.scale10 = sldec->pic->tmvp.scale10;
    ctudec->drv_ctx.inter_ctx.tmvp_ctx.scale01 = sldec->pic->tmvp.scale01;
    ctudec->drv_ctx.inter_ctx.tmvp_ctx.scale11 = sldec->pic->tmvp.scale11;

    memset(ctudec->drv_ctx.inter_ctx.tmvp_ctx.dir_map_v0, 0, 33 * sizeof(uint64_t));
    memset(ctudec->drv_ctx.inter_ctx.tmvp_ctx.dir_map_v1, 0, 33 * sizeof(uint64_t));

    /* FIXME tmp Reset DBF */
    memset(&ctudec->dbf_info, 0, sizeof(ctudec->dbf_info));

    /* FIXME entry might be check before attaching entry to CABAC so there
     * is no need for this check
     */
    ret = ovcabac_attach_entry(ctudec->cabac_ctx, einfo->entry_start, einfo->entry_end);
    if (ret < 0) {
        return OVVC_EINDATA;
    }

    /* FIXME Note cabac context tables could be initialised earlier
     * so we could only init once and recopy context tables to others
     * entries CABAC readers
     */
    ovcabac_init_slice_context_table(cabac_ctx.ctx_table, prms->sh->sh_slice_type, ctudec->slice_qp);

    reset_cabac_lines(sldec, prms);

    reset_drv_lines(sldec, prms);

    cabac_line_next_line(ctudec, sldec);

    drv_line_next_line(ctudec, sldec);

    slicedec_attach_frame_buff(ctudec, sldec, einfo);

    tmp_fbuff = ctudec->rcn_ctx.frame_buff;

    while (ctb_y < nb_ctu_h - 1) {
        uint8_t log2_ctb_s = ctudec->part_ctx->log2_ctu_s;

        ctudec->ctb_y = einfo->ctb_y + ctb_y;

        /* New ctu line */
        ret = decode_ctu_line(ctudec, sldec, prms, einfo, ctb_addr_rs);

        cabac_line_next_line(ctudec, sldec);

        drv_line_next_line(ctudec, sldec);

        /*TODO 
         * CLeaner Next CTU Line
         */
        fbuff_new_line(&tmp_fbuff, log2_ctb_s);
        ctudec->rcn_ctx.frame_buff = tmp_fbuff;
        ctudec->ctu_ngh_flags = CTU_UP_FLG|CTU_UPRGT_FLG;

        ctb_addr_rs += nb_ctu_w;
        ctb_y++;
    }

    ctudec->ctb_y = einfo->ctb_y + ctb_y;
    /* Last line */
    if (!einfo->implicit_h) {
        ret = decode_ctu_line(ctudec, sldec, prms, einfo, ctb_addr_rs);
    } else {
        ret = decode_ctu_last_line(ctudec, sldec, prms, einfo, ctb_addr_rs);
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

/* FIXME clean this init */
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
    ctudec->enabled_mip = sps->sps_mip_enabled_flag;

    ctudec->jcbcr_enabled = sps->sps_joint_cbcr_enabled_flag;
    ctudec->enable_lfnst  = sps->sps_lfnst_enabled_flag;
    ctudec->isp_enabled   = sps->sps_isp_enabled_flag;
    ctudec->enable_mrl    = sps->sps_mrl_enabled_flag;

    ctudec->transform_skip_enabled = sps->sps_transform_skip_enabled_flag;
    ctudec->max_num_merge_candidates = 6 - sps->sps_six_minus_max_num_merge_cand;

    ctudec->delta_qp_enabled = pps->pps_cu_qp_delta_enabled_flag;

#if 1
    ctudec->dbf_disable = sh->sh_deblocking_filter_disabled_flag |
                          ph->ph_deblocking_filter_disabled_flag |
                          pps->pps_deblocking_filter_disabled_flag;
#else
    ctudec->dbf_disable = 1;
#endif

    ctudec->lm_chroma_enabled = sps->sps_cclm_enabled_flag;
    if (ctudec->lm_chroma_enabled) {
        /* FIXME add support vertical */
        if (sps->sps_chroma_vertical_collocated_flag /*sps->sps_chroma_horizontal_collocated_flag*/) {
            rcn_init_cclm_functions_collocated(&ctudec->rcn_ctx.rcn_funcs);
        } else {
            rcn_init_cclm_functions(&ctudec->rcn_ctx.rcn_funcs);
        }
    }

    rcn_init_ict_functions(&ctudec->rcn_ctx.rcn_funcs, ict_type(ph));

    slice_init_qp_ctx(ctudec, prms);

    derive_dequant_ctx(ctudec, &ctudec->qp_ctx, 0);

    init_coding_coeff_coding_ctx(ctudec, prms);

    init_part_info(ctudec, prms);

    init_slice_tree_ctx(ctudec, prms);

    ctudec->drv_ctx.inter_ctx.tmvp_enabled = ph->ph_temporal_mvp_enabled_flag;

    /*FIXME move rcn functions pointers init */
    rcn_init_mc_functions(&ctudec->rcn_ctx.rcn_funcs);

    /* Note it is important here that part info has already been set before calling
     * this function since it will be used to set line sizes*/

    /* FIXME
     * move this somewhere else so we can handle lines at a higher level
     *     safer is allocated check
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
uninit_ctudec_list(OVSliceDec *const sldec)
{
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
        drv_lines_uninit(sldec);
    }

    ov_freep(sldec_p);

}

/* FIXME cleanup */
static void
derive_ctu_neighborhood(const OVSliceDec *const sldec,
                       OVCTUDec *const ctudec,
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

/* FIXME refactor dequant and remove this function */
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
