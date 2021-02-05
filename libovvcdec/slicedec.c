#include "nvcl_structures.h"

#include "ctudec.h"



static int
init_slice_tree_ctx()
{
    if (sh->slice_type == VVC_SLICE_I) {
        lc_ctx->part_ctx        = &vvc_ctx->intra_part_ctx;
        lc_ctx->part_ctx_chroma = &vvc_ctx->intra_part_ctx_chroma;
        lc_ctx->active_part_map = &lc_ctx->part_map;
        if (sps->qtbtt_dual_tree_intra_flag) {
            lc_ctx->coding_tree          = &dual_tree;
            lc_ctx->coding_tree_implicit = &dual_tree_implicit;
        } else {
            lc_ctx->coding_tree          = &coding_quadtree;
            lc_ctx->coding_tree_implicit = &coding_quadtree_implicit;
            lc_ctx->coding_unit          = &coding_unit_intra_st;
            lc_ctx->transform_unit       = &transform_unit_st;
        }
    } else {
        lc_ctx->coding_tree          = &coding_quadtree;
        lc_ctx->coding_tree_implicit = &coding_quadtree_implicit;
        lc_ctx->coding_unit          = &coding_unit_inter_st;
        lc_ctx->prediction_unit      = sh->slice_type == VVC_SLICE_B ? &prediction_unit_inter_b : &prediction_unit_inter_p;
        lc_ctx->transform_unit = &transform_unit_st;
        /*FIXME use inter version*/
        lc_ctx->part_ctx        = &vvc_ctx->inter_part_ctx;
        lc_ctx->part_ctx_chroma = &vvc_ctx->intra_part_ctx_chroma;
        lc_ctx->active_part_map = &lc_ctx->part_map;
    }
}


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
    lc_ctx->alf_num_chroma_alt = vvc_ctx->alf_num_alt_chroma;

    lc_ctx->loop_filter_state_flags = flag;
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

static int
init_coding_coeff_coding_ctx()
{
    uint8_t ph_joint_cbcr_sign_flag;

    uint8_t sh_dep_quant_used_flag;
    uint8_t sh_sign_data_hiding_used_flag;

    uint8_t sh_ts_residual_coding_disabled_flag;

    uint8_t ph_chroma_residual_scale_flag;
    uint8_t ph_explicit_scaling_list_enabled_flag;
    uint8_t ph_scaling_list_aps_id;

    lc_ctx->mts_enabled  = sps->sps_mts_enabled_flag && sps->sps_explicit_mts_intra_enabled_flag;
    lc_ctx->mts_implicit = sps->sps_mts_enabled_flag && !sps->sps_explicit_mts_intra_enabled_flag;

    if (sh->slice_dep_quant_enabled_flag) {
        lc_ctx->residual_coding = &residual_coding_dpq;
        lc_ctx->residual_coding_chroma = &residual_coding_chroma_dpq;
        lc_ctx->residual_coding_isp_h = &residual_coding_isp_h_dpq;
        lc_ctx->residual_coding_isp_v = &residual_coding_isp_v_dpq;
    } else {
        lc_ctx->residual_coding_isp_h = &residual_coding_isp_h_sdh;
        lc_ctx->residual_coding_isp_v = &residual_coding_isp_v_sdh;
        lc_ctx->residual_coding_chroma = &residual_coding_chroma_sdh;
        lc_ctx->residual_coding = &residual_coding_sdh;
        lc_ctx->enable_sdh = sh->slice_sign_data_hiding_enabled_flag;
    }
}

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

int
init_slice_tools(CTUDec *ctu_dec, const OVSH *const sh)
{
    lc_ctx->max_log2_transform_skip_size = sps->log2_transform_skip_max_size_minus2 + 2;

    lc_ctx->alf_num_chroma_alt = vvc_ctx->alf_num_alt_chroma;
    lc_ctx->lm_chroma_enabled = sps->sps_cclm_enabled_flag;
    lc_ctx->enabled_mip = sps->sps_mip_enabled_flag;

    lc_ctx->jcbcr_enabled = sps->sps_joint_cbcr_enabled_flag;
    lc_ctx->enable_lfnst  = sps->sps_lfnst_enabled_flag;
    lc_ctx->isp_enabled  = sps->sps_isp_enabled_flag;
    lc_ctx->enable_mrl = sps->sps_mrl_enabled_flag;
    lc_ctx->transform_skip_enabled = sps->sps_transform_skip_enabled_flag;
    lc_ctx->max_num_merge_candidates = 6 - sps->six_minus_max_num_merge_cand;
    lc_ctx->delta_qp_enabled = pps->cu_qp_delta_enabled_flag;

    lc_ctx->dbf_disable = sh->slice_deblocking_filter_disabled_flag | ph->ph_deblocking_filter_disabled_flag | pps->pps_deblocking_filter_disabled_flag;

    init_qp_ctx(lc_ctx, vvc_ctx);

    derive_dequant_ctx(lc_ctx, &lc_ctx->qp_ctx, 0);

    return 0;
}

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

static void
derive_ctu_neighborhood(const VVCContext *const vvc_ctx,
                       VVCLocalContext *const lc_ctx,
                       int ctb_address, int nb_ctu_w, int nb_ctu_h)
{
    int is_left_border = ((ctb_address) % nb_ctu_w)  ? 0 : 1;
    int is_up_border   = (ctb_address < nb_ctu_w) ? 1 : 0;
    int is_right_border= ((ctb_address + 1) % nb_ctu_w);

    uint8_t ctb_flags = 0;

    /* FIXME clean */
    if (!is_left_border) {
        flags |= VVC_CTU_LEFT_FLAG;
    }
    if (!is_up_border) {
        flags |= VVC_CTU_UP_FLAG;
    }
    if (!is_left_border && !is_up_border) {
        flags |= VVC_CTU_UPLEFT_FLAG;
    }
    if (!is_up_border && is_right_border) {
        flags |= VVC_CTU_UPRIGHT_FLAG;
    }

    lc_ctx->ctu_neighbour_flags = flags;
}


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
            lc_ctx->qp_ctx.current_qp = line_ctx.qp_x_map[0];

        /* Derive border ctx */
        load_ctu_ctx(lc_ctx, &line_map, &line_ctx, ctb_x, vvc_ctx->sh.slice_type);

        ovdec_decode_ctu(vvc_ctx, lc_ctx, ctb_x, ctb_y, 0);

        /* Post ctu reconstruction */
        if (!lc_ctx->dbf_disable)
            apply_dbf_ctu(lc_ctx, ctb_x, ctb_y, nb_ctu_w, nb_ctu_h);

        store_ctu_ctx(lc_ctx, &line_map, &line_ctx, ctb_x, vvc_ctx->sh.slice_type);

        if (vvc_ctx->ps.sps_data->sps_temporal_mvp_enabled_flag) {
            store_tmvp(vvc_ctx, lc_ctx, &lc_ctx->inter_ctx.tmvp_ctx);
        }
        ctb_addr_rs++;

        update_fd(&lc_ctx->frame_data, lc_ctx->part_ctx->log2_ctu_s);
        if (!(ctb_addr_rs % nb_ctu_w)) {
            update_fd_line(&lc_ctx->frame_data, nb_ctu_w, lc_ctx->part_ctx->log2_ctu_s);
        }
    }
}



static int
derive_dequant_ctx(VVCLocalContext *const lc_ctx, const VVCQPCTX *const qp_ctx,
                  int cu_qp_delta)
{
    /*FIXME avoid negative values especiallly in chroma_map derivation*/
    int base_qp = (qp_ctx->current_qp + cu_qp_delta + 64) & 63;
    lc_ctx->dequant_luma.qp = ((base_qp + 12) & 63);
    /*FIXME update transform skip ctx to VTM-10.0 */
    lc_ctx->dequant_luma_skip.qp = FFMAX(lc_ctx->dequant_luma.qp, qp_ctx->min_qp_prime_ts);
    lc_ctx->dequant_cb.qp = qp_ctx->chroma_qp_map_cb[(base_qp + qp_ctx->cb_offset + 64) & 63] + 12;
    lc_ctx->dequant_cr.qp = qp_ctx->chroma_qp_map_cr[(base_qp + qp_ctx->cr_offset + 64) & 63] + 12;
    lc_ctx->dequant_joint_cb_cr.qp = qp_ctx->chroma_qp_map_jcbcr[(base_qp + qp_ctx->jcbcr_offset + 64) & 63] + 12;
    lc_ctx->qp_ctx.current_qp = base_qp;
}

static void av_always_inline
load_depth_maps(const struct LineCtxs *const l, VVCLocalContext *const lc_ctx, int ctb_x)
{

    const uint8_t flags = lc_ctx->ctu_neighbour_flags;
    struct PartMap *part_map   = &lc_ctx->part_map;
    struct PartMap *part_map_c = &lc_ctx->part_map_c;
    int log2_nb_pb_ctu = 5 - (7 - lc_ctx->part_ctx->log2_ctu_s);
    int x_cb = ctb_x << log2_nb_pb_ctu;

    part_map->qt_depth_map_x  = &l->qt_depth_map_x_l [x_cb];
    part_map->log2_cu_w_map_x = &l->log2_cu_w_map_x_l[x_cb];

    part_map_c->qt_depth_map_x  = &l->qt_depth_map_x_c [x_cb];
    part_map_c->log2_cu_w_map_x = &l->log2_cu_w_map_x_c[x_cb];

    if (!(flags&VVC_CTU_LEFT_FLAG)) {
        memset(part_map->qt_depth_map_y    , VVC_NA_DEPTH, sizeof(uint8_t) << log2_nb_pb_ctu);
        memset(part_map->log2_cu_h_map_y   , VVC_NA_DEPTH, sizeof(uint8_t) << log2_nb_pb_ctu);

        memset(part_map_c->qt_depth_map_y  , VVC_NA_DEPTH, sizeof(uint8_t) << log2_nb_pb_ctu);
        memset(part_map_c->log2_cu_h_map_y , VVC_NA_DEPTH, sizeof(uint8_t) << log2_nb_pb_ctu);
    }

    if (!(flags&VVC_CTU_UP_FLAG)) {
        memset(part_map->qt_depth_map_x    , VVC_NA_DEPTH, sizeof(uint8_t) << log2_nb_pb_ctu);
        memset(part_map->log2_cu_w_map_x   , VVC_NA_DEPTH, sizeof(uint8_t) << log2_nb_pb_ctu);

        memset(part_map_c->qt_depth_map_x  , VVC_NA_DEPTH, sizeof(uint8_t) << log2_nb_pb_ctu);
        memset(part_map_c->log2_cu_w_map_x , VVC_NA_DEPTH, sizeof(uint8_t) << log2_nb_pb_ctu);
    }
}


static void
fill_ctu_borders(VVCLocalContext *const lc_ctx, const AVFrame *const frame,
                 uint16_t lim_w, uint16_t lim_h)
{
    const uint8_t log2_ctb_s = lc_ctx->part_ctx->log2_ctu_s;
    const struct VVCFrameData *const fd = &lc_ctx->frame_data;
    /*FIXME no need to copy from frame we can use line buffers and previous ctx*/

    const int frame_stride   = frame->linesize[0] >> 1;
    const int frame_stride_c = frame->linesize[1] >> 1;

    const uint8_t ctb_ngh_flg = lc_ctx->ctu_neighbour_flags;

    if ((ctb_ngh_flg & VVC_CTU_LEFT_FLAG)) {
        const uint16_t *src_y  = (uint16_t *) (fd->data_y)  - 4;
        const uint16_t *src_cb = (uint16_t *) (fd->data_cb) - 1;
        const uint16_t *src_cr = (uint16_t *) (fd->data_cr) - 1;

        uint16_t *dst_y  = &lc_ctx->ctu_data_y [VVC_CTB_OFFSET - 4];
        uint16_t *dst_cb = &lc_ctx->ctu_data_cb[VVC_CTB_OFFSET_CHROMA - 1];
        uint16_t *dst_cr = &lc_ctx->ctu_data_cr[VVC_CTB_OFFSET_CHROMA - 1];

        for (int i = 0; i < lim_h; ++i) {
            memcpy(dst_y, src_y, sizeof(uint16_t) * 4);
            dst_y += VVC_CTB_STRIDE;
            src_y += frame_stride;
        }

        for (int i = 0; i < (lim_h >> 1); i++) {
            memcpy(dst_cb, src_cb, sizeof(uint16_t));
            memcpy(dst_cr, src_cr, sizeof(uint16_t));
            dst_cb += VVC_CTB_STRIDE_CHROMA;
            dst_cr += VVC_CTB_STRIDE_CHROMA;
            src_cb += frame_stride_c;
            src_cr += frame_stride_c;
        }
    }

    if ((ctb_ngh_flg & VVC_CTU_UP_FLAG)) {
        const uint16_t *src_y  = (uint16_t *) (fd->data_y)  - frame_stride;
        const uint16_t *src_cb = (uint16_t *) (fd->data_cb) - frame_stride_c;
        const uint16_t *src_cr = (uint16_t *) (fd->data_cr) - frame_stride_c;

        uint16_t *dst_y  = &lc_ctx->ctu_data_y [VVC_CTB_OFFSET - VVC_CTB_STRIDE];
        uint16_t *dst_cb = &lc_ctx->ctu_data_cb[VVC_CTB_OFFSET_CHROMA - VVC_CTB_STRIDE_CHROMA];
        uint16_t *dst_cr = &lc_ctx->ctu_data_cr[VVC_CTB_OFFSET_CHROMA - VVC_CTB_STRIDE_CHROMA];
        if (ctb_ngh_flg & VVC_CTU_UPRIGHT_FLAG)
           lim_w <<= 1;

        if (lim_w > 196) lim_w = 196;
        memcpy(dst_y,  src_y,  sizeof(uint16_t) * lim_w);
        memcpy(dst_cb, src_cb, sizeof(uint16_t) * (lim_w >> 1));
        memcpy(dst_cr, src_cr, sizeof(uint16_t) * (lim_w >> 1));
    }

    if ((ctb_ngh_flg & (VVC_CTU_UP_FLAG | VVC_CTU_LEFT_FLAG)) == (VVC_CTU_UP_FLAG | VVC_CTU_LEFT_FLAG)) {
        const uint16_t *src_y  = (uint16_t *) (fd->data_y)  - frame_stride   - 4;
        const uint16_t *src_cb = (uint16_t *) (fd->data_cb) - frame_stride_c - 1;
        const uint16_t *src_cr = (uint16_t *) (fd->data_cr) - frame_stride_c - 1;

        uint16_t *dst_y  = &lc_ctx->ctu_data_y[VVC_CTB_OFFSET - VVC_CTB_STRIDE - 4];
        uint16_t *dst_cb = &lc_ctx->ctu_data_cb[VVC_CTB_OFFSET_CHROMA - VVC_CTB_STRIDE_CHROMA - 1];
        uint16_t *dst_cr = &lc_ctx->ctu_data_cr[VVC_CTB_OFFSET_CHROMA - VVC_CTB_STRIDE_CHROMA - 1];

        memcpy(dst_y, src_y, sizeof(uint16_t) * 4);
        *dst_cb = *src_cb;
        *dst_cr = *src_cr;
    }
}

static void
fill_ctu(VVCLocalContext *const lc_ctx, const AVFrame *const frame)
{
    const uint8_t log2_ctb_s = lc_ctx->part_ctx->log2_ctu_s;
    const struct VVCFrameData *const fd = &lc_ctx->frame_data;
    /*FIXME no need to copy from frame we can use line buffers and previous ctx*/

    const int frame_stride   = frame->linesize[0] >> 1;
    const int frame_stride_c = frame->linesize[1] >> 1;

    const uint8_t ctb_ngh_flg = lc_ctx->ctu_neighbour_flags;

    if ((ctb_ngh_flg & VVC_CTU_LEFT_FLAG)) {
        const int lim_h = 1 << log2_ctb_s;

        const uint16_t *src_y  = (uint16_t *) (fd->data_y)  - 4;
        const uint16_t *src_cb = (uint16_t *) (fd->data_cb) - 1;
        const uint16_t *src_cr = (uint16_t *) (fd->data_cr) - 1;

        uint16_t *dst_y  = &lc_ctx->ctu_data_y [VVC_CTB_OFFSET - 4];
        uint16_t *dst_cb = &lc_ctx->ctu_data_cb[VVC_CTB_OFFSET_CHROMA - 1];
        uint16_t *dst_cr = &lc_ctx->ctu_data_cr[VVC_CTB_OFFSET_CHROMA - 1];

        for (int i = 0; i < lim_h; ++i) {
            memcpy(dst_y, src_y, sizeof(uint16_t) * 4);
            dst_y += VVC_CTB_STRIDE;
            src_y += frame_stride;
        }

        for (int i = 0; i < (lim_h >> 1); i++) {
            memcpy(dst_cb, src_cb, sizeof(uint16_t));
            memcpy(dst_cr, src_cr, sizeof(uint16_t));
            dst_cb += VVC_CTB_STRIDE_CHROMA;
            dst_cr += VVC_CTB_STRIDE_CHROMA;
            src_cb += frame_stride_c;
            src_cr += frame_stride_c;
        }
    }

    if ((ctb_ngh_flg & VVC_CTU_UP_FLAG)) {
        int lim_w = 1 << (log2_ctb_s + 1);
        if (lim_w > 196) lim_w = 196;
        const uint16_t *src_y  = (uint16_t *) (fd->data_y)  - frame_stride;
        const uint16_t *src_cb = (uint16_t *) (fd->data_cb) - frame_stride_c;
        const uint16_t *src_cr = (uint16_t *) (fd->data_cr) - frame_stride_c;

        uint16_t *dst_y  = &lc_ctx->ctu_data_y [VVC_CTB_OFFSET - VVC_CTB_STRIDE];
        uint16_t *dst_cb = &lc_ctx->ctu_data_cb[VVC_CTB_OFFSET_CHROMA - VVC_CTB_STRIDE_CHROMA];
        uint16_t *dst_cr = &lc_ctx->ctu_data_cr[VVC_CTB_OFFSET_CHROMA - VVC_CTB_STRIDE_CHROMA];

        memcpy(dst_y,  src_y,  sizeof(uint16_t) * lim_w);
        memcpy(dst_cb, src_cb, sizeof(uint16_t) * (lim_w >> 1));
        memcpy(dst_cr, src_cr, sizeof(uint16_t) * (lim_w >> 1));
    }

    if ((ctb_ngh_flg & (VVC_CTU_UP_FLAG | VVC_CTU_LEFT_FLAG)) == (VVC_CTU_UP_FLAG | VVC_CTU_LEFT_FLAG)) {
        const uint16_t *src_y  = (uint16_t *) (fd->data_y)  - frame_stride   - 4;
        const uint16_t *src_cb = (uint16_t *) (fd->data_cb) - frame_stride_c - 1;
        const uint16_t *src_cr = (uint16_t *) (fd->data_cr) - frame_stride_c - 1;

        uint16_t *dst_y  = &lc_ctx->ctu_data_y[VVC_CTB_OFFSET - VVC_CTB_STRIDE - 4];
        uint16_t *dst_cb = &lc_ctx->ctu_data_cb[VVC_CTB_OFFSET_CHROMA - VVC_CTB_STRIDE_CHROMA - 1];
        uint16_t *dst_cr = &lc_ctx->ctu_data_cr[VVC_CTB_OFFSET_CHROMA - VVC_CTB_STRIDE_CHROMA - 1];

        memcpy(dst_y, src_y, sizeof(uint16_t) * 4);
        *dst_cb = *src_cb;
        *dst_cr = *src_cr;
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

    #if 0
    line_ctx->intra_modes_x_map += offset << 5;
    line_ctx->inter_modes_x_map += offset << 5;

    line_ctx->mv_row0 += offset << 5;
    line_ctx->mv_row1 += offset << 5;

    line_ctx->qp_x_map += offset << 5;
    line_ctx->qp_x_map_cb += offset << 5;
    line_ctx->qp_x_map_cr += offset << 5;

    line_ctx->dbf_qp_ver    += offset << 5;
    line_ctx->dbf_qp_ver_cb += offset << 5;
    line_ctx->dbf_qp_ver_cr += offset << 5;

    line_map->progress_map_l += offset;
    line_map->progress_map_c += offset;

    line_map->inter_map_rows_0 += offset;
    line_map->inter_map_rows_1 += offset;

    line_map->dbf_edge_ver += offset;
    line_map->dbf_edge_hor += offset;

    line_map->dbf_bs1_ver += offset;
    line_map->dbf_bs1_hor += offset;

    line_map->large_map_c += offset;

    line_map->dbf_bs1_ver_cb += offset;
    line_map->dbf_bs1_hor_cb += offset;

    line_map->dbf_bs1_ver_cr += offset;
    line_map->dbf_bs1_hor_cr += offset;

    line_map->dbf_bs2_ver += offset;
    line_map->dbf_bs2_hor += offset;
    #endif
}

static void
init_frame_data(const AVFrame *const frame, struct VVCFrameData *const fd,
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

static void
load_ctu_ctx(VVCLocalContext *const lc_ctx,
             struct LineMaps *const line_map,
             struct LineCtxs *const line_ctx,
             int ctb_x, int slice_type)
{
    fill_cu_availibility_map(lc_ctx, line_map, ctb_x);

    if (slice_type != VVC_SLICE_I) {
        struct VVCInterCtx *const inter_ctx = &lc_ctx->inter_ctx;
        struct VVCMVCtx *const mv_ctx0 = &inter_ctx->mv_ctx0;
        const struct LineMaps *const m = line_map;
        fill_inter_map(lc_ctx, mv_ctx0->map.rows, mv_ctx0->map.cols,
                (uint64_t)m->inter_map_rows_0[ctb_x],
                (uint64_t)m->inter_map_rows_0[ctb_x + 1],
                ctb_x);
    }

    if (slice_type == VVC_SLICE_B) {
        struct VVCInterCtx *const inter_ctx = &lc_ctx->inter_ctx;
        struct VVCMVCtx *const mv_ctx1 = &inter_ctx->mv_ctx1;
        const struct LineMaps *const m = line_map;
        fill_inter_map(lc_ctx, mv_ctx1->map.rows, mv_ctx1->map.cols,
                (uint64_t)m->inter_map_rows_1[ctb_x],
                (uint64_t)m->inter_map_rows_1[ctb_x + 1],
                ctb_x);
    }

    if (slice_type != VVC_SLICE_I)
        load_inter_ngh_ctx(lc_ctx, line_ctx, ctb_x);

    load_depth_maps(line_ctx, lc_ctx, ctb_x);

    load_intra_ngh_ctx(lc_ctx, line_ctx, ctb_x);

    if (!lc_ctx->dbf_disable) {
        load_dbf_ctx(&lc_ctx->dbf_info, line_ctx, lc_ctx, ctb_x);

        load_edge_map(&lc_ctx->dbf_info, line_map, lc_ctx->part_ctx, ctb_x);

        load_bs_map(&lc_ctx->dbf_info, line_map, lc_ctx->part_ctx, ctb_x);
    }
}

static void
store_ctu_ctx(VVCLocalContext *const lc_ctx,
              struct LineMaps *const line_map,
              struct LineCtxs *const line_ctx,
              int ctb_x, int slice_type)
{
    store_availability_maps(line_ctx, line_map, lc_ctx, ctb_x);

    if (slice_type != VVC_SLICE_I) {
        store_inter_maps(line_ctx, line_map, lc_ctx, ctb_x);
    }

    if (!lc_ctx->dbf_disable) {
        store_edge_map(&lc_ctx->dbf_info, line_map, lc_ctx->part_ctx, ctb_x);
        store_bs_map(&lc_ctx->dbf_info, line_map, lc_ctx->part_ctx, ctb_x);
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
    #if 0
    ov_freep(&l->intra_modes_x_map);
    ov_freep(&l->inter_modes_x_map);
    ov_freep(&l->mv_row0);
    ov_freep(&l->mv_row1);
    ov_freep(&l->qp_x_map);
    ov_freep(&l->qp_x_map_cb);
    ov_freep(&l->qp_x_map_cr);
    ov_freep(&l->dbf_qp_ver);
    ov_freep(&l->dbf_qp_ver_cb);
    ov_freep(&l->dbf_qp_ver_cr);
    ov_freep(&m->progress_map_l);
    ov_freep(&m->progress_map_c);
    ov_freep(&m->inter_map_rows_0);
    ov_freep(&m->inter_map_rows_1);
    ov_freep(&m->dbf_edge_ver);
    ov_freep(&m->dbf_edge_hor);
    ov_freep(&m->dbf_bs1_ver);
    ov_freep(&m->dbf_bs1_hor);
    ov_freep(&m->large_map_c);
    ov_freep(&m->dbf_bs1_ver_cb);
    ov_freep(&m->dbf_bs1_hor_cb);
    ov_freep(&m->dbf_bs1_ver_cr);
    ov_freep(&m->dbf_bs1_hor_cr);
    ov_freep(&m->dbf_bs2_ver);
    ov_freep(&m->dbf_bs2_hor);

    ov_freep(&s->ctb_alf_flag);
    #endif
}

