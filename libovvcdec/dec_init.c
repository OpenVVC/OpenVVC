#include "nvcl_structures.h"
#include "dec_structures.h"

/* This file is intended to translate Parameter Sets information to 
 * structures more usable by the decoder contexts
 * It contains translations from nvcl structures to decoder structures
 * and Parameters sets activation functions.
 */

static int
sps_init_partition_constraint_info(OVPartInfo *const pinfo, const OVSPS *const sps)
{
    uint8_t log2_ctu_s = sps->sps_log2_ctu_size_minus5 + 5;
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

#if 0
static int update_partition_context_inter(VVCPartSize *const part_ctx,
                                    const VVCSPSData *const sps){


    part_ctx->log2_ctu_s    = sps->sps_log2_ctu_size_minus5 + 5;
    part_ctx->log2_min_cb_s = sps->log2_min_luma_coding_block_size_minus2 + 2;

    part_ctx->log2_min_qt_s = sps->sps_log2_diff_min_qt_min_cb_inter_slice+
            part_ctx->log2_min_cb_s;

    part_ctx->log2_max_bt_s = part_ctx->log2_min_qt_s +
            sps->sps_log2_diff_max_bt_min_qt_inter_slice;
    part_ctx->log2_max_tt_s = part_ctx->log2_min_qt_s +
            sps->sps_log2_diff_max_tt_min_qt_inter_slice;

    part_ctx->max_mtt_depth = sps->sps_max_mtt_hierarchy_depth_inter_slice;

    part_ctx->log2_max_tb_s = 5 + sps->sps_max_luma_transform_size_64_flag;

    return 0;
}

static int update_partition_context_intra_chroma(VVCPartSize *const part_ctx,
                                                 const VVCSPSData *const sps)
{
    part_ctx->log2_ctu_s    = sps->sps_log2_ctu_size_minus5 + 5;
    part_ctx->log2_min_cb_s = sps->log2_min_luma_coding_block_size_minus2 + 2 - 1;

    part_ctx->log2_max_tb_s = 5 + sps->sps_max_luma_transform_size_64_flag - 1;

    part_ctx->log2_min_qt_s = sps->sps_log2_diff_min_qt_min_cb_intra_slice_chroma +
            part_ctx->log2_min_cb_s;

    part_ctx->log2_max_bt_s = part_ctx->log2_min_qt_s +
            sps->sps_log2_diff_max_bt_min_qt_intra_slice_chroma;
    part_ctx->log2_max_tt_s = part_ctx->log2_min_qt_s +
            sps->sps_log2_diff_max_tt_min_qt_intra_slice_chroma;

    part_ctx->max_mtt_depth = sps->sps_max_mtt_hierarchy_depth_intra_slice_chroma;

    return 0;
}
#endif

static void
sps_init_chroma_qp_tables(uint8_t *const qp_map_table, const OVSPS *const sps)
{
    uint8_t nb_qp_points = sps->sps_num_points_in_qp_table_minus1;
    int idx = sps->sps_qp_table_starts_minus26 + 26;
    int next_idx;
    int j;

    /*FIXME:
     *    -actual qp can be clipped to -qp_bitdepth_offset
     *    -fill first part of qp_map + check corner cases
     */

    qp_map_table[idx] = idx;
    for (j = idx - 1; j >= 0; --j){
        qp_map_table[j] = qp_map_table[j + 1] - 1;
    }

    for (j = 0; j <= nb_qp_points; ++j) {
        const int nb_step = sps->delta_qp_in_val_minus1[j] + 1;
        const int nb_qp_w = sps->delta_qp_in_val_minus1[j] ^ sps->delta_qp_diff_val[j];
        const int first_qp = qp_map_table[idx];
        const int round = nb_step >> 1;
        int qp_sum = nb_qp_w;
        int k;

        next_idx = idx + nb_step;

        for (k = idx + 1;  k <= next_idx; ++k) {
            qp_map_table[k]  = first_qp + (qp_sum + round) / nb_step;
            qp_sum += nb_qp_w;
        }
        idx = next_idx;
    }

    for (j = next_idx;  j < 64; ++j) {
        qp_map_table[j] = qp_map_table[j - 1] + 1;
    }
}

static void
sps_init_dpb_parameters(struct OVPS *prms)
{
    
}

static void
pps_init_qp(struct OVPS *prms)
{
    int base_qp = pps->pps_init_qp_minus26 + 26;
    /* FIXME chroma_qp_offsets */


  

}

static void
slice_init_qp_ctx(OVCTUDec *const ctudec, const struct OVPS *const prms)
{
    const uint8_t qp_bd_offset = 12;
    uint8_t pic_base_qp = pps->pps_init_qp_minus26 + 26;
    uint8_t pic_qp = pic_base_qp + ph->ph_qp_delta;
    int8_t slice_qp = pic_qp + sh->qp_delta;
    int8_t cb_qp_offset = sh->sh_cb_qp_offset + pps->pps_cb_qp_offset;
    int8_t cr_qp_offset = sh->sh_cr_qp_offset + pps->pps_cr_qp_offset;
    int8_t jcbcr_qp_offset = sh->sh_slice_joint_cbcr_qp_offset + pps->pps_joint_cbcr_qp_offset_value;

    ctudec->slice_qp = pic_qp + sh->sh->qp_delta;
    ctudec->dequant_luma.qp = slice_qp + qp_bd_offset;

    qp_ctx->chroma_qp_map_cb    = sps_info->chroma_qp_mapping_tables;
    qp_ctx->chroma_qp_map_cr    = vvc_ctx->chroma_qp_mapping_tables;
    qp_ctx->chroma_qp_map_jcbcr = vvc_ctx->chroma_qp_mapping_tables;

    qp_ctx->current_qp = slice_qp;
    qp_ctx->cb_offset = cb_qp_offset;
    qp_ctx->cr_offset = cr_qp_offset;
    qp_ctx->jcbcr_offset = jcbcr_qp_offset;

    lc_ctx->dequant_cb.qp = qp_ctx->chroma_qp_map_cb[slice_qp + cb_qp_offset] + qp_bd_offset;
    lc_ctx->dequant_cr.qp = qp_ctx->chroma_qp_map_cr[slice_qp + cr_qp_offset] + qp_bd_offset;
    lc_ctx->dequant_joint_cb_cr.qp = qp_ctx->chroma_qp_map_jcbcr[base_qp + jcbcr_qp_offset] + qp_bd_offset;
}



#if 0
static void
init_qp_ctx(VVCLocalContext *const lc_ctx, const VVCContext *const vvc_ctx)
{
    const VVCPPSData *const pps = vvc_ctx->ps.pps_data;
    const VVCSPSData *const sps = vvc_ctx->ps.sps_data;
    const VVCPictureHeaderData *const ph = vvc_ctx->ps.ph_data;
    const VVCSliceHeaderData *const sh = &vvc_ctx->sh;
		VVCQPCTX *const qp_ctx = &lc_ctx->qp_ctx;
    int base_qp = pps->init_qp_minus26 + 26;
    int cb_qp_offset = sh->slice_cb_qp_offset + pps->pps_cb_qp_offset;
    int cr_qp_offset = sh->slice_cr_qp_offset + pps->pps_cr_qp_offset;
    int jcbcr_qp_offset = sh->slice_joint_cbcr_qp_offset + pps->pps_joint_cbcr_qp_offset_value;

    base_qp += sh->slice_qp_delta;
    base_qp += ph->ph_qp_delta;


    lc_ctx->slice_qp = base_qp;

    lc_ctx->dequant_luma.qp = base_qp + 12;

    /*FIXME handle per components chroma_tables */
    qp_ctx->chroma_qp_map_cb =    vvc_ctx->chroma_qp_mapping_tables;
    qp_ctx->chroma_qp_map_cr =    vvc_ctx->chroma_qp_mapping_tables;
    qp_ctx->chroma_qp_map_jcbcr = vvc_ctx->chroma_qp_mapping_tables;

    qp_ctx->current_qp = base_qp;
    qp_ctx->cb_offset = cb_qp_offset;
    qp_ctx->cr_offset = cr_qp_offset;
    qp_ctx->jcbcr_offset = jcbcr_qp_offset;

    lc_ctx->dequant_cb.qp = qp_ctx->chroma_qp_map_cb[base_qp + cb_qp_offset] + 12;
    lc_ctx->dequant_cr.qp = qp_ctx->chroma_qp_map_cr[base_qp + cr_qp_offset] + 12;
    lc_ctx->dequant_joint_cb_cr.qp = qp_ctx->chroma_qp_map_jcbcr[base_qp + jcbcr_qp_offset] + 12;

}
#endif
