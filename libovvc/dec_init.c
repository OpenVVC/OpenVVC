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


static int
sps_init_partition_constraint_info_intra(OVPartInfo *const pinfo, const OVSPS *const sps)
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
    uint8_t nb_qp_points = nb_qp_min1;
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
        dst_qp_tab[j] = dst_qp_tab[j - 1] + 1;
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
        uint8_t *dst_qp_tab = sps_info->qp_tables_c[i].qp;

        sps_fill_qp_table(dst_qp_tab, dqp_input_val_min1, dqp_diff_val, start_qp,
                          nb_qp_min1);
    }

    /* Copy first table into others in case of same qp for all chroma tables */
    for (; i < 3 - sps->sps_joint_cbcr_enabled_flag; ++i) {
        sps_info->qp_tables_c[i] = sps_info->qp_tables_c[i-1];
    }
}

static void
sps_init_dpb_parameters(OVVCDec *const dec, struct OVPS *const prms)
{
    const OVSPS *const sps = prms->sps;
    const OVDPBParams *dpb_list = sps->dpb_parameters;

    #if 0
    dpb_init_params(dpb, &dpb_params[0]);
    #endif

}

#if 0
static void
pps_init_qp(struct OVPS *prms)
{
    int base_qp = pps->pps_init_qp_minus26 + 26;
    /* FIXME chroma_qp_offsets */




}

static int
pps_init_dpb(const OVPPS *const pps)
{
    pps_pic_width_in_luma_samples;
    pps_pic_height_in_luma_samples;
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

    qp_ctx->chroma_qp_map_cb    = sps_info->chroma_qp_mapping_tables;
    qp_ctx->chroma_qp_map_cr    = sps_info->chroma_qp_mapping_tables;
    qp_ctx->chroma_qp_map_jcbcr = sps_info->chroma_qp_mapping_tables;

    qp_ctx->current_qp = slice_qp;
    qp_ctx->cb_offset = cb_qp_offset;
    qp_ctx->cr_offset = cr_qp_offset;
    qp_ctx->jcbcr_offset = jcbcr_qp_offset;

    ctudec->dequant_luma.qp = slice_qp + qp_bd_offset;
    ctudec->dequant_cb.qp = qp_ctx->chroma_qp_map_cb[slice_qp + cb_qp_offset] + qp_bd_offset;
    ctudec->dequant_cr.qp = qp_ctx->chroma_qp_map_cr[slice_qp + cr_qp_offset] + qp_bd_offset;
    ctudec->dequant_joint_cb_cr.qp = qp_ctx->chroma_qp_map_jcbcr[base_qp + jcbcr_qp_offset] + qp_bd_offset;
}
#endif

static int
update_sps_info(struct SPSInfo *const sps_info, const OVSPS *const sps)
{
    sps_init_partition_constraint_info_intra(&sps_info->part_info[0], sps);
    sps_init_partition_constraint_info_inter(&sps_info->part_info[1], sps);
    sps_init_partition_constraint_info_chroma(&sps_info->part_info[2], sps);

    sps_init_chroma_qp_tables(sps_info, sps);

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
        sh_info->rbsp_entry[i] = nal->rbsp_data + rbsp_offset[i] + nb_sh_bytes;
    }

    /* Note this is so we can retrieve entry end by using rbsp_entry [i + 1] */
    sh_info->rbsp_entry[i] = nal->rbsp_data + nal->rbsp_size;

    /*FIXME check entries do not exceed rpbs size */
    return 0;
}

static void
init_tile_ctx(struct TileInfo *const tinfo, const OVPPS *const pps)
{
    int nb_cols = pps->pps_num_exp_tile_columns_minus1 + 1;
    int nb_rows = pps->pps_num_exp_tile_rows_minus1 + 1;

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

    /* FIXME review implicit last tile x and y
    */
    for (i = 1; i <=  nb_rows; ++i) {
        const int tile_nb_ctu_h = pps->pps_tile_row_height_minus1[i - 1] + 1;
        tinfo->ctu_y[i - 1] = nb_ctu_h - rem_ctu_h;
        tinfo->nb_ctu_h[i - 1] = tile_nb_ctu_h;
        rem_ctu_h -= tile_nb_ctu_h;
    }
    tinfo->ctu_y[i - 1] = nb_ctu_h - rem_ctu_h;
    tinfo->nb_ctu_h[i - 1] = rem_ctu_h;

    for (i = 1; i <= nb_cols; ++i) {
        const int tile_nb_ctu_w = pps->pps_tile_column_width_minus1[i - 1] + 1;
        tinfo->ctu_x[i - 1] = nb_ctu_w - rem_ctu_w;
        tinfo->nb_ctu_w[i - 1] = tile_nb_ctu_w;
        rem_ctu_w -= tile_nb_ctu_w;
    }

    tinfo->ctu_x[i - 1] = nb_ctu_w - rem_ctu_w;
    tinfo->nb_ctu_w[i - 1] = rem_ctu_w;

    tinfo->nb_tile_cols = nb_cols + !!rem_ctu_w;
    tinfo->nb_tile_rows = nb_rows + !!rem_ctu_h;
}

static int
update_pps_info(struct PPSInfo *const pps_info, const OVPPS *const pps,
                const OVSPS *const sps)
{
    struct TileInfo *const tinfo = &pps_info->tile_info;
    uint8_t have_tile = pps->pps_num_exp_tile_columns_minus1 + pps->pps_num_exp_tile_rows_minus1;

    if (have_tile) {
        init_tile_ctx(tinfo, pps);
    } else {
        /* Initialize tile context to 1 tile on whole picture based on
         * SPS informations
         */
        const int log2_ctu_s = sps->sps_log2_ctu_size_minus5 + 5;

        /* FIXME use SPS max pic size or pps ?  */
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

    return 0;
}

static int
update_ph_info(struct PHInfo *const ph_info, const OVPH *const ph)
{
    return 0;
}

static int
update_sh_info(struct SHInfo *const sh_info, const OVSH *const sh)
{
    return 0;
}

static const OVSPS *
retreive_sps(const OVNVCLCtx *const nvcl_ctx, const OVPPS *const pps)
{
    uint8_t sps_id = pps->pps_seq_parameter_set_id;
    const OVSPS *sps = NULL;
    if (sps_id < 16) {
        sps = nvcl_ctx->sps_list[sps_id];
    } else {
        ov_log(NULL, 3, "Invalid SPS ID  %d\n", sps_id);
    }

    return sps;
}

static const OVPPS *
retreive_pps(const OVNVCLCtx *const nvcl_ctx, const OVPH *const ph)
{
    uint8_t pps_id = ph->ph_pic_parameter_set_id;
    const OVPPS *pps = NULL;
    if (pps_id < 16) {
        pps = nvcl_ctx->pps_list[pps_id];
    } else {
        ov_log(NULL, 3, "Invalid PPS ID  %d\n", pps_id);
    }

    return pps;
}

int
decinit_update_params(OVVCDec *const dec, const OVNVCLCtx *const nvcl_ctx)
{
    /* FIXME assert nvcl_ctx params sets are not NULL*/
    struct OVPS *const ps = &dec->active_params;
    int ret;
    const OVSH *const sh = nvcl_ctx->sh;
    const OVPH *const ph = nvcl_ctx->ph;
    const OVPPS *const pps = retreive_pps(nvcl_ctx, ph);
    const OVSPS *const sps = retreive_sps(nvcl_ctx, pps);

    if (!sh || !ph || !pps || !sps) {
        ov_log(NULL, 3, "Missing Parameter sets for dec initialisation\n");
        return OVVC_EINDATA;
    }

    if (ps->sps != sps) {
        ret = update_sps_info(&ps->sps_info, sps);
        if (ret < 0) {
            goto failsps;
        }

        ps->sps = sps;
    }

    if (ps->pps != pps) {
        /* FIXME pps could trigger use a new sps */
        ret = update_pps_info(&ps->pps_info, pps, ps->sps);
        if (ret < 0) {
            goto failpps;
        }
        ps->pps = pps;
    }

    if (ps->ph != ph) {
        ret = update_ph_info(&ps->ph_info, ph);
        if (ret < 0) {
            goto failph;
        }
        ps->ph = ph;
    }

    /* We only check sh is present in ctx since we should be called
     * directely after the Slice Header is read so we are
     * reading a new slice
     */
    if (nvcl_ctx->sh) {
        ret = update_sh_info(&ps->sh_info, nvcl_ctx->sh);
        if (ret < 0) {
            goto failsh;
        }
        ps->sh = nvcl_ctx->sh;
    }


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
