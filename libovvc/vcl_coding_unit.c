#include <stdlib.h>
#include <string.h>
#include "vcl.h"
#include "ovutils.h"
#include "dec_structures.h"
#include "cabac_internal.h"
#include "ctudec.h"
#include "drv.h"
#include "drv_utils.h"
#include "rcn.h"
#include "dbf_utils.h"

/*FIXME find a more global location for these defintions */
enum CUMode {
    OV_NA = 0xFF,
    OV_INTER = 1,
    OV_INTRA = 2,
    OV_INTER_SKIP = 3,
    OV_MIP = 4,
};

/* FIXME refactor dequant*/
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

// /* FIXME only used by mip_idx */
// static inline uint8_t
// vvc_get_cabac_truncated(OVCABACCtx *const cabac_ctx, unsigned int max_symbol){
//     static const uint8_t threshold_lut[17] =
//     {
//         0, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 4
//     };
//     int threshold;
//     uint32_t ruiSymbol = 0;
//     /* MAX SYMBOL will not be > 16 */
//     #if 0
//     if( max_symbol > 256 ){
//         int thres_val = 1 << 8;
//         threshold = 8;
//         while( thres_val <= max_symbol ){
//             threshold++;
//             thres_val <<= 1;
//         }
//         threshold--;
//     }else{
//     #endif
//         threshold = threshold_lut[max_symbol];
//     #if 0
//     }
//     #endif

//     int val = 1 << threshold;
//     int b = max_symbol - val;

//     while(threshold--){
//         ruiSymbol <<= 1;
//         ruiSymbol |= ovcabac_bypass_read(cabac_ctx);
//     }

//     if( ruiSymbol >= val - b ){
//         uint32_t uiSymbol;
//         uiSymbol = ovcabac_bypass_read(cabac_ctx);
//         ruiSymbol <<= 1;
//         ruiSymbol += uiSymbol;
//         ruiSymbol -= ( val - b );
//     }

//     return ruiSymbol;
// }

/* skip_abv + skip_lft */
static uint8_t
ovcabac_read_ae_cu_skip_flag(OVCABACCtx *const cabac_ctx, uint8_t above_pu,
                             uint8_t left_pu)
{
    uint8_t cu_skip_flag;
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    int ctx_offset = (above_pu == OV_INTER_SKIP) + (left_pu == OV_INTER_SKIP);
    cu_skip_flag = ovcabac_ae_read(cabac_ctx, &cabac_state[SKIP_FLAG_CTX_OFFSET + ctx_offset]);
    return cu_skip_flag;
}

/* intra_abv | intra_lft */
static uint8_t
ovcabac_read_ae_pred_mode_flag(OVCABACCtx *const cabac_ctx,
                               uint8_t above_pu, uint8_t left_pu)
{
    uint8_t pred_mode_flag;
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    int ctx_offset = (above_pu == OV_INTRA) | (left_pu == OV_INTRA);
    pred_mode_flag = ovcabac_ae_read(cabac_ctx, &cabac_state[PRED_MODE_CTX_OFFSET + ctx_offset]);
    return pred_mode_flag;
}

static uint8_t
ovcabac_read_ae_cu_merge_flag(OVCABACCtx *const cabac_ctx)
{
    uint8_t merge_flag;
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    merge_flag = ovcabac_ae_read(cabac_ctx, &cabac_state[MERGE_FLAG_CTX_OFFSET]);
    return merge_flag;
}

/* FIXME check of max_nb_mrg_cand in function */
static uint8_t
ovcabac_read_ae_mvp_merge_idx(OVCABACCtx *const cabac_ctx,
                              uint8_t max_num_merge_cand)
{
    uint8_t merge_idx = 0;
    if (max_num_merge_cand) {
        uint64_t *const cabac_state = cabac_ctx->ctx_table;
        if (ovcabac_ae_read(cabac_ctx, &cabac_state[MERGE_IDX_CTX_OFFSET])) {
            merge_idx++;
            for (; merge_idx < max_num_merge_cand - 1; ++merge_idx) {
                if (!ovcabac_bypass_read(cabac_ctx)) {
                    break;
                }
            }
        }
    }
    return merge_idx;
}

static uint8_t
ovcabac_read_ae_inter_dir(OVCABACCtx *const cabac_ctx,
                          int log2_pb_w, int log2_pb_h)
{
    uint8_t inter_dir = 0;
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    /*FIXME Note no Bi predicition on 4x4 8x4 4x8*/
    if (log2_pb_w + log2_pb_h >= 5) {
        int ctx_id = 7 - ((log2_pb_w + log2_pb_h + 1) >> 1);
        inter_dir = ovcabac_ae_read(cabac_ctx, &cabac_state[INTER_DIR_CTX_OFFSET + ctx_id]);
        if (inter_dir) {
            return 3;
        }
    }

    return  1 + ovcabac_ae_read(cabac_ctx, &cabac_state[INTER_DIR_CTX_OFFSET + 5]);
}

/* FIXME optimize bypass CABAC reading for golomb + only used by mvd */
static int
vvc_exp_golomb_mv(OVCABACCtx *const cabac_ctx)
{
    unsigned prefix = 0;
    unsigned bit = 0;
    unsigned add_val = 0;
    unsigned length = 1, offset;

    do {
        prefix++;
        bit = ovcabac_bypass_read(cabac_ctx);
    } while (bit && prefix < (32 - 17));

    prefix -= 1 - bit;

    if (prefix < 0) {
        offset = prefix << 1;
    } else {
        offset = (((1 << prefix ) - 1) << 1);
        length += (prefix == (32 - 17) ? 17 - 1 : prefix);
    }

    while(length){
        add_val <<= 1;
        add_val |= ovcabac_bypass_read(cabac_ctx);
        length--;
    }

    return offset + add_val;
}

static OVMV
ovcabac_read_ae_mvd(OVCABACCtx *const cabac_ctx)
{
    int abs_x, abs_y;
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    OVMV mvd;

    abs_x = ovcabac_ae_read(cabac_ctx, &cabac_state[MVD_CTX_OFFSET]);
    abs_y = ovcabac_ae_read(cabac_ctx, &cabac_state[MVD_CTX_OFFSET]);

    if (abs_x) {
        abs_x += ovcabac_ae_read(cabac_ctx, &cabac_state[MVD_CTX_OFFSET + 1]);
    }

    if (abs_y) {
        abs_y += ovcabac_ae_read(cabac_ctx, &cabac_state[MVD_CTX_OFFSET + 1]);
    }

    if (abs_x) {
        uint8_t sign;

        if (abs_x > 1){
            abs_x += vvc_exp_golomb_mv(cabac_ctx);
        }

        sign = ovcabac_bypass_read(cabac_ctx);
        abs_x = sign ? -abs_x : abs_x;
    }

    if (abs_y) {
        uint8_t sign;

        if (abs_y > 1) {
            abs_y += vvc_exp_golomb_mv(cabac_ctx);
        }

        sign = ovcabac_bypass_read(cabac_ctx);
        abs_y = sign ? -abs_y : abs_y;
    }

    mvd.x = abs_x;
    mvd.y = abs_y;

    return mvd;
}

static uint8_t
ovcabac_read_ae_mvp_flag(OVCABACCtx *const cabac_ctx)
{
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    return ovcabac_ae_read(cabac_ctx, &cabac_state[MVP_IDX_CTX_OFFSET]);
}

/* mip_abv + mip_lft */
static uint8_t
ovcabac_read_ae_intra_mip(OVCABACCtx *const cabac_ctx,
                          int8_t log2_cb_w, int8_t log2_cb_h,
                          uint8_t mip_abv, uint8_t mip_lft)
{
    uint8_t ctx_offset;
    uint64_t *const cabac_state = cabac_ctx->ctx_table;

    if (abs(log2_cb_h - log2_cb_w) > 1){
        ctx_offset = 3;
    } else {
        ctx_offset  = mip_abv == OV_MIP;
        ctx_offset += mip_lft == OV_MIP;
    }

    return ovcabac_ae_read(cabac_ctx, &cabac_state[MIP_FLAG_CTX_OFFSET + ctx_offset]);
}

static uint8_t
ovcabac_read_ae_intra_mip_transpose_flag(OVCABACCtx *const cabac_ctx)
{
    return ovcabac_bypass_read(cabac_ctx);
}

static uint8_t
ovcabac_read_ae_intra_mip_mode(OVCABACCtx *const cabac_ctx, uint8_t log2_cb_w,
                               uint8_t log2_cb_h)
{
    #if 0
    int nb_mip_modes = 6;

    /* FIXME use LUT based on log2_sizes would be a better option */

    if ((log2_cb_h | log2_cb_w) == 2) {
        /* 4x4 */
        nb_mip_modes = 16;
    } else if (log2_cb_h == 2 || log2_cb_w == 2 || (log2_cb_h <= 3 && log2_cb_w <= 3)) {
        /* 8x8 || 4xX || Xx4 */
        nb_mip_modes = 8;
    }
    #else
    int nb_mip_modes = 6;
    if (log2_cb_h == log2_cb_w && log2_cb_h == 2) { //4x4
        nb_mip_modes = 16;
    } else if (log2_cb_h == 2 || log2_cb_w == 2 || (log2_cb_h == 3 && log2_cb_w == 3)) { //8x8
        nb_mip_modes = 8;
    }
    #endif

    return vvc_get_cabac_truncated(cabac_ctx, nb_mip_modes);
}

static uint8_t
ovcabac_read_ae_intra_multiref_flag(OVCABACCtx *const cabac_ctx)
{
    uint8_t value = 0;
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    value += ovcabac_ae_read(cabac_ctx, &cabac_state[MULTI_REF_LINE_IDX_CTX_OFFSET]);
    if (value) {
        value += (ovcabac_ae_read(cabac_ctx, &cabac_state[MULTI_REF_LINE_IDX_CTX_OFFSET + 1]));
    }
    return value;
}

/*FIXME check for allowed_isp before call */
static uint8_t
ovcabac_read_ae_intra_subpartition_flag(OVCABACCtx *const cabac_ctx,
                                        uint8_t allowed__direction_flags)
{
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    if (allowed__direction_flags && ovcabac_ae_read(cabac_ctx, &cabac_state[ISP_MODE_CTX_OFFSET])) {
        if (allowed__direction_flags == 3) {
            /* return 1 if hor 2 if ver */
            return 1 + ovcabac_ae_read(cabac_ctx, &cabac_state[ISP_MODE_CTX_OFFSET + 1]);
        } else {
            return allowed__direction_flags;
        }
    } else {
        return 0;
    }
}

static uint8_t
ovcabac_read_ae_intra_luma_mpm_flag(OVCABACCtx *const cabac_ctx)
{
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    return ovcabac_ae_read(cabac_ctx, &cabac_state[INTRA_LUMA_MPM_FLAG_CTX_OFFSET]);
}

static uint8_t
ovcabac_read_ae_intra_luma_mpm_idx(OVCABACCtx *const cabac_ctx,
                                   uint8_t is_multi_ref, uint8_t is_isp)
{
    uint8_t mpm_idx;
    /* FIXME use specific isp mrl cases ? */
    if (!is_multi_ref) {
        uint64_t *const cabac_state = cabac_ctx->ctx_table;
        mpm_idx = ovcabac_ae_read(cabac_ctx, &cabac_state[INTRA_LUMA_PLANAR_FLAG_CTX_OFFSET + !is_isp]);
    } else {
        mpm_idx = 1;
    }

    /* FIXME while + break loop instead */
    if (mpm_idx) {
        mpm_idx += ovcabac_bypass_read(cabac_ctx);
    }

    if (mpm_idx > 1) {
        mpm_idx += ovcabac_bypass_read(cabac_ctx);
    }

    if (mpm_idx > 2) {
        mpm_idx += ovcabac_bypass_read(cabac_ctx);
    }

    if (mpm_idx > 3) {
        mpm_idx += ovcabac_bypass_read(cabac_ctx);
    }

    return mpm_idx;
}

static uint8_t
ovcabac_read_ae_cclm_flag(OVCABACCtx *const cabac_ctx)
{
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    return ovcabac_ae_read(cabac_ctx, &cabac_state[CCLM_MODE_FLAG_CTX_OFFSET]);
}

static uint8_t
ovcabac_read_ae_intra_luma_mpm_remainder(OVCABACCtx *const cabac_ctx)
{
    uint8_t mpm_idx = ovcabac_bypass_read(cabac_ctx);
    int nb_bits = 5;//FIXME check behaviour

    while (--nb_bits) {
        mpm_idx = mpm_idx << 1;
        mpm_idx |= ovcabac_bypass_read(cabac_ctx);
    }

    if (mpm_idx >= (1 << 5) - (61 - 32)) {
        mpm_idx <<= 1;
        mpm_idx += ovcabac_bypass_read(cabac_ctx);
        mpm_idx -= (1 << 5) - (61 - 32);
    }

    return mpm_idx;
}

static uint8_t
ovcabac_read_ae_intra_chroma_mpm_flag(OVCABACCtx *const cabac_ctx)
{
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    return ovcabac_ae_read(cabac_ctx, &cabac_state[INTRA_CHROMA_PRED_MODE_CTX_OFFSET]);
}

static uint8_t
ovcabac_read_ae_intra_lm_chroma(OVCABACCtx *const cabac_ctx)
{
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    uint8_t lm_mode = ovcabac_ae_read(cabac_ctx, &cabac_state[CCLM_MODE_IDX_CTX_OFFSET]);

    if (lm_mode) {
        lm_mode += ovcabac_bypass_read(cabac_ctx);
    }

    return lm_mode;
}

static uint8_t
ovcabac_read_ae_intra_chroma_mpm_idx(OVCABACCtx *const cabac_ctx)
{
    uint8_t idx = ovcabac_bypass_read(cabac_ctx) << 1;

    return (idx | ovcabac_bypass_read(cabac_ctx));
}

int
coding_unit(OVCTUDec *const ctu_dec,
            const OVPartInfo *const part_ctx,
            uint8_t x0, uint8_t y0,
            uint8_t log2_cb_w, uint8_t log2_cb_h)
{
    unsigned int nb_cb_w = 1 << log2_cb_w >> part_ctx->log2_min_cb_s;
    unsigned int nb_cb_h = 1 << log2_cb_h >> part_ctx->log2_min_cb_s;

    int x_cb = x0 >> part_ctx->log2_min_cb_s;
    int y_cb = y0 >> part_ctx->log2_min_cb_s;

    VVCCU cu;

    int pred_qp = ((y0 ? ctu_dec->drv_ctx.qp_map_x[x_cb] : ctu_dec->qp_ctx.current_qp) +
                   (x0 ? ctu_dec->drv_ctx.qp_map_y[y_cb] : ctu_dec->qp_ctx.current_qp) + 1) >> 1;

    ctu_dec->qp_ctx.current_qp = pred_qp;

    cu = ctu_dec->coding_unit(ctu_dec, part_ctx, x0, y0, log2_cb_w, log2_cb_h);

    ctu_dec->dequant_chroma = &ctu_dec->dequant_cb;

    if (!(cu.cu_flags & flg_cu_skip_flag)) {
         /*TODO rename */
         transform_unit_wrap(ctu_dec, part_ctx,  x0, y0, log2_cb_w, log2_cb_h, cu);
    }

    /* FIXME delta qp clean */
    derive_dequant_ctx(ctu_dec, &ctu_dec->qp_ctx, 0);

    /* FIXME memset instead
     */
    /* update delta qp context */
    for (int i = 0; i < nb_cb_w; i++) {
        ctu_dec->drv_ctx.qp_map_x[x_cb + i] = ctu_dec->qp_ctx.current_qp;
    }

    for (int i = 0; i < nb_cb_h; i++) {
        ctu_dec->drv_ctx.qp_map_y[y_cb + i] = ctu_dec->qp_ctx.current_qp;
    }

    /* update dqp for deblocking filter usage */
    if (ctu_dec->coding_tree != dual_tree) {
        struct DBFInfo *dbf_info = &ctu_dec->dbf_info;
        fill_edge_map(&ctu_dec->dbf_info, x0, y0, log2_cb_w, log2_cb_h);
        fill_ctb_bound(&ctu_dec->dbf_info, x0, y0, log2_cb_w, log2_cb_h);


        fill_edge_map_c(&ctu_dec->dbf_info, x0, y0, log2_cb_w, log2_cb_h);
        fill_ctb_bound_c(&ctu_dec->dbf_info, x0, y0, log2_cb_w, log2_cb_h);

        dbf_fill_qp_map(&dbf_info->qp_map_y, x0, y0, log2_cb_w, log2_cb_h, ctu_dec->qp_ctx.current_qp);
        dbf_fill_qp_map(&dbf_info->qp_map_cb, x0, y0, log2_cb_w, log2_cb_h, ctu_dec->dequant_cb.qp - 12);
        dbf_fill_qp_map(&dbf_info->qp_map_cr, x0, y0, log2_cb_w, log2_cb_h, ctu_dec->dequant_cr.qp - 12);
    } else {
    if (!ctu_dec->dbf_disable && !(cu.cu_flags & flg_isp_flag)) {
        struct DBFInfo *dbf_info = &ctu_dec->dbf_info;
        dbf_fill_qp_map(&dbf_info->qp_map_y, x0, y0, log2_cb_w, log2_cb_h, ctu_dec->qp_ctx.current_qp);
        if (ctu_dec->coding_unit == &coding_unit_intra) {
            fill_edge_map(&ctu_dec->dbf_info, x0, y0, log2_cb_w, log2_cb_h);
            fill_ctb_bound(&ctu_dec->dbf_info, x0, y0, log2_cb_w, log2_cb_h);
        } else if (ctu_dec->coding_unit == &coding_unit_intra_c) {
            if (ctu_dec->dequant_chroma == &ctu_dec->dequant_joint_cb_cr) {
            dbf_fill_qp_map(&dbf_info->qp_map_cb, x0 << 1, y0 << 1, log2_cb_w + 1, log2_cb_h + 1, ctu_dec->dequant_joint_cb_cr.qp - 12);
            dbf_fill_qp_map(&dbf_info->qp_map_cr, x0 << 1, y0 << 1, log2_cb_w + 1, log2_cb_h + 1, ctu_dec->dequant_joint_cb_cr.qp - 12);
            } else {
            dbf_fill_qp_map(&dbf_info->qp_map_cb, x0 << 1, y0 << 1, log2_cb_w + 1, log2_cb_h + 1, ctu_dec->dequant_cb.qp - 12);
            dbf_fill_qp_map(&dbf_info->qp_map_cr, x0 << 1, y0 << 1, log2_cb_w + 1, log2_cb_h + 1, ctu_dec->dequant_cr.qp - 12);
            }
            fill_edge_map_c(&ctu_dec->dbf_info, x0 << 1, y0 << 1, log2_cb_w + 1, log2_cb_h + 1);
            fill_ctb_bound_c(&ctu_dec->dbf_info, x0 << 1, y0 << 1, log2_cb_w + 1, log2_cb_h + 1);

        } else {
            fill_edge_map(&ctu_dec->dbf_info, x0, y0, log2_cb_w, log2_cb_h);
            fill_ctb_bound(&ctu_dec->dbf_info, x0, y0, log2_cb_w, log2_cb_h);

            fill_edge_map_c(&ctu_dec->dbf_info, x0, y0, log2_cb_w, log2_cb_h);
            fill_ctb_bound_c(&ctu_dec->dbf_info, x0, y0, log2_cb_w, log2_cb_h);

            dbf_fill_qp_map(&dbf_info->qp_map_cb, x0, y0, log2_cb_w, log2_cb_h, ctu_dec->dequant_cb.qp - 12);
            dbf_fill_qp_map(&dbf_info->qp_map_cr, x0, y0, log2_cb_w, log2_cb_h, ctu_dec->dequant_cr.qp - 12);
            dbf_fill_qp_map(&dbf_info->qp_map_cb, x0, y0, log2_cb_w, log2_cb_h, ctu_dec->dequant_joint_cb_cr.qp - 12);
            dbf_fill_qp_map(&dbf_info->qp_map_cr, x0, y0, log2_cb_w, log2_cb_h, ctu_dec->dequant_joint_cb_cr.qp - 12);
        }
    } else {
        fill_ctb_bound(&ctu_dec->dbf_info, x0, y0, log2_cb_w, log2_cb_h);
    }
    }

    // Update depth_maps to selected depths
    struct PartMap *const part_map = ctu_dec->active_part_map;

    memset(&part_map->log2_cu_w_map_x[x_cb], log2_cb_w, sizeof(uint8_t) * nb_cb_w);
    memset(&part_map->log2_cu_h_map_y[y_cb], log2_cb_h, sizeof(uint8_t) * nb_cb_h);

    /* if we are in single tree and outside of separable tree we also need to
     * update chroma partition context
     */
    if (!ctu_dec->share &&
        ctu_dec->coding_tree != &dual_tree &&
        ctu_dec->coding_tree_implicit != &dual_tree_implicit
       ) {
        /* Note since we are in single tree the log2_diff_size of cu are the same no need
         * to recompute the nb_cb and cb coordinates
         */
         /*FIXME check if using log2_cb_size directely is OK for chroma tree splits*/
        memset(&ctu_dec->part_map_c.log2_cu_w_map_x[x_cb], log2_cb_w, sizeof(uint8_t) * nb_cb_w);
        memset(&ctu_dec->part_map_c.log2_cu_h_map_y[y_cb], log2_cb_h, sizeof(uint8_t) * nb_cb_h);
    }
    return 0;
}

static void
updt_cu_maps(OVCTUDec *const ctudec,
             const OVPartInfo *const part_ctx,
             unsigned int x0, unsigned int y0,
             unsigned int log2_cu_w, unsigned int log2_cu_h,
             uint8_t cu_mode)
{
    /*FIXME Separate non CABAC ctx and reconstruction buffers */

    struct PartMap *const part_map = &ctudec->part_map;

    uint8_t y_pu = y0 >> part_ctx->log2_min_cb_s;
    uint8_t x_pu = x0 >> part_ctx->log2_min_cb_s;
    uint8_t nb_pb_w = (1 << log2_cu_w) >> part_ctx->log2_min_cb_s;
    uint8_t nb_pb_h = (1 << log2_cu_h) >> part_ctx->log2_min_cb_s;

    memset(&part_map->cu_mode_x[x_pu], (uint8_t)cu_mode, sizeof(uint8_t) * nb_pb_w);
    memset(&part_map->cu_mode_y[y_pu], (uint8_t)cu_mode, sizeof(uint8_t) * nb_pb_h);
}

VVCCU
coding_unit_inter_st(OVCTUDec *const ctu_dec,
                     const OVPartInfo *const part_ctx,
                     uint8_t x0, uint8_t y0,
                     uint8_t log2_cu_w, uint8_t log2_cu_h)
{
    OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;
    uint8_t y_pu = y0 >> part_ctx->log2_min_cb_s;
    uint8_t x_pu = x0 >> part_ctx->log2_min_cb_s;
    #if 0
    VVCCTUPredContext *const pred_ctx = &ctu_dec->pred_ctx;
    #endif
    /*TODO fill */
    uint8_t cu_type_abv = ctu_dec->part_map.cu_mode_x[x_pu];
    uint8_t cu_type_lft = ctu_dec->part_map.cu_mode_y[y_pu];
    uint8_t cu_skip_flag;
    uint8_t cu_type;
    VVCCU cu = {0};

    cu_skip_flag = ovcabac_read_ae_cu_skip_flag(cabac_ctx, cu_type_abv,
                                                cu_type_lft);

    if (cu_skip_flag) {
        /* FIXME cu_skip_flag activation force merge_flag so we only need to read
           merge_idx */
        ctu_dec->prediction_unit(ctu_dec, part_ctx, x0, y0, log2_cu_w, log2_cu_h, 1);

        cu_type = OV_INTER_SKIP;

        FLG_STORE(cu_skip_flag, cu.cu_flags);

    } else {
        uint8_t pred_mode_flag;

        pred_mode_flag = ovcabac_read_ae_pred_mode_flag(cabac_ctx, cu_type_abv,
                                                        cu_type_lft);

        FLG_STORE(pred_mode_flag, cu.cu_flags);

        if (pred_mode_flag) {
            coding_unit_intra_st(ctu_dec, part_ctx, x0, y0, log2_cu_w, log2_cu_h);

            cu_type = OV_INTRA;

            fill_bs_map(&ctu_dec->dbf_info.bs2_map, x0, y0, log2_cu_w, log2_cu_h);
            fill_bs_map(&ctu_dec->dbf_info.bs2_map_c, x0, y0, log2_cu_w, log2_cu_h);

        } else {
            uint8_t merge_flag = ovcabac_read_ae_cu_merge_flag(cabac_ctx);

            ctu_dec->prediction_unit(ctu_dec, part_ctx, x0, y0, log2_cu_w, log2_cu_h, merge_flag);

            cu_type = OV_INTER;

            FLG_STORE(merge_flag, cu.cu_flags);

        }
    }

    updt_cu_maps(ctu_dec, part_ctx, x0, y0, log2_cu_w, log2_cu_h, cu_type);

    return cu;
}

VVCCU
coding_unit_intra_st(OVCTUDec *const ctu_dec,
                     const OVPartInfo *const part_ctx,
                     uint8_t x0, uint8_t y0,
                     uint8_t log2_cu_w, uint8_t log2_cu_h)
{
   VVCCU cu = {0};

   /* Force pred_mode_flag to 2 so we know cu was intra */
   cu.cu_flags = 2;

   coding_unit_intra(ctu_dec, part_ctx, x0, y0, log2_cu_w, log2_cu_h);

   fill_bs_map(&ctu_dec->dbf_info.bs2_map, x0, y0, log2_cu_w, log2_cu_h);
   /* if not in separable tree */
   if (!ctu_dec->share) {
       coding_unit_intra_c(ctu_dec, ctu_dec->part_ctx_c, x0 >> 1, y0 >> 1,
                           log2_cu_w - 1, log2_cu_h - 1);

       fill_bs_map(&ctu_dec->dbf_info.bs2_map_c, x0, y0, log2_cu_w, log2_cu_h);
   }


   return cu;
}

VVCCU
coding_unit_intra(OVCTUDec *const ctu_dec,
                  const OVPartInfo *const part_ctx,
                  uint8_t x0, uint8_t y0,
                  uint8_t log2_cb_w, uint8_t log2_cb_h)
{
    OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;
    #if 0
    uint8_t intra_mode;
    #endif
    uint8_t mip_flag = 0;
    uint8_t x_pu = x0 >> part_ctx->log2_min_cb_s;
    uint8_t y_pu = y0 >> part_ctx->log2_min_cb_s;
    uint8_t nb_pb_w = (1 << log2_cb_w) >> part_ctx->log2_min_cb_s;
    uint8_t nb_pb_h = (1 << log2_cb_h) >> part_ctx->log2_min_cb_s;
    #if 0
    uint8_t pu_shift = part_ctx->log2_min_cb_s - 2;
    #endif
    VVCCU cu = {0};

    cu.cu_flags |= flg_pred_mode_flag;

    if (ctu_dec->enabled_mip) {
        struct PartMap *part_map = &ctu_dec->part_map;
        uint8_t mip_abv = part_map->cu_mode_x[x_pu];
        uint8_t mip_lft = part_map->cu_mode_y[y_pu];

        mip_flag = ovcabac_read_ae_intra_mip(cabac_ctx, log2_cb_w, log2_cb_h,
                                             mip_abv, mip_lft);

        cu.cu_flags |= flg_mip_flag & (-(!!mip_flag));

        if (mip_flag) {
            uint8_t transpose_flag = ovcabac_read_ae_intra_mip_transpose_flag(cabac_ctx);
            uint8_t mip_mode;

            mip_mode = ovcabac_read_ae_intra_mip_mode(cabac_ctx, log2_cb_w, log2_cb_h);

            cu.cu_opaque  = transpose_flag << 7;
            cu.cu_opaque |= mip_mode;

            memset(&part_map->cu_mode_x[x_pu], OV_MIP, sizeof(uint8_t) * nb_pb_w);
            memset(&part_map->cu_mode_y[y_pu], OV_MIP, sizeof(uint8_t) * nb_pb_h);

            /* FIXME Check default to PLANAR for modes derivation */
        } else {
            memset(&part_map->cu_mode_x[x_pu], OV_INTRA, sizeof(uint8_t) * nb_pb_w);
            memset(&part_map->cu_mode_y[y_pu], OV_INTRA, sizeof(uint8_t) * nb_pb_h);
        }
    }

    /* FIXME missing IBC + PLT mode reading */
    if (!mip_flag) {

        uint8_t mrl_flag = ctu_dec->enable_mrl && y0 ? ovcabac_read_ae_intra_multiref_flag(cabac_ctx)
                                              : 0;
        uint8_t isp_mode = 0;
        uint8_t mpm_flag;
        uint8_t mrl_idx = mrl_flag;

        cu.cu_flags |= flg_mrl_flag & (-(!!mrl_flag));
        cu.cu_opaque = mrl_idx;

        if (!mrl_flag && ctu_dec->isp_enabled) {
            /* FIXME use LUTs based on sizes for split status */
            uint8_t isp_split_status = (log2_cb_w + log2_cb_h) > (2 << 1);

            if (isp_split_status) {
                isp_split_status  = (log2_cb_w <=  part_ctx->log2_max_tb_s) << 1;
                isp_split_status |= log2_cb_h <=  part_ctx->log2_max_tb_s;
            }

            isp_mode = ovcabac_read_ae_intra_subpartition_flag(cabac_ctx, isp_split_status);

            cu.cu_flags |= flg_isp_flag & (-(!!isp_mode));

            cu.cu_opaque = isp_mode;
        }

        mpm_flag = mrl_flag || ovcabac_read_ae_intra_luma_mpm_flag(cabac_ctx);

        cu.cu_flags |= flg_mpm_flag & (-(!!mpm_flag));

        if (mpm_flag) {
            uint8_t mpm_idx = ovcabac_read_ae_intra_luma_mpm_idx(cabac_ctx, mrl_flag, !!isp_mode);

            cu.cu_mode_info = mpm_idx;

        } else {
            uint8_t mpm_rem = ovcabac_read_ae_intra_luma_mpm_remainder(cabac_ctx);

            cu.cu_mode_info = mpm_rem;
        }
    }

    /* FIXME move after TU is read so we can reconstruct with or without
     * transform trees
     */
    cu = drv_intra_cu(ctu_dec, part_ctx, x0, y0, log2_cb_w, log2_cb_h, cu);
    fill_bs_map(&ctu_dec->dbf_info.bs2_map, x0, y0, log2_cb_w, log2_cb_h);

    return cu;
}

VVCCU
coding_unit_intra_c(OVCTUDec *const ctu_dec,
                    const OVPartInfo *const part_ctx,
                    uint8_t x0, uint8_t y0,
                    uint8_t log2_cb_w, uint8_t log2_cb_h)
{

    OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;
    struct IntraDRVInfo *const i_info = &ctu_dec->drv_ctx.intra_info;
    /* TODO set mode to default */
    uint8_t intra_mode;
    uint8_t mpm_idx = 0, cclm_idx = 1;
#if 1
    /* FIXME move to drv */
    uint8_t y_pu = y0 >> part_ctx->log2_min_cb_s;
    uint8_t x_pu = x0 >> part_ctx->log2_min_cb_s;
    uint8_t nb_pb_w = (1 << log2_cb_w) >> part_ctx->log2_min_cb_s;
    uint8_t nb_pb_h = (1 << log2_cb_h) >> part_ctx->log2_min_cb_s;
    uint8_t luma_mode = i_info->luma_modes[(x_pu + ((y_pu + (nb_pb_h >> 1)) << 5) + (nb_pb_w >> 1))];
    uint8_t pu_shift = ctu_dec->part_ctx->log2_min_cb_s - 2;
#endif
    uint8_t cclm_flag = 0;
    uint8_t mpm_flag  = 0;

    VVCCU cu = {0};

    cu.cu_flags = 2;

    /* FIXME CCLM luma partition constraints */
    if (ctu_dec->lm_chroma_enabled && !ctu_dec->tmp_disable_cclm &&
        ctu_dec->enable_cclm == 1) {

        cclm_flag = ovcabac_read_ae_cclm_flag(cabac_ctx);

        if (cclm_flag) {
            cclm_idx = ovcabac_read_ae_intra_lm_chroma(cabac_ctx);
        }
    }

    if (!cclm_flag) {
        mpm_flag = ovcabac_read_ae_intra_chroma_mpm_flag(cabac_ctx);
        if (mpm_flag) {
            mpm_idx = ovcabac_read_ae_intra_chroma_mpm_idx(cabac_ctx);
        }
    }

    //Note this is not required by cabac decoding.
    intra_mode = derive_intra_mode_c(cclm_flag, mpm_flag, mpm_idx,
                                     luma_mode, cclm_idx);


    /* FIXME move to RCN */
    ctu_field_set_rect_bitfield(&ctu_dec->rcn_ctx.progress_field_c, x_pu << pu_shift,
                                y_pu << pu_shift, nb_pb_w << pu_shift, nb_pb_h << pu_shift);


    ctu_dec->intra_mode_c = intra_mode;
    vvc_intra_pred_chroma(&ctu_dec->rcn_ctx, intra_mode, x0, y0, log2_cb_w, log2_cb_h);

    fill_bs_map(&ctu_dec->dbf_info.bs2_map_c, x0 << 1, y0 << 1, log2_cb_w + 1, log2_cb_h + 1);
    return cu;
}

int
prediction_unit_inter_p(OVCTUDec *const ctu_dec,
                        const OVPartInfo *const part_ctx,
                        uint8_t x0, uint8_t y0,
                        uint8_t log2_pb_w, uint8_t log2_pb_h,
                        uint8_t merge_flag)
{
    OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;
    struct IntraDRVInfo *const i_info = &ctu_dec->drv_ctx.intra_info;
#if 1
    struct InterDRVCtx *const inter_ctx = &ctu_dec->drv_ctx.inter_ctx;
    struct OVMVCtx *const mv_ctx0 = &inter_ctx->mv_ctx0;

#endif
    uint8_t y_pu = y0 >> part_ctx->log2_min_cb_s;
    uint8_t x_pu = x0 >> part_ctx->log2_min_cb_s;
    uint8_t nb_pb_w = (1 << log2_pb_w) >> part_ctx->log2_min_cb_s;
    uint8_t nb_pb_h = (1 << log2_pb_h) >> part_ctx->log2_min_cb_s;

    OVMV mv0;

    if (merge_flag) {
        uint8_t max_nb_cand = ctu_dec->max_num_merge_candidates;

        uint8_t merge_idx = ovcabac_read_ae_mvp_merge_idx(cabac_ctx, max_nb_cand);

        mv0 = drv_merge_mvp(inter_ctx, mv_ctx0,
                            x_pu, y_pu, nb_pb_w, nb_pb_h,
                            merge_idx, max_nb_cand);

    } else {
        /*FIXME add ref_idx*/
        OVMV mvd = ovcabac_read_ae_mvd(cabac_ctx);

        uint8_t mvp_idx = ovcabac_read_ae_mvp_flag(cabac_ctx);

        mv0 = drv_mvp_mvd(inter_ctx, mv_ctx0, mvd,
                          x_pu, y_pu, nb_pb_w, nb_pb_h,
                          mvp_idx, 1);
    }

    rcn_mcp(ctu_dec, x0, y0, log2_pb_w, log2_pb_h, mv0, 0);

    uint8_t pu_shift = part_ctx->log2_min_cb_s - 2;

    ctu_field_set_rect_bitfield(&ctu_dec->rcn_ctx.progress_field_c, x_pu << pu_shift,
                                y_pu << pu_shift, nb_pb_w << pu_shift, nb_pb_h << pu_shift);
    ctu_field_set_rect_bitfield(&ctu_dec->rcn_ctx.progress_field, x_pu << pu_shift,
                                y_pu << pu_shift, nb_pb_w << pu_shift, nb_pb_h << pu_shift);
#if 0
    fill_dbf_mv_map(&ctu_dec->dbf_info, mv_ctx0, mv0, x_pu, y_pu, nb_pb_w, nb_pb_h);
#endif

    /*FIXME this have to be moved to DRV */
    /* We need to reset Intra mode maps to PLANAR for correct MPM derivation */
    memset(&i_info->luma_mode_x[x_pu], OVINTRA_PLANAR, sizeof(uint8_t) * nb_pb_w);
    memset(&i_info->luma_mode_y[y_pu], OVINTRA_PLANAR, sizeof(uint8_t) * nb_pb_h);

    for (int i = 0; i < nb_pb_h; i++) {
        memset(&i_info->luma_modes[x_pu + (i << 5) + (y_pu << 5)], OVINTRA_PLANAR,
               sizeof(uint8_t) * nb_pb_w);
    }

    return merge_flag;
}

int
prediction_unit_inter_b(OVCTUDec *const ctu_dec,
                        const OVPartInfo *const part_ctx,
                        uint8_t x0, uint8_t y0,
                        uint8_t log2_pb_w, uint8_t log2_pb_h,
                        uint8_t merge_flag)
{
    OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;
    #if 1
    struct InterDRVCtx *const inter_ctx = &ctu_dec->drv_ctx.inter_ctx;
    struct IntraDRVInfo *const i_info = &ctu_dec->drv_ctx.intra_info;
    VVCMergeInfo mv_info;
    #endif

#if 1
    uint8_t y_pu = y0 >> part_ctx->log2_min_cb_s;
    uint8_t x_pu = x0 >> part_ctx->log2_min_cb_s;
    uint8_t nb_pb_w = (1 << log2_pb_w) >> part_ctx->log2_min_cb_s;
    uint8_t nb_pb_h = (1 << log2_pb_h) >> part_ctx->log2_min_cb_s;
#endif

    if (merge_flag) {
        uint8_t max_nb_cand = ctu_dec->max_num_merge_candidates;

        uint8_t merge_idx = ovcabac_read_ae_mvp_merge_idx(cabac_ctx, max_nb_cand);

        mv_info = drv_merge_mvp_b(inter_ctx, x_pu, y_pu,
                                  nb_pb_w, nb_pb_h, merge_idx,
                                  max_nb_cand);

    } else {
        OVMV mvd0, mvd1;
        uint8_t mvp_idx0 = 0;
        uint8_t mvp_idx1 = 0;

        uint8_t inter_dir = ovcabac_read_ae_inter_dir(cabac_ctx, log2_pb_w, log2_pb_h);

        if (inter_dir & 0x1) {
            /*FIXME add ref_idx*/
            mvd0 = ovcabac_read_ae_mvd(cabac_ctx);

            mvp_idx0 = ovcabac_read_ae_mvp_flag(cabac_ctx);
        }

        if (inter_dir & 0x2) {
            /*FIXME add ref_idx*/
            mvd1 = ovcabac_read_ae_mvd(cabac_ctx);

            mvp_idx1 = ovcabac_read_ae_mvp_flag(cabac_ctx);
        }

        mv_info = drv_mvp_b(inter_ctx, x_pu, y_pu, nb_pb_w, nb_pb_h,
                            mvd0, mvd1, mvp_idx0, mvp_idx1, inter_dir);
    }

    rcn_mcp_b(ctu_dec, inter_ctx, part_ctx, mv_info.mv0, mv_info.mv1, x0, y0,
              log2_pb_w, log2_pb_h, mv_info.inter_dir);

    uint8_t pu_shift = part_ctx->log2_min_cb_s - 2;

    ctu_field_set_rect_bitfield(&ctu_dec->rcn_ctx.progress_field_c, x_pu << pu_shift,
                                y_pu << pu_shift, nb_pb_w << pu_shift, nb_pb_h << pu_shift);
    ctu_field_set_rect_bitfield(&ctu_dec->rcn_ctx.progress_field, x_pu << pu_shift,
                                y_pu << pu_shift, nb_pb_w << pu_shift, nb_pb_h << pu_shift);

    /*FIXME this has to be moved to DRV */
    /* We need to reset Intra mode maps to PLANAR for correct MPM derivation */
    memset(&i_info->luma_mode_x[x_pu], OVINTRA_PLANAR, sizeof(uint8_t) * nb_pb_w);
    memset(&i_info->luma_mode_y[y_pu], OVINTRA_PLANAR, sizeof(uint8_t) * nb_pb_h);

    for (int i = 0; i < nb_pb_h; i++) {
        memset(&i_info->luma_modes[x_pu + (i << 5) + (y_pu << 5)], OVINTRA_PLANAR,
               sizeof(uint8_t) * nb_pb_w);
    }

    return merge_flag;
}
