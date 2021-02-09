#include <stdlib.h>
#include <string.h>
#include "vcl.h"
#include "dec_structures.h"
#include "cabac_internal.h"
#include "ctudec.h"

/*FIXME find a more global location for these defintions */
enum CUMode {
    OV_INTER_SKIP = 0,
    OV_INTER = 1,
    OV_INTRA = 2,
    OV_NA = 3,
    OV_MIP = 4,
};

struct OVMV {
    int32_t x;
    int32_t y;
};

/* FIXME shorten this LUT to min req size 
 * only used in truncated cabac reading
 */
static uint8_t g_tbMax[257] = { 0, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                                4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
                                5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
                                6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
                                6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7,
                                7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
                                7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
                                7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
                                7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
                                7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8 };

/* FIXME only used by mip_idx */
static inline uint8_t
vvc_get_cabac_truncated(OVCABACCtx *const cabac_ctx, unsigned int max_symbol){
    int threshold;
    uint32_t ruiSymbol = 0;
    /* MAX SYMBOL will not be > 16 */
    if( max_symbol > 256 ){
        int thres_val = 1 << 8;
        threshold = 8;
        while( thres_val <= max_symbol ){
            threshold++;
            thres_val <<= 1;
        }
        threshold--;
    }else{
        threshold = g_tbMax[max_symbol];
    }

    int val = 1 << threshold;
    int b = max_symbol - val;

    while(threshold--){
        ruiSymbol <<= 1;
        ruiSymbol |= ovcabac_bypass_read(cabac_ctx);
    }

    if( ruiSymbol >= val - b ){
        uint32_t uiSymbol;
        uiSymbol = ovcabac_bypass_read(cabac_ctx);
        ruiSymbol <<= 1;
        ruiSymbol += uiSymbol;
        ruiSymbol -= ( val - b );
    }

    return ruiSymbol;

}

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

static uint8_t
ovcabac_read_ae_root_cbf(OVCABACCtx *const cabac_ctx)
{
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    return ovcabac_ae_read(cabac_ctx, &cabac_state[QT_ROOT_CBF_CTX_OFFSET]);
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
        ctx_offset  = mip_abv == 75;
        ctx_offset += mip_lft == 75;
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
    int nb_mip_modes = 6;

    if ((log2_cb_h | log2_cb_w) == 2) {
        /* 4x4 */
        nb_mip_modes = 16;
    } else if (log2_cb_h == 2 || log2_cb_w == 2 || (log2_cb_h <= 3 && log2_cb_w <= 3)) {
        /* 8x8 || 4xX || Xx4 */
        nb_mip_modes = 8;
    }

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

#if 0
    /*TODO reset qp_map_up:left to current_qp at ctu start so there is no
     * testing required
     */
    int pred_qp = ((y0 ? ctu_dec->qp_map_up[x_cb]   : ctu_dec->qp_ctx.current_qp) +
                   (x0 ? ctu_dec->qp_map_left[y_cb] : ctu_dec->qp_ctx.current_qp) + 1) >> 1;

    ctu_dec->qp_ctx.current_qp = pred_qp;
#endif


    /* TODO set coding unit functor */
    cu = ctu_dec->coding_unit(ctu_dec, part_ctx, x0, y0, log2_cb_w, log2_cb_h);

    if (!(cu.cu_flags & flg_cu_skip_flag)) {
         /*TODO rename */
         transform_unit_wrap(ctu_dec, part_ctx,  x0, y0, log2_cb_w, log2_cb_h, cu);
    }

#if 0
    derive_dequant_ctx(ctu_dec, &ctu_dec->qp_ctx, 0);
#endif

#if 0
    /* update dqp for deblocking filter usage */
    if (!ctu_dec->dbf_disable) {
        fill_dbf_qp(&ctu_dec->dbf_info, x0, y0, log2_cb_w, log2_cb_h, ctu_dec->qp_ctx.current_qp);
        fill_dbf_qp_cb(&ctu_dec->dbf_info, x0, y0, log2_cb_w, log2_cb_h, ctu_dec->dequant_cb.qp - 12);
        fill_dbf_qp_cr(&ctu_dec->dbf_info, x0, y0, log2_cb_w, log2_cb_h, ctu_dec->dequant_cr.qp - 12);
    }
#endif

    // Update depth_maps to selected depths
#if 0
    /* update delta qp context */
    for (int i = 0; i < nb_cb_w; i++) {
        ctu_dec->qp_map_up[x_cb + i] = ctu_dec->qp_ctx.current_qp;
        ctu_dec->qp_map_up_cb[x_cb + i] = ctu_dec->dequant_cb.qp - 12;
        ctu_dec->qp_map_up_cr[x_cb + i] = ctu_dec->dequant_cr.qp - 12;
    }

    for (int i = 0; i < nb_cb_h; i++) {
        ctu_dec->qp_map_left[y_cb + i] = ctu_dec->qp_ctx.current_qp;
    }
#endif

#if 0
    if (!ctu_dec->dbf_disable) {
        fill_edge_map(&ctu_dec->dbf_info, x0, y0, log2_cb_w, log2_cb_h);
        fill_ctb_bound(&ctu_dec->dbf_info, x0, y0, log2_cb_w, log2_cb_h);
    }
#endif
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


        } else {
            uint8_t merge_flag = ovcabac_read_ae_cu_merge_flag(cabac_ctx);

            ctu_dec->prediction_unit(ctu_dec, part_ctx, x0, y0, log2_cu_w, log2_cu_h, merge_flag);

            cu_type = OV_INTER;

            FLG_STORE(merge_flag, cu.cu_flags);

        }
    }

    #if 0
    updt_cu_maps(ctu_dec, part_ctx, x0, y0, log2_cu_w, log2_cu_h, cu_type);
    #endif

    return cu;
}

VVCCU
coding_unit_intra_st(OVCTUDec *const ctu_dec,
                     const OVPartInfo *const part_ctx,
                     uint8_t x0, uint8_t y0,
                     uint8_t log2_cu_w, uint8_t log2_cu_h)
{
   VVCCU cu = {0};
   cu.cu_flags = 2;
   coding_unit_intra(ctu_dec, part_ctx, x0, y0, log2_cu_w, log2_cu_h);

   /* if not in separable tree */
   if (!ctu_dec->share) {
       coding_unit_intra_c(ctu_dec, ctu_dec->part_ctx_c, x0 >> 1, y0 >> 1,
               log2_cu_w - 1, log2_cu_h - 1);
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
    uint8_t intra_mode;
    uint8_t mip_flag = 0;
    uint8_t x_pu = x0 >> part_ctx->log2_min_cb_s;
    uint8_t y_pu = y0 >> part_ctx->log2_min_cb_s;
    uint8_t nb_pb_w = (1 << log2_cb_w) >> part_ctx->log2_min_cb_s;
    uint8_t nb_pb_h = (1 << log2_cb_h) >> part_ctx->log2_min_cb_s;
    uint8_t pu_shift = part_ctx->log2_min_cb_s - 2;
    VVCCU cu = {0};

    cu.cu_flags |= flg_pred_mode_flag;

    if (ctu_dec->enabled_mip) {
        #if 0
        VVCCTUPredContext *const pred_ctx = &ctu_dec->pred_ctx;
        #endif
        /* FIXME use other maps */
        uint8_t mip_abv = ctu_dec->part_map.cu_mode_x[x_pu] == OV_MIP;
        uint8_t mip_lft = ctu_dec->part_map.cu_mode_y[y_pu] == OV_MIP;

        mip_flag = ovcabac_read_ae_intra_mip(cabac_ctx, log2_cb_w, log2_cb_h,
                                             mip_abv, mip_lft);

        cu.cu_flags |= flg_mip_flag & (-(!!mip_flag));

        if (mip_flag) {
            struct PartMap *part_map = &ctu_dec->part_map;
            uint8_t transpose_flag = ovcabac_read_ae_intra_mip_transpose_flag(cabac_ctx);
            uint8_t mip_mode;

            /*FIXME use LUT based on log2_sizes would be a better option*/

            mip_mode = ovcabac_read_ae_intra_mip_mode(cabac_ctx, log2_cb_w, log2_cb_h);

            cu.cu_opaque  = transpose_flag << 7;
            cu.cu_opaque |= mip_mode;

            /*FIXME this is actually required by mip_flag CABAC context derivation */
            memset(&part_map->cu_mode_x[x_pu], OV_MIP, sizeof(uint8_t) * nb_pb_w);
            memset(&part_map->cu_mode_y[y_pu], OV_MIP, sizeof(uint8_t) * nb_pb_h);

            #if 0
            for (int i = 0; i < nb_pb_h; i++) {
                memset(&pred_ctx->cclm_intra_mode[x_pu + (i << 5) + (y_pu << 5)], OV_MIP,
                        sizeof(uint8_t) * nb_pb_w);
            }
            #endif
        }
    }

    /* FIXME missing IBC mode reading */
    if (!mip_flag) {

        uint8_t mrl_flag = ctu_dec->enable_mrl && y0 ? ovcabac_read_ae_intra_multiref_flag(cabac_ctx)
                                              : 0;
        uint8_t isp_mode = 0;
        uint8_t mpm_flag;
        uint8_t mrl_idx = mrl_flag;

        cu.cu_flags |= flg_mrl_flag & (-(!!mrl_flag));
        cu.cu_opaque = mrl_flag;

        if (!mrl_flag && ctu_dec->isp_enabled) {
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

#if 0
    cu = recon_intra_cu(ctu_dec, part_ctx, x0, y0, log2_cb_w, log2_cb_h, cu);
#endif

    return cu;
}

VVCCU
coding_unit_intra_c(OVCTUDec *const ctu_dec,
                    const OVPartInfo *const part_ctx,
                    uint8_t x0, uint8_t y0,
                    uint8_t log2_cb_w, uint8_t log2_cb_h)
{

    OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;
    #if 0
    VVCCTUPredContext *const pred_ctx = &ctu_dec->pred_ctx;
    #endif
    uint8_t intra_mode;
    uint8_t mpm_idx = 0, lm_idx = 1;
#if 0
    uint8_t y_pu = y0 >> part_ctx->log2_min_cb_s;
    uint8_t x_pu = x0 >> part_ctx->log2_min_cb_s;
    uint8_t nb_pb_w = (1 << log2_cb_w) >> part_ctx->log2_min_cb_s;
    uint8_t nb_pb_h = (1 << log2_cb_h) >> part_ctx->log2_min_cb_s;
    uint8_t mpm_flag = 0;
    uint8_t luma_mode = pred_ctx->cclm_intra_mode[(x_pu + ((y_pu + (nb_pb_h >> 1)) << 5) + (nb_pb_w >> 1))];
    uint8_t pu_shift = ctu_dec->part_ctx->log2_min_cb_s - 2;
#endif
    uint8_t cclm_flag = 0;

    VVCCU cu = {0};

    cu.cu_flags = 2;

    if (ctu_dec->lm_chroma_enabled && !ctu_dec->tmp_disable_cclm &&
            ctu_dec->enable_cclm == 1) {
        cclm_flag = ovcabac_read_ae_cclm_flag(cabac_ctx);
        if (cclm_flag) {
            lm_idx = ovcabac_read_ae_intra_lm_chroma(cabac_ctx);
        }
    }

    if (!cclm_flag) {
        uint8_t mpm_flag = ovcabac_read_ae_intra_chroma_mpm_flag(cabac_ctx);
        if (mpm_flag) {
            mpm_idx = ovcabac_read_ae_intra_chroma_mpm_idx(cabac_ctx);
        }
    }

#if 0
    //Note this is not required by cabac decoding.
    luma_mode = luma_mode != VVC_MIP_MODE ? luma_mode : VVC_PLANAR;
    intra_mode = ff_vvc_derive_intra_pred_mode_chroma(cclm_flag, mpm_flag, mpm_idx,
                                                      luma_mode, lm_idx);

    update_availability_maps(&ctu_dec->progress_map_c, x_pu << pu_shift, y_pu << pu_shift,
                             nb_pb_w << pu_shift, nb_pb_h << pu_shift);

    for (int i = 0; i < nb_pb_h; i++) {
        memset(&pred_ctx->intra_modes_chroma[x_pu + (y_pu << 5) + i*32], intra_mode,
                sizeof(uint8_t) * nb_pb_w);
    }

    if (intra_mode == VVC_DM_CHROMA) {//derive intra mode for use;
        uint16_t *const dst_cb = &ctu_dec->ctu_data_cb[VVC_CTB_OFFSET+(x0)+((y0)*VVC_CTB_STRIDE)];
        uint16_t *const dst_cr = &ctu_dec->ctu_data_cr[VVC_CTB_OFFSET+(x0)+((y0)*VVC_CTB_STRIDE)];
        vvc_intra_pred_chroma(ctu_dec, dst_cb, dst_cr,
                              VVC_CTB_STRIDE_CHROMA, luma_mode, x0, y0,
                              log2_cb_w, log2_cb_h);
    } else if (intra_mode == VVC_LM_CHROMA) {
        uint16_t *const dst_cb = &ctu_dec->ctu_data_cb[VVC_CTB_OFFSET+(x0)+((y0)*VVC_CTB_STRIDE)];
        uint16_t *const dst_cr = &ctu_dec->ctu_data_cr[VVC_CTB_OFFSET+(x0)+((y0)*VVC_CTB_STRIDE)];
        vvc_intra_pred_chroma(ctu_dec, dst_cb, dst_cr,
                              VVC_CTB_STRIDE_CHROMA, intra_mode, x0, y0,
                              log2_cb_w, log2_cb_h);
    } else {
        uint16_t *const dst_cb = &ctu_dec->ctu_data_cb[VVC_CTB_OFFSET+(x0)+((y0)*VVC_CTB_STRIDE)];
        uint16_t *const dst_cr = &ctu_dec->ctu_data_cr[VVC_CTB_OFFSET+(x0)+((y0)*VVC_CTB_STRIDE)];
        vvc_intra_pred_chroma(ctu_dec, dst_cb, dst_cr,
                              VVC_CTB_STRIDE_CHROMA, intra_mode, x0, y0,
                              log2_cb_w, log2_cb_h);
    }
#endif

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
#if 0
    struct VVCInterCtx *const inter_ctx = &ctu_dec->inter_ctx;
    struct VVCMVCtx *const mv_ctx0 = &inter_ctx->mv_ctx0;

    uint8_t y_pu = y0 >> part_ctx->log2_min_cb_s;
    uint8_t x_pu = x0 >> part_ctx->log2_min_cb_s;
    uint8_t nb_pb_w = (1 << log2_pb_w) >> part_ctx->log2_min_cb_s;
    uint8_t nb_pb_h = (1 << log2_pb_h) >> part_ctx->log2_min_cb_s;
#endif

    OVMV mv0;
    if (merge_flag) {
        uint8_t max_nb_cand = ctu_dec->max_num_merge_candidates;

        uint8_t merge_idx = ovcabac_read_ae_mvp_merge_idx(cabac_ctx, max_nb_cand);

#if 0
        mv0 = vvc_derive_merge_mvp(ctu_dec, inter_ctx, mv_ctx0,
                                   x_pu, y_pu, nb_pb_w, nb_pb_h,
                                   merge_idx, max_nb_cand);
#endif
    } else {
        /*FIXME add ref_idx*/
        OVMV mvd = ovcabac_read_ae_mvd(cabac_ctx);

        uint8_t mvp_idx = ovcabac_read_ae_mvp_flag(cabac_ctx);

#if 0
        mv0 = derive_mvp_candidates(ctu_dec, inter_ctx, mv_ctx0,
                                    x_pu, y_pu, nb_pb_w, nb_pb_h,
                                    mvp_idx, 1);
        mvd = scale_mvd(mvd);

        mv0.x += mvd.x;
        mv0.y += mvd.y;
#endif
    }

#if 0
    vvc_motion_compensation(ctu_dec, x0, y0, log2_pb_w, log2_pb_h, mv0, 0);

    fill_mvp_map(mv_ctx0, mv0, x_pu, y_pu, nb_pb_w, nb_pb_h);

    fill_dbf_mv_map(&ctu_dec->dbf_info, mv_ctx0, mv0, x_pu, y_pu, nb_pb_w, nb_pb_h);

    update_hmvp_lut(&inter_ctx->hmvp_lut, mv0);
#endif

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
    #if 0
    struct VVCInterCtx *const inter_ctx = &ctu_dec->inter_ctx;
    VVCMergeInfo mv_info;
    #endif

#if 0
    uint8_t y_pu = y0 >> part_ctx->log2_min_cb_s;
    uint8_t x_pu = x0 >> part_ctx->log2_min_cb_s;
    uint8_t nb_pb_w = (1 << log2_pb_w) >> part_ctx->log2_min_cb_s;
    uint8_t nb_pb_h = (1 << log2_pb_h) >> part_ctx->log2_min_cb_s;
#endif

    if (merge_flag) {
        uint8_t max_nb_cand = ctu_dec->max_num_merge_candidates;

        uint8_t merge_idx = ovcabac_read_ae_mvp_merge_idx(cabac_ctx, max_nb_cand);

#if 0
        mv_info = vvc_derive_merge_mvp_b(ctu_dec, inter_ctx, x_pu, y_pu,
                                         nb_pb_w, nb_pb_h, merge_idx,
                                         max_nb_cand);
#endif
    } else {
        uint8_t inter_dir;
        OVMV mvd0, mvd1;
        uint8_t mvp_idx0, mvp_idx1;

        inter_dir = ovcabac_read_ae_inter_dir(cabac_ctx, log2_pb_w, log2_pb_h);

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

#if 0
        mv_info = derive_mvp_b(ctu_dec, part_ctx, x0, y0, log2_pb_w, log2_pb_h,
                               mvd0, mvd1, mvp_idx0, mvp_idx1, inter_dir);
#endif
    }

#if 0
    update_mv_ctx_b(ctu_dec, inter_ctx, part_ctx, mv_info.mv0, mv_info.mv1, x0, y0,
                    log2_pb_w, log2_pb_h, mv_info.inter_dir);

    recon_mcp_b(ctu_dec, inter_ctx, part_ctx, mv_info.mv0, mv_info.mv1, x0, y0,
                log2_pb_w, log2_pb_h, mv_info.inter_dir);
#endif

    return merge_flag;
}
