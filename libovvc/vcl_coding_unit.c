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
    OV_AFFINE = 5,
    OV_INTER_SKIP_AFFINE = 6,
};

#define BCW_NUM                 5 ///< the number of weight options
#define BCW_DEFAULT             ((uint8_t)(BCW_NUM >> 1)) ///< Default weighting index representing for w=0.5
#define BCW_SIZE_CONSTRAINT     256 ///< disabling Bcw if cu size is smaller than 256

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
    ctudec->dequant_joint_cb_cr.qp = qp_ctx->chroma_qp_map_jcbcr[(base_qp + 64) & 63] + qp_ctx->jcbcr_offset + 12;
    ctudec->qp_ctx.current_qp = base_qp;
}

/* skip_abv + skip_lft */
static uint8_t
ovcabac_read_ae_cu_skip_flag(OVCABACCtx *const cabac_ctx, uint8_t above_pu,
                             uint8_t left_pu)
{
    uint8_t cu_skip_flag;
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    int ctx_offset = (above_pu == OV_INTER_SKIP || above_pu == OV_INTER_SKIP_AFFINE) + (left_pu == OV_INTER_SKIP || left_pu == OV_INTER_SKIP_AFFINE);
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
    int ctx_offset = (above_pu == OV_INTRA) | (left_pu == OV_INTRA) | (above_pu == OV_MIP) | (left_pu == OV_MIP);
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

static uint8_t
ovcabac_read_ae_sb_merge_flag(OVCABACCtx *const cabac_ctx, uint8_t lft_affine, uint8_t abv_affine)
{
    uint8_t sb_merge_flag;
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    uint8_t ctx_offset = abv_affine + lft_affine;
    sb_merge_flag = ovcabac_ae_read(cabac_ctx, &cabac_state[SUBBLOCK_MERGE_FLAG_CTX_OFFSET + ctx_offset]);
    return sb_merge_flag;
}


static uint8_t
ovcabac_read_ae_cu_affine_flag(OVCABACCtx *const cabac_ctx, uint8_t lft_affine, uint8_t abv_affine)
{
    uint8_t affine_flag;
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    uint8_t ctx_offset = abv_affine + lft_affine;
    affine_flag = ovcabac_ae_read(cabac_ctx, &cabac_state[AFFINE_FLAG_CTX_OFFSET + ctx_offset]);
    return affine_flag;
}

static uint8_t
ovcabac_read_ae_cu_affine_type(OVCABACCtx *const cabac_ctx)
{
    uint8_t affine_type;
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    affine_type = ovcabac_ae_read(cabac_ctx, &cabac_state[AFFINE_TYPE_CTX_OFFSET]);
    return affine_type;
}

static uint8_t
ovcabac_read_ae_affine_merge_idx(OVCABACCtx *const cabac_ctx, uint8_t nb_affine_merge_cand_min1)
{
    uint8_t merge_idx = 0;
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    if (nb_affine_merge_cand_min1 > 0) {
        if (ovcabac_ae_read(cabac_ctx, &cabac_state[AFF_MERGE_IDX_CTX_OFFSET])) {
            do {
                ++merge_idx;
            } while (--nb_affine_merge_cand_min1 && ovcabac_bypass_read(cabac_ctx));
        }
    }

    return merge_idx;
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

uint8_t
ovcabac_read_ae_mmvd_merge_idx(OVCABACCtx *const cabac_ctx,
                              uint8_t max_num_merge_cand)
{
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    uint8_t var0 = 0;
    if (max_num_merge_cand  > 1){
        var0 = ovcabac_ae_read(cabac_ctx, &cabac_state[MMVD_MERGE_IDX_CTX_OFFSET]);
    }

    int num_cand_minus1 = MMVD_REFINE_STEP - 1;
    int var1 = 0;
    if (ovcabac_ae_read(cabac_ctx, &cabac_state[MMVD_STEP_MVP_IDX_CTX_OFFSET])){
        var1++;
        for (; var1 < num_cand_minus1; var1++){
            if (!ovcabac_bypass_read(cabac_ctx)) {
                break;
            }
        }
    }
    int var2 = 0;
    if (ovcabac_bypass_read(cabac_ctx)){
    var2 += 2;
        if (ovcabac_bypass_read(cabac_ctx)){
            var2 += 1;
        }
    }
    else{
        var2 += 0;
        if (ovcabac_bypass_read(cabac_ctx)){
            var2 += 1;
        }
    }
    return (var0 * MMVD_MAX_REFINE_NUM + var1 * 4 + var2);
}

void
ovcabac_read_ae_gpm_merge_idx(OVCABACCtx *const cabac_ctx, struct VVCGPM* gpm_ctx,
                              uint8_t max_num_geo_cand)
{
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    gpm_ctx->split_dir = vvc_get_cabac_truncated(cabac_ctx, GEO_NUM_PARTITION_MODE);

    int num_cand_min2 = max_num_geo_cand  - 2;
    gpm_ctx->merge_idx0    = 0;
    gpm_ctx->merge_idx1    = 0;
    if (ovcabac_ae_read(cabac_ctx, &cabac_state[MERGE_IDX_CTX_OFFSET])){
        int max_symbol = num_cand_min2;
        for(int k = 0; k < max_symbol; k++ ){
            if(!ovcabac_bypass_read(cabac_ctx)){
                max_symbol = k;
                break;
            }
        }
        gpm_ctx->merge_idx0 += max_symbol + 1;
    }
    if (num_cand_min2 > 0){
        if (ovcabac_ae_read(cabac_ctx, &cabac_state[MERGE_IDX_CTX_OFFSET])){
            int max_symbol = num_cand_min2 - 1;
            for(int k = 0; k < max_symbol; k++ ){
                if(!ovcabac_bypass_read(cabac_ctx)){
                    max_symbol = k;
                    break;
                }
            }
            gpm_ctx->merge_idx1 += max_symbol + 1;
        }
    }
    gpm_ctx->merge_idx1 += (gpm_ctx->merge_idx1 >= gpm_ctx->merge_idx0) ? 1 : 0;
}

static uint8_t
ovcabac_read_ae_reg_merge_flag(OVCABACCtx *const cabac_ctx, uint8_t skip_flag){
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    uint8_t offset = skip_flag ? 0 : 1;
    return ovcabac_ae_read(cabac_ctx, &cabac_state[REGULAR_MERGE_FLAG_CTX_OFFSET + offset]);
}

static uint8_t
ovcabac_read_ae_mmvd_flag(OVCABACCtx *const cabac_ctx){
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    return ovcabac_ae_read(cabac_ctx, &cabac_state[MMVD_FLAG_CTX_OFFSET]);
}

static uint8_t
ovcabac_read_ae_ciip_flag(OVCABACCtx *const cabac_ctx){
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    return ovcabac_ae_read(cabac_ctx, &cabac_state[CIIP_FLAG_CTX_OFFSET]);
}

static uint8_t
ovcabac_read_ae_bcw_flag(OVCABACCtx *const cabac_ctx, uint8_t is_ldc){
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    uint32_t idx = 0;
    uint32_t symbol = ovcabac_ae_read(cabac_ctx, &cabac_state[BCW_IDX_CTX_OFFSET]);

    int32_t numBcw = is_ldc ? 5 : 3;
    if(symbol == 1){
        uint32_t prefixNumBits = numBcw - 2;
        uint32_t step = 1;

        idx = 1;
        for(int ui = 0; ui < prefixNumBits; ++ui){
            if(!ovcabac_bypass_read(cabac_ctx)){
                break;
            }
            idx += step;
        }
    }
    int parsing_order[BCW_NUM] =  { BCW_DEFAULT, BCW_DEFAULT+1, BCW_DEFAULT-1,
                                        BCW_DEFAULT+2, BCW_DEFAULT-2};
    return (uint8_t)parsing_order[idx];
}

static uint8_t
ovcabac_read_ae_amvr_precision(OVCABACCtx *const cabac_ctx, uint8_t ibc_flag)
{
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    uint8_t value = 0;
    uint8_t prec_idx = 0;
    if (ibc_flag){
        value = 1;
    }
    else{
        value = ovcabac_ae_read(cabac_ctx, &cabac_state[IMV_FLAG_CTX_OFFSET]);
    }

    if( value ){
        if (!ibc_flag){
            value = ovcabac_ae_read(cabac_ctx, &cabac_state[IMV_FLAG_CTX_OFFSET + 4]);
            prec_idx = value ? 1 : 3;
        }
        if (value){
            value = ovcabac_ae_read(cabac_ctx, &cabac_state[IMV_FLAG_CTX_OFFSET + 1]);
            prec_idx = value + 1;
        }
    }

    int amvr_precision[4] = { MV_PRECISION_QUARTER, MV_PRECISION_INT, MV_PRECISION_4PEL, MV_PRECISION_HALF };

    return amvr_precision[prec_idx];
}

static uint8_t
ovcabac_read_ae_affine_amvr_precision(OVCABACCtx *const cabac_ctx, uint8_t ibc_flag)
{
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    uint8_t value = 0;
    uint8_t prec_idx = 0;

    value = ovcabac_ae_read(cabac_ctx, &cabac_state[IMV_FLAG_CTX_OFFSET + 2]);
    if( value ){
        value = ovcabac_ae_read(cabac_ctx, &cabac_state[IMV_FLAG_CTX_OFFSET + 3]);
        prec_idx = value + 1;
    }

    int aff_amvr_precision[4] = { MV_PRECISION_QUARTER, MV_PRECISION_SIXTEENTH, MV_PRECISION_INT };
    return aff_amvr_precision[prec_idx];
}


static uint8_t
ovcabac_read_ae_inter_dir(OVCABACCtx *const cabac_ctx,
                          int log2_pb_w, int log2_pb_h)
{
    uint8_t inter_dir = 0;
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    /*FIXME Note no Bi predicition on 4x4 8x4 4x8*/
    if (log2_pb_w + log2_pb_h > 5) {
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
ovcabac_read_ae_smvd_flag(OVCABACCtx *const cabac_ctx)
{
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    return ovcabac_ae_read(cabac_ctx, &cabac_state[SMVD_FLAG_CTX_OFFSET]);
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

static uint8_t
ovcabac_read_ae_ref_idx(OVCABACCtx *const cabac_ctx, uint8_t nb_active_ref)
{
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    uint8_t ref_idx = 0;
    if (ovcabac_ae_read(cabac_ctx, &cabac_state[REF_PIC_CTX_OFFSET])) {
        ref_idx += 1;
        if (nb_active_ref > 2 && ovcabac_ae_read(cabac_ctx, &cabac_state[REF_PIC_CTX_OFFSET + 1])) {
            ref_idx += 1;
            while (nb_active_ref <= ref_idx && ovcabac_bypass_read(cabac_ctx)) {
                ref_idx++;
            }
        }
    }

    return ref_idx;
}

static uint8_t
drv_intra_mode_c(VVCCU cu, uint8_t luma_mode)
{
    uint8_t cclm_flag = !!(cu.cu_flags & flg_cclm_flag);
    uint8_t mpm_flag_c = !!(cu.cu_flags & flg_mpm_flag_c);
    uint8_t mpm_idx = cu.cu_mode_idx_c;
    uint8_t cclm_idx = cu.cu_mode_idx_c;
    uint8_t intra_mode = derive_intra_mode_c(cclm_flag, mpm_flag_c, mpm_idx,
                                             luma_mode, cclm_idx);

    return intra_mode;
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

    /* FIXME move after TU is read so we can reconstruct with or without
     * transform trees
     */
    if (cu.cu_flags & 0x2) {
        if (ctu_dec->coding_unit == &coding_unit_intra_st || ctu_dec->coding_unit == &coding_unit_inter_st) {
            uint8_t luma_mode;
            luma_mode = drv_intra_cu(ctu_dec, part_ctx, x0, y0, log2_cb_w, log2_cb_h, cu);
            cu.cu_mode_idx = luma_mode;

            ctu_dec->intra_mode_c = drv_intra_mode_c(cu, luma_mode);
            vvc_intra_pred_chroma(&ctu_dec->rcn_ctx, &ctu_dec->rcn_ctx.ctu_buff, ctu_dec->intra_mode_c, x0 >> 1, y0 >> 1, log2_cb_w - 1, log2_cb_h - 1);
        } else {
            /* FIXME inter */
            if (ctu_dec->coding_unit == &coding_unit_intra) {
                uint8_t luma_mode = drv_intra_cu(ctu_dec, part_ctx, x0, y0, log2_cb_w, log2_cb_h, cu);
                cu.cu_mode_idx = luma_mode;
            } else {
                struct IntraDRVInfo *const i_info = &ctu_dec->drv_ctx.intra_info;
                uint8_t y_pu = y0 >> part_ctx->log2_min_cb_s;
                uint8_t x_pu = x0 >> part_ctx->log2_min_cb_s;
                uint8_t nb_pb_w = (1 << log2_cb_w) >> part_ctx->log2_min_cb_s;
                uint8_t nb_pb_h = (1 << log2_cb_h) >> part_ctx->log2_min_cb_s;
                uint8_t pu_shift = ctu_dec->part_ctx->log2_min_cb_s - 2;
                uint8_t luma_mode = i_info->luma_modes[(x_pu + ((y_pu + (nb_pb_h >> 1)) << 5) + (nb_pb_w >> 1))];

                ctu_dec->intra_mode_c = drv_intra_mode_c(cu, luma_mode);

                ctu_field_set_rect_bitfield(&ctu_dec->rcn_ctx.progress_field_c, x_pu << pu_shift,
                                            y_pu << pu_shift, nb_pb_w << pu_shift, nb_pb_h << pu_shift);

                vvc_intra_pred_chroma(&ctu_dec->rcn_ctx, &ctu_dec->rcn_ctx.ctu_buff, ctu_dec->intra_mode_c, x0, y0, log2_cb_w, log2_cb_h);
            }
        }
    }

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
        /* FIXME separate tree */
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
    if (ctu_dec->share != 1  &&
        ctu_dec->coding_tree != &dual_tree &&
        ctu_dec->coding_tree_implicit != &dual_tree_implicit) {
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
    uint8_t cu_type = OV_INTER;
    VVCCU cu = {0};

    cu_skip_flag = ovcabac_read_ae_cu_skip_flag(cabac_ctx, cu_type_abv,
                                                cu_type_lft);

    if (cu_skip_flag) {
        /* FIXME cu_skip_flag activation force merge_flag so we only need to read
           merge_idx */
        uint8_t merge_flag = 1;
        cu_type = ctu_dec->prediction_unit(ctu_dec, part_ctx, x0, y0, log2_cu_w, log2_cu_h, cu_skip_flag, merge_flag);

        if (cu_type == OV_AFFINE) {
            cu_type = OV_INTER_SKIP_AFFINE;
        } else {
            cu_type = OV_INTER_SKIP;
        }

        FLG_STORE(cu_skip_flag, cu.cu_flags);
        ctu_dec->intra_mode_c = 0;

    } else {
        uint8_t pred_mode_flag = ctu_dec->share == 1 ? 1 : 0;

        if (!ctu_dec->share) {
            pred_mode_flag = ovcabac_read_ae_pred_mode_flag(cabac_ctx, cu_type_abv,
                                                            cu_type_lft);
        }

        FLG_STORE(pred_mode_flag, cu.cu_flags);

        if (pred_mode_flag) {
            cu = coding_unit_intra_st(ctu_dec, part_ctx, x0, y0, log2_cu_w, log2_cu_h);

            if (cu.cu_flags & flg_mip_flag) {
                cu_type = OV_MIP;
            } else {
                cu_type = OV_INTRA;
            }

        } else {
            uint8_t merge_flag = ovcabac_read_ae_cu_merge_flag(cabac_ctx);
            cu_type = ctu_dec->prediction_unit(ctu_dec, part_ctx, x0, y0, log2_cu_w, log2_cu_h, cu_skip_flag, merge_flag);

            FLG_STORE(merge_flag, cu.cu_flags);
            ctu_dec->intra_mode_c = 0;
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
   VVCCU cu_c = {0};

   /* Force pred_mode_flag to 2 so we know cu was intra */

   cu = coding_unit_intra(ctu_dec, part_ctx, x0, y0, log2_cu_w, log2_cu_h);

   /* if not in separable tree */
   if (!ctu_dec->share) {
       cu_c = coding_unit_intra_c(ctu_dec, ctu_dec->part_ctx_c, x0 >> 1, y0 >> 1,
                                  log2_cu_w - 1, log2_cu_h - 1);

       cu.cu_flags |= cu_c.cu_flags;
       cu.cu_mode_idx_c = cu_c.cu_mode_idx_c;
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
    uint8_t mip_flag = 0;
    VVCCU cu = {0};

    cu.cu_flags |= flg_pred_mode_flag;

    if (ctu_dec->enabled_mip) {
        struct PartMap *part_map = &ctu_dec->part_map;
        uint8_t x_pu = x0 >> part_ctx->log2_min_cb_s;
        uint8_t y_pu = y0 >> part_ctx->log2_min_cb_s;
        uint8_t nb_pb_w = (1 << log2_cb_w) >> part_ctx->log2_min_cb_s;
        uint8_t nb_pb_h = (1 << log2_cb_h) >> part_ctx->log2_min_cb_s;
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

    return cu;
}

VVCCU
coding_unit_intra_c(OVCTUDec *const ctu_dec,
                    const OVPartInfo *const part_ctx,
                    uint8_t x0, uint8_t y0,
                    uint8_t log2_cb_w, uint8_t log2_cb_h)
{

    OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;
    VVCCU cu = {0};
    uint8_t cclm_flag = 0;

    /* Force intra pred mode in case of separarate or dual tree */
    cu.cu_flags = 2;

    /* FIXME CCLM luma partition constraints */
    if (ctu_dec->lm_chroma_enabled && (!ctu_dec->tmp_disable_cclm &&
        ctu_dec->enable_cclm == 1 || ctu_dec->coding_tree != &dual_tree)) {

        cclm_flag = ovcabac_read_ae_cclm_flag(cabac_ctx);
        FLG_STORE(cclm_flag, cu.cu_flags);

        if (cclm_flag) {
            uint8_t cclm_idx = ovcabac_read_ae_intra_lm_chroma(cabac_ctx);
            cu.cu_mode_idx_c = cclm_idx;
        }
    }

    if (!cclm_flag) {
        uint8_t mpm_flag_c = ovcabac_read_ae_intra_chroma_mpm_flag(cabac_ctx);
        FLG_STORE(mpm_flag_c, cu.cu_flags);
        if (mpm_flag_c) {
            uint8_t mpm_idx = ovcabac_read_ae_intra_chroma_mpm_idx(cabac_ctx);
            cu.cu_mode_idx_c = mpm_idx;
        }
    }

    return cu;
}

int
prediction_unit_inter_p(OVCTUDec *const ctu_dec,
                        const OVPartInfo *const part_ctx,
                        uint8_t x0, uint8_t y0,
                        uint8_t log2_pb_w, uint8_t log2_pb_h,
                        uint8_t skip_flag, uint8_t merge_flag)
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
    uint8_t ref_idx = 0;

    OVMV mv0;
    uint8_t apply_ciip = 0;
    if (merge_flag) {
        uint8_t mmvd_flag  = 0;
        uint8_t max_nb_cand = ctu_dec->max_num_merge_candidates;

        /* FIXME missing affine in P */
        uint8_t sps_ciip_flag = inter_ctx->ciip_flag;
        uint8_t ciip_flag = sps_ciip_flag && !skip_flag && (1 << log2_pb_w) < 128
                                                        && (1 << log2_pb_h) < 128
                                                        && 1 << (log2_pb_w + log2_pb_h) >= 64;
        uint8_t  reg_merge_flag = 1;
        if (ciip_flag) {
            reg_merge_flag = ovcabac_read_ae_reg_merge_flag(cabac_ctx, skip_flag);
        }

        if (reg_merge_flag){
            if (inter_ctx->mmvd_flag){
                mmvd_flag = ovcabac_read_ae_mmvd_flag(cabac_ctx);
            }
        } else {
            apply_ciip = 1;
        }

        uint8_t merge_idx;
        if (mmvd_flag){
            merge_idx = ovcabac_read_ae_mmvd_merge_idx(cabac_ctx, max_nb_cand);
            mv0 = drv_mmvd_merge_mvp(inter_ctx, mv_ctx0,
                            x_pu, y_pu, nb_pb_w, nb_pb_h,
                            merge_idx, max_nb_cand);
        } else {
            uint8_t merge_idx = ovcabac_read_ae_mvp_merge_idx(cabac_ctx, max_nb_cand);
            mv0 = drv_merge_mvp(inter_ctx, mv_ctx0,
                                x_pu, y_pu, nb_pb_w, nb_pb_h,
                                merge_idx, max_nb_cand);
        }

    } else {
        /*FIXME add ref_idx*/
        OVMV mvd = ovcabac_read_ae_mvd(cabac_ctx);
        if (inter_ctx->nb_active_ref0 > 1) {
            ref_idx = ovcabac_read_ae_ref_idx(cabac_ctx, inter_ctx->nb_active_ref0);
        }

        uint8_t mvp_idx = ovcabac_read_ae_mvp_flag(cabac_ctx);

        uint8_t prec_amvr = MV_PRECISION_QUARTER;
        mv0 = drv_mvp_mvd(inter_ctx, mv_ctx0, mvd, prec_amvr,
                          x_pu, y_pu, nb_pb_w, nb_pb_h,
                          mvp_idx, 1, ref_idx, ref_idx);
    }

    if(apply_ciip) {
        rcn_ciip(ctu_dec, x0, y0, log2_pb_w, log2_pb_h, mv0, ref_idx);
    } else {
        rcn_mcp(ctu_dec, ctu_dec->rcn_ctx.ctu_buff, x0, y0,
                log2_pb_w, log2_pb_h, mv0, 0, ref_idx);
    }

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

static inline uint8_t
check_bdof(uint8_t log2_pu_w, uint8_t log2_pu_h, uint8_t ciip_flag, uint8_t smvd_flag, uint8_t bcw_flag)
{
    uint8_t bdof_enabled = (log2_pu_h >= 3) && (log2_pu_w >= 3) && ((log2_pu_w + log2_pu_h) >= 7);
    bdof_enabled &= !ciip_flag;
    bdof_enabled &= !smvd_flag;
    bdof_enabled &= !bcw_flag;
    return bdof_enabled;
}

static inline uint8_t
check_bdof_ref(struct InterDRVCtx *const inter_ctx, uint8_t ref_idx0, uint8_t ref_idx1)
{
    return inter_ctx->tmvp_ctx.dist_ref_0[ref_idx0] == -inter_ctx->tmvp_ctx.dist_ref_1[ref_idx1];
}

int
prediction_unit_inter_b(OVCTUDec *const ctu_dec,
                        const OVPartInfo *const part_ctx,
                        uint8_t x0, uint8_t y0,
                        uint8_t log2_pb_w, uint8_t log2_pb_h,
                        uint8_t skip_flag, uint8_t merge_flag)
{
    OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;
    #if 1
    struct InterDRVCtx *const inter_ctx = &ctu_dec->drv_ctx.inter_ctx;
    struct IntraDRVInfo *const i_info = &ctu_dec->drv_ctx.intra_info;
    VVCMergeInfo mv_info;
    #endif
    uint8_t ref_idx0 = 0;
    uint8_t ref_idx1 = 0;
    uint8_t cu_type = OV_INTER;

#if 1
    uint8_t y_pu = y0 >> part_ctx->log2_min_cb_s;
    uint8_t x_pu = x0 >> part_ctx->log2_min_cb_s;
    uint8_t nb_pb_w = (1 << log2_pb_w) >> part_ctx->log2_min_cb_s;
    uint8_t nb_pb_h = (1 << log2_pb_h) >> part_ctx->log2_min_cb_s;
    uint8_t pu_shift = part_ctx->log2_min_cb_s - 2;
#endif

    //TODO: use real ibc flag
    uint8_t ibc_flag = 0;
    uint8_t smvd_flag = 0;
    inter_ctx->prec_amvr = MV_PRECISION_QUARTER;
    uint8_t bcw_idx = BCW_DEFAULT;

    if (merge_flag) {
        uint8_t mmvd_flag = 0;
        uint8_t max_nb_cand = ctu_dec->max_num_merge_candidates;

        uint8_t ciip_enabled = inter_ctx->ciip_flag && !skip_flag &&  log2_pb_w < 7
                                                                  &&  log2_pb_h < 7
                                                                  && (log2_pb_w + log2_pb_h) >= 6;

        uint8_t gpm_enabled  = inter_ctx->gpm_flag && inter_ctx->max_gpm_cand > 1
                                                   && log2_pb_w >= 3 && log2_pb_h >= 3
                                                   && log2_pb_w <= 6 && log2_pb_h <= 6
                                                   && log2_pb_w < 3 + log2_pb_h
                                                   && log2_pb_h < 3 + log2_pb_w;

        if (ctu_dec->affine_enabled && log2_pb_w >= 3 && log2_pb_h >= 3) {
            uint8_t cu_type_abv = ctu_dec->part_map.cu_mode_x[x_pu];
            uint8_t cu_type_lft = ctu_dec->part_map.cu_mode_y[y_pu];

            uint8_t lft_affine = cu_type_lft == OV_AFFINE || cu_type_lft == OV_INTER_SKIP_AFFINE;
            uint8_t abv_affine = cu_type_abv == OV_AFFINE || cu_type_abv == OV_INTER_SKIP_AFFINE;

            uint8_t sb_merge_flag = ovcabac_read_ae_sb_merge_flag(cabac_ctx, lft_affine, abv_affine);
            if (sb_merge_flag) {
                uint8_t nb_affine_merge_cand_min1 = ctu_dec->affine_nb_merge_cand - 1;
                uint8_t merge_idx = ovcabac_read_ae_affine_merge_idx(cabac_ctx, nb_affine_merge_cand_min1);
                cu_type = OV_AFFINE;

                /* TODO call affine drv merge and rcn functions */
                drv_affine_merge_mvp_b(inter_ctx, x0, y0, log2_pb_w, log2_pb_h,
                                       merge_idx);

                goto end;
            }
        }

        uint8_t  reg_merge_flag = 1;
        if (gpm_enabled || ciip_enabled) {
            reg_merge_flag = ovcabac_read_ae_reg_merge_flag(cabac_ctx, skip_flag);
            if (!reg_merge_flag) {
                uint8_t ciip_flag = ciip_enabled;
                if (gpm_enabled && ciip_enabled) {
                    ciip_flag = ovcabac_read_ae_ciip_flag(cabac_ctx);
                }

                if (ciip_flag) {
                    uint8_t merge_idx = ovcabac_read_ae_mvp_merge_idx(cabac_ctx, max_nb_cand);

                    mv_info = drv_merge_mvp_b(inter_ctx, x_pu, y_pu,
                                              nb_pb_w, nb_pb_h, merge_idx,
                                              max_nb_cand, log2_pb_w + log2_pb_h <= 5);

                    inter_ctx->prec_amvr = mv_info.inter_dir & 0x1 ? mv_info.mv0.prec_amvr
                                                                   : mv_info.mv1.prec_amvr;

                    ref_idx0 = mv_info.mv0.ref_idx;
                    ref_idx1 = mv_info.mv1.ref_idx;

                    mv_info.mv0.bcw_idx_plus1 = 0;
                    mv_info.mv1.bcw_idx_plus1 = 0;

                    rcn_ciip_b(ctu_dec, mv_info.mv0, mv_info.mv1, x0, y0,
                               log2_pb_w, log2_pb_h, mv_info.inter_dir, ref_idx0, ref_idx1);

                    goto end;

                } else {
                    int max_num_gpm_cand = inter_ctx->max_gpm_cand;

                    ovcabac_read_ae_gpm_merge_idx(cabac_ctx, &inter_ctx->gpm_ctx, max_num_gpm_cand);

                    drv_gpm_merge_mvp_b(inter_ctx, x_pu, y_pu, nb_pb_w, nb_pb_h, max_nb_cand,
                                        log2_pb_w + log2_pb_h <= 5);

                    rcn_gpm_b(ctu_dec, &inter_ctx->gpm_ctx, x0, y0, log2_pb_w, log2_pb_h);

                    goto end;
                }
            }
        }

        if (reg_merge_flag) {
            if (inter_ctx->mmvd_flag) {
                mmvd_flag = ovcabac_read_ae_mmvd_flag(cabac_ctx);
                if (mmvd_flag) {
                    uint8_t merge_idx = ovcabac_read_ae_mmvd_merge_idx(cabac_ctx, max_nb_cand);

                    mv_info = drv_mmvd_merge_mvp_b(inter_ctx, x_pu, y_pu,
                                                   nb_pb_w, nb_pb_h, ctu_dec->cur_poc, merge_idx,
                                                   max_nb_cand, log2_pb_w + log2_pb_h <= 5);
                }
            }
        }

        if (!mmvd_flag) {
            uint8_t merge_idx = ovcabac_read_ae_mvp_merge_idx(cabac_ctx, max_nb_cand);
            mv_info = drv_merge_mvp_b(inter_ctx, x_pu, y_pu,
                                  nb_pb_w, nb_pb_h, merge_idx,
                                  max_nb_cand, log2_pb_w + log2_pb_h <= 5);
        }

        inter_ctx->prec_amvr = mv_info.inter_dir & 0x1 ? mv_info.mv0.prec_amvr
                                                       : mv_info.mv1.prec_amvr;
        ref_idx0 = mv_info.mv0.ref_idx;
        ref_idx1 = mv_info.mv1.ref_idx;

    } else {
        OVMV mvd0, mvd1 = {0};

        uint8_t bcw_idx = BCW_DEFAULT;

        uint8_t mvp_idx0 = 0;
        uint8_t mvp_idx1 = 0;

        uint8_t inter_dir = ovcabac_read_ae_inter_dir(cabac_ctx, log2_pb_w, log2_pb_h);

        uint8_t affine_flag = 0;
        if (ctu_dec->affine_enabled && log2_pb_w > 3 && log2_pb_h > 3) {
            uint8_t cu_type_abv = ctu_dec->part_map.cu_mode_x[x_pu];
            uint8_t cu_type_lft = ctu_dec->part_map.cu_mode_y[y_pu];

            uint8_t lft_affine = cu_type_lft == OV_AFFINE || cu_type_lft == OV_INTER_SKIP_AFFINE;
            uint8_t abv_affine = cu_type_abv == OV_AFFINE || cu_type_abv == OV_INTER_SKIP_AFFINE;

            affine_flag = ovcabac_read_ae_cu_affine_flag(cabac_ctx, lft_affine, abv_affine);
            if (affine_flag) {
                uint8_t six_affine_type = !!(ctu_dec->affine_status & 0x2);
                uint8_t affine_type = !six_affine_type ? 0 : ovcabac_read_ae_cu_affine_type(cabac_ctx);
                struct AffineControlInfo cp_mvd0;
                struct AffineControlInfo cp_mvd1;

                cu_type = OV_AFFINE;
                int32_t mvd_not_zero = 0;
                if (inter_dir & 0x1) {

                    if (inter_ctx->nb_active_ref0 > 1) {
                        ref_idx0 = ovcabac_read_ae_ref_idx(cabac_ctx, inter_ctx->nb_active_ref0);
                    }

                    cp_mvd0.lt = ovcabac_read_ae_mvd(cabac_ctx);
                    cp_mvd0.rt = ovcabac_read_ae_mvd(cabac_ctx);

                    if (affine_type) {
                        cp_mvd0.lb = ovcabac_read_ae_mvd(cabac_ctx);
                        mvd_not_zero |= (cp_mvd0.lb.x | cp_mvd0.lb.y);
                    }

                    mvp_idx0 = ovcabac_read_ae_mvp_flag(cabac_ctx);

                    mvd_not_zero |= (cp_mvd0.lt.x | cp_mvd0.lt.y);
                    mvd_not_zero |= (cp_mvd0.rt.x | cp_mvd0.rt.y);
                }

                if (inter_dir & 0x2) {
                    if (inter_ctx->nb_active_ref1 > 1) {
                        ref_idx1 = ovcabac_read_ae_ref_idx(cabac_ctx, inter_ctx->nb_active_ref1);
                    }

                    if (inter_dir & 0x1 && inter_ctx->mvd1_zero_flag) {
                        memset(&cp_mvd1, 0, sizeof(cp_mvd1));
                    } else {
                        cp_mvd1.lt = ovcabac_read_ae_mvd(cabac_ctx);
                        cp_mvd1.rt = ovcabac_read_ae_mvd(cabac_ctx);
                        if (affine_type) {
                            cp_mvd1.lb = ovcabac_read_ae_mvd(cabac_ctx);
                            mvd_not_zero |= (cp_mvd1.lb.x | cp_mvd1.lb.y);
                        }
                    }

                    mvp_idx1 = ovcabac_read_ae_mvp_flag(cabac_ctx);

                    mvd_not_zero |= (cp_mvd1.lt.x | cp_mvd1.lt.y);
                    mvd_not_zero |= (cp_mvd1.rt.x | cp_mvd1.rt.y);
                }

                if (inter_ctx->affine_amvr_flag && mvd_not_zero && affine_flag && !skip_flag) {
                    inter_ctx->prec_amvr = ovcabac_read_ae_affine_amvr_precision(cabac_ctx, ibc_flag);
                }

                if (inter_ctx->bcw_flag && !ibc_flag && (1 << (log2_pb_h + log2_pb_w) >= BCW_SIZE_CONSTRAINT)
                                                     && !skip_flag && inter_dir == 3) {
                    bcw_idx = ovcabac_read_ae_bcw_flag(cabac_ctx, inter_ctx->tmvp_ctx.ldc);
                }

                /* TODO call affine drv MVP and rcn functions */
                drv_affine_mvp_b(inter_ctx, x0, y0, log2_pb_w, log2_pb_h,
                                 &cp_mvd0, &cp_mvd1, mvp_idx0, mvp_idx1, bcw_idx,
                                 inter_dir, ref_idx0, ref_idx1,
                                 affine_type);

                goto end;
            }
        }


        if (inter_dir == 3 && !affine_flag && inter_ctx->bi_dir_pred_flag) {
            smvd_flag = ovcabac_read_ae_smvd_flag(cabac_ctx);
        }

        int32_t mvd_not_zero = 0;
        if (inter_dir & 0x1) {
            if (smvd_flag) {
                ref_idx0 = inter_ctx->ref_smvd_idx0;
            } else if (inter_ctx->nb_active_ref0 > 1) {
                ref_idx0 = ovcabac_read_ae_ref_idx(cabac_ctx, inter_ctx->nb_active_ref0);
            }

            mvd0 = ovcabac_read_ae_mvd(cabac_ctx);
            mvp_idx0 = ovcabac_read_ae_mvp_flag(cabac_ctx);

            mvd_not_zero |= (mvd0.x | mvd0.y);
        }

        if (inter_dir & 0x2) {
            if (!smvd_flag) {
                if (inter_ctx->nb_active_ref1 > 1) {
                    ref_idx1 = ovcabac_read_ae_ref_idx(cabac_ctx, inter_ctx->nb_active_ref1);
                }

                if (inter_dir & 0x1 && inter_ctx->mvd1_zero_flag) {
                } else {
                    mvd1 = ovcabac_read_ae_mvd(cabac_ctx);
                    mvd_not_zero |= (mvd1.x | mvd1.y);
                }
            }
            mvp_idx1 = ovcabac_read_ae_mvp_flag(cabac_ctx);
        }

        if (inter_ctx->amvr_flag && mvd_not_zero && !affine_flag && !skip_flag) {
            inter_ctx->prec_amvr = ovcabac_read_ae_amvr_precision(cabac_ctx, ibc_flag);
        }

        if (inter_ctx->bcw_flag && !ibc_flag && (1<<(log2_pb_h+log2_pb_w) >= BCW_SIZE_CONSTRAINT)
                && !skip_flag && inter_dir == 3) {
            bcw_idx = ovcabac_read_ae_bcw_flag( cabac_ctx, inter_ctx->tmvp_ctx.ldc);
        }

        if (smvd_flag) {
            mvd1.x       = -mvd0.x;
            mvd1.y       = -mvd0.y;
            mvd1.ref_idx = inter_ctx->ref_smvd_idx1;
            ref_idx1     = inter_ctx->ref_smvd_idx1;
        }

        mv_info = drv_mvp_b(inter_ctx, x_pu, y_pu, nb_pb_w, nb_pb_h,
                            mvd0, mvd1, inter_ctx->prec_amvr, mvp_idx0, mvp_idx1, bcw_idx,
                            inter_dir, ref_idx0, ref_idx1, log2_pb_w + log2_pb_h <= 5);
    }

    /* FIXME move all this to derivation
     * => nothing to read passed this line
     */
    {
        uint8_t bdof_enable = 0;
        uint8_t dmvr_enable = 0;
        if (ctu_dec->bdof_enabled && mv_info.inter_dir == 0x3) {
            /* Note ciip_flag is zero in this function */
            uint8_t ciip_flag = 0;
            uint8_t bcw_flag = (mv_info.mv0.bcw_idx_plus1 != 0 && mv_info.mv0.bcw_idx_plus1 != 3);
            bdof_enable = check_bdof(log2_pb_w, log2_pb_h, ciip_flag, bcw_flag, smvd_flag);

            bdof_enable = bdof_enable && check_bdof_ref(inter_ctx, ref_idx0, ref_idx1);
        }

        /* FIXME DMVR only enable when merge_flag is set
         * bcw_idx == default
         *
         */
        if (merge_flag && ctu_dec->dmvr_enabled && mv_info.inter_dir == 0x3) {
            dmvr_enable = check_bdof(log2_pb_w, log2_pb_h, 0, 0, 0);

            dmvr_enable = dmvr_enable && check_bdof_ref(inter_ctx, ref_idx0, ref_idx1);
        }

        if (!bdof_enable && !dmvr_enable) {
            rcn_mcp_b(ctu_dec, ctu_dec->rcn_ctx.ctu_buff, inter_ctx, part_ctx,
                      mv_info.mv0, mv_info.mv1, x0, y0,
                      log2_pb_w, log2_pb_h, mv_info.inter_dir, ref_idx0, ref_idx1);
        } else {
            uint8_t log2_w = OVMIN(log2_pb_w, 4);
            uint8_t log2_h = OVMIN(log2_pb_h, 4);
            uint8_t nb_sb_w = (1 << log2_pb_w) >> log2_w;
            uint8_t nb_sb_h = (1 << log2_pb_h) >> log2_h;
            int i, j;

            uint8_t disable_bdof = 0;
            OVMV mv0 = mv_info.mv0;
            OVMV mv1 = mv_info.mv1;
            for (i = 0; i < nb_sb_h; ++i) {
                for (j = 0; j < nb_sb_w; ++j) {
                    if (dmvr_enable) {
                        OVMV *tmvp_mv0 = inter_ctx->tmvp_mv[0].mvs;
                        OVMV *tmvp_mv1 = inter_ctx->tmvp_mv[1].mvs;
                        disable_bdof = rcn_dmvr_mv_refine(ctu_dec, ctu_dec->rcn_ctx.ctu_buff,
                                                          x0 + j * 16, y0 + i * 16,
                                                          log2_w, log2_h,
                                                          &mv0, &mv1,
                                                          ref_idx0, ref_idx1, bdof_enable);

                        /* FIXME temporary hack to override MVs on 8x8 grid for TMVP */
                        /* FIXME find an alternative in case x0 % 8  is != 0 */
                        tmvp_mv0[((x0 + 7 + j * 16) >> 3) + ((y0 + 7 + i * 16) >> 3) * 16] = mv0;
                        tmvp_mv1[((x0 + 7 + j * 16) >> 3) + ((y0 + 7 + i * 16) >> 3) * 16] = mv1;

                        if (log2_w > 3) {
                            tmvp_mv0[((x0 + 7 + j * 16) >> 3) + ((y0 + 7 + i * 16) >> 3) * 16 + 1] = mv0;
                            tmvp_mv1[((x0 + 7 + j * 16) >> 3) + ((y0 + 7 + i * 16) >> 3) * 16 + 1] = mv1;
                        }

                        if (log2_h > 3) {
                            tmvp_mv0[((x0 + 7 + j * 16) >> 3) + ((y0 + 7 + i * 16) >> 3) * 16 + 16] = mv0;
                            tmvp_mv1[((x0 + 7 + j * 16) >> 3) + ((y0 + 7 + i * 16) >> 3) * 16 + 16] = mv1;

                            if (log2_w > 3) {
                                tmvp_mv0[((x0 + 7 + j * 16) >> 3) + ((y0 + 7 + i * 16) >> 3) * 16 + 16+ 1] = mv0;
                                tmvp_mv1[((x0 + 7 + j * 16) >> 3) + ((y0 + 7 + i * 16) >> 3) * 16 + 16+ 1] = mv1;
                            }
                        }
                    }

                    if (bdof_enable) {
                        if (!dmvr_enable) {
                            rcn_bdof_mcp_l(ctu_dec, ctu_dec->rcn_ctx.ctu_buff,
                                           x0 + j * 16, y0 + i * 16, log2_w, log2_h,
                                           mv0, mv1, ref_idx0, ref_idx1);
                            rcn_mcp_b_c(ctu_dec, ctu_dec->rcn_ctx.ctu_buff, inter_ctx, part_ctx,
                                        mv0, mv1, x0 + j * 16, y0 + i * 16,
                                        log2_w, log2_h, mv_info.inter_dir, ref_idx0, ref_idx1);
                        }
                    } else {
                        if (!dmvr_enable) {
                            rcn_mcp_b(ctu_dec, ctu_dec->rcn_ctx.ctu_buff, inter_ctx, part_ctx,
                                      mv0, mv1, x0 + j * 16, y0 + i * 16,
                                      log2_w, log2_h, mv_info.inter_dir, ref_idx0, ref_idx1);
                        }
                    }
                    mv0 = mv_info.mv0;
                    mv1 = mv_info.mv1;
                }
            }
            #if 0
            if (!disable_bdof) {
                rcn_mcp_b_c(ctu_dec, ctu_dec->rcn_ctx.ctu_buff, inter_ctx, part_ctx,
                            mv_info.mv0, mv_info.mv1, x0, y0,
                            log2_pb_w, log2_pb_h, mv_info.inter_dir, ref_idx0, ref_idx1);
            }
            #endif

        }
    }

end:

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

    return cu_type;
}
