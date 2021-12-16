#include "ovutils.h"
#include "rcn_dequant.h"
#include "ctudec.h"
#include "bitdepth.h"

#define IQUANT_SHIFT 6
#define MAX_LOG2_TR_RANGE 15

static const int inverse_quant_scale_lut[2][6] =
{
    { 40, 45, 51, 57, 64,  72},
    { 57, 64, 72, 80, 90, 102}
};

#if BITDEPTH == 10
void
derive_dequant_ctx(OVCTUDec *const ctudec, const struct VVCQPCTX *const qp_ctx,
                  int cu_qp_delta)
{
    int qp_bd_offset = qp_ctx->qp_bd_offset;
    int base_qp = (qp_ctx->current_qp + cu_qp_delta + 64) & 63;
    ctudec->dequant_luma.qp = (base_qp + qp_bd_offset);

    ctudec->dequant_cb.qp = qp_ctx->chroma_qp_map_cb[(base_qp + 64) & 63] + qp_ctx->cb_offset + qp_bd_offset;
    ctudec->dequant_cr.qp = qp_ctx->chroma_qp_map_cr[(base_qp + 64) & 63] + qp_ctx->cr_offset + qp_bd_offset;
    ctudec->dequant_joint_cb_cr.qp = qp_ctx->chroma_qp_map_jcbcr[(base_qp + 64) & 63] + qp_ctx->jcbcr_offset + qp_bd_offset;
    ctudec->dequant_luma_skip.qp = OVMAX(ctudec->dequant_luma.qp, qp_ctx->min_qp_prime_ts);
    ctudec->dequant_cb_skip.qp   = OVMAX(ctudec->dequant_cb.qp, qp_ctx->min_qp_prime_ts);
    ctudec->dequant_cr_skip.qp   = OVMAX(ctudec->dequant_cr.qp, qp_ctx->min_qp_prime_ts);
    ctudec->dequant_jcbcr_skip.qp = OVMAX(ctudec->dequant_joint_cb_cr.qp, qp_ctx->min_qp_prime_ts);

    ctudec->qp_ctx.current_qp = base_qp;
}
#endif

static void
dequant_sb_neg(int16_t *const sb_coeffs, int scale, int shift)
{
    for( int i = 0; i < 16 ; i++ ){
        sb_coeffs[i] = ov_clip_intp2((int32_t)sb_coeffs[i] * (scale << shift),
                MAX_LOG2_TR_RANGE + 1);
    }
}

static void
dequant_sb(int16_t *const sb_coeffs, int scale, int shift)
{
    int add = (1 << shift) >> 1;
    for( int i = 0; i < 16 ; i++ ){
        sb_coeffs[i] = ov_clip_intp2((int32_t)(sb_coeffs[i] * scale + add) >> shift,
                MAX_LOG2_TR_RANGE + 1);
    }
}


static struct IQScale
derive_dequant_sdh(int qp, uint8_t log2_tb_w, uint8_t log2_tb_h)
{
    const uint8_t log2_tb_s = log2_tb_w + log2_tb_h;
    struct IQScale dequant_params;
    int shift = IQUANT_SHIFT - (MAX_LOG2_TR_RANGE - BITDEPTH)
        - (qp / 6) + (log2_tb_s >> 1) + (log2_tb_s & 1);

    int scale = inverse_quant_scale_lut[log2_tb_s & 1][qp % 6];

    if (shift >= 0){
        dequant_params.shift = shift;
        dequant_params.scale = scale;
        dequant_params.dequant_sb = &dequant_sb;
    } else {
        dequant_params.shift = -shift;
        dequant_params.scale = scale;
        dequant_params.dequant_sb = &dequant_sb_neg;
    }

    return dequant_params;
}

static struct IQScale
derive_dequant_dpq(int qp, uint8_t log2_tb_w, uint8_t log2_tb_h)
{
    const uint8_t log2_tb_s = log2_tb_w + log2_tb_h;
    struct IQScale dequant_params;

    int shift = IQUANT_SHIFT + 1 - (MAX_LOG2_TR_RANGE - BITDEPTH)
        - ((qp + 1) / 6) + (log2_tb_s >> 1) + (log2_tb_s & 1);

    int scale  = inverse_quant_scale_lut[log2_tb_s & 1][(qp + 1) % 6];

    if (shift >= 0){
        dequant_params.shift = shift;
        dequant_params.scale = scale;
        dequant_params.dequant_sb = &dequant_sb;
    } else {
        dequant_params.shift = -shift;
        dequant_params.scale = scale;
        dequant_params.dequant_sb = &dequant_sb_neg;
    }

    return dequant_params;
}

static struct IQScale
derive_dequant_ts(int qp, uint8_t log2_tb_w, uint8_t log2_tb_h)
{
    struct IQScale dequant_params;

    int shift = IQUANT_SHIFT - (qp / 6) ;
    int scale  = inverse_quant_scale_lut[0][qp % 6];

    if (shift >= 0){
        dequant_params.shift = shift;
        dequant_params.scale = scale;
        dequant_params.dequant_sb = &dequant_sb;
    } else {
        dequant_params.shift = -shift;
        dequant_params.scale = scale;
        dequant_params.dequant_sb = &dequant_sb_neg;
    }

    return dequant_params;
}

void
BD_DECL(rcn_init_dequant)(struct RCNFunctions *rcn_funcs)
{
     rcn_funcs->tmp.derive_dequant_sdh = &derive_dequant_sdh;
     rcn_funcs->tmp.derive_dequant_ts = &derive_dequant_ts;
     rcn_funcs->tmp.derive_dequant_dpq = &derive_dequant_dpq;
}
