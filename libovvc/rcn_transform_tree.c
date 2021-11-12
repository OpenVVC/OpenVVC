#include <stdint.h>
#include <string.h>

#include "ovmem.h"
#include "ctudec.h"
#include "dbf_utils.h"
#include "drv.h"
#include "vcl.h"

#define BITDEPTH 10

struct TBInfo {
   uint16_t last_pos;
   uint64_t sig_sb_map;
};

struct TUInfo {
   uint8_t is_sbt;
   uint8_t cbf_mask;
   uint16_t pos_offset;
   uint8_t tr_skip_mask;
   uint8_t cu_mts_flag;
   uint8_t cu_mts_idx;
   uint8_t lfnst_flag;
   uint8_t lfnst_idx;
   struct TBInfo tb_info[3];
};

struct ISPTUInfo {
   uint8_t cbf_mask;
   uint8_t tr_skip_mask;
   uint8_t cu_mts_flag;
   uint8_t cu_mts_idx;
   uint8_t lfnst_flag;
   uint8_t lfnst_idx;
   struct TBInfo tb_info[4];
};

static void
derive_dequant_ctx(OVCTUDec *const ctudec, const VVCQPCTX *const qp_ctx,
                  int cu_qp_delta)
{
    /*FIXME avoid negative values especiallly in chroma_map derivation*/
    int qp_bd_offset = 6 * (BITDEPTH - 8);
    int base_qp = (qp_ctx->current_qp + cu_qp_delta + 64) & 63;
    ctudec->dequant_luma.qp = ((base_qp + qp_bd_offset) & 63);

    /*FIXME update transform skip ctx to VTM-10.0 */
    ctudec->dequant_luma_skip.qp = OVMAX(ctudec->dequant_luma.qp, qp_ctx->min_qp_prime_ts);
    ctudec->dequant_cb.qp = qp_ctx->chroma_qp_map_cb[(base_qp + qp_ctx->cb_offset + 64) & 63] + qp_bd_offset;
    ctudec->dequant_cr.qp = qp_ctx->chroma_qp_map_cr[(base_qp + qp_ctx->cr_offset + 64) & 63] + qp_bd_offset;
    ctudec->dequant_joint_cb_cr.qp = qp_ctx->chroma_qp_map_jcbcr[(base_qp + 64) & 63] + qp_ctx->jcbcr_offset + qp_bd_offset;
    ctudec->qp_ctx.current_qp = base_qp;
}

static void
rcn_residual(OVCTUDec *const ctudec,
             int16_t *const dst, int16_t *src,
             uint8_t x0, uint8_t y0,
             unsigned int log2_tb_w, unsigned int log2_tb_h,
             unsigned int lim_cg_w,
             uint8_t cu_mts_flag, uint8_t cu_mts_idx,
             uint8_t is_dc, uint8_t lfnst_flag, uint8_t is_mip, uint8_t lfnst_idx, uint8_t sbt)
{
    struct TRFunctions *TRFunc = &ctudec->rcn_ctx.rcn_funcs.tr;
    fill_bs_map(&ctudec->dbf_info.bs1_map, x0, y0, log2_tb_w, log2_tb_h);
    int shift_v = 6 + 1;
    int shift_h = (6 + 15 - 1) - BITDEPTH;

    DECLARE_ALIGNED(32, int16_t, tmp)[64*64];

    int tb_w = 1 << log2_tb_w;
    int tb_h = 1 << log2_tb_h;

    memset(tmp, 0, sizeof(int16_t) << (log2_tb_w + log2_tb_h));

    if (lfnst_flag) {
        /* FIXME separate lfnst mode derivation from lfnst reconstruction */
        process_lfnst_luma(ctudec, src, ctudec->lfnst_subblock, log2_tb_w, log2_tb_h, x0, y0,
                           lfnst_idx);
        lim_cg_w = 8;
        is_dc = 0;
    }

    if (!is_mip && !cu_mts_flag && ctudec->mts_implicit && (log2_tb_w <= 4 || log2_tb_h <= 4) && !lfnst_flag) {
        /*FIXME condition on size in the if could be removed ?*/
        enum DCTType tr_h_idx = log2_tb_w <= 4 ? DST_VII : DCT_II;
        enum DCTType tr_v_idx = log2_tb_h <= 4 ? DST_VII : DCT_II;

        /* FIXME use coefficient zeroing in MTS */
        TRFunc->func[tr_v_idx][log2_tb_h](src, tmp, tb_w, tb_w, tb_h, shift_v);
        TRFunc->func[tr_h_idx][log2_tb_w](tmp, dst, tb_h, tb_h, tb_w, shift_h);

    } else if (!cu_mts_flag) {

        if (is_dc) {

            TRFunc->dc(dst, log2_tb_w, log2_tb_h, src[0]);

        } else {
            int nb_row = OVMIN(lim_cg_w, 1 << log2_tb_w);
            int nb_col = OVMIN(lim_cg_w, 1 << log2_tb_h);

            TRFunc->func[DCT_II][log2_tb_h](src, tmp, tb_w, nb_row, nb_col, shift_v);
            TRFunc->func[DCT_II][log2_tb_w](tmp, dst, tb_h, tb_h, nb_row, shift_h);
        }
    } else {
        enum DCTType tr_h_idx = cu_mts_idx  & 1;
        enum DCTType tr_v_idx = cu_mts_idx >> 1;

        TRFunc->func[tr_v_idx][log2_tb_h](src, tmp, tb_w, tb_w, tb_h, shift_v);
        TRFunc->func[tr_h_idx][log2_tb_w](tmp, dst, tb_h, tb_h, tb_w, shift_h);
    }
}

static void
rcn_residual_c(OVCTUDec *const ctudec,
               int16_t *const dst, int16_t *src,
               uint8_t x0, uint8_t y0,
               uint8_t log2_tb_w, uint8_t log2_tb_h,
               uint16_t last_pos,
               uint8_t lfnst_flag, uint8_t lfnst_idx)
{
    struct TRFunctions *TRFunc = &ctudec->rcn_ctx.rcn_funcs.tr;

    const int shift_v = 6 + 1;
    const int shift_h = (6 + 15 - 1) - BITDEPTH;

    DECLARE_ALIGNED(32, int16_t, tmp)[32*32];

    int tb_w = 1 << log2_tb_w;
    int tb_h = 1 << log2_tb_h;

    memset(tmp, 0, sizeof(int16_t) << (log2_tb_w + log2_tb_h));

    if (lfnst_flag && log2_tb_w > 1 && log2_tb_h > 1) {
        /* FIXME separate lfnst mode derivation from lfnst reconstruction */
        int16_t lfnst_sb[16];
        memcpy(lfnst_sb     , &src[0], sizeof(int16_t) * 4);
        memcpy(lfnst_sb +  4, &src[1 << log2_tb_w], sizeof(int16_t) * 4);
        memcpy(lfnst_sb +  8, &src[2 << log2_tb_w], sizeof(int16_t) * 4);
        memcpy(lfnst_sb + 12, &src[3 << log2_tb_w], sizeof(int16_t) * 4);

        process_lfnst(ctudec, src, lfnst_sb, log2_tb_w, log2_tb_h,
                      x0, y0, lfnst_idx);
    }

    if (!last_pos && !lfnst_flag) {

        TRFunc->dc(dst, log2_tb_w, log2_tb_h, src[0]);

    } else {
        int lim_sb_s = ((((last_pos >> 8)) >> 2) + (((last_pos & 0xFF))>> 2) + 1) << 2;
        if (lfnst_flag && log2_tb_w > 1 && log2_tb_h > 1) lim_sb_s = 8;
        int nb_row =  OVMIN(lim_sb_s, 1 << log2_tb_w);
        int nb_col =  OVMIN(lim_sb_s, 1 << log2_tb_h);
        /*FIXME might be transform SKIP */

        TRFunc->func[DCT_II][log2_tb_h](src, tmp, tb_w, nb_row, nb_col, shift_v);
        TRFunc->func[DCT_II][log2_tb_w](tmp, dst, tb_h, tb_h, nb_row, shift_h);
    }
}

void
rcn_res_c(OVCTUDec *const ctu_dec, const struct TUInfo *tu_info,
          uint8_t x0, uint8_t y0,
          uint8_t log2_tb_w, uint8_t log2_tb_h, uint8_t cbf_mask, uint8_t lfnst_flag)
{
    const struct RCNFunctions *const rcn_func = &ctu_dec->rcn_ctx.rcn_funcs;

    if (cbf_mask & 0x2) {
        uint16_t *const dst_cb = &ctu_dec->rcn_ctx.ctu_buff.cb[(x0) + (y0 * RCN_CTB_STRIDE)];
        int16_t scale  =  ctu_dec->lmcs_info.scale_c_flag ? ctu_dec->lmcs_info.lmcs_chroma_scale : 1<< 11;
        int16_t *const coeffs_cb = ctu_dec->residual_cb + tu_info->pos_offset;

        if (!(tu_info->tr_skip_mask & 0x2)) {
            const struct TBInfo *const tb_info_cb = &tu_info->tb_info[0];
            rcn_residual_c(ctu_dec, ctu_dec->transform_buff, coeffs_cb,
                           x0, y0, log2_tb_w, log2_tb_h,
                           tb_info_cb->last_pos, lfnst_flag, tu_info->lfnst_idx);
            if (log2_tb_w + log2_tb_h > 2) {
                rcn_func->ict.ict[log2_tb_w][0](ctu_dec->transform_buff, dst_cb, log2_tb_w, log2_tb_h, scale);
            } else {
                rcn_func->ict.add[log2_tb_w](ctu_dec->transform_buff, dst_cb, log2_tb_w, log2_tb_h, scale);
            }
        } else {
            if (log2_tb_w + log2_tb_h > 2) {
                rcn_func->ict.ict[log2_tb_w][0](coeffs_cb, dst_cb, log2_tb_w, log2_tb_h, scale);
            } else {
                rcn_func->ict.add[log2_tb_w](coeffs_cb, dst_cb, log2_tb_w, log2_tb_h, scale);
            }
        }


        fill_bs_map(&ctu_dec->dbf_info.bs1_map_cb, x0 << 1, y0 << 1, log2_tb_w + 1, log2_tb_h + 1);
        fill_ctb_bound_c(&ctu_dec->dbf_info, x0 << 1, y0 << 1, log2_tb_w + 1, log2_tb_h + 1);
    }

    if (cbf_mask & 0x1) {
        uint16_t *const dst_cr = &ctu_dec->rcn_ctx.ctu_buff.cr[(x0) + (y0 * RCN_CTB_STRIDE)];
        int16_t scale  =  ctu_dec->lmcs_info.scale_c_flag ? ctu_dec->lmcs_info.lmcs_chroma_scale : 1<< 11;
        int16_t *const coeffs_cr = ctu_dec->residual_cr + tu_info->pos_offset;

        if (!(tu_info->tr_skip_mask & 0x1)) {
            const struct TBInfo *const tb_info_cr = &tu_info->tb_info[1];
            rcn_residual_c(ctu_dec, ctu_dec->transform_buff, coeffs_cr,
                           x0, y0, log2_tb_w, log2_tb_h,
                           tb_info_cr->last_pos, lfnst_flag, tu_info->lfnst_idx);
            if (log2_tb_w + log2_tb_h > 2) {
                rcn_func->ict.ict[log2_tb_w][0](ctu_dec->transform_buff, dst_cr, log2_tb_w, log2_tb_h, scale);
            } else {
                rcn_func->ict.add[log2_tb_w](ctu_dec->transform_buff, dst_cr, log2_tb_w, log2_tb_h, scale);
            }
        } else {
            if (log2_tb_w + log2_tb_h > 2) {
                rcn_func->ict.ict[log2_tb_w][0](coeffs_cr, dst_cr, log2_tb_w, log2_tb_h, scale);
            } else {
                rcn_func->ict.add[log2_tb_w](coeffs_cr, dst_cr, log2_tb_w, log2_tb_h, scale);
            }
        }

        fill_bs_map(&ctu_dec->dbf_info.bs1_map_cr, x0 << 1, y0 << 1, log2_tb_w + 1, log2_tb_h + 1);
        fill_ctb_bound_c(&ctu_dec->dbf_info, x0 << 1, y0 << 1, log2_tb_w + 1, log2_tb_h + 1);
    }
}

void
rcn_jcbcr(OVCTUDec *const ctu_dec, const struct TUInfo *const tu_info,
          uint8_t x0, uint8_t y0, uint8_t log2_tb_w, uint8_t log2_tb_h,
          uint8_t cbf_mask, uint8_t lfnst_flag)
{
    const struct RCNFunctions *const rcn_func = &ctu_dec->rcn_ctx.rcn_funcs;
    uint16_t *const dst_cb = &ctu_dec->rcn_ctx.ctu_buff.cb[x0 + (y0 * RCN_CTB_STRIDE)];
    uint16_t *const dst_cr = &ctu_dec->rcn_ctx.ctu_buff.cr[x0 + (y0 * RCN_CTB_STRIDE)];
    if (!(tu_info->tr_skip_mask & 0x1)) {
        const struct TBInfo *const tb_info = &tu_info->tb_info[0];
        int16_t *const coeffs_jcbcr = ctu_dec->residual_cb + tu_info->pos_offset;
        rcn_residual_c(ctu_dec, ctu_dec->transform_buff, coeffs_jcbcr,
                       x0, y0, log2_tb_w, log2_tb_h,
                       tb_info->last_pos, lfnst_flag, tu_info->lfnst_idx);
    } else {
        int16_t *const coeffs_jcbcr = ctu_dec->residual_cb + tu_info->pos_offset;
        memcpy(ctu_dec->transform_buff, coeffs_jcbcr, sizeof(int16_t) << (log2_tb_w + log2_tb_h));
    }

    fill_bs_map(&ctu_dec->dbf_info.bs1_map_cb, x0 << 1, y0 << 1, log2_tb_w + 1, log2_tb_h + 1);
    fill_bs_map(&ctu_dec->dbf_info.bs1_map_cr, x0 << 1, y0 << 1, log2_tb_w + 1, log2_tb_h + 1);
    fill_ctb_bound_c(&ctu_dec->dbf_info, x0 << 1, y0 << 1, log2_tb_w + 1, log2_tb_h + 1);
    if ((cbf_mask & 0x3) == 0x3) {
        const uint8_t qp_bd_offset = 6 * (BITDEPTH - 8);
        uint8_t    qp = ctu_dec->dequant_joint_cb_cr.qp - qp_bd_offset;
        dbf_fill_qp_map(&ctu_dec->dbf_info.qp_map_cb, x0 << 1, y0 << 1, log2_tb_w + 1, log2_tb_h + 1, qp);
        dbf_fill_qp_map(&ctu_dec->dbf_info.qp_map_cr, x0 << 1, y0 << 1, log2_tb_w + 1, log2_tb_h + 1, qp);
    }
    /* FIXME set jcbcr dbf qp */

    /* FIXME better organisation based on cbf_mask */
    if (cbf_mask == 3) {
        int16_t scale  =  ctu_dec->lmcs_info.scale_c_flag ? ctu_dec->lmcs_info.lmcs_chroma_scale : 1<< 11;
        if (log2_tb_w + log2_tb_h == 2) scale = 1<<11;
        rcn_func->ict.ict[log2_tb_w][0](ctu_dec->transform_buff, dst_cb, log2_tb_w, log2_tb_h, scale);
        rcn_func->ict.ict[log2_tb_w][1](ctu_dec->transform_buff, dst_cr, log2_tb_w, log2_tb_h, scale);
    } else if (cbf_mask == 2) {
        int16_t scale  =  ctu_dec->lmcs_info.scale_c_flag ? ctu_dec->lmcs_info.lmcs_chroma_scale : 1<< 11;
        if (log2_tb_w + log2_tb_h == 2) scale = 1<<11;
        rcn_func->ict.ict[log2_tb_w][0](ctu_dec->transform_buff, dst_cb, log2_tb_w, log2_tb_h, scale);
        rcn_func->ict.ict[log2_tb_w][2](ctu_dec->transform_buff, dst_cr, log2_tb_w, log2_tb_h, scale);
    } else {
        int16_t scale  =  ctu_dec->lmcs_info.scale_c_flag ? ctu_dec->lmcs_info.lmcs_chroma_scale : 1<< 11;
        if (log2_tb_w + log2_tb_h == 2) scale = 1<<11;
        rcn_func->ict.ict[log2_tb_w][0](ctu_dec->transform_buff, dst_cr, log2_tb_w, log2_tb_h, scale);
        rcn_func->ict.ict[log2_tb_w][2](ctu_dec->transform_buff, dst_cb, log2_tb_w, log2_tb_h, scale);
    }

}

void
recon_isp_subtree_v(OVCTUDec *const ctudec,
                    unsigned int x0, unsigned int y0,
                    unsigned int log2_cb_w, unsigned int log2_cb_h,
                    uint8_t intra_mode,
                    const struct ISPTUInfo *const tu_info)
{
    #if 0
    struct OVDrvCtx *const pred_ctx = &ctudec->drv_ctx;
    #endif
    const struct TRFunctions *TRFunc = &ctudec->rcn_ctx.rcn_funcs.tr;
    const struct RCNFunctions *const rcn_func = &ctudec->rcn_ctx.rcn_funcs;
    uint8_t cbf_flags = tu_info->cbf_mask;
    uint8_t lfnst_flag = tu_info->lfnst_flag;

    int offset_x;
    int log2_pb_w = log2_cb_w - 2;
    int nb_pb;
    int pb_w;

    if (log2_cb_h < 4 && (log2_pb_w <= (4 - log2_cb_h))) {
        log2_pb_w = 4 - log2_cb_h;
    }
    nb_pb = (1 << log2_cb_w) >> log2_pb_w;
    pb_w =  1 << log2_pb_w;
    offset_x = 0;

    uint8_t type_h = ctudec->mts_enabled && log2_pb_w <= 4 && log2_pb_w > 1 ? DST_VII : DCT_II;
    uint8_t type_v = ctudec->mts_enabled && log2_cb_h <= 4 ? DST_VII : DCT_II;
    int i;

    for (i = 0; i < nb_pb; ++i) {
        uint8_t cbf = (cbf_flags >> (nb_pb - i - 1)) & 0x1;

        /* On vertical ISP the intra prediction can only be performed with a min width of 4
           thus requiring to perform prediction on 2 PBs for size 2xX and all PBs in the
           case of 1xX. Even when transforms are smaller
         */
         /*FIXME separate small cases */
        if (!(offset_x & 0x3)) {
            vvc_intra_pred_isp(ctudec, &ctudec->rcn_ctx.ctu_buff.y[0],
                               RCN_CTB_STRIDE, intra_mode, x0, y0,
                               log2_pb_w >= 2 ? log2_pb_w : 2, log2_cb_h,
                               log2_cb_w, log2_cb_h, offset_x, 0);
        }

        /* FIXME determine if DBF is applied on ISP subblocks */
        #if 1
        if (!(offset_x & 0x3) ) {
            dbf_fill_qp_map(&ctudec->dbf_info.qp_map_y, x0, y0,  log2_pb_w >= 2 ? log2_pb_w : 2, log2_cb_h, ctudec->qp_ctx.current_qp);
            fill_ctb_bound(&ctudec->dbf_info, x0 & ~0x3, y0, log2_pb_w >= 2 ? log2_pb_w : 2, log2_cb_h);
            fill_bs_map(&ctudec->dbf_info.bs2_map, x0 & ~0x3, y0, log2_pb_w >= 2 ? log2_pb_w : 2, log2_cb_h);
        }
        //fill_ctb_bound(&ctudec->dbf_info, x0 & ~0x3, y0, log2_pb_w, log2_cb_h);
        #endif

        if (cbf) {
            int16_t *coeffs_y = ctudec->residual_y + i * (1 << (log2_pb_w + log2_cb_h));

            //fill_bs_map(&ctudec->dbf_info.bs1_map, x0, y0, log2_pb_w, log2_cb_h);

            if (log2_pb_w) {
                int shift_v = 6 + 1;
                int shift_h = (6 + 15 - 1) - BITDEPTH;
                DECLARE_ALIGNED(32, int16_t, tmp)[64*64];
                int16_t *src = coeffs_y;
                int16_t *dst = ctudec->transform_buff;
                int cb_h = 1 << log2_cb_h;

                memset(tmp, 0, sizeof(int16_t) << (log2_pb_w + log2_cb_h));
                if (lfnst_flag) {
                    uint8_t lfnst_idx = tu_info->lfnst_idx;
                    int16_t lfnst_sb[16];
                    memcpy(lfnst_sb     , &coeffs_y[0], sizeof(int16_t) * 4);
                    memcpy(lfnst_sb +  4, &coeffs_y[1 << log2_pb_w], sizeof(int16_t) * 4);
                    memcpy(lfnst_sb +  8, &coeffs_y[2 << log2_pb_w], sizeof(int16_t) * 4);
                    memcpy(lfnst_sb + 12, &coeffs_y[3 << log2_pb_w], sizeof(int16_t) * 4);
                    process_lfnst_luma_isp(ctudec, coeffs_y, lfnst_sb, log2_pb_w, log2_cb_h,log2_cb_w, log2_cb_h, x0 -offset_x, y0,
                                       lfnst_idx);
                    /* lfnst forces IDCT II usage */
                    type_v = type_h = DCT_II;
                }

                TRFunc->func[type_v][OVMIN(log2_cb_h,6)](src, tmp, pb_w, pb_w, cb_h, shift_v);
                TRFunc->func[type_h][OVMIN(log2_pb_w,6)](tmp, dst, cb_h, cb_h, pb_w, shift_h);
            } else {
                int shift_h = (6 + 15 - 1) - BITDEPTH;
                int cb_h = 1 << log2_cb_h;
                DECLARE_ALIGNED(32, int16_t, tmp)[64];

                memset(tmp, 0, sizeof(int16_t) << (log2_pb_w + log2_cb_h));

                TRFunc->func[type_v][OVMIN(log2_cb_h,6)](coeffs_y, tmp, pb_w, pb_w, cb_h, shift_h + 1);

                memcpy(ctudec->transform_buff, tmp, sizeof(uint16_t) * (1 << log2_cb_h));
            }

            uint16_t *dst  = &ctudec->rcn_ctx.ctu_buff.y[x0 + y0 * RCN_CTB_STRIDE];
            int16_t *src  = ctudec->transform_buff;

            rcn_func->ict.add[log2_pb_w](src, dst, log2_pb_w, log2_cb_h, 0);
        }
        x0 += pb_w;
        offset_x += pb_w;
    }
}

void
recon_isp_subtree_h(OVCTUDec *const ctudec,
                    unsigned int x0, unsigned int y0,
                    unsigned int log2_cb_w, unsigned int log2_cb_h,
                    uint8_t intra_mode,
                    const struct ISPTUInfo *const tu_info)
{
    const struct TRFunctions *TRFunc = &ctudec->rcn_ctx.rcn_funcs.tr;
    const struct RCNFunctions *const rcn_func = &ctudec->rcn_ctx.rcn_funcs;

    int log2_pb_h = log2_cb_h - 2;
    int nb_pb;
    uint8_t cbf_flags = tu_info->cbf_mask;
    uint8_t lfnst_flag = tu_info->lfnst_flag;
    int pb_h, offset_y;

    // width < 16 imposes restrictions on split numbers
    if (log2_cb_w < 4 && (log2_pb_h <= (4 - log2_cb_w))) {
        log2_pb_h = 4 - log2_cb_w;
    }

    nb_pb = (1 << log2_cb_h) >> log2_pb_h;
    pb_h =  1 << log2_pb_h;
    offset_y = 0;

    uint8_t type_h = ctudec->mts_enabled && log2_cb_w <= 4 ? DST_VII : DCT_II;
    uint8_t type_v = ctudec->mts_enabled && log2_pb_h <= 4 && log2_pb_h > 1 ? DST_VII : DCT_II;
    int i;

    for (i = 0; i < nb_pb; ++i) {

        uint8_t cbf = (cbf_flags >> (nb_pb - i - 1)) & 0x1;

        /* FIXME determine if DBF is applied on ISP subblocks */
        #if 1
        if (!(offset_y & 0x3) ) {
            dbf_fill_qp_map(&ctudec->dbf_info.qp_map_y, x0, y0, log2_cb_w, log2_pb_h >=2 ? log2_pb_h : 2, ctudec->qp_ctx.current_qp);
            fill_ctb_bound(&ctudec->dbf_info, x0, y0, log2_cb_w, log2_pb_h >=2 ? log2_pb_h : 2);
            fill_bs_map(&ctudec->dbf_info.bs2_map, x0, y0, log2_cb_w, log2_pb_h >=2 ? log2_pb_h : 2);
        }
        //fill_ctb_bound(&ctudec->dbf_info, x0, y0, log2_cb_w, log2_pb_h);
        #endif

        vvc_intra_pred_isp(ctudec, &ctudec->rcn_ctx.ctu_buff.y[0],
                           RCN_CTB_STRIDE, intra_mode, x0, y0,
                           log2_cb_w, log2_pb_h, log2_cb_w, log2_cb_h, 0, offset_y);
        if (cbf) {
            int16_t *coeffs_y = ctudec->residual_y + i * (1 << (log2_cb_w + log2_pb_h));

            fill_bs_map(&ctudec->dbf_info.bs1_map, x0, y0, log2_cb_w, log2_pb_h);

            if (log2_pb_h) {
                DECLARE_ALIGNED(32, int16_t, tmp)[64*64];
                int shift_v = 6 + 1;
                int shift_h = (6 + 15 - 1) - BITDEPTH;
                int16_t *src = coeffs_y;
                int16_t *dst = ctudec->transform_buff;
                int cb_w = 1 << log2_cb_w;

                memset(tmp, 0, sizeof(int16_t) << (log2_cb_w + log2_pb_h));

                if (lfnst_flag) {
                    uint8_t lfnst_idx = tu_info->lfnst_idx;
                    /*FIXME avoid distinguishing isp case in lfnst */
                    int16_t lfnst_sb[16];
                    memcpy(lfnst_sb     , &coeffs_y[0], sizeof(int16_t) * 4);
                    memcpy(lfnst_sb +  4, &coeffs_y[1 << log2_cb_w], sizeof(int16_t) * 4);
                    memcpy(lfnst_sb +  8, &coeffs_y[2 << log2_cb_w], sizeof(int16_t) * 4);
                    memcpy(lfnst_sb + 12, &coeffs_y[3 << log2_cb_w], sizeof(int16_t) * 4);
                    process_lfnst_luma_isp(ctudec, coeffs_y, lfnst_sb, log2_cb_w, log2_pb_h,log2_cb_w, log2_cb_h, x0, y0-offset_y,
                                       lfnst_idx);

                    /* lfnst forces IDCT II usage */
                    type_v = type_h = DCT_II;
                }

                TRFunc->func[type_v][OVMIN(log2_pb_h,6)](src, tmp, cb_w, cb_w, pb_h, shift_v);
                TRFunc->func[type_h][OVMIN(log2_cb_w,6)](tmp, dst, pb_h, pb_h, cb_w, shift_h);
            } else {
                int shift_h = (6 + 15 - 1) - BITDEPTH;
                int cb_w = 1 << log2_cb_w;
                DECLARE_ALIGNED(32, int16_t, tmp)[64];

                memset(tmp, 0, sizeof(int16_t) << (log2_cb_w + log2_pb_h));

                TRFunc->func[type_h][OVMIN(log2_cb_w,6)](coeffs_y, tmp, pb_h, pb_h, cb_w, shift_h + 1);

                memcpy(ctudec->transform_buff, tmp, sizeof(uint16_t) * (1 << log2_cb_w));
            }

            uint16_t *dst  = &ctudec->rcn_ctx.ctu_buff.y[x0 + y0 * RCN_CTB_STRIDE];
            int16_t *src  = ctudec->transform_buff;

            rcn_func->ict.add[log2_cb_w](src, dst, log2_cb_w, log2_pb_h, 0);
        }
        y0 += pb_h;
        offset_y += pb_h;
    }
}

void
rcn_tu_st(OVCTUDec *const ctu_dec,
          uint8_t x0, uint8_t y0,
          uint8_t log2_tb_w, uint8_t log2_tb_h,
          uint8_t cu_flags, uint8_t cbf_mask,
          const struct TUInfo *const tu_info)
{
    uint8_t cbf_flag_l = cbf_mask & 0x10;
    uint8_t jcbcr_flag = cbf_mask & 0x8;
    uint8_t cbf_mask_c = cbf_mask & 0x3;

    if (cbf_flag_l) {
        const struct TBInfo *const tb_info = &tu_info->tb_info[2];
        const struct RCNFunctions *const rcn_func = &ctu_dec->rcn_ctx.rcn_funcs;

        if (!(tu_info->tr_skip_mask & 0x10)) {
            int lim_sb_s = ((((tb_info->last_pos >> 8)) >> 2) + (((tb_info->last_pos & 0xFF))>> 2) + 1) << 2;
            int16_t *const coeffs_y = ctu_dec->residual_y + tu_info->pos_offset;
            uint8_t is_mip = !!(cu_flags & flg_mip_flag);
            uint8_t is_intra = !!(cu_flags & 0x2);
            is_mip |= !is_intra;
            rcn_residual(ctu_dec, ctu_dec->transform_buff, coeffs_y, x0, y0, log2_tb_w, log2_tb_h,
                         lim_sb_s, tu_info->cu_mts_flag, tu_info->cu_mts_idx,
                         !tb_info->last_pos, tu_info->lfnst_flag, is_mip, tu_info->lfnst_idx, tu_info->is_sbt);

        } else {
            int16_t *const coeffs_y = ctu_dec->residual_y + tu_info->pos_offset;
            memcpy(ctu_dec->transform_buff, coeffs_y, sizeof(int16_t) << (log2_tb_w + log2_tb_h));
        }

        /* FIXME use transform add optimization */
        rcn_func->ict.add[log2_tb_w](ctu_dec->transform_buff, &ctu_dec->rcn_ctx.ctu_buff.y[x0 + y0 * RCN_CTB_STRIDE], log2_tb_w, log2_tb_h, 0);
        /* FIXME Avoid reprocessing CCLM from here by recontructing at the end of transform tree */
        if (ctu_dec->intra_mode_c >= 67 && ctu_dec->intra_mode_c < 70) {
            vvc_intra_pred_chroma(&ctu_dec->rcn_ctx, ctu_dec->intra_mode_c, x0 >> 1, y0 >> 1, log2_tb_w - 1, log2_tb_h - 1);
        }
        fill_bs_map(&ctu_dec->dbf_info.bs1_map, x0, y0, log2_tb_w, log2_tb_h);
        fill_ctb_bound(&ctu_dec->dbf_info, x0, y0, log2_tb_w, log2_tb_h);
    }

    if (jcbcr_flag) {

        rcn_jcbcr(ctu_dec, tu_info, x0 >> 1, y0 >> 1, log2_tb_w - 1, log2_tb_h - 1, cbf_mask_c, 0);

    } else if (cbf_mask_c) {

        rcn_res_c(ctu_dec, tu_info, x0 >> 1, y0 >> 1, log2_tb_w - 1, log2_tb_h - 1, cbf_mask_c, 0);

    }

    const uint8_t qp_bd_offset = 6 * (BITDEPTH - 8);
    derive_dequant_ctx(ctu_dec, &ctu_dec->qp_ctx, 0);
    struct DBFInfo *dbf_info = &ctu_dec->dbf_info;
    uint8_t qp_l  = ctu_dec->qp_ctx.current_qp;
    uint8_t qp_jc = ctu_dec->dequant_joint_cb_cr.qp - qp_bd_offset;
    uint8_t qp_cb = (cbf_mask & 0x3) == 0x3 && jcbcr_flag ? qp_jc : ctu_dec->dequant_cb.qp - qp_bd_offset;
    uint8_t qp_cr = (cbf_mask & 0x3) == 0x3 && jcbcr_flag ? qp_jc : ctu_dec->dequant_cr.qp - qp_bd_offset;

    dbf_fill_qp_map(&dbf_info->qp_map_y, x0, y0, log2_tb_w, log2_tb_h, qp_l);
    dbf_fill_qp_map(&dbf_info->qp_map_cb, x0, y0, log2_tb_w, log2_tb_h, qp_cb);
    dbf_fill_qp_map(&dbf_info->qp_map_cr, x0, y0, log2_tb_w, log2_tb_h, qp_cr);
}


static void
rcn_tu_l(OVCTUDec *const ctu_dec,
         uint8_t x0, uint8_t y0,
         uint8_t log2_tb_w, uint8_t log2_tb_h,
         uint8_t cu_flags, uint8_t cbf_mask,
         const struct TUInfo *const tu_info)
{
    const struct TBInfo *const tb_info = &tu_info->tb_info[2];
    if (cbf_mask) {
        const struct RCNFunctions *const rcn_func = &ctu_dec->rcn_ctx.rcn_funcs;

        if (!(tu_info->tr_skip_mask & 0x10)) {
            int lim_sb_s = ((((tb_info->last_pos >> 8)) >> 2) + (((tb_info->last_pos & 0xFF))>> 2) + 1) << 2;
            int16_t *const coeffs_y = ctu_dec->residual_y + tu_info->pos_offset;
            uint8_t is_mip = !!(cu_flags & flg_mip_flag);
            uint8_t is_intra = !!(cu_flags & 0x2);
            is_mip |= !is_intra;
            rcn_residual(ctu_dec, ctu_dec->transform_buff, coeffs_y, x0, y0, log2_tb_w, log2_tb_h,
                         lim_sb_s, tu_info->cu_mts_flag, tu_info->cu_mts_idx,
                         !tb_info->last_pos, tu_info->lfnst_flag, is_mip, tu_info->lfnst_idx, tu_info->is_sbt);

        } else {
            int16_t *const coeffs_y = ctu_dec->residual_y + tu_info->pos_offset;
            memcpy(ctu_dec->transform_buff, coeffs_y, sizeof(int16_t) << (log2_tb_w + log2_tb_h));
        }

        fill_bs_map(&ctu_dec->dbf_info.bs1_map, x0, y0, log2_tb_w, log2_tb_h);
        fill_ctb_bound(&ctu_dec->dbf_info, x0, y0, log2_tb_w, log2_tb_h);

        /* FIXME use transform add optimization */
        rcn_func->ict.add[log2_tb_w](ctu_dec->transform_buff, &ctu_dec->rcn_ctx.ctu_buff.y[x0 + y0 * RCN_CTB_STRIDE], log2_tb_w, log2_tb_h, 0);
    }
}

void
rcn_tu_c(OVCTUDec *const ctu_dec, uint8_t x0, uint8_t y0,
         uint8_t log2_tb_w, uint8_t log2_tb_h,
         uint8_t cu_flags, uint8_t cbf_mask,
         const struct TUInfo *const tu_info)
{
    uint8_t jcbcr_flag = cbf_mask & 0x8;
    uint8_t cbf_mask_c = cbf_mask & 0x3;

    if (jcbcr_flag) {

        rcn_jcbcr(ctu_dec, tu_info, x0, y0, log2_tb_w, log2_tb_h, cbf_mask_c, tu_info->lfnst_flag);

    } else if (cbf_mask_c) {

        rcn_res_c(ctu_dec, tu_info, x0, y0, log2_tb_w, log2_tb_h, cbf_mask_c, tu_info->lfnst_flag);

    }
}

static void
rcn_res_wrap(OVCTUDec *const ctu_dec, uint8_t x0, uint8_t y0,
             uint8_t log2_tb_w, uint8_t log2_tb_h, uint8_t cu_flags,
             const struct TUInfo *const tu_info)
{
    uint8_t cbf_mask = tu_info->cbf_mask;
    if (ctu_dec->transform_unit == &transform_unit_st && cbf_mask) {
        rcn_tu_st(ctu_dec, x0, y0, log2_tb_w, log2_tb_h, cu_flags, cbf_mask, tu_info);
    } else if (ctu_dec->transform_unit == &transform_unit_l && cbf_mask) {
        rcn_tu_l(ctu_dec, x0, y0, log2_tb_w, log2_tb_h, cu_flags, cbf_mask, tu_info);
    } else if (ctu_dec->transform_unit == &transform_unit_c && cbf_mask) {
        rcn_tu_c(ctu_dec, x0, y0, log2_tb_w, log2_tb_h, cu_flags, cbf_mask, tu_info);
    }
}

void
rcn_transform_tree(OVCTUDec *const ctu_dec, uint8_t x0, uint8_t y0,
                   uint8_t log2_tb_w, uint8_t log2_tb_h, uint8_t log2_max_tb_s,
                   uint8_t cu_flags, const struct TUInfo *const tu_info)
{
    uint8_t split_v = log2_tb_w > log2_max_tb_s;
    uint8_t split_h = log2_tb_h > log2_max_tb_s;

    if (split_v || split_h) {
        unsigned int tb_w1 = ((1 << log2_tb_w) >> split_v);
        unsigned int tb_h1 = ((1 << log2_tb_h) >> split_h);

        unsigned int log2_tb_w1 = log2_tb_w - split_v;
        unsigned int log2_tb_h1 = log2_tb_h - split_h;

        rcn_transform_tree(ctu_dec, x0, y0,
                           log2_tb_w1, log2_tb_h1,
                           log2_max_tb_s, cu_flags, &tu_info[0]);
        if (split_v) {
            rcn_transform_tree(ctu_dec, x0 + tb_w1, y0,
                               log2_tb_w1, log2_tb_h1,
                               log2_max_tb_s, cu_flags, &tu_info[1]);
        }

        if (split_h) {
            rcn_transform_tree(ctu_dec, x0, y0 + tb_h1,
                               log2_tb_w1, log2_tb_h1,
                               log2_max_tb_s, cu_flags, &tu_info[2]);
        }

        if (split_h && split_v) {
            rcn_transform_tree(ctu_dec, x0 + tb_w1, y0 + tb_h1,
                               log2_tb_w1, log2_tb_h1,
                               log2_max_tb_s, cu_flags, &tu_info[3]);
        }

    } else {
        rcn_res_wrap(ctu_dec, x0, y0, log2_tb_w, log2_tb_h, cu_flags, tu_info);
    }
}



