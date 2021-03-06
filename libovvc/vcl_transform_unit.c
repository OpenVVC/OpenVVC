#include <string.h>

#include "ovutils.h"
#include "cabac_internal.h"
#include "vcl.h"
#include "dec_structures.h"
#include "ctudec.h"
#include "rcn.h"
#include "dbf_utils.h"
#include "drv_utils.h"
#include "drv.h"

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

static uint8_t
ovcabac_read_ae_root_cbf(OVCABACCtx *const cabac_ctx)
{
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    return ovcabac_ae_read(cabac_ctx, &cabac_state[QT_ROOT_CBF_CTX_OFFSET]);
}

uint8_t
ovcabac_read_ae_tu_cbf_luma_isp(OVCABACCtx *const cabac_ctx,
                                uint8_t prev_cbf)
{
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    int offset = 2 + prev_cbf;
    return ovcabac_ae_read(cabac_ctx, &cabac_state[QT_CBF_CTX_OFFSET + offset] );
}

uint8_t
ovcabac_read_ae_tu_cbf_luma(OVCABACCtx *const cabac_ctx)
{
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    return ovcabac_ae_read(cabac_ctx, &cabac_state[QT_CBF_CTX_OFFSET] );
}

uint8_t
ovcabac_read_ae_tu_cbf_cb(OVCABACCtx *const cabac_ctx)
{
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    return ovcabac_ae_read(cabac_ctx, &cabac_state[QT_CBF_CB_CTX_OFFSET] );
}

uint8_t
ovcabac_read_ae_tu_cbf_cr(OVCABACCtx *const cabac_ctx,
                          uint8_t tu_cbf_cb)
{
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    return ovcabac_ae_read(cabac_ctx, &cabac_state[QT_CBF_CR_CTX_OFFSET + tu_cbf_cb] );
}

uint8_t
ovcabac_read_ae_joint_cb_cr_flag(OVCABACCtx *const cabac_ctx,
                                 uint8_t cbf_mask_minus1)
{
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    return ovcabac_ae_read(cabac_ctx, &cabac_state[JOINT_CB_CR_FLAG_CTX_OFFSET + cbf_mask_minus1] );
}

static int
vvc_exp_golomb(OVCABACCtx *const cabac_ctx)
{
    int symbol = 0;
    int bit = 1;
    int add_val = 0;
    int count = 0;
    /*FIXME determine max_num_bits to be read in bypass */
     while (bit && count <= 32){
        bit = ovcabac_bypass_read(cabac_ctx);
        symbol += bit << count++;
     }
     if (--count){
         add_val = 0;
         while (count){
             add_val <<= 1;
             add_val |= ovcabac_bypass_read(cabac_ctx);
             count--;
         }
     }
     return symbol + add_val;
}

int
ovcabac_read_ae_cu_delta_qp(OVCABACCtx *const cabac_ctx)
{
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
   int delta_qp = ovcabac_ae_read(cabac_ctx, &cabac_state[DELTA_QP_CTX_OFFSET]);
   if (delta_qp)
       while (delta_qp < 5 && ovcabac_ae_read(cabac_ctx, &cabac_state[DELTA_QP_CTX_OFFSET + 1]))
       delta_qp ++;

   if (delta_qp >= 5){
       delta_qp += vvc_exp_golomb(cabac_ctx);
   }

   if (delta_qp){
       delta_qp = ovcabac_bypass_read(cabac_ctx) ? -delta_qp : delta_qp;
   }

   return delta_qp;
}

uint8_t
ovcabac_read_ae_cu_mts_flag(OVCABACCtx *const cabac_ctx)
{
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    return ovcabac_ae_read(cabac_ctx,&cabac_state[MTS_IDX_CTX_OFFSET]);
}

uint8_t
ovcabac_read_ae_cu_mts_idx(OVCABACCtx *const cabac_ctx)
{
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    uint8_t symbol = ovcabac_ae_read(cabac_ctx, &cabac_state[MTS_IDX_CTX_OFFSET + 1]);
    if(symbol && ovcabac_ae_read(cabac_ctx, &cabac_state[MTS_IDX_CTX_OFFSET + 2])){
        symbol++;
        if( ovcabac_ae_read(cabac_ctx, &cabac_state[MTS_IDX_CTX_OFFSET + 3])){
            symbol++;
        }
    }
    return symbol;
}

uint8_t
ovcabac_read_ae_transform_skip_luma_flag(OVCABACCtx *const cabac_ctx)
{
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    return ovcabac_ae_read(cabac_ctx, &cabac_state[TRANSFORM_SKIP_FLAG_CTX_OFFSET]);
}

uint8_t
ovcabac_read_ae_transform_skip_flag_c(OVCABACCtx *const cabac_ctx)
{
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    return ovcabac_ae_read(cabac_ctx, &cabac_state[TRANSFORM_SKIP_FLAG_CTX_OFFSET + 1]);
}

static inline int
decode_last_sig_prefix(OVCABACCtx *const cabac_ctx,
                       unsigned int log2_tb_d, unsigned offset_ctx)
{
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    /* FIXME this tab could be adapted by adding some unused ctx to last sig ctx */
    static const int prefix_ctx[8]  = { 0, 0, 0, 3, 6, 10, 15, 21 };
    int pos = 0;
    int ctx_offset, ctx_shift;
    int max_symbol = OVMIN(log2_tb_d, 5) << 1;

    ctx_offset = prefix_ctx[log2_tb_d];
    ctx_shift  = (log2_tb_d + 1) >> 2 ;

    while(--max_symbol > 0 && ovcabac_ae_read(cabac_ctx, &cabac_state[offset_ctx + ctx_offset + (pos >> ctx_shift)])){
        ++pos;
    }
    return pos;
}

static inline int
decode_last_sig_prefix_c(OVCABACCtx *const cabac_ctx,
                         unsigned int log2_tb_h, unsigned offset_ctx)
{
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    int pos = 0;
    int ctx_shift;
    int max_symbol = log2_tb_h << 1;

    ctx_shift = (1 << log2_tb_h) >> 3 ;
    ctx_shift = ov_clip(ctx_shift,0,2);

    while(--max_symbol > 0 && ovcabac_ae_read(cabac_ctx, &cabac_state[offset_ctx + (pos >> ctx_shift)])){
        ++pos;
    }

    return pos;
}

static inline int
decode_last_sig_suffix(OVCABACCtx *const cabac_ctx, int prefix)
{
    int num_bins = (prefix - 2) >> 1;
    int val = 0;
    while (num_bins > 0){
        val = val << 1 | ovcabac_bypass_read(cabac_ctx) ;
        --num_bins;
    }
    val = (1 << ((prefix >> 1) - 1)) * (2 + (prefix & 1)) + val;
    return val;
}

uint16_t
ovcabac_read_ae_last_sig_pos(OVCABACCtx *const cabac_ctx,
                             uint8_t log2_tb_w, uint8_t log2_tb_h)
{
    uint8_t last_x;
    uint8_t last_y;

    last_x = decode_last_sig_prefix(cabac_ctx, log2_tb_w, LAST_X_CTX_OFFSET);

    last_y = decode_last_sig_prefix(cabac_ctx, log2_tb_h, LAST_Y_CTX_OFFSET);

    if (last_x > 3) {
        last_x = decode_last_sig_suffix(cabac_ctx, last_x);
    }

    if (last_y > 3) {
        last_y = decode_last_sig_suffix(cabac_ctx, last_y);
    }

    return ((uint16_t) last_y << 8) | (last_x & 0xFF);
}

uint16_t
ovcabac_read_ae_last_sig_pos_c(OVCABACCtx *const cabac_ctx,
                               uint8_t log2_tb_w, uint8_t log2_tb_h)
{
    uint8_t last_x;
    uint8_t last_y;

    last_x = decode_last_sig_prefix_c(cabac_ctx, log2_tb_w, LAST_X_C_CTX_OFFSET);
    last_y = decode_last_sig_prefix_c(cabac_ctx, log2_tb_h, LAST_Y_C_CTX_OFFSET);

    if (last_x > 3){
        last_x = decode_last_sig_suffix(cabac_ctx, last_x);
    }

    if(last_y > 3){
        last_y = decode_last_sig_suffix(cabac_ctx, last_y);
    }

    return ((uint16_t) last_y << 8) | (last_x & 0xFF);
}

uint8_t
ovcabac_read_ae_lfnst_flag(OVCABACCtx *const cabac_ctx, uint8_t is_dual_tree)
{
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    uint8_t lfnst_flag = ovcabac_ae_read(cabac_ctx, &cabac_state[LFNST_IDX_CTX_OFFSET + is_dual_tree]);
    return lfnst_flag;
}

uint8_t
ovcabac_read_ae_lfnst_idx(OVCABACCtx *const cabac_ctx)
{
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    return ovcabac_ae_read(cabac_ctx, &cabac_state[LFNST_IDX_CTX_OFFSET + 2]);
}

static uint8_t
decode_cbf_st(const OVCTUDec *const ctu_dec, uint8_t rqt_root_cbf)
{
    OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;
    uint8_t tu_cbf_cb = ovcabac_read_ae_tu_cbf_cb(cabac_ctx);
    uint8_t tu_cbf_cr = ovcabac_read_ae_tu_cbf_cr(cabac_ctx, tu_cbf_cb);
    uint8_t cbf_mask = (tu_cbf_cb << 1) | tu_cbf_cr;
    uint8_t tu_cbf_luma = rqt_root_cbf;

    if (!rqt_root_cbf || (cbf_mask && rqt_root_cbf)){
        tu_cbf_luma = ovcabac_read_ae_tu_cbf_luma(cabac_ctx);
    }

    if (ctu_dec->jcbcr_enabled && cbf_mask) {
        uint8_t joint_cb_cr = ovcabac_read_ae_joint_cb_cr_flag(cabac_ctx, (cbf_mask & 0x3) - 1);
        cbf_mask |= joint_cb_cr << 3;
    }

    return cbf_mask | (tu_cbf_luma << 4);
}

static uint8_t
decode_cbf_c(const OVCTUDec *const ctu_dec)
{
    OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;
    uint8_t tu_cbf_cb = ovcabac_read_ae_tu_cbf_cb(cabac_ctx);
    uint8_t tu_cbf_cr = ovcabac_read_ae_tu_cbf_cr(cabac_ctx, tu_cbf_cb);
    uint8_t cbf_mask = (tu_cbf_cb << 1) | tu_cbf_cr;

    if (ctu_dec->jcbcr_enabled && cbf_mask) {
        uint8_t joint_cb_cr = ovcabac_read_ae_joint_cb_cr_flag(cabac_ctx,
                (cbf_mask & 0x3) - 1);
        cbf_mask |= joint_cb_cr << 3;
    }
    return cbf_mask;
}

static void
recon_isp_subtree_v(OVCTUDec *const ctudec,
                    unsigned int x0, unsigned int y0,
                    unsigned int log2_cb_w, unsigned int log2_cb_h,
                    uint8_t intra_mode, uint8_t cbf_flags,
                    int16_t lfnst_sb[4][16], uint8_t lfnst_flag, uint8_t lfnst_idx)
{
    #if 0
    struct OVDrvCtx *const pred_ctx = &ctudec->drv_ctx;
    #endif
    const struct TRFunctions *TRFunc = &ctudec->rcn_ctx.rcn_funcs.tr;

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
            fill_edge_map(&ctudec->dbf_info, x0 & ~0x3, y0, log2_pb_w >= 2 ? log2_pb_w : 2, log2_cb_h);
            fill_ctb_bound(&ctudec->dbf_info, x0 & ~0x3, y0, log2_pb_w >= 2 ? log2_pb_w : 2, log2_cb_h);
            fill_bs_map(&ctudec->dbf_info.bs2_map, x0 & ~0x3, y0, log2_pb_w >= 2 ? log2_pb_w : 2, log2_cb_h);
        }
        //fill_ctb_bound(&ctudec->dbf_info, x0 & ~0x3, y0, log2_pb_w, log2_cb_h);
        //fill_edge_map(&ctudec->dbf_info, x0 & ~0x3, y0, log2_pb_w, log2_cb_h);
        #endif

        if (cbf) {
            int16_t *coeffs_y = ctudec->residual_y + i * (1 << (log2_pb_w + log2_cb_h));

            //fill_bs_map(&ctudec->dbf_info.bs1_map, x0, y0, log2_pb_w, log2_cb_h);

            if (log2_pb_w) {
                int shift_v = 6 + 1;
                int shift_h = (6 + 15 - 1) - 10;
                DECLARE_ALIGNED(32, int16_t, tmp)[64*64];
                int16_t *src = coeffs_y;
                int16_t *dst = ctudec->transform_buff;
                int cb_h = 1 << log2_cb_h;

                memset(tmp, 0, sizeof(int16_t) << (log2_pb_w + log2_cb_h));
                if (lfnst_flag) {
                    process_lfnst_luma_isp(ctudec, coeffs_y, lfnst_sb[i], log2_pb_w, log2_cb_h,log2_cb_w, log2_cb_h, x0 -offset_x, y0,
                                       lfnst_idx);
                    /* lfnst forces IDCT II usage */
                    type_v = type_h = DCT_II;
                }

                TRFunc->func[type_v][log2_cb_h](src, tmp, pb_w, pb_w, cb_h, shift_v);
                TRFunc->func[type_h][log2_pb_w](tmp, dst, cb_h, cb_h, pb_w, shift_h);
            } else {
                int shift_h = (6 + 15 - 1) - 10;
                int cb_h = 1 << log2_cb_h;
                DECLARE_ALIGNED(32, int16_t, tmp)[64];

                memset(tmp, 0, sizeof(int16_t) << (log2_pb_w + log2_cb_h));

                TRFunc->func[type_v][log2_cb_h](coeffs_y, tmp, pb_w, pb_w, cb_h, shift_h + 1);

                memcpy(ctudec->transform_buff, tmp, sizeof(uint16_t) * (1 << log2_cb_h));
            }

            uint16_t *dst  = &ctudec->rcn_ctx.ctu_buff.y[x0 + y0 * RCN_CTB_STRIDE];
            int16_t *src  = ctudec->transform_buff;

            #if 1
            vvc_add_residual(src, dst, log2_pb_w, log2_cb_h, 0);
            #else
            rcn_func->ict[0](src, dst, log2_pb_w, log2_cb_h, 512);
            #endif
        }
        x0 += pb_w;
        offset_x += pb_w;
    }
}

static void
recon_isp_subtree_h(OVCTUDec *const ctudec,
                    unsigned int x0, unsigned int y0,
                    unsigned int log2_cb_w, unsigned int log2_cb_h,
                uint8_t intra_mode, uint8_t cbf_flags,
                int16_t lfnst_sb[4][16], uint8_t lfnst_flag, uint8_t lfnst_idx)
{
    const struct TRFunctions *TRFunc = &ctudec->rcn_ctx.rcn_funcs.tr;
    int log2_pb_h = log2_cb_h - 2;
    int nb_pb;
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
            fill_edge_map(&ctudec->dbf_info, x0, y0, log2_cb_w, log2_pb_h >=2 ? log2_pb_h : 2);
            fill_ctb_bound(&ctudec->dbf_info, x0, y0, log2_cb_w, log2_pb_h >=2 ? log2_pb_h : 2);
            fill_bs_map(&ctudec->dbf_info.bs2_map, x0, y0, log2_cb_w, log2_pb_h >=2 ? log2_pb_h : 2);
        }
        //fill_ctb_bound(&ctudec->dbf_info, x0, y0, log2_cb_w, log2_pb_h);
        //fill_edge_map(&ctudec->dbf_info, x0, y0, log2_cb_w, log2_pb_h);
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
                int shift_h = (6 + 15 - 1) - 10;
                int16_t *src = coeffs_y;
                int16_t *dst = ctudec->transform_buff;
                int cb_w = 1 << log2_cb_w;

                memset(tmp, 0, sizeof(int16_t) << (log2_cb_w + log2_pb_h));

                if (lfnst_flag) {
                    /*FIXME avoid distinguishing isp case in lfnst */
                    process_lfnst_luma_isp(ctudec, coeffs_y, lfnst_sb[i], log2_cb_w, log2_pb_h,log2_cb_w, log2_cb_h, x0, y0-offset_y,
                                       lfnst_idx);

                    /* lfnst forces IDCT II usage */
                    type_v = type_h = DCT_II;
                }

                TRFunc->func[type_v][log2_pb_h](src, tmp, cb_w, cb_w, pb_h, shift_v);
                TRFunc->func[type_h][log2_cb_w](tmp, dst, pb_h, pb_h, cb_w, shift_h);
            } else {
                int shift_h = (6 + 15 - 1) - 10;
                int cb_w = 1 << log2_cb_w;
                DECLARE_ALIGNED(32, int16_t, tmp)[64];

                memset(tmp, 0, sizeof(int16_t) << (log2_cb_w + log2_pb_h));

                TRFunc->func[type_h][log2_cb_w](coeffs_y, tmp, pb_h, pb_h, cb_w, shift_h + 1);

                memcpy(ctudec->transform_buff, tmp, sizeof(uint16_t) * (1 << log2_cb_w));
            }

            uint16_t *dst  = &ctudec->rcn_ctx.ctu_buff.y[x0 + y0 * RCN_CTB_STRIDE];
            int16_t *src  = ctudec->transform_buff;

        #if 1
        vvc_add_residual(src, dst, log2_cb_w, log2_pb_h, 0);
        #else
            rcn_func->ict[0](src, dst, log2_cb_w, log2_pb_h, 0);
        #endif
        }
        y0 += pb_h;
        offset_y += pb_h;
    }
}

static int
transform_tree(OVCTUDec *const ctu_dec,
               const OVPartInfo *const part_ctx,
               unsigned int x0, unsigned int y0,
               unsigned int log2_tb_w, unsigned int log2_tb_h,
               unsigned int log2_max_tb_s, uint8_t rqt_root_cbf,
               uint8_t cu_flags)
{
    uint8_t split_v = log2_tb_w > log2_max_tb_s;
    uint8_t split_h = log2_tb_h > log2_max_tb_s;

    if (split_v || split_h) {
        unsigned int tb_w1 = ((1 << log2_tb_w) >> split_v);
        unsigned int tb_h1 = ((1 << log2_tb_h) >> split_h);

        unsigned int log2_tb_w1 = log2_tb_w - split_v;
        unsigned int log2_tb_h1 = log2_tb_h - split_h;

        transform_tree(ctu_dec, part_ctx, x0, y0,
                       log2_tb_w1, log2_tb_h1,
                       log2_max_tb_s, rqt_root_cbf, cu_flags);
        if (split_v) {
            transform_tree(ctu_dec, part_ctx, x0 + tb_w1, y0,
                           log2_tb_w1, log2_tb_h1,
                           log2_max_tb_s, rqt_root_cbf, cu_flags);
        }

        if (split_h) {
            transform_tree(ctu_dec, part_ctx, x0, y0 + tb_h1,
                           log2_tb_w1, log2_tb_h1,
                           log2_max_tb_s, rqt_root_cbf, cu_flags);
        }

        if (split_h && split_v) {
            transform_tree(ctu_dec, part_ctx, x0 + tb_w1, y0 + tb_h1,
                           log2_tb_w1, log2_tb_h1, log2_max_tb_s, rqt_root_cbf, cu_flags);
        }

    } else {
        ctu_dec->transform_unit(ctu_dec, x0, y0, log2_tb_w, log2_tb_h, rqt_root_cbf, cu_flags);
    }

    return 0;
}

static uint8_t
isp_subtree_v(OVCTUDec *const ctu_dec,
              unsigned int x0, unsigned int y0,
              unsigned int log2_cb_w, unsigned int log2_cb_h,
              uint8_t intra_mode)
{
    OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;

    uint8_t cbf_flags = 0;
    int cbf = 0;
    int log2_pb_w = log2_cb_w - 2;
    int nb_pb;
    uint64_t sig_sb_map[4] = {0};
    uint8_t nb_coeffs[4] = {0};
    int16_t lfnst_sb[4][16];
    int i;

    /* height < 16 imposes restrictions on split dimension */
    if (log2_cb_h < 4 && (log2_pb_w <= (4 - log2_cb_h))) {
        log2_pb_w = 4 - log2_cb_h;
    }
    memset(ctu_dec->residual_y, 0, sizeof(int16_t) << (log2_cb_w + log2_cb_h));

    nb_pb = (1 << log2_cb_w) >> log2_pb_w;

    for (i = 0; i < nb_pb-1; ++i) {
        cbf = ovcabac_read_ae_tu_cbf_luma_isp(cabac_ctx, cbf);
        cbf_flags <<= 1;
        cbf_flags |= cbf;

        if (cbf) {
            uint16_t last_pos = ovcabac_read_ae_last_sig_pos(cabac_ctx, log2_pb_w, log2_cb_h);
            int16_t *coeffs_y = ctu_dec->residual_y + i * (1 << (log2_pb_w + log2_cb_h));
            uint64_t scan_map = 0xFDA6EB73C8419520;
            int last_y = last_pos >> 8;
            int last_x = last_pos & 0xFF;
            nb_coeffs[i] = (scan_map >> ((last_x + (last_y << 2)) << 2)) & 0xF;

            if (log2_pb_w < 2) {
                sig_sb_map[i] = 2;
                ctu_dec->residual_coding_isp_v(ctu_dec, coeffs_y, log2_pb_w, log2_cb_h, last_pos);
            }else {
                sig_sb_map[i] = ctu_dec->residual_coding(ctu_dec, coeffs_y, log2_pb_w, log2_cb_h, last_pos);
                memcpy(lfnst_sb[i], ctu_dec->lfnst_subblock, sizeof(uint16_t) * 16);
            }
        }
    }

    cbf = !cbf_flags ? 1 : ovcabac_read_ae_tu_cbf_luma_isp(cabac_ctx, cbf);
    cbf_flags <<= 1;
    cbf_flags |= cbf;

    if (cbf) {
        uint16_t last_pos = ovcabac_read_ae_last_sig_pos(cabac_ctx, log2_pb_w, log2_cb_h);
        int16_t *coeffs_y = ctu_dec->residual_y + i * (1 << (log2_pb_w + log2_cb_h));

        if (log2_pb_w <= 1) {
        sig_sb_map[i] = 2;
            ctu_dec->residual_coding_isp_v(ctu_dec, coeffs_y, log2_pb_w, log2_cb_h, last_pos);
        } else {
            uint64_t scan_map = 0xFDA6EB73C8419520;
            int last_y = last_pos >> 8;
            int last_x = last_pos & 0xFF;
            nb_coeffs[i] = (scan_map >> ((last_x + (last_y << 2)) << 2)) & 0xF;
            sig_sb_map[i] = ctu_dec->residual_coding(ctu_dec, coeffs_y, log2_pb_w, log2_cb_h, last_pos);
            memcpy(lfnst_sb[i], ctu_dec->lfnst_subblock, sizeof(uint16_t) * 16);
        }
    }
    uint8_t lfnst_flag = 0;
    uint8_t lfnst_idx = 0;

    if (ctu_dec->enable_lfnst) {
        int max_lfnst_pos = (log2_cb_h == log2_pb_w) && (log2_pb_w <= 3) ? 7 : 15;
        uint8_t can_lfnst = (sig_sb_map[0] | sig_sb_map[1] | sig_sb_map[2] | sig_sb_map[3]) == 1;

        can_lfnst &= nb_coeffs[0] <= max_lfnst_pos;
        can_lfnst &= nb_coeffs[1] <= max_lfnst_pos;
        can_lfnst &= nb_coeffs[2] <= max_lfnst_pos;
        can_lfnst &= nb_coeffs[3] <= max_lfnst_pos;

        if (can_lfnst) {
            uint8_t is_dual = ctu_dec->transform_unit != transform_unit_st;
            lfnst_flag = ovcabac_read_ae_lfnst_flag(cabac_ctx, is_dual);
            if (lfnst_flag) {
                lfnst_idx = ovcabac_read_ae_lfnst_idx(cabac_ctx);
            }
        }
    }

#if 1
    recon_isp_subtree_v(ctu_dec, x0, y0, log2_cb_w, log2_cb_h, intra_mode, cbf_flags,
                        lfnst_sb, lfnst_flag, lfnst_idx);
#endif

    return cbf_flags;
}

static uint8_t
isp_subtree_h(OVCTUDec *const ctu_dec,
              unsigned int x0, unsigned int y0,
              unsigned int log2_cb_w, unsigned int log2_cb_h,
              uint8_t intra_mode)
{
    OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;
    int i;

    uint8_t cbf = 0;
    uint8_t cbf_flags = 0;
    int log2_pb_h = log2_cb_h - 2;
    int nb_pb;
    int16_t *coeffs_y = ctu_dec->residual_y;
    uint64_t sig_sb_map[4] = {0};
    uint8_t nb_coeffs[4] = {0};
    int16_t lfnst_sb[4][16];

    /* width < 16 imposes restrictions on split numbers */
    if (log2_cb_w < 4 && (log2_pb_h <= (4 - log2_cb_w))) {
        log2_pb_h = 4 - log2_cb_w;
    }

    int tb_s = 1 << (log2_cb_w + log2_pb_h);
    nb_pb = (1 << log2_cb_h) >> log2_pb_h;

    for (i = 0; i < nb_pb - 1; ++i) {
        cbf = ovcabac_read_ae_tu_cbf_luma_isp(cabac_ctx, cbf);
        cbf_flags <<= 1;
        cbf_flags |= cbf;
        if (cbf) {
            uint16_t last_pos = ovcabac_read_ae_last_sig_pos(cabac_ctx, log2_cb_w, log2_pb_h);

            if (log2_pb_h <= 1) {
                ctu_dec->residual_coding_isp_h(ctu_dec, coeffs_y, log2_cb_w, log2_pb_h, last_pos);
        sig_sb_map[i] = 2;
            } else {
                uint64_t scan_map = 0xFDA6EB73C8419520;
                int last_y = last_pos >> 8;
                int last_x = last_pos & 0xFF;
                nb_coeffs[i] = (scan_map >> ((last_x + (last_y << 2)) << 2)) & 0xF;
                sig_sb_map [i] = ctu_dec->residual_coding(ctu_dec, coeffs_y, log2_cb_w, log2_pb_h, last_pos);
                memcpy(lfnst_sb[i], ctu_dec->lfnst_subblock, sizeof(uint16_t) * 16);
            }
        }
        coeffs_y += tb_s;
    }

    cbf = !cbf_flags ? 1 : ovcabac_read_ae_tu_cbf_luma_isp(cabac_ctx, cbf);
    cbf_flags <<= 1;
    cbf_flags |= cbf;

    if (cbf) {
        uint16_t last_pos = ovcabac_read_ae_last_sig_pos(cabac_ctx, log2_cb_w, log2_pb_h);

        if (log2_pb_h <= 1) {
            ctu_dec->residual_coding_isp_h(ctu_dec, coeffs_y, log2_cb_w, log2_pb_h, last_pos);
        sig_sb_map[i] = 2;
        } else {
            uint64_t scan_map = 0xFDA6EB73C8419520;
            int last_y = last_pos >> 8;
            int last_x = last_pos & 0xFF;
            nb_coeffs[i] = (scan_map >> ((last_x + (last_y << 2)) << 2)) & 0xF;
            sig_sb_map[i] = ctu_dec->residual_coding(ctu_dec, coeffs_y, log2_cb_w, log2_pb_h, last_pos);
            memcpy(lfnst_sb[i], ctu_dec->lfnst_subblock, sizeof(uint16_t) * 16);
        }
    }

    uint8_t lfnst_flag = 0;
    uint8_t lfnst_idx = 0;
    if (ctu_dec->enable_lfnst) {
        uint8_t can_lfnst = (sig_sb_map[0] | sig_sb_map[1] | sig_sb_map[2] | sig_sb_map[3]) == 1;
        int max_lfnst_pos = (log2_pb_h == log2_cb_w) && (log2_cb_w <= 3) ? 7 : 15;

        can_lfnst &= nb_coeffs[0] <= max_lfnst_pos;
        can_lfnst &= nb_coeffs[1] <= max_lfnst_pos;
        can_lfnst &= nb_coeffs[2] <= max_lfnst_pos;
        can_lfnst &= nb_coeffs[3] <= max_lfnst_pos;

        if (can_lfnst) {
            uint8_t is_dual = ctu_dec->transform_unit != transform_unit_st;
            lfnst_flag = ovcabac_read_ae_lfnst_flag(cabac_ctx, is_dual);
            if (lfnst_flag) {
                lfnst_idx = ovcabac_read_ae_lfnst_idx(cabac_ctx);
            }
        }
    }

#if 1
    recon_isp_subtree_h(ctu_dec, x0, y0, log2_cb_w, log2_cb_h, intra_mode, cbf_flags,
                        lfnst_sb, lfnst_flag, lfnst_idx);
#endif

    return cbf_flags;
}

static int
transform_unit(OVCTUDec *const ctu_dec,
               unsigned int x0, unsigned int y0,
               unsigned int log2_tb_w, unsigned int log2_tb_h,
               uint8_t tu_cbf_luma, uint8_t cu_flags)
{

    if (tu_cbf_luma) {
        OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;
        uint8_t cu_mts_flag = 0;
        uint8_t cu_mts_idx = 0;
        uint8_t transform_skip_flag = 0;
        int16_t *const coeffs_y = ctu_dec->residual_y;
        #if 0
        int cu_qp_delta = 0;
        #endif


#if 1
        /*FIXME move bs map filling to to cbf_flag reading */
        fill_bs_map(&ctu_dec->dbf_info.bs1_map, x0, y0, log2_tb_w, log2_tb_h);
#endif

        if (ctu_dec->transform_skip_enabled && log2_tb_w <= ctu_dec->max_log2_transform_skip_size
                && log2_tb_h <= ctu_dec->max_log2_transform_skip_size && !(cu_flags &flg_isp_flag  )) {
            transform_skip_flag = ovcabac_read_ae_transform_skip_luma_flag(cabac_ctx);
        }


        if (!transform_skip_flag) {
            int lfnst_flag = 0;
            int lfnst_idx = 0;
            uint16_t last_pos = ovcabac_read_ae_last_sig_pos(cabac_ctx, log2_tb_w, log2_tb_h);
            int lim_cg_w = ((((last_pos >> 8)) >> 2) + (((last_pos & 0xFF))>> 2) + 1) << 2;
            uint8_t is_dc_c = !last_pos;
            uint64_t sig_sb_map;
            sig_sb_map = ctu_dec->residual_coding(ctu_dec, coeffs_y, log2_tb_w, log2_tb_h,
                                                 last_pos);

            if (ctu_dec->enable_lfnst && sig_sb_map == 0x1) {
                int max_lfnst_pos = (log2_tb_h == log2_tb_w) && (log2_tb_w <= 3) ? 7 : 15;
                int last_y = last_pos >> 8;
                int last_x = last_pos & 0xFF;
                uint8_t is_mip = !!(cu_flags & flg_mip_flag);
                uint8_t allow_mip_lfnst = !is_mip || (log2_tb_h >= 4 && log2_tb_w >= 4);
                uint64_t scan_map = 0xFDA6EB73C8419520;
                int nb_coeffs = (scan_map >> ((last_x + (last_y << 2)) << 2)) & 0xF;

                if (allow_mip_lfnst && nb_coeffs <= max_lfnst_pos && !is_dc_c) {
                    uint8_t is_dual = ctu_dec->transform_unit != transform_unit_st;
                    lfnst_flag = ovcabac_read_ae_lfnst_flag(cabac_ctx, is_dual);
                    if (lfnst_flag) {
                        lfnst_idx = ovcabac_read_ae_lfnst_idx(cabac_ctx);
                    }
                }
            }

            if (!lfnst_flag && !is_dc_c && ctu_dec->mts_enabled && (log2_tb_w < 6) && (log2_tb_h < 6)
                                    && !(sig_sb_map & (~0x000000000F0F0F0F))) {
                cu_mts_flag = ovcabac_read_ae_cu_mts_flag(cabac_ctx);
                if (cu_mts_flag) {
                    cu_mts_idx = ovcabac_read_ae_cu_mts_idx(cabac_ctx);
                }
            }

#if 1
            uint8_t is_mip = !!(cu_flags & flg_mip_flag);

            rcn_residual(ctu_dec, ctu_dec->transform_buff, coeffs_y, x0, y0, log2_tb_w, log2_tb_h,
                           lim_cg_w, cu_mts_flag, cu_mts_idx, is_dc_c, lfnst_flag, is_mip, lfnst_idx);
#endif

        } else { //Renorm transform skipped residuals
            #if 0
            int x_pu = x0 >> 2;
            int y_pu = y0 >> 2;
            #endif
            ctu_dec->dequant_skip = &ctu_dec->dequant_luma_skip;
            residual_coding_ts(ctu_dec, log2_tb_w, log2_tb_h);
            //FIXME transform residual is currently performed in the dequant function
        }

        #if 1
        vvc_add_residual( ctu_dec->transform_buff, &ctu_dec->rcn_ctx.ctu_buff.y[x0 + y0 * RCN_CTB_STRIDE], log2_tb_w, log2_tb_h, 0);
        #else
        rcn_func->ict[0](ctu_dec->transform_buff, &ctu_dec->rcn_ctx.ctu_buff.y[x0 + y0 * RCN_CTB_STRIDE],
                             log2_tb_w, log2_tb_h, 0);
        #endif
    }
    return 0;
}

static int
transform_unit_chroma(OVCTUDec *const ctu_dec,
                      unsigned int x0, unsigned int y0,
                      unsigned int log2_tb_w, unsigned int log2_tb_h,
                      uint8_t cbf_mask, uint8_t cu_flags)
{
    uint8_t joint_cb_cr = cbf_mask & (1 << 3);
    const struct RCNFunctions *const rcn_func = &ctu_dec->rcn_ctx.rcn_funcs;

    cbf_mask &= 0x3;

    if (!joint_cb_cr) {
        int16_t *const coeffs_cb = ctu_dec->residual_cb;
        int16_t *const coeffs_cr = ctu_dec->residual_cr;
        uint16_t last_pos_cb = (1 << (log2_tb_h + log2_tb_w)) - 1;
        uint16_t last_pos_cr = (1 << (log2_tb_h + log2_tb_w)) - 1;
        int lim_cg_w_cb = OVMAX(1 << log2_tb_w, 1 << log2_tb_h);
        int lim_cg_w_cr = OVMAX(1 << log2_tb_w, 1 << log2_tb_h);
        int16_t tmp_lfnst_cb[16];
        int16_t tmp_lfnst_cr[16];
        uint8_t lfnst_flag = 0;
        uint8_t lfnst_idx = 0;
        uint8_t transform_skip_flag = 0;

        /* FIXME move dequant to reconstruction this require modification in coeff
           reading */

        if (cbf_mask & 0x2) {
            OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;

            /* TODO use max_tr_skip_s = 0 to avoid using enabled test */
            if (ctu_dec->transform_skip_enabled && log2_tb_w <= ctu_dec->max_log2_transform_skip_size
                && log2_tb_h <= ctu_dec->max_log2_transform_skip_size) {
                OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;
                transform_skip_flag |= ovcabac_read_ae_transform_skip_flag_c(cabac_ctx) << 1;
            }

            if (!(transform_skip_flag & 0x2)) {
#if 1
                ctu_dec->dequant_chroma = &ctu_dec->dequant_cb;
#endif

                last_pos_cb = ovcabac_read_ae_last_sig_pos_c(cabac_ctx, log2_tb_w, log2_tb_h);

                lim_cg_w_cb = ((((last_pos_cb >> 8)) >> 2) + (((last_pos_cb & 0xFF))>> 2) + 1) << 2;

                last_pos_cb = ctu_dec->residual_coding_chroma(ctu_dec, coeffs_cb, log2_tb_w, log2_tb_h, last_pos_cb);

                /*FIXME avoid copy of lfnst_sb */
                memcpy(tmp_lfnst_cb, ctu_dec->lfnst_subblock, sizeof(int16_t) * 16);
            }  else {
                /* FIXME Chroma TS */
                ctu_dec->dequant_skip = &ctu_dec->dequant_cb_skip;
                residual_coding_ts(ctu_dec, log2_tb_w, log2_tb_h);
                memcpy(coeffs_cb ,ctu_dec->transform_buff, sizeof(uint16_t) << (log2_tb_h + log2_tb_w));
            }
        }

        if (cbf_mask & 0x1) {
            OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;

            /* TODO use max_tr_skip_s = 0 to avoid using enabled test */
            if (ctu_dec->transform_skip_enabled && log2_tb_w <= ctu_dec->max_log2_transform_skip_size
                                               && log2_tb_h <= ctu_dec->max_log2_transform_skip_size) {
                OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;
                transform_skip_flag |= ovcabac_read_ae_transform_skip_flag_c(cabac_ctx);
            }

            if (!(transform_skip_flag & 0x1)) {
#if 1
                ctu_dec->dequant_chroma = &ctu_dec->dequant_cr;
#endif

                last_pos_cr = ovcabac_read_ae_last_sig_pos_c(cabac_ctx, log2_tb_w, log2_tb_h);

                lim_cg_w_cr = ((((last_pos_cr >> 8)) >> 2) + (((last_pos_cr & 0xFF))>> 2) + 1) << 2;

                last_pos_cr = ctu_dec->residual_coding_chroma(ctu_dec, coeffs_cr, log2_tb_w, log2_tb_h, last_pos_cr);
            } else {
                ctu_dec->dequant_skip = &ctu_dec->dequant_cr_skip;
                residual_coding_ts(ctu_dec, log2_tb_w, log2_tb_h);
                memcpy(coeffs_cr, ctu_dec->transform_buff, sizeof(uint16_t) << (log2_tb_h + log2_tb_w));
            }

            /*FIXME avoid copy of lfnst_sb */
            memcpy(tmp_lfnst_cr, ctu_dec->lfnst_subblock, sizeof(int16_t) * 16);
        }

        if (ctu_dec->enable_lfnst && !transform_skip_flag) {
            uint8_t need_cb_chk = cbf_mask & 0x2;
            uint8_t need_cr_chk = cbf_mask & 0x1;

            int max_lfnst_pos = (log2_tb_h == log2_tb_w) && (log2_tb_w <= 3) ? 7 : 15;

            uint8_t can_lfnst = !!cbf_mask;

            /*FIXME better can_lfnst derivation */
            if (need_cb_chk && need_cr_chk) {
                can_lfnst &= last_pos_cr <= max_lfnst_pos && last_pos_cb <= max_lfnst_pos;
                can_lfnst &= !!(last_pos_cr | last_pos_cb);
            } else if (need_cb_chk) {
                can_lfnst &= last_pos_cb <= max_lfnst_pos;
                can_lfnst &= !!last_pos_cb;
            } else {
                can_lfnst &= last_pos_cr <= max_lfnst_pos;
                can_lfnst &= !!last_pos_cr;
            }

            if (can_lfnst) {
                OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;
                uint8_t is_dual = ctu_dec->transform_unit != transform_unit_st;
                lfnst_flag = ovcabac_read_ae_lfnst_flag(cabac_ctx, is_dual);
                if (lfnst_flag) {
                    lfnst_idx = ovcabac_read_ae_lfnst_idx(cabac_ctx);
                }
            }

        }


        #if 1
        if (cbf_mask & 0x2) {
            uint16_t *const dst_cb = &ctu_dec->rcn_ctx.ctu_buff.cb[(x0) + (y0 * RCN_CTB_STRIDE)];
            int16_t scale  = ctu_dec->lmcs_info.lmcs_chroma_scale;
            int16_t *const coeffs_cb = ctu_dec->residual_cb;

            if (!(transform_skip_flag & 0x2)) {
                rcn_residual_c(ctu_dec, ctu_dec->transform_buff, coeffs_cb, tmp_lfnst_cb, x0, y0, log2_tb_w, log2_tb_h,
                               lim_cg_w_cb, 0, 0, !last_pos_cb, lfnst_flag, 1, lfnst_idx);
            } else {
                memcpy(ctu_dec->transform_buff, coeffs_cb, sizeof(uint16_t) << (log2_tb_h + log2_tb_w));
            }

            rcn_func->ict[0](ctu_dec->transform_buff, dst_cb, log2_tb_w, log2_tb_h, scale);
#if 0
            (*ctu_dec->scale_addsub_residuals)[0](ctu_dec->transform_buff, dst_cb, log2_tb_w, log2_tb_h, scale);

#endif
            fill_bs_map(&ctu_dec->dbf_info.bs1_map_cb, x0 << 1, y0 << 1, log2_tb_w + 1, log2_tb_h + 1);
        }

        if (cbf_mask & 0x1) {
            uint16_t *const dst_cr = &ctu_dec->rcn_ctx.ctu_buff.cr[(x0) + (y0 * RCN_CTB_STRIDE)];
            int16_t scale  = ctu_dec->lmcs_info.lmcs_chroma_scale;

            if (!(transform_skip_flag & 0x1)) {
                rcn_residual_c(ctu_dec, ctu_dec->transform_buff, coeffs_cr, tmp_lfnst_cr, x0, y0, log2_tb_w, log2_tb_h,
                               lim_cg_w_cr, 0, 0, !last_pos_cr, lfnst_flag, 1, lfnst_idx);
            } else {
                memcpy(ctu_dec->transform_buff, coeffs_cr, sizeof(uint16_t) << (log2_tb_h + log2_tb_w));
            }

            rcn_func->ict[0](ctu_dec->transform_buff, dst_cr, log2_tb_w, log2_tb_h, scale);
#if 0
            (*ctu_dec->scale_addsub_residuals)[0](ctu_dec->transform_buff, dst_cr, log2_tb_w, log2_tb_h, scale);

#endif
            fill_bs_map(&ctu_dec->dbf_info.bs1_map_cr, x0 << 1, y0 << 1, log2_tb_w + 1, log2_tb_h + 1);
        }
        #endif

    } else { //Joint cb cr (ICT)
        OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;
        uint16_t last_pos_cbcr;
        uint16_t *const dst_cb = &ctu_dec->rcn_ctx.ctu_buff.cb[(x0) + (y0 * RCN_CTB_STRIDE)];
        uint16_t *const dst_cr = &ctu_dec->rcn_ctx.ctu_buff.cr[(x0) + (y0 * RCN_CTB_STRIDE)];
        int16_t scale = ctu_dec->lmcs_info.lmcs_chroma_scale;
        int16_t *const coeffs_jcbcr = ctu_dec->residual_cb;
        uint8_t lfnst_flag = 0;
        uint8_t lfnst_idx = 0;
        uint8_t transform_skip_flag = 0;
        int lim_cg_w_cbcr = 0;

        /* FIXME this is hackish joint cb cr involves a different delta qp from
          previous ones */

        #if 1
        if (cbf_mask == 3) {
            ctu_dec->dequant_chroma = &ctu_dec->dequant_joint_cb_cr;
            ctu_dec->dequant_skip = &ctu_dec->dequant_jcbcr_skip;
        } else if (cbf_mask == 1) {
            ctu_dec->dequant_chroma = &ctu_dec->dequant_cr;
            ctu_dec->dequant_skip = &ctu_dec->dequant_cr_skip;
        } else {
            ctu_dec->dequant_chroma = &ctu_dec->dequant_cb;
            ctu_dec->dequant_skip = &ctu_dec->dequant_cb_skip;
        }
        #endif

        if (ctu_dec->transform_skip_enabled && log2_tb_w <= ctu_dec->max_log2_transform_skip_size
            && log2_tb_h <= ctu_dec->max_log2_transform_skip_size) {
            OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;
            transform_skip_flag = ovcabac_read_ae_transform_skip_flag_c(cabac_ctx);
        }

        if (!transform_skip_flag) {
            last_pos_cbcr = ovcabac_read_ae_last_sig_pos_c(cabac_ctx, log2_tb_w, log2_tb_h);

            lim_cg_w_cbcr = ((((last_pos_cbcr >> 8)) >> 2) + (((last_pos_cbcr & 0xFF))>> 2) + 1) << 2;

            last_pos_cbcr = ctu_dec->residual_coding_chroma(ctu_dec, coeffs_jcbcr, log2_tb_w, log2_tb_h,
                                                            last_pos_cbcr);
        } else {
            residual_coding_ts(ctu_dec, log2_tb_w, log2_tb_h);
        }

        if (ctu_dec->enable_lfnst && !transform_skip_flag) {
            int max_lfnst_pos = (log2_tb_h == log2_tb_w) && (log2_tb_w <= 3) ? 7 : 15;

            uint8_t can_lfnst =  last_pos_cbcr <= max_lfnst_pos;

            can_lfnst &= !!last_pos_cbcr;

            if (can_lfnst) {
                uint8_t is_dual = ctu_dec->transform_unit != transform_unit_st;
                lfnst_flag = ovcabac_read_ae_lfnst_flag(cabac_ctx, is_dual);
                if (lfnst_flag) {
                    lfnst_idx = ovcabac_read_ae_lfnst_idx(cabac_ctx);
                }
            }
        }

        if (!transform_skip_flag) {
            rcn_residual_c(ctu_dec, ctu_dec->transform_buff, coeffs_jcbcr,
                           ctu_dec->lfnst_subblock, x0, y0, log2_tb_w, log2_tb_h,
                           lim_cg_w_cbcr, 0, 0, !last_pos_cbcr, lfnst_flag, 1, lfnst_idx);
        }

#if 0
        if (cbf_mask == 3) {
            (*ctu_dec->scale_addsub_residuals)[0](src, dst_cb, log2_tb_w, log2_tb_h, scale);
            (*ctu_dec->scale_addsub_residuals)[1](src, dst_cr, log2_tb_w, log2_tb_h, scale);
        } else if (cbf_mask == 2) {
            (*ctu_dec->scale_addsub_residuals)[0](src, dst_cb, log2_tb_w, log2_tb_h, scale);
            (*ctu_dec->scale_addsub_residuals)[2](src, dst_cr, log2_tb_w, log2_tb_h, scale);
        } else {
            (*ctu_dec->scale_addsub_residuals)[0](src, dst_cr, log2_tb_w, log2_tb_h, scale);
            (*ctu_dec->scale_addsub_residuals)[2](src, dst_cb, log2_tb_w, log2_tb_h, scale);
        }

#else
        fill_bs_map(&ctu_dec->dbf_info.bs1_map_cb, x0 << 1, y0 << 1, log2_tb_w + 1, log2_tb_h + 1);
        fill_bs_map(&ctu_dec->dbf_info.bs1_map_cr, x0 << 1, y0 << 1, log2_tb_w + 1, log2_tb_h + 1);
        /* FIXME better organisation based on cbf_mask */
        if (cbf_mask == 3) {
            rcn_func->ict[0](ctu_dec->transform_buff, dst_cb, log2_tb_w, log2_tb_h, scale);
            rcn_func->ict[1](ctu_dec->transform_buff, dst_cr, log2_tb_w, log2_tb_h, scale);
        } else if (cbf_mask == 2) {
            rcn_func->ict[0](ctu_dec->transform_buff, dst_cb, log2_tb_w, log2_tb_h, scale);
            rcn_func->ict[2](ctu_dec->transform_buff, dst_cr, log2_tb_w, log2_tb_h, scale);
        } else {
            rcn_func->ict[0](ctu_dec->transform_buff, dst_cr, log2_tb_w, log2_tb_h, scale);
            rcn_func->ict[2](ctu_dec->transform_buff, dst_cb, log2_tb_w, log2_tb_h, scale);
        }
#endif
    }
    return 0;
}

int
transform_unit_st(OVCTUDec *const ctu_dec,
                  unsigned int x0, unsigned int y0,
                  unsigned int log2_tb_w, unsigned int log2_tb_h,
                  uint8_t rqt_root_cbf, uint8_t cu_flags)
{
    uint8_t cbf_mask = decode_cbf_st(ctu_dec, rqt_root_cbf);

    /* FIXME check if delta_qp is read per cu or per tu */
    if (cbf_mask) {
        if (ctu_dec->delta_qp_enabled && cbf_mask) {
            OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;
            int cu_qp_delta = ovcabac_read_ae_cu_delta_qp(cabac_ctx);
            #if 1
            derive_dequant_ctx(ctu_dec, &ctu_dec->qp_ctx, cu_qp_delta);
            #endif
        }

        if (cbf_mask & (1 << 4)) {
            transform_unit(ctu_dec, x0, y0, log2_tb_w, log2_tb_h, 1, cu_flags);
        }

        if (cbf_mask & 0xF) {
            transform_unit_chroma(ctu_dec, x0 >> 1, y0 >> 1, log2_tb_w - 1,
                                  log2_tb_h - 1, cbf_mask, cu_flags);
        }
    }


    return 0;
}

int
transform_unit_l(OVCTUDec *const ctu_dec,
                  unsigned int x0, unsigned int y0,
                  unsigned int log2_tb_w, unsigned int log2_tb_h,
                  uint8_t rqt_root_cbf, uint8_t cu_flags)
{
    OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;
    uint8_t cbf_mask = ovcabac_read_ae_tu_cbf_luma(cabac_ctx);

    if (cbf_mask) {
        if (ctu_dec->delta_qp_enabled && cbf_mask) {
            int cu_qp_delta = ovcabac_read_ae_cu_delta_qp(cabac_ctx);
            #if 1
            derive_dequant_ctx(ctu_dec, &ctu_dec->qp_ctx, cu_qp_delta);
            #endif
        }

        transform_unit(ctu_dec, x0, y0, log2_tb_w, log2_tb_h, cbf_mask, cu_flags);
    }

    return 0;
}

int
transform_unit_c(OVCTUDec *const ctu_dec,
                  unsigned int x0, unsigned int y0,
                  unsigned int log2_tb_w, unsigned int log2_tb_h,
                  uint8_t rqt_root_cbf, uint8_t cu_flags)
{
    OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;
    uint8_t cbf_mask = decode_cbf_c(ctu_dec);

    if (cbf_mask) {
        if (ctu_dec->delta_qp_enabled && cbf_mask) {
            int cu_qp_delta = ovcabac_read_ae_cu_delta_qp(cabac_ctx);
            #if 1
            derive_dequant_ctx(ctu_dec, &ctu_dec->qp_ctx, cu_qp_delta);
            #endif
        }

        transform_unit_chroma(ctu_dec, x0, y0, log2_tb_w, log2_tb_h, cbf_mask, cu_flags);
    }

    return 0;
}

int
transform_unit_wrap(OVCTUDec *const ctu_dec,
                    const OVPartInfo *const part_ctx,
                    uint8_t x0, uint8_t y0,
                    uint8_t log2_cb_w, uint8_t log2_cb_h,
                    VVCCU cu)
{
    if (cu.cu_flags & flg_pred_mode_flag) {
        /* INTRA */
        if (!(cu.cu_flags & flg_isp_flag)) {
            /*FIXME check if part_ctx mandatory for transform_tree */
            transform_tree(ctu_dec, part_ctx, x0, y0, log2_cb_w, log2_cb_h,
                           part_ctx->log2_max_tb_s, 0, cu.cu_flags);
        } else {
            uint8_t isp_mode = cu.cu_opaque;
            uint8_t intra_mode = cu.cu_mode_idx;
            if (isp_mode) {
                /* Disable CCLM in 64x64 ISP CU*/
                ctu_dec->tmp_disable_cclm = log2_cb_h == log2_cb_w && log2_cb_h == 6;

                if (isp_mode == 2) {
                    isp_subtree_v(ctu_dec, x0, y0, log2_cb_w, log2_cb_h, intra_mode);

                    //return cu;

                } else if (isp_mode == 1) {
                    isp_subtree_h(ctu_dec, x0, y0, log2_cb_w, log2_cb_h, intra_mode);

                    ////return cu;
                }
            }
        }
    } else {
        /* INTER (should probably be default in a swicth*/
        /*FIXME move root_cbf_read into transform_tree */
        OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;
        uint8_t merge_flag   = !!(cu.cu_flags & flg_merge_flag);
        uint8_t rqt_root_cbf = merge_flag || ovcabac_read_ae_root_cbf(cabac_ctx);

        if (rqt_root_cbf) {
            transform_tree(ctu_dec, part_ctx, x0, y0, log2_cb_w, log2_cb_h,
                           part_ctx->log2_max_tb_s, 1, cu.cu_flags);
        }
    }
    return 0;
}
