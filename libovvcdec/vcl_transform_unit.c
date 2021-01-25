
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
    int max_symbol = FFMIN(log2_tb_d, 5) << 1;

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
    ctx_shift = av_clip(ctx_shift,0,2);

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
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    uint8_t last_x;
    uint8_t last_y;

    last_x = decode_last_sig_prefix(cabac_ctx, log2_tb_w, LAST_X_CTX_OFFSET);

    last_y = decode_last_sig_prefix(cabac_ctx, log2_tb_h, LAST_Y_CTX_OFFSET);

    if (last_x > 3){
        last_x = decode_last_sig_suffix(cabac_ctx, last_x);
    }

    if(last_y > 3){
        last_y = decode_last_sig_suffix(cabac_ctx, last_y);
    }

    return ((uint16_t) last_y << 8) | last_x & 0xFF;
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

    return ((uint16_t) last_y << 8) | last_x & 0xFF;
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

uint8_t
ovcabac_read_ae_significant_cg_flag(OVCABACCtx *const cabac_ctx,
                                    uint8_t got_significant_neighbour)
{
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    return ovcabac_ae_read(cabac_ctx, &ctx_table[SIG_COEFF_GROUP_CTX_OFFSET + got_significant_neighbour]);
}


uint8_t
ovcabac_read_ae_significant_cg_flag_chroma(OVCABACCtx *const cabac_ctx,
                                           uint8_t got_significant_neighbour)
{
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    return ovcabac_ae_read(cabac_ctx, &ctx_table[SIG_COEFF_GROUP_C_CTX_OFFSET + got_significant_neighbour]);
}

uint8_t
ovcabac_read_ae_significant_ts_cg_flag(OVCABACCtx *const cabac_ctx,
                                       uint8_t got_significant_neighbour)
{
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    return ovcabac_ae_read(cabac_ctx, &ctx_table[TS_SIG_COEFF_GROUP_CTX_OFFSET + got_significant_neighbour]);
}

static uint8_t
decode_cbf_st(const OVCTUDec *const ctu_dec, uint8_t rqt_root_cbf)
{
    VVCCABACContext *const cabac_ctx = ctu_dec->cabac_ctx;
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
    VVCCABACContext *const cabac_ctx = ctu_dec->cabac_ctx;
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

static int
transform_tree(OVCTUDec *const ctu_dec,
               const VVCPartSize *const part_ctx,
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
    VVCCTUPredContext *const pred_ctx = &ctu_dec->pred_ctx;
    VVCCABACContext *const cabac_ctx = ctu_dec->cabac_ctx;

    uint8_t cbf_flags = 0;
    int cbf = 0;
    int log2_pb_w = log2_cb_w - 2;
    int nb_pb;
    uint64_t sig_sb_map[4] = {0};
    uint8_t nb_coeffs[4] = {0};
    uint16_t lfnst_sb[4][16];
    int i;

    /* height < 16 imposes restrictions on split dimension */
    if (log2_cb_h < 4 && (log2_pb_w <= (4 - log2_cb_h))) {
        log2_pb_w = 4 - log2_cb_h;
    }

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
    uint8_t lfnst_idx;

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

#if 0
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
    VVCCABACContext *const cabac_ctx = ctu_dec->cabac_ctx;
    int i;

    uint8_t cbf = 0;
    uint8_t cbf_flags = 0;
    int log2_pb_h = log2_cb_h - 2;
    int nb_pb;
    int16_t *coeffs_y = ctu_dec->residual_y;
    uint64_t sig_sb_map[4] = {0};
    uint8_t nb_coeffs[4] = {0};
    uint16_t lfnst_sb[4][16];

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
    uint8_t lfnst_idx;
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

#if 0
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
        VVCCABACContext *const cabac_ctx = ctu_dec->cabac_ctx;
        uint8_t cu_mts_flag = 0;
        uint8_t cu_mts_idx = 0;
        uint8_t transform_skip_flag = 0;
        int16_t *const coeffs_y = ctu_dec->residual_y;
        int is_lfnst;
        int cu_qp_delta = 0;


#if 0
        /*FIXME move bs map filling to to cbf_flag reading */
        fill_bs_map(&ctu_dec->dbf_info.bs1_map, x0, y0, log2_tb_w, log2_tb_h);
#endif

        if (ctu_dec->transform_skip_enabled && log2_tb_w <= ctu_dec->max_log2_transform_skip_size
                && log2_tb_h <= ctu_dec->max_log2_transform_skip_size && !(cu_flags &flg_isp_flag  )) {
            transform_skip_flag = ovcabac_read_ae_transform_skip_luma_flag(cabac_ctx);
        }


        if (!transform_skip_flag) {
            int lfnst_flag = 0;
            int lfnst_idx;
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

#if 0
            uint8_t is_mip = !!(cu_flags & flg_mip_flag);

            recon_residual(ctu_dec, ctu_dec->transform_buff, coeffs_y, x0, y0, log2_tb_w, log2_tb_h,
                           lim_cg_w, cu_mts_flag, cu_mts_idx, is_dc_c, lfnst_flag, is_mip, lfnst_idx);
#endif

        } else { //Renorm transform skipped residuals
            int x_pu = x0 >> 2;
            int y_pu = y0 >> 2;
            is_lfnst = residual_coding_ts(ctu_dec, log2_tb_w, log2_tb_h);
            //FIXME transform residual is currently performed in the dequant function
        }
#if 0
        int16_t *dst  = &ctu_dec->ctu_data_y[VVC_CTB_OFFSET + x0 + VVC_CTB_STRIDE * y0];
        int16_t *src  = ctu_dec->transform_buff;

        vvc_dsp_context.vvc_transform_add(src, dst, log2_tb_w, log2_tb_h, 0);
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
    uint8_t force_lfnst = 0;

    cbf_mask &= 0x3;

    if (!joint_cb_cr) {
        int16_t *const coeffs_cb = ctu_dec->residual_cb;
        int16_t *const coeffs_cr = ctu_dec->residual_cr;
        uint16_t last_pos_cb; /* Max == 64 -> uint8_t */
        uint16_t last_pos_cr; /* Max == 64 -> uint8_t */
        int lim_cg_w_cb; /* Max == 64 -> uint8_t */
        int lim_cg_w_cr; /* Max == 64 -> uint8_t */
        int16_t tmp_lfnst_cb[16];
        int16_t tmp_lfnst_cr[16];
        uint8_t lfnst_flag = 0;
        uint8_t lfnst_idx; 

        /* FIXME move dequant to reconstruction this require modification in coeff
           reading */

        if (cbf_mask & 0x2) {
            VVCCABACContext *const cabac_ctx = ctu_dec->cabac_ctx;

            /* TODO use max_tr_skip_s = 0 to avoid using enabled test */
            if (ctu_dec->transform_skip_enabled && log2_tb_w <= ctu_dec->max_log2_transform_skip_size
                                               && log2_tb_h <= ctu_dec->max_log2_transform_skip_size) {
                VVCCABACContext *const cabac_ctx = ctu_dec->cabac_ctx;
                uint8_t transform_skip_flag = ovcabac_read_ae_transform_skip_flag_c(cabac_ctx);
            }

            #if 0
            ctu_dec->dequant_chroma = &ctu_dec->dequant_cb;
            #endif

            last_pos_cb = ovcabac_read_ae_last_sig_pos_c(cabac_ctx, log2_tb_w, log2_tb_h);

            lim_cg_w_cb = ((((last_pos_cb >> 8)) >> 2) + (((last_pos_cb & 0xFF))>> 2) + 1) << 2;

            last_pos_cb = ctu_dec->residual_coding_chroma(ctu_dec, coeffs_cb, log2_tb_w, log2_tb_h, last_pos_cb);

            /*FIXME avoid copy of lfnst_sb */
            memcpy(tmp_lfnst_cb, ctu_dec->lfnst_subblock, sizeof(int16_t) * 16);
        }

        if (cbf_mask & 0x1) {
            VVCCABACContext *const cabac_ctx = ctu_dec->cabac_ctx;

            /* TODO use max_tr_skip_s = 0 to avoid using enabled test */
            if (ctu_dec->transform_skip_enabled && log2_tb_w <= ctu_dec->max_log2_transform_skip_size
                                               && log2_tb_h <= ctu_dec->max_log2_transform_skip_size) {
                VVCCABACContext *const cabac_ctx = ctu_dec->cabac_ctx;
                uint8_t transform_skip_flag = ovcabac_read_ae_transform_skip_flag_c(cabac_ctx);
            }

            #if 0
            ctu_dec->dequant_chroma = &ctu_dec->dequant_cr;
            #endif

            last_pos_cr = ovcabac_read_ae_last_sig_pos_c(cabac_ctx, log2_tb_w, log2_tb_h);

            lim_cg_w_cr = ((((last_pos_cr >> 8)) >> 2) + (((last_pos_cr & 0xFF))>> 2) + 1) << 2;

            last_pos_cr = ctu_dec->residual_coding_chroma(ctu_dec, coeffs_cr, log2_tb_w, log2_tb_h, last_pos_cr);

            /*FIXME avoid copy of lfnst_sb */
            memcpy(tmp_lfnst_cr, ctu_dec->lfnst_subblock, sizeof(int16_t) * 16);
        }

        if (ctu_dec->enable_lfnst) {
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
                VVCCABACContext *const cabac_ctx = ctu_dec->cabac_ctx;
                uint8_t is_dual = ctu_dec->transform_unit != transform_unit_st;
                lfnst_flag = ovcabac_read_ae_lfnst_flag(cabac_ctx, is_dual);
                if (lfnst_flag) {
                    lfnst_idx = ovcabac_read_ae_lfnst_idx(cabac_ctx);
                }
            }

        }


        #if 0
        if (cbf_mask & 0x2) {
            uint16_t *dst_cb = &ctu_dec->ctu_data_cb[VVC_CTB_OFFSET_CHROMA + x0 + VVC_CTB_STRIDE_CHROMA * y0];
            int16_t scale  = ctu_dec->lmcs_chroma_scale;
            int16_t *const coeffs_cb = ctu_dec->residual_cb;

            recon_residual_c(ctu_dec, ctu_dec->transform_buff, coeffs_cb, tmp_lfnst_cb, x0, y0, log2_tb_w, log2_tb_h,
                           lim_cg_w_cb, 0, 0, !last_pos_cb, lfnst_flag, 1, lfnst_idx);

            (*ctu_dec->scale_addsub_residuals)[0](ctu_dec->transform_buff, dst_cb, log2_tb_w, log2_tb_h, scale);

            fill_bs_map(&ctu_dec->dbf_info.bs1_map_cb, x0 << 1, y0 << 1, log2_tb_w + 1, log2_tb_h + 1);
        }

        if (cbf_mask & 0x1) {
            uint16_t *dst_cr  = &ctu_dec->ctu_data_cr[VVC_CTB_OFFSET_CHROMA + x0 + VVC_CTB_STRIDE_CHROMA * y0];
            int16_t scale  = ctu_dec->lmcs_chroma_scale;

            recon_residual_c(ctu_dec, ctu_dec->transform_buff, coeffs_cr, tmp_lfnst_cr, x0, y0, log2_tb_w, log2_tb_h,
                             lim_cg_w_cr, 0, 0, !last_pos_cr, lfnst_flag, 1, lfnst_idx);

            (*ctu_dec->scale_addsub_residuals)[0](ctu_dec->transform_buff, dst_cr, log2_tb_w, log2_tb_h, scale);

            fill_bs_map(&ctu_dec->dbf_info.bs1_map_cr, x0 << 1, y0 << 1, log2_tb_w + 1, log2_tb_h + 1);
        }
        #endif

    } else { //Joint cb cr (ICT)
        VVCCABACContext *const cabac_ctx = ctu_dec->cabac_ctx;
        uint16_t last_pos_cbcr;
        #if 0
        int shift_v = 6 + 1;
        int shift_h = (6 + 15 - 1) - 10;
        int16_t *dst_cb  = &ctu_dec->ctu_data_cb[VVC_CTB_OFFSET_CHROMA + x0 + VVC_CTB_STRIDE_CHROMA * y0];
        int16_t *dst_cr  = &ctu_dec->ctu_data_cr[VVC_CTB_OFFSET_CHROMA + x0 + VVC_CTB_STRIDE_CHROMA * y0];
        int16_t *src  = ctu_dec->transform_buff;
        int16_t scale = ctu_dec->lmcs_chroma_scale;
        #endif
        uint16_t *const coeffs_jcbcr = ctu_dec->residual_cb;
        uint8_t lfnst_flag = 0;
        uint8_t lfnst_idx;


        /* FIXME this is hackish joint cb cr involves a different delta qp from
          previous ones */

        #if 0
        if (cbf_mask == 3) {
            ctu_dec->dequant_chroma = &ctu_dec->dequant_joint_cb_cr;
        } else if (cbf_mask == 1) {
            ctu_dec->dequant_chroma = &ctu_dec->dequant_cr;
        } else {
            ctu_dec->dequant_chroma = &ctu_dec->dequant_cb;
        }
        #endif

        if (ctu_dec->transform_skip_enabled && log2_tb_w <= ctu_dec->max_log2_transform_skip_size
            && log2_tb_h <= ctu_dec->max_log2_transform_skip_size) {
            VVCCABACContext *const cabac_ctx = ctu_dec->cabac_ctx;
            uint8_t transform_skip_flag = ovcabac_read_ae_transform_skip_flag_c(cabac_ctx);
        }

        last_pos_cbcr = ovcabac_read_ae_last_sig_pos_c(cabac_ctx, log2_tb_w, log2_tb_h);

        int lim_cg_w_cbcr = ((((last_pos_cbcr >> 8)) >> 2) + (((last_pos_cbcr & 0xFF))>> 2) + 1) << 2;

        last_pos_cbcr = ctu_dec->residual_coding_chroma(ctu_dec, coeffs_jcbcr, log2_tb_w, log2_tb_h,
                                                       last_pos_cbcr);

        if (ctu_dec->enable_lfnst) {
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

#if 0
        recon_residual_c(ctu_dec, ctu_dec->transform_buff, coeffs_jcbcr, ctu_dec->lfnst_subblock, x0, y0, log2_tb_w, log2_tb_h,
                         lim_cg_w_cbcr, 0, 0, !last_pos_cbcr, lfnst_flag, 1, lfnst_idx);

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

        fill_bs_map(&ctu_dec->dbf_info.bs1_map_cb, x0 << 1, y0 << 1, log2_tb_w + 1, log2_tb_h + 1);
        fill_bs_map(&ctu_dec->dbf_info.bs1_map_cr, x0 << 1, y0 << 1, log2_tb_w + 1, log2_tb_h + 1);
#endif
    }
    return 0;
}

static int
transform_unit_st(OVCTUDec *const ctu_dec,
                  unsigned int x0, unsigned int y0,
                  unsigned int log2_tb_w, unsigned int log2_tb_h,
                  uint8_t rqt_root_cbf, uint8_t cu_flags)
{
    uint8_t cbf_mask = decode_cbf_st(ctu_dec, rqt_root_cbf);

    /* FIXME check if delta_qp is read per cu or per tu */
    if (cbf_mask) {
        if (ctu_dec->delta_qp_enabled && cbf_mask) {
            VVCCABACContext *const cabac_ctx = ctu_dec->cabac_ctx;
            int cu_qp_delta = ovcabac_read_ae_cu_delta_qp(cabac_ctx);
            #if 0
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

static int
transform_unit_l(OVCTUDec *const ctu_dec,
                  unsigned int x0, unsigned int y0,
                  unsigned int log2_tb_w, unsigned int log2_tb_h,
                  uint8_t rqt_root_cbf, uint8_t cu_flags)
{
    VVCCABACContext *const cabac_ctx = ctu_dec->cabac_ctx;
    uint8_t cbf_mask = ovcabac_read_ae_tu_cbf_luma(cabac_ctx);

    if (cbf_mask) {
        if (ctu_dec->delta_qp_enabled && cbf_mask) {
            int cu_qp_delta = ovcabac_read_ae_cu_delta_qp(cabac_ctx);
            #if 0
            derive_dequant_ctx(ctu_dec, &ctu_dec->qp_ctx, cu_qp_delta);
            #endif
        }

        transform_unit(ctu_dec, x0, y0, log2_tb_w, log2_tb_h, cbf_mask, cu_flags);
    }

    return 0;
}

static int
transform_unit_c(OVCTUDec *const ctu_dec,
                  unsigned int x0, unsigned int y0,
                  unsigned int log2_tb_w, unsigned int log2_tb_h,
                  uint8_t rqt_root_cbf, uint8_t cu_flags)
{
    VVCCABACContext *const cabac_ctx = ctu_dec->cabac_ctx;
    uint8_t cbf_mask = decode_cbf_c(ctu_dec);

    if (cbf_mask) {
        if (ctu_dec->delta_qp_enabled && cbf_mask) {
            int cu_qp_delta = ovcabac_read_ae_cu_delta_qp(cabac_ctx);
            #if 0
            derive_dequant_ctx(ctu_dec, &ctu_dec->qp_ctx, cu_qp_delta);
            #endif
        }

        transform_unit_chroma(ctu_dec, x0, y0, log2_tb_w, log2_tb_h, cbf_mask, cu_flags);
    }

    return 0;
}
