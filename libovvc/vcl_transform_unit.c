/**
 *
 *   OpenVVC is open-source real time software decoder compliant with the 
 *   ITU-T H.266- MPEG-I - Part 3 VVC standard. OpenVVC is developed from 
 *   scratch in C as a library that provides consumers with real time and
 *   energy-aware decoding capabilities under different OS including MAC OS,
 *   Windows, Linux and Android targeting low energy real-time decoding of
 *   4K VVC videos on Intel x86 and ARM platforms.
 * 
 *   Copyright (C) 2020-2022  IETR-INSA Rennes :
 *   
 *   Pierre-Loup CABARAT
 *   Wassim HAMIDOUCHE
 *   Guillaume GAUTIER
 *   Thomas AMESTOY
 *
 *   This library is free software; you can redistribute it and/or
 *   modify it under the terms of the GNU Lesser General Public
 *   License as published by the Free Software Foundation; either
 *   version 2.1 of the License, or (at your option) any later version.
 *
 *   This library is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *   Lesser General Public License for more details.
 *
 *   You should have received a copy of the GNU Lesser General Public
 *   License along with this library; if not, write to the Free Software
 *   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301
 *   USA
 * 
 **/

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

#define LOG2_MIN_CU_S 2
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

static uint8_t
ovcabac_read_ae_root_cbf(OVCABACCtx *const cabac_ctx)
{
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    return ovcabac_ae_read(cabac_ctx, &cabac_state[QT_ROOT_CBF_CTX_OFFSET]);
}

static uint8_t
ovcabac_read_ae_tu_cbf_luma_isp(OVCABACCtx *const cabac_ctx,
                                uint8_t prev_cbf)
{
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    int offset = 2 + prev_cbf;
    return ovcabac_ae_read(cabac_ctx, &cabac_state[QT_CBF_CTX_OFFSET + offset] );
}

static uint8_t
ovcabac_read_ae_tu_cbf_luma(OVCABACCtx *const cabac_ctx, uint8_t bdpcm_flag)
{
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    return ovcabac_ae_read(cabac_ctx, &cabac_state[QT_CBF_CTX_OFFSET + bdpcm_flag]);
}

static uint8_t
ovcabac_read_ae_tu_cbf_cb(OVCABACCtx *const cabac_ctx, uint8_t bdpcm_flag)
{
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    return ovcabac_ae_read(cabac_ctx, &cabac_state[QT_CBF_CB_CTX_OFFSET + bdpcm_flag] );
}

static uint8_t
ovcabac_read_ae_tu_cbf_cr(OVCABACCtx *const cabac_ctx,
                          uint8_t tu_cbf_cb, uint8_t bdpcm_flag)
{
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    return ovcabac_ae_read(cabac_ctx, &cabac_state[QT_CBF_CR_CTX_OFFSET + (tu_cbf_cb | bdpcm_flag) + bdpcm_flag] );
}

static uint8_t
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

static int
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

static int
ovcabac_read_ae_cu_chroma_qp_offset(OVCABACCtx *const cabac_ctx, uint8_t length)
{
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    int delta_qp_idx = ovcabac_ae_read(cabac_ctx, &cabac_state[CHROMA_QP_ADJ_FLAG_CTX_OFFSET]);
    if (delta_qp_idx)
        while (delta_qp_idx <= length - 1 && ovcabac_ae_read(cabac_ctx, &cabac_state[CHROMA_QP_ADJ_IDC_CTX_OFFSET]))
            delta_qp_idx++;

    return delta_qp_idx;
}

static uint8_t
ovcabac_read_ae_cu_mts_flag(OVCABACCtx *const cabac_ctx)
{
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    return ovcabac_ae_read(cabac_ctx,&cabac_state[MTS_IDX_CTX_OFFSET]);
}

static uint8_t
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

static uint8_t
ovcabac_read_ae_transform_skip_luma_flag(OVCABACCtx *const cabac_ctx)
{
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    return ovcabac_ae_read(cabac_ctx, &cabac_state[TRANSFORM_SKIP_FLAG_CTX_OFFSET]);
}

static uint8_t
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
decode_last_sig_prefix_sbt_mts(OVCABACCtx *const cabac_ctx,
                               unsigned int log2_tb_d, unsigned int log2_tb_red,
                               unsigned offset_ctx)
{
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    /* FIXME this tab could be adapted by adding some unused ctx to last sig ctx */
    static const int prefix_ctx[8]  = { 0, 0, 0, 3, 6, 10, 15, 21 };
    int pos = 0;
    int ctx_offset, ctx_shift;
    int max_symbol = OVMIN(log2_tb_red, 5) << 1;

    ctx_offset = prefix_ctx[log2_tb_d];
    ctx_shift  = (log2_tb_red + 1) >> 2 ;

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

static uint16_t
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

static uint16_t
ovcabac_read_ae_last_sig_pos_red(OVCABACCtx *const cabac_ctx,
                                 uint8_t log2_tb_w, uint8_t log2_tb_h)
{
    uint8_t last_x;
    uint8_t last_y;
    uint8_t log2_red_w = log2_tb_w == 5 ? 4 : log2_tb_w;
    uint8_t log2_red_h = log2_tb_h == 5 ? 4 : log2_tb_h;

    last_x = decode_last_sig_prefix_sbt_mts(cabac_ctx, log2_tb_w, log2_red_w, LAST_X_CTX_OFFSET);

    last_y = decode_last_sig_prefix_sbt_mts(cabac_ctx, log2_tb_h, log2_red_h, LAST_Y_CTX_OFFSET);

    if (last_x > 3) {
        last_x = decode_last_sig_suffix(cabac_ctx, last_x);
    }

    if (last_y > 3) {
        last_y = decode_last_sig_suffix(cabac_ctx, last_y);
    }

    return ((uint16_t) last_y << 8) | (last_x & 0xFF);
}

static uint16_t
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

static uint8_t
ovcabac_read_ae_lfnst_flag(OVCABACCtx *const cabac_ctx, uint8_t is_dual_tree)
{
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    uint8_t lfnst_flag = ovcabac_ae_read(cabac_ctx, &cabac_state[LFNST_IDX_CTX_OFFSET + is_dual_tree]);
    return lfnst_flag;
}

static uint8_t
ovcabac_read_ae_lfnst_idx(OVCABACCtx *const cabac_ctx)
{
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    return ovcabac_ae_read(cabac_ctx, &cabac_state[LFNST_IDX_CTX_OFFSET + 2]);
}

static uint8_t
ovcabac_read_ae_sbt_flag(OVCABACCtx *const cabac_ctx, uint8_t log2_tb_w, uint8_t log2_tb_h)
{
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    uint8_t ctx_offset = log2_tb_w + log2_tb_h <= 8;
    return ovcabac_ae_read(cabac_ctx, &cabac_state[SBT_FLAG_CTX_OFFSET + ctx_offset]);
}

static uint8_t
ovcabac_read_ae_sbt_quad_flag(OVCABACCtx *const cabac_ctx)
{
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    return ovcabac_ae_read(cabac_ctx, &cabac_state[SBT_QUAD_FLAG_CTX_OFFSET]);
}

static uint8_t
ovcabac_read_ae_sbt_hor_flag(OVCABACCtx *const cabac_ctx, uint8_t log2_tb_w, uint8_t log2_tb_h)
{
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    uint8_t ctx_offset = (log2_tb_w == log2_tb_h) ? 0 :  (log2_tb_w < log2_tb_h) ? 1 : 2;
    return ovcabac_ae_read(cabac_ctx, &cabac_state[SBT_HOR_FLAG_CTX_OFFSET + ctx_offset]);
}

static uint8_t
ovcabac_read_ae_sbt_pos_flag(OVCABACCtx *const cabac_ctx)
{
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    return ovcabac_ae_read(cabac_ctx, &cabac_state[SBT_POS_FLAG_CTX_OFFSET]);
}


static uint8_t
decode_cbf_st(OVCTUDec *const ctu_dec, uint8_t rqt_root_cbf, uint8_t tr_depth, CUFlags cu_flags)
{
    OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;
    uint8_t intra_bdpcm_chroma_flag = !!(cu_flags & flg_intra_bdpcm_chroma_flag);
    uint8_t tu_cbf_cb = ovcabac_read_ae_tu_cbf_cb(cabac_ctx, intra_bdpcm_chroma_flag);
    uint8_t tu_cbf_cr = ovcabac_read_ae_tu_cbf_cr(cabac_ctx, tu_cbf_cb, intra_bdpcm_chroma_flag);
    uint8_t cbf_mask = (tu_cbf_cb << 1) | tu_cbf_cr;
    uint8_t tu_cbf_luma = rqt_root_cbf;

    if (!rqt_root_cbf || (cbf_mask && rqt_root_cbf) || tr_depth){
        uint8_t intra_bdpcm_luma_flag = !!(cu_flags & flg_intra_bdpcm_luma_flag);
        tu_cbf_luma = ovcabac_read_ae_tu_cbf_luma(cabac_ctx, intra_bdpcm_luma_flag);
    }

    /* FIXME delta_qp is only read on first significant TU in CU */
    if (ctu_dec->delta_qp_enabled && (rqt_root_cbf | cbf_mask | tu_cbf_luma) && ctu_dec->read_qp) {
        OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;
        int cu_qp_delta = ovcabac_read_ae_cu_delta_qp(cabac_ctx);
        derive_dequant_ctx(ctu_dec, &ctu_dec->qp_ctx, cu_qp_delta);
        ctu_dec->read_qp = 0;
    }

    #if 1
    if (ctu_dec->chroma_qp_offset_enabled && cbf_mask && ctu_dec->read_qp_c) {
        OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;
        int qp_bd_offset = ctu_dec->qp_ctx.qp_bd_offset;
        uint8_t length = ctu_dec->chroma_qp_offset_len;
        int cu_qp_delta = ovcabac_read_ae_cu_chroma_qp_offset(cabac_ctx, length);
        derive_dequant_ctx(ctu_dec, &ctu_dec->qp_ctx, 0);
        if (cu_qp_delta) {
        --cu_qp_delta;
        ctu_dec->qp_ctx.dqp_cb = ctu_dec->qp_ctx.cb_qp_offset_list[cu_qp_delta];
        ctu_dec->qp_ctx.dqp_cr = ctu_dec->qp_ctx.cr_qp_offset_list[cu_qp_delta];
        ctu_dec->qp_ctx.dqp_jcbcr = ctu_dec->qp_ctx.joint_cbcr_qp_offset_list[cu_qp_delta];
        ctu_dec->dequant_cb.qp = ov_clip(ctu_dec->dequant_cb.qp - qp_bd_offset + ctu_dec->qp_ctx.cb_qp_offset_list[cu_qp_delta], -qp_bd_offset, 63) + qp_bd_offset;
        ctu_dec->dequant_cr.qp = ov_clip(ctu_dec->dequant_cr.qp - qp_bd_offset + ctu_dec->qp_ctx.cr_qp_offset_list[cu_qp_delta], -qp_bd_offset, 63) + qp_bd_offset;
        ctu_dec->dequant_joint_cb_cr.qp = ov_clip(ctu_dec->dequant_joint_cb_cr.qp - qp_bd_offset + ctu_dec->qp_ctx.joint_cbcr_qp_offset_list[cu_qp_delta], -qp_bd_offset, 63) + qp_bd_offset;
        ctu_dec->dequant_cb_skip.qp   = OVMAX(ctu_dec->dequant_cb.qp, ctu_dec->qp_ctx.min_qp_prime_ts);
        ctu_dec->dequant_cr_skip.qp   = OVMAX(ctu_dec->dequant_cr.qp, ctu_dec->qp_ctx.min_qp_prime_ts);
        ctu_dec->dequant_jcbcr_skip.qp = OVMAX(ctu_dec->dequant_joint_cb_cr.qp, ctu_dec->qp_ctx.min_qp_prime_ts);
        }
        ctu_dec->read_qp_c = 0;
    }
    #endif

    /* FIXME intra if inter we only check for cbf_mask == 3*/
    if (ctu_dec->jcbcr_enabled && (!(cu_flags & flg_ibc_flag) && ((cu_flags & flg_pred_mode_flag) && cbf_mask) || cbf_mask == 3)) {
        uint8_t joint_cb_cr = ovcabac_read_ae_joint_cb_cr_flag(cabac_ctx, (cbf_mask & 0x3) - 1);
        cbf_mask |= joint_cb_cr << 3;
    }

    return cbf_mask | (tu_cbf_luma << 4);
}

static uint8_t
decode_cbf_c(OVCTUDec *const ctu_dec, CUFlags cu_flags)
{
    OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;
    uint8_t intra_bdpcm_chroma_flag = !!(cu_flags & flg_intra_bdpcm_chroma_flag);
    uint8_t tu_cbf_cb = ovcabac_read_ae_tu_cbf_cb(cabac_ctx, intra_bdpcm_chroma_flag);
    uint8_t tu_cbf_cr = ovcabac_read_ae_tu_cbf_cr(cabac_ctx, tu_cbf_cb, intra_bdpcm_chroma_flag);
    uint8_t cbf_mask = (tu_cbf_cb << 1) | tu_cbf_cr;

    if (ctu_dec->chroma_qp_offset_enabled && cbf_mask && ctu_dec->read_qp_c) {
        OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;
        int qp_bd_offset = ctu_dec->qp_ctx.qp_bd_offset;
        uint8_t length = ctu_dec->chroma_qp_offset_len;
        int cu_qp_delta = ovcabac_read_ae_cu_chroma_qp_offset(cabac_ctx, length);
        derive_dequant_ctx(ctu_dec, &ctu_dec->qp_ctx, 0);
        if (cu_qp_delta) {
        --cu_qp_delta;
        ctu_dec->qp_ctx.dqp_cb = ctu_dec->qp_ctx.cb_qp_offset_list[cu_qp_delta];
        ctu_dec->qp_ctx.dqp_cr = ctu_dec->qp_ctx.cr_qp_offset_list[cu_qp_delta];
        ctu_dec->qp_ctx.dqp_jcbcr = ctu_dec->qp_ctx.joint_cbcr_qp_offset_list[cu_qp_delta];
        ctu_dec->dequant_cb.qp = ov_clip(ctu_dec->dequant_cb.qp - qp_bd_offset + ctu_dec->qp_ctx.cb_qp_offset_list[cu_qp_delta], -qp_bd_offset, 63) + qp_bd_offset;
        ctu_dec->dequant_cr.qp = ov_clip(ctu_dec->dequant_cr.qp - qp_bd_offset + ctu_dec->qp_ctx.cr_qp_offset_list[cu_qp_delta], -qp_bd_offset, 63) + qp_bd_offset;
        ctu_dec->dequant_joint_cb_cr.qp = ov_clip(ctu_dec->dequant_joint_cb_cr.qp - qp_bd_offset + ctu_dec->qp_ctx.joint_cbcr_qp_offset_list[cu_qp_delta], -qp_bd_offset, 63) + qp_bd_offset;
        ctu_dec->dequant_cb_skip.qp   = OVMAX(ctu_dec->dequant_cb.qp, ctu_dec->qp_ctx.min_qp_prime_ts);
        ctu_dec->dequant_cr_skip.qp   = OVMAX(ctu_dec->dequant_cr.qp, ctu_dec->qp_ctx.min_qp_prime_ts);
        ctu_dec->dequant_jcbcr_skip.qp = OVMAX(ctu_dec->dequant_joint_cb_cr.qp, ctu_dec->qp_ctx.min_qp_prime_ts);
        }
        ctu_dec->read_qp_c = 0;
    }

    if (ctu_dec->jcbcr_enabled && cbf_mask) {
        uint8_t joint_cb_cr = ovcabac_read_ae_joint_cb_cr_flag(cabac_ctx,
                (cbf_mask & 0x3) - 1);
        cbf_mask |= joint_cb_cr << 3;
    }

    return cbf_mask;
}

static inline uint8_t
check_lfnst_nb_coeffs(uint16_t last_pos)
{
    static const uint64_t scan_map = 0xFDA6EB73C8419520;

    int last_y = last_pos >> 8;
    int last_x = last_pos & 0xFF;
    uint8_t ret_val = -!!((last_x >> 2) | (last_y >> 2));
    uint8_t last_sb_pos = ((last_x & 0x3) + ((last_y & 0x3) << 2));
    uint8_t nb_coeffs = (scan_map >> (last_sb_pos << 2)) & 0xF;

    return nb_coeffs | ret_val;
}

static uint8_t
jcbcr_lfnst_check(const struct TUInfo *tu_info, uint8_t log2_tb_w, uint8_t log2_tb_h)
{
    const struct TBInfo *const tb_info = &tu_info->tb_info[0];
    if (!tu_info->tr_skip_mask && tb_info->sig_sb_map <= 0x1 && log2_tb_h > 1 && log2_tb_w > 1) {
        int max_lfnst_pos = (log2_tb_h == log2_tb_w) && (log2_tb_w <= 3) ? 7 : 15;
        int nb_coeffs = check_lfnst_nb_coeffs(tb_info->last_pos);

        uint8_t can_lfnst = nb_coeffs <= max_lfnst_pos;

        can_lfnst &= !!nb_coeffs;

        return can_lfnst;
    }
    return 0;
}

static uint8_t
chroma_lfnst_check(const struct TUInfo *tu_info, uint8_t cbf_mask, uint8_t log2_tb_w, uint8_t log2_tb_h)
{
    uint8_t can_lfnst = 0;
    const struct TBInfo *const tb_info_cb = &tu_info->tb_info[0];
    const struct TBInfo *const tb_info_cr = &tu_info->tb_info[1];

    if (log2_tb_w > 1 && log2_tb_h > 1 && tb_info_cb->sig_sb_map <= 0x1 && tb_info_cr->sig_sb_map <= 0x1 && !tu_info->tr_skip_mask) {
        uint8_t need_cb_chk = cbf_mask & 0x2;
        uint8_t need_cr_chk = cbf_mask & 0x1;

        int max_lfnst_pos = (log2_tb_h == log2_tb_w) && (log2_tb_w <= 3) ? 7 : 15;

        can_lfnst = !!cbf_mask;

        /*fixme better can_lfnst derivation */
        if (need_cb_chk && need_cr_chk) {
            int nb_coeffs_cb = check_lfnst_nb_coeffs(tb_info_cb->last_pos);
            int nb_coeffs_cr = check_lfnst_nb_coeffs(tb_info_cr->last_pos);
            can_lfnst &= nb_coeffs_cr <= max_lfnst_pos && nb_coeffs_cb <= max_lfnst_pos;
            can_lfnst &= !!(nb_coeffs_cr | nb_coeffs_cb);
        } else if (need_cb_chk) {
            int nb_coeffs_cb = check_lfnst_nb_coeffs(tb_info_cb->last_pos);
            can_lfnst &= nb_coeffs_cb <= max_lfnst_pos;
            can_lfnst &= !!nb_coeffs_cb;
        } else {
            int nb_coeffs_cr = check_lfnst_nb_coeffs(tb_info_cr->last_pos);
            can_lfnst &= nb_coeffs_cr <= max_lfnst_pos;
            can_lfnst &= !!nb_coeffs_cr;
        }
    }

    return can_lfnst;
}

static uint8_t
lfnst_check_st(const struct TUInfo *const tu_info, uint8_t log2_tb_w, uint8_t log2_tb_h,
               uint8_t cbf_mask, CUFlags cu_flags)
{
    uint8_t cbf_flag_l  = cbf_mask & 0x10;
    uint8_t jcbcr_flag  = cbf_mask & 0x8;
    uint8_t cbf_mask_cb = cbf_mask & 0x2;
    uint8_t cbf_mask_cr = cbf_mask & 0x1;
    uint8_t non_only_dc = 0;
    /* FIXME chroma size ? */
    const uint8_t max_lfnst_pos   = (log2_tb_h == log2_tb_w) && (log2_tb_w <= 3) ? 7 : 15;
    const uint8_t max_lfnst_pos_c = (log2_tb_h == log2_tb_w) && (log2_tb_w <= 4) ? 7 : 15;

    /*FIXME MIP check on ly for luma TB*/
    uint8_t is_mip = !!(cu_flags & flg_mip_flag);
    uint8_t mip_lfnst = !is_mip || (log2_tb_h >= 4 && log2_tb_w >= 4);

    uint8_t can_lfnst = mip_lfnst && !(cu_flags & flg_ibc_flag);

    can_lfnst &= !(tu_info->tr_skip_mask);

    /* FIXME Note that sig_sb_map check is sufficient if max is 15 */
    if (can_lfnst) {
        if (cbf_flag_l) {
            const struct TBInfo *tb_info = &tu_info->tb_info[2];
            uint8_t nb_coeffs = check_lfnst_nb_coeffs(tb_info->last_pos);
            can_lfnst &= tb_info->sig_sb_map <= 0x1;
            can_lfnst &= nb_coeffs <= max_lfnst_pos;
            non_only_dc |= nb_coeffs;
        }

        if (jcbcr_flag && log2_tb_h > 2 && log2_tb_w > 2) {
            const struct TBInfo *const tb_info_cbcr = &tu_info->tb_info[0];
            uint8_t nb_coeffs_cbcr = check_lfnst_nb_coeffs(tb_info_cbcr->last_pos);
            can_lfnst &= tb_info_cbcr->sig_sb_map <= 0x1;
            can_lfnst &= nb_coeffs_cbcr <= max_lfnst_pos_c;
            non_only_dc |= nb_coeffs_cbcr;
        } else {

            if (cbf_mask_cb && log2_tb_h > 2 && log2_tb_w > 2) {
                const struct TBInfo *const tb_info_cb = &tu_info->tb_info[0];
                uint8_t nb_coeffs_cb = check_lfnst_nb_coeffs(tb_info_cb->last_pos);
                can_lfnst &= tb_info_cb->sig_sb_map <= 0x1;
                can_lfnst &= nb_coeffs_cb <= max_lfnst_pos_c;
                non_only_dc |= nb_coeffs_cb;
            }

            if (cbf_mask_cr && log2_tb_h > 2 && log2_tb_w > 2) {
                const struct TBInfo *const tb_info_cr = &tu_info->tb_info[1];
                uint8_t nb_coeffs_cr = check_lfnst_nb_coeffs(tb_info_cr->last_pos);
                can_lfnst &= tb_info_cr->sig_sb_map <= 0x1;
                can_lfnst &= nb_coeffs_cr <= max_lfnst_pos_c;
                non_only_dc |= nb_coeffs_cr;
            }
        }

        can_lfnst &= !!non_only_dc;

    }

    return can_lfnst;
}

static int
residual_coding_l(OVCTUDec *const ctu_dec,
                  unsigned int x0, unsigned int y0,
                  unsigned int log2_tb_w, unsigned int log2_tb_h,
                  CUFlags cu_flags, struct TUInfo *const tu_info)
{
    OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;
    uint8_t tr_skip_flag = 0;
    int16_t *const coeffs_y = ctu_dec->residual_y + tu_info->pos_offset;

    struct TBInfo *tb_info = &tu_info->tb_info[2];

    if (ctu_dec->transform_skip_enabled && log2_tb_w <= ctu_dec->max_log2_transform_skip_size
        && log2_tb_h <= ctu_dec->max_log2_transform_skip_size && !tu_info->is_sbt && !(cu_flags &flg_isp_flag)) {
        ctu_dec->dequant_skip = &ctu_dec->dequant_luma_skip;
        uint8_t intra_bdpcm_luma_flag = !!(cu_flags & flg_intra_bdpcm_luma_flag);
        tr_skip_flag = intra_bdpcm_luma_flag || ovcabac_read_ae_transform_skip_luma_flag(cabac_ctx);
        tu_info->tr_skip_mask |= tr_skip_flag << 4;
    }


    if (!tr_skip_flag || ctu_dec->sh_ts_disabled) {
        uint16_t last_pos;
        if (tu_info->is_sbt && tu_info->cu_mts_flag && log2_tb_w <= 5 && log2_tb_h <= 5) {
            last_pos = ovcabac_read_ae_last_sig_pos_red(cabac_ctx, log2_tb_w, log2_tb_h);
            uint8_t log2_red_w = log2_tb_w == 5 ? 4 : log2_tb_w;
            uint8_t log2_red_h = log2_tb_h == 5 ? 4 : log2_tb_h;

            /* FIXME recheck this */
            ctu_dec->tmp_red  =  log2_tb_w - log2_red_w;
            ctu_dec->tmp_red |= (log2_tb_h - log2_red_h) << 1;

            memset(coeffs_y, 0, sizeof(int16_t) << (log2_tb_h + log2_tb_w));
        } else {
            last_pos = ovcabac_read_ae_last_sig_pos(cabac_ctx, log2_tb_w, log2_tb_h);
        }

        uint64_t sig_sb_map = ctu_dec->residual_coding_l(ctu_dec, coeffs_y, log2_tb_w, log2_tb_h,
                                                         last_pos);

        ctu_dec->tmp_red  = 0;
        tb_info->sig_sb_map = sig_sb_map;
        tb_info->last_pos   = last_pos;


    } else {
        uint8_t intra_bdpcm_luma_flag = !!(cu_flags & flg_intra_bdpcm_luma_flag);
        ctu_dec->dequant_skip = &ctu_dec->dequant_luma_skip;
        residual_coding_ts(ctu_dec, coeffs_y, log2_tb_w, log2_tb_h,
                           intra_bdpcm_luma_flag);
    }


    return 0;
}

static int
residual_coding_c(OVCTUDec *const ctu_dec,
                  unsigned int x0, unsigned int y0,
                  unsigned int log2_tb_w, unsigned int log2_tb_h,
                  uint8_t cbf_mask, CUFlags cu_flags, struct TUInfo *tu_info)
{

    int16_t *const coeffs_cb = ctu_dec->residual_cb + tu_info->pos_offset;
    int16_t *const coeffs_cr = ctu_dec->residual_cr + tu_info->pos_offset;
    uint16_t last_pos_cb = (1 << (log2_tb_h + log2_tb_w)) - 1;
    uint16_t last_pos_cr = (1 << (log2_tb_h + log2_tb_w)) - 1;
    uint64_t sig_sb_map_cb = 0x1;
    uint64_t sig_sb_map_cr = 0x1;
    uint8_t transform_skip_flag = 0;
    struct TBInfo *const tb_info_cb = &tu_info->tb_info[0];
    struct TBInfo *const tb_info_cr = &tu_info->tb_info[1];


    /* FIXME move dequant to reconstruction this require modification in coeff
       reading */

    if (cbf_mask & 0x2) {
        OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;

        /* TODO use max_tr_skip_s = 0 to avoid using enabled test */
        if (ctu_dec->transform_skip_enabled && log2_tb_w <= ctu_dec->max_log2_transform_skip_size
            && log2_tb_h <= ctu_dec->max_log2_transform_skip_size && !tu_info->is_sbt) {
            OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;
            uint8_t intra_bdpcm_chroma_flag = !!(cu_flags & flg_intra_bdpcm_chroma_flag);
            ctu_dec->dequant_skip = &ctu_dec->dequant_cb_skip;
            transform_skip_flag |= (intra_bdpcm_chroma_flag || ovcabac_read_ae_transform_skip_flag_c(cabac_ctx)) << 1;
        }

        if (!(transform_skip_flag & 0x2) || ctu_dec->sh_ts_disabled) {

            last_pos_cb = ovcabac_read_ae_last_sig_pos_c(cabac_ctx, log2_tb_w, log2_tb_h);

            sig_sb_map_cb = ctu_dec->residual_coding_c(ctu_dec, coeffs_cb,
                                                       log2_tb_w, log2_tb_h,
                                                       last_pos_cb);

        }  else {
            /* FIXME Chroma TS */
            uint8_t intra_bdpcm_chroma_flag = !!(cu_flags & flg_intra_bdpcm_chroma_flag);
            ctu_dec->dequant_skip = &ctu_dec->dequant_cb_skip;
            residual_coding_ts(ctu_dec, coeffs_cb, log2_tb_w, log2_tb_h, intra_bdpcm_chroma_flag);
        }
    }

    if (cbf_mask & 0x1) {
        OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;

        /* TODO use max_tr_skip_s = 0 to avoid using enabled test */
        if (ctu_dec->transform_skip_enabled && log2_tb_w <= ctu_dec->max_log2_transform_skip_size
            && log2_tb_h <= ctu_dec->max_log2_transform_skip_size && !tu_info->is_sbt) {
            OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;
            uint8_t intra_bdpcm_chroma_flag = !!(cu_flags & flg_intra_bdpcm_chroma_flag);
            ctu_dec->dequant_skip = &ctu_dec->dequant_cr_skip;
            transform_skip_flag |= intra_bdpcm_chroma_flag || ovcabac_read_ae_transform_skip_flag_c(cabac_ctx);
        }

        if (!(transform_skip_flag & 0x1) || ctu_dec->sh_ts_disabled) {

            last_pos_cr = ovcabac_read_ae_last_sig_pos_c(cabac_ctx, log2_tb_w, log2_tb_h);

            sig_sb_map_cr = ctu_dec->residual_coding_c(ctu_dec, coeffs_cr,
                                                       log2_tb_w, log2_tb_h,
                                                       last_pos_cr);
        } else {
            uint8_t intra_bdpcm_chroma_flag = !!(cu_flags & flg_intra_bdpcm_chroma_flag);
            ctu_dec->dequant_skip = &ctu_dec->dequant_cr_skip;
            residual_coding_ts(ctu_dec, coeffs_cr, log2_tb_w, log2_tb_h, intra_bdpcm_chroma_flag);
        }
    }

    tu_info->tr_skip_mask  |= transform_skip_flag;
    tb_info_cb->last_pos   = last_pos_cb;
    tb_info_cr->last_pos   = last_pos_cr;
    tb_info_cb->sig_sb_map = sig_sb_map_cb;
    tb_info_cr->sig_sb_map = sig_sb_map_cr;

    return 0;
}

static int
residual_coding_jcbcr(OVCTUDec *const ctu_dec,
                      unsigned int x0, unsigned int y0,
                      unsigned int log2_tb_w, unsigned int log2_tb_h,
                      uint8_t cbf_mask, CUFlags cu_flags, struct TUInfo *tu_info)
{
    OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;

    uint16_t last_pos;
    int16_t *const coeffs_jcbcr = ctu_dec->residual_cb + tu_info->pos_offset;
    uint64_t sig_sb_map = 0x1;
    uint8_t transform_skip_flag = 0;
    struct TBInfo *const tb_info = &tu_info->tb_info[0];

    if (ctu_dec->transform_skip_enabled && log2_tb_w <= ctu_dec->max_log2_transform_skip_size
        && log2_tb_h <= ctu_dec->max_log2_transform_skip_size && !tu_info->is_sbt) {
        OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;
        uint8_t intra_bdpcm_chroma_flag = !!(cu_flags & flg_intra_bdpcm_chroma_flag);
        transform_skip_flag |= intra_bdpcm_chroma_flag || ovcabac_read_ae_transform_skip_flag_c(cabac_ctx);
        tu_info->tr_skip_mask |= transform_skip_flag;
    }

    if (!transform_skip_flag || ctu_dec->sh_ts_disabled) {
        last_pos = ovcabac_read_ae_last_sig_pos_c(cabac_ctx, log2_tb_w, log2_tb_h);

        sig_sb_map = ctu_dec->residual_coding_c(ctu_dec, coeffs_jcbcr, log2_tb_w, log2_tb_h,
                                                last_pos);
        tb_info->last_pos   = last_pos;
    } else {
        uint8_t intra_bdpcm_chroma_flag = !!(cu_flags & flg_intra_bdpcm_chroma_flag);

        /* FIXME this is hackish joint cb cr involves a different delta qp from
           previous ones */
        if (cbf_mask == 3) {
            ctu_dec->dequant_skip = &ctu_dec->dequant_jcbcr_skip;
        } else if (cbf_mask == 1) {
            ctu_dec->dequant_skip   = &ctu_dec->dequant_cr_skip;
        } else {
            ctu_dec->dequant_skip   = &ctu_dec->dequant_cb_skip;
        }

        residual_coding_ts(ctu_dec, ctu_dec->residual_cb + tu_info->pos_offset, log2_tb_w, log2_tb_h,
                           intra_bdpcm_chroma_flag);
    }

    tb_info->sig_sb_map = sig_sb_map;

    return 0;
}

int
transform_unit_st(OVCTUDec *const ctu_dec,
                  unsigned int x0, unsigned int y0,
                  unsigned int log2_tb_w, unsigned int log2_tb_h,
                  uint8_t rqt_root_cbf, CUFlags cu_flags, uint8_t tr_depth,
                  struct TUInfo *const tu_info)
{
    uint8_t cbf_mask = decode_cbf_st(ctu_dec, rqt_root_cbf, tr_depth, cu_flags);

    uint8_t cbf_flag_l = cbf_mask & 0x10;
    uint8_t jcbcr_flag = cbf_mask & 0x8;
    uint8_t cbf_mask_c = cbf_mask & 0x3;

    if (cbf_mask) {

        if (cbf_flag_l) {

            residual_coding_l(ctu_dec, x0, y0, log2_tb_w, log2_tb_h, cu_flags, tu_info);

        }

        if (jcbcr_flag) {

            residual_coding_jcbcr(ctu_dec, x0 >> 1, y0 >> 1, log2_tb_w - 1,
                                  log2_tb_h - 1, cbf_mask_c, cu_flags, tu_info);

        } else if (cbf_mask_c) {

            residual_coding_c(ctu_dec, x0 >> 1, y0 >> 1, log2_tb_w - 1,
                              log2_tb_h - 1, cbf_mask_c, cu_flags, tu_info);

        }

    }
    return cbf_mask;
}

int
transform_unit_l(struct OVCTUDec *const ctu_dec,
                 unsigned int x0, unsigned int y0,
                 unsigned int log2_tb_w, unsigned int log2_tb_h,
                 uint8_t rqt_root_cbf, CUFlags cu_flags, uint8_t tr_depth,
                 struct TUInfo *const tu_info)
{
    OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;
    uint8_t intra_bdpcm_luma_flag = !!(cu_flags & flg_intra_bdpcm_luma_flag);
    uint8_t cbf_mask = rqt_root_cbf || ovcabac_read_ae_tu_cbf_luma(cabac_ctx, intra_bdpcm_luma_flag);

    if (cbf_mask) {
        if (ctu_dec->delta_qp_enabled && cbf_mask && ctu_dec->read_qp) {
            int cu_qp_delta = ovcabac_read_ae_cu_delta_qp(cabac_ctx);
            derive_dequant_ctx(ctu_dec, &ctu_dec->qp_ctx, cu_qp_delta);
            ctu_dec->read_qp = 0;
        }

        residual_coding_l(ctu_dec, x0, y0, log2_tb_w, log2_tb_h, cu_flags, tu_info);

    }

    return cbf_mask;
}

int
transform_unit_c(OVCTUDec *const ctu_dec,
                 unsigned int x0, unsigned int y0,
                 unsigned int log2_tb_w, unsigned int log2_tb_h,
                 uint8_t rqt_root_cbf, CUFlags cu_flags, uint8_t tr_depth,
                 struct TUInfo *const tu_info)
{
    OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;
    uint8_t cbf_mask = decode_cbf_c(ctu_dec, cu_flags);
    uint8_t jcbcr_flag = cbf_mask & 0x8;
    uint8_t cbf_mask_c = cbf_mask & 0x3;

    if (cbf_mask) {

        if (jcbcr_flag) {

            residual_coding_jcbcr(ctu_dec, x0, y0, log2_tb_w, log2_tb_h,
                                  cbf_mask_c, cu_flags, tu_info);
        } else if (cbf_mask_c) {

            residual_coding_c(ctu_dec, x0, y0, log2_tb_w, log2_tb_h,
                              cbf_mask_c, cu_flags, tu_info);
        }
    }

    return cbf_mask;
}

static void
lfnst_mts(const OVCTUDec *const ctu_dec, uint8_t log2_tb_w, uint8_t log2_tb_h,
          CUFlags cu_flags, struct TUInfo *const tu_info)
{
    uint8_t cbf_mask = tu_info->cbf_mask;
    if (ctu_dec->transform_unit == &transform_unit_st && cbf_mask) {
                uint8_t is_mip = (cu_flags & flg_pred_mode_flag) && !!(cu_flags & flg_mip_flag);
                uint8_t allow_mip_lfnst = !is_mip || (log2_tb_h >= 4 && log2_tb_w >= 4);

        if (!(tu_info->tr_skip_mask))
        if ((cu_flags & flg_pred_mode_flag) && allow_mip_lfnst && ctu_dec->enable_lfnst && cu_flags & flg_pred_mode_flag && !(cu_flags & flg_ibc_flag)) {
            uint8_t can_lfnst = lfnst_check_st(tu_info, log2_tb_w, log2_tb_h,
                                               cbf_mask, cu_flags);

            if (can_lfnst) {
                OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;
                uint8_t is_dual = ctu_dec->transform_unit != &transform_unit_st;
                uint8_t lfnst_flag = ovcabac_read_ae_lfnst_flag(cabac_ctx, is_dual);
                tu_info->lfnst_flag = lfnst_flag;
                if (lfnst_flag) {
                    uint8_t lfnst_idx = ovcabac_read_ae_lfnst_idx(cabac_ctx);
                    tu_info->lfnst_idx = lfnst_idx;
                }
            }
        }

        if (!(tu_info->tr_skip_mask & 0x10)) {
        if (ctu_dec->mts_enabled && (((cu_flags & flg_pred_mode_flag) && ctu_dec->mts_explicit_intra)
                                 || (!(cu_flags & flg_pred_mode_flag) && ctu_dec->mts_explicit_inter))) {
            const struct TBInfo *tb_info = &tu_info->tb_info[2];
            if (!tu_info->lfnst_flag && !!tb_info->last_pos && (log2_tb_w < 6) && (log2_tb_h < 6)
                && !(tb_info->sig_sb_map & (~0x000000000F0F0F0F))) {
                OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;
                uint8_t cu_mts_flag = ovcabac_read_ae_cu_mts_flag(cabac_ctx);
                tu_info->cu_mts_flag = cu_mts_flag;
                if (cu_mts_flag) {
                    uint8_t cu_mts_idx = ovcabac_read_ae_cu_mts_idx(cabac_ctx);
                    tu_info->cu_mts_idx = cu_mts_idx;
                }
            }
        }
        }
    } else if (ctu_dec->transform_unit == &transform_unit_l && cbf_mask) {
        if (!(tu_info->tr_skip_mask & 0x10)) {
            const struct TBInfo *tb_info = &tu_info->tb_info[2];
            /* FIXME use sb_sig_map instead of last pos */
            if (ctu_dec->enable_lfnst && cu_flags & flg_pred_mode_flag && tb_info->sig_sb_map <= flg_pred_mode_flag&& !(cu_flags & flg_ibc_flag)) {
                int max_lfnst_pos = (log2_tb_h == log2_tb_w) && (log2_tb_w <= 3) ? 7 : 15;
                int nb_coeffs = check_lfnst_nb_coeffs(tb_info->last_pos);
                uint8_t is_mip = !!(cu_flags & flg_mip_flag);
                uint8_t allow_mip_lfnst = !is_mip || (log2_tb_h >= 4 && log2_tb_w >= 4);

                if (allow_mip_lfnst && nb_coeffs <= max_lfnst_pos && !!tb_info->last_pos) {
                    OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;
                    uint8_t is_dual = ctu_dec->transform_unit != &transform_unit_st;
                    uint8_t lfnst_flag = ovcabac_read_ae_lfnst_flag(cabac_ctx, is_dual);
                    tu_info->lfnst_flag = lfnst_flag;
                    if (lfnst_flag) {
                        uint8_t lfnst_idx = ovcabac_read_ae_lfnst_idx(cabac_ctx);
                        tu_info->lfnst_idx = lfnst_idx;
                    }
                }
            }

            if (!tu_info->lfnst_flag && !!tb_info->last_pos && (((cu_flags & flg_pred_mode_flag) && ctu_dec->mts_explicit_intra)
                                 || (!(cu_flags & flg_pred_mode_flag) && ctu_dec->mts_explicit_inter))
                && (log2_tb_w < 6) && (log2_tb_h < 6)
                && !(tb_info->sig_sb_map & (~0x000000000F0F0F0F))) {
                OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;

                uint8_t cu_mts_flag = ovcabac_read_ae_cu_mts_flag(cabac_ctx);
                tu_info->cu_mts_flag = cu_mts_flag;
                if (cu_mts_flag) {
                    uint8_t cu_mts_idx = ovcabac_read_ae_cu_mts_idx(cabac_ctx);
                    tu_info->cu_mts_idx = cu_mts_idx;
                }
            }
        }
    } else if (ctu_dec->transform_unit == &transform_unit_c && cbf_mask) {
        if (!(tu_info->tr_skip_mask))
        if (ctu_dec->enable_lfnst) {
            uint8_t jcbcr_flag = cbf_mask & 0x8;
            uint8_t cbf_mask_c = cbf_mask & 0x3;
            uint8_t can_lfnst = jcbcr_flag ? jcbcr_lfnst_check(tu_info, log2_tb_w, log2_tb_h)
                : chroma_lfnst_check(tu_info, cbf_mask_c, log2_tb_w, log2_tb_h);

            if (can_lfnst) {
                OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;
                uint8_t is_dual = ctu_dec->transform_unit != &transform_unit_st;
                uint8_t lfnst_flag = ovcabac_read_ae_lfnst_flag(cabac_ctx, is_dual);
                tu_info->lfnst_flag = lfnst_flag;
                if (lfnst_flag) {
                    uint8_t lfnst_idx = ovcabac_read_ae_lfnst_idx(cabac_ctx);
                    tu_info->lfnst_idx = lfnst_idx;
                }
            }
        }
    }
}

static int
transform_tree(OVCTUDec *const ctu_dec,
               unsigned int x0, unsigned int y0,
               unsigned int log2_tb_w, unsigned int log2_tb_h,
               unsigned int log2_max_tb_s, uint8_t rqt_root_cbf,
               CUFlags cu_flags, uint8_t tr_depth, struct TUInfo *tu_info)
{
    uint8_t split_v = log2_tb_w > log2_max_tb_s;
    uint8_t split_h = log2_tb_h > log2_max_tb_s;
    uint8_t nb_subtrees = tr_depth ? 1 : (1 << (split_v + split_h));
    if (log2_tb_w > 6 && log2_tb_h < 7) {

        transform_tree(ctu_dec, x0, y0, 6, log2_tb_h,
                       log2_max_tb_s, rqt_root_cbf, cu_flags,
                       tr_depth + 1, &tu_info[0]);

        tu_info[8].pos_offset = 1 << (log2_tb_h + 6);
        transform_tree(ctu_dec, x0 + 64, y0, 6, log2_tb_h,
                       log2_max_tb_s, rqt_root_cbf, cu_flags,
                       tr_depth + 1, &tu_info[8]);
        return 0;
    } else if (log2_tb_h > 6 && log2_tb_w < 7) {
        int residual_offset = tu_info[0].pos_offset;

        transform_tree(ctu_dec, x0, y0, log2_tb_w, 6,
                       log2_max_tb_s, rqt_root_cbf, cu_flags,
                       tr_depth + 1, &tu_info[0]);

        tu_info[8].pos_offset = 1 << (log2_tb_w + 6);
        transform_tree(ctu_dec, x0, y0 + 64, log2_tb_w, 6,
                       log2_max_tb_s, rqt_root_cbf, cu_flags,
                       tr_depth + 1, &tu_info[8]);
        return 0;
    }

    uint8_t cbf_mask = 0;

    if (split_v || split_h) {
        unsigned int tb_w1 = ((1 << log2_tb_w) >> split_v);
        unsigned int tb_h1 = ((1 << log2_tb_h) >> split_h);

        unsigned int log2_tb_w1 = log2_tb_w - split_v;
        unsigned int log2_tb_h1 = log2_tb_h - split_h;

        int residual_offset = tu_info[0].pos_offset;

        transform_tree(ctu_dec, x0, y0, log2_tb_w1, log2_tb_h1,
                       log2_max_tb_s, rqt_root_cbf, cu_flags,
                       tr_depth + 1, &tu_info[0 * nb_subtrees]);

        residual_offset += 1 << (log2_tb_w1 + log2_tb_h1);

        if (split_v) {
            tu_info[1 * nb_subtrees].pos_offset = residual_offset;
            transform_tree(ctu_dec, x0 + tb_w1, y0, log2_tb_w1, log2_tb_h1,
                           log2_max_tb_s, rqt_root_cbf, cu_flags,
                           tr_depth + 1, &tu_info[1 * nb_subtrees]);

            residual_offset += 1 << (log2_tb_w1 + log2_tb_h1);
        }

        if (split_h) {
            tu_info[2 * nb_subtrees].pos_offset = residual_offset;
            transform_tree(ctu_dec, x0, y0 + tb_h1, log2_tb_w1, log2_tb_h1,
                           log2_max_tb_s, rqt_root_cbf, cu_flags,
                           tr_depth + 1, &tu_info[2 * nb_subtrees]);

            residual_offset += 1 << (log2_tb_w1 + log2_tb_h1);
        }

        if (split_h && split_v) {
            tu_info[3 * nb_subtrees].pos_offset = residual_offset;
            transform_tree(ctu_dec, x0 + tb_w1, y0 + tb_h1, log2_tb_w1, log2_tb_h1,
                           log2_max_tb_s, rqt_root_cbf, cu_flags,
                           tr_depth + 1, &tu_info[3 * nb_subtrees]);
        }


    } else {

        cbf_mask = ctu_dec->transform_unit(ctu_dec, x0, y0, log2_tb_w, log2_tb_h, rqt_root_cbf, cu_flags, tr_depth, tu_info);
        tu_info->cbf_mask = cbf_mask;

    }

    return 0;
}

static int
sbt_half_ver(OVCTUDec *const ctu_dec,
             unsigned int x0, unsigned int y0,
             unsigned int log2_tb_w, unsigned int log2_tb_h,
             uint8_t sbt_pos, CUFlags cu_flags, struct TUInfo *tu_info)
{
    fill_ctb_bound(&ctu_dec->dbf_info, x0, y0, log2_tb_w, log2_tb_h);
    fill_ctb_bound_c(&ctu_dec->dbf_info, x0, y0, log2_tb_w, log2_tb_h);
    int qp_bd_offset = ctu_dec->qp_ctx.qp_bd_offset;
    if (!sbt_pos) {
        tu_info[0].is_sbt = 1;
        if (ctu_dec->mts_enabled && log2_tb_w - 1 <= 5 && log2_tb_h <= 5) {
            tu_info->cu_mts_flag = 1;
            tu_info->cu_mts_idx = 0x1;
        }
        uint8_t cbf_mask = ctu_dec->transform_unit(ctu_dec, x0, y0, log2_tb_w - 1, log2_tb_h,
                                                   1, cu_flags, 0, tu_info);

        struct DBFInfo *dbf_info = &ctu_dec->dbf_info;
        uint8_t qp_cb = ctu_dec->dequant_cb.qp - qp_bd_offset;
        uint8_t qp_cr = ctu_dec->dequant_cr.qp - qp_bd_offset;
        uint8_t qp_l  = ctu_dec->qp_ctx.current_qp;

        dbf_fill_qp_map(&dbf_info->qp_map_y, x0, y0, log2_tb_w, log2_tb_h, qp_l);
        dbf_fill_qp_map(&dbf_info->qp_map_cb, x0, y0, log2_tb_w, log2_tb_h, qp_cb);
        dbf_fill_qp_map(&dbf_info->qp_map_cr, x0, y0, log2_tb_w, log2_tb_h, qp_cr);
        ctu_dec->rcn_funcs.tmp.rcn_tu_st(ctu_dec, x0, y0, log2_tb_w - 1, log2_tb_h, cu_flags, cbf_mask, tu_info);

    } else {
        tu_info[0].is_sbt = 1;
        if (ctu_dec->mts_enabled && log2_tb_w - 1 <= 5 && log2_tb_h <= 5) {
            tu_info->cu_mts_flag = 1;
            tu_info->cu_mts_idx = 0x0;
        }
        uint8_t x1 = x0 + (1 << (log2_tb_w - 1));
        uint8_t cbf_mask = ctu_dec->transform_unit(ctu_dec, x0, y0, log2_tb_w - 1, log2_tb_h,
                                                   1, cu_flags, 0, tu_info);
        struct DBFInfo *dbf_info = &ctu_dec->dbf_info;
        uint8_t qp_cb = ctu_dec->dequant_cb.qp - qp_bd_offset;
        uint8_t qp_cr = ctu_dec->dequant_cr.qp - qp_bd_offset;
        uint8_t qp_l  = ctu_dec->qp_ctx.current_qp;

        dbf_fill_qp_map(&dbf_info->qp_map_y, x0, y0, log2_tb_w, log2_tb_h, qp_l);
        dbf_fill_qp_map(&dbf_info->qp_map_cb, x0, y0, log2_tb_w, log2_tb_h, qp_cb);
        dbf_fill_qp_map(&dbf_info->qp_map_cr, x0, y0, log2_tb_w, log2_tb_h, qp_cr);
        ctu_dec->rcn_funcs.tmp.rcn_tu_st(ctu_dec, x1, y0, log2_tb_w - 1, log2_tb_h, cu_flags, cbf_mask, tu_info);
    }

    return 0;
}

static int
sbt_half_hor(OVCTUDec *const ctu_dec,
             unsigned int x0, unsigned int y0,
             unsigned int log2_tb_w, unsigned int log2_tb_h,
             uint8_t sbt_pos, CUFlags cu_flags, struct TUInfo *tu_info)
{
    int qp_bd_offset = ctu_dec->qp_ctx.qp_bd_offset;
    fill_ctb_bound(&ctu_dec->dbf_info, x0, y0, log2_tb_w, log2_tb_h);
    fill_ctb_bound_c(&ctu_dec->dbf_info, x0, y0, log2_tb_w, log2_tb_h);
    if (!sbt_pos) {
        tu_info[0].is_sbt = 1;
        if (ctu_dec->mts_enabled && log2_tb_w <= 5 && log2_tb_h - 1 <= 5) {
            tu_info->cu_mts_flag = 1;
            tu_info->cu_mts_idx = 0x2;
        }

        uint8_t cbf_mask = ctu_dec->transform_unit(ctu_dec, x0, y0, log2_tb_w, log2_tb_h - 1,
                                                   1, cu_flags, 0, tu_info);

        struct DBFInfo *dbf_info = &ctu_dec->dbf_info;
        uint8_t qp_cb = ctu_dec->dequant_cb.qp - qp_bd_offset;
        uint8_t qp_cr = ctu_dec->dequant_cr.qp - qp_bd_offset;
        uint8_t qp_l  = ctu_dec->qp_ctx.current_qp;

        dbf_fill_qp_map(&dbf_info->qp_map_y, x0, y0, log2_tb_w, log2_tb_h, qp_l);
        dbf_fill_qp_map(&dbf_info->qp_map_cb, x0, y0, log2_tb_w, log2_tb_h, qp_cb);
        dbf_fill_qp_map(&dbf_info->qp_map_cr, x0, y0, log2_tb_w, log2_tb_h, qp_cr);
        ctu_dec->rcn_funcs.tmp.rcn_tu_st(ctu_dec, x0, y0, log2_tb_w, log2_tb_h - 1, cu_flags, cbf_mask, tu_info);

    } else {
        tu_info[0].is_sbt = 1;
        uint8_t y1 = y0 + (1 << (log2_tb_h - 1));
        uint8_t y3 = y0 + (1 << (log2_tb_h - 1));
        if (ctu_dec->mts_enabled && log2_tb_w <= 5 && log2_tb_h - 1 <= 5) {
            tu_info->cu_mts_flag = 1;
            tu_info->cu_mts_idx = 0x0;
        }
        uint8_t cbf_mask = ctu_dec->transform_unit(ctu_dec, x0, y0, log2_tb_w, log2_tb_h - 1,
                                                   1, cu_flags, 0, tu_info);


        struct DBFInfo *dbf_info = &ctu_dec->dbf_info;
        uint8_t qp_cb = ctu_dec->dequant_cb.qp - qp_bd_offset;
        uint8_t qp_cr = ctu_dec->dequant_cr.qp - qp_bd_offset;
        uint8_t qp_l  = ctu_dec->qp_ctx.current_qp;

        dbf_fill_qp_map(&dbf_info->qp_map_y, x0, y0, log2_tb_w, log2_tb_h, qp_l);
        dbf_fill_qp_map(&dbf_info->qp_map_cb, x0, y0, log2_tb_w, log2_tb_h, qp_cb);
        dbf_fill_qp_map(&dbf_info->qp_map_cr, x0, y0, log2_tb_w, log2_tb_h, qp_cr);
        ctu_dec->rcn_funcs.tmp.rcn_tu_st(ctu_dec, x0, y1, log2_tb_w, log2_tb_h - 1, cu_flags, cbf_mask, tu_info);
    }
    return 0;
}

static int
sbt_quad_ver(OVCTUDec *const ctu_dec,
             unsigned int x0, unsigned int y0,
             unsigned int log2_tb_w, unsigned int log2_tb_h,
             uint8_t sbt_pos, CUFlags cu_flags, struct TUInfo *tu_info)
{
    int qp_bd_offset = ctu_dec->qp_ctx.qp_bd_offset;
    fill_ctb_bound(&ctu_dec->dbf_info, x0, y0, log2_tb_w, log2_tb_h);
    fill_ctb_bound_c(&ctu_dec->dbf_info, x0, y0, log2_tb_w, log2_tb_h);
    if (!sbt_pos) {
        tu_info[0].is_sbt = 1;

        if (ctu_dec->mts_enabled && log2_tb_w - 2 <= 5 && log2_tb_h <= 5) {
            tu_info->cu_mts_flag = 1;
            tu_info->cu_mts_idx = 0x1;
        }

        uint8_t cbf_mask = ctu_dec->transform_unit(ctu_dec, x0, y0, log2_tb_w - 2, log2_tb_h,
                                                   1, cu_flags, 0, tu_info);


        struct DBFInfo *dbf_info = &ctu_dec->dbf_info;
        uint8_t qp_cb = ctu_dec->dequant_cb.qp - qp_bd_offset;
        uint8_t qp_cr = ctu_dec->dequant_cr.qp - qp_bd_offset;
        uint8_t qp_l  = ctu_dec->qp_ctx.current_qp;

        dbf_fill_qp_map(&dbf_info->qp_map_y, x0, y0, log2_tb_w, log2_tb_h, qp_l);
        dbf_fill_qp_map(&dbf_info->qp_map_cb, x0, y0, log2_tb_w, log2_tb_h, qp_cb);
        dbf_fill_qp_map(&dbf_info->qp_map_cr, x0, y0, log2_tb_w, log2_tb_h, qp_cr);
        ctu_dec->rcn_funcs.tmp.rcn_tu_st(ctu_dec, x0, y0, log2_tb_w - 2, log2_tb_h, cu_flags, cbf_mask, tu_info);

    } else {
        tu_info[0].is_sbt = 1;
        uint8_t x3 = x0 + (3 << (log2_tb_w - 2));

        if (ctu_dec->mts_enabled && log2_tb_w - 2 <= 5 && log2_tb_h <= 5) {
            tu_info->cu_mts_flag = 1;
            tu_info->cu_mts_idx = 0x0;
        }

        uint8_t cbf_mask = ctu_dec->transform_unit(ctu_dec, x0, y0, log2_tb_w - 2, log2_tb_h,
                                                   1, cu_flags, 0, tu_info);

        struct DBFInfo *dbf_info = &ctu_dec->dbf_info;
        uint8_t qp_cb = ctu_dec->dequant_cb.qp - qp_bd_offset;
        uint8_t qp_cr = ctu_dec->dequant_cr.qp - qp_bd_offset;
        uint8_t qp_l  = ctu_dec->qp_ctx.current_qp;

        dbf_fill_qp_map(&dbf_info->qp_map_y, x0, y0, log2_tb_w, log2_tb_h, qp_l);
        dbf_fill_qp_map(&dbf_info->qp_map_cb, x0, y0, log2_tb_w, log2_tb_h, qp_cb);
        dbf_fill_qp_map(&dbf_info->qp_map_cr, x0, y0, log2_tb_w, log2_tb_h, qp_cr);
        ctu_dec->rcn_funcs.tmp.rcn_tu_st(ctu_dec, x3, y0, log2_tb_w - 2, log2_tb_h, cu_flags, cbf_mask, tu_info);
    }
    return 0;
}

static int
sbt_quad_hor(OVCTUDec *const ctu_dec,
             unsigned int x0, unsigned int y0,
             unsigned int log2_tb_w, unsigned int log2_tb_h,
             uint8_t sbt_pos, CUFlags cu_flags, struct TUInfo *tu_info)
{
    int qp_bd_offset = ctu_dec->qp_ctx.qp_bd_offset;
    fill_ctb_bound(&ctu_dec->dbf_info, x0, y0, log2_tb_w, log2_tb_h);
    fill_ctb_bound_c(&ctu_dec->dbf_info, x0, y0, log2_tb_w, log2_tb_h);
    if (!sbt_pos) {
        tu_info[0].is_sbt = 1;

        if (ctu_dec->mts_enabled && log2_tb_w <= 5 && log2_tb_h - 2 <= 5) {
            tu_info->cu_mts_flag = 1;
            tu_info->cu_mts_idx = 0x2;
        }

        uint8_t cbf_mask = ctu_dec->transform_unit(ctu_dec, x0, y0, log2_tb_w, log2_tb_h - 2,
                                                   1, cu_flags, 0, tu_info);

        struct DBFInfo *dbf_info = &ctu_dec->dbf_info;
        uint8_t qp_cb = ctu_dec->dequant_cb.qp - qp_bd_offset;
        uint8_t qp_cr = ctu_dec->dequant_cr.qp - qp_bd_offset;
        uint8_t qp_l  = ctu_dec->qp_ctx.current_qp;

        dbf_fill_qp_map(&dbf_info->qp_map_y, x0, y0, log2_tb_w, log2_tb_h, qp_l);
        dbf_fill_qp_map(&dbf_info->qp_map_cb, x0, y0, log2_tb_w, log2_tb_h, qp_cb);
        dbf_fill_qp_map(&dbf_info->qp_map_cr, x0, y0, log2_tb_w, log2_tb_h, qp_cr);
        ctu_dec->rcn_funcs.tmp.rcn_tu_st(ctu_dec, x0, y0, log2_tb_w, log2_tb_h - 2, cu_flags, cbf_mask, tu_info);

    } else {
        tu_info[0].is_sbt = 1;
        uint8_t y3 = y0 + (3 << (log2_tb_h - 2));

        if (ctu_dec->mts_enabled && log2_tb_w <= 5 && log2_tb_h - 2 <= 5) {
            tu_info->cu_mts_flag = 1;
            tu_info->cu_mts_idx = 0x0;
        }

        uint8_t cbf_mask = ctu_dec->transform_unit(ctu_dec, x0, y0, log2_tb_w, log2_tb_h - 2,
                                                   1, cu_flags, 0, tu_info);


        struct DBFInfo *dbf_info = &ctu_dec->dbf_info;
        uint8_t qp_cb = ctu_dec->dequant_cb.qp - qp_bd_offset;
        uint8_t qp_cr = ctu_dec->dequant_cr.qp - qp_bd_offset;
        uint8_t qp_l  = ctu_dec->qp_ctx.current_qp;

        dbf_fill_qp_map(&dbf_info->qp_map_y, x0, y0, log2_tb_w, log2_tb_h, qp_l);
        dbf_fill_qp_map(&dbf_info->qp_map_cb, x0, y0, log2_tb_w, log2_tb_h, qp_cb);
        dbf_fill_qp_map(&dbf_info->qp_map_cr, x0, y0, log2_tb_w, log2_tb_h, qp_cr);
        ctu_dec->rcn_funcs.tmp.rcn_tu_st(ctu_dec, x0, y3, log2_tb_w, log2_tb_h - 2, cu_flags, cbf_mask, tu_info);

    }
    return 0;
}

static int
sbt_tree(OVCTUDec *const ctu_dec,
         unsigned int x0, unsigned int y0,
         unsigned int log2_tb_w, unsigned int log2_tb_h,
         uint8_t sbt_mode,
         uint8_t sbt_pos, CUFlags cu_flags, struct TUInfo *tu_info)
{
    /* FIXME remove this if zeroing is done by transform */
    switch (sbt_mode) {
    case 0x1:
         /*sbt_ver*/
         sbt_half_hor(ctu_dec, x0, y0, log2_tb_w, log2_tb_h,
                      sbt_pos, cu_flags, tu_info);
         break;
    case 0x2:
         /*sbt_hor*/
         sbt_half_ver(ctu_dec, x0, y0, log2_tb_w, log2_tb_h,
                      sbt_pos, cu_flags, tu_info);
         break;
    case 0x4:
         /*sbt_ver/quad*/
         sbt_quad_hor(ctu_dec, x0, y0, log2_tb_w, log2_tb_h,
                      sbt_pos, cu_flags, tu_info);
         break;
    case 0x8:
         /*sbt_hor/quad*/
         sbt_quad_ver(ctu_dec, x0, y0, log2_tb_w, log2_tb_h,
                      sbt_pos, cu_flags, tu_info);
         break;
    }
    struct DBFInfo *dbf_info = &ctu_dec->dbf_info;
    uint8_t qp_l  = ctu_dec->qp_ctx.current_qp;
    int qp_bd_offset = ctu_dec->qp_ctx.qp_bd_offset;
    dbf_fill_qp_map(&dbf_info->qp_map_y, x0, y0, log2_tb_w, log2_tb_h, qp_l);

    return 0;
}

static uint8_t
isp_subtree_v(OVCTUDec *const ctu_dec,
              unsigned int x0, unsigned int y0,
              unsigned int log2_cb_w, unsigned int log2_cb_h,
              CUFlags cu_flags,
              uint8_t intra_mode)
{
    OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;

    uint8_t cbf_flags = 0;
    int cbf = 0;
    int log2_pb_w = log2_cb_w - 2;
    int nb_pb;
    int i;
    struct ISPTUInfo tu_info = {0};
    struct TUInfo tu_info_c = {0};
    uint8_t cbf_mask_c = 0;

    /* height < 16 imposes restrictions on split dimension */
    if (log2_cb_h < 4 && (log2_pb_w <= (4 - log2_cb_h))) {
        log2_pb_w = 4 - log2_cb_h;
    }

    memset(ctu_dec->residual_y, 0, sizeof(int16_t) << (log2_cb_w + log2_cb_h));

    nb_pb = (1 << log2_cb_w) >> log2_pb_w;

    for (i = 0; i < nb_pb-1; ++i) {
        struct TBInfo *const tb_info = &tu_info.tb_info[i];
        cbf = ovcabac_read_ae_tu_cbf_luma_isp(cabac_ctx, cbf);
        cbf_flags <<= 1;
        cbf_flags |= cbf;

        if (cbf) {
            if (ctu_dec->delta_qp_enabled && ctu_dec->read_qp) {
                OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;
                int cu_qp_delta = ovcabac_read_ae_cu_delta_qp(cabac_ctx);
                derive_dequant_ctx(ctu_dec, &ctu_dec->qp_ctx, cu_qp_delta);
                ctu_dec->read_qp = 0;
            }
            uint16_t last_pos = ovcabac_read_ae_last_sig_pos(cabac_ctx, log2_pb_w, log2_cb_h);
            int16_t *coeffs_y = ctu_dec->residual_y + i * (1 << (log2_pb_w + log2_cb_h));

            tb_info->last_pos = last_pos;

            if (log2_pb_w < 2) {
                tb_info->sig_sb_map = ctu_dec->residual_coding_isp_v(ctu_dec, coeffs_y, log2_pb_w, log2_cb_h, last_pos);
            }else {
                tb_info->sig_sb_map = ctu_dec->residual_coding_l(ctu_dec, coeffs_y, log2_pb_w, log2_cb_h, last_pos);
            }
        }
    }

    if (ctu_dec->coding_tree != &dual_tree && ctu_dec->transform_unit == &transform_unit_st) {
        uint8_t intra_bdpcm_chroma_flag = 0;
        uint8_t tu_cbf_cb = ovcabac_read_ae_tu_cbf_cb(cabac_ctx, intra_bdpcm_chroma_flag);
        uint8_t tu_cbf_cr = ovcabac_read_ae_tu_cbf_cr(cabac_ctx, tu_cbf_cb, intra_bdpcm_chroma_flag);
        cbf_mask_c = (tu_cbf_cb << 1) | tu_cbf_cr;
    }

    cbf = !cbf_flags ? 1 : ovcabac_read_ae_tu_cbf_luma_isp(cabac_ctx, cbf);
    cbf_flags <<= 1;
    cbf_flags |= cbf;

    if (ctu_dec->delta_qp_enabled && cbf && ctu_dec->read_qp) {
        OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;
        int cu_qp_delta = ovcabac_read_ae_cu_delta_qp(cabac_ctx);
        derive_dequant_ctx(ctu_dec, &ctu_dec->qp_ctx, cu_qp_delta);
        ctu_dec->read_qp = 0;
    }

    #if 1
    if (ctu_dec->chroma_qp_offset_enabled && cbf_mask_c && ctu_dec->read_qp_c) {
        OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;
        int qp_bd_offset = ctu_dec->qp_ctx.qp_bd_offset;
        uint8_t length = ctu_dec->chroma_qp_offset_len;
        int cu_qp_delta = ovcabac_read_ae_cu_chroma_qp_offset(cabac_ctx, length);
        derive_dequant_ctx(ctu_dec, &ctu_dec->qp_ctx, 0);
        if (cu_qp_delta) {
            --cu_qp_delta;
            ctu_dec->qp_ctx.dqp_cb = ctu_dec->qp_ctx.cb_qp_offset_list[cu_qp_delta];
            ctu_dec->qp_ctx.dqp_cr = ctu_dec->qp_ctx.cr_qp_offset_list[cu_qp_delta];
            ctu_dec->qp_ctx.dqp_jcbcr = ctu_dec->qp_ctx.joint_cbcr_qp_offset_list[cu_qp_delta];
            ctu_dec->dequant_cb.qp = ov_clip(ctu_dec->dequant_cb.qp - qp_bd_offset + ctu_dec->qp_ctx.cb_qp_offset_list[cu_qp_delta], -qp_bd_offset, 63) + qp_bd_offset;
            ctu_dec->dequant_cr.qp = ov_clip(ctu_dec->dequant_cr.qp - qp_bd_offset + ctu_dec->qp_ctx.cr_qp_offset_list[cu_qp_delta], -qp_bd_offset, 63) + qp_bd_offset;
            ctu_dec->dequant_joint_cb_cr.qp = ov_clip(ctu_dec->dequant_joint_cb_cr.qp - qp_bd_offset + ctu_dec->qp_ctx.joint_cbcr_qp_offset_list[cu_qp_delta], -qp_bd_offset, 63) + qp_bd_offset;
            ctu_dec->dequant_cb_skip.qp   = OVMAX(ctu_dec->dequant_cb.qp, ctu_dec->qp_ctx.min_qp_prime_ts);
            ctu_dec->dequant_cr_skip.qp   = OVMAX(ctu_dec->dequant_cr.qp, ctu_dec->qp_ctx.min_qp_prime_ts);
            ctu_dec->dequant_jcbcr_skip.qp = OVMAX(ctu_dec->dequant_joint_cb_cr.qp, ctu_dec->qp_ctx.min_qp_prime_ts);
        }
        //derive_dequant_ctx(ctu_dec, &ctu_dec->qp_ctx, cu_qp_delta);
        //TODO update_chroma_qp
        ctu_dec->read_qp_c = 0;
    }
    #endif

    if (ctu_dec->jcbcr_enabled && cbf_mask_c) {
        uint8_t joint_cb_cr = ovcabac_read_ae_joint_cb_cr_flag(cabac_ctx,
                (cbf_mask_c & 0x3) - 1);
        cbf_mask_c |= joint_cb_cr << 3;
    }

    if (cbf) {
        struct TBInfo *const tb_info = &tu_info.tb_info[i];
        uint16_t last_pos = ovcabac_read_ae_last_sig_pos(cabac_ctx, log2_pb_w, log2_cb_h);
        int16_t *coeffs_y = ctu_dec->residual_y + i * (1 << (log2_pb_w + log2_cb_h));

        if (log2_pb_w <= 1) {
            tb_info->sig_sb_map = ctu_dec->residual_coding_isp_v(ctu_dec, coeffs_y, log2_pb_w, log2_cb_h, last_pos);
        } else {
            tb_info->sig_sb_map = ctu_dec->residual_coding_l(ctu_dec, coeffs_y, log2_pb_w, log2_cb_h, last_pos);
        }
        tb_info->last_pos = last_pos;
    }

    if (cbf_mask_c) {
        uint8_t jcbcr_flag = cbf_mask_c & 0x8;
        cbf_mask_c &= 0x3;
        if (jcbcr_flag) {

            residual_coding_jcbcr(ctu_dec, x0 >> 1, y0 >> 1, log2_cb_w - 1,
                                  log2_cb_h - 1, cbf_mask_c, 0, &tu_info_c);


        } else if (cbf_mask_c) {

            residual_coding_c(ctu_dec, x0 >> 1, y0 >> 1, log2_cb_w - 1,
                              log2_cb_h - 1, cbf_mask_c, 0, &tu_info_c);

        }
        cbf_mask_c |= jcbcr_flag;
    }

    if (ctu_dec->enable_lfnst && log2_pb_w > 1) {
        int max_lfnst_pos = (log2_cb_h == log2_pb_w) && (log2_pb_w <= 3) ? 7 : 15;
        uint8_t can_lfnst = (tu_info.tb_info[0].sig_sb_map | tu_info.tb_info[1].sig_sb_map | tu_info.tb_info[2].sig_sb_map | tu_info.tb_info[3].sig_sb_map) <= 1;

        can_lfnst &= check_lfnst_nb_coeffs(tu_info.tb_info[0].last_pos) <= max_lfnst_pos;
        can_lfnst &= check_lfnst_nb_coeffs(tu_info.tb_info[1].last_pos) <= max_lfnst_pos;
        can_lfnst &= check_lfnst_nb_coeffs(tu_info.tb_info[2].last_pos) <= max_lfnst_pos;
        can_lfnst &= check_lfnst_nb_coeffs(tu_info.tb_info[3].last_pos) <= max_lfnst_pos;

        if (cbf_mask_c) {
            uint8_t jcbcr_flag = cbf_mask_c & 0x8;
            int max_lfnst_pos_c = (log2_cb_h == log2_cb_w) && (log2_cb_w <= 4) ? 7 : 15;
            cbf_mask_c &= 0x3;
            if (log2_cb_w - 1 > 1 && log2_cb_h - 1 > 1) {
                if (jcbcr_flag) {
                    const struct TBInfo *const tb_info = &tu_info_c.tb_info[0];
                    uint8_t nb_coeffs = check_lfnst_nb_coeffs(tb_info->last_pos);

                    can_lfnst &= tb_info->sig_sb_map <= 0x1;
                    can_lfnst &= nb_coeffs <= max_lfnst_pos_c;
                } else {
                    uint8_t need_cb_chk = cbf_mask_c & 0x2;
                    uint8_t need_cr_chk = cbf_mask_c & 0x1;


                    if (need_cb_chk) {
                        const struct TBInfo *const tb_info_cb = &tu_info_c.tb_info[0];
                        uint8_t nb_coeffs_cb = check_lfnst_nb_coeffs(tb_info_cb->last_pos);
                        can_lfnst &= tb_info_cb->sig_sb_map <= 0x1;
                        can_lfnst &= nb_coeffs_cb <= max_lfnst_pos_c;
                    }

                    if (need_cr_chk) {
                        const struct TBInfo *const tb_info_cr = &tu_info_c.tb_info[1];
                        uint8_t nb_coeffs_cr = check_lfnst_nb_coeffs(tb_info_cr->last_pos);
                        can_lfnst &= tb_info_cr->sig_sb_map <= 0x1;
                        can_lfnst &= nb_coeffs_cr <= max_lfnst_pos_c;
                    }
                }
            }

            can_lfnst &= !tu_info_c.tr_skip_mask;
            cbf_mask_c |= jcbcr_flag;
        }

        if (can_lfnst) {
            uint8_t is_dual = ctu_dec->transform_unit != &transform_unit_st;
            uint8_t lfnst_flag = ovcabac_read_ae_lfnst_flag(cabac_ctx, is_dual);
            tu_info.lfnst_flag = lfnst_flag;
            if (lfnst_flag) {
                uint8_t lfnst_idx = ovcabac_read_ae_lfnst_idx(cabac_ctx);
                tu_info.lfnst_idx = lfnst_idx;
            }
        }
    }

    tu_info.cbf_mask = cbf_flags;

#if 1
    ctu_dec->rcn_funcs.tmp.recon_isp_subtree_v(ctu_dec, x0, y0, log2_cb_w, log2_cb_h, cu_flags, intra_mode, &tu_info);
#endif

    if (ctu_dec->transform_unit == &transform_unit_st) {

        int qp_bd_offset = ctu_dec->qp_ctx.qp_bd_offset;
        struct DBFInfo *dbf_info = &ctu_dec->dbf_info;
        uint8_t qp_cb = ctu_dec->dequant_cb.qp - qp_bd_offset;
        uint8_t qp_cr = ctu_dec->dequant_cr.qp - qp_bd_offset;

        fill_ctb_bound_c(&ctu_dec->dbf_info, x0, y0, log2_cb_w, log2_cb_h);

        uint8_t qp_l  = ctu_dec->qp_ctx.current_qp;
        dbf_fill_qp_map(&dbf_info->qp_map_y, x0, y0, log2_cb_w, log2_cb_h, qp_l);
        dbf_fill_qp_map(&dbf_info->qp_map_cb, x0, y0, log2_cb_w, log2_cb_h, qp_cb);
        dbf_fill_qp_map(&dbf_info->qp_map_cr, x0, y0, log2_cb_w, log2_cb_h, qp_cr);
        ctu_dec->rcn_funcs.tmp.rcn_tu_c(ctu_dec, x0 >> 1, y0 >> 1, log2_cb_w - 1, log2_cb_h - 1, 2, cbf_mask_c, &tu_info_c);

    } else {
        struct DBFInfo *dbf_info = &ctu_dec->dbf_info;
        uint8_t qp_l  = ctu_dec->qp_ctx.current_qp;
        dbf_fill_qp_map(&dbf_info->qp_map_y, x0, y0, log2_cb_w, log2_cb_h, qp_l);
    }

    return cbf_flags;
}

static uint8_t
isp_subtree_h(OVCTUDec *const ctu_dec,
              unsigned int x0, unsigned int y0,
              unsigned int log2_cb_w, unsigned int log2_cb_h,
              CUFlags cu_flags,
              uint8_t intra_mode)
{
    OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;
    int i;
    struct TUInfo tu_info_c = {0};

    uint8_t cbf = 0;
    uint8_t cbf_flags = 0;
    int log2_pb_h = log2_cb_h - 2;
    int nb_pb;
    int16_t *coeffs_y = ctu_dec->residual_y;
    struct ISPTUInfo tu_info = {0};
    uint8_t cbf_mask_c = 0;

    /* width < 16 imposes restrictions on split numbers */
    if (log2_cb_w < 4 && (log2_pb_h <= (4 - log2_cb_w))) {
        log2_pb_h = 4 - log2_cb_w;
    }

    int tb_s = 1 << (log2_cb_w + log2_pb_h);
    nb_pb = (1 << log2_cb_h) >> log2_pb_h;

    for (i = 0; i < nb_pb - 1; ++i) {
        struct TBInfo *const tb_info = &tu_info.tb_info[i];
        cbf = ovcabac_read_ae_tu_cbf_luma_isp(cabac_ctx, cbf);
        cbf_flags <<= 1;
        cbf_flags |= cbf;
        if (cbf) {
            if (ctu_dec->delta_qp_enabled && ctu_dec->read_qp) {
                OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;
                int cu_qp_delta = ovcabac_read_ae_cu_delta_qp(cabac_ctx);
                derive_dequant_ctx(ctu_dec, &ctu_dec->qp_ctx, cu_qp_delta);
                ctu_dec->read_qp = 0;
            }
            uint16_t last_pos = ovcabac_read_ae_last_sig_pos(cabac_ctx, log2_cb_w, log2_pb_h);

            if (log2_pb_h <= 1) {
                tb_info->sig_sb_map = ctu_dec->residual_coding_isp_h(ctu_dec, coeffs_y, log2_cb_w, log2_pb_h, last_pos);
            } else {
                tb_info->sig_sb_map = ctu_dec->residual_coding_l(ctu_dec, coeffs_y, log2_cb_w, log2_pb_h, last_pos);
            }
            tb_info->last_pos = last_pos;
        }
        coeffs_y += tb_s;
    }

    if (ctu_dec->coding_tree != &dual_tree && ctu_dec->transform_unit == &transform_unit_st) {
        uint8_t intra_bdpcm_chroma_flag = 0;
        uint8_t tu_cbf_cb = ovcabac_read_ae_tu_cbf_cb(cabac_ctx, intra_bdpcm_chroma_flag);
        uint8_t tu_cbf_cr = ovcabac_read_ae_tu_cbf_cr(cabac_ctx, tu_cbf_cb, intra_bdpcm_chroma_flag);
        cbf_mask_c = (tu_cbf_cb << 1) | tu_cbf_cr;
    }

    cbf = !cbf_flags ? 1 : ovcabac_read_ae_tu_cbf_luma_isp(cabac_ctx, cbf);
    cbf_flags <<= 1;
    cbf_flags |= cbf;

    if (ctu_dec->delta_qp_enabled && cbf && ctu_dec->read_qp) {
        OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;
        int cu_qp_delta = ovcabac_read_ae_cu_delta_qp(cabac_ctx);
        derive_dequant_ctx(ctu_dec, &ctu_dec->qp_ctx, cu_qp_delta);
        ctu_dec->read_qp = 0;
    }

    #if 1
    if (ctu_dec->chroma_qp_offset_enabled && cbf_mask_c && ctu_dec->read_qp_c) {
        OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;
        uint8_t length = ctu_dec->chroma_qp_offset_len;
        int cu_qp_delta = ovcabac_read_ae_cu_chroma_qp_offset(cabac_ctx, length);
        derive_dequant_ctx(ctu_dec, &ctu_dec->qp_ctx, 0);
        if (cu_qp_delta) {
            --cu_qp_delta;
            ctu_dec->qp_ctx.dqp_cb = ctu_dec->qp_ctx.cb_qp_offset_list[cu_qp_delta];
            ctu_dec->qp_ctx.dqp_cr = ctu_dec->qp_ctx.cr_qp_offset_list[cu_qp_delta];
            ctu_dec->qp_ctx.dqp_jcbcr = ctu_dec->qp_ctx.joint_cbcr_qp_offset_list[cu_qp_delta];
            int qp_bd_offset = ctu_dec->qp_ctx.qp_bd_offset;
            ctu_dec->dequant_cb.qp = ov_clip(ctu_dec->dequant_cb.qp - qp_bd_offset + ctu_dec->qp_ctx.cb_qp_offset_list[cu_qp_delta], -qp_bd_offset, 63) + qp_bd_offset;
            ctu_dec->dequant_cr.qp = ov_clip(ctu_dec->dequant_cr.qp - qp_bd_offset + ctu_dec->qp_ctx.cr_qp_offset_list[cu_qp_delta], -qp_bd_offset, 63) + qp_bd_offset;
            ctu_dec->dequant_joint_cb_cr.qp = ov_clip(ctu_dec->dequant_joint_cb_cr.qp - qp_bd_offset + ctu_dec->qp_ctx.joint_cbcr_qp_offset_list[cu_qp_delta], -qp_bd_offset, 63) + qp_bd_offset;
            ctu_dec->dequant_cb_skip.qp   = OVMAX(ctu_dec->dequant_cb.qp, ctu_dec->qp_ctx.min_qp_prime_ts);
            ctu_dec->dequant_cr_skip.qp   = OVMAX(ctu_dec->dequant_cr.qp, ctu_dec->qp_ctx.min_qp_prime_ts);
            ctu_dec->dequant_jcbcr_skip.qp = OVMAX(ctu_dec->dequant_joint_cb_cr.qp, ctu_dec->qp_ctx.min_qp_prime_ts);
            ctu_dec->read_qp = 0;
        }
        ctu_dec->read_qp_c = 0;
    }
    #endif

    if (ctu_dec->jcbcr_enabled && cbf_mask_c) {
        uint8_t joint_cb_cr = ovcabac_read_ae_joint_cb_cr_flag(cabac_ctx,
                (cbf_mask_c & 0x3) - 1);
        cbf_mask_c |= joint_cb_cr << 3;
    }

    if (cbf) {
        struct TBInfo *const tb_info = &tu_info.tb_info[i];
        uint16_t last_pos = ovcabac_read_ae_last_sig_pos(cabac_ctx, log2_cb_w, log2_pb_h);

        if (log2_pb_h <= 1) {
            tb_info->sig_sb_map = ctu_dec->residual_coding_isp_h(ctu_dec, coeffs_y, log2_cb_w, log2_pb_h, last_pos);
        } else {
            tb_info->sig_sb_map = ctu_dec->residual_coding_l(ctu_dec, coeffs_y, log2_cb_w, log2_pb_h, last_pos);
        }
        tb_info->last_pos = last_pos;
    }

    if (cbf_mask_c) {
        uint8_t jcbcr_flag = cbf_mask_c & 0x8;
        cbf_mask_c &= 0x3;
        if (jcbcr_flag) {

            residual_coding_jcbcr(ctu_dec, x0 >> 1, y0 >> 1, log2_cb_w - 1,
                                  log2_cb_h - 1, cbf_mask_c, 0, &tu_info_c);


        } else if (cbf_mask_c) {

            residual_coding_c(ctu_dec, x0 >> 1, y0 >> 1, log2_cb_w - 1,
                              log2_cb_h - 1, cbf_mask_c,  0, &tu_info_c);

        }
        cbf_mask_c |= jcbcr_flag;

    }

    if (ctu_dec->enable_lfnst && log2_pb_h > 1) {
        uint8_t can_lfnst = (tu_info.tb_info[0].sig_sb_map | tu_info.tb_info[1].sig_sb_map | tu_info.tb_info[2].sig_sb_map | tu_info.tb_info[3].sig_sb_map) <= 1;
        int max_lfnst_pos = (log2_pb_h == log2_cb_w) && (log2_cb_w <= 3) ? 7 : 15;

        can_lfnst &= check_lfnst_nb_coeffs(tu_info.tb_info[0].last_pos) <= max_lfnst_pos;
        can_lfnst &= check_lfnst_nb_coeffs(tu_info.tb_info[1].last_pos) <= max_lfnst_pos;
        can_lfnst &= check_lfnst_nb_coeffs(tu_info.tb_info[2].last_pos) <= max_lfnst_pos;
        can_lfnst &= check_lfnst_nb_coeffs(tu_info.tb_info[3].last_pos) <= max_lfnst_pos;

        if (cbf_mask_c) {
            int max_lfnst_pos_c = (log2_cb_h == log2_cb_w) && (log2_cb_w <= 4) ? 7 : 15;
            uint8_t jcbcr_flag = cbf_mask_c & 0x8;
            cbf_mask_c &= 0x3;
            if (log2_cb_w - 1 > 1 && log2_cb_h - 1 > 1) {
                if (jcbcr_flag) {
                    const struct TBInfo *const tb_info = &tu_info_c.tb_info[0];
                    uint8_t nb_coeffs = check_lfnst_nb_coeffs(tb_info->last_pos);
                    can_lfnst &= tb_info->sig_sb_map <= 0x1;
                    can_lfnst &= nb_coeffs <= max_lfnst_pos_c;
                } else {
                    uint8_t need_cb_chk = cbf_mask_c & 0x2;
                    uint8_t need_cr_chk = cbf_mask_c & 0x1;

                    if (need_cb_chk) {
                        const struct TBInfo *const tb_info_cb = &tu_info_c.tb_info[0];
                        uint8_t nb_coeffs_cb = check_lfnst_nb_coeffs(tb_info_cb->last_pos);
                        can_lfnst &= tb_info_cb->sig_sb_map <= 0x1;
                        can_lfnst &= nb_coeffs_cb <= max_lfnst_pos_c;
                    }

                    if (need_cr_chk) {
                        const struct TBInfo *const tb_info_cr = &tu_info_c.tb_info[1];
                        uint8_t nb_coeffs_cr = check_lfnst_nb_coeffs(tb_info_cr->last_pos);
                        can_lfnst &= tb_info_cr->sig_sb_map <= 0x1;
                        can_lfnst &= nb_coeffs_cr <= max_lfnst_pos_c;
                    }
                }
            }
            can_lfnst &= !tu_info_c.tr_skip_mask;
            cbf_mask_c |= jcbcr_flag;
        }

        if (can_lfnst) {
            uint8_t is_dual = ctu_dec->transform_unit != &transform_unit_st;
            uint8_t lfnst_flag = ovcabac_read_ae_lfnst_flag(cabac_ctx, is_dual);
            tu_info.lfnst_flag = lfnst_flag;
            if (lfnst_flag) {
                uint8_t lfnst_idx = ovcabac_read_ae_lfnst_idx(cabac_ctx);
                tu_info.lfnst_idx = lfnst_idx;
            }
        }
    }

    tu_info.cbf_mask = cbf_flags;

#if 1
    ctu_dec->rcn_funcs.tmp.recon_isp_subtree_h(ctu_dec, x0, y0, log2_cb_w, log2_cb_h, cu_flags, intra_mode, &tu_info);
#endif
    if (ctu_dec->transform_unit == &transform_unit_st) {
        int qp_bd_offset = ctu_dec->qp_ctx.qp_bd_offset;
        struct DBFInfo *dbf_info = &ctu_dec->dbf_info;
        uint8_t qp_cb = ctu_dec->dequant_cb.qp - qp_bd_offset;
        uint8_t qp_cr = ctu_dec->dequant_cr.qp - qp_bd_offset;

        fill_ctb_bound_c(&ctu_dec->dbf_info, x0, y0, log2_cb_w, log2_cb_h);

        uint8_t qp_l  = ctu_dec->qp_ctx.current_qp;
        dbf_fill_qp_map(&dbf_info->qp_map_y, x0, y0, log2_cb_w, log2_cb_h, qp_l);
        dbf_fill_qp_map(&dbf_info->qp_map_cb, x0, y0, log2_cb_w, log2_cb_h, qp_cb);
        dbf_fill_qp_map(&dbf_info->qp_map_cr, x0, y0, log2_cb_w, log2_cb_h, qp_cr);
        ctu_dec->rcn_funcs.tmp.rcn_tu_c(ctu_dec, x0 >> 1, y0 >> 1, log2_cb_w - 1, log2_cb_h - 1, 2, cbf_mask_c, &tu_info_c);

    } else {
        struct DBFInfo *dbf_info = &ctu_dec->dbf_info;
        uint8_t qp_l  = ctu_dec->qp_ctx.current_qp;
        dbf_fill_qp_map(&dbf_info->qp_map_y, x0, y0, log2_cb_w, log2_cb_h, qp_l);
    }

    return cbf_flags;
}

static uint8_t
sbt_allowed(uint8_t log2_tb_w, uint8_t log2_tb_h)
{
  uint8_t allow_half_v = log2_tb_w >= 3;
  uint8_t allow_half_h = log2_tb_h >= 3;
  uint8_t allow_quad_v = log2_tb_w >= 4;
  uint8_t allow_quad_h = log2_tb_h >= 4;

  uint8_t sbt_status = 0;

  sbt_status |= allow_half_h;
  sbt_status |= allow_half_v << 1;
  sbt_status |= allow_quad_h << 2;
  sbt_status |= allow_quad_v << 3;

  return sbt_status;
}

static uint8_t
sbt_mode(OVCABACCtx *const cabac_ctx, uint8_t log2_tb_w, uint8_t log2_tb_h, uint8_t sbt_mask)
{
    uint8_t sbt_mode;
    uint8_t sbt_quad_flag;
    uint8_t sbt_hor_flag;
    uint8_t sbt_pos_flag;

    if ((sbt_mask & 0xC) && (sbt_mask & 0x3)) {
        sbt_quad_flag = ovcabac_read_ae_sbt_quad_flag(cabac_ctx);
    } else {
        sbt_quad_flag = 0;
    }

    if ((sbt_quad_flag && (sbt_mask & 0x4) && (sbt_mask & 0x8)) ||
        (!sbt_quad_flag && (sbt_mask & 0x1) && (sbt_mask & 0x2))) {
        sbt_hor_flag = ovcabac_read_ae_sbt_hor_flag(cabac_ctx, log2_tb_w, log2_tb_h);
    } else {
        sbt_hor_flag = (sbt_quad_flag && (sbt_mask & 0x4)) || (!sbt_quad_flag && (sbt_mask & 0x1));
    }

    sbt_pos_flag = ovcabac_read_ae_sbt_pos_flag(cabac_ctx);

    sbt_mode = ((1 << (!sbt_hor_flag)) << (sbt_quad_flag << 1));
    sbt_mode |= sbt_pos_flag << 7;

    return sbt_mode;
}

int
transform_unit_wrap(OVCTUDec *const ctu_dec,
                    const OVPartInfo *const part_ctx,
                    uint8_t x0, uint8_t y0,
                    uint8_t log2_cb_w, uint8_t log2_cb_h,
                    VVCCU cu)
{
    uint8_t split_tu = ((log2_cb_w > part_ctx->log2_max_tb_s) | (log2_cb_h > part_ctx->log2_max_tb_s));
    if (cu.cu_flags & flg_pred_mode_flag && !(cu.cu_flags & flg_ibc_flag)) {
        /* INTRA */
        if (!(cu.cu_flags & flg_isp_flag)) {
            /*FIXME check if part_ctx mandatory for transform_tree */
            struct TUInfo tu_info[16] = {0};

            transform_tree(ctu_dec, x0, y0, log2_cb_w, log2_cb_h,
                           part_ctx->log2_max_tb_s, 0, cu.cu_flags, 0, tu_info);

            if (!split_tu) {
                lfnst_mts(ctu_dec, log2_cb_w, log2_cb_h, cu.cu_flags, tu_info);
            }

            if (ctu_dec->coding_unit == coding_unit_intra_c) {
                int qp_bd_offset = ctu_dec->qp_ctx.qp_bd_offset;
                struct DBFInfo *dbf_info = &ctu_dec->dbf_info;
                uint8_t qp_cb = ctu_dec->dequant_cb.qp - qp_bd_offset;
                uint8_t qp_cr = ctu_dec->dequant_cr.qp - qp_bd_offset;

                dbf_fill_qp_map(&dbf_info->qp_map_cb, x0 << 1, y0 << 1, log2_cb_w + 1, log2_cb_h + 1, qp_cb);
                dbf_fill_qp_map(&dbf_info->qp_map_cr, x0 << 1, y0 << 1, log2_cb_w + 1, log2_cb_h + 1, qp_cr);
            } else if (ctu_dec->coding_unit != &coding_unit_intra) {
                int qp_bd_offset = ctu_dec->qp_ctx.qp_bd_offset;
                struct DBFInfo *dbf_info = &ctu_dec->dbf_info;
                uint8_t qp_l  = ctu_dec->qp_ctx.current_qp;
                uint8_t qp_cb = ctu_dec->dequant_cb.qp - qp_bd_offset;
                uint8_t qp_cr = ctu_dec->dequant_cr.qp - qp_bd_offset;

                dbf_fill_qp_map(&dbf_info->qp_map_y, x0, y0, log2_cb_w, log2_cb_h, qp_l);
                dbf_fill_qp_map(&dbf_info->qp_map_cb, x0, y0, log2_cb_w, log2_cb_h, qp_cb);
                dbf_fill_qp_map(&dbf_info->qp_map_cr, x0, y0, log2_cb_w, log2_cb_h, qp_cr);
            } else {
                struct DBFInfo *dbf_info = &ctu_dec->dbf_info;
                uint8_t qp_l  = ctu_dec->qp_ctx.current_qp;
                dbf_fill_qp_map(&dbf_info->qp_map_y, x0, y0, log2_cb_w, log2_cb_h, qp_l);
            }

            ctu_dec->rcn_funcs.tmp.rcn_transform_tree(ctu_dec, x0, y0, log2_cb_w, log2_cb_h, part_ctx->log2_max_tb_s,
                                                      0, cu.cu_flags, tu_info);
        } else {
            uint8_t isp_mode = cu.cu_opaque;
            uint8_t intra_mode = cu.cu_mode_idx;

            if (isp_mode) {
                /* Disable CCLM in 64x64 ISP CU*/
                ctu_dec->tmp_disable_cclm |= log2_cb_h == log2_cb_w && log2_cb_h == 6;

                uint8_t x0_unit = x0 >> LOG2_MIN_CU_S;
                uint8_t y0_unit = y0 >> LOG2_MIN_CU_S;
                uint8_t nb_unit_w = (1 << log2_cb_w) >> LOG2_MIN_CU_S;
                uint8_t nb_unit_h = (1 << log2_cb_h) >> LOG2_MIN_CU_S;

                ctu_field_set_rect_bitfield(&ctu_dec->rcn_ctx.progress_field,
                                            x0_unit, y0_unit,
                                            nb_unit_w, nb_unit_h);

                if (isp_mode == 2) {
                    isp_subtree_v(ctu_dec, x0, y0, log2_cb_w, log2_cb_h, cu.cu_flags, intra_mode);
                } else if (isp_mode == 1) {
                    isp_subtree_h(ctu_dec, x0, y0, log2_cb_w, log2_cb_h, cu.cu_flags, intra_mode);
                }
            }
        }
    } else {
        OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;
        uint8_t cu_merge_flag = !!(cu.cu_flags & flg_merge_flag);
        uint8_t cu_skip_flag  = !!(cu.cu_flags & flg_cu_skip_flag);

        struct TUInfo tu_info[16] = {0};

        uint8_t rqt_root_cbf = !cu_skip_flag && (cu_merge_flag || ovcabac_read_ae_root_cbf(cabac_ctx));

        if (rqt_root_cbf) {
            uint8_t sbt_flag = 0;
            if (!split_tu && !(cu.cu_flags & flg_ibc_flag)) {
                if (ctu_dec->sbt_enabled && !ctu_dec->tmp_ciip && !(cu.cu_flags & flg_ibc_flag)) {
                    uint8_t sbt_mask = sbt_allowed(log2_cb_w, log2_cb_h);
                    if (sbt_mask) {
                        sbt_flag = ovcabac_read_ae_sbt_flag(cabac_ctx, log2_cb_w, log2_cb_h);
                        if (sbt_flag) {
                            uint8_t sbt_type = sbt_mode(cabac_ctx, log2_cb_w, log2_cb_h, sbt_mask);
                            uint8_t sbt_pos = !!(sbt_type & 0x80);

                            sbt_type &= 0x7F;

                            sbt_tree(ctu_dec, x0, y0, log2_cb_w, log2_cb_h,
                                     sbt_type, !!sbt_pos, cu.cu_flags, tu_info);
                        }
                    }
                }

                if (!sbt_flag) {
                    uint8_t cbf_mask = transform_unit_st(ctu_dec, x0, y0, log2_cb_w, log2_cb_h,
                                                         1, cu.cu_flags, 0, tu_info);

                    tu_info->cbf_mask = cbf_mask;

                    lfnst_mts(ctu_dec, log2_cb_w, log2_cb_h, cu.cu_flags, tu_info);
                }

            } else {
                transform_tree(ctu_dec, x0, y0, log2_cb_w, log2_cb_h,
                               part_ctx->log2_max_tb_s, 1, cu.cu_flags, 0, tu_info);
            }

            if (!sbt_flag) {

                int qp_bd_offset = ctu_dec->qp_ctx.qp_bd_offset;
                struct DBFInfo *dbf_info = &ctu_dec->dbf_info;
                uint8_t qp_l  = ctu_dec->qp_ctx.current_qp;
                uint8_t qp_cb = ctu_dec->dequant_cb.qp - qp_bd_offset;
                uint8_t qp_cr = ctu_dec->dequant_cr.qp - qp_bd_offset;

                dbf_fill_qp_map(&dbf_info->qp_map_y, x0, y0, log2_cb_w, log2_cb_h, qp_l);
                dbf_fill_qp_map(&dbf_info->qp_map_cb, x0, y0, log2_cb_w, log2_cb_h, qp_cb);
                dbf_fill_qp_map(&dbf_info->qp_map_cr, x0, y0, log2_cb_w, log2_cb_h, qp_cr);

                ctu_dec->rcn_funcs.tmp.rcn_transform_tree(ctu_dec, x0, y0, log2_cb_w, log2_cb_h, part_ctx->log2_max_tb_s,
                                                          0, cu.cu_flags, tu_info);
            }
        } else {
                int qp_bd_offset = ctu_dec->qp_ctx.qp_bd_offset;
                struct DBFInfo *dbf_info = &ctu_dec->dbf_info;
                uint8_t qp_l  = ctu_dec->qp_ctx.current_qp;
                uint8_t qp_cb = ctu_dec->dequant_cb.qp - qp_bd_offset;
                uint8_t qp_cr = ctu_dec->dequant_cr.qp - qp_bd_offset;

                dbf_fill_qp_map(&dbf_info->qp_map_y, x0, y0, log2_cb_w, log2_cb_h, qp_l);
                dbf_fill_qp_map(&dbf_info->qp_map_cb, x0, y0, log2_cb_w, log2_cb_h, qp_cb);
                dbf_fill_qp_map(&dbf_info->qp_map_cr, x0, y0, log2_cb_w, log2_cb_h, qp_cr);

                ctu_dec->rcn_funcs.tmp.rcn_transform_tree(ctu_dec, x0, y0, log2_cb_w, log2_cb_h, part_ctx->log2_max_tb_s,
                                                          0, cu.cu_flags, tu_info);
        }
    }
    return 0;
}

