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
    int min_qp_prime_ts = qp_ctx->min_qp_prime_ts;
    int base_qp = ((qp_ctx->current_qp + cu_qp_delta + 64 + 2 * qp_bd_offset) % (64 + qp_bd_offset)) - qp_bd_offset;
    int base_qp_c = ov_clip(base_qp, -qp_bd_offset, 63);
    ctudec->dequant_luma.qp = base_qp + qp_bd_offset;

    ctudec->dequant_cb.qp = ov_clip(qp_ctx->chroma_qp_map_cb[base_qp_c] + qp_ctx->cb_offset + qp_ctx->dqp_cb, -qp_bd_offset, 63) + qp_bd_offset;
    ctudec->dequant_cr.qp = ov_clip(qp_ctx->chroma_qp_map_cr[base_qp_c] + qp_ctx->cr_offset + qp_ctx->dqp_cr, -qp_bd_offset, 63) + qp_bd_offset;
    ctudec->dequant_joint_cb_cr.qp = ov_clip(qp_ctx->chroma_qp_map_jcbcr[base_qp_c] + qp_ctx->jcbcr_offset + qp_ctx->dqp_jcbcr, -qp_bd_offset, 63) + qp_bd_offset;
    ctudec->dequant_luma_skip.qp = OVMAX(ctudec->dequant_luma.qp, min_qp_prime_ts);
    ctudec->dequant_cb_skip.qp   = OVMAX(ctudec->dequant_cb.qp, min_qp_prime_ts);
    ctudec->dequant_cr_skip.qp   = OVMAX(ctudec->dequant_cr.qp, min_qp_prime_ts);
    ctudec->dequant_jcbcr_skip.qp = OVMAX(ctudec->dequant_joint_cb_cr.qp, min_qp_prime_ts);

    ctudec->qp_ctx.current_qp = base_qp;
}
#endif

static void
dequant_sb(int16_t *const sb_coeffs, int scale, int shift)
{
    int add = (1 << shift) >> 1;
    for( int i = 0; i < 16 ; i++ ){
        sb_coeffs[i] = ov_clip_intp2((int32_t)(sb_coeffs[i] * scale + add) >> shift,
                MAX_LOG2_TR_RANGE + 1);
    }
}

#define DQ_SHIFT(log2_tb_s, sh_dpq_flag) BITDEPTH - 5 + sh_dpq_flag + (log2_tb_s >> 1) + (log2_tb_s & 1)
static const uint8_t dequant_shift_lut_dpq[13] =
{
    DQ_SHIFT( 0, 1),
    DQ_SHIFT( 1, 1),
    DQ_SHIFT( 2, 1),
    DQ_SHIFT( 3, 1),
    DQ_SHIFT( 4, 1),
    DQ_SHIFT( 5, 1),
    DQ_SHIFT( 6, 1),
    DQ_SHIFT( 7, 1),
    DQ_SHIFT( 8, 1),
    DQ_SHIFT( 9, 1),
    DQ_SHIFT(10, 1),
    DQ_SHIFT(11, 1),
    DQ_SHIFT(12, 1)
};


static const uint8_t dequant_shift_lut_sdh[13] =
{
    DQ_SHIFT( 0, 0),
    DQ_SHIFT( 1, 0),
    DQ_SHIFT( 2, 0),
    DQ_SHIFT( 3, 0),
    DQ_SHIFT( 4, 0),
    DQ_SHIFT( 5, 0),
    DQ_SHIFT( 6, 0),
    DQ_SHIFT( 7, 0),
    DQ_SHIFT( 8, 0),
    DQ_SHIFT( 9, 0),
    DQ_SHIFT(10, 0),
    DQ_SHIFT(11, 0),
    DQ_SHIFT(12, 0)
};

static struct IQScale
derive_dequant_sdh(int qp, uint8_t log2_tb_w, uint8_t log2_tb_h)
{
    const uint8_t log2_tb_s = log2_tb_w + log2_tb_h;
    struct IQScale dequant_params;

    int scale = inverse_quant_scale_lut[log2_tb_s & 1][qp % 6];

    dequant_params.shift = dequant_shift_lut_sdh[log2_tb_s];
    dequant_params.scale = 16 * scale << (qp / 6);
    dequant_params.dequant_sb = &dequant_sb;

    return dequant_params;
}

static struct IQScale
derive_dequant_dpq(int qp, uint8_t log2_tb_w, uint8_t log2_tb_h)
{
    struct IQScale dequant_params;

    const uint8_t log2_tb_s = log2_tb_w + log2_tb_h;

    uint32_t scale = inverse_quant_scale_lut[log2_tb_s & 1][(qp + 1) % 6];

    dequant_params.shift = dequant_shift_lut_dpq[log2_tb_s];
    dequant_params.scale = 16 * scale << ((qp + 1) / 6);
    dequant_params.dequant_sb = &dequant_sb;

    return dequant_params;
}

static struct IQScale
derive_dequant_ts(int qp, uint8_t log2_tb_w, uint8_t log2_tb_h)
{
    struct IQScale dequant_params;

    int shift = IQUANT_SHIFT;
    int scale = inverse_quant_scale_lut[0][qp % 6];

    dequant_params.shift = shift;
    dequant_params.scale = scale << (qp / 6);
    dequant_params.dequant_sb = &dequant_sb;

    return dequant_params;
}

static void
dequant_tb_4x4(int16_t *dst, const int16_t *src, int scale, int shift,
               uint8_t log2_tb_w, uint8_t log2_tb_h, uint64_t sig_sb_map)
{
    int add = (1 << shift) >> 1;
    int nb_rows = 1 << OVMIN(5, log2_tb_h);//derive_nb_cols(sig_sb_map);
    int nb_cols = 1 << OVMIN(5, log2_tb_w);//derive_nb_rows(sig_sb_map);
    uint8_t src_stride = 1 << (OVMIN(5, log2_tb_w));
    uint8_t dst_stride = 1 << (OVMIN(5, log2_tb_w));

    /* Force sig_sb_map to one in case of DC coefficient */
    sig_sb_map |= !sig_sb_map;

    for (int i = 0; i < nb_rows/4 ; i++) {
        uint8_t sig_sb_row = sig_sb_map >> (i << 3);
        for (int j = 0; j < nb_cols/4 ; j++) {
            int16_t *_dst = dst + (j << 2);
            const int16_t *_src = src + (j << 4);
            if (sig_sb_row & 0x1) {
                _dst[0] = ov_clip_intp2((int32_t)(_src[ 0] * scale + add) >> shift, MAX_LOG2_TR_RANGE + 1);
                _dst[1] = ov_clip_intp2((int32_t)(_src[ 1] * scale + add) >> shift, MAX_LOG2_TR_RANGE + 1);
                _dst[2] = ov_clip_intp2((int32_t)(_src[ 2] * scale + add) >> shift, MAX_LOG2_TR_RANGE + 1);
                _dst[3] = ov_clip_intp2((int32_t)(_src[ 3] * scale + add) >> shift, MAX_LOG2_TR_RANGE + 1);

                _dst += dst_stride;

                _dst[0] = ov_clip_intp2((int32_t)(_src[ 4] * scale + add) >> shift, MAX_LOG2_TR_RANGE + 1);
                _dst[1] = ov_clip_intp2((int32_t)(_src[ 5] * scale + add) >> shift, MAX_LOG2_TR_RANGE + 1);
                _dst[2] = ov_clip_intp2((int32_t)(_src[ 6] * scale + add) >> shift, MAX_LOG2_TR_RANGE + 1);
                _dst[3] = ov_clip_intp2((int32_t)(_src[ 7] * scale + add) >> shift, MAX_LOG2_TR_RANGE + 1);

                _dst += dst_stride;

                _dst[0] = ov_clip_intp2((int32_t)(_src[ 8] * scale + add) >> shift, MAX_LOG2_TR_RANGE + 1);
                _dst[1] = ov_clip_intp2((int32_t)(_src[ 9] * scale + add) >> shift, MAX_LOG2_TR_RANGE + 1);
                _dst[2] = ov_clip_intp2((int32_t)(_src[10] * scale + add) >> shift, MAX_LOG2_TR_RANGE + 1);
                _dst[3] = ov_clip_intp2((int32_t)(_src[11] * scale + add) >> shift, MAX_LOG2_TR_RANGE + 1);

                _dst += dst_stride;

                _dst[0] = ov_clip_intp2((int32_t)(_src[12] * scale + add) >> shift, MAX_LOG2_TR_RANGE + 1);
                _dst[1] = ov_clip_intp2((int32_t)(_src[13] * scale + add) >> shift, MAX_LOG2_TR_RANGE + 1);
                _dst[2] = ov_clip_intp2((int32_t)(_src[14] * scale + add) >> shift, MAX_LOG2_TR_RANGE + 1);
                _dst[3] = ov_clip_intp2((int32_t)(_src[15] * scale + add) >> shift, MAX_LOG2_TR_RANGE + 1);
            } else {
                _dst[0] = 0;
                _dst[1] = 0;
                _dst[2] = 0;
                _dst[3] = 0;

                _dst += dst_stride;

                _dst[0] = 0;
                _dst[1] = 0;
                _dst[2] = 0;
                _dst[3] = 0;

                _dst += dst_stride;

                _dst[0] = 0;
                _dst[1] = 0;
                _dst[2] = 0;
                _dst[3] = 0;

                _dst += dst_stride;

                _dst[0] = 0;
                _dst[1] = 0;
                _dst[2] = 0;
                _dst[3] = 0;
            }
            sig_sb_row >>= 1;
        }
        src += src_stride << 2;
        dst += dst_stride << 2;
    }
}

static const uint16_t tb_lut_offset[64] =
{      0,    0,    2,    6,   14,   30,   62,   94,
      94,   96,  100,  108,  124,  156,  220,  284,
     284,  288,  296,  312,  344,  408,  536,  664,
     664,  672,  688,  720,  784,  912, 1168, 1424,
    1424, 1440, 1472, 1536, 1664, 1920, 2432, 2944,
    2944, 2976, 3040, 3168, 3424, 3936, 4960, 5984,
    5984, 6016, 6080, 6208, 6464, 6976, 8000, 9024,
    9024, 9024, 9024, 9024, 9024, 9024, 9024, 9024
};

static void
scale_dequant_tb_4x4(const int16_t *const luts, int16_t *dst, const int16_t *src, int scale, int shift,
               uint8_t log2_tb_w, uint8_t log2_tb_h, uint64_t sig_sb_map)
{
    int add = (1 << shift) >> 1;
    int nb_rows = 1 << OVMIN(5, log2_tb_h);//derive_nb_cols(sig_sb_map);
    int nb_cols = 1 << OVMIN(5, log2_tb_w);//derive_nb_rows(sig_sb_map);
    uint8_t src_stride = 1 << (OVMIN(5, log2_tb_w));
    uint8_t dst_stride = 1 << (OVMIN(5, log2_tb_w));

    /* Force sig_sb_map to one in case of DC coefficient */
    sig_sb_map |= !sig_sb_map;

    uint8_t lut_id = (log2_tb_w << 3) | log2_tb_h;
    const int16_t *sl = &luts[tb_lut_offset[lut_id]];
    scale /= 16;

    for (int i = 0; i < nb_rows/4 ; i++) {
        uint8_t sig_sb_row = sig_sb_map >> (i << 3);
        for (int j = 0; j < nb_cols/4 ; j++) {
            int16_t *_dst = dst + (j << 2);
            const int16_t *_src = src + (j << 4);
            if (sig_sb_row & 0x1) {
                _dst[0] = ov_clip_intp2((int32_t)((int32_t)_src[ 0] * sl[ 0] * scale + add) >> shift, MAX_LOG2_TR_RANGE + 1);
                _dst[1] = ov_clip_intp2((int32_t)((int32_t)_src[ 1] * sl[ 1] * scale + add) >> shift, MAX_LOG2_TR_RANGE + 1);
                _dst[2] = ov_clip_intp2((int32_t)((int32_t)_src[ 2] * sl[ 2] * scale + add) >> shift, MAX_LOG2_TR_RANGE + 1);
                _dst[3] = ov_clip_intp2((int32_t)((int32_t)_src[ 3] * sl[ 3] * scale + add) >> shift, MAX_LOG2_TR_RANGE + 1);

                _dst += dst_stride;

                _dst[0] = ov_clip_intp2((int32_t)((int32_t)_src[ 4] * sl[ 4] * scale + add) >> shift, MAX_LOG2_TR_RANGE + 1);
                _dst[1] = ov_clip_intp2((int32_t)((int32_t)_src[ 5] * sl[ 5] * scale + add) >> shift, MAX_LOG2_TR_RANGE + 1);
                _dst[2] = ov_clip_intp2((int32_t)((int32_t)_src[ 6] * sl[ 6] * scale + add) >> shift, MAX_LOG2_TR_RANGE + 1);
                _dst[3] = ov_clip_intp2((int32_t)((int32_t)_src[ 7] * sl[ 7] * scale + add) >> shift, MAX_LOG2_TR_RANGE + 1);

                _dst += dst_stride;

                _dst[0] = ov_clip_intp2((int32_t)((int32_t)_src[ 8] * sl[ 8] * scale + add) >> shift, MAX_LOG2_TR_RANGE + 1);
                _dst[1] = ov_clip_intp2((int32_t)((int32_t)_src[ 9] * sl[ 9] * scale + add) >> shift, MAX_LOG2_TR_RANGE + 1);
                _dst[2] = ov_clip_intp2((int32_t)((int32_t)_src[10] * sl[10] * scale + add) >> shift, MAX_LOG2_TR_RANGE + 1);
                _dst[3] = ov_clip_intp2((int32_t)((int32_t)_src[11] * sl[11] * scale + add) >> shift, MAX_LOG2_TR_RANGE + 1);

                _dst += dst_stride;

                _dst[0] = ov_clip_intp2((int32_t)((int32_t)_src[12] * sl[12] * scale + add) >> shift, MAX_LOG2_TR_RANGE + 1);
                _dst[1] = ov_clip_intp2((int32_t)((int32_t)_src[13] * sl[13] * scale + add) >> shift, MAX_LOG2_TR_RANGE + 1);
                _dst[2] = ov_clip_intp2((int32_t)((int32_t)_src[14] * sl[14] * scale + add) >> shift, MAX_LOG2_TR_RANGE + 1);
                _dst[3] = ov_clip_intp2((int32_t)((int32_t)_src[15] * sl[15] * scale + add) >> shift, MAX_LOG2_TR_RANGE + 1);
            } else {
                _dst[0] = 0;
                _dst[1] = 0;
                _dst[2] = 0;
                _dst[3] = 0;

                _dst += dst_stride;

                _dst[0] = 0;
                _dst[1] = 0;
                _dst[2] = 0;
                _dst[3] = 0;

                _dst += dst_stride;

                _dst[0] = 0;
                _dst[1] = 0;
                _dst[2] = 0;
                _dst[3] = 0;

                _dst += dst_stride;

                _dst[0] = 0;
                _dst[1] = 0;
                _dst[2] = 0;
                _dst[3] = 0;
            }
            sl += 16;
            sig_sb_row >>= 1;
        }
        src += src_stride << 2;
        dst += dst_stride << 2;
    }
}

void
BD_DECL(rcn_init_dequant_sdh)(struct RCNFunctions *rcn_funcs)
{
     rcn_funcs->tmp.derive_dequant     = &derive_dequant_sdh;
     rcn_funcs->tmp.derive_dequant_ts  = &derive_dequant_ts;
     rcn_funcs->tmp.dequant_tb_4x4 = &dequant_tb_4x4;
     rcn_funcs->tmp.scale_dequant_tb_4x4 = &scale_dequant_tb_4x4;
     //rcn_funcs->tmp.dequant_tb_4x4_neg = &dequant_tb_4x4_neg;
}

void
BD_DECL(rcn_init_dequant_dpq)(struct RCNFunctions *rcn_funcs)
{
     rcn_funcs->tmp.derive_dequant    = &derive_dequant_dpq;
     rcn_funcs->tmp.derive_dequant_ts = &derive_dequant_ts;
     rcn_funcs->tmp.dequant_tb_4x4 = &dequant_tb_4x4;
     rcn_funcs->tmp.scale_dequant_tb_4x4 = &scale_dequant_tb_4x4;
     //rcn_funcs->tmp.dequant_tb_4x4_neg = &dequant_tb_4x4_neg;
}
