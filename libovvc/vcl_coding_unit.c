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
#include "rcn_lmcs.h"
#include "dbf_utils.h"
#include "rcn_dequant.h"
#include "ovmv_structures.h"

#define LOG2_MIN_CU_S 2
/*FIXME find a more global location for these defintions */

#define BCW_NUM                 5 ///< the number of weight options
#define BCW_DEFAULT             ((uint8_t)(BCW_NUM >> 1)) ///< Default weighting index representing for w=0.5
#define BCW_SIZE_CONSTRAINT     256 ///< disabling Bcw if cu size is smaller than 256

/* skip_abv + skip_lft */
static uint8_t
ovcabac_read_ae_cu_skip_flag(OVCABACCtx *const cabac_ctx, uint8_t above_pu,
                             uint8_t left_pu)
{
    uint8_t cu_skip_flag;
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    int ctx_offset = ((above_pu == OV_IBC_SKIP) || (above_pu == OV_INTER_SKIP) || above_pu == OV_INTER_SKIP_AFFINE) + (left_pu == OV_INTER_SKIP || left_pu == OV_INTER_SKIP_AFFINE || (left_pu == OV_IBC_SKIP));
    cu_skip_flag = ovcabac_ae_read(cabac_ctx, &cabac_state[SKIP_FLAG_CTX_OFFSET + ctx_offset]);
    return cu_skip_flag;
}

static uint8_t
ovcabac_read_ae_cu_ibc_flag(OVCABACCtx *const cabac_ctx, uint8_t above_pu,
                             uint8_t left_pu)
{
    uint8_t cu_ibc_flag;
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    int ctx_offset = (above_pu == OV_IBC)  + (left_pu == OV_IBC) + (above_pu == OV_IBC_SKIP) + (left_pu == OV_IBC_SKIP);
    cu_ibc_flag = ovcabac_ae_read(cabac_ctx, &cabac_state[IBC_FLAG_CTX_OFFSET + ctx_offset]);
    return cu_ibc_flag;
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
    if (max_num_merge_cand > 1) {
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

static void
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

/* FIXME separate flag and idx */
static uint8_t
ovcabac_read_ae_bcw_flag(OVCABACCtx *const cabac_ctx)
{
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    uint8_t bcw_flag = ovcabac_ae_read(cabac_ctx, &cabac_state[BCW_IDX_CTX_OFFSET]);
    return bcw_flag;
}

static uint8_t
ovcabac_read_ae_bcw_idx(OVCABACCtx *const cabac_ctx, uint8_t is_ldc){
    static const int8_t parsing_order[BCW_NUM] = { 
        BCW_DEFAULT, BCW_DEFAULT + 1, BCW_DEFAULT - 1, BCW_DEFAULT + 2, BCW_DEFAULT - 2
    };

    int32_t nb_bcw_cand = is_ldc ? 5 : 3;
    uint32_t nb_idx_bits = nb_bcw_cand - 2;
    uint8_t bcw_idx = 1;
    int i;

    for (i = 0; i < nb_idx_bits; ++i) {
        if(!ovcabac_bypass_read(cabac_ctx)){
            break;
        }
        ++bcw_idx;
    }
    return parsing_order[bcw_idx];
}

static uint8_t
ovcabac_read_ae_ibc_amvr_precision(OVCABACCtx *const cabac_ctx)
{
    static const uint8_t amvr_precision[3] = {
        MV_PRECISION_INT, MV_PRECISION_INT, MV_PRECISION_4PEL,
    };

    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    uint8_t value = ovcabac_ae_read(cabac_ctx, &cabac_state[IMV_FLAG_CTX_OFFSET + 1]);
    uint8_t prec_idx = value + 1;
    return amvr_precision[prec_idx];
}

static uint8_t
ovcabac_read_ae_amvr_precision(OVCABACCtx *const cabac_ctx)
{
    static const uint8_t amvr_precision[4] = {
        MV_PRECISION_QUARTER, MV_PRECISION_INT, MV_PRECISION_4PEL, MV_PRECISION_HALF
    };

    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    uint8_t amvr_flag = ovcabac_ae_read(cabac_ctx, &cabac_state[IMV_FLAG_CTX_OFFSET]);
    if (amvr_flag) {
        uint8_t prec_idx = 3;
        uint8_t value = ovcabac_ae_read(cabac_ctx, &cabac_state[IMV_FLAG_CTX_OFFSET + 4]);

        if (value) {
            value = ovcabac_ae_read(cabac_ctx, &cabac_state[IMV_FLAG_CTX_OFFSET + 1]);
            prec_idx = value + 1;
        }
        return amvr_precision[prec_idx];
    }

    return amvr_precision[0];
}

static uint8_t
ovcabac_read_ae_affine_amvr_precision(OVCABACCtx *const cabac_ctx)
{
    static const uint8_t aff_amvr_precision[4] = {
        MV_PRECISION_QUARTER, MV_PRECISION_SIXTEENTH, MV_PRECISION_INT
    };

    uint64_t *const cabac_state = cabac_ctx->ctx_table;

    uint8_t prec_idx = ovcabac_ae_read(cabac_ctx, &cabac_state[IMV_FLAG_CTX_OFFSET + 2]);

    prec_idx += prec_idx && ovcabac_ae_read(cabac_ctx, &cabac_state[IMV_FLAG_CTX_OFFSET + 3]);

    return aff_amvr_precision[prec_idx];
}


static uint8_t
ovcabac_read_ae_inter_dir(OVCABACCtx *const cabac_ctx,
                          int log2_cb_w, int log2_cb_h)
{
    uint8_t inter_dir = 0;
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    /*FIXME Note no Bi predicition on 4x4 8x4 4x8*/
    if (log2_cb_w + log2_cb_h > 5) {
        int ctx_id = 7 - ((log2_cb_w + log2_cb_h + 1) >> 1);
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

    #if 0
    if (prefix < 0) {
        offset = prefix << 1;
    } else {
    #endif
        offset = (((1 << prefix ) - 1) << 1);
        length += (prefix == (32 - 17) ? 17 - 1 : prefix);
    #if 0
    }
    #endif

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

static uint8_t
ovcabac_read_ae_intra_bdpcm_flag(OVCABACCtx *const cabac_ctx)
{
    uint64_t *const cabac_state = cabac_ctx->ctx_table;

    return ovcabac_ae_read(cabac_ctx, &cabac_state[BDPCM_MODE_CTX_OFFSET]);
}

static uint8_t
ovcabac_read_ae_intra_bdpcm_dir(OVCABACCtx *const cabac_ctx)
{
    uint64_t *const cabac_state = cabac_ctx->ctx_table;

    return ovcabac_ae_read(cabac_ctx, &cabac_state[BDPCM_MODE_CTX_OFFSET + 1]);
}

static uint8_t
ovcabac_read_ae_intra_bdpcm_flag_c(OVCABACCtx *const cabac_ctx)
{
    uint64_t *const cabac_state = cabac_ctx->ctx_table;

    return ovcabac_ae_read(cabac_ctx, &cabac_state[BDPCM_MODE_CTX_OFFSET + 2]);
}

static uint8_t
ovcabac_read_ae_intra_bdpcm_dir_c(OVCABACCtx *const cabac_ctx)
{
    uint64_t *const cabac_state = cabac_ctx->ctx_table;

    return ovcabac_ae_read(cabac_ctx, &cabac_state[BDPCM_MODE_CTX_OFFSET + 3]);
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
            if (ref_idx + 1 < nb_active_ref)
            while (nb_active_ref > ref_idx + 1 && ovcabac_bypass_read(cabac_ctx)) {
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
    if (cu.cu_flags & flg_intra_bdpcm_chroma_flag) {
        intra_mode = (cu.cu_flags & flg_intra_bdpcm_chroma_dir) ? 50 : 18;
    }

    return intra_mode;
}

static int8_t
collocated_qp(const struct DBFQPMap *qp_map, int x0, int y0,
              int log2_cb_w, int log2_cb_h)
{
    uint8_t x0_u = x0 >> 2;
    uint8_t y0_u = y0 >> 2;
    int nb_units_w = (1 << (log2_cb_w - 1)) >> 2;
    int nb_units_h = (1 << (log2_cb_h - 1)) >> 2;

    x0_u += nb_units_w;
    y0_u += nb_units_h;

   return qp_map->hor[34 + 2 + x0_u + y0_u * 34];
}

int
coding_unit(OVCTUDec *const ctu_dec,
            const OVPartInfo *const part_ctx,
            uint8_t x0, uint8_t y0,
            uint8_t log2_cb_w, uint8_t log2_cb_h)
{
    uint8_t log2_min_cb_s = part_ctx->log2_min_cb_s;
    unsigned int nb_cb_w = (1 << log2_cb_w) >> log2_min_cb_s;
    unsigned int nb_cb_h = (1 << log2_cb_h) >> log2_min_cb_s;

    int x_cb = x0 >> log2_min_cb_s;
    int y_cb = y0 >> log2_min_cb_s;

    VVCCU cu;

    uint8_t compute_chr_scale = ((!(x0 & 0x3F) && !(y0 & 0x3F)) && ctu_dec->lmcs_info.lmcs_enabled_flag) ;
    if (compute_chr_scale){
        int16_t src_stride = ctu_dec->rcn_ctx.ctu_buff.stride;
        ctu_dec->rcn_funcs.rcn_lmcs_compute_chroma_scale(&ctu_dec->lmcs_info, src_stride, &ctu_dec->rcn_ctx.progress_field,
                                      ctu_dec->rcn_ctx.ctu_buff.y, x0, y0);
    }

    if (ctu_dec->coding_unit != &coding_unit_intra_c) {
        uint8_t cu_qp_delta_subdiv = ctu_dec->cu_qp_delta_subdiv;

        uint8_t qp_grp_msk = (1 << (part_ctx->log2_ctu_s + 1 - cu_qp_delta_subdiv))  - 1;
        uint8_t qp_grp_start = !((x0 & qp_grp_msk) | (y0 & qp_grp_msk));
        uint8_t x0_u = x0 >> 2;
        uint8_t y0_u = y0 >> 2;
        uint8_t nb_unit_w = (1 << log2_cb_w) >> 2;
        uint8_t nb_unit_h = (1 << log2_cb_h) >> 2;

        dbf_fill_cu_edge(&ctu_dec->dbf_info.cu_edge, x0_u, y0_u, nb_unit_w, nb_unit_h);
        if (qp_grp_start) {
            uint8_t qp_grp_x0 = x0 & (~qp_grp_msk);
            uint8_t qp_grp_y0 = y0 & (~qp_grp_msk);
            int pred_qp = ((y0 ? ctu_dec->drv_ctx.qp_map_x[qp_grp_x0 >> log2_min_cb_s] : ctu_dec->qp_ctx.current_qp) +
                           (x0 ? ctu_dec->drv_ctx.qp_map_y[qp_grp_y0 >> log2_min_cb_s] : ctu_dec->qp_ctx.current_qp) + 1) >> 1;
            ctu_dec->read_qp = 1;
            ctu_dec->qp_ctx.current_qp = pred_qp;
        }
        derive_dequant_ctx(ctu_dec, &ctu_dec->qp_ctx, 0);
    } else {
        int8_t qp_y = collocated_qp(&ctu_dec->dbf_info.qp_map_y, x0 << 1, y0 << 1, log2_cb_w + 1, log2_cb_h + 1);
        int base_qp_c = ov_clip(qp_y, 0, 63);
        int qp_bd_offset = ctu_dec->qp_ctx.qp_bd_offset;

        ctu_dec->dequant_cb.qp = ov_clip(ctu_dec->qp_ctx.chroma_qp_map_cb[base_qp_c] + ctu_dec->qp_ctx.dqp_cb + ctu_dec->qp_ctx.cb_offset, -qp_bd_offset, 63) + qp_bd_offset;
        ctu_dec->dequant_cr.qp = ov_clip(ctu_dec->qp_ctx.chroma_qp_map_cr[base_qp_c] + ctu_dec->qp_ctx.dqp_cr + ctu_dec->qp_ctx.cr_offset, -qp_bd_offset, 63) + qp_bd_offset;
        ctu_dec->dequant_joint_cb_cr.qp = ov_clip(ctu_dec->qp_ctx.chroma_qp_map_jcbcr[base_qp_c] + ctu_dec->qp_ctx.dqp_jcbcr + ctu_dec->qp_ctx.jcbcr_offset, -qp_bd_offset, 63) + qp_bd_offset;
        ctu_dec->dequant_cb_skip.qp   = OVMAX(ctu_dec->dequant_cb.qp, ctu_dec->qp_ctx.min_qp_prime_ts);
        ctu_dec->dequant_cr_skip.qp   = OVMAX(ctu_dec->dequant_cr.qp, ctu_dec->qp_ctx.min_qp_prime_ts);
        ctu_dec->dequant_jcbcr_skip.qp = OVMAX(ctu_dec->dequant_joint_cb_cr.qp, ctu_dec->qp_ctx.min_qp_prime_ts);
    }

    int qp_bd_offset = ctu_dec->qp_ctx.qp_bd_offset;

    ctu_dec->tmp_ciip = 0;

    cu = ctu_dec->coding_unit(ctu_dec, part_ctx, x0, y0, log2_cb_w, log2_cb_h);

    ctu_dec->dequant_chroma = &ctu_dec->dequant_cb;

    /* FIXME move after TU is read so we can reconstruct with or without
     * transform trees
     */

    if ((cu.cu_flags & flg_pred_mode_flag) && !(cu.cu_flags & flg_ibc_flag)) {
        if (ctu_dec->coding_unit == &coding_unit_intra_st || ctu_dec->coding_unit == &coding_unit_inter_st) {
            uint8_t luma_mode;
            luma_mode = drv_intra_cu(ctu_dec, part_ctx, x0, y0, log2_cb_w, log2_cb_h, cu);
            cu.cu_mode_idx = luma_mode;

            ctu_dec->intra_mode_c = drv_intra_mode_c(cu, luma_mode);
        } else {
            /* FIXME inter */
            if (ctu_dec->coding_unit == &coding_unit_intra) {
                uint8_t luma_mode = drv_intra_cu(ctu_dec, part_ctx, x0, y0, log2_cb_w, log2_cb_h, cu);
                cu.cu_mode_idx = luma_mode;
            } else {
                struct IntraDRVInfo *const i_info = &ctu_dec->drv_ctx.intra_info;
                uint8_t x0_unit = (x0 << 1) >> LOG2_MIN_CU_S;
                uint8_t y0_unit = (y0 << 1) >> LOG2_MIN_CU_S;
                uint8_t nb_unit_w = (2 << log2_cb_w) >> LOG2_MIN_CU_S;
                uint8_t nb_unit_h = (2 << log2_cb_h) >> LOG2_MIN_CU_S;

                uint8_t luma_mode = i_info->luma_modes[(x_cb + ((y_cb + (nb_cb_h >> 1)) << 5) + (nb_cb_w >> 1))];

                ctu_dec->intra_mode_c = drv_intra_mode_c(cu, luma_mode);
            }
        }
    }

    /*TODO rename */
    transform_unit_wrap(ctu_dec, part_ctx,  x0, y0, log2_cb_w, log2_cb_h, cu);

    /* FIXME memset instead
     */
    /* update delta qp context */
    if (ctu_dec->coding_unit != &coding_unit_intra_c) {

        for (int i = 0; i < nb_cb_w; i++) {
            ctu_dec->drv_ctx.qp_map_x[x_cb + i] = ctu_dec->qp_ctx.current_qp;
        }

        for (int i = 0; i < nb_cb_h; i++) {
            ctu_dec->drv_ctx.qp_map_y[y_cb + i] = ctu_dec->qp_ctx.current_qp;
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

static inline uint8_t
check_nz_mvd_p(const OVMV *const mvd)
{
    int32_t mvd_not_zero = 0;

    mvd_not_zero |= (mvd->x | mvd->y);

    return (uint8_t)!!mvd_not_zero;
}

static inline uint8_t
check_ibc_nz_mvd(const IBCMV *const mvd)
{
    int32_t mvd_not_zero = 0;

    mvd_not_zero |= (mvd->x | mvd->y);

    return (uint8_t)!!mvd_not_zero;
}

static struct IBCMVPData
inter_mvp_data_ibc(OVCTUDec *const ctu_dec, uint8_t nb_ibc_cand_min1)
{
    OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;
    struct IBCMVPData mvp_data;
    OVMV mvd;
    uint8_t mvp_idx = 0;

    mvd = ovcabac_read_ae_mvd(cabac_ctx);

    if (nb_ibc_cand_min1) {
        mvp_idx = ovcabac_read_ae_mvp_flag(cabac_ctx);
    }

    mvp_data.mvd.x   = mvd.x;
    mvp_data.mvd.y   = mvd.y;
    mvp_data.mvp_idx = mvp_idx;

    return mvp_data;
}

static struct MVPDataP
inter_mvp_data_p(OVCTUDec *const ctu_dec, uint8_t nb_active_ref_min1)
{
    OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;
    struct MVPDataP mvp_data;
    uint8_t ref_idx = nb_active_ref_min1;
    OVMV mvd;
    uint8_t mvp_idx;

    if (nb_active_ref_min1) {
        ref_idx = ovcabac_read_ae_ref_idx(cabac_ctx, nb_active_ref_min1 + 1);
    }

    mvd = ovcabac_read_ae_mvd(cabac_ctx);

    mvp_idx = ovcabac_read_ae_mvp_flag(cabac_ctx);


    mvp_data.ref_idx   = ref_idx;
    mvp_data.mvd       = mvd;
    mvp_data.mvp_idx   = mvp_idx;

    return mvp_data;
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

    uint8_t log2_min_cb_s = part_ctx->log2_min_cb_s;
    uint8_t x_cb = x0 >> log2_min_cb_s;
    uint8_t y_cb = y0 >> log2_min_cb_s;
    uint8_t nb_cb_w = (1 << log2_cu_w) >> log2_min_cb_s;
    uint8_t nb_cb_h = (1 << log2_cu_h) >> log2_min_cb_s;

    memset(&part_map->cu_mode_x[x_cb], (uint8_t)cu_mode, sizeof(uint8_t) * nb_cb_w);
    memset(&part_map->cu_mode_y[y_cb], (uint8_t)cu_mode, sizeof(uint8_t) * nb_cb_h);
}

static void
reset_intra_map_ibc(OVCTUDec *const ctudec, struct IntraDRVInfo *const intra_ctx,
                    uint8_t x0, uint8_t y0,
                    uint8_t log2_cb_w, uint8_t log2_cb_h,
                    uint8_t log2_min_cb_s)
{
    uint8_t x0_unit = x0 >> LOG2_MIN_CU_S;
    uint8_t y0_unit = y0 >> LOG2_MIN_CU_S;
    uint8_t nb_unit_w = (1 << log2_cb_w) >> LOG2_MIN_CU_S;
    uint8_t nb_unit_h = (1 << log2_cb_h) >> LOG2_MIN_CU_S;

    ctu_field_set_rect_bitfield(&ctudec->rcn_ctx.progress_field,
                                x0_unit, y0_unit, nb_unit_w, nb_unit_h);

    if (ctudec->coding_unit != coding_unit_intra) {
        ctu_field_set_rect_bitfield(&ctudec->rcn_ctx.progress_field_c,
                                    x0_unit, y0_unit, nb_unit_w, nb_unit_h);
    }

    uint8_t y_cb = y0 >> log2_min_cb_s;
    uint8_t x_cb = x0 >> log2_min_cb_s;
    uint8_t nb_cb_w = (1 << log2_cb_w) >> log2_min_cb_s;
    uint8_t nb_cb_h = (1 << log2_cb_h) >> log2_min_cb_s;

    /* We need to reset Intra mode maps to PLANAR for correct MPM derivation */
    memset(&intra_ctx->luma_mode_x[x_cb], OVINTRA_PLANAR, sizeof(uint8_t) * nb_cb_w);
    memset(&intra_ctx->luma_mode_y[y_cb], OVINTRA_PLANAR, sizeof(uint8_t) * nb_cb_h);

    for (int i = 0; i < nb_cb_h; i++) {
        memset(&intra_ctx->luma_modes[x_cb + (i << 5) + (y_cb << 5)], OVINTRA_DC,
               sizeof(uint8_t) * nb_cb_w);
    }
}

static void
reset_intra_map(OVCTUDec *const ctudec, struct IntraDRVInfo *const intra_ctx,
                uint8_t x0, uint8_t y0,
                uint8_t log2_cb_w, uint8_t log2_cb_h,
                uint8_t log2_min_cb_s)
{
    uint8_t x0_unit = x0 >> LOG2_MIN_CU_S;
    uint8_t y0_unit = y0 >> LOG2_MIN_CU_S;
    uint8_t nb_unit_w = (1 << log2_cb_w) >> LOG2_MIN_CU_S;
    uint8_t nb_unit_h = (1 << log2_cb_h) >> LOG2_MIN_CU_S;

    ctu_field_set_rect_bitfield(&ctudec->rcn_ctx.progress_field,
                                x0_unit, y0_unit, nb_unit_w, nb_unit_h);

    if (ctudec->coding_unit != coding_unit_intra) {
        ctu_field_set_rect_bitfield(&ctudec->rcn_ctx.progress_field_c,
                                    x0_unit, y0_unit, nb_unit_w, nb_unit_h);
    }

    uint8_t y_cb = y0 >> log2_min_cb_s;
    uint8_t x_cb = x0 >> log2_min_cb_s;
    uint8_t nb_cb_w = (1 << log2_cb_w) >> log2_min_cb_s;
    uint8_t nb_cb_h = (1 << log2_cb_h) >> log2_min_cb_s;

    /* We need to reset Intra mode maps to PLANAR for correct MPM derivation */
    memset(&intra_ctx->luma_mode_x[x_cb], OVINTRA_PLANAR, sizeof(uint8_t) * nb_cb_w);
    memset(&intra_ctx->luma_mode_y[y_cb], OVINTRA_PLANAR, sizeof(uint8_t) * nb_cb_h);

    for (int i = 0; i < nb_cb_h; i++) {
        memset(&intra_ctx->luma_modes[x_cb + (i << 5) + (y_cb << 5)], OVINTRA_PLANAR,
               sizeof(uint8_t) * nb_cb_w);
    }
}


VVCCU
coding_unit_inter_st(OVCTUDec *const ctu_dec,
                     const OVPartInfo *const part_ctx,
                     uint8_t x0, uint8_t y0,
                     uint8_t log2_cu_w, uint8_t log2_cu_h)
{
    OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;
    uint8_t log2_min_cb_s = part_ctx->log2_min_cb_s;
    uint8_t y_cb = y0 >> log2_min_cb_s;
    uint8_t x_cb = x0 >> log2_min_cb_s;
    /*TODO fill */
    uint8_t cu_type_abv = ctu_dec->part_map.cu_mode_x[x_cb];
    uint8_t cu_type_lft = ctu_dec->part_map.cu_mode_y[y_cb];
    uint8_t cu_type = OV_INTER;
    VVCCU cu = {0};

    uint8_t cu_skip_flag = ovcabac_read_ae_cu_skip_flag(cabac_ctx, cu_type_abv,
                                                cu_type_lft);

    if (cu_skip_flag) {
        FLG_STORE(cu_skip_flag, cu.cu_flags);
        if (ctu_dec->ibc_enabled && !ctu_dec->share && log2_cu_w < 7 && log2_cu_h < 7) {
            struct PartMap *part_map = &ctu_dec->part_map;
            uint8_t log2_min_cb_s = part_ctx->log2_min_cb_s;
            uint8_t x_cb = x0 >> log2_min_cb_s;
            uint8_t y_cb = y0 >> log2_min_cb_s;
            uint8_t nb_cb_w = (1 << log2_cu_w) >> log2_min_cb_s;
            uint8_t nb_cb_h = (1 << log2_cu_h) >> log2_min_cb_s;

            uint8_t ibc_abv = part_map->cu_mode_x[x_cb];
            uint8_t ibc_lft = part_map->cu_mode_y[y_cb];

            uint8_t ibc_flag = (log2_cu_w == 2 && log2_cu_h == 2) || ovcabac_read_ae_cu_ibc_flag(cabac_ctx, ibc_abv, ibc_lft);
            if (ibc_flag) {
                struct IntraDRVInfo *const i_info   = &ctu_dec->drv_ctx.intra_info;
                uint8_t nb_ibc_cand_min1 = ctu_dec->nb_ibc_cand_min1 - 1;
                uint8_t merge_flag = 1;
                IBCMV mv;
                uint8_t merge_idx = ovcabac_read_ae_mvp_merge_idx(cabac_ctx, nb_ibc_cand_min1 + 1);
                FLG_STORE(merge_flag, cu.cu_flags);

                FLG_STORE(ibc_flag, cu.cu_flags);

                if (!cu_skip_flag) {
                    updt_cu_maps(ctu_dec, part_ctx, x0, y0, log2_cu_w, log2_cu_h, OV_IBC);
                } else {
                    updt_cu_maps(ctu_dec, part_ctx, x0, y0, log2_cu_w, log2_cu_h, OV_IBC_SKIP);
                }
                mv = drv_ibc_merge_mv(&ctu_dec->drv_ctx.ibc_ctx, x0, y0, log2_cu_w, log2_cu_h, merge_idx, nb_ibc_cand_min1 + 1);
                ctu_dec->rcn_funcs.rcn_ibc_l(ctu_dec, x0, y0, log2_cu_w, log2_cu_h, part_ctx->log2_ctu_s, mv);
                if (ctu_dec->coding_unit != &coding_unit_intra && !ctu_dec->share) {
                    ctu_dec->rcn_funcs.rcn_ibc_c(ctu_dec, x0, y0, log2_cu_w, log2_cu_h, part_ctx->log2_ctu_s, mv);
                }

                reset_intra_map_ibc(ctu_dec, i_info, x0, y0, log2_cu_w, log2_cu_h, log2_min_cb_s);

                return cu;
            }
        }
        /* FIXME cu_skip_flag activation force merge_flag so we only need to read
           merge_idx */
        uint8_t merge_flag = 1;
        cu_type = ctu_dec->prediction_unit(ctu_dec, part_ctx, x0, y0, log2_cu_w, log2_cu_h, cu_skip_flag, 1);

        if (cu_type == OV_AFFINE) {
            cu_type = OV_INTER_SKIP_AFFINE;
        } else {
            cu_type = OV_INTER_SKIP;
        }

        FLG_STORE(cu_skip_flag, cu.cu_flags);
        ctu_dec->intra_mode_c = 0;

    } else {
        uint8_t pred_mode_flag = ctu_dec->share == 1 ? 1 : 0;

        if (ctu_dec->share == 1 || (log2_cu_w == 2 && log2_cu_h == 2)) {
            if (ctu_dec->ibc_enabled && log2_cu_w < 7 && log2_cu_h < 7) {
                struct PartMap *part_map = &ctu_dec->part_map;
                uint8_t log2_min_cb_s = part_ctx->log2_min_cb_s;
                uint8_t x_cb = x0 >> log2_min_cb_s;
                uint8_t y_cb = y0 >> log2_min_cb_s;
                uint8_t nb_cb_w = (1 << log2_cu_w) >> log2_min_cb_s;
                uint8_t nb_cb_h = (1 << log2_cu_h) >> log2_min_cb_s;

                uint8_t ibc_abv = part_map->cu_mode_x[x_cb];
                uint8_t ibc_lft = part_map->cu_mode_y[y_cb];

                uint8_t ibc_flag = ovcabac_read_ae_cu_ibc_flag(cabac_ctx, ibc_abv, ibc_lft);
                if (ibc_flag) {
                    IBCMV mv;
                    struct IntraDRVInfo *const i_info   = &ctu_dec->drv_ctx.intra_info;
                    uint8_t nb_ibc_cand_min1 = ctu_dec->nb_ibc_cand_min1 - 1;
                    uint8_t merge_flag = ovcabac_read_ae_cu_merge_flag(cabac_ctx);
                    if (merge_flag) {
                        uint8_t merge_idx = ovcabac_read_ae_mvp_merge_idx(cabac_ctx, nb_ibc_cand_min1 + 1);
                        FLG_STORE(merge_flag, cu.cu_flags);
                        mv = drv_ibc_merge_mv(&ctu_dec->drv_ctx.ibc_ctx, x0, y0, log2_cu_w, log2_cu_h, merge_idx, nb_ibc_cand_min1 + 1);
                    } else {
                        struct IBCMVPData mvp_data = inter_mvp_data_ibc(ctu_dec, nb_ibc_cand_min1);
                        uint8_t prec_amvr = MV_PRECISION_INT;

                        if (ctu_dec->amvr_enabled) {
                            uint8_t nz_mvd = check_ibc_nz_mvd(&mvp_data.mvd);
                            if (nz_mvd) {
                                prec_amvr = ovcabac_read_ae_ibc_amvr_precision(cabac_ctx);
                            }
                        }

                        mv = drv_ibc_mvp(&ctu_dec->drv_ctx.ibc_ctx, x0, y0, log2_cu_w, log2_cu_h, mvp_data.mvd, mvp_data.mvp_idx, prec_amvr);

                    }

                    ctu_dec->rcn_funcs.rcn_ibc_l(ctu_dec, x0, y0, log2_cu_w, log2_cu_h, part_ctx->log2_ctu_s, mv);
                    if (ctu_dec->coding_unit != &coding_unit_intra && !ctu_dec->share) {
                        ctu_dec->rcn_funcs.rcn_ibc_c(ctu_dec, x0, y0, log2_cu_w, log2_cu_h, part_ctx->log2_ctu_s, mv);
                    }

                    reset_intra_map_ibc(ctu_dec, i_info, x0, y0, log2_cu_w, log2_cu_h, log2_min_cb_s);

                    FLG_STORE(ibc_flag, cu.cu_flags);

                    updt_cu_maps(ctu_dec, part_ctx, x0, y0, log2_cu_w, log2_cu_h, OV_IBC);

                    return cu;
                }
            }
        } else if (!ctu_dec->share) {
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
            if (ctu_dec->ibc_enabled && !ctu_dec->share && log2_cu_w < 7 && log2_cu_h < 7) {
                struct PartMap *part_map = &ctu_dec->part_map;
                uint8_t log2_min_cb_s = part_ctx->log2_min_cb_s;
                uint8_t x_cb = x0 >> log2_min_cb_s;
                uint8_t y_cb = y0 >> log2_min_cb_s;
                uint8_t nb_cb_w = (1 << log2_cu_w) >> log2_min_cb_s;
                uint8_t nb_cb_h = (1 << log2_cu_h) >> log2_min_cb_s;

                uint8_t ibc_abv = part_map->cu_mode_x[x_cb];
                uint8_t ibc_lft = part_map->cu_mode_y[y_cb];

                uint8_t ibc_flag = ovcabac_read_ae_cu_ibc_flag(cabac_ctx, ibc_abv, ibc_lft);

                if (ibc_flag) {
                    uint8_t nb_ibc_cand_min1 = ctu_dec->nb_ibc_cand_min1 - 1;
                    uint8_t merge_flag = ovcabac_read_ae_cu_merge_flag(cabac_ctx);
                    struct IntraDRVInfo *const i_info   = &ctu_dec->drv_ctx.intra_info;
                    IBCMV mv;
                    if (merge_flag) {
                        uint8_t merge_idx = ovcabac_read_ae_mvp_merge_idx(cabac_ctx, nb_ibc_cand_min1 + 1);
                        FLG_STORE(merge_flag, cu.cu_flags);
                        mv = drv_ibc_merge_mv(&ctu_dec->drv_ctx.ibc_ctx, x0, y0, log2_cu_w, log2_cu_h, merge_idx, nb_ibc_cand_min1 + 1);
                    } else {
                        struct IBCMVPData mvp_data = inter_mvp_data_ibc(ctu_dec, nb_ibc_cand_min1);
                        uint8_t prec_amvr = MV_PRECISION_INT;
                        if (ctu_dec->amvr_enabled) {
                            uint8_t nz_mvd = check_ibc_nz_mvd(&mvp_data.mvd);
                            if (nz_mvd) {
                                prec_amvr = ovcabac_read_ae_ibc_amvr_precision(cabac_ctx);
                            }
                        }
                        mv = drv_ibc_mvp(&ctu_dec->drv_ctx.ibc_ctx, x0, y0, log2_cu_w, log2_cu_h, mvp_data.mvd, mvp_data.mvp_idx, prec_amvr);
                    }

                    ctu_dec->rcn_funcs.rcn_ibc_l(ctu_dec, x0, y0, log2_cu_w, log2_cu_h, part_ctx->log2_ctu_s, mv);
                    if (ctu_dec->coding_unit != &coding_unit_intra && !ctu_dec->share) {
                        ctu_dec->rcn_funcs.rcn_ibc_c(ctu_dec, x0, y0, log2_cu_w, log2_cu_h, part_ctx->log2_ctu_s, mv);
                    }

                    reset_intra_map_ibc(ctu_dec, i_info, x0, y0, log2_cu_w, log2_cu_h, log2_min_cb_s);

                    if (!cu_skip_flag) {
                        updt_cu_maps(ctu_dec, part_ctx, x0, y0, log2_cu_w, log2_cu_h, OV_IBC);
                    } else {
                        updt_cu_maps(ctu_dec, part_ctx, x0, y0, log2_cu_w, log2_cu_h, OV_IBC_SKIP);
                        FLG_STORE(cu_skip_flag, cu.cu_flags);
                    }
                    FLG_STORE(ibc_flag, cu.cu_flags);
                    return cu;
                }
            }
            uint8_t merge_flag = ovcabac_read_ae_cu_merge_flag(cabac_ctx);
            cu_type = ctu_dec->prediction_unit(ctu_dec, part_ctx, x0, y0, log2_cu_w, log2_cu_h, 0, merge_flag);

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
   if (!ctu_dec->share && !(cu.cu_flags & flg_ibc_flag)) {
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

    /* Check IBC only if we are in luma tree. This means we come either from dual tree or separate tree */
    if (ctu_dec->ibc_enabled && (ctu_dec->coding_unit == &coding_unit_intra ||
                                 ctu_dec->coding_unit == &coding_unit_intra_st) &&
                                 log2_cb_w < 7 && log2_cb_h < 7) {

        struct PartMap *part_map = &ctu_dec->part_map;
        uint8_t log2_min_cb_s = part_ctx->log2_min_cb_s;
        uint8_t x_cb = x0 >> log2_min_cb_s;
        uint8_t y_cb = y0 >> log2_min_cb_s;

        uint8_t ibc_abv = part_map->cu_mode_x[x_cb];
        uint8_t ibc_lft = part_map->cu_mode_y[y_cb];

        uint8_t cu_skip_flag = ovcabac_read_ae_cu_skip_flag(cabac_ctx, ibc_abv, ibc_lft);

        uint8_t ibc_flag = cu_skip_flag || ovcabac_read_ae_cu_ibc_flag(cabac_ctx, ibc_abv, ibc_lft);

        if (ibc_flag) {
            struct IntraDRVInfo *const i_info   = &ctu_dec->drv_ctx.intra_info;
            uint8_t nb_ibc_cand_min1 = ctu_dec->nb_ibc_cand_min1 - 1;
            uint8_t merge_flag = cu_skip_flag || ovcabac_read_ae_cu_merge_flag(cabac_ctx);
            IBCMV mv;
            if (merge_flag) {
                uint8_t merge_idx = ovcabac_read_ae_mvp_merge_idx(cabac_ctx, nb_ibc_cand_min1 + 1);
                FLG_STORE(merge_flag, cu.cu_flags);
                mv = drv_ibc_merge_mv(&ctu_dec->drv_ctx.ibc_ctx, x0, y0, log2_cb_w, log2_cb_h, merge_idx, nb_ibc_cand_min1 + 1);
            } else {
                struct IBCMVPData mvp_data = inter_mvp_data_ibc(ctu_dec, nb_ibc_cand_min1);
                uint8_t prec_amvr = MV_PRECISION_INT;
                if (ctu_dec->amvr_enabled) {
                    uint8_t nz_mvd = check_ibc_nz_mvd(&mvp_data.mvd);
                    if (nz_mvd) {
                        prec_amvr = ovcabac_read_ae_ibc_amvr_precision(cabac_ctx);
                    }
                }
                mv = drv_ibc_mvp(&ctu_dec->drv_ctx.ibc_ctx, x0, y0, log2_cb_w, log2_cb_h, mvp_data.mvd, mvp_data.mvp_idx, prec_amvr);
            }

            ctu_dec->rcn_funcs.rcn_ibc_l(ctu_dec, x0, y0, log2_cb_w, log2_cb_h, part_ctx->log2_ctu_s, mv);
            if (ctu_dec->coding_unit != &coding_unit_intra && !ctu_dec->share) {
                ctu_dec->rcn_funcs.rcn_ibc_c(ctu_dec, x0, y0, log2_cb_w, log2_cb_h, part_ctx->log2_ctu_s, mv);
            }

            reset_intra_map_ibc(ctu_dec, i_info, x0, y0, log2_cb_w, log2_cb_h, log2_min_cb_s);

            if (!cu_skip_flag) {
                updt_cu_maps(ctu_dec, part_ctx, x0, y0, log2_cb_w, log2_cb_h, OV_IBC);
            } else {
                updt_cu_maps(ctu_dec, part_ctx, x0, y0, log2_cb_w, log2_cb_h, OV_IBC_SKIP);
                FLG_STORE(cu_skip_flag, cu.cu_flags);
            }
            FLG_STORE(ibc_flag, cu.cu_flags);

            goto end;
        }
    }

    cu.cu_flags |= flg_pred_mode_flag;

    if (ctu_dec->bdpcm_enabled && log2_cb_w <= ctu_dec->max_log2_transform_skip_size
                               && log2_cb_h <= ctu_dec->max_log2_transform_skip_size) {
        uint8_t intra_bdpcm_luma_flag = ovcabac_read_ae_intra_bdpcm_flag(cabac_ctx);
        if (intra_bdpcm_luma_flag) {
            struct PartMap *part_map = &ctu_dec->part_map;
            uint8_t log2_min_cb_s = part_ctx->log2_min_cb_s;
            uint8_t x_cb = x0 >> log2_min_cb_s;
            uint8_t y_cb = y0 >> log2_min_cb_s;
            uint8_t nb_cb_w = (1 << log2_cb_w) >> log2_min_cb_s;
            uint8_t nb_cb_h = (1 << log2_cb_h) >> log2_min_cb_s;

            uint8_t intra_bdpcm_luma_dir = ovcabac_read_ae_intra_bdpcm_dir(cabac_ctx);

            memset(&part_map->cu_mode_x[x_cb], OV_INTRA, sizeof(uint8_t) * nb_cb_w);
            memset(&part_map->cu_mode_y[y_cb], OV_INTRA, sizeof(uint8_t) * nb_cb_h);

            FLG_STORE(intra_bdpcm_luma_flag, cu.cu_flags);
            FLG_STORE(intra_bdpcm_luma_dir, cu.cu_flags);
            updt_cu_maps(ctu_dec, part_ctx, x0, y0, log2_cb_w, log2_cb_h, OV_INTRA);
            goto end;
        }
    }

    if (ctu_dec->enabled_mip) {
        struct PartMap *part_map = &ctu_dec->part_map;
        uint8_t log2_min_cb_s = part_ctx->log2_min_cb_s;
        uint8_t x_cb = x0 >> log2_min_cb_s;
        uint8_t y_cb = y0 >> log2_min_cb_s;
        uint8_t nb_cb_w = (1 << log2_cb_w) >> log2_min_cb_s;
        uint8_t nb_cb_h = (1 << log2_cb_h) >> log2_min_cb_s;
        uint8_t mip_abv = part_map->cu_mode_x[x_cb];
        uint8_t mip_lft = part_map->cu_mode_y[y_cb];

        mip_flag = ovcabac_read_ae_intra_mip(cabac_ctx, log2_cb_w, log2_cb_h,
                                             mip_abv, mip_lft);

        cu.cu_flags |= flg_mip_flag & (-(!!mip_flag));

        if (mip_flag) {
            uint8_t transpose_flag = ovcabac_read_ae_intra_mip_transpose_flag(cabac_ctx);
            uint8_t mip_mode;

            mip_mode = ovcabac_read_ae_intra_mip_mode(cabac_ctx, log2_cb_w, log2_cb_h);

            cu.cu_opaque  = transpose_flag << 7;
            cu.cu_opaque |= mip_mode;

            memset(&part_map->cu_mode_x[x_cb], OV_MIP, sizeof(uint8_t) * nb_cb_w);
            memset(&part_map->cu_mode_y[y_cb], OV_MIP, sizeof(uint8_t) * nb_cb_h);

        } else {

            memset(&part_map->cu_mode_x[x_cb], OV_INTRA, sizeof(uint8_t) * nb_cb_w);
            memset(&part_map->cu_mode_y[y_cb], OV_INTRA, sizeof(uint8_t) * nb_cb_h);
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
            isp_split_status &= !(log2_cb_w > part_ctx->log2_max_tb_s || log2_cb_h > part_ctx->log2_max_tb_s);
            if (isp_split_status) {
                isp_split_status  = (log2_cb_w <= part_ctx->log2_max_tb_s) << 1;
                isp_split_status |=  log2_cb_h <= part_ctx->log2_max_tb_s;
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

        updt_cu_maps(ctu_dec, part_ctx, x0, y0, log2_cb_w, log2_cb_h, OV_INTRA);
    }
end:

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

    /* Force intra pred mode in case of separate or dual tree */
    cu.cu_flags |= flg_pred_mode_flag;

    if (ctu_dec->bdpcm_enabled && log2_cb_w <= ctu_dec->max_log2_transform_skip_size
                               && log2_cb_h <= ctu_dec->max_log2_transform_skip_size) {
        uint8_t intra_bdpcm_chroma_flag = ovcabac_read_ae_intra_bdpcm_flag_c(cabac_ctx);
        if (intra_bdpcm_chroma_flag) {
            struct PartMap *part_map = &ctu_dec->part_map;
            uint8_t log2_min_cb_s = part_ctx->log2_min_cb_s;
            uint8_t x_cb = x0 >> log2_min_cb_s;
            uint8_t y_cb = y0 >> log2_min_cb_s;
            uint8_t nb_cb_w = (1 << log2_cb_w) >> log2_min_cb_s;
            uint8_t nb_cb_h = (1 << log2_cb_h) >> log2_min_cb_s;

            uint8_t intra_bdpcm_chroma_dir = ovcabac_read_ae_intra_bdpcm_dir_c(cabac_ctx);

            FLG_STORE(intra_bdpcm_chroma_flag, cu.cu_flags);
            FLG_STORE(intra_bdpcm_chroma_dir, cu.cu_flags);

            goto end;
        }
    }

    /* FIXME CCLM luma partition constraints */
    if (ctu_dec->lm_chroma_enabled && ((!ctu_dec->tmp_disable_cclm &&
        ctu_dec->enable_cclm == 1 && ctu_dec->coding_tree == &dual_tree) ||
        ctu_dec->coding_tree != &dual_tree || part_ctx->log2_ctu_s == 5)) {

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

end:
    return cu;
}

static OVMV
drv_merge_motion_info_p(const struct MergeData *const mrg_data,
                        struct InterDRVCtx *const inter_ctx,
                        uint8_t x0, uint8_t y0,
                        uint8_t log2_cb_w, uint8_t log2_cb_h,
                        uint8_t log2_min_cb_s, uint8_t max_nb_cand)
{
    OVMV mv;
    enum MergeTypeP merge_type   = mrg_data->merge_type;
    uint8_t merge_idx = mrg_data->merge_idx;
    struct OVMVCtx *const mv_ctx = &inter_ctx->mv_ctx0;

    switch (merge_type) {
        case MMVD_MERGE:
            mv = drv_mmvd_merge_mvp(inter_ctx, mv_ctx,
                                    x0, y0, log2_cb_w, log2_cb_h,
                                    merge_idx, max_nb_cand);
            break;
        default:
            /* Note CIIP use classical MV derivation */
            mv = drv_merge_mvp(inter_ctx, mv_ctx,
                               x0, y0, log2_cb_w, log2_cb_h,
                               merge_idx, max_nb_cand);
            break;
    }
    return mv;
}

static struct MergeData
inter_skip_data_p(OVCTUDec *const ctu_dec,
                  const OVPartInfo *const part_ctx,
                  uint8_t x0, uint8_t y0,
                  uint8_t log2_cb_w, uint8_t log2_cb_h)
{
    const struct InterDRVCtx *const inter_ctx = &ctu_dec->drv_ctx.inter_ctx;
    OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;
    struct MergeData mrg_data;
    enum MergeTypeP mrg_type = DEFAULT_MERGE;

    uint8_t sb_merge_flag = 0;

    if ((ctu_dec->sbtmvp_enabled || ctu_dec->affine_enabled)
        && log2_cb_w >= 3 && log2_cb_h >= 3) {
        uint8_t log2_min_cb_s = part_ctx->log2_min_cb_s;
        uint8_t y_cb = y0 >> log2_min_cb_s;
        uint8_t x_cb = x0 >> log2_min_cb_s;
        uint8_t cu_type_abv = ctu_dec->part_map.cu_mode_x[x_cb];
        uint8_t cu_type_lft = ctu_dec->part_map.cu_mode_y[y_cb];

        uint8_t lft_affine = cu_type_lft == OV_AFFINE || cu_type_lft == OV_INTER_SKIP_AFFINE;
        uint8_t abv_affine = cu_type_abv == OV_AFFINE || cu_type_abv == OV_INTER_SKIP_AFFINE;

        sb_merge_flag = ovcabac_read_ae_sb_merge_flag(cabac_ctx, lft_affine, abv_affine);
    }

    uint8_t mmvd_enabled = ctu_dec->mmvd_enabled;
    uint8_t mmvd_flag = mmvd_enabled && !sb_merge_flag && ovcabac_read_ae_mmvd_flag(cabac_ctx);
    uint8_t merge_idx;

    if (sb_merge_flag) {
        if (ctu_dec->affine_enabled) {
            uint8_t nb_affine_merge_cand_min1 = ctu_dec->affine_nb_merge_cand - 1;
            merge_idx = ovcabac_read_ae_affine_merge_idx(cabac_ctx, nb_affine_merge_cand_min1);
        } else {
            merge_idx = 0;
        }
        mrg_type = SB_MERGE;
    } else if (mmvd_flag){
        uint8_t max_nb_cand = ctu_dec->max_num_merge_candidates;
        merge_idx = ovcabac_read_ae_mmvd_merge_idx(cabac_ctx, max_nb_cand);
        mrg_type = MMVD_MERGE;
    } else {
        uint8_t max_nb_cand = ctu_dec->max_num_merge_candidates;
        merge_idx = ovcabac_read_ae_mvp_merge_idx(cabac_ctx, max_nb_cand);
    }

    mrg_data.merge_type = mrg_type;
    mrg_data.merge_idx  = merge_idx;

    return mrg_data;
}

static struct MergeData
inter_merge_data_p(OVCTUDec *const ctu_dec,
                   const OVPartInfo *const part_ctx,
                   uint8_t x0, uint8_t y0,
                   uint8_t log2_cb_w, uint8_t log2_cb_h)
{
    OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;
    struct MergeData mrg_data;
    enum MergeTypeP mrg_type = DEFAULT_MERGE;
    uint8_t merge_idx;

    uint8_t max_nb_cand = ctu_dec->max_num_merge_candidates;

    /* FIXME missing affine in P */
    uint8_t sps_ciip_flag = ctu_dec->ciip_enabled;
    uint8_t ciip_enabled = sps_ciip_flag && log2_cb_w < 7
                                         && log2_cb_h < 7
                                         && (log2_cb_w + log2_cb_h) >= 6;

    uint8_t sb_merge_flag = 0;

    if ((ctu_dec->sbtmvp_enabled || ctu_dec->affine_enabled)
        && log2_cb_w >= 3 && log2_cb_h >= 3) {
        uint8_t log2_min_cb_s = part_ctx->log2_min_cb_s;
        uint8_t y_cb = y0 >> log2_min_cb_s;
        uint8_t x_cb = x0 >> log2_min_cb_s;
        uint8_t cu_type_abv = ctu_dec->part_map.cu_mode_x[x_cb];
        uint8_t cu_type_lft = ctu_dec->part_map.cu_mode_y[y_cb];

        uint8_t lft_affine = cu_type_lft == OV_AFFINE || cu_type_lft == OV_INTER_SKIP_AFFINE;
        uint8_t abv_affine = cu_type_abv == OV_AFFINE || cu_type_abv == OV_INTER_SKIP_AFFINE;

        sb_merge_flag = ovcabac_read_ae_sb_merge_flag(cabac_ctx, lft_affine, abv_affine);
    }

    uint8_t mmvd_flag  = 0;
    if (!sb_merge_flag) {
        uint8_t reg_merge_flag = !ciip_enabled || ovcabac_read_ae_reg_merge_flag(cabac_ctx, 0);

        if (reg_merge_flag) {
            uint8_t mmvd_enabled = ctu_dec->mmvd_enabled;
            mmvd_flag = mmvd_enabled && ovcabac_read_ae_mmvd_flag(cabac_ctx);
        } else {
            mrg_type = CIIP_MERGE;
            ctu_dec->tmp_ciip = 1;
        }
    }

    if (sb_merge_flag) {
        if (ctu_dec->affine_enabled) {
            uint8_t nb_affine_merge_cand_min1 = ctu_dec->affine_nb_merge_cand - 1;
            merge_idx = ovcabac_read_ae_affine_merge_idx(cabac_ctx, nb_affine_merge_cand_min1);
        } else {
            merge_idx = 0;
        }
        mrg_type = SB_MERGE;
    } else if (mmvd_flag){
        merge_idx = ovcabac_read_ae_mmvd_merge_idx(cabac_ctx, max_nb_cand);
        mrg_type = MMVD_MERGE;
    } else {
        merge_idx = ovcabac_read_ae_mvp_merge_idx(cabac_ctx, max_nb_cand);
    }

    mrg_data.merge_type = mrg_type;
    mrg_data.merge_idx  = merge_idx;

    return mrg_data;
}

static inline uint8_t
check_nz_affine_mvd_p(const struct AffineControlInfo *const cp_mvd, uint8_t affine_type)
{
    uint32_t mvd_not_zero = 0;
    if (affine_type) {
        mvd_not_zero |= (cp_mvd->lt.x | cp_mvd->lt.y);
        mvd_not_zero |= (cp_mvd->rt.x | cp_mvd->rt.y);
        mvd_not_zero |= (cp_mvd->lb.x | cp_mvd->lb.y);
    } else {
        mvd_not_zero |= (cp_mvd->lt.x | cp_mvd->lt.y);
        mvd_not_zero |= (cp_mvd->rt.x | cp_mvd->rt.y);
    }

    return (uint8_t)!!mvd_not_zero;
}

static struct AffineMVPDataP
inter_affine_mvp_data_p(OVCTUDec *const ctu_dec, uint8_t nb_active_ref_min1, uint8_t affine_type)
{
    struct AffineControlInfo cp_mvd;
    struct AffineMVPDataP mvp_data;
    OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;
    uint8_t ref_idx = nb_active_ref_min1;
    uint8_t mvp_idx;

    if (nb_active_ref_min1) {
        ref_idx = ovcabac_read_ae_ref_idx(cabac_ctx, nb_active_ref_min1 + 1);
    }

    cp_mvd.lt = ovcabac_read_ae_mvd(cabac_ctx);
    cp_mvd.rt = ovcabac_read_ae_mvd(cabac_ctx);

    if (affine_type) {
        cp_mvd.lb = ovcabac_read_ae_mvd(cabac_ctx);
    }

    mvp_idx = ovcabac_read_ae_mvp_flag(cabac_ctx);

    mvp_data.mvd       = cp_mvd;

    mvp_data.ref_idx   = ref_idx;
    mvp_data.mvp_idx   = mvp_idx;

    return mvp_data;
}

static struct AffineMVPDataB
inter_affine_mvp_data_b(OVCTUDec *const ctu_dec, uint8_t nb_active_ref0_min1,
                        uint8_t nb_active_ref1_min1, uint8_t affine_type)
{
    OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;
    struct AffineMVPDataB mvp_data;

    mvp_data.ref_idx0 = nb_active_ref0_min1;
    if (nb_active_ref0_min1) {
        mvp_data.ref_idx0 = ovcabac_read_ae_ref_idx(cabac_ctx, nb_active_ref0_min1 + 1);
    }

    mvp_data.mvd0.lt = ovcabac_read_ae_mvd(cabac_ctx);
    mvp_data.mvd0.rt = ovcabac_read_ae_mvd(cabac_ctx);
    if (affine_type) {
        mvp_data.mvd0.lb = ovcabac_read_ae_mvd(cabac_ctx);
    }

    mvp_data.mvp_idx0 = ovcabac_read_ae_mvp_flag(cabac_ctx);

    mvp_data.ref_idx1 = nb_active_ref1_min1;
    if (nb_active_ref1_min1) {
        mvp_data.ref_idx1 = ovcabac_read_ae_ref_idx(cabac_ctx, nb_active_ref1_min1 + 1);
    }

    if (!ctu_dec->mvd1_zero_enabled) {
        mvp_data.mvd1.lt = ovcabac_read_ae_mvd(cabac_ctx);
        mvp_data.mvd1.rt = ovcabac_read_ae_mvd(cabac_ctx);
        if (affine_type) {
            mvp_data.mvd1.lb = ovcabac_read_ae_mvd(cabac_ctx);
        }
    } else {
        memset(&mvp_data.mvd1, 0, sizeof(mvp_data.mvd1));
    }

    mvp_data.mvp_idx1 = ovcabac_read_ae_mvp_flag(cabac_ctx);

    return mvp_data;
}

static struct MVPDataB
inter_mvp_data_b(OVCTUDec *const ctu_dec, uint8_t nb_active_ref0_min1,
                 uint8_t nb_active_ref1_min1)
{
    OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;
    struct MVPDataB mvp_data;
    OVMV mvd0;
    OVMV mvd1 = {0};
    uint8_t ref_idx0 = 0;
    uint8_t mvp_idx0 = 0;
    uint8_t ref_idx1 = 0;
    uint8_t mvp_idx1 = 0;

    if (nb_active_ref0_min1) {
        ref_idx0 = ovcabac_read_ae_ref_idx(cabac_ctx, nb_active_ref0_min1 + 1);
    }

    mvd0 = ovcabac_read_ae_mvd(cabac_ctx);

    mvp_idx0 = ovcabac_read_ae_mvp_flag(cabac_ctx);

    if (nb_active_ref1_min1) {
        ref_idx1 = ovcabac_read_ae_ref_idx(cabac_ctx, nb_active_ref1_min1 + 1);
    }

    if (!ctu_dec->mvd1_zero_enabled) {
        mvd1 = ovcabac_read_ae_mvd(cabac_ctx);
    }

    mvp_idx1 = ovcabac_read_ae_mvp_flag(cabac_ctx);

    mvp_data.mvd0 = mvd0;
    mvp_data.mvd1 = mvd1;

    mvp_data.ref_idx0 = ref_idx0;
    mvp_data.mvp_idx0 = mvp_idx0;
    mvp_data.ref_idx1 = ref_idx1;
    mvp_data.mvp_idx1 = mvp_idx1;

    return mvp_data;
}

static struct MVPInfoP
inter_mvp_read_p(OVCTUDec *const ctu_dec,
                 uint8_t x0, uint8_t y0,
                 uint8_t log2_cb_w, uint8_t log2_cb_h,
                 uint8_t log2_min_cb_s,
                 uint8_t nb_active_ref_min1)
{
    OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;
    struct MVPInfoP mvp_info;
    uint8_t affine_flag = 0;

    if (ctu_dec->affine_enabled && log2_cb_w > 3 && log2_cb_h > 3) {
        uint8_t y_cb = y0 >> log2_min_cb_s;
        uint8_t x_cb = x0 >> log2_min_cb_s;
        uint8_t cu_type_abv = ctu_dec->part_map.cu_mode_x[x_cb];
        uint8_t cu_type_lft = ctu_dec->part_map.cu_mode_y[y_cb];

        uint8_t lft_affine = cu_type_lft == OV_AFFINE || cu_type_lft == OV_INTER_SKIP_AFFINE;
        uint8_t abv_affine = cu_type_abv == OV_AFFINE || cu_type_abv == OV_INTER_SKIP_AFFINE;

        affine_flag = ovcabac_read_ae_cu_affine_flag(cabac_ctx, lft_affine, abv_affine);
    }

    if (affine_flag) {
        struct AffineMVPDataP mvp_data;
        uint8_t prec_amvr = MV_PRECISION_QUARTER;
        uint8_t six_affine_type = !!(ctu_dec->affine_status & 0x2);
        uint8_t affine_type = six_affine_type && ovcabac_read_ae_cu_affine_type(cabac_ctx);

        mvp_data = inter_affine_mvp_data_p(ctu_dec, nb_active_ref_min1, affine_type);

        if (ctu_dec->affine_amvr_enabled) {
            int32_t nz_mvd = check_nz_affine_mvd_p(&mvp_data.mvd, affine_type);

            if (nz_mvd) {
                prec_amvr = ovcabac_read_ae_affine_amvr_precision(cabac_ctx);
            }
        }

        mvp_info.cu_type = OV_AFFINE;
        mvp_info.data.cp_mvd_info = mvp_data;
        mvp_info.prec_amvr = prec_amvr;
        mvp_info.affine_type = affine_type;

    } else {
        struct MVPDataP mvp_data;
        uint8_t prec_amvr = MV_PRECISION_QUARTER;

        mvp_data = inter_mvp_data_p(ctu_dec, nb_active_ref_min1);

        if (ctu_dec->amvr_enabled) {
            uint8_t nz_mvd = check_nz_mvd_p(&mvp_data.mvd);
            if (nz_mvd) {
                prec_amvr = ovcabac_read_ae_amvr_precision(cabac_ctx);
            }
        }

        mvp_info.cu_type = OV_INTER;
        mvp_info.data.mvd_info = mvp_data;
        mvp_info.prec_amvr = prec_amvr;
    }

    return mvp_info;
}

static int
drv_rcn_wrap_mvp_p(OVCTUDec *const ctu_dec,
                   const OVPartInfo *const part_ctx,
                   uint8_t x0, uint8_t y0,
                   uint8_t log2_cb_w, uint8_t log2_cb_h,
                   struct MVPInfoP mvp_info, uint8_t inter_dir)
{
    struct InterDRVCtx *const inter_ctx = &ctu_dec->drv_ctx.inter_ctx;
    struct IntraDRVInfo *const i_info   = &ctu_dec->drv_ctx.intra_info;
    uint8_t log2_min_cb_s = part_ctx->log2_min_cb_s;
    VVCMergeInfo mv_info;
    uint8_t ref_idx0 = 0;
    uint8_t ref_idx1 = 0;
    uint8_t cu_type = OV_INTER;
    if (mvp_info.cu_type == OV_AFFINE) {
        struct AffineMVPDataP *const mvp_data = &mvp_info.data.cp_mvd_info;

        drv_affine_mvp_p(inter_ctx, x0, y0, log2_cb_w, log2_cb_h,
                         &mvp_info, inter_dir);

        cu_type = OV_AFFINE;

        goto end;

    } else {
        const struct MVPDataP *const mvp_data = &mvp_info.data.mvd_info;

        if (inter_dir & 0x1) {
            uint8_t mvp_idx0 = mvp_data->mvp_idx;
            ref_idx0 = mvp_data->ref_idx;
            mv_info = drv_mvp_b(inter_ctx, x0, y0, log2_cb_w, log2_cb_h,
                                mvp_data->mvd, mvp_data->mvd,
                                mvp_info.prec_amvr,
                                mvp_idx0, -1,
                                BCW_DEFAULT,
                                0x1, ref_idx0, -1,
                                log2_cb_w + log2_cb_h <= 5);

            mv_info.mv0.prec_amvr = mvp_info.prec_amvr;
        } else {
            uint8_t mvp_idx1 = mvp_data->mvp_idx;
            ref_idx1 = mvp_data->ref_idx;
            mv_info = drv_mvp_b(inter_ctx, x0, y0, log2_cb_w, log2_cb_h,
                                mvp_data->mvd, mvp_data->mvd,
                                mvp_info.prec_amvr,
                                -1, mvp_idx1,
                                BCW_DEFAULT,
                                0x2, -1, ref_idx1,
                                log2_cb_w + log2_cb_h <= 5);

            mv_info.mv1.prec_amvr = mvp_info.prec_amvr;
        }
    }

    ctu_dec->rcn_funcs.rcn_mcp_b(ctu_dec, ctu_dec->rcn_ctx.ctu_buff, inter_ctx, part_ctx,
                                 mv_info.mv0, mv_info.mv1, x0, y0,
                                 log2_cb_w, log2_cb_h, inter_dir, ref_idx0, ref_idx1);

end:
    reset_intra_map(ctu_dec, i_info, x0, y0, log2_cb_w, log2_cb_h, log2_min_cb_s);

    return cu_type;
}

int
prediction_unit_inter_p(OVCTUDec *const ctu_dec,
                        const OVPartInfo *const part_ctx,
                        uint8_t x0, uint8_t y0,
                        uint8_t log2_cb_w, uint8_t log2_cb_h,
                        uint8_t skip_flag, uint8_t merge_flag)
{
    struct IntraDRVInfo *const i_info = &ctu_dec->drv_ctx.intra_info;
#if 1
    struct InterDRVCtx *const inter_ctx = &ctu_dec->drv_ctx.inter_ctx;
#endif

    uint8_t cu_type = OV_INTER;

    uint8_t log2_min_cb_s = part_ctx->log2_min_cb_s;
    uint8_t ref_idx = 0;

    ctu_dec->tmp_ciip = 0;

    OVMV mv0;
    if (merge_flag) {
        uint8_t max_nb_cand = ctu_dec->max_num_merge_candidates;
        struct MergeData mrg_data;
        if (skip_flag) {
            mrg_data = inter_skip_data_p(ctu_dec, part_ctx, x0, y0,
                                         log2_cb_w, log2_cb_h);
        } else {
            mrg_data = inter_merge_data_p(ctu_dec, part_ctx, x0, y0,
                                          log2_cb_w, log2_cb_h);
        }

    if (mrg_data.merge_type == SB_MERGE) {

        /* Note reconstruction is done in drv_function */
        drv_affine_merge_mvp_p(inter_ctx, x0, y0, log2_cb_w, log2_cb_h,
                               mrg_data.merge_idx);

        cu_type = OV_AFFINE;

        goto end;

    } else  {
        mv0 =  drv_merge_motion_info_p(&mrg_data, inter_ctx, x0, y0,
                                       log2_cb_w, log2_cb_h,
                                       log2_min_cb_s,
                                       max_nb_cand);
        ref_idx = mv0.ref_idx;

        if (mrg_data.merge_type == CIIP_MERGE) {
            ctu_dec->rcn_funcs.rcn_ciip(ctu_dec, x0, y0, log2_cb_w, log2_cb_h, mv0, ref_idx);
            goto end;
        }
    }

    } else {
        struct MVPInfoP mvp_info = inter_mvp_read_p(ctu_dec, x0, y0, log2_cb_w, log2_cb_h, log2_min_cb_s,
                                                    inter_ctx->nb_active_ref0 - 1);
        cu_type = drv_rcn_wrap_mvp_p(ctu_dec, part_ctx, x0, y0,
                                     log2_cb_w, log2_cb_h, mvp_info, 0x1);
        goto end;
    }

    ctu_dec->rcn_funcs.rcn_mcp(ctu_dec, ctu_dec->rcn_ctx.ctu_buff, x0, y0,
            log2_cb_w, log2_cb_h, mv0, 0, ref_idx);

end:
    /*FIXME this have to be moved to DRV */
    reset_intra_map(ctu_dec, i_info, x0, y0, log2_cb_w, log2_cb_h, log2_min_cb_s);

    return cu_type;
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
    uint8_t use_weights = inter_ctx->wp_info0[ref_idx0].flag || inter_ctx->wp_info1[ref_idx1].flag;
    return !use_weights && inter_ctx->dist_ref_0[ref_idx0] == -inter_ctx->dist_ref_1[ref_idx1];
}

static inline uint8_t
check_nz_affine_b(const struct AffineControlInfo *const cp_mvd0, const struct AffineControlInfo *const cp_mvd1, uint8_t affine_type)
{
    uint32_t mvd_not_zero = 0;
    if (affine_type) {
        mvd_not_zero |= (cp_mvd0->lt.x | cp_mvd0->lt.y);
        mvd_not_zero |= (cp_mvd0->rt.x | cp_mvd0->rt.y);
        mvd_not_zero |= (cp_mvd0->lb.x | cp_mvd0->lb.y);
        if (!mvd_not_zero) {
            mvd_not_zero |= (cp_mvd1->lt.x | cp_mvd1->lt.y);
            mvd_not_zero |= (cp_mvd1->rt.x | cp_mvd1->rt.y);
            mvd_not_zero |= (cp_mvd1->lb.x | cp_mvd1->lb.y);
        }
    } else {
        mvd_not_zero |= (cp_mvd0->lt.x | cp_mvd0->lt.y);
        mvd_not_zero |= (cp_mvd0->rt.x | cp_mvd0->rt.y);
        if (!mvd_not_zero) {
            mvd_not_zero |= (cp_mvd1->lt.x | cp_mvd1->lt.y);
            mvd_not_zero |= (cp_mvd1->rt.x | cp_mvd1->rt.y);
        }
    }

    return (uint8_t)!!mvd_not_zero;
}

static inline uint8_t
check_nz_mvd(const OVMV *const mvd0, const OVMV *const mvd1,
             uint8_t inter_dir, uint8_t smvd_flag, uint8_t mvd1_zero_flag)
{
    int32_t mvd_not_zero = 0;
    if (inter_dir & 0x1) {
        mvd_not_zero |= (mvd0->x | mvd0->y);
        if (inter_dir & 0x2 && !smvd_flag) {
            if (!mvd1_zero_flag) {
                mvd_not_zero |= (mvd1->x | mvd1->y);
            }
        }
    } else {
        mvd_not_zero |= (mvd1->x | mvd1->y);
    }

    return (uint8_t)!!mvd_not_zero;
}

static inline uint8_t
check_nz_mvd_b(const OVMV *const mvd0, const OVMV *const mvd1,
               uint8_t mvd1_zero_flag)
{
    uint32_t mvd_not_zero = 0;
    mvd_not_zero |= (mvd0->x | mvd0->y);
    if (!mvd1_zero_flag) {
        mvd_not_zero |= (mvd1->x | mvd1->y);
    }

    return (uint8_t)!!mvd_not_zero;
}

static inline uint8_t
check_nz_mvd_smvd(const OVMV *const mvd0)
{
    int32_t mvd_not_zero = 0;
    mvd_not_zero |= (mvd0->x | mvd0->y);
    return (uint8_t)!!mvd_not_zero;
}

static struct MergeData
inter_skip_data_b(OVCTUDec *const ctu_dec,
                  const OVPartInfo *const part_ctx,
                  uint8_t x0, uint8_t y0,
                  uint8_t log2_cb_w, uint8_t log2_cb_h)
{
    const struct InterDRVCtx *const inter_ctx = &ctu_dec->drv_ctx.inter_ctx;
    OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;
    struct MergeData mrg_data;
    enum MergeTypeP mrg_type = DEFAULT_MERGE;

    uint8_t sb_merge_flag = 0;
    uint8_t mmvd_flag = 0;
    uint8_t merge_idx;
    uint8_t reg_merge_flag;

    uint8_t gpm_enabled  = ctu_dec->gpm_enabled && ctu_dec->max_gpm_cand > 1
                                                && log2_cb_w > 2 && log2_cb_h > 2
                                                && log2_cb_w < 7 && log2_cb_h < 7
                                                && log2_cb_w < 3 + log2_cb_h
                                                && log2_cb_h < 3 + log2_cb_w;

    if ((ctu_dec->sbtmvp_enabled || ctu_dec->affine_enabled) &&
        log2_cb_w >= 3 && log2_cb_h >= 3) {
        uint8_t log2_min_cb_s = part_ctx->log2_min_cb_s;
        uint8_t y_cb = y0 >> log2_min_cb_s;
        uint8_t x_cb = x0 >> log2_min_cb_s;
        uint8_t cu_type_abv = ctu_dec->part_map.cu_mode_x[x_cb];
        uint8_t cu_type_lft = ctu_dec->part_map.cu_mode_y[y_cb];

        uint8_t lft_affine = cu_type_lft == OV_AFFINE || cu_type_lft == OV_INTER_SKIP_AFFINE;
        uint8_t abv_affine = cu_type_abv == OV_AFFINE || cu_type_abv == OV_INTER_SKIP_AFFINE;

        sb_merge_flag = ovcabac_read_ae_sb_merge_flag(cabac_ctx, lft_affine, abv_affine);
        if (sb_merge_flag) {
            goto idx;
        }
    }

    reg_merge_flag = !gpm_enabled || ovcabac_read_ae_reg_merge_flag(cabac_ctx, 1);

    if (reg_merge_flag) {
        uint8_t mmvd_enabled = ctu_dec->mmvd_enabled;
        mmvd_flag = mmvd_enabled && ovcabac_read_ae_mmvd_flag(cabac_ctx);
    }

idx:
    if (sb_merge_flag) {
        if (ctu_dec->affine_enabled) {
            uint8_t nb_affine_merge_cand_min1 = ctu_dec->affine_nb_merge_cand - 1;
            merge_idx = ovcabac_read_ae_affine_merge_idx(cabac_ctx, nb_affine_merge_cand_min1);
        } else {
            merge_idx = 0;
        }
        mrg_type = SB_MERGE;
    } else if (!reg_merge_flag) {
        int max_num_gpm_cand = ctu_dec->max_gpm_cand;
        struct VVCGPM *gpm_info = &ctu_dec->drv_ctx.inter_ctx.gpm_ctx;
        /* FIXME gpm idx type */
        ovcabac_read_ae_gpm_merge_idx(cabac_ctx, gpm_info, max_num_gpm_cand);
        mrg_type = GPM_MERGE;
        /* set to what was read in gpm merge idx to discard uninitialized value */
        merge_idx = gpm_info->merge_idx0;
    } else if (mmvd_flag){
        uint8_t max_nb_cand = ctu_dec->max_num_merge_candidates;
        merge_idx = ovcabac_read_ae_mmvd_merge_idx(cabac_ctx, max_nb_cand);
        mrg_type = MMVD_MERGE;
    } else {
        uint8_t max_nb_cand = ctu_dec->max_num_merge_candidates;
        merge_idx = ovcabac_read_ae_mvp_merge_idx(cabac_ctx, max_nb_cand);
    }

    mrg_data.merge_type = mrg_type;
    mrg_data.merge_idx  = merge_idx;

    return mrg_data;
}

static struct MergeData
inter_merge_data_b(OVCTUDec *const ctu_dec,
                   const OVPartInfo *const part_ctx,
                   uint8_t x0, uint8_t y0,
                   uint8_t log2_cb_w, uint8_t log2_cb_h)
{
    const struct InterDRVCtx *const inter_ctx = &ctu_dec->drv_ctx.inter_ctx;
    OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;
    struct MergeData mrg_data;
    enum MergeTypeP mrg_type = DEFAULT_MERGE;
    uint8_t merge_idx;

    uint8_t ciip_flag;

    uint8_t ciip_enabled = ctu_dec->ciip_enabled && log2_cb_w < 7 && log2_cb_h < 7
                                                 && (log2_cb_w + log2_cb_h) >= 6;

    uint8_t gpm_enabled  = ctu_dec->gpm_enabled && ctu_dec->max_gpm_cand > 1
                                                && log2_cb_w > 2 && log2_cb_h > 2
                                                && log2_cb_w < 7 && log2_cb_h < 7
                                                && log2_cb_w < 3 + log2_cb_h
                                                && log2_cb_h < 3 + log2_cb_w;

    uint8_t sb_merge_flag = 0;
    uint8_t mmvd_flag = 0;

    if ((ctu_dec->sbtmvp_enabled || ctu_dec->affine_enabled)
        && log2_cb_w >= 3 && log2_cb_h >= 3) {
        uint8_t log2_min_cb_s = part_ctx->log2_min_cb_s;
        uint8_t y_cb = y0 >> log2_min_cb_s;
        uint8_t x_cb = x0 >> log2_min_cb_s;
        uint8_t cu_type_abv = ctu_dec->part_map.cu_mode_x[x_cb];
        uint8_t cu_type_lft = ctu_dec->part_map.cu_mode_y[y_cb];

        uint8_t lft_affine = cu_type_lft == OV_AFFINE || cu_type_lft == OV_INTER_SKIP_AFFINE;
        uint8_t abv_affine = cu_type_abv == OV_AFFINE || cu_type_abv == OV_INTER_SKIP_AFFINE;

        sb_merge_flag = ovcabac_read_ae_sb_merge_flag(cabac_ctx, lft_affine, abv_affine);
        if (sb_merge_flag) {
            goto idx;
        }
    }

    uint8_t reg_merge_flag = !(gpm_enabled || ciip_enabled);

    reg_merge_flag = reg_merge_flag || ovcabac_read_ae_reg_merge_flag(cabac_ctx, 0);

    if (reg_merge_flag) {
        uint8_t mmvd_enabled = ctu_dec->mmvd_enabled;
        mmvd_flag = mmvd_enabled && ovcabac_read_ae_mmvd_flag(cabac_ctx);
    } else {
        ciip_flag = ciip_enabled;
        if (gpm_enabled && ciip_enabled) {
            ciip_flag = ovcabac_read_ae_ciip_flag(cabac_ctx);
        }
    }

idx:
    if (sb_merge_flag) {
        if (ctu_dec->affine_enabled) {
            uint8_t nb_affine_merge_cand_min1 = ctu_dec->affine_nb_merge_cand - 1;
            merge_idx = ovcabac_read_ae_affine_merge_idx(cabac_ctx, nb_affine_merge_cand_min1);
        } else {
            merge_idx = 0;
        }
        mrg_type = SB_MERGE;
    } else if (!reg_merge_flag && !ciip_flag) {
        int max_num_gpm_cand = ctu_dec->max_gpm_cand;
        struct VVCGPM *gpm_info = &ctu_dec->drv_ctx.inter_ctx.gpm_ctx;
        /* FIXME gpm idx type */
        ovcabac_read_ae_gpm_merge_idx(cabac_ctx, gpm_info, max_num_gpm_cand);
        mrg_type = GPM_MERGE;
    } else if (!reg_merge_flag && ciip_flag) {
        uint8_t max_nb_cand = ctu_dec->max_num_merge_candidates;
        merge_idx = ovcabac_read_ae_mvp_merge_idx(cabac_ctx, max_nb_cand);
        mrg_type = CIIP_MERGE;
    } else if (mmvd_flag){
        uint8_t max_nb_cand = ctu_dec->max_num_merge_candidates;
        merge_idx = ovcabac_read_ae_mmvd_merge_idx(cabac_ctx, max_nb_cand);
        mrg_type = MMVD_MERGE;
    } else {
        uint8_t max_nb_cand = ctu_dec->max_num_merge_candidates;
        merge_idx = ovcabac_read_ae_mvp_merge_idx(cabac_ctx, max_nb_cand);
    }

    mrg_data.merge_type = mrg_type;
    mrg_data.merge_idx  = merge_idx;

    return mrg_data;
}

uint8_t read_bidir_mvp(OVCTUDec *const ctu_dec,
                       uint8_t x0, uint8_t y0,
                       uint8_t log2_cb_w, uint8_t log2_cb_h,
                       uint8_t log2_min_cb_s,
                       uint8_t nb_active_ref0_min1, uint8_t nb_active_ref1_min1)
{
    OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;
    struct InterDRVCtx *const inter_ctx = &ctu_dec->drv_ctx.inter_ctx;
    struct MVPInfoB mvp_info;

    uint8_t cu_type = OV_INTER;

    uint8_t affine_flag = 0;
    uint8_t smvd_flag = 0;

    if (ctu_dec->affine_enabled && log2_cb_w > 3 && log2_cb_h > 3) {
        uint8_t y_cb = y0 >> log2_min_cb_s;
        uint8_t x_cb = x0 >> log2_min_cb_s;
        uint8_t cu_type_abv = ctu_dec->part_map.cu_mode_x[x_cb];
        uint8_t cu_type_lft = ctu_dec->part_map.cu_mode_y[y_cb];

        uint8_t lft_affine = cu_type_lft == OV_AFFINE || cu_type_lft == OV_INTER_SKIP_AFFINE;
        uint8_t abv_affine = cu_type_abv == OV_AFFINE || cu_type_abv == OV_INTER_SKIP_AFFINE;

        affine_flag = ovcabac_read_ae_cu_affine_flag(cabac_ctx, lft_affine, abv_affine);

        if (affine_flag) {
            uint8_t six_affine_type = !!(ctu_dec->affine_status & 0x2);
            uint8_t affine_type = six_affine_type && ovcabac_read_ae_cu_affine_type(cabac_ctx);
            uint8_t prec_amvr = MV_PRECISION_QUARTER;

            struct AffineMVPDataB mvp_data = inter_affine_mvp_data_b(ctu_dec, nb_active_ref0_min1, nb_active_ref1_min1, affine_type);

            /* Note affine is always be 1 here  skip_flag always 0 */
            if (ctu_dec->affine_amvr_enabled) {
                int32_t nz_mvd = check_nz_affine_b(&mvp_data.mvd0, &mvp_data.mvd1,
                                                   affine_type);

                if (nz_mvd) {
                    prec_amvr = ovcabac_read_ae_affine_amvr_precision(cabac_ctx);
                }
            }

            mvp_info.data.aff_mvp = mvp_data;
            mvp_info.prec_amvr = prec_amvr;
            mvp_info.affine_type = affine_type;
        }
    }

    if (!affine_flag) {

        uint8_t prec_amvr = MV_PRECISION_QUARTER;
        if (ctu_dec->smvd_enabled) {
            smvd_flag = ovcabac_read_ae_smvd_flag(cabac_ctx);
        }

        if (smvd_flag) {
            struct SMVDData smvd;

            smvd.mvd = ovcabac_read_ae_mvd(cabac_ctx);

            smvd.mvp_idx0 = ovcabac_read_ae_mvp_flag(cabac_ctx);
            smvd.mvp_idx1 = ovcabac_read_ae_mvp_flag(cabac_ctx);

            if (ctu_dec->amvr_enabled) {
                uint8_t nz_mvd = check_nz_mvd_smvd(&smvd.mvd);
                if (nz_mvd) {
                    prec_amvr = ovcabac_read_ae_amvr_precision(cabac_ctx);
                }
            }

            mvp_info.data.smvd = smvd;

        } else {

            struct MVPDataB mvp_b = inter_mvp_data_b(ctu_dec, nb_active_ref0_min1,
                                                     nb_active_ref1_min1);
            if (ctu_dec->amvr_enabled) {
                uint8_t nz_mvd = check_nz_mvd_b(&mvp_b.mvd0, &mvp_b.mvd1, ctu_dec->mvd1_zero_enabled);
                if (nz_mvd) {
                    prec_amvr = ovcabac_read_ae_amvr_precision(cabac_ctx);
                }
            }

            mvp_info.data.mvp = mvp_b;
        }

        mvp_info.prec_amvr = prec_amvr;
    }

    uint8_t bcw_idx = BCW_DEFAULT;
    if (ctu_dec->bcw_enabled && (1 << (log2_cb_h + log2_cb_w) >= BCW_SIZE_CONSTRAINT)) {
        uint8_t bcw_flag = ovcabac_read_ae_bcw_flag(cabac_ctx);
        if (bcw_flag) {
            bcw_idx = ovcabac_read_ae_bcw_idx(cabac_ctx, inter_ctx->tmvp_ctx.ldc);
        }
    }

    mvp_info.bcw_idx = bcw_idx;

    if (affine_flag) {
        const struct AffineMVPDataB *const mvp_data = &mvp_info.data.aff_mvp;

        drv_affine_mvp_b(inter_ctx, x0, y0, log2_cb_w, log2_cb_h,
                         &mvp_info, 0x3);

        cu_type = OV_AFFINE;

        return cu_type;

    } else {
        VVCMergeInfo mv_info;

        OVMV mvd0, mvd1 = {0};
        uint8_t mvp_idx0 = 0;
        uint8_t mvp_idx1 = 0;
        uint8_t ref_idx0 = 0;
        uint8_t ref_idx1 = 0;

        if (smvd_flag) {

            const struct SMVDData *const smvd = &mvp_info.data.smvd;
            mvd0 = smvd->mvd;
            mvp_idx0 = smvd->mvp_idx0;
            mvp_idx1 = smvd->mvp_idx1;

            ref_idx0     = inter_ctx->ref_smvd_idx0;
            ref_idx1     = inter_ctx->ref_smvd_idx1;
            mvd1.x       = -mvd0.x;
            mvd1.y       = -mvd0.y;
            /* FIXME check if necessary */
            mvd1.ref_idx = inter_ctx->ref_smvd_idx1;

        } else {
            const struct MVPDataB *const mvp_b = &mvp_info.data.mvp;

            mvd0 = mvp_b->mvd0;
            ref_idx0 = mvp_b->ref_idx0;
            mvp_idx0 = mvp_b->mvp_idx0;

            mvd1 = mvp_b->mvd1;
            ref_idx1 = mvp_b->ref_idx1;
            mvp_idx1 = mvp_b->mvp_idx1;
        }

        mv_info = drv_mvp_b(inter_ctx, x0, y0, log2_cb_w, log2_cb_h,
                            mvd0, mvd1, mvp_info.prec_amvr, mvp_idx0, mvp_idx1, bcw_idx,
                            0x3, ref_idx0, ref_idx1, log2_cb_w + log2_cb_h <= 5);

        uint8_t bdof_enable = 0;
        if (inter_ctx->bdof_enabled) {
            /* Note ciip_flag is zero in this function */
            uint8_t ciip_flag = 0;
            uint8_t bcw_flag = (mv_info.mv0.bcw_idx_plus1 != 0 && mv_info.mv0.bcw_idx_plus1 != 3);
            bdof_enable = check_bdof(log2_cb_w, log2_cb_h, ciip_flag, bcw_flag, smvd_flag);

            bdof_enable = bdof_enable && check_bdof_ref(inter_ctx, ref_idx0, ref_idx1);
        }

        if (!bdof_enable) {
            ctu_dec->rcn_funcs.rcn_mcp_b(ctu_dec, ctu_dec->rcn_ctx.ctu_buff, inter_ctx, ctu_dec->part_ctx,
                      mv_info.mv0, mv_info.mv1, x0, y0,
                      log2_cb_w, log2_cb_h, mv_info.inter_dir, ref_idx0, ref_idx1);
        } else {
            uint8_t log2_w = OVMIN(log2_cb_w, 4);
            uint8_t log2_h = OVMIN(log2_cb_h, 4);
            uint8_t nb_sb_w = (1 << log2_cb_w) >> log2_w;
            uint8_t nb_sb_h = (1 << log2_cb_h) >> log2_h;
            int i, j;

            OVMV mv0 = mv_info.mv0;
            OVMV mv1 = mv_info.mv1;
            for (i = 0; i < nb_sb_h; ++i) {
                for (j = 0; j < nb_sb_w; ++j) {
                    ctu_dec->rcn_funcs.rcn_bdof_mcp_l(ctu_dec, ctu_dec->rcn_ctx.ctu_buff,
                                   x0 + j * 16, y0 + i * 16, log2_w, log2_h,
                                   mv0, mv1, ref_idx0, ref_idx1);
                }
            }
            ctu_dec->rcn_funcs.rcn_mcp_b_c(ctu_dec, ctu_dec->rcn_ctx.ctu_buff, inter_ctx, ctu_dec->part_ctx,
                        mv0, mv1, x0, y0,
                        log2_cb_w, log2_cb_h, mv_info.inter_dir, ref_idx0, ref_idx1);
        }
    }
    return cu_type;
}

int
prediction_unit_inter_b(OVCTUDec *const ctu_dec,
                        const OVPartInfo *const part_ctx,
                        uint8_t x0, uint8_t y0,
                        uint8_t log2_cb_w, uint8_t log2_cb_h,
                        uint8_t skip_flag, uint8_t merge_flag)
{
    OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;
    #if 1
    struct InterDRVCtx *const inter_ctx = &ctu_dec->drv_ctx.inter_ctx;
    struct IntraDRVInfo *const i_info   = &ctu_dec->drv_ctx.intra_info;
    VVCMergeInfo mv_info;
    #endif
    uint8_t log2_min_cb_s = part_ctx->log2_min_cb_s;
    uint8_t ref_idx0 = 0;
    uint8_t ref_idx1 = 0;
    uint8_t cu_type = OV_INTER;

    uint8_t smvd_flag = 0;
    uint8_t mmvd_flag = 0;

    /* FIXME Move AMVR precision outside of inter_ctx */
    ctu_dec->tmp_ciip = 0;

    if (merge_flag) {
        struct MergeData mrg_data;
        if (skip_flag) {
            mrg_data = inter_skip_data_b(ctu_dec, part_ctx, x0, y0,
                                         log2_cb_w, log2_cb_h);
        } else {
            mrg_data = inter_merge_data_b(ctu_dec, part_ctx, x0, y0,
                                          log2_cb_w, log2_cb_h);
        }

        mmvd_flag = mrg_data.merge_type == MMVD_MERGE;

        if (mrg_data.merge_type == SB_MERGE) {

            /* Note reconstruction is done in drv_function */
            drv_affine_merge_mvp_b(inter_ctx, x0, y0, log2_cb_w, log2_cb_h,
                                   mrg_data.merge_idx);

            cu_type = OV_AFFINE;

            goto end;

        } else if (mrg_data.merge_type == GPM_MERGE) {
            uint8_t max_nb_cand = ctu_dec->max_num_merge_candidates;

            drv_gpm_merge_mvp_b(inter_ctx, x0, y0, log2_cb_w, log2_cb_h, max_nb_cand,
                                log2_cb_w + log2_cb_h <= 5);

            ctu_dec->rcn_funcs.rcn_gpm_b(ctu_dec, &inter_ctx->gpm_ctx, x0, y0, log2_cb_w, log2_cb_h);

            goto end;

        } else if (mrg_data.merge_type == CIIP_MERGE) {
            uint8_t max_nb_cand = ctu_dec->max_num_merge_candidates;

            mv_info = drv_merge_mvp_b(inter_ctx, x0, y0,
                                      log2_cb_w, log2_cb_h, mrg_data.merge_idx,
                                      max_nb_cand, log2_cb_w + log2_cb_h <= 5);

            fill_bs_map(&ctu_dec->dbf_info.bs2_map, x0, y0, log2_cb_w, log2_cb_h);
            fill_bs_map(&ctu_dec->dbf_info.bs2_map_c, x0, y0, log2_cb_w, log2_cb_h);

            ref_idx0 = mv_info.mv0.ref_idx;
            ref_idx1 = mv_info.mv1.ref_idx;

            mv_info.mv0.bcw_idx_plus1 = 0;
            mv_info.mv1.bcw_idx_plus1 = 0;

            ctu_dec->rcn_funcs.rcn_ciip_b(ctu_dec, mv_info.mv0, mv_info.mv1, x0, y0,
                       log2_cb_w, log2_cb_h, mv_info.inter_dir, ref_idx0, ref_idx1);

            ctu_dec->tmp_ciip = 1;

            goto end;

        } else if (mrg_data.merge_type == MMVD_MERGE){
            uint8_t max_nb_cand = ctu_dec->max_num_merge_candidates;

            mv_info = drv_mmvd_merge_mvp_b(inter_ctx, x0, y0,
                                           log2_cb_w, log2_cb_h, mrg_data.merge_idx,
                                           max_nb_cand, log2_cb_w + log2_cb_h <= 5);

        } else {
            uint8_t max_nb_cand = ctu_dec->max_num_merge_candidates;
            mv_info = drv_merge_mvp_b(inter_ctx, x0, y0,
                                      log2_cb_w, log2_cb_h, mrg_data.merge_idx,
                                      max_nb_cand, log2_cb_w + log2_cb_h <= 5);
        }

        ref_idx0 = mv_info.mv0.ref_idx;
        ref_idx1 = mv_info.mv1.ref_idx;

        {
            uint8_t bdof_enable = 0;
            uint8_t dmvr_enable = 0;
            if (inter_ctx->bdof_enabled && mv_info.inter_dir == 0x3) {
                /* Note ciip_flag is zero in this function */
                uint8_t ciip_flag = 0;
                uint8_t bcw_flag = (mv_info.mv0.bcw_idx_plus1 != 0 && mv_info.mv0.bcw_idx_plus1 != 3);
                bdof_enable = check_bdof(log2_cb_w, log2_cb_h, ciip_flag, bcw_flag, smvd_flag);

                bdof_enable = bdof_enable && check_bdof_ref(inter_ctx, ref_idx0, ref_idx1);
            }

            if (inter_ctx->dmvr_enabled && mv_info.inter_dir == 0x3) {
                /*FIXME check both mv in bir dir ?*/
                uint8_t bcw_flag = (mv_info.mv0.bcw_idx_plus1 != 0 && mv_info.mv0.bcw_idx_plus1 != 3);
                dmvr_enable = check_bdof(log2_cb_w, log2_cb_h, 0, mmvd_flag, bcw_flag);

                dmvr_enable = dmvr_enable && check_bdof_ref(inter_ctx, ref_idx0, ref_idx1);
            }

            if (!bdof_enable && !dmvr_enable) {
                ctu_dec->rcn_funcs.rcn_mcp_b(ctu_dec, ctu_dec->rcn_ctx.ctu_buff, inter_ctx, part_ctx,
                          mv_info.mv0, mv_info.mv1, x0, y0,
                          log2_cb_w, log2_cb_h, mv_info.inter_dir, ref_idx0, ref_idx1);
            } else if (dmvr_enable) {
                uint8_t log2_w = OVMIN(log2_cb_w, 4);
                uint8_t log2_h = OVMIN(log2_cb_h, 4);
                uint8_t nb_sb_w = (1 << log2_cb_w) >> log2_w;
                uint8_t nb_sb_h = (1 << log2_cb_h) >> log2_h;

                int i, j;

                OVMV mv0 = mv_info.mv0;
                OVMV mv1 = mv_info.mv1;
                for (i = 0; i < nb_sb_h; ++i) {
                    for (j = 0; j < nb_sb_w; ++j) {
                        if (dmvr_enable) {
                            struct TMVPMV *tmvp_mv0 = inter_ctx->tmvp_mv[0].mvs;
                            struct TMVPMV *tmvp_mv1 = inter_ctx->tmvp_mv[1].mvs;
                            ctu_dec->rcn_funcs.rcn_dmvr_mv_refine(ctu_dec, ctu_dec->rcn_ctx.ctu_buff,
                                               x0 + j * 16, y0 + i * 16,
                                               log2_w, log2_h,
                                               &mv0, &mv1,
                                               ref_idx0, ref_idx1, bdof_enable);

                            struct TMVPMV tmvpmv0 = {.mv.x=mv0.x, .mv.y=mv0.y, .z=inter_ctx->dist_ref_0[ref_idx0]};
                            struct TMVPMV tmvpmv1 = {.mv.x=mv1.x, .mv.y=mv1.y, .z=inter_ctx->dist_ref_1[ref_idx1]};
                            /* FIXME temporary hack to override MVs on 8x8 grid for TMVP */
                            /* FIXME find an alternative in case x0 % 8  is != 0 */
                            tmvp_mv0[((x0 + 7 + j * 16) >> 3) + ((y0 + 7 + i * 16) >> 3) * 16] = tmvpmv0;
                            tmvp_mv1[((x0 + 7 + j * 16) >> 3) + ((y0 + 7 + i * 16) >> 3) * 16] = tmvpmv1;

                            if (log2_w > 3) {
                                tmvp_mv0[((x0 + 7 + j * 16) >> 3) + ((y0 + 7 + i * 16) >> 3) * 16 + 1] = tmvpmv0;
                                tmvp_mv1[((x0 + 7 + j * 16) >> 3) + ((y0 + 7 + i * 16) >> 3) * 16 + 1] = tmvpmv1;
                            }

                            if (log2_h > 3) {
                                tmvp_mv0[((x0 + 7 + j * 16) >> 3) + ((y0 + 7 + i * 16) >> 3) * 16 + 16] = tmvpmv0;
                                tmvp_mv1[((x0 + 7 + j * 16) >> 3) + ((y0 + 7 + i * 16) >> 3) * 16 + 16] = tmvpmv1;

                                if (log2_w > 3) {
                                    tmvp_mv0[((x0 + 7 + j * 16) >> 3) + ((y0 + 7 + i * 16) >> 3) * 16 + 16+ 1] = tmvpmv0;
                                    tmvp_mv1[((x0 + 7 + j * 16) >> 3) + ((y0 + 7 + i * 16) >> 3) * 16 + 16+ 1] = tmvpmv1;
                                }
                            }
                        }
                        mv0 = mv_info.mv0;
                        mv1 = mv_info.mv1;
                    }
                }
            } else if (bdof_enable) {
                int i, j;
                uint8_t log2_w = OVMIN(log2_cb_w, 4);
                uint8_t log2_h = OVMIN(log2_cb_h, 4);
                uint8_t nb_sb_w = (1 << log2_cb_w) >> log2_w;
                uint8_t nb_sb_h = (1 << log2_cb_h) >> log2_h;
                OVMV mv0 = mv_info.mv0;
                OVMV mv1 = mv_info.mv1;
                for (i = 0; i < nb_sb_h; ++i) {
                    for (j = 0; j < nb_sb_w; ++j) {
                        ctu_dec->rcn_funcs.rcn_bdof_mcp_l(ctu_dec, ctu_dec->rcn_ctx.ctu_buff,
                                       x0 + j * 16, y0 + i * 16, log2_w, log2_h,
                                       mv0, mv1, ref_idx0, ref_idx1);
                    }
                }
                ctu_dec->rcn_funcs.rcn_mcp_b_c(ctu_dec, ctu_dec->rcn_ctx.ctu_buff, inter_ctx, part_ctx,
                            mv0, mv1, x0, y0,
                            log2_cb_w, log2_cb_h, mv_info.inter_dir, ref_idx0, ref_idx1);
            }
        }
    } else {
        uint8_t inter_dir = ovcabac_read_ae_inter_dir(cabac_ctx, log2_cb_w, log2_cb_h);

        if (inter_dir == 0x3) {
            uint8_t nb_active_ref0_min1 = inter_ctx->nb_active_ref0 - 1;
            uint8_t nb_active_ref1_min1 = inter_ctx->nb_active_ref1 - 1;

            cu_type = read_bidir_mvp(ctu_dec, x0, y0, log2_cb_w, log2_cb_h,
                                     log2_min_cb_s, nb_active_ref0_min1,
                                     nb_active_ref1_min1);
            goto end;

        } else {
            struct MVPInfoP mvp_info;
            if (ctu_dec->affine_enabled && log2_cb_w > 3 && log2_cb_h > 3) {
                uint8_t y_cb = y0 >> log2_min_cb_s;
                uint8_t x_cb = x0 >> log2_min_cb_s;
                uint8_t cu_type_abv = ctu_dec->part_map.cu_mode_x[x_cb];
                uint8_t cu_type_lft = ctu_dec->part_map.cu_mode_y[y_cb];

                uint8_t lft_affine = cu_type_lft == OV_AFFINE || cu_type_lft == OV_INTER_SKIP_AFFINE;
                uint8_t abv_affine = cu_type_abv == OV_AFFINE || cu_type_abv == OV_INTER_SKIP_AFFINE;

                uint8_t affine_flag = ovcabac_read_ae_cu_affine_flag(cabac_ctx, lft_affine, abv_affine);
                if (affine_flag) {
                    uint8_t six_affine_type = !!(ctu_dec->affine_status & 0x2);
                    uint8_t affine_type = six_affine_type && ovcabac_read_ae_cu_affine_type(cabac_ctx);
                    uint8_t nb_active_ref_min1 = inter_dir & 0x1 ? inter_ctx->nb_active_ref0 - 1
                                                                 : inter_ctx->nb_active_ref1 - 1;

                    uint8_t prec_amvr = MV_PRECISION_QUARTER;
                    struct AffineMVPDataP aff_mvp_data = inter_affine_mvp_data_p(ctu_dec, nb_active_ref_min1, affine_type);

                    if (ctu_dec->affine_amvr_enabled) {
                        int32_t nz_mvd = check_nz_affine_mvd_p(&aff_mvp_data.mvd, affine_type);

                        if (nz_mvd) {
                            prec_amvr = ovcabac_read_ae_affine_amvr_precision(cabac_ctx);
                        }
                    }

                    cu_type = OV_AFFINE;
                    mvp_info.data.cp_mvd_info = aff_mvp_data;
                    mvp_info.prec_amvr   = prec_amvr;
                    mvp_info.affine_type = affine_type;

                    drv_affine_mvp_p(inter_ctx, x0, y0, log2_cb_w, log2_cb_h,
                                     &mvp_info, inter_dir);

                    goto end;
                }
            }

            uint8_t prec_amvr = MV_PRECISION_QUARTER;
            uint8_t nb_active_ref_min1 = inter_dir & 0x1 ? inter_ctx->nb_active_ref0 - 1
                                                         : inter_ctx->nb_active_ref1 - 1;

            struct MVPDataP mvp_data = inter_mvp_data_p(ctu_dec, nb_active_ref_min1);

            if (ctu_dec->amvr_enabled) {
                uint8_t nz_mvd = check_nz_mvd_p(&mvp_data.mvd);
                if (nz_mvd) {
                    prec_amvr = ovcabac_read_ae_amvr_precision(cabac_ctx);
                }
            }

            mv_info = drv_mvp_b(inter_ctx, x0, y0, log2_cb_w, log2_cb_h,
                                mvp_data.mvd, mvp_data.mvd, prec_amvr, mvp_data.mvp_idx, mvp_data.mvp_idx, BCW_DEFAULT,
                                inter_dir, mvp_data.ref_idx, mvp_data.ref_idx, log2_cb_w + log2_cb_h <= 5);

            ctu_dec->rcn_funcs.rcn_mcp_b(ctu_dec, ctu_dec->rcn_ctx.ctu_buff, inter_ctx, part_ctx,
                      mv_info.mv0, mv_info.mv1, x0, y0,
                      log2_cb_w, log2_cb_h, inter_dir, mvp_data.ref_idx, mvp_data.ref_idx);
        }
    }

end:

    /*FIXME this have to be moved to DRV */
    reset_intra_map(ctu_dec, i_info, x0, y0, log2_cb_w, log2_cb_h, log2_min_cb_s);

    return cu_type;
}
