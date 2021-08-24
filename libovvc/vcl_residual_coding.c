/* FIXME All this file needs cleaning + refactoring :
 * Actions :
 *           - Move dequant processing outside of coeff coding funcitons
 *           - Procees by TU size, this will permit to remove ISP functions
 *         since when we know the size we also know the Sub Block size to use
 *         also it will be easier to reset neighbour context tables for contest
 *         offsets derivation
 *           - Find an other whay for sign_data hiding and dep quant specialization
 *           - Branch less neighbour context update
 */

#include <string.h>
#include "ovutils.h"
#include "data_scan_lut.h"
#include "cabac_internal.h"
#include "ctudec.h"
#include "vcl.h"

#define IQUANT_SHIFT 6
#define ADJ_QUANT_SHIFT 7
#define ADJ_DEQUANT_SHIFT ADJ_QUANT_SHIFT + 1

#define VVC_TR_CTX_STRIDE (64+2)
#define VVC_TR_CTX_OFFSET ((VVC_TR_CTX_STRIDE)*2+2)
#define VVC_TR_CTX_SIZE   (VVC_TR_CTX_STRIDE*VVC_TR_CTX_STRIDE)

enum TBSize
{
    TB_1x1     = 0,
    TB_1x2     = 1,
    TB_1x4     = 2,
    TB_1x8     = 3,
    TB_1x16    = 4,
    TB_1x32    = 5,
    TB_1x64    = 6,
    TB_1x128   = 7,
    TB_2x1     = 8,
    TB_2x2     = 9,
    TB_2x4     = 10,
    TB_2x8     = 11,
    TB_2x16    = 12,
    TB_2x32    = 13,
    TB_2x64    = 14,
    TB_2x128   = 15,
    TB_4x1     = 16,
    TB_4x2     = 17,
    TB_4x4     = 18,
    TB_4x8     = 19,
    TB_4x16    = 20,
    TB_4x32    = 21,
    TB_4x64    = 22,
    TB_4x128   = 23,
    TB_8x1     = 24,
    TB_8x2     = 25,
    TB_8x4     = 26,
    TB_8x8     = 27,
    TB_8x16    = 28,
    TB_8x32    = 29,
    TB_8x64    = 30,
    TB_8x128   = 31,
    TB_16x1    = 32,
    TB_16x2    = 33,
    TB_16x4    = 34,
    TB_16x8    = 35,
    TB_16x16   = 36,
    TB_16x32   = 37,
    TB_16x64   = 38,
    TB_16x128  = 39,
    TB_32x1    = 40,
    TB_32x2    = 41,
    TB_32x4    = 42,
    TB_32x8    = 43,
    TB_32x16   = 44,
    TB_32x32   = 45,
    TB_32x64   = 46,
    TB_32x128  = 47,
    TB_64x1    = 48,
    TB_64x2    = 49,
    TB_64x4    = 50,
    TB_64x8    = 51,
    TB_64x16   = 52,
    TB_64x32   = 53,
    TB_64x64   = 54,
    TB_64x128  = 55,
    TB_128x1   = 56,
    TB_128x2   = 57,
    TB_128x4   = 58,
    TB_128x8   = 59,
    TB_128x16  = 60,
    TB_128x32  = 61,
    TB_128x64  = 62,
    TB_128x128 = 63,
};

typedef struct VVCDepQuantCtx{
    const uint16_t state_trans_table;
    const uint8_t state_offset;
}VVCDepQuantCtx;

typedef struct VVCCoeffCodingCtx{
    uint8_t *sum_sig_nbs;
    uint8_t *sum_abs_lvl;
    uint8_t *sum_abs_lvl2;
    int16_t nb_remaining_bins;
    uint8_t enable_sdh;
}VVCCoeffCodingCtx;

typedef struct TSCoeffCodingCtx{
    uint8_t *nb_sig_ngh;
    uint8_t *sign_map;
    uint16_t *abs_coeffs;
    int16_t nb_remaining_bins;
}TSCoeffCodingCtx;

typedef struct VVCResidualStates{
    const uint16_t sig_flg_ctx_offset;
    const uint16_t abs_gt1_ctx_offset;
    const uint16_t par_lvl_ctx_offset;
    const uint16_t abs_gt2_ctx_offset;
}VVCSBStates;

typedef struct VVCSBScanContext{
   const uint64_t scan_map;
   const uint64_t scan_idx_map;
   const uint8_t log2_sb_w;
   const uint8_t log2_sb_h;
}VVCSBScanContext;

static const VVCDepQuantCtx luma_dep_quant_ctx = {
    0x7D28,
    12
};

static const VVCDepQuantCtx chroma_dep_quant_ctx = {
    0x7D28,
    8
};

#if 0
static const uint64_t inv_diag_map_4x4 = 0x041852C963DA7EBF;
static const uint64_t inv_diag_map_2x8 = 0x021436587A9CBEDF;
static const uint64_t inv_diag_map_8x2 = 0x08192A3B4C5D6E7F;
#endif

static const uint64_t diag_map_4X4 = 0xFBE7AD369C258140;
static const uint64_t diag_map_2x8 = 0xFDEBC9A785634120;
static const uint64_t diag_map_8x2 = 0xF7E6D5C4B3A29180;
static const uint64_t diag_map_2x2 = 0x0000000000003120;


static const VVCSBScanContext inv_diag_4x4_scan =
{
     0x041852C963DA7EBF,
     0x0259148C37BE6ADF,
     2,
     2,
};

static const VVCSBScanContext inv_diag_2x8_scan =
{
     0x021436587A9CBEDF,
     0x021436587A9CBEDF,
     1,
     3,
};

static const VVCSBScanContext inv_diag_8x2_scan =
{
     0x08192A3B4C5D6E7F,
     0x02468ACE13579BDF,
     3,
     1,
};

static const VVCSBScanContext inv_diag_2x4_scan =
{
     0x0213465700000000,
     0x0213465700000000,
     1,
     3,
};

static const VVCSBScanContext inv_diag_4x2_scan =
{
     0x08192A3B4C5D6E7F,
     0x02468ACE13579BDF,
     3,
     1,
};

static const VVCSBScanContext inv_diag_2x2_scan =
{
     0x0213,
     0x0213,
     1,
     1,
};

static const VVCSBScanContext inv_diag_1x16_scan =
{
     0x0123456789ABCDEF,
     0x0123456789ABCDEF,
     4,
     0,
};

/* Transform skip scan_ctx */
static const VVCSBScanContext diag_4x4_scan =
{
     0xFBE7AD369C258140,
     0x0,
     2,
     2,
};

static const VVCSBScanContext diag_2x8_scan =
{
     0xFDEBC9A785634120,
     0x0,
     1,
     3,
};

static const VVCSBScanContext diag_8x2_scan =
{
     0xF7E6D5C4B3A29180,
     0x0,
     3,
     1,
};

static const VVCSBScanContext diag_2x2_scan =
{
     0x3120,
     0x0,
     1,
     1,
};


static const uint64_t parity_flag_offset_map[3] =
{
    0xFAAAAA5555555555,
    0x5555555555555550,
    0x5550000000000000
};

static const uint64_t sig_flag_offset_map[3] =
{
    0x8884444444444000,
    0x4000000000000000,
    0x0000000000000000
};

static const VVCSBStates luma_ctx_offsets =
{
    SIG_FLAG_CTX_OFFSET,
    GT0_FLAG_CTX_OFFSET,
    PAR_FLAG_CTX_OFFSET,
    GT1_FLAG_CTX_OFFSET,
};

static const VVCSBStates chroma_ctx_offsets =
{
    SIG_FLAG_C_CTX_OFFSET,
    GT0_FLAG_C_CTX_OFFSET,
    PAR_FLAG_C_CTX_OFFSET,
    GT1_FLAG_C_CTX_OFFSET,
};

static uint8_t
ovcabac_read_ae_significant_sb_flag(OVCABACCtx *const cabac_ctx,
                                    uint8_t got_significant_neighbour)
{
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    return ovcabac_ae_read(cabac_ctx, &cabac_state[SIG_COEFF_GROUP_CTX_OFFSET + got_significant_neighbour]);
}

static uint8_t
ovcabac_read_ae_significant_sb_flag_chroma(OVCABACCtx *const cabac_ctx,
                                           uint8_t got_significant_neighbour)
{
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    return ovcabac_ae_read(cabac_ctx, &cabac_state[SIG_COEFF_GROUP_C_CTX_OFFSET + got_significant_neighbour]);
}

static uint8_t
ovcabac_read_ae_significant_ts_sb_flag(OVCABACCtx *const cabac_ctx,
                                       uint8_t got_significant_neighbour)
{
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    return ovcabac_ae_read(cabac_ctx, &cabac_state[TS_SIG_COEFF_GROUP_CTX_OFFSET + got_significant_neighbour]);
}

static void inline
set_implicit_coeff_ngbh(const VVCCoeffCodingCtx *const coef_nbh_ctx,
                        int tr_ctx_pos, int value)
{
    //int value2=value;//OVMIN(value,(4+(value&1)));
    coef_nbh_ctx->sum_abs_lvl[tr_ctx_pos                     - 1] += value;
    coef_nbh_ctx->sum_abs_lvl[tr_ctx_pos                     - 2] += value;
    coef_nbh_ctx->sum_abs_lvl[tr_ctx_pos - VVC_TR_CTX_STRIDE    ] += value;
    coef_nbh_ctx->sum_abs_lvl[tr_ctx_pos - VVC_TR_CTX_STRIDE - 1] += value;
    coef_nbh_ctx->sum_abs_lvl[tr_ctx_pos - VVC_TR_CTX_STRIDE * 2] += value;

    coef_nbh_ctx->sum_abs_lvl2[tr_ctx_pos                     - 1] += value;
    coef_nbh_ctx->sum_abs_lvl2[tr_ctx_pos                     - 2] += value;
    coef_nbh_ctx->sum_abs_lvl2[tr_ctx_pos - VVC_TR_CTX_STRIDE    ] += value;
    coef_nbh_ctx->sum_abs_lvl2[tr_ctx_pos - VVC_TR_CTX_STRIDE - 1] += value;
    coef_nbh_ctx->sum_abs_lvl2[tr_ctx_pos - VVC_TR_CTX_STRIDE * 2] += value;

    coef_nbh_ctx->sum_sig_nbs[tr_ctx_pos                     - 1] += value - 1;
    coef_nbh_ctx->sum_sig_nbs[tr_ctx_pos                     - 2] += value - 1;
    coef_nbh_ctx->sum_sig_nbs[tr_ctx_pos - VVC_TR_CTX_STRIDE    ] += value - 1;
    coef_nbh_ctx->sum_sig_nbs[tr_ctx_pos - VVC_TR_CTX_STRIDE - 1] += value - 1;
    coef_nbh_ctx->sum_sig_nbs[tr_ctx_pos - VVC_TR_CTX_STRIDE * 2] += value - 1;
}

#define updt_sat(x,val,sat) \
(x) = OVMIN((sat),(x)+(val));

static void inline
update_coeff_nbgh_bypassed(const VVCCoeffCodingCtx *const coef_nbh_ctx,
                           int tr_ctx_pos, int value)
{
    /* The value needs to be saturated to 5 or 4 for sum_abs_lvl1 */
    /*FIXME when in bypass we don't use sum_abs_lvl1 */

    updt_sat(coef_nbh_ctx->sum_abs_lvl2[tr_ctx_pos                     - 1], value, 51);
    updt_sat(coef_nbh_ctx->sum_abs_lvl2[tr_ctx_pos                     - 2], value, 51);
    updt_sat(coef_nbh_ctx->sum_abs_lvl2[tr_ctx_pos - VVC_TR_CTX_STRIDE    ], value, 51);
    updt_sat(coef_nbh_ctx->sum_abs_lvl2[tr_ctx_pos - VVC_TR_CTX_STRIDE - 1], value, 51);
    updt_sat(coef_nbh_ctx->sum_abs_lvl2[tr_ctx_pos - VVC_TR_CTX_STRIDE * 2], value, 51);
}

static void inline
update_coeff_nbgh_first_pass(const VVCCoeffCodingCtx *const coef_nbh_ctx,
                             int tr_ctx_pos, int value)
{
    coef_nbh_ctx->sum_abs_lvl[tr_ctx_pos                     - 1] += value;
    coef_nbh_ctx->sum_abs_lvl[tr_ctx_pos                     - 2] += value;
    coef_nbh_ctx->sum_abs_lvl[tr_ctx_pos - VVC_TR_CTX_STRIDE    ] += value;
    coef_nbh_ctx->sum_abs_lvl[tr_ctx_pos - VVC_TR_CTX_STRIDE - 1] += value;
    coef_nbh_ctx->sum_abs_lvl[tr_ctx_pos - VVC_TR_CTX_STRIDE * 2] += value;

    /* There is no need to sature value since it won't exceed 25*/
    coef_nbh_ctx->sum_abs_lvl2[tr_ctx_pos                     - 1] += value;
    coef_nbh_ctx->sum_abs_lvl2[tr_ctx_pos                     - 2] += value;
    coef_nbh_ctx->sum_abs_lvl2[tr_ctx_pos - VVC_TR_CTX_STRIDE    ] += value;
    coef_nbh_ctx->sum_abs_lvl2[tr_ctx_pos - VVC_TR_CTX_STRIDE - 1] += value;
    coef_nbh_ctx->sum_abs_lvl2[tr_ctx_pos - VVC_TR_CTX_STRIDE * 2] += value;

    coef_nbh_ctx->sum_sig_nbs[tr_ctx_pos                     - 1] += value - 1;
    coef_nbh_ctx->sum_sig_nbs[tr_ctx_pos                     - 2] += value - 1;
    coef_nbh_ctx->sum_sig_nbs[tr_ctx_pos - VVC_TR_CTX_STRIDE    ] += value - 1;
    coef_nbh_ctx->sum_sig_nbs[tr_ctx_pos - VVC_TR_CTX_STRIDE - 1] += value - 1;
    coef_nbh_ctx->sum_sig_nbs[tr_ctx_pos - VVC_TR_CTX_STRIDE * 2] += value - 1;
}

static void inline
update_coeff_nbgh_other_pass(const VVCCoeffCodingCtx *const coef_nbh_ctx,
                             int tr_ctx_pos, int value)
{
    updt_sat(coef_nbh_ctx->sum_abs_lvl2[tr_ctx_pos                     - 1], value, 51);
    updt_sat(coef_nbh_ctx->sum_abs_lvl2[tr_ctx_pos                     - 2], value, 51);
    updt_sat(coef_nbh_ctx->sum_abs_lvl2[tr_ctx_pos - VVC_TR_CTX_STRIDE    ], value, 51);
    updt_sat(coef_nbh_ctx->sum_abs_lvl2[tr_ctx_pos - VVC_TR_CTX_STRIDE - 1], value, 51);
    updt_sat(coef_nbh_ctx->sum_abs_lvl2[tr_ctx_pos - VVC_TR_CTX_STRIDE * 2], value, 51);
}

static int inline
decode_truncated_rice(OVCABACCtx *const cabac_ctx, unsigned int rice_param){

    unsigned int prefix = 0;
    unsigned int length = rice_param;
    int offset;
    int value = 0;
    int cut_off = 5;
    unsigned int code = 0;

    //FIXME check whether an other writing is more efficient
    do{
        prefix++;
        code = ovcabac_bypass_read(cabac_ctx);
    } while( code && prefix < 17 );
    prefix -= 1 - code;


    if(prefix < cut_off){
        offset  = prefix << rice_param;
    } else {
        offset  = ((( 1 << ( prefix - cut_off )) + cut_off - 1) << rice_param );
        length += ( prefix == 17 ? 15 - rice_param : prefix - 5 );
    }

    while(length--){
        value <<= 1;
        value |= ovcabac_bypass_read(cabac_ctx);
    }

    value += offset;
    value <<= 1;

    return value;
}

static const uint8_t rice_param_tab[32] =
{
  0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3
};

static inline void
decode_pass2_core(OVCABACCtx *const cabac_ctx,
                  int16_t  *const coeffs, int nb_pass2,
                  uint8_t  *const next_pass_idx_map,
                  const VVCSBScanContext *const scan_ctx,
                  const VVCCoeffCodingCtx *const c_coding_ctx)
{
    int scan_pos;
    #if 0
    const uint64_t inv_diag_map = scan_ctx->scan_map;
    #endif
    const uint8_t log2_sb_w     = scan_ctx->log2_sb_w;
    const uint8_t x_mask = (1 << log2_sb_w) - 1;

    for (scan_pos = 0; scan_pos < nb_pass2; ++scan_pos){
        int idx = next_pass_idx_map[scan_pos];
        int rem_abs_lvl;
        int sum_abs ;
        int rice_param ;
        int tr_ctx_pos = (idx & x_mask) + (idx >> log2_sb_w) * VVC_TR_CTX_STRIDE;

        sum_abs = OVMAX(OVMIN((int)c_coding_ctx->sum_abs_lvl2[tr_ctx_pos] - 20, 31), 0);

        rice_param = rice_param_tab [sum_abs];

        rem_abs_lvl = decode_truncated_rice(cabac_ctx, rice_param);

        if (rem_abs_lvl){
            update_coeff_nbgh_other_pass(c_coding_ctx, tr_ctx_pos, rem_abs_lvl);
        }

        coeffs[idx] += rem_abs_lvl;
    }
}

static inline void
decode_bypassed_coeff_core(OVCABACCtx *const cabac_ctx,
                      int16_t *const coeffs, int last_scan_pos,
                      uint8_t *const significant_map,
                      int *const nb_significant_coeff,
                      const VVCSBScanContext *const scan_ctx,
                      const VVCCoeffCodingCtx *const c_coding_ctx,
                      const VVCDepQuantCtx *const dep_quant,
                      int *const state,
                      uint32_t *const state_map)
{
    int scan_pos;
    int rice_param;

    const uint64_t inv_diag_map = scan_ctx->scan_map;
    const uint8_t log2_sb_w     = scan_ctx->log2_sb_w;
    const uint8_t x_mask = (1 << log2_sb_w) - 1;

    const uint16_t state_trans_tab = dep_quant->state_trans_table;
    int max_start = (1 << (scan_ctx->log2_sb_w + scan_ctx->log2_sb_h)) - 1;
    uint8_t pos_shift = ((max_start - (last_scan_pos & 0xF)) << 2);
    uint64_t scan_map = inv_diag_map >> pos_shift;

    for(scan_pos = last_scan_pos; scan_pos >= 0; --scan_pos){
        int idx = scan_map & 0xF;
        int value;

        int tr_ctx_pos = (idx & x_mask) + (idx >> log2_sb_w) * VVC_TR_CTX_STRIDE;

        int sum_abs = OVMIN(31, c_coding_ctx->sum_abs_lvl2[tr_ctx_pos]);

        rice_param = rice_param_tab[sum_abs];

        value = decode_truncated_rice(cabac_ctx, rice_param) >> 1;

        /*FIXME understand value saturation*/
        value = (value == ((*state < 2 ? 1: 2) << rice_param))
              ?  0 : (value  < ((*state < 2 ? 1: 2) << rice_param))
              ?  value + 1 : value;

        if(value){
            update_coeff_nbgh_bypassed(c_coding_ctx, tr_ctx_pos, value);

            coeffs[idx] = value;

            *state_map <<= 1;
            *state_map |= (*state >> 1);

            significant_map[(*nb_significant_coeff)++] = idx;
        }
        // update state transition context
        *state = ( state_trans_tab >> (( *state << 2) + ((value & 1) << 1)) ) & 3;
        scan_map >>= 4;
    }
}

static void inline
decode_signs(OVCABACCtx *const cabac_ctx, int16_t *const coeffs,
             uint32_t state_map, const uint8_t *const sig_c_idx_map,
             int nb_sig_c)
{
    uint32_t nb_signs = nb_sig_c;
    uint32_t signs_map = 0;

    while (nb_signs--){
        signs_map <<= 1;
        signs_map  |= ovcabac_bypass_read(cabac_ctx);
    }

    signs_map <<= 32 - nb_sig_c;
    state_map <<= 32 - nb_sig_c;

    for (unsigned  k = 0; k < nb_sig_c; k++ ){
        int idx = sig_c_idx_map[k];
        int add  = !!(state_map & (1u << 31));
        int sign = !!(signs_map & (1u << 31));
        int abs_coeff = (coeffs[idx] << 1) - add;
        coeffs[idx] = sign ? -abs_coeff : abs_coeff;
        signs_map <<= 1;
        state_map <<= 1;
    }
}

static inline int
residual_coding_first_subblock_4x4(OVCABACCtx *const cabac_ctx,
                                  int16_t  *const coeffs,
                                  int *const state, int start_pos,
                                  uint64_t par_flag_offset_map,
                                  uint64_t sig_flag_offset_map,
                                  VVCCoeffCodingCtx *const c_coding_ctx,
                                  const VVCSBStates *const ctx_offsets,
                                  const VVCDepQuantCtx *const dep_quant_ctx,
                                  const VVCSBScanContext *const scan_ctx)
{
    const uint64_t inv_diag_map = scan_ctx->scan_map;
    const uint8_t log2_sb_w     = scan_ctx->log2_sb_w;
    const uint8_t x_mask    = (1 << log2_sb_w) - 1;

    uint64_t *const sig_flg_ctx = &cabac_ctx->ctx_table[ctx_offsets->sig_flg_ctx_offset];
    uint64_t *const abs_gt1_ctx = &cabac_ctx->ctx_table[ctx_offsets->abs_gt1_ctx_offset];
    uint64_t *const par_lvl_ctx = &cabac_ctx->ctx_table[ctx_offsets->par_lvl_ctx_offset];
    uint64_t *const abs_gt2_ctx = &cabac_ctx->ctx_table[ctx_offsets->abs_gt2_ctx_offset];

    const uint16_t state_trans_tab = dep_quant_ctx->state_trans_table;
    const uint8_t state_offset     = dep_quant_ctx->state_offset;

    uint8_t sig_idx_map[16];
    uint8_t gt2_idx_map[16];
    int nb_sig_c = 0;
    int nb_pass2 = 0;

    uint32_t dep_quant_map = 0;

    uint8_t par_lvl_flag, abs_gt1_flag, abs_gt2_flag;
    int32_t coeff_val;
    uint16_t tr_ctx_pos;

    int nb_rem_bins = c_coding_ctx->nb_remaining_bins;
    int prev_state  = *state;

    // Implicit first coeff
    int scan_pos = start_pos;
    int max_start = (1 << (scan_ctx->log2_sb_w + scan_ctx->log2_sb_h)) - 1;
    uint8_t pos_shift = ((max_start - (scan_pos & 0xF)) << 2);
    uint64_t scan_map = inv_diag_map        >> pos_shift;
    uint64_t par_map  = par_flag_offset_map >> pos_shift;
    uint64_t sig_map  = sig_flag_offset_map >> pos_shift;

    uint8_t idx = scan_map & 0xF;

    abs_gt1_flag = ovcabac_ae_read(cabac_ctx, abs_gt1_ctx);

    --nb_rem_bins;
    coeff_val = 1 + abs_gt1_flag;

    if(abs_gt1_flag){
        par_lvl_flag = ovcabac_ae_read(cabac_ctx, par_lvl_ctx);
        abs_gt2_flag = ovcabac_ae_read(cabac_ctx, abs_gt2_ctx);
        coeff_val += par_lvl_flag;
        nb_rem_bins -= 2;
        if(abs_gt2_flag){
           coeff_val += 2; //(rem_abs_gt2_flag << 1)
           gt2_idx_map[nb_pass2++] = idx;
        }
    }

    coeffs[idx] = coeff_val;
    sig_idx_map[nb_sig_c++] = idx;

    tr_ctx_pos = (idx & x_mask) + (idx >> log2_sb_w) * VVC_TR_CTX_STRIDE;

    set_implicit_coeff_ngbh(c_coding_ctx, tr_ctx_pos, coeff_val);

    dep_quant_map |= (prev_state >> 1);
    prev_state = (state_trans_tab >> ((prev_state << 2) + ((coeff_val & 1) << 1))) & 3;

    --scan_pos;
    scan_map >>= 4;
    par_map  >>= 4;
    sig_map  >>= 4;

    // First pass
    for( ; scan_pos >= 0 && nb_rem_bins >= 4; --scan_pos ){

        uint8_t ctx_offset;
        uint8_t sig_coeff_flag;

        idx = scan_map & 0xF;
        tr_ctx_pos = (idx & x_mask) + (idx >> log2_sb_w) * VVC_TR_CTX_STRIDE;

        /*FIXME we could state ctx switch by same offset for chroma and luma
        */
        ctx_offset  = state_offset * OVMAX(0, prev_state - 1);
        ctx_offset += OVMIN(((c_coding_ctx->sum_abs_lvl[tr_ctx_pos] + 1) >> 1), 3);
        ctx_offset += sig_map & 0xF;

        sig_coeff_flag = ovcabac_ae_read(cabac_ctx, sig_flg_ctx + ctx_offset);

        coeff_val = sig_coeff_flag;
        --nb_rem_bins;

        if (sig_coeff_flag){

            ctx_offset  = 1;
            ctx_offset += OVMIN(c_coding_ctx->sum_sig_nbs[tr_ctx_pos], 4);
            ctx_offset += (par_map & 0xF);

            abs_gt1_flag = ovcabac_ae_read(cabac_ctx, abs_gt1_ctx + ctx_offset);

            if (abs_gt1_flag){
                par_lvl_flag = ovcabac_ae_read(cabac_ctx, par_lvl_ctx + ctx_offset);
                abs_gt2_flag = ovcabac_ae_read(cabac_ctx, abs_gt2_ctx + ctx_offset);
                coeff_val = 2 + par_lvl_flag;
                nb_rem_bins -= 2;
                if (abs_gt2_flag){
                    coeff_val += 2; /* (abs_gt2_flag << 1) */
                    gt2_idx_map[nb_pass2++] = idx;
                }
            }

            dep_quant_map <<= 1;
            dep_quant_map |= (prev_state >> 1);

            --nb_rem_bins;
            sig_idx_map[nb_sig_c++] = idx;

            coeffs[idx] = coeff_val;
            update_coeff_nbgh_first_pass(c_coding_ctx, tr_ctx_pos, coeff_val);
        }
        prev_state = (state_trans_tab >> (( prev_state << 2) + ((coeff_val & 1) << 1))) & 3;
        scan_map >>= 4;
        par_map  >>= 4;
        sig_map  >>= 4;
    }

    if (nb_pass2){
        decode_pass2_core(cabac_ctx, coeffs, nb_pass2, gt2_idx_map,
                          scan_ctx, c_coding_ctx);
    }

    if (scan_pos >= 0){
        decode_bypassed_coeff_core(cabac_ctx, coeffs, scan_pos, sig_idx_map,
                                   &nb_sig_c, scan_ctx, c_coding_ctx,
                                   &luma_dep_quant_ctx, &prev_state, &dep_quant_map);
    }

    decode_signs(cabac_ctx, coeffs, dep_quant_map, sig_idx_map, nb_sig_c);

    c_coding_ctx->nb_remaining_bins = nb_rem_bins;
    *state = prev_state;

    /*FIXME we could return state instead of nb_sig and avoid using pointers*/
    return nb_sig_c;
}

static inline int
residual_coding_subblock_4x4(OVCABACCtx *const cabac_ctx,
                             int16_t  *const coeffs,
                             int *const state, int start_pos,
                             uint64_t par_flag_offset_map,
                             uint64_t sig_flag_offset_map,
                             VVCCoeffCodingCtx *const c_coding_ctx,
                             const VVCSBStates *const ctx_offsets,
                             const VVCDepQuantCtx *const dep_quant_ctx,
                             const VVCSBScanContext *const scan_ctx)
{
    const uint64_t inv_diag_map = scan_ctx->scan_map;
    const uint8_t log2_sb_w     = scan_ctx->log2_sb_w;
    const uint8_t x_mask    = (1 << log2_sb_w) - 1;

    uint64_t *const sig_flg_ctx = &cabac_ctx->ctx_table[ctx_offsets->sig_flg_ctx_offset];
    uint64_t *const abs_gt1_ctx = &cabac_ctx->ctx_table[ctx_offsets->abs_gt1_ctx_offset];
    uint64_t *const par_lvl_ctx = &cabac_ctx->ctx_table[ctx_offsets->par_lvl_ctx_offset];
    uint64_t *const abs_gt2_ctx = &cabac_ctx->ctx_table[ctx_offsets->abs_gt2_ctx_offset];

    const uint16_t state_trans_tab = dep_quant_ctx->state_trans_table;
    const uint8_t state_offset     = dep_quant_ctx->state_offset;

    uint8_t sig_idx_map[16];
    uint8_t gt2_idx_map[16];
    int nb_sig_c = 0;
    int nb_pass2 = 0;

    uint32_t dep_quant_map = 0;

    uint8_t par_lvl_flag, abs_gt1_flag, abs_gt2_flag;
    int32_t coeff_val;
    uint16_t tr_ctx_pos;

    int nb_rem_bins = c_coding_ctx->nb_remaining_bins;
    int prev_state  = *state;

    int scan_pos = 15;
    uint64_t scan_map = inv_diag_map;
    uint64_t par_map  = par_flag_offset_map;
    uint64_t sig_map  = sig_flag_offset_map;

    uint8_t idx;

    // First pass
    for ( ; scan_pos > 0 && nb_rem_bins >= 4; --scan_pos ){

        uint8_t ctx_offset;
        uint8_t sig_coeff_flag;

        idx = scan_map & 0xF;
        tr_ctx_pos = (idx & x_mask) + (idx >> log2_sb_w) * VVC_TR_CTX_STRIDE;

        /*FIXME we could state ctx switch by same offset for chroma and luma
        */
        ctx_offset  = state_offset * OVMAX(0, prev_state - 1);
        ctx_offset += OVMIN(((c_coding_ctx->sum_abs_lvl[tr_ctx_pos] + 1) >> 1), 3);
        ctx_offset += sig_map & 0xF;

        sig_coeff_flag = ovcabac_ae_read(cabac_ctx, sig_flg_ctx + ctx_offset);

        coeff_val = sig_coeff_flag;
        --nb_rem_bins;

        if (sig_coeff_flag){

            ctx_offset  = 1;
            ctx_offset += OVMIN(c_coding_ctx->sum_sig_nbs[tr_ctx_pos], 4);
            ctx_offset += (par_map & 0xF);

            abs_gt1_flag = ovcabac_ae_read(cabac_ctx, abs_gt1_ctx + ctx_offset);

            if (abs_gt1_flag){
                par_lvl_flag = ovcabac_ae_read(cabac_ctx, par_lvl_ctx + ctx_offset);
                abs_gt2_flag = ovcabac_ae_read(cabac_ctx, abs_gt2_ctx + ctx_offset);
                coeff_val = 2 + par_lvl_flag;
                nb_rem_bins -= 2;
                if (abs_gt2_flag){
                    coeff_val += 2; /* (abs_gt2_flag << 1) */
                    gt2_idx_map[nb_pass2++] = idx;
                }
            }

            dep_quant_map <<= 1;
            dep_quant_map |= (prev_state >> 1);

            --nb_rem_bins;
            sig_idx_map[nb_sig_c++] = idx;

            coeffs[idx] = coeff_val;
            update_coeff_nbgh_first_pass(c_coding_ctx, tr_ctx_pos, coeff_val);
        }
        prev_state = (state_trans_tab >> (( prev_state << 2) + ((coeff_val & 1) << 1))) & 3;
        scan_map >>= 4;
        par_map  >>= 4;
        sig_map  >>= 4;
    }

    if (scan_pos == 0 && nb_rem_bins >= 4){

        uint8_t ctx_offset;
        uint8_t sig_coeff_flag;
        sig_coeff_flag= 1;
        idx = scan_map & 0xF;
        tr_ctx_pos = (idx & x_mask) + (idx >> log2_sb_w) * VVC_TR_CTX_STRIDE;

        //decrease scan_pos so we know last sig_coeff was read in first pass or not
        --scan_pos;

        if (nb_sig_c){
            ctx_offset  = state_offset * OVMAX(0, prev_state - 1);
            ctx_offset += OVMIN(((c_coding_ctx->sum_abs_lvl[tr_ctx_pos] + 1) >> 1), 3);
            ctx_offset += sig_map & 0xF;

            sig_coeff_flag = ovcabac_ae_read(cabac_ctx, sig_flg_ctx + ctx_offset);

            --nb_rem_bins;
        }

        coeff_val = sig_coeff_flag;

        if (sig_coeff_flag){

            ctx_offset  = 1;
            ctx_offset += OVMIN(c_coding_ctx->sum_sig_nbs[tr_ctx_pos], 4);
            ctx_offset += (par_map & 0xF);

            abs_gt1_flag = ovcabac_ae_read(cabac_ctx, abs_gt1_ctx + ctx_offset);

            if (abs_gt1_flag){
                par_lvl_flag = ovcabac_ae_read(cabac_ctx, par_lvl_ctx + ctx_offset);
                abs_gt2_flag = ovcabac_ae_read(cabac_ctx, abs_gt2_ctx + ctx_offset);
                coeff_val = 2 + par_lvl_flag;
                nb_rem_bins -= 2;
                if (abs_gt2_flag){
                    coeff_val += 2; /* (abs_gt2_flag << 1) */
                    gt2_idx_map[nb_pass2++] = idx;
                }
            }

            dep_quant_map <<= 1;
            dep_quant_map |= (prev_state >> 1);

            --nb_rem_bins;
            sig_idx_map[nb_sig_c++] = idx;

            coeffs[idx] = coeff_val;
            update_coeff_nbgh_first_pass(c_coding_ctx, tr_ctx_pos, coeff_val);
        }
        prev_state = (state_trans_tab >> (( prev_state << 2) + ((coeff_val & 1) << 1))) & 3;
        scan_map >>= 4;
        par_map  >>= 4;
        sig_map  >>= 4;
    }

    if (nb_pass2){
        decode_pass2_core(cabac_ctx, coeffs, nb_pass2, gt2_idx_map,
                          scan_ctx, c_coding_ctx);
    }

    if (scan_pos >= 0){
        decode_bypassed_coeff_core(cabac_ctx, coeffs, scan_pos, sig_idx_map,
                                   &nb_sig_c, scan_ctx, c_coding_ctx,
                                   &luma_dep_quant_ctx, &prev_state, &dep_quant_map);
    }

    decode_signs(cabac_ctx, coeffs, dep_quant_map, sig_idx_map, nb_sig_c);

    c_coding_ctx->nb_remaining_bins = nb_rem_bins;
    *state = prev_state;

    /*FIXME we could return state instead of nb_sig and avoid using pointers*/
    return nb_sig_c;
}

static inline int
residual_coding_subblock_dc(OVCABACCtx *const cabac_ctx,
                             int16_t  *const coeffs,
                             int *const state, int start_pos,
                             uint64_t par_flag_offset_map,
                             uint64_t sig_flag_offset_map,
                             VVCCoeffCodingCtx *const c_coding_ctx,
                             const VVCSBStates *const ctx_offsets,
                             const VVCDepQuantCtx *const dep_quant_ctx,
                             const VVCSBScanContext *const scan_ctx)
{
    const uint64_t inv_diag_map = scan_ctx->scan_map;
    const uint8_t log2_sb_w     = scan_ctx->log2_sb_w;
    const uint8_t x_mask    = (1 << log2_sb_w) - 1;

    uint64_t *const sig_flg_ctx = &cabac_ctx->ctx_table[ctx_offsets->sig_flg_ctx_offset];
    uint64_t *const abs_gt1_ctx = &cabac_ctx->ctx_table[ctx_offsets->abs_gt1_ctx_offset];
    uint64_t *const par_lvl_ctx = &cabac_ctx->ctx_table[ctx_offsets->par_lvl_ctx_offset];
    uint64_t *const abs_gt2_ctx = &cabac_ctx->ctx_table[ctx_offsets->abs_gt2_ctx_offset];

    const uint16_t state_trans_tab = dep_quant_ctx->state_trans_table;
    const uint8_t state_offset     = dep_quant_ctx->state_offset;

    uint8_t sig_idx_map[16];
    uint8_t gt2_idx_map[16];
    int nb_sig_c = 0;
    int nb_pass2 = 0;

    uint32_t dep_quant_map = 0;

    uint8_t par_lvl_flag, abs_gt1_flag, abs_gt2_flag;
    int32_t coeff_val;
    uint16_t tr_ctx_pos;

    int nb_rem_bins = c_coding_ctx->nb_remaining_bins;
    int prev_state  = *state;

    int scan_pos = start_pos;
    uint64_t scan_map = inv_diag_map;
    uint64_t par_map  = par_flag_offset_map;
    uint64_t sig_map  = sig_flag_offset_map;

    uint8_t idx;

    // First pass
    for ( ; scan_pos > 0 && nb_rem_bins >= 4; --scan_pos ){

        uint8_t ctx_offset;
        uint8_t sig_coeff_flag;

        idx = scan_map & 0xF;
        tr_ctx_pos = (idx & x_mask) + (idx >> log2_sb_w) * VVC_TR_CTX_STRIDE;

        /*FIXME we could state ctx switch by same offset for chroma and luma
        */
        ctx_offset  = state_offset * OVMAX(0, prev_state - 1);
        ctx_offset += OVMIN(((c_coding_ctx->sum_abs_lvl[tr_ctx_pos] + 1) >> 1), 3);
        ctx_offset += sig_map & 0xF;

        sig_coeff_flag = ovcabac_ae_read(cabac_ctx, sig_flg_ctx + ctx_offset);

        coeff_val = sig_coeff_flag;
        --nb_rem_bins;

        if (sig_coeff_flag){

            ctx_offset  = 1;
            ctx_offset += OVMIN(c_coding_ctx->sum_sig_nbs[tr_ctx_pos], 4);
            ctx_offset += (par_map & 0xF);

            abs_gt1_flag = ovcabac_ae_read(cabac_ctx, abs_gt1_ctx + ctx_offset);

            if (abs_gt1_flag){
                par_lvl_flag = ovcabac_ae_read(cabac_ctx, par_lvl_ctx + ctx_offset);
                abs_gt2_flag = ovcabac_ae_read(cabac_ctx, abs_gt2_ctx + ctx_offset);
                coeff_val = 2 + par_lvl_flag;
                nb_rem_bins -= 2;
                if (abs_gt2_flag){
                    coeff_val += 2; /* (abs_gt2_flag << 1) */
                    gt2_idx_map[nb_pass2++] = idx;
                }
            }

            dep_quant_map <<= 1;
            dep_quant_map |= (prev_state >> 1);

            --nb_rem_bins;
            sig_idx_map[nb_sig_c++] = idx;

            coeffs[idx] = coeff_val;
            update_coeff_nbgh_first_pass(c_coding_ctx, tr_ctx_pos, coeff_val);
        }
        prev_state = (state_trans_tab >> (( prev_state << 2) + ((coeff_val & 1) << 1))) & 3;
        scan_map >>= 4;
        par_map  >>= 4;
        sig_map  >>= 4;
    }

    if (scan_pos == 0 && nb_rem_bins >= 4){

        uint8_t ctx_offset;
        uint8_t sig_coeff_flag;
        sig_coeff_flag= 1;
        idx = scan_map & 0xF;
        tr_ctx_pos = (idx & x_mask) + (idx >> log2_sb_w) * VVC_TR_CTX_STRIDE;

        //decrease scan_pos so we know last sig_coeff was read in first pass or not
        --scan_pos;

        ctx_offset  = state_offset * OVMAX(0, prev_state - 1);
        ctx_offset += OVMIN(((c_coding_ctx->sum_abs_lvl[tr_ctx_pos] + 1) >> 1), 3);
        ctx_offset += sig_map & 0xF;

        sig_coeff_flag = ovcabac_ae_read(cabac_ctx, sig_flg_ctx + ctx_offset);

        --nb_rem_bins;

        coeff_val = sig_coeff_flag;

        if (sig_coeff_flag){

            ctx_offset  = 1;
            ctx_offset += OVMIN(c_coding_ctx->sum_sig_nbs[tr_ctx_pos], 4);
            ctx_offset += (par_map & 0xF);

            abs_gt1_flag = ovcabac_ae_read(cabac_ctx, abs_gt1_ctx + ctx_offset);

            if (abs_gt1_flag){
                par_lvl_flag = ovcabac_ae_read(cabac_ctx, par_lvl_ctx + ctx_offset);
                abs_gt2_flag = ovcabac_ae_read(cabac_ctx, abs_gt2_ctx + ctx_offset);
                coeff_val = 2 + par_lvl_flag;
                nb_rem_bins -= 2;
                if (abs_gt2_flag){
                    coeff_val += 2; /* (abs_gt2_flag << 1) */
                    gt2_idx_map[nb_pass2++] = idx;
                }
            }

            dep_quant_map <<= 1;
            dep_quant_map |= (prev_state >> 1);

            --nb_rem_bins;
            sig_idx_map[nb_sig_c++] = idx;

            coeffs[idx] = coeff_val;
            update_coeff_nbgh_first_pass(c_coding_ctx, tr_ctx_pos, coeff_val);
        }
        prev_state = (state_trans_tab >> (( prev_state << 2) + ((coeff_val & 1) << 1))) & 3;
        scan_map >>= 4;
        par_map  >>= 4;
        sig_map  >>= 4;
    }

    if (nb_pass2){
        decode_pass2_core(cabac_ctx, coeffs, nb_pass2, gt2_idx_map,
                          scan_ctx, c_coding_ctx);
    }

    if (scan_pos >= 0){
        decode_bypassed_coeff_core(cabac_ctx, coeffs, scan_pos, sig_idx_map,
                                   &nb_sig_c, scan_ctx, c_coding_ctx,
                                   &luma_dep_quant_ctx, &prev_state, &dep_quant_map);
    }

    decode_signs(cabac_ctx, coeffs, dep_quant_map, sig_idx_map, nb_sig_c);

    c_coding_ctx->nb_remaining_bins = nb_rem_bins;
    *state = prev_state;

    /*FIXME we could return state instead of nb_sig and avoid using pointers*/
    return nb_sig_c;
}

#if 0
static int
ovcabac_read_ae_sb_4x4_dc_dpq(OVCABACCtx *const cabac_ctx,
                         uint64_t *const ctx_table,
                         int16_t  *const coeffs,
                         int *const state, int start_coeff_idx,
                         VVCCoeffCodingCtx *const c_coding_ctx)
{
    residual_coding_first_subblock_4x4(cabac_ctx, ctx_table, coeffs, state,
                                       start_coeff_idx, 0xFAAAAA5555555555,
                                       0x8884444444444000, c_coding_ctx,
                                       &luma_ctx_offsets, &luma_dep_quant_ctx,
                                       &inv_diag_4x4_scan);
    return 0;
}
#endif

static int
ovcabac_read_ae_sb_4x4_first_dpq(OVCABACCtx *const cabac_ctx,
                            int16_t  *const coeffs,
                            int *const state, int start_coeff_idx,
                            uint8_t d_sb,
                            VVCCoeffCodingCtx *const c_coding_ctx)
{
    uint64_t par_flg_ofst_map = d_sb > 2 ? 0 : parity_flag_offset_map[d_sb];
    uint64_t sig_flg_ofst_map = d_sb > 2 ? 0 : sig_flag_offset_map[d_sb];

    residual_coding_first_subblock_4x4(cabac_ctx, coeffs, state,
                                       start_coeff_idx, par_flg_ofst_map,
                                       sig_flg_ofst_map, c_coding_ctx,
                                       &luma_ctx_offsets, &luma_dep_quant_ctx,
                                       &inv_diag_4x4_scan);
    return 0;
}

static int
ovcabac_read_ae_sb_4x4_dpq(OVCABACCtx *const cabac_ctx,
                           int16_t  *const coeffs,
                           int *const state,
                           uint8_t d_sb,
                           VVCCoeffCodingCtx *const c_coding_ctx)
{
    uint64_t par_flg_ofst_map = d_sb > 2 ? 0 : parity_flag_offset_map[d_sb];
    uint64_t sig_flg_ofst_map = d_sb > 2 ? 0 : sig_flag_offset_map[d_sb];

    residual_coding_subblock_4x4(cabac_ctx, coeffs, state,
                                 15, par_flg_ofst_map,
                                 sig_flg_ofst_map, c_coding_ctx,
                                 &luma_ctx_offsets, &luma_dep_quant_ctx,
                                 &inv_diag_4x4_scan);
    return 0;
}

static int
ovcabac_read_ae_sb_4x4_last_dc_dpq(OVCABACCtx *const cabac_ctx,
                              int16_t  *const coeffs,
                              int *const state,
                              VVCCoeffCodingCtx *const c_coding_ctx)
{
    residual_coding_subblock_dc(cabac_ctx, coeffs, state,
                                 15, 0xFAAAAA5555555555,
                                 0x8884444444444000, c_coding_ctx,
                                 &luma_ctx_offsets, &luma_dep_quant_ctx,
                                 &inv_diag_4x4_scan);
    return 0;
}

static int
ovcabac_read_ae_sb_dc_coeff_dpq(OVCABACCtx *const cabac_ctx,
                                int16_t *const coeffs)
{
    uint64_t *const ctx_table = cabac_ctx->ctx_table;
    uint8_t rem_abs_gt1_flag = ovcabac_ae_read(cabac_ctx, &ctx_table[GT0_FLAG_CTX_OFFSET]);

    uint32_t value = 1 + rem_abs_gt1_flag;

    uint8_t sign_flag;

    if (rem_abs_gt1_flag) {
        uint8_t par_level_flag   = ovcabac_ae_read(cabac_ctx, &ctx_table[PAR_FLAG_CTX_OFFSET]);
        uint8_t rem_abs_gt2_flag = ovcabac_ae_read(cabac_ctx, &ctx_table[GT1_FLAG_CTX_OFFSET]);

        value += (rem_abs_gt2_flag << 1) + par_level_flag;

        if (rem_abs_gt2_flag) {
            value += decode_truncated_rice(cabac_ctx, 0);
        }
    }

    sign_flag = ovcabac_bypass_read(cabac_ctx);

    /*FIXME might need a second pass*/
    //FIXME state = 0? find out shift
    //FIXME adapt this for int16_t and change coeff storage + apply inv quantif
    coeffs[0] = ( sign_flag ? -(int16_t)(value << 1) : (int16_t)(value << 1) );

    return 1;
}

static int
ovcabac_read_ae_sb_4x4_dc_c_dpq(OVCABACCtx *const cabac_ctx,
                                int16_t  *const coeffs,
                                int *const state, int start_coeff_idx,
                                VVCCoeffCodingCtx *const c_coding_ctx)
{
    residual_coding_first_subblock_4x4(cabac_ctx, coeffs, state,
                                       start_coeff_idx, 0x5000000000000000,
                                       0x4440000000000000, c_coding_ctx,
                                       &chroma_ctx_offsets, &chroma_dep_quant_ctx,
                                       &inv_diag_4x4_scan);
    return 0;
}

static int
ovcabac_read_ae_sb_4x4_first_c_dpq(OVCABACCtx *const cabac_ctx,
                                   int16_t  *const coeffs,
                                   int *const state, int start_coeff_idx,
                                   VVCCoeffCodingCtx *const c_coding_ctx)
{
    residual_coding_first_subblock_4x4(cabac_ctx, coeffs, state,
                                       start_coeff_idx, 0,
                                       0, c_coding_ctx,
                                       &chroma_ctx_offsets, &chroma_dep_quant_ctx,
                                       &inv_diag_4x4_scan);
    return 0;
}

static int
ovcabac_read_ae_sb_4x4_c_dpq(OVCABACCtx *const cabac_ctx,
                             int16_t  *const coeffs,
                             int *const state,
                             VVCCoeffCodingCtx *const c_coding_ctx)
{
    residual_coding_subblock_4x4(cabac_ctx, coeffs, state,
                                 0, 0,
                                 0, c_coding_ctx,
                                 &chroma_ctx_offsets, &chroma_dep_quant_ctx,
                                 &inv_diag_4x4_scan);
    return 0;
}

static int
ovcabac_read_ae_sb_4x4_last_dc_c_dpq(OVCABACCtx *const cabac_ctx,
                                     int16_t  *const coeffs,
                                     int *const state,
                                     VVCCoeffCodingCtx *const c_coding_ctx)
{
    residual_coding_subblock_dc(cabac_ctx, coeffs, state,
                                 15, 0x5000000000000000,
                                 0x4440000000000000, c_coding_ctx,
                                 &chroma_ctx_offsets, &chroma_dep_quant_ctx,
                                 &inv_diag_4x4_scan);
    return 0;
}

static int
ovcabac_read_ae_sb_8x2_dc_c_dpq(OVCABACCtx *const cabac_ctx,
                                int16_t  *const coeffs,
                                int *const state, int start_coeff_idx,
                                VVCCoeffCodingCtx *const c_coding_ctx)
{
    residual_coding_first_subblock_4x4(cabac_ctx, coeffs, state,
                                       start_coeff_idx, 0x5000000000000000,
                                       0x4440000000000000, c_coding_ctx,
                                       &chroma_ctx_offsets, &chroma_dep_quant_ctx,
                                       &inv_diag_8x2_scan);
    return 0;
}

static int
ovcabac_read_ae_sb_4x2_dc_c_dpq(OVCABACCtx *const cabac_ctx,
                                int16_t  *const coeffs,
                                int *const state, int start_coeff_idx,
                                VVCCoeffCodingCtx *const c_coding_ctx)
{
    residual_coding_first_subblock_4x4(cabac_ctx, coeffs, state,
                                       start_coeff_idx, 0x5000000000000000,
                                       0x4440000000000000, c_coding_ctx,
                                       &chroma_ctx_offsets, &chroma_dep_quant_ctx,
                                       &inv_diag_4x2_scan);
    return 0;
}

static int
ovcabac_read_ae_sb_2x8_dc_c_dpq(OVCABACCtx *const cabac_ctx,
                                int16_t  *const coeffs,
                                int *const state, int start_coeff_idx,
                                VVCCoeffCodingCtx *const c_coding_ctx)
{
    residual_coding_first_subblock_4x4(cabac_ctx, coeffs, state,
                                       start_coeff_idx, 0x5000000000000000,
                                       0x4440000000000000, c_coding_ctx,
                                       &chroma_ctx_offsets, &chroma_dep_quant_ctx,
                                       &inv_diag_2x8_scan);
    return 0;
}

static int
ovcabac_read_ae_sb_2x4_dc_c_dpq(OVCABACCtx *const cabac_ctx,
                                int16_t  *const coeffs,
                                int *const state, int start_coeff_idx,
                                VVCCoeffCodingCtx *const c_coding_ctx)
{
    residual_coding_first_subblock_4x4(cabac_ctx, coeffs, state,
                                       start_coeff_idx, 0x5000000000000000,
                                       0x4440000000000000, c_coding_ctx,
                                       &chroma_ctx_offsets, &chroma_dep_quant_ctx,
                                       &inv_diag_2x4_scan);
    return 0;
}

static int
ovcabac_read_ae_sb_2x2_c_dpq(OVCABACCtx *const cabac_ctx,
                                int16_t  *const coeffs,
                                int *const state, int start_coeff_idx,
                                VVCCoeffCodingCtx *const c_coding_ctx)
{
    residual_coding_first_subblock_4x4(cabac_ctx, coeffs, state,
                                       start_coeff_idx, 0x0,
                                       0x0, c_coding_ctx,
                                       &chroma_ctx_offsets, &chroma_dep_quant_ctx,
                                       &inv_diag_2x2_scan);
    return 0;
}

static int
ovcabac_read_ae_sb_2x2_dc_c_dpq(OVCABACCtx *const cabac_ctx,
                                int16_t  *const coeffs,
                                int *const state, int start_coeff_idx,
                                VVCCoeffCodingCtx *const c_coding_ctx)
{

    residual_coding_subblock_dc(cabac_ctx, coeffs, state,
                                3, 0x5000,
                                0x4440, c_coding_ctx,
                                &chroma_ctx_offsets, &chroma_dep_quant_ctx,
                                &inv_diag_2x2_scan);
    return 0;
}

static int
ovcabac_read_ae_sb_8x2_last_dc_c_dpq(OVCABACCtx *const cabac_ctx,
                                     int16_t  *const coeffs,
                                     int *const state,
                                     VVCCoeffCodingCtx *const c_coding_ctx)
{
    residual_coding_subblock_dc(cabac_ctx, coeffs, state,
                                 15, 0x5000000000000000,
                                 0x4440000000000000, c_coding_ctx,
                                 &chroma_ctx_offsets, &chroma_dep_quant_ctx,
                                 &inv_diag_8x2_scan);
    return 0;
}

static int
ovcabac_read_ae_sb_2x8_last_dc_c_dpq(OVCABACCtx *const cabac_ctx,
                                     int16_t  *const coeffs,
                                     int *const state,
                                     VVCCoeffCodingCtx *const c_coding_ctx)
{
    residual_coding_subblock_dc(cabac_ctx, coeffs, state,
                                 15, 0x5000000000000000,
                                 0x4440000000000000, c_coding_ctx,
                                 &chroma_ctx_offsets, &chroma_dep_quant_ctx,
                                 &inv_diag_2x8_scan);
    return 0;
}

static int
ovcabac_read_ae_sb_8x2_first_c_dpq(OVCABACCtx *const cabac_ctx,
                                   int16_t  *const coeffs,
                                   int *const state, int start_coeff_idx,
                                   VVCCoeffCodingCtx *const c_coding_ctx)
{
    residual_coding_first_subblock_4x4(cabac_ctx, coeffs, state,
                                       start_coeff_idx, 0,
                                       0, c_coding_ctx,
                                       &chroma_ctx_offsets, &chroma_dep_quant_ctx,
                                       &inv_diag_8x2_scan);
    return 0;
}

static int
ovcabac_read_ae_sb_2x8_first_c_dpq(OVCABACCtx *const cabac_ctx,
                                   int16_t  *const coeffs,
                                   int *const state, int start_coeff_idx,
                                   VVCCoeffCodingCtx *const c_coding_ctx)
{
    residual_coding_first_subblock_4x4(cabac_ctx, coeffs, state,
                                       start_coeff_idx, 0,
                                       0, c_coding_ctx,
                                       &chroma_ctx_offsets, &chroma_dep_quant_ctx,
                                       &inv_diag_2x8_scan);
    return 0;
}

static int
ovcabac_read_ae_sb_8x2_c_dpq(OVCABACCtx *const cabac_ctx,
                             int16_t  *const coeffs,
                             int *const state,
                             VVCCoeffCodingCtx *const c_coding_ctx)
{
    residual_coding_subblock_4x4(cabac_ctx, coeffs, state,
                                 0, 0,
                                 0, c_coding_ctx,
                                 &chroma_ctx_offsets, &chroma_dep_quant_ctx,
                                 &inv_diag_8x2_scan);
    return 0;
}

static int
ovcabac_read_ae_sb_2x8_c_dpq(OVCABACCtx *const cabac_ctx,
                             int16_t  *const coeffs,
                             int *const state,
                             VVCCoeffCodingCtx *const c_coding_ctx)
{
    residual_coding_subblock_4x4(cabac_ctx, coeffs, state,
                                 0, 0,
                                 0, c_coding_ctx,
                                 &chroma_ctx_offsets, &chroma_dep_quant_ctx,
                                 &inv_diag_2x8_scan);
    return 0;
}

static int
ovcabac_read_ae_sb_dc_coeff_c_dpq(OVCABACCtx *const cabac_ctx,
                                  int16_t *const coeffs)
{
    uint8_t sign_flag;
    uint64_t *const ctx_table = cabac_ctx->ctx_table;

    uint8_t rem_abs_gt1_flag = ovcabac_ae_read(cabac_ctx, &ctx_table[GT0_FLAG_C_CTX_OFFSET]);

    uint32_t value = 1 + rem_abs_gt1_flag;

    if (rem_abs_gt1_flag) {
        uint8_t par_level_flag   = ovcabac_ae_read(cabac_ctx, &ctx_table[PAR_FLAG_C_CTX_OFFSET]);
        uint8_t rem_abs_gt2_flag = ovcabac_ae_read(cabac_ctx, &ctx_table[GT1_FLAG_C_CTX_OFFSET]);

        value += (rem_abs_gt2_flag << 1) + par_level_flag;

        if (rem_abs_gt2_flag) {
            value += decode_truncated_rice(cabac_ctx, 0);
        }
    }

    sign_flag = ovcabac_bypass_read(cabac_ctx);

    //FIXME dep quant on dc
    //FIXME adapt this for int16_t and change coeff storage + apply inv quantif
    coeffs[0] = ( sign_flag ? -(int16_t)(value << 1) : (int16_t)(value << 1));

    return 1;
}

static int
ovcabac_read_ae_sb_8x2_dc_dpq(OVCABACCtx *const cabac_ctx,
                              int16_t  *const coeffs,
                              int *const state, int start_coeff_idx,
                              VVCCoeffCodingCtx *const c_coding_ctx)
{
    residual_coding_first_subblock_4x4(cabac_ctx, coeffs, state,
                                       start_coeff_idx, 0xFAAAA55555555555,
                                       0x8884444440000000, c_coding_ctx,
                                       &luma_ctx_offsets, &luma_dep_quant_ctx,
                                       &inv_diag_8x2_scan);
    return 0;
}

static int
ovcabac_read_ae_sb_2x8_dc_dpq(OVCABACCtx *const cabac_ctx,
                              int16_t  *const coeffs,
                              int *const state, int start_coeff_idx,
                              VVCCoeffCodingCtx *const c_coding_ctx)
{
    residual_coding_first_subblock_4x4(cabac_ctx, coeffs, state,
                                       start_coeff_idx, 0xFAAAA55555555555,
                                       0x8884444440000000, c_coding_ctx,
                                       &luma_ctx_offsets, &luma_dep_quant_ctx,
                                       &inv_diag_2x8_scan);
    return 0;
}

static int
ovcabac_read_ae_sb_1x16_dc_dpq(OVCABACCtx *const cabac_ctx,
                               int16_t  *const coeffs,
                               int *const state, int start_coeff_idx,
                               VVCCoeffCodingCtx *const c_coding_ctx)
{
    residual_coding_first_subblock_4x4(cabac_ctx, coeffs, state,
                                       start_coeff_idx, 0xFAA5555555000000,
                                       0x8844400000000000, c_coding_ctx,
                                       &luma_ctx_offsets, &luma_dep_quant_ctx,
                                       &inv_diag_1x16_scan);
    return 0;
}

static int
ovcabac_read_ae_sb_8x2_last_dc_dpq(OVCABACCtx *const cabac_ctx,
                                   int16_t  *const coeffs,
                                   int *const state,
                                   VVCCoeffCodingCtx *const c_coding_ctx)
{
    residual_coding_subblock_dc(cabac_ctx, coeffs, state,
                                 15, 0xFAAAA55555555555,
                                 0x8884444440000000, c_coding_ctx,
                                 &luma_ctx_offsets, &luma_dep_quant_ctx,
                                 &inv_diag_8x2_scan);
    return 0;
}

static int
ovcabac_read_ae_sb_2x8_last_dc_dpq(OVCABACCtx *const cabac_ctx,
                                   int16_t  *const coeffs,
                                   int *const state,
                                   VVCCoeffCodingCtx *const c_coding_ctx)
{
    residual_coding_subblock_dc(cabac_ctx, coeffs, state,
                                 15, 0xFAAAA55555555555,
                                 0x8884444440000000, c_coding_ctx,
                                 &luma_ctx_offsets, &luma_dep_quant_ctx,
                                 &inv_diag_2x8_scan);
    return 0;
}

static int
ovcabac_read_ae_sb_1x16_last_dc_dpq(OVCABACCtx *const cabac_ctx,
                                    int16_t  *const coeffs,
                                    int *const state,
                                    VVCCoeffCodingCtx *const c_coding_ctx)
{
    residual_coding_subblock_dc(cabac_ctx, coeffs, state,
                                15, 0xFAA5555555000000,
                                0x8844400000000000, c_coding_ctx,
                                &luma_ctx_offsets, &luma_dep_quant_ctx,
                                &inv_diag_1x16_scan);
    return 0;
}

static int
ovcabac_read_ae_sb_8x2_first_dpq(OVCABACCtx *const cabac_ctx,
                                 int16_t  *const coeffs,
                                 int *const state, int start_coeff_idx,
                                 VVCCoeffCodingCtx *const c_coding_ctx)
{
    residual_coding_first_subblock_4x4(cabac_ctx, coeffs, state,
                                       start_coeff_idx, 0x5550000000000000,
                                       0, c_coding_ctx,
                                       &luma_ctx_offsets, &luma_dep_quant_ctx,
                                       &inv_diag_8x2_scan);
    return 0;
}

static int
ovcabac_read_ae_sb_2x8_first_dpq(OVCABACCtx *const cabac_ctx,
                                 int16_t  *const coeffs,
                                 int *const state, int start_coeff_idx,
                                 VVCCoeffCodingCtx *const c_coding_ctx)
{
    residual_coding_first_subblock_4x4(cabac_ctx, coeffs, state,
                                       start_coeff_idx, 0x5550000000000000,
                                       0, c_coding_ctx,
                                       &luma_ctx_offsets, &luma_dep_quant_ctx,
                                       &inv_diag_2x8_scan);
    return 0;
}

static int
ovcabac_read_ae_sb_1x16_first_dpq(OVCABACCtx *const cabac_ctx,
                                  int16_t  *const coeffs,
                                  int *const state, int start_coeff_idx,
                                  VVCCoeffCodingCtx *const c_coding_ctx)
{
    residual_coding_first_subblock_4x4(cabac_ctx, coeffs, state,
                                       start_coeff_idx, 0,
                                       0, c_coding_ctx,
                                       &luma_ctx_offsets, &luma_dep_quant_ctx,
                                       &inv_diag_1x16_scan);
    return 0;
}

static int
ovcabac_read_ae_sb_8x2_first_far_dpq(OVCABACCtx *const cabac_ctx,
                                     int16_t  *const coeffs,
                                     int *const state, int start_coeff_idx,
                                     VVCCoeffCodingCtx *const c_coding_ctx)
{
    residual_coding_first_subblock_4x4(cabac_ctx, coeffs, state,
                                       start_coeff_idx, 0,
                                       0, c_coding_ctx,
                                       &luma_ctx_offsets, &luma_dep_quant_ctx,
                                       &inv_diag_8x2_scan);
    return 0;
}

static int
ovcabac_read_ae_sb_2x8_first_far_dpq(OVCABACCtx *const cabac_ctx,
                                     int16_t  *const coeffs,
                                     int *const state, int start_coeff_idx,
                                     VVCCoeffCodingCtx *const c_coding_ctx)
{
    residual_coding_first_subblock_4x4(cabac_ctx, coeffs, state,
                                       start_coeff_idx, 0,
                                       0, c_coding_ctx,
                                       &luma_ctx_offsets, &luma_dep_quant_ctx,
                                       &inv_diag_2x8_scan);
    return 0;
}

static int
ovcabac_read_ae_sb_8x2_dpq(OVCABACCtx *const cabac_ctx,
                           int16_t  *const coeffs,
                           int *const state,
                           VVCCoeffCodingCtx *const c_coding_ctx)
{
    residual_coding_subblock_4x4(cabac_ctx, coeffs, state,
                                 0, 0x5550000000000000,
                                 0, c_coding_ctx,
                                 &luma_ctx_offsets, &luma_dep_quant_ctx,
                                 &inv_diag_8x2_scan);
    return 0;
}

static int
ovcabac_read_ae_sb_2x8_dpq(OVCABACCtx *const cabac_ctx,
                           int16_t  *const coeffs,
                           int *const state,
                           VVCCoeffCodingCtx *const c_coding_ctx)
{
    residual_coding_subblock_4x4(cabac_ctx, coeffs, state,
                                 0, 0x5550000000000000,
                                 0, c_coding_ctx,
                                 &luma_ctx_offsets, &luma_dep_quant_ctx,
                                 &inv_diag_2x8_scan);
    return 0;
}

#if 0
static int
ovcabac_read_ae_sb_1x16_dpq(OVCABACCtx *const cabac_ctx,
                       int16_t  *const coeffs,
                       int *const state,
                       VVCCoeffCodingCtx *const c_coding_ctx)
{
    residual_coding_subblock_4x4(cabac_ctx, coeffs, state,
                                 0, 0,
                                 0, c_coding_ctx,
                                 &luma_ctx_offsets, &luma_dep_quant_ctx,
                                 &inv_diag_1x16_scan);
    return 0;
}
#endif

static int
ovcabac_read_ae_sb_8x2_far_dpq(OVCABACCtx *const cabac_ctx,
                               int16_t  *const coeffs,
                               int *const state,
                               VVCCoeffCodingCtx *const c_coding_ctx)
{
    residual_coding_subblock_4x4(cabac_ctx, coeffs, state,
                                 0, 0,
                                 0, c_coding_ctx,
                                 &luma_ctx_offsets, &luma_dep_quant_ctx,
                                 &inv_diag_8x2_scan);
    return 0;
}

static int
ovcabac_read_ae_sb_2x8_far_dpq(OVCABACCtx *const cabac_ctx,
                               int16_t  *const coeffs,
                               int *const state,
                               VVCCoeffCodingCtx *const c_coding_ctx)
{
    residual_coding_subblock_4x4(cabac_ctx, coeffs, state,
                                 0, 0,
                                 0, c_coding_ctx,
                                 &luma_ctx_offsets, &luma_dep_quant_ctx,
                                 &inv_diag_2x8_scan);
    return 0;
}

static void
update_ts_neighbourhood_first_pass(uint8_t   *const nb_significant,
                                   uint8_t   *const sign_map,
                                   uint16_t  *const abs_coeffs,
                                   int x, int y,
                                   int value, int sign){

    nb_significant[x + 1 +  y      * VVC_TR_CTX_STRIDE] += 1;
    nb_significant[x     + (y + 1) * VVC_TR_CTX_STRIDE] += 1;

    sign_map[x + 1 +  y      * VVC_TR_CTX_STRIDE] += sign;
    sign_map[x     + (y + 1) * VVC_TR_CTX_STRIDE] += sign;

    abs_coeffs[x +  y * VVC_TR_CTX_STRIDE] =  value;
}

static uint8_t
decode_pass2_ts(OVCABACCtx *const cabac_ctx,
                int16_t *const sb_coeffs,
                uint8_t nb_pass2,
                uint8_t *const pass2_idx_map,
                uint8_t *nb_pass3_p,
                uint8_t *const pass3_idx_map,
                uint16_t *const abs_coeffs,
                int16_t *const nb_remaining_bins,
                uint8_t log2_sb_w){
    int nb_read;
    int nb_pass3 = 0;
    const int x_mask = ((1 << log2_sb_w) - 1);
    uint64_t *const ts_ctx_table = &cabac_ctx->ctx_table[TS_GTX_FLAG_CTX_OFFSET];

    for (nb_read = 0; nb_read < nb_pass2 && *nb_remaining_bins >= 4; ++nb_read){
        int idx = pass2_idx_map[nb_read];
        int x = idx & x_mask;
        int y = idx >> log2_sb_w;
        int val = 0;

        uint8_t ts_gt2_flag = ovcabac_ae_read(cabac_ctx, ts_ctx_table + 1);
        --(*nb_remaining_bins);
        if (ts_gt2_flag) {
            uint8_t ts_gt3_flag = ovcabac_ae_read(cabac_ctx, ts_ctx_table + 2);
            --(*nb_remaining_bins);
            val += 2;
            if (ts_gt3_flag) {
                uint8_t ts_gt4_flag = ovcabac_ae_read(cabac_ctx, ts_ctx_table + 3);
                --(*nb_remaining_bins);
                val += 2;
                if (ts_gt4_flag) {
                    uint8_t ts_gt5_flag = ovcabac_ae_read(cabac_ctx, ts_ctx_table + 4);
                    --(*nb_remaining_bins);
                    val += 2;
                    if(ts_gt5_flag){
                        val += 2;
                        pass3_idx_map[nb_pass3++] = idx;
                    }
                }
            }
            sb_coeffs[idx] += val;
        }
        abs_coeffs[x +  y * VVC_TR_CTX_STRIDE] = sb_coeffs[idx];
    }
    *nb_pass3_p = nb_pass3;
    return nb_read;
}

static int
ovcabac_read_ae_sb_ts_core(OVCABACCtx *const cabac_ctx,
                           int16_t  *const coeffs,
                           const TSCoeffCodingCtx *cctx,
                           int16_t *const nb_remaining_bins,
                           const VVCSBScanContext *const scan_ctx)
{
    uint64_t *const ctx_table = cabac_ctx->ctx_table;
    const uint8_t log2_sb_w = scan_ctx->log2_sb_w;
    uint64_t scan_map = scan_ctx->scan_map;

    const int x_mask = (1 << log2_sb_w) - 1;
    const uint8_t max_scan_pos = (1 << (log2_sb_w + scan_ctx->log2_sb_h)) - 1;

    uint8_t sig_c_idx_map[16];
    uint8_t pass2_idx_map[16];
    uint8_t pass3_idx_map[16];

    uint8_t nb_sig_c = 0;
    uint8_t nb_pass2 = 0;
    uint8_t nb_pass3 = 0;

    uint8_t nb_read_pass2 = 0;

    uint16_t sign_map = 0;

    int coeff_idx;

    for (coeff_idx = 0; coeff_idx < max_scan_pos && *nb_remaining_bins >= 4; ++coeff_idx) {

        int idx = scan_map & 0xF;

        int x = idx & x_mask;
        int y = idx >> log2_sb_w;

        uint8_t nb_sig_c_ngh = cctx->nb_sig_ngh[x + y * VVC_TR_CTX_STRIDE];

        uint8_t ts_sig_c_flag = ovcabac_ae_read(cabac_ctx, &ctx_table[TS_SIG_FLAG_CTX_OFFSET + nb_sig_c_ngh]);

        --(*nb_remaining_bins);

        if (ts_sig_c_flag) {

            /* FIXME better computation of sign_offset */
            uint8_t nb_signs_ngh = cctx->sign_map[x + y * VVC_TR_CTX_STRIDE];
            int sign_offset = nb_sig_c_ngh != 2 ? nb_sig_c_ngh + nb_signs_ngh :
                (nb_signs_ngh == 2 ? 2 : nb_signs_ngh ^ 1);

            uint8_t ts_sign_flag = ovcabac_ae_read(cabac_ctx, &ctx_table[TS_RESIDUAL_SIGN_CTX_OFFSET + sign_offset]);
            uint8_t ts_gt1_flag  = ovcabac_ae_read(cabac_ctx, &ctx_table[TS_LRG1_FLAG_CTX_OFFSET + nb_sig_c_ngh]);
            int value = 1;

            sign_map |= ts_sign_flag << nb_sig_c;

            sig_c_idx_map[nb_sig_c++] = idx;

            --(*nb_remaining_bins);
            --(*nb_remaining_bins);

            if (ts_gt1_flag) {
                uint8_t ts_parity_flag = ovcabac_ae_read(cabac_ctx, &ctx_table[TS_PAR_FLAG_CTX_OFFSET]);

                value += 1 + ts_parity_flag;
                pass2_idx_map[nb_pass2++] = idx;

                --(*nb_remaining_bins);
            }

            coeffs[idx] = value;
            update_ts_neighbourhood_first_pass(cctx->nb_sig_ngh, cctx->sign_map,
                                               cctx->abs_coeffs, x, y, value, ts_sign_flag);
        }
        scan_map >>= 4;
    }

    if (*nb_remaining_bins >= 4) {
        int idx = scan_map & 0xF;

        int x = idx & x_mask;
        int y = idx >> log2_sb_w;

        uint8_t ts_sig_c_flag = !nb_sig_c;
        uint8_t nb_sig_c_ngh  = cctx->nb_sig_ngh[x + y * VVC_TR_CTX_STRIDE];

        if (nb_sig_c) {
            ts_sig_c_flag = ovcabac_ae_read(cabac_ctx, &ctx_table[TS_SIG_FLAG_CTX_OFFSET + nb_sig_c_ngh]);
            --(*nb_remaining_bins);
        }

        if (ts_sig_c_flag) {
            uint8_t nb_signs_ngh = cctx->sign_map[x + y * VVC_TR_CTX_STRIDE];
            int sign_offset = nb_sig_c_ngh != 2 ? nb_sig_c_ngh + nb_signs_ngh :
                                                 (nb_signs_ngh == 2 ? 2 : nb_signs_ngh ^ 1);

            uint8_t ts_sign_flag = ovcabac_ae_read(cabac_ctx, &ctx_table[TS_RESIDUAL_SIGN_CTX_OFFSET + sign_offset]);
            uint8_t ts_gt1_flag  = ovcabac_ae_read(cabac_ctx, &ctx_table[TS_LRG1_FLAG_CTX_OFFSET + nb_sig_c_ngh]);
            int value = 1;

            sign_map |= ts_sign_flag << nb_sig_c;

            sig_c_idx_map[nb_sig_c++] = max_scan_pos;

            --(*nb_remaining_bins);
            --(*nb_remaining_bins);

            if (ts_gt1_flag) {
                uint8_t ts_parity_flag = ovcabac_ae_read(cabac_ctx, &ctx_table[TS_PAR_FLAG_CTX_OFFSET]);
                value += 1 + ts_parity_flag;
                pass2_idx_map[nb_pass2++] = max_scan_pos;
                --(*nb_remaining_bins);
            }

            update_ts_neighbourhood_first_pass(cctx->nb_sig_ngh, cctx->sign_map, cctx->abs_coeffs,
                                               x, y, value, ts_sign_flag);

            coeffs[max_scan_pos] = value;
        }
        ++coeff_idx;
        scan_map >>= 4;
    }

    if (nb_pass2) {
        nb_read_pass2 = decode_pass2_ts(cabac_ctx, coeffs, nb_pass2, pass2_idx_map,
                                        &nb_pass3, pass3_idx_map, cctx->abs_coeffs, nb_remaining_bins,
                                        log2_sb_w);
    }

    /* Loop over already read pass2 requiring pass 3 */
    for (int i = 0; i < nb_pass3; i++) {
        int idx = pass3_idx_map[i];
        int remainder = decode_truncated_rice(cabac_ctx, 1);
        int x = idx & x_mask;
        int y = idx >> log2_sb_w;

        coeffs[idx] += remainder;
        cctx->abs_coeffs[x + y * VVC_TR_CTX_STRIDE] = coeffs[idx];
    }

    /* Loop over non already processed in pass 2 => bypass coded from pass 2 */
    for (int i = nb_read_pass2; i < nb_pass2; i++) {
        int idx = pass2_idx_map[i];
        int remainder = decode_truncated_rice(cabac_ctx, 1);
        int x = idx & x_mask;
        int y = idx >> log2_sb_w;

        coeffs[idx] += remainder;
        cctx->abs_coeffs[x + y * VVC_TR_CTX_STRIDE] = coeffs[idx];
    }

    /* FIXME we could probably move this to an other pass
     * based on the fact coeffs[idx] == 1 can only be true in first pass
     * when no greater flag is encountered
     */
    /* TS coeff Prediction */
    for (int i = 0; i < nb_sig_c; i++) {
        int idx = sig_c_idx_map[i];
        int x = idx & x_mask;
        int y = idx >> log2_sb_w;
        int max_abs_ngh = OVMAX(cctx->abs_coeffs[x     + ((y - 1) * VVC_TR_CTX_STRIDE)],
                                cctx->abs_coeffs[x - 1 + ( y      * VVC_TR_CTX_STRIDE)]);

        if (coeffs[idx] == 1 && max_abs_ngh) {
            coeffs[idx] = max_abs_ngh;
            cctx->abs_coeffs[x + y * VVC_TR_CTX_STRIDE] = coeffs[idx];
        } else {
            coeffs[idx] = coeffs[idx] - (coeffs[idx] <= max_abs_ngh);
            cctx->abs_coeffs[x + y * VVC_TR_CTX_STRIDE] = coeffs[idx];
        }
    }

    /* Bypass coded coefficients from pass 1 */
    for (int scan_idx = coeff_idx; scan_idx <= max_scan_pos; scan_idx++) {
        int idx = scan_map & 0xF;
        coeffs[idx] = decode_truncated_rice(cabac_ctx, 1) >> 1;

        if (coeffs[idx]) {
            uint8_t sign_flag = ovcabac_bypass_read(cabac_ctx);
            sign_map |= (sign_flag << nb_sig_c);
            sig_c_idx_map[nb_sig_c++] = idx;
        }
        scan_map >>= 4;
    }

    /* Apply coeff signs */
    for (int i = 0; i < nb_sig_c; i++) {
        int idx = sig_c_idx_map[i];

        coeffs[idx] = (sign_map & 0x1) ? -coeffs[idx] : coeffs[idx];
        sign_map = sign_map >> 1;
    }

    return 0;
}

static inline void
decode_bypassed_coeff_sdh(OVCABACCtx *const cabac_ctx,
                          int16_t *const coeffs, int last_scan_pos,
                          uint8_t *const sig_idx_map,
                          int *const nb_sig_c,
                          const VVCSBScanContext *const scan_ctx,
                          const VVCCoeffCodingCtx *const c_coding_ctx)
{
    int scan_pos;
    int rice_param;

    const uint64_t inv_diag_map = scan_ctx->scan_map;
    const uint8_t log2_sb_w     = scan_ctx->log2_sb_w;
    const uint8_t x_mask = (1 << log2_sb_w) - 1;

    int max_start = (1 << (scan_ctx->log2_sb_w + scan_ctx->log2_sb_h)) - 1;
    uint8_t pos_shift = ((max_start - (last_scan_pos & 0xF)) << 2);
    uint64_t scan_map = inv_diag_map >> pos_shift;

    for(scan_pos = last_scan_pos; scan_pos >= 0; --scan_pos){
        int idx = scan_map & 0xF;
        int value;

        int tr_ctx_pos = (idx & x_mask) + (idx >> log2_sb_w) * VVC_TR_CTX_STRIDE;

        int sum_abs = OVMIN(31, c_coding_ctx->sum_abs_lvl2[tr_ctx_pos]);

        rice_param = rice_param_tab[sum_abs];

        value = decode_truncated_rice(cabac_ctx, rice_param) >> 1;

        /*FIXME understand value saturation*/
        value = (value == ((0 < 2 ? 1: 2) << rice_param))
            ?  0 : (value  < ((0 < 2 ? 1: 2) << rice_param))
            ?  value + 1 : value;

        if(value){
            update_coeff_nbgh_bypassed(c_coding_ctx, tr_ctx_pos, value);

            coeffs[idx] = value;

            sig_idx_map[(*nb_sig_c)++] = idx;
        }
        scan_map >>= 4;
    }
}

static void inline
decode_signs_sdh(OVCABACCtx *const cabac_ctx, int16_t *const coeffs,
                 uint8_t *const sig_idx_map, int nb_sig_c,
                 uint8_t use_sdh)
{
    /*FIXME we could avoid sum_abs by xor on parity_flags */
    const uint32_t nb_signs = nb_sig_c - use_sdh;
    uint32_t signs_map  = 0;
    uint32_t sum_parity = 0;
    uint32_t nb_bins = nb_signs;

    while (nb_bins--){
        signs_map  = signs_map << 1;
        signs_map |= ovcabac_bypass_read(cabac_ctx);
    }

    signs_map  <<= 32 - nb_signs;

    for (unsigned k = 0; k < nb_signs; k++){
        int idx = sig_idx_map[k];
        int16_t abs_coeff = coeffs[idx];
        uint8_t sign = !!(signs_map & (1u << 31));
        sum_parity ^= abs_coeff;
        coeffs[idx] = sign ? -abs_coeff : abs_coeff ;
        signs_map <<= 1;
    }

    if (use_sdh){
        int idx = sig_idx_map[nb_signs];
        int16_t abs_coeff = coeffs[idx];
        sum_parity ^= abs_coeff;
        coeffs[idx] = (sum_parity & 1) ? -abs_coeff : abs_coeff;
    }
}

static inline int
residual_coding_first_subblock_sdh(OVCABACCtx *const cabac_ctx,
                                   int16_t  *const coeffs,
                                   int start_pos,
                                   uint64_t par_flag_offset_map,
                                   uint64_t sig_flag_offset_map,
                                   VVCCoeffCodingCtx *const c_coding_ctx,
                                   const VVCSBStates *const ctx_offsets,
                                   const VVCSBScanContext *const scan_ctx)
{
    const uint64_t inv_diag_map = scan_ctx->scan_map;
    const uint8_t log2_sb_w     = scan_ctx->log2_sb_w;
    const uint8_t x_mask    = (1 << log2_sb_w) - 1;

    uint64_t *const sig_flg_ctx = &cabac_ctx->ctx_table[ctx_offsets->sig_flg_ctx_offset];
    uint64_t *const abs_gt1_ctx = &cabac_ctx->ctx_table[ctx_offsets->abs_gt1_ctx_offset];
    uint64_t *const par_lvl_ctx = &cabac_ctx->ctx_table[ctx_offsets->par_lvl_ctx_offset];
    uint64_t *const abs_gt2_ctx = &cabac_ctx->ctx_table[ctx_offsets->abs_gt2_ctx_offset];

    uint8_t sig_idx_map[16];
    uint8_t gt2_idx_map[16];
    int nb_sig_c = 0;
    int nb_pass2 = 0;

    uint8_t par_lvl_flag, abs_gt1_flag, abs_gt2_flag;
    int32_t coeff_val;
    uint16_t tr_ctx_pos;

    int nb_rem_bins = c_coding_ctx->nb_remaining_bins;

    // Implicit first coeff
    int scan_pos = start_pos;
    int max_start = (1 << (scan_ctx->log2_sb_w + scan_ctx->log2_sb_h)) - 1;
    uint8_t pos_shift = ((max_start - (scan_pos & 0xF)) << 2);
    uint64_t scan_map = inv_diag_map        >> pos_shift;
    uint64_t par_map  = par_flag_offset_map >> pos_shift;
    uint64_t sig_map  = sig_flag_offset_map >> pos_shift;

    uint8_t idx = scan_map & 0xF;

    abs_gt1_flag = ovcabac_ae_read(cabac_ctx, abs_gt1_ctx);

    --nb_rem_bins;
    coeff_val = 1 + abs_gt1_flag;

    if(abs_gt1_flag){
        par_lvl_flag = ovcabac_ae_read(cabac_ctx, par_lvl_ctx);
        abs_gt2_flag = ovcabac_ae_read(cabac_ctx, abs_gt2_ctx);
        coeff_val += par_lvl_flag;
        nb_rem_bins -= 2;
        if(abs_gt2_flag){
           coeff_val += 2; //(rem_abs_gt2_flag << 1)
           gt2_idx_map[nb_pass2++] = idx;
        }
    }

    coeffs[idx] = coeff_val;
    sig_idx_map[nb_sig_c++] = idx;

    tr_ctx_pos = (idx & x_mask) + (idx >> log2_sb_w) * VVC_TR_CTX_STRIDE;

    set_implicit_coeff_ngbh(c_coding_ctx, tr_ctx_pos, coeff_val);

    --scan_pos;
    scan_map >>= 4;
    par_map  >>= 4;
    sig_map  >>= 4;

    // First pass
    for( ; scan_pos >= 0 && nb_rem_bins >= 4; --scan_pos ){

        uint8_t ctx_offset;
        uint8_t sig_coeff_flag;

        idx = scan_map & 0xF;
        tr_ctx_pos = (idx & x_mask) + (idx >> log2_sb_w) * VVC_TR_CTX_STRIDE;

        /*FIXME we could state ctx switch by same offset for chroma and luma
        */
        ctx_offset  = OVMIN(((c_coding_ctx->sum_abs_lvl[tr_ctx_pos] + 1) >> 1), 3);
        ctx_offset += sig_map & 0xF;

        sig_coeff_flag = ovcabac_ae_read(cabac_ctx, sig_flg_ctx + ctx_offset);

        coeff_val = sig_coeff_flag;
        --nb_rem_bins;

        if (sig_coeff_flag){

            ctx_offset  = 1;
            ctx_offset += OVMIN(c_coding_ctx->sum_sig_nbs[tr_ctx_pos], 4);
            ctx_offset += (par_map & 0xF);

            abs_gt1_flag = ovcabac_ae_read(cabac_ctx, abs_gt1_ctx + ctx_offset);

            if (abs_gt1_flag){
                par_lvl_flag = ovcabac_ae_read(cabac_ctx, par_lvl_ctx + ctx_offset);
                abs_gt2_flag = ovcabac_ae_read(cabac_ctx, abs_gt2_ctx + ctx_offset);
                coeff_val = 2 + par_lvl_flag;
                nb_rem_bins -= 2;
                if (abs_gt2_flag){
                    coeff_val += 2; /* (abs_gt2_flag << 1) */
                    gt2_idx_map[nb_pass2++] = idx;
                }
            }

            --nb_rem_bins;
            sig_idx_map[nb_sig_c++] = idx;

            coeffs[idx] = coeff_val;
            update_coeff_nbgh_first_pass(c_coding_ctx, tr_ctx_pos, coeff_val);
        }
        scan_map >>= 4;
        par_map  >>= 4;
        sig_map  >>= 4;
    }

    if (nb_pass2){
        decode_pass2_core(cabac_ctx, coeffs, nb_pass2, gt2_idx_map,
                          scan_ctx, c_coding_ctx);
    }

    if (scan_pos >= 0){
        decode_bypassed_coeff_sdh(cabac_ctx, coeffs, scan_pos, sig_idx_map,
                                   &nb_sig_c, scan_ctx, c_coding_ctx);
    }

    if (nb_sig_c){
        int max_start = (1 << (scan_ctx->log2_sb_w + scan_ctx->log2_sb_h)) - 1;
        uint8_t first_nz = (scan_ctx->scan_idx_map >> ((max_start - (sig_idx_map[0])) << 2)) & 0XF;
        uint8_t last_nz  = (scan_ctx->scan_idx_map >> ((max_start - (sig_idx_map[nb_sig_c - 1])) << 2)) & 0XF;
        uint8_t use_sdh =  c_coding_ctx->enable_sdh && (first_nz - last_nz) >= 4;
        decode_signs_sdh(cabac_ctx, coeffs, sig_idx_map, nb_sig_c, use_sdh);
    }

    c_coding_ctx->nb_remaining_bins = nb_rem_bins;

    /*FIXME we could return state instead of nb_sig and avoid using pointers*/
    return nb_sig_c;
}

static inline int
residual_coding_subblock_sdh(OVCABACCtx *const cabac_ctx,
                             int16_t  *const coeffs,
                             int start_pos,
                             uint64_t par_flag_offset_map,
                             uint64_t sig_flag_offset_map,
                             VVCCoeffCodingCtx *const c_coding_ctx,
                             const VVCSBStates *const ctx_offsets,
                             const VVCSBScanContext *const scan_ctx)
{
    const uint64_t inv_diag_map = scan_ctx->scan_map;
    const uint8_t log2_sb_w     = scan_ctx->log2_sb_w;
    const uint8_t x_mask    = (1 << log2_sb_w) - 1;

    uint64_t *const sig_flg_ctx = &cabac_ctx->ctx_table[ctx_offsets->sig_flg_ctx_offset];
    uint64_t *const abs_gt1_ctx = &cabac_ctx->ctx_table[ctx_offsets->abs_gt1_ctx_offset];
    uint64_t *const par_lvl_ctx = &cabac_ctx->ctx_table[ctx_offsets->par_lvl_ctx_offset];
    uint64_t *const abs_gt2_ctx = &cabac_ctx->ctx_table[ctx_offsets->abs_gt2_ctx_offset];

    uint8_t sig_idx_map[16];
    uint8_t gt2_idx_map[16];
    int nb_sig_c = 0;
    int nb_pass2 = 0;

    uint8_t par_lvl_flag, abs_gt1_flag, abs_gt2_flag;
    int32_t coeff_val;
    uint16_t tr_ctx_pos;

    int nb_rem_bins = c_coding_ctx->nb_remaining_bins;

    int max_start = (1 << (scan_ctx->log2_sb_w + scan_ctx->log2_sb_h)) - 1;
    int scan_pos = max_start;
    uint64_t scan_map = inv_diag_map;
    uint64_t par_map  = par_flag_offset_map;
    uint64_t sig_map  = sig_flag_offset_map;

    uint8_t idx;

    // First pass
    for ( ; scan_pos > 0 && nb_rem_bins >= 4; --scan_pos ){

        uint8_t ctx_offset;
        uint8_t sig_coeff_flag;

        idx = scan_map & 0xF;
        tr_ctx_pos = (idx & x_mask) + (idx >> log2_sb_w) * VVC_TR_CTX_STRIDE;

        /*FIXME we could state ctx switch by same offset for chroma and luma
        */
        ctx_offset  = OVMIN(((c_coding_ctx->sum_abs_lvl[tr_ctx_pos] + 1) >> 1), 3);
        ctx_offset += sig_map & 0xF;

        sig_coeff_flag = ovcabac_ae_read(cabac_ctx, sig_flg_ctx + ctx_offset);

        coeff_val = sig_coeff_flag;
        --nb_rem_bins;

        if (sig_coeff_flag){

            ctx_offset  = 1;
            ctx_offset += OVMIN(c_coding_ctx->sum_sig_nbs[tr_ctx_pos], 4);
            ctx_offset += (par_map & 0xF);

            abs_gt1_flag = ovcabac_ae_read(cabac_ctx, abs_gt1_ctx + ctx_offset);

            if (abs_gt1_flag){
                par_lvl_flag = ovcabac_ae_read(cabac_ctx, par_lvl_ctx + ctx_offset);
                abs_gt2_flag = ovcabac_ae_read(cabac_ctx, abs_gt2_ctx + ctx_offset);
                coeff_val = 2 + par_lvl_flag;
                nb_rem_bins -= 2;
                if (abs_gt2_flag){
                    coeff_val += 2; /* (abs_gt2_flag << 1) */
                    gt2_idx_map[nb_pass2++] = idx;
                }
            }

            --nb_rem_bins;
            sig_idx_map[nb_sig_c++] = idx;

            coeffs[idx] = coeff_val;
            update_coeff_nbgh_first_pass(c_coding_ctx, tr_ctx_pos, coeff_val);
        }
        scan_map >>= 4;
        par_map  >>= 4;
        sig_map  >>= 4;
    }

    if (scan_pos == 0 && nb_rem_bins >= 4){

        uint8_t ctx_offset;
        uint8_t sig_coeff_flag;
        sig_coeff_flag= 1;
        idx = scan_map & 0xF;
        tr_ctx_pos = (idx & x_mask) + (idx >> log2_sb_w) * VVC_TR_CTX_STRIDE;

        //decrease scan_pos so we know last sig_coeff was read in first pass or not
        --scan_pos;

        if (nb_sig_c){
            ctx_offset  = OVMIN(((c_coding_ctx->sum_abs_lvl[tr_ctx_pos] + 1) >> 1), 3);
            ctx_offset += sig_map & 0xF;

            sig_coeff_flag = ovcabac_ae_read(cabac_ctx, sig_flg_ctx + ctx_offset);

            --nb_rem_bins;
        }

        coeff_val = sig_coeff_flag;

        if (sig_coeff_flag){

            ctx_offset  = 1;
            ctx_offset += OVMIN(c_coding_ctx->sum_sig_nbs[tr_ctx_pos], 4);
            ctx_offset += (par_map & 0xF);

            abs_gt1_flag = ovcabac_ae_read(cabac_ctx, abs_gt1_ctx + ctx_offset);

            if (abs_gt1_flag){
                par_lvl_flag = ovcabac_ae_read(cabac_ctx, par_lvl_ctx + ctx_offset);
                abs_gt2_flag = ovcabac_ae_read(cabac_ctx, abs_gt2_ctx + ctx_offset);
                coeff_val = 2 + par_lvl_flag;
                nb_rem_bins -= 2;
                if (abs_gt2_flag){
                    coeff_val += 2; /* (abs_gt2_flag << 1) */
                    gt2_idx_map[nb_pass2++] = idx;
                }
            }

            --nb_rem_bins;
            sig_idx_map[nb_sig_c++] = idx;

            coeffs[idx] = coeff_val;
            update_coeff_nbgh_first_pass(c_coding_ctx, tr_ctx_pos, coeff_val);
        }
        scan_map >>= 4;
        par_map  >>= 4;
        sig_map  >>= 4;
    }

    if (nb_pass2){
        decode_pass2_core(cabac_ctx, coeffs, nb_pass2, gt2_idx_map,
                          scan_ctx, c_coding_ctx);
    }

    if (scan_pos >= 0){
        decode_bypassed_coeff_sdh(cabac_ctx, coeffs, scan_pos, sig_idx_map,
                                  &nb_sig_c, scan_ctx, c_coding_ctx);
    }

    if (nb_sig_c){
        int max_start = (1 << (scan_ctx->log2_sb_w + scan_ctx->log2_sb_h)) - 1;
        uint8_t first_nz = (scan_ctx->scan_idx_map >> ((max_start - (sig_idx_map[0])) << 2)) & 0XF;
        uint8_t last_nz  = (scan_ctx->scan_idx_map >> ((max_start - (sig_idx_map[nb_sig_c - 1])) << 2)) & 0XF;
        uint8_t use_sdh =  c_coding_ctx->enable_sdh && (first_nz - last_nz) >= 4;
        decode_signs_sdh(cabac_ctx, coeffs, sig_idx_map, nb_sig_c, use_sdh);
    }

    c_coding_ctx->nb_remaining_bins = nb_rem_bins;

    /*FIXME we could return state instead of nb_sig and avoid using pointers*/
    return nb_sig_c;
}

static inline int
residual_coding_subblock_dc_sdh(OVCABACCtx *const cabac_ctx,
                                int16_t  *const coeffs,
                                int start_pos,
                                uint64_t par_flag_offset_map,
                                uint64_t sig_flag_offset_map,
                                VVCCoeffCodingCtx *const c_coding_ctx,
                                const VVCSBStates *const ctx_offsets,
                                const VVCSBScanContext *const scan_ctx)
{
    const uint64_t inv_diag_map = scan_ctx->scan_map;
    const uint8_t log2_sb_w     = scan_ctx->log2_sb_w;
    const uint8_t x_mask    = (1 << log2_sb_w) - 1;

    uint64_t *const sig_flg_ctx = &cabac_ctx->ctx_table[ctx_offsets->sig_flg_ctx_offset];
    uint64_t *const abs_gt1_ctx = &cabac_ctx->ctx_table[ctx_offsets->abs_gt1_ctx_offset];
    uint64_t *const par_lvl_ctx = &cabac_ctx->ctx_table[ctx_offsets->par_lvl_ctx_offset];
    uint64_t *const abs_gt2_ctx = &cabac_ctx->ctx_table[ctx_offsets->abs_gt2_ctx_offset];

    uint8_t sig_idx_map[16];
    uint8_t gt2_idx_map[16];
    int nb_sig_c = 0;
    int nb_pass2 = 0;

    uint8_t par_lvl_flag, abs_gt1_flag, abs_gt2_flag;
    int32_t coeff_val;
    uint16_t tr_ctx_pos;

    int nb_rem_bins = c_coding_ctx->nb_remaining_bins;

    int max_start = (1 << (scan_ctx->log2_sb_w + scan_ctx->log2_sb_h)) - 1;
    int scan_pos = max_start;
    uint64_t scan_map = inv_diag_map;
    uint64_t par_map  = par_flag_offset_map;
    uint64_t sig_map  = sig_flag_offset_map;

    uint8_t idx;

    // First pass
    for ( ; scan_pos > 0 && nb_rem_bins >= 4; --scan_pos ){

        uint8_t ctx_offset;
        uint8_t sig_coeff_flag;

        idx = scan_map & 0xF;
        tr_ctx_pos = (idx & x_mask) + (idx >> log2_sb_w) * VVC_TR_CTX_STRIDE;

        /*FIXME we could state ctx switch by same offset for chroma and luma
        */
        ctx_offset  = OVMIN(((c_coding_ctx->sum_abs_lvl[tr_ctx_pos] + 1) >> 1), 3);
        ctx_offset += sig_map & 0xF;

        sig_coeff_flag = ovcabac_ae_read(cabac_ctx, sig_flg_ctx + ctx_offset);

        coeff_val = sig_coeff_flag;
        --nb_rem_bins;

        if (sig_coeff_flag){

            ctx_offset  = 1;
            ctx_offset += OVMIN(c_coding_ctx->sum_sig_nbs[tr_ctx_pos], 4);
            ctx_offset += (par_map & 0xF);

            abs_gt1_flag = ovcabac_ae_read(cabac_ctx, abs_gt1_ctx + ctx_offset);

            if (abs_gt1_flag){
                par_lvl_flag = ovcabac_ae_read(cabac_ctx, par_lvl_ctx + ctx_offset);
                abs_gt2_flag = ovcabac_ae_read(cabac_ctx, abs_gt2_ctx + ctx_offset);
                coeff_val = 2 + par_lvl_flag;
                nb_rem_bins -= 2;
                if (abs_gt2_flag){
                    coeff_val += 2; /* (abs_gt2_flag << 1) */
                    gt2_idx_map[nb_pass2++] = idx;
                }
            }

            --nb_rem_bins;
            sig_idx_map[nb_sig_c++] = idx;

            coeffs[idx] = coeff_val;
            update_coeff_nbgh_first_pass(c_coding_ctx, tr_ctx_pos, coeff_val);
        }
        scan_map >>= 4;
        par_map  >>= 4;
        sig_map  >>= 4;
    }

    if (scan_pos == 0 && nb_rem_bins >= 4){

        uint8_t ctx_offset;
        uint8_t sig_coeff_flag;
        sig_coeff_flag= 1;
        idx = scan_map & 0xF;
        tr_ctx_pos = (idx & x_mask) + (idx >> log2_sb_w) * VVC_TR_CTX_STRIDE;

        //decrease scan_pos so we know last sig_coeff was read in first pass or not
        --scan_pos;

        ctx_offset  = OVMIN(((c_coding_ctx->sum_abs_lvl[tr_ctx_pos] + 1) >> 1), 3);
        ctx_offset += sig_map & 0xF;

        sig_coeff_flag = ovcabac_ae_read(cabac_ctx, sig_flg_ctx + ctx_offset);

        --nb_rem_bins;

        coeff_val = sig_coeff_flag;

        if (sig_coeff_flag){

            ctx_offset  = 1;
            ctx_offset += OVMIN(c_coding_ctx->sum_sig_nbs[tr_ctx_pos], 4);
            ctx_offset += (par_map & 0xF);

            abs_gt1_flag = ovcabac_ae_read(cabac_ctx, abs_gt1_ctx + ctx_offset);

            if (abs_gt1_flag){
                par_lvl_flag = ovcabac_ae_read(cabac_ctx, par_lvl_ctx + ctx_offset);
                abs_gt2_flag = ovcabac_ae_read(cabac_ctx, abs_gt2_ctx + ctx_offset);
                coeff_val = 2 + par_lvl_flag;
                nb_rem_bins -= 2;
                if (abs_gt2_flag){
                    coeff_val += 2; /* (abs_gt2_flag << 1) */
                    gt2_idx_map[nb_pass2++] = idx;
                }
            }

            --nb_rem_bins;
            sig_idx_map[nb_sig_c++] = idx;

            coeffs[idx] = coeff_val;
            update_coeff_nbgh_first_pass(c_coding_ctx, tr_ctx_pos, coeff_val);
        }
        scan_map >>= 4;
        par_map  >>= 4;
        sig_map  >>= 4;
    }

    if (nb_pass2){
        decode_pass2_core(cabac_ctx, coeffs, nb_pass2, gt2_idx_map,
                          scan_ctx, c_coding_ctx);
    }

    if (scan_pos >= 0){
        decode_bypassed_coeff_sdh(cabac_ctx, coeffs, scan_pos, sig_idx_map,
                                   &nb_sig_c, scan_ctx, c_coding_ctx);
    }

    if (nb_sig_c){
        int max_start = (1 << (scan_ctx->log2_sb_w + scan_ctx->log2_sb_h)) - 1;
        uint8_t first_nz = (scan_ctx->scan_idx_map >> ((max_start - (sig_idx_map[0])) << 2)) & 0XF;
        uint8_t last_nz  = (scan_ctx->scan_idx_map >> ((max_start - (sig_idx_map[nb_sig_c - 1])) << 2)) & 0XF;
        uint8_t use_sdh =  c_coding_ctx->enable_sdh && (first_nz - last_nz) >= 4;
        decode_signs_sdh(cabac_ctx, coeffs, sig_idx_map, nb_sig_c, use_sdh);
    }

    c_coding_ctx->nb_remaining_bins = nb_rem_bins;

    /*FIXME we could return state instead of nb_sig and avoid using pointers*/
    return nb_sig_c;
}

#if 0
static int
ovcabac_read_ae_sb_4x4_dc_sdh(OVCABACCtx *const cabac_ctx,
                         int16_t  *const coeffs,
                         int *const state, int start_coeff_idx,
                         VVCCoeffCodingCtx *const c_coding_ctx)
{
    residual_coding_first_subblock_sdh(cabac_ctx, coeffs,
                                       start_coeff_idx, 0xFAAAAA5555555555,
                                       0x8884444444444000, c_coding_ctx,
                                       &luma_ctx_offsets,
                                       &inv_diag_4x4_scan);
    return 0;
}
#endif

static int
ovcabac_read_ae_sb_4x4_first_sdh(OVCABACCtx *const cabac_ctx,
                                 int16_t  *const coeffs,
                                 int start_coeff_idx,
                                 uint8_t d_sb,
                                 VVCCoeffCodingCtx *const c_coding_ctx)
{
    uint64_t par_flg_ofst_map = d_sb > 2 ? 0 : parity_flag_offset_map[d_sb];
    uint64_t sig_flg_ofst_map = d_sb > 2 ? 0 : sig_flag_offset_map[d_sb];

    residual_coding_first_subblock_sdh(cabac_ctx, coeffs,
                                       start_coeff_idx, par_flg_ofst_map,
                                       sig_flg_ofst_map, c_coding_ctx,
                                       &luma_ctx_offsets,
                                       &inv_diag_4x4_scan);
    return 0;
}

static int
ovcabac_read_ae_sb_4x4_sdh(OVCABACCtx *const cabac_ctx,
                           int16_t  *const coeffs,
                           uint8_t d_sb,
                           VVCCoeffCodingCtx *const c_coding_ctx)
{
    uint64_t par_flg_ofst_map = d_sb > 2 ? 0 : parity_flag_offset_map[d_sb];
    uint64_t sig_flg_ofst_map = d_sb > 2 ? 0 : sig_flag_offset_map[d_sb];

    residual_coding_subblock_sdh(cabac_ctx, coeffs,
                                 15, par_flg_ofst_map,
                                 sig_flg_ofst_map, c_coding_ctx,
                                 &luma_ctx_offsets,
                                 &inv_diag_4x4_scan);
    return 0;
}

static int
ovcabac_read_ae_sb_4x4_last_dc_sdh(OVCABACCtx *const cabac_ctx,
                                   int16_t  *const coeffs,
                                   VVCCoeffCodingCtx *const c_coding_ctx)
{
    residual_coding_subblock_dc_sdh(cabac_ctx, coeffs,
                                    15, 0xFAAAAA5555555555,
                                    0x8884444444444000, c_coding_ctx,
                                    &luma_ctx_offsets,
                                    &inv_diag_4x4_scan);
    return 0;
}

static int
ovcabac_read_ae_sb_dc_coeff_sdh(OVCABACCtx *const cabac_ctx,
                                int16_t *const coeffs)
{
    uint64_t *const ctx_table = cabac_ctx->ctx_table;
    uint8_t rem_abs_gt1_flag = ovcabac_ae_read(cabac_ctx, &ctx_table[GT0_FLAG_CTX_OFFSET]);

    uint32_t value = 1 + rem_abs_gt1_flag;
    uint8_t sign_flag;

    if (rem_abs_gt1_flag) {
        uint8_t par_level_flag   = ovcabac_ae_read(cabac_ctx, &ctx_table[PAR_FLAG_CTX_OFFSET]);
        uint8_t rem_abs_gt2_flag = ovcabac_ae_read(cabac_ctx, &ctx_table[GT1_FLAG_CTX_OFFSET]);

        value += (rem_abs_gt2_flag << 1) + par_level_flag;

        if (rem_abs_gt2_flag) {
            value += decode_truncated_rice(cabac_ctx,0);
        }
    }

    sign_flag = ovcabac_bypass_read(cabac_ctx);

    /*FIXME might need a second pass*/
    //FIXME state = 0? find out shift
    //FIXME adapt this for int16_t and change coeff storage + apply inv quantif
    coeffs[0] = ( sign_flag ? -(int16_t)(value) : (int16_t)(value) );

    return 1;
}

static int
ovcabac_read_ae_sb_4x4_dc_c_sdh(OVCABACCtx *const cabac_ctx,
                                int16_t  *const coeffs,
                                int start_coeff_idx,
                                VVCCoeffCodingCtx *const c_coding_ctx)
{
    residual_coding_first_subblock_sdh(cabac_ctx, coeffs,
                                       start_coeff_idx, 0x5000000000000000,
                                       0x4440000000000000, c_coding_ctx,
                                       &chroma_ctx_offsets,
                                       &inv_diag_4x4_scan);
    return 0;
}

static int
ovcabac_read_ae_sb_4x4_first_c_sdh(OVCABACCtx *const cabac_ctx,
                                   int16_t  *const coeffs,
                                   int start_coeff_idx,
                                   VVCCoeffCodingCtx *const c_coding_ctx)
{
    residual_coding_first_subblock_sdh(cabac_ctx, coeffs,
                                       start_coeff_idx, 0,
                                       0, c_coding_ctx,
                                       &chroma_ctx_offsets,
                                       &inv_diag_4x4_scan);
    return 0;
}

static int
ovcabac_read_ae_sb_4x4_c_sdh(OVCABACCtx *const cabac_ctx,
                             int16_t  *const coeffs,
                             VVCCoeffCodingCtx *const c_coding_ctx)
{
    residual_coding_subblock_sdh(cabac_ctx, coeffs,
                                 0, 0,
                                 0, c_coding_ctx,
                                 &chroma_ctx_offsets,
                                 &inv_diag_4x4_scan);
    return 0;
}

static int
ovcabac_read_ae_sb_4x4_last_dc_c_sdh(OVCABACCtx *const cabac_ctx,
                                     int16_t  *const coeffs,
                                     VVCCoeffCodingCtx *const c_coding_ctx)
{
    residual_coding_subblock_dc_sdh(cabac_ctx, coeffs,
                                    15, 0x5000000000000000,
                                    0x4440000000000000, c_coding_ctx,
                                    &chroma_ctx_offsets,
                                    &inv_diag_4x4_scan);
    return 0;
}

static int
ovcabac_read_ae_sb_8x2_dc_c_sdh(OVCABACCtx *const cabac_ctx,
                                int16_t  *const coeffs,
                                int start_coeff_idx,
                                VVCCoeffCodingCtx *const c_coding_ctx)
{
    residual_coding_first_subblock_sdh(cabac_ctx, coeffs,
                                       start_coeff_idx, 0x5000000000000000,
                                       0x4440000000000000, c_coding_ctx,
                                       &chroma_ctx_offsets,
                                       &inv_diag_8x2_scan);
    return 0;
}

static int
ovcabac_read_ae_sb_2x8_dc_c_sdh(OVCABACCtx *const cabac_ctx,
                                int16_t  *const coeffs,
                                int start_coeff_idx,
                                VVCCoeffCodingCtx *const c_coding_ctx)
{
    residual_coding_first_subblock_sdh(cabac_ctx, coeffs,
                                       start_coeff_idx, 0x5000000000000000,
                                       0x4440000000000000, c_coding_ctx,
                                       &chroma_ctx_offsets,
                                       &inv_diag_2x8_scan);
    return 0;
}

static int
ovcabac_read_ae_sb_8x2_last_dc_c_sdh(OVCABACCtx *const cabac_ctx,
                                     int16_t  *const coeffs,
                                     VVCCoeffCodingCtx *const c_coding_ctx)
{
    residual_coding_subblock_dc_sdh(cabac_ctx, coeffs,
                                    15, 0x5000000000000000,
                                    0x4440000000000000, c_coding_ctx,
                                    &chroma_ctx_offsets,
                                    &inv_diag_8x2_scan);
    return 0;
}

static int
ovcabac_read_ae_sb_2x8_last_dc_c_sdh(OVCABACCtx *const cabac_ctx,
                                     int16_t  *const coeffs,
                                     VVCCoeffCodingCtx *const c_coding_ctx)
{
    residual_coding_subblock_dc_sdh(cabac_ctx, coeffs,
                                    15, 0x5000000000000000,
                                    0x4440000000000000, c_coding_ctx,
                                    &chroma_ctx_offsets,
                                    &inv_diag_2x8_scan);
    return 0;
}

static int
ovcabac_read_ae_sb_8x2_first_c_sdh(OVCABACCtx *const cabac_ctx,
                                   int16_t  *const coeffs,
                                   int start_coeff_idx,
                                   VVCCoeffCodingCtx *const c_coding_ctx)
{
    residual_coding_first_subblock_sdh(cabac_ctx, coeffs,
                                       start_coeff_idx, 0,
                                       0, c_coding_ctx,
                                       &chroma_ctx_offsets,
                                       &inv_diag_8x2_scan);
    return 0;
}

static int
ovcabac_read_ae_sb_2x8_first_c_sdh(OVCABACCtx *const cabac_ctx,
                                   int16_t  *const coeffs,
                                   int start_coeff_idx,
                                   VVCCoeffCodingCtx *const c_coding_ctx)
{
    residual_coding_first_subblock_sdh(cabac_ctx, coeffs,
                                       start_coeff_idx, 0,
                                       0, c_coding_ctx,
                                       &chroma_ctx_offsets,
                                       &inv_diag_2x8_scan);
    return 0;
}

static int
ovcabac_read_ae_sb_8x2_c_sdh(OVCABACCtx *const cabac_ctx,
                             int16_t  *const coeffs,
                             VVCCoeffCodingCtx *const c_coding_ctx)
{
    residual_coding_subblock_sdh(cabac_ctx, coeffs,
                                 0, 0,
                                 0, c_coding_ctx,
                                 &chroma_ctx_offsets,
                                 &inv_diag_8x2_scan);
    return 0;
}

static int
ovcabac_read_ae_sb_2x8_c_sdh(OVCABACCtx *const cabac_ctx,
                             int16_t  *const coeffs,
                             VVCCoeffCodingCtx *const c_coding_ctx)
{
    residual_coding_subblock_sdh(cabac_ctx, coeffs,
                                 0, 0,
                                 0, c_coding_ctx,
                                 &chroma_ctx_offsets,
                                 &inv_diag_2x8_scan);
    return 0;
}

static int
ovcabac_read_ae_sb_dc_coeff_c_sdh(OVCABACCtx *const cabac_ctx,
                                  int16_t *const coeffs)
{
    uint64_t *const ctx_table = cabac_ctx->ctx_table;
    uint8_t rem_abs_gt1_flag = ovcabac_ae_read(cabac_ctx, &ctx_table[GT0_FLAG_C_CTX_OFFSET]);

    uint32_t value = 1 + rem_abs_gt1_flag;

    uint8_t sign_flag;

    if( rem_abs_gt1_flag ){
        uint8_t par_level_flag   = ovcabac_ae_read(cabac_ctx, &ctx_table[PAR_FLAG_C_CTX_OFFSET]);
        uint8_t rem_abs_gt2_flag = ovcabac_ae_read(cabac_ctx, &ctx_table[GT1_FLAG_C_CTX_OFFSET]);
        value += (rem_abs_gt2_flag << 1) + par_level_flag;
        if( rem_abs_gt2_flag ){
            //value += 2;//rem_abs_gt2_flag << 1
            value += decode_truncated_rice(cabac_ctx,0);
        }
    }

    sign_flag = ovcabac_bypass_read(cabac_ctx);

    //FIXME dep quant on dc
    //FIXME adapt this for int16_t and change coeff storage + apply inv quantif
    coeffs[0] = ( sign_flag ? -(int16_t)(value) : (int16_t)(value));

    return 1;
}

static     int
ovcabac_read_ae_sb_8x2_dc_sdh(OVCABACCtx *const cabac_ctx,
                              int16_t  *const coeffs,
                              int start_coeff_idx,
                              VVCCoeffCodingCtx *const c_coding_ctx)
{
    residual_coding_first_subblock_sdh(cabac_ctx, coeffs,
                                       start_coeff_idx, 0xFAAAA55555555555,
                                       0x8884444440000000, c_coding_ctx,
                                       &luma_ctx_offsets,
                                       &inv_diag_8x2_scan);
    return 0;
}
static int
ovcabac_read_ae_sb_2x8_dc_sdh(OVCABACCtx *const cabac_ctx,
                              int16_t  *const coeffs,
                              int start_coeff_idx,
                              VVCCoeffCodingCtx *const c_coding_ctx)
{
    residual_coding_first_subblock_sdh(cabac_ctx, coeffs,
                                       start_coeff_idx, 0xFAAAA55555555555,
                                       0x8884444440000000, c_coding_ctx,
                                       &luma_ctx_offsets,
                                       &inv_diag_2x8_scan);
    return 0;
}

static int
ovcabac_read_ae_sb_1x16_dc_sdh(OVCABACCtx *const cabac_ctx,
                               int16_t  *const coeffs,
                               int start_coeff_idx,
                               VVCCoeffCodingCtx *const c_coding_ctx)
{
    residual_coding_first_subblock_sdh(cabac_ctx, coeffs,
                                       start_coeff_idx, 0xFAA5555555000000,
                                       0x8844400000000000, c_coding_ctx,
                                       &luma_ctx_offsets,
                                       &inv_diag_1x16_scan);
    return 0;
}

static int
ovcabac_read_ae_sb_8x2_last_dc_sdh(OVCABACCtx *const cabac_ctx,
                                   int16_t  *const coeffs,
                                   VVCCoeffCodingCtx *const c_coding_ctx)
{
    residual_coding_subblock_dc_sdh(cabac_ctx, coeffs,
                                    15, 0xFAAAA55555555555,
                                    0x8884444440000000, c_coding_ctx,
                                    &luma_ctx_offsets,
                                    &inv_diag_8x2_scan);
    return 0;
}

static int
ovcabac_read_ae_sb_2x8_last_dc_sdh(OVCABACCtx *const cabac_ctx,
                                   int16_t  *const coeffs,
                                   VVCCoeffCodingCtx *const c_coding_ctx)
{
    residual_coding_subblock_dc_sdh(cabac_ctx, coeffs,
                                    15, 0xFAAAA55555555555,
                                    0x8884444440000000, c_coding_ctx,
                                    &luma_ctx_offsets,
                                    &inv_diag_2x8_scan);
    return 0;
}
static int
ovcabac_read_ae_sb_1x16_last_dc_sdh(OVCABACCtx *const cabac_ctx,
                                    int16_t  *const coeffs,
                                    VVCCoeffCodingCtx *const c_coding_ctx)
{
    residual_coding_subblock_dc_sdh(cabac_ctx, coeffs,
                                    15, 0xFAA5555555000000,
                                    0x8844400000000000, c_coding_ctx,
                                    &luma_ctx_offsets,
                                    &inv_diag_1x16_scan);
    return 0;
}
static int
ovcabac_read_ae_sb_8x2_first_sdh(OVCABACCtx *const cabac_ctx,
                                 int16_t  *const coeffs,
                                 int start_coeff_idx,
                                 VVCCoeffCodingCtx *const c_coding_ctx)
{
    residual_coding_first_subblock_sdh(cabac_ctx, coeffs,
                                       start_coeff_idx, 0x5550000000000000,
                                       0, c_coding_ctx,
                                       &luma_ctx_offsets,
                                       &inv_diag_8x2_scan);
    return 0;
}
static int
ovcabac_read_ae_sb_2x8_first_sdh(OVCABACCtx *const cabac_ctx,
                                 int16_t  *const coeffs,
                                 int start_coeff_idx,
                                 VVCCoeffCodingCtx *const c_coding_ctx)
{
    residual_coding_first_subblock_sdh(cabac_ctx, coeffs,
                                       start_coeff_idx, 0x5550000000000000,
                                       0, c_coding_ctx,
                                       &luma_ctx_offsets,
                                       &inv_diag_2x8_scan);
    return 0;
}
static int
ovcabac_read_ae_sb_1x16_first_sdh(OVCABACCtx *const cabac_ctx,
                                  int16_t  *const coeffs,
                                  int start_coeff_idx,
                                  VVCCoeffCodingCtx *const c_coding_ctx)
{
    residual_coding_first_subblock_sdh(cabac_ctx, coeffs,
                                       start_coeff_idx, 0,
                                       0, c_coding_ctx,
                                       &luma_ctx_offsets,
                                       &inv_diag_1x16_scan);
    return 0;
}
static int
ovcabac_read_ae_sb_8x2_first_far_sdh(OVCABACCtx *const cabac_ctx,
                                     int16_t  *const coeffs,
                                     int start_coeff_idx,
                                     VVCCoeffCodingCtx *const c_coding_ctx)
{
    residual_coding_first_subblock_sdh(cabac_ctx, coeffs,
                                       start_coeff_idx, 0,
                                       0, c_coding_ctx,
                                       &luma_ctx_offsets,
                                       &inv_diag_8x2_scan);
    return 0;
}

static int
ovcabac_read_ae_sb_2x8_first_far_sdh(OVCABACCtx *const cabac_ctx,
                                     int16_t  *const coeffs,
                                     int start_coeff_idx,
                                     VVCCoeffCodingCtx *const c_coding_ctx)
{
    residual_coding_first_subblock_sdh(cabac_ctx, coeffs,
                                       start_coeff_idx, 0,
                                       0, c_coding_ctx,
                                       &luma_ctx_offsets,
                                       &inv_diag_2x8_scan);
    return 0;
}

static int
ovcabac_read_ae_sb_8x2_sdh(OVCABACCtx *const cabac_ctx,
                           int16_t  *const coeffs,
                           VVCCoeffCodingCtx *const c_coding_ctx)
{
    residual_coding_subblock_sdh(cabac_ctx, coeffs,
                                 0, 0x5550000000000000,
                                 0, c_coding_ctx,
                                 &luma_ctx_offsets,
                                 &inv_diag_8x2_scan);
    return 0;
}

static int
ovcabac_read_ae_sb_2x8_sdh(OVCABACCtx *const cabac_ctx,
                           int16_t  *const coeffs,
                           VVCCoeffCodingCtx *const c_coding_ctx)
{
    residual_coding_subblock_sdh(cabac_ctx, coeffs,
                                 0, 0x5550000000000000,
                                 0, c_coding_ctx,
                                 &luma_ctx_offsets,
                                 &inv_diag_2x8_scan);
    return 0;
}

#if 0
static int
ovcabac_read_ae_sb_1x16_sdh(OVCABACCtx *const cabac_ctx,
                            int16_t  *const coeffs,
                            VVCCoeffCodingCtx *const c_coding_ctx)
{
    residual_coding_subblock_sdh(cabac_ctx, coeffs,
                                 0, 0,
                                 0, c_coding_ctx,
                                 &luma_ctx_offsets,
                                 &inv_diag_1x16_scan);
    return 0;
}
#endif

static int
ovcabac_read_ae_sb_8x2_far_sdh(OVCABACCtx *const cabac_ctx,
                               int16_t  *const coeffs,
                               VVCCoeffCodingCtx *const c_coding_ctx)
{
    residual_coding_subblock_sdh(cabac_ctx, coeffs,
                                 0, 0,
                                 0, c_coding_ctx,
                                 &luma_ctx_offsets,
                                 &inv_diag_8x2_scan);
    return 0;
}
static int
ovcabac_read_ae_sb_2x8_far_sdh(OVCABACCtx *const cabac_ctx,
                               int16_t  *const coeffs,
                               VVCCoeffCodingCtx *const c_coding_ctx)
{
    residual_coding_subblock_sdh(cabac_ctx, coeffs,
                                 0, 0,
                                 0, c_coding_ctx,
                                 &luma_ctx_offsets,
                                 &inv_diag_2x8_scan);
    return 0;
}

static const int inverse_quant_scale_lut[2][6] ={
    { 40, 45, 51, 57, 64,  72},
    { 57, 64, 72, 80, 90, 102}
};

struct IQScale{
    int scale;
    int shift;
    void (*dequant_sb)(int16_t *const sb_coeffs, int scale, int shift);
};

static void
dequant_sb_neg(int16_t *const sb_coeffs, int scale, int shift)
{
    const int     max_log2_tr_range = 15;
    const int32_t min_coeff_value   = -(1 << max_log2_tr_range);
    const int32_t max_coeff_value   =  (1 << max_log2_tr_range) - 1;

    for( int i = 0; i < 16 ; i++ ){
        sb_coeffs[i] = ov_clip((int32_t)(sb_coeffs[i] * scale) << shift ,
                min_coeff_value, max_coeff_value);
    }
}

static void
dequant_sb(int16_t *const sb_coeffs, int scale, int shift)
{
    const int     max_log2_tr_range = 15;
    const int32_t min_coeff_value   = -(1 << max_log2_tr_range);
    const int32_t max_coeff_value   =  (1 << max_log2_tr_range) - 1;

    int add = (1 << shift) >> 1;
    for( int i = 0; i < 16 ; i++ ){
        sb_coeffs[i] = ov_clip((int32_t)(sb_coeffs[i] * scale + add) >> shift ,
                min_coeff_value, max_coeff_value);
    }
}


static struct IQScale
derive_dequant_sdh(int qp, uint8_t log2_tb_w, uint8_t log2_tb_h)
{
    /*FIXME derive from ctx for range extentions*/
    const uint8_t max_log2_tr_range = 15;
    const int dep_quant_qp  = qp;
    const uint8_t log2_tb_s = log2_tb_w + log2_tb_h;
    struct IQScale dequant_params;
    int shift;
    int scale;
    /*FIXME non size dependent prefix could be derived from earlier ctx
      as soon as we know of bitdepth and tr range*/
    shift = IQUANT_SHIFT - (dep_quant_qp / 6)
              - (max_log2_tr_range - 10 - (log2_tb_s >> 1) - (log2_tb_s & 1));
    scale  = inverse_quant_scale_lut[log2_tb_s & 1][dep_quant_qp % 6];

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
    /*FIXME derive from ctx for range extentions*/
    const uint8_t max_log2_tr_range = 15;
    const int dep_quant_qp  = qp + 1;
    const uint8_t log2_tb_s = log2_tb_w + log2_tb_h;
    struct IQScale dequant_params;
    int shift;
    int scale;
    /*FIXME non size dependent prefix could be derived from earlier ctx
      as soon as we know of bitdepth and tr range*/
    shift = IQUANT_SHIFT + 1 - (dep_quant_qp / 6)
              - (max_log2_tr_range - 10 - (log2_tb_s >> 1) - (log2_tb_s & 1));
    scale  = inverse_quant_scale_lut[log2_tb_s & 1][dep_quant_qp % 6];

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
    /*FIXME derive from ctx for range extentions*/
    #if 0
    const uint8_t max_log2_tr_range = 15;
    #endif
    const int dep_quant_qp  = qp;
    struct IQScale dequant_params;
    int shift;
    int scale;
    /*FIXME non size dependent prefix could be derived from earlier ctx
      as soon as we know of bitdepth and tr range*/
    shift = IQUANT_SHIFT - (dep_quant_qp / 6) ;
    scale  = inverse_quant_scale_lut[0][dep_quant_qp % 6];

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

static void
reset_ctx_buffers (const VVCCoeffCodingCtx *ctx, int log2_w, int log2_h)
{
    uint8_t *nb_sig   = ctx->sum_sig_nbs;
    uint8_t *sum_abs1 = ctx->sum_abs_lvl;
    uint8_t *sum_abs2 = ctx->sum_abs_lvl2;
    for (int i = 0; i < (1 << log2_h); ++i){
        memset(nb_sig, 0, sizeof(*nb_sig) << log2_w);
        memset(sum_abs1, 0, sizeof(*sum_abs1) << log2_w);
        memset(sum_abs2, 0, sizeof(*sum_abs2) << log2_w);
        nb_sig   += VVC_TR_CTX_STRIDE;
        sum_abs1 += VVC_TR_CTX_STRIDE;
        sum_abs2 += VVC_TR_CTX_STRIDE;
    }
}

int
residual_coding_isp_h_sdh(OVCTUDec *const ctu_dec, int16_t *const dst,
                          unsigned int log2_tb_w, unsigned int log2_tb_h,
                          uint16_t last_pos)
{
    OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;
    int16_t *const _dst = dst;

    int last_x, last_y;
    int last_sb_x;
    int nb_sig_c;
    int start_coeff_idx;
    int x, y;
    int sb_offset;
    int sb_pos;
    int nb_sb;
    uint8_t sig_sb_flg = 1;
    int16_t sb_coeffs[16] = {0};
    uint8_t nb_significant [VVC_TR_CTX_SIZE];
    uint8_t sum_abs_level   [VVC_TR_CTX_SIZE];
    uint8_t sum_abs_level2 [VVC_TR_CTX_SIZE];

    uint16_t max_nb_bins = ((1 << (log2_tb_w + log2_tb_h))
                             * 28) >> 4;

    VVCCoeffCodingCtx c_coding_ctx = {
        .sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET],
        .sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET],
        .sum_sig_nbs  = &nb_significant[VVC_TR_CTX_OFFSET],
        .nb_remaining_bins = max_nb_bins,
        .enable_sdh = ctu_dec->enable_sdh
    };

    int qp = ctu_dec->dequant_luma.qp;

    struct IQScale deq_prms = derive_dequant_sdh(qp, log2_tb_w, log2_tb_h);

    memset(_dst, 0, sizeof(int16_t) * (1 << (log2_tb_w + log2_tb_h)));

    if (!last_pos) {
        ovcabac_read_ae_sb_dc_coeff_sdh(cabac_ctx, sb_coeffs);

        deq_prms.dequant_sb(sb_coeffs, deq_prms.scale, deq_prms.shift);

        _dst[0] = sb_coeffs[0];

        return 1;
    }

    reset_ctx_buffers(&c_coding_ctx, log2_tb_w, log2_tb_h);

    if(log2_tb_h){
        last_x =  last_pos       & 0x1F;
        last_y = (last_pos >> 8) & 0x1F;
        last_sb_x = last_x >> 3  ;

        if(!last_sb_x){
            int last_coeff_idx = last_x + (last_y << 3) ;
            int nb_coeffs = ff_vvc_diag_scan_8x2_num_cg [last_coeff_idx];
            nb_sig_c = ovcabac_read_ae_sb_8x2_dc_sdh(cabac_ctx, sb_coeffs,
                                                      nb_coeffs, &c_coding_ctx);

            deq_prms.dequant_sb(sb_coeffs, deq_prms.scale, deq_prms.shift);

            memcpy(&_dst[0             ], &sb_coeffs[0], sizeof(int16_t) * 8);
            memcpy(&_dst[1 << log2_tb_w], &sb_coeffs[8], sizeof(int16_t) * 8);

            return nb_sig_c;
        }

        x = last_x - (last_sb_x << 3) ;
        y = last_y;

        start_coeff_idx = ff_vvc_diag_scan_8x2_num_cg [x + y * (1 << 3)];

        sb_pos    = last_sb_x << 3;
        sb_offset = last_sb_x << 3;

        c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET + sb_offset];
        c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET + sb_offset];
        c_coding_ctx.sum_sig_nbs  = &nb_significant[VVC_TR_CTX_OFFSET + sb_offset];

        if(last_sb_x < 2){
            nb_sig_c = ovcabac_read_ae_sb_8x2_first_sdh(cabac_ctx, sb_coeffs,
                                                         start_coeff_idx, &c_coding_ctx);
        } else {
            nb_sig_c = ovcabac_read_ae_sb_8x2_first_far_sdh(cabac_ctx, sb_coeffs,
                                                             start_coeff_idx, &c_coding_ctx);
        }

        deq_prms.dequant_sb(sb_coeffs, deq_prms.scale, deq_prms.shift);

        memcpy(&_dst[sb_pos + 0               ], &sb_coeffs[0], sizeof(int16_t) * 8);
        memcpy(&_dst[sb_pos + (1 << log2_tb_w)], &sb_coeffs[8], sizeof(int16_t) * 8);

        nb_sb = last_sb_x;
        nb_sb--;

        for(;nb_sb > 0; --nb_sb){
            sig_sb_flg = ovcabac_read_ae_significant_sb_flag(cabac_ctx, sig_sb_flg);
            if(sig_sb_flg){

                memset(sb_coeffs, 0, sizeof(int16_t) * 16);

                sb_pos    = nb_sb << 3;
                sb_offset = nb_sb << 3;

                c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET + sb_offset];
                c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET + sb_offset];
                c_coding_ctx.sum_sig_nbs  = &nb_significant[VVC_TR_CTX_OFFSET + sb_offset];

                if(nb_sb < 2){
                    nb_sig_c = ovcabac_read_ae_sb_8x2_sdh(cabac_ctx, sb_coeffs,
                                                           &c_coding_ctx);
                } else {
                    nb_sig_c = ovcabac_read_ae_sb_8x2_far_sdh(cabac_ctx, sb_coeffs,
                                                               &c_coding_ctx);
                }

                deq_prms.dequant_sb(sb_coeffs, deq_prms.scale, deq_prms.shift);

                memcpy(&_dst[sb_pos + 0               ], &sb_coeffs[0],  sizeof(int16_t) * 8);
                memcpy(&_dst[sb_pos + (1 << log2_tb_w)], &sb_coeffs[8],  sizeof(int16_t) * 8);
            }
        }

        memset(sb_coeffs, 0, sizeof(uint16_t) * 16);

        c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET],
        c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET],
        c_coding_ctx.sum_sig_nbs  = &nb_significant[VVC_TR_CTX_OFFSET],

        nb_sig_c += ovcabac_read_ae_sb_8x2_last_dc_sdh(cabac_ctx, sb_coeffs,
                                                        &c_coding_ctx);

        deq_prms.dequant_sb(sb_coeffs, deq_prms.scale, deq_prms.shift);

        memcpy(&_dst[0             ], &sb_coeffs[0],  sizeof(int16_t) * 8);
        memcpy(&_dst[1 << log2_tb_w], &sb_coeffs[8],  sizeof(int16_t) * 8);

        return nb_sig_c;
    } else {
        last_x = last_pos & 0x1F;

        last_sb_x = last_x >> 4 ;

        if(!last_sb_x){
            int last_coeff_idx = last_x;
            int nb_coeffs = last_coeff_idx;

            nb_sig_c = ovcabac_read_ae_sb_1x16_dc_sdh(cabac_ctx, sb_coeffs,
                                                       nb_coeffs, &c_coding_ctx);

            deq_prms.dequant_sb(sb_coeffs, deq_prms.scale, deq_prms.shift);

            memcpy(&_dst[0] , &sb_coeffs[0], sizeof(int16_t) * 16);

            return nb_sig_c;
        }

        x = last_x - (last_sb_x << 4);

        sb_offset = last_sb_x << 4;
        sb_pos    = last_sb_x << 4;

        c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET + sb_offset];
        c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET + sb_offset];
        c_coding_ctx.sum_sig_nbs  = &nb_significant[VVC_TR_CTX_OFFSET + sb_offset];

        nb_sig_c = ovcabac_read_ae_sb_1x16_first_sdh(cabac_ctx, sb_coeffs,
                                                      x, &c_coding_ctx);

        deq_prms.dequant_sb(sb_coeffs, deq_prms.scale, deq_prms.shift);

        memcpy(&_dst[sb_pos] , &sb_coeffs[0], sizeof(int16_t) * 16);

        memset(sb_coeffs, 0, sizeof(uint16_t) * 16);

        c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET],
        c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET],
        c_coding_ctx.sum_sig_nbs  = &nb_significant[VVC_TR_CTX_OFFSET],

        nb_sig_c += ovcabac_read_ae_sb_1x16_last_dc_sdh(cabac_ctx, sb_coeffs,
                                                         &c_coding_ctx);

        deq_prms.dequant_sb(sb_coeffs, deq_prms.scale, deq_prms.shift);

        memcpy(&_dst[0] , &sb_coeffs[0], sizeof(int16_t) * 16);

        return nb_sig_c;
    }
}

int
residual_coding_isp_v_sdh(OVCTUDec *const ctu_dec, int16_t *const dst,
                          unsigned int log2_tb_w, unsigned int log2_tb_h,
                          uint16_t last_pos)
{
    OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;

    int16_t *const _dst = dst;
    int last_x, last_y;
    int last_sb_y;
    int nb_sig_c;
    int start_coeff_idx;
    int x, y;
    int sb_offset;
    int sb_pos;
    int nb_sb;
    uint8_t sig_sb_flg = 1;
    int16_t sb_coeffs[16] = {0};
    uint8_t  nb_significant[VVC_TR_CTX_SIZE];
    uint8_t  sum_abs_level  [VVC_TR_CTX_SIZE];
    uint8_t sum_abs_level2 [VVC_TR_CTX_SIZE];
    uint16_t max_nb_bins = ((1 << (log2_tb_w + log2_tb_h)) * 28) >> 4;

    VVCCoeffCodingCtx c_coding_ctx = {
        .sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET],
        .sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET],
        .sum_sig_nbs  = &nb_significant[VVC_TR_CTX_OFFSET],
        .nb_remaining_bins = max_nb_bins,
        .enable_sdh = ctu_dec->enable_sdh
    };

    int qp = ctu_dec->dequant_luma.qp;

    struct IQScale deq_prms = derive_dequant_sdh(qp, log2_tb_w, log2_tb_h);

    memset(_dst, 0, sizeof(int16_t) * (1 << (log2_tb_w + log2_tb_h)));

    if (!last_pos) {

        ovcabac_read_ae_sb_dc_coeff_sdh(cabac_ctx, sb_coeffs);

        deq_prms.dequant_sb(sb_coeffs, deq_prms.scale, deq_prms.shift);

        _dst[0] = sb_coeffs[0];
        return 1;
    }

    if (log2_tb_w) reset_ctx_buffers(&c_coding_ctx, log2_tb_w, log2_tb_h);
    else reset_ctx_buffers(&c_coding_ctx, log2_tb_h, log2_tb_w);

    if (log2_tb_w) {
        last_x =  last_pos       & 0x1F;
        last_y = (last_pos >> 8) & 0x1F;
        last_sb_y = last_y >> 3;

        if (!last_sb_y) {
            int last_coeff_idx = last_x + (last_y << 1);
            int nb_coeffs = ff_vvc_diag_scan_2x8_num_cg [last_coeff_idx];
            nb_sig_c = ovcabac_read_ae_sb_2x8_dc_sdh(cabac_ctx, sb_coeffs,
                                                      nb_coeffs, &c_coding_ctx);

            deq_prms.dequant_sb(sb_coeffs, deq_prms.scale, deq_prms.shift);

            memcpy(&_dst[0], &sb_coeffs[0],  sizeof(int16_t) * 16);

            return nb_sig_c;
        }

        sb_pos    =  last_sb_y << 4;
        sb_offset = (last_sb_y << 3) * VVC_TR_CTX_STRIDE;

        nb_sb = last_sb_y;

        x = last_x;
        y = last_y - (last_sb_y << 3);

        start_coeff_idx = ff_vvc_diag_scan_2x8_num_cg [x + (y << 1)];

        c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET + sb_offset];
        c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET + sb_offset];
        c_coding_ctx.sum_sig_nbs  = &nb_significant[VVC_TR_CTX_OFFSET + sb_offset];

        if(nb_sb < 2){
            nb_sig_c = ovcabac_read_ae_sb_2x8_first_sdh(cabac_ctx, sb_coeffs,
                                                         start_coeff_idx, &c_coding_ctx);
        } else {
            nb_sig_c = ovcabac_read_ae_sb_2x8_first_far_sdh(cabac_ctx, sb_coeffs,
                                                             start_coeff_idx, &c_coding_ctx);
        }

        deq_prms.dequant_sb(sb_coeffs, deq_prms.scale, deq_prms.shift);

        memcpy(&_dst[sb_pos + 0], &sb_coeffs[0], sizeof(int16_t) * 16);

        nb_sb--;

        for( ;nb_sb > 0; --nb_sb){
            sig_sb_flg = ovcabac_read_ae_significant_sb_flag(cabac_ctx, sig_sb_flg);
            if(sig_sb_flg){
                sb_pos    =  nb_sb << 4;
                sb_offset = (nb_sb << 3) * VVC_TR_CTX_STRIDE ;

                memset(sb_coeffs, 0, sizeof(int16_t) * 16);

                c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET + sb_offset];
                c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET + sb_offset];
                c_coding_ctx.sum_sig_nbs  = &nb_significant[VVC_TR_CTX_OFFSET + sb_offset];

                if(nb_sb < 2){
                    nb_sig_c = ovcabac_read_ae_sb_2x8_sdh(cabac_ctx, sb_coeffs,
                                                           &c_coding_ctx);
                } else {
                    nb_sig_c = ovcabac_read_ae_sb_2x8_far_sdh(cabac_ctx, sb_coeffs,
                                                               &c_coding_ctx);
                }

                deq_prms.dequant_sb(sb_coeffs, deq_prms.scale, deq_prms.shift);

                memcpy(&_dst[sb_pos], &sb_coeffs[0], sizeof(int16_t) * 16);
            }
        }

        memset(sb_coeffs, 0, sizeof(uint16_t) * 16);

        c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET];
        c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET];
        c_coding_ctx.sum_sig_nbs  = &nb_significant[VVC_TR_CTX_OFFSET];

        nb_sig_c += ovcabac_read_ae_sb_2x8_last_dc_sdh(cabac_ctx, sb_coeffs,
                                                        &c_coding_ctx);

        deq_prms.dequant_sb(sb_coeffs, deq_prms.scale, deq_prms.shift);

        memcpy(&_dst[0], &sb_coeffs[0],  sizeof(int16_t) * 16);

        return nb_sig_c;

    } else {
        last_y    = (last_pos >> 8) & 0x1F;
        last_sb_y = last_y >> 4;

        if(!last_sb_y){
            nb_sig_c = ovcabac_read_ae_sb_1x16_dc_sdh(cabac_ctx, sb_coeffs,
                                                       last_y, &c_coding_ctx);

            deq_prms.dequant_sb(sb_coeffs, deq_prms.scale, deq_prms.shift);

            memcpy(&_dst[0],  &sb_coeffs[0],  sizeof(int16_t) * 16);

            return nb_sig_c;
        }

        y = last_y - (last_sb_y << 4);

        sb_pos    = last_sb_y << 4;
        sb_offset = last_sb_y << 4;

        c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET + sb_offset];
        c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET + sb_offset];
        c_coding_ctx.sum_sig_nbs  = &nb_significant[VVC_TR_CTX_OFFSET + sb_offset];

        nb_sig_c = ovcabac_read_ae_sb_1x16_first_sdh(cabac_ctx, sb_coeffs,
                                                      y, &c_coding_ctx);

        deq_prms.dequant_sb(sb_coeffs, deq_prms.scale, deq_prms.shift);

        memcpy(&_dst[sb_pos + 0], &sb_coeffs[0],  sizeof(int16_t) * 16);

        memset(sb_coeffs, 0, sizeof(uint16_t) * 16);

        c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET];
        c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET];
        c_coding_ctx.sum_sig_nbs  = &nb_significant[VVC_TR_CTX_OFFSET];

        nb_sig_c += ovcabac_read_ae_sb_1x16_last_dc_sdh(cabac_ctx, sb_coeffs,
                                                         &c_coding_ctx);

        deq_prms.dequant_sb(sb_coeffs, deq_prms.scale, deq_prms.shift);

        memcpy(&_dst[0], &sb_coeffs[0], sizeof(int16_t) * 16);

        return nb_sig_c;
    }
}

uint64_t
residual_coding_sdh(OVCTUDec *const ctu_dec, int16_t *const dst,
                    unsigned int log2_tb_w, unsigned int log2_tb_h,
                    uint16_t last_pos)
{
    OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;
    int16_t *const _dst = dst;
    uint8_t lim_log2_w = OVMIN(log2_tb_w, 5);
    uint8_t lim_log2_h = OVMIN(log2_tb_h, 5);

    /*FIXME sort and rewrite scan map */
    const uint8_t *const sb_idx_2_sb_num = ff_vvc_idx_2_num[lim_log2_w  - 2]
                                                            [lim_log2_h - 2];

    const uint8_t *const scan_sb_x = ff_vvc_scan_x[lim_log2_w - 2][lim_log2_h - 2];
    const uint8_t *const scan_sb_y = ff_vvc_scan_y[lim_log2_w - 2][lim_log2_h - 2];
    int x, y;

    int last_x, last_y;
    int last_sb_x, last_sb_y;
    int nb_sb, nb_sig_c;
    int d_sb;

    int start_coeff_idx;
    int sb_offset;
    /*FIXME we could store coeff in smaller buffers and adapt transform
    to limited sizes*/
    uint16_t max_nb_bins = (((1 << (lim_log2_h + lim_log2_w)) << 5)
                          - ((1 << (lim_log2_h + lim_log2_w)) << 2)) >> 4;

    int16_t sb_coeffs[16] = {0};

    //TODO avoid offsets tabs
    uint8_t nb_significant[VVC_TR_CTX_SIZE];
    uint8_t sum_abs_level  [VVC_TR_CTX_SIZE];
    uint8_t sum_abs_level2 [VVC_TR_CTX_SIZE];

    VVCCoeffCodingCtx c_coding_ctx = {
        .sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET],
        .sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET],
        .sum_sig_nbs  = &nb_significant[VVC_TR_CTX_OFFSET],
        .nb_remaining_bins = max_nb_bins,
        .enable_sdh = ctu_dec->enable_sdh
    };

    uint64_t sig_sb_map = 0;
    int sb_pos;

    int qp = ctu_dec->dequant_luma.qp;

    struct IQScale deq_prms = derive_dequant_sdh(qp, log2_tb_w, log2_tb_h);

    memset(_dst, 0, sizeof(int16_t) * (1 << (log2_tb_w + log2_tb_h)));

    if (!last_pos){
        ovcabac_read_ae_sb_dc_coeff_sdh(cabac_ctx, sb_coeffs);

        deq_prms.dequant_sb(sb_coeffs, deq_prms.scale, deq_prms.shift);

        _dst[0] = sb_coeffs[0];
        memcpy(ctu_dec->lfnst_subblock, sb_coeffs, sizeof(int16_t) * 16);
        return 1;
    }

    reset_ctx_buffers(&c_coding_ctx, lim_log2_w, lim_log2_h);

    last_x =  last_pos       & 0x1F;
    last_y = (last_pos >> 8) & 0x1F;
    last_sb_x = last_x >> 2;
    last_sb_y = last_y >> 2;

    sig_sb_map |= 1llu << (last_sb_x + (last_sb_y << 3));

    if (!last_sb_x && !last_sb_y){
        int last_coeff_idx = last_x + (last_y << 2);
        int nb_coeffs = ff_vvc_diag_scan_4x4_num_cg [last_coeff_idx];

        nb_sig_c = ovcabac_read_ae_sb_4x4_first_sdh(cabac_ctx, sb_coeffs,
                                                     nb_coeffs, 0, &c_coding_ctx);

        deq_prms.dequant_sb(sb_coeffs, deq_prms.scale, deq_prms.shift);

        memcpy(&_dst[0]             , &sb_coeffs[ 0], sizeof(int16_t) * 4);
        memcpy(&_dst[1 << log2_tb_w], &sb_coeffs[ 4], sizeof(int16_t) * 4);
        memcpy(&_dst[2 << log2_tb_w], &sb_coeffs[ 8], sizeof(int16_t) * 4);
        memcpy(&_dst[3 << log2_tb_w], &sb_coeffs[12], sizeof(int16_t) * 4);
        memcpy(ctu_dec->lfnst_subblock, sb_coeffs, sizeof(int16_t) * 16);

        return sig_sb_map;
    }

    sb_pos    = (last_sb_x << 2) + ((last_sb_y << log2_tb_w) << 2);
    sb_offset = (last_sb_x << 2) + (last_sb_y << 2) * VVC_TR_CTX_STRIDE;

    d_sb = last_sb_x + last_sb_y;

    x = last_x & 3;
    y = last_y & 3;

    start_coeff_idx = ff_vvc_diag_scan_4x4_num_cg[x + (y << 2)];

    c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET + sb_offset];
    c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET + sb_offset];
    c_coding_ctx.sum_sig_nbs  = &nb_significant[VVC_TR_CTX_OFFSET + sb_offset];

    nb_sig_c = ovcabac_read_ae_sb_4x4_first_sdh(cabac_ctx, sb_coeffs,
                                                 start_coeff_idx, d_sb,
                                                 &c_coding_ctx);

    deq_prms.dequant_sb(sb_coeffs, deq_prms.scale, deq_prms.shift);

    memcpy(&_dst[sb_pos + (0)]             , &sb_coeffs[ 0], sizeof(int16_t) * 4);
    memcpy(&_dst[sb_pos + (1 << log2_tb_w)], &sb_coeffs[ 4], sizeof(int16_t) * 4);
    memcpy(&_dst[sb_pos + (2 << log2_tb_w)], &sb_coeffs[ 8], sizeof(int16_t) * 4);
    memcpy(&_dst[sb_pos + (3 << log2_tb_w)], &sb_coeffs[12], sizeof(int16_t) * 4);

    nb_sb = sb_idx_2_sb_num[last_sb_x + last_sb_y * ((1 << lim_log2_w) >> 2)];

    nb_sb--;

    for(; nb_sb > 0; --nb_sb){
        int x_sb = scan_sb_x[nb_sb];
        int y_sb = scan_sb_y[nb_sb];

        uint8_t sig_sb_flg;
        uint8_t sig_sb_blw = (((sig_sb_map >> ((y_sb + 1) << 3)) & 0xFF) >> x_sb) & 0x1;
        uint8_t sig_sb_rgt = (((sig_sb_map >> ( y_sb      << 3)) & 0xFF) >> (x_sb + 1)) & 0x1;

        sig_sb_flg = ovcabac_read_ae_significant_sb_flag(cabac_ctx, !!(sig_sb_rgt | sig_sb_blw));

        if(sig_sb_flg){

            memset(sb_coeffs, 0, sizeof(int16_t) * 16);

            sb_pos    = (x_sb << 2) + ((y_sb << log2_tb_w) << 2);
            sb_offset = (x_sb << 2) + (y_sb << 2) * (VVC_TR_CTX_STRIDE);

            sig_sb_map |= 1llu << (x_sb + (y_sb << 3));

            d_sb = x_sb + y_sb;

            c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET + sb_offset];
            c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET + sb_offset];
            c_coding_ctx.sum_sig_nbs  = &nb_significant[VVC_TR_CTX_OFFSET + sb_offset];

            nb_sig_c += ovcabac_read_ae_sb_4x4_sdh(cabac_ctx, sb_coeffs,
                                                    d_sb, &c_coding_ctx);

            deq_prms.dequant_sb(sb_coeffs, deq_prms.scale, deq_prms.shift);

            memcpy(&_dst[sb_pos + (0)]             , &sb_coeffs[ 0], sizeof(int16_t) * 4);
            memcpy(&_dst[sb_pos + (1 << log2_tb_w)], &sb_coeffs[ 4], sizeof(int16_t) * 4);
            memcpy(&_dst[sb_pos + (2 << log2_tb_w)], &sb_coeffs[ 8], sizeof(int16_t) * 4);
            memcpy(&_dst[sb_pos + (3 << log2_tb_w)], &sb_coeffs[12], sizeof(int16_t) * 4);
        }
    }

    memset(sb_coeffs, 0, sizeof(uint16_t) * 16);

    c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET];
    c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET];
    c_coding_ctx.sum_sig_nbs  = &nb_significant[VVC_TR_CTX_OFFSET];

    nb_sig_c += ovcabac_read_ae_sb_4x4_last_dc_sdh(cabac_ctx, sb_coeffs,
                                                    &c_coding_ctx);

    deq_prms.dequant_sb(sb_coeffs, deq_prms.scale, deq_prms.shift);

    memcpy(&_dst[0]             , &sb_coeffs[ 0], sizeof(int16_t) * 4);
    memcpy(&_dst[1 << log2_tb_w], &sb_coeffs[ 4], sizeof(int16_t) * 4);
    memcpy(&_dst[2 << log2_tb_w], &sb_coeffs[ 8], sizeof(int16_t) * 4);
    memcpy(&_dst[3 << log2_tb_w], &sb_coeffs[12], sizeof(int16_t) * 4);

    return sig_sb_map;
}

int
residual_coding_chroma_sdh(OVCTUDec *const ctu_dec, int16_t *const dst,
                           unsigned int log2_tb_w, unsigned int log2_tb_h,
                           uint16_t last_pos)
{
    OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;

    uint16_t max_nb_bins = ((1 << (log2_tb_w + log2_tb_h) << 5)
                          - ((1 << (log2_tb_w + log2_tb_h)) << 2)) >> 4;

    uint8_t nb_significant[VVC_TR_CTX_SIZE];
    uint8_t sum_abs_level  [VVC_TR_CTX_SIZE];
    uint8_t sum_abs_level2 [VVC_TR_CTX_SIZE];

    VVCCoeffCodingCtx c_coding_ctx = {
        .sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET],
        .sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET],
        .sum_sig_nbs  = &nb_significant[VVC_TR_CTX_OFFSET],
        .nb_remaining_bins = max_nb_bins,
        .enable_sdh = ctu_dec->enable_sdh
    };

    int16_t sb_coeffs[16] = {0}; //temporary table to store coeffs in process

    uint64_t sig_sb_map = 0;

    int qp = ctu_dec->dequant_chroma->qp;

    struct IQScale deq_prms = derive_dequant_sdh(qp, log2_tb_w, log2_tb_h);

    memset(dst, 0, sizeof(int16_t) * (1 << (log2_tb_w + log2_tb_h)));

    if (!last_pos) {
        ovcabac_read_ae_sb_dc_coeff_c_sdh(cabac_ctx, sb_coeffs);
        deq_prms.dequant_sb(sb_coeffs, deq_prms.scale, deq_prms.shift);
        dst[0] = sb_coeffs [0];
        return 0x1;
    }

    reset_ctx_buffers(&c_coding_ctx, log2_tb_w, log2_tb_h);

    if (log2_tb_w > 1 && log2_tb_h > 1) {
        const uint8_t *const sb_idx_2_sb_num = ff_vvc_idx_2_num[log2_tb_w - 2]
                                                               [log2_tb_h - 2];

        const uint8_t *const scan_sb_x = ff_vvc_scan_x[log2_tb_w - 2][log2_tb_h - 2];
        const uint8_t *const scan_sb_y = ff_vvc_scan_y[log2_tb_w - 2][log2_tb_h - 2];

        int nb_sig_coeff;
        int sb_offset;
        int sb_pos;

        int x =  last_pos       & 0x1F;
        int y = (last_pos >> 8) & 0x1F;

        int last_sb_x = x >> 2;
        int last_sb_y = y >> 2;

        int nb_coeffs = ff_vvc_diag_scan_4x4_num_cg[(x & 0x3) + ((y & 0x3) << 2)];

        int nb_sb = sb_idx_2_sb_num[last_sb_x + ((last_sb_y << log2_tb_w) >> 2)];

        if (!nb_sb) {
            ovcabac_read_ae_sb_4x4_dc_c_sdh(cabac_ctx, sb_coeffs,
                                            nb_coeffs, &c_coding_ctx);

            deq_prms.dequant_sb(sb_coeffs, deq_prms.scale, deq_prms.shift);

            memcpy(&dst[0]             , &sb_coeffs[ 0], sizeof(int16_t) * 4);
            memcpy(&dst[1 << log2_tb_w], &sb_coeffs[ 4], sizeof(int16_t) * 4);
            memcpy(&dst[2 << log2_tb_w], &sb_coeffs[ 8], sizeof(int16_t) * 4);
            memcpy(&dst[3 << log2_tb_w], &sb_coeffs[12], sizeof(int16_t) * 4);

            return 0x1;
        }

        sb_pos = (last_sb_x << 2) + (((last_sb_y << log2_tb_w) >> 2) << 4);

        sb_offset = (last_sb_x << 2) + (last_sb_y << 2) * (VVC_TR_CTX_STRIDE);

        c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET + sb_offset];
        c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET + sb_offset];
        c_coding_ctx.sum_sig_nbs  = &nb_significant[VVC_TR_CTX_OFFSET + sb_offset];

        nb_sig_coeff = ovcabac_read_ae_sb_4x4_first_c_sdh(cabac_ctx, sb_coeffs,
                                                           nb_coeffs,
                                                           &c_coding_ctx);

        deq_prms.dequant_sb(sb_coeffs, deq_prms.scale, deq_prms.shift);

        memcpy(&dst[sb_pos + (0)]             , &sb_coeffs[ 0], sizeof(int16_t) * 4);
        memcpy(&dst[sb_pos + (1 << log2_tb_w)], &sb_coeffs[ 4], sizeof(int16_t) * 4);
        memcpy(&dst[sb_pos + (2 << log2_tb_w)], &sb_coeffs[ 8], sizeof(int16_t) * 4);
        memcpy(&dst[sb_pos + (3 << log2_tb_w)], &sb_coeffs[12], sizeof(int16_t) * 4);

        sig_sb_map |= 1llu << (last_sb_x + (last_sb_y << 3));

        nb_sb--;

        for(int i = nb_sb; i > 0; --i){
            int x_sb = scan_sb_x[i];
            int y_sb = scan_sb_y[i];
            uint8_t sig_sb_blw = (((sig_sb_map >> ((y_sb + 1) << 3)) & 0xFF) >> x_sb) & 0x1;
            uint8_t sig_sb_rgt = (((sig_sb_map >> ( y_sb      << 3)) & 0xFF) >> (x_sb + 1)) & 0x1;

            uint8_t sig_sb_flg = ovcabac_read_ae_significant_sb_flag_chroma(cabac_ctx, !!(sig_sb_rgt | sig_sb_blw));

            if(sig_sb_flg){

                memset(sb_coeffs, 0, sizeof(int16_t) * 16);

                sb_pos = (x_sb << 2) + (((y_sb << log2_tb_w) >> 2) << 4);
                sb_offset = (x_sb << 2) + (y_sb << 2) * (VVC_TR_CTX_STRIDE);

                sig_sb_map |= 1llu << (x_sb + (y_sb << 3));

                c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET + sb_offset];
                c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET + sb_offset];
                c_coding_ctx.sum_sig_nbs  = &nb_significant[VVC_TR_CTX_OFFSET + sb_offset];

                nb_sig_coeff += ovcabac_read_ae_sb_4x4_c_sdh(cabac_ctx, sb_coeffs,
                                                              &c_coding_ctx);

                deq_prms.dequant_sb(sb_coeffs, deq_prms.scale, deq_prms.shift);

                memcpy(&dst[sb_pos + (0)]             , &sb_coeffs[ 0], sizeof(int16_t) * 4);
                memcpy(&dst[sb_pos + (1 << log2_tb_w)], &sb_coeffs[ 4], sizeof(int16_t) * 4);
                memcpy(&dst[sb_pos + (2 << log2_tb_w)], &sb_coeffs[ 8], sizeof(int16_t) * 4);
                memcpy(&dst[sb_pos + (3 << log2_tb_w)], &sb_coeffs[12], sizeof(int16_t) * 4);
            }
        }

        memset(sb_coeffs, 0, sizeof(int16_t) * 16);

        c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET];
        c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET];
        c_coding_ctx.sum_sig_nbs  = &nb_significant[VVC_TR_CTX_OFFSET];

        nb_sig_coeff += ovcabac_read_ae_sb_4x4_last_dc_c_sdh(cabac_ctx, sb_coeffs,
                                                              &c_coding_ctx);

        deq_prms.dequant_sb(sb_coeffs, deq_prms.scale, deq_prms.shift);

        memcpy(&dst[0]             , &sb_coeffs[ 0], sizeof(int16_t) * 4);
        memcpy(&dst[1 << log2_tb_w], &sb_coeffs[ 4], sizeof(int16_t) * 4);
        memcpy(&dst[2 << log2_tb_w], &sb_coeffs[ 8], sizeof(int16_t) * 4);
        memcpy(&dst[3 << log2_tb_w], &sb_coeffs[12], sizeof(int16_t) * 4);

        return 0xFFFF;

    } else if (log2_tb_h == 1) {
        int start_coeff_idx;
        int nb_sig_c;
        int sb_offset;
        int sb_pos;
        int nb_sb;
        uint8_t sig_sb_flg = 1;

        int x =  last_pos       & 0x1F;
        int y = (last_pos >> 8) & 0x1F;

        int last_sb_x = x >> 3;

        if (!last_sb_x) {
            int last_coeff_idx = x + (y << 3) ;
            int nb_coeffs = ff_vvc_diag_scan_8x2_num_cg [last_coeff_idx];
            nb_sig_c = ovcabac_read_ae_sb_8x2_dc_c_sdh(cabac_ctx, sb_coeffs,
                                                        nb_coeffs, &c_coding_ctx);

            deq_prms.dequant_sb(sb_coeffs, deq_prms.scale, deq_prms.shift);

            memcpy(&dst[0        ]     , &sb_coeffs[0], sizeof(int16_t) * 8);
            memcpy(&dst[1 << log2_tb_w], &sb_coeffs[8], sizeof(int16_t) * 8);
            /*FIXME determine whether or not we should read lfnst*/

            return 0x1;
        }

        x -= last_sb_x << 3;

        start_coeff_idx =  ff_vvc_diag_scan_8x2_num_cg [x + (y << 3)];

        sb_pos    = last_sb_x << 3;
        sb_offset = last_sb_x << 3;

        c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET + sb_offset];
        c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET + sb_offset];
        c_coding_ctx.sum_sig_nbs  = &nb_significant[VVC_TR_CTX_OFFSET + sb_offset];

        nb_sig_c = ovcabac_read_ae_sb_8x2_first_c_sdh(cabac_ctx, sb_coeffs,
                                                       start_coeff_idx, &c_coding_ctx);

        deq_prms.dequant_sb(sb_coeffs, deq_prms.scale, deq_prms.shift);

        memcpy(&dst[sb_pos + 0               ], &sb_coeffs[0], sizeof(int16_t) * 8);
        memcpy(&dst[sb_pos + (1 << log2_tb_w)], &sb_coeffs[8], sizeof(int16_t) * 8);

        nb_sb = last_sb_x;
        nb_sb--;

        for(;nb_sb > 0; --nb_sb){

            sig_sb_flg = ovcabac_read_ae_significant_sb_flag_chroma(cabac_ctx, sig_sb_flg);

            if(sig_sb_flg){

                memset(sb_coeffs, 0, sizeof(int16_t) * 16);

                sb_pos    = nb_sb << 3;
                sb_offset = nb_sb << 3;

                c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET + sb_offset];
                c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET + sb_offset];
                c_coding_ctx.sum_sig_nbs  = &nb_significant[VVC_TR_CTX_OFFSET + sb_offset];

                nb_sig_c = ovcabac_read_ae_sb_8x2_c_sdh(cabac_ctx, sb_coeffs,
                                                         &c_coding_ctx);

                deq_prms.dequant_sb(sb_coeffs, deq_prms.scale, deq_prms.shift);

                memcpy(&dst[sb_pos + 0               ], &sb_coeffs[0], sizeof(int16_t) * 8);
                memcpy(&dst[sb_pos + (1 << log2_tb_w)], &sb_coeffs[8], sizeof(int16_t) * 8);
            }
        }

        memset(sb_coeffs, 0, sizeof(uint16_t) * 16);

        c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET];
        c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET];
        c_coding_ctx.sum_sig_nbs  = &nb_significant[VVC_TR_CTX_OFFSET];

        nb_sig_c += ovcabac_read_ae_sb_8x2_last_dc_c_sdh(cabac_ctx, sb_coeffs,
                                                          &c_coding_ctx);

        deq_prms.dequant_sb(sb_coeffs, deq_prms.scale, deq_prms.shift);

        memcpy(&dst[0             ], &sb_coeffs[0], sizeof(int16_t) * 8);
        memcpy(&dst[1 << log2_tb_w], &sb_coeffs[8], sizeof(int16_t) * 8);

        return 0xFFFF;

    } else if(log2_tb_w == 1) {
        uint8_t sig_sb_flg = 1;
        int nb_sig_c;
        int sb_pos;
        int sb_offset;
        int start_coeff_idx;
        int nb_sb;

        int x =  last_pos       & 0x1F;
        int y = (last_pos >> 8) & 0x1F;
        int last_sb_y = y >> 3;

        if (!last_sb_y) {
            int last_coeff_idx = x + (y << 1);
            int nb_coeffs = ff_vvc_diag_scan_2x8_num_cg [last_coeff_idx];

            nb_sig_c = ovcabac_read_ae_sb_2x8_dc_c_sdh(cabac_ctx, sb_coeffs,
                                                        nb_coeffs, &c_coding_ctx);

            deq_prms.dequant_sb(sb_coeffs, deq_prms.scale, deq_prms.shift);

            memcpy(&dst[0], &sb_coeffs[0],  sizeof(int16_t) * 16);
            /*FIXME determine whether or not we should read lfnst*/

            return 0xFFFF;
        }

        sb_pos    =  last_sb_y << 4;
        sb_offset = (last_sb_y << 3) * VVC_TR_CTX_STRIDE;

        nb_sb = last_sb_y;

        x = x;
        y = y - (last_sb_y << 3);

        start_coeff_idx =  ff_vvc_diag_scan_2x8_num_cg [x + (y << 1)];

        c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET + sb_offset];
        c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET + sb_offset];
        c_coding_ctx.sum_sig_nbs  = &nb_significant[VVC_TR_CTX_OFFSET + sb_offset];

        nb_sig_c = ovcabac_read_ae_sb_2x8_first_c_sdh(cabac_ctx, sb_coeffs,
                                                       start_coeff_idx, &c_coding_ctx);

        deq_prms.dequant_sb(sb_coeffs, deq_prms.scale, deq_prms.shift);

        memcpy(&dst[sb_pos + 0], &sb_coeffs[0], sizeof(int16_t) * 16);

        nb_sb--;

        for( ;nb_sb > 0; --nb_sb){

            sig_sb_flg = ovcabac_read_ae_significant_sb_flag_chroma(cabac_ctx, sig_sb_flg);

            if(sig_sb_flg){
                sb_pos    = (nb_sb << 4);
                sb_offset = (nb_sb << 3) * VVC_TR_CTX_STRIDE;

                memset(sb_coeffs, 0, sizeof(int16_t) * 16);

                c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET + sb_offset];
                c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET + sb_offset];
                c_coding_ctx.sum_sig_nbs  = &nb_significant[VVC_TR_CTX_OFFSET + sb_offset];

                nb_sig_c = ovcabac_read_ae_sb_2x8_c_sdh(cabac_ctx, sb_coeffs,
                                                         &c_coding_ctx);

                deq_prms.dequant_sb(sb_coeffs, deq_prms.scale, deq_prms.shift);

                memcpy(&dst[sb_pos], &sb_coeffs[0], sizeof(int16_t) * 16);
            }
        }

        memset(sb_coeffs, 0, sizeof(uint16_t) * 16);

        c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET];
        c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET];
        c_coding_ctx.sum_sig_nbs  = &nb_significant[VVC_TR_CTX_OFFSET];

        nb_sig_c += ovcabac_read_ae_sb_2x8_last_dc_c_sdh(cabac_ctx, sb_coeffs,
                                                          &c_coding_ctx);

        deq_prms.dequant_sb(sb_coeffs, deq_prms.scale, deq_prms.shift);

        memcpy(&dst[0], &sb_coeffs[0],  sizeof(int16_t) * 16);

        return 0xFFFF;
    }
    return 0;
}

int
residual_coding_isp_h_dpq(OVCTUDec *const ctu_dec, int16_t *const dst,
                          unsigned int log2_tb_w, unsigned int log2_tb_h,
                          uint16_t last_pos)
{
    OVCABACCtx *const cabac_ctx =  ctu_dec->cabac_ctx;
    int16_t *const _dst = dst;

    int state = 0;
    int last_x, last_y;
    int last_sb_x;
    int nb_sig_c;
    int start_coeff_idx;
    int x, y;
    int sb_offset;
    int sb_pos;
    int nb_sb;
    uint8_t sig_sb_flg = 1;
    int16_t sb_coeffs[16] = {0};
    uint8_t nb_significant [VVC_TR_CTX_SIZE];
    uint8_t sum_abs_level   [VVC_TR_CTX_SIZE];
    uint8_t sum_abs_level2 [VVC_TR_CTX_SIZE];
    uint8_t log2_red_h = OVMIN(log2_tb_h, 5);
    uint8_t log2_red_w = OVMIN(log2_tb_w, 5);

    uint16_t max_nb_bins = ((1 << (log2_red_w + log2_red_h))
                             * 28) >> 4;

    int qp = ctu_dec->dequant_luma.qp;

    struct IQScale deq_prms = derive_dequant_dpq(qp, log2_tb_w, log2_tb_h);

    VVCCoeffCodingCtx c_coding_ctx = {
        .sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET],
        .sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET],
        .sum_sig_nbs  = &nb_significant[VVC_TR_CTX_OFFSET],
        .nb_remaining_bins = max_nb_bins,
        .enable_sdh = ctu_dec->enable_sdh
    };

    memset(_dst, 0, sizeof(int16_t) * (1 << (log2_tb_w + log2_tb_h)));

    if(!last_pos){
        ovcabac_read_ae_sb_dc_coeff_dpq(cabac_ctx, sb_coeffs);

        deq_prms.dequant_sb(sb_coeffs, deq_prms.scale, deq_prms.shift);

        _dst[0] = sb_coeffs[0];

        return 1;
    }

    reset_ctx_buffers(&c_coding_ctx, log2_tb_w, log2_tb_h);

    if(log2_tb_h){
        last_x =  last_pos       & 0x1F;
        last_y = (last_pos >> 8) & 0x1F;
        last_sb_x = last_x >> 3  ;

        if(!last_sb_x){
            int last_coeff_idx = last_x + (last_y << 3) ;
            int nb_coeffs = ff_vvc_diag_scan_8x2_num_cg [last_coeff_idx];
            nb_sig_c = ovcabac_read_ae_sb_8x2_dc_dpq(cabac_ctx, sb_coeffs,
                                                      &state, nb_coeffs,
                                                      &c_coding_ctx);

            deq_prms.dequant_sb(sb_coeffs, deq_prms.scale, deq_prms.shift);

            memcpy(&_dst[0             ], &sb_coeffs[0], sizeof(int16_t) * 8);
            memcpy(&_dst[1 << log2_tb_w], &sb_coeffs[8], sizeof(int16_t) * 8);

            return nb_sig_c;
        }

        x = last_x - (last_sb_x << 3) ;
        y = last_y;

        start_coeff_idx = ff_vvc_diag_scan_8x2_num_cg [x + y * (1 << 3)];

        sb_pos    = last_sb_x << 3;
        sb_offset = last_sb_x << 3;

        c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET + sb_offset];
        c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET + sb_offset];
        c_coding_ctx.sum_sig_nbs  = &nb_significant[VVC_TR_CTX_OFFSET + sb_offset];

        if(last_sb_x < 2){
            nb_sig_c = ovcabac_read_ae_sb_8x2_first_dpq(cabac_ctx, sb_coeffs,
                                                         &state, start_coeff_idx,
                                                         &c_coding_ctx);
        } else {
            nb_sig_c = ovcabac_read_ae_sb_8x2_first_far_dpq(cabac_ctx, sb_coeffs,
                                                             &state, start_coeff_idx,
                                                             &c_coding_ctx);
        }

        deq_prms.dequant_sb(sb_coeffs, deq_prms.scale, deq_prms.shift);

        memcpy(&_dst[sb_pos + 0               ], &sb_coeffs[0], sizeof(int16_t) * 8);
        memcpy(&_dst[sb_pos + (1 << log2_tb_w)], &sb_coeffs[8], sizeof(int16_t) * 8);

        nb_sb = last_sb_x;
        nb_sb--;

        for(;nb_sb > 0; --nb_sb){
            sig_sb_flg = ovcabac_read_ae_significant_sb_flag(cabac_ctx, sig_sb_flg);
            if(sig_sb_flg){

                memset(sb_coeffs, 0, sizeof(int16_t) * 16);

                sb_pos    = nb_sb << 3;
                sb_offset = nb_sb << 3;

                c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET + sb_offset];
                c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET + sb_offset];
                c_coding_ctx.sum_sig_nbs  = &nb_significant[VVC_TR_CTX_OFFSET + sb_offset];

                if(nb_sb < 2){
                    nb_sig_c = ovcabac_read_ae_sb_8x2_dpq(cabac_ctx, sb_coeffs,
                                                           &state, &c_coding_ctx);
                } else {
                    nb_sig_c = ovcabac_read_ae_sb_8x2_far_dpq(cabac_ctx, sb_coeffs,
                                                               &state, &c_coding_ctx);
                }

                deq_prms.dequant_sb(sb_coeffs, deq_prms.scale, deq_prms.shift);

                memcpy(&_dst[sb_pos + 0               ], &sb_coeffs[0],  sizeof(int16_t) * 8);
                memcpy(&_dst[sb_pos + (1 << log2_tb_w)], &sb_coeffs[8],  sizeof(int16_t) * 8);
            }
        }

        memset(sb_coeffs, 0, sizeof(uint16_t) * 16);

        c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET],
        c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET],
        c_coding_ctx.sum_sig_nbs  = &nb_significant[VVC_TR_CTX_OFFSET],

        nb_sig_c += ovcabac_read_ae_sb_8x2_last_dc_dpq(cabac_ctx,
                                                        sb_coeffs, &state,
                                                        &c_coding_ctx);

        deq_prms.dequant_sb(sb_coeffs, deq_prms.scale, deq_prms.shift);

        memcpy(&_dst[0             ], &sb_coeffs[0],  sizeof(int16_t) * 8);
        memcpy(&_dst[1 << log2_tb_w], &sb_coeffs[8],  sizeof(int16_t) * 8);

        return nb_sig_c;
    } else {
        last_x = last_pos & 0x1F;

        last_sb_x = last_x >> 4 ;

        if(!last_sb_x){
            int last_coeff_idx = last_x;
            int nb_coeffs = last_coeff_idx;

            nb_sig_c = ovcabac_read_ae_sb_1x16_dc_dpq(cabac_ctx, sb_coeffs,
                                                       &state, nb_coeffs,
                                                       &c_coding_ctx);

            deq_prms.dequant_sb(sb_coeffs, deq_prms.scale, deq_prms.shift);

            memcpy(&_dst[0] , &sb_coeffs[0], sizeof(int16_t) * 16);

            return nb_sig_c;
        }

        x = last_x - (last_sb_x << 4);

        sb_offset = last_sb_x << 4;
        sb_pos    = last_sb_x << 4;

        c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET + sb_offset];
        c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET + sb_offset];
        c_coding_ctx.sum_sig_nbs  = &nb_significant[VVC_TR_CTX_OFFSET + sb_offset];

        nb_sig_c = ovcabac_read_ae_sb_1x16_first_dpq(cabac_ctx, sb_coeffs,
                                                      &state, x, &c_coding_ctx);

        deq_prms.dequant_sb(sb_coeffs, deq_prms.scale, deq_prms.shift);

        memcpy(&_dst[sb_pos] , &sb_coeffs[0], sizeof(int16_t) * 16);

        memset(sb_coeffs, 0, sizeof(uint16_t) * 16);

        c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET],
        c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET],
        c_coding_ctx.sum_sig_nbs  = &nb_significant[VVC_TR_CTX_OFFSET],

        nb_sig_c += ovcabac_read_ae_sb_1x16_last_dc_dpq(cabac_ctx,
                                                         sb_coeffs, &state,
                                                         &c_coding_ctx);

        deq_prms.dequant_sb(sb_coeffs, deq_prms.scale, deq_prms.shift);

        memcpy(&_dst[0] , &sb_coeffs[0], sizeof(int16_t) * 16);

        return nb_sig_c;
    }
}

int
residual_coding_isp_v_dpq(OVCTUDec *const ctu_dec, int16_t *const dst,
                          unsigned int log2_tb_w, unsigned int log2_tb_h,
                          uint16_t last_pos)
{
    OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;

    int16_t *const _dst = dst;
    int state = 0;
    int last_x, last_y;
    int last_sb_y;
    int nb_sig_c;
    int start_coeff_idx;
    int x, y;
    int sb_offset;
    int sb_pos;
    int nb_sb;
    uint8_t sig_sb_flg = 1;
    int16_t sb_coeffs[16] = {0};
    uint8_t  nb_significant[VVC_TR_CTX_SIZE];
    uint8_t  sum_abs_level  [VVC_TR_CTX_SIZE];
    uint8_t sum_abs_level2 [VVC_TR_CTX_SIZE];
    uint8_t log2_red_h = OVMIN(log2_tb_h, 5);
    uint8_t log2_red_w = OVMIN(log2_tb_w, 5);
    uint16_t max_nb_bins = ((1 << (log2_red_w + log2_red_h)) * 28) >> 4;

    VVCCoeffCodingCtx c_coding_ctx = {
        .sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET],
        .sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET],
        .sum_sig_nbs  = &nb_significant[VVC_TR_CTX_OFFSET],
        .nb_remaining_bins = max_nb_bins,
        .enable_sdh = ctu_dec->enable_sdh
    };

    int qp = ctu_dec->dequant_luma.qp;

    struct IQScale deq_prms = derive_dequant_dpq(qp, log2_tb_w, log2_tb_h);

    memset(_dst, 0, sizeof(int16_t) * (1 << (log2_tb_w + log2_tb_h)));

    if(!last_pos){

        ovcabac_read_ae_sb_dc_coeff_dpq(cabac_ctx, sb_coeffs);

        deq_prms.dequant_sb(sb_coeffs, deq_prms.scale, deq_prms.shift);

        _dst[0] = sb_coeffs[0];
        return 1;
    }

    if (log2_tb_w) reset_ctx_buffers(&c_coding_ctx, log2_tb_w, log2_tb_h);
    else reset_ctx_buffers(&c_coding_ctx, log2_tb_h, log2_tb_w);

    if(log2_tb_w){
        last_x =  last_pos       & 0x1F;
        last_y = (last_pos >> 8) & 0x1F;
        last_sb_y = last_y >> 3;
        if(!last_sb_y){
            int last_coeff_idx = last_x + (last_y << 1);
            int nb_coeffs = ff_vvc_diag_scan_2x8_num_cg [last_coeff_idx];
            nb_sig_c = ovcabac_read_ae_sb_2x8_dc_dpq(cabac_ctx, sb_coeffs,
                                                      &state, nb_coeffs,
                                                      &c_coding_ctx);

            deq_prms.dequant_sb(sb_coeffs, deq_prms.scale, deq_prms.shift);

            memcpy(&_dst[0], &sb_coeffs[0],  sizeof(int16_t) * 16);

            return nb_sig_c;
        }
        sb_pos    =  last_sb_y << 4;
        sb_offset = (last_sb_y << 3) * VVC_TR_CTX_STRIDE;

        nb_sb = last_sb_y;

        x = last_x;
        y = last_y - (last_sb_y << 3);

        start_coeff_idx = ff_vvc_diag_scan_2x8_num_cg [x + (y << 1)];

        c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET + sb_offset];
        c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET + sb_offset];
        c_coding_ctx.sum_sig_nbs  = &nb_significant[VVC_TR_CTX_OFFSET + sb_offset];

        if(nb_sb < 2){
            nb_sig_c = ovcabac_read_ae_sb_2x8_first_dpq(cabac_ctx, sb_coeffs,
                                                         &state, start_coeff_idx,
                                                         &c_coding_ctx);
        } else {
            nb_sig_c = ovcabac_read_ae_sb_2x8_first_far_dpq(cabac_ctx, sb_coeffs,
                                                             &state, start_coeff_idx,
                                                             &c_coding_ctx);
        }

        deq_prms.dequant_sb(sb_coeffs, deq_prms.scale, deq_prms.shift);

        memcpy(&_dst[sb_pos + 0], &sb_coeffs[0], sizeof(int16_t) * 16);

        nb_sb--;

        for( ;nb_sb > 0; --nb_sb){
            sig_sb_flg = ovcabac_read_ae_significant_sb_flag(cabac_ctx, sig_sb_flg);
            if(sig_sb_flg){
                sb_pos    =  nb_sb << 4;
                sb_offset = (nb_sb << 3) * VVC_TR_CTX_STRIDE ;

                memset(sb_coeffs, 0, sizeof(int16_t) * 16);

                c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET + sb_offset];
                c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET + sb_offset];
                c_coding_ctx.sum_sig_nbs  = &nb_significant[VVC_TR_CTX_OFFSET + sb_offset];

                if(nb_sb < 2){
                    nb_sig_c = ovcabac_read_ae_sb_2x8_dpq(cabac_ctx, sb_coeffs,
                                                           &state, &c_coding_ctx);
                } else {
                    nb_sig_c = ovcabac_read_ae_sb_2x8_far_dpq(cabac_ctx, sb_coeffs,
                                                               &state, &c_coding_ctx);
                }

                deq_prms.dequant_sb(sb_coeffs, deq_prms.scale, deq_prms.shift);

                memcpy(&_dst[sb_pos], &sb_coeffs[0], sizeof(int16_t) * 16);
            }
        }

        memset(sb_coeffs, 0, sizeof(uint16_t) * 16);

        c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET];
        c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET];
        c_coding_ctx.sum_sig_nbs  = &nb_significant[VVC_TR_CTX_OFFSET];

        nb_sig_c += ovcabac_read_ae_sb_2x8_last_dc_dpq(cabac_ctx,
                                                        sb_coeffs, &state,
                                                        &c_coding_ctx);

        deq_prms.dequant_sb(sb_coeffs, deq_prms.scale, deq_prms.shift);

        memcpy(&_dst[0], &sb_coeffs[0],  sizeof(int16_t) * 16);

        return nb_sig_c;

    } else {
        last_y    = (last_pos >> 8) & 0x1F;
        last_sb_y = last_y >> 4;

        if(!last_sb_y){
            nb_sig_c = ovcabac_read_ae_sb_1x16_dc_dpq(cabac_ctx, sb_coeffs,
                                                       &state, last_y, &c_coding_ctx);

            deq_prms.dequant_sb(sb_coeffs, deq_prms.scale, deq_prms.shift);

            memcpy(&_dst[0],  &sb_coeffs[0],  sizeof(int16_t) * 16);

            return nb_sig_c;
        }

        y = last_y - (last_sb_y << 4);

        sb_pos    = last_sb_y << 4;
        sb_offset = last_sb_y << 4;

        c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET + sb_offset];
        c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET + sb_offset];
        c_coding_ctx.sum_sig_nbs  = &nb_significant[VVC_TR_CTX_OFFSET + sb_offset];

        nb_sig_c = ovcabac_read_ae_sb_1x16_first_dpq(cabac_ctx, sb_coeffs,
                                                      &state, y, &c_coding_ctx);

        deq_prms.dequant_sb(sb_coeffs, deq_prms.scale, deq_prms.shift);

        memcpy(&_dst[sb_pos + 0], &sb_coeffs[0],  sizeof(int16_t) * 16);

        memset(sb_coeffs, 0, sizeof(uint16_t) * 16);

        c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET];
        c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET];
        c_coding_ctx.sum_sig_nbs  = &nb_significant[VVC_TR_CTX_OFFSET];

        nb_sig_c += ovcabac_read_ae_sb_1x16_last_dc_dpq(cabac_ctx,
                                                         sb_coeffs, &state,
                                                         &c_coding_ctx);

        deq_prms.dequant_sb(sb_coeffs, deq_prms.scale, deq_prms.shift);

        memcpy(&_dst[0], &sb_coeffs[0], sizeof(int16_t) * 16);

        return nb_sig_c;
    }
}

static inline uint8_t
log2_size_2_idx(uint8_t log2_tb_w, uint8_t log2_tb_h)
{
    log2_tb_w &= 0x7;
    log2_tb_h &= 0x7;
    return (log2_tb_w << 3) | log2_tb_h;
}

static const VVCSBScanContext *ts_scan_ctx_lut[4] =
{
    &diag_2x2_scan, &diag_2x8_scan, &diag_8x2_scan, &diag_4x4_scan
};

static const VVCSBScanContext *const
select_scan_ctx(enum TBSize tb_size_idx)
{
    uint8_t idx;
    switch (tb_size_idx) {
        case TB_2x4:
        case TB_4x2:
            idx = 0;
            break;
        case TB_2x8:
        case TB_2x16:
        case TB_2x32:
        case TB_2x64:
            idx = 1;
            break;
        case TB_8x2:
        case TB_16x2:
        case TB_32x2:
        case TB_64x2:
            idx = 2;
            break;
        default:
            idx = 3;
            break;
    }
    return ts_scan_ctx_lut[idx];
}

static void store_sb_coeff(int16_t *dst, const int16_t *sb_coeffs, uint8_t log2_tb_w,
uint8_t log2_sb_w, uint8_t log2_tb_h)
{
    uint8_t cpy_s      = sizeof(*dst) << log2_sb_w;
    uint8_t dst_stride = 1 << log2_tb_w;
    uint8_t src_stride = 1 << log2_sb_w;

    if (log2_sb_w == 2) {
        memcpy(dst, sb_coeffs,  cpy_s);
        dst       += dst_stride;
        sb_coeffs += src_stride;

        memcpy(dst, sb_coeffs,  cpy_s);
        dst       += dst_stride;
        sb_coeffs += src_stride;

        memcpy(dst, sb_coeffs,  cpy_s);
        dst       += dst_stride;
        sb_coeffs += src_stride;

        memcpy(dst, sb_coeffs, cpy_s);
    } else {
        memcpy(dst, sb_coeffs,  cpy_s);
        dst       += dst_stride;
        sb_coeffs += src_stride;

        memcpy(dst, sb_coeffs,  cpy_s);

        if (log2_tb_h == 3) {
            dst       += dst_stride;
            sb_coeffs += src_stride;
            memcpy(dst, sb_coeffs,  cpy_s);

            dst       += dst_stride;
            sb_coeffs += src_stride;
            memcpy(dst, sb_coeffs,  cpy_s);

            dst       += dst_stride;
            sb_coeffs += src_stride;
            memcpy(dst, sb_coeffs,  cpy_s);

            dst       += dst_stride;
            sb_coeffs += src_stride;
            memcpy(dst, sb_coeffs,  cpy_s);

            dst       += dst_stride;
            sb_coeffs += src_stride;
            memcpy(dst, sb_coeffs,  cpy_s);

            dst       += dst_stride;
            sb_coeffs += src_stride;
            memcpy(dst, sb_coeffs,  cpy_s);
        }
    }
}

static const uint8_t *const
select_sb_scan_map_x(int8_t log2_tb_w, int8_t log2_tb_h, int8_t log2_sb_w, int8_t log2_sb_h)
{
    uint8_t idx_w = OVMAX(0, log2_tb_w - log2_sb_w);
    uint8_t idx_h = OVMAX(0, log2_tb_h - log2_sb_h);
    return ff_vvc_scan_x[idx_w][idx_h];
}

static const uint8_t *const
select_sb_scan_map_y(int8_t log2_tb_w, int8_t log2_tb_h, int8_t log2_sb_w, int8_t log2_sb_h)
{
    uint8_t idx_w = OVMAX(0, log2_tb_w - log2_sb_w);
    uint8_t idx_h = OVMAX(0, log2_tb_h - log2_sb_h);

    return ff_vvc_scan_y[idx_w][idx_h];
}

int
residual_coding_ts(OVCTUDec *const ctu_dec, int16_t *dst,
                   uint8_t log2_tb_w, uint8_t log2_tb_h)
{
    OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;
    /* FIXME smaller reset tables */
    uint8_t nb_significant[VVC_TR_CTX_SIZE]={0};
    uint8_t sign_map[VVC_TR_CTX_SIZE]={0};
    uint16_t abs_coeffs[VVC_TR_CTX_SIZE]={0};

    int16_t sb_coeffs[16] = {0};

    /* FIXME if called from chroma ? */
    int qp = ctu_dec->dequant_skip->qp;

    const struct IQScale deq_prms = derive_dequant_ts(qp, log2_tb_w, log2_tb_h);

    TSCoeffCodingCtx cctx = {
        .nb_sig_ngh = &nb_significant   [0],
        .sign_map   = &sign_map         [0],
        .abs_coeffs = &abs_coeffs  [VVC_TR_CTX_STRIDE + 0],
    };

    const enum TBSize tb_size = log2_size_2_idx(log2_tb_w, log2_tb_h);
    const VVCSBScanContext *const scan_ctx = select_scan_ctx(tb_size);

    uint8_t log2_tb_s = log2_tb_h + log2_tb_w;

    uint16_t max_nb_bins = (((1 << log2_tb_s) << 3) - (1 << log2_tb_s)) >> 2;
    uint8_t log2_sb_s = scan_ctx->log2_sb_w + scan_ctx->log2_sb_h;

    int nb_sb = (1 << log2_tb_s) >> (log2_sb_s);

    int nb_sig_c = 0;

    memset(dst, 0, sizeof(uint16_t) << log2_tb_s);

    if (nb_sb == 1) {
        int16_t *_dst = &dst[0];

        nb_sig_c += ovcabac_read_ae_sb_ts_core(cabac_ctx, sb_coeffs,
                                               &cctx,
                                               (int16_t*) &max_nb_bins,
                                               scan_ctx);

        deq_prms.dequant_sb(sb_coeffs, deq_prms.scale, deq_prms.shift);

        store_sb_coeff(_dst, sb_coeffs, log2_tb_w, scan_ctx->log2_sb_w,
                       scan_ctx->log2_sb_h);

        return 0xFFFF;

    } else {
        const uint8_t *const scan_sb_x = select_sb_scan_map_x(log2_tb_w, log2_tb_h,
                                                              scan_ctx->log2_sb_w,
                                                              scan_ctx->log2_sb_h);
        const uint8_t *const scan_sb_y = select_sb_scan_map_y(log2_tb_w, log2_tb_h,
                                                              scan_ctx->log2_sb_w,
                                                              scan_ctx->log2_sb_h);

        uint64_t sig_sb_map = 0;
        uint8_t  sig_sb_flg = 0;

        int i;

        for (i = 0; i < nb_sb - 1; ++i) {
            int x_sb = scan_sb_x[i];
            int y_sb = scan_sb_y[i];

            /* FIXME this could be simplified */
            uint8_t sig_sb_abv = (((sig_sb_map >> ((y_sb - 1) << 3)) & 0xFF) >> x_sb) & 0x1;
            uint8_t sig_sb_lft = (((sig_sb_map >> ( y_sb      << 3)) & 0xFF) >> (x_sb - 1)) & 0x1;

            int sig_sb_offset = !!sig_sb_abv + !!sig_sb_lft;


            sig_sb_flg = ovcabac_read_ae_significant_ts_sb_flag(cabac_ctx, sig_sb_offset);

            if (sig_sb_flg) {
                int16_t x_offset = x_sb << scan_ctx->log2_sb_w;
                int16_t y_offset = y_sb << scan_ctx->log2_sb_h;
                int16_t *_dst = &dst[(x_offset) + (y_offset << log2_tb_w)];

                int sb_offset = x_offset + y_offset * (VVC_TR_CTX_STRIDE);

                cctx.nb_sig_ngh = nb_significant + sb_offset;
                cctx.sign_map   = sign_map       + sb_offset;
                cctx.abs_coeffs = abs_coeffs     + sb_offset + VVC_TR_CTX_STRIDE;

                memset(sb_coeffs, 0, sizeof(int16_t) * 16);

                sig_sb_map |= 1llu << (x_sb + (y_sb << 3));

                nb_sig_c += ovcabac_read_ae_sb_ts_core(cabac_ctx, sb_coeffs,
                                                       &cctx,
                                                       (int16_t*) &max_nb_bins,
                                                       scan_ctx);

                deq_prms.dequant_sb(sb_coeffs, deq_prms.scale, deq_prms.shift);

                store_sb_coeff(_dst, sb_coeffs, log2_tb_w, scan_ctx->log2_sb_w,
                               scan_ctx->log2_sb_h);

            }
        }

        sig_sb_flg = !sig_sb_map;

        if (sig_sb_map) {
            int x_sb = scan_sb_x[i];
            int y_sb = scan_sb_y[i];
            uint8_t sig_sb_abv = (((sig_sb_map >> ((uint8_t)(y_sb - 1) << 3)) & 0xFF) >> x_sb) & 0x1;
            uint8_t sig_sb_lft = (((sig_sb_map >> ( y_sb      << 3)) & 0xFF) >> (uint8_t)(x_sb - 1)) & 0x1;
            int sig_sb_offset = !!sig_sb_abv + !!sig_sb_lft;

            sig_sb_flg = ovcabac_read_ae_significant_ts_sb_flag(cabac_ctx, sig_sb_offset);
        }

        if (sig_sb_flg) {
            int x_sb = scan_sb_x[i];
            int y_sb = scan_sb_y[i];
            int16_t x_offset = x_sb << scan_ctx->log2_sb_w;
            int16_t y_offset = y_sb << scan_ctx->log2_sb_h;
            int16_t *_dst = &dst[(x_offset) + (y_offset << log2_tb_w)];

            int sb_offset = x_offset + y_offset * (VVC_TR_CTX_STRIDE);

            cctx.nb_sig_ngh = nb_significant + sb_offset;
            cctx.sign_map   = sign_map       + sb_offset;
            cctx.abs_coeffs = abs_coeffs     + sb_offset + VVC_TR_CTX_STRIDE;

            memset(sb_coeffs, 0, sizeof(uint16_t) * 16);

            sig_sb_map |= 1llu << (x_sb + (y_sb << 3));

            nb_sig_c += ovcabac_read_ae_sb_ts_core(cabac_ctx, sb_coeffs,
                                                  &cctx,
                                                  (int16_t*) &max_nb_bins,
                                                  scan_ctx);

            deq_prms.dequant_sb(sb_coeffs, deq_prms.scale, deq_prms.shift);

            store_sb_coeff(_dst, sb_coeffs, log2_tb_w, scan_ctx->log2_sb_w,
                           scan_ctx->log2_sb_h);
        }
    }

    return 0xFFFF;
}

uint64_t
residual_coding_dpq(OVCTUDec *const ctu_dec, int16_t *const dst,
                    unsigned int log2_tb_w, unsigned int log2_tb_h,
                    uint16_t last_pos)
{
    OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;
    uint8_t tmp_red_w =   (ctu_dec->tmp_red & 0x1);
    uint8_t tmp_red_h = !!(ctu_dec->tmp_red & 0x2);
    int16_t *const _dst = dst;
    /*FIXME we can reduce LUT based on the fact it cannot be greater than 5 */
    uint8_t lim_log2_w = OVMIN(log2_tb_w - tmp_red_w, 5);
    uint8_t lim_log2_h = OVMIN(log2_tb_h - tmp_red_h, 5);

    const uint8_t *const sb_idx_2_sb_num = ff_vvc_idx_2_num[lim_log2_w - 2]
                                                           [lim_log2_h - 2];

    const uint8_t *const scan_sb_x = ff_vvc_scan_x[lim_log2_w - 2][lim_log2_h - 2];
    const uint8_t *const scan_sb_y = ff_vvc_scan_y[lim_log2_w - 2][lim_log2_h - 2];
    int x, y;

    int state = 0;

    int last_x, last_y;
    int last_sb_x, last_sb_y;
    int nb_sb, nb_sig_c;
    int d_sb;

    int start_coeff_idx;
    int sb_offset;
    uint16_t max_nb_bins = (((1 << (lim_log2_h + lim_log2_w)) << 5)
                          - ((1 << (lim_log2_h + lim_log2_w)) << 2)) >> 4;

    int16_t sb_coeffs[16] = {0};

    //TODO avoid offsets tabs
    uint8_t nb_significant[VVC_TR_CTX_SIZE];
    uint8_t sum_abs_level  [VVC_TR_CTX_SIZE];
    uint8_t sum_abs_level2 [VVC_TR_CTX_SIZE];

    VVCCoeffCodingCtx c_coding_ctx = {
        .sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET],
        .sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET],
        .sum_sig_nbs  = &nb_significant[VVC_TR_CTX_OFFSET],
        .nb_remaining_bins = max_nb_bins,
        .enable_sdh = ctu_dec->enable_sdh
    };

    uint64_t sig_sb_map = 0;
    int sb_pos;

    int qp = ctu_dec->dequant_luma.qp;

    struct IQScale deq_prms = derive_dequant_dpq(qp, log2_tb_w, log2_tb_h);

    memset(_dst, 0, sizeof(int16_t) * (1 << (log2_tb_w + log2_tb_h)));

    if (!last_pos){

        ovcabac_read_ae_sb_dc_coeff_dpq(cabac_ctx, sb_coeffs);

        deq_prms.dequant_sb(sb_coeffs, deq_prms.scale, deq_prms.shift);

        _dst[0] = sb_coeffs[0];
        memcpy(ctu_dec->lfnst_subblock, sb_coeffs, sizeof(int16_t) * 16);
        return 1;
    }

    reset_ctx_buffers(&c_coding_ctx, lim_log2_w + tmp_red_w, lim_log2_h + tmp_red_h);

    last_x =  last_pos       & 0x1F;
    last_y = (last_pos >> 8) & 0x1F;

    last_sb_x = last_x >> 2;
    last_sb_y = last_y >> 2;

    sig_sb_map |= 1llu << (last_sb_x + (last_sb_y << 3));

    if (!last_sb_x && !last_sb_y){
        int last_coeff_idx = last_x + (last_y << 2);
        int nb_coeffs = ff_vvc_diag_scan_4x4_num_cg [last_coeff_idx];

        nb_sig_c = ovcabac_read_ae_sb_4x4_first_dpq(cabac_ctx, sb_coeffs,
                                                     &state, nb_coeffs, 0,
                                                     &c_coding_ctx);

        deq_prms.dequant_sb(sb_coeffs, deq_prms.scale, deq_prms.shift);

        memcpy(&_dst[0]             , &sb_coeffs[ 0], sizeof(int16_t) * 4);
        memcpy(&_dst[1 << log2_tb_w], &sb_coeffs[ 4], sizeof(int16_t) * 4);
        memcpy(&_dst[2 << log2_tb_w], &sb_coeffs[ 8], sizeof(int16_t) * 4);
        memcpy(&_dst[3 << log2_tb_w], &sb_coeffs[12], sizeof(int16_t) * 4);

        memcpy(ctu_dec->lfnst_subblock, sb_coeffs, sizeof(int16_t) * 16);

        return sig_sb_map;
    }

    sb_pos    = (last_sb_x << 2) + ((last_sb_y << log2_tb_w) << 2);
    sb_offset = (last_sb_x << 2) + (last_sb_y << 2) * VVC_TR_CTX_STRIDE;

    d_sb = last_sb_x + last_sb_y;

    x = last_x & 0x3;
    y = last_y & 0x3;

    start_coeff_idx = ff_vvc_diag_scan_4x4_num_cg[x + (y << 2)];

    c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET + sb_offset];
    c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET + sb_offset];
    c_coding_ctx.sum_sig_nbs  = &nb_significant[VVC_TR_CTX_OFFSET + sb_offset];

    nb_sig_c = ovcabac_read_ae_sb_4x4_first_dpq(cabac_ctx, sb_coeffs,
                                                 &state, start_coeff_idx, d_sb,
                                                 &c_coding_ctx);

    deq_prms.dequant_sb(sb_coeffs, deq_prms.scale, deq_prms.shift);

    memcpy(&_dst[sb_pos + (0)]             , &sb_coeffs[ 0], sizeof(int16_t) * 4);
    memcpy(&_dst[sb_pos + (1 << log2_tb_w)], &sb_coeffs[ 4], sizeof(int16_t) * 4);
    memcpy(&_dst[sb_pos + (2 << log2_tb_w)], &sb_coeffs[ 8], sizeof(int16_t) * 4);
    memcpy(&_dst[sb_pos + (3 << log2_tb_w)], &sb_coeffs[12], sizeof(int16_t) * 4);

    nb_sb = sb_idx_2_sb_num[last_sb_x + last_sb_y * ((1 << lim_log2_w) >> 2)];

    nb_sb--;

    for(; nb_sb > 0; --nb_sb){
        int x_sb = scan_sb_x[nb_sb];
        int y_sb = scan_sb_y[nb_sb];

        uint8_t sig_sb_flg;
        uint8_t sig_sb_blw = ((uint8_t)((sig_sb_map >> ((y_sb + 1) << 3)) & 0xFF) >> x_sb) & 0x1;
        uint8_t sig_sb_rgt = ((uint8_t)((sig_sb_map >> ( y_sb      << 3)) & 0xFF) >> (x_sb + 1)) & 0x1;

        sig_sb_flg = ovcabac_read_ae_significant_sb_flag(cabac_ctx, (sig_sb_rgt | sig_sb_blw));

        if(sig_sb_flg){

            memset(sb_coeffs, 0, sizeof(int16_t) * 16);

            sb_pos    = (x_sb << 2) + ((y_sb << log2_tb_w) << 2);
            sb_offset = (x_sb << 2) + (y_sb << 2) * (VVC_TR_CTX_STRIDE);

            sig_sb_map |= 1llu << (x_sb + (y_sb << 3));

            d_sb = x_sb + y_sb;

            c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET + sb_offset];
            c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET + sb_offset];
            c_coding_ctx.sum_sig_nbs  = &nb_significant[VVC_TR_CTX_OFFSET + sb_offset];

            nb_sig_c += ovcabac_read_ae_sb_4x4_dpq(cabac_ctx, sb_coeffs,
                                                    &state, d_sb, &c_coding_ctx);

            deq_prms.dequant_sb(sb_coeffs, deq_prms.scale, deq_prms.shift);

            memcpy(&_dst[sb_pos + (0)]             , &sb_coeffs[ 0], sizeof(int16_t) * 4);
            memcpy(&_dst[sb_pos + (1 << log2_tb_w)], &sb_coeffs[ 4], sizeof(int16_t) * 4);
            memcpy(&_dst[sb_pos + (2 << log2_tb_w)], &sb_coeffs[ 8], sizeof(int16_t) * 4);
            memcpy(&_dst[sb_pos + (3 << log2_tb_w)], &sb_coeffs[12], sizeof(int16_t) * 4);
        }
    }

    memset(sb_coeffs, 0, sizeof(uint16_t) * 16);

    c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET];
    c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET];
    c_coding_ctx.sum_sig_nbs  = &nb_significant[VVC_TR_CTX_OFFSET];

    nb_sig_c += ovcabac_read_ae_sb_4x4_last_dc_dpq(cabac_ctx, sb_coeffs, &state,
                                                    &c_coding_ctx);

    deq_prms.dequant_sb(sb_coeffs, deq_prms.scale, deq_prms.shift);

    memcpy(&_dst[0]             , &sb_coeffs[ 0], sizeof(int16_t) * 4);
    memcpy(&_dst[1 << log2_tb_w], &sb_coeffs[ 4], sizeof(int16_t) * 4);
    memcpy(&_dst[2 << log2_tb_w], &sb_coeffs[ 8], sizeof(int16_t) * 4);
    memcpy(&_dst[3 << log2_tb_w], &sb_coeffs[12], sizeof(int16_t) * 4);

    return sig_sb_map | 1;
}

static void
position_cc_ctx(VVCCoeffCodingCtx *const cc_ctx, uint8_t* buff, int16_t buff_size,
                int16_t sb_offset)
{
    uint8_t *init_pos = buff + VVC_TR_CTX_OFFSET + sb_offset;
    cc_ctx->sum_abs_lvl  = init_pos;
    cc_ctx->sum_abs_lvl2 = init_pos +  buff_size;
    cc_ctx->sum_sig_nbs  = init_pos + (buff_size << 1);
}

static void
init_cc_ctx(VVCCoeffCodingCtx *const cc_ctx, uint8_t* buff,
            uint8_t log2_tb_w, uint8_t log2_tb_h)
{
    uint16_t nb_tb_smp = 1 << (log2_tb_w + log2_tb_h);
    uint16_t max_nb_bins = ((nb_tb_smp << 5) - (nb_tb_smp << 2)) >> 4;

    cc_ctx->nb_remaining_bins = max_nb_bins,

    position_cc_ctx(cc_ctx, buff, VVC_TR_CTX_SIZE, 0);
}


static int
decode_dpq_small_h_tu_c(OVCTUDec *const ctu_dec, int16_t *const dst,
                        unsigned int log2_tb_w, unsigned int log2_tb_h,
                        uint16_t last_pos)
{
    OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;
    int16_t sb_coeffs[16] = {0}; //temporary table to store coeffs in process
    int nb_sig_c;
    int16_t *const _dst = dst;
    int state = 0;

    int qp = ctu_dec->dequant_chroma->qp;
    const struct IQScale deq_prms = derive_dequant_dpq(qp, log2_tb_w, log2_tb_h);

    uint8_t buff[VVC_TR_CTX_SIZE * 3];

    uint8_t sig_sb_flg = 1;
    int16_t last_x =  last_pos       & 0x1F;
    int16_t last_y = (last_pos >> 8) & 0x1F;

    int8_t last_sb_x = last_x >> 3;

    int16_t sb_offset;

    VVCCoeffCodingCtx c_coding_ctx;

    init_cc_ctx(&c_coding_ctx, buff, log2_tb_w, log2_tb_h);
    c_coding_ctx.enable_sdh = ctu_dec->enable_sdh;

    reset_ctx_buffers(&c_coding_ctx, log2_tb_w, log2_tb_h);



    if (!last_sb_x) {
        int last_coeff_idx = last_x + (last_y << 3) ;
        int nb_coeffs = ff_vvc_diag_scan_8x2_num_cg [last_coeff_idx];
        if (log2_tb_w <= 2) {
            if (last_x >= 2) {

                sb_offset = 2;

                position_cc_ctx(&c_coding_ctx, buff, VVC_TR_CTX_SIZE, sb_offset);
                nb_coeffs -= 4;

                nb_sig_c = ovcabac_read_ae_sb_2x2_c_dpq(cabac_ctx, sb_coeffs,
                                                         &state, nb_coeffs,
                                                         &c_coding_ctx);

                deq_prms.dequant_sb(sb_coeffs, deq_prms.scale, deq_prms.shift);

                store_sb_coeff(_dst + 2, sb_coeffs, log2_tb_w, 1, 1);

                position_cc_ctx(&c_coding_ctx, buff, VVC_TR_CTX_SIZE, 0);

                nb_sig_c = ovcabac_read_ae_sb_2x2_dc_c_dpq(cabac_ctx, sb_coeffs + 4,
                                                            &state, 4,
                                                            &c_coding_ctx);

                deq_prms.dequant_sb(sb_coeffs, deq_prms.scale, deq_prms.shift);

                store_sb_coeff(_dst, sb_coeffs + 4, log2_tb_w, 1, 1);

            } else {
                nb_sig_c = ovcabac_read_ae_sb_4x2_dc_c_dpq(cabac_ctx, sb_coeffs,
                                                            &state, nb_coeffs,
                                                            &c_coding_ctx);

                deq_prms.dequant_sb(sb_coeffs, deq_prms.scale, deq_prms.shift);

                memcpy(&_dst[0        ]     , &sb_coeffs[0], sizeof(int16_t) * 2);
                memcpy(&_dst[1 << log2_tb_w], &sb_coeffs[8], sizeof(int16_t) * 2);
            }
        } else {
            nb_sig_c = ovcabac_read_ae_sb_8x2_dc_c_dpq(cabac_ctx, sb_coeffs,
                                                        &state, nb_coeffs,
                                                        &c_coding_ctx);

            deq_prms.dequant_sb(sb_coeffs, deq_prms.scale, deq_prms.shift);

            memcpy(&_dst[0        ]     , &sb_coeffs[0], sizeof(int16_t) * 8);
            memcpy(&_dst[1 << log2_tb_w], &sb_coeffs[8], sizeof(int16_t) * 8);
        }
        /*FIXME determine whether or not we should read lfnst*/

        return 0xFFFF;
    }

    int16_t x = last_x - (last_sb_x << 3);
    int16_t y = last_y;

    int16_t start_coeff_idx =  ff_vvc_diag_scan_8x2_num_cg [x + y * (1 << 3)];

    int16_t sb_pos    = last_sb_x << 3;
    sb_offset = last_sb_x << 3;

    position_cc_ctx(&c_coding_ctx, buff, VVC_TR_CTX_SIZE, sb_offset);

    nb_sig_c = ovcabac_read_ae_sb_8x2_first_c_dpq(cabac_ctx, sb_coeffs,
                                                   &state, start_coeff_idx,
                                                   &c_coding_ctx);

    deq_prms.dequant_sb(sb_coeffs, deq_prms.scale, deq_prms.shift);

    memcpy(&_dst[sb_pos + 0               ], &sb_coeffs[0], sizeof(int16_t) * 8);
    memcpy(&_dst[sb_pos + (1 << log2_tb_w)], &sb_coeffs[8], sizeof(int16_t) * 8);

    int16_t nb_sb = last_sb_x;
    nb_sb--;

    for ( ;nb_sb > 0; --nb_sb) {

        sig_sb_flg = ovcabac_read_ae_significant_sb_flag_chroma(cabac_ctx, sig_sb_flg);

        if(sig_sb_flg){

            memset(sb_coeffs, 0, sizeof(int16_t) * 16);

            sb_pos    = nb_sb << 3;
            sb_offset = nb_sb << 3;

            position_cc_ctx(&c_coding_ctx, buff, VVC_TR_CTX_SIZE, sb_offset);

            nb_sig_c = ovcabac_read_ae_sb_8x2_c_dpq(cabac_ctx, sb_coeffs, &state,
                                                     &c_coding_ctx);

            deq_prms.dequant_sb(sb_coeffs, deq_prms.scale, deq_prms.shift);

            memcpy(&_dst[sb_pos + 0               ], &sb_coeffs[0], sizeof(int16_t) * 8);
            memcpy(&_dst[sb_pos + (1 << log2_tb_w)], &sb_coeffs[8], sizeof(int16_t) * 8);
        }
    }

    memset(sb_coeffs, 0, sizeof(uint16_t) * 16);

    position_cc_ctx(&c_coding_ctx, buff, VVC_TR_CTX_SIZE, 0);

    nb_sig_c += ovcabac_read_ae_sb_8x2_last_dc_c_dpq(cabac_ctx, sb_coeffs, &state,
                                                      &c_coding_ctx);

    deq_prms.dequant_sb(sb_coeffs, deq_prms.scale, deq_prms.shift);

    memcpy(&_dst[0             ], &sb_coeffs[0], sizeof(int16_t) * 8);
    memcpy(&_dst[1 << log2_tb_w], &sb_coeffs[8], sizeof(int16_t) * 8);

    return 0xFFFF;

}

static int
decode_dpq_small_w_tu_c(OVCTUDec *const ctu_dec, int16_t *const dst,
                        unsigned int log2_tb_w, unsigned int log2_tb_h,
                        uint16_t last_pos)
{
    OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;
    int16_t sb_coeffs[16] = {0}; //temporary table to store coeffs in process
    uint8_t sig_sb_flg = 1;
    int state = 0;
    int nb_sig_c;
    int16_t *const _dst = dst;

    int qp = ctu_dec->dequant_chroma->qp;
    const struct IQScale deq_prms = derive_dequant_dpq(qp, log2_tb_w, log2_tb_h);

    /* FIXME reduce */
    uint8_t buff[VVC_TR_CTX_SIZE * 3];

    VVCCoeffCodingCtx c_coding_ctx;

    init_cc_ctx(&c_coding_ctx, buff, log2_tb_w, log2_tb_h);
    c_coding_ctx.enable_sdh = ctu_dec->enable_sdh;

    reset_ctx_buffers(&c_coding_ctx, log2_tb_w, log2_tb_h);

    int16_t last_x =  last_pos       & 0x1F;
    int16_t last_y = (last_pos >> 8) & 0x1F;

    int16_t last_sb_y = last_y >> 3;

    if(!last_sb_y){
        int last_coeff_idx = last_x + (last_y << 1);
        int nb_coeffs = ff_vvc_diag_scan_2x8_num_cg [last_coeff_idx];
        if (log2_tb_h <= 2) {
            static const uint8_t tmp_lut[4] = {
                0, 2, 1, 3
            };

            if (last_y >= 2) {

                nb_coeffs = tmp_lut[last_coeff_idx - 4];

                int16_t sb_offset = 2 * VVC_TR_CTX_STRIDE;

                position_cc_ctx(&c_coding_ctx, buff, VVC_TR_CTX_SIZE, sb_offset);

                nb_sig_c = ovcabac_read_ae_sb_2x2_c_dpq(cabac_ctx, sb_coeffs,
                                                         &state, nb_coeffs,
                                                         &c_coding_ctx);


                deq_prms.dequant_sb(sb_coeffs, deq_prms.scale, deq_prms.shift);

                store_sb_coeff(_dst + (2 << log2_tb_w), sb_coeffs, log2_tb_w, 1, 1);


                position_cc_ctx(&c_coding_ctx, buff, VVC_TR_CTX_SIZE, 0);

                nb_sig_c = ovcabac_read_ae_sb_2x2_dc_c_dpq(cabac_ctx, sb_coeffs + 4,
                                                            &state, 4,
                                                            &c_coding_ctx);

                deq_prms.dequant_sb(sb_coeffs, deq_prms.scale, deq_prms.shift);

                store_sb_coeff(_dst, sb_coeffs + 4, log2_tb_w, 1, 1);

            } else {
                nb_coeffs = tmp_lut[last_coeff_idx];
                nb_sig_c = ovcabac_read_ae_sb_2x4_dc_c_dpq(cabac_ctx, sb_coeffs,
                                                            &state, nb_coeffs,
                                                            &c_coding_ctx);

                deq_prms.dequant_sb(sb_coeffs, deq_prms.scale, deq_prms.shift);

                memcpy(&_dst[0], &sb_coeffs[0],  sizeof(int16_t) * 2);
                memcpy(&_dst[1<<log2_tb_w], &sb_coeffs[2],  sizeof(int16_t) * 2);
                memcpy(&_dst[2<<log2_tb_w], &sb_coeffs[4],  sizeof(int16_t) * 2);
                memcpy(&_dst[3<<log2_tb_w], &sb_coeffs[6],  sizeof(int16_t) * 2);
            }
        } else {
            nb_sig_c = ovcabac_read_ae_sb_2x8_dc_c_dpq(cabac_ctx, sb_coeffs,
                                                        &state, nb_coeffs,
                                                        &c_coding_ctx);
            deq_prms.dequant_sb(sb_coeffs, deq_prms.scale, deq_prms.shift);

            memcpy(&_dst[0], &sb_coeffs[0],  sizeof(int16_t) * 16);
        }

        /*FIXME determine whether or not we should read lfnst*/
        return 0xFFFF;
    }

    int16_t sb_pos    =  last_sb_y << 4;
    int16_t sb_offset = (last_sb_y << 3) * VVC_TR_CTX_STRIDE;

    int16_t nb_sb = last_sb_y;

    int16_t x = last_x;
    int16_t y = last_y - (last_sb_y << 3);

    int16_t start_coeff_idx =  ff_vvc_diag_scan_2x8_num_cg [x + (y << 1)];

    position_cc_ctx(&c_coding_ctx, buff, VVC_TR_CTX_SIZE, sb_offset);

    nb_sig_c =ovcabac_read_ae_sb_2x8_first_c_dpq(cabac_ctx, sb_coeffs,
                                                  &state, start_coeff_idx,
                                                  &c_coding_ctx);

    deq_prms.dequant_sb(sb_coeffs, deq_prms.scale, deq_prms.shift);

    memcpy(&_dst[sb_pos + 0], &sb_coeffs[0], sizeof(int16_t) * 16);

    nb_sb--;

    for ( ;nb_sb > 0; --nb_sb) {

        sig_sb_flg = ovcabac_read_ae_significant_sb_flag_chroma(cabac_ctx, sig_sb_flg);

        if(sig_sb_flg){
            sb_pos    = (nb_sb << 4);
            sb_offset = (nb_sb << 3) * VVC_TR_CTX_STRIDE;

            memset(sb_coeffs, 0, sizeof(int16_t) * 16);

            position_cc_ctx(&c_coding_ctx, buff, VVC_TR_CTX_SIZE, sb_offset);

            nb_sig_c = ovcabac_read_ae_sb_2x8_c_dpq(cabac_ctx, sb_coeffs, &state,
                                                     &c_coding_ctx);

            deq_prms.dequant_sb(sb_coeffs, deq_prms.scale, deq_prms.shift);

            memcpy(&_dst[sb_pos], &sb_coeffs[0], sizeof(int16_t) * 16);
        }
    }

    memset(sb_coeffs, 0, sizeof(uint16_t) * 16);

    position_cc_ctx(&c_coding_ctx, buff, VVC_TR_CTX_SIZE, 0);

    nb_sig_c += ovcabac_read_ae_sb_2x8_last_dc_c_dpq(cabac_ctx, sb_coeffs, &state,
                                                     &c_coding_ctx);

    deq_prms.dequant_sb(sb_coeffs, deq_prms.scale, deq_prms.shift);

    memcpy(&_dst[0], &sb_coeffs[0],  sizeof(int16_t) * 16);

    return 0xFFFF;
}

int
residual_coding_chroma_dpq(OVCTUDec *const ctu_dec, int16_t *const dst,
                           unsigned int log2_tb_w, unsigned int log2_tb_h,
                           uint16_t last_pos)
{
    OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;
    int16_t *const _dst = dst;
    //check for dependent quantization
    int state = 0;
    int16_t sb_coeffs[16] = {0}; //temporary table to store coeffs in process

    int qp = ctu_dec->dequant_chroma->qp;
    const struct IQScale deq_prms = derive_dequant_dpq(qp, log2_tb_w, log2_tb_h);

    uint8_t buff[VVC_TR_CTX_SIZE * 3];

    VVCCoeffCodingCtx c_coding_ctx;

    /* FIXME this is a bit wasteful if we only read a few sub blocks */
    memset(_dst, 0, sizeof(int16_t) * (1 << (log2_tb_w + log2_tb_h)));

    if (!last_pos){

        ovcabac_read_ae_sb_dc_coeff_c_dpq(cabac_ctx, sb_coeffs);

        deq_prms.dequant_sb(sb_coeffs, deq_prms.scale, deq_prms.shift);
        _dst[0] = sb_coeffs [0];
        return 0x1;
    }

    init_cc_ctx(&c_coding_ctx, buff, log2_tb_w, log2_tb_h);
    /* Sign data hiding is implicitly disabled when dependent quantization
     * is active.
     */
    c_coding_ctx.enable_sdh = 0;

    reset_ctx_buffers(&c_coding_ctx, log2_tb_w, log2_tb_h);

    if (log2_tb_w > 1 && log2_tb_h > 1) {
        const uint8_t log2_sb_w = 2;
        const uint8_t log2_sb_h = 2;
        int nb_sb_w  = (1 << log2_tb_w) >> log2_sb_w;
        const uint8_t *const sb_idx_2_sb_num = ff_vvc_idx_2_num[log2_tb_w - 2]
                                                               [log2_tb_h - 2];

        const uint8_t *const scan_sb_x = ff_vvc_scan_x[log2_tb_w - 2][log2_tb_h - 2];
        const uint8_t *const scan_sb_y = ff_vvc_scan_y[log2_tb_w - 2][log2_tb_h - 2];

        int nb_sig_coeff;

        int16_t last_x =  last_pos       & 0x1F;
        int16_t last_y = (last_pos >> 8) & 0x1F;

        int16_t x = last_x & ((1 << log2_sb_w) - 1);
        int16_t y = last_y & ((1 << log2_sb_h) - 1);

        int16_t sb_x = last_x >> log2_sb_w;
        int16_t sb_y = last_y >> log2_sb_h;

        int nb_c_first_sb = ff_vvc_diag_scan_4x4_num_cg[x + (y << log2_sb_h)];

        if(!sb_x && !sb_y){

            ovcabac_read_ae_sb_4x4_dc_c_dpq(cabac_ctx, sb_coeffs,
                                            &state, nb_c_first_sb, &c_coding_ctx);

            deq_prms.dequant_sb(sb_coeffs, deq_prms.scale, deq_prms.shift);

            memcpy(&_dst[0]             , &sb_coeffs[ 0], sizeof(int16_t) * 4);
            memcpy(&_dst[1 << log2_tb_w], &sb_coeffs[ 4], sizeof(int16_t) * 4);
            memcpy(&_dst[2 << log2_tb_w], &sb_coeffs[ 8], sizeof(int16_t) * 4);
            memcpy(&_dst[3 << log2_tb_w], &sb_coeffs[12], sizeof(int16_t) * 4);

            return 0x1;
        }

        uint8_t nb_sb = sb_idx_2_sb_num[sb_x + sb_y * nb_sb_w];
        uint64_t sig_sb_map = 0;

        int16_t sb_offset = (sb_x << log2_sb_w) + (sb_y << log2_sb_h) * (VVC_TR_CTX_STRIDE);
        position_cc_ctx(&c_coding_ctx, buff, VVC_TR_CTX_SIZE, sb_offset);

        nb_sig_coeff = ovcabac_read_ae_sb_4x4_first_c_dpq(cabac_ctx, sb_coeffs,
                                                           &state, nb_c_first_sb,
                                                           &c_coding_ctx);

        deq_prms.dequant_sb(sb_coeffs, deq_prms.scale, deq_prms.shift);

        int16_t sb_pos = (sb_x << log2_sb_w) + ((sb_y << log2_sb_h) * (nb_sb_w << log2_sb_w));
        memcpy(&_dst[sb_pos + (0)]             , &sb_coeffs[ 0], sizeof(int16_t) * 4);
        memcpy(&_dst[sb_pos + (1 << log2_tb_w)], &sb_coeffs[ 4], sizeof(int16_t) * 4);
        memcpy(&_dst[sb_pos + (2 << log2_tb_w)], &sb_coeffs[ 8], sizeof(int16_t) * 4);
        memcpy(&_dst[sb_pos + (3 << log2_tb_w)], &sb_coeffs[12], sizeof(int16_t) * 4);

        sig_sb_map |= 1llu << (sb_x + (sb_y << 3));

        nb_sb--;

        for(int i = nb_sb; i > 0; --i){
            int x_sb = scan_sb_x[i];
            int y_sb = scan_sb_y[i];
            uint8_t sig_sb_blw = (((sig_sb_map >> ((y_sb + 1) << 3)) & 0xFF) >> x_sb) & 0x1;
            uint8_t sig_sb_rgt = (((sig_sb_map >> ( y_sb      << 3)) & 0xFF) >> (x_sb + 1)) & 0x1;

            uint8_t sig_sb_flg = ovcabac_read_ae_significant_sb_flag_chroma(cabac_ctx, !!(sig_sb_rgt | sig_sb_blw));

            if(sig_sb_flg){

                memset(sb_coeffs, 0, sizeof(int16_t) * 16);

                sig_sb_map |= 1llu << (x_sb + (y_sb << 3));

                sb_offset = (x_sb << log2_sb_w) + (y_sb << log2_sb_h) * (VVC_TR_CTX_STRIDE);
                position_cc_ctx(&c_coding_ctx, buff, VVC_TR_CTX_SIZE, sb_offset);

                nb_sig_coeff += ovcabac_read_ae_sb_4x4_c_dpq(cabac_ctx, sb_coeffs,
                                                              &state, &c_coding_ctx);

                deq_prms.dequant_sb(sb_coeffs, deq_prms.scale, deq_prms.shift);

                sb_pos = (x_sb << log2_sb_w) + ((y_sb << log2_sb_h) * (nb_sb_w << log2_sb_w));
                memcpy(&_dst[sb_pos + (0)]             , &sb_coeffs[ 0], sizeof(int16_t) * 4);
                memcpy(&_dst[sb_pos + (1 << log2_tb_w)], &sb_coeffs[ 4], sizeof(int16_t) * 4);
                memcpy(&_dst[sb_pos + (2 << log2_tb_w)], &sb_coeffs[ 8], sizeof(int16_t) * 4);
                memcpy(&_dst[sb_pos + (3 << log2_tb_w)], &sb_coeffs[12], sizeof(int16_t) * 4);
            }
        }

        memset(sb_coeffs, 0, sizeof(int16_t) * 16);

        position_cc_ctx(&c_coding_ctx, buff, VVC_TR_CTX_SIZE, 0);

        nb_sig_coeff += ovcabac_read_ae_sb_4x4_last_dc_c_dpq(cabac_ctx, sb_coeffs,
                                                         &state, &c_coding_ctx);

        deq_prms.dequant_sb(sb_coeffs, deq_prms.scale, deq_prms.shift);

        memcpy(&_dst[0]             , &sb_coeffs[ 0], sizeof(int16_t) * 4);
        memcpy(&_dst[1 << log2_tb_w], &sb_coeffs[ 4], sizeof(int16_t) * 4);
        memcpy(&_dst[2 << log2_tb_w], &sb_coeffs[ 8], sizeof(int16_t) * 4);
        memcpy(&_dst[3 << log2_tb_w], &sb_coeffs[12], sizeof(int16_t) * 4);

       return 0xFFFF;

    } else if (log2_tb_h == 1) {
        return decode_dpq_small_h_tu_c(ctu_dec, dst, log2_tb_w, log2_tb_h, last_pos);
    } else if (log2_tb_w == 1) {
        return decode_dpq_small_w_tu_c(ctu_dec, dst, log2_tb_w, log2_tb_h, last_pos);
    }
    return 0;
}
