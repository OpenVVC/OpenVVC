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

#define VVC_TR_CTX_STRIDE (32+2)
#define VVC_TR_CTX_OFFSET ((VVC_TR_CTX_STRIDE)*2+2)
#define VVC_TR_CTX_SIZE   (VVC_TR_CTX_STRIDE*VVC_TR_CTX_STRIDE)

typedef struct VVCDepQuantCtx{
    const uint16_t state_trans_table;
    const uint8_t state_offset;
}VVCDepQuantCtx;

typedef struct VVCCoeffCodingCtx{
    uint8_t  *sum_sig_nbs;
    uint8_t  *sum_abs_lvl;
    uint8_t *sum_abs_lvl2;
    int16_t num_remaining_bins;
    uint8_t enable_sdh;
}VVCCoeffCodingCtx;

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
}VVCSBScanContext;

typedef struct VVCSBContext{
    VVCSBScanContext scan_ctx;
    VVCSBStates ctx_offsets;
    int num_rem_bins;
}VVCSBContext;


static const VVCDepQuantCtx luma_dep_quant_ctx = {
    0x7D28,
    12
};

static const VVCDepQuantCtx chroma_dep_quant_ctx = {
    0x7D28,
    8
};

static const uint64_t inv_diag_map_4x4 = 0x041852C963DA7EBF;
static const uint64_t inv_diag_map_2x8 = 0x021436587A9CBEDF;
static const uint64_t inv_diag_map_8x2 = 0x08192A3B4C5D6E7F;

static const uint64_t diag_map_4X4 = 0xFBE7AD369C258140;
static const uint64_t diag_map_2x8 = 0xFDEBC9A785634120;
static const uint64_t diag_map_8x2 = 0xF7E6D5C4B3A29180;


static const VVCSBScanContext inv_diag_4x4_scan = {
     0x041852C963DA7EBF,
     0x0259148C37BE6ADF,
     2,
};

static const VVCSBScanContext inv_diag_2x8_scan = {
     0x021436587A9CBEDF,
     0x021436587A9CBEDF,
     1,
};

static const VVCSBScanContext inv_diag_8x2_scan = {
     0x08192A3B4C5D6E7F,
     0x02468ACE13579BDF,
     3,
};

static const VVCSBScanContext inv_diag_1x16_scan = {
     0x0123456789ABCDEF,
     0x0123456789ABCDEF,
     4,
};

    static const uint64_t parity_flag_offset_map[3] ={
        0xFAAAAA5555555555,
        0x5555555555555550,
        0x5550000000000000
    };

static const uint64_t sig_flag_offset_map[3] ={
    0x8884444444444000,
    0x4000000000000000,
    0x0000000000000000
};

static const VVCSBStates luma_ctx_offsets = {
    SIG_FLAG_CTX_OFFSET,
    GT0_FLAG_CTX_OFFSET,
    PAR_FLAG_CTX_OFFSET,
    GT1_FLAG_CTX_OFFSET,
};

static const VVCSBStates chroma_ctx_offsets = {
    SIG_FLAG_C_CTX_OFFSET,
    GT0_FLAG_C_CTX_OFFSET,
    PAR_FLAG_C_CTX_OFFSET,
    GT1_FLAG_C_CTX_OFFSET,
};

static uint8_t
ovcabac_read_ae_significant_cg_flag(OVCABACCtx *const cabac_ctx,
                                    uint8_t got_significant_neighbour)
{
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    return ovcabac_ae_read(cabac_ctx, &cabac_state[SIG_COEFF_GROUP_CTX_OFFSET + got_significant_neighbour]);
}

static uint8_t
ovcabac_read_ae_significant_cg_flag_chroma(OVCABACCtx *const cabac_ctx,
                                           uint8_t got_significant_neighbour)
{
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    return ovcabac_ae_read(cabac_ctx, &cabac_state[SIG_COEFF_GROUP_C_CTX_OFFSET + got_significant_neighbour]);
}

static uint8_t
ovcabac_read_ae_significant_ts_cg_flag(OVCABACCtx *const cabac_ctx,
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

static int inline decode_truncated_rice(OVCABACCtx *const cabac_ctx,
                                                   unsigned int rice_param){

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
                  int16_t  *const coeffs, int num_pass2,
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

    for (scan_pos = 0; scan_pos < num_pass2; ++scan_pos){
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
                      int *const num_significant_coeff,
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
    uint8_t pos_shift = ((15 - (last_scan_pos & 0xF)) << 2);
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

            significant_map[(*num_significant_coeff)++] = idx;
        }
        // update state transition context
        *state = ( state_trans_tab >> (( *state << 2) + ((value & 1) << 1)) ) & 3;
        scan_map >>= 4;
    }
}

static void inline
decode_signs(OVCABACCtx *const cabac_ctx, int16_t *const coeffs,
             uint32_t state_map, const uint8_t *const sig_c_idx_map,
             int num_sig_c)
{
    uint32_t num_signs = num_sig_c;
    uint32_t signs_map = 0;

    while (num_signs--){
        signs_map <<= 1;
        signs_map  |= ovcabac_bypass_read(cabac_ctx);
    }

    signs_map <<= 32 - num_sig_c;
    state_map <<= 32 - num_sig_c;

    for (unsigned  k = 0; k < num_sig_c; k++ ){
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
    int num_sig_c = 0;
    int num_pass2 = 0;

    uint32_t dep_quant_map = 0;

    uint8_t par_lvl_flag, abs_gt1_flag, abs_gt2_flag;
    int32_t coeff_val;
    uint16_t tr_ctx_pos;

    int num_rem_bins = c_coding_ctx->num_remaining_bins;
    int prev_state  = *state;

    // Implicit first coeff
    int scan_pos = start_pos;
    uint8_t pos_shift = ((15 - (scan_pos & 0xF)) << 2);
    uint64_t scan_map = inv_diag_map        >> pos_shift;
    uint64_t par_map  = par_flag_offset_map >> pos_shift;
    uint64_t sig_map  = sig_flag_offset_map >> pos_shift;

    uint8_t idx = scan_map & 0xF;

    abs_gt1_flag = ovcabac_ae_read(cabac_ctx, abs_gt1_ctx);

    --num_rem_bins;
    coeff_val = 1 + abs_gt1_flag;

    if(abs_gt1_flag){
        par_lvl_flag = ovcabac_ae_read(cabac_ctx, par_lvl_ctx);
        abs_gt2_flag = ovcabac_ae_read(cabac_ctx, abs_gt2_ctx);
        coeff_val += par_lvl_flag;
        num_rem_bins -= 2;
        if(abs_gt2_flag){
           coeff_val += 2; //(rem_abs_gt2_flag << 1)
           gt2_idx_map[num_pass2++] = idx;
        }
    }

    coeffs[idx] = coeff_val;
    sig_idx_map[num_sig_c++] = idx;

    tr_ctx_pos = (idx & x_mask) + (idx >> log2_sb_w) * VVC_TR_CTX_STRIDE;

    set_implicit_coeff_ngbh(c_coding_ctx, tr_ctx_pos, coeff_val);

    dep_quant_map |= (prev_state >> 1);
    prev_state = (state_trans_tab >> ((prev_state << 2) + ((coeff_val & 1) << 1))) & 3;

    --scan_pos;
    scan_map >>= 4;
    par_map  >>= 4;
    sig_map  >>= 4;

    // First pass
    for( ; scan_pos >= 0 && num_rem_bins >= 4; --scan_pos ){

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
        --num_rem_bins;

        if (sig_coeff_flag){

            ctx_offset  = 1;
            ctx_offset += OVMIN(c_coding_ctx->sum_sig_nbs[tr_ctx_pos], 4);
            ctx_offset += (par_map & 0xF);

            abs_gt1_flag = ovcabac_ae_read(cabac_ctx, abs_gt1_ctx + ctx_offset);

            if (abs_gt1_flag){
                par_lvl_flag = ovcabac_ae_read(cabac_ctx, par_lvl_ctx + ctx_offset);
                abs_gt2_flag = ovcabac_ae_read(cabac_ctx, abs_gt2_ctx + ctx_offset);
                coeff_val = 2 + par_lvl_flag;
                num_rem_bins -= 2;
                if (abs_gt2_flag){
                    coeff_val += 2; /* (abs_gt2_flag << 1) */
                    gt2_idx_map[num_pass2++] = idx;
                }
            }

            dep_quant_map <<= 1;
            dep_quant_map |= (prev_state >> 1);

            --num_rem_bins;
            sig_idx_map[num_sig_c++] = idx;

            coeffs[idx] = coeff_val;
            update_coeff_nbgh_first_pass(c_coding_ctx, tr_ctx_pos, coeff_val);
        }
        prev_state = (state_trans_tab >> (( prev_state << 2) + ((coeff_val & 1) << 1))) & 3;
        scan_map >>= 4;
        par_map  >>= 4;
        sig_map  >>= 4;
    }

    if (num_pass2){
        decode_pass2_core(cabac_ctx, coeffs, num_pass2, gt2_idx_map,
                          scan_ctx, c_coding_ctx);
    }

    if (scan_pos >= 0){
        decode_bypassed_coeff_core(cabac_ctx, coeffs, scan_pos, sig_idx_map,
                                   &num_sig_c, scan_ctx, c_coding_ctx,
                                   &luma_dep_quant_ctx, &prev_state, &dep_quant_map);
    }

    decode_signs(cabac_ctx, coeffs, dep_quant_map, sig_idx_map, num_sig_c);

    c_coding_ctx->num_remaining_bins = num_rem_bins;
    *state = prev_state;

    /*FIXME we could return state instead of num_sig and avoid using pointers*/
    return num_sig_c;
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
    int num_sig_c = 0;
    int num_pass2 = 0;

    uint32_t dep_quant_map = 0;

    uint8_t par_lvl_flag, abs_gt1_flag, abs_gt2_flag;
    int32_t coeff_val;
    uint16_t tr_ctx_pos;

    int num_rem_bins = c_coding_ctx->num_remaining_bins;
    int prev_state  = *state;

    int scan_pos = 15;
    uint64_t scan_map = inv_diag_map;
    uint64_t par_map  = par_flag_offset_map;
    uint64_t sig_map  = sig_flag_offset_map;

    uint8_t idx;

    // First pass
    for ( ; scan_pos > 0 && num_rem_bins >= 4; --scan_pos ){

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
        --num_rem_bins;

        if (sig_coeff_flag){

            ctx_offset  = 1;
            ctx_offset += OVMIN(c_coding_ctx->sum_sig_nbs[tr_ctx_pos], 4);
            ctx_offset += (par_map & 0xF);

            abs_gt1_flag = ovcabac_ae_read(cabac_ctx, abs_gt1_ctx + ctx_offset);

            if (abs_gt1_flag){
                par_lvl_flag = ovcabac_ae_read(cabac_ctx, par_lvl_ctx + ctx_offset);
                abs_gt2_flag = ovcabac_ae_read(cabac_ctx, abs_gt2_ctx + ctx_offset);
                coeff_val = 2 + par_lvl_flag;
                num_rem_bins -= 2;
                if (abs_gt2_flag){
                    coeff_val += 2; /* (abs_gt2_flag << 1) */
                    gt2_idx_map[num_pass2++] = idx;
                }
            }

            dep_quant_map <<= 1;
            dep_quant_map |= (prev_state >> 1);

            --num_rem_bins;
            sig_idx_map[num_sig_c++] = idx;

            coeffs[idx] = coeff_val;
            update_coeff_nbgh_first_pass(c_coding_ctx, tr_ctx_pos, coeff_val);
        }
        prev_state = (state_trans_tab >> (( prev_state << 2) + ((coeff_val & 1) << 1))) & 3;
        scan_map >>= 4;
        par_map  >>= 4;
        sig_map  >>= 4;
    }

    if (scan_pos == 0 && num_rem_bins >= 4){

        uint8_t ctx_offset;
        uint8_t sig_coeff_flag;
        sig_coeff_flag= 1;
        idx = scan_map & 0xF;
        tr_ctx_pos = (idx & x_mask) + (idx >> log2_sb_w) * VVC_TR_CTX_STRIDE;

        //decrease scan_pos so we know last sig_coeff was read in first pass or not
        --scan_pos;

        if (num_sig_c){
            ctx_offset  = state_offset * OVMAX(0, prev_state - 1);
            ctx_offset += OVMIN(((c_coding_ctx->sum_abs_lvl[tr_ctx_pos] + 1) >> 1), 3);
            ctx_offset += sig_map & 0xF;

            sig_coeff_flag = ovcabac_ae_read(cabac_ctx, sig_flg_ctx + ctx_offset);

            --num_rem_bins;
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
                num_rem_bins -= 2;
                if (abs_gt2_flag){
                    coeff_val += 2; /* (abs_gt2_flag << 1) */
                    gt2_idx_map[num_pass2++] = idx;
                }
            }

            dep_quant_map <<= 1;
            dep_quant_map |= (prev_state >> 1);

            --num_rem_bins;
            sig_idx_map[num_sig_c++] = idx;

            coeffs[idx] = coeff_val;
            update_coeff_nbgh_first_pass(c_coding_ctx, tr_ctx_pos, coeff_val);
        }
        prev_state = (state_trans_tab >> (( prev_state << 2) + ((coeff_val & 1) << 1))) & 3;
        scan_map >>= 4;
        par_map  >>= 4;
        sig_map  >>= 4;
    }

    if (num_pass2){
        decode_pass2_core(cabac_ctx, coeffs, num_pass2, gt2_idx_map,
                          scan_ctx, c_coding_ctx);
    }

    if (scan_pos >= 0){
        decode_bypassed_coeff_core(cabac_ctx, coeffs, scan_pos, sig_idx_map,
                                   &num_sig_c, scan_ctx, c_coding_ctx,
                                   &luma_dep_quant_ctx, &prev_state, &dep_quant_map);
    }

    decode_signs(cabac_ctx, coeffs, dep_quant_map, sig_idx_map, num_sig_c);

    c_coding_ctx->num_remaining_bins = num_rem_bins;
    *state = prev_state;

    /*FIXME we could return state instead of num_sig and avoid using pointers*/
    return num_sig_c;
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
    int num_sig_c = 0;
    int num_pass2 = 0;

    uint32_t dep_quant_map = 0;

    uint8_t par_lvl_flag, abs_gt1_flag, abs_gt2_flag;
    int32_t coeff_val;
    uint16_t tr_ctx_pos;

    int num_rem_bins = c_coding_ctx->num_remaining_bins;
    int prev_state  = *state;

    int scan_pos = 15;
    uint64_t scan_map = inv_diag_map;
    uint64_t par_map  = par_flag_offset_map;
    uint64_t sig_map  = sig_flag_offset_map;

    uint8_t idx;

    // First pass
    for ( ; scan_pos > 0 && num_rem_bins >= 4; --scan_pos ){

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
        --num_rem_bins;

        if (sig_coeff_flag){

            ctx_offset  = 1;
            ctx_offset += OVMIN(c_coding_ctx->sum_sig_nbs[tr_ctx_pos], 4);
            ctx_offset += (par_map & 0xF);

            abs_gt1_flag = ovcabac_ae_read(cabac_ctx, abs_gt1_ctx + ctx_offset);

            if (abs_gt1_flag){
                par_lvl_flag = ovcabac_ae_read(cabac_ctx, par_lvl_ctx + ctx_offset);
                abs_gt2_flag = ovcabac_ae_read(cabac_ctx, abs_gt2_ctx + ctx_offset);
                coeff_val = 2 + par_lvl_flag;
                num_rem_bins -= 2;
                if (abs_gt2_flag){
                    coeff_val += 2; /* (abs_gt2_flag << 1) */
                    gt2_idx_map[num_pass2++] = idx;
                }
            }

            dep_quant_map <<= 1;
            dep_quant_map |= (prev_state >> 1);

            --num_rem_bins;
            sig_idx_map[num_sig_c++] = idx;

            coeffs[idx] = coeff_val;
            update_coeff_nbgh_first_pass(c_coding_ctx, tr_ctx_pos, coeff_val);
        }
        prev_state = (state_trans_tab >> (( prev_state << 2) + ((coeff_val & 1) << 1))) & 3;
        scan_map >>= 4;
        par_map  >>= 4;
        sig_map  >>= 4;
    }

    if (scan_pos == 0 && num_rem_bins >= 4){

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

        --num_rem_bins;

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
                num_rem_bins -= 2;
                if (abs_gt2_flag){
                    coeff_val += 2; /* (abs_gt2_flag << 1) */
                    gt2_idx_map[num_pass2++] = idx;
                }
            }

            dep_quant_map <<= 1;
            dep_quant_map |= (prev_state >> 1);

            --num_rem_bins;
            sig_idx_map[num_sig_c++] = idx;

            coeffs[idx] = coeff_val;
            update_coeff_nbgh_first_pass(c_coding_ctx, tr_ctx_pos, coeff_val);
        }
        prev_state = (state_trans_tab >> (( prev_state << 2) + ((coeff_val & 1) << 1))) & 3;
        scan_map >>= 4;
        par_map  >>= 4;
        sig_map  >>= 4;
    }

    if (num_pass2){
        decode_pass2_core(cabac_ctx, coeffs, num_pass2, gt2_idx_map,
                          scan_ctx, c_coding_ctx);
    }

    if (scan_pos >= 0){
        decode_bypassed_coeff_core(cabac_ctx, coeffs, scan_pos, sig_idx_map,
                                   &num_sig_c, scan_ctx, c_coding_ctx,
                                   &luma_dep_quant_ctx, &prev_state, &dep_quant_map);
    }

    decode_signs(cabac_ctx, coeffs, dep_quant_map, sig_idx_map, num_sig_c);

    c_coding_ctx->num_remaining_bins = num_rem_bins;
    *state = prev_state;

    /*FIXME we could return state instead of num_sig and avoid using pointers*/
    return num_sig_c;
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
                            uint8_t d_cg,
                            VVCCoeffCodingCtx *const c_coding_ctx)
{
    uint64_t par_flg_ofst_map = d_cg > 2 ? 0 : parity_flag_offset_map[d_cg];
    uint64_t sig_flg_ofst_map = d_cg > 2 ? 0 : sig_flag_offset_map[d_cg];

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
                           uint8_t d_cg,
                           VVCCoeffCodingCtx *const c_coding_ctx)
{
    uint64_t par_flg_ofst_map = d_cg > 2 ? 0 : parity_flag_offset_map[d_cg];
    uint64_t sig_flg_ofst_map = d_cg > 2 ? 0 : sig_flag_offset_map[d_cg];

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
update_ts_neighbourhood_first_pass(uint8_t   *const num_significant,
                                   uint8_t   *const sign_map,
                                   uint16_t  *const abs_coeffs,
                                   int x, int y,
                                   int value, int sign){

    num_significant[x + 1 +  y      * VVC_TR_CTX_STRIDE] += 1;
    num_significant[x     + (y + 1) * VVC_TR_CTX_STRIDE] += 1;

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
                int16_t *const num_remaining_bins,
                uint8_t log2_sb_w){
    int nb_read;
    int nb_pass3 = 0;
    uint64_t *const ts_ctx_table = &cabac_ctx->ctx_table[TS_GTX_FLAG_CTX_OFFSET];

    for (nb_read = 0; nb_read < nb_pass2 && *num_remaining_bins>= 4; ++nb_read){
        int idx = pass2_idx_map[nb_read];
        int x = idx & ((1 << log2_sb_w) - 1);
        int y = idx >> log2_sb_w;
        int val = 0;

        uint8_t ts_gt2_flag = ovcabac_ae_read(cabac_ctx, ts_ctx_table + 1);
        --(*num_remaining_bins);
        if (ts_gt2_flag) {
            uint8_t ts_gt3_flag = ovcabac_ae_read(cabac_ctx, ts_ctx_table + 2);
            --(*num_remaining_bins);
            val += 2;
            if (ts_gt3_flag) {
                uint8_t ts_gt4_flag = ovcabac_ae_read(cabac_ctx, ts_ctx_table + 3);
                --(*num_remaining_bins);
                val += 2;
                if (ts_gt4_flag) {
                    uint8_t ts_gt5_flag = ovcabac_ae_read(cabac_ctx, ts_ctx_table + 4);
                    --(*num_remaining_bins);
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
ovcabac_read_ae_sb_ts_4x4(OVCABACCtx *const cabac_ctx,
                          int16_t  *const coeffs,
                          uint8_t  *const nb_sig_ngh_map,
                          uint8_t  *const sign_table,
                          uint16_t *const abs_coeffs,
                          int16_t *const num_remaining_bins)
{
    uint64_t *const ctx_table = cabac_ctx->ctx_table;

    const uint8_t log2_sb_w = 2;
    
    uint8_t sig_c_idx_map[16];
    uint8_t pass2_idx_map[16];
    uint8_t pass3_idx_map[16];

    uint8_t nb_sig_c = 0;
    uint8_t nb_pass2 = 0;
    uint8_t nb_pass3 = 0;
    uint8_t nb_read_pass2 = 0;

    uint16_t sign_map = 0;

    uint64_t scan_map = diag_map_4X4;

    int coeff_idx;

    for (coeff_idx = 0; coeff_idx < 16 - 1 && *num_remaining_bins >= 4; ++coeff_idx) {

        int idx = scan_map & 0xF;

        int x = idx & ((1 << log2_sb_w) - 1);
        int y = idx >> log2_sb_w;

        uint8_t nb_sig_c_ngh = nb_sig_ngh_map[x + y * VVC_TR_CTX_STRIDE];

        uint8_t ts_sig_c_flag = ovcabac_ae_read(cabac_ctx, &ctx_table[TS_SIG_FLAG_CTX_OFFSET + nb_sig_c_ngh]);

        --(*num_remaining_bins);

        if (ts_sig_c_flag) {

            /* FIXME better computation of sign_offset */
            uint8_t nb_signs_ngh = sign_table[x + y * VVC_TR_CTX_STRIDE];
            int sign_offset = nb_sig_c_ngh != 2 ? nb_sig_c_ngh + nb_signs_ngh :
                (nb_signs_ngh == 2 ? 2 : nb_signs_ngh ^ 1);

            uint8_t ts_sign_flag = ovcabac_ae_read(cabac_ctx, &ctx_table[TS_RESIDUAL_SIGN_CTX_OFFSET + sign_offset]);
            uint8_t ts_gt1_flag  = ovcabac_ae_read(cabac_ctx, &ctx_table[TS_LRG1_FLAG_CTX_OFFSET + nb_sig_c_ngh]);
            int value = 1;

            sign_map |= ts_sign_flag << nb_sig_c;

            sig_c_idx_map[nb_sig_c++] = idx;

            --(*num_remaining_bins);
            --(*num_remaining_bins);

            if (ts_gt1_flag) {
                uint8_t ts_parity_flag = ovcabac_ae_read(cabac_ctx, &ctx_table[TS_PAR_FLAG_CTX_OFFSET]);

                value += 1 + ts_parity_flag;
                pass2_idx_map[nb_pass2++] = idx;

                --(*num_remaining_bins);
            }

            coeffs[idx] = value;
            update_ts_neighbourhood_first_pass(nb_sig_ngh_map, sign_table,
                                               abs_coeffs, x, y, value, ts_sign_flag);
        }
        scan_map >>= 4;
    }

    if (*num_remaining_bins >= 4) {
        int idx = scan_map & 0xF;

        int x = idx & ((1 << log2_sb_w) - 1);
        int y = idx >> log2_sb_w;

        uint8_t ts_sig_c_flag = !nb_sig_c;
        uint8_t nb_sig_c_ngh  = nb_sig_ngh_map[x + y * VVC_TR_CTX_STRIDE];

        if (nb_sig_c) {
            ts_sig_c_flag = ovcabac_ae_read(cabac_ctx, &ctx_table[TS_SIG_FLAG_CTX_OFFSET + nb_sig_c_ngh]);
            --(*num_remaining_bins);
        }

        if (ts_sig_c_flag) {
            uint8_t nb_signs_ngh = sign_table[x + y * VVC_TR_CTX_STRIDE];
            int sign_offset = nb_sig_c_ngh != 2 ? nb_sig_c_ngh + nb_signs_ngh :
                (nb_signs_ngh == 2 ? 2 : nb_signs_ngh ^ 1);

            uint8_t ts_sign_flag = ovcabac_ae_read(cabac_ctx, &ctx_table[TS_RESIDUAL_SIGN_CTX_OFFSET + sign_offset]);
            uint8_t ts_gt1_flag  = ovcabac_ae_read(cabac_ctx, &ctx_table[TS_LRG1_FLAG_CTX_OFFSET + nb_sig_c_ngh]);
            int value = 1;

            sign_map |= ts_sign_flag << nb_sig_c;

            sig_c_idx_map[nb_sig_c++] = 15;

            --(*num_remaining_bins);
            --(*num_remaining_bins);

            if (ts_gt1_flag) {
                uint8_t ts_parity_flag = ovcabac_ae_read(cabac_ctx, &ctx_table[TS_PAR_FLAG_CTX_OFFSET]);
                value += 1 + ts_parity_flag;
                pass2_idx_map[nb_pass2++] = 15;
                --(*num_remaining_bins);
            }

            update_ts_neighbourhood_first_pass(nb_sig_ngh_map, sign_table, abs_coeffs,
                                               x, y, value, ts_sign_flag);

            coeffs[15] = value;
        }
        ++coeff_idx;
        scan_map >>= 4;
    }

    if (nb_pass2) {
        nb_read_pass2 = decode_pass2_ts(cabac_ctx, coeffs, nb_pass2, pass2_idx_map,
                                        &nb_pass3, pass3_idx_map, abs_coeffs, num_remaining_bins,
                                        log2_sb_w);
    }

    /* Loop over already read pass2 requiring pass 3 */
    for (int i = 0; i < nb_pass3; i++) {
        int idx = pass3_idx_map[i];
        int remainder = decode_truncated_rice(cabac_ctx, 1);
        int x = idx & ((1 << log2_sb_w) - 1);
        int y = idx >> log2_sb_w;

        coeffs[idx] += remainder;
        abs_coeffs[x + y * VVC_TR_CTX_STRIDE] = coeffs[idx];
    }

    /* Loop over non already processed in pass 2 => bypass coded from pass 2 */
    for (int i = nb_read_pass2; i < nb_pass2; i++) {
        int idx = pass2_idx_map[i];
        int remainder = decode_truncated_rice(cabac_ctx, 1);
        int x = idx & ((1 << log2_sb_w) - 1);
        int y = idx >> log2_sb_w;

        coeffs[idx] += remainder;
        abs_coeffs[x + y * VVC_TR_CTX_STRIDE] = coeffs[idx];
    }

    /* FIXME we could probably move this to an other pass
     * based on the fact coeffs[idx] == 1 can only be true in first pass
     * when no greater flag is encountered
     */
    /* TS coeff Prediction */
    for (int i = 0; i < nb_sig_c; i++) {
        int idx = sig_c_idx_map[i];
        int x = idx & ((1 << log2_sb_w) - 1);
        int y = idx >> log2_sb_w;
        int max_abs_ngh = OVMAX(abs_coeffs[x     + ((y - 1) * VVC_TR_CTX_STRIDE)],
                                abs_coeffs[x - 1 + ( y      * VVC_TR_CTX_STRIDE)]);

        if (coeffs[idx] == 1 && max_abs_ngh) {
            coeffs[idx] = max_abs_ngh;
            abs_coeffs[x + y * VVC_TR_CTX_STRIDE] = coeffs[idx];
        } else {
            coeffs[idx] = coeffs[idx] - (coeffs[idx] <= max_abs_ngh);
            abs_coeffs[x + y * VVC_TR_CTX_STRIDE] = coeffs[idx];
        }
    }

    /* Bypass coded coefficients from pass 1 */
    for (int scan_idx = coeff_idx; scan_idx <= 15; scan_idx++) {
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
                          int *const num_sig_c,
                          const VVCSBScanContext *const scan_ctx,
                          const VVCCoeffCodingCtx *const c_coding_ctx)
{
    int scan_pos;
    int rice_param;

    const uint64_t inv_diag_map = scan_ctx->scan_map;
    const uint8_t log2_sb_w     = scan_ctx->log2_sb_w;
    const uint8_t x_mask = (1 << log2_sb_w) - 1;

    uint8_t pos_shift = ((15 - (last_scan_pos & 0xF)) << 2);
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

            sig_idx_map[(*num_sig_c)++] = idx;
        }
        scan_map >>= 4;
    }
}

static void inline
decode_signs_sdh(OVCABACCtx *const cabac_ctx, int16_t *const coeffs,
                 uint8_t *const sig_idx_map, int num_sig_c,
                 uint8_t use_sdh)
{
    /*FIXME we could avoid sum_abs by xor on parity_flags */
    const uint32_t num_signs = num_sig_c - use_sdh;
    uint32_t signs_map  = 0;
    uint32_t sum_parity = 0;
    uint32_t num_bins = num_signs;

    while (num_bins--){
        signs_map  = signs_map << 1;
        signs_map |= ovcabac_bypass_read(cabac_ctx);
    }

    signs_map  <<= 32 - num_signs;

    for (unsigned k = 0; k < num_signs; k++){
        int idx = sig_idx_map[k];
        int16_t abs_coeff = coeffs[idx];
        uint8_t sign = !!(signs_map & (1u << 31));
        sum_parity ^= abs_coeff;
        coeffs[idx] = sign ? -abs_coeff : abs_coeff ;
        signs_map <<= 1;
    }

    if (use_sdh){
        int idx = sig_idx_map[num_signs];
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
    int num_sig_c = 0;
    int num_pass2 = 0;

    uint8_t par_lvl_flag, abs_gt1_flag, abs_gt2_flag;
    int32_t coeff_val;
    uint16_t tr_ctx_pos;

    int num_rem_bins = c_coding_ctx->num_remaining_bins;

    // Implicit first coeff
    int scan_pos = start_pos;
    uint8_t pos_shift = ((15 - (scan_pos & 0xF)) << 2);
    uint64_t scan_map = inv_diag_map        >> pos_shift;
    uint64_t par_map  = par_flag_offset_map >> pos_shift;
    uint64_t sig_map  = sig_flag_offset_map >> pos_shift;

    uint8_t idx = scan_map & 0xF;

    abs_gt1_flag = ovcabac_ae_read(cabac_ctx, abs_gt1_ctx);

    --num_rem_bins;
    coeff_val = 1 + abs_gt1_flag;

    if(abs_gt1_flag){
        par_lvl_flag = ovcabac_ae_read(cabac_ctx, par_lvl_ctx);
        abs_gt2_flag = ovcabac_ae_read(cabac_ctx, abs_gt2_ctx);
        coeff_val += par_lvl_flag;
        num_rem_bins -= 2;
        if(abs_gt2_flag){
           coeff_val += 2; //(rem_abs_gt2_flag << 1)
           gt2_idx_map[num_pass2++] = idx;
        }
    }

    coeffs[idx] = coeff_val;
    sig_idx_map[num_sig_c++] = idx;

    tr_ctx_pos = (idx & x_mask) + (idx >> log2_sb_w) * VVC_TR_CTX_STRIDE;

    set_implicit_coeff_ngbh(c_coding_ctx, tr_ctx_pos, coeff_val);

    --scan_pos;
    scan_map >>= 4;
    par_map  >>= 4;
    sig_map  >>= 4;

    // First pass
    for( ; scan_pos >= 0 && num_rem_bins >= 4; --scan_pos ){

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
        --num_rem_bins;

        if (sig_coeff_flag){

            ctx_offset  = 1;
            ctx_offset += OVMIN(c_coding_ctx->sum_sig_nbs[tr_ctx_pos], 4);
            ctx_offset += (par_map & 0xF);

            abs_gt1_flag = ovcabac_ae_read(cabac_ctx, abs_gt1_ctx + ctx_offset);

            if (abs_gt1_flag){
                par_lvl_flag = ovcabac_ae_read(cabac_ctx, par_lvl_ctx + ctx_offset);
                abs_gt2_flag = ovcabac_ae_read(cabac_ctx, abs_gt2_ctx + ctx_offset);
                coeff_val = 2 + par_lvl_flag;
                num_rem_bins -= 2;
                if (abs_gt2_flag){
                    coeff_val += 2; /* (abs_gt2_flag << 1) */
                    gt2_idx_map[num_pass2++] = idx;
                }
            }

            --num_rem_bins;
            sig_idx_map[num_sig_c++] = idx;

            coeffs[idx] = coeff_val;
            update_coeff_nbgh_first_pass(c_coding_ctx, tr_ctx_pos, coeff_val);
        }
        scan_map >>= 4;
        par_map  >>= 4;
        sig_map  >>= 4;
    }

    if (num_pass2){
        decode_pass2_core(cabac_ctx, coeffs, num_pass2, gt2_idx_map,
                          scan_ctx, c_coding_ctx);
    }

    if (scan_pos >= 0){
        decode_bypassed_coeff_sdh(cabac_ctx, coeffs, scan_pos, sig_idx_map,
                                   &num_sig_c, scan_ctx, c_coding_ctx);
    }

    if (num_sig_c){
        uint8_t first_nz = (scan_ctx->scan_idx_map >> ((15 - (sig_idx_map[0])) << 2)) & 0XF;
        uint8_t last_nz  = (scan_ctx->scan_idx_map >> ((15 - (sig_idx_map[num_sig_c - 1])) << 2)) & 0XF;
        uint8_t use_sdh =  c_coding_ctx->enable_sdh && (first_nz - last_nz) >= 4;
        decode_signs_sdh(cabac_ctx, coeffs, sig_idx_map, num_sig_c, use_sdh);
    }

    c_coding_ctx->num_remaining_bins = num_rem_bins;

    /*FIXME we could return state instead of num_sig and avoid using pointers*/
    return num_sig_c;
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
    int num_sig_c = 0;
    int num_pass2 = 0;

    uint8_t par_lvl_flag, abs_gt1_flag, abs_gt2_flag;
    int32_t coeff_val;
    uint16_t tr_ctx_pos;

    int num_rem_bins = c_coding_ctx->num_remaining_bins;

    int scan_pos = 15;
    uint64_t scan_map = inv_diag_map;
    uint64_t par_map  = par_flag_offset_map;
    uint64_t sig_map  = sig_flag_offset_map;

    uint8_t idx;

    // First pass
    for ( ; scan_pos > 0 && num_rem_bins >= 4; --scan_pos ){

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
        --num_rem_bins;

        if (sig_coeff_flag){

            ctx_offset  = 1;
            ctx_offset += OVMIN(c_coding_ctx->sum_sig_nbs[tr_ctx_pos], 4);
            ctx_offset += (par_map & 0xF);

            abs_gt1_flag = ovcabac_ae_read(cabac_ctx, abs_gt1_ctx + ctx_offset);

            if (abs_gt1_flag){
                par_lvl_flag = ovcabac_ae_read(cabac_ctx, par_lvl_ctx + ctx_offset);
                abs_gt2_flag = ovcabac_ae_read(cabac_ctx, abs_gt2_ctx + ctx_offset);
                coeff_val = 2 + par_lvl_flag;
                num_rem_bins -= 2;
                if (abs_gt2_flag){
                    coeff_val += 2; /* (abs_gt2_flag << 1) */
                    gt2_idx_map[num_pass2++] = idx;
                }
            }

            --num_rem_bins;
            sig_idx_map[num_sig_c++] = idx;

            coeffs[idx] = coeff_val;
            update_coeff_nbgh_first_pass(c_coding_ctx, tr_ctx_pos, coeff_val);
        }
        scan_map >>= 4;
        par_map  >>= 4;
        sig_map  >>= 4;
    }

    if (scan_pos == 0 && num_rem_bins >= 4){

        uint8_t ctx_offset;
        uint8_t sig_coeff_flag;
        sig_coeff_flag= 1;
        idx = scan_map & 0xF;
        tr_ctx_pos = (idx & x_mask) + (idx >> log2_sb_w) * VVC_TR_CTX_STRIDE;

        //decrease scan_pos so we know last sig_coeff was read in first pass or not
        --scan_pos;

        if (num_sig_c){
            ctx_offset  = OVMIN(((c_coding_ctx->sum_abs_lvl[tr_ctx_pos] + 1) >> 1), 3);
            ctx_offset += sig_map & 0xF;

            sig_coeff_flag = ovcabac_ae_read(cabac_ctx, sig_flg_ctx + ctx_offset);

            --num_rem_bins;
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
                num_rem_bins -= 2;
                if (abs_gt2_flag){
                    coeff_val += 2; /* (abs_gt2_flag << 1) */
                    gt2_idx_map[num_pass2++] = idx;
                }
            }

            --num_rem_bins;
            sig_idx_map[num_sig_c++] = idx;

            coeffs[idx] = coeff_val;
            update_coeff_nbgh_first_pass(c_coding_ctx, tr_ctx_pos, coeff_val);
        }
        scan_map >>= 4;
        par_map  >>= 4;
        sig_map  >>= 4;
    }

    if (num_pass2){
        decode_pass2_core(cabac_ctx, coeffs, num_pass2, gt2_idx_map,
                          scan_ctx, c_coding_ctx);
    }

    if (scan_pos >= 0){
        decode_bypassed_coeff_sdh(cabac_ctx, coeffs, scan_pos, sig_idx_map,
                                  &num_sig_c, scan_ctx, c_coding_ctx);
    }

    if (num_sig_c){
        uint8_t first_nz = (scan_ctx->scan_idx_map >> ((15 - (sig_idx_map[0])) << 2)) & 0XF;
        uint8_t last_nz  = (scan_ctx->scan_idx_map >> ((15 - (sig_idx_map[num_sig_c - 1])) << 2)) & 0XF;
        uint8_t use_sdh =  c_coding_ctx->enable_sdh && (first_nz - last_nz) >= 4;
        decode_signs_sdh(cabac_ctx, coeffs, sig_idx_map, num_sig_c, use_sdh);
    }

    c_coding_ctx->num_remaining_bins = num_rem_bins;

    /*FIXME we could return state instead of num_sig and avoid using pointers*/
    return num_sig_c;
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
    int num_sig_c = 0;
    int num_pass2 = 0;

    uint8_t par_lvl_flag, abs_gt1_flag, abs_gt2_flag;
    int32_t coeff_val;
    uint16_t tr_ctx_pos;

    int num_rem_bins = c_coding_ctx->num_remaining_bins;

    int scan_pos = 15;
    uint64_t scan_map = inv_diag_map;
    uint64_t par_map  = par_flag_offset_map;
    uint64_t sig_map  = sig_flag_offset_map;

    uint8_t idx;

    // First pass
    for ( ; scan_pos > 0 && num_rem_bins >= 4; --scan_pos ){

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
        --num_rem_bins;

        if (sig_coeff_flag){

            ctx_offset  = 1;
            ctx_offset += OVMIN(c_coding_ctx->sum_sig_nbs[tr_ctx_pos], 4);
            ctx_offset += (par_map & 0xF);

            abs_gt1_flag = ovcabac_ae_read(cabac_ctx, abs_gt1_ctx + ctx_offset);

            if (abs_gt1_flag){
                par_lvl_flag = ovcabac_ae_read(cabac_ctx, par_lvl_ctx + ctx_offset);
                abs_gt2_flag = ovcabac_ae_read(cabac_ctx, abs_gt2_ctx + ctx_offset);
                coeff_val = 2 + par_lvl_flag;
                num_rem_bins -= 2;
                if (abs_gt2_flag){
                    coeff_val += 2; /* (abs_gt2_flag << 1) */
                    gt2_idx_map[num_pass2++] = idx;
                }
            }

            --num_rem_bins;
            sig_idx_map[num_sig_c++] = idx;

            coeffs[idx] = coeff_val;
            update_coeff_nbgh_first_pass(c_coding_ctx, tr_ctx_pos, coeff_val);
        }
        scan_map >>= 4;
        par_map  >>= 4;
        sig_map  >>= 4;
    }

    if (scan_pos == 0 && num_rem_bins >= 4){

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

        --num_rem_bins;

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
                num_rem_bins -= 2;
                if (abs_gt2_flag){
                    coeff_val += 2; /* (abs_gt2_flag << 1) */
                    gt2_idx_map[num_pass2++] = idx;
                }
            }

            --num_rem_bins;
            sig_idx_map[num_sig_c++] = idx;

            coeffs[idx] = coeff_val;
            update_coeff_nbgh_first_pass(c_coding_ctx, tr_ctx_pos, coeff_val);
        }
        scan_map >>= 4;
        par_map  >>= 4;
        sig_map  >>= 4;
    }

    if (num_pass2){
        decode_pass2_core(cabac_ctx, coeffs, num_pass2, gt2_idx_map,
                          scan_ctx, c_coding_ctx);
    }

    if (scan_pos >= 0){
        decode_bypassed_coeff_sdh(cabac_ctx, coeffs, scan_pos, sig_idx_map,
                                   &num_sig_c, scan_ctx, c_coding_ctx);
    }

    if (num_sig_c){
        uint8_t first_nz = (scan_ctx->scan_idx_map >> ((15 - (sig_idx_map[0])) << 2)) & 0XF;
        uint8_t last_nz  = (scan_ctx->scan_idx_map >> ((15 - (sig_idx_map[num_sig_c - 1])) << 2)) & 0XF;
        uint8_t use_sdh =  c_coding_ctx->enable_sdh && (first_nz - last_nz) >= 4;
        decode_signs_sdh(cabac_ctx, coeffs, sig_idx_map, num_sig_c, use_sdh);
    }

    c_coding_ctx->num_remaining_bins = num_rem_bins;

    /*FIXME we could return state instead of num_sig and avoid using pointers*/
    return num_sig_c;
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
                                 uint8_t d_cg,
                                 VVCCoeffCodingCtx *const c_coding_ctx)
{
    uint64_t par_flg_ofst_map = d_cg > 2 ? 0 : parity_flag_offset_map[d_cg];
    uint64_t sig_flg_ofst_map = d_cg > 2 ? 0 : sig_flag_offset_map[d_cg];

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
                           uint8_t d_cg,
                           VVCCoeffCodingCtx *const c_coding_ctx)
{
    uint64_t par_flg_ofst_map = d_cg > 2 ? 0 : parity_flag_offset_map[d_cg];
    uint64_t sig_flg_ofst_map = d_cg > 2 ? 0 : sig_flag_offset_map[d_cg];

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
    void (*dequant_sb)(int16_t *const cg_coeffs, int scale, int shift);
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
    int last_cg_x;
    int num_sig_c;
    int start_coeff_idx;
    int x, y;
    int cg_offset;
    int sb_pos;
    int num_cg;
    uint8_t sig_sb_flg = 1;
    int16_t cg_coeffs[16] = {0};
    uint8_t num_significant [VVC_TR_CTX_SIZE];
    uint8_t sum_abs_level   [VVC_TR_CTX_SIZE];
    uint8_t sum_abs_level2 [VVC_TR_CTX_SIZE];

    uint16_t max_num_bins = ((1 << (log2_tb_w + log2_tb_h))
                             * 28) >> 4;

    VVCCoeffCodingCtx c_coding_ctx = {
        .sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET],
        .sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET],
        .sum_sig_nbs  = &num_significant[VVC_TR_CTX_OFFSET],
        .num_remaining_bins = max_num_bins,
        .enable_sdh = ctu_dec->enable_sdh
    };

    int qp = ctu_dec->dequant_luma.qp;

    struct IQScale deq_prms = derive_dequant_sdh(qp, log2_tb_w, log2_tb_h);

    memset(_dst, 0, sizeof(int16_t) * (1 << (log2_tb_w + log2_tb_h)));

    if (!last_pos) {
        ovcabac_read_ae_sb_dc_coeff_sdh(cabac_ctx, cg_coeffs);

        deq_prms.dequant_sb(cg_coeffs, deq_prms.scale, deq_prms.shift);

        _dst[0] = cg_coeffs[0];

        return 1;
    }

    reset_ctx_buffers(&c_coding_ctx, log2_tb_w, log2_tb_h);

    if(log2_tb_h){
        last_x =  last_pos       & 0x1F;
        last_y = (last_pos >> 8) & 0x1F;
        last_cg_x = last_x >> 3  ;

        if(!last_cg_x){
            int last_coeff_idx = last_x + (last_y << 3) ;
            int num_coeffs = ff_vvc_diag_scan_8x2_num_cg [last_coeff_idx];
            num_sig_c = ovcabac_read_ae_sb_8x2_dc_sdh(cabac_ctx, cg_coeffs,
                                                      num_coeffs, &c_coding_ctx);

            deq_prms.dequant_sb(cg_coeffs, deq_prms.scale, deq_prms.shift);

            memcpy(&_dst[0             ], &cg_coeffs[0], sizeof(int16_t) * 8);
            memcpy(&_dst[1 << log2_tb_w], &cg_coeffs[8], sizeof(int16_t) * 8);

            return num_sig_c;
        }

        x = last_x - (last_cg_x << 3) ;
        y = last_y;

        start_coeff_idx = ff_vvc_diag_scan_8x2_num_cg [x + y * (1 << 3)];

        sb_pos    = last_cg_x << 3;
        cg_offset = last_cg_x << 3;

        c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET + cg_offset];
        c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET + cg_offset];
        c_coding_ctx.sum_sig_nbs  = &num_significant[VVC_TR_CTX_OFFSET + cg_offset];

        if(last_cg_x < 2){
            num_sig_c = ovcabac_read_ae_sb_8x2_first_sdh(cabac_ctx, cg_coeffs,
                                                         start_coeff_idx, &c_coding_ctx);
        } else {
            num_sig_c = ovcabac_read_ae_sb_8x2_first_far_sdh(cabac_ctx, cg_coeffs,
                                                             start_coeff_idx, &c_coding_ctx);
        }

        deq_prms.dequant_sb(cg_coeffs, deq_prms.scale, deq_prms.shift);

        memcpy(&_dst[sb_pos + 0               ], &cg_coeffs[0], sizeof(int16_t) * 8);
        memcpy(&_dst[sb_pos + (1 << log2_tb_w)], &cg_coeffs[8], sizeof(int16_t) * 8);

        num_cg = last_cg_x;
        num_cg--;

        for(;num_cg > 0; --num_cg){
            sig_sb_flg = ovcabac_read_ae_significant_cg_flag(cabac_ctx, sig_sb_flg);
            if(sig_sb_flg){

                memset(cg_coeffs, 0, sizeof(int16_t) * 16);

                sb_pos    = num_cg << 3;
                cg_offset = num_cg << 3;

                c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET + cg_offset];
                c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET + cg_offset];
                c_coding_ctx.sum_sig_nbs  = &num_significant[VVC_TR_CTX_OFFSET + cg_offset];

                if(num_cg < 2){
                    num_sig_c = ovcabac_read_ae_sb_8x2_sdh(cabac_ctx, cg_coeffs,
                                                           &c_coding_ctx);
                } else {
                    num_sig_c = ovcabac_read_ae_sb_8x2_far_sdh(cabac_ctx, cg_coeffs,
                                                               &c_coding_ctx);
                }

                deq_prms.dequant_sb(cg_coeffs, deq_prms.scale, deq_prms.shift);

                memcpy(&_dst[sb_pos + 0               ], &cg_coeffs[0],  sizeof(int16_t) * 8);
                memcpy(&_dst[sb_pos + (1 << log2_tb_w)], &cg_coeffs[8],  sizeof(int16_t) * 8);
            }
        }

        memset(cg_coeffs, 0, sizeof(uint16_t) * 16);

        c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET],
        c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET],
        c_coding_ctx.sum_sig_nbs  = &num_significant[VVC_TR_CTX_OFFSET],

        num_sig_c += ovcabac_read_ae_sb_8x2_last_dc_sdh(cabac_ctx, cg_coeffs,
                                                        &c_coding_ctx);

        deq_prms.dequant_sb(cg_coeffs, deq_prms.scale, deq_prms.shift);

        memcpy(&_dst[0             ], &cg_coeffs[0],  sizeof(int16_t) * 8);
        memcpy(&_dst[1 << log2_tb_w], &cg_coeffs[8],  sizeof(int16_t) * 8);

        return num_sig_c;
    } else {
        last_x = last_pos & 0x1F;

        last_cg_x = last_x >> 4 ;

        if(!last_cg_x){
            int last_coeff_idx = last_x;
            int num_coeffs = last_coeff_idx;

            num_sig_c = ovcabac_read_ae_sb_1x16_dc_sdh(cabac_ctx, cg_coeffs,
                                                       num_coeffs, &c_coding_ctx);

            deq_prms.dequant_sb(cg_coeffs, deq_prms.scale, deq_prms.shift);

            memcpy(&_dst[0] , &cg_coeffs[0], sizeof(int16_t) * 16);

            return num_sig_c;
        }

        x = last_x - (last_cg_x << 4);

        cg_offset = last_cg_x << 4;
        sb_pos    = last_cg_x << 4;

        c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET + cg_offset];
        c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET + cg_offset];
        c_coding_ctx.sum_sig_nbs  = &num_significant[VVC_TR_CTX_OFFSET + cg_offset];

        num_sig_c = ovcabac_read_ae_sb_1x16_first_sdh(cabac_ctx, cg_coeffs,
                                                      x, &c_coding_ctx);

        deq_prms.dequant_sb(cg_coeffs, deq_prms.scale, deq_prms.shift);

        memcpy(&_dst[sb_pos] , &cg_coeffs[0], sizeof(int16_t) * 16);

        memset(cg_coeffs, 0, sizeof(uint16_t) * 16);

        c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET],
        c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET],
        c_coding_ctx.sum_sig_nbs  = &num_significant[VVC_TR_CTX_OFFSET],

        num_sig_c += ovcabac_read_ae_sb_1x16_last_dc_sdh(cabac_ctx, cg_coeffs,
                                                         &c_coding_ctx);

        deq_prms.dequant_sb(cg_coeffs, deq_prms.scale, deq_prms.shift);

        memcpy(&_dst[0] , &cg_coeffs[0], sizeof(int16_t) * 16);

        return num_sig_c;
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
    int last_cg_y;
    int num_sig_c;
    int start_coeff_idx;
    int x, y;
    int cg_offset;
    int sb_pos;
    int num_cg;
    uint8_t sig_sb_flg = 1;
    int16_t cg_coeffs[16] = {0};
    uint8_t  num_significant[VVC_TR_CTX_SIZE];
    uint8_t  sum_abs_level  [VVC_TR_CTX_SIZE];
    uint8_t sum_abs_level2 [VVC_TR_CTX_SIZE];
    uint16_t max_num_bins = ((1 << (log2_tb_w + log2_tb_h)) * 28) >> 4;

    VVCCoeffCodingCtx c_coding_ctx = {
        .sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET],
        .sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET],
        .sum_sig_nbs  = &num_significant[VVC_TR_CTX_OFFSET],
        .num_remaining_bins = max_num_bins,
        .enable_sdh = ctu_dec->enable_sdh
    };

    int qp = ctu_dec->dequant_luma.qp;

    struct IQScale deq_prms = derive_dequant_sdh(qp, log2_tb_w, log2_tb_h);

    memset(_dst, 0, sizeof(int16_t) * (1 << (log2_tb_w + log2_tb_h)));

    if (!last_pos) {

        ovcabac_read_ae_sb_dc_coeff_sdh(cabac_ctx, cg_coeffs);

        deq_prms.dequant_sb(cg_coeffs, deq_prms.scale, deq_prms.shift);

        _dst[0] = cg_coeffs[0];
        return 1;
    }

    if (log2_tb_w) reset_ctx_buffers(&c_coding_ctx, log2_tb_w, log2_tb_h);
    else reset_ctx_buffers(&c_coding_ctx, log2_tb_h, log2_tb_w);

    if (log2_tb_w) {
        last_x =  last_pos       & 0x1F;
        last_y = (last_pos >> 8) & 0x1F;
        last_cg_y = last_y >> 3;

        if (!last_cg_y) {
            int last_coeff_idx = last_x + (last_y << 1);
            int num_coeffs = ff_vvc_diag_scan_2x8_num_cg [last_coeff_idx];
            num_sig_c = ovcabac_read_ae_sb_2x8_dc_sdh(cabac_ctx, cg_coeffs,
                                                      num_coeffs, &c_coding_ctx);

            deq_prms.dequant_sb(cg_coeffs, deq_prms.scale, deq_prms.shift);

            memcpy(&_dst[0], &cg_coeffs[0],  sizeof(int16_t) * 16);

            return num_sig_c;
        }

        sb_pos    =  last_cg_y << 4;
        cg_offset = (last_cg_y << 3) * VVC_TR_CTX_STRIDE;

        num_cg = last_cg_y;

        x = last_x;
        y = last_y - (last_cg_y << 3);

        start_coeff_idx = ff_vvc_diag_scan_2x8_num_cg [x + (y << 1)];

        c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET + cg_offset];
        c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET + cg_offset];
        c_coding_ctx.sum_sig_nbs  = &num_significant[VVC_TR_CTX_OFFSET + cg_offset];

        if(num_cg < 2){
            num_sig_c = ovcabac_read_ae_sb_2x8_first_sdh(cabac_ctx, cg_coeffs,
                                                         start_coeff_idx, &c_coding_ctx);
        } else {
            num_sig_c = ovcabac_read_ae_sb_2x8_first_far_sdh(cabac_ctx, cg_coeffs,
                                                             start_coeff_idx, &c_coding_ctx);
        }

        deq_prms.dequant_sb(cg_coeffs, deq_prms.scale, deq_prms.shift);

        memcpy(&_dst[sb_pos + 0], &cg_coeffs[0], sizeof(int16_t) * 16);

        num_cg--;

        for( ;num_cg > 0; --num_cg){
            sig_sb_flg = ovcabac_read_ae_significant_cg_flag(cabac_ctx, sig_sb_flg);
            if(sig_sb_flg){
                sb_pos    =  num_cg << 4;
                cg_offset = (num_cg << 3) * VVC_TR_CTX_STRIDE ;

                memset(cg_coeffs, 0, sizeof(int16_t) * 16);

                c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET + cg_offset];
                c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET + cg_offset];
                c_coding_ctx.sum_sig_nbs  = &num_significant[VVC_TR_CTX_OFFSET + cg_offset];

                if(num_cg < 2){
                    num_sig_c = ovcabac_read_ae_sb_2x8_sdh(cabac_ctx, cg_coeffs,
                                                           &c_coding_ctx);
                } else {
                    num_sig_c = ovcabac_read_ae_sb_2x8_far_sdh(cabac_ctx, cg_coeffs,
                                                               &c_coding_ctx);
                }

                deq_prms.dequant_sb(cg_coeffs, deq_prms.scale, deq_prms.shift);

                memcpy(&_dst[sb_pos], &cg_coeffs[0], sizeof(int16_t) * 16);
            }
        }

        memset(cg_coeffs, 0, sizeof(uint16_t) * 16);

        c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET];
        c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET];
        c_coding_ctx.sum_sig_nbs  = &num_significant[VVC_TR_CTX_OFFSET];

        num_sig_c += ovcabac_read_ae_sb_2x8_last_dc_sdh(cabac_ctx, cg_coeffs,
                                                        &c_coding_ctx);

        deq_prms.dequant_sb(cg_coeffs, deq_prms.scale, deq_prms.shift);

        memcpy(&_dst[0], &cg_coeffs[0],  sizeof(int16_t) * 16);

        return num_sig_c;

    } else {
        last_y    = (last_pos >> 8) & 0x1F;
        last_cg_y = last_y >> 4;

        if(!last_cg_y){
            num_sig_c = ovcabac_read_ae_sb_1x16_dc_sdh(cabac_ctx, cg_coeffs,
                                                       last_y, &c_coding_ctx);

            deq_prms.dequant_sb(cg_coeffs, deq_prms.scale, deq_prms.shift);

            memcpy(&_dst[0],  &cg_coeffs[0],  sizeof(int16_t) * 16);

            return num_sig_c;
        }

        y = last_y - (last_cg_y << 4);

        sb_pos    = last_cg_y << 4;
        cg_offset = last_cg_y << 4;

        c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET + cg_offset];
        c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET + cg_offset];
        c_coding_ctx.sum_sig_nbs  = &num_significant[VVC_TR_CTX_OFFSET + cg_offset];

        num_sig_c = ovcabac_read_ae_sb_1x16_first_sdh(cabac_ctx, cg_coeffs,
                                                      y, &c_coding_ctx);

        deq_prms.dequant_sb(cg_coeffs, deq_prms.scale, deq_prms.shift);

        memcpy(&_dst[sb_pos + 0], &cg_coeffs[0],  sizeof(int16_t) * 16);

        memset(cg_coeffs, 0, sizeof(uint16_t) * 16);

        c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET];
        c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET];
        c_coding_ctx.sum_sig_nbs  = &num_significant[VVC_TR_CTX_OFFSET];

        num_sig_c += ovcabac_read_ae_sb_1x16_last_dc_sdh(cabac_ctx, cg_coeffs,
                                                         &c_coding_ctx);

        deq_prms.dequant_sb(cg_coeffs, deq_prms.scale, deq_prms.shift);

        memcpy(&_dst[0], &cg_coeffs[0], sizeof(int16_t) * 16);

        return num_sig_c;
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
    const uint8_t *const cg_idx_2_cg_num = ff_vvc_idx_2_num[lim_log2_w  - 2]
                                                            [lim_log2_h - 2];

    const uint8_t *const scan_cg_x = ff_vvc_scan_x[lim_log2_w - 2][lim_log2_h - 2];
    const uint8_t *const scan_cg_y = ff_vvc_scan_y[lim_log2_w - 2][lim_log2_h - 2];
    int x, y;

    int last_x, last_y;
    int last_cg_x, last_cg_y;
    int num_cg, num_sig_c;
    int d_cg;

    int start_coeff_idx;
    int cg_offset;
    /*FIXME we could store coeff in smaller buffers and adapt transform
    to limited sizes*/
    uint16_t max_num_bins = (((1 << (lim_log2_h + lim_log2_w)) << 5)
                          - ((1 << (lim_log2_h + lim_log2_w)) << 2)) >> 4;

    int16_t cg_coeffs[16] = {0};

    //TODO avoid offsets tabs
    uint8_t num_significant[VVC_TR_CTX_SIZE];
    uint8_t sum_abs_level  [VVC_TR_CTX_SIZE];
    uint8_t sum_abs_level2 [VVC_TR_CTX_SIZE];

    VVCCoeffCodingCtx c_coding_ctx = {
        .sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET],
        .sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET],
        .sum_sig_nbs  = &num_significant[VVC_TR_CTX_OFFSET],
        .num_remaining_bins = max_num_bins,
        .enable_sdh = ctu_dec->enable_sdh
    };

    uint64_t sig_sb_map = 0;
    int sb_pos;

    int qp = ctu_dec->dequant_luma.qp;

    struct IQScale deq_prms = derive_dequant_sdh(qp, log2_tb_w, log2_tb_h);

    memset(_dst, 0, sizeof(int16_t) * (1 << (log2_tb_w + log2_tb_h)));

    if (!last_pos){
        ovcabac_read_ae_sb_dc_coeff_sdh(cabac_ctx, cg_coeffs);

        deq_prms.dequant_sb(cg_coeffs, deq_prms.scale, deq_prms.shift);

        _dst[0] = cg_coeffs[0];
        memcpy(ctu_dec->lfnst_subblock, cg_coeffs, sizeof(int16_t) * 16);
        return 1;
    }

    reset_ctx_buffers(&c_coding_ctx, lim_log2_w, lim_log2_h);

    last_x =  last_pos       & 0x1F;
    last_y = (last_pos >> 8) & 0x1F;
    last_cg_x = last_x >> 2;
    last_cg_y = last_y >> 2;

    sig_sb_map |= 1llu << (last_cg_x + (last_cg_y << 3));

    if (!last_cg_x && !last_cg_y){
        int last_coeff_idx = last_x + (last_y << 2);
        int num_coeffs = ff_vvc_diag_scan_4x4_num_cg [last_coeff_idx];

        num_sig_c = ovcabac_read_ae_sb_4x4_first_sdh(cabac_ctx, cg_coeffs,
                                                     num_coeffs, 0, &c_coding_ctx);

        deq_prms.dequant_sb(cg_coeffs, deq_prms.scale, deq_prms.shift);

        memcpy(&_dst[0]             , &cg_coeffs[ 0], sizeof(int16_t) * 4);
        memcpy(&_dst[1 << log2_tb_w], &cg_coeffs[ 4], sizeof(int16_t) * 4);
        memcpy(&_dst[2 << log2_tb_w], &cg_coeffs[ 8], sizeof(int16_t) * 4);
        memcpy(&_dst[3 << log2_tb_w], &cg_coeffs[12], sizeof(int16_t) * 4);
        memcpy(ctu_dec->lfnst_subblock, cg_coeffs, sizeof(int16_t) * 16);

        return sig_sb_map;
    }

    sb_pos    = (last_cg_x << 2) + ((last_cg_y << log2_tb_w) << 2);
    cg_offset = (last_cg_x << 2) + (last_cg_y << 2) * VVC_TR_CTX_STRIDE;

    d_cg = last_cg_x + last_cg_y;

    x = last_x & 3;
    y = last_y & 3;

    start_coeff_idx = ff_vvc_diag_scan_4x4_num_cg[x + (y << 2)];

    c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET + cg_offset];
    c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET + cg_offset];
    c_coding_ctx.sum_sig_nbs  = &num_significant[VVC_TR_CTX_OFFSET + cg_offset];

    num_sig_c = ovcabac_read_ae_sb_4x4_first_sdh(cabac_ctx, cg_coeffs,
                                                 start_coeff_idx, d_cg,
                                                 &c_coding_ctx);

    deq_prms.dequant_sb(cg_coeffs, deq_prms.scale, deq_prms.shift);

    memcpy(&_dst[sb_pos + (0)]             , &cg_coeffs[ 0], sizeof(int16_t) * 4);
    memcpy(&_dst[sb_pos + (1 << log2_tb_w)], &cg_coeffs[ 4], sizeof(int16_t) * 4);
    memcpy(&_dst[sb_pos + (2 << log2_tb_w)], &cg_coeffs[ 8], sizeof(int16_t) * 4);
    memcpy(&_dst[sb_pos + (3 << log2_tb_w)], &cg_coeffs[12], sizeof(int16_t) * 4);

    num_cg = cg_idx_2_cg_num[last_cg_x + last_cg_y * ((1 << lim_log2_w) >> 2)];

    num_cg--;

    for(; num_cg > 0; --num_cg){
        int x_cg = scan_cg_x[num_cg];
        int y_cg = scan_cg_y[num_cg];

        uint8_t sig_sb_flg;
        uint8_t sig_sb_blw = (((sig_sb_map >> ((y_cg + 1) << 3)) & 0xFF) >> x_cg) & 0x1;
        uint8_t sig_sb_rgt = (((sig_sb_map >> ( y_cg      << 3)) & 0xFF) >> (x_cg + 1)) & 0x1;

        sig_sb_flg = ovcabac_read_ae_significant_cg_flag(cabac_ctx, !!(sig_sb_rgt | sig_sb_blw));

        if(sig_sb_flg){

            memset(cg_coeffs, 0, sizeof(int16_t) * 16);

            sb_pos    = (x_cg << 2) + ((y_cg << log2_tb_w) << 2);
            cg_offset = (x_cg << 2) + (y_cg << 2) * (VVC_TR_CTX_STRIDE);

            sig_sb_map |= 1llu << (x_cg + (y_cg << 3));

            d_cg = x_cg + y_cg;

            c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET + cg_offset];
            c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET + cg_offset];
            c_coding_ctx.sum_sig_nbs  = &num_significant[VVC_TR_CTX_OFFSET + cg_offset];

            num_sig_c += ovcabac_read_ae_sb_4x4_sdh(cabac_ctx, cg_coeffs,
                                                    d_cg, &c_coding_ctx);

            deq_prms.dequant_sb(cg_coeffs, deq_prms.scale, deq_prms.shift);

            memcpy(&_dst[sb_pos + (0)]             , &cg_coeffs[ 0], sizeof(int16_t) * 4);
            memcpy(&_dst[sb_pos + (1 << log2_tb_w)], &cg_coeffs[ 4], sizeof(int16_t) * 4);
            memcpy(&_dst[sb_pos + (2 << log2_tb_w)], &cg_coeffs[ 8], sizeof(int16_t) * 4);
            memcpy(&_dst[sb_pos + (3 << log2_tb_w)], &cg_coeffs[12], sizeof(int16_t) * 4);
        }
    }

    memset(cg_coeffs, 0, sizeof(uint16_t) * 16);

    c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET];
    c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET];
    c_coding_ctx.sum_sig_nbs  = &num_significant[VVC_TR_CTX_OFFSET];

    num_sig_c += ovcabac_read_ae_sb_4x4_last_dc_sdh(cabac_ctx, cg_coeffs,
                                                    &c_coding_ctx);

    deq_prms.dequant_sb(cg_coeffs, deq_prms.scale, deq_prms.shift);

    memcpy(&_dst[0]             , &cg_coeffs[ 0], sizeof(int16_t) * 4);
    memcpy(&_dst[1 << log2_tb_w], &cg_coeffs[ 4], sizeof(int16_t) * 4);
    memcpy(&_dst[2 << log2_tb_w], &cg_coeffs[ 8], sizeof(int16_t) * 4);
    memcpy(&_dst[3 << log2_tb_w], &cg_coeffs[12], sizeof(int16_t) * 4);

    return sig_sb_map;
}

int
residual_coding_chroma_sdh(OVCTUDec *const ctu_dec, int16_t *const dst,
                           unsigned int log2_tb_w, unsigned int log2_tb_h,
                           uint16_t last_pos)
{
    OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;

    uint16_t max_num_bins = ((1 << (log2_tb_w + log2_tb_h) << 5)
                          - ((1 << (log2_tb_w + log2_tb_h)) << 2)) >> 4;

    uint8_t num_significant[VVC_TR_CTX_SIZE];
    uint8_t sum_abs_level  [VVC_TR_CTX_SIZE];
    uint8_t sum_abs_level2 [VVC_TR_CTX_SIZE];

    VVCCoeffCodingCtx c_coding_ctx = {
        .sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET],
        .sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET],
        .sum_sig_nbs  = &num_significant[VVC_TR_CTX_OFFSET],
        .num_remaining_bins = max_num_bins,
        .enable_sdh = ctu_dec->enable_sdh
    };

    int16_t cg_coeffs[16] = {0}; //temporary table to store coeffs in process

    uint64_t sig_sb_map = 0;

    int qp = ctu_dec->dequant_chroma->qp;

    struct IQScale deq_prms = derive_dequant_sdh(qp, log2_tb_w, log2_tb_h);

    memset(dst, 0, sizeof(int16_t) * (1 << (log2_tb_w + log2_tb_h)));

    if (!last_pos) {
        ovcabac_read_ae_sb_dc_coeff_c_sdh(cabac_ctx, cg_coeffs);
        deq_prms.dequant_sb(cg_coeffs, deq_prms.scale, deq_prms.shift);
        dst[0] = cg_coeffs [0];
        memcpy(ctu_dec->lfnst_subblock, cg_coeffs, sizeof(int16_t) * 16);
        return last_pos;
    }

    reset_ctx_buffers(&c_coding_ctx, log2_tb_w, log2_tb_h);

    if (log2_tb_w > 1 && log2_tb_h > 1) {
        const uint8_t *const cg_idx_2_cg_num = ff_vvc_idx_2_num[log2_tb_w - 2]
                                                               [log2_tb_h - 2];

        const uint8_t *const scan_cg_x = ff_vvc_scan_x[log2_tb_w - 2][log2_tb_h - 2];
        const uint8_t *const scan_cg_y = ff_vvc_scan_y[log2_tb_w - 2][log2_tb_h - 2];

        int num_sig_coeff;
        int cg_offset;
        int sb_pos;

        int x =  last_pos       & 0x1F;
        int y = (last_pos >> 8) & 0x1F;

        int last_cg_x = x >> 2;
        int last_cg_y = y >> 2;

        int num_coeffs = ff_vvc_diag_scan_4x4_num_cg[(x & 0x3) + ((y & 0x3) << 2)];

        int num_cg = cg_idx_2_cg_num[last_cg_x + ((last_cg_y << log2_tb_w) >> 2)];

        if (!num_cg) {
            ovcabac_read_ae_sb_4x4_dc_c_sdh(cabac_ctx, cg_coeffs,
                                            num_coeffs, &c_coding_ctx);

            deq_prms.dequant_sb(cg_coeffs, deq_prms.scale, deq_prms.shift);

            memcpy(&dst[0]             , &cg_coeffs[ 0], sizeof(int16_t) * 4);
            memcpy(&dst[1 << log2_tb_w], &cg_coeffs[ 4], sizeof(int16_t) * 4);
            memcpy(&dst[2 << log2_tb_w], &cg_coeffs[ 8], sizeof(int16_t) * 4);
            memcpy(&dst[3 << log2_tb_w], &cg_coeffs[12], sizeof(int16_t) * 4);

            memcpy(ctu_dec->lfnst_subblock, cg_coeffs, sizeof(int16_t) * 16);

            return num_coeffs;
        }

        sb_pos = (last_cg_x << 2) + (((last_cg_y << log2_tb_w) >> 2) << 4);

        cg_offset = (last_cg_x << 2) + (last_cg_y << 2) * (VVC_TR_CTX_STRIDE);

        c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET + cg_offset];
        c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET + cg_offset];
        c_coding_ctx.sum_sig_nbs  = &num_significant[VVC_TR_CTX_OFFSET + cg_offset];

        num_sig_coeff = ovcabac_read_ae_sb_4x4_first_c_sdh(cabac_ctx, cg_coeffs,
                                                           num_coeffs,
                                                           &c_coding_ctx);

        deq_prms.dequant_sb(cg_coeffs, deq_prms.scale, deq_prms.shift);

        memcpy(&dst[sb_pos + (0)]             , &cg_coeffs[ 0], sizeof(int16_t) * 4);
        memcpy(&dst[sb_pos + (1 << log2_tb_w)], &cg_coeffs[ 4], sizeof(int16_t) * 4);
        memcpy(&dst[sb_pos + (2 << log2_tb_w)], &cg_coeffs[ 8], sizeof(int16_t) * 4);
        memcpy(&dst[sb_pos + (3 << log2_tb_w)], &cg_coeffs[12], sizeof(int16_t) * 4);

        sig_sb_map |= 1llu << (last_cg_x + (last_cg_y << 3));

        num_cg--;

        for(int i = num_cg; i > 0; --i){
            int x_cg = scan_cg_x[i];
            int y_cg = scan_cg_y[i];
            uint8_t sig_sb_blw = (((sig_sb_map >> ((y_cg + 1) << 3)) & 0xFF) >> x_cg) & 0x1;
            uint8_t sig_sb_rgt = (((sig_sb_map >> ( y_cg      << 3)) & 0xFF) >> (x_cg + 1)) & 0x1;

            uint8_t sig_sb_flg = ovcabac_read_ae_significant_cg_flag_chroma(cabac_ctx, !!(sig_sb_rgt | sig_sb_blw));

            if(sig_sb_flg){

                memset(cg_coeffs, 0, sizeof(int16_t) * 16);

                sb_pos = (x_cg << 2) + (((y_cg << log2_tb_w) >> 2) << 4);
                cg_offset = (x_cg << 2) + (y_cg << 2) * (VVC_TR_CTX_STRIDE);

                sig_sb_map |= 1llu << (x_cg + (y_cg << 3));

                c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET + cg_offset];
                c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET + cg_offset];
                c_coding_ctx.sum_sig_nbs  = &num_significant[VVC_TR_CTX_OFFSET + cg_offset];

                num_sig_coeff += ovcabac_read_ae_sb_4x4_c_sdh(cabac_ctx, cg_coeffs,
                                                              &c_coding_ctx);

                deq_prms.dequant_sb(cg_coeffs, deq_prms.scale, deq_prms.shift);

                memcpy(&dst[sb_pos + (0)]             , &cg_coeffs[ 0], sizeof(int16_t) * 4);
                memcpy(&dst[sb_pos + (1 << log2_tb_w)], &cg_coeffs[ 4], sizeof(int16_t) * 4);
                memcpy(&dst[sb_pos + (2 << log2_tb_w)], &cg_coeffs[ 8], sizeof(int16_t) * 4);
                memcpy(&dst[sb_pos + (3 << log2_tb_w)], &cg_coeffs[12], sizeof(int16_t) * 4);
            }
        }

        memset(cg_coeffs, 0, sizeof(int16_t) * 16);

        c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET];
        c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET];
        c_coding_ctx.sum_sig_nbs  = &num_significant[VVC_TR_CTX_OFFSET];

        num_sig_coeff += ovcabac_read_ae_sb_4x4_last_dc_c_sdh(cabac_ctx, cg_coeffs,
                                                              &c_coding_ctx);

        deq_prms.dequant_sb(cg_coeffs, deq_prms.scale, deq_prms.shift);

        memcpy(&dst[0]             , &cg_coeffs[ 0], sizeof(int16_t) * 4);
        memcpy(&dst[1 << log2_tb_w], &cg_coeffs[ 4], sizeof(int16_t) * 4);
        memcpy(&dst[2 << log2_tb_w], &cg_coeffs[ 8], sizeof(int16_t) * 4);
        memcpy(&dst[3 << log2_tb_w], &cg_coeffs[12], sizeof(int16_t) * 4);

        return 0xFFFF;

    } else if (log2_tb_h == 1) {
        int start_coeff_idx;
        int num_sig_c;
        int cg_offset;
        int sb_pos;
        int num_cg;
        uint8_t sig_sb_flg = 1;

        int x =  last_pos       & 0x1F;
        int y = (last_pos >> 8) & 0x1F;

        int last_cg_x = x >> 3;

        if (!last_cg_x) {
            int last_coeff_idx = x + (y << 3) ;
            int num_coeffs = ff_vvc_diag_scan_8x2_num_cg [last_coeff_idx];
            num_sig_c = ovcabac_read_ae_sb_8x2_dc_c_sdh(cabac_ctx, cg_coeffs,
                                                        num_coeffs, &c_coding_ctx);

            deq_prms.dequant_sb(cg_coeffs, deq_prms.scale, deq_prms.shift);

            memcpy(&dst[0        ]     , &cg_coeffs[0], sizeof(int16_t) * 8);
            memcpy(&dst[1 << log2_tb_w], &cg_coeffs[8], sizeof(int16_t) * 8);
            /*FIXME determine whether or not we should read lfnst*/

            return 0xFFFF;
        }

        x -= last_cg_x << 3;

        start_coeff_idx =  ff_vvc_diag_scan_8x2_num_cg [x + (y << 3)];

        sb_pos    = last_cg_x << 3;
        cg_offset = last_cg_x << 3;

        c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET + cg_offset];
        c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET + cg_offset];
        c_coding_ctx.sum_sig_nbs  = &num_significant[VVC_TR_CTX_OFFSET + cg_offset];

        num_sig_c = ovcabac_read_ae_sb_8x2_first_c_sdh(cabac_ctx, cg_coeffs,
                                                       start_coeff_idx, &c_coding_ctx);

        deq_prms.dequant_sb(cg_coeffs, deq_prms.scale, deq_prms.shift);

        memcpy(&dst[sb_pos + 0               ], &cg_coeffs[0], sizeof(int16_t) * 8);
        memcpy(&dst[sb_pos + (1 << log2_tb_w)], &cg_coeffs[8], sizeof(int16_t) * 8);

        num_cg = last_cg_x;
        num_cg--;

        for(;num_cg > 0; --num_cg){

            sig_sb_flg = ovcabac_read_ae_significant_cg_flag_chroma(cabac_ctx, sig_sb_flg);

            if(sig_sb_flg){

                memset(cg_coeffs, 0, sizeof(int16_t) * 16);

                sb_pos    = num_cg << 3;
                cg_offset = num_cg << 3;

                c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET + cg_offset];
                c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET + cg_offset];
                c_coding_ctx.sum_sig_nbs  = &num_significant[VVC_TR_CTX_OFFSET + cg_offset];

                num_sig_c = ovcabac_read_ae_sb_8x2_c_sdh(cabac_ctx, cg_coeffs,
                                                         &c_coding_ctx);

                deq_prms.dequant_sb(cg_coeffs, deq_prms.scale, deq_prms.shift);

                memcpy(&dst[sb_pos + 0               ], &cg_coeffs[0], sizeof(int16_t) * 8);
                memcpy(&dst[sb_pos + (1 << log2_tb_w)], &cg_coeffs[8], sizeof(int16_t) * 8);
            }
        }

        memset(cg_coeffs, 0, sizeof(uint16_t) * 16);

        c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET];
        c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET];
        c_coding_ctx.sum_sig_nbs  = &num_significant[VVC_TR_CTX_OFFSET];

        num_sig_c += ovcabac_read_ae_sb_8x2_last_dc_c_sdh(cabac_ctx, cg_coeffs,
                                                          &c_coding_ctx);

        deq_prms.dequant_sb(cg_coeffs, deq_prms.scale, deq_prms.shift);

        memcpy(&dst[0             ], &cg_coeffs[0], sizeof(int16_t) * 8);
        memcpy(&dst[1 << log2_tb_w], &cg_coeffs[8], sizeof(int16_t) * 8);

        return 0xFFFF;

    } else if(log2_tb_w == 1) {
        uint8_t sig_sb_flg = 1;
        int num_sig_c;
        int sb_pos;
        int cg_offset;
        int start_coeff_idx;
        int num_cg;

        int x =  last_pos       & 0x1F;
        int y = (last_pos >> 8) & 0x1F;
        int last_cg_y = y >> 3;

        if (!last_cg_y) {
            int last_coeff_idx = x + (y << 1);
            int num_coeffs = ff_vvc_diag_scan_2x8_num_cg [last_coeff_idx];

            num_sig_c = ovcabac_read_ae_sb_2x8_dc_c_sdh(cabac_ctx, cg_coeffs,
                                                        num_coeffs, &c_coding_ctx);

            deq_prms.dequant_sb(cg_coeffs, deq_prms.scale, deq_prms.shift);

            memcpy(&dst[0], &cg_coeffs[0],  sizeof(int16_t) * 16);
            /*FIXME determine whether or not we should read lfnst*/

            return 0xFFFF;
        }

        sb_pos    =  last_cg_y << 4;
        cg_offset = (last_cg_y << 3) * VVC_TR_CTX_STRIDE;

        num_cg = last_cg_y;

        x = x;
        y = y - (last_cg_y << 3);

        start_coeff_idx =  ff_vvc_diag_scan_2x8_num_cg [x + (y << 1)];

        c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET + cg_offset];
        c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET + cg_offset];
        c_coding_ctx.sum_sig_nbs  = &num_significant[VVC_TR_CTX_OFFSET + cg_offset];

        num_sig_c = ovcabac_read_ae_sb_2x8_first_c_sdh(cabac_ctx, cg_coeffs,
                                                       start_coeff_idx, &c_coding_ctx);

        deq_prms.dequant_sb(cg_coeffs, deq_prms.scale, deq_prms.shift);

        memcpy(&dst[sb_pos + 0], &cg_coeffs[0], sizeof(int16_t) * 16);

        num_cg--;

        for( ;num_cg > 0; --num_cg){

            sig_sb_flg = ovcabac_read_ae_significant_cg_flag_chroma(cabac_ctx, sig_sb_flg);

            if(sig_sb_flg){
                sb_pos    = (num_cg << 4);
                cg_offset = (num_cg << 3) * VVC_TR_CTX_STRIDE;

                memset(cg_coeffs, 0, sizeof(int16_t) * 16);

                c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET + cg_offset];
                c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET + cg_offset];
                c_coding_ctx.sum_sig_nbs  = &num_significant[VVC_TR_CTX_OFFSET + cg_offset];

                num_sig_c = ovcabac_read_ae_sb_2x8_c_sdh(cabac_ctx, cg_coeffs,
                                                         &c_coding_ctx);

                deq_prms.dequant_sb(cg_coeffs, deq_prms.scale, deq_prms.shift);

                memcpy(&dst[sb_pos], &cg_coeffs[0], sizeof(int16_t) * 16);
            }
        }

        memset(cg_coeffs, 0, sizeof(uint16_t) * 16);

        c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET];
        c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET];
        c_coding_ctx.sum_sig_nbs  = &num_significant[VVC_TR_CTX_OFFSET];

        num_sig_c += ovcabac_read_ae_sb_2x8_last_dc_c_sdh(cabac_ctx, cg_coeffs,
                                                          &c_coding_ctx);

        deq_prms.dequant_sb(cg_coeffs, deq_prms.scale, deq_prms.shift);

        memcpy(&dst[0], &cg_coeffs[0],  sizeof(int16_t) * 16);

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
    int last_cg_x;
    int num_sig_c;
    int start_coeff_idx;
    int x, y;
    int cg_offset;
    int sb_pos;
    int num_cg;
    uint8_t sig_sb_flg = 1;
    int16_t cg_coeffs[16] = {0};
    uint8_t num_significant [VVC_TR_CTX_SIZE];
    uint8_t sum_abs_level   [VVC_TR_CTX_SIZE];
    uint8_t sum_abs_level2 [VVC_TR_CTX_SIZE];

    uint16_t max_num_bins = ((1 << (log2_tb_w + log2_tb_h))
                             * 28) >> 4;

    int qp = ctu_dec->dequant_luma.qp;

    struct IQScale deq_prms = derive_dequant_dpq(qp, log2_tb_w, log2_tb_h);

    VVCCoeffCodingCtx c_coding_ctx = {
        .sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET],
        .sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET],
        .sum_sig_nbs  = &num_significant[VVC_TR_CTX_OFFSET],
        .num_remaining_bins = max_num_bins,
        .enable_sdh = ctu_dec->enable_sdh
    };

    memset(_dst, 0, sizeof(int16_t) * (1 << (log2_tb_w + log2_tb_h)));

    if(!last_pos){
        ovcabac_read_ae_sb_dc_coeff_dpq(cabac_ctx, cg_coeffs);

        deq_prms.dequant_sb(cg_coeffs, deq_prms.scale, deq_prms.shift);

        _dst[0] = cg_coeffs[0];

        return 1;
    }

    reset_ctx_buffers(&c_coding_ctx, log2_tb_w, log2_tb_h);

    if(log2_tb_h){
        last_x =  last_pos       & 0x1F;
        last_y = (last_pos >> 8) & 0x1F;
        last_cg_x = last_x >> 3  ;

        if(!last_cg_x){
            int last_coeff_idx = last_x + (last_y << 3) ;
            int num_coeffs = ff_vvc_diag_scan_8x2_num_cg [last_coeff_idx];
            num_sig_c = ovcabac_read_ae_sb_8x2_dc_dpq(cabac_ctx, cg_coeffs,
                                                      &state, num_coeffs,
                                                      &c_coding_ctx);

            deq_prms.dequant_sb(cg_coeffs, deq_prms.scale, deq_prms.shift);

            memcpy(&_dst[0             ], &cg_coeffs[0], sizeof(int16_t) * 8);
            memcpy(&_dst[1 << log2_tb_w], &cg_coeffs[8], sizeof(int16_t) * 8);

            return num_sig_c;
        }

        x = last_x - (last_cg_x << 3) ;
        y = last_y;

        start_coeff_idx = ff_vvc_diag_scan_8x2_num_cg [x + y * (1 << 3)];

        sb_pos    = last_cg_x << 3;
        cg_offset = last_cg_x << 3;

        c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET + cg_offset];
        c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET + cg_offset];
        c_coding_ctx.sum_sig_nbs  = &num_significant[VVC_TR_CTX_OFFSET + cg_offset];

        if(last_cg_x < 2){
            num_sig_c = ovcabac_read_ae_sb_8x2_first_dpq(cabac_ctx, cg_coeffs,
                                                    &state, start_coeff_idx,
                                                    &c_coding_ctx);
        } else {
            num_sig_c = ovcabac_read_ae_sb_8x2_first_far_dpq(cabac_ctx, cg_coeffs,
                                                        &state, start_coeff_idx,
                                                        &c_coding_ctx);
        }

        deq_prms.dequant_sb(cg_coeffs, deq_prms.scale, deq_prms.shift);

        memcpy(&_dst[sb_pos + 0               ], &cg_coeffs[0], sizeof(int16_t) * 8);
        memcpy(&_dst[sb_pos + (1 << log2_tb_w)], &cg_coeffs[8], sizeof(int16_t) * 8);

        num_cg = last_cg_x;
        num_cg--;

        for(;num_cg > 0; --num_cg){
            sig_sb_flg = ovcabac_read_ae_significant_cg_flag(cabac_ctx, sig_sb_flg);
            if(sig_sb_flg){

                memset(cg_coeffs, 0, sizeof(int16_t) * 16);

                sb_pos    = num_cg << 3;
                cg_offset = num_cg << 3;

                c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET + cg_offset];
                c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET + cg_offset];
                c_coding_ctx.sum_sig_nbs  = &num_significant[VVC_TR_CTX_OFFSET + cg_offset];

                if(num_cg < 2){
                    num_sig_c = ovcabac_read_ae_sb_8x2_dpq(cabac_ctx, cg_coeffs,
                                                           &state, &c_coding_ctx);
                } else {
                    num_sig_c = ovcabac_read_ae_sb_8x2_far_dpq(cabac_ctx, cg_coeffs,
                                                               &state, &c_coding_ctx);
                }

                deq_prms.dequant_sb(cg_coeffs, deq_prms.scale, deq_prms.shift);

                memcpy(&_dst[sb_pos + 0               ], &cg_coeffs[0],  sizeof(int16_t) * 8);
                memcpy(&_dst[sb_pos + (1 << log2_tb_w)], &cg_coeffs[8],  sizeof(int16_t) * 8);
            }
        }

        memset(cg_coeffs, 0, sizeof(uint16_t) * 16);

        c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET],
        c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET],
        c_coding_ctx.sum_sig_nbs  = &num_significant[VVC_TR_CTX_OFFSET],

        num_sig_c += ovcabac_read_ae_sb_8x2_last_dc_dpq(cabac_ctx,
                                                        cg_coeffs, &state,
                                                        &c_coding_ctx);

        deq_prms.dequant_sb(cg_coeffs, deq_prms.scale, deq_prms.shift);

        memcpy(&_dst[0             ], &cg_coeffs[0],  sizeof(int16_t) * 8);
        memcpy(&_dst[1 << log2_tb_w], &cg_coeffs[8],  sizeof(int16_t) * 8);

        return num_sig_c;
    } else {
        last_x = last_pos & 0x1F;

        last_cg_x = last_x >> 4 ;

        if(!last_cg_x){
            int last_coeff_idx = last_x;
            int num_coeffs = last_coeff_idx;

            num_sig_c = ovcabac_read_ae_sb_1x16_dc_dpq(cabac_ctx, cg_coeffs,
                                                       &state, num_coeffs,
                                                       &c_coding_ctx);

            deq_prms.dequant_sb(cg_coeffs, deq_prms.scale, deq_prms.shift);

            memcpy(&_dst[0] , &cg_coeffs[0], sizeof(int16_t) * 16);

            return num_sig_c;
        }

        x = last_x - (last_cg_x << 4);

        cg_offset = last_cg_x << 4;
        sb_pos    = last_cg_x << 4;

        c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET + cg_offset];
        c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET + cg_offset];
        c_coding_ctx.sum_sig_nbs  = &num_significant[VVC_TR_CTX_OFFSET + cg_offset];

        num_sig_c = ovcabac_read_ae_sb_1x16_first_dpq(cabac_ctx, cg_coeffs,
                                                      &state, x, &c_coding_ctx);

        deq_prms.dequant_sb(cg_coeffs, deq_prms.scale, deq_prms.shift);

        memcpy(&_dst[sb_pos] , &cg_coeffs[0], sizeof(int16_t) * 16);

        memset(cg_coeffs, 0, sizeof(uint16_t) * 16);

        c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET],
        c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET],
        c_coding_ctx.sum_sig_nbs  = &num_significant[VVC_TR_CTX_OFFSET],

        num_sig_c += ovcabac_read_ae_sb_1x16_last_dc_dpq(cabac_ctx,
                                                         cg_coeffs, &state,
                                                         &c_coding_ctx);

        deq_prms.dequant_sb(cg_coeffs, deq_prms.scale, deq_prms.shift);

        memcpy(&_dst[0] , &cg_coeffs[0], sizeof(int16_t) * 16);

        return num_sig_c;
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
    int last_cg_y;
    int num_sig_c;
    int start_coeff_idx;
    int x, y;
    int cg_offset;
    int sb_pos;
    int num_cg;
    uint8_t sig_sb_flg = 1;
    int16_t cg_coeffs[16] = {0};
    uint8_t  num_significant[VVC_TR_CTX_SIZE];
    uint8_t  sum_abs_level  [VVC_TR_CTX_SIZE];
    uint8_t sum_abs_level2 [VVC_TR_CTX_SIZE];
    uint16_t max_num_bins = ((1 << (log2_tb_w + log2_tb_h)) * 28) >> 4;

    VVCCoeffCodingCtx c_coding_ctx = {
        .sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET],
        .sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET],
        .sum_sig_nbs  = &num_significant[VVC_TR_CTX_OFFSET],
        .num_remaining_bins = max_num_bins,
        .enable_sdh = ctu_dec->enable_sdh
    };

    int qp = ctu_dec->dequant_luma.qp;

    struct IQScale deq_prms = derive_dequant_dpq(qp, log2_tb_w, log2_tb_h);

    memset(_dst, 0, sizeof(int16_t) * (1 << (log2_tb_w + log2_tb_h)));

    if(!last_pos){

        ovcabac_read_ae_sb_dc_coeff_dpq(cabac_ctx, cg_coeffs);

        deq_prms.dequant_sb(cg_coeffs, deq_prms.scale, deq_prms.shift);

        _dst[0] = cg_coeffs[0];
        return 1;
    }

    if (log2_tb_w) reset_ctx_buffers(&c_coding_ctx, log2_tb_w, log2_tb_h);
    else reset_ctx_buffers(&c_coding_ctx, log2_tb_h, log2_tb_w);

    if(log2_tb_w){
        last_x =  last_pos       & 0x1F;
        last_y = (last_pos >> 8) & 0x1F;
        last_cg_y = last_y >> 3;
        if(!last_cg_y){
            int last_coeff_idx = last_x + (last_y << 1);
            int num_coeffs = ff_vvc_diag_scan_2x8_num_cg [last_coeff_idx];
            num_sig_c = ovcabac_read_ae_sb_2x8_dc_dpq(cabac_ctx, cg_coeffs,
                                                      &state, num_coeffs,
                                                      &c_coding_ctx);

            deq_prms.dequant_sb(cg_coeffs, deq_prms.scale, deq_prms.shift);

            memcpy(&_dst[0], &cg_coeffs[0],  sizeof(int16_t) * 16);

            return num_sig_c;
        }
        sb_pos    =  last_cg_y << 4;
        cg_offset = (last_cg_y << 3) * VVC_TR_CTX_STRIDE;

        num_cg = last_cg_y;

        x = last_x;
        y = last_y - (last_cg_y << 3);

        start_coeff_idx = ff_vvc_diag_scan_2x8_num_cg [x + (y << 1)];

        c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET + cg_offset];
        c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET + cg_offset];
        c_coding_ctx.sum_sig_nbs  = &num_significant[VVC_TR_CTX_OFFSET + cg_offset];

        if(num_cg < 2){
            num_sig_c = ovcabac_read_ae_sb_2x8_first_dpq(cabac_ctx, cg_coeffs,
                                                         &state, start_coeff_idx,
                                                         &c_coding_ctx);
        } else {
            num_sig_c = ovcabac_read_ae_sb_2x8_first_far_dpq(cabac_ctx, cg_coeffs,
                                                             &state, start_coeff_idx,
                                                             &c_coding_ctx);
        }

        deq_prms.dequant_sb(cg_coeffs, deq_prms.scale, deq_prms.shift);

        memcpy(&_dst[sb_pos + 0], &cg_coeffs[0], sizeof(int16_t) * 16);

        num_cg--;

        for( ;num_cg > 0; --num_cg){
            sig_sb_flg = ovcabac_read_ae_significant_cg_flag(cabac_ctx, sig_sb_flg);
            if(sig_sb_flg){
                sb_pos    =  num_cg << 4;
                cg_offset = (num_cg << 3) * VVC_TR_CTX_STRIDE ;

                memset(cg_coeffs, 0, sizeof(int16_t) * 16);

                c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET + cg_offset];
                c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET + cg_offset];
                c_coding_ctx.sum_sig_nbs  = &num_significant[VVC_TR_CTX_OFFSET + cg_offset];

                if(num_cg < 2){
                    num_sig_c = ovcabac_read_ae_sb_2x8_dpq(cabac_ctx, cg_coeffs,
                                                           &state, &c_coding_ctx);
                } else {
                    num_sig_c = ovcabac_read_ae_sb_2x8_far_dpq(cabac_ctx, cg_coeffs,
                                                               &state, &c_coding_ctx);
                }

                deq_prms.dequant_sb(cg_coeffs, deq_prms.scale, deq_prms.shift);

                memcpy(&_dst[sb_pos], &cg_coeffs[0], sizeof(int16_t) * 16);
            }
        }

        memset(cg_coeffs, 0, sizeof(uint16_t) * 16);

        c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET];
        c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET];
        c_coding_ctx.sum_sig_nbs  = &num_significant[VVC_TR_CTX_OFFSET];

        num_sig_c += ovcabac_read_ae_sb_2x8_last_dc_dpq(cabac_ctx,
                                                        cg_coeffs, &state,
                                                        &c_coding_ctx);

        deq_prms.dequant_sb(cg_coeffs, deq_prms.scale, deq_prms.shift);

        memcpy(&_dst[0], &cg_coeffs[0],  sizeof(int16_t) * 16);

        return num_sig_c;

    } else {
        last_y    = (last_pos >> 8) & 0x1F;
        last_cg_y = last_y >> 4;

        if(!last_cg_y){
            num_sig_c = ovcabac_read_ae_sb_1x16_dc_dpq(cabac_ctx, cg_coeffs,
                                                       &state, last_y, &c_coding_ctx);

            deq_prms.dequant_sb(cg_coeffs, deq_prms.scale, deq_prms.shift);

            memcpy(&_dst[0],  &cg_coeffs[0],  sizeof(int16_t) * 16);

            return num_sig_c;
        }

        y = last_y - (last_cg_y << 4);

        sb_pos    = last_cg_y << 4;
        cg_offset = last_cg_y << 4;

        c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET + cg_offset];
        c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET + cg_offset];
        c_coding_ctx.sum_sig_nbs  = &num_significant[VVC_TR_CTX_OFFSET + cg_offset];

        num_sig_c = ovcabac_read_ae_sb_1x16_first_dpq(cabac_ctx, cg_coeffs,
                                                      &state, y, &c_coding_ctx);

        deq_prms.dequant_sb(cg_coeffs, deq_prms.scale, deq_prms.shift);

        memcpy(&_dst[sb_pos + 0], &cg_coeffs[0],  sizeof(int16_t) * 16);

        memset(cg_coeffs, 0, sizeof(uint16_t) * 16);

        c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET];
        c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET];
        c_coding_ctx.sum_sig_nbs  = &num_significant[VVC_TR_CTX_OFFSET];

        num_sig_c += ovcabac_read_ae_sb_1x16_last_dc_dpq(cabac_ctx,
                                                         cg_coeffs, &state,
                                                         &c_coding_ctx);

        deq_prms.dequant_sb(cg_coeffs, deq_prms.scale, deq_prms.shift);

        memcpy(&_dst[0], &cg_coeffs[0], sizeof(int16_t) * 16);

        return num_sig_c;
    }
}


int
residual_coding_ts(OVCTUDec *const ctu_dec, unsigned int log2_tb_w, unsigned int log2_tb_h)
{
    OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;

    const uint8_t *const scan_cg_x = ff_vvc_scan_x[log2_tb_w  - 2][log2_tb_h - 2];
    const uint8_t *const scan_cg_y = ff_vvc_scan_y[log2_tb_w  - 2][log2_tb_h - 2];

    int i;
    int nb_cg, nb_sig_c;

    int nb_sig_cg = 0;
    int cg_offset;

    uint8_t log2_tb_s = log2_tb_h + log2_tb_w;

    uint16_t max_nb_bins = (((1 << log2_tb_s) << 3) - (1 << log2_tb_s)) >> 2;

    /* FIXME better tables */
    uint8_t nb_significant[VVC_TR_CTX_SIZE]={0};
    uint8_t sign_map[VVC_TR_CTX_SIZE]={0};
    uint16_t abs_coeffs[VVC_TR_CTX_SIZE]={0};

    uint8_t sig_cg_map[17*17] = {0};
    uint8_t sig_sb_flg = 0;

    /* FIXME use bit map of significant sb instead */
    uint8_t *significant_cg_map_2 = &sig_cg_map[0];

    /* FIXME if called from chroma ? */
    int qp = ctu_dec->dequant_luma_skip.qp;

    const struct IQScale deq_prms = derive_dequant_ts(qp, log2_tb_w, log2_tb_h);

    int16_t cg_coeffs[16] = {0};

    memset(&ctu_dec->transform_buff, 0, sizeof(uint16_t) << log2_tb_s);

    if (!(log2_tb_w - 2) && !(log2_tb_h - 2)) {
        int16_t *dst = &ctu_dec->transform_buff[0];

        nb_sig_c = ovcabac_read_ae_sb_ts_4x4(cabac_ctx, cg_coeffs,
                                             &nb_significant   [0],
                                             &sign_map         [0],
                                             &abs_coeffs  [36 + 0],
                                             (int16_t*)&max_nb_bins);

        deq_prms.dequant_sb(cg_coeffs, deq_prms.scale, deq_prms.shift);

        memcpy(&dst[0             ] , &cg_coeffs[0],  sizeof(int16_t) * 4);
        memcpy(&dst[1 << log2_tb_w] , &cg_coeffs[4],  sizeof(int16_t) * 4);
        memcpy(&dst[2 << log2_tb_w] , &cg_coeffs[8],  sizeof(int16_t) * 4);
        memcpy(&dst[3 << log2_tb_w] , &cg_coeffs[12], sizeof(int16_t) * 4);

        return nb_sig_c;
    }

    nb_cg = (1 << log2_tb_s) >> 4;

    for (i = 0; i < nb_cg - 1; ++i) {
        int x_cg = scan_cg_x[i];
        int y_cg = scan_cg_y[i];

        int significant_cg_offset = significant_cg_map_2 [x_cg + y_cg * (17)];

        sig_sb_flg = ovcabac_read_ae_significant_ts_cg_flag(cabac_ctx, significant_cg_offset);

        if (sig_sb_flg) {
            int16_t *dst = &ctu_dec->transform_buff[(x_cg << 2) + ((y_cg << log2_tb_w) << 2)];

            memset(cg_coeffs, 0, sizeof(int16_t) * 16);

            cg_offset = (x_cg << 2) + (y_cg << 2) * (VVC_TR_CTX_STRIDE);

            /* FIXME use bit map of significant sb instead */
            significant_cg_map_2[x_cg + 1 +  y_cg      * (17)] += 1;
            significant_cg_map_2[x_cg     + (y_cg + 1) * (17)] += 1;

            nb_sig_c += ovcabac_read_ae_sb_ts_4x4(cabac_ctx,
                                                  cg_coeffs,
                                                  &nb_significant[cg_offset],
                                                  &sign_map      [cg_offset],
                                                  &abs_coeffs [36+cg_offset],
                                                  (int16_t *)&max_nb_bins);

            deq_prms.dequant_sb(cg_coeffs, deq_prms.scale, deq_prms.shift);

            memcpy(&dst[0             ] , &cg_coeffs[0],  sizeof(int16_t) * 4);
            memcpy(&dst[1 << log2_tb_w] , &cg_coeffs[4],  sizeof(int16_t) * 4);
            memcpy(&dst[2 << log2_tb_w] , &cg_coeffs[8],  sizeof(int16_t) * 4);
            memcpy(&dst[3 << log2_tb_w] , &cg_coeffs[12], sizeof(int16_t) * 4);

            ++nb_sig_cg;
        }
    }

    sig_sb_flg = !nb_sig_cg;

    if (nb_sig_cg) {
        int x_cg = scan_cg_x[i];
        int y_cg = scan_cg_y[i];
        int significant_cg_offset = significant_cg_map_2 [x_cg + y_cg * (17)];

        sig_sb_flg = ovcabac_read_ae_significant_ts_cg_flag(cabac_ctx, significant_cg_offset);
    }

    if (sig_sb_flg) {
        int x_cg = scan_cg_x[i];
        int y_cg = scan_cg_y[i];
        int16_t *dst = &ctu_dec->transform_buff[(x_cg << 2) + ((y_cg << log2_tb_w) << 2)];

        cg_offset = (x_cg << 2) + (y_cg << 2) * (VVC_TR_CTX_STRIDE);

        memset(cg_coeffs, 0, sizeof(uint16_t) * 16);

        nb_sig_c += ovcabac_read_ae_sb_ts_4x4(cabac_ctx,
                                              cg_coeffs,
                                              &nb_significant[cg_offset],
                                              &sign_map      [cg_offset],
                                              &abs_coeffs   [36+cg_offset],
                                              (int16_t*) &max_nb_bins);

        deq_prms.dequant_sb(cg_coeffs, deq_prms.scale, deq_prms.shift);

        memcpy(&dst[0             ] , &cg_coeffs[0],  sizeof(int16_t) * 4);
        memcpy(&dst[1 << log2_tb_w] , &cg_coeffs[4],  sizeof(int16_t) * 4);
        memcpy(&dst[2 << log2_tb_w] , &cg_coeffs[8],  sizeof(int16_t) * 4);
        memcpy(&dst[3 << log2_tb_w] , &cg_coeffs[12], sizeof(int16_t) * 4);
    }

    return nb_sig_c;
}

uint64_t
residual_coding_dpq(OVCTUDec *const ctu_dec, int16_t *const dst,
                    unsigned int log2_tb_w, unsigned int log2_tb_h,
                    uint16_t last_pos)
{
    OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;
    int16_t *const _dst = dst;
    /*FIXME we can reduce LUT based on the fact it cannot be greater than 5 */
    uint8_t lim_log2_w = OVMIN(log2_tb_w, 5);
    uint8_t lim_log2_h = OVMIN(log2_tb_h, 5);

    const uint8_t *const cg_idx_2_cg_num = ff_vvc_idx_2_num[lim_log2_w  - 2]
                                                            [lim_log2_h - 2];

    const uint8_t *const scan_cg_x = ff_vvc_scan_x[lim_log2_w - 2][lim_log2_h - 2];
    const uint8_t *const scan_cg_y = ff_vvc_scan_y[lim_log2_w - 2][lim_log2_h - 2];
    int x, y;

    int state = 0;

    int last_x, last_y;
    int last_cg_x, last_cg_y;
    int num_cg, num_sig_c;
    int d_cg;

    int start_coeff_idx;
    int cg_offset;
    uint16_t max_num_bins = (((1 << (lim_log2_h + lim_log2_w)) << 5)
                          - ((1 << (lim_log2_h + lim_log2_w)) << 2)) >> 4;

    int16_t cg_coeffs[16] = {0};

    //TODO avoid offsets tabs
    uint8_t  num_significant[VVC_TR_CTX_SIZE];
    uint8_t  sum_abs_level  [VVC_TR_CTX_SIZE];
    uint8_t sum_abs_level2 [VVC_TR_CTX_SIZE];

    VVCCoeffCodingCtx c_coding_ctx = {
        .sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET],
        .sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET],
        .sum_sig_nbs  = &num_significant[VVC_TR_CTX_OFFSET],
        .num_remaining_bins = max_num_bins,
        .enable_sdh = ctu_dec->enable_sdh
    };

    uint64_t sig_sb_map = 0;
    int sb_pos;

    int qp = ctu_dec->dequant_luma.qp;

    struct IQScale deq_prms = derive_dequant_dpq(qp, log2_tb_w, log2_tb_h);

    memset(_dst, 0, sizeof(int16_t) * (1 << (log2_tb_w + log2_tb_h)));

    if (!last_pos){

        ovcabac_read_ae_sb_dc_coeff_dpq(cabac_ctx, cg_coeffs);

        deq_prms.dequant_sb(cg_coeffs, deq_prms.scale, deq_prms.shift);

        _dst[0] = cg_coeffs[0];
        memcpy(ctu_dec->lfnst_subblock, cg_coeffs, sizeof(int16_t) * 16);
        return 1;
    }

    reset_ctx_buffers(&c_coding_ctx, lim_log2_w, lim_log2_h);

    last_x =  last_pos       & 0x1F;
    last_y = (last_pos >> 8) & 0x1F;
    last_cg_x = last_x >> 2;
    last_cg_y = last_y >> 2;

    sig_sb_map |= 1llu << (last_cg_x + (last_cg_y << 3));

    if (!last_cg_x && !last_cg_y){
        int last_coeff_idx = last_x + (last_y << 2);
        int num_coeffs = ff_vvc_diag_scan_4x4_num_cg [last_coeff_idx];

        num_sig_c = ovcabac_read_ae_sb_4x4_first_dpq(cabac_ctx, cg_coeffs,
                                                     &state, num_coeffs, 0,
                                                     &c_coding_ctx);

        deq_prms.dequant_sb(cg_coeffs, deq_prms.scale, deq_prms.shift);

        memcpy(&_dst[0]             , &cg_coeffs[ 0], sizeof(int16_t) * 4);
        memcpy(&_dst[1 << log2_tb_w], &cg_coeffs[ 4], sizeof(int16_t) * 4);
        memcpy(&_dst[2 << log2_tb_w], &cg_coeffs[ 8], sizeof(int16_t) * 4);
        memcpy(&_dst[3 << log2_tb_w], &cg_coeffs[12], sizeof(int16_t) * 4);
        memcpy(ctu_dec->lfnst_subblock, cg_coeffs, sizeof(int16_t) * 16);

        return sig_sb_map;
    }

    sb_pos    = (last_cg_x << 2) + ((last_cg_y << log2_tb_w) << 2);
    cg_offset = (last_cg_x << 2) + (last_cg_y << 2) * VVC_TR_CTX_STRIDE;

    d_cg = last_cg_x + last_cg_y;

    x = last_x & 3;
    y = last_y & 3;

    start_coeff_idx = ff_vvc_diag_scan_4x4_num_cg[x + (y << 2)];

    c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET + cg_offset];
    c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET + cg_offset];
    c_coding_ctx.sum_sig_nbs  = &num_significant[VVC_TR_CTX_OFFSET + cg_offset];

    num_sig_c = ovcabac_read_ae_sb_4x4_first_dpq(cabac_ctx, cg_coeffs,
                                                 &state, start_coeff_idx, d_cg,
                                                 &c_coding_ctx);

    deq_prms.dequant_sb(cg_coeffs, deq_prms.scale, deq_prms.shift);

    memcpy(&_dst[sb_pos + (0)]             , &cg_coeffs[ 0], sizeof(int16_t) * 4);
    memcpy(&_dst[sb_pos + (1 << log2_tb_w)], &cg_coeffs[ 4], sizeof(int16_t) * 4);
    memcpy(&_dst[sb_pos + (2 << log2_tb_w)], &cg_coeffs[ 8], sizeof(int16_t) * 4);
    memcpy(&_dst[sb_pos + (3 << log2_tb_w)], &cg_coeffs[12], sizeof(int16_t) * 4);

    num_cg = cg_idx_2_cg_num[last_cg_x + last_cg_y * ((1 << lim_log2_w) >> 2)];

    num_cg--;

    for(; num_cg > 0; --num_cg){
        int x_cg = scan_cg_x[num_cg];
        int y_cg = scan_cg_y[num_cg];

        uint8_t sig_sb_flg;
        uint8_t sig_sb_blw = !!(((sig_sb_map >> ((y_cg + 1) << 3)) & 0xFF) & (1 << x_cg));
        uint8_t sig_sb_rgt = !!(((sig_sb_map >> ( y_cg      << 3)) & 0xFF) & (1 << (x_cg + 1)));

        sig_sb_flg = ovcabac_read_ae_significant_cg_flag(cabac_ctx, (sig_sb_rgt | sig_sb_blw));

        if(sig_sb_flg){

            memset(cg_coeffs, 0, sizeof(int16_t) * 16);

            sb_pos    = (x_cg << 2) + ((y_cg << log2_tb_w) << 2);
            cg_offset = (x_cg << 2) + (y_cg << 2) * (VVC_TR_CTX_STRIDE);

            sig_sb_map |= 1llu << (x_cg + (y_cg << 3));

            d_cg = x_cg + y_cg;

            c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET + cg_offset];
            c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET + cg_offset];
            c_coding_ctx.sum_sig_nbs  = &num_significant[VVC_TR_CTX_OFFSET + cg_offset];

            num_sig_c += ovcabac_read_ae_sb_4x4_dpq(cabac_ctx, cg_coeffs,
                                                    &state, d_cg, &c_coding_ctx);

            deq_prms.dequant_sb(cg_coeffs, deq_prms.scale, deq_prms.shift);

            memcpy(&_dst[sb_pos + (0)]             , &cg_coeffs[ 0], sizeof(int16_t) * 4);
            memcpy(&_dst[sb_pos + (1 << log2_tb_w)], &cg_coeffs[ 4], sizeof(int16_t) * 4);
            memcpy(&_dst[sb_pos + (2 << log2_tb_w)], &cg_coeffs[ 8], sizeof(int16_t) * 4);
            memcpy(&_dst[sb_pos + (3 << log2_tb_w)], &cg_coeffs[12], sizeof(int16_t) * 4);
        }
    }

    memset(cg_coeffs, 0, sizeof(uint16_t) * 16);

    c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET];
    c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET];
    c_coding_ctx.sum_sig_nbs  = &num_significant[VVC_TR_CTX_OFFSET];

    num_sig_c += ovcabac_read_ae_sb_4x4_last_dc_dpq(cabac_ctx, cg_coeffs, &state,
                                                    &c_coding_ctx);

    deq_prms.dequant_sb(cg_coeffs, deq_prms.scale, deq_prms.shift);

    memcpy(&_dst[0]             , &cg_coeffs[ 0], sizeof(int16_t) * 4);
    memcpy(&_dst[1 << log2_tb_w], &cg_coeffs[ 4], sizeof(int16_t) * 4);
    memcpy(&_dst[2 << log2_tb_w], &cg_coeffs[ 8], sizeof(int16_t) * 4);
    memcpy(&_dst[3 << log2_tb_w], &cg_coeffs[12], sizeof(int16_t) * 4);

    return sig_sb_map | 1;
}

int
residual_coding_chroma_dpq(OVCTUDec *const ctu_dec, int16_t *const dst,
                           unsigned int log2_tb_w, unsigned int log2_tb_h,
                           uint16_t last_pos)
{
    OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;
    int16_t *const _dst = dst;
    uint8_t tb_width  = 1 << log2_tb_w;
    //check for dependent quantization
    int state = 0;
    uint16_t max_num_bins = ((1 << (log2_tb_w + log2_tb_h) << 5)
                          - ((1 << (log2_tb_w + log2_tb_h)) << 2)) >> 4;
    int last_x, last_y;
    int last_cg_x, last_cg_y;
    int num_cg;
    int x, y;
    int cg_offset;
    int sb_pos;

    //TODO derive start coeff idx in cg
    int start_coeff_idx;
    int16_t cg_coeffs[16] = {0}; //temporary table to store coeffs in process

    uint8_t  num_significant[VVC_TR_CTX_SIZE];
    uint8_t  sum_abs_level  [VVC_TR_CTX_SIZE];
    uint8_t sum_abs_level2 [VVC_TR_CTX_SIZE];

    VVCCoeffCodingCtx c_coding_ctx = {
        .sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET],
        .sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET],
        .sum_sig_nbs  = &num_significant[VVC_TR_CTX_OFFSET],
        .num_remaining_bins = max_num_bins,
        .enable_sdh = ctu_dec->enable_sdh
    };

    uint64_t sig_sb_map = 0;

    int tb_width_in_cg  = tb_width  >> 2;

    int qp = ctu_dec->dequant_chroma->qp;

    const struct IQScale deq_prms = derive_dequant_dpq(qp, log2_tb_w, log2_tb_h);

    memset(_dst, 0, sizeof(int16_t) * (1 << (log2_tb_w + log2_tb_h)));

    if (!last_pos){

        ovcabac_read_ae_sb_dc_coeff_c_dpq(cabac_ctx, cg_coeffs);

        deq_prms.dequant_sb(cg_coeffs, deq_prms.scale, deq_prms.shift);
        _dst[0] = cg_coeffs [0];
        memcpy(ctu_dec->lfnst_subblock, cg_coeffs, sizeof(int16_t) * 16);
        return last_pos;
    }

    reset_ctx_buffers(&c_coding_ctx, log2_tb_w, log2_tb_h);

    if (log2_tb_w > 1 && log2_tb_h > 1) {
        const uint8_t *const cg_idx_2_cg_num = ff_vvc_idx_2_num[log2_tb_w - 2]
                                                               [log2_tb_h - 2];

        const uint8_t *const scan_cg_x = ff_vvc_scan_x[log2_tb_w - 2][log2_tb_h - 2];
        const uint8_t *const scan_cg_y = ff_vvc_scan_y[log2_tb_w - 2][log2_tb_h - 2];

        int num_sig_coeff;

        last_x =  last_pos       & 0x1F;
        last_y = (last_pos >> 8) & 0x1F;

        last_cg_x = last_x >> 2;
        last_cg_y = last_y >> 2;

        num_cg = cg_idx_2_cg_num[last_cg_x + last_cg_y * tb_width_in_cg];

        if(!num_cg){
            int last_coeff_idx = last_x + (last_y << 2);
            int num_coeffs = ff_vvc_diag_scan_4x4_num_cg [last_coeff_idx];

            ovcabac_read_ae_sb_4x4_dc_c_dpq(cabac_ctx, cg_coeffs,
                                            &state, num_coeffs, &c_coding_ctx);

            deq_prms.dequant_sb(cg_coeffs, deq_prms.scale, deq_prms.shift);

            memcpy(&_dst[0]             , &cg_coeffs[ 0], sizeof(int16_t) * 4);
            memcpy(&_dst[1 << log2_tb_w], &cg_coeffs[ 4], sizeof(int16_t) * 4);
            memcpy(&_dst[2 << log2_tb_w], &cg_coeffs[ 8], sizeof(int16_t) * 4);
            memcpy(&_dst[3 << log2_tb_w], &cg_coeffs[12], sizeof(int16_t) * 4);
            memcpy(ctu_dec->lfnst_subblock, cg_coeffs, sizeof(int16_t) * 16);

            return num_coeffs;
        }

        sb_pos = (last_cg_x << 2) + ((last_cg_y * tb_width_in_cg) << 4);

        x = last_x - (last_cg_x << 2);
        y = last_y - (last_cg_y << 2);

        start_coeff_idx = ff_vvc_diag_scan_4x4_num_cg[x + (y << 2)];

        cg_offset = (last_cg_x << 2) + (last_cg_y << 2) * (VVC_TR_CTX_STRIDE);

        c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET + cg_offset];
        c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET + cg_offset];
        c_coding_ctx.sum_sig_nbs  = &num_significant[VVC_TR_CTX_OFFSET + cg_offset];

        num_sig_coeff = ovcabac_read_ae_sb_4x4_first_c_dpq(cabac_ctx, cg_coeffs,
                                                           &state, start_coeff_idx,
                                                           &c_coding_ctx);

        deq_prms.dequant_sb(cg_coeffs, deq_prms.scale, deq_prms.shift);

        memcpy(&_dst[sb_pos + (0)]             , &cg_coeffs[ 0], sizeof(int16_t) * 4);
        memcpy(&_dst[sb_pos + (1 << log2_tb_w)], &cg_coeffs[ 4], sizeof(int16_t) * 4);
        memcpy(&_dst[sb_pos + (2 << log2_tb_w)], &cg_coeffs[ 8], sizeof(int16_t) * 4);
        memcpy(&_dst[sb_pos + (3 << log2_tb_w)], &cg_coeffs[12], sizeof(int16_t) * 4);

        sig_sb_map |= 1llu << (last_cg_x + (last_cg_y << 3));

        num_cg--;

        for(int i = num_cg; i > 0; --i){
            int x_cg = scan_cg_x[i];
            int y_cg = scan_cg_y[i];
            uint8_t sig_sb_blw = !!(((sig_sb_map >> ((y_cg + 1) << 3)) & 0xFF) & (1 << x_cg));
            uint8_t sig_sb_rgt = !!(((sig_sb_map >> ( y_cg      << 3)) & 0xFF) & (1 << (x_cg + 1)));

            uint8_t sig_sb_flg = ovcabac_read_ae_significant_cg_flag_chroma(cabac_ctx, !!(sig_sb_rgt | sig_sb_blw));

            if(sig_sb_flg){

                memset(cg_coeffs, 0, sizeof(int16_t) * 16);

                sb_pos = (x_cg << 2) + ((y_cg * tb_width_in_cg) << 4);
                cg_offset = (x_cg << 2) + (y_cg << 2) * (VVC_TR_CTX_STRIDE);

                sig_sb_map |= 1llu << (x_cg + (y_cg << 3));

                c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET + cg_offset];
                c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET + cg_offset];
                c_coding_ctx.sum_sig_nbs  = &num_significant[VVC_TR_CTX_OFFSET + cg_offset];

                num_sig_coeff += ovcabac_read_ae_sb_4x4_c_dpq(cabac_ctx, cg_coeffs,
                                                              &state, &c_coding_ctx);

                deq_prms.dequant_sb(cg_coeffs, deq_prms.scale, deq_prms.shift);

                memcpy(&_dst[sb_pos + (0)]             , &cg_coeffs[ 0], sizeof(int16_t) * 4);
                memcpy(&_dst[sb_pos + (1 << log2_tb_w)], &cg_coeffs[ 4], sizeof(int16_t) * 4);
                memcpy(&_dst[sb_pos + (2 << log2_tb_w)], &cg_coeffs[ 8], sizeof(int16_t) * 4);
                memcpy(&_dst[sb_pos + (3 << log2_tb_w)], &cg_coeffs[12], sizeof(int16_t) * 4);
            }
        }

        memset(cg_coeffs, 0, sizeof(int16_t) * 16);

        c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET];
        c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET];
        c_coding_ctx.sum_sig_nbs  = &num_significant[VVC_TR_CTX_OFFSET];

        num_sig_coeff += ovcabac_read_ae_sb_4x4_last_dc_c_dpq(cabac_ctx, cg_coeffs,
                                                         &state, &c_coding_ctx);

        deq_prms.dequant_sb(cg_coeffs, deq_prms.scale, deq_prms.shift);

        memcpy(&_dst[0]             , &cg_coeffs[ 0], sizeof(int16_t) * 4);
        memcpy(&_dst[1 << log2_tb_w], &cg_coeffs[ 4], sizeof(int16_t) * 4);
        memcpy(&_dst[2 << log2_tb_w], &cg_coeffs[ 8], sizeof(int16_t) * 4);
        memcpy(&_dst[3 << log2_tb_w], &cg_coeffs[12], sizeof(int16_t) * 4);

       return 0xFFFF;

    } else if (log2_tb_h == 1) {
        int num_sig_c;
        uint8_t sig_sb_flg = 1;

        last_x =  last_pos       & 0x1F;
        last_y = (last_pos >> 8) & 0x1F;

        last_cg_x = last_x >> 3;

        if (!last_cg_x) {
            int last_coeff_idx = last_x + (last_y << 3) ;
            int num_coeffs = ff_vvc_diag_scan_8x2_num_cg [last_coeff_idx];
            num_sig_c = ovcabac_read_ae_sb_8x2_dc_c_dpq(cabac_ctx, cg_coeffs,
                                                        &state, num_coeffs,
                                                        &c_coding_ctx);

            deq_prms.dequant_sb(cg_coeffs, deq_prms.scale, deq_prms.shift);

            memcpy(&_dst[0        ]     , &cg_coeffs[0], sizeof(int16_t) * 8);
            memcpy(&_dst[1 << log2_tb_w], &cg_coeffs[8], sizeof(int16_t) * 8);
            /*FIXME determine whether or not we should read lfnst*/

            return 0xFFFF;
        }

        x = last_x - (last_cg_x << 3);
        y = last_y;

        start_coeff_idx =  ff_vvc_diag_scan_8x2_num_cg [x + y * (1 << 3)];

        sb_pos    = last_cg_x << 3;
        cg_offset = last_cg_x << 3;

        c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET + cg_offset];
        c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET + cg_offset];
        c_coding_ctx.sum_sig_nbs  = &num_significant[VVC_TR_CTX_OFFSET + cg_offset];

        num_sig_c = ovcabac_read_ae_sb_8x2_first_c_dpq(cabac_ctx, cg_coeffs,
                                                       &state, start_coeff_idx,
                                                       &c_coding_ctx);

        deq_prms.dequant_sb(cg_coeffs, deq_prms.scale, deq_prms.shift);

        memcpy(&_dst[sb_pos + 0               ], &cg_coeffs[0], sizeof(int16_t) * 8);
        memcpy(&_dst[sb_pos + (1 << log2_tb_w)], &cg_coeffs[8], sizeof(int16_t) * 8);

        num_cg = last_cg_x;
        num_cg--;

        for ( ;num_cg > 0; --num_cg) {

            sig_sb_flg = ovcabac_read_ae_significant_cg_flag_chroma(cabac_ctx, sig_sb_flg);

            if(sig_sb_flg){

                memset(cg_coeffs, 0, sizeof(int16_t) * 16);

                sb_pos    = num_cg << 3;
                cg_offset = num_cg << 3;

                c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET + cg_offset];
                c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET + cg_offset];
                c_coding_ctx.sum_sig_nbs  = &num_significant[VVC_TR_CTX_OFFSET + cg_offset];

                num_sig_c = ovcabac_read_ae_sb_8x2_c_dpq(cabac_ctx, cg_coeffs, &state,
                                                         &c_coding_ctx);

                deq_prms.dequant_sb(cg_coeffs, deq_prms.scale, deq_prms.shift);

                memcpy(&_dst[sb_pos + 0               ], &cg_coeffs[0], sizeof(int16_t) * 8);
                memcpy(&_dst[sb_pos + (1 << log2_tb_w)], &cg_coeffs[8], sizeof(int16_t) * 8);
            }
        }

        memset(cg_coeffs, 0, sizeof(uint16_t) * 16);

        c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET];
        c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET];
        c_coding_ctx.sum_sig_nbs  = &num_significant[VVC_TR_CTX_OFFSET];

        num_sig_c += ovcabac_read_ae_sb_8x2_last_dc_c_dpq(cabac_ctx, cg_coeffs, &state,
                                                          &c_coding_ctx);

        deq_prms.dequant_sb(cg_coeffs, deq_prms.scale, deq_prms.shift);

        memcpy(&_dst[0             ], &cg_coeffs[0], sizeof(int16_t) * 8);
        memcpy(&_dst[1 << log2_tb_w], &cg_coeffs[8], sizeof(int16_t) * 8);

        return 0xFFFF;

    } else if (log2_tb_w == 1) {
        uint8_t sig_sb_flg = 1;
        int num_sig_c;

        last_x =  last_pos       & 0x1F;
        last_y = (last_pos >> 8) & 0x1F;
        last_cg_y = last_y >> 3;

        if(!last_cg_y){
            int last_coeff_idx = last_x + (last_y << 1);
            int num_coeffs = ff_vvc_diag_scan_2x8_num_cg [last_coeff_idx];
            num_sig_c = ovcabac_read_ae_sb_2x8_dc_c_dpq(cabac_ctx, cg_coeffs,
                                                        &state, num_coeffs,
                                                        &c_coding_ctx);

            deq_prms.dequant_sb(cg_coeffs, deq_prms.scale, deq_prms.shift);

            memcpy(&_dst[0], &cg_coeffs[0],  sizeof(int16_t) * 16);
            /*FIXME determine whether or not we should read lfnst*/

            return 0xFFFF;
        }

        sb_pos    =  last_cg_y << 4;
        cg_offset = (last_cg_y << 3) * VVC_TR_CTX_STRIDE;

        num_cg = last_cg_y;

        x = last_x;
        y = last_y - (last_cg_y << 3);

        start_coeff_idx =  ff_vvc_diag_scan_2x8_num_cg [x + (y << 1)];

        c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET + cg_offset];
        c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET + cg_offset];
        c_coding_ctx.sum_sig_nbs  = &num_significant[VVC_TR_CTX_OFFSET + cg_offset];

        num_sig_c =ovcabac_read_ae_sb_2x8_first_c_dpq(cabac_ctx, cg_coeffs,
                                                      &state, start_coeff_idx,
                                                      &c_coding_ctx);

        deq_prms.dequant_sb(cg_coeffs, deq_prms.scale, deq_prms.shift);

        memcpy(&_dst[sb_pos + 0], &cg_coeffs[0], sizeof(int16_t) * 16);

        num_cg--;

        for ( ;num_cg > 0; --num_cg) {

            sig_sb_flg = ovcabac_read_ae_significant_cg_flag_chroma(cabac_ctx, sig_sb_flg);

            if(sig_sb_flg){
                sb_pos    = (num_cg << 4);
                cg_offset = (num_cg << 3) * VVC_TR_CTX_STRIDE;

                memset(cg_coeffs, 0, sizeof(int16_t) * 16);

                c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET + cg_offset];
                c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET + cg_offset];
                c_coding_ctx.sum_sig_nbs  = &num_significant[VVC_TR_CTX_OFFSET + cg_offset];

                num_sig_c = ovcabac_read_ae_sb_2x8_c_dpq(cabac_ctx, cg_coeffs, &state,
                                                         &c_coding_ctx);

                deq_prms.dequant_sb(cg_coeffs, deq_prms.scale, deq_prms.shift);

                memcpy(&_dst[sb_pos], &cg_coeffs[0], sizeof(int16_t) * 16);
            }
        }

        memset(cg_coeffs, 0, sizeof(uint16_t) * 16);

        c_coding_ctx.sum_abs_lvl  = &sum_abs_level[VVC_TR_CTX_OFFSET];
        c_coding_ctx.sum_abs_lvl2 = &sum_abs_level2[VVC_TR_CTX_OFFSET];
        c_coding_ctx.sum_sig_nbs  = &num_significant[VVC_TR_CTX_OFFSET];

        num_sig_c += ovcabac_read_ae_sb_2x8_last_dc_c_dpq(cabac_ctx,
                                                          cg_coeffs, &state,
                                                          &c_coding_ctx);

        deq_prms.dequant_sb(cg_coeffs, deq_prms.scale, deq_prms.shift);

        memcpy(&_dst[0], &cg_coeffs[0],  sizeof(int16_t) * 16);

        return 0xFFFF;
    }
    return 0;
}
