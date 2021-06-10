#include <stdint.h>
#include <string.h>

#include "ovutils.h"
#include "rcn_lmcs.h"

void 
rcn_derive_lmcs_params(struct LMCSInfo *lmcs_info, uint16_t *const output_pivot, const OVLMCSData *const lmcs)
{
    int i;
    //BITDEPTH: only 10
    int bitdepth_luma = 10;
    uint8_t window_size = (1 << bitdepth_luma) >> 4; //(log2_cw_shift
    int16_t code_words[17] = {0};
    uint16_t   mapped_intervals[17];
    uint16_t unmapped_intervals[17];
    uint16_t interval_low =  lmcs->lmcs_min_bin_idx ? 0 : window_size;
    lmcs_info->min_idx = lmcs->lmcs_min_bin_idx;
    lmcs_info->max_idx = 16 - lmcs->lmcs_delta_max_bin_idx;

    memset(output_pivot, 0, sizeof(uint16_t)*16);

    code_words[0]         = 0;
    unmapped_intervals[0] = 0;
    mapped_intervals[0]   = 0;


    for (i = 0; i < 16; ++i){
        code_words[i] = lmcs->lmcs_delta_sign_cw_flag[i] ? -lmcs->lmcs_delta_abs_cw[i] + interval_low
                                                        : lmcs->lmcs_delta_abs_cw[i] + interval_low;
        interval_low = i + 1 >= lmcs->lmcs_min_bin_idx && i + 1 < lmcs_info->max_idx ? window_size : 0;
    }

    for (i = 0; i < 16; ++i){
        mapped_intervals[i+1] = code_words[i] + mapped_intervals[i];
        unmapped_intervals[i + 1]  = unmapped_intervals[i] + window_size;
        output_pivot[i] = mapped_intervals[i];
    }
    return;
}

void 
rcn_lmcs_compute_chroma_scale(struct OVCTUDec* ctudec, int x0, int y0)
{
    struct LMCSInfo* lmcs_info = &ctudec->lmcs_info;
    uint64_t upper_map = ctudec->rcn_ctx.progress_field.hfield[y0 >> 2];
    uint64_t left_map  = ctudec->rcn_ctx.progress_field.vfield[x0 >> 2];
    uint64_t needed_mask = (1 << 16) - 1;
    // uint32_t available_above_mask = ( upper_map & (needed_mask << ((x0 >> 2)+1)) >> ((x0 >> 2)+1));
    uint32_t available_above_mask = ( upper_map & (needed_mask << ((x0 >> 2)+1))) >> ((x0 >> 2)+1);
    // uint32_t available_left_mask  = ( left_map & (needed_mask << ((y0 >> 2)+1)) >> ((y0 >> 2)+1));
    uint32_t available_left_mask  = ( left_map & (needed_mask << ((y0 >> 2)+1)))  >> ((y0 >> 2)+1);
    // uint16_t *_src = &ctudec->ctu_data_y[RCN_CTB_PADDING + x0 + (y0 - 1) * RCN_CTB_STRIDE];
    uint16_t *_src = &ctudec->rcn_ctx.ctu_buff.y[x0 + (y0 - 1) * RCN_CTB_STRIDE];

    uint32_t luma_sum1 = 0;
    uint32_t luma_sum2 = 0;
    uint32_t luma_sum3 = 0;
    uint32_t luma_sum4 = 0;
    uint8_t num_luma_pu = 0;
    uint32_t log2_num_luma_pu = 0;
    uint32_t luma_avg=512;
    int idx = lmcs_info->min_idx;
    
    uint8_t num_luma_pu_above = 0;
    while (available_above_mask){
        luma_sum1 += _src[0];
        luma_sum2 += _src[1];
        luma_sum3 += _src[2];
        luma_sum4 += _src[3];
        _src += 4;
        ++num_luma_pu_above;
        available_above_mask >>= 1;
    }
    //When image border, num_luma_pum may be different from 16.
    if(num_luma_pu_above){
        _src -= 4;
        luma_sum1 += _src[3]*(16-num_luma_pu_above);
        luma_sum2 += _src[3]*(16-num_luma_pu_above);
        luma_sum3 += _src[3]*(16-num_luma_pu_above);
        luma_sum4 += _src[3]*(16-num_luma_pu_above);
        num_luma_pu_above = 16;
    }

    // _src = &ctudec->ctu_data_y[RCN_CTB_PADDING + x0 - 1 + y0  * RCN_CTB_STRIDE];
    _src = &ctudec->rcn_ctx.ctu_buff.y[x0 - 1 + y0  * RCN_CTB_STRIDE];
    uint8_t num_luma_pu_left = 0;
    while (available_left_mask){
        luma_sum1 += _src[0];
        luma_sum2 += _src[RCN_CTB_STRIDE];
        luma_sum3 += _src[RCN_CTB_STRIDE << 1];
        luma_sum4 += _src[RCN_CTB_STRIDE * 3];
        _src += RCN_CTB_STRIDE << 2;
        ++num_luma_pu_left;
        available_left_mask >>= 1;
    }
    //When image border, num_luma_pum may be different from 16.
    if(num_luma_pu_left){
        _src -= RCN_CTB_STRIDE << 2;
        #if 0
        luma_sum1 += _src[0]*(16-num_luma_pu_left);
        luma_sum2 += _src[RCN_CTB_STRIDE]*(16-num_luma_pu_left);
        luma_sum3 += _src[RCN_CTB_STRIDE << 1]*(16-num_luma_pu_left);
        #endif
        luma_sum1 += _src[RCN_CTB_STRIDE * 3]*(16-num_luma_pu_left);
        luma_sum2 += _src[RCN_CTB_STRIDE * 3]*(16-num_luma_pu_left);
        luma_sum3 += _src[RCN_CTB_STRIDE * 3]*(16-num_luma_pu_left);
        luma_sum4 += _src[RCN_CTB_STRIDE * 3]*(16-num_luma_pu_left);
        num_luma_pu_left = 16;
    }

    num_luma_pu = num_luma_pu_left + num_luma_pu_above;
    while(num_luma_pu){
        ++log2_num_luma_pu;
        num_luma_pu >>= 1;
    }
    //decrease of 1 if not zero
    log2_num_luma_pu -= !!log2_num_luma_pu;

    luma_avg = log2_num_luma_pu ? ((luma_sum1 + luma_sum2) + (luma_sum3 + luma_sum4) + (1 << (log2_num_luma_pu + 1))) >> (log2_num_luma_pu + 2) : 512;
    while (idx < lmcs_info->max_idx && luma_avg >= lmcs_info->lmcs_output_pivot[idx+1]){
        idx++;
    }

    uint32_t size_interval = (uint32_t)(lmcs_info->lmcs_output_pivot[idx+1] - lmcs_info->lmcs_output_pivot[idx]);
    lmcs_info->lmcs_chroma_scale = (size_interval == 0) ? 1 << 11 
                                                    : (1<<17)/ (size_interval + lmcs_info->lmcs_chroma_scaling_offset);
}


void 
rcn_lmcs_reshape_luma_blk(uint16_t *dst, ptrdiff_t stride_dst,
                            uint16_t* lmcs_output_pivot, int width, int height)
{
    int bitdepth = 10;
    uint8_t window_size = (1 << bitdepth) >> 4; 
    uint8_t idx = 0;

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++){
            idx = 1;
            while (idx < 16 && dst[x] >= lmcs_output_pivot[idx+1]){
                idx++;
            }

            int16_t map_low  = lmcs_output_pivot[idx];
            int16_t map_high = lmcs_output_pivot[idx+1];
            int16_t orig_low  = (idx) * window_size;
            int16_t orig_high = (idx+1) * window_size;
            int16_t factor    = (idx == 15) ? 0 : (orig_high - orig_low) * (1 << 11) / (map_high - map_low);
            int16_t luma_inv_reshaped = orig_low + (((dst[x] - map_low) * factor + (1 << (10))) >> 11);

            dst[x] = ov_clip_uintp2(luma_inv_reshaped, bitdepth);
        }
        dst += stride_dst;
    }
}

void 
rcn_lmcs_reshape_luma_blk_lut(uint16_t *dst, ptrdiff_t stride_dst, uint16_t* lmcs_lut_luma, int width, int height)
{
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++){ 
            dst[x] = lmcs_lut_luma[dst[x]];
        }
        dst += stride_dst;
    }
}

void 
rcn_lmcs_compute_lut_luma(struct LMCSInfo *lmcs_info, uint16_t* lmcs_lut_inv_luma, uint16_t* lmcs_lut_fwd_luma, 
                            uint16_t* lmcs_output_pivot)
{
    //BITDEPTH: only 10
    int bitdepth = 10;
    uint8_t window_size = (1 << bitdepth) >> 4; 
    uint16_t idx = lmcs_info->min_idx;
    int idx_fwd;
        
    int16_t map_high, map_low;
    int16_t orig_low, orig_high;
    int16_t factor_inv, luma_inv_reshaped;
    int16_t factor_fwd, luma_fwd_reshaped;

    for (uint16_t val = 0; val < (1<<bitdepth); val++){
        if (idx < lmcs_info->max_idx && val >= lmcs_output_pivot[idx+1]){
            idx++;
        }
        map_low  = lmcs_output_pivot[idx];
        map_high = lmcs_output_pivot[idx+1];
        orig_low  = (idx) * window_size;
        orig_high = (idx+1) * window_size;
        factor_inv    = (idx == 15) ? 0 : (orig_high - orig_low) * (1 << 11) / (map_high - map_low);
        luma_inv_reshaped = orig_low + (((val - map_low) * factor_inv + (1 << bitdepth)) >> 11);
        lmcs_lut_inv_luma[val] = ov_clip_uintp2(luma_inv_reshaped, bitdepth);
        
        idx_fwd = val / window_size;
        map_low  = lmcs_output_pivot[idx_fwd];
        map_high = lmcs_output_pivot[idx_fwd+1];
        orig_low  = (idx_fwd) * window_size;
        orig_high = (idx_fwd+1) * window_size;
        factor_fwd    = (idx_fwd == 15) ? 0 : (map_high - map_low) * (1 << 11) / (orig_high - orig_low);
        luma_fwd_reshaped = map_low + (((val - orig_low) * factor_fwd + (1 << bitdepth)) >> 11);
        lmcs_lut_fwd_luma[val] = ov_clip_uintp2(luma_fwd_reshaped, bitdepth);
    }
}


void rcn_lmcs_no_reshape(uint16_t *dst, ptrdiff_t stride_dst, uint16_t* lmcs_lut_luma, int width, int height){
    return;
}

void rcn_init_lmcs_function(struct RCNFunctions *rcn_func, uint8_t lmcs_flag){
    if(lmcs_flag){
        rcn_func->lmcs_reshape = &rcn_lmcs_reshape_luma_blk_lut;
    }
    else{
        rcn_func->lmcs_reshape = &rcn_lmcs_no_reshape;
    }
}