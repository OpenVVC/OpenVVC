#include <stdint.h>
#include <string.h>

#include "ovutils.h"
#include "rcn_lmcs.h"

void 
rcn_derive_lmcs_params(uint16_t *const output_pivot, const OVLMCSData *const lmcs)
{
    int i;
    // uint16_t *const output_pivot = vvc_ctx->lmcs_output_pivot;
    //TODO: why unused ?
    // uint8_t max_bin_idx = (PIC_CODE_CW_BINS - (lmcs->lmcs_delta_max_bin_idx + 1)) & 0x0F;
    // uint8_t min_bin_idx = lmcs->lmcs_min_bin_idx & 0x0F;
    //BITDEPTH: only 10
    int bitdepth_luma = 10;
    uint8_t window_size = (1 << bitdepth_luma) >> 4; //(log2_cw_shift
    int16_t code_words[17] = {0};
    uint16_t   mapped_intervals[17];
    uint16_t unmapped_intervals[17];
    // uint16_t scaled_unmapped_intervals[17];
    uint16_t interval_low = 0;

    memset(output_pivot, 0, sizeof(uint16_t)*16);

    code_words[0]         = 0;
    unmapped_intervals[0] = 0;
    mapped_intervals[0]   = 0;


    for (i = 0; i < 16; ++i){
        code_words[i] = lmcs->lmcs_delta_sign_cw_flag[i] ? -lmcs->lmcs_delta_abs_cw[i] + interval_low
                                                        : lmcs->lmcs_delta_abs_cw[i] + interval_low;
        interval_low = i + 1 < 15 ? window_size:0;
    }

    for (i = 0; i < 16; ++i){
        mapped_intervals[i+1] = code_words[i] + mapped_intervals[i];
        unmapped_intervals[i + 1]  = unmapped_intervals[i] + window_size;
        // scaled_unmapped_intervals[i] = (((int32_t)code_words[i] << FIXED_POINT_PREC)
        //                                 + ( 1 << (BITDEPTH_PREC - 1))) >> BITDEPTH_PREC;
        output_pivot[i] = mapped_intervals[i]/*(1 << (BITDEPTH_PREC + FIXED_POINT_PREC)) / (code_words[i] ? code_words[i] : 1)*/;
    }
    return;
}

void 
rcn_compute_lmcs_chroma_scale(struct OVCTUDec* ctudec, int x0, int y0)
{
    struct LMCSInfo* lmcs_info = &ctudec->lmcs_info;
    uint64_t upper_map = ctudec->rcn_ctx.progress_field.vfield[y0 >> 2];
    uint64_t left_map  = ctudec->rcn_ctx.progress_field.hfield[x0 >> 2];
    uint64_t needed_mask = (1 << 16) - 1;
    uint32_t available_above_mask = ( upper_map & (needed_mask << (63 - 16 - (x0 >> 2)))) >> (63 - 16 - (x0 >> 2));
    uint32_t available_left_mask  = ( left_map & (needed_mask << (63 - 16 - (y0 >> 2)))) >> (63 - 16 - (y0 >> 2));
    // uint16_t *_src = &ctudec->ctu_data_y[RCN_CTB_PADDING + x0 + (y0 - 1) * RCN_CTB_STRIDE];
    uint16_t *_src = &ctudec->rcn_ctx.ctu_buff.y[x0 + (y0 - 1) * RCN_CTB_STRIDE];

    uint32_t luma_sum1 = 0;
    uint32_t luma_sum2 = 0;
    uint32_t luma_sum3 = 0;
    uint32_t luma_sum4 = 0;
    uint8_t num_luma_pu = 0;
    uint32_t log2_num_luma_pu = 0;
    uint32_t luma_avg=512;
    int idx = 1;

    while (available_above_mask){
        luma_sum1 += _src[0];
        luma_sum2 += _src[1];
        luma_sum3 += _src[2];
        luma_sum4 += _src[3];
        _src += 4;
        ++num_luma_pu;
        available_above_mask >>= 1;
    }

    // _src = &ctudec->ctu_data_y[RCN_CTB_PADDING + x0 - 1 + y0  * RCN_CTB_STRIDE];
    _src = &ctudec->rcn_ctx.ctu_buff.y[x0 - 1 + y0  * RCN_CTB_STRIDE];

    while (available_left_mask){
        luma_sum1 += _src[0];
        luma_sum2 += _src[RCN_CTB_STRIDE];
        luma_sum3 += _src[RCN_CTB_STRIDE << 1];
        luma_sum4 += _src[RCN_CTB_STRIDE * 3];
        _src += RCN_CTB_STRIDE << 2;
        ++num_luma_pu;
        available_left_mask >>= 1;
    }

    while(num_luma_pu){
        ++log2_num_luma_pu;
        num_luma_pu >>= 1;
    }
    //decrease of 1 if not zero
    log2_num_luma_pu -= !!log2_num_luma_pu;

    luma_avg = log2_num_luma_pu ? ((luma_sum1 + luma_sum2) + (luma_sum3 + luma_sum4) + (1 << (log2_num_luma_pu + 1))) >> (log2_num_luma_pu + 2) : 512;
    while (idx < 16 && luma_avg >= lmcs_info->lmcs_output_pivot[idx+1]){
        idx++;
    }

    uint32_t size_interval = (uint32_t)(lmcs_info->lmcs_output_pivot[idx+1] - lmcs_info->lmcs_output_pivot[idx]);
    lmcs_info->lmcs_chroma_scale = (size_interval == 0) ? 1 << 11 
                                                    : (1<<17)/ (size_interval + lmcs_info->lmcs_chroma_scaling_offset);
}

void 
rcn_lmcs_reshape_luma_blk(uint8_t *_dst, ptrdiff_t stride_dst,
                            uint16_t* lmcs_output_pivot, int width, int height)
{
    int16_t *dst = (int16_t *)_dst;
    //BITDEPTH: uniquement pour bitdepth 10
    stride_dst /= sizeof(int16_t);

    int bitdepth = 10;
    uint8_t window_size = (1 << bitdepth) >> 4; 
    uint8_t idx = 0;

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++){
// int idxYInv = getPWLIdxInv(lumaSample);
// int invSample = m_inputPivot[idxYInv] + ((m_invScaleCoef[idxYInv] * (lumaSample - m_reshapePivot[idxYInv]) + (1 << (FP_PREC - 1))) >> FP_PREC); 
            idx = 1;
            while (idx < 16 && dst[x] >= lmcs_output_pivot[idx+1]){
                idx++;
            }
            //TODO: rename variables 
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