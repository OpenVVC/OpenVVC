#include <stdint.h>

#include "rcn_lmcs.h"

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
    uint16_t *_src = &ctudec->rcn_ctx.ctu_buff.y[RCN_CTB_PADDING + x0 + (y0 - 1) * RCN_CTB_STRIDE];

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
    _src = &ctudec->rcn_ctx.ctu_buff.y[RCN_CTB_PADDING + x0 - 1 + y0  * RCN_CTB_STRIDE];

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
    while (idx < 16 && luma_avg >= lmcs_info->lmcs_output_pivot[++idx]){
    }

    lmcs_info->lmcs_chroma_scale = (1<<17)/(uint32_t)(lmcs_info->lmcs_output_pivot[idx] - lmcs_info->lmcs_output_pivot[idx-1]);
}