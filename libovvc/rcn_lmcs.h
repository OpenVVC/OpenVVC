#ifndef RCN_LMCS_H
#define RCN_LMCS_H

#include <stdint.h>

struct CTUBitField;
struct RCNFunctions;
struct OVLMCSData;

struct LMCSInfo
{
    uint8_t  lmcs_enabled_flag;
    uint8_t  scale_c_flag;
    uint16_t lmcs_chroma_scale;
    int16_t  lmcs_chroma_scaling_offset;
    uint8_t min_idx;
    uint8_t max_idx;
    struct LMCSLUTs *luts;
};

void rcn_lmcs_compute_chroma_scale(struct LMCSInfo *const lmcs_info,
                                   const struct CTUBitField *const progress_field,
                                   const uint16_t *ctu_data_y, uint8_t x0, uint8_t y0);

void rcn_init_lmcs_function(struct RCNFunctions *rcn_func, uint8_t lmcs_flag);

void rcn_init_lmcs(struct LMCSInfo *lmcs_info, const struct OVLMCSData *const lmcs_data);

#endif //RCN_LMCS_H
