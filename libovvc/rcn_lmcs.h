#ifndef RCN_LMCS_H
#define RCN_LMCS_H

#include <stdint.h>


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


#endif //RCN_LMCS_H
