#ifndef RCN_LMCS_H
#define RCN_LMCS_H

#include "ctudec.h"

void rcn_lmcs_compute_chroma_scale(struct OVCTUDec* ctudec, int x0, int y0);

void rcn_init_lmcs_function(struct RCNFunctions *rcn_func, uint8_t lmcs_flag);

void rcn_init_lmcs(struct LMCSInfo *lmcs_info, const struct OVLMCSData *const lmcs_data);

#endif //RCN_LMCS_H
