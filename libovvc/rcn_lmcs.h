#ifndef RCN_LMCS_H
#define RCN_LMCS_H

#include "ctudec.h"

void rcn_derive_lmcs_params(uint16_t *const output_pivot, const OVLMCSData *const lmcs);

void rcn_compute_lmcs_chroma_scale(struct OVCTUDec* ctudec, int x0, int y0);

void rcn_lmcs_reshape_luma_blk(uint8_t *_dst, ptrdiff_t stride_dst,
                            uint16_t* lmcs_output_pivot, int width, int height);
#endif //RCN_LMCS_H
