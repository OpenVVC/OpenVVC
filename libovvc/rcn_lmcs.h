#ifndef RCN_LMCS_H
#define RCN_LMCS_H

#include "ctudec.h"

void rcn_derive_lmcs_params(uint16_t *const output_pivot, const OVLMCSData *const lmcs);

void rcn_lmcs_compute_chroma_scale(struct OVCTUDec* ctudec, int x0, int y0);

void rcn_lmcs_reshape_luma_blk(uint16_t *_dst, ptrdiff_t stride_dst,
                            uint16_t* lmcs_output_pivot, int width, int height);

void rcn_lmcs_reshape_luma_blk_lut(uint16_t *_dst, ptrdiff_t stride_dst, uint16_t* lmcs_lut_luma,
                            uint16_t* lmcs_output_pivot, int width, int height);

void rcn_lmcs_compute_lut_luma(uint16_t* lmcs_lut_luma, uint16_t* lmcs_output_pivot);

#endif //RCN_LMCS_H
