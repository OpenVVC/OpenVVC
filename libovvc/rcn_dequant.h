#ifndef RCN_DEQUANT_H
#define RCN_DEQUANT_H
#include <stdint.h>
#include "ovdefs.h"

struct IQScale
{
    int scale;
    int shift;
    void (*dequant_sb)(int16_t *const sb_coeffs, int scale, int shift);
};

struct VVCQPCTX;
void derive_dequant_ctx(OVCTUDec *const ctudec, const struct VVCQPCTX *const qp_ctx,
                        int cu_qp_delta);

#endif
