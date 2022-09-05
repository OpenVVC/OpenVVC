
#ifndef FG_BLEND_STRIPE_H
#define FG_BLEND_STRIPE_H

void asm_fg_blend_stripe_sse4(int16_t *dstSampleOffsetY, int16_t *srcSampleOffsetY, int32_t *grainStripe, uint32_t widthComp, uint32_t blockHeight, uint8_t bitDepth);

#endif 