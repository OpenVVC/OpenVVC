#ifndef CABAC_INTERNAL_H
#define CABAC_INTERNAL_H
#include "vcl_cabac.h"

#define NB_CABAC_BITS 16
#define CABAC_MASK ((1 << NB_CABAC_BITS) - 1)
#define ov_ctz(x) __builtin_ctz(x)

extern const uint8_t lps_renorm_table[64];
extern const uint8_t range_lps_lut[512];

inline uint8_t
ovcabac_ae_read(OVCABACCtx *const cabac_ctx, uint64_t *const cabac_state)
{
    uint16_t state_0 = (*cabac_state) >> 48;
    uint16_t state_1 = (*cabac_state) >> 32;
    uint8_t   rate_0 = (*cabac_state) >> 16;
    uint8_t   rate_1 = (*cabac_state) >>  0;
    int8_t state = (state_0 + state_1) >> 8;
    uint16_t symbol_mask;
    uint32_t range_lps;
    int32_t lps_mask;
    int log2_renorm;

    symbol_mask = (int16_t)state >> 7;
    state ^= symbol_mask;

    #if 0
    range_lps = ((state >> 2) * ((cabac_ctx->range /*& 0x1E0*/) >> 5) >> 1) + 4;
    #else
    range_lps = range_lps_lut[(cabac_ctx->range & 0x1E0) | (state >> 2)];
    #endif
    cabac_ctx->range -= range_lps;

    lps_mask   = (cabac_ctx->range << (NB_CABAC_BITS + 1)) - cabac_ctx->low_b - 1;
    lps_mask >>= 31;

    symbol_mask ^= lps_mask;

    state_0 -= (state_0 >> rate_0) & 0x7FE0;
    state_1 -= (state_1 >> rate_1) & 0x7FFE;
    state_0 += (0x7fffu >> rate_0) & 0x7FE0 & symbol_mask;
    state_1 += (0x7fffu >> rate_1) & 0x7FFE & symbol_mask;

    *cabac_state &= 0xFFFFFFFF;
    *cabac_state |= (uint64_t)state_0 << 48;
    *cabac_state |= (uint64_t)state_1 << 32;

    cabac_ctx->low_b -= (cabac_ctx->range << (NB_CABAC_BITS + 1)) & (lps_mask);
    cabac_ctx->range += (range_lps - cabac_ctx->range)         & (lps_mask);

    log2_renorm = lps_renorm_table[cabac_ctx->range >> 3];

    cabac_ctx->low_b <<= log2_renorm;
    cabac_ctx->range <<= log2_renorm;

    if (!(cabac_ctx->low_b & CABAC_MASK)){
        int num_bits = ov_ctz(cabac_ctx->low_b) - NB_CABAC_BITS;
        int tmp_fill = -CABAC_MASK;

        tmp_fill += cabac_ctx->bytestream[0] << 9;
        tmp_fill += cabac_ctx->bytestream[1] << 1;

        cabac_ctx->low_b += tmp_fill << num_bits;

        if (cabac_ctx->bytestream < cabac_ctx->bytestream_end){
            cabac_ctx->bytestream += NB_CABAC_BITS >> 3;
        }
    }
    return symbol_mask & 0x1;
}

inline uint8_t
ovcabac_bypass_read(OVCABACCtx *const cabac_ctx)
{
  int32_t range, lps_mask;

  cabac_ctx->low_b <<= 1;

    if (!(cabac_ctx->low_b & CABAC_MASK)){
        int num_bits = 0;
        int tmp_fill = -CABAC_MASK;
        tmp_fill += cabac_ctx->bytestream[0] << 9;
        tmp_fill += cabac_ctx->bytestream[1] << 1;
        cabac_ctx->low_b += tmp_fill << num_bits;
      if (cabac_ctx->bytestream < cabac_ctx->bytestream_end){
          cabac_ctx->bytestream += NB_CABAC_BITS >> 3;
      }
  }

  range = cabac_ctx->range << (NB_CABAC_BITS + 1);

  lps_mask = range - cabac_ctx->low_b - 1;
  lps_mask >>= 31;

  cabac_ctx->low_b -= range & lps_mask;

  return lps_mask & 0x1;
}
#endif
