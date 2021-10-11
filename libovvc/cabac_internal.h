#ifndef CABAC_INTERNAL_H
#define CABAC_INTERNAL_H
#include "ovutils.h"
#include "vcl_cabac.h"

#define NB_CABAC_BITS 16
#define CABAC_MASK ((1 << NB_CABAC_BITS) - 1)

extern const uint8_t lps_renorm_table[64];
extern const uint8_t range_lps_lut[512];

static inline uint8_t
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
        uint32_t tmp_fill = -CABAC_MASK;

        tmp_fill += cabac_ctx->bytestream[0] << 9;
        tmp_fill += cabac_ctx->bytestream[1] << 1;

        cabac_ctx->low_b += tmp_fill << num_bits;

        /* Last read will try to refill CABAC if last read occurs on
         * the bit just before alignment we will refill the CABAC
         * even if it was the last bit to be read in the entry
         * This is why we use <= instead of <.
         * Doing so permits to check for an error at the end
         * of each CTU line based on the position in the entry
         * Note that we can also check at the end of the entry
         * everything has been consumed.
         */
        if (cabac_ctx->bytestream <= cabac_ctx->bytestream_end){
            cabac_ctx->bytestream += NB_CABAC_BITS >> 3;
        } else {
            /* FIXME this permits to check if we needed to refill
             *  after end of entry
             */
#if 1
            cabac_ctx->bytestream = cabac_ctx->bytestream_end + 2;
            //printf("CABAC_EMPTY\n");
#endif
        }
    }
    return symbol_mask & 0x1;
}

static inline uint8_t
ovcabac_bypass_read(OVCABACCtx *const cabac_ctx)
{
  int32_t range, lps_mask;

  cabac_ctx->low_b <<= 1;

  if (!(cabac_ctx->low_b & CABAC_MASK)){
      int num_bits = 0;
      uint32_t tmp_fill = -CABAC_MASK;
      tmp_fill += cabac_ctx->bytestream[0] << 9;
      tmp_fill += cabac_ctx->bytestream[1] << 1;
      cabac_ctx->low_b += tmp_fill << num_bits;
      if (cabac_ctx->bytestream <= cabac_ctx->bytestream_end){
          cabac_ctx->bytestream += NB_CABAC_BITS >> 3;
      } else {
          /* FIXME this permits to check if we needed to refill
           *  after end of entry
           */
#if 1
          cabac_ctx->bytestream = cabac_ctx->bytestream_end + 2;
          //printf("CABAC_EMPTY\n");
#endif
      }
  }

  range = cabac_ctx->range << (NB_CABAC_BITS + 1);

  lps_mask = range - cabac_ctx->low_b - 1;
  lps_mask >>= 31;

  cabac_ctx->low_b -= range & lps_mask;

  return lps_mask & 0x1;
}
#endif



/* FIXME only used by mip_idx */
static inline uint8_t
vvc_get_cabac_truncated(OVCABACCtx *const cabac_ctx, unsigned int max_symbol){
    static const uint8_t threshold_lut[257] = {
        0, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
        5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
        6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
        6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
        6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
        6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
        7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
        7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
        7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
        7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
        7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
        7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
        7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
        7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
        8
    };
    int threshold;
    uint32_t ruiSymbol = 0;
    /* MAX SYMBOL will not be > 16 */
    #if 0
    if( max_symbol > 256 ){
        int thres_val = 1 << 8;
        threshold = 8;
        while( thres_val <= max_symbol ){
            threshold++;
            thres_val <<= 1;
        }
        threshold--;
    }else{
    #endif
        threshold = threshold_lut[max_symbol];
    #if 0
    }
    #endif

    int val = 1 << threshold;
    int b = max_symbol - val;

    while(threshold--){
        ruiSymbol <<= 1;
        ruiSymbol |= ovcabac_bypass_read(cabac_ctx);
    }

    if( ruiSymbol >= val - b ){
        uint32_t uiSymbol;
        uiSymbol = ovcabac_bypass_read(cabac_ctx);
        ruiSymbol <<= 1;
        ruiSymbol += uiSymbol;
        ruiSymbol -= ( val - b );
    }

    return ruiSymbol;
}
