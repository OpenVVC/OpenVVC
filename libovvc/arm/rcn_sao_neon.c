#include <stddef.h>
#include <stdint.h>
#include <arm_neon.h>

#include "dec_structures.h"
#include "rcn_structures.h"


static void
sao_band_filter_0_10_neon(uint8_t* _dst,
                          uint8_t* _src,
                          ptrdiff_t _stride_dst,
                          ptrdiff_t _stride_src,
                          struct SAOParamsCtu* sao,
                          int width,
                          int height,
                          int c_idx)
{
  int y, x;
  int shift = 10 - 5;
  int16_t* sao_offset_val = sao->offset_val[c_idx];
  uint8_t sao_left_class = sao->band_position[c_idx];
  int16x8_t r0, r1, r2, r3;
  int16x8_t x0, x1, x2, x3;
  int16x8_t sao1, sao2, sao3, sao4;
  int16x8_t src0, src2;
  uint16_t* dst = (uint16_t*)_dst;
  uint16_t* src = (uint16_t*)_src;
  ptrdiff_t stride_dst = _stride_dst >> 1;
  ptrdiff_t stride_src = _stride_src >> 1;

  r0 = vdupq_n_s16((sao_left_class)&31);
  r1 = vdupq_n_s16((sao_left_class + 1) & 31);
  r2 = vdupq_n_s16((sao_left_class + 2) & 31);
  r3 = vdupq_n_s16((sao_left_class + 3) & 31);
  sao1 = vdupq_n_s16(sao_offset_val[0]);
  sao2 = vdupq_n_s16(sao_offset_val[1]);
  sao3 = vdupq_n_s16(sao_offset_val[2]);
  sao4 = vdupq_n_s16(sao_offset_val[3]);
  for (y = 0; y < height; y++) {
    for (x = 0; x < width; x += 8) {
      src0 = vld1q_s16((const int16_t *)&src[x]);
      src2 = vshlq_s16(src0, vdupq_n_s16(-shift));
      x0 = (int16x8_t)vceqq_s16(src2, r0);
      x1 = (int16x8_t)vceqq_s16(src2, r1);
      x2 = (int16x8_t)vceqq_s16(src2, r2);
      x3 = (int16x8_t)vceqq_s16(src2, r3);
      x0 = vandq_s16(x0, sao1);
      x1 = vandq_s16(x1, sao2);
      x2 = vandq_s16(x2, sao3);
      x3 = vandq_s16(x3, sao4);
      x0 = vorrq_s16(x0, x1);
      x2 = vorrq_s16(x2, x3);
      x0 = vorrq_s16(x0, x2);
      src0 = vaddq_s16(src0, x0);
      src0 = vmaxq_s16(src0, vdupq_n_s16(0));
      src0 = vminq_s16(src0, vdupq_n_s16(0x03FF));
      vst1q_u16((uint16_t*)&dst[x],(uint16x8_t) src0);
    }
    dst += stride_dst;
    src += stride_src;
  }
}

static void
sao_edge_filter_10_neon(uint8_t* _dst,
                        uint8_t* _src,
                        ptrdiff_t _stride_dst,
                        ptrdiff_t _stride_src,
                        SAOParamsCtu* sao,
                        int width,
                        int height,
                        int c_idx)
{
  int x, y;
  int16_t* sao_offset_val = sao->offset_val[c_idx];
  int eo = sao->eo_class[c_idx];
  const int8_t pos[4][2][2] = {
    { { -1, 0 }, { 1, 0 } },
    { { 0, -1 }, { 0, 1 } },
    { { -1, -1 }, { 1, 1 } },
    { { 1, -1 }, { -1, 1 } },
  };
  int16x8_t x0, x1, x2, x3;
  int16x8_t offset0, offset1, offset2, offset3, offset4;
  int16x8_t r0, r1, r2, r3, r4;
  int16x8_t cmp0, cmp1;

  uint16_t* dst = (uint16_t*)_dst;
  uint16_t* src = (uint16_t*)_src;
  ptrdiff_t stride_dst = _stride_dst >> 1;
  ptrdiff_t stride_src = _stride_src >> 1;
  int a_stride = pos[eo][0][0] + pos[eo][0][1] * stride_src;
  int b_stride = pos[eo][1][0] + pos[eo][1][1] * stride_src;
  offset0 = vdupq_n_s16(sao_offset_val[0]);
  offset1 = vdupq_n_s16(sao_offset_val[1]);
  offset2 = vdupq_n_s16(sao_offset_val[2]);
  offset3 = vdupq_n_s16(sao_offset_val[3]);
  offset4 = vdupq_n_s16(sao_offset_val[4]);
  for (y = 0; y < height; y++) {
    for (x = 0; x < width; x += 8) {
      x0 = vld1q_s16((int16_t*)(src + x));
      cmp0 = vld1q_s16((int16_t*)(src + x + a_stride));
      cmp1 = vld1q_s16((int16_t*)(src + x + b_stride));
      r2 = vminq_s16(x0, cmp0);
      x1 = (int16x8_t)vceqq_s16(cmp0, r2);
      x2 = (int16x8_t)vceqq_s16(x0, r2);
      x1 = vsubq_s16(x2, x1);
      r2 = vminq_s16(x0, cmp1);
      x3 = (int16x8_t)vceqq_s16(cmp1, r2);
      x2 = (int16x8_t)vceqq_s16(x0, r2);
      x3 = (int16x8_t)vsubq_s16(x2, x3);
      x1 = (int16x8_t)vaddq_s16(x1, x3);
      r0 = (int16x8_t)vceqq_s16((int16x8_t)x1, vdupq_n_s16(-2));
      r1 = (int16x8_t)vceqq_s16((int16x8_t)x1, vdupq_n_s16(-1));
      r2 = (int16x8_t)vceqq_s16((int16x8_t)x1, vdupq_n_s16(0));
      r3 = (int16x8_t)vceqq_s16((int16x8_t)x1, vdupq_n_s16(1));
      r4 = (int16x8_t)vceqq_s16((int16x8_t)x1, vdupq_n_s16(2));
      r0 = vandq_s16(r0, offset0);
      r1 = vandq_s16(r1, offset1);
      r2 = vandq_s16(r2, offset2);
      r3 = vandq_s16(r3, offset3);
      r4 = vandq_s16(r4, offset4);
      r0 = vaddq_s16(r0, r1);
      r2 = vaddq_s16(r2, r3);
      r0 = vaddq_s16(r0, r4);
      r0 = vaddq_s16(r0, r2);
      r0 = vaddq_s16(r0, x0);
      r0 = vmaxq_s16(r0, vdupq_n_s16(0));
      r0 = vminq_s16(r0, vdupq_n_s16(0x03FF));
      vst1q_u16((uint16_t*)(dst + x),(uint16x8_t) r0);
    }
    src += stride_src;
    dst += stride_dst;
  }
}

static void
sao_edge_filter_7_10_neon(uint8_t* _dst,
                          uint8_t* _src,
                          ptrdiff_t _stride_dst,
                          ptrdiff_t _stride_src,
                          SAOParamsCtu* sao,
                          int width,
                          int height,
                          int c_idx)
{
  int x, y;
  int16_t* sao_offset_val = sao->offset_val[c_idx];
  int eo = sao->eo_class[c_idx];
  const int8_t pos[4][2][2] = {
    { { -1, 0 }, { 1, 0 } },
    { { 0, -1 }, { 0, 1 } },
    { { -1, -1 }, { 1, 1 } },
    { { 1, -1 }, { -1, 1 } },
  };
  int16x8_t x0, x1, x2, x3;
  int16x8_t offset0, offset1, offset2, offset3, offset4;
  int16x8_t r0, r1, r2, r3, r4;
  int16x8_t cmp0, cmp1;
  uint16_t* dst = (uint16_t*)_dst;
  uint16_t* src = (uint16_t*)_src;
  ptrdiff_t stride_dst = _stride_dst >> 1;
  ptrdiff_t stride_src = _stride_src >> 1;
  int a_stride = pos[eo][0][0] + pos[eo][0][1] * stride_src;
  int b_stride = pos[eo][1][0] + pos[eo][1][1] * stride_src;
  offset0 = vdupq_n_s16(sao_offset_val[0]);
  offset1 = vdupq_n_s16(sao_offset_val[1]);
  offset2 = vdupq_n_s16(sao_offset_val[2]);
  offset3 = vdupq_n_s16(sao_offset_val[3]);
  offset4 = vdupq_n_s16(sao_offset_val[4]);
  for (y = 0; y < height; y++) {
    for (x = 0; x < width - width%8; x += 8) {
      x0 = vld1q_s16((int16_t*)(src + x));
      cmp0 = vld1q_s16((int16_t*)(src + x + a_stride));
      cmp1 = vld1q_s16((int16_t*)(src + x + b_stride));
      r2 = vminq_s16(x0, cmp0);
      x1 = (int16x8_t)vceqq_s16(cmp0, r2);
      x2 = (int16x8_t)vceqq_s16(x0, r2);
      x1 = vsubq_s16(x2, x1);
      r2 = vminq_s16(x0, cmp1);
      x3 = (int16x8_t)vceqq_s16(cmp1, r2);
      x2 = (int16x8_t)vceqq_s16(x0, r2);
      x3 = (int16x8_t)vsubq_s16(x2, x3);
      x1 = (int16x8_t)vaddq_s16(x1, x3);
      r0 = (int16x8_t)vceqq_s16((int16x8_t)x1, vdupq_n_s16(-2));
      r1 = (int16x8_t)vceqq_s16((int16x8_t)x1, vdupq_n_s16(-1));
      r2 = (int16x8_t)vceqq_s16((int16x8_t)x1, vdupq_n_s16(0));
      r3 = (int16x8_t)vceqq_s16((int16x8_t)x1, vdupq_n_s16(1));
      r4 = (int16x8_t)vceqq_s16((int16x8_t)x1, vdupq_n_s16(2));
      r0 = vandq_s16(r0, offset0);
      r1 = vandq_s16(r1, offset1);
      r2 = vandq_s16(r2, offset2);
      r3 = vandq_s16(r3, offset3);
      r4 = vandq_s16(r4, offset4);
      r0 = vaddq_s16(r0, r1);
      r2 = vaddq_s16(r2, r3);
      r0 = vaddq_s16(r0, r4);
      r0 = vaddq_s16(r0, r2);
      r0 = vaddq_s16(r0, x0);
      r0 = vmaxq_s16(r0, vdupq_n_s16(0));
      r0 = vminq_s16(r0, vdupq_n_s16(0x03FF));
      vst1q_u16((uint16_t*)(dst + x),(uint16x8_t) r0);
    }
    x0 = vld1q_s16((int16_t*)(src + x));
    cmp0 = vld1q_s16((int16_t*)(src + x + a_stride));
    cmp1 = vld1q_s16((int16_t*)(src + x + b_stride));
    r2 = vminq_s16(x0, cmp0);
    x1 = (int16x8_t)vceqq_s16(cmp0, r2);
    x2 = (int16x8_t)vceqq_s16(x0, r2);
    x1 = vsubq_s16(x2, x1);
    r2 = vminq_s16(x0, cmp1);
    x3 = (int16x8_t)vceqq_s16(cmp1, r2);
    x2 = (int16x8_t)vceqq_s16(x0, r2);
    x3 = (int16x8_t)vsubq_s16(x2, x3);
    x1 = (int16x8_t)vaddq_s16(x1, x3);
    r0 = (int16x8_t)vceqq_s16((int16x8_t)x1, vdupq_n_s16(-2));
    r1 = (int16x8_t)vceqq_s16((int16x8_t)x1, vdupq_n_s16(-1));
    r2 = (int16x8_t)vceqq_s16((int16x8_t)x1, vdupq_n_s16(0));
    r3 = (int16x8_t)vceqq_s16((int16x8_t)x1, vdupq_n_s16(1));
    r4 = (int16x8_t)vceqq_s16((int16x8_t)x1, vdupq_n_s16(2));
    r0 = vandq_s16(r0, offset0);
    r1 = vandq_s16(r1, offset1);
    r2 = vandq_s16(r2, offset2);
    r3 = vandq_s16(r3, offset3);
    r4 = vandq_s16(r4, offset4);
    r0 = vaddq_s16(r0, r1);
    r2 = vaddq_s16(r2, r3);
    r0 = vaddq_s16(r0, r4);
    r0 = vaddq_s16(r0, r2);

    //mask to remove processing on last element
    r2 = vdupq_n_s16(0xFFFF);
    r2 = vextq_s16(r2, vdupq_n_s16(0), 1);
    r0 = vandq_s16(r0, r2);

    r0 = vaddq_s16(r0, x0);
    r0 = vmaxq_s16(r0, vdupq_n_s16(0));
    r0 = vminq_s16(r0, vdupq_n_s16(0x03FF));
    vst1q_u16((uint16_t*)(dst + x),(uint16x8_t) r0);
    src += stride_src;
    dst += stride_dst;
  }
}

void rcn_init_sao_functions_neon(struct RCNFunctions *const rcn_funcs){
    rcn_funcs->sao.band= &sao_band_filter_0_10_neon;
    rcn_funcs->sao.edge[0]= &sao_edge_filter_7_10_neon;
    rcn_funcs->sao.edge[1]= &sao_edge_filter_10_neon;
}
