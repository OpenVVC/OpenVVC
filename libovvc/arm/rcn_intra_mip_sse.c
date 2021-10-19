#include <stdint.h>
#include "rcn_neon.h"

#include "simde/x86/sse4.2.h"
#include "simde/x86/mmx.h"
#include "simde/x86/sse.h"
#include "simde/x86/sse2.h"
#include "simde/x86/sse3.h"
#include "simde/x86/ssse3.h"
#include "simde/x86/sse4.2.h"
#include "simde/x86/avx.h"
#include "simde/x86/avx2.h"
#include "simde/x86/avx512.h"

#include "rcn_structures.h"
#include "ovutils.h"

#define MIP_SHIFT 6

static inline void
mip_matmult_8_8(int16_t * bndy_line, uint16_t * mip_pred, const int stride_x, const uint8_t *matrix_mip, int16_t input_offset, const int rnd_mip, uint8_t log2_red_h)
  {
    uint8_t y;
    __m128i d[4], a[8], b, m[8], tmp[8];
    __m128i rnd = _mm_set1_epi32(rnd_mip);
    __m128i add = _mm_set1_epi32(input_offset);
    b = _mm_load_si128((__m128i *)bndy_line);
      for (y = 0; y < (1 << log2_red_h); y++) {
          d[0] = _mm_load_si128((__m128i *)&matrix_mip[0 * stride_x]);
          d[1] = _mm_load_si128((__m128i *)&matrix_mip[2 * stride_x]);
          d[2] = _mm_load_si128((__m128i *)&matrix_mip[4 * stride_x]);
          d[3] = _mm_load_si128((__m128i *)&matrix_mip[6 * stride_x]);

          tmp[0] = _mm_unpacklo_epi8(d[0], _mm_setzero_si128());
          tmp[1] = _mm_unpackhi_epi8(d[0], _mm_setzero_si128());
          tmp[2] = _mm_unpacklo_epi8(d[1], _mm_setzero_si128());
          tmp[3] = _mm_unpackhi_epi8(d[1], _mm_setzero_si128());

          tmp[4] = _mm_unpacklo_epi8(d[2], _mm_setzero_si128());
          tmp[5] = _mm_unpackhi_epi8(d[2], _mm_setzero_si128());
          tmp[6] = _mm_unpacklo_epi8(d[3], _mm_setzero_si128());
          tmp[7] = _mm_unpackhi_epi8(d[3], _mm_setzero_si128());

          m[0] = _mm_madd_epi16(b, tmp[0]);
          m[1] = _mm_madd_epi16(b, tmp[1]);
          m[2] = _mm_madd_epi16(b, tmp[2]);
          m[3] = _mm_madd_epi16(b, tmp[3]);
          m[4] = _mm_madd_epi16(b, tmp[4]);
          m[5] = _mm_madd_epi16(b, tmp[5]);
          m[6] = _mm_madd_epi16(b, tmp[6]);
          m[7] = _mm_madd_epi16(b, tmp[7]);

          a[0] = _mm_unpacklo_epi32(m[0], m[1]);
          a[1] = _mm_unpacklo_epi32(m[2], m[3]);
          a[2] = _mm_unpacklo_epi32(m[4], m[5]);
          a[3] = _mm_unpacklo_epi32(m[6], m[7]);

          a[4] = _mm_unpackhi_epi32(m[0], m[1]);
          a[5] = _mm_unpackhi_epi32(m[2], m[3]);
          a[6] = _mm_unpackhi_epi32(m[4], m[5]);
          a[7] = _mm_unpackhi_epi32(m[6], m[7]);

          m[0] = _mm_add_epi32(a[0], a[4]);
          m[1] = _mm_add_epi32(a[1], a[5]);
          m[2] = _mm_add_epi32(a[2], a[6]);
          m[3] = _mm_add_epi32(a[3], a[7]);

          a[0] = _mm_unpacklo_epi64(m[0], m[1]);
          a[1] = _mm_unpacklo_epi64(m[2], m[3]);

          a[2] = _mm_unpackhi_epi64(m[0], m[1]);
          a[3] = _mm_unpackhi_epi64(m[2], m[3]);

          m[0] = _mm_add_epi32(a[0], a[2]);
          m[1] = _mm_add_epi32(a[1], a[3]);

          m[0] = _mm_add_epi32(m[0], rnd);
          m[1] = _mm_add_epi32(m[1], rnd);

          m[0] = _mm_srai_epi32(m[0], MIP_SHIFT);
          m[1] = _mm_srai_epi32(m[1], MIP_SHIFT);

          m[0] = _mm_add_epi32(m[0], add);
          m[1] = _mm_add_epi32(m[1], add);

          m[0] = _mm_packs_epi32(m[0], m[1]);

          m[0] = _mm_min_epi16(m[0], _mm_set1_epi16(1023));
          m[0] = _mm_max_epi16(m[0], _mm_setzero_si128());

          _mm_store_si128((__m128i *)mip_pred, m[0]);
          mip_pred += 8;
          matrix_mip += 8*stride_x;
      }
  }

  static inline void
  mip_matmult_8_4(int16_t * bndy_line, uint16_t * mip_pred, const int stride_x, const uint8_t *matrix_mip, int16_t input_offset, const int rnd_mip, uint8_t log2_red_h)
    {
      uint8_t y;
      __m128i d[2], a[4], b, m[4], tmp[4];
      __m128i rnd = _mm_set1_epi32(rnd_mip);
      __m128i add = _mm_set1_epi32(input_offset);
      b = _mm_loadu_si128((__m128i *)bndy_line);
        for (y = 0; y < (1 << log2_red_h); y++) {
            d[0] = _mm_loadu_si128((__m128i *)&matrix_mip[0 * stride_x]);
            d[1] = _mm_loadu_si128((__m128i *)&matrix_mip[2 * stride_x]);

            tmp[0] = _mm_unpacklo_epi8(d[0], _mm_setzero_si128());
            tmp[1] = _mm_unpackhi_epi8(d[0], _mm_setzero_si128());
            tmp[2] = _mm_unpacklo_epi8(d[1], _mm_setzero_si128());
            tmp[3] = _mm_unpackhi_epi8(d[1], _mm_setzero_si128());

            m[0] = _mm_madd_epi16(b, tmp[0]);
            m[1] = _mm_madd_epi16(b, tmp[1]);
            m[2] = _mm_madd_epi16(b, tmp[2]);
            m[3] = _mm_madd_epi16(b, tmp[3]);

            a[0] = _mm_unpacklo_epi32(m[0], m[1]);
            a[1] = _mm_unpacklo_epi32(m[2], m[3]);

            a[2] = _mm_unpackhi_epi32(m[0], m[1]);
            a[3] = _mm_unpackhi_epi32(m[2], m[3]);

            m[0] = _mm_add_epi32(a[0], a[2]);
            m[1] = _mm_add_epi32(a[1], a[3]);

            a[0] = _mm_unpacklo_epi64(m[0], m[1]);
            a[1] = _mm_unpackhi_epi64(m[0], m[1]);

            m[0] = _mm_add_epi32(a[0], a[1]);

            m[0] = _mm_add_epi32(m[0], rnd);

            m[0] = _mm_srai_epi32(m[0], MIP_SHIFT);

            m[0] = _mm_add_epi32(m[0], add);

            m[0] = _mm_packs_epi32(m[0], m[1]);

            m[0] = _mm_min_epi16(m[0], _mm_set1_epi16(1023));
            m[0] = _mm_max_epi16(m[0], _mm_setzero_si128());

            _mm_storeu_si128((__m128i *)mip_pred, m[0]);
            mip_pred += 4;
            matrix_mip += 4*stride_x;
        }
    }

static inline void
mip_matmult_4_4(int16_t * bndy_line, uint16_t * mip_pred, const int stride_x, const uint8_t *matrix_mip, int16_t input_offset, const int rnd_mip, uint8_t log2_red_h)
  {
    uint8_t y;
    __m128i d, b, a, m[2], tmp[2];
    __m128i rnd = _mm_set1_epi32(rnd_mip);
    __m128i add = _mm_set1_epi32(input_offset);
    b = _mm_loadu_si128((__m128i *)bndy_line);
    b = _mm_unpacklo_epi64(b,b);
      for (y = 0; y < (1 << log2_red_h); y++) {
          d = _mm_loadu_si128((__m128i *)&matrix_mip[0 * stride_x]);

          tmp[0] = _mm_unpacklo_epi8(d, _mm_setzero_si128());
          tmp[1] = _mm_unpackhi_epi8(d, _mm_setzero_si128());

          m[0] = _mm_madd_epi16(b, tmp[0]);
          m[1] = _mm_madd_epi16(b, tmp[1]);


          a = _mm_hadd_epi32(m[0], m[1]);


          m[0] = _mm_add_epi32(a, rnd);

          m[0] = _mm_srai_epi32(m[0], MIP_SHIFT);

          m[0] = _mm_add_epi32(m[0], add);

          m[0] = _mm_packs_epi32(m[0], _mm_setzero_si128());

          m[0] = _mm_min_epi16(m[0], _mm_set1_epi16(1023));
          m[0] = _mm_max_epi16(m[0], _mm_setzero_si128());

          _mm_storel_epi64((__m128i *)mip_pred, m[0]);
          mip_pred += 4;
          matrix_mip += 4*stride_x;
      }
  }

static inline void
mip_matmult(int16_t * bndy_line, uint16_t * mip_pred, const int stride_x, const uint8_t *matrix_mip, int16_t input_offset, const int rnd_mip,
  uint8_t log2_bndy, uint8_t log2_red_w, uint8_t log2_red_h)
  {
    if (log2_bndy-1){
      if (log2_red_w == 3){
        mip_matmult_8_8(bndy_line, mip_pred, stride_x, matrix_mip, input_offset, rnd_mip, log2_red_h);
      }
      else if (log2_red_w == 2){
        mip_matmult_8_4(bndy_line, mip_pred, stride_x, matrix_mip, input_offset, rnd_mip, log2_red_h);
      }
    }
    else{
      mip_matmult_4_4(bndy_line, mip_pred, stride_x, matrix_mip, input_offset, rnd_mip, log2_red_h);
    }
}
static inline void
up_sample_h_4_2(uint16_t *const dst, const int16_t *const src,
          const uint16_t *ref,
          int log2_upsampled_size_src, int log2_opposite_size,
          int src_step, int src_stride,
          int dst_step, int dst_stride,
          int ref_step, int log2_scale)
{
    const int16_t *src_line   = src;
    const uint16_t *bndy_line = ref + ref_step;
    uint16_t *dst_line  = dst;
    __m128i a, b, c, r[4];
    __m128i rounding_offset = _mm_set1_epi16(1 << (log2_scale - 1));
    for (int i = 0; i < (1 << log2_opposite_size); i+=2) {
        int16_t *curr_dst = (int16_t*)dst_line;

        a = _mm_loadu_si128((__m128i *)src_line);
        b = _mm_slli_epi64(a, 16);
        c = _mm_set_epi64x (bndy_line[ref_step], bndy_line[0]);
        b = _mm_or_si128(b, c);

        r[0] = _mm_unpacklo_epi16(b,a);
        r[1] = _mm_unpackhi_epi16(b,a);
        r[2] = _mm_unpacklo_epi16(a,a);
        r[3] = _mm_unpackhi_epi16(a,a);

        r[0] = _mm_add_epi16(r[0], r[2]);
        r[1] = _mm_add_epi16(r[1], r[3]);

        r[0] = _mm_add_epi16(r[0], rounding_offset);
        r[1] = _mm_add_epi16(r[1], rounding_offset);

        r[0] = _mm_srai_epi16(r[0], log2_scale);
        r[1] = _mm_srai_epi16(r[1], log2_scale);

        _mm_storeu_si128((__m128i *)curr_dst      , r[0]);
        _mm_storeu_si128((__m128i *)(curr_dst + dst_stride), r[1]);

        src_line  += 2*src_stride;
        dst_line  += 2*dst_stride;
        bndy_line += 2*ref_step;
    }
}

static inline void
up_sample_h_8_2(uint16_t *const dst, const int16_t *const src,
          const uint16_t *ref,
          int log2_upsampled_size_src, int log2_opposite_size,
          int src_step, int src_stride,
          int dst_step, int dst_stride,
          int ref_step, int log2_scale)
{
    const int16_t *src_line   = src;
    const uint16_t *bndy_line = ref + ref_step;
    uint16_t *dst_line  = dst;
    __m128i a, b, c, r[4];
    __m128i rounding_offset = _mm_set1_epi16(1 << (log2_scale - 1));
    for (int i = 0; i < (1 << log2_opposite_size); i++) {
        int16_t *curr_dst = (int16_t*)dst_line;

        a = _mm_loadu_si128((__m128i *)src_line);
        b = _mm_bslli_si128(a, 2);
        c = _mm_set_epi64x (0, bndy_line[0]);
        b = _mm_or_si128(b, c);

        r[0] = _mm_unpacklo_epi16(b,a);
        r[1] = _mm_unpackhi_epi16(b,a);
        r[2] = _mm_unpacklo_epi16(a,a);
        r[3] = _mm_unpackhi_epi16(a,a);

        r[0] = _mm_add_epi16(r[0], r[2]);
        r[1] = _mm_add_epi16(r[1], r[3]);

        r[0] = _mm_add_epi16(r[0], rounding_offset);
        r[1] = _mm_add_epi16(r[1], rounding_offset);

        r[0] = _mm_srai_epi16(r[0], log2_scale);
        r[1] = _mm_srai_epi16(r[1], log2_scale);

        _mm_storeu_si128((__m128i *)curr_dst      , r[0]);
        _mm_storeu_si128((__m128i *)(curr_dst + 8), r[1]);

        src_line  += src_stride;
        dst_line  += dst_stride;
        bndy_line += ref_step;
    }
}

static inline void
up_sample_h_4_4(uint16_t *const dst, const int16_t *const src,
          const uint16_t *ref,
          int log2_upsampled_size_src, int log2_opposite_size,
          int src_step, int src_stride,
          int dst_step, int dst_stride,
          int ref_step, int log2_scale)
{
    const int16_t *src_line   = src;
    const uint16_t *bndy_line = ref + ref_step;
    uint16_t *dst_line  = dst;
    __m128i a, b, c, u[8], tmp[5], r[4];
    __m128i rounding_offset = _mm_set1_epi16(1 << (log2_scale - 1));
    for (int i = 0; i < (1 << log2_opposite_size); i+=2) {
        int16_t *curr_dst = (int16_t*)dst_line;

        a = _mm_loadu_si128((__m128i *)src_line);
        b = _mm_slli_epi64(a, 16);
        c = _mm_set_epi64x (bndy_line[ref_step], bndy_line[0]);
        b = _mm_or_si128(b, c);

        u[0] = _mm_unpacklo_epi16(b,b);
        u[2] = _mm_unpackhi_epi16(b,b);
        u[4] = _mm_unpacklo_epi16(a,a);
        u[6] = _mm_unpackhi_epi16(a,a);

        u[1] = _mm_unpackhi_epi32(u[0],u[0]);
        u[0] = _mm_unpacklo_epi32(u[0],u[0]);
        u[3] = _mm_unpackhi_epi32(u[2],u[2]);
        u[2] = _mm_unpacklo_epi32(u[2],u[2]);
        u[5] = _mm_unpackhi_epi32(u[4],u[4]);
        u[4] = _mm_unpacklo_epi32(u[4],u[4]);
        u[7] = _mm_unpackhi_epi32(u[6],u[6]);
        u[6] = _mm_unpacklo_epi32(u[6],u[6]);

        for (int8_t j = 0; j < 4; j++) {
          tmp[0] = _mm_and_si128(u[j], _mm_set1_epi32(0x0000FFFF));
          tmp[1] = _mm_slli_epi16(_mm_and_si128(u[j], _mm_set1_epi64x(0xFFFFFFFF)), 1);
          tmp[2] = _mm_and_si128(u[j + 4], _mm_set1_epi32(0x0000FFFF));
          tmp[3] = _mm_slli_epi16(_mm_and_si128(u[j + 4], _mm_set1_epi32(0xFFFF0000)), 1);
          tmp[4] = _mm_slli_epi16(_mm_and_si128(u[j + 4], _mm_set1_epi64x(0xFFFFFFFF00000000)), 1);

          tmp[0] = _mm_add_epi16(tmp[0], tmp[1]);
          tmp[1] = _mm_add_epi16(tmp[2], tmp[3]);
          tmp[2] = _mm_add_epi16(tmp[4], rounding_offset);

          tmp[0] = _mm_add_epi16(tmp[0], tmp[1]);
          tmp[0] = _mm_add_epi16(tmp[0], tmp[2]);

          r[j] = _mm_srai_epi16(tmp[0], log2_scale);
        }


        _mm_storeu_si128((__m128i *)curr_dst      , r[0]);
        _mm_storeu_si128((__m128i *)(curr_dst + 8), r[1]);
        _mm_storeu_si128((__m128i *)(curr_dst + dst_stride)    , r[2]);
        _mm_storeu_si128((__m128i *)(curr_dst + dst_stride + 8), r[3]);

        src_line  += 2*src_stride;
        dst_line  += 2*dst_stride;
        bndy_line += 2*ref_step;
    }
}

static inline void
up_sample_h_8_4(uint16_t *const dst, const int16_t *const src,
          const uint16_t *ref,
          int log2_upsampled_size_src, int log2_opposite_size,
          int src_step, int src_stride,
          int dst_step, int dst_stride,
          int ref_step, int log2_scale)
{
    const int16_t *src_line   = src;
    const uint16_t *bndy_line = ref + ref_step;
    uint16_t *dst_line  = dst;
    __m128i a, b, c, u[8], tmp[5], r[4];
    __m128i rounding_offset = _mm_set1_epi16(1 << (log2_scale - 1));
    for (int i = 0; i < (1 << log2_opposite_size); i++) {
        int16_t *curr_dst = (int16_t*)dst_line;

        a = _mm_loadu_si128((__m128i *)src_line);
        b = _mm_bslli_si128(a, 2);
        c = _mm_set_epi64x (0, bndy_line[0]);
        b = _mm_or_si128(b, c);

        u[0] = _mm_unpacklo_epi16(b,b);
        u[2] = _mm_unpackhi_epi16(b,b);
        u[4] = _mm_unpacklo_epi16(a,a);
        u[6] = _mm_unpackhi_epi16(a,a);

        u[1] = _mm_unpackhi_epi32(u[0],u[0]);
        u[0] = _mm_unpacklo_epi32(u[0],u[0]);
        u[3] = _mm_unpackhi_epi32(u[2],u[2]);
        u[2] = _mm_unpacklo_epi32(u[2],u[2]);
        u[5] = _mm_unpackhi_epi32(u[4],u[4]);
        u[4] = _mm_unpacklo_epi32(u[4],u[4]);
        u[7] = _mm_unpackhi_epi32(u[6],u[6]);
        u[6] = _mm_unpacklo_epi32(u[6],u[6]);

        for (int8_t j = 0; j < 4; j++) {
          tmp[0] = _mm_and_si128(u[j], _mm_set1_epi32(0x0000FFFF));
          tmp[1] = _mm_slli_epi16(_mm_and_si128(u[j], _mm_set1_epi64x(0xFFFFFFFF)), 1);
          tmp[2] = _mm_and_si128(u[j + 4], _mm_set1_epi32(0x0000FFFF));
          tmp[3] = _mm_slli_epi16(_mm_and_si128(u[j + 4], _mm_set1_epi32(0xFFFF0000)), 1);
          tmp[4] = _mm_slli_epi16(_mm_and_si128(u[j + 4], _mm_set1_epi64x(0xFFFFFFFF00000000)), 1);

          tmp[0] = _mm_add_epi16(tmp[0], tmp[1]);
          tmp[1] = _mm_add_epi16(tmp[2], tmp[3]);
          tmp[2] = _mm_add_epi16(tmp[4], rounding_offset);

          tmp[0] = _mm_add_epi16(tmp[0], tmp[1]);
          tmp[0] = _mm_add_epi16(tmp[0], tmp[2]);

          r[j] = _mm_srai_epi16(tmp[0], log2_scale);
        }


        _mm_storeu_si128((__m128i *)curr_dst       , r[0]);
        _mm_storeu_si128((__m128i *)(curr_dst + 8) , r[1]);
        _mm_storeu_si128((__m128i *)(curr_dst + 16), r[2]);
        _mm_storeu_si128((__m128i *)(curr_dst + 24), r[3]);

        src_line  += src_stride;
        dst_line  += dst_stride;
        bndy_line += ref_step;
    }
}

static inline void
up_sample_h_4_8(uint16_t *const dst, const int16_t *const src,
          const uint16_t *ref,
          int log2_upsampled_size_src, int log2_opposite_size,
          int src_step, int src_stride,
          int dst_step, int dst_stride,
          int ref_step, int log2_scale)
{
    const int16_t *src_line   = src;
    const uint16_t *bndy_line = ref + ref_step;
    uint16_t *dst_line  = dst;
    __m128i u[16], tmp[7], r[8];
    __m128i rounding_offset = _mm_set1_epi16(1 << (log2_scale - 1));
    for (int i = 0; i < (1 << log2_opposite_size); i+=2) {
        int16_t *curr_dst = (int16_t*)dst_line;

        u[0] = _mm_set1_epi16(bndy_line[0]);
        u[1] = _mm_set1_epi16(src_line[0]);
        u[2] = _mm_set1_epi16(src_line[1]);
        u[3] = _mm_set1_epi16(src_line[2]);
        u[4] = _mm_set1_epi16(bndy_line[ref_step]);
        u[5] = _mm_set1_epi16(src_line[4]);
        u[6] = _mm_set1_epi16(src_line[5]);
        u[7] = _mm_set1_epi16(src_line[6]);

        u[8] = _mm_set1_epi16(src_line[0]);
        u[9] = _mm_set1_epi16(src_line[1]);
        u[10] = _mm_set1_epi16(src_line[2]);
        u[11] = _mm_set1_epi16(src_line[3]);
        u[12] = _mm_set1_epi16(src_line[4]);
        u[13] = _mm_set1_epi16(src_line[5]);
        u[14] = _mm_set1_epi16(src_line[6]);
        u[15] = _mm_set1_epi16(src_line[7]);
        for (int8_t j = 0; j < 8; j++) {
          tmp[0] = _mm_and_si128(u[j], _mm_set1_epi32(0x0000FFFF));
          tmp[1] = _mm_slli_epi16(_mm_and_si128(u[j], _mm_set1_epi64x(0xFFFFFFFF)), 1);
          tmp[2] = _mm_slli_epi16(_mm_and_si128(u[j], _mm_set_epi64x(0, 0xFFFFFFFFFFFFFFFF)), 2);

          tmp[3] = _mm_and_si128(u[j + 8], _mm_set1_epi32(0x0000FFFF));
          tmp[4] = _mm_slli_epi16(_mm_and_si128(u[j + 8], _mm_set1_epi32(0xFFFF0000)), 1);
          tmp[5] = _mm_slli_epi16(_mm_and_si128(u[j + 8], _mm_set1_epi64x(0xFFFFFFFF00000000)), 1);
          tmp[6] = _mm_slli_epi16(_mm_and_si128(u[j + 8], _mm_set_epi64x(0xFFFFFFFFFFFFFFFF,0)), 2);

          tmp[0] = _mm_add_epi16(tmp[0], tmp[1]);
          tmp[1] = _mm_add_epi16(tmp[2], tmp[3]);
          tmp[2] = _mm_add_epi16(tmp[4], tmp[5]);
          tmp[3] = _mm_add_epi16(tmp[6], rounding_offset);

          tmp[0] = _mm_add_epi16(tmp[0], tmp[1]);
          tmp[1] = _mm_add_epi16(tmp[2], tmp[3]);

          tmp[0] = _mm_add_epi16(tmp[0], tmp[1]);

          r[j] = _mm_srai_epi16(tmp[0], log2_scale);
        }

        _mm_storeu_si128((__m128i *)curr_dst       , r[0]);
        _mm_storeu_si128((__m128i *)(curr_dst + 8) , r[1]);
        _mm_storeu_si128((__m128i *)(curr_dst + 16), r[2]);
        _mm_storeu_si128((__m128i *)(curr_dst + 24), r[3]);
        _mm_storeu_si128((__m128i *)(curr_dst  + dst_stride)    , r[4]);
        _mm_storeu_si128((__m128i *)(curr_dst + dst_stride + 8) , r[5]);
        _mm_storeu_si128((__m128i *)(curr_dst + dst_stride + 16), r[6]);
        _mm_storeu_si128((__m128i *)(curr_dst + dst_stride + 24), r[7]);

        src_line  += 2*src_stride;
        dst_line  += 2*dst_stride;
        bndy_line += 2*ref_step;
    }
}

static inline void
up_sample_h_8_8(uint16_t *const dst, const int16_t *const src,
          const uint16_t *ref,
          int log2_upsampled_size_src, int log2_opposite_size,
          int src_step, int src_stride,
          int dst_step, int dst_stride,
          int ref_step, int log2_scale)
{
    const int16_t *src_line   = src;
    const uint16_t *bndy_line = ref + ref_step;
    uint16_t *dst_line  = dst;
    __m128i u[16], tmp[7], r[8];
    __m128i rounding_offset = _mm_set1_epi16(1 << (log2_scale - 1));
    for (int i = 0; i < (1 << log2_opposite_size); i++) {
        int16_t *curr_dst = (int16_t*)dst_line;

        u[0] = _mm_set1_epi16(bndy_line[0]);
        u[1] = _mm_set1_epi16(src_line[0]);
        u[2] = _mm_set1_epi16(src_line[1]);
        u[3] = _mm_set1_epi16(src_line[2]);
        u[4] = _mm_set1_epi16(src_line[3]);
        u[5] = _mm_set1_epi16(src_line[4]);
        u[6] = _mm_set1_epi16(src_line[5]);
        u[7] = _mm_set1_epi16(src_line[6]);

        u[8] = _mm_set1_epi16(src_line[0]);
        u[9] = _mm_set1_epi16(src_line[1]);
        u[10] = _mm_set1_epi16(src_line[2]);
        u[11] = _mm_set1_epi16(src_line[3]);
        u[12] = _mm_set1_epi16(src_line[4]);
        u[13] = _mm_set1_epi16(src_line[5]);
        u[14] = _mm_set1_epi16(src_line[6]);
        u[15] = _mm_set1_epi16(src_line[7]);
        for (int8_t j = 0; j < 8; j++) {
          tmp[0] = _mm_and_si128(u[j], _mm_set1_epi32(0x0000FFFF));
          tmp[1] = _mm_slli_epi16(_mm_and_si128(u[j], _mm_set1_epi64x(0xFFFFFFFF)), 1);
          tmp[2] = _mm_slli_epi16(_mm_and_si128(u[j], _mm_set_epi64x(0, 0xFFFFFFFFFFFFFFFF)), 2);

          tmp[3] = _mm_and_si128(u[j + 8], _mm_set1_epi32(0x0000FFFF));
          tmp[4] = _mm_slli_epi16(_mm_and_si128(u[j + 8], _mm_set1_epi32(0xFFFF0000)), 1);
          tmp[5] = _mm_slli_epi16(_mm_and_si128(u[j + 8], _mm_set1_epi64x(0xFFFFFFFF00000000)), 1);
          tmp[6] = _mm_slli_epi16(_mm_and_si128(u[j + 8], _mm_set_epi64x(0xFFFFFFFFFFFFFFFF,0)), 2);

          tmp[0] = _mm_add_epi16(tmp[0], tmp[1]);
          tmp[1] = _mm_add_epi16(tmp[2], tmp[3]);
          tmp[2] = _mm_add_epi16(tmp[4], tmp[5]);
          tmp[3] = _mm_add_epi16(tmp[6], rounding_offset);

          tmp[0] = _mm_add_epi16(tmp[0], tmp[1]);
          tmp[1] = _mm_add_epi16(tmp[2], tmp[3]);

          tmp[0] = _mm_add_epi16(tmp[0], tmp[1]);

          r[j] = _mm_srai_epi16(tmp[0], log2_scale);
        }

        _mm_storeu_si128((__m128i *)curr_dst       , r[0]);
        _mm_storeu_si128((__m128i *)(curr_dst + 8) , r[1]);
        _mm_storeu_si128((__m128i *)(curr_dst + 16), r[2]);
        _mm_storeu_si128((__m128i *)(curr_dst + 24), r[3]);
        _mm_storeu_si128((__m128i *)(curr_dst + 32), r[4]);
        _mm_storeu_si128((__m128i *)(curr_dst + 40), r[5]);
        _mm_storeu_si128((__m128i *)(curr_dst + 48), r[6]);
        _mm_storeu_si128((__m128i *)(curr_dst + 56), r[7]);

        src_line  += src_stride;
        dst_line  += dst_stride;
        bndy_line += ref_step;
    }
}

static inline void
up_sample_v_4_2(uint16_t *const dst, const int16_t *const src,
          const uint16_t *ref,
          int log2_upsampled_size_src, int log2_opposite_size,
          int src_step, int src_stride,
          int dst_step, int dst_stride,
          int ref_step, int log2_scale)
{
    const int16_t *src_line   = src;
    const uint16_t *bndy_line = ref + ref_step;
    uint16_t *dst_line  = dst;
    __m128i s[5], tmp0, tmp1, r[8];
    __m128i rounding_offset = _mm_set1_epi16(1 << (log2_scale - 1));
    s[0] = _mm_loadu_si128((__m128i *)bndy_line);
    s[1] = _mm_loadu_si128((__m128i *)(src_line));
    s[2] = _mm_loadu_si128((__m128i *)(src_line+1*src_step));
    s[3] = _mm_loadu_si128((__m128i *)(src_line+2*src_step));
    s[4] = _mm_loadu_si128((__m128i *)(src_line+3*src_step));

    for (uint8_t i = 0; i < 4; i++) {
      tmp0 = _mm_add_epi16(s[i], s[i+1]);
      tmp1 = _mm_slli_epi16(s[i+1], 1);
      r[2*i] = _mm_srai_epi16(_mm_add_epi16(tmp0, rounding_offset), log2_scale);
      r[2*i+1] = _mm_srai_epi16(_mm_add_epi16(tmp1, rounding_offset), log2_scale);
    }
    _mm_storeu_si128((__m128i *)dst_line, r[0]);
    _mm_storeu_si128((__m128i *)(dst_line + 1 * dst_step), r[1]);
    _mm_storeu_si128((__m128i *)(dst_line + 2 * dst_step), r[2]);
    _mm_storeu_si128((__m128i *)(dst_line + 3 * dst_step), r[3]);
    _mm_storeu_si128((__m128i *)(dst_line + 4 * dst_step), r[4]);
    _mm_storeu_si128((__m128i *)(dst_line + 5 * dst_step), r[5]);
    _mm_storeu_si128((__m128i *)(dst_line + 6 * dst_step), r[6]);
    _mm_storeu_si128((__m128i *)(dst_line + 7 * dst_step), r[7]);
  }

  static inline void
  up_sample_v_4_4(uint16_t *const dst, const int16_t *const src,
            const uint16_t *ref,
            int log2_upsampled_size_src, int log2_opposite_size,
            int src_step, int src_stride,
            int dst_step, int dst_stride,
            int ref_step, int log2_scale)
  {
      const int16_t *src_line   = src;
      const uint16_t *bndy_line = ref + ref_step;
      uint16_t *dst_line  = dst;
      __m128i s[5], tmp0, tmp1, tmp2, tmp3, r[16], a2, b2, ab;
      __m128i rounding_offset = _mm_set1_epi16(1 << (log2_scale - 1));
      s[0] = _mm_loadu_si128((__m128i *)bndy_line);
      s[1] = _mm_loadu_si128((__m128i *)(src_line));
      s[2] = _mm_loadu_si128((__m128i *)(src_line+1*src_step));
      s[3] = _mm_loadu_si128((__m128i *)(src_line+2*src_step));
      s[4] = _mm_loadu_si128((__m128i *)(src_line+3*src_step));

      for (uint8_t i = 0; i < 4; i++) {
        b2 = _mm_slli_epi16(s[i], 1);
        a2 = _mm_slli_epi16(s[i+1], 1);
        ab = _mm_add_epi16(s[i], s[i+1]);
        tmp0 = _mm_add_epi16(b2, ab);
        tmp1 = _mm_add_epi16(b2, a2);
        tmp2 = _mm_add_epi16(a2, ab);
        tmp3 = _mm_slli_epi16(s[i+1], 2);

        r[4*i] = _mm_srai_epi16(_mm_add_epi16(tmp0, rounding_offset), log2_scale);
        r[4*i+1] = _mm_srai_epi16(_mm_add_epi16(tmp1, rounding_offset), log2_scale);
        r[4*i+2] = _mm_srai_epi16(_mm_add_epi16(tmp2, rounding_offset), log2_scale);
        r[4*i+3] = _mm_srai_epi16(_mm_add_epi16(tmp3, rounding_offset), log2_scale);
      }
      _mm_storeu_si128((__m128i *)dst_line, r[0]);
      _mm_storeu_si128((__m128i *)(dst_line + 1  * dst_step), r[1 ]);
      _mm_storeu_si128((__m128i *)(dst_line + 2  * dst_step), r[2 ]);
      _mm_storeu_si128((__m128i *)(dst_line + 3  * dst_step), r[3 ]);
      _mm_storeu_si128((__m128i *)(dst_line + 4  * dst_step), r[4 ]);
      _mm_storeu_si128((__m128i *)(dst_line + 5  * dst_step), r[5 ]);
      _mm_storeu_si128((__m128i *)(dst_line + 6  * dst_step), r[6 ]);
      _mm_storeu_si128((__m128i *)(dst_line + 7  * dst_step), r[7 ]);
      _mm_storeu_si128((__m128i *)(dst_line + 8  * dst_step), r[8 ]);
      _mm_storeu_si128((__m128i *)(dst_line + 9  * dst_step), r[9 ]);
      _mm_storeu_si128((__m128i *)(dst_line + 10 * dst_step), r[10]);
      _mm_storeu_si128((__m128i *)(dst_line + 11 * dst_step), r[11]);
      _mm_storeu_si128((__m128i *)(dst_line + 12 * dst_step), r[12]);
      _mm_storeu_si128((__m128i *)(dst_line + 13 * dst_step), r[13]);
      _mm_storeu_si128((__m128i *)(dst_line + 14 * dst_step), r[14]);
      _mm_storeu_si128((__m128i *)(dst_line + 15 * dst_step), r[15]);
    }

    static inline void
    up_sample_v_4_8(uint16_t *const dst, const int16_t *const src,
              const uint16_t *ref,
              int log2_upsampled_size_src, int log2_opposite_size,
              int src_step, int src_stride,
              int dst_step, int dst_stride,
              int ref_step, int log2_scale)
    {
        const int16_t *src_line   = src;
        const uint16_t *bndy_line = ref + ref_step;
        uint16_t *dst_line  = dst;
        __m128i s[5], tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, r[32], a2, b2, a4, b4, ab, a2b2;
        __m128i rounding_offset = _mm_set1_epi16(1 << (log2_scale - 1));
        s[0] = _mm_loadu_si128((__m128i *)bndy_line);
        s[1] = _mm_loadu_si128((__m128i *)(src_line));
        s[2] = _mm_loadu_si128((__m128i *)(src_line+1*src_step));
        s[3] = _mm_loadu_si128((__m128i *)(src_line+2*src_step));
        s[4] = _mm_loadu_si128((__m128i *)(src_line+3*src_step));

        for (uint8_t i = 0; i < 4; i++) {
          b2 = _mm_slli_epi16(s[i], 1);
          b4 = _mm_slli_epi16(s[i], 2);
          a2 = _mm_slli_epi16(s[i+1], 1);
          a4 = _mm_slli_epi16(s[i+1], 2);
          ab = _mm_add_epi16(s[i], s[i+1]);
          a2b2 = _mm_slli_epi16(ab, 1);
          tmp0 = _mm_add_epi16(b2, ab);
          tmp0 = _mm_add_epi16(tmp0, b4);
          tmp1 = _mm_add_epi16(b4, a2b2);
          tmp2 = _mm_add_epi16(a2, ab);
          tmp2 = _mm_add_epi16(tmp2, b4);
          tmp3 = _mm_slli_epi16(a2b2, 1);
          tmp4 = _mm_add_epi16(b2, ab);
          tmp4 = _mm_add_epi16(tmp4, a4);
          tmp5 = _mm_add_epi16(a4, a2b2);
          tmp6 = _mm_add_epi16(a2, ab);
          tmp6 = _mm_add_epi16(tmp6, a4);
          tmp7 = _mm_slli_epi16(a4, 1);

          r[8*i] = _mm_srai_epi16(_mm_add_epi16(tmp0, rounding_offset), log2_scale);
          r[8*i+1] = _mm_srai_epi16(_mm_add_epi16(tmp1, rounding_offset), log2_scale);
          r[8*i+2] = _mm_srai_epi16(_mm_add_epi16(tmp2, rounding_offset), log2_scale);
          r[8*i+3] = _mm_srai_epi16(_mm_add_epi16(tmp3, rounding_offset), log2_scale);
          r[8*i+4] = _mm_srai_epi16(_mm_add_epi16(tmp4, rounding_offset), log2_scale);
          r[8*i+5] = _mm_srai_epi16(_mm_add_epi16(tmp5, rounding_offset), log2_scale);
          r[8*i+6] = _mm_srai_epi16(_mm_add_epi16(tmp6, rounding_offset), log2_scale);
          r[8*i+7] = _mm_srai_epi16(_mm_add_epi16(tmp7, rounding_offset), log2_scale);
        }
        _mm_storeu_si128((__m128i *)dst_line, r[0]);
        _mm_storeu_si128((__m128i *)(dst_line + 1  * dst_step), r[1 ]);
        _mm_storeu_si128((__m128i *)(dst_line + 2  * dst_step), r[2 ]);
        _mm_storeu_si128((__m128i *)(dst_line + 3  * dst_step), r[3 ]);
        _mm_storeu_si128((__m128i *)(dst_line + 4  * dst_step), r[4 ]);
        _mm_storeu_si128((__m128i *)(dst_line + 5  * dst_step), r[5 ]);
        _mm_storeu_si128((__m128i *)(dst_line + 6  * dst_step), r[6 ]);
        _mm_storeu_si128((__m128i *)(dst_line + 7  * dst_step), r[7 ]);
        _mm_storeu_si128((__m128i *)(dst_line + 8  * dst_step), r[8 ]);
        _mm_storeu_si128((__m128i *)(dst_line + 9  * dst_step), r[9 ]);
        _mm_storeu_si128((__m128i *)(dst_line + 10 * dst_step), r[10]);
        _mm_storeu_si128((__m128i *)(dst_line + 11 * dst_step), r[11]);
        _mm_storeu_si128((__m128i *)(dst_line + 12 * dst_step), r[12]);
        _mm_storeu_si128((__m128i *)(dst_line + 13 * dst_step), r[13]);
        _mm_storeu_si128((__m128i *)(dst_line + 14 * dst_step), r[14]);
        _mm_storeu_si128((__m128i *)(dst_line + 15 * dst_step), r[15]);
        _mm_storeu_si128((__m128i *)(dst_line + 16 * dst_step), r[16]);
        _mm_storeu_si128((__m128i *)(dst_line + 17 * dst_step), r[17]);
        _mm_storeu_si128((__m128i *)(dst_line + 18 * dst_step), r[18]);
        _mm_storeu_si128((__m128i *)(dst_line + 19 * dst_step), r[19]);
        _mm_storeu_si128((__m128i *)(dst_line + 20 * dst_step), r[20]);
        _mm_storeu_si128((__m128i *)(dst_line + 21 * dst_step), r[21]);
        _mm_storeu_si128((__m128i *)(dst_line + 22 * dst_step), r[22]);
        _mm_storeu_si128((__m128i *)(dst_line + 23 * dst_step), r[23]);
        _mm_storeu_si128((__m128i *)(dst_line + 24 * dst_step), r[24]);
        _mm_storeu_si128((__m128i *)(dst_line + 25 * dst_step), r[25]);
        _mm_storeu_si128((__m128i *)(dst_line + 26 * dst_step), r[26]);
        _mm_storeu_si128((__m128i *)(dst_line + 27 * dst_step), r[27]);
        _mm_storeu_si128((__m128i *)(dst_line + 28 * dst_step), r[28]);
        _mm_storeu_si128((__m128i *)(dst_line + 29 * dst_step), r[29]);
        _mm_storeu_si128((__m128i *)(dst_line + 30 * dst_step), r[30]);
        _mm_storeu_si128((__m128i *)(dst_line + 31 * dst_step), r[31]);
      }
  static inline void
  up_sample_v_8_2(uint16_t *const dst, const int16_t *const src,
            const uint16_t *ref,
            int log2_upsampled_size_src, int log2_opposite_size,
            int src_step, int src_stride,
            int dst_step, int dst_stride,
            int ref_step, int log2_scale)
  {
    const int16_t *src_line   = src;
    const uint16_t *bndy_line = ref + ref_step;
    uint16_t *dst_line  = dst;
    for (uint8_t line = 0; line < (1<<(log2_opposite_size)); line+=8) {
      __m128i s[9], tmp0, tmp1, r[16];
      __m128i rounding_offset = _mm_set1_epi16(1 << (log2_scale - 1));
      s[0] = _mm_loadu_si128((__m128i *)bndy_line);
      s[1] = _mm_loadu_si128((__m128i *)(src_line));
      s[2] = _mm_loadu_si128((__m128i *)(src_line+1*src_step));
      s[3] = _mm_loadu_si128((__m128i *)(src_line+2*src_step));
      s[4] = _mm_loadu_si128((__m128i *)(src_line+3*src_step));
      s[5] = _mm_loadu_si128((__m128i *)(src_line+4*src_step));
      s[6] = _mm_loadu_si128((__m128i *)(src_line+5*src_step));
      s[7] = _mm_loadu_si128((__m128i *)(src_line+6*src_step));
      s[8] = _mm_loadu_si128((__m128i *)(src_line+7*src_step));

      for (uint8_t i = 0; i < 8; i++) {
        tmp0 = _mm_add_epi16(s[i], s[i+1]);
        tmp1 = _mm_slli_epi16(s[i+1], 1);
        r[2*i] = _mm_srai_epi16(_mm_add_epi16(tmp0, rounding_offset), log2_scale);
        r[2*i+1] = _mm_srai_epi16(_mm_add_epi16(tmp1, rounding_offset), log2_scale);
      }
      _mm_storeu_si128((__m128i *)dst_line, r[0]);
      _mm_storeu_si128((__m128i *)(dst_line + 1  * dst_step), r[1 ]);
      _mm_storeu_si128((__m128i *)(dst_line + 2  * dst_step), r[2 ]);
      _mm_storeu_si128((__m128i *)(dst_line + 3  * dst_step), r[3 ]);
      _mm_storeu_si128((__m128i *)(dst_line + 4  * dst_step), r[4 ]);
      _mm_storeu_si128((__m128i *)(dst_line + 5  * dst_step), r[5 ]);
      _mm_storeu_si128((__m128i *)(dst_line + 6  * dst_step), r[6 ]);
      _mm_storeu_si128((__m128i *)(dst_line + 7  * dst_step), r[7 ]);
      _mm_storeu_si128((__m128i *)(dst_line + 8  * dst_step), r[8 ]);
      _mm_storeu_si128((__m128i *)(dst_line + 9  * dst_step), r[9 ]);
      _mm_storeu_si128((__m128i *)(dst_line + 10 * dst_step), r[10]);
      _mm_storeu_si128((__m128i *)(dst_line + 11 * dst_step), r[11]);
      _mm_storeu_si128((__m128i *)(dst_line + 12 * dst_step), r[12]);
      _mm_storeu_si128((__m128i *)(dst_line + 13 * dst_step), r[13]);
      _mm_storeu_si128((__m128i *)(dst_line + 14 * dst_step), r[14]);
      _mm_storeu_si128((__m128i *)(dst_line + 15 * dst_step), r[15]);
      src_line  += 8*src_stride;
      dst_line  += 8*dst_stride;
      bndy_line += 8*ref_step;
    }
  }

  static inline void
  up_sample_v_8_4(uint16_t *const dst, const int16_t *const src,
            const uint16_t *ref,
            int log2_upsampled_size_src, int log2_opposite_size,
            int src_step, int src_stride,
            int dst_step, int dst_stride,
            int ref_step, int log2_scale)
  {
    const int16_t *src_line   = src;
    const uint16_t *bndy_line = ref + ref_step;
    uint16_t *dst_line  = dst;
    for (uint8_t line = 0; line < (1<<(log2_opposite_size)); line+=8) {
      __m128i s[9], tmp0, tmp1, tmp2, tmp3, r[32], a2, b2, ab;
      __m128i rounding_offset = _mm_set1_epi16(1 << (log2_scale - 1));
      s[0] = _mm_loadu_si128((__m128i *)bndy_line);
      s[1] = _mm_loadu_si128((__m128i *)(src_line));
      s[2] = _mm_loadu_si128((__m128i *)(src_line+1*src_step));
      s[3] = _mm_loadu_si128((__m128i *)(src_line+2*src_step));
      s[4] = _mm_loadu_si128((__m128i *)(src_line+3*src_step));
      s[5] = _mm_loadu_si128((__m128i *)(src_line+4*src_step));
      s[6] = _mm_loadu_si128((__m128i *)(src_line+5*src_step));
      s[7] = _mm_loadu_si128((__m128i *)(src_line+6*src_step));
      s[8] = _mm_loadu_si128((__m128i *)(src_line+7*src_step));

      for (uint8_t i = 0; i < 8; i++) {
        b2 = _mm_slli_epi16(s[i], 1);
        a2 = _mm_slli_epi16(s[i+1], 1);
        ab = _mm_add_epi16(s[i], s[i+1]);
        tmp0 = _mm_add_epi16(b2, ab);
        tmp1 = _mm_add_epi16(b2, a2);
        tmp2 = _mm_add_epi16(a2, ab);
        tmp3 = _mm_slli_epi16(s[i+1], 2);

        r[4*i] = _mm_srai_epi16(_mm_add_epi16(tmp0, rounding_offset), log2_scale);
        r[4*i+1] = _mm_srai_epi16(_mm_add_epi16(tmp1, rounding_offset), log2_scale);
        r[4*i+2] = _mm_srai_epi16(_mm_add_epi16(tmp2, rounding_offset), log2_scale);
        r[4*i+3] = _mm_srai_epi16(_mm_add_epi16(tmp3, rounding_offset), log2_scale);
      }
      _mm_storeu_si128((__m128i *)dst_line, r[0]);
      _mm_storeu_si128((__m128i *)(dst_line + 1  * dst_step), r[1 ]);
      _mm_storeu_si128((__m128i *)(dst_line + 2  * dst_step), r[2 ]);
      _mm_storeu_si128((__m128i *)(dst_line + 3  * dst_step), r[3 ]);
      _mm_storeu_si128((__m128i *)(dst_line + 4  * dst_step), r[4 ]);
      _mm_storeu_si128((__m128i *)(dst_line + 5  * dst_step), r[5 ]);
      _mm_storeu_si128((__m128i *)(dst_line + 6  * dst_step), r[6 ]);
      _mm_storeu_si128((__m128i *)(dst_line + 7  * dst_step), r[7 ]);
      _mm_storeu_si128((__m128i *)(dst_line + 8  * dst_step), r[8 ]);
      _mm_storeu_si128((__m128i *)(dst_line + 9  * dst_step), r[9 ]);
      _mm_storeu_si128((__m128i *)(dst_line + 10 * dst_step), r[10]);
      _mm_storeu_si128((__m128i *)(dst_line + 11 * dst_step), r[11]);
      _mm_storeu_si128((__m128i *)(dst_line + 12 * dst_step), r[12]);
      _mm_storeu_si128((__m128i *)(dst_line + 13 * dst_step), r[13]);
      _mm_storeu_si128((__m128i *)(dst_line + 14 * dst_step), r[14]);
      _mm_storeu_si128((__m128i *)(dst_line + 15 * dst_step), r[15]);
      _mm_storeu_si128((__m128i *)(dst_line + 16 * dst_step), r[16]);
      _mm_storeu_si128((__m128i *)(dst_line + 17 * dst_step), r[17]);
      _mm_storeu_si128((__m128i *)(dst_line + 18 * dst_step), r[18]);
      _mm_storeu_si128((__m128i *)(dst_line + 19 * dst_step), r[19]);
      _mm_storeu_si128((__m128i *)(dst_line + 20 * dst_step), r[20]);
      _mm_storeu_si128((__m128i *)(dst_line + 21 * dst_step), r[21]);
      _mm_storeu_si128((__m128i *)(dst_line + 22 * dst_step), r[22]);
      _mm_storeu_si128((__m128i *)(dst_line + 23 * dst_step), r[23]);
      _mm_storeu_si128((__m128i *)(dst_line + 24 * dst_step), r[24]);
      _mm_storeu_si128((__m128i *)(dst_line + 25 * dst_step), r[25]);
      _mm_storeu_si128((__m128i *)(dst_line + 26 * dst_step), r[26]);
      _mm_storeu_si128((__m128i *)(dst_line + 27 * dst_step), r[27]);
      _mm_storeu_si128((__m128i *)(dst_line + 28 * dst_step), r[28]);
      _mm_storeu_si128((__m128i *)(dst_line + 29 * dst_step), r[29]);
      _mm_storeu_si128((__m128i *)(dst_line + 30 * dst_step), r[30]);
      _mm_storeu_si128((__m128i *)(dst_line + 31 * dst_step), r[31]);
      src_line  += 8*src_stride;
      dst_line  += 8*dst_stride;
      bndy_line += 8*ref_step;
    }
  }

  static inline void
  up_sample_v_8_8(uint16_t *const dst, const int16_t *const src,
            const uint16_t *ref,
            int log2_upsampled_size_src, int log2_opposite_size,
            int src_step, int src_stride,
            int dst_step, int dst_stride,
            int ref_step, int log2_scale)
  {
    const int16_t *src_line   = src;
    const uint16_t *bndy_line = ref + ref_step;
    uint16_t *dst_line  = dst;
    for (uint8_t line = 0; line < (1<<(log2_opposite_size)); line+=8) {
      __m128i s[9], tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, r[64], a2, b2, a4, b4, ab, a2b2;
      __m128i rounding_offset = _mm_set1_epi16(1 << (log2_scale - 1));
      s[0] = _mm_loadu_si128((__m128i *)bndy_line);
      s[1] = _mm_loadu_si128((__m128i *)(src_line));
      s[2] = _mm_loadu_si128((__m128i *)(src_line+1*src_step));
      s[3] = _mm_loadu_si128((__m128i *)(src_line+2*src_step));
      s[4] = _mm_loadu_si128((__m128i *)(src_line+3*src_step));
      s[5] = _mm_loadu_si128((__m128i *)(src_line+4*src_step));
      s[6] = _mm_loadu_si128((__m128i *)(src_line+5*src_step));
      s[7] = _mm_loadu_si128((__m128i *)(src_line+6*src_step));
      s[8] = _mm_loadu_si128((__m128i *)(src_line+7*src_step));

      for (uint8_t i = 0; i < 8; i++) {
        b2 = _mm_slli_epi16(s[i], 1);
        b4 = _mm_slli_epi16(s[i], 2);
        a2 = _mm_slli_epi16(s[i+1], 1);
        a4 = _mm_slli_epi16(s[i+1], 2);
        ab = _mm_add_epi16(s[i], s[i+1]);
        a2b2 = _mm_slli_epi16(ab, 1);
        tmp0 = _mm_add_epi16(b2, ab);
        tmp0 = _mm_add_epi16(tmp0, b4);
        tmp1 = _mm_add_epi16(b4, a2b2);
        tmp2 = _mm_add_epi16(a2, ab);
        tmp2 = _mm_add_epi16(tmp2, b4);
        tmp3 = _mm_slli_epi16(a2b2, 1);
        tmp4 = _mm_add_epi16(b2, ab);
        tmp4 = _mm_add_epi16(tmp4, a4);
        tmp5 = _mm_add_epi16(a4, a2b2);
        tmp6 = _mm_add_epi16(a2, ab);
        tmp6 = _mm_add_epi16(tmp6, a4);
        tmp7 = _mm_slli_epi16(a4, 1);

        r[8*i] = _mm_srai_epi16(_mm_add_epi16(tmp0, rounding_offset), log2_scale);
        r[8*i+1] = _mm_srai_epi16(_mm_add_epi16(tmp1, rounding_offset), log2_scale);
        r[8*i+2] = _mm_srai_epi16(_mm_add_epi16(tmp2, rounding_offset), log2_scale);
        r[8*i+3] = _mm_srai_epi16(_mm_add_epi16(tmp3, rounding_offset), log2_scale);
        r[8*i+4] = _mm_srai_epi16(_mm_add_epi16(tmp4, rounding_offset), log2_scale);
        r[8*i+5] = _mm_srai_epi16(_mm_add_epi16(tmp5, rounding_offset), log2_scale);
        r[8*i+6] = _mm_srai_epi16(_mm_add_epi16(tmp6, rounding_offset), log2_scale);
        r[8*i+7] = _mm_srai_epi16(_mm_add_epi16(tmp7, rounding_offset), log2_scale);
      }
      _mm_storeu_si128((__m128i *)dst_line, r[0]);
      _mm_storeu_si128((__m128i *)(dst_line + 1  * dst_step), r[1 ]);
      _mm_storeu_si128((__m128i *)(dst_line + 2  * dst_step), r[2 ]);
      _mm_storeu_si128((__m128i *)(dst_line + 3  * dst_step), r[3 ]);
      _mm_storeu_si128((__m128i *)(dst_line + 4  * dst_step), r[4 ]);
      _mm_storeu_si128((__m128i *)(dst_line + 5  * dst_step), r[5 ]);
      _mm_storeu_si128((__m128i *)(dst_line + 6  * dst_step), r[6 ]);
      _mm_storeu_si128((__m128i *)(dst_line + 7  * dst_step), r[7 ]);
      _mm_storeu_si128((__m128i *)(dst_line + 8  * dst_step), r[8 ]);
      _mm_storeu_si128((__m128i *)(dst_line + 9  * dst_step), r[9 ]);
      _mm_storeu_si128((__m128i *)(dst_line + 10 * dst_step), r[10]);
      _mm_storeu_si128((__m128i *)(dst_line + 11 * dst_step), r[11]);
      _mm_storeu_si128((__m128i *)(dst_line + 12 * dst_step), r[12]);
      _mm_storeu_si128((__m128i *)(dst_line + 13 * dst_step), r[13]);
      _mm_storeu_si128((__m128i *)(dst_line + 14 * dst_step), r[14]);
      _mm_storeu_si128((__m128i *)(dst_line + 15 * dst_step), r[15]);
      _mm_storeu_si128((__m128i *)(dst_line + 16 * dst_step), r[16]);
      _mm_storeu_si128((__m128i *)(dst_line + 17 * dst_step), r[17]);
      _mm_storeu_si128((__m128i *)(dst_line + 18 * dst_step), r[18]);
      _mm_storeu_si128((__m128i *)(dst_line + 19 * dst_step), r[19]);
      _mm_storeu_si128((__m128i *)(dst_line + 20 * dst_step), r[20]);
      _mm_storeu_si128((__m128i *)(dst_line + 21 * dst_step), r[21]);
      _mm_storeu_si128((__m128i *)(dst_line + 22 * dst_step), r[22]);
      _mm_storeu_si128((__m128i *)(dst_line + 23 * dst_step), r[23]);
      _mm_storeu_si128((__m128i *)(dst_line + 24 * dst_step), r[24]);
      _mm_storeu_si128((__m128i *)(dst_line + 25 * dst_step), r[25]);
      _mm_storeu_si128((__m128i *)(dst_line + 26 * dst_step), r[26]);
      _mm_storeu_si128((__m128i *)(dst_line + 27 * dst_step), r[27]);
      _mm_storeu_si128((__m128i *)(dst_line + 28 * dst_step), r[28]);
      _mm_storeu_si128((__m128i *)(dst_line + 29 * dst_step), r[29]);
      _mm_storeu_si128((__m128i *)(dst_line + 30 * dst_step), r[30]);
      _mm_storeu_si128((__m128i *)(dst_line + 31 * dst_step), r[31]);
      _mm_storeu_si128((__m128i *)(dst_line + 32 * dst_step), r[32]);
      _mm_storeu_si128((__m128i *)(dst_line + 33 * dst_step), r[33]);
      _mm_storeu_si128((__m128i *)(dst_line + 34 * dst_step), r[34]);
      _mm_storeu_si128((__m128i *)(dst_line + 35 * dst_step), r[35]);
      _mm_storeu_si128((__m128i *)(dst_line + 36 * dst_step), r[36]);
      _mm_storeu_si128((__m128i *)(dst_line + 37 * dst_step), r[37]);
      _mm_storeu_si128((__m128i *)(dst_line + 38 * dst_step), r[38]);
      _mm_storeu_si128((__m128i *)(dst_line + 39 * dst_step), r[39]);
      _mm_storeu_si128((__m128i *)(dst_line + 40 * dst_step), r[40]);
      _mm_storeu_si128((__m128i *)(dst_line + 41 * dst_step), r[41]);
      _mm_storeu_si128((__m128i *)(dst_line + 42 * dst_step), r[42]);
      _mm_storeu_si128((__m128i *)(dst_line + 43 * dst_step), r[43]);
      _mm_storeu_si128((__m128i *)(dst_line + 44 * dst_step), r[44]);
      _mm_storeu_si128((__m128i *)(dst_line + 45 * dst_step), r[45]);
      _mm_storeu_si128((__m128i *)(dst_line + 46 * dst_step), r[46]);
      _mm_storeu_si128((__m128i *)(dst_line + 47 * dst_step), r[47]);
      _mm_storeu_si128((__m128i *)(dst_line + 48 * dst_step), r[48]);
      _mm_storeu_si128((__m128i *)(dst_line + 49 * dst_step), r[49]);
      _mm_storeu_si128((__m128i *)(dst_line + 50 * dst_step), r[50]);
      _mm_storeu_si128((__m128i *)(dst_line + 51 * dst_step), r[51]);
      _mm_storeu_si128((__m128i *)(dst_line + 52 * dst_step), r[52]);
      _mm_storeu_si128((__m128i *)(dst_line + 53 * dst_step), r[53]);
      _mm_storeu_si128((__m128i *)(dst_line + 54 * dst_step), r[54]);
      _mm_storeu_si128((__m128i *)(dst_line + 55 * dst_step), r[55]);
      _mm_storeu_si128((__m128i *)(dst_line + 56 * dst_step), r[56]);
      _mm_storeu_si128((__m128i *)(dst_line + 57 * dst_step), r[57]);
      _mm_storeu_si128((__m128i *)(dst_line + 58 * dst_step), r[58]);
      _mm_storeu_si128((__m128i *)(dst_line + 59 * dst_step), r[59]);
      _mm_storeu_si128((__m128i *)(dst_line + 60 * dst_step), r[60]);
      _mm_storeu_si128((__m128i *)(dst_line + 61 * dst_step), r[61]);
      _mm_storeu_si128((__m128i *)(dst_line + 62 * dst_step), r[62]);
      _mm_storeu_si128((__m128i *)(dst_line + 63 * dst_step), r[63]);
      src_line  += 8*src_stride;
      dst_line  += 8*dst_stride;
      bndy_line += 8*ref_step;
    }
  }


void
rcn_init_mip_functions_sse(struct RCNFunctions *const rcn_funcs)
{
  rcn_funcs->mip.upsample_h[0][0]= &up_sample_h_4_2;
  rcn_funcs->mip.upsample_h[0][1]= &up_sample_h_4_4;
  rcn_funcs->mip.upsample_h[0][2]= &up_sample_h_4_8;

  rcn_funcs->mip.upsample_h[1][0]= &up_sample_h_8_2;
  rcn_funcs->mip.upsample_h[1][1]= &up_sample_h_8_4;
  rcn_funcs->mip.upsample_h[1][2]= &up_sample_h_8_8;

  rcn_funcs->mip.upsample_v[0][0]= &up_sample_v_4_2;
  rcn_funcs->mip.upsample_v[0][1]= &up_sample_v_4_4;
  rcn_funcs->mip.upsample_v[0][2]= &up_sample_v_4_8;

  rcn_funcs->mip.upsample_v[1][0]= &up_sample_v_8_2;
  rcn_funcs->mip.upsample_v[1][1]= &up_sample_v_8_4;
  rcn_funcs->mip.upsample_v[1][2]= &up_sample_v_8_8;

    // rcn_funcs->mip.upsample_v= &up_sample_v;
    rcn_funcs->mip.matmult= &mip_matmult;
}
