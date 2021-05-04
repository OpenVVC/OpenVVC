#include <emmintrin.h>
#include <smmintrin.h>
#include <stdint.h>

#include "rcn_structures.h"
#include "ovutils.h"
#include "ctudec.h"

#define CLIP_10 ((1 << 10) - 1)

static void
compute_lm_subsample_4_sse(const uint16_t *lm_src, uint16_t *dst_cb, uint16_t *dst_cr,
                    ptrdiff_t lm_src_stride, ptrdiff_t dst_stride_c,
                    const struct CCLMParams *const lm_params,
                    int pb_w, int pb_h)
                    {
                      int i, j;
                      ptrdiff_t lm_src_stride2 = lm_src_stride << 1;
                      const __m128i scale_cb  = _mm_set1_epi16(lm_params->cb.a);
                      const int shift_cb  = lm_params->cb.shift;
                      const __m128i offset_cb = _mm_set1_epi32(lm_params->cb.b);
                      const __m128i scale_cr  = _mm_set1_epi16(lm_params->cr.a);
                      const int shift_cr  = lm_params->cr.shift;
                      const __m128i offset_cr = _mm_set1_epi32(lm_params->cr.b);

                      for (j = 0; j < pb_h; j++) {
                        __m128i x[4], a[2], t, r[2];
                        x[0] = _mm_loadu_si128((__m128i *) &lm_src[0]);
                        x[2] = _mm_loadu_si128((__m128i *) &lm_src[lm_src_stride]);
                        x[1] = _mm_or_si128(
                               _mm_slli_si128(x[0], 2),
                               _mm_setr_epi16(lm_src[0],0,0,0 ,0,0,0,0)
                             );
                        x[3] = _mm_or_si128(
                               _mm_slli_si128(x[2], 2),
                               _mm_setr_epi16(lm_src[lm_src_stride],0,0,0 ,0,0,0,0)
                             );

                        a[0] = _mm_add_epi16(x[0], x[2]);
                        a[1] = _mm_add_epi16(x[1], x[3]);

                        t = _mm_add_epi16(a[0], a[1]);
                        t = _mm_hadd_epi16(t, _mm_setzero_si128());
                        t = _mm_add_epi16(t, _mm_set1_epi16(4));
                        t = _mm_srai_epi16(t, 3);

                        x[0] = _mm_mullo_epi16(t, scale_cb);
                        x[1] = _mm_mulhi_epi16(t, scale_cb);
                        x[2] = _mm_mullo_epi16(t, scale_cr);
                        x[3] = _mm_mulhi_epi16(t, scale_cr);

                        r[0] = _mm_unpacklo_epi16(x[0],x[1]);
                        r[1] = _mm_unpacklo_epi16(x[2],x[3]);

                        r[0] = _mm_srai_epi32(r[0], shift_cb);
                        r[1] = _mm_srai_epi32(r[1], shift_cr);

                        r[0] = _mm_add_epi32(r[0], offset_cb);
                        r[1] = _mm_add_epi32(r[1], offset_cr);

                        r[0] = _mm_packs_epi32(r[0], _mm_setzero_si128());
                        r[1] = _mm_packs_epi32(r[1], _mm_setzero_si128());

                        r[0] = _mm_max_epi16(r[0], _mm_setzero_si128());
                        r[1] = _mm_max_epi16(r[1], _mm_setzero_si128());

                        r[0] = _mm_min_epi16(r[0], _mm_set1_epi16(CLIP_10));
                        r[1] = _mm_min_epi16(r[1], _mm_set1_epi16(CLIP_10));

                        _mm_storel_epi64((__m128i *)&dst_cb[0], r[0]);
                        _mm_storel_epi64((__m128i *)&dst_cr[0], r[1]);
                          for (i = 4; i < pb_w; i+=4) {
                              __m128i x[4], a[2], t, r[2];
                              x[0] = _mm_loadu_si128((__m128i *) &lm_src[2 * i - 1]);
                              x[1] = _mm_loadu_si128((__m128i *) &lm_src[2 * i]);
                              x[2] = _mm_loadu_si128((__m128i *) &lm_src[2 * i + lm_src_stride - 1]);
                              x[3] = _mm_loadu_si128((__m128i *) &lm_src[2 * i + lm_src_stride]);

                              a[0] = _mm_add_epi16(x[0], x[2]);
                              a[1] = _mm_add_epi16(x[1], x[3]);

                              t = _mm_add_epi16(a[0], a[1]);
                              t = _mm_hadd_epi16(t, _mm_setzero_si128());
                              t = _mm_add_epi16(t, _mm_set1_epi16(4));
                              t = _mm_srai_epi16(t, 3);

                              x[0] = _mm_mullo_epi16(t, scale_cb);
                              x[1] = _mm_mulhi_epi16(t, scale_cb);
                              x[2] = _mm_mullo_epi16(t, scale_cr);
                              x[3] = _mm_mulhi_epi16(t, scale_cr);

                              r[0] = _mm_unpacklo_epi16(x[0],x[1]);
                              r[1] = _mm_unpacklo_epi16(x[2],x[3]);

                              r[0] = _mm_srai_epi32(r[0], shift_cb);
                              r[1] = _mm_srai_epi32(r[1], shift_cr);

                              r[0] = _mm_add_epi32(r[0], offset_cb);
                              r[1] = _mm_add_epi32(r[1], offset_cr);

                              r[0] = _mm_packs_epi32(r[0], _mm_setzero_si128());
                              r[1] = _mm_packs_epi32(r[1], _mm_setzero_si128());

                              r[0] = _mm_max_epi16(r[0], _mm_setzero_si128());
                              r[1] = _mm_max_epi16(r[1], _mm_setzero_si128());

                              r[0] = _mm_min_epi16(r[0], _mm_set1_epi16(CLIP_10));
                              r[1] = _mm_min_epi16(r[1], _mm_set1_epi16(CLIP_10));

                              _mm_storel_epi64((__m128i *)&dst_cb[i], r[0]);
                              _mm_storel_epi64((__m128i *)&dst_cr[i], r[1]);
                          }
                          dst_cb += dst_stride_c;
                          dst_cr += dst_stride_c;
                          lm_src += lm_src_stride2;
                      }
                    }
static void
compute_lm_subsample_8_sse(const uint16_t *lm_src, uint16_t *dst_cb, uint16_t *dst_cr,
                    ptrdiff_t lm_src_stride, ptrdiff_t dst_stride_c,
                    const struct CCLMParams *const lm_params,
                    int pb_w, int pb_h)
                    {
                      int i, j;
                      ptrdiff_t lm_src_stride2 = lm_src_stride << 1;
                      const __m128i scale_cb  = _mm_set1_epi16(lm_params->cb.a);
                      const int shift_cb  = lm_params->cb.shift;
                      const __m128i offset_cb = _mm_set1_epi32(lm_params->cb.b);
                      const __m128i scale_cr  = _mm_set1_epi16(lm_params->cr.a);
                      const int shift_cr  = lm_params->cr.shift;
                      const __m128i offset_cr = _mm_set1_epi32(lm_params->cr.b);

                      for (j = 0; j < pb_h; j++) {
                          __m128i x[8], a[4], t, r[4];
                          x[0] = _mm_loadu_si128((__m128i *) &lm_src[0]);
                          x[2] = _mm_loadu_si128((__m128i *) &lm_src[lm_src_stride]);
                          x[1] = _mm_or_si128(
                                 _mm_slli_si128(x[0], 2),
                                 _mm_setr_epi16(lm_src[0],0,0,0 ,0,0,0,0)
                               );
                          x[3] = _mm_or_si128(
                                 _mm_slli_si128(x[2], 2),
                                 _mm_setr_epi16(lm_src[lm_src_stride],0,0,0 ,0,0,0,0)
                               );

                          x[4] = _mm_loadu_si128((__m128i *) &lm_src[7]);
                          x[5] = _mm_loadu_si128((__m128i *) &lm_src[8]);
                          x[6] = _mm_loadu_si128((__m128i *) &lm_src[7 + lm_src_stride]);
                          x[7] = _mm_loadu_si128((__m128i *) &lm_src[8 + lm_src_stride]);

                          a[0] = _mm_add_epi16(x[0], x[2]);
                          a[1] = _mm_add_epi16(x[1], x[3]);
                          a[2] = _mm_add_epi16(x[4], x[6]);
                          a[3] = _mm_add_epi16(x[5], x[7]);

                          t = _mm_hadd_epi16(_mm_add_epi16(a[0], a[1]), _mm_add_epi16(a[2], a[3]));
                          t = _mm_add_epi16(t, _mm_set1_epi16(4));
                          t = _mm_srai_epi16(t, 3);

                          x[0] = _mm_mullo_epi16(t, scale_cb);
                          x[1] = _mm_mulhi_epi16(t, scale_cb);
                          x[2] = _mm_mullo_epi16(t, scale_cr);
                          x[3] = _mm_mulhi_epi16(t, scale_cr);

                          r[0] = _mm_unpacklo_epi16(x[0],x[1]);
                          r[1] = _mm_unpacklo_epi16(x[2],x[3]);
                          r[2] = _mm_unpackhi_epi16(x[0],x[1]);
                          r[3] = _mm_unpackhi_epi16(x[2],x[3]);

                          r[0] = _mm_srai_epi32(r[0], shift_cb);
                          r[1] = _mm_srai_epi32(r[1], shift_cr);
                          r[2] = _mm_srai_epi32(r[2], shift_cb);
                          r[3] = _mm_srai_epi32(r[3], shift_cr);

                          r[0] = _mm_add_epi32(r[0], offset_cb);
                          r[1] = _mm_add_epi32(r[1], offset_cr);
                          r[2] = _mm_add_epi32(r[2], offset_cb);
                          r[3] = _mm_add_epi32(r[3], offset_cr);

                          r[0] = _mm_packs_epi32(r[0], r[2]);
                          r[1] = _mm_packs_epi32(r[1], r[3]);

                          r[0] = _mm_max_epi16(r[0], _mm_setzero_si128());
                          r[1] = _mm_max_epi16(r[1], _mm_setzero_si128());

                          r[0] = _mm_min_epi16(r[0], _mm_set1_epi16(CLIP_10));
                          r[1] = _mm_min_epi16(r[1], _mm_set1_epi16(CLIP_10));

                          _mm_storeu_si128((__m128i *)&dst_cb[0], r[0]);
                          _mm_storeu_si128((__m128i *)&dst_cr[0], r[1]);
                          for (i = 8; i < pb_w; i+=8) {
                              __m128i x[8], a[4], t, r[4];
                              x[0] = _mm_loadu_si128((__m128i *) &lm_src[2 * i - 1]);
                              x[1] = _mm_loadu_si128((__m128i *) &lm_src[2 * i]);
                              x[2] = _mm_loadu_si128((__m128i *) &lm_src[2 * i + lm_src_stride - 1]);
                              x[3] = _mm_loadu_si128((__m128i *) &lm_src[2 * i + lm_src_stride]);

                              x[4] = _mm_loadu_si128((__m128i *) &lm_src[2 * (i+4) - 1]);
                              x[5] = _mm_loadu_si128((__m128i *) &lm_src[2 * (i+4)]);
                              x[6] = _mm_loadu_si128((__m128i *) &lm_src[2 * (i+4) + lm_src_stride - 1]);
                              x[7] = _mm_loadu_si128((__m128i *) &lm_src[2 * (i+4) + lm_src_stride]);

                              a[0] = _mm_add_epi16(x[0], x[2]);
                              a[1] = _mm_add_epi16(x[1], x[3]);
                              a[2] = _mm_add_epi16(x[4], x[6]);
                              a[3] = _mm_add_epi16(x[5], x[7]);

                              t = _mm_hadd_epi16(_mm_add_epi16(a[0], a[1]), _mm_add_epi16(a[2], a[3]));
                              t = _mm_add_epi16(t, _mm_set1_epi16(4));
                              t = _mm_srai_epi16(t, 3);

                              x[0] = _mm_mullo_epi16(t, scale_cb);
                              x[1] = _mm_mulhi_epi16(t, scale_cb);
                              x[2] = _mm_mullo_epi16(t, scale_cr);
                              x[3] = _mm_mulhi_epi16(t, scale_cr);

                              r[0] = _mm_unpacklo_epi16(x[0],x[1]);
                              r[1] = _mm_unpacklo_epi16(x[2],x[3]);
                              r[2] = _mm_unpackhi_epi16(x[0],x[1]);
                              r[3] = _mm_unpackhi_epi16(x[2],x[3]);

                              r[0] = _mm_srai_epi32(r[0], shift_cb);
                              r[1] = _mm_srai_epi32(r[1], shift_cr);
                              r[2] = _mm_srai_epi32(r[2], shift_cb);
                              r[3] = _mm_srai_epi32(r[3], shift_cr);

                              r[0] = _mm_add_epi32(r[0], offset_cb);
                              r[1] = _mm_add_epi32(r[1], offset_cr);
                              r[2] = _mm_add_epi32(r[2], offset_cb);
                              r[3] = _mm_add_epi32(r[3], offset_cr);

                              r[0] = _mm_packs_epi32(r[0], r[2]);
                              r[1] = _mm_packs_epi32(r[1], r[3]);

                              r[0] = _mm_max_epi16(r[0], _mm_setzero_si128());
                              r[1] = _mm_max_epi16(r[1], _mm_setzero_si128());

                              r[0] = _mm_min_epi16(r[0], _mm_set1_epi16(CLIP_10));
                              r[1] = _mm_min_epi16(r[1], _mm_set1_epi16(CLIP_10));

                              _mm_storeu_si128((__m128i *)&dst_cb[i], r[0]);
                              _mm_storeu_si128((__m128i *)&dst_cr[i], r[1]);
                          }
                          dst_cb += dst_stride_c;
                          dst_cr += dst_stride_c;
                          lm_src += lm_src_stride2;
                      }
                    }

static void
compute_lm_subsample_4_lft_avail_sse(const uint16_t *lm_src, uint16_t *dst_cb, uint16_t *dst_cr,
                    ptrdiff_t lm_src_stride, ptrdiff_t dst_stride_c,
                    const struct CCLMParams *const lm_params,
                    int pb_w, int pb_h)
                    {
                      int i, j;
                      ptrdiff_t lm_src_stride2 = lm_src_stride << 1;
                      const __m128i scale_cb  = _mm_set1_epi16(lm_params->cb.a);
                      const int shift_cb  = lm_params->cb.shift;
                      const __m128i offset_cb = _mm_set1_epi32(lm_params->cb.b);
                      const __m128i scale_cr  = _mm_set1_epi16(lm_params->cr.a);
                      const int shift_cr  = lm_params->cr.shift;
                      const __m128i offset_cr = _mm_set1_epi32(lm_params->cr.b);

                      for (j = 0; j < pb_h; j++) {
                          for (i = 0; i < pb_w; i+=4) {
                              __m128i x[4], a[2], t, r[2];
                              x[0] = _mm_loadu_si128((__m128i *) &lm_src[2 * i - 1]);
                              x[1] = _mm_loadu_si128((__m128i *) &lm_src[2 * i]);
                              x[2] = _mm_loadu_si128((__m128i *) &lm_src[2 * i + lm_src_stride - 1]);
                              x[3] = _mm_loadu_si128((__m128i *) &lm_src[2 * i + lm_src_stride]);

                              a[0] = _mm_add_epi16(x[0], x[2]);
                              a[1] = _mm_add_epi16(x[1], x[3]);

                              t = _mm_add_epi16(a[0], a[1]);
                              t = _mm_hadd_epi16(t, _mm_setzero_si128());
                              t = _mm_add_epi16(t, _mm_set1_epi16(4));
                              t = _mm_srai_epi16(t, 3);

                              x[0] = _mm_mullo_epi16(t, scale_cb);
                              x[1] = _mm_mulhi_epi16(t, scale_cb);
                              x[2] = _mm_mullo_epi16(t, scale_cr);
                              x[3] = _mm_mulhi_epi16(t, scale_cr);

                              r[0] = _mm_unpacklo_epi16(x[0],x[1]);
                              r[1] = _mm_unpacklo_epi16(x[2],x[3]);

                              r[0] = _mm_srai_epi32(r[0], shift_cb);
                              r[1] = _mm_srai_epi32(r[1], shift_cr);

                              r[0] = _mm_add_epi32(r[0], offset_cb);
                              r[1] = _mm_add_epi32(r[1], offset_cr);

                              r[0] = _mm_packs_epi32(r[0], _mm_setzero_si128());
                              r[1] = _mm_packs_epi32(r[1], _mm_setzero_si128());

                              r[0] = _mm_max_epi16(r[0], _mm_setzero_si128());
                              r[1] = _mm_max_epi16(r[1], _mm_setzero_si128());

                              r[0] = _mm_min_epi16(r[0], _mm_set1_epi16(CLIP_10));
                              r[1] = _mm_min_epi16(r[1], _mm_set1_epi16(CLIP_10));

                              _mm_storel_epi64((__m128i *)&dst_cb[i], r[0]);
                              _mm_storel_epi64((__m128i *)&dst_cr[i], r[1]);
                          }
                          dst_cb += dst_stride_c;
                          dst_cr += dst_stride_c;
                          lm_src += lm_src_stride2;
                      }
                    }

static void
compute_lm_subsample_8_lft_avail_sse(const uint16_t *lm_src, uint16_t *dst_cb, uint16_t *dst_cr,
                    ptrdiff_t lm_src_stride, ptrdiff_t dst_stride_c,
                    const struct CCLMParams *const lm_params,
                    int pb_w, int pb_h)
                    {
                      int i, j;
                      ptrdiff_t lm_src_stride2 = lm_src_stride << 1;
                      const __m128i scale_cb  = _mm_set1_epi16(lm_params->cb.a);
                      const int shift_cb  = lm_params->cb.shift;
                      const __m128i offset_cb = _mm_set1_epi32(lm_params->cb.b);
                      const __m128i scale_cr  = _mm_set1_epi16(lm_params->cr.a);
                      const int shift_cr  = lm_params->cr.shift;
                      const __m128i offset_cr = _mm_set1_epi32(lm_params->cr.b);

                      for (j = 0; j < pb_h; j++) {
                          for (i = 0; i < pb_w; i+=8) {
                              __m128i x[8], a[4], t, r[4];
                              x[0] = _mm_loadu_si128((__m128i *) &lm_src[2 * i - 1]);
                              x[1] = _mm_loadu_si128((__m128i *) &lm_src[2 * i]);
                              x[2] = _mm_loadu_si128((__m128i *) &lm_src[2 * i + lm_src_stride - 1]);
                              x[3] = _mm_loadu_si128((__m128i *) &lm_src[2 * i + lm_src_stride]);

                              x[4] = _mm_loadu_si128((__m128i *) &lm_src[2 * (i+4) - 1]);
                              x[5] = _mm_loadu_si128((__m128i *) &lm_src[2 * (i+4)]);
                              x[6] = _mm_loadu_si128((__m128i *) &lm_src[2 * (i+4) + lm_src_stride - 1]);
                              x[7] = _mm_loadu_si128((__m128i *) &lm_src[2 * (i+4) + lm_src_stride]);

                              a[0] = _mm_add_epi16(x[0], x[2]);
                              a[1] = _mm_add_epi16(x[1], x[3]);
                              a[2] = _mm_add_epi16(x[4], x[6]);
                              a[3] = _mm_add_epi16(x[5], x[7]);

                              t = _mm_hadd_epi16(_mm_add_epi16(a[0], a[1]), _mm_add_epi16(a[2], a[3]));
                              t = _mm_add_epi16(t, _mm_set1_epi16(4));
                              t = _mm_srai_epi16(t, 3);

                              x[0] = _mm_mullo_epi16(t, scale_cb);
                              x[1] = _mm_mulhi_epi16(t, scale_cb);
                              x[2] = _mm_mullo_epi16(t, scale_cr);
                              x[3] = _mm_mulhi_epi16(t, scale_cr);

                              r[0] = _mm_unpacklo_epi16(x[0],x[1]);
                              r[1] = _mm_unpacklo_epi16(x[2],x[3]);
                              r[2] = _mm_unpackhi_epi16(x[0],x[1]);
                              r[3] = _mm_unpackhi_epi16(x[2],x[3]);

                              r[0] = _mm_srai_epi32(r[0], shift_cb);
                              r[1] = _mm_srai_epi32(r[1], shift_cr);
                              r[2] = _mm_srai_epi32(r[2], shift_cb);
                              r[3] = _mm_srai_epi32(r[3], shift_cr);

                              r[0] = _mm_add_epi32(r[0], offset_cb);
                              r[1] = _mm_add_epi32(r[1], offset_cr);
                              r[2] = _mm_add_epi32(r[2], offset_cb);
                              r[3] = _mm_add_epi32(r[3], offset_cr);

                              r[0] = _mm_packs_epi32(r[0], r[2]);
                              r[1] = _mm_packs_epi32(r[1], r[3]);

                              r[0] = _mm_max_epi16(r[0], _mm_setzero_si128());
                              r[1] = _mm_max_epi16(r[1], _mm_setzero_si128());

                              r[0] = _mm_min_epi16(r[0], _mm_set1_epi16(CLIP_10));
                              r[1] = _mm_min_epi16(r[1], _mm_set1_epi16(CLIP_10));

                              _mm_storeu_si128((__m128i *)&dst_cb[i], r[0]);
                              _mm_storeu_si128((__m128i *)&dst_cr[i], r[1]);
                          }
                          dst_cb += dst_stride_c;
                          dst_cr += dst_stride_c;
                          lm_src += lm_src_stride2;
                      }
                    }

static void
compute_lm_subsample_sse(const uint16_t *lm_src, uint16_t *dst_cb, uint16_t *dst_cr,
                     ptrdiff_t lm_src_stride, ptrdiff_t dst_stride_c,
                     const struct CCLMParams *const lm_params,
                     int pb_w, int pb_h, uint8_t lft_avail)
{
  if(lft_avail) {
    if (pb_w > 4) {
      compute_lm_subsample_8_lft_avail_sse(lm_src, dst_cb, dst_cr, lm_src_stride, dst_stride_c, lm_params, pb_w, pb_h);
    }
    else {
      compute_lm_subsample_4_lft_avail_sse(lm_src, dst_cb, dst_cr, lm_src_stride, dst_stride_c, lm_params, pb_w, pb_h);
    }
  }
  else {
    if (pb_w > 4) {
      compute_lm_subsample_8_sse(lm_src, dst_cb, dst_cr, lm_src_stride, dst_stride_c, lm_params, pb_w, pb_h);
    }
    else {
      compute_lm_subsample_4_sse(lm_src, dst_cb, dst_cr, lm_src_stride, dst_stride_c, lm_params, pb_w, pb_h);
    }
  }
}


void
rcn_init_cclm_functions_sse(struct RCNFunctions *rcn_func)
{
   struct CCLMFunctions *const cclm = &rcn_func->cclm;
   cclm->compute_subsample = &compute_lm_subsample_sse;
}
