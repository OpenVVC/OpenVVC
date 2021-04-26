#include <emmintrin.h>
#include <smmintrin.h>
#include <stdint.h>
#include <stddef.h>

#include "ctudec.h"
#include "rcn_transform.h"
#include "rcn.h"
#include "ovutils.h"

#define CLIP_10 ((1 << 10) - 1)

static inline void
ovvc_transform_add_sse_8_4_10(uint16_t *dst, ptrdiff_t dst_stride,
                              const int16_t *src, ptrdiff_t src_stride)
{
   __m128i l0, l1, l2, l3, r0, r1, r2, r3;

   l0 = _mm_loadu_si128((__m128i *)&dst[0 * dst_stride]);
   l1 = _mm_loadu_si128((__m128i *)&dst[1 * dst_stride]);
   l2 = _mm_loadu_si128((__m128i *)&dst[2 * dst_stride]);
   l3 = _mm_loadu_si128((__m128i *)&dst[3 * dst_stride]);
   r0 = _mm_loadu_si128((__m128i *)&src[0 * src_stride]);
   r1 = _mm_loadu_si128((__m128i *)&src[1 * src_stride]);
   r2 = _mm_loadu_si128((__m128i *)&src[2 * src_stride]);
   r3 = _mm_loadu_si128((__m128i *)&src[3 * src_stride]);

   l0 = _mm_add_epi16(l0, r0);
   l1 = _mm_add_epi16(l1, r1);
   l2 = _mm_add_epi16(l2, r2);
   l3 = _mm_add_epi16(l3, r3);

   l0 = _mm_max_epi16(l0, _mm_setzero_si128());
   l1 = _mm_max_epi16(l1, _mm_setzero_si128());
   l2 = _mm_max_epi16(l2, _mm_setzero_si128());
   l3 = _mm_max_epi16(l3, _mm_setzero_si128());

   l0 = _mm_min_epi16(l0, _mm_set1_epi16(CLIP_10));
   l1 = _mm_min_epi16(l1, _mm_set1_epi16(CLIP_10));
   l2 = _mm_min_epi16(l2, _mm_set1_epi16(CLIP_10));
   l3 = _mm_min_epi16(l3, _mm_set1_epi16(CLIP_10));

   _mm_storeu_si128((__m128i*)&dst[0 * dst_stride], l0);
   _mm_storeu_si128((__m128i*)&dst[1 * dst_stride], l1);
   _mm_storeu_si128((__m128i*)&dst[2 * dst_stride], l2);
   _mm_storeu_si128((__m128i*)&dst[3 * dst_stride], l3);
}

static inline void
ovvc_transform_add_sse_16_2_10(uint16_t *dst, ptrdiff_t dst_stride,
                               const int16_t *src, ptrdiff_t src_stride)
{
   __m128i l0, l1, l2, l3, r0, r1, r2, r3;

   l0 = _mm_loadu_si128((__m128i *)&dst[0 * dst_stride]);
   l1 = _mm_loadu_si128((__m128i *)&dst[1 * dst_stride]);
   l2 = _mm_loadu_si128((__m128i *)&dst[8 + 0 * dst_stride]);
   l3 = _mm_loadu_si128((__m128i *)&dst[8 + 1 * dst_stride]);
   r0 = _mm_loadu_si128((__m128i *)&src[0 * src_stride]);
   r1 = _mm_loadu_si128((__m128i *)&src[1 * src_stride]);
   r2 = _mm_loadu_si128((__m128i *)&src[8 + 0 * src_stride]);
   r3 = _mm_loadu_si128((__m128i *)&src[8 + 1 * src_stride]);

   l0 = _mm_add_epi16(l0, r0);
   l1 = _mm_add_epi16(l1, r1);
   l2 = _mm_add_epi16(l2, r2);
   l3 = _mm_add_epi16(l3, r3);

   l0 = _mm_max_epi16(l0, _mm_setzero_si128());
   l1 = _mm_max_epi16(l1, _mm_setzero_si128());
   l2 = _mm_max_epi16(l2, _mm_setzero_si128());
   l3 = _mm_max_epi16(l3, _mm_setzero_si128());

   l0 = _mm_min_epi16(l0, _mm_set1_epi16(CLIP_10));
   l1 = _mm_min_epi16(l1, _mm_set1_epi16(CLIP_10));
   l2 = _mm_min_epi16(l2, _mm_set1_epi16(CLIP_10));
   l3 = _mm_min_epi16(l3, _mm_set1_epi16(CLIP_10));

   _mm_storeu_si128((__m128i*)&dst[0 * dst_stride], l0);
   _mm_storeu_si128((__m128i*)&dst[1 * dst_stride], l1);
   _mm_storeu_si128((__m128i*)&dst[8 + 0 * dst_stride], l2);
   _mm_storeu_si128((__m128i*)&dst[8 + 1 * dst_stride], l3);
}

static inline void
ovvc_transform_add_sse_32_1_10(uint16_t *dst, ptrdiff_t dst_stride,
                               const int16_t *src, ptrdiff_t src_stride)
{
   __m128i l0, l1, l2, l3, r0, r1, r2, r3;

   l0 = _mm_loadu_si128((__m128i *)&dst[0 * dst_stride]);
   l1 = _mm_loadu_si128((__m128i *)&dst[8  + 0 * dst_stride]);
   l2 = _mm_loadu_si128((__m128i *)&dst[16 + 0 * dst_stride]);
   l3 = _mm_loadu_si128((__m128i *)&dst[24 + 0 * dst_stride]);
   r0 = _mm_loadu_si128((__m128i *)&src[0 * src_stride]);
   r1 = _mm_loadu_si128((__m128i *)&src[8  + 0 * src_stride]);
   r2 = _mm_loadu_si128((__m128i *)&src[16 + 0 * src_stride]);
   r3 = _mm_loadu_si128((__m128i *)&src[24 + 0 * src_stride]);

   l0 = _mm_add_epi16(l0, r0);
   l1 = _mm_add_epi16(l1, r1);
   l2 = _mm_add_epi16(l2, r2);
   l3 = _mm_add_epi16(l3, r3);

   l0 = _mm_max_epi16(l0, _mm_setzero_si128());
   l1 = _mm_max_epi16(l1, _mm_setzero_si128());
   l2 = _mm_max_epi16(l2, _mm_setzero_si128());
   l3 = _mm_max_epi16(l3, _mm_setzero_si128());

   l0 = _mm_min_epi16(l0, _mm_set1_epi16(CLIP_10));
   l1 = _mm_min_epi16(l1, _mm_set1_epi16(CLIP_10));
   l2 = _mm_min_epi16(l2, _mm_set1_epi16(CLIP_10));
   l3 = _mm_min_epi16(l3, _mm_set1_epi16(CLIP_10));

   _mm_storeu_si128((__m128i*)&dst[0 * dst_stride], l0);
   _mm_storeu_si128((__m128i*)&dst[8  + 0 * dst_stride], l1);
   _mm_storeu_si128((__m128i*)&dst[16 + 0 * dst_stride], l2);
   _mm_storeu_si128((__m128i*)&dst[24 + 0 * dst_stride], l3);
}


static void
vvc_add_residual_8_4_10_sse(const int16_t *const src, uint16_t *const dst,
                            int log2_tb_w, int log2_tb_h, int scale)
{
    int i;
    int tb_w = 1 << log2_tb_w;
    int tb_h = 1 << log2_tb_h;
    const int16_t *_src = (const int16_t *)src;
    uint16_t *_dst = dst;
    if (log2_tb_h > 1) {
        for (i = 0; i < tb_h >> 2; ++i){
            ovvc_transform_add_sse_8_4_10(_dst, RCN_CTB_STRIDE,
                                          _src, tb_w);
            _dst += RCN_CTB_STRIDE << 2;
            _src += tb_w << 2;
        }
    } else {
        vvc_add_residual(src, dst, log2_tb_w, log2_tb_h, 0);
    }
}

static void
vvc_add_residual_16_2_10_sse(const int16_t *const src, uint16_t *const dst,
                            int log2_tb_w, int log2_tb_h, int scale)
{
    int i;
    int tb_w = 1 << log2_tb_w;
    int tb_h = 1 << log2_tb_h;
    const int16_t *_src = (const int16_t *)src;
    uint16_t *_dst = dst;
    if (log2_tb_h > 1) {
      for (i = 0; i < tb_h >> 1; ++i){
          ovvc_transform_add_sse_16_2_10(_dst, RCN_CTB_STRIDE,
                                         _src, tb_w);
          _dst += RCN_CTB_STRIDE << 1;
          _src += tb_w << 1;
      }
    } else {
        vvc_add_residual(src, dst, log2_tb_w, log2_tb_h, 0);
    }
}

void
vvc_add_residual_32_1_10_sse(const int16_t *const src, uint16_t *const dst,
                     int log2_tb_w, int log2_tb_h, int scale)
{
    int i;
    int tb_w = 1 << log2_tb_w;
    int tb_h = 1 << log2_tb_h;
    const int16_t *_src = (const int16_t *)src;
    uint16_t *_dst = dst;
    for (i = 0; i < tb_h; ++i){
        ovvc_transform_add_sse_32_1_10(_dst, RCN_CTB_STRIDE,
                                       _src, tb_w);
        _dst += RCN_CTB_STRIDE;
        _src += tb_w;
    }
}

void
rcn_init_ict_functions_sse(struct RCNFunctions *rcn_func, uint8_t type)
{
    switch (type)
    {
        case 3:
            // rcn_func->ict[0][0] = &vvc_scale_add_residual;
            // rcn_func->ict[1][0] = &vvc_scale_add_residual;
            // rcn_func->ict[2][0] = &vvc_scale_add_residual;
            // rcn_func->ict[3][0] = &vvc_scale_add_residual;
            // rcn_func->ict[4][0] = &vvc_scale_add_residual;
            // rcn_func->ict[5][0] = &vvc_scale_add_residual;
            // rcn_func->ict[6][0] = &vvc_scale_add_residual;
            //
            // rcn_func->ict[0][1] = &vvc_scale_sub_residual;
            // rcn_func->ict[1][1] = &vvc_scale_sub_residual;
            // rcn_func->ict[2][1] = &vvc_scale_sub_residual;
            // rcn_func->ict[3][1] = &vvc_scale_sub_residual;
            // rcn_func->ict[4][1] = &vvc_scale_sub_residual;
            // rcn_func->ict[5][1] = &vvc_scale_sub_residual;
            // rcn_func->ict[6][1] = &vvc_scale_sub_residual;
            //
            // rcn_func->ict[0][2] = &vvc_scale_sub_half_residual;
            // rcn_func->ict[1][2] = &vvc_scale_sub_half_residual;
            // rcn_func->ict[2][2] = &vvc_scale_sub_half_residual;
            // rcn_func->ict[3][2] = &vvc_scale_sub_half_residual;
            // rcn_func->ict[4][2] = &vvc_scale_sub_half_residual;
            // rcn_func->ict[5][2] = &vvc_scale_sub_half_residual;
            // rcn_func->ict[6][2] = &vvc_scale_sub_half_residual;
            break;
        case 2:
            // rcn_func->ict[0][0] = &vvc_add_residual;
            // rcn_func->ict[1][0] = &vvc_add_residual;
            // rcn_func->ict[2][0] = &vvc_add_residual;
            rcn_func->ict[3][0] = &vvc_add_residual_8_4_10_sse;
            rcn_func->ict[4][0] = &vvc_add_residual_16_2_10_sse;
            rcn_func->ict[5][0] = &vvc_add_residual_32_1_10_sse;
            // rcn_func->ict[6][0] = &vvc_add_residual;

            // rcn_func->ict[0][1] = &vvc_sub_residual;
            // rcn_func->ict[1][1] = &vvc_sub_residual;
            // rcn_func->ict[2][1] = &vvc_sub_residual;
            // rcn_func->ict[3][1] = &vvc_sub_residual;
            // rcn_func->ict[4][1] = &vvc_sub_residual;
            // rcn_func->ict[5][1] = &vvc_sub_residual;
            // rcn_func->ict[6][1] = &vvc_sub_residual;
            //
            // rcn_func->ict[0][2] = &vvc_sub_half_residual;
            // rcn_func->ict[1][2] = &vvc_sub_half_residual;
            // rcn_func->ict[2][2] = &vvc_sub_half_residual;
            // rcn_func->ict[3][2] = &vvc_sub_half_residual;
            // rcn_func->ict[4][2] = &vvc_sub_half_residual;
            // rcn_func->ict[5][2] = &vvc_sub_half_residual;
            // rcn_func->ict[6][2] = &vvc_sub_half_residual;
            break;
        case 1:
            // rcn_func->ict[0][0] = &vvc_scale_add_residual;
            // rcn_func->ict[1][0] = &vvc_scale_add_residual;
            // rcn_func->ict[2][0] = &vvc_scale_add_residual;
            // rcn_func->ict[3][0] = &vvc_scale_add_residual;
            // rcn_func->ict[4][0] = &vvc_scale_add_residual;
            // rcn_func->ict[5][0] = &vvc_scale_add_residual;
            // rcn_func->ict[6][0] = &vvc_scale_add_residual;
            //
            // rcn_func->ict[0][1] = &vvc_scale_add_residual;
            // rcn_func->ict[1][1] = &vvc_scale_add_residual;
            // rcn_func->ict[2][1] = &vvc_scale_add_residual;
            // rcn_func->ict[3][1] = &vvc_scale_add_residual;
            // rcn_func->ict[4][1] = &vvc_scale_add_residual;
            // rcn_func->ict[5][1] = &vvc_scale_add_residual;
            // rcn_func->ict[6][1] = &vvc_scale_add_residual;
            //
            // rcn_func->ict[0][2] = &vvc_scale_add_half_residual;
            // rcn_func->ict[1][2] = &vvc_scale_add_half_residual;
            // rcn_func->ict[2][2] = &vvc_scale_add_half_residual;
            // rcn_func->ict[3][2] = &vvc_scale_add_half_residual;
            // rcn_func->ict[4][2] = &vvc_scale_add_half_residual;
            // rcn_func->ict[5][2] = &vvc_scale_add_half_residual;
            // rcn_func->ict[6][2] = &vvc_scale_add_half_residual;
            break;
        default:
            // rcn_func->ict[0][0] = &vvc_add_residual;
            // rcn_func->ict[1][0] = &vvc_add_residual;
            // rcn_func->ict[2][0] = &vvc_add_residual;
            rcn_func->ict[3][0] = &vvc_add_residual_8_4_10_sse;
            rcn_func->ict[4][0] = &vvc_add_residual_16_2_10_sse;
            rcn_func->ict[5][0] = &vvc_add_residual_32_1_10_sse;
            // rcn_func->ict[6][0] = &vvc_add_residual;

            // rcn_func->ict[0][1] = &vvc_add_residual;
            // rcn_func->ict[1][1] = &vvc_add_residual;
            // rcn_func->ict[2][1] = &vvc_add_residual;
            rcn_func->ict[3][1] = &vvc_add_residual_8_4_10_sse;
            rcn_func->ict[4][1] = &vvc_add_residual_16_2_10_sse;
            rcn_func->ict[5][1] = &vvc_add_residual_32_1_10_sse;
            // rcn_func->ict[6][1] = &vvc_add_residual;

            // rcn_func->ict[0][2] = &vvc_add_half_residual;
            // rcn_func->ict[1][2] = &vvc_add_half_residual;
            // rcn_func->ict[2][2] = &vvc_add_half_residual;
            // rcn_func->ict[3][2] = &vvc_add_half_residual;
            // rcn_func->ict[4][2] = &vvc_add_half_residual;
            // rcn_func->ict[5][2] = &vvc_add_half_residual;
            // rcn_func->ict[6][2] = &vvc_add_half_residual;
            break;
    }
}
