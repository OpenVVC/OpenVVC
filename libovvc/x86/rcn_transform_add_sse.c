#include <emmintrin.h>
#include <smmintrin.h>
#include <tmmintrin.h>
#include <stdint.h>
#include <stddef.h>

#include "ctudec.h"
#include "rcn_transform.h"
#include "rcn.h"
#include "ovutils.h"

#define CLIP_10 ((1 << 10) - 1)
#define SIGN_16 (int16_t)(1 << 15)

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

static inline void
ovvc_transform_add_half_sse_8_4_10(uint16_t *dst, ptrdiff_t dst_stride,
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

   r0 = _mm_srai_epi16(r0, 1);
   r1 = _mm_srai_epi16(r1, 1);
   r2 = _mm_srai_epi16(r2, 1);
   r3 = _mm_srai_epi16(r3, 1);

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
ovvc_transform_add_half_sse_16_2_10(uint16_t *dst, ptrdiff_t dst_stride,
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

   r0 = _mm_srai_epi16(r0, 1);
   r1 = _mm_srai_epi16(r1, 1);
   r2 = _mm_srai_epi16(r2, 1);
   r3 = _mm_srai_epi16(r3, 1);

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
ovvc_transform_add_half_sse_32_1_10(uint16_t *dst, ptrdiff_t dst_stride,
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

   r0 = _mm_srai_epi16(r0, 1);
   r1 = _mm_srai_epi16(r1, 1);
   r2 = _mm_srai_epi16(r2, 1);
   r3 = _mm_srai_epi16(r3, 1);

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

static inline void
ovvc_transform_sub_sse_8_4_10(uint16_t *dst, ptrdiff_t dst_stride,
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

   l0 = _mm_sub_epi16(l0, r0);
   l1 = _mm_sub_epi16(l1, r1);
   l2 = _mm_sub_epi16(l2, r2);
   l3 = _mm_sub_epi16(l3, r3);

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
ovvc_transform_sub_sse_16_2_10(uint16_t *dst, ptrdiff_t dst_stride,
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

   l0 = _mm_sub_epi16(l0, r0);
   l1 = _mm_sub_epi16(l1, r1);
   l2 = _mm_sub_epi16(l2, r2);
   l3 = _mm_sub_epi16(l3, r3);

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
ovvc_transform_sub_sse_32_1_10(uint16_t *dst, ptrdiff_t dst_stride,
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

   l0 = _mm_sub_epi16(l0, r0);
   l1 = _mm_sub_epi16(l1, r1);
   l2 = _mm_sub_epi16(l2, r2);
   l3 = _mm_sub_epi16(l3, r3);

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

static inline void
ovvc_transform_sub_half_sse_8_4_10(uint16_t *dst, ptrdiff_t dst_stride,
                              const int16_t *src, ptrdiff_t src_stride)
{
   __m128i l0, l1, l2, l3, r0, r1, r2, r3, s0, s1, s2, s3;

   l0 = _mm_loadu_si128((__m128i *)&dst[0 * dst_stride]);
   l1 = _mm_loadu_si128((__m128i *)&dst[1 * dst_stride]);
   l2 = _mm_loadu_si128((__m128i *)&dst[2 * dst_stride]);
   l3 = _mm_loadu_si128((__m128i *)&dst[3 * dst_stride]);
   r0 = _mm_loadu_si128((__m128i *)&src[0 * src_stride]);
   r1 = _mm_loadu_si128((__m128i *)&src[1 * src_stride]);
   r2 = _mm_loadu_si128((__m128i *)&src[2 * src_stride]);
   r3 = _mm_loadu_si128((__m128i *)&src[3 * src_stride]);

   //         sign  = value & (1 << 15);
   s0 = _mm_srai_epi16(_mm_and_si128(r0, _mm_set1_epi16(SIGN_16)), 15);
   s1 = _mm_srai_epi16(_mm_and_si128(r1, _mm_set1_epi16(SIGN_16)), 15);
   s2 = _mm_srai_epi16(_mm_and_si128(r2, _mm_set1_epi16(SIGN_16)), 15);
   s3 = _mm_srai_epi16(_mm_and_si128(r3, _mm_set1_epi16(SIGN_16)), 15);

   s0 = _mm_andnot_si128(_mm_abs_epi16(s0), _mm_set1_epi16(1));
   s1 = _mm_andnot_si128(_mm_abs_epi16(s1), _mm_set1_epi16(1));
   s2 = _mm_andnot_si128(_mm_abs_epi16(s2), _mm_set1_epi16(1));
   s3 = _mm_andnot_si128(_mm_abs_epi16(s3), _mm_set1_epi16(1));

   //         value = (-value);
   r0 = _mm_xor_si128(r0, _mm_set1_epi16((int16_t)0xFFFF));
   r1 = _mm_xor_si128(r1, _mm_set1_epi16((int16_t)0xFFFF));
   r2 = _mm_xor_si128(r2, _mm_set1_epi16((int16_t)0xFFFF));
   r3 = _mm_xor_si128(r3, _mm_set1_epi16((int16_t)0xFFFF));

   r0 = _mm_add_epi16(r0, _mm_set1_epi16(1));
   r1 = _mm_add_epi16(r1, _mm_set1_epi16(1));
   r2 = _mm_add_epi16(r2, _mm_set1_epi16(1));
   r3 = _mm_add_epi16(r3, _mm_set1_epi16(1));

   //         value = (abs(value);
   r0 = _mm_srai_epi16(r0, 1);
   r1 = _mm_srai_epi16(r1, 1);
   r2 = _mm_srai_epi16(r2, 1);
   r3 = _mm_srai_epi16(r3, 1);

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
ovvc_transform_sub_half_sse_16_2_10(uint16_t *dst, ptrdiff_t dst_stride,
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

   //         value = (-value);
   r0 = _mm_xor_si128(r0, _mm_set1_epi16((int16_t)0xFFFF));
   r1 = _mm_xor_si128(r1, _mm_set1_epi16((int16_t)0xFFFF));
   r2 = _mm_xor_si128(r2, _mm_set1_epi16((int16_t)0xFFFF));
   r3 = _mm_xor_si128(r3, _mm_set1_epi16((int16_t)0xFFFF));

   r0 = _mm_add_epi16(r0, _mm_set1_epi16(1));
   r1 = _mm_add_epi16(r1, _mm_set1_epi16(1));
   r2 = _mm_add_epi16(r2, _mm_set1_epi16(1));
   r3 = _mm_add_epi16(r3, _mm_set1_epi16(1));

   //         value = (abs(value);
   r0 = _mm_srai_epi16(r0, 1);
   r1 = _mm_srai_epi16(r1, 1);
   r2 = _mm_srai_epi16(r2, 1);
   r3 = _mm_srai_epi16(r3, 1);

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
ovvc_transform_sub_half_sse_32_1_10(uint16_t *dst, ptrdiff_t dst_stride,
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

   //         value = (-value);
   r0 = _mm_xor_si128(r0, _mm_set1_epi16((int16_t)0xFFFF));
   r1 = _mm_xor_si128(r1, _mm_set1_epi16((int16_t)0xFFFF));
   r2 = _mm_xor_si128(r2, _mm_set1_epi16((int16_t)0xFFFF));
   r3 = _mm_xor_si128(r3, _mm_set1_epi16((int16_t)0xFFFF));

   r0 = _mm_add_epi16(r0, _mm_set1_epi16(1));
   r1 = _mm_add_epi16(r1, _mm_set1_epi16(1));
   r2 = _mm_add_epi16(r2, _mm_set1_epi16(1));
   r3 = _mm_add_epi16(r3, _mm_set1_epi16(1));

   //         value = (abs(value);
   r0 = _mm_srai_epi16(r0, 1);
   r1 = _mm_srai_epi16(r1, 1);
   r2 = _mm_srai_epi16(r2, 1);
   r3 = _mm_srai_epi16(r3, 1);

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

static inline void
ovvc_transform_scale_add_sse_8_4_10(uint16_t *dst, ptrdiff_t dst_stride,
                               const int16_t *src, ptrdiff_t src_stride, __m128i scale)
{
   __m128i l[4], r[4], s[4], m[8];
   __m128i a = _mm_unpackhi_epi16(scale, _mm_set1_epi16(1));
   __m128i b = _mm_set1_epi16(1 << (11 - 1));

   l[0] = _mm_loadu_si128((__m128i *)&dst[0 * dst_stride]);
   l[1] = _mm_loadu_si128((__m128i *)&dst[1 * dst_stride]);
   l[2] = _mm_loadu_si128((__m128i *)&dst[2 * dst_stride]);
   l[3] = _mm_loadu_si128((__m128i *)&dst[3 * dst_stride]);
   r[0] = _mm_loadu_si128((__m128i *)&src[0 * src_stride]);
   r[1] = _mm_loadu_si128((__m128i *)&src[1 * src_stride]);
   r[2] = _mm_loadu_si128((__m128i *)&src[2 * src_stride]);
   r[3] = _mm_loadu_si128((__m128i *)&src[3 * src_stride]);

   //         sign  = value & (1 << 15);
   s[0] = _mm_srai_epi16(_mm_and_si128(r[0], _mm_set1_epi16(SIGN_16)), 15);
   s[1] = _mm_srai_epi16(_mm_and_si128(r[1], _mm_set1_epi16(SIGN_16)), 15);
   s[2] = _mm_srai_epi16(_mm_and_si128(r[2], _mm_set1_epi16(SIGN_16)), 15);
   s[3] = _mm_srai_epi16(_mm_and_si128(r[3], _mm_set1_epi16(SIGN_16)), 15);

   s[0] = _mm_abs_epi16(s[0]);
   s[1] = _mm_abs_epi16(s[1]);
   s[2] = _mm_abs_epi16(s[2]);
   s[3] = _mm_abs_epi16(s[3]);

   //         value = (abs(value);
   r[0] = _mm_abs_epi16(r[0]);
   r[1] = _mm_abs_epi16(r[1]);
   r[2] = _mm_abs_epi16(r[2]);
   r[3] = _mm_abs_epi16(r[3]);

   //         value = (value * scale + (1 << (11 - 1)));
   m[0] = _mm_unpacklo_epi16(r[0], b);
   m[1] = _mm_unpackhi_epi16(r[0], b);
   m[2] = _mm_unpacklo_epi16(r[1], b);
   m[3] = _mm_unpackhi_epi16(r[1], b);
   m[4] = _mm_unpacklo_epi16(r[2], b);
   m[5] = _mm_unpackhi_epi16(r[2], b);
   m[6] = _mm_unpacklo_epi16(r[3], b);
   m[7] = _mm_unpackhi_epi16(r[3], b);

   m[0] = _mm_madd_epi16(m[0], a);
   m[1] = _mm_madd_epi16(m[1], a);
   m[2] = _mm_madd_epi16(m[2], a);
   m[3] = _mm_madd_epi16(m[3], a);
   m[4] = _mm_madd_epi16(m[4], a);
   m[5] = _mm_madd_epi16(m[5], a);
   m[6] = _mm_madd_epi16(m[6], a);
   m[7] = _mm_madd_epi16(m[7], a);

   //         value = value >> 11;
   m[0] = _mm_srai_epi32(m[0], 11);
   m[1] = _mm_srai_epi32(m[1], 11);
   m[2] = _mm_srai_epi32(m[2], 11);
   m[3] = _mm_srai_epi32(m[3], 11);
   m[4] = _mm_srai_epi32(m[4], 11);
   m[5] = _mm_srai_epi32(m[5], 11);
   m[6] = _mm_srai_epi32(m[6], 11);
   m[7] = _mm_srai_epi32(m[7], 11);


   //         value = (ov_clip(value , 0,(1 << 16)-1));
   r[0] = _mm_packs_epi32(m[0], m[1]);
   r[1] = _mm_packs_epi32(m[2], m[3]);
   r[2] = _mm_packs_epi32(m[4], m[5]);
   r[3] = _mm_packs_epi32(m[6], m[7]);


   //         value = (sign ? -value : value);
   r[0] = _mm_xor_si128(r[0],_mm_cmpgt_epi16(s[0], _mm_setzero_si128()));
   r[1] = _mm_xor_si128(r[1],_mm_cmpgt_epi16(s[1], _mm_setzero_si128()));
   r[2] = _mm_xor_si128(r[2],_mm_cmpgt_epi16(s[2], _mm_setzero_si128()));
   r[3] = _mm_xor_si128(r[3],_mm_cmpgt_epi16(s[3], _mm_setzero_si128()));

   r[0] = _mm_add_epi16(r[0], s[0]);
   r[1] = _mm_add_epi16(r[1], s[1]);
   r[2] = _mm_add_epi16(r[2], s[2]);
   r[3] = _mm_add_epi16(r[3], s[3]);

   //         _dst[j] = ov_clip((int32_t)_dst[j] + value, 0, 1023);
   l[0] = _mm_add_epi16(l[0], r[0]);
   l[1] = _mm_add_epi16(l[1], r[1]);
   l[2] = _mm_add_epi16(l[2], r[2]);
   l[3] = _mm_add_epi16(l[3], r[3]);

   l[0] = _mm_max_epi16(l[0], _mm_setzero_si128());
   l[1] = _mm_max_epi16(l[1], _mm_setzero_si128());
   l[2] = _mm_max_epi16(l[2], _mm_setzero_si128());
   l[3] = _mm_max_epi16(l[3], _mm_setzero_si128());

   l[0] = _mm_min_epi16(l[0], _mm_set1_epi16(CLIP_10));
   l[1] = _mm_min_epi16(l[1], _mm_set1_epi16(CLIP_10));
   l[2] = _mm_min_epi16(l[2], _mm_set1_epi16(CLIP_10));
   l[3] = _mm_min_epi16(l[3], _mm_set1_epi16(CLIP_10));

   _mm_storeu_si128((__m128i*)&dst[0 * dst_stride], l[0]);
   _mm_storeu_si128((__m128i*)&dst[1 * dst_stride], l[1]);
   _mm_storeu_si128((__m128i*)&dst[2 * dst_stride], l[2]);
   _mm_storeu_si128((__m128i*)&dst[3 * dst_stride], l[3]);
}

static inline void
ovvc_transform_scale_add_sse_16_2_10(uint16_t *dst, ptrdiff_t dst_stride,
                               const int16_t *src, ptrdiff_t src_stride, __m128i scale)
{
   __m128i l[4], r[4], s[4], m[8];
   __m128i a = _mm_unpackhi_epi16(scale, _mm_set1_epi16(1));
   __m128i b = _mm_set1_epi16(1 << (11 - 1));

   l[0] = _mm_loadu_si128((__m128i *)&dst[0 * dst_stride]);
   l[1] = _mm_loadu_si128((__m128i *)&dst[1 * dst_stride]);
   l[2] = _mm_loadu_si128((__m128i *)&dst[8 + 0 * dst_stride]);
   l[3] = _mm_loadu_si128((__m128i *)&dst[8 + 1 * dst_stride]);
   r[0] = _mm_loadu_si128((__m128i *)&src[0 * src_stride]);
   r[1] = _mm_loadu_si128((__m128i *)&src[1 * src_stride]);
   r[2] = _mm_loadu_si128((__m128i *)&src[8 + 0 * src_stride]);
   r[3] = _mm_loadu_si128((__m128i *)&src[8 + 1 * src_stride]);

   //         sign  = value & (1 << 15);
   s[0] = _mm_srai_epi16(_mm_and_si128(r[0], _mm_set1_epi16(SIGN_16)), 15);
   s[1] = _mm_srai_epi16(_mm_and_si128(r[1], _mm_set1_epi16(SIGN_16)), 15);
   s[2] = _mm_srai_epi16(_mm_and_si128(r[2], _mm_set1_epi16(SIGN_16)), 15);
   s[3] = _mm_srai_epi16(_mm_and_si128(r[3], _mm_set1_epi16(SIGN_16)), 15);

   s[0] = _mm_abs_epi16(s[0]);
   s[1] = _mm_abs_epi16(s[1]);
   s[2] = _mm_abs_epi16(s[2]);
   s[3] = _mm_abs_epi16(s[3]);

   //         value = (abs(value);
   r[0] = _mm_abs_epi16(r[0]);
   r[1] = _mm_abs_epi16(r[1]);
   r[2] = _mm_abs_epi16(r[2]);
   r[3] = _mm_abs_epi16(r[3]);

   //         value = (value * scale + (1 << (11 - 1)));
   m[0] = _mm_unpacklo_epi16(r[0], b);
   m[1] = _mm_unpackhi_epi16(r[0], b);
   m[2] = _mm_unpacklo_epi16(r[1], b);
   m[3] = _mm_unpackhi_epi16(r[1], b);
   m[4] = _mm_unpacklo_epi16(r[2], b);
   m[5] = _mm_unpackhi_epi16(r[2], b);
   m[6] = _mm_unpacklo_epi16(r[3], b);
   m[7] = _mm_unpackhi_epi16(r[3], b);

   m[0] = _mm_madd_epi16(m[0], a);
   m[1] = _mm_madd_epi16(m[1], a);
   m[2] = _mm_madd_epi16(m[2], a);
   m[3] = _mm_madd_epi16(m[3], a);
   m[4] = _mm_madd_epi16(m[4], a);
   m[5] = _mm_madd_epi16(m[5], a);
   m[6] = _mm_madd_epi16(m[6], a);
   m[7] = _mm_madd_epi16(m[7], a);

   //         value = value >> 11;
   m[0] = _mm_srai_epi32(m[0], 11);
   m[1] = _mm_srai_epi32(m[1], 11);
   m[2] = _mm_srai_epi32(m[2], 11);
   m[3] = _mm_srai_epi32(m[3], 11);
   m[4] = _mm_srai_epi32(m[4], 11);
   m[5] = _mm_srai_epi32(m[5], 11);
   m[6] = _mm_srai_epi32(m[6], 11);
   m[7] = _mm_srai_epi32(m[7], 11);


   //         value = (ov_clip(value , 0,(1 << 16)-1));
   r[0] = _mm_packs_epi32(m[0], m[1]);
   r[1] = _mm_packs_epi32(m[2], m[3]);
   r[2] = _mm_packs_epi32(m[4], m[5]);
   r[3] = _mm_packs_epi32(m[6], m[7]);


   //         value = (sign ? -value : value);
   r[0] = _mm_xor_si128(r[0],_mm_cmpgt_epi16(s[0], _mm_setzero_si128()));
   r[1] = _mm_xor_si128(r[1],_mm_cmpgt_epi16(s[1], _mm_setzero_si128()));
   r[2] = _mm_xor_si128(r[2],_mm_cmpgt_epi16(s[2], _mm_setzero_si128()));
   r[3] = _mm_xor_si128(r[3],_mm_cmpgt_epi16(s[3], _mm_setzero_si128()));

   r[0] = _mm_add_epi16(r[0], s[0]);
   r[1] = _mm_add_epi16(r[1], s[1]);
   r[2] = _mm_add_epi16(r[2], s[2]);
   r[3] = _mm_add_epi16(r[3], s[3]);

   //         _dst[j] = ov_clip((int32_t)_dst[j] + value, 0, 1023);
   l[0] = _mm_add_epi16(l[0], r[0]);
   l[1] = _mm_add_epi16(l[1], r[1]);
   l[2] = _mm_add_epi16(l[2], r[2]);
   l[3] = _mm_add_epi16(l[3], r[3]);

   l[0] = _mm_max_epi16(l[0], _mm_setzero_si128());
   l[1] = _mm_max_epi16(l[1], _mm_setzero_si128());
   l[2] = _mm_max_epi16(l[2], _mm_setzero_si128());
   l[3] = _mm_max_epi16(l[3], _mm_setzero_si128());

   l[0] = _mm_min_epi16(l[0], _mm_set1_epi16(CLIP_10));
   l[1] = _mm_min_epi16(l[1], _mm_set1_epi16(CLIP_10));
   l[2] = _mm_min_epi16(l[2], _mm_set1_epi16(CLIP_10));
   l[3] = _mm_min_epi16(l[3], _mm_set1_epi16(CLIP_10));

   _mm_storeu_si128((__m128i*)&dst[0 * dst_stride], l[0]);
   _mm_storeu_si128((__m128i*)&dst[1 * dst_stride], l[1]);
   _mm_storeu_si128((__m128i*)&dst[8 + 0 * dst_stride], l[2]);
   _mm_storeu_si128((__m128i*)&dst[8 + 1 * dst_stride], l[3]);
}

static inline void
ovvc_transform_scale_add_sse_32_1_10(uint16_t *dst, ptrdiff_t dst_stride,
                               const int16_t *src, ptrdiff_t src_stride, __m128i scale)
{
   __m128i l[4], r[4], s[4], m[8];
   __m128i a = _mm_unpackhi_epi16(scale, _mm_set1_epi16(1));
   __m128i b = _mm_set1_epi16(1 << (11 - 1));

   l[0] = _mm_loadu_si128((__m128i *)&dst[0]);
   l[1] = _mm_loadu_si128((__m128i *)&dst[8]);
   l[2] = _mm_loadu_si128((__m128i *)&dst[16]);
   l[3] = _mm_loadu_si128((__m128i *)&dst[24]);
   //         value = _src[j];
   r[0] = _mm_loadu_si128((__m128i *)&src[0]);
   r[1] = _mm_loadu_si128((__m128i *)&src[8]);
   r[2] = _mm_loadu_si128((__m128i *)&src[16]);
   r[3] = _mm_loadu_si128((__m128i *)&src[24]);

   //         sign  = value & (1 << 15);
   s[0] = _mm_srai_epi16(_mm_and_si128(r[0], _mm_set1_epi16(SIGN_16)), 15);
   s[1] = _mm_srai_epi16(_mm_and_si128(r[1], _mm_set1_epi16(SIGN_16)), 15);
   s[2] = _mm_srai_epi16(_mm_and_si128(r[2], _mm_set1_epi16(SIGN_16)), 15);
   s[3] = _mm_srai_epi16(_mm_and_si128(r[3], _mm_set1_epi16(SIGN_16)), 15);

   s[0] = _mm_abs_epi16(s[0]);
   s[1] = _mm_abs_epi16(s[1]);
   s[2] = _mm_abs_epi16(s[2]);
   s[3] = _mm_abs_epi16(s[3]);

   //         value = (abs(value);
   r[0] = _mm_abs_epi16(r[0]);
   r[1] = _mm_abs_epi16(r[1]);
   r[2] = _mm_abs_epi16(r[2]);
   r[3] = _mm_abs_epi16(r[3]);

   //         value = (value * scale + (1 << (11 - 1)));
   m[0] = _mm_unpacklo_epi16(r[0], b);
   m[1] = _mm_unpackhi_epi16(r[0], b);
   m[2] = _mm_unpacklo_epi16(r[1], b);
   m[3] = _mm_unpackhi_epi16(r[1], b);
   m[4] = _mm_unpacklo_epi16(r[2], b);
   m[5] = _mm_unpackhi_epi16(r[2], b);
   m[6] = _mm_unpacklo_epi16(r[3], b);
   m[7] = _mm_unpackhi_epi16(r[3], b);

   m[0] = _mm_madd_epi16(m[0], a);
   m[1] = _mm_madd_epi16(m[1], a);
   m[2] = _mm_madd_epi16(m[2], a);
   m[3] = _mm_madd_epi16(m[3], a);
   m[4] = _mm_madd_epi16(m[4], a);
   m[5] = _mm_madd_epi16(m[5], a);
   m[6] = _mm_madd_epi16(m[6], a);
   m[7] = _mm_madd_epi16(m[7], a);

   //         value = value >> 11;
   m[0] = _mm_srai_epi32(m[0], 11);
   m[1] = _mm_srai_epi32(m[1], 11);
   m[2] = _mm_srai_epi32(m[2], 11);
   m[3] = _mm_srai_epi32(m[3], 11);
   m[4] = _mm_srai_epi32(m[4], 11);
   m[5] = _mm_srai_epi32(m[5], 11);
   m[6] = _mm_srai_epi32(m[6], 11);
   m[7] = _mm_srai_epi32(m[7], 11);


   //         value = (ov_clip(value , 0,(1 << 16)-1));
   r[0] = _mm_packs_epi32(m[0], m[1]);
   r[1] = _mm_packs_epi32(m[2], m[3]);
   r[2] = _mm_packs_epi32(m[4], m[5]);
   r[3] = _mm_packs_epi32(m[6], m[7]);


   //         value = (sign ? -value : value);
   r[0] = _mm_xor_si128(r[0],_mm_cmpgt_epi16(s[0], _mm_setzero_si128()));
   r[1] = _mm_xor_si128(r[1],_mm_cmpgt_epi16(s[1], _mm_setzero_si128()));
   r[2] = _mm_xor_si128(r[2],_mm_cmpgt_epi16(s[2], _mm_setzero_si128()));
   r[3] = _mm_xor_si128(r[3],_mm_cmpgt_epi16(s[3], _mm_setzero_si128()));

   r[0] = _mm_add_epi16(r[0], s[0]);
   r[1] = _mm_add_epi16(r[1], s[1]);
   r[2] = _mm_add_epi16(r[2], s[2]);
   r[3] = _mm_add_epi16(r[3], s[3]);

   //         _dst[j] = ov_clip((int32_t)_dst[j] + value, 0, 1023);
   l[0] = _mm_add_epi16(l[0], r[0]);
   l[1] = _mm_add_epi16(l[1], r[1]);
   l[2] = _mm_add_epi16(l[2], r[2]);
   l[3] = _mm_add_epi16(l[3], r[3]);

   l[0] = _mm_max_epi16(l[0], _mm_setzero_si128());
   l[1] = _mm_max_epi16(l[1], _mm_setzero_si128());
   l[2] = _mm_max_epi16(l[2], _mm_setzero_si128());
   l[3] = _mm_max_epi16(l[3], _mm_setzero_si128());

   l[0] = _mm_min_epi16(l[0], _mm_set1_epi16(CLIP_10));
   l[1] = _mm_min_epi16(l[1], _mm_set1_epi16(CLIP_10));
   l[2] = _mm_min_epi16(l[2], _mm_set1_epi16(CLIP_10));
   l[3] = _mm_min_epi16(l[3], _mm_set1_epi16(CLIP_10));

   _mm_storeu_si128((__m128i*)&dst[0],  l[0]);
   _mm_storeu_si128((__m128i*)&dst[8],  l[1]);
   _mm_storeu_si128((__m128i*)&dst[16], l[2]);
   _mm_storeu_si128((__m128i*)&dst[24], l[3]);
}

static inline void
ovvc_transform_scale_add_half_sse_8_4_10(uint16_t *dst, ptrdiff_t dst_stride,
                               const int16_t *src, ptrdiff_t src_stride, __m128i scale)
{
   __m128i l[4], r[4], s[4], m[8];
   __m128i a = _mm_unpackhi_epi16(scale, _mm_set1_epi16(1));
   __m128i b = _mm_set1_epi16(1 << (11 - 1));

   l[0] = _mm_loadu_si128((__m128i *)&dst[0 * dst_stride]);
   l[1] = _mm_loadu_si128((__m128i *)&dst[1 * dst_stride]);
   l[2] = _mm_loadu_si128((__m128i *)&dst[2 * dst_stride]);
   l[3] = _mm_loadu_si128((__m128i *)&dst[3 * dst_stride]);
   r[0] = _mm_loadu_si128((__m128i *)&src[0 * src_stride]);
   r[1] = _mm_loadu_si128((__m128i *)&src[1 * src_stride]);
   r[2] = _mm_loadu_si128((__m128i *)&src[2 * src_stride]);
   r[3] = _mm_loadu_si128((__m128i *)&src[3 * src_stride]);

   //         sign  = value & (1 << 15);
   s[0] = _mm_srai_epi16(_mm_and_si128(r[0], _mm_set1_epi16(SIGN_16)), 15);
   s[1] = _mm_srai_epi16(_mm_and_si128(r[1], _mm_set1_epi16(SIGN_16)), 15);
   s[2] = _mm_srai_epi16(_mm_and_si128(r[2], _mm_set1_epi16(SIGN_16)), 15);
   s[3] = _mm_srai_epi16(_mm_and_si128(r[3], _mm_set1_epi16(SIGN_16)), 15);

   s[0] = _mm_abs_epi16(s[0]);
   s[1] = _mm_abs_epi16(s[1]);
   s[2] = _mm_abs_epi16(s[2]);
   s[3] = _mm_abs_epi16(s[3]);

   //         value = (abs(value);
   r[0] = _mm_srai_epi16(r[0], 1);
   r[1] = _mm_srai_epi16(r[1], 1);
   r[2] = _mm_srai_epi16(r[2], 1);
   r[3] = _mm_srai_epi16(r[3], 1);

   r[0] = _mm_abs_epi16(r[0]);
   r[1] = _mm_abs_epi16(r[1]);
   r[2] = _mm_abs_epi16(r[2]);
   r[3] = _mm_abs_epi16(r[3]);

   //         value = (value * scale + (1 << (11 - 1)));
   m[0] = _mm_unpacklo_epi16(r[0], b);
   m[1] = _mm_unpackhi_epi16(r[0], b);
   m[2] = _mm_unpacklo_epi16(r[1], b);
   m[3] = _mm_unpackhi_epi16(r[1], b);
   m[4] = _mm_unpacklo_epi16(r[2], b);
   m[5] = _mm_unpackhi_epi16(r[2], b);
   m[6] = _mm_unpacklo_epi16(r[3], b);
   m[7] = _mm_unpackhi_epi16(r[3], b);

   m[0] = _mm_madd_epi16(m[0], a);
   m[1] = _mm_madd_epi16(m[1], a);
   m[2] = _mm_madd_epi16(m[2], a);
   m[3] = _mm_madd_epi16(m[3], a);
   m[4] = _mm_madd_epi16(m[4], a);
   m[5] = _mm_madd_epi16(m[5], a);
   m[6] = _mm_madd_epi16(m[6], a);
   m[7] = _mm_madd_epi16(m[7], a);

   //         value = value >> 11;
   m[0] = _mm_srai_epi32(m[0], 11);
   m[1] = _mm_srai_epi32(m[1], 11);
   m[2] = _mm_srai_epi32(m[2], 11);
   m[3] = _mm_srai_epi32(m[3], 11);
   m[4] = _mm_srai_epi32(m[4], 11);
   m[5] = _mm_srai_epi32(m[5], 11);
   m[6] = _mm_srai_epi32(m[6], 11);
   m[7] = _mm_srai_epi32(m[7], 11);


   //         value = (ov_clip(value , 0,(1 << 16)-1));
   r[0] = _mm_packs_epi32(m[0], m[1]);
   r[1] = _mm_packs_epi32(m[2], m[3]);
   r[2] = _mm_packs_epi32(m[4], m[5]);
   r[3] = _mm_packs_epi32(m[6], m[7]);


   //         value = (sign ? -value : value);
   r[0] = _mm_xor_si128(r[0],_mm_cmpgt_epi16(s[0], _mm_setzero_si128()));
   r[1] = _mm_xor_si128(r[1],_mm_cmpgt_epi16(s[1], _mm_setzero_si128()));
   r[2] = _mm_xor_si128(r[2],_mm_cmpgt_epi16(s[2], _mm_setzero_si128()));
   r[3] = _mm_xor_si128(r[3],_mm_cmpgt_epi16(s[3], _mm_setzero_si128()));

   r[0] = _mm_add_epi16(r[0], s[0]);
   r[1] = _mm_add_epi16(r[1], s[1]);
   r[2] = _mm_add_epi16(r[2], s[2]);
   r[3] = _mm_add_epi16(r[3], s[3]);

   //         _dst[j] = ov_clip((int32_t)_dst[j] + value, 0, 1023);
   l[0] = _mm_add_epi16(l[0], r[0]);
   l[1] = _mm_add_epi16(l[1], r[1]);
   l[2] = _mm_add_epi16(l[2], r[2]);
   l[3] = _mm_add_epi16(l[3], r[3]);

   l[0] = _mm_max_epi16(l[0], _mm_setzero_si128());
   l[1] = _mm_max_epi16(l[1], _mm_setzero_si128());
   l[2] = _mm_max_epi16(l[2], _mm_setzero_si128());
   l[3] = _mm_max_epi16(l[3], _mm_setzero_si128());

   l[0] = _mm_min_epi16(l[0], _mm_set1_epi16(CLIP_10));
   l[1] = _mm_min_epi16(l[1], _mm_set1_epi16(CLIP_10));
   l[2] = _mm_min_epi16(l[2], _mm_set1_epi16(CLIP_10));
   l[3] = _mm_min_epi16(l[3], _mm_set1_epi16(CLIP_10));

   _mm_storeu_si128((__m128i*)&dst[0 * dst_stride], l[0]);
   _mm_storeu_si128((__m128i*)&dst[1 * dst_stride], l[1]);
   _mm_storeu_si128((__m128i*)&dst[2 * dst_stride], l[2]);
   _mm_storeu_si128((__m128i*)&dst[3 * dst_stride], l[3]);
}

static inline void
ovvc_transform_scale_add_half_sse_16_2_10(uint16_t *dst, ptrdiff_t dst_stride,
                               const int16_t *src, ptrdiff_t src_stride, __m128i scale)
{
   __m128i l[4], r[4], s[4], m[8];
   __m128i a = _mm_unpackhi_epi16(scale, _mm_set1_epi16(1));
   __m128i b = _mm_set1_epi16(1 << (11 - 1));

   l[0] = _mm_loadu_si128((__m128i *)&dst[0 * dst_stride]);
   l[1] = _mm_loadu_si128((__m128i *)&dst[1 * dst_stride]);
   l[2] = _mm_loadu_si128((__m128i *)&dst[8 + 0 * dst_stride]);
   l[3] = _mm_loadu_si128((__m128i *)&dst[8 + 1 * dst_stride]);
   r[0] = _mm_loadu_si128((__m128i *)&src[0 * src_stride]);
   r[1] = _mm_loadu_si128((__m128i *)&src[1 * src_stride]);
   r[2] = _mm_loadu_si128((__m128i *)&src[8 + 0 * src_stride]);
   r[3] = _mm_loadu_si128((__m128i *)&src[8 + 1 * src_stride]);

   //         sign  = value & (1 << 15);
   s[0] = _mm_srai_epi16(_mm_and_si128(r[0], _mm_set1_epi16(SIGN_16)), 15);
   s[1] = _mm_srai_epi16(_mm_and_si128(r[1], _mm_set1_epi16(SIGN_16)), 15);
   s[2] = _mm_srai_epi16(_mm_and_si128(r[2], _mm_set1_epi16(SIGN_16)), 15);
   s[3] = _mm_srai_epi16(_mm_and_si128(r[3], _mm_set1_epi16(SIGN_16)), 15);

   s[0] = _mm_abs_epi16(s[0]);
   s[1] = _mm_abs_epi16(s[1]);
   s[2] = _mm_abs_epi16(s[2]);
   s[3] = _mm_abs_epi16(s[3]);

   //         value = (abs(value);
   r[0] = _mm_srai_epi16(r[0], 1);
   r[1] = _mm_srai_epi16(r[1], 1);
   r[2] = _mm_srai_epi16(r[2], 1);
   r[3] = _mm_srai_epi16(r[3], 1);

   r[0] = _mm_abs_epi16(r[0]);
   r[1] = _mm_abs_epi16(r[1]);
   r[2] = _mm_abs_epi16(r[2]);
   r[3] = _mm_abs_epi16(r[3]);

   //         value = (value * scale + (1 << (11 - 1)));
   m[0] = _mm_unpacklo_epi16(r[0], b);
   m[1] = _mm_unpackhi_epi16(r[0], b);
   m[2] = _mm_unpacklo_epi16(r[1], b);
   m[3] = _mm_unpackhi_epi16(r[1], b);
   m[4] = _mm_unpacklo_epi16(r[2], b);
   m[5] = _mm_unpackhi_epi16(r[2], b);
   m[6] = _mm_unpacklo_epi16(r[3], b);
   m[7] = _mm_unpackhi_epi16(r[3], b);

   m[0] = _mm_madd_epi16(m[0], a);
   m[1] = _mm_madd_epi16(m[1], a);
   m[2] = _mm_madd_epi16(m[2], a);
   m[3] = _mm_madd_epi16(m[3], a);
   m[4] = _mm_madd_epi16(m[4], a);
   m[5] = _mm_madd_epi16(m[5], a);
   m[6] = _mm_madd_epi16(m[6], a);
   m[7] = _mm_madd_epi16(m[7], a);

   //         value = value >> 11;
   m[0] = _mm_srai_epi32(m[0], 11);
   m[1] = _mm_srai_epi32(m[1], 11);
   m[2] = _mm_srai_epi32(m[2], 11);
   m[3] = _mm_srai_epi32(m[3], 11);
   m[4] = _mm_srai_epi32(m[4], 11);
   m[5] = _mm_srai_epi32(m[5], 11);
   m[6] = _mm_srai_epi32(m[6], 11);
   m[7] = _mm_srai_epi32(m[7], 11);


   //         value = (ov_clip(value , 0,(1 << 16)-1));
   r[0] = _mm_packs_epi32(m[0], m[1]);
   r[1] = _mm_packs_epi32(m[2], m[3]);
   r[2] = _mm_packs_epi32(m[4], m[5]);
   r[3] = _mm_packs_epi32(m[6], m[7]);


   //         value = (sign ? -value : value);
   r[0] = _mm_xor_si128(r[0],_mm_cmpgt_epi16(s[0], _mm_setzero_si128()));
   r[1] = _mm_xor_si128(r[1],_mm_cmpgt_epi16(s[1], _mm_setzero_si128()));
   r[2] = _mm_xor_si128(r[2],_mm_cmpgt_epi16(s[2], _mm_setzero_si128()));
   r[3] = _mm_xor_si128(r[3],_mm_cmpgt_epi16(s[3], _mm_setzero_si128()));

   r[0] = _mm_add_epi16(r[0], s[0]);
   r[1] = _mm_add_epi16(r[1], s[1]);
   r[2] = _mm_add_epi16(r[2], s[2]);
   r[3] = _mm_add_epi16(r[3], s[3]);

   //         _dst[j] = ov_clip((int32_t)_dst[j] + value, 0, 1023);
   l[0] = _mm_add_epi16(l[0], r[0]);
   l[1] = _mm_add_epi16(l[1], r[1]);
   l[2] = _mm_add_epi16(l[2], r[2]);
   l[3] = _mm_add_epi16(l[3], r[3]);

   l[0] = _mm_max_epi16(l[0], _mm_setzero_si128());
   l[1] = _mm_max_epi16(l[1], _mm_setzero_si128());
   l[2] = _mm_max_epi16(l[2], _mm_setzero_si128());
   l[3] = _mm_max_epi16(l[3], _mm_setzero_si128());

   l[0] = _mm_min_epi16(l[0], _mm_set1_epi16(CLIP_10));
   l[1] = _mm_min_epi16(l[1], _mm_set1_epi16(CLIP_10));
   l[2] = _mm_min_epi16(l[2], _mm_set1_epi16(CLIP_10));
   l[3] = _mm_min_epi16(l[3], _mm_set1_epi16(CLIP_10));

   _mm_storeu_si128((__m128i*)&dst[0 * dst_stride], l[0]);
   _mm_storeu_si128((__m128i*)&dst[1 * dst_stride], l[1]);
   _mm_storeu_si128((__m128i*)&dst[8 + 0 * dst_stride], l[2]);
   _mm_storeu_si128((__m128i*)&dst[8 + 1 * dst_stride], l[3]);
}

static inline void
ovvc_transform_scale_add_half_sse_32_1_10(uint16_t *dst, ptrdiff_t dst_stride,
                               const int16_t *src, ptrdiff_t src_stride, __m128i scale)
{
   __m128i l[4], r[4], s[4], m[8];
   __m128i a = _mm_unpackhi_epi16(scale, _mm_set1_epi16(1));
   __m128i b = _mm_set1_epi16(1 << (11 - 1));

   l[0] = _mm_loadu_si128((__m128i *)&dst[0]);
   l[1] = _mm_loadu_si128((__m128i *)&dst[8]);
   l[2] = _mm_loadu_si128((__m128i *)&dst[16]);
   l[3] = _mm_loadu_si128((__m128i *)&dst[24]);
   //         value = _src[j];
   r[0] = _mm_loadu_si128((__m128i *)&src[0]);
   r[1] = _mm_loadu_si128((__m128i *)&src[8]);
   r[2] = _mm_loadu_si128((__m128i *)&src[16]);
   r[3] = _mm_loadu_si128((__m128i *)&src[24]);

   //         sign  = value & (1 << 15);
   s[0] = _mm_srai_epi16(_mm_and_si128(r[0], _mm_set1_epi16(SIGN_16)), 15);
   s[1] = _mm_srai_epi16(_mm_and_si128(r[1], _mm_set1_epi16(SIGN_16)), 15);
   s[2] = _mm_srai_epi16(_mm_and_si128(r[2], _mm_set1_epi16(SIGN_16)), 15);
   s[3] = _mm_srai_epi16(_mm_and_si128(r[3], _mm_set1_epi16(SIGN_16)), 15);

   s[0] = _mm_abs_epi16(s[0]);
   s[1] = _mm_abs_epi16(s[1]);
   s[2] = _mm_abs_epi16(s[2]);
   s[3] = _mm_abs_epi16(s[3]);

   //         value = (abs(value);
   r[0] = _mm_srai_epi16(r[0], 1);
   r[1] = _mm_srai_epi16(r[1], 1);
   r[2] = _mm_srai_epi16(r[2], 1);
   r[3] = _mm_srai_epi16(r[3], 1);

   r[0] = _mm_abs_epi16(r[0]);
   r[1] = _mm_abs_epi16(r[1]);
   r[2] = _mm_abs_epi16(r[2]);
   r[3] = _mm_abs_epi16(r[3]);

   //         value = (value * scale + (1 << (11 - 1)));
   m[0] = _mm_unpacklo_epi16(r[0], b);
   m[1] = _mm_unpackhi_epi16(r[0], b);
   m[2] = _mm_unpacklo_epi16(r[1], b);
   m[3] = _mm_unpackhi_epi16(r[1], b);
   m[4] = _mm_unpacklo_epi16(r[2], b);
   m[5] = _mm_unpackhi_epi16(r[2], b);
   m[6] = _mm_unpacklo_epi16(r[3], b);
   m[7] = _mm_unpackhi_epi16(r[3], b);

   m[0] = _mm_madd_epi16(m[0], a);
   m[1] = _mm_madd_epi16(m[1], a);
   m[2] = _mm_madd_epi16(m[2], a);
   m[3] = _mm_madd_epi16(m[3], a);
   m[4] = _mm_madd_epi16(m[4], a);
   m[5] = _mm_madd_epi16(m[5], a);
   m[6] = _mm_madd_epi16(m[6], a);
   m[7] = _mm_madd_epi16(m[7], a);

   //         value = value >> 11;
   m[0] = _mm_srai_epi32(m[0], 11);
   m[1] = _mm_srai_epi32(m[1], 11);
   m[2] = _mm_srai_epi32(m[2], 11);
   m[3] = _mm_srai_epi32(m[3], 11);
   m[4] = _mm_srai_epi32(m[4], 11);
   m[5] = _mm_srai_epi32(m[5], 11);
   m[6] = _mm_srai_epi32(m[6], 11);
   m[7] = _mm_srai_epi32(m[7], 11);


   //         value = (ov_clip(value , 0,(1 << 16)-1));
   r[0] = _mm_packs_epi32(m[0], m[1]);
   r[1] = _mm_packs_epi32(m[2], m[3]);
   r[2] = _mm_packs_epi32(m[4], m[5]);
   r[3] = _mm_packs_epi32(m[6], m[7]);


   //         value = (sign ? -value : value);
   r[0] = _mm_xor_si128(r[0],_mm_cmpgt_epi16(s[0], _mm_setzero_si128()));
   r[1] = _mm_xor_si128(r[1],_mm_cmpgt_epi16(s[1], _mm_setzero_si128()));
   r[2] = _mm_xor_si128(r[2],_mm_cmpgt_epi16(s[2], _mm_setzero_si128()));
   r[3] = _mm_xor_si128(r[3],_mm_cmpgt_epi16(s[3], _mm_setzero_si128()));

   r[0] = _mm_add_epi16(r[0], s[0]);
   r[1] = _mm_add_epi16(r[1], s[1]);
   r[2] = _mm_add_epi16(r[2], s[2]);
   r[3] = _mm_add_epi16(r[3], s[3]);

   //         _dst[j] = ov_clip((int32_t)_dst[j] + value, 0, 1023);
   l[0] = _mm_add_epi16(l[0], r[0]);
   l[1] = _mm_add_epi16(l[1], r[1]);
   l[2] = _mm_add_epi16(l[2], r[2]);
   l[3] = _mm_add_epi16(l[3], r[3]);

   l[0] = _mm_max_epi16(l[0], _mm_setzero_si128());
   l[1] = _mm_max_epi16(l[1], _mm_setzero_si128());
   l[2] = _mm_max_epi16(l[2], _mm_setzero_si128());
   l[3] = _mm_max_epi16(l[3], _mm_setzero_si128());

   l[0] = _mm_min_epi16(l[0], _mm_set1_epi16(CLIP_10));
   l[1] = _mm_min_epi16(l[1], _mm_set1_epi16(CLIP_10));
   l[2] = _mm_min_epi16(l[2], _mm_set1_epi16(CLIP_10));
   l[3] = _mm_min_epi16(l[3], _mm_set1_epi16(CLIP_10));

   _mm_storeu_si128((__m128i*)&dst[0],  l[0]);
   _mm_storeu_si128((__m128i*)&dst[8],  l[1]);
   _mm_storeu_si128((__m128i*)&dst[16], l[2]);
   _mm_storeu_si128((__m128i*)&dst[24], l[3]);
}

static inline void
ovvc_transform_scale_sub_sse_8_4_10(uint16_t *dst, ptrdiff_t dst_stride,
                               const int16_t *src, ptrdiff_t src_stride, __m128i scale)
{
   __m128i l[4], r[4], s[4], m[8];
   __m128i a = _mm_unpackhi_epi16(scale, _mm_set1_epi16(1));
   __m128i b = _mm_set1_epi16(1 << (11 - 1));

   l[0] = _mm_loadu_si128((__m128i *)&dst[0 * dst_stride]);
   l[1] = _mm_loadu_si128((__m128i *)&dst[1 * dst_stride]);
   l[2] = _mm_loadu_si128((__m128i *)&dst[2 * dst_stride]);
   l[3] = _mm_loadu_si128((__m128i *)&dst[3 * dst_stride]);
   r[0] = _mm_loadu_si128((__m128i *)&src[0 * src_stride]);
   r[1] = _mm_loadu_si128((__m128i *)&src[1 * src_stride]);
   r[2] = _mm_loadu_si128((__m128i *)&src[2 * src_stride]);
   r[3] = _mm_loadu_si128((__m128i *)&src[3 * src_stride]);

   //         sign  = value & (1 << 15);
   s[0] = _mm_srai_epi16(_mm_and_si128(r[0], _mm_set1_epi16(SIGN_16)), 15);
   s[1] = _mm_srai_epi16(_mm_and_si128(r[1], _mm_set1_epi16(SIGN_16)), 15);
   s[2] = _mm_srai_epi16(_mm_and_si128(r[2], _mm_set1_epi16(SIGN_16)), 15);
   s[3] = _mm_srai_epi16(_mm_and_si128(r[3], _mm_set1_epi16(SIGN_16)), 15);

   s[0] = _mm_andnot_si128(_mm_abs_epi16(s[0]), _mm_set1_epi16(1));
   s[1] = _mm_andnot_si128(_mm_abs_epi16(s[1]), _mm_set1_epi16(1));
   s[2] = _mm_andnot_si128(_mm_abs_epi16(s[2]), _mm_set1_epi16(1));
   s[3] = _mm_andnot_si128(_mm_abs_epi16(s[3]), _mm_set1_epi16(1));

   //         value = (abs(value);
   r[0] = _mm_abs_epi16(r[0]);
   r[1] = _mm_abs_epi16(r[1]);
   r[2] = _mm_abs_epi16(r[2]);
   r[3] = _mm_abs_epi16(r[3]);

   //         value = (value * scale + (1 << (11 - 1)));
   m[0] = _mm_unpacklo_epi16(r[0], b);
   m[1] = _mm_unpackhi_epi16(r[0], b);
   m[2] = _mm_unpacklo_epi16(r[1], b);
   m[3] = _mm_unpackhi_epi16(r[1], b);
   m[4] = _mm_unpacklo_epi16(r[2], b);
   m[5] = _mm_unpackhi_epi16(r[2], b);
   m[6] = _mm_unpacklo_epi16(r[3], b);
   m[7] = _mm_unpackhi_epi16(r[3], b);

   m[0] = _mm_madd_epi16(m[0], a);
   m[1] = _mm_madd_epi16(m[1], a);
   m[2] = _mm_madd_epi16(m[2], a);
   m[3] = _mm_madd_epi16(m[3], a);
   m[4] = _mm_madd_epi16(m[4], a);
   m[5] = _mm_madd_epi16(m[5], a);
   m[6] = _mm_madd_epi16(m[6], a);
   m[7] = _mm_madd_epi16(m[7], a);

   //         value = value >> 11;
   m[0] = _mm_srai_epi32(m[0], 11);
   m[1] = _mm_srai_epi32(m[1], 11);
   m[2] = _mm_srai_epi32(m[2], 11);
   m[3] = _mm_srai_epi32(m[3], 11);
   m[4] = _mm_srai_epi32(m[4], 11);
   m[5] = _mm_srai_epi32(m[5], 11);
   m[6] = _mm_srai_epi32(m[6], 11);
   m[7] = _mm_srai_epi32(m[7], 11);


   //         value = (ov_clip(value , 0,(1 << 16)-1));
   r[0] = _mm_packs_epi32(m[0], m[1]);
   r[1] = _mm_packs_epi32(m[2], m[3]);
   r[2] = _mm_packs_epi32(m[4], m[5]);
   r[3] = _mm_packs_epi32(m[6], m[7]);


   //         value = (sign ? value : -value);
   r[0] = _mm_xor_si128(r[0],_mm_cmpgt_epi16(s[0], _mm_setzero_si128()));
   r[1] = _mm_xor_si128(r[1],_mm_cmpgt_epi16(s[1], _mm_setzero_si128()));
   r[2] = _mm_xor_si128(r[2],_mm_cmpgt_epi16(s[2], _mm_setzero_si128()));
   r[3] = _mm_xor_si128(r[3],_mm_cmpgt_epi16(s[3], _mm_setzero_si128()));

   r[0] = _mm_add_epi16(r[0], s[0]);
   r[1] = _mm_add_epi16(r[1], s[1]);
   r[2] = _mm_add_epi16(r[2], s[2]);
   r[3] = _mm_add_epi16(r[3], s[3]);

   //         _dst[j] = ov_clip((int32_t)_dst[j] + value, 0, 1023);
   l[0] = _mm_add_epi16(l[0], r[0]);
   l[1] = _mm_add_epi16(l[1], r[1]);
   l[2] = _mm_add_epi16(l[2], r[2]);
   l[3] = _mm_add_epi16(l[3], r[3]);

   l[0] = _mm_max_epi16(l[0], _mm_setzero_si128());
   l[1] = _mm_max_epi16(l[1], _mm_setzero_si128());
   l[2] = _mm_max_epi16(l[2], _mm_setzero_si128());
   l[3] = _mm_max_epi16(l[3], _mm_setzero_si128());

   l[0] = _mm_min_epi16(l[0], _mm_set1_epi16(CLIP_10));
   l[1] = _mm_min_epi16(l[1], _mm_set1_epi16(CLIP_10));
   l[2] = _mm_min_epi16(l[2], _mm_set1_epi16(CLIP_10));
   l[3] = _mm_min_epi16(l[3], _mm_set1_epi16(CLIP_10));

   _mm_storeu_si128((__m128i*)&dst[0 * dst_stride], l[0]);
   _mm_storeu_si128((__m128i*)&dst[1 * dst_stride], l[1]);
   _mm_storeu_si128((__m128i*)&dst[2 * dst_stride], l[2]);
   _mm_storeu_si128((__m128i*)&dst[3 * dst_stride], l[3]);
}

static inline void
ovvc_transform_scale_sub_sse_16_2_10(uint16_t *dst, ptrdiff_t dst_stride,
                               const int16_t *src, ptrdiff_t src_stride, __m128i scale)
{
   __m128i l[4], r[4], s[4], m[8];
   __m128i a = _mm_unpackhi_epi16(scale, _mm_set1_epi16(1));
   __m128i b = _mm_set1_epi16(1 << (11 - 1));

   l[0] = _mm_loadu_si128((__m128i *)&dst[0 * dst_stride]);
   l[1] = _mm_loadu_si128((__m128i *)&dst[1 * dst_stride]);
   l[2] = _mm_loadu_si128((__m128i *)&dst[8 + 0 * dst_stride]);
   l[3] = _mm_loadu_si128((__m128i *)&dst[8 + 1 * dst_stride]);
   r[0] = _mm_loadu_si128((__m128i *)&src[0 * src_stride]);
   r[1] = _mm_loadu_si128((__m128i *)&src[1 * src_stride]);
   r[2] = _mm_loadu_si128((__m128i *)&src[8 + 0 * src_stride]);
   r[3] = _mm_loadu_si128((__m128i *)&src[8 + 1 * src_stride]);

   //         sign  = value & (1 << 15);
   s[0] = _mm_srai_epi16(_mm_and_si128(r[0], _mm_set1_epi16(SIGN_16)), 15);
   s[1] = _mm_srai_epi16(_mm_and_si128(r[1], _mm_set1_epi16(SIGN_16)), 15);
   s[2] = _mm_srai_epi16(_mm_and_si128(r[2], _mm_set1_epi16(SIGN_16)), 15);
   s[3] = _mm_srai_epi16(_mm_and_si128(r[3], _mm_set1_epi16(SIGN_16)), 15);

   s[0] = _mm_andnot_si128(_mm_abs_epi16(s[0]), _mm_set1_epi16(1));
   s[1] = _mm_andnot_si128(_mm_abs_epi16(s[1]), _mm_set1_epi16(1));
   s[2] = _mm_andnot_si128(_mm_abs_epi16(s[2]), _mm_set1_epi16(1));
   s[3] = _mm_andnot_si128(_mm_abs_epi16(s[3]), _mm_set1_epi16(1));

   //         value = (abs(value);
   r[0] = _mm_abs_epi16(r[0]);
   r[1] = _mm_abs_epi16(r[1]);
   r[2] = _mm_abs_epi16(r[2]);
   r[3] = _mm_abs_epi16(r[3]);

   //         value = (value * scale + (1 << (11 - 1)));
   m[0] = _mm_unpacklo_epi16(r[0], b);
   m[1] = _mm_unpackhi_epi16(r[0], b);
   m[2] = _mm_unpacklo_epi16(r[1], b);
   m[3] = _mm_unpackhi_epi16(r[1], b);
   m[4] = _mm_unpacklo_epi16(r[2], b);
   m[5] = _mm_unpackhi_epi16(r[2], b);
   m[6] = _mm_unpacklo_epi16(r[3], b);
   m[7] = _mm_unpackhi_epi16(r[3], b);

   m[0] = _mm_madd_epi16(m[0], a);
   m[1] = _mm_madd_epi16(m[1], a);
   m[2] = _mm_madd_epi16(m[2], a);
   m[3] = _mm_madd_epi16(m[3], a);
   m[4] = _mm_madd_epi16(m[4], a);
   m[5] = _mm_madd_epi16(m[5], a);
   m[6] = _mm_madd_epi16(m[6], a);
   m[7] = _mm_madd_epi16(m[7], a);

   //         value = value >> 11;
   m[0] = _mm_srai_epi32(m[0], 11);
   m[1] = _mm_srai_epi32(m[1], 11);
   m[2] = _mm_srai_epi32(m[2], 11);
   m[3] = _mm_srai_epi32(m[3], 11);
   m[4] = _mm_srai_epi32(m[4], 11);
   m[5] = _mm_srai_epi32(m[5], 11);
   m[6] = _mm_srai_epi32(m[6], 11);
   m[7] = _mm_srai_epi32(m[7], 11);


   //         value = (ov_clip(value , 0,(1 << 16)-1));
   r[0] = _mm_packs_epi32(m[0], m[1]);
   r[1] = _mm_packs_epi32(m[2], m[3]);
   r[2] = _mm_packs_epi32(m[4], m[5]);
   r[3] = _mm_packs_epi32(m[6], m[7]);


   //         value = (sign ? -value : value);
   r[0] = _mm_xor_si128(r[0],_mm_cmpgt_epi16(s[0], _mm_setzero_si128()));
   r[1] = _mm_xor_si128(r[1],_mm_cmpgt_epi16(s[1], _mm_setzero_si128()));
   r[2] = _mm_xor_si128(r[2],_mm_cmpgt_epi16(s[2], _mm_setzero_si128()));
   r[3] = _mm_xor_si128(r[3],_mm_cmpgt_epi16(s[3], _mm_setzero_si128()));

   r[0] = _mm_add_epi16(r[0], s[0]);
   r[1] = _mm_add_epi16(r[1], s[1]);
   r[2] = _mm_add_epi16(r[2], s[2]);
   r[3] = _mm_add_epi16(r[3], s[3]);

   //         _dst[j] = ov_clip((int32_t)_dst[j] + value, 0, 1023);
   l[0] = _mm_add_epi16(l[0], r[0]);
   l[1] = _mm_add_epi16(l[1], r[1]);
   l[2] = _mm_add_epi16(l[2], r[2]);
   l[3] = _mm_add_epi16(l[3], r[3]);

   l[0] = _mm_max_epi16(l[0], _mm_setzero_si128());
   l[1] = _mm_max_epi16(l[1], _mm_setzero_si128());
   l[2] = _mm_max_epi16(l[2], _mm_setzero_si128());
   l[3] = _mm_max_epi16(l[3], _mm_setzero_si128());

   l[0] = _mm_min_epi16(l[0], _mm_set1_epi16(CLIP_10));
   l[1] = _mm_min_epi16(l[1], _mm_set1_epi16(CLIP_10));
   l[2] = _mm_min_epi16(l[2], _mm_set1_epi16(CLIP_10));
   l[3] = _mm_min_epi16(l[3], _mm_set1_epi16(CLIP_10));

   _mm_storeu_si128((__m128i*)&dst[0 * dst_stride], l[0]);
   _mm_storeu_si128((__m128i*)&dst[1 * dst_stride], l[1]);
   _mm_storeu_si128((__m128i*)&dst[8 + 0 * dst_stride], l[2]);
   _mm_storeu_si128((__m128i*)&dst[8 + 1 * dst_stride], l[3]);
}

static inline void
ovvc_transform_scale_sub_sse_32_1_10(uint16_t *dst, ptrdiff_t dst_stride,
                               const int16_t *src, ptrdiff_t src_stride, __m128i scale)
{
   __m128i l[4], r[4], s[4], m[8];
   __m128i a = _mm_unpackhi_epi16(scale, _mm_set1_epi16(1));
   __m128i b = _mm_set1_epi16(1 << (11 - 1));

   l[0] = _mm_loadu_si128((__m128i *)&dst[0]);
   l[1] = _mm_loadu_si128((__m128i *)&dst[8]);
   l[2] = _mm_loadu_si128((__m128i *)&dst[16]);
   l[3] = _mm_loadu_si128((__m128i *)&dst[24]);
   //         value = _src[j];
   r[0] = _mm_loadu_si128((__m128i *)&src[0]);
   r[1] = _mm_loadu_si128((__m128i *)&src[8]);
   r[2] = _mm_loadu_si128((__m128i *)&src[16]);
   r[3] = _mm_loadu_si128((__m128i *)&src[24]);

   //         sign  = value & (1 << 15);
   s[0] = _mm_srai_epi16(_mm_and_si128(r[0], _mm_set1_epi16(SIGN_16)), 15);
   s[1] = _mm_srai_epi16(_mm_and_si128(r[1], _mm_set1_epi16(SIGN_16)), 15);
   s[2] = _mm_srai_epi16(_mm_and_si128(r[2], _mm_set1_epi16(SIGN_16)), 15);
   s[3] = _mm_srai_epi16(_mm_and_si128(r[3], _mm_set1_epi16(SIGN_16)), 15);

   s[0] = _mm_andnot_si128(_mm_abs_epi16(s[0]), _mm_set1_epi16(1));
   s[1] = _mm_andnot_si128(_mm_abs_epi16(s[1]), _mm_set1_epi16(1));
   s[2] = _mm_andnot_si128(_mm_abs_epi16(s[2]), _mm_set1_epi16(1));
   s[3] = _mm_andnot_si128(_mm_abs_epi16(s[3]), _mm_set1_epi16(1));

   //         value = (abs(value);
   r[0] = _mm_abs_epi16(r[0]);
   r[1] = _mm_abs_epi16(r[1]);
   r[2] = _mm_abs_epi16(r[2]);
   r[3] = _mm_abs_epi16(r[3]);

   //         value = (value * scale + (1 << (11 - 1)));
   m[0] = _mm_unpacklo_epi16(r[0], b);
   m[1] = _mm_unpackhi_epi16(r[0], b);
   m[2] = _mm_unpacklo_epi16(r[1], b);
   m[3] = _mm_unpackhi_epi16(r[1], b);
   m[4] = _mm_unpacklo_epi16(r[2], b);
   m[5] = _mm_unpackhi_epi16(r[2], b);
   m[6] = _mm_unpacklo_epi16(r[3], b);
   m[7] = _mm_unpackhi_epi16(r[3], b);

   m[0] = _mm_madd_epi16(m[0], a);
   m[1] = _mm_madd_epi16(m[1], a);
   m[2] = _mm_madd_epi16(m[2], a);
   m[3] = _mm_madd_epi16(m[3], a);
   m[4] = _mm_madd_epi16(m[4], a);
   m[5] = _mm_madd_epi16(m[5], a);
   m[6] = _mm_madd_epi16(m[6], a);
   m[7] = _mm_madd_epi16(m[7], a);

   //         value = value >> 11;
   m[0] = _mm_srai_epi32(m[0], 11);
   m[1] = _mm_srai_epi32(m[1], 11);
   m[2] = _mm_srai_epi32(m[2], 11);
   m[3] = _mm_srai_epi32(m[3], 11);
   m[4] = _mm_srai_epi32(m[4], 11);
   m[5] = _mm_srai_epi32(m[5], 11);
   m[6] = _mm_srai_epi32(m[6], 11);
   m[7] = _mm_srai_epi32(m[7], 11);


   //         value = (ov_clip(value , 0,(1 << 16)-1));
   r[0] = _mm_packs_epi32(m[0], m[1]);
   r[1] = _mm_packs_epi32(m[2], m[3]);
   r[2] = _mm_packs_epi32(m[4], m[5]);
   r[3] = _mm_packs_epi32(m[6], m[7]);


   //         value = (sign ? -value : value);
   r[0] = _mm_xor_si128(r[0],_mm_cmpgt_epi16(s[0], _mm_setzero_si128()));
   r[1] = _mm_xor_si128(r[1],_mm_cmpgt_epi16(s[1], _mm_setzero_si128()));
   r[2] = _mm_xor_si128(r[2],_mm_cmpgt_epi16(s[2], _mm_setzero_si128()));
   r[3] = _mm_xor_si128(r[3],_mm_cmpgt_epi16(s[3], _mm_setzero_si128()));

   r[0] = _mm_add_epi16(r[0], s[0]);
   r[1] = _mm_add_epi16(r[1], s[1]);
   r[2] = _mm_add_epi16(r[2], s[2]);
   r[3] = _mm_add_epi16(r[3], s[3]);

   //         _dst[j] = ov_clip((int32_t)_dst[j] + value, 0, 1023);
   l[0] = _mm_add_epi16(l[0], r[0]);
   l[1] = _mm_add_epi16(l[1], r[1]);
   l[2] = _mm_add_epi16(l[2], r[2]);
   l[3] = _mm_add_epi16(l[3], r[3]);

   l[0] = _mm_max_epi16(l[0], _mm_setzero_si128());
   l[1] = _mm_max_epi16(l[1], _mm_setzero_si128());
   l[2] = _mm_max_epi16(l[2], _mm_setzero_si128());
   l[3] = _mm_max_epi16(l[3], _mm_setzero_si128());

   l[0] = _mm_min_epi16(l[0], _mm_set1_epi16(CLIP_10));
   l[1] = _mm_min_epi16(l[1], _mm_set1_epi16(CLIP_10));
   l[2] = _mm_min_epi16(l[2], _mm_set1_epi16(CLIP_10));
   l[3] = _mm_min_epi16(l[3], _mm_set1_epi16(CLIP_10));

   _mm_storeu_si128((__m128i*)&dst[0],  l[0]);
   _mm_storeu_si128((__m128i*)&dst[8],  l[1]);
   _mm_storeu_si128((__m128i*)&dst[16], l[2]);
   _mm_storeu_si128((__m128i*)&dst[24], l[3]);
}

static inline void
ovvc_transform_scale_sub_half_sse_8_4_10(uint16_t *dst, ptrdiff_t dst_stride,
                               const int16_t *src, ptrdiff_t src_stride, __m128i scale)
{
   __m128i l[4], r[4], s[4], m[8];
   __m128i a = _mm_unpackhi_epi16(scale, _mm_set1_epi16(1));
   __m128i b = _mm_set1_epi16(1 << (11 - 1));

   l[0] = _mm_loadu_si128((__m128i *)&dst[0 * dst_stride]);
   l[1] = _mm_loadu_si128((__m128i *)&dst[1 * dst_stride]);
   l[2] = _mm_loadu_si128((__m128i *)&dst[2 * dst_stride]);
   l[3] = _mm_loadu_si128((__m128i *)&dst[3 * dst_stride]);
   r[0] = _mm_loadu_si128((__m128i *)&src[0 * src_stride]);
   r[1] = _mm_loadu_si128((__m128i *)&src[1 * src_stride]);
   r[2] = _mm_loadu_si128((__m128i *)&src[2 * src_stride]);
   r[3] = _mm_loadu_si128((__m128i *)&src[3 * src_stride]);

   //         sign  = value & (1 << 15);
   s[0] = _mm_srai_epi16(_mm_and_si128(r[0], _mm_set1_epi16(SIGN_16)), 15);
   s[1] = _mm_srai_epi16(_mm_and_si128(r[1], _mm_set1_epi16(SIGN_16)), 15);
   s[2] = _mm_srai_epi16(_mm_and_si128(r[2], _mm_set1_epi16(SIGN_16)), 15);
   s[3] = _mm_srai_epi16(_mm_and_si128(r[3], _mm_set1_epi16(SIGN_16)), 15);

   s[0] = _mm_andnot_si128(_mm_abs_epi16(s[0]), _mm_set1_epi16(1));
   s[1] = _mm_andnot_si128(_mm_abs_epi16(s[1]), _mm_set1_epi16(1));
   s[2] = _mm_andnot_si128(_mm_abs_epi16(s[2]), _mm_set1_epi16(1));
   s[3] = _mm_andnot_si128(_mm_abs_epi16(s[3]), _mm_set1_epi16(1));

   //         value = (-value);
   r[0] = _mm_xor_si128(r[0], _mm_set1_epi16((int16_t)0xFFFF));
   r[1] = _mm_xor_si128(r[1], _mm_set1_epi16((int16_t)0xFFFF));
   r[2] = _mm_xor_si128(r[2], _mm_set1_epi16((int16_t)0xFFFF));
   r[3] = _mm_xor_si128(r[3], _mm_set1_epi16((int16_t)0xFFFF));

   r[0] = _mm_add_epi16(r[0], _mm_set1_epi16(1));
   r[1] = _mm_add_epi16(r[1], _mm_set1_epi16(1));
   r[2] = _mm_add_epi16(r[2], _mm_set1_epi16(1));
   r[3] = _mm_add_epi16(r[3], _mm_set1_epi16(1));

   //         value = (abs(value);
   r[0] = _mm_srai_epi16(r[0], 1);
   r[1] = _mm_srai_epi16(r[1], 1);
   r[2] = _mm_srai_epi16(r[2], 1);
   r[3] = _mm_srai_epi16(r[3], 1);

   r[0] = _mm_abs_epi16(r[0]);
   r[1] = _mm_abs_epi16(r[1]);
   r[2] = _mm_abs_epi16(r[2]);
   r[3] = _mm_abs_epi16(r[3]);

   //         value = (value * scale + (1 << (11 - 1)));
   m[0] = _mm_unpacklo_epi16(r[0], b);
   m[1] = _mm_unpackhi_epi16(r[0], b);
   m[2] = _mm_unpacklo_epi16(r[1], b);
   m[3] = _mm_unpackhi_epi16(r[1], b);
   m[4] = _mm_unpacklo_epi16(r[2], b);
   m[5] = _mm_unpackhi_epi16(r[2], b);
   m[6] = _mm_unpacklo_epi16(r[3], b);
   m[7] = _mm_unpackhi_epi16(r[3], b);

   m[0] = _mm_madd_epi16(m[0], a);
   m[1] = _mm_madd_epi16(m[1], a);
   m[2] = _mm_madd_epi16(m[2], a);
   m[3] = _mm_madd_epi16(m[3], a);
   m[4] = _mm_madd_epi16(m[4], a);
   m[5] = _mm_madd_epi16(m[5], a);
   m[6] = _mm_madd_epi16(m[6], a);
   m[7] = _mm_madd_epi16(m[7], a);

   //         value = value >> 11;
   m[0] = _mm_srai_epi32(m[0], 11);
   m[1] = _mm_srai_epi32(m[1], 11);
   m[2] = _mm_srai_epi32(m[2], 11);
   m[3] = _mm_srai_epi32(m[3], 11);
   m[4] = _mm_srai_epi32(m[4], 11);
   m[5] = _mm_srai_epi32(m[5], 11);
   m[6] = _mm_srai_epi32(m[6], 11);
   m[7] = _mm_srai_epi32(m[7], 11);


   //         value = (ov_clip(value , 0,(1 << 16)-1));
   r[0] = _mm_packs_epi32(m[0], m[1]);
   r[1] = _mm_packs_epi32(m[2], m[3]);
   r[2] = _mm_packs_epi32(m[4], m[5]);
   r[3] = _mm_packs_epi32(m[6], m[7]);

   //         value = (sign ? -value : value);
   r[0] = _mm_xor_si128(r[0],_mm_cmpgt_epi16(s[0], _mm_setzero_si128()));
   r[1] = _mm_xor_si128(r[1],_mm_cmpgt_epi16(s[1], _mm_setzero_si128()));
   r[2] = _mm_xor_si128(r[2],_mm_cmpgt_epi16(s[2], _mm_setzero_si128()));
   r[3] = _mm_xor_si128(r[3],_mm_cmpgt_epi16(s[3], _mm_setzero_si128()));

   r[0] = _mm_add_epi16(r[0], s[0]);
   r[1] = _mm_add_epi16(r[1], s[1]);
   r[2] = _mm_add_epi16(r[2], s[2]);
   r[3] = _mm_add_epi16(r[3], s[3]);

   //         _dst[j] = ov_clip((int32_t)_dst[j] + value, 0, 1023);
   l[0] = _mm_add_epi16(l[0], r[0]);
   l[1] = _mm_add_epi16(l[1], r[1]);
   l[2] = _mm_add_epi16(l[2], r[2]);
   l[3] = _mm_add_epi16(l[3], r[3]);

   l[0] = _mm_max_epi16(l[0], _mm_setzero_si128());
   l[1] = _mm_max_epi16(l[1], _mm_setzero_si128());
   l[2] = _mm_max_epi16(l[2], _mm_setzero_si128());
   l[3] = _mm_max_epi16(l[3], _mm_setzero_si128());

   l[0] = _mm_min_epi16(l[0], _mm_set1_epi16(CLIP_10));
   l[1] = _mm_min_epi16(l[1], _mm_set1_epi16(CLIP_10));
   l[2] = _mm_min_epi16(l[2], _mm_set1_epi16(CLIP_10));
   l[3] = _mm_min_epi16(l[3], _mm_set1_epi16(CLIP_10));

   _mm_storeu_si128((__m128i*)&dst[0 * dst_stride], l[0]);
   _mm_storeu_si128((__m128i*)&dst[1 * dst_stride], l[1]);
   _mm_storeu_si128((__m128i*)&dst[2 * dst_stride], l[2]);
   _mm_storeu_si128((__m128i*)&dst[3 * dst_stride], l[3]);
}

static inline void
ovvc_transform_scale_sub_half_sse_16_2_10(uint16_t *dst, ptrdiff_t dst_stride,
                               const int16_t *src, ptrdiff_t src_stride, __m128i scale)
{
   __m128i l[4], r[4], s[4], m[8];
   __m128i a = _mm_unpackhi_epi16(scale, _mm_set1_epi16(1));
   __m128i b = _mm_set1_epi16(1 << (11 - 1));

   l[0] = _mm_loadu_si128((__m128i *)&dst[0 * dst_stride]);
   l[1] = _mm_loadu_si128((__m128i *)&dst[1 * dst_stride]);
   l[2] = _mm_loadu_si128((__m128i *)&dst[8 + 0 * dst_stride]);
   l[3] = _mm_loadu_si128((__m128i *)&dst[8 + 1 * dst_stride]);
   r[0] = _mm_loadu_si128((__m128i *)&src[0 * src_stride]);
   r[1] = _mm_loadu_si128((__m128i *)&src[1 * src_stride]);
   r[2] = _mm_loadu_si128((__m128i *)&src[8 + 0 * src_stride]);
   r[3] = _mm_loadu_si128((__m128i *)&src[8 + 1 * src_stride]);

   //         sign  = value & (1 << 15);
   s[0] = _mm_srai_epi16(_mm_and_si128(r[0], _mm_set1_epi16(SIGN_16)), 15);
   s[1] = _mm_srai_epi16(_mm_and_si128(r[1], _mm_set1_epi16(SIGN_16)), 15);
   s[2] = _mm_srai_epi16(_mm_and_si128(r[2], _mm_set1_epi16(SIGN_16)), 15);
   s[3] = _mm_srai_epi16(_mm_and_si128(r[3], _mm_set1_epi16(SIGN_16)), 15);

   s[0] = _mm_andnot_si128(_mm_abs_epi16(s[0]), _mm_set1_epi16(1));
   s[1] = _mm_andnot_si128(_mm_abs_epi16(s[1]), _mm_set1_epi16(1));
   s[2] = _mm_andnot_si128(_mm_abs_epi16(s[2]), _mm_set1_epi16(1));
   s[3] = _mm_andnot_si128(_mm_abs_epi16(s[3]), _mm_set1_epi16(1));

   //         value = (-value);
   r[0] = _mm_xor_si128(r[0], _mm_set1_epi16((int16_t)0xFFFF));
   r[1] = _mm_xor_si128(r[1], _mm_set1_epi16((int16_t)0xFFFF));
   r[2] = _mm_xor_si128(r[2], _mm_set1_epi16((int16_t)0xFFFF));
   r[3] = _mm_xor_si128(r[3], _mm_set1_epi16((int16_t)0xFFFF));

   r[0] = _mm_add_epi16(r[0], _mm_set1_epi16(1));
   r[1] = _mm_add_epi16(r[1], _mm_set1_epi16(1));
   r[2] = _mm_add_epi16(r[2], _mm_set1_epi16(1));
   r[3] = _mm_add_epi16(r[3], _mm_set1_epi16(1));

   //         value = (abs(value);
   r[0] = _mm_srai_epi16(r[0], 1);
   r[1] = _mm_srai_epi16(r[1], 1);
   r[2] = _mm_srai_epi16(r[2], 1);
   r[3] = _mm_srai_epi16(r[3], 1);

   r[0] = _mm_abs_epi16(r[0]);
   r[1] = _mm_abs_epi16(r[1]);
   r[2] = _mm_abs_epi16(r[2]);
   r[3] = _mm_abs_epi16(r[3]);

   //         value = (value * scale + (1 << (11 - 1)));
   m[0] = _mm_unpacklo_epi16(r[0], b);
   m[1] = _mm_unpackhi_epi16(r[0], b);
   m[2] = _mm_unpacklo_epi16(r[1], b);
   m[3] = _mm_unpackhi_epi16(r[1], b);
   m[4] = _mm_unpacklo_epi16(r[2], b);
   m[5] = _mm_unpackhi_epi16(r[2], b);
   m[6] = _mm_unpacklo_epi16(r[3], b);
   m[7] = _mm_unpackhi_epi16(r[3], b);

   m[0] = _mm_madd_epi16(m[0], a);
   m[1] = _mm_madd_epi16(m[1], a);
   m[2] = _mm_madd_epi16(m[2], a);
   m[3] = _mm_madd_epi16(m[3], a);
   m[4] = _mm_madd_epi16(m[4], a);
   m[5] = _mm_madd_epi16(m[5], a);
   m[6] = _mm_madd_epi16(m[6], a);
   m[7] = _mm_madd_epi16(m[7], a);

   //         value = value >> 11;
   m[0] = _mm_srai_epi32(m[0], 11);
   m[1] = _mm_srai_epi32(m[1], 11);
   m[2] = _mm_srai_epi32(m[2], 11);
   m[3] = _mm_srai_epi32(m[3], 11);
   m[4] = _mm_srai_epi32(m[4], 11);
   m[5] = _mm_srai_epi32(m[5], 11);
   m[6] = _mm_srai_epi32(m[6], 11);
   m[7] = _mm_srai_epi32(m[7], 11);


   //         value = (ov_clip(value , 0,(1 << 16)-1));
   r[0] = _mm_packs_epi32(m[0], m[1]);
   r[1] = _mm_packs_epi32(m[2], m[3]);
   r[2] = _mm_packs_epi32(m[4], m[5]);
   r[3] = _mm_packs_epi32(m[6], m[7]);


   //         value = (sign ? -value : value);
   r[0] = _mm_xor_si128(r[0],_mm_cmpgt_epi16(s[0], _mm_setzero_si128()));
   r[1] = _mm_xor_si128(r[1],_mm_cmpgt_epi16(s[1], _mm_setzero_si128()));
   r[2] = _mm_xor_si128(r[2],_mm_cmpgt_epi16(s[2], _mm_setzero_si128()));
   r[3] = _mm_xor_si128(r[3],_mm_cmpgt_epi16(s[3], _mm_setzero_si128()));

   r[0] = _mm_add_epi16(r[0], s[0]);
   r[1] = _mm_add_epi16(r[1], s[1]);
   r[2] = _mm_add_epi16(r[2], s[2]);
   r[3] = _mm_add_epi16(r[3], s[3]);

   //         _dst[j] = ov_clip((int32_t)_dst[j] + value, 0, 1023);
   l[0] = _mm_add_epi16(l[0], r[0]);
   l[1] = _mm_add_epi16(l[1], r[1]);
   l[2] = _mm_add_epi16(l[2], r[2]);
   l[3] = _mm_add_epi16(l[3], r[3]);

   l[0] = _mm_max_epi16(l[0], _mm_setzero_si128());
   l[1] = _mm_max_epi16(l[1], _mm_setzero_si128());
   l[2] = _mm_max_epi16(l[2], _mm_setzero_si128());
   l[3] = _mm_max_epi16(l[3], _mm_setzero_si128());

   l[0] = _mm_min_epi16(l[0], _mm_set1_epi16(CLIP_10));
   l[1] = _mm_min_epi16(l[1], _mm_set1_epi16(CLIP_10));
   l[2] = _mm_min_epi16(l[2], _mm_set1_epi16(CLIP_10));
   l[3] = _mm_min_epi16(l[3], _mm_set1_epi16(CLIP_10));

   _mm_storeu_si128((__m128i*)&dst[0 * dst_stride], l[0]);
   _mm_storeu_si128((__m128i*)&dst[1 * dst_stride], l[1]);
   _mm_storeu_si128((__m128i*)&dst[8 + 0 * dst_stride], l[2]);
   _mm_storeu_si128((__m128i*)&dst[8 + 1 * dst_stride], l[3]);
}

static inline void
ovvc_transform_scale_sub_half_sse_32_1_10(uint16_t *dst, ptrdiff_t dst_stride,
                               const int16_t *src, ptrdiff_t src_stride, __m128i scale)
{
   __m128i l[4], r[4], s[4], m[8];
   __m128i a = _mm_unpackhi_epi16(scale, _mm_set1_epi16(1));
   __m128i b = _mm_set1_epi16(1 << (11 - 1));

   l[0] = _mm_loadu_si128((__m128i *)&dst[0]);
   l[1] = _mm_loadu_si128((__m128i *)&dst[8]);
   l[2] = _mm_loadu_si128((__m128i *)&dst[16]);
   l[3] = _mm_loadu_si128((__m128i *)&dst[24]);
   //         value = _src[j];
   r[0] = _mm_loadu_si128((__m128i *)&src[0]);
   r[1] = _mm_loadu_si128((__m128i *)&src[8]);
   r[2] = _mm_loadu_si128((__m128i *)&src[16]);
   r[3] = _mm_loadu_si128((__m128i *)&src[24]);

   //         sign  = value & (1 << 15);
   s[0] = _mm_srai_epi16(_mm_and_si128(r[0], _mm_set1_epi16(SIGN_16)), 15);
   s[1] = _mm_srai_epi16(_mm_and_si128(r[1], _mm_set1_epi16(SIGN_16)), 15);
   s[2] = _mm_srai_epi16(_mm_and_si128(r[2], _mm_set1_epi16(SIGN_16)), 15);
   s[3] = _mm_srai_epi16(_mm_and_si128(r[3], _mm_set1_epi16(SIGN_16)), 15);

   s[0] = _mm_andnot_si128(_mm_abs_epi16(s[0]), _mm_set1_epi16(1));
   s[1] = _mm_andnot_si128(_mm_abs_epi16(s[1]), _mm_set1_epi16(1));
   s[2] = _mm_andnot_si128(_mm_abs_epi16(s[2]), _mm_set1_epi16(1));
   s[3] = _mm_andnot_si128(_mm_abs_epi16(s[3]), _mm_set1_epi16(1));

   //         value = (-value);
   r[0] = _mm_xor_si128(r[0], _mm_set1_epi16((int16_t)0xFFFF));
   r[1] = _mm_xor_si128(r[1], _mm_set1_epi16((int16_t)0xFFFF));
   r[2] = _mm_xor_si128(r[2], _mm_set1_epi16((int16_t)0xFFFF));
   r[3] = _mm_xor_si128(r[3], _mm_set1_epi16((int16_t)0xFFFF));

   r[0] = _mm_add_epi16(r[0], _mm_set1_epi16(1));
   r[1] = _mm_add_epi16(r[1], _mm_set1_epi16(1));
   r[2] = _mm_add_epi16(r[2], _mm_set1_epi16(1));
   r[3] = _mm_add_epi16(r[3], _mm_set1_epi16(1));

   //         value = (abs(value);
   r[0] = _mm_srai_epi16(r[0], 1);
   r[1] = _mm_srai_epi16(r[1], 1);
   r[2] = _mm_srai_epi16(r[2], 1);
   r[3] = _mm_srai_epi16(r[3], 1);

   r[0] = _mm_abs_epi16(r[0]);
   r[1] = _mm_abs_epi16(r[1]);
   r[2] = _mm_abs_epi16(r[2]);
   r[3] = _mm_abs_epi16(r[3]);

   //         value = (value * scale + (1 << (11 - 1)));
   m[0] = _mm_unpacklo_epi16(r[0], b);
   m[1] = _mm_unpackhi_epi16(r[0], b);
   m[2] = _mm_unpacklo_epi16(r[1], b);
   m[3] = _mm_unpackhi_epi16(r[1], b);
   m[4] = _mm_unpacklo_epi16(r[2], b);
   m[5] = _mm_unpackhi_epi16(r[2], b);
   m[6] = _mm_unpacklo_epi16(r[3], b);
   m[7] = _mm_unpackhi_epi16(r[3], b);

   m[0] = _mm_madd_epi16(m[0], a);
   m[1] = _mm_madd_epi16(m[1], a);
   m[2] = _mm_madd_epi16(m[2], a);
   m[3] = _mm_madd_epi16(m[3], a);
   m[4] = _mm_madd_epi16(m[4], a);
   m[5] = _mm_madd_epi16(m[5], a);
   m[6] = _mm_madd_epi16(m[6], a);
   m[7] = _mm_madd_epi16(m[7], a);

   //         value = value >> 11;
   m[0] = _mm_srai_epi32(m[0], 11);
   m[1] = _mm_srai_epi32(m[1], 11);
   m[2] = _mm_srai_epi32(m[2], 11);
   m[3] = _mm_srai_epi32(m[3], 11);
   m[4] = _mm_srai_epi32(m[4], 11);
   m[5] = _mm_srai_epi32(m[5], 11);
   m[6] = _mm_srai_epi32(m[6], 11);
   m[7] = _mm_srai_epi32(m[7], 11);


   //         value = (ov_clip(value , 0,(1 << 16)-1));
   r[0] = _mm_packs_epi32(m[0], m[1]);
   r[1] = _mm_packs_epi32(m[2], m[3]);
   r[2] = _mm_packs_epi32(m[4], m[5]);
   r[3] = _mm_packs_epi32(m[6], m[7]);


   //         value = (sign ? -value : value);
   r[0] = _mm_xor_si128(r[0],_mm_cmpgt_epi16(s[0], _mm_setzero_si128()));
   r[1] = _mm_xor_si128(r[1],_mm_cmpgt_epi16(s[1], _mm_setzero_si128()));
   r[2] = _mm_xor_si128(r[2],_mm_cmpgt_epi16(s[2], _mm_setzero_si128()));
   r[3] = _mm_xor_si128(r[3],_mm_cmpgt_epi16(s[3], _mm_setzero_si128()));

   r[0] = _mm_add_epi16(r[0], s[0]);
   r[1] = _mm_add_epi16(r[1], s[1]);
   r[2] = _mm_add_epi16(r[2], s[2]);
   r[3] = _mm_add_epi16(r[3], s[3]);

   //         _dst[j] = ov_clip((int32_t)_dst[j] + value, 0, 1023);
   l[0] = _mm_add_epi16(l[0], r[0]);
   l[1] = _mm_add_epi16(l[1], r[1]);
   l[2] = _mm_add_epi16(l[2], r[2]);
   l[3] = _mm_add_epi16(l[3], r[3]);

   l[0] = _mm_max_epi16(l[0], _mm_setzero_si128());
   l[1] = _mm_max_epi16(l[1], _mm_setzero_si128());
   l[2] = _mm_max_epi16(l[2], _mm_setzero_si128());
   l[3] = _mm_max_epi16(l[3], _mm_setzero_si128());

   l[0] = _mm_min_epi16(l[0], _mm_set1_epi16(CLIP_10));
   l[1] = _mm_min_epi16(l[1], _mm_set1_epi16(CLIP_10));
   l[2] = _mm_min_epi16(l[2], _mm_set1_epi16(CLIP_10));
   l[3] = _mm_min_epi16(l[3], _mm_set1_epi16(CLIP_10));

   _mm_storeu_si128((__m128i*)&dst[0],  l[0]);
   _mm_storeu_si128((__m128i*)&dst[8],  l[1]);
   _mm_storeu_si128((__m128i*)&dst[16], l[2]);
   _mm_storeu_si128((__m128i*)&dst[24], l[3]);
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

static void
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

static void
vvc_add_residual_64_1_10_sse(const int16_t *const src, uint16_t *const dst,
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
        ovvc_transform_add_sse_32_1_10(_dst+32, RCN_CTB_STRIDE,
                                       _src+32, tb_w);
        _dst += RCN_CTB_STRIDE;
        _src += tb_w;
    }
}

static void
vvc_add_half_residual_8_4_10_sse(const int16_t *const src, uint16_t *const dst,
                            int log2_tb_w, int log2_tb_h, int scale)
{
    int i;
    int tb_w = 1 << log2_tb_w;
    int tb_h = 1 << log2_tb_h;
    const int16_t *_src = (const int16_t *)src;
    uint16_t *_dst = dst;
    if (log2_tb_h > 1) {
        for (i = 0; i < tb_h >> 2; ++i){
            ovvc_transform_add_half_sse_8_4_10(_dst, RCN_CTB_STRIDE,
                                          _src, tb_w);
            _dst += RCN_CTB_STRIDE << 2;
            _src += tb_w << 2;
        }
    } else {
        vvc_add_half_residual(src, dst, log2_tb_w, log2_tb_h, 0);
    }
}

static void
vvc_add_half_residual_16_2_10_sse(const int16_t *const src, uint16_t *const dst,
                            int log2_tb_w, int log2_tb_h, int scale)
{
    int i;
    int tb_w = 1 << log2_tb_w;
    int tb_h = 1 << log2_tb_h;
    const int16_t *_src = (const int16_t *)src;
    uint16_t *_dst = dst;
    if (log2_tb_h > 1) {
      for (i = 0; i < tb_h >> 1; ++i){
          ovvc_transform_add_half_sse_16_2_10(_dst, RCN_CTB_STRIDE,
                                         _src, tb_w);
          _dst += RCN_CTB_STRIDE << 1;
          _src += tb_w << 1;
      }
    } else {
        vvc_add_half_residual(src, dst, log2_tb_w, log2_tb_h, 0);
    }
}

static void
vvc_add_half_residual_32_1_10_sse(const int16_t *const src, uint16_t *const dst,
                     int log2_tb_w, int log2_tb_h, int scale)
{
    int i;
    int tb_w = 1 << log2_tb_w;
    int tb_h = 1 << log2_tb_h;
    const int16_t *_src = (const int16_t *)src;
    uint16_t *_dst = dst;
    for (i = 0; i < tb_h; ++i){
        ovvc_transform_add_half_sse_32_1_10(_dst, RCN_CTB_STRIDE,
                                       _src, tb_w);
        _dst += RCN_CTB_STRIDE;
        _src += tb_w;
    }
}

static void
vvc_sub_residual_8_4_10_sse(const int16_t *const src, uint16_t *const dst,
                            int log2_tb_w, int log2_tb_h, int scale)
{
    int i;
    int tb_w = 1 << log2_tb_w;
    int tb_h = 1 << log2_tb_h;
    const int16_t *_src = (const int16_t *)src;
    uint16_t *_dst = dst;
    if (log2_tb_h > 1) {
        for (i = 0; i < tb_h >> 2; ++i){
            ovvc_transform_sub_sse_8_4_10(_dst, RCN_CTB_STRIDE,
                                          _src, tb_w);
            _dst += RCN_CTB_STRIDE << 2;
            _src += tb_w << 2;
        }
    } else {
        vvc_sub_residual(src, dst, log2_tb_w, log2_tb_h, 0);
    }
}

static void
vvc_sub_residual_16_2_10_sse(const int16_t *const src, uint16_t *const dst,
                            int log2_tb_w, int log2_tb_h, int scale)
{
    int i;
    int tb_w = 1 << log2_tb_w;
    int tb_h = 1 << log2_tb_h;
    const int16_t *_src = (const int16_t *)src;
    uint16_t *_dst = dst;
    if (log2_tb_h > 1) {
      for (i = 0; i < tb_h >> 1; ++i){
          ovvc_transform_sub_sse_16_2_10(_dst, RCN_CTB_STRIDE,
                                         _src, tb_w);
          _dst += RCN_CTB_STRIDE << 1;
          _src += tb_w << 1;
      }
    } else {
        vvc_sub_residual(src, dst, log2_tb_w, log2_tb_h, 0);
    }
}

static void
vvc_sub_residual_32_1_10_sse(const int16_t *const src, uint16_t *const dst,
                     int log2_tb_w, int log2_tb_h, int scale)
{
    int i;
    int tb_w = 1 << log2_tb_w;
    int tb_h = 1 << log2_tb_h;
    const int16_t *_src = (const int16_t *)src;
    uint16_t *_dst = dst;
    for (i = 0; i < tb_h; ++i){
        ovvc_transform_sub_sse_32_1_10(_dst, RCN_CTB_STRIDE,
                                       _src, tb_w);
        _dst += RCN_CTB_STRIDE;
        _src += tb_w;
    }
}

static void
vvc_sub_half_residual_8_4_10_sse(const int16_t *const src, uint16_t *const dst,
                            int log2_tb_w, int log2_tb_h, int scale)
{
    int i;
    int tb_w = 1 << log2_tb_w;
    int tb_h = 1 << log2_tb_h;
    const int16_t *_src = (const int16_t *)src;
    uint16_t *_dst = dst;
    if (log2_tb_h > 1) {
        for (i = 0; i < tb_h >> 2; ++i){
            ovvc_transform_sub_half_sse_8_4_10(_dst, RCN_CTB_STRIDE,
                                          _src, tb_w);
            _dst += RCN_CTB_STRIDE << 2;
            _src += tb_w << 2;
        }
    } else {
        vvc_sub_half_residual(src, dst, log2_tb_w, log2_tb_h, 0);
    }
}

static void
vvc_sub_half_residual_16_2_10_sse(const int16_t *const src, uint16_t *const dst,
                            int log2_tb_w, int log2_tb_h, int scale)
{
    int i;
    int tb_w = 1 << log2_tb_w;
    int tb_h = 1 << log2_tb_h;
    const int16_t *_src = (const int16_t *)src;
    uint16_t *_dst = dst;
    if (log2_tb_h > 1) {
      for (i = 0; i < tb_h >> 1; ++i){
          ovvc_transform_sub_half_sse_16_2_10(_dst, RCN_CTB_STRIDE,
                                         _src, tb_w);
          _dst += RCN_CTB_STRIDE << 1;
          _src += tb_w << 1;
      }
    } else {
        vvc_sub_half_residual(src, dst, log2_tb_w, log2_tb_h, 0);
    }
}

static void
vvc_sub_half_residual_32_1_10_sse(const int16_t *const src, uint16_t *const dst,
                     int log2_tb_w, int log2_tb_h, int scale)
{
    int i;
    int tb_w = 1 << log2_tb_w;
    int tb_h = 1 << log2_tb_h;
    const int16_t *_src = (const int16_t *)src;
    uint16_t *_dst = dst;
    for (i = 0; i < tb_h; ++i){
        ovvc_transform_sub_half_sse_32_1_10(_dst, RCN_CTB_STRIDE,
                                       _src, tb_w);
        _dst += RCN_CTB_STRIDE;
        _src += tb_w;
    }
}

static void
vvc_scale_add_residual_8_4_10_sse(const int16_t *src, uint16_t *dst,
                      int log2_tb_w, int log2_tb_h,
                      int scale)
{

  int i;
  int tb_w = 1 << log2_tb_w;
  int tb_h = 1 << log2_tb_h;
  const int16_t *_src = (const int16_t *)src;
  uint16_t *_dst = dst;
  if (log2_tb_h > 1) {
    __m128i scale_vector = _mm_set1_epi16(scale);
    for (i = 0; i < tb_h >> 2; ++i){
        ovvc_transform_scale_add_sse_8_4_10(_dst, RCN_CTB_STRIDE,
                                       _src, tb_w, scale_vector);
        _dst += RCN_CTB_STRIDE << 2;
        _src += tb_w << 2;
    }
  } else {
      vvc_scale_add_residual(src, dst, log2_tb_w, log2_tb_h, scale);
  }
}

static void
vvc_scale_add_residual_16_2_10_sse(const int16_t *src, uint16_t *dst,
                      int log2_tb_w, int log2_tb_h,
                      int scale)
{
  int i;
  int tb_w = 1 << log2_tb_w;
  int tb_h = 1 << log2_tb_h;
  const int16_t *_src = (const int16_t *)src;
  uint16_t *_dst = dst;
  if (log2_tb_h > 1) {
    __m128i scale_vector = _mm_set1_epi16(scale);
    for (i = 0; i < tb_h >> 1; ++i){
        ovvc_transform_scale_add_sse_16_2_10(_dst, RCN_CTB_STRIDE,
                                       _src, tb_w, scale_vector);
        _dst += RCN_CTB_STRIDE << 1;
        _src += tb_w << 1;
    }
  } else {
      vvc_scale_add_residual(src, dst, log2_tb_w, log2_tb_h, scale);
  }
}

static void
vvc_scale_add_residual_32_1_10_sse(const int16_t *src, uint16_t *dst,
                      int log2_tb_w, int log2_tb_h,
                      int scale)
{
  int i;
  int tb_w = 1 << log2_tb_w;
  int tb_h = 1 << log2_tb_h;
  const int16_t *_src = (const int16_t *)src;
  uint16_t *_dst = dst;
  __m128i scale_vector = _mm_set1_epi16(scale);
  for (i = 0; i < tb_h; ++i){
      ovvc_transform_scale_add_sse_32_1_10(_dst, RCN_CTB_STRIDE,
                                     _src, tb_w, scale_vector);
      _dst += RCN_CTB_STRIDE;
      _src += tb_w;
  }
}

static void
vvc_scale_add_half_residual_8_4_10_sse(const int16_t *src, uint16_t *dst,
                      int log2_tb_w, int log2_tb_h,
                      int scale)
{
  int i;
  int tb_w = 1 << log2_tb_w;
  int tb_h = 1 << log2_tb_h;
  const int16_t *_src = (const int16_t *)src;
  uint16_t *_dst = dst;
  if (log2_tb_h > 1) {
    __m128i scale_vector = _mm_set1_epi16(scale);
    for (i = 0; i < tb_h >> 2; ++i){
        ovvc_transform_scale_add_half_sse_8_4_10(_dst, RCN_CTB_STRIDE,
                                       _src, tb_w, scale_vector);
        _dst += RCN_CTB_STRIDE << 2;
        _src += tb_w << 2;
    }
  } else {
      vvc_scale_add_half_residual(src, dst, log2_tb_w, log2_tb_h, scale);
  }
}

static void
vvc_scale_add_half_residual_16_2_10_sse(const int16_t *src, uint16_t *dst,
                      int log2_tb_w, int log2_tb_h,
                      int scale)
{
  int i;
  int tb_w = 1 << log2_tb_w;
  int tb_h = 1 << log2_tb_h;
  const int16_t *_src = (const int16_t *)src;
  uint16_t *_dst = dst;
  if (log2_tb_h > 1) {
    __m128i scale_vector = _mm_set1_epi16(scale);
    for (i = 0; i < tb_h >> 1; ++i){
        ovvc_transform_scale_add_half_sse_16_2_10(_dst, RCN_CTB_STRIDE,
                                       _src, tb_w, scale_vector);
        _dst += RCN_CTB_STRIDE << 1;
        _src += tb_w << 1;
    }
  } else {
      vvc_scale_add_half_residual(src, dst, log2_tb_w, log2_tb_h, scale);
  }
}

static void
vvc_scale_add_half_residual_32_1_10_sse(const int16_t *src, uint16_t *dst,
                      int log2_tb_w, int log2_tb_h,
                      int scale)
{
  int i;
  int tb_w = 1 << log2_tb_w;
  int tb_h = 1 << log2_tb_h;
  const int16_t *_src = (const int16_t *)src;
  uint16_t *_dst = dst;
  __m128i scale_vector = _mm_set1_epi16(scale);
  for (i = 0; i < tb_h; ++i){
      ovvc_transform_scale_add_half_sse_32_1_10(_dst, RCN_CTB_STRIDE,
                                     _src, tb_w, scale_vector);
      _dst += RCN_CTB_STRIDE;
      _src += tb_w;
  }
}

static void
vvc_scale_sub_residual_8_4_10_sse(const int16_t *src, uint16_t *dst,
                      int log2_tb_w, int log2_tb_h,
                      int scale)
{
  int i;
  int tb_w = 1 << log2_tb_w;
  int tb_h = 1 << log2_tb_h;
  const int16_t *_src = (const int16_t *)src;
  uint16_t *_dst = dst;
  if (log2_tb_h > 1) {
    __m128i scale_vector = _mm_set1_epi16(scale);
    for (i = 0; i < tb_h >> 2; ++i){
        ovvc_transform_scale_sub_sse_8_4_10(_dst, RCN_CTB_STRIDE,
                                       _src, tb_w, scale_vector);
        _dst += RCN_CTB_STRIDE << 2;
        _src += tb_w << 2;
    }
  } else {
      vvc_scale_sub_residual(src, dst, log2_tb_w, log2_tb_h, scale);
  }
}

static void
vvc_scale_sub_residual_16_2_10_sse(const int16_t *src, uint16_t *dst,
                      int log2_tb_w, int log2_tb_h,
                      int scale)
{
  int i;
  int tb_w = 1 << log2_tb_w;
  int tb_h = 1 << log2_tb_h;
  const int16_t *_src = (const int16_t *)src;
  uint16_t *_dst = dst;
  if (log2_tb_h > 1) {
    __m128i scale_vector = _mm_set1_epi16(scale);
    for (i = 0; i < tb_h >> 1; ++i){
        ovvc_transform_scale_sub_sse_16_2_10(_dst, RCN_CTB_STRIDE,
                                       _src, tb_w, scale_vector);
        _dst += RCN_CTB_STRIDE << 1;
        _src += tb_w << 1;
    }
  } else {
      vvc_scale_sub_residual(src, dst, log2_tb_w, log2_tb_h, scale);
  }
}

static void
vvc_scale_sub_residual_32_1_10_sse(const int16_t *src, uint16_t *dst,
                      int log2_tb_w, int log2_tb_h,
                      int scale)
{
  int i;
  int tb_w = 1 << log2_tb_w;
  int tb_h = 1 << log2_tb_h;
  const int16_t *_src = (const int16_t *)src;
  uint16_t *_dst = dst;
  __m128i scale_vector = _mm_set1_epi16(scale);
  for (i = 0; i < tb_h; ++i){
      ovvc_transform_scale_sub_sse_32_1_10(_dst, RCN_CTB_STRIDE,
                                     _src, tb_w, scale_vector);
      _dst += RCN_CTB_STRIDE;
      _src += tb_w;
  }
}

static void
vvc_scale_sub_half_residual_8_4_10_sse(const int16_t *src, uint16_t *dst,
                      int log2_tb_w, int log2_tb_h,
                      int scale)
{
  int i;
  int tb_w = 1 << log2_tb_w;
  int tb_h = 1 << log2_tb_h;
  const int16_t *_src = (const int16_t *)src;
  uint16_t *_dst = dst;
  if (log2_tb_h > 1) {
    __m128i scale_vector = _mm_set1_epi16(scale);
    for (i = 0; i < tb_h >> 2; ++i){
        ovvc_transform_scale_sub_half_sse_8_4_10(_dst, RCN_CTB_STRIDE,
                                       _src, tb_w, scale_vector);
        _dst += RCN_CTB_STRIDE << 2;
        _src += tb_w << 2;
    }
  } else {
      vvc_scale_sub_half_residual(src, dst, log2_tb_w, log2_tb_h, scale);
  }
}

static void
vvc_scale_sub_half_residual_16_2_10_sse(const int16_t *src, uint16_t *dst,
                      int log2_tb_w, int log2_tb_h,
                      int scale)
{
  int i;
  int tb_w = 1 << log2_tb_w;
  int tb_h = 1 << log2_tb_h;
  const int16_t *_src = (const int16_t *)src;
  uint16_t *_dst = dst;
  if (log2_tb_h > 1) {
    __m128i scale_vector = _mm_set1_epi16(scale);
    for (i = 0; i < tb_h >> 1; ++i){
        ovvc_transform_scale_sub_half_sse_16_2_10(_dst, RCN_CTB_STRIDE,
                                       _src, tb_w, scale_vector);
        _dst += RCN_CTB_STRIDE << 1;
        _src += tb_w << 1;
    }
  } else {
      vvc_scale_sub_half_residual(src, dst, log2_tb_w, log2_tb_h, scale);
  }
}

static void
vvc_scale_sub_half_residual_32_1_10_sse(const int16_t *src, uint16_t *dst,
                      int log2_tb_w, int log2_tb_h,
                      int scale)
{
  int i;
  int tb_w = 1 << log2_tb_w;
  int tb_h = 1 << log2_tb_h;
  const int16_t *_src = (const int16_t *)src;
  uint16_t *_dst = dst;
  __m128i scale_vector = _mm_set1_epi16(scale);
  for (i = 0; i < tb_h; ++i){
      ovvc_transform_scale_sub_half_sse_32_1_10(_dst, RCN_CTB_STRIDE,
                                     _src, tb_w, scale_vector);
      _dst += RCN_CTB_STRIDE;
      _src += tb_w;
  }
}

void
rcn_init_ict_functions_sse(struct RCNFunctions *rcn_func, uint8_t type)
{
    rcn_func->ict.add[3] = &vvc_add_residual_8_4_10_sse;
    rcn_func->ict.add[4] = &vvc_add_residual_16_2_10_sse;
    rcn_func->ict.add[5] = &vvc_add_residual_32_1_10_sse;
    rcn_func->ict.add[6] = &vvc_add_residual_64_1_10_sse;
    switch (type)
    {
        case 3:
            rcn_func->ict.ict[3][0] = &vvc_scale_add_residual_8_4_10_sse;
            rcn_func->ict.ict[4][0] = &vvc_scale_add_residual_16_2_10_sse;
            rcn_func->ict.ict[5][0] = &vvc_scale_add_residual_32_1_10_sse;

            rcn_func->ict.ict[3][1] = &vvc_scale_sub_residual_8_4_10_sse;
            rcn_func->ict.ict[4][1] = &vvc_scale_sub_residual_16_2_10_sse;
            rcn_func->ict.ict[5][1] = &vvc_scale_sub_residual_32_1_10_sse;

            rcn_func->ict.ict[3][2] = &vvc_scale_sub_half_residual_8_4_10_sse;
            rcn_func->ict.ict[4][2] = &vvc_scale_sub_half_residual_16_2_10_sse;
            rcn_func->ict.ict[5][2] = &vvc_scale_sub_half_residual_32_1_10_sse;
            break;
        case 2:
            rcn_func->ict.ict[3][0] = &vvc_add_residual_8_4_10_sse;
            rcn_func->ict.ict[4][0] = &vvc_add_residual_16_2_10_sse;
            rcn_func->ict.ict[5][0] = &vvc_add_residual_32_1_10_sse;

            rcn_func->ict.ict[3][1] = &vvc_sub_residual_8_4_10_sse;
            rcn_func->ict.ict[4][1] = &vvc_sub_residual_16_2_10_sse;
            rcn_func->ict.ict[5][1] = &vvc_sub_residual_32_1_10_sse;

            rcn_func->ict.ict[3][2] = &vvc_sub_half_residual_8_4_10_sse;
            rcn_func->ict.ict[4][2] = &vvc_sub_half_residual_16_2_10_sse;
            rcn_func->ict.ict[5][2] = &vvc_sub_half_residual_32_1_10_sse;
            break;
        case 1:
            rcn_func->ict.ict[3][0] = &vvc_scale_add_residual_8_4_10_sse;
            rcn_func->ict.ict[4][0] = &vvc_scale_add_residual_16_2_10_sse;
            rcn_func->ict.ict[5][0] = &vvc_scale_add_residual_32_1_10_sse;

            rcn_func->ict.ict[3][1] = &vvc_scale_add_residual_8_4_10_sse;
            rcn_func->ict.ict[4][1] = &vvc_scale_add_residual_16_2_10_sse;
            rcn_func->ict.ict[5][1] = &vvc_scale_add_residual_32_1_10_sse;

            rcn_func->ict.ict[3][2] = &vvc_scale_add_half_residual_8_4_10_sse;
            rcn_func->ict.ict[4][2] = &vvc_scale_add_half_residual_16_2_10_sse;
            rcn_func->ict.ict[5][2] = &vvc_scale_add_half_residual_32_1_10_sse;
            break;
        default:

            rcn_func->ict.ict[3][0] = &vvc_add_residual_8_4_10_sse;
            rcn_func->ict.ict[4][0] = &vvc_add_residual_16_2_10_sse;
            rcn_func->ict.ict[5][0] = &vvc_add_residual_32_1_10_sse;

            rcn_func->ict.ict[3][1] = &vvc_add_residual_8_4_10_sse;
            rcn_func->ict.ict[4][1] = &vvc_add_residual_16_2_10_sse;
            rcn_func->ict.ict[5][1] = &vvc_add_residual_32_1_10_sse;

            rcn_func->ict.ict[3][2] = &vvc_add_half_residual_8_4_10_sse;
            rcn_func->ict.ict[4][2] = &vvc_add_half_residual_16_2_10_sse;
            rcn_func->ict.ict[5][2] = &vvc_add_half_residual_32_1_10_sse;
            break;
    }
}
