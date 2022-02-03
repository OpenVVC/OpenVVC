#include <immintrin.h>
#include <stdint.h>
#include <stddef.h>

#include "ctudec.h"
#include "rcn_transform.h"
#include "rcn.h"
#include "ovutils.h"

#define CLIP_10 ((1 << 10) - 1)
#define SIGN_16 (int16_t)(1 << 15)

static inline void
ovvc_transform_add_avx2_8_4_10(uint16_t *dst, ptrdiff_t dst_stride,
                              const int16_t *src, ptrdiff_t src_stride)
{
    __m256i l0, l1, l2, l3, r0, r1, r2, r3;
    __m256i l01, l23, r01, r23;

    l0 = _mm256_loadu_si256((__m256i *)&dst[0 * dst_stride]);
    l1 = _mm256_loadu_si256((__m256i *)&dst[1 * dst_stride - 8]);
    l2 = _mm256_loadu_si256((__m256i *)&dst[2 * dst_stride]);
    l3 = _mm256_loadu_si256((__m256i *)&dst[3 * dst_stride - 8]);
    r0 = _mm256_loadu_si256((__m256i *)&src[0 * src_stride]);
    r1 = _mm256_loadu_si256((__m256i *)&src[1 * src_stride - 8]);
    r2 = _mm256_loadu_si256((__m256i *)&src[2 * src_stride]);
    r3 = _mm256_loadu_si256((__m256i *)&src[3 * src_stride - 8]);

    l01 = _mm256_blend_epi32(l0, l1, 0xF0);
    l23 = _mm256_blend_epi32(l2, l3, 0xF0);
    r01 = _mm256_blend_epi32(r0, r1, 0xF0);
    r23 = _mm256_blend_epi32(r2, r3, 0xF0);

    l01 = _mm256_add_epi16(l01, r01);
    l23 = _mm256_add_epi16(l23, r23);

    l01 = _mm256_max_epi16(l01, _mm256_setzero_si256());
    l23 = _mm256_max_epi16(l23, _mm256_setzero_si256());

    l01 = _mm256_min_epi16(l01, _mm256_set1_epi16(CLIP_10));
    l23 = _mm256_min_epi16(l23, _mm256_set1_epi16(CLIP_10));

    _mm_storeu_si128((__m128i*)&dst[0 * dst_stride], _mm256_extracti128_si256(l01, 0));
    _mm_storeu_si128((__m128i*)&dst[1 * dst_stride], _mm256_extracti128_si256(l01, 1));
    _mm_storeu_si128((__m128i*)&dst[2 * dst_stride], _mm256_extracti128_si256(l23, 0));
    _mm_storeu_si128((__m128i*)&dst[3 * dst_stride], _mm256_extracti128_si256(l23, 1));
}

static inline void
ovvc_transform_add_avx2_16_2_10(uint16_t *dst, ptrdiff_t dst_stride,
                               const int16_t *src, ptrdiff_t src_stride)
{
    __m256i l0, l1, r0, r1;

    l0 = _mm256_loadu_si256((__m256i *)&dst[0 * dst_stride]);
    l1 = _mm256_loadu_si256((__m256i *)&dst[1 * dst_stride]);
    r0 = _mm256_loadu_si256((__m256i *)&src[0 * src_stride]);
    r1 = _mm256_loadu_si256((__m256i *)&src[1 * src_stride]);

    l0 = _mm256_add_epi16(l0, r0);
    l1 = _mm256_add_epi16(l1, r1);

    l0 = _mm256_max_epi16(l0, _mm256_setzero_si256());
    l1 = _mm256_max_epi16(l1, _mm256_setzero_si256());

    l0 = _mm256_min_epi16(l0, _mm256_set1_epi16(CLIP_10));
    l1 = _mm256_min_epi16(l1, _mm256_set1_epi16(CLIP_10));

    _mm256_storeu_si256((__m256i*)&dst[0 * dst_stride], l0);
    _mm256_storeu_si256((__m256i*)&dst[1 * dst_stride], l1);
}

static inline void
ovvc_transform_add_avx2_32_1_10(uint16_t *dst, ptrdiff_t dst_stride,
                               const int16_t *src, ptrdiff_t src_stride)
{
    __m256i l0, l1, r0, r1;

    l0 = _mm256_loadu_si256((__m256i *)&dst[0]);
    l1 = _mm256_loadu_si256((__m256i *)&dst[16]);
    r0 = _mm256_loadu_si256((__m256i *)&src[0]);
    r1 = _mm256_loadu_si256((__m256i *)&src[16]);

    l0 = _mm256_add_epi16(l0, r0);
    l1 = _mm256_add_epi16(l1, r1);

    l0 = _mm256_max_epi16(l0, _mm256_setzero_si256());
    l1 = _mm256_max_epi16(l1, _mm256_setzero_si256());

    l0 = _mm256_min_epi16(l0, _mm256_set1_epi16(CLIP_10));
    l1 = _mm256_min_epi16(l1, _mm256_set1_epi16(CLIP_10));

    _mm256_storeu_si256((__m256i*)&dst[0], l0);
    _mm256_storeu_si256((__m256i*)&dst[16], l1);
}

static inline void
ovvc_transform_add_half_avx2_8_4_10(uint16_t *dst, ptrdiff_t dst_stride,
                              const int16_t *src, ptrdiff_t src_stride)
{
    __m256i l0, l1, l2, l3, r0, r1, r2, r3;
    __m256i l01, l23, r01, r23;
    
    l0 = _mm256_loadu_si256((__m256i *)&dst[0 * dst_stride]);
    l1 = _mm256_loadu_si256((__m256i *)&dst[1 * dst_stride - 8]);
    l2 = _mm256_loadu_si256((__m256i *)&dst[2 * dst_stride]);
    l3 = _mm256_loadu_si256((__m256i *)&dst[3 * dst_stride - 8]);
    r0 = _mm256_loadu_si256((__m256i *)&src[0 * src_stride]);
    r1 = _mm256_loadu_si256((__m256i *)&src[1 * src_stride - 8]);
    r2 = _mm256_loadu_si256((__m256i *)&src[2 * src_stride]);
    r3 = _mm256_loadu_si256((__m256i *)&src[3 * src_stride - 8]);

    l01 = _mm256_blend_epi32(l0, l1, 0xF0);
    l23 = _mm256_blend_epi32(l2, l3, 0xF0);
    r01 = _mm256_blend_epi32(r0, r1, 0xF0);
    r23 = _mm256_blend_epi32(r2, r3, 0xF0);

    r01 = _mm256_srai_epi16(r01, 1);
    r23 = _mm256_srai_epi16(r23, 1);

    l01 = _mm256_add_epi16(l01, r01);
    l23 = _mm256_add_epi16(l23, r23);

    l01 = _mm256_max_epi16(l01, _mm256_setzero_si256());
    l23 = _mm256_max_epi16(l23, _mm256_setzero_si256());

    l01 = _mm256_min_epi16(l01, _mm256_set1_epi16(CLIP_10));
    l23 = _mm256_min_epi16(l23, _mm256_set1_epi16(CLIP_10));

    _mm_storeu_si128((__m128i*)&dst[0 * dst_stride], _mm256_extracti128_si256(l01, 0));
    _mm_storeu_si128((__m128i*)&dst[1 * dst_stride], _mm256_extracti128_si256(l01, 1));
    _mm_storeu_si128((__m128i*)&dst[2 * dst_stride], _mm256_extracti128_si256(l23, 0));
    _mm_storeu_si128((__m128i*)&dst[3 * dst_stride], _mm256_extracti128_si256(l23, 1));
}

static inline void
ovvc_transform_add_half_avx2_16_2_10(uint16_t *dst, ptrdiff_t dst_stride,
                               const int16_t *src, ptrdiff_t src_stride)
{
    __m256i l0, l1, r0, r1;

    l0 = _mm256_loadu_si256((__m256i *)&dst[0 * dst_stride]);
    l1 = _mm256_loadu_si256((__m256i *)&dst[1 * dst_stride]);
    r0 = _mm256_loadu_si256((__m256i *)&src[0 * src_stride]);
    r1 = _mm256_loadu_si256((__m256i *)&src[1 * src_stride]);

    r0 = _mm256_srai_epi16(r0, 1);
    r1 = _mm256_srai_epi16(r1, 1);

    l0 = _mm256_add_epi16(l0, r0);
    l1 = _mm256_add_epi16(l1, r1);

    l0 = _mm256_max_epi16(l0, _mm256_setzero_si256());
    l1 = _mm256_max_epi16(l1, _mm256_setzero_si256());

    l0 = _mm256_min_epi16(l0, _mm256_set1_epi16(CLIP_10));
    l1 = _mm256_min_epi16(l1, _mm256_set1_epi16(CLIP_10));

    _mm256_storeu_si256((__m256i*)&dst[0 * dst_stride], l0);
    _mm256_storeu_si256((__m256i*)&dst[1 * dst_stride], l1);
}

static inline void
ovvc_transform_add_half_avx2_32_1_10(uint16_t *dst, ptrdiff_t dst_stride,
                               const int16_t *src, ptrdiff_t src_stride)
{
    __m256i l0, l1, r0, r1;

    l0 = _mm256_loadu_si256((__m256i *)&dst[0]);
    l1 = _mm256_loadu_si256((__m256i *)&dst[16]);
    r0 = _mm256_loadu_si256((__m256i *)&src[0]);
    r1 = _mm256_loadu_si256((__m256i *)&src[16]);

    r0 = _mm256_srai_epi16(r0, 1);
    r1 = _mm256_srai_epi16(r1, 1);

    l0 = _mm256_add_epi16(l0, r0);
    l1 = _mm256_add_epi16(l1, r1);

    l0 = _mm256_max_epi16(l0, _mm256_setzero_si256());
    l1 = _mm256_max_epi16(l1, _mm256_setzero_si256());

    l0 = _mm256_min_epi16(l0, _mm256_set1_epi16(CLIP_10));
    l1 = _mm256_min_epi16(l1, _mm256_set1_epi16(CLIP_10));

    _mm256_storeu_si256((__m256i*)&dst[0], l0);
    _mm256_storeu_si256((__m256i*)&dst[16], l1);
}

static inline void
ovvc_transform_sub_avx2_8_4_10(uint16_t *dst, ptrdiff_t dst_stride,
                              const int16_t *src, ptrdiff_t src_stride)
{
    __m256i l0, l1, l2, l3, r0, r1, r2, r3;
    __m256i l01, l23, r01, r23;

    l0 = _mm256_loadu_si256((__m256i *)&dst[0 * dst_stride]);
    l1 = _mm256_loadu_si256((__m256i *)&dst[1 * dst_stride - 8]);
    l2 = _mm256_loadu_si256((__m256i *)&dst[2 * dst_stride]);
    l3 = _mm256_loadu_si256((__m256i *)&dst[3 * dst_stride - 8]);
    r0 = _mm256_loadu_si256((__m256i *)&src[0 * src_stride]);
    r1 = _mm256_loadu_si256((__m256i *)&src[1 * src_stride - 8]);
    r2 = _mm256_loadu_si256((__m256i *)&src[2 * src_stride]);
    r3 = _mm256_loadu_si256((__m256i *)&src[3 * src_stride - 8]);

    l01 = _mm256_blend_epi32(l0, l1, 0xF0);
    l23 = _mm256_blend_epi32(l2, l3, 0xF0);
    r01 = _mm256_blend_epi32(r0, r1, 0xF0);
    r23 = _mm256_blend_epi32(r2, r3, 0xF0);

    l01 = _mm256_sub_epi16(l01, r01);
    l23 = _mm256_sub_epi16(l23, r23);

    l01 = _mm256_max_epi16(l01, _mm256_setzero_si256());
    l23 = _mm256_max_epi16(l23, _mm256_setzero_si256());

    l01 = _mm256_min_epi16(l01, _mm256_set1_epi16(CLIP_10));
    l23 = _mm256_min_epi16(l23, _mm256_set1_epi16(CLIP_10));

    _mm_storeu_si128((__m128i*)&dst[0 * dst_stride], _mm256_extracti128_si256(l01, 0));
    _mm_storeu_si128((__m128i*)&dst[1 * dst_stride], _mm256_extracti128_si256(l01, 1));
    _mm_storeu_si128((__m128i*)&dst[2 * dst_stride], _mm256_extracti128_si256(l23, 0));
    _mm_storeu_si128((__m128i*)&dst[3 * dst_stride], _mm256_extracti128_si256(l23, 1));
}

static inline void
ovvc_transform_sub_avx2_16_2_10(uint16_t *dst, ptrdiff_t dst_stride,
                               const int16_t *src, ptrdiff_t src_stride)
{
    __m256i l0, l1, r0, r1;

    l0 = _mm256_loadu_si256((__m256i *)&dst[0 * dst_stride]);
    l1 = _mm256_loadu_si256((__m256i *)&dst[1 * dst_stride]);
    r0 = _mm256_loadu_si256((__m256i *)&src[0 * src_stride]);
    r1 = _mm256_loadu_si256((__m256i *)&src[1 * src_stride]);

    l0 = _mm256_sub_epi16(l0, r0);
    l1 = _mm256_sub_epi16(l1, r1);

    l0 = _mm256_max_epi16(l0, _mm256_setzero_si256());
    l1 = _mm256_max_epi16(l1, _mm256_setzero_si256());

    l0 = _mm256_min_epi16(l0, _mm256_set1_epi16(CLIP_10));
    l1 = _mm256_min_epi16(l1, _mm256_set1_epi16(CLIP_10));

    _mm256_storeu_si256((__m256i*)&dst[0 * dst_stride], l0);
    _mm256_storeu_si256((__m256i*)&dst[1 * dst_stride], l1);
}

static inline void
ovvc_transform_sub_avx2_32_1_10(uint16_t *dst, ptrdiff_t dst_stride,
                               const int16_t *src, ptrdiff_t src_stride)
{
    __m256i l0, l1, r0, r1;

    l0 = _mm256_loadu_si256((__m256i *)&dst[0]);
    l1 = _mm256_loadu_si256((__m256i *)&dst[16]);
    r0 = _mm256_loadu_si256((__m256i *)&src[0]);
    r1 = _mm256_loadu_si256((__m256i *)&src[16]);

    l0 = _mm256_sub_epi16(l0, r0);
    l1 = _mm256_sub_epi16(l1, r1);

    l0 = _mm256_max_epi16(l0, _mm256_setzero_si256());
    l1 = _mm256_max_epi16(l1, _mm256_setzero_si256());

    l0 = _mm256_min_epi16(l0, _mm256_set1_epi16(CLIP_10));
    l1 = _mm256_min_epi16(l1, _mm256_set1_epi16(CLIP_10));

    _mm256_storeu_si256((__m256i*)&dst[0], l0);
    _mm256_storeu_si256((__m256i*)&dst[16], l1);
}

static inline void
ovvc_transform_sub_half_avx2_8_4_10(uint16_t *dst, ptrdiff_t dst_stride,
                              const int16_t *src, ptrdiff_t src_stride)
{
    __m256i l0, l1, l2, l3, r0, r1, r2, r3;
    __m256i l01, l23, r01, r23;

    l0 = _mm256_loadu_si256((__m256i *)&dst[0 * dst_stride]);
    l1 = _mm256_loadu_si256((__m256i *)&dst[1 * dst_stride - 8]);
    l2 = _mm256_loadu_si256((__m256i *)&dst[2 * dst_stride]);
    l3 = _mm256_loadu_si256((__m256i *)&dst[3 * dst_stride - 8]);
    r0 = _mm256_loadu_si256((__m256i *)&src[0 * src_stride]);
    r1 = _mm256_loadu_si256((__m256i *)&src[1 * src_stride - 8]);
    r2 = _mm256_loadu_si256((__m256i *)&src[2 * src_stride]);
    r3 = _mm256_loadu_si256((__m256i *)&src[3 * src_stride - 8]);

    l01 = _mm256_blend_epi32(l0, l1, 0xF0);
    l23 = _mm256_blend_epi32(l2, l3, 0xF0);
    r01 = _mm256_blend_epi32(r0, r1, 0xF0);
    r23 = _mm256_blend_epi32(r2, r3, 0xF0);

    //         value = (-value);
    r01 = _mm256_sub_epi16(_mm256_setzero_si256(), r01);
    r23 = _mm256_sub_epi16(_mm256_setzero_si256(), r23);

    //         value = (abs(value);
    r01 = _mm256_srai_epi16(r01, 1);
    r23 = _mm256_srai_epi16(r23, 1);

    l01 = _mm256_add_epi16(l01, r01);
    l23 = _mm256_add_epi16(l23, r23);

    l01 = _mm256_max_epi16(l01, _mm256_setzero_si256());
    l23 = _mm256_max_epi16(l23, _mm256_setzero_si256());

    l01 = _mm256_min_epi16(l01, _mm256_set1_epi16(CLIP_10));
    l23 = _mm256_min_epi16(l23, _mm256_set1_epi16(CLIP_10));

    _mm_storeu_si128((__m128i*)&dst[0 * dst_stride], _mm256_extracti128_si256(l01, 0));
    _mm_storeu_si128((__m128i*)&dst[1 * dst_stride], _mm256_extracti128_si256(l01, 1));
    _mm_storeu_si128((__m128i*)&dst[2 * dst_stride], _mm256_extracti128_si256(l23, 0));
    _mm_storeu_si128((__m128i*)&dst[3 * dst_stride], _mm256_extracti128_si256(l23, 1));
}

static inline void
ovvc_transform_sub_half_avx2_16_2_10(uint16_t *dst, ptrdiff_t dst_stride,
                               const int16_t *src, ptrdiff_t src_stride)
{
    __m256i l0, l1, r0, r1;

    l0 = _mm256_loadu_si256((__m256i *)&dst[0 * dst_stride]);
    l1 = _mm256_loadu_si256((__m256i *)&dst[1 * dst_stride]);
    r0 = _mm256_loadu_si256((__m256i *)&src[0 * src_stride]);
    r1 = _mm256_loadu_si256((__m256i *)&src[1 * src_stride]);

    //         value = (-value);
    r0 = _mm256_sub_epi16(_mm256_setzero_si256(), r0);
    r1 = _mm256_sub_epi16(_mm256_setzero_si256(), r1);

    //         value = (abs(value);
    r0 = _mm256_srai_epi16(r0, 1);
    r1 = _mm256_srai_epi16(r1, 1);

    l0 = _mm256_add_epi16(l0, r0);
    l1 = _mm256_add_epi16(l1, r1);

    l0 = _mm256_max_epi16(l0, _mm256_setzero_si256());
    l1 = _mm256_max_epi16(l1, _mm256_setzero_si256());

    l0 = _mm256_min_epi16(l0, _mm256_set1_epi16(CLIP_10));
    l1 = _mm256_min_epi16(l1, _mm256_set1_epi16(CLIP_10));

    _mm256_storeu_si256((__m256i*)&dst[0 * dst_stride], l0);
    _mm256_storeu_si256((__m256i*)&dst[1 * dst_stride], l1);
}

static inline void
ovvc_transform_sub_half_avx2_32_1_10(uint16_t *dst, ptrdiff_t dst_stride,
                               const int16_t *src, ptrdiff_t src_stride)
{
    __m256i l0, l1, r0, r1;

    l0 = _mm256_loadu_si256((__m256i *)&dst[0]);
    l1 = _mm256_loadu_si256((__m256i *)&dst[16]);
    r0 = _mm256_loadu_si256((__m256i *)&src[0]);
    r1 = _mm256_loadu_si256((__m256i *)&src[16]);

    //         value = (-value);
    r0 = _mm256_sub_epi16(_mm256_setzero_si256(), r0);
    r1 = _mm256_sub_epi16(_mm256_setzero_si256(), r1);

    //         value = (abs(value);
    r0 = _mm256_srai_epi16(r0, 1);
    r1 = _mm256_srai_epi16(r1, 1);

    l0 = _mm256_add_epi16(l0, r0);
    l1 = _mm256_add_epi16(l1, r1);

    l0 = _mm256_max_epi16(l0, _mm256_setzero_si256());
    l1 = _mm256_max_epi16(l1, _mm256_setzero_si256());

    l0 = _mm256_min_epi16(l0, _mm256_set1_epi16(CLIP_10));
    l1 = _mm256_min_epi16(l1, _mm256_set1_epi16(CLIP_10));

    _mm256_storeu_si256((__m256i*)&dst[0], l0);
    _mm256_storeu_si256((__m256i*)&dst[16], l1);
}

static inline void
ovvc_transform_scale_add_avx2_8_4_10(uint16_t *dst, ptrdiff_t dst_stride,
                               const int16_t *src, ptrdiff_t src_stride, __m256i scale)
{
    __m256i a = _mm256_unpackhi_epi16(scale, _mm256_set1_epi16(1));
    __m256i b = _mm256_set1_epi16(1 << (11 - 1));

    __m256i l0, l1, l2, l3, r0, r1, r2, r3;
    __m256i l01, l23, r01, r23, s01, s23;

    l0 = _mm256_loadu_si256((__m256i *)&dst[0 * dst_stride]);
    l1 = _mm256_loadu_si256((__m256i *)&dst[1 * dst_stride - 8]);
    l2 = _mm256_loadu_si256((__m256i *)&dst[2 * dst_stride]);
    l3 = _mm256_loadu_si256((__m256i *)&dst[3 * dst_stride - 8]);
    r0 = _mm256_loadu_si256((__m256i *)&src[0 * src_stride]);
    r1 = _mm256_loadu_si256((__m256i *)&src[1 * src_stride - 8]);
    r2 = _mm256_loadu_si256((__m256i *)&src[2 * src_stride]);
    r3 = _mm256_loadu_si256((__m256i *)&src[3 * src_stride - 8]);

    l01 = _mm256_blend_epi32(l0, l1, 0xF0);
    l23 = _mm256_blend_epi32(l2, l3, 0xF0);
    r01 = _mm256_blend_epi32(r0, r1, 0xF0);
    r23 = _mm256_blend_epi32(r2, r3, 0xF0);

    //         sign  = value & (1 << 15);
    s01 = _mm256_srai_epi16(_mm256_and_si256(r01, _mm256_set1_epi16(SIGN_16)), 15);
    s23 = _mm256_srai_epi16(_mm256_and_si256(r23, _mm256_set1_epi16(SIGN_16)), 15);


    s01 = _mm256_abs_epi16(s01);
    s23 = _mm256_abs_epi16(s23);

    //         value = (abs(value);
    r01 = _mm256_abs_epi16(r01);
    r23 = _mm256_abs_epi16(r23);

    //         value = (value * scale + (1 << (11 - 1)));
    __m256i m01l = _mm256_unpacklo_epi16(r01, b);
    __m256i m01h = _mm256_unpackhi_epi16(r01, b);
    __m256i m23l = _mm256_unpacklo_epi16(r23, b);
    __m256i m23h = _mm256_unpackhi_epi16(r23, b);


    m01l = _mm256_madd_epi16(m01l, a);
    m01h = _mm256_madd_epi16(m01h, a);
    m23l = _mm256_madd_epi16(m23l, a);
    m23h = _mm256_madd_epi16(m23h, a);

    //         value = value >> 11;
    m01l = _mm256_srai_epi32(m01l, 11);
    m01h = _mm256_srai_epi32(m01h, 11);
    m23l = _mm256_srai_epi32(m23l, 11);
    m23h = _mm256_srai_epi32(m23h, 11);


    //         value = (ov_clip(value , 0,(1 << 16)-1));
    r01 = _mm256_packs_epi32(m01l, m01h);
    r23 = _mm256_packs_epi32(m23l, m23h);

    //         value = (sign ? -value : value);
    r01 = _mm256_xor_si256(r01,_mm256_cmpgt_epi16(s01, _mm256_setzero_si256()));
    r23 = _mm256_xor_si256(r23,_mm256_cmpgt_epi16(s23, _mm256_setzero_si256()));

    r01 = _mm256_add_epi16(r01, s01);
    r23 = _mm256_add_epi16(r23, s23);

    //         _dst[j] = ov_clip((int32_t)_dst[j] + value, 0, 1023);
    l01 = _mm256_add_epi16(l01, r01);
    l23 = _mm256_add_epi16(l23, r23);

    l01 = _mm256_max_epi16(l01, _mm256_setzero_si256());
    l23 = _mm256_max_epi16(l23, _mm256_setzero_si256());

    l01 = _mm256_min_epi16(l01, _mm256_set1_epi16(CLIP_10));
    l23 = _mm256_min_epi16(l23, _mm256_set1_epi16(CLIP_10));

    _mm_storeu_si128((__m128i*)&dst[0 * dst_stride], _mm256_extracti128_si256(l01, 0));
    _mm_storeu_si128((__m128i*)&dst[1 * dst_stride], _mm256_extracti128_si256(l01, 1));
    _mm_storeu_si128((__m128i*)&dst[2 * dst_stride], _mm256_extracti128_si256(l23, 0));
    _mm_storeu_si128((__m128i*)&dst[3 * dst_stride], _mm256_extracti128_si256(l23, 1));
}

static inline void
ovvc_transform_scale_add_avx2_16_2_10(uint16_t *dst, ptrdiff_t dst_stride,
                               const int16_t *src, ptrdiff_t src_stride, __m256i scale)
{
    __m256i a = _mm256_unpackhi_epi16(scale, _mm256_set1_epi16(1));
    __m256i b = _mm256_set1_epi16(1 << (11 - 1));

    __m256i l0, l1, r0, r1, s0, s1, m0, m1, m2, m3;

    l0 = _mm256_loadu_si256((__m256i *)&dst[0 * dst_stride]);
    l1 = _mm256_loadu_si256((__m256i *)&dst[1 * dst_stride]);
    r0 = _mm256_loadu_si256((__m256i *)&src[0 * src_stride]);
    r1 = _mm256_loadu_si256((__m256i *)&src[1 * src_stride]);

    s0 = _mm256_srai_epi16(_mm256_and_si256(r0, _mm256_set1_epi16(SIGN_16)), 15);
    s1 = _mm256_srai_epi16(_mm256_and_si256(r1, _mm256_set1_epi16(SIGN_16)), 15);

    s0 = _mm256_abs_epi16(s0);
    s1 = _mm256_abs_epi16(s1);

    r0 = _mm256_abs_epi16(r0);
    r1 = _mm256_abs_epi16(r1);

    m0 = _mm256_unpacklo_epi16(r0, b);
    m1 = _mm256_unpackhi_epi16(r0, b);
    m2 = _mm256_unpacklo_epi16(r1, b);
    m3 = _mm256_unpackhi_epi16(r1, b);

    m0 = _mm256_madd_epi16(m0, a);
    m1 = _mm256_madd_epi16(m1, a);
    m2 = _mm256_madd_epi16(m2, a);
    m3 = _mm256_madd_epi16(m3, a);

    m0 = _mm256_srai_epi32(m0, 11);
    m1 = _mm256_srai_epi32(m1, 11);
    m2 = _mm256_srai_epi32(m2, 11);
    m3 = _mm256_srai_epi32(m3, 11);

    r0 = _mm256_packs_epi32(m0, m1);
    r1 = _mm256_packs_epi32(m2, m3);

    r0 = _mm256_xor_si256(r0,_mm256_cmpgt_epi16(s0, _mm256_setzero_si256()));
    r1 = _mm256_xor_si256(r1,_mm256_cmpgt_epi16(s1, _mm256_setzero_si256()));

    r0 = _mm256_add_epi16(r0, s0);
    r1 = _mm256_add_epi16(r1, s1);

    l0 = _mm256_add_epi16(l0, r0);
    l1 = _mm256_add_epi16(l1, r1);

    l0 = _mm256_max_epi16(l0, _mm256_setzero_si256());
    l1 = _mm256_max_epi16(l1, _mm256_setzero_si256());

    l0 = _mm256_min_epi16(l0, _mm256_set1_epi16(CLIP_10));
    l1 = _mm256_min_epi16(l1, _mm256_set1_epi16(CLIP_10));

    _mm256_storeu_si256((__m256i*)&dst[0 * dst_stride], l0);
    _mm256_storeu_si256((__m256i*)&dst[1 * dst_stride], l1);
}

static inline void
ovvc_transform_scale_add_avx2_32_1_10(uint16_t *dst, ptrdiff_t dst_stride,
                               const int16_t *src, ptrdiff_t src_stride, __m256i scale)
{
    __m256i a = _mm256_unpackhi_epi16(scale, _mm256_set1_epi16(1));
    __m256i b = _mm256_set1_epi16(1 << (11 - 1));

    __m256i l0, l1, r0, r1, s0, s1, m0, m1, m2, m3;

    l0 = _mm256_loadu_si256((__m256i *)&dst[0]);
    l1 = _mm256_loadu_si256((__m256i *)&dst[16]);
    r0 = _mm256_loadu_si256((__m256i *)&src[0]);
    r1 = _mm256_loadu_si256((__m256i *)&src[16]);

    s0 = _mm256_srai_epi16(_mm256_and_si256(r0, _mm256_set1_epi16(SIGN_16)), 15);
    s1 = _mm256_srai_epi16(_mm256_and_si256(r1, _mm256_set1_epi16(SIGN_16)), 15);

    s0 = _mm256_abs_epi16(s0);
    s1 = _mm256_abs_epi16(s1);

    r0 = _mm256_abs_epi16(r0);
    r1 = _mm256_abs_epi16(r1);

    m0 = _mm256_unpacklo_epi16(r0, b);
    m1 = _mm256_unpackhi_epi16(r0, b);
    m2 = _mm256_unpacklo_epi16(r1, b);
    m3 = _mm256_unpackhi_epi16(r1, b);

    m0 = _mm256_madd_epi16(m0, a);
    m1 = _mm256_madd_epi16(m1, a);
    m2 = _mm256_madd_epi16(m2, a);
    m3 = _mm256_madd_epi16(m3, a);

    m0 = _mm256_srai_epi32(m0, 11);
    m1 = _mm256_srai_epi32(m1, 11);
    m2 = _mm256_srai_epi32(m2, 11);
    m3 = _mm256_srai_epi32(m3, 11);

    r0 = _mm256_packs_epi32(m0, m1);
    r1 = _mm256_packs_epi32(m2, m3);

    r0 = _mm256_xor_si256(r0,_mm256_cmpgt_epi16(s0, _mm256_setzero_si256()));
    r1 = _mm256_xor_si256(r1,_mm256_cmpgt_epi16(s1, _mm256_setzero_si256()));

    r0 = _mm256_add_epi16(r0, s0);
    r1 = _mm256_add_epi16(r1, s1);

    l0 = _mm256_add_epi16(l0, r0);
    l1 = _mm256_add_epi16(l1, r1);

    l0 = _mm256_max_epi16(l0, _mm256_setzero_si256());
    l1 = _mm256_max_epi16(l1, _mm256_setzero_si256());

    l0 = _mm256_min_epi16(l0, _mm256_set1_epi16(CLIP_10));
    l1 = _mm256_min_epi16(l1, _mm256_set1_epi16(CLIP_10));

    _mm256_storeu_si256((__m256i*)&dst[0], l0);
    _mm256_storeu_si256((__m256i*)&dst[16], l1);
}

static inline void
ovvc_transform_scale_add_half_avx2_8_4_10(uint16_t *dst, ptrdiff_t dst_stride,
                               const int16_t *src, ptrdiff_t src_stride, __m256i scale)
{
    __m256i a = _mm256_unpackhi_epi16(scale, _mm256_set1_epi16(1));
   __m256i b = _mm256_set1_epi16(1 << (11 - 1));

    __m256i l0, l1, l2, l3, r0, r1, r2, r3;
    __m256i l01, l23, r01, r23, s01, s23;

    l0 = _mm256_loadu_si256((__m256i *)&dst[0 * dst_stride]);
    l1 = _mm256_loadu_si256((__m256i *)&dst[1 * dst_stride - 8]);
    l2 = _mm256_loadu_si256((__m256i *)&dst[2 * dst_stride]);
    l3 = _mm256_loadu_si256((__m256i *)&dst[3 * dst_stride - 8]);
    r0 = _mm256_loadu_si256((__m256i *)&src[0 * src_stride]);
    r1 = _mm256_loadu_si256((__m256i *)&src[1 * src_stride - 8]);
    r2 = _mm256_loadu_si256((__m256i *)&src[2 * src_stride]);
    r3 = _mm256_loadu_si256((__m256i *)&src[3 * src_stride - 8]);

    l01 = _mm256_blend_epi32(l0, l1, 0xF0);
    l23 = _mm256_blend_epi32(l2, l3, 0xF0);
    r01 = _mm256_blend_epi32(r0, r1, 0xF0);
    r23 = _mm256_blend_epi32(r2, r3, 0xF0);

   //         sign  = value & (1 << 15);
   s01 = _mm256_srai_epi16(_mm256_and_si256(r01, _mm256_set1_epi16(SIGN_16)), 15);
   s23 = _mm256_srai_epi16(_mm256_and_si256(r23, _mm256_set1_epi16(SIGN_16)), 15);


   s01 = _mm256_abs_epi16(s01);
   s23 = _mm256_abs_epi16(s23);

   //         value = (abs(value);
   r01 = _mm256_srai_epi16(r01, 1);
   r23 = _mm256_srai_epi16(r23, 1);

   r01 = _mm256_abs_epi16(r01);
   r23 = _mm256_abs_epi16(r23);

   //         value = (value * scale + (1 << (11 - 1)));
   __m256i m01l = _mm256_unpacklo_epi16(r01, b);
   __m256i m01h = _mm256_unpackhi_epi16(r01, b);
   __m256i m23l = _mm256_unpacklo_epi16(r23, b);
   __m256i m23h = _mm256_unpackhi_epi16(r23, b);


   m01l = _mm256_madd_epi16(m01l, a);
   m01h = _mm256_madd_epi16(m01h, a);
   m23l = _mm256_madd_epi16(m23l, a);
   m23h = _mm256_madd_epi16(m23h, a);

   //         value = value >> 11;
   m01l = _mm256_srai_epi32(m01l, 11);
   m01h = _mm256_srai_epi32(m01h, 11);
   m23l = _mm256_srai_epi32(m23l, 11);
   m23h = _mm256_srai_epi32(m23h, 11);


   //         value = (ov_clip(value , 0,(1 << 16)-1));
   r01 = _mm256_packs_epi32(m01l, m01h);
   r23 = _mm256_packs_epi32(m23l, m23h);

   //         value = (sign ? -value : value);
   r01 = _mm256_xor_si256(r01,_mm256_cmpgt_epi16(s01, _mm256_setzero_si256()));
   r23 = _mm256_xor_si256(r23,_mm256_cmpgt_epi16(s23, _mm256_setzero_si256()));

   r01 = _mm256_add_epi16(r01, s01);
   r23 = _mm256_add_epi16(r23, s23);

   //         _dst[j] = ov_clip((int32_t)_dst[j] + value, 0, 1023);
   l01 = _mm256_add_epi16(l01, r01);
   l23 = _mm256_add_epi16(l23, r23);

   l01 = _mm256_max_epi16(l01, _mm256_setzero_si256());
   l23 = _mm256_max_epi16(l23, _mm256_setzero_si256());

   l01 = _mm256_min_epi16(l01, _mm256_set1_epi16(CLIP_10));
   l23 = _mm256_min_epi16(l23, _mm256_set1_epi16(CLIP_10));

   _mm_storeu_si128((__m128i*)&dst[0 * dst_stride], _mm256_extracti128_si256(l01, 0));
   _mm_storeu_si128((__m128i*)&dst[1 * dst_stride], _mm256_extracti128_si256(l01, 1));
   _mm_storeu_si128((__m128i*)&dst[2 * dst_stride], _mm256_extracti128_si256(l23, 0));
   _mm_storeu_si128((__m128i*)&dst[3 * dst_stride], _mm256_extracti128_si256(l23, 1));
}

static inline void
ovvc_transform_scale_add_half_avx2_16_2_10(uint16_t *dst, ptrdiff_t dst_stride,
                               const int16_t *src, ptrdiff_t src_stride, __m256i scale)
{
    __m256i a = _mm256_unpackhi_epi16(scale, _mm256_set1_epi16(1));
    __m256i b = _mm256_set1_epi16(1 << (11 - 1));

    __m256i l0, l1, r0, r1, s0, s1, m0, m1, m2, m3;

    l0 = _mm256_loadu_si256((__m256i *)&dst[0 * dst_stride]);
    l1 = _mm256_loadu_si256((__m256i *)&dst[1 * dst_stride]);
    r0 = _mm256_loadu_si256((__m256i *)&src[0 * src_stride]);
    r1 = _mm256_loadu_si256((__m256i *)&src[1 * src_stride]);

    s0 = _mm256_srai_epi16(_mm256_and_si256(r0, _mm256_set1_epi16(SIGN_16)), 15);
    s1 = _mm256_srai_epi16(_mm256_and_si256(r1, _mm256_set1_epi16(SIGN_16)), 15);

    s0 = _mm256_abs_epi16(s0);
    s1 = _mm256_abs_epi16(s1);

    r0 = _mm256_srai_epi16(r0, 1);
    r1 = _mm256_srai_epi16(r1, 1);

    r0 = _mm256_abs_epi16(r0);
    r1 = _mm256_abs_epi16(r1);

    m0 = _mm256_unpacklo_epi16(r0, b);
    m1 = _mm256_unpackhi_epi16(r0, b);
    m2 = _mm256_unpacklo_epi16(r1, b);
    m3 = _mm256_unpackhi_epi16(r1, b);

    m0 = _mm256_madd_epi16(m0, a);
    m1 = _mm256_madd_epi16(m1, a);
    m2 = _mm256_madd_epi16(m2, a);
    m3 = _mm256_madd_epi16(m3, a);

    m0 = _mm256_srai_epi32(m0, 11);
    m1 = _mm256_srai_epi32(m1, 11);
    m2 = _mm256_srai_epi32(m2, 11);
    m3 = _mm256_srai_epi32(m3, 11);

    r0 = _mm256_packs_epi32(m0, m1);
    r1 = _mm256_packs_epi32(m2, m3);

    r0 = _mm256_xor_si256(r0,_mm256_cmpgt_epi16(s0, _mm256_setzero_si256()));
    r1 = _mm256_xor_si256(r1,_mm256_cmpgt_epi16(s1, _mm256_setzero_si256()));

    r0 = _mm256_add_epi16(r0, s0);
    r1 = _mm256_add_epi16(r1, s1);

    l0 = _mm256_add_epi16(l0, r0);
    l1 = _mm256_add_epi16(l1, r1);

    l0 = _mm256_max_epi16(l0, _mm256_setzero_si256());
    l1 = _mm256_max_epi16(l1, _mm256_setzero_si256());

    l0 = _mm256_min_epi16(l0, _mm256_set1_epi16(CLIP_10));
    l1 = _mm256_min_epi16(l1, _mm256_set1_epi16(CLIP_10));

    _mm256_storeu_si256((__m256i*)&dst[0 * dst_stride], l0);
    _mm256_storeu_si256((__m256i*)&dst[1 * dst_stride], l1);
}

static inline void
ovvc_transform_scale_add_half_avx2_32_1_10(uint16_t *dst, ptrdiff_t dst_stride,
                               const int16_t *src, ptrdiff_t src_stride, __m256i scale)
{
    __m256i a = _mm256_unpackhi_epi16(scale, _mm256_set1_epi16(1));
    __m256i b = _mm256_set1_epi16(1 << (11 - 1));

    __m256i l0, l1, r0, r1, s0, s1, m0, m1, m2, m3;

    l0 = _mm256_loadu_si256((__m256i *)&dst[0]);
    l1 = _mm256_loadu_si256((__m256i *)&dst[16]);
    r0 = _mm256_loadu_si256((__m256i *)&src[0]);
    r1 = _mm256_loadu_si256((__m256i *)&src[16]);

    s0 = _mm256_srai_epi16(_mm256_and_si256(r0, _mm256_set1_epi16(SIGN_16)), 15);
    s1 = _mm256_srai_epi16(_mm256_and_si256(r1, _mm256_set1_epi16(SIGN_16)), 15);

    s0 = _mm256_abs_epi16(s0);
    s1 = _mm256_abs_epi16(s1);

    r0 = _mm256_srai_epi16(r0, 1);
    r1 = _mm256_srai_epi16(r1, 1);

    r0 = _mm256_abs_epi16(r0);
    r1 = _mm256_abs_epi16(r1);

    m0 = _mm256_unpacklo_epi16(r0, b);
    m1 = _mm256_unpackhi_epi16(r0, b);
    m2 = _mm256_unpacklo_epi16(r1, b);
    m3 = _mm256_unpackhi_epi16(r1, b);

    m0 = _mm256_madd_epi16(m0, a);
    m1 = _mm256_madd_epi16(m1, a);
    m2 = _mm256_madd_epi16(m2, a);
    m3 = _mm256_madd_epi16(m3, a);

    m0 = _mm256_srai_epi32(m0, 11);
    m1 = _mm256_srai_epi32(m1, 11);
    m2 = _mm256_srai_epi32(m2, 11);
    m3 = _mm256_srai_epi32(m3, 11);

    r0 = _mm256_packs_epi32(m0, m1);
    r1 = _mm256_packs_epi32(m2, m3);

    r0 = _mm256_xor_si256(r0,_mm256_cmpgt_epi16(s0, _mm256_setzero_si256()));
    r1 = _mm256_xor_si256(r1,_mm256_cmpgt_epi16(s1, _mm256_setzero_si256()));

    r0 = _mm256_add_epi16(r0, s0);
    r1 = _mm256_add_epi16(r1, s1);

    l0 = _mm256_add_epi16(l0, r0);
    l1 = _mm256_add_epi16(l1, r1);

    l0 = _mm256_max_epi16(l0, _mm256_setzero_si256());
    l1 = _mm256_max_epi16(l1, _mm256_setzero_si256());

    l0 = _mm256_min_epi16(l0, _mm256_set1_epi16(CLIP_10));
    l1 = _mm256_min_epi16(l1, _mm256_set1_epi16(CLIP_10));

    _mm256_storeu_si256((__m256i*)&dst[0], l0);
    _mm256_storeu_si256((__m256i*)&dst[16], l1);
}

static inline void
ovvc_transform_scale_sub_avx2_8_4_10(uint16_t *dst, ptrdiff_t dst_stride,
                               const int16_t *src, ptrdiff_t src_stride, __m256i scale)
{
    __m256i a = _mm256_unpackhi_epi16(scale, _mm256_set1_epi16(1));
    __m256i b = _mm256_set1_epi16(1 << (11 - 1));

    __m256i l0, l1, l2, l3, r0, r1, r2, r3;
    __m256i l01, l23, r01, r23, s01, s23;

    l0 = _mm256_loadu_si256((__m256i *)&dst[0 * dst_stride]);
    l1 = _mm256_loadu_si256((__m256i *)&dst[1 * dst_stride - 8]);
    l2 = _mm256_loadu_si256((__m256i *)&dst[2 * dst_stride]);
    l3 = _mm256_loadu_si256((__m256i *)&dst[3 * dst_stride - 8]);
    r0 = _mm256_loadu_si256((__m256i *)&src[0 * src_stride]);
    r1 = _mm256_loadu_si256((__m256i *)&src[1 * src_stride - 8]);
    r2 = _mm256_loadu_si256((__m256i *)&src[2 * src_stride]);
    r3 = _mm256_loadu_si256((__m256i *)&src[3 * src_stride - 8]);

    l01 = _mm256_blend_epi32(l0, l1, 0xF0);
    l23 = _mm256_blend_epi32(l2, l3, 0xF0);
    r01 = _mm256_blend_epi32(r0, r1, 0xF0);
    r23 = _mm256_blend_epi32(r2, r3, 0xF0);

    //         sign  = value & (1 << 15);
    s01 = _mm256_srai_epi16(_mm256_and_si256(r01, _mm256_set1_epi16(SIGN_16)), 15);
    s23 = _mm256_srai_epi16(_mm256_and_si256(r23, _mm256_set1_epi16(SIGN_16)), 15);

    s01 = _mm256_andnot_si256(_mm256_abs_epi16(s01), _mm256_set1_epi16(1));
    s23 = _mm256_andnot_si256(_mm256_abs_epi16(s23), _mm256_set1_epi16(1));

    //         value = (abs(value);
    r01 = _mm256_abs_epi16(r01);
    r23 = _mm256_abs_epi16(r23);

    //         value = (value * scale + (1 << (11 - 1)));
    __m256i m01l = _mm256_unpacklo_epi16(r01, b);
    __m256i m01h = _mm256_unpackhi_epi16(r01, b);
    __m256i m23l = _mm256_unpacklo_epi16(r23, b);
    __m256i m23h = _mm256_unpackhi_epi16(r23, b);

    m01l = _mm256_madd_epi16(m01l, a);
    m01h = _mm256_madd_epi16(m01h, a);
    m23l = _mm256_madd_epi16(m23l, a);
    m23h = _mm256_madd_epi16(m23h, a);

    //         value = value >> 11;
    m01l = _mm256_srai_epi32(m01l, 11);
    m01h = _mm256_srai_epi32(m01h, 11);
    m23l = _mm256_srai_epi32(m23l, 11);
    m23h = _mm256_srai_epi32(m23h, 11);

    //         value = (ov_clip(value , 0,(1 << 16)-1));
    r01 = _mm256_packs_epi32(m01l, m01h);
    r23 = _mm256_packs_epi32(m23l, m23h);

    //         value = (sign ? value : -value);
    r01 = _mm256_xor_si256(r01,_mm256_cmpgt_epi16(s01, _mm256_setzero_si256()));
    r23 = _mm256_xor_si256(r23,_mm256_cmpgt_epi16(s23, _mm256_setzero_si256()));

    r01 = _mm256_add_epi16(r01, s01);
    r23 = _mm256_add_epi16(r23, s23);

    //         _dst[j] = ov_clip((int32_t)_dst[j] + value, 0, 1023);
    l01 = _mm256_add_epi16(l01, r01);
    l23 = _mm256_add_epi16(l23, r23);

    l01 = _mm256_max_epi16(l01, _mm256_setzero_si256());
    l23 = _mm256_max_epi16(l23, _mm256_setzero_si256());

    l01 = _mm256_min_epi16(l01, _mm256_set1_epi16(CLIP_10));
    l23 = _mm256_min_epi16(l23, _mm256_set1_epi16(CLIP_10));

    _mm_storeu_si128((__m128i*)&dst[0 * dst_stride], _mm256_extracti128_si256(l01, 0));
    _mm_storeu_si128((__m128i*)&dst[1 * dst_stride], _mm256_extracti128_si256(l01, 1));
    _mm_storeu_si128((__m128i*)&dst[2 * dst_stride], _mm256_extracti128_si256(l23, 0));
    _mm_storeu_si128((__m128i*)&dst[3 * dst_stride], _mm256_extracti128_si256(l23, 1));
}

static inline void
ovvc_transform_scale_sub_avx2_16_2_10(uint16_t *dst, ptrdiff_t dst_stride,
                               const int16_t *src, ptrdiff_t src_stride, __m256i scale)
{
    __m256i a = _mm256_unpackhi_epi16(scale, _mm256_set1_epi16(1));
    __m256i b = _mm256_set1_epi16(1 << (11 - 1));

    __m256i l0, l1, r0, r1, s0, s1, m0, m1, m2, m3;

    l0 = _mm256_loadu_si256((__m256i *)&dst[0 * dst_stride]);
    l1 = _mm256_loadu_si256((__m256i *)&dst[1 * dst_stride]);
    r0 = _mm256_loadu_si256((__m256i *)&src[0 * src_stride]);
    r1 = _mm256_loadu_si256((__m256i *)&src[1 * src_stride]);

    //         sign  = value & (1 << 15);
    s0 = _mm256_srai_epi16(_mm256_and_si256(r0, _mm256_set1_epi16(SIGN_16)), 15);
    s1 = _mm256_srai_epi16(_mm256_and_si256(r1, _mm256_set1_epi16(SIGN_16)), 15);

    s0 = _mm256_andnot_si256(_mm256_abs_epi16(s0), _mm256_set1_epi16(1));
    s1 = _mm256_andnot_si256(_mm256_abs_epi16(s1), _mm256_set1_epi16(1));

    //         value = (abs(value);
    r0 = _mm256_abs_epi16(r0);
    r1 = _mm256_abs_epi16(r1);

    //         value = (value * scale + (1 << (11 - 1)));
    m0 = _mm256_unpacklo_epi16(r0, b);
    m1 = _mm256_unpackhi_epi16(r0, b);
    m2 = _mm256_unpacklo_epi16(r1, b);
    m3 = _mm256_unpackhi_epi16(r1, b);

    m0 = _mm256_madd_epi16(m0, a);
    m1 = _mm256_madd_epi16(m1, a);
    m2 = _mm256_madd_epi16(m2, a);
    m3 = _mm256_madd_epi16(m3, a);

    //         value = value >> 11;
    m0 = _mm256_srai_epi32(m0, 11);
    m1 = _mm256_srai_epi32(m1, 11);
    m2 = _mm256_srai_epi32(m2, 11);
    m3 = _mm256_srai_epi32(m3, 11);

    //         value = (ov_clip(value , 0,(1 << 16)-1));
    r0 = _mm256_packs_epi32(m0, m1);
    r1 = _mm256_packs_epi32(m2, m3);

    //         value = (sign ? -value : value);
    r0 = _mm256_xor_si256(r0,_mm256_cmpgt_epi16(s0, _mm256_setzero_si256()));
    r1 = _mm256_xor_si256(r1,_mm256_cmpgt_epi16(s1, _mm256_setzero_si256()));

    r0 = _mm256_add_epi16(r0, s0);
    r1 = _mm256_add_epi16(r1, s1);

    //         _dst[j] = ov_clip((int32_t)_dst[j] + value, 0, 1023);
    l0 = _mm256_add_epi16(l0, r0);
    l1 = _mm256_add_epi16(l1, r1);

    l0 = _mm256_max_epi16(l0, _mm256_setzero_si256());
    l1 = _mm256_max_epi16(l1, _mm256_setzero_si256());

    l0 = _mm256_min_epi16(l0, _mm256_set1_epi16(CLIP_10));
    l1 = _mm256_min_epi16(l1, _mm256_set1_epi16(CLIP_10));

    _mm256_storeu_si256((__m256i*)&dst[0 * dst_stride], l0);
    _mm256_storeu_si256((__m256i*)&dst[1 * dst_stride], l1);
}

static inline void
ovvc_transform_scale_sub_avx2_32_1_10(uint16_t *dst, ptrdiff_t dst_stride,
                               const int16_t *src, ptrdiff_t src_stride, __m256i scale)
{
    __m256i a = _mm256_unpackhi_epi16(scale, _mm256_set1_epi16(1));
    __m256i b = _mm256_set1_epi16(1 << (11 - 1));

    __m256i l0, l1, r0, r1, s0, s1, m0, m1, m2, m3;

    l0 = _mm256_loadu_si256((__m256i *)&dst[0]);
    l1 = _mm256_loadu_si256((__m256i *)&dst[16]);
    r0 = _mm256_loadu_si256((__m256i *)&src[0]);
    r1 = _mm256_loadu_si256((__m256i *)&src[16]);

    //         sign  = value & (1 << 15);
    s0 = _mm256_srai_epi16(_mm256_and_si256(r0, _mm256_set1_epi16(SIGN_16)), 15);
    s1 = _mm256_srai_epi16(_mm256_and_si256(r1, _mm256_set1_epi16(SIGN_16)), 15);

    s0 = _mm256_andnot_si256(_mm256_abs_epi16(s0), _mm256_set1_epi16(1));
    s1 = _mm256_andnot_si256(_mm256_abs_epi16(s1), _mm256_set1_epi16(1));

    //         value = (abs(value);
    r0 = _mm256_abs_epi16(r0);
    r1 = _mm256_abs_epi16(r1);

    //         value = (value * scale + (1 << (11 - 1)));
    m0 = _mm256_unpacklo_epi16(r0, b);
    m1 = _mm256_unpackhi_epi16(r0, b);
    m2 = _mm256_unpacklo_epi16(r1, b);
    m3 = _mm256_unpackhi_epi16(r1, b);

    m0 = _mm256_madd_epi16(m0, a);
    m1 = _mm256_madd_epi16(m1, a);
    m2 = _mm256_madd_epi16(m2, a);
    m3 = _mm256_madd_epi16(m3, a);

    //         value = value >> 11;
    m0 = _mm256_srai_epi32(m0, 11);
    m1 = _mm256_srai_epi32(m1, 11);
    m2 = _mm256_srai_epi32(m2, 11);
    m3 = _mm256_srai_epi32(m3, 11);

    //         value = (ov_clip(value , 0,(1 << 16)-1));
    r0 = _mm256_packs_epi32(m0, m1);
    r1 = _mm256_packs_epi32(m2, m3);

    //         value = (sign ? -value : value);
    r0 = _mm256_xor_si256(r0,_mm256_cmpgt_epi16(s0, _mm256_setzero_si256()));
    r1 = _mm256_xor_si256(r1,_mm256_cmpgt_epi16(s1, _mm256_setzero_si256()));

    r0 = _mm256_add_epi16(r0, s0);
    r1 = _mm256_add_epi16(r1, s1);

    //         _dst[j] = ov_clip((int32_t)_dst[j] + value, 0, 1023);
    l0 = _mm256_add_epi16(l0, r0);
    l1 = _mm256_add_epi16(l1, r1);

    l0 = _mm256_max_epi16(l0, _mm256_setzero_si256());
    l1 = _mm256_max_epi16(l1, _mm256_setzero_si256());

    l0 = _mm256_min_epi16(l0, _mm256_set1_epi16(CLIP_10));
    l1 = _mm256_min_epi16(l1, _mm256_set1_epi16(CLIP_10));

    _mm256_storeu_si256((__m256i*)&dst[0], l0);
    _mm256_storeu_si256((__m256i*)&dst[16], l1);
}

static inline void
ovvc_transform_scale_sub_half_avx2_8_4_10(uint16_t *dst, ptrdiff_t dst_stride,
                               const int16_t *src, ptrdiff_t src_stride, __m256i scale)
{
    __m256i a = _mm256_unpackhi_epi16(scale, _mm256_set1_epi16(1));
    __m256i b = _mm256_set1_epi16(1 << (11 - 1));

    __m256i l0, l1, l2, l3, r0, r1, r2, r3;
    __m256i l01, l23, r01, r23, s01, s23;

    l0 = _mm256_loadu_si256((__m256i *)&dst[0 * dst_stride]);
    l1 = _mm256_loadu_si256((__m256i *)&dst[1 * dst_stride - 8]);
    l2 = _mm256_loadu_si256((__m256i *)&dst[2 * dst_stride]);
    l3 = _mm256_loadu_si256((__m256i *)&dst[3 * dst_stride - 8]);
    r0 = _mm256_loadu_si256((__m256i *)&src[0 * src_stride]);
    r1 = _mm256_loadu_si256((__m256i *)&src[1 * src_stride - 8]);
    r2 = _mm256_loadu_si256((__m256i *)&src[2 * src_stride]);
    r3 = _mm256_loadu_si256((__m256i *)&src[3 * src_stride - 8]);

    l01 = _mm256_blend_epi32(l0, l1, 0xF0);
    l23 = _mm256_blend_epi32(l2, l3, 0xF0);
    r01 = _mm256_blend_epi32(r0, r1, 0xF0);
    r23 = _mm256_blend_epi32(r2, r3, 0xF0);


    //         sign  = value & (1 << 15);
    s01 = _mm256_srai_epi16(_mm256_and_si256(r01, _mm256_set1_epi16(SIGN_16)), 15);
    s23 = _mm256_srai_epi16(_mm256_and_si256(r23, _mm256_set1_epi16(SIGN_16)), 15);

    s01 = _mm256_andnot_si256(_mm256_abs_epi16(s01), _mm256_set1_epi16(1));
    s23 = _mm256_andnot_si256(_mm256_abs_epi16(s23), _mm256_set1_epi16(1));

    //         value = (-value);
    r01 = _mm256_xor_si256(r01, _mm256_set1_epi16((int16_t)0xFFFF));
    r23 = _mm256_xor_si256(r23, _mm256_set1_epi16((int16_t)0xFFFF));

    r01 = _mm256_add_epi16(r01, _mm256_set1_epi16(1));
    r23 = _mm256_add_epi16(r23, _mm256_set1_epi16(1));

    //         value = (abs(value);
    r01 = _mm256_srai_epi16(r01, 1);
    r23 = _mm256_srai_epi16(r23, 1);

    r01 = _mm256_abs_epi16(r01);
    r23 = _mm256_abs_epi16(r23);

    //         value = (value * scale + (1 << (11 - 1)));
    __m256i m01l = _mm256_unpacklo_epi16(r01, b);
    __m256i m01h = _mm256_unpackhi_epi16(r01, b);
    __m256i m23l = _mm256_unpacklo_epi16(r23, b);
    __m256i m23h = _mm256_unpackhi_epi16(r23, b);

    m01l = _mm256_madd_epi16(m01l, a);
    m01h = _mm256_madd_epi16(m01h, a);
    m23l = _mm256_madd_epi16(m23l, a);
    m23h = _mm256_madd_epi16(m23h, a);

    //         value = value >> 11;
    m01l = _mm256_srai_epi32(m01l, 11);
    m01h = _mm256_srai_epi32(m01h, 11);
    m23l = _mm256_srai_epi32(m23l, 11);
    m23h = _mm256_srai_epi32(m23h, 11);

    //         value = (ov_clip(value , 0,(1 << 16)-1));
    r01 = _mm256_packs_epi32(m01l, m01h);
    r23 = _mm256_packs_epi32(m23l, m23h);

    //         value = (sign ? -value : value);
    r01 = _mm256_xor_si256(r01,_mm256_cmpgt_epi16(s01, _mm256_setzero_si256()));
    r23 = _mm256_xor_si256(r23,_mm256_cmpgt_epi16(s23, _mm256_setzero_si256()));

    r01 = _mm256_add_epi16(r01, s01);
    r23 = _mm256_add_epi16(r23, s23);

    //         _dst[j] = ov_clip((int32_t)_dst[j] + value, 0, 1023);
    l01 = _mm256_add_epi16(l01, r01);
    l23 = _mm256_add_epi16(l23, r23);

    l01 = _mm256_max_epi16(l01, _mm256_setzero_si256());
    l23 = _mm256_max_epi16(l23, _mm256_setzero_si256());

    l01 = _mm256_min_epi16(l01, _mm256_set1_epi16(CLIP_10));
    l23 = _mm256_min_epi16(l23, _mm256_set1_epi16(CLIP_10));

    _mm_storeu_si128((__m128i*)&dst[0 * dst_stride], _mm256_extracti128_si256(l01, 0));
    _mm_storeu_si128((__m128i*)&dst[1 * dst_stride], _mm256_extracti128_si256(l01, 1));
    _mm_storeu_si128((__m128i*)&dst[2 * dst_stride], _mm256_extracti128_si256(l23, 0));
    _mm_storeu_si128((__m128i*)&dst[3 * dst_stride], _mm256_extracti128_si256(l23, 1));
}

static inline void
ovvc_transform_scale_sub_half_avx2_16_2_10(uint16_t *dst, ptrdiff_t dst_stride,
                               const int16_t *src, ptrdiff_t src_stride, __m256i scale)
{
    __m256i a = _mm256_unpackhi_epi16(scale, _mm256_set1_epi16(1));
    __m256i b = _mm256_set1_epi16(1 << (11 - 1));

    __m256i l0, l1, r0, r1, s0, s1, m0, m1, m2, m3;

    l0 = _mm256_loadu_si256((__m256i *)&dst[0 * dst_stride]);
    l1 = _mm256_loadu_si256((__m256i *)&dst[1 * dst_stride]);
    r0 = _mm256_loadu_si256((__m256i *)&src[0 * src_stride]);
    r1 = _mm256_loadu_si256((__m256i *)&src[1 * src_stride]);

    //         sign  = value & (1 << 15);
    s0 = _mm256_srai_epi16(_mm256_and_si256(r0, _mm256_set1_epi16(SIGN_16)), 15);
    s1 = _mm256_srai_epi16(_mm256_and_si256(r1, _mm256_set1_epi16(SIGN_16)), 15);

    s0 = _mm256_andnot_si256(_mm256_abs_epi16(s0), _mm256_set1_epi16(1));
    s1 = _mm256_andnot_si256(_mm256_abs_epi16(s1), _mm256_set1_epi16(1));

    //         value = (-value);
    r0 = _mm256_xor_si256(r0, _mm256_set1_epi16((int16_t)0xFFFF));
    r1 = _mm256_xor_si256(r1, _mm256_set1_epi16((int16_t)0xFFFF));

    r0 = _mm256_add_epi16(r0, _mm256_set1_epi16(1));
    r1 = _mm256_add_epi16(r1, _mm256_set1_epi16(1));

    //         value = (abs(value);
    r0 = _mm256_srai_epi16(r0, 1);
    r1 = _mm256_srai_epi16(r1, 1);

    r0 = _mm256_abs_epi16(r0);
    r1 = _mm256_abs_epi16(r1);

    //         value = (value * scale + (1 << (11 - 1)));
    m0 = _mm256_unpacklo_epi16(r0, b);
    m1 = _mm256_unpackhi_epi16(r0, b);
    m2 = _mm256_unpacklo_epi16(r1, b);
    m3 = _mm256_unpackhi_epi16(r1, b);

    m0 = _mm256_madd_epi16(m0, a);
    m1 = _mm256_madd_epi16(m1, a);
    m2 = _mm256_madd_epi16(m2, a);
    m3 = _mm256_madd_epi16(m3, a);

    //         value = value >> 11;
    m0 = _mm256_srai_epi32(m0, 11);
    m1 = _mm256_srai_epi32(m1, 11);
    m2 = _mm256_srai_epi32(m2, 11);
    m3 = _mm256_srai_epi32(m3, 11);

    //         value = (ov_clip(value , 0,(1 << 16)-1));
    r0 = _mm256_packs_epi32(m0, m1);
    r1 = _mm256_packs_epi32(m2, m3);

    //         value = (sign ? -value : value);
    r0 = _mm256_xor_si256(r0,_mm256_cmpgt_epi16(s0, _mm256_setzero_si256()));
    r1 = _mm256_xor_si256(r1,_mm256_cmpgt_epi16(s1, _mm256_setzero_si256()));

    r0 = _mm256_add_epi16(r0, s0);
    r1 = _mm256_add_epi16(r1, s1);

    //         _dst[j] = ov_clip((int32_t)_dst[j] + value, 0, 1023);
    l0 = _mm256_add_epi16(l0, r0);
    l1 = _mm256_add_epi16(l1, r1);

    l0 = _mm256_max_epi16(l0, _mm256_setzero_si256());
    l1 = _mm256_max_epi16(l1, _mm256_setzero_si256());

    l0 = _mm256_min_epi16(l0, _mm256_set1_epi16(CLIP_10));
    l1 = _mm256_min_epi16(l1, _mm256_set1_epi16(CLIP_10));

    _mm256_storeu_si256((__m256i*)&dst[0 * dst_stride], l0);
    _mm256_storeu_si256((__m256i*)&dst[1 * dst_stride], l1);
}

static inline void
ovvc_transform_scale_sub_half_avx2_32_1_10(uint16_t *dst, ptrdiff_t dst_stride,
                               const int16_t *src, ptrdiff_t src_stride, __m256i scale)
{
    __m256i a = _mm256_unpackhi_epi16(scale, _mm256_set1_epi16(1));
    __m256i b = _mm256_set1_epi16(1 << (11 - 1));

    __m256i l0, l1, r0, r1, s0, s1, m0, m1, m2, m3;

    l0 = _mm256_loadu_si256((__m256i *)&dst[0]);
    l1 = _mm256_loadu_si256((__m256i *)&dst[16]);
    r0 = _mm256_loadu_si256((__m256i *)&src[0]);
    r1 = _mm256_loadu_si256((__m256i *)&src[16]);

    //         sign  = value & (1 << 15);
    s0 = _mm256_srai_epi16(_mm256_and_si256(r0, _mm256_set1_epi16(SIGN_16)), 15);
    s1 = _mm256_srai_epi16(_mm256_and_si256(r1, _mm256_set1_epi16(SIGN_16)), 15);

    s0 = _mm256_andnot_si256(_mm256_abs_epi16(s0), _mm256_set1_epi16(1));
    s1 = _mm256_andnot_si256(_mm256_abs_epi16(s1), _mm256_set1_epi16(1));

    //         value = (-value);
    r0 = _mm256_xor_si256(r0, _mm256_set1_epi16((int16_t)0xFFFF));
    r1 = _mm256_xor_si256(r1, _mm256_set1_epi16((int16_t)0xFFFF));

    r0 = _mm256_add_epi16(r0, _mm256_set1_epi16(1));
    r1 = _mm256_add_epi16(r1, _mm256_set1_epi16(1));

    //         value = (abs(value);
    r0 = _mm256_srai_epi16(r0, 1);
    r1 = _mm256_srai_epi16(r1, 1);

    r0 = _mm256_abs_epi16(r0);
    r1 = _mm256_abs_epi16(r1);

    //         value = (value * scale + (1 << (11 - 1)));
    m0 = _mm256_unpacklo_epi16(r0, b);
    m1 = _mm256_unpackhi_epi16(r0, b);
    m2 = _mm256_unpacklo_epi16(r1, b);
    m3 = _mm256_unpackhi_epi16(r1, b);

    m0 = _mm256_madd_epi16(m0, a);
    m1 = _mm256_madd_epi16(m1, a);
    m2 = _mm256_madd_epi16(m2, a);
    m3 = _mm256_madd_epi16(m3, a);

    //         value = value >> 11;
    m0 = _mm256_srai_epi32(m0, 11);
    m1 = _mm256_srai_epi32(m1, 11);
    m2 = _mm256_srai_epi32(m2, 11);
    m3 = _mm256_srai_epi32(m3, 11);

    //         value = (ov_clip(value , 0,(1 << 16)-1));
    r0 = _mm256_packs_epi32(m0, m1);
    r1 = _mm256_packs_epi32(m2, m3);

    //         value = (sign ? -value : value);
    r0 = _mm256_xor_si256(r0,_mm256_cmpgt_epi16(s0, _mm256_setzero_si256()));
    r1 = _mm256_xor_si256(r1,_mm256_cmpgt_epi16(s1, _mm256_setzero_si256()));

    r0 = _mm256_add_epi16(r0, s0);
    r1 = _mm256_add_epi16(r1, s1);

    //         _dst[j] = ov_clip((int32_t)_dst[j] + value, 0, 1023);
    l0 = _mm256_add_epi16(l0, r0);
    l1 = _mm256_add_epi16(l1, r1);

    l0 = _mm256_max_epi16(l0, _mm256_setzero_si256());
    l1 = _mm256_max_epi16(l1, _mm256_setzero_si256());

    l0 = _mm256_min_epi16(l0, _mm256_set1_epi16(CLIP_10));
    l1 = _mm256_min_epi16(l1, _mm256_set1_epi16(CLIP_10));

    _mm256_storeu_si256((__m256i*)&dst[0], l0);
    _mm256_storeu_si256((__m256i*)&dst[16], l1);
}

static void
vvc_add_residual_8_4_10_avx2(const int16_t *const src, uint16_t *const dst,
                            int log2_tb_w, int log2_tb_h, int scale)
{
    int i;
    int tb_w = 1 << log2_tb_w;
    int tb_h = 1 << log2_tb_h;
    const int16_t *_src = (const int16_t *)src;
    uint16_t *_dst = dst;
    if (log2_tb_h > 1) {
        for (i = 0; i < tb_h >> 2; ++i){
            ovvc_transform_add_avx2_8_4_10(_dst, RCN_CTB_STRIDE,
                                          _src, tb_w);
            _dst += RCN_CTB_STRIDE << 2;
            _src += tb_w << 2;
        }
    } else {
        vvc_add_residual(src, dst, log2_tb_w, log2_tb_h, 0);
    }
}

static void
vvc_add_residual_16_2_10_avx2(const int16_t *const src, uint16_t *const dst,
                            int log2_tb_w, int log2_tb_h, int scale)
{
    int i;
    int tb_w = 1 << log2_tb_w;
    int tb_h = 1 << log2_tb_h;
    const int16_t *_src = (const int16_t *)src;
    uint16_t *_dst = dst;
    if (log2_tb_h > 1) {
      for (i = 0; i < tb_h >> 1; ++i){
          ovvc_transform_add_avx2_16_2_10(_dst, RCN_CTB_STRIDE,
                                         _src, tb_w);
          _dst += RCN_CTB_STRIDE << 1;
          _src += tb_w << 1;
      }
    } else {
        vvc_add_residual(src, dst, log2_tb_w, log2_tb_h, 0);
    }
}

static void
vvc_add_residual_32_1_10_avx2(const int16_t *const src, uint16_t *const dst,
                     int log2_tb_w, int log2_tb_h, int scale)
{
    int i;
    int tb_w = 1 << log2_tb_w;
    int tb_h = 1 << log2_tb_h;
    const int16_t *_src = (const int16_t *)src;
    uint16_t *_dst = dst;
    for (i = 0; i < tb_h; ++i){
        ovvc_transform_add_avx2_32_1_10(_dst, RCN_CTB_STRIDE,
                                       _src, tb_w);
        _dst += RCN_CTB_STRIDE;
        _src += tb_w;
    }
}

static void
vvc_add_residual_64_1_10_avx2(const int16_t *const src, uint16_t *const dst,
                     int log2_tb_w, int log2_tb_h, int scale)
{
    int i;
    int tb_w = 1 << log2_tb_w;
    int tb_h = 1 << log2_tb_h;
    const int16_t *_src = (const int16_t *)src;
    uint16_t *_dst = dst;
    for (i = 0; i < tb_h; ++i){
        ovvc_transform_add_avx2_32_1_10(_dst, RCN_CTB_STRIDE,
                                       _src, tb_w);
        ovvc_transform_add_avx2_32_1_10(_dst+32, RCN_CTB_STRIDE,
                                       _src+32, tb_w);
        _dst += RCN_CTB_STRIDE;
        _src += tb_w;
    }
}

static void
vvc_add_half_residual_8_4_10_avx2(const int16_t *const src, uint16_t *const dst,
                            int log2_tb_w, int log2_tb_h, int scale)
{
    int i;
    int tb_w = 1 << log2_tb_w;
    int tb_h = 1 << log2_tb_h;
    const int16_t *_src = (const int16_t *)src;
    uint16_t *_dst = dst;
    if (log2_tb_h > 1) {
        for (i = 0; i < tb_h >> 2; ++i){
            ovvc_transform_add_half_avx2_8_4_10(_dst, RCN_CTB_STRIDE,
                                          _src, tb_w);
            _dst += RCN_CTB_STRIDE << 2;
            _src += tb_w << 2;
        }
    } else {
        vvc_add_half_residual(src, dst, log2_tb_w, log2_tb_h, 0);
    }
}

static void
vvc_add_half_residual_16_2_10_avx2(const int16_t *const src, uint16_t *const dst,
                            int log2_tb_w, int log2_tb_h, int scale)
{
    int i;
    int tb_w = 1 << log2_tb_w;
    int tb_h = 1 << log2_tb_h;
    const int16_t *_src = (const int16_t *)src;
    uint16_t *_dst = dst;
    if (log2_tb_h > 1) {
      for (i = 0; i < tb_h >> 1; ++i){
          ovvc_transform_add_half_avx2_16_2_10(_dst, RCN_CTB_STRIDE,
                                         _src, tb_w);
          _dst += RCN_CTB_STRIDE << 1;
          _src += tb_w << 1;
      }
    } else {
        vvc_add_half_residual(src, dst, log2_tb_w, log2_tb_h, 0);
    }
}

static void
vvc_add_half_residual_32_1_10_avx2(const int16_t *const src, uint16_t *const dst,
                     int log2_tb_w, int log2_tb_h, int scale)
{
    int i;
    int tb_w = 1 << log2_tb_w;
    int tb_h = 1 << log2_tb_h;
    const int16_t *_src = (const int16_t *)src;
    uint16_t *_dst = dst;
    for (i = 0; i < tb_h; ++i){
        ovvc_transform_add_half_avx2_32_1_10(_dst, RCN_CTB_STRIDE,
                                       _src, tb_w);
        _dst += RCN_CTB_STRIDE;
        _src += tb_w;
    }
}

static void
vvc_sub_residual_8_4_10_avx2(const int16_t *const src, uint16_t *const dst,
                            int log2_tb_w, int log2_tb_h, int scale)
{
    int i;
    int tb_w = 1 << log2_tb_w;
    int tb_h = 1 << log2_tb_h;
    const int16_t *_src = (const int16_t *)src;
    uint16_t *_dst = dst;
    if (log2_tb_h > 1) {
        for (i = 0; i < tb_h >> 2; ++i){
            ovvc_transform_sub_avx2_8_4_10(_dst, RCN_CTB_STRIDE,
                                          _src, tb_w);
            _dst += RCN_CTB_STRIDE << 2;
            _src += tb_w << 2;
        }
    } else {
        vvc_sub_residual(src, dst, log2_tb_w, log2_tb_h, 0);
    }
}

static void
vvc_sub_residual_16_2_10_avx2(const int16_t *const src, uint16_t *const dst,
                            int log2_tb_w, int log2_tb_h, int scale)
{
    int i;
    int tb_w = 1 << log2_tb_w;
    int tb_h = 1 << log2_tb_h;
    const int16_t *_src = (const int16_t *)src;
    uint16_t *_dst = dst;
    if (log2_tb_h > 1) {
      for (i = 0; i < tb_h >> 1; ++i){
          ovvc_transform_sub_avx2_16_2_10(_dst, RCN_CTB_STRIDE,
                                         _src, tb_w);
          _dst += RCN_CTB_STRIDE << 1;
          _src += tb_w << 1;
      }
    } else {
        vvc_sub_residual(src, dst, log2_tb_w, log2_tb_h, 0);
    }
}

static void
vvc_sub_residual_32_1_10_avx2(const int16_t *const src, uint16_t *const dst,
                     int log2_tb_w, int log2_tb_h, int scale)
{
    int i;
    int tb_w = 1 << log2_tb_w;
    int tb_h = 1 << log2_tb_h;
    const int16_t *_src = (const int16_t *)src;
    uint16_t *_dst = dst;
    for (i = 0; i < tb_h; ++i){
        ovvc_transform_sub_avx2_32_1_10(_dst, RCN_CTB_STRIDE,
                                       _src, tb_w);
        _dst += RCN_CTB_STRIDE;
        _src += tb_w;
    }
}

static void
vvc_sub_half_residual_8_4_10_avx2(const int16_t *const src, uint16_t *const dst,
                            int log2_tb_w, int log2_tb_h, int scale)
{
    int i;
    int tb_w = 1 << log2_tb_w;
    int tb_h = 1 << log2_tb_h;
    const int16_t *_src = (const int16_t *)src;
    uint16_t *_dst = dst;
    if (log2_tb_h > 1) {
        for (i = 0; i < tb_h >> 2; ++i){
            ovvc_transform_sub_half_avx2_8_4_10(_dst, RCN_CTB_STRIDE,
                                          _src, tb_w);
            _dst += RCN_CTB_STRIDE << 2;
            _src += tb_w << 2;
        }
    } else {
        vvc_sub_half_residual(src, dst, log2_tb_w, log2_tb_h, 0);
    }
}

static void
vvc_sub_half_residual_16_2_10_avx2(const int16_t *const src, uint16_t *const dst,
                            int log2_tb_w, int log2_tb_h, int scale)
{
    int i;
    int tb_w = 1 << log2_tb_w;
    int tb_h = 1 << log2_tb_h;
    const int16_t *_src = (const int16_t *)src;
    uint16_t *_dst = dst;
    if (log2_tb_h > 1) {
      for (i = 0; i < tb_h >> 1; ++i){
          ovvc_transform_sub_half_avx2_16_2_10(_dst, RCN_CTB_STRIDE,
                                         _src, tb_w);
          _dst += RCN_CTB_STRIDE << 1;
          _src += tb_w << 1;
      }
    } else {
        vvc_sub_half_residual(src, dst, log2_tb_w, log2_tb_h, 0);
    }
}

static void
vvc_sub_half_residual_32_1_10_avx2(const int16_t *const src, uint16_t *const dst,
                     int log2_tb_w, int log2_tb_h, int scale)
{
    int i;
    int tb_w = 1 << log2_tb_w;
    int tb_h = 1 << log2_tb_h;
    const int16_t *_src = (const int16_t *)src;
    uint16_t *_dst = dst;
    for (i = 0; i < tb_h; ++i){
        ovvc_transform_sub_half_avx2_32_1_10(_dst, RCN_CTB_STRIDE,
                                       _src, tb_w);
        _dst += RCN_CTB_STRIDE;
        _src += tb_w;
    }
}

static void
vvc_scale_add_residual_8_4_10_avx2(const int16_t *src, uint16_t *dst,
                      int log2_tb_w, int log2_tb_h,
                      int scale)
{

  int i;
  int tb_w = 1 << log2_tb_w;
  int tb_h = 1 << log2_tb_h;
  const int16_t *_src = (const int16_t *)src;
  uint16_t *_dst = dst;
  if (log2_tb_h > 1) {
    __m256i scale_vector = _mm256_set1_epi16(scale);
    for (i = 0; i < tb_h >> 2; ++i){
        ovvc_transform_scale_add_avx2_8_4_10(_dst, RCN_CTB_STRIDE,
                                       _src, tb_w, scale_vector);
        _dst += RCN_CTB_STRIDE << 2;
        _src += tb_w << 2;
    }
  } else {
      vvc_scale_add_residual(src, dst, log2_tb_w, log2_tb_h, scale);
  }
}

static void
vvc_scale_add_residual_16_2_10_avx2(const int16_t *src, uint16_t *dst,
                      int log2_tb_w, int log2_tb_h,
                      int scale)
{
  int i;
  int tb_w = 1 << log2_tb_w;
  int tb_h = 1 << log2_tb_h;
  const int16_t *_src = (const int16_t *)src;
  uint16_t *_dst = dst;
  if (log2_tb_h > 1) {
    __m256i scale_vector = _mm256_set1_epi16(scale);
    for (i = 0; i < tb_h >> 1; ++i){
        ovvc_transform_scale_add_avx2_16_2_10(_dst, RCN_CTB_STRIDE,
                                       _src, tb_w, scale_vector);
        _dst += RCN_CTB_STRIDE << 1;
        _src += tb_w << 1;
    }
  } else {
      vvc_scale_add_residual(src, dst, log2_tb_w, log2_tb_h, scale);
  }
}

static void
vvc_scale_add_residual_32_1_10_avx2(const int16_t *src, uint16_t *dst,
                      int log2_tb_w, int log2_tb_h,
                      int scale)
{
  int i;
  int tb_w = 1 << log2_tb_w;
  int tb_h = 1 << log2_tb_h;
  const int16_t *_src = (const int16_t *)src;
  uint16_t *_dst = dst;
  __m256i scale_vector = _mm256_set1_epi16(scale);
  for (i = 0; i < tb_h; ++i){
      ovvc_transform_scale_add_avx2_32_1_10(_dst, RCN_CTB_STRIDE,
                                     _src, tb_w, scale_vector);
      _dst += RCN_CTB_STRIDE;
      _src += tb_w;
  }
}

static void
vvc_scale_add_half_residual_8_4_10_avx2(const int16_t *src, uint16_t *dst,
                      int log2_tb_w, int log2_tb_h,
                      int scale)
{
  int i;
  int tb_w = 1 << log2_tb_w;
  int tb_h = 1 << log2_tb_h;
  const int16_t *_src = (const int16_t *)src;
  uint16_t *_dst = dst;
  if (log2_tb_h > 1) {
    __m256i scale_vector = _mm256_set1_epi16(scale);
    for (i = 0; i < tb_h >> 2; ++i){
        ovvc_transform_scale_add_half_avx2_8_4_10(_dst, RCN_CTB_STRIDE,
                                       _src, tb_w, scale_vector);
        _dst += RCN_CTB_STRIDE << 2;
        _src += tb_w << 2;
    }
  } else {
      vvc_scale_add_half_residual(src, dst, log2_tb_w, log2_tb_h, scale);
  }
}

static void
vvc_scale_add_half_residual_16_2_10_avx2(const int16_t *src, uint16_t *dst,
                      int log2_tb_w, int log2_tb_h,
                      int scale)
{
  int i;
  int tb_w = 1 << log2_tb_w;
  int tb_h = 1 << log2_tb_h;
  const int16_t *_src = (const int16_t *)src;
  uint16_t *_dst = dst;
  if (log2_tb_h > 1) {
    __m256i scale_vector = _mm256_set1_epi16(scale);
    for (i = 0; i < tb_h >> 1; ++i){
        ovvc_transform_scale_add_half_avx2_16_2_10(_dst, RCN_CTB_STRIDE,
                                       _src, tb_w, scale_vector);
        _dst += RCN_CTB_STRIDE << 1;
        _src += tb_w << 1;
    }
  } else {
      vvc_scale_add_half_residual(src, dst, log2_tb_w, log2_tb_h, scale);
  }
}

static void
vvc_scale_add_half_residual_32_1_10_avx2(const int16_t *src, uint16_t *dst,
                      int log2_tb_w, int log2_tb_h,
                      int scale)
{
  int i;
  int tb_w = 1 << log2_tb_w;
  int tb_h = 1 << log2_tb_h;
  const int16_t *_src = (const int16_t *)src;
  uint16_t *_dst = dst;
  __m256i scale_vector = _mm256_set1_epi16(scale);
  for (i = 0; i < tb_h; ++i){
      ovvc_transform_scale_add_half_avx2_32_1_10(_dst, RCN_CTB_STRIDE,
                                     _src, tb_w, scale_vector);
      _dst += RCN_CTB_STRIDE;
      _src += tb_w;
  }
}

static void
vvc_scale_sub_residual_8_4_10_avx2(const int16_t *src, uint16_t *dst,
                      int log2_tb_w, int log2_tb_h,
                      int scale)
{
  int i;
  int tb_w = 1 << log2_tb_w;
  int tb_h = 1 << log2_tb_h;
  const int16_t *_src = (const int16_t *)src;
  uint16_t *_dst = dst;
  if (log2_tb_h > 1) {
    __m256i scale_vector = _mm256_set1_epi16(scale);
    for (i = 0; i < tb_h >> 2; ++i){
        ovvc_transform_scale_sub_avx2_8_4_10(_dst, RCN_CTB_STRIDE,
                                       _src, tb_w, scale_vector);
        _dst += RCN_CTB_STRIDE << 2;
        _src += tb_w << 2;
    }
  } else {
      vvc_scale_sub_residual(src, dst, log2_tb_w, log2_tb_h, scale);
  }
}

static void
vvc_scale_sub_residual_16_2_10_avx2(const int16_t *src, uint16_t *dst,
                      int log2_tb_w, int log2_tb_h,
                      int scale)
{
  int i;
  int tb_w = 1 << log2_tb_w;
  int tb_h = 1 << log2_tb_h;
  const int16_t *_src = (const int16_t *)src;
  uint16_t *_dst = dst;
  if (log2_tb_h > 1) {
    __m256i scale_vector = _mm256_set1_epi16(scale);
    for (i = 0; i < tb_h >> 1; ++i){
        ovvc_transform_scale_sub_avx2_16_2_10(_dst, RCN_CTB_STRIDE,
                                       _src, tb_w, scale_vector);
        _dst += RCN_CTB_STRIDE << 1;
        _src += tb_w << 1;
    }
  } else {
      vvc_scale_sub_residual(src, dst, log2_tb_w, log2_tb_h, scale);
  }
}

static void
vvc_scale_sub_residual_32_1_10_avx2(const int16_t *src, uint16_t *dst,
                      int log2_tb_w, int log2_tb_h,
                      int scale)
{
  int i;
  int tb_w = 1 << log2_tb_w;
  int tb_h = 1 << log2_tb_h;
  const int16_t *_src = (const int16_t *)src;
  uint16_t *_dst = dst;
  __m256i scale_vector = _mm256_set1_epi16(scale);
  for (i = 0; i < tb_h; ++i){
      ovvc_transform_scale_sub_avx2_32_1_10(_dst, RCN_CTB_STRIDE,
                                     _src, tb_w, scale_vector);
      _dst += RCN_CTB_STRIDE;
      _src += tb_w;
  }
}

static void
vvc_scale_sub_half_residual_8_4_10_avx2(const int16_t *src, uint16_t *dst,
                      int log2_tb_w, int log2_tb_h,
                      int scale)
{
  int i;
  int tb_w = 1 << log2_tb_w;
  int tb_h = 1 << log2_tb_h;
  const int16_t *_src = (const int16_t *)src;
  uint16_t *_dst = dst;
  if (log2_tb_h > 1) {
    __m256i scale_vector = _mm256_set1_epi16(scale);
    for (i = 0; i < tb_h >> 2; ++i){
        ovvc_transform_scale_sub_half_avx2_8_4_10(_dst, RCN_CTB_STRIDE,
                                       _src, tb_w, scale_vector);
        _dst += RCN_CTB_STRIDE << 2;
        _src += tb_w << 2;
    }
  } else {
      vvc_scale_sub_half_residual(src, dst, log2_tb_w, log2_tb_h, scale);
  }
}

static void
vvc_scale_sub_half_residual_16_2_10_avx2(const int16_t *src, uint16_t *dst,
                      int log2_tb_w, int log2_tb_h,
                      int scale)
{
  int i;
  int tb_w = 1 << log2_tb_w;
  int tb_h = 1 << log2_tb_h;
  const int16_t *_src = (const int16_t *)src;
  uint16_t *_dst = dst;
  if (log2_tb_h > 1) {
    __m256i scale_vector = _mm256_set1_epi16(scale);
    for (i = 0; i < tb_h >> 1; ++i){
        ovvc_transform_scale_sub_half_avx2_16_2_10(_dst, RCN_CTB_STRIDE,
                                       _src, tb_w, scale_vector);
        _dst += RCN_CTB_STRIDE << 1;
        _src += tb_w << 1;
    }
  } else {
      vvc_scale_sub_half_residual(src, dst, log2_tb_w, log2_tb_h, scale);
  }
}

static void
vvc_scale_sub_half_residual_32_1_10_avx2(const int16_t *src, uint16_t *dst,
                      int log2_tb_w, int log2_tb_h,
                      int scale)
{
  int i;
  int tb_w = 1 << log2_tb_w;
  int tb_h = 1 << log2_tb_h;
  const int16_t *_src = (const int16_t *)src;
  uint16_t *_dst = dst;
  __m256i scale_vector = _mm256_set1_epi16(scale);
  for (i = 0; i < tb_h; ++i){
      ovvc_transform_scale_sub_half_avx2_32_1_10(_dst, RCN_CTB_STRIDE,
                                     _src, tb_w, scale_vector);
      _dst += RCN_CTB_STRIDE;
      _src += tb_w;
  }
}

void
rcn_init_ict_functions_avx2(struct RCNFunctions *rcn_func, uint8_t type)
{
    rcn_func->ict.add[3] = &vvc_add_residual_8_4_10_avx2;
    rcn_func->ict.add[4] = &vvc_add_residual_16_2_10_avx2;
    rcn_func->ict.add[5] = &vvc_add_residual_32_1_10_avx2;
    rcn_func->ict.add[6] = &vvc_add_residual_64_1_10_avx2;
    switch (type)
    {
        case 3:
            rcn_func->ict.ict[3][0] = &vvc_scale_add_residual_8_4_10_avx2;
            rcn_func->ict.ict[4][0] = &vvc_scale_add_residual_16_2_10_avx2;
            rcn_func->ict.ict[5][0] = &vvc_scale_add_residual_32_1_10_avx2;

            rcn_func->ict.ict[3][1] = &vvc_scale_sub_residual_8_4_10_avx2;
            rcn_func->ict.ict[4][1] = &vvc_scale_sub_residual_16_2_10_avx2;
            rcn_func->ict.ict[5][1] = &vvc_scale_sub_residual_32_1_10_avx2;

            rcn_func->ict.ict[3][2] = &vvc_scale_sub_half_residual_8_4_10_avx2;
            rcn_func->ict.ict[4][2] = &vvc_scale_sub_half_residual_16_2_10_avx2;
            rcn_func->ict.ict[5][2] = &vvc_scale_sub_half_residual_32_1_10_avx2;
            break;
        case 2:
            rcn_func->ict.ict[3][0] = &vvc_add_residual_8_4_10_avx2;
            rcn_func->ict.ict[4][0] = &vvc_add_residual_16_2_10_avx2;
            rcn_func->ict.ict[5][0] = &vvc_add_residual_32_1_10_avx2;

            rcn_func->ict.ict[3][1] = &vvc_sub_residual_8_4_10_avx2;
            rcn_func->ict.ict[4][1] = &vvc_sub_residual_16_2_10_avx2;
            rcn_func->ict.ict[5][1] = &vvc_sub_residual_32_1_10_avx2;

            rcn_func->ict.ict[3][2] = &vvc_sub_half_residual_8_4_10_avx2;
            rcn_func->ict.ict[4][2] = &vvc_sub_half_residual_16_2_10_avx2;
            rcn_func->ict.ict[5][2] = &vvc_sub_half_residual_32_1_10_avx2;
            break;
        case 1:
            rcn_func->ict.ict[3][0] = &vvc_scale_add_residual_8_4_10_avx2;
            rcn_func->ict.ict[4][0] = &vvc_scale_add_residual_16_2_10_avx2;
            rcn_func->ict.ict[5][0] = &vvc_scale_add_residual_32_1_10_avx2;

            rcn_func->ict.ict[3][1] = &vvc_scale_add_residual_8_4_10_avx2;
            rcn_func->ict.ict[4][1] = &vvc_scale_add_residual_16_2_10_avx2;
            rcn_func->ict.ict[5][1] = &vvc_scale_add_residual_32_1_10_avx2;

            rcn_func->ict.ict[3][2] = &vvc_scale_add_half_residual_8_4_10_avx2;
            rcn_func->ict.ict[4][2] = &vvc_scale_add_half_residual_16_2_10_avx2;
            rcn_func->ict.ict[5][2] = &vvc_scale_add_half_residual_32_1_10_avx2;
            break;
        default:

            rcn_func->ict.ict[3][0] = &vvc_add_residual_8_4_10_avx2;
            rcn_func->ict.ict[4][0] = &vvc_add_residual_16_2_10_avx2;
            rcn_func->ict.ict[5][0] = &vvc_add_residual_32_1_10_avx2;

            rcn_func->ict.ict[3][1] = &vvc_add_residual_8_4_10_avx2;
            rcn_func->ict.ict[4][1] = &vvc_add_residual_16_2_10_avx2;
            rcn_func->ict.ict[5][1] = &vvc_add_residual_32_1_10_avx2;

            rcn_func->ict.ict[3][2] = &vvc_add_half_residual_8_4_10_avx2;
            rcn_func->ict.ict[4][2] = &vvc_add_half_residual_16_2_10_avx2;
            rcn_func->ict.ict[5][2] = &vvc_add_half_residual_32_1_10_avx2;
            break;
    }
}
