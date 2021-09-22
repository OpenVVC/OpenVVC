#include <tmmintrin.h>

#include "rcn_structures.h"

#define GRAD_SHIFT 6
#define PROF_DELTA_LIMIT (1 << (BITDEPTH + 3))
static void
compute_prof_grad_4_sse(const uint16_t* src, int src_stride, int sb_w, int sb_h,
                  int grad_stride, int16_t* grad_x, int16_t* grad_y)
{
    int y, x;
    const int nb_smp_h = sb_h;
    const int nb_smp_w = sb_w;

    src += src_stride + 1;
    __m128i offset = _mm_set1_epi16(1 << 13);
    for (y = 0; y < nb_smp_h; ++y) {
        __m128i x1 = _mm_loadu_si64((__m128i *)&src[1]);
        __m128i x2 = _mm_loadu_si64((__m128i *)&src[-1]);
        __m128i y1 = _mm_loadu_si64((__m128i *)&src[src_stride]);
        __m128i y2 = _mm_loadu_si64((__m128i *)&src[-src_stride]);
        x1 = _mm_sub_epi16(x1, offset);
        x2 = _mm_sub_epi16(x2, offset);
        y1 = _mm_sub_epi16(y1, offset);
        y2 = _mm_sub_epi16(y2, offset);
        x1 = _mm_srai_epi16(x1, GRAD_SHIFT);
        x2 = _mm_srai_epi16(x2, GRAD_SHIFT);
        y1 = _mm_srai_epi16(y1, GRAD_SHIFT);
        y2 = _mm_srai_epi16(y2, GRAD_SHIFT);
        x1 = _mm_sub_epi16(x1, x2);
        y1 = _mm_sub_epi16(y1, y2);
        _mm_storeu_si64((__m128i *)grad_x, x1);
        _mm_storeu_si64((__m128i *)grad_y, y1);
        grad_x += grad_stride;
        grad_y += grad_stride;
        src += src_stride;
    }
}

static void
compute_prof_grad_8_sse(const uint16_t* src, int src_stride, int sb_w, int sb_h,
                  int grad_stride, int16_t* grad_x, int16_t* grad_y)
{
    int y, x;
    const int nb_smp_h = sb_h;
    const int nb_smp_w = sb_w;

    src += src_stride + 1;
    __m128i offset = _mm_set1_epi16(1 << 13);
    for (y = 0; y < nb_smp_h; ++y) {
        __m128i x1 = _mm_loadu_si128((__m128i *)&src[1]);
        __m128i x2 = _mm_loadu_si128((__m128i *)&src[-1]);
        __m128i y1 = _mm_loadu_si128((__m128i *)&src[src_stride]);
        __m128i y2 = _mm_loadu_si128((__m128i *)&src[-src_stride]);
        x1 = _mm_sub_epi16(x1, offset);
        x2 = _mm_sub_epi16(x2, offset);
        y1 = _mm_sub_epi16(y1, offset);
        y2 = _mm_sub_epi16(y2, offset);
        x1 = _mm_srai_epi16(x1, GRAD_SHIFT);
        x2 = _mm_srai_epi16(x2, GRAD_SHIFT);
        y1 = _mm_srai_epi16(y1, GRAD_SHIFT);
        y2 = _mm_srai_epi16(y2, GRAD_SHIFT);
        x1 = _mm_sub_epi16(x1, x2);
        y1 = _mm_sub_epi16(y1, y2);
        _mm_storeu_si128((__m128i *)grad_x, x1);
        _mm_storeu_si128((__m128i *)grad_y, y1);
        grad_x += grad_stride;
        grad_y += grad_stride;
        src += src_stride;
    }
}

static void
compute_prof_grad_16_sse(const uint16_t* src, int src_stride, int sb_w, int sb_h,
                  int grad_stride, int16_t* grad_x, int16_t* grad_y)
{
    int y, x;
    const int nb_smp_h = sb_h;
    const int nb_smp_w = sb_w;

    src += src_stride + 1;
    __m128i offset = _mm_set1_epi16(1 << 13);
    for (y = 0; y < nb_smp_h; ++y) {
        __m128i x1 = _mm_loadu_si128((__m128i *)&src[1]);
        __m128i x2 = _mm_loadu_si128((__m128i *)&src[-1]);
        __m128i x3 = _mm_loadu_si128((__m128i *)&src[1+8]);
        __m128i x4 = _mm_loadu_si128((__m128i *)&src[-1+8]);
        __m128i y1 = _mm_loadu_si128((__m128i *)&src[src_stride]);
        __m128i y2 = _mm_loadu_si128((__m128i *)&src[-src_stride]);
        __m128i y3 = _mm_loadu_si128((__m128i *)&src[src_stride+8]);
        __m128i y4 = _mm_loadu_si128((__m128i *)&src[-src_stride+8]);
        x1 = _mm_sub_epi16(x1, offset);
        x2 = _mm_sub_epi16(x2, offset);
        x3 = _mm_sub_epi16(x3, offset);
        x4 = _mm_sub_epi16(x4, offset);
        y1 = _mm_sub_epi16(y1, offset);
        y2 = _mm_sub_epi16(y2, offset);
        y3 = _mm_sub_epi16(y3, offset);
        y4 = _mm_sub_epi16(y4, offset);
        x1 = _mm_srai_epi16(x1, GRAD_SHIFT);
        x2 = _mm_srai_epi16(x2, GRAD_SHIFT);
        x3 = _mm_srai_epi16(x3, GRAD_SHIFT);
        x4 = _mm_srai_epi16(x4, GRAD_SHIFT);
        y1 = _mm_srai_epi16(y1, GRAD_SHIFT);
        y2 = _mm_srai_epi16(y2, GRAD_SHIFT);
        y3 = _mm_srai_epi16(y3, GRAD_SHIFT);
        y4 = _mm_srai_epi16(y4, GRAD_SHIFT);
        x1 = _mm_sub_epi16(x1, x2);
        x3 = _mm_sub_epi16(x3, x4);
        y1 = _mm_sub_epi16(y1, y2);
        y3 = _mm_sub_epi16(y3, y4);
        _mm_storeu_si128((__m128i *)grad_x, x1);
        _mm_storeu_si128((__m128i *)&grad_x[8], x3);
        _mm_storeu_si128((__m128i *)grad_y, y1);
        _mm_storeu_si128((__m128i *)&grad_y[8], y3);
        grad_x += grad_stride;
        grad_y += grad_stride;
        src += src_stride;
    }
}

static void
compute_prof_grad_sse(const uint16_t* src, int src_stride, int sb_w, int sb_h,
                  int grad_stride, int16_t* grad_x, int16_t* grad_y)
{
    if (sb_w == 16) {
      compute_prof_grad_16_sse(src, src_stride, sb_w, sb_h, grad_stride, grad_x, grad_y);
    } else if (sb_w == 8) {
      compute_prof_grad_8_sse(src, src_stride, sb_w, sb_h, grad_stride, grad_x, grad_y);
    } else {
      compute_prof_grad_4_sse(src, src_stride, sb_w, sb_h, grad_stride, grad_x, grad_y);
    }
}

void
rcn_prof_functions_sse(struct RCNFunctions *const rcn_funcs)
{
    rcn_funcs->prof.grad = &compute_prof_grad_sse;
}
