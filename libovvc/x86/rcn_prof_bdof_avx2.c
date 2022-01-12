#include <immintrin.h>

#include "rcn_structures.h"
#include "ovutils.h"
#include "bitdepth.h"

#define BITDEPTH 10
#define MAX_PB_SIZE 128

#define GRAD_SHIFT 6
#define PROF_DELTA_LIMIT (1 << (BITDEPTH + 3))
#define SB_H 4
#define SB_W 4
#define PROF_SMP_SHIFT (14 - BITDEPTH)

#define BDOF_SHIFT   (14 + 1 - BITDEPTH)
#define BDOF_OFFSET  ((1 << (BDOF_SHIFT - 1)))

static void
compute_prof_grad_4_avx2(const uint16_t* src, int src_stride, int sb_w, int sb_h,
                  int grad_stride, int16_t* grad_x, int16_t* grad_y)
{
    int y;
    const int nb_smp_h = sb_h;

    src += src_stride + 1;
    __m128i offset = _mm_set1_epi16(1 << 13);
    for (y = 0; y < nb_smp_h; ++y) {
        __m128i x1 = _mm_loadl_epi64((__m128i *)&src[1]);
        __m128i x2 = _mm_loadl_epi64((__m128i *)&src[-1]);
        __m128i y1 = _mm_loadl_epi64((__m128i *)&src[src_stride]);
        __m128i y2 = _mm_loadl_epi64((__m128i *)&src[-src_stride]);
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
        _mm_storel_epi64((__m128i *)grad_x, x1);
        _mm_storel_epi64((__m128i *)grad_y, y1);
        grad_x += grad_stride;
        grad_y += grad_stride;
        src += src_stride;
    }
}

static void
compute_prof_grad_8_avx2(const uint16_t* src, int src_stride, int sb_w, int sb_h,
                  int grad_stride, int16_t* grad_x, int16_t* grad_y)
{
    int y;
    const int nb_smp_h = sb_h;

    src += src_stride + 1;
    __m256i offset = _mm256_set1_epi16(1 << 13);
    for (y = 0; y < nb_smp_h; ++y) {
        __m128i x1 = _mm_loadu_si128((__m128i *)&src[1]);
        __m128i x2 = _mm_loadu_si128((__m128i *)&src[-1]);
        __m128i y1 = _mm_loadu_si128((__m128i *)&src[src_stride]);
        __m128i y2 = _mm_loadu_si128((__m128i *)&src[-src_stride]);

        __m256i xy1 = _mm256_inserti128_si256(_mm256_castsi128_si256(x1), y1, 1);
        __m256i xy2 = _mm256_inserti128_si256(_mm256_castsi128_si256(x2), y2, 1);
        xy1 = _mm256_sub_epi16(xy1, offset);
        xy2 = _mm256_sub_epi16(xy2, offset);
        xy1 = _mm256_srai_epi16(xy1, GRAD_SHIFT);
        xy2 = _mm256_srai_epi16(xy2, GRAD_SHIFT);
        xy1 = _mm256_sub_epi16(xy1, xy2);
        _mm_storeu_si128((__m128i *)grad_x, _mm256_extracti128_si256(xy1, 0));
        _mm_storeu_si128((__m128i *)grad_y, _mm256_extracti128_si256(xy1, 1));
        grad_x += grad_stride;
        grad_y += grad_stride;
        src += src_stride;
    }
}

static void
compute_prof_grad_16_avx2(const uint16_t* src, int src_stride, int sb_w, int sb_h,
                  int grad_stride, int16_t* grad_x, int16_t* grad_y)
{
    int y;
    const int nb_smp_h = sb_h;

    src += src_stride + 1;
    __m256i offset = _mm256_set1_epi16(1 << 13);
    for (y = 0; y < nb_smp_h; ++y) {
        __m256i x1 = _mm256_loadu_si256((__m256i *)&src[1]);
        __m256i x2 = _mm256_loadu_si256((__m256i *)&src[-1]);
        __m256i y1 = _mm256_loadu_si256((__m256i *)&src[src_stride]);
        __m256i y2 = _mm256_loadu_si256((__m256i *)&src[-src_stride]);
        x1 = _mm256_sub_epi16(x1, offset);
        x2 = _mm256_sub_epi16(x2, offset);

        y1 = _mm256_sub_epi16(y1, offset);
        y2 = _mm256_sub_epi16(y2, offset);

        x1 = _mm256_srai_epi16(x1, GRAD_SHIFT);
        x2 = _mm256_srai_epi16(x2, GRAD_SHIFT);

        y1 = _mm256_srai_epi16(y1, GRAD_SHIFT);
        y2 = _mm256_srai_epi16(y2, GRAD_SHIFT);

        x1 = _mm256_sub_epi16(x1, x2);
        y1 = _mm256_sub_epi16(y1, y2);
        _mm256_storeu_si256((__m256i *)grad_x, x1);
        _mm256_storeu_si256((__m256i *)grad_y, y1);
        grad_x += grad_stride;
        grad_y += grad_stride;
        src += src_stride;
    }
}

static void
compute_prof_grad_avx2(const int16_t* src, int src_stride, int sb_w, int sb_h,
                  int grad_stride, int16_t* grad_x, int16_t* grad_y)
{
    if (sb_w == 16) {
      compute_prof_grad_16_avx2(src, src_stride, sb_w, sb_h, grad_stride, grad_x, grad_y);
    } else if (sb_w == 8) {
      compute_prof_grad_8_avx2(src, src_stride, sb_w, sb_h, grad_stride, grad_x, grad_y);
    } else {
      compute_prof_grad_4_avx2(src, src_stride, sb_w, sb_h, grad_stride, grad_x, grad_y);
    }
}


void
rcn_init_prof_functions_avx2(struct RCNFunctions *const rcn_funcs)
{
    rcn_funcs->prof.grad = &compute_prof_grad_avx2;
}

void
rcn_init_bdof_functions_avx2(struct RCNFunctions *const rcn_funcs)
{
    rcn_funcs->bdof.grad = &compute_prof_grad_avx2;
}
