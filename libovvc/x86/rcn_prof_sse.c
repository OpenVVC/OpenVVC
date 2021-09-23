#include <smmintrin.h>

#include "rcn_structures.h"
#include "ovutils.h"

#define GRAD_SHIFT 6
#define PROF_DELTA_LIMIT (1 << (BITDEPTH + 3))
#define BITDEPTH 10
#define SB_H 4
#define SB_W 4
#define PROF_SMP_SHIFT (14 - BITDEPTH)

static void rcn_prof_sse(uint16_t* dst, int dst_stride, const uint16_t* src, int src_stride,
         const int16_t* grad_x, const int16_t* grad_y, int grad_stride,
         const int32_t* dmv_scale_h, const int32_t* dmv_scale_v,
         uint8_t bidir)
{
    //FIXME: Convert dmv_scale to int16_t to avoid _mm_packs_epi32
    int idx = 0;
    int x, y;
    if (!bidir) {
        __m128i sh1 = _mm_loadu_si128((__m128i *)&dmv_scale_h[0]);
        __m128i sh2 = _mm_loadu_si128((__m128i *)&dmv_scale_h[4]);
        __m128i sh3 = _mm_loadu_si128((__m128i *)&dmv_scale_h[8]);
        __m128i sh4 = _mm_loadu_si128((__m128i *)&dmv_scale_h[12]);

        __m128i sv1 = _mm_loadu_si128((__m128i *)&dmv_scale_v[0]);
        __m128i sv2 = _mm_loadu_si128((__m128i *)&dmv_scale_v[4]);
        __m128i sv3 = _mm_loadu_si128((__m128i *)&dmv_scale_v[8]);
        __m128i sv4 = _mm_loadu_si128((__m128i *)&dmv_scale_v[12]);

        __m128i x1 = _mm_loadu_si64((__m128i *)&grad_x[0*grad_stride]);
        __m128i x2 = _mm_loadu_si64((__m128i *)&grad_x[1*grad_stride]);
        __m128i x3 = _mm_loadu_si64((__m128i *)&grad_x[2*grad_stride]);
        __m128i x4 = _mm_loadu_si64((__m128i *)&grad_x[3*grad_stride]);

        __m128i y1 = _mm_loadu_si64((__m128i *)&grad_y[0*grad_stride]);
        __m128i y2 = _mm_loadu_si64((__m128i *)&grad_y[1*grad_stride]);
        __m128i y3 = _mm_loadu_si64((__m128i *)&grad_y[2*grad_stride]);
        __m128i y4 = _mm_loadu_si64((__m128i *)&grad_y[3*grad_stride]);

        __m128i srcV1 = _mm_loadu_si64((__m128i *)&src[0*src_stride]);
        __m128i srcV2 = _mm_loadu_si64((__m128i *)&src[1*src_stride]);
        __m128i srcV3 = _mm_loadu_si64((__m128i *)&src[2*src_stride]);
        __m128i srcV4 = _mm_loadu_si64((__m128i *)&src[3*src_stride]);

        sh1 = _mm_packs_epi32(sh1, _mm_setzero_si128());
        sh2 = _mm_packs_epi32(sh2, _mm_setzero_si128());
        sh3 = _mm_packs_epi32(sh3, _mm_setzero_si128());
        sh4 = _mm_packs_epi32(sh4, _mm_setzero_si128());

        sv1 = _mm_packs_epi32(sv1, _mm_setzero_si128());
        sv2 = _mm_packs_epi32(sv2, _mm_setzero_si128());
        sv3 = _mm_packs_epi32(sv3, _mm_setzero_si128());
        sv4 = _mm_packs_epi32(sv4, _mm_setzero_si128());

        srcV1 = _mm_unpacklo_epi16(srcV1, _mm_setzero_si128());
        srcV2 = _mm_unpacklo_epi16(srcV2, _mm_setzero_si128());
        srcV3 = _mm_unpacklo_epi16(srcV3, _mm_setzero_si128());
        srcV4 = _mm_unpacklo_epi16(srcV4, _mm_setzero_si128());

        // int32_t add = dmv_scale_h[idx] * grad_x[x] + dmv_scale_v[idx] * grad_y[x];
        x1 = _mm_unpacklo_epi16(x1, y1);
        x2 = _mm_unpacklo_epi16(x2, y2);
        x3 = _mm_unpacklo_epi16(x3, y3);
        x4 = _mm_unpacklo_epi16(x4, y4);

        sh1 = _mm_unpacklo_epi16(sh1, sv1);
        sh2 = _mm_unpacklo_epi16(sh2, sv2);
        sh3 = _mm_unpacklo_epi16(sh3, sv3);
        sh4 = _mm_unpacklo_epi16(sh4, sv4);

        x1 = _mm_madd_epi16(x1, sh1);
        x2 = _mm_madd_epi16(x2, sh2);
        x3 = _mm_madd_epi16(x3, sh3);
        x4 = _mm_madd_epi16(x4, sh4);

        // add = ov_clip(add, -PROF_DELTA_LIMIT, PROF_DELTA_LIMIT - 1);
        __m128i min_val = _mm_set1_epi32(-PROF_DELTA_LIMIT);
        __m128i max_val = _mm_set1_epi32(PROF_DELTA_LIMIT - 1);
        x1 = _mm_max_epi32(x1, min_val);
        x2 = _mm_max_epi32(x2, min_val);
        x3 = _mm_max_epi32(x3, min_val);
        x4 = _mm_max_epi32(x4, min_val);

        x1 = _mm_min_epi32(x1, max_val);
        x2 = _mm_min_epi32(x2, max_val);
        x3 = _mm_min_epi32(x3, max_val);
        x4 = _mm_min_epi32(x4, max_val);

        // val = (int16_t)src[x] -(1 << 13);
        // val += add;
        //
        // /* Clipping if not bi directional */
        // val = (val + 8200 /*+ PROF_SMP_OFFSET*/) >> PROF_SMP_SHIFT;
        // dst[x] = ov_clip(val, 0, 1023);

        __m128i offset = _mm_set1_epi32(-8200+(1 << 13));
        srcV1 = _mm_sub_epi32(srcV1, offset);
        srcV2 = _mm_sub_epi32(srcV2, offset);
        srcV3 = _mm_sub_epi32(srcV3, offset);
        srcV4 = _mm_sub_epi32(srcV4, offset);

        x1 = _mm_add_epi32(srcV1, x1);
        x2 = _mm_add_epi32(srcV2, x2);
        x3 = _mm_add_epi32(srcV3, x3);
        x4 = _mm_add_epi32(srcV4, x4);

        x1 = _mm_srai_epi16(x1, PROF_SMP_SHIFT);
        x2 = _mm_srai_epi16(x2, PROF_SMP_SHIFT);
        x3 = _mm_srai_epi16(x3, PROF_SMP_SHIFT);
        x4 = _mm_srai_epi16(x4, PROF_SMP_SHIFT);

        //pack without saturation
        x1 = _mm_shufflelo_epi16(x1, 0x88);
        x2 = _mm_shufflelo_epi16(x2, 0x88);
        x3 = _mm_shufflelo_epi16(x3, 0x88);
        x4 = _mm_shufflelo_epi16(x4, 0x88);

        x1 = _mm_shufflehi_epi16(x1, 0x88);
        x2 = _mm_shufflehi_epi16(x2, 0x88);
        x3 = _mm_shufflehi_epi16(x3, 0x88);
        x4 = _mm_shufflehi_epi16(x4, 0x88);

        x1 = _mm_shuffle_epi32(x1, 0x88);
        x2 = _mm_shuffle_epi32(x2, 0x88);
        x3 = _mm_shuffle_epi32(x3, 0x88);
        x4 = _mm_shuffle_epi32(x4, 0x88);

        x1 = _mm_max_epi16(x1, _mm_setzero_si128());
        x2 = _mm_max_epi16(x2, _mm_setzero_si128());
        x3 = _mm_max_epi16(x3, _mm_setzero_si128());
        x4 = _mm_max_epi16(x4, _mm_setzero_si128());

        x1 = _mm_min_epi16(x1, _mm_set1_epi16(1023));
        x2 = _mm_min_epi16(x2, _mm_set1_epi16(1023));
        x3 = _mm_min_epi16(x3, _mm_set1_epi16(1023));
        x4 = _mm_min_epi16(x4, _mm_set1_epi16(1023));

        _mm_storeu_si64((__m128i *)&dst[0*dst_stride], x1);
        _mm_storeu_si64((__m128i *)&dst[1*dst_stride], x2);
        _mm_storeu_si64((__m128i *)&dst[2*dst_stride], x3);
        _mm_storeu_si64((__m128i *)&dst[3*dst_stride], x4);
    } else {
      __m128i min_val = _mm_set1_epi32(-PROF_DELTA_LIMIT);
      __m128i max_val = _mm_set1_epi32(PROF_DELTA_LIMIT - 1);

      __m128i sh1 = _mm_loadu_si128((__m128i *)&dmv_scale_h[0]);
      __m128i sh2 = _mm_loadu_si128((__m128i *)&dmv_scale_h[4]);
      __m128i sh3 = _mm_loadu_si128((__m128i *)&dmv_scale_h[8]);
      __m128i sh4 = _mm_loadu_si128((__m128i *)&dmv_scale_h[12]);

      __m128i sv1 = _mm_loadu_si128((__m128i *)&dmv_scale_v[0]);
      __m128i sv2 = _mm_loadu_si128((__m128i *)&dmv_scale_v[4]);
      __m128i sv3 = _mm_loadu_si128((__m128i *)&dmv_scale_v[8]);
      __m128i sv4 = _mm_loadu_si128((__m128i *)&dmv_scale_v[12]);

      __m128i x1 = _mm_loadu_si64((__m128i *)&grad_x[0*grad_stride]);
      __m128i x2 = _mm_loadu_si64((__m128i *)&grad_x[1*grad_stride]);
      __m128i x3 = _mm_loadu_si64((__m128i *)&grad_x[2*grad_stride]);
      __m128i x4 = _mm_loadu_si64((__m128i *)&grad_x[3*grad_stride]);

      __m128i y1 = _mm_loadu_si64((__m128i *)&grad_y[0*grad_stride]);
      __m128i y2 = _mm_loadu_si64((__m128i *)&grad_y[1*grad_stride]);
      __m128i y3 = _mm_loadu_si64((__m128i *)&grad_y[2*grad_stride]);
      __m128i y4 = _mm_loadu_si64((__m128i *)&grad_y[3*grad_stride]);

      __m128i srcV1 = _mm_loadu_si64((__m128i *)&src[0*src_stride]);
      __m128i srcV2 = _mm_loadu_si64((__m128i *)&src[1*src_stride]);
      __m128i srcV3 = _mm_loadu_si64((__m128i *)&src[2*src_stride]);
      __m128i srcV4 = _mm_loadu_si64((__m128i *)&src[3*src_stride]);

      sh1 = _mm_packs_epi32(sh1, _mm_setzero_si128());
      sh2 = _mm_packs_epi32(sh2, _mm_setzero_si128());
      sh3 = _mm_packs_epi32(sh3, _mm_setzero_si128());
      sh4 = _mm_packs_epi32(sh4, _mm_setzero_si128());

      sv1 = _mm_packs_epi32(sv1, _mm_setzero_si128());
      sv2 = _mm_packs_epi32(sv2, _mm_setzero_si128());
      sv3 = _mm_packs_epi32(sv3, _mm_setzero_si128());
      sv4 = _mm_packs_epi32(sv4, _mm_setzero_si128());

      srcV1 = _mm_unpacklo_epi16(srcV1, _mm_setzero_si128());
      srcV2 = _mm_unpacklo_epi16(srcV2, _mm_setzero_si128());
      srcV3 = _mm_unpacklo_epi16(srcV3, _mm_setzero_si128());
      srcV4 = _mm_unpacklo_epi16(srcV4, _mm_setzero_si128());

      // int32_t add = dmv_scale_h[idx] * grad_x[x] + dmv_scale_v[idx] * grad_y[x];
      x1 = _mm_unpacklo_epi16(x1, y1);
      x2 = _mm_unpacklo_epi16(x2, y2);
      x3 = _mm_unpacklo_epi16(x3, y3);
      x4 = _mm_unpacklo_epi16(x4, y4);

      sh1 = _mm_unpacklo_epi16(sh1, sv1);
      sh2 = _mm_unpacklo_epi16(sh2, sv2);
      sh3 = _mm_unpacklo_epi16(sh3, sv3);
      sh4 = _mm_unpacklo_epi16(sh4, sv4);

      x1 = _mm_madd_epi16(x1, sh1);
      x2 = _mm_madd_epi16(x2, sh2);
      x3 = _mm_madd_epi16(x3, sh3);
      x4 = _mm_madd_epi16(x4, sh4);

      // add = ov_clip(add, -PROF_DELTA_LIMIT, PROF_DELTA_LIMIT - 1);
      x1 = _mm_max_epi32(x1, min_val);
      x2 = _mm_max_epi32(x2, min_val);
      x3 = _mm_max_epi32(x3, min_val);
      x4 = _mm_max_epi32(x4, min_val);

      x1 = _mm_min_epi32(x1, max_val);
      x2 = _mm_min_epi32(x2, max_val);
      x3 = _mm_min_epi32(x3, max_val);
      x4 = _mm_min_epi32(x4, max_val);

      // dst[x] = src[x] + add;
      x1 = _mm_add_epi32(srcV1, x1);
      x2 = _mm_add_epi32(srcV2, x2);
      x3 = _mm_add_epi32(srcV3, x3);
      x4 = _mm_add_epi32(srcV4, x4);

      //pack without saturation
      x1 = _mm_shufflelo_epi16(x1, 0x88);
      x2 = _mm_shufflelo_epi16(x2, 0x88);
      x3 = _mm_shufflelo_epi16(x3, 0x88);
      x4 = _mm_shufflelo_epi16(x4, 0x88);

      x1 = _mm_shufflehi_epi16(x1, 0x88);
      x2 = _mm_shufflehi_epi16(x2, 0x88);
      x3 = _mm_shufflehi_epi16(x3, 0x88);
      x4 = _mm_shufflehi_epi16(x4, 0x88);

      x1 = _mm_shuffle_epi32(x1, 0x88);
      x2 = _mm_shuffle_epi32(x2, 0x88);
      x3 = _mm_shuffle_epi32(x3, 0x88);
      x4 = _mm_shuffle_epi32(x4, 0x88);

      _mm_storeu_si64((__m128i *)&dst[0*dst_stride], x1);
      _mm_storeu_si64((__m128i *)&dst[1*dst_stride], x2);
      _mm_storeu_si64((__m128i *)&dst[2*dst_stride], x3);
      _mm_storeu_si64((__m128i *)&dst[3*dst_stride], x4);
    }
}


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
    rcn_funcs->prof.rcn = &rcn_prof_sse;
}
