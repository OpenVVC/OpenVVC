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


static void rcn_prof_sse(OVSample* dst, int dst_stride, const int16_t* src, int src_stride,
         const int16_t* grad_x, const int16_t* grad_y, int grad_stride,
         const int32_t* dmv_scale_h, const int32_t* dmv_scale_v,
         uint8_t bidir)
{
    //FIXME: Convert dmv_scale to int16_t to avoid _mm_packs_epi32
    if (!bidir) {
        __m128i sh1 = _mm_loadu_si128((__m128i *)&dmv_scale_h[0]);
        __m128i sh2 = _mm_loadu_si128((__m128i *)&dmv_scale_h[4]);
        __m128i sh3 = _mm_loadu_si128((__m128i *)&dmv_scale_h[8]);
        __m128i sh4 = _mm_loadu_si128((__m128i *)&dmv_scale_h[12]);

        __m128i sv1 = _mm_loadu_si128((__m128i *)&dmv_scale_v[0]);
        __m128i sv2 = _mm_loadu_si128((__m128i *)&dmv_scale_v[4]);
        __m128i sv3 = _mm_loadu_si128((__m128i *)&dmv_scale_v[8]);
        __m128i sv4 = _mm_loadu_si128((__m128i *)&dmv_scale_v[12]);

        __m128i x1 = _mm_loadl_epi64((__m128i *)&grad_x[0*grad_stride]);
        __m128i x2 = _mm_loadl_epi64((__m128i *)&grad_x[1*grad_stride]);
        __m128i x3 = _mm_loadl_epi64((__m128i *)&grad_x[2*grad_stride]);
        __m128i x4 = _mm_loadl_epi64((__m128i *)&grad_x[3*grad_stride]);

        __m128i y1 = _mm_loadl_epi64((__m128i *)&grad_y[0*grad_stride]);
        __m128i y2 = _mm_loadl_epi64((__m128i *)&grad_y[1*grad_stride]);
        __m128i y3 = _mm_loadl_epi64((__m128i *)&grad_y[2*grad_stride]);
        __m128i y4 = _mm_loadl_epi64((__m128i *)&grad_y[3*grad_stride]);

        __m128i srcV1 = _mm_loadl_epi64((__m128i *)&src[0*src_stride]);
        __m128i srcV2 = _mm_loadl_epi64((__m128i *)&src[1*src_stride]);
        __m128i srcV3 = _mm_loadl_epi64((__m128i *)&src[2*src_stride]);
        __m128i srcV4 = _mm_loadl_epi64((__m128i *)&src[3*src_stride]);

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

        _mm_storel_epi64((__m128i *)&dst[0*dst_stride], x1);
        _mm_storel_epi64((__m128i *)&dst[1*dst_stride], x2);
        _mm_storel_epi64((__m128i *)&dst[2*dst_stride], x3);
        _mm_storel_epi64((__m128i *)&dst[3*dst_stride], x4);
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

      __m128i x1 = _mm_loadl_epi64((__m128i *)&grad_x[0*grad_stride]);
      __m128i x2 = _mm_loadl_epi64((__m128i *)&grad_x[1*grad_stride]);
      __m128i x3 = _mm_loadl_epi64((__m128i *)&grad_x[2*grad_stride]);
      __m128i x4 = _mm_loadl_epi64((__m128i *)&grad_x[3*grad_stride]);

      __m128i y1 = _mm_loadl_epi64((__m128i *)&grad_y[0*grad_stride]);
      __m128i y2 = _mm_loadl_epi64((__m128i *)&grad_y[1*grad_stride]);
      __m128i y3 = _mm_loadl_epi64((__m128i *)&grad_y[2*grad_stride]);
      __m128i y4 = _mm_loadl_epi64((__m128i *)&grad_y[3*grad_stride]);

      __m128i srcV1 = _mm_loadl_epi64((__m128i *)&src[0*src_stride]);
      __m128i srcV2 = _mm_loadl_epi64((__m128i *)&src[1*src_stride]);
      __m128i srcV3 = _mm_loadl_epi64((__m128i *)&src[2*src_stride]);
      __m128i srcV4 = _mm_loadl_epi64((__m128i *)&src[3*src_stride]);

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

      _mm_storel_epi64((__m128i *)&dst[0*dst_stride], x1);
      _mm_storel_epi64((__m128i *)&dst[1*dst_stride], x2);
      _mm_storel_epi64((__m128i *)&dst[2*dst_stride], x3);
      _mm_storel_epi64((__m128i *)&dst[3*dst_stride], x4);
    }
}


static void
compute_prof_grad_4_sse(const uint16_t* src, int src_stride, int sb_w, int sb_h,
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
compute_prof_grad_8_sse(const uint16_t* src, int src_stride, int sb_w, int sb_h,
                  int grad_stride, int16_t* grad_x, int16_t* grad_y)
{
    int y;
    const int nb_smp_h = sb_h;

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
    int y;
    const int nb_smp_h = sb_h;

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
compute_prof_grad_sse(const int16_t* src, int src_stride, int sb_w, int sb_h,
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

static void
tmp_prof_mrg_sse(OVSample* _dst, ptrdiff_t _dststride,
             const int16_t* _src0, ptrdiff_t _srcstride,
             const int16_t* _src1, int height, intptr_t mx,
             intptr_t my, int width)
{
    int x, y;
    const int16_t* src0 = (int16_t *)_src0;
    const int16_t* src1 = (int16_t *)_src1;
    ptrdiff_t srcstride = _srcstride;
    OVSample* dst = (OVSample*)_dst;
    ptrdiff_t dststride = _dststride;
    int shift = 14 - BITDEPTH + 1;
    
    __m128i src00 = _mm_loadl_epi64((__m128i *)&src0[0*srcstride]);
    __m128i src01 = _mm_loadl_epi64((__m128i *)&src0[1*srcstride]);
    __m128i src02 = _mm_loadl_epi64((__m128i *)&src0[2*srcstride]);
    __m128i src03 = _mm_loadl_epi64((__m128i *)&src0[3*srcstride]);

    __m128i src10 = _mm_loadl_epi64((__m128i *)&src1[0*MAX_PB_SIZE]);
    __m128i src11 = _mm_loadl_epi64((__m128i *)&src1[1*MAX_PB_SIZE]);
    __m128i src12 = _mm_loadl_epi64((__m128i *)&src1[2*MAX_PB_SIZE]);
    __m128i src13 = _mm_loadl_epi64((__m128i *)&src1[3*MAX_PB_SIZE]);

    // increase to 32x4 to avoid overflow
    src00 = _mm_cvtepi16_epi32(src00);
    src01 = _mm_cvtepi16_epi32(src01);
    src02 = _mm_cvtepi16_epi32(src02);
    src03 = _mm_cvtepi16_epi32(src03);

    src10 = _mm_cvtepi16_epi32(src10);
    src11 = _mm_cvtepi16_epi32(src11);
    src12 = _mm_cvtepi16_epi32(src12);
    src13 = _mm_cvtepi16_epi32(src13);

    __m128i offset = _mm_set1_epi32(((1 << (13 - BITDEPTH + 1))));

    src00 = _mm_add_epi32(src00, src10);
    src01 = _mm_add_epi32(src01, src11);
    src02 = _mm_add_epi32(src02, src12);
    src03 = _mm_add_epi32(src03, src13);

    src00 = _mm_add_epi32(src00, offset);
    src01 = _mm_add_epi32(src01, offset);
    src02 = _mm_add_epi32(src02, offset);
    src03 = _mm_add_epi32(src03, offset);

    src00 = _mm_srai_epi32(src00, shift);
    src01 = _mm_srai_epi32(src01, shift);
    src02 = _mm_srai_epi32(src02, shift);
    src03 = _mm_srai_epi32(src03, shift);

    src00 = _mm_max_epi32(src00, _mm_setzero_si128());
    src01 = _mm_max_epi32(src01, _mm_setzero_si128());
    src02 = _mm_max_epi32(src02, _mm_setzero_si128());
    src03 = _mm_max_epi32(src03, _mm_setzero_si128());

    src00 = _mm_min_epi32(src00, _mm_set1_epi32(1023));
    src01 = _mm_min_epi32(src01, _mm_set1_epi32(1023));
    src02 = _mm_min_epi32(src02, _mm_set1_epi32(1023));
    src03 = _mm_min_epi32(src03, _mm_set1_epi32(1023));

    src00 = _mm_packs_epi32(src00, _mm_setzero_si128());
    src01 = _mm_packs_epi32(src01, _mm_setzero_si128());
    src02 = _mm_packs_epi32(src02, _mm_setzero_si128());
    src03 = _mm_packs_epi32(src03, _mm_setzero_si128());

    _mm_storel_epi64((__m128i *)&dst[0*dststride], src00);
    _mm_storel_epi64((__m128i *)&dst[1*dststride], src01);
    _mm_storel_epi64((__m128i *)&dst[2*dststride], src02);
    _mm_storel_epi64((__m128i *)&dst[3*dststride], src03);
}

static void
tmp_prof_mrg_w_sse(OVSample* _dst, ptrdiff_t _dststride,
               const int16_t* _src0, ptrdiff_t _srcstride,
               const int16_t* _src1, int height, intptr_t mx,
               intptr_t my, int width, int wt0, int wt1)
{
    int x, y;
    const int16_t* src0 = (int16_t *)_src0;
    const int16_t* src1 = (int16_t *)_src1;
    ptrdiff_t srcstride = _srcstride;
    OVSample* dst = (OVSample*)_dst;
    ptrdiff_t dststride = _dststride;
    int log_weights = floor_log2(wt0 + wt1);
    int shift = 14 - BITDEPTH + log_weights;

    __m128i offset = _mm_set1_epi32(2*((1 << (13 - BITDEPTH))) << (log_weights - 1));
    __m128i wt = _mm_set1_epi32(wt0&0xFFFF | (wt1<<16));

    __m128i src00 = _mm_loadl_epi64((__m128i *)&src0[0*srcstride]);
    __m128i src01 = _mm_loadl_epi64((__m128i *)&src0[1*srcstride]);
    __m128i src02 = _mm_loadl_epi64((__m128i *)&src0[2*srcstride]);
    __m128i src03 = _mm_loadl_epi64((__m128i *)&src0[3*srcstride]);

    __m128i src10 = _mm_loadl_epi64((__m128i *)&src1[0*MAX_PB_SIZE]);
    __m128i src11 = _mm_loadl_epi64((__m128i *)&src1[1*MAX_PB_SIZE]);
    __m128i src12 = _mm_loadl_epi64((__m128i *)&src1[2*MAX_PB_SIZE]);
    __m128i src13 = _mm_loadl_epi64((__m128i *)&src1[3*MAX_PB_SIZE]);

    src00 = _mm_unpacklo_epi16(src00, src10);
    src01 = _mm_unpacklo_epi16(src01, src11);
    src02 = _mm_unpacklo_epi16(src02, src12);
    src03 = _mm_unpacklo_epi16(src03, src13);

    src00 = _mm_madd_epi16(src00, wt);
    src01 = _mm_madd_epi16(src01, wt);
    src02 = _mm_madd_epi16(src02, wt);
    src03 = _mm_madd_epi16(src03, wt);

    src00 = _mm_add_epi32(src00, offset);
    src01 = _mm_add_epi32(src01, offset);
    src02 = _mm_add_epi32(src02, offset);
    src03 = _mm_add_epi32(src03, offset);

    src00 = _mm_srai_epi32(src00, shift);
    src01 = _mm_srai_epi32(src01, shift);
    src02 = _mm_srai_epi32(src02, shift);
    src03 = _mm_srai_epi32(src03, shift);

    src00 = _mm_max_epi32(src00, _mm_setzero_si128());
    src01 = _mm_max_epi32(src01, _mm_setzero_si128());
    src02 = _mm_max_epi32(src02, _mm_setzero_si128());
    src03 = _mm_max_epi32(src03, _mm_setzero_si128());

    src00 = _mm_min_epi32(src00, _mm_set1_epi32(1023));
    src01 = _mm_min_epi32(src01, _mm_set1_epi32(1023));
    src02 = _mm_min_epi32(src02, _mm_set1_epi32(1023));
    src03 = _mm_min_epi32(src03, _mm_set1_epi32(1023));

    src00 = _mm_packs_epi32(src00, _mm_setzero_si128());
    src01 = _mm_packs_epi32(src01, _mm_setzero_si128());
    src02 = _mm_packs_epi32(src02, _mm_setzero_si128());
    src03 = _mm_packs_epi32(src03, _mm_setzero_si128());

    _mm_storel_epi64((__m128i *)&dst[0*dststride], src00);
    _mm_storel_epi64((__m128i *)&dst[1*dststride], src01);
    _mm_storel_epi64((__m128i *)&dst[2*dststride], src02);
    _mm_storel_epi64((__m128i *)&dst[3*dststride], src03);
}

static void
rcn_apply_bdof_subblock_sse(const int16_t* src0, int src0_stride,
                        const int16_t* src1, int src1_stride,
                        OVSample *dst, int dst_stride,
                        const int16_t *gradX0, const int16_t *gradX1,
                        const int16_t *gradY0, const int16_t *gradY1, int grad_stride,
                        int wgt_x, int wgt_y)
{
    __m128i x01 = _mm_loadl_epi64((__m128i *)&gradX0[0*grad_stride]);
    __m128i x02 = _mm_loadl_epi64((__m128i *)&gradX0[1*grad_stride]);
    __m128i x03 = _mm_loadl_epi64((__m128i *)&gradX0[2*grad_stride]);
    __m128i x04 = _mm_loadl_epi64((__m128i *)&gradX0[3*grad_stride]);

    __m128i x11 = _mm_loadl_epi64((__m128i *)&gradX1[0*grad_stride]);
    __m128i x12 = _mm_loadl_epi64((__m128i *)&gradX1[1*grad_stride]);
    __m128i x13 = _mm_loadl_epi64((__m128i *)&gradX1[2*grad_stride]);
    __m128i x14 = _mm_loadl_epi64((__m128i *)&gradX1[3*grad_stride]);

    __m128i y01 = _mm_loadl_epi64((__m128i *)&gradY0[0*grad_stride]);
    __m128i y02 = _mm_loadl_epi64((__m128i *)&gradY0[1*grad_stride]);
    __m128i y03 = _mm_loadl_epi64((__m128i *)&gradY0[2*grad_stride]);
    __m128i y04 = _mm_loadl_epi64((__m128i *)&gradY0[3*grad_stride]);

    __m128i y11 = _mm_loadl_epi64((__m128i *)&gradY1[0*grad_stride]);
    __m128i y12 = _mm_loadl_epi64((__m128i *)&gradY1[1*grad_stride]);
    __m128i y13 = _mm_loadl_epi64((__m128i *)&gradY1[2*grad_stride]);
    __m128i y14 = _mm_loadl_epi64((__m128i *)&gradY1[3*grad_stride]);

    __m128i src01 = _mm_loadl_epi64((__m128i *)&src0[0*src0_stride]);
    __m128i src02 = _mm_loadl_epi64((__m128i *)&src0[1*src0_stride]);
    __m128i src03 = _mm_loadl_epi64((__m128i *)&src0[2*src0_stride]);
    __m128i src04 = _mm_loadl_epi64((__m128i *)&src0[3*src0_stride]);

    __m128i src11 = _mm_loadl_epi64((__m128i *)&src1[0*src1_stride]);
    __m128i src12 = _mm_loadl_epi64((__m128i *)&src1[1*src1_stride]);
    __m128i src13 = _mm_loadl_epi64((__m128i *)&src1[2*src1_stride]);
    __m128i src14 = _mm_loadl_epi64((__m128i *)&src1[3*src1_stride]);

    __m128i wx = _mm_set1_epi16(wgt_x);
    __m128i wy = _mm_set1_epi16(wgt_y);
    __m128i offset = _mm_set1_epi32(BDOF_OFFSET);

    __m128i w = _mm_unpacklo_epi16(wx, wy);

      x01 = _mm_unpacklo_epi16(x01, y01);
      x02 = _mm_unpacklo_epi16(x02, y02);
      x03 = _mm_unpacklo_epi16(x03, y03);
      x04 = _mm_unpacklo_epi16(x04, y04);

      x11 = _mm_unpacklo_epi16(x11, y11);
      x12 = _mm_unpacklo_epi16(x12, y12);
      x13 = _mm_unpacklo_epi16(x13, y13);
      x14 = _mm_unpacklo_epi16(x14, y14);

      src01 = _mm_cvtepi16_epi32(src01);
      src02 = _mm_cvtepi16_epi32(src02);
      src03 = _mm_cvtepi16_epi32(src03);
      src04 = _mm_cvtepi16_epi32(src04);

      src11 = _mm_cvtepi16_epi32(src11);
      src12 = _mm_cvtepi16_epi32(src12);
      src13 = _mm_cvtepi16_epi32(src13);
      src14 = _mm_cvtepi16_epi32(src14);

      x01 = _mm_sub_epi16(x01, x11);
      x02 = _mm_sub_epi16(x02, x12);
      x03 = _mm_sub_epi16(x03, x13);
      x04 = _mm_sub_epi16(x04, x14);

      x01 = _mm_madd_epi16(x01, w);
      x02 = _mm_madd_epi16(x02, w);
      x03 = _mm_madd_epi16(x03, w);
      x04 = _mm_madd_epi16(x04, w);

      src01 = _mm_add_epi32(src01, src11);
      src02 = _mm_add_epi32(src02, src12);
      src03 = _mm_add_epi32(src03, src13);
      src04 = _mm_add_epi32(src04, src14);

      x01 = _mm_add_epi32(x01, offset);
      x02 = _mm_add_epi32(x02, offset);
      x03 = _mm_add_epi32(x03, offset);
      x04 = _mm_add_epi32(x04, offset);

      x01 = _mm_add_epi32(x01, src01);
      x02 = _mm_add_epi32(x02, src02);
      x03 = _mm_add_epi32(x03, src03);
      x04 = _mm_add_epi32(x04, src04);

      x01 = _mm_srai_epi32(x01, BDOF_SHIFT);
      x02 = _mm_srai_epi32(x02, BDOF_SHIFT);
      x03 = _mm_srai_epi32(x03, BDOF_SHIFT);
      x04 = _mm_srai_epi32(x04, BDOF_SHIFT);

      x01 = _mm_shufflelo_epi16(x01, 0x88);
      x02 = _mm_shufflelo_epi16(x02, 0x88);
      x03 = _mm_shufflelo_epi16(x03, 0x88);
      x04 = _mm_shufflelo_epi16(x04, 0x88);

      x01 = _mm_shufflehi_epi16(x01, 0x88);
      x02 = _mm_shufflehi_epi16(x02, 0x88);
      x03 = _mm_shufflehi_epi16(x03, 0x88);
      x04 = _mm_shufflehi_epi16(x04, 0x88);

      x01 = _mm_shuffle_epi32(x01, 0x88);
      x02 = _mm_shuffle_epi32(x02, 0x88);
      x03 = _mm_shuffle_epi32(x03, 0x88);
      x04 = _mm_shuffle_epi32(x04, 0x88);

      x01 = _mm_max_epi16(x01, _mm_setzero_si128());
      x02 = _mm_max_epi16(x02, _mm_setzero_si128());
      x03 = _mm_max_epi16(x03, _mm_setzero_si128());
      x04 = _mm_max_epi16(x04, _mm_setzero_si128());

      x01 = _mm_min_epi16(x01, _mm_set1_epi16(1023));
      x02 = _mm_min_epi16(x02, _mm_set1_epi16(1023));
      x03 = _mm_min_epi16(x03, _mm_set1_epi16(1023));
      x04 = _mm_min_epi16(x04, _mm_set1_epi16(1023));

      _mm_storel_epi64((__m128i *)&dst[0*dst_stride], x01);
      _mm_storel_epi64((__m128i *)&dst[1*dst_stride], x02);
      _mm_storel_epi64((__m128i *)&dst[2*dst_stride], x03);
      _mm_storel_epi64((__m128i *)&dst[3*dst_stride], x04);
}

void
rcn_init_prof_functions_sse(struct RCNFunctions *const rcn_funcs)
{
    rcn_funcs->prof.grad = &compute_prof_grad_sse;
    rcn_funcs->prof.rcn = &rcn_prof_sse;
    rcn_funcs->prof.tmp_prof_mrg = &tmp_prof_mrg_sse;
    rcn_funcs->prof.tmp_prof_mrg_w = &tmp_prof_mrg_w_sse;
}

void
rcn_init_bdof_functions_sse(struct RCNFunctions *const rcn_funcs)
{
    rcn_funcs->bdof.grad = &compute_prof_grad_sse;
    rcn_funcs->bdof.subblock = &rcn_apply_bdof_subblock_sse;
}
