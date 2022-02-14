#include <emmintrin.h>
#include <smmintrin.h>

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
#define PROF_PREC_RND (1 << (14 - 1))

#define BDOF_SHIFT   (14 + 1 - BITDEPTH)
#define BDOF_OFFSET  ((1 << (BDOF_SHIFT - 1)))
#define BDOF_WGT_LIMIT ((1 << 4) - 1)


static void rcn_prof_sse(OVSample* dst, int dst_stride, const int16_t* src, int src_stride,
         const int16_t* grad_x, const int16_t* grad_y, int grad_stride,
         const int16_t* dmv_scale_h, const int16_t* dmv_scale_v,
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

        #if 0
        sh1 = _mm_packs_epi32(sh1, _mm_setzero_si128());
        sh2 = _mm_packs_epi32(sh2, _mm_setzero_si128());
        sh3 = _mm_packs_epi32(sh3, _mm_setzero_si128());
        sh4 = _mm_packs_epi32(sh4, _mm_setzero_si128());

        sv1 = _mm_packs_epi32(sv1, _mm_setzero_si128());
        sv2 = _mm_packs_epi32(sv2, _mm_setzero_si128());
        sv3 = _mm_packs_epi32(sv3, _mm_setzero_si128());
        sv4 = _mm_packs_epi32(sv4, _mm_setzero_si128());
        #endif

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
        __m128i min_val = _mm_set1_epi16(-PROF_DELTA_LIMIT);
        __m128i max_val = _mm_set1_epi16(PROF_DELTA_LIMIT - 1);
        x1 = _mm_packs_epi32(x1, x2);
        x3 = _mm_packs_epi32(x3, x4);

        x1 = _mm_max_epi16(x1, min_val);
        x3 = _mm_max_epi16(x3, min_val);

        x1 = _mm_min_epi16(x1, max_val);
        x3 = _mm_min_epi16(x3, max_val);

        // val = (int16_t)src[x] -(1 << 13);
        // val += add;
        //
        // /* Clipping if not bi directional */
        // val = (val + 8200 /*+ PROF_SMP_OFFSET*/) >> PROF_SMP_SHIFT;
        // dst[x] = ov_clip(val, 0, 1023);

        __m128i offset = _mm_set1_epi16((1 << (13 - BITDEPTH)));

        //srcV1 = _mm_cvtepi16_epi32(srcV1);
        //srcV2 = _mm_cvtepi16_epi32(srcV2);
        //srcV3 = _mm_cvtepi16_epi32(srcV3);
        //srcV4 = _mm_cvtepi16_epi32(srcV4);
        srcV1 = _mm_unpacklo_epi64(srcV1, srcV2);
        srcV3 = _mm_unpacklo_epi64(srcV3, srcV4);

        #if 0
        srcV1 = _mm_add_epi16(srcV1, offset);
        srcV2 = _mm_add_epi16(srcV2, offset);
        srcV3 = _mm_add_epi16(srcV3, offset);
        srcV4 = _mm_add_epi16(srcV4, offset);
        #endif
        srcV1 = _mm_add_epi16(srcV1, offset);
        srcV3 = _mm_add_epi16(srcV3, offset);

        x1 = _mm_add_epi16(srcV1, x1);
        x3 = _mm_add_epi16(srcV3, x3);

        x1 = _mm_srai_epi16(x1, PROF_SMP_SHIFT);
        x3 = _mm_srai_epi16(x3, PROF_SMP_SHIFT);

        x1 = _mm_max_epi16(x1, _mm_setzero_si128());
        x3 = _mm_max_epi16(x3, _mm_setzero_si128());

        x1 = _mm_min_epi16(x1, _mm_set1_epi16(1023));
        x3 = _mm_min_epi16(x3, _mm_set1_epi16(1023));

        _mm_storel_epi64((__m128i *)&dst[0*dst_stride], x1);
        _mm_storel_epi64((__m128i *)&dst[1*dst_stride], _mm_bsrli_si128(x1, 8));
        _mm_storel_epi64((__m128i *)&dst[2*dst_stride], x3);
        _mm_storel_epi64((__m128i *)&dst[3*dst_stride], _mm_bsrli_si128(x3, 8));
    } else {
      __m128i min_val = _mm_set1_epi16(-PROF_DELTA_LIMIT);
      __m128i max_val = _mm_set1_epi16(PROF_DELTA_LIMIT - 1);

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

      #if 0
      sh1 = _mm_packs_epi32(sh1, _mm_setzero_si128());
      sh2 = _mm_packs_epi32(sh2, _mm_setzero_si128());
      sh3 = _mm_packs_epi32(sh3, _mm_setzero_si128());
      sh4 = _mm_packs_epi32(sh4, _mm_setzero_si128());

      sv1 = _mm_packs_epi32(sv1, _mm_setzero_si128());
      sv2 = _mm_packs_epi32(sv2, _mm_setzero_si128());
      sv3 = _mm_packs_epi32(sv3, _mm_setzero_si128());
      sv4 = _mm_packs_epi32(sv4, _mm_setzero_si128());
      #endif

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
      x1 = _mm_packs_epi32(x1, x2);
      x3 = _mm_packs_epi32(x3, x4);

      x1 = _mm_max_epi16(x1, min_val);
      x3 = _mm_max_epi16(x3, min_val);

      x1 = _mm_min_epi16(x1, max_val);
      x3 = _mm_min_epi16(x3, max_val);

      srcV1 = _mm_unpacklo_epi64(srcV1, srcV2);
      srcV3 = _mm_unpacklo_epi64(srcV3, srcV4);

      x1 = _mm_add_epi16(srcV1, x1);
      x3 = _mm_add_epi16(srcV3, x3);

      _mm_storel_epi64((__m128i *)&dst[0*dst_stride], x1);
      _mm_storel_epi64((__m128i *)&dst[1*dst_stride], _mm_bsrli_si128(x1, 8));
      _mm_storel_epi64((__m128i *)&dst[2*dst_stride], x3);
      _mm_storel_epi64((__m128i *)&dst[3*dst_stride], _mm_bsrli_si128(x3, 8));
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
derive_bdof_weights(const int16_t* ref0, const int16_t* ref1,
                    const int16_t* grad_x0, const int16_t* grad_x1,
                    const int16_t* grad_y0, const int16_t* grad_y1,
                    const int src0_stride, const int src1_stride,
                    const int grad_stride,
                    int *weight_x, int *weight_y)
{
    int wgt_x = 0;
    int wgt_y = 0;
    int16_t var[8];

    int i, j;
    __m128i rnd = _mm_set1_epi16(PROF_PREC_RND);

    __m128i z = _mm_setzero_si128();
    __m128i accu_sgn_y_avg_x = z;
    __m128i accu_abs_x = z;
    __m128i accu_abs_y = z;
    __m128i accu_delta_x = z;
    __m128i accu_delta_y = z;

    for (i = 0; i < SB_H + 2; i++) {
        __m128i g_x0 = _mm_loadu_si128((__m128i *)grad_x0);
        __m128i g_x1 = _mm_loadu_si128((__m128i *)grad_x1);
        __m128i g_y0 = _mm_loadu_si128((__m128i *)grad_y0);
        __m128i g_y1 = _mm_loadu_si128((__m128i *)grad_y1);

        __m128i r0 = _mm_loadu_si128((__m128i *)ref0);
        __m128i r1 = _mm_loadu_si128((__m128i *)ref1);

        __m128i avg_g_x = _mm_add_epi16(g_x0, g_x1);
        __m128i avg_g_y = _mm_add_epi16(g_y0, g_y1);

        avg_g_x = _mm_srai_epi16(avg_g_x, 1);
        avg_g_y = _mm_srai_epi16(avg_g_y, 1);

        accu_abs_x = _mm_add_epi16(accu_abs_x, _mm_abs_epi16(avg_g_x));
        accu_abs_y = _mm_add_epi16(accu_abs_y, _mm_abs_epi16(avg_g_y));

        r0 = _mm_sub_epi16(r0, rnd);
        r1 = _mm_sub_epi16(r1, rnd);

        r0 = _mm_srai_epi16(r0, 4);
        r1 = _mm_srai_epi16(r1, 4);

        __m128i delta_ref = _mm_sub_epi16(r1, r0);

        __m128i sign_msk_y = _mm_cmplt_epi16(avg_g_y, z);
        __m128i z_msk_y    = _mm_cmpeq_epi16(avg_g_y, z);
        __m128i s_delta =_mm_xor_si128(delta_ref, sign_msk_y);

        s_delta = _mm_sub_epi16(s_delta, sign_msk_y);
        s_delta = _mm_andnot_si128(z_msk_y, s_delta);

        accu_delta_y = _mm_add_epi16(accu_delta_y, s_delta);

        __m128i sign_y_avg_x = _mm_xor_si128(avg_g_x, sign_msk_y); 

        sign_y_avg_x = _mm_sub_epi16(sign_y_avg_x, sign_msk_y);
        sign_y_avg_x = _mm_andnot_si128(z_msk_y, sign_y_avg_x);

        accu_sgn_y_avg_x = _mm_add_epi16(accu_sgn_y_avg_x, sign_y_avg_x);

        __m128i sign_msk_x = _mm_cmplt_epi16(avg_g_x, z);
        __m128i z_msk_x    = _mm_cmpeq_epi16(avg_g_x, z);

        s_delta =_mm_xor_si128(delta_ref, sign_msk_x);
        s_delta = _mm_sub_epi16(s_delta, sign_msk_x);
        s_delta = _mm_andnot_si128(z_msk_x, s_delta);

        accu_delta_x = _mm_add_epi16(accu_delta_x, s_delta);


        ref1 += src1_stride;
        ref0 += src0_stride;

        grad_x0 += grad_stride;
        grad_x1 += grad_stride;

        grad_y0 += grad_stride;
        grad_y1 += grad_stride;
    }

    __m128i msk = _mm_set1_epi16(-1);

    msk = _mm_bsrli_si128(msk, 4);


    accu_abs_x = _mm_and_si128(accu_abs_x, msk);
    accu_abs_y = _mm_and_si128(accu_abs_y, msk);
    accu_delta_x = _mm_and_si128(accu_delta_x, msk);
    accu_delta_y = _mm_and_si128(accu_delta_y, msk);
    accu_sgn_y_avg_x = _mm_and_si128(accu_sgn_y_avg_x, msk);

    __m128i sum_xy    = _mm_hadd_epi16(accu_abs_x, accu_abs_y);
    __m128i sum_delta = _mm_hadd_epi16(accu_delta_x, accu_delta_y);
    __m128i sum_sign  = _mm_hadd_epi16(accu_sgn_y_avg_x, z);

    sum_xy = _mm_hadd_epi16(sum_xy, sum_delta);
    sum_sign = _mm_hadd_epi16(sum_sign, z);
    sum_xy = _mm_hadd_epi16(sum_xy, sum_sign);

    _mm_storeu_si128((__m128i*)&var[0], sum_xy);

    int16_t sum_abs_x = var[0];
    int16_t sum_abs_y = var[1];
    int16_t sum_delta_x = var[2];
    int16_t sum_delta_y = var[3];
    int16_t sum_avg_x_y_signs = var[4];

    if (sum_abs_x) {
        int log2_renorm_x = floor_log2(sum_abs_x);

        wgt_x = (sum_delta_x << 2) >> log2_renorm_x;
        wgt_x = ov_clip(wgt_x, -BDOF_WGT_LIMIT, BDOF_WGT_LIMIT);
        *weight_x = wgt_x;
    }

    if (sum_abs_y) {
        int log2_renorm_y = floor_log2(sum_abs_y);
        int x_offset = 0;

        if (wgt_x) {
            /* FIXME understand this part */
            int high = sum_avg_x_y_signs >> 12;
            int low  = sum_avg_x_y_signs & ((1 << 12) - 1);
            x_offset = (((wgt_x * high) << 12) + (wgt_x * low)) >> 1;
        }

        wgt_y = ((sum_delta_y << 2) - x_offset) >> log2_renorm_y;
        wgt_y = ov_clip(wgt_y, -BDOF_WGT_LIMIT, BDOF_WGT_LIMIT);
        *weight_y = wgt_y;
    }
}

static void
rcn_bdof(struct BDOFFunctions *const bdof, OVSample *dst, int dst_stride,
         const int16_t *ref_bdof0, const int16_t *ref_bdof1, int ref_stride,
         const int16_t *grad_x0, const int16_t *grad_y0,
         const int16_t *grad_x1, const int16_t *grad_y1,
         int grad_stride, uint8_t pb_w, uint8_t pb_h)
{
    int nb_sb_w = (pb_w >> 2);
    int nb_sb_h = (pb_h >> 2);

    const int16_t *grad_x0_ln = grad_x0;
    const int16_t *grad_y0_ln = grad_y0;
    const int16_t *grad_x1_ln = grad_x1;
    const int16_t *grad_y1_ln = grad_y1;

    const int16_t *ref0_ln = ref_bdof0 - 128 - 1;
    const int16_t *ref1_ln = ref_bdof1 - 128 - 1;
    OVSample *dst_ln = dst;

    int i, j;

    for (i = 0; i < nb_sb_h; i++) {
        const int16_t *ref0_tmp = ref0_ln;
        const int16_t *ref1_tmp = ref1_ln;

        grad_x0 = grad_x0_ln;
        grad_y0 = grad_y0_ln;
        grad_x1 = grad_x1_ln;
        grad_y1 = grad_y1_ln;

        dst = dst_ln;

        for (j = 0; j < nb_sb_w; j++) {
            int wgt_x = 0;
            int wgt_y = 0;

            derive_bdof_weights(ref0_tmp, ref1_tmp,
                                grad_x0, grad_x1, grad_y0, grad_y1,
                                ref_stride, ref_stride,
                                grad_stride,
                                &wgt_x, &wgt_y);

            bdof->subblock(ref0_tmp + ref_stride + 1, ref_stride,
                           ref1_tmp + ref_stride + 1, ref_stride,
                           dst, dst_stride,
                           grad_x0 + grad_stride + 1, grad_x1 + grad_stride + 1,
                           grad_y0 + grad_stride + 1, grad_y1 + grad_stride + 1,
                           grad_stride,
                           wgt_x, wgt_y);

            grad_x0 += 1 << 2;
            grad_x1 += 1 << 2;
            grad_y0 += 1 << 2;
            grad_y1 += 1 << 2;
            ref0_tmp += 1 << 2;
            ref1_tmp += 1 << 2;
            dst    += 1 << 2;
        }

        grad_x0_ln += grad_stride << 2;
        grad_y0_ln += grad_stride << 2;
        grad_x1_ln += grad_stride << 2;
        grad_y1_ln += grad_stride << 2;

        ref0_ln += ref_stride << 2;
        ref1_ln += ref_stride << 2;

        dst_ln += dst_stride << 2;
    }
}

static void
rcn_apply_bdof_subblock_sse(const int16_t* src0, int src0_stride,
                            const int16_t* src1, int src1_stride,
                            OVSample *dst, int dst_stride,
                            const int16_t *gradX0, const int16_t *gradX1,
                            const int16_t *gradY0, const int16_t *gradY1, int grad_stride,
                            int wgt_x, int wgt_y)
{
    __m128i x01 = _mm_loadu_si128((__m128i *)&gradX0[0*grad_stride]);
    __m128i x02 = _mm_loadu_si128((__m128i *)&gradX0[1*grad_stride]);
    __m128i x03 = _mm_loadu_si128((__m128i *)&gradX0[2*grad_stride]);
    __m128i x04 = _mm_loadu_si128((__m128i *)&gradX0[3*grad_stride]);

    __m128i x11 = _mm_loadu_si128((__m128i *)&gradX1[0*grad_stride]);
    __m128i x12 = _mm_loadu_si128((__m128i *)&gradX1[1*grad_stride]);
    __m128i x13 = _mm_loadu_si128((__m128i *)&gradX1[2*grad_stride]);
    __m128i x14 = _mm_loadu_si128((__m128i *)&gradX1[3*grad_stride]);

    __m128i y01 = _mm_loadu_si128((__m128i *)&gradY0[0*grad_stride]);
    __m128i y02 = _mm_loadu_si128((__m128i *)&gradY0[1*grad_stride]);
    __m128i y03 = _mm_loadu_si128((__m128i *)&gradY0[2*grad_stride]);
    __m128i y04 = _mm_loadu_si128((__m128i *)&gradY0[3*grad_stride]);

    __m128i y11 = _mm_loadu_si128((__m128i *)&gradY1[0*grad_stride]);
    __m128i y12 = _mm_loadu_si128((__m128i *)&gradY1[1*grad_stride]);
    __m128i y13 = _mm_loadu_si128((__m128i *)&gradY1[2*grad_stride]);
    __m128i y14 = _mm_loadu_si128((__m128i *)&gradY1[3*grad_stride]);

    __m128i src01 = _mm_loadu_si128((__m128i *)&src0[0*src0_stride]);
    __m128i src02 = _mm_loadu_si128((__m128i *)&src0[1*src0_stride]);
    __m128i src03 = _mm_loadu_si128((__m128i *)&src0[2*src0_stride]);
    __m128i src04 = _mm_loadu_si128((__m128i *)&src0[3*src0_stride]);

    __m128i src11 = _mm_loadu_si128((__m128i *)&src1[0*src1_stride]);
    __m128i src12 = _mm_loadu_si128((__m128i *)&src1[1*src1_stride]);
    __m128i src13 = _mm_loadu_si128((__m128i *)&src1[2*src1_stride]);
    __m128i src14 = _mm_loadu_si128((__m128i *)&src1[3*src1_stride]);

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

      x01 = _mm_packs_epi32(x01, x02);
      x03 = _mm_packs_epi32(x03, x04);

      x01 = _mm_max_epi16(x01, _mm_setzero_si128());
      x03 = _mm_max_epi16(x03, _mm_setzero_si128());

      x01 = _mm_min_epi16(x01, _mm_set1_epi16(1023));
      x03 = _mm_min_epi16(x03, _mm_set1_epi16(1023));

      _mm_storel_epi64((__m128i *)&dst[0*dst_stride], x01);
      _mm_storel_epi64((__m128i *)&dst[1*dst_stride], _mm_bsrli_si128(x01, 8));
      _mm_storel_epi64((__m128i *)&dst[2*dst_stride], x03);
      _mm_storel_epi64((__m128i *)&dst[3*dst_stride], _mm_bsrli_si128(x03, 8));
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
    rcn_funcs->bdof.rcn_bdof = &rcn_bdof;
}
