#include <stdint.h>

#include "rcn_structures.h"
#include "ovutils.h"
#include "ctudec.h"

#include "bitdepth.h"

#define SWAP(type,a,b) do{type tmp= b; b = a; a = tmp;}while(0)


struct AVGMinMax{
    uint16_t min_l;
    uint16_t max_l;
    uint16_t min_cb;
    uint16_t max_cb;
    uint16_t min_cr;
    uint16_t max_cr;
};

static void
sub_sample_lm_ref_lft_collocated(const uint16_t *lm_src, const uint16_t *src_cb, const uint16_t *src_cr,
                      uint16_t *const lm_dst, uint16_t *dst_cb, uint16_t *dst_cr,
                      ptrdiff_t lm_src_stride, ptrdiff_t src_c_stride,
                      int lft_step, int nb_sample_lft, int abv_avail)
{
    int start_pos = lft_step >> 1;
    int lm_src_stride2 = lm_src_stride << 1;
    int stride_l = lm_src_stride2 * lft_step;
    int stride_c = src_c_stride * lft_step;
    const uint16_t *_src    = lm_src - 2 + start_pos * lm_src_stride2;
    const uint16_t *_src_cb = src_cb - 1 + start_pos * src_c_stride;
    const uint16_t *_src_cr = src_cr - 1 + start_pos * src_c_stride;
    int i;

    uint8_t padd_abv = start_pos == 0 && !abv_avail;
    for (i = 0; i < nb_sample_lft; i++) {
        int s = 4;

        s += _src[-(padd_abv ? 0 : lm_src_stride)];
        s += _src[0] * 4;
        s += _src[-1];
        s += _src[1];
        s += _src[lm_src_stride];

        lm_dst[i] = s >> 3;

        dst_cb[i] = _src_cb[0];
        dst_cr[i] = _src_cr[0];
        padd_abv = 0;

        _src    += stride_l;
        _src_cb += stride_c;
        _src_cr += stride_c;
    }
}

static void
sub_sample_lm_ref_abv_collocated(const uint16_t *lm_src, const uint16_t *src_cb, const uint16_t *src_cr,
                      uint16_t *const lm_dst, uint16_t *dst_cb, uint16_t *dst_cr,
                      ptrdiff_t lm_src_stride, ptrdiff_t src_c_stride,
                      int abv_step, int nb_sample_abv, uint8_t lft_avail)
{
    int start_pos = abv_step >> 1;
    const uint16_t *_src = lm_src - (lm_src_stride << 1) + (start_pos << 1);
    const uint16_t *_src_cb = src_cb - src_c_stride + start_pos;
    const uint16_t *_src_cr = src_cr - src_c_stride + start_pos;
    int i;

    int abv_step_l = abv_step << 1;

    uint8_t pad_left = start_pos == 0 && !lft_avail;

    for (i = 0; i < nb_sample_abv; i++) {
        int s = 4;

        s += _src[0 - lm_src_stride];
        s += _src[0] * 4;
        s += _src[0 - (!pad_left)];
        s += _src[0 + 1];
        s += _src[0 + lm_src_stride];

        lm_dst[i] = s >> 3;
        dst_cb[i] = _src_cb[0];
        dst_cr[i] = _src_cr[0];

        pad_left = 0;

        _src    += abv_step_l;
        _src_cb += abv_step;
        _src_cr += abv_step;
    }
}

static void
compute_lm_subsample_collocated(const uint16_t *lm_src, uint16_t *dst_cb, uint16_t *dst_cr,
                     ptrdiff_t lm_src_stride, ptrdiff_t dst_stride_c,
                     const struct CCLMParams *const lm_params,
                     int pb_w, int pb_h, uint8_t lft_avail, uint8_t abv_avail)
{
    int i, j;
    ptrdiff_t lm_src_stride2 = lm_src_stride << 1;

    const int scale_cb  = lm_params->cb.a;
    const int offset_cb = lm_params->cb.b;
    const int shift_cb  = lm_params->cb.shift;
    const int scale_cr  = lm_params->cr.a;
    const int offset_cr = lm_params->cr.b;
    const int shift_cr  = lm_params->cr.shift;

    for (j = 0; j < pb_h; j++) {
        const uint8_t padd_abv = j == 0 && !abv_avail;
        for (i = 0; i < pb_w; i++) {
            const uint8_t pad_left = i == 0 && !lft_avail;

            int lm_val = 4;
            lm_val += lm_src[2 * i - (padd_abv ? 0 : lm_src_stride)];
            lm_val += lm_src[2 * i] * 4;
            lm_val += lm_src[2 * i - (!pad_left)];
            lm_val += lm_src[2 * i + 1];
            lm_val += lm_src[2 * i + lm_src_stride];

            lm_val >>= 3;
            dst_cb[i] = ov_bdclip(((lm_val * scale_cb) >> shift_cb) + offset_cb);
            dst_cr[i] = ov_bdclip(((lm_val * scale_cr) >> shift_cr) + offset_cr);
        }
        dst_cb += dst_stride_c;
        dst_cr += dst_stride_c;
        lm_src += lm_src_stride2;
    }
}

static inline struct LMParams
compute_lm_params(int16_t avg_min_l, int16_t avg_min_c, int16_t avg_max_c, int16_t v, int8_t log2_rng_l)
{
    struct LMParams lm_params;
    int a, b, shift;
    int range_c = avg_max_c - avg_min_c;

    /* Note this is not max nor min of chroma samples, but
     * the corresponding samples to max and min luma samples
     * so we require an absolute value
     */

    int log2_rng_c_plus1 = range_c ? floor_log2(OVABS(range_c)) + 1 : 0;
    int add = (1 << log2_rng_c_plus1) >> 1;

    a = (range_c * v + add) >> log2_rng_c_plus1;

    shift = 3 + log2_rng_l - log2_rng_c_plus1;

    /* FIXME find branchless computation of shift and scale values */
    if (shift < 1) {
        shift = 1;
        a = -(!!a) & ((a < 0) ? -15 : 15);
    }

    b = avg_min_c - ((a * avg_min_l) >> shift);

    lm_params.shift = shift;
    lm_params.a = a;
    lm_params.b = b;

    return lm_params;
}

static inline struct CCLMParams
derive_cclm_params(const struct AVGMinMax *const avgs)
{
    struct CCLMParams lm_params = {
        .cb = {.a = 0, .b = avgs->min_cb, .shift = 0},
        .cr = {.a = 0, .b = avgs->min_cr, .shift = 0}
    };

    int range_l = avgs->max_l - avgs->min_l;

    if (range_l) {
        static const uint8_t div_lut[1 << 4] = {
            0,  7,  6,  5,  5,  4,  4,  3,  3,  2,  2,  1,  1,  1,  1,  0
        };

        int log2_rng_l = floor_log2(range_l);

        int norm_diff = ((range_l << 4) >> log2_rng_l) & 0xF;

        /*FIXME | 8 could be added to LUT */
        int v = div_lut[norm_diff] | 8;

        log2_rng_l += norm_diff != 0;

        lm_params.cb = compute_lm_params(avgs->min_l, avgs->min_cb, avgs->max_cb, v, log2_rng_l);
        lm_params.cr = compute_lm_params(avgs->min_l, avgs->min_cr, avgs->max_cr, v, log2_rng_l);
    }

    return lm_params;
}

static inline struct AVGMinMax
sort_average_lm_ref_samples(const uint16_t *const lm_smp, const uint16_t *const smp_cb,
                            const uint16_t *const smp_cr, int nb_samples)
{
    struct AVGMinMax avg_min_max;

    if (nb_samples == 2) {
        int min_idx = lm_smp[0] >= lm_smp[1];
        int max_idx = !min_idx;

        avg_min_max.min_l = lm_smp[min_idx];
        avg_min_max.max_l = lm_smp[max_idx];

        avg_min_max.min_cb = smp_cb[min_idx];
        avg_min_max.max_cb = smp_cb[max_idx];

        avg_min_max.min_cr = smp_cr[min_idx];
        avg_min_max.max_cr = smp_cr[max_idx];

    } else {
        /* FIXME better way of sorting LUT */

        int8_t idx[4] = { 0, 2, 1, 3 };

        int8_t *min_idx = &idx[0];
        int8_t *max_idx = &idx[2];

        if (lm_smp[0] > lm_smp[2]) SWAP(int8_t, min_idx[0], min_idx[1]);
        if (lm_smp[1] > lm_smp[3]) SWAP(int8_t, max_idx[0], max_idx[1]);

        if (lm_smp[min_idx[0]] > lm_smp[max_idx[1]]) SWAP(int8_t*, min_idx, max_idx);
        if (lm_smp[min_idx[1]] > lm_smp[max_idx[0]]) SWAP(int8_t, min_idx[1], max_idx[0]);

        avg_min_max.min_l = (lm_smp[min_idx[0]] + lm_smp[min_idx[1]] + 1) >> 1;
        avg_min_max.max_l = (lm_smp[max_idx[0]] + lm_smp[max_idx[1]] + 1) >> 1;

        avg_min_max.min_cb = (smp_cb[min_idx[0]] + smp_cb[min_idx[1]] + 1) >> 1;
        avg_min_max.max_cb = (smp_cb[max_idx[0]] + smp_cb[max_idx[1]] + 1) >> 1;

        avg_min_max.min_cr = (smp_cr[min_idx[0]] + smp_cr[min_idx[1]] + 1) >> 1;
        avg_min_max.max_cr = (smp_cr[max_idx[0]] + smp_cr[max_idx[1]] + 1) >> 1;
    }

    return avg_min_max;
}

static void
sub_sample_lm_ref_lft(const uint16_t *lm_src, const uint16_t *src_cb, const uint16_t *src_cr,
                      uint16_t *const lm_dst, uint16_t *dst_cb, uint16_t *dst_cr,
                      ptrdiff_t lm_src_stride, ptrdiff_t src_c_stride,
                      int lft_step, int nb_sample_lft)
{
    int start_pos = lft_step >> 1;
    int lm_src_stride2 = lm_src_stride << 1;
    int stride_l = lm_src_stride2 * lft_step;
    int stride_c = src_c_stride * lft_step;
    const uint16_t *_src    = lm_src - 2 + start_pos * lm_src_stride2;
    const uint16_t *_src_cb = src_cb - 1 + start_pos * src_c_stride;
    const uint16_t *_src_cr = src_cr - 1 + start_pos * src_c_stride;
    int i;

    for (i = 0; i < nb_sample_lft; i++) {
        int s = 4;
        s += _src[ 0] * 2;
        s += _src[ 1];
        s += _src[-1];
        s += _src[lm_src_stride] * 2;
        s += _src[lm_src_stride + 1];
        s += _src[lm_src_stride - 1];

        lm_dst[i] = s >> 3;

        dst_cb[i] = _src_cb[0];
        dst_cr[i] = _src_cr[0];

        _src    += stride_l;
        _src_cb += stride_c;
        _src_cr += stride_c;
    }
}

static void
sub_sample_lm_ref_abv0(const uint16_t *lm_src, const uint16_t *src_cb, const uint16_t *src_cr,
                       uint16_t *const lm_dst, uint16_t *dst_cb, uint16_t *dst_cr,
                       ptrdiff_t lm_src_stride, ptrdiff_t src_c_stride,
                       int abv_step, int nb_sample_abv, uint8_t lft_avail)
{
    int start_pos = abv_step >> 1;
    const uint16_t *_src = lm_src - lm_src_stride + (start_pos << 1);
    const uint16_t *_src_cb = src_cb - src_c_stride + start_pos;
    const uint16_t *_src_cr = src_cr - src_c_stride + start_pos;
    int i;

    int abv_step_l = abv_step << 1;
    uint8_t pad_left = start_pos == 0 && !lft_avail;

    for (i = 0; i < nb_sample_abv; i++) {

        int s = 2;
        s += _src[0              ] * 2;
        s += _src[0 - (!pad_left)];
        s += _src[0 + 1          ];

        lm_dst[i] = s >> 2;
        dst_cb[i] = _src_cb[0];
        dst_cr[i] = _src_cr[0];

        _src += abv_step_l;
        _src_cb += abv_step;
        _src_cr += abv_step;
        pad_left = 0;
    }
}

static void
sub_sample_lm_ref_abv(const uint16_t *lm_src, const uint16_t *src_cb, const uint16_t *src_cr,
                      uint16_t *const lm_dst, uint16_t *dst_cb, uint16_t *dst_cr,
                      ptrdiff_t lm_src_stride, ptrdiff_t src_c_stride,
                      int abv_step, int nb_sample_abv, uint8_t lft_avail)
{

    int start_pos = abv_step >> 1;
    const uint16_t *_src    = lm_src - (lm_src_stride << 1) + (start_pos << 1);
    const uint16_t *_src_cb = src_cb - src_c_stride + start_pos;
    const uint16_t *_src_cr = src_cr - src_c_stride + start_pos;
    int i;

    int abv_step_l = abv_step << 1;

    uint8_t pad_left = start_pos == 0 && !lft_avail;

    for (i = 0; i < nb_sample_abv; i++) {
        int s = 4;
        s += _src[0] * 2;
        s += _src[0 + 1];
        s += _src[0 - (!pad_left)];
        s += _src[0 + lm_src_stride] * 2;
        s += _src[0 + 1 + lm_src_stride];
        s += _src[0 + lm_src_stride - (!pad_left)];

        lm_dst[i] = s >> 3;
        dst_cb[i] = _src_cb[0];
        dst_cr[i] = _src_cr[0];

        pad_left = 0;

        _src += abv_step_l;
        _src_cb += abv_step;
        _src_cr += abv_step;
    }
}

static void
compute_lm_subsample(const uint16_t *lm_src, uint16_t *dst_cb, uint16_t *dst_cr,
                     ptrdiff_t lm_src_stride, ptrdiff_t dst_stride_c,
                     const struct CCLMParams *const lm_params,
                     int pb_w, int pb_h, uint8_t lft_avail)
{
    int i, j;
    ptrdiff_t lm_src_stride2 = lm_src_stride << 1;

    const int scale_cb  = lm_params->cb.a;
    const int offset_cb = lm_params->cb.b;
    const int shift_cb  = lm_params->cb.shift;
    const int scale_cr  = lm_params->cr.a;
    const int offset_cr = lm_params->cr.b;
    const int shift_cr  = lm_params->cr.shift;

    for (j = 0; j < pb_h; j++) {
        for (i = 0; i < pb_w; i++) {
            const uint8_t pad_left = i == 0 && !lft_avail;

            int lm_val = 4;
            lm_val += lm_src[2 * i + 1];
            lm_val += lm_src[2 * i - (!pad_left)];
            lm_val += lm_src[2 * i] * 2;
            lm_val += lm_src[2 * i + lm_src_stride] * 2;
            lm_val += lm_src[2 * i + 1 + lm_src_stride];
            lm_val += lm_src[2 * i + lm_src_stride - (!pad_left)];
            lm_val >>= 3;

            dst_cb[i] = ov_bdclip(((lm_val * scale_cb) >> shift_cb) + offset_cb);
            dst_cr[i] = ov_bdclip(((lm_val * scale_cr) >> shift_cr) + offset_cr);
        }
        dst_cb += dst_stride_c;
        dst_cr += dst_stride_c;
        lm_src += lm_src_stride2;
    }
}

static void
vvc_intra_cclm2(const uint16_t *lm_src, uint16_t *dst_cb, uint16_t *dst_cr,
               ptrdiff_t lm_src_stride, ptrdiff_t dst_stride_c,
               int log2_pb_w, int log2_pb_h, uint8_t lft_avail, uint8_t abv_avail,
               uint8_t y0, LMsubsampleFunc const compute_subsample)
{
    struct CCLMParams lm_params = {
        .cb = {.a = 0, .b = AVG_VAL, .shift = 0},
        .cr = {.a = 0, .b = AVG_VAL, .shift = 0}
    };

    int pb_w = 1 << log2_pb_w;
    int pb_h = 1 << log2_pb_h;

    if (abv_avail || lft_avail){
        struct AVGMinMax avg_min_max;
        uint16_t lm_smp[4];
        uint16_t smp_cb[4];
        uint16_t smp_cr[4];
        uint8_t log2_nb_smp_abv = !!abv_avail + !lft_avail;
        uint8_t log2_nb_smp_lft = !!lft_avail + !abv_avail;

        int nb_sample_abv = (abv_avail + !lft_avail) << 1;
        int nb_sample_lft = (lft_avail + !abv_avail) << 1;
        int nb_sample;

        if (abv_avail) {
            int abv_step = OVMAX(1, (pb_w >> log2_nb_smp_abv));
            uint8_t ctu_first_line = !y0;

            /* in case of abv only ref_length might be 2 while nb_sample_lft is 4
               We are forced to reduce nb_smp in this particular case*/
            nb_sample_abv = OVMIN(pb_w, nb_sample_abv);

            /*FIXME avoid checking for first_line */
            if (ctu_first_line) {
                sub_sample_lm_ref_abv0(lm_src, dst_cb, dst_cr, lm_smp, smp_cb, smp_cr,
                                       lm_src_stride, dst_stride_c,
                                       abv_step, nb_sample_abv, lft_avail);
            } else {
                sub_sample_lm_ref_abv(lm_src, dst_cb, dst_cr, lm_smp, smp_cb, smp_cr,
                                      lm_src_stride, dst_stride_c,
                                      abv_step, nb_sample_abv, lft_avail);
            }
        }

        if (lft_avail) {
            int lft_step = OVMAX(1, (pb_h >> log2_nb_smp_lft));

            /* in case of left only ref_length might be 2 while nb_sample_lft is 4
               We are forced to reduce nb_smp in this particular case */
            nb_sample_lft = OVMIN(pb_h, nb_sample_lft);

            sub_sample_lm_ref_lft(lm_src, dst_cb, dst_cr, &lm_smp[nb_sample_abv],
                                  &smp_cb[nb_sample_abv], &smp_cr[nb_sample_abv],
                                  lm_src_stride, dst_stride_c,
                                  lft_step, nb_sample_lft);
        }

        nb_sample = nb_sample_lft + nb_sample_abv;

        avg_min_max = sort_average_lm_ref_samples(lm_smp, smp_cb, smp_cr, nb_sample);

        lm_params = derive_cclm_params(&avg_min_max);
    }

    compute_subsample(lm_src, dst_cb, dst_cr, lm_src_stride, dst_stride_c,
                       &lm_params, pb_w, pb_h, lft_avail);
}

static void
vvc_intra_cclm(const uint16_t *const src_luma, uint16_t *const dst_cb,
               uint16_t *const dst_cr, int log2_pb_w, int log2_pb_h,
               int y0, int up_available, int left_available, LMsubsampleFunc const compute_subsample)
{
    vvc_intra_cclm2(src_luma, dst_cb, dst_cr, RCN_CTB_STRIDE, RCN_CTB_STRIDE,
               log2_pb_w, log2_pb_h, left_available, up_available, y0, compute_subsample);
}

static void
vvc_intra_cclm2_collocated(const uint16_t *lm_src, uint16_t *dst_cb, uint16_t *dst_cr,
                           ptrdiff_t lm_src_stride, ptrdiff_t dst_stride_c,
                           int log2_pb_w, int log2_pb_h, uint8_t lft_avail, uint8_t abv_avail,
                           uint8_t y0)
{
    struct CCLMParams lm_params = {
        .cb = {.a = 0, .b = AVG_VAL, .shift = 0},
        .cr = {.a = 0, .b = AVG_VAL, .shift = 0}
    };

    int pb_w = 1 << log2_pb_w;
    int pb_h = 1 << log2_pb_h;

    if (abv_avail || lft_avail){
        struct AVGMinMax avg_min_max;
        uint16_t lm_smp[4];
        uint16_t smp_cb[4];
        uint16_t smp_cr[4];
        uint8_t log2_nb_smp_abv = !!abv_avail + !lft_avail;
        uint8_t log2_nb_smp_lft = !!lft_avail + !abv_avail;

        int nb_sample_abv = (abv_avail + !lft_avail) << 1;
        int nb_sample_lft = (lft_avail + !abv_avail) << 1;
        int nb_sample;

        if (abv_avail) {
            int abv_step = OVMAX(1, (pb_w >> log2_nb_smp_abv));
            uint8_t ctu_first_line = !y0;

            /* in case of abv only ref_length might be 2 while nb_sample_lft is 4
               We are forced to reduce nb_smp in this particular case*/
            nb_sample_abv = OVMIN(pb_w, nb_sample_abv);

            /*FIXME avoid checking for first_line */
            if (ctu_first_line) {
                sub_sample_lm_ref_abv0(lm_src, dst_cb, dst_cr, lm_smp, smp_cb, smp_cr,
                                       lm_src_stride, dst_stride_c,
                                       abv_step, nb_sample_abv, lft_avail);
            } else {
                sub_sample_lm_ref_abv_collocated(lm_src, dst_cb, dst_cr, lm_smp, smp_cb, smp_cr,
                                                 lm_src_stride, dst_stride_c,
                                                 abv_step, nb_sample_abv, lft_avail);
            }
        }

        if (lft_avail) {
            int lft_step = OVMAX(1, (pb_h >> log2_nb_smp_lft));

            /* in case of left only ref_length might be 2 while nb_sample_lft is 4
               We are forced to reduce nb_smp in this particular case */
            nb_sample_lft = OVMIN(pb_h, nb_sample_lft);

            sub_sample_lm_ref_lft_collocated(lm_src, dst_cb, dst_cr, &lm_smp[nb_sample_abv],
                                             &smp_cb[nb_sample_abv], &smp_cr[nb_sample_abv],
                                             lm_src_stride, dst_stride_c,
                                             lft_step, nb_sample_lft, abv_avail);
        }

        nb_sample = nb_sample_lft + nb_sample_abv;

        avg_min_max = sort_average_lm_ref_samples(lm_smp, smp_cb, smp_cr, nb_sample);

        lm_params = derive_cclm_params(&avg_min_max);
    }

    compute_lm_subsample_collocated(lm_src, dst_cb, dst_cr, lm_src_stride, dst_stride_c,
                                    &lm_params, pb_w, pb_h, lft_avail, abv_avail);
}

static void
vvc_intra_cclm_cl(const uint16_t *const src_luma, uint16_t *const dst_cb,
                  uint16_t *const dst_cr, int log2_pb_w, int log2_pb_h,
                  int y0, int up_available, int left_available, LMsubsampleFunc const compute_subsample)
{
    vvc_intra_cclm2_collocated(src_luma, dst_cb, dst_cr, RCN_CTB_STRIDE, RCN_CTB_STRIDE,
               log2_pb_w, log2_pb_h, left_available, up_available, y0);
}


static void
vvc_intra_mdlm_t2(const uint16_t *lm_src, uint16_t *dst_cb, uint16_t *dst_cr,
                  ptrdiff_t lm_src_stride, ptrdiff_t dst_stride_c,
                  int log2_pb_w, int log2_pb_h, uint8_t lft_avail, uint8_t abv_avail,
                  uint8_t y0, uint8_t x0, uint64_t abv_map, LMsubsampleFunc const compute_subsample)
{
    struct CCLMParams lm_params = {
        .cb = {.a = 0, .b = AVG_VAL, .shift = 0},
        .cr = {.a = 0, .b = AVG_VAL, .shift = 0}
    };

    int pb_w = 1 << log2_pb_w;
    int pb_h = 1 << log2_pb_h;

    if (abv_avail) {
        struct AVGMinMax avg_min_max;
        uint16_t lm_smp[4];
        uint16_t smp_cb[4];
        uint16_t smp_cr[4];

        int ref_length_abv = pb_w + OVMIN(pb_h, pb_w);
        int nb_pb_ref_abv = ref_length_abv >> 1;

        int x_pb = x0 >> 1;

        uint64_t abv_msk = (1llu << nb_pb_ref_abv) - 1;
        uint64_t usable_mask = ((abv_map >> (x_pb + 1)) & abv_msk);
        int nb_avail_ref_pb_abv = __builtin_ctzll(~usable_mask);

        int avail_ref_length_abv = nb_avail_ref_pb_abv << 1;

        int nb_sample_abv = OVMIN(avail_ref_length_abv, 4);
        int abv_step = OVMAX(1, (avail_ref_length_abv >> 2));
        uint8_t ctu_first_line = !y0;

        /* in case of abv only ref_length might be 2 while nb_sample_lft is 4
           We are forced to reduce nb_smp in this particular case*/

        /*FIXME avoid checking for first_line */
        if (ctu_first_line) {
            sub_sample_lm_ref_abv0(lm_src, dst_cb, dst_cr, lm_smp, smp_cb, smp_cr,
                                   lm_src_stride, dst_stride_c,
                                   abv_step, nb_sample_abv, lft_avail);
        } else {
            sub_sample_lm_ref_abv(lm_src, dst_cb, dst_cr, lm_smp, smp_cb, smp_cr,
                                  lm_src_stride, dst_stride_c,
                                  abv_step, nb_sample_abv, lft_avail);
        }

        avg_min_max = sort_average_lm_ref_samples(lm_smp, smp_cb, smp_cr, nb_sample_abv);

        lm_params = derive_cclm_params(&avg_min_max);
    }

    compute_subsample(lm_src, dst_cb, dst_cr, lm_src_stride, dst_stride_c,
                         &lm_params, pb_w, pb_h, lft_avail);
}

static void
vvc_intra_mdlm_top(const uint16_t *const src_luma,
                   uint16_t *const dst_cb, uint16_t *const dst_cr,
                   uint64_t intra_map_rows, int log2_pb_w,
                   int log2_pb_h, int x0, int y0,
                   uint8_t left_available , uint8_t up_available, LMsubsampleFunc const compute_subsample)
{
    vvc_intra_mdlm_t2(src_luma, dst_cb, dst_cr, RCN_CTB_STRIDE, RCN_CTB_STRIDE,
                      log2_pb_w, log2_pb_h, left_available, up_available,
                      y0, x0, intra_map_rows, compute_subsample);
}

static void
vvc_intra_mdlm_t2_collocated(const uint16_t *lm_src, uint16_t *dst_cb, uint16_t *dst_cr,
                             ptrdiff_t lm_src_stride, ptrdiff_t dst_stride_c,
                             int log2_pb_w, int log2_pb_h, uint8_t lft_avail, uint8_t abv_avail,
                             uint8_t y0, uint8_t x0, uint64_t abv_map)
{
    struct CCLMParams lm_params = {
        .cb = {.a = 0, .b = AVG_VAL, .shift = 0},
        .cr = {.a = 0, .b = AVG_VAL, .shift = 0}
    };

    int pb_w = 1 << log2_pb_w;
    int pb_h = 1 << log2_pb_h;

    if (abv_avail) {
        struct AVGMinMax avg_min_max;
        uint16_t lm_smp[4];
        uint16_t smp_cb[4];
        uint16_t smp_cr[4];

        int ref_length_abv = pb_w + OVMIN(pb_h, pb_w);
        int nb_pb_ref_abv = ref_length_abv >> 1;

        int x_pb = x0 >> 1;

        uint64_t abv_msk = (1llu << nb_pb_ref_abv) - 1;
        uint64_t usable_mask = ((abv_map >> (x_pb + 1)) & abv_msk);
        int nb_avail_ref_pb_abv = __builtin_ctzll(~usable_mask);

        int avail_ref_length_abv = nb_avail_ref_pb_abv << 1;

        int nb_sample_abv = OVMIN(avail_ref_length_abv, 4);
        int abv_step = OVMAX(1, (avail_ref_length_abv >> 2));
        uint8_t ctu_first_line = !y0;

        /* in case of abv only ref_length might be 2 while nb_sample_lft is 4
           We are forced to reduce nb_smp in this particular case*/

        /*FIXME avoid checking for first_line */
        if (ctu_first_line) {
            sub_sample_lm_ref_abv0(lm_src, dst_cb, dst_cr, lm_smp, smp_cb, smp_cr,
                                   lm_src_stride, dst_stride_c,
                                   abv_step, nb_sample_abv, lft_avail);
        } else {
            sub_sample_lm_ref_abv_collocated(lm_src, dst_cb, dst_cr, lm_smp, smp_cb, smp_cr,
                                             lm_src_stride, dst_stride_c,
                                             abv_step, nb_sample_abv, lft_avail);
        }

        avg_min_max = sort_average_lm_ref_samples(lm_smp, smp_cb, smp_cr, nb_sample_abv);

        lm_params = derive_cclm_params(&avg_min_max);
    }

    compute_lm_subsample_collocated(lm_src, dst_cb, dst_cr, lm_src_stride, dst_stride_c,
                                    &lm_params, pb_w, pb_h, lft_avail, abv_avail);
}

static void
vvc_intra_mdlm_top_cl(const uint16_t *const src_luma,
                      uint16_t *const dst_cb, uint16_t *const dst_cr,
                      uint64_t intra_map_rows, int log2_pb_w,
                      int log2_pb_h, int x0, int y0,
                      uint8_t left_available , uint8_t up_available, LMsubsampleFunc const compute_subsample)
{
    vvc_intra_mdlm_t2_collocated(src_luma, dst_cb, dst_cr, RCN_CTB_STRIDE, RCN_CTB_STRIDE,
                      log2_pb_w, log2_pb_h, left_available, up_available,
                      y0, x0, intra_map_rows);
}

static void
vvc_intra_mdlm_l2(const uint16_t *lm_src, uint16_t *dst_cb, uint16_t *dst_cr,
                  ptrdiff_t lm_src_stride, ptrdiff_t dst_stride_c,
                  int log2_pb_w, int log2_pb_h, uint8_t lft_avail, uint8_t abv_avail,
                  uint8_t y0, uint8_t x0, uint64_t lft_map, LMsubsampleFunc const compute_subsample)
{
    struct CCLMParams lm_params = {
        .cb = {.a = 0, .b = AVG_VAL, .shift = 0},
        .cr = {.a = 0, .b = AVG_VAL, .shift = 0}
    };

    int pb_w = 1 << log2_pb_w;
    int pb_h = 1 << log2_pb_h;

    if (lft_avail) {
        struct AVGMinMax avg_min_max;
        uint16_t lm_smp[4];
        uint16_t smp_cb[4];
        uint16_t smp_cr[4];

        int ref_length_lft = pb_h + OVMIN(pb_h, pb_w);
        int nb_pb_ref_lft = ref_length_lft >> 1;

        int y_pb = y0 >> 1;

        uint64_t lft_msk = (1llu << nb_pb_ref_lft) - 1;
        uint64_t usable_mask = ((lft_map >> (y_pb + 1)) & lft_msk);
        int nb_avail_ref_pb_lft = __builtin_ctzll(~usable_mask);

        int avail_ref_length_lft = nb_avail_ref_pb_lft << 1;

        int nb_sample_lft = OVMIN(avail_ref_length_lft, 4);
        int lft_step = OVMAX(1, (avail_ref_length_lft >> 2));

        /* in case of lft only ref_length might be 2 while nb_sample_lft is 4
           We are forced to reduce nb_smp in this particular case*/

        sub_sample_lm_ref_lft(lm_src, dst_cb, dst_cr, lm_smp, smp_cb, smp_cr,
                              lm_src_stride, dst_stride_c,
                              lft_step, nb_sample_lft);

        avg_min_max = sort_average_lm_ref_samples(lm_smp, smp_cb, smp_cr, nb_sample_lft);

        lm_params = derive_cclm_params(&avg_min_max);
    }

    compute_subsample(lm_src, dst_cb, dst_cr, lm_src_stride, dst_stride_c,
                      &lm_params, pb_w, pb_h, lft_avail);
}

static void
vvc_intra_mdlm_left(const uint16_t *const src_luma,
                    uint16_t *const dst_cb, uint16_t *const dst_cr,
                    uint64_t intra_map_cols, int log2_pb_w,
                    int log2_pb_h, int x0, int y0,
                    uint8_t left_available, uint8_t up_available, LMsubsampleFunc const compute_subsample)
{
    vvc_intra_mdlm_l2(src_luma, dst_cb, dst_cr, RCN_CTB_STRIDE, RCN_CTB_STRIDE,
                      log2_pb_w, log2_pb_h, left_available, up_available,
                      y0, x0, intra_map_cols, compute_subsample);
}


static void
vvc_intra_mdlm_l2_collocated(const uint16_t *lm_src, uint16_t *dst_cb, uint16_t *dst_cr,
                             ptrdiff_t lm_src_stride, ptrdiff_t dst_stride_c,
                             int log2_pb_w, int log2_pb_h, uint8_t lft_avail, uint8_t abv_avail,
                             uint8_t y0, uint8_t x0, uint64_t lft_map)
{
    struct CCLMParams lm_params = {
        .cb = {.a = 0, .b = AVG_VAL, .shift = 0},
        .cr = {.a = 0, .b = AVG_VAL, .shift = 0}
    };

    int pb_w = 1 << log2_pb_w;
    int pb_h = 1 << log2_pb_h;

    if (lft_avail) {
        struct AVGMinMax avg_min_max;
        uint16_t lm_smp[4];
        uint16_t smp_cb[4];
        uint16_t smp_cr[4];

        int ref_length_lft = pb_h + OVMIN(pb_h, pb_w);
        int nb_pb_ref_lft = ref_length_lft >> 1;

        int y_pb = y0 >> 1;

        uint64_t lft_msk = (1llu << nb_pb_ref_lft) - 1;
        uint64_t usable_mask = ((lft_map >> (y_pb + 1)) & lft_msk);
        int nb_avail_ref_pb_lft = __builtin_ctzll(~usable_mask);

        int avail_ref_length_lft = nb_avail_ref_pb_lft << 1;

        int nb_sample_lft = OVMIN(avail_ref_length_lft, 4);
        int lft_step = OVMAX(1, (avail_ref_length_lft >> 2));

        /* in case of lft only ref_length might be 2 while nb_sample_lft is 4
           We are forced to reduce nb_smp in this particular case*/

        sub_sample_lm_ref_lft_collocated(lm_src, dst_cb, dst_cr, lm_smp, smp_cb, smp_cr,
                                         lm_src_stride, dst_stride_c,
                                         lft_step, nb_sample_lft, abv_avail);

        avg_min_max = sort_average_lm_ref_samples(lm_smp, smp_cb, smp_cr, nb_sample_lft);

        lm_params = derive_cclm_params(&avg_min_max);
    }

    compute_lm_subsample_collocated(lm_src, dst_cb, dst_cr, lm_src_stride, dst_stride_c,
                                    &lm_params, pb_w, pb_h, lft_avail, abv_avail);
}

static void
vvc_intra_mdlm_left_cl(const uint16_t *const src_luma,
                       uint16_t *const dst_cb, uint16_t *const dst_cr,
                       uint64_t intra_map_cols, int log2_pb_w,
                       int log2_pb_h, int x0, int y0,
                       uint8_t left_available, uint8_t up_available, LMsubsampleFunc const compute_subsample)
{
    vvc_intra_mdlm_l2_collocated(src_luma, dst_cb, dst_cr, RCN_CTB_STRIDE, RCN_CTB_STRIDE,
                      log2_pb_w, log2_pb_h, left_available, up_available,
                      y0, x0, intra_map_cols);
}

void
BD_DECL(rcn_init_cclm_functions_collocated)(struct RCNFunctions *rcn_func)
{
   struct CCLMFunctions *const cclm = &rcn_func->cclm; 
   cclm->cclm = &vvc_intra_cclm_cl;
   cclm->mdlm_left = &vvc_intra_mdlm_left_cl;
   cclm->mdlm_top  = &vvc_intra_mdlm_top_cl;
}

/* FIXME check vertical / horizontal */
void
BD_DECL(rcn_init_cclm_functions)(struct RCNFunctions *rcn_func)
{
   struct CCLMFunctions *const cclm = &rcn_func->cclm; 
   cclm->cclm              = &vvc_intra_cclm;
   cclm->mdlm_left         = &vvc_intra_mdlm_left;
   cclm->mdlm_top          = &vvc_intra_mdlm_top;
   cclm->compute_subsample = &compute_lm_subsample;
}

