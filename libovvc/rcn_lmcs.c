#include <stdint.h>
#include <string.h>

#include "ovutils.h"
#include "ovmem.h"

#include "nvcl_structures.h"
#include "rcn_structures.h"
#include "rcn_lmcs.h"


#define SMP_RNG (1 << BITDEPTH)

#include "bitdepth.h"

#define LOG2_NB_WND 4
#define NB_LMCS_WND (1 << LOG2_NB_WND)

#define NB_SMP_WND (SMP_RNG >> LOG2_NB_WND)

#define LOG2_WND_RNG (BITDEPTH - LOG2_NB_WND)
#define WND_RND (1 << (LOG2_WND_RNG - 1))

#define LMCS_PREC 11
#define LMCS_RND (1 << (LMCS_PREC - 1))

/* Window information */
struct WindowsInfo
{
    uint16_t scaled_fwd_step[NB_LMCS_WND];
    uint16_t scaled_bwd_step[NB_LMCS_WND];
    OVSample wnd_bnd[NB_LMCS_WND + 1];
    //int16_t wnd_sz[NB_LMCS_WND];
};

struct LMCSParams
{
    uint8_t min_bin_idx;
    uint8_t delta_max_bin_idx;
    int16_t cw_delta[NB_LMCS_WND];
};

/* Backward and forward LUTs */
struct LMCSLUTs
{
    OVSample fwd_lut[SMP_RNG];
    OVSample bwd_lut[SMP_RNG];
    OVSample wnd_bnd[NB_LMCS_WND + 1];
};

static uint8_t
get_bwd_idx(const OVSample *const wnd_bnd, OVSample val, uint8_t min_idx, uint8_t max_idx_plus1)
{
  uint8_t i = min_idx;
  for (; i < max_idx_plus1; i++) {
      if (val < wnd_bnd[i + 1]) {
          break;
      }
  }
  return (uint8_t)OVMIN(i, NB_LMCS_WND - 1);
}

static void
compute_windows_scale_steps(struct WindowsInfo *const wnd_info,
                            const int16_t *const cw_delta,
                            uint8_t min_idx, uint8_t max_idx_plus1)
{
  uint8_t min_idx_plus1 = min_idx + 1;

  uint16_t *const fwd_step = wnd_info->scaled_fwd_step;
  uint16_t *const bwd_step = wnd_info->scaled_bwd_step;
  OVSample *const wnd_bnd  = wnd_info->wnd_bnd;

  int i;

  /* Init left  with zero value until min_idx is reached */
  memset(wnd_bnd, 0, sizeof(*wnd_bnd) * min_idx_plus1);

  /* Init default to zero so padding is already done */
  memset(fwd_step, 0, sizeof(*fwd_step) << LOG2_NB_WND);
  memset(bwd_step, 0, sizeof(*bwd_step) << LOG2_NB_WND);

  /* Compute windows */
  for (i = min_idx; i < max_idx_plus1; i++) {
      int32_t wnd_sz = NB_SMP_WND + cw_delta[i];
      if (wnd_sz) {
          fwd_step[i] = ((wnd_sz << LMCS_PREC) + WND_RND) >> LOG2_WND_RNG;
          bwd_step[i] = (NB_SMP_WND << LMCS_PREC) / wnd_sz;
      }
      wnd_bnd[i + 1] = wnd_bnd[i] + wnd_sz;
  }

  /* Broadcast last value */
  for (i = max_idx_plus1; i < NB_LMCS_WND; i++) {
      wnd_bnd[i + 1] = wnd_bnd[i];
  }

}

static void
derive_forward_lut(OVSample *const fwd_lut, const struct WindowsInfo *const wnd_info)
{
    const uint16_t *const fwd_step = wnd_info->scaled_fwd_step;
    const OVSample *const wnd_bnd  = wnd_info->wnd_bnd;
    uint16_t val;

    for (val = 0; val < SMP_RNG; val++) {
        uint8_t wnd_idx = val >> LOG2_WND_RNG;
        int32_t wnd_lbnd = (int32_t)wnd_idx << LOG2_WND_RNG;

        int32_t nb_step = val - wnd_lbnd;

        int fwd_val = (int16_t)wnd_bnd[wnd_idx] + (((int32_t)fwd_step[wnd_idx] * nb_step + LMCS_RND) >> LMCS_PREC);

        fwd_lut[val] = ov_bdclip(fwd_val);
    }
}

static void
derive_backward_lut(OVSample *const bwd_lut, const struct WindowsInfo *const wnd_info,
                    uint8_t min_idx, uint8_t max_idx_plus1)
{
    const uint16_t *const bwd_step = wnd_info->scaled_bwd_step;
    const OVSample *const wnd_bnd  = wnd_info->wnd_bnd;
    uint16_t val;

    for (val = 0; val < SMP_RNG; val++) {
        uint8_t wnd_idx = get_bwd_idx(wnd_bnd, val, min_idx, max_idx_plus1);
        int32_t wnd_lbnd = (int32_t)wnd_idx << LOG2_WND_RNG;

        int32_t nb_step = val - wnd_bnd[wnd_idx];

        int bwd_val = wnd_lbnd + (((int32_t)bwd_step[wnd_idx] * nb_step + LMCS_RND) >> LMCS_PREC);

        bwd_lut[val] = ov_bdclip(bwd_val);
    }
}

static void
init_lmcs_lut(struct LMCSLUTs *const lmcs_luts, const struct LMCSParams *const params)
{
    struct WindowsInfo tmp_wnd;
    uint8_t min_idx = params->min_bin_idx;
    uint8_t max_idx_plus1 = NB_LMCS_WND - params->delta_max_bin_idx;

    compute_windows_scale_steps(&tmp_wnd, params->cw_delta, min_idx, max_idx_plus1);

    derive_forward_lut(lmcs_luts->fwd_lut, &tmp_wnd);

    derive_backward_lut(lmcs_luts->bwd_lut, &tmp_wnd, min_idx, max_idx_plus1);

    /* Only required for chroma scaling */
    memcpy(lmcs_luts->wnd_bnd, tmp_wnd.wnd_bnd, sizeof(lmcs_luts->wnd_bnd));
}

static void
lmcs_convert_data_to_info(struct LMCSParams *const dst, const struct OVLMCSData *const src)
{
    int i;

    dst->min_bin_idx       = src->lmcs_min_bin_idx;
    dst->delta_max_bin_idx = src->lmcs_delta_max_bin_idx;

    memset(dst->cw_delta, 0, sizeof(dst->cw_delta));

    for (i = dst->min_bin_idx; i < NB_LMCS_WND - dst->delta_max_bin_idx; ++i){
        dst->cw_delta[i] = src->lmcs_delta_sign_cw_flag[i] ? -src->lmcs_delta_abs_cw[i]
                                                           :  src->lmcs_delta_abs_cw[i];
    }
}

static uint32_t
lmcs_compute_luma_average(const OVSample *src, uint32_t abv_mask, uint32_t lft_mask)
{
    const OVSample *_src = src - RCN_CTB_STRIDE;

    uint32_t luma_sum1 = 0;
    uint32_t luma_sum2 = 0;
    uint32_t luma_sum3 = 0;
    uint32_t luma_sum4 = 0;

    uint8_t nb_units = 0;
    uint32_t log2_nb_units = 0;
    uint32_t luma_avg = AVG_VAL;

    uint8_t nb_units_abv = 0;
    while (abv_mask) {
        luma_sum1 += _src[0];
        luma_sum2 += _src[1];
        luma_sum3 += _src[2];
        luma_sum4 += _src[3];
        _src += 4;
        ++nb_units_abv;
        abv_mask >>= 1;
    }

    if (nb_units_abv) {
        uint32_t pad_val = _src[-1] * (16 - nb_units_abv);
        luma_sum1 += pad_val;
        luma_sum2 += pad_val;
        luma_sum3 += pad_val;
        luma_sum4 += pad_val;
        nb_units_abv = 16;
    }

    _src = src - 1;

    uint8_t nb_units_lft = 0;
    while (lft_mask) {
        luma_sum1 += _src[0];
        luma_sum2 += _src[RCN_CTB_STRIDE];
        luma_sum3 += _src[RCN_CTB_STRIDE << 1];
        luma_sum4 += _src[RCN_CTB_STRIDE * 3];
        _src += RCN_CTB_STRIDE << 2;
        ++nb_units_lft;
        lft_mask >>= 1;
    }

    if (nb_units_lft) {
        uint32_t pad_val = _src[-RCN_CTB_STRIDE] * (16 - nb_units_lft);
        luma_sum1 += pad_val;
        luma_sum2 += pad_val;
        luma_sum3 += pad_val;
        luma_sum4 += pad_val;
        nb_units_lft = 16;
    }

    nb_units = nb_units_abv + nb_units_lft;

    /* FIXME ctz */
    while (nb_units) {
        ++log2_nb_units;
        nb_units >>= 1;
    }

    luma_avg = log2_nb_units ? ((luma_sum1 + luma_sum2) + (luma_sum3 + luma_sum4) + (1 << log2_nb_units)) >> (log2_nb_units + 1) : AVG_VAL;

    return luma_avg;
}

static void
rcn_lmcs_reshape_luma_blk_lut(OVSample *dst, ptrdiff_t stride_dst, const OVSample *const lmcs_lut_luma,
                              int width, int height)
{
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++){
            dst[x] = lmcs_lut_luma[dst[x]];
        }
        dst += stride_dst;
    }
}

static void
rcn_lmcs_reshape_forward(OVSample *dst, ptrdiff_t stride_dst,
                         const struct LMCSLUTs *const luts,
                         int width, int height)
{
    rcn_lmcs_reshape_luma_blk_lut(dst, stride_dst, luts->fwd_lut, width, height);
}

static void
rcn_lmcs_reshape_backward(OVSample *dst, ptrdiff_t stride_dst,
                          const struct LMCSLUTs *const luts,
                          int width, int height)
{
    rcn_lmcs_reshape_luma_blk_lut(dst, stride_dst, luts->bwd_lut, width, height);
}

static void
rcn_lmcs_compute_lut_luma(struct LMCSInfo *lmcs_info, const struct OVLMCSData *const data)
{
    struct LMCSParams params;

    lmcs_convert_data_to_info(&params, data);

    init_lmcs_lut(lmcs_info->luts, &params);
}

static void
rcn_lmcs_no_reshape(OVSample *dst, ptrdiff_t stride_dst,
                    const struct LMCSLUTs *const luts,
                    int width, int height)
{
    return;
}

static void
rcn_lmcs_compute_chroma_scale(struct LMCSInfo *const lmcs_info,
                              const struct CTUBitField *const progress_field,
                              const OVSample *ctu_data_y, uint8_t x0, uint8_t y0)
{
    uint8_t x0_unit = x0 >> 2;
    uint8_t y0_unit = y0 >> 2;
    uint64_t abv_map = progress_field->hfield[y0_unit];
    uint64_t lft_map = progress_field->vfield[x0_unit];
    uint64_t needed_mask = (1 << 16) - 1;
    uint32_t abv_mask = (abv_map >> (x0_unit + 1)) & needed_mask;
    uint32_t lft_mask = (lft_map >> (y0_unit + 1)) & needed_mask;

    const OVSample *src = &ctu_data_y[x0 + y0 * RCN_CTB_STRIDE];

    uint32_t luma_avg = lmcs_compute_luma_average(src, abv_mask, lft_mask);

    int idx = get_bwd_idx(lmcs_info->luts->wnd_bnd, luma_avg, lmcs_info->min_idx, lmcs_info->max_idx);

    int32_t wnd_sz = (int32_t)(lmcs_info->luts->wnd_bnd[idx + 1] - lmcs_info->luts->wnd_bnd[idx]);

    lmcs_info->lmcs_chroma_scale = (wnd_sz == 0) ? 1 << LMCS_PREC
                                                 : (1 << (BITDEPTH - LOG2_NB_WND + LMCS_PREC)) / (wnd_sz + lmcs_info->lmcs_chroma_scaling_offset);
}

static void
rcn_init_lmcs(struct LMCSInfo *lmcs_info, const struct OVLMCSData *const lmcs_data)
{
    if (!lmcs_info->luts) {
        lmcs_info->luts = ov_malloc(sizeof(struct LMCSLUTs));
    }

    lmcs_info->min_idx = lmcs_data->lmcs_min_bin_idx;
    lmcs_info->max_idx = NB_LMCS_WND - lmcs_data->lmcs_delta_max_bin_idx;

    rcn_lmcs_compute_lut_luma(lmcs_info, lmcs_data);

    lmcs_info->lmcs_chroma_scaling_offset = lmcs_data->lmcs_delta_sign_crs_flag ?
        -lmcs_data->lmcs_delta_abs_crs
        : lmcs_data->lmcs_delta_abs_crs;

}

void
BD_DECL(rcn_init_lmcs_function)(struct RCNFunctions *rcn_func, uint8_t lmcs_flag)
{
    if(lmcs_flag){
        rcn_func->lmcs_reshape_forward  = &rcn_lmcs_reshape_forward;
        rcn_func->lmcs_reshape_backward = &rcn_lmcs_reshape_backward;
        rcn_func->rcn_lmcs_compute_chroma_scale = &rcn_lmcs_compute_chroma_scale;
        rcn_func->rcn_init_lmcs = &rcn_init_lmcs;
    } else {
        rcn_func->lmcs_reshape_forward  = &rcn_lmcs_no_reshape;
        rcn_func->lmcs_reshape_backward = &rcn_lmcs_no_reshape;
        rcn_func->rcn_lmcs_compute_chroma_scale = &rcn_lmcs_compute_chroma_scale;
        rcn_func->rcn_init_lmcs = &rcn_init_lmcs;
    }
}
