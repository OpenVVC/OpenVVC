#include <stdint.h>
#include <string.h>

#include "ovutils.h"

#include "rcn_lmcs.h"

#define SMP_RNG (1 << BITDEPTH)

#define BITDEPTH 10
#define ov_bdclip(val) ov_clip_uintp2(val, BITDEPTH);

#define LOG2_NB_WND 4
#define NB_LMCS_WND (1 << LOG2_NB_WND)

#define NB_SMP_WND (SMP_RNG >> LOG2_NB_WND)

#define LOG2_WND_RNG (BITDEPTH - LOG2_NB_WND)
#define WND_RND (1 << (LOG2_WND_RNG - 1))

#define LMCS_PREC 11
#define LMCS_RND (1 << (LMCS_PREC - 1))

#define AVG_VAL (1 << (BITDEPTH - 1))

/* Window information */
struct TMPWindowsInfo
{
    int16_t scaled_fwd_step[NB_LMCS_WND];
    int16_t scaled_bwd_step[NB_LMCS_WND];
    int16_t wnd_bnd[NB_LMCS_WND + 1];
    //int16_t wnd_sz[NB_LMCS_WND];
};

struct LMCSInfo2
{
    uint8_t lmcs_min_bin_idx;
    uint8_t lmcs_delta_max_bin_idx;
    int16_t lmcs_cw_delta[NB_LMCS_WND];
};

/* Backward and forward LUTs */
struct LMCSLUTs
{
    int16_t fwd_lut[SMP_RNG];
    int16_t bwd_lut[SMP_RNG];
    int16_t wnd_bnd[NB_LMCS_WND + 1];
};

static uint8_t
get_bwd_idx(const int16_t *const wnd_bnd, int16_t val, uint8_t min_idx, uint8_t max_idx_plus1)
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
compute_windows_scale_steps(struct TMPWindowsInfo *const tmp_wnd_info,
                            const int16_t *const lmcs_cw_delta,
                            uint8_t min_idx, uint8_t max_idx_plus1)
{
  uint8_t min_idx_plus1 = min_idx + 1;

  int16_t *const fwd_step = tmp_wnd_info->scaled_fwd_step;
  int16_t *const bwd_step = tmp_wnd_info->scaled_bwd_step;
  int16_t *const wnd_bnd = tmp_wnd_info->wnd_bnd;

  int i;

  /* Init left  with zero value until min_idx is reached */
  memset(wnd_bnd, 0, sizeof(*wnd_bnd) * min_idx_plus1);

  /* Init default to zero so padding is already done */
  memset(fwd_step, 0, sizeof(int16_t) << LOG2_NB_WND);
  memset(bwd_step, 0, sizeof(int16_t) << LOG2_NB_WND);

  /* Compute windows */
  for (i = min_idx; i < max_idx_plus1; i++) {
      int16_t wnd_sz = NB_SMP_WND + lmcs_cw_delta[i];
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
derive_forward_lut(int16_t *const fwd_lut, const struct TMPWindowsInfo *const tmp_wnd_info)
{
    const int16_t *const fwd_step = tmp_wnd_info->scaled_fwd_step;
    const int16_t *const wnd_bnd  = tmp_wnd_info->wnd_bnd;
    int16_t val;
    /* FIXME scaled_fwd/bwd_step[i] = 0 if outside of idx ranges
     * Find a way to saturate using min and max idx.
     */
    for (val = 0; val < SMP_RNG; val++) {
        uint8_t wnd_idx = val >> LOG2_WND_RNG;
        int32_t wnd_lbnd = (int32_t)wnd_idx << LOG2_WND_RNG;

        int32_t nb_step = val - wnd_lbnd;

        int fwd_val = wnd_bnd[wnd_idx] + ((fwd_step[wnd_idx] * nb_step + LMCS_RND) >> LMCS_PREC);

        fwd_lut[val] = ov_bdclip(fwd_val);
    }
}

static void
derive_backward_lut(int16_t *const bwd_lut, const struct TMPWindowsInfo *const tmp_wnd_info,
                    uint8_t min_idx, uint8_t max_idx_plus1)
{
    const int16_t *const bwd_step = tmp_wnd_info->scaled_bwd_step;
    const int16_t *const wnd_bnd  = tmp_wnd_info->wnd_bnd;
    int16_t val;
    /* FIXME scaled_fwd/bwd_step[i] = 0 if outside of idx ranges
     * Find a way to saturate using min and max idx.
     */
    for (val = 0; val < SMP_RNG; val++) {
        /* FIXME also construct an inverse idx LUT ? */
        uint8_t wnd_idx = get_bwd_idx(wnd_bnd, val, min_idx, max_idx_plus1);
        int32_t wnd_lbnd = (int32_t)wnd_idx << LOG2_WND_RNG;

        int32_t nb_step = val - wnd_bnd[wnd_idx];

        int bwd_val = wnd_lbnd + ((bwd_step[wnd_idx] * nb_step + LMCS_RND) >> LMCS_PREC);

        bwd_lut[val] = ov_bdclip(bwd_val);
    }
}

/* Note min_idx and max_bin_idx are supposed < 16 */
void
init_lmcs_lut(struct LMCSLUTs *const lmcs_luts, const struct LMCSInfo2 *const lmcs_info)
{
    struct TMPWindowsInfo tmp_wnd;
    uint8_t min_idx = lmcs_info->lmcs_min_bin_idx;
    uint8_t max_idx_plus1 = NB_LMCS_WND - lmcs_info->lmcs_delta_max_bin_idx;

    /* TODO alloc LUTs */
    compute_windows_scale_steps(&tmp_wnd, lmcs_info->lmcs_cw_delta, min_idx, max_idx_plus1);

    derive_forward_lut(lmcs_luts->fwd_lut, &tmp_wnd);

    derive_backward_lut(lmcs_luts->bwd_lut, &tmp_wnd, min_idx, max_idx_plus1);

    /* Only required for chroma scaling */
    memcpy(lmcs_luts->wnd_bnd, tmp_wnd.wnd_bnd, sizeof(lmcs_luts->wnd_bnd));
}

static void
lmcs_convert_data_to_info(struct LMCSInfo2 *const dst, const struct OVLMCSData *const src)
{
    int i;

    dst->lmcs_min_bin_idx = src->lmcs_min_bin_idx;
    dst->lmcs_delta_max_bin_idx = src->lmcs_delta_max_bin_idx;

    memset(dst->lmcs_cw_delta, 0, sizeof(dst->lmcs_cw_delta));

    for (i = dst->lmcs_min_bin_idx; i < NB_LMCS_WND - src->lmcs_delta_max_bin_idx; ++i){
        dst->lmcs_cw_delta[i] = src->lmcs_delta_sign_cw_flag[i] ? -src->lmcs_delta_abs_cw[i]
                                                                :  src->lmcs_delta_abs_cw[i];
    }
}

void
rcn_derive_lmcs_params(struct LMCSInfo *lmcs_info, uint16_t *const output_pivot, const OVLMCSData *const lmcs)
{
    lmcs_info->min_idx = lmcs->lmcs_min_bin_idx;
    lmcs_info->max_idx = 16 - lmcs->lmcs_delta_max_bin_idx;

    /* Keep track of associated parameters */
    lmcs_info->data = lmcs;
}

void 
rcn_lmcs_compute_chroma_scale(struct OVCTUDec* ctudec, int x0, int y0)
{
    struct LMCSInfo* lmcs_info = &ctudec->lmcs_info;
    uint64_t upper_map = ctudec->rcn_ctx.progress_field.hfield[y0 >> 2];
    uint64_t left_map  = ctudec->rcn_ctx.progress_field.vfield[x0 >> 2];
    uint64_t needed_mask = (1 << 16) - 1;
    // uint32_t available_above_mask = ( upper_map & (needed_mask << ((x0 >> 2)+1)) >> ((x0 >> 2)+1));
    uint32_t available_above_mask = ( upper_map & (needed_mask << ((x0 >> 2)+1))) >> ((x0 >> 2)+1);
    // uint32_t available_left_mask  = ( left_map & (needed_mask << ((y0 >> 2)+1)) >> ((y0 >> 2)+1));
    uint32_t available_left_mask  = ( left_map & (needed_mask << ((y0 >> 2)+1)))  >> ((y0 >> 2)+1);
    // uint16_t *_src = &ctudec->ctu_data_y[RCN_CTB_PADDING + x0 + (y0 - 1) * RCN_CTB_STRIDE];
    uint16_t *_src = &ctudec->rcn_ctx.ctu_buff.y[x0 + (y0 - 1) * RCN_CTB_STRIDE];

    uint32_t luma_sum1 = 0;
    uint32_t luma_sum2 = 0;
    uint32_t luma_sum3 = 0;
    uint32_t luma_sum4 = 0;
    uint8_t num_luma_pu = 0;
    uint32_t log2_num_luma_pu = 0;
    uint32_t luma_avg = AVG_VAL;
    int idx = lmcs_info->min_idx;
    
    uint8_t num_luma_pu_above = 0;
    while (available_above_mask){
        luma_sum1 += _src[0];
        luma_sum2 += _src[1];
        luma_sum3 += _src[2];
        luma_sum4 += _src[3];
        _src += 4;
        ++num_luma_pu_above;
        available_above_mask >>= 1;
    }
    //When image border, num_luma_pum may be different from 16.
    if(num_luma_pu_above){
        _src -= 4;
        luma_sum1 += _src[3]*(16-num_luma_pu_above);
        luma_sum2 += _src[3]*(16-num_luma_pu_above);
        luma_sum3 += _src[3]*(16-num_luma_pu_above);
        luma_sum4 += _src[3]*(16-num_luma_pu_above);
        num_luma_pu_above = 16;
    }

    // _src = &ctudec->ctu_data_y[RCN_CTB_PADDING + x0 - 1 + y0  * RCN_CTB_STRIDE];
    _src = &ctudec->rcn_ctx.ctu_buff.y[x0 - 1 + y0  * RCN_CTB_STRIDE];
    uint8_t num_luma_pu_left = 0;
    while (available_left_mask){
        luma_sum1 += _src[0];
        luma_sum2 += _src[RCN_CTB_STRIDE];
        luma_sum3 += _src[RCN_CTB_STRIDE << 1];
        luma_sum4 += _src[RCN_CTB_STRIDE * 3];
        _src += RCN_CTB_STRIDE << 2;
        ++num_luma_pu_left;
        available_left_mask >>= 1;
    }
    //When image border, num_luma_pum may be different from 16.
    if(num_luma_pu_left){
        _src -= RCN_CTB_STRIDE << 2;
        #if 0
        luma_sum1 += _src[0]*(16-num_luma_pu_left);
        luma_sum2 += _src[RCN_CTB_STRIDE]*(16-num_luma_pu_left);
        luma_sum3 += _src[RCN_CTB_STRIDE << 1]*(16-num_luma_pu_left);
        #endif
        luma_sum1 += _src[RCN_CTB_STRIDE * 3]*(16-num_luma_pu_left);
        luma_sum2 += _src[RCN_CTB_STRIDE * 3]*(16-num_luma_pu_left);
        luma_sum3 += _src[RCN_CTB_STRIDE * 3]*(16-num_luma_pu_left);
        luma_sum4 += _src[RCN_CTB_STRIDE * 3]*(16-num_luma_pu_left);
        num_luma_pu_left = 16;
    }

    num_luma_pu = num_luma_pu_left + num_luma_pu_above;
    while(num_luma_pu){
        ++log2_num_luma_pu;
        num_luma_pu >>= 1;
    }
    //decrease of 1 if not zero
    log2_num_luma_pu -= !!log2_num_luma_pu;

    luma_avg = log2_num_luma_pu ? ((luma_sum1 + luma_sum2) + (luma_sum3 + luma_sum4) + (1 << (log2_num_luma_pu + 1))) >> (log2_num_luma_pu + 2) : AVG_VAL;

    idx = get_bwd_idx((int16_t*)lmcs_info->lmcs_output_pivot, luma_avg,
                      lmcs_info->min_idx, lmcs_info->max_idx);

    /* FIXME use coded window size instead ? */
    uint32_t size_interval = (uint32_t)(lmcs_info->lmcs_output_pivot[idx + 1] - lmcs_info->lmcs_output_pivot[idx]);
    lmcs_info->lmcs_chroma_scale = (size_interval == 0) ? 1 << 11 
                                                    : (1<< (BITDEPTH - LOG2_NB_WND + 11))/ (size_interval + lmcs_info->lmcs_chroma_scaling_offset);
}

static void 
rcn_lmcs_reshape_luma_blk_lut(uint16_t *dst, ptrdiff_t stride_dst, uint16_t* lmcs_lut_luma,
                              int width, int height)
{
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++){ 
            dst[x] = lmcs_lut_luma[dst[x]];
        }
        dst += stride_dst;
    }
}

void 
rcn_lmcs_compute_lut_luma(struct LMCSInfo *lmcs_info, uint16_t* inverse_lut,
                          uint16_t* forward_lut, uint16_t* wnd_boundaries)
{
    struct LMCSInfo2 new_lmcs;
    struct LMCSLUTs new_luts;
    lmcs_convert_data_to_info(&new_lmcs, lmcs_info->data);

    init_lmcs_lut(&new_luts, &new_lmcs);

    memcpy(forward_lut, new_luts.fwd_lut, sizeof(new_luts.fwd_lut));
    memcpy(inverse_lut, new_luts.bwd_lut, sizeof(new_luts.bwd_lut));
    memcpy(wnd_boundaries, &new_luts.wnd_bnd, sizeof(new_luts.wnd_bnd) - 2);
}


void
rcn_lmcs_no_reshape(uint16_t *dst, ptrdiff_t stride_dst, uint16_t* lmcs_lut_luma,
                    int width, int height)
{
    return;
}

void
rcn_init_lmcs_function(struct RCNFunctions *rcn_func, uint8_t lmcs_flag)
{
    if(lmcs_flag){
        rcn_func->lmcs_reshape = &rcn_lmcs_reshape_luma_blk_lut;
    } else {
        rcn_func->lmcs_reshape = &rcn_lmcs_no_reshape;
    }
}
