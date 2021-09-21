#include <emmintrin.h>

#include "rcn_structures.h"

#define DMVR_NUM_ITERATION 2
#define DMVR_SAD_STRIDE  ((2 * DMVR_NUM_ITERATION) + 1)
#define DMVR_NB_IDX (DMVR_SAD_STRIDE * DMVR_SAD_STRIDE)

static const int8_t dmvr_mv_x[25 + 25] = {
    -2, -1, 0, 1, 2,
    -2, -1, 0, 1, 2,
    -2, -1, 0, 1, 2,
    -2, -1, 0, 1, 2,
    -2, -1, 0, 1, 2,
    2, 1, 0, -1, -2,
    2, 1, 0, -1, -2,
    2, 1, 0, -1, -2,
    2, 1, 0, -1, -2,
    2, 1, 0, -1, -2,
};

static const int8_t dmvr_mv_y[25 + 25] = {
    -2, -2, -2, -2, -2,
    -1, -1, -1, -1, -1,
     0,  0,  0,  0,  0,
     1,  1,  1,  1,  1,
     2,  2,  2,  2,  2,
     2,  2,  2,  2,  2,
     1,  1,  1,  1,  1,
     0,  0,  0,  0,  0,
    -1, -1, -1, -1, -1,
    -2, -2, -2, -2, -2,
};

static uint64_t
rcn_dmvr_sad_16(const int16_t *ref0, const int16_t *ref1,
             int16_t dmvr_stride, int16_t pb_w, int16_t pb_h)
{
  uint64_t sum = 0;
  int i;
  for (i = 0; i < (pb_h >> 1); ++i) {
      sum += abs(ref0[0]  - ref1[0]);
      sum += abs(ref0[1]  - ref1[1]);
      sum += abs(ref0[2]  - ref1[2]);
      sum += abs(ref0[3]  - ref1[3]);
      sum += abs(ref0[4]  - ref1[4]);
      sum += abs(ref0[5]  - ref1[5]);
      sum += abs(ref0[6]  - ref1[6]);
      sum += abs(ref0[7]  - ref1[7]);
      sum += abs(ref0[8]  - ref1[8]);
      sum += abs(ref0[9]  - ref1[9]);
      sum += abs(ref0[10] - ref1[10]);
      sum += abs(ref0[11] - ref1[11]);
      sum += abs(ref0[12] - ref1[12]);
      sum += abs(ref0[13] - ref1[13]);
      sum += abs(ref0[14] - ref1[14]);
      sum += abs(ref0[15] - ref1[15]);

      ref0 += dmvr_stride << 1;
      ref1 += dmvr_stride << 1;
  }
  return sum;
}

static uint64_t
rcn_dmvr_sad_8(const int16_t *ref0, const int16_t *ref1,
             int16_t dmvr_stride, int16_t pb_w, int16_t pb_h)
{
  uint64_t sum = 0;
  int i;
  for (i = 0; i < (pb_h >> 1); ++i) {
      sum += abs(ref0[0]  - ref1[0]);
      sum += abs(ref0[1]  - ref1[1]);
      sum += abs(ref0[2]  - ref1[2]);
      sum += abs(ref0[3]  - ref1[3]);
      sum += abs(ref0[4]  - ref1[4]);
      sum += abs(ref0[5]  - ref1[5]);
      sum += abs(ref0[6]  - ref1[6]);
      sum += abs(ref0[7]  - ref1[7]);

      ref0 += dmvr_stride << 1;
      ref1 += dmvr_stride << 1;
  }
  return sum;
}

/*FIXME return min_dmvr_idx; */
static uint8_t
dmvr_compute_sads_16(const int16_t *ref0, const int16_t *ref1,
                  uint64_t *sad_array, int sb_w, int sb_h)
{
    const int32_t stride_l0 = 128 + 4;
    const int32_t stride_l1 = 128 + 4;

    const int16_t *const ref0_start = ref0;
    const int16_t *const ref1_start = ref1;
    uint64_t min_cost = (uint64_t) -1;

    uint8_t idx;
    uint8_t dmvr_idx = 12;

    for (idx = 0; idx < 12; ++idx) {
        ref0 = ref0_start + (int16_t)dmvr_mv_x[idx]
                          + (int16_t)dmvr_mv_y[idx] * stride_l0;

        ref1 = ref1_start - (int16_t)dmvr_mv_x[idx]
                          - (int16_t)dmvr_mv_y[idx] * stride_l1;

        sad_array[idx] = rcn_dmvr_sad_16(ref0, ref1, stride_l1,
                                     sb_w, sb_h);
    }

    for (idx = 13; idx < DMVR_NB_IDX; ++idx) {
        ref0 = ref0_start + (int16_t)dmvr_mv_x[idx]
                          + (int16_t)dmvr_mv_y[idx] * stride_l0;

        ref1 = ref1_start - (int16_t)dmvr_mv_x[idx]
                          - (int16_t)dmvr_mv_y[idx] * stride_l1;

        sad_array[idx] = rcn_dmvr_sad_16(ref0, ref1, stride_l1,
                                     sb_w, sb_h);
    }
    for (idx = 0; idx < DMVR_NB_IDX; ++idx) {
        if (sad_array[idx] < min_cost || (idx == 12 && sad_array[idx] <= min_cost)) {
            min_cost = sad_array[idx];
            dmvr_idx = idx;
        }
    }

    return dmvr_idx;
}

/*FIXME return min_dmvr_idx; */
static uint8_t
dmvr_compute_sads_8(const int16_t *ref0, const int16_t *ref1,
                  uint64_t *sad_array, int sb_w, int sb_h)
{
    const int32_t stride_l0 = 128 + 4;
    const int32_t stride_l1 = 128 + 4;

    const int16_t *const ref0_start = ref0;
    const int16_t *const ref1_start = ref1;
    uint64_t min_cost = (uint64_t) -1;

    uint8_t idx;
    uint8_t dmvr_idx = 12;

    for (idx = 0; idx < 12; ++idx) {
        ref0 = ref0_start + (int16_t)dmvr_mv_x[idx]
                          + (int16_t)dmvr_mv_y[idx] * stride_l0;

        ref1 = ref1_start - (int16_t)dmvr_mv_x[idx]
                          - (int16_t)dmvr_mv_y[idx] * stride_l1;

        sad_array[idx] = rcn_dmvr_sad_8(ref0, ref1, stride_l1,
                                     sb_w, sb_h);
    }

    for (idx = 13; idx < DMVR_NB_IDX; ++idx) {
        ref0 = ref0_start + (int16_t)dmvr_mv_x[idx]
                          + (int16_t)dmvr_mv_y[idx] * stride_l0;

        ref1 = ref1_start - (int16_t)dmvr_mv_x[idx]
                          - (int16_t)dmvr_mv_y[idx] * stride_l1;

        sad_array[idx] = rcn_dmvr_sad_8(ref0, ref1, stride_l1,
                                     sb_w, sb_h);
    }
    for (idx = 0; idx < DMVR_NB_IDX; ++idx) {
        if (sad_array[idx] < min_cost || (idx == 12 && sad_array[idx] <= min_cost)) {
            min_cost = sad_array[idx];
            dmvr_idx = idx;
        }
    }

    return dmvr_idx;
}

void
rcn_dmvr_functions_sse(struct RCNFunctions *const rcn_funcs)
{
    rcn_funcs->dmvr.sad[0] = &rcn_dmvr_sad_8;
    rcn_funcs->dmvr.sad[1] = &rcn_dmvr_sad_16;

    rcn_funcs->dmvr.computeSB[0] = &dmvr_compute_sads_8;
    rcn_funcs->dmvr.computeSB[1] = &dmvr_compute_sads_16;
}
