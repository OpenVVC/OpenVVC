#include <stdint.h>

#include "rcn_structures.h"
#include "ovutils.h"
#include "ctudec.h"

#define SWAP(type,a,b) do{type tmp= b; b = a; a = tmp;}while(0)

#define W_SHIFT 1
#define H_SHIFT 1
#if 0
static void
vvc_intra_cclm_cl(const uint16_t* const src_luma, uint16_t* const dst_cb,
                  uint16_t* const dst_cr, int log2_pb_w, int log2_pb_h, int y0,
                  int up_available, int left_available)
{ // left and above
        int i, j, i_s;

        uint16_t lm_ref[4];
        uint16_t cb_ref[4];
        uint16_t cr_ref[4];
        VVCLMParams lm_params;
        const int src_stride = RCN_CTB_STRIDE;
        const int dst_stride = RCN_CTB_STRIDE;

        // ref above
        if (up_available && left_available) {
                const uint8_t first_row_in_ctu = !y0;
                const uint16_t *_src_cb, *_src_cr, *_src;
                int step_above = OVMAX(1, (1 << log2_pb_w) >> 1);
                int step_left = OVMAX(1, (1 << log2_pb_h) >> 1);

                int num_sample_l = OVMIN(1 << log2_pb_h, 2); // = 2
                int num_sample_a = OVMIN(1 << log2_pb_w, 2); // = 2

                int offset_pos_a = (1 << log2_pb_w) >> 2;
                int offset_pos_l = (1 << log2_pb_h) >> 2;
                if (first_row_in_ctu) {
                        _src = src_luma - src_stride;
                        _src_cb = dst_cb - src_stride;
                        _src_cr = dst_cr - src_stride;
                        for (i = 0, i_s = offset_pos_a << W_SHIFT;
                             i < num_sample_a;
                             ++i, i_s += step_above << W_SHIFT) {
                                lm_ref[i] = (_src[i_s - 1] + (_src[i_s] << 1) +
                                             _src[i_s + 1] + 2) >>
                                            2;
                                cb_ref[i] = _src_cb[i_s >> W_SHIFT];
                                cr_ref[i] = _src_cr[i_s >> W_SHIFT];
                        }
                } else {
                        _src = src_luma - (src_stride << H_SHIFT);
                        _src_cb = dst_cb - src_stride;
                        _src_cr = dst_cr - src_stride;
                        for (i = 0, i_s = offset_pos_a << W_SHIFT;
                             i < num_sample_a;
                             ++i, i_s += step_above << W_SHIFT) {
                                lm_ref[i] =
                                  (_src[i_s - src_stride] + _src[i_s - 1] +
                                   (_src[i_s] << 2) + _src[i_s + 1] +
                                   _src[i_s + src_stride] + 4) >>
                                  3;
                                cb_ref[i] = _src_cb[i_s >> W_SHIFT];
                                cr_ref[i] = _src_cr[i_s >> W_SHIFT];
                        }
                }

                // ref_left
                // FIXME check - 2 - W_SHIFT
                _src = src_luma - 2 - W_SHIFT +
                       (offset_pos_l << H_SHIFT) * src_stride;
                _src_cb = dst_cb - 1 + offset_pos_l * src_stride;
                _src_cr = dst_cr - 1 + offset_pos_l * src_stride;

                for (j = num_sample_a; j < num_sample_l + num_sample_a; ++j) {
                        lm_ref[j] =
                          (_src[1 - src_stride] + _src[0] + (_src[1] << 2) +
                           _src[2] + _src[1 + src_stride] + 4) >>
                          3;
                        cb_ref[j] = _src_cb[0];
                        cr_ref[j] = _src_cr[0];
                        _src += (src_stride * step_left) << H_SHIFT;
                        _src_cb += (src_stride * step_left);
                        _src_cr += (src_stride * step_left);
                }
        } else if (left_available) {
                const uint16_t *_src_cb, *_src_cr, *_src;
                int step_left = OVMAX(1, (1 << log2_pb_h) >> 2);
                int num_sample = OVMIN(1 << log2_pb_h, 4);
                int offset_pos_l = (1 << log2_pb_h) >> 3;

                _src = src_luma - 2 - W_SHIFT +
                       (offset_pos_l << H_SHIFT) * src_stride;
                _src_cb = dst_cb - 1 + offset_pos_l * src_stride;
                _src_cr = dst_cr - 1 + offset_pos_l * src_stride;

                if (!offset_pos_l) {
                        // FIXME first position

                        lm_ref[0] =
                          (_src[0] + (_src[1] << 1) + _src[2] + 2) >> 2;
                        cb_ref[0] = _src_cb[0];
                        cr_ref[0] = _src_cr[0];

                        _src += (src_stride * step_left) << H_SHIFT;
                        _src_cb += src_stride * step_left;
                        _src_cr += src_stride * step_left;
                }

                for (j = !offset_pos_l; j < num_sample; ++j) {
                        lm_ref[j] =
                          (_src[1 - src_stride] + _src[0] + (_src[1] << 2) +
                           _src[2] + _src[1 + src_stride] + 4) >>
                          3;
                        cb_ref[j] = _src_cb[0];
                        cr_ref[j] = _src_cr[0];
                        _src += (src_stride * step_left) << H_SHIFT;
                        _src_cb += src_stride * step_left;
                        _src_cr += src_stride * step_left;
                }
                if (num_sample == 2) {
                        lm_ref[3] = lm_ref[0];
                        cb_ref[3] = cb_ref[0];
                        cr_ref[3] = cr_ref[0];
                        lm_ref[2] = lm_ref[1];
                        cb_ref[2] = cb_ref[1];
                        cr_ref[2] = cr_ref[1];
                        lm_ref[0] = lm_ref[1];
                        cb_ref[0] = cb_ref[1];
                        cr_ref[0] = cr_ref[1];
                        lm_ref[1] = lm_ref[3];
                        cb_ref[1] = cb_ref[3];
                        cr_ref[1] = cr_ref[3];
                }
        } else if (up_available) {
                const uint8_t first_row_in_ctu = !y0;
                const uint16_t *_src_cb, *_src_cr, *_src;
                int step_above = OVMAX(1, (1 << log2_pb_w) >> 2);
                int num_sample = OVMIN(1 << log2_pb_w, 4);
                int offset_a = (1 << log2_pb_w) >> 3;

                // FIXME positions
                if (first_row_in_ctu) {
                        _src = src_luma - src_stride;
                        _src_cb = dst_cb - src_stride;
                        _src_cr = dst_cr - src_stride;
                        if (!offset_a) {
                                lm_ref[0] = _src[0];
                                cb_ref[0] = _src_cb[0];
                                cr_ref[0] = _src_cr[0];
                                _src += step_above << W_SHIFT;
                                _src_cb += step_above;
                                _src_cr += step_above;
                        }
                        for (i = !offset_a, i_s = offset_a << W_SHIFT;
                             i < num_sample;
                             ++i, i_s += step_above << W_SHIFT) {
                                lm_ref[i] = (_src[i_s - 1] + (_src[i_s] << 1) +
                                             _src[i_s + 1] + 2) >>
                                            2;
                                cb_ref[i] = _src_cb[i_s >> W_SHIFT];
                                cr_ref[i] = _src_cr[i_s >> W_SHIFT];
                        }
                } else {
                        _src = src_luma - (src_stride << H_SHIFT);
                        _src_cb = dst_cb - src_stride;
                        _src_cr = dst_cr - src_stride;
                        if (!offset_a) {
                                lm_ref[0] =
                                  (_src[0 - src_stride] + (_src[0] << 1) +
                                   _src[0 + src_stride] + 2) >>
                                  2;
                                cb_ref[0] = _src_cb[0];
                                cr_ref[0] = _src_cr[0];
                                _src += step_above << W_SHIFT;
                                _src_cb += step_above;
                                _src_cr += step_above;
                        }
                        for (i = !offset_a, i_s = offset_a << W_SHIFT;
                             i < num_sample;
                             ++i, i_s += step_above << W_SHIFT) {
                                lm_ref[i] =
                                  (_src[i_s - src_stride] + _src[i_s - 1] +
                                   (_src[i_s] << 2) + _src[i_s + 1] +
                                   _src[i_s + src_stride] + 4) >>
                                  3;
                                cb_ref[i] = _src_cb[i_s >> W_SHIFT];
                                cr_ref[i] = _src_cr[i_s >> W_SHIFT];
                        }
                }
                if (num_sample == 2) {
                        lm_ref[3] = lm_ref[0];
                        cb_ref[3] = cb_ref[0];
                        cr_ref[3] = cr_ref[0];
                        lm_ref[2] = lm_ref[1];
                        cb_ref[2] = cb_ref[1];
                        cr_ref[2] = cr_ref[1];
                        lm_ref[0] = lm_ref[1];
                        cb_ref[0] = cb_ref[1];
                        cr_ref[0] = cr_ref[1];
                        lm_ref[1] = lm_ref[3];
                        cb_ref[1] = cb_ref[3];
                        cr_ref[1] = cr_ref[3];
                }
        }

        // TODO compute lm_param for cb and cr components
        derive_lm_parameters(lm_ref, cb_ref, cr_ref, &lm_params);

        // inner part from reconstructed picture buffer
        // TODO directly scale cb and cr components;
        if (up_available && left_available) {
                int32_t value;
                uint16_t *_dst_cb, *_dst_cr;
                const uint16_t* _src;
                int offset_cb = lm_params.offset_cb;
                int offset_cr = lm_params.offset_cr;
                int scale_cb = lm_params.scale_cb;
                int scale_cr = lm_params.scale_cr;
                int shift_cb = lm_params.shift_cb;
                int shift_cr = lm_params.shift_cr;
                _dst_cb = dst_cb;
                _dst_cr = dst_cr;
                _src = src_luma;
                for (j = 0; j < 1 << log2_pb_h; ++j) {
                        for (i = 0, i_s = 0; i < 1 << log2_pb_w;
                             ++i, i_s = i << W_SHIFT) {
                                value =
                                  (_src[i_s - src_stride] + _src[i_s - 1] +
                                   (_src[i_s] << 2) + _src[i_s + 1] +
                                   _src[i_s + src_stride] + 4) >>
                                  3;
                                _dst_cb[i] = ov_clip(
                                  ((value * scale_cb) >> shift_cb) + offset_cb,
                                  0,
                                  1023);
                                _dst_cr[i] = ov_clip(
                                  ((value * scale_cr) >> shift_cr) + offset_cr,
                                  0,
                                  1023);
                        }
                        _dst_cb += dst_stride;
                        _dst_cr += dst_stride;
                        _src += src_stride << H_SHIFT;
                }
        } else if (left_available) {
                int32_t value;
                uint16_t *_dst_cb, *_dst_cr;
                const uint16_t* _src;
                int offset_cb = lm_params.offset_cb;
                int offset_cr = lm_params.offset_cr;
                int scale_cb = lm_params.scale_cb;
                int scale_cr = lm_params.scale_cr;
                int shift_cb = lm_params.shift_cb;
                int shift_cr = lm_params.shift_cr;
                _dst_cb = dst_cb;
                _dst_cr = dst_cr;
                _src = src_luma;
                for (i = 0, i_s = 0; i < 1 << log2_pb_w;
                     ++i, i_s = i << W_SHIFT) {
                        value = (_src[i_s - 1] + (_src[i_s] << 1) +
                                 _src[i_s + 1] + 2) >>
                                2;
                        _dst_cb[i] =
                          ov_clip(((value * scale_cb) >> shift_cb) + offset_cb,
                                  0,
                                  1023);
                        _dst_cr[i] =
                          ov_clip(((value * scale_cr) >> shift_cr) + offset_cr,
                                  0,
                                  1023);
                }

                _dst_cb += dst_stride;
                _dst_cr += dst_stride;
                _src += src_stride << H_SHIFT;

                for (j = 1; j < 1 << log2_pb_h; ++j) {
                        for (i = 0, i_s = 0; i < 1 << log2_pb_w;
                             ++i, i_s = i << W_SHIFT) {
                                value =
                                  (_src[i_s - src_stride] + _src[i_s - 1] +
                                   (_src[i_s] << 2) + _src[i_s + 1] +
                                   _src[i_s + src_stride] + 4) >>
                                  3;
                                _dst_cb[i] = ov_clip(
                                  ((value * scale_cb) >> shift_cb) + offset_cb,
                                  0,
                                  1023);
                                _dst_cr[i] = ov_clip(
                                  ((value * scale_cr) >> shift_cr) + offset_cr,
                                  0,
                                  1023);
                        }
                        _dst_cb += dst_stride;
                        _dst_cr += dst_stride;
                        _src += src_stride << H_SHIFT;
                }
        } else if (up_available) {
                int32_t value;
                uint16_t *_dst_cb, *_dst_cr;
                const uint16_t* _src;
                int offset_cb = lm_params.offset_cb;
                int offset_cr = lm_params.offset_cr;
                int scale_cb = lm_params.scale_cb;
                int scale_cr = lm_params.scale_cr;
                int shift_cb = lm_params.shift_cb;
                int shift_cr = lm_params.shift_cr;
                _dst_cb = dst_cb;
                _dst_cr = dst_cr;
                _src = src_luma;
                for (j = 0; j < 1 << log2_pb_h; ++j) {
                        value = (_src[0 - src_stride] + (_src[0] << 1) +
                                 _src[0 + src_stride] + 2) >>
                                2;
                        _dst_cb[0] =
                          ov_clip(((value * scale_cb) >> shift_cb) + offset_cb,
                                  0,
                                  1023);
                        _dst_cr[0] =
                          ov_clip(((value * scale_cr) >> shift_cr) + offset_cr,
                                  0,
                                  1023);

                        for (i = 1, i_s = 1 << W_SHIFT; i < 1 << log2_pb_w;
                             ++i, i_s = i << W_SHIFT) {
                                value =
                                  (_src[i_s - src_stride] + _src[i_s - 1] +
                                   (_src[i_s] << 2) + _src[i_s + 1] +
                                   _src[i_s + src_stride] + 4) >>
                                  3;
                                _dst_cb[i] = ov_clip(
                                  ((value * scale_cb) >> shift_cb) + offset_cb,
                                  0,
                                  1023);
                                _dst_cr[i] = ov_clip(
                                  ((value * scale_cr) >> shift_cr) + offset_cr,
                                  0,
                                  1023);
                        }
                        _dst_cb += dst_stride;
                        _dst_cr += dst_stride;
                        _src += src_stride << H_SHIFT;
                }
        }
}

static void
vvc_intra_mdlm_top_cl(const uint16_t* const src_luma, uint16_t* const dst_cb,
                      uint16_t* const dst_cr, uint64_t intra_map_rows,
                      int log2_pb_w, int log2_pb_h, int x0, int y0,
                      uint8_t left_available, uint8_t up_available)
{ // left and above
        int i, j, i_s;

        uint16_t lm_ref[4];
        uint16_t cb_ref[4];
        uint16_t cr_ref[4];
        VVCLMParams lm_params;
        const int src_stride = RCN_CTB_STRIDE;
        const int dst_stride = RCN_CTB_STRIDE;

        int x0_u = x0 >> 1;
        int above_ref_length_u = (1 << log2_pb_w) >> 1;

        above_ref_length_u +=
          OVMIN((1 << log2_pb_h) >> 1, (1 << log2_pb_w) >> 1);
        unsigned int mask_shift = 63 - (x0_u + above_ref_length_u);

        uint64_t mask_size = ((uint64_t)1 << (above_ref_length_u + 1)) - 1;
        uint64_t needed_mask = mask_size << mask_shift;
        uint64_t usable_mask = intra_map_rows & needed_mask;

        // We use PU width here since this won't change lm paramters
        // int sum_luma = 0, sum_chroma = 0, sum_xx = 0, sum_xy = 0;
        int min_log2_size;
        int min_size;
        if (!(usable_mask ^ needed_mask)) {
                min_size = above_ref_length_u << 1;
        } else {
                min_size = above_ref_length_u << 1;
                int remaining = 0;
                uint64_t mask_rem = (usable_mask ^ needed_mask) >> mask_shift;
                uint64_t mask_size2 =
                  ((uint64_t)1 << ((above_ref_length_u >> 1) + 1)) - 1;
                mask_rem &= mask_size2;
                while (mask_rem) {
                        remaining++;
                        mask_rem >>= 1;
                }
                min_size -= remaining << 1;
        }

        // ref above
        if (!up_available) {
                lm_params.offset_cb = 512;
                lm_params.offset_cr = 512;
                lm_params.scale_cb = 0;
                lm_params.scale_cr = 0;
                lm_params.shift_cb = 0;
                lm_params.shift_cr = 0;
        } else {
                if (left_available) {
                        const uint8_t first_row_in_ctu = !y0;
                        const uint16_t *_src_cb, *_src_cr, *_src;
                        int step_above = OVMAX(1, min_size >> 2);
                        int num_sample_a = OVMIN(min_size, 4); // = 2
                        int offset_pos_a = min_size >> 3;

                        if (first_row_in_ctu) {
                                _src = src_luma - src_stride;
                                _src_cb = dst_cb - src_stride;
                                _src_cr = dst_cr - src_stride;
                                for (i = 0, i_s = offset_pos_a << W_SHIFT;
                                     i < num_sample_a;
                                     ++i, i_s += step_above << W_SHIFT) {
                                        lm_ref[i] =
                                          (_src[i_s - 1] + (_src[i_s] << 1) +
                                           _src[i_s + 1] + 2) >>
                                          2;
                                        cb_ref[i] = _src_cb[i_s >> W_SHIFT];
                                        cr_ref[i] = _src_cr[i_s >> W_SHIFT];
                                }
                        } else {
                                _src = src_luma - (src_stride << H_SHIFT);
                                _src_cb = dst_cb - src_stride;
                                _src_cr = dst_cr - src_stride;
                                for (i = 0, i_s = offset_pos_a << W_SHIFT;
                                     i < num_sample_a;
                                     ++i, i_s += step_above << W_SHIFT) {
                                        lm_ref[i] =
                                          (_src[i_s - src_stride] +
                                           _src[i_s - 1] + (_src[i_s] << 2) +
                                           _src[i_s + 1] +
                                           _src[i_s + src_stride] + 4) >>
                                          3;
                                        cb_ref[i] = _src_cb[i_s >> W_SHIFT];
                                        cr_ref[i] = _src_cr[i_s >> W_SHIFT];
                                }
                        }
                } else {
                        const uint8_t first_row_in_ctu = !y0;
                        const uint16_t *_src_cb, *_src_cr, *_src;
                        int step_above = OVMAX(1, min_size >> 2);
                        int num_sample = OVMIN(min_size, 4); // = 2
                        int offset_a = min_size >> 3;

                        // FIXME positions
                        if (first_row_in_ctu) {
                                _src = src_luma - src_stride;
                                _src_cb = dst_cb - src_stride;
                                _src_cr = dst_cr - src_stride;
                                if (!offset_a) {
                                        lm_ref[0] = _src[0];
                                        cb_ref[0] = _src_cb[0];
                                        cr_ref[0] = _src_cr[0];
                                        _src += step_above << W_SHIFT;
                                        _src_cb += step_above;
                                        _src_cr += step_above;
                                }
                                for (i = !offset_a, i_s = offset_a << W_SHIFT;
                                     i < num_sample;
                                     ++i, i_s += step_above << W_SHIFT) {
                                        lm_ref[i] =
                                          (_src[i_s - 1] + (_src[i_s] << 1) +
                                           _src[i_s + 1] + 2) >>
                                          2;
                                        cb_ref[i] = _src_cb[i_s >> W_SHIFT];
                                        cr_ref[i] = _src_cr[i_s >> W_SHIFT];
                                }
                        } else {
                                _src = src_luma - (src_stride << H_SHIFT);
                                _src_cb = dst_cb - src_stride;
                                _src_cr = dst_cr - src_stride;
                                if (!offset_a) {
                                        lm_ref[0] =
                                          (_src[0 - src_stride] +
                                           (_src[0] << 1) +
                                           _src[0 + src_stride] + 2) >>
                                          2;
                                        cb_ref[0] = _src_cb[0];
                                        cr_ref[0] = _src_cr[0];
                                        _src += step_above << W_SHIFT;
                                        _src_cb += step_above;
                                        _src_cr += step_above;
                                }
                                for (i = !offset_a, i_s = offset_a << W_SHIFT;
                                     i < num_sample;
                                     ++i, i_s += step_above << W_SHIFT) {
                                        lm_ref[i] =
                                          (_src[i_s - src_stride] +
                                           _src[i_s - 1] + (_src[i_s] << 2) +
                                           _src[i_s + 1] +
                                           _src[i_s + src_stride] + 4) >>
                                          3;
                                        cb_ref[i] = _src_cb[i_s >> W_SHIFT];
                                        cr_ref[i] = _src_cr[i_s >> W_SHIFT];
                                }
                        }
                }

                if (OVMIN(min_size, 4) == 2) {
                        lm_ref[3] = lm_ref[0];
                        cb_ref[3] = cb_ref[0];
                        cr_ref[3] = cr_ref[0];
                        lm_ref[2] = lm_ref[1];
                        cb_ref[2] = cb_ref[1];
                        cr_ref[2] = cr_ref[1];
                        lm_ref[0] = lm_ref[1];
                        cb_ref[0] = cb_ref[1];
                        cr_ref[0] = cr_ref[1];
                        lm_ref[1] = lm_ref[3];
                        cb_ref[1] = cb_ref[3];
                        cr_ref[1] = cr_ref[3];
                }

                // TODO compute lm_param for cb and cr components
                derive_lm_parameters(lm_ref, cb_ref, cr_ref, &lm_params);
        }

        // inner part from reconstructed picture buffer
        // TODO directly scale cb and cr components;
        if (left_available) {
                int32_t value;
                uint16_t *_dst_cb, *_dst_cr;
                const uint16_t* _src;
                int offset_cb = lm_params.offset_cb;
                int offset_cr = lm_params.offset_cr;
                int scale_cb = lm_params.scale_cb;
                int scale_cr = lm_params.scale_cr;
                int shift_cb = lm_params.shift_cb;
                int shift_cr = lm_params.shift_cr;
                _dst_cb = dst_cb;
                _dst_cr = dst_cr;
                _src = src_luma;
                for (j = 0; j < 1 << log2_pb_h; ++j) {
                        for (i = 0, i_s = 0; i < 1 << log2_pb_w;
                             ++i, i_s = i << W_SHIFT) {
                                value =
                                  (_src[i_s - src_stride] + _src[i_s - 1] +
                                   (_src[i_s] << 2) + _src[i_s + 1] +
                                   _src[i_s + src_stride] + 4) >>
                                  3;
                                _dst_cb[i] = ov_clip(
                                  ((value * scale_cb) >> shift_cb) + offset_cb,
                                  0,
                                  1023);
                                _dst_cr[i] = ov_clip(
                                  ((value * scale_cr) >> shift_cr) + offset_cr,
                                  0,
                                  1023);
                        }
                        _dst_cb += dst_stride;
                        _dst_cr += dst_stride;
                        _src += src_stride << H_SHIFT;
                }
        } else {
                int32_t value;
                uint16_t *_dst_cb, *_dst_cr;
                const uint16_t* _src;
                int offset_cb = lm_params.offset_cb;
                int offset_cr = lm_params.offset_cr;
                int scale_cb = lm_params.scale_cb;
                int scale_cr = lm_params.scale_cr;
                int shift_cb = lm_params.shift_cb;
                int shift_cr = lm_params.shift_cr;
                _dst_cb = dst_cb;
                _dst_cr = dst_cr;
                _src = src_luma;
                for (j = 0; j < 1 << log2_pb_h; ++j) {
                        value = (_src[0 - src_stride] + (_src[0] << 1) +
                                 _src[0 + src_stride] + 2) >>
                                2;
                        _dst_cb[0] =
                          ov_clip(((value * scale_cb) >> shift_cb) + offset_cb,
                                  0,
                                  1023);
                        _dst_cr[0] =
                          ov_clip(((value * scale_cr) >> shift_cr) + offset_cr,
                                  0,
                                  1023);

                        for (i = 1, i_s = 1 << W_SHIFT; i < 1 << log2_pb_w;
                             ++i, i_s = i << W_SHIFT) {
                                value =
                                  (_src[i_s - src_stride] + _src[i_s - 1] +
                                   (_src[i_s] << 2) + _src[i_s + 1] +
                                   _src[i_s + src_stride] + 4) >>
                                  3;
                                _dst_cb[i] = ov_clip(
                                  ((value * scale_cb) >> shift_cb) + offset_cb,
                                  0,
                                  1023);
                                _dst_cr[i] = ov_clip(
                                  ((value * scale_cr) >> shift_cr) + offset_cr,
                                  0,
                                  1023);
                        }
                        _dst_cb += dst_stride;
                        _dst_cr += dst_stride;
                        _src += src_stride << H_SHIFT;
                }
        }
}

static void
vvc_intra_mdlm_left_cl(const uint16_t* const src_luma, uint16_t* const dst_cb,
                       uint16_t* const dst_cr, uint64_t intra_map_cols,
                       int log2_pb_w, int log2_pb_h, int x0, int y0,
                       uint8_t left_available, uint8_t up_available)
{ // left and above
        int i, j, i_s;

        uint16_t lm_ref[4];
        uint16_t cb_ref[4];
        uint16_t cr_ref[4];
        VVCLMParams lm_params;
        const int src_stride = RCN_CTB_STRIDE;
        const int dst_stride = RCN_CTB_STRIDE;

        int y0_u = y0 >> 1;

        int left_ref_length_u =
          /*(*/ (1 << (log2_pb_h /* + 1*/)) >> 1 /* + (1 << 1) - 1) >> 1*/;
        left_ref_length_u +=
          OVMIN((1 << log2_pb_w) >> 1, (1 << (log2_pb_h /* + 1*/)) >> 1);
        unsigned int mask_shift = 63 - (y0_u + left_ref_length_u);
        // (no need for top_left ref here so +1 instead of +2 for mask size)
        uint64_t mask_size = ((uint64_t)1 << (left_ref_length_u + 1)) - 1;
        uint64_t needed_mask = mask_size << mask_shift;
        uint64_t usable_mask = intra_map_cols & needed_mask;

        int min_log2_size;
        int min_size;
        if (!(usable_mask ^ needed_mask)) {
                min_size = left_ref_length_u << 1;
        } else {
                min_size = left_ref_length_u << 1;
                int remaining = 0;
                uint64_t mask_rem = (usable_mask ^ needed_mask) >> mask_shift;
                uint64_t mask_size2 =
                  ((uint64_t)1 << ((left_ref_length_u >> 1) + 1)) - 1;
                mask_rem &= mask_size2;
                while (mask_rem) {
                        remaining++;
                        mask_rem >>= 1;
                }
                min_size -= remaining << 1;
        }

        // ref above
        if (!left_available) {
                lm_params.offset_cb = 512;
                lm_params.offset_cr = 512;
                lm_params.scale_cb = 0;
                lm_params.scale_cr = 0;
                lm_params.shift_cb = 0;
                lm_params.shift_cr = 0;
        } else {
                if (up_available) {
                        const uint8_t first_row_in_ctu = !y0;
                        const uint16_t *_src_cb, *_src_cr, *_src;
                        int step_left = OVMAX(1, min_size >> 2);
                        int num_sample_l = OVMIN(min_size, 4); // = 2
                        int offset_pos_l = min_size >> 3;

                        // ref_left
                        // FIXME check - 2 - W_SHIFT
                        _src = src_luma - 2 - W_SHIFT +
                               (offset_pos_l << H_SHIFT) * src_stride;
                        _src_cb = dst_cb - 1 + offset_pos_l * src_stride;
                        _src_cr = dst_cr - 1 + offset_pos_l * src_stride;

                        for (j = 0; j < num_sample_l; ++j) {
                                lm_ref[j] = (_src[1 - src_stride] + _src[0] +
                                             (_src[1] << 2) + _src[2] +
                                             _src[1 + src_stride] + 4) >>
                                            3;
                                cb_ref[j] = _src_cb[0];
                                cr_ref[j] = _src_cr[0];
                                _src += (src_stride * step_left) << H_SHIFT;
                                _src_cb += (src_stride * step_left);
                                _src_cr += (src_stride * step_left);
                        }
                } else {
                        const uint16_t *_src_cb, *_src_cr, *_src;
                        int step_left = OVMAX(1, min_size >> 2);
                        int num_sample_l = OVMIN(min_size, 4); // = 2
                        int offset_pos_l = min_size >> 3;

                        _src = src_luma - 2 - W_SHIFT +
                               (offset_pos_l << H_SHIFT) * src_stride;
                        _src_cb = dst_cb - 1 + offset_pos_l * src_stride;
                        _src_cr = dst_cr - 1 + offset_pos_l * src_stride;

                        if (!offset_pos_l) {
                                // FIXME first position
                                _src = src_luma - 2 - W_SHIFT;

                                lm_ref[0] =
                                  (_src[0] + (_src[1] << 1) + _src[2] + 2) >> 2;
                                cb_ref[0] = _src_cb[0];
                                cr_ref[0] = _src_cr[0];

                                _src += (src_stride * step_left) << H_SHIFT;
                                _src_cb += src_stride * step_left;
                                _src_cr += src_stride * step_left;
                        }

                        for (j = !offset_pos_l; j < num_sample_l; ++j) {
                                lm_ref[j] = (_src[1 - src_stride] + _src[0] +
                                             (_src[1] << 2) + _src[2] +
                                             _src[1 + src_stride] + 4) >>
                                            3;
                                cb_ref[j] = _src_cb[0];
                                cr_ref[j] = _src_cr[0];
                                _src += (src_stride * step_left) << H_SHIFT;
                                _src_cb += src_stride * step_left;
                                _src_cr += src_stride * step_left;
                        }
                }

                if (OVMIN(min_size, 4) == 2) {
                        lm_ref[3] = lm_ref[0];
                        cb_ref[3] = cb_ref[0];
                        cr_ref[3] = cr_ref[0];
                        lm_ref[2] = lm_ref[1];
                        cb_ref[2] = cb_ref[1];
                        cr_ref[2] = cr_ref[1];
                        lm_ref[0] = lm_ref[1];
                        cb_ref[0] = cb_ref[1];
                        cr_ref[0] = cr_ref[1];
                        lm_ref[1] = lm_ref[3];
                        cb_ref[1] = cb_ref[3];
                        cr_ref[1] = cr_ref[3];
                }

                // TODO compute lm_param for cb and cr components
                derive_lm_parameters(lm_ref, cb_ref, cr_ref, &lm_params);
        }

        // inner part from reconstructed picture buffer
        // TODO directly scale cb and cr components;
        if (up_available) {
                int32_t value;
                uint16_t *_dst_cb, *_dst_cr;
                const uint16_t* _src;
                int offset_cb = lm_params.offset_cb;
                int offset_cr = lm_params.offset_cr;
                int scale_cb = lm_params.scale_cb;
                int scale_cr = lm_params.scale_cr;
                int shift_cb = lm_params.shift_cb;
                int shift_cr = lm_params.shift_cr;
                _dst_cb = dst_cb;
                _dst_cr = dst_cr;
                _src = src_luma;
                for (j = 0; j < 1 << log2_pb_h; ++j) {
                        for (i = 0, i_s = 0; i < 1 << log2_pb_w;
                             ++i, i_s = i << W_SHIFT) {
                                value =
                                  (_src[i_s - src_stride] + _src[i_s - 1] +
                                   (_src[i_s] << 2) + _src[i_s + 1] +
                                   _src[i_s + src_stride] + 4) >>
                                  3;
                                _dst_cb[i] = ov_clip(
                                  ((value * scale_cb) >> shift_cb) + offset_cb,
                                  0,
                                  1023);
                                _dst_cr[i] = ov_clip(
                                  ((value * scale_cr) >> shift_cr) + offset_cr,
                                  0,
                                  1023);
                        }
                        _dst_cb += dst_stride;
                        _dst_cr += dst_stride;
                        _src += src_stride << H_SHIFT;
                }
        } else {
                int32_t value;
                uint16_t *_dst_cb, *_dst_cr;
                const uint16_t* _src;
                int offset_cb = lm_params.offset_cb;
                int offset_cr = lm_params.offset_cr;
                int scale_cb = lm_params.scale_cb;
                int scale_cr = lm_params.scale_cr;
                int shift_cb = lm_params.shift_cb;
                int shift_cr = lm_params.shift_cr;
                _dst_cb = dst_cb;
                _dst_cr = dst_cr;
                _src = src_luma;
                for (i = 0, i_s = 0; i < 1 << log2_pb_w;
                     ++i, i_s = i << W_SHIFT) {
                        value = (_src[i_s - 1] + (_src[i_s] << 1) +
                                 _src[i_s + 1] + 2) >>
                                2;
                        _dst_cb[i] =
                          ov_clip(((value * scale_cb) >> shift_cb) + offset_cb,
                                  0,
                                  1023);
                        _dst_cr[i] =
                          ov_clip(((value * scale_cr) >> shift_cr) + offset_cr,
                                  0,
                                  1023);
                }

                _dst_cb += dst_stride;
                _dst_cr += dst_stride;
                _src += src_stride << H_SHIFT;

                for (j = 1; j < 1 << log2_pb_h; ++j) {
                        for (i = 0, i_s = 0; i < 1 << log2_pb_w;
                             ++i, i_s = i << W_SHIFT) {
                                value =
                                  (_src[i_s - src_stride] + _src[i_s - 1] +
                                   (_src[i_s] << 2) + _src[i_s + 1] +
                                   _src[i_s + src_stride] + 4) >>
                                  3;
                                _dst_cb[i] = ov_clip(
                                  ((value * scale_cb) >> shift_cb) + offset_cb,
                                  0,
                                  1023);
                                _dst_cr[i] = ov_clip(
                                  ((value * scale_cr) >> shift_cr) + offset_cr,
                                  0,
                                  1023);
                        }
                        _dst_cb += dst_stride;
                        _dst_cr += dst_stride;
                        _src += src_stride << H_SHIFT;
                }
        }
}
static void
vvc_intra_cclm(const uint16_t* const src_luma, uint16_t* const dst_cb,
               uint16_t* const dst_cr, int log2_pb_w, int log2_pb_h, int y0,
               int up_available, int left_available)
{ // left and above
        int i, j, i_s;

        uint16_t lm_ref[4];
        uint16_t cb_ref[4];
        uint16_t cr_ref[4];
        VVCLMParams lm_params;
        const int src_stride = RCN_CTB_STRIDE;
        const int dst_stride = RCN_CTB_STRIDE;

        // ref above
        if (up_available && left_available) {
                const uint8_t first_row_in_ctu = !y0;
                const uint16_t *_src_cb, *_src_cr, *_src;
                int step_above = OVMAX(1, (1 << log2_pb_w) >> 1);
                int step_left = OVMAX(1, (1 << log2_pb_h) >> 1);

                int num_sample_l = OVMIN(1 << log2_pb_h, 2); // = 2
                int num_sample_a = OVMIN(1 << log2_pb_w, 2); // = 2

                int offset_pos_a = (1 << log2_pb_w) >> 2;
                int offset_pos_l = (1 << log2_pb_h) >> 2;
                if (first_row_in_ctu) {
                        _src = src_luma - src_stride;
                        _src_cb = dst_cb - src_stride;
                        _src_cr = dst_cr - src_stride;
                        for (i = 0, i_s = offset_pos_a << W_SHIFT;
                             i < num_sample_a;
                             ++i, i_s += step_above << W_SHIFT) {
                                lm_ref[i] = (_src[i_s - 1] + (_src[i_s] << 1) +
                                             _src[i_s + 1] + 2) >>
                                            2;
                                cb_ref[i] = _src_cb[i_s >> W_SHIFT];
                                cr_ref[i] = _src_cr[i_s >> W_SHIFT];
                        }
                } else {
                        _src = src_luma - (src_stride << H_SHIFT);
                        _src_cb = dst_cb - src_stride;
                        _src_cr = dst_cr - src_stride;
                        for (i = 0, i_s = offset_pos_a << W_SHIFT;
                             i < num_sample_a;
                             ++i, i_s += step_above << W_SHIFT) {
                                lm_ref[i] =
                                  (_src[i_s - 1] + (_src[i_s] << 1) +
                                   _src[i_s + 1] + _src[i_s + src_stride - 1] +
                                   (_src[i_s + src_stride] << 1) +
                                   _src[i_s + src_stride + 1] + 4) >>
                                  3;
                                cb_ref[i] = _src_cb[i_s >> W_SHIFT];
                                cr_ref[i] = _src_cr[i_s >> W_SHIFT];
                        }
                }

                // ref_left
                // FIXME check - 2 - W_SHIFT
                _src = src_luma - 2 - W_SHIFT +
                       (offset_pos_l << H_SHIFT) * src_stride;
                _src_cb = dst_cb - 1 + offset_pos_l * src_stride;
                _src_cr = dst_cr - 1 + offset_pos_l * src_stride;

                for (j = num_sample_a; j < num_sample_l + num_sample_a; ++j) {
                        lm_ref[j] = (_src[1 - 1] + (_src[1] << 1) +
                                     _src[1 + 1] + _src[1 + src_stride - 1] +
                                     (_src[1 + src_stride] << 1) +
                                     _src[1 + src_stride + 1] + 4) >>
                                    3;
                        cb_ref[j] = _src_cb[0];
                        cr_ref[j] = _src_cr[0];
                        _src += (src_stride * step_left) << H_SHIFT;
                        _src_cb += (src_stride * step_left);
                        _src_cr += (src_stride * step_left);
                }
        } else if (left_available) {
                const uint16_t *_src_cb, *_src_cr, *_src;
                int step_left = OVMAX(1, (1 << log2_pb_h) >> 2);
                int num_sample = OVMIN(1 << log2_pb_h, 4);
                int offset_pos_l = (1 << log2_pb_h) >> 3;

                _src = src_luma - 2 - W_SHIFT +
                       (offset_pos_l << H_SHIFT) * src_stride;
                _src_cb = dst_cb - 1 + offset_pos_l * src_stride;
                _src_cr = dst_cr - 1 + offset_pos_l * src_stride;

                for (j = 0; j < num_sample; ++j) {
                        lm_ref[j] = (_src[1 - 1] + (_src[1] << 1) +
                                     _src[1 + 1] + _src[1 + src_stride - 1] +
                                     (_src[1 + src_stride] << 1) +
                                     _src[1 + src_stride + 1] + 4) >>
                                    3;
                        cb_ref[j] = _src_cb[0];
                        cr_ref[j] = _src_cr[0];
                        _src += (src_stride * step_left) << H_SHIFT;
                        _src_cb += src_stride * step_left;
                        _src_cr += src_stride * step_left;
                }
                if (num_sample == 2) {
                        lm_ref[3] = lm_ref[0];
                        cb_ref[3] = cb_ref[0];
                        cr_ref[3] = cr_ref[0];
                        lm_ref[2] = lm_ref[1];
                        cb_ref[2] = cb_ref[1];
                        cr_ref[2] = cr_ref[1];
                        lm_ref[0] = lm_ref[1];
                        cb_ref[0] = cb_ref[1];
                        cr_ref[0] = cr_ref[1];
                        lm_ref[1] = lm_ref[3];
                        cb_ref[1] = cb_ref[3];
                        cr_ref[1] = cr_ref[3];
                }
        } else if (up_available) {
                const uint8_t first_row_in_ctu = !y0;
                const uint16_t *_src_cb, *_src_cr, *_src;
                int step_above = OVMAX(1, (1 << log2_pb_w) >> 2);
                int num_sample = OVMIN(1 << log2_pb_w, 4);
                int offset_a = (1 << log2_pb_w) >> 3;

                // FIXME positions
                if (first_row_in_ctu) {
                        _src = src_luma - src_stride;
                        _src_cb = dst_cb - src_stride;
                        _src_cr = dst_cr - src_stride;
                        if (!offset_a) {
                                lm_ref[0] = _src[0];
                                cb_ref[0] = _src_cb[0];
                                cr_ref[0] = _src_cr[0];
                                _src += step_above << W_SHIFT;
                                _src_cb += step_above;
                                _src_cr += step_above;
                        }
                        for (i = !offset_a, i_s = offset_a << W_SHIFT;
                             i < num_sample;
                             ++i, i_s += step_above << W_SHIFT) {
                                lm_ref[i] = (_src[i_s - 1] + (_src[i_s] << 1) +
                                             _src[i_s + 1] + 2) >>
                                            2;
                                cb_ref[i] = _src_cb[i_s >> W_SHIFT];
                                cr_ref[i] = _src_cr[i_s >> W_SHIFT];
                        }
                } else {
                        _src = src_luma - (src_stride << H_SHIFT);
                        _src_cb = dst_cb - src_stride;
                        _src_cr = dst_cr - src_stride;

                        if (!offset_a) {
                                lm_ref[0] =
                                  (_src[0] + _src[0 + src_stride] + 1) >> 1;
                                cb_ref[0] = _src_cb[0];
                                cr_ref[0] = _src_cr[0];
                                _src += step_above << W_SHIFT;
                                _src_cb += step_above;
                                _src_cr += step_above;
                        }

                        for (i = !offset_a, i_s = offset_a << W_SHIFT;
                             i < num_sample;
                             ++i, i_s += step_above << W_SHIFT) {
                                lm_ref[i] =
                                  (_src[i_s - 1] + (_src[i_s] << 1) +
                                   _src[i_s + 1] + _src[i_s + src_stride - 1] +
                                   (_src[i_s + src_stride] << 1) +
                                   _src[i_s + src_stride + 1] + 4) >>
                                  3;
                                cb_ref[i] = _src_cb[i_s >> W_SHIFT];
                                cr_ref[i] = _src_cr[i_s >> W_SHIFT];
                        }
                }
                if (num_sample == 2) {
                        lm_ref[3] = lm_ref[0];
                        cb_ref[3] = cb_ref[0];
                        cr_ref[3] = cr_ref[0];
                        lm_ref[2] = lm_ref[1];
                        cb_ref[2] = cb_ref[1];
                        cr_ref[2] = cr_ref[1];
                        lm_ref[0] = lm_ref[1];
                        cb_ref[0] = cb_ref[1];
                        cr_ref[0] = cr_ref[1];
                        lm_ref[1] = lm_ref[3];
                        cb_ref[1] = cb_ref[3];
                        cr_ref[1] = cr_ref[3];
                }
        } else {
                lm_ref[0] = 512;
                lm_ref[1] = 512;
                lm_ref[2] = 512;
                lm_ref[3] = 512;
                cb_ref[0] = 512;
                cb_ref[1] = 512;
                cb_ref[2] = 512;
                cb_ref[3] = 512;
                cr_ref[0] = 512;
                cr_ref[1] = 512;
                cr_ref[2] = 512;
                cr_ref[3] = 512;
        }

        // TODO compute lm_param for cb and cr components
        derive_lm_parameters(lm_ref, cb_ref, cr_ref, &lm_params);

        // inner part from reconstructed picture buffer
        // TODO directly scale cb and cr components;
        if (up_available && left_available) {
                int32_t value;
                uint16_t *_dst_cb, *_dst_cr;
                const uint16_t* _src;
                int offset_cb = lm_params.offset_cb;
                int offset_cr = lm_params.offset_cr;
                int scale_cb = lm_params.scale_cb;
                int scale_cr = lm_params.scale_cr;
                int shift_cb = lm_params.shift_cb;
                int shift_cr = lm_params.shift_cr;
                _dst_cb = dst_cb;
                _dst_cr = dst_cr;
                _src = src_luma;
                for (j = 0; j < 1 << log2_pb_h; ++j) {
                        for (i = 0, i_s = 0; i < 1 << log2_pb_w;
                             ++i, i_s = i << W_SHIFT) {
                                value =
                                  (_src[i_s - 1] + (_src[i_s] << 1) +
                                   _src[i_s + 1] + _src[i_s + src_stride - 1] +
                                   (_src[i_s + src_stride] << 1) +
                                   _src[i_s + src_stride + 1] + 4) >>
                                  3;
                                _dst_cb[i] = ov_clip(
                                  ((value * scale_cb) >> shift_cb) + offset_cb,
                                  0,
                                  1023);
                                _dst_cr[i] = ov_clip(
                                  ((value * scale_cr) >> shift_cr) + offset_cr,
                                  0,
                                  1023);
                        }
                        _dst_cb += dst_stride;
                        _dst_cr += dst_stride;
                        _src += src_stride << H_SHIFT;
                }
        } else if (left_available) {
                int32_t value;
                uint16_t *_dst_cb, *_dst_cr;
                const uint16_t* _src;
                int offset_cb = lm_params.offset_cb;
                int offset_cr = lm_params.offset_cr;
                int scale_cb = lm_params.scale_cb;
                int scale_cr = lm_params.scale_cr;
                int shift_cb = lm_params.shift_cb;
                int shift_cr = lm_params.shift_cr;
                _dst_cb = dst_cb;
                _dst_cr = dst_cr;
                _src = src_luma;

                for (j = 0; j < 1 << log2_pb_h; ++j) {
                        for (i = 0, i_s = 0; i < 1 << log2_pb_w;
                             ++i, i_s = i << W_SHIFT) {
                                value =
                                  (_src[i_s - 1] + (_src[i_s] << 1) +
                                   _src[i_s + 1] + _src[i_s + src_stride - 1] +
                                   (_src[i_s + src_stride] << 1) +
                                   _src[i_s + src_stride + 1] + 4) >>
                                  3;
                                _dst_cb[i] = ov_clip(
                                  ((value * scale_cb) >> shift_cb) + offset_cb,
                                  0,
                                  1023);
                                _dst_cr[i] = ov_clip(
                                  ((value * scale_cr) >> shift_cr) + offset_cr,
                                  0,
                                  1023);
                        }
                        _dst_cb += dst_stride;
                        _dst_cr += dst_stride;
                        _src += src_stride << H_SHIFT;
                }
        } else if (up_available) {
                int32_t value;
                uint16_t *_dst_cb, *_dst_cr;
                const uint16_t* _src;
                int offset_cb = lm_params.offset_cb;
                int offset_cr = lm_params.offset_cr;
                int scale_cb = lm_params.scale_cb;
                int scale_cr = lm_params.scale_cr;
                int shift_cb = lm_params.shift_cb;
                int shift_cr = lm_params.shift_cr;
                _dst_cb = dst_cb;
                _dst_cr = dst_cr;
                _src = src_luma;
                for (j = 0; j < 1 << log2_pb_h; ++j) {
                        value = (_src[0] + _src[0 + src_stride] + 1) >> 1;
                        _dst_cb[0] =
                          ov_clip(((value * scale_cb) >> shift_cb) + offset_cb,
                                  0,
                                  1023);
                        _dst_cr[0] =
                          ov_clip(((value * scale_cr) >> shift_cr) + offset_cr,
                                  0,
                                  1023);

                        for (i = 1, i_s = 1 << W_SHIFT; i < 1 << log2_pb_w;
                             ++i, i_s = i << W_SHIFT) {
                                value =
                                  (_src[i_s - 1] + (_src[i_s] << 1) +
                                   _src[i_s + 1] + _src[i_s + src_stride - 1] +
                                   (_src[i_s + src_stride] << 1) +
                                   _src[i_s + src_stride + 1] + 4) >>
                                  3;
                                _dst_cb[i] = ov_clip(
                                  ((value * scale_cb) >> shift_cb) + offset_cb,
                                  0,
                                  1023);
                                _dst_cr[i] = ov_clip(
                                  ((value * scale_cr) >> shift_cr) + offset_cr,
                                  0,
                                  1023);
                        }
                        _dst_cb += dst_stride;
                        _dst_cr += dst_stride;
                        _src += src_stride << H_SHIFT;
                }
        } else {
                int32_t value;
                uint16_t *_dst_cb, *_dst_cr;
                const uint16_t* _src;
                int offset_cb = 512;
                int offset_cr = 512;
                int scale_cb = 0;
                int scale_cr = 0;
                int shift_cb = 0;
                int shift_cr = 0;
                _dst_cb = dst_cb;
                _dst_cr = dst_cr;
                _src = src_luma;
                for (j = 0; j < 1 << log2_pb_h; ++j) {
                        value = (_src[0] + _src[0 + src_stride] + 1) >> 1;
                        _dst_cb[0] =
                          ov_clip(((value * scale_cb) >> shift_cb) + offset_cb,
                                  0,
                                  1023);
                        _dst_cr[0] =
                          ov_clip(((value * scale_cr) >> shift_cr) + offset_cr,
                                  0,
                                  1023);

                        for (i = 1, i_s = 1 << W_SHIFT; i < 1 << log2_pb_w;
                             ++i, i_s = i << W_SHIFT) {
                                value =
                                  (_src[i_s - 1] + (_src[i_s] << 1) +
                                   _src[i_s + 1] + _src[i_s + src_stride - 1] +
                                   (_src[i_s + src_stride] << 1) +
                                   _src[i_s + src_stride + 1] + 4) >>
                                  3;
                                _dst_cb[i] = ov_clip(
                                  ((value * scale_cb) >> shift_cb) + offset_cb,
                                  0,
                                  1023);
                                _dst_cr[i] = ov_clip(
                                  ((value * scale_cr) >> shift_cr) + offset_cr,
                                  0,
                                  1023);
                        }
                        _dst_cb += dst_stride;
                        _dst_cr += dst_stride;
                        _src += src_stride << H_SHIFT;
                }
        }
}

static void
vvc_intra_mdlm_top(const uint16_t* const src_luma, uint16_t* const dst_cb,
                   uint16_t* const dst_cr, uint64_t intra_map_rows,
                   int log2_pb_w, int log2_pb_h, int x0, int y0,
                   uint8_t left_available, uint8_t up_available)
{ // left and above
        int i, j, i_s;

        uint16_t lm_ref[4];
        uint16_t cb_ref[4];
        uint16_t cr_ref[4];
        VVCLMParams lm_params;
        const int src_stride = RCN_CTB_STRIDE;
        const int dst_stride = RCN_CTB_STRIDE;

        int x0_u = x0 >> 1;
        int above_ref_length_u = (1 << log2_pb_w) >> 1;

        above_ref_length_u +=
          OVMIN((1 << log2_pb_h) >> 1, (1 << log2_pb_w) >> 1);
        unsigned int mask_shift = 63 - (x0_u + above_ref_length_u);

        uint64_t mask_size = ((uint64_t)1 << (above_ref_length_u + 1)) - 1;
        uint64_t needed_mask = mask_size << mask_shift;
        uint64_t usable_mask = intra_map_rows & needed_mask;

        // We use PU width here since this won't change lm paramters
        // int sum_luma = 0, sum_chroma = 0, sum_xx = 0, sum_xy = 0;
        int min_log2_size;
        int min_size;
        if (!(usable_mask ^ needed_mask)) {
                min_size = above_ref_length_u << 1;
        } else {
                min_size = above_ref_length_u << 1;
                int remaining = 0;
                uint64_t mask_rem = (usable_mask ^ needed_mask) >> mask_shift;
                uint64_t mask_size2 =
                  ((uint64_t)1 << ((above_ref_length_u >> 1) + 1)) - 1;
                mask_rem &= mask_size2;
                while (mask_rem) {
                        remaining++;
                        mask_rem >>= 1;
                }
                min_size -= remaining << 1;
        }

        // ref above
        if (!up_available) {
                lm_params.offset_cb = 512;
                lm_params.offset_cr = 512;
                lm_params.scale_cb = 0;
                lm_params.scale_cr = 0;
                lm_params.shift_cb = 0;
                lm_params.shift_cr = 0;
        } else {
                if (left_available) {
                        const uint8_t first_row_in_ctu = !y0;
                        const uint16_t *_src_cb, *_src_cr, *_src;
                        int step_above = OVMAX(1, min_size >> 2);
                        int num_sample_a = OVMIN(min_size, 4); // = 2
                        int offset_pos_a = min_size >> 3;

                        if (first_row_in_ctu) {
                                _src = src_luma - src_stride;
                                _src_cb = dst_cb - src_stride;
                                _src_cr = dst_cr - src_stride;
                                for (i = 0, i_s = offset_pos_a << W_SHIFT;
                                     i < num_sample_a;
                                     ++i, i_s += step_above << W_SHIFT) {
                                        lm_ref[i] =
                                          (_src[i_s - 1] + (_src[i_s] << 1) +
                                           _src[i_s + 1] + 2) >>
                                          2;
                                        cb_ref[i] = _src_cb[i_s >> W_SHIFT];
                                        cr_ref[i] = _src_cr[i_s >> W_SHIFT];
                                }
                        } else {
                                _src = src_luma - (src_stride << H_SHIFT);
                                _src_cb = dst_cb - src_stride;
                                _src_cr = dst_cr - src_stride;
                                for (i = 0, i_s = offset_pos_a << W_SHIFT;
                                     i < num_sample_a;
                                     ++i, i_s += step_above << W_SHIFT) {
                                        lm_ref[i] =
                                          (_src[i_s - 1] + (_src[i_s] << 1) +
                                           _src[i_s + 1] +
                                           _src[i_s + src_stride - 1] +
                                           (_src[i_s + src_stride] << 1) +
                                           _src[i_s + src_stride + 1] + 4) >>
                                          3;
                                        cb_ref[i] = _src_cb[i_s >> W_SHIFT];
                                        cr_ref[i] = _src_cr[i_s >> W_SHIFT];
                                }
                        }
                } else {
                        const uint8_t first_row_in_ctu = !y0;
                        const uint16_t *_src_cb, *_src_cr, *_src;
                        int step_above = OVMAX(1, min_size >> 2);
                        int num_sample = OVMIN(min_size, 4); // = 2
                        int offset_a = min_size >> 3;

                        // FIXME positions
                        if (first_row_in_ctu) {
                                _src = src_luma - src_stride;
                                _src_cb = dst_cb - src_stride;
                                _src_cr = dst_cr - src_stride;
                                if (!offset_a) {
                                        lm_ref[0] = _src[0];
                                        cb_ref[0] = _src_cb[0];
                                        cr_ref[0] = _src_cr[0];
                                        _src += step_above << W_SHIFT;
                                        _src_cb += step_above;
                                        _src_cr += step_above;
                                }
                                for (i = !offset_a, i_s = offset_a << W_SHIFT;
                                     i < num_sample;
                                     ++i, i_s += step_above << W_SHIFT) {
                                        lm_ref[i] =
                                          (_src[i_s - 1] + (_src[i_s] << 1) +
                                           _src[i_s + 1] + 2) >>
                                          2;
                                        cb_ref[i] = _src_cb[i_s >> W_SHIFT];
                                        cr_ref[i] = _src_cr[i_s >> W_SHIFT];
                                }
                        } else {
                                _src = src_luma - (src_stride << H_SHIFT);
                                _src_cb = dst_cb - src_stride;
                                _src_cr = dst_cr - src_stride;
                                if (!offset_a) {
                                        lm_ref[0] =
                                          (_src[0] + _src[0 + src_stride] +
                                           1) >>
                                          1;
                                        cb_ref[0] = _src_cb[0];
                                        cr_ref[0] = _src_cr[0];
                                        _src += step_above << W_SHIFT;
                                        _src_cb += step_above;
                                        _src_cr += step_above;
                                }
                                for (i = !offset_a, i_s = offset_a << W_SHIFT;
                                     i < num_sample;
                                     ++i, i_s += step_above << W_SHIFT) {
                                        lm_ref[i] =
                                          (_src[i_s - 1] + (_src[i_s] << 1) +
                                           _src[i_s + 1] +
                                           _src[i_s + src_stride - 1] +
                                           (_src[i_s + src_stride] << 1) +
                                           _src[i_s + src_stride + 1] + 4) >>
                                          3;
                                        cb_ref[i] = _src_cb[i_s >> W_SHIFT];
                                        cr_ref[i] = _src_cr[i_s >> W_SHIFT];
                                }
                        }
                }

                if (OVMIN(min_size, 4) == 2) {
                        lm_ref[3] = lm_ref[0];
                        cb_ref[3] = cb_ref[0];
                        cr_ref[3] = cr_ref[0];
                        lm_ref[2] = lm_ref[1];
                        cb_ref[2] = cb_ref[1];
                        cr_ref[2] = cr_ref[1];
                        lm_ref[0] = lm_ref[1];
                        cb_ref[0] = cb_ref[1];
                        cr_ref[0] = cr_ref[1];
                        lm_ref[1] = lm_ref[3];
                        cb_ref[1] = cb_ref[3];
                        cr_ref[1] = cr_ref[3];
                }

                // TODO compute lm_param for cb and cr components
                derive_lm_parameters(lm_ref, cb_ref, cr_ref, &lm_params);
        }

        // inner part from reconstructed picture buffer
        // TODO directly scale cb and cr components;
        if (left_available) {
                int32_t value;
                uint16_t *_dst_cb, *_dst_cr;
                const uint16_t* _src;
                int offset_cb = lm_params.offset_cb;
                int offset_cr = lm_params.offset_cr;
                int scale_cb = lm_params.scale_cb;
                int scale_cr = lm_params.scale_cr;
                int shift_cb = lm_params.shift_cb;
                int shift_cr = lm_params.shift_cr;
                _dst_cb = dst_cb;
                _dst_cr = dst_cr;
                _src = src_luma;
                for (j = 0; j < 1 << log2_pb_h; ++j) {
                        for (i = 0, i_s = 0; i < 1 << log2_pb_w;
                             ++i, i_s = i << W_SHIFT) {
                                value =
                                  (_src[i_s - 1] + (_src[i_s] << 1) +
                                   _src[i_s + 1] + _src[i_s + src_stride - 1] +
                                   (_src[i_s + src_stride] << 1) +
                                   _src[i_s + src_stride + 1] + 4) >>
                                  3;
                                _dst_cb[i] = ov_clip(
                                  ((value * scale_cb) >> shift_cb) + offset_cb,
                                  0,
                                  1023);
                                _dst_cr[i] = ov_clip(
                                  ((value * scale_cr) >> shift_cr) + offset_cr,
                                  0,
                                  1023);
                        }
                        _dst_cb += dst_stride;
                        _dst_cr += dst_stride;
                        _src += src_stride << H_SHIFT;
                }
        } else {
                int32_t value;
                uint16_t *_dst_cb, *_dst_cr;
                const uint16_t* _src;
                int offset_cb = lm_params.offset_cb;
                int offset_cr = lm_params.offset_cr;
                int scale_cb = lm_params.scale_cb;
                int scale_cr = lm_params.scale_cr;
                int shift_cb = lm_params.shift_cb;
                int shift_cr = lm_params.shift_cr;
                _dst_cb = dst_cb;
                _dst_cr = dst_cr;
                _src = src_luma;
                for (j = 0; j < 1 << log2_pb_h; ++j) {
                        value = (_src[0] + _src[0 + src_stride] + 1) >> 1;
                        _dst_cb[0] =
                          ov_clip(((value * scale_cb) >> shift_cb) + offset_cb,
                                  0,
                                  1023);
                        _dst_cr[0] =
                          ov_clip(((value * scale_cr) >> shift_cr) + offset_cr,
                                  0,
                                  1023);

                        for (i = 1, i_s = 1 << W_SHIFT; i < 1 << log2_pb_w;
                             ++i, i_s = i << W_SHIFT) {
                                value =
                                  (_src[i_s - 1] + (_src[i_s] << 1) +
                                   _src[i_s + 1] + _src[i_s + src_stride - 1] +
                                   (_src[i_s + src_stride] << 1) +
                                   _src[i_s + src_stride + 1] + 4) >>
                                  3;
                                _dst_cb[i] = ov_clip(
                                  ((value * scale_cb) >> shift_cb) + offset_cb,
                                  0,
                                  1023);
                                _dst_cr[i] = ov_clip(
                                  ((value * scale_cr) >> shift_cr) + offset_cr,
                                  0,
                                  1023);
                        }
                        _dst_cb += dst_stride;
                        _dst_cr += dst_stride;
                        _src += src_stride << H_SHIFT;
                }
        }
}

static void
vvc_intra_mdlm_left(const uint16_t* const src_luma, uint16_t* const dst_cb,
                    uint16_t* const dst_cr, uint64_t intra_map_cols,
                    int log2_pb_w, int log2_pb_h, int x0, int y0,
                    uint8_t left_available, uint8_t up_available)
{ // left and above
        int i, j, i_s;

        uint16_t lm_ref[4];
        uint16_t cb_ref[4];
        uint16_t cr_ref[4];
        VVCLMParams lm_params;
        const int src_stride = RCN_CTB_STRIDE;
        const int dst_stride = RCN_CTB_STRIDE;

        int y0_u = y0 >> 1;

        int left_ref_length_u = (1 << (log2_pb_h)) >> 1;
        left_ref_length_u +=
          OVMIN((1 << log2_pb_w) >> 1, (1 << (log2_pb_h)) >> 1);
        unsigned int mask_shift = 63 - (y0_u + left_ref_length_u);
        // (no need for top_left ref here so +1 instead of +2 for mask size)
        uint64_t mask_size = ((uint64_t)1 << (left_ref_length_u + 1)) - 1;
        uint64_t needed_mask = mask_size << mask_shift;
        uint64_t usable_mask = intra_map_cols & needed_mask;

        int min_log2_size;
        int min_size;
        if (!(usable_mask ^ needed_mask)) {
                min_size = left_ref_length_u << 1;
        } else {
                min_size = left_ref_length_u << 1;
                int remaining = 0;
                uint64_t mask_rem = (usable_mask ^ needed_mask) >> mask_shift;
                uint64_t mask_size2 =
                  ((uint64_t)1 << ((left_ref_length_u >> 1) + 1)) - 1;
                mask_rem &= mask_size2;
                while (mask_rem) {
                        remaining++;
                        mask_rem >>= 1;
                }
                min_size -= remaining << 1;
        }

        // ref above
        if (!left_available) {
                lm_params.offset_cb = 512;
                lm_params.offset_cr = 512;
                lm_params.scale_cb = 0;
                lm_params.scale_cr = 0;
                lm_params.shift_cb = 0;
                lm_params.shift_cr = 0;
        } else {
                if (up_available) {
                        const uint8_t first_row_in_ctu = !y0;
                        const uint16_t *_src_cb, *_src_cr, *_src;
                        int step_left = OVMAX(1, min_size >> 2);
                        int num_sample_l = OVMIN(min_size, 4); // = 2
                        int offset_pos_l = min_size >> 3;

                        // ref_left
                        // FIXME check - 2 - W_SHIFT
                        _src = src_luma - 2 - W_SHIFT +
                               (offset_pos_l << H_SHIFT) * src_stride;
                        _src_cb = dst_cb - 1 + offset_pos_l * src_stride;
                        _src_cr = dst_cr - 1 + offset_pos_l * src_stride;

                        for (j = 0; j < num_sample_l; ++j) {
                                lm_ref[j] =
                                  (_src[1 - 1] + (_src[1] << 1) + _src[1 + 1] +
                                   _src[1 + src_stride - 1] +
                                   (_src[1 + src_stride] << 1) +
                                   _src[1 + src_stride + 1] + 4) >>
                                  3;
                                cb_ref[j] = _src_cb[0];
                                cr_ref[j] = _src_cr[0];
                                _src += (src_stride * step_left) << H_SHIFT;
                                _src_cb += (src_stride * step_left);
                                _src_cr += (src_stride * step_left);
                        }
                } else {
                        const uint16_t *_src_cb, *_src_cr, *_src;
                        int step_left = OVMAX(1, min_size >> 2);
                        int num_sample_l = OVMIN(min_size, 4); // = 2
                        int offset_pos_l = min_size >> 3;

                        _src = src_luma - 2 - W_SHIFT +
                               (offset_pos_l << H_SHIFT) * src_stride;
                        _src_cb = dst_cb - 1 + offset_pos_l * src_stride;
                        _src_cr = dst_cr - 1 + offset_pos_l * src_stride;

                        for (j = 0; j < num_sample_l; ++j) {
                                lm_ref[j] =
                                  (_src[1 - 1] + (_src[1] << 1) + _src[1 + 1] +
                                   _src[1 + src_stride - 1] +
                                   (_src[1 + src_stride] << 1) +
                                   _src[1 + src_stride + 1] + 4) >>
                                  3;
                                cb_ref[j] = _src_cb[0];
                                cr_ref[j] = _src_cr[0];
                                _src += (src_stride * step_left) << H_SHIFT;
                                _src_cb += src_stride * step_left;
                                _src_cr += src_stride * step_left;
                        }
                }

                if (OVMIN(min_size, 4) == 2) {
                        lm_ref[3] = lm_ref[0];
                        cb_ref[3] = cb_ref[0];
                        cr_ref[3] = cr_ref[0];
                        lm_ref[2] = lm_ref[1];
                        cb_ref[2] = cb_ref[1];
                        cr_ref[2] = cr_ref[1];
                        lm_ref[0] = lm_ref[1];
                        cb_ref[0] = cb_ref[1];
                        cr_ref[0] = cr_ref[1];
                        lm_ref[1] = lm_ref[3];
                        cb_ref[1] = cb_ref[3];
                        cr_ref[1] = cr_ref[3];
                }

                // TODO compute lm_param for cb and cr components
                derive_lm_parameters(lm_ref, cb_ref, cr_ref, &lm_params);
        }

        // inner part from reconstructed picture buffer
        // TODO directly scale cb and cr components;
        if (left_available) {
                int32_t value;
                uint16_t *_dst_cb, *_dst_cr;
                const uint16_t* _src;
                int offset_cb = lm_params.offset_cb;
                int offset_cr = lm_params.offset_cr;
                int scale_cb = lm_params.scale_cb;
                int scale_cr = lm_params.scale_cr;
                int shift_cb = lm_params.shift_cb;
                int shift_cr = lm_params.shift_cr;
                _dst_cb = dst_cb;
                _dst_cr = dst_cr;
                _src = src_luma;
                for (j = 0; j < 1 << log2_pb_h; ++j) {
                        for (i = 0, i_s = 0; i < 1 << log2_pb_w;
                             ++i, i_s = i << W_SHIFT) {
                                value =
                                  (_src[i_s - 1] + (_src[i_s] << 1) +
                                   _src[i_s + 1] + _src[i_s + src_stride - 1] +
                                   (_src[i_s + src_stride] << 1) +
                                   _src[i_s + src_stride + 1] + 4) >>
                                  3;
                                _dst_cb[i] = ov_clip(
                                  ((value * scale_cb) >> shift_cb) + offset_cb,
                                  0,
                                  1023);
                                _dst_cr[i] = ov_clip(
                                  ((value * scale_cr) >> shift_cr) + offset_cr,
                                  0,
                                  1023);
                        }
                        _dst_cb += dst_stride;
                        _dst_cr += dst_stride;
                        _src += src_stride << H_SHIFT;
                }
        } else {
                int32_t value;
                uint16_t *_dst_cb, *_dst_cr;
                const uint16_t* _src;
                int offset_cb = lm_params.offset_cb;
                int offset_cr = lm_params.offset_cr;
                int scale_cb = lm_params.scale_cb;
                int scale_cr = lm_params.scale_cr;
                int shift_cb = lm_params.shift_cb;
                int shift_cr = lm_params.shift_cr;
                _dst_cb = dst_cb;
                _dst_cr = dst_cr;
                _src = src_luma;

                for (j = 0; j < 1 << log2_pb_h; ++j) {
                        value = (_src[0] + _src[0 + src_stride] + 1) >> 1;
                        _dst_cb[0] =
                          ov_clip(((value * scale_cb) >> shift_cb) + offset_cb,
                                  0,
                                  1023);
                        _dst_cr[0] =
                          ov_clip(((value * scale_cr) >> shift_cr) + offset_cr,
                                  0,
                                  1023);
                        for (i = 1, i_s = 1 << W_SHIFT; i < 1 << log2_pb_w;
                             ++i, i_s = i << W_SHIFT) {
                                value =
                                  (_src[i_s - 1] + (_src[i_s] << 1) +
                                   _src[i_s + 1] + _src[i_s + src_stride - 1] +
                                   (_src[i_s + src_stride] << 1) +
                                   _src[i_s + src_stride + 1] + 4) >>
                                  3;
                                _dst_cb[i] = ov_clip(
                                  ((value * scale_cb) >> shift_cb) + offset_cb,
                                  0,
                                  1023);
                                _dst_cr[i] = ov_clip(
                                  ((value * scale_cr) >> shift_cr) + offset_cr,
                                  0,
                                  1023);
                        }
                        _dst_cb += dst_stride;
                        _dst_cr += dst_stride;
                        _src += src_stride << H_SHIFT;
                }
        }
}
#endif


struct LMParams{
   int shift;
   int a;
   int b;
};

struct CCLMParams{
   struct LMParams cb;
   struct LMParams cr;
};

struct AVGMinMax{
    uint16_t min_l;
    uint16_t max_l;
    uint16_t min_cb;
    uint16_t max_cb;
    uint16_t min_cr;
    uint16_t max_cr;
};

static inline struct LMParams compute_lm_params(int16_t avg_min_l, int16_t avg_min_c, int16_t avg_max_c,
                                                int16_t v, int8_t log2_rng_l);

static inline struct CCLMParams derive_cclm_params(const struct AVGMinMax *const avg);

static inline struct AVGMinMax sort_average_lm_ref_samples(const uint16_t *const lm_spm, const uint16_t *const smp_cb,
                                                    const uint16_t *const smp_cr, int nb_samples);

static void sub_sample_lm_ref_lft(const uint16_t *lm_src, const uint16_t *src_cb, const uint16_t *src_cr,
                                  uint16_t *const lm_dst, uint16_t *dst_cb, uint16_t *dst_cr,
                                  ptrdiff_t lm_src_stride, ptrdiff_t src_c_stride,
                                  int lft_step, int nb_sample_lft);

static void sub_sample_lm_ref_abv0(const uint16_t *lm_src, const uint16_t *src_cb, const uint16_t *src_cr,
                                   uint16_t *const lm_dst, uint16_t *dst_cb, uint16_t *dst_cr,
                                   ptrdiff_t lm_src_stride, ptrdiff_t src_c_stride,
                                   int lft_step, int nb_sample_abv, uint8_t lft_avail);

static void sub_sample_lm_ref_abv(const uint16_t *lm_src, const uint16_t *src_cb, const uint16_t *src_cr,
                                  uint16_t *const lm_dst, uint16_t *dst_cb, uint16_t *dst_cr,
                                  ptrdiff_t lm_src_stride, ptrdiff_t src_c_stride,
                                  int abv_step, int nb_sample_abv, uint8_t lft_avail);



static void compute_lm_subsample(const uint16_t *lm_src, uint16_t *dst_cb, uint16_t *dst_cr,
                                 ptrdiff_t lm_src_stride, ptrdiff_t dst_stride_c,
                                 const struct CCLMParams *const lm_params,
                                 int pb_w, int pb_h, uint8_t lft_avail);

static void sub_sample_lm_ref_lft_collocated(const uint16_t *lm_src, const uint16_t *src_cb, const uint16_t *src_cr,
                                  uint16_t *const lm_dst, uint16_t *dst_cb, uint16_t *dst_cr,
                                  ptrdiff_t lm_src_stride, ptrdiff_t src_c_stride,
                                  int lft_step, int nb_sample_lft, int abv_avail);


static void sub_sample_lm_ref_abv_collocated(const uint16_t *lm_src, const uint16_t *src_cb, const uint16_t *src_cr,
                                  uint16_t *const lm_dst, uint16_t *dst_cb, uint16_t *dst_cr,
                                  ptrdiff_t lm_src_stride, ptrdiff_t src_c_stride,
                                  int abv_step, int nb_sample_abv, uint8_t lft_avail);



static void compute_lm_subsample_collocated(const uint16_t *lm_src, uint16_t *dst_cb, uint16_t *dst_cr,
                                            ptrdiff_t lm_src_stride, ptrdiff_t dst_stride_c,
                                            const struct CCLMParams *const lm_params,
                                            int pb_w, int pb_h, uint8_t lft_avail, uint8_t abv_avail);

static void vvc_intra_mdlm_l2_collocated(const uint16_t *lm_src, uint16_t *dst_cb, uint16_t *dst_cr,
                  ptrdiff_t lm_src_stride, ptrdiff_t dst_stride_c,
                  int log2_pb_w, int log2_pb_h, uint8_t lft_avail, uint8_t abv_avail,
                  uint8_t y0, uint8_t x0, uint64_t lft_map);

static void vvc_intra_mdlm_t2_collocated(const uint16_t *lm_src, uint16_t *dst_cb, uint16_t *dst_cr,
                  ptrdiff_t lm_src_stride, ptrdiff_t dst_stride_c,
                  int log2_pb_w, int log2_pb_h, uint8_t lft_avail, uint8_t abv_avail,
                  uint8_t y0, uint8_t x0, uint64_t lft_map);

void vvc_intra_cclm2_collocated(const uint16_t *lm_src, uint16_t *dst_cb, uint16_t *dst_cr,
                                ptrdiff_t lm_src_stride, ptrdiff_t dst_stride_c,
                                int log2_pb_w, int log2_pb_h, uint8_t lft_avail, uint8_t abv_avail,
                                uint8_t y0);

void
vvc_intra_cclm_cl(const uint16_t *const src_luma, uint16_t *const dst_cb,
                  uint16_t *const dst_cr, int log2_pb_w, int log2_pb_h,
                  int y0, int up_available, int left_available)
{
    vvc_intra_cclm2_collocated(src_luma, dst_cb, dst_cr, RCN_CTB_STRIDE, RCN_CTB_STRIDE,
               log2_pb_w, log2_pb_h, left_available, up_available, y0);
}

void
vvc_intra_mdlm_top_cl(const uint16_t *const src_luma,
                      uint16_t *const dst_cb, uint16_t *const dst_cr,
                      uint64_t intra_map_rows, int log2_pb_w,
                      int log2_pb_h, int x0, int y0,
                      uint8_t left_available , uint8_t up_available)
{
    vvc_intra_mdlm_t2_collocated(src_luma, dst_cb, dst_cr, RCN_CTB_STRIDE, RCN_CTB_STRIDE,
                      log2_pb_w, log2_pb_h, left_available, up_available,
                      y0, x0, intra_map_rows);
}

void
vvc_intra_mdlm_left_cl(const uint16_t *const src_luma,
                       uint16_t *const dst_cb, uint16_t *const dst_cr,
                       uint64_t intra_map_cols, int log2_pb_w,
                       int log2_pb_h, int x0, int y0,
                       uint8_t left_available, uint8_t up_available)
{
    vvc_intra_mdlm_l2_collocated(src_luma, dst_cb, dst_cr, RCN_CTB_STRIDE, RCN_CTB_STRIDE,
                      log2_pb_w, log2_pb_h, left_available, up_available,
                      y0, x0, intra_map_cols);
}

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
    int nb_sample;

    uint8_t padd_abv = start_pos == 0 && !abv_avail;
    for (nb_sample = 0; nb_sample < nb_sample_lft; nb_sample++) {
        int s = 4;

        s += _src[-(padd_abv ? 0 : lm_src_stride)];
        s += _src[0] * 4;
        s += _src[-1];
        s += _src[1];
        s += _src[lm_src_stride];

        lm_dst[nb_sample] = s >> 3;

        dst_cb[nb_sample] = _src_cb[0];
        dst_cr[nb_sample] = _src_cr[0];
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
    int nb_sample;

    int abv_step_l = abv_step << 1;

    uint8_t pad_left = start_pos == 0 && !lft_avail;

    for (nb_sample = 0; nb_sample < nb_sample_abv; nb_sample++) {
        int s = 4;

        s += _src[0 - lm_src_stride];
        s += _src[0] * 4;
        s += _src[0 - (!pad_left)];
        s += _src[0 + 1];
        s += _src[0 + lm_src_stride];

        lm_dst[nb_sample] = s >> 3;
        dst_cb[nb_sample] = _src_cb[0];
        dst_cr[nb_sample] = _src_cr[0];

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
            dst_cb[i] = ov_clip(((lm_val * scale_cb) >> shift_cb) + offset_cb, 0, 1023);
            dst_cr[i] = ov_clip(((lm_val * scale_cr) >> shift_cr) + offset_cr, 0, 1023);
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

    int log2_rng_c_plus1 = floor_log2(OVABS(range_c)) + 1;
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
        /*FIXME better way of sorting LUT ? ?*/

        int8_t idx[4] = { 0, 2, 1, 3 };

        int8_t *min_idx = &idx[0];
        int8_t *max_idx = &idx[2];

        #if 0
        uint8_t min_idxs = 0x20;
        uint8_t max_idxs = 0x31;
        #endif

        if (lm_smp[0] > lm_smp[2]) SWAP(int8_t, min_idx[0], min_idx[1]);
        if (lm_smp[1] > lm_smp[3]) SWAP(int8_t, max_idx[0], max_idx[1]);

        if (lm_smp[min_idx[0]] > lm_smp[max_idx[1]]) SWAP(int16_t*, min_idx, max_idx);
        if (lm_smp[min_idx[1]] > lm_smp[max_idx[0]]) SWAP(int8_t, min_idx[1], max_idx[0]);

        /* FIXME average should be performed in an other function */
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
    int nb_sample;

    for (nb_sample = 0; nb_sample < nb_sample_lft; nb_sample++) {
        int s = 4;
        s += _src[ 0] * 2;
        s += _src[ 1];
        s += _src[-1];
        s += _src[lm_src_stride] * 2;
        s += _src[lm_src_stride + 1];
        s += _src[lm_src_stride - 1];

        lm_dst[nb_sample] = s >> 3;

        dst_cb[nb_sample] = _src_cb[0];
        dst_cr[nb_sample] = _src_cr[0];

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
    int nb_sample;

    int abv_step_l = abv_step << 1;
    uint8_t pad_left = start_pos == 0 && !lft_avail;

    for (nb_sample = 0; nb_sample < nb_sample_abv; nb_sample++) {

        int s = 2;
        s += _src[0              ] * 2;
        s += _src[0 - (!pad_left)];
        s += _src[0 + 1          ];

        lm_dst[nb_sample] = s >> 2;
        dst_cb[nb_sample] = _src_cb[0];
        dst_cr[nb_sample] = _src_cr[0];
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
    const uint16_t *_src = lm_src - (lm_src_stride << 1) + (start_pos << 1);
    const uint16_t *_src_cb = src_cb - src_c_stride + start_pos;
    const uint16_t *_src_cr = src_cr - src_c_stride + start_pos;
    int nb_sample;

    int abv_step_l = abv_step << 1;

    uint8_t pad_left = start_pos == 0 && !lft_avail;

    for (nb_sample = 0; nb_sample < nb_sample_abv; nb_sample++) {
        int s = 4;
        s += _src[0] * 2;
        s += _src[0 + 1];
        s += _src[0 - (!pad_left)];
        s += _src[0 + lm_src_stride] * 2;
        s += _src[0 + 1 + lm_src_stride];
        s += _src[0 + lm_src_stride - (!pad_left)];

        lm_dst[nb_sample] = s >> 3;
        dst_cb[nb_sample] = _src_cb[0];
        dst_cr[nb_sample] = _src_cr[0];

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
            dst_cb[i] = ov_clip(((lm_val * scale_cb) >> shift_cb) + offset_cb, 0, 1023);
            dst_cr[i] = ov_clip(((lm_val * scale_cr) >> shift_cr) + offset_cr, 0, 1023);
        }
        dst_cb += dst_stride_c;
        dst_cr += dst_stride_c;
        lm_src += lm_src_stride2;
    }
}

void
vvc_intra_cclm2(const uint16_t *lm_src, uint16_t *dst_cb, uint16_t *dst_cr,
               ptrdiff_t lm_src_stride, ptrdiff_t dst_stride_c,
               int log2_pb_w, int log2_pb_h, uint8_t lft_avail, uint8_t abv_avail,
               uint8_t y0)
{
    struct CCLMParams lm_params = {
        .cb = {.a = 0, .b = 512, .shift = 0},
        .cr = {.a = 0, .b = 512, .shift = 0}
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
            int ref_length_abv = pb_w;
            int abv_step = OVMAX(1, (ref_length_abv >> log2_nb_smp_abv));
            uint8_t ctu_first_line = !y0;

            /* in case of abv only ref_length might be 2 while nb_sample_lft is 4
               We are forced to reduce nb_smp in this particular case*/
            nb_sample_abv = OVMIN(ref_length_abv, nb_sample_abv);

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
            int ref_length_lft = pb_h;
            int lft_step = OVMAX(1, (ref_length_lft >> log2_nb_smp_lft));

            /* in case of left only ref_length might be 2 while nb_sample_lft is 4
               We are forced to reduce nb_smp in this particular case */
            nb_sample_lft = OVMIN(ref_length_lft, nb_sample_lft);

            sub_sample_lm_ref_lft(lm_src, dst_cb, dst_cr, &lm_smp[nb_sample_abv],
                                  &smp_cb[nb_sample_abv], &smp_cr[nb_sample_abv],
                                  lm_src_stride, dst_stride_c,
                                  lft_step, nb_sample_lft);
        }

        nb_sample = nb_sample_lft + nb_sample_abv;

        avg_min_max = sort_average_lm_ref_samples(lm_smp, smp_cb, smp_cr, nb_sample);

        lm_params = derive_cclm_params(&avg_min_max);
    }

    compute_lm_subsample(lm_src, dst_cb, dst_cr, lm_src_stride, dst_stride_c,
                         &lm_params, pb_w, pb_h, lft_avail);
}

void
vvc_intra_cclm2_collocated(const uint16_t *lm_src, uint16_t *dst_cb, uint16_t *dst_cr,
                           ptrdiff_t lm_src_stride, ptrdiff_t dst_stride_c,
                           int log2_pb_w, int log2_pb_h, uint8_t lft_avail, uint8_t abv_avail,
                           uint8_t y0)
{
    struct CCLMParams lm_params = {
        .cb = {.a = 0, .b = 512, .shift = 0},
        .cr = {.a = 0, .b = 512, .shift = 0}
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
            int ref_length_abv = pb_w;
            int abv_step = OVMAX(1, (ref_length_abv >> log2_nb_smp_abv));
            uint8_t ctu_first_line = !y0;

            /* in case of abv only ref_length might be 2 while nb_sample_lft is 4
               We are forced to reduce nb_smp in this particular case*/
            nb_sample_abv = OVMIN(ref_length_abv, nb_sample_abv);

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
            int ref_length_lft = pb_h;
            int lft_step = OVMAX(1, (ref_length_lft >> log2_nb_smp_lft));

            /* in case of left only ref_length might be 2 while nb_sample_lft is 4
               We are forced to reduce nb_smp in this particular case */
            nb_sample_lft = OVMIN(ref_length_lft, nb_sample_lft);

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

void
vvc_intra_cclm(const uint16_t *const src_luma, uint16_t *const dst_cb,
               uint16_t *const dst_cr, int log2_pb_w, int log2_pb_h,
               int y0, int up_available, int left_available)
{
    vvc_intra_cclm2(src_luma, dst_cb, dst_cr, RCN_CTB_STRIDE, RCN_CTB_STRIDE,
               log2_pb_w, log2_pb_h, left_available, up_available, y0);
}

static void
vvc_intra_mdlm_t2(const uint16_t *lm_src, uint16_t *dst_cb, uint16_t *dst_cr,
                  ptrdiff_t lm_src_stride, ptrdiff_t dst_stride_c,
                  int log2_pb_w, int log2_pb_h, uint8_t lft_avail, uint8_t abv_avail,
                  uint8_t y0, uint8_t x0, uint64_t abv_map)
{
    struct CCLMParams lm_params = {
        .cb = {.a = 0, .b = 512, .shift = 0},
        .cr = {.a = 0, .b = 512, .shift = 0}
    };

    int pb_w = 1 << log2_pb_w;
    int pb_h = 1 << log2_pb_h;

    if (abv_avail) {
        struct AVGMinMax avg_min_max;
        uint16_t lm_smp[4];
        uint16_t smp_cb[4];
        uint16_t smp_cr[4];

        /* min 1 << log2 = 1 << ctz ((1 << log2) | (1 << log2)) */
        int ref_length_abv = pb_w + OVMIN((1 << log2_pb_h), (1 << log2_pb_w));
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

    compute_lm_subsample(lm_src, dst_cb, dst_cr, lm_src_stride, dst_stride_c,
                         &lm_params, pb_w, pb_h, lft_avail);
}

static void
vvc_intra_mdlm_t2_collocated(const uint16_t *lm_src, uint16_t *dst_cb, uint16_t *dst_cr,
                             ptrdiff_t lm_src_stride, ptrdiff_t dst_stride_c,
                             int log2_pb_w, int log2_pb_h, uint8_t lft_avail, uint8_t abv_avail,
                             uint8_t y0, uint8_t x0, uint64_t abv_map)
{
    struct CCLMParams lm_params = {
        .cb = {.a = 0, .b = 512, .shift = 0},
        .cr = {.a = 0, .b = 512, .shift = 0}
    };

    int pb_w = 1 << log2_pb_w;
    int pb_h = 1 << log2_pb_h;

    if (abv_avail) {
        struct AVGMinMax avg_min_max;
        uint16_t lm_smp[4];
        uint16_t smp_cb[4];
        uint16_t smp_cr[4];

        /* min 1 << log2 = 1 << ctz ((1 << log2) | (1 << log2)) */
        int ref_length_abv = pb_w + OVMIN((1 << log2_pb_h), (1 << log2_pb_w));
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

void
vvc_intra_mdlm_top(const uint16_t *const src_luma,
                   uint16_t *const dst_cb, uint16_t *const dst_cr,
                   uint64_t intra_map_rows, int log2_pb_w,
                   int log2_pb_h, int x0, int y0,
                   uint8_t left_available , uint8_t up_available)
{
    vvc_intra_mdlm_t2(src_luma, dst_cb, dst_cr, RCN_CTB_STRIDE, RCN_CTB_STRIDE,
                      log2_pb_w, log2_pb_h, left_available, up_available,
                      y0, x0, intra_map_rows);
}

static void
vvc_intra_mdlm_l2(const uint16_t *lm_src, uint16_t *dst_cb, uint16_t *dst_cr,
                  ptrdiff_t lm_src_stride, ptrdiff_t dst_stride_c,
                  int log2_pb_w, int log2_pb_h, uint8_t lft_avail, uint8_t abv_avail,
                  uint8_t y0, uint8_t x0, uint64_t lft_map)
{
    struct CCLMParams lm_params = {
        .cb = {.a = 0, .b = 512, .shift = 0},
        .cr = {.a = 0, .b = 512, .shift = 0}
    };

    int pb_w = 1 << log2_pb_w;
    int pb_h = 1 << log2_pb_h;

    if (lft_avail) {
        struct AVGMinMax avg_min_max;
        uint16_t lm_smp[4];
        uint16_t smp_cb[4];
        uint16_t smp_cr[4];

        /* min 1 << log2 = 1 << ctz ((1 << log2) | (1 << log2)) */
        int ref_length_lft = pb_h + OVMIN((1 << log2_pb_h), (1 << log2_pb_w));
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

    compute_lm_subsample(lm_src, dst_cb, dst_cr, lm_src_stride, dst_stride_c,
                         &lm_params, pb_w, pb_h, lft_avail);
}

static void
vvc_intra_mdlm_l2_collocated(const uint16_t *lm_src, uint16_t *dst_cb, uint16_t *dst_cr,
                             ptrdiff_t lm_src_stride, ptrdiff_t dst_stride_c,
                             int log2_pb_w, int log2_pb_h, uint8_t lft_avail, uint8_t abv_avail,
                             uint8_t y0, uint8_t x0, uint64_t lft_map)
{
    struct CCLMParams lm_params = {
        .cb = {.a = 0, .b = 512, .shift = 0},
        .cr = {.a = 0, .b = 512, .shift = 0}
    };

    int pb_w = 1 << log2_pb_w;
    int pb_h = 1 << log2_pb_h;

    if (lft_avail) {
        struct AVGMinMax avg_min_max;
        uint16_t lm_smp[4];
        uint16_t smp_cb[4];
        uint16_t smp_cr[4];

        /* min 1 << log2 = 1 << ctz ((1 << log2) | (1 << log2)) */
        int ref_length_lft = pb_h + OVMIN((1 << log2_pb_h), (1 << log2_pb_w));
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

void
vvc_intra_mdlm_left(const uint16_t *const src_luma,
                    uint16_t *const dst_cb, uint16_t *const dst_cr,
                    uint64_t intra_map_cols, int log2_pb_w,
                    int log2_pb_h, int x0, int y0,
                    uint8_t left_available, uint8_t up_available)
{
    vvc_intra_mdlm_l2(src_luma, dst_cb, dst_cr, RCN_CTB_STRIDE, RCN_CTB_STRIDE,
                      log2_pb_w, log2_pb_h, left_available, up_available,
                      y0, x0, intra_map_cols);
}

void
rcn_init_cclm_functions_collocated(struct RCNFunctions *rcn_func)
{
   struct CCLMFunctions *const cclm = &rcn_func->cclm; 
   cclm->cclm = &vvc_intra_cclm_cl;
   cclm->mdlm_left = &vvc_intra_mdlm_left_cl;
   cclm->mdlm_top  = &vvc_intra_mdlm_top_cl;
}

/* FIXME check vertical / horizontal */
void
rcn_init_cclm_functions(struct RCNFunctions *rcn_func)
{
   struct CCLMFunctions *const cclm = &rcn_func->cclm; 
   cclm->cclm = &vvc_intra_cclm;
   cclm->mdlm_left = &vvc_intra_mdlm_left;
   cclm->mdlm_top  = &vvc_intra_mdlm_top;
}

#undef W_SHIFT
#undef H_SHIFT
