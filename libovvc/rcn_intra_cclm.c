#define W_SHIFT 1
#define H_SHIFT 1
void
vvc_intra_cclm_cl(const uint16_t* const src_luma, uint16_t* const dst_cb,
                  uint16_t* const dst_cr, int log2_pb_w, int log2_pb_h, int y0,
                  int up_available, int left_available)
{ // left and above
        int i, j, i_s;

        uint16_t lm_ref[4];
        uint16_t cb_ref[4];
        uint16_t cr_ref[4];
        VVCLMParams lm_params;
        const int src_stride = VVC_CTB_STRIDE;
        const int dst_stride = VVC_CTB_STRIDE_CHROMA;

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

void
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
        const int src_stride = VVC_CTB_STRIDE;
        const int dst_stride = VVC_CTB_STRIDE_CHROMA;

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

void
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
        const int src_stride = VVC_CTB_STRIDE;
        const int dst_stride = VVC_CTB_STRIDE_CHROMA;

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
void
vvc_intra_cclm(const uint16_t* const src_luma, uint16_t* const dst_cb,
               uint16_t* const dst_cr, int log2_pb_w, int log2_pb_h, int y0,
               int up_available, int left_available)
{ // left and above
        int i, j, i_s;

        uint16_t lm_ref[4];
        uint16_t cb_ref[4];
        uint16_t cr_ref[4];
        VVCLMParams lm_params;
        const int src_stride = VVC_CTB_STRIDE;
        const int dst_stride = VVC_CTB_STRIDE_CHROMA;

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

void
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
        const int src_stride = VVC_CTB_STRIDE;
        const int dst_stride = VVC_CTB_STRIDE_CHROMA;

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

void
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
        const int src_stride = VVC_CTB_STRIDE;
        const int dst_stride = VVC_CTB_STRIDE_CHROMA;

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
#undef W_SHIFT
#undef H_SHIFT
