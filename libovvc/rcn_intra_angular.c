#include <stddef.h>
#include <stdint.h>
#include <string.h>

#include "ovutils.h"

#include "rcn_intra_angular.h"
#include "data_rcn_angular.h"



static const int8_t chroma_filter[4 * 32] = {
     0, 64,  0,  0,
    -1, 63,  2,  0,
    -2, 62,  4,  0,
    -2, 60,  7, -1,
    -2, 58, 10, -2,
    -3, 57, 12, -2,
    -4, 56, 14, -2,
    -4, 55, 15, -2,
    -4, 54, 16, -2,
    -5, 53, 18, -2,
    -6, 52, 20, -2,
    -6, 49, 24, -3,
    -6, 46, 28, -4,
    -5, 44, 29, -4,
    -4, 42, 30, -4,
    -4, 39, 33, -4,
    -4, 36, 36, -4,
    -4, 33, 39, -4,
    -4, 30, 42, -4,
    -4, 29, 44, -5,
    -4, 28, 46, -6,
    -3, 24, 49, -6,
    -2, 20, 52, -6,
    -2, 18, 53, -5,
    -2, 16, 54, -4,
    -2, 15, 55, -4,
    -2, 14, 56, -4,
    -2, 12, 57, -3,
    -2, 10, 58, -2,
    -1,  7, 60, -2,
     0,  4, 62, -2,
     0,  2, 63, -1
};



static const uint8_t vvc_pdpc_w[3][128] = {
        { 32, 8, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
        { 32, 16, 8, 4, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
        { 32, 32, 16, 16, 8, 8, 4, 4, 2, 2, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

void
vvc_intra_angular_hdia(const uint16_t* const ref_abv,
                       const uint16_t* const ref_lft, uint16_t* const dst,
                       ptrdiff_t dst_stride, int log2_pb_w,
                       int log2_pb_h)
{

        int16_t tmp_dst[128 * 128];
        const int tmp_stride = 128;
        int16_t* _tmp = tmp_dst;
        uint16_t* _dst = dst;
        int width = 1 << log2_pb_w;
        int height = 1 << log2_pb_h;
        int scale = OVMIN(2, log2_pb_w - (floor_log2(3 * 512 - 2) - 8));

        int delta_pos = 32; // angle val is 32 for diag
        for (int y = 0; y < width; y++) {
                const int delta_int = delta_pos >> 5;
                #if 0
                int wT = 16 >> OVMIN(31, ((y << 1) >> scale));
                #endif
                for (int x = 0; x < height; x++) {
                        _tmp[x] = ref_lft[x + delta_int + 1];
                }
                if (height >= 4 && scale >= 0)
                        for (int x = 0; x < OVMIN(3 << scale, height); x++) {
                                int wL = 32 >> (2 * x >> scale);
                                const int16_t above = ref_abv[y + x + 2];
                                _tmp[x] = ov_clip(
                                  _tmp[x] +
                                    ((wL * (above - _tmp[x]) + 32) >> 6),
                                  0,
                                  1023);
                        }
                delta_pos += 32;
                _tmp += tmp_stride;
        }

        _tmp = tmp_dst;
        // Transpose block
        for (int y = 0; y < width; y++) {
                _dst = &dst[y];
                for (int x = 0; x < height; x++) {
                        _dst[0] = _tmp[x];
                        _dst += dst_stride;
                }
                _tmp += tmp_stride;
        }
}

void
vvc_intra_angular_vdia(const uint16_t* const ref_abv,
                       const uint16_t* const ref_lft, uint16_t* const dst,
                       ptrdiff_t dst_stride, int log2_pb_w,
                       int log2_pb_h)
{
        int delta_pos = 32; // angle val = 32 for strict diag
        uint16_t* _dst = dst;
        int width = 1 << log2_pb_w;
        int height = 1 << log2_pb_h;
        int scale = OVMIN(2, log2_pb_h - (floor_log2(3 * 512 - 2) - 8));

        for (int y = 0; y < height; y++) {
                const int delta_int = delta_pos >> 5;
                #if 0
                int wT = 16 >> OVMIN(31, ((y << 1) >> scale));
                #endif

                for (int x = 0; x < width; x++) {
                        _dst[x] = ref_abv[x + delta_int + 1];
                }
                if (height >= 4 && scale >= 0)
                        for (int x = 0; x < OVMIN(3 << scale, width); x++) {
                                int wL = 32 >> (2 * x >> scale);
                                const int16_t left = ref_lft[y + x + 2];
                                _dst[x] = ov_clip(
                                  _dst[x] + ((wL * (left - _dst[x]) + 32) >> 6),
                                  0,
                                  1023);
                        }
                delta_pos += 32;
                _dst += dst_stride;
        }
}

void
vvc_intra_angular_h_c(const uint16_t* ref_lft, uint16_t* dst,
                      ptrdiff_t dst_stride, int log2_pb_w,
                      int log2_pb_h, int angle_val)
{
        int16_t tmp_dst[128 * 128];
        int16_t* _tmp = tmp_dst;
        uint16_t* _dst = dst;
        int width = 1 << log2_pb_w;
        int height = 1 << log2_pb_h;
        const int tmp_stride = 128;

        int delta_pos = angle_val;

        for (int y = 0; y < width; y++) {
                const int delta_int = delta_pos >> 5;
                const int delta_frac = delta_pos & 0x1F;
                // TODO for bit depth <= 10 we can use uint16_t for computation
                const uint16_t* pRM = ref_lft + delta_int + 1;
                int last_ref_val = *pRM++;
                for (int x = 0; x < height; x++) {
                        int curr_ref_val = *pRM;
                        int val;
                        val = (int16_t)last_ref_val +
                              (((delta_frac) * (curr_ref_val - last_ref_val) +
                                16) >>
                               5);
                        _tmp[x] = ov_clip(val, 0, 1023);
                        last_ref_val = curr_ref_val;
                        pRM++;
                }
                delta_pos += angle_val;
                _tmp += tmp_stride;
        }

        _tmp = tmp_dst;
        // Transpose block
        for (int y = 0; y < width; y++) {
                _dst = &dst[y];
                for (int x = 0; x < height; x++) {
                        _dst[0] = _tmp[x];
                        _dst += dst_stride;
                }
                _tmp += tmp_stride;
        }
}

void
vvc_intra_angular_v_c(const uint16_t* ref_abv, uint16_t* dst,
                      ptrdiff_t dst_stride, int log2_pb_w,
                      int log2_pb_h, int angle_val)
{
        uint16_t* _dst = dst;
        int width = 1 << log2_pb_w;
        int height = 1 << log2_pb_h;

        for (int y = 0, delta_pos = angle_val; y < height; y++) {
                const int delta_int = delta_pos >> 5;
                const int delta_frac = delta_pos & 0x1F;
                const uint16_t* pRM = ref_abv + delta_int + 1;
                int last_ref_val = *pRM++;
                for (int x = 0; x < width; pRM++, x++) {
                        int curr_ref_val = *pRM;
                        int val;
                        val = (int16_t)last_ref_val +
                              (((delta_frac) * (curr_ref_val - last_ref_val) +
                                16) >>
                               5);
                        _dst[x] = ov_clip(val, 0, 1023);
                        last_ref_val = curr_ref_val;
                }
                delta_pos += angle_val;
                _dst += dst_stride;
        }
}

void
vvc_intra_hor_pdpc(const uint16_t* const ref_abv,
                   const uint16_t* const ref_lft, uint16_t* const dst,
                   ptrdiff_t dst_stride, uint16_t log2_pb_w,
                   uint16_t log2_pb_h)
{
        uint16_t* _dst = dst;
        int pb_width = 1 << log2_pb_w;
        int pb_height = 1 << log2_pb_h;
        int pdpc_scale = (log2_pb_w + log2_pb_h - 2) >> 2;
        const uint8_t* pdpc_w = vvc_pdpc_w[pdpc_scale];

        const uint16_t tl_val = ref_abv[0];
        for (int y = 0; y < pb_height; y++) {
                int l_wgh = pdpc_w[y];
                int32_t l_val = ref_lft[y + 1];
                for (int x = 0; x < pb_width; x++) {
                        const int32_t t_val = ref_abv[x + 1];
                        int val =
                          (l_wgh * (t_val - tl_val) + (l_val << 6) + 32) >> 6;
                        _dst[x] = ov_clip(val, 0, 1023);
                }
                _dst += dst_stride;
        }
}

void
vvc_intra_ver_pdpc(const uint16_t* const ref_abv,
                   const uint16_t* const ref_lft, uint16_t* const dst,
                   ptrdiff_t dst_stride, uint16_t log2_pb_w,
                   uint16_t log2_pb_h)
{
        uint16_t* _dst = dst;
        const uint16_t tl_val = ref_abv[0];
        int pdpc_scale = (log2_pb_w + log2_pb_h - 2) >> 2;
        const uint8_t* pdpc_w = vvc_pdpc_w[pdpc_scale];

        int pb_width = 1 << log2_pb_w;
        int pb_height = 1 << log2_pb_h;

        for (int y = 0; y < pb_height; y++) {
                const uint16_t l_val = ref_lft[y + 1];
                for (int x = 0; x < pb_width; x++) {
                        const int32_t t_val = ref_abv[x + 1];
                        int l_wgh = pdpc_w[x];
                        int val =
                          (l_wgh * (l_val - tl_val) + (t_val << 6) + 32) >> 6;
                        _dst[x] = ov_clip(val, 0, 1023);
                }
                _dst += dst_stride;
        }
}

void
vvc_intra_hor(const uint16_t* const ref_abv, const uint16_t* const ref_lft,
              uint16_t* const dst, ptrdiff_t dst_stride, uint16_t log2_pb_w,
              uint16_t log2_pb_h)
{
        uint16_t* _dst = dst;
        int pb_width = 1 << log2_pb_w;
        int pb_height = 1 << log2_pb_h;

        for (int y = 0; y < pb_height; y++) {
                for (int j = 0; j < pb_width; j++) {
                        _dst[j] = ref_lft[y + 1];
                }
                //        memset(_dst, ref_lft[y+1], sizeof(uint16_t) *
                //        cb_width);
                _dst += dst_stride;
        }
}

void
vvc_intra_ver(const uint16_t* const ref_abv, const uint16_t* const ref_lft,
              uint16_t* const dst, ptrdiff_t dst_stride, uint16_t log2_pb_w,
              uint16_t log2_pb_h)
{
        uint16_t* _dst = dst;
        int pb_width = 1 << log2_pb_w;
        int pb_height = 1 << log2_pb_h;

        for (int y = 0; y < pb_height; y++) {
                memcpy(_dst, &ref_abv[1], sizeof(uint16_t) * pb_width);
                _dst += dst_stride;
        }
}

void
vvc_intra_angular_hpos_wide(const uint16_t* const ref_abv,
                            const uint16_t* const ref_lft, uint16_t* const dst,
                            ptrdiff_t dst_stride, int log2_pb_w,
                            int log2_pb_h, int mode_idx)
{
        int16_t tmp_dst[128 * 128];
        const int tmp_stride = 128;
        int16_t* _tmp = tmp_dst;
        uint16_t* _dst = dst;
        int width = 1 << log2_pb_w;
        int height = 1 << log2_pb_h;
        int angle_val = angle_table[mode_idx];
        int inv_angle = inverse_angle_table[mode_idx];
        int delta_pos = angle_val;
        int scale = OVMIN(2,
                          log2_pb_w -
                            (floor_log2(3 * inv_angle - 2) -
                             8)); //(log2_pb_w + log2_pb_h - 2) >> 2;
        #if 0
        int top_ref_length = 1 << (log2_pb_w + 1);
        #endif

        // swapped width/height for horizontal mode:
        for (int y = 0; y < width; y++) {
                const int delta_int = delta_pos >> 5;        // Integer part
                const int delta_frac = delta_pos & 0x1F; // Fractionnal part
                int inv_angle_sum = 256 + inv_angle;

                // Do linear filtering
                const uint16_t* pRM = ref_lft + delta_int + 1;
                int last_ref_val = *pRM++;
                for (int x = 0; x < height; pRM++, x++) {
                        int curr_ref_val = *pRM;
                        _tmp[x] =
                          (int16_t)last_ref_val +
                          (((delta_frac) * (curr_ref_val - last_ref_val) +
                            16) >>
                           5);
                        last_ref_val = curr_ref_val;
                }
                // TODO move this part to previous loop
                for (int x = 0; x < OVMIN(3 << scale, height); x++) {
                        // TODO check if we can use LUTs instead
                        int wL = 32 >> ((x << 1) >> scale);
                        const uint16_t* p =
                          ref_abv + y + (inv_angle_sum >> 9) + 1;

                        int32_t left = p[0];
                        _tmp[x] =
                          ov_clip(_tmp[x] + ((wL * (left - _tmp[x]) + 32) >> 6),
                                  0,
                                  1023);
                        inv_angle_sum += inv_angle;
                }
                delta_pos += angle_val;
                _tmp += tmp_stride;
        }

        _tmp = tmp_dst;
        // Transpose block
        for (int y = 0; y < width; y++) {
                _dst = &dst[y];
                for (int x = 0; x < height; x++) {
                        _dst[0] = _tmp[x];
                        _dst += dst_stride;
                }
                _tmp += tmp_stride;
        }
}

void
vvc_intra_angular_vpos_wide(const uint16_t* const ref_abv,
                            const uint16_t* const ref_lft, uint16_t* const dst,
                            ptrdiff_t dst_stride, int log2_pb_w,
                            int log2_pb_h, int mode_idx)
{

        uint16_t* _dst = dst;
        int angle_val = angle_table[mode_idx];
        int inv_angle = inverse_angle_table[mode_idx];
        int delta_pos = angle_val;
        int width = 1 << log2_pb_w;
        int height = 1 << log2_pb_h;
        int scale = OVMIN(2,
                          log2_pb_h -
                            (floor_log2(3 * inv_angle - 2) -
                             8)); //(log2_pb_w + log2_pb_h - 2) >> 2;
        for (int y = 0; y < height; y++) {
                const int delta_int = delta_pos >> 5;
                const int delta_frac = delta_pos & 0x1F;
                int inv_angle_sum = 256 + inv_angle;

                // Do linear filtering
                const uint16_t* pRM = ref_abv + delta_int + 1;
                int last_ref_val = *pRM++;
                for (int x = 0; x < width; pRM++, x++) {
                        int curr_ref_val = *pRM;
                        _dst[x] =
                          (int16_t)last_ref_val +
                          (((delta_frac) * (curr_ref_val - last_ref_val) +
                            16) >>
                           5);
                        last_ref_val = curr_ref_val;
                }
                for (int x = 0; x < OVMIN(3 << scale, width); x++) {
                        // TODO check if we can use LUTs instead
                        int wL = 32 >> ((x << 1) >> scale);
                        const uint16_t* p =
                          ref_lft + y + (inv_angle_sum >> 9) + 1;

                        int32_t left = p[0];
                        _dst[x] =
                          ov_clip(_dst[x] + ((wL * (left - _dst[x]) + 32) >> 6),
                                  0,
                                  1023);
                        inv_angle_sum += inv_angle;
                }
                delta_pos += angle_val;
                _dst += dst_stride;
        }
}

void
intra_angular_h_nofrac(const uint16_t* ref_lft, uint16_t* dst,
                       ptrdiff_t dst_stride, int log2_pb_w,
                       int log2_pb_h, int angle_val)
{
        uint16_t tmp_dst[128 * 128];
        const int tmp_stride = 128;
        uint16_t* _tmp = tmp_dst;
        uint16_t* _dst = dst;
        int width = 1 << log2_pb_w;
        int height = 1 << log2_pb_h;
        int delta_pos = angle_val >> 5;
        int y, x;

        for (y = 0; y < width; ++y) {
                for (x = 0; x < height; ++x) {
                        _tmp[x] = ref_lft[x + delta_pos + 1];
                }
                delta_pos += angle_val >> 5;
                _tmp += tmp_stride;
        }

        _tmp = tmp_dst;

        // Transpose block
        for (int y = 0; y < width; y++) {
                _dst = &dst[y];
                for (int x = 0; x < height; x++) {
                        _dst[0] = _tmp[x];
                        _dst += dst_stride;
                }
                _tmp += tmp_stride;
        }
}

void
intra_angular_h_nofrac_pdpc(const uint16_t* ref_abv, const uint16_t* ref_lft,
                            uint16_t* dst, ptrdiff_t dst_stride,
                            int log2_pb_w, int log2_pb_h, int mode_idx)
{
        int16_t tmp_dst[128 * 128];
        const int tmp_stride = 128;
        int16_t* _tmp = tmp_dst;
        uint16_t* _dst = dst;
        int angle_val = angle_table[mode_idx];
        int inv_angle = inverse_angle_table[mode_idx];
        int width = 1 << log2_pb_w;
        int height = 1 << log2_pb_h;
        int delta_pos = angle_val >> 5;
        int scale =
          OVMIN(2, log2_pb_w - (floor_log2(3 * inv_angle - 2) - 8));
        int y/*, x*/;

        for (y = 0; y < width; ++y) {
                int inv_angle_sum = 256 + inv_angle;
                for (int x = 0; x < height; x++) {
                        _tmp[x] = ref_lft[x + delta_pos + 1];
                }
                for (int x = 0; x < OVMIN(3 << scale, height); x++) {
                        int wL = 32 >> ((x << 1) >> scale);
                        const uint16_t* p =
                          ref_abv + y + (inv_angle_sum >> 9) + 1;

                        int16_t left = p[0];
                        _tmp[x] =
                          ov_clip(_tmp[x] + ((wL * (left - _tmp[x]) + 32) >> 6),
                                  0,
                                  1023);
                        inv_angle_sum += inv_angle;
                }
                delta_pos += angle_val >> 5;
                _tmp += tmp_stride;
        }

        _tmp = tmp_dst;

        // Transpose block
        for (int y = 0; y < width; y++) {
                _dst = &dst[y];
                for (int x = 0; x < height; x++) {
                        _dst[0] = _tmp[x];
                        _dst += dst_stride;
                }
                _tmp += tmp_stride;
        }
}

void
intra_angular_h_gauss_pdpc(const uint16_t* ref_abv, const uint16_t* ref_lft,
                           uint16_t* const dst, ptrdiff_t dst_stride,
                           int log2_pb_w, int log2_pb_h, int mode_idx)
{
        uint16_t tmp_dst[128 * 128];
        const int tmp_stride = 128;
        uint16_t* _tmp = tmp_dst;
        uint16_t* _dst = dst;
        int angle_val = angle_table[mode_idx];
        int inv_angle = inverse_angle_table[mode_idx];
        int width = 1 << log2_pb_w;
        int height = 1 << log2_pb_h;
        int delta_pos = angle_val;
        int scale =
          OVMIN(2, log2_pb_w - (floor_log2(3 * inv_angle - 2) - 8));

        for (int y = 0; y < width; y++) {
                const int delta_int = delta_pos >> 5;
                const int delta_frac = delta_pos & 0x1F;
                int inv_angle_sum = 256 + inv_angle;
                const int16_t* ref = (int16_t *)ref_lft + delta_int;
                for (int x = 0; x < height; x++) {
                        _tmp[x] =
                          ((int32_t)(ref[0] * (16 - (delta_frac >> 1))) +
                           (int32_t)(ref[1] * (32 - (delta_frac >> 1))) +
                           (int32_t)(ref[2] * (16 + (delta_frac >> 1))) +
                           (int32_t)(ref[3] * ((delta_frac >> 1))) + 32) >>
                          6;
                        ref++;
                }
                for (int x = 0; x < OVMIN(3 << scale, height); x++) {
                        int wL = 32 >> ((x << 1) >> scale);
                        const uint16_t* p =
                          ref_abv + y + (inv_angle_sum >> 9) + 1;

                        int16_t left = p[0];
                        _tmp[x] =
                          ov_clip(_tmp[x] + ((wL * (left - _tmp[x]) + 32) >> 6),
                                  0,
                                  1023);
                        inv_angle_sum += inv_angle;
                }

                delta_pos += angle_val;
                _tmp += tmp_stride;
        }

        _tmp = tmp_dst;

        for (int y = 0; y < width; y++) {
                _dst = &dst[y];
                for (int x = 0; x < height; x++) {
                        _dst[0] = _tmp[x];
                        _dst += dst_stride;
                }
                _tmp += tmp_stride;
        }
}

void
intra_angular_v_nofrac(const uint16_t* ref_abv, uint16_t* dst,
                       ptrdiff_t dst_stride, int log2_pb_w,
                       int log2_pb_h, int angle_val)
{
        uint16_t* _dst = dst;
        int width = 1 << log2_pb_w;
        int height = 1 << log2_pb_h;
        int delta_pos = angle_val;

        int y, x;
        for (y = 0; y < height; y++) {
                const int delta_int = delta_pos >> 5;
                for (x = 0; x < width; x++) {
                        _dst[x] = ref_abv[x + delta_int + 1];
                }
                delta_pos += angle_val;
                _dst += dst_stride;
        }
}

void
intra_angular_v_nofrac_pdpc(const uint16_t* ref_abv, const uint16_t* ref_lft,
                            uint16_t* const dst, ptrdiff_t dst_stride,
                            int log2_pb_w, int log2_pb_h, int mode_idx)
{
        uint16_t* _dst = dst;
        int angle_val = angle_table[mode_idx];
        int inv_angle = inverse_angle_table[mode_idx];
        int width = 1 << log2_pb_w;
        int height = 1 << log2_pb_h;
        int delta_pos = angle_val;
        int scale =
          OVMIN(2, log2_pb_h - (floor_log2(3 * inv_angle - 2) - 8));

        int y, x;
        for (y = 0; y < height; ++y) {
                const int delta_int = delta_pos >> 5;
                int inv_angle_sum = 256 + inv_angle;
                for (x = 0; x < width; x++) {
                        _dst[x] = ref_abv[x + delta_int + 1];
                }
                for (x = 0; x < OVMIN(3 << scale, width); ++x) {
                        int wL = 32 >> ((x << 1) >> scale);

                        const uint16_t* p =
                          ref_lft + y + (inv_angle_sum >> 9) + 1;

                        int16_t left = p[0];
                        _dst[x] =
                          ov_clip(_dst[x] + ((wL * (left - _dst[x]) + 32) >> 6),
                                  0,
                                  1023);
                        inv_angle_sum += inv_angle;
                }
                delta_pos += angle_val;
                _dst += dst_stride;
        }
}

void
intra_angular_v_gauss_pdpc(const uint16_t* ref_abv, const uint16_t* ref_lft,
                           uint16_t* const dst, ptrdiff_t dst_stride,
                           int log2_pb_w, int log2_pb_h, int mode_idx)
{
        uint16_t* _dst = dst;
        int angle_val = angle_table[mode_idx];
        int inv_angle = inverse_angle_table[mode_idx];
        int width = 1 << log2_pb_w;
        int height = 1 << log2_pb_h;
        int delta_pos = angle_val;
        int scale =
          OVMIN(2, log2_pb_h - (floor_log2(3 * inv_angle - 2) - 8));

        for (int y = 0; y < height; y++) {
                const int delta_int = delta_pos >> 5;
                const int delta_frac = delta_pos & 0x1F;
                int inv_angle_sum = 256 + inv_angle;
                const int16_t* ref = (int16_t*)ref_abv + delta_int;
                for (int x = 0; x < width; x++) {
                        _dst[x] =
                          ((int32_t)(ref[0] * (16 - (delta_frac >> 1))) +
                           (int32_t)(ref[1] * (32 - (delta_frac >> 1))) +
                           (int32_t)(ref[2] * (16 + (delta_frac >> 1))) +
                           (int32_t)(ref[3] * ((delta_frac >> 1))) + 32) >>
                          6;
                        ref++;
                }
                for (int x = 0; x < OVMIN(3 << scale, width); x++) {
                        int wL = 32 >> ((x << 1) >> scale);
                        const uint16_t* p =
                          ref_lft + y + (inv_angle_sum >> 9) + 1;

                        int16_t left = p[0];
                        _dst[x] =
                          ov_clip(_dst[x] + ((wL * (left - _dst[x]) + 32) >> 6),
                                  0,
                                  1023);
                        inv_angle_sum += inv_angle;
                }
                delta_pos += angle_val;
                _dst += dst_stride;
        }
}

void
intra_angular_h_cubic(const uint16_t* ref_lft, uint16_t* dst,
                      ptrdiff_t dst_stride, int log2_pb_w,
                      int log2_pb_h, int angle_val)
{
        uint16_t tmp_dst[128 * 128];
        const int tmp_stride = 128;
        int width = 1 << log2_pb_w;
        int height = 1 << log2_pb_h;
        uint16_t* _tmp = tmp_dst;
        uint16_t* _dst = dst;

        int delta_pos = angle_val;
        for (int y = 0; y < width; y++) {
                const int delta_int = delta_pos >> 5;
                const int delta_frac = delta_pos & 0x1F;
                const int16_t* ref = (int16_t*)ref_lft + delta_int;
                const int8_t* filter = &chroma_filter[delta_frac << 2];
                for (int x = 0; x < height; x++) {
                        int32_t val;
                        val = ((int32_t)(ref[0] * filter[0]) +
                               (int32_t)(ref[1] * filter[1]) +
                               (int32_t)(ref[2] * filter[2]) +
                               (int32_t)(ref[3] * filter[3]) + 32) >>
                              6;
                        _tmp[x] = ov_clip(val, 0, 1023);
                        ref++;
                }
                delta_pos += angle_val;
                _tmp += tmp_stride;
        }

        _tmp = tmp_dst;
        // Transpose block
        for (int y = 0; y < width; y++) {
                _dst = &dst[y];
                for (int x = 0; x < height; x++) {
                        _dst[0] = _tmp[x];
                        _dst += dst_stride;
                }
                _tmp += tmp_stride;
        }
}

void
intra_angular_h_gauss(const uint16_t* ref_lft, uint16_t* dst,
                      ptrdiff_t dst_stride, int log2_pb_w,
                      int log2_pb_h, int angle_val)
{
        uint16_t tmp_dst[128 * 128];
        const int tmp_stride = 128;
        int width = 1 << log2_pb_w;
        int height = 1 << log2_pb_h;
        uint16_t* _tmp = tmp_dst;
        uint16_t* _dst = dst;

        int delta_pos = angle_val;
        for (int y = 0; y < width; y++) {
                const int delta_int = delta_pos >> 5;
                const int delta_frac = delta_pos & 0x1F;
                const int16_t* ref = (int16_t*)ref_lft + delta_int;
                for (int x = 0; x < height; x++) {
                        _tmp[x] =
                          ((int32_t)(ref[0] * (16 - (delta_frac >> 1))) +
                           (int32_t)(ref[1] * (32 - (delta_frac >> 1))) +
                           (int32_t)(ref[2] * (16 + (delta_frac >> 1))) +
                           (int32_t)(ref[3] * ((delta_frac >> 1))) + 32) >>
                          6;
                        ref++;
                }
                delta_pos += angle_val;
                _tmp += tmp_stride;
        }

        _tmp = tmp_dst;
        // Transpose block
        for (int y = 0; y < width; y++) {
                _dst = &dst[y];
                for (int x = 0; x < height; x++) {
                        _dst[0] = _tmp[x];
                        _dst += dst_stride;
                }
                _tmp += tmp_stride;
        }
}

void
intra_angular_v_cubic(const uint16_t* ref_abv, uint16_t* dst,
                      ptrdiff_t dst_stride, int log2_pb_w,
                      int log2_pb_h, int angle_val)
{
        int delta_pos = angle_val;
        int width = 1 << log2_pb_w;
        int height = 1 << log2_pb_h;
        uint16_t* _dst = dst;

        for (int y = 0; y < height; y++) {
                const int delta_int = delta_pos >> 5;
                const int delta_frac = delta_pos & 0x1F;

                const int16_t* ref = (int16_t*)ref_abv + delta_int;
                const int8_t* filter = &chroma_filter[delta_frac << 2];

                for (int x = 0; x < width; x++) {
                        int32_t val;
                        val = ((int32_t)(ref[0] * filter[0]) +
                               (int32_t)(ref[1] * filter[1]) +
                               (int32_t)(ref[2] * filter[2]) +
                               (int32_t)(ref[3] * filter[3]) + 32) >>
                              6;
                        _dst[x] = ov_clip(val, 0, 1023);
                        ref++;
                }
                delta_pos += angle_val;
                _dst += dst_stride;
        }
}

void
intra_angular_v_gauss(const uint16_t* ref_abv, uint16_t* dst,
                      ptrdiff_t dst_stride, int log2_pb_w,
                      int log2_pb_h, int angle_val)
{
        int delta_pos = angle_val;
        int width = 1 << log2_pb_w;
        int height = 1 << log2_pb_h;
        uint16_t* _dst = dst;

        for (int y = 0; y < height; y++) {
                const int delta_int = delta_pos >> 5;
                const int delta_frac = delta_pos & 0x1F;

                const int16_t* ref = (int16_t*)ref_abv + delta_int;

                for (int x = 0; x < width; x++) {
                        _dst[x] =
                          ((int32_t)(ref[0] * (16 - (delta_frac >> 1))) +
                           (int32_t)(ref[1] * (32 - (delta_frac >> 1))) +
                           (int32_t)(ref[2] * (16 + (delta_frac >> 1))) +
                           (int32_t)(ref[3] * ((delta_frac >> 1))) + 32) >>
                          6;
                        ref++;
                }
                delta_pos += angle_val;
                _dst += dst_stride;
        }
}

void
intra_angular_h_cubic_pdpc(const uint16_t* ref_abv, const uint16_t* ref_lft,
                           uint16_t* const dst, ptrdiff_t dst_stride,
                           int log2_pb_w, int log2_pb_h, int mode_idx)
{
        uint16_t tmp_dst[128 * 128];
        const int tmp_stride = 128;
        uint16_t* _tmp = tmp_dst;
        uint16_t* _dst = dst;
        int angle_val = angle_table[mode_idx];
        int inv_angle = inverse_angle_table[mode_idx];
        int width = 1 << log2_pb_w;
        int height = 1 << log2_pb_h;
        int delta_pos = angle_val;
        int scale =
          OVMIN(2, log2_pb_w - (floor_log2(3 * inv_angle - 2) - 8));
        #if 0
        int top_ref_length = 1 << (log2_pb_w + 1);
        #endif

        for (int y = 0; y < width; y++) {
                const int delta_int = delta_pos >> 5;
                const int delta_frac = delta_pos & 0x1F;
                int inv_angle_sum = 256 + inv_angle;
                const int16_t* ref = (int16_t*)ref_lft + delta_int;
                const int8_t* filter = &chroma_filter[delta_frac << 2];

                for (int x = 0; x < height; x++) {
                        int val = ((int32_t)(ref[0] * filter[0]) +
                                   (int32_t)(ref[1] * filter[1]) +
                                   (int32_t)(ref[2] * filter[2]) +
                                   (int32_t)(ref[3] * filter[3]) + 32) >>
                                  6;
                        ref++;
                        _tmp[x] = ov_clip(val, 0, 1023);
                }

                for (int x = 0; x < OVMIN(3 << scale, height); x++) {
                        int wL = 32 >> ((x << 1) >> scale);
                        const uint16_t* p =
                          ref_abv + y + (inv_angle_sum >> 9) + 1;

                        int16_t left = p[0];
                        _tmp[x] =
                          ov_clip(_tmp[x] + ((wL * (left - _tmp[x]) + 32) >> 6),
                                  0,
                                  1023);
                        inv_angle_sum += inv_angle;
                }
                delta_pos += angle_val;
                _tmp += tmp_stride;
        }

        _tmp = tmp_dst;
        // Transpose block
        for (int y = 0; y < width; y++) {
                _dst = &dst[y];
                for (int x = 0; x < height; x++) {
                        _dst[0] = _tmp[x];
                        _dst += dst_stride;
                }
                _tmp += tmp_stride;
        }
}

void
intra_angular_v_cubic_pdpc(const uint16_t* ref_abv, const uint16_t* ref_lft,
                           uint16_t* const dst, ptrdiff_t dst_stride,
                           int log2_pb_w, int log2_pb_h, int mode_idx)
{
        uint16_t* _dst = dst;
        int width = 1 << log2_pb_w;
        int height = 1 << log2_pb_h;
        int angle_val = angle_table[mode_idx];
        int inv_angle = inverse_angle_table[mode_idx];
        int delta_pos = angle_val;
        int scale =
          OVMIN(2, log2_pb_h - (floor_log2(3 * inv_angle - 2) - 8));
        for (int y = 0; y < height; y++) {
                const int delta_int = delta_pos >> 5;
                const int delta_frac = delta_pos & 0x1F;
                int inv_angle_sum = 256 + inv_angle;
                const int16_t* ref = (int16_t *)ref_abv + delta_int;
                const int8_t* filter = &chroma_filter[delta_frac << 2];
                for (int x = 0; x < width; x++) {
                        int val = ((int32_t)(ref[0] * filter[0]) +
                                   (int32_t)(ref[1] * filter[1]) +
                                   (int32_t)(ref[2] * filter[2]) +
                                   (int32_t)(ref[3] * filter[3]) + 32) >>
                                  6;
                        ref++;
                        _dst[x] = ov_clip(val, 0, 1023);
                }

                for (int x = 0; x < OVMIN(3 << scale, width); x++) {
                        int wL = 32 >> ((x << 1) >> scale);
                        const uint16_t* p =
                          ref_lft + y + (inv_angle_sum >> 9) + 1;

                        int16_t left = p[0];
                        _dst[x] =
                          ov_clip(_dst[x] + ((wL * (left - _dst[x]) + 32) >> 6),
                                  0,
                                  1023);
                        inv_angle_sum += inv_angle;
                }
                delta_pos += angle_val;
                _dst += dst_stride;
        }
}

void
intra_angular_hdia_mref(const uint16_t* const ref_abv,
                        const uint16_t* const ref_lft, uint16_t* const dst,
                        ptrdiff_t dst_stride, int log2_pb_w,
                        int log2_pb_h, uint8_t mrl_idx)
{

    uint16_t tmp_dst[128 * 128];
    const int tmp_stride = 128;
    uint16_t* _tmp = tmp_dst;
    uint16_t* _dst = dst;
    int width = 1 << log2_pb_w;
    int height = 1 << log2_pb_h;

    int delta_pos = 32 * (mrl_idx + 1); // angle val is 32 for diag

    // and delta int is incremented of 1 (==y)
    for (int y = 0; y < width; y++) {
        const int delta_int = delta_pos >> 5;
        for (int x = 0; x < height; x++) {
            _tmp[x] = ref_lft[x + delta_int + 1];
        }
        delta_pos += 32;
        _tmp += tmp_stride;
    }

    _tmp = tmp_dst;
    // Transpose block
    for (int y = 0; y < width; y++) {
        _dst = &dst[y];
        for (int x = 0; x < height; x++) {
            _dst[0] = _tmp[x];
            _dst += dst_stride;
        }
        _tmp += tmp_stride;
    }
}

void
intra_angular_vdia_mref(const uint16_t* const ref_abv,
                        const uint16_t* const ref_lft, uint16_t* const dst,
                        ptrdiff_t dst_stride, int log2_pb_w,
                        int log2_pb_h, uint8_t mrl_idx)
{
    int delta_pos = 32 * (mrl_idx + 1); // angle val = 32 for strict diag
    uint16_t* _dst = dst;
    int width = 1 << log2_pb_w;
    int height = 1 << log2_pb_h;

    for (int y = 0; y < height; y++) {
        const int delta_int = delta_pos >> 5;

        for (int x = 0; x < width; x++) {
            _dst[x] = ref_abv[x + delta_int + 1];
        }

        delta_pos += 32;
        _dst += dst_stride;
    }
}

void
intra_angular_h_cubic_mref(const uint16_t* const ref_lft, uint16_t* const dst,
                           ptrdiff_t dst_stride,
                           int log2_pb_w, int log2_pb_h,
                           int angle_val, uint8_t mrl_idx)
{
    uint16_t tmp_dst[128 * 128];
    const int tmp_stride = 128;
    uint16_t* _tmp = tmp_dst;
    uint16_t* _dst = dst;
    int width = 1 << log2_pb_w;
    int height = 1 << log2_pb_h;
    int delta_pos = angle_val * (mrl_idx + 1);

    if ((angle_val & 0x1F)) {
        for (int y = 0; y < width; y++) {
            const int delta_int  = delta_pos >> 5;
            const int delta_frac = delta_pos & 0x1F;
            const int16_t* ref = (int16_t *)ref_lft + delta_int;
            const int8_t* filter = &chroma_filter[delta_frac << 2];
            for (int x = 0; x < height;
                 x++) { // FIXME check boundary
                int val =
                    ((int32_t)(ref[0] * filter[0]) +
                     (int32_t)(ref[1] * filter[1]) +
                     (int32_t)(ref[2] * filter[2]) +
                     (int32_t)(ref[3] * filter[3]) + 32) >>
                    6;
                ref++;
                _tmp[x] = ov_clip(val, 0, 1023);
            }
            delta_pos += angle_val;
            _tmp += tmp_stride;
        }
    } else {
        for (int y = 0; y < width; y++) {
            const int delta_int = delta_pos >> 5;
            for (int x = 0; x < height; x++) {
                _tmp[x] = ref_lft[x + delta_int + 1];
            }
            delta_pos += angle_val;
            _tmp += tmp_stride;
        }
    }

    _tmp = tmp_dst;
    // Transpose block
    for (int y = 0; y < width; y++) {
        _dst = &dst[y];
        for (int x = 0; x < height; x++) {
            _dst[0] = _tmp[x];
            _dst += dst_stride;
        }
        _tmp += tmp_stride;
    }
}

void
intra_angular_v_cubic_mref(const uint16_t* const ref_abv, uint16_t* const dst,
                           int dst_stride, int log2_pb_w,
                           int log2_pb_h, int angle_val,
                           uint8_t mrl_idx)
{
    int delta_pos = angle_val * (mrl_idx + 1);
    uint16_t* _dst = dst;
    int width = 1 << log2_pb_w;
    int height = 1 << log2_pb_h;
    // Note this can be vectorized
    if ((angle_val & 0x1F)) {
        for (int y = 0; y < height; y++) {
            const int delta_int = delta_pos >> 5;
            const int delta_frac = delta_pos & 0x1F;
            const int16_t* ref = (int16_t *)ref_abv + delta_int;
            const int8_t* filter = &chroma_filter[delta_frac << 2];
            for (int x = 0; x < width; x++) {
                int val =
                    ((int32_t)(ref[0] * filter[0]) +
                     (int32_t)(ref[1] * filter[1]) +
                     (int32_t)(ref[2] * filter[2]) +
                     (int32_t)(ref[3] * filter[3]) + 32) >>
                    6;
                ref++;
                _dst[x] = ov_clip(val, 0, 1023);
            }
            delta_pos += angle_val;
            _dst += dst_stride;
        }
    } else {
        for (int y = 0; y < height; y++) {
            const int delta_int = delta_pos >> 5;
            for (int x = 0; x < width; x++) {
                _dst[x] = ref_abv[x + delta_int + 1];
            }
            delta_pos += angle_val;
            _dst += dst_stride;
        }
    }
}
