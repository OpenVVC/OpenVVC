#include <stddef.h>
#include <stdint.h>

#include "ovutils.h"

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
vvc_intra_dc(const uint16_t* const src_above, const uint16_t* const src_left,
             uint16_t* const dst, ptrdiff_t dst_stride, int log2_pb_width,
             int log2_pb_height)
{
        int idx;
        uint16_t* _dst = dst;
        int dc_value = 0;

        const int denom = (log2_pb_width == log2_pb_height)
                            ? (log2_pb_width + 1)
                            : OVMAX(log2_pb_width, log2_pb_height);
        const int shift = denom;
        const int offset = ((1 << denom) >> 1);

        if (log2_pb_width >= log2_pb_height) {
                for (idx = 0; idx < (1 << log2_pb_width); idx++) {
                        dc_value += src_above[1 + idx];
                }
        }
        if (log2_pb_width <= log2_pb_height) {
                for (idx = 0; idx < (1 << log2_pb_height); idx++) {
                        dc_value += src_left[1 + idx];
                }
        }

        dc_value = (dc_value + offset) >> shift;

        for (int i = 0; i < (1 << log2_pb_height); ++i) {
                for (int j = 0; j < (1 << log2_pb_width); j++) {
                        _dst[j] = dc_value;
                }
                _dst += dst_stride;
        }
}

void
vvc_intra_planar(const uint16_t* const src_above,
                 const uint16_t* const src_left, uint16_t* const dst,
                 ptrdiff_t dst_stride, int log2_pb_width, int log2_pb_height)
{

    uint16_t* _dst = dst; // TODO template according to bitdepth;
    const uint32_t width = 1 << log2_pb_width;
    const uint32_t height = 1 << log2_pb_height;
    const uint32_t shift = 1 + log2_pb_width + log2_pb_height;
    const uint32_t offset = 1 << (log2_pb_width + log2_pb_height);
    int value;

    int top_row[128], bottom_row[128];
    int top_right = src_above[width + 1];

    for (int j = 0; j < width + 1; j++) {
        bottom_row[j] = src_left[height + 1];
    }

    for (int k = 0; k < width; k++) {
        value = src_above[k + 1];
        bottom_row[k] -= value; // bottom_left - val
        top_row[k] = value << log2_pb_height;
    }

    for (int y = 0; y < height; y++) {
        int value = (int)src_left[y + 1]; // left_value y
        int hor_pred = value << log2_pb_width;
        int right_pred = top_right - value;

        for (int x = 0; x < width; x++) {
            int vertPred;
            hor_pred += right_pred;
            top_row[x] += bottom_row[x];

            vertPred = top_row[x];

            _dst[x] = ((hor_pred << log2_pb_height) +
                       (vertPred << log2_pb_width) + offset) >>
                shift;
        }
        _dst += dst_stride;
    }
}

void
vvc_intra_dc_pdpc(const uint16_t* const src_above,
                  const uint16_t* const src_left, uint16_t* const dst,
                  ptrdiff_t dst_stride, int log2_pb_w, int log2_pb_h)
{
        int idx;
        uint16_t* _dst = dst;
        uint32_t dc_val = 0;

        const int denom = (log2_pb_w == log2_pb_h)
                            ? (log2_pb_w + 1)
                            : OVMAX(log2_pb_w, log2_pb_h);
        const int shift = denom;
        const int offset = ((1 << denom) >> 1);
        const int pdpc_scale = (log2_pb_w + log2_pb_h - 2) >> 2;
        const uint8_t* pdpc_w = vvc_pdpc_w[pdpc_scale];

        if (log2_pb_w >= log2_pb_h) {
                for (idx = 0; idx < (1 << log2_pb_w); idx++) {
                        dc_val += src_above[1 + idx];
                }
        }
        if (log2_pb_w <= log2_pb_h) {
                for (idx = 0; idx < (1 << log2_pb_h); idx++) {
                        dc_val += src_left[1 + idx];
                }
        }

        dc_val = (dc_val + offset) >> shift;

        for (int y = 0; y < (1 << log2_pb_h); y++) {
                int x;
                int l_wgh = pdpc_w[y];
                const int16_t l_val = src_left[y + 1];
                for (x = 0; x < (1 << log2_pb_w); x++) {
                        const int16_t t_val = src_above[x + 1];
                        int t_wgh = pdpc_w[x];
                        int val = ((t_wgh * l_val) + (l_wgh * t_val) +
                                   (64 - (t_wgh + l_wgh)) * dc_val + 32) >>
                                  6;
                        _dst[x] = ov_clip(val, 0, 1023);
                }
                _dst += dst_stride;
        }
}

// #define CUT_PDPC 0
void
vvc_intra_planar_pdpc(const uint16_t* const src_above,
                      const uint16_t* const src_left, uint16_t* const dst,
                      ptrdiff_t dst_stride, int log2_pb_w, int log2_pb_h)
{

        uint16_t* _dst = dst; // TODO template according to bitdepth;
        const uint32_t width = 1 << log2_pb_w;
        const uint32_t height = 1 << log2_pb_h;
        const uint32_t w_scale = OVMAX(1, log2_pb_w);
        const uint32_t h_scale = OVMAX(1, log2_pb_h);
        const uint32_t s_shift = w_scale + h_scale + 1;
        const uint32_t offset = 1 << (w_scale + h_scale);

        const uint8_t pdpc_scale = (log2_pb_w + log2_pb_h - 2) >> 2;

        int32_t t_row[128], b_row[128];
        const int8_t* pdpc_w = (int8_t *)vvc_pdpc_w[pdpc_scale];
        const int16_t bl_val = src_left[height + 1];
        const int16_t tr_val = src_above[width + 1];
        // #if CUT_PDPC
        //         const int pdpc_stop_w = OVMIN(3 << pdpc_scale, width);
        //         const int pdpc_stop_h = OVMIN(3 << pdpc_scale, height);
        // #endif
        int y;

        for (y = 0; y < width; ++y) {
                int32_t t_val = (int32_t)src_above[y + 1];
                b_row[y] = bl_val - t_val;
                t_row[y] = t_val << h_scale;
        }

        // #if CUT_PDPC
        //         for (y = 0; y < pdpc_stop_h; ++y) {
        // #else
        for (y = 0; y < height; ++y) {
                // #endif
                int x;
                int32_t l_val = (int32_t)src_left[y + 1];
                int y_wgh = pdpc_w[y];
                int32_t r_col_y = tr_val - l_val;
                int32_t l_col_y = l_val << w_scale;
                // #if CUT_PDPC
                //                 for (x = 0; x < pdpc_stop_w; ++x) {
                // #else
                for (x = 0; x < width; ++x) {
                        // #endif
                        int32_t val;
                        int x_wgh = pdpc_w[x];
                        int32_t t_val = (int32_t)src_above[x + 1];
                        l_col_y += r_col_y;
                        t_row[x] += b_row[x];
                        val = ((l_col_y << h_scale) + (t_row[x] << w_scale) +
                               offset) >>
                              s_shift;
                        val = ((x_wgh * l_val) + (y_wgh * t_val) +
                               (64 - (x_wgh + y_wgh)) * val + 32) >>
                              6;
                        _dst[x] = ov_clip(val, 0, 1023);
                }
                // #if CUT_PDPC
                //                 for (; x < width; ++x) {
                //                         int32_t val;
                //                         int32_t t_val = (int32_t)src_above[x
                //                         + 1]; l_col_y += r_col_y; t_row[x] +=
                //                         b_row[x]; val = ((l_col_y << h_scale)
                //                         + (t_row[x] << w_scale) +
                //                                offset) >>
                //                               s_shift;
                //                         val = ((y_wgh * t_val) + (64 - y_wgh)
                //                         * val + 32) >> 6; _dst[x] =
                //                         av_clip(val, 0, 1023);
                //                 }
                // #endif
                _dst += dst_stride;
        }

        // #if CUT_PDPC
        //         for (; y < height; ++y) {
        //                 int x;
        //                 int32_t l_val = (int32_t)src_left[y + 1];
        //                 int32_t r_col_y = tr_val - l_val;
        //                 int32_t l_col_y = l_val << w_scale;
        //                 for (x = 0; x < pdpc_stop_w; ++x) {
        //                         int32_t val;
        //                         int x_wgh = pdpc_w[x];
        //                         int32_t t_val = (int32_t)src_above[x + 1];
        //                         l_col_y += r_col_y;
        //                         t_row[x] += b_row[x];
        //                         val = ((l_col_y << h_scale) + (t_row[x] <<
        //                         w_scale) +
        //                                offset) >>
        //                               s_shift;
        //                         val = ((x_wgh * l_val) + (64 - x_wgh) * val +
        //                         32) >> 6; _dst[x] = av_clip(val, 0, 1023);
        //                 }
        //                 for (; x < width; ++x) {
        //                         int32_t val;
        //                         int32_t t_val = (int32_t)src_above[x + 1];
        //                         l_col_y += r_col_y;
        //                         t_row[x] += b_row[x];
        //                         val = ((l_col_y << h_scale) + (t_row[x] <<
        //                         w_scale) +
        //                                offset) >>
        //                               s_shift;
        //                         _dst[x] = av_clip(val, 0, 1023);
        //                 }
        //                 _dst += dst_stride;
        //         }
        // #endif
}

void
rcn_init_dc_planar_functions(struct RCNFunctions *const rcn_funcs)
{
    rcn_funcs->dc.func = &vvc_intra_dc;
    rcn_funcs->dc.pdpc = &vvc_intra_dc_pdpc;

    rcn_funcs->planar.func = &vvc_intra_planar;
    rcn_funcs->planar.pdpc[0] = &vvc_intra_planar_pdpc;
    rcn_funcs->planar.pdpc[1] = &vvc_intra_planar_pdpc;
}
