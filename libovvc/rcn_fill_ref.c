/**
 *
 *   OpenVVC is open-source real time software decoder compliant with the 
 *   ITU-T H.266- MPEG-I - Part 3 VVC standard. OpenVVC is developed from 
 *   scratch in C as a library that provides consumers with real time and
 *   energy-aware decoding capabilities under different OS including MAC OS,
 *   Windows, Linux and Android targeting low energy real-time decoding of
 *   4K VVC videos on Intel x86 and ARM platforms.
 * 
 *   Copyright (C) 2020-2022  IETR-INSA Rennes :
 *   
 *   Pierre-Loup CABARAT
 *   Wassim HAMIDOUCHE
 *   Guillaume GAUTIER
 *   Thomas AMESTOY
 *
 *   This library is free software; you can redistribute it and/or
 *   modify it under the terms of the GNU Lesser General Public
 *   License as published by the Free Software Foundation; either
 *   version 2.1 of the License, or (at your option) any later version.
 *
 *   This library is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *   Lesser General Public License for more details.
 *
 *   You should have received a copy of the GNU Lesser General Public
 *   License along with this library; if not, write to the Free Software
 *   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301
 *   USA
 * 
 **/

#include <stdint.h>
#include <string.h>

#include "ovutils.h"
#include "rcn_fill_ref.h"
#include "ctudec.h"

#include "bitdepth.h"

// note src2 is only usefull for top_left value filtering
// WARNING src and dst cannot be aliased
// FIXME we could probably merge this part with ref filling
// FIXME length should be an uint
static void
filter_ref_samples(const OVSample* const src, OVSample* const dst,
                   const OVSample* src2, int length)
{

    // Regular reference sample filter
    const OVSample* _src = src;
    OVSample* _dst = dst;

    // WARNING top left uses above and left values top-left
    // FIXME: implicit conversion
    *_dst = (src2[1] + (_src[0] << 1) + _src[1] + 2) >> 2;
    _dst++;
    _src++;
    // top row (left-to-right)
    for (uint32_t i = 1; i < length; i++, _dst++, _src++) {
        // FIXME: implicit conversion
        *_dst = (_src[1] + (_src[0] << 1) + _src[-1] + 2) >> 2;
    }
    *_dst = *_src;
}

// Filling a map of available neighbours for intra process
static void
fill_ref_left_0(const OVSample* const src, int src_stride,
                OVSample* const ref_left, uint64_t intra_map_cols,
                uint64_t intra_map_rows, int8_t x0, int8_t y0, int log2_pb_w,
                int log2_pb_h, int offset_y)
{
    const OVSample* src_col = &src[(x0 - 1)];
    uint64_t  avl_map_l =     available_units_map(intra_map_cols, y0, log2_pb_h);
    uint64_t navl_map_l = non_available_units_map(intra_map_cols, y0, log2_pb_h);
    const int ref_length_l = (1 << (log2_pb_h + 1)) + 1;
    int nb_pb_ref_l = ((1 << (log2_pb_h + 1)) >> 2) + 1;

    src_col += (y0 - 1) * src_stride;

    if (!navl_map_l) {
        int i;
        for (i = 0; i < ref_length_l; ++i) {
            ref_left[i] = *src_col;
            src_col += src_stride;
        }

    } else if (avl_map_l) {
        int nb_pb_avl = 64 - __builtin_clzll(avl_map_l);
        OVSample padding_val = AVG_VAL;
        OVSample* _dst = ref_left;
        int i;

        if (avl_map_l & 0x1) {
            *_dst = *(src_col + src_stride * offset_y);
        } else {
            *_dst = src_col[src_stride];
        }

        src_col += src_stride;
        ++_dst;

        for (i = 1; i < nb_pb_avl; ++i) {
            _dst[0] = src_col[0 * src_stride];
            _dst[1] = src_col[1 * src_stride];
            _dst[2] = src_col[2 * src_stride];
            _dst[3] = src_col[3 * src_stride];
            src_col += 4 * src_stride;
            _dst += 4;
        }

        padding_val = *(_dst - 1);

        for (; i < nb_pb_ref_l; ++i) {
            _dst[0] = padding_val;
            _dst[1] = padding_val;
            _dst[2] = padding_val;
            _dst[3] = padding_val;
            _dst += 4;
        }

    } else {
        /* Pad with first available sample in above ref */
        uint64_t avl_map_a = available_units_map(intra_map_rows, x0, log2_pb_w);

        OVSample padding_val = AVG_VAL;
        int i;

        if (avl_map_a) {
            padding_val = src_col[1 + src_stride * offset_y];
        }

        for (i = 0; i < ref_length_l; ++i) {
            ref_left[i] = padding_val;
        }
    }

    /* Padding for wide angle */
    for (int i = 0; i < 4 + offset_y; ++i) {
        ref_left[ref_length_l + i] = ref_left[ref_length_l - 1];
    }
}

static void
fill_ref_left_0_chroma(const OVSample* const src, int src_stride,
                       OVSample* const ref_left, uint64_t intra_map_cols,
                       uint64_t intra_map_rows, int8_t x0, int8_t y0,
                       int log2_pb_w, int log2_pb_h)
{
    const OVSample* src_col = &src[(x0 - 1)];
    int y_pb = y0 >> 1;
    int nb_pb_ref_l = ((1 << (log2_pb_h + 1)) >> 1) + 1;

    uint64_t ref_map_l = (1llu << (nb_pb_ref_l + 1)) - 1;
    uint64_t avl_map_l = (intra_map_cols >> y_pb) & ref_map_l;
    uint64_t navl_map_l = avl_map_l ^ ref_map_l;

    src_col += (y0 - 1) * src_stride;

    if (!navl_map_l) {
        const int ref_length_l = (1 << (log2_pb_h + 1)) + 1;
        int i;
        for (i = 0; i < ref_length_l; ++i) {
            ref_left[i] = *src_col;
            src_col += src_stride;
        }
    } else if (avl_map_l) {
        int nb_pb_avl = 64 - __builtin_clzll(avl_map_l);
        OVSample padding_val = AVG_VAL;
        OVSample* _dst = ref_left;
        int i;

        if (avl_map_l & 0x1) {
            *_dst = *src_col;
        } else {
            *_dst = src_col[src_stride];
        }

        src_col += src_stride;
        ++_dst;

        for (i = 1; i < nb_pb_avl; ++i) {
            _dst[0] = src_col[0 * src_stride];
            _dst[1] = src_col[1 * src_stride];
            padding_val = src_col[1 * src_stride];
            src_col += 2 * src_stride;
            _dst += 2;
        }

        for (; i < nb_pb_ref_l; ++i) {
            _dst[0] = padding_val;
            _dst[1] = padding_val;
            _dst += 2;
        }

    } else {
        /* Pad with first available sample in above ref */
        int x_pb = x0 >> 1;
        int nb_pb_ref_a = ((1 << (log2_pb_w + 1)) >> 1) + 1;

        uint64_t ref_map_a = (1llu << (nb_pb_ref_a + 1)) - 1;
        uint64_t avl_map_a = (intra_map_rows >> x_pb) & ref_map_a;

        const int ref_length_l = (1 << (log2_pb_h + 1)) + 1;
        OVSample padding_val = AVG_VAL;
        int i;

        if (avl_map_a) {
            padding_val = src_col[1];
        }

        for (i = 0; i < ref_length_l; ++i) {
            ref_left[i] = padding_val;
        }
    }

    /* Padding for wide angle */
    for (int i = 0; i < 2; ++i) {
        ref_left[(1 << (log2_pb_h + 1)) + 1 + i] =
            ref_left[(1 << (log2_pb_h + 1)) + i];
    }
}

static void
fill_ref_left_0_mref(const OVSample* const src, int src_stride,
                     OVSample* const ref_left, uint64_t intra_map_cols,
                     uint64_t intra_map_rows, int mref_idx, int8_t x0,
                     int8_t y0, int log2_pb_w, int log2_pb_h)
{
    const OVSample* src_col = &src[(x0 - (mref_idx + 1))];
    int y_pb = y0 >> 2;
    int nb_pb_ref_l = ((1 << (log2_pb_h + 1)) >> 2) + 1;
    int hw_ratio =
        OVMAX(1, (1 << log2_pb_h) >> log2_pb_w); // Note this is only

    uint64_t ref_map_l = (1llu << (nb_pb_ref_l + 1)) - 1;
    uint64_t avl_map_l = (intra_map_cols >> y_pb) & ref_map_l;
    uint64_t navl_map_l = avl_map_l ^ ref_map_l;

    src_col += (y0 - (1 + mref_idx)) * src_stride;

    if (!navl_map_l) {
        const int ref_length_l =
            (1 << (log2_pb_h + 1)) + 1 + (mref_idx + 1);
        int i;
        for (i = 0; i < ref_length_l; ++i) {
            ref_left[i] = *src_col;
            src_col += src_stride;
        }
    } else if (avl_map_l) {
        int nb_pb_avl = 64 - __builtin_clzll(avl_map_l);
        OVSample padding_val = AVG_VAL;
        OVSample* _dst = ref_left;
        int i;

        if (avl_map_l & 0x1) {
            for (i = 0; i < mref_idx + 1; ++i) {
                *_dst = *(src_col);
                _dst++;
                src_col += src_stride;
            }
        } else {
            src_col += src_stride;
            for (i = 0; i < mref_idx + 1; ++i) {
                *_dst = *src_col;
                _dst++;
            }
        }

        for (i = 1; i < nb_pb_avl; ++i) {
            _dst[0] = src_col[0 * src_stride];
            _dst[1] = src_col[1 * src_stride];
            _dst[2] = src_col[2 * src_stride];
            _dst[3] = src_col[3 * src_stride];
            padding_val = src_col[3 * src_stride];
            src_col += 4 * src_stride;
            _dst += 4;
        }

        for (; i < nb_pb_ref_l; ++i) {
            _dst[0] = padding_val;
            _dst[1] = padding_val;
            _dst[2] = padding_val;
            _dst[3] = padding_val;
            _dst += 4;
        }

    } else {
        /* Pad with first available sample in above ref */
        int x_pb = x0 >> 2;
        int nb_pb_ref_a = ((1 << (log2_pb_w + 1)) >> 2) + 1;

        uint64_t ref_map_a = (1llu << (nb_pb_ref_a + 1)) - 1;
        uint64_t avl_map_a = (intra_map_rows >> x_pb) & ref_map_a;

        const int ref_length_l =
            (1 << (log2_pb_h + 1)) + 1 + (mref_idx + 1);
        OVSample padding_val = AVG_VAL;
        int i;

        if (avl_map_a) {
            padding_val = src_col[1 + mref_idx];
        }

        for (i = 0; i < ref_length_l; ++i) {
            ref_left[i] = padding_val;
        }
    }

    /* Padding for wide angle */
    for (int i = 0; i < hw_ratio * (mref_idx + 1); ++i) {
        ref_left[(1 << (log2_pb_h + 1)) + (mref_idx + 1) + i] =
            ref_left[(1 << (log2_pb_h + 1)) + (mref_idx) + i];
    }
}

static void
fill_ref_above_0(const OVSample* const src, int src_stride,
                 OVSample* const ref_above, uint64_t intra_map_rows,
                 uint64_t intra_map_cols, int8_t x0, int8_t y0, int log2_pb_w,
                 int log2_pb_h, int offset_x)
{
    const OVSample *src_line = &src[(y0 - 1) * src_stride];

    uint64_t  avl_map_a =     available_units_map(intra_map_rows, x0, log2_pb_w);
    uint64_t navl_map_a = non_available_units_map(intra_map_rows, x0, log2_pb_w);
    const int ref_length_a = (1 << (log2_pb_w + 1)) + 1;

    src_line += (x0 - 1);

    if (!navl_map_a) {
        int i;
        for (i = 0; i < ref_length_a; ++i) {
            ref_above[i] = *src_line;
            ++src_line;
        }
    } else {
        OVSample padding_value = AVG_VAL;

        if (avl_map_a) {
            OVSample *_dst = ref_above + 1;
            int nb_pb_avl = 64 - __builtin_clzll(avl_map_a);

            memcpy(_dst, src_line + 1, (nb_pb_avl - 1) * (sizeof(*_dst) << LOG2_UNIT_S));
            _dst += (nb_pb_avl - 1) << LOG2_UNIT_S;

            if (avl_map_a & 0x1) {
                ref_above[0] = *(src_line + offset_x);
            } else {
                ref_above[0] = src_line[1];
            }

            padding_value = _dst[-1];

            navl_map_a >>= nb_pb_avl;

            while (navl_map_a) {
                _dst[0] = padding_value;
                _dst[1] = padding_value;
                _dst[2] = padding_value;
                _dst[3] = padding_value;
                navl_map_a >>= 1;
                _dst += 4;
            }

        } else {
            /* Pad with first available left ref sample value */
            uint64_t avl_map_l = available_units_map(intra_map_cols, y0, log2_pb_h);
            int i;

            padding_value = AVG_VAL;

            if (avl_map_l) {
                src_line += offset_x + src_stride;
                padding_value = *src_line;
            }

            for (i = 0; i < ref_length_a; ++i) {
                ref_above[i] = padding_value;
            }
        }
    }

    /* Padding for wide angle */
    for (int i = 0; i < 4 + offset_x; ++i) {
        ref_above[ref_length_a + i] = ref_above[ref_length_a - 1 + i];
    }
}

static void
fill_ref_above_0_chroma(const OVSample* const src, int src_stride,
                        OVSample* const ref_above, uint64_t intra_map_rows,
                        uint64_t intra_map_cols, int8_t x0, int8_t y0,
                        int log2_pb_w, int log2_pb_h)
{
    int x_pb = x0 >> 1;
    int nb_pb_ref_a = ((1 << (log2_pb_w + 1)) >> 1) + 1;

    uint64_t ref_map_a = (1llu << (nb_pb_ref_a + 1)) - 1;
    uint64_t avl_map_a = (intra_map_rows >> x_pb) & ref_map_a;
    uint64_t navl_map_a = avl_map_a ^ ref_map_a;

    const OVSample* src_line = &src[(y0 - 1) * src_stride];

    src_line += (x0 - 1);

    if (!navl_map_a) {
        const int ref_length_a = (1 << (log2_pb_w + 1)) + 1;
        int i;
        for (i = 0; i < ref_length_a; ++i) {
            ref_above[i] = *src_line;
            ++src_line;
        }
    } else {
        OVSample padding_value = AVG_VAL;
        if (avl_map_a) {

            // FIXME: int nb_pb_usable = 64 -
            // __builtin_clzll(avl_map_a);
            // FIXME: int nb_pb_missing = nb_pb_ref_a -
            // nb_pb_usable;
            OVSample* _dst = ref_above;

            if (avl_map_a & 0x1) {
                *_dst = src_line[0];
            } else {
                *_dst = src_line[1];
            }

            ++_dst;
            ++src_line;
            avl_map_a >>= 1;
            navl_map_a >>= 1;

            while (avl_map_a) {
                _dst[0] = src_line[0];
                _dst[1] = src_line[1];
                avl_map_a >>= 1;
                navl_map_a >>= 1;
               src_line += 2;
                _dst += 2;
            }

            padding_value = src_line[-1];

            while (navl_map_a) {
                _dst[0] = padding_value;
                _dst[1] = padding_value;
                navl_map_a >>= 1;
                _dst += 2;
            }

        } else {
            /* Pad with first available left ref sample value */
            int y_pb = y0 >> 1;
            const int ref_length_a = (1 << (log2_pb_w + 1)) + 1;
            int nb_pb_ref_l = ((1 << (log2_pb_h + 1)) >> 1) + 1;

            uint64_t needed_mask_l =
                (1llu << (nb_pb_ref_l + 1)) - 1;
            uint64_t usable_mask_l =
                (intra_map_cols >> y_pb) & needed_mask_l;

            int i;

            padding_value = AVG_VAL;

            if (usable_mask_l) {
                padding_value = src_line[src_stride];
            }

            for (i = 0; i < ref_length_a; ++i) {
                ref_above[i] = padding_value;
            }
        }
    }

    /* Padding for wide angle */
    for (int i = 0; i < 4; ++i) {
        ref_above[(1 << (log2_pb_w + 1)) + 1 + i] =
            ref_above[(1 << (log2_pb_w + 1)) + i];
    }
}

static void
fill_ref_above_0_mref(const OVSample* const src, int src_stride,
                      OVSample* const ref_above, uint64_t intra_map_rows,
                      uint64_t intra_map_cols, int mref_idx, int8_t x0,
                      int8_t y0, int log2_pb_w, int log2_pb_h)
{
    int x_pb = x0 >> 2;
    int nb_pb_ref_a = ((1 << (log2_pb_w + 1)) >> 2) + 1;
    int wh_ratio = OVMAX(1, (1 << log2_pb_w) >> log2_pb_h);

    uint64_t ref_map_a = (1llu << (nb_pb_ref_a + 1)) - 1;
    uint64_t avl_map_a = (intra_map_rows >> x_pb) & ref_map_a;
    uint64_t navl_map_a = avl_map_a ^ ref_map_a;

    const OVSample* src_line = &src[(y0 - (1 + mref_idx)) * src_stride];

    src_line += (x0 - (1 + mref_idx));

    if (!navl_map_a) {
        const int ref_length_a =
            (1 << (log2_pb_w + 1)) + 1 + (mref_idx + 1);
        int i;
        for (i = 0; i < ref_length_a; ++i) {
            ref_above[i] = *src_line;
            ++src_line;
        }
    } else {
        OVSample padding_value = AVG_VAL;
        if (avl_map_a) {

            // FIXME: int nb_pb_usable = 64 -
            // __builtin_clzll(avl_map_a);
            // FIXME: int nb_pb_missing = nb_pb_ref_a -
            // nb_pb_usable;
            OVSample* _dst = ref_above;
            int i;

            if (avl_map_a & 0x1) {
                for (i = 0; i < mref_idx + 1; ++i) {
                    *_dst = *(src_line);
                    ++src_line;
                    ++_dst;
                }
            } else {
                src_line += (mref_idx + 1);
                for (i = 0; i < mref_idx + 1; ++i) {
                    *_dst = *(src_line);
                    ++_dst;
                }
            }

            avl_map_a >>= 1;
            navl_map_a >>= 1;

            while (avl_map_a) {
                _dst[0] = src_line[0];
                _dst[1] = src_line[1];
                _dst[2] = src_line[2];
                _dst[3] = src_line[3];
                avl_map_a >>= 1;
                navl_map_a >>= 1;
                src_line += 4;
                _dst += 4;
            }

            padding_value = src_line[-1];

            while (navl_map_a) {
                _dst[0] = padding_value;
                _dst[1] = padding_value;
                _dst[2] = padding_value;
                _dst[3] = padding_value;
                navl_map_a >>= 1;
                _dst += 4;
            }

        } else {
            /* Pad with first available left ref sample value */
            int y_pb = y0 >> 2;
            const int ref_length_a = (1 << (log2_pb_w + 1)) + 1;
            // FIXME: +(mref_idx + 1);
            int nb_pb_ref_l = ((1 << (log2_pb_h + 1)) >> 2) + 1;

            uint64_t needed_mask_l =
                (1llu << (nb_pb_ref_l + 1)) - 1;
            uint64_t usable_mask_l =
                (intra_map_cols >> y_pb) & needed_mask_l;

            int i;
            const OVSample* src_line = &src[y0 * src_stride];

            src_line += (x0 - 1);

            padding_value = AVG_VAL;

            if (usable_mask_l) {
                padding_value = *src_line;
            }

            for (i = 0; i < ref_length_a; ++i) {
                ref_above[i] = padding_value;
            }
        }
    }

    /* Padding for wide angle */
    for (int i = 0; i < (mref_idx + 1) * wh_ratio; ++i) {
        ref_above[(1 << (log2_pb_w + 1)) + (mref_idx + 1) + i] =
            ref_above[(1 << (log2_pb_w + 1)) + (mref_idx) + i];
    }
}

void
BD_DECL(rcn_init_fill_ref)(struct RCNFunctions *rcn_funcs)
{
    rcn_funcs->tmp.filter_ref_samples      = &filter_ref_samples;
    rcn_funcs->tmp.fill_ref_left_0         = &fill_ref_left_0;
    rcn_funcs->tmp.fill_ref_left_0_chroma  = &fill_ref_left_0_chroma;
    rcn_funcs->tmp.fill_ref_left_0_mref    = &fill_ref_left_0_mref;
    rcn_funcs->tmp.fill_ref_above_0        = &fill_ref_above_0;
    rcn_funcs->tmp.fill_ref_above_0_chroma = &fill_ref_above_0_chroma;
    rcn_funcs->tmp.fill_ref_above_0_mref   = &fill_ref_above_0_mref;
}
