#include "ovutils.h"
#include <stdint.h>

#include "rcn_fill_ref.h"

// note src2 is only usefull for top_left value filtering
// WARNING src and dst cannot be aliased
// FIXME we could probably merge this part with ref filling
// FIXME length should be an uint
void
filter_ref_samples(const uint16_t* const src, uint16_t* const dst,
                   const uint16_t* src2, int length)
{

        // Regular reference sample filter
        const uint16_t* _src = src;
        uint16_t* _dst = dst;

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
void
fill_ref_left_0(const uint16_t* const src, int src_stride,
                uint16_t* const ref_left, uint64_t intra_map_cols,
                uint64_t intra_map_rows, int8_t x0, int8_t y0, int log2_pb_w,
                int log2_pb_h, int offset_y)
{
        const uint16_t* _src = &src[(x0 - 1) + (y0 - 1) * src_stride];
        int y_pb = y0 >> 2;
        int nb_pb_ref_l = ((1 << (log2_pb_h + 1)) >> 2) + 1;

        uint64_t ref_map_l = (1llu << (nb_pb_ref_l + 1)) - 1;
        uint64_t avl_map_l = (intra_map_cols >> y_pb) & ref_map_l;
        uint64_t navl_map_l = avl_map_l ^ ref_map_l;

        if (!navl_map_l) {
                const int ref_length_l = (1 << (log2_pb_h + 1)) + 1;
                int i;
                for (i = 0; i < ref_length_l; ++i) {
                        ref_left[i] = *_src;
                        _src += src_stride;
                }
        } else if (avl_map_l) {
                int nb_pb_avl = 64 - __builtin_clzll(avl_map_l);
                // int nb_pb_navl = nb_pb_ref_l - nb_pb_avl;
                uint16_t padding_val = 512;
                uint16_t* _dst = ref_left;
                int i;

                if (avl_map_l & 0x1) {
                        *_dst = *(_src + src_stride * offset_y);
                } else {
                        *_dst = _src[src_stride];
                }

                _src += src_stride;
                ++_dst;

                for (i = 1; i < nb_pb_avl; ++i) {
                        _dst[0] = _src[0 * src_stride];
                        _dst[1] = _src[1 * src_stride];
                        _dst[2] = _src[2 * src_stride];
                        _dst[3] = _src[3 * src_stride];
                        padding_val = _src[3 * src_stride];
                        _src += 4 * src_stride;
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

                const int ref_length_l = (1 << (log2_pb_h + 1)) + 1;
                uint16_t padding_val = 512;
                int i;

                if (avl_map_a) {
                        padding_val = _src[1 + src_stride * offset_y];
                }

                for (i = 0; i < ref_length_l; ++i) {
                        ref_left[i] = padding_val;
                }
        }

        /* Padding for wide angle */
        for (int i = 0; i < 4 + offset_y; ++i) {
                ref_left[(1 << (log2_pb_h + 1)) + 1 + i] =
                  ref_left[(1 << (log2_pb_h + 1)) + i];
        }
}

void
fill_ref_left_0_chroma(const uint16_t* const src, int src_stride,
                       uint16_t* const ref_left, uint64_t intra_map_cols,
                       uint64_t intra_map_rows, int8_t x0, int8_t y0,
                       int log2_pb_w, int log2_pb_h)
{
        const uint16_t* _src = &src[(x0 - 1) + (y0 - 1) * src_stride];
        int y_pb = y0 >> 1;
        int nb_pb_ref_l = ((1 << (log2_pb_h + 1)) >> 1) + 1;

        uint64_t ref_map_l = (1llu << (nb_pb_ref_l + 1)) - 1;
        uint64_t avl_map_l = (intra_map_cols >> y_pb) & ref_map_l;
        uint64_t navl_map_l = avl_map_l ^ ref_map_l;

        if (!navl_map_l) {
                const int ref_length_l = (1 << (log2_pb_h + 1)) + 1;
                int i;
                for (i = 0; i < ref_length_l; ++i) {
                        ref_left[i] = *_src;
                        _src += src_stride;
                }
        } else if (avl_map_l) {
                int nb_pb_avl = 64 - __builtin_clzll(avl_map_l);
                // int nb_pb_navl = nb_pb_ref_l - nb_pb_avl;
                uint16_t padding_val = 512;
                uint16_t* _dst = ref_left;
                int i;

                if (avl_map_l & 0x1) {
                        *_dst = *_src;
                } else {
                        *_dst = _src[src_stride];
                }

                _src += src_stride;
                ++_dst;

                for (i = 1; i < nb_pb_avl; ++i) {
                        _dst[0] = _src[0 * src_stride];
                        _dst[1] = _src[1 * src_stride];
                        padding_val = _src[1 * src_stride];
                        _src += 2 * src_stride;
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
                uint16_t padding_val = 512;
                int i;

                if (avl_map_a) {
                        padding_val = _src[1];
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

void
fill_ref_left_0_mref(const uint16_t* const src, int src_stride,
                     uint16_t* const ref_left, uint64_t intra_map_cols,
                     uint64_t intra_map_rows, int mref_idx, int8_t x0,
                     int8_t y0, int log2_pb_w, int log2_pb_h)
{
        const uint16_t* _src =
          &src[(x0 - (mref_idx + 1)) + (y0 - (mref_idx + 1)) * src_stride];
        int y_pb = y0 >> 2;
        int nb_pb_ref_l = ((1 << (log2_pb_h + 1)) >> 2) + 1;
        int hw_ratio =
          OVMAX(1, (1 << log2_pb_h) >> log2_pb_w); // Note this is only

        uint64_t ref_map_l = (1llu << (nb_pb_ref_l + 1)) - 1;
        uint64_t avl_map_l = (intra_map_cols >> y_pb) & ref_map_l;
        uint64_t navl_map_l = avl_map_l ^ ref_map_l;

        if (!navl_map_l) {
                const int ref_length_l =
                  (1 << (log2_pb_h + 1)) + 1 + (mref_idx + 1);
                int i;
                for (i = 0; i < ref_length_l; ++i) {
                        ref_left[i] = *_src;
                        _src += src_stride;
                }
        } else if (avl_map_l) {
                int nb_pb_avl = 64 - __builtin_clzll(avl_map_l);
                // int nb_pb_navl = nb_pb_ref_l - nb_pb_avl;
                uint16_t padding_val = 512;
                uint16_t* _dst = ref_left;
                int i;

                if (avl_map_l & 0x1) {
                        for (i = 0; i < mref_idx + 1; ++i) {
                                *_dst = *(_src);
                                _dst++;
                                _src += src_stride;
                        }
                } else {
                        _src += src_stride;
                        for (i = 0; i < mref_idx + 1; ++i) {
                                *_dst = *_src;
                                _dst++;
                        }
                }

                for (i = 1; i < nb_pb_avl; ++i) {
                        _dst[0] = _src[0 * src_stride];
                        _dst[1] = _src[1 * src_stride];
                        _dst[2] = _src[2 * src_stride];
                        _dst[3] = _src[3 * src_stride];
                        padding_val = _src[3 * src_stride];
                        _src += 4 * src_stride;
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
                uint16_t padding_val = 512;
                int i;

                if (avl_map_a) {
                        padding_val = _src[1 + mref_idx];
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

void
fill_ref_above_0(const uint16_t* const src, int src_stride,
                 uint16_t* const ref_above, uint64_t intra_map_rows,
                 uint64_t intra_map_cols, int8_t x0, int8_t y0, int log2_pb_w,
                 int log2_pb_h, int offset_x)
{
        int x_pb = x0 >> 2;
        int nb_pb_ref_a = ((1 << (log2_pb_w + 1)) >> 2) + 1;

        uint64_t ref_map_a = (1llu << (nb_pb_ref_a + 1)) - 1;
        uint64_t avl_map_a = (intra_map_rows >> x_pb) & ref_map_a;
        uint64_t navl_map_a = avl_map_a ^ ref_map_a;

        const uint16_t* _src = &src[(x0 - 1) + (y0 - 1) * src_stride];

        if (!navl_map_a) {
                const int ref_length_a = (1 << (log2_pb_w + 1)) + 1;
                int i;
                for (i = 0; i < ref_length_a; ++i) {
                        ref_above[i] = *_src;
                        ++_src;
                }
        } else {
                uint16_t padding_value = 512;
                if (avl_map_a) {

                        // FIXME: int nb_pb_usable = 64 -
                        // __builtin_clzll(avl_map_a);
                        // FIXME: int nb_pb_missing = nb_pb_ref_a -
                        // nb_pb_usable;
                        uint16_t* _dst = ref_above + 1;

                        if (avl_map_a & 0x1) {
                                ref_above[0] = *(_src + offset_x);
                        } else {
                                ref_above[0] = _src[1];
                        }

                        ++_src;
                        avl_map_a >>= 1;
                        navl_map_a >>= 1;

                        while (avl_map_a) {
                                _dst[0] = _src[0];
                                _dst[1] = _src[1];
                                _dst[2] = _src[2];
                                _dst[3] = _src[3];
                                avl_map_a >>= 1;
                                navl_map_a >>= 1;
                                _src += 4;
                                _dst += 4;
                        }

                        padding_value = _src[-1];

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
                        int nb_pb_ref_l = ((1 << (log2_pb_h + 1)) >> 2) + 1;

                        uint64_t needed_mask_l =
                          (1llu << (nb_pb_ref_l + 1)) - 1;
                        uint64_t usable_mask_l =
                          (intra_map_cols >> y_pb) & needed_mask_l;

                        int i;
                        // FIXME: Redeclaration de *_src
                        const uint16_t* _src =
                          &src[(x0 + offset_x - 1) + y0 * src_stride];

                        padding_value = 512;

                        if (usable_mask_l) {
                                padding_value = *_src;
                        }

                        for (i = 0; i < ref_length_a; ++i) {
                                ref_above[i] = padding_value;
                        }
                }
        }

        /* Padding for wide angle */
        for (int i = 0; i < 4 + offset_x; ++i) {
                ref_above[(1 << (log2_pb_w + 1)) + 1 + i] =
                  ref_above[(1 << (log2_pb_w + 1)) + i];
        }
}

void
fill_ref_above_0_chroma(const uint16_t* const src, int src_stride,
                        uint16_t* const ref_above, uint64_t intra_map_rows,
                        uint64_t intra_map_cols, int8_t x0, int8_t y0,
                        int log2_pb_w, int log2_pb_h)
{
        int x_pb = x0 >> 1;
        int nb_pb_ref_a = ((1 << (log2_pb_w + 1)) >> 1) + 1;

        uint64_t ref_map_a = (1llu << (nb_pb_ref_a + 1)) - 1;
        uint64_t avl_map_a = (intra_map_rows >> x_pb) & ref_map_a;
        uint64_t navl_map_a = avl_map_a ^ ref_map_a;

        const uint16_t* _src = &src[(x0 - 1) + (y0 - 1) * src_stride];

        if (!navl_map_a) {
                const int ref_length_a = (1 << (log2_pb_w + 1)) + 1;
                int i;
                for (i = 0; i < ref_length_a; ++i) {
                        ref_above[i] = *_src;
                        ++_src;
                }
        } else {
                uint16_t padding_value = 512;
                if (avl_map_a) {

                        // FIXME: int nb_pb_usable = 64 -
                        // __builtin_clzll(avl_map_a);
                        // FIXME: int nb_pb_missing = nb_pb_ref_a -
                        // nb_pb_usable;
                        uint16_t* _dst = ref_above;

                        if (avl_map_a & 0x1) {
                                *_dst = _src[0];
                        } else {
                                *_dst = _src[1];
                        }

                        ++_dst;
                        ++_src;
                        avl_map_a >>= 1;
                        navl_map_a >>= 1;

                        while (avl_map_a) {
                                _dst[0] = _src[0];
                                _dst[1] = _src[1];
                                avl_map_a >>= 1;
                                navl_map_a >>= 1;
                                _src += 2;
                                _dst += 2;
                        }

                        padding_value = _src[-1];

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

                        padding_value = 512;

                        if (usable_mask_l) {
                                padding_value = _src[src_stride];
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

void
fill_ref_above_0_mref(const uint16_t* const src, int src_stride,
                      uint16_t* const ref_above, uint64_t intra_map_rows,
                      uint64_t intra_map_cols, int mref_idx, int8_t x0,
                      int8_t y0, int log2_pb_w, int log2_pb_h)
{
        int x_pb = x0 >> 2;
        int nb_pb_ref_a = ((1 << (log2_pb_w + 1)) >> 2) + 1;
        int wh_ratio = OVMAX(1, (1 << log2_pb_w) >> log2_pb_h);

        uint64_t ref_map_a = (1llu << (nb_pb_ref_a + 1)) - 1;
        uint64_t avl_map_a = (intra_map_rows >> x_pb) & ref_map_a;
        uint64_t navl_map_a = avl_map_a ^ ref_map_a;

        const uint16_t* _src =
          &src[(x0 - (1 + mref_idx)) + (y0 - (1 + mref_idx)) * src_stride];

        if (!navl_map_a) {
                const int ref_length_a =
                  (1 << (log2_pb_w + 1)) + 1 + (mref_idx + 1);
                int i;
                for (i = 0; i < ref_length_a; ++i) {
                        ref_above[i] = *_src;
                        ++_src;
                }
        } else {
                uint16_t padding_value = 512;
                if (avl_map_a) {

                        // FIXME: int nb_pb_usable = 64 -
                        // __builtin_clzll(avl_map_a);
                        // FIXME: int nb_pb_missing = nb_pb_ref_a -
                        // nb_pb_usable;
                        uint16_t* _dst = ref_above;
                        int i;

                        if (avl_map_a & 0x1) {
                                for (i = 0; i < mref_idx + 1; ++i) {
                                        *_dst = *(_src);
                                        ++_src;
                                        ++_dst;
                                }
                        } else {
                                _src += (mref_idx + 1);
                                for (i = 0; i < mref_idx + 1; ++i) {
                                        *_dst = *(_src);
                                        ++_dst;
                                }
                        }

                        avl_map_a >>= 1;
                        navl_map_a >>= 1;

                        while (avl_map_a) {
                                _dst[0] = _src[0];
                                _dst[1] = _src[1];
                                _dst[2] = _src[2];
                                _dst[3] = _src[3];
                                avl_map_a >>= 1;
                                navl_map_a >>= 1;
                                _src += 4;
                                _dst += 4;
                        }

                        padding_value = _src[-1];

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
                        const uint16_t* _src = &src[(x0 - 1) + y0 * src_stride];

                        padding_value = 512;

                        if (usable_mask_l) {
                                padding_value = *_src;
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
