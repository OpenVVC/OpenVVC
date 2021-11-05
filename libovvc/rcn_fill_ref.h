#ifndef RCN_FILL_REF_H
#define RCN_FILL_REF_H

#include <stdint.h>

#define LOG2_UNIT_S 2

static inline uint64_t
non_available_units_map(uint64_t ref_map, uint8_t pos, uint8_t log2_pb_s)
{
    int pos_unit = pos >> LOG2_UNIT_S;
    int nb_ref_units = ((1 << (log2_pb_s + 1)) >> LOG2_UNIT_S) + 1;

    uint64_t req_unit_msk = (1llu << (nb_ref_units + 1)) - 1;
    uint64_t avl_map  = (ref_map >> pos_unit) & req_unit_msk;
    uint64_t navl_map = avl_map ^ req_unit_msk;

    return navl_map;
}

static inline uint64_t
available_units_map(uint64_t ref_map, uint8_t pos, uint8_t log2_pb_s)
{
    int pos_unit = pos >> LOG2_UNIT_S;
    int nb_ref_units = ((1 << (log2_pb_s + 1)) >> LOG2_UNIT_S) + 1;
    uint64_t req_unit_msk = (1llu << (nb_ref_units + 1)) - 1;

    uint64_t avl_map = (ref_map >> pos_unit) & req_unit_msk;

    return avl_map;
}

void
filter_ref_samples(const uint16_t* const src, uint16_t* const dst,
                   const uint16_t* src2, int length);

void
fill_ref_left_0(const uint16_t* const src, int src_stride,
                uint16_t* const ref_left, uint64_t intra_map_cols,
                uint64_t intra_map_rows, int8_t x0, int8_t y0, int log2_pb_w,
                int log2_pb_h, int offset_y);

void
fill_ref_left_0_chroma(const uint16_t* const src, int src_stride,
                       uint16_t* const ref_left, uint64_t intra_map_cols,
                       uint64_t intra_map_rows, int8_t x0, int8_t y0,
                       int log2_pb_w, int log2_pb_h);

void
fill_ref_left_0_mref(const uint16_t* const src, int src_stride,
                     uint16_t* const ref_left, uint64_t intra_map_cols,
                     uint64_t intra_map_rows, int mref_idx, int8_t x0,
                     int8_t y0, int log2_pb_w, int log2_pb_h);

void
fill_ref_above_0(const uint16_t* const src, int src_stride,
                 uint16_t* const ref_above, uint64_t intra_map_rows,
                 uint64_t intra_map_cols, int8_t x0, int8_t y0, int log2_pb_w,
                 int log2_pb_h, int offset_x);

void
fill_ref_above_0_chroma(const uint16_t* const src, int src_stride,
                        uint16_t* const ref_above, uint64_t intra_map_rows,
                        uint64_t intra_map_cols, int8_t x0, int8_t y0,
                        int log2_pb_w, int log2_pb_h);

void
fill_ref_above_0_mref(const uint16_t* const src, int src_stride,
                      uint16_t* const ref_above, uint64_t intra_map_rows,
                      uint64_t intra_map_cols, int mref_idx, int8_t x0,
                      int8_t y0, int log2_pb_w, int log2_pb_h);
#endif // RCN_FILL_REF_H
