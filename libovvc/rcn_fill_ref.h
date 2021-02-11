#ifndef RCN_FILL_REF_H
#define RCN_FILL_REF_H

#include <stdint.h>

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
