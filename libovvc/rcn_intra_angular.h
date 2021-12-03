#ifndef RCN_INTRA_ANGULAR_H
#define RCN_INTRA_ANGULAR_H

#include <stddef.h>
#include <stdint.h>
#include "bitdepth.h"

struct IntraAngularFunctions {
    void (*pure)(const OVSample* ref_lft, OVSample* dst,
                 ptrdiff_t dst_stride, int8_t log2_pb_w,
                 int8_t log2_pb_h);

    void (*diagonal)(const OVSample* ref_lft, OVSample* dst,
                     ptrdiff_t dst_stride, int8_t log2_pb_w,
                     int8_t log2_pb_h);

    void (*angular)(const OVSample* ref_lft, OVSample* dst,
                    ptrdiff_t dst_stride, int8_t log2_pb_w,
                    int8_t log2_pb_h, int angle_val);

    void (*pure_pdpc)(const OVSample* ref_abv, const OVSample* ref_lft,
                      OVSample* const dst, ptrdiff_t dst_stride,
                      int8_t log2_pb_w, int8_t log2_pb_h);

    void (*diagonal_pdpc)(const OVSample* ref_abv, const OVSample* ref_lft,
                         OVSample* const dst, ptrdiff_t dst_stride,
                         int8_t log2_pb_w, int8_t log2_pb_h);

    void (*angular_pdpc)(const OVSample* ref_abv, const OVSample* ref_lft,
                         OVSample* const dst, ptrdiff_t dst_stride,
                         int8_t log2_pb_w, int8_t log2_pb_h, int mode_idx);
};

struct IntraMRLFunctions
{
    void (*angular_h)(const OVSample* const ref_lft, OVSample* const dst,
                      ptrdiff_t dst_stride,
                      int8_t log2_pb_w, int8_t log2_pb_h,
                      int angle_val, uint8_t multi_ref_idx);

    void (*angular_v)(const OVSample* const ref_lft, OVSample* const dst,
                      ptrdiff_t dst_stride,
                      int8_t log2_pb_w, int8_t log2_pb_h,
                      int angle_val, uint8_t multi_ref_idx);
};

#endif // RCN_INTRA_ANGULAR_H
