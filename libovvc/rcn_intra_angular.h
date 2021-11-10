#ifndef RCN_INTRA_ANGULAR_H
#define RCN_INTRA_ANGULAR_H

#include <stddef.h>
#include <stdint.h>

struct IntraAngularFunctions {
    void (*pure)(const uint16_t* ref_lft, uint16_t* dst,
                 ptrdiff_t dst_stride, int8_t log2_pb_w,
                 int8_t log2_pb_h);

    void (*diagonal)(const uint16_t* ref_lft, uint16_t* dst,
                     ptrdiff_t dst_stride, int8_t log2_pb_w,
                     int8_t log2_pb_h);

    void (*angular)(const uint16_t* ref_lft, uint16_t* dst,
                    ptrdiff_t dst_stride, int8_t log2_pb_w,
                    int8_t log2_pb_h, int angle_val);

    void (*pure_pdpc)(const uint16_t* ref_abv, const uint16_t* ref_lft,
                      uint16_t* const dst, ptrdiff_t dst_stride,
                      int8_t log2_pb_w, int8_t log2_pb_h);

    void (*diagonal_pdpc)(const uint16_t* ref_abv, const uint16_t* ref_lft,
                         uint16_t* const dst, ptrdiff_t dst_stride,
                         int8_t log2_pb_w, int8_t log2_pb_h);

    void (*angular_pdpc)(const uint16_t* ref_abv, const uint16_t* ref_lft,
                         uint16_t* const dst, ptrdiff_t dst_stride,
                         int8_t log2_pb_w, int8_t log2_pb_h, int mode_idx);
};

struct IntraMRLFunctions
{
    void (*angular_h)(const uint16_t* const ref_lft, uint16_t* const dst,
                      ptrdiff_t dst_stride,
                      int8_t log2_pb_w, int8_t log2_pb_h,
                      int angle_val, uint8_t multi_ref_idx);

    void (*angular_v)(const uint16_t* const ref_lft, uint16_t* const dst,
                      ptrdiff_t dst_stride,
                      int8_t log2_pb_w, int8_t log2_pb_h,
                      int angle_val, uint8_t multi_ref_idx);
};

#endif // RCN_INTRA_ANGULAR_H
