#ifndef RCN_INTRA_ANGULAR_H
#define RCN_INTRA_ANGULAR_H

#include <stddef.h>
#include <stdint.h>

void
intra_angular_hdia_pdpc(const uint16_t* const ref_abv,
                        const uint16_t* const ref_lft, uint16_t* const dst,
                        ptrdiff_t dst_stride, int8_t log2_pb_w,
                        int8_t log2_pb_h);
void
intra_angular_vdia_pdpc(const uint16_t* const ref_abv,
                        const uint16_t* const ref_lft, uint16_t* const dst,
                        ptrdiff_t dst_stride, int8_t log2_pb_w,
                        int8_t log2_pb_h);

void
intra_angular_hdia(const uint16_t* const ref_lft, uint16_t* const dst,
                   ptrdiff_t dst_stride, int8_t log2_pb_w,
                   int8_t log2_pb_h);
void
intra_angular_vdia(const uint16_t* const ref_abv, uint16_t* const dst,
                   ptrdiff_t dst_stride, int8_t log2_pb_w,
                   int8_t log2_pb_h);
void
intra_angular_h_c(const uint16_t* ref_lft, uint16_t* dst,
                  ptrdiff_t dst_stride, int8_t log2_pb_w,
                  int8_t log2_pb_h, int angle_val);
void
intra_angular_v_c(const uint16_t* ref_abv, uint16_t* dst,
                  ptrdiff_t dst_stride, int8_t log2_pb_w,
                  int8_t log2_pb_h, int angle_val);
void
intra_angular_hor_pdpc(const uint16_t* const ref_abv,
                       const uint16_t* const ref_lft, uint16_t* const dst,
                       ptrdiff_t dst_stride, uint16_t log2_pb_w,
                       uint16_t log2_pb_h);
void
intra_angular_ver_pdpc(const uint16_t* const ref_abv,
                       const uint16_t* const ref_lft, uint16_t* const dst,
                       ptrdiff_t dst_stride, uint16_t log2_pb_w,
                       uint16_t log2_pb_h);
void
intra_angular_hor(const uint16_t* const ref_abv, const uint16_t* const ref_lft,
                  uint16_t* const dst, ptrdiff_t dst_stride, uint16_t log2_pb_w,
                  uint16_t log2_pb_h);
void
intra_angular_ver(const uint16_t* const ref_abv, const uint16_t* const ref_lft,
                  uint16_t* const dst, ptrdiff_t dst_stride, uint16_t log2_pb_w,
                  uint16_t log2_pb_h);
void
intra_angular_h_c_pdpc(const uint16_t* const ref_abv,
                       const uint16_t* const ref_lft, uint16_t* const dst,
                       ptrdiff_t dst_stride, int8_t log2_pb_w,
                       int8_t log2_pb_h, int mode_idx);
void
intra_angular_v_c_pdpc(const uint16_t* const ref_abv,
                       const uint16_t* const ref_lft, uint16_t* const dst,
                       ptrdiff_t dst_stride, int8_t log2_pb_w,
                       int8_t log2_pb_h, int mode_idx);

void
intra_angular_h_nofrac(const uint16_t* ref_lft, uint16_t* dst,
                       ptrdiff_t dst_stride, int8_t log2_pb_w,
                       int8_t log2_pb_h, int angle_val);
void
intra_angular_h_nofrac_pdpc(const uint16_t* ref_abv, const uint16_t* ref_lft,
                            uint16_t* dst, ptrdiff_t dst_stride,
                            int8_t log2_pb_w, int8_t log2_pb_h,
                            int mode_idx);
void
intra_angular_h_gauss_pdpc(const uint16_t* ref_abv, const uint16_t* ref_lft,
                           uint16_t* const dst, ptrdiff_t dst_stride,
                           int8_t log2_pb_w, int8_t log2_pb_h, int mode_idx);
void
intra_angular_v_nofrac(const uint16_t* ref_abv, uint16_t* dst,
                       ptrdiff_t dst_stride, int8_t log2_pb_w,
                       int8_t log2_pb_h, int angle_val);
void
intra_angular_v_nofrac_pdpc(const uint16_t* ref_abv, const uint16_t* ref_lft,
                            uint16_t* const dst, ptrdiff_t dst_stride,
                            int8_t log2_pb_w, int8_t log2_pb_h,
                            int mode_idx);
void
intra_angular_v_gauss_pdpc(const uint16_t* ref_abv, const uint16_t* ref_lft,
                           uint16_t* const dst, ptrdiff_t dst_stride,
                           int8_t log2_pb_w, int8_t log2_pb_h, int mode_idx);
void
intra_angular_h_cubic(const uint16_t* ref_lft, uint16_t* dst,
                      ptrdiff_t dst_stride, int8_t log2_pb_w,
                      int8_t log2_pb_h, int angle_val);
void
intra_angular_h_gauss(const uint16_t* ref_lft, uint16_t* dst,
                      ptrdiff_t dst_stride, int8_t log2_pb_w,
                      int8_t log2_pb_h, int angle_val);
void
intra_angular_v_cubic(const uint16_t* ref_abv, uint16_t* dst,
                      ptrdiff_t dst_stride, int8_t log2_pb_w,
                      int8_t log2_pb_h, int angle_val);
void
intra_angular_v_gauss(const uint16_t* ref_abv, uint16_t* dst,
                      ptrdiff_t dst_stride, int8_t log2_pb_w,
                      int8_t log2_pb_h, int angle_val);
void
intra_angular_h_cubic_pdpc(const uint16_t* ref_abv, const uint16_t* ref_lft,
                           uint16_t* const dst, ptrdiff_t dst_stride,
                           int8_t log2_pb_w, int8_t log2_pb_h, int mode_idx);
void
intra_angular_v_cubic_pdpc(const uint16_t* ref_abv, const uint16_t* ref_lft,
                           uint16_t* const dst, ptrdiff_t dst_stride,
                           int8_t log2_pb_w, int8_t log2_pb_h, int mode_idx);
void
intra_angular_h_cubic_mref(const uint16_t* const ref_lft, uint16_t* const dst,
                           ptrdiff_t dst_stride,
                           int8_t log2_pb_w, int8_t log2_pb_h,
                           int angle_val, uint8_t multi_ref_idx);
void
intra_angular_v_cubic_mref(const uint16_t* const ref_abv, uint16_t* const dst,
                           int dst_stride, int8_t log2_pb_w,
                           int8_t log2_pb_h,
                           int angle_val, uint8_t multi_ref_idx);

#endif // RCN_INTRA_ANGULAR_H
