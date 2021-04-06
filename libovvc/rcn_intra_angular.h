#ifndef RCN_INTRA_ANGULAR_H
#define RCN_INTRA_ANGULAR_H

#include <stddef.h>
#include <stdint.h>

void
vvc_intra_angular_hdia(const uint16_t* const ref_above,
                       const uint16_t* const ref_left, uint16_t* const dst,
                       ptrdiff_t dst_stride, int log2_pb_width,
                       int log2_pb_height);
void
vvc_intra_angular_vdia(const uint16_t* const ref_above,
                       const uint16_t* const ref_left, uint16_t* const dst,
                       ptrdiff_t dst_stride, int log2_pb_width,
                       int log2_pb_height);
void
vvc_intra_angular_h_c(const uint16_t* ref_left, uint16_t* dst,
                      ptrdiff_t dst_stride, int log2_pb_width,
                      int log2_pb_height, int angle_val);
void
vvc_intra_angular_v_c(const uint16_t* ref_above, uint16_t* dst,
                      ptrdiff_t dst_stride, int log2_pb_width,
                      int log2_pb_height, int angle_val);
void
vvc_intra_hor_pdpc(const uint16_t* const ref_above,
                   const uint16_t* const ref_left, uint16_t* const dst,
                   ptrdiff_t dst_stride, uint16_t log2_pb_width,
                   uint16_t log2_pb_height);
void
vvc_intra_ver_pdpc(const uint16_t* const ref_above,
                   const uint16_t* const ref_left, uint16_t* const dst,
                   ptrdiff_t dst_stride, uint16_t log2_pb_width,
                   uint16_t log2_pb_height);
void
vvc_intra_hor(const uint16_t* const ref_above, const uint16_t* const ref_left,
              uint16_t* const dst, ptrdiff_t dst_stride, uint16_t log2_pb_width,
              uint16_t log2_pb_height);
void
vvc_intra_ver(const uint16_t* const ref_above, const uint16_t* const ref_left,
              uint16_t* const dst, ptrdiff_t dst_stride, uint16_t log2_pb_width,
              uint16_t log2_pb_height);
void
vvc_intra_angular_hpos_wide(const uint16_t* const ref_above,
                            const uint16_t* const ref_left, uint16_t* const dst,
                            ptrdiff_t dst_stride, int log2_pb_width,
                            int log2_pb_height, int mode_idx);
void
vvc_intra_angular_vpos_wide(const uint16_t* const ref_above,
                            const uint16_t* const ref_left, uint16_t* const dst,
                            ptrdiff_t dst_stride, int log2_pb_width,
                            int log2_pb_height, int mode_idx);

void
intra_angular_h_nofrac(const uint16_t* ref_left, uint16_t* dst,
                       ptrdiff_t dst_stride, int log2_pb_width,
                       int log2_pb_height, int angle_val);
void
intra_angular_h_nofrac_pdpc(const uint16_t* ref_above, const uint16_t* ref_left,
                            uint16_t* dst, ptrdiff_t dst_stride,
                            int log2_pb_width, int log2_pb_height,
                            int mode_idx);
void
intra_angular_h_gauss_pdpc(const uint16_t* ref_above, const uint16_t* ref_left,
                           uint16_t* const dst, ptrdiff_t dst_stride,
                           int log2_pb_width, int log2_pb_height, int mode_idx);
void
intra_angular_v_nofrac(const uint16_t* ref_above, uint16_t* dst,
                       ptrdiff_t dst_stride, int log2_pb_width,
                       int log2_pb_height, int angle_val);
void
intra_angular_v_nofrac_pdpc(const uint16_t* ref_above, const uint16_t* ref_left,
                            uint16_t* const dst, ptrdiff_t dst_stride,
                            int log2_pb_width, int log2_pb_height,
                            int mode_idx);
void
intra_angular_v_gauss_pdpc(const uint16_t* ref_above, const uint16_t* ref_left,
                           uint16_t* const dst, ptrdiff_t dst_stride,
                           int log2_pb_width, int log2_pb_height, int mode_idx);
void
intra_angular_h_cubic(const uint16_t* ref_left, uint16_t* dst,
                      ptrdiff_t dst_stride, int log2_pb_width,
                      int log2_pb_height, int angle_val);
void
intra_angular_h_gauss(const uint16_t* ref_left, uint16_t* dst,
                      ptrdiff_t dst_stride, int log2_pb_width,
                      int log2_pb_height, int angle_val);
void
intra_angular_v_cubic(const uint16_t* ref_above, uint16_t* dst,
                      ptrdiff_t dst_stride, int log2_pb_width,
                      int log2_pb_height, int angle_val);
void
intra_angular_v_gauss(const uint16_t* ref_above, uint16_t* dst,
                      ptrdiff_t dst_stride, int log2_pb_width,
                      int log2_pb_height, int angle_val);
void
intra_angular_h_cubic_pdpc(const uint16_t* ref_above, const uint16_t* ref_left,
                           uint16_t* const dst, ptrdiff_t dst_stride,
                           int log2_pb_width, int log2_pb_height, int mode_idx);
void
intra_angular_v_cubic_pdpc(const uint16_t* ref_above, const uint16_t* ref_left,
                           uint16_t* const dst, ptrdiff_t dst_stride,
                           int log2_pb_width, int log2_pb_height, int mode_idx);
void
intra_angular_hdia_mref(const uint16_t* const ref_above,
                        const uint16_t* const ref_left, uint16_t* const dst,
                        ptrdiff_t dst_stride, int log2_pb_width,
                        int log2_pb_height, uint8_t multi_ref_idx);
void
intra_angular_vdia_mref(const uint16_t* const ref_above,
                        const uint16_t* const ref_left, uint16_t* const dst,
                        ptrdiff_t dst_stride, int log2_pb_width,
                        int log2_pb_height, uint8_t multi_ref_idx);
void
intra_angular_h_cubic_mref(const uint16_t* const ref_left, uint16_t* const dst,
                           ptrdiff_t dst_stride,
                           int log2_pb_width, int log2_pb_height,
                           int angle_val, uint8_t multi_ref_idx);
void
intra_angular_v_cubic_mref(const uint16_t* const ref_above, uint16_t* const dst,
                           int dst_stride, int log2_pb_width,
                           int log2_pb_height,
                           int angle_val, uint8_t multi_ref_idx);

#endif // RCN_INTRA_ANGULAR_H
