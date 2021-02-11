#ifndef DATA_SCAN_LUT_H
#define DATA_SCAN_LUT_H

#include <stdint.h>
extern const uint8_t ff_vvc_default_num_cg[16];
extern const uint8_t ff_vvc_diag_scan_2x2_num_cg [4];
extern const uint8_t ff_vvc_diag_scan_4x2_num_cg [8];
extern const uint8_t ff_vvc_diag_scan_2x4_num_cg [8];
extern const uint8_t ff_vvc_diag_scan_2x8_num_cg [16];
extern const uint8_t ff_vvc_diag_scan_4x4_num_cg [16];
extern const uint8_t ff_vvc_diag_scan_8x2_num_cg [16];
extern const uint8_t ff_vvc_diag_scan_4x8_num_cg [32];
extern const uint8_t ff_vvc_diag_scan_8x4_num_cg [32];
extern const uint8_t ff_vvc_diag_scan_8x8_num_cg [64];

extern const uint8_t *const ff_vvc_idx_2_num[6][6];

extern const uint8_t ff_vvc_default_num_cg[16];
extern const uint8_t ff_vvc_inv_diag_scan_2x2 [4];
extern const uint8_t ff_vvc_inv_diag_scan_4x2 [8];
extern const uint8_t ff_vvc_inv_diag_scan_2x4 [8];
extern const uint8_t ff_vvc_inv_diag_scan_2x8 [16];
extern const uint8_t ff_vvc_inv_diag_scan_4x4 [16];
extern const uint8_t ff_vvc_inv_diag_scan_8x2 [16];
extern const uint8_t ff_vvc_inv_diag_scan_4x8 [32];
extern const uint8_t ff_vvc_inv_diag_scan_8x4 [32];
extern const uint8_t ff_vvc_inv_diag_scan_8x8 [64];

extern const uint8_t ff_vvc_inv_diag_scan_2x2_x [4];
extern const uint8_t ff_vvc_inv_diag_scan_4x2_x [8];
extern const uint8_t ff_vvc_inv_diag_scan_2x4_x [8];
extern const uint8_t ff_vvc_inv_diag_scan_2x8_x [16];
extern const uint8_t ff_vvc_inv_diag_scan_4x4_x [16];
extern const uint8_t ff_vvc_inv_diag_scan_8x2_x [16];
extern const uint8_t ff_vvc_inv_diag_scan_4x8_x [32];
extern const uint8_t ff_vvc_inv_diag_scan_8x4_x [32];
extern const uint8_t ff_vvc_inv_diag_scan_8x8_x [64];

extern const uint8_t ff_vvc_inv_diag_scan_2x2_y [4];
extern const uint8_t ff_vvc_inv_diag_scan_4x2_y [8];
extern const uint8_t ff_vvc_inv_diag_scan_2x4_y [8];
extern const uint8_t ff_vvc_inv_diag_scan_2x8_y [16];
extern const uint8_t ff_vvc_inv_diag_scan_4x4_y [16];
extern const uint8_t ff_vvc_inv_diag_scan_8x2_y [16];
extern const uint8_t ff_vvc_inv_diag_scan_4x8_y [32];
extern const uint8_t ff_vvc_inv_diag_scan_8x4_y [32];
extern const uint8_t ff_vvc_inv_diag_scan_8x8_y [64];

extern const uint8_t *const ff_vvc_scan_2_idx[6][6];

extern const uint8_t ff_vvc_diag_scan4x4_x[16];
extern const uint8_t ff_vvc_diag_scan4x4_y[16];

extern const uint8_t ff_vvc_per_cg_dist_significant_flag_ctx_offset_0 [16];
extern const uint8_t ff_vvc_per_cg_dist_significant_flag_ctx_offset_1 [16];
extern const uint8_t ff_vvc_per_cg_dist_parity_flag_ctx_offset_0 [16];
extern const uint8_t ff_vvc_per_cg_dist_parity_flag_ctx_offset_1 [16];
extern const uint8_t ff_vvc_per_cg_dist_parity_flag_ctx_offset_2 [16];

extern const uint8_t ff_vvc_per_cg_dist_significant_flag_ctx_offset_isp_1x16 [16];
extern const uint8_t ff_vvc_per_cg_dist_significant_flag_ctx_offset_isp_2x8  [16];
extern const uint8_t ff_vvc_per_cg_dist_parity_flag_ctx_offset_isp_1x16 [16];
extern const uint8_t ff_vvc_per_cg_dist_parity_flag_ctx_offset_isp_2x8_0 [16];
extern const uint8_t ff_vvc_per_cg_dist_parity_flag_ctx_offset_isp_2x8_1 [16];

extern const uint8_t ff_vvc_parity_flag_ctx_offset_chroma [16];
extern const uint8_t ff_vvc_significant_flag_ctx_offset_chroma [16];

extern const uint8_t *const ff_vvc_parity_offset[3];
extern const uint8_t *const ff_vvc_significant_offset[3];

extern const uint8_t *const ff_vvc_scan_x[6][6];
extern const uint8_t *const ff_vvc_scan_y[6][6];

#endif
