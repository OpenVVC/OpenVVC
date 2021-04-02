#ifndef DATA_SCAN_LUT_H
#define DATA_SCAN_LUT_H

#include <stdint.h>

extern const uint8_t ff_vvc_diag_scan_2x8_num_cg [16];
extern const uint8_t ff_vvc_diag_scan_4x4_num_cg [16];
extern const uint8_t ff_vvc_diag_scan_8x2_num_cg [16];

extern const uint8_t ff_vvc_inv_diag_scan_2x8 [16];
extern const uint8_t ff_vvc_inv_diag_scan_4x4 [16];
extern const uint8_t ff_vvc_inv_diag_scan_8x2 [16];

extern const uint8_t *const ff_vvc_idx_2_num[6][6];

extern const uint8_t *const ff_vvc_scan_x[6][6];
extern const uint8_t *const ff_vvc_scan_y[6][6];

#endif
