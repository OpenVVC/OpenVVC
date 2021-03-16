#ifndef RCN_H
#define RCN_H

#include "ctudec.h"

/* FIXME
 * rename OVINTRA
 * MERGE with INTER modes ?
 * STORE on uint8_t
 * find meaningful place for definition
 */


/* FIXME remove some args and give RCNCTX instead of ctudec
 */

void rcn_residual(OVCTUDec *const ctudec,
             int16_t *const dst, int16_t *src,
             uint8_t x0, uint8_t y0,
             unsigned int log2_tb_w, unsigned int log2_tb_h,
             unsigned int lim_cg_w,
             uint8_t cu_mts_flag, uint8_t cu_mts_idx,
             uint8_t is_dc, uint8_t lfnst_flag, uint8_t is_mip, uint8_t lfnst_idx);

void rcn_residual_c(OVCTUDec *const ctudec,
               int16_t *const dst, int16_t *src,
               int16_t *const lfnst_sb,
               uint8_t x0, uint8_t y0,
               unsigned int log2_tb_w, unsigned int log2_tb_h,
               unsigned int lim_cg_w,
               uint8_t cu_mts_flag, uint8_t cu_mts_idx,
               uint8_t is_dc, uint8_t lfnst_flag, uint8_t is_mip, uint8_t lfnst_idx);

void rcn_frame_line_to_ctu(const struct OVRCNCtx *const rcn_ctx, uint8_t log2_ctb_s);

void ctu_copy_left_border(struct OVRCNCtx *rcn_ctx, uint8_t log2_ctb_s);

void rcn_update_ctu_border(struct OVRCNCtx *rcn_ctx, uint8_t log2_ctb_s);

void rcn_write_ctu_to_frame(const struct OVRCNCtx *const rcn_ctx, uint8_t log2_ctb_s);

void rcn_write_ctu_to_frame_border(const struct OVRCNCtx *const rcn_ctx,
                                   int last_ctu_w, int last_ctu_h);

void rcn_mcp_b(OVCTUDec*const lc_ctx, struct InterDRVCtx *const inter_ctx,
               const OVPartInfo *const part_ctx,
               const OVMV mv0, const OVMV mv1,
               unsigned int x0, unsigned int y0,
               unsigned int log2_pb_w, unsigned int log2_pb_h,
               uint8_t inter_dir);

void rcn_mcp(OVCTUDec *const ctudec, int x0, int y0, int log2_pu_w, int log2_pu_h,
             OVMV mv, uint8_t inter_dir);

void rcn_init_cclm_functions_collocated(struct RCNFunctions *rcn_func);

/* FIXME check vertical / horizontal */
void rcn_init_cclm_functions(struct RCNFunctions *rcn_func);

void rcn_init_ict_functions(struct RCNFunctions *rcn_func, uint8_t type);

void rcn_dbf_ctu(const struct OVRCNCtx  *const rcn_ctx, const struct DBFInfo *const dbf_info,
                 uint8_t log2_ctu_s, uint8_t last_x, uint8_t last_y);

void rcn_dbf_truncated_ctu(const struct OVRCNCtx  *const rcn_ctx, struct DBFInfo *const dbf_info,
                           uint8_t log2_ctu_s, uint8_t last_x, uint8_t last_y,
                           uint8_t ctu_w, uint8_t ctu_h);
#endif //RCN_H
