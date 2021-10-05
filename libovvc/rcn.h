#ifndef RCN_H
#define RCN_H

#include "ovdefs.h"


/* FIXME
 * rename OVINTRA
 * MERGE with INTER modes ?
 * STORE on uint8_t
 * find meaningful place for definition
 */


#if 0
extern const struct TrFunc tr_templates[NB_TR_TYPES][NB_TR_SIZES];
#endif
/* FIXME remove some args and give RCNCTX instead of ctudec
 */
struct OVRCNCtx;
struct InterDRVCtx;
struct RCNFunctions;
struct DBFInfo;
struct OVBuffInfo;
struct VVCGPM;

void rcn_residual(OVCTUDec *const ctudec,
             int16_t *const dst, int16_t *src,
             uint8_t x0, uint8_t y0,
             unsigned int log2_tb_w, unsigned int log2_tb_h,
             unsigned int lim_cg_w,
             uint8_t cu_mts_flag, uint8_t cu_mts_idx,
             uint8_t is_dc, uint8_t lfnst_flag, uint8_t is_mip, uint8_t lfnst_idx, uint8_t sbt);

void rcn_residual_c(OVCTUDec *const ctudec,
                    int16_t *const dst, int16_t *src,
                    uint8_t x0, uint8_t y0,
                    uint8_t log2_tb_w, uint8_t log2_tb_h,
                    uint16_t lim_cg_w,
                    uint8_t lfnst_flag, uint8_t lfnst_idx);

void rcn_frame_line_to_ctu(const struct OVRCNCtx *const rcn_ctx, uint8_t log2_ctb_s);

void rcn_intra_line_to_ctu(const struct OVRCNCtx *const rcn_ctx, int x_l, uint8_t log2_ctb_s);
void rcn_ctu_to_intra_line(OVCTUDec *const ctudec, int x_l);

void ctu_copy_left_border(struct OVRCNCtx *rcn_ctx, uint8_t log2_ctb_s);

void rcn_update_ctu_border(struct OVRCNCtx *rcn_ctx, uint8_t log2_ctb_s);

void rcn_write_ctu_to_frame(const struct OVRCNCtx *const rcn_ctx, uint8_t log2_ctb_s);

void rcn_write_ctu_to_frame_border(const struct OVRCNCtx *const rcn_ctx,
                                   int last_ctu_w, int last_ctu_h);

void rcn_mcp_b(OVCTUDec*const lc_ctx, struct OVBuffInfo dst, struct InterDRVCtx *const inter_ctx,
               const OVPartInfo *const part_ctx,
               const OVMV mv0, const OVMV mv1,
               unsigned int x0, unsigned int y0,
               unsigned int log2_pb_w, unsigned int log2_pb_h,
               uint8_t inter_dir, uint8_t ref_idx0, uint8_t ref_idx1);

void rcn_mcp_b_l(OVCTUDec*const lc_ctx, struct OVBuffInfo dst, struct InterDRVCtx *const inter_ctx,
                 const OVPartInfo *const part_ctx,
                 const OVMV mv0, const OVMV mv1,
                 unsigned int x0, unsigned int y0,
                 unsigned int log2_pb_w, unsigned int log2_pb_h,
                 uint8_t inter_dir, uint8_t ref_idx0, uint8_t ref_idx1);

struct PROFInfo;
void rcn_prof_mcp_b_l(OVCTUDec*const lc_ctx, struct OVBuffInfo dst, struct InterDRVCtx *const inter_ctx,
                      const OVPartInfo *const part_ctx,
                      const OVMV mv0, const OVMV mv1,
                      unsigned int x0, unsigned int y0,
                      unsigned int log2_pb_w, unsigned int log2_pb_h,
                      uint8_t inter_dir, uint8_t ref_idx0, uint8_t ref_idx1,
                      uint8_t prof_dir, const struct PROFInfo *const prof_info);

void rcn_mcp_b_c(OVCTUDec*const lc_ctx, struct OVBuffInfo dst, struct InterDRVCtx *const inter_ctx,
                 const OVPartInfo *const part_ctx,
                 const OVMV mv0, const OVMV mv1,
                 unsigned int x0, unsigned int y0,
                 unsigned int log2_pb_w, unsigned int log2_pb_h,
                 uint8_t inter_dir, uint8_t ref_idx0, uint8_t ref_idx1);

void rcn_bdof_mcp_l(OVCTUDec *const ctudec, struct OVBuffInfo dst,
                    uint8_t x0, uint8_t y0,
                    uint8_t log2_pu_w, uint8_t log2_pu_h,
                    OVMV mv0, OVMV mv1, uint8_t ref_idx0, uint8_t ref_idx1);

void rcn_mcp(OVCTUDec *const ctudec, struct OVBuffInfo dst, int x0, int y0, int log2_pu_w, int log2_pu_h,
             OVMV mv, uint8_t inter_dir, uint8_t ref_idx);

void rcn_ciip_b(OVCTUDec*const ctudec, const OVMV mv0, const OVMV mv1,
           unsigned int x0, unsigned int y0,
           unsigned int log2_pb_w, unsigned int log2_pb_h,
           uint8_t inter_dir, uint8_t ref_idx0, uint8_t ref_idx1);

void rcn_ciip(OVCTUDec *const ctudec,
         int x0, int y0, int log2_pb_w, int log2_pb_h,
         OVMV mv, uint8_t ref_idx);

void rcn_init_gpm_params();

void rcn_gpm_b(OVCTUDec *const ctudec, struct VVCGPM* gpm_ctx, int x0, int y0, int log2_pb_w, int log2_pb_h);

/* FIXME check vertical / horizontal */
void rcn_init_functions(struct RCNFunctions *rcn_func, uint8_t ict_type, uint8_t lm_chroma_enabled,
                        uint8_t sps_chroma_vertical_collocated_flag, uint8_t lmcs_flag);

void rcn_init_cclm_functions(struct RCNFunctions *rcn_func);

void rcn_init_cclm_functions_collocated(struct RCNFunctions *rcn_func);

void rcn_init_ict_functions(struct RCNFunctions *rcn_func, uint8_t type);

void rcn_init_tr_functions(struct RCNFunctions *const rcn_funcs);

void rcn_init_dc_planar_functions(struct RCNFunctions *const rcn_funcs);

void rcn_init_lfnst_functions(struct RCNFunctions *rcn_func);

void rcn_init_mip_functions(struct RCNFunctions *const rcn_funcs);

void rcn_init_sao_functions(struct RCNFunctions *const rcn_funcs);

void rcn_init_dmvr_functions(struct RCNFunctions *const rcn_funcs);

void rcn_init_prof_functions(struct RCNFunctions *const rcn_funcs);

void rcn_init_bdof_functions(struct RCNFunctions *const rcn_funcs);

void rcn_init_ciip_functions(struct RCNFunctions *const rcn_funcs);

void rcn_init_df_functions(struct RCNFunctions *const rcn_funcs);

void rcn_dbf_ctu(const struct OVRCNCtx  *const rcn_ctx, struct DBFInfo *const dbf_info,
                 uint8_t log2_ctu_s, uint8_t last_x, uint8_t last_y);

void rcn_dbf_truncated_ctu(const struct OVRCNCtx  *const rcn_ctx, struct DBFInfo *const dbf_info,
                           uint8_t log2_ctu_s, uint8_t last_x, uint8_t last_y,
                           uint8_t ctu_w, uint8_t ctu_h);

uint8_t rcn_dmvr_mv_refine(OVCTUDec *const ctudec, struct OVBuffInfo dst,
                           uint8_t x0, uint8_t y0,
                           uint8_t log2_pu_w, uint8_t log2_pu_h,
                           OVMV *mv0, OVMV *mv1, uint8_t ref_idx0, uint8_t ref_idx1, uint8_t
                           apply_bdof);
void rcn_lmcs_compute_chroma_scale(struct OVCTUDec* ctudec, int x0, int y0);
void
vvc_intra_pred_mip_tr(const struct OVRCNCtx *const rcn_ctx,
                     uint16_t *const dst,
                     int x0, int y0, int log2_pu_w, int log2_pu_h,
                     uint8_t mip_mode);

void
vvc_intra_pred_mip(const struct OVRCNCtx *const rcn_ctx,
                  uint16_t *const dst,
                  int x0, int y0, int log2_pu_w, int log2_pu_h,
                  uint8_t mip_mode);

void
vvc_add_residual(const int16_t *src, uint16_t *dst,
                 int log2_tb_w, int log2_tb_h,
                 int scale);

void
vvc_sub_residual(const int16_t *src, uint16_t *dst,
                int log2_tb_w, int log2_tb_h,
                int scale);

void
vvc_add_half_residual(const int16_t *src, uint16_t *dst,
                      int log2_tb_w, int log2_tb_h,
                      int scale);

void
vvc_sub_half_residual(const int16_t *src, uint16_t *dst,
                      int log2_tb_w, int log2_tb_h,
                      int scale);

void
vvc_scale_add_residual(const int16_t *src, uint16_t *dst,
                     int log2_tb_w, int log2_tb_h,
                     int scale);

void
vvc_scale_sub_residual(const int16_t *src, uint16_t *dst,
                     int log2_tb_w, int log2_tb_h,
                     int scale);

void
vvc_scale_add_half_residual(const int16_t *src, uint16_t *dst,
                          int log2_tb_w, int log2_tb_h,
                          int scale);

void
vvc_scale_sub_half_residual(const int16_t *src, uint16_t *dst,
                           int log2_tb_w, int log2_tb_h,
                           int scale);

#endif //RCN_H
