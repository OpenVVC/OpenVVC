#ifndef DRV_H
#define DRV_H
#include <stdint.h>

#include "ovdefs.h"
#include "ctudec.h"

enum OVIntraMode
{
    OVINTRA_PLANAR = 0,
    OVINTRA_DC = 1,
    OVINTRA_HOR = 18,
    OVINTRA_DIA = 34,
    OVINTRA_VER = 50,
    OVINTRA_VDIA = 66,
    OVINTRA_LM_CHROMA = 67,
    OVINTRA_NUM_INTRA_MODES = 67,
    OVINTRA_MDLM_LEFT = 68,
    OVINTRA_MDLM_TOP  = 69,
    OVINTRA_DM_CHROMA = 70,
    OVINTRA_MIP_MODE = 75,
    OVINTRA_NOT_AVAILABLE=128,
};

struct OVRCNCtx;
struct OVMVCtx;
struct InterDRVCtx;

uint8_t derive_intra_mode_c(uint8_t cclm_flag, uint8_t mpm_flag,
                            uint8_t mpm_idx, uint8_t luma_mode,
                            uint8_t cclm_idx);

/* FIXME vvc_intra functions deal with both reconstruction
 * and derivation this has to be clarified
 */
void vvc_intra_pred(const struct OVRCNCtx *const rcn_ctx, struct OVBuffInfo* ctu_buff,
                    uint8_t intra_mode, int x0, int y0,
                    int log2_pb_width, int log2_pb_height);

void vvc_intra_pred_chroma(const struct OVRCNCtx *const rcn_ctx, struct OVBuffInfo *ctu_buff,
                           uint8_t intra_mode, int x0, int y0,
                           int log2_pb_w, int log2_pb_h);

void vvc_intra_pred_isp(const OVCTUDec *const ctudec,
                        uint16_t *const src,
                        ptrdiff_t dst_stride,
                        uint8_t intra_mode,
                        int x0, int y0,
                        int log2_pb_width, int log2_pb_height,
                        int log2_cu_width,int log2_cu_height,
                        int offset_x, int offset_y);

void vvc_intra_pred_multi_ref( const OVCTUDec *const ctudec,
                               uint16_t *const src,
                               ptrdiff_t dst_stride,
                               uint8_t intra_mode, int x0, int y0,
                               int log2_pb_width, int log2_pb_height,
                               int multi_ref_idx);


uint8_t drv_intra_cu(OVCTUDec *const ctudec, const OVPartInfo *const part_ctx,
                     uint8_t x0, uint8_t y0, uint8_t log2_cb_w, uint8_t log2_cb_h,
                     VVCCU cu);

OVMV drv_mvp_mvd(struct InterDRVCtx *const inter_ctx,
                 const struct OVMVCtx *const mv_ctx,
                 OVMV mvd, uint8_t prec_amvr,
                 uint8_t pb_x, uint8_t pb_y,
                 uint8_t nb_pb_w, uint8_t nb_pb_h,
                 uint8_t mvp_idx, uint8_t inter_dir,
                 uint8_t ref_idx0, uint8_t ref_idx1);


OVMV drv_merge_mvp(struct InterDRVCtx *const inter_ctx,
                   const struct OVMVCtx *const mv_ctx,
                   uint8_t pb_x, uint8_t pb_y,
                   uint8_t nb_pb_w, uint8_t nb_pb_h,
                   uint8_t merge_idx, uint8_t max_nb_merge_cand);

OVMV drv_mmvd_merge_mvp(struct InterDRVCtx *const inter_ctx,
              const struct OVMVCtx *const mv_ctx,
              uint8_t pb_x, uint8_t pb_y,
              uint8_t nb_pb_w, uint8_t nb_pb_h,
              uint8_t merge_idx, uint8_t max_nb_merge_cand);

VVCMergeInfo drv_merge_mvp_b(const struct InterDRVCtx *const inter_ctx,
                             uint8_t pb_x, uint8_t pb_y,
                             uint8_t nb_pb_w, uint8_t nb_pb_h,
                             uint8_t merge_idx, uint8_t max_nb_merge_cand,
                             uint8_t is_small);

VVCMergeInfo drv_mmvd_merge_mvp_b(struct InterDRVCtx *const inter_ctx,
                                  uint8_t x0, uint8_t y0,
                                  uint8_t log2_cu_w, uint8_t log2_cu_h,
                                  uint8_t merge_idx,
                                  uint8_t max_nb_cand, uint8_t is_small);

VVCMergeInfo drv_mvp_b(struct InterDRVCtx *const inter_ctx,
                       uint8_t pb_x, uint8_t pb_y,
                       uint8_t nb_pb_w, uint8_t nb_pb_h,
                       OVMV mvd0, OVMV mvd1, int prec_amvr,
                       uint8_t mvp_idx0, uint8_t mvp_idx1, uint8_t bcw_idx,
                       uint8_t inter_dir, uint8_t ref_idx0, uint8_t ref_idx1,
                       uint8_t is_small);

void drv_affine_mvp_b(struct InterDRVCtx *const inter_ctx,
                      uint8_t x0, uint8_t y0,
                      uint8_t log2_cu_w, uint8_t log2_cu_h,
                      struct AffineControlInfo * cp_mvd0,
                      struct AffineControlInfo * cp_mvd1,
                      uint8_t mvp_idx0, uint8_t mvp_idx1, uint8_t bcw_idx,
                      uint8_t inter_dir, uint8_t ref_idx0, uint8_t ref_idx1,
                      uint8_t affine_type);

void drv_affine_merge_mvp_b(struct InterDRVCtx *const inter_ctx,
                            uint8_t x0, uint8_t y0,
                            uint8_t log2_cu_w, uint8_t log2_cu_h,
                            uint8_t merge_idx);

void drv_affine_mvp_p(struct InterDRVCtx *const inter_ctx,
                      uint8_t x0, uint8_t y0,
                      uint8_t log2_cu_w, uint8_t log2_cu_h,
                      struct AffineControlInfo * cp_mvd0,
                      struct AffineControlInfo * cp_mvd1,
                      uint8_t mvp_idx0, uint8_t mvp_idx1, uint8_t bcw_idx,
                      uint8_t inter_dir, uint8_t ref_idx0, uint8_t ref_idx1,
                      uint8_t affine_type);

void drv_affine_merge_mvp_p(struct InterDRVCtx *const inter_ctx,
                            uint8_t x0, uint8_t y0,
                            uint8_t log2_cu_w, uint8_t log2_cu_h,
                            uint8_t merge_idx);

void drv_gpm_merge_mvp_b(struct InterDRVCtx *const inter_ctx,
                         uint8_t x0, uint8_t y0,
                         uint8_t log2_cu_w, uint8_t log2_cu_h,
                         uint8_t max_nb_cand, uint8_t is_small);

OVMV drv_change_precision_mv(OVMV mv, int src, int dst);

OVMV drv_round_to_precision_mv(OVMV mv, int src, int dst);

void process_lfnst(OVCTUDec *const ctudec,
              int16_t *dst, const int16_t *src,
              int log2_tb_w, int log2_tb_h,
              int x0, int y0, uint8_t lfnst_idx);

void process_lfnst_luma(OVCTUDec *const ctudec,
                   int16_t *dst, const int16_t *src,
                   int log2_tb_w, int log2_tb_h,
                   int x0, int y0, uint8_t lfnst_idx);

void process_lfnst_luma_isp(OVCTUDec *const ctudec,
                            int16_t *dst, const int16_t *src,
                            int log2_tb_w, int log2_tb_h,
                            int log2_cb_w, int log2_cb_h,
                            int x0, int y0, uint8_t lfnst_idx);

void tmvp_inter_synchronization(OVPicture *ref_pic, int ctb_x, int ctb_y, int log2_ctu_s);

#endif
