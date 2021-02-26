#ifndef DRV_H
#define DRV_H
#include <stdint.h>
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


uint8_t derive_intra_mode_c(uint8_t cclm_flag, uint8_t mpm_flag,
                            uint8_t mpm_idx, uint8_t luma_mode,
                            uint8_t cclm_idx);

/* FIXME vvc_intra functions deal with both reconstruction
 * and derivation this has to be clarified
 */
void vvc_intra_pred(const struct OVRCNCtx *const rcn_ctx,
                    uint8_t intra_mode, int x0, int y0,
                    int log2_pb_width, int log2_pb_height);

void vvc_intra_pred_chroma(const struct OVRCNCtx *const rcn_ctx,
                           uint8_t intra_mode, int x0, int y0,
                           int log2_pb_w, int log2_pb_h);
VVCCU
drv_intra_cu(OVCTUDec *const ctudec, const OVPartInfo *const part_ctx,
             uint8_t x0, uint8_t y0, uint8_t log2_cb_w, uint8_t log2_cb_h,
             VVCCU cu);

VVCMergeInfo derive_mvp_b(struct InterDRVCtx *const inter_ctx,
                          const OVPartInfo *const part_ctx,
                          unsigned int x0, unsigned int y0,
                          unsigned int log2_pb_w, unsigned int log2_pb_h,
                          OVMV mvd0, OVMV mvd1,
                          uint8_t mvp_idx0, uint8_t mvp_idx1,
                          uint8_t inter_dir);

OVMV
derive_mvp_candidates(struct InterDRVCtx *const inter_ctx,
                      const struct OVMVCtx *const mv_ctx,
                      uint8_t pb_x, uint8_t pb_y,
                      uint8_t n_pb_w, uint8_t n_pb_h,
                      uint8_t mvp_idx, uint8_t inter_dir);
OVMV derive_mvp_mvd(struct InterDRVCtx *const inter_ctx,
                    const struct OVMVCtx *const mv_ctx,
                    OVMV mvd,
                    uint8_t pb_x, uint8_t pb_y,
                    uint8_t n_pb_w, uint8_t n_pb_h,
                    uint8_t mvp_idx, uint8_t inter_dir);

OVMV vvc_derive_merge_mvp(const struct InterDRVCtx *const inter_ctx,
                          const struct OVMVCtx *const mv_ctx,
                          uint8_t pb_x, uint8_t pb_y,
                          uint8_t n_pb_w, uint8_t n_pb_h,
                          uint8_t merge_idx, uint8_t max_nb_merge_cand);

VVCMergeInfo vvc_derive_merge_mvp_b(const struct InterDRVCtx *const inter_ctx,
                                    uint8_t pb_x, uint8_t pb_y,
                                    uint8_t n_pb_w, uint8_t n_pb_h,
                                    uint8_t merge_idx, uint8_t max_nb_merge_cand);
#endif
