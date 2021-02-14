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

void vvc_intra_pred(const OVCTUDec *const ctudec,
                    uint16_t *const src,
                    ptrdiff_t dst_stride,
                    uint8_t intra_mode, int x0, int y0,
                    int log2_pb_width, int log2_pb_height);

void vvc_intra_pred_chroma(const OVCTUDec *const ctudec,
                           uint16_t *const dst_NUU, uint16_t *const dst_NUk,
                           ptrdiff_t dst_stride,
                           uint8_t intra_mode, int x0, int y0,
                           int log2_pb_w, int log2_pb_h);
VVCCU
drv_intra_cu(OVCTUDec *const ctudec, const OVPartInfo *const part_ctx,
             uint8_t x0, uint8_t y0, uint8_t log2_cb_w, uint8_t log2_cb_h,
             VVCCU cu);
#endif
