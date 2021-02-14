#ifndef RCN_H
#define RCN_H

#include "ctudec.h"
#include "rcn_struct.h"

/* FIXME
 * rename OVINTRA
 * MERGE with INTER modes ?
 * STORE on uint8_t 
 * find meaningful place for definition
 */

// FIXED? :
#define VVC_CTB_STRIDE (128 + 16 + 64)
#define VVC_CTB_OFFSET (VVC_CTB_STRIDE * 4 + 16)

#define VVC_CTB_STRIDE_CHROMA (128 + 64 + 16)
#define VVC_CTB_CTX_SIZE   (VVC_CTB_STRIDE * (128 + 4))

#define VVC_CTU_UP_FLAG             (1 << 0)
#define VVC_CTU_LEFT_FLAG           (1 << 1)
#define VVC_CTU_UPLEFT_FLAG         (1 << 2)
#define VVC_CTU_UPRIGHT_FLAG        (1 << 3)

extern const struct TrFunc tr_templates[NB_TR_TYPES][NB_TR_SIZES];
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

void ctu_copy_left_border(struct OVRCNCtx *rcn_ctx, uint8_t log2_ctb_s);

void rcn_update_ctu_border(struct OVRCNCtx *rcn_ctx, uint8_t log2_ctb_s);

#endif //RCN_H
