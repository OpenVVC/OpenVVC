#include <stdint.h>
#include <string.h>

#include "ctudec.h"

void
rcn_ctu_copy_left_border(struct OVRCNCtx *rcn_ctx, uint8_t log2_ctb_s)
{
     uint8_t ctu_s = 1 << log2_ctb_s;
     struct OVBuffInfo *binfo = &rcn_ctx->ctu_buff;
     const uint16_t *src = binfo->y + ctu_s - 4;
     uint16_t *dst = binfo->y - 4;
     int i;

     /* Note we also copy border of above line for top left sample
      */
     src -= binfo->stride_c;
     dst -= binfo->stride_c;

     /* Copy 4 left samples because of MRL */
     for (i = 0; i < ctu_s + 1; ++i) {
         memcpy(dst, src, sizeof(*dst) * 4);
         dst += binfo->stride;
         src += binfo->stride;
     }

     src = binfo->cb + (ctu_s >> 1) - 1;
     dst = binfo->cb - 1;
     src -= binfo->stride_c;
     dst -= binfo->stride_c;

     for (i = 0; i < (ctu_s >> 1) + 1; ++i) {
         memcpy(dst, src, sizeof(*dst));
         dst += binfo->stride_c;
         src += binfo->stride_c;
     }

     src = binfo->cr + (ctu_s >> 1) - 1;
     dst = binfo->cr - 1;
     src -= binfo->stride_c;
     dst -= binfo->stride_c;

     for (i = 0; i < (ctu_s >> 1) + 1; ++i) {
         memcpy(dst, src, sizeof(*dst));
         dst += binfo->stride_c;
         src += binfo->stride_c;
     }
}

void
rcn_update_ctu_border(struct OVRCNCtx *rcn_ctx, uint8_t log2_ctb_s)
{
     uint8_t ctu_s = 1 << log2_ctb_s;
     struct OVBuffInfo *binfo = &rcn_ctx->ctu_buff;
     const uint16_t *src = binfo->y + ctu_s - 4;
     uint16_t *dst = binfo->y - 4;
     int i;

     struct OVBuffInfo *finfo = &rcn_ctx->frame_buff;
     finfo->y += ctu_s;
     finfo->cb += ctu_s >> 1;
     finfo->cr += ctu_s >> 1;
     /* Note we also copy border of above line for top left sample
      */
     src -= binfo->stride_c;
     dst -= binfo->stride_c;

     /* Copy 4 left samples because of MRL */
     for (i = 0; i < ctu_s + 1; ++i) {
         memcpy(dst, src, sizeof(*dst) * 4);
         dst += binfo->stride;
         src += binfo->stride;
     }

     src = binfo->cb + (ctu_s >> 1) - 1;
     dst = binfo->cb - 1;
     src -= binfo->stride_c;
     dst -= binfo->stride_c;

     for (i = 0; i < (ctu_s >> 1) + 1; ++i) {
         memcpy(dst, src, sizeof(*dst));
         dst += binfo->stride_c;
         src += binfo->stride_c;
     }

     src = binfo->cr + (ctu_s >> 1) - 1;
     dst = binfo->cr - 1;
     src -= binfo->stride_c;
     dst -= binfo->stride_c;

     for (i = 0; i < (ctu_s >> 1) + 1; ++i) {
         memcpy(dst, src, sizeof(*dst));
         dst += binfo->stride_c;
         src += binfo->stride_c;
     }
}

void
rcn_write_ctu_to_frame(const struct OVRCNCtx *const rcn_ctx, uint8_t log2_ctb_s)
{
    int i;
    const struct OVBuffInfo *const fd = &rcn_ctx->frame_buff;
    const uint16_t *src_y  = &rcn_ctx->ctu_buff.y [0];
    const uint16_t *src_cb = &rcn_ctx->ctu_buff.cb[0];
    const uint16_t *src_cr = &rcn_ctx->ctu_buff.cr[0];

    uint16_t *dst_y  = fd->y;
    uint16_t *dst_cb = fd->cb;
    uint16_t *dst_cr = fd->cr;

    for (i = 0; i < 1 << log2_ctb_s; i++) {
        memcpy(dst_y, src_y, sizeof(uint16_t) << log2_ctb_s);
        src_y += RCN_CTB_STRIDE;
        dst_y += fd->stride;
    }

    for (i = 0; i < 1 << (log2_ctb_s - 1); ++i) {
        memcpy(dst_cb, src_cb, sizeof(uint16_t) << (log2_ctb_s - 1));
        memcpy(dst_cr, src_cr, sizeof(uint16_t) << (log2_ctb_s - 1));
        src_cb += RCN_CTB_STRIDE;
        src_cr += RCN_CTB_STRIDE;
        dst_cb += fd->stride_c;
        dst_cr += fd->stride_c;
    }
}

void
rcn_frame_line_to_ctu(const struct OVRCNCtx *const rcn_ctx, uint8_t log2_ctb_s)
{
    const struct OVBuffInfo *const fd = &rcn_ctx->frame_buff;
    const uint16_t *src_y  =fd->y  - fd->stride;
    const uint16_t *src_cb =fd->cb - fd->stride_c;
    const uint16_t *src_cr =fd->cr - fd->stride_c;

    uint16_t *dst_y  =  rcn_ctx->ctu_buff.y  - RCN_CTB_STRIDE;
    uint16_t *dst_cb =  rcn_ctx->ctu_buff.cb - RCN_CTB_STRIDE;
    uint16_t *dst_cr =  rcn_ctx->ctu_buff.cr - RCN_CTB_STRIDE;

    memcpy(dst_y,  src_y , sizeof(uint16_t) * OVMIN(((1 << log2_ctb_s) + (1 << log2_ctb_s)), RCN_CTB_STRIDE - 16));
    memcpy(dst_cb, src_cb, sizeof(uint16_t) * ((1 << (log2_ctb_s - 1)) + (1 << (log2_ctb_s - 1))));
    memcpy(dst_cr, src_cr, sizeof(uint16_t) * ((1 << (log2_ctb_s - 1)) + (1 << (log2_ctb_s - 1))));
}

void
rcn_intra_line_to_ctu(const struct OVRCNCtx *const rcn_ctx, int x_l, uint8_t log2_ctb_s)
{
    const struct OVBuffInfo *const il = &rcn_ctx->intra_line_buff;
    const uint16_t *src_y  = il->y  + x_l;
    const uint16_t *src_cb = il->cb + (x_l>>1);
    const uint16_t *src_cr = il->cr + (x_l>>1);

    uint16_t *dst_y  =  rcn_ctx->ctu_buff.y  - RCN_CTB_STRIDE;
    uint16_t *dst_cb =  rcn_ctx->ctu_buff.cb - RCN_CTB_STRIDE;
    uint16_t *dst_cr =  rcn_ctx->ctu_buff.cr - RCN_CTB_STRIDE;

    memcpy(dst_y,  src_y , sizeof(uint16_t) * OVMIN(((1 << log2_ctb_s) + (1 << log2_ctb_s)), RCN_CTB_STRIDE - 16));
    memcpy(dst_cb, src_cb, sizeof(uint16_t) * ((1 << (log2_ctb_s - 1)) + (1 << (log2_ctb_s - 1))));
    memcpy(dst_cr, src_cr, sizeof(uint16_t) * ((1 << (log2_ctb_s - 1)) + (1 << (log2_ctb_s - 1))));
}

void rcn_ctu_to_intra_line(OVCTUDec *const ctudec, int x_l)
{
    struct OVRCNCtx *rcn_ctx = &ctudec->rcn_ctx;
    struct OVBuffInfo *intra_line_binfo = &rcn_ctx->intra_line_buff;

    const OVPartInfo *const pinfo = ctudec->part_ctx;
    uint8_t log2_ctb_size = pinfo->log2_ctu_s;
    int max_cu_width_l = 1 << log2_ctb_size;
    int max_cu_width_c = 1 << (log2_ctb_size-1);
    uint16_t *_src;

    _src = &ctudec->rcn_ctx.ctu_buff.y[(max_cu_width_l - 1) * RCN_CTB_STRIDE];
    memcpy(&intra_line_binfo->y[x_l], _src, sizeof(int16_t) << log2_ctb_size);

    _src = &ctudec->rcn_ctx.ctu_buff.cb[(max_cu_width_c - 1) * RCN_CTB_STRIDE ];
    memcpy(&intra_line_binfo->cb[x_l>>1], _src, sizeof(int16_t) << (log2_ctb_size-1));

    _src = &ctudec->rcn_ctx.ctu_buff.cr[(max_cu_width_c - 1) * RCN_CTB_STRIDE ];
    memcpy(&intra_line_binfo->cr[x_l>>1], _src, sizeof(int16_t) << (log2_ctb_size-1));
}

void
rcn_write_ctu_to_frame_border(const struct OVRCNCtx *const rcn_ctx,
                              int last_ctu_w, int last_ctu_h)
{
    const struct OVBuffInfo *const fd = &rcn_ctx->frame_buff;
    const uint16_t *src_y  = &rcn_ctx->ctu_buff.y [0];
    const uint16_t *src_cb = &rcn_ctx->ctu_buff.cb[0];
    const uint16_t *src_cr = &rcn_ctx->ctu_buff.cr[0];

    uint16_t *dst_y  = fd->y;
    uint16_t *dst_cb = fd->cb;
    uint16_t *dst_cr = fd->cr;

    for (int i = 0; i < last_ctu_h; ++i) {
        memcpy(dst_y, src_y, sizeof(uint16_t) * last_ctu_w);
        dst_y += fd->stride;
        src_y += RCN_CTB_STRIDE;
    }

    for (int i = 0; i < (last_ctu_h >> 1); i++) {
        memcpy(dst_cb, src_cb,  sizeof(uint16_t) * (last_ctu_w >> 1));
        memcpy(dst_cr, src_cr,  sizeof(uint16_t) * (last_ctu_w >> 1));
        dst_cb += fd->stride_c;
        dst_cr += fd->stride_c;
        src_cb += RCN_CTB_STRIDE;
        src_cr += RCN_CTB_STRIDE;
    } 
}
