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

    #if 1
    for (i = 0; i < 1 << (log2_ctb_s - 1); ++i) {
        #if 1
        memcpy(dst_cb, src_cb, sizeof(uint16_t) << (log2_ctb_s - 1));
        memcpy(dst_cr, src_cr, sizeof(uint16_t) << (log2_ctb_s - 1));
        #else
        memset(dst_cb, 1023, sizeof(uint16_t) << (log2_ctb_s - 1));
        memset(dst_cr, 0, sizeof(uint16_t) << (log2_ctb_s - 1));
        #endif
        src_cb += RCN_CTB_STRIDE;
        src_cr += RCN_CTB_STRIDE;
        dst_cb += fd->stride_c;
        dst_cr += fd->stride_c;
    }
    #endif
}

void
rcn_frame_line_to_ctu(const struct OVRCNCtx *const rcn_ctx, uint8_t log2_ctb_s)
{
    int i;
    const struct OVBuffInfo *const fd = &rcn_ctx->frame_buff;
    const uint16_t *src_y  =fd->y  - fd->stride;
    const uint16_t *src_cb =fd->cb - fd->stride_c;
    const uint16_t *src_cr =fd->cr - fd->stride_c;

    uint16_t *dst_y  =  rcn_ctx->ctu_buff.y  - RCN_CTB_STRIDE;
    uint16_t *dst_cb =  rcn_ctx->ctu_buff.cb - RCN_CTB_STRIDE;
    uint16_t *dst_cr =  rcn_ctx->ctu_buff.cr - RCN_CTB_STRIDE;

    memcpy(dst_y,  src_y , (2 << log2_ctb_s) + (1 << log2_ctb_s));
    #if 1
    memcpy(dst_cb, src_cb, (2 << (log2_ctb_s - 1)) + (1 << (log2_ctb_s - 1)));
    memcpy(dst_cr, src_cr, (2 << (log2_ctb_s - 1)) + (1 << (log2_ctb_s - 1)));
    #endif
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
