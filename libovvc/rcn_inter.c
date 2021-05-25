#include <string.h>
#include "ovdpb.h"
#include "dec_structures.h"
#include "ovdefs.h"
#include "ovutils.h"
#include "ctudec.h"
#include "rcn_structures.h"
#include "rcn_lmcs.h"
#include "drv.h"

#define MAX_PB_SIZE 128

#define REF_PADDING_C 1
#define EPEL_EXTRA_AFTER  2
#define EPEL_EXTRA REF_PADDING_C + EPEL_EXTRA_AFTER

#define REF_PADDING_L 3
#define QPEL_EXTRA_AFTER  4
#define QPEL_EXTRA REF_PADDING_L + QPEL_EXTRA_AFTER

static OVMV
clip_mv(int pos_x, int pos_y, int pic_w, int pic_h, int pb_w, int pb_h, OVMV mv)
{
    int x_max  = (pic_w + 2 - pos_x) << 4;
    int y_max  = (pic_h + 2 - pos_y) << 4;
    int x_min  = (-pb_w - 3 - pos_x) << 4;
    int y_min  = (-pb_h - 3 - pos_y) << 4;

    mv.x = ov_clip(mv.x, x_min, x_max);
    mv.y = ov_clip(mv.y, y_min, y_max);

    return mv;
}

static void
emulate_block_border(uint16_t *buf, const uint16_t *src,
                     ptrdiff_t buf_linesize,
                     ptrdiff_t src_linesize,
                     int block_w, int block_h,
                     int src_x, int src_y, int w, int h)
{
    int x, y;
    int start_y, start_x, end_y, end_x;

    if (!w || !h)
        return;

    if (src_y >= h) {
        src -= src_y * src_linesize;
        src += (h - 1) * src_linesize;
        src_y = h - 1;
    } else if (src_y <= -block_h) {
        src -= src_y * src_linesize;
        src += (1 - block_h) * src_linesize;
        src_y = 1 - block_h;
    }

    if (src_x >= w) {
        src  += (w - 1 - src_x);
        src_x = w - 1;
    } else if (src_x <= -block_w) {
        src  += (1 - block_w - src_x);
        src_x = 1 - block_w;
    }

    start_y = OVMAX(0, -src_y);
    start_x = OVMAX(0, -src_x);

    end_y = OVMIN(block_h, h-src_y);
    end_x = OVMIN(block_w, w-src_x);

    w    = end_x - start_x;
    src += start_y * src_linesize + start_x;
    buf += start_x;

    // top
    for (y = 0; y < start_y; y++) {
        memcpy(buf, src, w * sizeof(uint16_t));
        buf += buf_linesize;
    }

    // copy existing part
    for (; y < end_y; y++) {
        memcpy(buf, src, w * sizeof(uint16_t));
        src += src_linesize;
        buf += buf_linesize;
    }

    // bottom
    src -= src_linesize;
    for (; y < block_h; y++) {
        memcpy(buf, src, w * sizeof(uint16_t));
        buf += buf_linesize;
    }

    buf -= block_h * buf_linesize + start_x;

    while (block_h--) {
        uint16_t *bufp = (uint16_t *) buf;

        // left
        for(x = 0; x < start_x; x++) {
            bufp[x] = bufp[start_x];
        }

        // right
        for (x = end_x; x < block_w; x++) {
            bufp[x] = bufp[end_x - 1];
        }
        buf += buf_linesize;
    }
}

static uint8_t
test_for_edge_emulation_c(int pb_x, int pb_y, int pic_w, int pic_h,
                          int pb_w, int pb_h)
{
    uint8_t emulate_edge = 0;
    emulate_edge  =      pb_x - REF_PADDING_C < 0;
    emulate_edge |= 2 * (pb_y - REF_PADDING_C < 0);
    emulate_edge |= 4 * (pb_x >= pic_w);
    emulate_edge |= 8 * (pb_y >= pic_h);
    emulate_edge |= 4 * ((pb_x + pb_w + EPEL_EXTRA_AFTER) >= pic_w);
    emulate_edge |= 8 * ((pb_y + pb_h + EPEL_EXTRA_AFTER) >= pic_h);
    return emulate_edge;
}

static uint8_t
test_for_edge_emulation(int pb_x, int pb_y, int pic_w, int pic_h,
                        int pu_w, int pu_h)
{
    /* FIXME thi could be simplified */
    uint8_t emulate_edge = 0;
    emulate_edge =       pb_x - REF_PADDING_L < 0;
    emulate_edge |= 2 * (pb_y - REF_PADDING_L < 0);

    emulate_edge |= 4 * (pb_x >= pic_w);
    emulate_edge |= 8 * (pb_y >= pic_h);

    emulate_edge |= 4 * ((pb_x + pu_w + QPEL_EXTRA_AFTER) >= pic_w);
    emulate_edge |= 8 * ((pb_y + pu_h + QPEL_EXTRA_AFTER) >= pic_h);
    return emulate_edge;
}

static struct OVBuffInfo
derive_ref_buf_c(const OVPicture *const ref_pic, OVMV mv, int pos_x, int pos_y,
                 uint16_t *edge_buff0, uint16_t *edge_buff1,
                 int log2_pu_w, int log2_pu_h, int log2_ctu_s)
{
    struct OVBuffInfo ref_buff;
    uint16_t *const ref_cb  = (uint16_t *) ref_pic->frame->data[1];
    uint16_t *const ref_cr  = (uint16_t *) ref_pic->frame->data[2];

    int src_stride = ref_pic->frame->linesize[1] >> 1;

    /*FIXME check buff side derivation */
    int ref_pos_x = pos_x + (mv.x >> 5);
    int ref_pos_y = pos_y + (mv.y >> 5);

    const int pu_w = (1 << log2_pu_w) >> 1;
    const int pu_h = (1 << log2_pu_h) >> 1;

    const int pic_w = ref_pic->frame->width[0]  >> 1;
    const int pic_h = ref_pic->frame->height[0] >> 1;

    uint16_t *src_cb  = &ref_cb[ref_pos_x + ref_pos_y * src_stride];
    uint16_t *src_cr  = &ref_cr[ref_pos_x + ref_pos_y * src_stride];

    uint8_t emulate_edge = test_for_edge_emulation_c(ref_pos_x, ref_pos_y, pic_w, pic_h,
                                                     pu_w, pu_h);;

    if (emulate_edge){
        int src_off  = REF_PADDING_C * (src_stride) + (REF_PADDING_C);
        int buff_off = REF_PADDING_C * (RCN_CTB_STRIDE) + (REF_PADDING_C);
        int cpy_w = pu_w + EPEL_EXTRA;
        int cpy_h = pu_h + EPEL_EXTRA;
        /*FIXME clip to frame?*/
        int start_pos_x = ref_pos_x - REF_PADDING_C;
        int start_pos_y = ref_pos_y - REF_PADDING_C;

        emulate_block_border(edge_buff0, (src_cb - src_off),
                             RCN_CTB_STRIDE, src_stride,
                             cpy_w, cpy_h, start_pos_x, start_pos_y,
                             pic_w, pic_h);
        emulate_block_border(edge_buff1, (src_cr - src_off),
                             RCN_CTB_STRIDE, src_stride,
                             cpy_w, cpy_h, start_pos_x, start_pos_y,
                             pic_w, pic_h);

        ref_buff.cb = edge_buff0 + buff_off;
        ref_buff.cr = edge_buff1 + buff_off;
        ref_buff.stride_c = RCN_CTB_STRIDE;

    } else {
        ref_buff.cb = src_cb;
        ref_buff.cr = src_cr;
        ref_buff.stride_c = src_stride;
    }
    return ref_buff;
}

static struct OVBuffInfo
derive_ref_buf_y(const OVPicture *const ref_pic, OVMV mv, int pos_x, int pos_y,
                uint16_t *edge_buff, int log2_pu_w, int log2_pu_h, int log2_ctu_s)
{
    struct OVBuffInfo ref_buff;
    uint16_t *const ref_y  = (uint16_t *) ref_pic->frame->data[0];

    int src_stride = ref_pic->frame->linesize[0] >> 1;

    int ref_pos_x = pos_x + (mv.x >> 4);
    int ref_pos_y = pos_y + (mv.y >> 4);

    const int pu_w = 1 << log2_pu_w;
    const int pu_h = 1 << log2_pu_h;

    const int pic_w = ref_pic->frame->width[0];
    const int pic_h = ref_pic->frame->height[0];

    uint8_t emulate_edge = test_for_edge_emulation(ref_pos_x, ref_pos_y, pic_w, pic_h,
                                                   pu_w, pu_h);;

    /* FIXME Frame thread synchronization here to ensure data is available
     */

    if (emulate_edge){
        const uint16_t *src_y  = &ref_y[ref_pos_x + ref_pos_y * src_stride];
        int src_off  = REF_PADDING_L * (src_stride) + (REF_PADDING_L);
        int buff_off = REF_PADDING_L * (RCN_CTB_STRIDE) + (REF_PADDING_L);
        int cpy_w = pu_w + QPEL_EXTRA;
        int cpy_h = pu_h + QPEL_EXTRA;

        /*FIXME clip to frame?*/

        int start_pos_x = ref_pos_x - REF_PADDING_L;
        int start_pos_y = ref_pos_y - REF_PADDING_L;

        emulate_block_border(edge_buff, (src_y - src_off),
                             RCN_CTB_STRIDE, src_stride,
                             cpy_w, cpy_h, start_pos_x, start_pos_y,
                             pic_w, pic_h);
        ref_buff.y  = edge_buff + buff_off;
        ref_buff.stride = RCN_CTB_STRIDE;

    } else {

        ref_buff.y  = &ref_y[ref_pos_x + ref_pos_y * src_stride];
        ref_buff.stride = src_stride;
    }
    return ref_buff;
}

static void
rcn_motion_compensation_b(OVCTUDec *const ctudec, struct OVBuffInfo dst,
                          uint8_t x0, uint8_t y0,
                          uint8_t log2_pu_w, uint8_t log2_pu_h,
                          OVMV mv0, OVMV mv1, uint8_t ref_idx0, uint8_t ref_idx1)
{
    struct OVRCNCtx    *const rcn_ctx   = &ctudec->rcn_ctx;
    const struct InterDRVCtx *const inter_ctx = &ctudec->drv_ctx.inter_ctx;
    struct MCFunctions *mc_l = &rcn_ctx->rcn_funcs.mc_l;
    struct MCFunctions *mc_c = &rcn_ctx->rcn_funcs.mc_c;
    /* FIXME derive ref_idx */
    uint8_t ref_idx_0 = ref_idx0;
    uint8_t ref_idx_1 = ref_idx1;

    OVPicture *ref0 = inter_ctx->rpl0[ref_idx_0];
    OVPicture *ref1 = inter_ctx->rpl1[ref_idx_1];

    /* TMP buffers for edge emulation
     * FIXME use tmp buffers in local contexts
     */
    uint16_t edge_buff0[RCN_CTB_SIZE];
    uint16_t edge_buff1[RCN_CTB_SIZE];
    uint16_t edge_buff0_1[RCN_CTB_SIZE];
    uint16_t edge_buff1_1[RCN_CTB_SIZE];
    int16_t tmp_buff[RCN_CTB_SIZE];

    /*FIXME we suppose here both refs possess the same size*/

    const int log2_ctb_s = ctudec->part_ctx->log2_ctu_s;

    /* FIXME we should not need ctb_x/y
     * it could be retrieved from position in frame buff
     */
    int pos_x = (ctudec->ctb_x << log2_ctb_s) + x0;
    int pos_y = (ctudec->ctb_y << log2_ctb_s) + y0;

    mv0 = clip_mv(pos_x, pos_y, ref0->frame->width[0],
                  ref0->frame->height[0], 1 << log2_pu_w, 1 << log2_pu_h, mv0);

    mv1 = clip_mv(pos_x, pos_y, ref1->frame->width[0],
                  ref1->frame->height[0], 1 << log2_pu_w, 1 << log2_pu_h, mv1);


    const struct OVBuffInfo ref0_b = derive_ref_buf_y(ref0, mv0, pos_x, pos_y, edge_buff0,
                                                      log2_pu_w, log2_pu_h, log2_ctb_s);

    const struct OVBuffInfo ref1_b = derive_ref_buf_y(ref1, mv1, pos_x, pos_y, edge_buff1,
                                                      log2_pu_w, log2_pu_h, log2_ctb_s);

    const int pu_w = 1 << log2_pu_w;
    const int pu_h = 1 << log2_pu_h;

    uint8_t prec_x0 = (mv0.x) & 0xF;
    uint8_t prec_y0 = (mv0.y) & 0xF;

    uint8_t prec_x1 = (mv1.x) & 0xF;
    uint8_t prec_y1 = (mv1.y) & 0xF;

    uint8_t prec_0_mc_type = (prec_x0 > 0) + ((prec_y0 > 0) << 1);
    uint8_t prec_1_mc_type = (prec_x1 > 0) + ((prec_y1 > 0) << 1);

    dst.y  += x0 + y0 * dst.stride;
    dst.cb += (x0 >> 1) + (y0 >> 1) * dst.stride_c;
    dst.cr += (x0 >> 1) + (y0 >> 1) * dst.stride_c;

    mc_l->bidir0[prec_0_mc_type][log2_pu_w - 1](tmp_buff, ref0_b.y, ref0_b.stride, pu_h, prec_x0, prec_y0, pu_w);
    mc_l->bidir1[prec_1_mc_type][log2_pu_w - 1](dst.y, RCN_CTB_STRIDE, ref1_b.y, ref1_b.stride, tmp_buff, pu_h, prec_x1, prec_y1, pu_w);

    if (ctudec->lmcs_info.lmcs_enabled_flag){
        rcn_lmcs_reshape_luma_blk_lut(dst.y, RCN_CTB_STRIDE, ctudec->lmcs_info.lmcs_lut_fwd_luma, pu_w, pu_h);
    }

    const struct OVBuffInfo ref0_c = derive_ref_buf_c(ref0, mv0,
                                                      pos_x >> 1, pos_y >> 1,
                                                      edge_buff0, edge_buff0_1,
                                                      log2_pu_w, log2_pu_h, log2_ctb_s);

    const struct OVBuffInfo ref1_c = derive_ref_buf_c(ref1, mv1,
                                                      pos_x >> 1, pos_y >> 1,
                                                      edge_buff1, edge_buff1_1,
                                                      log2_pu_w, log2_pu_h, log2_ctb_s);
    prec_x0 = (mv0.x) & 0x1F;
    prec_y0 = (mv0.y) & 0x1F;

    prec_x1 = (mv1.x) & 0x1F;
    prec_y1 = (mv1.y) & 0x1F;

    prec_0_mc_type = (prec_x0 > 0) + ((prec_y0 > 0) << 1);
    prec_1_mc_type = (prec_x1 > 0) + ((prec_y1 > 0) << 1);

    int16_t* ref_data0 = tmp_buff;
    int16_t* ref_data1 = tmp_buff + MAX_PB_SIZE / 2;

    mc_c->bidir0[prec_0_mc_type][log2_pu_w - 1](ref_data0, ref0_c.cb, ref0_c.stride_c, pu_h >> 1, prec_x0, prec_y0, pu_w >> 1);
    mc_c->bidir0[prec_0_mc_type][log2_pu_w - 1](ref_data1, ref0_c.cr, ref0_c.stride_c, pu_h >> 1, prec_x0, prec_y0, pu_w >> 1);

    mc_c->bidir1[prec_1_mc_type][log2_pu_w - 1](dst.cb, RCN_CTB_STRIDE, ref1_c.cb, ref1_c.stride_c, ref_data0, pu_h >> 1, prec_x1, prec_y1, pu_w >> 1);
    mc_c->bidir1[prec_1_mc_type][log2_pu_w - 1](dst.cr, RCN_CTB_STRIDE, ref1_c.cr, ref1_c.stride_c, ref_data1, pu_h >> 1, prec_x1, prec_y1, pu_w >> 1);

}

void
rcn_mcp(OVCTUDec *const ctudec, struct OVBuffInfo dst, int x0, int y0, int log2_pu_w, int log2_pu_h,
        OVMV mv, uint8_t type, uint8_t ref_idx)
{
    struct OVRCNCtx    *const rcn_ctx   = &ctudec->rcn_ctx;
    struct InterDRVCtx *const inter_ctx = &ctudec->drv_ctx.inter_ctx;

    struct MCFunctions *mc_l = &rcn_ctx->rcn_funcs.mc_l;
    struct MCFunctions *mc_c = &rcn_ctx->rcn_funcs.mc_c;

    OVPicture *ref0 = inter_ctx->rpl0[ref_idx];
    OVPicture *ref1 = inter_ctx->rpl1[ref_idx];

    dst.y  += x0 + y0 * dst.stride;
    dst.cb += (x0 >> 1) + (y0 >> 1) * dst.stride_c;
    dst.cr += (x0 >> 1) + (y0 >> 1) * dst.stride_c;

    //TODOlmcs: do not allocate here, maybe use already allocated filter buffers ?
    //Why not use derive_ref_buf_y ?
    uint16_t tmp_buff [RCN_CTB_SIZE];

    const OVFrame *const frame0 =  type ? ref1->frame : ref0->frame;

    const uint16_t *const ref0_y  = (uint16_t *) frame0->data[0];
    const uint16_t *const ref0_cb = (uint16_t *) frame0->data[1];
    const uint16_t *const ref0_cr = (uint16_t *) frame0->data[2];

    int src_stride   = frame0->linesize[0] >> 1;
    int src_stride_c = frame0->linesize[1] >> 1;

    uint8_t log2_ctb_s = ctudec->part_ctx->log2_ctu_s;
    int pos_x = (ctudec->ctb_x << log2_ctb_s) + x0;
    int pos_y = (ctudec->ctb_y << log2_ctb_s) + y0;

    const int pu_w = 1 << log2_pu_w;
    const int pu_h = 1 << log2_pu_h;

    const int pic_w = frame0->width[0];
    const int pic_h = frame0->height[0];

    mv = clip_mv(pos_x, pos_y, pic_w, pic_h, pu_w, pu_h, mv);

    int ref_x = pos_x + (mv.x >> 4);
    int ref_y = pos_y + (mv.y >> 4);

    uint8_t prec_x   = (mv.x) & 0xF;
    uint8_t prec_y   = (mv.y) & 0xF;
    uint8_t prec_x_c = (mv.x) & 0x1F;
    uint8_t prec_y_c = (mv.y) & 0x1F;

    int prec_mc_type   = (prec_x  > 0) + ((prec_y > 0)   << 1);
    int prec_c_mc_type = (prec_x_c > 0) + ((prec_y_c > 0) << 1);

    uint8_t emulate_edge = test_for_edge_emulation(ref_x, ref_y, pic_w, pic_h,
                                                   pu_w, pu_h);;

    const uint16_t *src_y  = &ref0_y [ ref_x       + ref_y        * src_stride];
    const uint16_t *src_cb = &ref0_cb[(ref_x >> 1) + (ref_y >> 1) * src_stride_c];
    const uint16_t *src_cr = &ref0_cr[(ref_x >> 1) + (ref_y >> 1) * src_stride_c];

    /* FIXME
     * Thread synchronization to ensure data is available before usage
     */

    if (emulate_edge){
        int src_off  = REF_PADDING_L * (src_stride) + (REF_PADDING_L);
        int buff_off = REF_PADDING_L * (RCN_CTB_STRIDE) + (REF_PADDING_L);

        emulate_block_border(tmp_buff, (src_y - src_off),
                             RCN_CTB_STRIDE, src_stride,
                             pu_w + QPEL_EXTRA, pu_h + QPEL_EXTRA,
                             ref_x - REF_PADDING_L, ref_y - REF_PADDING_L,
                             pic_w, pic_h);

        src_y = tmp_buff + buff_off;
        src_stride = RCN_CTB_STRIDE;
    }

    mc_l->unidir[prec_mc_type][log2_pu_w](dst.y, RCN_CTB_STRIDE,
                                          src_y, src_stride, pu_h,
                                          prec_x, prec_y, pu_w);

    if (ctudec->lmcs_info.lmcs_enabled_flag){
        rcn_lmcs_reshape_luma_blk_lut(dst.y, RCN_CTB_STRIDE, ctudec->lmcs_info.lmcs_lut_fwd_luma, pu_w, pu_h);
    }


    emulate_edge = test_for_edge_emulation_c(ref_x >> 1, ref_y >> 1, pic_w >> 1, pic_h >> 1,
                                             pu_w >> 1, pu_h >> 1);;

    if (emulate_edge){
        int src_off  = REF_PADDING_C * (src_stride_c) + (REF_PADDING_C);
        int buff_off = REF_PADDING_C * (RCN_CTB_STRIDE) + (REF_PADDING_C);
        emulate_block_border(tmp_buff, (src_cb - src_off),
                             RCN_CTB_STRIDE, src_stride_c,
                             (pu_w >> 1)  + EPEL_EXTRA, (pu_h >> 1) + EPEL_EXTRA,
                             (pos_x >> 1) + (mv.x >> 5) - REF_PADDING_C, (pos_y >> 1) + (mv.y >> 5) - REF_PADDING_C,
                             (pic_w >> 1), (pic_h >> 1));
        src_cb = tmp_buff + buff_off;
        src_stride_c = RCN_CTB_STRIDE;
    }

    mc_c->unidir[prec_c_mc_type][log2_pu_w - 1](dst.cb, RCN_CTB_STRIDE,
                                                src_cb, src_stride_c,
                                                pu_h >> 1, prec_x_c, prec_y_c, pu_w >> 1);

    if (emulate_edge){
        int src_off  = REF_PADDING_C * (frame0->linesize[1] >> 1) + (REF_PADDING_C);
        int buff_off = REF_PADDING_C * (RCN_CTB_STRIDE) + (REF_PADDING_C);
        emulate_block_border(tmp_buff, (src_cr - src_off),
                             RCN_CTB_STRIDE, frame0->linesize[1] >> 1,
                             (pu_w >> 1) + EPEL_EXTRA, (pu_h >> 1) + EPEL_EXTRA,
                             (pos_x >> 1) + (mv.x >> 5) - REF_PADDING_C, (pos_y >> 1) + (mv.y >> 5) - REF_PADDING_C,
                             (pic_w >> 1), (pic_h >> 1));
        src_cr = tmp_buff + buff_off;
        src_stride_c = RCN_CTB_STRIDE;
    }
    mc_c->unidir[prec_c_mc_type][log2_pu_w - 1](dst.cr, RCN_CTB_STRIDE,
                                                src_cr, src_stride_c,
                                                pu_h >> 1, prec_x_c, prec_y_c, pu_w >> 1);
}

void
rcn_mcp_b(OVCTUDec*const lc_ctx, struct OVBuffInfo dst, struct InterDRVCtx *const inter_ctx,
          const OVPartInfo *const part_ctx,
          const OVMV mv0, const OVMV mv1,
          unsigned int x0, unsigned int y0,
          unsigned int log2_pb_w, unsigned int log2_pb_h,
          uint8_t inter_dir, uint8_t ref_idx0, uint8_t ref_idx1)
{
    if (inter_dir == 3) {

        rcn_motion_compensation_b(lc_ctx, dst, x0, y0, log2_pb_w, log2_pb_h, mv0, mv1, ref_idx0, ref_idx1);

    } else if (inter_dir & 0x2) {

        rcn_mcp(lc_ctx, dst, x0, y0, log2_pb_w, log2_pb_h, mv1, 1, ref_idx1);

    } else if (inter_dir & 0x1) {

        rcn_mcp(lc_ctx, dst, x0, y0, log2_pb_w, log2_pb_h, mv0, 0, ref_idx0);

    }
}

//TODOciip: do not define here
#define BIT_DEPTH 10
#define ov_clip_pixel(a) ov_clip_uintp2(a, BIT_DEPTH)

static void
put_weighted_ciip_pixels(uint16_t* dst, ptrdiff_t dststride,
                      const uint16_t* src_intra, const uint16_t* src_inter, ptrdiff_t srcstride,
                      int width, int height, int wt)
{   
    int x, y;
    int shift = 14 - BIT_DEPTH + 2;
    int offset = 2 + ( 1 << (shift - 1));
    for (y = 0; y < height; y++) {
        for (x = 0; x < width; ++x) {
            dst[x] = ov_clip_pixel( (src_intra[x] * wt + ((src_inter[x] * (4 - wt)) << (14 - BIT_DEPTH)) + offset) >> shift );
        }
        src_intra += srcstride;
        src_inter += srcstride;
        dst += dststride;
    }
}

/*Apply planar intra modes current
*/
void
rcn_ciip_b(OVCTUDec*const ctudec,
           const OVMV mv0, const OVMV mv1,
           unsigned int x0, unsigned int y0,
           unsigned int log2_pb_w, unsigned int log2_pb_h,
           uint8_t inter_dir, uint8_t ref_idx0, uint8_t ref_idx1)
{

    struct InterDRVCtx *const inter_ctx = &ctudec->drv_ctx.inter_ctx;
    const OVPartInfo *const part_ctx = ctudec->part_ctx;
    struct OVBuffInfo tmp_inter;
    uint16_t tmp_inter_l [RCN_CTB_SIZE], tmp_inter_cb[RCN_CTB_SIZE], tmp_inter_cr[RCN_CTB_SIZE] ;
    tmp_inter.y  = &tmp_inter_l [0];
    tmp_inter.cb = &tmp_inter_cb[0];
    tmp_inter.cr = &tmp_inter_cr[0];
    tmp_inter.stride   = RCN_CTB_STRIDE;
    tmp_inter.stride_c = RCN_CTB_STRIDE;
    rcn_mcp_b(ctudec, tmp_inter, inter_ctx, part_ctx, mv0, mv1, x0, y0, log2_pb_w, log2_pb_h, 
        inter_dir, ref_idx0, ref_idx1);

    struct OVBuffInfo tmp_intra;
    uint16_t tmp_intra_l [RCN_CTB_SIZE], tmp_intra_cb[RCN_CTB_SIZE], tmp_intra_cr[RCN_CTB_SIZE] ;
    tmp_intra.y  = &tmp_intra_l [0];
    tmp_intra.cb = &tmp_intra_cb[0];
    tmp_intra.cr = &tmp_intra_cr[0];
    tmp_intra.stride   = RCN_CTB_STRIDE;
    tmp_intra.stride_c = RCN_CTB_STRIDE;
    vvc_intra_pred(&ctudec->rcn_ctx, &tmp_intra, OVINTRA_PLANAR, x0, y0, log2_pb_w, log2_pb_w);
    vvc_intra_pred_chroma(&ctudec->rcn_ctx, &tmp_intra, OVINTRA_PLANAR, x0 >> 1, y0 >> 1, log2_pb_w - 1, log2_pb_w - 1);

    //TODOciip: chose correct value in function of neighbors
    int wt = 0;
    struct OVBuffInfo dst = ctudec->rcn_ctx.ctu_buff;
    dst.y  += x0 + y0 * dst.stride;
    tmp_intra.y  += x0 + y0 * tmp_intra.stride;
    tmp_inter.y  += x0 + y0 * tmp_inter.stride;
    put_weighted_ciip_pixels(dst.y, dst.stride, tmp_intra.y, tmp_inter.y, tmp_inter.stride,
                       1 << log2_pb_w, 1 << log2_pb_h, wt);

    dst.cb += (x0 >> 1) + (y0 >> 1) * dst.stride_c;
    tmp_intra.cb += (x0 >> 1) + (y0 >> 1) * tmp_intra.stride_c;
    tmp_inter.cb += (x0 >> 1) + (y0 >> 1) * tmp_inter.stride_c;
    put_weighted_ciip_pixels(dst.cb, dst.stride_c, tmp_intra.cb, tmp_inter.cb, tmp_inter.stride_c,
                       1 << log2_pb_w, 1 << log2_pb_h, wt);

    dst.cr += (x0 >> 1) + (y0 >> 1) * dst.stride_c;
    tmp_intra.cr += (x0 >> 1) + (y0 >> 1) * tmp_intra.stride_c;
    tmp_inter.cr += (x0 >> 1) + (y0 >> 1) * tmp_inter.stride_c;
    put_weighted_ciip_pixels(dst.cr, dst.stride_c, tmp_intra.cr, tmp_inter.cr, tmp_inter.stride_c,
                       1 << (log2_pb_w - 1), 1 << (log2_pb_h - 1), wt);
}
