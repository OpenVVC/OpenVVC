
#include <string.h>
#include "ctudec.h"
#include "rcn_structures.h"
#include "bitdepth.h"

static void
rcn_ibc_l(OVCTUDec *const ctu_dec,
        int16_t x0, int16_t y0,
        uint8_t log2_cu_w, uint8_t log2_cu_h,
        uint8_t log2_ctu_s,
        IBCMV mv)
{
    uint8_t ctb_msk = (256 * 128 >> (2*log2_ctu_s)) - 1;
    uint16_t ctb_pos = (ctu_dec->ctb_x & ctb_msk) << log2_ctu_s;
    uint16_t msk_h = ((1 << 8 + 7) >> log2_ctu_s) - 1;
    uint16_t msk_v = (1 << log2_ctu_s) - 1;

    int16_t ref_x = ctb_pos + x0 + mv.x;
    int16_t ref_y = y0 + mv.y;

    uint8_t cu_w = 1 << log2_cu_w;
    uint8_t cu_h = 1 << log2_cu_h;;

    ref_x &= msk_h;
    ref_y &= msk_v;

    uint16_t ibc_stride = 256 * 128 >> log2_ctu_s;
    uint8_t req_wrp = ref_x + cu_w > ibc_stride;

    if (!req_wrp) {
        const OVSample *src = (OVSample *) (ctu_dec->rcn_ctx.ctu_buff.y) + ref_x - ctb_pos + ref_y * ctu_dec->rcn_ctx.ctu_buff.stride;
        OVSample *dst = (OVSample *) (ctu_dec->rcn_ctx.ctu_buff.y) + x0 + y0 * ctu_dec->rcn_ctx.ctu_buff.stride;
        for(int i = 0; i < cu_h; ++i) {
            memcpy(dst, src, sizeof(OVSample) << log2_cu_w);
            dst += ctu_dec->rcn_ctx.ctu_buff.stride;
            src += ctu_dec->rcn_ctx.ctu_buff.stride;
        }
    } else {
        uint8_t size1 = ref_x + cu_w - ibc_stride;
        uint8_t size0 = cu_w - size1;
        const OVSample *src = (OVSample *) (ctu_dec->rcn_ctx.ctu_buff.y) + ref_x - ctb_pos + ref_y * ctu_dec->rcn_ctx.ctu_buff.stride;
        OVSample *dst = (OVSample *) (ctu_dec->rcn_ctx.ctu_buff.y) + x0 + y0 * ctu_dec->rcn_ctx.ctu_buff.stride;
        for(int i = 0; i < cu_h; ++i) {
            memcpy(dst, src, sizeof(OVSample) * size0);
            dst += ctu_dec->rcn_ctx.ctu_buff.stride;
            src += ctu_dec->rcn_ctx.ctu_buff.stride;
        }

        src = (OVSample *) (ctu_dec->rcn_ctx.ctu_buff.y) - ctb_pos + ref_y * ctu_dec->rcn_ctx.ctu_buff.stride;
        dst = (OVSample *) (ctu_dec->rcn_ctx.ctu_buff.y) + x0 + size0 + y0 * ctu_dec->rcn_ctx.ctu_buff.stride;
        for(int i = 0; i < cu_h; ++i) {
            memcpy(dst, src, sizeof(OVSample) * size1);
            dst += ctu_dec->rcn_ctx.ctu_buff.stride;
            src += ctu_dec->rcn_ctx.ctu_buff.stride;
        }
    }

}

static void
rcn_ibc_c(OVCTUDec *const ctu_dec,
          int16_t x0, int16_t y0,
          uint8_t log2_cu_w, uint8_t log2_cu_h,
          uint8_t log2_ctu_s,
          IBCMV mv)
{
    uint8_t ctb_msk = (256 * 128 >> (2*log2_ctu_s)) - 1;
    uint16_t ctb_pos = (ctu_dec->ctb_x & ctb_msk) << log2_ctu_s;
    uint16_t msk_h = ((1 << 8 + 7) >> log2_ctu_s) - 1;
    uint16_t msk_v = (1 << log2_ctu_s) - 1;

    int16_t ref_x = ctb_pos + x0 + mv.x;
    int16_t ref_y = y0 + mv.y;

    uint8_t cu_w = 1 << log2_cu_w;
    uint8_t cu_h = 1 << log2_cu_h;;

    ref_x &= msk_h;
    ref_y &= msk_v;

    uint16_t ibc_stride = 256 * 128 >> log2_ctu_s;
    uint8_t req_wrp = ref_x + cu_w > ibc_stride;

    cu_w >>= 1;
    cu_h >>= 1;;
    ref_x &= msk_h;
    ref_y &= msk_v;
    if (!req_wrp) {
        const OVSample *src = (OVSample *) (ctu_dec->rcn_ctx.ctu_buff.cb) + (ref_x - ctb_pos  >> 1)+ (ref_y >> 1) * ctu_dec->rcn_ctx.ctu_buff.stride_c;
        OVSample *dst = (OVSample *) (ctu_dec->rcn_ctx.ctu_buff.cb) + (x0 >> 1) + (y0 >> 1) * ctu_dec->rcn_ctx.ctu_buff.stride_c;
        for(int i = 0; i < cu_h; ++i) {
            memcpy(dst, src, sizeof(OVSample) << log2_cu_w >> 1);
            dst += ctu_dec->rcn_ctx.ctu_buff.stride_c;
            src += ctu_dec->rcn_ctx.ctu_buff.stride_c;
        }
        src = (OVSample *) (ctu_dec->rcn_ctx.ctu_buff.cr) + (ref_x- ctb_pos >> 1) + (ref_y >> 1) * ctu_dec->rcn_ctx.ctu_buff.stride_c;
        dst = (OVSample *) (ctu_dec->rcn_ctx.ctu_buff.cr) + (x0 >> 1) + (y0 >> 1) * ctu_dec->rcn_ctx.ctu_buff.stride_c;
        for(int i = 0; i < cu_h; ++i) {
            memcpy(dst, src, sizeof(OVSample) << log2_cu_w >> 1);
            dst += ctu_dec->rcn_ctx.ctu_buff.stride_c;
            src += ctu_dec->rcn_ctx.ctu_buff.stride_c;
        }
    } else {
        uint8_t size1 = (ref_x >> 1) + cu_w - (ibc_stride >> 1);
        uint8_t size0 = cu_w - size1;
        const OVSample *src = (OVSample *) (ctu_dec->rcn_ctx.ctu_buff.cb) + (ref_x - ctb_pos  >> 1)+ (ref_y >> 1) * ctu_dec->rcn_ctx.ctu_buff.stride_c;
        OVSample *dst = (OVSample *) (ctu_dec->rcn_ctx.ctu_buff.cb) + (x0 >> 1) + (y0 >> 1) * ctu_dec->rcn_ctx.ctu_buff.stride_c;
        for(int i = 0; i < cu_h; ++i) {
            memcpy(dst, src, sizeof(OVSample) * size0);
            dst += ctu_dec->rcn_ctx.ctu_buff.stride_c;
            src += ctu_dec->rcn_ctx.ctu_buff.stride_c;
        }

        src = (OVSample *) (ctu_dec->rcn_ctx.ctu_buff.cb) - (ctb_pos >> 1) + (ref_y >> 1) * ctu_dec->rcn_ctx.ctu_buff.stride_c;
        dst = (OVSample *) (ctu_dec->rcn_ctx.ctu_buff.cb) + (x0 >> 1) + size0 + (y0 >> 1) * ctu_dec->rcn_ctx.ctu_buff.stride_c;
        for(int i = 0; i < cu_h; ++i) {
            memcpy(dst, src, sizeof(OVSample) * size1);
            dst += ctu_dec->rcn_ctx.ctu_buff.stride_c;
            src += ctu_dec->rcn_ctx.ctu_buff.stride_c;
        }

        src = (OVSample *) (ctu_dec->rcn_ctx.ctu_buff.cr) + (ref_x - ctb_pos >> 1) + (ref_y >> 1) * ctu_dec->rcn_ctx.ctu_buff.stride_c;
        dst = (OVSample *) (ctu_dec->rcn_ctx.ctu_buff.cr) + (x0 >> 1) + (y0 >> 1) * ctu_dec->rcn_ctx.ctu_buff.stride_c;
        for(int i = 0; i < cu_h; ++i) {
            memcpy(dst, src, sizeof(OVSample) * size0);
            dst += ctu_dec->rcn_ctx.ctu_buff.stride_c;
            src += ctu_dec->rcn_ctx.ctu_buff.stride_c;
        }

        src = (OVSample *) (ctu_dec->rcn_ctx.ctu_buff.cr) - (ctb_pos >> 1) + (ref_y >> 1) * ctu_dec->rcn_ctx.ctu_buff.stride_c;
        dst = (OVSample *) (ctu_dec->rcn_ctx.ctu_buff.cr) + (x0 >> 1) + size0 + (y0 >> 1) * ctu_dec->rcn_ctx.ctu_buff.stride_c;
        for(int i = 0; i < cu_h; ++i) {
            memcpy(dst, src, sizeof(OVSample) * size1);
            dst += ctu_dec->rcn_ctx.ctu_buff.stride_c;
            src += ctu_dec->rcn_ctx.ctu_buff.stride_c;
        }
    }
}

void
BD_DECL(rcn_init_ibc)(struct RCNFunctions *const rcn_func)
{
     rcn_func->rcn_ibc_l = &rcn_ibc_l;
     rcn_func->rcn_ibc_c = &rcn_ibc_c;
}

