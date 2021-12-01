#include <stdint.h>
#include <string.h>

#include "data_rcn_transform.h"
#include "rcn.h"
#include "ovutils.h"
#include "drv.h"
#include "ctudec.h"

static const uint8_t lfnst_mode_map[67+28] =
{//0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94
   0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
};


static inline int
derive_wide_angular_mode2(int8_t log2_pb_w, int8_t log2_pb_h, int pred_mode)
{
    static const uint8_t mode_shift_lut[6] = {0, 6, 10, 12, 14, 15};
    int mode_shift = mode_shift_lut[OVABS(log2_pb_w - log2_pb_h)];
    if (log2_pb_w > log2_pb_h && pred_mode < 2 + mode_shift) {
        pred_mode += (OVINTRA_VDIA - 1);
    } else if (log2_pb_h > log2_pb_w && pred_mode > OVINTRA_VDIA - mode_shift) {
        pred_mode -= (OVINTRA_VDIA + 1);
    }
    return pred_mode;
}

/* FIXME move pare to mode derviation */
static int8_t
derive_lfnst_mode_l(OVCTUDec *const ctudec,
                   int log2_tb_w, int log2_tb_h,
                   int x0, int y0)
{
    const OVPartInfo *const part_ctx = ctudec->part_ctx;
    int log2_min_cb_s = part_ctx->log2_min_cb_s;
    int x_pu = x0 >> log2_min_cb_s;
    int y_pu = y0 >> log2_min_cb_s;
    int8_t intra_mode = ctudec->drv_ctx.intra_info.luma_modes[x_pu + (y_pu << 5)];

    if (intra_mode > OVINTRA_DC) {
        intra_mode = derive_wide_angular_mode2(log2_tb_w, log2_tb_h, intra_mode);
    }

    intra_mode = intra_mode < 0 ? intra_mode + 14 + 67: intra_mode >= 67
        ? intra_mode + 14 : intra_mode;

    return intra_mode;
}

static int8_t
derive_lfnst_mode_c(OVCTUDec *const ctudec,
                   int log2_tb_w, int log2_tb_h,
                   int x0, int y0, uint8_t intra_mode_c)
{
    const OVPartInfo *const part_ctx = ctudec->part_ctx_c;
    int log2_min_cb_s = part_ctx->log2_min_cb_s;
    int x_pu = x0 >> log2_min_cb_s;
    int y_pu = y0 >> log2_min_cb_s;
    /* FIXME chroma mode */
    int8_t intra_mode = intra_mode_c;//ctudec->pred_ctx.intra_modes_chroma[x_pu + (y_pu << 5)];
    uint8_t nb_pb_w = (1 << log2_tb_w) >> log2_min_cb_s;
    uint8_t nb_pb_h = (1 << log2_tb_h) >> log2_min_cb_s;
    uint8_t luma_mode = ctudec->drv_ctx.intra_info.luma_modes[(x_pu + ((y_pu + (nb_pb_h >> 1)) << 5) + (nb_pb_w >> 1))];

    intra_mode = intra_mode == OVINTRA_DM_CHROMA ||
        (intra_mode >= OVINTRA_LM_CHROMA && intra_mode <= OVINTRA_MDLM_TOP)
        ? luma_mode : intra_mode;

    /* We need to derive angular mode as used in intra prediction */
    if (intra_mode > OVINTRA_DC) {
        intra_mode = derive_wide_angular_mode2(log2_tb_w, log2_tb_h, intra_mode);
    }

    /* FIXME understatnd this check */
    intra_mode = intra_mode < 0 ? intra_mode + 14 + 67 : intra_mode >= 67
        ? intra_mode + 14 : intra_mode;

    return intra_mode;
}

void
process_lfnst(OVCTUDec *const ctudec,
              int16_t *dst, const int16_t *src,
              int log2_tb_w, int log2_tb_h,
              int x0, int y0, uint8_t lfnst_idx)
{
    struct RCNFunctions *const rcnFunc = &ctudec->rcn_funcs;

    int8_t intra_mode = derive_lfnst_mode_c(ctudec, log2_tb_w, log2_tb_h,
                                            x0, y0, ctudec->intra_mode_c);

    uint8_t lfnst_mode = lfnst_mode_map[intra_mode];

    uint8_t need_transpose = ((intra_mode < 67) && (intra_mode > OVINTRA_DIA)) || (intra_mode >= 67 + 14);

    uint64_t scan_map = 0xfbe7ad369c258140;

    int16_t tmp[16];

    for (int i = 0; i < 16; ++i) {
        tmp[i] = src[scan_map & 0xF];
        scan_map >>= 4;
    }

    /* FIXME reduce memset since limited transform */
    memset(dst, 0, sizeof(int16_t) << (log2_tb_w + log2_tb_h));

    const int8_t *lfnst_matrix = lfnst[(log2_tb_w >= 3 && log2_tb_h >= 3)][lfnst_mode][lfnst_idx];
    (rcnFunc->lfnst.func[need_transpose][(log2_tb_w >= 3 && log2_tb_h >= 3)])(tmp, dst, lfnst_matrix, log2_tb_w, log2_tb_h);
}

void
process_lfnst_luma(OVCTUDec *const ctudec,
                   int16_t *dst, const int16_t *src,
                   int log2_tb_w, int log2_tb_h,
                   int x0, int y0, uint8_t lfnst_idx)
{
    struct RCNFunctions *const rcnFunc = &ctudec->rcn_funcs;

    int8_t intra_mode = derive_lfnst_mode_l(ctudec, log2_tb_w, log2_tb_h, x0, y0);

    uint8_t need_transpose = (intra_mode < 67  && intra_mode > OVINTRA_DIA) || intra_mode >= 67 + 14;

    int lfnst_mode = lfnst_mode_map[intra_mode];
    uint64_t scan_map = 0xfbe7ad369c258140;
    int16_t tmp[16];

    for (int i = 0; i < 16; ++i) {
        tmp[i] = src[scan_map & 0xF];
        scan_map >>= 4;
    }

    /* FIXME reduce memset since limited transform */
    memset(dst, 0, sizeof(int16_t) << (log2_tb_w + log2_tb_h));

    const int8_t *lfnst_matrix = lfnst[(log2_tb_w >= 3 && log2_tb_h >= 3)][lfnst_mode][lfnst_idx];
    (rcnFunc->lfnst.func[need_transpose][(log2_tb_w >= 3 && log2_tb_h >= 3)])(tmp, dst, lfnst_matrix, log2_tb_w, log2_tb_h);
}

void
process_lfnst_luma_isp(OVCTUDec *const ctudec,
                   int16_t *dst, const int16_t *src,
                   int log2_tb_w, int log2_tb_h,
                   int log2_cb_w, int log2_cb_h,
                   int x0, int y0, uint8_t lfnst_idx)
{
    struct RCNFunctions *const rcnFunc = &ctudec->rcn_funcs;

    int8_t intra_mode = derive_lfnst_mode_l(ctudec, log2_cb_w, log2_cb_h, x0, y0);

    uint8_t need_transpose = (intra_mode < 67  && intra_mode > OVINTRA_DIA) || intra_mode >= 67 + 14;

    int lfnst_mode = lfnst_mode_map[intra_mode];
    uint64_t scan_map = 0xfbe7ad369c258140;
    int16_t tmp[16];

    for (int i = 0; i < 16; ++i) {
        tmp[i] = src[scan_map & 15];
        scan_map >>= 4;
    }

    /* FIXME reduce memset since limited transform */
    memset(dst, 0, sizeof(int16_t) << (log2_tb_w + log2_tb_h));

    const int8_t *lfnst_matrix = lfnst[(log2_tb_w >= 3 && log2_tb_h >= 3)][lfnst_mode][lfnst_idx];
    (rcnFunc->lfnst.func[need_transpose][(log2_tb_w >= 3 && log2_tb_h >= 3)])(tmp, dst, lfnst_matrix, log2_tb_w, log2_tb_h);
}
