#include <stddef.h>
#include <stdint.h>
#include <arm_neon.h>

#include "ovutils.h"
#include "rcn_structures.h"
#include "stdint.h"


static const uint16_t vvc_pdpc_w[3][128] = {
{32,  8,  2,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{32, 16,  8,  4, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{32, 32, 16, 16, 8, 8, 4, 4, 2, 2, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
};

static const uint8_t vvc_pdpc_w_sh[3][128] = {
{5,  3,  1,  0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F},
{5, 4, 3, 2, 1, 0, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F},
{5, 5, 4, 4, 3, 3, 2, 2, 1, 1, 0, 0, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F},
};

void
vvc_intra_dc_pdpc_neon(const uint16_t *const src_above,
                      const uint16_t *const src_left,
                      uint16_t *const dst, ptrdiff_t dst_stride,
                      int log2_pb_w, int log2_pb_h)
{
    int idx;
    uint16_t *_dst = dst;
    uint32_t dc_val = 0;

    const int denom  = (log2_pb_w == log2_pb_h) ? (log2_pb_w + 1)
                                                : OVMAX(log2_pb_w,log2_pb_h);
    const int shift  = denom;
    const int offset = ((1 << denom) >> 1);
    const int pdpc_scale = (log2_pb_w + log2_pb_h - 2) >> 2;
    #if 1
    const uint8_t *pdpc_sh = vvc_pdpc_w_sh[pdpc_scale];
    const uint16_t *pdpc_w = vvc_pdpc_w[pdpc_scale];
    #else
    const uint16_t *pdpc_w = vvc_pdpc_w[pdpc_scale];
    #endif
    int nb_neon_vec8_w = ((1 << log2_pb_w) + ((1 << 3) - 1)) >> 3;
    // int nb_neon_vec8_h = ((1 << log2_pb_h) + ((1 << 3) - 1)) >> 3;

    /* FIXME don't forget to optimize those loop when specializing by size*/
    if (log2_pb_w >= log2_pb_h){
        for (idx = 0; idx < (1 << log2_pb_w); idx++ ){
            dc_val += src_above[1 + idx];
        }
    }

    if (log2_pb_w <= log2_pb_h){
        for (idx = 0; idx < (1 << log2_pb_h); idx++ ){
            dc_val += src_left[1 + idx];
        }
    }

    dc_val = (dc_val + offset) >> shift;

    uint16x8_t dc_v  = vdupq_n_u16(dc_val);
    uint16x8_t add_v = vdupq_n_u16(32);

    for (int y = 0; y < (1 << log2_pb_h); ++y){
        int x;
        const int16_t l_val = src_left[y + 1];
        uint16x8_t l_v  = vdupq_n_u16(l_val);
        int y_wgh = pdpc_sh[y];

        /*FIXME PDPC weights are derived from log2 values all
           multiplications could be replaced by shift operations
           It will be easier to derive shift when specializing according
           to sizes*/
        for (x = 0; x < nb_neon_vec8_w; x++){
            uint16x8_t xl, yt, x_v, t_v, w_x, w_y;
            uint16x8_t pdpc_rnd, out_v;
            uint16x8_t tst;
            x_v = vld1q_u16((const uint16_t *) (pdpc_w    + 8 * x));
            t_v = vld1q_u16((const uint16_t *) (src_above + 8 * x + 1));
            tst = vshlq_n_u16(dc_v, 6);

            w_x = vshlq_n_u16(dc_v, y_wgh);
            w_y = vmulq_u16(dc_v, x_v);
            tst = vqsubq_u16(tst, w_x);
            tst = vqsubq_u16(tst, w_y);
            tst = vqaddq_u16(tst, add_v);

            xl = vmulq_u16(l_v, x_v);
            yt = vshlq_n_u16(t_v, y_wgh);

            pdpc_rnd = vqaddq_u16(xl, yt);

            tst = vqaddq_u16(tst, pdpc_rnd);
            out_v = vshrq_n_u16(tst,6);

            out_v = vminq_u16(out_v, vdupq_n_u16(1023));
            out_v = vmaxq_u16(out_v, vdupq_n_u16(0));

            vst1q_u16((_dst + 8 * x), out_v);
        }
        _dst += dst_stride;
    }
}

#define CUT_PDPC 0
void
vvc_intra_planar_pdpc_neon(const uint16_t *const src_above,
                          const uint16_t *const src_left,
                          uint16_t *const dst, ptrdiff_t dst_stride,
                          int log2_pb_w, int log2_pb_h)
{

    uint16_t *_dst = dst; //TODO template according to bitdepth;
    const uint32_t width  = 1 << log2_pb_w;
    const uint32_t height = 1 << log2_pb_h;
    //no need to check max since there is no pdpc call in isp
    // when these sizes are reached
    const uint32_t w_scale = OVMAX(1, log2_pb_w);
    const uint32_t h_scale = OVMAX(1, log2_pb_h);
    //const uint32_t w_scale = log2_pb_w;
    //const uint32_t h_scale = log2_pb_h;
    // const uint32_t s_shift = w_scale + h_scale + 1;
    // const uint32_t offset  = 1 << (w_scale + h_scale);

    const uint8_t pdpc_scale = (log2_pb_w + log2_pb_h - 2) >> 2;

    int16_t t_row[128],  r_col[128], l_col[128];//b_row[128],
    const uint16_t *pdpc_w = vvc_pdpc_w[pdpc_scale];
    // const int16_t bl_val = src_left[height + 1];
    // const int16_t tr_val = src_above[width + 1];
    #if CUT_PDPC
    const int pdpc_stop_w = OVMIN(3 << pdpc_scale, width);
    const int pdpc_stop_h = OVMIN(3 << pdpc_scale, height);
    #endif
    int x,y;
    int nb_neon_vec8_w = ((1 << log2_pb_w) + ((1 << 3) - 1)) >> 3;
    int nb_neon_vec8_h = ((1 << log2_pb_h) + ((1 << 3) - 1)) >> 3;

    uint16x8_t tr_val_v = vdupq_n_u16(src_above[width + 1]);
    uint16x8_t bl_val_v = vdupq_n_u16(src_left[height + 1]);

    for(x = 0; x < nb_neon_vec8_w; ++x){
        uint16x8_t tr_v;//, br_v;
        uint16x8_t src_a = vld1q_u16((uint16_t *)(src_above + x * 8 + 1));
        //br_v = _mm_sub_epi16(bl_val_v, src_a);
        tr_v = vshlq_n_u16(src_a, h_scale);
        //_mm_storeu_si128((__m128i*) &b_row[y*8], br_v);
        vst1q_u16((uint16_t*) &t_row[x*8], tr_v);
    }
    for(y = 0; y < nb_neon_vec8_h; ++y){
        uint16x8_t lc_v, rc_v;
        uint16x8_t src_l = vld1q_u16((uint16_t *)(src_left + y * 8 + 1));
        rc_v = vsubq_u16(tr_val_v, src_l);
        lc_v = vshlq_n_u16(src_l, w_scale);
        vst1q_u16((uint16_t*) &r_col[y*8], rc_v);
        vst1q_u16((uint16_t*) &l_col[y*8], lc_v);
    }
    int rc_scale = h_scale > w_scale ? h_scale - w_scale : 0;
    int tr_scale = w_scale > h_scale ? w_scale - h_scale : 0;
    int max_scale = OVMAX(h_scale,w_scale);
    uint16x8_t rnd_v = vdupq_n_u16(1 << max_scale);
    // static uint16_t vect[8] =  {8, 7, 6, 5, 4, 3, 2, 1};
    static uint16_t vect[8] =  {1, 2, 3, 4, 5, 6, 7, 8};
    uint16x8_t mul_col = vld1q_u16(vect);
    uint16x8_t add_v = vdupq_n_u16(32);
    for (y = 0; y < height; ++y){
        uint16x8_t rcol_v = vdupq_n_u16(r_col[y]);
        uint16x8_t lcol_v = vdupq_n_u16(l_col[y]);
        uint16x8_t rcol_mul_v = vmulq_u16(rcol_v, mul_col);
        uint16x8_t rcol_x ;
        /* pdpc var */
        int y_wgh = pdpc_w[y];
        int32_t l_val = (int32_t)src_left[y + 1];
        uint16x8_t l_v = vdupq_n_u16(l_val);
        uint16x8_t y_v = vdupq_n_u16(y_wgh);
        for (x = 0; x < nb_neon_vec8_w; ++x) {
            uint16x8_t out_lo;//, out_hi;
            uint16x8_t rc_v_lo;//, rc_v_hi;
            uint16x8_t tr_v_lo;//, tr_v_hi;
            uint16x8_t src_v_lo;//, src_v_hi;
            uint16x8_t str_v_lo;//, str_v_hi;

            uint16x8_t src_a = vld1q_u16(src_above + x * 8 + 1);
            uint16x8_t tr_v = vld1q_u16((uint16_t *)(t_row + x * 8));

            tr_v = vsubq_u16(tr_v, src_a);
            tr_v = vaddq_u16(tr_v, bl_val_v);

            vst1q_u16((uint16_t *)(&t_row[x*8]), tr_v);

            rcol_x = vaddq_u16(lcol_v, rcol_mul_v);
            rcol_mul_v = vaddq_u16(rcol_mul_v, vshlq_n_u16(rcol_v, 3));


            rc_v_lo = rcol_x;

            tr_v_lo = tr_v;

            src_v_lo = vshlq_n_u16(rc_v_lo, rc_scale);
            str_v_lo = vshlq_n_u16(tr_v_lo, tr_scale);

            out_lo = vaddq_u16(src_v_lo, str_v_lo);

            out_lo = vaddq_u16(out_lo, rnd_v);

            out_lo = vshrq_n_u16(out_lo, max_scale + 1);

            /* FIXME applying PDPC on the whole PU is useless
               since only max 12 pel cols from left and 12 pel
               rows from top require pdpc processing*/
            uint16x8_t xl, yt, x_v, w_x, w_y;//t_v,
            uint16x8_t pdpc_rnd, out_v;
            uint16x8_t tst;
            x_v = vld1q_u16(pdpc_w    + 8 * x);

            tst = vshlq_n_u16(out_lo, 6);

            w_x = vmulq_u16(out_lo, y_v);
            w_y = vmulq_u16(out_lo, x_v);

            tst = vsubq_u16(tst, w_x);
            tst = vsubq_u16(tst, w_y);
            tst = vaddq_u16(tst, add_v);

            xl = vmulq_u16(x_v, l_v);
            yt = vmulq_u16(y_v, src_a);

            pdpc_rnd = vaddq_u16(xl, yt);

            tst = vaddq_u16(tst, pdpc_rnd);

            out_v = vshrq_n_u16(tst,6);

            out_v = vminq_u16(out_v, vdupq_n_u16(1023));
            out_v = vmaxq_u16(out_v, vdupq_n_u16(0));

            vst1q_u16((_dst + 8 * x), out_v);
        }
        _dst += dst_stride;
    }
}

void
vvc_intra_planar_pdpc_2_neon(const uint16_t *const src_above,
                          const uint16_t *const src_left,
                          uint16_t *const dst, ptrdiff_t dst_stride,
                          int log2_pb_w, int log2_pb_h)
{

    uint16_t *_dst = dst; //TODO template according to bitdepth;
    const uint32_t width  = 1 << log2_pb_w;
    const uint32_t height = 1 << log2_pb_h;
    //no need to check max since there is no pdpc call in isp
    // when these sizes are reached
    const uint32_t w_scale = OVMAX(1, log2_pb_w);
    const uint32_t h_scale = OVMAX(1, log2_pb_h);
    //const uint32_t w_scale = log2_pb_w;
    //const uint32_t h_scale = log2_pb_h;
    const uint32_t s_shift = w_scale + h_scale + 1;
    const uint32_t offset  = 1 << (w_scale + h_scale);

    const uint8_t pdpc_scale = (log2_pb_w + log2_pb_h - 2) >> 2;

    int16_t t_row[128], r_col[128], l_col[128]; //b_row[128],
    const uint16_t *pdpc_w = vvc_pdpc_w[pdpc_scale];
    // const int16_t bl_val = src_left[height + 1];
    // const int16_t tr_val = src_above[width + 1];
    #if CUT_PDPC
    const int pdpc_stop_w = OVMIN(3 << pdpc_scale, width);
    const int pdpc_stop_h = OVMIN(3 << pdpc_scale, height);
    #endif
    int x,y;
    int nb_neon_vec8_w = ((1 << log2_pb_w) + ((1 << 3) - 1)) >> 3;
    int nb_neon_vec8_h = ((1 << log2_pb_h) + ((1 << 3) - 1)) >> 3;

    uint16x8_t tr_val_v = vdupq_n_u16(src_above[width + 1]);
    uint16x8_t bl_val_v = vdupq_n_u16(src_left[height + 1]);

    for(x = 0; x < nb_neon_vec8_w; ++x){
        uint16x8_t tr_v;//, br_v;
        uint16x8_t src_a = vld1q_u16((src_above + x * 8 + 1));
        //br_v = _mm_sub_epi16(bl_val_v, src_a);
        tr_v = vshlq_n_u16(src_a, h_scale);
        //_mm_storeu_si128((__m128i*) &b_row[y*8], br_v);
        vst1q_u16((uint16_t*) &t_row[x*8], tr_v);
    }
    for(y = 0; y < nb_neon_vec8_h; ++y){
        uint16x8_t lc_v, rc_v;
        uint16x8_t src_l = vld1q_u16((src_left + y * 8 + 1));
        rc_v = vsubq_u16(tr_val_v, src_l);
        lc_v = vshlq_n_u16(src_l, w_scale);
        vst1q_u16((uint16_t*) &r_col[y*8], rc_v);
        vst1q_u16((uint16_t*) &l_col[y*8], lc_v);
    }
    uint32x4_t rnd_v = vdupq_n_u32(offset);
    static uint16_t vect[8] =  {1, 2, 3, 4, 5, 6, 7, 8};
    uint16x8_t mul_col = vld1q_u16(vect);
    uint16x8_t add_v = vdupq_n_u16(32);
    for (y = 0; y < height; ++y){
        uint16x8_t rcol_v = vdupq_n_u16(r_col[y]);
        uint16x8_t lcol_v = vdupq_n_u16(l_col[y]);
        uint16x8_t rcol_mul_v = vmulq_u16(rcol_v, mul_col);
        uint16x8_t rcol_x ;
        /* pdpc var */
        int y_wgh = pdpc_w[y];
        int32_t l_val = (int32_t)src_left[y + 1];
        uint16x8_t l_v = vdupq_n_u16(l_val);
        uint16x8_t y_v = vdupq_n_u16(y_wgh);
        for (x = 0; x < nb_neon_vec8_w; ++x) {
            uint16x8_t out;
            uint32x4_t out_lo, out_hi;
            uint32x4_t rc_v_lo, rc_v_hi;
            uint32x4_t tr_v_lo, tr_v_hi;
            uint32x4_t src_v_lo, src_v_hi;
            uint32x4_t str_v_lo, str_v_hi;

            uint16x8_t src_a = vld1q_u16(src_above + x * 8 + 1);
            //__m128i br_v = _mm_loadu_si128((b_row + x * 8));
            uint16x8_t tr_v = vld1q_u16((uint16_t *)(t_row + x * 8));

            tr_v = vsubq_u16(tr_v, src_a);
            tr_v = vaddq_u16(tr_v, bl_val_v);

            vst1q_u16((uint16_t*) &t_row[x*8], tr_v);

            rcol_x = vaddq_u16(lcol_v, rcol_mul_v);
            rcol_mul_v = vaddq_u16(rcol_mul_v, vshlq_n_u16(rcol_v, 3));

            rc_v_lo = (uint32x4_t) vzip1q_u16(rcol_x, vdupq_n_u16(0));
            rc_v_hi = (uint32x4_t) vzip2q_u16(rcol_x, vdupq_n_u16(0));

            tr_v_lo = (uint32x4_t) vzip1q_u16(tr_v, vdupq_n_u16(0));
            tr_v_hi = (uint32x4_t) vzip2q_u16(tr_v, vdupq_n_u16(0));

            /*FIXME we could shift from diff from max
              and then use max to shift right to use smaller
              vectors */
            src_v_lo = vshlq_n_u32(rc_v_lo, h_scale);
            src_v_hi = vshlq_n_u32(rc_v_hi, h_scale);

            str_v_lo = vshlq_n_u32(tr_v_lo, w_scale);
            str_v_hi = vshlq_n_u32(tr_v_hi, w_scale);

            out_lo = vaddq_u32(src_v_lo, str_v_lo);
            out_hi = vaddq_u32(src_v_hi, str_v_hi);

            out_lo = vaddq_u32(out_lo, rnd_v);
            out_hi = vaddq_u32(out_hi, rnd_v);

            out_lo = vshrq_n_u32(out_lo, s_shift);
            out_hi = vshrq_n_u32(out_hi, s_shift);

            out = vcombine_u16(vqmovn_u32(out_lo), vqmovn_u32(out_hi));

            /* FIXME applying PDPC on the whole PU is useless
               since only max 12 pel cols from left and 12 pel
               rows from top require pdpc processing*/
            uint16x8_t xl, yt, x_v, w_x, w_y; //t_v,
            uint16x8_t pdpc_rnd, out_v;
            uint16x8_t tst;
            x_v = vld1q_u16(pdpc_w    + 8 * x);
            //t_v = _mm_loadu_si128( (src_above + 8 * x + 1));

            tst = vshlq_n_u16(out, 6);

            // x_v = _mm_unpacklo_epi8(x_v, _mm_setzero_si128());

            w_x = vmulq_u16(out, y_v);
            w_y = vmulq_u16(out, x_v);

            tst = vsubq_u16(tst, w_x);
            tst = vsubq_u16(tst, w_y);
            tst = vaddq_u16(tst, add_v);

            xl = vmulq_u16(x_v, l_v);
            yt = vmulq_u16(y_v, src_a);

            pdpc_rnd = vaddq_u16(xl, yt);

            tst = vaddq_u16(tst, pdpc_rnd);

            out_v = vshrq_n_u16(tst,6);

            out_v = vminq_u16(out_v, vdupq_n_u16(1023));
            out_v = vmaxq_u16(out_v, vdupq_n_u16(0));

            vst1q_u16((_dst + 8 * x), out_v);
        }
        _dst += dst_stride;
    }
}
void rcn_init_dc_planar_functions_neon(struct RCNFunctions *const rcn_funcs){
  rcn_funcs->dc.pdpc = &vvc_intra_dc_pdpc_neon;
  rcn_funcs->planar.pdpc[0] = &vvc_intra_planar_pdpc_neon;
  rcn_funcs->planar.pdpc[1] = &vvc_intra_planar_pdpc_2_neon;
}
