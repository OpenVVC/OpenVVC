#include "emmintrin.h"
#include "smmintrin.h"
#include "x86/rcn_intra_dc_planar_sse.h"
#include "ovutils.h"
#include "stdint.h"

static const uint8_t vvc_pdpc_w[3][128] = {
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
vvc_intra_dc_pdpc_sse(const uint16_t *const src_above,
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
    const uint8_t *pdpc_w = vvc_pdpc_w[pdpc_scale];
    #else
    const uint8_t *pdpc_w = vvc_pdpc_w[pdpc_scale];
    #endif
    int nb_sse_vec8_w = (1 << log2_pb_w) + ((1 << 3) - 1) >> 3;
    int nb_sse_vec8_h = (1 << log2_pb_h) + ((1 << 3) - 1) >> 3;

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

    __m128i dc_v  = _mm_set1_epi16(dc_val);
    __m128i add_v = _mm_set1_epi16(32);

    for (int y = 0; y < (1 << log2_pb_h); ++y){
        int x;
        __m128i l_v, y_v;
        #if 1
        int y_wgh = pdpc_sh[y];
        #else
        int y_wgh = pdpc_w[y];
        #endif
        const int16_t l_val = src_left[y + 1];
        l_v  = _mm_set1_epi16(l_val);
        y_v  = _mm_set1_epi16(y_wgh);
        /*FIXME PDPC weights are derived from log2 values all
           multiplications could be replaced by shift operations
           It will be easier to derive shift when specializing according
           to sizes*/
        for (x = 0; x < nb_sse_vec8_w; x++){
            __m128i xl, yt, x_v, t_v, w_x, w_y;
            __m128i pdpc_rnd, out_v;
            __m128i tst;
            x_v = _mm_loadu_si128((__m128i *) (pdpc_w    + 8 * x));
            t_v = _mm_loadu_si128((__m128i *) (src_above + 8 * x + 1));

            tst = _mm_slli_epi16(dc_v, 6);

            x_v = _mm_unpacklo_epi8(x_v, _mm_setzero_si128());

            #if 0
            w_x = _mm_mullo_epi16(dc_v, y_v);
            w_y = _mm_mullo_epi16(dc_v, x_v);
            #else
            w_x = _mm_slli_epi16(dc_v, y_wgh);
            w_y = _mm_mullo_epi16(dc_v, x_v);
            #endif

            tst = _mm_subs_epu16(tst, w_x);
            tst = _mm_subs_epu16(tst, w_y);
            tst = _mm_adds_epu16(tst, add_v);

            #if 0
            xl = _mm_mullo_epi16(l_v, x_v);
            yt = _mm_mullo_epi16(t_v, y_v);
            #else
            xl = _mm_mullo_epi16(l_v, x_v);
            yt = _mm_slli_epi16(t_v, y_wgh);
            #endif

            pdpc_rnd = _mm_adds_epu16(xl, yt);

            tst = _mm_adds_epu16(tst, pdpc_rnd);

            out_v = _mm_srli_epi16(tst,6);

            out_v = _mm_min_epi16(out_v, _mm_set1_epi16(1023));
            out_v = _mm_max_epi16(out_v, _mm_set1_epi16(0));

            _mm_storeu_si128((__m128i *)(_dst + 8 * x), out_v);
        }
        _dst += dst_stride;
    }
}

#define CUT_PDPC 0
void
vvc_intra_planar_pdpc_sse(const uint16_t *const src_above,
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

    int16_t t_row[128], b_row[128], r_col[128], l_col[128];
    const int8_t *pdpc_w = vvc_pdpc_w[pdpc_scale];
    const int16_t bl_val = src_left[height + 1];
    const int16_t tr_val = src_above[width + 1];
    #if CUT_PDPC
    const int pdpc_stop_w = OVMIN(3 << pdpc_scale, width);
    const int pdpc_stop_h = OVMIN(3 << pdpc_scale, height);
    #endif
    int x,y;
    int nb_sse_vec8_w = (1 << log2_pb_w) + ((1 << 3) - 1) >> 3;
    int nb_sse_vec8_h = (1 << log2_pb_h) + ((1 << 3) - 1) >> 3;

    __m128i tr_val_v = _mm_set1_epi16(src_above[width + 1]);
    __m128i bl_val_v = _mm_set1_epi16(src_left[height + 1]);

    for(x = 0; x < nb_sse_vec8_w; ++x){
        __m128i tr_v, br_v;
        __m128i src_a = _mm_loadu_si128((__m128i *)(src_above + x * 8 + 1));
        //br_v = _mm_sub_epi16(bl_val_v, src_a);
        tr_v = _mm_slli_epi16(src_a, h_scale);
        //_mm_storeu_si128((__m128i*) &b_row[y*8], br_v);
        _mm_storeu_si128((__m128i*) &t_row[x*8], tr_v);
    }
    for(y = 0; y < nb_sse_vec8_h; ++y){
        __m128i lc_v, rc_v;
        __m128i src_l = _mm_loadu_si128((__m128i *)(src_left + y * 8 + 1));
        rc_v = _mm_sub_epi16(tr_val_v, src_l);
        lc_v = _mm_slli_epi16(src_l, w_scale);
        _mm_storeu_si128((__m128i*) &r_col[y*8], rc_v);
        _mm_storeu_si128((__m128i*) &l_col[y*8], lc_v);
    }
    #if 1
    int rc_scale = h_scale > w_scale ? h_scale - w_scale : 0;
    int tr_scale = w_scale > h_scale ? w_scale - h_scale : 0;
    int max_scale = OVMAX(h_scale,w_scale);
    __m128i rnd_v = _mm_set1_epi16(1 << max_scale);
    #else
    __m128i rnd_v = _mm_set1_epi32(offset);
    #endif
    __m128i mul_col = _mm_set_epi16(8, 7, 6, 5, 4, 3, 2, 1);
    __m128i add_v = _mm_set1_epi16(32);
    for (y = 0; y < height; ++y){
        __m128i rcol_v = _mm_set1_epi16(r_col[y]);
        __m128i lcol_v = _mm_set1_epi16(l_col[y]);
        __m128i rcol_mul_v = _mm_mullo_epi16(rcol_v, mul_col);
        __m128i rcol_x ;
        /* pdpc var */
        int y_wgh = pdpc_w[y];
        int32_t l_val = (int32_t)src_left[y + 1];
        __m128i l_v = _mm_set1_epi16(l_val);
        __m128i y_v = _mm_set1_epi16(y_wgh);
        for (x = 0; x < nb_sse_vec8_w; ++x) {
            __m128i out_lo, out_hi;
            __m128i rc_v_lo, rc_v_hi;
            __m128i tr_v_lo, tr_v_hi;
            __m128i src_v_lo, src_v_hi;
            __m128i str_v_lo, str_v_hi;

            __m128i src_a = _mm_loadu_si128((__m128i *)(src_above + x * 8 + 1));
            //__m128i br_v = _mm_loadu_si128((__m128i *)(b_row + x * 8));
            __m128i tr_v = _mm_loadu_si128((__m128i *)(t_row + x * 8));

            tr_v = _mm_sub_epi16(tr_v, src_a);
            tr_v = _mm_add_epi16(tr_v, bl_val_v);

            _mm_storeu_si128((__m128i*) &t_row[x*8], tr_v);

            rcol_x = _mm_add_epi16(lcol_v, rcol_mul_v);
            rcol_mul_v = _mm_add_epi16(rcol_mul_v, _mm_slli_epi16(rcol_v, 3));

            #if 0
            rc_v_lo = _mm_unpacklo_epi16(rcol_x, _mm_setzero_si128());
            rc_v_hi = _mm_unpackhi_epi16(rcol_x, _mm_setzero_si128());

            tr_v_lo = _mm_unpacklo_epi16(tr_v, _mm_setzero_si128());
            tr_v_hi = _mm_unpackhi_epi16(tr_v, _mm_setzero_si128());

            /*FIXME we could shift from diff from max
              and then use max to shift right to use smaller
              vectors */
            src_v_lo = _mm_slli_epi32(rc_v_lo, h_scale);
            src_v_hi = _mm_slli_epi32(rc_v_hi, h_scale);

            str_v_lo = _mm_slli_epi32(tr_v_lo, w_scale);
            str_v_hi = _mm_slli_epi32(tr_v_hi, w_scale);

            out_lo = _mm_add_epi32(src_v_lo, str_v_lo);
            out_hi = _mm_add_epi32(src_v_hi, str_v_hi);

            out_lo = _mm_add_epi32(out_lo, rnd_v);
            out_hi = _mm_add_epi32(out_hi, rnd_v);

            out_lo = _mm_srli_epi32(out_lo, s_shift);
            out_hi = _mm_srli_epi32(out_hi, s_shift);

            out_lo = _mm_packs_epi32(out_lo, out_hi);
            #else
            rc_v_lo = rcol_x;

            tr_v_lo = tr_v;

            src_v_lo = _mm_slli_epi16(rc_v_lo, rc_scale);
            str_v_lo = _mm_slli_epi16(tr_v_lo, tr_scale);

            out_lo = _mm_add_epi16(src_v_lo, str_v_lo);

            out_lo = _mm_add_epi16(out_lo, rnd_v);

            out_lo = _mm_srli_epi16(out_lo, max_scale + 1);

            #endif

            //_mm_storeu_si128((__m128i *) (&tmp [8 * x]), tmp_lo);
        #if 0
        }
        for (x = 0; x < nb_sse_vec8_w; ++x) {
        #endif
            /* FIXME applying PDPC on the whole PU is useless
               since only max 12 pel cols from left and 12 pel
               rows from top require pdpc processing*/
            __m128i xl, yt, x_v, t_v, w_x, w_y;
            __m128i pdpc_rnd, out_v;
            __m128i tst;
            x_v = _mm_loadu_si128((__m128i *) (pdpc_w    + 8 * x));
            //t_v = _mm_loadu_si128((__m128i *) (src_above + 8 * x + 1));

            tst = _mm_slli_epi16(out_lo, 6);

            x_v = _mm_unpacklo_epi8(x_v, _mm_setzero_si128());

            w_x = _mm_mullo_epi16(out_lo, y_v);
            w_y = _mm_mullo_epi16(out_lo, x_v);

            tst = _mm_subs_epu16(tst, w_x);
            tst = _mm_subs_epu16(tst, w_y);
            tst = _mm_adds_epu16(tst, add_v);

            xl = _mm_mullo_epi16(x_v, l_v);
            yt = _mm_mullo_epi16(y_v, src_a);

            pdpc_rnd = _mm_adds_epu16(xl, yt);

            tst = _mm_adds_epu16(tst, pdpc_rnd);

            out_v = _mm_srli_epi16(tst,6);

            out_v = _mm_min_epi16(out_v, _mm_set1_epi16(1023));
            out_v = _mm_max_epi16(out_v, _mm_set1_epi16(0));

            _mm_storeu_si128((__m128i *)(_dst + 8 * x), out_v);
        }
        _dst += dst_stride;
    }
}

void
vvc_intra_planar_pdpc_2_sse(const uint16_t *const src_above,
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

    int16_t t_row[128], b_row[128], r_col[128], l_col[128];
    const int8_t *pdpc_w = vvc_pdpc_w[pdpc_scale];
    const int16_t bl_val = src_left[height + 1];
    const int16_t tr_val = src_above[width + 1];
    #if CUT_PDPC
    const int pdpc_stop_w = OVMIN(3 << pdpc_scale, width);
    const int pdpc_stop_h = OVMIN(3 << pdpc_scale, height);
    #endif
    int x,y;
    int nb_sse_vec8_w = (1 << log2_pb_w) + ((1 << 3) - 1) >> 3;
    int nb_sse_vec8_h = (1 << log2_pb_h) + ((1 << 3) - 1) >> 3;

    __m128i tr_val_v = _mm_set1_epi16(src_above[width + 1]);
    __m128i bl_val_v = _mm_set1_epi16(src_left[height + 1]);

    for(x = 0; x < nb_sse_vec8_w; ++x){
        __m128i tr_v, br_v;
        __m128i src_a = _mm_loadu_si128((__m128i *)(src_above + x * 8 + 1));
        //br_v = _mm_sub_epi16(bl_val_v, src_a);
        tr_v = _mm_slli_epi16(src_a, h_scale);
        //_mm_storeu_si128((__m128i*) &b_row[y*8], br_v);
        _mm_storeu_si128((__m128i*) &t_row[x*8], tr_v);
    }
    for(y = 0; y < nb_sse_vec8_h; ++y){
        __m128i lc_v, rc_v;
        __m128i src_l = _mm_loadu_si128((__m128i *)(src_left + y * 8 + 1));
        rc_v = _mm_sub_epi16(tr_val_v, src_l);
        lc_v = _mm_slli_epi16(src_l, w_scale);
        _mm_storeu_si128((__m128i*) &r_col[y*8], rc_v);
        _mm_storeu_si128((__m128i*) &l_col[y*8], lc_v);
    }
    #if 0
    int rc_scale = h_scale > w_scale ? h_scale - w_scale : 0;
    int tr_scale = w_scale > h_scale ? w_scale - h_scale : 0;
    int max_scale = OVMAX(h_scale,w_scale);
    __m128i rnd_v = _mm_set1_epi16(1 << max_scale);
    #else
    __m128i rnd_v = _mm_set1_epi32(offset);
    #endif
    __m128i mul_col = _mm_set_epi16(8, 7, 6, 5, 4, 3, 2, 1);
    __m128i add_v = _mm_set1_epi16(32);
    for (y = 0; y < height; ++y){
        __m128i rcol_v = _mm_set1_epi16(r_col[y]);
        __m128i lcol_v = _mm_set1_epi16(l_col[y]);
        __m128i rcol_mul_v = _mm_mullo_epi16(rcol_v, mul_col);
        __m128i rcol_x ;
        /* pdpc var */
        int y_wgh = pdpc_w[y];
        int32_t l_val = (int32_t)src_left[y + 1];
        __m128i l_v = _mm_set1_epi16(l_val);
        __m128i y_v = _mm_set1_epi16(y_wgh);
        for (x = 0; x < nb_sse_vec8_w; ++x) {
            __m128i out_lo, out_hi;
            __m128i rc_v_lo, rc_v_hi;
            __m128i tr_v_lo, tr_v_hi;
            __m128i src_v_lo, src_v_hi;
            __m128i str_v_lo, str_v_hi;

            __m128i src_a = _mm_loadu_si128((__m128i *)(src_above + x * 8 + 1));
            //__m128i br_v = _mm_loadu_si128((__m128i *)(b_row + x * 8));
            __m128i tr_v = _mm_loadu_si128((__m128i *)(t_row + x * 8));

            tr_v = _mm_sub_epi16(tr_v, src_a);
            tr_v = _mm_add_epi16(tr_v, bl_val_v);

            _mm_storeu_si128((__m128i*) &t_row[x*8], tr_v);

            rcol_x = _mm_add_epi16(lcol_v, rcol_mul_v);
            rcol_mul_v = _mm_add_epi16(rcol_mul_v, _mm_slli_epi16(rcol_v, 3));

            #if 1
            rc_v_lo = _mm_unpacklo_epi16(rcol_x, _mm_setzero_si128());
            rc_v_hi = _mm_unpackhi_epi16(rcol_x, _mm_setzero_si128());

            tr_v_lo = _mm_unpacklo_epi16(tr_v, _mm_setzero_si128());
            tr_v_hi = _mm_unpackhi_epi16(tr_v, _mm_setzero_si128());

            /*FIXME we could shift from diff from max
              and then use max to shift right to use smaller
              vectors */
            src_v_lo = _mm_slli_epi32(rc_v_lo, h_scale);
            src_v_hi = _mm_slli_epi32(rc_v_hi, h_scale);

            str_v_lo = _mm_slli_epi32(tr_v_lo, w_scale);
            str_v_hi = _mm_slli_epi32(tr_v_hi, w_scale);

            out_lo = _mm_add_epi32(src_v_lo, str_v_lo);
            out_hi = _mm_add_epi32(src_v_hi, str_v_hi);

            out_lo = _mm_add_epi32(out_lo, rnd_v);
            out_hi = _mm_add_epi32(out_hi, rnd_v);

            out_lo = _mm_srli_epi32(out_lo, s_shift);
            out_hi = _mm_srli_epi32(out_hi, s_shift);

            out_lo = _mm_packs_epi32(out_lo, out_hi);
            #else
            rc_v_lo = rcol_x;

            tr_v_lo = tr_v;

            src_v_lo = _mm_slli_epi16(rc_v_lo, rc_scale);
            str_v_lo = _mm_slli_epi16(tr_v_lo, tr_scale);

            out_lo = _mm_add_epi16(src_v_lo, str_v_lo);

            out_lo = _mm_add_epi16(out_lo, rnd_v);

            out_lo = _mm_srli_epi16(out_lo, max_scale + 1);

            #endif

            //_mm_storeu_si128((__m128i *) (&tmp [8 * x]), tmp_lo);
        #if 0
        }
        for (x = 0; x < nb_sse_vec8_w; ++x) {
        #endif
            /* FIXME applying PDPC on the whole PU is useless
               since only max 12 pel cols from left and 12 pel
               rows from top require pdpc processing*/
            __m128i xl, yt, x_v, t_v, w_x, w_y;
            __m128i pdpc_rnd, out_v;
            __m128i tst;
            x_v = _mm_loadu_si128((__m128i *) (pdpc_w    + 8 * x));
            //t_v = _mm_loadu_si128((__m128i *) (src_above + 8 * x + 1));

            tst = _mm_slli_epi16(out_lo, 6);

            x_v = _mm_unpacklo_epi8(x_v, _mm_setzero_si128());

            w_x = _mm_mullo_epi16(out_lo, y_v);
            w_y = _mm_mullo_epi16(out_lo, x_v);

            tst = _mm_subs_epu16(tst, w_x);
            tst = _mm_subs_epu16(tst, w_y);
            tst = _mm_adds_epu16(tst, add_v);

            xl = _mm_mullo_epi16(x_v, l_v);
            yt = _mm_mullo_epi16(y_v, src_a);

            pdpc_rnd = _mm_adds_epu16(xl, yt);

            tst = _mm_adds_epu16(tst, pdpc_rnd);

            out_v = _mm_srli_epi16(tst,6);

            out_v = _mm_min_epi16(out_v, _mm_set1_epi16(1023));
            out_v = _mm_max_epi16(out_v, _mm_set1_epi16(0));

            _mm_storeu_si128((__m128i *)(_dst + 8 * x), out_v);
        }
        _dst += dst_stride;
    }
}

void rcn_init_dc_planar_functions_sse(struct RCNFunctions *const rcn_funcs){
  rcn_funcs->dc.pdpc = &vvc_intra_dc_pdpc_sse;

  rcn_funcs->planar.pdpc[0] = &vvc_intra_planar_pdpc_sse;
  rcn_funcs->planar.pdpc[1] = &vvc_intra_planar_pdpc_2_sse;
}
