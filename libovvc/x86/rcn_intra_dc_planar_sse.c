#include "emmintrin.h"
#include "smmintrin.h"
#include "ovutils.h"
#include "rcn_structures.h"
#include "stdint.h"
#include "x86/vvc_utils_sse.h"

static const uint8_t vvc_pdpc_w[3][16] = {
{32,  8,  2,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{32, 16,  8,  4, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{32, 32, 16, 16, 8, 8, 4, 4, 2, 2, 1, 1, 0, 0, 0, 0},
};

static const uint8_t vvc_pdpc_w_sh[3][16] = {
{5,  3,  1,  0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F},
{5,  4,  3,     2,    1,    0, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F, 0x2F},
{5,  5,  4,     4,    3,    3,    2,    2,    1,    1,    0,    0, 0x2F, 0x2F, 0x2F, 0x2F},
};

#define VECTOR_DEFINITION_SCALE0()  const uint8_t *pdpc_sh = vvc_pdpc_w_sh[0];           \
                                    const uint8_t *pdpc_w = vvc_pdpc_w[0];               \
                                    __m128i x_v = _mm_loadl_epi64((__m128i *) (pdpc_w)); \
                                    x_v = _mm_unpacklo_epi8(x_v, _mm_setzero_si128());   \
                                    __m128i dc_v  = _mm_set1_epi16(dc_val);              \
                                    __m128i add_v = _mm_set1_epi16(32);                  \
                                    __m128i w_y = _mm_mullo_epi16(dc_v, x_v);            \
                                    __m128i tst = _mm_slli_epi16(dc_v, 6);               \
                                    tst = _mm_adds_epu16(tst, add_v);                    \
                                    __m128i tst0 = _mm_subs_epu16(tst, w_y);

#define VECTOR_DEFINITION_SCALE1()  const uint8_t *pdpc_sh = vvc_pdpc_w_sh[1];           \
                                    const uint8_t *pdpc_w = vvc_pdpc_w[1];               \
                                    __m128i x_v = _mm_loadu_si128((__m128i *) (pdpc_w)); \
                                    x_v = _mm_unpacklo_epi8(x_v, _mm_setzero_si128());   \
                                    __m128i dc_v  = _mm_set1_epi16(dc_val);              \
                                    __m128i add_v = _mm_set1_epi16(32);                  \
                                    __m128i w_y = _mm_mullo_epi16(dc_v, x_v);            \
                                    __m128i tst = _mm_slli_epi16(dc_v, 6);               \
                                    tst = _mm_adds_epu16(tst, add_v);                    \
                                    __m128i tst0 = _mm_subs_epu16(tst, w_y);

#define VECTOR_DEFINITION_SCALE2()  const uint8_t *pdpc_sh = vvc_pdpc_w_sh[2];           \
                                    const uint8_t *pdpc_w = vvc_pdpc_w[2];               \
                                    __m128i x_v = _mm_loadu_si128((__m128i *) (pdpc_w)); \
                                    __m128i x_v1, x_v2;                                  \
                                    x_v1 = _mm_unpacklo_epi8(x_v, _mm_setzero_si128());  \
                                    x_v2 = _mm_unpackhi_epi8(x_v, _mm_setzero_si128());  \
                                    __m128i dc_v  = _mm_set1_epi16(dc_val);              \
                                    __m128i add_v = _mm_set1_epi16(32);                  \
                                    __m128i w_y1 = _mm_mullo_epi16(dc_v, x_v1);          \
                                    __m128i w_y2 = _mm_mullo_epi16(dc_v, x_v2);          \
                                    __m128i tst = _mm_slli_epi16(dc_v, 6);               \
                                    tst = _mm_adds_epu16(tst, add_v);                    \
                                    __m128i tst01 = _mm_subs_epu16(tst, w_y1);           \
                                    __m128i tst02 = _mm_subs_epu16(tst, w_y2);

#define LOOP4_PDPC()    __m128i l_v;                                        \
                        int y_wgh = pdpc_sh[y];                             \
                                                                            \
                        const int16_t l_val = src_left[y + 1];              \
                        l_v  = _mm_set1_epi16(l_val);                       \
                                                                            \
                        __m128i xl, yt, w_x;                                \
                        __m128i pdpc_rnd, out_v;                            \
                        __m128i tst;                                        \
                                                                            \
                        w_x = _mm_slli_epi16(dc_v, y_wgh);                  \
                                                                            \
                        tst = _mm_subs_epu16(tst0, w_x);                    \
                                                                            \
                        xl = _mm_mullo_epi16(l_v, x_v);                     \
                        yt = _mm_slli_epi16(a1, y_wgh);                     \
                                                                            \
                        pdpc_rnd = _mm_adds_epu16(xl, yt);                  \
                                                                            \
                        tst = _mm_adds_epu16(tst, pdpc_rnd);                \
                                                                            \
                        out_v = _mm_srli_epi16(tst,6);                      \
                                                                            \
                        out_v = _mm_min_epi16(out_v, _mm_set1_epi16(1023)); \
                        out_v = _mm_max_epi16(out_v, _mm_set1_epi16(0));    \
                                                                            \
                        _mm_storel_epi64((__m128i *)(_dst), out_v);         \
                                                                            \
                        _dst += dst_stride;

#define LOOP4_NOPDPC()  __m128i l_v;                                        \
                        const int16_t l_val = src_left[y + 1];              \
                        l_v  = _mm_set1_epi16(l_val);                       \
                                                                            \
                        __m128i xl;                                         \
                        __m128i out_v;                                      \
                        __m128i tst;                                        \
                                                                            \
                        xl = _mm_mullo_epi16(l_v, x_v);                     \
                                                                            \
                        tst = _mm_adds_epu16(tst0, xl);                     \
                                                                            \
                        out_v = _mm_srli_epi16(tst,6);                      \
                                                                            \
                        out_v = _mm_min_epi16(out_v, _mm_set1_epi16(1023)); \
                        out_v = _mm_max_epi16(out_v, _mm_set1_epi16(0));    \
                                                                            \
                        _mm_storel_epi64((__m128i *)(_dst), out_v);         \
                                                                            \
                        _dst += dst_stride;

#define LOOP8_PDPC()    __m128i l_v;                                        \
                        int y_wgh = pdpc_sh[y];                             \
                                                                            \
                        const int16_t l_val = src_left[y + 1];              \
                        l_v  = _mm_set1_epi16(l_val);                       \
                                                                            \
                        __m128i xl, yt, w_x;                                \
                        __m128i pdpc_rnd, out_v;                            \
                        __m128i tst;                                        \
                                                                            \
                        w_x = _mm_slli_epi16(dc_v, y_wgh);                  \
                                                                            \
                        tst = _mm_subs_epu16(tst0, w_x);                    \
                                                                            \
                        xl = _mm_mullo_epi16(l_v, x_v);                     \
                        yt = _mm_slli_epi16(a1, y_wgh);                     \
                                                                            \
                        pdpc_rnd = _mm_adds_epu16(xl, yt);                  \
                                                                            \
                        tst = _mm_adds_epu16(tst, pdpc_rnd);                \
                                                                            \
                        out_v = _mm_srli_epi16(tst,6);                      \
                                                                            \
                        out_v = _mm_min_epi16(out_v, _mm_set1_epi16(1023)); \
                        out_v = _mm_max_epi16(out_v, _mm_set1_epi16(0));    \
                                                                            \
                        _mm_storeu_si128((__m128i *)(_dst), out_v);         \
                                                                            \
                        _dst += dst_stride;

#define LOOP8_NOPDPC()  __m128i l_v;                                        \
                        const int16_t l_val = src_left[y + 1];              \
                        l_v  = _mm_set1_epi16(l_val);                       \
                                                                            \
                        __m128i xl;                                         \
                        __m128i out_v;                                      \
                        __m128i tst;                                        \
                                                                            \
                        xl = _mm_mullo_epi16(l_v, x_v);                     \
                                                                            \
                        tst = _mm_adds_epu16(tst0, xl);                     \
                                                                            \
                        out_v = _mm_srli_epi16(tst,6);                      \
                                                                            \
                        out_v = _mm_min_epi16(out_v, _mm_set1_epi16(1023)); \
                        out_v = _mm_max_epi16(out_v, _mm_set1_epi16(0));    \
                                                                            \
                        _mm_storeu_si128((__m128i *)(_dst), out_v);         \
                                                                            \
                        _dst += dst_stride;

#define LOOP16_PDPC()   __m128i l_v;                                            \
                                                                                \
                        int y_wgh = pdpc_sh[y];                                 \
                                                                                \
                        const int16_t l_val = src_left[y + 1];                  \
                        l_v  = _mm_set1_epi16(l_val);                           \
                                                                                \
                        __m128i xl1, yt1, yt2, w_x;                             \
                        __m128i pdpc_rnd1;                                      \
                        __m128i out_v1, out_v2;                                 \
                        __m128i tst1, tst2;                                     \
                                                                                \
                        w_x = _mm_slli_epi16(dc_v, y_wgh);                      \
                                                                                \
                        tst1 = _mm_subs_epu16(tst0, w_x);                       \
                        tst2 = _mm_subs_epu16(tst, w_x);                        \
                                                                                \
                        xl1 = _mm_mullo_epi16(l_v, x_v);                        \
                        yt1 = _mm_slli_epi16(a1, y_wgh);                        \
                        yt2 = _mm_slli_epi16(a2, y_wgh);                        \
                                                                                \
                        pdpc_rnd1 = _mm_adds_epu16(xl1, yt1);                   \
                                                                                \
                        tst1 = _mm_adds_epu16(tst1, pdpc_rnd1);                 \
                        tst2 = _mm_adds_epu16(tst2, yt2);                       \
                                                                                \
                        out_v1 = _mm_srli_epi16(tst1, 6);                       \
                        out_v2 = _mm_srli_epi16(tst2, 6);                       \
                                                                                \
                        out_v1 = _mm_min_epi16(out_v1, _mm_set1_epi16(1023));   \
                        out_v1 = _mm_max_epi16(out_v1, _mm_set1_epi16(0));      \
                                                                                \
                        out_v2 = _mm_min_epi16(out_v2, _mm_set1_epi16(1023));   \
                        out_v2 = _mm_max_epi16(out_v2, _mm_set1_epi16(0));      \
                                                                                \
                        _mm_storeu_si128((__m128i *)(_dst), out_v1);            \
                        _mm_storeu_si128((__m128i *)(_dst + 8), out_v2);        \
                                                                                \
                        _dst += dst_stride;

#define LOOP16_PDPC2()  __m128i l_v;                                                    \
                                                                                        \
                        int y_wgh = pdpc_sh[y];                                         \
                                                                                        \
                        const int16_t l_val = src_left[y + 1];                          \
                        l_v  = _mm_set1_epi16(l_val);                                   \
                                                                                        \
                        __m128i xl1, xl2, yt1, yt2, w_x;                                \
                        __m128i pdpc_rnd1, pdpc_rnd2, out_v1, out_v2;                   \
                        __m128i tst1, tst2;                                             \
                                                                                        \
                        w_x = _mm_slli_epi16(dc_v, y_wgh);                              \
                                                                                        \
                        tst1 = _mm_subs_epu16(tst01, w_x);                              \
                        tst2 = _mm_subs_epu16(tst02, w_x);                              \
                                                                                        \
                        xl1 = _mm_mullo_epi16(l_v, x_v1);                               \
                        xl2 = _mm_mullo_epi16(l_v, x_v2);                               \
                        yt1 = _mm_slli_epi16(a1, y_wgh);                                \
                        yt2 = _mm_slli_epi16(a2, y_wgh);                                \
                                                                                        \
                        pdpc_rnd1 = _mm_adds_epu16(xl1, yt1);                           \
                        pdpc_rnd2 = _mm_adds_epu16(xl2, yt2);                           \
                                                                                        \
                        tst1 = _mm_adds_epu16(tst1, pdpc_rnd1);                         \
                        tst2 = _mm_adds_epu16(tst2, pdpc_rnd2);                         \
                                                                                        \
                        out_v1 = _mm_srli_epi16(tst1, 6);                               \
                        out_v2 = _mm_srli_epi16(tst2, 6);                               \
                                                                                        \
                        out_v1 = _mm_min_epi16(out_v1, _mm_set1_epi16(1023));           \
                        out_v2 = _mm_min_epi16(out_v2, _mm_set1_epi16(1023));           \
                        out_v1 = _mm_max_epi16(out_v1, _mm_set1_epi16(0));              \
                        out_v2 = _mm_max_epi16(out_v2, _mm_set1_epi16(0));              \
                                                                                        \
                        _mm_storeu_si128((__m128i *)(_dst), out_v1);                    \
                        _mm_storeu_si128((__m128i *)(_dst + 8), out_v2);                \
                                                                                        \
                        _dst += dst_stride;

#define LOOP16_NOPDPC() __m128i l_v;                                            \
                                                                                \
                        const int16_t l_val = src_left[y + 1];                  \
                        l_v  = _mm_set1_epi16(l_val);                           \
                                                                                \
                        __m128i xl1, yt1, yt2;                                  \
                        __m128i pdpc_rnd1;                                      \
                        __m128i out_v1, out_v2;                                 \
                        __m128i tst1, tst2;                                     \
                                                                                \
                        xl1 = _mm_mullo_epi16(l_v, x_v);                        \
                                                                                \
                        tst1 = _mm_adds_epu16(tst0, xl1);                       \
                                                                                \
                        out_v1 = _mm_srli_epi16(tst1, 6);                       \
                        out_v2 = _mm_srli_epi16(tst, 6);                        \
                                                                                \
                        out_v1 = _mm_min_epi16(out_v1, _mm_set1_epi16(1023));   \
                        out_v1 = _mm_max_epi16(out_v1, _mm_set1_epi16(0));      \
                                                                                \
                        out_v2 = _mm_min_epi16(out_v2, _mm_set1_epi16(1023));   \
                        out_v2 = _mm_max_epi16(out_v2, _mm_set1_epi16(0));      \
                                                                                \
                        _mm_storeu_si128((__m128i *)(_dst), out_v1);            \
                        _mm_storeu_si128((__m128i *)(_dst + 8), out_v2);        \
                                                                                \
                        _dst += dst_stride;

#define LOOP16_NOPDPC2()    __m128i l_v;                                           \
                            const int16_t l_val = src_left[y + 1];                 \
                            l_v  = _mm_set1_epi16(l_val);                          \
                                                                                   \
                            __m128i xl1, xl2;                                      \
                            __m128i pdpc_rnd1, pdpc_rnd2, out_v1, out_v2;          \
                            __m128i tst1, tst2;                                    \
                                                                                   \
                            xl1 = _mm_mullo_epi16(l_v, x_v1);                      \
                            xl2 = _mm_mullo_epi16(l_v, x_v2);                      \
                                                                                   \
                            tst1 = _mm_adds_epu16(tst01, xl1);                     \
                            tst2 = _mm_adds_epu16(tst02, xl2);                     \
                                                                                   \
                            out_v1 = _mm_srli_epi16(tst1, 6);                      \
                            out_v2 = _mm_srli_epi16(tst2, 6);                      \
                                                                                   \
                            out_v1 = _mm_min_epi16(out_v1, _mm_set1_epi16(1023));  \
                            out_v2 = _mm_min_epi16(out_v2, _mm_set1_epi16(1023));  \
                            out_v1 = _mm_max_epi16(out_v1, _mm_set1_epi16(0));     \
                            out_v2 = _mm_max_epi16(out_v2, _mm_set1_epi16(0));     \
                                                                                   \
                            _mm_storeu_si128((__m128i *)(_dst), out_v1);           \
                            _mm_storeu_si128((__m128i *)(_dst + 8), out_v2);       \
                                                                                   \
                            _dst += dst_stride;

#define LOOP32_PDPC()   __m128i l_v;                                            \
                                                                                \
                        int y_wgh = pdpc_sh[y];                                 \
                                                                                \
                        const int16_t l_val = src_left[y + 1];                  \
                        l_v  = _mm_set1_epi16(l_val);                           \
                                                                                \
                        __m128i xl1, yt1, yt2, yt3, yt4, w_x;                   \
                        __m128i pdpc_rnd1;                                      \
                        __m128i out_v1, out_v2, out_v3, out_v4;                 \
                        __m128i tst1, tst2, tst3, tst4;                         \
                                                                                \
                        w_x = _mm_slli_epi16(dc_v, y_wgh);                      \
                                                                                \
                        tst1 = _mm_subs_epu16(tst0, w_x);                       \
                        tst4 = _mm_subs_epu16(tst, w_x);                        \
                                                                                \
                        xl1 = _mm_mullo_epi16(l_v, x_v);                        \
                        yt1 = _mm_slli_epi16(a1, y_wgh);                        \
                        yt2 = _mm_slli_epi16(a2, y_wgh);                        \
                        yt3 = _mm_slli_epi16(a3, y_wgh);                        \
                        yt4 = _mm_slli_epi16(a4, y_wgh);                        \
                                                                                \
                        pdpc_rnd1 = _mm_adds_epu16(xl1, yt1);                   \
                                                                                \
                        tst1 = _mm_adds_epu16(tst1, pdpc_rnd1);                 \
                        tst2 = _mm_adds_epu16(tst4, yt2);                       \
                        tst3 = _mm_adds_epu16(tst4, yt3);                       \
                        tst4 = _mm_adds_epu16(tst4, yt4);                       \
                                                                                \
                        out_v1 = _mm_srli_epi16(tst1, 6);                       \
                        out_v2 = _mm_srli_epi16(tst2, 6);                       \
                        out_v3 = _mm_srli_epi16(tst3, 6);                       \
                        out_v4 = _mm_srli_epi16(tst4, 6);                       \
                                                                                \
                        out_v1 = _mm_min_epi16(out_v1, _mm_set1_epi16(1023));   \
                        out_v1 = _mm_max_epi16(out_v1, _mm_set1_epi16(0));      \
                                                                                \
                        out_v2 = _mm_min_epi16(out_v2, _mm_set1_epi16(1023));   \
                        out_v2 = _mm_max_epi16(out_v2, _mm_set1_epi16(0));      \
                                                                                \
                        out_v3 = _mm_min_epi16(out_v3, _mm_set1_epi16(1023));   \
                        out_v3 = _mm_max_epi16(out_v3, _mm_set1_epi16(0));      \
                                                                                \
                        out_v4 = _mm_min_epi16(out_v4, _mm_set1_epi16(1023));   \
                        out_v4 = _mm_max_epi16(out_v4, _mm_set1_epi16(0));      \
                                                                                \
                        _mm_storeu_si128((__m128i *)(_dst), out_v1);            \
                        _mm_storeu_si128((__m128i *)(_dst + 8), out_v2);        \
                        _mm_storeu_si128((__m128i *)(_dst + 16), out_v3);       \
                        _mm_storeu_si128((__m128i *)(_dst + 24), out_v4);       \
                                                                                \
                        _dst += dst_stride;

#define LOOP32_PDPC2()  __m128i l_v;                                                    \
                                                                                        \
                        int y_wgh = pdpc_sh[y];                                         \
                                                                                        \
                        const int16_t l_val = src_left[y + 1];                          \
                        l_v  = _mm_set1_epi16(l_val);                                   \
                                                                                        \
                        __m128i xl1, xl2, yt1, yt2, yt3, yt4, w_x;                      \
                        __m128i pdpc_rnd1, pdpc_rnd2, out_v1, out_v2, out_v3, out_v4;   \
                        __m128i tst1, tst2, tst3, tst4;                                 \
                                                                                        \
                        w_x = _mm_slli_epi16(dc_v, y_wgh);                              \
                                                                                        \
                        tst1 = _mm_subs_epu16(tst01, w_x);                              \
                        tst2 = _mm_subs_epu16(tst02, w_x);                              \
                        tst4 = _mm_subs_epu16(tst, w_x);                                \
                                                                                        \
                        xl1 = _mm_mullo_epi16(l_v, x_v1);                               \
                        xl2 = _mm_mullo_epi16(l_v, x_v2);                               \
                        yt1 = _mm_slli_epi16(a1, y_wgh);                                \
                        yt2 = _mm_slli_epi16(a2, y_wgh);                                \
                        yt3 = _mm_slli_epi16(a3, y_wgh);                                \
                        yt4 = _mm_slli_epi16(a4, y_wgh);                                \
                                                                                        \
                        pdpc_rnd1 = _mm_adds_epu16(xl1, yt1);                           \
                        pdpc_rnd2 = _mm_adds_epu16(xl2, yt2);                           \
                                                                                        \
                        tst1 = _mm_adds_epu16(tst1, pdpc_rnd1);                         \
                        tst2 = _mm_adds_epu16(tst2, pdpc_rnd2);                         \
                        tst3 = _mm_adds_epu16(tst4, yt3);                               \
                        tst4 = _mm_adds_epu16(tst4, yt4);                               \
                                                                                        \
                        out_v1 = _mm_srli_epi16(tst1, 6);                               \
                        out_v2 = _mm_srli_epi16(tst2, 6);                               \
                        out_v3 = _mm_srli_epi16(tst3, 6);                               \
                        out_v4 = _mm_srli_epi16(tst4, 6);                               \
                                                                                        \
                        out_v1 = _mm_min_epi16(out_v1, _mm_set1_epi16(1023));           \
                        out_v2 = _mm_min_epi16(out_v2, _mm_set1_epi16(1023));           \
                        out_v3 = _mm_min_epi16(out_v3, _mm_set1_epi16(1023));           \
                        out_v4 = _mm_min_epi16(out_v4, _mm_set1_epi16(1023));           \
                        out_v1 = _mm_max_epi16(out_v1, _mm_set1_epi16(0));              \
                        out_v2 = _mm_max_epi16(out_v2, _mm_set1_epi16(0));              \
                        out_v3 = _mm_max_epi16(out_v3, _mm_set1_epi16(0));              \
                        out_v4 = _mm_max_epi16(out_v4, _mm_set1_epi16(0));              \
                                                                                        \
                        _mm_storeu_si128((__m128i *)(_dst), out_v1);                    \
                        _mm_storeu_si128((__m128i *)(_dst + 8), out_v2);                \
                        _mm_storeu_si128((__m128i *)(_dst + 16), out_v3);               \
                        _mm_storeu_si128((__m128i *)(_dst + 24), out_v4);               \
                                                                                        \
                        _dst += dst_stride;

#define LOOP32_NOPDPC() __m128i l_v;                                            \
                                                                                \
                        const int16_t l_val = src_left[y + 1];                  \
                        l_v  = _mm_set1_epi16(l_val);                           \
                                                                                \
                        __m128i xl1;                                            \
                        __m128i out_v1, out_v2;                                 \
                        __m128i tst1;                                           \
                                                                                \
                        xl1 = _mm_mullo_epi16(l_v, x_v);                        \
                                                                                \
                        tst1 = _mm_adds_epu16(tst0, xl1);                       \
                                                                                \
                        out_v1 = _mm_srli_epi16(tst1, 6);                       \
                        out_v2 = _mm_srli_epi16(tst, 6);                        \
                                                                                \
                        out_v1 = _mm_min_epi16(out_v1, _mm_set1_epi16(1023));   \
                        out_v1 = _mm_max_epi16(out_v1, _mm_set1_epi16(0));      \
                                                                                \
                        out_v2 = _mm_min_epi16(out_v2, _mm_set1_epi16(1023));   \
                        out_v2 = _mm_max_epi16(out_v2, _mm_set1_epi16(0));      \
                                                                                \
                        _mm_storeu_si128((__m128i *)(_dst), out_v1);            \
                        _mm_storeu_si128((__m128i *)(_dst + 8), out_v2);        \
                        _mm_storeu_si128((__m128i *)(_dst + 16), out_v2);       \
                        _mm_storeu_si128((__m128i *)(_dst + 24), out_v2);       \
                                                                                \
                        _dst += dst_stride;

#define LOOP32_NOPDPC2()    __m128i l_v;                                           \
                            const int16_t l_val = src_left[y + 1];                 \
                            l_v  = _mm_set1_epi16(l_val);                          \
                                                                                   \
                            __m128i xl1, xl2;                                      \
                            __m128i pdpc_rnd1, pdpc_rnd2, out_v1, out_v2, out_v3;  \
                            __m128i tst1, tst2;                                    \
                                                                                   \
                            xl1 = _mm_mullo_epi16(l_v, x_v1);                      \
                            xl2 = _mm_mullo_epi16(l_v, x_v2);                      \
                                                                                   \
                            tst1 = _mm_adds_epu16(tst01, xl1);                     \
                            tst2 = _mm_adds_epu16(tst02, xl2);                     \
                                                                                   \
                            out_v1 = _mm_srli_epi16(tst1, 6);                      \
                            out_v2 = _mm_srli_epi16(tst2, 6);                      \
                            out_v3 = _mm_srli_epi16(tst, 6);                       \
                                                                                   \
                            out_v1 = _mm_min_epi16(out_v1, _mm_set1_epi16(1023));  \
                            out_v2 = _mm_min_epi16(out_v2, _mm_set1_epi16(1023));  \
                            out_v3 = _mm_min_epi16(out_v3, _mm_set1_epi16(1023));  \
                            out_v1 = _mm_max_epi16(out_v1, _mm_set1_epi16(0));     \
                            out_v2 = _mm_max_epi16(out_v2, _mm_set1_epi16(0));     \
                            out_v3 = _mm_max_epi16(out_v3, _mm_set1_epi16(0));     \
                                                                                   \
                            _mm_storeu_si128((__m128i *)(_dst), out_v1);           \
                            _mm_storeu_si128((__m128i *)(_dst + 8), out_v2);       \
                            _mm_storeu_si128((__m128i *)(_dst + 16), out_v3);      \
                            _mm_storeu_si128((__m128i *)(_dst + 24), out_v3);      \
                                                                                   \
                            _dst += dst_stride;

#define LOOP64_PDPC()   __m128i l_v;                                            \
                                                                                \
                        int y_wgh = pdpc_sh[y];                                 \
                                                                                \
                        const int16_t l_val = src_left[y + 1];                  \
                        l_v  = _mm_set1_epi16(l_val);                           \
                                                                                \
                        __m128i xl1, yt1, yt2, yt3, yt4, w_x;                   \
                        __m128i yt5, yt6, yt7, yt8;                             \
                        __m128i pdpc_rnd1;                                      \
                        __m128i out_v1, out_v2, out_v3, out_v4;                 \
                        __m128i out_v5, out_v6, out_v7, out_v8;                 \
                        __m128i tst1, tst2, tst3, tst4;                         \
                        __m128i tst5, tst6, tst7, tst8;                         \
                                                                                \
                        w_x = _mm_slli_epi16(dc_v, y_wgh);                      \
                                                                                \
                        tst1 = _mm_subs_epu16(tst0, w_x);                       \
                        tst8 = _mm_subs_epu16(tst, w_x);                        \
                                                                                \
                        xl1 = _mm_mullo_epi16(l_v, x_v);                        \
                        yt1 = _mm_slli_epi16(a1, y_wgh);                        \
                        yt2 = _mm_slli_epi16(a2, y_wgh);                        \
                        yt3 = _mm_slli_epi16(a3, y_wgh);                        \
                        yt4 = _mm_slli_epi16(a4, y_wgh);                        \
                        yt5 = _mm_slli_epi16(a5, y_wgh);                        \
                        yt6 = _mm_slli_epi16(a6, y_wgh);                        \
                        yt7 = _mm_slli_epi16(a7, y_wgh);                        \
                        yt8 = _mm_slli_epi16(a8, y_wgh);                        \
                                                                                \
                        pdpc_rnd1 = _mm_adds_epu16(xl1, yt1);                   \
                                                                                \
                        tst1 = _mm_adds_epu16(tst1, pdpc_rnd1);                 \
                        tst2 = _mm_adds_epu16(tst8, yt2);                       \
                        tst3 = _mm_adds_epu16(tst8, yt3);                       \
                        tst4 = _mm_adds_epu16(tst8, yt4);                       \
                        tst5 = _mm_adds_epu16(tst8, yt5);                       \
                        tst6 = _mm_adds_epu16(tst8, yt6);                       \
                        tst7 = _mm_adds_epu16(tst8, yt7);                       \
                        tst8 = _mm_adds_epu16(tst8, yt8);                       \
                                                                                \
                        out_v1 = _mm_srli_epi16(tst1, 6);                       \
                        out_v2 = _mm_srli_epi16(tst2, 6);                       \
                        out_v3 = _mm_srli_epi16(tst3, 6);                       \
                        out_v4 = _mm_srli_epi16(tst4, 6);                       \
                        out_v5 = _mm_srli_epi16(tst5, 6);                       \
                        out_v6 = _mm_srli_epi16(tst6, 6);                       \
                        out_v7 = _mm_srli_epi16(tst7, 6);                       \
                        out_v8 = _mm_srli_epi16(tst8, 6);                       \
                                                                                \
                        out_v1 = _mm_min_epi16(out_v1, _mm_set1_epi16(1023));   \
                        out_v1 = _mm_max_epi16(out_v1, _mm_set1_epi16(0));      \
                                                                                \
                        out_v2 = _mm_min_epi16(out_v2, _mm_set1_epi16(1023));   \
                        out_v2 = _mm_max_epi16(out_v2, _mm_set1_epi16(0));      \
                                                                                \
                        out_v3 = _mm_min_epi16(out_v3, _mm_set1_epi16(1023));   \
                        out_v3 = _mm_max_epi16(out_v3, _mm_set1_epi16(0));      \
                                                                                \
                        out_v4 = _mm_min_epi16(out_v4, _mm_set1_epi16(1023));   \
                        out_v4 = _mm_max_epi16(out_v4, _mm_set1_epi16(0));      \
                                                                                \
                        out_v5 = _mm_min_epi16(out_v5, _mm_set1_epi16(1023));   \
                        out_v5 = _mm_max_epi16(out_v5, _mm_set1_epi16(0));      \
                                                                                \
                        out_v6 = _mm_min_epi16(out_v6, _mm_set1_epi16(1023));   \
                        out_v6 = _mm_max_epi16(out_v6, _mm_set1_epi16(0));      \
                                                                                \
                        out_v7 = _mm_min_epi16(out_v7, _mm_set1_epi16(1023));   \
                        out_v7 = _mm_max_epi16(out_v7, _mm_set1_epi16(0));      \
                                                                                \
                        out_v8 = _mm_min_epi16(out_v8, _mm_set1_epi16(1023));   \
                        out_v8 = _mm_max_epi16(out_v8, _mm_set1_epi16(0));      \
                                                                                \
                        _mm_storeu_si128((__m128i *)(_dst), out_v1);            \
                        _mm_storeu_si128((__m128i *)(_dst + 8), out_v2);        \
                        _mm_storeu_si128((__m128i *)(_dst + 16), out_v3);       \
                        _mm_storeu_si128((__m128i *)(_dst + 24), out_v4);       \
                        _mm_storeu_si128((__m128i *)(_dst + 32), out_v5);       \
                        _mm_storeu_si128((__m128i *)(_dst + 40), out_v6);       \
                        _mm_storeu_si128((__m128i *)(_dst + 48), out_v7);       \
                        _mm_storeu_si128((__m128i *)(_dst + 56), out_v8);       \
                                                                                \
                        _dst += dst_stride;

#define LOOP64_PDPC2()  __m128i l_v;                                                    \
                                                                                        \
                        int y_wgh = pdpc_sh[y];                                         \
                                                                                        \
                        const int16_t l_val = src_left[y + 1];                          \
                        l_v  = _mm_set1_epi16(l_val);                                   \
                                                                                        \
                        __m128i xl1, xl2, yt1, yt2, yt3, yt4, w_x;                      \
                        __m128i yt5, yt6, yt7, yt8;                                     \
                        __m128i pdpc_rnd1, pdpc_rnd2, out_v1, out_v2, out_v3, out_v4;   \
                        __m128i out_v5, out_v6, out_v7, out_v8;                         \
                        __m128i tst1, tst2, tst3, tst4;                                 \
                        __m128i tst5, tst6, tst7, tst8;                                 \
                                                                                        \
                        w_x = _mm_slli_epi16(dc_v, y_wgh);                              \
                                                                                        \
                        tst1 = _mm_subs_epu16(tst01, w_x);                              \
                        tst2 = _mm_subs_epu16(tst02, w_x);                              \
                        tst8 = _mm_subs_epu16(tst, w_x);                                \
                                                                                        \
                        xl1 = _mm_mullo_epi16(l_v, x_v1);                               \
                        xl2 = _mm_mullo_epi16(l_v, x_v2);                               \
                        yt1 = _mm_slli_epi16(a1, y_wgh);                                \
                        yt2 = _mm_slli_epi16(a2, y_wgh);                                \
                        yt3 = _mm_slli_epi16(a3, y_wgh);                                \
                        yt4 = _mm_slli_epi16(a4, y_wgh);                                \
                        yt5 = _mm_slli_epi16(a5, y_wgh);                                \
                        yt6 = _mm_slli_epi16(a6, y_wgh);                                \
                        yt7 = _mm_slli_epi16(a7, y_wgh);                                \
                        yt8 = _mm_slli_epi16(a8, y_wgh);                                \
                                                                                        \
                        pdpc_rnd1 = _mm_adds_epu16(xl1, yt1);                           \
                        pdpc_rnd2 = _mm_adds_epu16(xl2, yt2);                           \
                                                                                        \
                        tst1 = _mm_adds_epu16(tst1, pdpc_rnd1);                         \
                        tst2 = _mm_adds_epu16(tst2, pdpc_rnd2);                         \
                        tst3 = _mm_adds_epu16(tst8, yt3);                               \
                        tst4 = _mm_adds_epu16(tst8, yt4);                               \
                        tst5 = _mm_adds_epu16(tst8, yt5);                               \
                        tst6 = _mm_adds_epu16(tst8, yt6);                               \
                        tst7 = _mm_adds_epu16(tst8, yt7);                               \
                        tst8 = _mm_adds_epu16(tst8, yt8);                               \
                                                                                        \
                        out_v1 = _mm_srli_epi16(tst1, 6);                               \
                        out_v2 = _mm_srli_epi16(tst2, 6);                               \
                        out_v3 = _mm_srli_epi16(tst3, 6);                               \
                        out_v4 = _mm_srli_epi16(tst4, 6);                               \
                        out_v5 = _mm_srli_epi16(tst5, 6);                               \
                        out_v6 = _mm_srli_epi16(tst6, 6);                               \
                        out_v7 = _mm_srli_epi16(tst7, 6);                               \
                        out_v8 = _mm_srli_epi16(tst8, 6);                               \
                                                                                        \
                        out_v1 = _mm_min_epi16(out_v1, _mm_set1_epi16(1023));           \
                        out_v1 = _mm_max_epi16(out_v1, _mm_set1_epi16(0));              \
                                                                                        \
                        out_v2 = _mm_min_epi16(out_v2, _mm_set1_epi16(1023));           \
                        out_v2 = _mm_max_epi16(out_v2, _mm_set1_epi16(0));              \
                                                                                        \
                        out_v3 = _mm_min_epi16(out_v3, _mm_set1_epi16(1023));           \
                        out_v3 = _mm_max_epi16(out_v3, _mm_set1_epi16(0));              \
                                                                                        \
                        out_v4 = _mm_min_epi16(out_v4, _mm_set1_epi16(1023));           \
                        out_v4 = _mm_max_epi16(out_v4, _mm_set1_epi16(0));              \
                                                                                        \
                        out_v5 = _mm_min_epi16(out_v5, _mm_set1_epi16(1023));           \
                        out_v5 = _mm_max_epi16(out_v5, _mm_set1_epi16(0));              \
                                                                                        \
                        out_v6 = _mm_min_epi16(out_v6, _mm_set1_epi16(1023));           \
                        out_v6 = _mm_max_epi16(out_v6, _mm_set1_epi16(0));              \
                                                                                        \
                        out_v7 = _mm_min_epi16(out_v7, _mm_set1_epi16(1023));           \
                        out_v7 = _mm_max_epi16(out_v7, _mm_set1_epi16(0));              \
                                                                                        \
                        out_v8 = _mm_min_epi16(out_v8, _mm_set1_epi16(1023));           \
                        out_v8 = _mm_max_epi16(out_v8, _mm_set1_epi16(0));              \
                                                                                        \
                        _mm_storeu_si128((__m128i *)(_dst), out_v1);                    \
                        _mm_storeu_si128((__m128i *)(_dst + 8), out_v2);                \
                        _mm_storeu_si128((__m128i *)(_dst + 16), out_v3);               \
                        _mm_storeu_si128((__m128i *)(_dst + 24), out_v4);               \
                        _mm_storeu_si128((__m128i *)(_dst + 32), out_v5);               \
                        _mm_storeu_si128((__m128i *)(_dst + 40), out_v6);               \
                        _mm_storeu_si128((__m128i *)(_dst + 48), out_v7);               \
                        _mm_storeu_si128((__m128i *)(_dst + 56), out_v8);               \
                                                                                        \
                        _dst += dst_stride;

#define LOOP64_NOPDPC() __m128i l_v;                                            \
                                                                                \
                        const int16_t l_val = src_left[y + 1];                  \
                        l_v  = _mm_set1_epi16(l_val);                           \
                                                                                \
                        __m128i xl1;                                            \
                        __m128i out_v1, out_v2;                                 \
                        __m128i tst1;                                           \
                                                                                \
                        xl1 = _mm_mullo_epi16(l_v, x_v);                        \
                                                                                \
                        tst1 = _mm_adds_epu16(tst0, xl1);                       \
                                                                                \
                        out_v1 = _mm_srli_epi16(tst1, 6);                       \
                        out_v2 = _mm_srli_epi16(tst, 6);                        \
                                                                                \
                        out_v1 = _mm_min_epi16(out_v1, _mm_set1_epi16(1023));   \
                        out_v1 = _mm_max_epi16(out_v1, _mm_set1_epi16(0));      \
                                                                                \
                        out_v2 = _mm_min_epi16(out_v2, _mm_set1_epi16(1023));   \
                        out_v2 = _mm_max_epi16(out_v2, _mm_set1_epi16(0));      \
                                                                                \
                        _mm_storeu_si128((__m128i *)(_dst), out_v1);            \
                        _mm_storeu_si128((__m128i *)(_dst + 8), out_v2);        \
                        _mm_storeu_si128((__m128i *)(_dst + 16), out_v2);       \
                        _mm_storeu_si128((__m128i *)(_dst + 24), out_v2);       \
                        _mm_storeu_si128((__m128i *)(_dst + 32), out_v2);       \
                        _mm_storeu_si128((__m128i *)(_dst + 40), out_v2);       \
                        _mm_storeu_si128((__m128i *)(_dst + 48), out_v2);       \
                        _mm_storeu_si128((__m128i *)(_dst + 56), out_v2);       \
                                                                                \
                        _dst += dst_stride;

#define LOOP64_NOPDPC2()    __m128i l_v;                                           \
                            const int16_t l_val = src_left[y + 1];                 \
                            l_v  = _mm_set1_epi16(l_val);                          \
                                                                                   \
                            __m128i xl1, xl2;                                      \
                            __m128i pdpc_rnd1, pdpc_rnd2, out_v1, out_v2, out_v3;  \
                            __m128i tst1, tst2;                                    \
                                                                                   \
                            xl1 = _mm_mullo_epi16(l_v, x_v1);                      \
                            xl2 = _mm_mullo_epi16(l_v, x_v2);                      \
                                                                                   \
                            tst1 = _mm_adds_epu16(tst01, xl1);                     \
                            tst2 = _mm_adds_epu16(tst02, xl2);                     \
                                                                                   \
                            out_v1 = _mm_srli_epi16(tst1, 6);                      \
                            out_v2 = _mm_srli_epi16(tst2, 6);                      \
                            out_v3 = _mm_srli_epi16(tst, 6);                       \
                                                                                   \
                            out_v1 = _mm_min_epi16(out_v1, _mm_set1_epi16(1023));  \
                            out_v2 = _mm_min_epi16(out_v2, _mm_set1_epi16(1023));  \
                            out_v3 = _mm_min_epi16(out_v3, _mm_set1_epi16(1023));  \
                            out_v1 = _mm_max_epi16(out_v1, _mm_set1_epi16(0));     \
                            out_v2 = _mm_max_epi16(out_v2, _mm_set1_epi16(0));     \
                            out_v3 = _mm_max_epi16(out_v3, _mm_set1_epi16(0));     \
                                                                                   \
                            _mm_storeu_si128((__m128i *)(_dst), out_v1);           \
                            _mm_storeu_si128((__m128i *)(_dst + 8), out_v2);       \
                            _mm_storeu_si128((__m128i *)(_dst + 16), out_v3);      \
                            _mm_storeu_si128((__m128i *)(_dst + 24), out_v3);      \
                            _mm_storeu_si128((__m128i *)(_dst + 32), out_v3);      \
                            _mm_storeu_si128((__m128i *)(_dst + 40), out_v3);      \
                            _mm_storeu_si128((__m128i *)(_dst + 48), out_v3);      \
                            _mm_storeu_si128((__m128i *)(_dst + 56), out_v3);      \
                                                                                   \
                            _dst += dst_stride;

static void
vvc_intra_dc_pdpc_w4_h4_sse(const uint16_t *const src_above,
                      const uint16_t *const src_left,
                      uint16_t *const dst, ptrdiff_t dst_stride,
                      int log2_pb_w, int log2_pb_h)
{
    int idx;
    uint16_t *_dst = dst;
    uint32_t dc_val = 0;
    const int shift  = 3;
    const int offset = 0x4;

    LOAD1(__m128i l, _mm_loadl_epi64, src_left + 1, 8)
    LOAD1(__m128i a, _mm_loadl_epi64, src_above + 1, 8)
    l1 = _mm_add_epi16(l1, a1);
    l1 = _mm_hadd_epi16(l1, l1);
    l1 = _mm_hadd_epi16(l1, l1);
    dc_val = _mm_extract_epi16(l1, 0);

    dc_val = (dc_val + offset) >> shift;

    VECTOR_DEFINITION_SCALE0();
    for (int y = 0; y < (1 << log2_pb_h); ++y){
        LOOP4_PDPC();
    }
}

static void
vvc_intra_dc_pdpc_w4_h8_sse(const uint16_t *const src_above,
                      const uint16_t *const src_left,
                      uint16_t *const dst, ptrdiff_t dst_stride,
                      int log2_pb_w, int log2_pb_h)
{
    int idx;
    uint16_t *_dst = dst;
    uint32_t dc_val = 0;
    const int shift  = 3;
    const int offset = 0x4;

    LOAD1(__m128i l, _mm_loadu_si128, src_left + 1, 8)
    LOAD1(__m128i a, _mm_loadl_epi64, src_above + 1, 8)
    l1 = _mm_hadd_epi16(l1, l1);
    l1 = _mm_hadd_epi16(l1, l1);
    l1 = _mm_hadd_epi16(l1, l1);
    dc_val = _mm_extract_epi16(l1, 0);

    dc_val = (dc_val + offset) >> shift;

    VECTOR_DEFINITION_SCALE0();
    for (int y = 0; y < 4; ++y){
        LOOP4_PDPC();
    }
    for (int y = 4; y < (1 << log2_pb_h); ++y){
        LOOP4_NOPDPC();
    }
}

static void
vvc_intra_dc_pdpc_w4_h16_sse(const uint16_t *const src_above,
                      const uint16_t *const src_left,
                      uint16_t *const dst, ptrdiff_t dst_stride,
                      int log2_pb_w, int log2_pb_h)
{
    int idx;
    uint16_t *_dst = dst;
    uint32_t dc_val = 0;
    const int shift  = 4;
    const int offset = 0x8;

    LOAD2(__m128i l, _mm_loadu_si128, src_left + 1, 8)
    LOAD1(__m128i a, _mm_loadl_epi64, src_above + 1, 8)
    l1 = _mm_add_epi16(l1, l2);
    l1 = _mm_hadd_epi16(l1, l1);
    l1 = _mm_hadd_epi16(l1, l1);
    l1 = _mm_hadd_epi16(l1, l1);
    dc_val = _mm_extract_epi16(l1, 0);

    dc_val = (dc_val + offset) >> shift;

    VECTOR_DEFINITION_SCALE1();
    for (int y = 0; y < 7; ++y){
        LOOP4_PDPC();
    }
    for (int y = 7; y < (1 << log2_pb_h); ++y){
        LOOP4_NOPDPC();
    }
}

static void
vvc_intra_dc_pdpc_w4_h32_sse(const uint16_t *const src_above,
                      const uint16_t *const src_left,
                      uint16_t *const dst, ptrdiff_t dst_stride,
                      int log2_pb_w, int log2_pb_h)
{
    int idx;
    uint16_t *_dst = dst;
    uint32_t dc_val = 0;
    const int shift  = 5;
    const int offset = 0x10;

    LOAD4(__m128i l, _mm_loadu_si128, src_left + 1, 8)
    LOAD1(__m128i a, _mm_loadl_epi64, src_above + 1, 8)
    l1 = _mm_add_epi16(l1, l2);
    l3 = _mm_add_epi16(l3, l4);
    l1 = _mm_add_epi16(l1, l3);
    l1 = _mm_hadd_epi16(l1, l1);
    l1 = _mm_hadd_epi16(l1, l1);
    l1 = _mm_hadd_epi16(l1, l1);
    dc_val = _mm_extract_epi16(l1, 0);

    dc_val = (dc_val + offset) >> shift;

    VECTOR_DEFINITION_SCALE1();
    for (int y = 0; y < 7; ++y){
        LOOP4_PDPC();
    }
    for (int y = 7; y < (1 << log2_pb_h); ++y){
        LOOP4_NOPDPC();
    }
}

static void
vvc_intra_dc_pdpc_w4_h64_sse(const uint16_t *const src_above,
                      const uint16_t *const src_left,
                      uint16_t *const dst, ptrdiff_t dst_stride,
                      int log2_pb_w, int log2_pb_h)
{
    int idx;
    uint16_t *_dst = dst;
    uint32_t dc_val = 0;
    const int shift  = 6;
    const int offset = 0x20;

    LOAD8(__m128i l, _mm_loadu_si128, src_left + 1, 8)
    LOAD1(__m128i a, _mm_loadl_epi64, src_above + 1, 8)
    l1 = _mm_add_epi16(l1, l2);
    l3 = _mm_add_epi16(l3, l4);
    l5 = _mm_add_epi16(l5, l6);
    l7 = _mm_add_epi16(l7, l8);
    l1 = _mm_add_epi16(l1, l3);
    l5 = _mm_add_epi16(l5, l7);
    l1 = _mm_add_epi16(l1, l5);
    l1 = _mm_hadd_epi16(l1, l1);
    l1 = _mm_hadd_epi16(l1, l1);
    l1 = _mm_hadd_epi16(l1, l1);
    dc_val = _mm_extract_epi16(l1, 0);

    dc_val = (dc_val + offset) >> shift;

    VECTOR_DEFINITION_SCALE1();
    for (int y = 0; y < 7; ++y){
        LOOP4_PDPC();
    }
    for (int y = 7; y < (1 << log2_pb_h); ++y){
        LOOP4_NOPDPC();
    }
}

static void
vvc_intra_dc_pdpc_w8_h4_sse(const uint16_t *const src_above,
                      const uint16_t *const src_left,
                      uint16_t *const dst, ptrdiff_t dst_stride,
                      int log2_pb_w, int log2_pb_h)
{
    int idx;
    uint16_t *_dst = dst;
    uint32_t dc_val = 0;
    const int shift  = 3;
    const int offset = 0x4;

    LOAD1(__m128i a, _mm_loadu_si128, src_above + 1, 8)
    __m128i l1 = _mm_hadd_epi16(a1, a1);
    l1 = _mm_hadd_epi16(l1, l1);
    l1 = _mm_hadd_epi16(l1, l1);
    dc_val = _mm_extract_epi16(l1, 0);

    dc_val = (dc_val + offset) >> shift;

    VECTOR_DEFINITION_SCALE0();
    for (int y = 0; y < (1 << log2_pb_h); ++y){
        LOOP8_PDPC();
    }
}

static void
vvc_intra_dc_pdpc_w8_h8_sse(const uint16_t *const src_above,
                      const uint16_t *const src_left,
                      uint16_t *const dst, ptrdiff_t dst_stride,
                      int log2_pb_w, int log2_pb_h)
{
    int idx;
    uint16_t *_dst = dst;
    uint32_t dc_val = 0;
    const int shift  = 4;
    const int offset = 0x8;

    LOAD1(__m128i l, _mm_loadu_si128, src_left + 1, 8)
    LOAD1(__m128i a, _mm_loadu_si128, src_above + 1, 8)
    l1 = _mm_add_epi16(l1, a1);
    l1 = _mm_hadd_epi16(l1, l1);
    l1 = _mm_hadd_epi16(l1, l1);
    l1 = _mm_hadd_epi16(l1, l1);
    dc_val = _mm_extract_epi16(l1, 0);

    dc_val = (dc_val + offset) >> shift;

    VECTOR_DEFINITION_SCALE1();
    for (int y = 0; y < 7; ++y){
        LOOP8_PDPC();
    }
    for (int y = 7; y < (1 << log2_pb_h); ++y){
        LOOP8_NOPDPC();
    }
}

static void
vvc_intra_dc_pdpc_w8_h16_sse(const uint16_t *const src_above,
                      const uint16_t *const src_left,
                      uint16_t *const dst, ptrdiff_t dst_stride,
                      int log2_pb_w, int log2_pb_h)
{
    int idx;
    uint16_t *_dst = dst;
    uint32_t dc_val = 0;
    const int shift  = 4;
    const int offset = 0x8;

    LOAD2(__m128i l, _mm_loadu_si128, src_left + 1, 8)
    LOAD1(__m128i a, _mm_loadu_si128, src_above + 1, 8)
    l1 = _mm_add_epi16(l1, l2);
    l1 = _mm_hadd_epi16(l1, l1);
    l1 = _mm_hadd_epi16(l1, l1);
    l1 = _mm_hadd_epi16(l1, l1);
    dc_val = _mm_extract_epi16(l1, 0);

    dc_val = (dc_val + offset) >> shift;

    VECTOR_DEFINITION_SCALE1();
    for (int y = 0; y < 7; ++y){
        LOOP8_PDPC();
    }
    for (int y = 7; y < (1 << log2_pb_h); ++y){
        LOOP8_NOPDPC();
    }
}

static void
vvc_intra_dc_pdpc_w8_h32_sse(const uint16_t *const src_above,
                      const uint16_t *const src_left,
                      uint16_t *const dst, ptrdiff_t dst_stride,
                      int log2_pb_w, int log2_pb_h)
{
    int idx;
    uint16_t *_dst = dst;
    uint32_t dc_val = 0;
    const int shift  = 5;
    const int offset = 0x10;

    LOAD4(__m128i l, _mm_loadu_si128, src_left + 1, 8)
    LOAD1(__m128i a, _mm_loadu_si128, src_above + 1, 8)
    l1 = _mm_add_epi16(l1, l2);
    l3 = _mm_add_epi16(l3, l4);
    l1 = _mm_add_epi16(l1, l3);
    l1 = _mm_hadd_epi16(l1, l1);
    l1 = _mm_hadd_epi16(l1, l1);
    l1 = _mm_hadd_epi16(l1, l1);
    dc_val = _mm_extract_epi16(l1, 0);

    dc_val = (dc_val + offset) >> shift;

    VECTOR_DEFINITION_SCALE1();
    for (int y = 0; y < 7; ++y){
        LOOP8_PDPC();
    }
    for (int y = 7; y < (1 << log2_pb_h); ++y){
        LOOP8_NOPDPC();
    }
}

static void
vvc_intra_dc_pdpc_w8_h64_sse(const uint16_t *const src_above,
                      const uint16_t *const src_left,
                      uint16_t *const dst, ptrdiff_t dst_stride,
                      int log2_pb_w, int log2_pb_h)
{
    int idx;
    uint16_t *_dst = dst;
    uint32_t dc_val = 0;
    const int shift  = 6;
    const int offset = 0x20;

    LOAD8(__m128i l, _mm_loadu_si128, src_left + 1, 8)
    LOAD1(__m128i a, _mm_loadu_si128, src_above + 1, 8)
    l1 = _mm_add_epi16(l1, l2);
    l3 = _mm_add_epi16(l3, l4);
    l5 = _mm_add_epi16(l5, l6);
    l7 = _mm_add_epi16(l7, l8);
    l1 = _mm_add_epi16(l1, l3);
    l5 = _mm_add_epi16(l5, l7);
    l1 = _mm_add_epi16(l1, l5);
    l1 = _mm_hadd_epi16(l1, l1);
    l1 = _mm_hadd_epi16(l1, l1);
    l1 = _mm_hadd_epi16(l1, l1);
    dc_val = _mm_extract_epi16(l1, 0);

    dc_val = (dc_val + offset) >> shift;

    VECTOR_DEFINITION_SCALE1();
    for (int y = 0; y < 7; ++y){
        LOOP8_PDPC();
    }
    for (int y = 7; y < (1 << log2_pb_h); ++y){
        LOOP8_NOPDPC();
    }
}

static void
vvc_intra_dc_pdpc_w16_h4_sse(const uint16_t *const src_above,
                      const uint16_t *const src_left,
                      uint16_t *const dst, ptrdiff_t dst_stride,
                      int log2_pb_w, int log2_pb_h)
{
    int idx;
    uint16_t *_dst = dst;
    uint32_t dc_val = 0;
    const int shift  = 4;
    const int offset = 0x8;


    LOAD2(__m128i a, _mm_loadu_si128, src_above + 1, 8)
    __m128i l1 = _mm_add_epi16(a1, a2);
    l1 = _mm_hadd_epi16(l1, l1);
    l1 = _mm_hadd_epi16(l1, l1);
    l1 = _mm_hadd_epi16(l1, l1);
    dc_val = _mm_extract_epi16(l1, 0);

    dc_val = (dc_val + offset) >> shift;

    VECTOR_DEFINITION_SCALE1();
    for (int y = 0; y < (1 << log2_pb_h); ++y){
        LOOP16_PDPC();
    }
}


static void
vvc_intra_dc_pdpc_w16_h8_sse(const uint16_t *const src_above,
                      const uint16_t *const src_left,
                      uint16_t *const dst, ptrdiff_t dst_stride,
                      int log2_pb_w, int log2_pb_h)
{
    int idx;
    uint16_t *_dst = dst;
    uint32_t dc_val = 0;
    const int shift  = 4;
    const int offset = 0x8;

    LOAD2(__m128i a, _mm_loadu_si128, src_above + 1, 8)
    __m128i l1 = _mm_add_epi16(a1, a2);
    l1 = _mm_hadd_epi16(l1, l1);
    l1 = _mm_hadd_epi16(l1, l1);
    l1 = _mm_hadd_epi16(l1, l1);
    dc_val = _mm_extract_epi16(l1, 0);

    dc_val = (dc_val + offset) >> shift;

    VECTOR_DEFINITION_SCALE1();
    for (int y = 0; y < 7; ++y){
        LOOP16_PDPC();
    }
    for (int y = 7; y < (1 << log2_pb_h); ++y){
        LOOP16_NOPDPC();
    }
}

static void
vvc_intra_dc_pdpc_w16_h16_sse(const uint16_t *const src_above,
                      const uint16_t *const src_left,
                      uint16_t *const dst, ptrdiff_t dst_stride,
                      int log2_pb_w, int log2_pb_h)
{
    int idx;
    uint16_t *_dst = dst;
    uint32_t dc_val = 0;
    const int shift  = 5;
    const int offset = 0x10;

    LOAD2(__m128i l, _mm_loadu_si128, src_left + 1, 8)
    LOAD2(__m128i a, _mm_loadu_si128, src_above + 1, 8)
    l1 = _mm_add_epi16(l1, a1);
    l2 = _mm_add_epi16(l2, a2);
    l1 = _mm_hadd_epi16(l1, l2);
    l1 = _mm_hadd_epi16(l1, l1);
    l1 = _mm_hadd_epi16(l1, l1);
    l1 = _mm_hadd_epi16(l1, l1);
    dc_val = _mm_extract_epi16(l1, 0);

    dc_val = (dc_val + offset) >> shift;

    VECTOR_DEFINITION_SCALE1();
    for (int y = 0; y < 7; ++y){
        LOOP16_PDPC();
    }
    for (int y = 7; y < (1 << log2_pb_h); ++y){
        LOOP16_NOPDPC();
    }
}

static void
vvc_intra_dc_pdpc_w16_h32_sse(const uint16_t *const src_above,
                      const uint16_t *const src_left,
                      uint16_t *const dst, ptrdiff_t dst_stride,
                      int log2_pb_w, int log2_pb_h)
{
    int idx;
    uint16_t *_dst = dst;
    uint32_t dc_val = 0;
    const int shift  = 5;
    const int offset = 0x10;

    LOAD4(__m128i l, _mm_loadu_si128, src_left + 1, 8)
    LOAD2(__m128i a, _mm_loadu_si128, src_above + 1, 8)
    l1 = _mm_add_epi16(l1, l2);
    l3 = _mm_add_epi16(l3, l4);
    l1 = _mm_add_epi16(l1, l3);
    l1 = _mm_hadd_epi16(l1, l1);
    l1 = _mm_hadd_epi16(l1, l1);
    l1 = _mm_hadd_epi16(l1, l1);
    dc_val = _mm_extract_epi16(l1, 0);

    dc_val = (dc_val + offset) >> shift;

    VECTOR_DEFINITION_SCALE1();
    for (int y = 0; y < 7; ++y){
        LOOP16_PDPC();
    }
    for (int y = 7; y < (1 << log2_pb_h); ++y){
        LOOP16_NOPDPC();
    }
}

static void
vvc_intra_dc_pdpc_w16_h64_sse(const uint16_t *const src_above,
                      const uint16_t *const src_left,
                      uint16_t *const dst, ptrdiff_t dst_stride,
                      int log2_pb_w, int log2_pb_h)
{
    int idx;
    uint16_t *_dst = dst;
    uint32_t dc_val = 0;
    const int shift  = 6;
    const int offset = 0x20;

    LOAD8(__m128i l, _mm_loadu_si128, src_left + 1, 8)
    LOAD2(__m128i a, _mm_loadu_si128, src_above + 1, 8)
    l1 = _mm_add_epi16(l1, l5);
    l2 = _mm_add_epi16(l2, l6);
    l3 = _mm_add_epi16(l3, l7);
    l4 = _mm_add_epi16(l4, l8);
    l1 = _mm_add_epi16(l1, l2);
    l3 = _mm_add_epi16(l3, l4);
    l1 = _mm_add_epi16(l1, l3);
    l1 = _mm_hadd_epi16(l1, l1);
    l1 = _mm_hadd_epi16(l1, l1);
    l1 = _mm_hadd_epi16(l1, l1);
    dc_val = _mm_extract_epi16(l1, 0);

    dc_val = (dc_val + offset) >> shift;

    VECTOR_DEFINITION_SCALE2();
    for (int y = 0; y < 13; ++y){
        LOOP16_PDPC2();
    }
    for (int y = 13; y < (1 << log2_pb_h); ++y){
        LOOP16_NOPDPC2();
    }
}

static void
vvc_intra_dc_pdpc_w32_h4_sse(const uint16_t *const src_above,
                      const uint16_t *const src_left,
                      uint16_t *const dst, ptrdiff_t dst_stride,
                      int log2_pb_w, int log2_pb_h)
{
    int idx;
    uint16_t *_dst = dst;
    uint32_t dc_val = 0;
    const int shift  = 5;
    const int offset = 0x10;

    LOAD4(__m128i a, _mm_loadu_si128, src_above + 1, 8)
    __m128i l1 = _mm_add_epi16(a1, a2);
    __m128i l2 = _mm_add_epi16(a3, a4);
    l1 = _mm_add_epi16(l1, l2);
    l1 = _mm_hadd_epi16(l1, l1);
    l1 = _mm_hadd_epi16(l1, l1);
    l1 = _mm_hadd_epi16(l1, l1);
    dc_val = _mm_extract_epi16(l1, 0);

    dc_val = (dc_val + offset) >> shift;

    VECTOR_DEFINITION_SCALE1();
    for (int y = 0; y < (1 << log2_pb_h); ++y){
        LOOP32_PDPC();
    }
}

static void
vvc_intra_dc_pdpc_w32_h8_sse(const uint16_t *const src_above,
                      const uint16_t *const src_left,
                      uint16_t *const dst, ptrdiff_t dst_stride,
                      int log2_pb_w, int log2_pb_h)
{
    int idx;
    uint16_t *_dst = dst;
    uint32_t dc_val = 0;
    const int shift  = 5;
    const int offset = 0x10;

    LOAD4(__m128i a, _mm_loadu_si128, src_above + 1, 8)
    __m128i l1 = _mm_add_epi16(a1, a2);
    __m128i l2 = _mm_add_epi16(a3, a4);
    l1 = _mm_add_epi16(l1, l2);
    l1 = _mm_hadd_epi16(l1, l1);
    l1 = _mm_hadd_epi16(l1, l1);
    l1 = _mm_hadd_epi16(l1, l1);
    dc_val = _mm_extract_epi16(l1, 0);

    dc_val = (dc_val + offset) >> shift;

    VECTOR_DEFINITION_SCALE1();
    for (int y = 0; y < 7; ++y){
        LOOP32_PDPC();
    }
    for (int y = 7; y < (1 << log2_pb_h); ++y){
        LOOP32_NOPDPC();
    }
}

static void
vvc_intra_dc_pdpc_w32_h32_sse(const uint16_t *const src_above,
                      const uint16_t *const src_left,
                      uint16_t *const dst, ptrdiff_t dst_stride,
                      int log2_pb_w, int log2_pb_h)
{
    int idx;
    uint16_t *_dst = dst;
    uint32_t dc_val = 0;
    const int shift  = 6;
    const int offset = 0x20;

    LOAD4(__m128i a, _mm_loadu_si128, src_above + 1, 8)
    LOAD4(__m128i l, _mm_loadu_si128, src_left + 1, 8)
    l1 = _mm_add_epi16(l1, a1);
    l2 = _mm_add_epi16(l2, a2);
    l3 = _mm_add_epi16(l3, a3);
    l4 = _mm_add_epi16(l4, a4);
    l1 = _mm_add_epi16(l1, l2);
    l3 = _mm_add_epi16(l3, l4);
    l1 = _mm_hadd_epi16(l1, l3);
    l1 = _mm_hadd_epi16(l1, l1);
    l1 = _mm_hadd_epi16(l1, l1);
    l1 = _mm_hadd_epi16(l1, l1);
    dc_val = _mm_extract_epi16(l1, 0);

    dc_val = (dc_val + offset) >> shift;

    VECTOR_DEFINITION_SCALE2();

    for (int y = 0; y < 13; ++y){
        LOOP32_PDPC2();
    }
    for (int y = 13; y < (1 << log2_pb_h); ++y){
        LOOP32_NOPDPC2();
    }
}

static void
vvc_intra_dc_pdpc_w32_h64_sse(const uint16_t *const src_above,
                      const uint16_t *const src_left,
                      uint16_t *const dst, ptrdiff_t dst_stride,
                      int log2_pb_w, int log2_pb_h)
{
    int idx;
    uint16_t *_dst = dst;
    uint32_t dc_val = 0;
    const int shift  = 6;
    const int offset = 0x20;

    LOAD4(__m128i a, _mm_loadu_si128, src_above + 1, 8)
    LOAD8(__m128i l, _mm_loadu_si128, src_left + 1, 8)
    l1 = _mm_add_epi16(l1, l5);
    l2 = _mm_add_epi16(l2, l6);
    l3 = _mm_add_epi16(l3, l7);
    l4 = _mm_add_epi16(l4, l8);
    l1 = _mm_add_epi16(l1, l2);
    l3 = _mm_add_epi16(l3, l4);
    l1 = _mm_hadd_epi16(l1, l3);
    l1 = _mm_hadd_epi16(l1, l1);
    l1 = _mm_hadd_epi16(l1, l1);
    l1 = _mm_hadd_epi16(l1, l1);
    dc_val = _mm_extract_epi16(l1, 0);

    dc_val = (dc_val + offset) >> shift;

    VECTOR_DEFINITION_SCALE2();

    for (int y = 0; y < 13; ++y){
        LOOP32_PDPC2();
    }
    for (int y = 13; y < (1 << log2_pb_h); ++y){
        LOOP32_NOPDPC2();
    }
}

static void
vvc_intra_dc_pdpc_w64_h4_sse(const uint16_t *const src_above,
                      const uint16_t *const src_left,
                      uint16_t *const dst, ptrdiff_t dst_stride,
                      int log2_pb_w, int log2_pb_h)
{
    int idx;
    uint16_t *_dst = dst;
    uint32_t dc_val = 0;
    const int shift  = 6;
    const int offset = 0x20;

    LOAD8(__m128i a, _mm_loadu_si128, src_above + 1, 8)
    __m128i l1 = _mm_add_epi16(a1, a5);
    __m128i l2 = _mm_add_epi16(a2, a6);
    __m128i l3 = _mm_add_epi16(a3, a7);
    __m128i l4 = _mm_add_epi16(a4, a8);
    l1 = _mm_add_epi16(l1, l2);
    l3 = _mm_add_epi16(l3, l4);
    l1 = _mm_add_epi16(l1, l3);
    l1 = _mm_hadd_epi16(l1, l1);
    l1 = _mm_hadd_epi16(l1, l1);
    l1 = _mm_hadd_epi16(l1, l1);
    dc_val = _mm_extract_epi16(l1, 0);

    dc_val = (dc_val + offset) >> shift;

    VECTOR_DEFINITION_SCALE1();

    for (int y = 0; y < (1 << log2_pb_h); ++y){
        LOOP64_PDPC();
    }
}

static void
vvc_intra_dc_pdpc_w64_h8_sse(const uint16_t *const src_above,
                      const uint16_t *const src_left,
                      uint16_t *const dst, ptrdiff_t dst_stride,
                      int log2_pb_w, int log2_pb_h)
{
    int idx;
    uint16_t *_dst = dst;
    uint32_t dc_val = 0;
    const int shift  = 6;
    const int offset = 0x20;

    LOAD8(__m128i a, _mm_loadu_si128, src_above + 1, 8)
    __m128i l1 = _mm_add_epi16(a1, a5);
    __m128i l2 = _mm_add_epi16(a2, a6);
    __m128i l3 = _mm_add_epi16(a3, a7);
    __m128i l4 = _mm_add_epi16(a4, a8);
    l1 = _mm_add_epi16(l1, l2);
    l3 = _mm_add_epi16(l3, l4);
    l1 = _mm_add_epi16(l1, l3);
    l1 = _mm_hadd_epi16(l1, l1);
    l1 = _mm_hadd_epi16(l1, l1);
    l1 = _mm_hadd_epi16(l1, l1);
    dc_val = _mm_extract_epi16(l1, 0);

    dc_val = (dc_val + offset) >> shift;

    VECTOR_DEFINITION_SCALE1();

    for (int y = 0; y < 7; ++y){
        LOOP64_PDPC();
    }
    for (int y = 7; y < (1 << log2_pb_h); ++y){
        LOOP64_NOPDPC();
    }
}

static void
vvc_intra_dc_pdpc_w64_h16_sse(const uint16_t *const src_above,
                      const uint16_t *const src_left,
                      uint16_t *const dst, ptrdiff_t dst_stride,
                      int log2_pb_w, int log2_pb_h)
{
    int idx;
    uint16_t *_dst = dst;
    uint32_t dc_val = 0;
    const int shift  = 6;
    const int offset = 0x20;

    LOAD8(__m128i a, _mm_loadu_si128, src_above + 1, 8)
    __m128i l1 = _mm_add_epi16(a1, a5);
    __m128i l2 = _mm_add_epi16(a2, a6);
    __m128i l3 = _mm_add_epi16(a3, a7);
    __m128i l4 = _mm_add_epi16(a4, a8);
    l1 = _mm_add_epi16(l1, l2);
    l3 = _mm_add_epi16(l3, l4);
    l1 = _mm_add_epi16(l1, l3);
    l1 = _mm_hadd_epi16(l1, l1);
    l1 = _mm_hadd_epi16(l1, l1);
    l1 = _mm_hadd_epi16(l1, l1);
    dc_val = _mm_extract_epi16(l1, 0);

    dc_val = (dc_val + offset) >> shift;

    VECTOR_DEFINITION_SCALE2();

    for (int y = 0; y < 13; ++y){
        LOOP64_PDPC2();
    }
    for (int y = 13; y < (1 << log2_pb_h); ++y){
        LOOP64_NOPDPC2();
    }
}

static void
vvc_intra_dc_pdpc_w64_h64_sse(const uint16_t *const src_above,
                      const uint16_t *const src_left,
                      uint16_t *const dst, ptrdiff_t dst_stride,
                      int log2_pb_w, int log2_pb_h)
{
    int idx;
    uint16_t *_dst = dst;
    uint32_t dc_val = 0;
    const int shift  = 7;
    const int offset = 0x40;

    LOAD8(__m128i a, _mm_loadu_si128, src_above + 1, 8)
    LOAD8(__m128i l, _mm_loadu_si128, src_left + 1, 8)
    l1 = _mm_add_epi16(l1, a1);
    l2 = _mm_add_epi16(l2, a2);
    l3 = _mm_add_epi16(l3, a3);
    l4 = _mm_add_epi16(l4, a4);
    l5 = _mm_add_epi16(l5, a5);
    l6 = _mm_add_epi16(l6, a6);
    l7 = _mm_add_epi16(l7, a7);
    l8 = _mm_add_epi16(l8, a8);
    l1 = _mm_add_epi16(l1, l2);
    l3 = _mm_add_epi16(l3, l4);
    l5 = _mm_add_epi16(l5, l6);
    l7 = _mm_add_epi16(l7, l8);
    l1 = _mm_add_epi16(l1, l3);
    l5 = _mm_add_epi16(l5, l7);
    l2 = _mm_unpackhi_epi16(l1, _mm_setzero_si128());
    l1 = _mm_unpacklo_epi16(l1, _mm_setzero_si128());
    l6 = _mm_unpackhi_epi16(l5, _mm_setzero_si128());
    l5 = _mm_unpacklo_epi16(l5, _mm_setzero_si128());
    l1 = _mm_add_epi32(l1, l5);
    l5 = _mm_add_epi32(l2, l6);
    l1 = _mm_hadd_epi32(l1, l5);
    l1 = _mm_hadd_epi32(l1, l1);
    l1 = _mm_hadd_epi32(l1, l1);
    dc_val = _mm_extract_epi32(l1, 0);

    dc_val = (dc_val + offset) >> shift;

    VECTOR_DEFINITION_SCALE2();

    for (int y = 0; y < 13; ++y){
        LOOP64_PDPC2();
    }
    for (int y = 13; y < (1 << log2_pb_h); ++y){
        LOOP64_NOPDPC2();
    }
}

#define CUT_PDPC 0

#define PLANAR_VECTOR_DEF_SCALE01(scale, size)  const uint8_t *pdpc_w = vvc_pdpc_w[scale];                    \
                                                __m128i tr_val_v = _mm_set1_epi16(src_above[width + 1]);      \
                                                __m128i bl_val_v = _mm_set1_epi16(src_left[height + 1]);      \
                                                __m128i rnd_v = _mm_set1_epi ## size(1 << max_scale);         \
                                                __m128i mul_col = _mm_set_epi16(8, 7, 6, 5, 4, 3, 2, 1);      \
                                                __m128i add_v = _mm_set1_epi16(32);                           \
                                                __m128i x_v = _mm_loadu_si128((__m128i *) (pdpc_w));          \
                                                __m128i x_v1 = _mm_unpacklo_epi8(x_v, _mm_setzero_si128());

#define PLANAR_VECTOR_DEF_SCALE2(size)      const uint8_t *pdpc_w = vvc_pdpc_w[2];                        \
                                            __m128i tr_val_v = _mm_set1_epi16(src_above[width + 1]);      \
                                            __m128i bl_val_v = _mm_set1_epi16(src_left[height + 1]);      \
                                            __m128i rnd_v = _mm_set1_epi ## size(1 << max_scale);         \
                                            __m128i mul_col = _mm_set_epi16(8, 7, 6, 5, 4, 3, 2, 1);      \
                                            __m128i add_v = _mm_set1_epi16(32);                           \
                                            __m128i x_v = _mm_loadu_si128((__m128i *) (pdpc_w));          \
                                            __m128i x_v1 = _mm_unpacklo_epi8(x_v, _mm_setzero_si128());   \
                                            __m128i x_v2 = _mm_unpackhi_epi8(x_v, _mm_setzero_si128());

#define INIT_LOOP_PDPC()    __m128i rcol_v = _mm_set1_epi16(r_col[y]);                            \
                            __m128i lcol_v = _mm_set1_epi16(l_col[y]);                            \
                            int y_wgh = pdpc_w[y];                                                \
                            int32_t l_val = (int32_t)src_left[y + 1];                             \
                            __m128i l_v = _mm_set1_epi16(l_val);                                  \
                            __m128i y_v = _mm_set1_epi16(y_wgh);

#define INIT_LOOP_NOPDPC()  __m128i rcol_v = _mm_set1_epi16(r_col[y]);                            \
                            __m128i lcol_v = _mm_set1_epi16(l_col[y]);                            \
                            int32_t l_val = (int32_t)src_left[y + 1];                             \
                            __m128i l_v = _mm_set1_epi16(l_val);


#define COMPUTE_RCOL_X_1()    __m128i rcol_mul_v = _mm_mullo_epi16(rcol_v, mul_col);                \
                            __m128i rcol_x1 = _mm_add_epi16(lcol_v, rcol_mul_v);

#define COMPUTE_RCOL_X_2()    __m128i rcol_mul_v = _mm_mullo_epi16(rcol_v, mul_col);                \
                            __m128i rcol_x1 = _mm_add_epi16(lcol_v, rcol_mul_v);                    \
                            rcol_mul_v = _mm_add_epi16(rcol_mul_v, _mm_slli_epi16(rcol_v, 3));      \
                            __m128i rcol_x2 = _mm_add_epi16(lcol_v, rcol_mul_v);

#define COMPUTE_RCOL_X_4()    __m128i rcol_mul_v = _mm_mullo_epi16(rcol_v, mul_col);                \
                            __m128i rcol_x1 = _mm_add_epi16(lcol_v, rcol_mul_v);                    \
                            rcol_mul_v = _mm_add_epi16(rcol_mul_v, _mm_slli_epi16(rcol_v, 3));      \
                            __m128i rcol_x2 = _mm_add_epi16(lcol_v, rcol_mul_v);                    \
                            rcol_mul_v = _mm_add_epi16(rcol_mul_v, _mm_slli_epi16(rcol_v, 3));      \
                            __m128i rcol_x3 = _mm_add_epi16(lcol_v, rcol_mul_v);                    \
                            rcol_mul_v = _mm_add_epi16(rcol_mul_v, _mm_slli_epi16(rcol_v, 3));      \
                            __m128i rcol_x4 = _mm_add_epi16(lcol_v, rcol_mul_v);

#define COMPUTE_RCOL_X_8()  __m128i rcol_mul_v = _mm_mullo_epi16(rcol_v, mul_col);                  \
                            __m128i rcol_x1 = _mm_add_epi16(lcol_v, rcol_mul_v);                    \
                            rcol_mul_v = _mm_add_epi16(rcol_mul_v, _mm_slli_epi16(rcol_v, 3));      \
                            __m128i rcol_x2 = _mm_add_epi16(lcol_v, rcol_mul_v);                    \
                            rcol_mul_v = _mm_add_epi16(rcol_mul_v, _mm_slli_epi16(rcol_v, 3));      \
                            __m128i rcol_x3 = _mm_add_epi16(lcol_v, rcol_mul_v);                    \
                            rcol_mul_v = _mm_add_epi16(rcol_mul_v, _mm_slli_epi16(rcol_v, 3));      \
                            __m128i rcol_x4 = _mm_add_epi16(lcol_v, rcol_mul_v);                    \
                            rcol_mul_v = _mm_add_epi16(rcol_mul_v, _mm_slli_epi16(rcol_v, 3));      \
                            __m128i rcol_x5 = _mm_add_epi16(lcol_v, rcol_mul_v);                    \
                            rcol_mul_v = _mm_add_epi16(rcol_mul_v, _mm_slli_epi16(rcol_v, 3));      \
                            __m128i rcol_x6 = _mm_add_epi16(lcol_v, rcol_mul_v);                    \
                            rcol_mul_v = _mm_add_epi16(rcol_mul_v, _mm_slli_epi16(rcol_v, 3));      \
                            __m128i rcol_x7 = _mm_add_epi16(lcol_v, rcol_mul_v);                    \
                            rcol_mul_v = _mm_add_epi16(rcol_mul_v, _mm_slli_epi16(rcol_v, 3));      \
                            __m128i rcol_x8 = _mm_add_epi16(lcol_v, rcol_mul_v);

#define COMPUTE_PLANAR_VECTOR(unroll_val)   UNROLL ## unroll_val(U111, tr_v, _mm_sub_epi16, tr_v, a)                        \
                                            UNROLL ## unroll_val(U110, tr_v, _mm_add_epi16, tr_v, bl_val_v)                 \
                                                                                                                            \
                                            COMPUTE_RCOL_X_## unroll_val()                                                  \
                                                                                                                            \
                                            UNROLL ## unroll_val(U110, __m128i src_v_lo, _mm_slli_epi16, rcol_x, rc_scale)  \
                                            UNROLL ## unroll_val(U110, __m128i str_v_lo, _mm_slli_epi16, tr_v, tr_scale)    \
                                                                                                                            \
                                            UNROLL ## unroll_val(U111, __m128i out_lo, _mm_add_epi16, src_v_lo, str_v_lo)   \
                                            UNROLL ## unroll_val(U110, out_lo, _mm_add_epi16, out_lo, rnd_v)                \
                                            UNROLL ## unroll_val(U110, out_lo, _mm_srli_epi16, out_lo, max_scale + 1)

#define COMPUTE_PLANAR_VECTOR2(unroll_val)  UNROLL ## unroll_val(U111, tr_v, _mm_sub_epi16, tr_v, a)                                        \
                                            UNROLL ## unroll_val(U110, tr_v, _mm_add_epi16, tr_v, bl_val_v)                                 \
                                                                                                                                            \
                                            COMPUTE_RCOL_X_## unroll_val()                                                                  \
                                            UNROLL ## unroll_val(U110, __m128i tr_v_lo, _mm_unpacklo_epi16, tr_v, _mm_setzero_si128())      \
                                            UNROLL ## unroll_val(U110, __m128i tr_v_hi, _mm_unpackhi_epi16, tr_v, _mm_setzero_si128())      \
                                                                                                                                            \
                                            UNROLL ## unroll_val(U110, __m128i rc_v_lo, _mm_unpacklo_epi16, rcol_x, _mm_setzero_si128())    \
                                            UNROLL ## unroll_val(U110, __m128i rc_v_hi, _mm_unpackhi_epi16, rcol_x, _mm_setzero_si128())    \
                                                                                                                                            \
                                            UNROLL ## unroll_val(U110, __m128i src_v_lo, _mm_slli_epi32, rc_v_lo, rc_scale)                 \
                                            UNROLL ## unroll_val(U110, __m128i src_v_hi, _mm_slli_epi32, rc_v_hi, rc_scale)                 \
                                            UNROLL ## unroll_val(U110, __m128i str_v_lo, _mm_slli_epi32, tr_v_lo, tr_scale)                 \
                                            UNROLL ## unroll_val(U110, __m128i str_v_hi, _mm_slli_epi32, tr_v_hi, tr_scale)                 \
                                                                                                                                            \
                                            UNROLL ## unroll_val(U111, __m128i out_lo, _mm_add_epi32, src_v_lo, str_v_lo)                   \
                                            UNROLL ## unroll_val(U111, __m128i out_hi, _mm_add_epi32, src_v_hi, str_v_hi)                   \
                                            UNROLL ## unroll_val(U110, out_lo, _mm_add_epi32, out_lo, rnd_v)                                \
                                            UNROLL ## unroll_val(U110, out_hi, _mm_add_epi32, out_hi, rnd_v)                                \
                                            UNROLL ## unroll_val(U110, out_lo, _mm_srli_epi32, out_lo, max_scale + 1)                       \
                                            UNROLL ## unroll_val(U110, out_hi, _mm_srli_epi32, out_hi, max_scale + 1)                       \
                                            UNROLL ## unroll_val(U111, out_lo, _mm_packs_epi32, out_lo, out_hi)

#define APPLY_PDPC(unroll_val, pdpc)    UNROLL ## unroll_val(U110, __m128i tst, _mm_slli_epi16, out_lo, 6)                  \
                                                                                                                            \
                                        UNROLL ## unroll_val(U110, __m128i w_x, _mm_mullo_epi16, out_lo, y_v)               \
                                        UNROLL ## pdpc(U111, __m128i w_y, _mm_mullo_epi16, out_lo, x_v)                     \
                                                                                                                            \
                                        UNROLL ## unroll_val(U111, tst, _mm_subs_epu16, tst, w_x)                           \
                                        UNROLL ## pdpc(U111, tst, _mm_subs_epu16, tst, w_y)                                 \
                                        UNROLL ## unroll_val(U110, tst, _mm_adds_epu16, tst, add_v)                         \
                                                                                                                            \
                                        UNROLL ## pdpc(U110, __m128i xl, _mm_mullo_epi16, x_v, l_v)                         \
                                        UNROLL ## unroll_val(U101, __m128i yt, _mm_mullo_epi16, y_v, a)                     \
                                                                                                                            \
                                        UNROLL ## pdpc(U111, tst, _mm_adds_epu16, tst, xl)                                  \
                                                                                                                            \
                                        UNROLL ## unroll_val(U111, tst, _mm_adds_epu16, tst, yt)                            \
                                                                                                                            \
                                        UNROLL ## unroll_val(U110, __m128i out_v, _mm_srli_epi16, tst, 6)                   \
                                                                                                                            \
                                        UNROLL ## unroll_val(U110, out_v, _mm_min_epi16, out_v, _mm_set1_epi16(1023))       \
                                        UNROLL ## unroll_val(U110, out_v, _mm_max_epi16, out_v, _mm_setzero_si128())

#define APPLY_NOPDPC(unroll_val, pdpc)  UNROLL ## unroll_val(U110, __m128i tst, _mm_slli_epi16, out_lo, 6)                  \
                                                                                                                            \
                                        UNROLL ## pdpc(U111, __m128i w_y, _mm_mullo_epi16, out_lo, x_v)                     \
                                                                                                                            \
                                        UNROLL ## pdpc(U111, tst, _mm_subs_epu16, tst, w_y)                                 \
                                        UNROLL ## unroll_val(U110, tst, _mm_adds_epu16, tst, add_v)                         \
                                                                                                                            \
                                        UNROLL ## pdpc(U110, __m128i xl, _mm_mullo_epi16, x_v, l_v)                         \
                                                                                                                            \
                                        UNROLL ## pdpc(U111, tst, _mm_adds_epu16, tst, xl)                                  \
                                                                                                                            \
                                        UNROLL ## unroll_val(U110, __m128i out_v, _mm_srli_epi16, tst, 6)                   \
                                                                                                                            \
                                        UNROLL ## unroll_val(U110, out_v, _mm_min_epi16, out_v, _mm_set1_epi16(1023))       \
                                        UNROLL ## unroll_val(U110, out_v, _mm_max_epi16, out_v, _mm_setzero_si128())

void
vvc_intra_planar_pdpc_w4_h4_sse(const uint16_t *const src_above,
                          const uint16_t *const src_left,
                          uint16_t *const dst, ptrdiff_t dst_stride,
                          int log2_pb_w, int log2_pb_h)
{
    uint16_t *_dst = dst;
    const uint32_t width  = 1 << log2_pb_w;
    const uint32_t height = 1 << log2_pb_h;
    int max_scale = 2;
    int rc_scale = 0;
    int tr_scale = 0;

    const uint8_t pdpc_scale = (log2_pb_w + log2_pb_h - 2) >> 2;

    int16_t r_col[4], l_col[4];

    PLANAR_VECTOR_DEF_SCALE01(0, 16)

    LOAD1(__m128i a, _mm_loadl_epi64, src_above + 1, 8)
    UNROLL1(U110, __m128i tr_v, _mm_slli_epi16, a, log2_pb_h)

    LOAD1(__m128i l, _mm_loadl_epi64, src_left + 1, 8)
    UNROLL1(U101, __m128i rc_v, _mm_sub_epi16, tr_val_v, l)
    UNROLL1(U110, __m128i lc_v, _mm_slli_epi16, l, log2_pb_w)

    STORE1(_mm_storel_epi64, (r_col),rc_v, 8)
    STORE1(_mm_storel_epi64, (l_col),lc_v, 8)

    for (int y = 0; y < 3; ++y){
        INIT_LOOP_PDPC()

        COMPUTE_PLANAR_VECTOR(1)

        APPLY_PDPC(1, 1)

        STORE1(_mm_storel_epi64, _dst, out_v, 8)
        _dst += dst_stride;
    }
    for (int y = 3; y < height; ++y){
        INIT_LOOP_NOPDPC()

        COMPUTE_PLANAR_VECTOR(1)

        APPLY_NOPDPC(1, 1)

        STORE1(_mm_storel_epi64, _dst, out_v, 8)
        _dst += dst_stride;
    }
}

void
vvc_intra_planar_pdpc_w4_h8_sse(const uint16_t *const src_above,
                          const uint16_t *const src_left,
                          uint16_t *const dst, ptrdiff_t dst_stride,
                          int log2_pb_w, int log2_pb_h)
{
    uint16_t *_dst = dst;
    const uint32_t width  = 1 << log2_pb_w;
    const uint32_t height = 1 << log2_pb_h;
    int max_scale = 3;
    int rc_scale = 1;
    int tr_scale = 0;

    int16_t r_col[8], l_col[8];

    PLANAR_VECTOR_DEF_SCALE01(0, 16)

    LOAD1(__m128i a, _mm_loadl_epi64, src_above + 1, 8)
    UNROLL1(U110, __m128i tr_v, _mm_slli_epi16, a, log2_pb_h)

    LOAD1(__m128i l, _mm_loadu_si128, src_left + 1, 8)
    UNROLL1(U101, __m128i rc_v, _mm_sub_epi16, tr_val_v, l)
    UNROLL1(U110, __m128i lc_v, _mm_slli_epi16, l, log2_pb_w)

    STORE1(_mm_storeu_si128, (r_col),rc_v, 8)
    STORE1(_mm_storeu_si128, (l_col),lc_v, 8)

    for (int y = 0; y < 3; ++y){
        INIT_LOOP_PDPC()

        COMPUTE_PLANAR_VECTOR(1)

        APPLY_PDPC(1, 1)

        STORE1(_mm_storel_epi64, _dst, out_v, 8)
        _dst += dst_stride;
    }
    for (int y = 3; y < height; ++y){
        INIT_LOOP_NOPDPC()

        COMPUTE_PLANAR_VECTOR(1)

        APPLY_NOPDPC(1, 1)

        STORE1(_mm_storel_epi64, _dst, out_v, 8)
        _dst += dst_stride;
    }
}

void
vvc_intra_planar_pdpc_w4_h16_sse(const uint16_t *const src_above,
                          const uint16_t *const src_left,
                          uint16_t *const dst, ptrdiff_t dst_stride,
                          int log2_pb_w, int log2_pb_h)
{
    uint16_t *_dst = dst;
    const uint32_t width  = 1 << log2_pb_w;
    const uint32_t height = 1 << log2_pb_h;
    int max_scale = 4;
    int rc_scale = 2;
    int tr_scale = 0;

    int16_t r_col[16], l_col[16];

    PLANAR_VECTOR_DEF_SCALE01(1, 16)

    LOAD1(__m128i a, _mm_loadl_epi64, src_above + 1, 8)
    UNROLL1(U110, __m128i tr_v, _mm_slli_epi16, a, log2_pb_h)

    LOAD2(__m128i l, _mm_loadu_si128, src_left + 1, 8)
    UNROLL2(U101, __m128i rc_v, _mm_sub_epi16, tr_val_v, l)
    UNROLL2(U110, __m128i lc_v, _mm_slli_epi16, l, log2_pb_w)

    STORE2(_mm_storeu_si128, r_col, rc_v, 8)
    STORE2(_mm_storeu_si128, l_col, lc_v, 8)

    for (int y = 0; y < 7; ++y){
        INIT_LOOP_PDPC()

        COMPUTE_PLANAR_VECTOR(1)

        APPLY_PDPC(1, 1)

        STORE1(_mm_storel_epi64, _dst, out_v, 8)
        _dst += dst_stride;
    }
    for (int y = 7; y < height; ++y){
        INIT_LOOP_NOPDPC()

        COMPUTE_PLANAR_VECTOR(1)

        APPLY_NOPDPC(1, 1)

        STORE1(_mm_storel_epi64, _dst, out_v, 8)
        _dst += dst_stride;
    }
}

void
vvc_intra_planar_pdpc_w4_h32_sse(const uint16_t *const src_above,
                          const uint16_t *const src_left,
                          uint16_t *const dst, ptrdiff_t dst_stride,
                          int log2_pb_w, int log2_pb_h)
{
    uint16_t *_dst = dst;
    const uint32_t width  = 1 << log2_pb_w;
    const uint32_t height = 1 << log2_pb_h;
    int max_scale = 5;
    int rc_scale = 3;
    int tr_scale = 0;

    int16_t r_col[32], l_col[32];

    PLANAR_VECTOR_DEF_SCALE01(1, 16)

    LOAD1(__m128i a, _mm_loadl_epi64, src_above + 1, 8)
    UNROLL1(U110, __m128i tr_v, _mm_slli_epi16, a, log2_pb_h)

    LOAD4(__m128i l, _mm_loadu_si128, src_left + 1, 8)
    UNROLL4(U101, __m128i rc_v, _mm_sub_epi16, tr_val_v, l)
    UNROLL4(U110, __m128i lc_v, _mm_slli_epi16, l, log2_pb_w)

    STORE4(_mm_storeu_si128, r_col, rc_v, 8)
    STORE4(_mm_storeu_si128, l_col, lc_v, 8)

    for (int y = 0; y < 7; ++y){
        INIT_LOOP_PDPC()

        COMPUTE_PLANAR_VECTOR(1)

        APPLY_PDPC(1, 1)

        STORE1(_mm_storel_epi64, _dst, out_v, 8)
        _dst += dst_stride;
    }
    for (int y = 7; y < height; ++y){
        INIT_LOOP_NOPDPC()

        COMPUTE_PLANAR_VECTOR(1)

        APPLY_NOPDPC(1, 1)

        STORE1(_mm_storel_epi64, _dst, out_v, 8)
        _dst += dst_stride;
    }
}

void
vvc_intra_planar_pdpc_w4_h64_sse(const uint16_t *const src_above,
                          const uint16_t *const src_left,
                          uint16_t *const dst, ptrdiff_t dst_stride,
                          int log2_pb_w, int log2_pb_h)
{   
    uint16_t *_dst = dst;
    const uint32_t width  = 1 << log2_pb_w;
    const uint32_t height = 1 << log2_pb_h;
    int max_scale = 6;
    int rc_scale = 4;
    int tr_scale = 0;

    int16_t r_col[64], l_col[64];

    PLANAR_VECTOR_DEF_SCALE01(1, 32)

    LOAD1(__m128i a, _mm_loadl_epi64, src_above + 1, 8)
    UNROLL1(U110, __m128i tr_v, _mm_slli_epi16, a, log2_pb_h)

    LOAD8(__m128i l, _mm_loadu_si128, src_left + 1, 8)
    UNROLL8(U101, __m128i rc_v, _mm_sub_epi16, tr_val_v, l)
    UNROLL8(U110, __m128i lc_v, _mm_slli_epi16, l, log2_pb_w)

    STORE8(_mm_storeu_si128, r_col, rc_v, 8)
    STORE8(_mm_storeu_si128, l_col, lc_v, 8)

    for (int y = 0; y < 7; ++y){
        INIT_LOOP_PDPC()

        COMPUTE_PLANAR_VECTOR2(1)

        APPLY_PDPC(1, 1)

        STORE1(_mm_storel_epi64, _dst, out_v, 8)
        _dst += dst_stride;
    }
    for (int y = 7; y < height; ++y){
        INIT_LOOP_NOPDPC()

        COMPUTE_PLANAR_VECTOR2(1)

        APPLY_NOPDPC(1, 1)

        STORE1(_mm_storel_epi64, _dst, out_v, 8)
        _dst += dst_stride;
    }
}

void
vvc_intra_planar_pdpc_w8_h4_sse(const uint16_t *const src_above,
                          const uint16_t *const src_left,
                          uint16_t *const dst, ptrdiff_t dst_stride,
                          int log2_pb_w, int log2_pb_h)
{
    uint16_t *_dst = dst;
    const uint32_t width  = 1 << log2_pb_w;
    const uint32_t height = 1 << log2_pb_h;
    int max_scale = 3;
    int rc_scale = 0;
    int tr_scale = 1;

    int16_t r_col[4], l_col[4];

    PLANAR_VECTOR_DEF_SCALE01(0, 16)

    LOAD1(__m128i a, _mm_loadu_si128, src_above + 1, 8)
    UNROLL1(U110, __m128i tr_v, _mm_slli_epi16, a, log2_pb_h)

    LOAD1(__m128i l, _mm_loadl_epi64, src_left + 1, 8)
    UNROLL1(U101, __m128i rc_v, _mm_sub_epi16, tr_val_v, l)
    UNROLL1(U110, __m128i lc_v, _mm_slli_epi16, l, log2_pb_w)

    STORE1(_mm_storel_epi64, (r_col),rc_v, 8)
    STORE1(_mm_storel_epi64, (l_col),lc_v, 8)

    for (int y = 0; y < 3; ++y){
        INIT_LOOP_PDPC()

        COMPUTE_PLANAR_VECTOR(1)

        APPLY_PDPC(1, 1)

        STORE1(_mm_storeu_si128, _dst, out_v, 8)
        _dst += dst_stride;
    }
    for (int y = 3; y < height; ++y){
        INIT_LOOP_NOPDPC()

        COMPUTE_PLANAR_VECTOR(1)

        APPLY_NOPDPC(1, 1)

        STORE1(_mm_storeu_si128, _dst, out_v, 8)
        _dst += dst_stride;
    }
}

void
vvc_intra_planar_pdpc_w8_h8_sse(const uint16_t *const src_above,
                          const uint16_t *const src_left,
                          uint16_t *const dst, ptrdiff_t dst_stride,
                          int log2_pb_w, int log2_pb_h)
{
    uint16_t *_dst = dst;
    const uint32_t width  = 1 << log2_pb_w;
    const uint32_t height = 1 << log2_pb_h;
    int max_scale = 3;
    int rc_scale = 0;
    int tr_scale = 0;

    int16_t r_col[8], l_col[8];

    PLANAR_VECTOR_DEF_SCALE01(1, 16)

    LOAD1(__m128i a, _mm_loadu_si128, src_above + 1, 8)
    UNROLL1(U110, __m128i tr_v, _mm_slli_epi16, a, log2_pb_h)

    LOAD1(__m128i l, _mm_loadu_si128, src_left + 1, 8)
    UNROLL1(U101, __m128i rc_v, _mm_sub_epi16, tr_val_v, l)
    UNROLL1(U110, __m128i lc_v, _mm_slli_epi16, l, log2_pb_w)

    STORE1(_mm_storeu_si128, (r_col),rc_v, 8)
    STORE1(_mm_storeu_si128, (l_col),lc_v, 8)

    for (int y = 0; y < 7; ++y){
        INIT_LOOP_PDPC()

        COMPUTE_PLANAR_VECTOR(1)

        APPLY_PDPC(1, 1)

        STORE1(_mm_storeu_si128, _dst, out_v, 8)
        _dst += dst_stride;
    }
    for (int y = 7; y < height; ++y){
        INIT_LOOP_NOPDPC()

        COMPUTE_PLANAR_VECTOR(1)

        APPLY_NOPDPC(1, 1)

        STORE1(_mm_storeu_si128, _dst, out_v, 8)
        _dst += dst_stride;
    }
}

void
vvc_intra_planar_pdpc_w8_h16_sse(const uint16_t *const src_above,
                          const uint16_t *const src_left,
                          uint16_t *const dst, ptrdiff_t dst_stride,
                          int log2_pb_w, int log2_pb_h)
{
    uint16_t *_dst = dst;
    const uint32_t width  = 1 << log2_pb_w;
    const uint32_t height = 1 << log2_pb_h;
    int max_scale = 4;
    int rc_scale = 1;
    int tr_scale = 0;

    int16_t r_col[16], l_col[16];

    PLANAR_VECTOR_DEF_SCALE01(1, 16)

    LOAD1(__m128i a, _mm_loadu_si128, src_above + 1, 8)
    UNROLL1(U110, __m128i tr_v, _mm_slli_epi16, a, log2_pb_h)

    LOAD2(__m128i l, _mm_loadu_si128, src_left + 1, 8)
    UNROLL2(U101, __m128i rc_v, _mm_sub_epi16, tr_val_v, l)
    UNROLL2(U110, __m128i lc_v, _mm_slli_epi16, l, log2_pb_w)

    STORE2(_mm_storeu_si128, (r_col),rc_v, 8)
    STORE2(_mm_storeu_si128, (l_col),lc_v, 8)

    for (int y = 0; y < 7; ++y){
        INIT_LOOP_PDPC()

        COMPUTE_PLANAR_VECTOR(1)

        APPLY_PDPC(1, 1)

        STORE1(_mm_storeu_si128, _dst, out_v, 8)
        _dst += dst_stride;
    }
    for (int y = 7; y < height; ++y){
        INIT_LOOP_NOPDPC()

        COMPUTE_PLANAR_VECTOR(1)

        APPLY_NOPDPC(1, 1)

        STORE1(_mm_storeu_si128, _dst, out_v, 8)
        _dst += dst_stride;
    }
}

void
vvc_intra_planar_pdpc_w8_h32_sse(const uint16_t *const src_above,
                          const uint16_t *const src_left,
                          uint16_t *const dst, ptrdiff_t dst_stride,
                          int log2_pb_w, int log2_pb_h)
{
    uint16_t *_dst = dst;
    const uint32_t width  = 1 << log2_pb_w;
    const uint32_t height = 1 << log2_pb_h;
    int max_scale = 5;
    int rc_scale = 2;
    int tr_scale = 0;

    int16_t r_col[32], l_col[32];

    PLANAR_VECTOR_DEF_SCALE01(1, 16)

    LOAD1(__m128i a, _mm_loadu_si128, src_above + 1, 8)
    UNROLL1(U110, __m128i tr_v, _mm_slli_epi16, a, log2_pb_h)

    LOAD4(__m128i l, _mm_loadu_si128, src_left + 1, 8)
    UNROLL4(U101, __m128i rc_v, _mm_sub_epi16, tr_val_v, l)
    UNROLL4(U110, __m128i lc_v, _mm_slli_epi16, l, log2_pb_w)

    STORE4(_mm_storeu_si128, (r_col),rc_v, 8)
    STORE4(_mm_storeu_si128, (l_col),lc_v, 8)

    for (int y = 0; y < 7; ++y){
        INIT_LOOP_PDPC()

        COMPUTE_PLANAR_VECTOR(1)

        APPLY_PDPC(1, 1)

        STORE1(_mm_storeu_si128, _dst, out_v, 8)
        _dst += dst_stride;
    }
    for (int y = 7; y < height; ++y){
        INIT_LOOP_NOPDPC()

        COMPUTE_PLANAR_VECTOR(1)

        APPLY_NOPDPC(1, 1)

        STORE1(_mm_storeu_si128, _dst, out_v, 8)
        _dst += dst_stride;
    }
}

void
vvc_intra_planar_pdpc_w8_h64_sse(const uint16_t *const src_above,
                          const uint16_t *const src_left,
                          uint16_t *const dst, ptrdiff_t dst_stride,
                          int log2_pb_w, int log2_pb_h)
{   
    uint16_t *_dst = dst;
    const uint32_t width  = 1 << log2_pb_w;
    const uint32_t height = 1 << log2_pb_h;
    int max_scale = 6;
    int rc_scale = 3;
    int tr_scale = 0;

    int16_t r_col[64], l_col[64];

    PLANAR_VECTOR_DEF_SCALE01(1, 32)

    LOAD1(__m128i a, _mm_loadu_si128, src_above + 1, 8)
    UNROLL1(U110, __m128i tr_v, _mm_slli_epi16, a, log2_pb_h)

    LOAD8(__m128i l, _mm_loadu_si128, src_left + 1, 8)
    UNROLL8(U101, __m128i rc_v, _mm_sub_epi16, tr_val_v, l)
    UNROLL8(U110, __m128i lc_v, _mm_slli_epi16, l, log2_pb_w)

    STORE8(_mm_storeu_si128, r_col, rc_v, 8)
    STORE8(_mm_storeu_si128, l_col, lc_v, 8)

    for (int y = 0; y < 7; ++y){
        INIT_LOOP_PDPC()

        COMPUTE_PLANAR_VECTOR2(1)

        APPLY_PDPC(1, 1)

        STORE1(_mm_storeu_si128, _dst, out_v, 8)
        _dst += dst_stride;
    }
    for (int y = 7; y < height; ++y){
        INIT_LOOP_NOPDPC()

        COMPUTE_PLANAR_VECTOR2(1)

        APPLY_NOPDPC(1, 1)

        STORE1(_mm_storeu_si128, _dst, out_v, 8)
        _dst += dst_stride;
    }
}

void
vvc_intra_planar_pdpc_w16_h4_sse(const uint16_t *const src_above,
                          const uint16_t *const src_left,
                          uint16_t *const dst, ptrdiff_t dst_stride,
                          int log2_pb_w, int log2_pb_h)
{
    uint16_t *_dst = dst;
    const uint32_t width  = 1 << log2_pb_w;
    const uint32_t height = 1 << log2_pb_h;
    int max_scale = 4;
    int rc_scale = 0;
    int tr_scale = 2;

    int16_t r_col[4], l_col[4];

    PLANAR_VECTOR_DEF_SCALE01(1, 16)

    LOAD2(__m128i a, _mm_loadu_si128, src_above + 1, 8)
    UNROLL2(U110, __m128i tr_v, _mm_slli_epi16, a, log2_pb_h)

    LOAD1(__m128i l, _mm_loadl_epi64, src_left + 1, 8)
    UNROLL1(U101, __m128i rc_v, _mm_sub_epi16, tr_val_v, l)
    UNROLL1(U110, __m128i lc_v, _mm_slli_epi16, l, log2_pb_w)

    STORE1(_mm_storel_epi64, (r_col),rc_v, 8)
    STORE1(_mm_storel_epi64, (l_col),lc_v, 8)

    for (int y = 0; y < height; ++y){
        INIT_LOOP_PDPC()

        COMPUTE_PLANAR_VECTOR(2)

        APPLY_PDPC(2, 1)

        STORE2(_mm_storeu_si128, _dst, out_v, 8)

        _dst += dst_stride;
    }
}

void
vvc_intra_planar_pdpc_w16_h8_sse(const uint16_t *const src_above,
                          const uint16_t *const src_left,
                          uint16_t *const dst, ptrdiff_t dst_stride,
                          int log2_pb_w, int log2_pb_h)
{
    uint16_t *_dst = dst;
    const uint32_t width  = 1 << log2_pb_w;
    const uint32_t height = 1 << log2_pb_h;
    int max_scale = 4;
    int rc_scale = 0;
    int tr_scale = 1;

    int16_t r_col[8], l_col[8];

    PLANAR_VECTOR_DEF_SCALE01(1, 16)

    LOAD2(__m128i a, _mm_loadu_si128, src_above + 1, 8)
    UNROLL2(U110, __m128i tr_v, _mm_slli_epi16, a, log2_pb_h)

    LOAD1(__m128i l, _mm_loadu_si128, src_left + 1, 8)
    UNROLL1(U101, __m128i rc_v, _mm_sub_epi16, tr_val_v, l)
    UNROLL1(U110, __m128i lc_v, _mm_slli_epi16, l, log2_pb_w)

    STORE1(_mm_storeu_si128, (r_col),rc_v, 8)
    STORE1(_mm_storeu_si128, (l_col),lc_v, 8)

    for (int y = 0; y < 7; ++y){
        INIT_LOOP_PDPC()

        COMPUTE_PLANAR_VECTOR(2)

        APPLY_PDPC(2, 1)

        STORE2(_mm_storeu_si128, _dst, out_v, 8)

        _dst += dst_stride;
    }
    for (int y = 7; y < height; ++y){
        INIT_LOOP_NOPDPC()

        COMPUTE_PLANAR_VECTOR(2)

        APPLY_NOPDPC(2, 1)

        STORE2(_mm_storeu_si128, _dst, out_v, 8)
        _dst += dst_stride;
    }
}

void
vvc_intra_planar_pdpc_w16_h16_sse(const uint16_t *const src_above,
                          const uint16_t *const src_left,
                          uint16_t *const dst, ptrdiff_t dst_stride,
                          int log2_pb_w, int log2_pb_h)
{
    uint16_t *_dst = dst;
    const uint32_t width  = 1 << log2_pb_w;
    const uint32_t height = 1 << log2_pb_h;
    int max_scale = 4;
    int rc_scale = 0;
    int tr_scale = 0;

    int16_t r_col[16], l_col[16];

    PLANAR_VECTOR_DEF_SCALE01(1, 16)

    LOAD2(__m128i a, _mm_loadu_si128, src_above + 1, 8)
    UNROLL2(U110, __m128i tr_v, _mm_slli_epi16, a, log2_pb_h)

    LOAD2(__m128i l, _mm_loadu_si128, src_left + 1, 8)
    UNROLL2(U101, __m128i rc_v, _mm_sub_epi16, tr_val_v, l)
    UNROLL2(U110, __m128i lc_v, _mm_slli_epi16, l, log2_pb_w)

    STORE2(_mm_storeu_si128, (r_col),rc_v, 8)
    STORE2(_mm_storeu_si128, (l_col),lc_v, 8)

    for (int y = 0; y < 7; ++y){
        INIT_LOOP_PDPC()

        COMPUTE_PLANAR_VECTOR(2)

        APPLY_PDPC(2, 1)

        STORE2(_mm_storeu_si128, _dst, out_v, 8)

        _dst += dst_stride;
    }
    for (int y = 7; y < height; ++y){
        INIT_LOOP_NOPDPC()

        COMPUTE_PLANAR_VECTOR(2)

        APPLY_NOPDPC(2, 1)

        STORE2(_mm_storeu_si128, _dst, out_v, 8)
        _dst += dst_stride;
    }
}

void
vvc_intra_planar_pdpc_w16_h32_sse(const uint16_t *const src_above,
                          const uint16_t *const src_left,
                          uint16_t *const dst, ptrdiff_t dst_stride,
                          int log2_pb_w, int log2_pb_h)
{
    uint16_t *_dst = dst;
    const uint32_t width  = 1 << log2_pb_w;
    const uint32_t height = 1 << log2_pb_h;
    int max_scale = 5;
    int rc_scale = 1;
    int tr_scale = 0;

    int16_t r_col[32], l_col[32];

    PLANAR_VECTOR_DEF_SCALE01(1, 16)

    LOAD2(__m128i a, _mm_loadu_si128, src_above + 1, 8)
    UNROLL2(U110, __m128i tr_v, _mm_slli_epi16, a, log2_pb_h)

    LOAD4(__m128i l, _mm_loadu_si128, src_left + 1, 8)
    UNROLL4(U101, __m128i rc_v, _mm_sub_epi16, tr_val_v, l)
    UNROLL4(U110, __m128i lc_v, _mm_slli_epi16, l, log2_pb_w)

    STORE4(_mm_storeu_si128, (r_col),rc_v, 8)
    STORE4(_mm_storeu_si128, (l_col),lc_v, 8)

    for (int y = 0; y < 7; ++y){
        INIT_LOOP_PDPC()

        COMPUTE_PLANAR_VECTOR(2)

        APPLY_PDPC(2, 1)

        STORE2(_mm_storeu_si128, _dst, out_v, 8)

        _dst += dst_stride;
    }
    for (int y = 7; y < height; ++y){
        INIT_LOOP_NOPDPC()

        COMPUTE_PLANAR_VECTOR(2)

        APPLY_NOPDPC(2, 1)

        STORE2(_mm_storeu_si128, _dst, out_v, 8)
        _dst += dst_stride;
    }
}

void
vvc_intra_planar_pdpc_w16_h64_sse(const uint16_t *const src_above,
                          const uint16_t *const src_left,
                          uint16_t *const dst, ptrdiff_t dst_stride,
                          int log2_pb_w, int log2_pb_h)
{   
    uint16_t *_dst = dst;
    const uint32_t width  = 1 << log2_pb_w;
    const uint32_t height = 1 << log2_pb_h;
    int max_scale = 6;
    int rc_scale = 2;
    int tr_scale = 0;

    int16_t r_col[64], l_col[64];

    PLANAR_VECTOR_DEF_SCALE2(32)

    LOAD2(__m128i a, _mm_loadu_si128, src_above + 1, 8)
    UNROLL2(U110, __m128i tr_v, _mm_slli_epi16, a, log2_pb_h)

    LOAD8(__m128i l, _mm_loadu_si128, src_left + 1, 8)
    UNROLL8(U101, __m128i rc_v, _mm_sub_epi16, tr_val_v, l)
    UNROLL8(U110, __m128i lc_v, _mm_slli_epi16, l, log2_pb_w)

    STORE8(_mm_storeu_si128, r_col, rc_v, 8)
    STORE8(_mm_storeu_si128, l_col, lc_v, 8)

    for (int y = 0; y < 12; ++y){
        INIT_LOOP_PDPC()

        COMPUTE_PLANAR_VECTOR2(2)

        APPLY_PDPC(2, 2)

        STORE2(_mm_storeu_si128, _dst, out_v, 8)
        _dst += dst_stride;
    }
    for (int y = 12; y < height; ++y){
        INIT_LOOP_NOPDPC()

        COMPUTE_PLANAR_VECTOR2(2)

        APPLY_NOPDPC(2, 2)

        STORE2(_mm_storeu_si128, _dst, out_v, 8)
        _dst += dst_stride;
    }
}

void
vvc_intra_planar_pdpc_w32_h4_sse(const uint16_t *const src_above,
                          const uint16_t *const src_left,
                          uint16_t *const dst, ptrdiff_t dst_stride,
                          int log2_pb_w, int log2_pb_h)
{
    uint16_t *_dst = dst;
    const uint32_t width  = 1 << log2_pb_w;
    const uint32_t height = 1 << log2_pb_h;
    int max_scale = 5;
    int rc_scale = 0;
    int tr_scale = 3;

    int16_t r_col[4], l_col[4];

    PLANAR_VECTOR_DEF_SCALE01(1, 16)

    LOAD4(__m128i a, _mm_loadu_si128, src_above + 1, 8)
    UNROLL4(U110, __m128i tr_v, _mm_slli_epi16, a, log2_pb_h)

    LOAD1(__m128i l, _mm_loadl_epi64, src_left + 1, 8)
    UNROLL1(U101, __m128i rc_v, _mm_sub_epi16, tr_val_v, l)
    UNROLL1(U110, __m128i lc_v, _mm_slli_epi16, l, log2_pb_w)

    STORE1(_mm_storel_epi64, (r_col),rc_v, 8)
    STORE1(_mm_storel_epi64, (l_col),lc_v, 8)

    for (int y = 0; y < height; ++y){
        INIT_LOOP_PDPC()

        COMPUTE_PLANAR_VECTOR(4)

        APPLY_PDPC(4, 1)

        STORE4(_mm_storeu_si128, _dst, out_v, 8)

        _dst += dst_stride;
    }
}

void
vvc_intra_planar_pdpc_w32_h8_sse(const uint16_t *const src_above,
                          const uint16_t *const src_left,
                          uint16_t *const dst, ptrdiff_t dst_stride,
                          int log2_pb_w, int log2_pb_h)
{
    uint16_t *_dst = dst;
    const uint32_t width  = 1 << log2_pb_w;
    const uint32_t height = 1 << log2_pb_h;
    int max_scale = 5;
    int rc_scale = 0;
    int tr_scale = 2;

    int16_t r_col[8], l_col[8];

    PLANAR_VECTOR_DEF_SCALE01(1, 16)

    LOAD4(__m128i a, _mm_loadu_si128, src_above + 1, 8)
    UNROLL4(U110, __m128i tr_v, _mm_slli_epi16, a, log2_pb_h)

    LOAD1(__m128i l, _mm_loadu_si128, src_left + 1, 8)
    UNROLL1(U101, __m128i rc_v, _mm_sub_epi16, tr_val_v, l)
    UNROLL1(U110, __m128i lc_v, _mm_slli_epi16, l, log2_pb_w)

    STORE1(_mm_storeu_si128, (r_col),rc_v, 8)
    STORE1(_mm_storeu_si128, (l_col),lc_v, 8)

    for (int y = 0; y < 7; ++y){
        INIT_LOOP_PDPC()

        COMPUTE_PLANAR_VECTOR(4)

        APPLY_PDPC(4, 1)

        STORE4(_mm_storeu_si128, _dst, out_v, 8)

        _dst += dst_stride;
    }
    for (int y = 7; y < height; ++y){
        INIT_LOOP_NOPDPC()

        COMPUTE_PLANAR_VECTOR(4)

        APPLY_NOPDPC(4, 1)

        STORE4(_mm_storeu_si128, _dst, out_v, 8)
        _dst += dst_stride;
    }
}

void
vvc_intra_planar_pdpc_w32_h16_sse(const uint16_t *const src_above,
                          const uint16_t *const src_left,
                          uint16_t *const dst, ptrdiff_t dst_stride,
                          int log2_pb_w, int log2_pb_h)
{
    uint16_t *_dst = dst;
    const uint32_t width  = 1 << log2_pb_w;
    const uint32_t height = 1 << log2_pb_h;
    int max_scale = 5;
    int rc_scale = 0;
    int tr_scale = 1;

    int16_t r_col[16], l_col[16];

    PLANAR_VECTOR_DEF_SCALE01(1, 16)

    LOAD4(__m128i a, _mm_loadu_si128, src_above + 1, 8)
    UNROLL4(U110, __m128i tr_v, _mm_slli_epi16, a, log2_pb_h)

    LOAD2(__m128i l, _mm_loadu_si128, src_left + 1, 8)
    UNROLL2(U101, __m128i rc_v, _mm_sub_epi16, tr_val_v, l)
    UNROLL2(U110, __m128i lc_v, _mm_slli_epi16, l, log2_pb_w)

    STORE2(_mm_storeu_si128, (r_col),rc_v, 8)
    STORE2(_mm_storeu_si128, (l_col),lc_v, 8)

    for (int y = 0; y < 7; ++y){
        INIT_LOOP_PDPC()

        COMPUTE_PLANAR_VECTOR(4)

        APPLY_PDPC(4, 1)

        STORE4(_mm_storeu_si128, _dst, out_v, 8)

        _dst += dst_stride;
    }
    for (int y = 7; y < height; ++y){
        INIT_LOOP_NOPDPC()

        COMPUTE_PLANAR_VECTOR(4)

        APPLY_NOPDPC(4, 1)

        STORE4(_mm_storeu_si128, _dst, out_v, 8)
        _dst += dst_stride;
    }
}

void
vvc_intra_planar_pdpc_w32_h32_sse(const uint16_t *const src_above,
                          const uint16_t *const src_left,
                          uint16_t *const dst, ptrdiff_t dst_stride,
                          int log2_pb_w, int log2_pb_h)
{
    uint16_t *_dst = dst;
    const uint32_t width  = 1 << log2_pb_w;
    const uint32_t height = 1 << log2_pb_h;
    int max_scale = 5;
    int rc_scale = 0;
    int tr_scale = 0;

    int16_t r_col[32], l_col[32];

    PLANAR_VECTOR_DEF_SCALE2(16);

    LOAD4(__m128i a, _mm_loadu_si128, src_above + 1, 8)
    UNROLL4(U110, __m128i tr_v, _mm_slli_epi16, a, log2_pb_h)

    LOAD4(__m128i l, _mm_loadu_si128, src_left + 1, 8)
    UNROLL4(U101, __m128i rc_v, _mm_sub_epi16, tr_val_v, l)
    UNROLL4(U110, __m128i lc_v, _mm_slli_epi16, l, log2_pb_w)

    STORE4(_mm_storeu_si128, (r_col),rc_v, 8)
    STORE4(_mm_storeu_si128, (l_col),lc_v, 8)

    for (int y = 0; y < 12; ++y){
        INIT_LOOP_PDPC()

        COMPUTE_PLANAR_VECTOR(4)

        APPLY_PDPC(4, 2)

        STORE4(_mm_storeu_si128, _dst, out_v, 8)

        _dst += dst_stride;
    }
    for (int y = 12; y < height; ++y){
        INIT_LOOP_NOPDPC()

        COMPUTE_PLANAR_VECTOR(4)

        APPLY_NOPDPC(4, 2)

        STORE4(_mm_storeu_si128, _dst, out_v, 8)
        _dst += dst_stride;
    }
}

void
vvc_intra_planar_pdpc_w32_h64_sse(const uint16_t *const src_above,
                          const uint16_t *const src_left,
                          uint16_t *const dst, ptrdiff_t dst_stride,
                          int log2_pb_w, int log2_pb_h)
{   
    uint16_t *_dst = dst;
    const uint32_t width  = 1 << log2_pb_w;
    const uint32_t height = 1 << log2_pb_h;
    int max_scale = 6;
    int rc_scale = 1;
    int tr_scale = 0;

    int16_t r_col[64], l_col[64];

    PLANAR_VECTOR_DEF_SCALE2(32)

    LOAD4(__m128i a, _mm_loadu_si128, src_above + 1, 8)
    UNROLL4(U110, __m128i tr_v, _mm_slli_epi16, a, log2_pb_h)

    LOAD8(__m128i l, _mm_loadu_si128, src_left + 1, 8)
    UNROLL8(U101, __m128i rc_v, _mm_sub_epi16, tr_val_v, l)
    UNROLL8(U110, __m128i lc_v, _mm_slli_epi16, l, log2_pb_w)

    STORE8(_mm_storeu_si128, r_col, rc_v, 8)
    STORE8(_mm_storeu_si128, l_col, lc_v, 8)

    for (int y = 0; y < 12; ++y){
        INIT_LOOP_PDPC()

        COMPUTE_PLANAR_VECTOR2(4)

        APPLY_PDPC(4, 2)

        STORE4(_mm_storeu_si128, _dst, out_v, 8)
        _dst += dst_stride;
    }
    for (int y = 12; y < height; ++y){
        INIT_LOOP_NOPDPC()

        COMPUTE_PLANAR_VECTOR2(4)

        APPLY_NOPDPC(4, 2)

        STORE4(_mm_storeu_si128, _dst, out_v, 8)
        _dst += dst_stride;
    }
}

void
vvc_intra_planar_pdpc_w64_h4_sse(const uint16_t *const src_above,
                          const uint16_t *const src_left,
                          uint16_t *const dst, ptrdiff_t dst_stride,
                          int log2_pb_w, int log2_pb_h)
{
    uint16_t *_dst = dst;
    const uint32_t width  = 1 << log2_pb_w;
    const uint32_t height = 1 << log2_pb_h;
    int max_scale = 6;
    int rc_scale = 0;
    int tr_scale = 4;

    int16_t r_col[4], l_col[4];

    PLANAR_VECTOR_DEF_SCALE01(1, 32)

    LOAD8(__m128i a, _mm_loadu_si128, src_above + 1, 8)
    UNROLL8(U110, __m128i tr_v, _mm_slli_epi16, a, log2_pb_h)

    LOAD1(__m128i l, _mm_loadl_epi64, src_left + 1, 8)
    UNROLL1(U101, __m128i rc_v, _mm_sub_epi16, tr_val_v, l)
    UNROLL1(U110, __m128i lc_v, _mm_slli_epi16, l, log2_pb_w)

    STORE1(_mm_storel_epi64, (r_col),rc_v, 8)
    STORE1(_mm_storel_epi64, (l_col),lc_v, 8)

    for (int y = 0; y < height; ++y){
        INIT_LOOP_PDPC()

        COMPUTE_PLANAR_VECTOR2(8)

        APPLY_PDPC(8, 1)

        STORE8(_mm_storeu_si128, _dst, out_v, 8)

        _dst += dst_stride;
    }
}

void
vvc_intra_planar_pdpc_w64_h8_sse(const uint16_t *const src_above,
                          const uint16_t *const src_left,
                          uint16_t *const dst, ptrdiff_t dst_stride,
                          int log2_pb_w, int log2_pb_h)
{
    uint16_t *_dst = dst;
    const uint32_t width  = 1 << log2_pb_w;
    const uint32_t height = 1 << log2_pb_h;
    int max_scale = 6;
    int rc_scale = 0;
    int tr_scale = 3;

    int16_t r_col[8], l_col[8];

    PLANAR_VECTOR_DEF_SCALE01(1, 32)

    LOAD8(__m128i a, _mm_loadu_si128, src_above + 1, 8)
    UNROLL8(U110, __m128i tr_v, _mm_slli_epi16, a, log2_pb_h)

    LOAD1(__m128i l, _mm_loadu_si128, src_left + 1, 8)
    UNROLL1(U101, __m128i rc_v, _mm_sub_epi16, tr_val_v, l)
    UNROLL1(U110, __m128i lc_v, _mm_slli_epi16, l, log2_pb_w)

    STORE1(_mm_storeu_si128, (r_col),rc_v, 8)
    STORE1(_mm_storeu_si128, (l_col),lc_v, 8)

    for (int y = 0; y < 7; ++y){
        INIT_LOOP_PDPC()

        COMPUTE_PLANAR_VECTOR2(8)

        APPLY_PDPC(8, 1)

        STORE8(_mm_storeu_si128, _dst, out_v, 8)

        _dst += dst_stride;
    }
    for (int y = 7; y < height; ++y){
        INIT_LOOP_NOPDPC()

        COMPUTE_PLANAR_VECTOR2(8)

        APPLY_NOPDPC(8, 1)

        STORE8(_mm_storeu_si128, _dst, out_v, 8)
        _dst += dst_stride;
    }
}

void
vvc_intra_planar_pdpc_w64_h16_sse(const uint16_t *const src_above,
                          const uint16_t *const src_left,
                          uint16_t *const dst, ptrdiff_t dst_stride,
                          int log2_pb_w, int log2_pb_h)
{
    uint16_t *_dst = dst;
    const uint32_t width  = 1 << log2_pb_w;
    const uint32_t height = 1 << log2_pb_h;
    int max_scale = 6;
    int rc_scale = 0;
    int tr_scale = 2;

    int16_t r_col[16], l_col[16];

    PLANAR_VECTOR_DEF_SCALE2(32)

    LOAD8(__m128i a, _mm_loadu_si128, src_above + 1, 8)
    UNROLL8(U110, __m128i tr_v, _mm_slli_epi16, a, log2_pb_h)

    LOAD2(__m128i l, _mm_loadu_si128, src_left + 1, 8)
    UNROLL2(U101, __m128i rc_v, _mm_sub_epi16, tr_val_v, l)
    UNROLL2(U110, __m128i lc_v, _mm_slli_epi16, l, log2_pb_w)

    STORE2(_mm_storeu_si128, (r_col),rc_v, 8)
    STORE2(_mm_storeu_si128, (l_col),lc_v, 8)

    for (int y = 0; y < 12; ++y){
        INIT_LOOP_PDPC()

        COMPUTE_PLANAR_VECTOR2(8)

        APPLY_PDPC(8, 2)

        STORE8(_mm_storeu_si128, _dst, out_v, 8)

        _dst += dst_stride;
    }
    for (int y = 12; y < height; ++y){
        INIT_LOOP_NOPDPC()

        COMPUTE_PLANAR_VECTOR2(8)

        APPLY_NOPDPC(8, 2)

        STORE8(_mm_storeu_si128, _dst, out_v, 8)
        _dst += dst_stride;
    }
}

void
vvc_intra_planar_pdpc_w64_h32_sse(const uint16_t *const src_above,
                          const uint16_t *const src_left,
                          uint16_t *const dst, ptrdiff_t dst_stride,
                          int log2_pb_w, int log2_pb_h)
{
    uint16_t *_dst = dst;
    const uint32_t width  = 1 << log2_pb_w;
    const uint32_t height = 1 << log2_pb_h;
    int max_scale = 6;
    int rc_scale = 0;
    int tr_scale = 1;

    int16_t r_col[32], l_col[32];

    PLANAR_VECTOR_DEF_SCALE2(32);

    LOAD8(__m128i a, _mm_loadu_si128, src_above + 1, 8)
    UNROLL8(U110, __m128i tr_v, _mm_slli_epi16, a, log2_pb_h)

    LOAD4(__m128i l, _mm_loadu_si128, src_left + 1, 8)
    UNROLL4(U101, __m128i rc_v, _mm_sub_epi16, tr_val_v, l)
    UNROLL4(U110, __m128i lc_v, _mm_slli_epi16, l, log2_pb_w)

    STORE4(_mm_storeu_si128, (r_col),rc_v, 8)
    STORE4(_mm_storeu_si128, (l_col),lc_v, 8)

    for (int y = 0; y < 12; ++y){
        INIT_LOOP_PDPC()

        COMPUTE_PLANAR_VECTOR2(8)

        APPLY_PDPC(8, 2)

        STORE8(_mm_storeu_si128, _dst, out_v, 8)

        _dst += dst_stride;
    }
    for (int y = 12; y < height; ++y){
        INIT_LOOP_NOPDPC()

        COMPUTE_PLANAR_VECTOR2(8)

        APPLY_NOPDPC(8, 2)

        STORE8(_mm_storeu_si128, _dst, out_v, 8)
        _dst += dst_stride;
    }
}

void
vvc_intra_planar_pdpc_w64_h64_sse(const uint16_t *const src_above,
                          const uint16_t *const src_left,
                          uint16_t *const dst, ptrdiff_t dst_stride,
                          int log2_pb_w, int log2_pb_h)
{   
    uint16_t *_dst = dst;
    const uint32_t width  = 1 << log2_pb_w;
    const uint32_t height = 1 << log2_pb_h;
    int max_scale = 6;
    int rc_scale = 0;
    int tr_scale = 0;

    int16_t r_col[64], l_col[64];

    PLANAR_VECTOR_DEF_SCALE2(32)

    LOAD8(__m128i a, _mm_loadu_si128, src_above + 1, 8)
    UNROLL8(U110, __m128i tr_v, _mm_slli_epi16, a, log2_pb_h)

    LOAD8(__m128i l, _mm_loadu_si128, src_left + 1, 8)
    UNROLL8(U101, __m128i rc_v, _mm_sub_epi16, tr_val_v, l)
    UNROLL8(U110, __m128i lc_v, _mm_slli_epi16, l, log2_pb_w)

    STORE8(_mm_storeu_si128, r_col, rc_v, 8)
    STORE8(_mm_storeu_si128, l_col, lc_v, 8)

    for (int y = 0; y < 12; ++y){
        INIT_LOOP_PDPC()

        COMPUTE_PLANAR_VECTOR2(8)

        APPLY_PDPC(8, 2)

        STORE8(_mm_storeu_si128, _dst, out_v, 8)
        _dst += dst_stride;
    }
    for (int y = 12; y < height; ++y){
        INIT_LOOP_NOPDPC()

        COMPUTE_PLANAR_VECTOR2(8)

        APPLY_NOPDPC(8, 2)

        STORE8(_mm_storeu_si128, _dst, out_v, 8)
        _dst += dst_stride;
    }
}

void rcn_init_dc_planar_functions_sse(struct RCNFunctions *const rcn_funcs){
    //DC
    rcn_funcs->dc.pdpc[0][0] = &vvc_intra_dc_pdpc_w4_h4_sse;
    rcn_funcs->dc.pdpc[0][1] = &vvc_intra_dc_pdpc_w4_h8_sse;
    rcn_funcs->dc.pdpc[0][2] = &vvc_intra_dc_pdpc_w4_h16_sse;
    rcn_funcs->dc.pdpc[0][3] = &vvc_intra_dc_pdpc_w4_h32_sse;
    rcn_funcs->dc.pdpc[0][4] = &vvc_intra_dc_pdpc_w4_h64_sse;

    rcn_funcs->dc.pdpc[1][0] = &vvc_intra_dc_pdpc_w8_h4_sse;
    rcn_funcs->dc.pdpc[1][1] = &vvc_intra_dc_pdpc_w8_h8_sse;
    rcn_funcs->dc.pdpc[1][2] = &vvc_intra_dc_pdpc_w8_h16_sse;
    rcn_funcs->dc.pdpc[1][3] = &vvc_intra_dc_pdpc_w8_h32_sse;
    rcn_funcs->dc.pdpc[1][4] = &vvc_intra_dc_pdpc_w8_h64_sse;

    rcn_funcs->dc.pdpc[2][0] = &vvc_intra_dc_pdpc_w16_h4_sse;
    rcn_funcs->dc.pdpc[2][1] = &vvc_intra_dc_pdpc_w16_h8_sse;
    rcn_funcs->dc.pdpc[2][2] = &vvc_intra_dc_pdpc_w16_h16_sse;
    rcn_funcs->dc.pdpc[2][3] = &vvc_intra_dc_pdpc_w16_h32_sse;
    rcn_funcs->dc.pdpc[2][4] = &vvc_intra_dc_pdpc_w16_h64_sse;

    rcn_funcs->dc.pdpc[3][0] = &vvc_intra_dc_pdpc_w32_h4_sse;
    rcn_funcs->dc.pdpc[3][1] = &vvc_intra_dc_pdpc_w32_h8_sse;
    rcn_funcs->dc.pdpc[3][2] = &vvc_intra_dc_pdpc_w32_h8_sse;//identical to h8
    rcn_funcs->dc.pdpc[3][3] = &vvc_intra_dc_pdpc_w32_h32_sse;
    rcn_funcs->dc.pdpc[3][4] = &vvc_intra_dc_pdpc_w32_h64_sse;

    rcn_funcs->dc.pdpc[4][0] = &vvc_intra_dc_pdpc_w64_h4_sse;
    rcn_funcs->dc.pdpc[4][1] = &vvc_intra_dc_pdpc_w64_h8_sse;
    rcn_funcs->dc.pdpc[4][2] = &vvc_intra_dc_pdpc_w64_h16_sse;
    rcn_funcs->dc.pdpc[4][3] = &vvc_intra_dc_pdpc_w64_h16_sse;//identical to h16
    rcn_funcs->dc.pdpc[4][4] = &vvc_intra_dc_pdpc_w64_h64_sse;

    //Planar
    rcn_funcs->planar.pdpc[0][0] = &vvc_intra_planar_pdpc_w4_h4_sse;
    rcn_funcs->planar.pdpc[0][1] = &vvc_intra_planar_pdpc_w4_h8_sse;
    rcn_funcs->planar.pdpc[0][2] = &vvc_intra_planar_pdpc_w4_h16_sse;
    rcn_funcs->planar.pdpc[0][3] = &vvc_intra_planar_pdpc_w4_h32_sse;
    rcn_funcs->planar.pdpc[0][4] = &vvc_intra_planar_pdpc_w4_h64_sse;

    rcn_funcs->planar.pdpc[1][0] = &vvc_intra_planar_pdpc_w8_h4_sse;
    rcn_funcs->planar.pdpc[1][1] = &vvc_intra_planar_pdpc_w8_h8_sse;
    rcn_funcs->planar.pdpc[1][2] = &vvc_intra_planar_pdpc_w8_h16_sse;
    rcn_funcs->planar.pdpc[1][3] = &vvc_intra_planar_pdpc_w8_h32_sse;
    rcn_funcs->planar.pdpc[1][4] = &vvc_intra_planar_pdpc_w8_h64_sse;

    rcn_funcs->planar.pdpc[2][0] = &vvc_intra_planar_pdpc_w16_h4_sse;
    rcn_funcs->planar.pdpc[2][1] = &vvc_intra_planar_pdpc_w16_h8_sse;
    rcn_funcs->planar.pdpc[2][2] = &vvc_intra_planar_pdpc_w16_h16_sse;
    rcn_funcs->planar.pdpc[2][3] = &vvc_intra_planar_pdpc_w16_h32_sse;
    rcn_funcs->planar.pdpc[2][4] = &vvc_intra_planar_pdpc_w16_h64_sse;

    rcn_funcs->planar.pdpc[3][0] = &vvc_intra_planar_pdpc_w32_h4_sse;
    rcn_funcs->planar.pdpc[3][1] = &vvc_intra_planar_pdpc_w32_h8_sse;
    rcn_funcs->planar.pdpc[3][2] = &vvc_intra_planar_pdpc_w32_h16_sse;
    rcn_funcs->planar.pdpc[3][3] = &vvc_intra_planar_pdpc_w32_h32_sse;
    rcn_funcs->planar.pdpc[3][4] = &vvc_intra_planar_pdpc_w32_h64_sse;

    rcn_funcs->planar.pdpc[4][0] = &vvc_intra_planar_pdpc_w64_h4_sse;
    rcn_funcs->planar.pdpc[4][1] = &vvc_intra_planar_pdpc_w64_h8_sse;
    rcn_funcs->planar.pdpc[4][2] = &vvc_intra_planar_pdpc_w64_h16_sse;
    rcn_funcs->planar.pdpc[4][3] = &vvc_intra_planar_pdpc_w64_h32_sse;
    rcn_funcs->planar.pdpc[4][4] = &vvc_intra_planar_pdpc_w64_h64_sse;
}
