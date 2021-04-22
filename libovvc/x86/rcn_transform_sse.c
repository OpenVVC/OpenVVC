#include <emmintrin.h>
#include <smmintrin.h>
#include <getopt.h>
#include <stdio.h>
#include <stddef.h>

#include "ovutils.h"
#include "data_rcn_transform.h"
#include "rcn_transform.h"
#include "x86/vvc_utils_sse.h"

#if 0
static inline void
multiply(__m128i *x, __m128i *y, __m128i *r)
{
    __m128i xylo = _mm_mullo_epi16(*x, *y);
    __m128i xyhi = _mm_mulhi_epi16(*x, *y);
    r[0] = _mm_unpacklo_epi16(xylo, xyhi);
    r[1] = _mm_unpackhi_epi16(xylo, xyhi);
}
#endif

static inline void
transpose4x4(__m128i *x, __m128i *r)
{
    __m128i tmp[4];

    tmp[0] = _mm_unpacklo_epi32(x[0], x[1]);
    tmp[1] = _mm_unpackhi_epi32(x[0], x[1]);
    tmp[2] = _mm_unpacklo_epi32(x[2], x[3]);
    tmp[3] = _mm_unpackhi_epi32(x[2], x[3]);

    r[0] = _mm_unpacklo_epi64(tmp[0], tmp[2]);
    r[1] = _mm_unpackhi_epi64(tmp[0], tmp[2]);
    r[2] = _mm_unpacklo_epi64(tmp[1], tmp[3]);
    r[3] = _mm_unpackhi_epi64(tmp[1], tmp[3]);
}

static inline void
transpose4x4_step2(__m128i *x, __m128i *r)
{
    __m128i tmp[4];

    tmp[0] = _mm_unpacklo_epi32(x[0], x[2]);
    tmp[1] = _mm_unpackhi_epi32(x[0], x[2]);
    tmp[2] = _mm_unpacklo_epi32(x[4], x[6]);
    tmp[3] = _mm_unpackhi_epi32(x[4], x[6]);

    r[0] = _mm_unpacklo_epi64(tmp[0], tmp[2]);
    r[1] = _mm_unpackhi_epi64(tmp[0], tmp[2]);
    r[2] = _mm_unpacklo_epi64(tmp[1], tmp[3]);
    r[3] = _mm_unpackhi_epi64(tmp[1], tmp[3]);
}

static inline void
matMult4x4_red1(__m128i *x, __m128i *d, __m128i *r)
{
    __m128i m[8], z;//a[4],
    z = _mm_setzero_si128();

    m[0] = _mm_unpacklo_epi16(x[0], z);
    m[1] = _mm_unpackhi_epi16(x[0], z);
    m[2] = _mm_unpacklo_epi16(d[0], z);
    m[3] = _mm_unpackhi_epi16(d[0], z);

    r[0] = _mm_madd_epi16(m[0], m[2]);
    r[1] = _mm_madd_epi16(m[1], m[3]);
}

static inline void
matMult4x4_red2(__m128i *x, __m128i *d, __m128i *r)
{
    __m128i m[8];//, a[4];

    m[0] = _mm_unpacklo_epi16(x[0], x[1]);
    m[1] = _mm_unpackhi_epi16(x[0], x[1]);
    m[2] = _mm_unpacklo_epi16(d[0], d[1]);
    m[3] = _mm_unpackhi_epi16(d[0], d[1]);

    r[0] = _mm_madd_epi16(m[0], m[2]);
    r[1] = _mm_madd_epi16(m[1], m[3]);
}

static inline void
matMult4x4(__m128i *x, __m128i *d, __m128i *r)
{
    __m128i m[8], a[4];

    #if 0
    multiply(x + 0, d + 0, m + 0);
    multiply(x + 1, d + 1, m + 2);
    multiply(x + 2, d + 2, m + 4);
    multiply(x + 3, d + 3, m + 6);

    a[0] = _mm_add_epi32(m[0], m[2]);
    a[1] = _mm_add_epi32(m[1], m[3]);
    a[2] = _mm_add_epi32(m[4], m[6]);
    a[3] = _mm_add_epi32(m[5], m[7]);
    #else
    m[0] = _mm_unpacklo_epi16(x[0], x[1]);
    m[1] = _mm_unpackhi_epi16(x[0], x[1]);
    m[2] = _mm_unpacklo_epi16(d[0], d[1]);
    m[3] = _mm_unpackhi_epi16(d[0], d[1]);
    m[4] = _mm_unpacklo_epi16(x[2], x[3]);
    m[5] = _mm_unpackhi_epi16(x[2], x[3]);
    m[6] = _mm_unpacklo_epi16(d[2], d[3]);
    m[7] = _mm_unpackhi_epi16(d[2], d[3]);

    a[0] = _mm_madd_epi16(m[0], m[2]);
    a[1] = _mm_madd_epi16(m[1], m[3]);
    a[2] = _mm_madd_epi16(m[4], m[6]);
    a[3] = _mm_madd_epi16(m[5], m[7]);
    #endif

    r[0] = _mm_add_epi32(a[0], a[2]);
    r[1] = _mm_add_epi32(a[1], a[3]);
}

static inline void
matMult8x8_red1(__m128i *x, __m128i *d, __m128i *r)
{
    __m128i m[16], z;//a[8], b[4],
    z = _mm_setzero_si128();

    m[ 0] = _mm_unpacklo_epi16(x[0], z);
    m[ 1] = _mm_unpackhi_epi16(x[0], z);
    m[ 2] = _mm_unpacklo_epi16(d[0], z);
    m[ 3] = _mm_unpackhi_epi16(d[0], z);

    r[0] = _mm_madd_epi16(m[ 0], m[ 2]);
    r[1] = _mm_madd_epi16(m[ 1], m[ 3]);
}

static inline void
matMult8x8_red2(__m128i *x, __m128i *d, __m128i *r)
{
    __m128i m[16];//, a[8], b[4];

    m[ 0] = _mm_unpacklo_epi16(x[0], x[1]);
    m[ 1] = _mm_unpackhi_epi16(x[0], x[1]);
    m[ 2] = _mm_unpacklo_epi16(d[0], d[1]);
    m[ 3] = _mm_unpackhi_epi16(d[0], d[1]);

    r[0] = _mm_madd_epi16(m[ 0], m[ 2]);
    r[1] = _mm_madd_epi16(m[ 1], m[ 3]);
}

static inline void
matMult8x8_red4(__m128i *x, __m128i *d, __m128i *r)
{
    __m128i m[16], a[8];//, b[4];

    m[ 0] = _mm_unpacklo_epi16(x[0], x[1]);
    m[ 1] = _mm_unpackhi_epi16(x[0], x[1]);
    m[ 2] = _mm_unpacklo_epi16(d[0], d[1]);
    m[ 3] = _mm_unpackhi_epi16(d[0], d[1]);
    m[ 4] = _mm_unpacklo_epi16(x[2], x[3]);
    m[ 5] = _mm_unpackhi_epi16(x[2], x[3]);
    m[ 6] = _mm_unpacklo_epi16(d[2], d[3]);
    m[ 7] = _mm_unpackhi_epi16(d[2], d[3]);

    a[0] = _mm_madd_epi16(m[ 0], m[ 2]);
    a[1] = _mm_madd_epi16(m[ 1], m[ 3]);
    a[2] = _mm_madd_epi16(m[ 4], m[ 6]);
    a[3] = _mm_madd_epi16(m[ 5], m[ 7]);

    r[0] = _mm_add_epi32(a[0], a[2]);
    r[1] = _mm_add_epi32(a[1], a[3]);
}

static inline void
matMult8x8(__m128i *x, __m128i *d, __m128i *r)
{
    __m128i m[16], a[8], b[4];

    #if 0
    multiply(x + 0, d + 0, m + 0);
    multiply(x + 1, d + 1, m + 2);
    multiply(x + 2, d + 2, m + 4);
    multiply(x + 3, d + 3, m + 6);
    multiply(x + 4, d + 4, m + 8);
    multiply(x + 5, d + 5, m + 10);
    multiply(x + 6, d + 6, m + 12);
    multiply(x + 7, d + 7, m + 14);

    a[0] = _mm_add_epi32(m[0], m[2]);
    a[1] = _mm_add_epi32(m[1], m[3]);
    a[2] = _mm_add_epi32(m[4], m[6]);
    a[3] = _mm_add_epi32(m[5], m[7]);
    a[4] = _mm_add_epi32(m[8], m[10]);
    a[5] = _mm_add_epi32(m[9], m[11]);
    a[6] = _mm_add_epi32(m[12], m[14]);
    a[7] = _mm_add_epi32(m[13], m[15]);

    b[0] = _mm_add_epi32(a[0], a[2]);
    b[1] = _mm_add_epi32(a[1], a[3]);
    b[2] = _mm_add_epi32(a[4], a[6]);
    b[3] = _mm_add_epi32(a[5], a[7]);
    #else
    m[ 0] = _mm_unpacklo_epi16(x[0], x[1]);
    m[ 1] = _mm_unpackhi_epi16(x[0], x[1]);
    m[ 2] = _mm_unpacklo_epi16(d[0], d[1]);
    m[ 3] = _mm_unpackhi_epi16(d[0], d[1]);
    m[ 4] = _mm_unpacklo_epi16(x[2], x[3]);
    m[ 5] = _mm_unpackhi_epi16(x[2], x[3]);
    m[ 6] = _mm_unpacklo_epi16(d[2], d[3]);
    m[ 7] = _mm_unpackhi_epi16(d[2], d[3]);

    a[0] = _mm_madd_epi16(m[ 0], m[ 2]);
    a[1] = _mm_madd_epi16(m[ 1], m[ 3]);
    a[2] = _mm_madd_epi16(m[ 4], m[ 6]);
    a[3] = _mm_madd_epi16(m[ 5], m[ 7]);

    b[0] = _mm_add_epi32(a[0], a[2]);
    b[1] = _mm_add_epi32(a[1], a[3]);

    m[ 8] = _mm_unpacklo_epi16(x[4], x[5]);
    m[ 9] = _mm_unpackhi_epi16(x[4], x[5]);
    m[10] = _mm_unpacklo_epi16(d[4], d[5]);
    m[11] = _mm_unpackhi_epi16(d[4], d[5]);
    m[12] = _mm_unpacklo_epi16(x[6], x[7]);
    m[13] = _mm_unpackhi_epi16(x[6], x[7]);
    m[14] = _mm_unpacklo_epi16(d[6], d[7]);
    m[15] = _mm_unpackhi_epi16(d[6], d[7]);

    a[4] = _mm_madd_epi16(m[ 8], m[10]);
    a[5] = _mm_madd_epi16(m[ 9], m[11]);
    a[6] = _mm_madd_epi16(m[12], m[14]);
    a[7] = _mm_madd_epi16(m[13], m[15]);

    b[2] = _mm_add_epi32(a[4], a[6]);
    b[3] = _mm_add_epi32(a[5], a[7]);
    #endif

    r[0] = _mm_add_epi32(b[0], b[2]);
    r[1] = _mm_add_epi32(b[1], b[3]);
}

static inline void
matMult16x16(__m128i *x, __m128i *d, __m128i *r)
{
    __m128i m[32], a[16], b[8], c[4];

    #if 0
    for (int k=0; k<16; k++){
        multiply(x + k, d + k, m + 2*k);
    }
    for (int k=0; k<8; k++){
        a[2*k+0] = _mm_add_epi32(m[4*k+0], m[4*k+2]);
        a[2*k+1] = _mm_add_epi32(m[4*k+1], m[4*k+3]);
    }
    for (int k=0; k<4; k++){
        b[2*k+0] = _mm_add_epi32(a[4*k+0], a[4*k+2]);
        b[2*k+1] = _mm_add_epi32(a[4*k+1], a[4*k+3]);
    }
    for (int k=0; k<2; k++){
        c[2*k+0] = _mm_add_epi32(b[4*k+0], b[4*k+2]);
        c[2*k+1] = _mm_add_epi32(b[4*k+1], b[4*k+3]);
    }
    #else
    m[ 0] = _mm_unpacklo_epi16(x[ 0], x[1]);
    m[ 1] = _mm_unpackhi_epi16(x[ 0], x[1]);
    m[ 2] = _mm_unpacklo_epi16(d[ 0], d[1]);
    m[ 3] = _mm_unpackhi_epi16(d[ 0], d[1]);
    m[ 4] = _mm_unpacklo_epi16(x[ 2], x[3]);
    m[ 5] = _mm_unpackhi_epi16(x[ 2], x[3]);
    m[ 6] = _mm_unpacklo_epi16(d[ 2], d[3]);
    m[ 7] = _mm_unpackhi_epi16(d[ 2], d[3]);

    a[ 0] = _mm_madd_epi16(m[ 0], m[ 2]);
    a[ 1] = _mm_madd_epi16(m[ 1], m[ 3]);
    a[ 2] = _mm_madd_epi16(m[ 4], m[ 6]);
    a[ 3] = _mm_madd_epi16(m[ 5], m[ 7]);

    b[ 0] = _mm_add_epi32(a[ 0], a[ 2]);
    b[ 1] = _mm_add_epi32(a[ 1], a[ 3]);

    m[ 8] = _mm_unpacklo_epi16(x[ 4], x[5]);
    m[ 9] = _mm_unpackhi_epi16(x[ 4], x[5]);
    m[10] = _mm_unpacklo_epi16(d[ 4], d[5]);
    m[11] = _mm_unpackhi_epi16(d[ 4], d[5]);
    m[12] = _mm_unpacklo_epi16(x[ 6], x[7]);
    m[13] = _mm_unpackhi_epi16(x[ 6], x[7]);
    m[14] = _mm_unpacklo_epi16(d[ 6], d[7]);
    m[15] = _mm_unpackhi_epi16(d[ 6], d[7]);

    a[ 4] = _mm_madd_epi16(m[ 8], m[10]);
    a[ 5] = _mm_madd_epi16(m[ 9], m[11]);
    a[ 6] = _mm_madd_epi16(m[12], m[14]);
    a[ 7] = _mm_madd_epi16(m[13], m[15]);

    b[ 2] = _mm_add_epi32(a[ 4], a[ 6]);
    b[ 3] = _mm_add_epi32(a[ 5], a[ 7]);

    c[ 0] = _mm_add_epi32(b[ 0], b[ 2]);
    c[ 1] = _mm_add_epi32(b[ 1], b[ 3]);

    m[16] = _mm_unpacklo_epi16(x[ 8], x[9]);
    m[17] = _mm_unpackhi_epi16(x[ 8], x[9]);
    m[18] = _mm_unpacklo_epi16(d[ 8], d[9]);
    m[19] = _mm_unpackhi_epi16(d[ 8], d[9]);
    m[20] = _mm_unpacklo_epi16(x[10], x[11]);
    m[21] = _mm_unpackhi_epi16(x[10], x[11]);
    m[22] = _mm_unpacklo_epi16(d[10], d[11]);
    m[23] = _mm_unpackhi_epi16(d[10], d[11]);

    a[ 8] = _mm_madd_epi16(m[16], m[18]);
    a[ 9] = _mm_madd_epi16(m[17], m[19]);
    a[10] = _mm_madd_epi16(m[20], m[22]);
    a[11] = _mm_madd_epi16(m[21], m[23]);

    b[ 4] = _mm_add_epi32(a[ 8], a[10]);
    b[ 5] = _mm_add_epi32(a[ 9], a[11]);

    m[24] = _mm_unpacklo_epi16(x[12], x[13]);
    m[25] = _mm_unpackhi_epi16(x[12], x[13]);
    m[26] = _mm_unpacklo_epi16(d[12], d[13]);
    m[27] = _mm_unpackhi_epi16(d[12], d[13]);
    m[28] = _mm_unpacklo_epi16(x[14], x[15]);
    m[29] = _mm_unpackhi_epi16(x[14], x[15]);
    m[30] = _mm_unpacklo_epi16(d[14], d[15]);
    m[31] = _mm_unpackhi_epi16(d[14], d[15]);

    a[12] = _mm_madd_epi16(m[24], m[26]);
    a[13] = _mm_madd_epi16(m[25], m[27]);
    a[14] = _mm_madd_epi16(m[28], m[30]);
    a[15] = _mm_madd_epi16(m[29], m[31]);

    b[ 6] = _mm_add_epi32(a[12], a[14]);
    b[ 7] = _mm_add_epi32(a[13], a[15]);

    c[ 2] = _mm_add_epi32(b[ 4], b[ 6]);
    c[ 3] = _mm_add_epi32(b[ 5], b[ 7]);
    #endif

    r[0] = _mm_add_epi32(c[0], c[2]);
    r[1] = _mm_add_epi32(c[1], c[3]);
}

static inline void
matMult16x16_red2(__m128i *x, __m128i *d, __m128i *r)
{
    __m128i m[32];//, a[16], b[8], c[4];

    m[ 0] = _mm_unpacklo_epi16(x[ 0], x[1]);
    m[ 2] = _mm_unpacklo_epi16(d[ 0], d[1]);

    r[ 0] = _mm_madd_epi16(m[ 0], m[ 2]);

    m[ 1] = _mm_unpackhi_epi16(x[ 0], x[1]);
    m[ 3] = _mm_unpackhi_epi16(d[ 0], d[1]);

    r[ 1] = _mm_madd_epi16(m[ 1], m[ 3]);
}

static inline void
matMult16x16_red4(__m128i *x, __m128i *d, __m128i *r)
{
    __m128i m[32], a[16];//, b[8], c[4];

    m[ 0] = _mm_unpacklo_epi16(x[ 0], x[1]);
    m[ 2] = _mm_unpacklo_epi16(d[ 0], d[1]);
    m[ 4] = _mm_unpacklo_epi16(x[ 2], x[3]);
    m[ 6] = _mm_unpacklo_epi16(d[ 2], d[3]);

    a[ 0] = _mm_madd_epi16(m[ 0], m[ 2]);
    a[ 2] = _mm_madd_epi16(m[ 4], m[ 6]);

    r[ 0] = _mm_add_epi32(a[ 0], a[ 2]);

    m[ 1] = _mm_unpackhi_epi16(x[ 0], x[1]);
    m[ 3] = _mm_unpackhi_epi16(d[ 0], d[1]);
    m[ 5] = _mm_unpackhi_epi16(x[ 2], x[3]);
    m[ 7] = _mm_unpackhi_epi16(d[ 2], d[3]);

    a[ 1] = _mm_madd_epi16(m[ 1], m[ 3]);
    a[ 3] = _mm_madd_epi16(m[ 5], m[ 7]);

    r[ 1] = _mm_add_epi32(a[ 1], a[ 3]);
}

static inline void
matMult16x16_red8(__m128i *x, __m128i *d, __m128i *r)
{
    __m128i m[32], a[16], b[8];//, c[4];

    m[ 0] = _mm_unpacklo_epi16(x[ 0], x[1]);
    m[ 2] = _mm_unpacklo_epi16(d[ 0], d[1]);
    m[ 4] = _mm_unpacklo_epi16(x[ 2], x[3]);
    m[ 6] = _mm_unpacklo_epi16(d[ 2], d[3]);

    a[ 0] = _mm_madd_epi16(m[ 0], m[ 2]);
    a[ 2] = _mm_madd_epi16(m[ 4], m[ 6]);

    b[ 0] = _mm_add_epi32(a[ 0], a[ 2]);

    m[ 1] = _mm_unpackhi_epi16(x[ 0], x[1]);
    m[ 3] = _mm_unpackhi_epi16(d[ 0], d[1]);
    m[ 5] = _mm_unpackhi_epi16(x[ 2], x[3]);
    m[ 7] = _mm_unpackhi_epi16(d[ 2], d[3]);

    a[ 1] = _mm_madd_epi16(m[ 1], m[ 3]);
    a[ 3] = _mm_madd_epi16(m[ 5], m[ 7]);

    b[ 1] = _mm_add_epi32(a[ 1], a[ 3]);

    m[ 8] = _mm_unpacklo_epi16(x[ 4], x[5]);
    m[ 9] = _mm_unpackhi_epi16(x[ 4], x[5]);
    m[10] = _mm_unpacklo_epi16(d[ 4], d[5]);
    m[11] = _mm_unpackhi_epi16(d[ 4], d[5]);
    m[12] = _mm_unpacklo_epi16(x[ 6], x[7]);
    m[13] = _mm_unpackhi_epi16(x[ 6], x[7]);
    m[14] = _mm_unpacklo_epi16(d[ 6], d[7]);
    m[15] = _mm_unpackhi_epi16(d[ 6], d[7]);

    a[ 4] = _mm_madd_epi16(m[ 8], m[10]);
    a[ 5] = _mm_madd_epi16(m[ 9], m[11]);
    a[ 6] = _mm_madd_epi16(m[12], m[14]);
    a[ 7] = _mm_madd_epi16(m[13], m[15]);

    b[ 2] = _mm_add_epi32(a[ 4], a[ 6]);
    b[ 3] = _mm_add_epi32(a[ 5], a[ 7]);

    r[ 0] = _mm_add_epi32(b[ 0], b[ 2]);
    r[ 1] = _mm_add_epi32(b[ 1], b[ 3]);
}

static inline void
transform_4x2(__m128i *x, __m128i *d, __m128i *r)
{
    #if 0
    __m128i m[4];

    multiply(x+0, d+0, m+0);
    multiply(x+1, d+1, m+2);

    r[0] = _mm_add_epi32(m[0], m[2]);
    r[1] = _mm_add_epi32(m[1], m[3]);
    #else
    __m128i m[4];
    m[0] = _mm_unpacklo_epi16(x[0], x[1]);
    m[1] = _mm_unpackhi_epi16(x[0], x[1]);
    m[2] = _mm_unpacklo_epi16(d[0], d[1]);
    m[3] = _mm_unpackhi_epi16(d[0], d[1]);

    r[0] = _mm_madd_epi16(m[0], m[2]);
    r[1] = _mm_madd_epi16(m[1], m[3]);
    #endif
}

static inline void
transform_4x4(__m128i *x, __m128i *d, __m128i *r)
{
    __m128i m[4];

    transform_4x2(x, d+0, m+0);
    transform_4x2(x, d+2, m+2);

    r[0] = _mm_add_epi32(m[0], m[1]);
    r[1] = _mm_add_epi32(m[2], m[3]);
    r[2] = _mm_sub_epi32(m[0], m[1]);
    r[3] = _mm_sub_epi32(m[2], m[3]);
}

static inline void
transform_4x8(__m128i *x, __m128i *d, __m128i *r)
{
    // __m128i m0[2], m1[2];
    __m128i a[8];

    #if 0
    multiply(x+0, d+0, m0);
    multiply(x+2, d+4, m1);
    a[0] = _mm_add_epi32(m0[0], m1[0]);
    a[1] = _mm_add_epi32(m0[1], m1[1]);

    multiply(x+0, d+1, m0);
    multiply(x+2, d+5, m1);
    a[2] = _mm_add_epi32(m0[0], m1[0]);
    a[3] = _mm_add_epi32(m0[1], m1[1]);

    multiply(x+1, d+2, m0);
    multiply(x+3, d+6, m1);
    a[4] = _mm_add_epi32(m0[0], m1[0]);
    a[5] = _mm_add_epi32(m0[1], m1[1]);

    multiply(x+1, d+3, m0);
    multiply(x+3, d+7, m1);
    a[6] = _mm_add_epi32(m0[0], m1[0]);
    a[7] = _mm_add_epi32(m0[1], m1[1]);
    #else

    __m128i m[4];
    m[0] = _mm_unpacklo_epi16(x[0], x[2]);
    m[1] = _mm_unpackhi_epi16(x[0], x[2]);
    m[2] = _mm_unpacklo_epi16(d[0], d[4]);
    m[3] = _mm_unpackhi_epi16(d[0], d[4]);

    a[0] = _mm_madd_epi16(m[0], m[2]);
    a[1] = _mm_madd_epi16(m[1], m[3]);

    m[2] = _mm_unpacklo_epi16(d[1], d[5]);
    m[3] = _mm_unpackhi_epi16(d[1], d[5]);

    a[2] = _mm_madd_epi16(m[0], m[2]);
    a[3] = _mm_madd_epi16(m[1], m[3]);

    m[0] = _mm_unpacklo_epi16(x[1], x[3]);
    m[1] = _mm_unpackhi_epi16(x[1], x[3]);
    m[2] = _mm_unpacklo_epi16(d[2], d[6]);
    m[3] = _mm_unpackhi_epi16(d[2], d[6]);

    a[4] = _mm_madd_epi16(m[0], m[2]);
    a[5] = _mm_madd_epi16(m[1], m[3]);

    m[2] = _mm_unpacklo_epi16(d[3], d[7]);
    m[3] = _mm_unpackhi_epi16(d[3], d[7]);

    a[6] = _mm_madd_epi16(m[0], m[2]);
    a[7] = _mm_madd_epi16(m[1], m[3]);

    #endif

    r[0] = _mm_add_epi32(a[0], a[4]);
    r[1] = _mm_add_epi32(a[1], a[5]);
    r[2] = _mm_add_epi32(a[2], a[6]);
    r[3] = _mm_add_epi32(a[3], a[7]);
    r[4] = _mm_sub_epi32(a[2], a[6]);
    r[5] = _mm_sub_epi32(a[3], a[7]);
    r[6] = _mm_sub_epi32(a[0], a[4]);
    r[7] = _mm_sub_epi32(a[1], a[5]);
}

static inline void
dct2_4x2(__m128i *x, __m128i *d, __m128i *r)
{
    __m128i m[2], a[2];

    transform_4x2(x,d,m);

    a[0] = _mm_add_epi32(m[0], m[1]);
    a[1] = _mm_sub_epi32(m[0], m[1]);

    r[0] = _mm_unpacklo_epi64(a[0], a[1]);
    r[1] = _mm_unpackhi_epi64(a[0], a[1]);
}

static inline void
dct2_4x4_red1(__m128i *x, __m128i *d, __m128i *r)
{
    __m128i m[4], t[4], z;

    z = _mm_setzero_si128();

    /*FIXME we could use swizzling to make more
     efficient use of madd here */
    t[0] = _mm_unpacklo_epi16(x[0], z);
    t[1] = _mm_unpackhi_epi16(x[0], z);
    t[2] = _mm_unpacklo_epi16(d[0], z);
    t[3] = _mm_unpackhi_epi16(d[0], z);

    m[0] = _mm_madd_epi16(t[0], t[2]);
    m[1] = _mm_madd_epi16(t[1], t[3]);

    t[2] = _mm_unpacklo_epi16(d[2], z);
    t[3] = _mm_unpackhi_epi16(d[2], z);

    m[2] = _mm_madd_epi16(t[0], t[2]);
    m[3] = _mm_madd_epi16(t[1], t[3]);

    t[0] = _mm_add_epi32(m[0], m[1]);
    t[1] = _mm_add_epi32(m[2], m[3]);
    t[2] = _mm_sub_epi32(m[0], m[1]);
    t[3] = _mm_sub_epi32(m[2], m[3]);

    /*FIXME transpose and shuffle are
     redundant blend might a better choice*/
    transpose4x4(t, r);

    r[0] = _mm_shuffle_epi32(r[0], 0xB4);
    r[1] = _mm_shuffle_epi32(r[1], 0xE1);
    r[2] = _mm_shuffle_epi32(r[2], 0xB4);
    r[3] = _mm_shuffle_epi32(r[3], 0xE1);
}

static inline void
dct2_4x4_red2(__m128i *x, __m128i *d, __m128i *r)
{
    __m128i m[4], t[4];

    /*FIXME we could use swizzling to make more
     efficient use of madd here */
    t[0] = _mm_unpacklo_epi16(x[0], x[1]);
    t[1] = _mm_unpackhi_epi16(x[0], x[1]);
    t[2] = _mm_unpacklo_epi16(d[0], d[1]);
    t[3] = _mm_unpackhi_epi16(d[0], d[1]);

    m[0] = _mm_madd_epi16(t[0], t[2]);
    m[1] = _mm_madd_epi16(t[1], t[3]);

    t[2] = _mm_unpacklo_epi16(d[2], d[3]);
    t[3] = _mm_unpackhi_epi16(d[2], d[3]);

    m[2] = _mm_madd_epi16(t[0], t[2]);
    m[3] = _mm_madd_epi16(t[1], t[3]);

    t[0] = _mm_add_epi32(m[0], m[1]);
    t[1] = _mm_add_epi32(m[2], m[3]);
    t[2] = _mm_sub_epi32(m[0], m[1]);
    t[3] = _mm_sub_epi32(m[2], m[3]);

    transpose4x4(t, r);

    r[0] = _mm_shuffle_epi32(r[0], 0xB4);
    r[1] = _mm_shuffle_epi32(r[1], 0xE1);
    r[2] = _mm_shuffle_epi32(r[2], 0xB4);
    r[3] = _mm_shuffle_epi32(r[3], 0xE1);
}

static inline void
dct2_4x4(__m128i *x, __m128i *d, __m128i *r)
{
    __m128i m[4];

    transform_4x4(x, d, m);

    transpose4x4(m, r);

    r[0] = _mm_shuffle_epi32(r[0], 0xB4);
    r[1] = _mm_shuffle_epi32(r[1], 0xE1);
    r[2] = _mm_shuffle_epi32(r[2], 0xB4);
    r[3] = _mm_shuffle_epi32(r[3], 0xE1);
}

static inline void
dct2_4x8(__m128i *x, __m128i *d, __m128i *r)
{
    __m128i m[8];

    transform_4x8(x, d, m);

    transpose4x4_step2(m  , r  );
    transpose4x4_step2(m+1, r+4);
}

static inline void
dct2_8x2(__m128i *x, __m128i *d, __m128i *r)
{
    __m128i e[2], o[2];

    dct2_4x2(x+0, d+0, e);

    matMult4x4(x+2, d+2, o);

    r[0] = _mm_add_epi32(e[0], o[0]);
    r[1] = _mm_sub_epi32(e[0], o[0]);
    r[2] = _mm_add_epi32(e[1], o[1]);
    r[3] = _mm_sub_epi32(e[1], o[1]);

    r[0] = _mm_shuffle_epi32(r[0], 0xB4);
    r[1] = _mm_shuffle_epi32(r[1], 0x1E);
    r[2] = _mm_shuffle_epi32(r[2], 0xB4);
    r[3] = _mm_shuffle_epi32(r[3], 0x1E);
}

static inline void
dct2_8x4_red1(__m128i *x, __m128i *d, __m128i *r)
{
    __m128i e[4];//, o[4];

    dct2_4x4_red1(x+0, d+0, e);

    r[0] = e[0];
    r[2] = e[1];
    r[4] = e[2];
    r[6] = e[3];

    r[1] = _mm_shuffle_epi32(e[0], 0x1B);
    r[3] = _mm_shuffle_epi32(e[1], 0x1B);
    r[5] = _mm_shuffle_epi32(e[2], 0x1B);
    r[7] = _mm_shuffle_epi32(e[3], 0x1B);
}

static inline void
dct2_8x4_red2(__m128i *x, __m128i *d, __m128i *r)
{
    __m128i e[4], o[4];

    dct2_4x4_red1(x+0, d+0, e);

    matMult4x4_red1(x+2, d+4, o+0);
    matMult4x4_red1(x+6, d+4, o+2);

    o[0] = _mm_shuffle_epi32(o[0], 0xB4);
    o[1] = _mm_shuffle_epi32(o[1], 0xB4);
    o[2] = _mm_shuffle_epi32(o[2], 0xB4);
    o[3] = _mm_shuffle_epi32(o[3], 0xB4);

    r[0] = _mm_add_epi32(e[0], o[0]);
    r[1] = _mm_sub_epi32(e[0], o[0]);
    r[2] = _mm_add_epi32(e[1], o[1]);
    r[3] = _mm_sub_epi32(e[1], o[1]);
    r[4] = _mm_add_epi32(e[2], o[2]);
    r[5] = _mm_sub_epi32(e[2], o[2]);
    r[6] = _mm_add_epi32(e[3], o[3]);
    r[7] = _mm_sub_epi32(e[3], o[3]);

    r[1] = _mm_shuffle_epi32(r[1], 0x1B);
    r[3] = _mm_shuffle_epi32(r[3], 0x1B);
    r[5] = _mm_shuffle_epi32(r[5], 0x1B);
    r[7] = _mm_shuffle_epi32(r[7], 0x1B);
}

static inline void
dct2_8x4_red4(__m128i *x, __m128i *d, __m128i *r)
{
    __m128i e[4], o[4];

    dct2_4x4_red2(x+0, d+0, e);

    matMult4x4_red2(x+2, d+4, o+0);
    matMult4x4_red2(x+6, d+4, o+2);

    o[0] = _mm_shuffle_epi32(o[0], 0xB4);
    o[1] = _mm_shuffle_epi32(o[1], 0xB4);
    o[2] = _mm_shuffle_epi32(o[2], 0xB4);
    o[3] = _mm_shuffle_epi32(o[3], 0xB4);

    r[0] = _mm_add_epi32(e[0], o[0]);
    r[1] = _mm_sub_epi32(e[0], o[0]);
    r[2] = _mm_add_epi32(e[1], o[1]);
    r[3] = _mm_sub_epi32(e[1], o[1]);
    r[4] = _mm_add_epi32(e[2], o[2]);
    r[5] = _mm_sub_epi32(e[2], o[2]);
    r[6] = _mm_add_epi32(e[3], o[3]);
    r[7] = _mm_sub_epi32(e[3], o[3]);

    r[1] = _mm_shuffle_epi32(r[1], 0x1B);
    r[3] = _mm_shuffle_epi32(r[3], 0x1B);
    r[5] = _mm_shuffle_epi32(r[5], 0x1B);
    r[7] = _mm_shuffle_epi32(r[7], 0x1B);
}

static inline void
dct2_8x4(__m128i *x, __m128i *d, __m128i *r)
{
    __m128i e[4], o[4];

    dct2_4x4(x+0, d+0, e);

    matMult4x4(x+2, d+4, o+0);
    matMult4x4(x+6, d+4, o+2);

    o[0] = _mm_shuffle_epi32(o[0], 0xB4);
    o[1] = _mm_shuffle_epi32(o[1], 0xB4);
    o[2] = _mm_shuffle_epi32(o[2], 0xB4);
    o[3] = _mm_shuffle_epi32(o[3], 0xB4);

    r[0] = _mm_add_epi32(e[0], o[0]);
    r[1] = _mm_sub_epi32(e[0], o[0]);
    r[2] = _mm_add_epi32(e[1], o[1]);
    r[3] = _mm_sub_epi32(e[1], o[1]);
    r[4] = _mm_add_epi32(e[2], o[2]);
    r[5] = _mm_sub_epi32(e[2], o[2]);
    r[6] = _mm_add_epi32(e[3], o[3]);
    r[7] = _mm_sub_epi32(e[3], o[3]);

    r[1] = _mm_shuffle_epi32(r[1], 0x1B);
    r[3] = _mm_shuffle_epi32(r[3], 0x1B);
    r[5] = _mm_shuffle_epi32(r[5], 0x1B);
    r[7] = _mm_shuffle_epi32(r[7], 0x1B);
}

static inline void
dct2_8x8(__m128i *x, __m128i *d, __m128i *r)
{
    __m128i e[8], o[8], oo[8];

    dct2_4x8(x+0, d+0, e);

    matMult4x4(x+4, d+ 8, oo+0);
    matMult4x4(x+4, d+12, oo+2);
    matMult4x4(x+4, d+16, oo+4);
    matMult4x4(x+4, d+20, oo+6);

    transpose4x4_step2(oo+0, o+0);
    transpose4x4_step2(oo+1, o+4);

    for(int k=0; k<8; k++){
        r[  k] = _mm_add_epi32(e[k], o[k]);
        r[8+k] = _mm_shuffle_epi32(_mm_sub_epi32(e[k],  o[k]),  0x1B);
    }
}

static inline void
dct2_16x2(__m128i *x, __m128i *d, __m128i *r)
{
    __m128i e[4], o[4];

    dct2_8x2(x+0, d+0, e);

    matMult8x8(x+ 6, d+6, o+0);
    matMult8x8(x+14, d+6, o+2);

    r[0] = _mm_add_epi32(e[0], o[0]);
    r[1] = _mm_add_epi32(e[1], o[1]);
    r[2] = _mm_sub_epi32(e[1], o[1]);
    r[3] = _mm_sub_epi32(e[0], o[0]);
    r[4] = _mm_add_epi32(e[2], o[2]);
    r[5] = _mm_add_epi32(e[3], o[3]);
    r[6] = _mm_sub_epi32(e[3], o[3]);
    r[7] = _mm_sub_epi32(e[2], o[2]);

    r[2] = _mm_shuffle_epi32(r[2], 0x1B);
    r[3] = _mm_shuffle_epi32(r[3], 0x1B);
    r[6] = _mm_shuffle_epi32(r[6], 0x1B);
    r[7] = _mm_shuffle_epi32(r[7], 0x1B);
}

static inline void
dct2_16x4_red2(__m128i *x, __m128i *d, __m128i *r)
{
    __m128i e[8], o[8];
    __m128i *odd_x = x + 10;

    dct2_8x4_red1(x+0, d+0, e);

    matMult8x8_red1(odd_x +  0, d+8, o+0);
    matMult8x8_red1(odd_x +  8, d+8, o+2);
    matMult8x8_red1(odd_x + 16, d+8, o+4);
    matMult8x8_red1(odd_x + 24, d+8, o+6);

    for(int k=0; k<8; k+=2){
        r[2*k+0] = _mm_add_epi32(e[k+0], o[k+0]);
        r[2*k+1] = _mm_add_epi32(e[k+1], o[k+1]);
        r[2*k+2] = _mm_shuffle_epi32(_mm_sub_epi32(e[k+1], o[k+1]), 0x1B);
        r[2*k+3] = _mm_shuffle_epi32(_mm_sub_epi32(e[k+0], o[k+0]), 0x1B);
    }
}

static inline void
dct2_16x4_red4(__m128i *x, __m128i *d, __m128i *r)
{
    __m128i e[8], o[8];
    __m128i *odd_x = x + 10;

    dct2_8x4_red2(x+0, d+0, e);

    matMult8x8_red2(odd_x +  0, d+8, o+0);
    matMult8x8_red2(odd_x +  8, d+8, o+2);
    matMult8x8_red2(odd_x + 16, d+8, o+4);
    matMult8x8_red2(odd_x + 24, d+8, o+6);

    for(int k=0; k<8; k+=2){
        r[2*k+0] = _mm_add_epi32(e[k+0], o[k+0]);
        r[2*k+1] = _mm_add_epi32(e[k+1], o[k+1]);
        r[2*k+2] = _mm_shuffle_epi32(_mm_sub_epi32(e[k+1], o[k+1]), 0x1B);
        r[2*k+3] = _mm_shuffle_epi32(_mm_sub_epi32(e[k+0], o[k+0]), 0x1B);
    }
}

static inline void
dct2_16x4_red8(__m128i *x, __m128i *d, __m128i *r)
{
    __m128i e[8], o[8];
    __m128i *odd_x = x + 10;

    dct2_8x4_red4(x+0, d+0, e);

    matMult8x8_red4(odd_x +  0, d+8, o+0);
    matMult8x8_red4(odd_x +  8, d+8, o+2);
    matMult8x8_red4(odd_x + 16, d+8, o+4);
    matMult8x8_red4(odd_x + 24, d+8, o+6);

    for(int k=0; k<8; k+=2){
        r[2*k+0] = _mm_add_epi32(e[k+0], o[k+0]);
        r[2*k+1] = _mm_add_epi32(e[k+1], o[k+1]);
        r[2*k+2] = _mm_shuffle_epi32(_mm_sub_epi32(e[k+1], o[k+1]), 0x1B);
        r[2*k+3] = _mm_shuffle_epi32(_mm_sub_epi32(e[k+0], o[k+0]), 0x1B);
    }
}

static inline void
dct2_16x4(__m128i *x, __m128i *d, __m128i *r)
{
    __m128i e[8], o[8];

    dct2_8x4(x+0, d+0, e);

    matMult8x8(x+10, d+8, o+0);
    matMult8x8(x+18, d+8, o+2);
    matMult8x8(x+26, d+8, o+4);
    matMult8x8(x+34, d+8, o+6);

    for(int k=0; k<8; k+=2){
        r[2*k+0] = _mm_add_epi32(e[k+0], o[k+0]);
        r[2*k+1] = _mm_add_epi32(e[k+1], o[k+1]);
        r[2*k+2] = _mm_shuffle_epi32(_mm_sub_epi32(e[k+1], o[k+1]), 0x1B);
        r[2*k+3] = _mm_shuffle_epi32(_mm_sub_epi32(e[k+0], o[k+0]), 0x1B);
    }
}

static inline void
dct2_16x8(__m128i *x, __m128i *d, __m128i *r)
{
    __m128i e[16], o[16];

    dct2_8x8(x+0, d+0, e);

    for(int k=0; k<8; k++){
        matMult8x8(x+8+8*k, d+24, o+2*k);
    }

    for(int k=0; k<8; k++){
        r[4*k+0] = _mm_add_epi32(e[k+0], o[2*k+0]);
        r[4*k+1] = _mm_add_epi32(e[k+8], o[2*k+1]);
        r[4*k+2] = _mm_shuffle_epi32(_mm_sub_epi32(e[k+8], o[2*k+1]), 0x1B);
        r[4*k+3] = _mm_shuffle_epi32(_mm_sub_epi32(e[k+0], o[2*k+0]), 0x1B);
    }
}

static inline void
dct2_32x2(__m128i *x, __m128i *d, __m128i *r)
{
    __m128i e[8], o[8];

    dct2_16x2(x+0, d+0, e);

    for(int k=0; k<2; k++){
        matMult16x16(x+24+16*k, d+14, o+4*k+0);
        matMult16x16(x+24+16*k, d+30, o+4*k+2);
    }

    for(int k=0; k<2; k++){
        r[8*k+0] = _mm_add_epi32(e[4*k+0], o[4*k+0]);
        r[8*k+1] = _mm_add_epi32(e[4*k+1], o[4*k+1]);
        r[8*k+2] = _mm_add_epi32(e[4*k+2], o[4*k+2]);
        r[8*k+3] = _mm_add_epi32(e[4*k+3], o[4*k+3]);
        r[8*k+4] = _mm_shuffle_epi32(_mm_sub_epi32(e[4*k+3], o[4*k+3]), 0x1B);
        r[8*k+5] = _mm_shuffle_epi32(_mm_sub_epi32(e[4*k+2], o[4*k+2]), 0x1B);
        r[8*k+6] = _mm_shuffle_epi32(_mm_sub_epi32(e[4*k+1], o[4*k+1]), 0x1B);
        r[8*k+7] = _mm_shuffle_epi32(_mm_sub_epi32(e[4*k+0], o[4*k+0]), 0x1B);
    }
}

static inline void
dct2_32x4_red4(__m128i *x, __m128i *d, __m128i *r)
{
    __m128i e[16], o[16];

    dct2_16x4_red2(x+0, d+0, e);

    for(int k=0; k<4; k++){
        matMult16x16_red2(x+74+16*k, d+16, o+4*k+0);
        matMult16x16_red2(x+74+16*k, d+32, o+4*k+2);
    }

    for(int k=0; k<4; k++){
        r[8*k+0] = _mm_add_epi32(e[4*k+0], o[4*k+0]);
        r[8*k+1] = _mm_add_epi32(e[4*k+1], o[4*k+1]);
        r[8*k+2] = _mm_add_epi32(e[4*k+2], o[4*k+2]);
        r[8*k+3] = _mm_add_epi32(e[4*k+3], o[4*k+3]);
        r[8*k+4] = _mm_shuffle_epi32(_mm_sub_epi32(e[4*k+3], o[4*k+3]), 0x1B);
        r[8*k+5] = _mm_shuffle_epi32(_mm_sub_epi32(e[4*k+2], o[4*k+2]), 0x1B);
        r[8*k+6] = _mm_shuffle_epi32(_mm_sub_epi32(e[4*k+1], o[4*k+1]), 0x1B);
        r[8*k+7] = _mm_shuffle_epi32(_mm_sub_epi32(e[4*k+0], o[4*k+0]), 0x1B);
    }
}

static inline void
dct2_32x4_red8(__m128i *x, __m128i *d, __m128i *r)
{
    __m128i e[16], o[16];

    dct2_16x4_red4(x+0, d+0, e);

    for(int k=0; k<4; k++){
        matMult16x16_red4(x+74+16*k, d+16, o+4*k+0);
        matMult16x16_red4(x+74+16*k, d+32, o+4*k+2);
    }

    for(int k=0; k<4; k++){
        r[8*k+0] = _mm_add_epi32(e[4*k+0], o[4*k+0]);
        r[8*k+1] = _mm_add_epi32(e[4*k+1], o[4*k+1]);
        r[8*k+2] = _mm_add_epi32(e[4*k+2], o[4*k+2]);
        r[8*k+3] = _mm_add_epi32(e[4*k+3], o[4*k+3]);
        r[8*k+4] = _mm_shuffle_epi32(_mm_sub_epi32(e[4*k+3], o[4*k+3]), 0x1B);
        r[8*k+5] = _mm_shuffle_epi32(_mm_sub_epi32(e[4*k+2], o[4*k+2]), 0x1B);
        r[8*k+6] = _mm_shuffle_epi32(_mm_sub_epi32(e[4*k+1], o[4*k+1]), 0x1B);
        r[8*k+7] = _mm_shuffle_epi32(_mm_sub_epi32(e[4*k+0], o[4*k+0]), 0x1B);
    }
}

static inline void
dct2_32x4_red16(__m128i *x, __m128i *d, __m128i *r)
{
    __m128i e[16], o[16];

    dct2_16x4_red8(x+0, d+0, e);

    for(int k=0; k<4; k++){
        matMult16x16_red8(x+74+16*k, d+16, o+4*k+0);
        matMult16x16_red8(x+74+16*k, d+32, o+4*k+2);
    }

    for(int k=0; k<4; k++){
        r[8*k+0] = _mm_add_epi32(e[4*k+0], o[4*k+0]);
        r[8*k+1] = _mm_add_epi32(e[4*k+1], o[4*k+1]);
        r[8*k+2] = _mm_add_epi32(e[4*k+2], o[4*k+2]);
        r[8*k+3] = _mm_add_epi32(e[4*k+3], o[4*k+3]);
        r[8*k+4] = _mm_shuffle_epi32(_mm_sub_epi32(e[4*k+3], o[4*k+3]), 0x1B);
        r[8*k+5] = _mm_shuffle_epi32(_mm_sub_epi32(e[4*k+2], o[4*k+2]), 0x1B);
        r[8*k+6] = _mm_shuffle_epi32(_mm_sub_epi32(e[4*k+1], o[4*k+1]), 0x1B);
        r[8*k+7] = _mm_shuffle_epi32(_mm_sub_epi32(e[4*k+0], o[4*k+0]), 0x1B);
    }
}

static inline void
dct2_32x4(__m128i *x, __m128i *d, __m128i *r)
{
    __m128i e[16], o[16];

    dct2_16x4(x+0, d+0, e);

    for(int k=0; k<4; k++){
        matMult16x16(x+74+16*k, d+16, o+4*k+0);
        matMult16x16(x+74+16*k, d+32, o+4*k+2);
    }

    for(int k=0; k<4; k++){
        r[8*k+0] = _mm_add_epi32(e[4*k+0], o[4*k+0]);
        r[8*k+1] = _mm_add_epi32(e[4*k+1], o[4*k+1]);
        r[8*k+2] = _mm_add_epi32(e[4*k+2], o[4*k+2]);
        r[8*k+3] = _mm_add_epi32(e[4*k+3], o[4*k+3]);
        r[8*k+4] = _mm_shuffle_epi32(_mm_sub_epi32(e[4*k+3], o[4*k+3]), 0x1B);
        r[8*k+5] = _mm_shuffle_epi32(_mm_sub_epi32(e[4*k+2], o[4*k+2]), 0x1B);
        r[8*k+6] = _mm_shuffle_epi32(_mm_sub_epi32(e[4*k+1], o[4*k+1]), 0x1B);
        r[8*k+7] = _mm_shuffle_epi32(_mm_sub_epi32(e[4*k+0], o[4*k+0]), 0x1B);
    }
}

static inline void
dct2_32x8(__m128i *x, __m128i *d, __m128i *r){
    __m128i e[32], o[32];

    dct2_16x8(x+0, d+0, e+0);

    for(int k=0; k<8; k++){
        matMult16x16(x+72+16*k, d+32, o+4*k+0);
        matMult16x16(x+72+16*k, d+48, o+4*k+2);
    }

    for(int k=0; k<8; k++){
        r[8*k+0] = _mm_add_epi32(e[4*k+0], o[4*k+0]);
        r[8*k+1] = _mm_add_epi32(e[4*k+1], o[4*k+1]);
        r[8*k+2] = _mm_add_epi32(e[4*k+2], o[4*k+2]);
        r[8*k+3] = _mm_add_epi32(e[4*k+3], o[4*k+3]);
        r[8*k+4] = _mm_shuffle_epi32(_mm_sub_epi32(e[4*k+3], o[4*k+3]), 0x1B);
        r[8*k+5] = _mm_shuffle_epi32(_mm_sub_epi32(e[4*k+2], o[4*k+2]), 0x1B);
        r[8*k+6] = _mm_shuffle_epi32(_mm_sub_epi32(e[4*k+1], o[4*k+1]), 0x1B);
        r[8*k+7] = _mm_shuffle_epi32(_mm_sub_epi32(e[4*k+0], o[4*k+0]), 0x1B);
    }
}

void vvc_inverse_dct_ii_2_sse(const int16_t *src, int16_t *dst, ptrdiff_t src_stride,
                              int num_lines, int line_brk, int shift)
{
    __m128i x1, x2;
    __m128i e, o;
    __m128i outl, outh;

    __m128i add = _mm_set1_epi32(1<<(shift-6-1));
    __m128i val1 = _mm_set1_epi16(1);
    for(int i = 0; i < num_lines / 8; i++){
        //Même chose d'utiliser unpacklo et unpackhi ici
        x1 = _mm_load_si128((__m128i *)&src[0]);
        x2 = _mm_load_si128((__m128i *)&src[src_stride]);

        __m128i xylo = _mm_mullo_epi16(x1, val1);
        __m128i xyhi = _mm_mulhi_epi16(x1, val1);
        __m128i x1lo = _mm_unpacklo_epi16(xylo, xyhi);
        __m128i x1hi = _mm_unpackhi_epi16(xylo, xyhi);

        xylo = _mm_mullo_epi16(x2, val1);
        xyhi = _mm_mulhi_epi16(x2, val1);
        __m128i x2lo = _mm_unpacklo_epi16(xylo, xyhi);
        __m128i x2hi = _mm_unpackhi_epi16(xylo, xyhi);

        __m128i elo = _mm_add_epi32(x1lo,x2lo);
        __m128i ehi = _mm_add_epi32(x1hi,x2hi);
        __m128i olo = _mm_sub_epi32(x1lo,x2lo);
        __m128i ohi = _mm_sub_epi32(x1hi,x2hi);

        elo = _mm_add_epi32(elo, add);
        ehi = _mm_add_epi32(ehi, add);
        olo = _mm_add_epi32(olo, add);
        ohi = _mm_add_epi32(ohi, add);

        elo = _mm_srai_epi32(elo, shift-6);
        ehi = _mm_srai_epi32(ehi, shift-6);
        olo = _mm_srai_epi32(olo, shift-6);
        ohi = _mm_srai_epi32(ohi, shift-6);

        __m128i outllo = _mm_unpacklo_epi32(elo,olo);
        __m128i outhlo = _mm_unpackhi_epi32(elo,olo);
        __m128i outlhi = _mm_unpacklo_epi32(ehi,ohi);
        __m128i outhhi = _mm_unpackhi_epi32(ehi,ohi);

        outl = _mm_packs_epi32(outllo, outhlo);
        outh = _mm_packs_epi32(outlhi, outhhi);

        _mm_store_si128((__m128i *) &(dst[0]), outl);
        _mm_store_si128((__m128i *) &(dst[8]), outh);

        src += 8;
        dst += 16;
    }

    if (!(num_lines & 0x7)) return;

    if (num_lines & 0x4){
        __m128i add = _mm_set1_epi32(1<<(shift-6-1));
        __m128i val1 = _mm_set1_epi16(1);

        x1 = _mm_loadl_epi64((__m128i *)&src[0]);
        x2 = _mm_loadl_epi64((__m128i *)&src[src_stride]);

        __m128i xylo = _mm_mullo_epi16(x1, val1);
        __m128i xyhi = _mm_mulhi_epi16(x1, val1);
        x1 = _mm_unpacklo_epi16(xylo, xyhi);

        xylo = _mm_mullo_epi16(x2, val1);
        xyhi = _mm_mulhi_epi16(x2, val1);
        x2 = _mm_unpacklo_epi16(xylo, xyhi);

        e = _mm_add_epi32(x1,x2);
        o = _mm_sub_epi32(x1,x2);

        e = _mm_add_epi32(e, add);
        o = _mm_add_epi32(o, add);

        e = _mm_srai_epi32(e, shift-6);
        o = _mm_srai_epi32(o, shift-6);

        outl = _mm_packs_epi32(_mm_unpacklo_epi32(e,o),_mm_unpackhi_epi32(e,o));
        _mm_store_si128((__m128i *) &(dst[0]), outl);
    }

    if (num_lines & 0x2){
        x1 = _mm_unpacklo_epi64(_mm_set1_epi32(src[0]),_mm_set1_epi32(src[1]));
        x2 = _mm_unpacklo_epi64(_mm_set1_epi32(src[src_stride]),_mm_set1_epi32(src[src_stride+1]));

        __m128i xo = _mm_set1_epi64x(0xFFFFFFFF00000000);
        x2 = _mm_add_epi32(_mm_xor_si128(x2, xo),_mm_set1_epi64x(0x0000000100000000));
        e = _mm_add_epi32(x1,x2);
        //        eo = _mm_slli_epi32(eo, 6);

        __m128i add = _mm_set1_epi32(1<<(shift-6-1));
        e = _mm_add_epi32(e, add);
        e = _mm_srai_epi32(e, shift-6);
        outl = _mm_packs_epi32(e, e); //clip pour repasser en 16
        _mm_storel_epi64((__m128i *) dst, outl);
    }

    if (num_lines & 0x1){
        vvc_inverse_dct_ii_2(src, dst, src_stride, num_lines & 0x1, line_brk, shift);
    }

}

void
vvc_inverse_dct_ii_4_sse(const int16_t *src, int16_t *dst, ptrdiff_t src_stride,
                        int num_lines, int line_brk, int shift)
{
    __m128i add = _mm_set1_epi32(1 << (shift - 1));

    for (int j = 0; j < num_lines / 8; j++) {
        __m128i x[4], d[8], r[8];

        x[0] = _mm_load_si128((__m128i*)(src + 0 * src_stride));
        x[1] = _mm_load_si128((__m128i*)(src + 1 * src_stride));
        x[2] = _mm_load_si128((__m128i*)(src + 2 * src_stride));
        x[3] = _mm_load_si128((__m128i*)(src + 3 * src_stride));

        d[0] = _mm_set1_epi16(DCT_II_4[0]);
        d[1] = _mm_set1_epi16(DCT_II_4[1]);
        d[2] = _mm_set1_epi16(DCT_II_4[4]);
        d[3] = _mm_set1_epi16(DCT_II_4[5]);
        d[4] = _mm_set1_epi16(DCT_II_4[8]);
        d[5] = _mm_set1_epi16(DCT_II_4[9]);
        d[6] = _mm_set1_epi16(DCT_II_4[12]);
        d[7] = _mm_set1_epi16(DCT_II_4[13]);

        dct2_4x8(x,d,r);

        for(int i=0; i<8; i+=2){
            __m128i o;
            r[i+0] = _mm_add_epi32(r[i+0], add);
            r[i+1] = _mm_add_epi32(r[i+1], add);

            r[i+0] = _mm_srai_epi32(r[i+0], shift);
            r[i+1] = _mm_srai_epi32(r[i+1], shift);

            o = _mm_packs_epi32(r[i+0], r[i+1]);

            _mm_store_si128((__m128i *) (dst+i/2*8), o);
        }
        src += 8;
        dst += 32;
    }

    if (!(num_lines & 0x7)) return;

    if (num_lines & 0x4){
        __m128i x[2], d[4], r[4];

        x[0] = _mm_unpacklo_epi64(
                _mm_loadl_epi64((__m128i*)(src + 0 * src_stride)),
                _mm_loadl_epi64((__m128i*)(src + 1 * src_stride))
                );
        x[1] = _mm_unpacklo_epi64(
                _mm_loadl_epi64((__m128i*)(src + 2 * src_stride)),
                _mm_loadl_epi64((__m128i*)(src + 3 * src_stride))
                );

        static const int16_t DCT_II_4_4_sse[8 * 4] = {
            64,  64,  64,  64,  83,  36,  83,  36,
            64, -64,  64, -64,  36, -83,  36, -83,
            64,  64,  64,  64,  36,  83,  36,  83,
            -64,  64, -64,  64, -83,  36, -83,  36
        };

        d[0] = _mm_load_si128((__m128i*)(DCT_II_4_4_sse + 0));
        d[1] = _mm_load_si128((__m128i*)(DCT_II_4_4_sse + 8));
        d[2] = _mm_load_si128((__m128i*)(DCT_II_4_4_sse + 16));
        d[3] = _mm_load_si128((__m128i*)(DCT_II_4_4_sse + 24));

        dct2_4x4(x,d,r);

        __m128i add = _mm_set1_epi32(1 << (shift - 1));
        r[0] = _mm_add_epi32(r[0], add);
        r[1] = _mm_add_epi32(r[1], add);
        r[2] = _mm_add_epi32(r[2], add);
        r[3] = _mm_add_epi32(r[3], add);

        r[0] = _mm_srai_epi32(r[0], shift);
        r[1] = _mm_srai_epi32(r[1], shift);
        r[2] = _mm_srai_epi32(r[2], shift);
        r[3] = _mm_srai_epi32(r[3], shift);

        __m128i o0 = _mm_packs_epi32(r[0], r[1]);
        __m128i o1 = _mm_packs_epi32(r[2], r[3]);

        _mm_store_si128((__m128i *) (dst+0), o0);
        _mm_store_si128((__m128i *) (dst+8), o1);
    }

    if (num_lines & 0x2){
        __m128i x[2], d[2], r[2];

        x[0] = _mm_unpacklo_epi32(
                _mm_loadl_epi64((__m128i*)(src + 0 * src_stride)),
                _mm_loadl_epi64((__m128i*)(src + 1 * src_stride))
                );
        x[0] = _mm_unpacklo_epi16(x[0],x[0]);
        x[1] = _mm_unpacklo_epi32(
                _mm_loadl_epi64((__m128i*)(src + 2 * src_stride)),
                _mm_loadl_epi64((__m128i*)(src + 3 * src_stride))
                );
        x[1] = _mm_unpacklo_epi16(x[1],x[1]);

        static const int16_t DCT_II_4_2_sse[8 * 2]= {
            64, 64, 64, 64, 83, 36, 83, 36,
            64,-64, 64,-64, 36,-83, 36,-83
        };

        d[0] = _mm_load_si128((__m128i*)(DCT_II_4_2_sse+0));
        d[1] = _mm_load_si128((__m128i*)(DCT_II_4_2_sse+8));

        dct2_4x2(x,d,r);

        r[0] = _mm_shuffle_epi32(r[0], 0xB4);
        r[1] = _mm_shuffle_epi32(r[1], 0xB4);

        __m128i add = _mm_set1_epi32(1 << (shift - 1));
        r[0] = _mm_add_epi32(r[0], add);
        r[1] = _mm_add_epi32(r[1], add);

        r[0] = _mm_srai_epi32(r[0], shift);
        r[1] = _mm_srai_epi32(r[1], shift);

        __m128i out = _mm_packs_epi32(r[0], r[1]);

        _mm_store_si128((__m128i *) dst, out);
    }

    if (num_lines & 0x1){
        vvc_inverse_dct_ii_4(src, dst, src_stride, num_lines & 0x1, line_brk, shift);
    }
}

void
vvc_inverse_dct_ii_8_sse(const int16_t *src, int16_t *dst, ptrdiff_t src_stride,
                         int num_lines, int line_brk, int shift)
{
    __m128i add = _mm_set1_epi32(1 << (shift - 1));

    for (int j = 0; j < num_lines / 8; j++) {
        __m128i x[8], d[24], r[16];

        x[0] = _mm_load_si128((__m128i*)(src + 0 * src_stride));
        x[1] = _mm_load_si128((__m128i*)(src + 2 * src_stride));
        x[2] = _mm_load_si128((__m128i*)(src + 4 * src_stride));
        x[3] = _mm_load_si128((__m128i*)(src + 6 * src_stride));
        x[4] = _mm_load_si128((__m128i*)(src + 1 * src_stride));
        x[5] = _mm_load_si128((__m128i*)(src + 3 * src_stride));
        x[6] = _mm_load_si128((__m128i*)(src + 5 * src_stride));
        x[7] = _mm_load_si128((__m128i*)(src + 7 * src_stride));

        d[0] = _mm_set1_epi16(DCT_II_8[0]);
        d[1] = _mm_set1_epi16(DCT_II_8[1]);
        d[2] = _mm_set1_epi16(DCT_II_8[16]);
        d[3] = _mm_set1_epi16(DCT_II_8[17]);
        d[4] = _mm_set1_epi16(DCT_II_8[32]);
        d[5] = _mm_set1_epi16(DCT_II_8[33]);
        d[6] = _mm_set1_epi16(DCT_II_8[48]);
        d[7] = _mm_set1_epi16(DCT_II_8[49]);

        for(int k=0; k<4; k++){
            d[8+4*k+0] = _mm_set1_epi16(DCT_II_8[8+k]);
            d[8+4*k+1] = _mm_set1_epi16(DCT_II_8[24+k]);
            d[8+4*k+2] = _mm_set1_epi16(DCT_II_8[40+k]);
            d[8+4*k+3] = _mm_set1_epi16(DCT_II_8[56+k]);
        }

        dct2_8x8(x,d,r);

        for(int i=0; i<8; i++){
            __m128i o;

            r[i+0] = _mm_add_epi32(r[i+0], add);
            r[i+8] = _mm_add_epi32(r[i+8], add);

            r[i+0] = _mm_srai_epi32(r[i+0], shift);
            r[i+8] = _mm_srai_epi32(r[i+8], shift);

            o = _mm_packs_epi32(r[i+0], r[i+8]);

            _mm_store_si128((__m128i *) (dst+i*8), o);
        }
        src += 8;
        dst += 64;
    }

    if (num_lines & 0x4){
        __m128i x[10], d[8], r[8];

        x[0] = _mm_unpacklo_epi64(
                _mm_loadl_epi64((__m128i*)(src + 0 * src_stride)),
                _mm_loadl_epi64((__m128i*)(src + 2 * src_stride))
                );
        x[1] = _mm_unpacklo_epi64(
                _mm_loadl_epi64((__m128i*)(src + 4 * src_stride)),
                _mm_loadl_epi64((__m128i*)(src + 6 * src_stride))
                );

        x[2] = _mm_unpacklo_epi64(
                _mm_set1_epi16(src[1 * src_stride + 0]),
                _mm_set1_epi16(src[1 * src_stride + 1])
                );
        x[3] = _mm_unpacklo_epi64(
                _mm_set1_epi16(src[3 * src_stride + 0]),
                _mm_set1_epi16(src[3 * src_stride + 1])
                );
        x[4] = _mm_unpacklo_epi64(
                _mm_set1_epi16(src[5 * src_stride + 0]),
                _mm_set1_epi16(src[5 * src_stride + 1])
                );
        x[5] = _mm_unpacklo_epi64(
                _mm_set1_epi16(src[7 * src_stride + 0]),
                _mm_set1_epi16(src[7 * src_stride + 1])
                );

        x[6] = _mm_unpacklo_epi64(
                _mm_set1_epi16(src[1 * src_stride + 2]),
                _mm_set1_epi16(src[1 * src_stride + 3])
                );
        x[7] = _mm_unpacklo_epi64(
                _mm_set1_epi16(src[3 * src_stride + 2]),
                _mm_set1_epi16(src[3 * src_stride + 3])
                );
        x[8] = _mm_unpacklo_epi64(
                _mm_set1_epi16(src[5 * src_stride + 2]),
                _mm_set1_epi16(src[5 * src_stride + 3])
                );
        x[9] = _mm_unpacklo_epi64(
                _mm_set1_epi16(src[7 * src_stride + 2]),
                _mm_set1_epi16(src[7 * src_stride + 3])
                );

        static const int16_t DCT_II_8_4_sse[8 * 8] = {
            64,  64,  64,  64,  83,  36,  83,  36,
            64, -64,  64, -64,  36, -83,  36, -83,
            64,  64,  64,  64,  36,  83,  36,  83,
            -64, 64, -64,  64, -83,  36, -83,  36,
            89,  75,  18,  50,  89,  75,  18,  50,
            75, -18, -50, -89,  75, -18, -50, -89,
            50, -89,  75,  18,  50, -89,  75,  18,
            18, -50, -89,  75,  18, -50, -89,  75
        };

        d[0] = _mm_load_si128((__m128i*)(DCT_II_8_4_sse + 0));
        d[1] = _mm_load_si128((__m128i*)(DCT_II_8_4_sse + 8));
        d[2] = _mm_load_si128((__m128i*)(DCT_II_8_4_sse + 16));
        d[3] = _mm_load_si128((__m128i*)(DCT_II_8_4_sse + 24));
        d[4] = _mm_load_si128((__m128i*)(DCT_II_8_4_sse + 32));
        d[5] = _mm_load_si128((__m128i*)(DCT_II_8_4_sse + 40));
        d[6] = _mm_load_si128((__m128i*)(DCT_II_8_4_sse + 48));
        d[7] = _mm_load_si128((__m128i*)(DCT_II_8_4_sse + 56));

        dct2_8x4(x, d, r);

        __m128i add = _mm_set1_epi32(1 << (shift - 1));
        for (int i=0; i<8; i+=2){
            __m128i o;

            r[i+0] = _mm_add_epi32(r[i+0], add);
            r[i+1] = _mm_add_epi32(r[i+1], add);

            r[i+0] = _mm_srai_epi32(r[i+0], shift);
            r[i+1] = _mm_srai_epi32(r[i+1], shift);

            o = _mm_packs_epi32(r[i+0], r[i+1]);

            _mm_store_si128((__m128i *) (dst+i/2*8), o);
        }
    }

    if (num_lines & 0x2){
        __m128i x[6];
        __m128i d[6];
        __m128i r[4];

        x[0] = _mm_unpacklo_epi32(
                _mm_loadl_epi64((__m128i*)(src + 0 * src_stride)),
                _mm_loadl_epi64((__m128i*)(src + 2 * src_stride))
                );
        x[0] = _mm_unpacklo_epi16(x[0],x[0]);
        x[1] = _mm_unpacklo_epi32(
                _mm_loadl_epi64((__m128i*)(src + 4 * src_stride)),
                _mm_loadl_epi64((__m128i*)(src + 6 * src_stride))
                );
        x[1] = _mm_unpacklo_epi16(x[1],x[1]);

        x[2] = _mm_unpacklo_epi64(
                _mm_set1_epi16(src[1 * src_stride    ]),
                _mm_set1_epi16(src[1 * src_stride + 1])
                );
        x[3] = _mm_unpacklo_epi64(
                _mm_set1_epi16(src[3 * src_stride    ]),
                _mm_set1_epi16(src[3 * src_stride + 1])
                );
        x[4] = _mm_unpacklo_epi64(
                _mm_set1_epi16(src[5 * src_stride    ]),
                _mm_set1_epi16(src[5 * src_stride + 1])
                );
        x[5] = _mm_unpacklo_epi64(
                _mm_set1_epi16(src[7 * src_stride    ]),
                _mm_set1_epi16(src[7 * src_stride + 1])
                );

        static const int16_t DCT_II_8_2_sse[8 * 6] = {
            64,  64,  64,  64,  83,  36,  83,  36,
            64, -64,  64, -64,  36, -83,  36, -83,
            89,  75,  18,  50,  89,  75,  18,  50,
            75, -18, -50, -89,  75, -18, -50, -89,
            50, -89,  75,  18,  50, -89,  75,  18,
            18, -50, -89,  75,  18, -50, -89,  75
        };

        d[0] = _mm_load_si128((__m128i*)(DCT_II_8_2_sse +  0));
        d[1] = _mm_load_si128((__m128i*)(DCT_II_8_2_sse +  8));
        d[2] = _mm_load_si128((__m128i*)(DCT_II_8_2_sse + 16));
        d[3] = _mm_load_si128((__m128i*)(DCT_II_8_2_sse + 24));
        d[4] = _mm_load_si128((__m128i*)(DCT_II_8_2_sse + 32));
        d[5] = _mm_load_si128((__m128i*)(DCT_II_8_2_sse + 40));

        dct2_8x2(x, d, r);

        __m128i add = _mm_set1_epi32(1 << (shift - 1));
        r[0] = _mm_add_epi32(r[0], add);
        r[1] = _mm_add_epi32(r[1], add);
        r[2] = _mm_add_epi32(r[2], add);
        r[3] = _mm_add_epi32(r[3], add);

        r[0] = _mm_srai_epi32(r[0], shift);
        r[1] = _mm_srai_epi32(r[1], shift);
        r[2] = _mm_srai_epi32(r[2], shift);
        r[3] = _mm_srai_epi32(r[3], shift);

        __m128i o0 = _mm_packs_epi32(r[0], r[1]);
        __m128i o1 = _mm_packs_epi32(r[2], r[3]);

        _mm_store_si128((__m128i *) (dst+0), o0);
        _mm_store_si128((__m128i *) (dst+8), o1);
    }

    if (num_lines & 0x1){
        vvc_inverse_dct_ii_8(src, dst, src_stride, num_lines & 0x1, line_brk, shift);
    }
}

void
vvc_inverse_dct_ii_16_sse(const int16_t *src, int16_t *dst, ptrdiff_t src_stride,
                          int num_lines, int line_brk, int shift)
{
    __m128i add = _mm_set1_epi32(1 << (shift - 1));

    for (int j = 0; j < num_lines / 8; j++) {
        __m128i x[72], d[32], r[32];

        x[0] = _mm_load_si128((__m128i*)(src +  0 * src_stride));
        x[1] = _mm_load_si128((__m128i*)(src +  4 * src_stride));
        x[2] = _mm_load_si128((__m128i*)(src +  8 * src_stride));
        x[3] = _mm_load_si128((__m128i*)(src + 12 * src_stride));
        x[4] = _mm_load_si128((__m128i*)(src +  2 * src_stride));
        x[5] = _mm_load_si128((__m128i*)(src +  6 * src_stride));
        x[6] = _mm_load_si128((__m128i*)(src + 10 * src_stride));
        x[7] = _mm_load_si128((__m128i*)(src + 14 * src_stride));

        for (int k=0; k<8; k++){
            for (int l=0; l<8; l++)
                x[8+8*k+l] = _mm_set1_epi16(src[(2*l+1) * src_stride + k]);
        }

        d[0] = _mm_set1_epi16(DCT_II_16[  0]);
        d[1] = _mm_set1_epi16(DCT_II_16[  1]);
        d[2] = _mm_set1_epi16(DCT_II_16[ 64]);
        d[3] = _mm_set1_epi16(DCT_II_16[ 65]);
        d[4] = _mm_set1_epi16(DCT_II_16[128]);
        d[5] = _mm_set1_epi16(DCT_II_16[129]);
        d[6] = _mm_set1_epi16(DCT_II_16[192]);
        d[7] = _mm_set1_epi16(DCT_II_16[193]);

        for (int k = 0; k < 4; k++) {
            d[8 + 4 * k + 0] = _mm_set1_epi16(DCT_II_16[ 32 + k]);
            d[8 + 4 * k + 1] = _mm_set1_epi16(DCT_II_16[ 96 + k]);
            d[8 + 4 * k + 2] = _mm_set1_epi16(DCT_II_16[160 + k]);
            d[8 + 4 * k + 3] = _mm_set1_epi16(DCT_II_16[224 + k]);
        }

        static const int16_t DCT_II_16_8_sse[8 * 8] = {
            90,  87,  80,  70,  57,  43,  25,   9,
            87,  57,   9, -43, -80, -90, -70, -25,
            80,   9, -70, -87, -25,  57,  90,  43,
            70, -43, -87,   9,  90,  25, -80, -57,
            57, -80, -25,  90,  -9, -87,  43,  70,
            43, -90,  57,  25, -87,  70,   9, -80,
            25, -70,  90, -80,  43,   9, -57,  87,
            9, -25,  43, -57,  70, -80,  87, -90
        };

        for (int k=0; k<8; k++) {
            d[24+k] = _mm_load_si128((__m128i*)(DCT_II_16_8_sse + 8*k));
        }

        dct2_16x8(x, d, r);

        for (int i = 0; i < 32; i+=2) {
            __m128i o;

            r[i + 0] = _mm_add_epi32(r[i + 0], add);
            r[i + 1] = _mm_add_epi32(r[i + 1], add);

            r[i + 0] = _mm_srai_epi32(r[i + 0], shift);
            r[i + 1] = _mm_srai_epi32(r[i + 1], shift);

            o = _mm_packs_epi32(r[i + 0], r[i + 1]);

            _mm_store_si128((__m128i *) (dst + i * 8/2), o);
        }
        src += 8;
        dst += 128;
    }

    if (!(num_lines & 0x7)) return;

    if (num_lines & 0x4){
        __m128i x[42], d[16], r[16];

        x[0] = _mm_unpacklo_epi64(
                _mm_loadl_epi64((__m128i*)(src + 0 * src_stride)),
                _mm_loadl_epi64((__m128i*)(src + 4 * src_stride))
                );
        x[1] = _mm_unpacklo_epi64(
                _mm_loadl_epi64((__m128i*)(src +  8 * src_stride)),
                _mm_loadl_epi64((__m128i*)(src + 12 * src_stride))
                );

        x[2] = _mm_unpacklo_epi64(
                _mm_set1_epi16(src[2 * src_stride + 0]),
                _mm_set1_epi16(src[2 * src_stride + 1])
                );
        x[3] = _mm_unpacklo_epi64(
                _mm_set1_epi16(src[6 * src_stride + 0]),
                _mm_set1_epi16(src[6 * src_stride + 1])
                );
        x[4] = _mm_unpacklo_epi64(
                _mm_set1_epi16(src[10 * src_stride + 0]),
                _mm_set1_epi16(src[10 * src_stride + 1])
                );
        x[5] = _mm_unpacklo_epi64(
                _mm_set1_epi16(src[14 * src_stride + 0]),
                _mm_set1_epi16(src[14 * src_stride + 1])
                );

        x[6] = _mm_unpacklo_epi64(
                _mm_set1_epi16(src[2 * src_stride + 2]),
                _mm_set1_epi16(src[2 * src_stride + 3])
                );
        x[7] = _mm_unpacklo_epi64(
                _mm_set1_epi16(src[6 * src_stride + 2]),
                _mm_set1_epi16(src[6 * src_stride + 3])
                );
        x[8] = _mm_unpacklo_epi64(
                _mm_set1_epi16(src[10 * src_stride + 2]),
                _mm_set1_epi16(src[10 * src_stride + 3])
                );
        x[9] = _mm_unpacklo_epi64(
                _mm_set1_epi16(src[14 * src_stride + 2]),
                _mm_set1_epi16(src[14 * src_stride + 3])
                );

        for (int k=0; k<4; k++){
            for (int l=0; l<8; l++)
                x[10+8*k+l] = _mm_set1_epi16(src[(2*l+1) * src_stride + k]);
        }

        static const int16_t DCT_II_16_4_sse[8 * 16] = {
            64, 64, 64, 64, 83, 36, 83, 36,
            64, -64, 64, -64, 36, -83, 36, -83,
            64, 64, 64, 64, 36, 83, 36, 83,
            -64, 64, -64, 64, -83, 36, -83, 36,
            89, 75, 18, 50, 89, 75, 18, 50,
            75, -18, -50, -89, 75, -18, -50, -89,
            50, -89, 75, 18, 50, -89, 75, 18,
            18, -50, -89, 75, 18, -50, -89, 75,
            90,  87,  80,  70,  57,  43,  25,   9,
            87,  57,   9, -43, -80, -90, -70, -25,
            80,   9, -70, -87, -25,  57,  90,  43,
            70, -43, -87,   9,  90,  25, -80, -57,
            57, -80, -25,  90,  -9, -87,  43,  70,
            43, -90,  57,  25, -87,  70,   9, -80,
            25, -70,  90, -80,  43,   9, -57,  87,
            9, -25,  43, -57,  70, -80,  87, -90
        };

        for (int k=0; k<16; k++) {
            d[k] = _mm_load_si128((__m128i*)(DCT_II_16_4_sse + 8*k));
        }

        dct2_16x4(x, d, r);

        __m128i add = _mm_set1_epi32(1 << (shift - 1));
        for (int i = 0; i < 16; i += 2) {
            __m128i o;

            r[i + 0] = _mm_add_epi32(r[i + 0], add);
            r[i + 1] = _mm_add_epi32(r[i + 1], add);

            r[i + 0] = _mm_srai_epi32(r[i + 0], shift);
            r[i + 1] = _mm_srai_epi32(r[i + 1], shift);

            o = _mm_packs_epi32(r[i + 0], r[i + 1]);

            _mm_store_si128((__m128i *) (dst + i / 2 * 8), o);
        }
    }

    if (num_lines & 0x2){
        __m128i x[24], d[14], r[8];

        x[0] = _mm_unpacklo_epi32(
                _mm_loadl_epi64((__m128i*)(src + 0 * src_stride)),
                _mm_loadl_epi64((__m128i*)(src + 4 * src_stride))
                );
        x[0] = _mm_unpacklo_epi16(x[0], x[0]);
        x[1] = _mm_unpacklo_epi32(
                _mm_loadl_epi64((__m128i*)(src + 8 * src_stride)),
                _mm_loadl_epi64((__m128i*)(src + 12 * src_stride))
                );
        x[1] = _mm_unpacklo_epi16(x[1], x[1]);

        x[2] = _mm_unpacklo_epi64(
                _mm_set1_epi16(src[2 * src_stride]),
                _mm_set1_epi16(src[2 * src_stride + 1])
                );
        x[3] = _mm_unpacklo_epi64(
                _mm_set1_epi16(src[6 * src_stride]),
                _mm_set1_epi16(src[6 * src_stride + 1])
                );
        x[4] = _mm_unpacklo_epi64(
                _mm_set1_epi16(src[10 * src_stride]),
                _mm_set1_epi16(src[10 * src_stride + 1])
                );
        x[5] = _mm_unpacklo_epi64(
                _mm_set1_epi16(src[14 * src_stride]),
                _mm_set1_epi16(src[14 * src_stride + 1])
                );

        for (int k=0; k<2; k++){
            x[6+8*k+0] = _mm_set1_epi16(src[ 1 * src_stride + k]);
            x[6+8*k+1] = _mm_set1_epi16(src[ 3 * src_stride + k]);
            x[6+8*k+2] = _mm_set1_epi16(src[ 5 * src_stride + k]);
            x[6+8*k+3] = _mm_set1_epi16(src[ 7 * src_stride + k]);
            x[6+8*k+4] = _mm_set1_epi16(src[ 9 * src_stride + k]);
            x[6+8*k+5] = _mm_set1_epi16(src[11 * src_stride + k]);
            x[6+8*k+6] = _mm_set1_epi16(src[13 * src_stride + k]);
            x[6+8*k+7] = _mm_set1_epi16(src[15 * src_stride + k]);
        }

        static const int16_t DCT_II_16_2_sse[8 * 14] = {
            64,  64,  64,  64,  83,  36,  83,  36,
            64, -64,  64, -64,  36, -83,  36, -83,
            89,  75,  18,  50,  89,  75,  18,  50,
            75, -18, -50, -89,  75, -18, -50, -89,
            50, -89,  75,  18,  50, -89,  75,  18,
            18, -50, -89,  75,  18, -50, -89,  75,
            90,  87,  80,  70,  57,  43,  25,   9,
            87,  57,   9, -43, -80, -90, -70, -25,
            80,   9, -70, -87, -25,  57,  90,  43,
            70, -43, -87,   9,  90,  25, -80, -57,
            57, -80, -25,  90,  -9, -87,  43,  70,
            43, -90,  57,  25, -87,  70,   9, -80,
            25, -70,  90, -80,  43,   9, -57,  87,
            9, -25,  43, -57,  70, -80,  87, -90
        };

        for (int k=0; k<14; k++) {
            d[k] = _mm_load_si128((__m128i*)(DCT_II_16_2_sse + 8*k));
        }

        dct2_16x2(x, d, r);

        __m128i add = _mm_set1_epi32(1 << (shift - 1));

        for (int i = 0; i < 8; i += 2) {
            __m128i o;

            r[i + 0] =  _mm_add_epi32(r[i + 0], add);
            r[i + 1] =  _mm_add_epi32(r[i + 1], add);

            r[i + 0] = _mm_srai_epi32(r[i + 0], shift);
            r[i + 1] = _mm_srai_epi32(r[i + 1], shift);

            o = _mm_packs_epi32(r[i + 0], r[i + 1]);

            _mm_store_si128((__m128i *) (dst + i / 2 * 8), o);
        }
    }

    if (num_lines & 0x1){
        vvc_inverse_dct_ii_16(src, dst, src_stride, num_lines & 0x1, line_brk, shift);
    }
}

static inline void
dct2_4x8_red1(__m128i *x, __m128i *d, __m128i *r)
{
    __m128i a[4], tmp[4], m[4];
    __m128i z = _mm_setzero_si128();

    m[0] = _mm_unpacklo_epi16(x[0], z);
    m[1] = _mm_unpackhi_epi16(x[0], z);
    m[2] = _mm_unpacklo_epi16(d[0], z);
    m[3] = _mm_unpackhi_epi16(d[0], z);

    a[0] = _mm_madd_epi16(m[0], m[2]);
    a[1] = _mm_madd_epi16(m[1], m[3]);

    m[2] = _mm_unpacklo_epi16(d[1], z);
    m[3] = _mm_unpackhi_epi16(d[1], z);

    a[2] = _mm_madd_epi16(m[0], m[2]);
    a[3] = _mm_madd_epi16(m[1], m[3]);

    tmp[0] = _mm_unpacklo_epi32(a[0], a[2]);
    tmp[1] = _mm_unpackhi_epi32(a[0], a[2]);
    tmp[2] = _mm_unpacklo_epi32(a[2], a[0]);
    tmp[3] = _mm_unpackhi_epi32(a[2], a[0]);

    r[0] = _mm_unpacklo_epi64(tmp[0], tmp[2]);
    r[1] = _mm_unpackhi_epi64(tmp[0], tmp[2]);
    r[2] = _mm_unpacklo_epi64(tmp[1], tmp[3]);
    r[3] = _mm_unpackhi_epi64(tmp[1], tmp[3]);

    tmp[0] = _mm_unpacklo_epi32(a[1], a[3]);
    tmp[1] = _mm_unpackhi_epi32(a[1], a[3]);
    tmp[2] = _mm_unpacklo_epi32(a[3], a[1]);
    tmp[3] = _mm_unpackhi_epi32(a[3], a[1]);

    r[4] = _mm_unpacklo_epi64(tmp[0], tmp[2]);
    r[5] = _mm_unpackhi_epi64(tmp[0], tmp[2]);
    r[6] = _mm_unpacklo_epi64(tmp[1], tmp[3]);
    r[7] = _mm_unpackhi_epi64(tmp[1], tmp[3]);
}

static inline void
dct2_8x8_red2(__m128i *x, __m128i *d, __m128i *r)
{
    __m128i e[8], o[8], oo[8];
    __m128i m[8], z;//a[4],

    dct2_4x8_red1(x+0, d+0, e);

    __m128i tmp[4];
    z = _mm_setzero_si128();

    m[0] = _mm_unpacklo_epi16(x[4], z);
    m[1] = _mm_unpackhi_epi16(x[4], z);
    m[2] = _mm_unpacklo_epi16(d[8], z);
    m[3] = _mm_unpackhi_epi16(d[8], z);

    oo[0] = _mm_madd_epi16(m[0], m[2]);
    oo[1] = _mm_madd_epi16(m[1], m[3]);

    m[2] = _mm_unpacklo_epi16(d[12], z);
    m[3] = _mm_unpackhi_epi16(d[12], z);

    oo[2] = _mm_madd_epi16(m[0], m[2]);
    oo[3] = _mm_madd_epi16(m[1], m[3]);

    m[2] = _mm_unpacklo_epi16(d[16], z);
    m[3] = _mm_unpackhi_epi16(d[16], z);

    oo[4] = _mm_madd_epi16(m[0], m[2]);
    oo[5] = _mm_madd_epi16(m[1], m[3]);

    m[2] = _mm_unpacklo_epi16(d[20], z);
    m[3] = _mm_unpackhi_epi16(d[20], z);

    oo[6] = _mm_madd_epi16(m[0], m[2]);
    oo[7] = _mm_madd_epi16(m[1], m[3]);

    tmp[0] = _mm_unpacklo_epi32(oo[0], oo[2]);
    tmp[1] = _mm_unpackhi_epi32(oo[0], oo[2]);
    tmp[2] = _mm_unpacklo_epi32(oo[4], oo[6]);
    tmp[3] = _mm_unpackhi_epi32(oo[4], oo[6]);

    o[0] = _mm_unpacklo_epi64(tmp[0], tmp[2]);
    o[1] = _mm_unpackhi_epi64(tmp[0], tmp[2]);
    o[2] = _mm_unpacklo_epi64(tmp[1], tmp[3]);
    o[3] = _mm_unpackhi_epi64(tmp[1], tmp[3]);

    tmp[0] = _mm_unpacklo_epi32(oo[1], oo[3]);
    tmp[1] = _mm_unpackhi_epi32(oo[1], oo[3]);
    tmp[2] = _mm_unpacklo_epi32(oo[5], oo[7]);
    tmp[3] = _mm_unpackhi_epi32(oo[5], oo[7]);

    o[4] = _mm_unpacklo_epi64(tmp[0], tmp[2]);
    o[5] = _mm_unpackhi_epi64(tmp[0], tmp[2]);
    o[6] = _mm_unpacklo_epi64(tmp[1], tmp[3]);
    o[7] = _mm_unpackhi_epi64(tmp[1], tmp[3]);

    for(int k=0; k<8; k++){
        r[  k] = _mm_add_epi32(e[k], o[k]);
        r[8+k] = _mm_shuffle_epi32(_mm_sub_epi32(e[k],  o[k]),  0x1B);
    }
}

static inline void
dct2_16x8_red4(__m128i *x, __m128i *d, __m128i *r)
{
    __m128i e[16], o[16];
    __m128i *evn_x = x;
    __m128i *odd_x = x + 8;
    __m128i *odd_d = d + 24;

    dct2_8x8_red2(evn_x, d+0, e);

    for(int k=0; k<8; k++){
        __m128i m[4], a[2];
        m[ 0] = _mm_unpacklo_epi16(odd_x[2*k], odd_x[2*k + 1]);
        m[ 1] = _mm_unpackhi_epi16(odd_x[2*k], odd_x[2*k + 1]);
        m[ 2] = _mm_unpacklo_epi16(odd_d[0], odd_d[1]);
        m[ 3] = _mm_unpackhi_epi16(odd_d[0], odd_d[1]);

        a[0] = _mm_madd_epi16(m[0], m[2]);
        a[1] = _mm_madd_epi16(m[1], m[3]);

        o[2*k    ] = a[0];
        o[2*k + 1] = a[1];
    #if 0
    }

    for(int k=0; k<8; k++){
    #endif
        r[4*k+0] = _mm_add_epi32(e[k+0], o[2*k+0]);
        r[4*k+1] = _mm_add_epi32(e[k+8], o[2*k+1]);
        r[4*k+2] = _mm_shuffle_epi32(_mm_sub_epi32(e[k+8], o[2*k+1]), 0x1B);
        r[4*k+3] = _mm_shuffle_epi32(_mm_sub_epi32(e[k+0], o[2*k+0]), 0x1B);
    }
}

static inline void
dct2_32x8_red8(__m128i *x, __m128i *d, __m128i *r){
    __m128i e[32], o[32];
    __m128i *evn_x = x;
    __m128i *odd_x = x + 24;
    __m128i *odd_d1 = d + 32;
    __m128i *odd_d2 = d + 48;

    /*FIXME we can also use a reduced version of dct2_16x8*/
    dct2_16x8_red4(evn_x, d+0, e+0);

    #if 0
    for(int k=0; k<8; k++){
        matMult16x16_red4(odd_x + 4*k, odd_d1, o+4*k+0);
        matMult16x16_red4(odd_x + 4*k, odd_d2, o+4*k+2);
    }

    /*FIXME find out which e and o results are null*/
    for(int k=0; k<8; k++){
    #else
    for(int k=0; k<8; k++){
        matMult16x16_red4(odd_x + 4*k, odd_d1, o+4*k+0);
        matMult16x16_red4(odd_x + 4*k, odd_d2, o+4*k+2);
    #endif
        r[8*k+0] = _mm_add_epi32(e[4*k+0], o[4*k+0]);
        r[8*k+1] = _mm_add_epi32(e[4*k+1], o[4*k+1]);
        r[8*k+2] = _mm_add_epi32(e[4*k+2], o[4*k+2]);
        r[8*k+3] = _mm_add_epi32(e[4*k+3], o[4*k+3]);
        r[8*k+4] = _mm_shuffle_epi32(_mm_sub_epi32(e[4*k+3], o[4*k+3]), 0x1b);
        r[8*k+5] = _mm_shuffle_epi32(_mm_sub_epi32(e[4*k+2], o[4*k+2]), 0x1b);
        r[8*k+6] = _mm_shuffle_epi32(_mm_sub_epi32(e[4*k+1], o[4*k+1]), 0x1b);
        r[8*k+7] = _mm_shuffle_epi32(_mm_sub_epi32(e[4*k+0], o[4*k+0]), 0x1b);
    }
}

static inline void
load_src_32x8_red8(const int16_t *src, ptrdiff_t src_stride, __m128i *x)
{
    __m128i *evn_xl = x;
    __m128i *odd_xl = x + 24;
    int k;

    evn_xl[0] = _mm_load_si128((__m128i*)(src +  0 * src_stride));
    #if 0
    evn_xl[1] = _mm_setzero_si128();
    evn_xl[2] = _mm_setzero_si128();
    evn_xl[3] = _mm_setzero_si128();
    #endif
    evn_xl[4] = _mm_load_si128((__m128i*)(src +  4 * src_stride));
    #if 0
    evn_xl[5] = _mm_setzero_si128();
    evn_xl[6] = _mm_setzero_si128();
    evn_xl[7] = _mm_setzero_si128();
    #endif

    for (k = 0; k < 8; ++k){
        const int16_t *src_col = src + k;
        evn_xl[ 8] = _mm_set1_epi16(src_col[ 2 * src_stride]);
        evn_xl[ 9] = _mm_set1_epi16(src_col[ 6 * src_stride]);
        evn_xl += 2;
    }

    for (k = 0; k < 8; ++k){
        const int16_t *src_col = src + k;
        odd_xl[ 0] = _mm_set1_epi16(src_col[ 1 * src_stride]);
        odd_xl[ 1] = _mm_set1_epi16(src_col[ 3 * src_stride]);
        odd_xl[ 2] = _mm_set1_epi16(src_col[ 5 * src_stride]);
        odd_xl[ 3] = _mm_set1_epi16(src_col[ 7 * src_stride]);
        odd_xl += 4;
    }
}

static inline void
dct2_4x8_red2(__m128i *x, __m128i *d, __m128i *r)
{
    __m128i a[8], tmp[4], m[8];
    __m128i z = _mm_setzero_si128();

    m[0] = _mm_unpacklo_epi16(x[0], z);
    m[1] = _mm_unpackhi_epi16(x[0], z);
    m[2] = _mm_unpacklo_epi16(d[0], z);
    m[3] = _mm_unpackhi_epi16(d[0], z);

    a[0] = _mm_madd_epi16(m[0], m[2]);
    a[1] = _mm_madd_epi16(m[1], m[3]);

    m[2] = _mm_unpacklo_epi16(d[1], z);
    m[3] = _mm_unpackhi_epi16(d[1], z);

    a[2] = _mm_madd_epi16(m[0], m[2]);
    a[3] = _mm_madd_epi16(m[1], m[3]);

    m[0] = _mm_unpacklo_epi16(x[1], z);
    m[1] = _mm_unpackhi_epi16(x[1], z);
    m[2] = _mm_unpacklo_epi16(d[2], z);
    m[3] = _mm_unpackhi_epi16(d[2], z);

    a[4] = _mm_madd_epi16(m[0], m[2]);
    a[5] = _mm_madd_epi16(m[1], m[3]);

    m[2] = _mm_unpacklo_epi16(d[3], z);
    m[3] = _mm_unpackhi_epi16(d[3], z);

    a[6] = _mm_madd_epi16(m[0], m[2]);
    a[7] = _mm_madd_epi16(m[1], m[3]);

    m[0] = _mm_add_epi32(a[0], a[4]);
    m[1] = _mm_add_epi32(a[1], a[5]);
    m[2] = _mm_add_epi32(a[2], a[6]);
    m[3] = _mm_add_epi32(a[3], a[7]);
    m[4] = _mm_sub_epi32(a[2], a[6]);
    m[5] = _mm_sub_epi32(a[3], a[7]);
    m[6] = _mm_sub_epi32(a[0], a[4]);
    m[7] = _mm_sub_epi32(a[1], a[5]);

    tmp[0] = _mm_unpacklo_epi32(m[0], m[2]);
    tmp[1] = _mm_unpackhi_epi32(m[0], m[2]);
    tmp[2] = _mm_unpacklo_epi32(m[4], m[6]);
    tmp[3] = _mm_unpackhi_epi32(m[4], m[6]);

    r[0] = _mm_unpacklo_epi64(tmp[0], tmp[2]);
    r[1] = _mm_unpackhi_epi64(tmp[0], tmp[2]);
    r[2] = _mm_unpacklo_epi64(tmp[1], tmp[3]);
    r[3] = _mm_unpackhi_epi64(tmp[1], tmp[3]);

    tmp[0] = _mm_unpacklo_epi32(m[1], m[3]);
    tmp[1] = _mm_unpackhi_epi32(m[1], m[3]);
    tmp[2] = _mm_unpacklo_epi32(m[5], m[7]);
    tmp[3] = _mm_unpackhi_epi32(m[5], m[7]);

    r[4] = _mm_unpacklo_epi64(tmp[0], tmp[2]);
    r[5] = _mm_unpackhi_epi64(tmp[0], tmp[2]);
    r[6] = _mm_unpacklo_epi64(tmp[1], tmp[3]);
    r[7] = _mm_unpackhi_epi64(tmp[1], tmp[3]);
}

static inline void
dct2_8x8_red4(__m128i *x, __m128i *d, __m128i *r)
{
    __m128i e[8], o[8], oo[8];
    __m128i m[8];

    dct2_4x8_red2(x+0, d+0, e);

    __m128i tmp[4];

    m[0] = _mm_unpacklo_epi16(x[4], x[5]);
    m[1] = _mm_unpackhi_epi16(x[4], x[5]);
    m[2] = _mm_unpacklo_epi16(d[8], d[8+1]);
    m[3] = _mm_unpackhi_epi16(d[8], d[8+1]);

    oo[0] = _mm_madd_epi16(m[0], m[2]);
    oo[1] = _mm_madd_epi16(m[1], m[3]);

    m[2] = _mm_unpacklo_epi16(d[12], d[12+1]);
    m[3] = _mm_unpackhi_epi16(d[12], d[12+1]);

    oo[2] = _mm_madd_epi16(m[0], m[2]);
    oo[3] = _mm_madd_epi16(m[1], m[3]);

    m[2] = _mm_unpacklo_epi16(d[16], d[16 + 1]);
    m[3] = _mm_unpackhi_epi16(d[16], d[16 + 1]);

    oo[4] = _mm_madd_epi16(m[0], m[2]);
    oo[5] = _mm_madd_epi16(m[1], m[3]);

    m[2] = _mm_unpacklo_epi16(d[20], d[20 + 1]);
    m[3] = _mm_unpackhi_epi16(d[20], d[20 + 1]);

    oo[6] = _mm_madd_epi16(m[0], m[2]);
    oo[7] = _mm_madd_epi16(m[1], m[3]);

    tmp[0] = _mm_unpacklo_epi32(oo[0], oo[2]);
    tmp[1] = _mm_unpackhi_epi32(oo[0], oo[2]);
    tmp[2] = _mm_unpacklo_epi32(oo[4], oo[6]);
    tmp[3] = _mm_unpackhi_epi32(oo[4], oo[6]);

    o[0] = _mm_unpacklo_epi64(tmp[0], tmp[2]);
    o[1] = _mm_unpackhi_epi64(tmp[0], tmp[2]);
    o[2] = _mm_unpacklo_epi64(tmp[1], tmp[3]);
    o[3] = _mm_unpackhi_epi64(tmp[1], tmp[3]);

    tmp[0] = _mm_unpacklo_epi32(oo[1], oo[3]);
    tmp[1] = _mm_unpackhi_epi32(oo[1], oo[3]);
    tmp[2] = _mm_unpacklo_epi32(oo[5], oo[7]);
    tmp[3] = _mm_unpackhi_epi32(oo[5], oo[7]);

    o[4] = _mm_unpacklo_epi64(tmp[0], tmp[2]);
    o[5] = _mm_unpackhi_epi64(tmp[0], tmp[2]);
    o[6] = _mm_unpacklo_epi64(tmp[1], tmp[3]);
    o[7] = _mm_unpackhi_epi64(tmp[1], tmp[3]);

    for(int k=0; k<8; k++){
        r[  k] = _mm_add_epi32(e[k], o[k]);
        r[8+k] = _mm_shuffle_epi32(_mm_sub_epi32(e[k],  o[k]),  0x1B);
    }
}

static inline void
dct2_16x8_red8(__m128i *x, __m128i *d, __m128i *r)
{
    __m128i e[16], o[16];
    __m128i *evn_x = x;
    __m128i *odd_x = x + 8;
    __m128i *odd_d = d + 24;

    dct2_8x8_red4(evn_x, d+0, e);

    for(int k=0; k<8; k++){
        __m128i m[8], a[4];
        m[ 0] = _mm_unpacklo_epi16(odd_x[4 * k], odd_x[4 * k + 1]);
        m[ 1] = _mm_unpackhi_epi16(odd_x[4 * k], odd_x[4 * k + 1]);
        m[ 2] = _mm_unpacklo_epi16(odd_d[0], odd_d[1]);
        m[ 3] = _mm_unpackhi_epi16(odd_d[0], odd_d[1]);
        m[ 4] = _mm_unpacklo_epi16(odd_x[4 * k + 2], odd_x[4 * k + 3]);
        m[ 5] = _mm_unpackhi_epi16(odd_x[4 * k + 2], odd_x[4 * k + 3]);
        m[ 6] = _mm_unpacklo_epi16(odd_d[2], odd_d[3]);
        m[ 7] = _mm_unpackhi_epi16(odd_d[2], odd_d[3]);

        a[0] = _mm_madd_epi16(m[ 0], m[ 2]);
        a[1] = _mm_madd_epi16(m[ 1], m[ 3]);
        a[2] = _mm_madd_epi16(m[ 4], m[ 6]);
        a[3] = _mm_madd_epi16(m[ 5], m[ 7]);

        o[2*k    ] = _mm_add_epi32(a[0], a[2]);
        o[2*k + 1] = _mm_add_epi32(a[1], a[3]);
    #if 0
    }

    for(int k=0; k<8; k++){
    #endif
        r[4*k+0] = _mm_add_epi32(e[k+0], o[2*k+0]);
        r[4*k+1] = _mm_add_epi32(e[k+8], o[2*k+1]);
        r[4*k+2] = _mm_shuffle_epi32(_mm_sub_epi32(e[k+8], o[2*k+1]), 0x1B);
        r[4*k+3] = _mm_shuffle_epi32(_mm_sub_epi32(e[k+0], o[2*k+0]), 0x1B);
    }
}

static inline void
dct2_32x8_red16(__m128i *x, __m128i *d, __m128i *r){
    __m128i e[32], o[32];
    __m128i *evn_x = x;
    __m128i *odd_x = x + 40;
    __m128i *odd_d1 = d + 32;
    __m128i *odd_d2 = d + 48;

    dct2_16x8_red8(evn_x, d+0, e+0);

    for(int k=0; k<8; k++){
        matMult16x16_red8(odd_x + 8*k, odd_d1, o+4*k+0);
        matMult16x16_red8(odd_x + 8*k, odd_d2, o+4*k+2);
    #if 0
    }

    for(int k=0; k<8; k++){
    #endif
        r[8*k+0] = _mm_add_epi32(e[4*k+0], o[4*k+0]);
        r[8*k+1] = _mm_add_epi32(e[4*k+1], o[4*k+1]);
        r[8*k+2] = _mm_add_epi32(e[4*k+2], o[4*k+2]);
        r[8*k+3] = _mm_add_epi32(e[4*k+3], o[4*k+3]);
        r[8*k+4] = _mm_shuffle_epi32(_mm_sub_epi32(e[4*k+3], o[4*k+3]), 0x1b);
        r[8*k+5] = _mm_shuffle_epi32(_mm_sub_epi32(e[4*k+2], o[4*k+2]), 0x1b);
        r[8*k+6] = _mm_shuffle_epi32(_mm_sub_epi32(e[4*k+1], o[4*k+1]), 0x1b);
        r[8*k+7] = _mm_shuffle_epi32(_mm_sub_epi32(e[4*k+0], o[4*k+0]), 0x1b);
    }
}

static inline void
load_src_32x8_red16(const int16_t *src, ptrdiff_t src_stride, __m128i *x)
{
    __m128i *evn_xl = x;
    __m128i *odd_xl = x + 40;
    int k;

    evn_xl[0] = _mm_load_si128((__m128i*)(src +  0 * src_stride));
    evn_xl[1] = _mm_load_si128((__m128i*)(src +  8 * src_stride));
    evn_xl[2] = _mm_setzero_si128();
    evn_xl[3] = _mm_setzero_si128();
    evn_xl[4] = _mm_load_si128((__m128i*)(src +  4 * src_stride));
    evn_xl[5] = _mm_load_si128((__m128i*)(src + 12 * src_stride));
    evn_xl[6] = _mm_setzero_si128();
    evn_xl[7] = _mm_setzero_si128();

    for (k = 0; k < 8; ++k){
        const int16_t *src_col = src + k;
        evn_xl[ 8] = _mm_set1_epi16(src_col[ 2 * src_stride]);
        evn_xl[ 9] = _mm_set1_epi16(src_col[ 6 * src_stride]);
        evn_xl[10] = _mm_set1_epi16(src_col[10 * src_stride]);
        evn_xl[11] = _mm_set1_epi16(src_col[14 * src_stride]);
        evn_xl += 4;
    }

    for (k = 0; k < 8; ++k){
        const int16_t *src_col = src + k;
        odd_xl[ 0] = _mm_set1_epi16(src_col[ 1 * src_stride]);
        odd_xl[ 1] = _mm_set1_epi16(src_col[ 3 * src_stride]);
        odd_xl[ 2] = _mm_set1_epi16(src_col[ 5 * src_stride]);
        odd_xl[ 3] = _mm_set1_epi16(src_col[ 7 * src_stride]);
        odd_xl[ 4] = _mm_set1_epi16(src_col[ 9 * src_stride]);
        odd_xl[ 5] = _mm_set1_epi16(src_col[11 * src_stride]);
        odd_xl[ 6] = _mm_set1_epi16(src_col[13 * src_stride]);
        odd_xl[ 7] = _mm_set1_epi16(src_col[15 * src_stride]);
        odd_xl += 8;
    }
}

#if 0
static inline void
load_src_32x8_24(const int16_t *src, ptrdiff_t src_stride, __m128i *x)
{
    __m128i *evn_xl = x;
    __m128i *odd_xl = x + 72;
    int k;

    evn_xl[0] = _mm_load_si128((__m128i*)(src +  0 * src_stride));
    evn_xl[1] = _mm_load_si128((__m128i*)(src +  8 * src_stride));
    evn_xl[2] = _mm_load_si128((__m128i*)(src + 16 * src_stride));
    evn_xl[3] = _mm_setzero_si128();
    evn_xl[4] = _mm_load_si128((__m128i*)(src +  4 * src_stride));
    evn_xl[5] = _mm_load_si128((__m128i*)(src + 12 * src_stride));
    evn_xl[6] = _mm_load_si128((__m128i*)(src + 20 * src_stride));
    evn_xl[7] = _mm_setzero_si128();

    for (k = 0; k < 8; ++k){
        const uint16_t *src_col = src + k;
        evn_xl[ 8] = _mm_set1_epi16(src_col[ 2 * src_stride]);
        evn_xl[ 9] = _mm_set1_epi16(src_col[ 6 * src_stride]);
        evn_xl[10] = _mm_set1_epi16(src_col[10 * src_stride]);
        evn_xl[11] = _mm_set1_epi16(src_col[14 * src_stride]);
        evn_xl[12] = _mm_set1_epi16(src_col[18 * src_stride]);
        evn_xl[13] = _mm_set1_epi16(src_col[22 * src_stride]);
        evn_xl[14] = _mm_setzero_si128();
        evn_xl[15] = _mm_setzero_si128();
        evn_xl += 8;
    }

    for (k = 0; k < 8; ++k){
        const uint16_t *src_col = src + k;
        odd_xl[ 0] = _mm_set1_epi16(src_col[ 1 * src_stride]);
        odd_xl[ 1] = _mm_set1_epi16(src_col[ 3 * src_stride]);
        odd_xl[ 2] = _mm_set1_epi16(src_col[ 5 * src_stride]);
        odd_xl[ 3] = _mm_set1_epi16(src_col[ 7 * src_stride]);
        odd_xl[ 4] = _mm_set1_epi16(src_col[ 9 * src_stride]);
        odd_xl[ 5] = _mm_set1_epi16(src_col[11 * src_stride]);
        odd_xl[ 6] = _mm_set1_epi16(src_col[13 * src_stride]);
        odd_xl[ 7] = _mm_set1_epi16(src_col[15 * src_stride]);
        odd_xl[ 8] = _mm_set1_epi16(src_col[17 * src_stride]);
        odd_xl[ 9] = _mm_set1_epi16(src_col[19 * src_stride]);
        odd_xl[10] = _mm_set1_epi16(src_col[21 * src_stride]);
        odd_xl[11] = _mm_set1_epi16(src_col[23 * src_stride]);
        odd_xl[12] = _mm_setzero_si128();
        odd_xl[13] = _mm_setzero_si128();
        odd_xl[14] = _mm_setzero_si128();
        odd_xl[15] = _mm_setzero_si128();
        odd_xl += 16;
    }
}
#endif

static inline void
load_src_32x8(const int16_t *src, ptrdiff_t src_stride, __m128i *x)
{
    __m128i *evn_xl = x;
    __m128i *odd_xl = x + 72;
    int k;

    evn_xl[0] = _mm_load_si128((__m128i*)(src +  0 * src_stride));
    evn_xl[1] = _mm_load_si128((__m128i*)(src +  8 * src_stride));
    evn_xl[2] = _mm_load_si128((__m128i*)(src + 16 * src_stride));
    evn_xl[3] = _mm_load_si128((__m128i*)(src + 24 * src_stride));
    evn_xl[4] = _mm_load_si128((__m128i*)(src +  4 * src_stride));
    evn_xl[5] = _mm_load_si128((__m128i*)(src + 12 * src_stride));
    evn_xl[6] = _mm_load_si128((__m128i*)(src + 20 * src_stride));
    evn_xl[7] = _mm_load_si128((__m128i*)(src + 28 * src_stride));

    for (k = 0; k < 8; ++k){
        const int16_t *src_col = src + k;
        evn_xl[ 8] = _mm_set1_epi16(src_col[ 2 * src_stride]);
        evn_xl[ 9] = _mm_set1_epi16(src_col[ 6 * src_stride]);
        evn_xl[10] = _mm_set1_epi16(src_col[10 * src_stride]);
        evn_xl[11] = _mm_set1_epi16(src_col[14 * src_stride]);
        evn_xl[12] = _mm_set1_epi16(src_col[18 * src_stride]);
        evn_xl[13] = _mm_set1_epi16(src_col[22 * src_stride]);
        evn_xl[14] = _mm_set1_epi16(src_col[26 * src_stride]);
        evn_xl[15] = _mm_set1_epi16(src_col[30 * src_stride]);
        evn_xl += 8;
    }

    for (k = 0; k < 8; ++k){
        const int16_t *src_col = src + k;
        odd_xl[ 0] = _mm_set1_epi16(src_col[ 1 * src_stride]);
        odd_xl[ 1] = _mm_set1_epi16(src_col[ 3 * src_stride]);
        odd_xl[ 2] = _mm_set1_epi16(src_col[ 5 * src_stride]);
        odd_xl[ 3] = _mm_set1_epi16(src_col[ 7 * src_stride]);
        odd_xl[ 4] = _mm_set1_epi16(src_col[ 9 * src_stride]);
        odd_xl[ 5] = _mm_set1_epi16(src_col[11 * src_stride]);
        odd_xl[ 6] = _mm_set1_epi16(src_col[13 * src_stride]);
        odd_xl[ 7] = _mm_set1_epi16(src_col[15 * src_stride]);
        odd_xl[ 8] = _mm_set1_epi16(src_col[17 * src_stride]);
        odd_xl[ 9] = _mm_set1_epi16(src_col[19 * src_stride]);
        odd_xl[10] = _mm_set1_epi16(src_col[21 * src_stride]);
        odd_xl[11] = _mm_set1_epi16(src_col[23 * src_stride]);
        odd_xl[12] = _mm_set1_epi16(src_col[25 * src_stride]);
        odd_xl[13] = _mm_set1_epi16(src_col[27 * src_stride]);
        odd_xl[14] = _mm_set1_epi16(src_col[29 * src_stride]);
        odd_xl[15] = _mm_set1_epi16(src_col[31 * src_stride]);
        odd_xl += 16;
    }
}

static inline void
load_dct_32x8(__m128i *d)
{
    int k;
    d[0] = _mm_set1_epi16(DCT_II_32[0]);
    d[1] = _mm_set1_epi16(DCT_II_32[1]);
    d[2] = _mm_set1_epi16(DCT_II_32[256]);
    d[3] = _mm_set1_epi16(DCT_II_32[257]);
    d[4] = _mm_set1_epi16(DCT_II_32[512]);
    d[5] = _mm_set1_epi16(DCT_II_32[513]);
    d[6] = _mm_set1_epi16(DCT_II_32[768]);
    d[7] = _mm_set1_epi16(DCT_II_32[769]);

    for (k = 0; k < 4; k++) {
        d[8 + 4 * k + 0] = _mm_set1_epi16(DCT_II_32[128 + k]);
        d[8 + 4 * k + 1] = _mm_set1_epi16(DCT_II_32[384 + k]);
        d[8 + 4 * k + 2] = _mm_set1_epi16(DCT_II_32[640 + k]);
        d[8 + 4 * k + 3] = _mm_set1_epi16(DCT_II_32[896 + k]);
    }

    static const int16_t DCT_II_32_8_sse[8 * 40] = {
        90,  87,  80,  70,  57,  43,  25,   9,
        87,  57,   9, -43, -80, -90, -70, -25,
        80,   9, -70, -87, -25,  57,  90,  43,
        70, -43, -87,   9,  90,  25, -80, -57,
        57, -80, -25,  90,  -9, -87,  43,  70,
        43, -90,  57,  25, -87,  70,   9, -80,
        25, -70,  90, -80,  43,   9, -57,  87,
        9, -25,  43, -57,  70, -80,  87, -90,

        90, 90, 88, 85, 82, 78, 73, 67,
        90, 82, 67, 46, 22, -4, -31, -54,
        88, 67, 31, -13, -54, -82, -90, -78,
        85, 46, -13, -67, -90, -73, -22, 38,
        82, 22, -54, -90, -61, 13, 78, 85,
        78, -4, -82, -73, 13, 85, 67, -22,
        73, -31, -90, -22, 78, 67, -38, -90,
        67, -54, -78, 38, 85, -22, -90, 4,
        61, -73, -46, 82, 31, -88, -13, 90,
        54, -85, -4, 88, -46, -61, 82, 13,
        46, -90, 38, 54, -90, 31, 61, -88,
        38, -88, 73, -4, -67, 90, -46, -31,
        31, -78, 90, -61, 4, 54, -88, 82,
        22, -61, 85, -90, 73, -38, -4, 46,
        13, -38, 61, -78, 88, -90, 85, -73,
        4, -13, 22, -31, 38, -46, 54, -61,

        61, 54, 46, 38, 31, 22, 13, 4,
        -73, -85, -90, -88, -78, -61, -38, -13,
        -46, -4, 38, 73, 90, 85, 61, 22,
        82, 88, 54, -4, -61, -90, -78, -31,
        31, -46, -90, -67, 4, 73, 88, 38,
        -88, -61, 31, 90, 54, -38, -90, -46,
        -13, 82, 61, -46, -88, -4, 85, 54,
        90, 13, -88, -31, 82, 46, -73, -61,
        -4, -90, 22, 85, -38, -78, 54, 67,
        -90, 38, 67, -78, -22, 90, -31, -73,
        22, 67, -85, 13, 73, -82, 4, 78,
        85, -78, 13, 61, -90, 54, 22, -82,
        -38, -22, 73, -90, 67, -13, -46, 85,
        -78, 90, -82, 54, -13, -31, 67, -88,
        54, -31, 4, 22, -46, 67, -82, 90,
        67, -73, 78, -82, 85, -88, 90, -90
    };

    for (int k=0; k<40; k++) {
        d[24+k] = _mm_load_si128((__m128i*)(DCT_II_32_8_sse + 8*k));
    }

}

static void
dct2_32_8lines_red8(const int16_t *src, int16_t *dst, ptrdiff_t src_stride,
               int num_lines, int line_brk, int shift)
{
    int rnd_add = 1 << (shift - 1);
    __m128i d[64];
    /*FIXME we can still reduce dct coefficient loading operation
      however it is less intensive than src_load since it is done
      outside of loop */
    load_dct_32x8(d);
    for (int j = 0; j < num_lines / 8; j++) {
        __m128i x[56], r[64], rnd_add_v;

        load_src_32x8_red8(src, src_stride, x);

        dct2_32x8_red8(x, d, r);

        rnd_add_v = _mm_set1_epi32(rnd_add);

        for (int i = 0; i < 64; i+=2) {
            __m128i o;

            r[i + 0] = _mm_add_epi32(r[i + 0], rnd_add_v);
            r[i + 1] = _mm_add_epi32(r[i + 1], rnd_add_v);

            r[i + 0] = _mm_srai_epi32(r[i + 0], shift);
            r[i + 1] = _mm_srai_epi32(r[i + 1], shift);

            o = _mm_packs_epi32(r[i + 0], r[i + 1]);

            _mm_store_si128((__m128i *) (dst + i * 8/2), o);
        }
        src += 8;
        dst += 256;
    }
}

static void
dct2_32_8lines_red16(const int16_t *src, int16_t *dst, ptrdiff_t src_stride,
                    int num_lines, int line_brk, int shift)
{
    int rnd_add = 1 << (shift - 1);
    __m128i d[64];
    load_dct_32x8(d);
    for (int j = 0; j < num_lines / 8; j++) {
        __m128i x[104], r[64], rnd_add_v;

        load_src_32x8_red16(src, src_stride, x);

        dct2_32x8_red16(x, d, r);

        rnd_add_v = _mm_set1_epi32(rnd_add);

        for (int i = 0; i < 64; i+=2) {
            __m128i o;

            r[i + 0] = _mm_add_epi32(r[i + 0], rnd_add_v);
            r[i + 1] = _mm_add_epi32(r[i + 1], rnd_add_v);

            r[i + 0] = _mm_srai_epi32(r[i + 0], shift);
            r[i + 1] = _mm_srai_epi32(r[i + 1], shift);

            o = _mm_packs_epi32(r[i + 0], r[i + 1]);

            _mm_store_si128((__m128i *) (dst + i * 8/2), o);
        }
        src += 8;
        dst += 256;
    }
}

static void
dct2_32_8lines(const int16_t *src, int16_t *dst, ptrdiff_t src_stride,
               int num_lines, int line_brk, int shift)
{
    int rnd_add = 1 << (shift - 1);
    __m128i d[64];
    load_dct_32x8(d);
    for (int j = 0; j < num_lines / 8; j++) {
        __m128i x[200], r[64], rnd_add_v;

        load_src_32x8(src, src_stride, x);

        dct2_32x8(x, d, r);

        rnd_add_v = _mm_set1_epi32(rnd_add);

        for (int i = 0; i < 64; i+=2) {
            __m128i o;

            r[i + 0] = _mm_add_epi32(r[i + 0], rnd_add_v);
            r[i + 1] = _mm_add_epi32(r[i + 1], rnd_add_v);

            r[i + 0] = _mm_srai_epi32(r[i + 0], shift);
            r[i + 1] = _mm_srai_epi32(r[i + 1], shift);

            o = _mm_packs_epi32(r[i + 0], r[i + 1]);

            _mm_store_si128((__m128i *) (dst + i * 8/2), o);
        }
        src += 8;
        dst += 256;
    }
}

static inline void
load_dct_32x4(__m128i *d)
{
    static const int16_t DCT_II_32_4_sse[8 * 48] = {
        64, 64, 64, 64, 83, 36, 83, 36,
        64, -64, 64, -64, 36, -83, 36, -83,
        64, 64, 64, 64, 36, 83, 36, 83,
        -64, 64, -64, 64, -83, 36, -83, 36,

        89, 75, 18, 50, 89, 75, 18, 50,
        75, -18, -50, -89, 75, -18, -50, -89,
        50, -89, 75, 18, 50, -89, 75, 18,
        18, -50, -89, 75, 18, -50, -89, 75,

        90,  87,  80,  70,  57,  43,  25,   9,
        87,  57,   9, -43, -80, -90, -70, -25,
        80,   9, -70, -87, -25,  57,  90,  43,
        70, -43, -87,   9,  90,  25, -80, -57,
        57, -80, -25,  90,  -9, -87,  43,  70,
        43, -90,  57,  25, -87,  70,   9, -80,
        25, -70,  90, -80,  43,   9, -57,  87,
        9, -25,  43, -57,  70, -80,  87, -90,

        90, 90, 88, 85, 82, 78, 73, 67,
        90, 82, 67, 46, 22, -4, -31, -54,
        88, 67, 31, -13, -54, -82, -90, -78,
        85, 46, -13, -67, -90, -73, -22, 38,
        82, 22, -54, -90, -61, 13, 78, 85,
        78, -4, -82, -73, 13, 85, 67, -22,
        73, -31, -90, -22, 78, 67, -38, -90,
        67, -54, -78, 38, 85, -22, -90, 4,
        61, -73, -46, 82, 31, -88, -13, 90,
        54, -85, -4, 88, -46, -61, 82, 13,
        46, -90, 38, 54, -90, 31, 61, -88,
        38, -88, 73, -4, -67, 90, -46, -31,
        31, -78, 90, -61, 4, 54, -88, 82,
        22, -61, 85, -90, 73, -38, -4, 46,
        13, -38, 61, -78, 88, -90, 85, -73,
        4, -13, 22, -31, 38, -46, 54, -61,

        61, 54, 46, 38, 31, 22, 13, 4,
        -73, -85, -90, -88, -78, -61, -38, -13,
        -46, -4, 38, 73, 90, 85, 61, 22,
        82, 88, 54, -4, -61, -90, -78, -31,
        31, -46, -90, -67, 4, 73, 88, 38,
        -88, -61, 31, 90, 54, -38, -90, -46,
        -13, 82, 61, -46, -88, -4, 85, 54,
        90, 13, -88, -31, 82, 46, -73, -61,
        -4, -90, 22, 85, -38, -78, 54, 67,
        -90, 38, 67, -78, -22, 90, -31, -73,
        22, 67, -85, 13, 73, -82, 4, 78,
        85, -78, 13, 61, -90, 54, 22, -82,
        -38, -22, 73, -90, 67, -13, -46, 85,
        -78, 90, -82, 54, -13, -31, 67, -88,
        54, -31, 4, 22, -46, 67, -82, 90,
        67, -73, 78, -82, 85, -88, 90, -90
    };

    for (int k=0; k<48; k++) {
        d[k] = _mm_load_si128((__m128i*)(DCT_II_32_4_sse + 8*k));
    }
}

static inline void
load_src_32x4_red4(const int16_t *src, ptrdiff_t src_stride, __m128i *x)
{
    __m128i *evn_xl = x;
    __m128i *odd_xl = x + 74;
    int k;

    x[0] = _mm_unpacklo_epi64(
            _mm_loadl_epi64((__m128i*)(src + 0 * src_stride)),
            _mm_setzero_si128()
            );
    #if 0
    x[1] = _mm_setzero_si128();
    x[2] = _mm_setzero_si128();
    x[3] = _mm_setzero_si128();
    x[4] = _mm_setzero_si128();
    x[5] = _mm_setzero_si128();
    x[6] = _mm_setzero_si128();
    x[7] = _mm_setzero_si128();
    x[8] = _mm_setzero_si128();
    x[9] = _mm_setzero_si128();
    #endif

    for (k = 0; k < 8; ++k){
        const int16_t *src_col = src + k;
        evn_xl[10] = _mm_set1_epi16(src_col[ 2 * src_stride]);
        #if 0
        evn_xl[11] = _mm_setzero_si128();
        evn_xl[12] = _mm_setzero_si128();
        evn_xl[13] = _mm_setzero_si128();
        evn_xl[14] = _mm_setzero_si128();
        evn_xl[15] = _mm_setzero_si128();
        evn_xl[16] = _mm_setzero_si128();
        evn_xl[17] = _mm_setzero_si128();
        #endif
        evn_xl += 8;
    }

    for (k = 0; k < 8; ++k){
        const int16_t *src_col = src + k;
        odd_xl[ 0] = _mm_set1_epi16(src_col[ 1 * src_stride]);
        odd_xl[ 1] = _mm_set1_epi16(src_col[ 3 * src_stride]);
        #if 0
        odd_xl[ 2] = _mm_setzero_si128();
        odd_xl[ 3] = _mm_setzero_si128();
        odd_xl[ 4] = _mm_setzero_si128();
        odd_xl[ 5] = _mm_setzero_si128();
        odd_xl[ 6] = _mm_setzero_si128();
        odd_xl[ 7] = _mm_setzero_si128();
        odd_xl[ 8] = _mm_setzero_si128();
        odd_xl[ 9] = _mm_setzero_si128();
        odd_xl[10] = _mm_setzero_si128();
        odd_xl[11] = _mm_setzero_si128();
        odd_xl[12] = _mm_setzero_si128();
        odd_xl[13] = _mm_setzero_si128();
        odd_xl[14] = _mm_setzero_si128();
        odd_xl[15] = _mm_setzero_si128();
        #endif
        odd_xl += 16;
    }
}

static inline void
load_src_32x4_red8(const int16_t *src, ptrdiff_t src_stride, __m128i *x)
{
    __m128i *evn_xl = x;
    __m128i *odd_xl = x + 74;
    int k;

    x[0] = _mm_unpacklo_epi64(
            _mm_loadl_epi64((__m128i*)(src + 0 * src_stride)),
            _mm_setzero_si128()
            );
    x[1] = _mm_setzero_si128();
    x[2] = _mm_unpacklo_epi64(
            _mm_set1_epi16(src[4 * src_stride + 0]),
            _mm_set1_epi16(src[4 * src_stride + 1])
            );
    x[3] = _mm_setzero_si128();
    x[4] = _mm_setzero_si128();
    x[5] = _mm_setzero_si128();
    x[6] = _mm_unpacklo_epi64(
            _mm_set1_epi16(src[4 * src_stride + 2]),
            _mm_set1_epi16(src[4 * src_stride + 3])
            );
    x[7] = _mm_setzero_si128();
    x[8] = _mm_setzero_si128();
    x[9] = _mm_setzero_si128();

    for (k = 0; k < 8; ++k){
        const int16_t *src_col = src + k;
        evn_xl[10] = _mm_set1_epi16(src_col[ 2 * src_stride]);
        evn_xl[11] = _mm_set1_epi16(src_col[ 6 * src_stride]);
        #if 0
        evn_xl[12] = _mm_setzero_si128();
        evn_xl[13] = _mm_setzero_si128();
        evn_xl[14] = _mm_setzero_si128();
        evn_xl[15] = _mm_setzero_si128();
        evn_xl[16] = _mm_setzero_si128();
        evn_xl[17] = _mm_setzero_si128();
        #endif
        evn_xl += 8;
    }

    for (k = 0; k < 8; ++k){
        const int16_t *src_col = src + k;
        odd_xl[ 0] = _mm_set1_epi16(src_col[ 1 * src_stride]);
        odd_xl[ 1] = _mm_set1_epi16(src_col[ 3 * src_stride]);
        odd_xl[ 2] = _mm_set1_epi16(src_col[ 5 * src_stride]);
        odd_xl[ 3] = _mm_set1_epi16(src_col[ 7 * src_stride]);
        #if 0
        odd_xl[ 4] = _mm_setzero_si128();
        odd_xl[ 5] = _mm_setzero_si128();
        odd_xl[ 6] = _mm_setzero_si128();
        odd_xl[ 7] = _mm_setzero_si128();
        odd_xl[ 8] = _mm_setzero_si128();
        odd_xl[ 9] = _mm_setzero_si128();
        odd_xl[10] = _mm_setzero_si128();
        odd_xl[11] = _mm_setzero_si128();
        odd_xl[12] = _mm_setzero_si128();
        odd_xl[13] = _mm_setzero_si128();
        odd_xl[14] = _mm_setzero_si128();
        odd_xl[15] = _mm_setzero_si128();
        #endif
        odd_xl += 16;
    }
}

static inline void
load_src_32x4_red16(const int16_t *src, ptrdiff_t src_stride, __m128i *x)
{
    __m128i *evn_xl = x;
    __m128i *odd_xl = x + 74;
    int k;

    x[0] = _mm_unpacklo_epi64(
            _mm_loadl_epi64((__m128i*)(src + 0 * src_stride)),
            _mm_loadl_epi64((__m128i*)(src + 8 * src_stride))
            );
    x[1] = _mm_unpacklo_epi64(
            _mm_loadl_epi64((__m128i*)(src + 16 * src_stride)),
            _mm_setzero_si128()
            );
    x[2] = _mm_unpacklo_epi64(
            _mm_set1_epi16(src[4 * src_stride + 0]),
            _mm_set1_epi16(src[4 * src_stride + 1])
            );
    x[3] = _mm_unpacklo_epi64(
            _mm_set1_epi16(src[12 * src_stride + 0]),
            _mm_set1_epi16(src[12 * src_stride + 1])
            );
    x[4] = _mm_setzero_si128();
    x[5] = _mm_setzero_si128();
    x[6] = _mm_unpacklo_epi64(
            _mm_set1_epi16(src[4 * src_stride + 2]),
            _mm_set1_epi16(src[4 * src_stride + 3])
            );
    x[7] = _mm_unpacklo_epi64(
            _mm_set1_epi16(src[12 * src_stride + 2]),
            _mm_set1_epi16(src[12 * src_stride + 3])
            );
    x[8] = _mm_setzero_si128();
    x[9] = _mm_setzero_si128();

    for (k = 0; k < 8; ++k){
        const int16_t *src_col = src + k;
        evn_xl[10] = _mm_set1_epi16(src_col[ 2 * src_stride]);
        evn_xl[11] = _mm_set1_epi16(src_col[ 6 * src_stride]);
        evn_xl[12] = _mm_set1_epi16(src_col[10 * src_stride]);
        evn_xl[13] = _mm_set1_epi16(src_col[14 * src_stride]);
        #if 0
        evn_xl[14] = _mm_setzero_si128();
        evn_xl[15] = _mm_setzero_si128();
        evn_xl[16] = _mm_setzero_si128();
        evn_xl[17] = _mm_setzero_si128();
        #endif
        evn_xl += 8;
    }

    for (k = 0; k < 8; ++k){
        const int16_t *src_col = src + k;
        odd_xl[ 0] = _mm_set1_epi16(src_col[ 1 * src_stride]);
        odd_xl[ 1] = _mm_set1_epi16(src_col[ 3 * src_stride]);
        odd_xl[ 2] = _mm_set1_epi16(src_col[ 5 * src_stride]);
        odd_xl[ 3] = _mm_set1_epi16(src_col[ 7 * src_stride]);
        odd_xl[ 4] = _mm_set1_epi16(src_col[ 9 * src_stride]);
        odd_xl[ 5] = _mm_set1_epi16(src_col[11 * src_stride]);
        odd_xl[ 6] = _mm_set1_epi16(src_col[13 * src_stride]);
        odd_xl[ 7] = _mm_set1_epi16(src_col[15 * src_stride]);
        #if 0
        odd_xl[ 8] = _mm_setzero_si128();
        odd_xl[ 9] = _mm_setzero_si128();
        odd_xl[10] = _mm_setzero_si128();
        odd_xl[11] = _mm_setzero_si128();
        odd_xl[12] = _mm_setzero_si128();
        odd_xl[13] = _mm_setzero_si128();
        odd_xl[14] = _mm_setzero_si128();
        odd_xl[15] = _mm_setzero_si128();
        #endif
        odd_xl += 16;
    }
}

static inline void
load_src_32x4(const int16_t *src, ptrdiff_t src_stride, __m128i *x)
{
    __m128i *evn_xl = x;
    __m128i *odd_xl = x + 74;
    int k;

    x[0] = _mm_unpacklo_epi64(
            _mm_loadl_epi64((__m128i*)(src + 0 * src_stride)),
            _mm_loadl_epi64((__m128i*)(src + 8 * src_stride))
            );
    x[1] = _mm_unpacklo_epi64(
            _mm_loadl_epi64((__m128i*)(src + 16 * src_stride)),
            _mm_loadl_epi64((__m128i*)(src + 24 * src_stride))
            );
    x[2] = _mm_unpacklo_epi64(
            _mm_set1_epi16(src[4 * src_stride + 0]),
            _mm_set1_epi16(src[4 * src_stride + 1])
            );
    x[3] = _mm_unpacklo_epi64(
            _mm_set1_epi16(src[12 * src_stride + 0]),
            _mm_set1_epi16(src[12 * src_stride + 1])
            );
    x[4] = _mm_unpacklo_epi64(
            _mm_set1_epi16(src[20 * src_stride + 0]),
            _mm_set1_epi16(src[20 * src_stride + 1])
            );
    x[5] = _mm_unpacklo_epi64(
            _mm_set1_epi16(src[28 * src_stride + 0]),
            _mm_set1_epi16(src[28 * src_stride + 1])
            );
    x[6] = _mm_unpacklo_epi64(
            _mm_set1_epi16(src[4 * src_stride + 2]),
            _mm_set1_epi16(src[4 * src_stride + 3])
            );
    x[7] = _mm_unpacklo_epi64(
            _mm_set1_epi16(src[12 * src_stride + 2]),
            _mm_set1_epi16(src[12 * src_stride + 3])
            );
    x[8] = _mm_unpacklo_epi64(
            _mm_set1_epi16(src[20 * src_stride + 2]),
            _mm_set1_epi16(src[20 * src_stride + 3])
            );
    x[9] = _mm_unpacklo_epi64(
            _mm_set1_epi16(src[28 * src_stride + 2]),
            _mm_set1_epi16(src[28 * src_stride + 3])
            );

    for (k = 0; k < 8; ++k){
        const int16_t *src_col = src + k;
        evn_xl[10] = _mm_set1_epi16(src_col[ 2 * src_stride]);
        evn_xl[11] = _mm_set1_epi16(src_col[ 6 * src_stride]);
        evn_xl[12] = _mm_set1_epi16(src_col[10 * src_stride]);
        evn_xl[13] = _mm_set1_epi16(src_col[14 * src_stride]);
        evn_xl[14] = _mm_set1_epi16(src_col[18 * src_stride]);
        evn_xl[15] = _mm_set1_epi16(src_col[22 * src_stride]);
        evn_xl[16] = _mm_set1_epi16(src_col[26 * src_stride]);
        evn_xl[17] = _mm_set1_epi16(src_col[30 * src_stride]);
        evn_xl += 8;
    }

    for (k = 0; k < 8; ++k){
        const int16_t *src_col = src + k;
        odd_xl[ 0] = _mm_set1_epi16(src_col[ 1 * src_stride]);
        odd_xl[ 1] = _mm_set1_epi16(src_col[ 3 * src_stride]);
        odd_xl[ 2] = _mm_set1_epi16(src_col[ 5 * src_stride]);
        odd_xl[ 3] = _mm_set1_epi16(src_col[ 7 * src_stride]);
        odd_xl[ 4] = _mm_set1_epi16(src_col[ 9 * src_stride]);
        odd_xl[ 5] = _mm_set1_epi16(src_col[11 * src_stride]);
        odd_xl[ 6] = _mm_set1_epi16(src_col[13 * src_stride]);
        odd_xl[ 7] = _mm_set1_epi16(src_col[15 * src_stride]);
        odd_xl[ 8] = _mm_set1_epi16(src_col[17 * src_stride]);
        odd_xl[ 9] = _mm_set1_epi16(src_col[19 * src_stride]);
        odd_xl[10] = _mm_set1_epi16(src_col[21 * src_stride]);
        odd_xl[11] = _mm_set1_epi16(src_col[23 * src_stride]);
        odd_xl[12] = _mm_set1_epi16(src_col[25 * src_stride]);
        odd_xl[13] = _mm_set1_epi16(src_col[27 * src_stride]);
        odd_xl[14] = _mm_set1_epi16(src_col[29 * src_stride]);
        odd_xl[15] = _mm_set1_epi16(src_col[31 * src_stride]);
        odd_xl += 16;
    }
}

static void
dct2_32_4lines_red4(const int16_t *src, int16_t *dst, ptrdiff_t src_stride,
               int num_lines, int line_brk, int shift)
{
    __m128i x[202], d[48], r[32];

    load_src_32x4_red4(src, src_stride, x);

    load_dct_32x4(d);

    dct2_32x4_red4(x, d, r);

    __m128i add = _mm_set1_epi32(1 << (shift - 1));
    for (int i = 0; i < 32; i += 2) {
        __m128i o;
        r[i + 0] = _mm_add_epi32(r[i + 0], add);
        r[i + 1] = _mm_add_epi32(r[i + 1], add);

        r[i + 0] = _mm_srai_epi32(r[i + 0], shift);
        r[i + 1] = _mm_srai_epi32(r[i + 1], shift);

        o = _mm_packs_epi32(r[i + 0], r[i + 1]);

        _mm_store_si128((__m128i *) (dst + i / 2 * 8), o);
    }
}

static void
dct2_32_4lines_red8(const int16_t *src, int16_t *dst, ptrdiff_t src_stride,
               int num_lines, int line_brk, int shift)
{
    __m128i x[202], d[48], r[32];

    load_src_32x4_red8(src, src_stride, x);

    load_dct_32x4(d);

    dct2_32x4_red8(x, d, r);

    __m128i add = _mm_set1_epi32(1 << (shift - 1));
    for (int i = 0; i < 32; i += 2) {
        __m128i o;
        r[i + 0] = _mm_add_epi32(r[i + 0], add);
        r[i + 1] = _mm_add_epi32(r[i + 1], add);

        r[i + 0] = _mm_srai_epi32(r[i + 0], shift);
        r[i + 1] = _mm_srai_epi32(r[i + 1], shift);

        o = _mm_packs_epi32(r[i + 0], r[i + 1]);

        _mm_store_si128((__m128i *) (dst + i / 2 * 8), o);
    }
}

static void
dct2_32_4lines_red16(const int16_t *src, int16_t *dst, ptrdiff_t src_stride,
               int num_lines, int line_brk, int shift)
{
    __m128i x[202], d[48], r[32];

    load_src_32x4_red16(src, src_stride, x);

    load_dct_32x4(d);

    dct2_32x4_red16(x, d, r);

    __m128i add = _mm_set1_epi32(1 << (shift - 1));
    for (int i = 0; i < 32; i += 2) {
        __m128i o;
        r[i + 0] = _mm_add_epi32(r[i + 0], add);
        r[i + 1] = _mm_add_epi32(r[i + 1], add);

        r[i + 0] = _mm_srai_epi32(r[i + 0], shift);
        r[i + 1] = _mm_srai_epi32(r[i + 1], shift);

        o = _mm_packs_epi32(r[i + 0], r[i + 1]);

        _mm_store_si128((__m128i *) (dst + i / 2 * 8), o);
    }
}

static void
dct2_32_4lines(const int16_t *src, int16_t *dst, ptrdiff_t src_stride,
               int num_lines, int line_brk, int shift)
{
    __m128i x[202], d[48], r[32];

    load_src_32x4(src, src_stride, x);

    load_dct_32x4(d);

    dct2_32x4(x, d, r);

    __m128i add = _mm_set1_epi32(1 << (shift - 1));
    for (int i = 0; i < 32; i += 2) {
        __m128i o;
        r[i + 0] = _mm_add_epi32(r[i + 0], add);
        r[i + 1] = _mm_add_epi32(r[i + 1], add);

        r[i + 0] = _mm_srai_epi32(r[i + 0], shift);
        r[i + 1] = _mm_srai_epi32(r[i + 1], shift);

        o = _mm_packs_epi32(r[i + 0], r[i + 1]);

        _mm_store_si128((__m128i *) (dst + i / 2 * 8), o);
    }
}

void
vvc_inverse_dct_ii_32_sse(const int16_t *src, int16_t *dst, ptrdiff_t src_stride,
                          int num_lines, int line_brk, int shift)
{

    /*FIXME it would be more efficient to get per reduced lines speialised
      functions instead off checking */
    if (num_lines/8){
        if (line_brk > 16)
            dct2_32_8lines(src, dst, src_stride, num_lines, line_brk, shift);
        else if (line_brk > 8)
            dct2_32_8lines_red16(src, dst, src_stride, num_lines, line_brk, shift);
        else
            dct2_32_8lines_red8(src, dst, src_stride, num_lines, line_brk, shift);
        src += (num_lines & 0xF8);
        dst += (num_lines >> 3) << 8;
    }

    if (!(num_lines & 0x7)) return;

    if (num_lines & 0x4){
        if (line_brk > 16)
            dct2_32_4lines(src, dst, src_stride, num_lines, line_brk, shift);
        else if (line_brk > 8)
            dct2_32_4lines_red16(src, dst, src_stride, num_lines, line_brk, shift);
        else if (line_brk > 4)
            dct2_32_4lines_red8(src, dst, src_stride, num_lines, line_brk, shift);
        else
            dct2_32_4lines_red4(src, dst, src_stride, num_lines, line_brk, shift);

        src += 4;
        dst += 128;
    }

    if (num_lines & 0x2){
        __m128i x[56], d[46], r[16];

        x[0] = _mm_unpacklo_epi32(
                _mm_loadl_epi64((__m128i*)(src + 0 * src_stride)),
                _mm_loadl_epi64((__m128i*)(src + 8 * src_stride))
                );
        x[0] = _mm_unpacklo_epi16(x[0], x[0]);
        x[1] = _mm_unpacklo_epi32(
                _mm_loadl_epi64((__m128i*)(src + 16 * src_stride)),
                _mm_loadl_epi64((__m128i*)(src + 24 * src_stride))
                );
        x[1] = _mm_unpacklo_epi16(x[1], x[1]);

        x[2] = _mm_unpacklo_epi64(
                _mm_set1_epi16(src[ 4 * src_stride    ]),
                _mm_set1_epi16(src[ 4 * src_stride + 1])
                );
        x[3] = _mm_unpacklo_epi64(
                _mm_set1_epi16(src[12 * src_stride    ]),
                _mm_set1_epi16(src[12 * src_stride + 1])
                );
        x[4] = _mm_unpacklo_epi64(
                _mm_set1_epi16(src[20 * src_stride    ]),
                _mm_set1_epi16(src[20 * src_stride + 1])
                );
        x[5] = _mm_unpacklo_epi64(
                _mm_set1_epi16(src[28 * src_stride    ]),
                _mm_set1_epi16(src[28 * src_stride + 1])
                );

        for (int k=0; k<2; k++){
            x[6+8*k+0] = _mm_set1_epi16(src[ 2 * src_stride + k]);
            x[6+8*k+1] = _mm_set1_epi16(src[ 6 * src_stride + k]);
            x[6+8*k+2] = _mm_set1_epi16(src[10 * src_stride + k]);
            x[6+8*k+3] = _mm_set1_epi16(src[14 * src_stride + k]);
            x[6+8*k+4] = _mm_set1_epi16(src[18 * src_stride + k]);
            x[6+8*k+5] = _mm_set1_epi16(src[22 * src_stride + k]);
            x[6+8*k+6] = _mm_set1_epi16(src[26 * src_stride + k]);
            x[6+8*k+7] = _mm_set1_epi16(src[30 * src_stride + k]);
        }

        for (int k=0; k<2; k++) {
            for (int l = 0; l < 16; l++) {
                x[24+16*k+l] = _mm_set1_epi16(src[(2*l+1) * src_stride + k]);
            }
        }

        static const int16_t DCT_II_32_2_sse[8 * 46] = {
            64,  64,  64,  64,  83,  36,  83,  36,
            64, -64,  64, -64,  36, -83,  36, -83,

            89,  75,  18,  50,  89,  75,  18,  50,
            75, -18, -50, -89,  75, -18, -50, -89,
            50, -89,  75,  18,  50, -89,  75,  18,
            18, -50, -89,  75,  18, -50, -89,  75,

            90,  87,  80,  70,  57,  43,  25,   9,
            87,  57,   9, -43, -80, -90, -70, -25,
            80,   9, -70, -87, -25,  57,  90,  43,
            70, -43, -87,   9,  90,  25, -80, -57,
            57, -80, -25,  90,  -9, -87,  43,  70,
            43, -90,  57,  25, -87,  70,   9, -80,
            25, -70,  90, -80,  43,   9, -57,  87,
            9, -25,  43, -57,  70, -80,  87, -90,

            90, 90, 88, 85, 82, 78, 73, 67,
            90, 82, 67, 46, 22, -4, -31, -54,
            88, 67, 31, -13, -54, -82, -90, -78,
            85, 46, -13, -67, -90, -73, -22, 38,
            82, 22, -54, -90, -61, 13, 78, 85,
            78, -4, -82, -73, 13, 85, 67, -22,
            73, -31, -90, -22, 78, 67, -38, -90,
            67, -54, -78, 38, 85, -22, -90, 4,
            61, -73, -46, 82, 31, -88, -13, 90,
            54, -85, -4, 88, -46, -61, 82, 13,
            46, -90, 38, 54, -90, 31, 61, -88,
            38, -88, 73, -4, -67, 90, -46, -31,
            31, -78, 90, -61, 4, 54, -88, 82,
            22, -61, 85, -90, 73, -38, -4, 46,
            13, -38, 61, -78, 88, -90, 85, -73,
            4, -13, 22, -31, 38, -46, 54, -61,

            61, 54, 46, 38, 31, 22, 13, 4,
            -73, -85, -90, -88, -78, -61, -38, -13,
            -46, -4, 38, 73, 90, 85, 61, 22,
            82, 88, 54, -4, -61, -90, -78, -31,
            31, -46, -90, -67, 4, 73, 88, 38,
            -88, -61, 31, 90, 54, -38, -90, -46,
            -13, 82, 61, -46, -88, -4, 85, 54,
            90, 13, -88, -31, 82, 46, -73, -61,
            -4, -90, 22, 85, -38, -78, 54, 67,
            -90, 38, 67, -78, -22, 90, -31, -73,
            22, 67, -85, 13, 73, -82, 4, 78,
            85, -78, 13, 61, -90, 54, 22, -82,
            -38, -22, 73, -90, 67, -13, -46, 85,
            -78, 90, -82, 54, -13, -31, 67, -88,
            54, -31, 4, 22, -46, 67, -82, 90,
            67, -73, 78, -82, 85, -88, 90, -90
        };

        for (int k=0; k<46; k++) {
            d[k] = _mm_load_si128((const __m128i*) (DCT_II_32_2_sse + 8*k));
        }

        dct2_32x2(x, d, r);

        __m128i add = _mm_set1_epi32(1 << (shift - 1));

        for (int i = 0; i < 16; i += 2) {
            __m128i o;

            r[i + 0] = _mm_add_epi32(r[i + 0], add);
            r[i + 1] = _mm_add_epi32(r[i + 1], add);

            r[i + 0] = _mm_srai_epi32(r[i + 0], shift);
            r[i + 1] = _mm_srai_epi32(r[i + 1], shift);

            o = _mm_packs_epi32(r[i + 0], r[i + 1]);

            _mm_store_si128((__m128i *) (dst + i / 2 * 8), o);
        }
    }

    if (num_lines & 0x1){
        vvc_inverse_dct_ii_32(src, dst, src_stride, num_lines & 0x1, line_brk, shift);
    }
}

void
vvc_inverse_dct_ii_64_sse_8lines_red16(const int16_t* src, int16_t* dst, ptrdiff_t src_stride,
                      int num_lines, int num_columns, int shift)
                      {

                      }
void
vvc_inverse_dct_ii_64_sse_8lines(const int16_t* src, int16_t* dst, ptrdiff_t src_stride,
                      int num_lines, int num_columns, int shift)
{
        __m128i add = _mm_set1_epi32(1 << (shift - 1));

        int j, k;

        for (j = 0; j < num_lines>>3; j++) {
                /* Utilizing symmetry properties to the maximum to minimize
                the
                 * number of multiplications */
                 __m128i x[16], d[16], di[8];
                 __m128i m[16], a[16], b[8], c[4], r[128], r2[128], r3[128];
                 //EEEEO
                 x[0 ]=_mm_load_si128((__m128i*)(src +  16 * src_stride));
                 d[0] = _mm_set1_epi16(DCT_II_64_EEEEO[0]);
                 d[1] = _mm_set1_epi16(DCT_II_64_EEEEO[1]);

                 m[0] = _mm_mullo_epi16(x[0], d[0]);
                 m[1] = _mm_mulhi_epi16(x[0], d[0]);
                 m[2] = _mm_mullo_epi16(x[0], d[1]);
                 m[3] = _mm_mulhi_epi16(x[0], d[1]);

                 b[0] = _mm_unpacklo_epi16(m[0], m[1]);
                 b[1] = _mm_unpackhi_epi16(m[0], m[1]);

                 b[2] = _mm_unpacklo_epi16(m[2], m[3]);
                 b[3] = _mm_unpackhi_epi16(m[2], m[3]);

                 //EEEEE
                 x[0]=_mm_load_si128((__m128i*)(src));
                 d[0] = _mm_set1_epi16(DCT_II_64_EEEEE[0]);
                 d[1] = _mm_set1_epi16(DCT_II_64_EEEEE[1]);

                 m[0] = _mm_mullo_epi16(x[0], d[0]);
                 m[1] = _mm_mulhi_epi16(x[0], d[0]);
                 m[2] = _mm_mullo_epi16(x[0], d[1]);
                 m[3] = _mm_mulhi_epi16(x[0], d[1]);

                 c[0] = _mm_unpacklo_epi16(m[0], m[1]);
                 c[1] = _mm_unpackhi_epi16(m[0], m[1]);

                 c[2] = _mm_unpacklo_epi16(m[2], m[3]);
                 c[3] = _mm_unpackhi_epi16(m[2], m[3]);

                 /* Combining even and odd terms at each hierarchy levels to
                  * calculate the final spatial domain vector */
                 //EEEE
                 for (k = 0; k < 2; k++) {
                   a[k*2] = _mm_add_epi32(b[k*2], c[k*2]);
                   a[k*2+1] = _mm_add_epi32(b[k*2+1], c[k*2+1]);

                   a[k*2+4] = _mm_sub_epi32(c[2-k*2], b[2-k*2]);
                   a[k*2+5] = _mm_sub_epi32(c[3-k*2], b[3-k*2]);
                 }

                 //EEEO
                 x[0 ]=_mm_load_si128((__m128i*)(src +  8 * src_stride));
                 x[1 ]=_mm_load_si128((__m128i*)(src +  24 * src_stride));
                 for (k = 0; k < 4; k++) {
                   d[ 0] = _mm_set1_epi16(DCT_II_64_EEEOT[k * 2 + 0]);
                   d[ 1] = _mm_set1_epi16(DCT_II_64_EEEOT[k * 2 + 1]);

                   m[ 0] = _mm_unpacklo_epi16(x[0],  x[1]);
                   m[ 1] = _mm_unpackhi_epi16(x[0],  x[1]);

                   di[0] = _mm_unpacklo_epi16(d[0],  d[1]);

                   r3[k*2] = _mm_madd_epi16(m[0], di[0]);
                   r3[k*2+1] = _mm_madd_epi16(m[1], di[0]);
                 }
                 //EEE
                 /* Combining even and odd terms at each hierarchy levels to
                  * calculate the final spatial domain vector */
                  for (k = 0; k < 4; k++) {
                    r[k*2] = _mm_add_epi32(r3[k*2], a[k*2]);
                    r[k*2+1] = _mm_add_epi32(r3[k*2+1], a[k*2+1]);

                    r[k*2+8] = _mm_sub_epi32(a[6-k*2], r3[6-k*2]);
                    r[k*2+9] = _mm_sub_epi32(a[7-k*2], r3[7-k*2]);
                  }

                 //EEO
                 x[0 ]=_mm_load_si128((__m128i*)(src +  4 * src_stride));
                 x[1 ]=_mm_load_si128((__m128i*)(src +  12 * src_stride));
                 x[2 ]=_mm_load_si128((__m128i*)(src +  20 * src_stride));
                 x[3 ]=_mm_load_si128((__m128i*)(src +  28 * src_stride));
                 for (k = 0; k < 8; k++) {
                   d[ 0] = _mm_set1_epi16(DCT_II_64_EEOT[k * 4 + 0]);
                   d[ 1] = _mm_set1_epi16(DCT_II_64_EEOT[k * 4 + 1]);
                   d[ 2] = _mm_set1_epi16(DCT_II_64_EEOT[k * 4 + 2]);
                   d[ 3] = _mm_set1_epi16(DCT_II_64_EEOT[k * 4 + 3]);

                   m[ 0] = _mm_unpacklo_epi16(x[0],  x[1]);
                   m[ 1] = _mm_unpacklo_epi16(x[2],  x[3]);

                   m[ 2] = _mm_unpackhi_epi16(x[0],  x[1]);
                   m[ 3] = _mm_unpackhi_epi16(x[2],  x[3]);

                   di[0] = _mm_unpacklo_epi16(d[0],  d[1]);
                   di[1] = _mm_unpacklo_epi16(d[2],  d[3]);


                   a[0] = _mm_madd_epi16(m[0], di[0]);
                   a[1] = _mm_madd_epi16(m[1], di[1]);
                   a[2] = _mm_madd_epi16(m[2], di[0]);
                   a[3] = _mm_madd_epi16(m[3], di[1]);


                   r2[k*2] = _mm_add_epi32(a[0], a[1]);
                   r2[k*2+1] = _mm_add_epi32(a[2], a[3]);
                 }
                 //EE
                 /* Combining even and odd terms at each hierarchy levels to
                  * calculate the final spatial domain vector */
                 for (k = 0; k < 8; k++) {
                   r3[k*2] = _mm_add_epi32(r2[k*2], r[k*2]);
                   r3[k*2+1] = _mm_add_epi32(r2[k*2+1], r[k*2+1]);

                   r3[k*2+16] = _mm_sub_epi32(r[14-k*2], r2[14-k*2]);
                   r3[k*2+17] = _mm_sub_epi32(r[15-k*2], r2[15-k*2]);
                 }

                 //EO
                 x[0 ]=_mm_load_si128((__m128i*)(src +  2 * src_stride));
                 x[1 ]=_mm_load_si128((__m128i*)(src +  6 * src_stride));
                 x[2 ]=_mm_load_si128((__m128i*)(src +  10 * src_stride));
                 x[3 ]=_mm_load_si128((__m128i*)(src +  14 * src_stride));
                 x[4 ]=_mm_load_si128((__m128i*)(src +  18 * src_stride));
                 x[5 ]=_mm_load_si128((__m128i*)(src +  22 * src_stride));
                 x[6 ]=_mm_load_si128((__m128i*)(src +  26 * src_stride));
                 x[7 ]=_mm_load_si128((__m128i*)(src +  30 * src_stride));
                 for (k = 0; k < 16; k++) {
                   d[ 0] = _mm_set1_epi16(DCT_II_64_EOT[k * 8 + 0]);
                   d[ 1] = _mm_set1_epi16(DCT_II_64_EOT[k * 8 + 1]);
                   d[ 2] = _mm_set1_epi16(DCT_II_64_EOT[k * 8 + 2]);
                   d[ 3] = _mm_set1_epi16(DCT_II_64_EOT[k * 8 + 3]);
                   d[ 4] = _mm_set1_epi16(DCT_II_64_EOT[k * 8 + 4]);
                   d[ 5] = _mm_set1_epi16(DCT_II_64_EOT[k * 8 + 5]);
                   d[ 6] = _mm_set1_epi16(DCT_II_64_EOT[k * 8 + 6]);
                   d[ 7] = _mm_set1_epi16(DCT_II_64_EOT[k * 8 + 7]);

                   m[ 0] = _mm_unpacklo_epi16(x[0],  x[1]);
                   m[ 1] = _mm_unpacklo_epi16(x[2],  x[3]);
                   m[ 2] = _mm_unpacklo_epi16(x[4],  x[5]);
                   m[ 3] = _mm_unpacklo_epi16(x[6],  x[7]);

                   m[ 4] = _mm_unpackhi_epi16(x[0],  x[1]);
                   m[ 5] = _mm_unpackhi_epi16(x[2],  x[3]);
                   m[ 6] = _mm_unpackhi_epi16(x[4],  x[5]);
                   m[ 7] = _mm_unpackhi_epi16(x[6],  x[7]);

                   di[0] = _mm_unpacklo_epi16(d[0],  d[1]);
                   di[1] = _mm_unpacklo_epi16(d[2],  d[3]);
                   di[2] = _mm_unpacklo_epi16(d[4],  d[5]);
                   di[3] = _mm_unpacklo_epi16(d[6],  d[7]);

                   a[0] = _mm_madd_epi16(m[0], di[0]);
                   a[1] = _mm_madd_epi16(m[1], di[1]);
                   a[2] = _mm_madd_epi16(m[2], di[2]);
                   a[3] = _mm_madd_epi16(m[3], di[3]);

                   a[4] = _mm_madd_epi16(m[4], di[0]);
                   a[5] = _mm_madd_epi16(m[5], di[1]);
                   a[6] = _mm_madd_epi16(m[6], di[2]);
                   a[7] = _mm_madd_epi16(m[7], di[3]);

                   b[0] = _mm_add_epi32(a[0], a[1]);
                   b[1] = _mm_add_epi32(a[2], a[3]);

                   b[2] = _mm_add_epi32(a[4], a[5]);
                   b[3] = _mm_add_epi32(a[6], a[7]);


                   r[k*2] = _mm_add_epi32(b[0], b[1]);

                   r[k*2+1] = _mm_add_epi32(b[2], b[3]);
                 }

                 //E
                 /* Combining even and odd terms at each hierarchy levels to
                  * calculate the final spatial domain vector */
                 for (k = 0; k < 16; k++) {
                   r2[k*2] = _mm_add_epi32(r[k*2], r3[k*2]);
                   r2[k*2+1] = _mm_add_epi32(r[k*2+1], r3[k*2+1]);

                   r2[k*2+32] = _mm_sub_epi32(r3[30-k*2], r[30-k*2]);
                   r2[k*2+33] = _mm_sub_epi32(r3[31-k*2], r[31-k*2]);
                 }
                 //O
                 x[0 ]=_mm_load_si128((__m128i*)(src +  src_stride));
                 x[1 ]=_mm_load_si128((__m128i*)(src +  3  * src_stride));
                 x[2 ]=_mm_load_si128((__m128i*)(src +  5  * src_stride));
                 x[3 ]=_mm_load_si128((__m128i*)(src +  7  * src_stride));
                 x[4 ]=_mm_load_si128((__m128i*)(src +  9  * src_stride));
                 x[5 ]=_mm_load_si128((__m128i*)(src +  11 * src_stride));
                 x[6 ]=_mm_load_si128((__m128i*)(src +  13 * src_stride));
                 x[7 ]=_mm_load_si128((__m128i*)(src +  15 * src_stride));
                 x[8 ]=_mm_load_si128((__m128i*)(src +  17 * src_stride));
                 x[9 ]=_mm_load_si128((__m128i*)(src +  19 * src_stride));
                 x[10]=_mm_load_si128((__m128i*)(src +  21 * src_stride));
                 x[11]=_mm_load_si128((__m128i*)(src +  23 * src_stride));
                 x[12]=_mm_load_si128((__m128i*)(src +  25 * src_stride));
                 x[13]=_mm_load_si128((__m128i*)(src +  27 * src_stride));
                 x[14]=_mm_load_si128((__m128i*)(src +  29 * src_stride));
                 x[15]=_mm_load_si128((__m128i*)(src +  31 * src_stride));
                 for (k = 0; k < 32; k++) {
                   d[ 0] = _mm_set1_epi16(DCT_II_64_OT[k * 16 + 0 ]);
                   d[ 1] = _mm_set1_epi16(DCT_II_64_OT[k * 16 + 1 ]);
                   d[ 2] = _mm_set1_epi16(DCT_II_64_OT[k * 16 + 2 ]);
                   d[ 3] = _mm_set1_epi16(DCT_II_64_OT[k * 16 + 3 ]);
                   d[ 4] = _mm_set1_epi16(DCT_II_64_OT[k * 16 + 4 ]);
                   d[ 5] = _mm_set1_epi16(DCT_II_64_OT[k * 16 + 5 ]);
                   d[ 6] = _mm_set1_epi16(DCT_II_64_OT[k * 16 + 6 ]);
                   d[ 7] = _mm_set1_epi16(DCT_II_64_OT[k * 16 + 7 ]);
                   d[ 8] = _mm_set1_epi16(DCT_II_64_OT[k * 16 + 8 ]);
                   d[ 9] = _mm_set1_epi16(DCT_II_64_OT[k * 16 + 9 ]);
                   d[10] = _mm_set1_epi16(DCT_II_64_OT[k * 16 + 10]);
                   d[11] = _mm_set1_epi16(DCT_II_64_OT[k * 16 + 11]);
                   d[12] = _mm_set1_epi16(DCT_II_64_OT[k * 16 + 12]);
                   d[13] = _mm_set1_epi16(DCT_II_64_OT[k * 16 + 13]);
                   d[14] = _mm_set1_epi16(DCT_II_64_OT[k * 16 + 14]);
                   d[15] = _mm_set1_epi16(DCT_II_64_OT[k * 16 + 15]);

                   m[ 0] = _mm_unpacklo_epi16(x[0],  x[1]);
                   m[ 1] = _mm_unpacklo_epi16(x[2],  x[3]);
                   m[ 2] = _mm_unpacklo_epi16(x[4],  x[5]);
                   m[ 3] = _mm_unpacklo_epi16(x[6],  x[7]);
                   m[ 4] = _mm_unpacklo_epi16(x[8],  x[9]);
                   m[ 5] = _mm_unpacklo_epi16(x[10], x[11]);
                   m[ 6] = _mm_unpacklo_epi16(x[12], x[13]);
                   m[ 7] = _mm_unpacklo_epi16(x[14], x[15]);

                   m[ 8] = _mm_unpackhi_epi16(x[0],  x[1]);
                   m[ 9] = _mm_unpackhi_epi16(x[2],  x[3]);
                   m[10] = _mm_unpackhi_epi16(x[4],  x[5]);
                   m[11] = _mm_unpackhi_epi16(x[6],  x[7]);
                   m[12] = _mm_unpackhi_epi16(x[8],  x[9]);
                   m[13] = _mm_unpackhi_epi16(x[10], x[11]);
                   m[14] = _mm_unpackhi_epi16(x[12], x[13]);
                   m[15] = _mm_unpackhi_epi16(x[14], x[15]);

                   di[0] = _mm_unpacklo_epi16(d[0],  d[1]);
                   di[1] = _mm_unpacklo_epi16(d[2],  d[3]);
                   di[2] = _mm_unpacklo_epi16(d[4],  d[5]);
                   di[3] = _mm_unpacklo_epi16(d[6],  d[7]);
                   di[4] = _mm_unpacklo_epi16(d[8],  d[9]);
                   di[5] = _mm_unpacklo_epi16(d[10], d[11]);
                   di[6] = _mm_unpacklo_epi16(d[12], d[13]);
                   di[7] = _mm_unpacklo_epi16(d[14], d[15]);

                   a[0] = _mm_madd_epi16(m[0], di[0]);
                   a[1] = _mm_madd_epi16(m[1], di[1]);
                   a[2] = _mm_madd_epi16(m[2], di[2]);
                   a[3] = _mm_madd_epi16(m[3], di[3]);
                   a[4] = _mm_madd_epi16(m[4], di[4]);
                   a[5] = _mm_madd_epi16(m[5], di[5]);
                   a[6] = _mm_madd_epi16(m[6], di[6]);
                   a[7] = _mm_madd_epi16(m[7], di[7]);

                   a[ 8] = _mm_madd_epi16(m[ 8], di[0]);
                   a[ 9] = _mm_madd_epi16(m[ 9], di[1]);
                   a[10] = _mm_madd_epi16(m[10], di[2]);
                   a[11] = _mm_madd_epi16(m[11], di[3]);
                   a[12] = _mm_madd_epi16(m[12], di[4]);
                   a[13] = _mm_madd_epi16(m[13], di[5]);
                   a[14] = _mm_madd_epi16(m[14], di[6]);
                   a[15] = _mm_madd_epi16(m[15], di[7]);

                   b[0] = _mm_add_epi32(a[0], a[1]);
                   b[1] = _mm_add_epi32(a[2], a[3]);
                   b[2] = _mm_add_epi32(a[4], a[5]);
                   b[3] = _mm_add_epi32(a[6], a[7]);

                   b[4] = _mm_add_epi32(a[ 8], a[9]);
                   b[5] = _mm_add_epi32(a[10], a[11]);
                   b[6] = _mm_add_epi32(a[12], a[13]);
                   b[7] = _mm_add_epi32(a[14], a[15]);

                   c[0] = _mm_add_epi32(b[0], b[1]);
                   c[1] = _mm_add_epi32(b[2], b[3]);

                   c[2] = _mm_add_epi32(b[4], b[5]);
                   c[3] = _mm_add_epi32(b[6], b[7]);

                   r3[k*2] = _mm_add_epi32(c[0], c[1]);
                   r3[k*2+1] = _mm_add_epi32(c[2], c[3]);
                 }

                 //Result
                 for (k = 0; k < 32; k++) {
                   r[k*2] =   _mm_add_epi32(r3[k*2], r2[k*2]);
                   r[k*2] =   _mm_add_epi32(r[k*2],add);
                   r[k*2] =   _mm_srai_epi32(r[k*2],shift);
                   r[k*2+1] = _mm_add_epi32(r3[k*2+1], r2[k*2+1]);
                   r[k*2+1] = _mm_add_epi32(r[k*2+1],add);
                   r[k*2+1] = _mm_srai_epi32(r[k*2+1],shift);
                   r[k*2] =_mm_packs_epi32(r[k*2], r[k*2+1]);

                   r[k*2+64] = _mm_sub_epi32(r2[62-k*2], r3[62-k*2]);
                   r[k*2+64] = _mm_add_epi32(r[k*2+64],add);
                   r[k*2+64] = _mm_srai_epi32(r[k*2+64],shift);
                   r[k*2+65] = _mm_sub_epi32(r2[63-k*2], r3[63-k*2]);
                   r[k*2+65] = _mm_add_epi32(r[k*2+65],add);
                   r[k*2+65] = _mm_srai_epi32(r[k*2+65],shift);
                   r[k*2+64] =_mm_packs_epi32(r[k*2+64], r[k*2+65]);
                 }

                for (k = 0; k < 8; k++) {
                      __m128i tmp[8]; __m128i tmp2[8];
                      tmp[0] = _mm_unpacklo_epi16(r[k*16 + 0], r[k*16 + 2]);
                      tmp[1] = _mm_unpackhi_epi16(r[k*16 + 0], r[k*16 + 2]);
                      tmp[2] = _mm_unpacklo_epi16(r[k*16 + 4], r[k*16 + 6]);
                      tmp[3] = _mm_unpackhi_epi16(r[k*16 + 4], r[k*16 + 6]);
                      tmp[4] = _mm_unpacklo_epi16(r[k*16 + 8], r[k*16 + 10]);
                      tmp[5] = _mm_unpackhi_epi16(r[k*16 + 8], r[k*16 + 10]);
                      tmp[6] = _mm_unpacklo_epi16(r[k*16 + 12], r[k*16 + 14]);
                      tmp[7] = _mm_unpackhi_epi16(r[k*16 + 12], r[k*16 + 14]);

                      tmp2[0] = _mm_unpacklo_epi32(tmp[0], tmp[2]);
                      tmp2[1] = _mm_unpackhi_epi32(tmp[0], tmp[2]);
                      tmp2[2] = _mm_unpacklo_epi32(tmp[1], tmp[3]);
                      tmp2[3] = _mm_unpackhi_epi32(tmp[1], tmp[3]);
                      tmp2[4] = _mm_unpacklo_epi32(tmp[4], tmp[6]);
                      tmp2[5] = _mm_unpackhi_epi32(tmp[4], tmp[6]);
                      tmp2[6] = _mm_unpacklo_epi32(tmp[5], tmp[7]);
                      tmp2[7] = _mm_unpackhi_epi32(tmp[5], tmp[7]);

                      tmp[0] = _mm_unpacklo_epi64(tmp2[0], tmp2[4]);
                      tmp[1] = _mm_unpackhi_epi64(tmp2[0], tmp2[4]);
                      tmp[2] = _mm_unpacklo_epi64(tmp2[1], tmp2[5]);
                      tmp[3] = _mm_unpackhi_epi64(tmp2[1], tmp2[5]);
                      tmp[4] = _mm_unpacklo_epi64(tmp2[2], tmp2[6]);
                      tmp[5] = _mm_unpackhi_epi64(tmp2[2], tmp2[6]);
                      tmp[6] = _mm_unpacklo_epi64(tmp2[3], tmp2[7]);
                      tmp[7] = _mm_unpackhi_epi64(tmp2[3], tmp2[7]);

                      _mm_store_si128((__m128i *) (dst+k*8+0), tmp[0]);
                      _mm_store_si128((__m128i *) (dst+k*8+64), tmp[1]);
                      _mm_store_si128((__m128i *) (dst+k*8+128), tmp[2]);
                      _mm_store_si128((__m128i *) (dst+k*8+192), tmp[3]);
                      _mm_store_si128((__m128i *) (dst+k*8+256), tmp[4]);
                      _mm_store_si128((__m128i *) (dst+k*8+320), tmp[5]);
                      _mm_store_si128((__m128i *) (dst+k*8+384), tmp[6]);
                      _mm_store_si128((__m128i *) (dst+k*8+448), tmp[7]);
                    }
                      src+=8;
                      dst += 512;
                }
}


void
vvc_inverse_dct_ii_64_sse(const int16_t *src, int16_t *dst, ptrdiff_t src_stride,
                          int num_lines, int line_brk, int shift)
{
    /*FIXME it would be more efficient to get per reduced lines speialised
      functions instead off checking */
    if (num_lines>>3){
      if (line_brk > 16) {
          vvc_inverse_dct_ii_64_sse_8lines(src, dst, src_stride, num_lines&0xF8, line_brk, shift);
      }
      else if (line_brk > 8){
          vvc_inverse_dct_ii_64_sse_8lines(src, dst, src_stride, num_lines&0xF8, line_brk, shift);
      }
      else {
          vvc_inverse_dct_ii_64_sse_8lines(src, dst, src_stride, num_lines&0xF8, line_brk, shift);
      }
        src += (num_lines & 0xF8);
        dst += (num_lines >> 3) << 9;
    }

    // if (!(num_lines & 0x7)) return;

    if (num_lines & 0x4){
      vvc_inverse_dct_ii_64(src, dst, src_stride, num_lines & 0x4, line_brk, shift);
      src += 4;
      dst += 256;
    }

    if (num_lines & 0x2){
      vvc_inverse_dct_ii_64(src, dst, src_stride, num_lines & 0x2, line_brk, shift);
      src += 2;
      dst += 128;
    }

    if (num_lines & 0x1){
        vvc_inverse_dct_ii_64(src, dst, src_stride, num_lines & 0x1, line_brk, shift);
    }
}


void vvc_inverse_dst_vii_4_sse(const int16_t *src, int16_t *dst, ptrdiff_t src_stride,
                               int num_lines, int line_brk, int shift){
    inverse_sse2_B4(src, dst, src_stride, shift, num_lines, DST_VII_4);
}

void vvc_inverse_dst_vii_8_sse(const int16_t *src, int16_t *dst, ptrdiff_t src_stride,
                               int num_lines, int line_brk, int shift){
    inverse_sse2_B8(src, dst, src_stride, shift, num_lines, DST_VII_8);
}

void vvc_inverse_dst_vii_16_sse(const int16_t *src, int16_t *dst, ptrdiff_t src_stride,
                                int num_lines, int line_brk, int shift){
    inverse_sse2_B16(src, dst, src_stride, shift, num_lines, DST_VII_16);
}

void vvc_inverse_dst_vii_32_sse(const int16_t *src, int16_t *dst, ptrdiff_t src_stride,
                                int num_lines, int line_brk, int shift){
    inverse_sse2_B32(src, dst, src_stride, shift, num_lines, DST_VII_32);
}

void vvc_inverse_dct_viii_4_sse(const int16_t *src, int16_t *dst, ptrdiff_t src_stride,
                                int num_lines, int line_brk, int shift){
    inverse_sse2_B4(src, dst, src_stride, shift, num_lines, DCT_VIII_4);
}

void vvc_inverse_dct_viii_8_sse(const int16_t *src, int16_t *dst, ptrdiff_t src_stride,
                                int num_lines, int line_brk, int shift){
    inverse_sse2_B8(src, dst, src_stride, shift, num_lines, DCT_VIII_8);
}

void vvc_inverse_dct_viii_16_sse(const int16_t *src, int16_t *dst, ptrdiff_t src_stride,
                                 int num_lines, int line_brk, int shift){
    inverse_sse2_B16(src, dst, src_stride, shift, num_lines, DCT_VIII_16);
}

void vvc_inverse_dct_viii_32_sse(const int16_t *src, int16_t *dst, ptrdiff_t src_stride,
                                 int num_lines, int line_brk, int shift){
    inverse_sse2_B32(src, dst, src_stride, shift, num_lines, DCT_VIII_32);
}


void
vvc_inverse_dct_ii_dc_sse(int16_t *const dst, int log2_tb_w, int log2_tb_h,
                          int dc_val)
{

    int i;//, j;
    int tb_w = 1 << log2_tb_w;
    int tb_h = 1 << log2_tb_h;
    int clip_min = -(1 << 15);
    int clip_max = (1 << 15)-1;
    int16_t * _dst = (int16_t *)dst;
    int value = (((dc_val + 1) >> 1) + 8) >> 4;
    value = ov_clip(value, clip_min, clip_max);
    __m128i x0 = _mm_set1_epi16(value);
    switch (log2_tb_w){
        case 3:{
        for (i = 0; i < tb_h; ++i){
            _mm_store_si128((__m128i *)_dst, x0);
            _dst += tb_w;
        }
                   break;
               }
        #if 1
        case 4:{
        for (i = 0; i < tb_h; ++i){
            _mm_store_si128((__m128i *)_dst, x0);
            _mm_store_si128((__m128i *)&_dst[8], x0);
            _dst += tb_w;
        }
                   break;
               }
        case 5:{
        for (i = 0; i < tb_h; ++i){
            _mm_store_si128((__m128i *)_dst, x0);
            _mm_store_si128((__m128i *)&_dst[8], x0);
            _mm_store_si128((__m128i *)&_dst[16], x0);
            _mm_store_si128((__m128i *)&_dst[24], x0);
            _dst += tb_w;
        }
                   break;
               }
               #endif
        default:
            vvc_inverse_dct_ii_dc(dst, log2_tb_w, log2_tb_h, dc_val);
    }
}


void rcn_init_tr_functions_sse(struct RCNFunctions *const rcn_funcs){
  rcn_funcs->tr.func[DST_VII][2] = &vvc_inverse_dst_vii_4_sse;
  rcn_funcs->tr.func[DST_VII][3] = &vvc_inverse_dst_vii_8_sse;
  rcn_funcs->tr.func[DST_VII][4] = &vvc_inverse_dst_vii_16_sse;
  rcn_funcs->tr.func[DST_VII][5] = &vvc_inverse_dst_vii_32_sse;
  //
  rcn_funcs->tr.func[DCT_VIII][2] = &vvc_inverse_dct_viii_4_sse;
  rcn_funcs->tr.func[DCT_VIII][3] = &vvc_inverse_dct_viii_8_sse;
  rcn_funcs->tr.func[DCT_VIII][4] = &vvc_inverse_dct_viii_16_sse;
  rcn_funcs->tr.func[DCT_VIII][5] = &vvc_inverse_dct_viii_32_sse;
  //
  rcn_funcs->tr.func[DCT_II][1] = &vvc_inverse_dct_ii_2_sse;
  rcn_funcs->tr.func[DCT_II][2] = &vvc_inverse_dct_ii_4_sse;
  rcn_funcs->tr.func[DCT_II][3] = &vvc_inverse_dct_ii_8_sse;
  rcn_funcs->tr.func[DCT_II][4] = &vvc_inverse_dct_ii_16_sse;
  rcn_funcs->tr.func[DCT_II][5] = &vvc_inverse_dct_ii_32_sse;
  rcn_funcs->tr.func[DCT_II][6] = &vvc_inverse_dct_ii_64_sse;

  rcn_funcs->tr.dc = &vvc_inverse_dct_ii_dc_sse;
}
