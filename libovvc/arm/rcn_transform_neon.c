#include <stdio.h>
#include <stddef.h>
#include <arm_neon.h>

#include "ovutils.h"
#include "data_rcn_transform.h"
#include "rcn_transform.h"

static void afficherVecteur16NEON128(int8x16_t cible)
{
    int8_t test[16];
    vst1q_s8(test, cible);
//    printf(" ________ ________ ________ ________ ________ ________ ________ ________\n|        |        |        |        |        |        |        |        |\n");
    printf("| %06d | %06d | %06d | %06d | %06d | %06d | %06d | %06d | %06d | %06d | %06d | %06d | %06d | %06d | %06d | %06d |\n", test[0], test[1], test[2], test[3], test[4], test[5], test[6], test[7], test[8], test[9], test[10], test[11], test[12], test[13], test[14], test[15]);
//    printf("|________|________|________|________|________|________|________|________|\n");
}

static void afficherVecteur8NEON128(int16x8_t cible)
{
    int16_t test[8];
    vst1q_s16(test, cible);
//    printf(" ________ ________ ________ ________ ________ ________ ________ ________\n|        |        |        |        |        |        |        |        |\n");
    printf("| %06d | %06d | %06d | %06d | %06d | %06d | %06d | %06d |\n", test[0], test[1], test[2], test[3], test[4], test[5], test[6], test[7]);
//    printf("|________|________|________|________|________|________|________|________|\n");
}

static void afficherVecteur4NEON128(int32x4_t cible)
{
    int32_t test[4];
    vst1q_s32(test, cible);
//    printf(" ________ ________ ________ ________ ________ ________ ________ ________\n|        |        |        |        |        |        |        |        |\n");
    printf("| %06d | %06d | %06d | %06d |\n", test[0], test[1], test[2], test[3]);
//    printf("|________|________|________|________|________|________|________|________|\n");
}

static inline void
transpose4x4(int32x4_t *x, int32x4_t *r)
{
    int64x2_t tmp[4];

    tmp[0] = (int64x2_t) vzip1q_s32(x[0], x[1]);
    tmp[1] = (int64x2_t) vzip2q_s32(x[0], x[1]);
    tmp[2] = (int64x2_t) vzip1q_s32(x[2], x[3]);
    tmp[3] = (int64x2_t) vzip2q_s32(x[2], x[3]);

    r[0] = (int32x4_t) vzip1q_s64(tmp[0], tmp[2]);
    r[1] = (int32x4_t) vzip2q_s64(tmp[0], tmp[2]);
    r[2] = (int32x4_t) vzip1q_s64(tmp[1], tmp[3]);
    r[3] = (int32x4_t) vzip2q_s64(tmp[1], tmp[3]);
}

static inline void
transpose4x4_step2(int32x4_t *x, int32x4_t *r)
{
    int64x2_t tmp[4];

    tmp[0] = (int64x2_t) vzip1q_s32(x[0], x[2]);
    tmp[1] = (int64x2_t) vzip2q_s32(x[0], x[2]);
    tmp[2] = (int64x2_t) vzip1q_s32(x[4], x[6]);
    tmp[3] = (int64x2_t) vzip2q_s32(x[4], x[6]);

    r[0] = (int32x4_t) vzip1q_s64(tmp[0], tmp[2]);
    r[1] = (int32x4_t) vzip2q_s64(tmp[0], tmp[2]);
    r[2] = (int32x4_t) vzip1q_s64(tmp[1], tmp[3]);
    r[3] = (int32x4_t) vzip2q_s64(tmp[1], tmp[3]);
}

static inline void
matMult4x4_red1(int16x8_t *x, int16x8_t *d, int32x4_t *r)
{
    r[0] = vmull_s16(vget_low_s16(x[0]), vget_low_s16(d[0]));
    r[1] = vmull_s16(vget_high_s16(x[0]), vget_high_s16(d[0]));
}

static inline void
matMult4x4_red2(int16x8_t *x, int16x8_t *d, int32x4_t *r)
{
    int32x4_t m[4];//, a[4];
    m[0] = vmull_s16(vget_low_s16(x[0]), vget_low_s16(d[0]));
    m[1] = vmull_s16(vget_high_s16(x[0]), vget_high_s16(d[0]));
    m[2] = vmull_s16(vget_low_s16(x[1]), vget_low_s16(d[1]));
    m[3] = vmull_s16(vget_high_s16(x[1]), vget_high_s16(d[1]));

    r[0] = vaddq_s32(m[0], m[2]);
    r[1] = vaddq_s32(m[1], m[3]);
}

static inline void
matMult4x4(int16x8_t *x, int16x8_t *d, int32x4_t *r)
{
    int32x4_t m[8], a[4];

    m[0] = vmull_s16(vget_low_s16(x[0]), vget_low_s16(d[0]));
    m[1] = vmull_s16(vget_high_s16(x[0]), vget_high_s16(d[0]));
    m[2] = vmull_s16(vget_low_s16(x[1]), vget_low_s16(d[1]));
    m[3] = vmull_s16(vget_high_s16(x[1]), vget_high_s16(d[1]));
    m[4] = vmull_s16(vget_low_s16(x[2]), vget_low_s16(d[2]));
    m[5] = vmull_s16(vget_high_s16(x[2]), vget_high_s16(d[2]));
    m[6] = vmull_s16(vget_low_s16(x[3]), vget_low_s16(d[3]));
    m[7] = vmull_s16(vget_high_s16(x[3]), vget_high_s16(d[3]));

    a[0] = vaddq_s32(m[0], m[2]);
    a[1] = vaddq_s32(m[1], m[3]);
    a[2] = vaddq_s32(m[4], m[6]);
    a[3] = vaddq_s32(m[5], m[7]);

    r[0] = vaddq_s32(a[0], a[2]);
    r[1] = vaddq_s32(a[1], a[3]);
}

static inline void
matMult8x8_red1(int16x8_t *x, int16x8_t *d, int32x4_t *r)
{
    r[0] = vmull_s16(vget_low_s16(x[0]), vget_low_s16(d[0]));
    r[1] = vmull_s16(vget_high_s16(x[0]), vget_high_s16(d[0]));
}

static inline void
matMult8x8_red2(int16x8_t *x, int16x8_t *d, int32x4_t *r)
{
    int32x4_t m[16];//, a[8], b[4];

    m[ 0] = vmull_s16(vget_low_s16(x[0]), vget_low_s16(d[0]));
    m[ 1] = vmull_s16(vget_high_s16(x[0]), vget_high_s16(d[0]));
    m[ 2] = vmull_s16(vget_low_s16(x[1]), vget_low_s16(d[1]));
    m[ 3] = vmull_s16(vget_high_s16(x[1]), vget_high_s16(d[1]));

    r[0] = vaddq_s32(m[ 0], m[ 2]);
    r[1] = vaddq_s32(m[ 1], m[ 3]);
}

static inline void
matMult8x8_red4(int16x8_t *x, int16x8_t *d, int32x4_t *r)
{
    int32x4_t m[16], a[8];//, b[4];

    m[ 0] = vmull_s16(vget_low_s16(x[0]), vget_low_s16(d[0]));
    m[ 1] = vmull_s16(vget_high_s16(x[0]), vget_high_s16(d[0]));
    m[ 2] = vmull_s16(vget_low_s16(x[1]), vget_low_s16(d[1]));
    m[ 3] = vmull_s16(vget_high_s16(x[1]), vget_high_s16(d[1]));
    m[ 4] = vmull_s16(vget_low_s16(x[2]), vget_low_s16(d[2]));
    m[ 5] = vmull_s16(vget_high_s16(x[2]), vget_high_s16(d[2]));
    m[ 6] = vmull_s16(vget_low_s16(x[3]), vget_low_s16(d[3]));
    m[ 7] = vmull_s16(vget_high_s16(x[3]), vget_high_s16(d[3]));

    a[0] = vaddq_s32(m[ 0], m[ 2]);
    a[1] = vaddq_s32(m[ 1], m[ 3]);
    a[2] = vaddq_s32(m[ 4], m[ 6]);
    a[3] = vaddq_s32(m[ 5], m[ 7]);

    r[0] = vaddq_s32(a[0], a[2]);
    r[1] = vaddq_s32(a[1], a[3]);
}

static inline void
matMult8x8(int16x8_t *x, int16x8_t *d, int32x4_t *r)
{
    int32x4_t m[16], a[8], b[4];

    m[0] = vmull_s16(vget_low_s16(x[0]), vget_low_s16(d[0]));
    m[1] = vmull_s16(vget_high_s16(x[0]), vget_high_s16(d[0]));
    m[2] = vmull_s16(vget_low_s16(x[1]), vget_low_s16(d[1]));
    m[3] = vmull_s16(vget_high_s16(x[1]), vget_high_s16(d[1]));
    m[4] = vmull_s16(vget_low_s16(x[2]), vget_low_s16(d[2]));
    m[5] = vmull_s16(vget_high_s16(x[2]), vget_high_s16(d[2]));
    m[6] = vmull_s16(vget_low_s16(x[3]), vget_low_s16(d[3]));
    m[7] = vmull_s16(vget_high_s16(x[3]), vget_high_s16(d[3]));

    a[0] = vaddq_s32(m[ 0], m[ 2]);
    a[1] = vaddq_s32(m[ 1], m[ 3]);
    a[2] = vaddq_s32(m[ 4], m[ 6]);
    a[3] = vaddq_s32(m[ 5], m[ 7]);

    b[0] = vaddq_s32(a[0], a[2]);
    b[1] = vaddq_s32(a[1], a[3]);

    m[ 8] = vmull_s16(vget_low_s16(x[4]), vget_low_s16(d[4]));
    m[ 9] = vmull_s16(vget_high_s16(x[4]), vget_high_s16(d[4]));
    m[10] = vmull_s16(vget_low_s16(x[5]), vget_low_s16(d[5]));
    m[11] = vmull_s16(vget_high_s16(x[5]), vget_high_s16(d[5]));
    m[12] = vmull_s16(vget_low_s16(x[6]), vget_low_s16(d[6]));
    m[13] = vmull_s16(vget_high_s16(x[6]), vget_high_s16(d[6]));
    m[14] = vmull_s16(vget_low_s16(x[7]), vget_low_s16(d[7]));
    m[15] = vmull_s16(vget_high_s16(x[7]), vget_high_s16(d[7]));

    a[4] = vaddq_s32(m[ 8], m[10]);
    a[5] = vaddq_s32(m[ 9], m[11]);
    a[6] = vaddq_s32(m[12], m[14]);
    a[7] = vaddq_s32(m[13], m[15]);

    b[2] = vaddq_s32(a[4], a[6]);
    b[3] = vaddq_s32(a[5], a[7]);

    r[0] = vaddq_s32(b[0], b[2]);
    r[1] = vaddq_s32(b[1], b[3]);
}
static inline void
matMult16x16(int16x8_t *x, int16x8_t *d, int32x4_t *r)
{
    int32x4_t m[32], a[16], b[8], c[4];

    m[0] = vmull_s16(vget_low_s16(x[0]), vget_low_s16(d[0]));
    m[1] = vmull_s16(vget_high_s16(x[0]), vget_high_s16(d[0]));
    m[2] = vmull_s16(vget_low_s16(x[1]), vget_low_s16(d[1]));
    m[3] = vmull_s16(vget_high_s16(x[1]), vget_high_s16(d[1]));
    m[4] = vmull_s16(vget_low_s16(x[2]), vget_low_s16(d[2]));
    m[5] = vmull_s16(vget_high_s16(x[2]), vget_high_s16(d[2]));
    m[6] = vmull_s16(vget_low_s16(x[3]), vget_low_s16(d[3]));
    m[7] = vmull_s16(vget_high_s16(x[3]), vget_high_s16(d[3]));

    a[0] = vaddq_s32(m[ 0], m[ 2]);
    a[1] = vaddq_s32(m[ 1], m[ 3]);
    a[2] = vaddq_s32(m[ 4], m[ 6]);
    a[3] = vaddq_s32(m[ 5], m[ 7]);

    b[0] = vaddq_s32(a[ 0], a[ 2]);
    b[1] = vaddq_s32(a[ 1], a[ 3]);

    m[ 8] = vmull_s16(vget_low_s16(x[4]), vget_low_s16(d[4]));
    m[ 9] = vmull_s16(vget_high_s16(x[4]), vget_high_s16(d[4]));
    m[10] = vmull_s16(vget_low_s16(x[5]), vget_low_s16(d[5]));
    m[11] = vmull_s16(vget_high_s16(x[5]), vget_high_s16(d[5]));
    m[12] = vmull_s16(vget_low_s16(x[6]), vget_low_s16(d[6]));
    m[13] = vmull_s16(vget_high_s16(x[6]), vget_high_s16(d[6]));
    m[14] = vmull_s16(vget_low_s16(x[7]), vget_low_s16(d[7]));
    m[15] = vmull_s16(vget_high_s16(x[7]), vget_high_s16(d[7]));

    a[ 4] = vaddq_s32(m[ 8], m[10]);
    a[ 5] = vaddq_s32(m[ 9], m[11]);
    a[ 6] = vaddq_s32(m[12], m[14]);
    a[ 7] = vaddq_s32(m[13], m[15]);

    b[ 2] = vaddq_s32(a[ 4], a[ 6]);
    b[ 3] = vaddq_s32(a[ 5], a[ 7]);

    c[ 0] = vaddq_s32(b[ 0], b[ 2]);
    c[ 1] = vaddq_s32(b[ 1], b[ 3]);

    m[16] = vmull_s16(vget_low_s16(x[8]), vget_low_s16(d[8]));
    m[17] = vmull_s16(vget_high_s16(x[8]), vget_high_s16(d[8]));
    m[18] = vmull_s16(vget_low_s16(x[9]), vget_low_s16(d[9]));
    m[19] = vmull_s16(vget_high_s16(x[9]), vget_high_s16(d[9]));
    m[20] = vmull_s16(vget_low_s16(x[10]), vget_low_s16(d[10]));
    m[21] = vmull_s16(vget_high_s16(x[10]), vget_high_s16(d[10]));
    m[22] = vmull_s16(vget_low_s16(x[11]), vget_low_s16(d[11]));
    m[23] = vmull_s16(vget_high_s16(x[11]), vget_high_s16(d[11]));

    a[ 8] = vaddq_s32(m[16], m[18]);
    a[ 9] = vaddq_s32(m[17], m[19]);
    a[10] = vaddq_s32(m[20], m[22]);
    a[11] = vaddq_s32(m[21], m[23]);

    b[ 4] = vaddq_s32(a[ 8], a[10]);
    b[ 5] = vaddq_s32(a[ 9], a[11]);

    m[24] = vmull_s16(vget_low_s16(x[12]), vget_low_s16(d[12]));
    m[25] = vmull_s16(vget_high_s16(x[12]), vget_high_s16(d[12]));
    m[26] = vmull_s16(vget_low_s16(x[13]), vget_low_s16(d[13]));
    m[27] = vmull_s16(vget_high_s16(x[13]), vget_high_s16(d[13]));
    m[28] = vmull_s16(vget_low_s16(x[14]), vget_low_s16(d[14]));
    m[29] = vmull_s16(vget_high_s16(x[14]), vget_high_s16(d[14]));
    m[30] = vmull_s16(vget_low_s16(x[15]), vget_low_s16(d[15]));
    m[31] = vmull_s16(vget_high_s16(x[15]), vget_high_s16(d[15]));

    a[12] = vaddq_s32(m[24], m[26]);
    a[13] = vaddq_s32(m[25], m[27]);
    a[14] = vaddq_s32(m[28], m[30]);
    a[15] = vaddq_s32(m[29], m[31]);

    b[ 6] = vaddq_s32(a[12], a[14]);
    b[ 7] = vaddq_s32(a[13], a[15]);

    c[ 2] = vaddq_s32(b[ 4], b[ 6]);
    c[ 3] = vaddq_s32(b[ 5], b[ 7]);

    r[0] = vaddq_s32(c[0], c[2]);
    r[1] = vaddq_s32(c[1], c[3]);
}

static inline void
matMult16x16_red2(int16x8_t *x, int16x8_t *d, int32x4_t *r)
{
    int32x4_t m[4];

    m[0] = vmull_s16(vget_low_s16(x[0]), vget_low_s16(d[0]));
    m[1] = vmull_s16(vget_high_s16(x[0]), vget_high_s16(d[0]));
    m[2] = vmull_s16(vget_low_s16(x[1]), vget_low_s16(d[1]));
    m[3] = vmull_s16(vget_high_s16(x[1]), vget_high_s16(d[1]));

    r[ 0] = vaddq_s32(m[ 0], m[ 2]);
    r[ 1] = vaddq_s32(m[ 1], m[ 3]);
}

static inline void
matMult16x16_red4(int16x8_t *x, int16x8_t *d, int32x4_t *r)
{
    int32x4_t m[8], a[4];

    m[0] = vmull_s16(vget_low_s16(x[0]), vget_low_s16(d[0]));
    m[1] = vmull_s16(vget_high_s16(x[0]), vget_high_s16(d[0]));
    m[2] = vmull_s16(vget_low_s16(x[1]), vget_low_s16(d[1]));
    m[3] = vmull_s16(vget_high_s16(x[1]), vget_high_s16(d[1]));
    m[4] = vmull_s16(vget_low_s16(x[2]), vget_low_s16(d[2]));
    m[5] = vmull_s16(vget_high_s16(x[2]), vget_high_s16(d[2]));
    m[6] = vmull_s16(vget_low_s16(x[3]), vget_low_s16(d[3]));
    m[7] = vmull_s16(vget_high_s16(x[3]), vget_high_s16(d[3]));

    a[0] = vaddq_s32(m[ 0], m[ 2]);
    a[1] = vaddq_s32(m[ 1], m[ 3]);
    a[2] = vaddq_s32(m[ 4], m[ 6]);
    a[3] = vaddq_s32(m[ 5], m[ 7]);

    r[0] = vaddq_s32(a[ 0], a[ 2]);
    r[1] = vaddq_s32(a[ 1], a[ 3]);
}

static inline void
matMult16x16_red8(int16x8_t *x, int16x8_t *d, int32x4_t *r)
{
    int32x4_t m[16], a[8], b[4];

    m[0] = vmull_s16(vget_low_s16(x[0]), vget_low_s16(d[0]));
    m[1] = vmull_s16(vget_high_s16(x[0]), vget_high_s16(d[0]));
    m[2] = vmull_s16(vget_low_s16(x[1]), vget_low_s16(d[1]));
    m[3] = vmull_s16(vget_high_s16(x[1]), vget_high_s16(d[1]));
    m[4] = vmull_s16(vget_low_s16(x[2]), vget_low_s16(d[2]));
    m[5] = vmull_s16(vget_high_s16(x[2]), vget_high_s16(d[2]));
    m[6] = vmull_s16(vget_low_s16(x[3]), vget_low_s16(d[3]));
    m[7] = vmull_s16(vget_high_s16(x[3]), vget_high_s16(d[3]));

    a[0] = vaddq_s32(m[ 0], m[ 2]);
    a[1] = vaddq_s32(m[ 1], m[ 3]);
    a[2] = vaddq_s32(m[ 4], m[ 6]);
    a[3] = vaddq_s32(m[ 5], m[ 7]);

    b[0] = vaddq_s32(a[ 0], a[ 2]);
    b[1] = vaddq_s32(a[ 1], a[ 3]);

    m[ 8] = vmull_s16(vget_low_s16(x[4]), vget_low_s16(d[4]));
    m[ 9] = vmull_s16(vget_high_s16(x[4]), vget_high_s16(d[4]));
    m[10] = vmull_s16(vget_low_s16(x[5]), vget_low_s16(d[5]));
    m[11] = vmull_s16(vget_high_s16(x[5]), vget_high_s16(d[5]));
    m[12] = vmull_s16(vget_low_s16(x[6]), vget_low_s16(d[6]));
    m[13] = vmull_s16(vget_high_s16(x[6]), vget_high_s16(d[6]));
    m[14] = vmull_s16(vget_low_s16(x[7]), vget_low_s16(d[7]));
    m[15] = vmull_s16(vget_high_s16(x[7]), vget_high_s16(d[7]));

    a[ 4] = vaddq_s32(m[ 8], m[10]);
    a[ 5] = vaddq_s32(m[ 9], m[11]);
    a[ 6] = vaddq_s32(m[12], m[14]);
    a[ 7] = vaddq_s32(m[13], m[15]);

    b[ 2] = vaddq_s32(a[ 4], a[ 6]);
    b[ 3] = vaddq_s32(a[ 5], a[ 7]);

    r[ 0] = vaddq_s32(b[ 0], b[ 2]);
    r[ 1] = vaddq_s32(b[ 1], b[ 3]);
}

static inline void
transform_4x2(int16x8_t *x, int16x8_t *d, int32x4_t *r)
{
    int32x4_t m[4];
    m[0] = vmull_s16(vget_low_s16(x[0]), vget_low_s16(d[0]));
    m[1] = vmull_s16(vget_high_s16(x[0]), vget_high_s16(d[0]));
    m[2] = vmull_s16(vget_low_s16(x[1]), vget_low_s16(d[1]));
    m[3] = vmull_s16(vget_high_s16(x[1]), vget_high_s16(d[1]));

    r[0] = vaddq_s32(m[0], m[2]);
    r[1] = vaddq_s32(m[1], m[3]);
}

static inline void
transform_4x4(int16x8_t *x, int16x8_t *d, int32x4_t *r)
{
    int32x4_t m[4];

    transform_4x2(x, d+0, m+0);
    transform_4x2(x, d+2, m+2);

    r[0] = vaddq_s32(m[0], m[1]);
    r[1] = vaddq_s32(m[2], m[3]);
    r[2] = vsubq_s32(m[0], m[1]);
    r[3] = vsubq_s32(m[2], m[3]);
}
static inline void
transform_4x8(int16x8_t *x, int16x8_t *d, int32x4_t *r)
{
    int32x4_t a[8];
    int32x4_t m[4];

    m[0] = vmull_s16(vget_low_s16(x[0]), vget_low_s16(d[0]));
    m[1] = vmull_s16(vget_high_s16(x[0]), vget_high_s16(d[0]));
    m[2] = vmull_s16(vget_low_s16(x[2]), vget_low_s16(d[4]));
    m[3] = vmull_s16(vget_high_s16(x[2]), vget_high_s16(d[4]));

    a[0] = vaddq_s32(m[0], m[2]);
    a[1] = vaddq_s32(m[1], m[3]);

    m[0] = vmull_s16(vget_low_s16(x[0]), vget_low_s16(d[1]));
    m[1] = vmull_s16(vget_high_s16(x[0]), vget_high_s16(d[1]));
    m[2] = vmull_s16(vget_low_s16(x[2]), vget_low_s16(d[5]));
    m[3] = vmull_s16(vget_high_s16(x[2]), vget_high_s16(d[5]));

    a[2] = vaddq_s32(m[0], m[2]);
    a[3] = vaddq_s32(m[1], m[3]);

    m[0] = vmull_s16(vget_low_s16(x[1]), vget_low_s16(d[2]));
    m[1] = vmull_s16(vget_high_s16(x[1]), vget_high_s16(d[2]));
    m[2] = vmull_s16(vget_low_s16(x[3]), vget_low_s16(d[6]));
    m[3] = vmull_s16(vget_high_s16(x[3]), vget_high_s16(d[6]));

    a[4] = vaddq_s32(m[0], m[2]);
    a[5] = vaddq_s32(m[1], m[3]);

    m[0] = vmull_s16(vget_low_s16(x[1]), vget_low_s16(d[3]));
    m[1] = vmull_s16(vget_high_s16(x[1]), vget_high_s16(d[3]));
    m[2] = vmull_s16(vget_low_s16(x[3]), vget_low_s16(d[7]));
    m[3] = vmull_s16(vget_high_s16(x[3]), vget_high_s16(d[7]));

    a[6] = vaddq_s32(m[0], m[2]);
    a[7] = vaddq_s32(m[1], m[3]);

    r[0] = vaddq_s32(a[0], a[4]);
    r[1] = vaddq_s32(a[1], a[5]);
    r[2] = vaddq_s32(a[2], a[6]);
    r[3] = vaddq_s32(a[3], a[7]);
    r[4] = vsubq_s32(a[2], a[6]);
    r[5] = vsubq_s32(a[3], a[7]);
    r[6] = vsubq_s32(a[0], a[4]);
    r[7] = vsubq_s32(a[1], a[5]);
}

static inline void
dct2_4x2(int16x8_t *x, int16x8_t *d, int32x4_t *r)
{
    int32x4_t m[2], a[2];

    transform_4x2(x,d,m);

    a[0] = vaddq_s32(m[0], m[1]);
    a[1] = vsubq_s32(m[0], m[1]);

    r[0] = (int32x4_t)vzip1q_s64((int64x2_t) a[0], (int64x2_t) a[1]);
    r[1] = (int32x4_t)vzip2q_s64((int64x2_t) a[0], (int64x2_t) a[1]);
}

static inline void
dct2_4x4_red1(int16x8_t *x, int16x8_t *d, int32x4_t *r)
{
    int32x4_t m[4], t[4];

    /*FIXME we could use swizzling to make more
     efficient use of madd here */
    m[0] = vmull_s16(vget_low_s16(x[0]), vget_low_s16(d[0]));
    m[1] = vmull_s16(vget_high_s16(x[0]), vget_high_s16(d[0]));
    m[2] = vmull_s16(vget_low_s16(x[0]), vget_low_s16(d[2]));
    m[3] = vmull_s16(vget_high_s16(x[0]), vget_high_s16(d[2]));

    t[0] = vaddq_s32(m[0], m[1]);
    t[1] = vaddq_s32(m[2], m[3]);
    t[2] = vsubq_s32(m[0], m[1]);
    t[3] = vsubq_s32(m[2], m[3]);

    /*FIXME transpose and shuffle are
     redundant blend might a better choice*/
    transpose4x4(t, r);

    r[0] = vcombine_s32(vget_low_s32(r[0]), vrev64_s32(vget_high_s32(r[0])));//_mm_shuffle_epi32(r[0], 0xB4);
    r[1] = vcombine_s32(vrev64_s32(vget_low_s32(r[1])), vget_high_s32(r[1]));//_mm_shuffle_epi32(r[1], 0xE1);
    r[2] = vcombine_s32(vget_low_s32(r[2]), vrev64_s32(vget_high_s32(r[2])));//_mm_shuffle_epi32(r[2], 0xB4);
    r[3] = vcombine_s32(vrev64_s32(vget_low_s32(r[3])), vget_high_s32(r[3]));//_mm_shuffle_epi32(r[3], 0xE1);
}

static inline void
dct2_4x4_red2(int16x8_t *x, int16x8_t *d, int32x4_t *r)
{
    int32x4_t m[4], t[4];

    /*FIXME we could use swizzling to make more
     efficient use of madd here */
    t[0] = vmull_s16(vget_low_s16(x[0]), vget_low_s16(d[0]));
    t[1] = vmull_s16(vget_high_s16(x[0]), vget_high_s16(d[0]));
    t[2] = vmull_s16(vget_low_s16(x[1]), vget_low_s16(d[1]));
    t[3] = vmull_s16(vget_high_s16(x[1]), vget_high_s16(d[1]));

    m[0] = vaddq_s32(t[0], t[2]);
    m[1] = vaddq_s32(t[1], t[3]);

    t[0] = vmull_s16(vget_low_s16(x[0]), vget_low_s16(d[2]));
    t[1] = vmull_s16(vget_high_s16(x[0]), vget_high_s16(d[2]));
    t[2] = vmull_s16(vget_low_s16(x[1]), vget_low_s16(d[3]));
    t[3] = vmull_s16(vget_high_s16(x[1]), vget_high_s16(d[3]));

    m[2] = vaddq_s32(t[0], t[2]);
    m[3] = vaddq_s32(t[1], t[3]);

    t[0] = vaddq_s32(m[0], m[1]);
    t[1] = vaddq_s32(m[2], m[3]);
    t[2] = vsubq_s32(m[0], m[1]);
    t[3] = vsubq_s32(m[2], m[3]);

    transpose4x4(t, r);

    r[0] = vcombine_s32(vget_low_s32(r[0]), vrev64_s32(vget_high_s32(r[0])));//_mm_shuffle_epi32(r[0], 0xB4);
    r[1] = vcombine_s32(vrev64_s32(vget_low_s32(r[1])), vget_high_s32(r[1]));//_mm_shuffle_epi32(r[1], 0xE1);
    r[2] = vcombine_s32(vget_low_s32(r[2]), vrev64_s32(vget_high_s32(r[2])));//_mm_shuffle_epi32(r[2], 0xB4);
    r[3] = vcombine_s32(vrev64_s32(vget_low_s32(r[3])), vget_high_s32(r[3]));//_mm_shuffle_epi32(r[3], 0xE1);
}

static inline void
dct2_4x4(int16x8_t *x, int16x8_t *d, int32x4_t *r)
{
    int32x4_t m[4];

    transform_4x4(x, d, m);

    transpose4x4(m, r);

    r[0] = vcombine_s32(vget_low_s32(r[0]), vrev64_s32(vget_high_s32(r[0])));//_mm_shuffle_epi32(r[0], 0xB4);
    r[1] = vcombine_s32(vrev64_s32(vget_low_s32(r[1])), vget_high_s32(r[1]));//_mm_shuffle_epi32(r[1], 0xE1);
    r[2] = vcombine_s32(vget_low_s32(r[2]), vrev64_s32(vget_high_s32(r[2])));//_mm_shuffle_epi32(r[2], 0xB4);
    r[3] = vcombine_s32(vrev64_s32(vget_low_s32(r[3])), vget_high_s32(r[3]));//_mm_shuffle_epi32(r[3], 0xE1);
}

static inline void
dct2_4x8(int16x8_t *x, int16x8_t *d, int32x4_t *r)
{
    int32x4_t m[8];

    transform_4x8(x, d, m);

    transpose4x4_step2(m  , r  );
    transpose4x4_step2(m+1, r+4);
}

static inline void
dct2_8x2(int16x8_t *x, int16x8_t *d, int32x4_t *r)
{
    int32x4_t e[2], o[2];

    dct2_4x2(x+0, d+0, e);

    matMult4x4(x+2, d+2, o);

    r[0] = vaddq_s32(e[0], o[0]);
    r[1] = vsubq_s32(e[0], o[0]);
    r[2] = vaddq_s32(e[1], o[1]);
    r[3] = vsubq_s32(e[1], o[1]);

    r[0] = vcombine_s32(vget_low_s32(r[0]), vrev64_s32(vget_high_s32(r[0])));//_mm_shuffle_epi32(r[0], 0xB4);
    r[1] = vcombine_s32(vget_high_s32(r[1]), vrev64_s32(vget_low_s32(r[1])));//_mm_shuffle_epi32(r[1], 0x1E);
    r[2] = vcombine_s32(vget_low_s32(r[2]), vrev64_s32(vget_high_s32(r[2])));//_mm_shuffle_epi32(r[2], 0xB4);
    r[3] = vcombine_s32(vget_high_s32(r[3]), vrev64_s32(vget_low_s32(r[3])));//_mm_shuffle_epi32(r[3], 0x1E);
}

static inline void
dct2_8x4_red1(int16x8_t *x, int16x8_t *d, int32x4_t *r)
{
    int32x4_t e[4];//, o[4];

    dct2_4x4_red1(x+0, d+0, e);

    r[0] = e[0];
    r[2] = e[1];
    r[4] = e[2];
    r[6] = e[3];

    r[1] = vcombine_s32(vrev64_s32(vget_high_s32(e[0])), vrev64_s32(vget_low_s32(e[0])));//_mm_shuffle_epi32(e[1], 0x1B);
    r[3] = vcombine_s32(vrev64_s32(vget_high_s32(e[1])), vrev64_s32(vget_low_s32(e[1])));//_mm_shuffle_epi32(e[3], 0x1B);
    r[5] = vcombine_s32(vrev64_s32(vget_high_s32(e[2])), vrev64_s32(vget_low_s32(e[2])));//_mm_shuffle_epi32(e[5], 0x1B);
    r[7] = vcombine_s32(vrev64_s32(vget_high_s32(e[3])), vrev64_s32(vget_low_s32(e[3])));//_mm_shuffle_epi32(e[7], 0x1B);
}

static inline void
dct2_8x4_red2(int16x8_t *x, int16x8_t *d, int32x4_t *r)
{
    int32x4_t e[4], o[4];

    dct2_4x4_red1(x+0, d+0, e);

    matMult4x4_red1(x+2, d+4, o+0);
    matMult4x4_red1(x+6, d+4, o+2);

    o[0] = vcombine_s32(vget_low_s32(o[0]), vrev64_s32(vget_high_s32(o[0])));//_mm_shuffle_epi32(o[0], 0xB4);
    o[1] = vcombine_s32(vget_low_s32(o[1]), vrev64_s32(vget_high_s32(o[1])));//_mm_shuffle_epi32(o[1], 0xB4);
    o[2] = vcombine_s32(vget_low_s32(o[2]), vrev64_s32(vget_high_s32(o[2])));//_mm_shuffle_epi32(o[2], 0xB4);
    o[3] = vcombine_s32(vget_low_s32(o[3]), vrev64_s32(vget_high_s32(o[3])));//_mm_shuffle_epi32(o[3], 0xB4);

    r[0] = vaddq_s32(e[0], o[0]);
    r[1] = vsubq_s32(e[0], o[0]);
    r[2] = vaddq_s32(e[1], o[1]);
    r[3] = vsubq_s32(e[1], o[1]);
    r[4] = vaddq_s32(e[2], o[2]);
    r[5] = vsubq_s32(e[2], o[2]);
    r[6] = vaddq_s32(e[3], o[3]);
    r[7] = vsubq_s32(e[3], o[3]);

    r[1] = vcombine_s32(vrev64_s32(vget_high_s32(r[1])), vrev64_s32(vget_low_s32(r[1])));//_mm_shuffle_epi32(r[1], 0x1B);
    r[3] = vcombine_s32(vrev64_s32(vget_high_s32(r[3])), vrev64_s32(vget_low_s32(r[3])));//_mm_shuffle_epi32(r[3], 0x1B);
    r[5] = vcombine_s32(vrev64_s32(vget_high_s32(r[5])), vrev64_s32(vget_low_s32(r[5])));//_mm_shuffle_epi32(r[5], 0x1B);
    r[7] = vcombine_s32(vrev64_s32(vget_high_s32(r[7])), vrev64_s32(vget_low_s32(r[7])));//_mm_shuffle_epi32(r[7], 0x1B);
}

static inline void
dct2_8x4_red4(int16x8_t *x, int16x8_t *d, int32x4_t *r)
{
    int32x4_t e[4], o[4];

    dct2_4x4_red2(x+0, d+0, e);

    matMult4x4_red2(x+2, d+4, o+0);
    matMult4x4_red2(x+6, d+4, o+2);

    o[0] = vcombine_s32(vget_low_s32(o[0]), vrev64_s32(vget_high_s32(o[0])));//_mm_shuffle_epi32(o[0], 0xB4);
    o[1] = vcombine_s32(vget_low_s32(o[1]), vrev64_s32(vget_high_s32(o[1])));//_mm_shuffle_epi32(o[1], 0xB4);
    o[2] = vcombine_s32(vget_low_s32(o[2]), vrev64_s32(vget_high_s32(o[2])));//_mm_shuffle_epi32(o[2], 0xB4);
    o[3] = vcombine_s32(vget_low_s32(o[3]), vrev64_s32(vget_high_s32(o[3])));//_mm_shuffle_epi32(o[3], 0xB4);

    r[0] = vaddq_s32(e[0], o[0]);
    r[1] = vsubq_s32(e[0], o[0]);
    r[2] = vaddq_s32(e[1], o[1]);
    r[3] = vsubq_s32(e[1], o[1]);
    r[4] = vaddq_s32(e[2], o[2]);
    r[5] = vsubq_s32(e[2], o[2]);
    r[6] = vaddq_s32(e[3], o[3]);
    r[7] = vsubq_s32(e[3], o[3]);

    r[1] = vcombine_s32(vrev64_s32(vget_high_s32(r[1])), vrev64_s32(vget_low_s32(r[1])));//_mm_shuffle_epi32(r[1], 0x1B);
    r[3] = vcombine_s32(vrev64_s32(vget_high_s32(r[3])), vrev64_s32(vget_low_s32(r[3])));//_mm_shuffle_epi32(r[3], 0x1B);
    r[5] = vcombine_s32(vrev64_s32(vget_high_s32(r[5])), vrev64_s32(vget_low_s32(r[5])));//_mm_shuffle_epi32(r[5], 0x1B);
    r[7] = vcombine_s32(vrev64_s32(vget_high_s32(r[7])), vrev64_s32(vget_low_s32(r[7])));//_mm_shuffle_epi32(r[7], 0x1B);
}

static inline void
dct2_8x4(int16x8_t *x, int16x8_t *d, int32x4_t *r)
{
    int32x4_t e[4], o[4];

    dct2_4x4(x+0, d+0, e);

    matMult4x4(x+2, d+4, o+0);
    matMult4x4(x+6, d+4, o+2);

    o[0] = vcombine_s32(vget_low_s32(o[0]), vrev64_s32(vget_high_s32(o[0])));//_mm_shuffle_epi32(o[0], 0xB4);
    o[1] = vcombine_s32(vget_low_s32(o[1]), vrev64_s32(vget_high_s32(o[1])));//_mm_shuffle_epi32(o[1], 0xB4);
    o[2] = vcombine_s32(vget_low_s32(o[2]), vrev64_s32(vget_high_s32(o[2])));//_mm_shuffle_epi32(o[2], 0xB4);
    o[3] = vcombine_s32(vget_low_s32(o[3]), vrev64_s32(vget_high_s32(o[3])));//_mm_shuffle_epi32(o[3], 0xB4);

    r[0] = vaddq_s32(e[0], o[0]);
    r[1] = vsubq_s32(e[0], o[0]);
    r[2] = vaddq_s32(e[1], o[1]);
    r[3] = vsubq_s32(e[1], o[1]);
    r[4] = vaddq_s32(e[2], o[2]);
    r[5] = vsubq_s32(e[2], o[2]);
    r[6] = vaddq_s32(e[3], o[3]);
    r[7] = vsubq_s32(e[3], o[3]);

    r[1] = vcombine_s32(vrev64_s32(vget_high_s32(r[1])), vrev64_s32(vget_low_s32(r[1])));//_mm_shuffle_epi32(r[1], 0x1B);
    r[3] = vcombine_s32(vrev64_s32(vget_high_s32(r[3])), vrev64_s32(vget_low_s32(r[3])));//_mm_shuffle_epi32(r[3], 0x1B);
    r[5] = vcombine_s32(vrev64_s32(vget_high_s32(r[5])), vrev64_s32(vget_low_s32(r[5])));//_mm_shuffle_epi32(r[5], 0x1B);
    r[7] = vcombine_s32(vrev64_s32(vget_high_s32(r[7])), vrev64_s32(vget_low_s32(r[7])));//_mm_shuffle_epi32(r[7], 0x1B);
}

static inline void
dct2_8x8(int16x8_t *x, int16x8_t *d, int32x4_t *r)
{
    int32x4_t e[8], o[8], oo[8];

    dct2_4x8(x+0, d+0, e);

    matMult4x4(x+4, d+ 8, oo+0);
    matMult4x4(x+4, d+12, oo+2);
    matMult4x4(x+4, d+16, oo+4);
    matMult4x4(x+4, d+20, oo+6);

    transpose4x4_step2(oo+0, o+0);
    transpose4x4_step2(oo+1, o+4);

    for(int k=0; k<8; k++){
        r[  k] = vaddq_s32(e[k], o[k]);
        int32x4_t tmp = vsubq_s32(e[k],  o[k]);
        r[8+k] = vcombine_s32(vrev64_s32(vget_high_s32(tmp)), vrev64_s32(vget_low_s32(tmp)));//_mm_shuffle_epi32(tmp,  0x1B);
    }
}
static inline void
dct2_16x2(int16x8_t *x, int16x8_t *d, int32x4_t *r)
{
    int32x4_t e[4], o[4];

    dct2_8x2(x+0, d+0, e);

    matMult8x8(x+ 6, d+6, o+0);
    matMult8x8(x+14, d+6, o+2);

    r[0] = vaddq_s32(e[0], o[0]);
    r[1] = vaddq_s32(e[1], o[1]);
    r[2] = vsubq_s32(e[1], o[1]);
    r[3] = vsubq_s32(e[0], o[0]);
    r[4] = vaddq_s32(e[2], o[2]);
    r[5] = vaddq_s32(e[3], o[3]);
    r[6] = vsubq_s32(e[3], o[3]);
    r[7] = vsubq_s32(e[2], o[2]);

    r[2] = vcombine_s32(vrev64_s32(vget_high_s32(r[2])), vrev64_s32(vget_low_s32(r[2])));//_mm_shuffle_epi32(r[2], 0x1B);
    r[3] = vcombine_s32(vrev64_s32(vget_high_s32(r[3])), vrev64_s32(vget_low_s32(r[3])));//_mm_shuffle_epi32(r[3], 0x1B);
    r[6] = vcombine_s32(vrev64_s32(vget_high_s32(r[6])), vrev64_s32(vget_low_s32(r[6])));//_mm_shuffle_epi32(r[6], 0x1B);
    r[7] = vcombine_s32(vrev64_s32(vget_high_s32(r[7])), vrev64_s32(vget_low_s32(r[7])));//_mm_shuffle_epi32(r[7], 0x1B);
}

static inline void
dct2_16x4_red2(int16x8_t *x, int16x8_t *d, int32x4_t *r)
{
    int32x4_t e[8], o[8];
    int16x8_t *odd_x = x + 10;

    dct2_8x4_red1(x+0, d+0, e);

    matMult8x8_red1(odd_x +  0, d+8, o+0);
    matMult8x8_red1(odd_x +  8, d+8, o+2);
    matMult8x8_red1(odd_x + 16, d+8, o+4);
    matMult8x8_red1(odd_x + 24, d+8, o+6);

    for(int k=0; k<8; k+=2){
        r[2*k+0] = vaddq_s32(e[k+0], o[k+0]);
        r[2*k+1] = vaddq_s32(e[k+1], o[k+1]);
        int32x4_t temp = vsubq_s32(e[k+1], o[k+1]);
        r[2*k+2] = vcombine_s32(vrev64_s32(vget_high_s32(temp)), vrev64_s32(vget_low_s32(temp)));//_mm_shuffle_epi32(temp, 0x1B);
        temp = vsubq_s32(e[k+0], o[k+0]);
        r[2*k+3] = vcombine_s32(vrev64_s32(vget_high_s32(temp)), vrev64_s32(vget_low_s32(temp)));//_mm_shuffle_epi32(temp, 0x1B);
    }
}

static inline void
dct2_16x4_red4(int16x8_t *x, int16x8_t *d, int32x4_t *r)
{
    int32x4_t e[8], o[8];
    int16x8_t *odd_x = x + 10;

    dct2_8x4_red2(x+0, d+0, e);

    matMult8x8_red2(odd_x +  0, d+8, o+0);
    matMult8x8_red2(odd_x +  8, d+8, o+2);
    matMult8x8_red2(odd_x + 16, d+8, o+4);
    matMult8x8_red2(odd_x + 24, d+8, o+6);

    for(int k=0; k<8; k+=2){
        r[2*k+0] = vaddq_s32(e[k+0], o[k+0]);
        r[2*k+1] = vaddq_s32(e[k+1], o[k+1]);
        int32x4_t temp = vsubq_s32(e[k+1], o[k+1]);
        r[2*k+2] = vcombine_s32(vrev64_s32(vget_high_s32(temp)), vrev64_s32(vget_low_s32(temp)));//_mm_shuffle_epi32(temp, 0x1B);
        temp = vsubq_s32(e[k+0], o[k+0]);
        r[2*k+3] = vcombine_s32(vrev64_s32(vget_high_s32(temp)), vrev64_s32(vget_low_s32(temp)));//_mm_shuffle_epi32(temp, 0x1B);
    }
}

static inline void
dct2_16x4_red8(int16x8_t *x, int16x8_t *d, int32x4_t *r)
{
    int32x4_t e[8], o[8];
    int16x8_t *odd_x = x + 10;

    dct2_8x4_red4(x+0, d+0, e);

    matMult8x8_red4(odd_x +  0, d+8, o+0);
    matMult8x8_red4(odd_x +  8, d+8, o+2);
    matMult8x8_red4(odd_x + 16, d+8, o+4);
    matMult8x8_red4(odd_x + 24, d+8, o+6);

    for(int k=0; k<8; k+=2){
        r[2*k+0] = vaddq_s32(e[k+0], o[k+0]);
        r[2*k+1] = vaddq_s32(e[k+1], o[k+1]);
        int32x4_t temp = vsubq_s32(e[k+1], o[k+1]);
        r[2*k+2] = vcombine_s32(vrev64_s32(vget_high_s32(temp)), vrev64_s32(vget_low_s32(temp)));//_mm_shuffle_epi32(temp, 0x1B);
        temp = vsubq_s32(e[k+0], o[k+0]);
        r[2*k+3] = vcombine_s32(vrev64_s32(vget_high_s32(temp)), vrev64_s32(vget_low_s32(temp)));//_mm_shuffle_epi32(temp, 0x1B);
    }
}

static inline void
dct2_16x4(int16x8_t *x, int16x8_t *d, int32x4_t *r)
{
    int32x4_t e[8], o[8];

    dct2_8x4(x+0, d+0, e);

    matMult8x8(x+10, d+8, o+0);
    matMult8x8(x+18, d+8, o+2);
    matMult8x8(x+26, d+8, o+4);
    matMult8x8(x+34, d+8, o+6);

    for(int k=0; k<8; k+=2){
        r[2*k+0] = vaddq_s32(e[k+0], o[k+0]);
        r[2*k+1] = vaddq_s32(e[k+1], o[k+1]);
        int32x4_t temp = vsubq_s32(e[k+1], o[k+1]);
        r[2*k+2] = vcombine_s32(vrev64_s32(vget_high_s32(temp)), vrev64_s32(vget_low_s32(temp)));//_mm_shuffle_epi32(temp, 0x1B);
        temp = vsubq_s32(e[k+0], o[k+0]);
        r[2*k+3] = vcombine_s32(vrev64_s32(vget_high_s32(temp)), vrev64_s32(vget_low_s32(temp)));//_mm_shuffle_epi32(temp, 0x1B);
    }
}
static inline void
dct2_16x8(int16x8_t *x, int16x8_t *d, int32x4_t *r)
{
    int32x4_t e[16], o[16];

    dct2_8x8(x+0, d+0, e);

    for(int k=0; k<8; k++){
        matMult8x8(x+8+8*k, d+24, o+2*k);
    }

    for(int k=0; k<8; k++){
        r[4*k+0] = vaddq_s32(e[k+0], o[2*k+0]);
        r[4*k+1] = vaddq_s32(e[k+8], o[2*k+1]);
        int32x4_t temp = vsubq_s32(e[k+8], o[2*k+1]);
        r[4*k+2] = vcombine_s32(vrev64_s32(vget_high_s32(temp)), vrev64_s32(vget_low_s32(temp)));//_mm_shuffle_epi32(temp, 0x1B);
        temp = vsubq_s32(e[k+0], o[2*k+0]);
        r[4*k+3] = vcombine_s32(vrev64_s32(vget_high_s32(temp)), vrev64_s32(vget_low_s32(temp)));//_mm_shuffle_epi32(temp, 0x1B);
    }
}

static inline void
dct2_32x2(int16x8_t *x, int16x8_t *d, int32x4_t *r)
{
    int32x4_t e[8], o[8];

    dct2_16x2(x+0, d+0, e);

    for(int k=0; k<2; k++){
        matMult16x16(x+24+16*k, d+14, o+4*k+0);
        matMult16x16(x+24+16*k, d+30, o+4*k+2);
    }

    for(int k=0; k<2; k++){
        r[8*k+0] = vaddq_s32(e[4*k+0], o[4*k+0]);
        r[8*k+1] = vaddq_s32(e[4*k+1], o[4*k+1]);
        r[8*k+2] = vaddq_s32(e[4*k+2], o[4*k+2]);
        r[8*k+3] = vaddq_s32(e[4*k+3], o[4*k+3]);
        int32x4_t temp = vsubq_s32(e[4*k+3], o[4*k+3]);
        r[8*k+4] = vcombine_s32(vrev64_s32(vget_high_s32(temp)), vrev64_s32(vget_low_s32(temp)));//_mm_shuffle_epi32(temp, 0x1B);
        temp = vsubq_s32(e[4*k+2], o[4*k+2]);
        r[8*k+5] = vcombine_s32(vrev64_s32(vget_high_s32(temp)), vrev64_s32(vget_low_s32(temp)));//_mm_shuffle_epi32(temp, 0x1B);
        temp = vsubq_s32(e[4*k+1], o[4*k+1]);
        r[8*k+6] = vcombine_s32(vrev64_s32(vget_high_s32(temp)), vrev64_s32(vget_low_s32(temp)));//_mm_shuffle_epi32(temp, 0x1B);
        temp = vsubq_s32(e[4*k+0], o[4*k+0]);
        r[8*k+7] = vcombine_s32(vrev64_s32(vget_high_s32(temp)), vrev64_s32(vget_low_s32(temp)));//_mm_shuffle_epi32(temp, 0x1B);
    }
}

static inline void
dct2_32x4_red4(int16x8_t *x, int16x8_t *d, int32x4_t *r)
{
    int32x4_t e[16], o[16];

    dct2_16x4_red2(x+0, d+0, e);

    for(int k=0; k<4; k++){
        matMult16x16_red2(x+74+16*k, d+16, o+4*k+0);
        matMult16x16_red2(x+74+16*k, d+32, o+4*k+2);
    }

    for(int k=0; k<4; k++){
        r[8*k+0] = vaddq_s32(e[4*k+0], o[4*k+0]);
        r[8*k+1] = vaddq_s32(e[4*k+1], o[4*k+1]);
        r[8*k+2] = vaddq_s32(e[4*k+2], o[4*k+2]);
        r[8*k+3] = vaddq_s32(e[4*k+3], o[4*k+3]);
        int32x4_t temp = vsubq_s32(e[4*k+3], o[4*k+3]);
        r[8*k+4] = vcombine_s32(vrev64_s32(vget_high_s32(temp)), vrev64_s32(vget_low_s32(temp)));//_mm_shuffle_epi32(temp, 0x1B);
        temp = vsubq_s32(e[4*k+2], o[4*k+2]);
        r[8*k+5] = vcombine_s32(vrev64_s32(vget_high_s32(temp)), vrev64_s32(vget_low_s32(temp)));//_mm_shuffle_epi32(temp, 0x1B);
        temp = vsubq_s32(e[4*k+1], o[4*k+1]);
        r[8*k+6] = vcombine_s32(vrev64_s32(vget_high_s32(temp)), vrev64_s32(vget_low_s32(temp)));//_mm_shuffle_epi32(temp, 0x1B);
        temp = vsubq_s32(e[4*k+0], o[4*k+0]);
        r[8*k+7] = vcombine_s32(vrev64_s32(vget_high_s32(temp)), vrev64_s32(vget_low_s32(temp)));//_mm_shuffle_epi32(temp, 0x1B);
    }
}

static inline void
dct2_32x4_red8(int16x8_t *x, int16x8_t *d, int32x4_t *r)
{
    int32x4_t e[16], o[16];

    dct2_16x4_red4(x+0, d+0, e);

    for(int k=0; k<4; k++){
        matMult16x16_red4(x+74+16*k, d+16, o+4*k+0);
        matMult16x16_red4(x+74+16*k, d+32, o+4*k+2);
    }

    for(int k=0; k<4; k++){
        r[8*k+0] = vaddq_s32(e[4*k+0], o[4*k+0]);
        r[8*k+1] = vaddq_s32(e[4*k+1], o[4*k+1]);
        r[8*k+2] = vaddq_s32(e[4*k+2], o[4*k+2]);
        r[8*k+3] = vaddq_s32(e[4*k+3], o[4*k+3]);
        int32x4_t temp = vsubq_s32(e[4*k+3], o[4*k+3]);
        r[8*k+4] = vcombine_s32(vrev64_s32(vget_high_s32(temp)), vrev64_s32(vget_low_s32(temp)));//_mm_shuffle_epi32(temp, 0x1B);
        temp = vsubq_s32(e[4*k+2], o[4*k+2]);
        r[8*k+5] = vcombine_s32(vrev64_s32(vget_high_s32(temp)), vrev64_s32(vget_low_s32(temp)));//_mm_shuffle_epi32(temp, 0x1B);
        temp = vsubq_s32(e[4*k+1], o[4*k+1]);
        r[8*k+6] = vcombine_s32(vrev64_s32(vget_high_s32(temp)), vrev64_s32(vget_low_s32(temp)));//_mm_shuffle_epi32(temp, 0x1B);
        temp = vsubq_s32(e[4*k+0], o[4*k+0]);
        r[8*k+7] = vcombine_s32(vrev64_s32(vget_high_s32(temp)), vrev64_s32(vget_low_s32(temp)));//_mm_shuffle_epi32(temp, 0x1B);
    }
}

static inline void
dct2_32x4_red16(int16x8_t *x, int16x8_t *d, int32x4_t *r)
{
    int32x4_t e[16], o[16];

    dct2_16x4_red8(x+0, d+0, e);

    for(int k=0; k<4; k++){
        matMult16x16_red8(x+74+16*k, d+16, o+4*k+0);
        matMult16x16_red8(x+74+16*k, d+32, o+4*k+2);
    }

    for(int k=0; k<4; k++){
        r[8*k+0] = vaddq_s32(e[4*k+0], o[4*k+0]);
        r[8*k+1] = vaddq_s32(e[4*k+1], o[4*k+1]);
        r[8*k+2] = vaddq_s32(e[4*k+2], o[4*k+2]);
        r[8*k+3] = vaddq_s32(e[4*k+3], o[4*k+3]);
        int32x4_t temp = vsubq_s32(e[4*k+3], o[4*k+3]);
        r[8*k+4] = vcombine_s32(vrev64_s32(vget_high_s32(temp)), vrev64_s32(vget_low_s32(temp)));//_mm_shuffle_epi32(temp, 0x1B);
        temp = vsubq_s32(e[4*k+2], o[4*k+2]);
        r[8*k+5] = vcombine_s32(vrev64_s32(vget_high_s32(temp)), vrev64_s32(vget_low_s32(temp)));//_mm_shuffle_epi32(temp, 0x1B);
        temp = vsubq_s32(e[4*k+1], o[4*k+1]);
        r[8*k+6] = vcombine_s32(vrev64_s32(vget_high_s32(temp)), vrev64_s32(vget_low_s32(temp)));//_mm_shuffle_epi32(temp, 0x1B);
        temp = vsubq_s32(e[4*k+0], o[4*k+0]);
        r[8*k+7] = vcombine_s32(vrev64_s32(vget_high_s32(temp)), vrev64_s32(vget_low_s32(temp)));//_mm_shuffle_epi32(temp, 0x1B);
    }
}

static inline void
dct2_32x4(int16x8_t *x, int16x8_t *d, int32x4_t *r)
{
    int32x4_t e[16], o[16];

    dct2_16x4(x+0, d+0, e);

    for(int k=0; k<4; k++){
        matMult16x16(x+74+16*k, d+16, o+4*k+0);
        matMult16x16(x+74+16*k, d+32, o+4*k+2);
    }

    for(int k=0; k<4; k++){
        r[8*k+0] = vaddq_s32(e[4*k+0], o[4*k+0]);
        r[8*k+1] = vaddq_s32(e[4*k+1], o[4*k+1]);
        r[8*k+2] = vaddq_s32(e[4*k+2], o[4*k+2]);
        r[8*k+3] = vaddq_s32(e[4*k+3], o[4*k+3]);
        int32x4_t temp = vsubq_s32(e[4*k+3], o[4*k+3]);
        r[8*k+4] = vcombine_s32(vrev64_s32(vget_high_s32(temp)), vrev64_s32(vget_low_s32(temp)));//_mm_shuffle_epi32(temp, 0x1B);
        temp = vsubq_s32(e[4*k+2], o[4*k+2]);
        r[8*k+5] = vcombine_s32(vrev64_s32(vget_high_s32(temp)), vrev64_s32(vget_low_s32(temp)));//_mm_shuffle_epi32(temp, 0x1B);
        temp = vsubq_s32(e[4*k+1], o[4*k+1]);
        r[8*k+6] = vcombine_s32(vrev64_s32(vget_high_s32(temp)), vrev64_s32(vget_low_s32(temp)));//_mm_shuffle_epi32(temp, 0x1B);
        temp = vsubq_s32(e[4*k+0], o[4*k+0]);
        r[8*k+7] = vcombine_s32(vrev64_s32(vget_high_s32(temp)), vrev64_s32(vget_low_s32(temp)));//_mm_shuffle_epi32(temp, 0x1B);
    }
}

static inline void
dct2_32x8(int16x8_t *x, int16x8_t *d, int32x4_t *r){
    int32x4_t e[32], o[32];

    dct2_16x8(x+0, d+0, e+0);

    for(int k=0; k<8; k++){
        matMult16x16(x+72+16*k, d+32, o+4*k+0);
        matMult16x16(x+72+16*k, d+48, o+4*k+2);
    }

    for(int k=0; k<8; k++){
        r[8*k+0] = vaddq_s32(e[4*k+0], o[4*k+0]);
        r[8*k+1] = vaddq_s32(e[4*k+1], o[4*k+1]);
        r[8*k+2] = vaddq_s32(e[4*k+2], o[4*k+2]);
        r[8*k+3] = vaddq_s32(e[4*k+3], o[4*k+3]);
        int32x4_t temp = vsubq_s32(e[4*k+3], o[4*k+3]);
        r[8*k+4] = vcombine_s32(vrev64_s32(vget_high_s32(temp)), vrev64_s32(vget_low_s32(temp)));//_mm_shuffle_epi32(temp, 0x1B);
        temp = vsubq_s32(e[4*k+2], o[4*k+2]);
        r[8*k+5] = vcombine_s32(vrev64_s32(vget_high_s32(temp)), vrev64_s32(vget_low_s32(temp)));//_mm_shuffle_epi32(temp, 0x1B);
        temp = vsubq_s32(e[4*k+1], o[4*k+1]);
        r[8*k+6] = vcombine_s32(vrev64_s32(vget_high_s32(temp)), vrev64_s32(vget_low_s32(temp)));//_mm_shuffle_epi32(temp, 0x1B);
        temp = vsubq_s32(e[4*k+0], o[4*k+0]);
        r[8*k+7] = vcombine_s32(vrev64_s32(vget_high_s32(temp)), vrev64_s32(vget_low_s32(temp)));//_mm_shuffle_epi32(temp, 0x1B);
    }
}
void vvc_inverse_dct_ii_2_neon(const int16_t *src, int16_t *dst, ptrdiff_t src_stride,
                              int num_lines, int line_brk, int shift)
{
    int16x8_t x1, x2;
    int16x8_t outl, outh;

    int32x4_t add = vdupq_n_s32(1<<(shift-6-1));
    for(int i = 0; i < num_lines / 8; i++){
        //Même chose d'utiliser unpacklo et unpackhi ici
        x1 = vld1q_s16(&src[0]);
        x2 = vld1q_s16(&src[src_stride]);

        int32x4_t x1lo = vmull_n_s16(vget_low_s16((int16x8_t)x1), 1);
        int32x4_t x1hi = vmull_n_s16(vget_high_s16((int16x8_t)x1), 1);

        int32x4_t x2lo = vmull_n_s16(vget_low_s16((int16x8_t)x2), 1);
        int32x4_t x2hi = vmull_n_s16(vget_high_s16((int16x8_t)x2), 1);

        int32x4_t elo = vaddq_s32(x1lo,x2lo);
        int32x4_t ehi = vaddq_s32(x1hi,x2hi);
        int32x4_t olo = vsubq_s32(x1lo,x2lo);
        int32x4_t ohi = vsubq_s32(x1hi,x2hi);

        elo = vaddq_s32(elo, add);
        ehi = vaddq_s32(ehi, add);
        olo = vaddq_s32(olo, add);
        ohi = vaddq_s32(ohi, add);

        elo = vshrq_n_s32(elo, shift-6);
        ehi = vshrq_n_s32(ehi, shift-6);
        olo = vshrq_n_s32(olo, shift-6);
        ohi = vshrq_n_s32(ohi, shift-6);

        int32x4_t outllo = vzip1q_s32(elo,olo);
        int32x4_t outhlo = vzip2q_s32(elo,olo);
        int32x4_t outlhi = vzip1q_s32(ehi,ohi);
        int32x4_t outhhi = vzip2q_s32(ehi,ohi);

        outl = vcombine_s16(vqmovn_s32(outllo), vqmovn_s32(outhlo));
        outh = vcombine_s16(vqmovn_s32(outlhi), vqmovn_s32(outhhi));

        vst1q_s16(&(dst[0]), outl);
        vst1q_s16(&(dst[8]), outh);

        src += 8;
        dst += 16;
    }

    if (!(num_lines & 0x7)) return;

    if (num_lines & 0x4){
        x1 = vcombine_s16(vld1_s16(&src[0]), vcreate_s16(0));
        x2 = vcombine_s16(vld1_s16(&src[src_stride]), vcreate_s16(0));

        int32x4_t x1l = vmull_n_s16(vget_low_s16(x1), 1);
        int32x4_t x2l = vmull_n_s16(vget_low_s16(x2), 1);

        int32x4_t e = vaddq_s32(x1l,x2l);
        int32x4_t o = vsubq_s32(x1l,x2l);

        e = vaddq_s32(e, add);
        o = vaddq_s32(o, add);

        e = vshrq_n_s32(e, shift-6);
        o = vshrq_n_s32(o, shift-6);

        int32x4_t outllo = vzip1q_s32(e,o);
        int32x4_t outhlo = vzip2q_s32(e,o);

        outl = vcombine_s16(vqmovn_s32(outllo), vqmovn_s32(outhlo));
        vst1q_s16(&(dst[0]), outl);
    }

    if (num_lines & 0x2){
      //FIXME Untestable with CTC AI
      int32_t vect1[4] = {src[0], src[0], src[1], src[1]};
      int32_t vect2[4] = {src[src_stride], src[src_stride], src[src_stride + 1], src[src_stride + 1]};
      int32x4_t x1l = vld1q_s32(vect1);
      int32x4_t x2l = vld1q_s32(vect2);

      int32x4_t xo = (int32x4_t)vdupq_n_s64(0xFFFFFFFF00000000);
      x2l = vaddq_s32(veorq_s32(x2l, xo),(int32x4_t)vdupq_n_s64(0x0000000100000000));
      int32x4_t e = vaddq_s32(x1l,x2l);

      int32x4_t add = vdupq_n_s32(1<<(shift-6-1));
      e = vaddq_s32(e, add);
      e = vshrq_n_s32(e, shift-6);
      outl = vcombine_s16(vqmovn_s32(e), vqmovn_s32(e));
      vst1q_s16(&(dst[0]), outl);
    }

    if (num_lines & 0x1){
        vvc_inverse_dct_ii_2(src, dst, src_stride, num_lines & 0x1, line_brk, shift);
    }
}

void
vvc_inverse_dct_ii_4_neon(const int16_t *src, int16_t *dst, ptrdiff_t src_stride,
                        int num_lines, int line_brk, int shift)
{
    int32x4_t add = vdupq_n_s32(1 << (shift - 1));

    for (int j = 0; j < num_lines / 8; j++) {
        int16x8_t x[4], d[8];
        int32x4_t r[8];

        x[0] = vld1q_s16(src + 0 * src_stride);
        x[1] = vld1q_s16(src + 1 * src_stride);
        x[2] = vld1q_s16(src + 2 * src_stride);
        x[3] = vld1q_s16(src + 3 * src_stride);

        d[0] = vdupq_n_s16(DCT_II_4[0]);
        d[1] = vdupq_n_s16(DCT_II_4[1]);
        d[2] = vdupq_n_s16(DCT_II_4[4]);
        d[3] = vdupq_n_s16(DCT_II_4[5]);
        d[4] = vdupq_n_s16(DCT_II_4[8]);
        d[5] = vdupq_n_s16(DCT_II_4[9]);
        d[6] = vdupq_n_s16(DCT_II_4[12]);
        d[7] = vdupq_n_s16(DCT_II_4[13]);

        dct2_4x8(x,d,r);

        for(int i=0; i<8; i+=2){
            int16x8_t o;
            r[i+0] = vaddq_s32(r[i+0], add);
            r[i+1] = vaddq_s32(r[i+1], add);

            r[i+0] = vshrq_n_s32(r[i+0], shift);
            r[i+1] = vshrq_n_s32(r[i+1], shift);

            o = vcombine_s16(vqmovn_s32(r[i+0]),
                     vqmovn_s32(r[i+1]));

            vst1q_s16((dst+i/2*8), o);
        }
        src += 8;
        dst += 32;
    }

    if (!(num_lines & 0x7)) return;

    if (num_lines & 0x4){
        int16x8_t x[2], d[4];
        int32x4_t r[4];

        x[0] = vcombine_s16(
                	vld1_s16(src + 0 * src_stride),
                	vld1_s16(src + 1 * src_stride)
                );
        x[1] = vcombine_s16(
                	vld1_s16(src + 2 * src_stride),
                	vld1_s16(src + 3 * src_stride)
                );

        static const int16_t DCT_II_4_4_neon[8 * 4] = {
            64,  64,  64,  64,  83,  36,  83,  36,
            64, -64,  64, -64,  36, -83,  36, -83,
            64,  64,  64,  64,  36,  83,  36,  83,
            -64,  64, -64,  64, -83,  36, -83,  36
        };

        d[0] = vld1q_s16(DCT_II_4_4_neon + 0);
        d[1] = vld1q_s16(DCT_II_4_4_neon + 8);
        d[2] = vld1q_s16(DCT_II_4_4_neon + 16);
        d[3] = vld1q_s16(DCT_II_4_4_neon + 24);

        dct2_4x4(x,d,r);

        int32x4_t add = vdupq_n_s32(1 << (shift - 1));
        r[0] = vaddq_s32(r[0], add);
        r[1] = vaddq_s32(r[1], add);
        r[2] = vaddq_s32(r[2], add);
        r[3] = vaddq_s32(r[3], add);

        r[0] = vshrq_n_s32(r[0], shift);
        r[1] = vshrq_n_s32(r[1], shift);
        r[2] = vshrq_n_s32(r[2], shift);
        r[3] = vshrq_n_s32(r[3], shift);

        int16x8_t o0 = vcombine_s16(vqmovn_s32(r[0]), vqmovn_s32(r[1]));
        int16x8_t o1 = vcombine_s16(vqmovn_s32(r[2]), vqmovn_s32(r[3]));

        vst1q_s16(dst+0, o0);
        vst1q_s16(dst+8, o1);
    }

    if (num_lines & 0x2){
      int16x8_t x[2], d[2];
      int32x4_t r[2];
        //FIXME: Not use on CTC AI
        x[0] = (int16x8_t) vzip1q_s32(
                	(int32x4_t)vcombine_s16(vld1_s16(src + 0 * src_stride), vdup_n_s16(0)),
                	(int32x4_t)vcombine_s16(vld1_s16(src + 1 * src_stride), vdup_n_s16(0))
                );
        x[0] = vzip1q_s16(x[0], x[0]);
        x[1] = (int16x8_t) vzip1q_s32(
                	(int32x4_t)vcombine_s16(vld1_s16(src + 2 * src_stride), vdup_n_s16(0)),
                	(int32x4_t)vcombine_s16(vld1_s16(src + 3 * src_stride), vdup_n_s16(0))
                );
        x[1] = vzip1q_s16(x[1], x[1]);
        static const int16_t DCT_II_4_2_neon[8 * 2]= {
            64, 64, 64, 64, 83, 36, 83, 36,
            64,-64, 64,-64, 36,-83, 36,-83
        };

        d[0] = vld1q_s16(DCT_II_4_2_neon + 0);
        d[1] = vld1q_s16(DCT_II_4_2_neon + 8);

        dct2_4x2(x,d,r);

        r[0] = vcombine_s32(vget_low_s32(r[0]), vrev64_s32(vget_high_s32(r[0])));//_mm_shuffle_epi32(r[0], 0xB4);
        r[1] = vcombine_s32(vget_low_s32(r[1]), vrev64_s32(vget_high_s32(r[1])));//_mm_shuffle_epi32(r[1], 0xB4);

        int32x4_t add = vdupq_n_s32(1 << (shift - 1));
        r[0] = vaddq_s32(r[0], add);
        r[1] = vaddq_s32(r[1], add);

        r[0] = vshrq_n_s32(r[0], shift);
        r[1] = vshrq_n_s32(r[1], shift);

        int16x8_t out = vcombine_s16(vqmovn_s32(r[0]), vqmovn_s32(r[1]));

        vst1q_s16(dst, out);
    }

    if (num_lines & 0x1){
        vvc_inverse_dct_ii_4(src, dst, src_stride, num_lines & 0x1, line_brk, shift);
    }
}

void
vvc_inverse_dct_ii_8_neon(const int16_t *src, int16_t *dst, ptrdiff_t src_stride,
                         int num_lines, int line_brk, int shift)
{
    int32x4_t add = vdupq_n_s32(1 << (shift - 1));

    for (int j = 0; j < num_lines / 8; j++) {
        int16x8_t x[8], d[24];
        int32x4_t r[16];

        x[0] = vld1q_s16(src + 0 * src_stride);
        x[1] = vld1q_s16(src + 2 * src_stride);
        x[2] = vld1q_s16(src + 4 * src_stride);
        x[3] = vld1q_s16(src + 6 * src_stride);
        x[4] = vld1q_s16(src + 1 * src_stride);
        x[5] = vld1q_s16(src + 3 * src_stride);
        x[6] = vld1q_s16(src + 5 * src_stride);
        x[7] = vld1q_s16(src + 7 * src_stride);

        d[0] = vdupq_n_s16(DCT_II_8[0]);
        d[1] = vdupq_n_s16(DCT_II_8[1]);
        d[2] = vdupq_n_s16(DCT_II_8[16]);
        d[3] = vdupq_n_s16(DCT_II_8[17]);
        d[4] = vdupq_n_s16(DCT_II_8[32]);
        d[5] = vdupq_n_s16(DCT_II_8[33]);
        d[6] = vdupq_n_s16(DCT_II_8[48]);
        d[7] = vdupq_n_s16(DCT_II_8[49]);

        for(int k=0; k<4; k++){
            d[8+4*k+0] = vdupq_n_s16(DCT_II_8[8+k]);
            d[8+4*k+1] = vdupq_n_s16(DCT_II_8[24+k]);
            d[8+4*k+2] = vdupq_n_s16(DCT_II_8[40+k]);
            d[8+4*k+3] = vdupq_n_s16(DCT_II_8[56+k]);
        }

        dct2_8x8(x,d,r);

        for(int i=0; i<8; i++){
            int16x8_t o;

            r[i+0] = vaddq_s32(r[i+0], add);
            r[i+8] = vaddq_s32(r[i+8], add);

            r[i+0] = vshrq_n_s32(r[i+0], shift);
            r[i+8] = vshrq_n_s32(r[i+8], shift);

            o = vcombine_s16(vqmovn_s32(r[i+0]), vqmovn_s32(r[i+8]));

            vst1q_s16(dst+i*8, o);
        }
        src += 8;
        dst += 64;
    }

    if (num_lines & 0x4){
        int16x8_t x[10], d[8];
        int32x4_t r[8];

        x[0] = vcombine_s16(
                	vld1_s16(src + 0 * src_stride),
                	vld1_s16(src + 2 * src_stride)
                );
        x[1] = vcombine_s16(
                	vld1_s16(src + 4 * src_stride),
                	vld1_s16(src + 6 * src_stride)
                );

        x[2] = vcombine_s16(
                vdup_n_s16(src[1 * src_stride + 0]),
                vdup_n_s16(src[1 * src_stride + 1])
                );
        x[3] = vcombine_s16(
                vdup_n_s16(src[3 * src_stride + 0]),
                vdup_n_s16(src[3 * src_stride + 1])
                );
        x[4] = vcombine_s16(
                vdup_n_s16(src[5 * src_stride + 0]),
                vdup_n_s16(src[5 * src_stride + 1])
                );
        x[5] = vcombine_s16(
                vdup_n_s16(src[7 * src_stride + 0]),
                vdup_n_s16(src[7 * src_stride + 1])
                );

        x[6] = vcombine_s16(
                vdup_n_s16(src[1 * src_stride + 2]),
                vdup_n_s16(src[1 * src_stride + 3])
                );
        x[7] = vcombine_s16(
                vdup_n_s16(src[3 * src_stride + 2]),
                vdup_n_s16(src[3 * src_stride + 3])
                );
        x[8] = vcombine_s16(
                vdup_n_s16(src[5 * src_stride + 2]),
                vdup_n_s16(src[5 * src_stride + 3])
                );
        x[9] = vcombine_s16(
                vdup_n_s16(src[7 * src_stride + 2]),
                vdup_n_s16(src[7 * src_stride + 3])
                );

        static const int16_t DCT_II_8_4_neon[8 * 8] = {
            64,  64,  64,  64,  83,  36,  83,  36,
            64, -64,  64, -64,  36, -83,  36, -83,
            64,  64,  64,  64,  36,  83,  36,  83,
            -64, 64, -64,  64, -83,  36, -83,  36,
            89,  75,  18,  50,  89,  75,  18,  50,
            75, -18, -50, -89,  75, -18, -50, -89,
            50, -89,  75,  18,  50, -89,  75,  18,
            18, -50, -89,  75,  18, -50, -89,  75
        };

        d[0] = vld1q_s16(DCT_II_8_4_neon + 0);
        d[1] = vld1q_s16(DCT_II_8_4_neon + 8);
        d[2] = vld1q_s16(DCT_II_8_4_neon + 16);
        d[3] = vld1q_s16(DCT_II_8_4_neon + 24);
        d[4] = vld1q_s16(DCT_II_8_4_neon + 32);
        d[5] = vld1q_s16(DCT_II_8_4_neon + 40);
        d[6] = vld1q_s16(DCT_II_8_4_neon + 48);
        d[7] = vld1q_s16(DCT_II_8_4_neon + 56);

        dct2_8x4(x, d, r);

        int32x4_t add = vdupq_n_s32(1 << (shift - 1));
        for (int i=0; i<8; i+=2){
            int16x8_t o;

            r[i+0] = vaddq_s32(r[i+0], add);
            r[i+1] = vaddq_s32(r[i+1], add);

            r[i+0] = vshrq_n_s32(r[i+0], shift);
            r[i+1] = vshrq_n_s32(r[i+1], shift);

            o = vcombine_s16(vqmovn_s32(r[i+0]), vqmovn_s32(r[i+1]));

            vst1q_s16(dst+i/2*8, o);
        }
    }

    if (num_lines & 0x2){
        int16x8_t x[6];
        int16x8_t d[6];
        int32x4_t r[4];

        x[0] = (int16x8_t) vzip1q_s32(
                	(int32x4_t)vcombine_s16(vld1_s16(src + 0 * src_stride), vdup_n_s16(0)),
                	(int32x4_t)vcombine_s16(vld1_s16(src + 2 * src_stride), vdup_n_s16(0))
                );
        x[0] = vzip1q_s16(x[0], x[0]);
        x[1] = (int16x8_t) vzip1q_s32(
                	(int32x4_t)vcombine_s16(vld1_s16(src + 4 * src_stride), vdup_n_s16(0)),
                	(int32x4_t)vcombine_s16(vld1_s16(src + 6 * src_stride), vdup_n_s16(0))
                );
        x[1] = vzip1q_s16(x[1], x[1]);

        x[2] = vcombine_s16(
                vdup_n_s16(src[1 * src_stride    ]),
                vdup_n_s16(src[1 * src_stride + 1])
                );
        x[3] = vcombine_s16(
                vdup_n_s16(src[3 * src_stride    ]),
                vdup_n_s16(src[3 * src_stride + 1])
                );
        x[4] = vcombine_s16(
                vdup_n_s16(src[5 * src_stride    ]),
                vdup_n_s16(src[5 * src_stride + 1])
                );
        x[5] = vcombine_s16(
                vdup_n_s16(src[7 * src_stride    ]),
                vdup_n_s16(src[7 * src_stride + 1])
                );

        static const int16_t DCT_II_8_2_neon[8 * 6] = {
            64,  64,  64,  64,  83,  36,  83,  36,
            64, -64,  64, -64,  36, -83,  36, -83,
            89,  75,  18,  50,  89,  75,  18,  50,
            75, -18, -50, -89,  75, -18, -50, -89,
            50, -89,  75,  18,  50, -89,  75,  18,
            18, -50, -89,  75,  18, -50, -89,  75
        };

        d[0] = vld1q_s16(DCT_II_8_2_neon +  0);
        d[1] = vld1q_s16(DCT_II_8_2_neon +  8);
        d[2] = vld1q_s16(DCT_II_8_2_neon + 16);
        d[3] = vld1q_s16(DCT_II_8_2_neon + 24);
        d[4] = vld1q_s16(DCT_II_8_2_neon + 32);
        d[5] = vld1q_s16(DCT_II_8_2_neon + 40);

        dct2_8x2(x, d, r);
        int32x4_t add = vdupq_n_s32(1 << (shift - 1));
        r[0] = vaddq_s32(r[0], add);
        r[1] = vaddq_s32(r[1], add);
        r[2] = vaddq_s32(r[2], add);
        r[3] = vaddq_s32(r[3], add);

        r[0] = vshrq_n_s32(r[0], shift);
        r[1] = vshrq_n_s32(r[1], shift);
        r[2] = vshrq_n_s32(r[2], shift);
        r[3] = vshrq_n_s32(r[3], shift);

        int16x8_t o0 = vcombine_s16(vqmovn_s32(r[0]), vqmovn_s32(r[1]));
        int16x8_t o1 = vcombine_s16(vqmovn_s32(r[2]), vqmovn_s32(r[3]));

        vst1q_s16(dst+0, o0);
        vst1q_s16(dst+8, o1);
    }

    if (num_lines & 0x1){
        vvc_inverse_dct_ii_8(src, dst, src_stride, num_lines & 0x1, line_brk, shift);
    }
}

void
vvc_inverse_dct_ii_16_neon(const int16_t *src, int16_t *dst, ptrdiff_t src_stride,
                          int num_lines, int line_brk, int shift)
{
    int32x4_t add = vdupq_n_s32(1 << (shift - 1));

    for (int j = 0; j < num_lines / 8; j++) {
        int16x8_t x[72], d[32];
        int32x4_t r[32];

        x[0] = vld1q_s16(src +  0 * src_stride);
        x[1] = vld1q_s16(src +  4 * src_stride);
        x[2] = vld1q_s16(src +  8 * src_stride);
        x[3] = vld1q_s16(src + 12 * src_stride);
        x[4] = vld1q_s16(src +  2 * src_stride);
        x[5] = vld1q_s16(src +  6 * src_stride);
        x[6] = vld1q_s16(src + 10 * src_stride);
        x[7] = vld1q_s16(src + 14 * src_stride);

        for (int k=0; k<8; k++){
            for (int l=0; l<8; l++)
                x[8+8*k+l] = vdupq_n_s16(src[(2*l+1) * src_stride + k]);
        }

        d[0] = vdupq_n_s16(DCT_II_16[  0]);
        d[1] = vdupq_n_s16(DCT_II_16[  1]);
        d[2] = vdupq_n_s16(DCT_II_16[ 64]);
        d[3] = vdupq_n_s16(DCT_II_16[ 65]);
        d[4] = vdupq_n_s16(DCT_II_16[128]);
        d[5] = vdupq_n_s16(DCT_II_16[129]);
        d[6] = vdupq_n_s16(DCT_II_16[192]);
        d[7] = vdupq_n_s16(DCT_II_16[193]);

        for (int k = 0; k < 4; k++) {
            d[8 + 4 * k + 0] = vdupq_n_s16(DCT_II_16[ 32 + k]);
            d[8 + 4 * k + 1] = vdupq_n_s16(DCT_II_16[ 96 + k]);
            d[8 + 4 * k + 2] = vdupq_n_s16(DCT_II_16[160 + k]);
            d[8 + 4 * k + 3] = vdupq_n_s16(DCT_II_16[224 + k]);
        }

        static const int16_t DCT_II_16_8_neon[8 * 8] = {
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
            d[24+k] = vld1q_s16(DCT_II_16_8_neon + 8*k);
        }

        dct2_16x8(x, d, r);

        for (int i = 0; i < 32; i+=2) {
            int16x8_t o;

            r[i + 0] = vaddq_s32(r[i + 0], add);
            r[i + 1] = vaddq_s32(r[i + 1], add);

            r[i + 0] = vshrq_n_s32(r[i + 0], shift);
            r[i + 1] = vshrq_n_s32(r[i + 1], shift);

            o = vcombine_s16(vqmovn_s32(r[i+0]), vqmovn_s32(r[i+1]));

            vst1q_s16(dst + i * 8/2, o);
        }
        src += 8;
        dst += 128;
    }

    if (!(num_lines & 0x7)) return;

    if (num_lines & 0x4){
        int16x8_t x[42], d[16];
        int32x4_t r[16];

        x[0] = vcombine_s16(
                vld1_s16(src + 0 * src_stride),
                vld1_s16(src + 4 * src_stride)
                );
        x[1] = vcombine_s16(
                vld1_s16(src +  8 * src_stride),
                vld1_s16(src + 12 * src_stride)
                );

        x[2] = vcombine_s16(
                vdup_n_s16(src[2 * src_stride + 0]),
                vdup_n_s16(src[2 * src_stride + 1])
                );
        x[3] = vcombine_s16(
                vdup_n_s16(src[6 * src_stride + 0]),
                vdup_n_s16(src[6 * src_stride + 1])
                );
        x[4] = vcombine_s16(
                vdup_n_s16(src[10 * src_stride + 0]),
                vdup_n_s16(src[10 * src_stride + 1])
                );
        x[5] = vcombine_s16(
                vdup_n_s16(src[14 * src_stride + 0]),
                vdup_n_s16(src[14 * src_stride + 1])
                );

        x[6] = vcombine_s16(
                vdup_n_s16(src[2 * src_stride + 2]),
                vdup_n_s16(src[2 * src_stride + 3])
                );
        x[7] = vcombine_s16(
                vdup_n_s16(src[6 * src_stride + 2]),
                vdup_n_s16(src[6 * src_stride + 3])
                );
        x[8] = vcombine_s16(
                vdup_n_s16(src[10 * src_stride + 2]),
                vdup_n_s16(src[10 * src_stride + 3])
                );
        x[9] = vcombine_s16(
                vdup_n_s16(src[14 * src_stride + 2]),
                vdup_n_s16(src[14 * src_stride + 3])
                );

        for (int k=0; k<4; k++){
            for (int l=0; l<8; l++)
                x[10+8*k+l] = vdupq_n_s16(src[(2*l+1) * src_stride + k]);
        }

        static const int16_t DCT_II_16_4_neon[8 * 16] = {
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
            d[k] = vld1q_s16(DCT_II_16_4_neon + 8*k);
        }

        dct2_16x4(x, d, r);

        int32x4_t add = vdupq_n_s32(1 << (shift - 1));
        for (int i = 0; i < 16; i += 2) {
            int16x8_t o;

            r[i + 0] = vaddq_s32(r[i + 0], add);
            r[i + 1] = vaddq_s32(r[i + 1], add);

            r[i + 0] = vshrq_n_s32(r[i + 0], shift);
            r[i + 1] = vshrq_n_s32(r[i + 1], shift);

            o = vcombine_s16(vqmovn_s32(r[i+0]), vqmovn_s32(r[i+1]));

            vst1q_s16(dst + i / 2 * 8, o);
        }
    }

    if (num_lines & 0x2){
        int16x8_t x[24], d[14];
        int32x4_t r[8];

        x[0] = (int16x8_t) vzip1q_s32(
                	(int32x4_t)vcombine_s16(vld1_s16(src + 0 * src_stride), vdup_n_s16(0)),
                	(int32x4_t)vcombine_s16(vld1_s16(src + 4 * src_stride), vdup_n_s16(0))
                );
        x[0] = vzip1q_s16(x[0], x[0]);
        x[1] = (int16x8_t) vzip1q_s32(
                	(int32x4_t)vcombine_s16(vld1_s16(src + 8 * src_stride), vdup_n_s16(0)),
                	(int32x4_t)vcombine_s16(vld1_s16(src + 12 * src_stride), vdup_n_s16(0))
                );
        x[1] = vzip1q_s16(x[1], x[1]);

        x[2] = vcombine_s16(
                vdup_n_s16(src[2 * src_stride]),
                vdup_n_s16(src[2 * src_stride + 1])
                );
        x[3] = vcombine_s16(
                vdup_n_s16(src[6 * src_stride]),
                vdup_n_s16(src[6 * src_stride + 1])
                );
        x[4] = vcombine_s16(
                vdup_n_s16(src[10 * src_stride]),
                vdup_n_s16(src[10 * src_stride + 1])
                );
        x[5] = vcombine_s16(
                vdup_n_s16(src[14 * src_stride]),
                vdup_n_s16(src[14 * src_stride + 1])
                );

        for (int k=0; k<2; k++){
            x[6+8*k+0] = vdupq_n_s16(src[ 1 * src_stride + k]);
            x[6+8*k+1] = vdupq_n_s16(src[ 3 * src_stride + k]);
            x[6+8*k+2] = vdupq_n_s16(src[ 5 * src_stride + k]);
            x[6+8*k+3] = vdupq_n_s16(src[ 7 * src_stride + k]);
            x[6+8*k+4] = vdupq_n_s16(src[ 9 * src_stride + k]);
            x[6+8*k+5] = vdupq_n_s16(src[11 * src_stride + k]);
            x[6+8*k+6] = vdupq_n_s16(src[13 * src_stride + k]);
            x[6+8*k+7] = vdupq_n_s16(src[15 * src_stride + k]);
        }

        static const int16_t DCT_II_16_2_neon[8 * 14] = {
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
            d[k] = vld1q_s16(DCT_II_16_2_neon + 8*k);
        }

        dct2_16x2(x, d, r);

        int32x4_t add = vdupq_n_s32(1 << (shift - 1));

        for (int i = 0; i < 8; i += 2) {
            int16x8_t o;

            r[i + 0] =  vaddq_s32(r[i + 0], add);
            r[i + 1] =  vaddq_s32(r[i + 1], add);

            r[i + 0] = vshrq_n_s32(r[i + 0], shift);
            r[i + 1] = vshrq_n_s32(r[i + 1], shift);

            o = vcombine_s16(vqmovn_s32(r[i+0]), vqmovn_s32(r[i+1]));

            vst1q_s16(dst + i / 2 * 8, o);
        }
    }

    if (num_lines & 0x1){
        vvc_inverse_dct_ii_16(src, dst, src_stride, num_lines & 0x1, line_brk, shift);
    }
}

static inline void
dct2_4x8_red1(int16x8_t *x, int16x8_t *d, int32x4_t *r)
{
    int32x4_t a[4], tmp[4];

    a[0] = vmull_s16(vget_low_s16(x[0]), vget_low_s16(d[0]));
    a[1] = vmull_s16(vget_high_s16(x[0]), vget_high_s16(d[0]));
    a[2] = vmull_s16(vget_low_s16(x[0]), vget_low_s16(d[1]));
    a[3] = vmull_s16(vget_high_s16(x[0]), vget_high_s16(d[1]));

    tmp[0] = vzip1q_s32(a[0], a[2]);
    tmp[1] = vzip2q_s32(a[0], a[2]);
    tmp[2] = vzip1q_s32(a[2], a[0]);
    tmp[3] = vzip2q_s32(a[2], a[0]);

    r[0] = (int32x4_t) vzip1q_s64((int64x2_t) tmp[0], (int64x2_t) tmp[2]);
    r[1] = (int32x4_t) vzip2q_s64((int64x2_t) tmp[0], (int64x2_t) tmp[2]);
    r[2] = (int32x4_t) vzip1q_s64((int64x2_t) tmp[1], (int64x2_t) tmp[3]);
    r[3] = (int32x4_t) vzip2q_s64((int64x2_t) tmp[1], (int64x2_t) tmp[3]);

    tmp[0] = vzip1q_s32(a[1], a[3]);
    tmp[1] = vzip2q_s32(a[1], a[3]);
    tmp[2] = vzip1q_s32(a[3], a[1]);
    tmp[3] = vzip2q_s32(a[3], a[1]);

    r[4] = (int32x4_t) vzip1q_s64((int64x2_t) tmp[0], (int64x2_t) tmp[2]);
    r[5] = (int32x4_t) vzip2q_s64((int64x2_t) tmp[0], (int64x2_t) tmp[2]);
    r[6] = (int32x4_t) vzip1q_s64((int64x2_t) tmp[1], (int64x2_t) tmp[3]);
    r[7] = (int32x4_t) vzip2q_s64((int64x2_t) tmp[1], (int64x2_t) tmp[3]);
}

static inline void
dct2_8x8_red2(int16x8_t *x, int16x8_t *d, int32x4_t *r)
{
    int32x4_t e[8], o[8], oo[8];

    dct2_4x8_red1(x+0, d+0, e);

    int32x4_t tmp[4];

    oo[0] = vmull_s16(vget_low_s16(x[4]), vget_low_s16(d[8]));
    oo[1] = vmull_s16(vget_high_s16(x[4]), vget_high_s16(d[8]));
    oo[2] = vmull_s16(vget_low_s16(x[4]), vget_low_s16(d[12]));
    oo[3] = vmull_s16(vget_high_s16(x[4]), vget_high_s16(d[12]));

    oo[4] = vmull_s16(vget_low_s16(x[4]), vget_low_s16(d[16]));
    oo[5] = vmull_s16(vget_high_s16(x[4]), vget_high_s16(d[16]));
    oo[6] = vmull_s16(vget_low_s16(x[4]), vget_low_s16(d[20]));
    oo[7] = vmull_s16(vget_high_s16(x[4]), vget_high_s16(d[20]));

    tmp[0] = vzip1q_s32(oo[0], oo[2]);
    tmp[1] = vzip2q_s32(oo[0], oo[2]);
    tmp[2] = vzip1q_s32(oo[4], oo[6]);
    tmp[3] = vzip2q_s32(oo[4], oo[6]);

    o[0] = (int32x4_t) vzip1q_s64((int64x2_t) tmp[0], (int64x2_t) tmp[2]);
    o[1] = (int32x4_t) vzip2q_s64((int64x2_t) tmp[0], (int64x2_t) tmp[2]);
    o[2] = (int32x4_t) vzip1q_s64((int64x2_t) tmp[1], (int64x2_t) tmp[3]);
    o[3] = (int32x4_t) vzip2q_s64((int64x2_t) tmp[1], (int64x2_t) tmp[3]);

    tmp[0] = vzip1q_s32(oo[1], oo[3]);
    tmp[1] = vzip2q_s32(oo[1], oo[3]);
    tmp[2] = vzip1q_s32(oo[5], oo[7]);
    tmp[3] = vzip2q_s32(oo[5], oo[7]);

    o[4] = (int32x4_t) vzip1q_s64((int64x2_t) tmp[0], (int64x2_t) tmp[2]);
    o[5] = (int32x4_t) vzip2q_s64((int64x2_t) tmp[0], (int64x2_t) tmp[2]);
    o[6] = (int32x4_t) vzip1q_s64((int64x2_t) tmp[1], (int64x2_t) tmp[3]);
    o[7] = (int32x4_t) vzip2q_s64((int64x2_t) tmp[1], (int64x2_t) tmp[3]);

    for(int k=0; k<8; k++){
        r[  k] = vaddq_s32(e[k], o[k]);
        int32x4_t tmp = vsubq_s32(e[k],  o[k]);
        r[8+k] = vcombine_s32(vrev64_s32(vget_high_s32(tmp)), vrev64_s32(vget_low_s32(tmp)));//_mm_shuffle_epi32(tmp, 0x1B)
    }
}

static inline void
dct2_16x8_red4(int16x8_t *x, int16x8_t *d, int32x4_t *r)
{
    int32x4_t e[16], o[16];
    int16x8_t *evn_x = x;
    int16x8_t *odd_x = x + 8;
    int16x8_t *odd_d = d + 24;

    dct2_8x8_red2(evn_x, d+0, e);

    for(int k=0; k<8; k++){
        int32x4_t m[4], a[2];
        m[0] = vmull_s16(vget_low_s16(odd_x[2*k]), vget_low_s16(odd_d[0]));
        m[1] = vmull_s16(vget_high_s16(odd_x[2*k]), vget_high_s16(odd_d[0]));
        m[2] = vmull_s16(vget_low_s16(odd_x[2*k + 1]), vget_low_s16(odd_d[1]));
        m[3] = vmull_s16(vget_high_s16(odd_x[2*k + 1]), vget_high_s16(odd_d[1]));

        a[0] = vaddq_s32(m[0], m[2]);
        a[1] = vaddq_s32(m[1], m[3]);

        o[2*k    ] = a[0];
        o[2*k + 1] = a[1];

        r[4*k+0] = vaddq_s32(e[k+0], o[2*k+0]);
        r[4*k+1] = vaddq_s32(e[k+8], o[2*k+1]);
        int32x4_t tmp = vsubq_s32(e[k+8], o[2*k+1]);
        r[4*k+2] = vcombine_s32(vrev64_s32(vget_high_s32(tmp)), vrev64_s32(vget_low_s32(tmp)));//_mm_shuffle_epi32(tmp, 0x1B);
        tmp = vsubq_s32(e[k+0], o[2*k+0]);
        r[4*k+3] = vcombine_s32(vrev64_s32(vget_high_s32(tmp)), vrev64_s32(vget_low_s32(tmp)));//_mm_shuffle_epi32(tmp, 0x1B);
    }
}

static inline void
dct2_32x8_red8(int16x8_t *x, int16x8_t *d, int32x4_t *r){
    int32x4_t e[32], o[32];
    int16x8_t *evn_x = x;
    int16x8_t *odd_x = x + 24;
    int16x8_t *odd_d1 = d + 32;
    int16x8_t *odd_d2 = d + 48;

    /*FIXME we can also use a reduced version of dct2_16x8*/
    dct2_16x8_red4(evn_x, d+0, e+0);

    for(int k=0; k<8; k++){
        matMult16x16_red4(odd_x + 4*k, odd_d1, o+4*k+0);
        matMult16x16_red4(odd_x + 4*k, odd_d2, o+4*k+2);
        r[8*k+0] = vaddq_s32(e[4*k+0], o[4*k+0]);
        r[8*k+1] = vaddq_s32(e[4*k+1], o[4*k+1]);
        r[8*k+2] = vaddq_s32(e[4*k+2], o[4*k+2]);
        r[8*k+3] = vaddq_s32(e[4*k+3], o[4*k+3]);
        int32x4_t tmp = vsubq_s32(e[4*k+3], o[4*k+3]);
        r[8*k+4] = vcombine_s32(vrev64_s32(vget_high_s32(tmp)), vrev64_s32(vget_low_s32(tmp)));//_mm_shuffle_epi32(tmp, 0x1B);
        tmp = vsubq_s32(e[4*k+2], o[4*k+2]);
        r[8*k+5] = vcombine_s32(vrev64_s32(vget_high_s32(tmp)), vrev64_s32(vget_low_s32(tmp)));//_mm_shuffle_epi32(tmp, 0x1B);
        tmp = vsubq_s32(e[4*k+1], o[4*k+1]);
        r[8*k+6] = vcombine_s32(vrev64_s32(vget_high_s32(tmp)), vrev64_s32(vget_low_s32(tmp)));//_mm_shuffle_epi32(tmp, 0x1B);
        tmp = vsubq_s32(e[4*k+0], o[4*k+0]);
        r[8*k+7] = vcombine_s32(vrev64_s32(vget_high_s32(tmp)), vrev64_s32(vget_low_s32(tmp)));//_mm_shuffle_epi32(tmp, 0x1B);
    }
}

static inline void
load_src_32x8_red8(const int16_t *src, ptrdiff_t src_stride, int16x8_t *x)
{
    int16x8_t *evn_xl = x;
    int16x8_t *odd_xl = x + 24;
    int k;

    evn_xl[0] = vld1q_s16(src +  0 * src_stride);
    evn_xl[4] = vld1q_s16(src +  4 * src_stride);

    for (k = 0; k < 8; ++k){
        const int16_t *src_col = src + k;
        evn_xl[ 8] = vdupq_n_s16(src_col[ 2 * src_stride]);
        evn_xl[ 9] = vdupq_n_s16(src_col[ 6 * src_stride]);
        evn_xl += 2;
    }

    for (k = 0; k < 8; ++k){
        const int16_t *src_col = src + k;
        odd_xl[ 0] = vdupq_n_s16(src_col[ 1 * src_stride]);
        odd_xl[ 1] = vdupq_n_s16(src_col[ 3 * src_stride]);
        odd_xl[ 2] = vdupq_n_s16(src_col[ 5 * src_stride]);
        odd_xl[ 3] = vdupq_n_s16(src_col[ 7 * src_stride]);
        odd_xl += 4;
    }
}

static inline void
dct2_4x8_red2(int16x8_t *x, int16x8_t *d, int32x4_t *r)
{
    int32x4_t a[8], tmp[4], m[8];

    a[0] = vmull_s16(vget_low_s16(x[0]), vget_low_s16(d[0]));
    a[1] = vmull_s16(vget_high_s16(x[0]), vget_high_s16(d[0]));
    a[2] = vmull_s16(vget_low_s16(x[0]), vget_low_s16(d[1]));
    a[3] = vmull_s16(vget_high_s16(x[0]), vget_high_s16(d[1]));

    a[4] = vmull_s16(vget_low_s16(x[1]), vget_low_s16(d[2]));
    a[5] = vmull_s16(vget_high_s16(x[1]), vget_high_s16(d[2]));
    a[6] = vmull_s16(vget_low_s16(x[1]), vget_low_s16(d[3]));
    a[7] = vmull_s16(vget_high_s16(x[1]), vget_high_s16(d[3]));

    m[0] = vaddq_s32(a[0], a[4]);
    m[1] = vaddq_s32(a[1], a[5]);
    m[2] = vaddq_s32(a[2], a[6]);
    m[3] = vaddq_s32(a[3], a[7]);
    m[4] = vsubq_s32(a[2], a[6]);
    m[5] = vsubq_s32(a[3], a[7]);
    m[6] = vsubq_s32(a[0], a[4]);
    m[7] = vsubq_s32(a[1], a[5]);

    tmp[0] = vzip1q_s32(m[0], m[2]);
    tmp[1] = vzip2q_s32(m[0], m[2]);
    tmp[2] = vzip1q_s32(m[4], m[6]);
    tmp[3] = vzip2q_s32(m[4], m[6]);

    r[0] = (int32x4_t) vzip1q_s64((int64x2_t) tmp[0], (int64x2_t) tmp[2]);
    r[1] = (int32x4_t) vzip2q_s64((int64x2_t) tmp[0], (int64x2_t) tmp[2]);
    r[2] = (int32x4_t) vzip1q_s64((int64x2_t) tmp[1], (int64x2_t) tmp[3]);
    r[3] = (int32x4_t) vzip2q_s64((int64x2_t) tmp[1], (int64x2_t) tmp[3]);

    tmp[0] = vzip1q_s32(m[1], m[3]);
    tmp[1] = vzip2q_s32(m[1], m[3]);
    tmp[2] = vzip1q_s32(m[5], m[7]);
    tmp[3] = vzip2q_s32(m[5], m[7]);

    r[4] = (int32x4_t) vzip1q_s64((int64x2_t) tmp[0], (int64x2_t) tmp[2]);
    r[5] = (int32x4_t) vzip2q_s64((int64x2_t) tmp[0], (int64x2_t) tmp[2]);
    r[6] = (int32x4_t) vzip1q_s64((int64x2_t) tmp[1], (int64x2_t) tmp[3]);
    r[7] = (int32x4_t) vzip2q_s64((int64x2_t) tmp[1], (int64x2_t) tmp[3]);
}

static inline void
dct2_8x8_red4(int16x8_t *x, int16x8_t *d, int32x4_t *r)
{
    int32x4_t e[8], o[8], oo[8];
    int32x4_t m[4];

    dct2_4x8_red2(x+0, d+0, e);

    int32x4_t tmp[4];

    m[0] = vmull_s16(vget_low_s16(x[4]), vget_low_s16(d[8]));
    m[1] = vmull_s16(vget_high_s16(x[4]), vget_high_s16(d[8]));
    m[2] = vmull_s16(vget_low_s16(x[5]), vget_low_s16(d[9]));
    m[3] = vmull_s16(vget_high_s16(x[5]), vget_high_s16(d[9]));

    oo[0] = vaddq_s32(m[0], m[2]);
    oo[1] = vaddq_s32(m[1], m[3]);

    m[0] = vmull_s16(vget_low_s16(x[4]), vget_low_s16(d[12]));
    m[1] = vmull_s16(vget_high_s16(x[4]), vget_high_s16(d[12]));
    m[2] = vmull_s16(vget_low_s16(x[5]), vget_low_s16(d[13]));
    m[3] = vmull_s16(vget_high_s16(x[5]), vget_high_s16(d[13]));

    oo[2] = vaddq_s32(m[0], m[2]);
    oo[3] = vaddq_s32(m[1], m[3]);

    m[0] = vmull_s16(vget_low_s16(x[4]), vget_low_s16(d[16]));
    m[1] = vmull_s16(vget_high_s16(x[4]), vget_high_s16(d[16]));
    m[2] = vmull_s16(vget_low_s16(x[5]), vget_low_s16(d[17]));
    m[3] = vmull_s16(vget_high_s16(x[5]), vget_high_s16(d[17]));

    oo[4] = vaddq_s32(m[0], m[2]);
    oo[5] = vaddq_s32(m[1], m[3]);

    m[0] = vmull_s16(vget_low_s16(x[4]), vget_low_s16(d[20]));
    m[1] = vmull_s16(vget_high_s16(x[4]), vget_high_s16(d[20]));
    m[2] = vmull_s16(vget_low_s16(x[5]), vget_low_s16(d[21]));
    m[3] = vmull_s16(vget_high_s16(x[5]), vget_high_s16(d[21]));

    oo[6] = vaddq_s32(m[0], m[2]);
    oo[7] = vaddq_s32(m[1], m[3]);

    tmp[0] = vzip1q_s32(oo[0], oo[2]);
    tmp[1] = vzip2q_s32(oo[0], oo[2]);
    tmp[2] = vzip1q_s32(oo[4], oo[6]);
    tmp[3] = vzip2q_s32(oo[4], oo[6]);

    o[0] = (int32x4_t) vzip1q_s64((int64x2_t) tmp[0], (int64x2_t) tmp[2]);
    o[1] = (int32x4_t) vzip2q_s64((int64x2_t) tmp[0], (int64x2_t) tmp[2]);
    o[2] = (int32x4_t) vzip1q_s64((int64x2_t) tmp[1], (int64x2_t) tmp[3]);
    o[3] = (int32x4_t) vzip2q_s64((int64x2_t) tmp[1], (int64x2_t) tmp[3]);

    tmp[0] = vzip1q_s32(oo[1], oo[3]);
    tmp[1] = vzip2q_s32(oo[1], oo[3]);
    tmp[2] = vzip1q_s32(oo[5], oo[7]);
    tmp[3] = vzip2q_s32(oo[5], oo[7]);

    o[4] = (int32x4_t) vzip1q_s64((int64x2_t) tmp[0], (int64x2_t) tmp[2]);
    o[5] = (int32x4_t) vzip2q_s64((int64x2_t) tmp[0], (int64x2_t) tmp[2]);
    o[6] = (int32x4_t) vzip1q_s64((int64x2_t) tmp[1], (int64x2_t) tmp[3]);
    o[7] = (int32x4_t) vzip2q_s64((int64x2_t) tmp[1], (int64x2_t) tmp[3]);

    for(int k=0; k<8; k++){
        r[  k] = vaddq_s32(e[k], o[k]);
        int32x4_t tmp = vsubq_s32(e[k],  o[k]);
        r[8+k] = vcombine_s32(vrev64_s32(vget_high_s32(tmp)), vrev64_s32(vget_low_s32(tmp)));//_mm_shuffle_epi32(tmp, 0x1B);
    }
}

static inline void
dct2_16x8_red8(int16x8_t *x, int16x8_t *d, int32x4_t *r)
{
    int32x4_t e[16], o[16];
    int16x8_t *evn_x = x;
    int16x8_t *odd_x = x + 8;
    int16x8_t *odd_d = d + 24;

    dct2_8x8_red4(evn_x, d+0, e);

    for(int k=0; k<8; k++){
        int32x4_t m[8], a[4];
        m[ 0] = vmull_s16(vget_low_s16(odd_x[4 * k]), vget_low_s16(odd_d[0]));
        m[ 1] = vmull_s16(vget_high_s16(odd_x[4 * k]), vget_high_s16(odd_d[0]));
        m[ 2] = vmull_s16(vget_low_s16(odd_x[4 * k + 1]), vget_low_s16(odd_d[1]));
        m[ 3] = vmull_s16(vget_high_s16(odd_x[4 * k + 1]), vget_high_s16(odd_d[1]));
        m[ 4] = vmull_s16(vget_low_s16(odd_x[4 * k + 2]), vget_low_s16(odd_d[2]));
        m[ 5] = vmull_s16(vget_high_s16(odd_x[4 * k + 2]), vget_high_s16(odd_d[2]));
        m[ 6] = vmull_s16(vget_low_s16(odd_x[4 * k + 3]), vget_low_s16(odd_d[3]));
        m[ 7] = vmull_s16(vget_high_s16(odd_x[4 * k + 3]), vget_high_s16(odd_d[3]));

        a[0] = vaddq_s32(m[ 0], m[ 2]);
        a[1] = vaddq_s32(m[ 1], m[ 3]);
        a[2] = vaddq_s32(m[ 4], m[ 6]);
        a[3] = vaddq_s32(m[ 5], m[ 7]);

        o[2*k    ] = vaddq_s32(a[0], a[2]);
        o[2*k + 1] = vaddq_s32(a[1], a[3]);

        r[4*k+0] = vaddq_s32(e[k+0], o[2*k+0]);
        r[4*k+1] = vaddq_s32(e[k+8], o[2*k+1]);
        int32x4_t tmp = vsubq_s32(e[k+8], o[2*k+1]);
        r[4*k+2] = vcombine_s32(vrev64_s32(vget_high_s32(tmp)), vrev64_s32(vget_low_s32(tmp)));//_mm_shuffle_epi32(tmp, 0x1B);
        tmp = vsubq_s32(e[k+0], o[2*k+0]);
        r[4*k+3] = vcombine_s32(vrev64_s32(vget_high_s32(tmp)), vrev64_s32(vget_low_s32(tmp)));//_mm_shuffle_epi32(tmp, 0x1B);
    }
}

static inline void
dct2_32x8_red16(int16x8_t *x, int16x8_t *d, int32x4_t *r){
    int32x4_t e[32], o[32];
    int16x8_t *evn_x = x;
    int16x8_t *odd_x = x + 40;
    int16x8_t *odd_d1 = d + 32;
    int16x8_t *odd_d2 = d + 48;

    dct2_16x8_red8(evn_x, d+0, e+0);

    for(int k=0; k<8; k++){
        matMult16x16_red8(odd_x + 8*k, odd_d1, o+4*k+0);
        matMult16x16_red8(odd_x + 8*k, odd_d2, o+4*k+2);

        r[8*k+0] = vaddq_s32(e[4*k+0], o[4*k+0]);
        r[8*k+1] = vaddq_s32(e[4*k+1], o[4*k+1]);
        r[8*k+2] = vaddq_s32(e[4*k+2], o[4*k+2]);
        r[8*k+3] = vaddq_s32(e[4*k+3], o[4*k+3]);
        int32x4_t tmp = vsubq_s32(e[4*k+3], o[4*k+3]);
        r[8*k+4] = vcombine_s32(vrev64_s32(vget_high_s32(tmp)), vrev64_s32(vget_low_s32(tmp)));//_mm_shuffle_epi32(tmp, 0x1B);
        tmp = vsubq_s32(e[4*k+2], o[4*k+2]);
        r[8*k+5] = vcombine_s32(vrev64_s32(vget_high_s32(tmp)), vrev64_s32(vget_low_s32(tmp)));//_mm_shuffle_epi32(tmp, 0x1B);
        tmp = vsubq_s32(e[4*k+1], o[4*k+1]);
        r[8*k+6] = vcombine_s32(vrev64_s32(vget_high_s32(tmp)), vrev64_s32(vget_low_s32(tmp)));//_mm_shuffle_epi32(tmp, 0x1B);
        tmp = vsubq_s32(e[4*k+0], o[4*k+0]);
        r[8*k+7] = vcombine_s32(vrev64_s32(vget_high_s32(tmp)), vrev64_s32(vget_low_s32(tmp)));//_mm_shuffle_epi32(tmp, 0x1B);
    }
}

static inline void
load_src_32x8_red16(const int16_t *src, ptrdiff_t src_stride, int16x8_t *x)
{
    int16x8_t *evn_xl = x;
    int16x8_t *odd_xl = x + 40;
    int k;

    evn_xl[0] = vld1q_s16(src +  0 * src_stride);
    evn_xl[1] = vld1q_s16(src +  8 * src_stride);
    evn_xl[2] = vdupq_n_s16(0);
    evn_xl[3] = vdupq_n_s16(0);
    evn_xl[4] = vld1q_s16(src +  4 * src_stride);
    evn_xl[5] = vld1q_s16(src + 12 * src_stride);
    evn_xl[6] = vdupq_n_s16(0);
    evn_xl[7] = vdupq_n_s16(0);

    for (k = 0; k < 8; ++k){
        const int16_t *src_col = src + k;
        evn_xl[ 8] = vdupq_n_s16(src_col[ 2 * src_stride]);
        evn_xl[ 9] = vdupq_n_s16(src_col[ 6 * src_stride]);
        evn_xl[10] = vdupq_n_s16(src_col[10 * src_stride]);
        evn_xl[11] = vdupq_n_s16(src_col[14 * src_stride]);
        evn_xl += 4;
    }

    for (k = 0; k < 8; ++k){
        const int16_t *src_col = src + k;
        odd_xl[ 0] = vdupq_n_s16(src_col[ 1 * src_stride]);
        odd_xl[ 1] = vdupq_n_s16(src_col[ 3 * src_stride]);
        odd_xl[ 2] = vdupq_n_s16(src_col[ 5 * src_stride]);
        odd_xl[ 3] = vdupq_n_s16(src_col[ 7 * src_stride]);
        odd_xl[ 4] = vdupq_n_s16(src_col[ 9 * src_stride]);
        odd_xl[ 5] = vdupq_n_s16(src_col[11 * src_stride]);
        odd_xl[ 6] = vdupq_n_s16(src_col[13 * src_stride]);
        odd_xl[ 7] = vdupq_n_s16(src_col[15 * src_stride]);
        odd_xl += 8;
    }
}

static inline void
load_src_32x8(const int16_t *src, ptrdiff_t src_stride, int16x8_t *x)
{
    int16x8_t *evn_xl = x;
    int16x8_t *odd_xl = x + 72;
    int k;

    evn_xl[0] = vld1q_s16(src +  0 * src_stride);
    evn_xl[1] = vld1q_s16(src +  8 * src_stride);
    evn_xl[2] = vld1q_s16(src + 16 * src_stride);
    evn_xl[3] = vld1q_s16(src + 24 * src_stride);
    evn_xl[4] = vld1q_s16(src +  4 * src_stride);
    evn_xl[5] = vld1q_s16(src + 12 * src_stride);
    evn_xl[6] = vld1q_s16(src + 20 * src_stride);
    evn_xl[7] = vld1q_s16(src + 28 * src_stride);

    for (k = 0; k < 8; ++k){
        const int16_t *src_col = src + k;
        evn_xl[ 8] = vdupq_n_s16(src_col[ 2 * src_stride]);
        evn_xl[ 9] = vdupq_n_s16(src_col[ 6 * src_stride]);
        evn_xl[10] = vdupq_n_s16(src_col[10 * src_stride]);
        evn_xl[11] = vdupq_n_s16(src_col[14 * src_stride]);
        evn_xl[12] = vdupq_n_s16(src_col[18 * src_stride]);
        evn_xl[13] = vdupq_n_s16(src_col[22 * src_stride]);
        evn_xl[14] = vdupq_n_s16(src_col[26 * src_stride]);
        evn_xl[15] = vdupq_n_s16(src_col[30 * src_stride]);
        evn_xl += 8;
    }

    for (k = 0; k < 8; ++k){
        const int16_t *src_col = src + k;
        odd_xl[ 0] = vdupq_n_s16(src_col[ 1 * src_stride]);
        odd_xl[ 1] = vdupq_n_s16(src_col[ 3 * src_stride]);
        odd_xl[ 2] = vdupq_n_s16(src_col[ 5 * src_stride]);
        odd_xl[ 3] = vdupq_n_s16(src_col[ 7 * src_stride]);
        odd_xl[ 4] = vdupq_n_s16(src_col[ 9 * src_stride]);
        odd_xl[ 5] = vdupq_n_s16(src_col[11 * src_stride]);
        odd_xl[ 6] = vdupq_n_s16(src_col[13 * src_stride]);
        odd_xl[ 7] = vdupq_n_s16(src_col[15 * src_stride]);
        odd_xl[ 8] = vdupq_n_s16(src_col[17 * src_stride]);
        odd_xl[ 9] = vdupq_n_s16(src_col[19 * src_stride]);
        odd_xl[10] = vdupq_n_s16(src_col[21 * src_stride]);
        odd_xl[11] = vdupq_n_s16(src_col[23 * src_stride]);
        odd_xl[12] = vdupq_n_s16(src_col[25 * src_stride]);
        odd_xl[13] = vdupq_n_s16(src_col[27 * src_stride]);
        odd_xl[14] = vdupq_n_s16(src_col[29 * src_stride]);
        odd_xl[15] = vdupq_n_s16(src_col[31 * src_stride]);
        odd_xl += 16;
    }
}

static inline void
load_dct_32x8(int16x8_t *d)
{
    int k;
    d[0] = vdupq_n_s16(DCT_II_32[0]);
    d[1] = vdupq_n_s16(DCT_II_32[1]);
    d[2] = vdupq_n_s16(DCT_II_32[256]);
    d[3] = vdupq_n_s16(DCT_II_32[257]);
    d[4] = vdupq_n_s16(DCT_II_32[512]);
    d[5] = vdupq_n_s16(DCT_II_32[513]);
    d[6] = vdupq_n_s16(DCT_II_32[768]);
    d[7] = vdupq_n_s16(DCT_II_32[769]);

    for (k = 0; k < 4; k++) {
        d[8 + 4 * k + 0] = vdupq_n_s16(DCT_II_32[128 + k]);
        d[8 + 4 * k + 1] = vdupq_n_s16(DCT_II_32[384 + k]);
        d[8 + 4 * k + 2] = vdupq_n_s16(DCT_II_32[640 + k]);
        d[8 + 4 * k + 3] = vdupq_n_s16(DCT_II_32[896 + k]);
    }

    static const int16_t DCT_II_32_8_neon[8 * 40] = {
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
        d[24+k] = vld1q_s16(DCT_II_32_8_neon + 8*k);
    }

}

static void
dct2_32_8lines_red8(const int16_t *src, int16_t *dst, ptrdiff_t src_stride,
               int num_lines, int line_brk, int shift)
{
    int rnd_add = 1 << (shift - 1);
    int16x8_t d[64];
    /*FIXME we can still reduce dct coefficient loading operation
      however it is less intensive than src_load since it is done
      outside of loop */
    load_dct_32x8(d);
    int32x4_t rnd_add_v = vdupq_n_s32(rnd_add);
    for (int j = 0; j < num_lines / 8; j++) {
        int16x8_t x[56];
        int32x4_t r[64];

        load_src_32x8_red8(src, src_stride, x);

        dct2_32x8_red8(x, d, r);


        for (int i = 0; i < 64; i+=2) {
            int16x8_t o;

            r[i + 0] = vaddq_s32(r[i + 0], rnd_add_v);
            r[i + 1] = vaddq_s32(r[i + 1], rnd_add_v);

            r[i + 0] = vshrq_n_s32(r[i + 0], shift);
            r[i + 1] = vshrq_n_s32(r[i + 1], shift);

            o = vcombine_s16(vqmovn_s32(r[i + 0]), vqmovn_s32(r[i + 1]));

            vst1q_s16(dst + i * 8/2, o);
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
    int16x8_t d[64];
    load_dct_32x8(d);
    int32x4_t rnd_add_v = vdupq_n_s32(rnd_add);
    for (int j = 0; j < num_lines / 8; j++) {
        int16x8_t x[104];
        int32x4_t r[64];

        load_src_32x8_red16(src, src_stride, x);

        dct2_32x8_red16(x, d, r);


        for (int i = 0; i < 64; i+=2) {
            int16x8_t o;

            r[i + 0] = vaddq_s32(r[i + 0], rnd_add_v);
            r[i + 1] = vaddq_s32(r[i + 1], rnd_add_v);

            r[i + 0] = vshrq_n_s32(r[i + 0], shift);
            r[i + 1] = vshrq_n_s32(r[i + 1], shift);

            o = vcombine_s16(vqmovn_s32(r[i + 0]), vqmovn_s32(r[i + 1]));

            vst1q_s16(dst + i * 8/2, o);
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
    int16x8_t d[64];
    load_dct_32x8(d);
    int32x4_t rnd_add_v = vdupq_n_s32(rnd_add);
    for (int j = 0; j < num_lines / 8; j++) {
        int16x8_t x[200];
        int32x4_t r[64];

        load_src_32x8(src, src_stride, x);

        dct2_32x8(x, d, r);


        for (int i = 0; i < 64; i+=2) {
            int16x8_t o;

            r[i + 0] = vaddq_s32(r[i + 0], rnd_add_v);
            r[i + 1] = vaddq_s32(r[i + 1], rnd_add_v);

            r[i + 0] = vshrq_n_s32(r[i + 0], shift);
            r[i + 1] = vshrq_n_s32(r[i + 1], shift);

            o = vcombine_s16(vqmovn_s32(r[i + 0]), vqmovn_s32(r[i + 1]));

            vst1q_s16(dst + i * 8/2, o);
        }
        src += 8;
        dst += 256;
    }
}
static inline void
load_dct_32x4(int16x8_t *d)
{
    static const int16_t DCT_II_32_4_neon[8 * 48] = {
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
        d[k] = vld1q_s16(DCT_II_32_4_neon + 8*k);
    }
}

static inline void
load_src_32x4_red4(const int16_t *src, ptrdiff_t src_stride, int16x8_t *x)
{
    int16x8_t *evn_xl = x;
    int16x8_t *odd_xl = x + 74;
    int k;

    x[0] = vcombine_s16(
            vld1_s16(src + 0 * src_stride),
            vdup_n_s16(0)
            );

    for (k = 0; k < 8; ++k){
        const int16_t *src_col = src + k;
        evn_xl[10] = vdupq_n_s16(src_col[ 2 * src_stride]);
        evn_xl += 8;
    }

    for (k = 0; k < 8; ++k){
        const int16_t *src_col = src + k;
        odd_xl[ 0] = vdupq_n_s16(src_col[ 1 * src_stride]);
        odd_xl[ 1] = vdupq_n_s16(src_col[ 3 * src_stride]);
        odd_xl += 16;
    }
}

static inline void
load_src_32x4_red8(const int16_t *src, ptrdiff_t src_stride, int16x8_t *x)
{
    int16x8_t *evn_xl = x;
    int16x8_t *odd_xl = x + 74;
    int k;

    x[0] = vcombine_s16(
            vld1_s16(src + 0 * src_stride),
            vdup_n_s16(0)
            );
    x[1] = vdupq_n_s16(0);
    x[2] = vcombine_s16(
            vdup_n_s16(src[4 * src_stride + 0]),
            vdup_n_s16(src[4 * src_stride + 1])
            );
    x[3] = vdupq_n_s16(0);
    x[4] = vdupq_n_s16(0);
    x[5] = vdupq_n_s16(0);
    x[6] = vcombine_s16(
            vdup_n_s16(src[4 * src_stride + 2]),
            vdup_n_s16(src[4 * src_stride + 3])
            );
    x[7] = vdupq_n_s16(0);
    x[8] = vdupq_n_s16(0);
    x[9] = vdupq_n_s16(0);

    for (k = 0; k < 8; ++k){
        const int16_t *src_col = src + k;
        evn_xl[10] = vdupq_n_s16(src_col[ 2 * src_stride]);
        evn_xl[11] = vdupq_n_s16(src_col[ 6 * src_stride]);
        evn_xl += 8;
    }

    for (k = 0; k < 8; ++k){
        const int16_t *src_col = src + k;
        odd_xl[ 0] = vdupq_n_s16(src_col[ 1 * src_stride]);
        odd_xl[ 1] = vdupq_n_s16(src_col[ 3 * src_stride]);
        odd_xl[ 2] = vdupq_n_s16(src_col[ 5 * src_stride]);
        odd_xl[ 3] = vdupq_n_s16(src_col[ 7 * src_stride]);
        odd_xl += 16;
    }
}

static inline void
load_src_32x4_red16(const int16_t *src, ptrdiff_t src_stride, int16x8_t *x)
{
    int16x8_t *evn_xl = x;
    int16x8_t *odd_xl = x + 74;
    int k;

    x[0] = vcombine_s16(
            vld1_s16(src + 0 * src_stride),
            vld1_s16(src + 8 * src_stride)
            );
    x[1] = vcombine_s16(
            vld1_s16(src + 16 * src_stride),
            vdup_n_s16(0)
            );
    x[2] = vcombine_s16(
            vdup_n_s16(src[4 * src_stride + 0]),
            vdup_n_s16(src[4 * src_stride + 1])
            );
    x[3] = vcombine_s16(
            vdup_n_s16(src[12 * src_stride + 0]),
            vdup_n_s16(src[12 * src_stride + 1])
            );
    x[4] = vdupq_n_s16(0);
    x[5] = vdupq_n_s16(0);
    x[6] = vcombine_s16(
            vdup_n_s16(src[4 * src_stride + 2]),
            vdup_n_s16(src[4 * src_stride + 3])
            );
    x[7] = vcombine_s16(
            vdup_n_s16(src[12 * src_stride + 2]),
            vdup_n_s16(src[12 * src_stride + 3])
            );
    x[8] = vdupq_n_s16(0);
    x[9] = vdupq_n_s16(0);

    for (k = 0; k < 8; ++k){
        const int16_t *src_col = src + k;
        evn_xl[10] = vdupq_n_s16(src_col[ 2 * src_stride]);
        evn_xl[11] = vdupq_n_s16(src_col[ 6 * src_stride]);
        evn_xl[12] = vdupq_n_s16(src_col[10 * src_stride]);
        evn_xl[13] = vdupq_n_s16(src_col[14 * src_stride]);
        evn_xl += 8;
    }

    for (k = 0; k < 8; ++k){
        const int16_t *src_col = src + k;
        odd_xl[ 0] = vdupq_n_s16(src_col[ 1 * src_stride]);
        odd_xl[ 1] = vdupq_n_s16(src_col[ 3 * src_stride]);
        odd_xl[ 2] = vdupq_n_s16(src_col[ 5 * src_stride]);
        odd_xl[ 3] = vdupq_n_s16(src_col[ 7 * src_stride]);
        odd_xl[ 4] = vdupq_n_s16(src_col[ 9 * src_stride]);
        odd_xl[ 5] = vdupq_n_s16(src_col[11 * src_stride]);
        odd_xl[ 6] = vdupq_n_s16(src_col[13 * src_stride]);
        odd_xl[ 7] = vdupq_n_s16(src_col[15 * src_stride]);
        odd_xl += 16;
    }
}

static inline void
load_src_32x4(const int16_t *src, ptrdiff_t src_stride, int16x8_t *x)
{
    int16x8_t *evn_xl = x;
    int16x8_t *odd_xl = x + 74;
    int k;

    x[0] = vcombine_s16(
            vld1_s16(src + 0 * src_stride),
            vld1_s16(src + 8 * src_stride)
            );
    x[1] = vcombine_s16(
            vld1_s16(src + 16 * src_stride),
            vld1_s16(src + 24 * src_stride)
            );
    x[2] = vcombine_s16(
            vdup_n_s16(src[4 * src_stride + 0]),
            vdup_n_s16(src[4 * src_stride + 1])
            );
    x[3] = vcombine_s16(
            vdup_n_s16(src[12 * src_stride + 0]),
            vdup_n_s16(src[12 * src_stride + 1])
            );
    x[4] = vcombine_s16(
            vdup_n_s16(src[20 * src_stride + 0]),
            vdup_n_s16(src[20 * src_stride + 1])
            );
    x[5] = vcombine_s16(
            vdup_n_s16(src[28 * src_stride + 0]),
            vdup_n_s16(src[28 * src_stride + 1])
            );
    x[6] = vcombine_s16(
            vdup_n_s16(src[4 * src_stride + 2]),
            vdup_n_s16(src[4 * src_stride + 3])
            );
    x[7] = vcombine_s16(
            vdup_n_s16(src[12 * src_stride + 2]),
            vdup_n_s16(src[12 * src_stride + 3])
            );
    x[8] = vcombine_s16(
            vdup_n_s16(src[20 * src_stride + 2]),
            vdup_n_s16(src[20 * src_stride + 3])
            );
    x[9] = vcombine_s16(
            vdup_n_s16(src[28 * src_stride + 2]),
            vdup_n_s16(src[28 * src_stride + 3])
            );

    for (k = 0; k < 8; ++k){
        const int16_t *src_col = src + k;
        evn_xl[10] = vdupq_n_s16(src_col[ 2 * src_stride]);
        evn_xl[11] = vdupq_n_s16(src_col[ 6 * src_stride]);
        evn_xl[12] = vdupq_n_s16(src_col[10 * src_stride]);
        evn_xl[13] = vdupq_n_s16(src_col[14 * src_stride]);
        evn_xl[14] = vdupq_n_s16(src_col[18 * src_stride]);
        evn_xl[15] = vdupq_n_s16(src_col[22 * src_stride]);
        evn_xl[16] = vdupq_n_s16(src_col[26 * src_stride]);
        evn_xl[17] = vdupq_n_s16(src_col[30 * src_stride]);
        evn_xl += 8;
    }

    for (k = 0; k < 8; ++k){
        const int16_t *src_col = src + k;
        odd_xl[ 0] = vdupq_n_s16(src_col[ 1 * src_stride]);
        odd_xl[ 1] = vdupq_n_s16(src_col[ 3 * src_stride]);
        odd_xl[ 2] = vdupq_n_s16(src_col[ 5 * src_stride]);
        odd_xl[ 3] = vdupq_n_s16(src_col[ 7 * src_stride]);
        odd_xl[ 4] = vdupq_n_s16(src_col[ 9 * src_stride]);
        odd_xl[ 5] = vdupq_n_s16(src_col[11 * src_stride]);
        odd_xl[ 6] = vdupq_n_s16(src_col[13 * src_stride]);
        odd_xl[ 7] = vdupq_n_s16(src_col[15 * src_stride]);
        odd_xl[ 8] = vdupq_n_s16(src_col[17 * src_stride]);
        odd_xl[ 9] = vdupq_n_s16(src_col[19 * src_stride]);
        odd_xl[10] = vdupq_n_s16(src_col[21 * src_stride]);
        odd_xl[11] = vdupq_n_s16(src_col[23 * src_stride]);
        odd_xl[12] = vdupq_n_s16(src_col[25 * src_stride]);
        odd_xl[13] = vdupq_n_s16(src_col[27 * src_stride]);
        odd_xl[14] = vdupq_n_s16(src_col[29 * src_stride]);
        odd_xl[15] = vdupq_n_s16(src_col[31 * src_stride]);
        odd_xl += 16;
    }
}

static void
dct2_32_4lines_red4(const int16_t *src, int16_t *dst, ptrdiff_t src_stride,
               int num_lines, int line_brk, int shift)
{
    int16x8_t x[202], d[48];
    int32x4_t r[32];

    load_src_32x4_red4(src, src_stride, x);

    load_dct_32x4(d);

    dct2_32x4_red4(x, d, r);

    int32x4_t add = vdupq_n_s32(1 << (shift - 1));
    for (int i = 0; i < 32; i += 2) {
        int16x8_t o;
        r[i + 0] = vaddq_s32(r[i + 0], add);
        r[i + 1] = vaddq_s32(r[i + 1], add);

        r[i + 0] = vshrq_n_s32(r[i + 0], shift);
        r[i + 1] = vshrq_n_s32(r[i + 1], shift);

        o = vcombine_s16(vqmovn_s32(r[i + 0]), vqmovn_s32(r[i + 1]));

        vst1q_s16(dst + i / 2 * 8, o);
    }
}

static void
dct2_32_4lines_red8(const int16_t *src, int16_t *dst, ptrdiff_t src_stride,
               int num_lines, int line_brk, int shift)
{
    int16x8_t x[202], d[48];
    int32x4_t r[32];

    load_src_32x4_red8(src, src_stride, x);

    load_dct_32x4(d);

    dct2_32x4_red8(x, d, r);

    int32x4_t add = vdupq_n_s32(1 << (shift - 1));
    for (int i = 0; i < 32; i += 2) {
        int16x8_t o;
        r[i + 0] = vaddq_s32(r[i + 0], add);
        r[i + 1] = vaddq_s32(r[i + 1], add);

        r[i + 0] = vshrq_n_s32(r[i + 0], shift);
        r[i + 1] = vshrq_n_s32(r[i + 1], shift);

        o = vcombine_s16(vqmovn_s32(r[i + 0]), vqmovn_s32(r[i + 1]));

        vst1q_s16(dst + i / 2 * 8, o);
    }
}

static void
dct2_32_4lines_red16(const int16_t *src, int16_t *dst, ptrdiff_t src_stride,
               int num_lines, int line_brk, int shift)
{
    int16x8_t x[202], d[48];
    int32x4_t r[32];

    load_src_32x4_red16(src, src_stride, x);

    load_dct_32x4(d);

    dct2_32x4_red16(x, d, r);

    int32x4_t add = vdupq_n_s32(1 << (shift - 1));
    for (int i = 0; i < 32; i += 2) {
        int16x8_t o;
        r[i + 0] = vaddq_s32(r[i + 0], add);
        r[i + 1] = vaddq_s32(r[i + 1], add);

        r[i + 0] = vshrq_n_s32(r[i + 0], shift);
        r[i + 1] = vshrq_n_s32(r[i + 1], shift);

        o = vcombine_s16(vqmovn_s32(r[i + 0]), vqmovn_s32(r[i + 1]));

        vst1q_s16(dst + i / 2 * 8, o);
    }
}

static void
dct2_32_4lines(const int16_t *src, int16_t *dst, ptrdiff_t src_stride,
               int num_lines, int line_brk, int shift)
{
    int16x8_t x[202], d[48];
    int32x4_t r[32];

    load_src_32x4(src, src_stride, x);

    load_dct_32x4(d);

    dct2_32x4(x, d, r);

    int32x4_t add = vdupq_n_s32(1 << (shift - 1));
    for (int i = 0; i < 32; i += 2) {
        int16x8_t o;
        r[i + 0] = vaddq_s32(r[i + 0], add);
        r[i + 1] = vaddq_s32(r[i + 1], add);

        r[i + 0] = vshrq_n_s32(r[i + 0], shift);
        r[i + 1] = vshrq_n_s32(r[i + 1], shift);

        o = vcombine_s16(vqmovn_s32(r[i + 0]), vqmovn_s32(r[i + 1]));

        vst1q_s16(dst + i / 2 * 8, o);
    }
}

void
vvc_inverse_dct_ii_32_neon(const int16_t *src, int16_t *dst, ptrdiff_t src_stride,
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
        int16x8_t x[56], d[46];
        int32x4_t r[16];

        x[0] = (int16x8_t) vzip1q_s32(
                	(int32x4_t)vcombine_s16(vld1_s16(src + 0 * src_stride), vdup_n_s16(0)),
                	(int32x4_t)vcombine_s16(vld1_s16(src + 8 * src_stride), vdup_n_s16(0))
                );
        x[0] = vzip1q_s16(x[0], x[0]);
        x[1] = (int16x8_t) vzip1q_s32(
                	(int32x4_t)vcombine_s16(vld1_s16(src + 16 * src_stride), vdup_n_s16(0)),
                	(int32x4_t)vcombine_s16(vld1_s16(src + 24 * src_stride), vdup_n_s16(0))
                );
        x[1] = vzip1q_s16(x[1], x[1]);

        x[2] = vcombine_s16(
                vdup_n_s16(src[ 4 * src_stride    ]),
                vdup_n_s16(src[ 4 * src_stride + 1])
                );
        x[3] = vcombine_s16(
                vdup_n_s16(src[12 * src_stride    ]),
                vdup_n_s16(src[12 * src_stride + 1])
                );
        x[4] = vcombine_s16(
                vdup_n_s16(src[20 * src_stride    ]),
                vdup_n_s16(src[20 * src_stride + 1])
                );
        x[5] = vcombine_s16(
                vdup_n_s16(src[28 * src_stride    ]),
                vdup_n_s16(src[28 * src_stride + 1])
                );

        for (int k=0; k<2; k++){
            x[6+8*k+0] = vdupq_n_s16(src[ 2 * src_stride + k]);
            x[6+8*k+1] = vdupq_n_s16(src[ 6 * src_stride + k]);
            x[6+8*k+2] = vdupq_n_s16(src[10 * src_stride + k]);
            x[6+8*k+3] = vdupq_n_s16(src[14 * src_stride + k]);
            x[6+8*k+4] = vdupq_n_s16(src[18 * src_stride + k]);
            x[6+8*k+5] = vdupq_n_s16(src[22 * src_stride + k]);
            x[6+8*k+6] = vdupq_n_s16(src[26 * src_stride + k]);
            x[6+8*k+7] = vdupq_n_s16(src[30 * src_stride + k]);
        }

        for (int k=0; k<2; k++) {
            for (int l = 0; l < 16; l++) {
                x[24+16*k+l] = vdupq_n_s16(src[(2*l+1) * src_stride + k]);
            }
        }

        static const int16_t DCT_II_32_2_neon[8 * 46] = {
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
            d[k] = vld1q_s16(DCT_II_32_2_neon + 8*k);
        }

        dct2_32x2(x, d, r);

        int32x4_t add = vdupq_n_s32(1 << (shift - 1));

        for (int i = 0; i < 16; i += 2) {
            int16x8_t o;

            r[i + 0] = vaddq_s32(r[i + 0], add);
            r[i + 1] = vaddq_s32(r[i + 1], add);

            r[i + 0] = vshrq_n_s32(r[i + 0], shift);
            r[i + 1] = vshrq_n_s32(r[i + 1], shift);

            o = vcombine_s16(vqmovn_s32(r[i + 0]), vqmovn_s32(r[i + 1]));

            vst1q_s16(dst + i / 2 * 8, o);
        }
    }

    if (num_lines & 0x1){
      vvc_inverse_dct_ii_32(src, dst, src_stride, num_lines & 0x1, line_brk, shift);
    }
}
void
vvc_inverse_dct_ii_64_neon_8lines_red8(const int16_t* src, int16_t* dst, ptrdiff_t src_stride,
                      int num_lines, int num_columns, int shift)
{
  int32x4_t add = vdupq_n_s32(1 << (shift - 1));

  int j, k;

  for (j = 0; j < num_lines>>3; j++) {
        /* Utilizing symmetry properties to the maximum to minimize
        the
         * number of multiplications */
         int16x8_t x[4], d[4];
         int32x4_t m[4], a[4], r[128], r2[128], r3[128];

         //EEE
         uint16_t DCT_II_64_EEE[8] ={64,  64,  64,  64, 64,  64,  64,  64};
         x[0 ]=vld1q_s16(src);
         for (k = 0; k < 8; k++) {
           d[ 0] = vdupq_n_s16(DCT_II_64_EEE[k]);

           r[k*2] = vmull_s16(vget_low_s16(x[0]), vget_low_s16(d[0]));
           r[k*2+1] = vmull_s16(vget_high_s16(x[0]), vget_high_s16(d[0]));
         }

         //EEO
         x[0 ]=vld1q_s16(src +  4 * src_stride);
         for (k = 0; k < 8; k++) {
           d[ 0] = vdupq_n_s16(DCT_II_64_EEOT[k * 4 + 0]);

           r2[k*2] = vmull_s16(vget_low_s16(x[0]), vget_low_s16(d[0]));
           r2[k*2+1] = vmull_s16(vget_high_s16(x[0]), vget_high_s16(d[0]));
         }
         //EE
         /* Combining even and odd terms at each hierarchy levels to
          * calculate the final spatial domain vector */
         for (k = 0; k < 8; k++) {
           r3[k*2] = vaddq_s32(r2[k*2], r[k*2]);
           r3[k*2+1] = vaddq_s32(r2[k*2+1], r[k*2+1]);

           r3[k*2+16] = vsubq_s32(r[14-k*2], r2[14-k*2]);
           r3[k*2+17] = vsubq_s32(r[15-k*2], r2[15-k*2]);
         }

/**/
         //EO
         x[0]=vld1q_s16(src +  2 * src_stride);
         x[1]=vld1q_s16(src +  6 * src_stride);
         for (k = 0; k < 16; k++) {
           d[0] = vdupq_n_s16(DCT_II_64_EOT[k * 8 + 0]);
           d[1] = vdupq_n_s16(DCT_II_64_EOT[k * 8 + 1]);

           m[0] = vmull_s16(vget_low_s16(x[0]), vget_low_s16(d[0]));
           m[1] = vmull_s16(vget_high_s16(x[0]), vget_high_s16(d[0]));
           m[2] = vmull_s16(vget_low_s16(x[1]), vget_low_s16(d[1]));
           m[3] = vmull_s16(vget_high_s16(x[1]), vget_high_s16(d[1]));

           r[k*2] = vaddq_s32(m[0], m[2]);
           r[k*2+1] = vaddq_s32(m[1], m[3]);
         }

         //E
         /* Combining even and odd terms at each hierarchy levels to
          * calculate the final spatial domain vector */
         for (k = 0; k < 16; k++) {
           r2[k*2] = vaddq_s32(r[k*2], r3[k*2]);
           r2[k*2+1] = vaddq_s32(r[k*2+1], r3[k*2+1]);

           r2[k*2+32] = vsubq_s32(r3[30-k*2], r[30-k*2]);
           r2[k*2+33] = vsubq_s32(r3[31-k*2], r[31-k*2]);
         }
         //O
         x[0]=vld1q_s16(src +  src_stride);
         x[1]=vld1q_s16(src +  3  * src_stride);
         x[2]=vld1q_s16(src +  5  * src_stride);
         x[3]=vld1q_s16(src +  7  * src_stride);
         for (k = 0; k < 32; k++) {
           d[0] = vdupq_n_s16(DCT_II_64_OT[k * 16 + 0 ]);
           d[1] = vdupq_n_s16(DCT_II_64_OT[k * 16 + 1 ]);
           d[2] = vdupq_n_s16(DCT_II_64_OT[k * 16 + 2 ]);
           d[3] = vdupq_n_s16(DCT_II_64_OT[k * 16 + 3 ]);

           // m[0] = vzip1q_s16(x[0],  x[1]);
           // m[1] = vzip1q_s16(x[2],  x[3]);
           //
           // m[2] = vzip2q_s16(x[0],  x[1]);
           // m[3] = vzip2q_s16(x[2],  x[3]);
           //
           // di[0] = vzip1q_s16(d[0],  d[1]);
           // di[1] = vzip1q_s16(d[2],  d[3]);
           //
           // a[0] = _mm_madd_epi16(m[0], di[0]);
           // a[1] = _mm_madd_epi16(m[1], di[1]);
           //
           // a[2] = _mm_madd_epi16(m[2], di[0]);
           // a[3] = _mm_madd_epi16(m[3], di[1]);

           m[0] = vmull_s16(vget_low_s16(x[0]), vget_low_s16(d[0]));
           m[1] = vmull_s16(vget_high_s16(x[0]), vget_high_s16(d[0]));
           m[2] = vmull_s16(vget_low_s16(x[1]), vget_low_s16(d[1]));
           m[3] = vmull_s16(vget_high_s16(x[1]), vget_high_s16(d[1]));
           m[4] = vmull_s16(vget_low_s16(x[2]), vget_low_s16(d[2]));
           m[5] = vmull_s16(vget_high_s16(x[2]), vget_high_s16(d[2]));
           m[6] = vmull_s16(vget_low_s16(x[3]), vget_low_s16(d[3]));
           m[7] = vmull_s16(vget_high_s16(x[3]), vget_high_s16(d[3]));

           a[0] = vaddq_s32(m[0], m[2]);
           a[1] = vaddq_s32(m[4], m[6]);

           a[2] = vaddq_s32(m[1], m[3]);
           a[3] = vaddq_s32(m[5], m[7]);

           r3[k*2] = vaddq_s32(a[0], a[1]);

           r3[k*2+1] = vaddq_s32(a[2], a[3]);
         }

         //Result
         for (k = 0; k < 32; k++) {
           r[k*2] =   vaddq_s32(r3[k*2], r2[k*2]);
           r[k*2] =   vaddq_s32(r[k*2],add);
           r[k*2] =   vshrq_n_s32(r[k*2],shift);
           r[k*2+1] = vaddq_s32(r3[k*2+1], r2[k*2+1]);
           r[k*2+1] = vaddq_s32(r[k*2+1],add);
           r[k*2+1] = vshrq_n_s32(r[k*2+1],shift);
           r[k*2] = (int32x4_t) vcombine_s16(vqmovn_s32(r[k*2]), vqmovn_s32(r[k*2+1]));

           r[k*2+64] = vsubq_s32(r2[62-k*2], r3[62-k*2]);
           r[k*2+64] = vaddq_s32(r[k*2+64],add);
           r[k*2+64] = vshrq_n_s32(r[k*2+64],shift);
           r[k*2+65] = vsubq_s32(r2[63-k*2], r3[63-k*2]);
           r[k*2+65] = vaddq_s32(r[k*2+65],add);
           r[k*2+65] = vshrq_n_s32(r[k*2+65],shift);
           r[k*2+64] = (int32x4_t) vcombine_s16(vqmovn_s32(r[k*2+64]), vqmovn_s32(r[k*2+65]));
         }

        for (k = 0; k < 8; k++) {
              int16x8_t tmp[8], tmp2[8];
              tmp[0] = vzip1q_s16((int16x8_t)r[k*16 + 0], (int16x8_t) r[k*16 + 2]);
              tmp[1] = vzip2q_s16((int16x8_t)r[k*16 + 0], (int16x8_t) r[k*16 + 2]);
              tmp[2] = vzip1q_s16((int16x8_t)r[k*16 + 4], (int16x8_t) r[k*16 + 6]);
              tmp[3] = vzip2q_s16((int16x8_t)r[k*16 + 4], (int16x8_t) r[k*16 + 6]);
              tmp[4] = vzip1q_s16((int16x8_t)r[k*16 + 8], (int16x8_t) r[k*16 + 10]);
              tmp[5] = vzip2q_s16((int16x8_t)r[k*16 + 8], (int16x8_t) r[k*16 + 10]);
              tmp[6] = vzip1q_s16((int16x8_t)r[k*16 + 12], (int16x8_t) r[k*16 + 14]);
              tmp[7] = vzip2q_s16((int16x8_t)r[k*16 + 12], (int16x8_t) r[k*16 + 14]);

              tmp2[0] = (int16x8_t) vzip1q_s32((int32x4_t)tmp[0], (int32x4_t)tmp[2]);
              tmp2[1] = (int16x8_t) vzip2q_s32((int32x4_t)tmp[0], (int32x4_t)tmp[2]);
              tmp2[2] = (int16x8_t) vzip1q_s32((int32x4_t)tmp[1], (int32x4_t)tmp[3]);
              tmp2[3] = (int16x8_t) vzip2q_s32((int32x4_t)tmp[1], (int32x4_t)tmp[3]);
              tmp2[4] = (int16x8_t) vzip1q_s32((int32x4_t)tmp[4], (int32x4_t)tmp[6]);
              tmp2[5] = (int16x8_t) vzip2q_s32((int32x4_t)tmp[4], (int32x4_t)tmp[6]);
              tmp2[6] = (int16x8_t) vzip1q_s32((int32x4_t)tmp[5], (int32x4_t)tmp[7]);
              tmp2[7] = (int16x8_t) vzip2q_s32((int32x4_t)tmp[5], (int32x4_t)tmp[7]);

              tmp[0] = (int16x8_t) vzip1q_s64((int64x2_t) tmp2[0], (int64x2_t) tmp2[4]);
              tmp[1] = (int16x8_t) vzip2q_s64((int64x2_t) tmp2[0], (int64x2_t) tmp2[4]);
              tmp[2] = (int16x8_t) vzip1q_s64((int64x2_t) tmp2[1], (int64x2_t) tmp2[5]);
              tmp[3] = (int16x8_t) vzip2q_s64((int64x2_t) tmp2[1], (int64x2_t) tmp2[5]);
              tmp[4] = (int16x8_t) vzip1q_s64((int64x2_t) tmp2[2], (int64x2_t) tmp2[6]);
              tmp[5] = (int16x8_t) vzip2q_s64((int64x2_t) tmp2[2], (int64x2_t) tmp2[6]);
              tmp[6] = (int16x8_t) vzip1q_s64((int64x2_t) tmp2[3], (int64x2_t) tmp2[7]);
              tmp[7] = (int16x8_t) vzip2q_s64((int64x2_t) tmp2[3], (int64x2_t) tmp2[7]);

              vst1q_s16(dst+k*8+0, tmp[0]);
              vst1q_s16(dst+k*8+64, tmp[1]);
              vst1q_s16(dst+k*8+128, tmp[2]);
              vst1q_s16(dst+k*8+192, tmp[3]);
              vst1q_s16(dst+k*8+256, tmp[4]);
              vst1q_s16(dst+k*8+320, tmp[5]);
              vst1q_s16(dst+k*8+384, tmp[6]);
              vst1q_s16(dst+k*8+448, tmp[7]);
            }
              src+=8;
              dst += 512;
        }
}

void
vvc_inverse_dct_ii_64_neon_8lines_red16(const int16_t* src, int16_t* dst, ptrdiff_t src_stride,
                      int num_lines, int num_columns, int shift)
{
  int32x4_t add = vdupq_n_s32(1 << (shift - 1));

  int j, k;

  for (j = 0; j < num_lines>>3; j++) {
        /* Utilizing symmetry properties to the maximum to minimize
        the
         * number of multiplications */
         int16x8_t x[8], d[8];
         int32x4_t m[8], a[8], b[4], r[128], r2[128], r3[128];

         //EEEE
         uint16_t DCT_II_64_EEEE[4] ={64,  64,  64,  64};
         x[0 ]=vld1q_s16(src);
         for (k = 0; k < 4; k++) {
           d[ 0] = vdupq_n_s16(DCT_II_64_EEEE[k]);

           a[k*2] = vmull_s16(vget_low_s16(x[0]), vget_low_s16(d[0]));
           a[k*2+1] = vmull_s16(vget_high_s16(x[0]), vget_high_s16(d[0]));
         }

         //EEEO
         x[0 ]=vld1q_s16(src +  8 * src_stride);
         for (k = 0; k < 4; k++) {
           d[ 0] = vdupq_n_s16(DCT_II_64_EEEOT[k * 2 + 0]);

           r3[k*2] = vmull_s16(vget_low_s16(x[0]), vget_low_s16(d[0]));
           r3[k*2+1] = vmull_s16(vget_high_s16(x[0]), vget_high_s16(d[0]));
         }
         //EEE
         /* Combining even and odd terms at each hierarchy levels to
          * calculate the final spatial domain vector */
          for (k = 0; k < 4; k++) {
            r[k*2] = vaddq_s32(r3[k*2], a[k*2]);
            r[k*2+1] = vaddq_s32(r3[k*2+1], a[k*2+1]);

            r[k*2+8] = vsubq_s32(a[6-k*2], r3[6-k*2]);
            r[k*2+9] = vsubq_s32(a[7-k*2], r3[7-k*2]);
          }
         //EEO
         x[0 ]=vld1q_s16(src +  4 * src_stride);
         x[1 ]=vld1q_s16(src +  12 * src_stride);
         for (k = 0; k < 8; k++) {
           d[ 0] = vdupq_n_s16(DCT_II_64_EEOT[k * 4 + 0]);
           d[ 1] = vdupq_n_s16(DCT_II_64_EEOT[k * 4 + 1]);

           m[0] = vmull_s16(vget_low_s16(x[0]), vget_low_s16(d[0]));
           m[1] = vmull_s16(vget_high_s16(x[0]), vget_high_s16(d[0]));
           m[2] = vmull_s16(vget_low_s16(x[1]), vget_low_s16(d[1]));
           m[3] = vmull_s16(vget_high_s16(x[1]), vget_high_s16(d[1]));

           r2[k*2] = vaddq_s32(m[0], m[2]);
           r2[k*2+1] = vaddq_s32(m[1], m[3]);

         }
         //EE
         /* Combining even and odd terms at each hierarchy levels to
          * calculate the final spatial domain vector */
         for (k = 0; k < 8; k++) {
           r3[k*2] = vaddq_s32(r2[k*2], r[k*2]);
           r3[k*2+1] = vaddq_s32(r2[k*2+1], r[k*2+1]);

           r3[k*2+16] = vsubq_s32(r[14-k*2], r2[14-k*2]);
           r3[k*2+17] = vsubq_s32(r[15-k*2], r2[15-k*2]);
         }

         //EO
         x[0]=vld1q_s16(src +  2 * src_stride);
         x[1]=vld1q_s16(src +  6 * src_stride);
         x[2]=vld1q_s16(src +  10 * src_stride);
         x[3]=vld1q_s16(src +  14 * src_stride);
         for (k = 0; k < 16; k++) {
           d[0] = vdupq_n_s16(DCT_II_64_EOT[k * 8 + 0]);
           d[1] = vdupq_n_s16(DCT_II_64_EOT[k * 8 + 1]);
           d[2] = vdupq_n_s16(DCT_II_64_EOT[k * 8 + 2]);
           d[3] = vdupq_n_s16(DCT_II_64_EOT[k * 8 + 3]);

           m[0] = vmull_s16(vget_low_s16(x[0]), vget_low_s16(d[0]));
           m[1] = vmull_s16(vget_high_s16(x[0]), vget_high_s16(d[0]));
           m[2] = vmull_s16(vget_low_s16(x[1]), vget_low_s16(d[1]));
           m[3] = vmull_s16(vget_high_s16(x[1]), vget_high_s16(d[1]));
           m[4] = vmull_s16(vget_low_s16(x[2]), vget_low_s16(d[2]));
           m[5] = vmull_s16(vget_high_s16(x[2]), vget_high_s16(d[2]));
           m[6] = vmull_s16(vget_low_s16(x[3]), vget_low_s16(d[3]));
           m[7] = vmull_s16(vget_high_s16(x[3]), vget_high_s16(d[3]));

           a[0] = vaddq_s32(m[0], m[2]);
           a[1] = vaddq_s32(m[4], m[6]);

           a[2] = vaddq_s32(m[1], m[3]);
           a[3] = vaddq_s32(m[5], m[7]);

           r[k*2] = vaddq_s32(a[0], a[1]);

           r[k*2+1] = vaddq_s32(a[2], a[3]);
         }

         //E
         /* Combining even and odd terms at each hierarchy levels to
          * calculate the final spatial domain vector */
         for (k = 0; k < 16; k++) {
           r2[k*2] = vaddq_s32(r[k*2], r3[k*2]);
           r2[k*2+1] = vaddq_s32(r[k*2+1], r3[k*2+1]);

           r2[k*2+32] = vsubq_s32(r3[30-k*2], r[30-k*2]);
           r2[k*2+33] = vsubq_s32(r3[31-k*2], r[31-k*2]);
         }
         //O
         x[0]=vld1q_s16(src +  src_stride);
         x[1]=vld1q_s16(src +  3  * src_stride);
         x[2]=vld1q_s16(src +  5  * src_stride);
         x[3]=vld1q_s16(src +  7  * src_stride);
         x[4]=vld1q_s16(src +  9  * src_stride);
         x[5]=vld1q_s16(src +  11 * src_stride);
         x[6]=vld1q_s16(src +  13 * src_stride);
         x[7]=vld1q_s16(src +  15 * src_stride);
         for (k = 0; k < 32; k++) {
           d[0] = vdupq_n_s16(DCT_II_64_OT[k * 16 + 0 ]);
           d[1] = vdupq_n_s16(DCT_II_64_OT[k * 16 + 1 ]);
           d[2] = vdupq_n_s16(DCT_II_64_OT[k * 16 + 2 ]);
           d[3] = vdupq_n_s16(DCT_II_64_OT[k * 16 + 3 ]);
           d[4] = vdupq_n_s16(DCT_II_64_OT[k * 16 + 4 ]);
           d[5] = vdupq_n_s16(DCT_II_64_OT[k * 16 + 5 ]);
           d[6] = vdupq_n_s16(DCT_II_64_OT[k * 16 + 6 ]);
           d[7] = vdupq_n_s16(DCT_II_64_OT[k * 16 + 7 ]);

           m[0] = vmull_s16(vget_low_s16(x[0]), vget_low_s16(d[0]));
           m[1] = vmull_s16(vget_high_s16(x[0]), vget_high_s16(d[0]));
           m[2] = vmull_s16(vget_low_s16(x[1]), vget_low_s16(d[1]));
           m[3] = vmull_s16(vget_high_s16(x[1]), vget_high_s16(d[1]));
           m[4] = vmull_s16(vget_low_s16(x[2]), vget_low_s16(d[2]));
           m[5] = vmull_s16(vget_high_s16(x[2]), vget_high_s16(d[2]));
           m[6] = vmull_s16(vget_low_s16(x[3]), vget_low_s16(d[3]));
           m[7] = vmull_s16(vget_high_s16(x[3]), vget_high_s16(d[3]));

           m[8]  = vmull_s16(vget_low_s16(x[4]), vget_low_s16(d[4]));
           m[9]  = vmull_s16(vget_high_s16(x[4]), vget_high_s16(d[4]));
           m[10] = vmull_s16(vget_low_s16(x[5]), vget_low_s16(d[5]));
           m[11] = vmull_s16(vget_high_s16(x[5]), vget_high_s16(d[5]));
           m[12] = vmull_s16(vget_low_s16(x[6]), vget_low_s16(d[6]));
           m[13] = vmull_s16(vget_high_s16(x[6]), vget_high_s16(d[6]));
           m[14] = vmull_s16(vget_low_s16(x[7]), vget_low_s16(d[7]));
           m[15] = vmull_s16(vget_high_s16(x[7]), vget_high_s16(d[7]));

           a[0] = vaddq_s32(m[0 ], m[2 ]);
           a[1] = vaddq_s32(m[4 ], m[6 ]);
           a[2] = vaddq_s32(m[8 ], m[10]);
           a[3] = vaddq_s32(m[12], m[14]);

           a[4] = vaddq_s32(m[1 ], m[3 ]);
           a[5] = vaddq_s32(m[5 ], m[7 ]);
           a[6] = vaddq_s32(m[9 ], m[11]);
           a[7] = vaddq_s32(m[13], m[15]);

           b[0] = vaddq_s32(a[0], a[1]);
           b[1] = vaddq_s32(a[2], a[3]);

           b[2] = vaddq_s32(a[4], a[5]);
           b[3] = vaddq_s32(a[6], a[7]);

           r3[k*2] = vaddq_s32(b[0], b[1]);

           r3[k*2+1] = vaddq_s32(b[2], b[3]);
         }

         //Result
         for (k = 0; k < 32; k++) {
           r[k*2] =   vaddq_s32(r3[k*2], r2[k*2]);
           r[k*2] =   vaddq_s32(r[k*2],add);
           r[k*2] =   vshrq_n_s32(r[k*2],shift);
           r[k*2+1] = vaddq_s32(r3[k*2+1], r2[k*2+1]);
           r[k*2+1] = vaddq_s32(r[k*2+1],add);
           r[k*2+1] = vshrq_n_s32(r[k*2+1],shift);
           r[k*2] = (int32x4_t) vcombine_s16(vqmovn_s32(r[k*2]), vqmovn_s32(r[k*2+1]));

           r[k*2+64] = vsubq_s32(r2[62-k*2], r3[62-k*2]);
           r[k*2+64] = vaddq_s32(r[k*2+64],add);
           r[k*2+64] = vshrq_n_s32(r[k*2+64],shift);
           r[k*2+65] = vsubq_s32(r2[63-k*2], r3[63-k*2]);
           r[k*2+65] = vaddq_s32(r[k*2+65],add);
           r[k*2+65] = vshrq_n_s32(r[k*2+65],shift);
           r[k*2+64] = (int32x4_t) vcombine_s16(vqmovn_s32(r[k*2+64]), vqmovn_s32(r[k*2+65]));
         }

         for (k = 0; k < 8; k++) {
               int16x8_t tmp[8], tmp2[8];
               tmp[0] = vzip1q_s16((int16x8_t) r[k*16 + 0], (int16x8_t) r[k*16 + 2]);
               tmp[1] = vzip2q_s16((int16x8_t) r[k*16 + 0], (int16x8_t) r[k*16 + 2]);
               tmp[2] = vzip1q_s16((int16x8_t) r[k*16 + 4], (int16x8_t) r[k*16 + 6]);
               tmp[3] = vzip2q_s16((int16x8_t) r[k*16 + 4], (int16x8_t) r[k*16 + 6]);
               tmp[4] = vzip1q_s16((int16x8_t) r[k*16 + 8], (int16x8_t) r[k*16 + 10]);
               tmp[5] = vzip2q_s16((int16x8_t) r[k*16 + 8], (int16x8_t) r[k*16 + 10]);
               tmp[6] = vzip1q_s16((int16x8_t) r[k*16 + 12], (int16x8_t) r[k*16 + 14]);
               tmp[7] = vzip2q_s16((int16x8_t) r[k*16 + 12], (int16x8_t) r[k*16 + 14]);

               tmp2[0] = (int16x8_t) vzip1q_s32((int32x4_t)tmp[0], (int32x4_t)tmp[2]);
               tmp2[1] = (int16x8_t) vzip2q_s32((int32x4_t)tmp[0], (int32x4_t)tmp[2]);
               tmp2[2] = (int16x8_t) vzip1q_s32((int32x4_t)tmp[1], (int32x4_t)tmp[3]);
               tmp2[3] = (int16x8_t) vzip2q_s32((int32x4_t)tmp[1], (int32x4_t)tmp[3]);
               tmp2[4] = (int16x8_t) vzip1q_s32((int32x4_t)tmp[4], (int32x4_t)tmp[6]);
               tmp2[5] = (int16x8_t) vzip2q_s32((int32x4_t)tmp[4], (int32x4_t)tmp[6]);
               tmp2[6] = (int16x8_t) vzip1q_s32((int32x4_t)tmp[5], (int32x4_t)tmp[7]);
               tmp2[7] = (int16x8_t) vzip2q_s32((int32x4_t)tmp[5], (int32x4_t)tmp[7]);

               tmp[0] = (int16x8_t) vzip1q_s64((int64x2_t) tmp2[0], (int64x2_t) tmp2[4]);
               tmp[1] = (int16x8_t) vzip2q_s64((int64x2_t) tmp2[0], (int64x2_t) tmp2[4]);
               tmp[2] = (int16x8_t) vzip1q_s64((int64x2_t) tmp2[1], (int64x2_t) tmp2[5]);
               tmp[3] = (int16x8_t) vzip2q_s64((int64x2_t) tmp2[1], (int64x2_t) tmp2[5]);
               tmp[4] = (int16x8_t) vzip1q_s64((int64x2_t) tmp2[2], (int64x2_t) tmp2[6]);
               tmp[5] = (int16x8_t) vzip2q_s64((int64x2_t) tmp2[2], (int64x2_t) tmp2[6]);
               tmp[6] = (int16x8_t) vzip1q_s64((int64x2_t) tmp2[3], (int64x2_t) tmp2[7]);
               tmp[7] = (int16x8_t) vzip2q_s64((int64x2_t) tmp2[3], (int64x2_t) tmp2[7]);

               vst1q_s16(dst+k*8+0, tmp[0]);
               vst1q_s16(dst+k*8+64, tmp[1]);
               vst1q_s16(dst+k*8+128, tmp[2]);
               vst1q_s16(dst+k*8+192, tmp[3]);
               vst1q_s16(dst+k*8+256, tmp[4]);
               vst1q_s16(dst+k*8+320, tmp[5]);
               vst1q_s16(dst+k*8+384, tmp[6]);
               vst1q_s16(dst+k*8+448, tmp[7]);
             }
              src+=8;
              dst += 512;
        }
}
void
vvc_inverse_dct_ii_64_neon_8lines(const int16_t* src, int16_t* dst, ptrdiff_t src_stride,
                      int num_lines, int num_columns, int shift)
{
        int32x4_t add = vdupq_n_s32(1 << (shift - 1));

        int j, k;

        for (j = 0; j < num_lines>>3; j++) {
                /* Utilizing symmetry properties to the maximum to minimize
                the
                 * number of multiplications */
                 int16x8_t x[16], d[16];
                 int32x4_t m[16], a[16], b[8], c[4], r[128], r2[128], r3[128];
                 //EEEEO
                 x[0] = vld1q_s16(src +  16 * src_stride);
                 d[0] = vdupq_n_s16(DCT_II_64_EEEEO[0]);
                 d[1] = vdupq_n_s16(DCT_II_64_EEEEO[1]);

                 b[0] = vmull_s16(vget_low_s16(x[0]), vget_low_s16(d[0]));
                 b[1] = vmull_s16(vget_high_s16(x[0]), vget_high_s16(d[0]));
                 b[2] = vmull_s16(vget_low_s16(x[0]), vget_low_s16(d[1]));
                 b[3] = vmull_s16(vget_high_s16(x[0]), vget_high_s16(d[1]));

                 //EEEEE
                 x[0] = vld1q_s16(src);
                 d[0] = vdupq_n_s16(DCT_II_64_EEEEE[0]);
                 d[1] = vdupq_n_s16(DCT_II_64_EEEEE[1]);

                 c[0] = vmull_s16(vget_low_s16(x[0]), vget_low_s16(d[0]));
                 c[1] = vmull_s16(vget_high_s16(x[0]), vget_high_s16(d[0]));
                 c[2] = vmull_s16(vget_low_s16(x[0]), vget_low_s16(d[1]));
                 c[3] = vmull_s16(vget_high_s16(x[0]), vget_high_s16(d[1]));

                 /* Combining even and odd terms at each hierarchy levels to
                  * calculate the final spatial domain vector */
                 //EEEE
                 for (k = 0; k < 2; k++) {
                   a[k*2] = vaddq_s32(b[k*2], c[k*2]);
                   a[k*2+1] = vaddq_s32(b[k*2+1], c[k*2+1]);

                   a[k*2+4] = vsubq_s32(c[2-k*2], b[2-k*2]);
                   a[k*2+5] = vsubq_s32(c[3-k*2], b[3-k*2]);
                 }

                 //EEEO
                 x[0 ]=vld1q_s16(src +  8 * src_stride);
                 x[1 ]=vld1q_s16(src +  24 * src_stride);
                 for (k = 0; k < 4; k++) {
                   d[ 0] = vdupq_n_s16(DCT_II_64_EEEOT[k * 2 + 0]);
                   d[ 1] = vdupq_n_s16(DCT_II_64_EEEOT[k * 2 + 1]);

                   m[0] = vmull_s16(vget_low_s16(x[0]), vget_low_s16(d[0]));
                   m[1] = vmull_s16(vget_high_s16(x[0]), vget_high_s16(d[0]));
                   m[2] = vmull_s16(vget_low_s16(x[1]), vget_low_s16(d[1]));
                   m[3] = vmull_s16(vget_high_s16(x[1]), vget_high_s16(d[1]));

                   r3[k*2] = vaddq_s32(m[0], m[2]);
                   r3[k*2+1] = vaddq_s32(m[1], m[3]);
                 }
                 //EEE
                 /* Combining even and odd terms at each hierarchy levels to
                  * calculate the final spatial domain vector */
                  for (k = 0; k < 4; k++) {
                    r[k*2] = vaddq_s32(r3[k*2], a[k*2]);
                    r[k*2+1] = vaddq_s32(r3[k*2+1], a[k*2+1]);

                    r[k*2+8] = vsubq_s32(a[6-k*2], r3[6-k*2]);
                    r[k*2+9] = vsubq_s32(a[7-k*2], r3[7-k*2]);
                  }

                 //EEO
                 x[0 ]=vld1q_s16(src +  4 * src_stride);
                 x[1 ]=vld1q_s16(src +  12 * src_stride);
                 x[2 ]=vld1q_s16(src +  20 * src_stride);
                 x[3 ]=vld1q_s16(src +  28 * src_stride);
                 for (k = 0; k < 8; k++) {
                   d[ 0] = vdupq_n_s16(DCT_II_64_EEOT[k * 4 + 0]);
                   d[ 1] = vdupq_n_s16(DCT_II_64_EEOT[k * 4 + 1]);
                   d[ 2] = vdupq_n_s16(DCT_II_64_EEOT[k * 4 + 2]);
                   d[ 3] = vdupq_n_s16(DCT_II_64_EEOT[k * 4 + 3]);

                   m[0] = vmull_s16(vget_low_s16(x[0]), vget_low_s16(d[0]));
                   m[1] = vmull_s16(vget_high_s16(x[0]), vget_high_s16(d[0]));
                   m[2] = vmull_s16(vget_low_s16(x[1]), vget_low_s16(d[1]));
                   m[3] = vmull_s16(vget_high_s16(x[1]), vget_high_s16(d[1]));
                   m[4] = vmull_s16(vget_low_s16(x[2]), vget_low_s16(d[2]));
                   m[5] = vmull_s16(vget_high_s16(x[2]), vget_high_s16(d[2]));
                   m[6] = vmull_s16(vget_low_s16(x[3]), vget_low_s16(d[3]));
                   m[7] = vmull_s16(vget_high_s16(x[3]), vget_high_s16(d[3]));

                   a[0] = vaddq_s32(m[0], m[2]);
                   a[1] = vaddq_s32(m[4], m[6]);

                   a[2] = vaddq_s32(m[1], m[3]);
                   a[3] = vaddq_s32(m[5], m[7]);

                   r2[k*2] = vaddq_s32(a[0], a[1]);
                   r2[k*2+1] = vaddq_s32(a[2], a[3]);
                 }
                 //EE
                 /* Combining even and odd terms at each hierarchy levels to
                  * calculate the final spatial domain vector */
                 for (k = 0; k < 8; k++) {
                   r3[k*2] = vaddq_s32(r2[k*2], r[k*2]);
                   r3[k*2+1] = vaddq_s32(r2[k*2+1], r[k*2+1]);

                   r3[k*2+16] = vsubq_s32(r[14-k*2], r2[14-k*2]);
                   r3[k*2+17] = vsubq_s32(r[15-k*2], r2[15-k*2]);
                 }

                 //EO
                 x[0 ]=vld1q_s16(src +  2 * src_stride);
                 x[1 ]=vld1q_s16(src +  6 * src_stride);
                 x[2 ]=vld1q_s16(src +  10 * src_stride);
                 x[3 ]=vld1q_s16(src +  14 * src_stride);
                 x[4 ]=vld1q_s16(src +  18 * src_stride);
                 x[5 ]=vld1q_s16(src +  22 * src_stride);
                 x[6 ]=vld1q_s16(src +  26 * src_stride);
                 x[7 ]=vld1q_s16(src +  30 * src_stride);
                 for (k = 0; k < 16; k++) {
                   d[ 0] = vdupq_n_s16(DCT_II_64_EOT[k * 8 + 0]);
                   d[ 1] = vdupq_n_s16(DCT_II_64_EOT[k * 8 + 1]);
                   d[ 2] = vdupq_n_s16(DCT_II_64_EOT[k * 8 + 2]);
                   d[ 3] = vdupq_n_s16(DCT_II_64_EOT[k * 8 + 3]);
                   d[ 4] = vdupq_n_s16(DCT_II_64_EOT[k * 8 + 4]);
                   d[ 5] = vdupq_n_s16(DCT_II_64_EOT[k * 8 + 5]);
                   d[ 6] = vdupq_n_s16(DCT_II_64_EOT[k * 8 + 6]);
                   d[ 7] = vdupq_n_s16(DCT_II_64_EOT[k * 8 + 7]);

                   m[0] = vmull_s16(vget_low_s16(x[0]), vget_low_s16(d[0]));
                   m[1] = vmull_s16(vget_high_s16(x[0]), vget_high_s16(d[0]));
                   m[2] = vmull_s16(vget_low_s16(x[1]), vget_low_s16(d[1]));
                   m[3] = vmull_s16(vget_high_s16(x[1]), vget_high_s16(d[1]));
                   m[4] = vmull_s16(vget_low_s16(x[2]), vget_low_s16(d[2]));
                   m[5] = vmull_s16(vget_high_s16(x[2]), vget_high_s16(d[2]));
                   m[6] = vmull_s16(vget_low_s16(x[3]), vget_low_s16(d[3]));
                   m[7] = vmull_s16(vget_high_s16(x[3]), vget_high_s16(d[3]));

                   m[8]  = vmull_s16(vget_low_s16(x[4]), vget_low_s16(d[4]));
                   m[9]  = vmull_s16(vget_high_s16(x[4]), vget_high_s16(d[4]));
                   m[10] = vmull_s16(vget_low_s16(x[5]), vget_low_s16(d[5]));
                   m[11] = vmull_s16(vget_high_s16(x[5]), vget_high_s16(d[5]));
                   m[12] = vmull_s16(vget_low_s16(x[6]), vget_low_s16(d[6]));
                   m[13] = vmull_s16(vget_high_s16(x[6]), vget_high_s16(d[6]));
                   m[14] = vmull_s16(vget_low_s16(x[7]), vget_low_s16(d[7]));
                   m[15] = vmull_s16(vget_high_s16(x[7]), vget_high_s16(d[7]));

                   a[0] = vaddq_s32(m[0 ], m[2 ]);
                   a[1] = vaddq_s32(m[4 ], m[6 ]);
                   a[2] = vaddq_s32(m[8 ], m[10]);
                   a[3] = vaddq_s32(m[12], m[14]);

                   a[4] = vaddq_s32(m[1 ], m[3 ]);
                   a[5] = vaddq_s32(m[5 ], m[7 ]);
                   a[6] = vaddq_s32(m[9 ], m[11]);
                   a[7] = vaddq_s32(m[13], m[15]);

                   b[0] = vaddq_s32(a[0], a[1]);
                   b[1] = vaddq_s32(a[2], a[3]);

                   b[2] = vaddq_s32(a[4], a[5]);
                   b[3] = vaddq_s32(a[6], a[7]);

                   r[k*2] = vaddq_s32(b[0], b[1]);

                   r[k*2+1] = vaddq_s32(b[2], b[3]);
                 }

                 //E
                 /* Combining even and odd terms at each hierarchy levels to
                  * calculate the final spatial domain vector */
                 for (k = 0; k < 16; k++) {
                   r2[k*2] = vaddq_s32(r[k*2], r3[k*2]);
                   r2[k*2+1] = vaddq_s32(r[k*2+1], r3[k*2+1]);

                   r2[k*2+32] = vsubq_s32(r3[30-k*2], r[30-k*2]);
                   r2[k*2+33] = vsubq_s32(r3[31-k*2], r[31-k*2]);
                 }
                 //O
                 x[0 ]=vld1q_s16(src +  src_stride);
                 x[1 ]=vld1q_s16(src +  3  * src_stride);
                 x[2 ]=vld1q_s16(src +  5  * src_stride);
                 x[3 ]=vld1q_s16(src +  7  * src_stride);
                 x[4 ]=vld1q_s16(src +  9  * src_stride);
                 x[5 ]=vld1q_s16(src +  11 * src_stride);
                 x[6 ]=vld1q_s16(src +  13 * src_stride);
                 x[7 ]=vld1q_s16(src +  15 * src_stride);
                 x[8 ]=vld1q_s16(src +  17 * src_stride);
                 x[9 ]=vld1q_s16(src +  19 * src_stride);
                 x[10]=vld1q_s16(src +  21 * src_stride);
                 x[11]=vld1q_s16(src +  23 * src_stride);
                 x[12]=vld1q_s16(src +  25 * src_stride);
                 x[13]=vld1q_s16(src +  27 * src_stride);
                 x[14]=vld1q_s16(src +  29 * src_stride);
                 x[15]=vld1q_s16(src +  31 * src_stride);
                 for (k = 0; k < 32; k++) {
                   d[ 0] = vdupq_n_s16(DCT_II_64_OT[k * 16 + 0 ]);
                   d[ 1] = vdupq_n_s16(DCT_II_64_OT[k * 16 + 1 ]);
                   d[ 2] = vdupq_n_s16(DCT_II_64_OT[k * 16 + 2 ]);
                   d[ 3] = vdupq_n_s16(DCT_II_64_OT[k * 16 + 3 ]);
                   d[ 4] = vdupq_n_s16(DCT_II_64_OT[k * 16 + 4 ]);
                   d[ 5] = vdupq_n_s16(DCT_II_64_OT[k * 16 + 5 ]);
                   d[ 6] = vdupq_n_s16(DCT_II_64_OT[k * 16 + 6 ]);
                   d[ 7] = vdupq_n_s16(DCT_II_64_OT[k * 16 + 7 ]);
                   d[ 8] = vdupq_n_s16(DCT_II_64_OT[k * 16 + 8 ]);
                   d[ 9] = vdupq_n_s16(DCT_II_64_OT[k * 16 + 9 ]);
                   d[10] = vdupq_n_s16(DCT_II_64_OT[k * 16 + 10]);
                   d[11] = vdupq_n_s16(DCT_II_64_OT[k * 16 + 11]);
                   d[12] = vdupq_n_s16(DCT_II_64_OT[k * 16 + 12]);
                   d[13] = vdupq_n_s16(DCT_II_64_OT[k * 16 + 13]);
                   d[14] = vdupq_n_s16(DCT_II_64_OT[k * 16 + 14]);
                   d[15] = vdupq_n_s16(DCT_II_64_OT[k * 16 + 15]);

                   m[0] = vmull_s16(vget_low_s16(x[0]), vget_low_s16(d[0]));
                   m[1] = vmull_s16(vget_high_s16(x[0]), vget_high_s16(d[0]));
                   m[2] = vmull_s16(vget_low_s16(x[1]), vget_low_s16(d[1]));
                   m[3] = vmull_s16(vget_high_s16(x[1]), vget_high_s16(d[1]));
                   m[4] = vmull_s16(vget_low_s16(x[2]), vget_low_s16(d[2]));
                   m[5] = vmull_s16(vget_high_s16(x[2]), vget_high_s16(d[2]));
                   m[6] = vmull_s16(vget_low_s16(x[3]), vget_low_s16(d[3]));
                   m[7] = vmull_s16(vget_high_s16(x[3]), vget_high_s16(d[3]));

                   m[8]  = vmull_s16(vget_low_s16(x[4]), vget_low_s16(d[4]));
                   m[9]  = vmull_s16(vget_high_s16(x[4]), vget_high_s16(d[4]));
                   m[10] = vmull_s16(vget_low_s16(x[5]), vget_low_s16(d[5]));
                   m[11] = vmull_s16(vget_high_s16(x[5]), vget_high_s16(d[5]));
                   m[12] = vmull_s16(vget_low_s16(x[6]), vget_low_s16(d[6]));
                   m[13] = vmull_s16(vget_high_s16(x[6]), vget_high_s16(d[6]));
                   m[14] = vmull_s16(vget_low_s16(x[7]), vget_low_s16(d[7]));
                   m[15] = vmull_s16(vget_high_s16(x[7]), vget_high_s16(d[7]));

                   m[16] = vmull_s16(vget_low_s16(x[8]), vget_low_s16(d[8]));
                   m[17] = vmull_s16(vget_high_s16(x[8]), vget_high_s16(d[8]));
                   m[18] = vmull_s16(vget_low_s16(x[9]), vget_low_s16(d[9]));
                   m[19] = vmull_s16(vget_high_s16(x[9]), vget_high_s16(d[9]));
                   m[20] = vmull_s16(vget_low_s16(x[10]), vget_low_s16(d[10]));
                   m[21] = vmull_s16(vget_high_s16(x[10]), vget_high_s16(d[10]));
                   m[22] = vmull_s16(vget_low_s16(x[11]), vget_low_s16(d[11]));
                   m[23] = vmull_s16(vget_high_s16(x[11]), vget_high_s16(d[11]));

                   m[24]  = vmull_s16(vget_low_s16(x[12]), vget_low_s16(d[12]));
                   m[25]  = vmull_s16(vget_high_s16(x[12]), vget_high_s16(d[12]));
                   m[26] = vmull_s16(vget_low_s16(x[13]), vget_low_s16(d[13]));
                   m[27] = vmull_s16(vget_high_s16(x[13]), vget_high_s16(d[13]));
                   m[28] = vmull_s16(vget_low_s16(x[14]), vget_low_s16(d[14]));
                   m[29] = vmull_s16(vget_high_s16(x[14]), vget_high_s16(d[14]));
                   m[30] = vmull_s16(vget_low_s16(x[15]), vget_low_s16(d[15]));
                   m[31] = vmull_s16(vget_high_s16(x[15]), vget_high_s16(d[15]));

                   a[0] = vaddq_s32(m[0 ], m[2 ]);
                   a[1] = vaddq_s32(m[4 ], m[6 ]);
                   a[2] = vaddq_s32(m[8 ], m[10]);
                   a[3] = vaddq_s32(m[12], m[14]);
                   a[4] = vaddq_s32(m[16], m[18]);
                   a[5] = vaddq_s32(m[20], m[22]);
                   a[6] = vaddq_s32(m[24], m[26]);
                   a[7] = vaddq_s32(m[28], m[30]);

                   a[ 8] = vaddq_s32(m[1 ], m[3 ]);
                   a[ 9] = vaddq_s32(m[5 ], m[7 ]);
                   a[10] = vaddq_s32(m[9 ], m[11]);
                   a[11] = vaddq_s32(m[13], m[15]);
                   a[12] = vaddq_s32(m[17], m[19]);
                   a[13] = vaddq_s32(m[21], m[23]);
                   a[14] = vaddq_s32(m[25], m[27]);
                   a[15] = vaddq_s32(m[29], m[31]);

                   b[0] = vaddq_s32(a[0], a[1]);
                   b[1] = vaddq_s32(a[2], a[3]);
                   b[2] = vaddq_s32(a[4], a[5]);
                   b[3] = vaddq_s32(a[6], a[7]);

                   b[4] = vaddq_s32(a[ 8], a[9]);
                   b[5] = vaddq_s32(a[10], a[11]);
                   b[6] = vaddq_s32(a[12], a[13]);
                   b[7] = vaddq_s32(a[14], a[15]);

                   c[0] = vaddq_s32(b[0], b[1]);
                   c[1] = vaddq_s32(b[2], b[3]);

                   c[2] = vaddq_s32(b[4], b[5]);
                   c[3] = vaddq_s32(b[6], b[7]);

                   r3[k*2] = vaddq_s32(c[0], c[1]);
                   r3[k*2+1] = vaddq_s32(c[2], c[3]);
                 }

                 //Result
                 for (k = 0; k < 32; k++) {
                   r[k*2] =   vaddq_s32(r3[k*2], r2[k*2]);
                   r[k*2] =   vaddq_s32(r[k*2],add);
                   r[k*2] =   vshrq_n_s32(r[k*2],shift);
                   r[k*2+1] = vaddq_s32(r3[k*2+1], r2[k*2+1]);
                   r[k*2+1] = vaddq_s32(r[k*2+1],add);
                   r[k*2+1] = vshrq_n_s32(r[k*2+1],shift);
                   r[k*2] = (int32x4_t)vcombine_s16(vqmovn_s32(r[k*2]), vqmovn_s32(r[k*2+1]));

                   r[k*2+64] = vsubq_s32(r2[62-k*2], r3[62-k*2]);
                   r[k*2+64] = vaddq_s32(r[k*2+64],add);
                   r[k*2+64] = vshrq_n_s32(r[k*2+64],shift);
                   r[k*2+65] = vsubq_s32(r2[63-k*2], r3[63-k*2]);
                   r[k*2+65] = vaddq_s32(r[k*2+65],add);
                   r[k*2+65] = vshrq_n_s32(r[k*2+65],shift);
                   r[k*2+64] = (int32x4_t)vcombine_s16(vqmovn_s32(r[k*2+64]), vqmovn_s32(r[k*2+65]));
                 }

                 for (k = 0; k < 8; k++) {
                       int16x8_t tmp[8], tmp2[8];
                       tmp[0] = vzip1q_s16((int16x8_t) r[k*16 + 0], (int16x8_t) r[k*16 + 2]);
                       tmp[1] = vzip2q_s16((int16x8_t) r[k*16 + 0], (int16x8_t) r[k*16 + 2]);
                       tmp[2] = vzip1q_s16((int16x8_t) r[k*16 + 4], (int16x8_t) r[k*16 + 6]);
                       tmp[3] = vzip2q_s16((int16x8_t) r[k*16 + 4], (int16x8_t) r[k*16 + 6]);
                       tmp[4] = vzip1q_s16((int16x8_t) r[k*16 + 8], (int16x8_t) r[k*16 + 10]);
                       tmp[5] = vzip2q_s16((int16x8_t) r[k*16 + 8], (int16x8_t) r[k*16 + 10]);
                       tmp[6] = vzip1q_s16((int16x8_t) r[k*16 + 12], (int16x8_t) r[k*16 + 14]);
                       tmp[7] = vzip2q_s16((int16x8_t) r[k*16 + 12], (int16x8_t) r[k*16 + 14]);

                       tmp2[0] = (int16x8_t) vzip1q_s32((int32x4_t)tmp[0], (int32x4_t)tmp[2]);
                       tmp2[1] = (int16x8_t) vzip2q_s32((int32x4_t)tmp[0], (int32x4_t)tmp[2]);
                       tmp2[2] = (int16x8_t) vzip1q_s32((int32x4_t)tmp[1], (int32x4_t)tmp[3]);
                       tmp2[3] = (int16x8_t) vzip2q_s32((int32x4_t)tmp[1], (int32x4_t)tmp[3]);
                       tmp2[4] = (int16x8_t) vzip1q_s32((int32x4_t)tmp[4], (int32x4_t)tmp[6]);
                       tmp2[5] = (int16x8_t) vzip2q_s32((int32x4_t)tmp[4], (int32x4_t)tmp[6]);
                       tmp2[6] = (int16x8_t) vzip1q_s32((int32x4_t)tmp[5], (int32x4_t)tmp[7]);
                       tmp2[7] = (int16x8_t) vzip2q_s32((int32x4_t)tmp[5], (int32x4_t)tmp[7]);

                       tmp[0] = (int16x8_t) vzip1q_s64((int64x2_t) tmp2[0], (int64x2_t) tmp2[4]);
                       tmp[1] = (int16x8_t) vzip2q_s64((int64x2_t) tmp2[0], (int64x2_t) tmp2[4]);
                       tmp[2] = (int16x8_t) vzip1q_s64((int64x2_t) tmp2[1], (int64x2_t) tmp2[5]);
                       tmp[3] = (int16x8_t) vzip2q_s64((int64x2_t) tmp2[1], (int64x2_t) tmp2[5]);
                       tmp[4] = (int16x8_t) vzip1q_s64((int64x2_t) tmp2[2], (int64x2_t) tmp2[6]);
                       tmp[5] = (int16x8_t) vzip2q_s64((int64x2_t) tmp2[2], (int64x2_t) tmp2[6]);
                       tmp[6] = (int16x8_t) vzip1q_s64((int64x2_t) tmp2[3], (int64x2_t) tmp2[7]);
                       tmp[7] = (int16x8_t) vzip2q_s64((int64x2_t) tmp2[3], (int64x2_t) tmp2[7]);

                       vst1q_s16(dst+k*8+0, tmp[0]);
                       vst1q_s16(dst+k*8+64, tmp[1]);
                       vst1q_s16(dst+k*8+128, tmp[2]);
                       vst1q_s16(dst+k*8+192, tmp[3]);
                       vst1q_s16(dst+k*8+256, tmp[4]);
                       vst1q_s16(dst+k*8+320, tmp[5]);
                       vst1q_s16(dst+k*8+384, tmp[6]);
                       vst1q_s16(dst+k*8+448, tmp[7]);
                     }
                      src+=8;
                      dst += 512;
                }
}


void
vvc_inverse_dct_ii_64_neon(const int16_t *src, int16_t *dst, ptrdiff_t src_stride,
                          int num_lines, int line_brk, int shift)
{
    /*FIXME it would be more efficient to get per reduced lines speialised
      functions instead off checking */
    if (num_lines>>3){
      if (line_brk > 16) {
          vvc_inverse_dct_ii_64_neon_8lines(src, dst, src_stride, num_lines&0xF8, line_brk, shift);
      }
      else if (line_brk > 8){
          vvc_inverse_dct_ii_64_neon_8lines_red16(src, dst, src_stride, num_lines&0xF8, line_brk, shift);
      }
      else {
          vvc_inverse_dct_ii_64_neon_8lines_red8(src, dst, src_stride, num_lines&0xF8, line_brk, shift);
      }
        src += (num_lines & 0xF8);
        dst += (num_lines >> 3) << 9;
    }

    if (!(num_lines & 0x7)) return;
    vvc_inverse_dct_ii_64(src, dst, src_stride, num_lines & 0x7, line_brk, shift);

    // if (num_lines & 0x4){
    //   vvc_inverse_dct_ii_64(src, dst, src_stride, num_lines & 0x4, line_brk, shift);
    //   src += 4;
    //   dst += 256;
    // }
    //
    // if (num_lines & 0x2){
    //   vvc_inverse_dct_ii_64(src, dst, src_stride, num_lines & 0x2, line_brk, shift);
    //   src += 2;
    //   dst += 128;
    // }
    //
    // if (num_lines & 0x1){
    //     vvc_inverse_dct_ii_64(src, dst, src_stride, num_lines & 0x1, line_brk, shift);
    // }
}
#if 0
void vvc_inverse_dst_vii_4_neon(const int16_t *src, int16_t *dst, ptrdiff_t src_stride,
                               int num_lines, int line_brk, int shift){
    inverse_neon2_B4(src, dst, src_stride, shift, num_lines, DST_VII_4);
}

void vvc_inverse_dst_vii_8_neon(const int16_t *src, int16_t *dst, ptrdiff_t src_stride,
                               int num_lines, int line_brk, int shift){
    inverse_neon2_B8(src, dst, src_stride, shift, num_lines, DST_VII_8);
}

void vvc_inverse_dst_vii_16_neon(const int16_t *src, int16_t *dst, ptrdiff_t src_stride,
                                int num_lines, int line_brk, int shift){
    inverse_neon2_B16(src, dst, src_stride, shift, num_lines, DST_VII_16);
}

void vvc_inverse_dst_vii_32_neon(const int16_t *src, int16_t *dst, ptrdiff_t src_stride,
                                int num_lines, int line_brk, int shift){
    inverse_neon2_B32(src, dst, src_stride, shift, num_lines, DST_VII_32);
}

void vvc_inverse_dct_viii_4_neon(const int16_t *src, int16_t *dst, ptrdiff_t src_stride,
                                int num_lines, int line_brk, int shift){
    inverse_neon2_B4(src, dst, src_stride, shift, num_lines, DCT_VIII_4);
}

void vvc_inverse_dct_viii_8_neon(const int16_t *src, int16_t *dst, ptrdiff_t src_stride,
                                int num_lines, int line_brk, int shift){
    inverse_neon2_B8(src, dst, src_stride, shift, num_lines, DCT_VIII_8);
}

void vvc_inverse_dct_viii_16_neon(const int16_t *src, int16_t *dst, ptrdiff_t src_stride,
                                 int num_lines, int line_brk, int shift){
    inverse_neon2_B16(src, dst, src_stride, shift, num_lines, DCT_VIII_16);
}

void vvc_inverse_dct_viii_32_neon(const int16_t *src, int16_t *dst, ptrdiff_t src_stride,
                                 int num_lines, int line_brk, int shift){
    inverse_neon2_B32(src, dst, src_stride, shift, num_lines, DCT_VIII_32);
}
#endif


void
vvc_inverse_dct_ii_dc_neon(int16_t *const dst, int log2_tb_w, int log2_tb_h,
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
    int16x8_t x0 = vdupq_n_s16(value);
    switch (log2_tb_w){
        case 3:{
        for (i = 0; i < tb_h; ++i){
            vst1q_s16(_dst, x0);
            _dst += tb_w;
        }
                   break;
               }
        case 4:{
        for (i = 0; i < tb_h; ++i){
            vst1q_s16(_dst, x0);
            vst1q_s16(&_dst[8], x0);
            _dst += tb_w;
        }
                   break;
               }
        case 5:{
        for (i = 0; i < tb_h; ++i){
            vst1q_s16(_dst, x0);
            vst1q_s16(&_dst[8], x0);
            vst1q_s16(&_dst[16], x0);
            vst1q_s16(&_dst[24], x0);
            _dst += tb_w;
        }
                   break;
               }
        default:
            vvc_inverse_dct_ii_dc(dst, log2_tb_w, log2_tb_h, dc_val);
    }
}


void rcn_init_tr_functions_neon(struct RCNFunctions *const rcn_funcs){
  // rcn_funcs->tr.func[DST_VII][2] = &vvc_inverse_dst_vii_4_neon;
  // rcn_funcs->tr.func[DST_VII][3] = &vvc_inverse_dst_vii_8_neon;
  // rcn_funcs->tr.func[DST_VII][4] = &vvc_inverse_dst_vii_16_neon;
  // rcn_funcs->tr.func[DST_VII][5] = &vvc_inverse_dst_vii_32_neon;
  
  // rcn_funcs->tr.func[DCT_VIII][2] = &vvc_inverse_dct_viii_4_neon;
  // rcn_funcs->tr.func[DCT_VIII][3] = &vvc_inverse_dct_viii_8_neon;
  // rcn_funcs->tr.func[DCT_VIII][4] = &vvc_inverse_dct_viii_16_neon;
  // rcn_funcs->tr.func[DCT_VIII][5] = &vvc_inverse_dct_viii_32_neon;

  rcn_funcs->tr.func[DCT_II][1] = &vvc_inverse_dct_ii_2_neon;
  rcn_funcs->tr.func[DCT_II][2] = &vvc_inverse_dct_ii_4_neon;
  rcn_funcs->tr.func[DCT_II][3] = &vvc_inverse_dct_ii_8_neon;
  rcn_funcs->tr.func[DCT_II][4] = &vvc_inverse_dct_ii_16_neon;
  rcn_funcs->tr.func[DCT_II][5] = &vvc_inverse_dct_ii_32_neon;
  rcn_funcs->tr.func[DCT_II][6] = &vvc_inverse_dct_ii_64_neon;

  rcn_funcs->tr.dc = &vvc_inverse_dct_ii_dc_neon;
}
