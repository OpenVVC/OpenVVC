#if defined _MSC_VER
#include <tmmintrin.h>
#else
#include <x86intrin.h>
#endif

#include <stdint.h>
#include "rcn_alf.h"
#include "rcn_structures.h"

#define process2coeffs5x5(i, ptr0, ptr1, ptr2, ptr3) { \
                    const __m128i val00 = _mm_sub_epi16(_mm_loadu_si128((const __m128i *) (ptr0)), cur);\
                    const __m128i val10 = _mm_sub_epi16(_mm_loadu_si128((const __m128i *) (ptr2)), cur);\
                    const __m128i val01 = _mm_sub_epi16(_mm_loadu_si128((const __m128i *) (ptr1)), cur);\
                    const __m128i val11 = _mm_sub_epi16(_mm_loadu_si128((const __m128i *) (ptr3)), cur);\
\
                    __m128i val01A = _mm_unpacklo_epi16(val00, val10);\
                    __m128i val01B = _mm_unpackhi_epi16(val00, val10);\
                    __m128i val01C = _mm_unpacklo_epi16(val01, val11);\
                    __m128i val01D = _mm_unpackhi_epi16(val01, val11);\
\
                    __m128i limit01A = params[1][i];\
\
                    val01A = _mm_min_epi16(val01A, limit01A);\
                    val01B = _mm_min_epi16(val01B, limit01A);\
                    val01C = _mm_min_epi16(val01C, limit01A);\
                    val01D = _mm_min_epi16(val01D, limit01A);\
\
                    limit01A = _mm_sub_epi16(_mm_setzero_si128(), limit01A);\
\
                    val01A = _mm_max_epi16(val01A, limit01A);\
                    val01B = _mm_max_epi16(val01B, limit01A);\
                    val01C = _mm_max_epi16(val01C, limit01A);\
                    val01D = _mm_max_epi16(val01D, limit01A);\
\
                    val01A = _mm_add_epi16(val01A, val01C);\
                    val01B = _mm_add_epi16(val01B, val01D);\
\
                    __m128i coeff01A = params[0][i];\
\
                    accumA = _mm_add_epi32(accumA, _mm_madd_epi16(val01A, coeff01A));\
                    accumB = _mm_add_epi32(accumB, _mm_madd_epi16(val01B, coeff01A));\
                };\

#define process2coeffs7x7(i, ptr0, ptr1, ptr2, ptr3) { \
        const __m128i val00 = _mm_sub_epi16(_mm_loadu_si128((const __m128i *) (ptr0)), cur); \
        const __m128i val10 = _mm_sub_epi16(_mm_loadu_si128((const __m128i *) (ptr2)), cur); \
        const __m128i val01 = _mm_sub_epi16(_mm_loadu_si128((const __m128i *) (ptr1)), cur); \
        const __m128i val11 = _mm_sub_epi16(_mm_loadu_si128((const __m128i *) (ptr3)), cur); \
         \
        __m128i val01A = _mm_unpacklo_epi16(val00, val10); \
        __m128i val01B = _mm_unpackhi_epi16(val00, val10); \
        __m128i val01C = _mm_unpacklo_epi16(val01, val11); \
        __m128i val01D = _mm_unpackhi_epi16(val01, val11); \
         \
        __m128i limit01A = params[0][1][i]; \
        __m128i limit01B = params[1][1][i]; \
         \
        val01A = _mm_min_epi16(val01A, limit01A); \
        val01B = _mm_min_epi16(val01B, limit01B); \
        val01C = _mm_min_epi16(val01C, limit01A); \
        val01D = _mm_min_epi16(val01D, limit01B); \
         \
        limit01A = _mm_sub_epi16(_mm_setzero_si128(), limit01A); \
        limit01B = _mm_sub_epi16(_mm_setzero_si128(), limit01B); \
         \
        val01A = _mm_max_epi16(val01A, limit01A); \
        val01B = _mm_max_epi16(val01B, limit01B); \
        val01C = _mm_max_epi16(val01C, limit01A); \
        val01D = _mm_max_epi16(val01D, limit01B); \
         \
        val01A = _mm_add_epi16(val01A, val01C); \
        val01B = _mm_add_epi16(val01B, val01D); \
         \
        const __m128i coeff01A = params[0][0][i]; \
        const __m128i coeff01B = params[1][0][i]; \
         \
        accumA = _mm_add_epi32(accumA, _mm_madd_epi16(val01A, coeff01A)); \
        accumB = _mm_add_epi32(accumB, _mm_madd_epi16(val01B, coeff01B)); \
}; \


static void
simdFilter5x5Blk(int16_t *const dst, const int16_t *const src,
                 const int dstStride, const int srcStride,
                 Area blk_dst,
                 const int16_t *const filter_set, const int16_t *const clip_set,
                 const int ctu_height, int virbnd_pos)
{
    const int SHIFT = NUM_BITS - 1;
    const int ROUND = 1 << (SHIFT - 1);

    const size_t STEP_X = 8;
    const size_t STEP_Y = 4;

    const int clpMin = 0;
    const int clpMax = (1<<10) - 1;

    int16_t * _src = (int16_t *) src;
    int16_t * _dst = dst;

    const __m128i mmOffset = _mm_set1_epi32(ROUND);
    const __m128i mmMin = _mm_set1_epi16( clpMin );
    const __m128i mmMax = _mm_set1_epi16( clpMax );

    __m128i params[2][3];
    __m128i fs   = _mm_loadu_si128((__m128i *) filter_set);
    params[0][0] = _mm_shuffle_epi32(fs, 0x00);
    params[0][1] = _mm_shuffle_epi32(fs, 0x55);
    params[0][2] = _mm_shuffle_epi32(fs, 0xaa);
    __m128i fc   = _mm_loadu_si128((__m128i *) clip_set);
    params[1][0] = _mm_shuffle_epi32(fc, 0x00);
    params[1][1] = _mm_shuffle_epi32(fc, 0x55);
    params[1][2] = _mm_shuffle_epi32(fc, 0xaa);

    for (size_t i = 0; i < blk_dst.height; i += STEP_Y) {
        for (size_t j = 0; j < blk_dst.width; j += STEP_X) {
            for (size_t ii = 0; ii < STEP_Y; ii++) {
                const uint16_t *pImg0, *pImg1, *pImg2, *pImg3, *pImg4;

                pImg0 = (uint16_t*)_src + j + ii * srcStride ;
                pImg1 = pImg0 + srcStride;
                pImg2 = pImg0 - srcStride;
                pImg3 = pImg1 + srcStride;
                pImg4 = pImg2 - srcStride;

                __m128i cur = _mm_loadu_si128((const __m128i *) pImg0);

                __m128i accumA, accumB;

                accumA = mmOffset;
                accumB = mmOffset;

                process2coeffs5x5(0, pImg3 + 0, pImg4 + 0, pImg1 + 1, pImg2 - 1);
                process2coeffs5x5(1, pImg1 + 0, pImg2 + 0, pImg1 - 1, pImg2 + 1);
                process2coeffs5x5(2, pImg0 + 2, pImg0 - 2, pImg0 + 1, pImg0 - 1);

                accumA = _mm_srai_epi32(accumA, SHIFT);
                accumB = _mm_srai_epi32(accumB, SHIFT);

                accumA = _mm_packs_epi32(accumA, accumB);
                accumA = _mm_add_epi16(accumA, cur);
                accumA = _mm_min_epi16(mmMax, _mm_max_epi16(accumA, mmMin));

                _mm_storeu_si128((__m128i *) (_dst + ii * dstStride + j), accumA);
            }
        }

        _src += srcStride * STEP_Y;
        _dst += dstStride * STEP_Y;
    }
}

static void
simdFilter5x5BlkVB(int16_t *const dst, const int16_t *const src,
                 const int dstStride, const int srcStride,
                 Area blk_dst,
                 const int16_t *const filter_set, const int16_t *const clip_set,
                 const int ctu_height, int virbnd_pos)
{
    const int SHIFT = NUM_BITS - 1;
    const int ROUND = 1 << (SHIFT - 1);

    const size_t STEP_X = 8;
    const size_t STEP_Y = 4;

    const int clpMin = 0;
    const int clpMax = (1<<10) - 1;

    int16_t * _src = (int16_t *) src;
    int16_t * _dst = dst;

    const __m128i mmOffset = _mm_set1_epi32(ROUND);
    const __m128i mmMin = _mm_set1_epi16( clpMin );
    const __m128i mmMax = _mm_set1_epi16( clpMax );

    const __m128i mmOffsetborder = _mm_set1_epi32(1 << ((SHIFT + 3) - 1));
    const int SHIFTborder = SHIFT+3;


    __m128i params[2][3];
    __m128i fs   = _mm_loadu_si128((__m128i *) filter_set);
    params[0][0] = _mm_shuffle_epi32(fs, 0x00);
    params[0][1] = _mm_shuffle_epi32(fs, 0x55);
    params[0][2] = _mm_shuffle_epi32(fs, 0xaa);
    __m128i fc   = _mm_loadu_si128((__m128i *) clip_set);
    params[1][0] = _mm_shuffle_epi32(fc, 0x00);
    params[1][1] = _mm_shuffle_epi32(fc, 0x55);
    params[1][2] = _mm_shuffle_epi32(fc, 0xaa);

    for (size_t i = 0; i < blk_dst.height; i += STEP_Y) {
        for (size_t j = 0; j < blk_dst.width; j += STEP_X) {
            for (size_t ii = 0; ii < STEP_Y; ii++) {
                const uint16_t *pImg0, *pImg1, *pImg2, *pImg3, *pImg4;

                pImg0 = (uint16_t*)_src + j + ii * srcStride ;
                pImg1 = pImg0 + srcStride;
                pImg2 = pImg0 - srcStride;
                pImg3 = pImg1 + srcStride;
                pImg4 = pImg2 - srcStride;

                const int yVb = (blk_dst.y + i + ii) & (ctu_height - 1);

                if (yVb < virbnd_pos && (yVb >= virbnd_pos - 2)) {
                    pImg1 = (yVb == virbnd_pos - 1) ? pImg0 : pImg1;
                    pImg3 = (yVb >= virbnd_pos - 2) ? pImg1 : pImg3;

                    pImg2 = (yVb == virbnd_pos - 1) ? pImg0 : pImg2;
                    pImg4 = (yVb >= virbnd_pos - 2) ? pImg2 : pImg4;
                } else if (yVb >= virbnd_pos && (yVb <= virbnd_pos + 1)) {
                    pImg2 = (yVb == virbnd_pos) ? pImg0 : pImg2;
                    pImg4 = (yVb <= virbnd_pos + 1) ? pImg2 : pImg4;

                    pImg1 = (yVb == virbnd_pos) ? pImg0 : pImg1;
                    pImg3 = (yVb <= virbnd_pos + 1) ? pImg1 : pImg3;
                }

                __m128i cur = _mm_loadu_si128((const __m128i *) pImg0);

                uint8_t isNearVBabove = yVb < virbnd_pos && (yVb >= virbnd_pos - 1);
                uint8_t isNearVBbelow = yVb >= virbnd_pos && (yVb <= virbnd_pos);

                __m128i accumA, accumB;
                if (!(isNearVBabove || isNearVBbelow)) {
                    accumA = mmOffset;
                    accumB = mmOffset;
                } else {
                    //Rounding offset fix
                    accumA = mmOffsetborder;
                    accumB = mmOffsetborder;
                }

                process2coeffs5x5(0, pImg3 + 0, pImg4 + 0, pImg1 + 1, pImg2 - 1);
                process2coeffs5x5(1, pImg1 + 0, pImg2 + 0, pImg1 - 1, pImg2 + 1);
                process2coeffs5x5(2, pImg0 + 2, pImg0 - 2, pImg0 + 1, pImg0 - 1);

                if (!(isNearVBabove || isNearVBbelow)) {
                    accumA = _mm_srai_epi32(accumA, SHIFT);
                    accumB = _mm_srai_epi32(accumB, SHIFT);
                } else {
                    //Rounding offset fix
                    accumA = _mm_srai_epi32(accumA, SHIFTborder);
                    accumB = _mm_srai_epi32(accumB, SHIFTborder);
                }

                accumA = _mm_packs_epi32(accumA, accumB);
                accumA = _mm_add_epi16(accumA, cur);
                accumA = _mm_min_epi16(mmMax, _mm_max_epi16(accumA, mmMin));

                if (j + STEP_X <= blk_dst.width) {
                  _mm_storeu_si128((__m128i *) (_dst + ii * dstStride + j), accumA);
                } else {
                  _mm_storel_epi64((__m128i *) (_dst + ii * dstStride + j), accumA);
                }
            }
        }

        _src += srcStride * STEP_Y;
        _dst += dstStride * STEP_Y;
    }
}

static void
simdFilter7x7Blk(uint8_t * class_idx_arr, uint8_t * transpose_idx_arr, int16_t *const dst, int16_t *const src, const int dstStride, const int srcStride,
                         Area blk_dst, const int16_t *filter_set, const int16_t *clip_set,
                         const int ctu_height, int virbnd_pos)
{
    const int SHIFT = NUM_BITS - 1;
    const int ROUND = 1 << (SHIFT - 1);

    const size_t STEP_X = 8;
    const size_t STEP_Y = 4;

    const int clpMin = 0;
    const int clpMax = (1<<10) - 1;

    int16_t * _src = src;
    int16_t * _dst = dst;

    int transpose_idx = 0;
    int class_idx = 0;

    const __m128i mmOffset = _mm_set1_epi32(ROUND);
    const __m128i mmMin = _mm_set1_epi16( clpMin );
    const __m128i mmMax = _mm_set1_epi16( clpMax );

    for (size_t i = 0; i < blk_dst.height; i += STEP_Y)
    {
        for (size_t j = 0; j < blk_dst.width; j += STEP_X)
        {
            __m128i params[2][2][6];

            for (int k = 0; k < 2; ++k)
            {
                transpose_idx = transpose_idx_arr[(i>>2) * CLASSIFICATION_BLK_SIZE + (j>>2) + k];
                class_idx = class_idx_arr[(i>>2) * CLASSIFICATION_BLK_SIZE + (j>>2) + k];

                const __m128i rawCoeffLo = _mm_loadu_si128((const __m128i *) (filter_set + transpose_idx * MAX_NUM_ALF_CLASSES * MAX_NUM_ALF_LUMA_COEFF + class_idx * MAX_NUM_ALF_LUMA_COEFF));
                const __m128i rawCoeffHi = _mm_loadl_epi64((const __m128i *) (filter_set + transpose_idx * MAX_NUM_ALF_CLASSES * MAX_NUM_ALF_LUMA_COEFF + class_idx * MAX_NUM_ALF_LUMA_COEFF + 8));
                const __m128i rawClipLo  = _mm_loadu_si128((const __m128i *) (clip_set + transpose_idx * MAX_NUM_ALF_CLASSES * MAX_NUM_ALF_LUMA_COEFF + class_idx * MAX_NUM_ALF_LUMA_COEFF));
                const __m128i rawClipHi  = _mm_loadl_epi64((const __m128i *) (clip_set + transpose_idx * MAX_NUM_ALF_CLASSES * MAX_NUM_ALF_LUMA_COEFF + class_idx * MAX_NUM_ALF_LUMA_COEFF + 8));

                params[k][0][0] = _mm_shuffle_epi32(rawCoeffLo, 0x00);
                params[k][0][1] = _mm_shuffle_epi32(rawCoeffLo, 0x55);
                params[k][0][2] = _mm_shuffle_epi32(rawCoeffLo, 0xaa);
                params[k][0][3] = _mm_shuffle_epi32(rawCoeffLo, 0xff);
                params[k][0][4] = _mm_shuffle_epi32(rawCoeffHi, 0x00);
                params[k][0][5] = _mm_shuffle_epi32(rawCoeffHi, 0x55);
                params[k][1][0] = _mm_shuffle_epi32(rawClipLo, 0x00);
                params[k][1][1] = _mm_shuffle_epi32(rawClipLo, 0x55);
                params[k][1][2] = _mm_shuffle_epi32(rawClipLo, 0xaa);
                params[k][1][3] = _mm_shuffle_epi32(rawClipLo, 0xff);
                params[k][1][4] = _mm_shuffle_epi32(rawClipHi, 0x00);
                params[k][1][5] = _mm_shuffle_epi32(rawClipHi, 0x55);
            }

            for (size_t ii = 0; ii < STEP_Y; ii++)
            {
                const uint16_t *pImg0, *pImg1, *pImg2, *pImg3, *pImg4, *pImg5, *pImg6;

                pImg0 = (uint16_t *)_src + j + ii * srcStride;
                pImg1 = pImg0 + srcStride;
                pImg2 = pImg0 - srcStride;
                pImg3 = pImg1 + srcStride;
                pImg4 = pImg2 - srcStride;
                pImg5 = pImg3 + srcStride;
                pImg6 = pImg4 - srcStride;

                __m128i cur = _mm_loadu_si128((const __m128i *) pImg0);

                __m128i accumA = mmOffset;
                __m128i accumB = mmOffset;

                process2coeffs7x7(0, pImg5 + 0, pImg6 + 0, pImg3 + 1, pImg4 - 1);
                process2coeffs7x7(1, pImg3 + 0, pImg4 + 0, pImg3 - 1, pImg4 + 1);
                process2coeffs7x7(2, pImg1 + 2, pImg2 - 2, pImg1 + 1, pImg2 - 1);
                process2coeffs7x7(3, pImg1 + 0, pImg2 + 0, pImg1 - 1, pImg2 + 1);
                process2coeffs7x7(4, pImg1 - 2, pImg2 + 2, pImg0 + 3, pImg0 - 3);
                process2coeffs7x7(5, pImg0 + 2, pImg0 - 2, pImg0 + 1, pImg0 - 1);

                accumA = _mm_srai_epi32(accumA, SHIFT);
                accumB = _mm_srai_epi32(accumB, SHIFT);

                accumA = _mm_packs_epi32(accumA, accumB);
                accumA = _mm_add_epi16(accumA, cur);
                accumA = _mm_min_epi16(mmMax, _mm_max_epi16(accumA, mmMin));

                _mm_storeu_si128((__m128i *) (_dst + ii * dstStride + j), accumA);
            }
        }

        _src += srcStride * STEP_Y;
        _dst += dstStride * STEP_Y;
    }
}

static void
simdFilter7x7BlkVB(uint8_t * class_idx_arr, uint8_t * transpose_idx_arr, int16_t *const dst, int16_t *const src, const int dstStride, const int srcStride,
                         Area blk_dst, const int16_t *filter_set, const int16_t *clip_set,
                         const int ctu_height, int virbnd_pos)
{
    const int SHIFT = NUM_BITS - 1;
    const int ROUND = 1 << (SHIFT - 1);

    const size_t STEP_X = 8;
    const size_t STEP_Y = 4;

    const int clpMin = 0;
    const int clpMax = (1<<10) - 1;

    int16_t * _src = src;
    int16_t * _dst = dst;

    int transpose_idx = 0;
    int class_idx = 0;

    const __m128i mmOffset = _mm_set1_epi32(ROUND);
    const __m128i mmMin = _mm_set1_epi16( clpMin );
    const __m128i mmMax = _mm_set1_epi16( clpMax );

    const __m128i mmOffsetborder = _mm_set1_epi32(1 << ((SHIFT + 3) - 1));
    const int SHIFTborder = SHIFT+3;

    for (size_t i = 0; i < blk_dst.height; i += STEP_Y)
    {
        for (size_t j = 0; j < blk_dst.width; j += STEP_X)
        {
            __m128i params[2][2][6];

            for (int k = 0; k < 2; ++k)
            {
                transpose_idx = transpose_idx_arr[(i>>2) * CLASSIFICATION_BLK_SIZE + (j>>2) + k];
                class_idx = class_idx_arr[(i>>2) * CLASSIFICATION_BLK_SIZE + (j>>2) + k];

                const __m128i rawCoeffLo = _mm_loadu_si128((const __m128i *) (filter_set + transpose_idx * MAX_NUM_ALF_CLASSES * MAX_NUM_ALF_LUMA_COEFF + class_idx * MAX_NUM_ALF_LUMA_COEFF));
                const __m128i rawCoeffHi = _mm_loadl_epi64((const __m128i *) (filter_set + transpose_idx * MAX_NUM_ALF_CLASSES * MAX_NUM_ALF_LUMA_COEFF + class_idx * MAX_NUM_ALF_LUMA_COEFF + 8));
                const __m128i rawClipLo  = _mm_loadu_si128((const __m128i *) (clip_set + transpose_idx * MAX_NUM_ALF_CLASSES * MAX_NUM_ALF_LUMA_COEFF + class_idx * MAX_NUM_ALF_LUMA_COEFF));
                const __m128i rawClipHi  = _mm_loadl_epi64((const __m128i *) (clip_set + transpose_idx * MAX_NUM_ALF_CLASSES * MAX_NUM_ALF_LUMA_COEFF + class_idx * MAX_NUM_ALF_LUMA_COEFF + 8));

                params[k][0][0] = _mm_shuffle_epi32(rawCoeffLo, 0x00);
                params[k][0][1] = _mm_shuffle_epi32(rawCoeffLo, 0x55);
                params[k][0][2] = _mm_shuffle_epi32(rawCoeffLo, 0xaa);
                params[k][0][3] = _mm_shuffle_epi32(rawCoeffLo, 0xff);
                params[k][0][4] = _mm_shuffle_epi32(rawCoeffHi, 0x00);
                params[k][0][5] = _mm_shuffle_epi32(rawCoeffHi, 0x55);
                params[k][1][0] = _mm_shuffle_epi32(rawClipLo, 0x00);
                params[k][1][1] = _mm_shuffle_epi32(rawClipLo, 0x55);
                params[k][1][2] = _mm_shuffle_epi32(rawClipLo, 0xaa);
                params[k][1][3] = _mm_shuffle_epi32(rawClipLo, 0xff);
                params[k][1][4] = _mm_shuffle_epi32(rawClipHi, 0x00);
                params[k][1][5] = _mm_shuffle_epi32(rawClipHi, 0x55);
            }

            for (size_t ii = 0; ii < STEP_Y; ii++)
            {
                const uint16_t *pImg0, *pImg1, *pImg2, *pImg3, *pImg4, *pImg5, *pImg6;

                pImg0 = (uint16_t *)_src + j + ii * srcStride;
                pImg1 = pImg0 + srcStride;
                pImg2 = pImg0 - srcStride;
                pImg3 = pImg1 + srcStride;
                pImg4 = pImg2 - srcStride;
                pImg5 = pImg3 + srcStride;
                pImg6 = pImg4 - srcStride;

                const int yVb = (blk_dst.y + i + ii) & (ctu_height - 1);
                if (yVb < virbnd_pos && (yVb >= virbnd_pos - 4))   // above
                {
                    pImg1 = (yVb == virbnd_pos - 1) ? pImg0 : pImg1;
                    pImg3 = (yVb >= virbnd_pos - 2) ? pImg1 : pImg3;
                    pImg5 = (yVb >= virbnd_pos - 3) ? pImg3 : pImg5;

                    pImg2 = (yVb == virbnd_pos - 1) ? pImg0 : pImg2;
                    pImg4 = (yVb >= virbnd_pos - 2) ? pImg2 : pImg4;
                    pImg6 = (yVb >= virbnd_pos - 3) ? pImg4 : pImg6;
                }
                else if (yVb >= virbnd_pos && (yVb <= virbnd_pos + 3))   // bottom
                {
                    pImg2 = (yVb == virbnd_pos) ? pImg0 : pImg2;
                    pImg4 = (yVb <= virbnd_pos + 1) ? pImg2 : pImg4;
                    pImg6 = (yVb <= virbnd_pos + 2) ? pImg4 : pImg6;

                    pImg1 = (yVb == virbnd_pos) ? pImg0 : pImg1;
                    pImg3 = (yVb <= virbnd_pos + 1) ? pImg1 : pImg3;
                    pImg5 = (yVb <= virbnd_pos + 2) ? pImg3 : pImg5;
                }

                __m128i cur = _mm_loadu_si128((const __m128i *) pImg0);

                uint8_t isNearVBabove = yVb < virbnd_pos && (yVb >= virbnd_pos - 1);
                uint8_t isNearVBbelow = yVb >= virbnd_pos && (yVb <= virbnd_pos);

                __m128i accumA, accumB;
                if (!(isNearVBabove || isNearVBbelow))
                {
                    accumA = mmOffset;
                    accumB = mmOffset;
                }
                else
                {
                    //Rounding offset fix
                    accumA = mmOffsetborder;
                    accumB = mmOffsetborder;
                }

                process2coeffs7x7(0, pImg5 + 0, pImg6 + 0, pImg3 + 1, pImg4 - 1);
                process2coeffs7x7(1, pImg3 + 0, pImg4 + 0, pImg3 - 1, pImg4 + 1);
                process2coeffs7x7(2, pImg1 + 2, pImg2 - 2, pImg1 + 1, pImg2 - 1);
                process2coeffs7x7(3, pImg1 + 0, pImg2 + 0, pImg1 - 1, pImg2 + 1);
                process2coeffs7x7(4, pImg1 - 2, pImg2 + 2, pImg0 + 3, pImg0 - 3);
                process2coeffs7x7(5, pImg0 + 2, pImg0 - 2, pImg0 + 1, pImg0 - 1);

                if (!(isNearVBabove || isNearVBbelow))
                {
                    accumA = _mm_srai_epi32(accumA, SHIFT);
                    accumB = _mm_srai_epi32(accumB, SHIFT);
                }
                else
                {
                    //Rounding offset fix
                    accumA = _mm_srai_epi32(accumA, SHIFTborder);
                    accumB = _mm_srai_epi32(accumB, SHIFTborder);
                }
                accumA = _mm_packs_epi32(accumA, accumB);
                accumA = _mm_add_epi16(accumA, cur);
                accumA = _mm_min_epi16(mmMax, _mm_max_epi16(accumA, mmMin));

                // if (j>=x0 && i+ii>=y0 && j<=w_max && i+ii<=h_max)
                _mm_storeu_si128((__m128i *) (_dst + ii * dstStride + j), accumA);
            }
        }

        _src += srcStride * STEP_Y;
        _dst += dstStride * STEP_Y;
    }
}

#define selectEvenValues(dest, src0, src1) {\
  __m128i a0 = _mm_shufflelo_epi16(src0, 0xD8);\
  __m128i a1 = _mm_shufflelo_epi16(src1, 0xD8);\
  a0 = _mm_shufflehi_epi16(a0, 0xD8);\
  a1 = _mm_shufflehi_epi16(a1, 0xD8);\
  __m128i b0 = _mm_unpacklo_epi32(a0, a1);\
  __m128i b1 = _mm_unpackhi_epi32(a0, a1);\
  dest = _mm_unpacklo_epi32(b0, b1);\
}

#define selectOddValues(dest, src0, src1) {\
  __m128i a0 = _mm_shufflelo_epi16(src0, 0xD8);\
  __m128i a1 = _mm_shufflelo_epi16(src1, 0xD8);\
  a0 = _mm_shufflehi_epi16(a0, 0xD8);\
  a1 = _mm_shufflehi_epi16(a1, 0xD8);\
  afficherVecteur8SSE128(a0);\
  afficherVecteur8SSE128(a1);\
  printf("\n");\
  __m128i b0 = _mm_unpacklo_epi32(a0, a1);\
  __m128i b1 = _mm_unpackhi_epi32(a0, a1);\
  dest = _mm_unpackhi_epi32(b0, b1);\
}

#define selectEvenOddValues(even, odd, src0, src1) {\
  __m128i a0 = _mm_shufflelo_epi16(src0, 0xD8);\
  __m128i a1 = _mm_shufflelo_epi16(src1, 0xD8);\
  a0 = _mm_shufflehi_epi16(a0, 0xD8);\
  a1 = _mm_shufflehi_epi16(a1, 0xD8);\
  __m128i b0 = _mm_unpacklo_epi32(a0, a1);\
  __m128i b1 = _mm_unpackhi_epi32(a0, a1);\
  even = _mm_unpacklo_epi32(b0, b1);\
  odd = _mm_unpackhi_epi32(b0, b1);\
}

void cc_alf_filterBlkVB_sse(int16_t * chroma_dst, int16_t * luma_src, const int chr_stride, const int luma_stride,
                            const Area blk_dst, const uint8_t c_id, const int16_t *filt_coeff,
                            const int vbCTUHeight, int vbPos)
{
  const size_t STEP_X = 8;
  const size_t STEP_Y = 4;

  //ATTENTION: scaleX et Y fixed to 1 (en 4 2 0)
  const int scaleX             = 1;
  const int scaleY             = 1;

  __m128i filter01 = _mm_set1_epi32((filt_coeff[0] & 0xFFFF) | ((filt_coeff[1] & 0xFFFF)<<16));
  __m128i filter23 = _mm_set1_epi32((filt_coeff[2] & 0xFFFF) | ((filt_coeff[3] & 0xFFFF)<<16));
  __m128i filter45 = _mm_set1_epi32((filt_coeff[4] & 0xFFFF) | ((filt_coeff[5] & 0xFFFF)<<16));
  __m128i filter6  = _mm_set1_epi16(filt_coeff[6]);

  const int scale_bits = 7;
  __m128i scale_offset = _mm_set1_epi32((1 << scale_bits ) >> 1);

  //BITDEPTH: uniquement pour bitdepth 10
  const int bit_depth = 10;
  __m128i offset = _mm_set1_epi16((1 << bit_depth) >> 1);

  __m128i clip_h = _mm_set1_epi16((1<<bit_depth) - 1 );
  __m128i clip_l = _mm_setzero_si128();

  for( int i = 0; i < blk_dst.height; i += STEP_Y )
  {
    for( int j = 0; j < blk_dst.width; j += STEP_X )
    {
      for( int ii = 0; ii < STEP_Y; ii++ )
      {
        int row       = ii;
        int col       = j;
        int16_t *srcSelf  = chroma_dst + col + row * chr_stride;

        int offset1 = luma_stride;
        int offset2 = -luma_stride;
        int offset3 = 2 * luma_stride;
        row <<= scaleY;
        col <<= scaleX;
        const int16_t *srcCross = luma_src + col + row * luma_stride;

        int pos = ((blk_dst.y + i + ii) << scaleY) & (vbCTUHeight - 1);
        if (!(scaleY == 0 && (pos == vbPos || pos == vbPos + 1)))
        {
          if (pos == (vbPos - 2) || pos == (vbPos + 1))
          {
            offset3 = offset1;
          }
          else if (pos == (vbPos - 1) || pos == vbPos)
          {
            offset1 = 0;
            offset2 = 0;
            offset3 = 0;
          }
          __m128i val0, val1, val2, val3, val4, val5, val6;
          __m128i curr;

          __m128i self = _mm_loadu_si128((const __m128i *) (srcSelf));
          __m128i x00 = _mm_loadu_si128((const __m128i *) (srcCross + offset2));
          __m128i x01 = _mm_loadu_si128((const __m128i *) (srcCross + offset2 + 8));
          __m128i x10 = _mm_loadu_si128((const __m128i *) (srcCross));
          __m128i x11 = _mm_loadu_si128((const __m128i *) (srcCross + 8));
          __m128i x20 = _mm_loadu_si128((const __m128i *) (srcCross + offset1));
          __m128i x21 = _mm_loadu_si128((const __m128i *) (srcCross + offset1 + 8));
          __m128i x30 = _mm_loadu_si128((const __m128i *) (srcCross + offset3));
          __m128i x31 = _mm_loadu_si128((const __m128i *) (srcCross + offset3 + 8));

          selectEvenValues(val0, x00, x01);
          selectEvenOddValues(curr, val2, x10, x11);
          selectEvenOddValues(val4, val5, x20, x21);
          selectEvenValues(val6, x30, x31);

          val1 = _mm_setr_epi16(srcCross[-1], 0, 0, 0, 0, 0, 0, 0);
          val3 = _mm_setr_epi16(srcCross[offset1 -1], 0, 0, 0, 0, 0, 0, 0);

          val1 = _mm_add_epi16(val1, _mm_bslli_si128(val2, 2));
          val3 = _mm_add_epi16(val3, _mm_bslli_si128(val5, 2));

          val0 = _mm_sub_epi16(val0, curr);
          val1 = _mm_sub_epi16(val1, curr);
          val2 = _mm_sub_epi16(val2, curr);
          val3 = _mm_sub_epi16(val3, curr);
          val4 = _mm_sub_epi16(val4, curr);
          val5 = _mm_sub_epi16(val5, curr);
          val6 = _mm_sub_epi16(val6, curr);

          __m128i val01l = _mm_unpacklo_epi16(val0, val1);
          __m128i val01h = _mm_unpackhi_epi16(val0, val1);
          __m128i val23l = _mm_unpacklo_epi16(val2, val3);
          __m128i val23h = _mm_unpackhi_epi16(val2, val3);
          __m128i val45l = _mm_unpacklo_epi16(val4, val5);
          __m128i val45h = _mm_unpackhi_epi16(val4, val5);

          val01l = _mm_madd_epi16(val01l, filter01);
          val01h = _mm_madd_epi16(val01h, filter01);
          val23l = _mm_madd_epi16(val23l, filter23);
          val23h = _mm_madd_epi16(val23h, filter23);
          val45l = _mm_madd_epi16(val45l, filter45);
          val45h = _mm_madd_epi16(val45h, filter45);
          __m128i val6lo = _mm_mullo_epi16(val6, filter6);
          __m128i val6hi = _mm_mulhi_epi16(val6, filter6);
          __m128i val6l = _mm_unpacklo_epi16(val6lo, val6hi);
          __m128i val6h = _mm_unpackhi_epi16(val6lo, val6hi);

          __m128i a0 = _mm_add_epi32(val01l, val23l);
          __m128i a1 = _mm_add_epi32(val45l, val6l);
          __m128i a2 = _mm_add_epi32(val01h, val23h);
          __m128i a3 = _mm_add_epi32(val45h, val6h);

          a0 = _mm_add_epi32(a0, a1);
          a1 = _mm_add_epi32(a2, a3);

          a0 = _mm_add_epi32(a0, scale_offset);
          a0 = _mm_srai_epi32(a0, scale_bits);

          a1 = _mm_add_epi32(a1, scale_offset);
          a1 = _mm_srai_epi32(a1, scale_bits);

          a0 = _mm_packs_epi32(a0, a1);
          a0 = _mm_add_epi16(a0, offset);

          a0 = _mm_min_epi16(a0, clip_h);
          a0 = _mm_max_epi16(a0, clip_l);

          a0 = _mm_sub_epi16(a0, offset);
          a0 = _mm_add_epi16(a0, self);

          a0 = _mm_min_epi16(a0, clip_h);
          a0 = _mm_max_epi16(a0, clip_l);

          _mm_storeu_si128((__m128i *) (srcSelf), a0);
        }
      }
    }
    chroma_dst += chr_stride * STEP_Y;
    luma_src += luma_stride * STEP_Y << scaleY;
  }
}

static void simdDeriveClassificationBlk(uint8_t * class_idx_arr, uint8_t * transpose_idx_arr,
                                        int16_t *const src, const int stride, const Area blk,
                                        const int shift, const int ctu_s, int virbnd_pos)
{
    int blk_h = blk.height;
    int blk_w = blk.width;

    uint16_t colSums[18][40];
    int i;
    const uint32_t ctb_msk = ctu_s - 1;

    for (i = 0; i < blk_h + 4; i += 2) {
        int yoffset = (i - 3) * stride - 3;
        const int16_t *src0 = &src[yoffset];
        const int16_t *src1 = &src[yoffset + stride];
        const int16_t *src2 = &src[yoffset + stride * 2];
        const int16_t *src3 = &src[yoffset + stride * 3];

        const int y = blk.y - 2 + i;
        int j;

        if (y > 0 && (y & ctb_msk) == virbnd_pos - 2) {
            src3 = src2;
        } else if (y > 0 && (y & ctb_msk) == virbnd_pos) {
            src0 = src1;
        }

        __m128i prev = _mm_setzero_si128();

        for (j = 0; j < blk_w + 4; j += 8) {
            const __m128i x0 = _mm_loadu_si128((const __m128i *) (src0 + j));
            const __m128i x1 = _mm_loadu_si128((const __m128i *) (src1 + j));
            const __m128i x2 = _mm_loadu_si128((const __m128i *) (src2 + j));
            const __m128i x3 = _mm_loadu_si128((const __m128i *) (src3 + j));

            const __m128i x4 = _mm_loadu_si128((const __m128i *) (src0 + j + 2));
            const __m128i x5 = _mm_loadu_si128((const __m128i *) (src1 + j + 2));
            const __m128i x6 = _mm_loadu_si128((const __m128i *) (src2 + j + 2));
            const __m128i x7 = _mm_loadu_si128((const __m128i *) (src3 + j + 2));

            const __m128i nw = _mm_blend_epi16(x0, x1, 0xaa);
            const __m128i n  = _mm_blend_epi16(x0, x5, 0x55);
            const __m128i ne = _mm_blend_epi16(x4, x5, 0xaa);
            const __m128i w  = _mm_blend_epi16(x1, x2, 0xaa);
            const __m128i e  = _mm_blend_epi16(x5, x6, 0xaa);
            const __m128i sw = _mm_blend_epi16(x2, x3, 0xaa);
            const __m128i s  = _mm_blend_epi16(x2, x7, 0x55);
            const __m128i se = _mm_blend_epi16(x6, x7, 0xaa);

            __m128i c = _mm_blend_epi16(x1, x6, 0x55);
            c         = _mm_add_epi16(c, c);
            __m128i d = _mm_shuffle_epi8(c, _mm_setr_epi8(2, 3, 0, 1, 6, 7, 4, 5, 10, 11, 8, 9, 14, 15, 12, 13));

            const __m128i ver = _mm_abs_epi16(_mm_sub_epi16(c, _mm_add_epi16(n, s)));
            const __m128i hor = _mm_abs_epi16(_mm_sub_epi16(d, _mm_add_epi16(w, e)));
            const __m128i di0 = _mm_abs_epi16(_mm_sub_epi16(d, _mm_add_epi16(nw, se)));
            const __m128i di1 = _mm_abs_epi16(_mm_sub_epi16(d, _mm_add_epi16(ne, sw)));

            const __m128i hv  = _mm_hadd_epi16(ver, hor);
            const __m128i di  = _mm_hadd_epi16(di0, di1);
            const __m128i all = _mm_hadd_epi16(hv, di);

            const __m128i t = _mm_blend_epi16(all, prev, 0xaa);
            _mm_storeu_si128((__m128i *) &colSums[i >> 1][j], _mm_hadd_epi16(t, all));
            prev = all;
        }
    }

    for (i = 0; i < (blk_h >> 1); i += 4) {
        __m128i class_idx[4], transpose_idx[4];
        for (size_t k = 0; k < 4; k++) {
            __m128i x0, x1, x2, x3, x4, x5, x6, x7;

            const uint32_t z = (2 * i + blk.y) & ctb_msk;
            const uint32_t z2 = (2 * i + 4 + blk.y) & ctb_msk;

            x0 = (z == virbnd_pos) ? _mm_setzero_si128() : _mm_loadu_si128((__m128i *) &colSums[i + 0][(k*8) + 4]);
            x1 = _mm_loadu_si128((__m128i *) &colSums[i + 1][(k*8) + 4]);
            x2 = _mm_loadu_si128((__m128i *) &colSums[i + 2][(k*8) + 4]);
            x3 = (z == virbnd_pos - 4) ? _mm_setzero_si128() : _mm_loadu_si128((__m128i *) &colSums[i + 3][(k*8) + 4]);

            x4 = (z2 == virbnd_pos) ? _mm_setzero_si128() : _mm_loadu_si128((__m128i *) &colSums[i + 2][(k*8) + 4]);
            x5 = _mm_loadu_si128((__m128i *) &colSums[i + 3][(k*8) + 4]);
            x6 = _mm_loadu_si128((__m128i *) &colSums[i + 4][(k*8) + 4]);
            x7 = (z2 == virbnd_pos - 4) ? _mm_setzero_si128() : _mm_loadu_si128((__m128i *) &colSums[i + 5][(k*8) + 4]);

            __m128i x0l = _mm_cvtepu16_epi32(x0);
            __m128i x0h = _mm_unpackhi_epi16(x0, _mm_setzero_si128());
            __m128i x1l = _mm_cvtepu16_epi32(x1);
            __m128i x1h = _mm_unpackhi_epi16(x1, _mm_setzero_si128());
            __m128i x2l = _mm_cvtepu16_epi32(x2);
            __m128i x2h = _mm_unpackhi_epi16(x2, _mm_setzero_si128());
            __m128i x3l = _mm_cvtepu16_epi32(x3);
            __m128i x3h = _mm_unpackhi_epi16(x3, _mm_setzero_si128());
            __m128i x4l = _mm_cvtepu16_epi32(x4);
            __m128i x4h = _mm_unpackhi_epi16(x4, _mm_setzero_si128());
            __m128i x5l = _mm_cvtepu16_epi32(x5);
            __m128i x5h = _mm_unpackhi_epi16(x5, _mm_setzero_si128());
            __m128i x6l = _mm_cvtepu16_epi32(x6);
            __m128i x6h = _mm_unpackhi_epi16(x6, _mm_setzero_si128());
            __m128i x7l = _mm_cvtepu16_epi32(x7);
            __m128i x7h = _mm_unpackhi_epi16(x7, _mm_setzero_si128());

            x0l = _mm_add_epi32(x0l, x1l);
            x2l = _mm_add_epi32(x2l, x3l);
            x4l = _mm_add_epi32(x4l, x5l);
            x6l = _mm_add_epi32(x6l, x7l);
            x0h = _mm_add_epi32(x0h, x1h);
            x2h = _mm_add_epi32(x2h, x3h);
            x4h = _mm_add_epi32(x4h, x5h);
            x6h = _mm_add_epi32(x6h, x7h);

            x0l = _mm_add_epi32(x0l, x2l);
            x4l = _mm_add_epi32(x4l, x6l);
            x0h = _mm_add_epi32(x0h, x2h);
            x4h = _mm_add_epi32(x4h, x6h);

            x2l = _mm_unpacklo_epi32(x0l, x4l);
            x2h = _mm_unpackhi_epi32(x0l, x4l);
            x6l = _mm_unpacklo_epi32(x0h, x4h);
            x6h = _mm_unpackhi_epi32(x0h, x4h);

            __m128i sumV  = _mm_unpacklo_epi32(x2l, x6l);
            __m128i sumH  = _mm_unpackhi_epi32(x2l, x6l);
            __m128i sumD0 = _mm_unpacklo_epi32(x2h, x6h);
            __m128i sumD1 = _mm_unpackhi_epi32(x2h, x6h);

            __m128i tempAct = _mm_add_epi32(sumV, sumH);

            const uint32_t scale  = (z == virbnd_pos - 4 || z == virbnd_pos) ? 96 : 64;
            const uint32_t scale2 = (z2 == virbnd_pos - 4 || z2 == virbnd_pos) ? 96 : 64;
            __m128i activity = _mm_mullo_epi32(tempAct, _mm_unpacklo_epi64(_mm_set1_epi32(scale), _mm_set1_epi32(scale2)));
            activity         = _mm_srl_epi32(activity, _mm_cvtsi32_si128(shift));
            activity         = _mm_min_epi32(activity, _mm_set1_epi32(15));
            __m128i classIdx = _mm_shuffle_epi8(_mm_setr_epi8(0, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 4), activity);

            __m128i dirTempHVMinus1 = _mm_cmpgt_epi32(sumV, sumH);
            __m128i hv1             = _mm_max_epi32(sumV, sumH);
            __m128i hv0             = _mm_min_epi32(sumV, sumH);

            __m128i dirTempDMinus1 = _mm_cmpgt_epi32(sumD0, sumD1);
            __m128i d1             = _mm_max_epi32(sumD0, sumD1);
            __m128i d0             = _mm_min_epi32(sumD0, sumD1);

            __m128i a      = _mm_xor_si128(_mm_mullo_epi32(d1, hv0), _mm_set1_epi32(0x80000000));
            __m128i b      = _mm_xor_si128(_mm_mullo_epi32(hv1, d0), _mm_set1_epi32(0x80000000));
            __m128i dirIdx = _mm_cmpgt_epi32(a, b);
            __m128i hvd1   = _mm_blendv_epi8(hv1, d1, dirIdx);
            __m128i hvd0   = _mm_blendv_epi8(hv0, d0, dirIdx);

            __m128i strength1 = _mm_cmpgt_epi32(hvd1, _mm_add_epi32(hvd0, hvd0));
            __m128i strength2 = _mm_cmpgt_epi32(_mm_add_epi32(hvd1, hvd1), _mm_add_epi32(hvd0, _mm_slli_epi32(hvd0, 3)));
            __m128i offset    = _mm_and_si128(strength1, _mm_set1_epi32(5));
            classIdx          = _mm_add_epi32(classIdx, offset);
            classIdx          = _mm_add_epi32(classIdx, _mm_and_si128(strength2, _mm_set1_epi32(5)));
            offset            = _mm_andnot_si128(dirIdx, offset);
            offset            = _mm_add_epi32(offset, offset);
            classIdx          = _mm_add_epi32(classIdx, offset);

            __m128i transposeIdx = _mm_set1_epi32(3);
            transposeIdx         = _mm_add_epi32(transposeIdx, dirTempHVMinus1);
            transposeIdx         = _mm_add_epi32(transposeIdx, dirTempDMinus1);
            transposeIdx         = _mm_add_epi32(transposeIdx, dirTempDMinus1);

            class_idx[k] = _mm_shuffle_epi8(classIdx, _mm_setr_epi8(0, 4, 8, 12, 0, 4, 8, 12, 0, 4, 8, 12, 0, 4, 8, 12));
            transpose_idx[k] = _mm_shuffle_epi8(transposeIdx, _mm_setr_epi8(0, 4, 8, 12, 0, 4, 8, 12, 0, 4, 8, 12, 0, 4, 8, 12));
        }


        __m128i c1, c2, t1, t2;

        c1 = _mm_unpacklo_epi16(class_idx[0], class_idx[1]);
        c2 = _mm_unpacklo_epi16(class_idx[2], class_idx[3]);
        c1 = _mm_unpacklo_epi32(c1, c2);

        t1 = _mm_unpacklo_epi16(transpose_idx[0], transpose_idx[1]);
        t2 = _mm_unpacklo_epi16(transpose_idx[2], transpose_idx[3]);
        t1 = _mm_unpacklo_epi32(t1, t2);


        int yOffset = (2*i + blk.y) & ctb_msk;
        int xOffset = (blk.x) & ctb_msk;
        _mm_storel_epi64((__m128i *) (class_idx_arr + (yOffset>>2) * CLASSIFICATION_BLK_SIZE + (xOffset>>2)), c1);
        _mm_storel_epi64((__m128i *) (class_idx_arr + ((yOffset>>2)+1) * CLASSIFICATION_BLK_SIZE + (xOffset>>2)), _mm_bsrli_si128(c1, 8));


        _mm_storel_epi64((__m128i *) (transpose_idx_arr + (yOffset>>2) * CLASSIFICATION_BLK_SIZE + (xOffset>>2)), t1);
        _mm_storel_epi64((__m128i *) (transpose_idx_arr + ((yOffset>>2)+1) * CLASSIFICATION_BLK_SIZE + (xOffset>>2)), _mm_bsrli_si128(t1, 8));
    }
}

void rcn_init_alf_functions_sse(struct RCNFunctions *rcn_func){
  rcn_func->alf.classif=&simdDeriveClassificationBlk;
  rcn_func->alf.luma[0]=&simdFilter7x7Blk;
  rcn_func->alf.luma[1]=&simdFilter7x7BlkVB;
  rcn_func->alf.chroma[0]=&simdFilter5x5Blk;
  rcn_func->alf.chroma[1]=&simdFilter5x5BlkVB;
  rcn_func->alf.ccalf[0]=&cc_alf_filterBlkVB_sse;
  rcn_func->alf.ccalf[1]=&cc_alf_filterBlkVB_sse;
}
