#if defined _MSC_VER
#include <tmmintrin.h>
#else
#include <x86intrin.h>
#endif

#include <stdint.h>
#include "rcn_alf.h"
#include "rcn_structures.h"

#define sh(x) (0x0202 * (x & 7) + 0x0100 + 0x1010 * (x & 8))

static const uint16_t shuffleTab[4][2][8] = {
        {
                { sh(0), sh(1), sh(2), sh(3), sh(4), sh(5), sh(6), sh(7) },
                { sh(8), sh(9), sh(10), sh(11), sh(12), sh(13), sh(14), sh(15) },
        },
        {
                { sh(9), sh(4), sh(10), sh(8), sh(1), sh(5), sh(11), sh(7) },
                { sh(3), sh(0), sh(2), sh(6), sh(12), sh(13), sh(14), sh(15) },
        },
        {
                { sh(0), sh(3), sh(2), sh(1), sh(8), sh(7), sh(6), sh(5) },
                { sh(4), sh(9), sh(10), sh(11), sh(12), sh(13), sh(14), sh(15) },
        },
        {
                { sh(9), sh(8), sh(10), sh(4), sh(3), sh(7), sh(11), sh(5) },
                { sh(1), sh(0), sh(2), sh(6), sh(12), sh(13), sh(14), sh(15) },
        },
};

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


// void simdFilter5x5Blk(ALFClassifier **classifier, AVFrame *const recDst, int16_t **const recSrc,
//                              const Area blk_dst, const Area blk, const ComponentID compId,
//                              const short *filter_set, const short *clip_set,
//                              const int ctu_height, int virbnd_pos, AlfFilterType filtType, int margin)

void simdFilter5x5Blk(ALFClassifier **classifier, int16_t *const dst, int16_t *const src, const int dstStride, const int srcStride,
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

    for (size_t i = 0; i < blk_dst.height; i += STEP_Y)
    {
        for (size_t j = 0; j < blk_dst.width; j += STEP_X)
        {
            for (size_t ii = 0; ii < STEP_Y; ii++)
            {
                const uint16_t *pImg0, *pImg1, *pImg2, *pImg3, *pImg4;

                pImg0 = _src + j + ii * srcStride ;
                pImg1 = pImg0 + srcStride;
                pImg2 = pImg0 - srcStride;
                pImg3 = pImg1 + srcStride;
                pImg4 = pImg2 - srcStride;

                const int yVb = (blk_dst.y + i + ii) & (ctu_height - 1);
                if (yVb < virbnd_pos && (yVb >= virbnd_pos - 2))   // above
                {
                    pImg1 = (yVb == virbnd_pos - 1) ? pImg0 : pImg1;
                    pImg3 = (yVb >= virbnd_pos - 2) ? pImg1 : pImg3;

                    pImg2 = (yVb == virbnd_pos - 1) ? pImg0 : pImg2;
                    pImg4 = (yVb >= virbnd_pos - 2) ? pImg2 : pImg4;
                }
                else if (yVb >= virbnd_pos && (yVb <= virbnd_pos + 1))   // bottom
                {
                    pImg2 = (yVb == virbnd_pos) ? pImg0 : pImg2;
                    pImg4 = (yVb <= virbnd_pos + 1) ? pImg2 : pImg4;

                    pImg1 = (yVb == virbnd_pos) ? pImg0 : pImg1;
                    pImg3 = (yVb <= virbnd_pos + 1) ? pImg1 : pImg3;
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

                process2coeffs5x5(0, pImg3 + 0, pImg4 + 0, pImg1 + 1, pImg2 - 1);
                process2coeffs5x5(1, pImg1 + 0, pImg2 + 0, pImg1 - 1, pImg2 + 1);
                process2coeffs5x5(2, pImg0 + 2, pImg0 - 2, pImg0 + 1, pImg0 - 1);

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

                if (j + STEP_X <= blk_dst.width)
                {
                    // if (j>=x0 && i+ii>=y0 && j<=w_max && i+ii<=h_max)
                        _mm_storeu_si128((__m128i *) (_dst + ii * dstStride + j), accumA);
                }
                else
                {
                    // if (j>=x0 && i+ii>=y0 && j<=w_max && i+ii<=h_max)
                        _mm_storel_epi64((__m128i *) (_dst + ii * dstStride + j), accumA);
                }
            }
        }

        _src += srcStride * STEP_Y;
        _dst += dstStride * STEP_Y;
    }
}

// void simdFilter7x7Blk(ALFClassifier **classifier, AVFrame *const recDst, int16_t **const recSrc,
//                              const Area blk_dst, const Area blk, const ComponentID compId,
//                              const short *filter_set, const short *clip_set,
//                              const int ctu_height, int virbnd_pos, AlfFilterType filtType, int margin)
void simdFilter7x7Blk(ALFClassifier **classifier, int16_t *const dst, int16_t *const src, const int dstStride, const int srcStride,
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

    const __m128i mmOffset = _mm_set1_epi32(ROUND);
    const __m128i mmMin = _mm_set1_epi16( clpMin );
    const __m128i mmMax = _mm_set1_epi16( clpMax );

    const __m128i mmOffsetborder = _mm_set1_epi32(1 << ((SHIFT + 3) - 1));
    const int SHIFTborder = SHIFT+3;

    for (size_t i = 0; i < blk_dst.height; i += STEP_Y)
    {
        // const ALFClassifier *pClass = classifier[blk_dst.y + i] + blk_dst.x;
        const ALFClassifier *pClass = classifier[i/4];

        for (size_t j = 0; j < blk_dst.width; j += STEP_X)
        {
            __m128i params[2][2][6];

            for (int k = 0; k < 2; ++k)
            {
                const ALFClassifier cl = pClass[j/4 * k];

                const int transpose_idx = cl.transpose_idx;
                const int class_idx     = cl.class_idx;

                __m128i rawCoeff0, rawCoeff1;
                __m128i rawClip0, rawClip1;

                rawCoeff0 = _mm_loadu_si128((const __m128i *) (filter_set + class_idx * MAX_NUM_ALF_LUMA_COEFF));
                rawCoeff1 = _mm_loadl_epi64((const __m128i *) (filter_set + class_idx * MAX_NUM_ALF_LUMA_COEFF + 8));

                rawClip0 = _mm_loadu_si128((const __m128i *) (clip_set + class_idx * MAX_NUM_ALF_LUMA_COEFF));
                rawClip1 = _mm_loadl_epi64((const __m128i *) (clip_set + class_idx * MAX_NUM_ALF_LUMA_COEFF + 8));

                const __m128i s0 = _mm_loadu_si128((const __m128i *) shuffleTab[transpose_idx][0]);
                const __m128i s1 = _mm_xor_si128(s0, _mm_set1_epi8((char) 0x80));
                const __m128i s2 = _mm_loadu_si128((const __m128i *) shuffleTab[transpose_idx][1]);
                const __m128i s3 = _mm_xor_si128(s2, _mm_set1_epi8((char) 0x80));

                const __m128i rawCoeffLo = _mm_or_si128(_mm_shuffle_epi8(rawCoeff0, s0), _mm_shuffle_epi8(rawCoeff1, s1));
                const __m128i rawCoeffHi = _mm_or_si128(_mm_shuffle_epi8(rawCoeff0, s2), _mm_shuffle_epi8(rawCoeff1, s3));
                const __m128i rawClipLo  = _mm_or_si128(_mm_shuffle_epi8(rawClip0, s0), _mm_shuffle_epi8(rawClip1, s1));
                const __m128i rawClipHi  = _mm_or_si128(_mm_shuffle_epi8(rawClip0, s2), _mm_shuffle_epi8(rawClip1, s3));

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

                pImg0 = _src + j + ii * srcStride;
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
                int16_t* pRec1 = _dst + ii * dstStride + j;
            }
        }

        _src += srcStride * STEP_Y;
        _dst += dstStride * STEP_Y;
    }
}

void rcn_init_alf_functions_sse(struct RCNFunctions *rcn_func){
  rcn_func->alf.luma=&simdFilter7x7Blk;
  rcn_func->alf.chroma=&simdFilter5x5Blk;
}
