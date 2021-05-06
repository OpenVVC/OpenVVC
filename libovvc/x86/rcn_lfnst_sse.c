#include <emmintrin.h>

#include "rcn_structures.h"

static void
lfnst_4x4_16(const int16_t* const src, __m128i* r,
             const int8_t* const lfnst_matrix)
             {
               int32_t* const _src = (int32_t*) src;
                    __m128i x[8], c[16], l[32], d[32];
                    x[0] = _mm_set1_epi32(_src[0]);
                    x[1] = _mm_set1_epi32(_src[1]);
                    x[2] = _mm_set1_epi32(_src[2]);
                    x[3] = _mm_set1_epi32(_src[3]);
                    x[4] = _mm_set1_epi32(_src[4]);
                    x[5] = _mm_set1_epi32(_src[5]);
                    x[6] = _mm_set1_epi32(_src[6]);
                    x[7] = _mm_set1_epi32(_src[7]);

                    c[0]  = _mm_loadu_si128((__m128i *)&lfnst_matrix[0 * 16]);
                    c[1]  = _mm_loadu_si128((__m128i *)&lfnst_matrix[1 * 16]);
                    c[2]  = _mm_loadu_si128((__m128i *)&lfnst_matrix[2 * 16]);
                    c[3]  = _mm_loadu_si128((__m128i *)&lfnst_matrix[3 * 16]);
                    c[4]  = _mm_loadu_si128((__m128i *)&lfnst_matrix[4 * 16]);
                    c[5]  = _mm_loadu_si128((__m128i *)&lfnst_matrix[5 * 16]);
                    c[6]  = _mm_loadu_si128((__m128i *)&lfnst_matrix[6 * 16]);
                    c[7]  = _mm_loadu_si128((__m128i *)&lfnst_matrix[7 * 16]);
                    c[8]  = _mm_loadu_si128((__m128i *)&lfnst_matrix[8 * 16]);
                    c[9]  = _mm_loadu_si128((__m128i *)&lfnst_matrix[9 * 16]);
                    c[10] = _mm_loadu_si128((__m128i *)&lfnst_matrix[10 * 16]);
                    c[11] = _mm_loadu_si128((__m128i *)&lfnst_matrix[11 * 16]);
                    c[12] = _mm_loadu_si128((__m128i *)&lfnst_matrix[12 * 16]);
                    c[13] = _mm_loadu_si128((__m128i *)&lfnst_matrix[13 * 16]);
                    c[14] = _mm_loadu_si128((__m128i *)&lfnst_matrix[14 * 16]);
                    c[15] = _mm_loadu_si128((__m128i *)&lfnst_matrix[15 * 16]);


                    l[0]  = _mm_unpacklo_epi8(c[0] , _mm_cmplt_epi8(c[0] , _mm_setzero_si128()));
                    l[1]  = _mm_unpacklo_epi8(c[1] , _mm_cmplt_epi8(c[1] , _mm_setzero_si128()));
                    l[2]  = _mm_unpacklo_epi8(c[2] , _mm_cmplt_epi8(c[2] , _mm_setzero_si128()));
                    l[3]  = _mm_unpacklo_epi8(c[3] , _mm_cmplt_epi8(c[3] , _mm_setzero_si128()));
                    l[4]  = _mm_unpacklo_epi8(c[4] , _mm_cmplt_epi8(c[4] , _mm_setzero_si128()));
                    l[5]  = _mm_unpacklo_epi8(c[5] , _mm_cmplt_epi8(c[5] , _mm_setzero_si128()));
                    l[6]  = _mm_unpacklo_epi8(c[6] , _mm_cmplt_epi8(c[6] , _mm_setzero_si128()));
                    l[7]  = _mm_unpacklo_epi8(c[7] , _mm_cmplt_epi8(c[7] , _mm_setzero_si128()));
                    l[8]  = _mm_unpacklo_epi8(c[8] , _mm_cmplt_epi8(c[8] , _mm_setzero_si128()));
                    l[9]  = _mm_unpacklo_epi8(c[9] , _mm_cmplt_epi8(c[9] , _mm_setzero_si128()));
                    l[10] = _mm_unpacklo_epi8(c[10], _mm_cmplt_epi8(c[10], _mm_setzero_si128()));
                    l[11] = _mm_unpacklo_epi8(c[11], _mm_cmplt_epi8(c[11], _mm_setzero_si128()));
                    l[12] = _mm_unpacklo_epi8(c[12], _mm_cmplt_epi8(c[12], _mm_setzero_si128()));
                    l[13] = _mm_unpacklo_epi8(c[13], _mm_cmplt_epi8(c[13], _mm_setzero_si128()));
                    l[14] = _mm_unpacklo_epi8(c[14], _mm_cmplt_epi8(c[14], _mm_setzero_si128()));
                    l[15] = _mm_unpacklo_epi8(c[15], _mm_cmplt_epi8(c[15], _mm_setzero_si128()));

                    l[16] = _mm_unpackhi_epi8(c[0] , _mm_cmplt_epi8(c[0] , _mm_setzero_si128()));
                    l[17] = _mm_unpackhi_epi8(c[1] , _mm_cmplt_epi8(c[1] , _mm_setzero_si128()));
                    l[18] = _mm_unpackhi_epi8(c[2] , _mm_cmplt_epi8(c[2] , _mm_setzero_si128()));
                    l[19] = _mm_unpackhi_epi8(c[3] , _mm_cmplt_epi8(c[3] , _mm_setzero_si128()));
                    l[20] = _mm_unpackhi_epi8(c[4] , _mm_cmplt_epi8(c[4] , _mm_setzero_si128()));
                    l[21] = _mm_unpackhi_epi8(c[5] , _mm_cmplt_epi8(c[5] , _mm_setzero_si128()));
                    l[22] = _mm_unpackhi_epi8(c[6] , _mm_cmplt_epi8(c[6] , _mm_setzero_si128()));
                    l[23] = _mm_unpackhi_epi8(c[7] , _mm_cmplt_epi8(c[7] , _mm_setzero_si128()));
                    l[24] = _mm_unpackhi_epi8(c[8] , _mm_cmplt_epi8(c[8] , _mm_setzero_si128()));
                    l[25] = _mm_unpackhi_epi8(c[9] , _mm_cmplt_epi8(c[9] , _mm_setzero_si128()));
                    l[26] = _mm_unpackhi_epi8(c[10], _mm_cmplt_epi8(c[10], _mm_setzero_si128()));
                    l[27] = _mm_unpackhi_epi8(c[11], _mm_cmplt_epi8(c[11], _mm_setzero_si128()));
                    l[28] = _mm_unpackhi_epi8(c[12], _mm_cmplt_epi8(c[12], _mm_setzero_si128()));
                    l[29] = _mm_unpackhi_epi8(c[13], _mm_cmplt_epi8(c[13], _mm_setzero_si128()));
                    l[30] = _mm_unpackhi_epi8(c[14], _mm_cmplt_epi8(c[14], _mm_setzero_si128()));
                    l[31] = _mm_unpackhi_epi8(c[15], _mm_cmplt_epi8(c[15], _mm_setzero_si128()));

                    d[0]  = _mm_unpacklo_epi16(l[0], l[1]);
                    d[1]  = _mm_unpackhi_epi16(l[0], l[1]);
                    d[2]  = _mm_unpacklo_epi16(l[2], l[3]);
                    d[3]  = _mm_unpackhi_epi16(l[2], l[3]);
                    d[4]  = _mm_unpacklo_epi16(l[4], l[5]);
                    d[5]  = _mm_unpackhi_epi16(l[4], l[5]);
                    d[6]  = _mm_unpacklo_epi16(l[6], l[7]);
                    d[7]  = _mm_unpackhi_epi16(l[6], l[7]);
                    d[8]  = _mm_unpacklo_epi16(l[8], l[9]);
                    d[9]  = _mm_unpackhi_epi16(l[8], l[9]);
                    d[10] = _mm_unpacklo_epi16(l[10], l[11]);
                    d[11] = _mm_unpackhi_epi16(l[10], l[11]);
                    d[12] = _mm_unpacklo_epi16(l[12], l[13]);
                    d[13] = _mm_unpackhi_epi16(l[12], l[13]);
                    d[14] = _mm_unpacklo_epi16(l[14], l[15]);
                    d[15] = _mm_unpackhi_epi16(l[14], l[15]);

                    d[16] = _mm_unpacklo_epi16(l[16], l[17]);
                    d[17] = _mm_unpackhi_epi16(l[16], l[17]);
                    d[18] = _mm_unpacklo_epi16(l[18], l[19]);
                    d[19] = _mm_unpackhi_epi16(l[18], l[19]);
                    d[20] = _mm_unpacklo_epi16(l[20], l[21]);
                    d[21] = _mm_unpackhi_epi16(l[20], l[21]);
                    d[22] = _mm_unpacklo_epi16(l[22], l[23]);
                    d[23] = _mm_unpackhi_epi16(l[22], l[23]);
                    d[24] = _mm_unpacklo_epi16(l[24], l[25]);
                    d[25] = _mm_unpackhi_epi16(l[24], l[25]);
                    d[26] = _mm_unpacklo_epi16(l[26], l[27]);
                    d[27] = _mm_unpackhi_epi16(l[26], l[27]);
                    d[28] = _mm_unpacklo_epi16(l[28], l[29]);
                    d[29] = _mm_unpackhi_epi16(l[28], l[29]);
                    d[30] = _mm_unpacklo_epi16(l[30], l[31]);
                    d[31] = _mm_unpackhi_epi16(l[30], l[31]);

                    l[0]  = _mm_madd_epi16(x[0],d[0]);
                    l[1]  = _mm_madd_epi16(x[0],d[1]);
                    l[2]  = _mm_madd_epi16(x[1],d[2]);
                    l[3]  = _mm_madd_epi16(x[1],d[3]);
                    l[4]  = _mm_madd_epi16(x[2],d[4]);
                    l[5]  = _mm_madd_epi16(x[2],d[5]);
                    l[6]  = _mm_madd_epi16(x[3],d[6]);
                    l[7]  = _mm_madd_epi16(x[3],d[7]);
                    l[8]  = _mm_madd_epi16(x[4],d[8]);
                    l[9]  = _mm_madd_epi16(x[4],d[9]);
                    l[10] = _mm_madd_epi16(x[5],d[10]);
                    l[11] = _mm_madd_epi16(x[5],d[11]);
                    l[12] = _mm_madd_epi16(x[6],d[12]);
                    l[13] = _mm_madd_epi16(x[6],d[13]);
                    l[14] = _mm_madd_epi16(x[7],d[14]);
                    l[15] = _mm_madd_epi16(x[7],d[15]);

                    l[16] = _mm_madd_epi16(x[0],d[16]);
                    l[17] = _mm_madd_epi16(x[0],d[17]);
                    l[18] = _mm_madd_epi16(x[1],d[18]);
                    l[19] = _mm_madd_epi16(x[1],d[19]);
                    l[20] = _mm_madd_epi16(x[2],d[20]);
                    l[21] = _mm_madd_epi16(x[2],d[21]);
                    l[22] = _mm_madd_epi16(x[3],d[22]);
                    l[23] = _mm_madd_epi16(x[3],d[23]);
                    l[24] = _mm_madd_epi16(x[4],d[24]);
                    l[25] = _mm_madd_epi16(x[4],d[25]);
                    l[26] = _mm_madd_epi16(x[5],d[26]);
                    l[27] = _mm_madd_epi16(x[5],d[27]);
                    l[28] = _mm_madd_epi16(x[6],d[28]);
                    l[29] = _mm_madd_epi16(x[6],d[29]);
                    l[30] = _mm_madd_epi16(x[7],d[30]);
                    l[31] = _mm_madd_epi16(x[7],d[31]);

                    d[0]  = _mm_add_epi32(l[0],l[8]);
                    d[1]  = _mm_add_epi32(l[1],l[9]);
                    d[2]  = _mm_add_epi32(l[2],l[10]);
                    d[3]  = _mm_add_epi32(l[3],l[11]);
                    d[4]  = _mm_add_epi32(l[4],l[12]);
                    d[5]  = _mm_add_epi32(l[5],l[13]);
                    d[6]  = _mm_add_epi32(l[6],l[14]);
                    d[7]  = _mm_add_epi32(l[7],l[15]);

                    d[8]  = _mm_add_epi32(l[16],l[24]);
                    d[9]  = _mm_add_epi32(l[17],l[25]);
                    d[10] = _mm_add_epi32(l[18],l[26]);
                    d[11] = _mm_add_epi32(l[19],l[27]);
                    d[12] = _mm_add_epi32(l[20],l[28]);
                    d[13] = _mm_add_epi32(l[21],l[29]);
                    d[14] = _mm_add_epi32(l[22],l[30]);
                    d[15] = _mm_add_epi32(l[23],l[31]);

                    l[0]  = _mm_add_epi32(d[0],d[4]);
                    l[1]  = _mm_add_epi32(d[1],d[5]);
                    l[2]  = _mm_add_epi32(d[2],d[6]);
                    l[3]  = _mm_add_epi32(d[3],d[7]);

                    l[4]  = _mm_add_epi32(d[8],d[12]);
                    l[5]  = _mm_add_epi32(d[9],d[13]);
                    l[6]  = _mm_add_epi32(d[10],d[14]);
                    l[7]  = _mm_add_epi32(d[11],d[15]);

                    d[0]  = _mm_add_epi32(l[0],l[2]);
                    d[1]  = _mm_add_epi32(l[1],l[3]);

                    d[2]  = _mm_add_epi32(l[4],l[6]);
                    d[3]  = _mm_add_epi32(l[5],l[7]);


                    d[0]  = _mm_add_epi32(d[0],_mm_set1_epi32(64));
                    d[1]  = _mm_add_epi32(d[1],_mm_set1_epi32(64));
                    d[2]  = _mm_add_epi32(d[2],_mm_set1_epi32(64));
                    d[3]  = _mm_add_epi32(d[3],_mm_set1_epi32(64));

                    d[0]  = _mm_srai_epi32(d[0],7);
                    d[1]  = _mm_srai_epi32(d[1],7);
                    d[2]  = _mm_srai_epi32(d[2],7);
                    d[3]  = _mm_srai_epi32(d[3],7);

                    r[0]  = _mm_packs_epi32(d[0],d[1]);
                    r[1]  = _mm_packs_epi32(d[2],d[3]);
                  }

static void
lfnst_4x4_8(const int16_t* const src, __m128i* r,
             const int8_t* const lfnst_matrix)
             {
               int32_t* const _src = (int32_t*) src;
               __m128i x[8], c[16], l[16], d[16];
               x[0] = _mm_set1_epi32(_src[0]);
               x[1] = _mm_set1_epi32(_src[1]);
               x[2] = _mm_set1_epi32(_src[2]);
               x[3] = _mm_set1_epi32(_src[3]);


               c[0]  = _mm_loadu_si128((__m128i *)&lfnst_matrix[0 * 16]);
               c[1]  = _mm_loadu_si128((__m128i *)&lfnst_matrix[1 * 16]);
               c[2]  = _mm_loadu_si128((__m128i *)&lfnst_matrix[2 * 16]);
               c[3]  = _mm_loadu_si128((__m128i *)&lfnst_matrix[3 * 16]);
               c[4]  = _mm_loadu_si128((__m128i *)&lfnst_matrix[4 * 16]);
               c[5]  = _mm_loadu_si128((__m128i *)&lfnst_matrix[5 * 16]);
               c[6]  = _mm_loadu_si128((__m128i *)&lfnst_matrix[6 * 16]);
               c[7]  = _mm_loadu_si128((__m128i *)&lfnst_matrix[7 * 16]);


               l[0]  = _mm_unpacklo_epi8(c[0], _mm_cmplt_epi8(c[0], _mm_setzero_si128()));
               l[1]  = _mm_unpacklo_epi8(c[1], _mm_cmplt_epi8(c[1], _mm_setzero_si128()));
               l[2]  = _mm_unpacklo_epi8(c[2], _mm_cmplt_epi8(c[2], _mm_setzero_si128()));
               l[3]  = _mm_unpacklo_epi8(c[3], _mm_cmplt_epi8(c[3], _mm_setzero_si128()));
               l[4]  = _mm_unpacklo_epi8(c[4], _mm_cmplt_epi8(c[4], _mm_setzero_si128()));
               l[5]  = _mm_unpacklo_epi8(c[5], _mm_cmplt_epi8(c[5], _mm_setzero_si128()));
               l[6]  = _mm_unpacklo_epi8(c[6], _mm_cmplt_epi8(c[6], _mm_setzero_si128()));
               l[7]  = _mm_unpacklo_epi8(c[7], _mm_cmplt_epi8(c[7], _mm_setzero_si128()));

               l[8]  = _mm_unpackhi_epi8(c[0], _mm_cmplt_epi8(c[0], _mm_setzero_si128()));
               l[9]  = _mm_unpackhi_epi8(c[1], _mm_cmplt_epi8(c[1], _mm_setzero_si128()));
               l[10] = _mm_unpackhi_epi8(c[2], _mm_cmplt_epi8(c[2], _mm_setzero_si128()));
               l[11] = _mm_unpackhi_epi8(c[3], _mm_cmplt_epi8(c[3], _mm_setzero_si128()));
               l[12] = _mm_unpackhi_epi8(c[4], _mm_cmplt_epi8(c[4], _mm_setzero_si128()));
               l[13] = _mm_unpackhi_epi8(c[5], _mm_cmplt_epi8(c[5], _mm_setzero_si128()));
               l[14] = _mm_unpackhi_epi8(c[6], _mm_cmplt_epi8(c[6], _mm_setzero_si128()));
               l[15] = _mm_unpackhi_epi8(c[7], _mm_cmplt_epi8(c[7], _mm_setzero_si128()));

               d[0]  = _mm_unpacklo_epi16(l[0], l[1]);
               d[1]  = _mm_unpackhi_epi16(l[0], l[1]);
               d[2]  = _mm_unpacklo_epi16(l[2], l[3]);
               d[3]  = _mm_unpackhi_epi16(l[2], l[3]);
               d[4]  = _mm_unpacklo_epi16(l[4], l[5]);
               d[5]  = _mm_unpackhi_epi16(l[4], l[5]);
               d[6]  = _mm_unpacklo_epi16(l[6], l[7]);
               d[7]  = _mm_unpackhi_epi16(l[6], l[7]);

               d[8]  = _mm_unpacklo_epi16(l[8], l[9]);
               d[9]  = _mm_unpackhi_epi16(l[8], l[9]);
               d[10] = _mm_unpacklo_epi16(l[10], l[11]);
               d[11] = _mm_unpackhi_epi16(l[10], l[11]);
               d[12] = _mm_unpacklo_epi16(l[12], l[13]);
               d[13] = _mm_unpackhi_epi16(l[12], l[13]);
               d[14] = _mm_unpacklo_epi16(l[14], l[15]);
               d[15] = _mm_unpackhi_epi16(l[14], l[15]);

               l[0]  = _mm_madd_epi16(x[0],d[0]);
               l[1]  = _mm_madd_epi16(x[0],d[1]);
               l[2]  = _mm_madd_epi16(x[1],d[2]);
               l[3]  = _mm_madd_epi16(x[1],d[3]);
               l[4]  = _mm_madd_epi16(x[2],d[4]);
               l[5]  = _mm_madd_epi16(x[2],d[5]);
               l[6]  = _mm_madd_epi16(x[3],d[6]);
               l[7]  = _mm_madd_epi16(x[3],d[7]);

               l[8]  = _mm_madd_epi16(x[0],d[8]);
               l[9]  = _mm_madd_epi16(x[0],d[9]);
               l[10] = _mm_madd_epi16(x[1],d[10]);
               l[11] = _mm_madd_epi16(x[1],d[11]);
               l[12] = _mm_madd_epi16(x[2],d[12]);
               l[13] = _mm_madd_epi16(x[2],d[13]);
               l[14] = _mm_madd_epi16(x[3],d[14]);
               l[15] = _mm_madd_epi16(x[3],d[15]);

               d[0]  = _mm_add_epi32(l[0],l[4]);
               d[1]  = _mm_add_epi32(l[1],l[5]);
               d[2]  = _mm_add_epi32(l[2],l[6]);
               d[3]  = _mm_add_epi32(l[3],l[7]);

               d[4]  = _mm_add_epi32(l[8],l[12]);
               d[5]  = _mm_add_epi32(l[9],l[13]);
               d[6]  = _mm_add_epi32(l[10],l[14]);
               d[7]  = _mm_add_epi32(l[11],l[15]);

               l[0]  = _mm_add_epi32(d[0],d[2]);
               l[1]  = _mm_add_epi32(d[1],d[3]);

               l[2]  = _mm_add_epi32(d[4],d[6]);
               l[3]  = _mm_add_epi32(d[5],d[7]);

               l[0]  = _mm_add_epi32(l[0],_mm_set1_epi32(64));
               l[1]  = _mm_add_epi32(l[1],_mm_set1_epi32(64));
               l[2]  = _mm_add_epi32(l[2],_mm_set1_epi32(64));
               l[3]  = _mm_add_epi32(l[3],_mm_set1_epi32(64));

               l[0]  = _mm_srai_epi32(l[0],7);
               l[1]  = _mm_srai_epi32(l[1],7);
               l[2]  = _mm_srai_epi32(l[2],7);
               l[3]  = _mm_srai_epi32(l[3],7);

               r[0]  = _mm_packs_epi32(l[0],l[1]);
               r[1]  = _mm_packs_epi32(l[2],l[3]);
             }

static void
lfnst_8x8_16(const int16_t* const src, __m128i* r,
            const int8_t* const lfnst_matrix)
            {
              int32_t* const _src = (int32_t*) src;
                   __m128i x[8], c[16], l[32], d[32];
                   x[0] = _mm_set1_epi32(_src[0]);
                   x[1] = _mm_set1_epi32(_src[1]);
                   x[2] = _mm_set1_epi32(_src[2]);
                   x[3] = _mm_set1_epi32(_src[3]);
                   x[4] = _mm_set1_epi32(_src[4]);
                   x[5] = _mm_set1_epi32(_src[5]);
                   x[6] = _mm_set1_epi32(_src[6]);
                   x[7] = _mm_set1_epi32(_src[7]);

                   c[0]  = _mm_loadu_si128((__m128i *)&lfnst_matrix[0 * 48]);
                   c[1]  = _mm_loadu_si128((__m128i *)&lfnst_matrix[1 * 48]);
                   c[2]  = _mm_loadu_si128((__m128i *)&lfnst_matrix[2 * 48]);
                   c[3]  = _mm_loadu_si128((__m128i *)&lfnst_matrix[3 * 48]);
                   c[4]  = _mm_loadu_si128((__m128i *)&lfnst_matrix[4 * 48]);
                   c[5]  = _mm_loadu_si128((__m128i *)&lfnst_matrix[5 * 48]);
                   c[6]  = _mm_loadu_si128((__m128i *)&lfnst_matrix[6 * 48]);
                   c[7]  = _mm_loadu_si128((__m128i *)&lfnst_matrix[7 * 48]);
                   c[8]  = _mm_loadu_si128((__m128i *)&lfnst_matrix[8 * 48]);
                   c[9]  = _mm_loadu_si128((__m128i *)&lfnst_matrix[9 * 48]);
                   c[10] = _mm_loadu_si128((__m128i *)&lfnst_matrix[10 * 48]);
                   c[11] = _mm_loadu_si128((__m128i *)&lfnst_matrix[11 * 48]);
                   c[12] = _mm_loadu_si128((__m128i *)&lfnst_matrix[12 * 48]);
                   c[13] = _mm_loadu_si128((__m128i *)&lfnst_matrix[13 * 48]);
                   c[14] = _mm_loadu_si128((__m128i *)&lfnst_matrix[14 * 48]);
                   c[15] = _mm_loadu_si128((__m128i *)&lfnst_matrix[15 * 48]);


                   l[0]  = _mm_unpacklo_epi8(c[0] , _mm_cmplt_epi8(c[0] , _mm_setzero_si128()));
                   l[1]  = _mm_unpacklo_epi8(c[1] , _mm_cmplt_epi8(c[1] , _mm_setzero_si128()));
                   l[2]  = _mm_unpacklo_epi8(c[2] , _mm_cmplt_epi8(c[2] , _mm_setzero_si128()));
                   l[3]  = _mm_unpacklo_epi8(c[3] , _mm_cmplt_epi8(c[3] , _mm_setzero_si128()));
                   l[4]  = _mm_unpacklo_epi8(c[4] , _mm_cmplt_epi8(c[4] , _mm_setzero_si128()));
                   l[5]  = _mm_unpacklo_epi8(c[5] , _mm_cmplt_epi8(c[5] , _mm_setzero_si128()));
                   l[6]  = _mm_unpacklo_epi8(c[6] , _mm_cmplt_epi8(c[6] , _mm_setzero_si128()));
                   l[7]  = _mm_unpacklo_epi8(c[7] , _mm_cmplt_epi8(c[7] , _mm_setzero_si128()));
                   l[8]  = _mm_unpacklo_epi8(c[8] , _mm_cmplt_epi8(c[8] , _mm_setzero_si128()));
                   l[9]  = _mm_unpacklo_epi8(c[9] , _mm_cmplt_epi8(c[9] , _mm_setzero_si128()));
                   l[10] = _mm_unpacklo_epi8(c[10], _mm_cmplt_epi8(c[10], _mm_setzero_si128()));
                   l[11] = _mm_unpacklo_epi8(c[11], _mm_cmplt_epi8(c[11], _mm_setzero_si128()));
                   l[12] = _mm_unpacklo_epi8(c[12], _mm_cmplt_epi8(c[12], _mm_setzero_si128()));
                   l[13] = _mm_unpacklo_epi8(c[13], _mm_cmplt_epi8(c[13], _mm_setzero_si128()));
                   l[14] = _mm_unpacklo_epi8(c[14], _mm_cmplt_epi8(c[14], _mm_setzero_si128()));
                   l[15] = _mm_unpacklo_epi8(c[15], _mm_cmplt_epi8(c[15], _mm_setzero_si128()));

                   l[16] = _mm_unpackhi_epi8(c[0] , _mm_cmplt_epi8(c[0] , _mm_setzero_si128()));
                   l[17] = _mm_unpackhi_epi8(c[1] , _mm_cmplt_epi8(c[1] , _mm_setzero_si128()));
                   l[18] = _mm_unpackhi_epi8(c[2] , _mm_cmplt_epi8(c[2] , _mm_setzero_si128()));
                   l[19] = _mm_unpackhi_epi8(c[3] , _mm_cmplt_epi8(c[3] , _mm_setzero_si128()));
                   l[20] = _mm_unpackhi_epi8(c[4] , _mm_cmplt_epi8(c[4] , _mm_setzero_si128()));
                   l[21] = _mm_unpackhi_epi8(c[5] , _mm_cmplt_epi8(c[5] , _mm_setzero_si128()));
                   l[22] = _mm_unpackhi_epi8(c[6] , _mm_cmplt_epi8(c[6] , _mm_setzero_si128()));
                   l[23] = _mm_unpackhi_epi8(c[7] , _mm_cmplt_epi8(c[7] , _mm_setzero_si128()));
                   l[24] = _mm_unpackhi_epi8(c[8] , _mm_cmplt_epi8(c[8] , _mm_setzero_si128()));
                   l[25] = _mm_unpackhi_epi8(c[9] , _mm_cmplt_epi8(c[9] , _mm_setzero_si128()));
                   l[26] = _mm_unpackhi_epi8(c[10], _mm_cmplt_epi8(c[10], _mm_setzero_si128()));
                   l[27] = _mm_unpackhi_epi8(c[11], _mm_cmplt_epi8(c[11], _mm_setzero_si128()));
                   l[28] = _mm_unpackhi_epi8(c[12], _mm_cmplt_epi8(c[12], _mm_setzero_si128()));
                   l[29] = _mm_unpackhi_epi8(c[13], _mm_cmplt_epi8(c[13], _mm_setzero_si128()));
                   l[30] = _mm_unpackhi_epi8(c[14], _mm_cmplt_epi8(c[14], _mm_setzero_si128()));
                   l[31] = _mm_unpackhi_epi8(c[15], _mm_cmplt_epi8(c[15], _mm_setzero_si128()));

                   d[0]  = _mm_unpacklo_epi16(l[0], l[1]);
                   d[1]  = _mm_unpackhi_epi16(l[0], l[1]);
                   d[2]  = _mm_unpacklo_epi16(l[2], l[3]);
                   d[3]  = _mm_unpackhi_epi16(l[2], l[3]);
                   d[4]  = _mm_unpacklo_epi16(l[4], l[5]);
                   d[5]  = _mm_unpackhi_epi16(l[4], l[5]);
                   d[6]  = _mm_unpacklo_epi16(l[6], l[7]);
                   d[7]  = _mm_unpackhi_epi16(l[6], l[7]);
                   d[8]  = _mm_unpacklo_epi16(l[8], l[9]);
                   d[9]  = _mm_unpackhi_epi16(l[8], l[9]);
                   d[10] = _mm_unpacklo_epi16(l[10], l[11]);
                   d[11] = _mm_unpackhi_epi16(l[10], l[11]);
                   d[12] = _mm_unpacklo_epi16(l[12], l[13]);
                   d[13] = _mm_unpackhi_epi16(l[12], l[13]);
                   d[14] = _mm_unpacklo_epi16(l[14], l[15]);
                   d[15] = _mm_unpackhi_epi16(l[14], l[15]);

                   d[16] = _mm_unpacklo_epi16(l[16], l[17]);
                   d[17] = _mm_unpackhi_epi16(l[16], l[17]);
                   d[18] = _mm_unpacklo_epi16(l[18], l[19]);
                   d[19] = _mm_unpackhi_epi16(l[18], l[19]);
                   d[20] = _mm_unpacklo_epi16(l[20], l[21]);
                   d[21] = _mm_unpackhi_epi16(l[20], l[21]);
                   d[22] = _mm_unpacklo_epi16(l[22], l[23]);
                   d[23] = _mm_unpackhi_epi16(l[22], l[23]);
                   d[24] = _mm_unpacklo_epi16(l[24], l[25]);
                   d[25] = _mm_unpackhi_epi16(l[24], l[25]);
                   d[26] = _mm_unpacklo_epi16(l[26], l[27]);
                   d[27] = _mm_unpackhi_epi16(l[26], l[27]);
                   d[28] = _mm_unpacklo_epi16(l[28], l[29]);
                   d[29] = _mm_unpackhi_epi16(l[28], l[29]);
                   d[30] = _mm_unpacklo_epi16(l[30], l[31]);
                   d[31] = _mm_unpackhi_epi16(l[30], l[31]);

                   l[0]  = _mm_madd_epi16(x[0],d[0]);
                   l[1]  = _mm_madd_epi16(x[0],d[1]);
                   l[2]  = _mm_madd_epi16(x[1],d[2]);
                   l[3]  = _mm_madd_epi16(x[1],d[3]);
                   l[4]  = _mm_madd_epi16(x[2],d[4]);
                   l[5]  = _mm_madd_epi16(x[2],d[5]);
                   l[6]  = _mm_madd_epi16(x[3],d[6]);
                   l[7]  = _mm_madd_epi16(x[3],d[7]);
                   l[8]  = _mm_madd_epi16(x[4],d[8]);
                   l[9]  = _mm_madd_epi16(x[4],d[9]);
                   l[10] = _mm_madd_epi16(x[5],d[10]);
                   l[11] = _mm_madd_epi16(x[5],d[11]);
                   l[12] = _mm_madd_epi16(x[6],d[12]);
                   l[13] = _mm_madd_epi16(x[6],d[13]);
                   l[14] = _mm_madd_epi16(x[7],d[14]);
                   l[15] = _mm_madd_epi16(x[7],d[15]);

                   l[16] = _mm_madd_epi16(x[0],d[16]);
                   l[17] = _mm_madd_epi16(x[0],d[17]);
                   l[18] = _mm_madd_epi16(x[1],d[18]);
                   l[19] = _mm_madd_epi16(x[1],d[19]);
                   l[20] = _mm_madd_epi16(x[2],d[20]);
                   l[21] = _mm_madd_epi16(x[2],d[21]);
                   l[22] = _mm_madd_epi16(x[3],d[22]);
                   l[23] = _mm_madd_epi16(x[3],d[23]);
                   l[24] = _mm_madd_epi16(x[4],d[24]);
                   l[25] = _mm_madd_epi16(x[4],d[25]);
                   l[26] = _mm_madd_epi16(x[5],d[26]);
                   l[27] = _mm_madd_epi16(x[5],d[27]);
                   l[28] = _mm_madd_epi16(x[6],d[28]);
                   l[29] = _mm_madd_epi16(x[6],d[29]);
                   l[30] = _mm_madd_epi16(x[7],d[30]);
                   l[31] = _mm_madd_epi16(x[7],d[31]);

                   d[0]  = _mm_add_epi32(l[0],l[8]);
                   d[1]  = _mm_add_epi32(l[1],l[9]);
                   d[2]  = _mm_add_epi32(l[2],l[10]);
                   d[3]  = _mm_add_epi32(l[3],l[11]);
                   d[4]  = _mm_add_epi32(l[4],l[12]);
                   d[5]  = _mm_add_epi32(l[5],l[13]);
                   d[6]  = _mm_add_epi32(l[6],l[14]);
                   d[7]  = _mm_add_epi32(l[7],l[15]);

                   d[8]  = _mm_add_epi32(l[16],l[24]);
                   d[9]  = _mm_add_epi32(l[17],l[25]);
                   d[10] = _mm_add_epi32(l[18],l[26]);
                   d[11] = _mm_add_epi32(l[19],l[27]);
                   d[12] = _mm_add_epi32(l[20],l[28]);
                   d[13] = _mm_add_epi32(l[21],l[29]);
                   d[14] = _mm_add_epi32(l[22],l[30]);
                   d[15] = _mm_add_epi32(l[23],l[31]);

                   l[0]  = _mm_add_epi32(d[0],d[4]);
                   l[1]  = _mm_add_epi32(d[1],d[5]);
                   l[2]  = _mm_add_epi32(d[2],d[6]);
                   l[3]  = _mm_add_epi32(d[3],d[7]);

                   l[4]  = _mm_add_epi32(d[8],d[12]);
                   l[5]  = _mm_add_epi32(d[9],d[13]);
                   l[6]  = _mm_add_epi32(d[10],d[14]);
                   l[7]  = _mm_add_epi32(d[11],d[15]);

                   d[0]  = _mm_add_epi32(l[0],l[2]);
                   d[1]  = _mm_add_epi32(l[1],l[3]);

                   d[2]  = _mm_add_epi32(l[4],l[6]);
                   d[3]  = _mm_add_epi32(l[5],l[7]);


                   d[0]  = _mm_add_epi32(d[0],_mm_set1_epi32(64));
                   d[1]  = _mm_add_epi32(d[1],_mm_set1_epi32(64));
                   d[2]  = _mm_add_epi32(d[2],_mm_set1_epi32(64));
                   d[3]  = _mm_add_epi32(d[3],_mm_set1_epi32(64));

                   d[0]  = _mm_srai_epi32(d[0],7);
                   d[1]  = _mm_srai_epi32(d[1],7);
                   d[2]  = _mm_srai_epi32(d[2],7);
                   d[3]  = _mm_srai_epi32(d[3],7);

                   r[0]  = _mm_packs_epi32(d[0],d[1]);
                   r[1]  = _mm_packs_epi32(d[2],d[3]);
                 }

static void
lfnst_8x8_8(const int16_t* const src, __m128i* r,
            const int8_t* const lfnst_matrix)
            {
              int32_t* const _src = (int32_t*) src;
              __m128i x[8], c[16], l[16], d[16];
              x[0] = _mm_set1_epi32(_src[0]);
              x[1] = _mm_set1_epi32(_src[1]);
              x[2] = _mm_set1_epi32(_src[2]);
              x[3] = _mm_set1_epi32(_src[3]);


              c[0]  = _mm_loadu_si128((__m128i *)&lfnst_matrix[0 * 48]);
              c[1]  = _mm_loadu_si128((__m128i *)&lfnst_matrix[1 * 48]);
              c[2]  = _mm_loadu_si128((__m128i *)&lfnst_matrix[2 * 48]);
              c[3]  = _mm_loadu_si128((__m128i *)&lfnst_matrix[3 * 48]);
              c[4]  = _mm_loadu_si128((__m128i *)&lfnst_matrix[4 * 48]);
              c[5]  = _mm_loadu_si128((__m128i *)&lfnst_matrix[5 * 48]);
              c[6]  = _mm_loadu_si128((__m128i *)&lfnst_matrix[6 * 48]);
              c[7]  = _mm_loadu_si128((__m128i *)&lfnst_matrix[7 * 48]);


              l[0]  = _mm_unpacklo_epi8(c[0], _mm_cmplt_epi8(c[0], _mm_setzero_si128()));
              l[1]  = _mm_unpacklo_epi8(c[1], _mm_cmplt_epi8(c[1], _mm_setzero_si128()));
              l[2]  = _mm_unpacklo_epi8(c[2], _mm_cmplt_epi8(c[2], _mm_setzero_si128()));
              l[3]  = _mm_unpacklo_epi8(c[3], _mm_cmplt_epi8(c[3], _mm_setzero_si128()));
              l[4]  = _mm_unpacklo_epi8(c[4], _mm_cmplt_epi8(c[4], _mm_setzero_si128()));
              l[5]  = _mm_unpacklo_epi8(c[5], _mm_cmplt_epi8(c[5], _mm_setzero_si128()));
              l[6]  = _mm_unpacklo_epi8(c[6], _mm_cmplt_epi8(c[6], _mm_setzero_si128()));
              l[7]  = _mm_unpacklo_epi8(c[7], _mm_cmplt_epi8(c[7], _mm_setzero_si128()));

              l[8]  = _mm_unpackhi_epi8(c[0], _mm_cmplt_epi8(c[0], _mm_setzero_si128()));
              l[9]  = _mm_unpackhi_epi8(c[1], _mm_cmplt_epi8(c[1], _mm_setzero_si128()));
              l[10] = _mm_unpackhi_epi8(c[2], _mm_cmplt_epi8(c[2], _mm_setzero_si128()));
              l[11] = _mm_unpackhi_epi8(c[3], _mm_cmplt_epi8(c[3], _mm_setzero_si128()));
              l[12] = _mm_unpackhi_epi8(c[4], _mm_cmplt_epi8(c[4], _mm_setzero_si128()));
              l[13] = _mm_unpackhi_epi8(c[5], _mm_cmplt_epi8(c[5], _mm_setzero_si128()));
              l[14] = _mm_unpackhi_epi8(c[6], _mm_cmplt_epi8(c[6], _mm_setzero_si128()));
              l[15] = _mm_unpackhi_epi8(c[7], _mm_cmplt_epi8(c[7], _mm_setzero_si128()));

              d[0]  = _mm_unpacklo_epi16(l[0], l[1]);
              d[1]  = _mm_unpackhi_epi16(l[0], l[1]);
              d[2]  = _mm_unpacklo_epi16(l[2], l[3]);
              d[3]  = _mm_unpackhi_epi16(l[2], l[3]);
              d[4]  = _mm_unpacklo_epi16(l[4], l[5]);
              d[5]  = _mm_unpackhi_epi16(l[4], l[5]);
              d[6]  = _mm_unpacklo_epi16(l[6], l[7]);
              d[7]  = _mm_unpackhi_epi16(l[6], l[7]);

              d[8]  = _mm_unpacklo_epi16(l[8], l[9]);
              d[9]  = _mm_unpackhi_epi16(l[8], l[9]);
              d[10] = _mm_unpacklo_epi16(l[10], l[11]);
              d[11] = _mm_unpackhi_epi16(l[10], l[11]);
              d[12] = _mm_unpacklo_epi16(l[12], l[13]);
              d[13] = _mm_unpackhi_epi16(l[12], l[13]);
              d[14] = _mm_unpacklo_epi16(l[14], l[15]);
              d[15] = _mm_unpackhi_epi16(l[14], l[15]);

              l[0]  = _mm_madd_epi16(x[0],d[0]);
              l[1]  = _mm_madd_epi16(x[0],d[1]);
              l[2]  = _mm_madd_epi16(x[1],d[2]);
              l[3]  = _mm_madd_epi16(x[1],d[3]);
              l[4]  = _mm_madd_epi16(x[2],d[4]);
              l[5]  = _mm_madd_epi16(x[2],d[5]);
              l[6]  = _mm_madd_epi16(x[3],d[6]);
              l[7]  = _mm_madd_epi16(x[3],d[7]);

              l[8]  = _mm_madd_epi16(x[0],d[8]);
              l[9]  = _mm_madd_epi16(x[0],d[9]);
              l[10] = _mm_madd_epi16(x[1],d[10]);
              l[11] = _mm_madd_epi16(x[1],d[11]);
              l[12] = _mm_madd_epi16(x[2],d[12]);
              l[13] = _mm_madd_epi16(x[2],d[13]);
              l[14] = _mm_madd_epi16(x[3],d[14]);
              l[15] = _mm_madd_epi16(x[3],d[15]);

              d[0]  = _mm_add_epi32(l[0],l[4]);
              d[1]  = _mm_add_epi32(l[1],l[5]);
              d[2]  = _mm_add_epi32(l[2],l[6]);
              d[3]  = _mm_add_epi32(l[3],l[7]);

              d[4]  = _mm_add_epi32(l[8],l[12]);
              d[5]  = _mm_add_epi32(l[9],l[13]);
              d[6]  = _mm_add_epi32(l[10],l[14]);
              d[7]  = _mm_add_epi32(l[11],l[15]);

              l[0]  = _mm_add_epi32(d[0],d[2]);
              l[1]  = _mm_add_epi32(d[1],d[3]);

              l[2]  = _mm_add_epi32(d[4],d[6]);
              l[3]  = _mm_add_epi32(d[5],d[7]);

              l[0]  = _mm_add_epi32(l[0],_mm_set1_epi32(64));
              l[1]  = _mm_add_epi32(l[1],_mm_set1_epi32(64));
              l[2]  = _mm_add_epi32(l[2],_mm_set1_epi32(64));
              l[3]  = _mm_add_epi32(l[3],_mm_set1_epi32(64));

              l[0]  = _mm_srai_epi32(l[0],7);
              l[1]  = _mm_srai_epi32(l[1],7);
              l[2]  = _mm_srai_epi32(l[2],7);
              l[3]  = _mm_srai_epi32(l[3],7);

              r[0]  = _mm_packs_epi32(l[0],l[1]);
              r[1]  = _mm_packs_epi32(l[2],l[3]);
            }

static void
compute_lfnst_4x4_sse(const int16_t* const src, int16_t* const dst,
                  const int8_t* const lfnst_matrix, int log2_tb_w,
                  int log2_tb_h)
{
  __m128i r[2];
  if (!(log2_tb_w == log2_tb_h)) {
    lfnst_4x4_16(src, r, lfnst_matrix);
  }
  else {
    lfnst_4x4_8(src, r, lfnst_matrix);
  }
  _mm_storel_epi64((__m128i*) &dst[0], r[0]);
  _mm_storel_epi64((__m128i*) &dst[1<<log2_tb_w], _mm_bsrli_si128(r[0], 8));
  _mm_storel_epi64((__m128i*) &dst[2<<log2_tb_w], r[1]);
  _mm_storel_epi64((__m128i*) &dst[3<<log2_tb_w], _mm_bsrli_si128(r[1], 8));
}

static void
compute_lfnst_8x8(const int16_t* const src, int16_t* const dst,
                  const int8_t* const lfnst_matrix, int log2_tb_w,
                  int log2_tb_h)
{
  __m128i r[6];
  if (!(log2_tb_w == log2_tb_h && log2_tb_w == 8)) {
    lfnst_8x8_16(src, &r[0], &lfnst_matrix[0]);
    lfnst_8x8_16(src, &r[2], &lfnst_matrix[16]);
    lfnst_8x8_16(src, &r[4], &lfnst_matrix[32]);
  }
  else {
    lfnst_8x8_8(src, &r[0], &lfnst_matrix[0]);
    lfnst_8x8_8(src, &r[2], &lfnst_matrix[16]);
    lfnst_8x8_8(src, &r[4], &lfnst_matrix[32]);
  }

  _mm_storeu_si128((__m128i *) &dst[0], r[0]);
  _mm_storeu_si128((__m128i *) &dst[(1 << log2_tb_w)], r[1]);
  _mm_storeu_si128((__m128i *) &dst[(2 << log2_tb_w)], r[2]);
  _mm_storeu_si128((__m128i *) &dst[(3 << log2_tb_w)], r[3]);

  _mm_storel_epi64((__m128i*) &dst[4<<log2_tb_w], r[4]);
  _mm_storel_epi64((__m128i*) &dst[5<<log2_tb_w], _mm_bsrli_si128(r[4], 8));
  _mm_storel_epi64((__m128i*) &dst[6<<log2_tb_w], r[5]);
  _mm_storel_epi64((__m128i*) &dst[7<<log2_tb_w], _mm_bsrli_si128(r[5], 8));
}

static void
compute_lfnst_4x4_tr(const int16_t* const src, int16_t* const dst,
                     const int8_t* const lfnst_matrix, int log2_tb_w,
                     int log2_tb_h)
{
        __m128i r[2];
        __m128i t[2];
        if (!(log2_tb_w == log2_tb_h)) {
          lfnst_4x4_16(src, r, lfnst_matrix);
        }
        else {
          lfnst_4x4_8(src, r, lfnst_matrix);
        }

        t[0] = _mm_unpacklo_epi16(r[0], r[1]);
        t[1] = _mm_unpackhi_epi16(r[0], r[1]);

        r[0] = _mm_unpacklo_epi32(t[0], t[1]);
        r[1] = _mm_unpackhi_epi32(t[0], t[1]);

        r[0] = _mm_shufflelo_epi16(r[0], 0xD8);
        r[0] = _mm_shufflehi_epi16(r[0], 0xD8);
        r[1] = _mm_shufflelo_epi16(r[1], 0xD8);
        r[1] = _mm_shufflehi_epi16(r[1], 0xD8);

        _mm_storel_epi64((__m128i*) &dst[0], r[0]);
        _mm_storel_epi64((__m128i*) &dst[1<<log2_tb_w], _mm_bsrli_si128(r[0], 8));
        _mm_storel_epi64((__m128i*) &dst[2<<log2_tb_w], r[1]);
        _mm_storel_epi64((__m128i*) &dst[3<<log2_tb_w], _mm_bsrli_si128(r[1], 8));
}

static void
compute_lfnst_8x8_tr(const int16_t* const src, int16_t* const dst,
                     const int8_t* const lfnst_matrix, int log2_tb_w,
                     int log2_tb_h)
{
        __m128i r[6], t[4];
        if (!(log2_tb_w == log2_tb_h && log2_tb_w == 8)) {
          lfnst_8x8_16(src, &r[0], &lfnst_matrix[0]);
          lfnst_8x8_16(src, &r[2], &lfnst_matrix[16]);
          lfnst_8x8_16(src, &r[4], &lfnst_matrix[32]);
        }
        else {
          lfnst_8x8_8(src, &r[0], &lfnst_matrix[0]);
          lfnst_8x8_8(src, &r[2], &lfnst_matrix[16]);
          lfnst_8x8_8(src, &r[4], &lfnst_matrix[32]);
        }

        t[0] = _mm_unpacklo_epi16(r[0], r[1]);
        t[1] = _mm_unpackhi_epi16(r[0], r[1]);
        t[2] = _mm_unpacklo_epi16(r[2], r[3]);
        t[3] = _mm_unpackhi_epi16(r[2], r[3]);

        r[0] = _mm_unpacklo_epi32(t[0], t[2]);
        r[1] = _mm_unpackhi_epi32(t[0], t[2]);
        r[2] = _mm_unpacklo_epi32(t[1], t[3]);
        r[3] = _mm_unpackhi_epi32(t[1], t[3]);

        t[0] = _mm_unpacklo_epi16(r[4], r[5]);
        t[1] = _mm_unpackhi_epi16(r[4], r[5]);

        r[4] = _mm_unpacklo_epi32(t[0], t[1]);
        r[5] = _mm_unpackhi_epi32(t[0], t[1]);

        r[4] = _mm_shufflelo_epi16(r[4], 0xD8);
        r[4] = _mm_shufflehi_epi16(r[4], 0xD8);
        r[5] = _mm_shufflelo_epi16(r[5], 0xD8);
        r[5] = _mm_shufflehi_epi16(r[5], 0xD8);

        _mm_storel_epi64((__m128i*) &dst[(0<<log2_tb_w)], r[0]);
        _mm_storel_epi64((__m128i*) &dst[(1<<log2_tb_w)], _mm_bsrli_si128(r[0], 8));
        _mm_storel_epi64((__m128i*) &dst[(2<<log2_tb_w)], r[1]);
        _mm_storel_epi64((__m128i*) &dst[(3<<log2_tb_w)], _mm_bsrli_si128(r[1], 8));

        _mm_storel_epi64((__m128i*) &dst[(4<<log2_tb_w)], r[2]);
        _mm_storel_epi64((__m128i*) &dst[(5<<log2_tb_w)], _mm_bsrli_si128(r[2], 8));
        _mm_storel_epi64((__m128i*) &dst[(6<<log2_tb_w)], r[3]);
        _mm_storel_epi64((__m128i*) &dst[(7<<log2_tb_w)], _mm_bsrli_si128(r[3], 8));

        _mm_storel_epi64((__m128i*) &dst[(0<<log2_tb_w) + 4], r[4]);
        _mm_storel_epi64((__m128i*) &dst[(1<<log2_tb_w) + 4], _mm_bsrli_si128(r[4], 8));
        _mm_storel_epi64((__m128i*) &dst[(2<<log2_tb_w) + 4], r[5]);
        _mm_storel_epi64((__m128i*) &dst[(3<<log2_tb_w) + 4], _mm_bsrli_si128(r[5], 8));
}

void rcn_init_lfnst_functions_sse(struct RCNFunctions *rcn_func){
  rcn_func->lfnst.func[0][0] = &compute_lfnst_4x4_sse;
  rcn_func->lfnst.func[0][1] = &compute_lfnst_8x8;
  rcn_func->lfnst.func[1][0] = &compute_lfnst_4x4_tr;
  rcn_func->lfnst.func[1][1] = &compute_lfnst_8x8_tr;
}
