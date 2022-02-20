/**
 *
 *   OpenVVC is open-source real time software decoder compliant with the 
 *   ITU-T H.266- MPEG-I - Part 3 VVC standard. OpenVVC is developed from 
 *   scratch in C as a library that provides consumers with real time and
 *   energy-aware decoding capabilities under different OS including MAC OS,
 *   Windows, Linux and Android targeting low energy real-time decoding of
 *   4K VVC videos on Intel x86 and ARM platforms.
 * 
 *   Copyright (C) 2020-2022  IETR-INSA Rennes :
 *   
 *   Pierre-Loup CABARAT
 *   Wassim HAMIDOUCHE
 *   Guillaume GAUTIER
 *   Thomas AMESTOY
 *
 *   This library is free software; you can redistribute it and/or
 *   modify it under the terms of the GNU Lesser General Public
 *   License as published by the Free Software Foundation; either
 *   version 2.1 of the License, or (at your option) any later version.
 *
 *   This library is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *   Lesser General Public License for more details.
 *
 *   You should have received a copy of the GNU Lesser General Public
 *   License along with this library; if not, write to the Free Software
 *   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301
 *   USA
 * 
 **/

#include <stddef.h>
#include <stdint.h>

#include "data_rcn_transform.h"
#include "rcn_transform.h"
#include "ovutils.h"
#include "rcn_structures.h"

#include "bitdepth.h"

static void
matrix_multiplication(const int16_t* src, const int16_t* const tr_matrix,
                      int16_t* dst, ptrdiff_t src_stride, int tr_size,
                      int num_lines, int num_columns, int shift)
{
    const int round_factor = 1 << (shift - 1);
#if 0
    const int cut_off = tr_size - num_columns;
#endif
    int clip_min = -(1 << 15);
    int clip_max = (1 << 15) - 1;

    for (int i = 0; i < num_lines; i++) {
        for (int j = 0; j < tr_size; j++) {
            int sum = 0;
            for (int k = 0; k < tr_size; k++) {
                sum += (int)src[k * src_stride + i] *
                    tr_matrix[k * tr_size + j];
            }
            dst[i * tr_size + j] =
                ov_clip((int)(sum + round_factor) >> shift,
                        clip_min,
                        clip_max);
        }
    }
}

static void
inverse_dct_ii_2(const int16_t* src, int16_t* dst, ptrdiff_t src_stride,
                     int num_lines, int num_columns, int shift)
{
    int j;
    int E, O;
    int add = 1 << (shift - 1);
    int clip_min = -(1 << 15);
    int clip_max = (1 << 15) - 1;

    // const int  num_lines = line - iSkipLine;
    for (j = 0; j < num_lines; j++) {
        /* Utilizing symmetry properties to the maximum to minimize
           the
         * number of multiplications */
        E = DCT_II_2[0] * (src[0] + src[src_stride]);
        O = DCT_II_2[2] * (src[0] - src[src_stride]);

        /* Combining even and odd terms at each hierarchy levels to
         * calculate the final spatial domain vector */
        dst[0] = ov_clip((E + add) >> shift, clip_min, clip_max);
        dst[1] = ov_clip((O + add) >> shift, clip_min, clip_max);

        src++;
        dst += 2;
    }
}

static void
inverse_dct_ii_4(const int16_t* src, int16_t* dst, ptrdiff_t src_stride,
                     int num_lines, int num_columns, int shift)
{
    int j;
    int E[2], O[2];
    int add = 1 << (shift - 1);
    int clip_min = -(1 << 15);
    int clip_max = (1 << 15) - 1;

    for (j = 0; j < num_lines; j++) {
        /* Utilizing symmetry properties to the maximum to minimize
           the
         * number of multiplications */
        O[0] = DCT_II_4[1 * 4 + 0] * src[src_stride] +
            DCT_II_4[3 * 4 + 0] * src[3 * src_stride];
        O[1] = DCT_II_4[1 * 4 + 1] * src[src_stride] +
            DCT_II_4[3 * 4 + 1] * src[3 * src_stride];
        E[0] = DCT_II_4[0 * 4 + 0] * src[0] +
            DCT_II_4[2 * 4 + 0] * src[2 * src_stride];
        E[1] = DCT_II_4[0 * 4 + 1] * src[0] +
            DCT_II_4[2 * 4 + 1] * src[2 * src_stride];

        /* Combining even and odd terms at each hierarchy levels to
         * calculate the final spatial domain vector */
        dst[0] = ov_clip((E[0] + O[0] + add) >> shift, clip_min, clip_max);
        dst[1] = ov_clip((E[1] + O[1] + add) >> shift, clip_min, clip_max);
        dst[2] = ov_clip((E[1] - O[1] + add) >> shift, clip_min, clip_max);
        dst[3] = ov_clip((E[0] - O[0] + add) >> shift, clip_min, clip_max);

        src++;
        dst += 4;
    }
}

static void
inverse_dct_ii_8(const int16_t* src, int16_t* dst, ptrdiff_t src_stride,
                     int num_lines, int num_columns, int shift)
{
    int j, k;
    int E[4], O[4];
    int EE[2], EO[2];
    int add = 1 << (shift - 1);
    int clip_min = -(1 << 15);
    int clip_max = (1 << 15) - 1;

    for (j = 0; j < num_lines; j++) {
        /* Utilizing symmetry properties to the maximum to minimize
           the
         * number of multiplications */
        for (k = 0; k < 4; k++) {
            O[k] = DCT_II_8[1 * 8 + k] * src[src_stride] +
                DCT_II_8[3 * 8 + k] * src[3 * src_stride] +
                DCT_II_8[5 * 8 + k] * src[5 * src_stride] +
                DCT_II_8[7 * 8 + k] * src[7 * src_stride];
        }

        EO[0] = DCT_II_8[2 * 8 + 0] * src[2 * src_stride] +
            DCT_II_8[6 * 8 + 0] * src[6 * src_stride];
        EO[1] = DCT_II_8[2 * 8 + 1] * src[2 * src_stride] +
            DCT_II_8[6 * 8 + 1] * src[6 * src_stride];
        EE[0] = DCT_II_8[0 * 8 + 0] * src[0] +
            DCT_II_8[4 * 8 + 0] * src[4 * src_stride];
        EE[1] = DCT_II_8[0 * 8 + 1] * src[0] +
            DCT_II_8[4 * 8 + 1] * src[4 * src_stride];

        /* Combining even and odd terms at each hierarchy levels to
         * calculate the final spatial domain vector */
        E[0] = EE[0] + EO[0];
        E[3] = EE[0] - EO[0];
        E[1] = EE[1] + EO[1];
        E[2] = EE[1] - EO[1];

        for (k = 0; k < 4; k++) {
            dst[k] = ov_clip((E[k] + O[k] + add) >> shift, clip_min, clip_max);
            dst[k + 4] = ov_clip((E[3 - k] - O[3 - k] + add) >> shift,
                                 clip_min,
                                 clip_max);
        }
        src++;
        dst += 8;
    }
}

static void
inverse_dct_ii_16(const int16_t* src, int16_t* dst, ptrdiff_t src_stride,
                      int num_lines, int num_columns, int shift)
{
    int j, k;
    int E[8], O[8];
    int EE[4], EO[4];
    int EEE[2], EEO[2];
    int add = 1 << (shift - 1);
    int clip_min = -(1 << 15);
    int clip_max = (1 << 15) - 1;

    for (j = 0; j < num_lines; j++) {
        /* Utilizing symmetry properties to the maximum to minimize
           the
         * number of multiplications */
        for (k = 0; k < 8; k++) {
            O[k] = DCT_II_16[1 * 16 + k] * src[src_stride] +
                DCT_II_16[3 * 16 + k] * src[3 * src_stride] +
                DCT_II_16[5 * 16 + k] * src[5 * src_stride] +
                DCT_II_16[7 * 16 + k] * src[7 * src_stride] +
                DCT_II_16[9 * 16 + k] * src[9 * src_stride] +
                DCT_II_16[11 * 16 + k] * src[11 * src_stride] +
                DCT_II_16[13 * 16 + k] * src[13 * src_stride] +
                DCT_II_16[15 * 16 + k] * src[15 * src_stride];
        }
        for (k = 0; k < 4; k++) {
            EO[k] = DCT_II_16[2 * 16 + k] * src[2 * src_stride] +
                DCT_II_16[6 * 16 + k] * src[6 * src_stride] +
                DCT_II_16[10 * 16 + k] * src[10 * src_stride] +
                DCT_II_16[14 * 16 + k] * src[14 * src_stride];
        }
        EEO[0] = DCT_II_16[4 * 16] * src[4 * src_stride] +
            DCT_II_16[12 * 16] * src[12 * src_stride];
        EEE[0] = DCT_II_16[0] * src[0] +
            DCT_II_16[8 * 16] * src[8 * src_stride];
        EEO[1] = DCT_II_16[4 * 16 + 1] * src[4 * src_stride] +
            DCT_II_16[12 * 16 + 1] * src[12 * src_stride];
        EEE[1] = DCT_II_16[0 * 16 + 1] * src[0] +
            DCT_II_16[8 * 16 + 1] * src[8 * src_stride];

        /* Combining even and odd terms at each hierarchy levels to
         * calculate the final spatial domain vector */
        for (k = 0; k < 2; k++) {
            EE[k] = EEE[k] + EEO[k];
            EE[k + 2] = EEE[1 - k] - EEO[1 - k];
        }
        for (k = 0; k < 4; k++) {
            E[k] = EE[k] + EO[k];
            E[k + 4] = EE[3 - k] - EO[3 - k];
        }
        for (k = 0; k < 8; k++) {
            dst[k] = ov_clip(
                             (E[k] + O[k] + add) >> shift, clip_min, clip_max);
            dst[k + 8] =
                ov_clip((E[7 - k] - O[7 - k] + add) >> shift,
                        clip_min,
                        clip_max);
        }
        src++;
        dst += 16;
    }
}

static void
inverse_dct_ii_32(const int16_t* src, int16_t* dst, ptrdiff_t src_stride,
                      int num_lines, int num_columns, int shift)
{
    int j, k;
    int E[16], O[16];
    int EE[8], EO[8];
    int EEE[4], EEO[4];
    int EEEE[2], EEEO[2];
    int add = 1 << (shift - 1);
    int clip_min = -(1 << 15);
    int clip_max = (1 << 15) - 1;

    for (j = 0; j < num_lines; j++) {
        /* Utilizing symmetry properties to the maximum to minimize
           the
         * number of multiplications */
        for (k = 0; k < 16; k++) {
            O[k] = DCT_II_32[1 * 32 + k] * src[src_stride] +
                DCT_II_32[3 * 32 + k] * src[3 * src_stride] +
                DCT_II_32[5 * 32 + k] * src[5 * src_stride] +
                DCT_II_32[7 * 32 + k] * src[7 * src_stride] +
                DCT_II_32[9 * 32 + k] * src[9 * src_stride] +
                DCT_II_32[11 * 32 + k] * src[11 * src_stride] +
                DCT_II_32[13 * 32 + k] * src[13 * src_stride] +
                DCT_II_32[15 * 32 + k] * src[15 * src_stride] +
                DCT_II_32[17 * 32 + k] * src[17 * src_stride] +
                DCT_II_32[19 * 32 + k] * src[19 * src_stride] +
                DCT_II_32[21 * 32 + k] * src[21 * src_stride] +
                DCT_II_32[23 * 32 + k] * src[23 * src_stride] +
                DCT_II_32[25 * 32 + k] * src[25 * src_stride] +
                DCT_II_32[27 * 32 + k] * src[27 * src_stride] +
                DCT_II_32[29 * 32 + k] * src[29 * src_stride] +
                DCT_II_32[31 * 32 + k] * src[31 * src_stride];
        }
        for (k = 0; k < 8; k++) {
            EO[k] = DCT_II_32[2 * 32 + k] * src[2 * src_stride] +
                DCT_II_32[6 * 32 + k] * src[6 * src_stride] +
                DCT_II_32[10 * 32 + k] * src[10 * src_stride] +
                DCT_II_32[14 * 32 + k] * src[14 * src_stride] +
                DCT_II_32[18 * 32 + k] * src[18 * src_stride] +
                DCT_II_32[22 * 32 + k] * src[22 * src_stride] +
                DCT_II_32[26 * 32 + k] * src[26 * src_stride] +
                DCT_II_32[30 * 32 + k] * src[30 * src_stride];
        }
        for (k = 0; k < 4; k++) {
            EEO[k] = DCT_II_32[4 * 32 + k] * src[4 * src_stride] +
                DCT_II_32[12 * 32 + k] * src[12 * src_stride] +
                DCT_II_32[20 * 32 + k] * src[20 * src_stride] +
                DCT_II_32[28 * 32 + k] * src[28 * src_stride];
        }
        EEEO[0] = DCT_II_32[8 * 32 + 0] * src[8 * src_stride] +
            DCT_II_32[24 * 32 + 0] * src[24 * src_stride];
        EEEO[1] = DCT_II_32[8 * 32 + 1] * src[8 * src_stride] +
            DCT_II_32[24 * 32 + 1] * src[24 * src_stride];
        EEEE[0] = DCT_II_32[0 * 32 + 0] * src[0] +
            DCT_II_32[16 * 32 + 0] * src[16 * src_stride];
        EEEE[1] = DCT_II_32[0 * 32 + 1] * src[0] +
            DCT_II_32[16 * 32 + 1] * src[16 * src_stride];

        /* Combining even and odd terms at each hierarchy levels to
         * calculate the final spatial domain vector */
        EEE[0] = EEEE[0] + EEEO[0];
        EEE[3] = EEEE[0] - EEEO[0];
        EEE[1] = EEEE[1] + EEEO[1];
        EEE[2] = EEEE[1] - EEEO[1];
        for (k = 0; k < 4; k++) {
            EE[k] = EEE[k] + EEO[k];
            EE[k + 4] = EEE[3 - k] - EEO[3 - k];
        }
        for (k = 0; k < 8; k++) {
            E[k] = EE[k] + EO[k];
            E[k + 8] = EE[7 - k] - EO[7 - k];
        }
        for (k = 0; k < 16; k++) {
            dst[k] = ov_clip(
                             (E[k] + O[k] + add) >> shift, clip_min, clip_max);
            dst[k + 16] =
                ov_clip((E[15 - k] - O[15 - k] + add) >> shift,
                        clip_min,
                        clip_max);
        }
        src++;
        dst += 32;
    }
}

static void
inverse_dct_ii_64(const int16_t* src, int16_t* dst, ptrdiff_t src_stride,
                      int num_lines, int num_columns, int shift)
{
    int round_factor = 1 << (shift - 1);

    int clip_min = -(1 << 15);
    int clip_max = (1 << 15) - 1;
    int j, k;
    int E[32], O[32];
    int EE[16], EO[16];
    int EEE[8], EEO[8];
    int EEEE[4], EEEO[4];
    int EEEEE[2], EEEEO[2];
    for (j = 0; j < num_lines; j++) {
        /* Utilizing symmetry properties to the maximum to minimize
           the
         * number of multiplications */
        for (k = 0; k < 32; k++) {
            O[k] =
                DCT_II_64_OT[k * 16 + 0 ] * src[src_stride] +
                DCT_II_64_OT[k * 16 + 1 ] * src[3 * src_stride] +
                DCT_II_64_OT[k * 16 + 2 ] * src[5 * src_stride] +
                DCT_II_64_OT[k * 16 + 3 ] * src[7 * src_stride] +
                DCT_II_64_OT[k * 16 + 4 ] * src[9 * src_stride] +
                DCT_II_64_OT[k * 16 + 5 ] * src[11 * src_stride] +
                DCT_II_64_OT[k * 16 + 6 ] * src[13 * src_stride] +
                DCT_II_64_OT[k * 16 + 7 ] * src[15 * src_stride] +
                DCT_II_64_OT[k * 16 + 8 ] * src[17 * src_stride] +
                DCT_II_64_OT[k * 16 + 9 ] * src[19 * src_stride] +
                DCT_II_64_OT[k * 16 + 10] * src[21 * src_stride] +
                DCT_II_64_OT[k * 16 + 11] * src[23 * src_stride] +
                DCT_II_64_OT[k * 16 + 12] * src[25 * src_stride] +
                DCT_II_64_OT[k * 16 + 13] * src[27 * src_stride] +
                DCT_II_64_OT[k * 16 + 14] * src[29 * src_stride] +
                DCT_II_64_OT[k * 16 + 15] * src[31 * src_stride];
        }
        for (k = 0; k < 16; k++) {
            EO[k] =
                DCT_II_64_EOT[k * 8 + 0] * src[2 * src_stride] +
                DCT_II_64_EOT[k * 8 + 1] * src[6 * src_stride] +
                DCT_II_64_EOT[k * 8 + 2] * src[10 * src_stride] +
                DCT_II_64_EOT[k * 8 + 3] * src[14 * src_stride] +
                DCT_II_64_EOT[k * 8 + 4] * src[18 * src_stride] +
                DCT_II_64_EOT[k * 8 + 5] * src[22 * src_stride] +
                DCT_II_64_EOT[k * 8 + 6] * src[26 * src_stride] +
                DCT_II_64_EOT[k * 8 + 7] * src[30 * src_stride];
        }
        for (k = 0; k < 8; k++) {
            EEO[k] =
                DCT_II_64_EEOT[k * 4 + 0] * src[4 * src_stride] +
                DCT_II_64_EEOT[k * 4 + 1] * src[12 * src_stride] +
                DCT_II_64_EEOT[k * 4 + 2] * src[20 * src_stride] +
                DCT_II_64_EEOT[k * 4 + 3] * src[28 * src_stride];
        }
        for (k = 0; k < 4; k++) {
            EEEO[k] =
                DCT_II_64_EEEOT[k * 2 + 0] * src[8 * src_stride] +
                DCT_II_64_EEEOT[k * 2 + 1] * src[24 * src_stride];
        }
        EEEEO[0] =
            DCT_II_64_EEEEO[0] * src[16 * src_stride];
        EEEEO[1] =
            DCT_II_64_EEEEO[1] * src[16 * src_stride];
        EEEEE[0] =
            DCT_II_64_EEEEE[0] * src[0];
        EEEEE[1] =
            DCT_II_64_EEEEE[1] * src[0];

        /* Combining even and odd terms at each hierarchy levels to
         * calculate the final spatial domain vector */
        for (k = 0; k < 2; k++) {
            EEEE[k] = EEEEE[k] + EEEEO[k];
            EEEE[k + 2] = EEEEE[1 - k] - EEEEO[1 - k];
        }
        for (k = 0; k < 4; k++) {
            EEE[k] = EEEE[k] + EEEO[k];
            EEE[k + 4] = EEEE[3 - k] - EEEO[3 - k];
        }
        for (k = 0; k < 8; k++) {
            EE[k] = EEE[k] + EEO[k];
            EE[k + 8] = EEE[7 - k] - EEO[7 - k];
        }
        for (k = 0; k < 16; k++) {
            E[k] = EE[k] + EO[k];
            E[k + 16] = EE[15 - k] - EO[15 - k];
        }
        for (k = 0; k < 32; k++) {
            dst[k] = ov_clip((E[k] + O[k] + round_factor) >> shift,
                             clip_min,
                             clip_max);
            dst[k + 32] = ov_clip(
                                  (E[31 - k] - O[31 - k] + round_factor) >> shift,
                                  clip_min,
                                  clip_max);
        }
        src++;
        dst += 64;
    }
}

static void
inverse_dct_viii_4(const int16_t* src, int16_t* dst, ptrdiff_t src_stride,
                       int num_lines, int num_columns, int shift)
{
    int i;
    int clip_min = -(1 << 15);
    int clip_max = (1 << 15) - 1;
    int round_factor = 1 << (shift - 1);

    int c[4];
    //    const int  reducedLine = line - iSkipLine;
    for (i = 0; i < num_lines; i++) {
        // Intermediate Variables
        c[0] = src[0 * src_stride] + src[3 * src_stride];
        c[1] = src[2 * src_stride] + src[0 * src_stride];
        c[2] = src[3 * src_stride] - src[2 * src_stride];
        c[3] = (int)DCT_VIII_4[1] * src[1 * src_stride];

        dst[0] =
            ov_clip(((int)DCT_VIII_4[3] * c[0] +
                     (int)DCT_VIII_4[2] * c[1] + c[3] + round_factor) >>
                    shift,
                    clip_min,
                    clip_max);
        dst[1] = ov_clip(((int)DCT_VIII_4[1] *
                          (src[0 * src_stride] - src[2 * src_stride] -
                           src[3 * src_stride]) +
                          round_factor) >>
                         shift,
                         clip_min,
                         clip_max);
        dst[2] =
            ov_clip(((int)DCT_VIII_4[3] * c[2] +
                     (int)DCT_VIII_4[2] * c[0] - c[3] + round_factor) >>
                    shift,
                    clip_min,
                    clip_max);
        dst[3] =
            ov_clip(((int)DCT_VIII_4[3] * c[1] -
                     (int)DCT_VIII_4[2] * c[2] - c[3] + round_factor) >>
                    shift,
                    clip_min,
                    clip_max);

        dst += 4;
        src++;
    }
}

static void
inverse_dct_viii_8(const int16_t* src, int16_t* dst, ptrdiff_t src_stride,
                       int num_lines, int num_columns, int shift)
{
    matrix_multiplication(
                          src, DCT_VIII_8, dst, src_stride, 8, num_lines, num_columns, shift);
}

static void
inverse_dct_viii_16(const int16_t* src, int16_t* dst, ptrdiff_t src_stride,
                        int num_lines, int num_columns, int shift)
{
    matrix_multiplication(
                          src, DCT_VIII_16, dst, src_stride, 16, num_lines, num_columns, shift);
}

static void
inverse_dct_viii_32(const int16_t* src, int16_t* dst, ptrdiff_t src_stride,
                        int num_lines, int num_columns, int shift)
{
    matrix_multiplication(
                          src, DCT_VIII_32, dst, src_stride, 32, num_lines, num_columns, shift);
}

static void
inverse_dst_vii_4(const int16_t* src, int16_t* dst, ptrdiff_t src_stride,
                      int num_lines, int num_columns, int shift)
{
    int i;
    int clip_min = -(1 << 15);
    int clip_max = (1 << 15) - 1;
    int round_factor = 1 << (shift - 1);

    int c[4];
    //    const int  reducedLine = line - iSkipLine;
    for (i = 0; i < num_lines; i++) {
        // Intermediate Variables
        c[0] = src[0 * src_stride] + src[2 * src_stride];
        c[1] = src[2 * src_stride] + src[3 * src_stride];
        c[2] = src[0 * src_stride] - src[3 * src_stride];
        c[3] = DST_VII_4[2] * src[1 * src_stride];

        dst[0] = ov_clip((DST_VII_4[0] * c[0] + DST_VII_4[1] * c[1] +
                          c[3] + round_factor) >>
                         shift,
                         clip_min,
                         clip_max);
        dst[1] = ov_clip((DST_VII_4[1] * c[2] - DST_VII_4[0] * c[1] +
                          c[3] + round_factor) >>
                         shift,
                         clip_min,
                         clip_max);
        dst[2] = ov_clip(
                         (DST_VII_4[2] * (src[0 * src_stride] - src[2 * src_stride] +
                                          src[3 * src_stride]) +
                          round_factor) >>
                         shift,
                         clip_min,
                         clip_max);
        dst[3] = ov_clip((DST_VII_4[1] * c[0] + DST_VII_4[0] * c[2] -
                          c[3] + round_factor) >>
                         shift,
                         clip_min,
                         clip_max);

        dst += 4;
        src++;
    }
}

static void
inverse_dst_vii_8(const int16_t* src, int16_t* dst, ptrdiff_t src_stride,
                      int num_lines, int num_columns, int shift)
{

    matrix_multiplication(
                          src, DST_VII_8, dst, src_stride, 8, num_lines, num_columns, shift);
}

static void
inverse_dst_vii_16(const int16_t* src, int16_t* dst, ptrdiff_t src_stride,
                       int num_lines, int num_columns, int shift)
{
    matrix_multiplication(
                          src, DST_VII_16, dst, src_stride, 16, num_lines, num_columns, shift);
}

static void
inverse_dst_vii_32(const int16_t* src, int16_t* dst, ptrdiff_t src_stride,
                       int num_lines, int num_columns, int shift)
{
    matrix_multiplication(
                          src, DST_VII_32, dst, src_stride, 32, num_lines, num_columns, shift);
}

static void
inverse_dct_ii_dc(int16_t* const dst, int log2_tb_w, int log2_tb_h,
                      int dc_val)
{
    int i, j;
    int16_t* _dst = dst;
    const int tb_w = 1 << log2_tb_w;
    const int tb_h = 1 << log2_tb_h;
    int clip_min = -(1 << 15);
    int clip_max = (1 << 15) - 1;
    const int add = 1 << ((14 - BITDEPTH) - 1);
    int value = (((dc_val + 1) >> 1) + add) >> (14 - BITDEPTH);
    value = ov_clip(value, clip_min, clip_max);
    for (i = 0; i < tb_h; ++i) {
        for (j = 0; j < tb_w; ++j) {
            _dst[j] = value;
        }
        _dst += tb_w;
    }
}

void
BD_DECL(rcn_init_tr_functions)(struct RCNFunctions *const rcn_funcs)
{
    rcn_funcs->tr.func[DST_VII][0] = NULL;
    rcn_funcs->tr.func[DST_VII][1] = NULL;
    rcn_funcs->tr.func[DST_VII][2] = &inverse_dst_vii_4;
    rcn_funcs->tr.func[DST_VII][3] = &inverse_dst_vii_8;
    rcn_funcs->tr.func[DST_VII][4] = &inverse_dst_vii_16;
    rcn_funcs->tr.func[DST_VII][5] = &inverse_dst_vii_32;
    rcn_funcs->tr.func[DST_VII][6] = NULL;

    rcn_funcs->tr.func[DCT_VIII][0] = NULL;
    rcn_funcs->tr.func[DCT_VIII][1] = NULL;
    rcn_funcs->tr.func[DCT_VIII][2] = &inverse_dct_viii_4;
    rcn_funcs->tr.func[DCT_VIII][3] = &inverse_dct_viii_8;
    rcn_funcs->tr.func[DCT_VIII][4] = &inverse_dct_viii_16;
    rcn_funcs->tr.func[DCT_VIII][5] = &inverse_dct_viii_32;
    rcn_funcs->tr.func[DCT_VIII][6] = NULL;

    rcn_funcs->tr.func[DCT_II][0] = NULL;
    rcn_funcs->tr.func[DCT_II][1] = &inverse_dct_ii_2;
    rcn_funcs->tr.func[DCT_II][2] = &inverse_dct_ii_4;
    rcn_funcs->tr.func[DCT_II][3] = &inverse_dct_ii_8;
    rcn_funcs->tr.func[DCT_II][4] = &inverse_dct_ii_16;
    rcn_funcs->tr.func[DCT_II][5] = &inverse_dct_ii_32;
    rcn_funcs->tr.func[DCT_II][6] = &inverse_dct_ii_64;

    rcn_funcs->tr.dc = &inverse_dct_ii_dc;
}
