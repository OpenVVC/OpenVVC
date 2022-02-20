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

#include "rcn_neon.h"

#include "simde/x86/sse4.2.h"
#include "simde/x86/mmx.h"
#include "simde/x86/sse.h"
#include "simde/x86/sse2.h"
#include "simde/x86/sse3.h"
#include "simde/x86/ssse3.h"
#include "simde/x86/sse4.2.h"
#include "simde/x86/avx.h"
#include "simde/x86/avx2.h"
#include "simde/x86/avx512.h"
#include <stdio.h>
#include <stddef.h>
#include <stdint.h>

#include "rcn_transform.h"

typedef int16_t TCoeff;
typedef int16_t TMatrixCoeff;

#define SHIFTV 7
#define SHIFTH 10
#define SIZEV 32
#define SIZEH 2

void afficherVecteur4SSE128(__m128i cible)
{
    int32_t test[4];
    _mm_store_si128((__m128i *)test, cible);

//    printf(" ________ ________ ________ ________\n|        |        |        |        |\n");
    printf("| %06d | %06d | %06d | %06d |\n",  test[0], test[1], test[2], test[3]);
//    printf("|________|________|________|________|\n");
}

void afficherVecteur8SSE128(__m128i cible)
{
    int16_t test[8];
    _mm_store_si128((__m128i *)test, cible);

//    printf(" ________ ________ ________ ________ ________ ________ ________ ________\n|        |        |        |        |        |        |        |        |\n");
    printf("| %06d | %06d | %06d | %06d | %06d | %06d | %06d | %06d |\n", test[0], test[1], test[2], test[3], test[4], test[5], test[6], test[7]);
//    printf("|________|________|________|________|________|________|________|________|\n");
}

void afficherTableau2D(int v, int h, TCoeff * m){
    int i,j;
    for (i = 0; i < v; ++i) {
        for (j = 0; j < h; ++j) {
            printf("%i ",m[i*h+j]);
        }
        printf("\n");
    }
    printf("\n");
}
static void
transpose8x8(__m128i *in, __m128i *out){
  __m128i tmp1[8], tmp2[8];
  tmp1[0] = _mm_unpacklo_epi16(in[0], in[1]);
  tmp1[1] = _mm_unpackhi_epi16(in[0], in[1]);
  tmp1[2] = _mm_unpacklo_epi16(in[2], in[3]);
  tmp1[3] = _mm_unpackhi_epi16(in[2], in[3]);
  tmp1[4] = _mm_unpacklo_epi16(in[4], in[5]);
  tmp1[5] = _mm_unpackhi_epi16(in[4], in[5]);
  tmp1[6] = _mm_unpacklo_epi16(in[6], in[7]);
  tmp1[7] = _mm_unpackhi_epi16(in[6], in[7]);

  tmp2[0] = _mm_unpacklo_epi32(tmp1[0], tmp1[2]);
  tmp2[1] = _mm_unpackhi_epi32(tmp1[0], tmp1[2]);
  tmp2[2] = _mm_unpacklo_epi32(tmp1[1], tmp1[3]);
  tmp2[3] = _mm_unpackhi_epi32(tmp1[1], tmp1[3]);
  tmp2[4] = _mm_unpacklo_epi32(tmp1[4], tmp1[6]);
  tmp2[5] = _mm_unpackhi_epi32(tmp1[4], tmp1[6]);
  tmp2[6] = _mm_unpacklo_epi32(tmp1[5], tmp1[7]);
  tmp2[7] = _mm_unpackhi_epi32(tmp1[5], tmp1[7]);

  out[0] = _mm_unpacklo_epi64(tmp2[0], tmp2[4]);
  out[1] = _mm_unpackhi_epi64(tmp2[0], tmp2[4]);
  out[2] = _mm_unpacklo_epi64(tmp2[1], tmp2[5]);
  out[3] = _mm_unpackhi_epi64(tmp2[1], tmp2[5]);
  out[4] = _mm_unpacklo_epi64(tmp2[2], tmp2[6]);
  out[5] = _mm_unpackhi_epi64(tmp2[2], tmp2[6]);
  out[6] = _mm_unpacklo_epi64(tmp2[3], tmp2[7]);
  out[7] = _mm_unpackhi_epi64(tmp2[3], tmp2[7]);
}

void inverse_sse2_B4(const TCoeff *src, TCoeff *dst, int src_stride, int shift, int line, const TMatrixCoeff* iT){
    __m128i x1, x2; //Contient le vecteur à transformer
    __m128i d, d2, d3, d4;	//Contient coefficient DCT ou DST
    __m128i result[line]; //Variables pour calculs (result[] il faudrait mettre line ou la taille max 32)
    __m128i add = _mm_set1_epi32(1 << (shift - 1));

    int nbstore = line/2;

    d = _mm_load_si128((__m128i *)&(iT[0]));
    d2 = _mm_load_si128((__m128i *)&(iT[8])); //8 car [taille coef DCT ou DST 4] x [2 (le milieu de la hauteur)]
    d3 = _mm_unpacklo_epi16(d,d2);
    d4 = _mm_unpackhi_epi16(d,d2);
    for(int i = 0; i < line; ++i){
        x1 = _mm_unpacklo_epi16(_mm_set1_epi16(src[i]),_mm_set1_epi16(src[2*src_stride+i]));
        x2 = _mm_unpacklo_epi16(_mm_set1_epi16(src[src_stride+i]),_mm_set1_epi16(src[3*src_stride+i]));

        result[i] = _mm_add_epi32(_mm_madd_epi16(x1, d3),_mm_madd_epi16(x2, d4));

        result[i] = _mm_add_epi32(result[i], add);

        result[i] = _mm_srai_epi32(result[i], shift);
    }

    for(int i = 0; i<nbstore; i++){
        result[i] = _mm_packs_epi32(result[2*i], result[2*i+1]); //clip pour repasser en 16
        _mm_store_si128((__m128i *)&(dst[i*8]), result[i]);	//dst[i*8] car result contient 8 résultat
    }
}


void inverse_sse2_B8(const TCoeff *src, TCoeff *dst, int src_stride, int shift, int line, const TMatrixCoeff* iT){
  __m128i add = _mm_set1_epi32(1 << (shift - 1));
  TCoeff *_src = (TCoeff *)src;
  TCoeff *_dst = dst;
  uint8_t i, k;
  __m128i x[8], d[8], m[8], di[8], a[8], b[4], r[2], res[8];
  for (i = 0; i < line>>3; i++) {
    x[0]=_mm_load_si128((__m128i*)(_src + 0 * src_stride));
    x[1]=_mm_load_si128((__m128i*)(_src + 1 * src_stride));
    x[2]=_mm_load_si128((__m128i*)(_src + 2 * src_stride));
    x[3]=_mm_load_si128((__m128i*)(_src + 3 * src_stride));
    x[4]=_mm_load_si128((__m128i*)(_src + 4 * src_stride));
    x[5]=_mm_load_si128((__m128i*)(_src + 5 * src_stride));
    x[6]=_mm_load_si128((__m128i*)(_src + 6 * src_stride));
    x[7]=_mm_load_si128((__m128i*)(_src + 7 * src_stride));

    for (k = 0; k < 8; k++) {
      d[0] = _mm_set1_epi16(iT[k + 0 * 8]);
      d[1] = _mm_set1_epi16(iT[k + 1 * 8]);
      d[2] = _mm_set1_epi16(iT[k + 2 * 8]);
      d[3] = _mm_set1_epi16(iT[k + 3 * 8]);
      d[4] = _mm_set1_epi16(iT[k + 4 * 8]);
      d[5] = _mm_set1_epi16(iT[k + 5 * 8]);
      d[6] = _mm_set1_epi16(iT[k + 6 * 8]);
      d[7] = _mm_set1_epi16(iT[k + 7 * 8]);

      m[0] = _mm_unpacklo_epi16(x[0],  x[1]);
      m[1] = _mm_unpacklo_epi16(x[2],  x[3]);
      m[2] = _mm_unpacklo_epi16(x[4],  x[5]);
      m[3] = _mm_unpacklo_epi16(x[6],  x[7]);

      m[4] = _mm_unpackhi_epi16(x[0],  x[1]);
      m[5] = _mm_unpackhi_epi16(x[2],  x[3]);
      m[6] = _mm_unpackhi_epi16(x[4],  x[5]);
      m[7] = _mm_unpackhi_epi16(x[6],  x[7]);

      di[0] = _mm_unpacklo_epi16(d[0],  d[1]);
      di[1] = _mm_unpacklo_epi16(d[2],  d[3]);
      di[2] = _mm_unpacklo_epi16(d[4],  d[5]);
      di[3] = _mm_unpacklo_epi16(d[6],  d[7]);

      di[4] = _mm_unpackhi_epi16(d[0],  d[1]);
      di[5] = _mm_unpackhi_epi16(d[2],  d[3]);
      di[6] = _mm_unpackhi_epi16(d[4],  d[5]);
      di[7] = _mm_unpackhi_epi16(d[6],  d[7]);


      a[0] = _mm_madd_epi16(m[0], di[0]);
      a[1] = _mm_madd_epi16(m[1], di[1]);
      a[2] = _mm_madd_epi16(m[2], di[2]);
      a[3] = _mm_madd_epi16(m[3], di[3]);

      a[4] = _mm_madd_epi16(m[4], di[4]);
      a[5] = _mm_madd_epi16(m[5], di[5]);
      a[6] = _mm_madd_epi16(m[6], di[6]);
      a[7] = _mm_madd_epi16(m[7], di[7]);

      b[0] = _mm_add_epi32(a[0], a[2]);
      b[1] = _mm_add_epi32(a[1], a[3]);

      b[2] = _mm_add_epi32(a[4], a[6]);
      b[3] = _mm_add_epi32(a[5], a[7]);

      r[0] = _mm_add_epi32(b[0], b[1]);
      r[1] = _mm_add_epi32(b[2], b[3]);

      r[0] = _mm_add_epi32(r[0], add);
      r[1] = _mm_add_epi32(r[1], add);

      r[0] = _mm_srai_epi32(r[0], shift);
      r[1] = _mm_srai_epi32(r[1], shift);

      res[k] = _mm_packs_epi32(r[0], r[1]);
    }
        transpose8x8(res,m);

        _mm_store_si128((__m128i *)&_dst[0],  m[0]);
        _mm_store_si128((__m128i *)&_dst[8],  m[1]);
        _mm_store_si128((__m128i *)&_dst[16], m[2]);
        _mm_store_si128((__m128i *)&_dst[24], m[3]);
        _mm_store_si128((__m128i *)&_dst[32], m[4]);
        _mm_store_si128((__m128i *)&_dst[40], m[5]);
        _mm_store_si128((__m128i *)&_dst[48], m[6]);
        _mm_store_si128((__m128i *)&_dst[56], m[7]);
        _src+=8;
        _dst+=64;
    }

    if (!(line&0x7)) {
      return ;
    }

    d[0] = _mm_load_si128((__m128i *)&(iT[0*8]));
    d[1] = _mm_load_si128((__m128i *)&(iT[1*8]));
    d[2] = _mm_load_si128((__m128i *)&(iT[2*8]));
    d[3] = _mm_load_si128((__m128i *)&(iT[3*8]));
    d[4] = _mm_load_si128((__m128i *)&(iT[4*8]));
    d[5] = _mm_load_si128((__m128i *)&(iT[5*8]));
    d[6] = _mm_load_si128((__m128i *)&(iT[6*8]));
    d[7] = _mm_load_si128((__m128i *)&(iT[7*8]));

    for(i = 0; i < (line&0x7); ++i){
      x[0] = _mm_set1_epi16(_src[0*src_stride+i]);
      x[1] = _mm_set1_epi16(_src[1*src_stride+i]);
      x[2] = _mm_set1_epi16(_src[2*src_stride+i]);
      x[3] = _mm_set1_epi16(_src[3*src_stride+i]);
      x[4] = _mm_set1_epi16(_src[4*src_stride+i]);
      x[5] = _mm_set1_epi16(_src[5*src_stride+i]);
      x[6] = _mm_set1_epi16(_src[6*src_stride+i]);
      x[7] = _mm_set1_epi16(_src[7*src_stride+i]);

      m[0] = _mm_unpacklo_epi16(x[0],  x[1]);
      m[1] = _mm_unpacklo_epi16(x[2],  x[3]);
      m[2] = _mm_unpacklo_epi16(x[4],  x[5]);
      m[3] = _mm_unpacklo_epi16(x[6],  x[7]);

      m[4] = _mm_unpackhi_epi16(x[0],  x[1]);
      m[5] = _mm_unpackhi_epi16(x[2],  x[3]);
      m[6] = _mm_unpackhi_epi16(x[4],  x[5]);
      m[7] = _mm_unpackhi_epi16(x[6],  x[7]);

      di[0] = _mm_unpacklo_epi16(d[0],  d[1]);
      di[1] = _mm_unpacklo_epi16(d[2],  d[3]);
      di[2] = _mm_unpacklo_epi16(d[4],  d[5]);
      di[3] = _mm_unpacklo_epi16(d[6],  d[7]);

      di[4] = _mm_unpackhi_epi16(d[0],  d[1]);
      di[5] = _mm_unpackhi_epi16(d[2],  d[3]);
      di[6] = _mm_unpackhi_epi16(d[4],  d[5]);
      di[7] = _mm_unpackhi_epi16(d[6],  d[7]);

      a[0] = _mm_madd_epi16(m[0], di[0]);
      a[1] = _mm_madd_epi16(m[1], di[1]);
      a[2] = _mm_madd_epi16(m[2], di[2]);
      a[3] = _mm_madd_epi16(m[3], di[3]);

      a[4] = _mm_madd_epi16(m[4], di[4]);
      a[5] = _mm_madd_epi16(m[5], di[5]);
      a[6] = _mm_madd_epi16(m[6], di[6]);
      a[7] = _mm_madd_epi16(m[7], di[7]);

      b[0] = _mm_add_epi32(a[0], a[2]);
      b[1] = _mm_add_epi32(a[1], a[3]);

      b[2] = _mm_add_epi32(a[4], a[6]);
      b[3] = _mm_add_epi32(a[5], a[7]);

      r[0] = _mm_add_epi32(b[0], b[1]);
      r[1] = _mm_add_epi32(b[2], b[3]);

      r[0] = _mm_add_epi32(r[0], add);
      r[1] = _mm_add_epi32(r[1], add);

      r[0] = _mm_srai_epi32(r[0], shift);
      r[1] = _mm_srai_epi32(r[1], shift);

      _mm_store_si128((__m128i *)&(_dst[i*8]), _mm_packs_epi32(r[0], r[1])); //dst[i*8] car result contient 8 résultat
    }
}

void inverse_sse2_B16(const TCoeff *src, TCoeff *dst, int src_stride, int shift, int line, const TMatrixCoeff* iT){
    int i, k;
    TCoeff *_src = (TCoeff *)src;
    TCoeff *_dst = dst;
    __m128i x[16], m[16];
    __m128i d[16], di[8], a[16], b[8], c[4], r[2], res[16];
    __m128i add = _mm_set1_epi32(1 << (shift - 1));

    for (i = 0; i < line>>3; i++) {
      x[0 ]=_mm_load_si128((__m128i*)(_src +  0  * src_stride));
      x[1 ]=_mm_load_si128((__m128i*)(_src +  1  * src_stride));
      x[2 ]=_mm_load_si128((__m128i*)(_src +  2  * src_stride));
      x[3 ]=_mm_load_si128((__m128i*)(_src +  3  * src_stride));
      x[4 ]=_mm_load_si128((__m128i*)(_src +  4  * src_stride));
      x[5 ]=_mm_load_si128((__m128i*)(_src +  5  * src_stride));
      x[6 ]=_mm_load_si128((__m128i*)(_src +  6  * src_stride));
      x[7 ]=_mm_load_si128((__m128i*)(_src +  7  * src_stride));
      x[8 ]=_mm_load_si128((__m128i*)(_src +  8  * src_stride));
      x[9 ]=_mm_load_si128((__m128i*)(_src +  9  * src_stride));
      x[10]=_mm_load_si128((__m128i*)(_src +  10 * src_stride));
      x[11]=_mm_load_si128((__m128i*)(_src +  11 * src_stride));
      x[12]=_mm_load_si128((__m128i*)(_src +  12 * src_stride));
      x[13]=_mm_load_si128((__m128i*)(_src +  13 * src_stride));
      x[14]=_mm_load_si128((__m128i*)(_src +  14 * src_stride));
      x[15]=_mm_load_si128((__m128i*)(_src +  15 * src_stride));
      for (k = 0; k < 16; k++) {
        d[0 ] = _mm_set1_epi16(iT[k + 0  * 16]);
        d[1 ] = _mm_set1_epi16(iT[k + 1  * 16]);
        d[2 ] = _mm_set1_epi16(iT[k + 2  * 16]);
        d[3 ] = _mm_set1_epi16(iT[k + 3  * 16]);
        d[4 ] = _mm_set1_epi16(iT[k + 4  * 16]);
        d[5 ] = _mm_set1_epi16(iT[k + 5  * 16]);
        d[6 ] = _mm_set1_epi16(iT[k + 6  * 16]);
        d[7 ] = _mm_set1_epi16(iT[k + 7  * 16]);
        d[8 ] = _mm_set1_epi16(iT[k + 8  * 16]);
        d[9 ] = _mm_set1_epi16(iT[k + 9  * 16]);
        d[10] = _mm_set1_epi16(iT[k + 10 * 16]);
        d[11] = _mm_set1_epi16(iT[k + 11 * 16]);
        d[12] = _mm_set1_epi16(iT[k + 12 * 16]);
        d[13] = _mm_set1_epi16(iT[k + 13 * 16]);
        d[14] = _mm_set1_epi16(iT[k + 14 * 16]);
        d[15] = _mm_set1_epi16(iT[k + 15 * 16]);

        m[0 ] = _mm_unpacklo_epi16(x[0 ],  x[1 ]);
        m[1 ] = _mm_unpacklo_epi16(x[2 ],  x[3 ]);
        m[2 ] = _mm_unpacklo_epi16(x[4 ],  x[5 ]);
        m[3 ] = _mm_unpacklo_epi16(x[6 ],  x[7 ]);
        m[4 ] = _mm_unpacklo_epi16(x[8 ],  x[9 ]);
        m[5 ] = _mm_unpacklo_epi16(x[10],  x[11]);
        m[6 ] = _mm_unpacklo_epi16(x[12],  x[13]);
        m[7 ] = _mm_unpacklo_epi16(x[14],  x[15]);

        m[8 ] = _mm_unpackhi_epi16(x[0 ],  x[1 ]);
        m[9 ] = _mm_unpackhi_epi16(x[2 ],  x[3 ]);
        m[10] = _mm_unpackhi_epi16(x[4 ],  x[5 ]);
        m[11] = _mm_unpackhi_epi16(x[6 ],  x[7 ]);
        m[12] = _mm_unpackhi_epi16(x[8 ],  x[9 ]);
        m[13] = _mm_unpackhi_epi16(x[10],  x[11]);
        m[14] = _mm_unpackhi_epi16(x[12],  x[13]);
        m[15] = _mm_unpackhi_epi16(x[14],  x[15]);

        di[0 ] = _mm_unpacklo_epi16(d[0 ],  d[1 ]);
        di[1 ] = _mm_unpacklo_epi16(d[2 ],  d[3 ]);
        di[2 ] = _mm_unpacklo_epi16(d[4 ],  d[5 ]);
        di[3 ] = _mm_unpacklo_epi16(d[6 ],  d[7 ]);
        di[4 ] = _mm_unpacklo_epi16(d[8 ],  d[9 ]);
        di[5 ] = _mm_unpacklo_epi16(d[10],  d[11]);
        di[6 ] = _mm_unpacklo_epi16(d[12],  d[13]);
        di[7 ] = _mm_unpacklo_epi16(d[14],  d[15]);

        a[0 ] = _mm_madd_epi16(m[0 ], di[0 ]);
        a[1 ] = _mm_madd_epi16(m[1 ], di[1 ]);
        a[2 ] = _mm_madd_epi16(m[2 ], di[2 ]);
        a[3 ] = _mm_madd_epi16(m[3 ], di[3 ]);
        a[4 ] = _mm_madd_epi16(m[4 ], di[4 ]);
        a[5 ] = _mm_madd_epi16(m[5 ], di[5 ]);
        a[6 ] = _mm_madd_epi16(m[6 ], di[6 ]);
        a[7 ] = _mm_madd_epi16(m[7 ], di[7 ]);

        a[8 ] = _mm_madd_epi16(m[8 ], di[0 ]);
        a[9 ] = _mm_madd_epi16(m[9 ], di[1 ]);
        a[10] = _mm_madd_epi16(m[10], di[2 ]);
        a[11] = _mm_madd_epi16(m[11], di[3 ]);
        a[12] = _mm_madd_epi16(m[12], di[4 ]);
        a[13] = _mm_madd_epi16(m[13], di[5 ]);
        a[14] = _mm_madd_epi16(m[14], di[6 ]);
        a[15] = _mm_madd_epi16(m[15], di[7 ]);

        b[0 ] = _mm_add_epi32(a[0 ], a[1 ]);
        b[1 ] = _mm_add_epi32(a[2 ], a[3 ]);
        b[2 ] = _mm_add_epi32(a[4 ], a[5 ]);
        b[3 ] = _mm_add_epi32(a[6 ], a[7 ]);

        b[4 ] = _mm_add_epi32(a[8 ], a[9 ]);
        b[5 ] = _mm_add_epi32(a[10], a[11]);
        b[6 ] = _mm_add_epi32(a[12], a[13]);
        b[7 ] = _mm_add_epi32(a[14], a[15]);

        c[0 ] = _mm_add_epi32(b[0 ], b[1 ]);
        c[1 ] = _mm_add_epi32(b[2 ], b[3 ]);

        c[2 ] = _mm_add_epi32(b[4 ], b[5 ]);
        c[3 ] = _mm_add_epi32(b[6 ], b[7 ]);

        r[0] = _mm_add_epi32(c[0], c[1]);
        r[1] = _mm_add_epi32(c[2], c[3]);

        r[0] = _mm_add_epi32(r[0], add);
        r[1] = _mm_add_epi32(r[1], add);

        r[0] = _mm_srai_epi32(r[0], shift);
        r[1] = _mm_srai_epi32(r[1], shift);

        res[k] = _mm_packs_epi32(r[0], r[1]);
      }
      transpose8x8(&res[0], &m[0]);
      transpose8x8(&res[8], &m[8]);

      _mm_store_si128((__m128i *)&_dst[0],   m[0 ]);
      _mm_store_si128((__m128i *)&_dst[16],  m[1 ]);
      _mm_store_si128((__m128i *)&_dst[32],  m[2 ]);
      _mm_store_si128((__m128i *)&_dst[48],  m[3 ]);
      _mm_store_si128((__m128i *)&_dst[64],  m[4 ]);
      _mm_store_si128((__m128i *)&_dst[80],  m[5 ]);
      _mm_store_si128((__m128i *)&_dst[96],  m[6 ]);
      _mm_store_si128((__m128i *)&_dst[112], m[7 ]);

      _mm_store_si128((__m128i *)&_dst[8],   m[8 ]);
      _mm_store_si128((__m128i *)&_dst[24],  m[9 ]);
      _mm_store_si128((__m128i *)&_dst[40],  m[10]);
      _mm_store_si128((__m128i *)&_dst[56],  m[11]);
      _mm_store_si128((__m128i *)&_dst[72],  m[12]);
      _mm_store_si128((__m128i *)&_dst[88],  m[13]);
      _mm_store_si128((__m128i *)&_dst[104], m[14]);
      _mm_store_si128((__m128i *)&_dst[120], m[15]);
      _src+=8;
      _dst+=128;
    }

    if (!(line&0x7)) {
      return ;
    }
    line = line&0x7;
    __m128i vhi[8], vlo[8], vh[16], vl[16], result[line][4];

    for(int l = 0; l < 2; ++l){
        for(int i = 0; i < line; ++i){
            result[i][2*l] = _mm_set1_epi32(0);
            result[i][2*l+1] = _mm_set1_epi32(0);
        }
    }

    for(int l = 0; l < 2; ++l){
        for(int k = 0; k < 2; ++k){
            for(int i = 0; i < 8; ++i){
                d[i] = _mm_load_si128((__m128i *)&(iT[i*16+k*128+l*8]));
            }

            for(int i = 0; i < line; ++i){

                for(int j = 0; j < 8; ++j){
                    x[j] = _mm_set1_epi16(_src[j*src_stride+i+k*8*src_stride]);

                    vhi[j] = _mm_mulhi_epi16(x[j], d[j]);
                    vlo[j] = _mm_mullo_epi16(x[j], d[j]);
                    vl[j+8*k] = _mm_unpacklo_epi16(vlo[j], vhi[j]);
                    vh[j+8*k] = _mm_unpackhi_epi16(vlo[j], vhi[j]);
                }

                for(int j = 0; j < 4; ++j){
                    result[i][2*l] = _mm_add_epi32(result[i][2*l], _mm_add_epi32(vl[2*j+8*k],vl[2*j+1+8*k]));
                    result[i][2*l+1] = _mm_add_epi32(result[i][2*l+1], _mm_add_epi32(vh[2*j+8*k],vh[2*j+1+8*k]));
                }
            }
        }
    }

    for(int l = 0; l < 2; ++l){
        for(int i = 0; i < line; ++i){
            result[i][2*l] = _mm_add_epi32(result[i][2*l], add);
            result[i][2*l+1] = _mm_add_epi32(result[i][2*l+1], add);
            result[i][2*l] = _mm_srai_epi32(result[i][2*l], shift);
            result[i][2*l+1] = _mm_srai_epi32(result[i][2*l+1], shift);
        }
    }

    for(int l = 0; l < 2; ++l){
        for(int i = 0; i < line; i++){
            result[i][l] = _mm_packs_epi32(result[i][2*l], result[i][2*l+1]); //clip pour repasser en 16
            _mm_store_si128((__m128i *)&(_dst[i*16+l*8]), result[i][l]);
        }
    }
}

void inverse_sse2_B32(const TCoeff *src, TCoeff *dst, int src_stride, int shift, int line, const TMatrixCoeff* iT){
    int i, k;
    TCoeff *_src = (TCoeff *)src;
    TCoeff *_dst = dst;
    __m128i x[32], m[32];
    __m128i d[32], di[16], a[32], b[16], c[8], r[2], res[32];
    __m128i add = _mm_set1_epi32(1 << (shift - 1));

    for (i = 0; i < line>>3; i++) {
      x[0 ]=_mm_load_si128((__m128i*)(_src +  0  * src_stride));
      x[1 ]=_mm_load_si128((__m128i*)(_src +  1  * src_stride));
      x[2 ]=_mm_load_si128((__m128i*)(_src +  2  * src_stride));
      x[3 ]=_mm_load_si128((__m128i*)(_src +  3  * src_stride));
      x[4 ]=_mm_load_si128((__m128i*)(_src +  4  * src_stride));
      x[5 ]=_mm_load_si128((__m128i*)(_src +  5  * src_stride));
      x[6 ]=_mm_load_si128((__m128i*)(_src +  6  * src_stride));
      x[7 ]=_mm_load_si128((__m128i*)(_src +  7  * src_stride));
      x[8 ]=_mm_load_si128((__m128i*)(_src +  8  * src_stride));
      x[9 ]=_mm_load_si128((__m128i*)(_src +  9  * src_stride));
      x[10]=_mm_load_si128((__m128i*)(_src +  10 * src_stride));
      x[11]=_mm_load_si128((__m128i*)(_src +  11 * src_stride));
      x[12]=_mm_load_si128((__m128i*)(_src +  12 * src_stride));
      x[13]=_mm_load_si128((__m128i*)(_src +  13 * src_stride));
      x[14]=_mm_load_si128((__m128i*)(_src +  14 * src_stride));
      x[15]=_mm_load_si128((__m128i*)(_src +  15 * src_stride));
      x[16]=_mm_load_si128((__m128i*)(_src +  16 * src_stride));
      x[17]=_mm_load_si128((__m128i*)(_src +  17 * src_stride));
      x[18]=_mm_load_si128((__m128i*)(_src +  18 * src_stride));
      x[19]=_mm_load_si128((__m128i*)(_src +  19 * src_stride));
      x[20]=_mm_load_si128((__m128i*)(_src +  20 * src_stride));
      x[21]=_mm_load_si128((__m128i*)(_src +  21 * src_stride));
      x[22]=_mm_load_si128((__m128i*)(_src +  22 * src_stride));
      x[23]=_mm_load_si128((__m128i*)(_src +  23 * src_stride));
      x[24]=_mm_load_si128((__m128i*)(_src +  24 * src_stride));
      x[25]=_mm_load_si128((__m128i*)(_src +  25 * src_stride));
      x[26]=_mm_load_si128((__m128i*)(_src +  26 * src_stride));
      x[27]=_mm_load_si128((__m128i*)(_src +  27 * src_stride));
      x[28]=_mm_load_si128((__m128i*)(_src +  28 * src_stride));
      x[29]=_mm_load_si128((__m128i*)(_src +  29 * src_stride));
      x[30]=_mm_load_si128((__m128i*)(_src +  30 * src_stride));
      x[31]=_mm_load_si128((__m128i*)(_src +  31 * src_stride));
      for (k = 0; k < 32; k++) {
        d[0 ] = _mm_set1_epi16(iT[k + 0  * 32]);
        d[1 ] = _mm_set1_epi16(iT[k + 1  * 32]);
        d[2 ] = _mm_set1_epi16(iT[k + 2  * 32]);
        d[3 ] = _mm_set1_epi16(iT[k + 3  * 32]);
        d[4 ] = _mm_set1_epi16(iT[k + 4  * 32]);
        d[5 ] = _mm_set1_epi16(iT[k + 5  * 32]);
        d[6 ] = _mm_set1_epi16(iT[k + 6  * 32]);
        d[7 ] = _mm_set1_epi16(iT[k + 7  * 32]);
        d[8 ] = _mm_set1_epi16(iT[k + 8  * 32]);
        d[9 ] = _mm_set1_epi16(iT[k + 9  * 32]);
        d[10] = _mm_set1_epi16(iT[k + 10 * 32]);
        d[11] = _mm_set1_epi16(iT[k + 11 * 32]);
        d[12] = _mm_set1_epi16(iT[k + 12 * 32]);
        d[13] = _mm_set1_epi16(iT[k + 13 * 32]);
        d[14] = _mm_set1_epi16(iT[k + 14 * 32]);
        d[15] = _mm_set1_epi16(iT[k + 15 * 32]);
        d[16] = _mm_set1_epi16(iT[k + 16 * 32]);
        d[17] = _mm_set1_epi16(iT[k + 17 * 32]);
        d[18] = _mm_set1_epi16(iT[k + 18 * 32]);
        d[19] = _mm_set1_epi16(iT[k + 19 * 32]);
        d[20] = _mm_set1_epi16(iT[k + 20 * 32]);
        d[21] = _mm_set1_epi16(iT[k + 21 * 32]);
        d[22] = _mm_set1_epi16(iT[k + 22 * 32]);
        d[23] = _mm_set1_epi16(iT[k + 23 * 32]);
        d[24] = _mm_set1_epi16(iT[k + 24 * 32]);
        d[25] = _mm_set1_epi16(iT[k + 25 * 32]);
        d[26] = _mm_set1_epi16(iT[k + 26 * 32]);
        d[27] = _mm_set1_epi16(iT[k + 27 * 32]);
        d[28] = _mm_set1_epi16(iT[k + 28 * 32]);
        d[29] = _mm_set1_epi16(iT[k + 29 * 32]);
        d[30] = _mm_set1_epi16(iT[k + 30 * 32]);
        d[31] = _mm_set1_epi16(iT[k + 31 * 32]);


        m[0 ] = _mm_unpacklo_epi16(x[0 ],  x[1 ]);
        m[1 ] = _mm_unpacklo_epi16(x[2 ],  x[3 ]);
        m[2 ] = _mm_unpacklo_epi16(x[4 ],  x[5 ]);
        m[3 ] = _mm_unpacklo_epi16(x[6 ],  x[7 ]);
        m[4 ] = _mm_unpacklo_epi16(x[8 ],  x[9 ]);
        m[5 ] = _mm_unpacklo_epi16(x[10],  x[11]);
        m[6 ] = _mm_unpacklo_epi16(x[12],  x[13]);
        m[7 ] = _mm_unpacklo_epi16(x[14],  x[15]);
        m[8 ] = _mm_unpacklo_epi16(x[16],  x[17]);
        m[9 ] = _mm_unpacklo_epi16(x[18],  x[19]);
        m[10] = _mm_unpacklo_epi16(x[20],  x[21]);
        m[11] = _mm_unpacklo_epi16(x[22],  x[23]);
        m[12] = _mm_unpacklo_epi16(x[24],  x[25]);
        m[13] = _mm_unpacklo_epi16(x[26],  x[27]);
        m[14] = _mm_unpacklo_epi16(x[28],  x[29]);
        m[15] = _mm_unpacklo_epi16(x[30],  x[31]);

        m[16] = _mm_unpackhi_epi16(x[0 ],  x[1 ]);
        m[17] = _mm_unpackhi_epi16(x[2 ],  x[3 ]);
        m[18] = _mm_unpackhi_epi16(x[4 ],  x[5 ]);
        m[19] = _mm_unpackhi_epi16(x[6 ],  x[7 ]);
        m[20] = _mm_unpackhi_epi16(x[8 ],  x[9 ]);
        m[21] = _mm_unpackhi_epi16(x[10],  x[11]);
        m[22] = _mm_unpackhi_epi16(x[12],  x[13]);
        m[23] = _mm_unpackhi_epi16(x[14],  x[15]);
        m[24] = _mm_unpackhi_epi16(x[16],  x[17]);
        m[25] = _mm_unpackhi_epi16(x[18],  x[19]);
        m[26] = _mm_unpackhi_epi16(x[20],  x[21]);
        m[27] = _mm_unpackhi_epi16(x[22],  x[23]);
        m[28] = _mm_unpackhi_epi16(x[24],  x[25]);
        m[29] = _mm_unpackhi_epi16(x[26],  x[27]);
        m[30] = _mm_unpackhi_epi16(x[28],  x[29]);
        m[31] = _mm_unpackhi_epi16(x[30],  x[31]);

        di[0 ] = _mm_unpacklo_epi16(d[0 ],  d[1 ]);
        di[1 ] = _mm_unpacklo_epi16(d[2 ],  d[3 ]);
        di[2 ] = _mm_unpacklo_epi16(d[4 ],  d[5 ]);
        di[3 ] = _mm_unpacklo_epi16(d[6 ],  d[7 ]);
        di[4 ] = _mm_unpacklo_epi16(d[8 ],  d[9 ]);
        di[5 ] = _mm_unpacklo_epi16(d[10],  d[11]);
        di[6 ] = _mm_unpacklo_epi16(d[12],  d[13]);
        di[7 ] = _mm_unpacklo_epi16(d[14],  d[15]);
        di[8 ] = _mm_unpacklo_epi16(d[16],  d[17]);
        di[9 ] = _mm_unpacklo_epi16(d[18],  d[19]);
        di[10] = _mm_unpacklo_epi16(d[20],  d[21]);
        di[11] = _mm_unpacklo_epi16(d[22],  d[23]);
        di[12] = _mm_unpacklo_epi16(d[24],  d[25]);
        di[13] = _mm_unpacklo_epi16(d[26],  d[27]);
        di[14] = _mm_unpacklo_epi16(d[28],  d[29]);
        di[15] = _mm_unpacklo_epi16(d[30],  d[31]);

        a[0 ] = _mm_madd_epi16(m[0 ], di[0 ]);
        a[1 ] = _mm_madd_epi16(m[1 ], di[1 ]);
        a[2 ] = _mm_madd_epi16(m[2 ], di[2 ]);
        a[3 ] = _mm_madd_epi16(m[3 ], di[3 ]);
        a[4 ] = _mm_madd_epi16(m[4 ], di[4 ]);
        a[5 ] = _mm_madd_epi16(m[5 ], di[5 ]);
        a[6 ] = _mm_madd_epi16(m[6 ], di[6 ]);
        a[7 ] = _mm_madd_epi16(m[7 ], di[7 ]);
        a[8 ] = _mm_madd_epi16(m[8 ], di[8 ]);
        a[9 ] = _mm_madd_epi16(m[9 ], di[9 ]);
        a[10] = _mm_madd_epi16(m[10], di[10]);
        a[11] = _mm_madd_epi16(m[11], di[11]);
        a[12] = _mm_madd_epi16(m[12], di[12]);
        a[13] = _mm_madd_epi16(m[13], di[13]);
        a[14] = _mm_madd_epi16(m[14], di[14]);
        a[15] = _mm_madd_epi16(m[15], di[15]);

        a[16] = _mm_madd_epi16(m[16], di[0 ]);
        a[17] = _mm_madd_epi16(m[17], di[1 ]);
        a[18] = _mm_madd_epi16(m[18], di[2 ]);
        a[19] = _mm_madd_epi16(m[19], di[3 ]);
        a[20] = _mm_madd_epi16(m[20], di[4 ]);
        a[21] = _mm_madd_epi16(m[21], di[5 ]);
        a[22] = _mm_madd_epi16(m[22], di[6 ]);
        a[23] = _mm_madd_epi16(m[23], di[7 ]);
        a[24] = _mm_madd_epi16(m[24], di[8 ]);
        a[25] = _mm_madd_epi16(m[25], di[9 ]);
        a[26] = _mm_madd_epi16(m[26], di[10]);
        a[27] = _mm_madd_epi16(m[27], di[11]);
        a[28] = _mm_madd_epi16(m[28], di[12]);
        a[29] = _mm_madd_epi16(m[29], di[13]);
        a[30] = _mm_madd_epi16(m[30], di[14]);
        a[31] = _mm_madd_epi16(m[31], di[15]);

        b[0 ] = _mm_add_epi32(a[0 ], a[1 ]);
        b[1 ] = _mm_add_epi32(a[2 ], a[3 ]);
        b[2 ] = _mm_add_epi32(a[4 ], a[5 ]);
        b[3 ] = _mm_add_epi32(a[6 ], a[7 ]);
        b[4 ] = _mm_add_epi32(a[8 ], a[9 ]);
        b[5 ] = _mm_add_epi32(a[10], a[11]);
        b[6 ] = _mm_add_epi32(a[12], a[13]);
        b[7 ] = _mm_add_epi32(a[14], a[15]);

        b[8 ] = _mm_add_epi32(a[16], a[17]);
        b[9 ] = _mm_add_epi32(a[18], a[19]);
        b[10] = _mm_add_epi32(a[20], a[21]);
        b[11] = _mm_add_epi32(a[22], a[23]);
        b[12] = _mm_add_epi32(a[24], a[25]);
        b[13] = _mm_add_epi32(a[26], a[27]);
        b[14] = _mm_add_epi32(a[28], a[29]);
        b[15] = _mm_add_epi32(a[30], a[31]);

        c[0 ] = _mm_add_epi32(b[0 ], b[1 ]);
        c[1 ] = _mm_add_epi32(b[2 ], b[3 ]);
        c[2 ] = _mm_add_epi32(b[4 ], b[5 ]);
        c[3 ] = _mm_add_epi32(b[6 ], b[7 ]);

        c[4 ] = _mm_add_epi32(b[8 ], b[9 ]);
        c[5 ] = _mm_add_epi32(b[10], b[11]);
        c[6 ] = _mm_add_epi32(b[12], b[13]);
        c[7 ] = _mm_add_epi32(b[14], b[15]);

        d[0 ] = _mm_add_epi32(c[0 ], c[1 ]);
        d[1 ] = _mm_add_epi32(c[2 ], c[3 ]);
        d[2 ] = _mm_add_epi32(c[4 ], c[5 ]);
        d[3 ] = _mm_add_epi32(c[6 ], c[7 ]);

        r[0] = _mm_add_epi32(d[0], d[1]);
        r[1] = _mm_add_epi32(d[2], d[3]);

        r[0] = _mm_add_epi32(r[0], add);
        r[1] = _mm_add_epi32(r[1], add);

        r[0] = _mm_srai_epi32(r[0], shift);
        r[1] = _mm_srai_epi32(r[1], shift);

        res[k] = _mm_packs_epi32(r[0], r[1]);
      }
      transpose8x8(&res[0], &m[0]);
      transpose8x8(&res[8], &m[8]);
      transpose8x8(&res[16], &m[16]);
      transpose8x8(&res[24], &m[24]);

      _mm_store_si128((__m128i *)&_dst[0],   m[0 ]);
      _mm_store_si128((__m128i *)&_dst[32],  m[1 ]);
      _mm_store_si128((__m128i *)&_dst[64],  m[2 ]);
      _mm_store_si128((__m128i *)&_dst[96],  m[3 ]);
      _mm_store_si128((__m128i *)&_dst[128], m[4 ]);
      _mm_store_si128((__m128i *)&_dst[160], m[5 ]);
      _mm_store_si128((__m128i *)&_dst[192], m[6 ]);
      _mm_store_si128((__m128i *)&_dst[224], m[7 ]);

      _mm_store_si128((__m128i *)&_dst[8],   m[8 ]);
      _mm_store_si128((__m128i *)&_dst[40],  m[9 ]);
      _mm_store_si128((__m128i *)&_dst[72],  m[10]);
      _mm_store_si128((__m128i *)&_dst[104], m[11]);
      _mm_store_si128((__m128i *)&_dst[136], m[12]);
      _mm_store_si128((__m128i *)&_dst[168], m[13]);
      _mm_store_si128((__m128i *)&_dst[200], m[14]);
      _mm_store_si128((__m128i *)&_dst[232], m[15]);

      _mm_store_si128((__m128i *)&_dst[16],  m[16]);
      _mm_store_si128((__m128i *)&_dst[48],  m[17]);
      _mm_store_si128((__m128i *)&_dst[80],  m[18]);
      _mm_store_si128((__m128i *)&_dst[112], m[19]);
      _mm_store_si128((__m128i *)&_dst[144], m[20]);
      _mm_store_si128((__m128i *)&_dst[176], m[21]);
      _mm_store_si128((__m128i *)&_dst[208], m[22]);
      _mm_store_si128((__m128i *)&_dst[240], m[23]);

      _mm_store_si128((__m128i *)&_dst[24],  m[24]);
      _mm_store_si128((__m128i *)&_dst[56],  m[25]);
      _mm_store_si128((__m128i *)&_dst[88],  m[26]);
      _mm_store_si128((__m128i *)&_dst[120], m[27]);
      _mm_store_si128((__m128i *)&_dst[152], m[28]);
      _mm_store_si128((__m128i *)&_dst[184], m[29]);
      _mm_store_si128((__m128i *)&_dst[216], m[30]);
      _mm_store_si128((__m128i *)&_dst[248], m[31]);

      _src+=8;
      _dst+=256;
    }

    if (!(line&0x7)) {
      return ;
    }
    line = line&0x7;
    __m128i vhi[8], vlo[8], vh[32], vl[32], result[line][8];

    for(int l = 0; l < 4; ++l){
        for(int i = 0; i < line; ++i){
            result[i][2*l] = _mm_set1_epi32(0);
            result[i][2*l+1] = _mm_set1_epi32(0);
        }
    }

    for(int l = 0; l < 4; ++l){
        for(int k = 0; k < 4; ++k){
            for(int i = 0; i < 8; ++i){
                d[i] = _mm_load_si128((__m128i *)&(iT[i*32+k*256+l*8]));
            }

            for(int i = 0; i < line; ++i){
                for(int j = 0; j < 8; ++j){
                    x[j] = _mm_set1_epi16(_src[j*src_stride+i+k*8*src_stride]);

                    vhi[j] = _mm_mulhi_epi16(x[j], d[j]);
                    vlo[j] = _mm_mullo_epi16(x[j], d[j]);
                    vl[j+8*k] = _mm_unpacklo_epi16(vlo[j], vhi[j]);
                    vh[j+8*k] = _mm_unpackhi_epi16(vlo[j], vhi[j]);
                }
                for(int j = 0; j < 4; ++j){
                    result[i][2*l] = _mm_add_epi32(result[i][2*l], _mm_add_epi32(vl[2*j+8*k],vl[2*j+1+8*k]));
                    result[i][2*l+1] = _mm_add_epi32(result[i][2*l+1], _mm_add_epi32(vh[2*j+8*k],vh[2*j+1+8*k]));
                }
            }
        }
    }

    for(int l = 0; l < 4; ++l){
        for(int i = 0; i < line; ++i){
            result[i][2*l] = _mm_add_epi32(result[i][2*l], add);
            result[i][2*l+1] = _mm_add_epi32(result[i][2*l+1], add);
            result[i][2*l] = _mm_srai_epi32(result[i][2*l], shift);
            result[i][2*l+1] = _mm_srai_epi32(result[i][2*l+1], shift);
        }
    }

    for(int l = 0; l < 4; ++l){
        for(int i = 0; i < line; i++){
            result[i][l] = _mm_packs_epi32(result[i][2*l], result[i][2*l+1]);
            _mm_store_si128((__m128i *)&(_dst[i*32+l*8]), result[i][l]);
        }
    }
}
