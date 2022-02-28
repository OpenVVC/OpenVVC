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

#include <stdint.h>
#include <tmmintrin.h>
#include <smmintrin.h>

#include "ovutils.h"
#include "rcn_structures.h"

static void
filter_h_7_7(OVSample *src, const int stride, const int tc)
{
    static const int16_t db7_p[8] = { 0, 5, 14, 23, 32, 41, 50, 59 };
    static const int16_t tc7_p[8] = { 0, 1, 1, 2, 3, 4, 5, 6, };
    static const int16_t db7_q[8] = { 59, 50, 41, 32, 23, 14, 5, 0 };
    static const int16_t tc7_q[8] = { 6, 5, 4, 3, 2, 1, 1, 0 };

    __m128i tmp = _mm_set1_epi16(tc);
    __m128i clp_p = _mm_loadu_si128((__m128i *)tc7_p);
    clp_p = _mm_mullo_epi16(clp_p, tmp);
    clp_p = _mm_srli_epi16(clp_p, 1);

    __m128i clp_q = _mm_loadu_si128((__m128i *)tc7_q);
    clp_q = _mm_mullo_epi16(clp_q, tmp);
    clp_q = _mm_srli_epi16(clp_q, 1);

    __m128i db7p = _mm_loadu_si128((__m128i *)db7_p);
    __m128i db7q = _mm_loadu_si128((__m128i *)db7_q);


    int16_t ref_p[4];
    int16_t ref_q[4];
    int16_t db_ref[4];

    OVSample* _src = src;
    for (int i = 0; i < 4; ++i) {
        OVSample* srcP = _src - 1;
        OVSample* srcQ = _src;
        ref_p[i] = (srcP[-6 * 1] + srcP[-7 * 1] + 1) >> 1;
        ref_q[i] = (srcQ[ 6 * 1] + srcQ[ 7 * 1] + 1) >> 1;

        db_ref[i] = (2 * (srcP[0] + srcQ[0])
                     + srcP[-1] + srcP[-2 * 1] + srcP[-3 * 1] + srcP[-4 * 1] + srcP[-5 * 1] + srcP[-6 * 1]
                     + srcQ[ 1] + srcQ[ 2 * 1] + srcQ[ 3 * 1] + srcQ[ 4 * 1] + srcQ[ 5 * 1] + srcQ[ 6 * 1] + 8) >> 4;
        _src += stride;
    }
    __m128i add_32 = _mm_set1_epi16(32);

    for (int i = 0; i < 4; i++) {
        OVSample* _src = src - 8;
        __m128i p = _mm_load_si128((__m128i *) _src);
        __m128i q = _mm_load_si128((__m128i *)  src);
        __m128i clp_min_p = _mm_sub_epi16(p, clp_p);
        __m128i clp_max_p = _mm_add_epi16(p, clp_p);
        __m128i clp_min_q = _mm_sub_epi16(q, clp_q);
        __m128i clp_max_q = _mm_add_epi16(q, clp_q);
        __m128i dref = _mm_set1_epi16(db_ref[i]);
        __m128i rf_p = _mm_set1_epi16(ref_p[i]);
        __m128i rf_q = _mm_set1_epi16(ref_q[i]);

        __m128i dref_p = _mm_mullo_epi16(dref, db7p);
        __m128i dref_q = _mm_mullo_epi16(dref, db7q);
        __m128i x1 = _mm_mullo_epi16(rf_p, db7p);
        __m128i x2 = _mm_mullo_epi16(rf_q, db7q);

        __m128i x0 = _mm_subs_epu16(_mm_slli_epi16(rf_p, 6), x1);
        __m128i x3 = _mm_subs_epu16(_mm_slli_epi16(rf_q, 6), x2);

        x0 = _mm_adds_epu16(x0, dref_p);
        x0 = _mm_adds_epu16(x0, add_32);
        x3 = _mm_adds_epu16(x3, dref_q);
        x3 = _mm_adds_epu16(x3, add_32);

        x3 = _mm_srli_epi16(x3, 6);
        x0 = _mm_srli_epi16(x0, 6);

        x0 = _mm_max_epi16(x0, clp_min_p);
        x3 = _mm_max_epi16(x3, clp_min_q);
        x0 = _mm_min_epi16(x0, clp_max_p);
        x3 = _mm_min_epi16(x3, clp_max_q);

        _mm_store_si128((__m128i *) _src, x0);
        _mm_store_si128((__m128i *) src, x3);

        src += stride;
    }
}

static void
filter_v_7_7(OVSample *src, const int stride, const int tc)
{
    static const int dbCoeffs7[7] = { 59, 50, 41, 32, 23, 14, 5 };
    static const int8_t tc7[7] = { 6, 5, 4, 3, 2, 1, 1};

    int16_t* srcP = src - stride;
    int16_t* srcQ = src;

    __m128i srcP0  = _mm_loadl_epi64( ( const __m128i* ) &srcP[0] );
    __m128i srcP1  = _mm_loadl_epi64( ( const __m128i* ) &srcP[-stride] );
    __m128i srcP2  = _mm_loadl_epi64( ( const __m128i* ) &srcP[stride * -2] );
    __m128i srcP3  = _mm_loadl_epi64( ( const __m128i* ) &srcP[stride * -3] );
    __m128i srcP4  = _mm_loadl_epi64( ( const __m128i* ) &srcP[stride * -4] );
    __m128i srcP5  = _mm_loadl_epi64( ( const __m128i* ) &srcP[stride * -5] );
    __m128i srcP6  = _mm_loadl_epi64( ( const __m128i* ) &srcP[stride * -6] );
    __m128i srcP7  = _mm_loadl_epi64( ( const __m128i* ) &srcP[stride * -7] );

    __m128i srcQ0  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[0] );
    __m128i srcQ1  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride] );
    __m128i srcQ2  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * 2] );
    __m128i srcQ3  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * 3] );
    __m128i srcQ4  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * 4] );
    __m128i srcQ5  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * 5] );
    __m128i srcQ6  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * 6] );
    __m128i srcQ7  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * 7] );

    __m128i refPx4 = _mm_add_epi16(srcP6, srcP7);
    __m128i refQx4 = _mm_add_epi16(srcQ6, srcQ7);
    refPx4 = _mm_add_epi16(refPx4, _mm_set1_epi16(1));
    refQx4 = _mm_add_epi16(refQx4, _mm_set1_epi16(1));
    refPx4 = _mm_srai_epi16(refPx4, 1);
    refQx4 = _mm_srai_epi16(refQx4, 1);

    srcP0  = _mm_cvtepi16_epi32(srcP0);
    srcP1  = _mm_cvtepi16_epi32(srcP1);
    srcP2  = _mm_cvtepi16_epi32(srcP2);
    srcP3  = _mm_cvtepi16_epi32(srcP3);
    srcP4  = _mm_cvtepi16_epi32(srcP4);
    srcP5  = _mm_cvtepi16_epi32(srcP5);
    srcP6  = _mm_cvtepi16_epi32(srcP6);

    srcQ0  = _mm_cvtepi16_epi32(srcQ0);
    srcQ1  = _mm_cvtepi16_epi32(srcQ1);
    srcQ2  = _mm_cvtepi16_epi32(srcQ2);
    srcQ3  = _mm_cvtepi16_epi32(srcQ3);
    srcQ4  = _mm_cvtepi16_epi32(srcQ4);
    srcQ5  = _mm_cvtepi16_epi32(srcQ5);
    srcQ6  = _mm_cvtepi16_epi32(srcQ6);

    srcP0 = _mm_add_epi32(srcP0, srcQ0);
    srcP1 = _mm_add_epi32(srcP1, srcQ1);
    srcP2 = _mm_add_epi32(srcP2, srcQ2);
    srcP3 = _mm_add_epi32(srcP3, srcQ3);
    srcP4 = _mm_add_epi32(srcP4, srcQ4);
    srcP5 = _mm_add_epi32(srcP5, srcQ5);
    srcP6 = _mm_add_epi32(srcP6, srcQ6);

    srcP0 = _mm_slli_epi32(srcP0, 1);
    srcP1 = _mm_add_epi32(srcP1, srcP2);
    srcP3 = _mm_add_epi32(srcP3, srcP4);
    srcP5 = _mm_add_epi32(srcP5, srcP6);

    srcP0 = _mm_add_epi32(srcP0, srcP1);
    srcP3 = _mm_add_epi32(srcP3, srcP5);

    srcP0 = _mm_add_epi32(srcP0, srcP3);
    __m128i refMiddlex4 = _mm_add_epi32(srcP0, _mm_set1_epi32(8));
    refMiddlex4 = _mm_srai_epi32(refMiddlex4, 4);
    refMiddlex4 = _mm_packs_epi32(refMiddlex4, _mm_setzero_si128());

   for( int pos = 0; pos < 7; pos++ )
   {
     __m128i vref1 = refPx4;
     __m128i vref0 = refMiddlex4;
     __m128i vsrc  = _mm_loadl_epi64( ( const __m128i* ) &srcP[-stride * pos] );
     __m128i vmax  = _mm_set1_epi16( ( tc * tc7[pos] ) >> 1 );
     __m128i vmin  = _mm_sub_epi16( vsrc, vmax );
     vmax          = _mm_add_epi16( vsrc, vmax );
     vref0         = _mm_unpacklo_epi16( vref0, vref1 );
     __m128i vtmp  = _mm_set1_epi32( dbCoeffs7[pos] | ( ( 64 - dbCoeffs7[pos] ) << 16 ) );
     vtmp          = _mm_madd_epi16( vref0, vtmp );
     vtmp          = _mm_add_epi32( vtmp, _mm_set1_epi32( 32 ) );
     vtmp          = _mm_srli_epi32( vtmp, 6 );
     vtmp          = _mm_packs_epi32( vtmp, vtmp );
     vtmp          = _mm_min_epi16( _mm_max_epi16( vtmp, vmin ), vmax );
     _mm_storel_epi64( ( __m128i* ) &srcP[-stride * pos], vtmp );
   }

   for( int pos = 0; pos < 7; pos++ )
   {
     __m128i vref1 = refQx4;
     __m128i vref0 = refMiddlex4;
     __m128i vsrc  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * pos] );
     __m128i vmax  = _mm_set1_epi16( ( tc * tc7[pos] ) >> 1 );
     __m128i vmin  = _mm_sub_epi16( vsrc, vmax );
     vmax          = _mm_add_epi16( vsrc, vmax );
     vref0         = _mm_unpacklo_epi16( vref0, vref1 );
     __m128i vtmp  = _mm_set1_epi32( dbCoeffs7[pos] | ( ( 64 - dbCoeffs7[pos] ) << 16 ) );
     vtmp          = _mm_madd_epi16( vref0, vtmp );
     vtmp          = _mm_add_epi32( vtmp, _mm_set1_epi32( 32 ) );
     vtmp          = _mm_srli_epi32( vtmp, 6 );
     vtmp          = _mm_packs_epi32( vtmp, vtmp );
     vtmp          = _mm_min_epi16( _mm_max_epi16( vtmp, vmin ), vmax );
     _mm_storel_epi64( ( __m128i* ) &srcQ[stride * pos], vtmp );
   }
}

static void
filter_v_7_5(OVSample *src, const int stride, const int tc)
{
      static const int dbCoeffs7[7] = { 59, 50, 41, 32, 23, 14, 5 };
      static const int dbCoeffs5[5] = { 58, 45, 32, 19, 6};
      static const int8_t tc7[7] = { 6, 5, 4, 3, 2, 1, 1};

      int16_t* srcP = src - stride;
      int16_t* srcQ = src;

      __m128i srcP0  = _mm_loadl_epi64( ( const __m128i* ) &srcP[0] );
      __m128i srcP1  = _mm_loadl_epi64( ( const __m128i* ) &srcP[-stride] );
      __m128i srcP2  = _mm_loadl_epi64( ( const __m128i* ) &srcP[stride * -2] );
      __m128i srcP3  = _mm_loadl_epi64( ( const __m128i* ) &srcP[stride * -3] );
      __m128i srcP4  = _mm_loadl_epi64( ( const __m128i* ) &srcP[stride * -4] );
      __m128i srcP5  = _mm_loadl_epi64( ( const __m128i* ) &srcP[stride * -5] );
      __m128i srcP6  = _mm_loadl_epi64( ( const __m128i* ) &srcP[stride * -6] );
      __m128i srcP7  = _mm_loadl_epi64( ( const __m128i* ) &srcP[stride * -7] );

      __m128i srcQ0  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[0] );
      __m128i srcQ1  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride] );
      __m128i srcQ2  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * 2] );
      __m128i srcQ3  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * 3] );
      __m128i srcQ4  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * 4] );
      __m128i srcQ5  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * 5] );

      __m128i refPx4 = _mm_add_epi16(srcP6, srcP7);
      __m128i refQx4 = _mm_add_epi16(srcQ4, srcQ5);
      refPx4 = _mm_add_epi16(refPx4, _mm_set1_epi16(1));
      refQx4 = _mm_add_epi16(refQx4, _mm_set1_epi16(1));
      refPx4 = _mm_srai_epi16(refPx4, 1);
      refQx4 = _mm_srai_epi16(refQx4, 1);

      srcP0  = _mm_cvtepi16_epi32(srcP0);
      srcP1  = _mm_cvtepi16_epi32(srcP1);
      srcP2  = _mm_cvtepi16_epi32(srcP2);
      srcP3  = _mm_cvtepi16_epi32(srcP3);
      srcP4  = _mm_cvtepi16_epi32(srcP4);
      srcP5  = _mm_cvtepi16_epi32(srcP5);

      srcQ0  = _mm_cvtepi16_epi32(srcQ0);
      srcQ1  = _mm_cvtepi16_epi32(srcQ1);
      srcQ2  = _mm_cvtepi16_epi32(srcQ2);
      srcQ3  = _mm_cvtepi16_epi32(srcQ3);
      srcQ4  = _mm_cvtepi16_epi32(srcQ4);
      srcQ5  = _mm_cvtepi16_epi32(srcQ5);

      srcP0 = _mm_add_epi32(srcP0, srcQ0);
      srcP1 = _mm_add_epi32(srcP1, srcQ1);
      srcP2 = _mm_add_epi32(srcP2, srcQ2);
      srcP3 = _mm_add_epi32(srcP3, srcQ3);
      srcP4 = _mm_add_epi32(srcP4, srcQ4);
      srcP5 = _mm_add_epi32(srcP5, srcQ5);

      srcP0 = _mm_add_epi32(srcP0, srcP1);
      srcP2 = _mm_add_epi32(srcP2, srcP3);
      srcP4 = _mm_add_epi32(srcP4, srcP5);
      srcP0 = _mm_slli_epi32(srcP0, 1);

      srcP0 = _mm_add_epi32(srcP0, srcP2);

      srcP0 = _mm_add_epi32(srcP0, srcP4);
      __m128i refMiddlex4 = _mm_add_epi32(srcP0, _mm_set1_epi32(8));
      refMiddlex4 = _mm_srai_epi32(refMiddlex4, 4);
      refMiddlex4 = _mm_packs_epi32(refMiddlex4, _mm_setzero_si128());

     for( int pos = 0; pos < 7; pos++ )
     {
       __m128i vref1 = refPx4;
       __m128i vref0 = refMiddlex4;
       __m128i vsrc  = _mm_loadl_epi64( ( const __m128i* ) &srcP[-stride * pos] );
       __m128i vmax  = _mm_set1_epi16( ( tc * tc7[pos] ) >> 1 );
       __m128i vmin  = _mm_sub_epi16( vsrc, vmax );
       vmax          = _mm_add_epi16( vsrc, vmax );
       vref0         = _mm_unpacklo_epi16( vref0, vref1 );
       __m128i vtmp  = _mm_set1_epi32( dbCoeffs7[pos] | ( ( 64 - dbCoeffs7[pos] ) << 16 ) );
       vtmp          = _mm_madd_epi16( vref0, vtmp );
       vtmp          = _mm_add_epi32( vtmp, _mm_set1_epi32( 32 ) );
       vtmp          = _mm_srli_epi32( vtmp, 6 );
       vtmp          = _mm_packs_epi32( vtmp, vtmp );
       vtmp          = _mm_min_epi16( _mm_max_epi16( vtmp, vmin ), vmax );
       _mm_storel_epi64( ( __m128i* ) &srcP[-stride * pos], vtmp );
     }

     for( int pos = 0; pos < 5; pos++ )
     {
       __m128i vref1 = refQx4;
       __m128i vref0 = refMiddlex4;
       __m128i vsrc  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * pos] );
       __m128i vmax  = _mm_set1_epi16( ( tc * tc7[pos] ) >> 1 );
       __m128i vmin  = _mm_sub_epi16( vsrc, vmax );
       vmax          = _mm_add_epi16( vsrc, vmax );
       vref0         = _mm_unpacklo_epi16( vref0, vref1 );
       __m128i vtmp  = _mm_set1_epi32( dbCoeffs5[pos] | ( ( 64 - dbCoeffs5[pos] ) << 16 ) );
       vtmp          = _mm_madd_epi16( vref0, vtmp );
       vtmp          = _mm_add_epi32( vtmp, _mm_set1_epi32( 32 ) );
       vtmp          = _mm_srli_epi32( vtmp, 6 );
       vtmp          = _mm_packs_epi32( vtmp, vtmp );
       vtmp          = _mm_min_epi16( _mm_max_epi16( vtmp, vmin ), vmax );
       _mm_storel_epi64( ( __m128i* ) &srcQ[stride * pos], vtmp );
     }
}

static void
filter_v_5_7(OVSample *src, const int stride, const int tc)
{
  static const int dbCoeffs7[7] = { 59, 50, 41, 32, 23, 14, 5 };
  static const int dbCoeffs5[5] = { 58, 45, 32, 19, 6};
  static const int8_t tc7[7] = { 6, 5, 4, 3, 2, 1, 1};

  int16_t* srcP = src - stride;
  int16_t* srcQ = src;

  __m128i srcP0  = _mm_loadl_epi64( ( const __m128i* ) &srcP[0] );
  __m128i srcP1  = _mm_loadl_epi64( ( const __m128i* ) &srcP[-stride] );
  __m128i srcP2  = _mm_loadl_epi64( ( const __m128i* ) &srcP[stride * -2] );
  __m128i srcP3  = _mm_loadl_epi64( ( const __m128i* ) &srcP[stride * -3] );
  __m128i srcP4  = _mm_loadl_epi64( ( const __m128i* ) &srcP[stride * -4] );
  __m128i srcP5  = _mm_loadl_epi64( ( const __m128i* ) &srcP[stride * -5] );

  __m128i srcQ0  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[0] );
  __m128i srcQ1  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride] );
  __m128i srcQ2  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * 2] );
  __m128i srcQ3  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * 3] );
  __m128i srcQ4  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * 4] );
  __m128i srcQ5  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * 5] );
  __m128i srcQ6  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * 6] );
  __m128i srcQ7  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * 7] );

  __m128i refPx4 = _mm_add_epi16(srcP4, srcP5);
  __m128i refQx4 = _mm_add_epi16(srcQ6, srcQ7);
  refPx4 = _mm_add_epi16(refPx4, _mm_set1_epi16(1));
  refQx4 = _mm_add_epi16(refQx4, _mm_set1_epi16(1));
  refPx4 = _mm_srai_epi16(refPx4, 1);
  refQx4 = _mm_srai_epi16(refQx4, 1);

  srcP0  = _mm_cvtepi16_epi32(srcP0);
  srcP1  = _mm_cvtepi16_epi32(srcP1);
  srcP2  = _mm_cvtepi16_epi32(srcP2);
  srcP3  = _mm_cvtepi16_epi32(srcP3);
  srcP4  = _mm_cvtepi16_epi32(srcP4);
  srcP5  = _mm_cvtepi16_epi32(srcP5);

  srcQ0  = _mm_cvtepi16_epi32(srcQ0);
  srcQ1  = _mm_cvtepi16_epi32(srcQ1);
  srcQ2  = _mm_cvtepi16_epi32(srcQ2);
  srcQ3  = _mm_cvtepi16_epi32(srcQ3);
  srcQ4  = _mm_cvtepi16_epi32(srcQ4);
  srcQ5  = _mm_cvtepi16_epi32(srcQ5);

  srcP0 = _mm_add_epi32(srcP0, srcQ0);
  srcP1 = _mm_add_epi32(srcP1, srcQ1);
  srcP2 = _mm_add_epi32(srcP2, srcQ2);
  srcP3 = _mm_add_epi32(srcP3, srcQ3);
  srcP4 = _mm_add_epi32(srcP4, srcQ4);
  srcP5 = _mm_add_epi32(srcP5, srcQ5);

  srcP0 = _mm_add_epi32(srcP0, srcP1);
  srcP2 = _mm_add_epi32(srcP2, srcP3);
  srcP4 = _mm_add_epi32(srcP4, srcP5);
  srcP0 = _mm_slli_epi32(srcP0, 1);

  srcP0 = _mm_add_epi32(srcP0, srcP2);

  srcP0 = _mm_add_epi32(srcP0, srcP4);
  __m128i refMiddlex4 = _mm_add_epi32(srcP0, _mm_set1_epi32(8));
  refMiddlex4 = _mm_srai_epi32(refMiddlex4, 4);
  refMiddlex4 = _mm_packs_epi32(refMiddlex4, _mm_setzero_si128());

 for( int pos = 0; pos < 5; pos++ )
 {
   __m128i vref1 = refPx4;
   __m128i vref0 = refMiddlex4;
   __m128i vsrc  = _mm_loadl_epi64( ( const __m128i* ) &srcP[-stride * pos] );
   __m128i vmax  = _mm_set1_epi16( ( tc * tc7[pos] ) >> 1 );
   __m128i vmin  = _mm_sub_epi16( vsrc, vmax );
   vmax          = _mm_add_epi16( vsrc, vmax );
   vref0         = _mm_unpacklo_epi16( vref0, vref1 );
   __m128i vtmp  = _mm_set1_epi32( dbCoeffs5[pos] | ( ( 64 - dbCoeffs5[pos] ) << 16 ) );
   vtmp          = _mm_madd_epi16( vref0, vtmp );
   vtmp          = _mm_add_epi32( vtmp, _mm_set1_epi32( 32 ) );
   vtmp          = _mm_srli_epi32( vtmp, 6 );
   vtmp          = _mm_packs_epi32( vtmp, vtmp );
   vtmp          = _mm_min_epi16( _mm_max_epi16( vtmp, vmin ), vmax );
   _mm_storel_epi64( ( __m128i* ) &srcP[-stride * pos], vtmp );
 }

 for( int pos = 0; pos < 7; pos++ )
 {
   __m128i vref1 = refQx4;
   __m128i vref0 = refMiddlex4;
   __m128i vsrc  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * pos] );
   __m128i vmax  = _mm_set1_epi16( ( tc * tc7[pos] ) >> 1 );
   __m128i vmin  = _mm_sub_epi16( vsrc, vmax );
   vmax          = _mm_add_epi16( vsrc, vmax );
   vref0         = _mm_unpacklo_epi16( vref0, vref1 );
   __m128i vtmp  = _mm_set1_epi32( dbCoeffs7[pos] | ( ( 64 - dbCoeffs7[pos] ) << 16 ) );
   vtmp          = _mm_madd_epi16( vref0, vtmp );
   vtmp          = _mm_add_epi32( vtmp, _mm_set1_epi32( 32 ) );
   vtmp          = _mm_srli_epi32( vtmp, 6 );
   vtmp          = _mm_packs_epi32( vtmp, vtmp );
   vtmp          = _mm_min_epi16( _mm_max_epi16( vtmp, vmin ), vmax );
   _mm_storel_epi64( ( __m128i* ) &srcQ[stride * pos], vtmp );
 }

}

static void
filter_v_5_5(OVSample *src, const int stride, const int tc)
{
  static const int dbCoeffs5[5] = { 58, 45, 32, 19, 6};
  static const int8_t tc7[7] = { 6, 5, 4, 3, 2, 1, 1};

  int16_t* srcP = src - stride;
  int16_t* srcQ = src;

  __m128i srcP0  = _mm_loadl_epi64( ( const __m128i* ) &srcP[0] );
  __m128i srcP1  = _mm_loadl_epi64( ( const __m128i* ) &srcP[-stride] );
  __m128i srcP2  = _mm_loadl_epi64( ( const __m128i* ) &srcP[stride * -2] );
  __m128i srcP3  = _mm_loadl_epi64( ( const __m128i* ) &srcP[stride * -3] );
  __m128i srcP4  = _mm_loadl_epi64( ( const __m128i* ) &srcP[stride * -4] );
  __m128i srcP5  = _mm_loadl_epi64( ( const __m128i* ) &srcP[stride * -5] );

  __m128i srcQ0  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[0] );
  __m128i srcQ1  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride] );
  __m128i srcQ2  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * 2] );
  __m128i srcQ3  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * 3] );
  __m128i srcQ4  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * 4] );
  __m128i srcQ5  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * 5] );

  __m128i refPx4 = _mm_add_epi16(srcP4, srcP5);
  __m128i refQx4 = _mm_add_epi16(srcQ4, srcQ5);
  refPx4 = _mm_add_epi16(refPx4, _mm_set1_epi16(1));
  refQx4 = _mm_add_epi16(refQx4, _mm_set1_epi16(1));
  refPx4 = _mm_srai_epi16(refPx4, 1);
  refQx4 = _mm_srai_epi16(refQx4, 1);

  srcP0  = _mm_cvtepi16_epi32(srcP0);
  srcP1  = _mm_cvtepi16_epi32(srcP1);
  srcP2  = _mm_cvtepi16_epi32(srcP2);
  srcP3  = _mm_cvtepi16_epi32(srcP3);
  srcP4  = _mm_cvtepi16_epi32(srcP4);

  srcQ0  = _mm_cvtepi16_epi32(srcQ0);
  srcQ1  = _mm_cvtepi16_epi32(srcQ1);
  srcQ2  = _mm_cvtepi16_epi32(srcQ2);
  srcQ3  = _mm_cvtepi16_epi32(srcQ3);
  srcQ4  = _mm_cvtepi16_epi32(srcQ4);

  srcP0 = _mm_add_epi32(srcP0, srcQ0);
  srcP1 = _mm_add_epi32(srcP1, srcQ1);
  srcP2 = _mm_add_epi32(srcP2, srcQ2);
  srcP3 = _mm_add_epi32(srcP3, srcQ3);
  srcP4 = _mm_add_epi32(srcP4, srcQ4);

  srcP0 = _mm_add_epi32(srcP0, srcP1);
  srcP0 = _mm_add_epi32(srcP0, srcP2);
  srcP3 = _mm_add_epi32(srcP3, srcP4);
  srcP0 = _mm_slli_epi32(srcP0, 1);

  srcP0 = _mm_add_epi32(srcP0, srcP3);

  __m128i refMiddlex4 = _mm_add_epi32(srcP0, _mm_set1_epi32(8));
  refMiddlex4 = _mm_srai_epi32(refMiddlex4, 4);
  refMiddlex4 = _mm_packs_epi32(refMiddlex4, _mm_setzero_si128());

 for( int pos = 0; pos < 5; pos++ )
 {
   __m128i vref1 = refPx4;
   __m128i vref0 = refMiddlex4;
   __m128i vsrc  = _mm_loadl_epi64( ( const __m128i* ) &srcP[-stride * pos] );
   __m128i vmax  = _mm_set1_epi16( ( tc * tc7[pos] ) >> 1 );
   __m128i vmin  = _mm_sub_epi16( vsrc, vmax );
   vmax          = _mm_add_epi16( vsrc, vmax );
   vref0         = _mm_unpacklo_epi16( vref0, vref1 );
   __m128i vtmp  = _mm_set1_epi32( dbCoeffs5[pos] | ( ( 64 - dbCoeffs5[pos] ) << 16 ) );
   vtmp          = _mm_madd_epi16( vref0, vtmp );
   vtmp          = _mm_add_epi32( vtmp, _mm_set1_epi32( 32 ) );
   vtmp          = _mm_srli_epi32( vtmp, 6 );
   vtmp          = _mm_packs_epi32( vtmp, vtmp );
   vtmp          = _mm_min_epi16( _mm_max_epi16( vtmp, vmin ), vmax );
   _mm_storel_epi64( ( __m128i* ) &srcP[-stride * pos], vtmp );
 }

 for( int pos = 0; pos < 5; pos++ )
 {
   __m128i vref1 = refQx4;
   __m128i vref0 = refMiddlex4;
   __m128i vsrc  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * pos] );
   __m128i vmax  = _mm_set1_epi16( ( tc * tc7[pos] ) >> 1 );
   __m128i vmin  = _mm_sub_epi16( vsrc, vmax );
   vmax          = _mm_add_epi16( vsrc, vmax );
   vref0         = _mm_unpacklo_epi16( vref0, vref1 );
   __m128i vtmp  = _mm_set1_epi32( dbCoeffs5[pos] | ( ( 64 - dbCoeffs5[pos] ) << 16 ) );
   vtmp          = _mm_madd_epi16( vref0, vtmp );
   vtmp          = _mm_add_epi32( vtmp, _mm_set1_epi32( 32 ) );
   vtmp          = _mm_srli_epi32( vtmp, 6 );
   vtmp          = _mm_packs_epi32( vtmp, vtmp );
   vtmp          = _mm_min_epi16( _mm_max_epi16( vtmp, vmin ), vmax );
   _mm_storel_epi64( ( __m128i* ) &srcQ[stride * pos], vtmp );
 }
}

static void
filter_v_7_3(OVSample *src, const int stride, const int tc)
{
  static const int dbCoeffs7[7] = { 59, 50, 41, 32, 23, 14, 5 };
  static const int dbCoeffs3[3] = { 53, 32, 11 };
  static const int8_t tc7[7] = { 6, 5, 4, 3, 2, 1, 1};
  static const int8_t tc3[3] = { 6, 4, 2 };
  int16_t* srcP = src - stride;
  int16_t* srcQ = src;

  __m128i srcP0  = _mm_loadl_epi64( ( const __m128i* ) &srcP[0] );
  __m128i srcP1  = _mm_loadl_epi64( ( const __m128i* ) &srcP[-stride] );
  __m128i srcP2  = _mm_loadl_epi64( ( const __m128i* ) &srcP[stride * -2] );
  __m128i srcP3  = _mm_loadl_epi64( ( const __m128i* ) &srcP[stride * -3] );
  __m128i srcP4  = _mm_loadl_epi64( ( const __m128i* ) &srcP[stride * -4] );
  __m128i srcP5  = _mm_loadl_epi64( ( const __m128i* ) &srcP[stride * -5] );
  __m128i srcP6  = _mm_loadl_epi64( ( const __m128i* ) &srcP[stride * -6] );
  __m128i srcP7  = _mm_loadl_epi64( ( const __m128i* ) &srcP[stride * -7] );

  __m128i srcQ0  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[0] );
  __m128i srcQ1  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride] );
  __m128i srcQ2  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * 2] );
  __m128i srcQ3  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * 3] );

  __m128i refPx4 = _mm_add_epi16(srcP6, srcP7);
  __m128i refQx4 = _mm_add_epi16(srcQ2, srcQ3);
  refPx4 = _mm_add_epi16(refPx4, _mm_set1_epi16(1));
  refQx4 = _mm_add_epi16(refQx4, _mm_set1_epi16(1));
  refPx4 = _mm_srai_epi16(refPx4, 1);
  refQx4 = _mm_srai_epi16(refQx4, 1);

  srcP0  = _mm_cvtepi16_epi32(srcP0);
  srcP1  = _mm_cvtepi16_epi32(srcP1);
  srcP2  = _mm_cvtepi16_epi32(srcP2);
  srcP3  = _mm_cvtepi16_epi32(srcP3);
  srcP4  = _mm_cvtepi16_epi32(srcP4);
  srcP5  = _mm_cvtepi16_epi32(srcP5);
  srcP6  = _mm_cvtepi16_epi32(srcP6);

  srcQ0  = _mm_cvtepi16_epi32(srcQ0);
  srcQ1  = _mm_cvtepi16_epi32(srcQ1);
  srcQ2  = _mm_cvtepi16_epi32(srcQ2);


  srcP0 = _mm_add_epi32(srcP0, srcQ0);
  srcP1 = _mm_add_epi32(srcP1, srcQ0);
  srcP2 = _mm_add_epi32(srcP2, srcQ1);
  srcP3 = _mm_add_epi32(srcP3, srcQ1);
  srcP4 = _mm_add_epi32(srcP4, srcQ1);
  srcP5 = _mm_add_epi32(srcP5, srcQ2);
  srcP6 = _mm_add_epi32(srcP6, srcQ2);

  srcP0 = _mm_slli_epi32(srcP0, 1);
  srcP1 = _mm_add_epi32(srcP1, srcP2);
  srcP3 = _mm_add_epi32(srcP3, srcP4);
  srcP5 = _mm_add_epi32(srcP5, srcP6);

  srcP0 = _mm_add_epi32(srcP0, srcP1);
  srcP3 = _mm_add_epi32(srcP3, srcP5);

  srcP0 = _mm_add_epi32(srcP0, srcP3);
  __m128i refMiddlex4 = _mm_add_epi32(srcP0, _mm_set1_epi32(8));
  refMiddlex4 = _mm_srai_epi32(refMiddlex4, 4);
  refMiddlex4 = _mm_packs_epi32(refMiddlex4, _mm_setzero_si128());

 for( int pos = 0; pos < 7; pos++ )
 {
   __m128i vref1 = refPx4;
   __m128i vref0 = refMiddlex4;
   __m128i vsrc  = _mm_loadl_epi64( ( const __m128i* ) &srcP[-stride * pos] );
   __m128i vmax  = _mm_set1_epi16( ( tc * tc7[pos] ) >> 1 );
   __m128i vmin  = _mm_sub_epi16( vsrc, vmax );
   vmax          = _mm_add_epi16( vsrc, vmax );
   vref0         = _mm_unpacklo_epi16( vref0, vref1 );
   __m128i vtmp  = _mm_set1_epi32( dbCoeffs7[pos] | ( ( 64 - dbCoeffs7[pos] ) << 16 ) );
   vtmp          = _mm_madd_epi16( vref0, vtmp );
   vtmp          = _mm_add_epi32( vtmp, _mm_set1_epi32( 32 ) );
   vtmp          = _mm_srli_epi32( vtmp, 6 );
   vtmp          = _mm_packs_epi32( vtmp, vtmp );
   vtmp          = _mm_min_epi16( _mm_max_epi16( vtmp, vmin ), vmax );
   _mm_storel_epi64( ( __m128i* ) &srcP[-stride * pos], vtmp );
 }

 for( int pos = 0; pos < 3; pos++ )
 {
   __m128i vref1 = refQx4;
   __m128i vref0 = refMiddlex4;
   __m128i vsrc  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * pos] );
   __m128i vmax  = _mm_set1_epi16( ( tc * tc3[pos] ) >> 1 );
   __m128i vmin  = _mm_sub_epi16( vsrc, vmax );
   vmax          = _mm_add_epi16( vsrc, vmax );
   vref0         = _mm_unpacklo_epi16( vref0, vref1 );
   __m128i vtmp  = _mm_set1_epi32( dbCoeffs3[pos] | ( ( 64 - dbCoeffs3[pos] ) << 16 ) );
   vtmp          = _mm_madd_epi16( vref0, vtmp );
   vtmp          = _mm_add_epi32( vtmp, _mm_set1_epi32( 32 ) );
   vtmp          = _mm_srli_epi32( vtmp, 6 );
   vtmp          = _mm_packs_epi32( vtmp, vtmp );
   vtmp          = _mm_min_epi16( _mm_max_epi16( vtmp, vmin ), vmax );
   _mm_storel_epi64( ( __m128i* ) &srcQ[stride * pos], vtmp );
 }
}

static void
filter_v_3_7(OVSample *src, const int stride, const int tc)
{
  static const int dbCoeffs7[7] = { 59, 50, 41, 32, 23, 14, 5 };
  static const int dbCoeffs3[3] = { 53, 32, 11 };
  static const int8_t tc7[7] = { 6, 5, 4, 3, 2, 1, 1};
  static const int8_t tc3[3] = { 6, 4, 2 };

  int16_t* srcP = src - stride;
  int16_t* srcQ = src;

  __m128i srcP0  = _mm_loadl_epi64( ( const __m128i* ) &srcP[0] );
  __m128i srcP1  = _mm_loadl_epi64( ( const __m128i* ) &srcP[-stride] );
  __m128i srcP2  = _mm_loadl_epi64( ( const __m128i* ) &srcP[stride * -2] );
  __m128i srcP3  = _mm_loadl_epi64( ( const __m128i* ) &srcP[stride * -3] );

  __m128i srcQ0  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[0] );
  __m128i srcQ1  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride] );
  __m128i srcQ2  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * 2] );
  __m128i srcQ3  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * 3] );
  __m128i srcQ4  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * 4] );
  __m128i srcQ5  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * 5] );
  __m128i srcQ6  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * 6] );
  __m128i srcQ7  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * 7] );

  __m128i refPx4 = _mm_add_epi16(srcP2, srcP3);
  __m128i refQx4 = _mm_add_epi16(srcQ6, srcQ7);
  refPx4 = _mm_add_epi16(refPx4, _mm_set1_epi16(1));
  refQx4 = _mm_add_epi16(refQx4, _mm_set1_epi16(1));
  refPx4 = _mm_srai_epi16(refPx4, 1);
  refQx4 = _mm_srai_epi16(refQx4, 1);

  srcP0  = _mm_cvtepi16_epi32(srcP0);
  srcP1  = _mm_cvtepi16_epi32(srcP1);
  srcP2  = _mm_cvtepi16_epi32(srcP2);

  srcQ0  = _mm_cvtepi16_epi32(srcQ0);
  srcQ1  = _mm_cvtepi16_epi32(srcQ1);
  srcQ2  = _mm_cvtepi16_epi32(srcQ2);
  srcQ3  = _mm_cvtepi16_epi32(srcQ3);
  srcQ4  = _mm_cvtepi16_epi32(srcQ4);
  srcQ5  = _mm_cvtepi16_epi32(srcQ5);
  srcQ6  = _mm_cvtepi16_epi32(srcQ6);

  srcQ0 = _mm_add_epi32(srcP0, srcQ0);
  srcQ1 = _mm_add_epi32(srcP0, srcQ1);
  srcQ2 = _mm_add_epi32(srcP1, srcQ2);
  srcQ3 = _mm_add_epi32(srcP1, srcQ3);
  srcQ4 = _mm_add_epi32(srcP1, srcQ4);
  srcQ5 = _mm_add_epi32(srcP2, srcQ5);
  srcQ6 = _mm_add_epi32(srcP2, srcQ6);

  srcQ0 = _mm_slli_epi32(srcQ0, 1);
  srcQ1 = _mm_add_epi32(srcQ1, srcQ2);
  srcQ3 = _mm_add_epi32(srcQ3, srcQ4);
  srcQ5 = _mm_add_epi32(srcQ5, srcQ6);

  srcQ0 = _mm_add_epi32(srcQ0, srcQ1);
  srcQ3 = _mm_add_epi32(srcQ3, srcQ5);

  srcQ0 = _mm_add_epi32(srcQ0, srcQ3);
  __m128i refMiddlex4 = _mm_add_epi32(srcQ0, _mm_set1_epi32(8));
  refMiddlex4 = _mm_srai_epi32(refMiddlex4, 4);
  refMiddlex4 = _mm_packs_epi32(refMiddlex4, _mm_setzero_si128());

 for( int pos = 0; pos < 3; pos++ )
 {
   __m128i vref1 = refPx4;
   __m128i vref0 = refMiddlex4;
   __m128i vsrc  = _mm_loadl_epi64( ( const __m128i* ) &srcP[-stride * pos] );
   __m128i vmax  = _mm_set1_epi16( ( tc * tc3[pos] ) >> 1 );
   __m128i vmin  = _mm_sub_epi16( vsrc, vmax );
   vmax          = _mm_add_epi16( vsrc, vmax );
   vref0         = _mm_unpacklo_epi16( vref0, vref1 );
   __m128i vtmp  = _mm_set1_epi32( dbCoeffs3[pos] | ( ( 64 - dbCoeffs3[pos] ) << 16 ) );
   vtmp          = _mm_madd_epi16( vref0, vtmp );
   vtmp          = _mm_add_epi32( vtmp, _mm_set1_epi32( 32 ) );
   vtmp          = _mm_srli_epi32( vtmp, 6 );
   vtmp          = _mm_packs_epi32( vtmp, vtmp );
   vtmp          = _mm_min_epi16( _mm_max_epi16( vtmp, vmin ), vmax );
   _mm_storel_epi64( ( __m128i* ) &srcP[-stride * pos], vtmp );
 }

 for( int pos = 0; pos < 7; pos++ )
 {
   __m128i vref1 = refQx4;
   __m128i vref0 = refMiddlex4;
   __m128i vsrc  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * pos] );
   __m128i vmax  = _mm_set1_epi16( ( tc * tc7[pos] ) >> 1 );
   __m128i vmin  = _mm_sub_epi16( vsrc, vmax );
   vmax          = _mm_add_epi16( vsrc, vmax );
   vref0         = _mm_unpacklo_epi16( vref0, vref1 );
   __m128i vtmp  = _mm_set1_epi32( dbCoeffs7[pos] | ( ( 64 - dbCoeffs7[pos] ) << 16 ) );
   vtmp          = _mm_madd_epi16( vref0, vtmp );
   vtmp          = _mm_add_epi32( vtmp, _mm_set1_epi32( 32 ) );
   vtmp          = _mm_srli_epi32( vtmp, 6 );
   vtmp          = _mm_packs_epi32( vtmp, vtmp );
   vtmp          = _mm_min_epi16( _mm_max_epi16( vtmp, vmin ), vmax );
   _mm_storel_epi64( ( __m128i* ) &srcQ[stride * pos], vtmp );
 }
}

static void
filter_v_5_3(OVSample *src, const int stride, const int tc)
{
  static const int dbCoeffs3[3] = { 53, 32, 11 };
  static const int dbCoeffs5[5] = { 58, 45, 32, 19, 6};
  static const int8_t tc7[7] = { 6, 5, 4, 3, 2, 1, 1};
  static const int8_t tc3[3] = { 6, 4, 2 };

  int16_t* srcP = src - stride;
  int16_t* srcQ = src;

  __m128i srcP0  = _mm_loadl_epi64( ( const __m128i* ) &srcP[0] );
  __m128i srcP1  = _mm_loadl_epi64( ( const __m128i* ) &srcP[-stride] );
  __m128i srcP2  = _mm_loadl_epi64( ( const __m128i* ) &srcP[stride * -2] );
  __m128i srcP3  = _mm_loadl_epi64( ( const __m128i* ) &srcP[stride * -3] );
  __m128i srcP4  = _mm_loadl_epi64( ( const __m128i* ) &srcP[stride * -4] );
  __m128i srcP5  = _mm_loadl_epi64( ( const __m128i* ) &srcP[stride * -5] );

  __m128i srcQ0  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[0] );
  __m128i srcQ1  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride] );
  __m128i srcQ2  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * 2] );
  __m128i srcQ3  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * 3] );


  __m128i refPx4 = _mm_add_epi16(srcP4, srcP5);
  __m128i refQx4 = _mm_add_epi16(srcQ2, srcQ3);
  refPx4 = _mm_add_epi16(refPx4, _mm_set1_epi16(1));
  refQx4 = _mm_add_epi16(refQx4, _mm_set1_epi16(1));
  refPx4 = _mm_srai_epi16(refPx4, 1);
  refQx4 = _mm_srai_epi16(refQx4, 1);

  srcP0  = _mm_cvtepi16_epi32(srcP0);
  srcP1  = _mm_cvtepi16_epi32(srcP1);
  srcP2  = _mm_cvtepi16_epi32(srcP2);
  srcP3  = _mm_cvtepi16_epi32(srcP3);

  srcQ0  = _mm_cvtepi16_epi32(srcQ0);
  srcQ1  = _mm_cvtepi16_epi32(srcQ1);
  srcQ2  = _mm_cvtepi16_epi32(srcQ2);
  srcQ3  = _mm_cvtepi16_epi32(srcQ3);

  srcP0 = _mm_add_epi32(srcP0, srcQ0);
  srcP1 = _mm_add_epi32(srcP1, srcQ1);
  srcP2 = _mm_add_epi32(srcP2, srcQ2);
  srcP3 = _mm_add_epi32(srcP3, srcQ3);

  srcP0 = _mm_add_epi32(srcP0, srcP1);
  srcP2 = _mm_add_epi32(srcP2, srcP3);

  srcP0 = _mm_add_epi32(srcP0, srcP2);

  __m128i refMiddlex4 = _mm_add_epi32(srcP0, _mm_set1_epi32(4));
  refMiddlex4 = _mm_srai_epi32(refMiddlex4, 3);
  refMiddlex4 = _mm_packs_epi32(refMiddlex4, _mm_setzero_si128());

 for( int pos = 0; pos < 5; pos++ )
 {
   __m128i vref1 = refPx4;
   __m128i vref0 = refMiddlex4;
   __m128i vsrc  = _mm_loadl_epi64( ( const __m128i* ) &srcP[-stride * pos] );
   __m128i vmax  = _mm_set1_epi16( ( tc * tc7[pos] ) >> 1 );
   __m128i vmin  = _mm_sub_epi16( vsrc, vmax );
   vmax          = _mm_add_epi16( vsrc, vmax );
   vref0         = _mm_unpacklo_epi16( vref0, vref1 );
   __m128i vtmp  = _mm_set1_epi32( dbCoeffs5[pos] | ( ( 64 - dbCoeffs5[pos] ) << 16 ) );
   vtmp          = _mm_madd_epi16( vref0, vtmp );
   vtmp          = _mm_add_epi32( vtmp, _mm_set1_epi32( 32 ) );
   vtmp          = _mm_srli_epi32( vtmp, 6 );
   vtmp          = _mm_packs_epi32( vtmp, vtmp );
   vtmp          = _mm_min_epi16( _mm_max_epi16( vtmp, vmin ), vmax );
   _mm_storel_epi64( ( __m128i* ) &srcP[-stride * pos], vtmp );
 }

 for( int pos = 0; pos < 3; pos++ )
 {
   __m128i vref1 = refQx4;
   __m128i vref0 = refMiddlex4;
   __m128i vsrc  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * pos] );
   __m128i vmax  = _mm_set1_epi16( ( tc * tc3[pos] ) >> 1 );
   __m128i vmin  = _mm_sub_epi16( vsrc, vmax );
   vmax          = _mm_add_epi16( vsrc, vmax );
   vref0         = _mm_unpacklo_epi16( vref0, vref1 );
   __m128i vtmp  = _mm_set1_epi32( dbCoeffs3[pos] | ( ( 64 - dbCoeffs3[pos] ) << 16 ) );
   vtmp          = _mm_madd_epi16( vref0, vtmp );
   vtmp          = _mm_add_epi32( vtmp, _mm_set1_epi32( 32 ) );
   vtmp          = _mm_srli_epi32( vtmp, 6 );
   vtmp          = _mm_packs_epi32( vtmp, vtmp );
   vtmp          = _mm_min_epi16( _mm_max_epi16( vtmp, vmin ), vmax );
   _mm_storel_epi64( ( __m128i* ) &srcQ[stride * pos], vtmp );
 }
}

static void
filter_v_3_5(OVSample *src, const int stride, const int tc)
{
  static const int dbCoeffs3[3] = { 53, 32, 11 };
  static const int dbCoeffs5[5] = { 58, 45, 32, 19, 6};
  static const int8_t tc7[7] = { 6, 5, 4, 3, 2, 1, 1};
  static const int8_t tc3[3] = { 6, 4, 2 };

  int16_t* srcP = src - stride;
  int16_t* srcQ = src;

  __m128i srcP0  = _mm_loadl_epi64( ( const __m128i* ) &srcP[0] );
  __m128i srcP1  = _mm_loadl_epi64( ( const __m128i* ) &srcP[-stride] );
  __m128i srcP2  = _mm_loadl_epi64( ( const __m128i* ) &srcP[stride * -2] );
  __m128i srcP3  = _mm_loadl_epi64( ( const __m128i* ) &srcP[stride * -3] );

  __m128i srcQ0  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[0] );
  __m128i srcQ1  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride] );
  __m128i srcQ2  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * 2] );
  __m128i srcQ3  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * 3] );
  __m128i srcQ4  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * 4] );
  __m128i srcQ5  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * 5] );


  __m128i refPx4 = _mm_add_epi16(srcP2, srcP3);
  __m128i refQx4 = _mm_add_epi16(srcQ4, srcQ5);
  refPx4 = _mm_add_epi16(refPx4, _mm_set1_epi16(1));
  refQx4 = _mm_add_epi16(refQx4, _mm_set1_epi16(1));
  refPx4 = _mm_srai_epi16(refPx4, 1);
  refQx4 = _mm_srai_epi16(refQx4, 1);

  srcP0  = _mm_cvtepi16_epi32(srcP0);
  srcP1  = _mm_cvtepi16_epi32(srcP1);
  srcP2  = _mm_cvtepi16_epi32(srcP2);
  srcP3  = _mm_cvtepi16_epi32(srcP3);

  srcQ0  = _mm_cvtepi16_epi32(srcQ0);
  srcQ1  = _mm_cvtepi16_epi32(srcQ1);
  srcQ2  = _mm_cvtepi16_epi32(srcQ2);
  srcQ3  = _mm_cvtepi16_epi32(srcQ3);

  srcP0 = _mm_add_epi32(srcP0, srcQ0);
  srcP1 = _mm_add_epi32(srcP1, srcQ1);
  srcP2 = _mm_add_epi32(srcP2, srcQ2);
  srcP3 = _mm_add_epi32(srcP3, srcQ3);

  srcP0 = _mm_add_epi32(srcP0, srcP1);
  srcP2 = _mm_add_epi32(srcP2, srcP3);

  srcP0 = _mm_add_epi32(srcP0, srcP2);

  __m128i refMiddlex4 = _mm_add_epi32(srcP0, _mm_set1_epi32(4));
  refMiddlex4 = _mm_srai_epi32(refMiddlex4, 3);
  refMiddlex4 = _mm_packs_epi32(refMiddlex4, _mm_setzero_si128());

 for( int pos = 0; pos < 3; pos++ )
 {
   __m128i vref1 = refPx4;
   __m128i vref0 = refMiddlex4;
   __m128i vsrc  = _mm_loadl_epi64( ( const __m128i* ) &srcP[-stride * pos] );
   __m128i vmax  = _mm_set1_epi16( ( tc * tc3[pos] ) >> 1 );
   __m128i vmin  = _mm_sub_epi16( vsrc, vmax );
   vmax          = _mm_add_epi16( vsrc, vmax );
   vref0         = _mm_unpacklo_epi16( vref0, vref1 );
   __m128i vtmp  = _mm_set1_epi32( dbCoeffs3[pos] | ( ( 64 - dbCoeffs3[pos] ) << 16 ) );
   vtmp          = _mm_madd_epi16( vref0, vtmp );
   vtmp          = _mm_add_epi32( vtmp, _mm_set1_epi32( 32 ) );
   vtmp          = _mm_srli_epi32( vtmp, 6 );
   vtmp          = _mm_packs_epi32( vtmp, vtmp );
   vtmp          = _mm_min_epi16( _mm_max_epi16( vtmp, vmin ), vmax );
   _mm_storel_epi64( ( __m128i* ) &srcP[-stride * pos], vtmp );
 }

 for( int pos = 0; pos < 5; pos++ )
 {
   __m128i vref1 = refQx4;
   __m128i vref0 = refMiddlex4;
   __m128i vsrc  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * pos] );
   __m128i vmax  = _mm_set1_epi16( ( tc * tc7[pos] ) >> 1 );
   __m128i vmin  = _mm_sub_epi16( vsrc, vmax );
   vmax          = _mm_add_epi16( vsrc, vmax );
   vref0         = _mm_unpacklo_epi16( vref0, vref1 );
   __m128i vtmp  = _mm_set1_epi32( dbCoeffs5[pos] | ( ( 64 - dbCoeffs5[pos] ) << 16 ) );
   vtmp          = _mm_madd_epi16( vref0, vtmp );
   vtmp          = _mm_add_epi32( vtmp, _mm_set1_epi32( 32 ) );
   vtmp          = _mm_srli_epi32( vtmp, 6 );
   vtmp          = _mm_packs_epi32( vtmp, vtmp );
   vtmp          = _mm_min_epi16( _mm_max_epi16( vtmp, vmin ), vmax );
   _mm_storel_epi64( ( __m128i* ) &srcQ[stride * pos], vtmp );
 }
}

void
rcn_init_df_functions_sse(struct RCNFunctions *const rcn_funcs)
{
  rcn_funcs->df.filter_v[0] = NULL;
  rcn_funcs->df.filter_v[1] = &filter_v_3_5;
  rcn_funcs->df.filter_v[2] = &filter_v_3_7;
  rcn_funcs->df.filter_v[3] = NULL;
  rcn_funcs->df.filter_v[4] = &filter_v_5_3;
  rcn_funcs->df.filter_v[5] = &filter_v_5_5;
  rcn_funcs->df.filter_v[6] = &filter_v_5_7;
  rcn_funcs->df.filter_v[7] = NULL;
  rcn_funcs->df.filter_v[8] = &filter_v_7_3;
  rcn_funcs->df.filter_v[9] = &filter_v_7_5;
  rcn_funcs->df.filter_v[10]= &filter_v_7_7;
  rcn_funcs->df.filter_h[10]= &filter_h_7_7;
}
