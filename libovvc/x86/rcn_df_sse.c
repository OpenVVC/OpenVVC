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
}
