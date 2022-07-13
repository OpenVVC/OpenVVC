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
*   Ibrahim FARHAT
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

#include "ovutils.h"
#include "data_rcn_transform.h"
#include "rcn_transform.h"
#include "arm/vvc_utils_sse.h"



//void ov_idct_ii_2_neon(const int16_t *src, int16_t *dst, ptrdiff_t src_stride, int num_lines, int line_brk, int shift);
void ov_idct_ii_4_neon(const int16_t *src, int16_t *dst, ptrdiff_t src_stride, int shift, int shift_add, const int16_t *Mat);

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
dct2_4x2(__m128i *x, __m128i *r)
{
  __m128i m[2], a[2], d[2];

  static const int16_t DCT_II_4_2_sse[8 * 2]= {
    64, 64, 64, 64, 83, 36, 83, 36,
    64,-64, 64,-64, 36,-83, 36,-83
  };

  d[0] = _mm_load_si128((__m128i*)(DCT_II_4_2_sse+0));
  d[1] = _mm_load_si128((__m128i*)(DCT_II_4_2_sse+8));

  transform_4x2(x,d,m);

  a[0] = _mm_add_epi32(m[0], m[1]);
  a[1] = _mm_sub_epi32(m[0], m[1]);

  r[0] = _mm_unpacklo_epi64(a[0], a[1]);
  r[1] = _mm_unpackhi_epi64(a[0], a[1]);
}

static inline void
dct2_4x8(__m128i *x, __m128i *r)
{
  __m128i m[8], d[8];

  d[0] = _mm_set1_epi16(DCT_II_8[0]);
  d[1] = _mm_set1_epi16(DCT_II_8[1]);
  d[2] = _mm_set1_epi16(DCT_II_8[16]);
  d[3] = _mm_set1_epi16(DCT_II_8[17]);
  d[4] = _mm_set1_epi16(DCT_II_8[32]);
  d[5] = _mm_set1_epi16(DCT_II_8[33]);
  d[6] = _mm_set1_epi16(DCT_II_8[48]);
  d[7] = _mm_set1_epi16(DCT_II_8[49]);

  transform_4x8(x, d, m);

  transpose4x4_step2(m  , r  );
  transpose4x4_step2(m+1, r+4);
}

void vvc_inverse_dct_ii_2_neon(const int16_t *src, int16_t *dst, ptrdiff_t src_stride,
                          int num_lines, int line_brk, int shift)
{
  const int16_t *_src = (const int16_t *)src;
  uint16_t *_dst = dst;
  //ov_idct_ii_2_neon(_dst, _src, src_stride, num_lines, line_brk, shift);

}


void
vvc_inverse_dct_ii_4_neon(const int16_t *src, int16_t *dst, ptrdiff_t src_stride,
                          int num_lines, int line_brk, int shift)
{
  static const int16_t DCT_II_4_4_neon[4 * 4] = {
    64, 64, 64, 64, 83, 36, -36,-83,
    64, -64, -64, 64, 36, -83, 83, -36
  };
  static const int16_t DCT_II_8_4_neon[4 * 8] = {
    64,  64,  64,  64,  83,  36,  83,  36,
    64, -64,  64, -64,  36, -83,  36, -83,
    64,  64,  64,  64,  36,  83,  36,  83,
    -64, 64, -64,  64, -83,  36, -83,  36,
  };

  __m128i add = _mm_set1_epi32(1 << (shift - 1));

  for (int j = 0; j < num_lines / 8; j++) {
    __m128i x[4], d[8], r[8];

    x[0] = _mm_load_si128((__m128i*)(src + 0 * src_stride));
    x[1] = _mm_load_si128((__m128i*)(src + 1 * src_stride));
    x[2] = _mm_load_si128((__m128i*)(src + 2 * src_stride));
    x[3] = _mm_load_si128((__m128i*)(src + 3 * src_stride));

    dct2_4x8(x, r);

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

    ov_idct_ii_4_neon(src, dst, src_stride<<1, -shift,(1<<shift-1), DCT_II_4_4_neon); // functional

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

    dct2_4x2(x, r);

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


void rcn_init_tr_functions_neon(struct RCNFunctions *const rcn_funcs){
  /* rcn_funcs->tr.func[DST_VII][2] = &vvc_inverse_dst_vii_4_sse;
  rcn_funcs->tr.func[DST_VII][3] = &vvc_inverse_dst_vii_8_sse;
  rcn_funcs->tr.func[DST_VII][4] = &vvc_inverse_dst_vii_16_sse;
  rcn_funcs->tr.func[DST_VII][5] = &vvc_inverse_dst_vii_32_sse;
  */
  /*rcn_funcs->tr.func[DCT_VIII][2] = &vvc_inverse_dct_viii_4_sse;
  rcn_funcs->tr.func[DCT_VIII][3] = &vvc_inverse_dct_viii_8_sse;
  rcn_funcs->tr.func[DCT_VIII][4] = &vvc_inverse_dct_viii_16_sse;
  rcn_funcs->tr.func[DCT_VIII][5] = &vvc_inverse_dct_viii_32_sse;
   */
  //rcn_funcs->tr.func[DCT_II][1] = &vvc_inverse_dct_ii_2_neon;
  rcn_funcs->tr.func[DCT_II][2] = &vvc_inverse_dct_ii_4_neon;
  /*rcn_funcs->tr.func[DCT_II][3] = &vvc_inverse_dct_ii_8_sse;
  rcn_funcs->tr.func[DCT_II][4] = &vvc_inverse_dct_ii_16_sse;
  rcn_funcs->tr.func[DCT_II][5] = &vvc_inverse_dct_ii_32_sse;
  rcn_funcs->tr.func[DCT_II][6] = &vvc_inverse_dct_ii_64_sse;

  rcn_funcs->tr.dc = &vvc_inverse_dct_ii_dc_sse;*/
}
