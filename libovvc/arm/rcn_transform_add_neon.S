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

#include "arm/asm.S"

.macro CLIP, val, min_val, max_val
   SMAX  \val, \val, \min_val
   SMIN  \val, \val, \max_val
.endm

function transform_add_8_4_10_neon, export=1

  mov x15, x0
  mov w10, #0x03FF
  dup v30.8h, w10
  dup v31.8h, wzr

  ld1 {v13.8h}, [x15], x1
  ld1 {v14.8h}, [x15], x1
  ld1 {v15.8h}, [x15], x1
  ld1 {v16.8h}, [x15]

  ld1 {v17.8h}, [x2], x3
  ld1 {v18.8h}, [x2], x3
  ld1 {v19.8h}, [x2], x3
  ld1 {v20.8h}, [x2]

  add v13.8h, v13.8h, v17.8h
  add v14.8h, v14.8h, v18.8h
  add v15.8h, v15.8h, v19.8h
  add v16.8h, v16.8h, v20.8h

  CLIP v13.8h, v31.8h, v30.8h
  CLIP v14.8h, v31.8h, v30.8h
  CLIP v15.8h, v31.8h, v30.8h
  CLIP v16.8h, v31.8h, v30.8h

  st1 {v13.8h}, [x0], x1
  st1 {v14.8h}, [x0], x1
  st1 {v15.8h}, [x0], x1
  st1 {v16.8h}, [x0]

  ret lr
endfunc

function transform_add_16_2_10_neon, export=1

  mov x15, x0
  mov x16, x2
  mov w10, #0x03FF
  dup v30.8h, w10
  dup v31.8h, wzr

  ld1 {v13.8h}, [x15], x1
  ld1 {v14.8h}, [x15]
  add x15, x0, #16
  ld1 {v15.8h}, [x15], x1
  ld1 {v16.8h}, [x15]

  ld1 {v17.8h}, [x16], x3
  ld1 {v18.8h}, [x16]
  add x16, x2, #16
  ld1 {v19.8h}, [x16], x3
  ld1 {v20.8h}, [x16]

  add v13.8h, v13.8h, v17.8h
  add v14.8h, v14.8h, v18.8h
  add v15.8h, v15.8h, v19.8h
  add v16.8h, v16.8h, v20.8h

  CLIP v13.8h, v31.8h, v30.8h
  CLIP v14.8h, v31.8h, v30.8h
  CLIP v15.8h, v31.8h, v30.8h
  CLIP v16.8h, v31.8h, v30.8h

  mov x15, x0
  st1 {v13.8h}, [x15], x1
  st1 {v14.8h}, [x15]
  add x15, x0, #16
  st1 {v15.8h}, [x15], x1
  st1 {v16.8h}, [x15]

  ret lr
endfunc

function transform_add_32_1_10_neon, export=1

  mov x15, x0
  mov w10, #0x03FF
  dup v30.8h, w10
  dup v31.8h, wzr

  ld1 {v13.8h}, [x15], #16
  ld1 {v14.8h}, [x15], #16
  ld1 {v15.8h}, [x15], #16
  ld1 {v16.8h}, [x15]

  ld1 {v17.8h}, [x2], #16
  ld1 {v18.8h}, [x2], #16
  ld1 {v19.8h}, [x2], #16
  ld1 {v20.8h}, [x2]

  add v13.8h, v13.8h, v17.8h
  add v14.8h, v14.8h, v18.8h
  add v15.8h, v15.8h, v19.8h
  add v16.8h, v16.8h, v20.8h

  CLIP v13.8h, v31.8h, v30.8h
  CLIP v14.8h, v31.8h, v30.8h
  CLIP v15.8h, v31.8h, v30.8h
  CLIP v16.8h, v31.8h, v30.8h

  st1 {v13.8h}, [x0], #16
  st1 {v14.8h}, [x0], #16
  st1 {v15.8h}, [x0], #16
  st1 {v16.8h}, [x0]

  ret lr
endfunc

function transform_add_half_8_4_10_neon, export=1

  mov x15, x0
  mov w10, #0x03FF
  dup v30.8h, w10
  dup v31.8h, wzr

  ld1 {v13.8h}, [x15], x1
  ld1 {v14.8h}, [x15], x1
  ld1 {v15.8h}, [x15], x1
  ld1 {v16.8h}, [x15]

  ld1 {v17.8h}, [x2], x3
  ld1 {v18.8h}, [x2], x3
  ld1 {v19.8h}, [x2], x3
  ld1 {v20.8h}, [x2]

  sshr v17.8h, v17.8h, #1
  sshr v18.8h, v18.8h, #1
  sshr v19.8h, v19.8h, #1
  sshr v20.8h, v20.8h, #1

  add v13.8h, v13.8h, v17.8h
  add v14.8h, v14.8h, v18.8h
  add v15.8h, v15.8h, v19.8h
  add v16.8h, v16.8h, v20.8h

  CLIP v13.8h, v31.8h, v30.8h
  CLIP v14.8h, v31.8h, v30.8h
  CLIP v15.8h, v31.8h, v30.8h
  CLIP v16.8h, v31.8h, v30.8h

  st1 {v13.8h}, [x0], x1
  st1 {v14.8h}, [x0], x1
  st1 {v15.8h}, [x0], x1
  st1 {v16.8h}, [x0]

  ret lr
endfunc


function transform_add_half_16_2_10_neon, export=1

   mov x15, x0
   mov x16, x2
   mov w10, #0x03FF
   dup v30.8h, w10
   dup v31.8h, wzr

   ld1 {v13.8h}, [x15], x1
   ld1 {v14.8h}, [x15]
   add x15, x0, #16
   ld1 {v15.8h}, [x15], x1
   ld1 {v16.8h}, [x15]

   ld1 {v17.8h}, [x16], x3
   ld1 {v18.8h}, [x16]
   add x16, x2, #16
   ld1 {v19.8h}, [x16], x3
   ld1 {v20.8h}, [x16]

   sshr v17.8h, v17.8h, #1
   sshr v18.8h, v18.8h, #1
   sshr v19.8h, v19.8h, #1
   sshr v20.8h, v20.8h, #1

   add v13.8h, v13.8h, v17.8h
   add v14.8h, v14.8h, v18.8h
   add v15.8h, v15.8h, v19.8h
   add v16.8h, v16.8h, v20.8h

   CLIP v13.8h, v31.8h, v30.8h
   CLIP v14.8h, v31.8h, v30.8h
   CLIP v15.8h, v31.8h, v30.8h
   CLIP v16.8h, v31.8h, v30.8h

   mov x15, x0
   st1 {v13.8h}, [x15], x1
   st1 {v14.8h}, [x15]
   add x15, x0, #16
   st1 {v15.8h}, [x15], x1
   st1 {v16.8h}, [x15]

   ret lr
endfunc

function transform_add_half_32_1_10_neon, export=1

  mov x15, x0
  mov w10, #0x03FF
  dup v30.8h, w10
  dup v31.8h, wzr

  ld1 {v13.8h}, [x15], #16
  ld1 {v14.8h}, [x15], #16
  ld1 {v15.8h}, [x15], #16
  ld1 {v16.8h}, [x15]

  ld1 {v17.8h}, [x2], #16
  ld1 {v18.8h}, [x2], #16
  ld1 {v19.8h}, [x2], #16
  ld1 {v20.8h}, [x2]

  sshr v17.8h, v17.8h, #1
  sshr v18.8h, v18.8h, #1
  sshr v19.8h, v19.8h, #1
  sshr v20.8h, v20.8h, #1

  add v13.8h, v13.8h, v17.8h
  add v14.8h, v14.8h, v18.8h
  add v15.8h, v15.8h, v19.8h
  add v16.8h, v16.8h, v20.8h

  CLIP v13.8h, v31.8h, v30.8h
  CLIP v14.8h, v31.8h, v30.8h
  CLIP v15.8h, v31.8h, v30.8h
  CLIP v16.8h, v31.8h, v30.8h

  st1 {v13.8h}, [x0], #16
  st1 {v14.8h}, [x0], #16
  st1 {v15.8h}, [x0], #16
  st1 {v16.8h}, [x0]

  ret lr
endfunc

// sub residuals
function transform_sub_8_4_10_neon, export=1

  mov x15, x0
  mov w10, #0x03FF
  dup v30.8h, w10
  dup v31.8h, wzr

  ld1 {v13.8h}, [x15], x1
  ld1 {v14.8h}, [x15], x1
  ld1 {v15.8h}, [x15], x1
  ld1 {v16.8h}, [x15]

  ld1 {v17.8h}, [x2], x3
  ld1 {v18.8h}, [x2], x3
  ld1 {v19.8h}, [x2], x3
  ld1 {v20.8h}, [x2]

  sub v13.8h, v13.8h, v17.8h
  sub v14.8h, v14.8h, v18.8h
  sub v15.8h, v15.8h, v19.8h
  sub v16.8h, v16.8h, v20.8h

  CLIP v13.8h, v31.8h, v30.8h
  CLIP v14.8h, v31.8h, v30.8h
  CLIP v15.8h, v31.8h, v30.8h
  CLIP v16.8h, v31.8h, v30.8h

  st1 {v13.8h}, [x0], x1
  st1 {v14.8h}, [x0], x1
  st1 {v15.8h}, [x0], x1
  st1 {v16.8h}, [x0]

  ret lr
endfunc

function transform_sub_16_2_10_neon, export=1

  mov x15, x0
  mov x16, x2
  mov w10, #0x03FF
  dup v30.8h, w10
  dup v31.8h, wzr

  ld1 {v13.8h}, [x15], x1
  ld1 {v14.8h}, [x15]
  add x15, x0, #16
  ld1 {v15.8h}, [x15], x1
  ld1 {v16.8h}, [x15]

  ld1 {v17.8h}, [x16], x3
  ld1 {v18.8h}, [x16]
  add x16, x2, #16
  ld1 {v19.8h}, [x16], x3
  ld1 {v20.8h}, [x16]

  sub v13.8h, v13.8h, v17.8h
  sub v14.8h, v14.8h, v18.8h
  sub v15.8h, v15.8h, v19.8h
  sub v16.8h, v16.8h, v20.8h

  CLIP v13.8h, v31.8h, v30.8h
  CLIP v14.8h, v31.8h, v30.8h
  CLIP v15.8h, v31.8h, v30.8h
  CLIP v16.8h, v31.8h, v30.8h

  mov x15, x0
  st1 {v13.8h}, [x15], x1
  st1 {v14.8h}, [x15]
  add x15, x0, #16
  st1 {v15.8h}, [x15], x1
  st1 {v16.8h}, [x15]

  ret lr
endfunc

function transform_sub_32_1_10_neon, export=1

  mov x15, x0
  mov w10, #0x03FF
  dup v30.8h, w10
  dup v31.8h, wzr

  ld1 {v13.8h}, [x15], #16
  ld1 {v14.8h}, [x15], #16
  ld1 {v15.8h}, [x15], #16
  ld1 {v16.8h}, [x15]

  ld1 {v17.8h}, [x2], #16
  ld1 {v18.8h}, [x2], #16
  ld1 {v19.8h}, [x2], #16
  ld1 {v20.8h}, [x2]

  sub v13.8h, v13.8h, v17.8h
  sub v14.8h, v14.8h, v18.8h
  sub v15.8h, v15.8h, v19.8h
  sub v16.8h, v16.8h, v20.8h

  CLIP v13.8h, v31.8h, v30.8h
  CLIP v14.8h, v31.8h, v30.8h
  CLIP v15.8h, v31.8h, v30.8h
  CLIP v16.8h, v31.8h, v30.8h

  st1 {v13.8h}, [x0], #16
  st1 {v14.8h}, [x0], #16
  st1 {v15.8h}, [x0], #16
  st1 {v16.8h}, [x0]

  ret lr
endfunc

// sub half
function transform_sub_half_8_4_10_neon, export=1

  mov x15, x0
  mov w10, #0x03FF
  dup v30.8h, w10
  dup v31.8h, wzr

  ld1 {v13.8h}, [x15], x1
  ld1 {v14.8h}, [x15], x1
  ld1 {v15.8h}, [x15], x1
  ld1 {v16.8h}, [x15]

  ld1 {v17.8h}, [x2], x3
  ld1 {v18.8h}, [x2], x3
  ld1 {v19.8h}, [x2], x3
  ld1 {v20.8h}, [x2]

  sub v17.8h, v31.8h, v17.8h
  sub v18.8h, v31.8h, v18.8h
  sub v19.8h, v31.8h, v19.8h
  sub v20.8h, v31.8h, v20.8h

  sshr v17.8h, v17.8h, #1
  sshr v18.8h, v18.8h, #1
  sshr v19.8h, v19.8h, #1
  sshr v20.8h, v20.8h, #1

  add v13.8h, v13.8h, v17.8h
  add v14.8h, v14.8h, v18.8h
  add v15.8h, v15.8h, v19.8h
  add v16.8h, v16.8h, v20.8h

  CLIP v13.8h, v31.8h, v30.8h
  CLIP v14.8h, v31.8h, v30.8h
  CLIP v15.8h, v31.8h, v30.8h
  CLIP v16.8h, v31.8h, v30.8h

  st1 {v13.8h}, [x0], x1
  st1 {v14.8h}, [x0], x1
  st1 {v15.8h}, [x0], x1
  st1 {v16.8h}, [x0]

  ret lr
endfunc

function transform_sub_half_16_2_10_neon, export=1

  mov x15, x0
  mov x16, x2
  mov w10, #0x03FF
  dup v30.8h, w10
  dup v31.8h, wzr

  ld1 {v13.8h}, [x15], x1
  ld1 {v14.8h}, [x15]
  add x15, x0, #16
  ld1 {v15.8h}, [x15], x1
  ld1 {v16.8h}, [x15]

  ld1 {v17.8h}, [x16], x3
  ld1 {v18.8h}, [x16]
  add x16, x2, #16
  ld1 {v19.8h}, [x16], x3
  ld1 {v20.8h}, [x16]

  sub v17.8h, v31.8h, v17.8h
  sub v18.8h, v31.8h, v18.8h
  sub v19.8h, v31.8h, v19.8h
  sub v20.8h, v31.8h, v20.8h

  sshr v17.8h, v17.8h, #1
  sshr v18.8h, v18.8h, #1
  sshr v19.8h, v19.8h, #1
  sshr v20.8h, v20.8h, #1

  add v13.8h, v13.8h, v17.8h
  add v14.8h, v14.8h, v18.8h
  add v15.8h, v15.8h, v19.8h
  add v16.8h, v16.8h, v20.8h

  CLIP v13.8h, v31.8h, v30.8h
  CLIP v14.8h, v31.8h, v30.8h
  CLIP v15.8h, v31.8h, v30.8h
  CLIP v16.8h, v31.8h, v30.8h

  mov x15, x0
  st1 {v13.8h}, [x15], x1
  st1 {v14.8h}, [x15]
  add x15, x0, #16
  st1 {v15.8h}, [x15], x1
  st1 {v16.8h}, [x15]

  ret lr
endfunc

function transform_sub_half_32_1_10_neon, export=1

  mov x15, x0
  mov w10, #0x03FF
  dup v30.8h, w10
  dup v31.8h, wzr

  ld1 {v13.8h}, [x15], #16
  ld1 {v14.8h}, [x15], #16
  ld1 {v15.8h}, [x15], #16
  ld1 {v16.8h}, [x15]

  ld1 {v17.8h}, [x2], #16
  ld1 {v18.8h}, [x2], #16
  ld1 {v19.8h}, [x2], #16
  ld1 {v20.8h}, [x2]

  sub v17.8h, v31.8h, v17.8h
  sub v18.8h, v31.8h, v18.8h
  sub v19.8h, v31.8h, v19.8h
  sub v20.8h, v31.8h, v20.8h

  sshr v17.8h, v17.8h, #1
  sshr v18.8h, v18.8h, #1
  sshr v19.8h, v19.8h, #1
  sshr v20.8h, v20.8h, #1

  add v13.8h, v13.8h, v17.8h
  add v14.8h, v14.8h, v18.8h
  add v15.8h, v15.8h, v19.8h
  add v16.8h, v16.8h, v20.8h

  CLIP v13.8h, v31.8h, v30.8h
  CLIP v14.8h, v31.8h, v30.8h
  CLIP v15.8h, v31.8h, v30.8h
  CLIP v16.8h, v31.8h, v30.8h

  st1 {v13.8h}, [x0], #16
  st1 {v14.8h}, [x0], #16
  st1 {v15.8h}, [x0], #16
  st1 {v16.8h}, [x0]

  ret lr
endfunc

// add scale
function transform_scale_add_8_4_10_neon, export=1

  mov x15, x0
  mov w10, #0x03FF
  mov w11, #1024
  dup v9.2d, x4
  dup v28.4s, w11
  dup v30.8h, w10
  dup v31.8h, wzr

  ld1 {v13.8h}, [x15], x1
  ld1 {v14.8h}, [x15], x1
  ld1 {v15.8h}, [x15], x1
  ld1 {v16.8h}, [x15]

  ld1 {v17.8h}, [x2], x3
  ld1 {v18.8h}, [x2], x3
  ld1 {v19.8h}, [x2], x3
  ld1 {v20.8h}, [x2]

  sshr v21.8h, v17.8h, #15
  sshr v22.8h, v18.8h, #15
  sshr v23.8h, v19.8h, #15
  sshr v24.8h, v20.8h, #15

  abs v17.8h, v17.8h
  abs v18.8h, v18.8h
  abs v19.8h, v19.8h
  abs v20.8h, v20.8h

  CLIP v17.8h, v31.8h, v30.8h
  CLIP v18.8h, v31.8h, v30.8h
  CLIP v19.8h, v31.8h, v30.8h
  CLIP v20.8h, v31.8h, v30.8h

  smull v1.4s, v17.4h, v9.h[0]
  smull v2.4s, v18.4h, v9.h[0]
  smull v3.4s, v19.4h, v9.h[0]
  smull v4.4s, v20.4h, v9.h[0]

  smull2 v5.4s, v17.8h, v9.h[0]
  smull2 v6.4s, v18.8h, v9.h[0]
  smull2 v7.4s, v19.8h, v9.h[0]
  smull2 v8.4s, v20.8h, v9.h[0]

  add v1.4s, v1.4s, v28.4s
  add v2.4s, v2.4s, v28.4s
  add v3.4s, v3.4s, v28.4s
  add v4.4s, v4.4s, v28.4s
  add v5.4s, v5.4s, v28.4s
  add v6.4s, v6.4s, v28.4s
  add v7.4s, v7.4s, v28.4s
  add v8.4s, v8.4s, v28.4s

  ushr v1.4s, v1.4s, #11
  ushr v2.4s, v2.4s, #11
  ushr v3.4s, v3.4s, #11
  ushr v4.4s, v4.4s, #11
  ushr v5.4s, v5.4s, #11
  ushr v6.4s, v6.4s, #11
  ushr v7.4s, v7.4s, #11
  ushr v8.4s, v8.4s, #11

  uzp1 v17.8h, v1.8h, v5.8h
  uzp1 v18.8h, v2.8h, v6.8h
  uzp1 v19.8h, v3.8h, v7.8h
  uzp1 v20.8h, v4.8h, v8.8h

  eor v17.16b, v17.16b, v21.16b
  eor v18.16b, v18.16b, v22.16b
  eor v19.16b, v19.16b, v23.16b
  eor v20.16b, v20.16b, v24.16b

  sub v17.8h, v17.8h, v21.8h
  sub v18.8h, v18.8h, v22.8h
  sub v19.8h, v19.8h, v23.8h
  sub v20.8h, v20.8h, v24.8h

  add v13.8h, v13.8h, v17.8h
  add v14.8h, v14.8h, v18.8h
  add v15.8h, v15.8h, v19.8h
  add v16.8h, v16.8h, v20.8h

  CLIP v13.8h, v31.8h, v30.8h
  CLIP v14.8h, v31.8h, v30.8h
  CLIP v15.8h, v31.8h, v30.8h
  CLIP v16.8h, v31.8h, v30.8h

  st1 {v13.8h}, [x0], x1
  st1 {v14.8h}, [x0], x1
  st1 {v15.8h}, [x0], x1
  st1 {v16.8h}, [x0]

  ret lr
endfunc

function transform_scale_add_16_2_10_neon, export=1

  mov x15, x0
  mov x16, x2
  mov w10, #0x03FF
  mov w11, #1024

  dup v30.8h, w10
  dup v31.8h, wzr
  dup v9.2d, x4
  dup v28.4s, w11


  ld1 {v13.8h}, [x15], x1
  ld1 {v14.8h}, [x15]
  add x15, x0, #16
  ld1 {v15.8h}, [x15], x1
  ld1 {v16.8h}, [x15]

  ld1 {v17.8h}, [x16], x3
  ld1 {v18.8h}, [x16]
  add x16, x2, #16
  ld1 {v19.8h}, [x16], x3
  ld1 {v20.8h}, [x16]

  sshr v21.8h, v17.8h, #15
  sshr v22.8h, v18.8h, #15
  sshr v23.8h, v19.8h, #15
  sshr v24.8h, v20.8h, #15

  abs v17.8h, v17.8h
  abs v18.8h, v18.8h
  abs v19.8h, v19.8h
  abs v20.8h, v20.8h

  CLIP v17.8h, v31.8h, v30.8h
  CLIP v18.8h, v31.8h, v30.8h
  CLIP v19.8h, v31.8h, v30.8h
  CLIP v20.8h, v31.8h, v30.8h

  smull v1.4s, v17.4h, v9.h[0]
  smull v2.4s, v18.4h, v9.h[0]
  smull v3.4s, v19.4h, v9.h[0]
  smull v4.4s, v20.4h, v9.h[0]

  smull2 v5.4s, v17.8h, v9.h[0]
  smull2 v6.4s, v18.8h, v9.h[0]
  smull2 v7.4s, v19.8h, v9.h[0]
  smull2 v8.4s, v20.8h, v9.h[0]

  add v1.4s, v1.4s, v28.4s
  add v2.4s, v2.4s, v28.4s
  add v3.4s, v3.4s, v28.4s
  add v4.4s, v4.4s, v28.4s
  add v5.4s, v5.4s, v28.4s
  add v6.4s, v6.4s, v28.4s
  add v7.4s, v7.4s, v28.4s
  add v8.4s, v8.4s, v28.4s

  ushr v1.4s, v1.4s, #11
  ushr v2.4s, v2.4s, #11
  ushr v3.4s, v3.4s, #11
  ushr v4.4s, v4.4s, #11
  ushr v5.4s, v5.4s, #11
  ushr v6.4s, v6.4s, #11
  ushr v7.4s, v7.4s, #11
  ushr v8.4s, v8.4s, #11

  uzp1 v17.8h, v1.8h, v5.8h
  uzp1 v18.8h, v2.8h, v6.8h
  uzp1 v19.8h, v3.8h, v7.8h
  uzp1 v20.8h, v4.8h, v8.8h

  eor v17.16b, v17.16b, v21.16b
  eor v18.16b, v18.16b, v22.16b
  eor v19.16b, v19.16b, v23.16b
  eor v20.16b, v20.16b, v24.16b

  sub v17.8h, v17.8h, v21.8h
  sub v18.8h, v18.8h, v22.8h
  sub v19.8h, v19.8h, v23.8h
  sub v20.8h, v20.8h, v24.8h

  add v13.8h, v13.8h, v17.8h
  add v14.8h, v14.8h, v18.8h
  add v15.8h, v15.8h, v19.8h
  add v16.8h, v16.8h, v20.8h

  CLIP v13.8h, v31.8h, v30.8h
  CLIP v14.8h, v31.8h, v30.8h
  CLIP v15.8h, v31.8h, v30.8h
  CLIP v16.8h, v31.8h, v30.8h

  mov x15, x0
  st1 {v13.8h}, [x15], x1
  st1 {v14.8h}, [x15]
  add x15, x0, #16
  st1 {v15.8h}, [x15], x1
  st1 {v16.8h}, [x15]

  ret lr
endfunc

function transform_scale_add_32_1_10_neon, export=1

  mov x15, x0
  mov w10, #0x03FF
  mov w11, #1024
  dup v9.2d, x4
  dup v28.4s, w11
  dup v30.8h, w10
  dup v31.8h, wzr

  ld1 {v13.8h}, [x15], #16
  ld1 {v14.8h}, [x15], #16
  ld1 {v15.8h}, [x15], #16
  ld1 {v16.8h}, [x15]

  ld1 {v17.8h}, [x2], #16
  ld1 {v18.8h}, [x2], #16
  ld1 {v19.8h}, [x2], #16
  ld1 {v20.8h}, [x2]

  sshr v21.8h, v17.8h, #15
  sshr v22.8h, v18.8h, #15
  sshr v23.8h, v19.8h, #15
  sshr v24.8h, v20.8h, #15

  abs v17.8h, v17.8h
  abs v18.8h, v18.8h
  abs v19.8h, v19.8h
  abs v20.8h, v20.8h

  CLIP v17.8h, v31.8h, v30.8h
  CLIP v18.8h, v31.8h, v30.8h
  CLIP v19.8h, v31.8h, v30.8h
  CLIP v20.8h, v31.8h, v30.8h

  smull v1.4s, v17.4h, v9.h[0]
  smull v2.4s, v18.4h, v9.h[0]
  smull v3.4s, v19.4h, v9.h[0]
  smull v4.4s, v20.4h, v9.h[0]

  smull2 v5.4s, v17.8h, v9.h[0]
  smull2 v6.4s, v18.8h, v9.h[0]
  smull2 v7.4s, v19.8h, v9.h[0]
  smull2 v8.4s, v20.8h, v9.h[0]

  add v1.4s, v1.4s, v28.4s
  add v2.4s, v2.4s, v28.4s
  add v3.4s, v3.4s, v28.4s
  add v4.4s, v4.4s, v28.4s
  add v5.4s, v5.4s, v28.4s
  add v6.4s, v6.4s, v28.4s
  add v7.4s, v7.4s, v28.4s
  add v8.4s, v8.4s, v28.4s

  ushr v1.4s, v1.4s, #11
  ushr v2.4s, v2.4s, #11
  ushr v3.4s, v3.4s, #11
  ushr v4.4s, v4.4s, #11
  ushr v5.4s, v5.4s, #11
  ushr v6.4s, v6.4s, #11
  ushr v7.4s, v7.4s, #11
  ushr v8.4s, v8.4s, #11

  uzp1 v17.8h, v1.8h, v5.8h
  uzp1 v18.8h, v2.8h, v6.8h
  uzp1 v19.8h, v3.8h, v7.8h
  uzp1 v20.8h, v4.8h, v8.8h

  eor v17.16b, v17.16b, v21.16b
  eor v18.16b, v18.16b, v22.16b
  eor v19.16b, v19.16b, v23.16b
  eor v20.16b, v20.16b, v24.16b

  sub v17.8h, v17.8h, v21.8h
  sub v18.8h, v18.8h, v22.8h
  sub v19.8h, v19.8h, v23.8h
  sub v20.8h, v20.8h, v24.8h

  add v13.8h, v13.8h, v17.8h
  add v14.8h, v14.8h, v18.8h
  add v15.8h, v15.8h, v19.8h
  add v16.8h, v16.8h, v20.8h

  CLIP v13.8h, v31.8h, v30.8h
  CLIP v14.8h, v31.8h, v30.8h
  CLIP v15.8h, v31.8h, v30.8h
  CLIP v16.8h, v31.8h, v30.8h

  st1 {v13.8h}, [x0], #16
  st1 {v14.8h}, [x0], #16
  st1 {v15.8h}, [x0], #16
  st1 {v16.8h}, [x0]

  ret lr
endfunc

  // scale add half
function transform_scale_add_half_8_4_10_neon, export=1

  mov x15, x0
  mov w10, #0x03FF
  mov w11, #1024
  dup v9.2d, x4
  dup v28.4s, w11
  dup v30.8h, w10
  dup v31.8h, wzr

  ld1 {v13.8h}, [x15], x1
  ld1 {v14.8h}, [x15], x1
  ld1 {v15.8h}, [x15], x1
  ld1 {v16.8h}, [x15]

  ld1 {v17.8h}, [x2], x3
  ld1 {v18.8h}, [x2], x3
  ld1 {v19.8h}, [x2], x3
  ld1 {v20.8h}, [x2]

  sshr v21.8h, v17.8h, #15
  sshr v22.8h, v18.8h, #15
  sshr v23.8h, v19.8h, #15
  sshr v24.8h, v20.8h, #15

  sshr v17.8h, v17.8h, #1
  sshr v18.8h, v18.8h, #1
  sshr v19.8h, v19.8h, #1
  sshr v20.8h, v20.8h, #1

  abs v17.8h, v17.8h
  abs v18.8h, v18.8h
  abs v19.8h, v19.8h
  abs v20.8h, v20.8h

  CLIP v17.8h, v31.8h, v30.8h
  CLIP v18.8h, v31.8h, v30.8h
  CLIP v19.8h, v31.8h, v30.8h
  CLIP v20.8h, v31.8h, v30.8h

  smull v1.4s, v17.4h, v9.h[0]
  smull v2.4s, v18.4h, v9.h[0]
  smull v3.4s, v19.4h, v9.h[0]
  smull v4.4s, v20.4h, v9.h[0]

  smull2 v5.4s, v17.8h, v9.h[0]
  smull2 v6.4s, v18.8h, v9.h[0]
  smull2 v7.4s, v19.8h, v9.h[0]
  smull2 v8.4s, v20.8h, v9.h[0]

  add v1.4s, v1.4s, v28.4s
  add v2.4s, v2.4s, v28.4s
  add v3.4s, v3.4s, v28.4s
  add v4.4s, v4.4s, v28.4s
  add v5.4s, v5.4s, v28.4s
  add v6.4s, v6.4s, v28.4s
  add v7.4s, v7.4s, v28.4s
  add v8.4s, v8.4s, v28.4s

  ushr v1.4s, v1.4s, #11
  ushr v2.4s, v2.4s, #11
  ushr v3.4s, v3.4s, #11
  ushr v4.4s, v4.4s, #11
  ushr v5.4s, v5.4s, #11
  ushr v6.4s, v6.4s, #11
  ushr v7.4s, v7.4s, #11
  ushr v8.4s, v8.4s, #11

  uzp1 v17.8h, v1.8h, v5.8h
  uzp1 v18.8h, v2.8h, v6.8h
  uzp1 v19.8h, v3.8h, v7.8h
  uzp1 v20.8h, v4.8h, v8.8h

  eor v17.16b, v17.16b, v21.16b
  eor v18.16b, v18.16b, v22.16b
  eor v19.16b, v19.16b, v23.16b
  eor v20.16b, v20.16b, v24.16b

  sub v17.8h, v17.8h, v21.8h
  sub v18.8h, v18.8h, v22.8h
  sub v19.8h, v19.8h, v23.8h
  sub v20.8h, v20.8h, v24.8h

  add v13.8h, v13.8h, v17.8h
  add v14.8h, v14.8h, v18.8h
  add v15.8h, v15.8h, v19.8h
  add v16.8h, v16.8h, v20.8h

  CLIP v13.8h, v31.8h, v30.8h
  CLIP v14.8h, v31.8h, v30.8h
  CLIP v15.8h, v31.8h, v30.8h
  CLIP v16.8h, v31.8h, v30.8h

  st1 {v13.8h}, [x0], x1
  st1 {v14.8h}, [x0], x1
  st1 {v15.8h}, [x0], x1
  st1 {v16.8h}, [x0]

  ret lr
endfunc

function transform_scale_add_half_16_2_10_neon, export=1

  mov x15, x0
  mov x16, x2
  mov w10, #0x03FF
  mov w11, #1024

  dup v30.8h, w10
  dup v31.8h, wzr
  dup v9.2d, x4
  dup v28.4s, w11

  ld1 {v13.8h}, [x15], x1
  ld1 {v14.8h}, [x15]
  add x15, x0, #16
  ld1 {v15.8h}, [x15], x1
  ld1 {v16.8h}, [x15]

  ld1 {v17.8h}, [x16], x3
  ld1 {v18.8h}, [x16]
  add x16, x2, #16
  ld1 {v19.8h}, [x16], x3
  ld1 {v20.8h}, [x16]

  sshr v21.8h, v17.8h, #15
  sshr v22.8h, v18.8h, #15
  sshr v23.8h, v19.8h, #15
  sshr v24.8h, v20.8h, #15

  sshr v17.8h, v17.8h, #1
  sshr v18.8h, v18.8h, #1
  sshr v19.8h, v19.8h, #1
  sshr v20.8h, v20.8h, #1

  abs v17.8h, v17.8h
  abs v18.8h, v18.8h
  abs v19.8h, v19.8h
  abs v20.8h, v20.8h

  CLIP v17.8h, v31.8h, v30.8h
  CLIP v18.8h, v31.8h, v30.8h
  CLIP v19.8h, v31.8h, v30.8h
  CLIP v20.8h, v31.8h, v30.8h

  smull v1.4s, v17.4h, v9.h[0]
  smull v2.4s, v18.4h, v9.h[0]
  smull v3.4s, v19.4h, v9.h[0]
  smull v4.4s, v20.4h, v9.h[0]

  smull2 v5.4s, v17.8h, v9.h[0]
  smull2 v6.4s, v18.8h, v9.h[0]
  smull2 v7.4s, v19.8h, v9.h[0]
  smull2 v8.4s, v20.8h, v9.h[0]

  add v1.4s, v1.4s, v28.4s
  add v2.4s, v2.4s, v28.4s
  add v3.4s, v3.4s, v28.4s
  add v4.4s, v4.4s, v28.4s
  add v5.4s, v5.4s, v28.4s
  add v6.4s, v6.4s, v28.4s
  add v7.4s, v7.4s, v28.4s
  add v8.4s, v8.4s, v28.4s

  ushr v1.4s, v1.4s, #11
  ushr v2.4s, v2.4s, #11
  ushr v3.4s, v3.4s, #11
  ushr v4.4s, v4.4s, #11
  ushr v5.4s, v5.4s, #11
  ushr v6.4s, v6.4s, #11
  ushr v7.4s, v7.4s, #11
  ushr v8.4s, v8.4s, #11

  uzp1 v17.8h, v1.8h, v5.8h
  uzp1 v18.8h, v2.8h, v6.8h
  uzp1 v19.8h, v3.8h, v7.8h
  uzp1 v20.8h, v4.8h, v8.8h

  eor v17.16b, v17.16b, v21.16b
  eor v18.16b, v18.16b, v22.16b
  eor v19.16b, v19.16b, v23.16b
  eor v20.16b, v20.16b, v24.16b

  sub v17.8h, v17.8h, v21.8h
  sub v18.8h, v18.8h, v22.8h
  sub v19.8h, v19.8h, v23.8h
  sub v20.8h, v20.8h, v24.8h

  add v13.8h, v13.8h, v17.8h
  add v14.8h, v14.8h, v18.8h
  add v15.8h, v15.8h, v19.8h
  add v16.8h, v16.8h, v20.8h

  CLIP v13.8h, v31.8h, v30.8h
  CLIP v14.8h, v31.8h, v30.8h
  CLIP v15.8h, v31.8h, v30.8h
  CLIP v16.8h, v31.8h, v30.8h

  mov x15, x0
  st1 {v13.8h}, [x15], x1
  st1 {v14.8h}, [x15]
  add x15, x0, #16
  st1 {v15.8h}, [x15], x1
  st1 {v16.8h}, [x15]

  ret lr
endfunc

function transform_scale_add_half_32_1_10_neon, export=1

  mov x15, x0
  mov w10, #0x03FF
  mov w11, #1024
  dup v9.2d, x4
  dup v28.4s, w11
  dup v30.8h, w10
  dup v31.8h, wzr

  ld1 {v13.8h}, [x15], #16
  ld1 {v14.8h}, [x15], #16
  ld1 {v15.8h}, [x15], #16
  ld1 {v16.8h}, [x15]

  ld1 {v17.8h}, [x2], #16
  ld1 {v18.8h}, [x2], #16
  ld1 {v19.8h}, [x2], #16
  ld1 {v20.8h}, [x2]

  sshr v21.8h, v17.8h, #15
  sshr v22.8h, v18.8h, #15
  sshr v23.8h, v19.8h, #15
  sshr v24.8h, v20.8h, #15

  sshr v17.8h, v17.8h, #1
  sshr v18.8h, v18.8h, #1
  sshr v19.8h, v19.8h, #1
  sshr v20.8h, v20.8h, #1

  abs v17.8h, v17.8h
  abs v18.8h, v18.8h
  abs v19.8h, v19.8h
  abs v20.8h, v20.8h

  CLIP v17.8h, v31.8h, v30.8h
  CLIP v18.8h, v31.8h, v30.8h
  CLIP v19.8h, v31.8h, v30.8h
  CLIP v20.8h, v31.8h, v30.8h

  smull v1.4s, v17.4h, v9.h[0]
  smull v2.4s, v18.4h, v9.h[0]
  smull v3.4s, v19.4h, v9.h[0]
  smull v4.4s, v20.4h, v9.h[0]

  smull2 v5.4s, v17.8h, v9.h[0]
  smull2 v6.4s, v18.8h, v9.h[0]
  smull2 v7.4s, v19.8h, v9.h[0]
  smull2 v8.4s, v20.8h, v9.h[0]

  add v1.4s, v1.4s, v28.4s
  add v2.4s, v2.4s, v28.4s
  add v3.4s, v3.4s, v28.4s
  add v4.4s, v4.4s, v28.4s
  add v5.4s, v5.4s, v28.4s
  add v6.4s, v6.4s, v28.4s
  add v7.4s, v7.4s, v28.4s
  add v8.4s, v8.4s, v28.4s

  ushr v1.4s, v1.4s, #11
  ushr v2.4s, v2.4s, #11
  ushr v3.4s, v3.4s, #11
  ushr v4.4s, v4.4s, #11
  ushr v5.4s, v5.4s, #11
  ushr v6.4s, v6.4s, #11
  ushr v7.4s, v7.4s, #11
  ushr v8.4s, v8.4s, #11

  uzp1 v17.8h, v1.8h, v5.8h
  uzp1 v18.8h, v2.8h, v6.8h
  uzp1 v19.8h, v3.8h, v7.8h
  uzp1 v20.8h, v4.8h, v8.8h

  eor v17.16b, v17.16b, v21.16b
  eor v18.16b, v18.16b, v22.16b
  eor v19.16b, v19.16b, v23.16b
  eor v20.16b, v20.16b, v24.16b

  sub v17.8h, v17.8h, v21.8h
  sub v18.8h, v18.8h, v22.8h
  sub v19.8h, v19.8h, v23.8h
  sub v20.8h, v20.8h, v24.8h

  add v13.8h, v13.8h, v17.8h
  add v14.8h, v14.8h, v18.8h
  add v15.8h, v15.8h, v19.8h
  add v16.8h, v16.8h, v20.8h

  CLIP v13.8h, v31.8h, v30.8h
  CLIP v14.8h, v31.8h, v30.8h
  CLIP v15.8h, v31.8h, v30.8h
  CLIP v16.8h, v31.8h, v30.8h

  st1 {v13.8h}, [x0], #16
  st1 {v14.8h}, [x0], #16
  st1 {v15.8h}, [x0], #16
  st1 {v16.8h}, [x0]

  ret lr
endfunc

// sub scale
function transform_scale_sub_8_4_10_neon, export=1

  mov x15, x0
  mov w10, #0x03FF
  mov w11, #1024
  mov w12, -1
  dup v30.8h, w10
  dup v31.8h, wzr
  dup v9.2d, x4
  dup v10.8h, w12
  dup v28.4s, w11

  ld1 {v13.8h}, [x15], x1
  ld1 {v14.8h}, [x15], x1
  ld1 {v15.8h}, [x15], x1
  ld1 {v16.8h}, [x15]

  ld1 {v17.8h}, [x2], x3
  ld1 {v18.8h}, [x2], x3
  ld1 {v19.8h}, [x2], x3
  ld1 {v20.8h}, [x2]

  sshr v21.8h, v17.8h, #15
  sshr v22.8h, v18.8h, #15
  sshr v23.8h, v19.8h, #15
  sshr v24.8h, v20.8h, #15

  abs v17.8h, v17.8h
  abs v18.8h, v18.8h
  abs v19.8h, v19.8h
  abs v20.8h, v20.8h

  CLIP v17.8h, v31.8h, v30.8h
  CLIP v18.8h, v31.8h, v30.8h
  CLIP v19.8h, v31.8h, v30.8h
  CLIP v20.8h, v31.8h, v30.8h

  smull v1.4s, v17.4h, v9.h[0]
  smull v2.4s, v18.4h, v9.h[0]
  smull v3.4s, v19.4h, v9.h[0]
  smull v4.4s, v20.4h, v9.h[0]

  smull2 v5.4s, v17.8h, v9.h[0]
  smull2 v6.4s, v18.8h, v9.h[0]
  smull2 v7.4s, v19.8h, v9.h[0]
  smull2 v8.4s, v20.8h, v9.h[0]

  add v1.4s, v1.4s, v28.4s
  add v2.4s, v2.4s, v28.4s
  add v3.4s, v3.4s, v28.4s
  add v4.4s, v4.4s, v28.4s
  add v5.4s, v5.4s, v28.4s
  add v6.4s, v6.4s, v28.4s
  add v7.4s, v7.4s, v28.4s
  add v8.4s, v8.4s, v28.4s

  ushr v1.4s, v1.4s, #11
  ushr v2.4s, v2.4s, #11
  ushr v3.4s, v3.4s, #11
  ushr v4.4s, v4.4s, #11
  ushr v5.4s, v5.4s, #11
  ushr v6.4s, v6.4s, #11
  ushr v7.4s, v7.4s, #11
  ushr v8.4s, v8.4s, #11

  uzp1 v17.8h, v1.8h, v5.8h
  uzp1 v18.8h, v2.8h, v6.8h
  uzp1 v19.8h, v3.8h, v7.8h
  uzp1 v20.8h, v4.8h, v8.8h

  eor v21.16b, v21.16b, v10.16b
  eor v22.16b, v22.16b, v10.16b
  eor v23.16b, v23.16b, v10.16b
  eor v24.16b, v24.16b, v10.16b

  eor v17.16b, v17.16b, v21.16b
  eor v18.16b, v18.16b, v22.16b
  eor v19.16b, v19.16b, v23.16b
  eor v20.16b, v20.16b, v24.16b

  sub v17.8h, v17.8h, v21.8h
  sub v18.8h, v18.8h, v22.8h
  sub v19.8h, v19.8h, v23.8h
  sub v20.8h, v20.8h, v24.8h

  add v13.8h, v13.8h, v17.8h
  add v14.8h, v14.8h, v18.8h
  add v15.8h, v15.8h, v19.8h
  add v16.8h, v16.8h, v20.8h

  CLIP v13.8h, v31.8h, v30.8h
  CLIP v14.8h, v31.8h, v30.8h
  CLIP v15.8h, v31.8h, v30.8h
  CLIP v16.8h, v31.8h, v30.8h

  st1 {v13.8h}, [x0], x1
  st1 {v14.8h}, [x0], x1
  st1 {v15.8h}, [x0], x1
  st1 {v16.8h}, [x0]

  ret lr
endfunc

function transform_scale_sub_16_2_10_neon, export=1

  mov x15, x0
  mov x16, x2
  mov w10, #0x03FF
  mov w11, #1024
  mov w12, -1
  dup v30.8h, w10
  dup v31.8h, wzr
  dup v9.2d, x4
  dup v10.8h, w12
  dup v28.4s, w11

  ld1 {v13.8h}, [x15], x1
  ld1 {v14.8h}, [x15]
  add x15, x0, #16
  ld1 {v15.8h}, [x15], x1
  ld1 {v16.8h}, [x15]

  ld1 {v17.8h}, [x16], x3
  ld1 {v18.8h}, [x16]
  add x16, x2, #16
  ld1 {v19.8h}, [x16], x3
  ld1 {v20.8h}, [x16]

  sshr v21.8h, v17.8h, #15
  sshr v22.8h, v18.8h, #15
  sshr v23.8h, v19.8h, #15
  sshr v24.8h, v20.8h, #15

  abs v17.8h, v17.8h
  abs v18.8h, v18.8h
  abs v19.8h, v19.8h
  abs v20.8h, v20.8h

  CLIP v17.8h, v31.8h, v30.8h
  CLIP v18.8h, v31.8h, v30.8h
  CLIP v19.8h, v31.8h, v30.8h
  CLIP v20.8h, v31.8h, v30.8h

  smull v1.4s, v17.4h, v9.h[0]
  smull v2.4s, v18.4h, v9.h[0]
  smull v3.4s, v19.4h, v9.h[0]
  smull v4.4s, v20.4h, v9.h[0]

  smull2 v5.4s, v17.8h, v9.h[0]
  smull2 v6.4s, v18.8h, v9.h[0]
  smull2 v7.4s, v19.8h, v9.h[0]
  smull2 v8.4s, v20.8h, v9.h[0]

  add v1.4s, v1.4s, v28.4s
  add v2.4s, v2.4s, v28.4s
  add v3.4s, v3.4s, v28.4s
  add v4.4s, v4.4s, v28.4s
  add v5.4s, v5.4s, v28.4s
  add v6.4s, v6.4s, v28.4s
  add v7.4s, v7.4s, v28.4s
  add v8.4s, v8.4s, v28.4s

  ushr v1.4s, v1.4s, #11
  ushr v2.4s, v2.4s, #11
  ushr v3.4s, v3.4s, #11
  ushr v4.4s, v4.4s, #11
  ushr v5.4s, v5.4s, #11
  ushr v6.4s, v6.4s, #11
  ushr v7.4s, v7.4s, #11
  ushr v8.4s, v8.4s, #11

  uzp1 v17.8h, v1.8h, v5.8h
  uzp1 v18.8h, v2.8h, v6.8h
  uzp1 v19.8h, v3.8h, v7.8h
  uzp1 v20.8h, v4.8h, v8.8h

  eor v21.16b, v21.16b, v10.16b
  eor v22.16b, v22.16b, v10.16b
  eor v23.16b, v23.16b, v10.16b
  eor v24.16b, v24.16b, v10.16b

  eor v17.16b, v17.16b, v21.16b
  eor v18.16b, v18.16b, v22.16b
  eor v19.16b, v19.16b, v23.16b
  eor v20.16b, v20.16b, v24.16b

  sub v17.8h, v17.8h, v21.8h
  sub v18.8h, v18.8h, v22.8h
  sub v19.8h, v19.8h, v23.8h
  sub v20.8h, v20.8h, v24.8h

  add v13.8h, v13.8h, v17.8h
  add v14.8h, v14.8h, v18.8h
  add v15.8h, v15.8h, v19.8h
  add v16.8h, v16.8h, v20.8h

  CLIP v13.8h, v31.8h, v30.8h
  CLIP v14.8h, v31.8h, v30.8h
  CLIP v15.8h, v31.8h, v30.8h
  CLIP v16.8h, v31.8h, v30.8h

  mov x15, x0
  st1 {v13.8h}, [x15], x1
  st1 {v14.8h}, [x15]
  add x15, x0, #16
  st1 {v15.8h}, [x15], x1
  st1 {v16.8h}, [x15]

  ret lr
endfunc

function transform_scale_sub_32_1_10_neon, export=1

  mov x15, x0
  mov w10, #0x03FF
  mov w11, #1024
  mov w12, -1
  dup v30.8h, w10
  dup v31.8h, wzr
  dup v9.2d, x4
  dup v10.8h, w12
  dup v28.4s, w11

  ld1 {v13.8h}, [x15], #16
  ld1 {v14.8h}, [x15], #16
  ld1 {v15.8h}, [x15], #16
  ld1 {v16.8h}, [x15]

  ld1 {v17.8h}, [x2], #16
  ld1 {v18.8h}, [x2], #16
  ld1 {v19.8h}, [x2], #16
  ld1 {v20.8h}, [x2]

  sshr v21.8h, v17.8h, #15
  sshr v22.8h, v18.8h, #15
  sshr v23.8h, v19.8h, #15
  sshr v24.8h, v20.8h, #15

  abs v17.8h, v17.8h
  abs v18.8h, v18.8h
  abs v19.8h, v19.8h
  abs v20.8h, v20.8h

  CLIP v17.8h, v31.8h, v30.8h
  CLIP v18.8h, v31.8h, v30.8h
  CLIP v19.8h, v31.8h, v30.8h
  CLIP v20.8h, v31.8h, v30.8h

  smull v1.4s, v17.4h, v9.h[0]
  smull v2.4s, v18.4h, v9.h[0]
  smull v3.4s, v19.4h, v9.h[0]
  smull v4.4s, v20.4h, v9.h[0]

  smull2 v5.4s, v17.8h, v9.h[0]
  smull2 v6.4s, v18.8h, v9.h[0]
  smull2 v7.4s, v19.8h, v9.h[0]
  smull2 v8.4s, v20.8h, v9.h[0]

  add v1.4s, v1.4s, v28.4s
  add v2.4s, v2.4s, v28.4s
  add v3.4s, v3.4s, v28.4s
  add v4.4s, v4.4s, v28.4s
  add v5.4s, v5.4s, v28.4s
  add v6.4s, v6.4s, v28.4s
  add v7.4s, v7.4s, v28.4s
  add v8.4s, v8.4s, v28.4s

  ushr v1.4s, v1.4s, #11
  ushr v2.4s, v2.4s, #11
  ushr v3.4s, v3.4s, #11
  ushr v4.4s, v4.4s, #11
  ushr v5.4s, v5.4s, #11
  ushr v6.4s, v6.4s, #11
  ushr v7.4s, v7.4s, #11
  ushr v8.4s, v8.4s, #11

  uzp1 v17.8h, v1.8h, v5.8h
  uzp1 v18.8h, v2.8h, v6.8h
  uzp1 v19.8h, v3.8h, v7.8h
  uzp1 v20.8h, v4.8h, v8.8h

  eor v21.16b, v21.16b, v10.16b
  eor v22.16b, v22.16b, v10.16b
  eor v23.16b, v23.16b, v10.16b
  eor v24.16b, v24.16b, v10.16b

  eor v17.16b, v17.16b, v21.16b
  eor v18.16b, v18.16b, v22.16b
  eor v19.16b, v19.16b, v23.16b
  eor v20.16b, v20.16b, v24.16b

  sub v17.8h, v17.8h, v21.8h
  sub v18.8h, v18.8h, v22.8h
  sub v19.8h, v19.8h, v23.8h
  sub v20.8h, v20.8h, v24.8h

  add v13.8h, v13.8h, v17.8h
  add v14.8h, v14.8h, v18.8h
  add v15.8h, v15.8h, v19.8h
  add v16.8h, v16.8h, v20.8h

  CLIP v13.8h, v31.8h, v30.8h
  CLIP v14.8h, v31.8h, v30.8h
  CLIP v15.8h, v31.8h, v30.8h
  CLIP v16.8h, v31.8h, v30.8h

  st1 {v13.8h}, [x0], #16
  st1 {v14.8h}, [x0], #16
  st1 {v15.8h}, [x0], #16
  st1 {v16.8h}, [x0]

  ret lr
endfunc

// scale sub half
function transform_scale_sub_half_8_4_10_neon, export=1
  mov x15, x0
  mov w10, #0x03FF
  mov w11, #1024
  mov w12, -1
  dup v30.8h, w10
  dup v31.8h, wzr
  dup v9.2d, x4
  dup v10.8h, w12
  dup v28.4s, w11

  ld1 {v13.8h}, [x15], x1
  ld1 {v14.8h}, [x15], x1
  ld1 {v15.8h}, [x15], x1
  ld1 {v16.8h}, [x15]

  ld1 {v17.8h}, [x2], x3
  ld1 {v18.8h}, [x2], x3
  ld1 {v19.8h}, [x2], x3
  ld1 {v20.8h}, [x2]

  sshr v21.8h, v17.8h, #15
  sshr v22.8h, v18.8h, #15
  sshr v23.8h, v19.8h, #15
  sshr v24.8h, v20.8h, #15

  sub v17.8h, v31.8h, v17.8h
  sub v18.8h, v31.8h, v18.8h
  sub v19.8h, v31.8h, v19.8h
  sub v20.8h, v31.8h, v20.8h

  sshr v17.8h, v17.8h, #1
  sshr v18.8h, v18.8h, #1
  sshr v19.8h, v19.8h, #1
  sshr v20.8h, v20.8h, #1

  abs v17.8h, v17.8h
  abs v18.8h, v18.8h
  abs v19.8h, v19.8h
  abs v20.8h, v20.8h

  CLIP v17.8h, v31.8h, v30.8h
  CLIP v18.8h, v31.8h, v30.8h
  CLIP v19.8h, v31.8h, v30.8h
  CLIP v20.8h, v31.8h, v30.8h

  smull v1.4s, v17.4h, v9.h[0]
  smull v2.4s, v18.4h, v9.h[0]
  smull v3.4s, v19.4h, v9.h[0]
  smull v4.4s, v20.4h, v9.h[0]

  smull2 v5.4s, v17.8h, v9.h[0]
  smull2 v6.4s, v18.8h, v9.h[0]
  smull2 v7.4s, v19.8h, v9.h[0]
  smull2 v8.4s, v20.8h, v9.h[0]

  add v1.4s, v1.4s, v28.4s
  add v2.4s, v2.4s, v28.4s
  add v3.4s, v3.4s, v28.4s
  add v4.4s, v4.4s, v28.4s
  add v5.4s, v5.4s, v28.4s
  add v6.4s, v6.4s, v28.4s
  add v7.4s, v7.4s, v28.4s
  add v8.4s, v8.4s, v28.4s

  ushr v1.4s, v1.4s, #11
  ushr v2.4s, v2.4s, #11
  ushr v3.4s, v3.4s, #11
  ushr v4.4s, v4.4s, #11
  ushr v5.4s, v5.4s, #11
  ushr v6.4s, v6.4s, #11
  ushr v7.4s, v7.4s, #11
  ushr v8.4s, v8.4s, #11

  uzp1 v17.8h, v1.8h, v5.8h
  uzp1 v18.8h, v2.8h, v6.8h
  uzp1 v19.8h, v3.8h, v7.8h
  uzp1 v20.8h, v4.8h, v8.8h

  eor v21.16b, v21.16b, v10.16b
  eor v22.16b, v22.16b, v10.16b
  eor v23.16b, v23.16b, v10.16b
  eor v24.16b, v24.16b, v10.16b

  eor v17.16b, v17.16b, v21.16b
  eor v18.16b, v18.16b, v22.16b
  eor v19.16b, v19.16b, v23.16b
  eor v20.16b, v20.16b, v24.16b

  sub v17.8h, v17.8h, v21.8h
  sub v18.8h, v18.8h, v22.8h
  sub v19.8h, v19.8h, v23.8h
  sub v20.8h, v20.8h, v24.8h

  add v13.8h, v13.8h, v17.8h
  add v14.8h, v14.8h, v18.8h
  add v15.8h, v15.8h, v19.8h
  add v16.8h, v16.8h, v20.8h

  CLIP v13.8h, v31.8h, v30.8h
  CLIP v14.8h, v31.8h, v30.8h
  CLIP v15.8h, v31.8h, v30.8h
  CLIP v16.8h, v31.8h, v30.8h

  st1 {v13.8h}, [x0], x1
  st1 {v14.8h}, [x0], x1
  st1 {v15.8h}, [x0], x1
  st1 {v16.8h}, [x0]

  ret lr
endfunc

function transform_scale_sub_half_16_2_10_neon, export=1

  mov x15, x0
  mov x16, x2
  mov w10, #0x03FF
  mov w11, #1024
  mov w12, -1
  dup v30.8h, w10
  dup v31.8h, wzr
  dup v9.2d, x4
  dup v10.8h, w12
  dup v28.4s, w11

  ld1 {v13.8h}, [x15], x1
  ld1 {v14.8h}, [x15]
  add x15, x0, #16
  ld1 {v15.8h}, [x15], x1
  ld1 {v16.8h}, [x15]

  ld1 {v17.8h}, [x16], x3
  ld1 {v18.8h}, [x16]
  add x16, x2, #16
  ld1 {v19.8h}, [x16], x3
  ld1 {v20.8h}, [x16]

  sshr v21.8h, v17.8h, #15
  sshr v22.8h, v18.8h, #15
  sshr v23.8h, v19.8h, #15
  sshr v24.8h, v20.8h, #15

  sub v17.8h, v31.8h, v17.8h
  sub v18.8h, v31.8h, v18.8h
  sub v19.8h, v31.8h, v19.8h
  sub v20.8h, v31.8h, v20.8h

  sshr v17.8h, v17.8h, #1
  sshr v18.8h, v18.8h, #1
  sshr v19.8h, v19.8h, #1
  sshr v20.8h, v20.8h, #1

  abs v17.8h, v17.8h
  abs v18.8h, v18.8h
  abs v19.8h, v19.8h
  abs v20.8h, v20.8h

  CLIP v17.8h, v31.8h, v30.8h
  CLIP v18.8h, v31.8h, v30.8h
  CLIP v19.8h, v31.8h, v30.8h
  CLIP v20.8h, v31.8h, v30.8h

  smull v1.4s, v17.4h, v9.h[0]
  smull v2.4s, v18.4h, v9.h[0]
  smull v3.4s, v19.4h, v9.h[0]
  smull v4.4s, v20.4h, v9.h[0]

  smull2 v5.4s, v17.8h, v9.h[0]
  smull2 v6.4s, v18.8h, v9.h[0]
  smull2 v7.4s, v19.8h, v9.h[0]
  smull2 v8.4s, v20.8h, v9.h[0]

  add v1.4s, v1.4s, v28.4s
  add v2.4s, v2.4s, v28.4s
  add v3.4s, v3.4s, v28.4s
  add v4.4s, v4.4s, v28.4s
  add v5.4s, v5.4s, v28.4s
  add v6.4s, v6.4s, v28.4s
  add v7.4s, v7.4s, v28.4s
  add v8.4s, v8.4s, v28.4s

  ushr v1.4s, v1.4s, #11
  ushr v2.4s, v2.4s, #11
  ushr v3.4s, v3.4s, #11
  ushr v4.4s, v4.4s, #11
  ushr v5.4s, v5.4s, #11
  ushr v6.4s, v6.4s, #11
  ushr v7.4s, v7.4s, #11
  ushr v8.4s, v8.4s, #11

  uzp1 v17.8h, v1.8h, v5.8h
  uzp1 v18.8h, v2.8h, v6.8h
  uzp1 v19.8h, v3.8h, v7.8h
  uzp1 v20.8h, v4.8h, v8.8h

  eor v21.16b, v21.16b, v10.16b
  eor v22.16b, v22.16b, v10.16b
  eor v23.16b, v23.16b, v10.16b
  eor v24.16b, v24.16b, v10.16b

  eor v17.16b, v17.16b, v21.16b
  eor v18.16b, v18.16b, v22.16b
  eor v19.16b, v19.16b, v23.16b
  eor v20.16b, v20.16b, v24.16b

  sub v17.8h, v17.8h, v21.8h
  sub v18.8h, v18.8h, v22.8h
  sub v19.8h, v19.8h, v23.8h
  sub v20.8h, v20.8h, v24.8h

  add v13.8h, v13.8h, v17.8h
  add v14.8h, v14.8h, v18.8h
  add v15.8h, v15.8h, v19.8h
  add v16.8h, v16.8h, v20.8h

  CLIP v13.8h, v31.8h, v30.8h
  CLIP v14.8h, v31.8h, v30.8h
  CLIP v15.8h, v31.8h, v30.8h
  CLIP v16.8h, v31.8h, v30.8h

  mov x15, x0
  st1 {v13.8h}, [x15], x1
  st1 {v14.8h}, [x15]
  add x15, x0, #16
  st1 {v15.8h}, [x15], x1
  st1 {v16.8h}, [x15]

  ret lr
endfunc

function transform_scale_sub_half_32_1_10_neon, export=1

  mov x15, x0
  mov w10, #0x03FF
  mov w11, #1024
  mov w12, -1
  dup v30.8h, w10
  dup v31.8h, wzr
  dup v9.2d, x4
  dup v10.8h, w12
  dup v28.4s, w11

  ld1 {v13.8h}, [x15], #16
  ld1 {v14.8h}, [x15], #16
  ld1 {v15.8h}, [x15], #16
  ld1 {v16.8h}, [x15]

  ld1 {v17.8h}, [x2], #16
  ld1 {v18.8h}, [x2], #16
  ld1 {v19.8h}, [x2], #16
  ld1 {v20.8h}, [x2]

  sshr v21.8h, v17.8h, #15
  sshr v22.8h, v18.8h, #15
  sshr v23.8h, v19.8h, #15
  sshr v24.8h, v20.8h, #15

  sub v17.8h, v31.8h, v17.8h
  sub v18.8h, v31.8h, v18.8h
  sub v19.8h, v31.8h, v19.8h
  sub v20.8h, v31.8h, v20.8h

  sshr v17.8h, v17.8h, #1
  sshr v18.8h, v18.8h, #1
  sshr v19.8h, v19.8h, #1
  sshr v20.8h, v20.8h, #1

  abs v17.8h, v17.8h
  abs v18.8h, v18.8h
  abs v19.8h, v19.8h
  abs v20.8h, v20.8h

  CLIP v17.8h, v31.8h, v30.8h
  CLIP v18.8h, v31.8h, v30.8h
  CLIP v19.8h, v31.8h, v30.8h
  CLIP v20.8h, v31.8h, v30.8h

  smull v1.4s, v17.4h, v9.h[0]
  smull v2.4s, v18.4h, v9.h[0]
  smull v3.4s, v19.4h, v9.h[0]
  smull v4.4s, v20.4h, v9.h[0]

  smull2 v5.4s, v17.8h, v9.h[0]
  smull2 v6.4s, v18.8h, v9.h[0]
  smull2 v7.4s, v19.8h, v9.h[0]
  smull2 v8.4s, v20.8h, v9.h[0]

  add v1.4s, v1.4s, v28.4s
  add v2.4s, v2.4s, v28.4s
  add v3.4s, v3.4s, v28.4s
  add v4.4s, v4.4s, v28.4s
  add v5.4s, v5.4s, v28.4s
  add v6.4s, v6.4s, v28.4s
  add v7.4s, v7.4s, v28.4s
  add v8.4s, v8.4s, v28.4s

  ushr v1.4s, v1.4s, #11
  ushr v2.4s, v2.4s, #11
  ushr v3.4s, v3.4s, #11
  ushr v4.4s, v4.4s, #11
  ushr v5.4s, v5.4s, #11
  ushr v6.4s, v6.4s, #11
  ushr v7.4s, v7.4s, #11
  ushr v8.4s, v8.4s, #11

  uzp1 v17.8h, v1.8h, v5.8h
  uzp1 v18.8h, v2.8h, v6.8h
  uzp1 v19.8h, v3.8h, v7.8h
  uzp1 v20.8h, v4.8h, v8.8h

  eor v21.16b, v21.16b, v10.16b
  eor v22.16b, v22.16b, v10.16b
  eor v23.16b, v23.16b, v10.16b
  eor v24.16b, v24.16b, v10.16b

  eor v17.16b, v17.16b, v21.16b
  eor v18.16b, v18.16b, v22.16b
  eor v19.16b, v19.16b, v23.16b
  eor v20.16b, v20.16b, v24.16b

  sub v17.8h, v17.8h, v21.8h
  sub v18.8h, v18.8h, v22.8h
  sub v19.8h, v19.8h, v23.8h
  sub v20.8h, v20.8h, v24.8h

  add v13.8h, v13.8h, v17.8h
  add v14.8h, v14.8h, v18.8h
  add v15.8h, v15.8h, v19.8h
  add v16.8h, v16.8h, v20.8h

  CLIP v13.8h, v31.8h, v30.8h
  CLIP v14.8h, v31.8h, v30.8h
  CLIP v15.8h, v31.8h, v30.8h
  CLIP v16.8h, v31.8h, v30.8h

  st1 {v13.8h}, [x0], #16
  st1 {v14.8h}, [x0], #16
  st1 {v15.8h}, [x0], #16
  st1 {v16.8h}, [x0]

  ret lr
endfunc