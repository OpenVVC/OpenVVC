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

  ushr v17.8h, v17.8h, #1
  ushr v18.8h, v18.8h, #1
  ushr v19.8h, v19.8h, #1
  ushr v20.8h, v20.8h, #1

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

   ushr v17.8h, v17.8h, #1
   ushr v18.8h, v18.8h, #1
   ushr v19.8h, v19.8h, #1
   ushr v20.8h, v20.8h, #1

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

  ushr v17.8h, v17.8h, #1
  ushr v18.8h, v18.8h, #1
  ushr v19.8h, v19.8h, #1
  ushr v20.8h, v20.8h, #1

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