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
  dup v31.8h, wzr   // store all '0'

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
  add v16.8h, v16.8h, v20.8h    // adding L and R vectors

  CLIP v13.8h, v31.8h, v30.8h
  CLIP v14.8h, v31.8h, v30.8h
  CLIP v15.8h, v31.8h, v30.8h
  CLIP v16.8h, v31.8h, v30.8h

  st1 {v13.8h}, [x0], x1
  st1 {v14.8h}, [x0], x1
  st1 {v15.8h}, [x0], x1
  st1 {v16.8h}, [x0] // store result to the dst ptr

  ret lr
endfunc
