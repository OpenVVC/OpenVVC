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
 
#include "arm/asm.S"

#define MAX_PB_SIZE #64
#define MAX_PB_SIZE_8 #128
#define MAX_PB_SIZE_16 #256
#define SHIFT_VAL #4 //14 - bitdepth, for 10 -> 4, for 8 -> 6

.macro CLIP, val, min_val, max_val
      SMAX  \val, \val, \min_val
      SMIN  \val, \val, \max_val
.endm

function put_vvc_bi0_pel_pixels_10_4_neon, export=1
      mov    x12, MAX_PB_SIZE_16
      lsl    x2, x2, #1
0:
      subs x3, x3, #2
      LD1   {V0.4H}, [x1], x2
      LD1   {V1.4H}, [x1], x2
      SHL V0.4H, V0.4H, SHIFT_VAL
      SHL V1.4H, V1.4H, SHIFT_VAL
      ST1   {V0.4H}, [x0], x12
      ST1   {V1.4H}, [x0], x12
      b.ne 0b
      ret lr
endfunc

function put_vvc_bi0_pel_pixels_10_8_neon, export=1
      mov    x12, MAX_PB_SIZE_16
      lsl    x2, x2, #1
0:
      subs x3, x3, #2
      LD1   {V0.8H}, [x1], x2
      LD1   {V1.8H}, [x1], x2
      SHL V0.8H, V0.8H, SHIFT_VAL
      SHL V1.8H, V1.8H, SHIFT_VAL
      ST1   {V0.8H}, [x0], x12
      ST1   {V1.8H}, [x0], x12
      b.ne 0b
      ret lr
endfunc

function put_vvc_bi0_pel_pixels_10_16_neon, export=1
      lsl    x2, x2, #1
      mov    x12, MAX_PB_SIZE_16
0:
      subs x3, x3, #2
      LD1   {V0.8H-V1.8H}, [x1], x2
      LD1   {V2.8H-V3.8H}, [x1], x2
      SHL V0.8H, V0.8H, SHIFT_VAL
      SHL V1.8H, V1.8H, SHIFT_VAL
      SHL V2.8H, V2.8H, SHIFT_VAL
      SHL V3.8H, V3.8H, SHIFT_VAL
      ST1   {V0.8H-V1.8H}, [x0], x12
      ST1   {V2.8H-V3.8H}, [x0], x12
      b.ne 0b
      ret lr
endfunc

function put_vvc_bi0_pel_pixels_10_32_neon, export=1
      lsl    x2, x2, #1
      mov    x12, MAX_PB_SIZE_16
0:
      subs x3, x3, #2
      LD1   {V0.8H-V3.8H}, [x1], x2
      LD1   {V4.8H-V7.8H}, [x1], x2
      SHL V0.8H, V0.8H, SHIFT_VAL
      SHL V1.8H, V1.8H, SHIFT_VAL
      SHL V2.8H, V2.8H, SHIFT_VAL
      SHL V3.8H, V3.8H, SHIFT_VAL
      SHL V4.8H, V4.8H, SHIFT_VAL
      SHL V5.8H, V5.8H, SHIFT_VAL
      SHL V6.8H, V6.8H, SHIFT_VAL
      SHL V7.8H, V7.8H, SHIFT_VAL
      ST1   {V0.8H-V3.8H}, [x0], x12
      ST1   {V4.8H-V7.8H}, [x0], x12
      b.ne 0b
      ret lr
endfunc

function put_vvc_bi1_pel_pixels_10_4_neon, export=1
      mov    x12, MAX_PB_SIZE_16 //src2stride
      mov    w11, #0x1 //offset val
      mov    w10, #0x03FF //max val
      lsl    x1, x1, #1 //_dststride
      lsl    x3, x3, #1 //_srcstride
      DUP    V0.4H, w11 //offset vector
      DUP    V1.4H, wzr //min val vector
      DUP    V2.4H, w10 //max val vector
0:
      subs x5, x5, #2
      LD1   {V3.4H}, [x2], x3
      LD1   {V4.4H}, [x2], x3
      LD1   {V7.4H}, [x4], x12
      LD1   {V8.4H}, [x4], x12
      ADD   V3.4H,   V3.4H, V0.4H
      ADD   V4.4H,   V4.4H, V0.4H
      SHL   V3.4H, V3.4H, SHIFT_VAL
      SHL   V4.4H, V4.4H, SHIFT_VAL
      SHADD   V3.4H, V3.4H, V7.4H
      SHADD   V4.4H, V4.4H, V8.4H
      SSHR  V3.4H, V3.4H, SHIFT_VAL
      SSHR  V4.4H, V4.4H, SHIFT_VAL
      CLIP  V3.4H, V1.4H, V2.4H
      CLIP  V4.4H, V1.4H, V2.4H
      ST1   {V3.4H}, [x0], x1
      ST1   {V4.4H}, [x0], x1
      b.ne 0b
      ret lr
endfunc

function put_vvc_bi1_pel_pixels_10_8_neon, export=1
      mov    x12, MAX_PB_SIZE_16 //src2stride
      mov    w11, #0x1 //offset val
      mov    w10, #0x03FF //max val
      lsl    x1, x1, #1 //_dststride
      lsl    x3, x3, #1 //_srcstride
      DUP    V0.8H, w11 //offset vector
      DUP    V1.8H, wzr //min val vector
      DUP    V2.8H, w10 //max val vector
0:
      subs x5, x5, #2
      LD1   {V3.8H}, [x2], x3
      LD1   {V4.8H}, [x2], x3
      LD1   {V7.8H}, [x4], x12
      LD1   {V8.8H}, [x4], x12
      ADD   V3.8H,   V3.8H, V0.8H
      ADD   V4.8H,   V4.8H, V0.8H
      SHL   V3.8H, V3.8H, SHIFT_VAL
      SHL   V4.8H, V4.8H, SHIFT_VAL
      SHADD   V3.8H, V3.8H, V7.8H
      SHADD   V4.8H, V4.8H, V8.8H
      SSHR  V3.8H, V3.8H, SHIFT_VAL
      SSHR  V4.8H, V4.8H, SHIFT_VAL
      CLIP  V3.8H, V1.8H, V2.8H
      CLIP  V4.8H, V1.8H, V2.8H
      ST1   {V3.8H}, [x0], x1
      ST1   {V4.8H}, [x0], x1
      b.ne 0b
      ret lr
endfunc

function put_vvc_bi1_pel_pixels_10_16_neon, export=1
      mov    x12, MAX_PB_SIZE_16 //src2stride
      mov    w11, #0x1 //offset val
      mov    w10, #0x03FF //max val
      lsl    x1, x1, #1 //_dststride
      lsl    x3, x3, #1 //_srcstride
      DUP    V0.8H, w11 //offset vector
      DUP    V1.8H, wzr //min val vector
      DUP    V2.8H, w10 //max val vector
0:
      subs x5, x5, #1
      LD1   {V3.8H-V4.8H}, [x2], x3
      LD1   {V7.8H-V8.8H}, [x4], x12
      ADD   V3.8H,   V3.8H, V0.8H
      ADD   V4.8H,   V4.8H, V0.8H
      SHL   V3.8H, V3.8H, SHIFT_VAL
      SHL   V4.8H, V4.8H, SHIFT_VAL
      SHADD   V3.8H, V3.8H, V7.8H
      SHADD   V4.8H, V4.8H, V8.8H
      SSHR  V3.8H, V3.8H, SHIFT_VAL
      SSHR  V4.8H, V4.8H, SHIFT_VAL
      CLIP  V3.8H, V1.8H, V2.8H
      CLIP  V4.8H, V1.8H, V2.8H
      ST1   {V3.8H-V4.8H}, [x0], x1
      b.ne 0b
      ret lr
endfunc

function put_vvc_bi1_pel_pixels_10_32_neon, export=1
      mov    x12, MAX_PB_SIZE_16 //src2stride
      mov    w11, #0x1 //offset val
      mov    w10, #0x03FF //max val
      lsl    x1, x1, #1 //_dststride
      lsl    x3, x3, #1 //_srcstride
      DUP    V0.8H, w11 //offset vector
      DUP    V1.8H, wzr //min val vector
      DUP    V2.8H, w10 //max val vector
0:
      subs x5, x5, #1
      LD1   {V3.8H-V6.8H}, [x2], x3
      LD1   {V7.8H-V10.8H}, [x4], x12
      ADD   V3.8H,   V3.8H, V0.8H
      ADD   V4.8H,   V4.8H, V0.8H
      ADD   V5.8H,   V5.8H, V0.8H
      ADD   V6.8H,   V6.8H, V0.8H
      SHL   V3.8H, V3.8H, SHIFT_VAL
      SHL   V4.8H, V4.8H, SHIFT_VAL
      SHL   V5.8H, V5.8H, SHIFT_VAL
      SHL   V6.8H, V6.8H, SHIFT_VAL
      SHADD   V3.8H, V3.8H, V7.8H
      SHADD   V4.8H, V4.8H, V8.8H
      SHADD   V5.8H, V5.8H, V9.8H
      SHADD   V6.8H, V6.8H, V10.8H
      SSHR  V3.8H, V3.8H, SHIFT_VAL
      SSHR  V4.8H, V4.8H, SHIFT_VAL
      SSHR  V5.8H, V5.8H, SHIFT_VAL
      SSHR  V6.8H, V6.8H, SHIFT_VAL
      CLIP  V3.8H, V1.8H, V2.8H
      CLIP  V4.8H, V1.8H, V2.8H
      CLIP  V5.8H, V1.8H, V2.8H
      CLIP  V6.8H, V1.8H, V2.8H
      ST1   {V3.8H-V6.8H}, [x0], x1
      b.ne 0b
      ret lr
endfunc

.macro uni_qpel_filter, out, coef=V0, l0=V1, l1=V2, l2=V3, l3=V4, l4=V5, l5=V6, l6=V7, l7=V8, shift=#6, s1=.4H, s2=.8H, s3=.H
      SMULL2 V11.4S,  \l0\s2, \coef\s3[0]
      SMLAL2 V11.4S,  \l1\s2, \coef\s3[1]
      SMLAL2 V11.4S, \l2\s2, \coef\s3[2]
      SMLAL2 V11.4S, \l3\s2, \coef\s3[3]
      SMLAL2 V11.4S, \l4\s2, \coef\s3[4]
      SMLAL2 V11.4S, \l5\s2, \coef\s3[5]
      SMLAL2 V11.4S, \l6\s2, \coef\s3[6]
      SMLAL2 V11.4S, \l7\s2, \coef\s3[7]
      SMULL V10.4S, \l0\s1, \coef\s3[0]
      SMLAL V10.4S, \l1\s1, \coef\s3[1]
      SMLAL V10.4S, \l2\s1, \coef\s3[2]
      SMLAL V10.4S, \l3\s1, \coef\s3[3]
      SMLAL V10.4S, \l4\s1, \coef\s3[4]
      SMLAL V10.4S, \l5\s1, \coef\s3[5]
      SMLAL V10.4S, \l6\s1, \coef\s3[6]
      SMLAL V10.4S, \l7\s1, \coef\s3[7]
      ADD   V10.4S, V10.4S, V29.4S
      ADD   V11.4S, V11.4S, V29.4S
      SSHR  V10.4S, V10.4S, \shift
      SSHR  V11.4S, V11.4S, \shift
      UZP1  \out, V10.8H, V11.8H
      CLIP  \out, V30.8H, V31.8H
.endm

.macro bi0_qpel_filter, out, coef=V0, l0=V1, l1=V2, l2=V3, l3=V4, l4=V5, l5=V6, l6=V7, l7=V8, shift=#2,s1=.4H, s2=.8H, s3=.H
      SMULL2 V11.4S,  \l0\s2, \coef\s3[0]
      SMLAL2 V11.4S,  \l1\s2, \coef\s3[1]
      SMLAL2 V11.4S, \l2\s2, \coef\s3[2]
      SMLAL2 V11.4S, \l3\s2, \coef\s3[3]
      SMLAL2 V11.4S, \l4\s2, \coef\s3[4]
      SMLAL2 V11.4S, \l5\s2, \coef\s3[5]
      SMLAL2 V11.4S, \l6\s2, \coef\s3[6]
      SMLAL2 V11.4S, \l7\s2, \coef\s3[7]
      SMULL V10.4S, \l0\s1, \coef\s3[0]
      SMLAL V10.4S, \l1\s1, \coef\s3[1]
      SMLAL V10.4S, \l2\s1, \coef\s3[2]
      SMLAL V10.4S, \l3\s1, \coef\s3[3]
      SMLAL V10.4S, \l4\s1, \coef\s3[4]
      SMLAL V10.4S, \l5\s1, \coef\s3[5]
      SMLAL V10.4S, \l6\s1, \coef\s3[6]
      SMLAL V10.4S, \l7\s1, \coef\s3[7]
      SSHR  V10.4S, V10.4S, \shift
      SSHR  V11.4S, V11.4S, \shift
      UZP1  \out, V10.8H, V11.8H
.endm

.macro bi1_qpel_filter, out, src1, coef=V0, l0=V1, l1=V2, l2=V3, l3=V4, l4=V5, l5=V6, l6=V7, l7=V8, shift=#2, s1=.4H, s2=.8H, s3=.H
      SMULL2 V11.4S,  \l0\s2, \coef\s3[0]
      SMLAL2 V11.4S,  \l1\s2, \coef\s3[1]
      SMLAL2 V11.4S, \l2\s2, \coef\s3[2]
      SMLAL2 V11.4S, \l3\s2, \coef\s3[3]
      SMLAL2 V11.4S, \l4\s2, \coef\s3[4]
      SMLAL2 V11.4S, \l5\s2, \coef\s3[5]
      SMLAL2 V11.4S, \l6\s2, \coef\s3[6]
      SMLAL2 V11.4S, \l7\s2, \coef\s3[7]
      SMULL V10.4S, \l0\s1, \coef\s3[0]
      SMLAL V10.4S, \l1\s1, \coef\s3[1]
      SMLAL V10.4S, \l2\s1, \coef\s3[2]
      SMLAL V10.4S, \l3\s1, \coef\s3[3]
      SMLAL V10.4S, \l4\s1, \coef\s3[4]
      SMLAL V10.4S, \l5\s1, \coef\s3[5]
      SMLAL V10.4S, \l6\s1, \coef\s3[6]
      SMLAL V10.4S, \l7\s1, \coef\s3[7]
      ADD   V10.4S, V10.4S, V29.4S
      ADD   V11.4S, V11.4S, V29.4S
      SSHR  V10.4S, V10.4S, \shift
      SSHR  V11.4S, V11.4S, \shift
      SADDW V10.4S, V10.4S,  \src1\s1
      SADDW2 V11.4S, V11.4S, \src1\s2
      SSHR  V10.4S, V10.4S, #5
      SSHR  V11.4S, V11.4S, #5
      UZP1  \out, V10.8H, V11.8H
      CLIP  \out, V30.8H, V31.8H
.endm

.macro uni_qpel_h_filter, out, last_jump
      LD1   {V1.8H}, [x2], x12
      LD1   {V2.8H}, [x2], x12
      LD1   {V3.8H}, [x2], x12
      LD1   {V4.8H}, [x2], x12
      LD1   {V5.8H}, [x2], x12
      LD1   {V6.8H}, [x2], x12
      LD1   {V7.8H}, [x2], x12
      LD1   {V8.8H}, [x2], \last_jump
      uni_qpel_filter \out, V0
.endm

.macro bi0_qpel_h_filter, out, last_jump, coef=V0, jump=x12, in=x1, l0=V1, l1=V2, l2=V3, l3=V4, l4=V5, l5=V6, l6=V7, l7=V8, s1=.4H, s2=.8H, s3=.H
      LD1   {\l0\s2}, [\in], \jump
      LD1   {\l1\s2}, [\in], \jump
      LD1   {\l2\s2}, [\in], \jump
      LD1   {\l3\s2}, [\in], \jump
      LD1   {\l4\s2}, [\in], \jump
      LD1   {\l5\s2}, [\in], \jump
      LD1   {\l6\s2}, [\in], \jump
      LD1   {\l7\s2}, [\in], \last_jump
      bi0_qpel_filter \out, \coef, \l0, \l1, \l2, \l3, \l4, \l5, \l6, \l7
.endm

function put_vvc_uni_qpel_h_10_8_neon, export=1
      mov    x12, #2
      mov    w11, #32 //offset val
      mov    w10, #0x03FF //max val
      sub    x2, x2, #6
      sub    x3, x3, #7
      lsl    x1, x1, #1 //_dststride
      lsl    x3, x3, #1 //_srcstride
      DUP    V29.4S, w11 //offset vector
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      LD1    {V0.8B}, [x5]
      SADDW   V0.8H, V30.8H, V0.8B
0:
      subs x4, x4, #2
      uni_qpel_h_filter V17.8H, x3
      uni_qpel_h_filter V18.8H, x3
      ST1   {V17.8H}, [x0], x1
      ST1   {V18.8H}, [x0], x1
      b.ne 0b
      ret lr
endfunc

function put_vvc_uni_qpel_h_10_16_neon, export=1
      mov    x12, #2
      mov    w11, #32 //offset val
      mov    w10, #0x03FF //max val
      sub    x2, x2, #6
      sub    x3, x3, #15
      lsl    x1, x1, #1 //_dststride
      lsl    x3, x3, #1 //_srcstride
      DUP    V29.4S, w11 //offset vector
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      LD1    {V0.8B}, [x5]
      SADDW   V0.8H, V30.8H, V0.8B
0:
      subs x4, x4, #1
      uni_qpel_h_filter V17.8H, x12
      uni_qpel_h_filter V18.8H, x3
      ST1   {V17.8H-V18.8H}, [x0], x1
      b.ne 0b
      ret lr
endfunc

function put_vvc_uni_qpel_h_10_32_neon, export=1
      mov    x12, #2
      mov    w11, #32 //offset val
      mov    w10, #0x03FF //max val
      sub    x2, x2, #6
      sub    x3, x3, #31
      lsl    x1, x1, #1 //_dststride
      lsl    x3, x3, #1 //_srcstride
      DUP    V29.4S, w11 //offset vector
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      LD1    {V0.8B}, [x5]
      SADDW   V0.8H, V30.8H, V0.8B
0:
      subs x4, x4, #1
      uni_qpel_h_filter V17.8H, x12
      uni_qpel_h_filter V18.8H, x12
      uni_qpel_h_filter V19.8H, x12
      uni_qpel_h_filter V20.8H, x3
      ST1   {V17.8H-V20.8H}, [x0], x1
      b.ne 0b
      ret lr
endfunc

.macro uni_qpel_v_filter_8lines, out, in
      LD1    {V1.8H}, [\in], x3
      LD1    {V2.8H}, [\in], x3
      LD1    {V3.8H}, [\in], x3
      LD1    {V4.8H}, [\in], x3
      LD1    {V5.8H}, [\in], x3
      LD1    {V6.8H}, [\in], x3
      LD1    {V7.8H}, [\in], x3
0:              subs x4, x4, #8
      LD1    {V8.8H}, [\in], x3
      uni_qpel_filter V17.8H, V0, V1, V2, V3, V4, V5, V6, V7, V8
      LD1    {V1.8H}, [\in], x3
      uni_qpel_filter V18.8H, V0, V2, V3, V4, V5, V6, V7, V8, V1
      LD1    {V2.8H}, [\in], x3
      uni_qpel_filter V19.8H, V0, V3, V4, V5, V6, V7, V8, V1, V2
      LD1    {V3.8H}, [\in], x3
      uni_qpel_filter V20.8H, V0, V4, V5, V6, V7, V8, V1, V2, V3
      LD1    {V4.8H}, [\in], x3
      uni_qpel_filter V21.8H, V0, V5, V6, V7, V8, V1, V2, V3, V4
      LD1    {V5.8H}, [\in], x3
      uni_qpel_filter V22.8H, V0, V6, V7, V8, V1, V2, V3, V4, V5
      LD1    {V6.8H}, [\in], x3
      uni_qpel_filter V23.8H, V0, V7, V8, V1, V2, V3, V4, V5, V6
      LD1    {V7.8H}, [\in], x3
      uni_qpel_filter V24.8H, V0, V8, V1, V2, V3, V4, V5, V6, V7
      ST1   {V17.8H}, [\out], x1
      ST1   {V18.8H}, [\out], x1
      ST1   {V19.8H}, [\out], x1
      ST1   {V20.8H}, [\out], x1
      ST1   {V21.8H}, [\out], x1
      ST1   {V22.8H}, [\out], x1
      ST1   {V23.8H}, [\out], x1
      ST1   {V24.8H}, [\out], x1
      b.ne 0b
.endm

.macro uni_qpel_v_filter_4lines, out, in
      LD1    {V1.8H}, [\in], x3
      LD1    {V2.8H}, [\in], x3
      LD1    {V3.8H}, [\in], x3
      LD1    {V4.8H}, [\in], x3
      LD1    {V5.8H}, [\in], x3
      LD1    {V6.8H}, [\in], x3
      LD1    {V7.8H}, [\in], x3
0:              subs x4, x4, #4
      LD1    {V8.8H}, [\in], x3
      uni_qpel_filter V17.8H, V0, V1, V2, V3, V4, V5, V6, V7, V8
      LD1    {V1.8H}, [\in], x3
      uni_qpel_filter V18.8H, V0, V2, V3, V4, V5, V6, V7, V8, V1
      LD1    {V2.8H}, [\in], x3
      uni_qpel_filter V19.8H, V0, V3, V4, V5, V6, V7, V8, V1, V2
      LD1    {V3.8H}, [\in], x3
      uni_qpel_filter V20.8H, V0, V4, V5, V6, V7, V8, V1, V2, V3
      ST1   {V17.8H}, [\out], x1
      ST1   {V18.8H}, [\out], x1
      ST1   {V19.8H}, [\out], x1
      ST1   {V20.8H}, [\out], x1
      b.ne 0b
.endm

function put_vvc_uni_qpel_v_10_8_neon_8lines, export=1
      mov    x12, #2
      mov    w11, #32 //offset val
      mov    w10, #0x03FF //max val
      lsl    x1, x1, #1 //_dststride
      lsl    x3, x3, #1 //_srcstride
      lsl    x13, x3, #1 //2*_srcstride
      sub    x2, x2, x3
      sub    x2, x2, x13
      DUP    V29.4S, w11 //offset vector
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      LD1    {V0.8B}, [x5]
      SADDW   V0.8H, V30.8H, V0.8B
      uni_qpel_v_filter_8lines x0, x2
      ret lr
endfunc

function put_vvc_uni_qpel_v_10_8_neon_4lines, export=1
      mov    x12, #2
      mov    w11, #32 //offset val
      mov    w10, #0x03FF //max val
      lsl    x1, x1, #1 //_dststride
      lsl    x3, x3, #1 //_srcstride
      lsl    x13, x3, #1 //2*_srcstride
      sub    x2, x2, x3
      sub    x2, x2, x13
      DUP    V29.4S, w11 //offset vector
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      LD1    {V0.8B}, [x5]
      SADDW   V0.8H, V30.8H, V0.8B
      uni_qpel_v_filter_4lines x0, x2
      ret lr
endfunc

function put_vvc_uni_qpel_v_10_16_neon_8lines, export=1
      mov    x12, #2
      mov    w11, #32 //offset val
      mov    w10, #0x03FF //max val
      lsl    x1, x1, #1 //_dststride
      lsl    x3, x3, #1 //_srcstride
      lsl    x13, x3, #1
      sub    x2, x2, x3
      sub    x2, x2, x13
      mov    x9, x2  //save x2[-3*_srcstride]
      mov    x8, x0  //save x0[0]
      mov    x14, x4  //save height
      DUP    V29.4S, w11 //offset vector
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      LD1    {V0.8B}, [x5]
      SADDW   V0.8H, V30.8H, V0.8B
      uni_qpel_v_filter_8lines x0, x2
      add x0, x8, #16
      add x2, x9, #16
      mov    x4, x14
      uni_qpel_v_filter_8lines x0, x2
      ret lr
endfunc

function put_vvc_uni_qpel_v_10_16_neon_4lines, export=1
      mov    x12, #2
      mov    w11, #32 //offset val
      mov    w10, #0x03FF //max val
      lsl    x1, x1, #1 //_dststride
      lsl    x3, x3, #1 //_srcstride
      lsl    x13, x3, #1
      sub    x2, x2, x3
      sub    x2, x2, x13
      mov    x9, x2  //save x2[-3*_srcstride]
      mov    x8, x0  //save x0[0]
      mov    x14, x4  //save height
      DUP    V29.4S, w11 //offset vector
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      LD1    {V0.8B}, [x5]
      SADDW   V0.8H, V30.8H, V0.8B
      uni_qpel_v_filter_4lines x0, x2
      add x0, x8, #16
      add x2, x9, #16
      mov    x4, x14
      uni_qpel_v_filter_4lines x0, x2
      ret lr
endfunc

function put_vvc_uni_qpel_v_10_32_neon_8lines, export=1
      mov    x12, #2
      mov    w11, #32 //offset val
      mov    w10, #0x03FF //max val
      lsl    x1, x1, #1 //_dststride
      lsl    x3, x3, #1 //_srcstride
      lsl    x13, x3, #1
      sub    x2, x2, x3
      sub    x2, x2, x13
      mov    x9, x2  //save x2[-3*_srcstride]
      mov    x8, x0  //save x0[0]
      mov    x14, x4  //save height
      DUP    V29.4S, w11 //offset vector
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      LD1    {V0.8B}, [x5]
      SADDW   V0.8H, V30.8H, V0.8B
      uni_qpel_v_filter_8lines x0, x2
      add x0, x8, #16
      add x2, x9, #16
      mov    x4, x14
      uni_qpel_v_filter_8lines x0, x2
      add x0, x8, #32
      add x2, x9, #32
      mov    x4, x14
      uni_qpel_v_filter_8lines x0, x2
      add x0, x8, #48
      add x2, x9, #48
      mov    x4, x14
      uni_qpel_v_filter_8lines x0, x2
      ret lr
endfunc

function put_vvc_uni_qpel_v_10_32_neon_4lines, export=1
      mov    x12, #2
      mov    w11, #32 //offset val
      mov    w10, #0x03FF //max val
      lsl    x1, x1, #1 //_dststride
      lsl    x3, x3, #1 //_srcstride
      lsl    x13, x3, #1
      sub    x2, x2, x3
      sub    x2, x2, x13
      mov    x9, x2  //save x2[-3*_srcstride]
      mov    x8, x0  //save x0[0]
      mov    x14, x4  //save height
      DUP    V29.4S, w11 //offset vector
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      LD1    {V0.8B}, [x5]
      SADDW   V0.8H, V30.8H, V0.8B
      uni_qpel_v_filter_4lines x0, x2
      add x0, x8, #16
      add x2, x9, #16
      mov    x4, x14
      uni_qpel_v_filter_4lines x0, x2
      add x0, x8, #32
      add x2, x9, #32
      mov    x4, x14
      uni_qpel_v_filter_4lines x0, x2
      add x0, x8, #48
      add x2, x9, #48
      mov    x4, x14
      uni_qpel_v_filter_4lines x0, x2
      ret lr
endfunc

.macro uni_qpel_hv_filter_8lines, out, in, coefx=V0, coefy=V1, last_jump=x3, jump=x12
      bi0_qpel_h_filter V2.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
      bi0_qpel_h_filter V3.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
      bi0_qpel_h_filter V4.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
      bi0_qpel_h_filter V5.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
      bi0_qpel_h_filter V6.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
      bi0_qpel_h_filter V7.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
      bi0_qpel_h_filter V8.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
0:              subs x4, x4, #8
      bi0_qpel_h_filter V9.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
      uni_qpel_filter V20.8H, \coefy, V2, V3, V4, V5, V6, V7, V8, V9, shift=#10
      bi0_qpel_h_filter V2.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
      uni_qpel_filter V21.8H, \coefy, V3, V4, V5, V6, V7, V8, V9, V2, shift=#10
      bi0_qpel_h_filter V3.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
      uni_qpel_filter V22.8H, \coefy, V4, V5, V6, V7, V8, V9, V2, V3, shift=#10
      bi0_qpel_h_filter V4.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
      uni_qpel_filter V23.8H, \coefy, V5, V6, V7, V8, V9, V2, V3, V4, shift=#10
      bi0_qpel_h_filter V5.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
      uni_qpel_filter V24.8H, \coefy, V6, V7, V8, V9, V2, V3, V4, V5, shift=#10
      bi0_qpel_h_filter V6.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
      uni_qpel_filter V25.8H, \coefy, V7, V8, V9, V2, V3, V4, V5, V6, shift=#10
      bi0_qpel_h_filter V7.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
      uni_qpel_filter V26.8H, \coefy, V8, V9, V2, V3, V4, V5, V6, V7, shift=#10
      bi0_qpel_h_filter V8.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
      uni_qpel_filter V27.8H, \coefy, V9, V2, V3, V4, V5, V6, V7, V8, shift=#10
      ST1   {V20.8H}, [\out], x1
      ST1   {V21.8H}, [\out], x1
      ST1   {V22.8H}, [\out], x1
      ST1   {V23.8H}, [\out], x1
      ST1   {V24.8H}, [\out], x1
      ST1   {V25.8H}, [\out], x1
      ST1   {V26.8H}, [\out], x1
      ST1   {V27.8H}, [\out], x1
      b.ne 0b
.endm

.macro uni_qpel_hv_filter_4lines, out, in, coefx=V0, coefy=V1, last_jump=x3, jump=x12
      bi0_qpel_h_filter V2.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
      bi0_qpel_h_filter V3.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
      bi0_qpel_h_filter V4.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
      bi0_qpel_h_filter V5.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
      bi0_qpel_h_filter V6.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
      bi0_qpel_h_filter V7.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
      bi0_qpel_h_filter V8.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
0:              subs x4, x4, #4
      bi0_qpel_h_filter V9.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
      uni_qpel_filter V20.8H, \coefy, V2, V3, V4, V5, V6, V7, V8, V9, shift=#10
      bi0_qpel_h_filter V2.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
      uni_qpel_filter V21.8H, \coefy, V3, V4, V5, V6, V7, V8, V9, V2, shift=#10
      bi0_qpel_h_filter V3.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
      uni_qpel_filter V22.8H, \coefy, V4, V5, V6, V7, V8, V9, V2, V3, shift=#10
      bi0_qpel_h_filter V4.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
      uni_qpel_filter V23.8H, \coefy, V5, V6, V7, V8, V9, V2, V3, V4, shift=#10
      ST1   {V20.8H}, [\out], x1
      ST1   {V21.8H}, [\out], x1
      ST1   {V22.8H}, [\out], x1
      ST1   {V23.8H}, [\out], x1
      b.ne 0b
.endm

function put_vvc_uni_qpel_hv_10_8_neon_8lines, export=1
      mov    w11, #512 //offset val
      mov    w10, #0x03FF //max val
      DUP    V29.4S, w11 //offset vector
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      mov    x12, #2
      sub    x2, x2, #6//start h src[]
      sub    x13, x3, #7
      lsl    x13, x13, #1 // last jump h
      lsl    x1, x1, #1 //_dststride
      lsl    x14, x3, #1 //_srcstride
      lsl    x15, x14, #1 //2*_srcstride
      sub    x2, x2, x14
      sub    x2, x2, x15//start hv src[]
      LD1    {V0.8B}, [x5]
      LD1    {V1.8B}, [x6]
      SADDW   V0.8H, V30.8H, V0.8B
      SADDW   V1.8H, V30.8H, V1.8B
      uni_qpel_hv_filter_8lines x0, x2, V0, V1, x13, x12
      ret lr
endfunc

function put_vvc_uni_qpel_hv_10_8_neon_4lines, export=1
      mov    w11, #512 //offset val
      mov    w10, #0x03FF //max val
      DUP    V29.4S, w11 //offset vector
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      mov    x12, #2
      sub    x2, x2, #6//start h src[]
      sub    x13, x3, #7
      lsl    x13, x13, #1 // last jump h
      lsl    x1, x1, #1 //_dststride
      lsl    x14, x3, #1 //_srcstride
      lsl    x15, x14, #1 //2*_srcstride
      sub    x2, x2, x14
      sub    x2, x2, x15//start hv src[]
      LD1    {V0.8B}, [x5]
      LD1    {V1.8B}, [x6]
      SADDW   V0.8H, V30.8H, V0.8B
      SADDW   V1.8H, V30.8H, V1.8B
      uni_qpel_hv_filter_4lines x0, x2, V0, V1, x13, x12
      ret lr
endfunc

function put_vvc_uni_qpel_hv_10_16_neon_8lines, export=1
      mov    w11, #512 //offset val
      mov    w10, #0x03FF //max val
      DUP    V29.4S, w11 //offset vector
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      mov    x12, #2
      sub    x2, x2, #6//start h src[]
      sub    x13, x3, #7
      lsl    x13, x13, #1 // last jump h
      lsl    x1, x1, #1 //_dststride
      lsl    x14, x3, #1 //_srcstride
      lsl    x15, x14, #1 //2*_srcstride
      sub    x2, x2, x14
      sub    x2, x2, x15//start hv src[]
      mov    x16, x4//save height
      mov    x17, x2//save start src
      mov    x15, x0//save start dst
      LD1    {V0.8B}, [x5]
      LD1    {V1.8B}, [x6]
      SADDW   V0.8H, V30.8H, V0.8B
      SADDW   V1.8H, V30.8H, V1.8B
      uni_qpel_hv_filter_8lines x0, x2, V0, V1, x13, x12
      add x0, x15, #16
      add x2, x17, #16
      mov x4, x16
      uni_qpel_hv_filter_8lines x0, x2, V0, V1, x13, x12
      ret lr
endfunc

function put_vvc_uni_qpel_hv_10_16_neon_4lines, export=1
      mov    w11, #512 //offset val
      mov    w10, #0x03FF //max val
      DUP    V29.4S, w11 //offset vector
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      mov    x12, #2
      sub    x2, x2, #6//start h src[]
      sub    x13, x3, #7
      lsl    x13, x13, #1 // last jump h
      lsl    x1, x1, #1 //_dststride
      lsl    x14, x3, #1 //_srcstride
      lsl    x15, x14, #1 //2*_srcstride
      sub    x2, x2, x14
      sub    x2, x2, x15//start hv src[]
      mov    x16, x4//save height
      mov    x17, x2//save start src
      mov    x15, x0//save start dst
      LD1    {V0.8B}, [x5]
      LD1    {V1.8B}, [x6]
      SADDW   V0.8H, V30.8H, V0.8B
      SADDW   V1.8H, V30.8H, V1.8B
      uni_qpel_hv_filter_4lines x0, x2, V0, V1, x13, x12
      add x0, x15, #16
      add x2, x17, #16
      mov x4, x16
      uni_qpel_hv_filter_4lines x0, x2, V0, V1, x13, x12
      ret lr
endfunc

function put_vvc_uni_qpel_hv_10_32_neon_8lines, export=1
      mov    w11, #512 //offset val
      mov    w10, #0x03FF //max val
      DUP    V29.4S, w11 //offset vector
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      mov    x12, #2
      sub    x2, x2, #6//start h src[]
      sub    x13, x3, #7
      lsl    x13, x13, #1 // last jump h
      lsl    x1, x1, #1 //_dststride
      lsl    x14, x3, #1 //_srcstride
      lsl    x15, x14, #1 //2*_srcstride
      sub    x2, x2, x14
      sub    x2, x2, x15//start hv src[]
      mov    x16, x4//save height
      mov    x17, x2//save start src
      mov    x15, x0//save start dst
      LD1    {V0.8B}, [x5]
      LD1    {V1.8B}, [x6]
      SADDW   V0.8H, V30.8H, V0.8B
      SADDW   V1.8H, V30.8H, V1.8B
      uni_qpel_hv_filter_8lines x0, x2, V0, V1, x13, x12
      add x0, x15, #16
      add x2, x17, #16
      mov x4, x16
      uni_qpel_hv_filter_8lines x0, x2, V0, V1, x13, x12
      add x0, x15, #32
      add x2, x17, #32
      mov x4, x16
      uni_qpel_hv_filter_8lines x0, x2, V0, V1, x13, x12
      add x0, x15, #48
      add x2, x17, #48
      mov x4, x16
      uni_qpel_hv_filter_8lines x0, x2, V0, V1, x13, x12
      ret lr
endfunc

function put_vvc_uni_qpel_hv_10_32_neon_4lines, export=1
      mov    w11, #512 //offset val
      mov    w10, #0x03FF //max val
      DUP    V29.4S, w11 //offset vector
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      mov    x12, #2
      sub    x2, x2, #6//start h src[]
      sub    x13, x3, #7
      lsl    x13, x13, #1 // last jump h
      lsl    x1, x1, #1 //_dststride
      lsl    x14, x3, #1 //_srcstride
      lsl    x15, x14, #1 //2*_srcstride
      sub    x2, x2, x14
      sub    x2, x2, x15//start hv src[]
      mov    x16, x4//save height
      mov    x17, x2//save start src
      mov    x15, x0//save start dst
      LD1    {V0.8B}, [x5]
      LD1    {V1.8B}, [x6]
      SADDW   V0.8H, V30.8H, V0.8B
      SADDW   V1.8H, V30.8H, V1.8B
      uni_qpel_hv_filter_4lines x0, x2, V0, V1, x13, x12
      add x0, x15, #16
      add x2, x17, #16
      mov x4, x16
      uni_qpel_hv_filter_4lines x0, x2, V0, V1, x13, x12
      add x0, x15, #32
      add x2, x17, #32
      mov x4, x16
      uni_qpel_hv_filter_4lines x0, x2, V0, V1, x13, x12
      add x0, x15, #48
      add x2, x17, #48
      mov x4, x16
      uni_qpel_hv_filter_4lines x0, x2, V0, V1, x13, x12
      ret lr
endfunc

function put_vvc_bi0_qpel_h_10_8_neon, export=1
      mov    x12, #2
      mov    x11, MAX_PB_SIZE_16 //offset val
      mov    w10, #0x03FF //max val
      sub    x1, x1, #6
      sub    x2, x2, #7
      lsl    x2, x2, #1 //_srcstride
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      LD1    {V0.8B}, [x4]
      SADDW   V0.8H, V30.8H, V0.8B
0:
      subs x3, x3, #2
      bi0_qpel_h_filter V17.8H, x2
      bi0_qpel_h_filter V18.8H, x2
      ST1   {V17.8H}, [x0], x11
      ST1   {V18.8H}, [x0], x11
      b.ne 0b
      ret lr
endfunc

function put_vvc_bi0_qpel_h_10_16_neon, export=1
      mov    x12, #2
      mov    x11, MAX_PB_SIZE_16 //offset val
      mov    w10, #0x03FF //max val
      sub    x1, x1, #6
      sub    x2, x2, #15
      lsl    x2, x2, #1 //_srcstride
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      LD1    {V0.8B}, [x4]
      SADDW   V0.8H, V30.8H, V0.8B
0:
      subs x3, x3, #1
      bi0_qpel_h_filter V17.8H, x12
      bi0_qpel_h_filter V18.8H, x2
      ST1   {V17.8H-V18.8H}, [x0], x11
      b.ne 0b
      ret lr
endfunc

function put_vvc_bi0_qpel_h_10_32_neon, export=1
      mov    x12, #2
      mov    x11, MAX_PB_SIZE_16 //offset val
      mov    w10, #0x03FF //max val
      sub    x1, x1, #6
      sub    x2, x2, #31
      lsl    x2, x2, #1 //_srcstride
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      LD1    {V0.8B}, [x4]
      SADDW   V0.8H, V30.8H, V0.8B
0:
      subs x3, x3, #1
      bi0_qpel_h_filter V17.8H, x12
      bi0_qpel_h_filter V18.8H, x12
      bi0_qpel_h_filter V19.8H, x12
      bi0_qpel_h_filter V20.8H, x2
      ST1   {V17.8H-V20.8H}, [x0], x11
      b.ne 0b
      ret lr
endfunc

.macro bi0_qpel_v_filter_8lines, out, in
      LD1    {V1.8H}, [\in], x2
      LD1    {V2.8H}, [\in], x2
      LD1    {V3.8H}, [\in], x2
      LD1    {V4.8H}, [\in], x2
      LD1    {V5.8H}, [\in], x2
      LD1    {V6.8H}, [\in], x2
      LD1    {V7.8H}, [\in], x2
0:              subs x3, x3, #8
      LD1    {V8.8H}, [\in], x2
      bi0_qpel_filter V17.8H, V0, V1, V2, V3, V4, V5, V6, V7, V8
      LD1    {V1.8H}, [\in], x2
      bi0_qpel_filter V18.8H, V0, V2, V3, V4, V5, V6, V7, V8, V1
      LD1    {V2.8H}, [\in], x2
      bi0_qpel_filter V19.8H, V0, V3, V4, V5, V6, V7, V8, V1, V2
      LD1    {V3.8H}, [\in], x2
      bi0_qpel_filter V20.8H, V0, V4, V5, V6, V7, V8, V1, V2, V3
      LD1    {V4.8H}, [\in], x2
      bi0_qpel_filter V21.8H, V0, V5, V6, V7, V8, V1, V2, V3, V4
      LD1    {V5.8H}, [\in], x2
      bi0_qpel_filter V22.8H, V0, V6, V7, V8, V1, V2, V3, V4, V5
      LD1    {V6.8H}, [\in], x2
      bi0_qpel_filter V23.8H, V0, V7, V8, V1, V2, V3, V4, V5, V6
      LD1    {V7.8H}, [\in], x2
      bi0_qpel_filter V24.8H, V0, V8, V1, V2, V3, V4, V5, V6, V7
      ST1   {V17.8H}, [\out], x11
      ST1   {V18.8H}, [\out], x11
      ST1   {V19.8H}, [\out], x11
      ST1   {V20.8H}, [\out], x11
      ST1   {V21.8H}, [\out], x11
      ST1   {V22.8H}, [\out], x11
      ST1   {V23.8H}, [\out], x11
      ST1   {V24.8H}, [\out], x11
      b.ne 0b
.endm

.macro bi0_qpel_v_filter_4lines, out, in
      LD1    {V1.8H}, [\in], x2
      LD1    {V2.8H}, [\in], x2
      LD1    {V3.8H}, [\in], x2
      LD1    {V4.8H}, [\in], x2
      LD1    {V5.8H}, [\in], x2
      LD1    {V6.8H}, [\in], x2
      LD1    {V7.8H}, [\in], x2
0:              subs x3, x3, #4
      LD1    {V8.8H}, [\in], x2
      bi0_qpel_filter V17.8H, V0, V1, V2, V3, V4, V5, V6, V7, V8
      LD1    {V1.8H}, [\in], x2
      bi0_qpel_filter V18.8H, V0, V2, V3, V4, V5, V6, V7, V8, V1
      LD1    {V2.8H}, [\in], x2
      bi0_qpel_filter V19.8H, V0, V3, V4, V5, V6, V7, V8, V1, V2
      LD1    {V3.8H}, [\in], x2
      bi0_qpel_filter V20.8H, V0, V4, V5, V6, V7, V8, V1, V2, V3
      ST1   {V17.8H}, [\out], x11
      ST1   {V18.8H}, [\out], x11
      ST1   {V19.8H}, [\out], x11
      ST1   {V20.8H}, [\out], x11
      b.ne 0b
.endm

function put_vvc_bi0_qpel_v_10_8_neon_8lines, export=1
      mov    x11, MAX_PB_SIZE_16 //offset val
      mov    w10, #0x03FF //max val
      lsl    x2, x2, #1 //_srcstride
      lsl    x12, x2, #1 //2*_srcstride
      sub    x1, x1, x2
      sub    x1, x1, x12
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      LD1    {V0.8B}, [x4]
      SADDW   V0.8H, V30.8H, V0.8B
      bi0_qpel_v_filter_8lines x0, x1
      ret lr
endfunc

function put_vvc_bi0_qpel_v_10_8_neon_4lines, export=1
      mov    x11, MAX_PB_SIZE_16 //offset val
      mov    w10, #0x03FF //max val
      lsl    x2, x2, #1 //_srcstride
      lsl    x12, x2, #1 //2*_srcstride
      sub    x1, x1, x2
      sub    x1, x1, x12
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      LD1    {V0.8B}, [x4]
      SADDW   V0.8H, V30.8H, V0.8B
      bi0_qpel_v_filter_4lines x0, x1
      ret lr
endfunc

function put_vvc_bi0_qpel_v_10_16_neon_8lines, export=1
      mov    x11, MAX_PB_SIZE_16 //offset val
      mov    w10, #0x03FF //max val
      lsl    x2, x2, #1 //_srcstride
      lsl    x12, x2, #1 //2*_srcstride
      sub    x1, x1, x2
      sub    x1, x1, x12
      mov    x9, x1  //save x1[-3*_srcstride]
      mov    x8, x0  //save x0[0]
      mov    x14, x3  //save height
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      LD1    {V0.8B}, [x4]
      SADDW   V0.8H, V30.8H, V0.8B
      bi0_qpel_v_filter_8lines x0, x1
      add x0, x8, #16
      add x1, x9, #16
      mov x3, x14
      bi0_qpel_v_filter_8lines x0, x1
      ret lr
endfunc

function put_vvc_bi0_qpel_v_10_16_neon_4lines, export=1
      mov    x11, MAX_PB_SIZE_16 //offset val
      mov    w10, #0x03FF //max val
      lsl    x2, x2, #1 //_srcstride
      lsl    x12, x2, #1 //2*_srcstride
      sub    x1, x1, x2
      sub    x1, x1, x12
      mov    x9, x1  //save x1[-3*_srcstride]
      mov    x8, x0  //save x0[0]
      mov    x14, x3  //save height
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      LD1    {V0.8B}, [x4]
      SADDW   V0.8H, V30.8H, V0.8B
      bi0_qpel_v_filter_4lines x0, x1
      add x0, x8, #16
      add x1, x9, #16
      mov x3, x14
      bi0_qpel_v_filter_4lines x0, x1
      ret lr
endfunc

function put_vvc_bi0_qpel_v_10_32_neon_8lines, export=1
      mov    x11, MAX_PB_SIZE_16 //offset val
      mov    w10, #0x03FF //max val
      lsl    x2, x2, #1 //_srcstride
      lsl    x12, x2, #1 //2*_srcstride
      sub    x1, x1, x2
      sub    x1, x1, x12
      mov    x9, x1  //save x1[-3*_srcstride]
      mov    x8, x0  //save x0[0]
      mov    x14, x3  //save height
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      LD1    {V0.8B}, [x4]
      SADDW   V0.8H, V30.8H, V0.8B
      bi0_qpel_v_filter_8lines x0, x1
      add x0, x8, #16
      add x1, x9, #16
      mov x3, x14
      bi0_qpel_v_filter_8lines x0, x1
      add x0, x8, #32
      add x1, x9, #32
      mov x3, x14
      bi0_qpel_v_filter_8lines x0, x1
      add x0, x8, #48
      add x1, x9, #48
      mov x3, x14
      bi0_qpel_v_filter_8lines x0, x1
      ret lr
endfunc

function put_vvc_bi0_qpel_v_10_32_neon_4lines, export=1
      mov    x11, MAX_PB_SIZE_16 //offset val
      mov    w10, #0x03FF //max val
      lsl    x2, x2, #1 //_srcstride
      lsl    x12, x2, #1 //2*_srcstride
      sub    x1, x1, x2
      sub    x1, x1, x12
      mov    x9, x1  //save x1[-3*_srcstride]
      mov    x8, x0  //save x0[0]
      mov    x14, x3  //save height
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      LD1    {V0.8B}, [x4]
      SADDW   V0.8H, V30.8H, V0.8B
      bi0_qpel_v_filter_4lines x0, x1
      add x0, x8, #16
      add x1, x9, #16
      mov x3, x14
      bi0_qpel_v_filter_4lines x0, x1
      add x0, x8, #32
      add x1, x9, #32
      mov x3, x14
      bi0_qpel_v_filter_4lines x0, x1
      add x0, x8, #48
      add x1, x9, #48
      mov x3, x14
      bi0_qpel_v_filter_4lines x0, x1
      ret lr
endfunc

.macro bi0_qpel_hv_filter_8lines, out, in, coefx=V0, coefy=V1, last_jump=x3, jump=x12
      bi0_qpel_h_filter V2.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
      bi0_qpel_h_filter V3.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
      bi0_qpel_h_filter V4.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
      bi0_qpel_h_filter V5.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
      bi0_qpel_h_filter V6.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
      bi0_qpel_h_filter V7.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
      bi0_qpel_h_filter V8.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
0:              subs x3, x3, #8
      bi0_qpel_h_filter V9.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
      bi0_qpel_filter V20.8H, \coefy, V2, V3, V4, V5, V6, V7, V8, V9, shift=#6
      bi0_qpel_h_filter V2.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
      bi0_qpel_filter V21.8H, \coefy, V3, V4, V5, V6, V7, V8, V9, V2, shift=#6
      bi0_qpel_h_filter V3.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
      bi0_qpel_filter V22.8H, \coefy, V4, V5, V6, V7, V8, V9, V2, V3, shift=#6
      bi0_qpel_h_filter V4.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
      bi0_qpel_filter V23.8H, \coefy, V5, V6, V7, V8, V9, V2, V3, V4, shift=#6
      bi0_qpel_h_filter V5.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
      bi0_qpel_filter V24.8H, \coefy, V6, V7, V8, V9, V2, V3, V4, V5, shift=#6
      bi0_qpel_h_filter V6.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
      bi0_qpel_filter V25.8H, \coefy, V7, V8, V9, V2, V3, V4, V5, V6, shift=#6
      bi0_qpel_h_filter V7.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
      bi0_qpel_filter V26.8H, \coefy, V8, V9, V2, V3, V4, V5, V6, V7, shift=#6
      bi0_qpel_h_filter V8.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
      bi0_qpel_filter V27.8H, \coefy, V9, V2, V3, V4, V5, V6, V7, V8, shift=#6
      ST1   {V20.8H}, [\out], x11
      ST1   {V21.8H}, [\out], x11
      ST1   {V22.8H}, [\out], x11
      ST1   {V23.8H}, [\out], x11
      ST1   {V24.8H}, [\out], x11
      ST1   {V25.8H}, [\out], x11
      ST1   {V26.8H}, [\out], x11
      ST1   {V27.8H}, [\out], x11
      b.ne 0b
.endm

.macro bi0_qpel_hv_filter_4lines, out, in, coefx=V0, coefy=V1, last_jump=x3, jump=x12
      bi0_qpel_h_filter V2.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
      bi0_qpel_h_filter V3.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
      bi0_qpel_h_filter V4.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
      bi0_qpel_h_filter V5.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
      bi0_qpel_h_filter V6.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
      bi0_qpel_h_filter V7.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
      bi0_qpel_h_filter V8.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
0:              subs x3, x3, #4
      bi0_qpel_h_filter V9.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
      bi0_qpel_filter V20.8H, \coefy, V2, V3, V4, V5, V6, V7, V8, V9, shift=#6
      bi0_qpel_h_filter V2.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
      bi0_qpel_filter V21.8H, \coefy, V3, V4, V5, V6, V7, V8, V9, V2, shift=#6
      bi0_qpel_h_filter V3.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
      bi0_qpel_filter V22.8H, \coefy, V4, V5, V6, V7, V8, V9, V2, V3, shift=#6
      bi0_qpel_h_filter V4.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
      bi0_qpel_filter V23.8H, \coefy, V5, V6, V7, V8, V9, V2, V3, V4, shift=#6
      ST1   {V20.8H}, [\out], x11
      ST1   {V21.8H}, [\out], x11
      ST1   {V22.8H}, [\out], x11
      ST1   {V23.8H}, [\out], x11
      b.ne 0b
.endm

function put_vvc_bi0_qpel_hv_10_8_neon_8lines, export=1
      mov    x11, MAX_PB_SIZE_16
      mov    w10, #0x03FF //max val
      DUP    V29.4S, w11 //offset vector
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      mov    x12, #2
      sub    x1, x1, #6//start h src[]
      sub    x13, x2, #7
      lsl    x13, x13, #1 // last jump h
      lsl    x14, x2, #1 //_srcstride
      lsl    x15, x14, #1 //2*_srcstride
      sub    x1, x1, x14
      sub    x1, x1, x15//start hv src[]
      LD1    {V0.8B}, [x4]
      LD1    {V1.8B}, [x5]
      SADDW   V0.8H, V30.8H, V0.8B
      SADDW   V1.8H, V30.8H, V1.8B
      bi0_qpel_hv_filter_8lines x0, x1, V0, V1, x13, x12
      ret lr
endfunc

function put_vvc_bi0_qpel_hv_10_8_neon_4lines, export=1
      mov    x11, MAX_PB_SIZE_16
      mov    w10, #0x03FF //max val
      DUP    V29.4S, w11 //offset vector
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      mov    x12, #2
      sub    x1, x1, #6//start h src[]
      sub    x13, x2, #7
      lsl    x13, x13, #1 // last jump h
      lsl    x14, x2, #1 //_srcstride
      lsl    x15, x14, #1 //2*_srcstride
      sub    x1, x1, x14
      sub    x1, x1, x15//start hv src[]
      LD1    {V0.8B}, [x4]
      LD1    {V1.8B}, [x5]
      SADDW   V0.8H, V30.8H, V0.8B
      SADDW   V1.8H, V30.8H, V1.8B
      bi0_qpel_hv_filter_4lines x0, x1, V0, V1, x13, x12
      ret lr
endfunc

function put_vvc_bi0_qpel_hv_10_16_neon_8lines, export=1
      mov    x11, MAX_PB_SIZE_16
      mov    w10, #0x03FF //max val
      DUP    V29.4S, w11 //offset vector
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      mov    x12, #2
      sub    x1, x1, #6//start h src[]
      sub    x13, x2, #7
      lsl    x13, x13, #1 // last jump h
      lsl    x14, x2, #1 //_srcstride
      lsl    x15, x14, #1 //2*_srcstride
      sub    x1, x1, x14
      sub    x1, x1, x15//start hv src[]
      mov    x16, x3//save height
      mov    x17, x1//save start src
      mov    x15, x0//save start dst
      LD1    {V0.8B}, [x4]
      LD1    {V1.8B}, [x5]
      SADDW   V0.8H, V30.8H, V0.8B
      SADDW   V1.8H, V30.8H, V1.8B
      bi0_qpel_hv_filter_8lines x0, x1, V0, V1, x13, x12
      add x0, x15, #16
      add x1, x17, #16
      mov x3, x16
      bi0_qpel_hv_filter_8lines x0, x1, V0, V1, x13, x12
      ret lr
endfunc

function put_vvc_bi0_qpel_hv_10_16_neon_4lines, export=1
      mov    x11, MAX_PB_SIZE_16
      mov    w10, #0x03FF //max val
      DUP    V29.4S, w11 //offset vector
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      mov    x12, #2
      sub    x1, x1, #6//start h src[]
      sub    x13, x2, #7
      lsl    x13, x13, #1 // last jump h
      lsl    x14, x2, #1 //_srcstride
      lsl    x15, x14, #1 //2*_srcstride
      sub    x1, x1, x14
      sub    x1, x1, x15//start hv src[]
      mov    x16, x3//save height
      mov    x17, x1//save start src
      mov    x15, x0//save start dst
      LD1    {V0.8B}, [x4]
      LD1    {V1.8B}, [x5]
      SADDW   V0.8H, V30.8H, V0.8B
      SADDW   V1.8H, V30.8H, V1.8B
      bi0_qpel_hv_filter_4lines x0, x1, V0, V1, x13, x12
      add x0, x15, #16
      add x1, x17, #16
      mov x3, x16
      bi0_qpel_hv_filter_4lines x0, x1, V0, V1, x13, x12
      ret lr
endfunc

function put_vvc_bi0_qpel_hv_10_32_neon_8lines, export=1
      mov    x11, MAX_PB_SIZE_16
      mov    w10, #0x03FF //max val
      DUP    V29.4S, w11 //offset vector
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      mov    x12, #2
      sub    x1, x1, #6//start h src[]
      sub    x13, x2, #7
      lsl    x13, x13, #1 // last jump h
      lsl    x14, x2, #1 //_srcstride
      lsl    x15, x14, #1 //2*_srcstride
      sub    x1, x1, x14
      sub    x1, x1, x15//start hv src[]
      mov    x16, x3//save height
      mov    x17, x1//save start src
      mov    x15, x0//save start dst
      LD1    {V0.8B}, [x4]
      LD1    {V1.8B}, [x5]
      SADDW   V0.8H, V30.8H, V0.8B
      SADDW   V1.8H, V30.8H, V1.8B
      bi0_qpel_hv_filter_8lines x0, x1, V0, V1, x13, x12
      add x0, x15, #16
      add x1, x17, #16
      mov x3, x16
      bi0_qpel_hv_filter_8lines x0, x1, V0, V1, x13, x12
      add x0, x15, #32
      add x1, x17, #32
      mov x3, x16
      bi0_qpel_hv_filter_8lines x0, x1, V0, V1, x13, x12
      add x0, x15, #48
      add x1, x17, #48
      mov x3, x16
      bi0_qpel_hv_filter_8lines x0, x1, V0, V1, x13, x12
      ret lr
endfunc

function put_vvc_bi0_qpel_hv_10_32_neon_4lines, export=1
      mov    x11, MAX_PB_SIZE_16
      mov    w10, #0x03FF //max val
      DUP    V29.4S, w11 //offset vector
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      mov    x12, #2
      sub    x1, x1, #6//start h src[]
      sub    x13, x2, #7
      lsl    x13, x13, #1 // last jump h
      lsl    x14, x2, #1 //_srcstride
      lsl    x15, x14, #1 //2*_srcstride
      sub    x1, x1, x14
      sub    x1, x1, x15//start hv src[]
      mov    x16, x3//save height
      mov    x17, x1//save start src
      mov    x15, x0//save start dst
      LD1    {V0.8B}, [x4]
      LD1    {V1.8B}, [x5]
      SADDW   V0.8H, V30.8H, V0.8B
      SADDW   V1.8H, V30.8H, V1.8B
      bi0_qpel_hv_filter_4lines x0, x1, V0, V1, x13, x12
      add x0, x15, #16
      add x1, x17, #16
      mov x3, x16
      bi0_qpel_hv_filter_4lines x0, x1, V0, V1, x13, x12
      add x0, x15, #32
      add x1, x17, #32
      mov x3, x16
      bi0_qpel_hv_filter_4lines x0, x1, V0, V1, x13, x12
      add x0, x15, #48
      add x1, x17, #48
      mov x3, x16
      bi0_qpel_hv_filter_4lines x0, x1, V0, V1, x13, x12
      ret lr
endfunc

.macro bi1_qpel_h_filter, out, last_jump, src1, suffix1=.4H, suffix2=.8H
      LD1   {V1.8H}, [x2], x12
      LD1   {V2.8H}, [x2], x12
      LD1   {V3.8H}, [x2], x12
      LD1   {V4.8H}, [x2], x12
      LD1   {V5.8H}, [x2], x12
      LD1   {V6.8H}, [x2], x12
      LD1   {V7.8H}, [x2], x12
      LD1   {V8.8H}, [x2], \last_jump
      bi1_qpel_filter \out, \src1, V0
.endm

function put_vvc_bi1_qpel_h_10_8_neon, export=1
      mov    x12, #2
      mov    w11, #64 //offset val
      mov    w10, #0x03FF //max val
      mov    x9, MAX_PB_SIZE_16
      sub    x2, x2, #6
      sub    x3, x3, #7
      lsl    x1, x1, #1 //_dststride
      lsl    x3, x3, #1 //_srcstride
      DUP    V29.4S, w11 //offset vector
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      LD1    {V0.8B}, [x6]
      SADDW   V0.8H, V30.8H, V0.8B
0:
      subs x5, x5, #2
      LD1   {V27.8H}, [x4], x9
      LD1   {V28.8H}, [x4], x9
      bi1_qpel_h_filter V17.8H, x3, V27
      bi1_qpel_h_filter V18.8H, x3, V28
      ST1   {V17.8H}, [x0], x1
      ST1   {V18.8H}, [x0], x1
      b.ne 0b
      ret lr
endfunc

function put_vvc_bi1_qpel_h_10_16_neon, export=1
      mov    x12, #2
      mov    w11, #64 //offset val
      mov    w10, #0x03FF //max val
      mov    x9, MAX_PB_SIZE_16
      sub    x2, x2, #6
      sub    x3, x3, #15
      lsl    x1, x1, #1 //_dststride
      lsl    x3, x3, #1 //_srcstride
      DUP    V29.4S, w11 //offset vector
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      LD1    {V0.8B}, [x6]
      SADDW   V0.8H, V30.8H, V0.8B
0:
      subs x5, x5, #1
      LD1   {V27.8H-V28.8H}, [x4], x9
      bi1_qpel_h_filter V17.8H, x12, V27
      bi1_qpel_h_filter V18.8H, x3, V28
      ST1   {V17.8H-V18.8H}, [x0], x1
      b.ne 0b
      ret lr
endfunc

function put_vvc_bi1_qpel_h_10_32_neon, export=1
      mov    x12, #2
      mov    w11, #64 //offset val
      mov    w10, #0x03FF //max val
      mov    x9, MAX_PB_SIZE_16
      sub    x2, x2, #6
      sub    x3, x3, #31
      lsl    x1, x1, #1 //_dststride
      lsl    x3, x3, #1 //_srcstride
      DUP    V29.4S, w11 //offset vector
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      LD1    {V0.8B}, [x6]
      SADDW   V0.8H, V30.8H, V0.8B
0:
      subs x5, x5, #1
      LD1   {V25.8H-V28.8H}, [x4], x9
      bi1_qpel_h_filter V17.8H, x12, V25
      bi1_qpel_h_filter V18.8H, x12, V26
      bi1_qpel_h_filter V19.8H, x12, V27
      bi1_qpel_h_filter V20.8H, x3, V28
      ST1   {V17.8H-V20.8H}, [x0], x1
      b.ne 0b
      ret lr
endfunc

.macro bi1_qpel_v_filter_8lines, out, in
      LD1    {V1.8H}, [\in], x3
      LD1    {V2.8H}, [\in], x3
      LD1    {V3.8H}, [\in], x3
      LD1    {V4.8H}, [\in], x3
      LD1    {V5.8H}, [\in], x3
      LD1    {V6.8H}, [\in], x3
      LD1    {V7.8H}, [\in], x3
0:              subs   x5, x5, #8
      LD1    {V8.8H}, [\in], x3
      LD1    {V28.8H}, [x4], x9
      bi1_qpel_filter V17.8H, v28, V0, V1, V2, V3, V4, V5, V6, V7, V8
      LD1    {V1.8H}, [\in], x3
      LD1    {V28.8H}, [x4], x9
      bi1_qpel_filter V18.8H, v28, V0, V2, V3, V4, V5, V6, V7, V8, V1
      LD1    {V2.8H}, [\in], x3
      LD1    {V28.8H}, [x4], x9
      bi1_qpel_filter V19.8H, v28, V0, V3, V4, V5, V6, V7, V8, V1, V2
      LD1    {V3.8H}, [\in], x3
      LD1    {V28.8H}, [x4], x9
      bi1_qpel_filter V20.8H, v28, V0, V4, V5, V6, V7, V8, V1, V2, V3
      LD1    {V4.8H}, [\in], x3
      LD1    {V28.8H}, [x4], x9
      bi1_qpel_filter V21.8H, v28, V0, V5, V6, V7, V8, V1, V2, V3, V4
      LD1    {V5.8H}, [\in], x3
      LD1    {V28.8H}, [x4], x9
      bi1_qpel_filter V22.8H, v28, V0, V6, V7, V8, V1, V2, V3, V4, V5
      LD1    {V6.8H}, [\in], x3
      LD1    {V28.8H}, [x4], x9
      bi1_qpel_filter V23.8H, v28, V0, V7, V8, V1, V2, V3, V4, V5, V6
      LD1    {V7.8H}, [\in], x3
      LD1    {V28.8H}, [x4], x9
      bi1_qpel_filter V24.8H, v28, V0, V8, V1, V2, V3, V4, V5, V6, V7
      ST1   {V17.8H}, [\out], x1
      ST1   {V18.8H}, [\out], x1
      ST1   {V19.8H}, [\out], x1
      ST1   {V20.8H}, [\out], x1
      ST1   {V21.8H}, [\out], x1
      ST1   {V22.8H}, [\out], x1
      ST1   {V23.8H}, [\out], x1
      ST1   {V24.8H}, [\out], x1
      b.ne 0b
.endm

.macro bi1_qpel_v_filter_4lines, out, in
      LD1    {V1.8H}, [\in], x3
      LD1    {V2.8H}, [\in], x3
      LD1    {V3.8H}, [\in], x3
      LD1    {V4.8H}, [\in], x3
      LD1    {V5.8H}, [\in], x3
      LD1    {V6.8H}, [\in], x3
      LD1    {V7.8H}, [\in], x3
0:              subs   x5, x5, #4
      LD1    {V8.8H}, [\in], x3
      LD1    {V28.8H}, [x4], x9
      bi1_qpel_filter V17.8H, v28, V0, V1, V2, V3, V4, V5, V6, V7, V8
      LD1    {V1.8H}, [\in], x3
      LD1    {V28.8H}, [x4], x9
      bi1_qpel_filter V18.8H, v28, V0, V2, V3, V4, V5, V6, V7, V8, V1
      LD1    {V2.8H}, [\in], x3
      LD1    {V28.8H}, [x4], x9
      bi1_qpel_filter V19.8H, v28, V0, V3, V4, V5, V6, V7, V8, V1, V2
      LD1    {V3.8H}, [\in], x3
      LD1    {V28.8H}, [x4], x9
      bi1_qpel_filter V20.8H, v28, V0, V4, V5, V6, V7, V8, V1, V2, V3
      ST1   {V17.8H}, [\out], x1
      ST1   {V18.8H}, [\out], x1
      ST1   {V19.8H}, [\out], x1
      ST1   {V20.8H}, [\out], x1
      b.ne 0b
.endm

function put_vvc_bi1_qpel_v_10_8_neon_8lines, export=1
      mov    x12, #2
      mov    w11, #64 //offset val
      mov    w10, #0x03FF //max val
      mov    x9, MAX_PB_SIZE_16
      lsl    x1, x1, #1 //_dststride
      lsl    x3, x3, #1 //_srcstride
      lsl    x13, x3, #1 //2*_srcstride
      sub    x2, x2, x3
      sub    x2, x2, x13
      DUP    V29.4S, w11 //offset vector
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      LD1    {V0.8B}, [x6]
      SADDW   V0.8H, V30.8H, V0.8B
      bi1_qpel_v_filter_8lines x0, x2
      ret lr
endfunc

function put_vvc_bi1_qpel_v_10_8_neon_4lines, export=1
      mov    x12, #2
      mov    w11, #64 //offset val
      mov    w10, #0x03FF //max val
      mov    x9, MAX_PB_SIZE_16
      lsl    x1, x1, #1 //_dststride
      lsl    x3, x3, #1 //_srcstride
      lsl    x13, x3, #1 //2*_srcstride
      sub    x2, x2, x3
      sub    x2, x2, x13
      DUP    V29.4S, w11 //offset vector
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      LD1    {V0.8B}, [x6]
      SADDW   V0.8H, V30.8H, V0.8B
      bi1_qpel_v_filter_4lines x0, x2
      ret lr
endfunc

function put_vvc_bi1_qpel_v_10_16_neon_8lines, export=1
      mov    w11, #64 //offset val
      mov    w10, #0x03FF //max val
      mov    x9, MAX_PB_SIZE_16
      lsl    x1, x1, #1 //_dststride
      lsl    x3, x3, #1 //_srcstride
      lsl    x13, x3, #1 //2*_srcstride
      sub    x2, x2, x3
      sub    x2, x2, x13
      mov    x13, x0  //save x0[0]
      mov    x12, x2  //save x2[-3*_srcstride]
      mov    x14, x4  //save x4[0]
      mov    x15, x5  //save height
      DUP    V29.4S, w11 //offset vector
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      LD1    {V0.8B}, [x6]
      SADDW   V0.8H, V30.8H, V0.8B
      bi1_qpel_v_filter_8lines x0, x2
      add x2, x12, #16
      add x0, x13, #16
      add x4, x14, #16
      mov x5, x15
      bi1_qpel_v_filter_8lines x0, x2
      ret lr
endfunc

function put_vvc_bi1_qpel_v_10_16_neon_4lines, export=1
      mov    w11, #64 //offset val
      mov    w10, #0x03FF //max val
      mov    x9, MAX_PB_SIZE_16
      lsl    x1, x1, #1 //_dststride
      lsl    x3, x3, #1 //_srcstride
      lsl    x13, x3, #1 //2*_srcstride
      sub    x2, x2, x3
      sub    x2, x2, x13
      mov    x13, x0  //save x0[0]
      mov    x12, x2  //save x2[-3*_srcstride]
      mov    x14, x4  //save x4[0]
      mov    x15, x5  //save height
      DUP    V29.4S, w11 //offset vector
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      LD1    {V0.8B}, [x6]
      SADDW   V0.8H, V30.8H, V0.8B
      bi1_qpel_v_filter_4lines x0, x2
      add x2, x12, #16
      add x0, x13, #16
      add x4, x14, #16
      mov x5, x15
      bi1_qpel_v_filter_4lines x0, x2
      ret lr
endfunc

function put_vvc_bi1_qpel_v_10_32_neon_8lines, export=1
      mov    w11, #64 //offset val
      mov    w10, #0x03FF //max val
      mov    x9, MAX_PB_SIZE_16
      lsl    x1, x1, #1 //_dststride
      lsl    x3, x3, #1 //_srcstride
      lsl    x13, x3, #1 //2*_srcstride
      sub    x2, x2, x3
      sub    x2, x2, x13
      mov    x12, x2  //save x2[-3*_srcstride]
      mov    x13, x0  //save x0[0]
      mov    x14, x4  //save x4[0]
      mov    x15, x5  //save height
      DUP    V29.4S, w11 //offset vector
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      LD1    {V0.8B}, [x6]
      SADDW   V0.8H, V30.8H, V0.8B
      bi1_qpel_v_filter_8lines x0, x2
      add x2, x12, #16
      add x0, x13, #16
      add x4, x14, #16
      mov x5, x15
      bi1_qpel_v_filter_8lines x0, x2
      add x2, x12, #32
      add x0, x13, #32
      add x4, x14, #32
      mov x5, x15
      bi1_qpel_v_filter_8lines x0, x2
      add x2, x12, #48
      add x0, x13, #48
      add x4, x14, #48
      mov x5, x15
      bi1_qpel_v_filter_8lines x0, x2
      ret lr
endfunc

function put_vvc_bi1_qpel_v_10_32_neon_4lines, export=1
      mov    w11, #64 //offset val
      mov    w10, #0x03FF //max val
      mov    x9, MAX_PB_SIZE_16
      lsl    x1, x1, #1 //_dststride
      lsl    x3, x3, #1 //_srcstride
      lsl    x13, x3, #1 //2*_srcstride
      sub    x2, x2, x3
      sub    x2, x2, x13
      mov    x12, x2  //save x2[-3*_srcstride]
      mov    x13, x0  //save x0[0]
      mov    x14, x4  //save x4[0]
      mov    x15, x5  //save height
      DUP    V29.4S, w11 //offset vector
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      LD1    {V0.8B}, [x6]
      SADDW   V0.8H, V30.8H, V0.8B
      bi1_qpel_v_filter_4lines x0, x2
      add x2, x12, #16
      add x0, x13, #16
      add x4, x14, #16
      mov x5, x15
      bi1_qpel_v_filter_4lines x0, x2
      add x2, x12, #32
      add x0, x13, #32
      add x4, x14, #32
      mov x5, x15
      bi1_qpel_v_filter_4lines x0, x2
      add x2, x12, #48
      add x0, x13, #48
      add x4, x14, #48
      mov x5, x15
      bi1_qpel_v_filter_4lines x0, x2
      ret lr
endfunc

.macro bi1_qpel_hv_filter_8lines, out, in, coefx=V0, coefy=V1, last_jump=x3, jump=x12
      bi0_qpel_h_filter V2.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
      bi0_qpel_h_filter V3.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
      bi0_qpel_h_filter V4.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
      bi0_qpel_h_filter V5.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
      bi0_qpel_h_filter V6.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
      bi0_qpel_h_filter V7.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
      bi0_qpel_h_filter V8.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
0:              subs x5, x5, #8
      bi0_qpel_h_filter V9.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
      LD1    {V28.8H}, [x4], x9
      bi1_qpel_filter V20.8H, V28, \coefy, V2, V3, V4, V5, V6, V7, V8, V9, shift=#6
      bi0_qpel_h_filter V2.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
      LD1    {V28.8H}, [x4], x9
      bi1_qpel_filter V21.8H, V28, \coefy, V3, V4, V5, V6, V7, V8, V9, V2, shift=#6
      bi0_qpel_h_filter V3.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
      LD1    {V28.8H}, [x4], x9
      bi1_qpel_filter V22.8H, V28, \coefy, V4, V5, V6, V7, V8, V9, V2, V3, shift=#6
      bi0_qpel_h_filter V4.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
      LD1    {V28.8H}, [x4], x9
      bi1_qpel_filter V23.8H, V28, \coefy, V5, V6, V7, V8, V9, V2, V3, V4, shift=#6
      bi0_qpel_h_filter V5.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
      LD1    {V28.8H}, [x4], x9
      bi1_qpel_filter V24.8H, V28, \coefy, V6, V7, V8, V9, V2, V3, V4, V5, shift=#6
      bi0_qpel_h_filter V6.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
      LD1    {V28.8H}, [x4], x9
      bi1_qpel_filter V25.8H, V28, \coefy, V7, V8, V9, V2, V3, V4, V5, V6, shift=#6
      bi0_qpel_h_filter V7.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
      LD1    {V28.8H}, [x4], x9
      bi1_qpel_filter V26.8H, V28, \coefy, V8, V9, V2, V3, V4, V5, V6, V7, shift=#6
      bi0_qpel_h_filter V8.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
      LD1    {V28.8H}, [x4], x9
      bi1_qpel_filter V27.8H, V28, \coefy, V9, V2, V3, V4, V5, V6, V7, V8, shift=#6
      ST1   {V20.8H}, [\out], x1
      ST1   {V21.8H}, [\out], x1
      ST1   {V22.8H}, [\out], x1
      ST1   {V23.8H}, [\out], x1
      ST1   {V24.8H}, [\out], x1
      ST1   {V25.8H}, [\out], x1
      ST1   {V26.8H}, [\out], x1
      ST1   {V27.8H}, [\out], x1
      b.ne 0b
.endm

.macro bi1_qpel_hv_filter_4lines, out, in, coefx=V0, coefy=V1, last_jump=x3, jump=x12
      bi0_qpel_h_filter V2.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
      bi0_qpel_h_filter V3.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
      bi0_qpel_h_filter V4.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
      bi0_qpel_h_filter V5.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
      bi0_qpel_h_filter V6.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
      bi0_qpel_h_filter V7.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
      bi0_qpel_h_filter V8.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
0:              subs x5, x5, #4
      bi0_qpel_h_filter V9.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
      LD1    {V28.8H}, [x4], x9
      bi1_qpel_filter V20.8H, V28, \coefy, V2, V3, V4, V5, V6, V7, V8, V9, shift=#6
      bi0_qpel_h_filter V2.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
      LD1    {V28.8H}, [x4], x9
      bi1_qpel_filter V21.8H, V28, \coefy, V3, V4, V5, V6, V7, V8, V9, V2, shift=#6
      bi0_qpel_h_filter V3.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
      LD1    {V28.8H}, [x4], x9
      bi1_qpel_filter V22.8H, V28, \coefy, V4, V5, V6, V7, V8, V9, V2, V3, shift=#6
      bi0_qpel_h_filter V4.8H, \last_jump, \coefx, \jump, \in, V12, V13, V14, V15, V16, V17, V18, V19
      LD1    {V28.8H}, [x4], x9
      bi1_qpel_filter V23.8H, V28, \coefy, V5, V6, V7, V8, V9, V2, V3, V4, shift=#6
      ST1   {V20.8H}, [\out], x1
      ST1   {V21.8H}, [\out], x1
      ST1   {V22.8H}, [\out], x1
      ST1   {V23.8H}, [\out], x1
      b.ne 0b
.endm

function put_vvc_bi1_qpel_hv_10_8_neon_8lines, export=1
      mov    x9, MAX_PB_SIZE_16
      mov    w11, #1024 //offset val
      mov    w10, #0x03FF //max val
      DUP    V29.4S, w11 //offset vector
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      mov    x12, #2
      sub    x2, x2, #6//start h src[]
      sub    x13, x3, #7
      lsl    x13, x13, #1 // last jump h
      lsl    x1, x1, #1 //_dststride
      lsl    x14, x3, #1 //_srcstride
      lsl    x15, x14, #1 //2*_srcstride
      sub    x2, x2, x14
      sub    x2, x2, x15//start hv src[]
      LD1    {V0.8B}, [x6]
      LD1    {V1.8B}, [x7]
      SADDW   V0.8H, V30.8H, V0.8B
      SADDW   V1.8H, V30.8H, V1.8B
      bi1_qpel_hv_filter_8lines x0, x2, V0, V1, x13, x12
      ret lr
endfunc

function put_vvc_bi1_qpel_hv_10_8_neon_4lines, export=1
      mov    x9, MAX_PB_SIZE_16
      mov    w11, #1024 //offset val
      mov    w10, #0x03FF //max val
      DUP    V29.4S, w11 //offset vector
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      mov    x12, #2
      sub    x2, x2, #6//start h src[]
      sub    x13, x3, #7
      lsl    x13, x13, #1 // last jump h
      lsl    x1, x1, #1 //_dststride
      lsl    x14, x3, #1 //_srcstride
      lsl    x15, x14, #1 //2*_srcstride
      sub    x2, x2, x14
      sub    x2, x2, x15//start hv src[]
      LD1    {V0.8B}, [x6]
      LD1    {V1.8B}, [x7]
      SADDW   V0.8H, V30.8H, V0.8B
      SADDW   V1.8H, V30.8H, V1.8B
      bi1_qpel_hv_filter_4lines x0, x2, V0, V1, x13, x12
      ret lr
endfunc

function put_vvc_bi1_qpel_hv_10_16_neon_8lines, export=1
      mov    x9, MAX_PB_SIZE_16
      mov    w11, #1024 //offset val
      mov    w10, #0x03FF //max val
      DUP    V29.4S, w11 //offset vector
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      mov    x12, #2
      sub    x2, x2, #6//start h src[]
      sub    x13, x3, #7
      lsl    x13, x13, #1 // last jump h
      lsl    x1, x1, #1 //_dststride
      lsl    x14, x3, #1 //_srcstride
      lsl    x15, x14, #1 //2*_srcstride
      sub    x2, x2, x14
      sub    x2, x2, x15//start hv src[]
      mov    x16, x5//save height
      mov    x17, x2//save start src
      mov    x15, x0//save start dst
      mov    x14, x4//save start dst
      LD1    {V0.8B}, [x6]
      LD1    {V1.8B}, [x7]
      SADDW   V0.8H, V30.8H, V0.8B
      SADDW   V1.8H, V30.8H, V1.8B
      bi1_qpel_hv_filter_8lines x0, x2, V0, V1, x13, x12
      add x0, x15, #16
      add x2, x17, #16
      add x4, x14, #16
      mov x5, x16
      bi1_qpel_hv_filter_8lines x0, x2, V0, V1, x13, x12
      ret lr
endfunc

function put_vvc_bi1_qpel_hv_10_16_neon_4lines, export=1
      mov    x9, MAX_PB_SIZE_16
      mov    w11, #1024 //offset val
      mov    w10, #0x03FF //max val
      DUP    V29.4S, w11 //offset vector
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      mov    x12, #2
      sub    x2, x2, #6//start h src[]
      sub    x13, x3, #7
      lsl    x13, x13, #1 // last jump h
      lsl    x1, x1, #1 //_dststride
      lsl    x14, x3, #1 //_srcstride
      lsl    x15, x14, #1 //2*_srcstride
      sub    x2, x2, x14
      sub    x2, x2, x15//start hv src[]
      mov    x16, x5//save height
      mov    x17, x2//save start src
      mov    x15, x0//save start dst
      mov    x14, x4//save start dst
      LD1    {V0.8B}, [x6]
      LD1    {V1.8B}, [x7]
      SADDW   V0.8H, V30.8H, V0.8B
      SADDW   V1.8H, V30.8H, V1.8B
      bi1_qpel_hv_filter_4lines x0, x2, V0, V1, x13, x12
      add x0, x15, #16
      add x2, x17, #16
      add x4, x14, #16
      mov x5, x16
      bi1_qpel_hv_filter_4lines x0, x2, V0, V1, x13, x12
      ret lr
endfunc

function put_vvc_bi1_qpel_hv_10_32_neon_8lines, export=1
      mov    x9, MAX_PB_SIZE_16
      mov    w11, #1024 //offset val
      mov    w10, #0x03FF //max val
      DUP    V29.4S, w11 //offset vector
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      mov    x12, #2
      sub    x2, x2, #6//start h src[]
      sub    x13, x3, #7
      lsl    x13, x13, #1 // last jump h
      lsl    x1, x1, #1 //_dststride
      lsl    x14, x3, #1 //_srcstride
      lsl    x15, x14, #1 //2*_srcstride
      sub    x2, x2, x14
      sub    x2, x2, x15//start hv src[]
      mov    x16, x5//save height
      mov    x17, x2//save start src
      mov    x15, x0//save start dst
      mov    x14, x4//save start dst
      LD1    {V0.8B}, [x6]
      LD1    {V1.8B}, [x7]
      SADDW   V0.8H, V30.8H, V0.8B
      SADDW   V1.8H, V30.8H, V1.8B
      bi1_qpel_hv_filter_8lines x0, x2, V0, V1, x13, x12
      add x0, x15, #16
      add x2, x17, #16
      add x4, x14, #16
      mov x5, x16
      bi1_qpel_hv_filter_8lines x0, x2, V0, V1, x13, x12
      add x0, x15, #32
      add x2, x17, #32
      add x4, x14, #32
      mov x5, x16
      bi1_qpel_hv_filter_8lines x0, x2, V0, V1, x13, x12
      add x0, x15, #48
      add x2, x17, #48
      add x4, x14, #48
      mov x5, x16
      bi1_qpel_hv_filter_8lines x0, x2, V0, V1, x13, x12
      ret lr
endfunc

function put_vvc_bi1_qpel_hv_10_32_neon_4lines, export=1
      mov    x9, MAX_PB_SIZE_16
      mov    w11, #1024 //offset val
      mov    w10, #0x03FF //max val
      DUP    V29.4S, w11 //offset vector
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      mov    x12, #2
      sub    x2, x2, #6//start h src[]
      sub    x13, x3, #7
      lsl    x13, x13, #1 // last jump h
      lsl    x1, x1, #1 //_dststride
      lsl    x14, x3, #1 //_srcstride
      lsl    x15, x14, #1 //2*_srcstride
      sub    x2, x2, x14
      sub    x2, x2, x15//start hv src[]
      mov    x16, x5//save height
      mov    x17, x2//save start src
      mov    x15, x0//save start dst
      mov    x14, x4//save start dst
      LD1    {V0.8B}, [x6]
      LD1    {V1.8B}, [x7]
      SADDW   V0.8H, V30.8H, V0.8B
      SADDW   V1.8H, V30.8H, V1.8B
      bi1_qpel_hv_filter_4lines x0, x2, V0, V1, x13, x12
      add x0, x15, #16
      add x2, x17, #16
      add x4, x14, #16
      mov x5, x16
      bi1_qpel_hv_filter_4lines x0, x2, V0, V1, x13, x12
      add x0, x15, #32
      add x2, x17, #32
      add x4, x14, #32
      mov x5, x16
      bi1_qpel_hv_filter_4lines x0, x2, V0, V1, x13, x12
      add x0, x15, #48
      add x2, x17, #48
      add x4, x14, #48
      mov x5, x16
      bi1_qpel_hv_filter_4lines x0, x2, V0, V1, x13, x12
      ret lr
endfunc

.macro uni_epel_filter, out, coef=V0, l0=V1, l1=V2, l2=V3, l3=V4, shift=#6, s1=.4H, s2=.8H, s3=.H
      SMULL2 V11.4S,  \l0\s2, \coef\s3[0]
      SMLAL2 V11.4S,  \l1\s2, \coef\s3[1]
      SMLAL2 V11.4S, \l2\s2, \coef\s3[2]
      SMLAL2 V11.4S, \l3\s2, \coef\s3[3]
      SMULL V10.4S, \l0\s1, \coef\s3[0]
      SMLAL V10.4S, \l1\s1, \coef\s3[1]
      SMLAL V10.4S, \l2\s1, \coef\s3[2]
      SMLAL V10.4S, \l3\s1, \coef\s3[3]
      ADD   V10.4S, V10.4S, V29.4S
      ADD   V11.4S, V11.4S, V29.4S
      SSHR  V10.4S, V10.4S, \shift
      SSHR  V11.4S, V11.4S, \shift
      UZP1  \out, V10.8H, V11.8H
      CLIP  \out, V30.8H, V31.8H
.endm

.macro bi0_epel_filter, out, coef=V0, l0=V1, l1=V2, l2=V3, l3=V4, shift=#2,s1=.4H, s2=.8H, s3=.H
      SMULL2 V11.4S,  \l0\s2, \coef\s3[0]
      SMLAL2 V11.4S,  \l1\s2, \coef\s3[1]
      SMLAL2 V11.4S, \l2\s2, \coef\s3[2]
      SMLAL2 V11.4S, \l3\s2, \coef\s3[3]
      SMULL V10.4S, \l0\s1, \coef\s3[0]
      SMLAL V10.4S, \l1\s1, \coef\s3[1]
      SMLAL V10.4S, \l2\s1, \coef\s3[2]
      SMLAL V10.4S, \l3\s1, \coef\s3[3]
      SSHR  V10.4S, V10.4S, \shift
      SSHR  V11.4S, V11.4S, \shift
      UZP1  \out, V10.8H, V11.8H
.endm

.macro bi1_epel_filter, out, src1, coef=V0, l0=V1, l1=V2, l2=V3, l3=V4, shift=#2, s1=.4H, s2=.8H, s3=.H
      SMULL2 V11.4S,  \l0\s2, \coef\s3[0]
      SMLAL2 V11.4S,  \l1\s2, \coef\s3[1]
      SMLAL2 V11.4S, \l2\s2, \coef\s3[2]
      SMLAL2 V11.4S, \l3\s2, \coef\s3[3]
      SMULL V10.4S, \l0\s1, \coef\s3[0]
      SMLAL V10.4S, \l1\s1, \coef\s3[1]
      SMLAL V10.4S, \l2\s1, \coef\s3[2]
      SMLAL V10.4S, \l3\s1, \coef\s3[3]
      ADD   V10.4S, V10.4S, V29.4S
      ADD   V11.4S, V11.4S, V29.4S
      SSHR  V10.4S, V10.4S, \shift
      SSHR  V11.4S, V11.4S, \shift
      SADDW V10.4S, V10.4S,  \src1\s1
      SADDW2 V11.4S, V11.4S, \src1\s2
      SSHR  V10.4S, V10.4S, #5
      SSHR  V11.4S, V11.4S, #5
      UZP1  \out, V10.8H, V11.8H
      CLIP  \out, V30.8H, V31.8H
.endm

.macro uni_epel_h_filter, out, last_jump
      LD1   {V1.8H}, [x2], x12
      LD1   {V2.8H}, [x2], x12
      LD1   {V3.8H}, [x2], x12
      LD1   {V4.8H}, [x2], \last_jump
      uni_epel_filter \out
.endm

.macro bi0_epel_h_filter, out, last_jump, coef=V0, jump=x12, in=x1, l0=V1, l1=V2, l2=V3, l3=V4, s1=.4H, s2=.8H, s3=.H
      LD1   {\l0\s2}, [\in], \jump
      LD1   {\l1\s2}, [\in], \jump
      LD1   {\l2\s2}, [\in], \jump
      LD1   {\l3\s2}, [\in], \last_jump
      bi0_epel_filter \out, \coef, \l0, \l1, \l2, \l3
.endm

function put_vvc_uni_epel_h_10_8_neon, export=1
      mov    x12, #2
      mov    w11, #32 //offset val
      mov    w10, #0x03FF //max val
      sub    x2, x2, #2
      sub    x3, x3, #3
      lsl    x1, x1, #1 //_dststride
      lsl    x3, x3, #1 //_srcstride
      DUP    V29.4S, w11 //offset vector
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      LD1    {V0.8B}, [x5]
      SADDW   V0.8H, V30.8H, V0.8B
0:
      subs x4, x4, #2
      uni_epel_h_filter V14.8H, x3
      uni_epel_h_filter V15.8H, x3
      ST1   {V14.8H}, [x0], x1
      ST1   {V15.8H}, [x0], x1
      b.ne 0b
      ret lr
endfunc

function put_vvc_uni_epel_h_10_16_neon, export=1
      mov    x12, #2
      mov    x13, #10
      mov    w11, #32 //offset val
      mov    w10, #0x03FF //max val
      sub    x2, x2, #2
      sub    x3, x3, #11
      lsl    x1, x1, #1 //_dststride
      lsl    x3, x3, #1 //_srcstride
      DUP    V29.4S, w11 //offset vector
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      LD1    {V0.8B}, [x5]
      SADDW   V0.8H, V30.8H, V0.8B
0:
      subs x4, x4, #1
      uni_epel_h_filter V14.8H, x13
      uni_epel_h_filter V15.8H, x3
      ST1   {V14.8H-V15.8H}, [x0], x1
      b.ne 0b
      ret lr
endfunc

function put_vvc_uni_epel_h_10_32_neon, export=1
      mov    x12, #2
      mov    x13, #10
      mov    w11, #32 //offset val
      mov    w10, #0x03FF //max val
      sub    x2, x2, #2
      sub    x3, x3, #27
      lsl    x1, x1, #1 //_dststride
      lsl    x3, x3, #1 //_srcstride
      DUP    V29.4S, w11 //offset vector
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      LD1    {V0.8B}, [x5]
      SADDW   V0.8H, V30.8H, V0.8B
0:
      subs x4, x4, #1
      uni_epel_h_filter V14.8H, x13
      uni_epel_h_filter V15.8H, x13
      uni_epel_h_filter V16.8H, x13
      uni_epel_h_filter V17.8H, x3
      ST1   {V14.8H-V17.8H}, [x0], x1
      b.ne 0b
      ret lr
endfunc

.macro uni_epel_v_filter_4lines, out, in
      LD1    {V1.8H}, [\in], x3
      LD1    {V2.8H}, [\in], x3
      LD1    {V3.8H}, [\in], x3
0:              subs x4, x4, #4
      LD1    {V4.8H}, [\in], x3
      uni_epel_filter V17.8H, V0, V1, V2, V3, V4
      LD1    {V1.8H}, [\in], x3
      uni_epel_filter V18.8H, V0, V2, V3, V4, V1
      LD1    {V2.8H}, [\in], x3
      uni_epel_filter V19.8H, V0, V3, V4, V1, V2
      LD1    {V3.8H}, [\in], x3
      uni_epel_filter V20.8H, V0, V4, V1, V2, V3
      ST1   {V17.8H}, [\out], x1
      ST1   {V18.8H}, [\out], x1
      ST1   {V19.8H}, [\out], x1
      ST1   {V20.8H}, [\out], x1
      b.ne 0b
.endm

.macro uni_epel_v_filter_2lines, out, in
      LD1    {V1.8H}, [\in], x3
      LD1    {V2.8H}, [\in], x3
      LD1    {V3.8H}, [\in], x3
0:              subs x4, x4, #2
      LD1    {V4.8H}, [\in], x3
      uni_epel_filter V17.8H, V0, V1, V2, V3, V4
      LD1    {V1.8H}, [\in], x3
      uni_epel_filter V18.8H, V0, V2, V3, V4, V1
      ST1   {V17.8H}, [\out], x1
      ST1   {V18.8H}, [\out], x1
      b.ne 0b
.endm

function put_vvc_uni_epel_v_10_8_neon_4lines, export=1
      mov    x12, #2
      mov    w11, #32 //offset val
      mov    w10, #0x03FF //max val
      lsl    x1, x1, #1 //_dststride
      lsl    x3, x3, #1 //_srcstride
      sub    x2, x2, x3
      DUP    V29.4S, w11 //offset vector
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      LD1    {V0.8B}, [x5]
      SADDW   V0.8H, V30.8H, V0.8B
      uni_epel_v_filter_4lines x0, x2
      ret lr
endfunc

function put_vvc_uni_epel_v_10_8_neon_2lines, export=1
      mov    x12, #2
      mov    w11, #32 //offset val
      mov    w10, #0x03FF //max val
      lsl    x1, x1, #1 //_dststride
      lsl    x3, x3, #1 //_srcstride
      sub    x2, x2, x3
      DUP    V29.4S, w11 //offset vector
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      LD1    {V0.8B}, [x5]
      SADDW   V0.8H, V30.8H, V0.8B
      uni_epel_v_filter_2lines x0, x2
      ret lr
endfunc

function put_vvc_uni_epel_v_10_16_neon_4lines, export=1
      mov    x12, #2
      mov    w11, #32 //offset val
      mov    w10, #0x03FF //max val
      lsl    x1, x1, #1 //_dststride
      lsl    x3, x3, #1 //_srcstride
      sub    x2, x2, x3
      mov    x9, x2  //save x2[-_srcstride]
      mov    x8, x0  //save x0[0]
      mov    x14, x4  //save height
      DUP    V29.4S, w11 //offset vector
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      LD1    {V0.8B}, [x5]
      SADDW   V0.8H, V30.8H, V0.8B
      uni_epel_v_filter_4lines x0, x2
      add x0, x8, #16
      add x2, x9, #16
      mov x4, x14
      uni_epel_v_filter_4lines x0, x2
      ret lr
endfunc

function put_vvc_uni_epel_v_10_16_neon_2lines, export=1
      mov    x12, #2
      mov    w11, #32 //offset val
      mov    w10, #0x03FF //max val
      lsl    x1, x1, #1 //_dststride
      lsl    x3, x3, #1 //_srcstride
      sub    x2, x2, x3
      mov    x9, x2  //save x2[-_srcstride]
      mov    x8, x0  //save x0[0]
      mov    x14, x4  //save height
      DUP    V29.4S, w11 //offset vector
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      LD1    {V0.8B}, [x5]
      SADDW   V0.8H, V30.8H, V0.8B
      uni_epel_v_filter_2lines x0, x2
      add x0, x8, #16
      add x2, x9, #16
      mov x4, x14
      uni_epel_v_filter_2lines x0, x2
      ret lr
endfunc

function put_vvc_uni_epel_v_10_32_neon_4lines, export=1
      mov    x12, #2
      mov    w11, #32 //offset val
      mov    w10, #0x03FF //max val
      lsl    x1, x1, #1 //_dststride
      lsl    x3, x3, #1 //_srcstride
      sub    x2, x2, x3
      mov    x9, x2  //save x2[-_srcstride]
      mov    x8, x0  //save x0[0]
      mov    x14, x4  //save height
      DUP    V29.4S, w11 //offset vector
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      LD1    {V0.8B}, [x5]
      SADDW   V0.8H, V30.8H, V0.8B
      uni_epel_v_filter_4lines x0, x2
      add x0, x8, #16
      add x2, x9, #16
      mov x4, x14
      uni_epel_v_filter_4lines x0, x2
      add x0, x8, #32
      add x2, x9, #32
      mov x4, x14
      uni_epel_v_filter_4lines x0, x2
      add x0, x8, #48
      add x2, x9, #48
      mov x4, x14
      uni_epel_v_filter_4lines x0, x2
      ret lr
endfunc

function put_vvc_uni_epel_v_10_32_neon_2lines, export=1
      mov    x12, #2
      mov    w11, #32 //offset val
      mov    w10, #0x03FF //max val
      lsl    x1, x1, #1 //_dststride
      lsl    x3, x3, #1 //_srcstride
      sub    x2, x2, x3
      mov    x9, x2  //save x2[-_srcstride]
      mov    x8, x0  //save x0[0]
      mov    x14, x4  //save height
      DUP    V29.4S, w11 //offset vector
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      LD1    {V0.8B}, [x5]
      SADDW   V0.8H, V30.8H, V0.8B
      uni_epel_v_filter_2lines x0, x2
      add x0, x8, #16
      add x2, x9, #16
      mov x4, x14
      uni_epel_v_filter_2lines x0, x2
      add x0, x8, #32
      add x2, x9, #32
      mov x4, x14
      uni_epel_v_filter_2lines x0, x2
      add x0, x8, #48
      add x2, x9, #48
      mov x4, x14
      uni_epel_v_filter_2lines x0, x2
      ret lr
endfunc

.macro uni_epel_hv_filter_4lines, out, in, coefx=V0, coefy=V1, last_jump=x3, jump=x12
      bi0_epel_h_filter V2.8H, \last_jump, \coefx, \jump, \in, V6, V7, V8, V9
      bi0_epel_h_filter V3.8H, \last_jump, \coefx, \jump, \in, V6, V7, V8, V9
      bi0_epel_h_filter V4.8H, \last_jump, \coefx, \jump, \in, V6, V7, V8, V9
0:              subs x4, x4, #4
      bi0_epel_h_filter V5.8H, \last_jump, \coefx, \jump, \in, V6, V7, V8, V9
      uni_epel_filter V20.8H, \coefy, V2, V3, V4, V5, shift=#10
      bi0_epel_h_filter V2.8H, \last_jump, \coefx, \jump, \in, V6, V7, V8, V9
      uni_epel_filter V21.8H, \coefy, V3, V4, V5, V2, shift=#10
      bi0_epel_h_filter V3.8H, \last_jump, \coefx, \jump, \in, V6, V7, V8, V9
      uni_epel_filter V22.8H, \coefy, V4, V5, V2, V3, shift=#10
      bi0_epel_h_filter V4.8H, \last_jump, \coefx, \jump, \in, V6, V7, V8, V9
      uni_epel_filter V23.8H, \coefy, V5, V2, V3, V4, shift=#10
      ST1   {V20.8H}, [\out], x1
      ST1   {V21.8H}, [\out], x1
      ST1   {V22.8H}, [\out], x1
      ST1   {V23.8H}, [\out], x1
      b.ne 0b
.endm

.macro uni_epel_hv_filter_2lines, out, in, coefx=V0, coefy=V1, last_jump=x3, jump=x12
      bi0_epel_h_filter V2.8H, \last_jump, \coefx, \jump, \in, V6, V7, V8, V9
      bi0_epel_h_filter V3.8H, \last_jump, \coefx, \jump, \in, V6, V7, V8, V9
      bi0_epel_h_filter V4.8H, \last_jump, \coefx, \jump, \in, V6, V7, V8, V9
0:              subs x4, x4, #2
      bi0_epel_h_filter V5.8H, \last_jump, \coefx, \jump, \in, V6, V7, V8, V9
      uni_epel_filter V20.8H, \coefy, V2, V3, V4, V5, shift=#10
      bi0_epel_h_filter V2.8H, \last_jump, \coefx, \jump, \in, V6, V7, V8, V9
      uni_epel_filter V21.8H, \coefy, V3, V4, V5, V2, shift=#10
      ST1   {V20.8H}, [\out], x1
      ST1   {V21.8H}, [\out], x1
      b.ne 0b
.endm

function put_vvc_uni_epel_hv_10_8_neon_4lines, export=1
      mov    w11, #512 //offset val
      mov    w10, #0x03FF //max val
      DUP    V29.4S, w11 //offset vector
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      mov    x12, #2
      sub    x2, x2, #2//start h src[]
      sub    x13, x3, #3
      lsl    x13, x13, #1 // last jump h
      lsl    x1, x1, #1 //_dststride
      lsl    x14, x3, #1 //_srcstride
      sub    x2, x2, x14//start hv src[]
      LD1    {V0.8B}, [x5]
      LD1    {V1.8B}, [x6]
      SADDW   V0.8H, V30.8H, V0.8B
      SADDW   V1.8H, V30.8H, V1.8B
      uni_epel_hv_filter_4lines x0, x2, V0, V1, x13, x12
      ret lr
endfunc

function put_vvc_uni_epel_hv_10_8_neon_2lines, export=1
      mov    w11, #512 //offset val
      mov    w10, #0x03FF //max val
      DUP    V29.4S, w11 //offset vector
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      mov    x12, #2
      sub    x2, x2, #2//start h src[]
      sub    x13, x3, #3
      lsl    x13, x13, #1 // last jump h
      lsl    x1, x1, #1 //_dststride
      lsl    x14, x3, #1 //_srcstride
      sub    x2, x2, x14//start hv src[]
      LD1    {V0.8B}, [x5]
      LD1    {V1.8B}, [x6]
      SADDW   V0.8H, V30.8H, V0.8B
      SADDW   V1.8H, V30.8H, V1.8B
      uni_epel_hv_filter_2lines x0, x2, V0, V1, x13, x12
      ret lr
endfunc

function put_vvc_uni_epel_hv_10_16_neon_4lines, export=1
      mov    w11, #512 //offset val
      mov    w10, #0x03FF //max val
      DUP    V29.4S, w11 //offset vector
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      mov    x12, #2
      sub    x2, x2, #2//start h src[]
      sub    x13, x3, #3
      lsl    x13, x13, #1 // last jump h
      lsl    x1, x1, #1 //_dststride
      lsl    x14, x3, #1 //_srcstride
      sub    x2, x2, x14//start hv src[]
      mov    x16, x4//save height
      mov    x17, x2//save start src
      mov    x15, x0//save start dst
      LD1    {V0.8B}, [x5]
      LD1    {V1.8B}, [x6]
      SADDW   V0.8H, V30.8H, V0.8B
      SADDW   V1.8H, V30.8H, V1.8B
      uni_epel_hv_filter_4lines x0, x2, V0, V1, x13, x12
      add x0, x15, #16
      add x2, x17, #16
      mov x4, x16
      uni_epel_hv_filter_4lines x0, x2, V0, V1, x13, x12
      ret lr
endfunc

function put_vvc_uni_epel_hv_10_16_neon_2lines, export=1
      mov    w11, #512 //offset val
      mov    w10, #0x03FF //max val
      DUP    V29.4S, w11 //offset vector
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      mov    x12, #2
      sub    x2, x2, #2//start h src[]
      sub    x13, x3, #3
      lsl    x13, x13, #1 // last jump h
      lsl    x1, x1, #1 //_dststride
      lsl    x14, x3, #1 //_srcstride
      sub    x2, x2, x14//start hv src[]
      mov    x16, x4//save height
      mov    x17, x2//save start src
      mov    x15, x0//save start dst
      LD1    {V0.8B}, [x5]
      LD1    {V1.8B}, [x6]
      SADDW   V0.8H, V30.8H, V0.8B
      SADDW   V1.8H, V30.8H, V1.8B
      uni_epel_hv_filter_2lines x0, x2, V0, V1, x13, x12
      add x0, x15, #16
      add x2, x17, #16
      mov x4, x16
      uni_epel_hv_filter_2lines x0, x2, V0, V1, x13, x12
      ret lr
endfunc

function put_vvc_uni_epel_hv_10_32_neon_4lines, export=1
      mov    w11, #512 //offset val
      mov    w10, #0x03FF //max val
      DUP    V29.4S, w11 //offset vector
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      mov    x12, #2
      sub    x2, x2, #2//start h src[]
      sub    x13, x3, #3
      lsl    x13, x13, #1 // last jump h
      lsl    x1, x1, #1 //_dststride
      lsl    x14, x3, #1 //_srcstride
      sub    x2, x2, x14//start hv src[]
      mov    x16, x4//save height
      mov    x17, x2//save start src
      mov    x15, x0//save start dst
      LD1    {V0.8B}, [x5]
      LD1    {V1.8B}, [x6]
      SADDW   V0.8H, V30.8H, V0.8B
      SADDW   V1.8H, V30.8H, V1.8B
      uni_epel_hv_filter_4lines x0, x2, V0, V1, x13, x12
      add x0, x15, #16
      add x2, x17, #16
      mov x4, x16
      uni_epel_hv_filter_4lines x0, x2, V0, V1, x13, x12
      add x0, x15, #32
      add x2, x17, #32
      mov x4, x16
      uni_epel_hv_filter_4lines x0, x2, V0, V1, x13, x12
      add x0, x15, #48
      add x2, x17, #48
      mov x4, x16
      uni_epel_hv_filter_4lines x0, x2, V0, V1, x13, x12
      ret lr
endfunc

function put_vvc_uni_epel_hv_10_32_neon_2lines, export=1
      mov    w11, #512 //offset val
      mov    w10, #0x03FF //max val
      DUP    V29.4S, w11 //offset vector
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      mov    x12, #2
      sub    x2, x2, #2//start h src[]
      sub    x13, x3, #3
      lsl    x13, x13, #1 // last jump h
      lsl    x1, x1, #1 //_dststride
      lsl    x14, x3, #1 //_srcstride
      sub    x2, x2, x14//start hv src[]
      mov    x16, x4//save height
      mov    x17, x2//save start src
      mov    x15, x0//save start dst
      LD1    {V0.8B}, [x5]
      LD1    {V1.8B}, [x6]
      SADDW   V0.8H, V30.8H, V0.8B
      SADDW   V1.8H, V30.8H, V1.8B
      uni_epel_hv_filter_2lines x0, x2, V0, V1, x13, x12
      add x0, x15, #16
      add x2, x17, #16
      mov x4, x16
      uni_epel_hv_filter_2lines x0, x2, V0, V1, x13, x12
      add x0, x15, #32
      add x2, x17, #32
      mov x4, x16
      uni_epel_hv_filter_2lines x0, x2, V0, V1, x13, x12
      add x0, x15, #48
      add x2, x17, #48
      mov x4, x16
      uni_epel_hv_filter_2lines x0, x2, V0, V1, x13, x12
      ret lr
endfunc

function put_vvc_bi0_epel_h_10_8_neon, export=1
      mov    x12, #2
      mov    x11, MAX_PB_SIZE_16 //offset val
      mov    w10, #0x03FF //max val
      sub    x1, x1, #2
      sub    x2, x2, #3
      lsl    x2, x2, #1 //_srcstride
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      LD1    {V0.8B}, [x4]
      SADDW   V0.8H, V30.8H, V0.8B
0:
      subs x3, x3, #2
      bi0_epel_h_filter V14.8H, x2
      bi0_epel_h_filter V15.8H, x2
      ST1   {V14.8H}, [x0], x11
      ST1   {V15.8H}, [x0], x11
      b.ne 0b
      ret lr
endfunc

function put_vvc_bi0_epel_h_10_16_neon, export=1
      mov    x12, #2
      mov    x13, #10
      mov    x11, MAX_PB_SIZE_16 //offset val
      mov    w10, #0x03FF //max val
      sub    x1, x1, #2
      sub    x2, x2, #11
      lsl    x2, x2, #1 //_srcstride
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      LD1    {V0.8B}, [x4]
      SADDW   V0.8H, V30.8H, V0.8B
0:
      subs x3, x3, #1
      bi0_epel_h_filter V14.8H, x13
      bi0_epel_h_filter V15.8H, x2
      ST1   {V14.8H-V15.8H}, [x0], x11
      b.ne 0b
      ret lr
endfunc

function put_vvc_bi0_epel_h_10_32_neon, export=1
      mov    x12, #2
      mov    x13, #10
      mov    x11, MAX_PB_SIZE_16 //offset val
      mov    w10, #0x03FF //max val
      sub    x1, x1, #2
      sub    x2, x2, #27
      lsl    x2, x2, #1 //_srcstride
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      LD1    {V0.8B}, [x4]
      SADDW   V0.8H, V30.8H, V0.8B
0:
      subs x3, x3, #1
      bi0_epel_h_filter V14.8H, x13
      bi0_epel_h_filter V15.8H, x13
      bi0_epel_h_filter V16.8H, x13
      bi0_epel_h_filter V17.8H, x2
      ST1   {V14.8H-V17.8H}, [x0], x11
      b.ne 0b
      ret lr
endfunc

.macro bi0_epel_v_filter_4lines, out, in
      LD1    {V1.8H}, [\in], x2
      LD1    {V2.8H}, [\in], x2
      LD1    {V3.8H}, [\in], x2

0:              subs x3, x3, #4
      LD1    {V4.8H}, [\in], x2
      bi0_epel_filter V17.8H, V0, V1, V2, V3, V4
      LD1    {V1.8H}, [\in], x2
      bi0_epel_filter V18.8H, V0, V2, V3, V4, V1
      LD1    {V2.8H}, [\in], x2
      bi0_epel_filter V19.8H, V0, V3, V4, V1, V2
      LD1    {V3.8H}, [\in], x2
      bi0_epel_filter V20.8H, V0, V4, V1, V2, V3
      ST1   {V17.8H}, [\out], x11
      ST1   {V18.8H}, [\out], x11
      ST1   {V19.8H}, [\out], x11
      ST1   {V20.8H}, [\out], x11
      b.ne 0b
.endm

.macro bi0_epel_v_filter_2lines, out, in
      LD1    {V1.8H}, [\in], x2
      LD1    {V2.8H}, [\in], x2
      LD1    {V3.8H}, [\in], x2

0:              subs x3, x3, #2
      LD1    {V4.8H}, [\in], x2
      bi0_epel_filter V17.8H, V0, V1, V2, V3, V4
      LD1    {V1.8H}, [\in], x2
      bi0_epel_filter V18.8H, V0, V2, V3, V4, V1
      ST1   {V17.8H}, [\out], x11
      ST1   {V18.8H}, [\out], x11
      b.ne 0b
.endm

function put_vvc_bi0_epel_v_10_8_neon_4lines, export=1
      mov    x11, MAX_PB_SIZE_16 //offset val
      mov    w10, #0x03FF //max val
      lsl    x2, x2, #1 //_srcstride
      sub    x1, x1, x2
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      LD1    {V0.8B}, [x4]
      SADDW   V0.8H, V30.8H, V0.8B
      bi0_epel_v_filter_4lines x0, x1
      ret lr
endfunc

function put_vvc_bi0_epel_v_10_8_neon_2lines, export=1
      mov    x11, MAX_PB_SIZE_16 //offset val
      mov    w10, #0x03FF //max val
      lsl    x2, x2, #1 //_srcstride
      sub    x1, x1, x2
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      LD1    {V0.8B}, [x4]
      SADDW   V0.8H, V30.8H, V0.8B
      bi0_epel_v_filter_2lines x0, x1
      ret lr
endfunc

function put_vvc_bi0_epel_v_10_16_neon_4lines, export=1
      mov    x11, MAX_PB_SIZE_16 //offset val
      mov    w10, #0x03FF //max val
      lsl    x2, x2, #1 //_srcstride
      sub    x1, x1, x2
      mov    x9, x1  //save x1[-3*_srcstride]
      mov    x8, x0  //save x0[0]
      mov    x14, x3  //save height
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      LD1    {V0.8B}, [x4]
      SADDW   V0.8H, V30.8H, V0.8B
      bi0_epel_v_filter_4lines x0, x1
      add x0, x8, #16
      add x1, x9, #16
      mov x3, x14
      bi0_epel_v_filter_4lines x0, x1
      ret lr
endfunc

function put_vvc_bi0_epel_v_10_16_neon_2lines, export=1
      mov    x11, MAX_PB_SIZE_16 //offset val
      mov    w10, #0x03FF //max val
      lsl    x2, x2, #1 //_srcstride
      sub    x1, x1, x2
      mov    x9, x1  //save x1[-3*_srcstride]
      mov    x8, x0  //save x0[0]
      mov    x14, x3  //save height
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      LD1    {V0.8B}, [x4]
      SADDW   V0.8H, V30.8H, V0.8B
      bi0_epel_v_filter_2lines x0, x1
      add x0, x8, #16
      add x1, x9, #16
      mov x3, x14
      bi0_epel_v_filter_2lines x0, x1
      ret lr
endfunc

function put_vvc_bi0_epel_v_10_32_neon_4lines, export=1
      mov    x11, MAX_PB_SIZE_16 //offset val
      mov    w10, #0x03FF //max val
      lsl    x2, x2, #1 //_srcstride
      sub    x1, x1, x2
      mov    x9, x1  //save x1[-3*_srcstride]
      mov    x8, x0  //save x0[0]
      mov    x14, x3  //save height
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      LD1    {V0.8B}, [x4]
      SADDW   V0.8H, V30.8H, V0.8B
      bi0_epel_v_filter_4lines x0, x1
      add x0, x8, #16
      add x1, x9, #16
      mov x3, x14
      bi0_epel_v_filter_4lines x0, x1
      add x0, x8, #32
      add x1, x9, #32
      mov x3, x14
      bi0_epel_v_filter_4lines x0, x1
      add x0, x8, #48
      add x1, x9, #48
      mov x3, x14
      bi0_epel_v_filter_4lines x0, x1
      ret lr
endfunc

function put_vvc_bi0_epel_v_10_32_neon_2lines, export=1
      mov    x11, MAX_PB_SIZE_16 //offset val
      mov    w10, #0x03FF //max val
      lsl    x2, x2, #1 //_srcstride
      sub    x1, x1, x2
      mov    x9, x1  //save x1[-3*_srcstride]
      mov    x8, x0  //save x0[0]
      mov    x14, x3  //save height
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      LD1    {V0.8B}, [x4]
      SADDW   V0.8H, V30.8H, V0.8B
      bi0_epel_v_filter_2lines x0, x1
      add x0, x8, #16
      add x1, x9, #16
      mov x3, x14
      bi0_epel_v_filter_2lines x0, x1
      add x0, x8, #32
      add x1, x9, #32
      mov x3, x14
      bi0_epel_v_filter_2lines x0, x1
      add x0, x8, #48
      add x1, x9, #48
      mov x3, x14
      bi0_epel_v_filter_2lines x0, x1
      ret lr
endfunc

.macro bi0_epel_hv_filter_4lines, out, in, coefx=V0, coefy=V1, last_jump=x3, jump=x12
      bi0_epel_h_filter V2.8H, \last_jump, \coefx, \jump, \in, V6, V7, V8, V9
      bi0_epel_h_filter V3.8H, \last_jump, \coefx, \jump, \in, V6, V7, V8, V9
      bi0_epel_h_filter V4.8H, \last_jump, \coefx, \jump, \in, V6, V7, V8, V9
0:              subs x3, x3, #4
      bi0_epel_h_filter V5.8H, \last_jump, \coefx, \jump, \in, V6, V7, V8, V9
      bi0_epel_filter V20.8H, \coefy, V2, V3, V4, V5, shift=#6
      bi0_epel_h_filter V2.8H, \last_jump, \coefx, \jump, \in, V6, V7, V8, V9
      bi0_epel_filter V21.8H, \coefy, V3, V4, V5, V2, shift=#6
      bi0_epel_h_filter V3.8H, \last_jump, \coefx, \jump, \in, V6, V7, V8, V9
      bi0_epel_filter V22.8H, \coefy, V4, V5, V2, V3, shift=#6
      bi0_epel_h_filter V4.8H, \last_jump, \coefx, \jump, \in, V6, V7, V8, V9
      bi0_epel_filter V23.8H, \coefy, V5, V2, V3, V4, shift=#6
      ST1   {V20.8H}, [\out], x11
      ST1   {V21.8H}, [\out], x11
      ST1   {V22.8H}, [\out], x11
      ST1   {V23.8H}, [\out], x11
      b.ne 0b
.endm

.macro bi0_epel_hv_filter_2lines, out, in, coefx=V0, coefy=V1, last_jump=x3, jump=x12
      bi0_epel_h_filter V2.8H, \last_jump, \coefx, \jump, \in, V6, V7, V8, V9
      bi0_epel_h_filter V3.8H, \last_jump, \coefx, \jump, \in, V6, V7, V8, V9
      bi0_epel_h_filter V4.8H, \last_jump, \coefx, \jump, \in, V6, V7, V8, V9
0:              subs x3, x3, #2
      bi0_epel_h_filter V5.8H, \last_jump, \coefx, \jump, \in, V6, V7, V8, V9
      bi0_epel_filter V20.8H, \coefy, V2, V3, V4, V5, shift=#6
      bi0_epel_h_filter V2.8H, \last_jump, \coefx, \jump, \in, V6, V7, V8, V9
      bi0_epel_filter V21.8H, \coefy, V3, V4, V5, V2, shift=#6
      ST1   {V20.8H}, [\out], x11
      ST1   {V21.8H}, [\out], x11
      b.ne 0b
.endm

function put_vvc_bi0_epel_hv_10_8_neon_4lines, export=1
      mov    x11, MAX_PB_SIZE_16
      mov    w10, #0x03FF //max val
      DUP    V29.4S, w11 //offset vector
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      mov    x12, #2
      sub    x1, x1, #2//start h src[]
      sub    x13, x2, #3
      lsl    x13, x13, #1 // last jump h
      lsl    x14, x2, #1 //_srcstride
      sub    x1, x1, x14
      LD1    {V0.8B}, [x4]
      LD1    {V1.8B}, [x5]
      SADDW   V0.8H, V30.8H, V0.8B
      SADDW   V1.8H, V30.8H, V1.8B
      bi0_epel_hv_filter_4lines x0, x1, V0, V1, x13, x12
      ret lr
endfunc

function put_vvc_bi0_epel_hv_10_8_neon_2lines, export=1
      mov    x11, MAX_PB_SIZE_16
      mov    w10, #0x03FF //max val
      DUP    V29.4S, w11 //offset vector
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      mov    x12, #2
      sub    x1, x1, #2//start h src[]
      sub    x13, x2, #3
      lsl    x13, x13, #1 // last jump h
      lsl    x14, x2, #1 //_srcstride
      sub    x1, x1, x14
      LD1    {V0.8B}, [x4]
      LD1    {V1.8B}, [x5]
      SADDW   V0.8H, V30.8H, V0.8B
      SADDW   V1.8H, V30.8H, V1.8B
      bi0_epel_hv_filter_2lines x0, x1, V0, V1, x13, x12
      ret lr
endfunc

function put_vvc_bi0_epel_hv_10_16_neon_4lines, export=1
      mov    x11, MAX_PB_SIZE_16
      mov    w10, #0x03FF //max val
      DUP    V29.4S, w11 //offset vector
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      mov    x12, #2
      sub    x1, x1, #2//start h src[]
      sub    x13, x2, #3
      lsl    x13, x13, #1 // last jump h
      lsl    x14, x2, #1 //_srcstride
      sub    x1, x1, x14
      mov    x16, x3//save height
      mov    x17, x1//save start src
      mov    x15, x0//save start dst
      LD1    {V0.8B}, [x4]
      LD1    {V1.8B}, [x5]
      SADDW   V0.8H, V30.8H, V0.8B
      SADDW   V1.8H, V30.8H, V1.8B
      bi0_epel_hv_filter_4lines x0, x1, V0, V1, x13, x12
      add x0, x15, #16
      add x1, x17, #16
      mov x3, x16
      bi0_epel_hv_filter_4lines x0, x1, V0, V1, x13, x12
      ret lr
endfunc

function put_vvc_bi0_epel_hv_10_16_neon_2lines, export=1
      mov    x11, MAX_PB_SIZE_16
      mov    w10, #0x03FF //max val
      DUP    V29.4S, w11 //offset vector
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      mov    x12, #2
      sub    x1, x1, #2//start h src[]
      sub    x13, x2, #3
      lsl    x13, x13, #1 // last jump h
      lsl    x14, x2, #1 //_srcstride
      sub    x1, x1, x14
      mov    x16, x3//save height
      mov    x17, x1//save start src
      mov    x15, x0//save start dst
      LD1    {V0.8B}, [x4]
      LD1    {V1.8B}, [x5]
      SADDW   V0.8H, V30.8H, V0.8B
      SADDW   V1.8H, V30.8H, V1.8B
      bi0_epel_hv_filter_2lines x0, x1, V0, V1, x13, x12
      add x0, x15, #16
      add x1, x17, #16
      mov x3, x16
      bi0_epel_hv_filter_2lines x0, x1, V0, V1, x13, x12
      ret lr
endfunc

function put_vvc_bi0_epel_hv_10_32_neon_4lines, export=1
      mov    x11, MAX_PB_SIZE_16
      mov    w10, #0x03FF //max val
      DUP    V29.4S, w11 //offset vector
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      mov    x12, #2
      sub    x1, x1, #2//start h src[]
      sub    x13, x2, #3
      lsl    x13, x13, #1 // last jump h
      lsl    x14, x2, #1 //_srcstride
      sub    x1, x1, x14
      mov    x16, x3//save height
      mov    x17, x1//save start src
      mov    x15, x0//save start dst
      LD1    {V0.8B}, [x4]
      LD1    {V1.8B}, [x5]
      SADDW   V0.8H, V30.8H, V0.8B
      SADDW   V1.8H, V30.8H, V1.8B
      bi0_epel_hv_filter_4lines x0, x1, V0, V1, x13, x12
      add x0, x15, #16
      add x1, x17, #16
      mov x3, x16
      bi0_epel_hv_filter_4lines x0, x1, V0, V1, x13, x12
      add x0, x15, #32
      add x1, x17, #32
      mov x3, x16
      bi0_epel_hv_filter_4lines x0, x1, V0, V1, x13, x12
      add x0, x15, #48
      add x1, x17, #48
      mov x3, x16
      bi0_epel_hv_filter_4lines x0, x1, V0, V1, x13, x12
      ret lr
endfunc

function put_vvc_bi0_epel_hv_10_32_neon_2lines, export=1
      mov    x11, MAX_PB_SIZE_16
      mov    w10, #0x03FF //max val
      DUP    V29.4S, w11 //offset vector
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      mov    x12, #2
      sub    x1, x1, #2//start h src[]
      sub    x13, x2, #3
      lsl    x13, x13, #1 // last jump h
      lsl    x14, x2, #1 //_srcstride
      sub    x1, x1, x14
      mov    x16, x3//save height
      mov    x17, x1//save start src
      mov    x15, x0//save start dst
      LD1    {V0.8B}, [x4]
      LD1    {V1.8B}, [x5]
      SADDW   V0.8H, V30.8H, V0.8B
      SADDW   V1.8H, V30.8H, V1.8B
      bi0_epel_hv_filter_2lines x0, x1, V0, V1, x13, x12
      add x0, x15, #16
      add x1, x17, #16
      mov x3, x16
      bi0_epel_hv_filter_2lines x0, x1, V0, V1, x13, x12
      add x0, x15, #32
      add x1, x17, #32
      mov x3, x16
      bi0_epel_hv_filter_2lines x0, x1, V0, V1, x13, x12
      add x0, x15, #48
      add x1, x17, #48
      mov x3, x16
      bi0_epel_hv_filter_2lines x0, x1, V0, V1, x13, x12
      ret lr
endfunc

.macro bi1_epel_h_filter, out, last_jump, src1, suffix1=.4H, suffix2=.8H
      LD1   {V1.8H}, [x2], x12
      LD1   {V2.8H}, [x2], x12
      LD1   {V3.8H}, [x2], x12
      LD1   {V4.8H}, [x2], \last_jump
      bi1_epel_filter \out, \src1, V0
.endm

function put_vvc_bi1_epel_h_10_8_neon, export=1
      mov    x12, #2
      mov    w11, #64 //offset val
      mov    w10, #0x03FF //max val
      mov    x9, MAX_PB_SIZE_16
      sub    x2, x2, #2
      sub    x3, x3, #3
      lsl    x1, x1, #1 //_dststride
      lsl    x3, x3, #1 //_srcstride
      DUP    V29.4S, w11 //offset vector
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      LD1    {V0.8B}, [x6]
      SADDW   V0.8H, V30.8H, V0.8B
0:
      subs x5, x5, #2
      LD1   {V27.8H}, [x4], x9
      LD1   {V28.8H}, [x4], x9
      bi1_epel_h_filter V14.8H, x3, V27
      bi1_epel_h_filter V15.8H, x3, V28
      ST1   {V14.8H}, [x0], x1
      ST1   {V15.8H}, [x0], x1
      b.ne 0b
      ret lr
endfunc

function put_vvc_bi1_epel_h_10_16_neon, export=1
      mov    x12, #2
      mov    x13, #10
      mov    w11, #64 //offset val
      mov    w10, #0x03FF //max val
      mov    x9, MAX_PB_SIZE_16
      sub    x2, x2, #2
      sub    x3, x3, #11
      lsl    x1, x1, #1 //_dststride
      lsl    x3, x3, #1 //_srcstride
      DUP    V29.4S, w11 //offset vector
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      LD1    {V0.8B}, [x6]
      SADDW   V0.8H, V30.8H, V0.8B
0:
      subs x5, x5, #1
      LD1   {V27.8H-V28.8H}, [x4], x9
      bi1_epel_h_filter V14.8H, x13, V27
      bi1_epel_h_filter V15.8H, x3, V28
      ST1   {V14.8H-V15.8H}, [x0], x1
      b.ne 0b
      ret lr
endfunc

function put_vvc_bi1_epel_h_10_32_neon, export=1
      mov    x12, #2
      mov    x13, #10
      mov    w11, #64 //offset val
      mov    w10, #0x03FF //max val
      mov    x9, MAX_PB_SIZE_16
      sub    x2, x2, #2
      sub    x3, x3, #27
      lsl    x1, x1, #1 //_dststride
      lsl    x3, x3, #1 //_srcstride
      DUP    V29.4S, w11 //offset vector
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      LD1    {V0.8B}, [x6]
      SADDW   V0.8H, V30.8H, V0.8B
0:
      subs x5, x5, #1
      LD1   {V25.8H-V28.8H}, [x4], x9
      bi1_epel_h_filter V14.8H, x13, V25
      bi1_epel_h_filter V15.8H, x13, V26
      bi1_epel_h_filter V16.8H, x13, V27
      bi1_epel_h_filter V17.8H, x3, V28
      ST1   {V14.8H-V17.8H}, [x0], x1
      b.ne 0b
      ret lr
endfunc

.macro bi1_epel_v_filter_4lines, out, in
      LD1    {V1.8H}, [\in], x3
      LD1    {V2.8H}, [\in], x3
      LD1    {V3.8H}, [\in], x3
0:              subs   x5, x5, #4
      LD1    {V4.8H}, [\in], x3
      LD1    {V28.8H}, [x4], x9
      bi1_epel_filter V17.8H, v28, V0, V1, V2, V3, V4
      LD1    {V1.8H}, [\in], x3
      LD1    {V28.8H}, [x4], x9
      bi1_epel_filter V18.8H, v28, V0, V2, V3, V4, V1
      LD1    {V2.8H}, [\in], x3
      LD1    {V28.8H}, [x4], x9
      bi1_epel_filter V19.8H, v28, V0, V3, V4, V1, V2
      LD1    {V3.8H}, [\in], x3
      LD1    {V28.8H}, [x4], x9
      bi1_epel_filter V20.8H, v28, V0, V4, V1, V2, V3
      ST1   {V17.8H}, [\out], x1
      ST1   {V18.8H}, [\out], x1
      ST1   {V19.8H}, [\out], x1
      ST1   {V20.8H}, [\out], x1
      b.ne 0b
.endm

.macro bi1_epel_v_filter_2lines, out, in
      LD1    {V1.8H}, [\in], x3
      LD1    {V2.8H}, [\in], x3
      LD1    {V3.8H}, [\in], x3
0:              subs   x5, x5, #2
      LD1    {V4.8H}, [\in], x3
      LD1    {V28.8H}, [x4], x9
      bi1_epel_filter V17.8H, v28, V0, V1, V2, V3, V4
      LD1    {V1.8H}, [\in], x3
      LD1    {V28.8H}, [x4], x9
      bi1_epel_filter V18.8H, v28, V0, V2, V3, V4, V1
      ST1   {V17.8H}, [\out], x1
      ST1   {V18.8H}, [\out], x1
      b.ne 0b
.endm

function put_vvc_bi1_epel_v_10_8_neon_4lines, export=1
      mov    x12, #2
      mov    w11, #64 //offset val
      mov    w10, #0x03FF //max val
      mov    x9, MAX_PB_SIZE_16
      lsl    x1, x1, #1 //_dststride
      lsl    x3, x3, #1 //_srcstride
      sub    x2, x2, x3
      DUP    V29.4S, w11 //offset vector
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      LD1    {V0.8B}, [x6]
      SADDW   V0.8H, V30.8H, V0.8B
      bi1_epel_v_filter_4lines x0, x2
      ret lr
endfunc

function put_vvc_bi1_epel_v_10_8_neon_2lines, export=1
      mov    x12, #2
      mov    w11, #64 //offset val
      mov    w10, #0x03FF //max val
      mov    x9, MAX_PB_SIZE_16
      lsl    x1, x1, #1 //_dststride
      lsl    x3, x3, #1 //_srcstride
      sub    x2, x2, x3
      DUP    V29.4S, w11 //offset vector
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      LD1    {V0.8B}, [x6]
      SADDW   V0.8H, V30.8H, V0.8B
      bi1_epel_v_filter_2lines x0, x2
      ret lr
endfunc

function put_vvc_bi1_epel_v_10_16_neon_4lines, export=1
      mov    w11, #64 //offset val
      mov    w10, #0x03FF //max val
      mov    x9, MAX_PB_SIZE_16
      lsl    x1, x1, #1 //_dststride
      lsl    x3, x3, #1 //_srcstride
      sub    x2, x2, x3
      mov    x13, x0  //save x0[0]
      mov    x12, x2  //save x2[-_srcstride]
      mov    x14, x4  //save x4[0]
      mov    x15, x5  //save height
      DUP    V29.4S, w11 //offset vector
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      LD1    {V0.8B}, [x6]
      SADDW   V0.8H, V30.8H, V0.8B
      bi1_epel_v_filter_4lines x0, x2
      add x2, x12, #16
      add x0, x13, #16
      add x4, x14, #16
      mov x5, x15
      bi1_epel_v_filter_4lines x0, x2
      ret lr
endfunc

function put_vvc_bi1_epel_v_10_16_neon_2lines, export=1
      mov    w11, #64 //offset val
      mov    w10, #0x03FF //max val
      mov    x9, MAX_PB_SIZE_16
      lsl    x1, x1, #1 //_dststride
      lsl    x3, x3, #1 //_srcstride
      sub    x2, x2, x3
      mov    x13, x0  //save x0[0]
      mov    x12, x2  //save x2[-_srcstride]
      mov    x14, x4  //save x4[0]
      mov    x15, x5  //save height
      DUP    V29.4S, w11 //offset vector
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      LD1    {V0.8B}, [x6]
      SADDW   V0.8H, V30.8H, V0.8B
      bi1_epel_v_filter_2lines x0, x2
      add x2, x12, #16
      add x0, x13, #16
      add x4, x14, #16
      mov x5, x15
      bi1_epel_v_filter_2lines x0, x2
      ret lr
endfunc

function put_vvc_bi1_epel_v_10_32_neon_4lines, export=1
      mov    w11, #64 //offset val
      mov    w10, #0x03FF //max val
      mov    x9, MAX_PB_SIZE_16
      lsl    x1, x1, #1 //_dststride
      lsl    x3, x3, #1 //_srcstride
      sub    x2, x2, x3
      mov    x13, x0  //save x0[0]
      mov    x12, x2  //save x2[-_srcstride]
      mov    x14, x4  //save x4[0]
      mov    x15, x5  //save height
      DUP    V29.4S, w11 //offset vector
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      LD1    {V0.8B}, [x6]
      SADDW   V0.8H, V30.8H, V0.8B
      bi1_epel_v_filter_4lines x0, x2
      add x2, x12, #16
      add x0, x13, #16
      add x4, x14, #16
      mov x5, x15
      bi1_epel_v_filter_4lines x0, x2
      add x2, x12, #32
      add x0, x13, #32
      add x4, x14, #32
      mov x5, x15
      bi1_epel_v_filter_4lines x0, x2
      add x2, x12, #48
      add x0, x13, #48
      add x4, x14, #48
      mov x5, x15
      bi1_epel_v_filter_4lines x0, x2
      ret lr
endfunc

function put_vvc_bi1_epel_v_10_32_neon_2lines, export=1
      mov    w11, #64 //offset val
      mov    w10, #0x03FF //max val
      mov    x9, MAX_PB_SIZE_16
      lsl    x1, x1, #1 //_dststride
      lsl    x3, x3, #1 //_srcstride
      sub    x2, x2, x3
      mov    x13, x0  //save x0[0]
      mov    x12, x2  //save x2[-_srcstride]
      mov    x14, x4  //save x4[0]
      mov    x15, x5  //save height
      DUP    V29.4S, w11 //offset vector
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      LD1    {V0.8B}, [x6]
      SADDW   V0.8H, V30.8H, V0.8B
      bi1_epel_v_filter_2lines x0, x2
      add x2, x12, #16
      add x0, x13, #16
      add x4, x14, #16
      mov x5, x15
      bi1_epel_v_filter_2lines x0, x2
      add x2, x12, #32
      add x0, x13, #32
      add x4, x14, #32
      mov x5, x15
      bi1_epel_v_filter_2lines x0, x2
      add x2, x12, #48
      add x0, x13, #48
      add x4, x14, #48
      mov x5, x15
      bi1_epel_v_filter_2lines x0, x2
      ret lr
endfunc

.macro bi1_epel_hv_filter_4lines, out, in, coefx=V0, coefy=V1, last_jump=x3, jump=x12
      bi0_epel_h_filter V2.8H, \last_jump, \coefx, \jump, \in, V6, V7, V8, V9
      bi0_epel_h_filter V3.8H, \last_jump, \coefx, \jump, \in, V6, V7, V8, V9
      bi0_epel_h_filter V4.8H, \last_jump, \coefx, \jump, \in, V6, V7, V8, V9
0:              subs x5, x5, #4
      bi0_epel_h_filter V5.8H, \last_jump, \coefx, \jump, \in, V6, V7, V8, V9
      LD1    {V28.8H}, [x4], x9
      bi1_epel_filter V20.8H, V28, \coefy, V2, V3, V4, V5, shift=#6
      bi0_epel_h_filter V2.8H, \last_jump, \coefx, \jump, \in, V6, V7, V8, V9
      LD1    {V28.8H}, [x4], x9
      bi1_epel_filter V21.8H, V28, \coefy, V3, V4, V5, V2, shift=#6
      bi0_epel_h_filter V3.8H, \last_jump, \coefx, \jump, \in, V6, V7, V8, V9
      LD1    {V28.8H}, [x4], x9
      bi1_epel_filter V22.8H, V28, \coefy, V4, V5, V2, V3, shift=#6
      bi0_epel_h_filter V4.8H, \last_jump, \coefx, \jump, \in, V6, V7, V8, V9
      LD1    {V28.8H}, [x4], x9
      bi1_epel_filter V23.8H, V28, \coefy, V5, V2, V3, V4, shift=#6
      ST1   {V20.8H}, [\out], x1
      ST1   {V21.8H}, [\out], x1
      ST1   {V22.8H}, [\out], x1
      ST1   {V23.8H}, [\out], x1
      b.ne 0b
.endm

.macro bi1_epel_hv_filter_2lines, out, in, coefx=V0, coefy=V1, last_jump=x3, jump=x12
      bi0_epel_h_filter V2.8H, \last_jump, \coefx, \jump, \in, V6, V7, V8, V9
      bi0_epel_h_filter V3.8H, \last_jump, \coefx, \jump, \in, V6, V7, V8, V9
      bi0_epel_h_filter V4.8H, \last_jump, \coefx, \jump, \in, V6, V7, V8, V9
0:              subs x5, x5, #2
      bi0_epel_h_filter V5.8H, \last_jump, \coefx, \jump, \in, V6, V7, V8, V9
      LD1    {V28.8H}, [x4], x9
      bi1_epel_filter V20.8H, V28, \coefy, V2, V3, V4, V5, shift=#6
      bi0_epel_h_filter V2.8H, \last_jump, \coefx, \jump, \in, V6, V7, V8, V9
      LD1    {V28.8H}, [x4], x9
      bi1_epel_filter V21.8H, V28, \coefy, V3, V4, V5, V2, shift=#6
      ST1   {V20.8H}, [\out], x1
      ST1   {V21.8H}, [\out], x1
      b.ne 0b
.endm

function put_vvc_bi1_epel_hv_10_8_neon_4lines, export=1
      mov    x9, MAX_PB_SIZE_16
      mov    w11, #1024 //offset val
      mov    w10, #0x03FF //max val
      DUP    V29.4S, w11 //offset vector
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      mov    x12, #2
      sub    x2, x2, #2//start h src[]
      sub    x13, x3, #3
      lsl    x13, x13, #1 // last jump h
      lsl    x1, x1, #1 //_dststride
      lsl    x14, x3, #1 //_srcstride
      sub    x2, x2, x14//start hv src[]
      LD1    {V0.8B}, [x6]
      LD1    {V1.8B}, [x7]
      SADDW   V0.8H, V30.8H, V0.8B
      SADDW   V1.8H, V30.8H, V1.8B
      bi1_epel_hv_filter_4lines x0, x2, V0, V1, x13, x12
      ret lr
endfunc

function put_vvc_bi1_epel_hv_10_8_neon_2lines, export=1
      mov    x9, MAX_PB_SIZE_16
      mov    w11, #1024 //offset val
      mov    w10, #0x03FF //max val
      DUP    V29.4S, w11 //offset vector
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      mov    x12, #2
      sub    x2, x2, #2//start h src[]
      sub    x13, x3, #3
      lsl    x13, x13, #1 // last jump h
      lsl    x1, x1, #1 //_dststride
      lsl    x14, x3, #1 //_srcstride
      sub    x2, x2, x14//start hv src[]
      LD1    {V0.8B}, [x6]
      LD1    {V1.8B}, [x7]
      SADDW   V0.8H, V30.8H, V0.8B
      SADDW   V1.8H, V30.8H, V1.8B
      bi1_epel_hv_filter_2lines x0, x2, V0, V1, x13, x12
      ret lr
endfunc

function put_vvc_bi1_epel_hv_10_16_neon_4lines, export=1
      mov    x9, MAX_PB_SIZE_16
      mov    w11, #1024 //offset val
      mov    w10, #0x03FF //max val
      DUP    V29.4S, w11 //offset vector
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      mov    x12, #2
      sub    x2, x2, #2//start h src[]
      sub    x13, x3, #3
      lsl    x13, x13, #1 // last jump h
      lsl    x1, x1, #1 //_dststride
      lsl    x14, x3, #1 //_srcstride
      sub    x2, x2, x14//start hv src[]
      mov    x16, x5//save height
      mov    x17, x2//save start src
      mov    x15, x0//save start dst
      mov    x14, x4//save start dst
      LD1    {V0.8B}, [x6]
      LD1    {V1.8B}, [x7]
      SADDW   V0.8H, V30.8H, V0.8B
      SADDW   V1.8H, V30.8H, V1.8B
      bi1_epel_hv_filter_4lines x0, x2, V0, V1, x13, x12
      add x0, x15, #16
      add x2, x17, #16
      add x4, x14, #16
      mov x5, x16
      bi1_epel_hv_filter_4lines x0, x2, V0, V1, x13, x12
      ret lr
endfunc

function put_vvc_bi1_epel_hv_10_16_neon_2lines, export=1
      mov    x9, MAX_PB_SIZE_16
      mov    w11, #1024 //offset val
      mov    w10, #0x03FF //max val
      DUP    V29.4S, w11 //offset vector
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      mov    x12, #2
      sub    x2, x2, #2//start h src[]
      sub    x13, x3, #3
      lsl    x13, x13, #1 // last jump h
      lsl    x1, x1, #1 //_dststride
      lsl    x14, x3, #1 //_srcstride
      sub    x2, x2, x14//start hv src[]
      mov    x16, x5//save height
      mov    x17, x2//save start src
      mov    x15, x0//save start dst
      mov    x14, x4//save start dst
      LD1    {V0.8B}, [x6]
      LD1    {V1.8B}, [x7]
      SADDW   V0.8H, V30.8H, V0.8B
      SADDW   V1.8H, V30.8H, V1.8B
      bi1_epel_hv_filter_2lines x0, x2, V0, V1, x13, x12
      add x0, x15, #16
      add x2, x17, #16
      add x4, x14, #16
      mov x5, x16
      bi1_epel_hv_filter_2lines x0, x2, V0, V1, x13, x12
      ret lr
endfunc

function put_vvc_bi1_epel_hv_10_32_neon_4lines, export=1
      mov    x9, MAX_PB_SIZE_16
      mov    w11, #1024 //offset val
      mov    w10, #0x03FF //max val
      DUP    V29.4S, w11 //offset vector
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      mov    x12, #2
      sub    x2, x2, #2//start h src[]
      sub    x13, x3, #3
      lsl    x13, x13, #1 // last jump h
      lsl    x1, x1, #1 //_dststride
      lsl    x14, x3, #1 //_srcstride
      sub    x2, x2, x14//start hv src[]
      mov    x16, x5//save height
      mov    x17, x2//save start src
      mov    x15, x0//save start dst
      mov    x14, x4//save start dst
      LD1    {V0.8B}, [x6]
      LD1    {V1.8B}, [x7]
      SADDW   V0.8H, V30.8H, V0.8B
      SADDW   V1.8H, V30.8H, V1.8B
      bi1_epel_hv_filter_4lines x0, x2, V0, V1, x13, x12
      add x0, x15, #16
      add x2, x17, #16
      add x4, x14, #16
      mov x5, x16
      bi1_epel_hv_filter_4lines x0, x2, V0, V1, x13, x12
      add x0, x15, #32
      add x2, x17, #32
      add x4, x14, #32
      mov x5, x16
      bi1_epel_hv_filter_4lines x0, x2, V0, V1, x13, x12
      add x0, x15, #48
      add x2, x17, #48
      add x4, x14, #48
      mov x5, x16
      bi1_epel_hv_filter_4lines x0, x2, V0, V1, x13, x12
      ret lr
endfunc

function put_vvc_bi1_epel_hv_10_32_neon_2lines, export=1
      mov    x9, MAX_PB_SIZE_16
      mov    w11, #1024 //offset val
      mov    w10, #0x03FF //max val
      DUP    V29.4S, w11 //offset vector
      DUP    V30.8H, wzr //min val vector
      DUP    V31.8H, w10 //max val vector
      mov    x12, #2
      sub    x2, x2, #2//start h src[]
      sub    x13, x3, #3
      lsl    x13, x13, #1 // last jump h
      lsl    x1, x1, #1 //_dststride
      lsl    x14, x3, #1 //_srcstride
      sub    x2, x2, x14//start hv src[]
      mov    x16, x5//save height
      mov    x17, x2//save start src
      mov    x15, x0//save start dst
      mov    x14, x4//save start dst
      LD1    {V0.8B}, [x6]
      LD1    {V1.8B}, [x7]
      SADDW   V0.8H, V30.8H, V0.8B
      SADDW   V1.8H, V30.8H, V1.8B
      bi1_epel_hv_filter_2lines x0, x2, V0, V1, x13, x12
      add x0, x15, #16
      add x2, x17, #16
      add x4, x14, #16
      mov x5, x16
      bi1_epel_hv_filter_2lines x0, x2, V0, V1, x13, x12
      add x0, x15, #32
      add x2, x17, #32
      add x4, x14, #32
      mov x5, x16
      bi1_epel_hv_filter_2lines x0, x2, V0, V1, x13, x12
      add x0, x15, #48
      add x2, x17, #48
      add x4, x14, #48
      mov x5, x16
      bi1_epel_hv_filter_2lines x0, x2, V0, V1, x13, x12
      ret lr
endfunc
