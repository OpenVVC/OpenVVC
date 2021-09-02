#include "arm/asm.S"

#define MAX_PB_SIZE #64
#define MAX_PB_SIZE_8 #128
#define MAX_PB_SIZE_16 #256
#define SHIFT_VAL #4 //14 - bitdepth, for 10 -> 4, for 8 -> 6

.macro clip, val, min_val, max_val
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
                DUP    V3.4H, wzr //min val vector
                DUP    V4.4H, w10 //max val vector
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

.macro uni_qpel_h_filter, out, last_jump
      LD1   {V1.8H}, [x2], x12
      LD1   {V2.8H}, [x2], x12
      LD1   {V3.8H}, [x2], x12
      LD1   {V4.8H}, [x2], x12
      LD1   {V5.8H}, [x2], x12
      LD1   {V6.8H}, [x2], x12
      LD1   {V7.8H}, [x2], x12
      LD1   {V8.8H}, [x2], \last_jump
      SMULL2 V9.4S, V1.8H, V0.H[0]
      SMLAL2 V9.4S, V2.8H, V0.H[1]
      SMULL2 V10.4S, V3.8H, V0.H[2]
      SMLAL2 V10.4S, V4.8H, V0.H[3]
      SMULL2 V11.4S, V5.8H, V0.H[4]
      SMLAL2 V11.4S, V6.8H, V0.H[5]
      SMULL2 V12.4S, V7.8H, V0.H[6]
      SMLAL2 V12.4S, V8.8H, V0.H[7]
      SMULL V1.4S, V1.4H, V0.H[0]
      SMLAL V1.4S, V2.4H, V0.H[1]
      SMULL V2.4S, V3.4H, V0.H[2]
      SMLAL V2.4S, V4.4H, V0.H[3]
      SMULL V3.4S, V5.4H, V0.H[4]
      SMLAL V3.4S, V6.4H, V0.H[5]
      SMULL V4.4S, V7.4H, V0.H[6]
      SMLAL V4.4S, V8.4H, V0.H[7]
      ADD   V9.4S, V9.4S, V10.4S
      ADD   V10.4S, V11.4S, V12.4S
      ADD   V1.4S, V1.4S, V2.4S
      ADD   V2.4S, V3.4S, V4.4S
      ADD   V1.4S, V1.4S, V2.4S
      ADD   V2.4S, V9.4S, V10.4S
      ADD   V1.4S, V1.4S, V29.4S
      ADD   V2.4S, V2.4S, V29.4S
      SSHR  V1.4S, V1.4S, #6
      SSHR  V2.4S, V2.4S, #6
      UZP1  \out, V1.8H, V2.8H
      CLIP  \out, V30.8H, V31.8H
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
                uni_qpel_h_filter V14.8H, x3
                uni_qpel_h_filter V15.8H, x3
                ST1   {V14.8H}, [x0], x1
                ST1   {V15.8H}, [x0], x1
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
                uni_qpel_h_filter V14.8H, x12
                uni_qpel_h_filter V15.8H, x3
                ST1   {V14.8H-V15.8H}, [x0], x1
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
                uni_qpel_h_filter V14.8H, x12
                uni_qpel_h_filter V15.8H, x12
                uni_qpel_h_filter V16.8H, x12
                uni_qpel_h_filter V17.8H, x3
                ST1   {V14.8H-V17.8H}, [x0], x1
                b.ne 0b
                ret lr
endfunc

.macro bi0_qpel_h_filter, out, last_jump
      LD1   {V1.8H}, [x1], x12
      LD1   {V2.8H}, [x1], x12
      LD1   {V3.8H}, [x1], x12
      LD1   {V4.8H}, [x1], x12
      LD1   {V5.8H}, [x1], x12
      LD1   {V6.8H}, [x1], x12
      LD1   {V7.8H}, [x1], x12
      LD1   {V8.8H}, [x1], \last_jump
      SMULL2 V9.4S, V1.8H, V0.H[0]
      SMLAL2 V9.4S, V2.8H, V0.H[1]
      SMULL2 V10.4S, V3.8H, V0.H[2]
      SMLAL2 V10.4S, V4.8H, V0.H[3]
      SMULL2 V11.4S, V5.8H, V0.H[4]
      SMLAL2 V11.4S, V6.8H, V0.H[5]
      SMULL2 V12.4S, V7.8H, V0.H[6]
      SMLAL2 V12.4S, V8.8H, V0.H[7]
      SMULL V1.4S, V1.4H, V0.H[0]
      SMLAL V1.4S, V2.4H, V0.H[1]
      SMULL V2.4S, V3.4H, V0.H[2]
      SMLAL V2.4S, V4.4H, V0.H[3]
      SMULL V3.4S, V5.4H, V0.H[4]
      SMLAL V3.4S, V6.4H, V0.H[5]
      SMULL V4.4S, V7.4H, V0.H[6]
      SMLAL V4.4S, V8.4H, V0.H[7]
      ADD   V9.4S, V9.4S, V10.4S
      ADD   V10.4S, V11.4S, V12.4S
      ADD   V1.4S, V1.4S, V2.4S
      ADD   V2.4S, V3.4S, V4.4S
      ADD   V1.4S, V1.4S, V2.4S
      ADD   V2.4S, V9.4S, V10.4S
      SSHR  V1.4S, V1.4S, #2
      SSHR  V2.4S, V2.4S, #2
      UZP1  \out, V1.8H, V2.8H
.endm

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
                bi0_qpel_h_filter V14.8H, x2
                bi0_qpel_h_filter V15.8H, x2
                ST1   {V14.8H}, [x0], x11
                ST1   {V15.8H}, [x0], x11
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
                bi0_qpel_h_filter V14.8H, x12
                bi0_qpel_h_filter V15.8H, x2
                ST1   {V14.8H-V15.8H}, [x0], x11
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
                bi0_qpel_h_filter V14.8H, x12
                bi0_qpel_h_filter V15.8H, x12
                bi0_qpel_h_filter V16.8H, x12
                bi0_qpel_h_filter V17.8H, x2
                ST1   {V14.8H-V17.8H}, [x0], x11
                b.ne 0b
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
      SMULL2 V9.4S, V1.8H, V0.H[0]
      SMLAL2 V9.4S, V2.8H, V0.H[1]
      SMULL2 V10.4S, V3.8H, V0.H[2]
      SMLAL2 V10.4S, V4.8H, V0.H[3]
      SMULL2 V11.4S, V5.8H, V0.H[4]
      SMLAL2 V11.4S, V6.8H, V0.H[5]
      SMULL2 V12.4S, V7.8H, V0.H[6]
      SMLAL2 V12.4S, V8.8H, V0.H[7]
      SMULL V1.4S, V1.4H, V0.H[0]
      SMLAL V1.4S, V2.4H, V0.H[1]
      SMULL V2.4S, V3.4H, V0.H[2]
      SMLAL V2.4S, V4.4H, V0.H[3]
      SMULL V3.4S, V5.4H, V0.H[4]
      SMLAL V3.4S, V6.4H, V0.H[5]
      SMULL V4.4S, V7.4H, V0.H[6]
      SMLAL V4.4S, V8.4H, V0.H[7]
      ADD   V9.4S, V9.4S, V10.4S
      ADD   V10.4S, V11.4S, V12.4S
      ADD   V1.4S, V1.4S, V2.4S
      ADD   V2.4S, V3.4S, V4.4S
      ADD   V1.4S, V1.4S, V2.4S
      ADD   V2.4S, V9.4S, V10.4S
      ADD   V1.4S, V1.4S, V29.4S
      ADD   V2.4S, V2.4S, V29.4S
      SSHR  V1.4S, V1.4S, #2
      SSHR  V2.4S, V2.4S, #2
      SADDW V1.4S, V1.4S,  \src1\suffix1
      SADDW2 V2.4S, V2.4S, \src1\suffix2
      SSHR  V1.4S, V1.4S, #5
      SSHR  V2.4S, V2.4S, #5
      UZP1  \out, V1.8H, V2.8H
      CLIP  \out, V30.8H, V31.8H
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
                bi1_qpel_h_filter V14.8H, x3, V27
                bi1_qpel_h_filter V15.8H, x3, V28
                ST1   {V14.8H}, [x0], x1
                ST1   {V15.8H}, [x0], x1
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
                bi1_qpel_h_filter V14.8H, x12, V27
                bi1_qpel_h_filter V15.8H, x3, V28
                ST1   {V14.8H-V15.8H}, [x0], x1
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
                bi1_qpel_h_filter V14.8H, x12, V25
                bi1_qpel_h_filter V15.8H, x12, V26
                bi1_qpel_h_filter V16.8H, x12, V27
                bi1_qpel_h_filter V17.8H, x3, V28
                ST1   {V14.8H-V17.8H}, [x0], x1
                b.ne 0b
                ret lr
endfunc

.macro uni_epel_h_filter, out, last_jump
      LD1   {V1.8H}, [x2], x12
      LD1   {V2.8H}, [x2], x12
      LD1   {V3.8H}, [x2], x12
      LD1   {V4.8H}, [x2], \last_jump
      SMULL2 V9.4S, V1.8H, V0.H[0]
      SMLAL2 V9.4S, V2.8H, V0.H[1]
      SMULL2 V10.4S, V3.8H, V0.H[2]
      SMLAL2 V10.4S, V4.8H, V0.H[3]
      SMULL V1.4S, V1.4H, V0.H[0]
      SMLAL V1.4S, V2.4H, V0.H[1]
      SMULL V2.4S, V3.4H, V0.H[2]
      SMLAL V2.4S, V4.4H, V0.H[3]
      ADD   V1.4S, V1.4S, V2.4S
      ADD   V2.4S, V9.4S, V10.4S
      ADD   V1.4S, V1.4S, V29.4S
      ADD   V2.4S, V2.4S, V29.4S
      SSHR  V1.4S, V1.4S, #6
      SSHR  V2.4S, V2.4S, #6
      UZP1  \out, V1.8H, V2.8H
      CLIP  \out, V30.8H, V31.8H
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
                mov    w11, #32 //offset val
                mov    w10, #0x03FF //max val
                sub    x2, x2, #2
                sub    x3, x3, #7
                lsl    x1, x1, #1 //_dststride
                lsl    x3, x3, #1 //_srcstride
                DUP    V29.4S, w11 //offset vector
                DUP    V30.8H, wzr //min val vector
                DUP    V31.8H, w10 //max val vector
                LD1    {V0.8B}, [x5]
                SADDW   V0.8H, V30.8H, V0.8B
0:
                subs x4, x4, #1
                uni_epel_h_filter V14.8H, x12
                uni_epel_h_filter V15.8H, x3
                ST1   {V14.8H-V15.8H}, [x0], x1
                b.ne 0b
                ret lr
endfunc

.macro bi0_epel_h_filter, out, last_jump
      LD1   {V1.8H}, [x1], x12
      LD1   {V2.8H}, [x1], x12
      LD1   {V3.8H}, [x1], x12
      LD1   {V4.8H}, [x1], \last_jump
      SMULL2 V9.4S, V1.8H, V0.H[0]
      SMLAL2 V9.4S, V2.8H, V0.H[1]
      SMULL2 V10.4S, V3.8H, V0.H[2]
      SMLAL2 V10.4S, V4.8H, V0.H[3]
      SMULL V1.4S, V1.4H, V0.H[0]
      SMLAL V1.4S, V2.4H, V0.H[1]
      SMULL V2.4S, V3.4H, V0.H[2]
      SMLAL V2.4S, V4.4H, V0.H[3]
      ADD   V1.4S, V1.4S, V2.4S
      ADD   V2.4S, V9.4S, V10.4S
      SSHR  V1.4S, V1.4S, #2
      SSHR  V2.4S, V2.4S, #2
      UZP1  \out, V1.8H, V2.8H
.endm

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
                mov    x11, MAX_PB_SIZE_16 //offset val
                mov    w10, #0x03FF //max val
                sub    x1, x1, #2
                sub    x2, x2, #7
                lsl    x2, x2, #1 //_srcstride
                DUP    V30.8H, wzr //min val vector
                DUP    V31.8H, w10 //max val vector
                LD1    {V0.8B}, [x4]
                SADDW   V0.8H, V30.8H, V0.8B
0:
                subs x3, x3, #1
                bi0_epel_h_filter V14.8H, x12
                bi0_epel_h_filter V15.8H, x2
                ST1   {V14.8H-V15.8H}, [x0], x11
                b.ne 0b
                ret lr
endfunc

.macro bi1_epel_h_filter, out, last_jump, src1, suffix1=.4H, suffix2=.8H
      LD1   {V1.8H}, [x2], x12
      LD1   {V2.8H}, [x2], x12
      LD1   {V3.8H}, [x2], x12
      LD1   {V4.8H}, [x2], \last_jump
      SMULL2 V9.4S, V1.8H, V0.H[0]
      SMLAL2 V9.4S, V2.8H, V0.H[1]
      SMULL2 V10.4S, V3.8H, V0.H[2]
      SMLAL2 V10.4S, V4.8H, V0.H[3]
      SMULL V1.4S, V1.4H, V0.H[0]
      SMLAL V1.4S, V2.4H, V0.H[1]
      SMULL V2.4S, V3.4H, V0.H[2]
      SMLAL V2.4S, V4.4H, V0.H[3]
      ADD   V1.4S, V1.4S, V2.4S
      ADD   V2.4S, V9.4S, V10.4S
      ADD   V1.4S, V1.4S, V29.4S
      ADD   V2.4S, V2.4S, V29.4S
      SSHR  V1.4S, V1.4S, #2
      SSHR  V2.4S, V2.4S, #2
      SADDW V1.4S, V1.4S,  \src1\suffix1
      SADDW2 V2.4S, V2.4S, \src1\suffix2
      SSHR  V1.4S, V1.4S, #5
      SSHR  V2.4S, V2.4S, #5
      UZP1  \out, V1.8H, V2.8H
      CLIP  \out, V30.8H, V31.8H
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
                mov    w11, #64 //offset val
                mov    w10, #0x03FF //max val
                mov    x9, MAX_PB_SIZE_16
                sub    x2, x2, #2
                sub    x3, x3, #7
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
                bi1_epel_h_filter V14.8H, x12, V27
                bi1_epel_h_filter V15.8H, x3, V28
                ST1   {V14.8H-V15.8H}, [x0], x1
                b.ne 0b
                ret lr
endfunc
