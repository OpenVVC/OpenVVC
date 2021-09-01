#include "arm/asm.S"

#define MAX_PB_SIZE #64
#define MAX_PB_SIZE_8 #128
#define MAX_PB_SIZE_16 #256
#define OH_ONE_VECTOR #16
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
