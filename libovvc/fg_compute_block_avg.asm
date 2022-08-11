%define ARCH_X86_64 1
%define private_prefix asm

%include "x86inc.asm"

    section .text
INIT_XMM sse4

; line 576 of x86inc.asm 
;DECLARE_REG 0,  rdi
;DECLARE_REG 1,  rsi
;DECLARE_REG 2,  rdx
;DECLARE_REG 3,  rcx
;DECLARE_REG 4,  R8
;DECLARE_REG 5,  R9
;DECLARE_REG 6,  rax, 8
;DECLARE_REG 7,  R10, 16
;DECLARE_REG 8,  R11, 24
;DECLARE_REG 9,  rbx, 32
;DECLARE_REG 10, rbp, 40
;DECLARE_REG 11, R14, 48
;DECLARE_REG 12, R15, 56
;DECLARE_REG 13, R12, 64
;DECLARE_REG 14, R13, 72


; dstSampeBlk8 => r0 (rdi)
; widthComp    => r1 (rsi)
; pNumSamples  => r2 (rdx)
; ySize        => r3 (rcx)
; xSize        => r4 (R8)
; bitDepth     => r5 (R9)
; blockAvg     => r6 (rax)
; numSamples   => r7 (R10)
; ???          => r8 (R11)
; ???          => r9 (rbx)
	

; function_name, #args (%1), #regs (%2), #xmm_regs (%3), [stack_size,] (%4) arg_names... (%5-*)
cglobal fg_compute_block_avg, 6, 8, 7, dstSampleBlk8, widthComp, pNumSamples, ySize, xSize, bitDepth \
                                     , blockAvg, numSamples
; declares a function  that:
; - automatically loads 6 arguments (dstSampleBlk8, widthComp, pNumSamples, ySize, xSize, bitDepth) into registers
; - uses 2 additional register (blockAvg, numSamples)
; - 7 vector registers (m0-m7)
; - perform no allocation of stack space.
    
    .begin_global:
        ; Compute total number of samples in the block 
        xor  blockAvgq, blockAvgq
        xor numSamplesq, numSamplesq
        
        ; From https://github.com/OpenVVC/OpenVVC/blob/master/libovvc/pp_film_grain.c#L910,
        ; xSize & ySize can either be 0 or 8.
        ; Thus, all possible combinations of xSize & ySize leads to 0, 64
        movzx   r6q, ySizeb
        mul     xSizeb
        mov     numSamplesw, r6w

        ; Clear registers
        pxor m0, m0
        pxor m1, m1
        pxor m2, m2
        pxor m3, m3
        pxor m4, m4
        pxor m5, m5
        pxor m6, m6
        pxor m7, m7

    .core_global:        
        ; Fetch 8 values at the time and stored them in vector registers.
        mova  m0, [dstSampleBlk8q + widthCompq * 0]
        mova  m1, [dstSampleBlk8q + widthCompq * 2]
        mova  m2, [dstSampleBlk8q + widthCompq * 4]

        lea  r6q, [widthCompq * 4]
        lea  r6q, [r6q + widthCompq * 2]
        mova  m3, [dstSampleBlk8q + r6q] ;[dstSampleBlk8q + widthCompq * 6]
        
        lea  r6q, [r6q + widthCompq * 2]
        mova  m4, [dstSampleBlk8q + r6q] ;[dstSampleBlk8q + widthCompq * 8] 
        
        lea  r6q, [r6q + widthCompq * 2]
        mova  m5, [dstSampleBlk8q + r6q] ;[dstSampleBlk8q + widthCompq * 10]

        lea  r6q, [r6q + widthCompq * 2]    
        mova  m6, [dstSampleBlk8q + r6q] ;[dstSampleBlk8q + widthCompq * 12]
        
        lea  r6q, [r6q + widthCompq * 2]
        mova  m7, [dstSampleBlk8q + r6q] ;[dstSampleBlk8q + widthCompq * 14]

        ; Add all the values in the vector registers
        paddd  m0, m1
        paddd  m2, m3
        paddd  m4, m5
        paddd  m6, m7

        paddd m0, m2
        paddd m4, m6
        paddd m0, m4
        
        ; Perfom horizontal sum 3 times. After that, every entries of m0 will have the same value.
        phaddw  m0, m0
        phaddw  m0, m0
        phaddw  m0, m0

        ; Store value for return value function (blockAvgw = rax)
        movd  r8d, m0
        mov   blockAvgw, r8w
        
        test   numSamplesq, numSamplesq
        jz     .end_global

        ; perfor; division
        shr   blockAvgq, 6 ; blockAvgq /= 64
        ; handle high bit depths
        push   r3q
        push  blockAvgq
            mov   r6b, bitDepthb
            sub   r6b, 8
            mov   r3b, r6b    
        pop   blockAvgq
        shr   blockAvgw, r3b
        pop   r3q

        call clip_uintp2

    .end_global:
        ; Update pNumSamples
        mov  WORD [pNumSamplesq], numSamplesw
        RET

clip_uintp2:
    .begin_func:
        ; val = blockAvg = rax
        push    r9
        push    r3
        push    r2
        
        push    r10
        mov     r10, rsp
        
    .core_func:

        test    r6d, r6d
        jle     .else

        .if:
            mov     r2d, 1
            mov     r3d, 8
            sal     r2d, r3b ; mask  = (1 << a);
            ; We need to substract 1 to mask later
            ; Otherwise we won't get the expected result when we negate it (ex: ~256 = -256)
            
            mov     r9d, r2d
            neg     r9d
            test    r9d, r6d
            setne   r9b
            sub     r2d, 1 ; mask -= 1;
            movzx   r3d, r9b ; overflow = (val & (~mask));
        
            neg     r3d
            and     r3d, r2d
            and     r2d, r6d
            or      r6d, r3d ; ((-overflow) & mask) | (val & mask)
            jmp .end_func
        .else:
            xor  r6d, r6d

    .end_func:
        mov rsp, r10
        pop r10

        pop    r2
        pop    r3
        pop    r9

        ret