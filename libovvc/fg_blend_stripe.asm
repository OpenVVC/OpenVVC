%define ARCH_X86_64 1
%define private_prefix asm

%include "x86inc.asm"

    section .rodata align=16
zeros_vec:  dd  0, 0, 0, 0
ones_vec:   dd  1, 1, 1, 1
neg_ones_vec:   dd  -1, -1, -1, -1

    section .data align=16
mask:       dd  1, 1, 1, 1

    section .bss  align=16
not_mask:   resd  4         ; reserve 4 dword (128 bits)


    section .text
INIT_XMM sse4

; dstSampleOffsetY  => r0 (rdi)
; srcSampleOffsetY  => r1 (rsi)
; grainStripe       => r2 (rdx)
; widthComp         => r3 (rcx)
; blockHeight       => r4 (R8)
; bitDepth          => r5 (R9)
; X                 => r6 (rax)
; ???               => r7 (R10)
; ???               => r8 (R11)
; X                 => r9 (rbx)
; X                 => r10 (rbp)
; X                 => r11 (r14)
; ???               => r12 (r15)
; ???               => r13 (r12)
; ???               => r14 (r13)
; ???               => r15 (rsp)


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
	

; function_name, #args (%1), #regs (%2), #xmm_regs (%3), [stack_size,] (%4) arg_names... (%5-*)
cglobal fg_blend_stripe, 6, 6, 7, dstSampleOffsetY, srcSampleOffsetY, grainStripe, widthComp, blockHeight, bitDepth, row, col
; declares a function  that:
; - automatically loads 6 arguments (dstSampleOffsetY, srcSampleOffsetY, grainStripe, widthComp, blockHeight, bitDepth)
; - uses 0 additional register ()
; - 7 vector registers (m0-m7)
; - perform no allocation of stack space.
    
    .begin_global:
        ; Clear registers
        pxor m0, m0
        pxor m1, m1
        pxor m2, m2
        pxor m3, m3
        pxor m4, m4
        pxor m5, m5
        pxor m6, m6
        pxor m7, m7
 
        ;  Prepare SIMD SSE4 ov_clip_uintp2: mask & not_mask
        ; for loop
        push r3q
        mov r3q, bitDepthq

        ; Restore after each call values of mask
        mov DWORD [mask], 1
        mov DWORD [mask + 4], 1
        mov DWORD [mask + 8], 1
        mov DWORD [mask + 12], 1
        
        shl DWORD [mask], r3b
        shl DWORD [mask + 4], r3b
        shl DWORD [mask + 8], r3b
        shl DWORD [mask + 12], r3b
        
        mov  r3w, [mask]
        mov  [not_mask], r3q
        xor  [not_mask], r3q
        sub [not_mask], r3q
        mov  r3q, 1
        sub [mask], r3q

        mov  r3w, [mask + 4]
        mov  [not_mask + 4], r3q
        xor  [not_mask + 4], r3q
        sub [not_mask + 4], r3q
        mov  r3q, 1
        sub [mask + 4], r3q
        
        mov  r3w, [mask + 8]
        mov  [not_mask + 8], r3q
        xor  [not_mask + 8], r3q
        sub [not_mask + 8], r3q
        mov  r3q, 1
        sub [mask + 8], r3q
        
        mov  r3w, [mask + 12]
        mov  [not_mask + 12], r3q
        xor  [not_mask + 12], r3q
        sub [not_mask + 12], r3q
        mov  r3q, 1
        sub [mask + 12], r3q

        pop r3q

        push    r11q
        push    r9q
        
        ; left shift
        mov r11q, bitDepthq 
        sub r11q, 8

        push r11q ; make sure it is on 16 byte boundary

        push    r6q
        push    r10q
        mov     r10q, rsp
        sub     r10q, 96 ; Make space for row (8 bytes), col (8 bytes), double push r9q for alignment (16 bytes), 16 chunks of int32_t vals (64 bytes)
       
    .core_global:        

        mov  QWORD [rsp-32], 0 ; row
        
        .L1:
            ; 1 iteration = 4 rows
            cmp [rsp-32], blockHeightq ; row < blockHeight
            jge  .end_global

            mov  QWORD [rsp-24], 0 ; col
            mov r9q, 0
            ; row offset: We need to skip (#row * widthComp * 4 bytes)
            add r9q, [rsp-32]
            imul r9q, widthCompq
            imul r9q, 2 ; int32_t = 4 bytes (we later multiply by 2 to make it 4)

            .L2_begin:
                cmp [rsp-24], widthCompq ; col < widthComp
                jge  .L1_end
                
                ; Push twice to keep 16 bytes alignment
                push r9q
                push r9q

                ; 2 registers = total load of 8 values on the same row
                ; row 1
                ; get the first 4 chunks of int32_t
                mova  m0, [grainStripeq + r9q * 2] ; multiply by the remaining 2 to make it 4
                mova  m1, [grainStripeq + r9q * 2 + 16]
                ; Left shift
                pslld m0, [rsp + 32]
                pslld m1, [rsp + 32]
                
                ; add srcSampleOffsetY
                ; packed move with sign extension
                pmovsxwd m2, [srcSampleOffsetYq + r9q]
                pmovsxwd m3, [srcSampleOffsetYq + r9q + 8]
                paddd m0, m2
                paddd m1, m3

                ; SIMD SSE4 ov_clip_uintp2
                ;  Set to 0 all negative values.
                ; grainSample = _mm_max_epi32(grainSample, _mm_setzero_si128());
                pmaxsd m0, [zeros_vec]
                pmaxsd m1, [zeros_vec]

                ; Save grainSample
                movdqu  [rbp], m0
                movdqu  [rbp + 16], m1
           
                ; int32_t overflow = !!(val & (~mask));
                ; __m128i overflow = _mm_and_si128(grainSample, not_mask);
                pand m0, [not_mask]
                pand m1, [not_mask]
                ; overflow = _mm_min_epi32(overflow, _mm_set1_epi32(1));
                pminsd m0, [ones_vec] 
                pminsd m1, [ones_vec]
                ; overflow = _mm_sub_epi32(_mm_set1_epi32(0), overflow);
                pmulld m0, [neg_ones_vec]
                pmulld m1, [neg_ones_vec]

                ; ((-overflow) & mask) | (val & mask);
                ; __m128i lhs = _mm_and_si128(overflow, mask);
                pand m0, [mask]
                pand m1, [mask]

                ; Save lhs
                movdqu  [rbp + 32], m0
                movdqu [rbp + 48], m1

                ; Retrieve grainSample
                movdqu  m0, [rbp]
                movdqu  m1, [rbp + 16]
                
                ;__m128i rhs = _mm_and_si128(grainSample, mask);
                pand m0, [mask]
                pand m1, [mask]

                ; __m128i clipped_val = _mm_or_si128(lhs, rhs);
                ; rhs OR lhs
                por m0, [rbp + 32]
                por m1, [rbp + 48]

                ; Update dstSampleOffsetY
                pextrd  eax, m0, 0
                mov     WORD [dstSampleOffsetYq + r9q + 0], ax
                pextrd  eax, m0, 1
                mov     WORD [dstSampleOffsetYq + r9q + 2], ax
                pextrd  eax, m0, 2
                mov     WORD [dstSampleOffsetYq + r9q + 4], ax
                pextrd  eax, m0, 3
                mov     WORD [dstSampleOffsetYq + r9q + 6], ax
                
                pextrd  eax, m1, 0
                mov     WORD [dstSampleOffsetYq + r9q + 8], ax
                pextrd  eax, m1, 1
                mov     WORD [dstSampleOffsetYq + r9q + 10], ax
                pextrd  eax, m1, 2
                mov     WORD [dstSampleOffsetYq + r9q + 12], ax
                pextrd  eax, m1, 3
                mov     WORD [dstSampleOffsetYq + r9q + 14], ax
                
                ; row 2
                lea  r9q, [r9q + widthCompq * 2]
                mova  m2, [grainStripeq + r9q * 2]
                mova  m3, [grainStripeq + r9q * 2 + 16]
                
                pslld m2, [rsp + 32]
                pslld m3, [rsp + 32]

                pmovsxwd m0, [srcSampleOffsetYq + r9q]
                pmovsxwd m1, [srcSampleOffsetYq + r9q + 8]
                paddd m2, m0
                paddd m3, m1

                pmaxsd m2, [zeros_vec]
                pmaxsd m3, [zeros_vec]

                movdqu  [rbp], m2
                movdqu  [rbp + 16], m3

                pand m2, [not_mask]
                pand m3, [not_mask]

                pminsd m2, [ones_vec] 
                pminsd m3, [ones_vec]

                pmulld m2, [neg_ones_vec]
                pmulld m3, [neg_ones_vec]

                pand m2, [mask]
                pand m3, [mask]

                movdqu  [rbp + 32], m2
                movdqu [rbp + 48], m3

                movdqu  m2, [rbp]
                movdqu  m3, [rbp + 16]
                
                pand m2, [mask]
                pand m3, [mask]

                por m2, [rbp + 32]
                por m3, [rbp + 48]

                pextrd  eax, m2, 0
                mov     WORD [dstSampleOffsetYq + r9q + 0], ax
                pextrd  eax, m2, 1
                mov     WORD [dstSampleOffsetYq + r9q + 2], ax
                pextrd  eax, m2, 2
                mov     WORD [dstSampleOffsetYq + r9q + 4], ax
                pextrd  eax, m2, 3
                mov     WORD [dstSampleOffsetYq + r9q + 6], ax
                
                pextrd  eax, m3, 0
                mov     WORD [dstSampleOffsetYq + r9q + 8], ax
                pextrd  eax, m3, 1
                mov     WORD [dstSampleOffsetYq + r9q + 10], ax
                pextrd  eax, m3, 2
                mov     WORD [dstSampleOffsetYq + r9q + 12], ax
                pextrd  eax, m3, 3
                mov     WORD [dstSampleOffsetYq + r9q + 14], ax

                ; row 3
                lea  r9q, [r9q + widthCompq * 2]
                mova  m4, [grainStripeq + r9q * 2]
                mova  m5, [grainStripeq + r9q * 2 + 16]
                
                pslld m4, [rsp + 32]
                pslld m5, [rsp + 32]

                pmovsxwd m0, [srcSampleOffsetYq + r9q]
                pmovsxwd m1, [srcSampleOffsetYq + r9q + 8]
                paddd m4, m0
                paddd m5, m1

                pmaxsd m4, [zeros_vec]
                pmaxsd m5, [zeros_vec]

                movdqu  [rbp], m4
                movdqu  [rbp + 16], m5

                pand m4, [not_mask]
                pand m5, [not_mask]

                pminsd m4, [ones_vec] 
                pminsd m5, [ones_vec]

                pmulld m4, [neg_ones_vec]
                pmulld m5, [neg_ones_vec]

                pand m4, [mask]
                pand m5, [mask]

                movdqu  [rbp + 32], m4
                movdqu [rbp + 48], m5

                movdqu  m4, [rbp]
                movdqu  m5, [rbp + 16]
                
                pand m4, [mask]
                pand m5, [mask]

                por m4, [rbp + 32]
                por m5, [rbp + 48]

                pextrd  eax, m4, 0
                mov     WORD [dstSampleOffsetYq + r9q + 0], ax
                pextrd  eax, m4, 1
                mov     WORD [dstSampleOffsetYq + r9q + 2], ax
                pextrd  eax, m4, 2
                mov     WORD [dstSampleOffsetYq + r9q + 4], ax
                pextrd  eax, m4, 3
                mov     WORD [dstSampleOffsetYq + r9q + 6], ax
                
                pextrd  eax, m5, 0
                mov     WORD [dstSampleOffsetYq + r9q + 8], ax
                pextrd  eax, m5, 1
                mov     WORD [dstSampleOffsetYq + r9q + 10], ax
                pextrd  eax, m5, 2
                mov     WORD [dstSampleOffsetYq + r9q + 12], ax
                pextrd  eax, m5, 3
                mov     WORD [dstSampleOffsetYq + r9q + 14], ax

                ; row 4
                lea  r9q, [r9q + widthCompq * 2]
                mova  m6, [grainStripeq + r9q * 2]
                mova  m7, [grainStripeq + r9q * 2 + 16]

                pslld m6, [rsp + 32]
                pslld m7, [rsp + 32]

                pmovsxwd m0, [srcSampleOffsetYq + r9q]
                pmovsxwd m1, [srcSampleOffsetYq + r9q + 8]
                paddd m6, m0
                paddd m7, m1

                pmaxsd m6, [zeros_vec]
                pmaxsd m7, [zeros_vec]

                movdqu  [rbp], m6
                movdqu  [rbp + 16], m7

                pand m6, [not_mask]
                pand m7, [not_mask]

                pminsd m6, [ones_vec] 
                pminsd m7, [ones_vec]

                pmulld m6, [neg_ones_vec]
                pmulld m7, [neg_ones_vec]

                pand m6, [mask]
                pand m7, [mask]

                movdqu  [rbp + 32], m6
                movdqu [rbp + 48], m7

                movdqu  m6, [rbp]
                movdqu  m7, [rbp + 16]
                
                pand m6, [mask]
                pand m7, [mask]

                por m6, [rbp + 32]
                por m7, [rbp + 48]

                pextrd  eax, m6, 0
                mov     WORD [dstSampleOffsetYq + r9q + 0], ax
                pextrd  eax, m6, 1
                mov     WORD [dstSampleOffsetYq + r9q + 2], ax
                pextrd  eax, m6, 2
                mov     WORD [dstSampleOffsetYq + r9q + 4], ax
                pextrd  eax, m6, 3
                mov     WORD [dstSampleOffsetYq + r9q + 6], ax
                
                pextrd  eax, m7, 0
                mov     WORD [dstSampleOffsetYq + r9q + 8], ax
                pextrd  eax, m7, 1
                mov     WORD [dstSampleOffsetYq + r9q + 10], ax
                pextrd  eax, m7, 2
                mov     WORD [dstSampleOffsetYq + r9q + 12], ax
                pextrd  eax, m7, 3
                mov     WORD [dstSampleOffsetYq + r9q + 14], ax

                ; Prepare next iteration
                pop r9q
                pop r9q
                
                add QWORD [rsp-24], 8
                add  r9q, 16

                jmp .L2_begin

        .L1_end:
            add QWORD [rsp-32], 4
            jmp .L1

    .end_global:
        add r10q, 96
        mov rsp, r10q
        pop r10q
        pop r6q
        pop r11q
        pop r9q
        pop r11q

        RET