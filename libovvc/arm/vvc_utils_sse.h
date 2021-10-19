#ifndef VVC_UTILS_SSE_H
#define VVC_UTILS_SSE_H

#include <stdint.h>
#include <stddef.h>
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

typedef int16_t TCoeff;
typedef int16_t TMatrixCoeff;

void inverse_sse2_B4(const TCoeff *src, TCoeff *dst, int src_stride, int shift, int line, const TMatrixCoeff* iT);
void inverse_sse2_B8(const TCoeff *src, TCoeff *dst, int src_stride, int shift, int line, const TMatrixCoeff* iT);
void inverse_sse2_B16(const TCoeff *src, TCoeff *dst, int src_stride, int shift, int line, const TMatrixCoeff* iT);
void inverse_sse2_B32(const TCoeff *src, TCoeff *dst, int src_stride, int shift, int line, const TMatrixCoeff* iT);

void afficherVecteur4SSE128(__m128i cible);
void afficherVecteur8SSE128(__m128i cible);
void afficherTableau2D(int v, int h, TCoeff * m);

#endif//VVC_UTILS_SSE_H
