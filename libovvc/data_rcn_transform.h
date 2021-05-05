#ifndef RCN_TRANSFORM_DATA_H
#define RCN_TRANSFORM_DATA_H

#include <stdint.h>

extern const int16_t DCT_II_2[4];
extern const int16_t DCT_II_4[16];
extern const int16_t DCT_II_8[64];
extern const int16_t DCT_II_16[256];
extern const int16_t DCT_II_32[1024];
extern const int16_t DCT_II_64_OT[16*32];
extern const int16_t DCT_II_64_EOT[8*16];
extern const int16_t DCT_II_64_EEOT[4*8];
extern const int16_t DCT_II_64_EEEOT[2*4];
extern const int16_t DCT_II_64_EEEEO[2];
extern const int16_t DCT_II_64_EEEEE[2];

extern const int16_t DCT_VIII_4[16];
extern const int16_t DCT_VIII_8[64];
extern const int16_t DCT_VIII_16[256];
extern const int16_t DCT_VIII_32[1024];

extern const int16_t DST_VII_4[16];
extern const int16_t DST_VII_8[64];
extern const int16_t DST_VII_16[256];
extern const int16_t DST_VII_32[1024];

extern const int8_t* const lfnst[2][4][2];

#endif // RCN_TRANSFORM_DATA_H
