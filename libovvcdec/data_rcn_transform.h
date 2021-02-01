#ifndef RCN_TRANSFORM_DATA_H
#define RCN_TRANSFORM_DATA_H

#include <stdint.h>

extern const int16_t DCT_II_2[4];
extern const int16_t DCT_II_4[16];
extern const int16_t DCT_II_8[64];
extern const int16_t DCT_II_16[256];
extern const int16_t DCT_II_32[1024];
extern const int16_t DCT_II_64[4096];

extern const int16_t DCT_VIII_4[16];
extern const int16_t DCT_VIII_8[64];
extern const int16_t DCT_VIII_16[256];
extern const int16_t DCT_VIII_32[1024];

extern const int16_t DST_VII_4[16];
extern const int16_t DST_VII_8[64];
extern const int16_t DST_VII_16[256];
extern const int16_t DST_VII_32[1024];

extern const int8_t* const lfnst_8x8[4][2];
extern const int8_t lfnst_0_0_8x8[16 * 48];
extern const int8_t lfnst_0_1_8x8[16 * 48];
extern const int8_t lfnst_1_0_8x8[16 * 48];
extern const int8_t lfnst_1_1_8x8[16 * 48];
extern const int8_t lfnst_2_0_8x8[16 * 48];
extern const int8_t lfnst_2_1_8x8[16 * 48];
extern const int8_t lfnst_3_0_8x8[16 * 48];
extern const int8_t lfnst_3_1_8x8[16 * 48];

extern const int8_t* const lfnst_4x4[4][2];
extern const int8_t lfnst_0_0_4x4[16 * 16];
extern const int8_t lfnst_0_1_4x4[16 * 16];
extern const int8_t lfnst_1_0_4x4[16 * 16];
extern const int8_t lfnst_1_1_4x4[16 * 16];
extern const int8_t lfnst_2_0_4x4[16 * 16];
extern const int8_t lfnst_2_1_4x4[16 * 16];
extern const int8_t lfnst_3_0_4x4[16 * 16];
extern const int8_t lfnst_3_1_4x4[16 * 16];
#endif // RCN_TRANSFORM_DATA_H
