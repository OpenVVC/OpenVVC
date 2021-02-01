#ifndef OVVCUTILS_H
#define OVVCUTILS_H

#include <stdint.h>

#define OVMAX(a, b) (((a) > (b)) ? (a) : (b))
#define OVMIN(a, b) (((a) < (b)) ? (a) : (b))

void
ov_log(void* ctx, int log_level, const char* log_content, ...);

uint32_t
ov_clip(uint32_t val, uint32_t a, uint32_t b);

uint32_t
ov_clip_uintp2(uint32_t val, uint32_t a);

int
floor_log2(unsigned x);

#endif
