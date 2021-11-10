#ifndef OVVCUTILS_H
#define OVVCUTILS_H

#include <stdarg.h>
#include <stdint.h>

#include "ovlog.h"


#define OVMAX(a, b) (((a) > (b)) ? (a) : (b))
#define OVMIN(a, b) (((a) < (b)) ? (a) : (b))
#define OVABS(a) (((a) < (0)) ? -(a) : (a))

#define ov_clz(x) __builtin_clz(x)
#define ov_ctz(x) __builtin_ctz(x)

#define ov_clz64(x) __builtin_clzl(x)
#define ov_ctz64(x) __builtin_ctzl(x)

#define ov_ceil_log2(x) 32 - __builtin_clz(x)

/* FIXME
 * Add specific clip for unsigned */
static inline int32_t
ov_clip(int32_t val, int32_t a, int32_t b)
{
    return OVMIN(OVMAX(val, a), b);
}

uint32_t
ov_clip_uintp2(int32_t val, uint32_t a);

int32_t
ov_clip_intp2(int32_t val, uint32_t a);

static inline int
floor_log2(unsigned x)
{
#if 0
    int bits = -1;
    while (x > 0) {
        bits++;
        x >>= 1;
    }
    return bits;
#else
    return 32 - ov_clz(x) - 1;
#endif
}

int get_number_of_cores();

#endif
