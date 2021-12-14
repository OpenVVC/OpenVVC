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

#define ov_ceil_log2(x) 32 - __builtin_clz((x - !!x) + !(x - !!x))

/* FIXME
 * Add specific clip for unsigned */
static inline int32_t
ov_clip(int32_t val, int32_t a, int32_t b)
{
    return OVMIN(OVMAX(val, a), b);
}

static inline uint32_t
ov_clip_uintp2(int32_t val, uint32_t a)
{
    if (val > 0) {
        int32_t mask  = (1 << a) - 1;
        int32_t overflow = !!(val & (~mask));
        return ((-overflow) & mask) | (val & mask);
    } else {
        return 0;
    }
    #if 0
    return OVMIN(OVMAX(0, val), (1 << a) - 1);
    #endif
}

static inline int32_t
ov_clip_intp2(int32_t val, uint32_t a)
{
    int b = a - 1;
    if (val > 0) {
        int32_t mask  = (1 << b) - 1;
        int32_t overflow = !!(val & (~mask));
        return ((-overflow) & mask) | (val & mask);
    } else {
        val = -val;
        int32_t mask  = (1 << b) - 1;
        int32_t overflow = !!(val & (~mask));
        return -(((-overflow) & mask) | (val & mask));
    }
}

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
    return 31 - ov_clz(x + !x);
#endif
}

int get_number_of_cores();

#endif
