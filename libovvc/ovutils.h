/**
 *
 *   OpenVVC is open-source real time software decoder compliant with the 
 *   ITU-T H.266- MPEG-I - Part 3 VVC standard. OpenVVC is developed from 
 *   scratch in C as a library that provides consumers with real time and
 *   energy-aware decoding capabilities under different OS including MAC OS,
 *   Windows, Linux and Android targeting low energy real-time decoding of
 *   4K VVC videos on Intel x86 and ARM platforms.
 * 
 *   Copyright (C) 2020-2022  IETR-INSA Rennes :
 *   
 *   Pierre-Loup CABARAT
 *   Wassim HAMIDOUCHE
 *   Guillaume GAUTIER
 *   Thomas AMESTOY
 *
 *   This library is free software; you can redistribute it and/or
 *   modify it under the terms of the GNU Lesser General Public
 *   License as published by the Free Software Foundation; either
 *   version 2.1 of the License, or (at your option) any later version.
 *
 *   This library is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *   Lesser General Public License for more details.
 *
 *   You should have received a copy of the GNU Lesser General Public
 *   License along with this library; if not, write to the Free Software
 *   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301
 *   USA
 * 
 **/

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
    int32_t lim = (1 << b);
    uint32_t overflow_msk = ~((lim << 1) - 1);
    uint32_t clip = ((uint32_t)val + lim) & overflow_msk;
    if (clip) {
        uint32_t msk = (1 << b) - 1;
        return (val >> 31) ^ msk;
    } else {
        return val;
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
