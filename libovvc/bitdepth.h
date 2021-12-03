#ifndef BITDEPTH_H
#define BITDEPTH_H

#if BITDEPTH > 8
typedef uint8_t OVSample;
#else
typedef uint8_t OVSample;
#endif

#define ov_bdclip(a) ov_clip_uintp2(a, BITDEPTH)
#define AVG_VAL ((1 << (BITDEPTH)) >> 1)

#define FUNC3(a, b, c)  a ## _ ## b ##  c
#define FUNC2(a, b, c)  FUNC3(a, b, c)
#define BD_DECL(a)  FUNC2(a, BITDEPTH, )

#endif
