#ifndef OVVCUTILS_H
#define OVVCUTILS_H

#include <stdint.h>

#define OVMAX(a, b) (((a) > (b)) ? (a) : (b))
#define OVMIN(a, b) (((a) < (b)) ? (a) : (b))
#define OVABS(a) (((a) < (0)) ? (-a) : (a))

enum OVLOG_TYPE {OVLOG_ERROR, OVLOG_WARNING, OVLOG_INFO, OVLOG_VERBOSE, OVLOG_DEBUG, OVLOG_TRACE };

void
ov_log(void* ctx, int log_level, const char* log_content, ...);

/*FIXME use separate clipping functions for
 * signed and unsigned clipping
 */
int32_t
ov_clip(int32_t val, int32_t a, int32_t b);

uint32_t
ov_clip_uintp2(uint32_t val, uint32_t a);

int
floor_log2(unsigned x);

#endif
