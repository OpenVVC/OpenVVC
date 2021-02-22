#ifndef OVVCUTILS_H
#define OVVCUTILS_H

#include <stdint.h>

#define RED   "\x1B[31m"
#define GRN   "\x1B[32m"
#define YEL   "\x1B[33m"
#define BLU   "\x1B[34m"
#define MAG   "\x1B[35m"
#define CYN   "\x1B[36m"
#define WHT   "\x1B[37m"
#define RST "\x1B[0m"

#define OVMAX(a, b) (((a) > (b)) ? (a) : (b))
#define OVMIN(a, b) (((a) < (b)) ? (a) : (b))
#define OVABS(a) (((a) < (0)) ? -(a) : (a))

typedef enum
{
    OVLOG_ERROR,
    OVLOG_WARNING,
    OVLOG_INFO,
    OVLOG_VERBOSE,
    OVLOG_DEBUG,
    OVLOG_TRACE
} OVLOG_TYPE;

extern OVLOG_TYPE ov_log_level;



void
set_ov_log_level(OVLOG_TYPE log_level);

void
ov_log(void* ctx, int log_level, const char* log_content, ...);

/* FIXME
 * Add specific clip for unsigned */
int32_t
ov_clip(int32_t val, int32_t a, int32_t b);

uint32_t
ov_clip_uintp2(uint32_t val, uint32_t a);

int
floor_log2(unsigned x);

#endif
