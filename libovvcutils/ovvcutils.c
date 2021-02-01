#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>

#include "ovvcutils.h"

static int ov_log_level = 1;
static const char* vvctype = "VVCDec";

void
ov_log(void* ctx, int log_level, const char* log_content, ...)
{
        va_list args;

        va_start(args, log_content);

        if (log_level > ov_log_level) {
                const char* type = "NULL";
                if (ctx != NULL) {
                        type = vvctype;
                }
                printf("[%s @ Ox%.16lx] : ", type, (long unsigned)ctx);
                vprintf(log_content, args);
        }

        va_end(args);
}

uint32_t
ov_clip(uint32_t val, uint32_t a, uint32_t b)
{
        return OVMIN(OVMAX(val, a), b);
}

uint32_t
ov_clip_uintp2(uint32_t val, uint32_t a)
{
        return OVMIN(val, a);
}

int
floor_log2(unsigned x)
{
        int bits = -1;
        while (x > 0) {
                bits++;
                x >>= 1;
        }
        return bits;
}
