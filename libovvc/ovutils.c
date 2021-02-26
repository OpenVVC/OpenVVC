#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>

#include "ovutils.h"

OVLOG_TYPE ov_log_level = OVLOG_INFO;
#define ov_clz(x) __builtin_clz(x)

static const char* vvctype = "VVCDec";

static const char *OVLOG_COLORIFY[6] = { RED, YEL, BLU, CYN, GRN, MAG};

void
set_ov_log_level(OVLOG_TYPE log_level){
        ov_log_level = log_level;
}

void
ov_log(void* ctx, int log_level, const char* log_content, ...)
{
        va_list args;

        va_start(args, log_content);

        if (log_level <= ov_log_level) {
                const char* type = "NULL";
                if (ctx != NULL) {
                        type = vvctype;
                }
                fprintf(stderr, "%s", OVLOG_COLORIFY[log_level]);
                fprintf(stderr, "[%s @ Ox%.16lx] : ", type, (long unsigned)ctx);
                vfprintf(stderr, log_content, args);
                fprintf(stderr, "%s", RST);
        }

        va_end(args);
}

int32_t
ov_clip(int32_t val, int32_t a, int32_t b)
{
        return OVMIN(OVMAX(val, a), b);
}

/* FIXME check if used on negative numbers */
uint32_t
ov_clip_uintp2(int32_t val, uint32_t a)
{
        return OVMIN(OVMAX(0, val), 1 << a);

}

int
floor_log2(unsigned x)
{
        #if 0
        int bits = -1;
        while (x > 0) {
                bits++;
                x >>= 1;
        }
        return bits;
        #endif
        return 31 - ov_clz(x);
}
