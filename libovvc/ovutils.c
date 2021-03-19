#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>

#include "ovutils.h"

OVLOG_TYPE ov_log_level = OVLOG_INFO;

static const char* vvctype = "VVCDec";

static const char *OVLOG_COLORIFY[6] = { RED, YEL, BLU, CYN, GRN, MAG};

void
set_ov_log_level(OVLOG_TYPE log_level)
{
        ov_log_level = log_level;
}

void
ov_log(void* ctx, int log_level, const char* log_content, ...)
{
        va_list args;

        va_start(args, log_content);

        #if 0
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
        #endif

        va_end(args);
}

/* FIXME check if used on negative numbers */
uint32_t
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
