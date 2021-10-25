#include <stdio.h>
#include "ovlog.h"

OVLogLevel ov_log_level = OVLOG_INFO;

static const char* vvctype = "VVCDec";

#define RED   "\x1B[31m"
#define GRN   "\x1B[32m"
#define YEL   "\x1B[33m"
#define BLU   "\x1B[34m"
#define MAG   "\x1B[35m"
#define CYN   "\x1B[36m"
#define WHT   "\x1B[37m"
#define RST "\x1B[0m"

static const char *OVLOG_COLORIFY[6] = { RED, YEL, BLU, CYN, GRN, MAG};

void
ovlog_set_log_level(OVLogLevel log_level)
{
    ov_log_level = log_level;
}

static void
ov_log_default(void* ctx, int log_level, const char* log_content, va_list vl)
{
    if (log_level <= ov_log_level) {
        const char* type = "NULL";
        if (ctx != NULL) {
            type = vvctype;
        }
        fprintf(stderr, "%s", OVLOG_COLORIFY[log_level]);
        fprintf(stderr, "[%s @ %p] : ", type, ctx);
        vfprintf(stderr, log_content, vl);
        fprintf(stderr, "%s", RST);
    }
}

static void (*ov_log_callback)(void* ctx, int log_level, const char* log_content, va_list vl) = ov_log_default;

void
ovlog_set_callback(void (*log_function)(void* ctx, int log_level, const char* log_content, va_list vl))
{
    ov_log_callback = log_function;
}

void
ov_log(void* ctx, int log_level, const char* log_content, ...)
{
        va_list args;

        va_start(args, log_content);

        ov_log_callback(ctx, log_level, log_content, args);

        va_end(args);
}

