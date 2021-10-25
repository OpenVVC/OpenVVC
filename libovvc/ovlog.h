#ifndef OVLOG_H
#define OVLOG_H

#include <stdarg.h>

typedef enum
{
    OVLOG_ERROR,
    OVLOG_WARNING,
    OVLOG_INFO,
    OVLOG_VERBOSE,
    OVLOG_TRACE,
    OVLOG_DEBUG
} OVLogLevel;

void ovlog_set_log_level(OVLogLevel log_level);

void ov_log(void* ctx, int log_level, const char* log_content, ...);

void ovlog_set_callback(void (*log_function)(void* ctx, int log_level, const char* log_content, va_list vl));

#endif

