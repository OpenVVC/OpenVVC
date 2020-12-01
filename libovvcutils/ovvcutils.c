#include <stdio.h>
#include <stdarg.h>

#include "ovvcutils.h"

static int ov_log_level = 1;
static const char* vvctype = "VVCDec";

void
ov_log(void* ctx, int log_level, const char* log_content, ... )
{
    va_list args;

    va_start(args, log_content);

    if (log_level > ov_log_level){
        const char *type = "NULL";
        if (ctx != NULL){
            type = vvctype;
        }
        printf("[%s @ Ox%.16lx] : ", type, (long unsigned)ctx);
        vprintf(log_content, args);
    }

    va_end(args);
}

