#include "ovvcutils.h"
#include "stdio.h"

static int ov_log_level = 1;

/*TODO proper structures initialisation etc. */
void
ov_log(void* ctx, int log_level, const char* log_content)
{
    if (log_level > ov_log_level){
        printf("[NULL @ Ox%.16lx] : %s", (long unsigned)ctx, log_content);
    }
}

