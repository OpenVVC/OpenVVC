#include <stdint.h>

#include "ovutils.h"
#include "ovconfig.h"
#if _WIN32
#include <windows.h>
#else
#include <unistd.h>
#endif

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

int32_t
ov_clip_intp2(int32_t val, uint32_t a)
{   
    int b = a-1;
    if (val > 0) {
        int32_t mask  = (1 << b) - 1;
        int32_t overflow = !!(val & (~mask));
        return ((-overflow) & mask) | (val & mask); 
    } else {
        val = -val;
        int32_t mask  = (1 << b) - 1;
        int32_t overflow = !!(val & (~mask));
        return -(((-overflow) & mask) | (val & mask));
    }
}

int get_number_of_cores() {
#if _WIN32
    SYSTEM_INFO sysinfo;
    GetSystemInfo(&sysinfo);
    return sysinfo.dwNumberOfProcessors;
#else
    return sysconf(_SC_NPROCESSORS_ONLN);
#endif
}
