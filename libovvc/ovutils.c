#include <stdint.h>

#include "ovutils.h"
#include "ovconfig.h"
#if _WIN32
#include <windows.h>
#else
#include <unistd.h>
#endif

/* FIXME check if used on negative numbers */

int
get_number_of_cores() {
#if _WIN32
    SYSTEM_INFO sysinfo;
    GetSystemInfo(&sysinfo);
    return sysinfo.dwNumberOfProcessors;
#else
    return sysconf(_SC_NPROCESSORS_ONLN);
#endif
}
