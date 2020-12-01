#include <stdlib.h>

#include "ovvcdmx.h"

struct OVVCDmx {
    int val;
    struct {
        int val;
    }my_struct;
};

int
ovdmx_init(OVVCDmx **vvcdmx)
{
    *vvcdmx = malloc(sizeof(**vvcdmx));

    return 0;
}

void
ovdmx_close(OVVCDmx **vvcdmx)
{
    if (*vvcdmx) free(*vvcdmx);

    *vvcdmx = NULL;
}

