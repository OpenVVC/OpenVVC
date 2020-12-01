#include <stdlib.h>

#include "ovvcdec.h"

struct OVVCDec{
    int val;
    struct {
        int opt1;
    }options;
};


int
ovdec_init(OVVCDec **vvcdec)
{
    *vvcdec = malloc(sizeof(OVVCDec));
    if (*vvcdec == NULL) goto fail;

    return 0;

fail:
    /* TODO proper error management (ENOMEM)*/
    return -1;
}

void
ovdec_close(OVVCDec **vvcdec)
{
    if (*vvcdec != NULL) free(*vvcdec);

    *vvcdec = NULL;
}

