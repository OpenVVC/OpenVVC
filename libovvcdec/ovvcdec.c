#include <stdlib.h>

#include "ovvcdec.h"
#include "libovvcutils/ovvcutils.h"

static const char *const decname = "Open VVC Decoder";

struct OVVCDec{
    const char *name;
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

    (*vvcdec)->name = decname;

    return 0;

fail:
    /* TODO proper error management (ENOMEM)*/
    return -1;
}

int
ovdec_close(OVVCDec *vvcdec)
{
    int not_dec;
    if (vvcdec != NULL) {

        not_dec = vvcdec->name != decname;

        if (not_dec) goto fail;

        free(vvcdec);

        return 0;
    }

fail:
    ov_log(vvcdec, 3, "Trying to close a something not a decoder.\n");
    return -1;
}
