#include <stdlib.h>

#include "libovvcutils/ovmem.h"
#include "libovvcutils/ovvcutils.h"

#include "ovvcdmx.h"

static const char *const demux_name = "Open VVC Annex B demuxer";

struct OVVCDmx {
    const char *name;
    int val;

    struct {
        int val;
    }options;
};

int
ovdmx_init(OVVCDmx **vvcdmx)
{
    *vvcdmx = malloc(sizeof(**vvcdmx));
    
    (*vvcdmx)->name = demux_name;

    if (*vvcdmx == NULL) return -1;

    return 0;
}

int
ovdmx_close(OVVCDmx *vvcdmx)
{
    int not_dmx;
    if (vvcdmx != NULL) {

        not_dmx = vvcdmx->name != demux_name;

        if (not_dmx) goto fail;

        free(vvcdmx);

        return 0;
    }

fail:
    ov_log(vvcdmx, 3, "Trying to close a something not a demuxer.\n");
    return -1;
}

