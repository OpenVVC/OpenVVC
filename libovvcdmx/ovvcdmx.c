#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <stdint.h>

#include "libovvcutils/ovmem.h"
#include "libovvcutils/ovvcutils.h"

#include "ovvcdmx.h"


/* Use Fixed size of demux read cache buffer to 64K */

#define OVVCDMX_IO_BUFF_SIZE (1 << 16)

#define OVVCDMX_IO_BUFF_MASK ((1 << 16) - 1)

static const char *const demux_name = "Open VVC Annex B demuxer";

typedef struct OVReadBuff{
    uint8_t *data;
}OVReadBuff;

struct OVVCDmx{
    const char *name;

    FILE *fstream;
    OVReadBuff cache_buffer;

    struct{
        int val;
    }options;
};

static int ovread_buff_alloc(OVVCDmx *const dmx, size_t buff_size);

static void ovread_buff_free(OVVCDmx *const dmx);

int
ovdmx_init(OVVCDmx **vvcdmx)
{
    *vvcdmx = ov_mallocz(sizeof(**vvcdmx));

    (*vvcdmx)->name = demux_name;

    if (*vvcdmx == NULL) return -1;

    return 0;
}

int
ovdmx_close(OVVCDmx *vvcdmx)
{
    int not_dmx = 0;
    if (vvcdmx != NULL) {

        not_dmx = vvcdmx->name != demux_name;

        if (not_dmx) goto fail;

        ov_free(vvcdmx);

        return 0;
    }

fail:
    ov_log(vvcdmx, 3, "Trying to close a something not a demuxer.\n");
    return -1;
}

int
ovdmx_attach_stream(OVVCDmx *const dmx, FILE *fstream)
{
    /* FiXME is this check necessary this function should not
       be called if stream is not allocated.
       Maybe we should open file ourselves / and use a wrapper around
       I/Os */
    int ret;
    if (fstream == NULL) {
        return -1;
    }

    dmx->fstream = fstream;

    ret = ovread_buff_alloc(dmx, OVVCDMX_IO_BUFF_SIZE);

    return ret;
}

void
ovdmx_detach_stream(OVVCDmx *const dmx)
{
    dmx->fstream = NULL;
    ovread_buff_free(dmx);
}

static int
ovread_buff_alloc(OVVCDmx *const dmx, size_t buff_size)
{
    uint8_t *byte_stream;

    /* Note we keep a 8 bytes left and 8 bytes padding at the right of
       our buffer in order to check,to store the last bytes of the
       previous chunk of the bytestream 8 bytes is because we might
       want to use 64bit types in order to quickly probe for successive
       zero bytes */
    byte_stream = ov_mallocz(buff_size + 8 + 8);

    dmx->cache_buffer.data = byte_stream;

    if (byte_stream == NULL) {
        return -1;
    }

    return 1;
}

void
ovread_buff_free(OVVCDmx *const dmx)
{
    ov_freep(&dmx->cache_buffer.data);
}

