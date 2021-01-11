#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <stdint.h>

#include "libovvcutils/ovvcutils.h"
#include "libovvcutils/ovmem.h"

#include "ovvcdmx.h"
#include "ovio.h"
#include "ovannexb.h"
#include "ovunits.h"


/* Use Fixed size of demux read cache buffer to 64K */

#define OVVCDMX_IO_BUFF_SIZE (1 << 16)

#define OVVCDMX_IO_BUFF_MASK ((1 << 16) - 1)

static const char *const demux_name = "Open VVC Annex B demuxer";


struct OVVCDmx
{
    const char *name;

    FILE *fstream;
    /*OVReadBuff cache_buffer;*/
    int nb_stc;
    int nb_epb;

    /* Points to a read only IO context */
    OVIOStream *io_str;

    /* Demuxer options to be passed at init */
    struct{
        int val;
    }options;
};

static int process_chunk(OVVCDmx *const dmx, const uint8_t *byte, uint64_t byte_pos);

static int process_last_chunk(OVVCDmx *const dmx, const uint8_t *byte, uint64_t byte_pos,
                              int nb_bytes_last);

int
ovdmx_init(OVVCDmx **vvcdmx)
{
    *vvcdmx = ov_mallocz(sizeof(**vvcdmx));

    if (*vvcdmx == NULL) return -1;

    (*vvcdmx)->name = demux_name;
    (*vvcdmx)->io_str = NULL;

    return 0;
}

int
ovdmx_close(OVVCDmx *vvcdmx)
{
    int not_dmx = 0;
    if (vvcdmx != NULL) {

        not_dmx = vvcdmx->name != demux_name;

        if (not_dmx) goto fail;

        ovdmx_detach_stream(vvcdmx);

        ov_free(vvcdmx);

        return 0;
    }

fail:
    ov_log(vvcdmx, 3, "Trying to close a something not a demuxer.\n");
    return -1;
}

/*FIXME do not use FILE use a wrapper to something more general instead
  so we can support other IO types */
int
ovdmx_attach_stream(OVVCDmx *const dmx, FILE *fstream)
{
    /* FiXME is this check necessary this function should not
       be called if stream is not allocated.
       Maybe we should open file ourselves / and use a wrapper around
       I/Os */
    int ret = 0;
    if (fstream == NULL) {
        return -1;
    }

    dmx->fstream = fstream;

    /* TODO distinguish init and open / attach */
    dmx->io_str = ovio_stream_open(fstream);
    if (dmx->io_str == NULL) {
        /* no mem*/
        return -1;
    }

    return ret;
}

/*FIXME share bytestream mem with OVIOStream so we avoid copying from one to
  another*/
void
ovdmx_detach_stream(OVVCDmx *const dmx)
{
    dmx->fstream = NULL;

    /* FIXME decide if it should free OVIOStream cache buff */
    if (dmx->io_str != NULL) {
        ovio_stream_close(dmx->io_str);
    }

    dmx->io_str = NULL;
}


/*
   FIXME there might be a need to find a fast way to remove obvious
   zero filling data between RBPS and next NALU
   FIXME we should move thi part to other file for annexb demux
   if we plan to support more demuxers format*/
int
ovdmx_read_stream(OVVCDmx *const dmx)
{
    uint64_t nb_chunks = 0;
    uint64_t byte_pos = 0;
    long int nb_bytes_last = 0;
    int nb_nalus = 0;
    OVIOStream *const io_str = dmx->io_str;
    const uint8_t *io_cache;

    if (io_str == NULL) {
       return -1;
    }

    /* Note if we need to read more than one chunk or
       if we uee a reentrant function we need to plan on not
       reseting this value since we might need the last few bytes
       of the cache used in previous call in cases */

    while (!ovio_stream_eof(io_str)){
       const uint8_t *byte;
       int read_in_buf = 0;
       int ret;

       /* request new read from io layer */
       read_in_buf += ovio_stream_read(&io_cache, OVVCDMX_IO_BUFF_SIZE, io_str);
       byte = io_cache - 8;

       if (!read_in_buf) {
           goto last_chunk;
       }

       nb_chunks += read_in_buf;

       ret = process_chunk(dmx, byte, byte_pos);
       if (ret < 0) {
           goto readfail;
       }

       printf("num_stc %d\n", dmx->nb_stc);
       printf("num_prev %d\n", dmx->nb_epb);
       nb_nalus += dmx->nb_stc;
    }

    return 1;

last_chunk:
    if (ovio_stream_error(io_str)) {
        return -1;
    } else {
       const uint64_t mask = OVVCDMX_IO_BUFF_MASK;
       int ret;

        /* we do not check return si error was already reported ?*/

        nb_bytes_last = ovio_stream_tell(io_str) & mask;

        if (ovio_stream_eof(io_str) && (nb_bytes_last > 0)) {

            const uint8_t *byte = io_cache - 8;

            ret = process_last_chunk(dmx, byte, byte_pos, nb_bytes_last);
            if (ret < 0) {
                goto readfail;
            }

            printf("num_stc last %d\n", dmx->nb_stc);
            printf("num_prev last %d\n", dmx->nb_epb);
            printf("EOF reached\n");
            nb_nalus += dmx->nb_stc;
        }
    }

    printf("Num bytes read %ld\n", (nb_chunks << 16) + nb_bytes_last);
    printf("byte_pos %ld\n", (nb_chunks << 16) + nb_bytes_last);
    printf("Num NALU read %d\n", nb_nalus);

    return 1;

readfail:
    /* free current NALU memory + return error
    */
    printf("FAILED NAL read!\n");
    return -1;
}

/**
 * returns: -1 Invalid data
 *          byte_pos in chunk if stc
 *          0 if nothing found and needs a new read
 */
static int
process_chunk(OVVCDmx *const dmx, const uint8_t *byte, uint64_t byte_pos)
{
    const uint64_t mask = OVVCDMX_IO_BUFF_MASK;
    do {
        /*TODO bintricks for start code detection */

        /*FIXME we will actually loop over this more than once even if a start
          code has been detected. This is a bit inefficient */

       /* WARNING We need to be careful on endianness here if we plan
          to use bigger read sizes */

        if (byte[(byte_pos) & mask] == 0) {
            int stc_or_epb;
            stc_or_epb = check_stc_or_epb(&byte[(byte_pos) & mask]);
            if (stc_or_epb < 0) {
                printf("Invalid\n");
                return -1;
            }

            if (stc_or_epb) {
                /* TODO handle what is to be done with it here */
                dmx_process_elem(dmx, &byte[(byte_pos) & mask], byte_pos, stc_or_epb);
                dmx->nb_stc += stc_or_epb == 1;
                dmx->nb_epb += stc_or_epb == 2;
                #if 0
                return (byte_pos & mask);
                #endif
            }
        }
    } while ((++byte_pos) & mask);

    return (byte_pos & mask);
}

static int
process_last_chunk(OVVCDmx *const dmx, const uint8_t *byte, uint64_t byte_pos,
                   int nb_bytes_last)
{
    const uint64_t mask = OVVCDMX_IO_BUFF_MASK;
    do {
        /*TODO bintricks for start code detection */
        if (byte[(byte_pos) & mask] == 0) {
            int stc_or_epb;
            stc_or_epb = check_stc_or_epb(&byte[(byte_pos) & mask]);

            if (stc_or_epb < 0) {
                return -1;
            }

            if (stc_or_epb) {
                /* TODO handle what is to be done with it here */
                dmx_process_elem(dmx, &byte[(byte_pos) & mask], byte_pos, stc_or_epb);
                dmx->nb_stc += stc_or_epb == 1;
                dmx->nb_epb += stc_or_epb == 2;
            }
        }
    } while (((++byte_pos) & mask) <= nb_bytes_last);
    /* FIXME handle EOF EOB stuff */
    return (byte_pos & mask);
}
