#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <stdint.h>

#include "libovvcutils/ovvcutils.h"
#include "libovvcutils/ovmem.h"

#include "ovvcdmx.h"

/* Use Fixed size of demux read cache buffer to 64K */

#define OVVCDMX_IO_BUFF_SIZE (1 << 16)

#define OVVCDMX_IO_BUFF_MASK ((1 << 16) - 1)

static const char *const demux_name = "Open VVC Annex B demuxer";

typedef struct OVReadBuff{
    uint8_t *bytestream;
    void *opaque;
}OVReadBuff;

struct OVVCDmx{
    const char *name;

    FILE *fstream;
    OVReadBuff cache_buffer;

    struct{
        int val;
    }options;
};


static int ovread_buff_init(struct OVReadBuff *const cache_buff,
                             size_t buff_size);

static void ovread_buff_close(struct OVReadBuff *const cache_buff);

static inline int check_stc_or_epb(const uint8_t *byte);

static int dmx_process_elem(OVVCDmx *const dmx,
                            const uint8_t *const bytestream,
                            uint64_t byte_pos,
                            int stc_or_epb);

/* TODO create a wrapper for  IO functions names  should
   mimic stdio functions since they are basically offering
   the same interface */
static size_t io_read_stream(OVVCDmx *const dmx);

static int io_eof_stream(OVVCDmx *const dmx);

static int io_error_stream(OVVCDmx *const dmx);

static long int io_tell_stream(OVVCDmx *const dmx);

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

/*FIXME do not use FILE use a wrapper to something more general instead
  so we can support other IO types */
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

    ret = ovread_buff_init(&dmx->cache_buffer, OVVCDMX_IO_BUFF_SIZE);

    return ret;
}

void
ovdmx_detach_stream(OVVCDmx *const dmx)
{
    dmx->fstream = NULL;
    ovread_buff_close(&dmx->cache_buffer);
}

static int
ovread_buff_init(struct OVReadBuff *const cache_buff, size_t buff_size)
{
    uint8_t *byte_stream;

    /* Note we keep a 8 bytes left and 8 bytes padding at the right of
     * our buffer in order to check,to store the last bytes of the
     * previous chunk of the bytestream 8 bytes is because we might
     * want to use 64bit types in order to quickly probe for successive
     * zero bytes
     */
    byte_stream = ov_mallocz(8 + buff_size + 8);

    if (byte_stream == NULL) {
        return -1;
    }

    cache_buff->opaque = (void *)byte_stream;

    /* Bytestream will be cached from this position in allocated memory
     * the padding byte will be used to keep a copy of the last 8 bytes
     * in order to check for start codes overlapping between two chunks
     */
    cache_buff->bytestream = byte_stream + 8;

    /* last 16 bytes are set to 0xFF so we do not detect any zero byte
     * past the actual available data when checking for a start or emulation
     * prevention code.
     * Using. 0xFF shoul prevent patterns such as 0x000003 or 0x000001
     * we used 16 bytes so first copy from read function will copy 0XFF
     * bytes at the averlapping area reader will then ignore them since
     * it cannot be taken as start code.
     */
    memset(byte_stream + 8 + buff_size, 0XFF, sizeof(*byte_stream) * 8);

    return 1;
}

void
ovread_buff_close(struct OVReadBuff *const cache_buff)
{
    cache_buff->bytestream = NULL;

    ov_freep(&cache_buff->opaque);
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
    uint8_t *const cache_buffer = dmx->cache_buffer.bytestream;

    if (cache_buffer == NULL) {
       return -1;
    }

    /* Note if we need to read more than one chunk or
       if we uee a reentrant function we need to plan on not
       reseting this value since we might need the last few bytes
       of the cache used in previous call in cases */

    while (!io_eof_stream(dmx)){
       const unsigned char *byte = cache_buffer - 8;
       const uint64_t mask = OVVCDMX_IO_BUFF_MASK;
       int nb_stc = 0;
       int nb_epb = 0;
       int read_in_buf = 0;

       /* request new read from io layer */
       read_in_buf += io_read_stream(dmx);

       if (!read_in_buf) {
       /*FIXME check the case of eof is aligned with buffer*/
           goto last_chunk;
       }

       nb_chunks += read_in_buf;

       /* WARNING We need to be careful on endianness here if we plan
          to use bigger read sizes */
       do {
           /* Mask is cosmetic here only to signal the
              use of a circular buffer in case we want to
              change default cachesize into something more dynamic*/

           /*TODO bintricks for start code detection */

           /*FIXME we will actually loop over this more than once even if a start
             code has been detected. This is a bit inefficient */
           if (byte[(byte_pos) & mask] == 0) {
               int stc_or_epb;
               stc_or_epb = check_stc_or_epb(&byte[(byte_pos) & mask]);
               if (stc_or_epb < 0) {
                   return -1;
               }

               if (stc_or_epb) {
                   /* TODO handle what is to be done with it here */
                   dmx_process_elem(dmx, byte, byte_pos, stc_or_epb);
                   nb_stc += stc_or_epb == 1;
                   nb_epb += stc_or_epb == 2;
               }
           }
       } while ((++byte_pos) & mask);

       /* FIXME copy more bytes */
       /* copy last_two byte to padding area */
       printf("num_stc %d\n", nb_stc);
       printf("num_prev %d\n", nb_epb);
       nb_nalus += nb_stc;
    }

last_chunk:
    if (io_error_stream(dmx)) {
        return -1;
    } else {

        const uint64_t mask = OVVCDMX_IO_BUFF_MASK;
        /* we do not check return si error was already reported ?*/
        nb_bytes_last = io_tell_stream(dmx) & mask;

        if (io_eof_stream(dmx) && (nb_bytes_last > 0)) {
            int nb_stc = 0;
            int nb_epb = 0;
            uint8_t *byte = cache_buffer - 8;
            /*FIXME write 0xFF to prevent returning stc or emu
            at the end of the buffer*/

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
                        dmx_process_elem(dmx, byte, byte_pos, stc_or_epb);
                        nb_stc += stc_or_epb == 1;
                        nb_epb += stc_or_epb == 2;
                    }
                }
            } while (((++byte_pos) & mask) <= nb_bytes_last);
            printf("num_stc last %d\n", nb_stc);
            printf("num_prev last %d\n", nb_epb);
            printf("EOF reached\n");
            nb_nalus += nb_stc;
        }
    }

    printf("Num bytes read %ld\n", (nb_chunks << 16) + nb_bytes_last);
    printf("byte_pos %ld\n", (nb_chunks << 16) + nb_bytes_last);
    printf("Num NALU read %d\n", nb_nalus);

    return 1;
}

/* We use size_t trying to mimic fread return type */
static size_t
io_read_stream(OVVCDmx *const dmx)
{
    const size_t i_buff_size = OVVCDMX_IO_BUFF_SIZE;
    FILE *fstream = dmx->fstream;
    uint8_t *cache_start = dmx->cache_buffer.bytestream;
    uint8_t *cache_end = cache_start + i_buff_size;
    size_t read_in_buf;

    /* FIXME this might depend on the demux maybe this should
       be done somewhere else this force cache buffer */
    memcpy(cache_start - 8, cache_end - 8, sizeof(*cache_start) * 8);

    read_in_buf = fread(cache_start, i_buff_size, 1, fstream);

    return read_in_buf;
}

/* FIXME decide of an EOF value check what is usual*/
static int
io_eof_stream(OVVCDmx *const dmx)
{
    FILE *fstream = dmx->fstream;

    return feof(fstream);
}

static int
io_error_stream(OVVCDmx *const dmx)
{
    FILE *fstream = dmx->fstream;

    return ferror(fstream);
}

static long int
io_tell_stream(OVVCDmx *const dmx)
{
    FILE *fstream = dmx->fstream;

    return ftell(fstream);
}

/* Check for start_code or emulation prevention byte
 * returns:
 *    0 if no start code was found
      1 if a valid start code was found
      2 if a valid zero byte emulation prevention was found
      FIXME: should we return 4 if byte2 was zero ? (might be useful to
      abort detection on filling data or detect a new sequence without EOS
      or EOB )
      negative value if something invalid emulation prevention
      or start code with an invalid NAL header
 */
static inline int
check_stc_or_epb(const uint8_t *byte)
{
    /* we consider first byte has already been detected ?*/
    uint8_t byte1 = byte[1];
    uint8_t byte2 = byte[2];
    uint8_t byte3 = byte[3];

    if ((byte1 == 0) && !(byte2 & (~0x3))) {
        if (byte2 == 0x01) {
            uint8_t fbdn_zbit = byte3 & 0x01;
            if (fbdn_zbit) {
                goto invalid_data;
            }
            return 1;
        }

        /* Note this is only useful when removing rbsp
           and for probing*/
        if (byte2 == 0x03) {
            uint8_t invalid_emu = byte3 & (~0x3);
            if (invalid_emu) {
                goto invalid_data;
            }
            return 2;
        }
    }

    return 0;

invalid_data:
    /*TODO detect log what went wrong and return correct value
      Note this could be handled at a higher level based on return
      in case of failure we do not need to be fast*/
    return -1;
}

static int
dmx_process_elem(OVVCDmx *const dmx, const uint8_t *const bytestream,
        uint64_t byte_pos, int stc_or_epb)
{
    /* TODO decide what to do with this this is supposed to become
       an internal function this could be used as a callback to do
       various things according to the demuxer context*/
    const uint64_t mask = OVVCDMX_IO_BUFF_MASK;
    if (stc_or_epb == 1) {
        uint8_t byte4 = bytestream[((byte_pos) & mask) + 4];
        uint8_t nalu_type = (byte4 >> 3) & 0x1F;
        /* Start code */
        printf("STC at pos : %ld, NAL :%d\n", byte_pos - 8, nalu_type);
    } else {
        /* Emulation zero byte prevention */
        printf("EPB at pos : %ld\n", byte_pos - 8);
    }

    return 1;
}

