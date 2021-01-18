#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <stdint.h>

#include "libovvcutils/ovvcutils.h"
#include "libovvcutils/ovmem.h"
#include "libovvcutils/mempool.h"
#include "libovvcutils/mempool_internal.h"

#include "ovvcdmx.h"
#include "ovio.h"
#include "ovannexb.h"
#include "ovunits.h"


/* Use Fixed size of demux read cache buffer to 64K */

#define OVVCDMX_IO_BUFF_SIZE (1 << 16)

#define OVVCDMX_IO_BUFF_MASK ((1 << 16) - 1)

#define OVRBSP_CACHE_SIZE (1 << 16)

#define OVEPB_CACHE_SIZE (16 * sizeof(uint32_t))

static const char *const demux_name = "Open VVC Annex B demuxer";

enum RBSPSegmentDelimiter
{
    ANNEXB_STC = 1,
    ANNEXB_EPB = 2,
    END_OF_CACHE = 3
};

struct RBSPSegment
{
    /* position of segment first byte in reader cache */
    uint32_t start;

    /* position of segment last byte in cache */
    uint32_t end;
};

struct RBSPCacheData
{
    /* Cache buffer used to catenate RBSP chunks while
       extracting RBSP_data.
       Its size is initialised at 64kB and will grow
       to the max RBSP size encountered in the stream */
    uint8_t *start;

    uint8_t *end;

    /* Size of rbsp cache buffer used to check
       if we need to grow the cahce buffer when
       appending a chunk of RBSP buffer to the cache */
    size_t cache_size;

    /* Size of current catenate rbsp data
       this correspond to the last position where
       a chunk of RBSP data can be appended to the cache*/
    size_t rbsp_size;
};

struct EPBCacheInfo
{
    /* A table containing the locations of Emulation Prevention
       Bytes encountered in current NAL Unit RBSP
       this table is allocated to 16 * sizeof(*epb_pos) and will
       grow by the same amount of memory each time it is needed */
    uint32_t *epb_pos;

    /* Size in bytes of EPB encountered.
       Since encountering EBP is unlikely in compressed data
       the table will only grow of 16 * sizeof(*epb_pos) */
    size_t cache_size;

    /* The number of Emulation Prevention Bytes encountered
       int the current SODB while extracting the current
       NAL Unit RBSP Data*/
    int nb_epb;
};

struct NALUnitListElem
{
    struct NALUnitListElem *prev_nalu;
    struct NALUnitListElem *next_nalu;
    OVNALUnit nalu;
    struct {
        MemPoolElem *pool_ref;
    } private;
};

struct NALUnitsList
{
    struct NALUnitListElem *first_nalu;
    struct NALUnitListElem *last_nalu;

    int nb_nalus;
};

struct ReaderCache
{
    /* Pointer to io_cached buffer */
    const uint8_t *data_start;

    /* Cursor position in cache */
    uint32_t first_pos;

    /* FIXME the nb_chunk_read might overflow on really
     * long streams when the decoder is used for a long time
     * this will not affect decoding process but might result
     * in wrong informations reporting about the demuxer status
     */

    /* Number of 64 kB chunks already processed by
     * the demuxer */
    uint64_t nb_chunk_read;
};

struct OVVCDmx
{
    const char *name;

    FILE *fstream;
    /*OVReadBuff cache_buffer;*/
    int nb_stc;
    int nb_epb;

    /* Points to a read only IO context */
    OVIOStream *io_str;

    /* Information on current Stream Of Data Bytes (SODB)
     */
    struct ReaderCache cache_ctx;

    /* Cache used for RBSP extraction */
    struct RBSPCacheData rbsp_ctx;

    /* Cache uses to store Emulation Prevention Bytes (EPB)
       when extracting RBSP Data from current NAL Unit */
    struct EPBCacheInfo epb_info;

    /* A chained list containing current Acces Unit (AU)
       or Picture Unit (PU) NAL Units */
    struct NALUnitsList nalu_list;

    /* Memory pool for NALUListElem */
    MemPool *nalu_elem_pool;

    /* Demuxer options to be passed at init */
    struct{
        int val;
    }options;
};

static int extract_cache_segments(OVVCDmx *const dmx, struct ReaderCache *const cache_ctx);

static int process_last_chunk(OVVCDmx *const dmx, const uint8_t *byte, uint64_t byte_pos,
                              int nb_bytes_last);

static int init_rbsp_cache(struct RBSPCacheData *const rbsp_ctx);

/* Realloc rbsp_cache adding an extra 64KB to previously allocated size
   and copy previous content */
static int extend_rbsp_cache(struct RBSPCacheData *const rbsp_ctx);

static void free_rbsp_cache(struct RBSPCacheData *const rbsp_ctx);

static int init_epb_cache(struct EPBCacheInfo *const epb_info);

/* Realloc EPB locations caching adding an extra 16 elements to previously
 * allocated size and copy previous content
 */
static int extend_epb_cache(struct EPBCacheInfo *const epb_info);

static void free_epb_cache(struct EPBCacheInfo *const epb_info);

static void free_nalu_list(struct NALUnitsList *list);

static int refill_reader_cache(struct ReaderCache *const cache_ctx,
                               OVIOStream *const io_str);

int
ovdmx_init(OVVCDmx **vvcdmx)
{
    int ret;
    *vvcdmx = ov_mallocz(sizeof(**vvcdmx));

    if (*vvcdmx == NULL) return -1;

    (*vvcdmx)->name = demux_name;
    (*vvcdmx)->io_str = NULL;

    (*vvcdmx)->nalu_elem_pool = ovmempool_init(sizeof(struct NALUnitListElem));

    if ((*vvcdmx)->nalu_elem_pool == NULL) {
        goto fail_pool_init;
    }

    ret = init_rbsp_cache(&(*vvcdmx)->rbsp_ctx);
    if (ret < 0) {
        goto fail_rbsp_cache;
    }

    ret = init_epb_cache(&(*vvcdmx)->epb_info);
    if (ret < 0) {
        goto fail_epb_cache;
    }

    return 0;

fail_epb_cache:
    free_rbsp_cache(&(*vvcdmx)->rbsp_ctx);

fail_rbsp_cache:
    ovmempool_uninit(&(*vvcdmx)->nalu_elem_pool);

fail_pool_init:
    ov_freep(vvcdmx);

    return -1;
}

int
ovdmx_close(OVVCDmx *vvcdmx)
{
    int not_dmx = 0;
    if (vvcdmx != NULL) {

        not_dmx = vvcdmx->name != demux_name;

        if (not_dmx) goto fail;

        ovdmx_detach_stream(vvcdmx);

        ovmempool_uninit(&vvcdmx->nalu_elem_pool);

        free_rbsp_cache(&vvcdmx->rbsp_ctx);

        free_epb_cache(&vvcdmx->epb_info);

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
    int ret = 0;

    /* FiXME is this check necessary this function should not
       be called if stream is not allocated.
       Maybe we should open file ourselves / and use a wrapper around
       I/Os */
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

    /* TODO init SODB ctx by first read */

    if (!ovio_stream_eof(dmx->io_str)) {
        int read_in_buf;
        struct ReaderCache *const cache_ctx = &dmx->cache_ctx;

        read_in_buf += ovio_stream_read(&cache_ctx->data_start, OVVCDMX_IO_BUFF_SIZE,
                                        dmx->io_str);
        cache_ctx->data_start -= 8;
        /* pos is init at 8 so we do not overread first data chunk */
        cache_ctx->first_pos = 0;

        /* FIXME Process first chunk of data ? */
        ret = extract_cache_segments(dmx, cache_ctx);

        if (!read_in_buf) {
            /* TODO error handling if end of file is encountered on first read */
        }
        cache_ctx->nb_chunk_read = 1;
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
    /* FIXME ReaderCache  should be reset */
}

static int
refill_reader_cache(struct ReaderCache *const cache_ctx, OVIOStream *const io_str)
{
    int read_in_buf;

    read_in_buf = ovio_stream_read(&cache_ctx->data_start, OVVCDMX_IO_BUFF_SIZE,
                                   io_str);

    cache_ctx->data_start -= 8;

    cache_ctx->first_pos   = 0;

    cache_ctx->nb_chunk_read += read_in_buf;

    if (!read_in_buf) {
        return -1;
    }

    return 0;
}

static int
extract_nal_unit(OVVCDmx *dmx)
{
    return 0;
}

static uint8_t
is_access_unit_delimiter(struct NALUnitListElem *elem)
{
    /* FIXME add other rules based on NAL Unit types ordering
     * Since some AU delimitations rules involve POC computation
     * this require reading until slice header
     */
    return elem->nalu.type == OVNALU_AUD;
}

static struct NALUnitListElem *pop_nalu_elem(struct NALUnitsList *list)
{
    struct NALUnitListElem *elem = NULL;
    elem = list->first_nalu;
    if (elem) {
        list->first_nalu = elem->next_nalu;
        if (elem->next_nalu) {
            elem->next_nalu->prev_nalu = NULL;
        }
    } else {
        /*FIXME ensure list is emptied */
        list->last_nalu = NULL;
    }
    elem->prev_nalu = NULL;
    elem->next_nalu = NULL;

    return elem;
}

static void
push_nalu_elem(struct NALUnitsList *list, struct NALUnitListElem *elem)
{
    if (list->last_nalu) {
        elem->prev_nalu = list->last_nalu;
        list->last_nalu = elem;
    } else {
       /*FIXME check no first elem */
       list->first_nalu = elem;
       list->last_nalu = elem;
    }
}

static int
extract_access_unit(OVVCDmx *const dmx)
{
    struct NALUnitsList *nalu_list = &dmx->nalu_list;
    struct NALUnitListElem *current_nalu = nalu_list->first_nalu;
    struct NALUnitsList pending_nalu_list = {0};
    uint8_t au_end_found = 0;
    do {
        if (!current_nalu) {
            int ret;
            struct ReaderCache *const cache_ctx = &dmx->cache_ctx;
            current_nalu = nalu_list->last_nalu;
            /* No pending NALU in dmx try to extract NALU units from next
               chunk */
            ret = refill_reader_cache(cache_ctx, dmx->io_str);
            if (ret < 0) {
                goto last_chunk;
            }

            ret = extract_cache_segments(dmx, cache_ctx);

            if (!current_nalu) {
                current_nalu = nalu_list->first_nalu;
            }
        }

        if (current_nalu) {

            do {
               /* Move NALU from NALU list to pending list */
               push_nalu_elem(&pending_nalu_list, current_nalu);

               current_nalu = pop_nalu_elem(nalu_list);

            } while (current_nalu && !is_access_unit_delimiter(current_nalu));

            au_end_found = current_nalu != NULL;
        }
    } while (!au_end_found);

    free_nalu_list(&pending_nalu_list);

    /* TODO NALUList to packet */

    return 0;

last_chunk:
    {
    int nb_chunks = dmx->cache_ctx.nb_chunk_read;
    int nb_bytes_last = 0;
    if (ovio_stream_error(dmx->io_str)) {
        return -1;
    } else {
       const uint64_t mask = OVVCDMX_IO_BUFF_MASK;
       int ret;

        /* we do not check return si error was already reported ?*/

        nb_bytes_last = ovio_stream_tell(dmx->io_str) & mask;

        if (ovio_stream_eof(dmx->io_str) && (nb_bytes_last > 0)) {

            const uint8_t *byte = dmx->cache_ctx.data_start;

            ret = process_last_chunk(dmx, byte, 0, nb_bytes_last);
            if (ret < 0) {
                goto readfail;
            }

            printf("num_stc last %d\n", dmx->nb_stc);
            printf("num_prev last %d\n", dmx->nb_epb);
            printf("EOF reached\n");
        }
    }

    printf("Num bytes read %ld\n", (nb_chunks << 16) + nb_bytes_last);
    printf("byte_pos %ld\n", (nb_chunks << 16) + nb_bytes_last);

    return 0;
    }

readfail:
    /* free current NALU memory + return error
    */
    printf("FAILED NAL read!\n");
    return -1;
}

static void free_nalu_elem(struct NALUnitListElem *nalu_elem);

static void
free_nalu_list(struct NALUnitsList *list)
{
    struct NALUnitListElem *elem = list->first_nalu;
    while (elem) {
        struct NALUnitListElem *to_free = elem;
        elem = elem->next_nalu;
        free_nalu_elem(to_free);
    }
    list->first_nalu = NULL;
    list->last_nalu = NULL;
}

int
ovdmx_extract_picture_unit(OVVCDmx *const dmx, OVPictureUnit **dst_pu)
{
    int ret;
    OVPictureUnit *pu = ov_mallocz(sizeof(*pu));
    if (!pu) {
        *dst_pu = NULL;
        return -1;
    }

    ret = extract_access_unit(dmx);
    if (ret < 0) {
        ov_free(pu);
        return -1;
    }

    *dst_pu = pu;

    return 0;
}

/*
   FIXME there might be a need to find a fast way to remove obvious
   zero filling data between RBPS and next NALU
   FIXME we should move thi part to other file for annexb demux
   if we plan to support more demuxers format*/
int
ovdmx_read_stream(OVVCDmx *const dmx)
{
    uint64_t nb_chunks = 1;
    uint64_t byte_pos = 0;
    long int nb_bytes_last = 0;
    int nb_nalus = 0;
    OVIOStream *const io_str = dmx->io_str;
    struct ReaderCache *const cache_ctx = &dmx->cache_ctx;
    const uint8_t **io_cache = &cache_ctx->data_start;

    if (io_str == NULL) {
       return -1;
    }

    /* Note if we need to read more than one chunk or
       if we uee a reentrant function we need to plan on not
       reseting this value since we might need the last few bytes
       of the cache used in previous call in cases */

    if (!ovio_stream_eof(io_str)) {

        while (!ovio_stream_eof(io_str)) {
            extract_access_unit(dmx);

            free_nalu_list(&dmx->nalu_list);
        }

    } else {
        goto last_chunk;
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

            const uint8_t *byte = *io_cache - 8;

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

static struct NALUnitListElem *
create_nalu_elem(OVVCDmx *const dmx)
{
    MemPool *const mempool = dmx->nalu_elem_pool;
    MemPoolElem *elem = ovmempool_popelem(mempool);
    struct NALUnitListElem *nalu_elem = NULL;
    if (!elem) {
        return NULL;
    }
    nalu_elem = elem->data;
    nalu_elem->private.pool_ref = elem;

    nalu_elem->nalu.rbsp_data = NULL;
    nalu_elem->nalu.rbsp_size = 0;

    nalu_elem->nalu.epb_pos = NULL;
    nalu_elem->nalu.nb_epb = 0;

    return nalu_elem;
}

static void
free_nalu_elem(struct NALUnitListElem *nalu_elem)
{
    /* TODO unref NALU instead of free */
    ov_freep(&nalu_elem->nalu.rbsp_data);

    if (nalu_elem->nalu.epb_pos){
        ov_freep(&nalu_elem->nalu.epb_pos);
    }

    /* clean up before returning to pool */
    nalu_elem->nalu.nb_epb    = 0;
    nalu_elem->nalu.rbsp_size = 0;
    ovmempool_pushelem(nalu_elem->private.pool_ref);
}

static void
append_nalu_elem(struct NALUnitsList *const list, struct NALUnitListElem *elem)
{
    if (!list->last_nalu) {
        /* list is empty */
        elem->prev_nalu = NULL;
        list->first_nalu = elem;
    } else {
        elem->prev_nalu = list->last_nalu;
        list->last_nalu->next_nalu = elem;
    }
    elem->next_nalu = NULL;
    list->last_nalu = elem;
}

static int
append_rbsp_segment_to_cache(struct ReaderCache *const cache_ctx,
                             struct RBSPCacheData *rbsp_cache,
                             struct RBSPSegment *sgmt_ctx)
{
    int sgmt_size = sgmt_ctx->end - sgmt_ctx->start;
    if (rbsp_cache->cache_size < rbsp_cache->rbsp_size + sgmt_size) {
        int ret;
         ret = extend_rbsp_cache(rbsp_cache);
         if (ret < 0) {
             return -1;
         }
    }

    memcpy(rbsp_cache->end, &cache_ctx->data_start[sgmt_ctx->start], sgmt_size);

    rbsp_cache->end       += sgmt_size;
    rbsp_cache->rbsp_size += sgmt_size;

    sgmt_ctx->start = sgmt_ctx->end + 3;

    return 0;
}

static void
empty_rbsp_cache(struct RBSPCacheData *rbsp_cache)
{
    rbsp_cache->end        = rbsp_cache->start;
    rbsp_cache->rbsp_size = 0;
}

static int
process_start_code(OVVCDmx *const dmx, struct ReaderCache *const cache_ctx,
                   uint64_t byte_pos, struct RBSPSegment *sgmt_ctx)
{
    const uint64_t mask = OVVCDMX_IO_BUFF_MASK;
    const uint8_t *bytestream = &cache_ctx->data_start[byte_pos & mask];
    struct NALUnitsList *nalu_list = &dmx->nalu_list;
    struct NALUnitListElem *nalu_elem = create_nalu_elem(dmx);

    enum OVNALUType nalu_type = (bytestream[4] >> 3) & 0x1F;

    if (!nalu_elem) {
        printf("NALU alloc fail\n");
        return -1;
    }

    sgmt_ctx->end = byte_pos;
    append_rbsp_segment_to_cache(cache_ctx, &dmx->rbsp_ctx, sgmt_ctx);
    sgmt_ctx->start = byte_pos + 3;

    if (nalu_list->last_nalu) {
        /*TODO attach rbsp_data to nalu*/
        uint8_t *rbsp_data = ov_malloc(dmx->rbsp_ctx.rbsp_size);
        if (!rbsp_data) {
            free_nalu_elem(nalu_elem);
            return -1;
        }

        if (dmx->epb_info.nb_epb) {
            uint32_t *epb_pos = NULL;
            epb_pos = ov_malloc(dmx->epb_info.nb_epb * sizeof(*epb_pos));
            if (!epb_pos) {
                free_nalu_elem(nalu_elem);
                ov_free(rbsp_data);
                return -1;
            }
            memcpy(epb_pos, dmx->epb_info.epb_pos, dmx->epb_info.nb_epb * sizeof(*epb_pos));
            nalu_list->last_nalu->nalu.epb_pos = epb_pos;
            nalu_list->last_nalu->nalu.nb_epb = dmx->epb_info.nb_epb;
        }
        dmx->epb_info.nb_epb = 0;

        nalu_list->last_nalu->nalu.rbsp_data = rbsp_data;
        nalu_list->last_nalu->nalu.rbsp_size = dmx->rbsp_ctx.rbsp_size;

        memcpy(rbsp_data, dmx->rbsp_ctx.start, dmx->rbsp_ctx.rbsp_size);

        empty_rbsp_cache(&dmx->rbsp_ctx);
    }

    nalu_elem->nalu.type = nalu_type;

    append_nalu_elem(nalu_list, nalu_elem);

    dmx->nb_stc++;

    printf("STC at pos : %ld, NAL :%d\n", byte_pos - 8, nalu_type);
    return 0;
}

static int
process_emulation_prevention_byte(OVVCDmx *const dmx, struct ReaderCache *const cache_ctx,
                                  uint64_t byte_pos, struct RBSPSegment *sgmt_ctx)
{
    struct EPBCacheInfo *const epb_info = &dmx->epb_info;
    sgmt_ctx->end = byte_pos + 2;
    append_rbsp_segment_to_cache(cache_ctx, &dmx->rbsp_ctx, sgmt_ctx);
    sgmt_ctx->start = sgmt_ctx->end + 1;

    if (epb_info->nb_epb + 1 > (epb_info->cache_size)/sizeof(*epb_info->epb_pos)) {
        int ret = extend_epb_cache(epb_info);
        if (ret < 0) {
            printf("ERROR extending cache\n");
            return -1;
        }
    }

    epb_info->epb_pos[epb_info->nb_epb++];

    dmx->nb_epb++;

    printf("EPB at pos : %ld\n", byte_pos - 8);
    return 0;
}

/**
 * returns: -1 Invalid data
 *          byte_pos in chunk if stc
 *          0 if nothing found and needs a new read
 */
/* WARNING We need to be careful on endianness here if we plan
   to use bigger read sizes */
static int
extract_cache_segments(OVVCDmx *const dmx, struct ReaderCache *const cache_ctx)
{
    const uint64_t mask = OVVCDMX_IO_BUFF_MASK;
    const uint8_t *byte = &cache_ctx->data_start[cache_ctx->first_pos];
    uint32_t byte_pos = cache_ctx->first_pos & mask;
    uint8_t end_of_cache;
    struct RBSPSegment sgmt_ctx = {0};

    do {
        const uint8_t *bytestream = &byte[byte_pos & mask];

        /*FIXME we will actually loop over this more than once even if a start
          code has been detected. This is a bit inefficient */
        /*TODO bintricks for two fast zero bytes detection*/
        if (*bytestream == 0) {
            int ret;
            ret = ovannexb_check_stc_or_epb(bytestream);
            if (ret < 0) {
                printf("Invalid\n");
                ret = -1;
            }

            if (ret) {
                enum RBSPSegmentDelimiter dlm = ret;

                switch (dlm) {
                case ANNEXB_STC:
                    ret = process_start_code(dmx, cache_ctx, byte_pos, &sgmt_ctx);
                    break;
                case ANNEXB_EPB:
                    ret = process_emulation_prevention_byte(dmx, cache_ctx, byte_pos, &sgmt_ctx);
                    break;
                }

                if (ret < 0) {
                    return ret;
                }
            }
        }
        end_of_cache = ((++byte_pos) & mask) == 0;
    } while (!end_of_cache);

    /* End of cache */
    sgmt_ctx.end = byte_pos;
    append_rbsp_segment_to_cache(cache_ctx, &dmx->rbsp_ctx, &sgmt_ctx);

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
            stc_or_epb = ovannexb_check_stc_or_epb(&byte[(byte_pos) & mask]);

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

static int
init_rbsp_cache(struct RBSPCacheData *const rbsp_ctx)
{
    rbsp_ctx->start = ov_mallocz(OVRBSP_CACHE_SIZE);
    if (rbsp_ctx->start == NULL) {
        return -1;
    }

    rbsp_ctx->end = rbsp_ctx->start;
    rbsp_ctx->cache_size = OVRBSP_CACHE_SIZE;

    return 0;
}

static void
free_rbsp_cache(struct RBSPCacheData *const rbsp_ctx)
{
    ov_freep(&rbsp_ctx->start);
}

static int
extend_rbsp_cache(struct RBSPCacheData *const rbsp_ctx)
{
    uint8_t *old_cache = rbsp_ctx->start;
    uint8_t *new_cache;
    size_t new_size = rbsp_ctx->cache_size + OVRBSP_CACHE_SIZE;
    new_cache = ov_malloc(new_size);
    if (!new_cache) {
        return -1;
    }

    memcpy(new_cache, old_cache, rbsp_ctx->rbsp_size);

    ov_free(old_cache);

    rbsp_ctx->start = new_cache;
    rbsp_ctx->end = rbsp_ctx->start + rbsp_ctx->rbsp_size;
    rbsp_ctx->cache_size = new_size;
    return 0;
}

static int
init_epb_cache(struct EPBCacheInfo *const epb_info)
{
    epb_info->epb_pos = ov_mallocz(OVEPB_CACHE_SIZE);
    if (epb_info->epb_pos == NULL) {
        return -1;
    }

    epb_info->cache_size = OVEPB_CACHE_SIZE;

    return 0;
}

/* Realloc EPB locations caching adding an extra 16 elements to previously
 * allocated size and copy previous content
 */
static int
extend_epb_cache(struct EPBCacheInfo *const epb_info)
{
    uint32_t *old_cache = epb_info->epb_pos;
    uint32_t *new_cache;
    size_t new_size = epb_info->cache_size + OVEPB_CACHE_SIZE;
    new_cache = ov_malloc(new_size);
    if (!new_cache) {
        return -1;
    }

    memcpy(new_cache, old_cache, epb_info->nb_epb * sizeof(*epb_info->epb_pos));

    ov_free(old_cache);

    epb_info->epb_pos    = new_cache;
    epb_info->cache_size = new_size;
    return 0;
}

static void
free_epb_cache(struct EPBCacheInfo *const epb_info)
{
    ov_freep(&epb_info->epb_pos);
}
