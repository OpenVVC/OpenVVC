/**
 *
 *   OpenVVC is open-source real time software decoder compliant with the 
 *   ITU-T H.266- MPEG-I - Part 3 VVC standard. OpenVVC is developed from 
 *   scratch in C as a library that provides consumers with real time and
 *   energy-aware decoding capabilities under different OS including MAC OS,
 *   Windows, Linux and Android targeting low energy real-time decoding of
 *   4K VVC videos on Intel x86 and ARM platforms.
 * 
 *   Copyright (C) 2020-2022  IETR-INSA Rennes :
 *   
 *   Pierre-Loup CABARAT
 *   Wassim HAMIDOUCHE
 *   Guillaume GAUTIER
 *   Thomas AMESTOY
 *
 *   This library is free software; you can redistribute it and/or
 *   modify it under the terms of the GNU Lesser General Public
 *   License as published by the Free Software Foundation; either
 *   version 2.1 of the License, or (at your option) any later version.
 *
 *   This library is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *   Lesser General Public License for more details.
 *
 *   You should have received a copy of the GNU Lesser General Public
 *   License along with this library; if not, write to the Free Software
 *   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301
 *   USA
 * 
 **/

#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <stdint.h>

#include "ovutils.h"
#include "overror.h"
#include "ovmem.h"
#include "mempool.h"
#include "mempool_internal.h"

#include "ovdmx.h"
#include "ovio.h"
#include "ovannexb.h"
#include "ovunits.h"


#define OVRBSP_CACHE_SIZE (1 << 16)

#define OVEPB_CACHE_SIZE (16 * sizeof(uint32_t))

#define OV_RBSP_PADDING 8

enum DMXReturn
{
    OV_INVALID_DATA = -1,
    OV_ENOMEM = -2,
};

static const char *const demux_name = "Open VVC Annex B demuxer";

enum RBSPSegmentDelimiter
{
    /* NAL Unit Start Code 0x000001 */
    ANNEXB_STC = 1,
    /* NAL Unit Start Code 0x000003 */
    ANNEXB_EPB = 2,
    /* Reached end of cache or EOF */
    END_OF_CACHE = 3
};

struct RBSPSegment
{
    /* position of segment first byte in reader cache */
    const uint8_t *start_p;

    /* position of segment last byte in cache */
    const uint8_t *end_p;
};

struct RBSPCacheData
{
    /* Cache buffer used to catenate RBSP chunks while extracting RBSP_data.
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

    const uint8_t *cache_start;
    const uint8_t *cache_end;

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

struct OVDemux
{
    const char *name;

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

    /* Current NAL Unit to be added to list when end is found */
    struct NALUnitListElem *nalu_pending;

    /* Memory pool for NALUListElem */
    MemPool *nalu_elem_pool;

    uint8_t eof;

    /* Demuxer options to be passed at init */
    struct{
        int val;
    }options;
};

static int extract_cache_segments(OVDemux *const dmx, struct ReaderCache *const cache_ctx);

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

static void append_nalu_elem(struct NALUnitsList *const list, struct NALUnitListElem *elem);

static void free_nalu_elem(struct NALUnitListElem *nalu_elem);

int
ovdmx_init(OVDemux **dmx_p)
{
    int ret = 0;
    OVDemux *dmx;

    dmx = ov_mallocz(sizeof(*dmx));

    if (!dmx) return OV_ENOMEM;

    *dmx_p = dmx;

    dmx->name = demux_name;
    dmx->io_str = NULL;

    dmx->nalu_elem_pool = ovmempool_init(sizeof(struct NALUnitListElem));

    if (dmx->nalu_elem_pool == NULL) {
        goto fail_pool_init;
    }

    ret = init_rbsp_cache(&dmx->rbsp_ctx);
    if (ret < 0) {
        goto fail_rbsp_cache;
    }

    ret = init_epb_cache(&dmx->epb_info);
    if (ret < 0) {
        goto fail_epb_cache;
    }

    return 0;

fail_epb_cache:
    free_rbsp_cache(&dmx->rbsp_ctx);

fail_rbsp_cache:
    ovmempool_uninit(&dmx->nalu_elem_pool);

fail_pool_init:
    ov_freep(dmx_p);

    return ret;
}

int
ovdmx_close(OVDemux *dmx)
{
    int not_dmx = 0;
    if (dmx != NULL) {

        not_dmx = dmx->name != demux_name;

        if (not_dmx) goto fail;

        ovdmx_detach_stream(dmx);

        ovmempool_uninit(&dmx->nalu_elem_pool);

        free_rbsp_cache(&dmx->rbsp_ctx);

        free_epb_cache(&dmx->epb_info);

        free_nalu_list(&dmx->nalu_list);

        if (dmx->nalu_pending) {
            free_nalu_elem(dmx->nalu_pending);
        }

        ov_free(dmx);

        return 0;
    }

fail:
    ov_log(dmx, OVLOG_ERROR, "Trying to close a something not a demuxer.\n");
    return -1;
}

int
ovdmx_attach_stream(OVDemux *const dmx, OVIO *io)
{
    int ret = 0;

    /* FiXME is this check necessary this function should not
       be called if stream is not allocated.
       Maybe we should open file ourselves / and use a wrapper around
       I/Os */
    if (io == NULL) {
        ov_log(dmx, OVLOG_ERROR, "No stream to attach.\n");
        return OVVC_EINDATA;
    }

    /* TODO distinguish init and open / attach */
    dmx->io_str = ovio_stream_open(io);
    if (dmx->io_str == NULL) {
        ov_log(dmx, OVLOG_ERROR, "Failed to open stream.\n");
        return OVVC_EINDATA;
    }

    /* Initialise reader cache by first read */
    if (!ovio_stream_eof(dmx->io_str)) {
        struct ReaderCache *const cache_ctx = &dmx->cache_ctx;
        int read_in_buf;

        read_in_buf = ovio_stream_read(&cache_ctx->data_start, dmx->io_str);
        cache_ctx->first_pos = 0;
        cache_ctx->cache_start = cache_ctx->data_start;

        cache_ctx->cache_end = cache_ctx->cache_start + read_in_buf;

        if (read_in_buf < io->size) {
            dmx->eof = 1;
        }
        else {
            /* Buffer end is set to size minus 8 so we do not overread first data chunk */
            cache_ctx->cache_end -= 8;
        }

        /* FIXME Process first chunk of data ? */
        ret = extract_cache_segments(dmx, cache_ctx);

        cache_ctx->nb_chunk_read = 1;
    }

    return ret;
}

/*FIXME share bytestream mem with OVIOStream so we avoid copying from one to
  another*/
void
ovdmx_detach_stream(OVDemux *const dmx)
{
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
    read_in_buf = ovio_stream_read(&cache_ctx->data_start, io_str);
    cache_ctx->data_start -= 8;

    cache_ctx->cache_start = cache_ctx->data_start;
    cache_ctx->cache_end   = cache_ctx->data_start + read_in_buf;

    cache_ctx->nb_chunk_read += read_in_buf;

    if (read_in_buf != ovio_stream_buff_size(io_str)) {
        cache_ctx->cache_end += 8;
        return 1;
    }

    return 0;
}

#if 0
static uint8_t
is_access_unit_delimiter(struct NALUnitListElem *elem)
{
    /* FIXME add other rules based on NAL Unit types ordering
     * Since some AU delimitations rules involve POC computation
     * this require reading until slice header
     */
    return elem->nalu.type == OVNALU_AUD || elem->nalu.type == OVNALU_PPS;
}
#endif

static struct NALUnitListElem *pop_nalu_elem(struct NALUnitsList *list)
{
    struct NALUnitListElem *elem = NULL;
    elem = list->first_nalu;
    if (elem) {
        list->first_nalu = elem->next_nalu;
        if (elem->next_nalu) {
            elem->next_nalu->prev_nalu = NULL;
        }
        elem->prev_nalu = NULL;
        elem->next_nalu = NULL;
    } else {
        /*FIXME ensure list is emptied */
        list->last_nalu = NULL;
    }

    return elem;
}

static int
extract_nal_unit(OVDemux *const dmx, struct NALUnitsList *const dst_list)
{
    struct NALUnitsList *nalu_list = &dmx->nalu_list;
    struct NALUnitListElem *current_nalu = pop_nalu_elem(nalu_list);

    do {
        if (!current_nalu && !dmx->eof) {
            struct ReaderCache *const cache_ctx = &dmx->cache_ctx;

            /* FIXME error handling from demux + use return values */
            dmx->eof = refill_reader_cache(cache_ctx, dmx->io_str);

            extract_cache_segments(dmx, cache_ctx);

            current_nalu = pop_nalu_elem(nalu_list);
        }

        if (current_nalu) {
            append_nalu_elem(dst_list, current_nalu);
        }

    } while (current_nalu == NULL && !dmx->eof);

    return -(current_nalu == NULL && dmx->eof);
}

#if 0
static int
extract_access_unit(OVDemux *const dmx, struct NALUnitsList *const dst_list)
{
    struct NALUnitsList *nalu_list = &dmx->nalu_list;
    struct NALUnitListElem *current_nalu = pop_nalu_elem(nalu_list);
    uint8_t au_end_found = 0;
    do {
        if (!current_nalu) {
            /* No NALU in dmx try to extract NALU units from next
               chunk */
            struct ReaderCache *const cache_ctx = &dmx->cache_ctx;
            int ret = -1;

            if (!eof)
            ret = refill_reader_cache(cache_ctx, dmx->io_str);

            if (ret < 0) {
                eof = 1;
                #if 0
                goto last_chunk;
                #endif
            }

            ret = extract_cache_segments(dmx, cache_ctx);

            current_nalu = pop_nalu_elem(nalu_list);
        }

        if (current_nalu) {
            do {
                /* Move NALU from NALU list to pending list */
                append_nalu_elem(dst_list, current_nalu);
                printf("%d\n",current_nalu->nalu.type);

                current_nalu = pop_nalu_elem(nalu_list);

            } while (current_nalu && !is_access_unit_delimiter(current_nalu));

            au_end_found = current_nalu != NULL || (current_nalu == NULL && eof);
        }
    } while (!au_end_found && !eof);

    if (current_nalu)
        append_nalu_elem(dst_list, current_nalu);

    /* TODO NALUList to packet */

    return 0;

#if 0
readfail:
    /* free current NALU memory + return error
    */
    printf("FAILED NAL read!\n");
    return -1;
#endif
}
#endif

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

static void
move_nalu_elem_to_ovnalu(struct NALUnitListElem *lelem, OVNALUnit *nalu)
{
    memcpy(nalu, &lelem->nalu, sizeof(lelem->nalu)) ;
    memset(&lelem->nalu, 0, sizeof(lelem->nalu));
}

static int
count_list_nal_units(struct NALUnitsList *const src)
{
    /* FIXME Avoid counting list elements */
    int nb_nalus = 0;
    struct NALUnitListElem *lelem = src->first_nalu;
    while (lelem) {
        lelem = lelem->next_nalu;
        nb_nalus++;
    }
    return nb_nalus;
}

int
ovdmx_init_pu_from_list(OVPictureUnit **ovpu_p, struct NALUnitsList *const src)
{
    OVPictureUnit *ovpu;
    struct NALUnitListElem *lelem = src->first_nalu;
    int nb_nalus = count_list_nal_units(src);
    int ret = ovpu_init(ovpu_p, nb_nalus);
    int i;
    if (ret < 0) {
        return ret;
    }

    ovpu = *ovpu_p;

    for (i = 0; i < nb_nalus; i++) {
        ovpu->nalus[i] = ov_mallocz(sizeof(*ovpu->nalus[i]));
        if (!ovpu->nalus[i]) {
            goto fail_nalu_alloc;
        }
        move_nalu_elem_to_ovnalu(lelem, ovpu->nalus[i]);
        lelem = lelem->next_nalu;
    }

    ovpu->nb_nalus = nb_nalus;

    *ovpu_p = ovpu;

    return 0;

fail_nalu_alloc:
    ovpu_unref(ovpu_p);
    return OVVC_ENOMEM;
}

int
ovdmx_extract_picture_unit(OVDemux *const dmx, OVPictureUnit **dst_pu_p)
{
    int ret;
    struct NALUnitsList pending_nalu_list = {0};

    #if 0
    if (/*!dmx->eof &&*/ dmx->nalu_list.first_nalu) {
        ret = extract_access_unit(dmx, &pending_nalu_list);

        /* FIXME return */


        if (!dmx->eof && ret < 0) {
            ov_log(dmx, OVLOG_ERROR, "No valid Access Unit found \n");
            free_nalu_list(&pending_nalu_list);
            ov_free(pu);
            return ret;
        }
    } else {
        ov_free(pu);
        *dst_pu = NULL;
        return -1;
    }
    #else
    ret = extract_nal_unit(dmx, &pending_nalu_list);
    if (!dmx->eof && ret < 0) {
        ov_log(dmx, OVLOG_ERROR, "No valid Access Unit found \n");
        free_nalu_list(&pending_nalu_list);
        return ret;
    }
    #endif

    int ret2 = ovdmx_init_pu_from_list(dst_pu_p, &pending_nalu_list);
    if (ret2 < 0) {
        free_nalu_list(&pending_nalu_list);
        return ret2;
    }

    free_nalu_list(&pending_nalu_list);

    return ret;
}

static struct NALUnitListElem *
create_nalu_elem(OVDemux *const dmx)
{
    MemPool *const mempool = dmx->nalu_elem_pool;
    MemPoolElem *elem = ovmempool_popelem(mempool);

    struct NALUnitListElem *nalu_elem = NULL;
    if (!elem) {
        return NULL;
    }

    nalu_elem = elem->data;
    nalu_elem->private.pool_ref = elem;

    ov_nalu_init(&nalu_elem->nalu);

    return nalu_elem;
}

static void
free_nalu_elem(struct NALUnitListElem *nalu_elem)
{
    /* TODO unref NALU instead of free */
    if (nalu_elem->nalu.rbsp_data) {
        ov_freep(&nalu_elem->nalu.rbsp_data);
    }

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

/* FIXME remove unused cache_ctx */
static int
append_rbsp_segment_to_cache(struct ReaderCache *const cache_ctx,
                             struct RBSPCacheData *rbsp_cache,
                             const struct RBSPSegment *sgmt_ctx)
{
    ptrdiff_t sgmt_size = sgmt_ctx->end_p - sgmt_ctx->start_p;
    /* FIXME use an assert instead this is not supposed to happen */
    if (sgmt_size < 0) {
        ov_log(NULL, OVLOG_ERROR, "Invalid segment\n");
        return -1;
    }

    if (rbsp_cache->cache_size < rbsp_cache->rbsp_size + sgmt_size) {
        int ret;
         ret = extend_rbsp_cache(rbsp_cache);
         if (ret < 0) {
             return ret;
         }
    }

    memcpy(rbsp_cache->end, sgmt_ctx->start_p, (size_t)sgmt_size);

    rbsp_cache->end       += sgmt_size;
    rbsp_cache->rbsp_size += sgmt_size;

    return 0;
}

static void
empty_rbsp_cache(struct RBSPCacheData *rbsp_cache)
{
    rbsp_cache->end        = rbsp_cache->start;
    rbsp_cache->rbsp_size = 0;
}

static int
process_start_code(OVDemux *const dmx, struct ReaderCache *const cache_ctx,
                   const struct RBSPSegment *sgmt_ctx)
{
    const uint8_t *bytestream = sgmt_ctx->end_p;
    struct NALUnitsList *nalu_list = &dmx->nalu_list;
    struct NALUnitListElem *nalu_elem = create_nalu_elem(dmx);
    struct NALUnitListElem *nalu_pending = dmx->nalu_pending;

    enum OVNALUType nalu_type = (bytestream[4] >> 3) & 0x1F;

    if (!nalu_elem) {
        ov_log(dmx, OVLOG_ERROR, "Could not alloc NALU element\n");
        return OV_ENOMEM;
    }

    /* New NAL Unit start code found we end so we can process previous
     * NAL Unit data
     */
    append_rbsp_segment_to_cache(cache_ctx, &dmx->rbsp_ctx, sgmt_ctx);

    if (nalu_pending) {
        /* FIXME Using of mallocz is to prevent padding to be not zero */
        uint8_t *rbsp_data = ov_mallocz(dmx->rbsp_ctx.rbsp_size + OV_RBSP_PADDING);
        if (!rbsp_data) {
            free_nalu_elem(nalu_elem);
            return OV_ENOMEM;
        }

        if (dmx->epb_info.nb_epb) {
            uint32_t *epb_pos = NULL;
            epb_pos = ov_malloc(dmx->epb_info.nb_epb * sizeof(*epb_pos));
            if (!epb_pos) {
                free_nalu_elem(nalu_elem);
                ov_free(rbsp_data);
                return OV_ENOMEM;
            }

            memcpy(epb_pos, dmx->epb_info.epb_pos, dmx->epb_info.nb_epb * sizeof(*epb_pos));

            nalu_pending->nalu.epb_pos = epb_pos;
            nalu_pending->nalu.nb_epb = dmx->epb_info.nb_epb;
        }

        dmx->epb_info.nb_epb = 0;

        nalu_pending->nalu.rbsp_data = rbsp_data;
        nalu_pending->nalu.rbsp_size = dmx->rbsp_ctx.rbsp_size;

        memcpy(rbsp_data, dmx->rbsp_ctx.start, dmx->rbsp_ctx.rbsp_size);

        empty_rbsp_cache(&dmx->rbsp_ctx);

        append_nalu_elem(nalu_list, nalu_pending);
    } else {
        ov_log(dmx, OVLOG_TRACE, "No pending nalu when processing start_code, skipping.\n");
        empty_rbsp_cache(&dmx->rbsp_ctx);
    }

    nalu_elem->nalu.type = nalu_type;

    dmx->nalu_pending = nalu_elem;

    return 0;
}

static int
process_emulation_prevention_byte(OVDemux *const dmx, struct ReaderCache *const cache_ctx,
                                  const struct RBSPSegment *sgmt_ctx)
{
    struct EPBCacheInfo *const epb_info = &dmx->epb_info;

    append_rbsp_segment_to_cache(cache_ctx, &dmx->rbsp_ctx, sgmt_ctx);

    if (epb_info->nb_epb + 1 > (epb_info->cache_size)/sizeof(*epb_info->epb_pos)) {
        int ret = extend_epb_cache(epb_info);
        if (ret < 0) {
            ov_log(dmx, OVLOG_ERROR, "ERROR extending cache\n");
            return ret;
        }
    }

    /* FIXME new computation of epb position */
    epb_info->epb_pos[epb_info->nb_epb] = dmx->rbsp_ctx.rbsp_size - 1 /*+ epb_info->nb_epb*/;
    epb_info->nb_epb++;

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
extract_cache_segments(OVDemux *const dmx, struct ReaderCache *const cache_ctx)
{
    const uint8_t *byte = cache_ctx->cache_start;
    const uint8_t *const cache_end = cache_ctx->cache_end;
    uint32_t byte_pos = cache_ctx->first_pos;
    uint8_t end_of_cache;
    struct RBSPSegment sgmt_ctx = {0};

    sgmt_ctx.start_p = byte + byte_pos;
    sgmt_ctx.end_p   = byte + byte_pos;

    do {
        const uint8_t *bytestream = &byte[byte_pos];

        /* FIXME we will actually loop over this more than once even if a start
         * code has been detected. This is a bit inefficient
         * TODO bin tricks for two fast zero bytes detection
         */
        if (*bytestream == 0) {
            int ret;
            ret = ovannexb_check_stc_or_epb(bytestream);
            if (ret < 0) {
                ov_log(dmx, OVLOG_ERROR, "Invalid raw VVC data\n");
                ret = OV_INVALID_DATA;
            }

            if (ret) {
                enum RBSPSegmentDelimiter dlm = ret;

                switch (dlm) {
                case ANNEXB_STC:
                    sgmt_ctx.end_p = bytestream;

                    ret = process_start_code(dmx, cache_ctx, &sgmt_ctx);

                    /* Next segment start is located after start code three bytes */
                    sgmt_ctx.end_p = sgmt_ctx.start_p = bytestream + 3;
                    if (sgmt_ctx.end_p > cache_end) {
                        ov_log(dmx, OVLOG_DEBUG, "STC over cache end\n");
                    }
                    break;
                case ANNEXB_EPB:
                    /* Keep the two zero bytes of emulation prevention three bytes */
                    sgmt_ctx.end_p = bytestream + 2;

                    ret = process_emulation_prevention_byte(dmx, cache_ctx, &sgmt_ctx);

                    /* Remove the emulation prevention 0x03 byte */
                    sgmt_ctx.end_p = sgmt_ctx.start_p = bytestream + 3;
                    if (sgmt_ctx.end_p > cache_end) {
                        ov_log(dmx, OVLOG_DEBUG, "EBP over cache end\n");
                    }

                    break;
                default:
                    /* FIXME we should not have something different from STC or
                     * EPB here
                     */
                    ov_log(dmx, OVLOG_ERROR, "Invalid raw VVC data\n");
                    ret = OV_INVALID_DATA;
                    break;
                }

                if (ret < 0) {
                    return ret;
                }
                byte_pos += 2;
            }
        }
        /* Note we actually mean >= here since the last bytes reside
         * into padded area sometimes we might read up to 6 bytes ahead
         * of cache end and the next segment start might be located
         * at cache end + 2 in case an Emulation prevention three bytes
         * overlap cache end
         */
        end_of_cache = &byte[++byte_pos] >= cache_end;

    } while (!end_of_cache);

    if (dmx->eof) {
        ov_log(dmx, OVLOG_TRACE, "EOF reached\n");
        sgmt_ctx.end_p = cache_end;
        return process_start_code(dmx, cache_ctx, &sgmt_ctx);
    }

    /* Keep track of overlapping removed start code or EBP*/
    cache_ctx->first_pos = &byte[byte_pos] - cache_end;

    /* Recopy cache to RBSP cache before refill */
    if (sgmt_ctx.start_p < byte + byte_pos) {
        sgmt_ctx.end_p = byte + byte_pos;
        append_rbsp_segment_to_cache(cache_ctx, &dmx->rbsp_ctx, &sgmt_ctx);
    }

    return 0;
}

static int
init_rbsp_cache(struct RBSPCacheData *const rbsp_ctx)
{
    rbsp_ctx->start = ov_mallocz(OVRBSP_CACHE_SIZE);
    if (rbsp_ctx->start == NULL) {
        return OV_ENOMEM;
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
        return OV_ENOMEM;
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
        return OV_ENOMEM;
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
        return OV_ENOMEM;
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
