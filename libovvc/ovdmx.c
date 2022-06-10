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

static const char *const demux_name = "Open VVC Annex B demuxer";

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
    OVNALUnit *nalu;
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
    const uint8_t *start;
    const uint8_t *end;

    /* Number of bytes to be skipped at start */
    uint32_t nb_skip;
};

struct NALUListStatus
{
    struct NALUnitsList nalu_list;
    int nb_nalus;
    uint8_t got_vcl;
    uint8_t got_nvcl_suffix;
    uint8_t got_ph;
};

struct OVDemux
{
    const char *name;

    /* Points to a read only IO context */
    OVIOStream *io_str;

    /* Information on current Stream Of Data Bytes (SODB)
     */
    struct ReaderCache rdr_cache;

    /* Cache used for RBSP extraction */
    struct RBSPCacheData rbsp_cache;

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

static int extract_cache_segments(OVDemux *const dmx, struct ReaderCache *const rdr_cache);

static int init_rbsp_cache(struct RBSPCacheData *const rbsp_cache);

/* Realloc rbsp_cache adding an extra 64KB to previously allocated size
   and copy previous content */
static int extend_rbsp_cache(struct RBSPCacheData *const rbsp_cache);

static void free_rbsp_cache(struct RBSPCacheData *const rbsp_cache);

static int init_epb_cache(struct EPBCacheInfo *const epb_info);

/* Realloc EPB locations caching adding an extra 16 elements to previously
 * allocated size and copy previous content
 */
static int extend_epb_cache(struct EPBCacheInfo *const epb_info);

static void free_epb_cache(struct EPBCacheInfo *const epb_info);

static void clear_nalu_list(struct NALUnitsList *list);

static int refill_reader_cache(struct ReaderCache *const rdr_cache,
                               OVIOStream *const io_str);

static void append_nalu_elem(struct NALUnitsList *const list, struct NALUnitListElem *elem);
static void prepend_nalu_elem(struct NALUnitsList *const list, struct NALUnitListElem *elem);

static void release_nalu_elem(struct NALUnitListElem *nalu_elem);

int
ovdmx_init(OVDemux **dmx_p)
{
    int ret = 0;
    OVDemux *dmx;

    dmx = ov_mallocz(sizeof(*dmx));

    if (!dmx) return OVVC_ENOMEM;

    *dmx_p = dmx;

    dmx->name = demux_name;
    dmx->io_str = NULL;

    dmx->nalu_elem_pool = ovmempool_init(sizeof(struct NALUnitListElem));

    if (dmx->nalu_elem_pool == NULL) {
        goto fail_pool_init;
    }

    ret = init_rbsp_cache(&dmx->rbsp_cache);
    if (ret < 0) {
        goto fail_rbsp_cache;
    }

    ret = init_epb_cache(&dmx->epb_info);
    if (ret < 0) {
        goto fail_epb_cache;
    }

    return 0;

fail_epb_cache:
    free_rbsp_cache(&dmx->rbsp_cache);

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

        free_rbsp_cache(&dmx->rbsp_cache);

        free_epb_cache(&dmx->epb_info);

        clear_nalu_list(&dmx->nalu_list);

        if (dmx->nalu_pending) {
            release_nalu_elem(dmx->nalu_pending);
        }

        ovmempool_uninit(&dmx->nalu_elem_pool);

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

    dmx->io_str = ovio_stream_open(io);

    if (dmx->io_str == NULL) {
        ov_log(dmx, OVLOG_ERROR, "Failed to open stream.\n");
        return OVVC_EINDATA;
    }

    dmx->rdr_cache.nb_skip = 0;

    return ret;
}

static void
empty_rbsp_cache(struct RBSPCacheData *rbsp_cache)
{
    rbsp_cache->end       = rbsp_cache->start;
    rbsp_cache->rbsp_size = 0;
}

static void
empty_epb_cache(struct EPBCacheInfo *const epb_info)
{
    epb_info->nb_epb = 0;
}

void
ovdmx_detach_stream(OVDemux *const dmx)
{
    if (dmx->io_str != NULL) {
        ovio_stream_close(dmx->io_str);
    }

    dmx->io_str = NULL;

    empty_rbsp_cache(&dmx->rbsp_cache);
    empty_epb_cache(&dmx->epb_info);
}

static int
refill_reader_cache(struct ReaderCache *const rdr_cache, OVIOStream *const io_str)
{
    int bytes_read = ovio_stream_read(&rdr_cache->start, io_str);

    rdr_cache->end   = rdr_cache->start + bytes_read;

    return bytes_read;
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
        elem->prev_nalu = NULL;
        elem->next_nalu = NULL;
    } else {
        /*FIXME ensure list is emptied */
        list->last_nalu = NULL;
    }

    return elem;
}

static int process_start_code(OVDemux *const dmx);

static int
extract_nal_unit(OVDemux *const dmx, struct NALUnitListElem **dst_nalup)
{
    struct NALUnitsList *nalu_list = &dmx->nalu_list;

    do {

        *dst_nalup = pop_nalu_elem(nalu_list);

        if (*dst_nalup) {

            return 1;

        } else if (!ovio_stream_eof(dmx->io_str)) {

            struct ReaderCache *const rdr_cache = &dmx->rdr_cache;

            int nb_bytes_read = refill_reader_cache(rdr_cache, dmx->io_str);

            int ret = extract_cache_segments(dmx, rdr_cache);

            if (ret < 0) {
                return ret;
            }


            if (ovio_stream_eof(dmx->io_str)) {
                ov_log(dmx, OVLOG_TRACE, "EOF reached\n");
                dmx->eof = 1;
                ret = process_start_code(dmx);
                if (ret < 0) {
                    return ret;
                }
            }

        } else {

            ov_log(dmx, OVLOG_TRACE, "IO EOF reached and list empty\n");

            return 0;
        }

    } while (1);
}

static void
clear_nalu_list(struct NALUnitsList *list)
{
    struct NALUnitListElem *elem = list->first_nalu;
    while (elem) {
        struct NALUnitListElem *to_release = elem;
        elem = elem->next_nalu;
        release_nalu_elem(to_release);
    }
    list->first_nalu = NULL;
    list->last_nalu = NULL;
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

static int
ovdmx_init_pu_from_list(OVPictureUnit **ovpu_p, struct NALUListStatus *const status)
{
    OVPictureUnit *ovpu;
    struct NALUnitsList *nalu_list = &status->nalu_list;
    struct NALUnitListElem *lelem = nalu_list->first_nalu;

    int ret = ovpu_init(ovpu_p, status->nb_nalus);

    if (ret < 0) {
        return ret;
    }

    ovpu = *ovpu_p;

    ov_log(NULL, OVLOG_ERROR, "Picture Unit containing %d NAL Units\n", status->nb_nalus);

    for (int i = 0; i < status->nb_nalus; i++) {
        ov_log(NULL, OVLOG_ERROR, "/t %s\n", nalu_name[lelem->nalu->type]);

        ov_nalu_new_ref(&ovpu->nalus[i], lelem->nalu);

        lelem = lelem->next_nalu;
    }

    *ovpu_p = ovpu;

    return 0;
}

static inline uint8_t
is_vcl(const struct OVNALUnit *const nalu) {
    return (nalu->type < OVNALU_OPI);
}

static inline uint8_t
is_nvcl_delimiter(const struct OVNALUnit *const nalu) {
    uint8_t nalu_type = nalu->type;

    uint8_t is_nvcl = (nalu->type >= OVNALU_OPI);

    uint8_t is_nvcl_inside  = nalu_type == OVNALU_SUFFIX_SEI;
    is_nvcl_inside |= nalu_type == OVNALU_SUFFIX_APS;
    is_nvcl_inside |= nalu_type == OVNALU_AUD;
    is_nvcl_inside |= nalu_type == OVNALU_EOS;
    is_nvcl_inside |= nalu_type == OVNALU_EOB;
    is_nvcl_inside |= nalu_type == OVNALU_FD;

    return  is_nvcl && !is_nvcl_inside;
}

static uint8_t
is_next_pu_start(struct NALUListStatus *const status,
                 const struct OVNALUnit *const nalu)
{
    status->got_ph |= nalu->type == OVNALU_PH;

    if (!status->got_vcl) {

        status->got_vcl |= is_vcl(nalu);
        goto no;

    } else if (!is_vcl(nalu) && !is_nvcl_delimiter(nalu)) {

        status->got_nvcl_suffix = 1;
        goto no;

    } else {
        if (is_vcl(nalu) && status->got_ph) {
            goto no;
        } else {
            goto yes;
        }
    }

no:
    return 0;
yes:
    return 1;

}

int
ovdmx_extract_picture_unit(OVDemux *const dmx, OVPictureUnit **dst_pu_p)
{
    struct NALUListStatus status = {0};
    int ret;

    do {
        struct NALUnitListElem *nalu;
        ret = extract_nal_unit(dmx, &nalu);
        if (ret < 0) {
            ov_log(dmx, OVLOG_ERROR, "Error extracting NAL Unit.\n");
            goto extraction_error;
        }

        if (!nalu) {
            ov_log(dmx, OVLOG_DEBUG, "No NALU available.\n");
            break;
        }

        if (!is_next_pu_start(&status, nalu->nalu)) {
            status.nb_nalus++;
            append_nalu_elem(&status.nalu_list, nalu);
        } else {
            /* Move NALU element back into demux list */
            prepend_nalu_elem(&dmx->nalu_list, nalu);
            break;
        }

    } while (1);

    ret = ovdmx_init_pu_from_list(dst_pu_p, &status);

extraction_error:
    /* We could also try to build a Picture Unit however this behaviour is safer*/
    clear_nalu_list(&status.nalu_list);

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

    ovnalu_init2(&nalu_elem->nalu);

    return nalu_elem;
}

static void
release_nalu_elem(struct NALUnitListElem *nalu_elem)
{
    if (nalu_elem->nalu) {
        ov_nalu_unref(&nalu_elem->nalu);
    }
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

static void
prepend_nalu_elem(struct NALUnitsList *const list, struct NALUnitListElem *elem)
{
    elem->prev_nalu = NULL;
    elem->next_nalu = list->first_nalu;
    if (list->first_nalu) {
        list->first_nalu->prev_nalu = elem;
    } else {
        list->last_nalu = elem;
    }
    list->first_nalu = elem;
}

static int
append_rbsp_segment_to_cache(struct RBSPCacheData *const rbsp_cache,
                             const struct RBSPSegment *const sgmt_ctx)
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
free_nalu_data(struct OVNALUnit **nalu_p)
{
    OVNALUnit *nalu = *nalu_p;
    ov_freep(&nalu->rbsp_data);
    ov_freep(nalu_p);
}

static int
allocate_nalu_data(struct OVNALUnit *const nalu,
                   const struct EPBCacheInfo *const epb_info,
                   const struct RBSPCacheData *const rbsp_cache)
{
    size_t nalu_data_size = rbsp_cache->rbsp_size + OV_RBSP_PADDING;
    size_t epb_data_size  = epb_info->nb_epb * sizeof(*nalu->epb_pos);
    size_t alloc_size = nalu_data_size + epb_data_size;

    uint8_t *data = ov_malloc(alloc_size);
    if (!data) {
        return OVVC_ENOMEM;
    }

    memcpy(data, rbsp_cache->start, rbsp_cache->rbsp_size);

    nalu->rbsp_data = data;
    nalu->rbsp_size = rbsp_cache->rbsp_size;

    /* FIXME temporary compat for NALU free */
    nalu->release = &free_nalu_data;

    /* Set padding area to zero */
    memset(data + rbsp_cache->rbsp_size, 0, OV_RBSP_PADDING);

    nalu->nb_epb  = epb_info->nb_epb;

    if (epb_data_size) {

        data += nalu_data_size;

        memcpy(data, epb_info->epb_pos, epb_data_size);

        nalu->epb_pos = (uint32_t*) data;

        return 0;
    }

    nalu->epb_pos = NULL;

    return 0;
}

static int
process_start_code(OVDemux *const dmx)
{
    struct NALUnitListElem *nalu_pending = dmx->nalu_pending;

    if (nalu_pending) {
        enum OVNALUType nalu_type = (dmx->rbsp_cache.start[1] >> 3) & 0x1F;

        struct OVNALUnit *const nalu = nalu_pending->nalu;

        int ret = allocate_nalu_data(nalu, &dmx->epb_info, &dmx->rbsp_cache);

        if (ret < 0) {
            empty_rbsp_cache(&dmx->rbsp_cache);
            empty_epb_cache(&dmx->epb_info);
            return OVVC_ENOMEM;
        }

        nalu->type = nalu_type;

        append_nalu_elem(&dmx->nalu_list, nalu_pending);
    } else {
        ov_log(dmx, OVLOG_TRACE, "No pending nalu when processing start_code, skipping.\n");
    }

    empty_rbsp_cache(&dmx->rbsp_cache);
    empty_epb_cache(&dmx->epb_info);

    struct NALUnitListElem *nalu_elem = create_nalu_elem(dmx);
    if (!nalu_elem) {
        ov_log(dmx, OVLOG_ERROR, "Could not alloc NALU element\n");
        return OVVC_ENOMEM;
    }

    dmx->nalu_pending = nalu_elem;

    return 0;
}

static int
process_emulation_prevention_byte(OVDemux *const dmx)
{
    struct EPBCacheInfo *const epb_info = &dmx->epb_info;

    if (epb_info->nb_epb + 1 > (epb_info->cache_size) / sizeof(*epb_info->epb_pos)) {
        int ret = extend_epb_cache(epb_info);
        if (ret < 0) {
            ov_log(dmx, OVLOG_ERROR, "ERROR extending cache\n");
            return ret;
        }
    }

    epb_info->epb_pos[epb_info->nb_epb] = dmx->rbsp_cache.rbsp_size - 1;
    epb_info->nb_epb++;

    return 0;
}

static int
process_rbsp_delimiter(OVDemux *const dmx, enum RBSPSegmentDelimiter dlm)
{
    switch (dlm) {
        case ANNEXB_STC:

            return process_start_code(dmx);

            break;
        case ANNEXB_EPB:

            return process_emulation_prevention_byte(dmx);

            break;
        default:
            ov_log(dmx, OVLOG_ERROR, "Invalid raw VVC data\n");
            return OVVC_EINDATA;
    }
}


/**
 * returns: OVVC_EINDATA on invalid data found
 *          OVVC_ENOMEM on allocation error from either RBSP/EPB caches or data buffer
 *          for NALU,
 *          0 otherwise;
 */
/* WARNING We need to be careful on endianness here if we plan
   to use bigger read sizes */
static int
extract_cache_segments(OVDemux *const dmx, struct ReaderCache *const rdr_cache)
{
    const uint8_t *cursor = rdr_cache->start + rdr_cache->nb_skip;
    struct RBSPSegment sgmt_ctx = {.start_p = cursor, .end_p = cursor};
    int ret;

    do {

        /* TODO bin tricks for two fast zero bytes detection */
        if (*cursor == 0) {

            enum RBSPSegmentDelimiter dlm = ovannexb_check_stc_or_epb(cursor);

            if (dlm) {

                /* Keep the delimiter first two bytes inside the segment */
                sgmt_ctx.end_p = cursor += 2;

                ret = append_rbsp_segment_to_cache(&dmx->rbsp_cache, &sgmt_ctx);
                if (ret < 0) goto error;

                ret = process_rbsp_delimiter(dmx, dlm);
                if (ret < 0) goto error;

                /* Next segment start is located after delimiter 3 bytes */
                sgmt_ctx.start_p = sgmt_ctx.end_p += 1;
            }
        }

    } while (++cursor < rdr_cache->end);

    /* Keep track of overlapping start code or EBP */
    rdr_cache->nb_skip = cursor - rdr_cache->end;

    if (sgmt_ctx.start_p < cursor) {

        sgmt_ctx.end_p = cursor;

        /* Recopy cache to RBSP cache before refill */
        ret = append_rbsp_segment_to_cache(&dmx->rbsp_cache, &sgmt_ctx);
        if (ret < 0) goto error;

    }

    return 0;

error:
    empty_rbsp_cache(&dmx->rbsp_cache);
    empty_epb_cache(&dmx->epb_info);
    return ret;
}

static int
init_rbsp_cache(struct RBSPCacheData *const rbsp_cache)
{
    rbsp_cache->start = ov_malloc(OVRBSP_CACHE_SIZE);
    if (!rbsp_cache->start) {
        return OVVC_ENOMEM;
    }

    rbsp_cache->end = rbsp_cache->start;
    rbsp_cache->cache_size = OVRBSP_CACHE_SIZE;

    return 0;
}

static void
free_rbsp_cache(struct RBSPCacheData *const rbsp_cache)
{
    ov_freep(&rbsp_cache->start);
}

static int
extend_rbsp_cache(struct RBSPCacheData *const rbsp_cache)
{
    uint8_t *old_cache = rbsp_cache->start;
    uint8_t *new_cache;
    size_t new_size = rbsp_cache->cache_size + OVRBSP_CACHE_SIZE;

    new_cache = ov_malloc(new_size);
    if (!new_cache) {
        return OVVC_ENOMEM;
    }

    memcpy(new_cache, old_cache, rbsp_cache->rbsp_size);

    ov_free(old_cache);

    rbsp_cache->start = new_cache;
    rbsp_cache->end = rbsp_cache->start + rbsp_cache->rbsp_size;
    rbsp_cache->cache_size = new_size;
    return 0;
}

static int
init_epb_cache(struct EPBCacheInfo *const epb_info)
{
    epb_info->epb_pos = ov_malloc(OVEPB_CACHE_SIZE);
    if (epb_info->epb_pos == NULL) {
        return OVVC_ENOMEM;
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
        return OVVC_ENOMEM;
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
