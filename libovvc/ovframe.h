#ifndef OV_fRAME_H
#define OV_fRAME_H

#include <stddef.h>
#include <stdint.h>
#include "ovdefs.h"
#include <stdatomic.h>

/* TODO Decide on values */
enum ChromaFmt
{
    OV_YUV_420_P8  = 0,
    OV_YUV_420_P10 = 1,
};

struct MemPool;
/* Miscelaneous information on Picture */
struct FrameInfo
{
    enum ChromaFmt chromat_format;

    /* TODO add VUI info
     *      add RefPicList info
     */
};

struct PlaneProp
{
    uint16_t stride;
    uint16_t width;
    uint16_t height;
    uint16_t depth;
};

struct FramePool
{
    struct MemPool *plane_pool[4];
    struct PlaneProp plane_prop[4];
};

struct FrameInternal
{
    /* reference counter */
    atomic_uint ref_count;

    struct frame_pool *frame_pool;
    void *pool_elem[4];
};

struct Frame
{
    /* Pointer to Picture data planes per component */
    uint8_t *data[3];

    /* Per component line size in bytes */
    size_t linesize[3];
    size_t width[3];
    size_t height[3];

    /* Picture Order Count */
    uint32_t poc;

    struct FrameInfo frame_info;

    /* private data */
    struct FrameInternal internal;
    /* TODO we could attach associated NAL Units to
     * the frame so the user might process things as SEI
     * outside of the actual decoding process
     */
};


int ovframe_new_ref(OVFrame **dst, OVFrame *src);

void ovframe_unref(OVFrame **frame);

int framepool_request_planes(OVFrame *const frame, struct FramePool *const fp);
#endif
