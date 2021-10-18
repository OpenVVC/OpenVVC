#ifndef OVDPB_INTERNAL_H
#define OVDPB_INTERNAL_H

#include "ovdpb.h"
#include "ovdefs.h"
#include "ovframe.h"

struct FramePool;

struct DPBInternal
{
    struct FramePool *frame_pool;
};

void dpb_uninit_framepool(struct DPBInternal *dpb_priv);

int dpb_init_framepool(struct DPBInternal *dpb_priv, const OVSPS *const sps);

int dpbpriv_request_frame(struct DPBInternal *dpb_priv, OVFrame **frame);

void dpbpriv_uninit_framepool(struct DPBInternal *dpb_priv);

int dpbpriv_init_framepool(struct DPBInternal *dpb_priv, const OVSPS *const sps);

#endif
