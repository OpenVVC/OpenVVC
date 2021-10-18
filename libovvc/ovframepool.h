#ifndef OVFRAMEPOOL_H
#define OVFRAMEPOOL_H
#include "ovframe.h"

struct FramePool;

void ovframepool_uninit(struct FramePool **fpool_p);

int ovframepool_init(struct FramePool **fpool_p, uint8_t fmt, uint16_t pic_w, uint16_t pic_h);

int ovframepool_request_planes(OVFrame *const frame, struct FramePool *const fp);

void ovframepool_release_planes(OVFrame *const frame);

#endif

