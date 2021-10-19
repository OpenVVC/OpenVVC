#ifndef OVFRAMEPOOL_H
#define OVFRAMEPOOL_H
#include "ovframe.h"

struct FramePool;

void ovframepool_uninit(struct FramePool **fpool_p);

int ovframepool_init(struct FramePool **fpool_p, uint8_t fmt, uint16_t pic_w, uint16_t pic_h);

OVFrame *ovframepool_request_frame(struct FramePool *fpool);

void ovframepool_release_frame(OVFrame **frame_p);
#endif

