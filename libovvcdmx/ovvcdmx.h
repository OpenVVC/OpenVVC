#ifndef OVVCDMX_H
#define OVVCDMX_H

#include <stdio.h>

typedef struct OVVCDmx OVVCDmx;

int ovdmx_init(OVVCDmx **vvcdmx);

int ovdmx_close(OVVCDmx *vvcdmx);

int ovdmx_attach_stream(OVVCDmx *const dmx, FILE *fstream);

void ovdmx_detach_stream(OVVCDmx *const dmx);

int ovdmx_read_stream(OVVCDmx *const dmx);

#endif

