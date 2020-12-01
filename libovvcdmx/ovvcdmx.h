#ifndef OVVCDMX_H
#define OVVCDMX_H

#include "ovvcdmx.h"

typedef struct OVVCDmx OVVCDmx;

int ovdmx_init(OVVCDmx **vvcdmx);

void ovdmx_close(OVVCDmx **vvcdmx);
#endif

