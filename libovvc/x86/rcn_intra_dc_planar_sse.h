#ifndef VVC_INTRA_PRED_SSE_H
#define VVC_INTRA_PRED_SSE_H

#include "stdint.h"
#include "stddef.h"
#include "rcn_structures.h"

void rcn_init_dc_planar_functions_sse(struct RCNFunctions *const rcn_funcs);

#endif
