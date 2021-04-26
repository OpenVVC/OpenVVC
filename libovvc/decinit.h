/* This file contains utilities to initialise the decoder
 */

#ifndef DECINIT_H
#define DECINIT_H
#include "ovdefs.h"

int decinit_update_params(struct OVPS *const ps, const OVNVCLCtx *const nvcl_ctx);

int decinit_set_entry_points(OVPS *const prms, const OVNALUnit *nal, uint32_t nb_sh_bytes);

#endif
