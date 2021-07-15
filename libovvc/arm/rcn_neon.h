#ifndef RCN_SSE_H
#define RCN_SSE_H
#include "rcn_structures.h"

void rcn_init_dc_planar_functions_neon(struct RCNFunctions *const rcn_funcs);
void rcn_init_sao_functions_neon(struct RCNFunctions *const rcn_funcs);


#endif//RCN_SSE_H
