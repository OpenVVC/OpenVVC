#ifndef RCN_SSE_H
#define RCN_SSE_H
#include "rcn_structures.h"

void rcn_init_mc_functions_neon(struct RCNFunctions *const rcn_funcs);
void rcn_init_tr_functions_neon(struct RCNFunctions *const rcn_funcs);
void rcn_init_dc_planar_functions_neon(struct RCNFunctions *const rcn_funcs);
// void rcn_init_ict_functions_neon(struct RCNFunctions *rcn_func, uint8_t type);
// void rcn_init_alf_functions_neon(struct RCNFunctions *rcn_func);
// void rcn_init_cclm_functions_neon(struct RCNFunctions *rcn_func);
// void rcn_init_lfnst_functions_neon(struct RCNFunctions *rcn_func);
// void rcn_init_mip_functions_neon(struct RCNFunctions *const rcn_funcs);
void rcn_init_sao_functions_neon(struct RCNFunctions *const rcn_funcs);


#endif//RCN_SSE_H
