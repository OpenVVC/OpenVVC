#ifndef RCN_MC_H
#define RCN_MC_H

struct RCNFunctions;

void rcn_init_mc_functions(struct RCNFunctions *const rcn_funcs);

void put_weighted_ciip_pixels(uint16_t* dst, ptrdiff_t dststride,
                      const uint16_t* src_intra, const uint16_t* src_inter, ptrdiff_t srcstride,
                      int width, int height, int wt);

#endif // RCN_MC_H
