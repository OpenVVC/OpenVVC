#ifndef RCN_TRANSFORM_TREE
#define RCN_TRANSFORM_TREE

#include <stdint.h>
#include <string.h>

#include "ovdefs.h"

struct ISPTUInfo;
struct TUInfo;

void rcn_tu_st(OVCTUDec *const ctu_dec,
               uint8_t x0, uint8_t y0,
               uint8_t log2_tb_w, uint8_t log2_tb_h,
               uint8_t cu_flags, uint8_t cbf_mask,
               const struct TUInfo *const tu_info);

void rcn_tu_c(OVCTUDec *const ctu_dec, uint8_t x0, uint8_t y0,
              uint8_t log2_tb_w, uint8_t log2_tb_h,
              uint8_t cu_flags, uint8_t cbf_mask,
              const struct TUInfo *const tu_info);

void rcn_transform_tree(OVCTUDec *const ctu_dec, uint8_t x0, uint8_t y0,
                        uint8_t log2_tb_w, uint8_t log2_tb_h, uint8_t log2_max_tb_s,
                        uint8_t cu_flags, const struct TUInfo *const tu_info);

void recon_isp_subtree_v(OVCTUDec *const ctudec,
                         unsigned int x0, unsigned int y0,
                         unsigned int log2_cb_w, unsigned int log2_cb_h,
                         uint8_t intra_mode,
                         const struct ISPTUInfo *const tu_info);

void recon_isp_subtree_h(OVCTUDec *const ctudec,
                         unsigned int x0, unsigned int y0,
                         unsigned int log2_cb_w, unsigned int log2_cb_h,
                         uint8_t intra_mode,
                         const struct ISPTUInfo *const tu_info);

#endif

