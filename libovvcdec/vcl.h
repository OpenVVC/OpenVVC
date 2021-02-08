#ifndef VCL_H
#define VCL_H
#include <stdint.h>

#include "ovdefs.h"

int transform_unit_st(OVCTUDec *const ctu_dec,
                      unsigned int x0, unsigned int y0,
                      unsigned int log2_tb_w, unsigned int log2_tb_h,
                      uint8_t rqt_root_cbf, uint8_t cu_flags);

int transform_unit_l(OVCTUDec *const ctu_dec,
                     unsigned int x0, unsigned int y0,
                     unsigned int log2_tb_w, unsigned int log2_tb_h,
                     uint8_t rqt_root_cbf, uint8_t cu_flags);

int transform_unit_c(OVCTUDec *const ctu_dec,
                     unsigned int x0, unsigned int y0,
                     unsigned int log2_tb_w, unsigned int log2_tb_h,
                     uint8_t rqt_root_cbf, uint8_t cu_flags);

VVCCU coding_unit_intra_st(OVCTUDec *const ctu_dec,
                           const OVPartInfo *const part_ctx,
                           uint8_t x0, uint8_t y0,
                           uint8_t log2_cu_w, uint8_t log2_cu_h);

VVCCU coding_unit_intra(OVCTUDec *const ctu_dec,
                        const OVPartInfo *const part_ctx,
                        uint8_t x0, uint8_t y0,
                        uint8_t log2_cb_w, uint8_t log2_cb_h);

VVCCU coding_unit_intra_c(OVCTUDec *const ctu_dec,
                          const OVPartInfo *const part_ctx,
                          uint8_t x0, uint8_t y0,
                          uint8_t log2_cb_w, uint8_t log2_cb_h);

int coding_quadtree(OVCTUDec *const ctu_dec,
                    const OVPartInfo *const part_ctx,
                    unsigned int x0, unsigned int y0,
                    unsigned int log2_cb_s, unsigned int qt_depth);

int coding_quadtree_implicit(OVCTUDec *const ctu_dec,
                             const OVPartInfo *const part_ctx,
                             unsigned int x0, unsigned int y0,
                             unsigned int log2_cb_s, unsigned int qt_depth,
                             unsigned int rem_w, unsigned int rem_h);

int dual_tree(OVCTUDec *const ctu_dec,
              const OVPartInfo *const part_ctx,
              unsigned int x0, unsigned int y0,
              unsigned int log2_cb_s, unsigned int qt_depth);

int dual_tree_implicit(OVCTUDec *const ctu_dec,
                       const OVPartInfo *const part_ctx,
                       unsigned int x0, unsigned int y0,
                       unsigned int log2_cb_s, unsigned int qt_depth,
                       unsigned int rem_w,
                       unsigned int rem_h);

#endif
