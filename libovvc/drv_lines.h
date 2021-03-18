#ifndef DRV_LINES_H
#define DRV_LINES_H
#include "ovdefs.h"

int init_drv_lines(OVSliceDec *sldec, const OVPS *const prms);

void reset_drv_lines(OVSliceDec *sldec, const OVPS *const prms);

void drv_line_next_line(OVCTUDec *const ctudec, const struct DRVLines *const lns);

#if 0
void drv_line_next_ctu(OVCTUDec *const ctudec, OVSliceDec *sldec, struct DRVLines *drv_line,
                       const OVPS *const prms, uint16_t ctb_x);
#endif

void drv_lines_uninit(OVSliceDec *sldec);

void store_inter_maps(const struct DRVLines *const l,
                      OVCTUDec *const ctudec,
                      unsigned int ctb_x, uint8_t is_last);

void dbf_load_info(struct DBFInfo *const dbf_info,
                   const struct DBFLines *const dbf_lines,
                   uint8_t log2_ctu_s, int ctb_x);

void dbf_store_info(struct DBFInfo *const dbf_info,
                    const struct DBFLines *const dbf_lines,
                    uint8_t log2_ctu_s, int ctb_x);
void offset_drv_lines(struct DRVLines *const lns, uint8_t tile_x, uint8_t tile_y,
                      uint8_t ctb_x,
                      uint8_t log2_ctb_s, uint8_t log2_min_cb_s,
                      uint8_t  nb_tile_cols, uint16_t nb_ctb_pic_w);
#endif
