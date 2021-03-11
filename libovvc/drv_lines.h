#ifndef DRV_LINES_H
#define DRV_LINES_H
#include "ovdefs.h"

int init_drv_lines(OVSliceDec *sldec, const OVPS *const prms);

void reset_drv_lines(OVSliceDec *sldec, const OVPS *const prms);

void drv_line_next_line(OVCTUDec *const ctudec, const OVSliceDec *const sldec);

void drv_line_next_ctu(OVCTUDec *const ctudec, OVSliceDec *sldec, struct DRVLines *drv_line,
                       const OVPS *const prms, uint16_t ctb_x);

void drv_lines_uninit(OVSliceDec *sldec);

void dbf_load_info(struct DBFInfo *const dbf_info,
                   const struct DBFLines *const dbf_lines,
                   uint8_t log2_ctu_s, int ctb_x);

void dbf_store_info(struct DBFInfo *const dbf_info,
                    const struct DBFLines *const dbf_lines,
                    uint8_t log2_ctu_s, int ctb_x);
#endif
