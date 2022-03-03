/**
 *
 *   OpenVVC is open-source real time software decoder compliant with the 
 *   ITU-T H.266- MPEG-I - Part 3 VVC standard. OpenVVC is developed from 
 *   scratch in C as a library that provides consumers with real time and
 *   energy-aware decoding capabilities under different OS including MAC OS,
 *   Windows, Linux and Android targeting low energy real-time decoding of
 *   4K VVC videos on Intel x86 and ARM platforms.
 * 
 *   Copyright (C) 2020-2022  IETR-INSA Rennes :
 *   
 *   Pierre-Loup CABARAT
 *   Wassim HAMIDOUCHE
 *   Guillaume GAUTIER
 *   Thomas AMESTOY
 *
 *   This library is free software; you can redistribute it and/or
 *   modify it under the terms of the GNU Lesser General Public
 *   License as published by the Free Software Foundation; either
 *   version 2.1 of the License, or (at your option) any later version.
 *
 *   This library is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *   Lesser General Public License for more details.
 *
 *   You should have received a copy of the GNU Lesser General Public
 *   License along with this library; if not, write to the Free Software
 *   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301
 *   USA
 * 
 **/

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

void store_ibc_maps(const struct DRVLines *const l,
                    OVCTUDec *const ctudec,
                    unsigned int ctb_x, uint8_t is_last);
#endif
