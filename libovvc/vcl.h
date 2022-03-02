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

#ifndef VCL_H
#define VCL_H
#include <stdint.h>
#include <string.h>

#include "cu_utils.h"
#include "dec_structures.h"
#include "ovdefs.h"
#include "vcl_cabac.h"

struct TUInfo;

int
transform_unit_wrap(OVCTUDec *const ctu_dec,
                    const OVPartInfo *const part_ctx,
                    uint8_t x0, uint8_t y0,
                    uint8_t log2_cb_w, uint8_t log2_cb_h,
                    VVCCU cu);

int transform_unit_st(OVCTUDec *const ctu_dec,
                      unsigned int x0, unsigned int y0,
                      unsigned int log2_tb_w, unsigned int log2_tb_h,
                      uint8_t rqt_root_cbf, CUFlags cu_flags, uint8_t tr_depth,
                      struct TUInfo *const tu_info);

int transform_unit_l(OVCTUDec *const ctu_dec,
                     unsigned int x0, unsigned int y0,
                     unsigned int log2_tb_w, unsigned int log2_tb_h,
                     uint8_t rqt_root_cbf, CUFlags cu_flags, uint8_t tr_depth,
                     struct TUInfo *const tu_info);

int transform_unit_c(OVCTUDec *const ctu_dec,
                     unsigned int x0, unsigned int y0,
                     unsigned int log2_tb_w, unsigned int log2_tb_h,
                     uint8_t rqt_root_cbf, CUFlags cu_flags, uint8_t tr_depth,
                     struct TUInfo *const tu_info);

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

VVCCU
coding_unit_inter_st(OVCTUDec *const ctu_dec,
                     const OVPartInfo *const part_ctx,
                     uint8_t x0, uint8_t y0,
                     uint8_t log2_cu_w, uint8_t log2_cu_h);

int
prediction_unit_inter_b(OVCTUDec *const ctu_dec,
                        const OVPartInfo *const part_ctx,
                        uint8_t x0, uint8_t y0,
                        uint8_t log2_pb_w, uint8_t log2_pb_h,
                        uint8_t skip_flag, uint8_t merge_flag);
int
prediction_unit_inter_p(OVCTUDec *const ctu_dec,
                        const OVPartInfo *const part_ctx,
                        uint8_t x0, uint8_t y0,
                        uint8_t log2_pb_w, uint8_t log2_pb_h,
                        uint8_t skip_flag, uint8_t merge_flag);

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

uint64_t residual_coding_isp_h_sdh(OVCTUDec *const ctu_dec, int16_t *const dst,
                              uint8_t log2_tb_w, uint8_t log2_tb_h,
                              uint16_t last_pos);

uint64_t residual_coding_isp_v_sdh(OVCTUDec *const ctu_dec, int16_t *const dst,
                              uint8_t log2_tb_w, uint8_t log2_tb_h,
                              uint16_t last_pos);

uint64_t residual_coding_isp_h_dpq(OVCTUDec *const ctu_dec, int16_t *const dst,
                              uint8_t log2_tb_w, uint8_t log2_tb_h,
                              uint16_t last_pos);

uint64_t residual_coding_isp_v_dpq(OVCTUDec *const ctu_dec, int16_t *const dst,
                              uint8_t log2_tb_w, uint8_t log2_tb_h,
                              uint16_t last_pos);

uint64_t residual_coding_sdh(OVCTUDec *const ctu_dec, int16_t *const dst,
                             uint8_t log2_tb_w, uint8_t log2_tb_h,
                             uint16_t last_pos);

uint64_t residual_coding_chroma_sdh(OVCTUDec *const ctu_dec, int16_t *const dst,
                                    uint8_t log2_tb_w, uint8_t log2_tb_h,
                                    uint16_t last_pos);

uint64_t residual_coding_dpq(OVCTUDec *const ctu_dec, int16_t *const dst,
                             uint8_t log2_tb_w, uint8_t log2_tb_h,
                             uint16_t last_pos);

uint64_t residual_coding_chroma_dpq(OVCTUDec *const ctu_dec, int16_t *const dst,
                                    uint8_t log2_tb_w, uint8_t log2_tb_h,
                                    uint16_t last_pos);

uint64_t residual_coding_ts(OVCTUDec *const ctu_dec, int16_t *dst, uint8_t log2_tb_w, uint8_t log2_tb_h, uint8_t bdpcm);

int coding_unit(OVCTUDec *const ctu_dec,
                const OVPartInfo *const part_ctx,
                uint8_t x0, uint8_t y0,
                uint8_t log2_cb_w, uint8_t log2_cb_h);


void ovcabac_read_ae_sao_ctu( OVCTUDec *const ctudec, int ctb_rs, uint16_t nb_ctu_w);

void ovcabac_read_ae_alf_ctu( OVCTUDec *const ctudec, uint16_t ctb_rs, uint16_t nb_ctu_w);

void ovcabac_read_ae_cc_alf_ctu(OVCTUDec *const ctudec, uint16_t ctb_rs, uint16_t nb_ctu_w);

#endif
