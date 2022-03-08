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

#include <stdlib.h>
#include <stdint.h>
#include "ovutils.h"
#include "rcn.h"
#include "rcn_fill_ref.h"
#include "rcn_intra_angular.h"
#include "data_rcn_angular.h"
#include "ctudec.h"
#include "drv.h"
#include "dbf_utils.h"


/* Modify angular mode for non square according to w / h
 * ratio.
 * WARNING: do not call if DC or PLANAR, or LM
 *          return value is not the mode to be used
 *          for derivation but for reconstruction.
 * FIXME clean return unsigned and smaller sizes
 * FIXME remove the + 2 if specialized angular modes
 */
static int
derive_wide_angular_mode(int log2_pb_w, int log2_pb_h, int pred_mode)
{
    static const uint8_t mode_shift_tab[6] = {0, 6, 10, 12, 14, 15};
    int mode_shift = mode_shift_tab[OVABS(log2_pb_w - log2_pb_h)];

    if (log2_pb_w > log2_pb_h && pred_mode < 2 + mode_shift) {
        pred_mode += (OVINTRA_VDIA - 1);
    } else if (log2_pb_h > log2_pb_w && pred_mode > OVINTRA_VDIA - mode_shift) {
        pred_mode -= (OVINTRA_VDIA - 1);
    }
    return pred_mode;
}

static void
intra_angular_gauss_v(const struct OVRCNCtx *rcn_ctx, OVSample *ref1, OVSample *ref2, OVSample *dst, int dst_stride,
                      uint8_t log2_pb_w, uint8_t log2_pb_h, int8_t mode_idx)
{
    const struct IntraAngularFunctions *gauss_v = rcn_ctx->ctudec->rcn_funcs.intra_angular_gauss_v;
    OVSample filtered_ref_abv[(128 << 1) + 128];
    OVSample filtered_ref_lft[(128 << 1) + 128];
    switch (mode_idx) {
        case 0:
            gauss_v->pure_pdpc(ref1, ref2, dst, dst_stride, log2_pb_w, log2_pb_h);
            break;
        case (16):
            {
                int abv_ref_length = 1 << (log2_pb_w + 1);
                int lft_ref_length = 1 << (log2_pb_h + 1);

                rcn_ctx->ctudec->rcn_funcs.tmp.filter_ref_samples(ref1, filtered_ref_abv, ref2,
                                                          abv_ref_length);
                rcn_ctx->ctudec->rcn_funcs.tmp.filter_ref_samples(ref2, filtered_ref_lft, ref1,
                                                          lft_ref_length);

                ref1 = filtered_ref_abv;
                ref2 = filtered_ref_lft;

                gauss_v->diagonal_pdpc(ref1, ref2, dst, dst_stride, log2_pb_w, log2_pb_h);
            }
            break;
        default:
            {
                int abs_angle_val = angle_table[OVABS(mode_idx)];
                int inv_angle = inverse_angle_table[OVABS(mode_idx)];
                uint8_t req_frac = !!(abs_angle_val & 0x1F);
                int8_t pdpc_scale = OVMIN(2, log2_pb_h - (floor_log2(3 * inverse_angle_table[OVABS(mode_idx)] - 2) - 8));

                if (!req_frac) {
                    const struct IntraAngularFunctions *nofrac_v = rcn_ctx->ctudec->rcn_funcs.intra_angular_nofrac_v;
                    if (mode_idx < 0) {
                        int abv_ref_length = 1 << (log2_pb_w + 1);
                        int lft_ref_length = 1 << (log2_pb_h + 1);
                        int pu_w = 1 << log2_pb_w;
                        int pu_h = 1 << log2_pb_h;
                        int inv_angle_sum = 256;

                        rcn_ctx->ctudec->rcn_funcs.tmp.filter_ref_samples(ref1, filtered_ref_abv + pu_h,
                                                                  ref2, abv_ref_length);
                        rcn_ctx->ctudec->rcn_funcs.tmp.filter_ref_samples(ref2, filtered_ref_lft + pu_w,
                                                                  ref1, lft_ref_length);
                        ref1 = filtered_ref_abv + pu_h;
                        ref2 = filtered_ref_lft + pu_w;

                        for (int k = -1; k >= -pu_h; k--) {
                            inv_angle_sum += inv_angle;
                            ref1[k] = ref2[OVMIN(inv_angle_sum >> 9,pu_h)];
                        }

                        nofrac_v->angular(ref1, dst, dst_stride, log2_pb_w, log2_pb_h,
                                          -abs_angle_val);

                    } else if (pdpc_scale < 0) {
                        int abv_ref_length = 1 << (log2_pb_w + 1);
                        rcn_ctx->ctudec->rcn_funcs.tmp.filter_ref_samples(ref1, filtered_ref_abv, ref2,
                                                                  abv_ref_length);
                        ref1 = filtered_ref_abv;

                        nofrac_v->angular(ref1, dst, dst_stride, log2_pb_w, log2_pb_h,
                                          abs_angle_val);

                    } else {
                        int abv_ref_length = 1 << (log2_pb_w + 1);
                        int lft_ref_length = 1 << (log2_pb_h + 1);

                        rcn_ctx->ctudec->rcn_funcs.tmp.filter_ref_samples(ref1, filtered_ref_abv, ref2,
                                           abv_ref_length);
                        rcn_ctx->ctudec->rcn_funcs.tmp.filter_ref_samples(ref2, filtered_ref_lft, ref1,
                                           lft_ref_length);

                        ref1 = filtered_ref_abv;
                        ref2 = filtered_ref_lft;

                        nofrac_v->angular_pdpc(ref1, ref2, dst, dst_stride, log2_pb_w, log2_pb_h,
                                               mode_idx);
                    }
                } else {
                    if (mode_idx < 0) {
                        int pu_h = 1 << log2_pb_h;
                        int inv_angle_sum = 256;
                        for (int k = -1; k >= -pu_h; k-- ){
                            inv_angle_sum += inv_angle;
                            ref1[k] = ref2[OVMIN(inv_angle_sum >> 9,pu_h)];
                        }

                        gauss_v->angular(ref1, dst, dst_stride, log2_pb_w, log2_pb_h,
                                         -abs_angle_val);

                    } else if (pdpc_scale < 0) {

                        gauss_v->angular(ref1, dst, dst_stride, log2_pb_w, log2_pb_h,
                                         abs_angle_val);

                    } else {

                        gauss_v->angular_pdpc(ref1, ref2, dst, dst_stride,
                                              log2_pb_w, log2_pb_h,
                                              mode_idx);
                    }
                }
            }
            break;
    }
}

static void
intra_angular_gauss_h(const struct OVRCNCtx *rcn_ctx, OVSample *ref1, OVSample *ref2, OVSample *dst, int dst_stride,
                      uint8_t log2_pb_w, uint8_t log2_pb_h, int8_t mode_idx)
{
    const struct IntraAngularFunctions *gauss_h = rcn_ctx->ctudec->rcn_funcs.intra_angular_gauss_h;
    OVSample filtered_ref_abv[(128 << 1) + 128];
    OVSample filtered_ref_lft[(128 << 1) + 128];
    switch (mode_idx) {
        case 0:
            gauss_h->pure_pdpc(ref1, ref2, dst, dst_stride, log2_pb_w, log2_pb_h);
            break;

        case (16):
            {
                int abv_ref_length = 1 << (log2_pb_w + 1);
                int lft_ref_length = 1 << (log2_pb_h + 1);
                rcn_ctx->ctudec->rcn_funcs.tmp.filter_ref_samples(ref1, filtered_ref_abv, ref2,
                                                          abv_ref_length);
                rcn_ctx->ctudec->rcn_funcs.tmp.filter_ref_samples(ref2, filtered_ref_lft, ref1,
                                                          lft_ref_length);
                ref1 = filtered_ref_abv;
                ref2 = filtered_ref_lft;

                gauss_h->diagonal_pdpc(ref1, ref2, dst, dst_stride, log2_pb_w, log2_pb_h);
            }
            break;
        default:
            {
                int abs_angle_val = angle_table[OVABS(mode_idx)];
                int inv_angle = inverse_angle_table[OVABS(mode_idx)];
                uint8_t req_frac = !!(abs_angle_val & 0x1F);
                int8_t pdpc_scale = OVMIN(2, log2_pb_w - (floor_log2(3 * inverse_angle_table[OVABS(mode_idx)] - 2) - 8));

                if (!req_frac) {
                    const struct IntraAngularFunctions *nofrac_h = rcn_ctx->ctudec->rcn_funcs.intra_angular_nofrac_h;
                    if (mode_idx < 0){
                        int pu_w = 1 << log2_pb_w;
                        int pu_h = 1 << log2_pb_h;

                        int abv_ref_length = 1 << (log2_pb_w + 1);
                        int lft_ref_length = 1 << (log2_pb_h + 1);
                        int inv_angle_sum  = 256;

                        rcn_ctx->ctudec->rcn_funcs.tmp.filter_ref_samples(ref1, filtered_ref_abv + pu_h,
                                                                  ref2, abv_ref_length);
                        rcn_ctx->ctudec->rcn_funcs.tmp.filter_ref_samples(ref2, filtered_ref_lft + pu_w,
                                                                  ref1, lft_ref_length);
                        ref1 = filtered_ref_abv + pu_h;
                        ref2 = filtered_ref_lft + pu_w;

                        for (int k = -1; k >= -pu_w; k--) {
                            inv_angle_sum += inv_angle;
                            ref2[k] = ref1[OVMIN(inv_angle_sum >> 9,pu_w)];
                        }

                        nofrac_h->angular(ref2, dst, dst_stride, log2_pb_w, log2_pb_h,
                                          -abs_angle_val);

                    } else if (pdpc_scale < 0) {
                        int lft_ref_length = 1 << (log2_pb_h + 1);
                        rcn_ctx->ctudec->rcn_funcs.tmp.filter_ref_samples(ref2, filtered_ref_lft,
                                                                  ref1, lft_ref_length);
                        ref2 = filtered_ref_lft;
                        nofrac_h->angular(ref2, dst, dst_stride, log2_pb_w, log2_pb_h,
                                          abs_angle_val);
                    } else {
                        int abv_ref_length = 1 << (log2_pb_w + 1);
                        int lft_ref_length = 1 << (log2_pb_h + 1);

                        rcn_ctx->ctudec->rcn_funcs.tmp.filter_ref_samples(ref1, filtered_ref_abv,
                                                                  ref2, abv_ref_length);
                        rcn_ctx->ctudec->rcn_funcs.tmp.filter_ref_samples(ref2, filtered_ref_lft,
                                                                  ref1, lft_ref_length);

                        ref1 = filtered_ref_abv;
                        ref2 = filtered_ref_lft;

                        nofrac_h->angular_pdpc(ref1, ref2, dst, dst_stride,
                                               log2_pb_w, log2_pb_h,
                                               mode_idx);
                    }
                } else {
                    if (mode_idx < 0){
                        int pu_w = 1 << log2_pb_w;

                        int inv_angle_sum    = 256;
                        for (int k = -1; k >= -pu_w; k--) {
                            inv_angle_sum += inv_angle;
                            ref2[k] = ref1[OVMIN(inv_angle_sum >> 9,pu_w)];
                        }
                        gauss_h->angular(ref2, dst, dst_stride, log2_pb_w, log2_pb_h,
                                         -abs_angle_val);

                    } else if (pdpc_scale < 0) {

                        gauss_h->angular(ref2, dst, dst_stride, log2_pb_w, log2_pb_h,
                                         abs_angle_val);

                    } else {
                        gauss_h->angular_pdpc(ref1, ref2, dst, dst_stride,
                                              log2_pb_w, log2_pb_h,
                                              mode_idx);
                    }
                }
            }
            break;
    }
}

static void
intra_angular_cubic_v(const struct OVRCNCtx *rcn_ctx, OVSample *ref1, OVSample *ref2, OVSample *dst, int dst_stride,
                      uint8_t log2_pb_w, uint8_t log2_pb_h, int8_t mode_idx)
{
    const struct IntraAngularFunctions *cubic_v = rcn_ctx->ctudec->rcn_funcs.intra_angular_cubic_v;

    switch (mode_idx) {
        case 0:
            if (log2_pb_h > 1) {
                cubic_v->pure_pdpc(ref1, ref2, dst, dst_stride, log2_pb_w,
                                          log2_pb_h);
            } else {
                cubic_v->pure(ref1, dst, dst_stride, log2_pb_w, log2_pb_h);
            }
            break;
        case (16):
            if (log2_pb_h > 1) {
                cubic_v->diagonal_pdpc(ref1, ref2, dst, dst_stride, log2_pb_w, log2_pb_h);
            } else {
                cubic_v->diagonal(ref1, dst, dst_stride, log2_pb_w, log2_pb_h);
            }
            break;
        default:
            {
                const struct IntraAngularFunctions *nofrac_v = rcn_ctx->ctudec->rcn_funcs.intra_angular_nofrac_v;
                int abs_angle_val = angle_table[OVABS(mode_idx)];
                int inv_angle = inverse_angle_table[OVABS(mode_idx)];
                uint8_t req_frac = !!(abs_angle_val & 0x1F);
                int8_t pdpc_scale = OVMIN(2, log2_pb_h - (floor_log2(3 * inverse_angle_table[OVABS(mode_idx)] - 2) - 8));
                if (!req_frac) {
                    if (mode_idx < 0) {
                        int pu_h = 1 << log2_pb_h;
                        int inv_angle_sum = 256;

                        for (int k = -1; k >= -pu_h; k--){
                            inv_angle_sum += inv_angle;
                            ref1[k] = ref2[OVMIN(inv_angle_sum >> 9,pu_h)];
                        }

                        nofrac_v->angular(ref1, dst, dst_stride, log2_pb_w, log2_pb_h,
                                          -abs_angle_val);

                    } else if (pdpc_scale < 0 || log2_pb_h < 2){

                        nofrac_v->angular(ref1, dst, dst_stride, log2_pb_w, log2_pb_h,
                                          abs_angle_val);

                    } else {

                        nofrac_v->angular_pdpc(ref1, ref2, dst, dst_stride,
                                               log2_pb_w, log2_pb_h,
                                               mode_idx);
                    }
                } else {
                    if (mode_idx < 0) {
                        int pu_h = 1 << log2_pb_h;
                        int inv_angle_sum = 256;

                        for (int k = -1; k >= -pu_h; k--){
                            inv_angle_sum += inv_angle;
                            ref1[k] = ref2[OVMIN(inv_angle_sum >> 9,pu_h)];
                        }

                        cubic_v->angular(ref1, dst, dst_stride, log2_pb_w, log2_pb_h,
                                         -abs_angle_val);

                    } else if (pdpc_scale < 0 || log2_pb_h < 2){

                        cubic_v->angular(ref1, dst, dst_stride, log2_pb_w, log2_pb_h,
                                         abs_angle_val);

                    } else {

                        cubic_v->angular_pdpc(ref1, ref2, dst, dst_stride,
                                              log2_pb_w, log2_pb_h,
                                              mode_idx);
                    }
                }
            }
            break;
    }
}

static void
intra_angular_cubic_h(const struct OVRCNCtx *rcn_ctx, OVSample *ref1, OVSample *ref2, OVSample *dst, int dst_stride,
                      uint8_t log2_pb_w, uint8_t log2_pb_h, int8_t mode_idx)
{
    const struct IntraAngularFunctions *cubic_h = rcn_ctx->ctudec->rcn_funcs.intra_angular_cubic_h;
    switch (mode_idx) {
        case 0:
            if (log2_pb_h > 1){
                cubic_h->pure_pdpc(ref1, ref2, dst, dst_stride, log2_pb_w, log2_pb_h);
            } else {
                cubic_h->pure(ref2, dst, dst_stride, log2_pb_w, log2_pb_h);
            }
            break;

        case (16):
            if (log2_pb_h > 1){
                cubic_h->diagonal_pdpc(ref1, ref2, dst, dst_stride, log2_pb_w, log2_pb_h);
            } else {
                cubic_h->diagonal(ref2, dst, dst_stride, log2_pb_w, log2_pb_h);
            }
            break;
        default:
            {
                const struct IntraAngularFunctions *nofrac_h = rcn_ctx->ctudec->rcn_funcs.intra_angular_nofrac_h;
                int abs_angle_val = angle_table[OVABS(mode_idx)];
                int inv_angle = inverse_angle_table[OVABS(mode_idx)];
                uint8_t req_frac = !!(abs_angle_val & 0x1F);
                int8_t pdpc_scale = OVMIN(2, log2_pb_w - (floor_log2(3 * inverse_angle_table[OVABS(mode_idx)] - 2) - 8));
                if (!req_frac) {
                    if (mode_idx < 0){
                        int pu_w  = 1 << log2_pb_w;
                        int inv_angle_sum = 256;

                        for (int k = -1; k >= -pu_w; k--) {
                            inv_angle_sum += inv_angle;
                            ref2[k] = ref1[OVMIN(inv_angle_sum >> 9, pu_w)];
                        }

                        nofrac_h->angular(ref2, dst, dst_stride, log2_pb_w, log2_pb_h,
                                          -abs_angle_val);

                    } else if (pdpc_scale < 0 || log2_pb_h < 2) {

                        nofrac_h->angular(ref2, dst, dst_stride, log2_pb_w, log2_pb_h,
                                                 abs_angle_val);

                    } else {

                        nofrac_h->angular_pdpc(ref1, ref2, dst, dst_stride,
                                               log2_pb_w, log2_pb_h,
                                               mode_idx);
                    }

                } else {
                    if (mode_idx < 0){
                        int pu_w  = 1 << log2_pb_w;
                        int inv_angle_sum = 256;

                        for (int k = -1; k >= -pu_w; k--) {
                            inv_angle_sum += inv_angle;
                            ref2[k] = ref1[OVMIN(inv_angle_sum >> 9, pu_w)];
                        }

                        cubic_h->angular(ref2, dst, dst_stride, log2_pb_w, log2_pb_h,
                                         -abs_angle_val);

                    } else if (pdpc_scale < 0 || log2_pb_h < 2) {

                        cubic_h->angular(ref2, dst, dst_stride, log2_pb_w, log2_pb_h,
                                         abs_angle_val);
                    } else {

                        cubic_h->angular_pdpc(ref1, ref2, dst, dst_stride,
                                              log2_pb_w, log2_pb_h,
                                              mode_idx);
                    }
                }
            }
            break;
    }
}

static void
intra_angular_v(const struct OVRCNCtx *rcn_ctx, OVSample *ref1, OVSample *ref2, OVSample *dst, int dst_stride,
                uint8_t log2_pb_w, uint8_t log2_pb_h, int8_t pred_mode)
{
    int mode_idx = pred_mode - (int)OVINTRA_VER;
    uint8_t log2_nb_smp = log2_pb_w + log2_pb_h;
    uint8_t use_gauss_filter = log2_nb_smp > 5 && OVABS(mode_idx) > intra_filter[log2_nb_smp >> 1];

    if (use_gauss_filter) {
        intra_angular_gauss_v(rcn_ctx, ref1, ref2, dst, dst_stride, log2_pb_w, log2_pb_h, mode_idx);
    } else {
        intra_angular_cubic_v(rcn_ctx, ref1, ref2, dst, dst_stride, log2_pb_w, log2_pb_h, mode_idx);
    }
}

static void
intra_angular_h(const struct OVRCNCtx *rcn_ctx, OVSample *ref1, OVSample *ref2, OVSample *dst, int dst_stride,
                uint8_t log2_pb_w, uint8_t log2_pb_h, int8_t pred_mode)
{
    int mode_idx = -(pred_mode - (int)OVINTRA_HOR);
    uint8_t log2_nb_smp = log2_pb_w + log2_pb_h;

    uint8_t use_gauss_filter = log2_nb_smp > 5 && OVABS(mode_idx) > intra_filter[log2_nb_smp >> 1];
    if (use_gauss_filter) {
        intra_angular_gauss_h(rcn_ctx, ref1, ref2, dst, dst_stride, log2_pb_w, log2_pb_h, mode_idx);
    } else {
        intra_angular_cubic_h(rcn_ctx, ref1, ref2, dst, dst_stride, log2_pb_w, log2_pb_h, mode_idx);
    }
}

static void
vvc_intra_pred(const struct OVRCNCtx *const rcn_ctx, const struct OVBuffInfo* ctu_buff,
               uint8_t intra_mode, int x0, int y0,
               int log2_pb_w, int log2_pb_h, CUFlags cu_flags)
{
    const struct DCFunctions *dc = &rcn_ctx->ctudec->rcn_funcs.dc;
    const struct PlanarFunctions *planar = &rcn_ctx->ctudec->rcn_funcs.planar;

    OVSample ref_above[(128 << 1) + 128];
    OVSample ref_left [(128 << 1) + 128];

    ptrdiff_t dst_stride = ctu_buff->stride;

    const OVSample *src = ctu_buff->y;
    OVSample *dst = &ctu_buff->y[x0 + (y0 * dst_stride)];

    OVSample *ref1 = ref_above + (1 << log2_pb_h);
    OVSample *ref2 = ref_left + (1 << log2_pb_w);

    rcn_ctx->ctudec->rcn_funcs.tmp.fill_ref_left_0(src,dst_stride,ref2,
                                           rcn_ctx->progress_field.vfield[x0 >> 2],
                                           rcn_ctx->progress_field.hfield[y0 >> 2],
                                           x0, y0, log2_pb_w, log2_pb_h, 0);

    rcn_ctx->ctudec->rcn_funcs.tmp.fill_ref_above_0(src, dst_stride, ref1,
                                            rcn_ctx->progress_field.hfield[y0 >> 2],
                                            rcn_ctx->progress_field.vfield[x0 >> 2],
                                            x0, y0, log2_pb_w, log2_pb_h, 0);

    if (cu_flags & flg_intra_bdpcm_luma_flag) {

        if (cu_flags & flg_intra_bdpcm_luma_dir) {
            const struct IntraAngularFunctions *cubic_v = rcn_ctx->ctudec->rcn_funcs.intra_angular_cubic_v;
            cubic_v->pure(ref1, dst, dst_stride, log2_pb_w, log2_pb_h);
        } else {
            const struct IntraAngularFunctions *cubic_h = rcn_ctx->ctudec->rcn_funcs.intra_angular_cubic_h;
            cubic_h->pure(ref2, dst, dst_stride, log2_pb_w, log2_pb_h);
        }

        return;
    }

    switch (intra_mode) {
    case OVINTRA_PLANAR:
    {
        OVSample filtered_ref_abv[(128 << 1) + 128];
        OVSample filtered_ref_lft[(128 << 1) + 128];

        if ((log2_pb_h + log2_pb_w) > 5) {
            rcn_ctx->ctudec->rcn_funcs.tmp.filter_ref_samples(ref1, filtered_ref_abv, ref2,
                                                      (1 << log2_pb_w) + 4);
            rcn_ctx->ctudec->rcn_funcs.tmp.filter_ref_samples(ref2, filtered_ref_lft, ref1,
                                                      (1 << log2_pb_h) + 4);
            ref1 = filtered_ref_abv;
            ref2 = filtered_ref_lft;
        }

        planar->pdpc[log2_pb_w-2][log2_pb_h-2](ref1, ref2, dst, dst_stride,
                                                     log2_pb_w, log2_pb_h);
        break;
    }
    case OVINTRA_DC:
    {
        dc->pdpc[log2_pb_w-2][log2_pb_h-2](ref1, ref2, dst, dst_stride, log2_pb_w, log2_pb_h);

        break;
    }
    default:
    {
        int pred_mode = derive_wide_angular_mode(log2_pb_w, log2_pb_h, intra_mode);

        int is_vertical = pred_mode >= OVINTRA_DIA;

        if (is_vertical) {
            intra_angular_v(rcn_ctx, ref1, ref2, dst, dst_stride, log2_pb_w, log2_pb_h, pred_mode);
        } else {
            intra_angular_h(rcn_ctx, ref1, ref2, dst, dst_stride, log2_pb_w, log2_pb_h, pred_mode);
        }
        break;
    }
    }
}

static void
vvc_intra_pred_isp(const OVCTUDec *const ctudec,
                   OVSample *const src,
                   ptrdiff_t dst_stride,
                   uint8_t intra_mode,
                   int x0, int y0,
                   int log2_pb_w, int log2_pb_h,
                   int log2_cb_w, int log2_cb_h,
                   int offset_x, int offset_y)
{

    OVSample ref_abv[(128<<1) + 128];
    OVSample ref_lft[(128<<1) + 128];
    OVSample *dst = &src[x0 + (y0 * dst_stride)];
    OVSample *ref1 = ref_abv + (1 << log2_pb_h);
    OVSample *ref2 = ref_lft + (1 << log2_pb_w);
    const struct OVRCNCtx *const rcn_ctx = &ctudec->rcn_ctx;
    const struct DCFunctions *dc = &ctudec->rcn_funcs.dc;
    const struct PlanarFunctions *planar = &ctudec->rcn_funcs.planar;

    ctudec->rcn_funcs.tmp.fill_ref_left_0(src, dst_stride, ref2,
                                                  ctudec->rcn_ctx.progress_field.vfield[(x0 >> 2) + !!(offset_x % 4)],
                    ctudec->rcn_ctx.progress_field.hfield[((y0 ) >> 2) + !!(y0 % 4)],
                    x0, y0 - offset_y, log2_cb_w, log2_cb_h, offset_y);

    ctudec->rcn_funcs.tmp.fill_ref_above_0(src, dst_stride, ref1,
                     ctudec->rcn_ctx.progress_field.hfield[(y0 >> 2) + !!(offset_y % 4)],
                     ctudec->rcn_ctx.progress_field.vfield[((x0 ) >> 2) + !!(x0 % 4)],
                     x0 - offset_x, y0, log2_cb_w, log2_cb_h, offset_x);

    ref1 += offset_x;
    ref2 += offset_y;

    for (int i = 0; i < 4; ++i){
        ref2[(1 << log2_cb_h) + (1 << log2_pb_h) + 1 + i] = ref2[(1 << log2_cb_h) + (1 << log2_pb_h) + i];
    }

    for (int i = 0; i < 4 ; ++i){
        ref1[(1 << log2_cb_w) + (1 << log2_pb_w) + 1 + i] = ref1[(1 << log2_cb_w) + (1 << log2_pb_w) + i];
    }

    switch (intra_mode) {
    case OVINTRA_PLANAR:
    {
        if (log2_pb_h > 1) {
            planar->pdpc[log2_pb_w-2][log2_pb_h-2](ref1, ref2, dst, dst_stride,
                                                         log2_pb_w, log2_pb_h);
        } else {
            planar->func(ref1, ref2, dst, dst_stride, log2_pb_w, log2_pb_h);
        }
        break;
    }
    case OVINTRA_DC:
    {
        if (log2_pb_h > 1) {
            dc->pdpc[log2_pb_w-2][log2_pb_h-2](ref1, ref2, dst, dst_stride, log2_pb_w, log2_pb_h);
        } else {
            dc->func(ref1, ref2, dst, dst_stride, log2_pb_w, log2_pb_h);
        }
        break;
    }
    default:
    {
        int pred_mode = derive_wide_angular_mode(log2_cb_w, log2_cb_h, intra_mode);

        int is_vertical = pred_mode >= OVINTRA_DIA ? 1 : 0;

        if (is_vertical) {
            int mode_idx = pred_mode - (int)OVINTRA_VER;
            intra_angular_cubic_v(rcn_ctx, ref1, ref2, dst, dst_stride, log2_pb_w, log2_pb_h, mode_idx);
        } else {
            int mode_idx = -(pred_mode - (int)OVINTRA_HOR);
            intra_angular_cubic_h(rcn_ctx, ref1, ref2, dst, dst_stride, log2_pb_w, log2_pb_h, mode_idx);
        }
        break;
    }
    }
}

static void
vvc_intra_pred_multi_ref(const OVCTUDec *const ctudec,
                         OVSample *const src,
                         ptrdiff_t dst_stride,
                         uint8_t intra_mode, int x0, int y0,
                         int log2_pb_w, int log2_pb_h,
                         int mrl_idx)
{
    OVSample ref_abv[(128<<1) + 128];
    OVSample ref_lft [(128<<1) + 128];
    OVSample *dst = &src[x0 + (y0 * dst_stride)];
    OVSample *ref1 = ref_abv + (1 << log2_pb_h);
    OVSample *ref2 = ref_lft  + (1 << log2_pb_w);
    const struct DCFunctions *dc = &ctudec->rcn_funcs.dc;
    const struct PlanarFunctions *planar = &ctudec->rcn_funcs.planar;
    const struct IntraMRLFunctions *mrl_func = ctudec->rcn_funcs.intra_mrl;
    const struct IntraAngularFunctions *nofrac_v = ctudec->rcn_funcs.intra_angular_nofrac_v;
    const struct IntraAngularFunctions *nofrac_h = ctudec->rcn_funcs.intra_angular_nofrac_h;

    ctudec->rcn_funcs.tmp.fill_ref_left_0_mref(src, dst_stride, ref2,
                         ctudec->rcn_ctx.progress_field.vfield[x0 >> 2],
                         ctudec->rcn_ctx.progress_field.hfield[y0 >> 2],
                         mrl_idx, x0, y0,
                         log2_pb_w, log2_pb_h);

    ctudec->rcn_funcs.tmp.fill_ref_above_0_mref(src, dst_stride, ref1,
                          ctudec->rcn_ctx.progress_field.hfield[y0 >> 2],
                          ctudec->rcn_ctx.progress_field.vfield[x0 >> 2],
                          mrl_idx, x0 , y0,
                          log2_pb_w, log2_pb_h);

    ref1 += mrl_idx;
    ref2 += mrl_idx;

    switch (intra_mode) {
    case OVINTRA_PLANAR:
        planar->func(ref1, ref2, dst, dst_stride, log2_pb_w, log2_pb_h);
        break;
    case OVINTRA_DC:
        dc->func(ref1, ref2, dst, dst_stride, log2_pb_w, log2_pb_h);
        break;
    default:
    {
        int pred_mode = derive_wide_angular_mode(log2_pb_w, log2_pb_h, intra_mode);

        int is_vertical = pred_mode >= OVINTRA_DIA ? 1 : 0;

        if(is_vertical){
            int mode_idx = pred_mode - OVINTRA_VER ;
            switch (mode_idx) {
            case 0:
                nofrac_v->pure(ref1, dst, dst_stride, log2_pb_w, log2_pb_h);
                break;
            case 16:
                ref1 += mrl_idx;
                nofrac_v->diagonal(ref1, dst, dst_stride, log2_pb_w, log2_pb_h);
                break;
            default:
            {
                int angle_val;
                if (mode_idx < 0){
                    int inv_angle = inverse_angle_table[-mode_idx];
                    int pb_h  = 1 << log2_pb_h;
                    int inv_angle_sum = 256;

                    OVSample *dst_ref = ref1 - mrl_idx;
                    OVSample *tmp_lft = ref2 - mrl_idx;

                    angle_val = -angle_table[-mode_idx];

                    /* FIXME two stage fill and broadcast last value */
                    for (int k = -1; k >= -pb_h; k--) {
                        inv_angle_sum += inv_angle;
                        dst_ref[k] = tmp_lft[OVMIN(inv_angle_sum >> 9, pb_h)];
                    }

                } else {
                    angle_val = angle_table[mode_idx];
                }

                if (angle_val & 0x1F) {
                    mrl_func->angular_v(ref1, dst, dst_stride, log2_pb_w, log2_pb_h,
                                        angle_val, mrl_idx);
                } else {
                    ref1 += (angle_val * mrl_idx) >> 5;
                    nofrac_v->angular(ref1, dst, dst_stride, log2_pb_w, log2_pb_h,
                                      angle_val);
                }

                break;
            }
            }
        } else {
            int mode_idx = -(pred_mode - OVINTRA_HOR);
            switch (mode_idx) {
            case 0:
                nofrac_h->pure(ref2, dst, dst_stride, log2_pb_w, log2_pb_h);
                break;
            case 16:
                ref2 += mrl_idx;
                nofrac_h->diagonal(ref2, dst, dst_stride, log2_pb_w, log2_pb_h);
                break;
            default:
            {
                int angle_val;
                if (mode_idx < 0) {
                    int inv_angle = inverse_angle_table[-mode_idx];
                    int inv_angle_sum = 256;

                    OVSample *dst_ref = ref2 - mrl_idx;
                    OVSample *tmp_abv = ref1 - mrl_idx;

                    int pb_w = 1 << log2_pb_w;

                    angle_val = -angle_table[-mode_idx];

                    /* FIXME two stage fill and broadcast last value */
                    for (int k = -1; k >= -pb_w; k--) {
                        inv_angle_sum += inv_angle;
                        dst_ref[k] = tmp_abv[OVMIN((inv_angle_sum >> 9), pb_w)];
                    }

                } else {
                    angle_val = angle_table[mode_idx];
                }

                if (angle_val & 0x1F) {
                    mrl_func->angular_h(ref2, dst, dst_stride, log2_pb_w, log2_pb_h,
                                        angle_val, mrl_idx);
                } else {
                    ref2 += (angle_val * mrl_idx) >> 5;
                    nofrac_h->angular(ref2, dst, dst_stride, log2_pb_w, log2_pb_h,
                                      angle_val);
                }
                break;
            }
            }
        }
        break;
    }
    }
}

static void
intra_angular_chroma_v(const struct OVRCNCtx *rcn_ctx, OVSample *ref1, OVSample *ref2, OVSample *dst, int dst_stride,
                       uint8_t log2_pb_w, uint8_t log2_pb_h, int8_t mode_idx)
{
    int inv_angle = inverse_angle_table[OVABS(mode_idx)];
    int abs_angle = angle_table[OVABS(mode_idx)];
    uint8_t req_frac = !!(abs_angle& 0x1F);
    int8_t pdpc_scale = OVMIN(2, log2_pb_h - (floor_log2(3 * inverse_angle_table[OVABS(mode_idx)] - 2) - 8));
    const struct IntraAngularFunctions *nofrac_v = rcn_ctx->ctudec->rcn_funcs.intra_angular_nofrac_v;
    const struct IntraAngularFunctions *c_v = rcn_ctx->ctudec->rcn_funcs.intra_angular_c_v;

    if (mode_idx < 0) {
        int pb_h = 1 << log2_pb_h;
        int inv_angle_sum    = 256;

        for( int k = -1; k >= -pb_h; k-- ) {
            inv_angle_sum += inv_angle;
            ref1[k] = ref2[OVMIN(inv_angle_sum >> 9, pb_h)];
        }

        if (!req_frac) {

            nofrac_v->angular(ref1, dst, dst_stride, log2_pb_w, log2_pb_h,
                              -abs_angle);

        } else {
            c_v->angular(ref1, dst, dst_stride, log2_pb_w, log2_pb_h,
                         -abs_angle);
        }
    } else if (pdpc_scale < 0) {

        if (!req_frac) {

            nofrac_v->angular(ref1, dst, dst_stride, log2_pb_w, log2_pb_h,
                              abs_angle);

        } else {
            c_v->angular(ref1, dst, dst_stride, log2_pb_w, log2_pb_h,
                         abs_angle);
        }

    } else {
        if (!req_frac) {
            if (log2_pb_h > 1 && log2_pb_w > 1 && pdpc_scale >= 0) {
                c_v->angular_pdpc(ref1, ref2, dst, dst_stride,
                                  log2_pb_w, log2_pb_h,
                                  mode_idx);
            } else {

                nofrac_v->angular(ref1, dst, dst_stride, log2_pb_w, log2_pb_h,
                                  abs_angle);

            }
        } else {
            if (log2_pb_h > 1 && log2_pb_w > 1 && pdpc_scale >= 0) {
                c_v->angular_pdpc(ref1, ref2, dst, dst_stride,
                                         log2_pb_w, log2_pb_h,
                                         mode_idx);
            } else {
                c_v->angular(ref1, dst, dst_stride, log2_pb_w, log2_pb_h,
                             abs_angle);
            }
        }
    }
}

static void
intra_angular_chroma_h(const struct OVRCNCtx *rcn_ctx, OVSample *ref1, OVSample *ref2, OVSample *dst, int dst_stride,
                       uint8_t log2_pb_w, uint8_t log2_pb_h, int8_t mode_idx)
{
    int inv_angle = inverse_angle_table[OVABS(mode_idx)];
    int abs_angle = angle_table[OVABS(mode_idx)];
    uint8_t req_frac = !!(abs_angle& 0x1F);
    int8_t pdpc_scale = OVMIN(2, log2_pb_w - (floor_log2(3 * inverse_angle_table[OVABS(mode_idx)] - 2) - 8));
    const struct IntraAngularFunctions *nofrac_h = rcn_ctx->ctudec->rcn_funcs.intra_angular_nofrac_h;
    const struct IntraAngularFunctions *c_h = rcn_ctx->ctudec->rcn_funcs.intra_angular_c_h;

    if (mode_idx < 0) {
        int pb_w = 1 << log2_pb_w;
        int inv_angle_sum    = 256;

        for( int k = -1; k >= -pb_w; k-- ) {
            inv_angle_sum += inv_angle;
            ref2[k] = ref1[OVMIN(inv_angle_sum >> 9, pb_w)];
        }

        if (!req_frac) {

            nofrac_h->angular(ref2, dst, dst_stride, log2_pb_w, log2_pb_h,
                              -abs_angle);

        } else {

            c_h->angular(ref2, dst, dst_stride, log2_pb_w, log2_pb_h,
                              -abs_angle);

        }

    } else if (mode_idx < 8 && pdpc_scale < 0) {
        if (!req_frac) {

            nofrac_h->angular(ref2, dst, dst_stride, log2_pb_w, log2_pb_h,
                              abs_angle);

        } else {

            c_h->angular(ref2, dst, dst_stride, log2_pb_w, log2_pb_h,
                         abs_angle);
        }

    } else {
        if (!req_frac) {
            if (log2_pb_h > 1 && log2_pb_w > 1 && pdpc_scale >= 0) {

                nofrac_h->angular_pdpc(ref1, ref2, dst, dst_stride,
                                       log2_pb_w, log2_pb_h,
                                       mode_idx);

            } else {

                nofrac_h->angular(ref2, dst, dst_stride, log2_pb_w, log2_pb_h,
                                  abs_angle);

            }
        } else {
            if (log2_pb_h > 1 && log2_pb_w > 1 && pdpc_scale >= 0) {

                c_h->angular_pdpc(ref1, ref2, dst, dst_stride,
                                  log2_pb_w, log2_pb_h,
                                  mode_idx);

            } else {

                c_h->angular(ref2, dst, dst_stride, log2_pb_w, log2_pb_h,
                             abs_angle);

            }

        }
    }
}

static void
vvc_intra_chroma_angular(const struct OVRCNCtx *rcn_ctx, const OVSample *const src, OVSample *const dst,
                         OVSample *ref_left, OVSample *ref_above,
                         uint64_t left_col_map, uint64_t top_row_map,
                         int8_t log2_pb_w, int8_t log2_pb_h,
                         int8_t x0, int8_t y0,
                         int8_t intra_mode)
{
    int pred_mode = derive_wide_angular_mode(log2_pb_w, log2_pb_h,
                                             intra_mode);

    int dst_stride = rcn_ctx->ctu_buff.stride_c;
    OVSample *ref1 = ref_above + (1 << log2_pb_h);
    OVSample *ref2 = ref_left + (1 << log2_pb_w);
    int is_vertical = pred_mode >= OVINTRA_DIA ? 1 : 0;

    rcn_ctx->ctudec->rcn_funcs.tmp.fill_ref_left_0_chroma(src, dst_stride, ref2,
                                                  left_col_map, top_row_map,
                                                  x0, y0, log2_pb_w, log2_pb_h);

    rcn_ctx->ctudec->rcn_funcs.tmp.fill_ref_above_0_chroma(src, dst_stride, ref1,
                                                   top_row_map, left_col_map,
                                                   x0, y0, log2_pb_w, log2_pb_h);

    if (is_vertical) {
        const struct IntraAngularFunctions *c_v = rcn_ctx->ctudec->rcn_funcs.intra_angular_c_v;
        int mode_idx = pred_mode - (int)OVINTRA_VER;
        switch (mode_idx) {
            case 0:
                if (log2_pb_h > 1 && log2_pb_w > 1) {
                    c_v->pure_pdpc(ref1, ref2, dst, dst_stride, log2_pb_w, log2_pb_h);
                } else {
                    c_v->pure(ref1, dst, dst_stride, log2_pb_w, log2_pb_h);
                }

                break;
            case (16):
            {
                if (log2_pb_h > 1 && log2_pb_w > 1) {
                    c_v->diagonal_pdpc(ref1, ref2, dst, dst_stride,
                                       log2_pb_w, log2_pb_h);
                } else {
                    c_v->diagonal(ref1, dst, dst_stride, log2_pb_w, log2_pb_h);
                }
            }
                break;
            default:
            {
                intra_angular_chroma_v(rcn_ctx, ref1, ref2, dst, dst_stride, log2_pb_w, log2_pb_h, mode_idx);
            }
            break;
        }
    } else {
        const struct IntraAngularFunctions *c_h = rcn_ctx->ctudec->rcn_funcs.intra_angular_c_h;
        int mode_idx = -(pred_mode - (int)OVINTRA_HOR);
        switch (mode_idx) {
            case 0:
                if (log2_pb_h > 1 && log2_pb_w > 1) {
                    c_h->pure_pdpc(ref1, ref2, dst, dst_stride,
                                           log2_pb_w, log2_pb_h);
                } else {
                    c_h->pure(ref2, dst, dst_stride, log2_pb_w, log2_pb_h);
                }

                break;
            case (16):
            {
                if (log2_pb_h > 1 && log2_pb_w > 1) {
                    c_h->diagonal_pdpc(ref1, ref2, dst, dst_stride,
                                              log2_pb_w, log2_pb_h);
                } else {
                    c_h->diagonal(ref2, dst, dst_stride, log2_pb_w, log2_pb_h);
                }
            }
                break;
            default:
                {
                    intra_angular_chroma_h(rcn_ctx, ref1, ref2, dst, dst_stride, log2_pb_w, log2_pb_h, mode_idx);
                }
                break;
        }
    }
}

static void
vvc_intra_pred_chroma(const struct OVRCNCtx *const rcn_ctx,
                      uint8_t intra_mode, int x0, int y0,
                      int log2_pb_w, int log2_pb_h, CUFlags cu_flags)
{
    const struct RCNFunctions *rcn_func = &rcn_ctx->ctudec->rcn_funcs;
    const struct DCFunctions *dc = &rcn_ctx->ctudec->rcn_funcs.dc;
    const struct PlanarFunctions *planar = &rcn_ctx->ctudec->rcn_funcs.planar;
    OVCTUDec *const ctudec = rcn_ctx->ctudec;

    const struct OVBuffInfo *ctu_buff = &rcn_ctx->ctu_buff;

    OVSample *const dst_cb = &ctu_buff->cb[(x0) + (y0 * ctu_buff->stride_c)];
    OVSample *const dst_cr = &ctu_buff->cr[(x0) + (y0 * ctu_buff->stride_c)];

    const OVSample *const src_cb = &rcn_ctx->ctu_buff.cb[0];
    const OVSample *const src_cr = &rcn_ctx->ctu_buff.cr[0];

    ptrdiff_t dst_stride = ctu_buff->stride_c;

    /*TODO load ref_sample for cb and cr in same function*/
    OVSample ref_above[(128<<1) + 128];
    OVSample ref_left [(128<<1) + 128];
    OVSample *ref1 = ref_above;
    OVSample *ref2 = ref_left;

    uint64_t left_col_map = rcn_ctx->progress_field_c.vfield[x0 >> 1];
    uint64_t top_row_map  = rcn_ctx->progress_field_c.hfield[y0 >> 1];

    if (cu_flags & flg_intra_bdpcm_chroma_flag) {
        if (cu_flags & flg_intra_bdpcm_chroma_dir) {
            const struct IntraAngularFunctions *cubic_v = rcn_ctx->ctudec->rcn_funcs.intra_angular_cubic_v;
            ctudec->rcn_funcs.tmp.fill_ref_above_0_chroma(src_cb, dst_stride, ref_above,
                                                          top_row_map, left_col_map,
                                                          x0, y0, log2_pb_w, log2_pb_h);

            cubic_v->pure(ref1, dst_cb, dst_stride, log2_pb_w, log2_pb_h);

            ctudec->rcn_funcs.tmp.fill_ref_above_0_chroma(src_cr, dst_stride, ref_above,
                                                          top_row_map, left_col_map,
                                                          x0, y0, log2_pb_w, log2_pb_h);

            cubic_v->pure(ref1, dst_cr, dst_stride, log2_pb_w, log2_pb_h);
        } else  {
            const struct IntraAngularFunctions *cubic_h = rcn_ctx->ctudec->rcn_funcs.intra_angular_cubic_h;

            ctudec->rcn_funcs.tmp.fill_ref_left_0_chroma(src_cb, dst_stride, ref_left,
                                                         left_col_map, top_row_map,
                                                         x0, y0, log2_pb_w, log2_pb_h);

            cubic_h->pure(ref2, dst_cb, dst_stride, log2_pb_w, log2_pb_h);

            ctudec->rcn_funcs.tmp.fill_ref_left_0_chroma(src_cr, dst_stride, ref_left,
                                                         left_col_map, top_row_map,
                                                         x0, y0, log2_pb_w, log2_pb_h);


            cubic_h->pure(ref2, dst_cr, dst_stride, log2_pb_w, log2_pb_h);
        }

        return;
    }

    switch (intra_mode) {
    case OVINTRA_PLANAR://PLANAR
    {
        ctudec->rcn_funcs.tmp.fill_ref_left_0_chroma(src_cb, dst_stride, ref_left,
                               left_col_map, top_row_map,
                               x0, y0, log2_pb_w, log2_pb_h);

        ctudec->rcn_funcs.tmp.fill_ref_above_0_chroma(src_cb, dst_stride, ref_above,
                                top_row_map, left_col_map,
                                x0, y0, log2_pb_w, log2_pb_h);

        if (log2_pb_h > 1 && log2_pb_w > 1) {
            planar->pdpc[log2_pb_w-2][log2_pb_h-2](ref1, ref2, dst_cb, dst_stride,
                            log2_pb_w, log2_pb_h);
        } else {
            planar->func(ref1, ref2, dst_cb, dst_stride,
                         log2_pb_w, log2_pb_h);

        }

        ctudec->rcn_funcs.tmp.fill_ref_left_0_chroma(src_cr, dst_stride, ref_left,
                               left_col_map, top_row_map,
                               x0, y0, log2_pb_w, log2_pb_h);

        ctudec->rcn_funcs.tmp.fill_ref_above_0_chroma(src_cr, dst_stride, ref_above,
                                top_row_map, left_col_map,
                                x0, y0, log2_pb_w, log2_pb_h);

        if (log2_pb_h > 1 && log2_pb_w > 1) {
            planar->pdpc[log2_pb_w-2][log2_pb_h-2](ref1, ref2, dst_cr, dst_stride,
                            log2_pb_w, log2_pb_h);
        } else {
            planar->func(ref1, ref2, dst_cr, dst_stride,
                         log2_pb_w, log2_pb_h);
        }
        break;
    }
    case OVINTRA_DC://DC
    {
        ctudec->rcn_funcs.tmp.fill_ref_left_0_chroma(src_cb, dst_stride, ref_left,
                               left_col_map, top_row_map,
                               x0, y0, log2_pb_w, log2_pb_h);

        ctudec->rcn_funcs.tmp.fill_ref_above_0_chroma(src_cb, dst_stride, ref_above,
                                top_row_map, left_col_map,
                                x0, y0, log2_pb_w, log2_pb_h);

        /* PDPC disable for 4xX and Xx4 blocks */
        if (log2_pb_h > 1 && log2_pb_w > 1) {
            dc->pdpc[log2_pb_w-2][log2_pb_h-2](ref1, ref2, dst_cb, dst_stride, log2_pb_w,
                              log2_pb_h);
        } else {
            dc->func(ref1, ref2, dst_cb, dst_stride, log2_pb_w,
                         log2_pb_h);
        }

        ctudec->rcn_funcs.tmp.fill_ref_left_0_chroma(src_cr, dst_stride, ref_left,
                               left_col_map, top_row_map,
                               x0, y0, log2_pb_w, log2_pb_h);

        ctudec->rcn_funcs.tmp.fill_ref_above_0_chroma(src_cr, dst_stride, ref_above,
                                top_row_map, left_col_map,
                                x0, y0, log2_pb_w, log2_pb_h);

        /* PDPC disable for 4xX and Xx4 blocks */
        if (log2_pb_h > 1 && log2_pb_w > 1) {
            dc->pdpc[log2_pb_w-2][log2_pb_h-2](ref1, ref2, dst_cr, dst_stride, log2_pb_w,
                              log2_pb_h);
        } else {
            dc->func(ref1, ref2, dst_cr, dst_stride, log2_pb_w, log2_pb_h);
        }

        break;
    }
    case OVINTRA_LM_CHROMA:
    {
        const OVSample  *const src_luma = &ctu_buff->y[(x0<<1)+((y0<<1)*ctu_buff->stride)];
        /* FIXME to be replaced by progress fields */
        uint8_t neighbour = rcn_ctx->ctudec->ctu_ngh_flags;
        uint8_t got_left_ctu = neighbour & CTU_LFT_FLG;
        uint8_t got_top_ctu  = neighbour & CTU_UP_FLG;

        rcn_func->cclm.cclm(rcn_ctx, log2_pb_w, log2_pb_h, x0, y0,
                            got_top_ctu || y0, got_left_ctu || x0);
        break;
    }
    case OVINTRA_MDLM_LEFT:
    {
        uint8_t neighbour = rcn_ctx->ctudec->ctu_ngh_flags;
        uint8_t got_left_ctu = neighbour & CTU_LFT_FLG;
        uint8_t got_top_ctu  = neighbour & CTU_UP_FLG;

        rcn_func->cclm.mdlm_left(rcn_ctx,
                                 left_col_map, log2_pb_w, log2_pb_h,
                                 x0, y0, x0 || got_left_ctu, y0 || got_top_ctu);
        break;
    }
    case OVINTRA_MDLM_TOP:
    {
        uint8_t neighbour = rcn_ctx->ctudec->ctu_ngh_flags;
        uint8_t got_left_ctu = neighbour & CTU_LFT_FLG;
        uint8_t got_top_ctu  = neighbour & CTU_UP_FLG;

        rcn_func->cclm.mdlm_top(rcn_ctx,
                                top_row_map, log2_pb_w, log2_pb_h,
                                x0, y0, x0 || got_left_ctu, y0 || got_top_ctu);
        break;
    }
    default://angular
    {
        vvc_intra_chroma_angular(rcn_ctx, src_cb, dst_cb, ref_left, ref_above, left_col_map,
                                 top_row_map, log2_pb_w, log2_pb_h,
                                 x0, y0, intra_mode);

        vvc_intra_chroma_angular(rcn_ctx, src_cr, dst_cr, ref_left, ref_above, left_col_map,
                                 top_row_map, log2_pb_w, log2_pb_h,
                                 x0, y0, intra_mode);
        break;

    }
    }

    fill_bs_map(&ctudec->dbf_info.bs2_map_c, x0 << 1, y0 << 1, log2_pb_w + 1, log2_pb_h + 1);

}

void
BD_DECL(rcn_init_intra_functions)(struct RCNFunctions *const rcn_funcs)
{
    rcn_funcs->intra_pred   = &vvc_intra_pred;
    rcn_funcs->intra_pred_c = &vvc_intra_pred_chroma;

    rcn_funcs->intra_pred_isp = &vvc_intra_pred_isp;
    rcn_funcs->intra_pred_mrl = &vvc_intra_pred_multi_ref;
}
