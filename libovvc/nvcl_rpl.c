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

#include <string.h>
#include "nvcl.h"
#include "nvcl_utils.h"
#include "nvcl_structures.h"
#include "nvcl_private.h"

#if 0
struct RPL{
    /* Storage for convenience variables non signaled as syntax
     * elements but mandatory to derive other syntax elements
     * without recomputing it from scratch
     */
    struct {
        /* Number of long term ref */
        uint8_t nb_ltrp;

        /* Number of short term ref */
        uint8_t nb_strp;

        /* Number of inter layer ref pictures*/
        uint8_t nb_ilrp;
    } helper;
};
#endif

static int
ref_pic_list_ilrp_ltrp(OVNVCLReader *const rdr, struct OVRPL *const rpl,
                       const OVSPS *const sps)
{
    if (!rpl->ltrp_in_header_flag) {
        int i;
        for (i = 0; i <= ((rpl->num_ref_entries - 1) & 0xF); i++) {
            struct RefPic *rp = &rpl->rp_list[i];
            rp->inter_layer_ref_pic_flag = nvcl_read_flag(rdr);
            if (rp->inter_layer_ref_pic_flag) {
                rp->ilrp_idx = nvcl_read_u_expgolomb(rdr);
            } else {
                rp->st_ref_pic_flag = nvcl_read_flag(rdr);
                if (rp->st_ref_pic_flag) {
                    rp->abs_delta_poc_st = nvcl_read_u_expgolomb(rdr);
                    if (rp->abs_delta_poc_st > 0 || (i == 0)) {
                        rp->strp_entry_sign_flag = nvcl_read_flag(rdr);
                    }
                } else {
                    const uint8_t nb_bits = sps->sps_log2_max_pic_order_cnt_lsb_minus4 + 4;
                    rp->rpls_poc_lsb_lt = nvcl_read_bits(rdr, nb_bits);
                }
            }
        }
    } else {
        int i;
        for (i = 0; i <= ((rpl->num_ref_entries - 1) & 0xF); i++) {
            struct RefPic *rp = &rpl->rp_list[i];
            rp->inter_layer_ref_pic_flag = nvcl_read_flag(rdr);
            if (rp->inter_layer_ref_pic_flag) {
                rp->ilrp_idx = nvcl_read_u_expgolomb(rdr);
            } else {
                rp->st_ref_pic_flag = nvcl_read_flag(rdr);
                if (rp->st_ref_pic_flag) {
                    rp->abs_delta_poc_st = nvcl_read_u_expgolomb(rdr);
                    if (rp->abs_delta_poc_st > 0 || (i == 0)) {
                        rp->strp_entry_sign_flag = nvcl_read_flag(rdr);
                    }
                }
            }
        }
    }
    return 0;
}

static int
ref_pic_list_ilrp(OVNVCLReader *const rdr, struct OVRPL *const rpl)
{
    int i;
    for (i = 0; i <= ((rpl->num_ref_entries - 1) & 0xF); i++) {
        struct RefPic *rp = &rpl->rp_list[i];
        rp->inter_layer_ref_pic_flag = nvcl_read_flag(rdr);
        if (rp->inter_layer_ref_pic_flag) {
            rp->ilrp_idx = nvcl_read_u_expgolomb(rdr);
        } else {
            rp->abs_delta_poc_st = nvcl_read_u_expgolomb(rdr);
            if (rp->abs_delta_poc_st > 0 || (i == 0)) {
                rp->strp_entry_sign_flag = nvcl_read_flag(rdr);
            }
        }
    }
    return 0;
}

static int
ref_pic_list_ltrp(OVNVCLReader *const rdr, struct OVRPL *const rpl,
                  const OVSPS *const sps)
{
    if (rpl->ltrp_in_header_flag) {
        /* Read only short term pic info
         * Long term POC LSB will be read outside of
         * the ref pic list
         */
         int i;
        for (i = 0; i <= ((rpl->num_ref_entries - 1) & 0xF); i++) {
            struct RefPic *rp = &rpl->rp_list[i];
            rp->st_ref_pic_flag = nvcl_read_flag(rdr);
            if (rp->st_ref_pic_flag) {
                rp->abs_delta_poc_st = nvcl_read_u_expgolomb(rdr);
                if (rp->abs_delta_poc_st > 0 || (i == 0)) {
                    rp->strp_entry_sign_flag = nvcl_read_flag(rdr);
                }
            }
        }
    } else {
        /* Long term POC LSB info is present in structure
         * Read it if ref pic is long term.
         * Note we can pass here only if ref pic list is
         * read from sps
         */
         int i;
        for (i = 0; i <= ((rpl->num_ref_entries - 1) & 0xF); i++) {
            struct RefPic *rp = &rpl->rp_list[i];
            rp->st_ref_pic_flag = nvcl_read_flag(rdr);
            if (rp->st_ref_pic_flag) {
                rp->abs_delta_poc_st = nvcl_read_u_expgolomb(rdr);
                if (rp->abs_delta_poc_st > 0 || (i == 0)) {
                    rp->strp_entry_sign_flag = nvcl_read_flag(rdr);
                }
            } else {
                const uint8_t nb_bits = sps->sps_log2_max_pic_order_cnt_lsb_minus4 + 4;
                rp->rpls_poc_lsb_lt = nvcl_read_bits(rdr, nb_bits);
            }
        }
    }
    return 0;
}

static int
ref_pic_list_strp(OVNVCLReader *const rdr, struct OVRPL *const rpl, uint8_t use_weighted_pred)
{
    int i;
    for (i = 0; i <= ((rpl->num_ref_entries - 1) & 0xF); i++) {
        struct RefPic *rp = &rpl->rp_list[i];
        /* Default value to 1 to avoid confusion between
         * short and long reference pictures when reading
         * long term data from a header (SH or PH)
         */
        rp->st_ref_pic_flag = 1;
        rp->abs_delta_poc_st = nvcl_read_u_expgolomb(rdr);
        if (rp->abs_delta_poc_st > 0 || (!use_weighted_pred) || (i == 0)) {
            rp->strp_entry_sign_flag = nvcl_read_flag(rdr);
        }
    }
    return 0;
}

#if 0
static int
ref_pic_list_struct2(OVNVCLReader *const rdr,listIdx, rplsIdx)
{
    int i;

    OVRPL *rpl = ref_pic_list_set[listIdx][rplsIdx];

    rpl->num_ref_entries = nvcl_read_u_expgolomb(rdr);

    /* FIXME use swithc and status based on sps flags? */
    if (sps->sps_inter_layer_prediction_enabled_flag && sps->sps_long_term_ref_pics_flag) {
        /* default to 1 */
        rpl->ltrp_in_header_flag = 1;

        /* FIXME my guess is this means we are reading list from SPS
         * since rplsIdx is < to the number of ref pic lists in sps
         * and this function has been called (meaning we either come
         * directly from SPS reader or if called from the PH header
         * rpl_sps_flag was OFF so we did not read an rpls_idx
         */
        if (rplsIdx < sps->sps_num_ref_pic_lists[listIdx] && rpl->num_ref_entries > 0) {
            rpl->ltrp_in_header_flag = nvcl_read_flag(rdr);
        }

        ref_pic_list_ilrp_ltrp();
    } else if (sps->sps_inter_layer_prediction_enabled_flag) {
        ref_pic_list_ilrp();
    } else if (sps->sps_long_term_ref_pics_flag) {
        /* default to 1 */
        rpl->ltrp_in_header_flag = 1;

        if (rplsIdx < sps->sps_num_ref_pic_lists[listIdx] && rpl->num_ref_entries > 0) {
            rpl->ltrp_in_header_flag = nvcl_read_flag(rdr);
        }

        ref_pic_list_ltrp();
    } else {
        ref_pic_list_strp();
    }
}
#endif

static int
ref_pic_list_header(OVNVCLReader *const rdr, const OVSPS *const sps,
                    struct OVRPL *const rpl)
{
    uint8_t use_weighted_pred = sps->sps_weighted_pred_flag || sps->sps_weighted_bipred_flag;

    rpl->num_ref_entries = nvcl_read_u_expgolomb(rdr);

    if (rpl->num_ref_entries > 0) {
        /* FIXME use swithc and status based on sps flags? */
        if (sps->sps_long_term_ref_pics_flag) {
            rpl->ltrp_in_header_flag = 1;
            if (sps->sps_inter_layer_prediction_enabled_flag) {
                ref_pic_list_ilrp_ltrp(rdr, rpl, sps);
            } else {
                ref_pic_list_ltrp(rdr, rpl, sps);
            }
        } else {
            /* No need of ltrp_in_header_flag */
            if (sps->sps_inter_layer_prediction_enabled_flag) {
                ref_pic_list_ilrp(rdr, rpl);
            } else {
                ref_pic_list_strp(rdr, rpl, use_weighted_pred);
            }
        }
    }
    return 0;
}

/* Read Additional info for long term ref pictures */
/* FIXME :
 *   -Only read when some LTRP exist in the selected RPL
 *   -Retrieve correct RPL when used from SPS
 *   -Decide if we should separate LT ST and ILRP Ref Pic Lists
 */
static int
header_read_long_term_info(OVNVCLReader *const rdr, const struct OVRPL *const rpl,
                           struct RPLHeader *rpl_h, const OVSPS *const sps)
{
    if (rpl->ltrp_in_header_flag) {
        int j;
        const uint8_t nb_bits = sps->sps_log2_max_pic_order_cnt_lsb_minus4 + 4;
        for (j = 0; j <= ((rpl->num_ref_entries - 1) & 0xF); j++) {
            const struct RefPic *rp = &rpl->rp_list[j];
            if (!rp->st_ref_pic_flag && !rp->inter_layer_ref_pic_flag) {
                struct LTInfo *lti = &rpl_h->lt_info[j];

                lti->poc_lsb_lt = nvcl_read_bits(rdr, nb_bits);

                lti->delta_poc_msb_cycle_present_flag = nvcl_read_flag(rdr);

                if (lti->delta_poc_msb_cycle_present_flag) {
                    lti->delta_poc_msb_cycle_lt = nvcl_read_u_expgolomb(rdr);
                }
            }
        }
    } else {
        /* This can only be called when rpl comes from sps */
        int j;
        for (j = 0; j <= ((rpl->num_ref_entries - 1) & 0xF); j++) {
            const struct RefPic *rp = &rpl->rp_list[j];
            /* FIXME avoid recount of lt_ref */
            if (!rp->st_ref_pic_flag && !rp->inter_layer_ref_pic_flag) {
                struct LTInfo *lti = &rpl_h->lt_info[j];

                lti->delta_poc_msb_cycle_present_flag = nvcl_read_flag(rdr);

                if (lti->delta_poc_msb_cycle_present_flag) {
                    lti->delta_poc_msb_cycle_lt = nvcl_read_u_expgolomb(rdr);
                }
            }
        }
    }
    return 0;
}

/* This one is called by PH/SH reader */
int
nvcl_read_header_ref_pic_lists(OVNVCLReader *const rdr, OVHRPL *const rpl_h,
                               const OVSPS *const sps, const OVPPS *pps)
{
    /* TODO move outside */
    struct RPLHeader *const rpl_h0 = &rpl_h->rpl_h0;
    struct RPLHeader *const rpl_h1 = &rpl_h->rpl_h1;

    if (sps->sps_num_ref_pic_lists0 > 0) {
        rpl_h0->rpl_sps_flag = nvcl_read_flag(rdr);
    }

    if (rpl_h0->rpl_sps_flag) {
        const struct OVRPL *rpl0;
        if (sps->sps_num_ref_pic_lists0 > 1) {
            int nb_bits = ov_ceil_log2(sps->sps_num_ref_pic_lists0);
            rpl_h0->rpl_idx = nvcl_read_bits(rdr, nb_bits);
        }

        rpl0 = &sps->rpl_s0[rpl_h0->rpl_idx & 0x3F];

        if (sps->sps_long_term_ref_pics_flag) {
            /* Call long term post function */
            header_read_long_term_info(rdr, rpl0, rpl_h0, sps);
        }
        /* FIXME update long term informations */
        memcpy(&rpl_h0->rpl_data, rpl0, sizeof(*rpl0));
    } else {
        struct OVRPL *rpl0 = &rpl_h0->rpl_data;

        ref_pic_list_header(rdr, sps, rpl0);
        if (sps->sps_long_term_ref_pics_flag) {
            /* Call long term post function with lt_header*/
            header_read_long_term_info(rdr, rpl0, rpl_h0, sps);
        }
        //memcpy(&rpl_h0->rpl_data, rpl0, sizeof(*rpl0));
    }
    rpl_h->rpl0 = &rpl_h0->rpl_data;

    /* l1 list */
    if (sps->sps_num_ref_pic_lists1 > 0 && pps->pps_rpl1_idx_present_flag) {
        rpl_h1->rpl_sps_flag = nvcl_read_flag(rdr);
    }

    if (rpl_h1->rpl_sps_flag) {
        /* FIXME pps_rpl_idx_presetn already tested by rpl_sps_flag */
        const struct OVRPL *rpl1;
        if (sps->sps_num_ref_pic_lists1 > 1) {
            int v = ov_ceil_log2(sps->sps_num_ref_pic_lists1);
            rpl_h1->rpl_idx = nvcl_read_bits(rdr, v);
        }

        rpl1 = &sps->rpl_s1[rpl_h1->rpl_idx & 0x3F];

        if (sps->sps_long_term_ref_pics_flag) {
            /* Call long term post function */
            header_read_long_term_info(rdr, rpl1, rpl_h1, sps);
        }
        /* FIXME update long term informations */
        memcpy(&rpl_h1->rpl_data, rpl1, sizeof(*rpl1));

    } else if (rpl_h0->rpl_sps_flag && sps->sps_num_ref_pic_lists1 > 0) {
        /*FIXME check nb_rpl1  + long term poc complement */
        const struct OVRPL *rpl1 = &sps->rpl_s1[rpl_h0->rpl_idx & 0x3F];

        if (sps->sps_long_term_ref_pics_flag) {
            /* Call long term post function with lt_header*/
            header_read_long_term_info(rdr, rpl1, rpl_h1, sps);
        }
        memcpy(&rpl_h1->rpl_data, rpl1, sizeof(*rpl1));
    } else if (rpl_h0->rpl_sps_flag && sps->sps_rpl1_same_as_rpl0_flag) {
        memcpy(&rpl_h1->rpl_data, rpl_h->rpl0, sizeof(*rpl_h->rpl0));
    } else {
        struct OVRPL *rpl1 = &rpl_h1->rpl_data;

        ref_pic_list_header(rdr, sps, rpl1);

        if (sps->sps_long_term_ref_pics_flag) {
            /* Call long term post function with lt_header*/
            header_read_long_term_info(rdr, rpl1, rpl_h1, sps);
        }
    }
        rpl_h->rpl1 = &rpl_h1->rpl_data;

    return 0;
}

int
nvcl_read_sps_ref_pic_list(OVNVCLReader *const rdr, const OVSPS *const sps,
                           OVRPL *const rpl)
{
    uint8_t use_weighted_pred = sps->sps_weighted_pred_flag || sps->sps_weighted_bipred_flag;
    rpl->num_ref_entries = nvcl_read_u_expgolomb(rdr);

    if (rpl->num_ref_entries > 0) {
        /* FIXME use swithc and status based on sps flags? */
        if (sps->sps_long_term_ref_pics_flag) {
            rpl->ltrp_in_header_flag = nvcl_read_flag(rdr);
            if (sps->sps_inter_layer_prediction_enabled_flag) {
                ref_pic_list_ilrp_ltrp(rdr, rpl, sps);
            } else {
                ref_pic_list_ltrp(rdr, rpl, sps);
            }
        } else {
            /* No need of ltrp_in_header_flag */
            if (sps->sps_inter_layer_prediction_enabled_flag) {
                ref_pic_list_ilrp(rdr, rpl);
            } else {
                ref_pic_list_strp(rdr, rpl, use_weighted_pred);
            }
        }
    }
    return 0;
}

