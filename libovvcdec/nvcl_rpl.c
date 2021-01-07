#include "nvcl.h"
#include "nvcl_utils.h"

typedef struct OVRefPicInfo
{
    uint8_t inter_layer_ref_pic_flag;
    uint8_t st_ref_pic_flag;
    uint8_t abs_delta_poc_st;
    uint8_t strp_entry_sign_flag;
    uint8_t rpls_poc_lsb_lt;
    uint8_t ilrp_idx;
} OVRefPicInfo;

typedef struct OVRPLSet 
{
    uint8_t rpl_sps_flag[i];
    uint8_t rpl_idx[i];
    uint8_t poc_lsb_lt[i][j];
    uint8_t delta_poc_msb_cycle_present_flag[i][j];
    uint8_t delta_poc_msb_cycle_lt[i][j];
} OVRPLSet;

typedef struct OVRPL
{
    uint8_t num_ref_entries;
    uint8_t ltrp_in_header_flag;
    uint8_t inter_layer_ref_pic_flag[i];
    uint8_t st_ref_pic_flag[i];
    uint8_t abs_delta_poc_st[i];
    uint8_t strp_entry_sign_flag[i];
    uint8_t rpls_poc_lsb_lt[j++];/*fixme*/
    uint8_t ilrp_idx[i];
} OVRPL;

typedef struct OVRefPic
{
  uint8_t rp_type;
  union {
    uint8_t poc_lsb_lt;
    uint8_t ilrp_idx;
    uint8_t 
  } desc;
} OVRefPic;

/*FIXME separate reading to make parsing more readable */
/* This one is called by ph/sh */
ref_pic_lists(OVNVCLReader *const rdr)
{
    int i;
    for (i = 0; i < 2; i++) {
        if (sps->sps_num_ref_pic_lists[i] > 0 && (i == 0 || (i == 1 && pps->pps_rpl1_idx_present_flag))) {
            rpls->rpl_sps_flag[i] = nvcl_read_flag(rdr);
        }

        if (rpls->rpl_sps_flag[i]) {
            if (sps->sps_num_ref_pic_lists[i] > 1 && (i == 0 || (i == 1 && pps->pps_rpl1_idx_present_flag))) {
                /*TODO v = ceil log2 sps_num-ref_pic_lists*/
                int v = 0;
                rpls->rpl_idx[i] = nvcl_read_bits(rdr, v);
            }
        } else {
            ref_pic_list_struct(i, sps->sps_num_ref_pic_lists[i]);
        }

        for (j = 0; j < NumLtrpEntries[i][RplsIdx[i]]; j++) {
            if (ltrp_in_header_flag[i][RplsIdx[i]]) {
                rpls->poc_lsb_lt[i][j];
            }

            rpls->delta_poc_msb_cycle_present_flag[i][j] = nvcl_read_flag(rdr);
            if (delta_poc_msb_cycle_present_flag[i][j]) {
                rpls->delta_poc_msb_cycle_lt[i][j] = nvcl_read_u_expgolomb(rdr);
            }
        }
    }
}

#if 0
/*TODO separate ilrp  lt and st cases */
ref_pic_list_struct(listIdx, rplsIdx)
{
    int i, j;
    rpl->num_ref_entries[listIdx][rplsIdx] = nvcl_read_u_expgolomb(rdr);
    if (sps->sps_long_term_ref_pics_flag && rplsIdx < sps->sps_num_ref_pic_lists[listIdx] && num_ref_entries[listIdx][rplsIdx] > 0) {
        rpl->ltrp_in_header_flag[listIdx][rplsIdx] = nvcl_read_flag(rdr);
    }

    for (i = 0, j = 0; i < num_ref_entries[listIdx][rplsIdx]; i++) {

        if (sps->sps_inter_layer_prediction_enabled_flag) {
            rpl->inter_layer_ref_pic_flag[listIdx][rplsIdx][i] = nvcl_read_flag(rdr);
        }

        if (!inter_layer_ref_pic_flag[listIdx][rplsIdx][i]) {

            if (sps->sps_long_term_ref_pics_flag) {
                rpl->st_ref_pic_flag[listIdx][rplsIdx][i] = nvcl_read_flag(rdr);
            }

            if (st_ref_pic_flag[listIdx][rplsIdx][i]) {

                rpl->abs_delta_poc_st[listIdx][rplsIdx][i] = nvcl_read_u_expgolomb(rdr);
                if (AbsDeltaPocSt[listIdx][rplsIdx][i] > 0) {
                    rpl->strp_entry_sign_flag[listIdx][rplsIdx][i];
                }

            } else if (!ltrp_in_header_flag[listIdx][rplsIdx]) {
                rpl->rpls_poc_lsb_lt[listIdx][rplsIdx][j++];
            }
        } else {
            rpl->ilrp_idx[listIdx][rplsIdx][i] = nvcl_read_u_expgolomb(rdr);
        }
    }
}
#endif

static int ref_pic_list_ilrp_ltrp(OVNVCLReader *const rdr)
{
    if (rplsIdx < sps->sps_num_ref_pic_lists[listIdx] && rpl->num_ref_entries > 0) {
        rpl->ltrp_in_header_flag = nvcl_read_flag(rdr);
    }

    for (i = 0, j = 0; i < rpl->num_ref_entries; i++) {

        rpl->inter_layer_ref_pic_flag[i] = nvcl_read_flag(rdr);
        if (rpl->inter_layer_ref_pic_flag[i]) {
            rpl->ilrp_idx[i] = nvcl_read_u_expgolomb(rdr);
        } else {
            if (sps->sps_long_term_ref_pics_flag) {
                rpl->st_ref_pic_flag[i] = nvcl_read_flag(rdr);
            }

            if (st_ref_pic_flag[i]) {

                rpl->abs_delta_poc_st[i] = nvcl_read_u_expgolomb(rdr);
                if (AbsDeltaPocSt[i] > 0) {
                    rpl->strp_entry_sign_flag[i];
                }

            } else if (!rpl->ltrp_in_header_flag) {
                rpl->rpls_poc_lsb_lt[j++];
            }
        }
    }
}

static int ref_pic_list_ilrp(OVNVCLReader *const rdr)
{
    for (i = 0, j = 0; i < rpl->num_ref_entries; i++) {

        rpl->inter_layer_ref_pic_flag[i] = nvcl_read_flag(rdr);
        if (rpl->inter_layer_ref_pic_flag[i]) {
            rpl->ilrp_idx[i] = nvcl_read_u_expgolomb(rdr);
        } else {
            rpl->abs_delta_poc_st[i] = nvcl_read_u_expgolomb(rdr);
            if (AbsDeltaPocSt[i] > 0) {
                rpl->strp_entry_sign_flag[i];
            }
        }
    }
}

static int ref_pic_list_ltrp(OVNVCLReader *const rdr)
{
    if (rplsIdx < sps->sps_num_ref_pic_lists[listIdx] && rpl->num_ref_entries > 0) {
        rpl->ltrp_in_header_flag = nvcl_read_flag(rdr);
        if (rpl->ltrp_in_header_flag) {
            for (i = 0, j = 0; i < num_ref_entries; i++) {
                VVCRefPic *ref_pic = rpl->ref_pic[i];
                rpl->st_ref_pic_flag[i] = nvcl_read_flag(rdr);
                if (st_ref_pic_flag[i]) {
                    rpl->abs_delta_poc_st[i] = nvcl_read_u_expgolomb(rdr);
                    if (AbsDeltaPocSt[i] > 0) {
                        rpl->strp_entry_sign_flag[i];
                    }
                }
            }
        } else {
            for (i = 0, j = 0; i < num_ref_entries; i++) {
                rpl->st_ref_pic_flag[i] = nvcl_read_flag(rdr);
                if (st_ref_pic_flag[i]) {
                    rpl->abs_delta_poc_st[i] = nvcl_read_u_expgolomb(rdr);
                    if (AbsDeltaPocSt[i] > 0) {
                        rpl->strp_entry_sign_flag[i];
                    }
                } else {
                    rpl->rpls_poc_lsb_lt[j++];
                }
            }
        }
    } else {
        for (i = 0, j = 0; i < rpl->num_ref_entries; i++) {
            rpl->st_ref_pic_flag[i] = nvcl_read_flag(rdr);
            if (st_ref_pic_flag[i]) {
                rpl->abs_delta_poc_st[i] = nvcl_read_u_expgolomb(rdr);
                if (AbsDeltaPocSt[i] > 0) {
                    rpl->strp_ent[i];
                }
            } else {
                rpl->rpls_poc_lsb_lt[j++];
            }
        }
    }
}

static int ref_pic_list_strp(OVNVCLReader *const rdr)
{
    for (i = 0, j = 0; i < rpl->num_ref_entries; i++) {
        rpl->abs_delta_poc_st[i] = nvcl_read_u_expgolomb(rdr);
        if (AbsDeltaPocSt[i] > 0) {
            rpl->strp_entry_sign_flag[i];
        }
    }
}

static int ref_pic_list_struct2(OVNVCLReader *const rdr,listIdx, rplsIdx)
{
    int i, j;

    OVRPL *rpl = ref_pic_list_set[listIdx][rplsIdx];

    rpl->num_ref_entries = nvcl_read_u_expgolomb(rdr);

    if (sps->sps_inter_layer_prediction_enabled_flag && sps->sps_long_term_ref_pics_flag) {
        ref_pic_list_ilrp_ltrp();
    } else if (sps->sps_inter_layer_prediction_enabled_flag) {
        ref_pic_list_ilrp();
    } else if (sps->sps_long_term_ref_pics_flag) {
        ref_pic_list_ltrp();
    } else {
        ref_pic_list_strp();
    }
}

static int ref_pic_lists2(OVNVCLReader *const rdr)
{
    int i;
    /* l0 list */
    if (sps->sps_num_ref_pic_lists[0] > 0) {
        rpls->rpl_sps_flag[0] = nvcl_read_flag(rdr);
    }

    if (rpls->rpl_sps_flag[0]) {
        if (sps->sps_num_ref_pic_lists[0] > 1) {
            /*TODO v = ceil log2 sps_num-ref_pic_lists*/
            int v = 0;
            rpls->rpl_idx[0] = nvcl_read_bits(rdr, v);
        }
    } else {
        ref_pic_list_struct(0, sps->sps_num_ref_pic_lists[0]);
    }

    for (j = 0; j < NumLtrpEntries[0][RplsIdx[0]]; j++) {
        if (ltrp_in_header_flag[0][RplsIdx[0]]) {
            rpls->poc_lsb_lt[0][j];
        }

        rpls->delta_poc_msb_cycle_present_flag[0][j] = nvcl_read_flag(rdr);
        if (delta_poc_msb_cycle_present_flag[0][j]) {
            rpls->delta_poc_msb_cycle_lt[0][j] = nvcl_read_u_expgolomb(rdr);
        }
    }

    /* l1 list */
    if (sps->sps_num_ref_pic_lists[1] > 0 && ((pps->pps_rpl1_idx_present_flag))) {
        rpls->rpl_sps_flag[1] = nvcl_read_flag(rdr);
    }

    if (rpls->rpl_sps_flag[1]) {
        if (sps->sps_num_ref_pic_lists[1] > 1 && ((pps->pps_rpl1_idx_present_flag))) {
            /*TODO v = ceil log2 sps_num-ref_pic_lists*/
            int v = 0;
            rpls->rpl_idx[1] = nvcl_read_bits(rdr, v);
        }
    } else {
        ref_pic_list_struct(1, sps->sps_num_ref_pic_lists[1]);
    }

    for (j = 0; j < NumLtrpEntries[1][RplsIdx[1]]; j++) {
        if (ltrp_in_header_flag[1][RplsIdx[1]]) {
            rpls->poc_lsb_lt[1][j];
        }

        rpls->delta_poc_msb_cycle_present_flag[1][j] = nvcl_read_flag(rdr);
        if (delta_poc_msb_cycle_present_flag[1][j]) {
            rpls->delta_poc_msb_cycle_lt[1][j] = nvcl_read_u_expgolomb(rdr);
        }
    }
}
