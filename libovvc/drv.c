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


static void
fill_mvp_map(struct VVCMVCtx *const mv_ctx, VVCMV mv,
             int pb_x, int pb_y, int n_pb_w, int n_pb_h)
{
    int i, j;
    for (j = 0; j < n_pb_h; ++j) {
        for (i = 0; i < n_pb_w; ++i) {
            memcpy(&mv_ctx->mvs[PB_POS_IN_BUF(pb_x + i, pb_y + j)], &mv, sizeof(VVCMV));
        }
    }
    update_availability_maps(&mv_ctx->map, pb_x, pb_y, n_pb_w, n_pb_h);
}

#define LF_MV_THRESHOLD 8
static void
fill_dbf_mv_map_b(struct DBFInfo *const dbf_info, struct VVCMVCtx *const mv_ctx, struct VVCMVCtx *const mv_ctx1, VVCMV mv,
                  int pb_x, int pb_y, int n_pb_w, int n_pb_h)
{
    int i, j;
    int log2_diff_min_cu = 1;
    int mask = (1 << (log2_diff_min_cu + 1)) - 1;
    int shift_v = 1 + (pb_y << log2_diff_min_cu);
    int shift_h = 2 + (pb_x << log2_diff_min_cu);

    uint64_t val = dbf_info->bs1_map.hor[(pb_y << log2_diff_min_cu)];

    uint64_t tmp_mask_h = (uint64_t)mask << shift_h;
    uint64_t tmp_mask_v = (uint64_t)mask << shift_v;

    for (j = 0; j < n_pb_w; ++j) {
        VVCMV mv_above = mv_ctx->mvs[PB_POS_IN_BUF(pb_x + j, pb_y - 1)];
        int64_t above_avail = -((!!(mv_ctx->map.rows[pb_y]  & POS_MASK(pb_x + j, 0))
                                & !(mv_ctx1->map.rows[pb_y] & POS_MASK(pb_x + j, 0))));
        int64_t abv_th = -((FFABS(mv_above.x - mv.x) >= LF_MV_THRESHOLD) |
                           (FFABS(mv_above.y - mv.y) >= LF_MV_THRESHOLD));
        val |= (tmp_mask_h & abv_th & above_avail) | (tmp_mask_h & (-(!above_avail)));
        tmp_mask_h  <<= (1 << log2_diff_min_cu);
    }
    dbf_info->bs1_map.hor[(pb_y << log2_diff_min_cu)] |= val;

    val = dbf_info->bs1_map.ver[(pb_x << log2_diff_min_cu)];

    for (i = 0; i < n_pb_h; ++i) {
        VVCMV mv_left = mv_ctx->mvs[PB_POS_IN_BUF(pb_x - 1, pb_y + i)];
        int64_t left_avail = -(!!(mv_ctx->map.cols[pb_x]  & POS_MASK(pb_y + i, 0))
                              & !(mv_ctx1->map.cols[pb_x] & POS_MASK(pb_y + i, 0)));
        int64_t abv_th = -((FFABS(mv_left.x - mv.x) >= LF_MV_THRESHOLD) |
                           (FFABS(mv_left.y - mv.y) >= LF_MV_THRESHOLD));
        val |= (tmp_mask_v & abv_th & left_avail) | (tmp_mask_v & (-(!left_avail)));
        tmp_mask_v <<= (1 << log2_diff_min_cu);
    }
    dbf_info->bs1_map.ver[(pb_x << log2_diff_min_cu)] |= val;
}

static void
fill_dbf_mv_map(struct DBFInfo *const dbf_info, struct VVCMVCtx *const mv_ctx, VVCMV mv,
                int pb_x, int pb_y, int n_pb_w, int n_pb_h)
{
    int i, j;
    int log2_diff_min_cu = 1;
    int mask = (1 << (log2_diff_min_cu + 1)) - 1;
    int shift_v = 1 + (pb_y << log2_diff_min_cu);
    int shift_h = 2 + (pb_x << log2_diff_min_cu);

    uint64_t val = dbf_info->bs1_map.hor[(pb_y << log2_diff_min_cu)];

    uint64_t tmp_mask_h = (uint64_t)mask << shift_h;
    uint64_t tmp_mask_v = (uint64_t)mask << shift_v;
    for (j = 0; j < n_pb_w; ++j) {
        VVCMV mv_above = mv_ctx->mvs[PB_POS_IN_BUF(pb_x + j, pb_y - 1)];
        int64_t above_avail = -(!!(mv_ctx->map.rows[pb_y] & POS_MASK(pb_x + j, 0)));
        int64_t abv_th = -((FFABS(mv_above.x - mv.x) >= LF_MV_THRESHOLD) |
                           (FFABS(mv_above.y - mv.y) >= LF_MV_THRESHOLD));
        val |= (tmp_mask_h & abv_th & above_avail) | (tmp_mask_h & (-(!above_avail)));
        tmp_mask_h  <<= (1 << log2_diff_min_cu);
    }
    dbf_info->bs1_map.hor[(pb_y << log2_diff_min_cu)] |= val;

    val = dbf_info->bs1_map.ver[(pb_x << log2_diff_min_cu)];

    for (i = 0; i < n_pb_h; ++i) {
        VVCMV mv_left = mv_ctx->mvs[PB_POS_IN_BUF(pb_x - 1, pb_y + i)];
        int64_t left_avail = -(!!(mv_ctx->map.cols[pb_x] & POS_MASK(pb_y + i, 0)));
        int64_t abv_th = -((FFABS(mv_left.x - mv.x) >= LF_MV_THRESHOLD) |
                           (FFABS(mv_left.y - mv.y) >= LF_MV_THRESHOLD));
        val |= (tmp_mask_v & abv_th & left_avail) | (tmp_mask_v & (-(!left_avail)));
        tmp_mask_v <<= (1 << log2_diff_min_cu);
    }
    dbf_info->bs1_map.ver[(pb_x << log2_diff_min_cu)] |= val;
}

static VVCMergeInfo
derive_mvp_b(VVCLocalContext *const lc_ctx,
             const VVCPartSize *const part_ctx,
             unsigned int x0, unsigned int y0,
             unsigned int log2_pb_w, unsigned int log2_pb_h,
             VVCMV mvd0, VVCMV mvd1,
             uint8_t mvp_idx0, uint8_t mvp_idx1,
             uint8_t inter_dir)
{
    struct VVCInterCtx *const inter_ctx = &lc_ctx->inter_ctx;
    uint8_t y_pu = y0 >> part_ctx->log2_min_cb_s;
    uint8_t x_pu = x0 >> part_ctx->log2_min_cb_s;
    uint8_t nb_pb_w = (1 << log2_pb_w) >> part_ctx->log2_min_cb_s;
    uint8_t nb_pb_h = (1 << log2_pb_h) >> part_ctx->log2_min_cb_s;
    VVCMV mv0 = {0}, mv1 = {0};
    VVCMergeInfo mv_info;

    if (inter_dir & 0x1) {
        struct VVCMVCtx *const mv_ctx0 = &inter_ctx->mv_ctx0;

        mv0 = derive_mvp_candidates(lc_ctx, inter_ctx, mv_ctx0,
                                    x_pu, y_pu, nb_pb_w, nb_pb_h,
                                    mvp_idx0, inter_dir & 0x1);

        mvd0 = scale_mvd(mvd0);

        mv0.x += mvd0.x;
        mv0.y += mvd0.y;

    }

    if (inter_dir & 0x2) {
        struct VVCMVCtx *const mv_ctx1 = &inter_ctx->mv_ctx1;

        mv1 = derive_mvp_candidates(lc_ctx, inter_ctx, mv_ctx1,
                                    x_pu, y_pu,
                                    nb_pb_w, nb_pb_h,
                                    mvp_idx1, inter_dir & 0x2);

        mvd1 = scale_mvd(mvd1);

        mv1.x += mvd1.x;
        mv1.y += mvd1.y;

    }

    mv_info.inter_dir = inter_dir;
    mv_info.mv0 = mv0;
    mv_info.mv1 = mv1;

    return mv_info;
}

static void
update_mv_ctx_b(VVCLocalContext *const lc_ctx, struct VVCInterCtx *const inter_ctx,
                const VVCPartSize *const part_ctx,
                const VVCMV mv0, const VVCMV mv1,
                unsigned int x0, unsigned int y0,
                unsigned int log2_pb_w, unsigned int log2_pb_h,
                uint8_t inter_dir)
{
    /*FIXME Use specific DBF update function if DBF is disabled */
    struct DBFInfo *const dbf_info = &lc_ctx->dbf_info;
    uint8_t y_pu = y0 >> part_ctx->log2_min_cb_s;
    uint8_t x_pu = x0 >> part_ctx->log2_min_cb_s;
    uint8_t nb_pb_w = (1 << log2_pb_w) >> part_ctx->log2_min_cb_s;
    uint8_t nb_pb_h = (1 << log2_pb_h) >> part_ctx->log2_min_cb_s;

    if (inter_dir == 3) {
        struct VVCMVCtx *const mv_ctx0 = &inter_ctx->mv_ctx0;
        struct VVCMVCtx *const mv_ctx1 = &inter_ctx->mv_ctx1;

        fill_mvp_map(mv_ctx0, mv0, x_pu, y_pu, nb_pb_w, nb_pb_h);

        fill_mvp_map(mv_ctx1, mv1, x_pu, y_pu, nb_pb_w, nb_pb_h);

        fill_dbf_mv_map(dbf_info, mv_ctx0, mv0, x_pu, y_pu, nb_pb_w, nb_pb_h);

        fill_dbf_mv_map(dbf_info, mv_ctx1, mv1, x_pu, y_pu, nb_pb_w, nb_pb_h);

    } else if (inter_dir & 0x2) {
        struct VVCMVCtx *const mv_ctx0 = &inter_ctx->mv_ctx0;
        struct VVCMVCtx *const mv_ctx1 = &inter_ctx->mv_ctx1;

        fill_mvp_map(mv_ctx1, mv1, x_pu, y_pu, nb_pb_w, nb_pb_h);

        fill_dbf_mv_map_b(dbf_info, mv_ctx1, mv_ctx0, mv1, x_pu, y_pu, nb_pb_w, nb_pb_h);

    } else if (inter_dir & 0x1) {
        struct VVCMVCtx *const mv_ctx0 = &inter_ctx->mv_ctx0;
        struct VVCMVCtx *const mv_ctx1 = &inter_ctx->mv_ctx1;

        fill_mvp_map(mv_ctx0, mv0, x_pu, y_pu, nb_pb_w, nb_pb_h);

        fill_dbf_mv_map_b(dbf_info, mv_ctx0, mv_ctx1, mv0, x_pu, y_pu, nb_pb_w, nb_pb_h);

    }

    update_hmvp_lut_b(&inter_ctx->hmvp_lut, mv0, mv1, inter_dir);
}
