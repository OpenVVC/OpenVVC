#include <string.h>
#include "slicedec.h"
#include "ctudec.h"
#include "nvcl_structures.h"
#include "dec_structures.h"
#include "drv_lines.h"
#include "ovutils.h"
#include "ovmem.h"
#include "overror.h"
#include "ovdec_internal.h"

/* TODO define in a header */
enum SliceType {
     SLICE_B = 0,
     SLICE_P = 1,
     SLICE_I = 2
};

static void
free_inter_drv_lines(struct DRVLines *const drv_lns)
{
    struct InterLines *const lns = &drv_lns->inter_lines;

    if (lns->mv0) {
        ov_freep(&lns->mv0);
    }

    if (lns->mv1) {
        ov_freep(&lns->mv1);
    }

    if (lns->dir0) {
        ov_freep(&lns->dir0);
    }

    if (lns->dir1) {
        ov_freep(&lns->dir1);
    }

    if (lns->affine) {
        ov_freep(&lns->affine);
    }

    if (lns->aff_info) {
        ov_freep(&lns->aff_info);
    }
}

static int
init_inter_drv_lines(struct DRVLines *const drv_lns, int nb_pb_ctb,
                     int nb_ctb_pic_w)
{
    struct InterLines *const lns = &drv_lns->inter_lines;

    lns->mv0  = ov_mallocz(sizeof(*lns->mv0) * nb_pb_ctb * (nb_ctb_pic_w + 2));
    lns->mv1  = ov_mallocz(sizeof(*lns->mv1) * nb_pb_ctb * (nb_ctb_pic_w + 2));

    lns->dir0  = ov_mallocz(sizeof(*lns->dir0) * (nb_ctb_pic_w + 2));
    lns->dir1  = ov_mallocz(sizeof(*lns->dir1) * (nb_ctb_pic_w + 2));

    lns->affine = ov_mallocz(sizeof(*lns->affine) * (nb_ctb_pic_w + 2));
    lns->aff_info = ov_mallocz(sizeof(*lns->aff_info) * (nb_ctb_pic_w + 2));

    if (!lns->mv0 || !lns->mv1 || !lns->dir0 || !lns->dir1) {
        free_inter_drv_lines(drv_lns);
        return OVVC_ENOMEM;
    }

    return 0;
}

static void
offset_inter_drv_lines(struct DRVLines *const drv_lns, int ctb_offset,
                       uint8_t log2_ctb_s, uint8_t log2_min_cb_s)
{
    struct InterLines *const lns = &drv_lns->inter_lines;

    lns->mv0 += ctb_offset;
    lns->mv1 += ctb_offset;

    lns->dir0 += ctb_offset;
    lns->dir1 += ctb_offset;

    lns->affine   += ctb_offset;
    lns->aff_info += ctb_offset;
}

#if 0
/* Reset inter direction maps according to CTU neighbourhood
 */
static void
fill_inter_map(OVCTUDec *const ctudec, uint64_t above_map, uint64_t tr_map)
{
    const int nb_ctb_pb =  (1 << ((ctudec->part_ctx->log2_ctu_s) & 7)) >> ctudec->part_ctx->log2_min_cb_s;
    struct InterDRVCtx *const inter_ctx = &ctudec->drv_ctx.inter_ctx;

    uint64_t *const rows_map = inter_ctx->mv_ctx0.map.hfield;
    uint64_t *const cols_map = inter_ctx->mv_ctx0.map.vfield;

    const uint8_t ctb_ngh_flg = ctudec->ctu_ngh_flags;
    int i;

    rows_map[0] = 0;
    cols_map[0] = 0;

    if (ctb_ngh_flg & CTU_LFT_FLG) {
        cols_map[0] = cols_map[nb_ctb_pb];
        for (i = 1; i < nb_ctb_pb + 1; ++i) {
            uint8_t left_available = !!(cols_map[nb_ctb_pb] & (1llu << i));
            rows_map[i] = (uint64_t)left_available;
        }
    } else {
        for (int i = 1; i < nb_ctb_pb + 1; ++i) {
            rows_map[i] = 0;
        }
    }

    if (ctb_ngh_flg & CTU_LFT_FLG) {

        rows_map[0] = above_map << 1;

        if (ctb_ngh_flg & CTU_UPRGT_FLG) {
            rows_map[0] |= tr_map << (nb_ctb_pb + 1);
        }

        if (ctb_ngh_flg & CTU_UPLFT_FLG) {
            rows_map[0] |= cols_map[0] & 0x1;
        }

        for (i = 1; i < nb_ctb_pb + 1; i++) {
            uint8_t top_available = !!(above_map & (1llu << (i - 1)));
            cols_map[i] = (uint64_t)top_available;
        }
    } else {
        for (i = 1; i < nb_ctb_pb + 1; i++) {
            cols_map[i] = 0;
        }
    }
}
#endif

static void
tmvp_store_mv(OVCTUDec *ctudec)
{
    uint16_t ctb_x = ctudec->ctb_x;
    uint16_t ctb_y = ctudec->ctb_y;

    uint8_t log2_ctb_s = ctudec->part_ctx->log2_ctu_s;
    uint8_t log2_min_cb_s = ctudec->part_ctx->log2_min_cb_s;

    int nb_pb_ctb_w = (1 << log2_ctb_s) >> log2_min_cb_s;
    struct InterDRVCtx *const inter_ctx = &ctudec->drv_ctx.inter_ctx;

    const struct MVPlane *plane0 = inter_ctx->tmvp_ctx.plane0;
    const struct MVPlane *plane1 = inter_ctx->tmvp_ctx.plane1;

    int nb_ctb_w = ctudec->nb_ctb_pic_w;
    uint16_t ctb_addr_rs = ctb_x + ctb_y * nb_ctb_w;

    if (plane0->dirs) {
        uint64_t *dst_dirs = plane0->dirs + ctb_addr_rs * nb_pb_ctb_w;
        int32_t nb_tmvp_unit = nb_pb_ctb_w >> 1;

        int32_t pln_stride = nb_tmvp_unit * nb_ctb_w;
        int32_t ctb_offset = ctb_x * nb_tmvp_unit + (ctb_y * nb_tmvp_unit * pln_stride);
        OVMV *dst_mv = plane0->mvs + ctb_offset;
        struct OVMVCtx *mv_ctx = &inter_ctx->mv_ctx0;
        const OVMV *src2 = mv_ctx->mvs + 35;
        const OVMV *src = inter_ctx->tmvp_mv[0].mvs;
        int i, j;

        memcpy(dst_dirs, &mv_ctx->map.vfield[1], sizeof(uint64_t) * nb_pb_ctb_w);
        for (i = 0; i < nb_pb_ctb_w; i += 2) {
            for (j = 0; j < nb_pb_ctb_w; j += 2) {
                dst_mv[j >> 1] = src[j >> 1];
            }
            src2 += 34 * 2;
            src += 16;
            dst_mv += pln_stride;
        }
    }

    if (plane1->dirs) {
        struct OVMVCtx *mv_ctx = &inter_ctx->mv_ctx1;
        //const OVMV *src = mv_ctx->mvs + 35;
        const OVMV *src2 = mv_ctx->mvs + 35;
        const OVMV *src = inter_ctx->tmvp_mv[1].mvs;
        uint64_t *dst_dirs = plane1->dirs + ctb_addr_rs * nb_pb_ctb_w;
        int i, j;
        int32_t nb_tmvp_unit = nb_pb_ctb_w >> 1;
        int32_t pln_stride = nb_tmvp_unit * nb_ctb_w;
        int32_t ctb_offset = ctb_x * nb_tmvp_unit + (ctb_y * nb_tmvp_unit * pln_stride);
        OVMV *dst_mv = plane1->mvs + ctb_offset;

        /*FIXME memory could be spared with smaller map size when possible */
        memcpy(dst_dirs, &mv_ctx->map.vfield[1], sizeof(uint64_t) * nb_pb_ctb_w);
        for (i = 0; i < nb_pb_ctb_w; i += 2) {
            for (j = 0; j < nb_pb_ctb_w; j += 2) {
                dst_mv[j >> 1] = src[j >> 1];
            }
            src2 += 34 * 2;
            src += 16;
            dst_mv += pln_stride;
        }
    }
}

static void
rotate_affine_cp(struct AffineInfo *const aff_info, struct AffineInfo *const lns, uint64_t msk,
                 OVMV *mv0, OVMV *mv1, uint16_t nb_pb_ctb)
{
    //memcpy(&lns->aff_info[ctb_x * nb_ctb_pb], &aff_info[1 + nb_ctb_pb * 34], sizeof(struct AffineInfo) * nb_ctb_pb);
    int i = 0;
    while (msk) {
        if (msk & 0x1) {
            int x_pb = aff_info[i].pb.x_pb;
            int nb_pb_w = aff_info[i].pb.nb_pb_w;
            lns[i].cps[0].lt = mv0[x_pb + (nb_pb_ctb * 34)+ 1];
            lns[i].cps[0].rt = mv0[x_pb + (nb_pb_ctb * 34)+ nb_pb_w];
            lns[i].cps[1].lt = mv1[x_pb + (nb_pb_ctb * 34)+ 1];
            lns[i].cps[1].rt = mv1[x_pb + (nb_pb_ctb * 34)+ nb_pb_w];
            lns[i].type = aff_info[i].type;
            lns[i].pb = aff_info[i].pb;
        }

        i++;
        msk >>= 1;
    }
}

/* Copy last Motion Vector from CTU to corresponding line in 
 * DRVLine
 */
void
store_inter_maps(const struct DRVLines *const l,
                 OVCTUDec *const ctudec,
                 unsigned int ctb_x, uint8_t is_last)
{
    uint8_t nb_ctb_pb =  (1 << ((ctudec->part_ctx->log2_ctu_s) & 7)) >> ctudec->part_ctx->log2_min_cb_s;
    struct InterDRVCtx *const inter_ctx = &ctudec->drv_ctx.inter_ctx;
    struct AffineDRVInfo *const aff_ctx = &inter_ctx->affine_ctx;
    const struct InterLines  *const lns = &l->inter_lines;

    struct OVMVCtx *const mv_ctx0 = &inter_ctx->mv_ctx0;
    struct OVMVCtx *const mv_ctx1 = &inter_ctx->mv_ctx1;

    uint64_t *const rows_map0 = mv_ctx0->map.hfield;
    uint64_t *const cols_map0 = mv_ctx0->map.vfield;

    uint64_t *const rows_map1 = mv_ctx1->map.hfield;
    uint64_t *const cols_map1 = mv_ctx1->map.vfield;

    uint64_t *const rows_affn = aff_ctx->map.hfield;
    uint64_t *const cols_affn = aff_ctx->map.vfield;
    struct AffineInfo *const aff_info = aff_ctx->affine_info;

    int i;

    const uint64_t lst_row0 = rows_map0[nb_ctb_pb];
    const uint64_t lst_row1 = rows_map1[nb_ctb_pb];

    const uint64_t lst_row_aff = rows_affn[nb_ctb_pb];

    uint64_t above_map0 = (uint64_t)lns->dir0[ctb_x + 1] | ((uint64_t)lns->dir0 [ctb_x + 2] << nb_ctb_pb);
    uint64_t above_map1 = (uint64_t)lns->dir1[ctb_x + 1] | ((uint64_t)lns->dir1 [ctb_x + 2] << nb_ctb_pb);
    uint64_t affine_map = (uint64_t)lns->affine[ctb_x + 1] | ((uint64_t)lns->affine [ctb_x + 2] << nb_ctb_pb);

    if (is_last) {
       above_map0 = 0;
       above_map1 = 0;
       affine_map = 0;
    }

    tmvp_store_mv(ctudec);

    /* Copy last MV column to next CTU left MV column
     */
    for (i = 0; i < nb_ctb_pb + 1; i++) {
        uint64_t left_available0 = !!(cols_map0[nb_ctb_pb] & (1llu << i));
        uint64_t left_available1 = !!(cols_map1[nb_ctb_pb] & (1llu << i));
        uint64_t lft_aff = !!(cols_affn[nb_ctb_pb] & (1llu << i));
        rows_map0[i] = left_available0;
        rows_map1[i] = left_available1;
        rows_affn[i] = lft_aff;
        mv_ctx0->mvs[i * 34] = mv_ctx0->mvs[i* 34 + nb_ctb_pb];
        mv_ctx1->mvs[i * 34] = mv_ctx1->mvs[i* 34 + nb_ctb_pb];
        aff_info[i * 34] = aff_info[i* 34 + nb_ctb_pb];
    }

    cols_map0[0] = cols_map0[nb_ctb_pb];
    cols_map1[0] = cols_map1[nb_ctb_pb];
    cols_affn[0] = cols_affn[nb_ctb_pb];

    /* Replace CTU above MV line by line MV at ctb_x + 1*/
    memcpy(&mv_ctx0->mvs[1], &lns->mv0[(ctb_x + 1) * nb_ctb_pb], sizeof(OVMV) * (nb_ctb_pb));
    memcpy(&mv_ctx1->mvs[1], &lns->mv1[(ctb_x + 1) * nb_ctb_pb], sizeof(OVMV) * (nb_ctb_pb));
    memcpy(&aff_info[1], &lns->aff_info[(ctb_x + 1) * nb_ctb_pb], sizeof(struct AffineInfo) * (nb_ctb_pb));

    #if 1
    mv_ctx0->mvs[1 + nb_ctb_pb] = lns->mv0[(ctb_x + 2) * nb_ctb_pb];
    mv_ctx1->mvs[1 + nb_ctb_pb] = lns->mv1[(ctb_x + 2) * nb_ctb_pb];
    aff_info[1 + nb_ctb_pb] = lns->aff_info[(ctb_x + 2) * nb_ctb_pb];
    #endif

    for (i = 1; i < nb_ctb_pb + 1; i++) {
        uint64_t top_available0 = !!(above_map0 & (1llu << (i - 1)));
        uint64_t top_available1 = !!(above_map1 & (1llu << (i - 1)));
        uint64_t abv_affn = !!(affine_map & (1llu << (i - 1)));
        cols_map0[i] = (uint64_t)top_available0;
        cols_map1[i] = (uint64_t)top_available1;
        cols_affn[i] = abv_affn;
    }

    rows_map0[0] |= above_map0 << 1;
    rows_map1[0] |= above_map1 << 1;
    rows_affn[0] |= affine_map << 1;

    /* Save last CTU MV line to line at ctb_x */
    memcpy(&lns->mv0[ctb_x * nb_ctb_pb], &mv_ctx0->mvs[1 + nb_ctb_pb * 34], sizeof(OVMV) * nb_ctb_pb);
    memcpy(&lns->mv1[ctb_x * nb_ctb_pb], &mv_ctx1->mvs[1 + nb_ctb_pb * 34], sizeof(OVMV) * nb_ctb_pb);
    rotate_affine_cp(&aff_info[1 + nb_ctb_pb * 34], &lns->aff_info[ctb_x * nb_ctb_pb],
                     lst_row_aff >> 1, mv_ctx0->mvs, mv_ctx1->mvs, nb_ctb_pb);

    /* Store last inter dir info onto line */
    lns->dir0[ctb_x] = (uint32_t)(lst_row0 >> 1);
    lns->dir1[ctb_x] = (uint32_t)(lst_row1 >> 1);
    lns->affine[ctb_x] = (uint32_t)(lst_row_aff >> 1);
}

static void
free_dbf_lines(struct DBFLines *const l)
{
    /* FIXME check before free */
    ov_freep(&l->qp_x_map);
    ov_freep(&l->qp_x_map_cb);
    ov_freep(&l->qp_x_map_cr);

    ov_freep(&l->dbf_edge_hor);
    ov_freep(&l->dbf_edge_hor_c);

    ov_freep(&l->large_map_c);
    ov_freep(&l->small_map);

    ov_freep(&l->dbf_bs2_hor);
    ov_freep(&l->dbf_bs2_hor_c);

    ov_freep(&l->dbf_bs1_hor);
    ov_freep(&l->dbf_bs1_hor_cb);
    ov_freep(&l->dbf_bs1_hor_cr);
}


static int
init_dbf_lines(struct DBFLines *const l, int nb_ctu_line, int nb_pu_line)
{
    uint8_t malloc_chk = 0;

    l->qp_x_map       = ov_mallocz((nb_pu_line + 1) * sizeof(int8_t));
    l->qp_x_map_cb    = ov_mallocz((nb_pu_line + 1) * sizeof(int8_t));
    l->qp_x_map_cr    = ov_mallocz((nb_pu_line + 1) * sizeof(int8_t));

    l->dbf_edge_hor   = ov_mallocz((nb_ctu_line + 1) * sizeof(uint64_t));
    l->dbf_edge_hor_c = ov_mallocz((nb_ctu_line + 1) * sizeof(uint64_t));

    l->small_map      = ov_mallocz((nb_ctu_line + 1) * sizeof(uint64_t));
    l->large_map_c    = ov_mallocz((nb_ctu_line + 1) * sizeof(uint64_t));

    l->dbf_bs1_hor    = ov_mallocz((nb_ctu_line + 1) * sizeof(uint64_t));
    l->dbf_bs1_hor_cb = ov_mallocz((nb_ctu_line + 1) * sizeof(uint64_t));
    l->dbf_bs1_hor_cr = ov_mallocz((nb_ctu_line + 1) * sizeof(uint64_t));

    l->dbf_bs2_hor    = ov_mallocz((nb_ctu_line + 1) * sizeof(uint64_t));
    l->dbf_bs2_hor_c  = ov_mallocz((nb_ctu_line + 1) * sizeof(uint64_t));

    malloc_chk |= l->qp_x_map       == NULL;
    malloc_chk |= l->qp_x_map_cb    == NULL;
    malloc_chk |= l->qp_x_map_cr    == NULL;

    malloc_chk |= l->dbf_edge_hor   == NULL;
    malloc_chk |= l->dbf_edge_hor_c == NULL;

    malloc_chk |= l->large_map_c    == NULL;

    malloc_chk |= l->dbf_bs1_hor    == NULL;
    malloc_chk |= l->dbf_bs1_hor_cb == NULL;
    malloc_chk |= l->dbf_bs1_hor_cr == NULL;

    malloc_chk |= l->dbf_bs2_hor    == NULL;
    malloc_chk |= l->dbf_bs2_hor_c  == NULL;

    if (malloc_chk) {
        free_dbf_lines(l);
        return OVVC_ENOMEM;
    }

    return 0;
}

static void
offset_dbf_lines(struct DBFLines *const l, int ctb_offset,
                 uint8_t log2_ctb_s, uint8_t log2_min_cb_s)
{
    int pb_offset = (ctb_offset << log2_ctb_s) >> log2_min_cb_s;
    pb_offset <<= 5;

    l->qp_x_map       += pb_offset;
    l->qp_x_map_cb    += pb_offset;
    l->qp_x_map_cr    += pb_offset;

    l->dbf_edge_hor   += ctb_offset;
    l->dbf_edge_hor_c += ctb_offset;

    l->large_map_c    += ctb_offset;
    l->small_map      += ctb_offset;

    l->dbf_bs1_hor    += ctb_offset;
    l->dbf_bs1_hor_cb += ctb_offset;
    l->dbf_bs1_hor_cr += ctb_offset;

    l->dbf_bs2_hor    += ctb_offset;
    l->dbf_bs2_hor_c  += ctb_offset;

}

static void
dbf_clear_lines(const struct DBFLines *const l, int nb_ctu_line, int nb_pu_line)
{
    memset(l->qp_x_map,    0, (nb_pu_line + 1) * sizeof(int8_t));
    memset(l->qp_x_map_cb, 0, (nb_pu_line + 1) * sizeof(int8_t));
    memset(l->qp_x_map_cr, 0, (nb_pu_line + 1) * sizeof(int8_t));

    memset(l->dbf_edge_hor,   0, (nb_ctu_line + 1) * sizeof(uint64_t));
    memset(l->dbf_edge_hor_c, 0, (nb_ctu_line + 1) * sizeof(uint64_t));

    memset(l->small_map,   0, (nb_ctu_line + 1) * sizeof(uint64_t));
    memset(l->large_map_c, 0, (nb_ctu_line + 1) * sizeof(uint64_t));

    memset(l->dbf_bs1_hor,    0, (nb_ctu_line + 1) * sizeof(uint64_t));
    memset(l->dbf_bs1_hor_cb, 0, (nb_ctu_line + 1) * sizeof(uint64_t));
    memset(l->dbf_bs1_hor_cr, 0, (nb_ctu_line + 1) * sizeof(uint64_t));

    memset(l->dbf_bs2_hor,   0, (nb_ctu_line + 1) * sizeof(uint64_t));
    memset(l->dbf_bs2_hor_c, 0, (nb_ctu_line + 1) * sizeof(uint64_t));
}

static void
dbf_load_qp_map(struct DBFInfo *const dbf_info, const struct DBFLines *const l,
                struct DRVLines *l2,
                uint8_t log2_ctu_s, int ctb_x)
{
    uint8_t nb_units_ctb = 1 << (log2_ctu_s & 7) >> 2;
    int i;

    for (i = 0; i < nb_units_ctb + 1; ++i) {
        dbf_info->qp_map_y.hor [i * 34] = dbf_info->qp_map_y.hor [i * 34 + nb_units_ctb];
        dbf_info->qp_map_cb.hor[i * 34] = dbf_info->qp_map_cb.hor[i * 34 + nb_units_ctb];
        dbf_info->qp_map_cr.hor[i * 34] = dbf_info->qp_map_cr.hor[i * 34 + nb_units_ctb];

        dbf_info->qp_map_y.hor [1 + i * 34] = dbf_info->qp_map_y.hor [i * 34 + nb_units_ctb + 1];
        dbf_info->qp_map_cb.hor[1 + i * 34] = dbf_info->qp_map_cb.hor[i * 34 + nb_units_ctb + 1];
        dbf_info->qp_map_cr.hor[1 + i * 34] = dbf_info->qp_map_cr.hor[i * 34 + nb_units_ctb + 1];

        dbf_info->qp_map_y.ver [i] = dbf_info->qp_map_y.ver [i + 34 * nb_units_ctb];
        dbf_info->qp_map_cb.ver[i] = dbf_info->qp_map_cb.ver[i + 34 * nb_units_ctb];
        dbf_info->qp_map_cr.ver[i] = dbf_info->qp_map_cr.ver[i + 34 * nb_units_ctb];
    }

    for (i = 0; i < nb_units_ctb; ++i) {
        dbf_info->qp_map_y.hor [2 + i] = l->qp_x_map   [(ctb_x << 5) + i];
        dbf_info->qp_map_cb.hor[2 + i] = l->qp_x_map_cb[(ctb_x << 5) + i];
        dbf_info->qp_map_cr.hor[2 + i] = l->qp_x_map_cr[(ctb_x << 5) + i];
    }
}

static void
dbf_store_qp_map(const struct DBFInfo *const dbf_info,
                 const struct DBFLines *const l,
                 struct DRVLines *l2,
                 uint8_t log2_ctu_s,
                 unsigned int ctb_x)
{
    uint8_t nb_units_ctb =  1 << (log2_ctu_s & 7) >> 2;

    memcpy(&l->qp_x_map   [ctb_x << 5], &dbf_info->qp_map_y.hor[2 + 34*nb_units_ctb],  sizeof(uint8_t) * nb_units_ctb);
    memcpy(&l->qp_x_map_cb[ctb_x << 5], &dbf_info->qp_map_cb.hor[2 + 34*nb_units_ctb], sizeof(uint8_t) * nb_units_ctb);
    memcpy(&l->qp_x_map_cr[ctb_x << 5], &dbf_info->qp_map_cr.hor[2 + 34*nb_units_ctb], sizeof(uint8_t) * nb_units_ctb);
}

static void
dbf_load_edge_map(struct DBFInfo *const dbf_info, const struct DBFLines *const l,
                  uint8_t log2_ctu_s, int ctb_x)
{
    int nb_pb_s = (1 << log2_ctu_s) >> 2;
    uint64_t ctb_lft_msk = (uint64_t)-(!!ctb_x);

    if (ctb_x) {
        /* copy previous 8 vertical edges used for filter length derivation */
        uint64_t *ctb_bnd = dbf_info->ctb_bound_ver;
        uint64_t *ctb_bnd_c = dbf_info->ctb_bound_ver_c;
        memcpy(ctb_bnd, &ctb_bnd[nb_pb_s], 8 * sizeof(uint64_t));
        memcpy(ctb_bnd_c + 8 - 3, &ctb_bnd_c[nb_pb_s + 8 - 3], 3 * sizeof(uint64_t));
    }

    for (int i = 0; i < nb_pb_s; ++i) {
        dbf_info->edge_map_hor[i] = ctb_lft_msk & (dbf_info->edge_map_hor[i] >> (nb_pb_s)) & 0x3;
        dbf_info->ctb_bound_hor[i + 8] = ctb_lft_msk & (dbf_info->ctb_bound_hor[i + 8] >> (nb_pb_s)) & 0x3;
        dbf_info->edge_map_ver[i] = 0;
        dbf_info->ctb_bound_ver[i + 8] = 0;
    }

    for (int i = 0; i < nb_pb_s; ++i) {
        dbf_info->edge_map_hor_c[i] = ctb_lft_msk & (dbf_info->edge_map_hor_c[i] >> (nb_pb_s)) & 0x3;
        dbf_info->ctb_bound_hor_c[i + 8] = ctb_lft_msk & (dbf_info->ctb_bound_hor_c[i + 8] >> (nb_pb_s)) & 0x3;
        dbf_info->edge_map_ver_c[i] = 0;
        dbf_info->ctb_bound_ver_c[i + 8] = 0;
    }

    dbf_info->ctb_bound_hor[7] = l->small_map[ctb_x];
    /* FIXME move at decoder intialisation */
    dbf_info->ctb_bound_hor[6] = (uint64_t) - 1ll;
    /* ctb_bound_hor requires its two first bits set to one since we need 2 last units
     * from left CTU.
     */
    dbf_info->ctb_bound_hor[8]           = 0x3;
    dbf_info->ctb_bound_hor[8 + nb_pb_s] = 0x3;

    dbf_info->edge_map_hor  [0] |= l->dbf_edge_hor[ctb_x];
    dbf_info->edge_map_hor_c[0] |= l->dbf_edge_hor_c[ctb_x];
    /* FIXME try to remove this */
    dbf_info->edge_map_ver[0]    = ctb_lft_msk & dbf_info->edge_map_ver[nb_pb_s];
    dbf_info->edge_map_ver_c[0]  = ctb_lft_msk & dbf_info->edge_map_ver_c[nb_pb_s];
}

static void
dbf_store_edge_map(struct DBFInfo *const dbf_info, const struct DBFLines *const l,
                   uint8_t log2_ctu_s, int ctb_x)
{
    int nb_pb_s = (1 << log2_ctu_s) >> 2;

    /* Store last edge line info and reset first column */
    l->dbf_edge_hor[ctb_x]   = dbf_info->edge_map_hor[nb_pb_s];
    l->dbf_edge_hor_c[ctb_x] = dbf_info->edge_map_hor_c[nb_pb_s];

    /* Use last horizontal edges to determine if last CUs were large
     * if there were no edges large in the last three lines will be 0 
     */
    l->large_map_c[ctb_x]  = dbf_info->edge_map_hor_c[nb_pb_s - 1];
    l->large_map_c[ctb_x] |= dbf_info->edge_map_hor_c[nb_pb_s - 2];
    l->large_map_c[ctb_x] |= dbf_info->edge_map_hor_c[nb_pb_s - 3];

    l->small_map[ctb_x] = dbf_info->edge_map_hor[nb_pb_s - 1];
}

static void
dbf_load_bs_map(struct DBFInfo *const dbf_info, const struct DBFLines *const l,
            uint8_t log2_ctu_s, int ctb_x)
{
    int nb_pb_s = (1 << log2_ctu_s) >> 2;
    struct DBFMap *const bs2_map    = &dbf_info->bs2_map;
    struct DBFMap *const bs2_map_c  = &dbf_info->bs2_map_c;
    struct DBFMap *const bs1_map    = &dbf_info->bs1_map;
    struct DBFMap *const bs1_map_cb = &dbf_info->bs1_map_cb;
    struct DBFMap *const bs1_map_cr = &dbf_info->bs1_map_cr;
    int i;

    /* FIXME check cast */
    uint64_t ctb_lft_msk = -(!!ctb_x);

    dbf_info->ctb_bound_hor_c[8 - 1] = l->large_map_c[ctb_x];

    bs2_map->ver[0]   = ctb_lft_msk & bs2_map->ver[nb_pb_s];
    bs2_map_c->ver[0] = ctb_lft_msk & bs2_map_c->ver[nb_pb_s];

    bs1_map->ver[0]    = ctb_lft_msk & bs1_map->ver[nb_pb_s];
    bs1_map_cb->ver[0] = ctb_lft_msk & bs1_map_cb->ver[nb_pb_s];
    bs1_map_cr->ver[0] = ctb_lft_msk & bs1_map_cr->ver[nb_pb_s];

    bs2_map->hor[0]   = ctb_lft_msk & (bs2_map->hor[0]   >> nb_pb_s) & 0x3;
    bs2_map_c->hor[0] = ctb_lft_msk & (bs2_map_c->hor[0] >> nb_pb_s) & 0x3;

    bs1_map->hor[0]    = ctb_lft_msk & (bs1_map->hor[0]    >> nb_pb_s) & 0x3;
    bs1_map_cb->hor[0] = ctb_lft_msk & (bs1_map_cb->hor[0] >> nb_pb_s) & 0x3;
    bs1_map_cr->hor[0] = ctb_lft_msk & (bs1_map_cr->hor[0] >> nb_pb_s) & 0x3;

    bs2_map->hor[0]   |= l->dbf_bs2_hor[ctb_x];
    bs2_map_c->hor[0] |= l->dbf_bs2_hor_c[ctb_x];

    bs1_map->hor[0]    |= l->dbf_bs1_hor[ctb_x];
    bs1_map_cb->hor[0] |= l->dbf_bs1_hor_cb[ctb_x];
    bs1_map_cr->hor[0] |= l->dbf_bs1_hor_cr[ctb_x];

    /*FIXME separate hor and ver loop  and merge all loops*/

    for (i = 1; i < nb_pb_s + 1; ++i) {
        bs1_map->hor[i] = ctb_lft_msk & (bs1_map->hor[i] >> nb_pb_s) & 0x3;
        bs1_map->ver[i] = 0;
    }

    for (i = 1; i < nb_pb_s + 1; ++i) {
        bs2_map->hor[i] = ctb_lft_msk & (bs2_map->hor[i] >> nb_pb_s) & 0x3;
        bs2_map->ver[i] = 0;
    }

    for (i = 1; i < nb_pb_s + 1; ++i) {
        bs2_map_c->hor[i] = ctb_lft_msk & (bs2_map_c->hor[i] >> nb_pb_s) & 0x3;
        bs2_map_c->ver[i] = 0;
    }

    for (i = 1; i < nb_pb_s + 1; ++i) {
        bs1_map_cb->hor[i] = ctb_lft_msk & (bs1_map_cb->hor[i] >> nb_pb_s) & 0x3;
        bs1_map_cb->ver[i] = 0;
    }

    for (i = 1; i < nb_pb_s + 1; ++i) {
        bs1_map_cr->hor[i] = ctb_lft_msk & (bs1_map_cr->hor[i] >> nb_pb_s) & 0x3;
        bs1_map_cr->ver[i] = 0;
    }
}

static void
dbf_store_bs_map(struct DBFInfo *const dbf_info, const struct DBFLines *const l,
              uint8_t log2_ctu_s, int ctb_x)
{
    int nb_pb_s = (1 << log2_ctu_s) >> 2;

    struct DBFMap *const bs2_map = &dbf_info->bs2_map;
    struct DBFMap *const bs2_map_c = &dbf_info->bs2_map_c;

    struct DBFMap *const bs1_map = &dbf_info->bs1_map;
    struct DBFMap *const bs1_map_cb = &dbf_info->bs1_map_cb;
    struct DBFMap *const bs1_map_cr = &dbf_info->bs1_map_cr;

    l->dbf_bs2_hor[ctb_x]   =  bs2_map->hor[nb_pb_s];
    l->dbf_bs2_hor_c[ctb_x] =  bs2_map_c->hor[nb_pb_s];

    l->dbf_bs1_hor[ctb_x]    =  bs1_map->hor[nb_pb_s];
    l->dbf_bs1_hor_cb[ctb_x] =  bs1_map_cb->hor[nb_pb_s];
    l->dbf_bs1_hor_cr[ctb_x] =  bs1_map_cr->hor[nb_pb_s];
}

void
dbf_load_info(struct DBFInfo *const dbf_info,
              const struct DBFLines *const dbf_lines,
              uint8_t log2_ctu_s, int ctb_x)
{
    dbf_load_edge_map(dbf_info, dbf_lines, log2_ctu_s, ctb_x);
    dbf_load_bs_map(dbf_info, dbf_lines, log2_ctu_s, ctb_x);
    dbf_load_qp_map(dbf_info, dbf_lines, NULL, log2_ctu_s, ctb_x);
}

void
dbf_store_info(struct DBFInfo *const dbf_info,
               const struct DBFLines *const dbf_lines,
               uint8_t log2_ctu_s, int ctb_x)
{
    dbf_store_edge_map(dbf_info, dbf_lines, log2_ctu_s, ctb_x);
    dbf_store_bs_map(dbf_info, dbf_lines, log2_ctu_s, ctb_x);
    dbf_store_qp_map(dbf_info, dbf_lines, NULL, log2_ctu_s, ctb_x);
}

int
init_drv_lines(OVSliceDec *sldec, const OVPS *const prms)
{
     int ret;
     uint8_t slice_type = sldec->slice_type;
    const OVPartInfo *const pinfo = slice_type == SLICE_I ? &prms->sps_info.part_info[0]
                                                          : &prms->sps_info.part_info[1];
     const OVSPS *const sps = prms->sps;
     const struct TileInfo *tinfo = &prms->pps_info.tile_info;

     struct DRVLines *const lns = &sldec->drv_lines;

     uint8_t log2_ctb_s    = pinfo->log2_ctu_s;
     uint8_t log2_min_cb_s = pinfo->log2_min_cb_s;

     /* TODO use active parameters such as generic pic info
      * or something instead of this since this could be
      * overridden by PPS in case of sub_pic etc.
      * we could compute those values once for all earlier
      * in the decoding process
      */

     uint16_t pic_w = sps->sps_pic_width_max_in_luma_samples;

     uint16_t nb_ctb_pic_w = (pic_w + ((1 << log2_ctb_s) - 1)) >> log2_ctb_s;

     nb_ctb_pic_w += tinfo->nb_tile_cols * 2;
     nb_ctb_pic_w *= tinfo->nb_tile_rows;

     uint16_t nb_pb_pic_w = nb_ctb_pic_w << (log2_ctb_s - log2_min_cb_s);
     uint8_t nb_pb_ctb = (1 << log2_ctb_s) >> log2_min_cb_s;

     ret = init_inter_drv_lines(lns, nb_pb_ctb, 32*nb_ctb_pic_w);

     /* FIXME return */
     ret = init_dbf_lines(&lns->dbf_lines, nb_ctb_pic_w, 32 * nb_pb_pic_w);

     lns->intra_luma_x  = ov_mallocz(sizeof(*lns->intra_luma_x) * nb_pb_pic_w);

     if (!lns->intra_luma_x || ret < 0) {
         drv_lines_uninit(sldec);
         return OVVC_ENOMEM;
     }

     return 0;
}

void
reset_drv_lines(OVSliceDec *sldec, const OVPS *const prms)
{
     uint8_t slice_type = sldec->slice_type;
    const OVPartInfo *const pinfo = slice_type == SLICE_I ? &prms->sps_info.part_info[0]
                                                          : &prms->sps_info.part_info[1];
    const OVSPS *const sps = prms->sps;
     const struct TileInfo *tinfo = &prms->pps_info.tile_info;

    struct DRVLines *const lns = &sldec->drv_lines;

    uint8_t log2_ctb_s = pinfo->log2_ctu_s;
    uint8_t log2_min_cb_s = pinfo->log2_min_cb_s;
    /* TODO use active parameters such as generic pic info
     * see init_cabac_lines
     */
    uint16_t pic_w = sps->sps_pic_width_max_in_luma_samples;
    uint16_t nb_ctb_pic_w = (pic_w + ((1 << log2_ctb_s) - 1)) >> log2_ctb_s;

     nb_ctb_pic_w += tinfo->nb_tile_cols;
     nb_ctb_pic_w *= tinfo->nb_tile_rows;

    uint16_t nb_pb_pic_w = nb_ctb_pic_w << (log2_ctb_s - log2_min_cb_s);

    struct InterLines *const i_lns = &lns->inter_lines;

    /*FIXME check 32 */
    memset(i_lns->dir0, 0, (sizeof(*i_lns->dir0) * 32*(nb_ctb_pic_w + 1)));
    memset(i_lns->dir1, 0, (sizeof(*i_lns->dir1) * 32*(nb_ctb_pic_w + 1)));
    memset(i_lns->affine, 0, (sizeof(*i_lns->affine) * 32*(nb_ctb_pic_w + 1)));

    /* PLANAR  = 0 value is used if absent so we use it as reset value
    */
    memset(lns->intra_luma_x,     0,  sizeof(*lns->intra_luma_x) * nb_pb_pic_w);

    dbf_clear_lines(&lns->dbf_lines, nb_ctb_pic_w,  32*nb_pb_pic_w);
}

static void
load_first_ctu_inter(const struct DRVLines *const l,
                 OVCTUDec *const ctudec,
                 unsigned int ctb_x)
{
    uint8_t nb_ctb_pb =  (1 << ((ctudec->part_ctx->log2_ctu_s) & 7)) >> ctudec->part_ctx->log2_min_cb_s;

    struct InterDRVCtx *const inter_ctx = &ctudec->drv_ctx.inter_ctx;
    struct AffineDRVInfo *const aff_ctx = &inter_ctx->affine_ctx;
    const struct InterLines *const lns = &l->inter_lines;

    struct OVMVCtx *const mv_ctx0 = &inter_ctx->mv_ctx0;
    struct OVMVCtx *const mv_ctx1 = &inter_ctx->mv_ctx1;
    struct AffineInfo *const aff_info = aff_ctx->affine_info;

    uint64_t *const rows_map0 = mv_ctx0->map.hfield;
    uint64_t *const cols_map0 = mv_ctx0->map.vfield;

    uint64_t *const rows_map1 = mv_ctx1->map.hfield;
    uint64_t *const cols_map1 = mv_ctx1->map.vfield;

    uint64_t *const rows_affn = aff_ctx->map.hfield;
    uint64_t *const cols_affn = aff_ctx->map.vfield;

    int i;

    uint64_t above_map0 = (uint64_t)lns->dir0[0] | ((uint64_t)lns->dir0[1] << nb_ctb_pb);
    uint64_t above_map1 = (uint64_t)lns->dir1[0] | ((uint64_t)lns->dir1[1] << nb_ctb_pb);
    uint64_t affine_map = (uint64_t)lns->affine[0] | ((uint64_t)lns->affine[1] << nb_ctb_pb);

    for (i = 0; i < nb_ctb_pb + 1; i++) {
        #if 0
        uint64_t left_available0 = !!(cols_map0[nb_ctb_pb] & (1llu << i));
        uint64_t left_available1 = !!(cols_map1[nb_ctb_pb] & (1llu << i));
        #endif
        rows_map0[i] = 0;
        rows_map1[i] = 0;
        rows_affn[i] = 0;
        #if 0
        mv_ctx0->mvs[i * 34] = mv_ctx0->mvs[i* 34 + nb_ctb_pb];
        mv_ctx1->mvs[i * 34] = mv_ctx1->mvs[i* 34 + nb_ctb_pb];
        #endif
    }
    cols_map0[0] = 0;
    cols_map1[0] = 0;
    cols_affn[0] = 0;

    memcpy(&mv_ctx0->mvs[1], &lns->mv0[0], sizeof(OVMV) * nb_ctb_pb);
    memcpy(&mv_ctx1->mvs[1], &lns->mv1[0], sizeof(OVMV) * nb_ctb_pb);

    memcpy(&aff_info[1], &lns->aff_info[0], sizeof(struct AffineInfo) * nb_ctb_pb);

    mv_ctx0->mvs[1 + nb_ctb_pb] = lns->mv0[nb_ctb_pb];
    mv_ctx1->mvs[1 + nb_ctb_pb] = lns->mv1[nb_ctb_pb];
    aff_info[1 + nb_ctb_pb] = lns->aff_info[nb_ctb_pb];

    for (i = 1; i < nb_ctb_pb + 1; i++) {
        uint64_t top_available0 = !!(above_map0 & (1llu << (i - 1)));
        uint64_t top_available1 = !!(above_map1 & (1llu << (i - 1)));
        uint64_t abv_affine = !!(affine_map & (1llu << (i - 1)));
        cols_map0[i] = top_available0;
        cols_map1[i] = top_available1;
        cols_affn[i] = abv_affine;
    }
    rows_map0[0] |= above_map0 << 1;
    rows_map1[0] |= above_map1 << 1;
    rows_affn[0] |= affine_map << 1;

    /* Reset HMVP Look Up table */
    inter_ctx->hmvp_lut.nb_mv = 0;
}


void
drv_line_next_line(OVCTUDec *const ctudec, const struct DRVLines *const lns)
{
    #if 0
    struct OVDrvCtx *const drv_ctx        = &ctudec->drv_ctx;
    #endif
    struct IntraDRVInfo *const intra_info = &ctudec->drv_ctx.intra_info;
    const OVPartInfo *const pinfo = ctudec->part_ctx;

    uint8_t log2_ctb_s    = pinfo->log2_ctu_s;
    uint8_t log2_min_cb_s = pinfo->log2_min_cb_s;

    uint16_t nb_pb_ctb_w = (1 << log2_ctb_s) >> log2_min_cb_s;

    /* FIXME
     *     done twice on new entry see (reset lines function)
     *     use partition limitations for reset
     */
    /* Reset to 0 == PLANAR */
    struct OVDrvCtx *const drv_ctx = &ctudec->drv_ctx;
    int8_t qp_val = ctudec->qp_ctx.current_qp;

    /* FIXME remove if unused */
    intra_info->luma_mode_x = lns->intra_luma_x;

    load_first_ctu_inter(lns, ctudec, 0);

    memset(intra_info->luma_mode_y, 0, sizeof(*intra_info->luma_mode_y) * nb_pb_ctb_w);

    memset(drv_ctx->qp_map_x, qp_val, sizeof(*drv_ctx->qp_map_x) * nb_pb_ctb_w);
    memset(drv_ctx->qp_map_y, qp_val, sizeof(*drv_ctx->qp_map_y) * nb_pb_ctb_w);

    #if 0
    memset(&ctudec->dbf_info, 0, sizeof(ctudec->dbf_info));
    #endif
    dbf_load_info(&ctudec->dbf_info, &lns->dbf_lines, log2_ctb_s, 0);
}


#if 0
void
drv_line_next_ctu(OVCTUDec *const ctudec, OVSliceDec *sldec, struct DRVLines *drv_lns,
                  const OVPS *const prms, uint16_t ctb_x)
{
    const OVPartInfo *const pinfo = ctudec->part_ctx;
    #if 0
    struct IntraDRVInfo *const intra_info = &ctudec->drv_ctx.intra_info;
    #endif
    struct OVDrvCtx *const drv_ctx = &ctudec->drv_ctx;

    uint8_t log2_ctb_s    = pinfo->log2_ctu_s;
    uint8_t log2_min_cb_s = pinfo->log2_min_cb_s;

    uint16_t nb_pb_ctb_w = (1 << log2_ctb_s) >> log2_min_cb_s;

    int8_t qp_val = ctudec->qp_ctx.current_qp;

    /* FIXME Unecessary ? */
    #if 0
    memset(intra_info->luma_mode_x, 0, sizeof(*intra_info->luma_mode_x) * nb_pb_ctb_w);
    #endif

    /* Reset QP prediction map to previous current QP prediction
     * value
     */
    #if 1

    memset(drv_ctx->qp_map_x, qp_val, sizeof(*drv_ctx->qp_map_x) * nb_pb_ctb_w);
    memset(drv_ctx->qp_map_y, qp_val, sizeof(*drv_ctx->qp_map_y) * nb_pb_ctb_w);

    #endif
}
#endif

void
offset_drv_lines(struct DRVLines *const lns, uint8_t tile_x, uint8_t tile_y,
                 uint8_t ctb_x,
                 uint8_t log2_ctb_s, uint8_t log2_min_cb_s,
                 uint8_t  nb_tile_cols, uint16_t nb_ctb_pic_w)
{
     int ctb_offset = (uint32_t) (nb_ctb_pic_w + nb_tile_cols * 1) * tile_y;

     ctb_offset += tile_x;
     ctb_offset += ctb_x;

     offset_inter_drv_lines(lns, (ctb_offset << log2_ctb_s) >> log2_min_cb_s, log2_ctb_s, log2_min_cb_s);

     ctb_offset -= tile_x;
     ctb_offset -= tile_y * nb_tile_cols;
     /* FIXME return */
     offset_dbf_lines(&lns->dbf_lines, ctb_offset, log2_ctb_s, log2_min_cb_s);
}

void
drv_lines_uninit(OVSliceDec *sldec)
{
     struct DRVLines *const lns = &sldec->drv_lines;

     ov_freep(&lns->intra_luma_x);
     free_inter_drv_lines(lns);
     free_dbf_lines(&lns->dbf_lines);
}

