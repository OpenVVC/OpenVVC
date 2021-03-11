#include <string.h>
#include "slicedec.h"
#include "ctudec.h"
#include "nvcl_structures.h"
#include "dec_structures.h"
#include "drv_lines.h"
#include "ovutils.h"
#include "ovmem.h"
#include "overror.h"

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
}

static void
inter_lines_next_ctu()
{
}

static int
init_inter_drv_lines(struct DRVLines *const drv_lns, int nb_pb_pic_w,
                     int nb_ctb_pic_w)
{
    struct InterLines *const lns = &drv_lns->inter_lines;

    lns->mv0  = ov_mallocz(sizeof(*lns->mv0) * nb_ctb_pic_w * 32);
    lns->mv1  = ov_mallocz(sizeof(*lns->mv1) * nb_ctb_pic_w * 32);

    lns->dir0  = ov_mallocz(sizeof(*lns->dir0) * (nb_ctb_pic_w + 1));
    lns->dir1  = ov_mallocz(sizeof(*lns->dir1) * (nb_ctb_pic_w + 1));

    if (!lns->mv0 || !lns->mv1 || !lns->dir0 || !lns->dir1) {
        free_inter_drv_lines(drv_lns);
         return OVVC_ENOMEM;
    }

    return 0;
}

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

/* Copy last Motion Vector from CTU to corresponding line in 
 * DRVLine
 */
void
store_inter_maps(struct DRVLines *const l,
                 const OVCTUDec *const ctudec,
                 unsigned int ctb_x)
{
    uint8_t nb_ctb_pb =  (1 << ((ctudec->part_ctx->log2_ctu_s) & 7)) >> ctudec->part_ctx->log2_min_cb_s;
    struct InterDRVCtx *const inter_ctx = &ctudec->drv_ctx.inter_ctx;
    struct InterLines  *const lns = &l->inter_lines;

    struct OVMVCtx *const mv_ctx0 = &inter_ctx->mv_ctx0;
    struct OVMVCtx *const mv_ctx1 = &inter_ctx->mv_ctx1;

    uint64_t *const rows_map0 = mv_ctx0->map.hfield;
    uint64_t *const cols_map0 = mv_ctx0->map.vfield;

    uint64_t *const rows_map1 = mv_ctx1->map.hfield;
    uint64_t *const cols_map1 = mv_ctx1->map.vfield;

    int i;

    const uint64_t lst_row0 = rows_map0[nb_ctb_pb];
    #if 0
    const uint64_t lst_col0 = cols_map0[nb_ctb_pb];
    #endif

    const uint64_t lst_row1 = rows_map1[nb_ctb_pb];
    #if 0
    const uint64_t lst_col1 = cols_map1[nb_ctb_pb];
    #endif

    uint64_t above_map0 = (uint64_t)lns->dir0[ctb_x + 1] | ((uint64_t)lns->dir0 [ctb_x + 2] << nb_ctb_pb);
    uint64_t above_map1 = (uint64_t)lns->dir1[ctb_x + 1] | ((uint64_t)lns->dir1 [ctb_x + 2] << nb_ctb_pb);


    /* Copy last MV column to next CTU left MV column
     */
    for (i = 0; i < nb_ctb_pb + 1; i++) {
        uint64_t left_available0 = !!(cols_map0[nb_ctb_pb] & (1llu << i));
        uint64_t left_available1 = !!(cols_map1[nb_ctb_pb] & (1llu << i));
        rows_map0[i] = left_available0;
        rows_map1[i] = left_available1;
        mv_ctx0->mvs[i * 34] = mv_ctx0->mvs[i* 34 + nb_ctb_pb];
        mv_ctx1->mvs[i * 34] = mv_ctx1->mvs[i* 34 + nb_ctb_pb];
    }

    cols_map0[0] = cols_map0[nb_ctb_pb];
    cols_map1[0] = cols_map1[nb_ctb_pb];

    /* Replace CTU above MV line by line MV at ctb_x + 1*/
    memcpy(&mv_ctx0->mvs[1], &lns->mv0[(ctb_x + 1) << 5], sizeof(OVMV) * nb_ctb_pb);
    memcpy(&mv_ctx1->mvs[1], &lns->mv1[(ctb_x + 1) << 5], sizeof(OVMV) * nb_ctb_pb);
    mv_ctx0->mvs[1 + nb_ctb_pb] = lns->mv0[(ctb_x + 2) << 5];
    mv_ctx1->mvs[1 + nb_ctb_pb] = lns->mv1[(ctb_x + 2) << 5];

    for (i = 1; i < nb_ctb_pb + 1; i++) {
        uint64_t top_available0 = !!(above_map0 & (1llu << (i - 1)));
        uint64_t top_available1 = !!(above_map1 & (1llu << (i - 1)));
        cols_map0[i] = (uint64_t)top_available0;
        cols_map1[i] = (uint64_t)top_available1;
    }

    rows_map0[0] |= above_map0 << 1;
    rows_map1[0] |= above_map1 << 1;

    /* Save last CTU MV line to line at ctb_x */
    memcpy(&lns->mv0[ctb_x << 5], &mv_ctx0->mvs[1 + nb_ctb_pb * 34], sizeof(OVMV) * nb_ctb_pb);
    memcpy(&lns->mv1[ctb_x << 5], &mv_ctx1->mvs[1 + nb_ctb_pb * 34], sizeof(OVMV) * nb_ctb_pb);

    /* Store last inter dir info onto line */
    lns->dir0[ctb_x] = (uint32_t)(lst_row0 >> 1);
    lns->dir1[ctb_x] = (uint32_t)(lst_row1 >> 1);
}

static void
free_dbf_lines(struct DBFLines *const l)
{
    /* FIXME check before free */
    ov_freep(&l->qp_x_map);
    ov_freep(&l->qp_x_map_cb);
    ov_freep(&l->qp_x_map_cr);
    ov_freep(&l->dbf_qp_ver);
    ov_freep(&l->dbf_qp_ver_cb);
    ov_freep(&l->dbf_qp_ver_cr);

    ov_freep(&l->dbf_edge_ver);
    ov_freep(&l->dbf_edge_hor);

    ov_freep(&l->dbf_bs2_ver);
    ov_freep(&l->dbf_bs2_hor);

    ov_freep(&l->dbf_bs1_ver);
    ov_freep(&l->dbf_bs1_hor);
    ov_freep(&l->dbf_bs1_ver_cb);
    ov_freep(&l->dbf_bs1_hor_cb);
    ov_freep(&l->dbf_bs1_ver_cr);
    ov_freep(&l->dbf_bs1_hor_cr);

    ov_freep(&l->large_map_c);
}


static int
init_dbf_lines(struct DBFLines *const l, int nb_ctu_line, int nb_pu_line)
{
    uint8_t malloc_chk = 0;

    l->qp_x_map       = ov_mallocz((nb_pu_line + 1) * sizeof(int8_t));
    l->qp_x_map_cb    = ov_mallocz((nb_pu_line + 1) * sizeof(int8_t));
    l->qp_x_map_cr    = ov_mallocz((nb_pu_line + 1) * sizeof(int8_t));
    l->dbf_qp_ver     = ov_mallocz((nb_pu_line + 1) * sizeof(int8_t));
    l->dbf_qp_ver_cb  = ov_mallocz((nb_pu_line + 1) * sizeof(int8_t));
    l->dbf_qp_ver_cr  = ov_mallocz((nb_pu_line + 1) * sizeof(int8_t));

    l->dbf_edge_ver   = ov_mallocz((nb_ctu_line + 1) * sizeof(uint64_t));
    l->dbf_edge_hor   = ov_mallocz((nb_ctu_line + 1) * sizeof(uint64_t));

    l->dbf_bs1_ver    = ov_mallocz((nb_ctu_line + 1) * sizeof(uint64_t));
    l->dbf_bs1_hor    = ov_mallocz((nb_ctu_line + 1) * sizeof(uint64_t));

    l->large_map_c    = ov_mallocz((nb_ctu_line + 1) * sizeof(uint64_t));

    l->dbf_bs1_ver_cb = ov_mallocz((nb_ctu_line + 1) * sizeof(uint64_t));
    l->dbf_bs1_hor_cb = ov_mallocz((nb_ctu_line + 1) * sizeof(uint64_t));

    l->dbf_bs1_ver_cr = ov_mallocz((nb_ctu_line + 1) * sizeof(uint64_t));
    l->dbf_bs1_hor_cr = ov_mallocz((nb_ctu_line + 1) * sizeof(uint64_t));

    l->dbf_bs2_ver    = ov_mallocz((nb_ctu_line + 1) * sizeof(uint64_t));
    l->dbf_bs2_hor    = ov_mallocz((nb_ctu_line + 1) * sizeof(uint64_t));

    malloc_chk |= l->qp_x_map       == NULL;
    malloc_chk |= l->qp_x_map_cb    == NULL;
    malloc_chk |= l->qp_x_map_cr    == NULL;
    malloc_chk |= l->dbf_qp_ver     == NULL;
    malloc_chk |= l->dbf_qp_ver_cb  == NULL;
    malloc_chk |= l->dbf_qp_ver_cr  == NULL;

    malloc_chk |= l->dbf_edge_ver   == NULL;
    malloc_chk |= l->dbf_edge_hor   == NULL;

    malloc_chk |= l->dbf_bs1_ver    == NULL;
    malloc_chk |= l->dbf_bs1_hor    == NULL;

    malloc_chk |= l->large_map_c    == NULL;

    malloc_chk |= l->dbf_bs1_ver_cb == NULL;
    malloc_chk |= l->dbf_bs1_hor_cb == NULL;

    malloc_chk |= l->dbf_bs1_ver_cr == NULL;
    malloc_chk |= l->dbf_bs1_hor_cr == NULL;

    malloc_chk |= l->dbf_bs2_ver    == NULL;
    malloc_chk |= l->dbf_bs2_hor    == NULL;

    if (malloc_chk) {
        free_dbf_lines(l);
        return OVVC_ENOMEM;
    }
}

static void
dbf_clear_lines(const struct DBFLines *const l, int nb_ctu_line, int nb_pu_line)
{
    memset(l->qp_x_map, 0, (nb_pu_line + 1) * sizeof(int8_t));
    memset(l->qp_x_map_cb, 0, (nb_pu_line + 1) * sizeof(int8_t));
    memset(l->qp_x_map_cr, 0, (nb_pu_line + 1) * sizeof(int8_t));
    memset(l->dbf_qp_ver, 0, (nb_pu_line + 1) * sizeof(int8_t));
    memset(l->dbf_qp_ver_cb, 0, (nb_pu_line + 1) * sizeof(int8_t));
    memset(l->dbf_qp_ver_cr, 0, (nb_pu_line + 1) * sizeof(int8_t));

    memset(l->dbf_edge_ver, 0, (nb_ctu_line + 1) * sizeof(uint64_t));
    memset(l->dbf_edge_hor, 0, (nb_ctu_line + 1) * sizeof(uint64_t));

    memset(l->dbf_bs1_ver, 0, (nb_ctu_line + 1) * sizeof(uint64_t));
    memset(l->dbf_bs1_hor, 0, (nb_ctu_line + 1) * sizeof(uint64_t));

    memset(l->large_map_c, 0, (nb_ctu_line + 1) * sizeof(uint64_t));

    memset(l->dbf_bs1_ver_cb, 0, (nb_ctu_line + 1) * sizeof(uint64_t));
    memset(l->dbf_bs1_hor_cb, 0, (nb_ctu_line + 1) * sizeof(uint64_t));

    memset(l->dbf_bs1_ver_cr, 0, (nb_ctu_line + 1) * sizeof(uint64_t));
    memset(l->dbf_bs1_hor_cr, 0, (nb_ctu_line + 1) * sizeof(uint64_t));

    memset(l->dbf_bs2_ver, 0, (nb_ctu_line + 1) * sizeof(uint64_t));
    memset(l->dbf_bs2_hor, 0, (nb_ctu_line + 1) * sizeof(uint64_t));
}

static void
dbf_load_qp_map(struct DBFInfo *const dbf_info, const struct DBFLines *const l,
                struct DRVLines *l2,
                uint8_t log2_ctu_s, int ctb_x)
{
    uint8_t nb_units_ctb =  1 << (log2_ctu_s & 7) >> 2;
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
        /* FIXME copied from qp ctx*/
        #if 0
        dbf_info->qp_map.hor[2 + i]    = l2->qp_map_x[i >> 1];
        dbf_info->qp_map_cb.hor[2 + i] = l2->qp_map_cb_up[i >> 1];
        dbf_info->qp_map_cr.hor[2 + i] = l2->qp_map_cr_up[i >> 1];
        #else
        dbf_info->qp_map_y.hor [2 + i] = l->qp_x_map   [(ctb_x << 5) + i];
        dbf_info->qp_map_cb.hor[2 + i] = l->qp_x_map_cb[(ctb_x << 5) + i];
        dbf_info->qp_map_cr.hor[2 + i] = l->qp_x_map_cr[(ctb_x << 5) + i];
        #endif
    }

    for (i = 0; i < nb_units_ctb; ++i) {
        dbf_info->qp_map_y.ver [34 * i] = l->dbf_qp_ver   [(ctb_x << 5) + i];
        dbf_info->qp_map_cb.ver[34 * i] = l->dbf_qp_ver_cb[(ctb_x << 5) + i];
        dbf_info->qp_map_cr.ver[34 * i] = l->dbf_qp_ver_cr[(ctb_x << 5) + i];
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

    /* Lst line copied from qp ctx */
    #if 0
    memcpy(&l2->qp_x_map   [ctb_x << 5], dbf_info->qp_map_up,    sizeof(uint8_t) * nb_units_ctb);
    memcpy(&l2->qp_x_map_cb[ctb_x << 5], dbf_info->qp_map_up_cb, sizeof(uint8_t) * nb_units_ctb);
    memcpy(&l2->qp_x_map_cr[ctb_x << 5], dbf_info->qp_map_up_cr, sizeof(uint8_t) * nb_units_ctb);
    #else
    memcpy(&l->qp_x_map   [ctb_x << 5], &dbf_info->qp_map_y.hor[2 + 34*nb_units_ctb],  sizeof(uint8_t) * nb_units_ctb);
    memcpy(&l->qp_x_map_cb[ctb_x << 5], &dbf_info->qp_map_cb.hor[2 + 34*nb_units_ctb], sizeof(uint8_t) * nb_units_ctb);
    memcpy(&l->qp_x_map_cr[ctb_x << 5], &dbf_info->qp_map_cr.hor[2 + 34*nb_units_ctb], sizeof(uint8_t) * nb_units_ctb);

    #if 0
    memset(&dbf_info->qp_map_y.hor[2 + 34*nb_units_ctb], 1, sizeof(uint8_t) * nb_units_ctb);
    memset(&dbf_info->qp_map_cb.hor[2 + 34*nb_units_ctb], 1, sizeof(uint8_t) * nb_units_ctb);
    memset(&dbf_info->qp_map_cr.hor[2 + 34*nb_units_ctb], -1, sizeof(uint8_t) * nb_units_ctb);
    #endif
    #endif

    #if 0
    if (!dbf_info->dbf_disable) {
    #endif
        for (int i = 0; i < nb_units_ctb + 1; i++) {
            l->dbf_qp_ver[(ctb_x << 5) + i]    = dbf_info->qp_map_y.ver [nb_units_ctb + 34 * i];
            l->dbf_qp_ver_cb[(ctb_x << 5) + i] = dbf_info->qp_map_cb.ver[nb_units_ctb + 34 * i];
            l->dbf_qp_ver_cr[(ctb_x << 5) + i] = dbf_info->qp_map_cr.ver[nb_units_ctb + 34 * i];
        }
    #if 0
    }
    #endif
}

static void
dbf_load_edge_map(struct DBFInfo *const dbf_info, const struct DBFLines *const l,
              uint8_t log2_ctu_s, int ctb_x)
{
    int nb_pb_s = (1 << log2_ctu_s) >> 2;
    uint64_t tmp;

    if (ctb_x) {
        /* copy previous 8 vertical edges used for filter length derivation */
        uint64_t *ctb_bnd = dbf_info->ctb_bound_ver;
        memcpy(ctb_bnd, &ctb_bnd[nb_pb_s], 8 * sizeof(uint64_t));
    }

    /*FIXME move to other function */
    tmp = l->dbf_edge_ver[ctb_x];
    dbf_info->ctb_bound_ver[8] |= tmp & 0x1;
    tmp >>= 1;

    for (int i = 1; i < nb_pb_s; ++i) {
        dbf_info->ctb_bound_hor[i + 8] = (-(!!ctb_x)) & (dbf_info->ctb_bound_hor[i + 8] >> (nb_pb_s)) & 0x3;
        dbf_info->ctb_bound_ver[i + 8] = tmp & 0x1;
        tmp >>= 1;
    }

    dbf_info->ctb_bound_ver[8 + nb_pb_s] = (uint64_t) - 1ll;
    dbf_info->ctb_bound_hor[8 + nb_pb_s] = (uint64_t) - 1ll;
    dbf_info->ctb_bound_ver[8] = (uint64_t) - 1ll;
    dbf_info->ctb_bound_hor[8] = (uint64_t) - 1ll;
    dbf_info->ctb_bound_hor[6] = (uint64_t) - 1ll;

    dbf_info->edge_map_ver[0] = (-(!!ctb_x)) & dbf_info->edge_map_ver[nb_pb_s];
    dbf_info->edge_map_hor[0] = (-(!!ctb_x)) & (dbf_info->edge_map_hor[0] >> (nb_pb_s)) & 0x3;
    dbf_info->edge_map_hor[0] |= l->dbf_edge_hor[ctb_x];

    tmp = l->dbf_edge_ver[ctb_x];
    dbf_info->edge_map_ver[0] |= tmp & 0x1;
    tmp >>= 1;

    for (int i = 1; i < nb_pb_s; ++i) {
        dbf_info->edge_map_hor[i] = (-(!!ctb_x)) & (dbf_info->edge_map_hor[i] >> (nb_pb_s)) & 0x3;
        dbf_info->edge_map_ver[i] = tmp & 0x1;
        tmp >>= 1;
    }
}

static void
dbf_store_edge_map(struct DBFInfo *const dbf_info, const struct DBFLines *const l,
                   uint8_t log2_ctu_s, int ctb_x)
{
    int nb_pb_s = (1 << log2_ctu_s) >> 2;
    uint64_t mask = (uint64_t)1 << nb_pb_s;
    int i;

    l->dbf_edge_hor[ctb_x] = dbf_info->edge_map_hor[nb_pb_s];
    l->dbf_edge_ver[ctb_x] = 0;

    for (i = 0; i < nb_pb_s + 1; ++i) {
        l->dbf_edge_ver[ctb_x] |= (!!(dbf_info->edge_map_ver[i] & mask)) << i;
    }

    /* Use last horizontal edges to determine if last CUs were large
     * if there were no edges large in the last three lines will be 0 
     */
    l->large_map_c[ctb_x]  = dbf_info->edge_map_hor[nb_pb_s - 1];
    l->large_map_c[ctb_x] |= dbf_info->edge_map_hor[nb_pb_s - 2];
    l->large_map_c[ctb_x] |= dbf_info->edge_map_hor[nb_pb_s - 3];
}

static void
dbf_load_bs_map(struct DBFInfo *const dbf_info, const struct DBFLines *const l,
            uint8_t log2_ctu_s, int ctb_x)
{
    int nb_pb_s = (1 << log2_ctu_s) >> 2;
    uint64_t tmp;
    struct DBFMap *const bs2_map    = &dbf_info->bs2_map;
    struct DBFMap *const bs1_map    = &dbf_info->bs1_map;
    struct DBFMap *const bs1_map_cb = &dbf_info->bs1_map_cb;
    struct DBFMap *const bs1_map_cr = &dbf_info->bs1_map_cr;
    int i;

    /* FIXME check cast */
    uint64_t ctb_lft_msk = -(!!ctb_x);

    dbf_info->large_map_c = l->large_map_c[ctb_x];

    bs2_map->ver[0] = ctb_lft_msk & bs2_map->ver[nb_pb_s];
    bs1_map->ver[0] = ctb_lft_msk & bs1_map->ver[nb_pb_s];
    bs1_map_cb->ver[0] = ctb_lft_msk & bs1_map_cb->ver[nb_pb_s];
    bs1_map_cr->ver[0] = ctb_lft_msk & bs1_map_cr->ver[nb_pb_s];

    bs2_map->hor[0] = ctb_lft_msk & (bs2_map->hor[0] >> nb_pb_s) & 0x3;
    bs1_map->hor[0] = ctb_lft_msk & (bs1_map->hor[0] >> nb_pb_s) & 0x3;
    bs1_map_cb->hor[0] = ctb_lft_msk & (bs1_map_cb->hor[0] >> nb_pb_s) & 0x3;
    bs1_map_cr->hor[0] = ctb_lft_msk & (bs1_map_cr->hor[0] >> nb_pb_s) & 0x3;

    bs2_map->hor[0] |= l->dbf_bs2_hor[ctb_x];
    bs1_map->hor[0] |= l->dbf_bs1_hor[ctb_x];
    bs1_map_cb->hor[0] |= l->dbf_bs1_hor_cb[ctb_x];
    bs1_map_cr->hor[0] |= l->dbf_bs1_hor_cr[ctb_x];

    /*FIXME separate hor and ver loop  and merge all loops*/
    tmp = l->dbf_bs1_ver[ctb_x];
    bs1_map->ver[0] |= tmp & 0x1;
    tmp >>= 1;

    for (i = 1; i < nb_pb_s + 1; ++i) {
        bs1_map->hor[i] = ctb_lft_msk & (bs1_map->hor[i] >> nb_pb_s) & 0x3;
        bs1_map->ver[i] = tmp & 0x1;
        tmp >>= 1;
    }


    tmp = l->dbf_bs2_ver[ctb_x];
    bs2_map->ver[0] |= tmp & 0x1;
    tmp >>= 1;

    for (i = 1; i < nb_pb_s + 1; ++i) {
        bs2_map->hor[i] = ctb_lft_msk & (bs2_map->hor[i] >> nb_pb_s) & 0x3;
        bs2_map->ver[i] = tmp & 0x1;
        tmp >>= 1;
    }

    tmp = l->dbf_bs1_ver_cb[ctb_x];
    bs1_map_cb->ver[0] |= tmp & 0x1;
    tmp >>= 1;

    for (i = 1; i < nb_pb_s + 1; ++i) {
        bs1_map_cb->hor[i] = ctb_lft_msk & (bs1_map_cb->hor[i] >> nb_pb_s) & 0x3;
        bs1_map_cb->ver[i] = tmp & 0x1;
        tmp >>= 1;
    }


    tmp = l->dbf_bs1_ver_cr[ctb_x];
    bs1_map_cr->ver[0] |= tmp & 0x1;
    tmp >>= 1;

    for (i = 1; i < nb_pb_s + 1; ++i) {
        bs1_map_cr->hor[i] = ctb_lft_msk & (bs1_map_cr->hor[i] >> nb_pb_s) & 0x3;
        bs1_map_cr->ver[i] = tmp & 0x1;
        tmp >>= 1;
    }
}

static void
dbf_store_bs_map(struct DBFInfo *const dbf_info, const struct DBFLines *const l,
              uint8_t log2_ctu_s, int ctb_x)
{
    int nb_pb_s = (1 << log2_ctu_s) >> 2;
    uint64_t mask = (uint64_t)1 << nb_pb_s;
    int i;

    struct DBFMap *const bs2_map = &dbf_info->bs2_map;
    struct DBFMap *const bs1_map = &dbf_info->bs1_map;
    struct DBFMap *const bs1_map_cb = &dbf_info->bs1_map_cb;
    struct DBFMap *const bs1_map_cr = &dbf_info->bs1_map_cr;

    l->dbf_bs2_hor[ctb_x] =  bs2_map->hor[nb_pb_s];
    l->dbf_bs1_hor[ctb_x] =  bs1_map->hor[nb_pb_s];
    l->dbf_bs1_hor_cb[ctb_x] =  bs1_map_cb->hor[nb_pb_s];
    l->dbf_bs1_hor_cr[ctb_x] =  bs1_map_cr->hor[nb_pb_s];

    /* Reset vertical edge maps on line storage*/
    l->dbf_bs2_ver[ctb_x] = 0;
    l->dbf_bs1_ver[ctb_x] = 0;
    l->dbf_bs1_ver_cb[ctb_x] = 0;
    l->dbf_bs1_ver_cr[ctb_x] = 0;

    /* if the i got an vertical edge on the last pb store it */
    for (i = 0; i < nb_pb_s; ++i) {
        l->dbf_bs2_ver[ctb_x] |= (!!(bs2_map->ver[i] & mask)) << i;
        l->dbf_bs1_ver[ctb_x] |= (!!(bs1_map->ver[i] & mask)) << i;
        l->dbf_bs1_ver_cb[ctb_x] |= (!!(bs1_map_cb->ver[i] & mask)) << i;
        l->dbf_bs1_ver_cr[ctb_x] |= (!!(bs1_map_cr->ver[i] & mask)) << i;
    }
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
     const OVPartInfo *const pinfo = sldec->ctudec_list->part_ctx;
     const OVSPS *const sps = prms->sps;

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
     uint16_t nb_pb_pic_w = nb_ctb_pic_w << (log2_ctb_s - log2_min_cb_s);

     ret = init_inter_drv_lines(lns, nb_ctb_pic_w, nb_ctb_pic_w);

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
    const OVPartInfo *pinfo = sldec->ctudec_list->part_ctx;
    const OVSPS *const sps = prms->sps;

    struct DRVLines *const lns = &sldec->drv_lines;

    uint8_t log2_ctb_s = pinfo->log2_ctu_s;
    uint8_t log2_min_cb_s = pinfo->log2_min_cb_s;
    /* TODO use active parameters such as generic pic info
     * see init_cabac_lines
     */
    uint16_t pic_w = sps->sps_pic_width_max_in_luma_samples;
    uint16_t nb_ctb_pic_w = (pic_w + ((1 << log2_ctb_s) - 1)) >> log2_ctb_s;
    uint16_t nb_pb_pic_w = nb_ctb_pic_w << (log2_ctb_s - log2_min_cb_s);

    struct InterLines *const i_lns = &lns->inter_lines;

    memset(i_lns->dir0, 0, (sizeof(*i_lns->dir0) * (nb_ctb_pic_w + 1)));
    memset(i_lns->dir1, 0, (sizeof(*i_lns->dir1) * (nb_ctb_pic_w + 1)));

     /* PLANAR  = 0 value is used if absent so we use it as reset value
      */
     memset(lns->intra_luma_x,     0,  sizeof(*lns->intra_luma_x) * nb_pb_pic_w);
     dbf_clear_lines(&lns->dbf_lines, nb_ctb_pic_w,  nb_pb_pic_w);
}

static void
load_first_ctu_inter(struct DRVLines *const l,
                 const OVCTUDec *const ctudec,
                 unsigned int ctb_x)
{
    uint8_t nb_ctb_pb =  (1 << ((ctudec->part_ctx->log2_ctu_s) & 7)) >> ctudec->part_ctx->log2_min_cb_s;

    struct InterDRVCtx *const inter_ctx = &ctudec->drv_ctx.inter_ctx;
    const struct InterLines *const lns = &l->inter_lines;

    struct OVMVCtx *const mv_ctx0 = &inter_ctx->mv_ctx0;
    struct OVMVCtx *const mv_ctx1 = &inter_ctx->mv_ctx1;

    uint64_t *const rows_map0 = inter_ctx->mv_ctx0.map.hfield;
    uint64_t *const cols_map0 = inter_ctx->mv_ctx0.map.vfield;

    uint64_t *const rows_map1 = inter_ctx->mv_ctx1.map.hfield;
    uint64_t *const cols_map1 = inter_ctx->mv_ctx1.map.vfield;

    int i;

    uint64_t above_map0 = (uint64_t)lns->dir0[0] | ((uint64_t)lns->dir0[1] << nb_ctb_pb);
    uint64_t above_map1 = (uint64_t)lns->dir1[0] | ((uint64_t)lns->dir1[1] << nb_ctb_pb);

    for (i = 0; i < nb_ctb_pb + 1; i++) {
        uint64_t left_available0 = !!(cols_map0[nb_ctb_pb] & (1llu << i));
        uint64_t left_available1 = !!(cols_map1[nb_ctb_pb] & (1llu << i));
        rows_map0[i] = 0;
        rows_map1[i] = 0;
        mv_ctx0->mvs[i * 34] = mv_ctx0->mvs[i* 34 + nb_ctb_pb];
        mv_ctx1->mvs[i * 34] = mv_ctx1->mvs[i* 34 + nb_ctb_pb];
    }
    cols_map0[0] = 0;
    cols_map1[0] = 0;

    /* Replace CTU above MV line by line MV at ctb_x + 1*/
    memcpy(&mv_ctx0->mvs[1], &lns->mv0[0], sizeof(OVMV) * nb_ctb_pb);
    memcpy(&mv_ctx1->mvs[1], &lns->mv1[0], sizeof(OVMV) * nb_ctb_pb);

    mv_ctx0->mvs[1 + nb_ctb_pb] = lns->mv0[32];
    mv_ctx1->mvs[1 + nb_ctb_pb] = lns->mv1[32];

    for (i = 1; i < nb_ctb_pb + 1; i++) {
        uint64_t top_available0 = !!(above_map0 & (1llu << (i - 1)));
        uint64_t top_available1 = !!(above_map1 & (1llu << (i - 1)));
        cols_map0[i] = top_available0;
        cols_map1[i] = top_available1;
    }
    rows_map0[0] |= above_map0 << 1;
    rows_map1[0] |= above_map1 << 1;

    /* Reset HMVP Look Up table */
    inter_ctx->hmvp_lut.nb_mv = 0;
}


void
drv_line_next_line(OVCTUDec *const ctudec, const OVSliceDec *const sldec)
{
    #if 0
    struct OVDrvCtx *const drv_ctx        = &ctudec->drv_ctx;
    #endif
    struct IntraDRVInfo *const intra_info = &ctudec->drv_ctx.intra_info;
    const struct DRVLines *const lns      = &sldec->drv_lines;
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

    #if 1
    memset(&ctudec->dbf_info, 0, sizeof(ctudec->dbf_info));
    #endif
    dbf_load_info(&ctudec->dbf_info, &sldec->drv_lines.dbf_lines, log2_ctb_s, 0);
}


void
drv_line_next_ctu(OVCTUDec *const ctudec, OVSliceDec *sldec, struct DRVLines *drv_lns,
                  const OVPS *const prms, uint16_t ctb_x)
{
    const OVPartInfo *const pinfo = ctudec->part_ctx;
    struct IntraDRVInfo *const intra_info = &ctudec->drv_ctx.intra_info;
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

void
drv_lines_uninit(OVSliceDec *sldec)
{
     struct DRVLines *const lns = &sldec->drv_lines;

     ov_freep(&lns->intra_luma_x);
     free_inter_drv_lines(lns);
     free_dbf_lines(&lns->dbf_lines);
}
