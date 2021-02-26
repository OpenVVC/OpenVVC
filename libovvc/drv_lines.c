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

    uint64_t *const rows_map0 = inter_ctx->mv_ctx0.map.hfield;
    uint64_t *const cols_map0 = inter_ctx->mv_ctx0.map.vfield;

    uint64_t *const rows_map1 = inter_ctx->mv_ctx1.map.hfield;
    uint64_t *const cols_map1 = inter_ctx->mv_ctx1.map.vfield;

    int i;

    const uint64_t lst_row0 = rows_map0[nb_ctb_pb];
    const uint64_t lst_col0 = cols_map0[nb_ctb_pb];

    const uint64_t lst_row1 = rows_map1[nb_ctb_pb];
    const uint64_t lst_col1 = cols_map1[nb_ctb_pb];

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
    #if 0
    rows_map0[0] |= cols_map0[0] & 0x1;
    rows_map1[0] |= cols_map1[0] & 0x1;
    #else
    rows_map0[0] |= above_map0 << 1;
    rows_map1[0] |= above_map1 << 1;
    #endif

    /* Save last CTU MV line to line at ctb_x */
    memcpy(&lns->mv0[ctb_x << 5], &mv_ctx0->mvs[1 + nb_ctb_pb * 34], sizeof(OVMV) * nb_ctb_pb);
    memcpy(&lns->mv1[ctb_x << 5], &mv_ctx1->mvs[1 + nb_ctb_pb * 34], sizeof(OVMV) * nb_ctb_pb);

    /* Store last inter dir info onto line */
    lns->dir0[ctb_x] = (uint32_t)(lst_row0 >> 1);
    lns->dir1[ctb_x] = (uint32_t)(lst_row1 >> 1);
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
    intra_info->luma_mode_x = lns->intra_luma_x;

    /* Reset HMVP Look Up table */
    ctudec->drv_ctx.inter_ctx.hmvp_lut.nb_mv = 0;
    load_first_ctu_inter(lns, ctudec, 0);

    memset(intra_info->luma_mode_y, 0, sizeof(*intra_info->luma_mode_y) * nb_pb_ctb_w);
    memset(drv_ctx->qp_map_x, qp_val, sizeof(*drv_ctx->qp_map_x) * nb_pb_ctb_w);
    memset(drv_ctx->qp_map_y, qp_val, sizeof(*drv_ctx->qp_map_y) * nb_pb_ctb_w);
}

void
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
}

void
drv_line_next_ctu(OVCTUDec *const ctudec, struct DRVLines *drv_lns,
                  const OVPS *const prms, uint16_t ctb_x)
{
    const OVPartInfo *const pinfo = ctudec->part_ctx;
    struct IntraDRVInfo *const intra_info = &ctudec->drv_ctx.intra_info;
    struct OVDrvCtx *const drv_ctx = &ctudec->drv_ctx;

    uint8_t log2_ctb_s    = pinfo->log2_ctu_s;
    uint8_t log2_min_cb_s = pinfo->log2_min_cb_s;

    uint16_t nb_pb_ctb_w = (1 << log2_ctb_s) >> log2_min_cb_s;

    int8_t qp_val = ctudec->qp_ctx.current_qp;

        #if 0
        store_inter_maps(drv_lns, ctudec, ctb_x);
        #endif

#if 0
    update_inter_ctu_dec(ctudec, &drv_lns->inter_lines, 0, ctb_x);
#endif

    /* FIXME Unecessary ? */
    memset(intra_info->luma_mode_x, 0, sizeof(*intra_info->luma_mode_x) * nb_pb_ctb_w);

    /* Reset QP prediction map to previous current QP prediction
     * value
     */
    memset(drv_ctx->qp_map_x, qp_val, sizeof(*drv_ctx->qp_map_x) * nb_pb_ctb_w);
    memset(drv_ctx->qp_map_y, qp_val, sizeof(*drv_ctx->qp_map_y) * nb_pb_ctb_w);
}

void
drv_lines_uninit(OVSliceDec *sldec)
{
     struct DRVLines *const lns = &sldec->drv_lines;

     ov_freep(&lns->intra_luma_x);
     free_inter_drv_lines(lns);
}
