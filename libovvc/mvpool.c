#include "ovdefs.h"
#include "ovdec_internal.h"
#include "ovmem.h"
#include "ovutils.h"
#include "mempool.h"
#include "mempool_internal.h"
#include "overror.h"

/*FIXME only for MV declaration which should be
 * declared somewhere else
 */
#include "ctudec.h"


#define LOG2_MIN_CU_S 2

static int
init_mv_pool(struct MVPool *const mv_pool, const struct PicPartInfo *const pinfo)
{
    size_t   nb_ctb_pic = (size_t) pinfo->nb_ctb_w * pinfo->nb_ctb_h;
    uint16_t nb_pb_ctb_w = (1 << pinfo->log2_ctu_s) >> LOG2_MIN_CU_S;

    /*FIXME tmp size to ben reduced to 8x8*/
    size_t elem_size = nb_ctb_pic * sizeof(OVMV) * nb_pb_ctb_w * nb_pb_ctb_w;

    mv_pool->mv_pool = ovmempool_init(elem_size);

    if (!mv_pool->dir_pool) {
       return OVVC_ENOMEM;
    }

    return 0;
}

static int
init_dir_field_pool(struct MVPool *const mv_pool, const struct PicPartInfo *const pinfo)
{
    size_t nb_ctb_pic = (size_t) pinfo->nb_ctb_w * pinfo->nb_ctb_h;
    uint8_t nb_pb_ctb_w = (1 << pinfo->log2_ctu_s) >> LOG2_MIN_CU_S;
    size_t elem_size = nb_ctb_pic * sizeof(uint64_t) * nb_pb_ctb_w;

    mv_pool->dir_pool = ovmempool_init(elem_size);

    if (!mv_pool->dir_pool) {
       return OVVC_ENOMEM;
    }

    return 0;
}

int
mvpool_init(struct MVPool **mv_pool_p, const struct PicPartInfo *const pinfo)
{
    struct MVPool *mv_pool;
    int ret;

    mv_pool = ov_mallocz(sizeof(*mv_pool));
    if (!mv_pool) {
        return OVVC_ENOMEM;
    }

    *mv_pool_p = mv_pool;

    ret = init_dir_field_pool(mv_pool, pinfo);
    if (ret < 0) {
        goto fail_field;
    }

    ret = init_mv_pool(mv_pool, pinfo);
    if (ret < 0) {
        goto fail_mv;
    }

    return 0;

fail_mv :
    ovmempool_uninit(&mv_pool->dir_pool);

fail_field :
    ov_freep(mv_pool_p);
    ov_log(NULL, OVLOG_ERROR, "MV pool intialisation failed\n");
    return OVVC_ENOMEM;
}

void
mvpool_uninit(struct MVPool **mv_pool_p)
{
    struct MVPool *mv_pool = *mv_pool_p;

    ovmempool_uninit(&mv_pool->dir_pool);

    ovmempool_uninit(&mv_pool->mv_pool);

    ov_freep(mv_pool_p);
}

int
mvpool_request_mv_plane(struct MVPool *mv_pool, struct MVPlane *mv_plane)
{
    MemPoolElem *mv_elem;
    MemPoolElem *dir_elem;

    mv_elem  = ovmempool_popelem(mv_pool->mv_pool);
    if (!mv_elem) {
        goto fail_mv;
    }

    dir_elem  = ovmempool_popelem(mv_pool->dir_pool);
    if (!dir_elem) {
        goto fail_dir;
    }

    mv_plane->dir_elem = dir_elem;
    mv_plane->mv_elem  = mv_elem;

    mv_plane->mvs  = mv_elem->data;
    mv_plane->dirs = dir_elem->data;

    return 0;

fail_dir:
    ovmempool_pushelem(mv_plane->mv_elem);
fail_mv:
    return OVVC_ENOMEM;
}

void
mvpool_release_mv_plane(struct MVPlane *mv_plane)
{
    ovmempool_pushelem(mv_plane->mv_elem);

    ovmempool_pushelem(mv_plane->dir_elem);

    mv_plane->mvs  = NULL;

    mv_plane->dirs = NULL;
}
