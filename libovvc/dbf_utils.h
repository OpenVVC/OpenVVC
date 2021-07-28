#ifndef DBF_UTILS_H
#define DBF_UTILS_H

/*FIXME harmonize notations */
/* Store average qp between two CB into line */
/*FIXME should we really inline this one */
static inline void
dbf_fill_qp_map(struct DBFQPMap *qp_map, int x0, int y0,
                int log2_cb_w, int log2_cb_h, int8_t qp)
{
    /* FIXME use border storage */
    uint8_t x0_u = x0 >> 2;
    uint8_t y0_u = y0 >> 2;
    int nb_cb_w = (1 << log2_cb_w) >> 2;
    int nb_cb_h = (1 << log2_cb_h) >> 2;

    int first_pos_hor = 34 + 2 + x0_u + y0_u * 34;
    int i;

    for (i = 0; i < nb_cb_h; i++) {
        memset(&qp_map->hor[first_pos_hor + 34 * i], qp, sizeof(uint8_t) * nb_cb_w);
    }
}

static inline void
dbf_fill_aff_map(struct DBFMap *aff_map, int x0_u, int y0_u,
                int nb_unit_w, int nb_unit_h)
{
    /* FIXME avoid duplicate storage with MVs */
    uint64_t mask_ver = (uint64_t)((uint64_t)1 << nb_unit_h) - 1;
    uint64_t mask_hor = (uint64_t)((uint64_t)1 << nb_unit_w) - 1;
    int i;

    mask_ver <<= y0_u;
    mask_hor <<= x0_u + 2;

    for (i = 0; i < nb_unit_w; i++) {
        aff_map->ver[x0_u + i + 1] |= mask_ver;
    }

    for (i = 0; i < nb_unit_h; i++) {
        aff_map->hor[y0_u + i + 1] |= mask_hor;
    }
}

/* Set nb_unit_w/h bits to 1 onto left/above and right/bottom parts of a CB in BS maps*/
static inline void
fill_bs_map(struct DBFMap *const dbf_map, int x0, int y0, int log2_cu_w, int log2_cu_h)
{
    uint16_t nb_unit_w = (1 << log2_cu_w) >> 2;
    uint16_t nb_unit_h = (1 << log2_cu_h) >> 2;
    const uint64_t mask_ver = (uint64_t)((uint64_t)1 << nb_unit_h) - 1;
    const uint64_t mask_hor = (uint64_t)((uint64_t)1 << nb_unit_w) - 1;
    uint8_t x0_u = x0 >> 2;
    uint8_t y0_u = y0 >> 2;

    dbf_map->ver[x0_u + nb_unit_w] |= mask_ver << y0_u;
    dbf_map->hor[y0_u + nb_unit_h] |= mask_hor << (2 + x0_u);

    dbf_map->ver[x0_u] |= mask_ver << y0_u;
    dbf_map->hor[y0_u] |= mask_hor << (2 + x0_u);
}

/* Set nb_unit_w/h bits to 1 onto right/bottom part of a CB in edge_map  */
static inline void
fill_ctb_bound(struct DBFInfo *const dbf_info, int x0, int y0, int log2_cu_w, int log2_cu_h)
{
    uint16_t nb_unit_w = (1 << log2_cu_w) >> 2;
    uint16_t nb_unit_h = (1 << log2_cu_h) >> 2;

    const uint64_t mask_ver = (uint64_t)((uint64_t)1 << nb_unit_h) - 1;
    const uint64_t mask_hor = (uint64_t)((uint64_t)1 << nb_unit_w) - 1;
    uint8_t x0_u = x0 >> 2;
    uint8_t y0_u = y0 >> 2;

    dbf_info->ctb_bound_ver[8 + x0_u + nb_unit_w] |= mask_ver << y0_u;
    dbf_info->ctb_bound_hor[8 + y0_u + nb_unit_h] |= mask_hor << (2 + x0_u);

    dbf_info->ctb_bound_ver[8 + x0_u] |= mask_ver << y0_u;
    dbf_info->ctb_bound_hor[8 + y0_u] |= mask_hor << (2 + x0_u);
}

static inline void
fill_ctb_bound_c(struct DBFInfo *const dbf_info, int x0, int y0, int log2_cu_w, int log2_cu_h)
{
    uint16_t nb_unit_w = (1 << log2_cu_w) >> 2;
    uint16_t nb_unit_h = (1 << log2_cu_h) >> 2;
    const uint64_t mask_ver = (uint64_t)((uint64_t)1 << nb_unit_h) - 1;
    const uint64_t mask_hor = (uint64_t)((uint64_t)1 << nb_unit_w) - 1;
    uint8_t x0_u = x0 >> 2;
    uint8_t y0_u = y0 >> 2;

    dbf_info->ctb_bound_ver_c[8 + x0_u + nb_unit_w] |= mask_ver << y0_u;
    dbf_info->ctb_bound_hor_c[8 + y0_u + nb_unit_h] |= mask_hor << (2 + x0_u);

    dbf_info->ctb_bound_ver_c[8 + x0_u] |= mask_ver << y0_u;
    dbf_info->ctb_bound_hor_c[8 + y0_u] |= mask_hor << (2 + x0_u);
}

#endif
