#ifndef DBF_UTILS_H
#define DBF_UTILS_H

/*FIXME harmonize notations */
/* Store average qp between two CB into line */
/*FIXME should we really inline this one */
static inline void
dbf_fill_qp_map(struct DBFQPMap *qp_map, int x0, int y0,
                int log2_cb_w, int log2_cb_h, int8_t qp)
{
    uint8_t x0_u = x0 >> 2;
    uint8_t y0_u = y0 >> 2;
    int nb_cb_w = (1 << log2_cb_w) >> 2;
    int nb_cb_h = (1 << log2_cb_h) >> 2;

    int first_pos_hor = 2 + x0_u + y0_u * 34;
    int first_pos_ver = y0_u + x0_u * 34;
    int i;

    for (i = 0; i < nb_cb_w; i++) {
        qp_map->hor[first_pos_hor + i] += qp + 1;
        qp_map->hor[first_pos_hor + i] >>= 1;
        qp_map->hor[first_pos_hor + i + 34 * nb_cb_h] = qp;
    }

    /* we actually need 2 maps because conflict between vertical
     * and horizontal average QP could occur in map corners.
     */
    for (i = 0; i < nb_cb_h; i++) {
        qp_map->ver[first_pos_ver + i] += qp + 1;
        qp_map->ver[first_pos_ver + i] >>= 1;
        qp_map->ver[first_pos_ver + i + 34 * nb_cb_w] = qp;
    }
}

/* Set nb_unit_w/h bits to 1 onto left/above and right/bottom parts of a CB in BS maps*/
static inline void
fill_bs_map(struct DBFMap *const dbf_map, int x0, int y0, int log2_cu_w, int log2_cu_h)
{
    const uint64_t mask_ver = (1 << ((1 << log2_cu_h) >> 2)) - 1;
    const uint64_t mask_hor = (1 << ((1 << log2_cu_w) >> 2)) - 1;
    uint8_t x0_u = x0 >> 2;
    uint8_t y0_u = y0 >> 2;

    dbf_map->ver[(x0 + (1 << log2_cu_w)) >> 2] |= mask_ver << y0_u;
    dbf_map->hor[(y0 + (1 << log2_cu_h)) >> 2] |= mask_hor << (2 + x0_u);

    dbf_map->ver[x0_u] |= mask_ver << y0_u;
    dbf_map->hor[y0_u] |= mask_hor << (2 + x0_u);
}

/* Set nb_unit_w/h bits to 1 onto right/bottom part of a CB in edge_map  */
static inline void
fill_ctb_bound(struct DBFInfo *const dbf_info, int x0, int y0, int log2_cu_w, int log2_cu_h)
{
    const uint64_t mask_ver = (1 << ((1 << log2_cu_h) >> 2)) - 1;
    const uint64_t mask_hor = (1 << ((1 << log2_cu_w) >> 2)) - 1;
    uint8_t x0_u = x0 >> 2;
    uint8_t y0_u = y0 >> 2;

    dbf_info->ctb_bound_ver[8 + ((x0 + (1 << log2_cu_w)) >> 2)] |= mask_ver << y0_u;
    dbf_info->ctb_bound_hor[8 + ((y0 + (1 << log2_cu_h)) >> 2)] |= mask_hor << (2 + x0_u);

    dbf_info->ctb_bound_ver[8 + x0_u] |= mask_ver << y0_u;
    dbf_info->ctb_bound_hor[8 + y0_u] |= mask_hor << (2 + x0_u);
}

static inline void
fill_ctb_bound_c(struct DBFInfo *const dbf_info, int x0, int y0, int log2_cu_w, int log2_cu_h)
{
    const uint64_t mask_ver = (1 << ((1 << log2_cu_h) >> 2)) - 1;
    const uint64_t mask_hor = (1 << ((1 << log2_cu_w) >> 2)) - 1;
    uint8_t x0_u = x0 >> 2;
    uint8_t y0_u = y0 >> 2;

    dbf_info->ctb_bound_ver_c[8 + ((x0 + (1 << log2_cu_w)) >> 2)] |= mask_ver << y0_u;
    dbf_info->ctb_bound_hor_c[8 + ((y0 + (1 << log2_cu_h)) >> 2)] |= mask_hor << (2 + x0_u);

    dbf_info->ctb_bound_ver_c[8 + x0_u] |= mask_ver << y0_u;
    dbf_info->ctb_bound_hor_c[8 + y0_u] |= mask_hor << (2 + x0_u);
}

#endif
