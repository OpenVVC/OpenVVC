#include "ovunits.h"
#include "ovmem.h"
#include "ovutils.h"

int
ov_nalu_init(OVNALUnit *nalu)
{
    nalu->rbsp_data = NULL;
    nalu->rbsp_size = 0;

    nalu->epb_pos = NULL;
    nalu->nb_epb = 0;

    atomic_init(&nalu->ref_count, 0);
    return 1;
}

int
ov_nalu_new_ref(OVNALUnit **nalu_p, OVNALUnit *nalu)
{
    if (!nalu) {
        return -1;
    }

    atomic_fetch_add_explicit(&nalu->ref_count, 1, memory_order_acq_rel);

    *nalu_p = nalu;

    return 0;
}

static void
ovnalu_free(OVNALUnit *nalu)
{
    ov_freep(&nalu->rbsp_data);
    if (nalu->epb_pos) {
        ov_freep(&nalu->epb_pos);
    }
    ov_free(nalu);
}

void
ov_nalu_unref(OVNALUnit **nalu_p)
{
    if (!nalu_p)
        return;

    OVNALUnit *nalu = *nalu_p;

    if (!nalu){
        ov_log(NULL, OVLOG_ERROR, "Trying to unref NULL nalu\n");
        return;
    }

    unsigned ref_count = atomic_fetch_add_explicit(&nalu->ref_count, -1, memory_order_acq_rel);

    if (!ref_count) {
        ovnalu_free(nalu);
    }

    *nalu_p = NULL;
}

/*FIXME Add and reference counting */
void
ov_free_pu(OVPictureUnit **pu)
{
    OVPictureUnit *to_free = *pu;
    if (to_free) {
        int i;
        for (i = 0; i < to_free->nb_nalus; ++i) {
            OVNALUnit *nalu = to_free->nalus[i];

            ov_nalu_unref(&nalu);
        }
        ov_free(to_free->nalus);
    }
    ov_freep(pu);
}
