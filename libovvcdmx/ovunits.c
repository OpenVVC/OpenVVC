#include "ovunits.h"
#include "libovvcutils/ovmem.h"


int
ov_init_nalu()
{
    return 1;
}

/*FIXME Add and reference counting */
void
ov_free_pu(OVPictureUnit **pu)
{
    OVPictureUnit *to_free = *pu;
    if (to_free) {
        int i;
        for (i = 0; i < to_free->nb_nalus; ++i) {
            OVNALUnit *nalu = &to_free->nalus[i];

            ov_freep(&nalu->rbsp_data);

            if (nalu->epb_pos) {
                ov_freep(&nalu->epb_pos);
            }
        }
        ov_free(to_free->nalus);
    }
    ov_freep(pu);
}
