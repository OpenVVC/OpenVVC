#include <stddef.h>

#include "libovvcdec/ovvcdec.h"
#include "libovvcdmx/ovvcdmx.h"
#include "libovvcutils/ovvcutils.h"


int
main(int argc, char** argv)
{
    OVVCDec *vvcdec;
    OVVCDmx *vvcdmx;
    int log_level = 2;
    int ret;

    ret = ovdec_init(&vvcdec);

    if (ret < 0) goto faildec;

    ov_log(vvcdec, log_level, "Success at decoder init.\n");

    ret = ovdmx_init(&vvcdmx);

    if (ret < 0) goto faildmx;

    ov_log(vvcdmx, log_level, "Success at demuxer init.\n");

    /* Do stuff here */

    ovdec_close(&vvcdec);
    ovdmx_close(&vvcdmx);

    return 0;

faildec:
#if 1
    ov_log(NULL, log_level, "Error at decoder init.\n");
#endif
    return ret;

faildmx:
#if 1
    ov_log(NULL, log_level, "Error at demuxer init.\n");
#endif
    ovdec_close(&vvcdec);
    return ret;
}
