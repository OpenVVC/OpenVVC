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

    ov_log(vvcdec, log_level, "Decoder init.\n");

    ret = ovdmx_init(&vvcdmx);

    if (ret < 0) goto faildmx;

    ov_log(vvcdmx, log_level, "Demuxer init.\n");

    /* Do stuff here */

    ret = ovdec_close(vvcdec);

    if (ret < 0) goto faildecclose;

    ret = ovdmx_close(vvcdmx);

    if (ret < 0) goto faildmxclose;

    return 0;

faildec:
    ov_log(NULL, log_level, "Decoder failed at init.\n");
    return ret;

faildmx:
    ov_log(NULL, log_level, "Demuxer failed at init.\n");
    ovdec_close(vvcdec);
    return ret;

faildecclose:
    ov_log(NULL, log_level, "Demuxer failed at cloture.\n");
    ovdmx_close(vvcdmx);

faildmxclose:
    return ret;
}
