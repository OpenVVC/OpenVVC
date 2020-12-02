#include <stddef.h>

#include "libovvcdec/ovvcdec.h"
#include "libovvcdmx/ovvcdmx.h"
#include "libovvcutils/ovvcutils.h"

typedef struct OVVCHdl{
    OVVCDmx *dmx;
    OVVCDec *dec;
}OVVCHdl;

/* Use global static value so we can easily change log level on
   the whole object this is useful when debugging */
static int log_level = 2;

static int init_openvvc_hdl(OVVCHdl *const ovvc_hdl);

static int close_openvvc_hdl(OVVCHdl *const ovvc_hdl);

int
main(int argc, char** argv)
{
    OVVCHdl ovvc_hdl;
    int ret;


    ret = init_openvvc_hdl(&ovvc_hdl);

    if (ret < 0) goto failinit;


    if (ret < 0) goto failattach;

    /* Do stuff here */


failattach:
    ret = close_openvvc_hdl(&ovvc_hdl);

failinit:

    return ret;
}

static int
init_openvvc_hdl(OVVCHdl *const ovvc_hdl)
{
    OVVCDec **vvcdec = &ovvc_hdl->dec;
    OVVCDmx **vvcdmx = &ovvc_hdl->dmx;
    int ret;

    ret = ovdec_init(vvcdec);

    if (ret < 0) goto faildec;

    ov_log(vvcdec, log_level, "Decoder init.\n");

    ret = ovdmx_init(vvcdmx);

    if (ret < 0) goto faildmx;

    ov_log(vvcdmx, log_level, "Demuxer init.\n");

    return 0;

faildec:
    ov_log(NULL, log_level, "Decoder failed at init.\n");
    return ret;

faildmx:
    ov_log(NULL, log_level, "Demuxer failed at init.\n");
    ovdec_close(*vvcdec);
    return ret;
}

static int
close_openvvc_hdl(OVVCHdl *const ovvc_hdl)
{
    OVVCDec *vvcdec = ovvc_hdl->dec;
    OVVCDmx *vvcdmx = ovvc_hdl->dmx;
    int ret;

    ret = ovdec_close(vvcdec);

    if (ret < 0) goto faildecclose;

    ret = ovdmx_close(vvcdmx);

    if (ret < 0) goto faildmxclose;

    return 0;

faildecclose:
    /* Do not check for dmx failure  since it might override
       return value to a correct one in either success or 
       failure we already raised an error*/
    ov_log(NULL, log_level, "Decoder failed at cloture.\n");
    ovdmx_close(vvcdmx);

faildmxclose:
    return ret;
}
