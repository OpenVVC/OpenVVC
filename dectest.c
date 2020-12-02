#include <stddef.h>
#include <stdio.h>

#include "libovvcdec/ovvcdec.h"
#include "libovvcdmx/ovvcdmx.h"
#include "libovvcutils/ovvcutils.h"

typedef struct OVVCHdl{
    OVVCDmx *dmx;
    OVVCDec *dec;
    /* TODO decide whether or not file pointer must be given
       to dmx or not  it is only given to the hadle so when can
       close the file later since we opened it from here*/
    FILE *fp;
}OVVCHdl;

/* Use global static value so we can easily change log level on
   the whole object this is useful when debugging */
static int log_level = 2;

static int dmx_attach_file(OVVCHdl *const vvc_hdl, const char *const file_name);

static int init_openvvc_hdl(OVVCHdl *const ovvc_hdl);

static int close_openvvc_hdl(OVVCHdl *const ovvc_hdl);

int
main(int argc, char** argv)
{
    OVVCHdl ovvc_hdl;
    const char *file_name = "test.266";
    int ret;

    /* TODO add basic opitons parser and assign
       filenames intoa functions*/

    if (argc > 2) {
        fprintf(stderr, "Too many arguments.\n");
        /* TODO show_usage();*/
        return -1;
    } else if (argc == 2) {
        file_name = argv[1];
    }

    ret = init_openvvc_hdl(&ovvc_hdl);

    if (ret < 0) goto failinit;

    ret = dmx_attach_file(&ovvc_hdl, file_name);

    if (ret < 0) goto failattach;

    /* Do stuff here */

failattach:
    ret = close_openvvc_hdl(&ovvc_hdl);

failinit:

    return ret;
}

static int
dmx_attach_file(OVVCHdl *const vvc_hdl, const char *const file_name)
{
    FILE *file = fopen(file_name,"rb");

    if (file == NULL) {
        perror(file_name);
       vvc_hdl->fp = NULL;
       return -1;
    }

    vvc_hdl->fp = file;

    /*TOOO Call ret = ovvc_dmx_attach */

    return 0;
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

    if (ovvc_hdl->fp != NULL) {
        fclose(ovvc_hdl->fp);
    }

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

