#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>

#include "ovdec.h"
#include "ovio.h"
#include "ovdefs.h"
#include "ovdmx.h"
#include "ovframe.h"
#include "ovdpb.h"
#include "ovutils.h"
#include "ovversion.h"
#include "ovmem.h"

typedef struct OVVCHdl
{
    OVVCDmx *dmx;
    OVVCDec *dec;
    OVIO *io;
} OVVCHdl;

static int dmx_attach_file(OVVCHdl *const vvc_hdl, const char *const input_file_name);

static int init_openvvc_hdl(OVVCHdl *const ovvc_hdl, const char *output_file_name, int nb_frame_th, int nb_entry_th, int upscale_flag);

static int close_openvvc_hdl(OVVCHdl *const ovvc_hdl);

static int read_write_stream(OVVCHdl *const hdl, FILE *fout);

static int write_decoded_frame_to_file(OVFrame *const frame, FILE *fp, int* max_frame_h);

static void print_version(void);

static void print_usage(void);


int
main(int argc, char** argv)
{
    /* basic options parser and assign
       filenames into a functions*/
    int c;
    const char *input_file_name = NULL, *output_file_name = NULL;
    int ov_log_level = OVLOG_INFO;
    FILE *fout = NULL;
    int nb_frame_th = 0;
    int nb_entry_th = 0;
    int upscale_flag = 0;

    uint8_t options_flag=0;

    OVVCHdl ovvc_hdl;
    int ret = 0;

    while (1) {

        const static struct option long_options[] =
        {
            {"version",   no_argument,       0, 'v'},
            {"help",      no_argument,       0, 'h'},
            {"log-level", required_argument, 0, 'l'},
            {"infile",    required_argument, 0, 'i'},
            {"outfile",   required_argument, 0, 'o'},
            {"framethr",  required_argument, 0, 't'},
            {"entrythr",  required_argument, 0, 'e'},
            {"upscale",   required_argument, 0, 'u'},
        };

        int option_index = 0;

        c = getopt_long(argc, argv, "vhl:i:o:t:e:u:", long_options,
                        &option_index);
        if (c == -1){
            break;
        }

        switch (c)
        {
            case 'v':
                options_flag += 0x01;
                break;

            case 'h':
                options_flag += 0x10;
                break;

            case 'l':
                ov_log_level = optarg[0]-'0';
                break;

            case 'i':
                /*TODO: Sanitize filename*/
                input_file_name = optarg;
                break;

            case 'o':
                /*TODO: Sanitize filename*/
                output_file_name = optarg;
                break;

            case 'u':
                upscale_flag = atoi(optarg);
                break;

            case 't':
                nb_frame_th = atoi(optarg);
                break;

            case 'e':
                nb_entry_th = atoi(optarg);
                break;

            case '?':
                options_flag += 0x10;
                break;
            default:
                abort();
        }
    }

    if (OVLOG_ERROR <= ov_log_level && ov_log_level <= OVLOG_DEBUG){
        ovlog_set_log_level(ov_log_level);
    }

    if (input_file_name == NULL){
        input_file_name ="test.266";
    }

    if (output_file_name == NULL){
        output_file_name ="test.yuv";
    }

    if (options_flag){

        if (options_flag & 0x01) {
            print_version();
        }

        if (options_flag & 0x10) {
            print_usage();
        }

        return 0;
    }

    fout = fopen(output_file_name, "wb");

    if (fout == NULL) {
        ov_log(NULL, OVLOG_ERROR, "Failed to open output file '%s'.\n", output_file_name);
        goto failinit;
    } else {
        ov_log(NULL, OVLOG_INFO, "Decoded stream will be written to '%s'.\n", output_file_name);
    }

    ret = init_openvvc_hdl(&ovvc_hdl, output_file_name, nb_frame_th, nb_entry_th, upscale_flag);

    if (ret < 0) goto failinit;

    ret = dmx_attach_file(&ovvc_hdl, input_file_name);

    if (ret < 0) goto failattach;

    read_write_stream(&ovvc_hdl, fout);

    ovdmx_detach_stream(ovvc_hdl.dmx);

failattach:
    ret = close_openvvc_hdl(&ovvc_hdl);
    fclose(fout);

failinit:
    return ret;
}

static int
dmx_attach_file(OVVCHdl *const vvc_hdl, const char *const input_file_name)
{
    OVFileIO *file_io = ovio_new_fileio(input_file_name, "rb");
    int ret;

    if (file_io == NULL) {
        perror(input_file_name);
        vvc_hdl->io = NULL;
        return -1;
    }

    vvc_hdl->io = (OVIO*) file_io;

    ret = ovdmx_attach_stream(vvc_hdl->dmx, (OVIO*) vvc_hdl->io);

    return ret;
}

static int
init_openvvc_hdl(OVVCHdl *const ovvc_hdl, const char *output_file_name, int nb_frame_th, int nb_entry_th, int upscale_flag)
{
    OVVCDec **vvcdec = &ovvc_hdl->dec;
    OVVCDmx **vvcdmx = &ovvc_hdl->dmx;
    int ret;

    int display_output = !!strcmp(output_file_name, "/dev/null");

    ret = ovdec_init(vvcdec, display_output, nb_frame_th, nb_entry_th, upscale_flag);

    if (ret < 0) goto faildec;

    ov_log(vvcdec, OVLOG_TRACE, "Decoder init.\n");

    ret = ovdmx_init(vvcdmx);

    if (ret < 0) goto faildmx;

    ov_log(vvcdmx, OVLOG_TRACE, "Demuxer init.\n");

    return 0;

faildec:
    ov_log(NULL, OVLOG_ERROR, "Decoder failed at init.\n");
    return ret;

faildmx:
    ov_log(NULL, OVLOG_ERROR, "Demuxer failed at init.\n");
    ovdec_close(*vvcdec);
    return ret;
}

static int
close_openvvc_hdl(OVVCHdl *const ovvc_hdl)
{
    OVVCDec *vvcdec = ovvc_hdl->dec;
    OVVCDmx *vvcdmx = ovvc_hdl->dmx;
    OVIO* io = ovvc_hdl->io;
    int ret;

    if (io != NULL) {
        io->close(io);
    }

    ret = ovdec_close(vvcdec);

    if (ret < 0) goto faildecclose;

    ret = ovdmx_close(vvcdmx);

    if (ret < 0) goto faildmxclose;

    return 0;

faildecclose:
    /* Do not check for dmx failure since it might override
       return value to a correct one in either success or
       failure we already raised an error */
    ov_log(NULL, OVLOG_ERROR, "Decoder failed at cloture.\n");
    ovdmx_close(vvcdmx);

faildmxclose:
    return ret;
}

static int
read_write_stream(OVVCHdl *const hdl, FILE *fout)
{
    OVVCDec *const dec = hdl->dec;
    int nb_pic = 0;
    int ret;
    int max_frame_h[3] = {0};

    do {
        OVPictureUnit *pu = NULL;
        ret = ovdmx_extract_picture_unit(hdl->dmx, &pu);
        if (ret < 0) {
            break;
        }

        if (pu) {
            int nb_pic2;
            ret = ovdec_submit_picture_unit(dec, pu);

            ovpu_unref(&pu);

            do {
                OVFrame *frame = NULL;
                nb_pic2 = ovdec_receive_picture(dec, &frame);

                if (frame) {
                    write_decoded_frame_to_file(frame, fout, max_frame_h);

                    ov_log(NULL, OVLOG_DEBUG, "Got output picture with POC %d.\n", frame->poc);

                    ovframe_unref(&frame);

                    ++nb_pic;
                }
            } while (nb_pic2 > 0);

        } else {
            break;
        }


    } while (ret >= 0);

    ret = 1;

    while (ret > 0) {
        OVFrame *frame = NULL;
        ret = ovdec_drain_picture(dec, &frame);

        if (frame) {
            write_decoded_frame_to_file(frame, fout, max_frame_h);

            ov_log(NULL, OVLOG_DEBUG, "Drain picture with POC %d.\n", frame->poc);

            ovframe_unref(&frame);

            ++nb_pic;
        }
    }

    ov_log(NULL, OVLOG_INFO, "Decoded %d pictures\n", nb_pic);

    return 1;
}

static int
write_decoded_frame_to_file(OVFrame *const frame, FILE *fp, int* max_frame_h)
{
    uint8_t component = 0;
    int ret = 0;
    struct ScalingInfo scale_info = frame->scale_info;
    int bd_shift = (frame->frame_info.chroma_format == OV_YUV_420_P8) ? 0: 1;

    uint8_t * zeros = ov_mallocz(frame->linesize[0] * sizeof(uint8_t));
    for (component = 0; component < 3; component++) {
        uint16_t add_w = (scale_info.scaling_win_left + scale_info.scaling_win_right);
        uint16_t add_h = (scale_info.scaling_win_top  + scale_info.scaling_win_bottom);
        add_w = component ? add_w : add_w << 1; 
        add_h = component ? add_h : add_h << 1;
        uint16_t win_top  =  component ? scale_info.scaling_win_top  : scale_info.scaling_win_top  << 1;
        uint16_t win_left =  component ? scale_info.scaling_win_left : scale_info.scaling_win_left << 1;
        int frame_h = frame->height[component] - add_h;
        int frame_w = frame->width[component]  - add_w;
        max_frame_h[component] = OVMAX(max_frame_h[component], frame_h);

        for (int j = win_top; j < frame_h + win_top; j++){
            int offset_h = j * frame->linesize[component] ;
            int offset   = offset_h + (win_left << 1) ;
            ret += fwrite(&frame->data[component][offset], frame_w << bd_shift, sizeof(uint8_t), fp);
            ret += fwrite(zeros, frame->linesize[component] - (frame_w << bd_shift), sizeof(uint8_t), fp);
        }
        for (int j = frame_h ; j < max_frame_h[component]; j++){
            ret += fwrite(zeros, frame->linesize[component], sizeof(uint8_t), fp);
        }
    }
    ov_freep(&zeros);
    // for (component = 0; component < 3; component++) {
    //     max_frame_h[component] = OVMAX(max_frame_h[component], frame->height[component]);
    //     int frame_size = max_frame_h[component] * frame->linesize[component];
    //     ret += fwrite(frame->data[component], frame_size, sizeof(uint8_t), fp);
    // }

    return ret;
}

static void
print_version(){
    printf("libovvc version %s\n", ovdec_version());
}

static void print_usage(){
  printf("usage: dectest [options]\n");
  printf("options:\n");
  printf("\t-h, --help\t\t\t\tShow this message.\n");
  printf("\t-v, --version\t\t\t\tShow version information.\n");
  printf("\t-l <level>, --log-level=<level>\t\tDefine the level of verbosity. Value between 0 and 6. (Default: 2)\n");
  printf("\t-i <file>, --infile=<file>\t\tPath to the file to be decoded (Default: test.266).\n");
  printf("\t-o <file>, --outfile=<file>\t\tPath to the output file (Default: test.yuv).\n");
  printf("\t-t <nbthreads>, --framethr=<nbthreads>\t\tNumber of simultaneous frames decoded (Default: 0).\n");
  printf("\t-e <nbthreads>, --entrythr=<nbthreads>\t\tNumber of simultaneous entries decoded per frame (Default: 0).\n");
}
