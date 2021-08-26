#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>

#include "ovdec.h"
#include "ovdefs.h"
#include "ovdmx.h"
#include "ovframe.h"
#include "ovdpb.h"
#include "ovutils.h"
#include "ovversion.h"

typedef struct OVVCHdl{
    OVVCDmx *dmx;
    OVVCDec *dec;
    /* TODO decide whether or not file pointer must be given
       to dmx or not  it is only given to the hadle so when can
       close the file later since we opened it from here*/
    FILE *fp;
}OVVCHdl;

static int dmx_attach_file(OVVCHdl *const vvc_hdl, const char *const input_file_name);

static int init_openvvc_hdl(OVVCHdl *const ovvc_hdl, const char *output_file_name, int nb_threads);

static int close_openvvc_hdl(OVVCHdl *const ovvc_hdl);

static int read_write_stream(OVVCHdl *const hdl, FILE *fp, FILE *fout);
static int read_stream(OVVCHdl *const hdl, FILE *fp);

static uint32_t write_decoded_frame_to_file(OVFrame *const frame, FILE *fp);

static void print_version(void);

static void print_usage(void);


int
main(int argc, char** argv)
{
  /* basic options parser and assign
     filenames into a functions*/
  int c;
  const char *input_file_name = NULL, *output_file_name = NULL;
  int ov_log_level=OVLOG_INFO;
  FILE *fout = NULL;
  int nb_threads = 1;

  uint8_t options_flag=0;

  OVVCHdl ovvc_hdl;
  int ret = 0;
  while (1)
    {
      static struct option long_options[] =
        {
          {"version", no_argument,      0, 'v'},
          {"help",    no_argument,       0, 'h'},
          {"log-level", required_argument, 0, 'l'},
          {"infile",      required_argument, 0, 'i'},
          {"outfile",      required_argument, 0, 'o'},
          {"threads",      required_argument, 0, 't'},
        };
      int option_index = 0;

      c = getopt_long (argc, argv, "vhl:i:o:t:",
                       long_options, &option_index);
      if (c == -1){
        break;
      }
      switch (c)
        {
        case 'v':
          options_flag+=0x01;
          break;

        case 'h':
          options_flag+=0x10;
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

        case 't':
          nb_threads = atoi(optarg);
          break;

        case '?':
          options_flag+=0x10;
          break;
        default:
          abort ();
        }
    }

    if (OVLOG_ERROR <= ov_log_level && ov_log_level <= OVLOG_TRACE){
      set_ov_log_level(ov_log_level);
    }

    if (input_file_name == NULL){
      input_file_name ="test.266";
    }

    if (output_file_name == NULL){
      output_file_name ="test.yuv";
    }
    if (output_file_name == NULL){
      output_file_name ="test.yuv";
    }
    if (options_flag){
      if (options_flag & 0x01) {print_version();}
      if (options_flag & 0x10) {print_usage();}
      return 0;
    }

    fout = fopen(output_file_name, "wb");
    if (fout == NULL) {
      ov_log(NULL, OVLOG_ERROR, "Failed to open output file '%s'.\n", output_file_name);
      goto failinit;
    } else {
      ov_log(NULL, OVLOG_INFO, "Decoded stream will be written to '%s'.\n", output_file_name);
    }
    ret = init_openvvc_hdl(&ovvc_hdl, output_file_name, nb_threads);

    if (ret < 0) goto failinit;

    ret = dmx_attach_file(&ovvc_hdl, input_file_name);

    if (ret < 0) goto failattach;

    read_write_stream(&ovvc_hdl, ovvc_hdl.fp, fout);
    // read_stream(&ovvc_hdl, ovvc_hdl.fp);

    ovdmx_detach_stream(ovvc_hdl.dmx);


    /* Do stuff here */

failattach:
    ret = close_openvvc_hdl(&ovvc_hdl);
    fclose(fout);

failinit:
    return ret;
}

static int
dmx_attach_file(OVVCHdl *const vvc_hdl, const char *const input_file_name)
{
    int ret;
    FILE *file = fopen(input_file_name,"rb");

    if (file == NULL) {
        perror(input_file_name);
       vvc_hdl->fp = NULL;
       return -1;
    }

    vvc_hdl->fp = file;

    ret = ovdmx_attach_stream(vvc_hdl->dmx, vvc_hdl->fp);

    return ret;
}

static int
init_openvvc_hdl(OVVCHdl *const ovvc_hdl, const char *output_file_name, int nb_threads)
{
    OVVCDec **vvcdec = &ovvc_hdl->dec;
    OVVCDmx **vvcdmx = &ovvc_hdl->dmx;
    int ret;

    ret = ovdec_init(vvcdec, output_file_name, nb_threads);

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
    ov_log(NULL, OVLOG_ERROR, "Decoder failed at cloture.\n");
    ovdmx_close(vvcdmx);

faildmxclose:
    return ret;
}


static int
read_stream(OVVCHdl *const hdl, FILE *fp)
{
    int ret;
    OVVCDmx *const dmx = hdl->dmx;
    OVVCDec *const dec = hdl->dec;
    OVPictureUnit *pu = NULL;
    
    do {
        ret = ovdmx_extract_picture_unit(dmx, &pu);
        if (ret < 0) {
            ov_log(NULL, OVLOG_INFO, "Picture unit not extracted\n");
            goto end_out;
        }

        if (pu){
            ret = ovdec_submit_picture_unit(dec, pu);
            if (ret < 0) {
                ov_log(NULL, OVLOG_INFO, "Picture unit not submitted\n");
                ov_free_pu(&pu);
                goto end_out;
            }
            /* FIXME Picture unit freeing be inside the decoder
             * use ref_counted buffer and call unref here instead
             */
            ov_free_pu(&pu);
        }
    } while (ret >= 0);

//TODOpar: stop cleanly program in case of error
end_out:
    return 1;
}

static int
read_write_stream(OVVCHdl *const hdl, FILE *fp, FILE *fout)
{
    #if 0
    ovdmx_read_stream(dmx);
    #else
    int ret;
    OVVCDmx *const dmx = hdl->dmx;
    OVVCDec *const dec = hdl->dec;
    OVPictureUnit *pu = NULL;
    OVFrame *frame = NULL;

    int nb_pic = 0;
    do {
        ret = ovdmx_extract_picture_unit(dmx, &pu);
        if (ret < 0) {
            break;
        }

        if (pu){
            ret = ovdec_submit_picture_unit(dec, pu);
            if (ret < 0) {
                ov_free_pu(&pu);
                break;
            }

            do {
                frame = NULL;
                ovdec_receive_picture(dec, &frame);

                /* FIXME use ret instead of frame */
                if (frame) {
                    write_decoded_frame_to_file(frame, fout);
                    ++nb_pic;
                    ov_log(NULL, OVLOG_DEBUG, "Got ouput picture with POC %d.\n", frame->poc);

                    ovframe_unref(&frame);
                }
            } while (frame);

            /* FIXME Picture unit freeing be inside the decoder
             * use ref_counted buffer and call unref here instead
             */
            ov_free_pu(&pu);
        }


    } while (ret >= 0);

    //Wait for all the sub decoders to finish before draining dpb
    ovdec_uninit_subdec_list(dec);

    ret = 1;
    while (ret > 0) {
        frame = NULL;
        ret = ovdec_drain_picture(dec, &frame);
        /* FIXME use ret instead of frame */
        if (frame) {
            write_decoded_frame_to_file(frame, fout);
            ++nb_pic;
            ov_log(NULL, OVLOG_DEBUG, "Drain picture with POC %d.\n", frame->poc);

            ovframe_unref(&frame);
        }
    }

    ov_log(NULL, OVLOG_INFO, "Decoded %d pictures\n", nb_pic);

    #endif
    return 1;
}


static uint32_t write_decoded_frame_to_file(OVFrame *const frame, FILE *fp){
  uint8_t component = 0;
  uint32_t ret = 0;
  for(component=0; component<3; component++){
    uint32_t frame_size = frame->height[component] * frame->linesize[component];
    ret +=fwrite(frame->data[component], frame_size, sizeof(uint8_t), fp);
  }
  return ret;
}



static void print_version(){
  print_ov_lib_version();
}

static void print_usage(){
  printf("usage: dectest [options]\n");
  printf("options:\n");
  printf("\t-h, --help\t\t\t\tShow this message.\n");
  printf("\t-v, --version\t\t\t\tShow version information.\n");
  printf("\t-l <level>, --log-level=<level>\t\tDefine the level of verbosity. Value between 0 and 6. (Default: 2)\n");
  printf("\t-i <file>, --infile=<file>\t\tPath to the file to be decoded (Default: test.266).\n");
  printf("\t-o <file>, --outfile=<file>\t\tPath to the output file (Default: test.yuv).\n");
  printf("\t-t <nbthreads>, --threads=<nbthreads>\t\tNumber of decoding threads (Default: 1).\n");
}
