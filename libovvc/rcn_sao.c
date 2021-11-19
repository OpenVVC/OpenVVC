#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "ovutils.h"
#include "ovdefs.h"
#include "ovframe.h"
#include "ovdpb.h"
#include "slicedec.h"

#include "bitdepth.h"

static void
sao_band_filter(uint8_t *_dst, uint8_t *_src,
                ptrdiff_t stride_dst, ptrdiff_t stride_src,
                SAOParamsCtu *sao,
                int width, int height,
                int c_idx)
{
    int16_t *dst = (int16_t *)_dst;
    int16_t *src = (int16_t *)_src;
    int offset_table[32] = { 0 };
    int k, y, x;
    int shift  = BITDEPTH - 5;

    int16_t *sao_offset_val = sao->offset_val[c_idx];
    uint8_t sao_left_class  = sao->band_position[c_idx];

    stride_src /= sizeof(int16_t);
    stride_dst /= sizeof(int16_t);

    for (k = 0; k < 4; k++) {
        offset_table[(k + sao_left_class) & 31] = sao_offset_val[k];
    }

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = ov_bdclip(src[x] + offset_table[src[x] >> shift]);
        }
        dst += stride_dst;
        src += stride_src;
    }
}

#define CMP(a, b) ((a) > (b) ? 1 : ((a) == (b) ? 0 : -1))
static void
sao_edge_filter(uint8_t *_dst, uint8_t *_src,
                ptrdiff_t stride_dst, ptrdiff_t stride_src,
                SAOParamsCtu *sao, int width, int height,
                int c_idx)
{
    // static const uint8_t edge_idx[] = { 1, 2, 0, 3, 4 };
    static const int8_t pos[4][2][2] = {
        { { -1,  0 }, {  1, 0 } }, // horizontal
        { {  0, -1 }, {  0, 1 } }, // vertical
        { { -1, -1 }, {  1, 1 } }, // 45 degree
        { {  1, -1 }, { -1, 1 } }, // 135 degree
    };

    int16_t *sao_offset_val = sao->offset_val[c_idx];
    uint8_t eo = sao->eo_class[c_idx];
    int16_t *dst = (int16_t *)_dst;
    int16_t *src = (int16_t *)_src;

    stride_src /= sizeof(int16_t);
    stride_dst /= sizeof(int16_t);
    int a_stride, b_stride;
    int src_offset = 0;
    int dst_offset = 0;
    int x, y;

    a_stride = pos[eo][0][0] + pos[eo][0][1] * stride_src;
    b_stride = pos[eo][1][0] + pos[eo][1][1] * stride_src;
    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            int diff0         = CMP(src[x + src_offset], src[x + src_offset + a_stride]);
            int diff1         = CMP(src[x + src_offset], src[x + src_offset + b_stride]);
            int offset_val    = 2 + diff0 + diff1;
            dst[x + dst_offset] = ov_bdclip(src[x + src_offset] + sao_offset_val[offset_val]);
        }
        src_offset += stride_src;
        dst_offset += stride_dst;
    }
}

static void
rcn_sao_ctu(OVCTUDec *const ctudec, SAOParamsCtu *sao, int x_start_pic, int y_start_pic, int y_end_pic, int fb_offset, uint8_t is_border)
{   
    struct OVFilterBuffers* fb   = &ctudec->filter_buffers;
    const OVPartInfo *const pinfo = ctudec->part_ctx;
    uint8_t log2_ctb_s = pinfo->log2_ctu_s;

    const struct SAOFunctions *saofunc = &ctudec->rcn_ctx.rcn_funcs.sao;
    OVFrame *frame = fb->pic_frame;

    for (int c_idx = 0; c_idx < (ctudec->sao_info.chroma_format_idc ? 3 : 1); c_idx++) {
        int shift_chr = c_idx == 0 ? 0 : 1;
        int x0       = x_start_pic >> shift_chr;
        int y0       = y_start_pic >> shift_chr;
        int f_width  = (ctudec->pic_w) >> shift_chr;
        int f_height = (ctudec->pic_h) >> shift_chr;

        int ctb_size_h = (1 << log2_ctb_s) >> shift_chr;
        int ctb_size_v = (y_end_pic - y_start_pic) >> shift_chr;
        int width    = OVMIN(ctb_size_h, f_width - x0);
        int height   = OVMIN(ctb_size_v, f_height - y0);

        int int16_t_shift = 1;
        ptrdiff_t stride_out_pic = frame->linesize[c_idx];
        uint8_t *out_pic = frame->data[c_idx];
        out_pic = &out_pic[ y0 * stride_out_pic + (x0<<int16_t_shift)];

        uint8_t *filtered = (uint8_t *) fb->filter_region[c_idx];
        int stride_filtered = fb->filter_region_stride[c_idx]<<int16_t_shift;
        filtered = &filtered[(fb->filter_region_offset[c_idx]<<int16_t_shift) + stride_filtered*(fb_offset>> shift_chr)];  

        switch (sao->type_idx[c_idx]) {
            case SAO_BAND:
                saofunc->band(out_pic, filtered, stride_out_pic, stride_filtered,
                               sao, width, height, c_idx);

                break;

            case SAO_EDGE:
            {

                int x_start = 0;
                int y_start = 0;
                if ((is_border & OV_BOUNDARY_LEFT_RECT) && sao->eo_class[c_idx] != 1){
                    x_start = 1;
                    width   = width-1;
                }
                if ((is_border & OV_BOUNDARY_UPPER_RECT) && sao->eo_class[c_idx] != 0){
                    y_start = 1;
                    height  = height-1;
                }
                if ((is_border & OV_BOUNDARY_RIGHT_RECT) && sao->eo_class[c_idx] != 1){
                    width   = width-1;
                }
                if ((is_border & OV_BOUNDARY_BOTTOM_RECT) && sao->eo_class[c_idx] != 0){
                    height  = height-1;
                }

                int src_offset = y_start*stride_out_pic + x_start*sizeof(uint16_t);
                int dst_offset = y_start*stride_filtered + x_start*sizeof(uint16_t);

                saofunc->edge[!(width % 8)](out_pic + src_offset, filtered + dst_offset, stride_out_pic, stride_filtered, sao, width, height, c_idx);

                break;
            }
        }
    }
}

void
rcn_sao_filter_line(OVCTUDec *const ctudec, const struct RectEntryInfo *const einfo, uint16_t ctb_y)
{
    if (!ctudec->sao_info.sao_luma_flag && !ctudec->sao_info.sao_chroma_flag){
        return;
    }

    const OVPartInfo *const pinfo = ctudec->part_ctx;
    uint8_t log2_ctb_s = pinfo->log2_ctu_s;
    int ctu_width  = 1 << log2_ctb_s;

    struct OVFilterBuffers* fb = &ctudec->filter_buffers;
    int margin = 2*fb->margin;

    for (int ctb_x = 0; ctb_x < einfo->nb_ctu_w; ctb_x++) {
        int ctb_x_pic = ctb_x + einfo->ctb_x;
        int ctb_y_pic = ctb_y + einfo->ctb_y;
        int x_pos = ctu_width * ctb_x;
        int x_pos_pic = ctu_width * ctb_x_pic;
        int y_pos_pic = ctu_width * ctb_y_pic;
        
        uint8_t is_border = 0;
        is_border = (ctb_x == 0)                   ? is_border | OV_BOUNDARY_LEFT_RECT: is_border;
        is_border = (ctb_x == einfo->nb_ctu_w - 1) ? is_border | OV_BOUNDARY_RIGHT_RECT: is_border;
        is_border = (ctb_y == einfo->nb_ctu_h - 1) ? is_border | OV_BOUNDARY_BOTTOM_RECT: is_border;

        int y_start_pic = y_pos_pic + margin;

        ctudec_extend_filter_region(ctudec, fb->saved_rows_sao, x_pos, x_pos_pic, y_start_pic, is_border);

        int fb_offset = 0;
        int ctb_addr_rs = ctb_y * einfo->nb_ctu_w + ctb_x;
        SAOParamsCtu *sao  = &ctudec->sao_info.sao_params[ctb_addr_rs];
        rcn_sao_ctu(ctudec, sao, x_pos_pic, y_start_pic, y_pos_pic + ctu_width, fb_offset, is_border);

        if (!(is_border & OV_BOUNDARY_BOTTOM_RECT)) {
            fb_offset = ctu_width - margin;
            y_pos_pic += ctu_width;
            ctb_addr_rs += einfo->nb_ctu_w;
            sao         =  &ctudec->sao_info.sao_params[ctb_addr_rs];
            rcn_sao_ctu(ctudec, sao, x_pos_pic, y_pos_pic, y_pos_pic + margin, fb_offset, is_border);
        }

        ctudec_save_last_rows(ctudec, fb->saved_rows_sao, x_pos, x_pos_pic, y_start_pic, is_border);
        ctudec_save_last_cols(ctudec, x_pos_pic, y_start_pic, is_border);
    }
}

void
rcn_sao_first_pix_rows(OVCTUDec *const ctudec, const struct RectEntryInfo *const einfo, uint16_t ctb_y)
{
    if (!ctudec->sao_info.sao_luma_flag && !ctudec->sao_info.sao_chroma_flag){
        return;
    }

    const OVPartInfo *const pinfo = ctudec->part_ctx;
    uint8_t log2_ctb_s = pinfo->log2_ctu_s;
    int ctu_width  = 1 << log2_ctb_s;
    int ctb_y_pic = ctb_y + einfo->ctb_y;

    struct OVFilterBuffers* fb = &ctudec->filter_buffers;
    int margin = 2 * fb->margin;

    for (int ctb_x = 0; ctb_x < einfo->nb_ctu_w; ctb_x++) {
        int ctb_x_pic = ctb_x + einfo->ctb_x;
        int x_pos = ctu_width * ctb_x;
        int x_pos_pic = ctu_width * ctb_x_pic;
        
        //left | right | up | down
        uint8_t is_border = 0; 
        is_border = (ctb_x==0)          ? is_border | OV_BOUNDARY_LEFT_RECT: is_border;
        is_border = (ctb_x == einfo->nb_ctu_w - 1) ? is_border | OV_BOUNDARY_RIGHT_RECT: is_border;
        is_border = (ctb_y==0)          ? is_border | OV_BOUNDARY_UPPER_RECT: is_border;

        //Apply SAO of previous ctu line
        int y_start_pic = ctu_width * ctb_y_pic;
        ctudec_extend_filter_region(ctudec, fb->saved_rows_sao, x_pos, x_pos_pic, y_start_pic, is_border);

        int fb_offset = 0;
        int ctb_addr_rs    = ctb_y * einfo->nb_ctu_w + ctb_x;
        SAOParamsCtu *sao  = &ctudec->sao_info.sao_params[ctb_addr_rs];
        rcn_sao_ctu(ctudec, sao, x_pos_pic, y_start_pic, y_start_pic + margin, fb_offset, is_border);

        const int width_l = fb->filter_region_w[0];
        for (int comp = 0; comp < 3; comp++) {
            int16_t* saved_rows = fb->saved_rows_sao[comp];
            int16_t* filter_region = fb->filter_region[comp];
            int stride_filter = fb->filter_region_stride[comp];

            int ratio_luma_chroma = 2;
            int ratio = comp==0 ? 1 : ratio_luma_chroma;        
            const int width = width_l/ratio;
            const int x = x_pos/ratio;
            int stride_rows = fb->saved_rows_stride[comp];
            int offset_y = comp==0 ? 2 * fb->margin : fb->margin;
   
            for (int ii=0; ii < fb->margin; ii++) {
                memcpy(&saved_rows[ii*stride_rows + x], &filter_region[(offset_y+ii)*stride_filter + fb->margin], width * sizeof(int16_t));
            }
        }

        ctudec_save_last_cols(ctudec, x_pos_pic, y_start_pic, is_border);
    }
}

void
BD_DECL(rcn_init_sao_functions)(struct RCNFunctions *const rcn_funcs)
{
    rcn_funcs->sao.band= &sao_band_filter;
    rcn_funcs->sao.edge[0]= &sao_edge_filter;
    rcn_funcs->sao.edge[1]= &sao_edge_filter;
}
