#include <stdint.h>
#include <stdlib.h>

#include "ovutils.h"
#include "ovdefs.h"
#include "ovframe.h"
#include "ovdpb.h"
#include "slicedec.h"

// #include "dec_structures.h"
// #include "ctudec.h"

#define BIT_DEPTH 10

static void sao_band_filter(uint8_t *_dst, uint8_t *_src,
        ptrdiff_t stride_dst, ptrdiff_t stride_src,
        SAOParamsCtu *sao,
         int width, int height,
        int c_idx)
{
    int16_t *dst = (int16_t *)_dst;
    int16_t *src = (int16_t *)_src;
    int offset_table[32] = { 0 };
    int k, y, x;
    //BITDEPTH: uniquement pour bitdepth 10
    int shift  = BIT_DEPTH - 5;

    int16_t *sao_offset_val = sao->offset_val[c_idx];
    uint8_t sao_left_class  = sao->band_position[c_idx];

    stride_src /= sizeof(int16_t);
    stride_dst /= sizeof(int16_t);

    for (k = 0; k < 4; k++)
        offset_table[(k + sao_left_class) & 31] = sao_offset_val[k];
        // offset_table[(k + sao_left_class) & 31] = sao_offset_val[k + 1];
    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++)
            #if 0
            dst[x] = ov_clip_uintp2(src[x] + offset_table[src[x] >> shift], BIT_DEPTH);
            #else
            dst[x] = ov_clip(src[x] + offset_table[src[x] >> shift], 0, 1023);
            #endif
        dst += stride_dst;
        src += stride_src;
    }
}

////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
#define CMP(a, b) ((a) > (b) ? 1 : ((a) == (b) ? 0 : -1))

static void sao_edge_filter(uint8_t *_dst, uint8_t *_src,
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
            #if 0
            dst[x + dst_offset] = ov_clip_uintp2( src[x + src_offset] + sao_offset_val[offset_val], BIT_DEPTH );
            #else
            dst[x + dst_offset] = ov_clip( src[x + src_offset] + sao_offset_val[offset_val], 0, 1023 );
            #endif
        }
        src_offset += stride_src;
        dst_offset += stride_dst;
    }
}


void rcn_sao_ctu(OVCTUDec *const ctudec, int ctb_x_pic, int ctb_y_pic, int start_h, int end_h, int fb_offset, int nb_ctu_w, uint8_t is_border) 
{   
    struct OVFilterBuffers fb   = ctudec->filter_buffers;
    const OVPartInfo *const pinfo = ctudec->part_ctx;
    uint8_t log2_ctb_size = pinfo->log2_ctu_s;

    int x           = ctb_x_pic << log2_ctb_size;
    // int y           = ctb_y_pic << log2_ctb_size;
    int y           = start_h;
    int ctb_addr_rs = ctb_y_pic * nb_ctu_w + ctb_x_pic;
    SAOParamsCtu *sao  = &ctudec->sao_info.sao_params[ctb_addr_rs];
    const struct SAOFunctions *saofunc = &ctudec->rcn_ctx.rcn_funcs.sao;

    OVFrame *frame = fb.pic_frame;

    for (int c_idx = 0; c_idx < (ctudec->sao_info.chroma_format_idc ? 3 : 1); c_idx++) {
        int shift_chr = c_idx==0 ? 0 : 1;
        int x0       = x >> shift_chr;
        int y0       = y >> shift_chr;
        int f_width  = (ctudec->pic_w) >> shift_chr;
        int f_height = (ctudec->pic_h) >> shift_chr;
        int ctb_size_h = (1 << log2_ctb_size) >> shift_chr;
        // int ctb_size_v = (1 << log2_ctb_size) >> shift_chr;
        int ctb_size_v = (end_h - start_h) >> shift_chr;
        int width    = OVMIN(ctb_size_h, f_width - x0);
        int height   = OVMIN(ctb_size_v, f_height - y0);

        //BITDEPTH: uniquement pour bitdepth 10
        int int16_t_shift = 1;
        ptrdiff_t stride_out_pic = frame->linesize[c_idx];
        uint8_t *out_pic = frame->data[c_idx];
        out_pic = &out_pic[ y0 * stride_out_pic + (x0<<int16_t_shift)];

        uint8_t *filtered = (uint8_t *) fb.filter_region[c_idx];
        int stride_filtered = fb.filter_region_stride[c_idx]<<int16_t_shift;
        filtered = &filtered[(fb.filter_region_offset[c_idx]<<int16_t_shift) + stride_filtered*(fb_offset>> shift_chr)];  

        switch (sao->type_idx[c_idx]) {
            case SAO_BAND:
                // printf("BO %i %i %i\n", c_idx, x, y);
                saofunc->band(out_pic, filtered, stride_out_pic, stride_filtered,
                               sao, width, height, c_idx);

                break;

            case SAO_EDGE:
            {
                // printf("EO %i %i %i\n", c_idx, x, y);

                //Do not apply filters on image borders, compliant with VTM
                int x_start = 0;
                int y_start = 0;
                if ( (is_border & OV_BOUNDARY_LEFT_RECT) && sao->eo_class[c_idx] != 1){
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

                //parameters: buffer
                int src_offset = y_start*stride_out_pic + x_start*sizeof(uint16_t);
                int dst_offset = y_start*stride_filtered + x_start*sizeof(uint16_t);

                saofunc->edge[!(width % 8)](out_pic + src_offset, filtered + dst_offset, stride_out_pic, stride_filtered, sao, width, height, c_idx);

                break;
            }
        }
    }
}

void rcn_sao_filter_line(OVCTUDec *const ctudec, int nb_ctu_w, uint16_t ctb_y_pic)
{
    if (!ctudec->sao_info.sao_luma_flag && !ctudec->sao_info.sao_chroma_flag){
        return;
    }
    const OVPartInfo *const pinfo = ctudec->part_ctx;
    uint8_t log2_ctb_size = pinfo->log2_ctu_s;
    int ctu_width  = 1 << log2_ctb_size;

    struct OVFilterBuffers* fb = &ctudec->filter_buffers;
    int margin = 2*fb->margin;

    for (int ctb_x = 0; ctb_x < nb_ctu_w; ctb_x++) 
    {
        // int ctb_x_pic = tile_ctx->ctu_x[tile_x] + ctb_x;
        int ctb_x_pic   = ctb_x;
        int x_pos_ctu = ctu_width * ctb_x_pic;
        int y_pos_ctu = ctu_width * ctb_y_pic;
        
        //left | right | up | down
        uint8_t is_border = 0;
        is_border = (ctb_x==0)          ? is_border | OV_BOUNDARY_LEFT_RECT: is_border;
        is_border = (ctb_x==nb_ctu_w-1) ? is_border | OV_BOUNDARY_RIGHT_RECT: is_border;
        is_border = (ctb_y_pic==0)          ? is_border | OV_BOUNDARY_UPPER_RECT: is_border;
        // is_border = (ctb_y==nb_ctu_h-1) ? is_border | OV_BOUNDARY_BOTTOM_RECT: is_border;
        is_border = (y_pos_ctu + ctu_width >= ctudec->pic_h) ? is_border | OV_BOUNDARY_BOTTOM_RECT: is_border;

        //TODOfilt: leave y_start_h negative and treat in ctudec extend and 
        //ctudec_save_last_rows and ctudec_save_last_cols

        //Apply SAO of previous ctu line
        int y_start = y_pos_ctu + margin;
        ctudec_extend_filter_region(ctudec, x_pos_ctu, y_start, is_border);

        int fb_offset = 0;
        rcn_sao_ctu(ctudec, ctb_x_pic, ctb_y_pic, y_start, y_pos_ctu + ctu_width, fb_offset, nb_ctu_w, is_border);
        if ( ! (is_border & OV_BOUNDARY_BOTTOM_RECT)){
            fb_offset = ctu_width - margin;
            y_pos_ctu += ctu_width;
            rcn_sao_ctu(ctudec, ctb_x_pic, ctb_y_pic+1, y_pos_ctu, y_pos_ctu + margin, fb_offset, nb_ctu_w, is_border);
        }

        ctudec_save_last_rows(ctudec, x_pos_ctu, y_start, is_border);
        ctudec_save_last_cols(ctudec, x_pos_ctu, y_start, is_border);
    }
}

void rcn_init_sao_functions(struct RCNFunctions *const rcn_funcs){
    rcn_funcs->sao.band= &sao_band_filter;
    rcn_funcs->sao.edge[0]= &sao_edge_filter;
    rcn_funcs->sao.edge[1]= &sao_edge_filter;
}
