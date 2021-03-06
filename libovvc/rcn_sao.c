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

void sao_band_filter(uint8_t *_dst, uint8_t *_src,
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

void sao_edge_filter(uint8_t *_dst, uint8_t *_src,
        ptrdiff_t stride_dst, ptrdiff_t stride_src,
        SAOParamsCtu *sao,
         int width, int height, int x_start, int y_start, 
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
    int src_offset = y_start * stride_src;
    int dst_offset = y_start * stride_dst;
    int x, y;

    a_stride = pos[eo][0][0] + pos[eo][0][1] * stride_src;
    b_stride = pos[eo][1][0] + pos[eo][1][1] * stride_src;
    for (y = y_start; y < y_start + height; y++) {
        for (x = x_start; x < x_start + width; x++) {
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


void rcn_sao_ctu(OVCTUDec *const ctudec, int ctb_x_pic, int ctb_y_pic, int nb_ctu_w, uint8_t is_border) 
{   
    struct OVFilterBuffers fb   = ctudec->filter_buffers;
    const OVPartInfo *const pinfo = ctudec->part_ctx;
    uint8_t log2_ctb_size = pinfo->log2_ctu_s;

    int x           = ctb_x_pic << log2_ctb_size;
    int y           = ctb_y_pic << log2_ctb_size;
    int ctb_addr_rs = ctb_y_pic * nb_ctu_w + ctb_x_pic;
    SAOParamsCtu *sao  = &ctudec->sao_info.sao_params[ctb_addr_rs];

    OVFrame *frame = fb.pic_frame;
    // struct OVBuffInfo frame_buff = ctudec->rcn_ctx.frame_buff;

    for (int c_idx = 0; c_idx < (ctudec->sao_info.chroma_format_idc ? 3 : 1); c_idx++) {
        int shift_chr = c_idx==0 ? 0 : 1;
        int x0       = x >> shift_chr;
        int y0       = y >> shift_chr;
        int f_width  = (ctudec->pic_w) >> shift_chr;
        int f_height = (ctudec->pic_h) >> shift_chr;
        int ctb_size_h = (1 << log2_ctb_size) >> shift_chr;
        int ctb_size_v = (1 << log2_ctb_size) >> shift_chr;
        int width    = OVMIN(ctb_size_h, f_width - x0);
        int height   = OVMIN(ctb_size_v, f_height - y0);

        //BITDEPTH: uniquement pour bitdepth 10
        int int16_t_shift = 1;
        ptrdiff_t stride_out_pic = frame->linesize[c_idx];
        uint8_t *out_pic = frame->data[c_idx];
        // ptrdiff_t stride_out_pic = c_idx==0 ? frame_buff.stride : frame_buff.stride_c;
        // uint8_t *out_pic = frame_buff.y;
        // if (c_idx != 0) 
        //     c_idx==1 ? out_pic = frame_buff.cb : frame_buff.cr;
        out_pic = &out_pic[ y0 * stride_out_pic + (x0<<int16_t_shift)];
 
        uint8_t *filtered = (uint8_t *) fb.filter_region[c_idx];
        int stride_filtered = fb.filter_region_stride[c_idx]<<int16_t_shift;
        filtered = &filtered[(fb.filter_region_offset[c_idx]<<int16_t_shift)];  
        // filtered = &filtered[ (y0 * stride_filtered) + (x0<<int16_t_shift) + (fb.filter_region_offset[c_idx]<<int16_t_shift)];  

        switch (sao->type_idx[c_idx]) {
            case SAO_BAND:
                // printf("BO %i %i %i\n", c_idx, x, y);

                sao_band_filter(out_pic, filtered, stride_out_pic, stride_filtered,
                               sao, width, height, c_idx);

                sao->type_idx[c_idx] = SAO_APPLIED;
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
                sao_edge_filter(out_pic, filtered, stride_out_pic, stride_filtered, sao, width, height, x_start, y_start, c_idx);

                sao->type_idx[c_idx] = SAO_APPLIED;
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

    for (int ctb_x = 0; ctb_x < nb_ctu_w; ctb_x++) 
    {
        // int ctb_x_pic = tile_ctx->ctu_x[tile_x] + ctb_x;
        int ctb_x_pic   = ctb_x;
        int xPos = ctu_width * ctb_x_pic;
        int yPos = ctu_width * ctb_y_pic;
        
        //left | right | up | down
        uint8_t is_border = 0; 
        is_border = (ctb_x==0)          ? is_border | OV_BOUNDARY_LEFT_RECT: is_border;
        is_border = (ctb_x==nb_ctu_w-1) ? is_border | OV_BOUNDARY_RIGHT_RECT: is_border;
        is_border = (ctb_y_pic==0)          ? is_border | OV_BOUNDARY_UPPER_RECT: is_border;
        // is_border = (ctb_y==nb_ctu_h-1) ? is_border | OV_BOUNDARY_BOTTOM_RECT: is_border;
        is_border = (yPos + ctu_width >= ctudec->pic_h) ? is_border | OV_BOUNDARY_BOTTOM_RECT: is_border;


        ctudec_extend_filter_region(ctudec, xPos, yPos, is_border);

        rcn_sao_ctu(ctudec, ctb_x_pic, ctb_y_pic, nb_ctu_w, is_border);

        ctudec_save_last_rows(ctudec, xPos, yPos, is_border);
        ctudec_save_last_cols(ctudec, xPos, yPos, is_border);
    }
}
