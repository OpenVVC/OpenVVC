#include <stdint.h>
#include <string.h>

#include "ctudec.h"

void
rcn_ctu_copy_left_border(struct OVRCNCtx *rcn_ctx, uint8_t log2_ctb_s)
{
     uint8_t ctu_s = 1 << log2_ctb_s;
     struct OVBuffInfo *binfo = &rcn_ctx->ctu_buff;
     const OVSample *src = binfo->y + ctu_s - 4;
     OVSample *dst = binfo->y - 4;
     int i;

     /* Note we also copy border of above line for top left sample
      */
     src -= binfo->stride_c;
     dst -= binfo->stride_c;

     /* Copy 4 left samples because of MRL */
     for (i = 0; i < ctu_s + 1; ++i) {
         memcpy(dst, src, sizeof(*dst) * 4);
         dst += binfo->stride;
         src += binfo->stride;
     }

     src = binfo->cb + (ctu_s >> 1) - 1;
     dst = binfo->cb - 1;
     src -= binfo->stride_c;
     dst -= binfo->stride_c;

     for (i = 0; i < (ctu_s >> 1) + 1; ++i) {
         memcpy(dst, src, sizeof(*dst));
         dst += binfo->stride_c;
         src += binfo->stride_c;
     }

     src = binfo->cr + (ctu_s >> 1) - 1;
     dst = binfo->cr - 1;
     src -= binfo->stride_c;
     dst -= binfo->stride_c;

     for (i = 0; i < (ctu_s >> 1) + 1; ++i) {
         memcpy(dst, src, sizeof(*dst));
         dst += binfo->stride_c;
         src += binfo->stride_c;
     }
}

void
rcn_update_ctu_border(struct OVRCNCtx *rcn_ctx, uint8_t log2_ctb_s)
{
     uint8_t ctu_s = 1 << log2_ctb_s;
     struct OVBuffInfo *binfo = &rcn_ctx->ctu_buff;
     const OVSample *src = binfo->y + ctu_s - 4;
     OVSample *dst = binfo->y - 4;
     int i;

     struct OVBuffInfo *finfo = &rcn_ctx->frame_buff;
     finfo->y += ctu_s;
     finfo->cb += ctu_s >> 1;
     finfo->cr += ctu_s >> 1;
     /* Note we also copy border of above line for top left sample
      */
     src -= binfo->stride_c;
     dst -= binfo->stride_c;

     /* Copy 4 left samples because of MRL */
     for (i = 0; i < ctu_s + 1; ++i) {
         memcpy(dst, src, sizeof(*dst) * 4);
         dst += binfo->stride;
         src += binfo->stride;
     }

     src = binfo->cb + (ctu_s >> 1) - 1;
     dst = binfo->cb - 1;
     src -= binfo->stride_c;
     dst -= binfo->stride_c;

     for (i = 0; i < (ctu_s >> 1) + 1; ++i) {
         memcpy(dst, src, sizeof(*dst));
         dst += binfo->stride_c;
         src += binfo->stride_c;
     }

     src = binfo->cr + (ctu_s >> 1) - 1;
     dst = binfo->cr - 1;
     src -= binfo->stride_c;
     dst -= binfo->stride_c;

     for (i = 0; i < (ctu_s >> 1) + 1; ++i) {
         memcpy(dst, src, sizeof(*dst));
         dst += binfo->stride_c;
         src += binfo->stride_c;
     }
}

void
rcn_write_ctu_to_frame(const struct OVRCNCtx *const rcn_ctx, uint8_t log2_ctb_s)
{
    int i;
    const struct OVBuffInfo *const fd = &rcn_ctx->frame_buff;
    const OVSample *src_y  = &rcn_ctx->ctu_buff.y [0];
    const OVSample *src_cb = &rcn_ctx->ctu_buff.cb[0];
    const OVSample *src_cr = &rcn_ctx->ctu_buff.cr[0];

    OVSample *dst_y  = fd->y;
    OVSample *dst_cb = fd->cb;
    OVSample *dst_cr = fd->cr;

    for (i = 0; i < 1 << log2_ctb_s; i++) {
        memcpy(dst_y, src_y, sizeof(OVSample) << log2_ctb_s);
        src_y += RCN_CTB_STRIDE;
        dst_y += fd->stride;
    }

    for (i = 0; i < 1 << (log2_ctb_s - 1); ++i) {
        memcpy(dst_cb, src_cb, sizeof(OVSample) << (log2_ctb_s - 1));
        memcpy(dst_cr, src_cr, sizeof(OVSample) << (log2_ctb_s - 1));
        src_cb += RCN_CTB_STRIDE;
        src_cr += RCN_CTB_STRIDE;
        dst_cb += fd->stride_c;
        dst_cr += fd->stride_c;
    }
}

void
rcn_frame_line_to_ctu(const struct OVRCNCtx *const rcn_ctx, uint8_t log2_ctb_s)
{
    const struct OVBuffInfo *const fd = &rcn_ctx->frame_buff;
    const OVSample *src_y  =fd->y  - fd->stride;
    const OVSample *src_cb =fd->cb - fd->stride_c;
    const OVSample *src_cr =fd->cr - fd->stride_c;

    OVSample *dst_y  =  rcn_ctx->ctu_buff.y  - RCN_CTB_STRIDE;
    OVSample *dst_cb =  rcn_ctx->ctu_buff.cb - RCN_CTB_STRIDE;
    OVSample *dst_cr =  rcn_ctx->ctu_buff.cr - RCN_CTB_STRIDE;

    memcpy(dst_y,  src_y , sizeof(OVSample) * OVMIN(((1 << log2_ctb_s) + (1 << log2_ctb_s)), RCN_CTB_STRIDE - 16));
    memcpy(dst_cb, src_cb, sizeof(OVSample) * ((1 << (log2_ctb_s - 1)) + (1 << (log2_ctb_s - 1))));
    memcpy(dst_cr, src_cr, sizeof(OVSample) * ((1 << (log2_ctb_s - 1)) + (1 << (log2_ctb_s - 1))));
}

void
rcn_intra_line_to_ctu(const struct OVRCNCtx *const rcn_ctx, int x_l, uint8_t log2_ctb_s)
{
    const struct OVBuffInfo *const il = &rcn_ctx->intra_line_buff;
    const OVSample *src_y  = il->y  + x_l;
    const OVSample *src_cb = il->cb + (x_l>>1);
    const OVSample *src_cr = il->cr + (x_l>>1);

    OVSample *dst_y  =  rcn_ctx->ctu_buff.y  - RCN_CTB_STRIDE;
    OVSample *dst_cb =  rcn_ctx->ctu_buff.cb - RCN_CTB_STRIDE;
    OVSample *dst_cr =  rcn_ctx->ctu_buff.cr - RCN_CTB_STRIDE;

    memcpy(dst_y,  src_y , sizeof(OVSample) * OVMIN(((1 << log2_ctb_s) + (1 << log2_ctb_s)), RCN_CTB_STRIDE - 16));
    memcpy(dst_cb, src_cb, sizeof(OVSample) * ((1 << (log2_ctb_s - 1)) + (1 << (log2_ctb_s - 1))));
    memcpy(dst_cr, src_cr, sizeof(OVSample) * ((1 << (log2_ctb_s - 1)) + (1 << (log2_ctb_s - 1))));
}

void
rcn_ctu_to_intra_line(const struct OVRCNCtx *const rcn_ctx, int x_l, uint8_t log2_ctu_s)
{
    struct OVBuffInfo *intra_line_binfo = &rcn_ctx->intra_line_buff;

    int max_cu_width_l = 1 << log2_ctu_s;
    int max_cu_width_c = 1 << (log2_ctu_s - 1);

    OVSample *_src;

    _src = &rcn_ctx->ctu_buff.y[(max_cu_width_l - 1) * RCN_CTB_STRIDE];
    memcpy(&intra_line_binfo->y[x_l], _src, sizeof(*_src) << log2_ctu_s);

    _src = &rcn_ctx->ctu_buff.cb[(max_cu_width_c - 1) * RCN_CTB_STRIDE ];
    memcpy(&intra_line_binfo->cb[x_l>>1], _src, sizeof(*_src) << (log2_ctu_s - 1));

    _src = &rcn_ctx->ctu_buff.cr[(max_cu_width_c - 1) * RCN_CTB_STRIDE ];
    memcpy(&intra_line_binfo->cr[x_l>>1], _src, sizeof(*_src) << (log2_ctu_s - 1));
}

void
rcn_write_ctu_to_frame_border(const struct OVRCNCtx *const rcn_ctx,
                              int last_ctu_w, int last_ctu_h)
{
    const struct OVBuffInfo *const fd = &rcn_ctx->frame_buff;
    const OVSample *src_y  = &rcn_ctx->ctu_buff.y [0];
    const OVSample *src_cb = &rcn_ctx->ctu_buff.cb[0];
    const OVSample *src_cr = &rcn_ctx->ctu_buff.cr[0];

    OVSample *dst_y  = fd->y;
    OVSample *dst_cb = fd->cb;
    OVSample *dst_cr = fd->cr;

    for (int i = 0; i < last_ctu_h; ++i) {
        memcpy(dst_y, src_y, sizeof(*src_y) * last_ctu_w);
        dst_y += fd->stride;
        src_y += RCN_CTB_STRIDE;
    }

    for (int i = 0; i < (last_ctu_h >> 1); i++) {
        memcpy(dst_cb, src_cb,  sizeof(*src_cb) * (last_ctu_w >> 1));
        memcpy(dst_cr, src_cr,  sizeof(*src_cr) * (last_ctu_w >> 1));
        dst_cb += fd->stride_c;
        dst_cr += fd->stride_c;
        src_cb += RCN_CTB_STRIDE;
        src_cr += RCN_CTB_STRIDE;
    }
}

static void
ctudec_free_intra_line_buff(struct OVRCNCtx *const rcn_ctx)
{
    struct OVBuffInfo *intra_line_b = &rcn_ctx->intra_line_buff;

    if (intra_line_b->y) {
        ov_freep(&intra_line_b->y);
        ov_freep(&intra_line_b->cb);
        ov_freep(&intra_line_b->cr);
    }
}

static void
ctudec_free_filter_buffers(struct OVRCNCtx *const rcn_ctx)
{
    OVSample** saved_rows_sao    = rcn_ctx->filter_buffers.saved_rows_sao;
    OVSample** saved_rows_alf    = rcn_ctx->filter_buffers.saved_rows_alf;
    OVSample** saved_cols        = rcn_ctx->filter_buffers.saved_cols;
    OVSample** filter_region     = rcn_ctx->filter_buffers.filter_region;

    for(int comp = 0; comp < 3; comp++)
    {
        if(filter_region[comp])     ov_freep(&filter_region[comp]);
        if(saved_rows_sao[comp])    ov_freep(&saved_rows_sao[comp]);
        if(saved_rows_alf[comp])    ov_freep(&saved_rows_alf[comp]);
        if(saved_cols[comp])        ov_freep(&saved_cols[comp]);
    }
}

void
rcn_buff_uninit(struct OVRCNCtx *const rcn_ctx)
{
    ctudec_free_filter_buffers(rcn_ctx);
    ctudec_free_intra_line_buff(rcn_ctx);
}

void
ctudec_alloc_intra_line_buff(struct OVRCNCtx *const rcn_ctx, int nb_ctu_w, uint8_t log2_ctb_s)
{
    struct OVBuffInfo *intra_line_b = &rcn_ctx->intra_line_buff;

    intra_line_b->stride   = nb_ctu_w << log2_ctb_s;
    intra_line_b->stride_c = (nb_ctu_w << log2_ctb_s) >> 1;

    //Free and re-alloc when new ctu width for rectangular entry.
    ctudec_free_intra_line_buff(rcn_ctx);

    intra_line_b->y  = ov_malloc(intra_line_b->stride   * sizeof(*intra_line_b->y));
    intra_line_b->cb = ov_malloc(intra_line_b->stride_c * sizeof(*intra_line_b->cb));
    intra_line_b->cr = ov_malloc(intra_line_b->stride_c * sizeof(*intra_line_b->cr));
}

void
ctudec_save_last_cols(struct OVRCNCtx *const rcn_ctx, int x_pic_l, int y_pic_l, uint8_t is_border_rect)
{
    if (is_border_rect & OV_BOUNDARY_RIGHT_RECT)
        return;

    struct OVFilterBuffers* fb = &rcn_ctx->filter_buffers;
    const int width_l = ( x_pic_l + fb->filter_region_w[0] > rcn_ctx->frame_start->width[0] ) ? ( rcn_ctx->frame_start->width[0] - x_pic_l ) : fb->filter_region_w[0];
    const int height_l = ( y_pic_l + fb->filter_region_h[0] > rcn_ctx->frame_start->height[0] ) ? ( rcn_ctx->frame_start->height[0] - y_pic_l ) : fb->filter_region_h[0];
    const int margin = fb->margin;

    for(int comp = 0; comp < 3; comp++) {
        OVSample* saved_cols = fb->saved_cols[comp];
        OVSample* filter_region = fb->filter_region[comp];
        int stride_filter = fb->filter_region_stride[comp];

        int ratio_luma_chroma = 2;
        int ratio = comp==0 ? 1 : ratio_luma_chroma;
        const int width = width_l/ratio;
        const int height = height_l/ratio;

        for(int ii=0; ii < height; ii++) {
            for(int jj=0; jj < margin; jj++) {
                saved_cols[ii*margin + jj] = filter_region[(ii+margin)*stride_filter + width + jj];
            }
        }
    }
}

void
ctudec_save_last_rows(struct OVRCNCtx *const rcn_ctx, OVSample** saved_rows, int x_l, int x_pic_l, int y_pic_l, uint8_t is_border_rect)
{
    struct OVFilterBuffers* fb = &rcn_ctx->filter_buffers;
    const int width_l = ( x_pic_l + fb->filter_region_w[0] > rcn_ctx->frame_start->width[0] ) ? ( rcn_ctx->frame_start->width[0] - x_pic_l ) : fb->filter_region_w[0];
    const int height_l = ( y_pic_l + fb->filter_region_h[0] > rcn_ctx->frame_start->height[0] ) ? ( rcn_ctx->frame_start->height[0] - y_pic_l ) : fb->filter_region_h[0];
    const int margin = fb->margin;

    for(int comp = 0; comp < 3; comp++) {
        OVSample* saved_rows_comp = saved_rows[comp];
        OVSample* filter_region = fb->filter_region[comp];
        int stride_filter = fb->filter_region_stride[comp];

        int ratio_luma_chroma = 2;
        int ratio = comp==0 ? 1 : ratio_luma_chroma;
        const int width = width_l/ratio;
        const int height = height_l/ratio;
        const int x = x_l/ratio;

        int stride_rows = fb->saved_rows_stride[comp];
        //save pixels in top left corner of ctu filter
        for(int ii=0; ii < margin; ii++) {
            for(int jj=0; jj < margin; jj++) {
                // if ( is_border_rect & VVC_BOUNDARY_RIGHT_TILE)
                if ( 0 )
                    filter_region[ii*stride_filter + jj] = saved_rows_comp[ii*stride_rows];
                else
                    filter_region[ii*stride_filter + jj] = saved_rows_comp[ii*stride_rows + x + width - margin + jj];
            }
        }

        if ( is_border_rect & OV_BOUNDARY_BOTTOM_RECT)
            continue;

        for(int ii=0 ; ii < margin; ii++) {
            memcpy(&saved_rows_comp[ii*stride_rows + x], &filter_region[(height+ii)*stride_filter + margin], width * sizeof(OVSample));
        }
    }
}

void
ctudec_extend_filter_region(struct OVRCNCtx *const rcn_ctx, OVSample** saved_rows, int x_l,
                            int x_pic_l, int y_pic_l, uint8_t bnd_msk)
{

    struct OVFilterBuffers* fb = &rcn_ctx->filter_buffers;

    const int width_l = (x_pic_l + fb->filter_region_w[0] > rcn_ctx->frame_start->width[0]) ? (rcn_ctx->frame_start->width[0] - x_pic_l)
                                                                           : fb->filter_region_w[0];
    const int height_l = (y_pic_l + fb->filter_region_h[0] > rcn_ctx->frame_start->height[0]) ? (rcn_ctx->frame_start->height[0] - y_pic_l)
                                                                            : fb->filter_region_h[0];
    const int margin = fb->margin;

    for (int comp = 0; comp < 3; comp++) {
        int ratio_luma_chroma = 2;
        int ratio = comp==0 ? 1 : ratio_luma_chroma;

        const int width = width_l/ratio;
        const int height = height_l/ratio;
        const int x = x_l/ratio;
        const int x_pic = x_pic_l/ratio;
        const int y_pic = y_pic_l/ratio;

        OVSample* saved_rows_comp = saved_rows[comp];
        OVSample* saved_cols = fb->saved_cols[comp];
        OVSample* filter_region = fb->filter_region[comp];
        int stride_filter = fb->filter_region_stride[comp];

        int stride_pic = rcn_ctx->frame_start->linesize[comp]/sizeof(OVSample);
        OVSample* frame = (OVSample*) rcn_ctx->frame_start->data[comp] + y_pic*stride_pic + x_pic;
        int stride_rows = fb->saved_rows_stride[comp];
        int filter_region_offset = fb->filter_region_offset[comp];

        //*******************************************************/
        //Copy of entire CTU from frame, before border extension
        if (1) {
            OVSample *dst = &filter_region[filter_region_offset];
            const OVSample *src = &frame[0];
            int cpy_s = sizeof(OVSample) * width;
            for (int ii=0; ii < height; ii++) {
                memcpy(dst, src, cpy_s);
                dst += stride_filter;
                src += stride_pic;
            }
        }

        // //*******************************************************/
        //Left margins
        if (bnd_msk & OV_BOUNDARY_LEFT_RECT) {
            OVSample *dst = &filter_region[margin * stride_filter];
            for (int ii=0; ii < height; ii++) {
                const OVSample pad = dst[margin];
                for (int jj=0; jj < margin; jj++) {
                    dst[jj] = pad;
                }
                dst += stride_filter;
            }
        } else {
            int cpy_s = sizeof(OVSample) * margin;
            OVSample *dst = &filter_region[margin * stride_filter];
            const OVSample *src = saved_cols;
            for (int ii=0; ii < height; ii++) {
                for (int jj=0; jj < margin; jj++) {
                    memcpy(dst, src, cpy_s);
                }
                src += margin;
                dst += stride_filter;
            }
        }

        //Right margins
        int h = margin + height;
        int w = margin + width;
        if (bnd_msk & OV_BOUNDARY_RIGHT_RECT) {
            OVSample *dst = &filter_region[margin * stride_filter + w];
            for (int ii=0; ii < height; ii++) {
                const OVSample pad = dst[-1];
                for (int jj=0; jj < margin; jj++) {
                    dst[jj] = pad;
                }
                dst += stride_filter;
            }
        } else {
            int cpy_s = sizeof(OVSample) * margin;
            const OVSample *src = &frame[width];
                  OVSample *dst = &filter_region[margin * stride_filter + w];
            for (int ii=0; ii < height; ii++) {
                for (int jj=0; jj < margin; jj++) {
                    memcpy(dst, src, cpy_s);
                }
                dst += stride_filter;
                src += stride_pic;
            }
        }

        //*******************************************************/
        //Upper margins
        int x_offset_end = 0;

        if (!(bnd_msk & OV_BOUNDARY_RIGHT_RECT))
            x_offset_end = margin;

        if (bnd_msk & OV_BOUNDARY_UPPER_RECT) {
            int cpy_s = sizeof(OVSample) * (width + x_offset_end);
            OVSample *dst = &filter_region[margin];
            const OVSample *src = dst + (margin * stride_filter);
            for (int ii=0; ii < margin; ii++) {
                memcpy(dst, src, cpy_s);
                dst += stride_filter;
            }
        } else {
            int cpy_s = sizeof(OVSample) * (width + x_offset_end);
            OVSample *dst = &filter_region[margin];
            const OVSample *src = &saved_rows_comp[x];
            for (int ii=0; ii < margin; ii++) {
                memcpy(dst, src, cpy_s);
                dst += stride_filter;
                src += stride_rows;
            }
        }

        //Bottom margins
        if (bnd_msk & OV_BOUNDARY_BOTTOM_RECT) {
            int cpy_s = sizeof(OVSample) * (width + 2 * margin);
                  OVSample *dst = &filter_region[h * stride_filter];
            const OVSample *src = dst - stride_filter;
            for (int ii=0; ii < margin; ii++) {
                memcpy(dst, src, cpy_s);
                dst += stride_filter;
            }
        } else {
            int cpy_s = sizeof(OVSample) * (width + 2 * margin);
                  OVSample *dst = &filter_region[h * stride_filter];
            const OVSample *src = &frame[height * stride_pic - margin];
            for (int ii=0; ii < margin; ii++) {
                memcpy(dst, src, cpy_s);
                dst += stride_filter;
                src += stride_pic;
            }
        }

        //*******************************************************/
        //Fill all corners on boudaries
        if (bnd_msk & OV_BOUNDARY_UPPER_RECT) {
            int cpy_s = sizeof(OVSample) * margin;
                  OVSample *dst = &filter_region[0];
            const OVSample *src = &filter_region[margin * stride_filter];
            for (int ii = 0; ii < margin; ii++) {
                memcpy(dst    , src    , cpy_s);
                memcpy(dst + w, src + w, cpy_s);
                dst += stride_filter;
            }
        }

        if (bnd_msk & OV_BOUNDARY_BOTTOM_RECT) {
            int cpy_s = sizeof(OVSample) * margin;
                  OVSample *dst = &filter_region[h * stride_filter];
            const OVSample *src = dst - stride_filter;
            for (int ii = 0; ii < margin; ii++) {
                memcpy(dst    , src    , cpy_s);
                memcpy(dst + w, src + w, cpy_s);
                dst += stride_filter;
            }
        }

        if (bnd_msk & OV_BOUNDARY_LEFT_RECT) {
            OVSample *src  = &filter_region[0];
            OVSample *src2 = &filter_region[h * stride_filter];
            for (int ii = 0; ii < margin; ii++) {
                const OVSample pad  = src [margin];
                const OVSample pad2 = src2[margin];
                for (int jj = 0; jj < margin; jj++) {
                    src [jj] = pad;
                    src2[jj] = pad2;
                }
                src  += stride_filter;
                src2 += stride_filter;
            }
        }

        if (bnd_msk & OV_BOUNDARY_RIGHT_RECT) {
            OVSample *src  = &filter_region[w];
            OVSample *src2 = &filter_region[w + (h * stride_filter)];
            for (int ii=0; ii < margin; ii++) {
                const OVSample pad  = src [-1];
                const OVSample pad2 = src2[-1];
                for (int jj=0; jj < margin; jj++) {
                    src [jj] = pad;
                    src2[jj] = pad2;
                }
                src  += stride_filter;
                src2 += stride_filter;
            }
        }
    }
}

void
ctudec_alloc_filter_buffers(struct OVRCNCtx *const rcn_ctx, int nb_ctu_w, int margin, uint8_t log2_ctb_s)
{
    int ctu_s = 1 << log2_ctb_s;

    struct OVFilterBuffers* fb = &rcn_ctx->filter_buffers;
    OVSample** saved_rows_sao = fb->saved_rows_sao;
    OVSample** saved_rows_alf = fb->saved_rows_alf;
    OVSample** saved_cols    = fb->saved_cols;
    OVSample** filter_region = fb->filter_region;
    fb->margin = margin;

    for(int comp = 0; comp < 3; comp++) {
        int ratio_luma_chroma = 2;
        int ratio = comp==0 ? 1 : ratio_luma_chroma;

        fb->filter_region_w[comp]        = ctu_s / ratio ;
        fb->filter_region_h[comp]        = ctu_s / ratio ;
        fb->filter_region_stride[comp]   = ctu_s / ratio + 2 * margin ;
        fb->filter_region_offset[comp]   = margin * fb->filter_region_stride[comp] + margin;

        if (!filter_region[comp]) {
            int ext_size = fb->filter_region_stride[comp] * (fb->filter_region_h[comp] + 2 * margin + 1);
            filter_region[comp] = ov_malloc(ext_size * sizeof(OVSample));
            saved_cols[comp]    = ov_malloc(fb->filter_region_h[comp] * margin * sizeof(OVSample));
        }

        fb->saved_rows_stride[comp] = nb_ctu_w * ctu_s / ratio; ;
        if (saved_rows_sao[comp]) {
            ov_freep(&saved_rows_sao[comp]);
        }
        saved_rows_sao[comp] = ov_mallocz(margin * fb->saved_rows_stride[comp] * sizeof(OVSample));

        if(saved_rows_alf[comp]){
            ov_freep(&saved_rows_alf[comp]);
        }
        saved_rows_alf[comp] = ov_mallocz(margin * fb->saved_rows_stride[comp] * sizeof(OVSample));
    }

}

