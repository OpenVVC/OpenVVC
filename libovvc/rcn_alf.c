#include <stdint.h>
#include <stdlib.h>

#include "ovutils.h"
#include "ovmem.h"
#include "nvcl_structures.h"

#include "ctudec.h"
#include "rcn_alf.h"
#include "rcn_structures.h"


static inline int clipALF(const int clip, const int32_t ref, const int32_t val0, const int32_t val1)
  {
    // return Clip3<int>(-clip, +clip, val0-ref) + Clip3<int>(-clip, +clip, val1-ref);
    int clip1 = (int) OVMIN( OVMAX(-clip, val0 - ref) , clip);
    int clip2 = (int) OVMIN( OVMAX(-clip, val1 - ref) , clip);
    return clip1 + clip2;
  }


static const int16_t fixed_filter_coeff[ALF_FIXED_FILTER_NUM][MAX_NUM_ALF_LUMA_COEFF] =
{
    { 0,   0,   2,  -3,   1,  -4,   1,   7,  -1,   1,  -1,   5, 0 },
    { 0,   0,   0,   0,   0,  -1,   0,   1,   0,   0,  -1,   2, 0 },
    { 0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   0,   0, 0 },
    { 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -1,   1, 0 },
    { 2,   2,  -7,  -3,   0,  -5,  13,  22,  12,  -3,  -3,  17,  0 },
    { -1,   0,   6,  -8,   1,  -5,   1,  23,   0,   2,  -5,  10,  0 },
    { 0,   0,  -1,  -1,   0,  -1,   2,   1,   0,   0,  -1,   4, 0 },
    { 0,   0,   3, -11,   1,   0,  -1,  35,   5,   2,  -9,   9,  0 },
    { 0,   0,   8,  -8,  -2,  -7,   4,   4,   2,   1,  -1,  25,  0 },
    { 0,   0,   1,  -1,   0,  -3,   1,   3,  -1,   1,  -1,   3, 0 },
    { 0,   0,   3,  -3,   0,  -6,   5,  -1,   2,   1,  -4,  21,  0 },
    { -7,   1,   5,   4,  -3,   5,  11,  13,  12,  -8,  11,  12,  0 },
    { -5,  -3,   6,  -2,  -3,   8,  14,  15,   2,  -7,  11,  16,  0 },
    { 2,  -1,  -6,  -5,  -2,  -2,  20,  14,  -4,   0,  -3,  25,  0 },
    { 3,   1,  -8,  -4,   0,  -8,  22,   5,  -3,   2, -10,  29,  0 },
    { 2,   1,  -7,  -1,   2, -11,  23,  -5,   0,   2, -10,  29,  0 },
    { -6,  -3,   8,   9,  -4,   8,   9,   7,  14,  -2,   8,   9,  0 },
    { 2,   1,  -4,  -7,   0,  -8,  17,  22,   1,  -1,  -4,  23,  0 },
    { 3,   0,  -5,  -7,   0,  -7,  15,  18,  -5,   0,  -5,  27,  0 },
    { 2,   0,   0,  -7,   1, -10,  13,  13,  -4,   2,  -7,  24,  0 },
    { 3,   3, -13,   4,  -2,  -5,   9,  21,  25,  -2,  -3,  12,  0 },
    { -5,  -2,   7,  -3,  -7,   9,   8,   9,  16,  -2,  15,  12,  0 },
    { 0,  -1,   0,  -7,  -5,   4,  11,  11,   8,  -6,  12,  21,  0 },
    { 3,  -2,  -3,  -8,  -4,  -1,  16,  15,  -2,  -3,   3,  26,  0 },
    { 2,   1,  -5,  -4,  -1,  -8,  16,   4,  -2,   1,  -7,  33,  0 },
    { 2,   1,  -4,  -2,   1, -10,  17,  -2,   0,   2, -11,  33,  0 },
    { 1,  -2,   7, -15, -16,  10,   8,   8,  20,  11,  14,  11,  0 },
    { 2,   2,   3, -13, -13,   4,   8,  12,   2,  -3,  16,  24,  0 },
    { 1,   4,   0,  -7,  -8,  -4,   9,   9,  -2,  -2,   8,  29,  0 },
    { 1,   1,   2,  -4,  -1,  -6,   6,   3,  -1,  -1,  -3,  30,  0 },
    { -7,   3,   2,  10,  -2,   3,   7,  11,  19,  -7,   8,  10, 0 },
    { 0,  -2,  -5,  -3,  -2,   4,  20,  15,  -1,  -3,  -1,  22,  0 },
    { 3,  -1,  -8,  -4,  -1,  -4,  22,   8,  -4,   2,  -8,  28,  0 },
    { 0,   3, -14,   3,   0,   1,  19,  17,   8,  -3,  -7,  20,  0 },
    { 0,   2,  -1,  -8,   3,  -6,   5,  21,   1,   1,  -9,  13,  0 },
    { -4,  -2,   8,  20,  -2,   2,   3,   5,  21,   4,   6,   1, 0 },
    { 2,  -2,  -3,  -9,  -4,   2,  14,  16,   3,  -6,   8,  24,  0 },
    { 2,   1,   5, -16,  -7,   2,   3,  11,  15,  -3,  11,  22,  0 },
    { 1,   2,   3, -11,  -2,  -5,   4,   8,   9,  -3,  -2,  26,  0 },
    { 0,  -1,  10,  -9,  -1,  -8,   2,   3,   4,   0,   0,  29,  0 },
    { 1,   2,   0,  -5,   1,  -9,   9,   3,   0,   1,  -7,  20,  0 },
    { -2,   8,  -6,  -4,   3,  -9,  -8,  45,  14,   2, -13,   7, 0 },
    { 1,  -1,  16, -19,  -8,  -4,  -3,   2,  19,   0,   4,  30,  0 },
    { 1,   1,  -3,   0,   2, -11,  15,  -5,   1,   2,  -9,  24,  0 },
    { 0,   1,  -2,   0,   1,  -4,   4,   0,   0,   1,  -4,   7,  0 },
    { 0,   1,   2,  -5,   1,  -6,   4,  10,  -2,   1,  -4,  10,  0 },
    { 3,   0,  -3,  -6,  -2,  -6,  14,   8,  -1,  -1,  -3,  31,  0 },
    { 0,   1,   0,  -2,   1,  -6,   5,   1,   0,   1,  -5,  13,  0 },
    { 3,   1,   9, -19, -21,   9,   7,   6,  13,   5,  15,  21,  0 },
    { 2,   4,   3, -12, -13,   1,   7,   8,   3,   0,  12,  26,  0 },
    { 3,   1,  -8,  -2,   0,  -6,  18,   2,  -2,   3, -10,  23,  0 },
    { 1,   1,  -4,  -1,   1,  -5,   8,   1,  -1,   2,  -5,  10,  0 },
    { 0,   1,  -1,   0,   0,  -2,   2,   0,   0,   1,  -2,   3,  0 },
    { 1,   1,  -2,  -7,   1,  -7,  14,  18,   0,   0,  -7,  21,  0 },
    { 0,   1,   0,  -2,   0,  -7,   8,   1,  -2,   0,  -3,  24,  0 },
    { 0,   1,   1,  -2,   2, -10,  10,   0,  -2,   1,  -7,  23,  0 },
    { 0,   2,   2, -11,   2,  -4,  -3,  39,   7,   1, -10,   9,  0 },
    { 1,   0,  13, -16,  -5,  -6,  -1,   8,   6,   0,   6,  29,  0 },
    { 1,   3,   1,  -6,  -4,  -7,   9,   6,  -3,  -2,   3,  33,  0 },
    { 4,   0, -17,  -1,  -1,   5,  26,   8,  -2,   3, -15,  30,  0 },
    { 0,   1,  -2,   0,   2,  -8,  12,  -6,   1,   1,  -6,  16,  0 },
    { 0,   0,   0,  -1,   1,  -4,   4,   0,   0,   0,  -3,  11,  0 },
    { 0,   1,   2,  -8,   2,  -6,   5,  15,   0,   2,  -7,   9,  0 },
    { 1,  -1,  12, -15,  -7,  -2,   3,   6,   6,  -1,   7,  30,  0 },
};

static const int16_t class_to_filter_mapping[NUM_FIXED_FILTER_SETS][MAX_NUM_ALF_CLASSES] =
{
    { 8,   2,   2,   2,   3,   4,  53,   9,   9,  52,   4,   4,   5,   9,   2,   8,  10,   9,   1,   3,  39,  39,  10,   9,  52 },
    { 11,  12,  13,  14,  15,  30,  11,  17,  18,  19,  16,  20,  20,   4,  53,  21,  22,  23,  14,  25,  26,  26,  27,  28,  10 },
    { 16,  12,  31,  32,  14,  16,  30,  33,  53,  34,  35,  16,  20,   4,   7,  16,  21,  36,  18,  19,  21,  26,  37,  38,  39 },
    { 35,  11,  13,  14,  43,  35,  16,   4,  34,  62,  35,  35,  30,  56,   7,  35,  21,  38,  24,  40,  16,  21,  48,  57,  39 },
    { 11,  31,  32,  43,  44,  16,   4,  17,  34,  45,  30,  20,  20,   7,   5,  21,  22,  46,  40,  47,  26,  48,  63,  58,  10 },
    { 12,  13,  50,  51,  52,  11,  17,  53,  45,   9,  30,   4,  53,  19,   0,  22,  23,  25,  43,  44,  37,  27,  28,  10,  55 },
    { 30,  33,  62,  51,  44,  20,  41,  56,  34,  45,  20,  41,  41,  56,   5,  30,  56,  38,  40,  47,  11,  37,  42,  57,   8 },
    { 35,  11,  23,  32,  14,  35,  20,   4,  17,  18,  21,  20,  20,  20,   4,  16,  21,  36,  46,  25,  41,  26,  48,  49,  58 },
    { 12,  31,  59,  59,   3,  33,  33,  59,  59,  52,   4,  33,  17,  59,  55,  22,  36,  59,  59,  60,  22,  36,  59,  25,  55 },
    { 31,  25,  15,  60,  60,  22,  17,  19,  55,  55,  20,  20,  53,  19,  55,  22,  46,  25,  43,  60,  37,  28,  10,  55,  52 },
    { 12,  31,  32,  50,  51,  11,  33,  53,  19,  45,  16,   4,   4,  53,   5,  22,  36,  18,  25,  43,  26,  27,  27,  28,  10 },
    { 5,   2,  44,  52,   3,   4,  53,  45,   9,   3,   4,  56,   5,   0,   2,   5,  10,  47,  52,   3,  63,  39,  10,   9,  52 },
    { 12,  34,  44,  44,   3,  56,  56,  62,  45,   9,  56,  56,   7,   5,   0,  22,  38,  40,  47,  52,  48,  57,  39,  10,   9 },
    { 35,  11,  23,  14,  51,  35,  20,  41,  56,  62,  16,  20,  41,  56,   7,  16,  21,  38,  24,  40,  26,  26,  42,  57,  39 },
    { 33,  34,  51,  51,  52,  41,  41,  34,  62,   0,  41,  41,  56,   7,   5,  56,  38,  38,  40,  44,  37,  42,  57,  39,  10 },
    { 16,  31,  32,  15,  60,  30,   4,  17,  19,  25,  22,  20,   4,  53,  19,  21,  22,  46,  25,  55,  26,  48,  63,  58,  55 },
};


void
rcn_alf_create(OVCTUDec *const ctudec, RCNALF* alf)
{
    const OVPartInfo *const pinfo = ctudec->part_ctx;
    uint8_t log2_ctb_size = pinfo->log2_ctu_s;
    int ctu_width  = 1 << log2_ctb_size;  

    //BITDEPTH: uniquement pour bitdepth 10
    int bit_depth = 10; 
    int shift_luma = bit_depth - 8;
    int shift_chroma = bit_depth - 8;

    /* FIXME would be better to get (1 << bitdepth) - 1 for clipping */
    alf->alf_clipping_values[CHANNEL_TYPE_LUMA][0] = 1 << bit_depth;
    alf->alf_clipping_values[CHANNEL_TYPE_CHROMA][0] = 1 << bit_depth;
    for(int i = 1; i < MAX_ALF_NUM_CLIP_VAL; ++i) {
        alf->alf_clipping_values[CHANNEL_TYPE_LUMA][i] = 1 << (7 - 2 * i + shift_luma);
        alf->alf_clipping_values[CHANNEL_TYPE_CHROMA][i] = 1 << (7 - 2 * i + shift_chroma);
    }

    /* FIXME Use two fixed 32x32 tables instead of alloc (one for tr_idx one for cl_idx
     * We could also use a joint idx for both.
     */
    if (alf->classifier == 0) {
        alf->classifier = ov_malloc(sizeof(ALFClassifier*)*ctu_width >> 2);
        for (int i = 0; i < ctu_width/4; i++) {
            alf->classifier[i] = ov_malloc(sizeof(ALFClassifier)*ctu_width >> 2 );
        }
    }

    for (int filter_set_idx = 0; filter_set_idx < NUM_FIXED_FILTER_SETS; filter_set_idx++) {
        for (int class_idx = 0; class_idx < MAX_NUM_ALF_CLASSES; class_idx++) {
            int fixed_filter_idx = class_to_filter_mapping[filter_set_idx][class_idx];
            for (int i = 0; i < MAX_NUM_ALF_LUMA_COEFF - 1; i++) {
                alf->fixed_filter_coeff_dec[filter_set_idx][class_idx * MAX_NUM_ALF_LUMA_COEFF + i] = fixed_filter_coeff[fixed_filter_idx][i];
            }
            alf->fixed_filter_coeff_dec[filter_set_idx][class_idx * MAX_NUM_ALF_LUMA_COEFF + MAX_NUM_ALF_LUMA_COEFF - 1] = (1 << (NUM_BITS - 1));
        }
    }

    for (int i = 0; i < MAX_NUM_ALF_LUMA_COEFF * MAX_NUM_ALF_CLASSES; i++) {
        alf->clip_default[i] = alf->alf_clipping_values[CHANNEL_TYPE_LUMA][0];
    }
}


void alf_reconstructCoeff_luma(RCNALF* alf, const struct OVALFData* alf_data)
{
    int factor = 1 << (NUM_BITS - 1);
    int num_classes =  MAX_NUM_ALF_CLASSES;
    int num_coeff = 13;
    int num_coeff_minus1 = num_coeff - 1;
    
    int num_filters = alf_data->alf_luma_num_filters_signalled_minus1 + 1 ;
    int16_t* coeff ;
    int16_t* clip;

    coeff = (int16_t*) alf_data->alf_luma_coeff;
    clip = (int16_t*) alf_data->alf_luma_clip_idx;
    for( int filter_idx = 0; filter_idx < num_filters; filter_idx++ )
    {
        coeff[filter_idx* MAX_NUM_ALF_LUMA_COEFF + num_coeff_minus1] = factor;
    }

    for( int class_idx = 0; class_idx < num_classes; class_idx++ )
    {
        int filter_idx = alf_data->alf_luma_coeff_delta_idx[class_idx];
        for (int coeffIdx = 0; coeffIdx < num_coeff_minus1; ++coeffIdx)
        {
            alf->coeff_final[class_idx * MAX_NUM_ALF_LUMA_COEFF + coeffIdx] = coeff[filter_idx * MAX_NUM_ALF_LUMA_COEFF + coeffIdx];
        }
        alf->coeff_final[class_idx* MAX_NUM_ALF_LUMA_COEFF + num_coeff_minus1] = factor;
        alf->clip_final[class_idx* MAX_NUM_ALF_LUMA_COEFF + num_coeff_minus1] = alf->alf_clipping_values[0][0];
        for( int coeffIdx = 0; coeffIdx < num_coeff_minus1; ++coeffIdx )
        {
            int clipIdx = alf_data->alf_luma_clip_flag ? (clip + filter_idx * MAX_NUM_ALF_LUMA_COEFF)[coeffIdx] : 0;
            (alf->clip_final + class_idx * MAX_NUM_ALF_LUMA_COEFF)[coeffIdx] = alf->alf_clipping_values[0][clipIdx];
        }
        alf->clip_final[class_idx* MAX_NUM_ALF_LUMA_COEFF + num_coeff_minus1] = alf->alf_clipping_values[0][0];
    }
}

void alf_reconstructCoeff_chroma(RCNALF* alf, const struct OVALFData* alf_data)
{
    int factor = 1 << (NUM_BITS - 1);
    int num_classes = 1;
    int num_coeff = 7 ;
    int num_coeff_minus1 = num_coeff - 1;
    const int num_alts = alf_data->alf_chroma_num_alt_filters_minus1 + 1;
    
    int16_t* coeff ;
    int16_t* clip;

    for( int alt_idx = 0; alt_idx < num_alts; ++ alt_idx )
    {
        coeff = (int16_t*) alf_data->alf_chroma_coeff[alt_idx];
        clip = (int16_t*) alf_data->alf_chroma_clip_idx[alt_idx];
        for( int coeffIdx = 0; coeffIdx < num_coeff_minus1; ++coeffIdx )
        {
            alf->chroma_coeff_final[alt_idx][coeffIdx] = coeff[coeffIdx];
            int clipIdx = alf_data->alf_chroma_clip_flag ? clip[coeffIdx] : 0;
            alf->chroma_clip_final[alt_idx][coeffIdx] = alf->alf_clipping_values[1][clipIdx];
        }
        alf->chroma_coeff_final[alt_idx][num_coeff_minus1] = factor;
        alf->chroma_clip_final[alt_idx][num_coeff_minus1] = alf->alf_clipping_values[1][0];

        for( int class_idx = 0; class_idx < num_classes; class_idx++ )
        {
            int filter_idx = alf_data->alf_luma_coeff_delta_idx[class_idx];
            for (int coeffIdx = 0; coeffIdx < num_coeff_minus1; ++coeffIdx)
            {
                alf->coeff_final[class_idx * MAX_NUM_ALF_LUMA_COEFF + coeffIdx] = coeff[filter_idx * MAX_NUM_ALF_LUMA_COEFF + coeffIdx];
            }
            alf->coeff_final[class_idx* MAX_NUM_ALF_LUMA_COEFF + num_coeff_minus1] = factor;
            alf->clip_final[class_idx* MAX_NUM_ALF_LUMA_COEFF + num_coeff_minus1] = alf->alf_clipping_values[1][0];
            for( int coeffIdx = 0; coeffIdx < num_coeff_minus1; ++coeffIdx )
            {
                int clipIdx = alf_data->alf_chroma_clip_flag ? (clip + filter_idx * MAX_NUM_ALF_LUMA_COEFF)[coeffIdx] : 0;
                (alf->clip_final + class_idx * MAX_NUM_ALF_LUMA_COEFF)[coeffIdx] = alf->alf_clipping_values[1][clipIdx];
            }
            alf->clip_final[class_idx* MAX_NUM_ALF_LUMA_COEFF + num_coeff_minus1] = alf->alf_clipping_values[1][0];
        }
    }
}

void rcn_alf_reconstruct_coeff_APS(RCNALF* alf, OVCTUDec *const ctudec, uint8_t luma_flag, uint8_t chroma_flag)
{ 
    if (luma_flag){
        for (int i = 0; i < ctudec->alf_info.num_alf_aps_ids_luma; i++)
        {
            const struct OVALFData* alf_data = ctudec->alf_info.aps_alf_data[i];

            alf_reconstructCoeff_luma(alf, alf_data);

            for(int j = 0; j < MAX_NUM_ALF_CLASSES * MAX_NUM_ALF_LUMA_COEFF; j++){
                alf->coeff_aps_luma[0][j] = alf->coeff_final[j];
                alf->clip_aps_luma[0][j] = alf->clip_final[j];
            }
        }
    }

    if (chroma_flag){
        const struct OVALFData* alf_data_c = ctudec->alf_info.aps_alf_data_c;
        alf_reconstructCoeff_chroma(alf, alf_data_c);
    }
}



void
rcn_alf_derive_classificationBlk(ALFClassifier **classifier,
                                 int16_t *const src, const int stride, const Area blk,
                                 const int shift, const int ctu_height, int virbnd_pos)
{

    static const int th[16] = { 0, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 4 };
    const int maxActivity = 15;
    int laplacian[NUM_DIRECTIONS][CLASSIFICATION_BLK_SIZE + 5][CLASSIFICATION_BLK_SIZE + 5];

    int fl = 2;
    int flP1 = fl + 1;
    int fl2 = 2 * fl;

    int main_dir, secondary_dir, dir_temp_hv, dir_temp_d;

    int pixY;
    int height = blk.height + fl2;
    int width  = blk.width + fl2;

    for( int i = 0; i < height; i += 2 ) {
        int yoffset = ( i + 1 - flP1 ) * stride - flP1;
        const int16_t *src0 = &src[yoffset - stride];
        const int16_t *src1 = &src[yoffset];
        const int16_t *src2 = &src[yoffset + stride];
        const int16_t *src3 = &src[yoffset + stride * 2];

        const int y = blk.y - 2 + i;
        if (y > 0 && (y & (ctu_height - 1)) == virbnd_pos - 2) {
            src3 = &src[yoffset + stride];
        } else if (y > 0 && (y & (ctu_height - 1)) == virbnd_pos) {
            src0 = &src[yoffset];
        }

        int* pYver = laplacian[VER][i];
        int* pYhor = laplacian[HOR][i];
        int* pYdig0 = laplacian[DIAG0][i];
        int* pYdig1 = laplacian[DIAG1][i];

        for( int j = 0; j < width; j += 2 ) {
            // pixY = j + 1 + posX;
            pixY = j + 1 ;
            const int16_t *pY = src1 + pixY;
            const int16_t *pYdown = src0 + pixY;
            const int16_t *pYup = src2 + pixY;
            const int16_t *pYup2 = src3 + pixY;

            const int16_t y0 = pY[0] << 1;
            const int16_t yup1 = pYup[1] << 1;

            //modification des buffers laplacian ici 
            pYver[j]  = abs(y0 - pYdown[0]  - pYup[0])   + abs(yup1 - pY[1]    - pYup2[1]);
            pYhor[j]  = abs(y0 - pY[1]      - pY[-1])    + abs(yup1 - pYup[2]  - pYup[0]);
            pYdig0[j] = abs(y0 - pYdown[-1] - pYup[1])   + abs(yup1 - pY[0]    - pYup2[2]);
            pYdig1[j] = abs(y0 - pYup[-1]   - pYdown[1]) + abs(yup1 - pYup2[0] - pY[2]);
        }

        for( int j = 6; j < width; j += 4 ) {
            int jM6 = j - 6;
            int jM4 = j - 4;
            int jM2 = j - 2;

            pYver[jM6]  += pYver[jM4]  + pYver[jM2]  + pYver[j];
            pYhor[jM6]  += pYhor[jM4]  + pYhor[jM2]  + pYhor[j];
            pYdig0[jM6] += pYdig0[jM4] + pYdig0[jM2] + pYdig0[j];
            pYdig1[jM6] += pYdig1[jM4] + pYdig1[jM2] + pYdig1[j];
        }
    }

    // classification block size
    const int clsSizeY = 4;
    const int clsSizeX = 4;

    for( int i = 0; i < blk.height; i += clsSizeY ) {
        int* pYver = laplacian[VER][i];
        int* pYver2 = laplacian[VER][i + 2];
        int* pYver4 = laplacian[VER][i + 4];
        int* pYver6 = laplacian[VER][i + 6];

        int* pYhor = laplacian[HOR][i];
        int* pYhor2 = laplacian[HOR][i + 2];
        int* pYhor4 = laplacian[HOR][i + 4];
        int* pYhor6 = laplacian[HOR][i + 6];

        int* pYdig0 = laplacian[DIAG0][i];
        int* pYdig02 = laplacian[DIAG0][i + 2];
        int* pYdig04 = laplacian[DIAG0][i + 4];
        int* pYdig06 = laplacian[DIAG0][i + 6];

        int* pYdig1 = laplacian[DIAG1][i];
        int* pYdig12 = laplacian[DIAG1][i + 2];
        int* pYdig14 = laplacian[DIAG1][i + 4];
        int* pYdig16 = laplacian[DIAG1][i + 6];

        for( int j = 0; j < blk.width; j += clsSizeX ) {
            int sumV = 0; int sumH = 0; int sumD0 = 0; int sumD1 = 0;

            if (((i + blk.y) % ctu_height) == (virbnd_pos - 4)) {
                sumV = pYver[j] + pYver2[j] + pYver4[j];
                sumH = pYhor[j] + pYhor2[j] + pYhor4[j];
                sumD0 = pYdig0[j] + pYdig02[j] + pYdig04[j];
                sumD1 = pYdig1[j] + pYdig12[j] + pYdig14[j];
            } else if (((i + blk.y) % ctu_height) == virbnd_pos) {
                sumV = pYver2[j] + pYver4[j] + pYver6[j];
                sumH = pYhor2[j] + pYhor4[j] + pYhor6[j];
                sumD0 = pYdig02[j] + pYdig04[j] + pYdig06[j];
                sumD1 = pYdig12[j] + pYdig14[j] + pYdig16[j];
            } else {
                sumV = pYver[j] + pYver2[j] + pYver4[j] + pYver6[j];
                sumH = pYhor[j] + pYhor2[j] + pYhor4[j] + pYhor6[j];
                sumD0 = pYdig0[j] + pYdig02[j] + pYdig04[j] + pYdig06[j];
                sumD1 = pYdig1[j] + pYdig12[j] + pYdig14[j] + pYdig16[j];
            }

            int tempAct = sumV + sumH;
            int activity = 0;

            const int y = (i + blk.y) & (ctu_height - 1);

            if (y == virbnd_pos - 4 || y == virbnd_pos) {
                activity = (int16_t) OVMIN( OVMAX(0, (tempAct * 96) >> shift) , maxActivity);
            } else {
                activity = (int16_t) OVMIN( OVMAX(0, (tempAct * 64) >> shift) , maxActivity);
            }

            int class_idx = th[activity];

            int hv1, hv0, d1, d0, hvd1, hvd0;

            if( sumV > sumH ) {
                hv1 = sumV;
                hv0 = sumH;
                dir_temp_hv = 1;
            } else {
                hv1 = sumH;
                hv0 = sumV;
                dir_temp_hv = 3;
            }

            if( sumD0 > sumD1 ) {
                d1 = sumD0;
                d0 = sumD1;
                dir_temp_d = 0;
            } else {
                d1 = sumD1;
                d0 = sumD0;
                dir_temp_d = 2;
            }

            if( (uint32_t)d1 * (uint32_t)hv0 > (uint32_t)hv1 * (uint32_t)d0 ) {
                hvd1 = d1;
                hvd0 = d0;
                main_dir = dir_temp_d;
                secondary_dir = dir_temp_hv;
            } else {
                hvd1 = hv1;
                hvd0 = hv0;
                main_dir = dir_temp_hv;
                secondary_dir = dir_temp_d;
            }

            int directionStrength = 0;

            if( hvd1 > 2 * hvd0 ) {
                directionStrength = 1;

            }

            if( hvd1 * 2 > 9 * hvd0 ) {
                directionStrength = 2;
            }

            if( directionStrength ) {
                class_idx += ( ( ( main_dir & 0x1 ) << 1 ) + directionStrength ) * 5;
            }

            static const int transposeTable[8] = { 0, 1, 0, 2, 2, 3, 1, 3 };
            int transpose_idx = transposeTable[main_dir * 2 + ( secondary_dir >> 1 )];

            int yOffset = (i + blk.y) % ctu_height;
            int xOffset = (j + blk.x) % ctu_height;
            ALFClassifier alf_class;
            alf_class.class_idx = class_idx;
            alf_class.transpose_idx = transpose_idx;
            classifier[yOffset>>2][xOffset>>2] = alf_class;
        }
    }
}


void rcn_alf_derive_classification(RCNALF *alf, int16_t *const rcn_img, const int stride, Area blk, int ctu_width, int pic_h )
{
    ALFClassifier** classifier = alf->classifier;
    int height = blk.y + blk.height;
    int width = blk.x + blk.width;
    //BITDEPTH: uniquement pour bitdepth 10
    int bit_depth = 10;

    for( int i = blk.y; i < height; i += CLASSIFICATION_BLK_SIZE )
    {
        int nHeight = OVMIN( i + CLASSIFICATION_BLK_SIZE, height ) - i;

        for( int j = blk.x; j < width; j += CLASSIFICATION_BLK_SIZE )
        {
            int nWidth = OVMIN( j + CLASSIFICATION_BLK_SIZE, width ) - j;
            Area blk_class;
            blk_class.x = j;
            blk_class.y = i;
            blk_class.width = nWidth;
            blk_class.height = nHeight;

            int16_t* rcn_img_class = rcn_img + (i - blk.y) * stride + (j - blk.x);

            rcn_alf_derive_classificationBlk(classifier, rcn_img_class, stride, blk_class,
                                             bit_depth + 4, ctu_width,
                                             (blk.height<ctu_width) ? pic_h : blk.height - ALF_VB_POS_ABOVE_CTUROW_LUMA);
        }
    }
}

void cc_alf_filterBlk(int16_t * chroma_dst, int16_t * luma_src, const int chr_stride, const int luma_stride, 
                        const Area blk_dst, const uint8_t c_id, const int16_t *filt_coeff, 
                        const int vbCTUHeight, int vbPos)
{
  const int clsSizeY           = 4;
  const int clsSizeX           = 4;
  //ATTENTION: scaleX et Y fixed to 1 (en 4 2 0)
  const int scaleX             = 1;
  const int scaleY             = 1;

  for( int i = 0; i < blk_dst.height; i += clsSizeY )
  {
    for( int j = 0; j < blk_dst.width; j += clsSizeX )
    {
      for( int ii = 0; ii < clsSizeY; ii++ )
      {
        int row       = ii;
        int col       = j;
        int16_t *srcSelf  = chroma_dst + col + row * chr_stride;

        int offset1 = luma_stride;
        int offset2 = -luma_stride;
        int offset3 = 2 * luma_stride;
        row <<= scaleY;
        col <<= scaleX;
        const int16_t *srcCross = luma_src + col + row * luma_stride;

        int pos = ((blk_dst.y + i + ii) << scaleY) & (vbCTUHeight - 1);
        if (scaleY == 0 && (pos == vbPos || pos == vbPos + 1))
        {
          continue;
        }
        if (pos == (vbPos - 2) || pos == (vbPos + 1))
        {
          offset3 = offset1;
        }
        else if (pos == (vbPos - 1) || pos == vbPos)
        {
          offset1 = 0;
          offset2 = 0;
          offset3 = 0;
        }

        for (int jj = 0; jj < clsSizeX; jj++)
        {
          const int jj2     = (jj << scaleX);
          const int offset0 = 0;

          int sum = 0;
          const int16_t currSrcCross = srcCross[offset0 + jj2];
          sum += filt_coeff[0] * (srcCross[offset2 + jj2    ] - currSrcCross);
          sum += filt_coeff[1] * (srcCross[offset0 + jj2 - 1] - currSrcCross);
          sum += filt_coeff[2] * (srcCross[offset0 + jj2 + 1] - currSrcCross);
          sum += filt_coeff[3] * (srcCross[offset1 + jj2 - 1] - currSrcCross);
          sum += filt_coeff[4] * (srcCross[offset1 + jj2    ] - currSrcCross);
          sum += filt_coeff[5] * (srcCross[offset1 + jj2 + 1] - currSrcCross);
          sum += filt_coeff[6] * (srcCross[offset3 + jj2    ] - currSrcCross);

          const int scale_bits = 7;
          sum = (sum + ((1 << scale_bits ) >> 1)) >> scale_bits;    

          //BITDEPTH: uniquement pour bitdepth 10
          const int bit_depth = 10;
          const int offset = 1 << bit_depth >> 1;
          sum = OVMAX( OVMIN( sum + offset, (1<<bit_depth) - 1 ), 0) - offset;

          sum += srcSelf[jj];
          srcSelf[jj] = OVMAX( OVMIN( sum, (1<<bit_depth) - 1 ), 0) ;
        }
      }
    }
    chroma_dst += chr_stride * clsSizeY;
    luma_src += luma_stride * clsSizeY << scaleY;
  }
}


// dst   : buffer Frame output, pointing to the begining of CTU.
// src   : filter buffer pre-ALF (of size CTU)
// blk_dst: location and dimension of destination block in frame
// blk   : location and dimension of destination block in filter buffer
void alf_filter_c(int16_t *const dst, const int16_t *const src,
                  const int dst_stride, const int src_stride,
                  Area blk_dst,
                  const int16_t *const filter_set, const int16_t *const clip_set,
                  const int ctu_height, int virbnd_pos)
{
    const int shift = NUM_BITS - 1;
    const int offset = 1 << ( shift - 1 );

    const int blk_h = 4;
    const int blk_w = 4;

    int dst_blk_stride = dst_stride * blk_h;
    int src_blk_stride = src_stride * blk_h;

    int16_t* dst0 = dst ;
    int16_t* dst1 = dst + dst_stride;
    int i;

    const int16_t *lm_src0 = src;
    const int16_t *lm_src1 = lm_src0 + src_stride;
    const int16_t *lm_src2 = lm_src0 - src_stride;
    const int16_t *lm_src3 = lm_src1 + src_stride;
    const int16_t *lm_src4 = lm_src2 - src_stride;
    const int16_t *lm_src5 = lm_src3 + src_stride;
    const int16_t *lm_src6 = lm_src4 - src_stride;

    for (i = 0; i < blk_dst.height; i += blk_h) {
        int j;
        for (j = 0; j < blk_dst.width; j += blk_w) {
            int k;
            for (k = 0; k < blk_h; k++) {
                int l;
                const int16_t *src_0 = lm_src0 + j + k * src_stride;
                const int16_t *src_1 = lm_src1 + j + k * src_stride;
                const int16_t *src_2 = lm_src2 + j + k * src_stride;
                const int16_t *src_3 = lm_src3 + j + k * src_stride;
                const int16_t *src_4 = lm_src4 + j + k * src_stride;
                const int16_t *src_5 = lm_src5 + j + k * src_stride;
                const int16_t *src_6 = lm_src6 + j + k * src_stride;

                dst1 = dst0 + j + k * dst_stride;

                int yVb = (blk_dst.y + i + k) & (ctu_height - 1);
                if (yVb < virbnd_pos && (yVb >= virbnd_pos - 2)) {
                    src_1 = (yVb == virbnd_pos - 1) ? src_0 : src_1;
                    src_3 = (yVb >= virbnd_pos - 2) ? src_1 : src_3;
                    src_5 = (yVb >= virbnd_pos - 3) ? src_3 : src_5;

                    src_2 = (yVb == virbnd_pos - 1) ? src_0 : src_2;
                    src_4 = (yVb >= virbnd_pos - 2) ? src_2 : src_4;
                    src_6 = (yVb >= virbnd_pos - 3) ? src_4 : src_6;
                } else if (yVb >= virbnd_pos && (yVb <= virbnd_pos + 1)) {
                    src_2 = (yVb == virbnd_pos    ) ? src_0 : src_2;
                    src_4 = (yVb <= virbnd_pos + 1) ? src_2 : src_4;
                    src_6 = (yVb <= virbnd_pos + 2) ? src_4 : src_6;

                    src_1 = (yVb == virbnd_pos    ) ? src_0 : src_1;
                    src_3 = (yVb <= virbnd_pos + 1) ? src_1 : src_3;
                    src_5 = (yVb <= virbnd_pos + 2) ? src_3 : src_5;
                }

                uint8_t isNearVBabove = yVb < virbnd_pos && (yVb >= virbnd_pos - 1);
                uint8_t isNearVBbelow = yVb >= virbnd_pos && (yVb <= virbnd_pos);

                for (l = 0; l < blk_w; l++) {
                    int sum = 0;
                    const int16_t curr = src_0[0];
                    sum += filter_set[0] * clipALF(clip_set[0], curr, src_3[ 0], src_4[ 0]);
                    sum += filter_set[1] * clipALF(clip_set[1], curr, src_1[ 1], src_2[-1]);
                    sum += filter_set[2] * clipALF(clip_set[2], curr, src_1[ 0], src_2[ 0]);
                    sum += filter_set[3] * clipALF(clip_set[3], curr, src_1[-1], src_2[ 1]);
                    sum += filter_set[4] * clipALF(clip_set[4], curr, src_0[ 2], src_0[-2]);
                    sum += filter_set[5] * clipALF(clip_set[5], curr, src_0[ 1], src_0[-1]);

                    if (!(isNearVBabove || isNearVBbelow)) {
                        sum = (sum + offset) >> shift;
                    } else {
                        //Rounding offset fix
                        sum = (sum + (1 << ((shift + 3) - 1))) >> (shift + 3);
                    }

                    sum += curr;
                    dst1[l] = OVMAX( OVMIN( sum, (1<<10) - 1 ), 0);

                    src_0++;
                    src_1++;
                    src_2++;
                    src_3++;
                    src_4++;
                    src_5++;
                    src_6++;
                }
            }
        }

        dst0 += dst_blk_stride;
        dst1 += dst_blk_stride;

        lm_src0 += src_blk_stride;
        lm_src1 += src_blk_stride;
        lm_src2 += src_blk_stride;
        lm_src3 += src_blk_stride;
        lm_src4 += src_blk_stride;
        lm_src5 += src_blk_stride;
        lm_src6 += src_blk_stride;
    }
}

void alf_filterBlkLuma(ALFClassifier **classifier, int16_t *const dst, int16_t *const src, const int dstStride, const int srcStride,
                         Area blk_dst, const int16_t *filter_set, const int16_t *clip_set,
                         const int ctu_height, int virbnd_pos)
{
  const int16_t *pImg0, *pImg1, *pImg2, *pImg3, *pImg4, *pImg5, *pImg6;

  const int16_t *coef = filter_set;
  const int16_t *clip = clip_set;

  const int shift = NUM_BITS - 1;
  const int offset = 1 << ( shift - 1 );

  int transpose_idx = 0;
  const int clsSizeY = 4;
  const int clsSizeX = 4;
  ALFClassifier *pClass = 0;

  int dstStride2 = dstStride * clsSizeY;
  int srcStride2 = srcStride * clsSizeY;

  int16_t * _src = src;
  int16_t * _dst = dst;

  int16_t* pRec0 = dst ;
  int16_t* pRec1 = pRec0 + dstStride;

  for( int i = 0; i < blk_dst.height; i += clsSizeY )
  {
    pClass = classifier[i>>2];
    for( int j = 0; j < blk_dst.width; j += clsSizeX )
    {
      ALFClassifier cl = pClass[j>>2];
      transpose_idx = cl.transpose_idx;
      coef = filter_set + cl.class_idx * MAX_NUM_ALF_LUMA_COEFF;
      clip = clip_set + cl.class_idx * MAX_NUM_ALF_LUMA_COEFF;

      static const uint8_t shuffle_lut[4][12] = {
          {  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11},
          {  9,  4, 10,  8,  1,  5, 11,  7,  3,  0,  2,  6},
          {  0,  3,  2,  1,  8,  7,  6,  5,  4,  9, 10, 11},
          {  9,  8, 10,  4,  3,  7, 11,  5,  1,  0,  2,  6}
      };

      const int filt_coeff[MAX_NUM_ALF_LUMA_COEFF] = {
          coef[shuffle_lut[transpose_idx][0]],
          coef[shuffle_lut[transpose_idx][1]],
          coef[shuffle_lut[transpose_idx][2]],
          coef[shuffle_lut[transpose_idx][3]],
          coef[shuffle_lut[transpose_idx][4]],
          coef[shuffle_lut[transpose_idx][5]],
          coef[shuffle_lut[transpose_idx][6]],
          coef[shuffle_lut[transpose_idx][7]],
          coef[shuffle_lut[transpose_idx][8]],
          coef[shuffle_lut[transpose_idx][9]],
          coef[shuffle_lut[transpose_idx][10]],
          coef[shuffle_lut[transpose_idx][11]],
          coef[12]
      };

      const int filt_clip[MAX_NUM_ALF_LUMA_COEFF] = {
          clip[shuffle_lut[transpose_idx][0]],
          clip[shuffle_lut[transpose_idx][1]],
          clip[shuffle_lut[transpose_idx][2]],
          clip[shuffle_lut[transpose_idx][3]],
          clip[shuffle_lut[transpose_idx][4]],
          clip[shuffle_lut[transpose_idx][5]],
          clip[shuffle_lut[transpose_idx][6]],
          clip[shuffle_lut[transpose_idx][7]],
          clip[shuffle_lut[transpose_idx][8]],
          clip[shuffle_lut[transpose_idx][9]],
          clip[shuffle_lut[transpose_idx][10]],
          clip[shuffle_lut[transpose_idx][11]],
          clip[12]
      };

      for( int ii = 0; ii < clsSizeY; ii++ ) {
        pImg0 = _src + j + ii * srcStride;
        pImg1 = pImg0 + srcStride;
        pImg2 = pImg0 - srcStride;
        pImg3 = pImg1 + srcStride;
        pImg4 = pImg2 - srcStride;
        pImg5 = pImg3 + srcStride;
        pImg6 = pImg4 - srcStride;

        pRec1 = pRec0 + j + ii * dstStride;

        int yVb = (blk_dst.y + i + ii) & (ctu_height - 1);
        if (yVb < virbnd_pos && (yVb >= virbnd_pos - 4)) {
          pImg1 = (yVb == virbnd_pos - 1) ? pImg0 : pImg1;
          pImg3 = (yVb >= virbnd_pos - 2) ? pImg1 : pImg3;
          pImg5 = (yVb >= virbnd_pos - 3) ? pImg3 : pImg5;

          pImg2 = (yVb == virbnd_pos - 1) ? pImg0 : pImg2;
          pImg4 = (yVb >= virbnd_pos - 2) ? pImg2 : pImg4;
          pImg6 = (yVb >= virbnd_pos - 3) ? pImg4 : pImg6;
        } else if (yVb >= virbnd_pos && (yVb <= virbnd_pos + 3)) {
          pImg2 = (yVb == virbnd_pos) ? pImg0 : pImg2;
          pImg4 = (yVb <= virbnd_pos + 1) ? pImg2 : pImg4;
          pImg6 = (yVb <= virbnd_pos + 2) ? pImg4 : pImg6;

          pImg1 = (yVb == virbnd_pos) ? pImg0 : pImg1;
          pImg3 = (yVb <= virbnd_pos + 1) ? pImg1 : pImg3;
          pImg5 = (yVb <= virbnd_pos + 2) ? pImg3 : pImg5;
        }


        uint8_t isNearVBabove = yVb < virbnd_pos && (yVb >= virbnd_pos - 1);
        uint8_t isNearVBbelow = yVb >= virbnd_pos && (yVb <= virbnd_pos);
        for( int jj = 0; jj < clsSizeX; jj++ ) {
          int sum = 0;
          const int16_t curr = pImg0[+0];
          sum += filt_coeff[0] * ( clipALF(filt_clip[0], curr, pImg5[+0], pImg6[+0]) );
          sum += filt_coeff[1] * ( clipALF(filt_clip[1], curr, pImg3[+1], pImg4[-1]) );

          sum += filt_coeff[2] * ( clipALF(filt_clip[2], curr, pImg3[+0], pImg4[+0]) );
          sum += filt_coeff[3] * ( clipALF(filt_clip[3], curr, pImg3[-1], pImg4[+1]) );

          sum += filt_coeff[4] * ( clipALF(filt_clip[4], curr, pImg1[+2], pImg2[-2]) );
          sum += filt_coeff[5] * ( clipALF(filt_clip[5], curr, pImg1[+1], pImg2[-1]) );

          sum += filt_coeff[6] * ( clipALF(filt_clip[6], curr, pImg1[+0], pImg2[+0]) );
          sum += filt_coeff[7] * ( clipALF(filt_clip[7], curr, pImg1[-1], pImg2[+1]) );

          sum += filt_coeff[8] * ( clipALF(filt_clip[8], curr, pImg1[-2], pImg2[+2]) );
          sum += filt_coeff[9] * ( clipALF(filt_clip[9], curr, pImg0[+3], pImg0[-3]) );

          sum += filt_coeff[10] * ( clipALF(filt_clip[10], curr, pImg0[+2], pImg0[-2]) );
          sum += filt_coeff[11] * ( clipALF(filt_clip[11], curr, pImg0[+1], pImg0[-1]) );

          if (!(isNearVBabove || isNearVBbelow)) {
            sum = (sum + offset) >> shift;
          } else {
            sum = (sum + (1 << ((shift + 3) - 1))) >> (shift + 3);
          }

          sum += curr;
          pRec1[jj] = OVMAX( OVMIN( sum, (1<<10) - 1 ), 0);

          pImg0++;
          pImg1++;
          pImg2++;
          pImg3++;
          pImg4++;
          pImg5++;
          pImg6++;
        }
      }
    }

    pRec0 += dstStride2;
    pRec1 += dstStride2;

    _src += srcStride2;
    _dst += dstStride2;
  }
}



void rcn_alf_filter_line(OVCTUDec *const ctudec, int nb_ctu_w, uint16_t ctb_y_pic)
{
    struct ALFInfo* alf_info = &ctudec->alf_info;
    if (!alf_info->alf_luma_enabled_flag && !alf_info->alf_cb_enabled_flag && !alf_info->alf_cr_enabled_flag){
    return;
    }
    struct OVFilterBuffers fb = ctudec->filter_buffers;
    OVFrame *frame = fb.pic_frame;

    const OVPartInfo *const pinfo = ctudec->part_ctx;
    uint8_t log2_ctb_size = pinfo->log2_ctu_s;
    int ctu_width  = 1 << log2_ctb_size;

    RCNALF* alf = &alf_info->rcn_alf;
    for (int ctb_x = 0; ctb_x < nb_ctu_w; ctb_x++) {
        //TODO: change when applied on rectangular region
        int ctb_y = ctb_y_pic;
        int xPos = ctu_width * ctb_x;
        int yPos = ctu_width * ctb_y;
        int width = ( xPos + ctu_width > ctudec->pic_w ) ? ( ctudec->pic_w - xPos ) : ctu_width;
        int height = ( yPos + ctu_width > ctudec->pic_h ) ? ( ctudec->pic_h - yPos ) : ctu_width;

        //left | right | up | down
        uint8_t is_border = 0;
        is_border = (ctb_x==0)          ? is_border | OV_BOUNDARY_LEFT_RECT: is_border;
        is_border = (ctb_x==nb_ctu_w-1) ? is_border | OV_BOUNDARY_RIGHT_RECT: is_border;
        is_border = (ctb_y==0)          ? is_border | OV_BOUNDARY_UPPER_RECT: is_border;
        // is_border = (ctb_y==nb_ctu_h-1) ? is_border | OV_BOUNDARY_BOTTOM_RECT: is_border;
        is_border = (yPos + ctu_width >= ctudec->pic_h) ? is_border | OV_BOUNDARY_BOTTOM_RECT: is_border;

        int ctu_rs_addr = ctb_x + ctb_y * nb_ctu_w ;
        ALFParamsCtu alf_params_ctu = alf_info->ctb_alf_params[ctu_rs_addr];

        int16_t **src = fb.filter_region;
        ctudec_extend_filter_region(ctudec, xPos, yPos, is_border);

        if (alf_params_ctu.ctb_alf_flag & 0x4) {
            uint8_t c_idx = 0;
            Area blk_dst;
            //Source block in the filter buffers image
            int stride_src = fb.filter_region_stride[c_idx];
            int16_t*  src_luma = &src[c_idx][fb.filter_region_offset[c_idx]];

            //Destination block in the final image
            blk_dst.x=xPos; blk_dst.y=yPos;
            blk_dst.width=width; blk_dst.height=height;
            //BITDEPTH
            int stride_dst = frame->linesize[c_idx]/2;
            int16_t*  dst_luma = (int16_t*) frame->data[c_idx] + blk_dst.y*stride_dst + blk_dst.x;

            rcn_alf_derive_classification(alf, src_luma, stride_src, blk_dst, ctu_width, ctudec->pic_h);

            int16_t filter_idx = alf_params_ctu.ctb_alf_idx;
            int16_t *coeff;
            int16_t *clip;

            if (filter_idx >= NUM_FIXED_FILTER_SETS) {
                coeff = alf->coeff_aps_luma[filter_idx - NUM_FIXED_FILTER_SETS];
                clip  = alf->clip_aps_luma[filter_idx - NUM_FIXED_FILTER_SETS];
            } else {
                coeff = alf->fixed_filter_coeff_dec[filter_idx];
                clip  = alf->clip_default;
            }

            (ctudec->rcn_ctx.rcn_funcs.alf.luma)(alf->classifier, dst_luma, src_luma, stride_dst, stride_src,
                blk_dst, coeff, clip,
                ctu_width, (yPos + ctu_width >= ctudec->pic_h) ? ctudec->pic_h : height - ALF_VB_POS_ABOVE_CTUROW_LUMA);
        }

        for( uint8_t c_idx = 1; c_idx < 3; c_idx++ )
        {
            const int chr_scale = frame->linesize[0] / frame->linesize[c_idx];

            if( (c_idx==1 && (alf_params_ctu.ctb_alf_flag & 2)) || (c_idx==2 && (alf_params_ctu.ctb_alf_flag & 1)))
            {
                Area blk_dst;
                int stride_src = fb.filter_region_stride[c_idx];
                int16_t*  src_chroma = &src[c_idx][fb.filter_region_offset[c_idx]];
                //Destination block in the final image
                blk_dst.x = xPos/chr_scale;
                blk_dst.y = yPos/chr_scale;

                blk_dst.width  = width /chr_scale;
                blk_dst.height = height/chr_scale;

                int stride_dst = frame->linesize[c_idx]/2;

                int16_t*  dst_chroma = (int16_t*) frame->data[c_idx] + blk_dst.y*stride_dst + blk_dst.x;

                uint8_t alt_num = (c_idx == 1) ? alf_params_ctu.cb_alternative : alf_params_ctu.cr_alternative;

                ctudec->rcn_ctx.rcn_funcs.alf.chroma(dst_chroma, src_chroma,
                                                     stride_dst, stride_src, blk_dst,
                                                     alf->chroma_coeff_final[alt_num],
                                                     alf->chroma_clip_final[alt_num],
                                                     ctu_width/chr_scale,
                                                     (( yPos + ctu_width >= ctudec->pic_h) ? ctudec->pic_h/chr_scale
                                                      : (ctu_width - ALF_VB_POS_ABOVE_CTUROW_LUMA)/chr_scale));
            }

            if ((c_idx==1 && alf_info->cc_alf_cb_enabled_flag) || (c_idx==2 && alf_info->cc_alf_cr_enabled_flag))
            {
                const OVALFData* alf_data = (c_idx==1) ? alf_info->aps_cc_alf_data_cb : alf_info->aps_cc_alf_data_cr;

                const int filt_idx = alf_info->ctb_cc_alf_filter_idx[c_idx - 1][ctu_rs_addr];
                if (filt_idx != 0)
                {
                    //TODO: maybe reverse buffer use, the alf reconstructed pixels are in the pic frame.
                    Area blk_dst;
                    //Source block in the filter buffers image
                    int stride_src = fb.filter_region_stride[0];
                    int16_t*  src_chroma = &src[0][fb.filter_region_offset[0]];

                    //Destination block in the final image
                    blk_dst.x=xPos/chr_scale; blk_dst.y=yPos/chr_scale;
                    blk_dst.width=width/chr_scale; blk_dst.height=height/chr_scale;
                    //BITDEPTH
                    int stride_dst = frame->linesize[c_idx]/2;
                    int16_t*  dst_chroma = (int16_t*) frame->data[c_idx] + blk_dst.y*stride_dst + blk_dst.x;

                    // const int16_t *filt_coeff = alf_data.alf_cc_mapped_coeff[c_idx - 1][filt_idx];
                    const int16_t *filt_coeff = alf_data->alf_cc_mapped_coeff[c_idx - 1][filt_idx - 1];
                    cc_alf_filterBlk(dst_chroma, src_chroma, stride_dst, stride_src, blk_dst, c_idx, filt_coeff,
                    ctu_width, (( yPos + ctu_width >= ctudec->pic_h) ? ctudec->pic_h/chr_scale : (ctu_width - ALF_VB_POS_ABOVE_CTUROW_LUMA)));
                }
            }

        }
        ctudec_save_last_rows(ctudec, xPos, yPos, is_border);
        ctudec_save_last_cols(ctudec, xPos, yPos, is_border);
    }
}


void rcn_alf_destroy(RCNALF* alf, int16_t ctu_width)
{
    if( alf->classifier )
    {
        for (int i = 0; i < ctu_width>>2; i++){
            ov_free(alf->classifier[i]);
        }
        ov_free(alf->classifier);
        alf->classifier = 0;
    }
}


void rcn_init_alf_functions(struct RCNFunctions *rcn_func){
  rcn_func->alf.luma=&alf_filterBlkLuma;
  rcn_func->alf.chroma=&alf_filter_c;
}
