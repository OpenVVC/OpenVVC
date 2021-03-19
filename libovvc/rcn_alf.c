#include <stdint.h>
#include <stdlib.h>

#include "ovutils.h"
#include "ovmem.h"
#include "nvcl_structures.h"

#include "ctudec.h"
#include "rcn_alf.h"


static inline int clipALF(const int clip, const int16_t ref, const int16_t val0, const int16_t val1)
  {
    // return Clip3<int>(-clip, +clip, val0-ref) + Clip3<int>(-clip, +clip, val1-ref);
    int clip1 = (int) OVMIN( OVMAX(-clip, val0-ref) , +clip);
    int clip2 = (int) OVMIN( OVMAX(-clip, val1-ref) , +clip);
    return clip1 + clip2;
  }


//TODO: (on all functions) use the pointer to RcnALF, or directly the structure ?
void alf_create_arrays(RcnALF* alf)
{
    int tab[ALF_FIXED_FILTER_NUM][MAX_NUM_ALF_LUMA_COEFF] =
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
    for(int i=0; i<ALF_FIXED_FILTER_NUM; i++) 
    {
        for(int j=0; j<MAX_NUM_ALF_LUMA_COEFF; j++) 
            alf->fixed_filter_coeff[i][j] = tab[i][j];
    }

      //*********************************************************
    const int tab2[NUM_FIXED_FILTER_SETS][MAX_NUM_ALF_CLASSES] =
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
    for(int i=0; i<NUM_FIXED_FILTER_SETS; i++) 
    {
        for(int j=0; j<MAX_NUM_ALF_CLASSES; j++) 
            alf->class_to_filter_mapping[i][j] = tab2[i][j];
    }
}


void alf_create(OVCTUDec *const ctudec, RcnALF* alf)
{
    const OVPartInfo *const pinfo = ctudec->part_ctx;
    uint8_t log2_ctb_size = pinfo->log2_ctu_s;
    int ctu_width  = 1 << log2_ctb_size;  

    alf->classifier =  0 ;
    for (int i = 0; i < NUM_DIRECTIONS; i++)
    {
        alf->laplacian[i] = alf->laplacianPtr[i];
        for (int j = 0; j < sizeof(alf->laplacianPtr[i]) / sizeof(alf->laplacianPtr[i][0]); j++)
        {
            alf->laplacianPtr[i][j] = alf->laplacianData[i][j];
        }
    }
    // alf->numCTUsInWidth = ( ctudec->pic_w / ctu_width ) + ( ( ctudec->pic_w % ctu_width ) ? 1 : 0 );
    // alf->numCTUsInHeight = ( ctudec->pic_h / ctu_width ) + ( ( ctudec->pic_h % ctu_width ) ? 1 : 0 );
    // alf->numCTUsInPic = alf->numCTUsInHeight * alf->numCTUsInWidth;

    //BITDEPTH: uniquement pour bitdepth 10
    int bit_depth = 10; 
    alf->alfClippingValues[CHANNEL_TYPE_LUMA][0] = 1 << bit_depth;  
    int shift_luma = bit_depth - 8;
    for( int i = 1; i < MAX_ALF_NUM_CLIP_VAL; ++i )
    {
        alf->alfClippingValues[CHANNEL_TYPE_LUMA][i] = 1 << (7 - 2 * i + shift_luma);
    }

    alf->alfClippingValues[CHANNEL_TYPE_CHROMA][0] = 1 << bit_depth;
    int shift_chroma = bit_depth - 8;
    alf->alfClippingValues[CHANNEL_TYPE_CHROMA][0] = 1 << bit_depth;
    for( int i = 1; i < MAX_ALF_NUM_CLIP_VAL; ++i )
    {
        alf->alfClippingValues[CHANNEL_TYPE_CHROMA][i] = 1 << (7 - 2 * i + shift_chroma);
    }

    alf_create_arrays(alf);

    // Classification
    if ( alf->classifier == 0 )
    {
        alf->classifier = ov_malloc(sizeof(ALFClassifier*)*ctu_width);
        for (int i = 0; i < ctu_width; i++)
        {
            alf->classifier[i] = ov_malloc(sizeof(ALFClassifier)*ctu_width );
        }
    }

    for (int filter_set_idx = 0; filter_set_idx < NUM_FIXED_FILTER_SETS; filter_set_idx++)
    {
        for (int class_idx = 0; class_idx < MAX_NUM_ALF_CLASSES; class_idx++)
        {
            int fixed_filter_idx = alf->class_to_filter_mapping[filter_set_idx][class_idx];
            for (int i = 0; i < MAX_NUM_ALF_LUMA_COEFF - 1; i++)
            {
                alf->fixedFilterSetCoeffDec[filter_set_idx][class_idx * MAX_NUM_ALF_LUMA_COEFF + i] = alf->fixed_filter_coeff[fixed_filter_idx][i];
            }
            alf->fixedFilterSetCoeffDec[filter_set_idx][class_idx * MAX_NUM_ALF_LUMA_COEFF + MAX_NUM_ALF_LUMA_COEFF - 1] = (1 << (NUM_BITS - 1));
        }
    }
    for (int i = 0; i < MAX_NUM_ALF_LUMA_COEFF * MAX_NUM_ALF_CLASSES; i++)
    {
        alf->clipDefault[i] = alf->alfClippingValues[CHANNEL_TYPE_LUMA][0];
    }
}


void alf_reconstructCoeff(RcnALF* alf, const struct OVALFData* alf_data, ChannelType channel)
{
    int factor = 1 << (NUM_BITS - 1);
    int num_classes = ( channel == CHANNEL_TYPE_LUMA ) ? MAX_NUM_ALF_CLASSES : 1;
    int num_coeff = ( channel == CHANNEL_TYPE_CHROMA ) ? 7 : 13;
    int num_coeff_minus1 = num_coeff - 1;
    const int num_alts = ( channel == CHANNEL_TYPE_LUMA ) ? 1 : alf_data->alf_chroma_num_alt_filters_minus1 + 1;
    
    int numFilters = ( channel == CHANNEL_TYPE_LUMA ) ? alf_data->alf_luma_num_filters_signalled_minus1 + 1 : 1;
    int16_t* coeff ;
    int16_t* clipp;

    for( int alt_idx = 0; alt_idx < num_alts; ++ alt_idx )
    {
        if( channel == CHANNEL_TYPE_CHROMA )
        {
            coeff = alf_data->alf_chroma_coeff[alt_idx];
            clipp = alf_data->alf_chroma_clip_idx[alt_idx];
            for( int coeffIdx = 0; coeffIdx < num_coeff_minus1; ++coeffIdx )
            {
                alf->chroma_coeffFinal[alt_idx][coeffIdx] = coeff[coeffIdx];
                int clipIdx = alf_data->alf_chroma_clip_flag ? clipp[coeffIdx] : 0;
                alf->chroma_clippFinal[alt_idx][coeffIdx] = alf->alfClippingValues[channel][clipIdx];
            }
            alf->chroma_coeffFinal[alt_idx][num_coeff_minus1] = factor;
            alf->chroma_clippFinal[alt_idx][num_coeff_minus1] = alf->alfClippingValues[channel][0];
            continue;
        }
        else{
            coeff = alf_data->alf_luma_coeff;
            clipp = alf_data->alf_luma_clip_idx;
            for( int filter_idx = 0; filter_idx < numFilters; filter_idx++ )
            {
                coeff[filter_idx* MAX_NUM_ALF_LUMA_COEFF + num_coeff_minus1] = factor;
            }
        }

        for( int class_idx = 0; class_idx < num_classes; class_idx++ )
        {
            int filter_idx = alf_data->alf_luma_coeff_delta_idx[class_idx];
            for (int coeffIdx = 0; coeffIdx < num_coeff_minus1; ++coeffIdx)
            {
                alf->coeffFinal[class_idx * MAX_NUM_ALF_LUMA_COEFF + coeffIdx] = coeff[filter_idx * MAX_NUM_ALF_LUMA_COEFF + coeffIdx];
            }
            alf->coeffFinal[class_idx* MAX_NUM_ALF_LUMA_COEFF + num_coeff_minus1] = factor;
            alf->clippFinal[class_idx* MAX_NUM_ALF_LUMA_COEFF + num_coeff_minus1] = alf->alfClippingValues[channel][0];
            for( int coeffIdx = 0; coeffIdx < num_coeff_minus1; ++coeffIdx )
            {
                int clipIdx = alf_data->alf_chroma_clip_flag ? (clipp + filter_idx * MAX_NUM_ALF_LUMA_COEFF)[coeffIdx] : 0;
                (alf->clippFinal + class_idx * MAX_NUM_ALF_LUMA_COEFF)[coeffIdx] = alf->alfClippingValues[channel][clipIdx];
            }
            alf->clippFinal[class_idx* MAX_NUM_ALF_LUMA_COEFF + num_coeff_minus1] = alf->alfClippingValues[channel][0];
        }
    }
}

void alf_reconstruct_coeff_APS(RcnALF* alf, OVCTUDec *const ctudec, uint8_t luma_flag, uint8_t chroma_flag)
{
    const struct OVALFData* alf_data = ctudec->alf_info.aps_alf_data;

    if (luma_flag){
        // for (int i = 0; i < sh.slice_num_alf_aps_ids_luma; i++)
        // {
        //     //ATTENTION: pourquoi pas un tableau comme VTM -> parce que toujours une seule boucle
        //     int apsIdx = sh.slice_alf_aps_id_luma;
        //     int ret = ff_vvc_activate_alf(vvc_ctx, apsIdx);

            alf_reconstructCoeff(alf, alf_data, CHANNEL_TYPE_LUMA);

            for(int j = 0; j < MAX_NUM_ALF_CLASSES * MAX_NUM_ALF_LUMA_COEFF; j++){
                //TODO: m_coeffApsLuma may be only 1 dimensional
                alf->coeffApsLuma[0][j] = alf->coeffFinal[j];
                alf->clippApsLuma[0][j] = alf->clippFinal[j];
            }
        // }
    }

    if (chroma_flag){
        // int apsIdx = sh.slice_alf_aps_id_chroma;
        // int ret = ff_vvc_activate_alf(vvc_ctx, apsIdx);
        alf_reconstructCoeff(alf, alf_data, CHANNEL_TYPE_CHROMA);
    }
}



void alf_deriveClassificationBlk(ALFClassifier **classifier, int **laplacian[NUM_DIRECTIONS],
                                                 int16_t *const src, const int stride, const Area blk,
                                                 const int shift, const int vbCTUHeight, int vbPos)
{
 
  static const int th[16] = { 0, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 4 };
  const int maxActivity = 15;

  int fl = 2;
  int flP1 = fl + 1;
  int fl2 = 2 * fl;

  int mainDirection, secondaryDirection, dirTempHV, dirTempD;

  int pixY;
  int height = blk.height + fl2;
  int width = blk.width + fl2;

  for( int i = 0; i < height; i += 2 )
  {
    int yoffset = ( i + 1 - flP1 ) * stride - flP1;
    const int16_t *src0 = &src[yoffset - stride];
    const int16_t *src1 = &src[yoffset];
    const int16_t *src2 = &src[yoffset + stride];
    const int16_t *src3 = &src[yoffset + stride * 2];

    const int y = blk.y - 2 + i;
    if (y > 0 && (y & (vbCTUHeight - 1)) == vbPos - 2)
    {
      src3 = &src[yoffset + stride];
    }
    else if (y > 0 && (y & (vbCTUHeight - 1)) == vbPos)
    {
      src0 = &src[yoffset];
    }
    int* pYver = laplacian[VER][i];
    int* pYhor = laplacian[HOR][i];
    int* pYdig0 = laplacian[DIAG0][i];
    int* pYdig1 = laplacian[DIAG1][i];

    for( int j = 0; j < width; j += 2 )
    {
      // pixY = j + 1 + posX;
      pixY = j + 1 ;
      const int16_t *pY = src1 + pixY;
      const int16_t* pYdown = src0 + pixY;
      const int16_t* pYup = src2 + pixY;
      const int16_t* pYup2 = src3 + pixY;

      const int16_t y0 = pY[0] << 1;
      const int16_t yup1 = pYup[1] << 1;

      //modification des buffers laplacian ici 
      pYver[j] = abs( y0 - pYdown[0] - pYup[0] ) + abs( yup1 - pY[1] - pYup2[1] );
      pYhor[j] = abs( y0 - pY[1] - pY[-1] ) + abs( yup1 - pYup[2] - pYup[0] );
      pYdig0[j] = abs( y0 - pYdown[-1] - pYup[1] ) + abs( yup1 - pY[0] - pYup2[2] );
      pYdig1[j] = abs( y0 - pYup[-1] - pYdown[1] ) + abs( yup1 - pYup2[0] - pY[2] );

      if( j > 4 && ( j - 6 ) % 4 == 0 )
      {
        int jM6 = j - 6;
        int jM4 = j - 4;
        int jM2 = j - 2;

        pYver[jM6] += pYver[jM4] + pYver[jM2] + pYver[j];
        pYhor[jM6] += pYhor[jM4] + pYhor[jM2] + pYhor[j];
        pYdig0[jM6] += pYdig0[jM4] + pYdig0[jM2] + pYdig0[j];
        pYdig1[jM6] += pYdig1[jM4] + pYdig1[jM2] + pYdig1[j];
      }
    }
  }

  // classification block size
  const int clsSizeY = 4;
  const int clsSizeX = 4;

  for( int i = 0; i < blk.height; i += clsSizeY )
  {
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

    for( int j = 0; j < blk.width; j += clsSizeX )
    {
      int sumV = 0; int sumH = 0; int sumD0 = 0; int sumD1 = 0;
      if (((i + blk.y) % vbCTUHeight) == (vbPos - 4))
      {
        sumV = pYver[j] + pYver2[j] + pYver4[j];
        sumH = pYhor[j] + pYhor2[j] + pYhor4[j];
        sumD0 = pYdig0[j] + pYdig02[j] + pYdig04[j];
        sumD1 = pYdig1[j] + pYdig12[j] + pYdig14[j];
      }
      else if (((i + blk.y) % vbCTUHeight) == vbPos)
      {
        sumV = pYver2[j] + pYver4[j] + pYver6[j];
        sumH = pYhor2[j] + pYhor4[j] + pYhor6[j];
        sumD0 = pYdig02[j] + pYdig04[j] + pYdig06[j];
        sumD1 = pYdig12[j] + pYdig14[j] + pYdig16[j];
      }
      else
      {
        sumV = pYver[j] + pYver2[j] + pYver4[j] + pYver6[j];
        sumH = pYhor[j] + pYhor2[j] + pYhor4[j] + pYhor6[j];
        sumD0 = pYdig0[j] + pYdig02[j] + pYdig04[j] + pYdig06[j];
        sumD1 = pYdig1[j] + pYdig12[j] + pYdig14[j] + pYdig16[j];
      }

      int tempAct = sumV + sumH;
      int activity = 0;

      const int y = (i + blk.y) & (vbCTUHeight - 1);
      if (y == vbPos - 4 || y == vbPos)
      {
        activity = (int16_t) OVMIN( OVMAX(0, (tempAct * 96) >> shift) , maxActivity);
      }
      else
      {
        activity = (int16_t) OVMIN( OVMAX(0, (tempAct * 64) >> shift) , maxActivity);
      }
      int class_idx = th[activity];

      int hv1, hv0, d1, d0, hvd1, hvd0;

      if( sumV > sumH )
      {
        hv1 = sumV;
        hv0 = sumH;
        dirTempHV = 1;
      }
      else
      {
        hv1 = sumH;
        hv0 = sumV;
        dirTempHV = 3;
      }
      if( sumD0 > sumD1 )
      {
        d1 = sumD0;
        d0 = sumD1;
        dirTempD = 0;
      }
      else
      {
        d1 = sumD1;
        d0 = sumD0;
        dirTempD = 2;
      }
      if( (uint32_t)d1 * (uint32_t)hv0 > (uint32_t)hv1 * (uint32_t)d0 )
      {
        hvd1 = d1;
        hvd0 = d0;
        mainDirection = dirTempD;
        secondaryDirection = dirTempHV;
      }
      else
      {
        hvd1 = hv1;
        hvd0 = hv0;
        mainDirection = dirTempHV;
        secondaryDirection = dirTempD;
      }

      int directionStrength = 0;
      if( hvd1 > 2 * hvd0 )
      {
        directionStrength = 1;
      }
      if( hvd1 * 2 > 9 * hvd0 )
      {
        directionStrength = 2;
      }

      if( directionStrength )
      {
        class_idx += ( ( ( mainDirection & 0x1 ) << 1 ) + directionStrength ) * 5;
      }

      static const int transposeTable[8] = { 0, 1, 0, 2, 2, 3, 1, 3 };
      int transpose_idx = transposeTable[mainDirection * 2 + ( secondaryDirection >> 1 )];

      int yOffset = (i + blk.y) % vbCTUHeight;
      int xOffset = (j + blk.x) % vbCTUHeight;

      //TODO: divide by 16 the size of alf.classifier
      ALFClassifier *cl0 = &classifier[yOffset][xOffset];
      ALFClassifier *cl1 = &classifier[yOffset + 1][xOffset];
      ALFClassifier *cl2 = &classifier[yOffset + 2][xOffset];
      ALFClassifier *cl3 = &classifier[yOffset + 3][xOffset];
      ALFClassifier alf_class;
      alf_class.class_idx = class_idx;
      alf_class.transpose_idx = transpose_idx;
      cl0[0] = alf_class;
      cl0[1] = alf_class;
      cl0[2] = alf_class;
      cl0[3] = alf_class;
      cl1[0] = alf_class;
      cl1[1] = alf_class;
      cl1[2] = alf_class;
      cl1[3] = alf_class;
      cl2[0] = alf_class;
      cl2[1] = alf_class;
      cl2[2] = alf_class;
      cl2[3] = alf_class;
      cl3[0] = alf_class;
      cl3[1] = alf_class;
      cl3[2] = alf_class;
      cl3[3] = alf_class;
    }
  }
}


void alf_deriveClassification(RcnALF alf, int16_t *const rcn_img, const int stride, Area blk, int ctu_width, int pic_h )
{
    ALFClassifier** classifier = alf.classifier;  
    int height = blk.y + blk.height;
    int width = blk.x + blk.width;
    //BITDEPTH: uniquement pour bitdepth 10
    int bit_depth = 10;

    for( int i = blk.y; i < height; i += CLASSIFICATION_BLK_SIZE )
    // for( int i = 0; i < blk.height; i += CLASSIFICATION_BLK_SIZE )
    {
        //TODO: handle alf.alfVBLumaPos
        int nHeight = OVMIN( i + CLASSIFICATION_BLK_SIZE, height ) - i;
        // int nHeight = OVMIN( i + CLASSIFICATION_BLK_SIZE, blk.height ) - i;

        for( int j = blk.x; j < width; j += CLASSIFICATION_BLK_SIZE )
        // for( int j = 0; j < blk.width; j += CLASSIFICATION_BLK_SIZE )
        {
            int nWidth = OVMIN( j + CLASSIFICATION_BLK_SIZE, width ) - j;
            // int nWidth = OVMIN( j + CLASSIFICATION_BLK_SIZE, blk.width ) - j;
            Area blk_class;
            blk_class.x = j; blk_class.y = i;
            blk_class.width = nWidth; blk_class.height = nHeight; 

            int16_t* rcn_img_class = rcn_img + (i-blk.y)*stride + (j-blk.x);
            alf_deriveClassificationBlk(classifier, alf.laplacian, rcn_img_class, stride, blk_class, bit_depth + 4
            , ctu_width, (blk.height<ctu_width) ? pic_h : blk.height - ALF_VB_POS_ABOVE_CTUROW_LUMA
            );
        }
    }
}


// dst   : buffer Frame output, pointing to the begining of CTU.
// src   : filter buffer pre-ALF (of size CTU)
// blkDst: location and dimension of destination block in frame  
// blk   : location and dimension of destination block in filter buffer  
void alf_filterBlk(ALFClassifier **classifier, int16_t *const dst, int16_t *const src, const int dstStride, const int srcStride, 
                         Area blkDst, Area blk, const ComponentID compId,
                         const int16_t *filterSet, const int16_t *fClipSet,
                         const int vbCTUHeight, int vbPos)
{
  const int16_t *pImgYPad0, *pImgYPad1, *pImgYPad2, *pImgYPad3, *pImgYPad4, *pImgYPad5, *pImgYPad6;
  const int16_t *pImg0, *pImg1, *pImg2, *pImg3, *pImg4, *pImg5, *pImg6;

  const int16_t *coef = filterSet;
  const int16_t *clip = fClipSet;

  const int shift = NUM_BITS - 1;
  const int offset = 1 << ( shift - 1 );

  int transpose_idx = 0;
  const int clsSizeY = 4;
  const int clsSizeX = 4;
  ALFClassifier *pClass = 0;

  int dstStride2 = dstStride * clsSizeY;
  int srcStride2 = srcStride * clsSizeY;

  pImgYPad0 = src ;
  pImgYPad1 = pImgYPad0 + srcStride;
  pImgYPad2 = pImgYPad0 - srcStride;
  pImgYPad3 = pImgYPad1 + srcStride;
  pImgYPad4 = pImgYPad2 - srcStride;
  pImgYPad5 = pImgYPad3 + srcStride;
  pImgYPad6 = pImgYPad4 - srcStride;


  //Partie a changer pour enlever les mallocs
  //creer une shuffleTab comme dans alf_sse
  int16_t* filterCoeff;
  int16_t* filterClipp;
    const uint8_t bChroma = compId > 0;
  if( !bChroma )
  {
    filterCoeff = ov_malloc(MAX_NUM_ALF_LUMA_COEFF*sizeof(int));
    filterClipp = ov_malloc(MAX_NUM_ALF_LUMA_COEFF*sizeof(int));
  }
  else
  {
    filterCoeff = ov_malloc(MAX_NUM_ALF_CHROMA_COEFF*sizeof(int));
    filterClipp = ov_malloc(MAX_NUM_ALF_CHROMA_COEFF*sizeof(int));
  }

  int16_t* pRec0 = dst ;
  int16_t* pRec1 = pRec0 + dstStride;

  for( int i = 0; i < blkDst.height; i += clsSizeY )
  {
    if( !bChroma )
    {
      pClass = classifier[i];
    }

    for( int j = 0; j < blkDst.width; j += clsSizeX )
    {
      if( !bChroma )
      {
        ALFClassifier cl = pClass[j];
        transpose_idx = cl.transpose_idx;
        coef = filterSet + cl.class_idx * MAX_NUM_ALF_LUMA_COEFF;
        clip = fClipSet + cl.class_idx * MAX_NUM_ALF_LUMA_COEFF;
        // printf("%i %i : %i %i\n", j, i, cl.class_idx, cl.transpose_idx);

        if( transpose_idx == 1 )
        {
          int filterCoeff_tmp[MAX_NUM_ALF_LUMA_COEFF] = { coef[9], coef[4], coef[10], coef[8], coef[1], coef[5], coef[11], coef[7], coef[3], coef[0], coef[2], coef[6], coef[12] };
          int filterClipp_tmp[MAX_NUM_ALF_LUMA_COEFF] = { clip[9], clip[4], clip[10], clip[8], clip[1], clip[5], clip[11], clip[7], clip[3], clip[0], clip[2], clip[6], clip[12] };
          for (int i=0; i<MAX_NUM_ALF_LUMA_COEFF; i++)
          {
            filterCoeff[i] = filterCoeff_tmp[i];
            filterClipp[i] = filterClipp_tmp[i];
          }

        }
        else if( transpose_idx == 2 )
        {
          int filterCoeff_tmp[MAX_NUM_ALF_LUMA_COEFF] = { coef[0], coef[3], coef[2], coef[1], coef[8], coef[7], coef[6], coef[5], coef[4], coef[9], coef[10], coef[11], coef[12] };
          int filterClipp_tmp[MAX_NUM_ALF_LUMA_COEFF] = { clip[0], clip[3], clip[2], clip[1], clip[8], clip[7], clip[6], clip[5], clip[4], clip[9], clip[10], clip[11], clip[12] };
          for (int i=0; i<MAX_NUM_ALF_LUMA_COEFF; i++)
          {
            filterCoeff[i] = filterCoeff_tmp[i];
            filterClipp[i] = filterClipp_tmp[i];
          }
        }
        else if( transpose_idx == 3 )
        {
          int filterCoeff_tmp[MAX_NUM_ALF_LUMA_COEFF] = { coef[9], coef[8], coef[10], coef[4], coef[3], coef[7], coef[11], coef[5], coef[1], coef[0], coef[2], coef[6], coef[12] };
          int filterClipp_tmp[MAX_NUM_ALF_LUMA_COEFF] = { clip[9], clip[8], clip[10], clip[4], clip[3], clip[7], clip[11], clip[5], clip[1], clip[0], clip[2], clip[6], clip[12] };
          for (int i=0; i<MAX_NUM_ALF_LUMA_COEFF; i++)
          {
            filterCoeff[i] = filterCoeff_tmp[i];
            filterClipp[i] = filterClipp_tmp[i];
          }
        }
        else
        {
          int filterCoeff_tmp[MAX_NUM_ALF_LUMA_COEFF] = { coef[0], coef[1], coef[2], coef[3], coef[4], coef[5], coef[6], coef[7], coef[8], coef[9], coef[10], coef[11], coef[12] };
          int filterClipp_tmp[MAX_NUM_ALF_LUMA_COEFF] = { clip[0], clip[1], clip[2], clip[3], clip[4], clip[5], clip[6], clip[7], clip[8], clip[9], clip[10], clip[11], clip[12] };
          for (int i=0; i<MAX_NUM_ALF_LUMA_COEFF; i++)
          {
            filterCoeff[i] = filterCoeff_tmp[i];
            filterClipp[i] = filterClipp_tmp[i];
          }
        }
      }
      else
      {
        if( transpose_idx == 1 )
        {
          int filterCoeff_tmp[MAX_NUM_ALF_CHROMA_COEFF] = { coef[4], coef[1], coef[5], coef[3], coef[0], coef[2], coef[6] };
          int filterClipp_tmp[MAX_NUM_ALF_CHROMA_COEFF] = { clip[4], clip[1], clip[5], clip[3], clip[0], clip[2], clip[6] };
          for (int i=0; i<MAX_NUM_ALF_CHROMA_COEFF; i++)
          {
            filterCoeff[i] = filterCoeff_tmp[i];
            filterClipp[i] = filterClipp_tmp[i];
          }
        }
        else if( transpose_idx == 2 )
        {
          int filterCoeff_tmp[MAX_NUM_ALF_CHROMA_COEFF] = { coef[0], coef[3], coef[2], coef[1], coef[4], coef[5], coef[6] };
          int filterClipp_tmp[MAX_NUM_ALF_CHROMA_COEFF] = { clip[0], clip[3], clip[2], clip[1], clip[4], clip[5], clip[6] };
          for (int i=0; i<MAX_NUM_ALF_CHROMA_COEFF; i++)
          {
            filterCoeff[i] = filterCoeff_tmp[i];
            filterClipp[i] = filterClipp_tmp[i];
          }
        }
        else if( transpose_idx == 3 )
        {
          int filterCoeff_tmp[MAX_NUM_ALF_CHROMA_COEFF] = { coef[4], coef[3], coef[5], coef[1], coef[0], coef[2], coef[6] };
          int filterClipp_tmp[MAX_NUM_ALF_CHROMA_COEFF] = { clip[4], clip[3], clip[5], clip[1], clip[0], clip[2], clip[6] };
          for (int i=0; i<MAX_NUM_ALF_CHROMA_COEFF; i++)
          {
            filterCoeff[i] = filterCoeff_tmp[i];
            filterClipp[i] = filterClipp_tmp[i];
          }
        }
        else
        {
          int filterCoeff_tmp[MAX_NUM_ALF_CHROMA_COEFF] = { coef[0], coef[1], coef[2], coef[3], coef[4], coef[5], coef[6] };
          int filterClipp_tmp[MAX_NUM_ALF_CHROMA_COEFF] = { clip[0], clip[1], clip[2], clip[3], clip[4], clip[5], clip[6] };
          for (int i=0; i<MAX_NUM_ALF_CHROMA_COEFF; i++)
          {
            filterCoeff[i] = filterCoeff_tmp[i];
            filterClipp[i] = filterClipp_tmp[i];
          }
        }
      }

      for( int ii = 0; ii < clsSizeY; ii++ )
      {
        pImg0 = pImgYPad0 + j + ii * srcStride;
        pImg1 = pImgYPad1 + j + ii * srcStride;
        pImg2 = pImgYPad2 + j + ii * srcStride;
        pImg3 = pImgYPad3 + j + ii * srcStride;
        pImg4 = pImgYPad4 + j + ii * srcStride;
        pImg5 = pImgYPad5 + j + ii * srcStride;
        pImg6 = pImgYPad6 + j + ii * srcStride;

        pRec1 = pRec0 + j + ii * dstStride;

        int yVb = (blkDst.y + i + ii) & (vbCTUHeight - 1);
        if (yVb < vbPos && (yVb >= vbPos - (bChroma ? 2 : 4)))   // above
        {
          pImg1 = (yVb == vbPos - 1) ? pImg0 : pImg1;
          pImg3 = (yVb >= vbPos - 2) ? pImg1 : pImg3;
          pImg5 = (yVb >= vbPos - 3) ? pImg3 : pImg5;

          pImg2 = (yVb == vbPos - 1) ? pImg0 : pImg2;
          pImg4 = (yVb >= vbPos - 2) ? pImg2 : pImg4;
          pImg6 = (yVb >= vbPos - 3) ? pImg4 : pImg6;
        }
        else if (yVb >= vbPos && (yVb <= vbPos + (bChroma ? 1 : 3)))   // bottom
        {
          pImg2 = (yVb == vbPos) ? pImg0 : pImg2;
          pImg4 = (yVb <= vbPos + 1) ? pImg2 : pImg4;
          pImg6 = (yVb <= vbPos + 2) ? pImg4 : pImg6;

          pImg1 = (yVb == vbPos) ? pImg0 : pImg1;
          pImg3 = (yVb <= vbPos + 1) ? pImg1 : pImg3;
          pImg5 = (yVb <= vbPos + 2) ? pImg3 : pImg5;
        }


        uint8_t isNearVBabove = yVb < vbPos && (yVb >= vbPos - 1);
        uint8_t isNearVBbelow = yVb >= vbPos && (yVb <= vbPos);
        for( int jj = 0; jj < clsSizeX; jj++ )
        {
          int sum = 0;
          const int16_t curr = pImg0[+0];
          if( !bChroma )
          {
            sum += filterCoeff[0] * ( clipALF(filterClipp[0], curr, pImg5[+0], pImg6[+0]) );

            sum += filterCoeff[1] * ( clipALF(filterClipp[1], curr, pImg3[+1], pImg4[-1]) );
            sum += filterCoeff[2] * ( clipALF(filterClipp[2], curr, pImg3[+0], pImg4[+0]) );
            sum += filterCoeff[3] * ( clipALF(filterClipp[3], curr, pImg3[-1], pImg4[+1]) );

            sum += filterCoeff[4] * ( clipALF(filterClipp[4], curr, pImg1[+2], pImg2[-2]) );
            sum += filterCoeff[5] * ( clipALF(filterClipp[5], curr, pImg1[+1], pImg2[-1]) );
            sum += filterCoeff[6] * ( clipALF(filterClipp[6], curr, pImg1[+0], pImg2[+0]) );
            sum += filterCoeff[7] * ( clipALF(filterClipp[7], curr, pImg1[-1], pImg2[+1]) );
            sum += filterCoeff[8] * ( clipALF(filterClipp[8], curr, pImg1[-2], pImg2[+2]) );

            sum += filterCoeff[9] * ( clipALF(filterClipp[9], curr, pImg0[+3], pImg0[-3]) );
            sum += filterCoeff[10] * ( clipALF(filterClipp[10], curr, pImg0[+2], pImg0[-2]) );
            sum += filterCoeff[11] * ( clipALF(filterClipp[11], curr, pImg0[+1], pImg0[-1]) );
          }
          else
          {
            sum += filterCoeff[0] * ( clipALF(filterClipp[0], curr, pImg3[+0], pImg4[+0]) );

            sum += filterCoeff[1] * ( clipALF(filterClipp[1], curr, pImg1[+1], pImg2[-1]) );
            sum += filterCoeff[2] * ( clipALF(filterClipp[2], curr, pImg1[+0], pImg2[+0]) );
            sum += filterCoeff[3] * ( clipALF(filterClipp[3], curr, pImg1[-1], pImg2[+1]) );

            sum += filterCoeff[4] * ( clipALF(filterClipp[4], curr, pImg0[+2], pImg0[-2]) );
            sum += filterCoeff[5] * ( clipALF(filterClipp[5], curr, pImg0[+1], pImg0[-1]) );
          }

          if (!(isNearVBabove || isNearVBbelow))
          {
            sum = (sum + offset) >> shift;
          }
          else
          {
            //Rounding offset fix
            sum = (sum + (1 << ((shift + 3) - 1))) >> (shift + 3);
          }

          sum += curr;
          //BITDEPTH: pour etre generique utiliser VVCSPSData: bit_depth_luma_minus8; bit_depth_chroma_minus8;
          // pRec1[jj] = ClipPel( sum, clpRng );
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

    pImgYPad0 += srcStride2;
    pImgYPad1 += srcStride2;
    pImgYPad2 += srcStride2;
    pImgYPad3 += srcStride2;
    pImgYPad4 += srcStride2;
    pImgYPad5 += srcStride2;
    pImgYPad6 += srcStride2;
  }
  ov_free(filterCoeff);
  ov_free(filterClipp);
}


void rcn_alf_filter_line(OVCTUDec *const ctudec, int nb_ctu_w, uint16_t ctb_y_pic)
{
    struct ALFInfo alf_info = ctudec->alf_info;
    if (!alf_info.alf_luma_enabled_flag && !alf_info.alf_cb_enabled_flag && !alf_info.alf_cr_enabled_flag){
    return;
    }
    struct OVFilterBuffers fb = ctudec->filter_buffers;
    OVFrame *frame = fb.pic_frame;

    const OVPartInfo *const pinfo = ctudec->part_ctx;
    uint8_t log2_ctb_size = pinfo->log2_ctu_s;
    int ctu_width  = 1 << log2_ctb_size;

    //BITDEPTH: uniquement pour bitdepth 10
    int int16_t_shift = 1;

    //Initialization of ALF reconstruction structures
    RcnALF alf;
    alf_create(ctudec, &alf);

    uint8_t luma_flag = 1;
    uint8_t chroma_flag = alf_info.alf_cb_enabled_flag || alf_info.alf_cr_enabled_flag;
    alf_reconstruct_coeff_APS(&alf, ctudec, luma_flag, chroma_flag);

    uint8_t clipTop = 0, clipBottom = 0, clipLeft = 0, clipRight = 0;
    int numHorVirBndry = 0, numVerVirBndry = 0;
    int horVirBndryPos[] = { 0, 0, 0 };
    int verVirBndryPos[] = { 0, 0, 0 };

    for (int ctb_x = 0; ctb_x < nb_ctu_w; ctb_x++) {
        //TODO: change when applied on rectangular region
        int ctb_y = ctb_y_pic;

        //left | right | up | down
        uint8_t is_border = 0; 
        // is_border = (ctb_x==0)          ? is_border | VVC_BOUNDARY_LEFT_TILE: is_border;
        // is_border = (ctb_x==nb_ctu_w-1) ? is_border | VVC_BOUNDARY_RIGHT_TILE: is_border;
        // is_border = (ctb_y==0)          ? is_border | VVC_BOUNDARY_UPPER_TILE: is_border;
        // is_border = (ctb_y==nb_ctu_h-1) ? is_border | VVC_BOUNDARY_BOTTOM_TILE: is_border;
        int ctu_rs_addr = ctb_x + ctb_y * nb_ctu_w ;
        ALFParamsCtu alf_params_ctu = ctudec->alf_info.alf_params[ctu_rs_addr];

        int xPos = ctu_width * ctb_x;
        int yPos = ctu_width * ctb_y;

        int width = ( xPos + ctu_width > ctudec->pic_w ) ? ( ctudec->pic_w - xPos ) : ctu_width;
        int height = ( yPos + ctu_width > ctudec->pic_h ) ? ( ctudec->pic_h - yPos ) : ctu_width;

        // extend_ctu_filter_buffer(frame, lc_ctx, tile_ctx, tile_idx,
        //       ctb_x, ctb_y, ctu_width, vvc_ctx->margin, is_border);

        //TODO: int16_t** or uint8_t**??
        int margin = fb.margin;
        int16_t **src = fb.filter_region;

        if( alf_params_ctu.ctb_alf_flag ){
            Area blk,blk_dst;
            //Source block in the filter buffers image
            // blk.x=0; blk.y=0;
            blk.x=xPos+margin; blk.y=yPos+margin;
            blk.width=width; blk.height=height;
            int stride_src = fb.filter_region_stride[0];
            int16_t*  src_luma = &src[0][blk.y*stride_src + blk.x];

            //Destination block in the final image
            blk_dst.x=xPos; blk_dst.y=yPos;
            blk_dst.width=width; blk_dst.height=height;
            //BITDEPTH
            int stride_dst = frame->linesize[0]/2;
            int16_t*  dst_luma = (int16_t*) frame->data[0] + blk_dst.y*stride_dst + blk_dst.x;
            // alf_deriveClassification( alf, dst_luma, stride_dst, blk_dst );
            alf_deriveClassification( alf, src_luma, stride_src, blk_dst, ctu_width, ctudec->pic_h);

            int16_t filterSetIndex = alf_params_ctu.ctb_alf_idx;
            int16_t *coeff;
            int16_t *clip;
            if (filterSetIndex >= NUM_FIXED_FILTER_SETS){
                coeff = alf.coeffApsLuma[filterSetIndex - NUM_FIXED_FILTER_SETS];
                clip = alf.clippApsLuma[filterSetIndex - NUM_FIXED_FILTER_SETS];
            }
            else{
                coeff = alf.fixedFilterSetCoeffDec[filterSetIndex];
                clip = alf.clipDefault;
            }

            alf_filterBlk(alf.classifier, dst_luma, src_luma, stride_dst, stride_src,
                blk_dst, blk, COMPONENT_Y, coeff, clip, 
                ctu_width, (yPos + ctu_width >= ctudec->pic_h) ? ctudec->pic_h : blk.height - ALF_VB_POS_ABOVE_CTUROW_LUMA);
        }

        for( int compIdx = 1; compIdx < MAX_NUM_COMPONENT; compIdx++ )
        {
            ComponentID compID = (ComponentID) compIdx ;
            const int chr_scale = frame->linesize[0] / frame->linesize[compID];

            if( (compID==1 && (alf_params_ctu.ctb_alf_flag & 2)) || (compID==2 && (alf_params_ctu.ctb_alf_flag & 1)))
            {
                Area blk,blk_dst;
                //Source block in the filter buffers image
                // blk.x=0; blk.y=0;
                blk.x=xPos/chr_scale+margin; blk.y=yPos/chr_scale+margin;
                blk.width=width/chr_scale; blk.height=height/chr_scale;
                int stride_src = fb.filter_region_stride[compID];
                int16_t*  src_chroma = &src[compID][blk.y*stride_src + blk.x];

                //Destination block in the final image
                blk_dst.x=xPos/chr_scale; blk_dst.y=yPos/chr_scale;
                blk_dst.width=width/chr_scale; blk_dst.height=height/chr_scale;
                //BITDEPTH
                int stride_dst = frame->linesize[compID]/2;
                int16_t*  dst_chroma = (int16_t*) frame->data[compID] + blk_dst.y*stride_dst + blk_dst.x;

                uint8_t alt_num = (compID == COMPONENT_Cb) ? alf_params_ctu.cb_alternative : alf_params_ctu.cr_alternative;

                alf_filterBlk(alf.classifier, dst_chroma, src_chroma, stride_dst, stride_src, blk_dst, blk, compID,
                    alf.chroma_coeffFinal[alt_num], alf.chroma_clippFinal[alt_num], 
                    ctu_width/chr_scale, (( yPos + ctu_width >= ctudec->pic_h) ? ctudec->pic_h/chr_scale : (ctu_width - ALF_VB_POS_ABOVE_CTUROW_LUMA)/chr_scale));
            }
        }
        // save_last_rows_ctu(lc_ctx, nb_ctu_w, ctb_x, ctu_width, vvc_ctx->margin, is_border);

        // //func save_last_cols_ctu
        // save_last_cols_ctu(lc_ctx, ctu_width, vvc_ctx->margin, is_border);
    }
}


void alf_destroy(RcnALF* alf, int16_t ctu_width)
{
    if( alf->classifier )
    {
        for (int i = 0; i < ctu_width; i++)
            ov_free(alf->classifier[i]);
        ov_free(alf->classifier);
        alf->classifier = 0;
    }
}

