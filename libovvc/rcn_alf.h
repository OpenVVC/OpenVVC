#ifndef RCN_ALF_H
#define RCN_ALF_H

struct OVCTUDec;
struct RCNFunctions;
struct RectEntryInfo;

#define MAX_ALF_NUM_CLIP_VAL  4
#define NUM_BITS  8
#define CLASSIFICATION_BLK_SIZE  32  //non-normative, local buffer size
#define ALF_UNUSED_CLASSIDX  255
#define ALF_UNUSED_TRANSPOSIDX  255
#define ALF_VB_POS_ABOVE_CTUROW_LUMA   4
#define ALF_VB_POS_ABOVE_CTUROW_CHMA   2

#define MAX_NUM_ALF_CLASSES                            25
#define MAX_NUM_ALF_LUMA_COEFF                         13
#define MAX_NUM_ALF_CHROMA_COEFF                        7
#define MAX_ALF_FILTER_LENGTH                           7
#define MAX_NUM_ALF_COEFF                              MAX_ALF_FILTER_LENGTH * MAX_ALF_FILTER_LENGTH / 2 + 1
#define MAX_ALF_PADDING_SIZE                            4

#define ALF_FIXED_FILTER_NUM                           64
#define ALF_CTB_MAX_NUM_APS                             8
#define ALF_CTB_MAX_NUM_TRANSPOSE                       4
#define NUM_FIXED_FILTER_SETS                          16
#define NUM_TOTAL_FILTER_SETS                          NUM_FIXED_FILTER_SETS + ALF_CTB_MAX_NUM_APS
#define MAX_NUM_ALF_ALTERNATIVES_CHROMA                8
#define MAX_NUM_CC_ALF_FILTERS                         4
#define MAX_NUM_CC_ALF_CHROMA_COEFF                    8
#define CCALF_DYNAMIC_RANGE                            6
#define CCALF_BITS_PER_COEFF_LEVEL                     3


typedef struct Area
{
  int x,y;
  int width,height;
} Area;

enum Direction
{
  HOR,
  VER,
  DIAG0,
  DIAG1,
  NUM_DIRECTIONS
};

typedef enum
{
  CHROMA_400        = 0,
  CHROMA_420        = 1,
  CHROMA_422        = 2,
  CHROMA_444        = 3,
  NUM_CHROMA_FORMAT = 4
}ChromaFormat;


typedef enum
{
  CHANNEL_TYPE_LUMA    = 0,
  CHANNEL_TYPE_CHROMA  = 1,
  MAX_NUM_CHANNEL_TYPE = 2
}ChannelType;

typedef struct
{
  int16_t           filter_coeff_dec[NUM_FIXED_FILTER_SETS+ALF_CTB_MAX_NUM_APS][ALF_CTB_MAX_NUM_TRANSPOSE*MAX_NUM_ALF_CLASSES * MAX_NUM_ALF_LUMA_COEFF];
  int16_t           filter_clip_dec[NUM_FIXED_FILTER_SETS+ALF_CTB_MAX_NUM_APS][ALF_CTB_MAX_NUM_TRANSPOSE*MAX_NUM_ALF_CLASSES * MAX_NUM_ALF_LUMA_COEFF];
  uint8_t           created;
  uint8_t           class_idx[CLASSIFICATION_BLK_SIZE*CLASSIFICATION_BLK_SIZE];
  uint8_t           transpose_idx[CLASSIFICATION_BLK_SIZE*CLASSIFICATION_BLK_SIZE];
  int16_t           coeff_final[MAX_NUM_ALF_CLASSES * MAX_NUM_ALF_LUMA_COEFF];
  int16_t           clip_final[MAX_NUM_ALF_CLASSES * MAX_NUM_ALF_LUMA_COEFF];
  int16_t           chroma_coeff_final[MAX_NUM_ALF_ALTERNATIVES_CHROMA][MAX_NUM_ALF_CHROMA_COEFF];
  int16_t           chroma_clip_final[MAX_NUM_ALF_ALTERNATIVES_CHROMA][MAX_NUM_ALF_CHROMA_COEFF];
  int16_t           alf_clipping_values[MAX_NUM_CHANNEL_TYPE][MAX_ALF_NUM_CLIP_VAL];
}RCNALF;

void rcn_alf_create(RCNALF* alf);

void rcn_alf_filter_line(struct OVCTUDec *const ctudec, const struct RectEntryInfo *const einfo, uint16_t ctb_y);

void rcn_alf_reconstruct_coeff_APS(RCNALF* alf, struct OVCTUDec *const ctudec, uint8_t luma_flag, uint8_t chroma_flag);

void rcn_init_alf_functions(struct RCNFunctions *rcn_func);

#endif
