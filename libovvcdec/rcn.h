#ifndef RCN_H
#define RCN_H

#include "ctudec.h"

enum VVCIntraPredMode{
    VVC_PLANAR = 0,
    VVC_DC = 1,
    VVC_HOR = 18,
    VVC_DIA = 34,
    VVC_VER = 50,
    VVC_VDIA = 66,
    VVC_LM_CHROMA = 67,
    VVC_NUM_INTRA_MODES = 67,
    VVC_MDLM_LEFT = 68,
    VVC_MDLM_TOP  = 69,
    VVC_DM_CHROMA = 70,
    VVC_MIP_MODE = 75,
    VVC_NOT_AVAILABLE=128,
};

// FIXED? :
#define VVC_CTB_STRIDE (128 + 16 + 64)
#define VVC_CTB_OFFSET (VVC_CTB_STRIDE * 4 + 16)

#define VVC_CTB_STRIDE_CHROMA (128 + 64 + 16)

#define VVC_CTU_UP_FLAG             (1 << 0)
#define VVC_CTU_LEFT_FLAG           (1 << 1)
#define VVC_CTU_UPLEFT_FLAG         (1 << 2)
#define VVC_CTU_UPRIGHT_FLAG        (1 << 3)

#endif //RCN_H
