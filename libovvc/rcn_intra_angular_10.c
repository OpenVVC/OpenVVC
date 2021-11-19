#define BITDEPTH 10

#define ov_bdclip(val) ov_clip_uintp2(val, BITDEPTH);
#include "rcn_intra_angular.c"
