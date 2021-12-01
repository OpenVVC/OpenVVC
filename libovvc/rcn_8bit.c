#define BITDEPTH 8

#include "rcn_alf.c"
#include "rcn_df.c"
#include "rcn_intra.c"
#include "rcn_intra_cclm.c"
#include "rcn_prof_bdof.c"
#include "rcn_transform.c"
#include "rcn_transform_tree.c"

#include "rcn_fill_ref.c"
#include "rcn_intra_dc_planar.c"
#include "rcn_lmcs.c"
#include "rcn_residuals.c"

#include "rcn_dequant.c"
#include "rcn_intra_mip.c"
#include "rcn_mc.c"
#include "rcn_sao.c"

#undef BITDEPTH
