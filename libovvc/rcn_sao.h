#ifndef RCN_SAO_H
#define RCN_SAO_H

#include "ctudec.h"

void rcn_sao_filter_line(OVCTUDec *const ctudec, const struct RectEntryInfo *const einfo, uint16_t ctb_y);
void rcn_sao_first_pix_rows(OVCTUDec *const ctudec, const struct RectEntryInfo *const einfo, uint16_t ctb_y);

#endif

