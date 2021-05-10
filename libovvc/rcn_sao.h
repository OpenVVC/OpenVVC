#ifndef RCN_SAO_H
#define RCN_SAO_H

#include "ctudec.h"

void rcn_sao_filter_line(OVCTUDec *const ctudec, int nb_ctu_w, uint16_t ctb_y_pic) ;
void rcn_sao_first_pix_rows(OVCTUDec *const ctudec, int nb_ctu_w, uint16_t ctb_y_pic) ; 

#endif

