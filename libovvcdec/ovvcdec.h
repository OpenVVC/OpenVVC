#ifndef OPENVVC_H
#define OPENVVC_H

typedef struct OVVCDec OVVCDec;

int ovdec_init(OVVCDec **ovvcdec);

void ovdec_close(OVVCDec **ovvcdec);

#endif
