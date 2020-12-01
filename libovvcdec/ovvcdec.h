#ifndef OPENVVC_H
#define OPENVVC_H

typedef struct OVVCDec OVVCDec;

int ovdec_init(OVVCDec **ovvcdec);

int ovdec_close(OVVCDec *ovvcdec);

#endif
