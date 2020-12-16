#ifndef OVVC_ERROR_H
#define OVVC_ERROR_H

#include <stdint.h>

/*List of error*/
#define OVVC_ENOMEM           OVVC_ERROR_TAG('N','M','E','M')
#define OVVC_EINDATA          OVVC_ERROR_TAG('I','N','D','A')
#define OVVC_EUNSUPPORTED     OVVC_ERROR_TAG('U','N','S','P')
#define OVVC_EAGAIN           OVVC_ERROR_TAG('E','A','G','N')

#define OVVC_ERROR_TAG(a,b,c,d) -(((a)<<24)+((b)<<16)+((c)<<8)+(d))

#define OVVC_MAX_ERR_STRLEN 30

const char * ovvc_error_stringify(uint32_t error_code);

#endif/*OVVC_ERROR_H*/
