#include <stdint.h>

#include "ovvcerror.h"

const char * ovvc_error_stringify(uint32_t error_code){
  switch (error_code) {
    case OVVC_EAGAIN:
      return "No data available. Try again later.";
      break;
    case OVVC_EINDATA:
      return "Invalide data.";
      break;
    case OVVC_ENOMEM:
      return "Missing memory.";
      break;
    case OVVC_EUNSUPPORTED:
      return "Unsupported tool.";
      break;
    default:
      return "Error code not recognized.";
  }
}
