
#ifndef _IMAGE_COMMON_H_
#define _IMAGE_COMMON_H_


#include "../basic/basic_common.h"


# ifdef IMAGE_EXPORTS
#   define IMAGE_API  EXPORT_LIBRARY
# else
#   define IMAGE_API  IMPORT_LIBRARY
# endif


#endif