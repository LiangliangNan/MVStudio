
#ifndef _OPENGL_COMMON_COMMON_H_
#define _OPENGL_COMMON_COMMON_H_


#include "../basic/basic_common.h"


# ifdef OPENGL_EXPORTS
#   define OPENGL_API  EXPORT_LIBRARY
# else
#   define OPENGL_API  IMPORT_LIBRARY
# endif


#endif
