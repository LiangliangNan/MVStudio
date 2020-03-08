
#ifndef _SFM_EXPORT_H_
#define _SFM_EXPORT_H_


#include "../basic/basic_common.h"


# ifdef SFM_EXPORTS
#   define SFM_API  EXPORT_LIBRARY
# else
#   define SFM_API  IMPORT_LIBRARY
# endif


#endif