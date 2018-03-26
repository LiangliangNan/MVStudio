
#ifndef _OPENGL_GLINFORMATION_H_
#define _OPENGL_GLINFORMATION_H_

#include "opengl_common.h"
#include "glew.h"

#include <string>



#ifndef NDEBUG
#define ogf_check_gl {\
	GLInfo::check_gl(__FILE__, __LINE__) ;\
}
#else
#define ogf_check_gl
#endif

class OPENGL_API GLInfo
{
public:
	/**
	* Prints the last GL error to the Logger.
	*/
	static void check_gl(const std::string& file, int line) ;

	static std::string gl_vendor() ;
	static std::string gl_renderer() ;
	static std::string gl_version() ;
	static std::string gl_extensions() ;

	static std::string glew_version() ;
	static std::string glsl_version() ;
} ;


#endif


