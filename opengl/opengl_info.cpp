#include "opengl_info.h"
#include <iostream>


static const std::string err_msg = "error(null_string)";

void GLInfo::check_gl(const std::string& file, int line) {
	GLenum error_code = glGetError() ;
	if(error_code != GL_NO_ERROR) {
		std::string str = (char*)(gluErrorString(error_code));
		std::cerr << "GL error in file \'" << file << "\' @ line " << line << ": " << str << std::endl ;
	}
}

std::string GLInfo::gl_vendor() {
	const GLubyte* str = glGetString(GL_VENDOR) ;
	return std::string(reinterpret_cast<const char*>(str)) ;
}

std::string GLInfo::gl_renderer() {
	const GLubyte* str = glGetString(GL_RENDERER) ;
	if (str)
		return std::string(reinterpret_cast<const char*>(str)) ;
	else
		return err_msg;
}

std::string GLInfo::gl_version() {
	const GLubyte* str = glGetString(GL_VERSION) ;
	if (str)
		return std::string(reinterpret_cast<const char*>(str)) ;
	else
		return err_msg;
}

std::string GLInfo::gl_extensions() {
	const GLubyte* str = glGetString(GL_EXTENSIONS) ;
	if (str)
		return std::string(reinterpret_cast<const char*>(str)) ;
	else
		return err_msg;
}

std::string GLInfo::glew_version() {
	const GLubyte* str = glewGetString(GLEW_VERSION) ;
	if (str)
		return std::string(reinterpret_cast<const char*>(str)) ;
	else
		return err_msg;
}

std::string GLInfo::glsl_version() {
	const GLubyte* str = glGetString(GL_SHADING_LANGUAGE_VERSION) ;
	if (str)
		return std::string(reinterpret_cast<const char*>(str)) ;
	else
		return "not supported";
}

