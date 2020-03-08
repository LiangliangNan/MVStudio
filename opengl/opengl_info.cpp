#include "opengl_info.h"
#include <iostream>


static const std::string err_msg = "error(null_string)";


// Convert an OpenGL error code into a descriptive string.
inline const char* gl_error_string(GLenum code) {
    switch (code)
    {
        case GL_NO_ERROR:
            return "No error";
        case GL_INVALID_ENUM:
            return "Invalid enum (An unacceptable value is specified for an enumerated argument)";
        case GL_INVALID_VALUE:
            return "Invalid value (A numeric argument is out of range)";
        case GL_INVALID_OPERATION:
            return "Invalid operation (The specified operation is not allowed in the current state)";
        case GL_OUT_OF_MEMORY:
            return "Out of memory (There is not enough memory left to execute the command)";
#ifdef GL_STACK_OVERFLOW
        case GL_STACK_OVERFLOW:
                return "Stack overflow (An attempt has been made to perform an operation that would cause an internal stack to overflow)";
#endif
#ifdef GL_STACK_UNDERFLOW
        case GL_STACK_UNDERFLOW:
                return "Stack underflow (An attempt has been made to perform an operation that would cause an internal stack to underflow)";
#endif
        case GL_TABLE_TOO_LARGE:
            return "Table too large";
        case GL_INVALID_FRAMEBUFFER_OPERATION:
            return "Invalid framebuffer operation (The framebuffer object is not complete)";
    }
    return "Unknown error";
}

void GLInfo::check_gl(const std::string& file, int line) {
	GLenum error_code = glGetError() ;
	if(error_code != GL_NO_ERROR) {
		std::string str = (char*)(gl_error_string(error_code));
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

