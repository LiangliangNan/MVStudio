
#include "shader_factory.h"
#include "shader.h"
#include "opengl_info.h"
#include "../basic/logger.h"
#include "../basic/basic_types.h"


bool ShaderFactory::read_shader_souce_file(const std::string& file_name, std::string& source) {
	std::ifstream in(file_name.c_str()) ;
	if(!in) {
		Logger::err(title()) << "could not open file:"
			<< "\'" << file_name << "\'"
			<< std::endl ;
		return false ;
	}
	while(in) {
		char line[1024] ;
		in.getline(line,1024) ; 
		source += std::string(line) ;
		source += "\n" ;
	}

	// If there is a END in the middle of the program,
	// terminate the string (so that everything after
	// it will be ignored)
	char* p = (char*)strstr(source.c_str(), "\nEND\n") ;
	if(p != nil) {
		*(p+4) = '\0' ;
	}
	return true ;
}


VertexShader* ShaderFactory::create_vertex_shader_from_file(const std::string& file_name) 
{
	std::string source ;
	if(read_shader_souce_file(file_name, source)) {
		return create_vertex_shader_from_string(source) ;
	} else {
		return nil ;
	}
}

FragmentShader* ShaderFactory::create_fragment_shader_from_file(const std::string& file_name) 
{
	std::string source ;
	if(read_shader_souce_file(file_name, source)) {
		return create_fragment_shader_from_string(source) ;
	} else {
		return nil ;
	}
}

GeometryShader* ShaderFactory::create_geometry_shader_from_file(const std::string& file_name) 
{
	std::string source ;
	if(read_shader_souce_file(file_name, source)) {
		return create_geometry_shader_from_string(source) ;
	} else {
		return nil ;
	}
}

VertexShader* ShaderFactory::create_vertex_shader_from_string(const std::string& code) {
	// Had problems with ATI (the rasterpos does not seem to vary within a GL_POINT)
	std::string vendor = GLInfo::gl_vendor() ;
	if(vendor.length() < 6 || vendor.substr(0,6) != "NVIDIA") {
		return nil ;
	}

	if (!Shader::is_supported()) {
		Logger::err(title()) << "GLSL vertex shader unsupported" << std::endl ;
		return nil;	
	}

	VertexShader* result = new VertexShader() ;
	ogf_assert(result->is_supported());

	if(!result->compile(code.c_str())) {
		Logger::err(title()) << "error occurred while parsing vertex shader" << std::endl ;
		delete result ;
		return nil ;
	}

	return result ;
}

FragmentShader* ShaderFactory::create_fragment_shader_from_string(const std::string& code) {
	// Had problems with ATI (the rasterpos does not seem to vary within a GL_POINT)
	std::string vendor = GLInfo::gl_vendor() ;
	if(vendor.length() < 6 || vendor.substr(0,6) != "NVIDIA") {
		return nil ;
	}

	if (!FragmentShader::is_supported()) {
		Logger::err(title()) << "GLSL fragment shader unsupported" << std::endl ;
		return nil;	
	}

	FragmentShader* result = new FragmentShader() ;
	ogf_assert(result->is_supported());

	if(!result->compile(code.c_str())) {
		Logger::err(title()) << "error occurred while parsing fragment shader" << std::endl ;
		delete result ;
		return nil ;
	}

	return result ;
}

GeometryShader* ShaderFactory::create_geometry_shader_from_string(const std::string& code) {
	// Had problems with ATI (the rasterpos does not seem to vary within a GL_POINT)
	std::string vendor = GLInfo::gl_vendor() ;
	if(vendor.length() < 6 || vendor.substr(0,6) != "NVIDIA") {
		return nil ;
	}

	if (!GeometryShader::is_supported()) {
		Logger::err(title()) << "GLSL geometry shader unsupported" << std::endl ;
		return nil;	
	}

	GeometryShader* result = new GeometryShader() ;
	ogf_assert(result->is_supported());

	if(!result->compile(code.c_str())) {
		Logger::err(title()) << "error occurred while parsing geometry shader" << std::endl ;
		delete result ;
		return nil ;
	}

	return result ;
}

