
#ifndef _OPENGL_SHADER_FACTORY_H_
#define _OPENGL_SHADER_FACTORY_H_

#include "opengl_common.h"

#include <string>


// class Canvas;
class VertexShader;
class FragmentShader;
class GeometryShader;

class OPENGL_API ShaderFactory
{
public:
	static std::string title() { return "ShaderFactory"; }

	static bool read_shader_souce_file(const std::string& file_name, std::string& source) ;

	static VertexShader*	create_vertex_shader_from_file(const std::string& file_name) ;
	static VertexShader*	create_vertex_shader_from_string(const std::string& code) ;

	static FragmentShader*  create_fragment_shader_from_file(const std::string& file_name) ;
	static FragmentShader*  create_fragment_shader_from_string(const std::string& code) ;
	
	static GeometryShader*  create_geometry_shader_from_file(const std::string& file_name) ;
	static GeometryShader*  create_geometry_shader_from_string(const std::string& code) ;
} ;



#endif
