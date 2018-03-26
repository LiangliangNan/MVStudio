
#ifndef _OPENGL_SHADER_H_
#define _OPENGL_SHADER_H_

#include "opengl_common.h"
#include "../basic/counted.h"
#include "../basic/smart_pointer.h"

#include "glew.h"

#include <string>
#include <map>



class OPENGL_API Shader : public Counted
{
public:
	static bool is_supported() ;

	Shader() ;
	~Shader() ;

	void add_uniform(const std::string name) ;

	virtual void activate() = 0 ;
	virtual void deactivate() = 0 ;
	virtual bool compile(const char* source) = 0 ;

	virtual void set_program_parameter(const std::string& name, double x) ;
	virtual void set_program_parameter(const std::string& name, double x, double y) ;
	virtual void set_program_parameter(const std::string& name, double x, double y, double z) ;
	virtual void set_program_parameter(const std::string& name, double x, double y, double z, double w) ;

protected:
	static bool init() ;

	bool compile_shader(GLhandleARB shader, const char* source) ;
	void set_id(unsigned int x) { id_ = x ; }
	unsigned int id() const  { return id_ ; }

	void update() ;

protected:
	std::map<std::string, int>  uniform_ids_;

	static GLhandleARB program_id_ ;
	static bool initialized_ ;
	static bool error_ ;
	static bool	activated_;

private:
	unsigned int id_ ;
} ;


//_________________________________________________________

class OPENGL_API VertexShader : public Shader
{
public:
	VertexShader() ;
	virtual ~VertexShader() ;

	virtual void activate() ;
	virtual void deactivate() ;

	virtual bool compile(const char* code) ;

	virtual void set_vertex_parameter(int i, double x1, double x2, double x3, double x4) ;

protected:
	static VertexShader* current() { return current_ ; }
	friend class Shader ;

private:
	static VertexShader* current_ ;
} ;

//_________________________________________________________

class OPENGL_API FragmentShader : public Shader
{
public:
	FragmentShader() ;
	virtual ~FragmentShader() ;

	virtual void activate() ;
	virtual void deactivate() ;

	virtual bool compile(const char* code) ;

protected:
	static FragmentShader* current() { return current_ ; }
	friend class Shader ;

private:
	static FragmentShader* current_ ;
} ;

//_________________________________________________________

class OPENGL_API GeometryShader : public Shader
{
public:
	GeometryShader() ;
	virtual ~GeometryShader() ;

	virtual void activate() ;
	virtual void deactivate() ;

	virtual bool compile(const char* code) ;

protected:
	static GeometryShader* current() { return current_ ; }
	friend class Shader ;

private:
	static GeometryShader* current_ ;
} ;


//_________________________________________________________


typedef SmartPointer<VertexShader>		VertexShader_var ;
typedef SmartPointer<FragmentShader>	FragmentShader_var ;
typedef SmartPointer<GeometryShader>	GeometryShader_var ;



#endif
