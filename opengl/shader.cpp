
#include "shader.h"
#include "../basic/logger.h"
#include "../basic/basic_types.h"


GLhandleARB Shader::program_id_ = 0 ;
bool Shader::initialized_ = false ;
bool Shader::error_       = false ;
bool Shader::activated_   = false ;


Shader::Shader() 
: id_(0)
{
	init() ;
}

Shader::~Shader() {
}

bool Shader::init() {
	if(initialized_) {
		return !error_ ;
	}
	initialized_ = true ;
	bool ok = glewIsSupported("GL_ARB_shader_objects") != 0 ;
    ok = ok && (glewIsSupported("GL_ARB_vertex_shader") != 0) ;
    ok = ok && (glewIsSupported("GL_ARB_fragment_shader") != 0) ;
	// ok = ok && glewIsSupported("GL_ARB_shading_language_100") ;
	// Note: It seems that GL_ARB_shading_language_100 does not have any associated
	// function prototype, so we do not need initializing it (and trying to
	// initialize it returns an error in glew)
	error_ = !ok ;

	if(ok) {
		program_id_ = glCreateProgramObjectARB() ;
	} else {
		Logger::err("Shader") << "could not initialize OpenGL Shading Language" << std::endl ;
	}

	return ok ;
}


bool Shader::is_supported() {
	return init() ;
}


void Shader::update() {
	if(VertexShader::current() == nil || FragmentShader::current() == nil) {
		glUseProgramObjectARB(0) ;
		activated_ = false;
	} else {
		glLinkProgramARB(program_id_) ;

		GLint status ;
		glGetObjectParameterivARB(program_id_, GL_OBJECT_LINK_STATUS_ARB, &status) ;

		bool ok = (status != GL_FALSE) ;
		if(!ok) {
			// Note: we only display errors related with incompatibility between vertex and fragment shader
			// (when we use the RenderingPipeline API, we declare the shaders one by one, and the system 
			// may bark in the transient state if non-uniform variables are used)
			if(VertexShader::current() != nil && FragmentShader::current()  != nil) {
				static char log_buffer[4096] ;
				GLsizei length ;
				glGetInfoLogARB(program_id_, 4096, &length, log_buffer) ;
				Logger::err("Shader") << log_buffer << std::endl ;
			}
		}  else {
			glUseProgramObjectARB(program_id_) ;
			activated_ = true;
		}
	}
}

bool Shader::compile_shader(GLhandleARB shader, const char* code) {
	GLint length = strlen(code) ;
	glShaderSourceARB(shader, 1, &code, &length) ;
	glCompileShaderARB(shader) ;

	GLint status ;
	glGetObjectParameterivARB(shader, GL_OBJECT_COMPILE_STATUS_ARB, &status) ;

	bool ok = (status != GL_FALSE) ;
	if(!ok) {
		static char log_buffer[4096] ;
		GLsizei length ;
		glGetInfoLogARB(shader, 4096, &length, log_buffer) ;
		Logger::err("Shader") << log_buffer << std::endl ;
	}
	return ok ;
}

void Shader::add_uniform(const std::string name) {
	activate();
	uniform_ids_[name] = glGetUniformLocationARB(program_id_, name.c_str()) ;
	deactivate();
}


void Shader::set_program_parameter(const std::string& name, double x1) {
	ogf_assert(activated_);
	glUniform1fARB(uniform_ids_[name], float(x1)) ;
}

void Shader::set_program_parameter(const std::string& name, double x1, double x2) {
	ogf_assert(activated_);
	glUniform2fARB(uniform_ids_[name], float(x1), float(x2)) ;
}

void Shader::set_program_parameter(const std::string& name, double x1, double x2, double x3) {
	ogf_assert(activated_);
	glUniform3fARB(uniform_ids_[name], float(x1), float(x2), float(x3)) ;
}


void Shader::set_program_parameter(const std::string& name, double x1, double x2, double x3, double x4) {
	ogf_assert(activated_);
	glUniform4fARB(uniform_ids_[name], float(x1), float(x2), float(x3), float(x4)) ;
}

//________________________________________________________________________

VertexShader* VertexShader::current_ = nil ;

VertexShader::VertexShader() : Shader() {
	if(initialized_) {
		set_id(glCreateShaderObjectARB(GL_VERTEX_SHADER_ARB)) ;
	} else {
		set_id(0) ;
	}
}

VertexShader::~VertexShader() {
	if(id() != 0) {
		glDeleteObjectARB(id()) ;
	}
	set_id(0) ;
}

void VertexShader::activate() {
	glEnable(GL_VERTEX_PROGRAM_POINT_SIZE_ARB) ;
	glEnable(GL_POINT_SPRITE_ARB);

	glAttachObjectARB(program_id_, id()) ;

	current_ = this ;
	update() ;
}

void VertexShader::deactivate() {
	glDisable(GL_VERTEX_PROGRAM_POINT_SIZE_ARB) ;  
	glDisable(GL_POINT_SPRITE_ARB);

	glDetachObjectARB(program_id_, id()) ;

	current_ = nil ;
	update() ;
}

bool VertexShader::compile(const char* code) {
	return compile_shader(id(), code) ;
}


void VertexShader::set_vertex_parameter(int i, double x1, double x2, double x3, double x4) {
	glVertexAttrib4dARB(i, x1, x2, x3, x4) ;
}


//_________________________________________________________

FragmentShader* FragmentShader::current_ = nil ;

FragmentShader::FragmentShader() : Shader() {
	if(initialized_) {
		set_id(glCreateShaderObjectARB(GL_FRAGMENT_SHADER_ARB)) ;
	} else {
		set_id(0) ;
	}
}

FragmentShader::~FragmentShader() {
	if(id() != 0) {
		glDeleteObjectARB(id()) ;
	}
	set_id(0) ;
}

void FragmentShader::activate() {
	glAttachObjectARB(program_id_, id()) ;

	current_ = this ;
	update() ;
}

void FragmentShader::deactivate() {
	glDetachObjectARB(program_id_, id()) ;

	current_ = nil ;
	update() ;
}

bool FragmentShader::compile(const char* code) {
	return compile_shader(id(), code) ;
}


//_________________________________________________________

GeometryShader* GeometryShader::current_ = nil ;

GeometryShader::GeometryShader() : Shader() {
	if(initialized_) {
		set_id(glCreateShaderObjectARB(GL_GEOMETRY_SHADER_ARB)) ;
	} else {
		set_id(0) ;
	}
}

GeometryShader::~GeometryShader() {
	if(id() != 0) {
		glDeleteObjectARB(id()) ;
	}
	set_id(0) ;
}

void GeometryShader::activate() {
	glAttachObjectARB(program_id_, id()) ;

	current_ = this ;
	update() ;
}

void GeometryShader::deactivate() {
	glDetachObjectARB(program_id_, id()) ;

	current_ = nil ;
	update() ;
}

bool GeometryShader::compile(const char* code) {
	return compile_shader(id(), code) ;
}
