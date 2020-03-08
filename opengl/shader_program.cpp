#include "shader_program.h"
#include "../basic/file_utils.h"
#include "../basic/logger.h"


bool ShaderProgram::load_shader_file(GLuint& shader_id, GLuint shader_type, std::string const& filename) {
    std::string shader_code;
    FileUtils::read_file_to_string(filename, shader_code);

	if(shader_code.empty()) {
		Logger::err(title()) << "failed reading shader file" << std::endl;
		return false;
	}

	if(!load_shader_code(shader_id, shader_type, shader_code)) {
		Logger::err(title()) << "failed reading shader file" << std::endl;
		return false;
	}

	return true;
}

bool ShaderProgram::load_shader_code(GLuint& shader_id, GLuint shader_type, std::string const& code) {
    if(shader_id == 0) {
        shader_id = glCreateShader(shader_type);
        ogf_check_gl;

        glAttachShader(prog_id, shader_id);
        ogf_check_gl;
    }

	if (compile_shader(shader_id, code)) {
		need_to_link = true;
		return true;
	} else
		return false;
}

/* ---------------------------------------------------------------- */

void ShaderProgram::unload_shader(GLuint& shader_id) {
    if(shader_id != 0)
    {
        glDetachShader(prog_id, shader_id);
        ogf_check_gl;
        glDeleteShader(shader_id);
        ogf_check_gl;
        shader_id = 0;
    }
}

/* ---------------------------------------------------------------- */

bool ShaderProgram::compile_shader(GLuint shader_id, std::string const& code) {
    /* Pass code to OpenGL. */
    char const* data[1] = { code.c_str() };
    glShaderSource(shader_id, 1, data, NULL);
    ogf_check_gl;

    /* Compile shader. */
    glCompileShader(shader_id);
    ogf_check_gl;
    if(get_shader_property(shader_id, GL_COMPILE_STATUS) == GL_FALSE) {
        GLint log_size = get_shader_property(shader_id, GL_INFO_LOG_LENGTH);
		if(log_size == 0) 
			Logger::err(title()) << "shader compilation failed" << std::endl;

        std::string log;
        log.append(log_size + 1, '\0');
        glGetShaderInfoLog(shader_id, log_size + 1, NULL, &log[0]);
		Logger::err(title()) << log << std::endl;
		return false;
    } else
		return true;
}

/* ---------------------------------------------------------------- */

bool ShaderProgram::try_load_all(std::string const& basename) {
    std::string vert_filename = basename + ".vert";
    std::string geom_filename = basename + ".geom";
    std::string frag_filename = basename + ".frag";

	if(!FileUtils::is_file(vert_filename) || !FileUtils::is_file(frag_filename)) {
        std::cerr << "Skipping shaders from " << basename << ".*" << std::endl;
        return false;
    }

    std::cerr << "Loading shaders from " << basename << ".*" << std::endl;

    load_vert_file(vert_filename);

	if(FileUtils::is_file(geom_filename))
        load_geom_file(geom_filename);

    load_frag_file(frag_filename);

    return true;
}


ShaderProgram::ShaderProgram (void) {
	vert_id = 0;
	geom_id = 0;
	frag_id = 0;
	prog_id = glCreateProgram();
	ogf_check_gl;
	need_to_link = false;
}

ShaderProgram::~ShaderProgram (void) {
	glDeleteProgram(prog_id);
	ogf_check_gl;
	glDeleteShader(vert_id);
	ogf_check_gl;
	glDeleteShader(geom_id);
	ogf_check_gl;
	glDeleteShader(frag_id);
	ogf_check_gl;
}

ShaderProgram::Ptr ShaderProgram::create (void) {
	return Ptr(new ShaderProgram);
}

bool ShaderProgram::load_vert_file (std::string const& filename) {
	return load_shader_file(vert_id, GL_VERTEX_SHADER, filename);
}

bool ShaderProgram::load_geom_file (std::string const& filename) {
	return load_shader_file(geom_id, GL_GEOMETRY_SHADER, filename);
}

bool ShaderProgram::load_frag_file (std::string const& filename) {
	return load_shader_file(frag_id, GL_FRAGMENT_SHADER, filename);
}

bool ShaderProgram::load_vert_code (std::string const& code) {
	return load_shader_code(vert_id, GL_VERTEX_SHADER, code);
}

bool ShaderProgram::load_geom_code (std::string const& code) {
	return load_shader_code(geom_id, GL_GEOMETRY_SHADER, code);
}

bool ShaderProgram::load_frag_code (std::string const& code) {
	return load_shader_code(frag_id, GL_FRAGMENT_SHADER, code);
}

void ShaderProgram::unload_vert (void) {
	glDetachShader(prog_id, vert_id);
	ogf_check_gl;

	glDeleteShader(vert_id);
	ogf_check_gl;

	vert_id = 0;
}

void ShaderProgram::unload_geom (void) {
	glDetachShader(prog_id, geom_id);
	ogf_check_gl;

	glDeleteShader(geom_id);
	ogf_check_gl;

	geom_id = 0;
}

void ShaderProgram::unload_frag (void) {
	glDetachShader(prog_id, frag_id);
	ogf_check_gl;

	glDeleteShader(frag_id);
	ogf_check_gl;

	frag_id = 0;
}

GLint ShaderProgram::get_attrib_location (char const* name) {
	ensure_linked();
	return glGetAttribLocation(prog_id, name);
}

GLint ShaderProgram::get_uniform_location (char const* name) {
	ensure_linked();
	GLint loc = glGetUniformLocation(prog_id, name);
	return loc;
}

void ShaderProgram::send_uniform (char const* name, Vec3f const& v) {
	GLint loc = get_uniform_location(name);
	if (loc < 0)
		return;
	glUniform3fv(loc, 1, v.data());
}

void ShaderProgram::send_uniform (char const* name, Vec4f const& v) {
	GLint loc = get_uniform_location(name);
	if (loc < 0)
		return;
	glUniform4fv(loc, 1, v.data());
}

void ShaderProgram::send_uniform (const char* name, Matrix4f const& m) {
	GLint loc = get_uniform_location(name);
	if (loc < 0)
		return;
	glUniformMatrix4fv(loc, 1, true, m.data());
}

void ShaderProgram::send_uniform (const char* name, GLint val) {
	GLint loc = get_uniform_location(name);
	if (loc < 0)
		return;
	glUniform1i(loc, val);
}

void ShaderProgram::send_uniform (const char* name, GLfloat val) {
	GLint loc = get_uniform_location(name);
	if (loc < 0)
		return;
	glUniform1f(loc, val);
}

void ShaderProgram::bind (void) {
	ensure_linked();
	glUseProgram(prog_id);
	ogf_check_gl;
}

void ShaderProgram::unbind (void) const {
	glUseProgram(0);
	ogf_check_gl;
}

void ShaderProgram::ensure_linked (void) {
	if (need_to_link)
	{
		glLinkProgram(prog_id);
		ogf_check_gl;
		need_to_link = false;
	}
}

GLint ShaderProgram::get_program_property (int pname) {
	GLint ret;
	glGetProgramiv(prog_id, pname, &ret);
	ogf_check_gl;
	return ret;
}

GLint ShaderProgram::get_shader_property (GLuint shader_id, int pname) {
	GLint ret;
	glGetShaderiv(shader_id, pname, &ret);
	ogf_check_gl;
	return ret;
}
