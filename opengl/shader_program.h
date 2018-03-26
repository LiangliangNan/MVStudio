
#ifndef _OPENGL_SHADER_PROGRAM_H_
#define _OPENGL_SHADER_PROGRAM_H_

#include "opengl_common.h"

#include <string>
#include "../basic/smart_pointer.h"
#include "../math/math_types.h"
#include "../basic/counted.h"
#include "opengl_info.h"

#define OGL_ATTRIB_POSITION		"pos"
#define OGL_ATTRIB_NORMAL		"normal"
#define OGL_ATTRIB_COLOR		"color"
#define OGL_ATTRIB_TEXCOORD		"texuv"


typedef		GeometricTypes::Vector3d_float32	Vec3f;
typedef		GeometricTypes::Vector4d_float32	Vec4f;
typedef		GeometricTypes::Matrix4_float32		Matrix4f;

/**
 * Abstraction for OpenGL Shader Programs.
 */
class OPENGL_API ShaderProgram : public Counted
{
public:
    typedef SmartPointer<ShaderProgram>			Ptr;
    typedef SmartPointer<const ShaderProgram>	ConstPtr;

public:
    ~ShaderProgram (void);
    
	static std::string title() { return "ShaderProgram"; }

	static Ptr create(void);

    /**
     * Try loading all shaders by appending ".vert", ".geom"
     * and ".frag" to basename.
     */
    bool try_load_all(std::string const& basename);

    /** Loads a vertex shader from file. */
    bool load_vert_file(std::string const& filename);
    /** Loads optional geometry shader from file. */
    bool load_geom_file(std::string const& filename);
    /** Load fragment shader from file. */
    bool load_frag_file(std::string const& filename);

    /** Loads a vertex shader from code in memory. */
    bool load_vert_code(std::string const& code);
    /** Loads optional geometry shader from code in memory. */
    bool load_geom_code(std::string const& code);
    /** Load fragment shader from code in memory. */
    bool load_frag_code(std::string const& code);

    /** Unloads vertex shader. */
    void unload_vert(void);
    /** Unloads geometry shader. */
    void unload_geom(void);
    /** Unloads fragment shader. */
    void unload_frag(void);

    /**
     * Returns attribute location for the program. If the program
     * has not yet been linked, it is linked first. If there is no
     * attribute by that name, -1 is returned.
     */
    GLint get_attrib_location(char const* name);

    /**
     * Returns the uniform location of the program. If the program
     * has not yet been linked, it is linked first. If there is no
     * uniform variable by that name, -1 is returned.
     */
    GLint get_uniform_location(char const* name);

    /** Sends 3-vector 'v' to uniform location 'name'. */
    void send_uniform(char const* name, Vec3f const& v);
    /** Sends 4-vector 'v' to uniform location 'name'. */
    void send_uniform(char const* name, Vec4f const& v);
    /** Sends 4x4-matrx 'm' to uniform location 'name'. */
    void send_uniform(char const* name, Matrix4f const& m);
    /** Sends integer 'val' to uniform location 'name'. */
    void send_uniform(char const* name, GLint val);
    /** Sends float 'val' to uniform location 'name'. */
    void send_uniform(char const* name, GLfloat val);

    /** Selects the shader program for rendering. */
    void bind(void);

    /** Deselects the currend shader program. */
    void unbind(void) const;

private:
    ShaderProgram (void);

    bool load_shader_file(GLuint& shader_id, GLuint shader_type, std::string const& filename);
    bool load_shader_code(GLuint& shader_id, GLuint shader_type, std::string const& code);
    void unload_shader(GLuint& shader_id);
    bool compile_shader(GLuint shader_id, std::string const& code);

    GLint get_program_property(int pname);
    GLint get_shader_property(GLuint shader_id, int pname);

    void  ensure_linked(void);

private:
    GLuint prog_id;
    GLuint vert_id;
    GLuint geom_id;
    GLuint frag_id;

    bool need_to_link;
};


#endif /* _OPENGL_SHADER_PROGRAM_H_ */
