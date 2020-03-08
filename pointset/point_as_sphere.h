
#ifndef _POINT_AS_SPHERE_H_
#define _POINT_AS_SPHERE_H_


#include "../opengl/shader.h"


class PointAsSphere
{
public:
	PointAsSphere();
	~PointAsSphere(void);

	bool points_as_spheres() const ;
	void set_points_as_spheres(bool x) ;

protected:
	void activate_points_shaders();
	void deactivate_points_shaders();

protected:
	bool				points_as_spheres_ ;
	static bool			has_points_shaders_ ;

	VertexShader_var	points_vertex_shader_;
	FragmentShader_var	points_fragment_shader_;
};



#endif