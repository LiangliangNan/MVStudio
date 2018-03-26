#include "point_as_sphere.h"
#include "../opengl/shader_factory.h"
#include "../basic/logger.h"
#include "../basic/file_utils.h"



bool PointAsSphere::has_points_shaders_ = true;

PointAsSphere::PointAsSphere()
: points_as_spheres_(true)
{
}

PointAsSphere::~PointAsSphere(void)
{
}


bool PointAsSphere::points_as_spheres() const { 
	return points_as_spheres_ ; 
}

void PointAsSphere::set_points_as_spheres(bool x) {
	points_as_spheres_ = x;
}


void PointAsSphere::activate_points_shaders() {
	if (has_points_shaders_) {
		if(points_vertex_shader_.is_nil()) {
			std::string dir = FileUtils::MVStudio_resource_directory() + "/gpu/";

			std::string vertex_source = dir + "sphere.vs";
			points_vertex_shader_ = ShaderFactory::create_vertex_shader_from_file(vertex_source);
			if (points_vertex_shader_.is_nil()) {
				Logger::err("PointAsSphere") << "failed creating vertex program" << std::endl; 
				has_points_shaders_ = false;
				return;
			} else {
				//points_vertex_shader_->add_uniform(xxx);
			}
		}

		if (points_fragment_shader_.is_nil()) {
			std::string dir = FileUtils::MVStudio_resource_directory() + "/gpu/";

			std::string pixel_source  = dir + "sphere.fs";
			points_fragment_shader_ = ShaderFactory::create_fragment_shader_from_file(pixel_source);
			if (points_fragment_shader_.is_nil()) {
				Logger::err("PointAsSphere") << "failed creating pixel program" << std::endl; 
				has_points_shaders_ = false;
				return;
			} else {
				//points_fragment_shader_->add_uniform("lightDir");
			}
		}
	}

	if(!points_vertex_shader_.is_nil()) {
		points_vertex_shader_->activate() ;
		//points_vertex_shader_->set_program_parameter(xxx);
		//points_vertex_shader_->set_vertex_parameter(xxx);
	}
	if (!points_fragment_shader_.is_nil()) {
		points_fragment_shader_->activate() ;
		// 		float pos[4];
		// 		glGetLightfv(GL_LIGHT0, GL_POSITION, pos);
		// 		points_fragment_shader_->set_program_parameter("lightDir", pos[0], pos[1], pos[2]);
	}
}

void PointAsSphere::deactivate_points_shaders() {
	if(!points_vertex_shader_.is_nil()) {
		points_vertex_shader_->deactivate() ;
	}

	if(!points_fragment_shader_.is_nil()) {
		points_fragment_shader_->deactivate() ;
	}
}
