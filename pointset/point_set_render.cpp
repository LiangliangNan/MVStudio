#include "point_set_render.h"
#include "../pointset/point_set.h"

#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif




PointSetRender::PointSetRender() {
	point_set_ = nil;
	points_as_spheres_ = false;

	point_size_ = 4 ;
	lighting_ = true;
}


PointSetRender::~PointSetRender(void)
{
}


void PointSetRender::draw() {
	if (!point_set_)
		return;

	int num = target()->num_points();
	if (num < 1)
		return;

	if (points_as_spheres_ && has_points_shaders_)
		activate_points_shaders();

	float* points = &(point_set_->points()[0].x);
	float* colors = 0;
	if (target()->has_colors())
		colors = &(point_set_->colors()[0].x);
	else {
		colors = new float[num * 3];
		memset(colors, 0, num * 3 * 4);
	}

	glPointSize(point_size_);

	bool has_normals = point_set_->has_normals();
	if (has_normals) {
		float* normals = &(point_set_->normals()[0].x);
		if (lighting_)
			glEnable(GL_LIGHTING);
		else
			glDisable(GL_LIGHTING);

		glEnableClientState(GL_VERTEX_ARRAY);
		glEnableClientState(GL_NORMAL_ARRAY);
		glEnableClientState(GL_COLOR_ARRAY);
		glVertexPointer(3, GL_FLOAT, 0, points);
		glNormalPointer(GL_FLOAT, 0, normals);
		glColorPointer(3, GL_FLOAT, 0, colors);
		glDrawArrays(GL_POINTS, 0, num);
		glDisableClientState(GL_COLOR_ARRAY);
		glDisableClientState(GL_NORMAL_ARRAY);
		glDisableClientState(GL_VERTEX_ARRAY);
	}
	else {
		glDisable(GL_LIGHTING); // always off for points without normals

		glEnableClientState(GL_VERTEX_ARRAY);
		glEnableClientState(GL_COLOR_ARRAY);
		glVertexPointer(3, GL_FLOAT, 0, points);
		glColorPointer(3, GL_FLOAT, 0, colors);
		glDrawArrays(GL_POINTS, 0, num);
		glDisableClientState(GL_COLOR_ARRAY);
		glDisableClientState(GL_VERTEX_ARRAY);
		glEnable(GL_LIGHTING);
	}

	if (points_as_spheres_ && has_points_shaders_)
		deactivate_points_shaders();
}

