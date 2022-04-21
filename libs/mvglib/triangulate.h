
/* Triangulate two image points */

#ifndef _MVGlIB_TRIANGULATE_H_
#define _MVGlIB_TRIANGULATE_H_

#include "../math/math_types.h"



/* Project a point onto an image */
vec2d project(double *R, double *t0, double *P);

/* Find the point with the smallest squared projection error */
vec3d triangulate(vec2d p, vec2d q,
	double *R0, double *t0,
	double *R1, double *t1, double *error);

/* Find the point with the smallest squared projection error */
vec3d triangulate_n(int num_points,
	vec2d *p, double *R, double *t, double *error_out);

/* Find the point with the smallest squared projection error */
vec3d triangulate_n_refine(vec3d pt, int num_points,
	vec2d *p, double *R, double *t, double *error_out);

/* Given an F matrix, two calibration matrices, and a point
* correspondence, find R and t */
void find_extrinsics(double *F, double *K1, double *K2,
	vec2d p1, vec2d p2, double *R, double *t);

/* Given an E matrix and a point correspondence, find R and t */
int find_extrinsics_essential(double *E, vec2d p1, vec2d p2,
	double *R, double *t);

int find_extrinsics_essential_multipt(double *E, int n,
	vec2d *p1, vec2d *p2,
	double *R, double *t);

/* Solve for a 3x4 projection matrix, given a set of 3D points and 2D
* projections */
int find_projection_3x4(int num_pts, vec3d *points, vec2d *projs, double *P);

/* Solve for a 3x4 projection matrix using RANSAC, given a set of 3D
* points and 2D projections */
int find_projection_3x4_ransac(int num_pts, vec3d *points, vec2d *projs,
	double *P,
	int ransac_rounds, double ransac_threshold);


#endif /* _MVGlIB_TRIANGULATE_H_ */
