
/* Triangulate two image points */

#ifndef _MVGlIB_TRIANGULATE_H_
#define _MVGlIB_TRIANGULATE_H_

#include "easy3d/core/types.h"



/* Project a point onto an image */
easy3d::dvec2 project(double *R, double *t0, double *P);

/* Find the point with the smallest squared projection error */
easy3d::dvec3 triangulate(easy3d::dvec2 p, easy3d::dvec2 q,
	double *R0, double *t0,
	double *R1, double *t1, double *error);

/* Find the point with the smallest squared projection error */
easy3d::dvec3 triangulate_n(int num_points,
	easy3d::dvec2 *p, double *R, double *t, double *error_out);

/* Find the point with the smallest squared projection error */
easy3d::dvec3 triangulate_n_refine(easy3d::dvec3 pt, int num_points,
	easy3d::dvec2 *p, double *R, double *t, double *error_out);

/* Given an F matrix, two calibration matrices, and a point
* correspondence, find R and t */
void find_extrinsics(double *F, double *K1, double *K2,
	easy3d::dvec2 p1, easy3d::dvec2 p2, double *R, double *t);

/* Given an E matrix and a point correspondence, find R and t */
int find_extrinsics_essential(double *E, easy3d::dvec2 p1, easy3d::dvec2 p2,
	double *R, double *t);

int find_extrinsics_essential_multipt(double *E, int n,
	easy3d::dvec2 *p1, easy3d::dvec2 *p2,
	double *R, double *t);

/* Solve for a 3x4 projection matrix, given a set of 3D points and 2D
* projections */
int find_projection_3x4(int num_pts, easy3d::dvec3 *points, easy3d::dvec2 *projs, double *P);

/* Solve for a 3x4 projection matrix using RANSAC, given a set of 3D
* points and 2D projections */
int find_projection_3x4_ransac(int num_pts, easy3d::dvec3 *points, easy3d::dvec2 *projs,
	double *P,
	int ransac_rounds, double ransac_threshold);


#endif /* _MVGlIB_TRIANGULATE_H_ */
