
/* Solve the 5-point relative pose problem */

#ifndef _MVGLIB_5POINT_RELATIVE_POSE_H_
#define _MVGLIB_5POINT_RELATIVE_POSE_H_


#include "easy3d/core/types.h"


int compute_pose_ransac(
	int n, easy3d::dvec2* r_pts, easy3d::dvec2* l_pts,
	double *K1, double *K2,
	double ransac_threshold, int ransac_rounds,
	double *R_out, double *t_out
	);



#endif /* _MVGLIB_5POINT_RELATIVE_POSE_H_ */
