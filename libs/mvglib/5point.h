
/* Solve the 5-point relative pose problem */

#ifndef _MVGLIB_5POINT_RELATIVE_POSE_H_
#define _MVGLIB_5POINT_RELATIVE_POSE_H_


#include "../math/math_types.h"


int compute_pose_ransac(
	int n, vec2d* r_pts, vec2d* l_pts,
	double *K1, double *K2,
	double ransac_threshold, int ransac_rounds,
	double *R_out, double *t_out
	);



#endif /* _MVGLIB_5POINT_RELATIVE_POSE_H_ */
