
/* homography.h */
/* Computes a homography */

#ifndef _MVGLIB_HOMOGRAPHY_H_
#define _MVGLIB_HOMOGRAPHY_H_


#include "easy3d/core/types.h"



/* Computes the homography that, when applied to the points in l_pts,
* minimizes the least-squares error between the result and the
* corresponding points in r_pts.
*
* n -- number of points
* r_pts -- matches
* l_pts -- initial points
* Tout -- on return, contains the 3x3 transformation matrix */
void align_homography(
	int num_pts, easy3d::dvec3* r_pts, easy3d::dvec3* l_pts, double* Tout, int refine
	);

/* Use non-linear least squares to refine a homography */
void align_homography_non_linear(
	int num_pts, easy3d::dvec3* r_pts, easy3d::dvec3* l_pts, double* Tin, double* Tout
	);



#endif /* _SFM_HOMOGRAPHY_H_ */
