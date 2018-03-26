
/* Compute the closed-form least-squares solution to a rigid body alignment */

#ifndef _MVGLIB_HORN_H_
#define _MVGLIB_HORN_H_

#include "../math/math_types.h"


/* Computes the closed-form least-squares solution to a rigid
 * body alignment.
 *
 * n: the number of points
 * right_pts: One set of n points (the set of matches)
 * left_pts:  The other set of n points (the set of input points)
 * R: on return, contains the optimal (2x2) rotation 
 * T: on return, contains the optimal (3x3) translation
 * weight: array of weights assigned to each pair of points.  Leave
 *         parameter NULL to weight all points equally */
double align_horn(int n, vec3d *right_pts, vec3d *left_pts, 
		  double *R, double *T, double *Tout, 
		  double *scale, double *weight);

/* Align two sets of points with a 3D rotation */
double align_3D_rotation(int n, vec3d *r_pts, vec3d *l_pts, double *R);

/* Computes the closed-form least-squares solution to a rigid
 * body alignment.
 *
 * n: the number of points
 * right_pts: Target set of n points 
 * left_pts:  Source set of n points */
double align_horn_3D(int n, vec3d *right_pts, vec3d *left_pts, int scale_xform, 
		     double *Tout);    

double align_horn_3D_2(int n, vec3d *right_pts, vec3d *left_pts, int scale_xform,
                       double *Tout);

int align_2D_ransac(int n, vec3d *r_pts, vec3d *l_pts, 
                    int num_ransac_rounds, double ransac_thresh,
                    double *Tret);

/* Align two sets of points with a 2D similarity transform */
int align_horn_ransac(int n, vec3d *r_pts, vec3d *l_pts, 
                      int num_ransac_rounds, double ransac_thresh,
                      double *Tret);

/* Align two sets of points with a 3D similarity transform */
int align_horn_3D_ransac(int n, vec3d *r_pts, vec3d *l_pts, 
                         int num_ransac_rounds, double ransac_thresh,
                         double *Tret);
    
/* Align two sets of points with a 3D rotation */
int align_3D_rotation_ransac(int n, vec3d *r_pts, vec3d *l_pts, 
			     int num_ransac_rounds, double ransac_thresh,
			     double *R);
    

#endif /* __horn_h__ */
