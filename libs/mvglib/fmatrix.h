
/* Routines for estimating the fundamental matrix of a pair of images */

#ifndef _MVGLIB_FMATRIX_H_
#define _MVGLIB_FMATRIX_H_

#include "../math/math_types.h"




/* Compute the epipoles of an F-matrix */
void fmatrix_compute_epipoles(double *F, double *e1, double *e2);

/* Compute the distance from l to the epipolar line of r under F */
double fmatrix_compute_residual(double *F, vec3d r, vec3d l);

/* Use RANSAC to estimate an F-matrix */
int estimate_fmatrix_ransac_matches(int num_pts, vec3d *a_pts, vec3d *b_pts, 
                                    int num_trials, double threshold, 
                                    double success_ratio, 
                                    int essential, double *F);

/* Use linear least-squares to estimate the fundamantal matrix.  The
 * F-matrix is returned in Fout, and the two epipoles in e1 and e2 */
int estimate_fmatrix_linear(int num_pts, vec3d *r_pts, vec3d *l_pts, 
                            int essential, 
                            double *Fout, double *e1, double *e2);


/* Estimate the essential matrix from an F-matrix, assuming 
 * same focal lengths */
void estimate_essential_same_focal_lengths(double *F, double *alpha, double *E);

/* Estimate the essential matrix from an F-matrix, assuming 
 * different focal lengths */
void estimate_essential_different_focal_lengths(double *F, 
						double *calib1, double *calib2,
						double *E);

/* Find the closest rank 2 matrix to the given 3x3 matrix */
int closest_rank2_matrix(double *Fin, double *Fout, double *U, double *VT);

/* Refine an F-matrix estimate using LM */
void refine_fmatrix_nonlinear_matches(int num_pts, vec3d *r_pts, vec3d *l_pts, 
				      double *F0, double *Fout);

/* Compute an F-matrix from two sets of camera parameters */
void fmatrix_from_parameters(double *i0, double *R0, double *t0, double *i1, double *R1, double *t1, double *F);



#endif /* _MVGLIB_FMATRIX_H_ */
