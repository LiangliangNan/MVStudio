
#ifndef _MATH_SVD_H_
#define _MATH_SVD_H_


#include "math_types.h"


int svd(int m, int n, int withu, int withv, double eps, double tol,
	double *a, double *q, double *u, double *v, double *vt);

/* Compute the mean of a set of vectors */
void vec_svd(int n, const vec3d *v, double *U, double *S, double *VT);

#endif /* _MATH_SVD_H_ */
