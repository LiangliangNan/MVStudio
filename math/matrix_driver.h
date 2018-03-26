
#ifndef _SFM_OPTIMIZATION_H_
#define _SFM_OPTIMIZATION_H_

#include "math_common.h"


#ifdef __cplusplus
extern "C" {
#endif

	/* Fill a given matrix with an n x n identity matrix */
	void MATH_API matrix_ident(int n, double *A);

	/* Fill a given matrix with an m x n matrix of zeroes */
	void MATH_API matrix_zeroes(int m, int n, double *A);

	/* Compute the cross product of two 3 x 1 vectors */
	void MATH_API matrix_cross(const double *u, const double *v, double *w);

	/* Create the 3x3 cross product matrix from a 3-vector */
	void MATH_API matrix_cross_matrix(double *v, double *v_cross);

	/* Get the norm of the matrix */
	double MATH_API matrix_norm(int m, int n, double *A);

	/* Get the [squared] norm of the matrix */
	double MATH_API matrix_norm2(int m, int n, double *A);

	/* Transpose the m x n matrix A and put the result in the n x m matrix AT */
	void MATH_API matrix_transpose(int m, int n, double *A, double *AT);

	/* Compute the matrix product R = AB */
	void MATH_API matrix_product(int Am, int An, int Bm, int Bn,
		const double *A, const double *B, double *R);

	/* Compute the matrix product R = A^T B */
	void MATH_API matrix_transpose_product(int Am, int An, int Bm, int Bn, double *A, double *B, double *R);
	/* Compute the matrix product R = A B^T */
	void MATH_API matrix_transpose_product2(int Am, int An, int Bm, int Bn, double *A, double *B, double *R);

	void MATH_API matrix_product33(double *A, double *B, double *R);
	void MATH_API matrix_product121(double *A, double *b, double *r);
	void MATH_API matrix_product131(double *A, double *b, double *r);
	void MATH_API matrix_product331(double *A, double *b, double *r);
	void MATH_API matrix_product341(double *A, double *b, double *r);
	void MATH_API matrix_product44(double *A, double *B, double *R);
	void MATH_API matrix_product441(double *A, double *b, double *r);

	/* Compute the power of a matrix */
	void MATH_API matrix_power(int n, double *A, int pow, double *R);

	/* Compute the matrix sum R = A + B */
	void MATH_API matrix_sum(int Am, int An, int Bm, int Bn,
		double *A, double *B, double *R);

	/* Compute the matrix difference R = A - B */
	void MATH_API matrix_diff(int Am, int An, int Bm, int Bn, double *A, double *B, double *R);

	/* Compute the determinant of a 3x3 matrix */
	double MATH_API matrix_determinant3(double *A);

	/* Scale a matrix by a scalar */
	void MATH_API matrix_scale(int m, int n, double *A, double s, double *R);

	/* Print the given m x n matrix */
	void MATH_API matrix_print(int m, int n, double *A);

	/* Read a matrix from a file */
	void MATH_API matrix_read_file(int m, int n, double *matrix, char *fname);

	/* Write a matrix to a file */
	void MATH_API matrix_write_file(int m, int n, double *matrix, char *fname);

	/* Compute (transpose of) LU decomposition of A */
	void MATH_API matrix_lu(int n, double *A, double *LU, int *ipiv);
	void MATH_API matrix_lu_no_transpose(int n, double *A, double *LU, int *ipiv);

	/* Solve a system of equations using a precomputed LU decomposition */
	void MATH_API matrix_solve_lu(int n, double *LU, int *ipiv, double *b, double *x);

	/* Invert the n-by-n matrix A, storing the result in Ainv */
	void MATH_API matrix_invert(int n, double *A, double *Ainv);
	void MATH_API matrix_invert_inplace(int n, double *A);

	/* Convert a rotation matrix to axis and angle representation */
	void MATH_API matrix_to_axis_angle(double *R, double *axis, double *angle);
	void MATH_API axis_angle_to_matrix(double *axis, double angle, double *R);
	void MATH_API axis_angle_to_matrix4(double *axis, double angle, double *R);

	/* Convert a matrix to a normalize quaternion */
	void MATH_API matrix_to_quaternion(double *R, double *q);
	/* Convert a normalized quaternion to a matrix */
	void MATH_API quaternion_to_matrix(double *q, double *R);

	/* Decompose a square matrix into an orthogonal matrix and a symmetric
	* positive semidefinite matrix */
	void MATH_API matrix_polar_decomposition(int n, double *A, double *Q, double *S);

	/* Find the unit vector that minimizes ||Ax|| */
	void MATH_API matrix_minimum_unit_norm_solution(int m, int n, double *A, double *x);

	/* Driver for the minpack function lmdif, which uses
	* Levenberg-Marquardt for non-linear least squares minimization */
	void MATH_API lmdif_driver(void *fcn, int m, int n, double *xvec, double tol);
	void MATH_API lmdif_driver2(void *fcn, int m, int n, double *xvec, double tol);

	/* Driver for the lapack function dgelss, which finds x to minimize
	* norm(b - A * x) */
	void MATH_API dgelss_driver(double *A, double *b, double *x, int m, int n, int nrhs);
	void MATH_API dgelsy_driver(double *A, double *b, double *x, int m, int n, int nrhs);

	/* Version of above where matrix is already in column-major order */
	void MATH_API dgelsy_driver_transpose(double *A, double *b, double *x, int m, int n, int nrhs);

	/* Solve an n x n system */
	void MATH_API dgesv_driver(int n, double *A, double *b, double *x);

	/* n: the order of matrix A
	* A: matrix for which the eigenvectors/values are to be computed
	* evec: output array containing the eigenvectors
	* eval: output array containing the eigenvalues
	*
	* Note: Assumes the results are real! */
	int MATH_API dgeev_driver(int n, double *A, double *evec, double *eval);

	/* Compute singular value decomposition of an m x n matrix A */
	int MATH_API dgesvd_driver(int m, int n, double *A, double *U, double *S, double *VT);
	/* Compute singular value decomposition of an m x n matrix A
	* (only compute S and VT) */
	int MATH_API dgesvd_driver_vt(int m, int n, double *A, double *S, double *VT);

	/* Compute Cholesky decomposition of an nxn matrix */
	void MATH_API dpotrf_driver(int n, double *A, double *U);

	/* Compute a QR factorization of an m by n matrix A */
	void MATH_API dgeqrf_driver(int m, int n, double *A, double *Q, double *R);

	/* Compute an RQ factorization of an m by n matrix A */
	void MATH_API dgerqf_driver(int m, int n, double *A, double *R, double *Q);

#ifdef __cplusplus
}

#endif

#endif /* __matrix_h__ */
