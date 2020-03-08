#include "horn.h"
#include "qsort.h"
#include "../math/matrix_driver.h"
#include "../math/svd.h"

#include <assert.h>



#ifdef WIN32
#define isnan _isnan
//#define isinf _isinf
#endif

using namespace easy3d;

/* Computes the closed-form least-squares solution to a rigid
* body alignment.
*
* n: the number of points
* right_pts: Target set of n points
* left_pts:  Source set of n points */
double align_horn(int n, dvec3 *right_pts, dvec3 *left_pts,
	double *R, double *T,
	double *Tout, double *scale, double *weight) {
	int i;
	dvec3 right_centroid(0.0, 0.0, 0.0);
	dvec3 left_centroid(0.0, 0.0, 0.0);
	double M[2][2] = { { 0.0, 0.0 },
	{ 0.0, 0.0 } };
	double MT[2][2];
	double MTM[2][2];
	double eval[2], sqrteval[2];
	double evec[2][2];
	double S[2][2], Sinv[2][2], U[2][2];
	double Tcenter[3][3] = { { 1.0, 0.0, 0.0 },
	{ 0.0, 1.0, 0.0 },
	{ 0.0, 0.0, 1.0 } };

	double Ttmp[3][3];

	double sum_num, sum_den, RMS_sum;

#if 1
	double weight_sum = 0.0;

	if (weight == NULL) {
		weight_sum = n;

		for (i = 0; i < n; i++) {
			right_centroid =
				(right_centroid + right_pts[i]);
			left_centroid =
				(left_centroid + left_pts[i]);
		}

		right_centroid = (1.0 / weight_sum * right_centroid);
		left_centroid = (1.0 / weight_sum * left_centroid);
	}
	else {
		/* Compute the weighted centroid of both point sets */
		for (i = 0; i < n; i++) {
			right_centroid = right_centroid + (weight[i] * right_pts[i]);
			left_centroid = left_centroid + (weight[i] * left_pts[i]);
			weight_sum += weight[i];
		}

		right_centroid = 1.0 / weight_sum * right_centroid;
		left_centroid = 1.0 / weight_sum * left_centroid;
	}
#else
	/* Calculate the centroid of both sets of points */
	for (i = 0; i < n; i++) {
		right_centroid = v3_add(right_centroid, right_pts[i]);
		left_centroid = v3_add(left_centroid, left_pts[i]);
	}

	right_centroid = v3_scale(1.0 / n, right_centroid);
	left_centroid = v3_scale(1.0 / n, left_centroid);
#endif

	/* Compute the scale */
	sum_num = sum_den = 0.0;

	for (i = 0; i < n; i++) {
		dvec3 r = right_centroid - right_pts[i];
		dvec3 l = left_centroid - left_pts[i];

		sum_num = r.length2();
		sum_den = l.length2();
	}

	*scale = sqrt(sum_num / sum_den);

	/* Fill in the matrix M */
	for (i = 0; i < n; i++) {
		dvec3 r = right_centroid - right_pts[i];
		dvec3 l = left_centroid - left_pts[i];

		if (weight != NULL) {
			M[0][0] += r.x * l.x;
			M[0][1] += r.x * l.y;
			M[1][0] += r.y * l.x;
			M[1][1] += r.y * l.y;
		}
		else {
			M[0][0] += r.x * l.x;
			M[0][1] += r.x * l.y;
			M[1][0] += r.y * l.x;
			M[1][1] += r.y * l.y;
		}
	}

	/* Compute MTM */
	matrix_transpose(2, 2, (double *)M, (double *)MT);
	matrix_product(2, 2, 2, 2, (double *)MT, (double *)M, (double *)MTM);

	/* Calculate Sinv, the inverse of the square root of MTM */
	dgeev_driver(2, (double *)MTM, (double *)evec, eval);

	/* MTM = eval[0] * evec[0]T * evec[0] + eval[1] * evec[1]T * evec[1] */
	/* S = sqrt(eval[0]) * evec[0]T * evec[0] + sqrt(eval[1]) * evec[1]T * evec[1] */
	sqrteval[0] = sqrt(eval[0]);
	sqrteval[1] = sqrt(eval[1]);

	S[0][0] =
		(sqrteval[0]) * evec[0][0] * evec[0][0] +
		(sqrteval[1]) * evec[1][0] * evec[1][0];
	S[0][1] =
		(sqrteval[0]) * evec[0][0] * evec[0][1] +
		(sqrteval[1]) * evec[1][0] * evec[1][1];
	S[1][0] =
		(sqrteval[0]) * evec[0][1] * evec[0][0] +
		(sqrteval[1]) * evec[1][1] * evec[1][0];
	S[1][1] =
		(sqrteval[0]) * evec[0][1] * evec[0][1] +
		(sqrteval[1]) * evec[1][1] * evec[1][1];

	Sinv[0][0] =
		(1.0 / sqrteval[0]) * evec[0][0] * evec[0][0] +
		(1.0 / sqrteval[1]) * evec[1][0] * evec[1][0];
	Sinv[0][1] =
		(1.0 / sqrteval[0]) * evec[0][0] * evec[0][1] +
		(1.0 / sqrteval[1]) * evec[1][0] * evec[1][1];
	Sinv[1][0] =
		(1.0 / sqrteval[0]) * evec[0][1] * evec[0][0] +
		(1.0 / sqrteval[1]) * evec[1][1] * evec[1][0];
	Sinv[1][1] =
		(1.0 / sqrteval[0]) * evec[0][1] * evec[0][1] +
		(1.0 / sqrteval[1]) * evec[1][1] * evec[1][1];

	// matrix_product(2, 2, 2, 2, (double *)S, (double *)Sinv, (double *)U);

	/* U = M * Sinv */
	matrix_product(2, 2, 2, 2, (double *)M, (double *)Sinv, (double *)U);

	/* Fill in the rotation matrix */
	R[0] = U[0][0]; R[1] = U[0][1]; R[2] = 0.0;
	R[3] = U[1][0], R[4] = U[1][1]; R[5] = 0.0;
	R[6] = 0.0;     R[7] = 0.0;     R[8] = 1.0;

	// memcpy(R, U, sizeof(double) * 4);

	/* Fill in the translation matrix */
	T[0] = T[4] = T[8] = 1.0;
	T[1] = T[3] = T[6] = T[7] = 0.0;
	T[2] = right_centroid.x;
	T[5] = right_centroid.y;

	Tcenter[0][0] = *scale;
	Tcenter[1][1] = *scale;
	Tcenter[0][2] = -*scale * left_centroid.x;
	Tcenter[1][2] = -*scale * left_centroid.y;

	matrix_product(3, 3, 3, 3, T, R, (double *)Ttmp);


	matrix_product(3, 3, 3, 3, (double *)Ttmp, (double *)Tcenter, Tout);

	T[2] = (right_centroid - left_centroid).x;
	T[5] = (right_centroid - left_centroid).y;


	/* Now compute the RMS error between the points */
	RMS_sum = 0.0;

	for (i = 0; i < n; i++) {
		dvec3 r = (right_centroid - right_pts[i]);
		dvec3 l = (left_centroid - left_pts[i]);
		dvec3 resid;

		/* Rotate, scale l */
		dvec3 Rl, SRl;

		Rl.x = R[0] * l.x + R[1] * l.y + R[2] * l.z;
		Rl.y = R[3] * l.x + R[4] * l.y + R[5] * l.z;
		Rl.z = R[6] * l.x + R[7] * l.y + R[8] * l.z;

		SRl = (*scale) * Rl;

		resid = r - SRl;
		RMS_sum += resid.length2();
	}

	return sqrt(RMS_sum / n);
}

/* Computes the closed-form least-squares solution to a rigid
* body alignment.
*
* n: the number of points
* right_pts: Target set of n points
* left_pts:  Source set of n points */
double align_horn_3D(int n, dvec3 *right_pts, dvec3 *left_pts, int scale_xform,
	double *Tout) {
	int i;
	dvec3 right_centroid(0.0, 0.0, 0.0);
	dvec3 left_centroid(0.0, 0.0, 0.0);
	double M[3][3] = { { 0.0, 0.0, 0.0, },
	{ 0.0, 0.0, 0.0, },
	{ 0.0, 0.0, 0.0, } };
	double MT[3][3];
	double MTM[3][3];
	double eval[3], sqrteval_inv[3];
	double evec[3][3], evec_tmp[3][3];
	double Sinv[3][3], U[3][3];
	double Tcenter[4][4] = { { 1.0, 0.0, 0.0, 0.0 },
	{ 0.0, 1.0, 0.0, 0.0 },
	{ 0.0, 0.0, 1.0, 0.0 },
	{ 0.0, 0.0, 0.0, 1.0 } };

	double Ttmp[4][4];
	double T[16], R[16];

	double sum_num, sum_den, scale, RMS_sum;

	int perm[3];

	/* Compute the centroid of both point sets */
	right_centroid = geom::vec_mean(n, right_pts);
	left_centroid = geom::vec_mean(n, left_pts);

	/* Compute the scale */
	sum_num = sum_den = 0.0;

	for (i = 0; i < n; i++) {
		dvec3 r = (right_centroid - right_pts[i]);
		dvec3 l = (left_centroid - left_pts[i]);

		sum_num += r.length2();
		sum_den += l.length2();
	}

	scale = sqrt(sum_num / sum_den);

	/* Fill in the matrix M */
	for (i = 0; i < n; i++) {
		dvec3 r = (right_centroid - right_pts[i]);
		dvec3 l = (left_centroid - left_pts[i]);

		M[0][0] += r.x * l.x;
		M[0][1] += r.x * l.y;
		M[0][2] += r.x * l.z;

		M[1][0] += r.y * l.x;
		M[1][1] += r.y * l.y;
		M[1][2] += r.y * l.z;

		M[2][0] += r.z * l.x;
		M[2][1] += r.z * l.y;
		M[2][2] += r.z * l.z;
	}

	/* Compute MTM */
	matrix_transpose(3, 3, (double *)M, (double *)MT);
	matrix_product(3, 3, 3, 3, (double *)MT, (double *)M, (double *)MTM);

	/* Calculate Sinv, the inverse of the square root of MTM */
	dgeev_driver(3, (double *)MTM, (double *)evec, eval);

	/* Sort the eigenvalues */
	qsort_descending();
	qsort_perm(3, eval, perm);

	memcpy(evec_tmp[0], evec[perm[0]], sizeof(double) * 3);
	memcpy(evec_tmp[1], evec[perm[1]], sizeof(double) * 3);
	memcpy(evec_tmp[2], evec[perm[2]], sizeof(double) * 3);
	memcpy(evec, evec_tmp, sizeof(double) * 9);

	sqrteval_inv[0] = 1.0 / sqrt(eval[0]);
	sqrteval_inv[1] = 1.0 / sqrt(eval[1]);

	if (eval[2] < 1.0e-8 * eval[0]) {
		sqrteval_inv[2] = 0.0;
	}
	else {
		sqrteval_inv[2] = 1.0 / sqrt(eval[2]);
	}

	Sinv[0][0] =
		sqrteval_inv[0] * evec[0][0] * evec[0][0] +
		sqrteval_inv[1] * evec[1][0] * evec[1][0] +
		sqrteval_inv[2] * evec[2][0] * evec[2][0];
	Sinv[0][1] =
		sqrteval_inv[0] * evec[0][0] * evec[0][1] +
		sqrteval_inv[1] * evec[1][0] * evec[1][1] +
		sqrteval_inv[2] * evec[2][0] * evec[2][1];
	Sinv[0][2] =
		sqrteval_inv[0] * evec[0][0] * evec[0][2] +
		sqrteval_inv[1] * evec[1][0] * evec[1][2] +
		sqrteval_inv[2] * evec[2][0] * evec[2][2];

	Sinv[1][0] =
		sqrteval_inv[0] * evec[0][1] * evec[0][0] +
		sqrteval_inv[1] * evec[1][1] * evec[1][0] +
		sqrteval_inv[2] * evec[2][1] * evec[2][0];
	Sinv[1][1] =
		sqrteval_inv[0] * evec[0][1] * evec[0][1] +
		sqrteval_inv[1] * evec[1][1] * evec[1][1] +
		sqrteval_inv[2] * evec[2][1] * evec[2][1];
	Sinv[1][2] =
		sqrteval_inv[0] * evec[0][1] * evec[0][2] +
		sqrteval_inv[1] * evec[1][1] * evec[1][2] +
		sqrteval_inv[2] * evec[2][1] * evec[2][2];

	Sinv[2][0] =
		sqrteval_inv[0] * evec[0][2] * evec[0][0] +
		sqrteval_inv[1] * evec[1][2] * evec[1][0] +
		sqrteval_inv[2] * evec[2][2] * evec[2][0];
	Sinv[2][1] =
		sqrteval_inv[0] * evec[0][2] * evec[0][1] +
		sqrteval_inv[1] * evec[1][2] * evec[1][1] +
		sqrteval_inv[2] * evec[2][2] * evec[2][1];
	Sinv[2][2] =
		sqrteval_inv[0] * evec[0][2] * evec[0][2] +
		sqrteval_inv[1] * evec[1][2] * evec[1][2] +
		sqrteval_inv[2] * evec[2][2] * evec[2][2];

	/* U = M * Sinv */
	matrix_product(3, 3, 3, 3, (double *)M, (double *)Sinv, (double *)U);

	if (eval[2] < 1.0e-8 * eval[0]) {
		double u3u3[9], Utmp[9];
		//matrix_transpose_product2(3, 1, 3, 1, evec[2], evec[2], u3u3);
		matrix_transpose_product(3, 1, 3, 1, evec[2], evec[2], u3u3);

		matrix_sum(3, 3, 3, 3, (double *)U, u3u3, Utmp);

		if (matrix_determinant3(Utmp) < 0.0) {
			printf("[align_horn_3D] Recomputing matrix...\n");
			matrix_diff(3, 3, 3, 3, (double *)U, u3u3, Utmp);
		}

		memcpy(U, Utmp, 9 * sizeof(double));
	}

	/* Fill in the rotation matrix */
	R[0] = U[0][0]; R[1] = U[0][1]; R[2] = U[0][2]; R[3] = 0.0;
	R[4] = U[1][0]; R[5] = U[1][1]; R[6] = U[1][2]; R[7] = 0.0;
	R[8] = U[2][0]; R[9] = U[2][1]; R[10] = U[2][2]; R[11] = 0.0;
	R[12] = 0.0;     R[13] = 0.0;     R[14] = 0.0;     R[15] = 1.0;

	/* Fill in the translation matrix */
	matrix_ident(4, T);
	T[3] = right_centroid.x;
	T[7] = right_centroid.y;
	T[11] = right_centroid.z;

	if (scale_xform == 0)
		scale = 1.0;

	Tcenter[0][0] = scale;
	Tcenter[1][1] = scale;
	Tcenter[2][2] = scale;

	Tcenter[0][3] = -scale * left_centroid.x;
	Tcenter[1][3] = -scale * left_centroid.y;
	Tcenter[2][3] = -scale * left_centroid.z;

	matrix_product(4, 4, 4, 4, T, R, (double *)Ttmp);
	matrix_product(4, 4, 4, 4, (double *)Ttmp, (double *)Tcenter, Tout);

	/* Now compute the RMS error between the points */
	RMS_sum = 0.0;

	for (i = 0; i < n; i++) {
		double left[4] = { 
			left_pts[i].x,
			left_pts[i].y,
			left_pts[i].z, 1.0 };
		double left_prime[3];
		double dx, dy, dz;

		matrix_product(4, 4, 4, 1, Tout, left, left_prime);

		dx = left_prime[0] - right_pts[i].x;
		dy = left_prime[1] - right_pts[i].y;
		dz = left_prime[2] - right_pts[i].z;

		RMS_sum += dx * dx + dy * dy + dz * dz;
	}

	return sqrt(RMS_sum / n);
}


/* Computes the closed-form least-squares solution to a rigid
* body alignment.
*
* n: the number of points
* right_pts: Target set of n points
* left_pts:  Source set of n points */
double align_horn_3D_2(int n, dvec3 *right_pts, dvec3 *left_pts, int scale_xform,
	double *Tout)
{
	int i;
	dvec3 right_centroid(0.0, 0.0, 0.0);
	dvec3 left_centroid(0.0, 0.0, 0.0);
	double Tcenter[16] = { 1.0, 0.0, 0.0, 0.0,
		0.0, 1.0, 0.0, 0.0,
		0.0, 0.0, 1.0, 0.0,
		0.0, 0.0, 0.0, 1.0 };

	double Ttmp[4][4];
	double T[16], R[16], R3x3[9];

	double sum_num, sum_den, scale, RMS_sum;

	dvec3 *left_pts_zm = new dvec3[n];
	dvec3 *right_pts_zm = new dvec3[n];

	double error = 0.0;

	/* Compute the centroid of both point sets */
	right_centroid = geom::vec_mean(n, right_pts);
	left_centroid = geom::vec_mean(n, left_pts);

	/* Compute the scale */
	sum_num = sum_den = 0.0;

	for (i = 0; i < n; i++) {
		dvec3 r = right_centroid - right_pts[i];
		dvec3 l = left_centroid - left_pts[i];

		sum_num += r.length2();
		sum_den += l.length2();
	}

	scale = sqrt(sum_num / sum_den);

	for (i = 0; i < n; i++) {
		dvec3 r = right_centroid - right_pts[i];
		dvec3 l = left_centroid - left_pts[i];

		right_pts_zm[i] = r;
		left_pts_zm[i] = scale * l;
	}

	/* Compute the rotation */
	error = align_3D_rotation(n, right_pts_zm, left_pts_zm, R3x3);
	// printf("error[%d]: %0.3f\n", n, error);
	// matrix_print(3, 3, R3x3);

	/* Fill in the rotation matrix */
	R[0] = R3x3[0]; R[1] = R3x3[1]; R[2] = R3x3[2]; R[3] = 0.0;
	R[4] = R3x3[3]; R[5] = R3x3[4]; R[6] = R3x3[5]; R[7] = 0.0;
	R[8] = R3x3[6]; R[9] = R3x3[7]; R[10] = R3x3[8]; R[11] = 0.0;
	R[12] = 0.0;     R[13] = 0.0;     R[14] = 0.0;     R[15] = 1.0;

	/* Fill in the translation matrix */
	// matrix_ident(4, T);
	T[0] = 1.0;  T[1] = 0.0;  T[2] = 0.0;  T[3] = right_centroid.x;
	T[4] = 0.0;  T[5] = 1.0;  T[6] = 0.0;  T[7] = right_centroid.y;
	T[8] = 0.0;  T[9] = 0.0;  T[10] = 1.0; T[11] = right_centroid.z;
	T[12] = 0.0; T[13] = 0.0; T[14] = 0.0; T[15] = 1.0;

	if (scale_xform == 0)
		scale = 1.0;

	Tcenter[0] = scale;
	Tcenter[5] = scale;
	Tcenter[10] = scale;

	Tcenter[3] = -scale * left_centroid.x;
	Tcenter[7] = -scale * left_centroid.y;
	Tcenter[11] = -scale * left_centroid.z;

	matrix_product44(T, R, (double *)Ttmp);
	matrix_product44((double *)Ttmp, (double *)Tcenter, Tout);

	/* Now compute the RMS error between the points */
	RMS_sum = 0.0;

	for (i = 0; i < n; i++) {
		double left[4] = { 
			left_pts[i].x,
			left_pts[i].y,
			left_pts[i].z, 1.0 };
		double left_prime[3];
		double dx, dy, dz;

		matrix_product441(Tout, left, left_prime);

		dx = left_prime[0] - right_pts[i].x;
		dy = left_prime[1] - right_pts[i].y;
		dz = left_prime[2] - right_pts[i].z;

		RMS_sum += dx * dx + dy * dy + dz * dz;
	}

	delete[] (left_pts_zm);
	delete[] (right_pts_zm);

	return sqrt(RMS_sum / n);
}

/* Align two sets of points with a 3D rotation */
double align_3D_rotation(int n, dvec3 *r_pts, dvec3 *l_pts, double *R)
{
	double A[9];
	double U[9], S[3], V[9], VT[9], RT[9];
	int i;
	double error;

	for (i = 0; i < 9; i++)
		A[i] = 0.0;

	for (i = 0; i < n; i++) {
		double *a = l_pts[i].data(), *b = r_pts[i].data();
		// matrix_product(3, 1, 1, 3, l_pts[i].p, r_pts[i].p, tensor);
		A[0] += a[0] * b[0];
		A[1] += a[0] * b[1];
		A[2] += a[0] * b[2];

		A[3] += a[1] * b[0];
		A[4] += a[1] * b[1];
		A[5] += a[1] * b[2];

		A[6] += a[2] * b[0];
		A[7] += a[2] * b[1];
		A[8] += a[2] * b[2];
	}

	svd(3, 3, 1, 1, 1.0e-12, 1.0e-12, A, S, U, V, VT);

	matrix_product33(U, VT, RT);
	matrix_transpose(3, 3, RT, R);

	// printf("R:\n");
	// matrix_print(3, 3, R);

	if (matrix_determinant3(R) < 0.0) {
		/* We're dealing with a reflection */
		double tmp[9];
		double reflectZ[9] = { 1.0, 0.0, 0.0,
			0.0, 1.0, 0.0,
			0.0, 0.0, -1.0 };

		matrix_product33(U, reflectZ, tmp);
		matrix_product33(tmp, VT, RT);
		matrix_transpose(3, 3, RT, R);
	}

	/* Compute error */
	error = 0.0;
	for (i = 0; i < n; i++) {
		double rot[3];
		double diff[3];
		double dist;
		matrix_product331(R, l_pts[i].data(), rot);
		matrix_diff(3, 1, 3, 1, rot, r_pts[i].data(), diff);
		dist = matrix_norm(3, 1, diff);

		// printf("d[%d] = %0.6f\n", i, dist);
		error += dist;
	}

	return error / n;
}

double align_2D(int n, dvec3 *right_pts, dvec3 *left_pts,
	double *R, double *T,
	double *Tout, double *scale, double *weight)
{
	int i;
	dvec3 right_centroid = geom::vec_mean(n, right_pts);
	dvec3 left_centroid = geom::vec_mean(n, left_pts);
	double Ttmp[3][3], Tcenter[3][3];

	/* Setup the matrix */
	int num_vars = 2;
	int num_eqns = 2 * n;
	double *A = (double *)malloc(sizeof(double) * num_vars * num_eqns);
	double *b = (double *)malloc(sizeof(double) * num_eqns);
	double x[2];

	double RMS_sum = 0.0;

	assert(n >= 2);
	for (i = 0; i < n; i++) {
		double *row1 = A + 2 * i * num_vars;
		double *row2 = A + (2 * i + 1) * num_vars;

		row1[0] = left_pts[i].x - left_centroid.x;
		row1[1] = left_pts[i].y - left_centroid.y;
		b[2 * i + 0] = right_pts[i].x - right_centroid.x;

		row2[0] = left_pts[i].y - left_centroid.y;
		row2[1] = -(left_pts[i].x - left_centroid.x);
		b[2 * i + 1] = right_pts[i].y - right_centroid.y;
	}

	/* Solve the system */
	dgelsy_driver(A, b, x, num_eqns, num_vars, 1);

	/* Fill in the rotation matrix */
	R[0] = x[0];  R[1] = x[1]; R[2] = 0.0;
	R[3] = -x[1]; R[4] = x[0]; R[5] = 0.0;
	R[6] = 0.0;   R[7] = 0.0;  R[8] = 1.0;

	/* Fill in the translation matrix */
	matrix_ident(3, T);
	T[2] = right_centroid.x;
	T[5] = right_centroid.y;

	matrix_ident(3, (double *)Tcenter);
	Tcenter[0][2] = -left_centroid.x;
	Tcenter[1][2] = -left_centroid.y;

	matrix_product(3, 3, 3, 3, T, R, (double *)Ttmp);
	matrix_product(3, 3, 3, 3, (double *)Ttmp, (double *)Tcenter, Tout);

	T[2] = (right_centroid - left_centroid).x;
	T[5] = (right_centroid - left_centroid).y;

	/* Now compute the RMS error between the points */
	RMS_sum = 0.0;

	for (i = 0; i < n; i++) {
		dvec3 r = (right_pts[i] - right_centroid);
		dvec3 l = (left_pts[i] - left_centroid);
		dvec3 resid;

		/* Rotate, scale l */
		dvec3 SRl;

		SRl.x = R[0] * l.x + R[1] * l.y + R[2] * l.z;
		SRl.y = R[3] * l.x + R[4] * l.y + R[5] * l.z;
		SRl.z = R[6] * l.x + R[7] * l.y + R[8] * l.z;

		// SRl = v3_scale(*scale, Rl);

		resid = r - SRl;
		RMS_sum += resid.length2();
	}

	free(A);
	free(b);

	return sqrt(RMS_sum / n);
}

/* Align two sets of points with a 2D similarity transform */
int align_2D_ransac(int n, dvec3 *r_pts, dvec3 *l_pts,
	int num_ransac_rounds, double ransac_thresh,
	double *Tret)
{
	int round;
#define MIN_SUPPORT 2
	dvec3 *l_inliers, *r_inliers;
	int num_inliers, max_inliers = 0;
	double Tbest[9];

	if (n < 3) {
		printf("[align_2D_ransac] Error: need at least 3 points!\n");
		return 0;
	}

	l_inliers = (dvec3 *)malloc(sizeof(dvec3) * n);
	r_inliers = (dvec3 *)malloc(sizeof(dvec3) * n);

	for (round = 0; round < num_ransac_rounds; round++) {
		int support[MIN_SUPPORT];
		int i, j;
		dvec3 r_mean, l_mean, r0, l0;
		dvec3 r_pts_small[MIN_SUPPORT], l_pts_small[MIN_SUPPORT];
		double Rtmp[9], T1tmp[9], T2tmp[9], tmp[9], Tout[9];
		double a, b;

		for (i = 0; i < MIN_SUPPORT; i++) {
			/* Select an index from 0 to n-1 */
			int idx, reselect;

			do {
				reselect = 0;
				idx = rand() % n;
				for (j = 0; j < i; j++) {
					if (support[j] == idx) {
						reselect = 1;
						break;
					}
				}
			} while (reselect);

			support[i] = idx;
			r_pts_small[i] = r_pts[idx];
			l_pts_small[i] = l_pts[idx];
		}

		r_mean = 0.5 * (r_pts_small[0] + r_pts_small[1]);
		l_mean = 0.5 * (l_pts_small[0] + l_pts_small[1]);

		r0 = r_pts_small[0] - r_mean;
		l0 = l_pts_small[0] - l_mean;

		a = (r0.y + r0.x * l0.x / l0.y) /
			(l0.y + l0.x * l0.x / l0.y);
		b = (r0.x - a * l0.x) / l0.y;

		Rtmp[0] = a;  Rtmp[1] = b;  Rtmp[2] = 0.0;
		Rtmp[3] = -b; Rtmp[4] = a;  Rtmp[5] = 0.0;
		Rtmp[6] = 0;  Rtmp[7] = 0;  Rtmp[8] = 1.0;

		matrix_ident(3, T1tmp);
		T1tmp[2] = -l_mean.x;
		T1tmp[5] = -l_mean.y;

		matrix_ident(3, T2tmp);
		T2tmp[2] = r_mean.x;
		T2tmp[5] = r_mean.y;

		matrix_product(3, 3, 3, 3, Rtmp, T1tmp, tmp);
		matrix_product(3, 3, 3, 3, T2tmp, tmp, Tout);

		/* Count inliers */
		num_inliers = 0;
		for (i = 0; i < n; i++) {
			double Tp[3];
			double diff[3];
			double dist;
			matrix_product(3, 3, 3, 1, Tout, l_pts[i].data(), Tp);
			matrix_diff(3, 1, 3, 1, Tp, r_pts[i].data(), diff);
			dist = matrix_norm(3, 1, diff);

			if (dist < ransac_thresh) {
				num_inliers++;
			}

			if (num_inliers > max_inliers) {
				max_inliers = num_inliers;
				memcpy(Tbest, Tout, sizeof(double) * 9);
				// printf(" inliers_new: %d\n", num_inliers);
			}
		}
	}

	memcpy(Tret, Tbest, 9 * sizeof(double));

	free(r_inliers);
	free(l_inliers);

	return num_inliers;
#undef MIN_SUPPORT
}

/* Align two sets of points with a 2D similarity transform */
int align_horn_ransac(int n, dvec3 *r_pts, dvec3 *l_pts,
	int num_ransac_rounds, double ransac_thresh,
	double *Tret)
{
	int round;
#define MIN_SUPPORT 3
	dvec3 *l_inliers, *r_inliers;
	int num_inliers, max_inliers = 0;
	double Tbest[9];

	if (n < 3) {
		printf("[align_horn_ransac] Error: need at least 3 points!\n");
		return 0;
	}

	l_inliers = (dvec3 *)malloc(sizeof(dvec3) * n);
	r_inliers = (dvec3 *)malloc(sizeof(dvec3) * n);

	for (round = 0; round < num_ransac_rounds; round++) {
		int support[MIN_SUPPORT];
		int i, j;
		dvec3 r_pts_small[MIN_SUPPORT], l_pts_small[MIN_SUPPORT];
		double Rtmp[9], Ttmp[9], Tout[9], scale_tmp;

		for (i = 0; i < MIN_SUPPORT; i++) {
			/* Select an index from 0 to n-1 */
			int idx, reselect;
			do {
				reselect = 0;
				idx = rand() % n;
				for (j = 0; j < i; j++) {
					if (support[j] == idx) {
						reselect = 1;
						break;
					}
				}
			} while (reselect);

			support[i] = idx;
			r_pts_small[i] = r_pts[idx];
			l_pts_small[i] = l_pts[idx];
		}

		align_horn(MIN_SUPPORT, r_pts_small, l_pts_small, Rtmp, Ttmp,
			Tout, &scale_tmp, NULL);

		/* Count inliers */
		num_inliers = 0;
		for (i = 0; i < n; i++) {
			double Tp[3];
			double diff[3];
			double dist;
			matrix_product(3, 3, 3, 1, Tout, l_pts[i].data(), Tp);
			matrix_diff(3, 1, 3, 1, Tp, r_pts[i].data(), diff);
			dist = matrix_norm(3, 1, diff);

			if (dist < ransac_thresh) {
				num_inliers++;
			}

			if (num_inliers > max_inliers) {
				max_inliers = num_inliers;
				memcpy(Tbest, Tout, sizeof(double) * 9);
				// printf(" inliers_new: %d\n", num_inliers);
			}
		}
	}

#if 0
	/* Reestimate using all inliers */
	num_inliers = 0;
	for (i = 0; i < n; i++) {
		double Tp[3];
		double diff[3];
		double dist;
		matrix_product(3, 3, 3, 1, Tbest, l_pts[i].p, Tp);
		matrix_diff(3, 1, 3, 1, Tp, r_pts[i].p, diff);
		dist = matrix_norm(3, 1, diff);

		if (dist < ransac_thresh) {
			r_inliers[num_inliers] = r_pts[i];
			l_inliers[num_inliers] = l_pts[i];
			num_inliers++;
		}
	}

	// printf(" inliers: %d\n", num_inliers);

	align_horn(num_inliers, r_inliers, l_inliers, R, T, Tret, &scale, NULL);
#else
	memcpy(Tret, Tbest, 9 * sizeof(double));
#endif

	free(r_inliers);
	free(l_inliers);

	return num_inliers;
#undef MIN_SUPPORT
}

/* Align two sets of points with a 3D similarity transform */
int align_horn_3D_ransac(int n, dvec3 *r_pts, dvec3 *l_pts,
	int num_ransac_rounds, double ransac_thresh,
	double *Tret)
{
	int round;
#define MIN_SUPPORT 3
	dvec3 *l_inliers, *r_inliers;
	double *Vp, *TVp;
	int num_inliers, max_inliers = 0;
	double Tbest[16];
	int i;
	double *ptr;

	double ransac_threshsq = ransac_thresh * ransac_thresh;

	if (n < MIN_SUPPORT) {
		printf("[align_horn_3D_ransac] Error: need at least %d points!\n",
			MIN_SUPPORT);
		return 0;
	}

	l_inliers = (dvec3 *)malloc(sizeof(dvec3) * n);
	r_inliers = (dvec3 *)malloc(sizeof(dvec3) * n);

	Vp = (double *)malloc(sizeof(double) * 4 * n);
	TVp = (double *)malloc(sizeof(double) * 4 * n);

	for (i = 0; i < n; i++) {
		memcpy(Vp + 4 * i, l_pts[i].data(), 3 * sizeof(double));
		Vp[4 * i + 3] = 1.0;
	}

	for (round = 0; round < num_ransac_rounds; round++) {
		int support[MIN_SUPPORT];
		int i, j;
		dvec3 r_pts_small[MIN_SUPPORT], l_pts_small[MIN_SUPPORT];
		double Tout[16], ToutT[16];
		int nan = 0;

		for (i = 0; i < MIN_SUPPORT; i++) {
			/* Select an index from 0 to n-1 */
			int idx, reselect;
			do {
				reselect = 0;
				idx = rand() % n;
				for (j = 0; j < i; j++) {
					if (support[j] == idx) {
						reselect = 1;
						break;
					}
				}
			} while (reselect);

			support[i] = idx;
			r_pts_small[i] = r_pts[idx];
			l_pts_small[i] = l_pts[idx];
		}

		align_horn_3D_2(MIN_SUPPORT, r_pts_small, l_pts_small, 1, Tout);

#if 1
		for (i = 0; i < 16; i++) {
			if (isnan(Tout[i]) || Tout[i] != Tout[i]) {
				nan = 1;
				break;
			}
		}

		if (nan == 1)
			continue;
#endif

		/* Count inliers */
		num_inliers = 0;

#if 0
		for (i = 0; i < n; i++) {
			double Tp[4];
			double diff[3];
			double dist;

			double p[4] = { l_pts[i].p[0], l_pts[i].p[1], l_pts[i].p[2], 1.0 };

			matrix_product(4, 4, 4, 1, Tout, p, Tp);
			matrix_diff(3, 1, 3, 1, Tp, r_pts[i].p, diff);
			dist = matrix_norm(3, 1, diff);

			if (dist < ransac_thresh) {
				num_inliers++;
			}
		}
#else
		matrix_transpose(4, 4, Tout, ToutT);
		matrix_product(n, 4, 4, 4, Vp, ToutT, TVp);
		//matrix_product_old(n, 4, 4, 4, Vp, ToutT, TVp);

		ptr = TVp;

		for (i = 0; i < n; i++) {
			// double diff[3], dist;
			double dx, dy, dz, dist;
			// matrix_diff(3, 1, 3, 1, TVp + 4 * i, r_pts[i].p, diff);
			dx = ptr[0] - r_pts[i].x;
			dy = ptr[1] - r_pts[i].y;
			dz = ptr[2] - r_pts[i].z;

			dist = dx * dx + dy * dy + dz * dz; // matrix_normsq(3, 1, diff);

			if (dist < ransac_threshsq)
				num_inliers++;

			ptr += 4;
		}
#endif

		if (num_inliers > max_inliers) {
			max_inliers = num_inliers;
			memcpy(Tbest, Tout, sizeof(double) * 16);
		}
	}

	/* Reestimate using all inliers */
#if 0
	matrix_transpose(4, 4, Tbest, TbestT);
	matrix_product(n, 4, 4, 4, Vp, TbestT, TVp);

	num_inliers = 0;
	for (i = 0; i < n; i++) {
		// double Tp[4];
		double diff[3];
		double dist;

		matrix_diff(3, 1, 3, 1, TVp + 4 * i, r_pts[i].p, diff);

		// double p[4] = { l_pts[i].p[0], l_pts[i].p[1], l_pts[i].p[2], 1.0 };

		// matrix_product(4, 4, 4, 1, Tbest, p, Tp);
		// matrix_diff(3, 1, 3, 1, Tp, r_pts[i].p, diff);
		dist = matrix_normsq(3, 1, diff);

		if (dist < ransac_threshsq) {
			r_inliers[num_inliers] = r_pts[i];
			l_inliers[num_inliers] = l_pts[i];
			num_inliers++;
		}
	}

	align_horn_3D_2(num_inliers, r_inliers, l_inliers, 1, Tret);

	// memcpy(Tret, Tbest, 16 * sizeof(double));

	if (isnan(Tret[0]) || Tret[0] != Tret[0]) {
		printf("[align_horn_3D_ransac] nan at end [num_inliers: %d], "
			"restoring old matrix\n", num_inliers);
		memcpy(Tret, Tbest, sizeof(double) * 16);
	}
#else
	memcpy(Tret, Tbest, sizeof(double) * 16);
#endif

	free(r_inliers);
	free(l_inliers);
	free(Vp);
	free(TVp);

	return max_inliers;
#undef MIN_SUPPORT
}

/* Align two sets of points with a 3D rotation */
int align_3D_rotation_ransac(int n, dvec3 *r_pts, dvec3 *l_pts,
	int num_ransac_rounds, double ransac_thresh,
	double *R)
{
	int round;
	double error = 0.0;
#define MIN_SUPPORT 3
	// const int min_support = 3;
	dvec3 *l_inliers, *r_inliers;
	int num_inliers, max_inliers = 0;
	double Rbest[9];

	if (n < 3) {
		printf("[align_3D_rotation_ransac] Error: need at least 3 points!\n");
		return 0;
	}

	l_inliers = (dvec3 *)malloc(sizeof(dvec3) * n);
	r_inliers = (dvec3 *)malloc(sizeof(dvec3) * n);

	for (round = 0; round < num_ransac_rounds; round++) {
		int support[MIN_SUPPORT];
		int i, j;
		dvec3 r_pts_small[MIN_SUPPORT], l_pts_small[MIN_SUPPORT];
		double Rtmp[9];

		for (i = 0; i < MIN_SUPPORT; i++) {
			/* Select an index from 0 to n-1 */
			int idx, reselect;
			do {
				reselect = 0;
				idx = rand() % n;
				for (j = 0; j < i; j++) {
					if (support[j] == idx) {
						reselect = 1;
						break;
					}
				}
			} while (reselect);

			support[i] = idx;
			r_pts_small[i] = r_pts[idx];
			l_pts_small[i] = l_pts[idx];
		}

		align_3D_rotation(MIN_SUPPORT, r_pts_small, l_pts_small, Rtmp);

		/* Count inliers */
		num_inliers = 0;

		for (i = 0; i < n; i++) {
			double rot[3];
			double diff[3];
			double dist;
			matrix_product(3, 3, 3, 1, Rtmp, l_pts[i].data(), rot);
			matrix_diff(3, 1, 3, 1, rot, r_pts[i].data(), diff);
			dist = matrix_norm(3, 1, diff);

			if (dist < ransac_thresh) {
				num_inliers++;
			}

			if (num_inliers > max_inliers) {
				max_inliers = num_inliers;
				memcpy(Rbest, Rtmp, sizeof(double) * 9);
			}
		}
	}

#if 0
	/* Reestimate using all inliers */
	num_inliers = 0;
	for (i = 0; i < n; i++) {
		double rot[3];
		double diff[3];
		double dist;
		matrix_product(3, 3, 3, 1, Rbest, l_pts[i].p, rot);
		matrix_diff(3, 1, 3, 1, rot, r_pts[i].p, diff);
		dist = matrix_norm(3, 1, diff);

		if (dist < ransac_thresh) {
			r_inliers[num_inliers] = r_pts[i];
			l_inliers[num_inliers] = l_pts[i];
			num_inliers++;
		}
	}

	error = align_3D_rotation(num_inliers, r_inliers, l_inliers, R);

	printf("[align_3D_rotation] Error: %0.3f\n", error);

	free(r_inliers);
	free(l_inliers);

	return num_inliers;
#else
	memcpy(R, Rbest, 9 * sizeof(double));

	free(r_inliers);
	free(l_inliers);

	return max_inliers;
#endif

#undef MIN_SUPPORT
}
