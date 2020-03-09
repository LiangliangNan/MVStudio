#include "triangulate.h"

#include "../math/matrix_driver.h"
#include "../math/svd.h"


static vec2d global_p, global_q;
static double *global_R0, *global_t0, *global_R1, *global_t1;

void quick_svd(double *E, double *U, double *S, double *VT) {
	double e1[3] = { E[0], E[3], E[6] };
	double e2[3] = { E[1], E[4], E[7] };
	double e3[3] = { E[2], E[5], E[8] };

	double e1e2[3], e1e3[3], e2e3[3];
	double e1e2_mag, e1e3_mag, e2e3_mag;

	double u_a[3], u_b[3], u_c[3];
	double v_a[3], v_b[3], v_c[3];
	double tmp[3], mag;

	double e1_mag = matrix_norm(3, 1, e1);
	double e2_mag = matrix_norm(3, 1, e2);

	matrix_cross(e1, e2, e1e2);
	matrix_cross(e1, e3, e1e3);
	matrix_cross(e2, e3, e2e3);

	e1e2_mag = matrix_norm(3, 1, e1e2);
	e1e3_mag = matrix_norm(3, 1, e1e3);
	e2e3_mag = matrix_norm(3, 1, e2e3);

	if (e1e2_mag >= e1e3_mag && e1e2_mag >= e1e3_mag) {
		matrix_scale(3, 1, e1, 1.0 / e1_mag, v_a);
		matrix_scale(3, 1, e1e2, 1.0 / e1e2_mag, v_c);
		matrix_cross(v_c, v_a, v_b);
	}
	else if (e1e3_mag > e1e2_mag && e1e3_mag > e2e3_mag) {
		matrix_scale(3, 1, e1, 1.0 / e1_mag, v_a);
		matrix_scale(3, 1, e1e3, 1.0 / e1e3_mag, v_b);
		matrix_cross(v_b, v_a, v_c);
	}
	else if (e2e3_mag > e1e2_mag && e2e3_mag > e1e3_mag) {
		matrix_scale(3, 1, e2, 1.0 / e2_mag, v_a);
		matrix_scale(3, 1, e2e3, 1.0 / e2e3_mag, v_a);
		matrix_cross(v_b, v_c, v_a);
	}

	matrix_product331(E, v_a, tmp);
	mag = matrix_norm(3, 1, tmp);
	matrix_scale(3, 1, tmp, 1.0 / mag, u_a);

	matrix_product331(E, v_b, tmp);
	mag = matrix_norm(3, 1, tmp);
	matrix_scale(3, 1, tmp, 1.0 / mag, u_b);

	matrix_cross(u_a, u_b, u_c);
}

/* Project a point onto an image */
vec2d project(double *R, double *t0, double *P) {
	double tmp[3], tmp2[3];
	vec2d result;

	/* Rigid transform */
	matrix_product331(R, P, tmp);
	matrix_sum(3, 1, 3, 1, tmp, t0, tmp2);

	/* Perspective division */
	result.x = tmp2[0] / tmp2[2];
	result.y = tmp2[1] / tmp2[2];

	return result;
}

void triangulation_residual(const int *m, const int *n, double *x,
	double *fvec, double *iflag)
{
	/* Project the point into the two views */
	vec2d p = project(global_R0, global_t0, x);
	vec2d q = project(global_R1, global_t1, x);

	fvec[0] = global_p.x - p.x;
	fvec[1] = global_p.y - p.y;
	fvec[2] = global_q.x - q.x;
	fvec[3] = global_q.y - q.y;
}

static int global_num_points;
static double *global_Rs = NULL;
static double *global_ts = NULL;
static vec2d *global_ps;

void triangulate_n_residual(const int *m, const int *n,
	double *x, double *fvec, double *iflag)
{
	int i;

	for (i = 0; i < global_num_points; i++) {
		int Roff = 9 * i;
		int toff = 3 * i;

		/* Project the point into the view */
		vec2d p = project(global_Rs + Roff, global_ts + toff, x);

		fvec[2 * i + 0] = global_ps[i].x - p.x;
		fvec[2 * i + 1] = global_ps[i].y - p.y;
	}
}


/* Find the point with the smallest squared projection error */
vec3d triangulate_n_refine(vec3d pt, int num_points,
	vec2d *p, double *R, double *t, double *error_out)
{
	int num_eqs = 2 * num_points;
	int num_vars = 3;

	double x[3] = { pt.x, pt.y, pt.z };
	double error;

	int i;

	/* Run a non-linear optimization to polish the result */
	global_num_points = num_points;
	global_ps = p;
	global_Rs = R;  global_ts = t;
	lmdif_driver((void*)triangulate_n_residual, num_eqs, num_vars, x, 1.0e-5);

	error = 0.0;
	for (i = 0; i < num_points; i++) {
		double dx, dy;
		int Roff = 9 * i;
		int toff = 3 * i;
		double pp[3];

		/* Compute projection error */
		matrix_product331(R + Roff, x, pp);
		pp[0] += t[toff + 0];
		pp[1] += t[toff + 1];
		pp[2] += t[toff + 2];

		dx = pp[0] / pp[2] - p[i].x;
		dy = pp[1] / pp[2] - p[i].y;
		error += dx * dx + dy * dy;
	}

	error = sqrt(error / num_points);

	// printf("[triangulate_n] Error [after polishing]: %0.3e\n", error);

	if (error_out != NULL) {
		*error_out = error;
	}

	return vec3d(x[0], x[1], x[2]);
}


/* Find the point with the smallest squared projection error */
vec3d triangulate_n(int num_points,
	vec2d *p, double *R, double *t, double *error_out)
{
	int num_eqs = 2 * num_points;
	int num_vars = 3;

	double *A = (double *)malloc(sizeof(double) * num_eqs * num_vars);
	double *b = (double *)malloc(sizeof(double) * num_eqs);
	double *x = (double *)malloc(sizeof(double) * num_vars);

	int i;
	double error;

	vec3d r;

	for (i = 0; i < num_points; i++) {
		int Roff = 9 * i;
		int row = 6 * i;
		int brow = 2 * i;
		int toff = 3 * i;

		A[row + 0] = R[Roff + 0] - p[i].x * R[Roff + 6];
		A[row + 1] = R[Roff + 1] - p[i].x * R[Roff + 7];
		A[row + 2] = R[Roff + 2] - p[i].x * R[Roff + 8];

		A[row + 3] = R[Roff + 3] - p[i].y * R[Roff + 6];
		A[row + 4] = R[Roff + 4] - p[i].y * R[Roff + 7];
		A[row + 5] = R[Roff + 5] - p[i].y * R[Roff + 8];

		b[brow + 0] = t[toff + 2] * p[i].x - t[toff + 0];
		b[brow + 1] = t[toff + 2] * p[i].y - t[toff + 1];
	}

	/* Find the least squares result */
	dgelsy_driver(A, b, x, num_eqs, num_vars, 1);

	error = 0.0;
	for (i = 0; i < num_points; i++) {
		double dx, dy;
		int Roff = 9 * i;
		int toff = 3 * i;
		double pp[3];

		/* Compute projection error */
		matrix_product331(R + Roff, x, pp);
		pp[0] += t[toff + 0];
		pp[1] += t[toff + 1];
		pp[2] += t[toff + 2];

		dx = pp[0] / pp[2] - p[i].x;
		dy = pp[1] / pp[2] - p[i].y;
		error += dx * dx + dy * dy;
	}

	error = sqrt(error / num_points);

	// printf("[triangulate_n] Error [before polishing]: %0.3e\n", error);

	/* Run a non-linear optimization to refine the result */
	global_num_points = num_points;
	global_ps = p;
	global_Rs = R;  global_ts = t;
	lmdif_driver((void*)triangulate_n_residual, num_eqs, num_vars, x, 1.0e-5);

	error = 0.0;
	for (i = 0; i < num_points; i++) {
		double dx, dy;
		int Roff = 9 * i;
		int toff = 3 * i;
		double pp[3];

		/* Compute projection error */
		matrix_product331(R + Roff, x, pp);
		pp[0] += t[toff + 0];
		pp[1] += t[toff + 1];
		pp[2] += t[toff + 2];

		dx = pp[0] / pp[2] - p[i].x;
		dy = pp[1] / pp[2] - p[i].y;
		error += dx * dx + dy * dy;
	}

	error = sqrt(error / num_points);

	// printf("[triangulate_n] Error [after polishing]: %0.3e\n", error);

	if (error_out != NULL) {
		*error_out = error;
	}

	r = vec3d(x[0], x[1], x[2]);

	free(A);
	free(b);
	free(x);

	return r;
}


/* Find the point with the smallest squared projection error */
vec3d triangulate(vec2d p, vec2d q,
	double *R0, double *t0,
	double *R1, double *t1, double *error)
{
	double A[12];
	double b[4];
	double x[3];

	double dx1, dx2, dy1, dy2;

	A[0] = R0[0] - p.x * R0[6];
	A[1] = R0[1] - p.x * R0[7];
	A[2] = R0[2] - p.x * R0[8];

	A[3] = R0[3] - p.y * R0[6];
	A[4] = R0[4] - p.y * R0[7];
	A[5] = R0[5] - p.y * R0[8];

	A[6] = R1[0] - q.x * R1[6];
	A[7] = R1[1] - q.x * R1[7];
	A[8] = R1[2] - q.x * R1[8];

	A[9] = R1[3]  - q.y * R1[6];
	A[10] = R1[4] - q.y * R1[7];
	A[11] = R1[5] - q.y * R1[8];

	b[0] = t0[2] * p.x - t0[0];
	b[1] = t0[2] * p.y - t0[1];
	b[2] = t1[2] * q.x - t1[0];
	b[3] = t1[2] * q.y - t1[1];

	/* Find the least squares result */
	dgelsy_driver(A, b, x, 4, 3, 1);

	/* Run a non-linear optimization to refine the result */
	global_p = p;
	global_q = q;
	global_R0 = R0;  global_t0 = t0;
	global_R1 = R1;  global_t1 = t1;
	lmdif_driver((void*)triangulation_residual, 4, 3, x, 1.0e-10);

	if (error != NULL) {
		double pp[3], qp[3];

		/* Compute projection error */
		matrix_product331(R0, x, pp);
		pp[0] += t0[0];
		pp[1] += t0[1];
		pp[2] += t0[2];

		matrix_product331(R1, x, qp);
		qp[0] += t1[0];
		qp[1] += t1[1];
		qp[2] += t1[2];

		dx1 = pp[0] / pp[2] - p.x;
		dy1 = pp[1] / pp[2] - p.y;

		dx2 = qp[0] / qp[2] - q.x;
		dy2 = qp[1] / qp[2] - q.y;

		*error = dx1 * dx1 + dy1 * dy1 + dx2 * dx2 + dy2 * dy2;
	}

	return vec3d(x[0], x[1], x[2]);
}

/* Given an F matrix, two calibration matrices, and a point correspondence, find R and t */
void find_extrinsics(double *F, double *K1, double *K2,
	vec2d p1, vec2d p2, double *R, double *t) {
	double E[9];
	double K1_inv[9], K2_inv[9], tmp[9];

	/* Find the essential matrix */
	matrix_invert(3, K1, K1_inv);
	matrix_invert(3, K2, K2_inv);

	matrix_product33(F, K1, tmp);
	matrix_transpose_product(3, 3, 3, 3, K2, tmp, E);

	find_extrinsics_essential(E, p1, p2, R, t);
}

/* Given an E matrix and a point correspondence, find R and t */
int find_extrinsics_essential(double *E, vec2d p1, vec2d p2,
	double *R, double *t)
{
	double tmp[9], tmp2[3], Qv[3];
	double U[9], S[3], VT[9];
	double tu[3], Ra[9], Rb[9];

	double D[9] =
	{ 0.0, 1.0, 0.0,
	-1.0, 0.0, 0.0,
	0.0, 0.0, 1.0 };

	double DT[9] =
	{ 0.0, -1.0, 0.0,
	1.0, 0.0, 0.0,
	0.0, 0.0, 1.0 };

	double I[9] =
	{ 1.0, 0.0, 0.0,
	0.0, 1.0, 0.0,
	0.0, 0.0, 1.0 };

	double t0[3] = { 0.0, 0.0, 0.0 };

	vec3d Q, PQ;
	double c1, c2;
	double error;

	/* Now find the SVD of E */
	dgesvd_driver(3, 3, E, U, S, VT);

	/* Now find R and t */
	tu[0] = U[2];  tu[1] = U[5];  tu[2] = U[8];

	matrix_product33(U, D, tmp);
	matrix_product33(tmp, VT, Ra);
	matrix_product33(U, DT, tmp);

	matrix_product33(tmp, VT, Rb);

	if (matrix_determinant3(Ra) < 0.0) {
		// printf("flipping...\n");
		matrix_scale(3, 3, Ra, -1.0, Ra);
	}

	if (matrix_determinant3(Rb) < 0.0) {
		// printf("flopping...\n");
		matrix_scale(3, 3, Rb, -1.0, Rb);
	}

	/* Figure out which configuration is correct using the supplied
	* point */

	Q = triangulate(p1, p2, I, t0, Ra, tu, &error);
	Qv[0] = Q.x, Qv[1] = Q.y, Qv[2] = Q.z;
	matrix_product331(Ra, Qv, tmp);
	matrix_sum(3, 1, 3, 1, tmp, tu, tmp2);
	PQ = vec3d(tmp2[0], tmp2[1], tmp2[2]);

	c1 = Q.z;
	c2 = PQ.z;

	// EDIT!!!
	if (c1 < 0 && c2 < 0) {
		memcpy(R, Ra, 9 * sizeof(double));
		t[0] = tu[0]; t[1] = tu[1]; t[2] = tu[2];
		// printf("[find_extrinsics] Case 1\n");
	}
	else if (c1 > 0 && c2 > 0) { // EDIT!!!
		memcpy(R, Ra, 9 * sizeof(double));
		// // t[0] = -tu[0]; t[1] = -tu[1]; t[2] = -tu[2];
		t[0] = -tu[0]; t[1] = -tu[1]; t[2] = -tu[2];
		// printf("[find_extrinsics] Case 2\n");

		Q = triangulate(p1, p2, I, t0, R, t, &error);
	}
	else {
		/* Triangulate again */
		Q = triangulate(p1, p2, I, t0, Rb, tu, &error);
		Qv[0] = Q.x, Qv[1] = Q.y, Qv[2] = Q.z;
		matrix_product331(Rb, Qv, tmp);
		matrix_sum(3, 1, 3, 1, tmp, tu, tmp2);
		PQ = vec3d(tmp2[0], tmp2[1], tmp2[2]);

		c1 = Q.z;
		c2 = PQ.z;

		// EDIT!!!
		if (c1 < 0 && c2 < 0) {
			memcpy(R, Rb, 9 * sizeof(double));
			t[0] = tu[0]; t[1] = tu[1]; t[2] = tu[2];
			// printf("[find_extrinsics] Case 3\n");
		}
		else if (c1 > 0 && c2 > 0) { // EDIT!!!
			memcpy(R, Rb, 9 * sizeof(double));
			// // t[0] = -tu[0]; t[1] = -tu[1]; t[2] = -tu[2];
			t[0] = -tu[0]; t[1] = -tu[1]; t[2] = -tu[2];
			// printf("[find_extrinsics] Case 4\n");

			Q = triangulate(p1, p2, I, t0, R, t, &error);
		}
		else {
			printf("[find_extrinsics] Error: no case found!\n");
			return 0;
		}
	}

	matrix_product331(R, Q.data(), tmp2);
	matrix_sum(3, 1, 3, 1, tmp2, t, tmp2);

	printf("[find_extrinsics] error: %0.3e\n", error);

	return 1;
}

/* Given an E matrix and a point correspondence, find R and t */
int find_extrinsics_essential_multipt(double *E, int n,
	vec2d *p1, vec2d *p2,
	double *R, double *t)
{
	double tmp[9], tmp2[3], Qv[3];
	double U[9], S[3], VT[9];
	double tu[3], Ra[9], Rb[9];

	double D[9] =
	{ 0.0, 1.0, 0.0,
	-1.0, 0.0, 0.0,
	0.0, 0.0, 1.0 };

	double DT[9] =
	{ 0.0, -1.0, 0.0,
	1.0, 0.0, 0.0,
	0.0, 0.0, 1.0 };

	double I[9] =
	{ 1.0, 0.0, 0.0,
	0.0, 1.0, 0.0,
	0.0, 0.0, 1.0 };

	double t0[3] = { 0.0, 0.0, 0.0 };

	vec3d Q, PQ;
	double c1, c2;
	double error;

	int i;
	int c1_pos = 0, c1_neg = 0;
	int c2_pos = 0, c2_neg = 0;

	/* Now find the SVD of E */
	dgesvd_driver(3, 3, E, U, S, VT);

	/* Now find R and t */
	tu[0] = U[2];  tu[1] = U[5];  tu[2] = U[8];

	matrix_product33(U, D, tmp);
	matrix_product33(tmp, VT, Ra);
	matrix_product33(U, DT, tmp);

	matrix_product33(tmp, VT, Rb);

	if (matrix_determinant3(Ra) < 0.0) {
		// printf("flipping...\n");
		matrix_scale(3, 3, Ra, -1.0, Ra);
	}

	if (matrix_determinant3(Rb) < 0.0) {
		// printf("flopping...\n");
		matrix_scale(3, 3, Rb, -1.0, Rb);
	}

	/* Figure out which configuration is correct using the supplied
	* points */

	for (i = 0; i < n; i++) {
		Q = triangulate(p1[i], p2[i], I, t0, Ra, tu, &error);
		Qv[0] = Q.x, Qv[1] = Q.y, Qv[2] = Q.z;
		matrix_product331(Ra, Qv, tmp);
		matrix_sum(3, 1, 3, 1, tmp, tu, tmp2);
		PQ = vec3d(tmp2[0], tmp2[1], tmp2[2]);

		c1 = Q.z;
		c2 = PQ.z;

		if (c1 > 0)
			c1_pos++;
		else
			c1_neg++;

		if (c2 > 0)
			c2_pos++;
		else
			c2_neg++;
	}

	// printf("[1] c1_pos: %d, c1_neg: %d\n", c1_pos, c1_neg);
	// printf("[1] c2_pos: %d, c2_neg: %d\n", c2_pos, c2_neg);

	// EDIT!!!
	if (c1_pos < c1_neg && c2_pos < c2_neg) {
		memcpy(R, Ra, 9 * sizeof(double));
		t[0] = tu[0]; t[1] = tu[1]; t[2] = tu[2];
		// printf("[find_extrinsics] Case 1\n");
	}
	else if (c1_pos > c1_neg && c2_pos > c2_neg) { // EDIT!!!
		memcpy(R, Ra, 9 * sizeof(double));
		// // t[0] = -tu[0]; t[1] = -tu[1]; t[2] = -tu[2];
		t[0] = -tu[0]; t[1] = -tu[1]; t[2] = -tu[2];
		// printf("[find_extrinsics] Case 2\n");
		// Q = triangulate(p1, p2, I, t0, R, t, &error);
	}
	else {
		/* Triangulate again */
		c1_pos = c1_neg = c2_pos = c2_neg = 0;

		for (i = 0; i < n; i++) {
			Q = triangulate(p1[i], p2[i], I, t0, Rb, tu, &error);
			Qv[0] = Q.x, Qv[1] = Q.y, Qv[2] = Q.z;
			matrix_product331(Rb, Qv, tmp);
			matrix_sum(3, 1, 3, 1, tmp, tu, tmp2);
			PQ = vec3d(tmp2[0], tmp2[1], tmp2[2]);

			c1 = Q.z;
			c2 = PQ.z;

			if (c1 > 0)
				c1_pos++;
			else
				c1_neg++;

			if (c2 > 0)
				c2_pos++;
			else
				c2_neg++;
		}

		// printf("[2] c1_pos: %d, c1_neg: %d\n", c1_pos, c1_neg);
		// printf("[2] c2_pos: %d, c2_neg: %d\n", c2_pos, c2_neg);

		if (c1_pos < c1_neg && c2_pos < c2_neg) { // EDIT!!!
			memcpy(R, Rb, 9 * sizeof(double));
			t[0] = tu[0]; t[1] = tu[1]; t[2] = tu[2];
			// printf("[find_extrinsics] Case 3\n");
		}
		else if (c1_pos > c1_neg && c2_pos > c2_neg) { // EDIT!!!
			memcpy(R, Rb, 9 * sizeof(double));
			// // t[0] = -tu[0]; t[1] = -tu[1]; t[2] = -tu[2];
			t[0] = -tu[0]; t[1] = -tu[1]; t[2] = -tu[2];
			// printf("[find_extrinsics] Case 4\n");

			// Q = triangulate(p1, p2, I, t0, R, t, &error);
		}
		else {
			fprintf(stderr, "[find_extrinsics] Error: no case found!\n");
			return 0;
		}
	}

	return 1;
}


static vec3d *condition_points_3D(int num_points, vec3d *pts, double *T) {
	vec3d *pts_new = (vec3d *)malloc(sizeof(vec3d) * num_points);
	vec3d *pts_zero_mean = (vec3d *)malloc(sizeof(vec3d) * num_points);
	double U[9], S[3], VT[9], U16[9], Sinv[16], Ttrans[16], action[16];

	vec3d mean = Geom::vec_mean(num_points, pts);

	double total_dist = 0.0;
	double avg_dist;
	double factor;
	int i;

	for (i = 0; i < num_points; i++) {
		double dx = pts[i].x - mean.x;
		double dy = pts[i].y - mean.y;
		double dz = pts[i].z - mean.z;

		pts_zero_mean[i] = vec3d(
			pts[i].x - mean.x,
			pts[i].y - mean.y,
			pts[i].z - mean.z);

		total_dist += sqrt(dx * dx + dy * dy + dz * dz);
	}

	vec_svd(num_points, pts_zero_mean, U, S, VT);

	Sinv[0] = sqrt(3.0) / sqrt(S[0]);  Sinv[1] = 0.0;  Sinv[2] = 0.0;  Sinv[3] = 0.0;
	Sinv[4] = 0.0;  Sinv[5] = sqrt(3.0) / sqrt(S[1]);  Sinv[6] = 0.0;  Sinv[7] = 0.0;
	Sinv[8] = 0.0;  Sinv[9] = 0.0;  Sinv[10] = sqrt(3.0) / sqrt(S[2]); Sinv[11] = 0.0;
	Sinv[12] = 0.0; Sinv[13] = 0.0; Sinv[14] = 0.0; Sinv[15] = 1.0;

	matrix_ident(4, U16);
	memcpy(U16 + 0, U + 0, 3 * sizeof(double));
	memcpy(U16 + 4, U + 3, 3 * sizeof(double));
	memcpy(U16 + 8, U + 6, 3 * sizeof(double));

	matrix_ident(4, Ttrans);
	Ttrans[3]  = -mean.x;
	Ttrans[7]  = -mean.y;
	Ttrans[11] = -mean.z;

	avg_dist = total_dist / num_points;
	factor = sqrt(3.0) / avg_dist;

	//matrix_transpose_product2(4, 4, 4, 4, Sinv, U16, action);
	matrix_transpose_product(4, 4, 4, 4, Sinv, U16, action);
	matrix_product44(action, Ttrans, T);

	for (i = 0; i < num_points; i++) {
		double x = factor * (pts[i].x - mean.x);
		double y = factor * (pts[i].y - mean.y);
		double z = factor * (pts[i].z - mean.z);

		double pt[4] = { pts[i].x, pts[i].y, pts[i].z, 1.0 };
		double Tpt[4];

		matrix_product441(T, pt, Tpt);

		pts_new[i] = vec3d(Tpt[0], Tpt[1], Tpt[2]);
	}

	free(pts_zero_mean);

	return pts_new;
}

static vec2d *condition_points_2D(int num_points, vec2d *pts, double *T) {
	vec2d *pts_new = (vec2d *)malloc(sizeof(vec2d) * num_points);

	vec2d mean = Geom::vec_mean(num_points, pts);
	double total_dist = 0.0;
	double avg_dist;
	double factor;
	int i;

	for (i = 0; i < num_points; i++) {
		double dx = pts[i].x - mean.x;
		double dy = pts[i].y - mean.y;
		total_dist += sqrt(dx * dx + dy * dy);
	}

	avg_dist = total_dist / num_points;
	factor = sqrt(2.0) / avg_dist;

	for (i = 0; i < num_points; i++) {
		double x = factor * (pts[i].x - mean.x);
		double y = factor * (pts[i].y - mean.y);
		pts_new[i] = vec2d(x, y);
	}

	T[0] = factor;  T[1] = 0.0;     T[2] = -factor * mean.x;
	T[3] = 0.0;     T[4] = factor;  T[5] = -factor * mean.y;
	T[6] = 0.0;      T[7] = 0.0;    T[8] = 1.0;

	return pts_new;
}

/* Solve for a 3x4 projection matrix, given a set of 3D points and 2D
* projections */
int find_projection_3x4(int num_pts, vec3d *points, vec2d *projs, double *P) {
	if (num_pts < 6) {
		printf("[find_projection_3x4] Need at least 6 points!\n");
		return -1;
	}
	else {

		// #define _CONDITION_
#ifdef _CONDITION_
		double Tpoints[16];
		vec3d *points_new = condition_points_3D(num_pts, points, Tpoints);

		double Tprojs[9];
		vec2d *projs_new = condition_points_2D(num_pts, projs, Tprojs);

		double Tprojs_inv[9];
		double Ptmp[12];
#else
		vec3d *points_new = points;
		vec2d *projs_new = projs;
#endif

		int num_eqns = 2 * num_pts;
		int num_vars = 11;

		double *A = new double[num_eqns * num_vars];
		double *b = new double[num_eqns];
		double X[11];

		double error = 0.0;

		int i;

		for (i = 0; i < num_pts; i++) {
			double *row1 = A + 2 * i * num_vars;
			double *row2 = A + (2 * i + 1) * num_vars;

			row1[0] = points_new[i].x;
			row1[1] = points_new[i].y;
			row1[2] = points_new[i].z;
			row1[3] = 1.0;

			row1[4] = 0.0;
			row1[5] = 0.0;
			row1[6] = 0.0;
			row1[7] = 0.0;

			// EDIT!!!
			row1[8] = projs_new[i].x * points_new[i].x;
			row1[9] = projs_new[i].x * points_new[i].y;
			row1[10] = projs_new[i].x * points_new[i].z;

			b[2 * i] = -projs_new[i].x;


			row2[0] = 0.0;
			row2[1] = 0.0;
			row2[2] = 0.0;
			row2[3] = 0.0;

			row2[4] = points_new[i].x;
			row2[5] = points_new[i].y;
			row2[6] = points_new[i].z;
			row2[7] = 1.0;

			// EDIT!!!
			row2[8]  = projs_new[i].y * points_new[i].x;
			row2[9]  = projs_new[i].y * points_new[i].y;
			row2[10] = projs_new[i].y * points_new[i].z;

			b[2 * i + 1] = -projs_new[i].y;
		}

		dgelsy_driver(A, b, X, num_eqns, num_vars, 1);

		memcpy(P, X, sizeof(double) * 11);
		P[11] = 1.0;

#ifdef _CONDITION_
		matrix_invert(3, Tprojs, Tprojs_inv);
		matrix_product(3, 3, 3, 4, Tprojs_inv, P, Ptmp);
		matrix_product(3, 4, 4, 4, Ptmp, Tpoints, P);

		matrix_scale(3, 4, P, 1.0 / P[11], P);
#endif

		for (i = 0; i < num_pts; i++) {
			double pt[4] = { 
				points[i].x,
				points[i].y,
				points[i].z, 1.0 };
			double pr[3];
			double dx, dy, dist;

			matrix_product341(P, pt, pr);
			// EDIT!!!
			pr[0] /= -pr[2];
			pr[1] /= -pr[2];

			dx = pr[0] - projs[i].x;
			dy = pr[1] - projs[i].y;

			dist = dx * dx + dy * dy;

			error += dist;
		}

		// printf("[find_projection_3x4] Average error is %0.3f\n", 
		//       error / num_pts);

		delete[] (A);
		delete[] (b);

#ifdef _CONDITION_
		free(points_new);
		free(projs_new);
#endif

		return 0;
	}
}

static int global_num_pts;
static vec3d *global_points;
static vec2d *global_projs;

static void projection_residual(const int *m, const int *n, double *x,
	double *fvec, double *iflag)
{
	int i;

	double P[12];
	memcpy(P, x, sizeof(double) * 11);
	P[11] = 1.0;

	for (i = 0; i < global_num_pts; i++) {
		double pt[4] = { 
			global_points[i].x,
			global_points[i].y,
			global_points[i].z, 1.0 };

		double pr[3];
		double dx, dy;

		matrix_product341(P, pt, pr);
		// EDIT!!
		pr[0] /= -pr[2];
		pr[1] /= -pr[2];

		dx = pr[0] - global_projs[i].x;
		dy = pr[1] - global_projs[i].y;

		fvec[2 * i + 0] = dx;
		fvec[2 * i + 1] = dy;
	}
}

/* Solve for a 3x4 projection matrix, given a set of 3D points and 2D
* projections using non-linear optimization */
int find_projection_3x4_nonlinear(int num_pts, vec3d *points, vec2d *projs,
	double *Pin, double *Pout)
{
	if (num_pts < 6) {
		printf("[find_projection_3x4_nonlinear] Need at least 6 points!\n");
		return -1;
	}
	else {
		int num_eqns = 2 * num_pts;
		int num_vars = 11;
		double x[11];

		global_num_pts = num_pts;
		global_points = points;
		global_projs = projs;

		memcpy(x, Pin, sizeof(double) * 11);
		lmdif_driver((void*)projection_residual, num_eqns, num_vars, x, 1.0e-5);

		memcpy(Pout, x, sizeof(double) * 11);
		Pout[11] = 1.0;

		return 0;
	}
}


/* Solve for a 3x4 projection matrix using RANSAC, given a set of 3D
* points and 2D projections */
int find_projection_3x4_ransac(int num_pts, vec3d *points, vec2d *projs,
	double *P,
	int ransac_rounds, double ransac_threshold)
{
	if (num_pts < 6) {
		printf("[find_projection_3x4_ransac] Error: need at least 6 points!\n");
		return -1;
	}
	else {
#define MIN_PTS 6
		// const int min_pts = 6;
		int *inliers = (int *)malloc(sizeof(int) * num_pts);
		int indices[MIN_PTS];
		int round, i, j;
		int max_inliers = 0;
		double max_error = 0.0;
		double Pbest[12];
		int num_inliers = 0, num_inliers_new = 0;
		vec3d *pts_final = NULL;
		vec2d *projs_final = NULL;
		double Plinear[12];

		double Rinit[9];
		double triangular[9], orthogonal[9];
		int neg, sign;

		double thresh_sq = ransac_threshold * ransac_threshold;
		double error = 0.0;

		int num_inliers_polished = 0;

		for (round = 0; round < ransac_rounds; round++) {
			vec3d pts_inner[MIN_PTS];
			vec2d projs_inner[MIN_PTS];
			double Ptmp[12];

			num_inliers = 0;
			for (i = 0; i < MIN_PTS; i++) {
				int redo = 0;
				int idx;
				int redo_count = 0;

				do {
					if (redo_count > 10000) {
						free(inliers);
						return -1;
					}

					idx = rand() % num_pts;

					redo = 0;
					for (j = 0; j < i; j++) {
						if (idx == indices[j]) {
							redo = 1;
							break;
						}
						else if (projs[idx].x == projs[indices[j]].x &&
							projs[idx].y == projs[indices[j]].y) {
							redo = 1;
						}
					}

					redo_count++;
				} while (redo);

				indices[i] = idx;
				pts_inner[i] = points[idx];
				projs_inner[i] = projs[idx];
			}

			/* Solve for the parameters */
			find_projection_3x4(MIN_PTS, pts_inner, projs_inner, Ptmp);

#if 1
			/* Fix the sign on the P matrix */
			memcpy(Rinit + 0, Ptmp + 0, 3 * sizeof(double));
			memcpy(Rinit + 3, Ptmp + 4, 3 * sizeof(double));
			memcpy(Rinit + 6, Ptmp + 8, 3 * sizeof(double));

			dgerqf_driver(3, 3, Rinit, triangular, orthogonal);

			/* Check the parity along the diagonal */
			neg =
				(triangular[0] < 0.0) +
				(triangular[4] < 0.0) +
				(triangular[8] < 0.0);

			if ((neg % 2) == 1) {
				sign = -1;
			}
			else {
				sign = 1;
			}
#endif

			/* Count the number of inliers */
			error = 0.0;
			for (i = 0; i < num_pts; i++) {
				double pt[4] = { 
					points[i].x,
					points[i].y,
					points[i].z, 1.0 };
				double pr[3];
				double dx, dy, dist;

				matrix_product341(Ptmp, pt, pr);

				/* Check cheirality */
				// EDIT!!!
				if (sign * pr[2] > 0.0)
					continue;

				// EDIT!!!
				pr[0] /= -pr[2];
				pr[1] /= -pr[2];

				dx = pr[0] - projs[i].x;
				dy = pr[1] - projs[i].y;

				dist = dx * dx + dy * dy;

				if (dist < thresh_sq) {
					inliers[num_inliers] = i;
					num_inliers++;
					error += dist;
				}
			}

			if (num_inliers > max_inliers) {
				memcpy(Pbest, Ptmp, sizeof(double) * 12);
				max_error = error;
				max_inliers = num_inliers;
			}
		}

		memcpy(P, Pbest, sizeof(double) * 12);

		printf("[find_projection_3x4_ransac] num_inliers = %d (out of %d)\n",
			max_inliers, num_pts);
		printf("[find_projection_3x4_ransac] error = %0.3f\n",
			sqrt(max_error / max_inliers));

		if (max_inliers < 6) {
			printf("[find_projection_3x4_ransac] "
				"Too few inliers to continue.\n");

			free(inliers);

			return -1;
		}

		/* Do the final least squares minimization */

#if 1
		/* Fix the sign on the P matrix */
		memcpy(Rinit + 0, Pbest + 0, 3 * sizeof(double));
		memcpy(Rinit + 3, Pbest + 4, 3 * sizeof(double));
		memcpy(Rinit + 6, Pbest + 8, 3 * sizeof(double));

		dgerqf_driver(3, 3, Rinit, triangular, orthogonal);

		/* Check the parity along the diagonal */
		neg =
			(triangular[0] < 0.0) +
			(triangular[4] < 0.0) +
			(triangular[8] < 0.0);

		if ((neg % 2) == 1) {
			sign = -1;
		}
		else {
			sign = 1;
		}
#endif

		num_inliers = 0;
		pts_final = (vec3d *)malloc(sizeof(vec3d) * max_inliers);
		projs_final = (vec2d *)malloc(sizeof(vec2d) * max_inliers);

		for (i = 0; i < num_pts; i++) {
			double pt[4] = { 
				points[i].x,
				points[i].y,
				points[i].z, 1.0 };

			double pr[3];
			double dx, dy, dist;

			matrix_product341(Pbest, pt, pr);

			/* Check cheirality */
			// EDIT!!!
			if (sign * pr[2] > 0.0)
				continue;

			// EDIT!!!
			pr[0] /= -pr[2];
			pr[1] /= -pr[2];

			dx = pr[0] - projs[i].x;
			dy = pr[1] - projs[i].y;

			dist = dx * dx + dy * dy;

			if (dist < thresh_sq) {
				pts_final[num_inliers] = points[i];
				projs_final[num_inliers] = projs[i];
				num_inliers++;
			}
		}

		if (num_inliers != max_inliers) {
			printf("[find_projection_3x4_ransac] Error! There was a miscount "
				"somewhere: (%d != %d)\n", num_inliers, max_inliers);
		}

		find_projection_3x4(max_inliers, pts_final, projs_final, Plinear);

#if 1
		/* Fix the sign on the P matrix */
		memcpy(Rinit + 0, Plinear + 0, 3 * sizeof(double));
		memcpy(Rinit + 3, Plinear + 4, 3 * sizeof(double));
		memcpy(Rinit + 6, Plinear + 8, 3 * sizeof(double));

		dgerqf_driver(3, 3, Rinit, triangular, orthogonal);

		/* Check the parity along the diagonal */
		neg =
			(triangular[0] < 0.0) +
			(triangular[4] < 0.0) +
			(triangular[8] < 0.0);

		if ((neg % 2) == 1) {
			sign = -1;
		}
		else {
			sign = 1;
		}
#endif

		for (i = 0; i < num_pts; i++) {
			double pt[4] =
			{ points[i].x, points[i].y, points[i].z, 1.0 };
			double pr[3];
			double dx, dy, dist;

			matrix_product341(Plinear, pt, pr);

			// EDIT!!!
			if (sign * pr[2] > 0.0)
				continue;

			// EDIT!!!
			pr[0] /= -pr[2];
			pr[1] /= -pr[2];

			dx = pr[0] - projs[i].x;
			dy = pr[1] - projs[i].y;

			dist = dx * dx + dy * dy;

			if (dist < thresh_sq) {
				num_inliers_new++;
			}
		}

		if (num_inliers_new < max_inliers) {
			printf("[find_projection_3x4_ransac] Reverting to old solution\n");
			memcpy(Plinear, Pbest, 12 * sizeof(double));
		}

// 		printf("Best matrix (pre-opt):\n");
// 		matrix_print(3, 4, Plinear);

		error = 0.0;
		for (i = 0; i < max_inliers; i++) {
			double pt[4] =
			{ pts_final[i].x, pts_final[i].y, pts_final[i].z, 1.0 };
			double pr[3];
			double dx, dy, dist;

			matrix_product341(Plinear, pt, pr);
			pr[0] /= pr[2];
			pr[1] /= pr[2];

			dx = pr[0] - projs_final[i].x;
			dy = pr[1] - projs_final[i].y;

			dist = dx * dx + dy * dy;

			error += dist;
		}

// 		printf("Old error: %0.3e\n", sqrt(error / max_inliers));

		/* Polish the result */
		if (max_inliers >= 6) {
			int num_inliers_polished = 0;
			find_projection_3x4_nonlinear(max_inliers, pts_final, projs_final,
				Plinear, P);

#if 1
			/* Fix the sign on the P matrix */
			memcpy(Rinit + 0, P + 0, 3 * sizeof(double));
			memcpy(Rinit + 3, P + 4, 3 * sizeof(double));
			memcpy(Rinit + 6, P + 8, 3 * sizeof(double));

			dgerqf_driver(3, 3, Rinit, triangular, orthogonal);

			/* Check the parity along the diagonal */
			neg =
				(triangular[0] < 0.0) +
				(triangular[4] < 0.0) +
				(triangular[8] < 0.0);

			if ((neg % 2) == 1) {
				sign = -1;
			}
			else {
				sign = 1;
			}
#endif

			/* Check that the number of inliers hasn't gone down */
			num_inliers_polished = 0;
			for (i = 0; i < num_pts; i++) {
				double pt[4] =
				{ points[i].x, points[i].y, points[i].z, 1.0 };
				double pr[3];
				double dx, dy, dist;

				matrix_product341(P, pt, pr);

				// EDIT!!!
				if (sign * pr[2] > 0.0)
					continue;

				// EDIT!!!
				pr[0] /= -pr[2];
				pr[1] /= -pr[2];

				dx = pr[0] - projs[i].x;
				dy = pr[1] - projs[i].y;

				dist = dx * dx + dy * dy;

				if (dist < thresh_sq) {
					num_inliers_polished++;
				}
			}

			if (num_inliers_polished < max_inliers) {
				printf("Decreased number of inliers (%d < %d), reverting\n",
					num_inliers_polished, max_inliers);

				memcpy(P, Plinear, sizeof(double) * 12);
			}
		}
		else {
			memcpy(P, Plinear, sizeof(double) * 12);
		}

// 		printf("Best matrix (post-opt):\n");
// 		matrix_print(3, 4, P);

		error = 0.0;
		for (i = 0; i < max_inliers; i++) {
			double pt[4] =
			{ pts_final[i].x, pts_final[i].y, pts_final[i].z, 1.0 };
			double pr[3];
			double dx, dy, dist;

			matrix_product341(P, pt, pr);

			// EDIT!!!
			pr[0] /= -pr[2];
			pr[1] /= -pr[2];

			dx = pr[0] - projs_final[i].x;
			dy = pr[1] - projs_final[i].y;

			dist = dx * dx + dy * dy;

			error += dist;
		}

		printf("New error: %0.3e\n", sqrt(error / max_inliers));

		free(inliers);
		free(pts_final);
		free(projs_final);

		return max_inliers;
	}
#undef MIN_PTS
}

