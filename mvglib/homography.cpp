
#include "homography.h"
#include "../math/matrix_driver.h"



static vec3d* condition_points(int num_points, vec3d* pts, double *T) {
	vec3d* pts_new = (vec3d*)malloc(sizeof(vec3d) * num_points);

	// vec3d mean = v3_mean(num_points, pts);
	vec3d mean(0, 0, 0);
	for (int i = 0; i < num_points; ++i) {
		mean = mean + pts[i];
	}
	mean /= num_points;

	double total_dist = 0.0;
	double avg_dist;
	double factor;
	int i;

	for (i = 0; i < num_points; i++) {
		double dx = pts[i].x - mean.x;
		double dy = pts[i].y - mean.y;
		total_dist += std::sqrt(dx * dx + dy * dy);
	}

	avg_dist = total_dist / num_points;
	factor = std::sqrt(2.0) / avg_dist;

	for (i = 0; i < num_points; i++) {
		double x = factor * (pts[i].x - mean.x);
		double y = factor * (pts[i].y - mean.y);
		pts_new[i] = vec3d(x, y, 1.0);
	}

	T[0] = factor;  T[1] = 0.0;     T[2] = -factor * mean.x;
	T[3] = 0.0;     T[4] = factor;  T[5] = -factor * mean.y;
	T[6] = 0.0;      T[7] = 0.0;    T[8] = 1.0;

	return pts_new;
}

/* Computes the homography that, when applied to the points in l_pts,
* minimizes the least-squares error between the result and the
* corresponding points in r_pts.
*
* n -- number of points
* r_pts -- matches
* l_pts -- initial points
* Tout -- on return, contains the 3x3 transformation matrix */
void align_homography(int num_pts, vec3d* r_pts, vec3d* l_pts,
	double *Tout, int refine)
{
	int m = num_pts * 2;
	int n = 8;
	int nrhs = 1;
	int i, base;

	double *A = new double[m * n];    /* Left-hand matrix */
	double *B = new double[m * nrhs]; /* Right-hand matrix */

	double Ttmp[9];
	double T1[9], T2[9];

#define _CONDITION_
#ifdef _CONDITION_
	/* Normalize the points */
	vec3d* r_pts_norm = condition_points(num_pts, r_pts, T1);
	vec3d* l_pts_norm = condition_points(num_pts, l_pts, T2);
	double T1inv[9];
#else
	vec3d* r_pts_norm = r_pts;
	vec3d* l_pts_norm = l_pts;
#endif

	for (i = 0; i < num_pts; i++) {
		base = 2 * i * n;
		A[base + 0] = l_pts_norm[i].x;
		A[base + 1] = l_pts_norm[i].y;
		A[base + 2] = 1.0;
		A[base + 3] = A[base + 4] = A[base + 5] = 0.0;
		A[base + 6] = -l_pts_norm[i].x * r_pts_norm[i].x;
		A[base + 7] = -l_pts_norm[i].y * r_pts_norm[i].x;

		base = (2 * i + 1) * n;
		A[base + 0] = A[base + 1] = A[base + 2] = 0.0;
		A[base + 3] = l_pts_norm[i].x;
		A[base + 4] = l_pts_norm[i].y;
		A[base + 5] = 1.0;
		A[base + 6] = -l_pts_norm[i].x * r_pts_norm[i].y;
		A[base + 7] = -l_pts_norm[i].y * r_pts_norm[i].y;

		B[2 * i + 0] = r_pts_norm[i].x;
		B[2 * i + 1] = r_pts_norm[i].y;
	}

	/* Make the call to dgelsy */
	dgelsy_driver(A, B, Tout, m, n, nrhs);

	Tout[8] = 1.0;

#ifdef _CONDITION_
	/* Undo normalization */
	matrix_invert(3, T1, T1inv);

	matrix_product(3, 3, 3, 3, T1inv, Tout, Ttmp);
	matrix_product(3, 3, 3, 3, Ttmp, T2, Tout);

	matrix_scale(3, 3, Tout, 1.0 / Tout[8], Tout);
#endif

	if (refine) {
		memcpy(Ttmp, Tout, sizeof(double) * 9);
		align_homography_non_linear(num_pts, r_pts, l_pts, Ttmp, Tout);
	}

	free(A);
	free(B);

#ifdef _CONDITION_
	free(r_pts_norm);
	free(l_pts_norm);
#endif
}

static int global_num_pts;
static vec3d* global_r_pts;
static vec3d* global_l_pts;
static int global_round;

static void homography_resids(int *m, int *n, double *x, double *fvec, int *iflag)
{
	int i;

	double resids = 0.0;

	double H[9];
	memcpy(H, x, 8 * sizeof(double));
	H[8] = 1.0;

	if (*iflag == 0 && global_num_pts > 4) {
		printf("[Round %d]\n", global_round);
		printf("  H=(%0.5f, %0.5f, %0.5f, %0.5f, %0.5f, %0.5f, %0.5f, %0.5f, %0.1f)\n",
			H[0], H[1], H[2], H[3], H[4], H[5], H[6], H[7], H[8]);
		global_round++;
	}

	for (i = 0; i < global_num_pts; i++) {
		double p[3], q[3];

		p[0] = global_l_pts[i].x;
		p[1] = global_l_pts[i].y;
		p[2] = global_l_pts[i].z;

		if (*iflag == 0 && global_num_pts > 4)
			printf("    p=(%0.3f, %0.3f, %0.3f)\n", p[0], p[1], p[2]);

		matrix_product(3, 3, 3, 1, H, p, q);

		if (*iflag == 0 && global_num_pts > 4)
			printf("    q=(%0.3f, %0.3f, %0.3f)\n", q[0], q[1], q[2]);

		q[0] /= q[2];
		q[1] /= q[2];

		fvec[2 * i + 0] = q[0] - global_r_pts[i].x;
		fvec[2 * i + 1] = q[1] - global_r_pts[i].y;

		if (*iflag == 0 && global_num_pts > 4)
			printf("    (%0.3f, %0.3f) ==> (%0.3f, %0.3f)\n", q[0], q[1], global_r_pts[i].x, global_r_pts[i].y);
	}

	for (i = 0; i < 2 * global_num_pts; i++) {
		resids += fvec[i] * fvec[i];
	}

	if (*iflag == 0 && global_num_pts > 4)
		printf("resids = %0.3f\n", resids);
}

/* Use non-linear least squares to refine a homography */
void align_homography_non_linear(int num_pts, vec3d* r_pts, vec3d* l_pts,
	double *Tin, double *Tout)
{
	double x[8];

#if 0
	if (num_pts > 4) {
		printf("pre: ");
		matrix_print(3, 3, Tin);
	}
#endif

	memcpy(x, Tin, 8 * sizeof(double));

	global_num_pts = num_pts;
	global_r_pts = r_pts;
	global_l_pts = l_pts;
	global_round = 0;

	lmdif_driver((void*)homography_resids, 2 * num_pts, 8, x, 1.0e-4);

	memcpy(Tout, x, 8 * sizeof(double));
	Tout[8] = 1.0;

#if 0
	if (num_pts > 4) {
		printf("post: ");
		matrix_print(3, 3, Tout);
	}
#endif
}
