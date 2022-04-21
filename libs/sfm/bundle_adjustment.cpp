#include "bundle_adjustment.h"
#include "../math/matrix_driver.h"
#include "../../3rd_party/sba/sba.h"



#define TEST_FOCAL
#define NUM_STDDEV			2.0 // 3.0 // 6.0
#define NUM_ERROR_BINS			10
#define MIN_INLIERS_EST_PROJECTION	6 /* 7 */ /* 30 */ /* This constant needs adjustment */



namespace sfm {


	typedef struct {
		int num_cameras;               /* Number of cameras */
		int num_points;                /* Number of points */
		int num_params_per_camera;     /* Number of parameters for each camera */

		int est_focal_length;          /* Should the focal length be estimated? */
		int const_focal_length;        /* Is the focal length constant for all
							 * cameras? */
		int explicit_camera_centers;   /* Are the camera centers explicit? */
		int estimate_distortion;       /* Apply undistortion? */

		camera_params_t global_params;
		camera_params_t *init_params;  /* Initial camera parameters */

		vec3d *points;
	} sfm_global_t;




	static int global_num_points = 0;
	static sfm_global_t *global_params = NULL;
	static vec3d *global_points = NULL;
	static vec2d *global_projections = NULL;
	static int global_constrain_focal = 0;
	static double global_init_focal = 0.0;
	static double global_constrain_focal_weight = 0.0;
	static double global_constrain_rd_weight = 0.0;
	static int global_round = 0;



	static void sfm_project_point(int j, int i, double *aj, double *bi,
		double *xij, void *adata)
	{
		sfm_global_t *globs = (sfm_global_t *)adata;

		double K[9] =
		{ 1.0, 0.0, 0.0,
		0.0, 1.0, 0.0,
		0.0, 0.0, 1.0 };

		double *w, *dt;

		/* Compute intrinsics */
		if (!globs->est_focal_length) {
			K[0] = K[4] = globs->init_params[j].f; // globs->global_params.f;
		}
		else if (globs->const_focal_length) {
			printf("Error: case of constant focal length "
				"has not been implemented.\n");
			K[0] = K[4] = globs->global_params.f;
		}
		else {
			K[0] = K[4] = aj[6];
		}

		/* Compute translation, rotation update */
		dt = aj + 0;
		w = aj + 3;

#ifdef COLIN_HACK
		w[0] = w[1] = w[2] = 0.0;
		dt[2] = 0.0;
#endif

		if (globs->estimate_distortion == 0) {
			sfm_project(globs->init_params + j, K, w, dt, bi, xij,
				globs->explicit_camera_centers);
		}
		else {
			double Rnew[9];
			rot_update(globs->init_params->R, w, Rnew);
			sfm_project_rd(globs->init_params + j, K, aj + 7,
				Rnew, dt, bi, xij, 1, globs->explicit_camera_centers);
		}
	}


	void camera_refine_residual(const int *m, const int *n,
		double *x, double *fvec, int *iflag)
	{
		int i;
		double error = 0.0, error2 = 0.0;

		for (i = 0; i < global_num_points; i++) {
			double pt[3] = { 
				global_points[i].x,
				global_points[i].y,
				global_points[i].z };

			double proj[2], dx, dy;

			sfm_project_point(0, i, x, pt, proj, (void *)global_params);

			dx = global_projections[i].x - proj[0];
			dy = global_projections[i].y - proj[1];

			fvec[2 * i + 0] = dx;
			fvec[2 * i + 1] = dy;

			if (*iflag == 0) {
				error += dx * dx + dy * dy;
				error2 += sqrt(dx * dx + dy * dy);
			}
		}

		if (global_constrain_focal == 1) {
			double focal_diff = global_init_focal - x[6];
			fvec[2 * global_num_points] =
				global_constrain_focal_weight * focal_diff;

			if (global_params->estimate_distortion) {
				fvec[2 * global_num_points + 1] =
					-global_constrain_rd_weight * x[7];
				fvec[2 * global_num_points + 2] =
					-global_constrain_rd_weight * x[8];
			}
		}
		else if (global_params->estimate_distortion) {
			fvec[2 * global_num_points + 0] =
				-global_constrain_rd_weight * x[7];
			fvec[2 * global_num_points + 1] =
				-global_constrain_rd_weight * x[8];
		}

		if (*iflag == 0) {
			if (global_params->estimate_distortion) {
				// 			printf("  Round[%d]: RMS error = %0.8f [%0.8f], "
				// 				"f = %0.3f; %0.3e %0.3e\n",
				// 				global_round, sqrt(error / global_num_points),
				// 				error2 / global_num_points, x[6], x[7], x[8]);
			}
			else {
				if (global_params->est_focal_length) {
					// 				printf("  Round[%d]: RMS error = %0.8f [%0.8f], f = %0.3f\n",
					// 					global_round, sqrt(error / global_num_points),
					// 					error2 / global_num_points, x[6]);
				}
				else {
					// 				printf("  Round[%d]: RMS error = %0.8f [%0.8f]\n",
					// 					global_round, sqrt(error / global_num_points),
					// 					error2 / global_num_points);
				}
			}

			if (global_constrain_focal == 1) {
				printf("  Round[%d]: df = %0.3f\n",
					global_round, global_init_focal - x[6]);
			}

			global_round++;
		}
	}


	/* Refine the position of a single camera */
	void camera_refine(int num_points, vec3d *points, vec2d *projs,
		camera_params_t *params, int adjust_focal,
		int estimate_distortion)
	{
		if (adjust_focal) {
			int num_camera_params = 7;
			double x[9] = { params->t[0], params->t[1], params->t[2],
				0.0, 0.0, 0.0, params->f, params->k[0], params->k[1] };
			double Rnew[9];

			sfm_global_t globs;
			int focal_constraint = 0;

			if (estimate_distortion)
				num_camera_params += 2;

			globs.num_cameras = 1;
			globs.num_points = num_points;
			globs.est_focal_length = 1;
			globs.const_focal_length = 0;
			globs.explicit_camera_centers = 1;
			globs.global_params.f = params->f;
			globs.init_params = params;
			globs.estimate_distortion = estimate_distortion;

			global_num_points = num_points;
			global_params = &globs;
			global_points = points;
			global_projections = projs;
			global_round = 0;

			if (params->constrained[6]) {
				printf("[camera_refine] Constraining focal length to %0.3f "
					"(weight: %0.3f)\n",
					params->constraints[6],
					num_points * params->weights[6]);
				focal_constraint = 1;
				global_init_focal = params->constraints[6];
				global_constrain_focal = 1;
				global_constrain_focal_weight =
					1.0e0 /*1.0e1*/ * num_points * params->weights[6];
			}
			else {
				focal_constraint = 0;
				global_init_focal = 0.0;
				global_constrain_focal = 0;
				global_constrain_focal_weight = 0.0;
			}

			if (estimate_distortion) {
				global_constrain_rd_weight = 0.05 * num_points;
				// 1.0e-1 * num_points;
			}

			lmdif_driver2((void*)camera_refine_residual,
				2 * num_points + focal_constraint +
				2 * estimate_distortion,
				num_camera_params, x, 1.0e-12);

			/* Copy out the parameters */
			memcpy(params->t, x + 0, 3 * sizeof(double));
			rot_update(params->R, x + 3, Rnew);
			memcpy(params->R, Rnew, 9 * sizeof(double));
			params->f = x[6];

			if (estimate_distortion) {
				params->k[0] = x[7];
				params->k[1] = x[8];
			}
		}
		else {
			double x[6] = { params->t[0], params->t[1], params->t[2],
				0.0, 0.0, 0.0 };
			double Rnew[9];

			sfm_global_t globs;

			globs.num_cameras = 1;
			globs.num_points = num_points;
			globs.est_focal_length = 0;
			globs.const_focal_length = 1;
			globs.explicit_camera_centers = 1;
			globs.global_params.f = params->f;
			globs.init_params = params;
			globs.estimate_distortion = estimate_distortion;

			global_num_points = num_points;
			global_params = &globs;
			global_points = points;
			global_projections = projs;
			global_round = 0;

			global_init_focal = 0.0;
			global_constrain_focal = 0;
			global_constrain_focal_weight = 0.0;

			lmdif_driver2((void*)camera_refine_residual, 2 * num_points, 6, x, 1.0e-12);

			/* Copy out the parameters */
			memcpy(params->t, x + 0, 3 * sizeof(double));
			rot_update(params->R, x + 3, Rnew);
			memcpy(params->R, Rnew, 9 * sizeof(double));
		}
	}


	vec2d sfm_project_final(camera_params_t *params, vec3d pt,
		int explicit_camera_centers, int undistort) {
		double b_cam[3], b_proj[3];
		double b[3] = { pt.x, pt.y, pt.z };
		vec2d proj;

		/* Project! */
		if (!explicit_camera_centers) {
			matrix_product331(params->R, b, b_cam);
			b_cam[0] += params->t[0];
			b_cam[1] += params->t[1];
			b_cam[2] += params->t[2];
		}
		else {
			double b2[3];
			b2[0] = b[0] - params->t[0];
			b2[1] = b[1] - params->t[1];
			b2[2] = b[2] - params->t[2];
			matrix_product331(params->R, b2, b_cam);
		}

		if (!params->known_intrinsics) {
			b_proj[0] = b_cam[0] * params->f;
			b_proj[1] = b_cam[1] * params->f;
			b_proj[2] = b_cam[2];

			b_proj[0] /= -b_proj[2];
			b_proj[1] /= -b_proj[2];

			if (undistort) {
				double rsq =
					(b_proj[0] * b_proj[0] + b_proj[1] * b_proj[1]) /
					(params->f * params->f);
				double factor = 1.0 + params->k[0] * rsq + params->k[1] * rsq * rsq;

				b_proj[0] *= factor;
				b_proj[1] *= factor;
			}
		}
		else {
			/* Apply intrinsics */
			double x_n = -b_cam[0] / b_cam[2];
			double y_n = -b_cam[1] / b_cam[2];

			double *k = params->k_known;
			double rsq = x_n * x_n + y_n * y_n;
			double factor = 1.0 + k[0] * rsq +
				k[1] * rsq * rsq + k[4] * rsq * rsq * rsq;

			double dx_x = 2 * k[2] * x_n * y_n + k[3] * (rsq + 2 * x_n * x_n);
			double dx_y = k[2] * (rsq + 2 * y_n * y_n) + 2 * k[3] * x_n * y_n;

			double x_d = x_n * factor + dx_x;
			double y_d = y_n * factor + dx_y;

			double *K = params->K_known;
			b_proj[0] = K[0] * x_d + K[1] * y_d + K[2];
			b_proj[1] = K[4] * y_d + K[5];
		}

		proj = vec2d(b_proj[0], b_proj[1]);
		return proj;
	}



	/* Compute an updated rotation matrix given the initial rotation (R)
	* and the correction (w) */
	void rot_update(double *R, double *w, double *Rnew)
	{
		double theta, sinth, costh, n[3];
		double nx[9], nxsq[9];
		double term2[9], term3[9];
		double tmp[9], dR[9];

		double ident[9] =
		{ 1.0, 0.0, 0.0,
		0.0, 1.0, 0.0,
		0.0, 0.0, 1.0 };

		theta = sqrt(w[0] * w[0] + w[1] * w[1] + w[2] * w[2]);

		if (theta == 0.0) {
			memcpy(Rnew, R, sizeof(double) * 9);
			return;
		}

		n[0] = w[0] / theta;
		n[1] = w[1] / theta;
		n[2] = w[2] / theta;

		nx[0] = 0.0;   nx[1] = -n[2];  nx[2] = n[1];
		nx[3] = n[2];  nx[4] = 0.0;    nx[5] = -n[0];
		nx[6] = -n[1]; nx[7] = n[0];   nx[8] = 0.0;

		matrix_product33(nx, nx, nxsq);

		sinth = sin(theta);
		costh = cos(theta);

		matrix_scale(3, 3, nx, sinth, term2);
		matrix_scale(3, 3, nxsq, 1.0 - costh, term3);

		matrix_sum(3, 3, 3, 3, ident, term2, tmp);
		matrix_sum(3, 3, 3, 3, tmp, term3, dR);

		matrix_product33(dR, R, Rnew);
	}


	void sfm_project(camera_params_t *init, double *K,
		double *w, double *dt, double *b, double *p,
		int explicit_camera_centers)
	{
		double *R, *t;

		double Rnew[9];
		double tnew[3];

		double b_cam[3], b_proj[3];

		R = init->R;
		t = init->t;

		rot_update(R, w, Rnew);

		tnew[0] = dt[0]; // t[0] + dt[0];
		tnew[1] = dt[1]; // t[1] + dt[1];  // 0.0
		tnew[2] = dt[2]; // t[2] + dt[2];  // 0.0

		/* Project! */
		if (!explicit_camera_centers) {
			matrix_product331(Rnew, b, b_cam);
			b_cam[0] += tnew[0];
			b_cam[1] += tnew[1];
			b_cam[2] += tnew[2];
		}
		else {
			double b2[3];
			b2[0] = b[0] - tnew[0];
			b2[1] = b[1] - tnew[1];
			b2[2] = b[2] - tnew[2];
			matrix_product331(Rnew, b2, b_cam);
		}

		if (!init->known_intrinsics) {
			matrix_product331(K, b_cam, b_proj);
			p[0] = -b_proj[0] / b_proj[2];
			p[1] = -b_proj[1] / b_proj[2];
		}
		else {
			/* Apply intrinsics */
			double x_n = -b_cam[0] / b_cam[2];
			double y_n = -b_cam[1] / b_cam[2];

			double *k = init->k_known;
			double rsq = x_n * x_n + y_n * y_n;
			double factor = 1.0 + k[0] * rsq +
				k[1] * rsq * rsq + k[4] * rsq * rsq * rsq;

			double dx_x = 2 * k[2] * x_n * y_n + k[3] * (rsq + 2 * x_n * x_n);
			double dx_y = k[2] * (rsq + 2 * y_n * y_n) + 2 * k[3] * x_n * y_n;

			double x_d = x_n * factor + dx_x;
			double y_d = y_n * factor + dx_y;

			double *K = init->K_known;
			p[0] = K[0] * x_d + K[1] * y_d + K[2];
			p[1] = K[4] * y_d + K[5];
		}
	}



	void sfm_project_rd(camera_params_t *init, double *K, double *k,
		double *R, double *dt, double *b, double *p,
		int undistort, int explicit_camera_centers)
	{
		double *t;

		double tnew[3];
		double b_cam[3];

		t = init->t;

		tnew[0] = dt[0];
		tnew[1] = dt[1];
		tnew[2] = dt[2];

		/* Project! */
		if (!explicit_camera_centers) {
			matrix_product331(R, b, b_cam);
			b_cam[0] += tnew[0];
			b_cam[1] += tnew[1];
			b_cam[2] += tnew[2];
		}
		else {
			double b2[3];
			b2[0] = b[0] - tnew[0];
			b2[1] = b[1] - tnew[1];
			b2[2] = b[2] - tnew[2];
			matrix_product331(R, b2, b_cam);
		}

		if (!init->known_intrinsics) {
			p[0] = -b_cam[0] * K[0] / b_cam[2];
			p[1] = -b_cam[1] * K[0] / b_cam[2];
		}
		else {
			/* Apply intrinsics */
			double x_n = -b_cam[0] / b_cam[2];
			double y_n = -b_cam[1] / b_cam[2];

			double *k = init->k_known;
			double rsq = x_n * x_n + y_n * y_n;
			double factor = 1.0 + k[0] * rsq +
				k[1] * rsq * rsq + k[4] * rsq * rsq * rsq;

			double dx_x = 2 * k[2] * x_n * y_n + k[3] * (rsq + 2 * x_n * x_n);
			double dx_y = k[2] * (rsq + 2 * y_n * y_n) + 2 * k[3] * x_n * y_n;

			double x_d = x_n * factor + dx_x;
			double y_d = y_n * factor + dx_y;

			double *K = init->K_known;
			p[0] = K[0] * x_d + K[1] * y_d + K[2];
			p[1] = K[4] * y_d + K[5];
		}

		// p[0] = b_cam[0] * K[0] / b_cam[2];
		// p[1] = b_cam[1] * K[0] / b_cam[2];

		/* Apply radial distortion */
		if (undistort) {
#ifndef TEST_FOCAL
			double k1 = k[0], k2 = k[1];
#else
			double k1 = k[0] / init->k_scale;
			double k2 = k[1] / init->k_scale;
#endif

			double rsq = (p[0] * p[0] + p[1] * p[1]) / (K[0] * K[0]);
			double factor = 1.0 + k1 * rsq + k2 * rsq * rsq;

			p[0] *= factor;
			p[1] *= factor;
		}
	}





	static void *safe_malloc(int n, char *where)
	{
		void *mem = malloc(n);

		if (mem == NULL) {
			if (where) {
				printf("[safe_malloc] Error allocating %d bytes "
					"of memory at %s\n", n, where);
			}
			else {
				printf("[safe_malloc] Error allocating %d bytes of memory\n", n);
			}

			fflush(stdout);
			exit(1);
		}

		return mem;
	}

	static double *global_last_ws = NULL;
	static double *global_last_Rs = NULL;

	static void sfm_project_point3(int j, int i, double *aj, double *bi,
		double *xij, void *adata)
	{
		sfm_global_t *globs = (sfm_global_t *)adata;

		double K[9] =
		{ 1.0, 0.0, 0.0,
		0.0, 1.0, 0.0,
		0.0, 0.0, 1.0 };

		double *w, *dt, *k;

		/* Compute intrinsics */
		if (!globs->est_focal_length) {
			K[0] = K[4] = globs->init_params[j].f; // globs->global_params.f;
		}
		else if (globs->const_focal_length) {
			printf("Error: case of constant focal length "
				"has not been implemented.\n");
			K[0] = K[4] = globs->global_params.f;
		}
		else {
#ifndef TEST_FOCAL
			K[0] = K[4] = aj[6];
#else
			K[0] = K[4] = aj[6] / globs->init_params[j].f_scale;
#endif
		}

		/* Compute translation, rotation update */
		dt = aj + 0;
		w = aj + 3;

		if (globs->est_focal_length)
			k = aj + 7;
		else
			k = aj + 6;

		if (w[0] != global_last_ws[3 * j + 0] ||
			w[1] != global_last_ws[3 * j + 1] ||
			w[2] != global_last_ws[3 * j + 2]) {

			rot_update(globs->init_params[j].R, w, global_last_Rs + 9 * j);
			global_last_ws[3 * j + 0] = w[0];
			global_last_ws[3 * j + 1] = w[1];
			global_last_ws[3 * j + 2] = w[2];
		}

		sfm_project_rd(globs->init_params + j, K, k, global_last_Rs + 9 * j,
			dt, bi, xij, globs->estimate_distortion,
			globs->explicit_camera_centers);
	}

	static void sfm_project_point3_mot(int j, int i, double *aj,
		double *xij, void *adata)
	{
		sfm_global_t *globs = (sfm_global_t *)adata;
		double *b = globs->points[i].data();

		sfm_project_point3(j, i, aj, b, xij, adata);
	}


#define SBA_V121

	void bundle_adjustment(int num_pts, int num_cameras, int ncons,
		char *vmask,
		double *projections,
		int est_focal_length,
		int const_focal_length,
		int undistort,
		int explicit_camera_centers,
		camera_params_t *init_camera_params,
		vec3d *init_pts,
		int use_constraints,
		int use_point_constraints,
		vec3d *pt_constraints,
		double pt_constraint_weight,
		int fix_points,
		double eps2,
		double *Vout,
		double *Sout,
		double *Uout, double *Wout
		/* size num_cameras ** 2 * cnp * cnp */)
	{
		int cnp;
		double *params;

#ifdef SBA_V121
		double opts[6]; // opts[5];
#else
		double opts[3];
#endif
		double info[10];

		int i, j, base;
		int num_camera_params, num_pt_params, num_params;

		sfm_global_t global_params;

		// #ifndef SBA_V121
		camera_constraints_t *constraints = NULL;
		// #endif

		point_constraints_t *point_constraints = NULL;

		const double f_scale = 0.001;
		const double k_scale = 5.0;

		if (est_focal_length)
			cnp = 7;
		else
			cnp = 6;

		if (undistort)
			cnp += 2;

		num_camera_params = cnp * num_cameras;
		num_pt_params = 3 * num_pts;
		num_params = num_camera_params + num_pt_params;

		params = (double *)safe_malloc(sizeof(double) * num_params, "params");

		/* Fill parameters */
		for (j = 0; j < num_cameras; j++) {
			int c = 0;

#ifdef TEST_FOCAL
			init_camera_params[j].f_scale = f_scale;
			init_camera_params[j].k_scale = k_scale;
#else
			init_camera_params[j].f_scale = 1.0;
			init_camera_params[j].k_scale = 1.0;
#endif

			/* Translation is zero */
			params[cnp * j + 0] = init_camera_params[j].t[0]; // 0.0;
			params[cnp * j + 1] = init_camera_params[j].t[1]; // 0.0;
			params[cnp * j + 2] = init_camera_params[j].t[2]; // 0.0;

			/* Rotation is zero */
			params[cnp * j + 3] = 0.0;
			params[cnp * j + 4] = 0.0;
			params[cnp * j + 5] = 0.0;

			if (est_focal_length) {
				/* Focal length is initial estimate */
#ifndef TEST_FOCAL
				params[cnp * j + 6] = init_camera_params[j].f;
#else
				params[cnp * j + 6] =
					init_camera_params[j].f * init_camera_params[j].f_scale;
#endif
				c = 7;
			}
			else {
				c = 6;
			}

			if (undistort) {
#ifndef TEST_FOCAL
				params[cnp * j + c] = init_camera_params[j].k[0];
				params[cnp * j + c + 1] = init_camera_params[j].k[1];
#else
				double scale = init_camera_params[j].k_scale;
				params[cnp * j + c] = init_camera_params[j].k[0] * scale;
				params[cnp * j + c + 1] = init_camera_params[j].k[1] * scale;
#endif
			}
		}

		base = num_camera_params;
		for (i = 0; i < num_pts; i++) {
			params[base + 3 * i + 0] = init_pts[i].x;
			params[base + 3 * i + 1] = init_pts[i].y;
			params[base + 3 * i + 2] = init_pts[i].z;
		}

		opts[0] = 1.0e-3;
		opts[1] = 1.0e-10; // 1.0e-15;
		opts[2] = eps2; // 0.0;  // 1.0e-10; // 1.0e-15;

#ifdef SBA_V121
		opts[3] = 1.0e-12;
		// opts[4] = 4.0e-2;
		opts[4] = 0.0;
		opts[5] = 4.0e-2; // change this back to opts[4] for sba v1.2.1
#endif

		// opts[1] = 1.0e-8;
		// opts[2] = 1.0e-8;

		// #ifndef SBA_V121
		/* Create the constraints */
		if (use_constraints) {
			constraints =
				(camera_constraints_t *)
				malloc(num_cameras * sizeof(camera_constraints_t));

			for (i = 0; i < num_cameras; i++) {
				constraints[i].constrained = (char *)malloc(cnp);
				constraints[i].constraints =
					(double *)malloc(sizeof(double) * cnp);
				constraints[i].weights = (double *)malloc(sizeof(double) * cnp);

				memcpy(constraints[i].constrained,
					init_camera_params[i].constrained, cnp * sizeof(char));
				memcpy(constraints[i].constraints,
					init_camera_params[i].constraints, cnp * sizeof(double));
				memcpy(constraints[i].weights,
					init_camera_params[i].weights, cnp * sizeof(double));

#ifdef TEST_FOCAL
				if (est_focal_length) {
					constraints[i].constraints[6] *= f_scale;
					constraints[i].weights[6] *= (1.0 / (f_scale * f_scale));
				}

				if (undistort) {
					constraints[i].constraints[7] *= k_scale;
					constraints[i].weights[7] *= (1.0 / (k_scale * k_scale));

					constraints[i].constraints[8] *= k_scale;
					constraints[i].weights[8] *= (1.0 / (k_scale * k_scale));
				}
#endif
			}
		}
		// #endif

		if (use_point_constraints) {
			point_constraints =
				(point_constraints_t *)
				malloc(num_pts * sizeof(point_constraints_t));

			for (i = 0; i < num_pts; i++) {
				if (pt_constraints[i].x == 0.0 && pt_constraints[i].y == 0.0 && pt_constraints[i].z == 0.0) {
					point_constraints[i].constrained = 0;
					point_constraints[i].constraints[0] = 0.0;
					point_constraints[i].constraints[1] = 0.0;
					point_constraints[i].constraints[2] = 0.0;
					point_constraints[i].weight = 0.0;
				}
				else {
					// printf("[run_sfm] Constraining point %d\n", i);
					point_constraints[i].constrained = 1;
					point_constraints[i].weight = pt_constraint_weight;
					point_constraints[i].constraints[0] = pt_constraints[i].x;
					point_constraints[i].constraints[1] = pt_constraints[i].y;
					point_constraints[i].constraints[2] = pt_constraints[i].z;
				}
			}
		}

		/* Fill global param struct */
		global_params.num_cameras = num_cameras;
		global_params.num_points = num_pts;
		global_params.num_params_per_camera = cnp;

		global_params.est_focal_length = est_focal_length;
		global_params.const_focal_length = const_focal_length;
		global_params.estimate_distortion = undistort;
		global_params.explicit_camera_centers = explicit_camera_centers,

			global_params.global_params.f = 1.0;
		global_params.init_params = init_camera_params;

		global_last_ws =
			(double*)safe_malloc(3 * num_cameras * sizeof(double), "global_last_ws");

		global_last_Rs =
			(double*)safe_malloc(9 * num_cameras * sizeof(double), "global_last_ws");

		global_params.points = init_pts;

		for (i = 0; i < num_cameras; i++) {
			global_last_ws[3 * i + 0] = 0.0;
			global_last_ws[3 * i + 1] = 0.0;
			global_last_ws[3 * i + 2] = 0.0;

			memcpy(global_last_Rs + 9 * i,
				init_camera_params[i].R, 9 * sizeof(double));
		}

		/* Run sparse bundle adjustment */
#define MAX_ITERS 150 // 256
#define VERBOSITY 1

#ifdef SBA_V121
		if (fix_points == 0) {
				sba_motstr_levmar(num_pts, num_cameras, ncons,
					vmask, params, cnp, 3, projections, NULL, 2,
					//remove NULL in prev line for sba v1.2.1
					sfm_project_point3, NULL,
					(void *)(&global_params),
					MAX_ITERS, VERBOSITY, opts, info,
					use_constraints, constraints,
					use_point_constraints,
					point_constraints, Vout, Sout, Uout, Wout);
		}
		else {
				sba_mot_levmar(num_pts, num_cameras, ncons,
					vmask, params, cnp, projections, NULL, 2,
					sfm_project_point3_mot, NULL,
					(void *)(&global_params),
					MAX_ITERS, VERBOSITY, opts, info,
					use_constraints, constraints);
		}
#else
		if (fix_points == 0) {
			sba_motstr_levmar(num_pts, num_cameras, ncons,
				vmask, params, cnp, 3, projections, 2,
				sfm_project_point2, NULL, (void *)(&global_params),
				MAX_ITERS, VERBOSITY, opts, info,
				use_constraints, constraints,
				Vout, Sout, Uout, Wout);
		}
		else {
			sba_mot_levmar(num_pts, num_cameras, ncons,
				vmask, params, cnp, projections, 2,
				sfm_mot_project_point, NULL, (void *)(&global_params),
				MAX_ITERS, VERBOSITY, opts, info);
		}
#endif

		printf("[run_sfm] Number of iterations: %d\n", (int)info[5]);
		printf("info[6] = %0.3f\n", info[6]);

		/* Copy out the params */
		for (j = 0; j < num_cameras; j++) {
			double *dt = params + cnp * j + 0;
			double *w = params + cnp * j + 3;
			double Rnew[9];
			int c;

			/* Translation */
			init_camera_params[j].t[0] = dt[0];
			init_camera_params[j].t[1] = dt[1];
			init_camera_params[j].t[2] = dt[2];
			// init_camera_params[j].t[0] += dt[0];
			// init_camera_params[j].t[1] += dt[1];
			// init_camera_params[j].t[2] += dt[2];

			/* Rotation */
			rot_update(init_camera_params[j].R, w, Rnew);
			memcpy(init_camera_params[j].R, Rnew, 9 * sizeof(double));

			/* Focal length */
			if (est_focal_length) {
				c = 7;
#ifndef TEST_FOCAL
				init_camera_params[j].f = params[cnp * j + 6];
#else
				init_camera_params[j].f =
					params[cnp * j + 6] / init_camera_params[j].f_scale;
#endif
			}
			else {
				c = 6;
			}

			if (undistort) {
#ifndef TEST_FOCAL
				init_camera_params[j].k[0] = params[cnp * j + c];
				init_camera_params[j].k[1] = params[cnp * j + c + 1];
#else
				double scale = init_camera_params[j].k_scale;
				init_camera_params[j].k[0] = params[cnp * j + c] / scale;
				init_camera_params[j].k[1] = params[cnp * j + c + 1] / scale;
#endif
			}

#ifdef TEST_FOCAL
			init_camera_params[j].f_scale = 1.0;
			init_camera_params[j].k_scale = 1.0;
#endif
		}

		base = num_camera_params;
		for (i = 0; i < num_pts; i++) {
			init_pts[i].x = params[base + 3 * i + 0];
			init_pts[i].y = params[base + 3 * i + 1];
			init_pts[i].z = params[base + 3 * i + 2];
		}

		// #define DEBUG_SFM
#ifdef DEBUG_SFM
		for (i = 0; i < num_cameras; i++) {
			int num_projs = 0;
			double error = 0.0;
			double error_max = 0.0;
			int idx_max = 0;
			double px_max = 0.0, py_max = 0.0;

			double K[9] = { init_camera_params[i].f, 0.0, 0.0,
				0.0, init_camera_params[i].f, 0.0,
				0.0, 0.0, 1.0 };
			double w[3] = { 0.0, 0.0, 0.0 };
			double dt[3] = { init_camera_params[i].t[0],
				init_camera_params[i].t[1],
				init_camera_params[i].t[2] };

			// double dt[3] = { 0.0, 0.0, 0.0 };

			for (j = 0; j < num_pts; j++) {
				double b[3], pr[2];
				double dx, dy, dist;

				if (!vmask[j * num_cameras + i])
					continue;

				b[0] = Vx(init_pts[j]);
				b[1] = Vy(init_pts[j]);
				b[2] = Vz(init_pts[j]);

				sfm_project(&(init_camera_params[i]), K, w, dt, b, pr,
					global_params.explicit_camera_centers);

				dx = pr[0] - Vx(projections[j * num_cameras + i]);
				dy = pr[1] - Vy(projections[j * num_cameras + i]);

				dist = dx * dx + dy * dy;
				error += dist;

				if (dist > error_max) {
					idx_max = j;
					error_max = dist;
					px_max = Vx(projections[j * num_cameras + i]);
					py_max = Vy(projections[j * num_cameras + i]);
				}

				num_projs++;
			}

			printf("Camera %d:  error = %0.3f (%0.3f)\n", i,
				error, sqrt(error / num_projs));
			printf("           error_max = %0.3f (%d)\n", sqrt(error_max), idx_max);
			printf("           proj = %0.3f, %0.3f\n", px_max, py_max);
		}
#endif /* DEBUG_SFM */

		free(params);

		// #ifndef SBA_V121
		if (use_constraints) {
			for (i = 0; i < num_cameras; i++) {
				free(constraints[i].constraints);
				free(constraints[i].constrained);
				free(constraints[i].weights);
			}
			free(constraints);
		}

		free(global_last_ws);
		free(global_last_Rs);

		// #endif
	}


}













#include "../mvglib/qsort.h"
#include "../mvglib/triangulate.h"
namespace sfm {

	void initialize_camera_parameters(const ImageData &data, camera_params_t &camera)
	{
		matrix_ident(3, camera.R);
		camera.t[0] = camera.t[1] = camera.t[2] = 0.0;
		camera.f = 0.0;
		camera.k[0] = camera.k[1] = 0.0;

		camera.k_inv[0] = camera.k_inv[2] = camera.k_inv[3] = 0.0;
		camera.k_inv[4] = camera.k_inv[5] = 0.0;
		camera.k_inv[1] = 1.0;

		camera.f_scale = 1.0;
		camera.k_scale = 1.0;

		for (int i = 0; i < NUM_CAMERA_PARAMS; i++) {
			camera.constrained[i] = 0;
			camera.constraints[i] = 0.0;
			camera.weights[i] = 0.0;
		}

		if (data.known_intrinsics) {
			camera.known_intrinsics = 1;
			memcpy(camera.K_known, data.K, 9 * sizeof(double));
			memcpy(camera.k_known, data.k, 5 * sizeof(double));
		}
		else {
			camera.known_intrinsics = 0;
		}
	}



	static int compare_doubles(const void *d1, const void *d2)
	{
		double a = *(double *)d1;
		double b = *(double *)d2;

		if (a < b) return -1;
		if (a > b) return 1;
		return 0;
	}

	std::vector<int> refine_camera_parameters(const ImageData &data,
		int num_points,
		vec3d *points, vec2d *projs,
		int *pt_idxs, camera_params_t *camera,
		double *error_out,
		bool adjust_focal,
		bool remove_outliers,
		//                                         bool optimize_for_fisheye,
		bool estimate_distortion,
		double min_proj_error_threshold,
		double max_proj_error_threshold)
	{
		int num_points_curr = num_points;
		vec3d *points_curr = new vec3d[num_points];
		vec2d *projs_curr = new vec2d[num_points];

		memcpy(points_curr, points, num_points * sizeof(vec3d));
		memcpy(projs_curr, projs, num_points * sizeof(vec2d));

		std::vector<int> inliers;

		for (int i = 0; i < num_points; i++)
			inliers.push_back(i);

		int round = 0;

		/* First refine with the focal length fixed */
		camera_refine(num_points_curr, points_curr, projs_curr, camera, 0, 0);

		while (1) {
			printf("[RefineCameraParameters] Calling with %d points\n",
				num_points_curr);

			camera_refine(num_points_curr, points_curr, projs_curr, camera,
				adjust_focal ? 1 : 0, estimate_distortion ? 1 : 0);

			if (!remove_outliers)
				break;

			vec3d *points_next = new vec3d[num_points];
			vec2d *projs_next = new vec2d[num_points];

			int count = 0;
			double error = 0.0;
			std::vector<int> inliers_next;

			double *errors = new double[num_points_curr];

			for (int i = 0; i < num_points_curr; i++) {
				vec2d pr = sfm_project_final(camera, points_curr[i], 1,
					estimate_distortion ? 1 : 0);

				double dx = pr.x - projs_curr[i].x;
				double dy = pr.y - projs_curr[i].y;
				double diff = sqrt(dx * dx + dy * dy);

				errors[i] = diff;
				error += diff;
			}

			printf("[RefineCameraParameters] Error: %0.3f\n",
				error / num_points_curr);

			/* Sort and histogram errors */
			double med = kth_element_copy(num_points_curr,
				iround(0.95 * num_points_curr),
				errors);

			/* We will tolerate any match with projection error < 8.0 */
			double threshold = 1.2 * NUM_STDDEV * med; /* k * stddev */
			ogf_clamp(threshold, min_proj_error_threshold,
				max_proj_error_threshold);

			printf("[RefineCameraParameters] Threshold = %0.3f\n", threshold);
			for (int i = 0; i < num_points_curr; i++) {
				if (errors[i] < threshold) {
					inliers_next.push_back(inliers[i]);

					points_next[count] = points_curr[i];
					projs_next[count] = projs_curr[i];

					count++;
				}
				else {
					if (pt_idxs != NULL) {
						printf("[RefineCameraParameters] Removing point [%d] with "
							"reprojection error %0.3f\n", pt_idxs[i],
							errors[i]);
					}
					else {
						printf("[RefineCameraParameters] Removing point with "
							"reprojection error %0.3f\n", errors[i]);
					}
				}
			}

#if 1
			qsort(errors, num_points_curr, sizeof(double), compare_doubles);

			double pr_min = errors[0];
			double pr_max = errors[num_points_curr - 1];
			double pr_step = (pr_max - pr_min) / NUM_ERROR_BINS;

			/* Break histogram into 10 bins */
			int idx_count = 0;
			for (int i = 0; i < NUM_ERROR_BINS; i++) {
				double max = pr_min + (i + 1) * pr_step;
				int start = idx_count;

				while (idx_count < num_points_curr && errors[idx_count] <= max)
					idx_count++;

				int bin_size = idx_count - start;
				// 			printf("   E[%0.3e--%0.3e]: %d [%0.3f]\n",
				// 				max - pr_step, max, bin_size,
				// 				bin_size / (double)num_points_curr);
			}
#endif

			delete[] points_curr;
			delete[] projs_curr;
			delete[] errors;

			points_curr = points_next;
			projs_curr = projs_next;

			if (count == num_points_curr)
				break;  /* We're done */

			num_points_curr = count;

			inliers = inliers_next;

			if (count == 0) /* Out of measurements */
				break;

			round++;

			if (error_out != NULL) {
				*error_out = error;
			}
		}

		printf("[RefineCameraParameters] Exiting after %d rounds "
			"with %d / %d points\n", round + 1, num_points_curr, num_points);

		delete[] points_curr;
		delete[] projs_curr;

		return inliers;
	}


	void fix_intrinsics(double *P, double *K, double *R, double *t) {
		/* Check the parity along the diagonal */
		int neg = (K[0] < 0.0) + (K[4] < 0.0) + (K[8] < 0.0);

		/* If odd parity, negate the instrinsic matrix */
		if ((neg % 2) == 1) {
			matrix_scale(3, 3, K, -1.0, K);
			matrix_scale(3, 4, P, -1.0, P);
		}

		/* Now deal with case of even parity */
		double fix[9];
		matrix_ident(3, fix);
		double tmp[9], tmp2[12];

		if (K[0] < 0.0 && K[4] < 0.0) {
			fix[0] = -1.0;
			fix[4] = -1.0;
		}
		else if (K[0] < 0.0) {
			fix[0] = -1.0;
			fix[8] = -1.0;
		}
		else if (K[4] < 0.0) {
			fix[4] = -1.0;
			fix[8] = -1.0;
		}
		else {
			/* No change needed */
		}

		matrix_product(3, 3, 3, 3, K, fix, tmp);
		memcpy(K, tmp, sizeof(double) * 3 * 3);

		double Kinv[9];
		matrix_invert(3, K, Kinv);

		matrix_product(3, 3, 3, 4, Kinv, P, tmp2);

		memcpy(R + 0, tmp2 + 0, sizeof(double) * 3);
		memcpy(R + 3, tmp2 + 4, sizeof(double) * 3);
		memcpy(R + 6, tmp2 + 8, sizeof(double) * 3);

		t[0] = tmp2[3];
		t[1] = tmp2[7];
		t[2] = tmp2[11];
	}



	bool find_and_verify_camera(int num_points, vec3d *points_solve, vec2d *projs_solve,
		int *idxs_solve,
		double *K, double *R, double *t,
		double proj_estimation_threshold,
		double proj_estimation_threshold_weak,
		std::vector<int> &inliers,
		std::vector<int> &inliers_weak,
		std::vector<int> &outliers)
	{
		/* First, find the projection matrix */
		double P[12];
		int r = -1;

		if (num_points >= 9) {
			r = find_projection_3x4_ransac(num_points,
				points_solve, projs_solve,
				P, /* 2048 */ 4096 /* 100000 */,
				proj_estimation_threshold);
		}

		if (r == -1) {
			printf("[FindAndVerifyCamera] Couldn't find projection matrix\n");
			return false;
		}

		/* If number of inliers is too low, fail */
		if (r <= MIN_INLIERS_EST_PROJECTION) {
			printf("[FindAndVerifyCamera] Too few inliers to use "
				"projection matrix\n");
			return false;
		}

		double KRinit[9], Kinit[9], Rinit[9], tinit[3];
		memcpy(KRinit + 0, P + 0, 3 * sizeof(double));
		memcpy(KRinit + 3, P + 4, 3 * sizeof(double));
		memcpy(KRinit + 6, P + 8, 3 * sizeof(double));

		dgerqf_driver(3, 3, KRinit, Kinit, Rinit);

		/* We want our intrinsics to have a certain form */
		fix_intrinsics(P, Kinit, Rinit, tinit);
		matrix_scale(3, 3, Kinit, 1.0 / Kinit[8], Kinit);

		//printf("[FindAndVerifyCamera] Estimated intrinsics:\n");
		//matrix_print(3, 3, Kinit);
		//printf("[FindAndVerifyCamera] Estimated extrinsics:\n");
		//matrix_print(3, 3, Rinit);
		//matrix_print(1, 3, tinit);
		//fflush(stdout);

		/* Check cheirality constraint */
		//printf("[FindAndVerifyCamera] Checking consistency...\n");

		double Rigid[12] =
		{ Rinit[0], Rinit[1], Rinit[2], tinit[0],
		Rinit[3], Rinit[4], Rinit[5], tinit[1],
		Rinit[6], Rinit[7], Rinit[8], tinit[2] };

		int num_behind = 0;
		for (int j = 0; j < num_points; j++) {
			double p[4] = {
				points_solve[j].x,
				points_solve[j].y,
				points_solve[j].z, 1.0 };
			double q[3], q2[3];

			matrix_product(3, 4, 4, 1, Rigid, p, q);
			matrix_product331(Kinit, q, q2);

			double pimg[2] = { -q2[0] / q2[2], -q2[1] / q2[2] };
			double diff =
				(pimg[0] - projs_solve[j].x) *
				(pimg[0] - projs_solve[j].x) +
				(pimg[1] - projs_solve[j].y) *
				(pimg[1] - projs_solve[j].y);

			diff = sqrt(diff);

			if (diff < proj_estimation_threshold)
				inliers.push_back(j);

			if (diff < proj_estimation_threshold_weak) {
				inliers_weak.push_back(j);
			}
			else {
				printf("[FindAndVerifyCamera] Removing point [%d] "
					"(reproj. error = %0.3f)\n", idxs_solve[j], diff);
				outliers.push_back(j);
			}

			// EDIT!!!
			if (q[2] > 0.0)
				num_behind++;  /* Cheirality constraint violated */
		}

		if (num_behind >= 0.9 * num_points) {
			printf("[FindAndVerifyCamera] Error: camera is pointing "
				"away from scene\n");
			return false;
		}

		memcpy(K, Kinit, sizeof(double) * 9);
		memcpy(R, Rinit, sizeof(double) * 9);
		memcpy(t, tinit, sizeof(double) * 3);

		// #define COLIN_HACK
#ifdef COLIN_HACK
		matrix_ident(3, R);
		t[0] = t[1] = t[2] = 0.0;
#endif

		return true;
	}


	/* Triangulate two points */
	vec3d triangulate(vec2d p, vec2d q,
		camera_params_t c1, camera_params_t c2,
		double &proj_error, bool &in_front, double &angle,
		bool explicit_camera_centers)
	{
		double K1[9], K2[9];
		double K1inv[9], K2inv[9];

		get_intrinsics(c1, K1);
		get_intrinsics(c2, K2);

		matrix_invert(3, K1, K1inv);
		matrix_invert(3, K2, K2inv);

		/* Set up the 3D point */
		// EDIT!!!
		double proj1[3] = { p.x, p.y, -1.0 };
		double proj2[3] = { q.x, q.y, -1.0 };

		double proj1_norm[3], proj2_norm[3];

		matrix_product(3, 3, 3, 1, K1inv, proj1, proj1_norm);
		matrix_product(3, 3, 3, 1, K2inv, proj2, proj2_norm);

		vec2d p_norm(proj1_norm[0] / proj1_norm[2], proj1_norm[1] / proj1_norm[2]);
		vec2d q_norm(proj2_norm[0] / proj2_norm[2], proj2_norm[1] / proj2_norm[2]);

		/* Undo radial distortion */
		p_norm = undistort_normalized_point(p_norm, c1);
		q_norm = undistort_normalized_point(q_norm, c2);

		/* Compute the angle between the rays */
		angle = compute_ray_angle(p, q, c1, c2);

		/* Triangulate the point */
		vec3d pt;
		if (!explicit_camera_centers) {
			pt = triangulate(p_norm, q_norm, c1.R, c1.t, c2.R, c2.t, &proj_error);
		}
		else {
			double t1[3];
			double t2[3];

			/* Put the translation in standard form */
			matrix_product(3, 3, 3, 1, c1.R, c1.t, t1);
			matrix_scale(3, 1, t1, -1.0, t1);
			matrix_product(3, 3, 3, 1, c2.R, c2.t, t2);
			matrix_scale(3, 1, t2, -1.0, t2);

			pt = triangulate(p_norm, q_norm, c1.R, t1, c2.R, t2, &proj_error);
		}

		proj_error = (c1.f + c2.f) * 0.5 * sqrt(proj_error * 0.5);

		/* Check cheirality */
		bool cc1 = check_cheirality(pt, c1);
		bool cc2 = check_cheirality(pt, c2);

		in_front = (cc1 && cc2);

		return pt;
	}



}