#include "epipolar.h"

#include "../mvglib/5point.h"
#include "../math/math_types.h"
#include "../math/matrix_driver.h"

#include "../mvglib/fmatrix.h"

#include <cassert>



namespace sfm {

	/* Estimate an E-matrix from a given set of point matches */
	std::vector<int> estimate_e_matrix(const std::vector<Keypoint> &k1,
		const std::vector<Keypoint> &k2,
		const std::vector<KeypointMatch>& matches,
		int num_trials, double threshold,
		double f1, double f2,
		double *E, double *F)
	{
		int num_keys1 = (int)k1.size();
		int num_keys2 = (int)k2.size();

		std::vector<Keypoint> k1_norm, k2_norm;
		k1_norm.resize(num_keys1);
		k2_norm.resize(num_keys2);

		for (int i = 0; i < num_keys1; i++) {
			Keypoint k;
			k.x = float(k1[i].x / f1);
			k.y = float(k1[i].y / f1);
			k1_norm[i] = k;
		}

		for (int i = 0; i < num_keys2; i++) {
			Keypoint k;
			k.x = float(k2[i].x / f2);
			k.y = float(k2[i].y / f2);
			k2_norm[i] = k;
		}

		double scale = 0.5 * (f1 + f2);

		std::vector<int> inliers =
			estimate_f_Matrix(k1_norm, k2_norm, matches, num_trials, threshold / (scale * scale), E, true);

		double K1_inv[9] = { 1.0 / f1, 0.0, 0.0,
			0.0, 1.0 / f1, 0.0,
			0.0, 0.0, 1.0 };
		double K2_inv[9] = { 1.0 / f2, 0.0, 0.0,
			0.0, 1.0 / f2, 0.0,
			0.0, 0.0, 1.0 };

		double tmp[9];
		matrix_product(3, 3, 3, 3, K1_inv, E, tmp);
		matrix_product(3, 3, 3, 3, K2_inv, tmp, F);

		return inliers;
	}


	/* Estimate an F-matrix from a given set of point matches */
	std::vector<int> estimate_f_Matrix(
		const std::vector<Keypoint> &k1,
		const std::vector<Keypoint> &k2,
		const std::vector<KeypointMatch>& matches,
		int num_trials, double threshold,
		double *F, bool essential
		)
	{
		int num_pts = (int)matches.size();

		/* num_pts should be greater than a threshold */
		if (num_pts < 20) {
			std::vector<int> inliers;
			return inliers;
		}

		vec3d* k1_pts = new vec3d[num_pts];
		vec3d* k2_pts = new vec3d[num_pts];

		vec3d* k1_pts_in = new vec3d[num_pts];
		vec3d* k2_pts_in = new vec3d[num_pts];

		for (int i = 0; i < num_pts; i++) {
			int idx1 = matches[i].key_idx1;
			int idx2 = matches[i].key_idx2;

			assert(idx1 < (int)k1.size());
			assert(idx2 < (int)k2.size());

			k1_pts[i] = vec3d(k1[idx1].x, k1[idx1].y, 1.0);
			k2_pts[i] = vec3d(k2[idx2].x, k2[idx2].y, 1.0);
		}

		estimate_fmatrix_ransac_matches(
			num_pts, k2_pts, k1_pts, num_trials, threshold, 0.95, (essential ? 1 : 0), F
			);

		/* Find the inliers */
		std::vector<int> inliers;

		for (int i = 0; i < num_pts; i++) {
			double dist = fmatrix_compute_residual(F, k2_pts[i], k1_pts[i]);
			if (dist < threshold) {
				inliers.push_back(i);
			}
		}

		/* Re-estimate using inliers */
		int num_inliers = (int)inliers.size();

		for (int i = 0; i < num_inliers; i++) {
			k1_pts_in[i] = k1_pts[inliers[i]]; // v3_new(k1[idx1]->m_x, k1[idx1]->m_y, 1.0);
			k2_pts_in[i] = k2_pts[inliers[i]]; // v3_new(k2[idx2]->m_x, k2[idx2]->m_y, 1.0);
		}

		double F0[9];
		memcpy(F0, F, sizeof(double) * 9);

		if (!essential) {
			/* Refine using NLLS */
			for (int i = 0; i < num_inliers; i++) {
				k1_pts_in[i] = k1_pts[inliers[i]];
				k2_pts_in[i] = k2_pts[inliers[i]];
			}

			refine_fmatrix_nonlinear_matches(
				num_inliers, k2_pts_in, k1_pts_in, F0, F
				);
		}
		else {
			memcpy(F, F0, sizeof(double) * 9);
		}

		inliers.clear();
		for (int i = 0; i < num_pts; i++) {
			double dist = fmatrix_compute_residual(F, k2_pts[i], k1_pts[i]);
			if (dist < threshold) {
				inliers.push_back(i);
			}
		}
		num_inliers = (int)inliers.size();

		delete[] k1_pts;
		delete[] k2_pts;
		delete[] k1_pts_in;
		delete[] k2_pts_in;

		return inliers;
	}


	/* Estimate relative pose from a given set of point matches */
	int estimate_pose_5_point(
		const std::vector<Keypoint> &k1,
		const std::vector<Keypoint> &k2,
		const std::vector<KeypointMatch>& matches,
		int num_trials, double threshold,
		double *K1, double *K2,
		double *R, double *t
		)
	{
		int num_pts = (int)matches.size();

		vec2d* k1_pts = new vec2d[num_pts];
		vec2d* k2_pts = new vec2d[num_pts];

		for (int i = 0; i < num_pts; i++) {
			int idx1 = matches[i].key_idx1;
			int idx2 = matches[i].key_idx2;

			k1_pts[i] = vec2d(k1[idx1].x, k1[idx1].y);
			k2_pts[i] = vec2d(k2[idx2].x, k2[idx2].y);
		}

		int num_inliers = compute_pose_ransac(
			num_pts, k1_pts, k2_pts, K1, K2, threshold, num_trials, R, t
			);

		delete[] k1_pts;
		delete[] k2_pts;

		return num_inliers;
	}


}