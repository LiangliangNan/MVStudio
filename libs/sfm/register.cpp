
#include "register.h"
#include "../math/math_types.h"
#include "../math/matrix_driver.h"
#include "../mvglib/homography.h"
#include "../mvglib/horn.h"




namespace sfm {


	static int CountInliers(const std::vector<Keypoint> &k1,
		const std::vector<Keypoint> &k2,
		std::vector<KeypointMatch> matches,
		double *M, double thresh, std::vector<int> &inliers);

	static int LeastSquaresFit(const std::vector<Keypoint> &k1,
		const std::vector<Keypoint> &k2,
		std::vector<KeypointMatch> matches, MotionModel mm,
		const std::vector<int> &inliers, double *M);

	/* Estimate a transform between two sets of keypoints */
	std::vector<int> estimate_transform(
		const std::vector<Keypoint> &k1,
		const std::vector<Keypoint> &k2,
		const std::vector<KeypointMatch> &matches,
		MotionModel mm,
		int nRANSAC, double RANSACthresh,
		double *Mout)
	{
		int min_matches = -1;
		switch (mm) {
		case MotionRigid:
			min_matches = 3;
			break;
		case MotionHomography:
			min_matches = 4;
			break;
		}

		int *match_idxs = new int[min_matches];

		int num_matches = (int)matches.size();
		int max_inliers = 0;
		double Mbest[9];

		if (num_matches < min_matches) {
			std::vector<int> empty;
			printf("Cannot estimate rigid transform\n");
			return empty;
		}

		vec3d* r_pts = new vec3d[min_matches];
		vec3d* l_pts = new vec3d[min_matches];
		double *weight = new double[min_matches];

		for (int round = 0; round < nRANSAC; round++) {
			for (int i = 0; i < min_matches; i++) {
				bool found;
				int idx;

				do {
					found = true;
					idx = rand() % num_matches;

					for (int j = 0; j < i; j++) {
						if (match_idxs[j] == idx) {
							found = false;
							break;
						}
					}
				} while (!found);

				match_idxs[i] = idx;
			}

			/* Solve for the motion */

			for (int i = 0; i < min_matches; i++) {
				int idx1 = matches[match_idxs[i]].key_idx1;
				int idx2 = matches[match_idxs[i]].key_idx2;

				l_pts[i].x = k1[idx1].x;
				l_pts[i].y = k1[idx1].y;
				l_pts[i].z = 1.0;

				r_pts[i].x = k2[idx2].x;
				r_pts[i].y = k2[idx2].y;
				r_pts[i].z = 1.0;

				weight[i] = 1.0;
			}

			double Mcurr[9];

			switch (mm)	{
			case MotionRigid:	{
				double R[9], T[9], Tout[9], scale;
				align_horn(min_matches, r_pts, l_pts, R, T, Tout, &scale, weight);
				memcpy(Mcurr, Tout, 9 * sizeof(double));
				break;
			}

			case MotionHomography: {
				align_homography(min_matches, r_pts, l_pts, Mcurr, 0);
				break;
			}
			}


			std::vector<int> inliers;
			int num_inliers = CountInliers(k1, k2, matches, Mcurr,
				RANSACthresh, inliers);

			if (num_inliers > max_inliers) {
				max_inliers = num_inliers;
				memcpy(Mbest, Mcurr, 9 * sizeof(double));
			}
		}

		std::vector<int> inliers;
		CountInliers(k1, k2, matches, Mbest, RANSACthresh, inliers);
		memcpy(Mout, Mbest, 9 * sizeof(double));
		LeastSquaresFit(k1, k2, matches, mm, inliers, Mout);

		// memcpy(Mout, Mbest, 9 * sizeof(double));

		delete[] match_idxs;
		delete[] r_pts;
		delete[] l_pts;
		delete[] weight;

		return inliers;
	}

	static int CountInliers(const std::vector<Keypoint> &k1,
		const std::vector<Keypoint> &k2,
		std::vector<KeypointMatch> matches,
		double *M, double thresh, std::vector<int> &inliers)
	{
		inliers.clear();
		int count = 0;

		for (unsigned int i = 0; i < matches.size(); i++) {
			/* Determine if the ith feature in f1, when transformed by M,
			 * is within RANSACthresh of its match in f2 (if one exists)
			 *
			 * if so, increment count and append i to inliers */

			double p[3];

			p[0] = k1[matches[i].key_idx1].x;
			p[1] = k1[matches[i].key_idx1].y;
			p[2] = 1.0;

			double q[3];
			matrix_product(3, 3, 3, 1, M, p, q);

			double qx = q[0] / q[2];
			double qy = q[1] / q[2];

			double dx = qx - k2[matches[i].key_idx2].x;
			double dy = qy - k2[matches[i].key_idx2].y;

			double dist = sqrt(dx * dx + dy * dy);

			if (dist <= thresh) {
				count++;
				inliers.push_back(i);
			}
		}

		return count;
	}

	static int LeastSquaresFit(const std::vector<Keypoint> &k1,
		const std::vector<Keypoint> &k2,
		std::vector<KeypointMatch> matches, MotionModel mm,
		const std::vector<int> &inliers, double *M)
	{
		vec3d* r_pts = new vec3d[inliers.size()];
		vec3d* l_pts = new vec3d[inliers.size()];
		double *weight = new double[inliers.size()];

		/* Compute residual */
		double error = 0.0;
		for (int i = 0; i < (int)inliers.size(); i++) {
			int idx1 = matches[inliers[i]].key_idx1;
			int idx2 = matches[inliers[i]].key_idx2;

			double r[3], l[3];
			l[0] = k1[idx1].x;
			l[1] = k1[idx1].y;
			l[2] = 1.0;

			r[0] = k2[idx2].x;
			r[1] = k2[idx2].y;
			r[2] = 1.0;

			double rp[3];
			matrix_product(3, 3, 3, 1, M, l, rp);

			rp[0] /= rp[2];
			rp[1] /= rp[2];

			double dx = rp[0] - r[0];
			double dy = rp[1] - r[1];

			error += dx * dx + dy * dy;
		}

		//printf("[LeastSquaresFit] Residual error (before) is %0.3e\n", error);


		for (int i = 0; i < (int)inliers.size(); i++) {
			int idx1 = matches[inliers[i]].key_idx1;
			int idx2 = matches[inliers[i]].key_idx2;

			l_pts[i].x = k1[idx1].x;
			l_pts[i].y = k1[idx1].y;
			l_pts[i].z = 1.0;

			r_pts[i].x = k2[idx2].x;
			r_pts[i].y = k2[idx2].y;
			r_pts[i].z = 1.0;

			weight[i] = 1.0;
		}

		switch (mm) {
		case MotionRigid: {
			double R[9], T[9], Tout[9], scale;
			align_horn((int)inliers.size(), r_pts, l_pts, R, T, Tout, &scale, weight);
			memcpy(M, Tout, 9 * sizeof(double));
			break;
		}

		case MotionHomography: {
			align_homography((int)inliers.size(), r_pts, l_pts, M, 1);
			break;
		}
		}

		/* Compute residual */
		error = 0.0;
		for (int i = 0; i < (int)inliers.size(); i++) {
			int idx1 = matches[inliers[i]].key_idx1;
			int idx2 = matches[inliers[i]].key_idx2;

			double r[3], l[3];
			l[0] = k1[idx1].x;
			l[1] = k1[idx1].y;
			l[2] = 1.0;

			r[0] = k2[idx2].x;
			r[1] = k2[idx2].y;
			r[2] = 1.0;

			double rp[3];
			matrix_product(3, 3, 3, 1, M, l, rp);

			rp[0] /= rp[2];
			rp[1] /= rp[2];

			double dx = rp[0] - r[0];
			double dy = rp[1] - r[1];

			error += dx * dx + dy * dy;
		}

		//printf("[LeastSquaresFit] Residual error (after) is %0.3e\n", error);

		delete[] r_pts;
		delete[] l_pts;
		delete[] weight;

		return 0;
	}


}
