
/* Routines for computing epipolar geometry */

#ifndef _SFM_EPIPOLAR_H_
#define _SFM_EPIPOLAR_H_


#include "keys.h"
#include <vector>


namespace sfm {


	/* Estimate an E-matrix from a given set of point matches */
	std::vector<int> estimate_e_matrix(
		const std::vector<Keypoint> &k1,
		const std::vector<Keypoint> &k2,
		const std::vector<KeypointMatch>& matches,
		int num_trials, double threshold,
		double f1, double f2,
		double *E, double *F
		);

	/* Estimate an F-matrix from a given set of point matches */
	std::vector<int> estimate_f_Matrix(
		const std::vector<Keypoint> &k1,
		const std::vector<Keypoint> &k2,
		const std::vector<KeypointMatch>& matches,
		int num_trials, double threshold,
		double *F, bool essential = false
		);

	/* Estimate relative pose from a given set of point matches */
	int estimate_pose_5_point(
		const std::vector<Keypoint> &k1,
		const std::vector<Keypoint> &k2,
		const std::vector<KeypointMatch>& matches,
		int num_trials, double threshold,
		double *K1, double *K2,
		double *R, double *t
		);


}


#endif /* _SFM_EPIPOLAR_H_ */
