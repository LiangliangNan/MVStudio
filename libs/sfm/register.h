
/* Compute relationships between images */

#ifndef _SFM_REGISTER_H_
#define _SFM_REGISTER_H_

#include "keys.h"
#include <vector>


namespace sfm {

	enum MotionModel {
		MotionRigid,
		MotionHomography
	};

	/* Estimate a transform between two sets of key points */
	std::vector<int> estimate_transform(
		const std::vector<Keypoint> &k1,
		const std::vector<Keypoint> &k2,
		const std::vector<KeypointMatch> &matches,
		MotionModel mm,
		int nRANSAC, double RANSACthresh,
		double *Mout
		);

}

#endif /* _SFM_REGISTER_H_ */
