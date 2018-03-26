/*
 *  Copyright (c) 2008-2010  Noah Snavely (snavely (at) cs.cornell.edu)
 *    and the University of Washington
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 */

/* RelativePose.cpp */
/* Functions for computing the relative pose of two images */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "BundlerApp.h"

#include "Epipolar.h"
#include "Register.h"
#include "Bundle.h"

#include "sfm.h"
#include "fmatrix.h"
#include "matrix.h"
#include "triangulate.h"



bool BundlerApp::EstimateRelativePose2(int i1, int i2,
	camera_params_t &camera1,
	camera_params_t &camera2)
{
	// int num_images = GetNumImages();
	MatchIndex list_idx;

	if (i1 < i2)
		list_idx = GetMatchIndex(i1, i2); // i1 * num_images + i2;
	else
		list_idx = GetMatchIndex(i2, i1); // i2 * num_images + i1;

	std::vector<KeypointMatch> &matches = m_matches.GetMatchList(list_idx);
	// int num_matches = (int) m_match_lists[list_idx].size();
	int num_matches = (int)matches.size();

	double K1[9], K2[9];
	GetIntrinsics(camera1, K1);
	GetIntrinsics(camera2, K2);

	double R0[9], t0[3];
	int num_inliers = 0;

	num_inliers =
		EstimatePose5Point(m_image_data[i1].m_keys,
		m_image_data[i2].m_keys,
		matches,
		512, /* m_fmatrix_rounds, 8 * m_fmatrix_rounds */
		0.25 * m_fmatrix_threshold, // 0.003, // 0.004 /*0.001,*/ // /*0.5 **/ m_fmatrix_threshold, 
		K1, K2, R0, t0);

	if (num_inliers == 0)
		return false;

	printf("  Found %d / %d inliers (%0.3f%%)\n", num_inliers, num_matches,
		100.0 * num_inliers / num_matches);

	bool initialized = false;
	if (!initialized) {
		memcpy(camera2.R, R0, sizeof(double) * 9);

		matrix_transpose_product(3, 3, 3, 1, R0, t0, camera2.t);
		matrix_scale(3, 1, camera2.t, -1.0, camera2.t);
	}

	return true;
}
