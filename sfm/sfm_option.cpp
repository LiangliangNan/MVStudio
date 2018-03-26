#include "sfm_option.h"


SfmOption::SfmOption() {
	output_directory = "./sfm";
	image_directory = "./images";
	key_directory = "./keys";

	match_table_file = "matche_table.txt";
	bundle_output_file = "bundle.out";
	bundle_output_base = "bundle_";

	skip_full_bundle = false;

	// Should we fix the focal length for all cameras (set to init_focal_length)?
	fixed_focal_length = false;
	init_focal_length = 532.0;

	// Should we use camera constraints?
	use_camera_constraints = false;

	// Should we add soft constraint on focal lengths to stay near their estimated values
	constrain_focal = true;
	// The weight for the focal length constraint
	constrain_focal_weight = 0.0001;

	// Should we estimate radial distortion parameters (2 coefficients) for each camera?
	estimate_distortion = true;

	/* Point constraint fields */
	use_point_constraints = false;
	point_constraint_weight = 0.0;

	// Ignore key points too close to the border of an image
	keypoint_border_width = 0;
	// Ignore key points too close to the bottom of an image
	keypoint_border_bottom = 0;

	// Number of features matches for an image pair to be considered a match
	min_num_feat_matches = 16;
	// Ray angle threshold
	ray_angle_threshold = 2.0;

	// Homography RANSAC params
	homography_rounds = 256;
	homography_threshold = 6.0;

	// F-matrix RANSAC params
	fmatrix_rounds = 2048;
	fmatrix_threshold = 9.0;

	// RANSAC threshold for estimating projection matrix
	projection_estimation_threshold = 4.0; 

	min_proj_error_threshold = 8.0;
	max_proj_error_threshold = 16.0;

}
