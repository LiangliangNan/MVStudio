#ifndef _SFM_OPTION_H_
#define _SFM_OPTION_H_


#include <string>


struct SfmOption
{
	SfmOption();

	std::string		output_directory;
	std::string		image_directory;
	std::string		key_directory;

	std::string		list_file;
	std::string		match_table_file;
	std::string		bundle_output_file;
	std::string		bundle_output_base;


	// Should we fix the focal length for all cameras (set to init_focal_length)?
	bool		fixed_focal_length;		
	double		init_focal_length;

	// Should we use camera constraints?
	bool		use_camera_constraints;  
	std::string	camera_constraint_file;

	// Should we add soft constraint on focal lengths to stay near their estimated values
	bool		constrain_focal;		
	// The weight for the focal length constraint
	double		constrain_focal_weight;

	// Should we estimate radial distortion parameters (2 coefficients) for each camera?
	bool		estimate_distortion;

	/* Point constraint fields */
	bool		use_point_constraints;
	double		point_constraint_weight;
	std::string	point_constraint_file;

	// Skip full optimization stages
	bool		skip_full_bundle;

	// Ignore key points too close to the border of an image
	int		keypoint_border_width; 
	// Ignore key points too close to the bottom of an image
	int		keypoint_border_bottom; 

	// Number of features matches for an image pair to be considered a match
	// (Minimum number of matches needed to register an image)
	int		min_num_feat_matches;   

	// Ray angle threshold
	double	ray_angle_threshold;    

	// Homography RANSAC params
	int		homography_rounds;  
	double	homography_threshold;

	// F-matrix RANSAC params
	int		fmatrix_rounds;        
	double	fmatrix_threshold;

	// RANSAC threshold for estimating projection matrix
	double	projection_estimation_threshold;

	double	min_proj_error_threshold;
	double	max_proj_error_threshold;

};


#endif
