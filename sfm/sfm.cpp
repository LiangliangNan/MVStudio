#include "sfm.h"
#include "../basic/logger.h"


namespace sfm {


	SfM::SfM(SfmOption opt)
		: option_(opt)
		, point_constraints_(nil)
	{
		initial_pair_[0] = -1;
		initial_pair_[1] = -1;
	}

	SfM::~SfM()
	{
		if (point_constraints_)
			delete point_constraints_;
	}

	void SfM::apply(PointSet* pset) {
		if (option_.fixed_focal_length && option_.estimate_distortion) {
			Logger::err(title()) << "\'fixed_focal_length\' and \'estimate_distortion\' are incompatible" << std::endl;
			return;
		}

		Logger::out(title()) << "loading image names..." << std::endl;
		load_image_names();
		if (image_data_.empty()) {
			return;
		}
		int num_images = num_of_images();
		Logger::out(title()) << "done. " << num_images << " images." << std::endl;

		if (option_.use_camera_constraints)
			read_camera_constraints(option_.camera_constraint_file);

// 		if (has_geometric_constraints) {
// 			read_geometric_constraints(filename);
// 			return;
// 		}
// 		else 
 		{
			Logger::out(title()) << "loading match table..." << std::endl;
			load_match_table();
			prune_multiple_matches();

			Logger::out(title()) << "loading keys..." << std::endl;
			load_keys();

			if (option_.keypoint_border_width > 0)
				remove_border_matches();
			if (option_.keypoint_border_bottom > 0)
				remove_bottom_matches();

			compute_geometric_constraints();		// Compute initial image information

// 			save the geometric constraints into file, so you don't need to compute it anymore.
// 			write_geometric_constraints(filename);
		}

		run_sfm(pset);

		//if (m_bundle_version < 0.3)
		//	FixReflectionBug();

	}

	/* Return the number of images */
	int  SfM::num_of_images() const { 
		return (int)image_data_.size(); 
	}

	/* Return the number of matches for image pair (img1, img2) */
	int  SfM::num_of_matches(int img1, int img2) const {
		int img_min = std::min(img1, img2);
		int img_max = std::max(img1, img2);

		MatchIndex idx(img_min, img_max);
		return match_table_.num_of_matches(idx);
	}

	Keypoint &SfM::get_key(int img, int key) {
		return image_data_[img].keys[key];
	}


	int SfM::num_of_keys(int img) {
		return (int)image_data_[img].keys.size();
	}



}