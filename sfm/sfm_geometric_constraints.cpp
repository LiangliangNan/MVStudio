#include "sfm.h"
#include "epipolar.h"
#include "register.h"

#include "../basic/logger.h"
#include "../math/matrix_driver.h"

#include <cassert>


namespace sfm {

	void SfM::compute_geometric_constraints(int new_image_start) {

		compute_epipolar_geometry(true, new_image_start);
		compute_transforms(false, new_image_start);

		make_match_lists_symmetric();

		compute_tracks(new_image_start);
	}


	/* Compute epipolar geometry between all matching images */
	void SfM::compute_epipolar_geometry(bool remove_bad_matches, int new_image_start) {
		transforms_.clear();

		std::vector<MatchIndex> remove;

		int num_images = num_of_images();
		for (int i = 0; i < num_images; i++) {
			MatchTable::MatchAdjList::iterator iter;
			for (iter = match_table_.begin(i); iter != match_table_.end(i); iter++) {
				unsigned int j = iter->index; // first;

				assert(images_match(i, j));

				MatchIndex idx(i, j);
				MatchIndex idx_rev(j, i);

				bool connect12 = compute_epipolar_geometry(i, j, remove_bad_matches);
				if (!connect12) {
					if (remove_bad_matches) {
						remove.push_back(idx);
						remove.push_back(idx_rev);

						transforms_.erase(idx);
						transforms_.erase(idx_rev);
					}
				}
				else {
					matrix_transpose(3, 3,
						transforms_[idx].fmatrix,
						transforms_[idx_rev].fmatrix);
				}
			}
		}

		int num_removed = (int)remove.size();
		for (int i = 0; i < num_removed; i++) {
			int img1 = remove[i].first;
			int img2 = remove[i].second;
			remove_match(img1, img2);
		}
	}


	/* Compute epipolar geometry between a given pair of images */
	bool SfM::compute_epipolar_geometry(int idx1, int idx2, bool remove_bad_matches) {
		if (idx1 == idx2) {
			Logger::warn(title()) << "computing epipolar geometry for identical images (ignored)" << std::endl;
			return false;
		}

		assert(image_data_[idx1].keys_loaded);
		assert(image_data_[idx2].keys_loaded);

		MatchIndex offset(idx1, idx2);
		MatchIndex offset_rev(idx2, idx1);
		std::vector<KeypointMatch> &list = match_table_.match_list(offset);

		double F[9];
		std::vector<int> inliers = estimate_f_Matrix(
			image_data_[idx1].keys,	
			image_data_[idx2].keys,
			list,
			option_.fmatrix_rounds,
			option_.fmatrix_threshold,
			F
			);

		int num_inliers = (int)inliers.size();
#ifdef _DEBUG
		Logger::out(title()) << "inliers [" << idx1 << ", " << idx2 << "] = " << num_inliers << " out of " << list.size() << std::endl;
#endif // _DEBUG
		if (remove_bad_matches) {
			/* Refine the matches */
			std::vector<KeypointMatch> new_match_list;

			for (int i = 0; i < num_inliers; i++) {
				new_match_list.push_back(list[inliers[i]]);
			}

			list.clear();
			list = new_match_list;
		}

		if (num_inliers >= option_.min_num_feat_matches) {
			transforms_[offset] = Transform();
			transforms_[offset_rev] = Transform();

			memcpy(transforms_[offset].fmatrix, F, 9 * sizeof(double));
#ifdef _DEBUG
			Logger::out(title()) << "inliers [" << idx1 << ", " << idx2 << "] = " << num_inliers << " out of " << list.size() << std::endl;
#endif // _DEBUG
			return true;
		}
		else {
			return false;
		}
	}


	/* Compute transforms between all matching images */
	void SfM::compute_transforms(bool remove_bad_matches, int new_image_start) {
		transforms_.clear();

		int num_images = num_of_images();
		for (int i = 0; i < num_images; i++) {
			MatchTable::MatchAdjList::iterator iter;
			for (iter = match_table_.begin(i); iter != match_table_.end(i); iter++) {
				unsigned int j = iter->index;

				assert(images_match(i, j));

				MatchIndex idx(i, j);
				MatchIndex idx_rev(j, i);

				transforms_[idx] = Transform();
				transforms_[idx_rev] = Transform();

				bool connect12 = compute_transform(i, j, remove_bad_matches);

				if (!connect12) {
					if (remove_bad_matches) {
						match_table_.remove_match(idx);
						match_table_.remove_match(idx_rev);

						transforms_.erase(idx);
						transforms_.erase(idx_rev);
					}
				}
				else {
					matrix_invert(3, transforms_[idx].H, transforms_[idx_rev].H);
				}
			}
		}

#ifdef _DEBUG
		/* Print the inlier ratios */
		FILE *f = fopen("pairwise_scores2.txt", "w");

		for (int i = 0; i < num_images; i++) {
			MatchTable::MatchAdjList::iterator iter;

			for (iter = match_table_.begin(i); iter != match_table_.end(i); iter++) {
				unsigned int j = iter->index; // first;
				assert(images_match(i, j));

				MatchIndex idx(i, j);
				fprintf(f, "%d %d %0.5f\n", i, j, transforms_[idx].inlier_ratio);
			}
		}
		fclose(f);
#endif // _DEBUG
	}


	/* Compute a transform between a given pair of images */
	bool SfM::compute_transform(int idx1, int idx2, bool remove_bad_matches) {
		if (idx1 == idx2) {
			Logger::warn(title()) << "computing transform for identical images (ignored)" << std::endl;
			return false;
		}

		assert(image_data_[idx1].keys_loaded);
		assert(image_data_[idx2].keys_loaded);

		MatchIndex offset(idx1, idx2);
		std::vector<KeypointMatch> &list = match_table_.match_list(offset);

		double M[9];
		std::vector<int> inliers = estimate_transform(
			image_data_[idx1].keys,
			image_data_[idx2].keys,
			list, MotionHomography,
			option_.homography_rounds,
			option_.homography_threshold,
			M
			);

		int num_inliers = (int)inliers.size();
#ifdef _DEBUG
		Logger::out(title()) << "inliers [" << idx1 << ", " << idx2 << "] = " << num_inliers << " out of " << list.size() << std::endl;
#endif // _DEBUG
		if (remove_bad_matches) {
			/* Refine the matches */
			std::vector<KeypointMatch> new_match_list;

			for (int i = 0; i < num_inliers; i++) {
				new_match_list.push_back(list[inliers[i]]);
			}

			list.clear();
			list = new_match_list;
		}

#define MIN_INLIERS 10
		if (num_inliers >= MIN_INLIERS) {
			transforms_[offset].num_inliers = num_inliers;
			transforms_[offset].inlier_ratio = ((double)num_inliers) / ((double)list.size());

			memcpy(transforms_[offset].H, M, 9 * sizeof(double));
			return true;
		}
		else {
			return false;
		}
	}



}