#include "sfm.h"
#include <easy3d/util/logging.h>


using namespace easy3d;
namespace sfm {



	void SfM::remove_all_matches() {
		match_table_.remove_all();
	}

	void SfM::set_match(int img1, int img2) {
		match_table_.set_match(MatchIndex(img1, img2));
	}

	void SfM::remove_match(int img1, int img2) {
		match_table_.remove_match(MatchIndex(img1, img2));
	}

	bool SfM::images_match(int img1, int img2) {
		return match_table_.contains(MatchIndex(img1, img2));
	}

	void SfM::prune_multiple_matches() {
		unsigned int num_images = num_of_images();

		for (unsigned int i = 0; i < num_images; i++) {
			std::vector<unsigned int> remove;

			MatchTable::MatchAdjList::iterator iter;
			for (iter = match_table_.begin(i); iter != match_table_.end(i); iter++) {
				HashSetInt seen;

				int num_pruned = 0;
				std::vector<KeypointMatch> &list = iter->match_list;

				/* Unmark keys */
				int num_matches = (int)list.size();
				for (int k = 0; k < num_matches; k++) {
					int idx2 = list[k].key_idx2;

					// if (GetKey(j,idx2).m_extra != -1) {
					if (seen.find(idx2) != seen.end()) {
						/* This is a repeat */
						// printf("[%d] Pruning repeat %d\n", i, idx2);
						list.erase(list.begin() + k);
						num_matches--;
						k--;

						num_pruned++;
					}
					else {
						/* Mark this key as matched */
						seen.insert(idx2);
					}
				}

				unsigned int j = iter->index; // first;
				//LOG(INFO) << "Pruned [" << i << ", " << j << "] = " << num_pruned << " / " << num_matches + num_pruned << std::endl;

				if (num_matches < option_.min_num_feat_matches) {
					/* Get rid of... */
					remove.push_back(iter->index); // first);
				}
			}

			for (unsigned int j = 0; j < remove.size(); j++) {
				int idx2 = remove[j];
				match_table_.remove_match(MatchIndex(i, idx2));
				LOG(INFO) << "Removing [" << i << ", " << idx2 << "]" << std::endl;
			}
		}
	}

	void SfM::remove_border_matches() {
		int num = num_of_images();
		for (int i = 0; i < num; i++) {
			for (int j = i + 1; j < num; j++) {
 				if (!images_match(i, j))
 					continue;
 
				remove_border_matches(i, j, option_.keypoint_border_width);
 			}
 		}
	}

	void SfM::remove_bottom_matches() {
		int num = num_of_images();
		for (int i = 0; i < num; i++) {
			for (int j = i + 1; j < num; j++) {
				if (!images_match(i, j))
					continue;

				remove_bottom_matches(i, j, option_.keypoint_border_bottom);
			}
		}
	}


	/* Remove matches close to the edges of the two given images */
	void SfM::remove_border_matches(int i1, int i2, int border_width) {
		MatchIndex idx(i1, i2);

		assert(match_table_.contains(idx));
		assert(image_data_[i1].keys_loaded);
		assert(image_data_[i2].keys_loaded);

		std::vector<KeypointMatch>& list = match_table_.match_list(idx);
		int num_matches = (int)list.size();

		int w1 = image_data_[i1].width();
		int h1 = image_data_[i1].height();
		int w2 = image_data_[i2].width();
		int h2 = image_data_[i2].height();

		double w1_min = -0.5 * w1 + border_width;
		double w1_max = 0.5 * w1 - border_width;
		double h1_min = -0.5 * h1 + border_width;
		double h1_max = 0.5 * h1 - border_width;

		double w2_min = -0.5 * w2 + border_width;
		double w2_max = 0.5 * w2 - border_width;
		double h2_min = -0.5 * h2 + border_width;
		double h2_max = 0.5 * h2 - border_width;

		int num_removed = 0;
		for (int k = 0; k < num_matches; k++) {
			KeypointMatch &m = list[k];

			const Keypoint &k1 = image_data_[i1].keys[m.key_idx1];
			const Keypoint &k2 = image_data_[i2].keys[m.key_idx2];

			if (k1.x < w1_min || k1.x > w1_max ||
				k1.y < h1_min || k1.y > h1_max ||
				k2.x < w2_min || k2.x > w2_max ||
				k2.y < h2_min || k2.y > h2_max) {

				/* Erase this match */
				list.erase(list.begin() + k);
				k--;
				num_matches--;

				num_removed++;
			}
		}

		LOG(INFO) << "Removed " << num_removed << " matches from pair (" << i1 << ", " << i2 << ")" << std::endl;
	}

	/* Remove matches close to the bottom edge of the two given images */
	void SfM::remove_bottom_matches(int i1, int i2, int border_width)
	{
		MatchIndex idx(i1, i2);

		// assert(m_match_lists.find(idx) != m_match_lists.end());
		assert(match_table_.contains(idx));
		assert(image_data_[i1].keys_loaded);
		assert(image_data_[i2].keys_loaded);

		std::vector<KeypointMatch> &list = match_table_.match_list(idx);
		// m_match_lists[idx];
		int num_matches = (int)list.size();

		int h1 = image_data_[i1].height();
		int h2 = image_data_[i2].height();

		double h1_min = -0.5 * h1 + border_width;
		double h2_min = -0.5 * h2 + border_width;

		int num_removed = 0;
		for (int k = 0; k < num_matches; k++) {
			KeypointMatch &m = list[k];

			const Keypoint &k1 = image_data_[i1].keys[m.key_idx1];
			const Keypoint &k2 = image_data_[i2].keys[m.key_idx2];

			if (k1.y < h1_min || k2.y < h2_min) {
				/* Erase this match */
				list.erase(list.begin() + k);
				k--;
				num_matches--;
				num_removed++;
			}
		}

		LOG(INFO) << "Removed " << num_removed << " matches from pair (" << i1 << ", " << i2 << ")" << std::endl;
	}


	/* Make match lists symmetric */
	void SfM::make_match_lists_symmetric()
	{
		unsigned int num_images = num_of_images();

		std::vector<MatchIndex> matches;
		for (unsigned int i = 0; i < num_images; i++) {
			MatchTable::MatchAdjList::const_iterator iter;
			for (iter = match_table_.begin(i); iter != match_table_.end(i); iter++) {
				unsigned int j = iter->index; // iter->second;

				if (j <= i)
					continue;

				assert(images_match(i, j));

				MatchIndex idx(i, j);
				MatchIndex idx_rev(j, i);

				const std::vector<KeypointMatch> &list = iter->match_list;
				int num_matches = (int)list.size();

				match_table_.set_match(idx_rev);
				match_table_.clear_match(idx_rev);

				for (int k = 0; k < num_matches; k++) {
					KeypointMatch m1, m2;

					m1 = list[k];

					m2.key_idx1 = m1.key_idx2;
					m2.key_idx2 = m1.key_idx1;

					match_table_.add_match(idx_rev, m2);
				}

				matches.push_back(idx);
			}
		}

		int num_matches = (int)matches.size();
		for (int i = 0; i < num_matches; i++) {
			unsigned int img1 = matches[i].first;
			unsigned int img2 = matches[i].second;
			set_match(img2, img1);
		}

		matches.clear();
	}


	/* Use the bundle-adjusted points to create a new set of matches */
	void SfM::set_matches_from_points(int threshold) {
		/* Clear all matches */
		remove_all_matches();

		int num_points = (int)point_data_.size();
		for (int i = 0; i < num_points; i++) {
			int num_views = (int)point_data_[i].views.size();
			if (num_views < threshold)
				continue;

			for (int j = 0; j < num_views; j++) {
				for (int k = 0; k < num_views; k++) {
					if (j == k) continue;

					ImageKey view1 = point_data_[i].views[j];
					ImageKey view2 = point_data_[i].views[k];

					int v1 = view1.first;
					int v2 = view2.first;

					int k1 = view1.second;
					int k2 = view2.second;

					KeypointMatch m;
					m.key_idx1 = k1;
					m.key_idx2 = k2;

					set_match(v1, v2);
					MatchIndex idx(v1, v2);
					match_table_.add_match(idx, m);
				}
			}
		}
	}


}