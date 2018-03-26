#include "sfm.h"
#include "../basic/logger.h"

#include <queue>
#include <cassert>


namespace sfm {

	static bool CompareFirst(const KeypointMatch &k1, const KeypointMatch &k2) {
		return (k1.key_idx1 < k2.key_idx1);
	}


	/* Compute a set of tracks that explain the matches */
	void SfM::compute_tracks(int new_image_start)
	{
		int num_images = num_of_images();

		/* Clear all marks for new images */
		for (int i = 0; i < num_images; i++) {
			/* If this image has no neighbors, don't worry about its keys */
			int num_nbrs = (int)match_table_.num_of_neighbors(i);
			if (num_nbrs == 0)
				continue;

			int num_features = image_data_[i].num_of_keys();
			image_data_[i].key_flags.resize(num_features);
		}

		for (int i = 0; i < num_images; i++) {
			MatchTable::MatchAdjList::iterator iter;
			for (iter = match_table_.begin(i); iter != match_table_.end(i); iter++) {
				// MatchIndex idx = *iter;
				std::vector<KeypointMatch> &list = iter->match_list; // iter->second; // m_match_lists[idx];
				sort(list.begin(), list.end(), CompareFirst);
			}
		}

		bool *img_marked = new bool[num_images];
		memset(img_marked, 0, num_images * sizeof(bool));

		std::vector<int> touched;
		touched.reserve(num_images);

		int pt_idx = 0;
		std::vector<TrackData> tracks;
		for (int i = 0; i < num_images; i++) {
			int num_features = image_data_[i].num_of_keys();

			/* If this image has no neighbors, skip it */
			int num_nbrs = (int)match_table_.num_of_neighbors(i);
			if (num_nbrs == 0)
				continue;

			for (int j = 0; j < num_features; j++) {
				ImageKeyVector features;
				std::queue<ImageKey> features_queue;

				/* Check if this feature was visited */
				if (image_data_[i].key_flags[j])
					continue; // already visited this feature

				/* Reset flags */
				int num_touched = (int)touched.size();
				for (int k = 0; k < num_touched; k++)
					img_marked[touched[k]] = false;
				touched.clear();

				/* Do a breadth first search given this feature */
				image_data_[i].key_flags[j] = true;

				features.push_back(ImageKey(i, j));
				features_queue.push(ImageKey(i, j));

				img_marked[i] = true;
				touched.push_back(i);

				int num_rounds = 0;
				while (!features_queue.empty()) {
					num_rounds++;

					ImageKey feature = features_queue.front();
					features_queue.pop();

					int img1 = feature.first;
					int f1 = feature.second;
					KeypointMatch dummy;
					dummy.key_idx1 = f1;

					int start_idx;
					/* Limit new images to point only to other new images */
					if (img1 >= new_image_start) {
						start_idx = new_image_start;
					}
					else {
						start_idx = 0;
					}

					MatchTable::MatchAdjList &nbrs = match_table_.neighbors(img1);

					MatchTable::MatchAdjList::iterator iter;
					for (iter = nbrs.begin(); iter != nbrs.end(); iter++) {
						unsigned int k = iter->index; 

						if (img_marked[k])
							continue;

						MatchIndex base(img1, k);
						std::vector<KeypointMatch> &list = match_table_.match_list(base);

						/* Do a binary search for the feature */
						std::pair<std::vector<KeypointMatch>::iterator,
							std::vector<KeypointMatch>::iterator> p;

						p = equal_range(list.begin(), list.end(),
							dummy, CompareFirst);

						if (p.first == p.second)
							continue;  /* not found */

						assert((p.first)->key_idx1 == f1);
						int idx2 = (p.first)->key_idx2;

						/* Check if we visited this point already */
						assert(idx2 < image_data_[k].num_of_keys());

						if (image_data_[k].key_flags[idx2])
							continue;

						/* Mark and push the point */
						image_data_[k].key_flags[idx2] = true;
						features.push_back(ImageKey(k, idx2));
						features_queue.push(ImageKey(k, idx2));

						img_marked[k] = true;
						touched.push_back(k);
					}
				} /* while loop */

				if (features.size() >= 2) {
					//Logger::out(title()) << "point with " << features.size() << " projections found" << std::endl;
					tracks.push_back(TrackData(features));
					pt_idx++;
				}
				else {
					// printf("Feature only has %d points (%d inconsistent)\n", 
					//       (int) features.size(), num_inconsistent);
				}

			} /* for loop over features */
		} /* for loop over images */

		Logger::out(title()) << "found " << pt_idx << " points" << std::endl;

		if (pt_idx != (int)tracks.size()) {
			Logger::err(title()) << "point count inconsistent" << std::endl;
		}

		/* Clear match lists */
		Logger::out(title()) << "clearing match lists..." << std::endl;


		/* Create the new consistent match lists */
		Logger::out(title()) << "creating consistent match lists..." << std::endl;

		int num_pts = pt_idx;
		for (int i = 0; i < num_pts; i++) {
			int num_features = (int)tracks[i].views.size();

			for (int j = 0; j < num_features; j++) {
				int img1 = tracks[i].views[j].first;
				int key1 = tracks[i].views[j].second;

				image_data_[img1].visible_points.push_back(i);
				image_data_[img1].visible_keys.push_back(key1);
			}
		}

		/* Save the tracks */
		track_data_ = tracks;

		Logger::out(title()) << "done" << std::endl;

		//////////////////////////////////////////////////////////////////////////

		/* Set match flags */
		remove_all_matches();

		int num_tracks = (int)track_data_.size();
		for (int i = 0; i < num_tracks; i++) {
			TrackData &track = track_data_[i];
			int num_views = (int)track.views.size();

			for (int j = 0; j < num_views; j++) {
				int img1 = track.views[j].first;
				assert(img1 >= 0 && img1 < num_images);

				for (int k = j + 1; k < num_views; k++) {
					int img2 = track.views[k].first;

					assert(img2 >= 0 && img2 < num_images);

					set_match(img1, img2);
					set_match(img2, img1);
				}
			}
		}
	}


	void SfM::set_tracks(int image) {
		Logger::out(title()) << "setting tracks for image " << image << std::endl;

		ImageData &img_data = image_data_[image];
		assert(img_data.keys_loaded);

		int num_tracks = (int)img_data.visible_points.size();

		for (int i = 0; i < num_tracks; i++) {
			int tr = img_data.visible_points[i];
			int key = img_data.visible_keys[i];

			assert(key < (int)img_data.keys.size());

			img_data.keys[key].track = tr;
		}
	}




	/* Return the intersection of two int vectors */
	static std::vector<int> vector_intersection(const std::vector<int> &v1, const std::vector<int> &v2)
	{
#ifndef WIN32
		__gnu_cxx::hash_set<int> seen;
#else
		std::unordered_set<int> seen;
#endif

		int v1_size = (int)v1.size();
		int v2_size = (int)v2.size();

		std::vector<int> intersection;

		for (int i = 0; i < v1_size; i++)
			seen.insert(v1[i]);

		for (int i = 0; i < v2_size; i++) {
			if (seen.find(v2[i]) != seen.end())
				intersection.push_back(v2[i]);
		}

		seen.clear();

		return intersection;
	}


	void SfM::set_matches_from_tracks(int img1, int img2) {
		std::vector<int> &tracks1 = image_data_[img1].visible_points;
		std::vector<int> &tracks2 = image_data_[img2].visible_points;

		std::vector<int> isect = vector_intersection(tracks1, tracks2);

		int num_isect = (int)isect.size();

		if (num_isect == 0)
			return;

		MatchIndex idx(img1, img2);

		std::vector<KeypointMatch> &matches = match_table_.match_list(idx);
		// m_match_lists[idx];

		matches.clear();
		matches.resize(num_isect);

		for (int i = 0; i < num_isect; i++) {
			int tr = isect[i];

			std::pair < std::vector<int>::const_iterator,
				std::vector<int>::const_iterator > p;
			const std::vector<int> &pt1 = image_data_[img1].visible_points;
			p = equal_range(pt1.begin(), pt1.end(), tr);
			assert(p.first != p.second);
			int offset = int(p.first - pt1.begin());
			int k1 = image_data_[img1].visible_keys[offset];

			const std::vector<int> &pt2 = image_data_[img2].visible_points;
			p = equal_range(pt2.begin(), pt2.end(), tr);
			assert(p.first != p.second);
			offset = int(p.first - pt2.begin());
			int k2 = image_data_[img2].visible_keys[offset];

			matches[i] = KeypointMatch(k1, k2);
		}
	}


	int SfM::get_num_track_matches(int img1, int img2)
	{
		const std::vector<int> &tracks1 = image_data_[img1].visible_points;
		const std::vector<int> &tracks2 = image_data_[img2].visible_points;

		// std::vector<int> isect = GetVectorIntersection(tracks1, tracks2);
		// int num_isect = (int) isect.size();

		std::vector<int>::const_iterator iter;
		for (iter = tracks2.begin(); iter != tracks2.end(); iter++) {
			int track_idx = *iter;
			track_data_[track_idx].extra = 0;
		}

		for (iter = tracks1.begin(); iter != tracks1.end(); iter++) {
			int track_idx = *iter;
			track_data_[track_idx].extra = 1;
		}

		int num_isect = 0;
		for (iter = tracks2.begin(); iter != tracks2.end(); iter++) {
			int track_idx = *iter;
			num_isect += track_data_[track_idx].extra;
		}

		return num_isect;
	}



}