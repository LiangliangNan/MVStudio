#include "sfm.h"
#include "geometry.h"
#include "epipolar.h"
#include "bundle_adjustment.h"
#include "../basic/logger.h"
#include "../basic/stop_watch.h"
#include "../math/matrix_driver.h"
#include "../mvglib/qsort.h"
#include "../mvglib/triangulate.h"
#include "../pointset/point_set.h"
#include "../basic/file_utils.h"
#include "../basic/progress.h"

#include <cassert>
#include <iomanip>



namespace sfm {

	void SfM::run_sfm(PointSet* pset) {
		StopWatch w;

		/* Set track pointers to -1 */
		for (int i = 0; i < (int)track_data_.size(); i++) {
			track_data_[i].extra = -1;
		}

		/* For now, assume all images form one connected component */
		int num_images = num_of_images();
		int *added_order = new int[num_images];
		int *added_order_inv = new int[num_images];

		/* **** Run bundle adjustment! **** */

		camera_params_t *cameras = new camera_params_t[num_images];
		int max_pts = (int)track_data_.size(); // 1243742; /* HACK! */
		vec3d *points = new vec3d[max_pts];
		vec3i *colors = new vec3i[max_pts];
		std::vector<ImageKeyVector> pt_views;


		/* Initialize the bundle adjustment */
		int num_init_cams = 0;
		initialize_bundle_adjust(num_init_cams, added_order, added_order_inv,
			cameras, points, colors, pt_views);

		int i_best = -1, j_best = -1, max_matches = 0;
		double max_score = 0.0;
		int curr_num_cameras, curr_num_pts;
		int pt_count;

		if (num_init_cams == 0) {
			bundle_pick_initial_pair(i_best, j_best, true);

			added_order[0] = i_best;
			added_order[1] = j_best;

			Logger::out(title()) << "Adjusting cameras " << i_best << " and " << j_best << " (score = " << max_score << ")" << std::endl;

			/* **** Set up the initial cameras **** */
			double init_focal_length_0 = 0.0, init_focal_length_1 = 0.0;
			pt_count = curr_num_pts =
				setup_initial_camera_pair(i_best, j_best,
				init_focal_length_0, init_focal_length_1,
				cameras, points, colors, pt_views);
#ifdef _DEBUG
			std::string tmp_file = option_.output_directory + "/bundle.init.out";
			dump_output_file(tmp_file, num_images, 2, curr_num_pts, added_order, cameras, points, colors, pt_views);
#endif // _DEBUG

			/* Run sfm for the first time */
			double error0;
			error0 = run_bundle_adjustment(curr_num_pts, 2, 0, false, cameras, points, added_order, colors, pt_views);
			Logger::out(title()) << "focal lengths: " << cameras[0].f << ", " << cameras[1].f << std::endl;

			if (pset) {
				std::vector<vec3f>& pts = pset->points();	pts.resize(curr_num_pts);
				std::vector<vec3f>& cls = pset->colors();	cls.resize(curr_num_pts);
				for (std::size_t i = 0; i < curr_num_pts; ++i) {
					pts[i] = vec3f(points[i].x, points[i].y, points[i].z);
					cls[i] = vec3f(colors[i].x / 255.0f, colors[i].y / 255.0f, colors[i].z / 255.0f);
				}
				pset->fit();
				pset->update_all();
			}
#ifdef _DEBUG	
			tmp_file = option_.output_directory + "/points_0001.ply";
			dump_points_to_ply(tmp_file, curr_num_pts, 2, points, colors, cameras);

			std::ostringstream string_stream;
			string_stream << option_.output_directory << "/" << option_.bundle_output_base << "0001.out";
			tmp_file = string_stream.str();
			dump_output_file(tmp_file, num_images, 2, curr_num_pts, added_order, cameras, points, colors, pt_views);
#endif // _DEBUG

			curr_num_cameras = 2;
		}
		else {
			curr_num_cameras = num_init_cams;
			pt_count = curr_num_pts = (int)point_data_.size();
		}

		int round = 0;
		ProgressLogger progress(num_images);
		while (curr_num_cameras < num_images) {
			if (progress.is_canceled())
				break;
			progress.notify(curr_num_cameras);

			int parent_idx;
			int max_cam =
				find_camera_with_most_matches(curr_num_cameras, curr_num_pts,
				added_order, parent_idx,
				max_matches, pt_views);

			printf("[SifterApp::BundleAdjust] max_matches = %d\n", max_matches);

			if (max_matches < option_.min_num_feat_matches)
				break; /* No more connections */

			/* Find all images with 90% of the matches of the maximum */
			std::vector<ImagePair> image_set;

			if (false && max_matches < 48) {
				image_set.push_back(ImagePair(max_cam, parent_idx));
			}
			else {
				// int nMatches = MIN(100, iround(0.75 * max_matches));
				int nMatches = iround(0.75 * max_matches);
				image_set =
					find_cameras_with_n_matches(nMatches,
					curr_num_cameras, curr_num_pts,
					added_order, pt_views);
			}

			int num_added_images = (int)image_set.size();

			printf("[SifterApp::BundleAdjustFast] Registering %d images\n",
				num_added_images);

			for (int i = 0; i < num_added_images; i++)
				printf("[SifterApp::BundleAdjustFast] Adjusting camera %d\n",
				image_set[i].first);

			/* Now, throw the new cameras into the mix */
			int image_count = 0;
			for (int i = 0; i < num_added_images; i++) {
				int next_idx = image_set[i].first;
				int parent_idx = image_set[i].second;

				added_order[curr_num_cameras + image_count] = next_idx;

				printf("[SifterApp::BundleAdjust[%d]] Adjusting camera %d "
					"(parent = %d)\n",
					round, next_idx,
					(parent_idx == -1 ? -1 : added_order[parent_idx]));

				/* **** Set up the new camera **** */
				bool success = false;
				camera_params_t camera_new =
					bundle_initialize_image(
					image_data_[next_idx],
					next_idx, curr_num_cameras + image_count,
					curr_num_cameras, curr_num_pts,
					added_order, points,
					NULL /*cameras + parent_idx*/, cameras,
					pt_views, &success);

				if (success) {
					cameras[curr_num_cameras + image_count] = camera_new;
					image_count++;
				}
				else {
					printf("[BundleAdjust] Couldn't initialize image %d\n",
						next_idx);
					image_data_[next_idx].ignore_in_bundle = true;
				}
			}


			/* Compute the distance between the first pair of cameras */
#if 0
			double dist0 = camera_distance(...);
#else
			double dist0 = 0.0;
#endif

			printf("[SifterApp::BundleAdjust] Adding new matches\n");

			pt_count = curr_num_pts;

			curr_num_cameras += image_count;

			pt_count =
				bundle_adjust_add_all_new_points(pt_count, curr_num_cameras,
				added_order, cameras,
				points, colors,
				dist0, pt_views);

			curr_num_pts = pt_count;

			printf("[SifterApp::BundleAdjust] Number of points = %d\n", pt_count);
			fflush(stdout);

			if (!option_.skip_full_bundle) {
				/* Run sfm again to update parameters */
				run_bundle_adjustment(curr_num_pts, curr_num_cameras, 0, false,
					cameras, points, added_order, colors, pt_views);

				/* Remove bad points and cameras */
				remove_bad_points_and_cameras(curr_num_pts, curr_num_cameras + 1,
					added_order, cameras, points, colors,
					pt_views);

				//printf("  focal lengths:\n");

				for (int i = 0; i < curr_num_cameras; i++) {
					Logger::out(title()) << "focal length for image [" << i << "](" << FileUtils::simple_name(image_data_[added_order[i]].image_file) << "): "
						<< cameras[i].f << " (" << image_data_[added_order[i]].init_focal << ")" << std::endl;
// 					printf("   [%03d] %0.3f (%0.3f) %s %d; %0.3e %0.3e\n",
// 						i, cameras[i].f,
// 						image_data_[added_order[i]].init_focal,
// 						image_data_[added_order[i]].image_file.c_str(),
// 						added_order[i], cameras[i].k[0], cameras[i].k[1]);
				}
			}

			round++;

			if (pset) {
				std::vector<vec3f>& pts = pset->points();	pts.resize(curr_num_pts);
				std::vector<vec3f>& cls = pset->colors();	cls.resize(curr_num_pts);
				for (std::size_t i = 0; i < curr_num_pts; ++i) {
					pts[i] = vec3f(points[i].x, points[i].y, points[i].z);
					cls[i] = vec3f(colors[i].x / 255.0f, colors[i].y / 255.0f, colors[i].z / 255.0f);
				}				
				pset->update_all();
			}

#ifdef _DEBUG
			/* Dump output for this round */
			std::ostringstream string_stream;
			string_stream << option_.output_directory << "/points_" << std::setfill('0') << std::setw(4) << curr_num_cameras << ".ply";
			std::string tmp_file = string_stream.str();
			dump_points_to_ply(tmp_file, curr_num_pts, curr_num_cameras, points, colors, cameras);

			string_stream.str("");
			string_stream.clear();
			string_stream << option_.output_directory << "/" << option_.bundle_output_base << std::setfill('0') << std::setw(4) << curr_num_cameras << ".out";
			tmp_file = string_stream.str();
			dump_output_file(tmp_file, num_images, curr_num_cameras, curr_num_pts, added_order, cameras, points, colors, pt_views);
#endif // _DEBUG

		}

		Logger::out(title()) << "bundle adjustment took " << w.elapsed() << " seconds" << std::endl;

		/* Dump output */
		dump_output_file(option_.bundle_output_file,
			num_images, curr_num_cameras, curr_num_pts,
			added_order, cameras, points, colors, pt_views);

		/* Save the camera parameters and points */

		/* Cameras */
		for (int i = 0; i < num_images; i++) {
			image_data_[i].camera.adjusted = false;
		}

		for (int i = 0; i < curr_num_cameras; i++) {
			int img = added_order[i];

			image_data_[img].camera.adjusted = true;
			memcpy(image_data_[img].camera.R, cameras[i].R,
				9 * sizeof(double));

			matrix_product(3, 3, 3, 1,
				cameras[i].R, cameras[i].t,
				image_data_[img].camera.t);

			matrix_scale(3, 1,
				image_data_[img].camera.t, -1.0,
				image_data_[img].camera.t);

			image_data_[img].camera.focal = cameras[i].f;

			image_data_[img].camera.Finalize();
		}

		/* Points */
		for (int i = 0; i < curr_num_pts; i++) {
			/* Check if the point is visible in any view */
			if ((int)pt_views[i].size() == 0)
				continue; /* Invisible */

			PointData pdata;
			pdata.pos = points[i];
			pdata.color = colors[i];

			for (int j = 0; j < (int)pt_views[i].size(); j++) {
				int v = pt_views[i][j].first;
				int vnew = added_order[v];
				pdata.views.push_back(ImageKey(vnew, pt_views[i][j].second));
			}

			point_data_.push_back(pdata);
		}

		delete[] added_order;
		delete[] added_order_inv;
		delete[] points;
		delete[] colors;

		set_matches_from_points();

		// 		for (int i = 0; i < num_images; i++) {
		// 			if (image_data_[i].camera.adjusted)
		// 				image_data_[i] was used during bundle adjustment
		// 			else
		// 				image_data_[i] was NOT used during bundle adjustment
		// 		}
	}



	static int compare_doubles(const void *d1, const void *d2)
	{
		double a = *(double *)d1;
		double b = *(double *)d2;

		if (a < b) return -1;
		if (a > b) return 1;
		return 0;
	}



	double SfM::run_bundle_adjustment(int num_pts, int num_cameras, int start_camera,
		bool fix_points, camera_params_t *init_camera_params,
		vec3d *init_pts, int *added_order, vec3i *colors,
		std::vector<ImageKeyVector> &pt_views, double eps2,
		double *S, double *U, double *V,
		double *W, bool remove_outliers)
	{
#define MIN_POINTS 20
		int num_outliers = 0;
		int total_outliers = 0;
		double dist_total = 0.0;
		int num_dists = 0;

		int *remap = new int[num_pts];
		vec3d *nz_pts = new vec3d[num_pts];

		do {
			if (num_pts - total_outliers < MIN_POINTS) {
				printf("[RunSFM] Too few points remaining, exiting!\n");
				fflush(stdout);

				dist_total = DBL_MAX;
				break;
			}

			/* Set up the vmask and projections */
			char *vmask = NULL;
			double *projections = NULL;

			int num_projections = 0;
			for (int i = 0; i < num_pts; i++) {
				num_projections += (int)pt_views[i].size();
			}

			vmask = new char[num_pts * num_cameras];
			projections = new double[2 * num_projections];

			for (int i = 0; i < num_pts * num_cameras; i++)
				vmask[i] = 0;

			int arr_idx = 0;
			int nz_count = 0;
			for (int i = 0; i < num_pts; i++) {
				int num_views = (int)pt_views[i].size();

				if (num_views > 0) {
					for (int j = 0; j < num_views; j++) {
						int c = pt_views[i][j].first;
						int v = added_order[c];
						int k = pt_views[i][j].second;

						vmask[nz_count * num_cameras + c] = 1;

						projections[2 * arr_idx + 0] = get_key(v, k).x;
						projections[2 * arr_idx + 1] = get_key(v, k).y;

						arr_idx++;
					}

					remap[i] = nz_count;
					nz_pts[nz_count] = init_pts[i];
					nz_count++;
				}
				else {
					remap[i] = -1;
				}
			}

			dist_total = 0.0;
			num_dists = 0;

			bool fixed_focal = option_.fixed_focal_length;

			StopWatch w;
			bundle_adjustment(nz_count, num_cameras, start_camera, vmask, projections,
				fixed_focal ? 0 : 1, 0,
				option_.estimate_distortion ? 1 : 0, 1,
				init_camera_params, nz_pts,
				(option_.constrain_focal ? 1 : 0),
				(option_.use_point_constraints) ? 1 : 0,
				point_constraints_, option_.point_constraint_weight,
				fix_points ? 1 : 0, eps2, V, S, U, W);

			Logger::out(title()) << "run_bundle_adjustment took " << w.elapsed() << " seconds" << std::endl;

			/* Check for outliers */

			std::vector<int> outliers;
			std::vector<double> reproj_errors;

			for (int i = 0; i < num_cameras; i++) {
				ImageData &data = image_data_[added_order[i]];

				double K[9] = { init_camera_params[i].f, 0.0, 0.0,
					0.0, init_camera_params[i].f, 0.0,
					0.0, 0.0, 1.0 };
				// double w[3] = { 0.0, 0.0, 0.0 };
				double dt[3] = { init_camera_params[i].t[0],
					init_camera_params[i].t[1],
					init_camera_params[i].t[2] };

				/* Compute inverse distortion parameters */
				if (option_.estimate_distortion) {
					double *k = init_camera_params[i].k;
					double k_dist[6] = { 0.0, 1.0, 0.0, k[0], 0.0, k[1] };
					double w_2 = 0.5 * data.width();
					double h_2 = 0.5 * data.height();
					double max_radius =
						sqrt(w_2 * w_2 + h_2 * h_2) / init_camera_params[i].f;

					invert_distortion(6, 6, 0.0, max_radius,
						k_dist, init_camera_params[i].k_inv);
				}

				if (data.known_intrinsics) {
					double *k = init_camera_params[i].k_known;
					double k_dist[8] =
					{ 0.0, 1.0, 0.0, k[0], 0.0, k[1], 0.0, k[4] };
					double w_2 = 0.5 * data.width();
					double h_2 = 0.5 * data.height();
					double max_radius =
						sqrt(w_2 * w_2 + h_2 * h_2) /
						init_camera_params[i].K_known[0];

					invert_distortion(8, 6, 0.0, max_radius, k_dist,
						init_camera_params[i].k_inv);
				}

				int num_keys = num_of_keys(added_order[i]);

				int num_pts_proj = 0;
				for (int j = 0; j < num_keys; j++) {
					if (get_key(added_order[i], j).extra >= 0) {
						num_pts_proj++;
					}
				}

				double *dists = new double[num_pts_proj];
				int pt_count = 0;

				std::vector<Keypoint>::iterator iter;
				// for (int j = 0; j < num_keys; j++) {
				for (iter = image_data_[added_order[i]].keys.begin();
					iter != image_data_[added_order[i]].keys.end();
					iter++) {

					const Keypoint &key = *iter;

					if (key.extra >= 0) {
						double b[3], pr[2];
						double dx, dy, dist;
						int pt_idx = key.extra;

						b[0] = nz_pts[remap[pt_idx]].x;
						b[1] = nz_pts[remap[pt_idx]].y;
						b[2] = nz_pts[remap[pt_idx]].z;

						sfm_project_rd(&(init_camera_params[i]), K,
							init_camera_params[i].k,
							init_camera_params[i].R, dt, b, pr,
							option_.estimate_distortion, true);

						dx = pr[0] - key.x;
						dy = pr[1] - key.y;

						dist = sqrt(dx * dx + dy * dy);
						dist_total += dist;
						num_dists++;

						dists[pt_count] = dist;

						pt_count++;
					}
				}

				/* Estimate the median of the distances */
				double med = kth_element_copy(num_pts_proj,
					iround(0.8 /* 0.9 */ * num_pts_proj),
					dists);

				median_copy(num_pts_proj, dists);

#define NUM_STDDEV 2.0 // 3.0 // 6.0
				double thresh = 1.2 * NUM_STDDEV * med; /* k * stddev */
				ogf_clamp(thresh, option_.min_proj_error_threshold, option_.max_proj_error_threshold);

				/* Compute the average reprojection error for this
				* camera */

				double sum = 0.0;
				for (int j = 0; j < num_pts_proj; j++) {
					sum += dists[j];
				}

				double avg = sum / num_pts_proj;
				// 			printf("[RunSFM] Mean error cam %d[%d] [%d pts]: %0.3e "
				// 				"[med: %0.3e, %0.3e]\n",
				// 				i, added_order[i], num_pts_proj, avg,
				// 				kth_element_copy(num_pts_proj,
				// 				iround(0.5 * num_pts_proj), dists),
				// 				thresh);

				// printf("Outlier threshold is %0.3f\n", thresh);

				pt_count = 0;
				for (int j = 0; j < num_keys; j++) {
					int pt_idx = get_key(added_order[i], j).extra;

					if (pt_idx < 0)
						continue;

					/* Don't remove constrained points */
					if (option_.use_point_constraints && point_constraints_[pt_idx].x != 0.0) {
						pt_count++;
						continue;
					}

					if (dists[pt_count] > thresh) {
						/* Remove this point from consideration */
						bool found = false;
						for (int k = 0; k < (int)outliers.size(); k++) {
							if (outliers[k] == pt_idx) {
								found = true;
								break;
							}
						}

						if (!found) {
							outliers.push_back(pt_idx);
							reproj_errors.push_back(dists[pt_count]);
						}
					}
					pt_count++;
				}

#define OUTPUT_VERBOSE_STATS
#ifdef OUTPUT_VERBOSE_STATS
#define NUM_ERROR_BINS 10
				qsort(dists, num_pts_proj, sizeof(double), compare_doubles);

				double pr_min = dists[0];
				double pr_max = dists[num_pts_proj - 1];
				double pr_step = (pr_max - pr_min) / NUM_ERROR_BINS;

				/* Break histogram into 10 bins */
				int idx_count = 0;
				for (int i = 0; i < NUM_ERROR_BINS; i++) {
					double max = pr_min + (i + 1) * pr_step;
					int start = idx_count;

					while (idx_count < num_pts_proj && dists[idx_count] <= max)
						idx_count++;

					int bin_size = idx_count - start;
					// 				printf("   E[%0.3e--%0.3e]: %d [%0.3f]\n",
					// 					max - pr_step, max, bin_size,
					// 					bin_size / (double)num_pts_proj);
				}
#endif

				delete[] dists;
			}

			/* Remove outlying points */
			if (remove_outliers) {
				for (int i = 0; i < (int)outliers.size(); i++) {
					int idx = outliers[i];

					// 					printf("[RunSFM] Removing outlier %d "
					// 						"(reproj error: %0.3f)\n", idx, reproj_errors[i]);

					if (colors != NULL) {
						colors[idx] = vec3i(0x0, 0x0, 0xff);
					}

					int num_views = (int)pt_views[idx].size();

					for (int j = 0; j < num_views; j++) {
						int v = pt_views[idx][j].first;
						int k = pt_views[idx][j].second;

						vmask[idx * num_cameras + v] = 0;

						/* Sanity check */
						if (get_key(added_order[v], k).extra != idx)
							printf("Error!  Entry for (%d,%d) "
							"should be %d, but is %d\n",
							added_order[v], k,
							idx, get_key(added_order[v], k).extra);

						get_key(added_order[v], k).extra = -2;
					}

					pt_views[idx].clear();
				}

				num_outliers = (int)outliers.size();
				total_outliers += num_outliers;

				Logger::out(title()) << "removing " << num_outliers << " outliers" << std::endl;
			}

			delete[] vmask;
			delete[] projections;

			for (int i = 0; i < num_pts; i++) {
				if (remap[i] != -1) {
					init_pts[i] = nz_pts[remap[i]];
				}
			}

			if (!remove_outliers) break;

		} while (num_outliers > 0);

		delete[] remap;
		delete[] nz_pts;

		return dist_total / num_dists;
	}



	/* Initialize the bundle adjustment procedure (loading an existing
	* model if one exists) */
	void SfM::initialize_bundle_adjust(int &num_init_cams,
		int *added_order,
		int *added_order_inv,
		camera_params_t *cameras,
		vec3d *points, vec3i *colors,
		std::vector<ImageKeyVector> &pt_views)
	{
		int num_images = num_of_images();

		/* Initialize all keypoints to have not been matched */
		for (int i = 0; i < num_images; i++) {
			std::vector<Keypoint>::iterator iter;
			for (iter = image_data_[i].keys.begin();
				iter != image_data_[i].keys.end();
				iter++) {

				iter->extra = -1;
			}
		}

		/* Initialize the bundle adjustment with the existing model (if
		* there is one) */
		/* Cameras */
		num_init_cams = 0;

		for (int i = 0; i < num_images; i++) {
			if (added_order_inv != NULL)
				added_order_inv[i] = -1;
		}

		/* Points */
		int n = 0;
		double error = 0.0;
		for (int i = 0; i < (int)point_data_.size(); i++) {
			points[i] = vec3d(point_data_[i].pos[0],
				point_data_[i].pos[1],
				point_data_[i].pos[2]);
			// -point_data_[i].pos[2]);

			colors[i] = vec3i(point_data_[i].color[0],
				point_data_[i].color[1],
				point_data_[i].color[2]);

			ImageKeyVector views;
			int num_views = (int)point_data_[i].views.size();
			double *views_arr = new double[num_views];
			int *perm = new int[num_views];

			for (int j = 0; j < num_views; j++) {
				ImageKey p = point_data_[i].views[j];

				if (added_order_inv[p.first] == -1) {
					printf("Error: added_order_inv[%d] doesn't exist!\n", p.first);
				}

				if (p.first < 0 || p.first >= num_images) {
					printf("Error: p.first[%d] out of range\n", p.first);
				}

				if (p.second < 0 ||
					p.second >= (int)image_data_[p.first].keys.size()) {

					printf("Error: p.second[%d,%d] out of range\n", p.second,
						(int)image_data_[p.first].keys.size());
				}

				int image_idx = p.first;
				int key_idx = p.second;

				image_data_[p.first].keys[p.second].extra = i;
				p.first = added_order_inv[p.first];

				views_arr[j] = (double)p.first;
				views.push_back(p);

				int track = image_data_[image_idx].keys[key_idx].track;

				if (track != -1) { // assert(track != -1);
					track_data_[track].extra = i;
				}
			}

			/* Sort the views */
			qsort_ascending();
			qsort_perm(num_views, views_arr, perm);

			ImageKeyVector views_sorted;
			for (int j = 0; j < num_views; j++) {
				views_sorted.push_back(views[perm[j]]);
			}

			pt_views.push_back(views_sorted);

			delete[] views_arr;
			delete[] perm;

			// pt_views.push_back(views);
		}

		printf("  Avg. proj error [%d projections] = %0.3e\n", n,
			sqrt(error / n));
	}

	/* Find the camera with the most matches to existing points */
	int SfM::find_camera_with_most_matches(int num_cameras, int num_points,
		int *added_order,
		int &parent_idx, int &max_matches,
		const std::vector<ImageKeyVector> &pt_views)
	{
		max_matches = 0;

		int i_best = -1;
		double top_score = 0.0;

		parent_idx = -1;

		int num_images = num_of_images();
		for (int i = 0; i < num_images; i++) {
			if (image_data_[i].ignore_in_bundle)
				continue;

			/* Check if we added this image already */
			bool added = false;
			for (int j = 0; j < num_cameras; j++) {
				if (added_order[j] == i) {
					added = true;
					break;
				}
			}

			if (added)
				continue;

			int num_existing_matches = 0;
			int parent_idx_best = -1;

			/* Find the tracks seen by this image */
			const std::vector<int> &tracks = image_data_[i].visible_points;
			int num_tracks = (int)tracks.size();

			for (int j = 0; j < num_tracks; j++) {
				int tr = tracks[j];
				if (track_data_[tr].extra < 0)
					continue;

				/* This tracks corresponds to a point */
				int pt = track_data_[tr].extra;
				if ((int)pt_views[pt].size() == 0)
					continue;

				num_existing_matches++;
			}

			if (num_existing_matches > 0)
				printf("  existing_matches[%d] = %d\n", i, num_existing_matches);

			double score = num_existing_matches;

			if (score > 0.0)
				printf("  score[%d]   = %0.3f\n", i, score);

			if (score > top_score) {
				i_best = i;
				parent_idx = parent_idx_best;
				max_matches = num_existing_matches;
				top_score = score;
			}

			// delete [] saw;
		}

		if (parent_idx == -1) {
			printf("Error: parent not found\n");
		}

		return i_best;
	}

	/* Find all cameras with at least N matches to existing points */
	std::vector<ImagePair> SfM::find_cameras_with_n_matches(int n,
		int num_cameras,
		int num_points,
		int *added_order,
		const std::vector<ImageKeyVector> &pt_views)
	{
		std::vector<ImagePair> image_pairs;

		int num_images = num_of_images();
		for (int i = 0; i < num_images; i++) {
			if (image_data_[i].ignore_in_bundle)
				continue;

			/* Check if we added this image already */
			bool added = false;
			for (int j = 0; j < num_cameras; j++) {
				if (added_order[j] == i) {
					added = true;
					break;
				}
			}

			if (added)
				continue;

			int num_existing_matches = 0;
			int parent_idx_best = -1;

			/* Find the tracks seen by this image */
			const std::vector<int> &tracks = image_data_[i].visible_points;
			int num_tracks = (int)tracks.size();

			for (int j = 0; j < num_tracks; j++) {
				int tr = tracks[j];
				if (track_data_[tr].extra < 0)
					continue;

				/* This tracks corresponds to a point */
				int pt = track_data_[tr].extra;
				if ((int)pt_views[pt].size() == 0)
					continue;

				num_existing_matches++;
			}

			if (num_existing_matches >= n)
				image_pairs.push_back(ImagePair(i, parent_idx_best));

			// delete [] saw;
		}

		return image_pairs;
	}


	/* Triangulate a subtrack */
	vec3d SfM::triangulate_n_views(const ImageKeyVector &views,
		int *added_order, camera_params_t *cameras,
		double &error, bool explicit_camera_centers)
	{
		int num_views = (int)views.size();

		vec2d *pv = new vec2d[num_views];
		double *Rs = new double[9 * num_views];
		double *ts = new double[3 * num_views];

		for (int i = 0; i < num_views; i++) {
			camera_params_t *cam = NULL;

			int camera_idx = views[i].first;
			int image_idx = added_order[camera_idx];
			int key_idx = views[i].second;
			Keypoint &key = get_key(image_idx, key_idx);

			double p3[3] = { key.x, key.y, 1.0 };

			double K[9], Kinv[9];
			get_intrinsics(cameras[camera_idx], K);
			matrix_invert(3, K, Kinv);

			double p_n[3];
			matrix_product(3, 3, 3, 1, Kinv, p3, p_n);

			// EDIT!!!
			pv[i] = vec2d(-p_n[0], -p_n[1]);
			pv[i] = undistort_normalized_point(pv[i], cameras[camera_idx]);

			cam = cameras + camera_idx;

			memcpy(Rs + 9 * i, cam->R, 9 * sizeof(double));
			if (!explicit_camera_centers) {
				memcpy(ts + 3 * i, cam->t, 3 * sizeof(double));
			}
			else {
				matrix_product(3, 3, 3, 1, cam->R, cam->t, ts + 3 * i);
				matrix_scale(3, 1, ts + 3 * i, -1.0, ts + 3 * i);
			}
		}

		vec3d pt = triangulate_n(num_views, pv, Rs, ts, &error);

		error = 0.0;
		for (int i = 0; i < num_views; i++) {
			int camera_idx = views[i].first;
			int image_idx = added_order[camera_idx];
			int key_idx = views[i].second;
			Keypoint &key = get_key(image_idx, key_idx);

			vec2d pr = sfm_project_final(cameras + camera_idx, pt,
				explicit_camera_centers ? 1 : 0,
				option_.estimate_distortion ? 1 : 0);

			double dx = pr.x - key.x;
			double dy = pr.y - key.y;

			error += dx * dx + dy * dy;
		}

		error = sqrt(error / num_views);

		delete[] pv;
		delete[] Rs;
		delete[] ts;

		return pt;
	}


	/* Add new points to the bundle adjustment */
	int SfM::bundle_adjust_add_all_new_points(int num_points, int num_cameras,
		int *added_order,
		camera_params_t *cameras,
		vec3d *points, vec3i *colors,
		double reference_baseline,
		std::vector<ImageKeyVector> &pt_views,
		double max_reprojection_error,
		int min_views)
	{
		std::vector<int> track_idxs;
		std::vector<ImageKeyVector> new_tracks;

		// __gnu_cxx::hash_map<int,bool> tracks_seen;
		int num_tracks_total = (int)track_data_.size();
		int *tracks_seen = new int[num_tracks_total];
		for (int i = 0; i < num_tracks_total; i++) {
			tracks_seen[i] = -1;
		}

		/* Gather up the projections of all the new tracks */
		for (int i = 0; i < num_cameras; i++) {
			int image_idx1 = added_order[i];

			int num_keys = num_of_keys(image_idx1);

			for (int j = 0; j < num_keys; j++) {
				Keypoint &key = get_key(image_idx1, j);

				if (key.track == -1)
					continue;  /* Key belongs to no track */

				if (key.extra != -1)
					continue;  /* Key is outlier or has already been added */

				int track_idx = key.track;

				/* Check if this track is already associated with a point */
				if (track_data_[track_idx].extra != -1)
					continue;

				/* Check if we've seen this track */
				int seen = tracks_seen[track_idx];

				if (seen == -1) {
					/* We haven't yet seen this track, create a new track */
					tracks_seen[track_idx] = (int)new_tracks.size();

					ImageKeyVector track;
					track.push_back(ImageKey(i, j));
					new_tracks.push_back(track);
					track_idxs.push_back(track_idx);
				}
				else {
					new_tracks[seen].push_back(ImageKey(i, j));
				}
			}
		}

		delete[] tracks_seen;

		/* Now for each (sub) track, triangulate to see if the track is
		* consistent */
		int pt_count = num_points;

		int num_ill_conditioned = 0;
		int num_high_reprojection = 0;
		int num_cheirality_failed = 0;
		int num_added = 0;

		int num_tracks = (int)new_tracks.size();
		for (int i = 0; i < num_tracks; i++) {
			int num_views = (int)new_tracks[i].size();

			if (num_views < min_views) continue;  /* Not enough views */


			/* Check if at least two cameras fix the position of the point */
			bool conditioned = false;
			bool good_distance = false;
			double max_angle = 0.0;
			for (int j = 0; j < num_views; j++) {
				for (int k = j + 1; k < num_views; k++) {
					int camera_idx1 = new_tracks[i][j].first;
					int image_idx1 = added_order[camera_idx1];
					int key_idx1 = new_tracks[i][j].second;

					int camera_idx2 = new_tracks[i][k].first;
					int image_idx2 = added_order[camera_idx2];
					int key_idx2 = new_tracks[i][k].second;

					Keypoint &key1 = get_key(image_idx1, key_idx1);
					Keypoint &key2 = get_key(image_idx2, key_idx2);

					vec2d p(key1.x, key1.y);
					vec2d q(key2.x, key2.y);

					double angle = compute_ray_angle(p, q,
						cameras[camera_idx1],
						cameras[camera_idx2]);

					if (angle > max_angle)
						max_angle = angle;

					/* Check that the angle between the rays is large enough */
					if (RAD2DEG(angle) >= option_.ray_angle_threshold) {
						conditioned = true;
					}

					good_distance = true;
				}
			}

			if (!conditioned || !good_distance) {
				num_ill_conditioned++;
				continue;
			}

			double error;
			vec3d pt = triangulate_n_views(new_tracks[i], added_order, cameras,
				error, true);

			if (isnan(error) || error > max_reprojection_error) {
				num_high_reprojection++;
				continue;
			}

			bool all_in_front = true;
			for (int j = 0; j < num_views; j++) {
				int camera_idx = new_tracks[i][j].first;
				bool in_front = check_cheirality(pt, cameras[camera_idx]);

				if (!in_front) {
					all_in_front = false;
					break;
				}
			}

			if (!all_in_front) {
				num_cheirality_failed++;

				continue;
			}

			/* All tests succeeded, so let's add the point */

			fflush(stdout);

			points[pt_count] = pt;

			int camera_idx = new_tracks[i][0].first;
			int image_idx = added_order[camera_idx];
			int key_idx = new_tracks[i][0].second;

			unsigned char r = get_key(image_idx, key_idx).r;
			unsigned char g = get_key(image_idx, key_idx).g;
			unsigned char b = get_key(image_idx, key_idx).b;
			colors[pt_count] = vec3i((int)r, (int)g, (int)b);

			pt_views.push_back(new_tracks[i]);

			/* Set the point index on the keys */
			for (int j = 0; j < num_views; j++) {
				int camera_idx = new_tracks[i][j].first;
				int image_idx = added_order[camera_idx];
				int key_idx = new_tracks[i][j].second;
				get_key(image_idx, key_idx).extra = pt_count;
			}

			int track_idx = track_idxs[i];
			track_data_[track_idx].extra = pt_count;

			pt_count++;
			num_added++;
		}

		printf("[AddAllNewPoints] Added %d new points\n", num_added);
		printf("[AddAllNewPoints] Ill-conditioned tracks: %d\n",
			num_ill_conditioned);
		printf("[AddAllNewPoints] Bad reprojections: %d\n", num_high_reprojection);
		printf("[AddAllNewPoints] Failed cheirality checks: %d\n",
			num_cheirality_failed);

		return pt_count;
	}


	int SfM::remove_bad_points_and_cameras(int num_points, int num_cameras,
		int *added_order,
		camera_params_t *cameras,
		vec3d *points, vec3i *colors,
		std::vector<ImageKeyVector> &pt_views)
	{
		int num_pruned = 0;

		for (int i = 0; i < num_points; i++) {
			double *pos = points[i].data();
			int num_views = (int)pt_views[i].size();

			if (num_views == 0)
				continue;

			double max_angle = 0.0;
			for (int j = 0; j < num_views; j++) {
				int v1 = pt_views[i][j].first;

				double r1[3];
				matrix_diff(3, 1, 3, 1, pos, cameras[v1].t, r1);
				double norm = matrix_norm(3, 1, r1);
				matrix_scale(3, 1, r1, 1.0 / norm, r1);

				for (int k = j + 1; k < num_views; k++) {
					int v2 = pt_views[i][k].first;

					double r2[3];
					matrix_diff(3, 1, 3, 1, pos, cameras[v2].t, r2);
					double norm = matrix_norm(3, 1, r2);
					matrix_scale(3, 1, r2, 1.0 / norm, r2);

					double dot;
					matrix_product(1, 3, 3, 1, r1, r2, &dot);

					ogf_clamp(dot, -1.0 + 1.0e-8, 1.0 - 1.0e-8);
					double angle = acos(dot);

					if (angle > max_angle) {
						max_angle = angle;
					}
				}
			}

			if (RAD2DEG(max_angle) < 0.5 * option_.ray_angle_threshold) {
				printf("[RemoveBadPointsAndCamera] "
					"Removing point %d with angle %0.3f\n",
					i, RAD2DEG(max_angle));

				for (int j = 0; j < num_views; j++) {
					// Set extra flag back to 0
					int v = pt_views[i][j].first;
					int k = pt_views[i][j].second;
					get_key(added_order[v], k).extra = -1;
				}

				pt_views[i].clear();

				if (colors != NULL) {
					colors[i] = vec3i(0x0, 0x0, 0xff);
				}

				num_pruned++;
			}
		}

		printf("[RemoveBadPointsAndCameras] Pruned %d points\n", num_pruned);

		return num_pruned;
	}


	/* Pick a good initial pair of cameras to bootstrap the bundle
	* adjustment */
	void SfM::bundle_pick_initial_pair(int &i_best, int &j_best,
		bool use_init_focal_only)
	{
		/* Compute the match matrix */
		int num_images = num_of_images();
		int max_matches = 0;
		double max_score = 0.0;

		int max_matches_2 = 0;
		double max_score_2 = 0.0;

		i_best = j_best = -1;
		int i_best_2 = -1;
		int j_best_2 = -1;

		int i_best_3 = -1;
		int j_best_3 = -1;

		if (initial_pair_[0] != -1 && initial_pair_[1] != -1) {
			i_best = initial_pair_[0];
			j_best = initial_pair_[1];

			printf("[BundleAdjust] Setting initial pair to "
				"%d and %d\n", i_best, j_best);

			return;
		}

		double SCORE_THRESHOLD = 2.0; // 1.0 // 2.0 // 1.0

		/* Compute score for each image pair */
		int max_pts = 0;
		for (int i = 0; i < num_images; i++) {
			if (image_data_[i].ignore_in_bundle)
				continue;

			for (int j = i + 1; j < num_images; j++) {
				if (image_data_[j].ignore_in_bundle)
					continue;

				MatchIndex idx(i, j);

				int num_matches = get_num_track_matches(i, j);
				max_pts += num_matches;

#define MATCH_THRESHOLD 32
#define MIN_SCORE 1.0e-1
#define MIN_MATCHES 80

				if (num_matches <= MATCH_THRESHOLD) {
					continue;
				}

				double score = 0.0;
				double ratio = transforms_[idx].inlier_ratio;
				if (ratio == 0.0) {
					score = MIN_SCORE;
				}
				else {
					score = 1.0 / transforms_[idx].inlier_ratio;
				}

				/* Compute the primary score */
				if (num_matches > max_matches && score > SCORE_THRESHOLD) {
					max_matches = num_matches;
					max_score = score;

					i_best = i;
					j_best = j;
				}

				/* Compute the backup score */
				if (num_matches > MIN_MATCHES && score > max_score_2) {
					max_matches_2 = num_matches;
					max_score_2 = score;
					i_best_2 = i;
					j_best_2 = j;
				}
			}
		}

		/* Set track pointers to -1 (GetNumTrackMatches alters these
		* values) */
		for (int i = 0; i < (int)track_data_.size(); i++) {
			track_data_[i].extra = -1;
		}

		if (i_best == -1 && j_best == -1) {
			if (i_best_2 == -1 && j_best_2 == -1) {
				printf("[BundleAdjust] Error: no good camera pairs found!\n");

				if (use_init_focal_only) {
					printf("[BundleAdjust] Trying a backup approach...\n");
					bundle_pick_initial_pair(i_best, j_best, false);
				}
				else {
					printf("[BundleAdjust] Picking first two cameras...\n");

					i_best = 0;
					j_best = 1;
				}
			}
			else {
				i_best = i_best_2;
				j_best = j_best_2;
			}
		}
	}


	/* Setup the initial camera pair for bundle adjustment */
	int SfM::setup_initial_camera_pair(int i_best, int j_best,
		double &init_focal_length_0,
		double &init_focal_length_1,
		camera_params_t *cameras,
		vec3d *points, vec3i *colors,
		std::vector<ImageKeyVector> &pt_views)
	{
		/* Load the keys for the images */
		image_data_[i_best].load_keys();
		image_data_[j_best].load_keys();
		image_data_[i_best].load_key_colors();
		image_data_[j_best].load_key_colors();

		set_matches_from_tracks(i_best, j_best);

		set_tracks(i_best);
		set_tracks(j_best);

		initialize_camera_parameters(image_data_[i_best], cameras[0]);
		initialize_camera_parameters(image_data_[j_best], cameras[1]);
		set_camera_constraints(i_best, cameras + 0);
		set_camera_constraints(j_best, cameras + 1);

		/* Put first camera at origin */
		cameras[0].R[0] = 1.0;  cameras[0].R[1] = 0.0;  cameras[0].R[2] = 0.0;
		cameras[0].R[3] = 0.0;  cameras[0].R[4] = 1.0;  cameras[0].R[5] = 0.0;
		cameras[0].R[6] = 0.0;  cameras[0].R[7] = 0.0;  cameras[0].R[8] = 1.0;

		/* Initialize the positions of the cameras (using constraints,
		* if provided) */
		if (image_data_[i_best].camera.constrained[0])
			cameras[0].t[0] = image_data_[i_best].camera.constraints[0];
		else
			cameras[0].t[0] = 0.0;

		if (image_data_[i_best].camera.constrained[1])
			cameras[0].t[1] = image_data_[i_best].camera.constraints[1];
		else
			cameras[0].t[1] = 0.0;

		if (image_data_[i_best].camera.constrained[2])
			cameras[0].t[2] = image_data_[i_best].camera.constraints[2];
		else
			cameras[0].t[2] = 0.0;


		init_focal_length_0 = cameras[0].f = image_data_[i_best].init_focal;
		init_focal_length_1 = cameras[1].f = image_data_[j_best].init_focal;

		bool solved_for_extrinsics = false;
		/* Solve for the initial locations */
		if (estimate_relative_pose2(i_best, j_best, cameras[0], cameras[1])) {
			solved_for_extrinsics = true;
		}

		else {
#define INITIAL_DEPTH 3.0 // 1000.0 // 3.0 /* was 3.0, fixme! */

			/* Put second camera at origin too */
			cameras[1].R[0] = 1.0;  cameras[1].R[1] = 0.0;  cameras[1].R[2] = 0.0;
			cameras[1].R[3] = 0.0;  cameras[1].R[4] = 1.0;  cameras[1].R[5] = 0.0;
			cameras[1].R[6] = 0.0;  cameras[1].R[7] = 0.0;  cameras[1].R[8] = 1.0;

			if (image_data_[j_best].camera.constrained[0])
				cameras[1].t[0] = image_data_[j_best].camera.constraints[0];
			else
				cameras[1].t[0] = 0.0;

			if (image_data_[j_best].camera.constrained[1])
				cameras[1].t[1] = image_data_[j_best].camera.constraints[1];
			else
				cameras[1].t[1] = 0.0;

			if (image_data_[j_best].camera.constrained[2])
				cameras[1].t[2] = image_data_[j_best].camera.constraints[2];
			else
				cameras[1].t[2] = 0.0;
		}

		if (option_.constrain_focal) {
			set_focal_constraint(image_data_[i_best], cameras + 0);
			set_focal_constraint(image_data_[j_best], cameras + 1);
		}

		/* **** Set up the initial 3D points **** */
		printf("[BundleAdjust] Adding initial matches...\n");

		int pt_count = 0;
		MatchIndex list_idx(i_best, j_best);
		std::vector<KeypointMatch> &list = match_table_.match_list(list_idx);

		int num_matches = (int)list.size();
		for (int i = 0; i < num_matches; i++) {
			int key_idx1 = list[i].key_idx1;
			int key_idx2 = list[i].key_idx2;

			double x_proj = get_key(i_best, key_idx1).x;
			double y_proj = get_key(i_best, key_idx1).y;

			/* Back project the point to a constant depth */
			if (!solved_for_extrinsics) {
				double x_pt = (x_proj / option_.init_focal_length) * INITIAL_DEPTH;
				double y_pt = (y_proj / option_.init_focal_length) * INITIAL_DEPTH;
				double z_pt = INITIAL_DEPTH + cameras[0].t[2];

				points[pt_count] = vec3d(x_pt, y_pt, z_pt);
			}
			else {
				double x_proj1 = get_key(i_best, key_idx1).x;
				double y_proj1 = get_key(i_best, key_idx1).y;
				double x_proj2 = get_key(j_best, key_idx2).x;
				double y_proj2 = get_key(j_best, key_idx2).y;

				double error;

				vec2d p(x_proj1, y_proj1);
				vec2d q(x_proj2, y_proj2);

				bool in_front = true;
				double angle = 0.0;
				points[pt_count] = triangulate(p, q, cameras[0], cameras[1],
					error, in_front, angle, true);

				// 			printf(" tri.error[%d] = %0.3f\n", i, error);

				if (error > option_.projection_estimation_threshold) {
					printf(" skipping point\n");
					continue;
				}
			}

			/* Get the color of the point */
			unsigned char r = get_key(i_best, key_idx1).r;
			unsigned char g = get_key(i_best, key_idx1).g;
			unsigned char b = get_key(i_best, key_idx1).b;
			colors[pt_count] = vec3i((int)r, (int)g, (int)b);

			get_key(i_best, key_idx1).extra = pt_count;
			get_key(j_best, key_idx2).extra = pt_count;

			int track_idx = get_key(i_best, key_idx1).track;
			track_data_[track_idx].extra = pt_count;

			ImageKeyVector views;
			views.push_back(ImageKey(0, key_idx1));
			views.push_back(ImageKey(1, key_idx2));
			pt_views.push_back(views);

			pt_count++;
		}

		// m_match_lists[list_idx].clear();

		match_table_.clear_match(list_idx);

		return pt_count;
	}



	camera_params_t SfM::bundle_initialize_image(ImageData &data,
		int image_idx, int camera_idx,
		int num_cameras, int num_points,
		int *added_order,
		vec3d *points,
		camera_params_t *parent,
		camera_params_t *cameras,
		std::vector<ImageKeyVector> &pt_views,
		bool *success_out,
		bool refine_cameras_and_points)
	{
		if (success_out != NULL)
			*success_out = true;

		/* Load the keys */
		data.load_keys();
		set_tracks(image_idx);

		/* **** Connect the new camera to any existing points **** */
		int num_pts_solve = 0;
		int num_keys = (int)data.keys.size();
		vec3d *points_solve = new vec3d[num_keys];
		vec2d *projs_solve = new vec2d[num_keys];
		vec2d *projs_solve_orig = new vec2d[num_keys];
		int *idxs_solve = new int[num_keys];
		int *keys_solve = new int[num_keys];

		printf("[BundleInitializeImage] "
			"Connecting existing matches...\n");

		/* Find the tracks seen by this image */
		std::vector<int> &tracks = data.visible_points;
		int num_tracks = (int)tracks.size();

		for (int i = 0; i < num_tracks; i++) {
			int tr = tracks[i];
			if (track_data_[tr].extra < 0)
				continue;

			/* This tracks corresponds to a point */
			int pt = track_data_[tr].extra;
			if ((int)pt_views[pt].size() == 0)
				continue;

			int key = data.visible_keys[i];

			// printf("  Connecting existing point [%d] (cam: %d)\n", 
			//        pt, image_idx);

			/* Add the point to the set we'll use to solve for
			* the camera position */
			points_solve[num_pts_solve] = points[pt];

			projs_solve[num_pts_solve] = vec2d(data.keys[key].x, data.keys[key].y);

			idxs_solve[num_pts_solve] = pt;
			keys_solve[num_pts_solve] = key;

			num_pts_solve++;
		}

		if (num_pts_solve < option_.min_num_feat_matches) {
			printf("[BundleInitializeImage] Couldn't initialize\n");

			if (success_out != NULL)
				*success_out = false;

			camera_params_t dummy;

			delete[] points_solve;
			delete[] projs_solve;
			delete[] projs_solve_orig;
			delete[] idxs_solve;
			delete[] keys_solve;

			image_data_[image_idx].unload_keys();

			return dummy;
		}

		/* **** Solve for the camera position **** */
		printf("[BundleInitializeImage] Initializing camera...\n");
		fflush(stdout);
		double Kinit[9], Rinit[9], tinit[3];
		std::vector<int> inliers, inliers_weak, outliers;
		bool success =
			find_and_verify_camera(num_pts_solve, points_solve, projs_solve,
			idxs_solve, Kinit, Rinit, tinit,
			option_.projection_estimation_threshold,
			16.0 * option_.projection_estimation_threshold, /*4.0*/
			inliers, inliers_weak, outliers);

		if (!success) {
			printf("[BundleInitializeImage] Couldn't initialize\n");

			if (success_out != NULL)
				*success_out = false;

			camera_params_t dummy;

			delete[] points_solve;
			delete[] projs_solve;
			delete[] projs_solve_orig;
			delete[] idxs_solve;
			delete[] keys_solve;

			image_data_[image_idx].unload_keys();

			return dummy;
		}

		camera_params_t camera_new;
		initialize_camera_parameters(data, camera_new);

		/* Start with the new camera at same place as the best match */
		if (success) {
			/* Set up the new camera */
			memcpy(camera_new.R, Rinit, 9 * sizeof(double));

			matrix_transpose_product(3, 3, 3, 1, Rinit, tinit, camera_new.t);
			matrix_scale(3, 1, camera_new.t, -1.0, camera_new.t);

			/* Set up the new focal length */
			set_camera_constraints(added_order[num_cameras], &camera_new);

			if (option_.fixed_focal_length) {
				camera_new.f = option_.init_focal_length;
			}
			else {
				double ratio;
				double init = data.init_focal;
				double obs = 0.5 * (Kinit[0] + Kinit[4]);

				printf("[BundleInitializeImage] "
					"Camera has initial focal length of %0.3f\n", init);

				if (init > obs) ratio = init / obs;
				else            ratio = obs / init;

				if (ratio < 1.4) {
					camera_new.f = data.init_focal;
					if (option_.constrain_focal)
						set_focal_constraint(image_data_[image_idx],
						&camera_new);
				}
				else {
					printf("[BundleInitializeImage] "
						"Estimated focal length of %0.3f "
						"is too different\n", obs);
					camera_new.f = 0.5 * (Kinit[0] + Kinit[4]);
				}
			}
		}
		else {
			printf("[BundleInitializeImage] Error!  "
				"Pose estimation failed!\n");
		}

		/* **** Finally, start the bundle adjustment **** */
		printf("[BundleInitializeImage] Adjusting...\n");
		fflush(stdout);

		int num_inliers = (int)inliers_weak.size();

		vec3d *points_final = new vec3d[num_inliers];
		vec2d *projs_final = new vec2d[num_inliers];
		int *idxs_final = new int[num_inliers];
		int *keys_final = new int[num_inliers];
		int num_points_final = num_inliers;

		for (int i = 0; i < num_inliers; i++) {
			points_final[i] = points_solve[inliers_weak[i]];
			projs_final[i] = projs_solve[inliers_weak[i]];

			idxs_final[i] = idxs_solve[inliers_weak[i]];
			keys_final[i] = keys_solve[inliers_weak[i]];
		}

		if (refine_cameras_and_points) {
			inliers =
				refine_camera_and_points(data, num_points_final,
				points_final, projs_final, idxs_final,
				cameras, added_order, pt_views, &camera_new,
				true);
		}
		else {
			inliers = refine_camera_parameters(data, num_points_final,
				points_final, projs_final,
				idxs_final, &camera_new,
				NULL, !option_.fixed_focal_length, true,
				option_.estimate_distortion,
				option_.min_proj_error_threshold,
				option_.max_proj_error_threshold);
		}

		if ((int)inliers.size() < 8 || camera_new.f < 0.1 * data.width()) {
			printf("[BundleInitializeImage] Bad camera\n");
			if (success_out)
				*success_out = false;

			delete[] points_final;
			delete[] projs_final;
			delete[] idxs_final;
			delete[] keys_final;

			delete[] points_solve;
			delete[] projs_solve;
			delete[] projs_solve_orig;
			delete[] idxs_solve;
			delete[] keys_solve;

			camera_params_t dummy;
			return dummy;
		}

		/* Point the keys to their corresponding points */
		num_inliers = (int)inliers.size();
		for (int i = 0; i < num_inliers; i++) {
			int inlier_idx = inliers[i];
			// printf("[BundleInitializeImage] Connecting point [%d]\n",
			//        idxs_final[inlier_idx]);
			data.keys[keys_final[inlier_idx]].extra = idxs_final[inlier_idx];
			pt_views[idxs_final[inlier_idx]].
				push_back(ImageKey(camera_idx, keys_final[inlier_idx]));
		}
		fflush(stdout);

		delete[] points_final;
		delete[] projs_final;
		delete[] idxs_final;
		delete[] keys_final;

		delete[] points_solve;
		delete[] projs_solve;
		delete[] projs_solve_orig;
		delete[] idxs_solve;
		delete[] keys_solve;

		data.load_key_colors();
		data.camera.adjusted = true;

		return camera_new;
	}


	/* Set constraints on cameras */
	void SfM::set_camera_constraints(int cam_idx, camera_params_t *params) {

	}


	void SfM::set_focal_constraint(const ImageData &data, camera_params_t *params) {

	}


	/* Refine a given camera and the points it observes */
	std::vector<int> SfM::refine_camera_and_points(const ImageData &data,
		int num_points,
		vec3d *points, vec2d *projs,
		int *pt_idxs,
		camera_params_t *cameras,
		int *added_order,
		const std::vector < ImageKeyVector >
		&pt_views,
		camera_params_t *camera_out,
		bool remove_outliers)
	{
		// double error_thresh = 1.0e-6;
		double error_old = DBL_MAX;
		double derror;

		bool removed;
		std::vector<int> inliers_out;

		for (int i = 0; i < num_points; i++)
			inliers_out.push_back(i);

		int num_points_curr = num_points;
		vec3d *points_curr = new vec3d[num_points];
		vec2d *projs_curr = new vec2d[num_points];
		int *pt_idxs_curr = new int[num_points];

		memcpy(points_curr, points, sizeof(vec3d) * num_points);
		memcpy(projs_curr, projs, sizeof(vec2d) * num_points);
		memcpy(pt_idxs_curr, pt_idxs, sizeof(int) * num_points);

		do {
			removed = false;
			// do {
			{
				double error = 0.0;

				/* Refine the camera */
				refine_camera_parameters(data, num_points_curr,
					points_curr, projs_curr,
					pt_idxs_curr, camera_out, &error,
					!option_.fixed_focal_length, false, // true,
					option_.estimate_distortion,
					option_.min_proj_error_threshold,
					option_.max_proj_error_threshold);


#if 1 /* Change to 0 to skip point polishing */
				/* Refine the points */
				error = refine_points(num_points_curr, points_curr, projs_curr,
					pt_idxs_curr, cameras, added_order,
					pt_views, camera_out);

				printf("[RefineCameraAndPoints] "
					"Error (after point polishing): %0.3f\n", error);

				derror = error_old - error;
				error_old = error;
#endif
			}
			// } while (derror > error_thresh);

			if (remove_outliers) {
				/* Refine the camera once more and remove outliers */
				std::vector<int> inliers;
				inliers = refine_camera_parameters(data, num_points_curr,
					points_curr, projs_curr,
					pt_idxs_curr, camera_out,
					NULL, !option_.fixed_focal_length, true,
					option_.estimate_distortion,
					option_.min_proj_error_threshold,
					option_.max_proj_error_threshold);

				int num_inliers = (int)inliers.size();

				std::vector<int> inliers_out_next;
				if (num_inliers < num_points_curr) {
					removed = true;

					for (int i = 0; i < num_inliers; i++) {
						points_curr[i] = points_curr[inliers[i]];
						projs_curr[i] = projs_curr[inliers[i]];
						pt_idxs_curr[i] = pt_idxs_curr[inliers[i]];
						inliers_out_next.push_back(inliers_out[i]);
					}

					num_points_curr = num_inliers;
					inliers_out = inliers_out_next;
				}
			}
			else {
				refine_camera_parameters(data, num_points_curr,
					points_curr, projs_curr,
					pt_idxs_curr, camera_out,
					NULL, !option_.fixed_focal_length, false,
					option_.estimate_distortion,
					option_.min_proj_error_threshold,
					option_.max_proj_error_threshold);
			}

		} while (removed);

		delete[] points_curr;
		delete[] projs_curr;
		delete[] pt_idxs_curr;

		return inliers_out;
	}


	bool SfM::estimate_relative_pose2(int i1, int i2,
		camera_params_t &camera1,
		camera_params_t &camera2)
	{
		MatchIndex list_idx(i1, i2);
		if (i1 > i2)
			list_idx = MatchIndex(i2, i1);

		std::vector<KeypointMatch> &matches = match_table_.match_list(list_idx);
		int num_matches = (int)matches.size();

		double K1[9], K2[9];
		get_intrinsics(camera1, K1);
		get_intrinsics(camera2, K2);

		double R0[9], t0[3];
		int num_inliers = 0;

		num_inliers =
			estimate_pose_5_point(image_data_[i1].keys,
			image_data_[i2].keys,
			matches,
			512, /* m_fmatrix_rounds, 8 * m_fmatrix_rounds */
			0.25 * option_.fmatrix_threshold, // 0.003, // 0.004 /*0.001,*/ // /*0.5 **/ m_fmatrix_threshold, 
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


	double SfM::refine_points(int num_points, vec3d *points, vec2d *projs,
		int *pt_idxs, camera_params_t *cameras,
		int *added_order,
		const std::vector<ImageKeyVector> &pt_views,
		camera_params_t *camera_out)
	{
		double error = 0.0;

		/* Triangulate each of the points */
		for (int i = 0; i < num_points; i++) {
			int pt_idx = pt_idxs[i];

			int num_views = (int)pt_views[pt_idx].size() + 1;

			if (num_views < 2) continue;

			vec2d *pv = new vec2d[num_views];
			double *Rs = new double[9 * num_views];
			double *ts = new double[3 * num_views];

			for (int j = 0; j < num_views; j++) {
				camera_params_t *cam = NULL;

				if (j < num_views - 1) {
					int camera_idx = pt_views[pt_idx][j].first;
					int image_idx = added_order[camera_idx];
					int key_idx = pt_views[pt_idx][j].second;
					Keypoint &key = get_key(image_idx, key_idx);

					double p3[3] = { key.x, key.y, 1.0 };
					double K[9], Kinv[9];
					get_intrinsics(cameras[camera_idx], K);
					matrix_invert(3, K, Kinv);

					double p_n[3];
					matrix_product(3, 3, 3, 1, Kinv, p3, p_n);

					pv[j] = vec2d(p_n[0], p_n[1]);

					cam = cameras + camera_idx;
				}
				else {
					double p3[3] = { projs[i].x, projs[i].y, 1.0 };
					double K[9], Kinv[9];
					get_intrinsics(*camera_out, K);
					matrix_invert(3, K, Kinv);

					double p_n[3];
					matrix_product(3, 3, 3, 1, Kinv, p3, p_n);

					pv[j] = vec2d(p_n[0], p_n[1]);
					cam = camera_out;
				}

				memcpy(Rs + 9 * j, cam->R, 9 * sizeof(double));

				matrix_product(3, 3, 3, 1, cam->R, cam->t, ts + 3 * j);
				matrix_scale(3, 1, ts + 3 * j, -1.0, ts + 3 * j);
			}

			// points[i] = triangulate_n(num_views, pv, Rs, ts, &error_curr);
			double error_curr = 0.0;
			points[i] = triangulate_n_refine(points[i], num_views, pv, Rs, ts,
				&error_curr);

			vec2d pr = sfm_project_final(camera_out, points[i], 1,
				option_.estimate_distortion ? 1 : 0);

			double dx = pr.x - projs[i].x;
			double dy = pr.y - projs[i].y;

			error += dx * dx + dy * dy;

			delete[] pv;
			delete[] Rs;
			delete[] ts;
		}

		return sqrt(error / num_points);
	}
}