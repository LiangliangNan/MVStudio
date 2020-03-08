#ifndef _SFM_H_
#define _SFM_H_

#include "sfm_option.h"
#include "image_data.h"
#include "match_table.h"
#include "geometry.h"

#include <string>
#include <vector>




enum SfMMethod {
	SFM_BUNDLER,
	SFM_PBA
};


namespace easy3d {
    class PointCloud;
}


namespace sfm {

	typedef std::pair<int, int> ImagePair;

	
	class SfM
	{
	public:
		SfM(SfmOption opt);
		~SfM();

		static std::string title() { return "SfM"; }

		// On return, 'pset' will also be filled with the computed sparse point cloud 
		// if pset is not NULL. And meanwhile, the point cloud will also be saved in 
		// the file specified by 'SfmOption.bundle_output_file'.  
		void apply(easy3d::PointCloud* pset = nullptr);

	private:

		/* Return the number of images */
		int  num_of_images() const;

		/* Return the number of matches for image pair (img1, img2) */
		int  num_of_matches(int img1, int img2) const;

		/* Get keys */
		Keypoint& get_key(int img, int key);
		int num_of_keys(int img);

		// ----------- file management ------------

		void load_image_names();
		void load_match_table();
		void load_keys();

		void read_camera_constraints(const std::string& file);
		void read_point_constraints(const std::string& file);

		void read_geometric_constraints(const std::string& file);
		void write_geometric_constraints(const std::string& file);

		/* Dump an output file containing information about the current state of the scene */
		void dump_output_file(const std::string& filename,
			int num_images, int num_cameras, int num_points, int* added_order,
			camera_params_t* cameras, easy3d::dvec3* points, easy3d::ivec3* colors,
			std::vector<ImageKeyVector>& pt_views);

		/* Write points to a ply file */
		void dump_points_to_ply(const std::string& filename,
			int num_points, int num_cameras,
			easy3d::dvec3* points, easy3d::ivec3* colors, camera_params_t* cameras);


		// ----------- matches  ------------

		void remove_all_matches();

		void set_match(int img1, int img2);
		void remove_match(int img1, int img2);
		bool images_match(int img1, int img2);

		/* Prune points that match to multiple targets */
		void prune_multiple_matches();

		/* Make match lists symmetric */
		void make_match_lists_symmetric();

		/* Remove matches close to the edges of the images */
		void remove_border_matches();
		/* Remove matches close to the bottom edge of the images */
		void remove_bottom_matches();
		/* Remove matches close to the edges of the two given images */
		void remove_border_matches(int img1, int img2, int border_width);
		/* Remove matches close to the bottom edge of the two given images */
		void remove_bottom_matches(int img1, int img2, int border_width);

		/* Use the bundle-adjusted points to create a new set of matches */
		void set_matches_from_points(int threshold = 0);

		// ----------- geometric constraints ------------

		/* Compute geometric information about image pairs */
		void compute_geometric_constraints(int new_image_start = 0);

		/* Compute epipolar geometry between all matching images */
		void compute_epipolar_geometry(bool remove_bad_matches, int new_image_start = 0);
		/* Compute epipolar geometry between a given pair of images */
		bool compute_epipolar_geometry(int idx1, int idx2, bool remove_bad_matches);

		/* Compute transforms between all matching images */
		void compute_transforms(bool remove_bad_matches, int new_image_start = 0);
		/* Compute a transform between a given pair of images */
		bool compute_transform(int idx1, int idx2, bool remove_bad_matches);

		/* Compute a set of tracks that explain the matches */
		void compute_tracks(int new_image_start = 0);

		void set_tracks(int image);
		void set_matches_from_tracks(int img1, int img2);		/* Use tracks to setup matches */
		int  get_num_track_matches(int img1, int img2);


		// ----------- bundle adjustment ------------
		void  run_sfm(easy3d::PointCloud* pset);

		double run_bundle_adjustment(int num_pts, int num_cameras, int start_camera,
			bool fix_points, camera_params_t *init_camera_params,
			easy3d::dvec3 *init_pts, int *added_order, easy3d::ivec3 *colors,
			std::vector<ImageKeyVector> &pt_views, double eps2 = 1.0e-12,
			double *S = NULL, double *U = NULL, double *V = NULL,
			double *W = NULL, bool remove_outliers = true);

		/* Initialize the bundle adjustment procedure (loading an existing
		* model if one exists) */
		void initialize_bundle_adjust(int &num_init_cams,
			int *added_order,
			int *added_order_inv,
			camera_params_t *cameras,
			easy3d::dvec3 *points, easy3d::ivec3 *colors,
			std::vector<ImageKeyVector> &pt_views);

		/* Find the camera with the most matches to existing points */
		int find_camera_with_most_matches(int num_cameras, int num_points,
			int *added_order,
			int &parent_idx, int &max_matches,
			const std::vector<ImageKeyVector> &pt_views);

		/* Find all cameras with at least N matches to existing points */
		std::vector<ImagePair> find_cameras_with_n_matches(int n,
			int num_cameras,
			int num_points,
			int *added_order,
			const std::vector<ImageKeyVector> &pt_views);

		/* Triangulate a subtrack */
		easy3d::dvec3 triangulate_n_views(const ImageKeyVector &views,
			int *added_order, camera_params_t *cameras,
			double &error, bool explicit_camera_centers);

		/* Add new points to the bundle adjustment */
		int bundle_adjust_add_all_new_points(int num_points, int num_cameras,
			int *added_order,
			camera_params_t *cameras,
			easy3d::dvec3 *points, easy3d::ivec3 *colors,
			double reference_baseline,
			std::vector<ImageKeyVector> &pt_views,
			double max_reprojection_error = 16.0,
			int min_views = 2);

		/* Remove bad points and cameras from a reconstruction */
		int remove_bad_points_and_cameras(int num_points, int num_cameras,
			int *added_order,
			camera_params_t *cameras,
			easy3d::dvec3 *points, easy3d::ivec3 *colors,
			std::vector<ImageKeyVector> &pt_views);

		/* Pick a good initial pair of cameras to bootstrap the bundle
		* adjustment */
		void bundle_pick_initial_pair(int &i_best, int &j_best,
			bool use_init_focal_only);

		/* Setup the initial camera pair for bundle adjustment */
		int setup_initial_camera_pair(int i_best, int j_best,
			double &init_focal_length_0,
			double &init_focal_length_1,
			camera_params_t *cameras,
			easy3d::dvec3 *points, easy3d::ivec3 *colors,
			std::vector<ImageKeyVector> &pt_views);

		/* Initialize an image for bundle adjustment */
		camera_params_t bundle_initialize_image(ImageData &data,
			int image_idx, int camera_idx,
			int num_cameras, int num_points,
			int *added_order, easy3d::dvec3 *points,
			camera_params_t *parent,
			camera_params_t *cameras,
			std::vector<ImageKeyVector> &pt_views,
			bool *success_out = NULL,
			bool refine_cameras_and_points = false);

		/* Set constraints on cameras */
		void set_camera_constraints(int cam_idx, camera_params_t *params);
		void set_focal_constraint(const ImageData &data, camera_params_t *params);

		/* Refine a set of 3D points */
		double refine_points(int num_points, easy3d::dvec3 *points, easy3d::dvec2 *projs,
			int *pt_idxs, camera_params_t *cameras,
			int *added_order,
			const std::vector<ImageKeyVector> &pt_views,
			camera_params_t *camera_out);

		/* Refine a given camera and the points it observes */
		std::vector<int> refine_camera_and_points(const ImageData &data,
			int num_points,
			easy3d::dvec3 *points, easy3d::dvec2 *projs,
			int *pt_idxs,
			camera_params_t *cameras,
			int *added_order,
			const std::vector<ImageKeyVector>
			&pt_views,
			camera_params_t *camera_out,
			bool remove_outliers);

		bool estimate_relative_pose2(int i1, int i2,
			camera_params_t &camera1,
			camera_params_t &camera2);

	private:
		SfmOption			option_;

		std::vector<ImageData>	image_data_;	/* Image data */
		std::vector<PointData>  point_data_;	/* Information about 3D points in the scene */
		std::vector<TrackData>  track_data_;	/* Information about the detected 3D tracks */

		MatchTable			match_table_;
		TransformHashMap		transforms_;

		/* Images to use as the initial pair during bundle adjustment */
		int				initial_pair_[2];   

		/* Point constraint fields */
		easy3d::dvec3*	point_constraints_;

	};

}

#endif
