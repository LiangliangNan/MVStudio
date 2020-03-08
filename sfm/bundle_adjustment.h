#ifndef _SFM_UTIL_H_
#define _SFM_UTIL_H_


#include "camera.h"
#include "../math/math_types.h"

#include <vector>


namespace sfm {


	/* Compute an updated rotation matrix given the initial rotation (R)
	* and the correction (w) */
	void rot_update(double *R, double *w, double *Rnew);

	void sfm_project(camera_params_t *init, double *K,
		double *w, double *dt, double *b, double *p,
		int explicit_camera_centers);

	void sfm_project_rd(camera_params_t *init, double *K, double *k,
		double *R, double *dt, double *b, double *p,
		int undistort, int explicit_camera_centers);

	vec2d sfm_project_final(camera_params_t *params, vec3d pt,
		int explicit_camera_centers, int undistort);

	void bundle_adjustment(int num_pts, int num_cameras, int ncons,
		char *vmask,
		double *projections,
		int est_focal_length,
		int const_focal_length,
		int undistort,
		int explicit_camera_centers,
		camera_params_t *init_camera_params,
		vec3d *init_pts,
		int use_constraints,
		int use_point_constraints,
		vec3d *points_constraints,
		double point_constraint_weight,
		int fix_points,
		double eps2,
		double *Vout,
		double *Sout,
		double *Uout, double *Wout);


	/* Refine the position of a single camera */
	void camera_refine(int num_points, vec3d *points, vec2d *projs,
		camera_params_t *params, int adjust_focal,
		int estimate_distortion);
}



















// functions that should be put somewhere else
#include "image_data.h"
namespace sfm {
	void initialize_camera_parameters(const ImageData &data, camera_params_t &camera);

	std::vector<int> refine_camera_parameters(const ImageData &data,
		int num_points,
		vec3d *points, vec2d *projs,
		int *pt_idxs, camera_params_t *camera,
		double *error_out,
		bool adjust_focal,
		bool remove_outliers,
		//                                         bool optimize_for_fisheye,
		bool estimate_distortion,
		double min_proj_error_threshold,
		double max_proj_error_threshold
		);

	/* Use a 180 rotation to fix up the intrinsic matrix */
	void fix_intrinsics(double *P, double *K, double *R, double *t);

	bool find_and_verify_camera(int num_points, vec3d *points_solve, vec2d *projs_solve,
		int *idxs_solve,
		double *K, double *R, double *t,
		double proj_estimation_threshold,
		double proj_estimation_threshold_weak,
		std::vector<int> &inliers,
		std::vector<int> &inliers_weak,
		std::vector<int> &outliers);

	/* Triangulate two points */
	vec3d triangulate(vec2d p, vec2d q,
		camera_params_t c1, camera_params_t c2,
		double &proj_error, bool &in_front, double &angle,
		bool explicit_camera_centers);
}
#endif
