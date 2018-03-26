#ifndef _SFM_SFM_UTIL_H_
#define _SFM_SFM_UTIL_H_


#include <vector>
#include "../sfm/camera.h"


namespace sfm {

	typedef struct {
		int image;
		int key;
		double x;
		double y;
	} view_t;

	typedef struct {
		double pos[3];
		double color[3];
		std::vector<view_t> views;
	} point_t;



	void read_bundle_file(const std::string& bundle_file,
		std::vector<camera_params_t> &cameras,
		std::vector<point_t> &points
		);


	void write_bundle_file(const std::string& bundle_file,
		const std::vector<camera_params_t> &cameras,
		const std::vector<point_t> &points
		);


	void read_list_file(const std::string& list_file, std::vector<std::string> &files);


	void write_pmvs(const std::string& pmvs_path,
		const std::string& txt_path,
		const std::string& list_file,
		const std::string& bundle_file,
		const std::string& pmvs_option_file,
		std::vector<std::string> images,
		std::vector<camera_params_t> &cameras
		);

	void write_vis_file(const std::string& vis_file,
		std::vector<camera_params_t> &cameras,
		std::vector<point_t> &points
		);

	void undistort_image(const std::string &in_file,
		const camera_params_t &camera,
		const std::string &out_file
		);
};


#endif
