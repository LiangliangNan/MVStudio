
#ifndef _DIRECTORIES_H
#define _DIRECTORIES_H


#include "../basic/file_utils.h"
#include "../sfm/sfm.h"

#include <string>
#include <vector>



struct ProjectImage {
	std::string	file;
	bool		ignored;
};


class Project
{
public:
	Project() {}
	~Project() {}

	static std::string title() { return "Project"; }

	// --------------- create -------------------
	// create a project (providing images)
	bool create(const std::string& proj_file, const std::string& img_dir);

	// create empty project (without images)
	bool create(const std::string& proj_file);
	// set images. Existing images will be replaced
	bool add_images(const std::string& img_dir);	

	void set_ignore_image(const std::string& image, bool ign);

	// --------------- read / write -------------------
	bool load(const std::string& file);
	bool save(const std::string& file) const;

	// ------------------------------------------------
	bool is_valid();

public:
	std::string  image_dir;
	std::string  project_dir;

	std::vector<ProjectImage> images;

	// sfm
	std::string  sfm_keys_dir;
	std::string  sfm_output_dir;
	std::string  sfm_list_file;
	std::string  sfm_match_table_file;
	std::string  sfm_out_file;
	std::string  sfm_dump_file_base;

	// pmvs
	std::string  pmvs_dir;
	std::string  pmvs_models_dir;
	std::string  pmvs_txt_dir;
	std::string  pmvs_visualize_dir;
	std::string  pmvs_vis_data_file;
	std::string  pmvs_option_file;

private:
	void initialize();
};


#endif