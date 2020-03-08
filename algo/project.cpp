#include "project.h"
#include <easy3d/util/logging.h>
#include <easy3d/util/file_system.h>
#include <fstream>

using namespace easy3d;

bool Project::create(const std::string& proj_file, const std::string& img_dir, SfMMethod sfm) {
	project_dir = file_system::parent_directory(proj_file);
	image_dir = img_dir;
	sfm_method_ = sfm;

	std::vector<std::string> entries;
	file_system::get_directory_entries(image_dir, entries, false);

	images.clear();
	for (std::size_t i = 0; i < entries.size(); ++i) {
		const std::string& name = entries[i];
		std::string ext = file_system::extension(name);
		if (ext == "jpg" || ext == "JPG") {
			ProjectImage img;
			img.file = image_dir + "/" + name;
			img.ignored = false;
			images.push_back(img);
		}
	}

	if (images.size() < 3) {
		LOG(WARNING) << "too few (only " << images.size() << ") images provided";
		return false;
	}

	initialize();

	if (is_valid())
		return save(proj_file);

	return false;
}


bool Project::create(const std::string& proj_file, SfMMethod sfm) {
	project_dir = file_system::parent_directory(proj_file);
	sfm_method_ = sfm;

	images.clear();
	initialize();

	return save(proj_file);
}


bool Project::add_images(const std::string& img_dir) {
	image_dir = img_dir;

	std::vector<std::string> entries;
	file_system::get_directory_entries(image_dir, entries, false);

	images.clear();
	for (std::size_t i = 0; i < entries.size(); ++i) {
		const std::string& name = entries[i];
		std::string ext = file_system::extension(name);
		if (ext == "jpg" || ext == "JPG") {
			ProjectImage img;
			img.file = image_dir + "/" + name;
			img.ignored = false;
			images.push_back(img);
		}
	}

	if (images.size() < 3) {
		LOG(WARNING) << "too few (only " << images.size() << ") images provided";
		return false;
	}

	return true;
}


void Project::set_ignore_image(const std::string& image, bool ign) {
	for (std::size_t i = 0; i < images.size(); ++i) {
		const std::string&  file = images[i].file;
		if (image == file) {
			images[i].ignored = ign;
		}
	}
}


void Project::initialize() {
	// sift features and image matching
	sfm_keys_dir = project_dir + "/keys";
	if (!file_system::is_directory(sfm_keys_dir))
		file_system::create_directory(sfm_keys_dir);
	sfm_list_file = sfm_keys_dir + "/list_focal.txt";
	sfm_match_table_file = sfm_keys_dir + "/match_table.txt";

	// sfm
	sfm_output_dir = project_dir + "/sfm";	
	if (!file_system::is_directory(sfm_output_dir))
		file_system::create_directory(sfm_output_dir);
	sfm_out_file = sfm_output_dir + "/bundle.out";
	sfm_dump_file_base = "/bundle_";

	// pmvs
	pmvs_dir = project_dir + "/pmvs";		if (!file_system::is_directory(pmvs_dir))		file_system::create_directory(pmvs_dir);
	pmvs_models_dir = pmvs_dir + "/models";	if (!file_system::is_directory(pmvs_models_dir))	file_system::create_directory(pmvs_models_dir);
	pmvs_txt_dir = pmvs_dir + "/txt";		if (!file_system::is_directory(pmvs_txt_dir))	file_system::create_directory(pmvs_txt_dir);
	pmvs_visualize_dir = pmvs_dir + "/visualize";	if (!file_system::is_directory(pmvs_visualize_dir))	file_system::create_directory(pmvs_visualize_dir);
	pmvs_vis_data_file = pmvs_dir + "/vis.dat";
	pmvs_option_file = pmvs_dir + "/pmvs_options.txt";
}


void Project::set_sfm_method(SfMMethod sfm) {
	sfm_method_ = sfm;
}


bool Project::is_valid() {
	// check if all directories exist
	if (!file_system::is_directory(project_dir)) {
		if (!file_system::create_directory(project_dir))
			return false;
	}
	if (!file_system::is_directory(sfm_keys_dir)) {
		if (!file_system::create_directory(sfm_keys_dir))
			return false;
	}	
	if (!file_system::is_directory(sfm_output_dir)) {
		if (!file_system::create_directory(sfm_output_dir))
			return false;
	}
	if (!file_system::is_directory(pmvs_dir)) {
		if (!file_system::create_directory(pmvs_dir))
			return false;
	}
	if (!file_system::is_directory(pmvs_models_dir)) {
		if (!file_system::create_directory(pmvs_models_dir))
			return false;
	}
	if (!file_system::is_directory(pmvs_txt_dir)) {
		if (!file_system::create_directory(pmvs_txt_dir))
			return false;
	}
	if (!file_system::is_directory(pmvs_visualize_dir)) {
		if (!file_system::create_directory(pmvs_visualize_dir))
			return false;
	}

	if (!file_system::is_directory(image_dir)) {
		if (!file_system::create_directory(image_dir))
			return false;
	}

	if (images.size() < 3) {
		LOG(WARNING) << "too few (only " << images.size() << ") images provided";
		return false;
	}

	//////////////////////////////////////////////////////////////////////////
	
	// check if all images exist
	unsigned int count = 0;
	for (std::size_t i = 0; i < images.size(); ++i) {
		const std::string&  file = images[i].file;
		if (!file_system::is_file(file)) {
			LOG(WARNING) << "missing image \'" << file_system::simple_name(file) << "\'";
			images[i].ignored = true;
			count++;
		}
	}

	return true;
}


bool Project::load(const std::string& file) {
	project_dir = file_system::parent_directory(file);

	std::ifstream input(file.c_str());
	if (input.fail()) {
		LOG(WARNING) << "could not load project \'" << file << "\'";
		return false;
	}

	std::string dummy;
	getline(input, dummy);  // skip the comments

	std::string sfm;
	input >> dummy >> sfm;

	if (sfm == "BUNDLER")
		sfm_method_ = SFM_BUNDLER;
	else if (sfm == "PBA")
		sfm_method_ = SFM_PBA;
	else {
		LOG(WARNING) << "unknown SfM method \'" << sfm << "\'. Use default \'BUNDLER\'";
	}

	//////////////////////////////////////////////////////////////////////////

	input >> dummy >> image_dir;

	images.clear();
	int num;
	input >> dummy >> num;
	for (int i = 0; i < num; ++i) {
		std::string name;
		bool ignored;
		input >> name >> ignored;
		std::string ext = file_system::extension(name);
		if (ext != "jpg" && ext != "JPG") {
			LOG(WARNING) << "unsupported image format \'" << ext << "\' (ignored)" << std::endl;
			continue;
		}

		if (!file_system::is_file(name)) {
			LOG(WARNING) << "missing image \'" << file_system::simple_name(name) << "\'";
			continue;
		}
		ProjectImage img;
		img.file = name;
		img.ignored = ignored;
		images.push_back(img);
	}

	initialize();

	return true;
}


bool Project::save(const std::string& file) const {
	std::ofstream output(file.c_str());
	if (output.fail()) {
		LOG(WARNING) << "could not save project \'" << file << "\'";
		return false;
	}

	output << "# MVStudio Project File. Saved by liangliang.nan@gmail.com" << std::endl;

	if (sfm_method_ == SFM_PBA)
		output << std::endl << "[SfMMethod]: " << "PBA" << std::endl;
	else
		output << std::endl << "[SfMMethod]: " << "BUNDLER" << std::endl;
	
	output << std::endl << "[ImageDirectory]: " << image_dir << std::endl;
	output << "[ImageNumber]: " << images.size() << std::endl;
	for (std::size_t i = 0; i < images.size(); ++i) {
		output << "     " << images[i].file << "     " << images[i].ignored << std::endl;
	}

	return true;
}