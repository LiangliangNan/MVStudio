#include <iostream>
#include "../basic/logger.h"
#include "../basic/file_utils.h"
#include "../basic/basic_types.h"


#include "../algo/project.h"
#include "../sfm/sfm.h"


int main(int argc, char* argv[]) {

	//std::string project_file = "../../data/ET/ET(small).mvsproj";
	std::string project_file = "../../data/building/building.mvsproj";
	Project* project_ = new Project;
	project_->load(project_file);


	if (!project_ || !project_->is_valid()) {
		std::cerr << "invalid project" << std::endl;
		return 0;
	}

	std::string listfocal = project_->sfm_list_file;
	std::string match_table = project_->sfm_match_table_file;
	if (!FileUtils::is_file(listfocal) || !FileUtils::is_file(match_table)) {
		std::cerr << "please run image matching first" << std::endl;
		return 0;
	}

	std::vector<ProjectImage>& images = project_->images;
	for (std::size_t i = 0; i < images.size(); ++i) {
		if (images[i].ignored)
			continue;
		const std::string& basename = FileUtils::base_name(images[i].file);
		std::string key_file = project_->sfm_keys_dir + "/" + basename + ".key";
		if (!FileUtils::is_file(key_file)) {
			std::cerr << basename + ".key" << " doesn\'t exist" << std::endl;
			return 0;
		}
	}

	SfmOption option;
	option.list_file = project_->sfm_list_file;
	option.match_table_file = project_->sfm_match_table_file;
	option.bundle_output_file = project_->sfm_out_file;
	option.bundle_output_base = project_->sfm_dump_file_base;
	option.output_directory = project_->sfm_output_dir;
	option.image_directory = project_->image_dir;
	option.key_directory = project_->sfm_keys_dir;
	option.fixed_focal_length = false;
	option.constrain_focal = true;
	option.constrain_focal_weight = 0.0001;
	option.estimate_distortion = true;

	sfm::SfM s(option);
	s.apply();


    return 0;
}

