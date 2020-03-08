
#include "sparse_reconstruction.h"
#include "project.h"

#include <easy3d/util/logging.h>
#include <easy3d/util/file_system.h>
#include <easy3d/util/stop_watch.h>


using namespace easy3d;

SparseReconstruction::SparseReconstruction(Project* proj) 
	: project_(proj)
{
}

SparseReconstruction::~SparseReconstruction() {
}


bool SparseReconstruction::apply(PointCloud* pset /* = nil */) {
	if (!project_ || !project_->is_valid()) {
		LOG(WARNING) << "invalid project" << std::endl;
		return false;
	}

	std::string listfocal = project_->sfm_list_file;
	std::string match_table = project_->sfm_match_table_file;
	if (!file_system::is_file(listfocal) || !file_system::is_file(match_table)) {
		LOG(WARNING) << "please run image matching first" << std::endl;
		return false;
	}

	std::vector<ProjectImage>& images = project_->images;
	for (std::size_t i = 0; i < images.size(); ++i) {
		if (images[i].ignored)
			continue;
		const std::string& basename = file_system::base_name(images[i].file);
		std::string key_file = project_->sfm_keys_dir + "/" + basename + ".key";
		if (!file_system::is_file(key_file)) {
			LOG(WARNING) << basename + ".key" << " doesn\'t exist" << std::endl;
			return false;
		}
	}

	StopWatch w;
	run_sfm(pset);

	return true;
}



void SparseReconstruction::run_sfm(PointCloud* pset) {
	if (!project_ || !project_->is_valid()) {
		LOG(WARNING) << "invalid project" << std::endl;
		return;
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
	//option.skip_full_bundle = true;

	sfm::SfM s(option);
	s.apply(pset);
}

