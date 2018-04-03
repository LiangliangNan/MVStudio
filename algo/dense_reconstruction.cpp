
#include "dense_reconstruction.h"
#include "project.h"

#include "../basic/file_utils.h"
#include "../basic/progress.h"
#include "../basic/logger.h"
#include "../basic/stop_watch.h"
#include "../pointset/point_set_io.h"
#include "../sfm/sfm_util.h"
#include "../3rd_party/cmvs-pmvs/pmvs/findMatch.h"
#include "../3rd_party/cmvs-pmvs/pmvs/option.h"

#include <cassert>
#include <iomanip>



using namespace sfm;

DenseReconstruction::DenseReconstruction(Project* proj)
: project_(proj)
{
}

DenseReconstruction::~DenseReconstruction() {
}


bool DenseReconstruction::convert_bundler_to_pmvs() {
	if (!project_ || !project_->is_valid()) {
		Logger::warn(title()) << "invalid project" << std::endl;
		return false;
	}

	if (!FileUtils::is_file(project_->sfm_out_file)) {
		Logger::warn(title()) << "please run sparse reconstruction first" << std::endl;
		return false;
	}

	// read the list file
	const std::string& list_file = project_->sfm_list_file;
	std::vector<std::string> images;
	sfm::read_list_file(list_file, images);
	assert(images.size() > 0);

	// read the bundle file
	std::vector<camera_params_t> cameras;
	std::vector<point_t> points;
	const std::string& bundle_file = project_->sfm_out_file;
	sfm::read_bundle_file(bundle_file, cameras, points);

	//////////////////////////////////////////////////////////////////////////
	// write camera geometry in the PMVS file format 
	const std::string& pmvs_path = project_->pmvs_dir;
	const std::string& txt_path = project_->pmvs_txt_dir;
	sfm::write_pmvs(pmvs_path, txt_path, list_file, bundle_file, project_->pmvs_option_file, images, cameras);

	//////////////////////////////////////////////////////////////////////////

	// write new files
	std::string rd_list_file = pmvs_path + "/list.rd.txt";
	std::ofstream output(rd_list_file.c_str());
	if (output.fail()) {
		Logger::warn(title()) << "could not open file: \'" << rd_list_file << "\'" << std::endl;
		return false;
	}
	for (int i = 0; i < images.size(); i++) {
		if (cameras[i].f == 0.0)
			continue;
		output << images[i] << std::endl;
	}


	std::string rd_bundle_file = pmvs_path + "/bundle.rd.out";
	sfm::write_bundle_file(rd_bundle_file, cameras, points);

	//////////////////////////////////////////////////////////////////////////
	// write vis.dat file
	const std::string& vis_file = project_->pmvs_vis_data_file;
	write_vis_file(vis_file, cameras, points);

	//////////////////////////////////////////////////////////////////////////

	//undistort images

	Logger::out(title()) << "undistorting images..." << std::endl;

	const std::string& visualize_path = project_->pmvs_visualize_dir;
	int count = 0;

	ProgressLogger progress(images.size());
	for (int i = 0; i < images.size(); i++) {
		if (cameras[i].f == 0.0)
			continue;

		std::ostringstream string_stream;
		string_stream << visualize_path << "/" << std::setfill('0') << std::setw(4) << count << ".jpg";
		std::string out = string_stream.str();
		sfm::undistort_image(images[i], cameras[i], out);

		++ count;
		progress.next();
	}

	Logger::out(title()) << "undistorting images done" << std::endl;

	return true;
}

bool DenseReconstruction::convert_pba_to_pmvs() {
	if (!project_ || !project_->is_valid()) {
		Logger::warn(title()) << "invalid project" << std::endl;
		return false;
	}

	return true;
}

bool DenseReconstruction::run_pmvs(PointSet* pset) {
	if (!project_ || !project_->is_valid()) {
		Logger::warn(title()) << "invalid project" << std::endl;
		return false;
	}

	bool ready = false;
	if (project_->sfm_method() == SFM_BUNDLER)
		ready = convert_bundler_to_pmvs();
	else
		ready = convert_pba_to_pmvs();

	if (ready) {
		Logger::out(title()) << "running PMVS..." << std::endl;
		StopWatch w;

		PMVS3::Soption option;
		option.init_from_file(project_->pmvs_dir, project_->pmvs_option_file);

		PMVS3::CfindMatch findMatch;
		findMatch.init(option);
		findMatch.run();

		Logger::out(title()) << "PMVS done. Time: " << w.elapsed() << std::endl;

		Logger::out(title()) << "saving results..." << std::endl;
		w.start();
		//const std::vector<Ppatch>& patches = findMatch.m_pos.m_ppatches;

		findMatch.fillPointSet(pset);

		//Liangliang: change if you want to save intermediate results 
		//bool export_patch = false;
		//bool export_pSet = false;
		//findMatch.write(project_->pmvs_models_dir, export_patch, export_pSet);
		
		std::string ply_file = project_->pmvs_models_dir + "/dense.ply";
		PointSetIO::save(ply_file, pset);
		Logger::out(title()) << "saving results done. Time: " << w.elapsed() << std::endl;
	}

	return true;
}