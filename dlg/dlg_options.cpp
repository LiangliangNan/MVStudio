#include "dlg_options.h"
#include "paint_canvas.h"
#include "main_window.h"
#include "../basic/logger.h"
#include "../basic/file_utils.h"
#include "../pointset/point_set.h"

#include <QDoubleValidator>


DlgOptions::DlgOptions(PaintCanvas* cvs, QWidget* parent)
: QDialog(parent)
, canvas_(cvs)
{
	setupUi(this);
	//setWindowFlags(Qt::WindowStaysOnTopHint);

	density_attr_name_ = "density";

	// default value
	default_octree_depth_ = 8;
	default_samples_per_node_ = "1.0";
	default_trim_value_ = "5.0";
	default_area_ratio_ = "0.001";
	
	lineEditSamplesPerNode->setValidator(new QDoubleValidator(1.0, 20, 3, this));
	lineEditTrimValue->setValidator(new QDoubleValidator(1.0, 20, 3, this));
	lineEditAreaRatio->setValidator(new QDoubleValidator(0.0, 1.0, 3, this));
	resetParameters();

	connect(hintButton, SIGNAL(toggled(bool)), this, SLOT(showHint(bool)));
	connect(defaultPrametersButton, SIGNAL(clicked()), this, SLOT(resetParameters()));
	connect(reconstructButton, SIGNAL(clicked()), this, SLOT(reconstruct()));
	connect(trimButton, SIGNAL(clicked()), this, SLOT(trim()));

	showHint(false);
}

DlgOptions::~DlgOptions()
{
}

void DlgOptions::resetParameters() {
	spinBoxOctreeDepth->setValue(default_octree_depth_);
	lineEditSamplesPerNode->setText(default_samples_per_node_);

	lineEditTrimValue->setText(default_trim_value_);
	lineEditAreaRatio->setText(default_area_ratio_);
}

void DlgOptions::showHint(bool b) {
	if (b) 
		widgetHint->show();
	else
		widgetHint->hide();
}

void DlgOptions::reconstruct() {
	PaintCanvas* cvs = dynamic_cast<PaintCanvas*>(canvas_);
	PointSet* pset = cvs->pointSet();
	if (pset) {
// 		int		octree_depth = spinBoxOctreeDepth->value();
// 		float	sampers_per_node = lineEditSamplesPerNode->text().toFloat();
// 
// 		PoissonReconstruction recon;
// 		recon.set_octree_depth(octree_depth);
// 		recon.set_sampers_per_node(sampers_per_node);
// 
// 		Map* mesh = recon.apply(pset, density_attr_name_);
// 		if (mesh) {
// 			std::string name = FileUtils::dir_name(pset->name()) + "/" + FileUtils::base_name(pset->name()) + "_poisson_reconstruction";
// 			mesh->set_name(name);
// 			cvs->mainWindow()->addObject(mesh, false, false);
// 			cur_model_name_ = mesh->name();
// 		}
// 		show();
	}
}


void DlgOptions::trim() {
// 	PaintCanvas* cvs = dynamic_cast<PaintCanvas*>(canvas_);
// 
// 	Map* mesh = canvas_->mesh();
// 	if (mesh) {
// 		float	trim_value = lineEditTrimValue->text().toFloat();
// 		float   area_ratio = lineEditAreaRatio->text().toFloat();
// 		bool    triangulate = false; // I can do it using my algorithm :-)
// 		int		smooth = 0;
// 
// 		if (MapVertexAttribute<float>::is_defined(mesh, density_attr_name_)) {
// 			MapVertexAttribute<float> density(mesh, density_attr_name_);
// 			float min_density = FLT_MAX;
// 			float max_density = -FLT_MAX;
// 			FOR_EACH_VERTEX_CONST(Map, mesh, it) {
// 				float value = density[it];
// 				min_density = std::min(min_density, value);
// 				max_density = std::max(max_density, value);
// 			}
// 
// 			if (trim_value <= min_density || trim_value >= max_density) {
// 				Logger::err("PoissonRecon") << "trim value (" << trim_value 
// 					<< ") out of density range [" << clip_precision(min_density, 2) << ", " 
// 					<< clip_precision(max_density, 2) << "]" << std::endl; 
// 				return;
// 			}
// 		} 
// 
// 		Map* trimmed_mesh = PoissonReconstruction::trim(mesh, density_attr_name_, trim_value, area_ratio, triangulate, smooth);
// 		if (trimmed_mesh) {
// 			std::string name = FileUtils::dir_name(mesh->name()) + "/" + FileUtils::base_name(mesh->name()) + "_trimmed";
// 			trimmed_mesh->set_name(name);
// 			cvs->mainWindow()->addObject(trimmed_mesh, false, false);
// 			cur_model_name_ = trimmed_mesh->name();
// 		}
// 		show();
// 	}
}
