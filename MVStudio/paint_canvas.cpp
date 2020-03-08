#include "paint_canvas.h"
#include "main_window.h"
#include "../3rd_party/QGLViewer/manipulatedCameraFrame.h"
#include <easy3d/util/file_system.h>
#include "../pointCloud/point_set.h"
#include "../pointCloud/point_set_io.h"
#include "../pointCloud/point_set_render.h"
#include "../opengl/opengl_info.h"
#include "../algo/image_matching.h"
#include "../algo/sparse_reconstruction.h"
#include "../algo/dense_reconstruction.h"
#include "../sfm/sfm_util.h"

#include <QFileDialog>
#include <QMouseEvent>
#include <QMessageBox>
#include <QColorMap>
#include <QColorDialog>

#include <cassert>
#include <fstream>

using namespace sfm;


PaintCanvas::PaintCanvas(QWidget *parent)
	: QGLViewer(parent)
	, point_set_(nil)
	, coord_system_region_size_(150)
	, show_coord_sys_(true)
	, project_(nil)
{
	main_window_ = dynamic_cast<MainWindow*>(parent);

	render_ = new PointSetRender;

	light_pos_ = dvec3(0.27f, 0.27f, 0.92f);

	//////////////////////////////////////////////////////////////////////////

	// Move camera according to viewer type (on X, Y or Z axis)
	// e.g. 'GL_FRONT' equals to 
	camera()->setViewDirection(qglviewer::Vec(0.0, 1.0, 0.0));

	camera()->lookAt(sceneCenter());
	camera()->setType(qglviewer::Camera::PERSPECTIVE);
	camera()->showEntireScene();
}


PaintCanvas::~PaintCanvas() {
	delete render_;
}


void PaintCanvas::init()
{
	LOG(INFO) << "initializing..." << std::endl ;

 	//////////////////////////////////////////////////////////////////////////
 	GLenum err = glewInit();
 	if (GLEW_OK != err) {
 		// Problem: glewInit failed, something is seriously wrong. 
 		LOG(ERROR) << glewGetErrorString(err) << std::endl ;
 	}
	LOG(INFO) << "Using GLEW: " << GLInfo::glew_version() << std::endl ;
	LOG(INFO) << "GL Vendor: " << GLInfo::gl_vendor() << std::endl ;
	LOG(INFO) << "GL Renderer: " << GLInfo::gl_renderer() << std::endl ;
	LOG(INFO) << "GL Version: " << GLInfo::gl_version() << std::endl ;

	LOG(INFO) << "GLSL Version: " << GLInfo::glsl_version() << std::endl ;
	//LOG(INFO) << "GL_Extensions: " << GLInfo::gl_extensions() << std::endl;

	//////////////////////////////////////////////////////////////////////////

	setStateFileName("");

	// Default value is 0.005, which is appropriate for most applications. 
	// A lower value will prevent clipping of very close objects at the 
	// expense of a worst Z precision.
	camera()->setZNearCoefficient(0.005f);

	// Default value is square root of 3.0 (so that a cube of size 
	// sceneRadius() is not clipped).
	camera()->setZClippingCoefficient(std::sqrt(3.0f));

	camera()->setViewDirection(qglviewer::Vec(0.0, 1.0, 0.0));
	camera()->setType(qglviewer::Camera::PERSPECTIVE);
	showEntireScene();

	camera()->frame()->setSpinningSensitivity(/*1.0f*/100.0f);
	setMouseTracking(true);

	//////////////////////////////////////////////////////////////////////////
	// I like the inverse direction
	setShortcut(MOVE_CAMERA_LEFT, Qt::Key_Right);
	setShortcut(MOVE_CAMERA_RIGHT, Qt::Key_Left);
	setShortcut(MOVE_CAMERA_UP,	Qt::Key_Down);
	setShortcut(MOVE_CAMERA_DOWN, Qt::Key_Up);

  	setMouseBinding(Qt::ShiftModifier, Qt::LeftButton, CAMERA, SCREEN_ROTATE);
  	setMouseBinding(Qt::ControlModifier, Qt::LeftButton, CAMERA, ZOOM_ON_REGION);

	//////////////////////////////////////////////////////////////////////////

	glEnable(GL_DEPTH_TEST);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);  //GL_ONE_MINUS_DST_ALPHA

	//////////////////////////////////////////////////////////////////////////

	QColor bkgrd_color = Qt::white;
	setBackgroundColor(bkgrd_color);
	
	//////////////////////////////////////////////////////////////////////////

	//float pos[] = { 0.5f, 0.5f, 0.5f, 1.0f };
	float pos[] = {static_cast<float>(light_pos_.x), static_cast<float>(light_pos_.y), static_cast<float>(light_pos_.z), 0.0f };
	glLightfv(GL_LIGHT0, GL_POSITION, pos);

	glEnable(GL_LIGHT0);		
	glEnable(GL_LIGHTING);
	glEnable(GL_NORMALIZE);
}

void PaintCanvas::setLightPosition(const dvec3& pos) {
	light_pos_ = pos;

	makeCurrent();

	glPushMatrix();
	glLoadIdentity();
	// 0 = infinite light 
	// 1 = local light
	float p[] = {static_cast<float>(light_pos_.x), static_cast<float>(light_pos_.y), static_cast<float>(light_pos_.z), 0.0f };
	glLightfv(GL_LIGHT0, GL_POSITION, p);
	glPopMatrix();	
	
	update_graphics(); 
}


void PaintCanvas::showPointAsSphere(bool b) {
	render_->set_points_as_spheres(b);
	update_graphics();
}


void PaintCanvas::increasePointSize() {
	float size = render_->point_size();
	++size;
	render_->set_point_size(size);
	update_graphics();
}


void PaintCanvas::decreasePointSize() {
	float size = render_->point_size();
	if (size >= 2)
		--size;
	render_->set_point_size(size);
	update_graphics();
}

void PaintCanvas::draw() {
	if (show_coord_sys_)  {
		glEnable(GL_MULTISAMPLE);
		drawCornerAxis();
	}
	
	if (point_set_) {
		glDisable(GL_MULTISAMPLE);
		render_->draw();
	}

	// drawCameras();
}


void PaintCanvas::mousePressEvent(QMouseEvent* e)
{
	QGLViewer::mousePressEvent(e);
}


void PaintCanvas::mouseMoveEvent(QMouseEvent* e)
{
	bool found = false;
	qglviewer::Vec p = camera()->pointUnderPixel(e->pos(), found);
	main_window_->showCoordinateUnderMouse(dvec3(p.x, p.y, p.z), found);

	QGLViewer::mouseMoveEvent(e); 
}


void PaintCanvas::mouseReleaseEvent(QMouseEvent* e)
{	
	QGLViewer::mouseReleaseEvent(e);
}


void PaintCanvas::keyPressEvent( QKeyEvent *e ) {
	switch (e->key())
	{
	case Qt::Key_Escape: // I don't want to close the window
		break;

	default:
		QGLViewer::keyPressEvent(e);
	}
}


void PaintCanvas::snapshotScreen() {
	bool need_hide = show_coord_sys_;
	if (need_hide)
		show_coord_sys_ = false;  // hide the coord system temporally

	saveSnapshot(false);

	if (need_hide) {
		show_coord_sys_ = true;
		update_graphics();
	}
}


void PaintCanvas::copyCamera() {
	if (camera()->frame()) {
		const qglviewer::Vec pos = camera()->frame()->position();
		const qglviewer::Quaternion q = camera()->frame()->orientation();

		QString cam_str = QString("%1 %2 %3 %4 %5 %6 %7")
			.arg(pos[0])
			.arg(pos[1])
			.arg(pos[2])
			.arg(q[0])
			.arg(q[1])
			.arg(q[2])
			.arg(q[3]);
		qApp->clipboard()->setText(cam_str);
	} 
}


void PaintCanvas::pasteCamera() {
	// get the camera parameters from clipboard string
	QString str = qApp->clipboard()->text();
	QStringList list = str.split(" ", QString::SkipEmptyParts);
	if(list.size() != 7) {
		LOG(ERROR) << "camera not available in clipboard" << std::endl;
		return;
	}

	float vec[3];
	for(int i=0; i<3; ++i){
		bool ok;
		vec[i] = list[i].toFloat(&ok);
		if(!ok) {
			LOG(ERROR) << "camera not available in clipboard" << std::endl;
			return;
		}
	}

	double orient[4];
	for(int i=0; i< 4; ++i) {
		bool ok;
		orient[i] = list[i+3].toDouble(&ok);
		if(!ok) {
			LOG(ERROR) << "camera not available in clipboard" << std::endl;
			return;
		}
	}

	//////////////////////////////////////////////////////////////////////////

	// change view 
	qglviewer::Frame new_frame;
	new_frame.setPosition(qglviewer::Vec(vec[0], vec[1], vec[2]));
	new_frame.setOrientation(orient[0], orient[1], orient[2], orient[3]);

	float duration = 0.5f;
	camera()->interpolateTo(new_frame, duration); 
	update_graphics();
}


void PaintCanvas::setLighting(bool b) {
	render_->set_lighting(b);
	update_graphics();
}

void PaintCanvas::fit() {
	fitScreen();
}

void PaintCanvas::update_graphics() {
	update();  
}

void PaintCanvas::update_all() {
	main_window_->updateStatusBar();
	update_graphics();

	// This approach has significant drawbacks. For example, imagine you wanted to perform two such loops 
	// in parallel-calling one of them would effectively halt the other until the first one is finished 
	// (so you can't distribute computing power among different tasks). It also makes the application react
	// with delays to events. Furthermore the code is difficult to read and analyze, therefore this solution
	// is only suited for short and simple problems that are to be processed in a single thread, such as 
	// splash screens and the monitoring of short operations.
	QCoreApplication::processEvents();  
}


void PaintCanvas::setProjectionMode(bool b) {
	if (b)
		camera()->setType(qglviewer::Camera::PERSPECTIVE);
	else
		camera()->setType(qglviewer::Camera::ORTHOGRAPHIC);

	update_graphics();
}

void PaintCanvas::showCoordinateSystem(bool b) { 
	show_coord_sys_= b;
	update_graphics();
}


void PaintCanvas::fitScreen() {
	if (point_set_) {
		Box3f box = point_set_->bounding_box();

		qglviewer::Vec vmin(box.x_min(), box.y_min(), box.z_min());
		qglviewer::Vec vmax(box.x_max(), box.y_max(), box.z_max());

		//Vec vmin(-2.0f, -2.0f, -2.0f);
		//Vec vmax(2.0f, 2.0f, 2.0f);
		setSceneBoundingBox(vmin, vmax);
		showEntireScene();
		update_graphics();
	}
}


void PaintCanvas::drawCornerAxis()  
{	
	glEnable(GL_MULTISAMPLE);

	// The viewport and the scissor are changed to fit the lower left
	// corner. Original values are saved.
	int viewport[4];
	int scissor[4];
	glGetIntegerv(GL_VIEWPORT, viewport);
	glGetIntegerv(GL_SCISSOR_BOX, scissor);

	//////////////////////////////////////////////////////////////////////////

	// Axis viewport size, in pixels
	glViewport(0, 0, coord_system_region_size_, coord_system_region_size_);
	glScissor(0, 0, coord_system_region_size_, coord_system_region_size_);

	// The Z-buffer is cleared to make the axis appear over the
	// original image.
	glClear(GL_DEPTH_BUFFER_BIT);

	// Tune for best line rendering
	glLineWidth(3.0);

	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	glOrtho(-1, 1, -1, 1, -1, 1);

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();
	glMultMatrixd(camera()->orientation().inverse().matrix());

	float axis_size = 0.9f; // other 0.2 space for drawing the x, y, z labels
	drawAxis(axis_size); 

	// Draw text id
	glColor3f(0, 0, 0);

	// Liangliang: It seems the renderText() func disables multi-sample.
	// Is this a bug in Qt ?
	GLboolean anti_alias = glIsEnabled(GL_MULTISAMPLE);
	const_cast<PaintCanvas*>(this)->renderText(axis_size, 0, 0, "X");
	const_cast<PaintCanvas*>(this)->renderText(0, axis_size, 0, "Y");
	const_cast<PaintCanvas*>(this)->renderText(0, 0, axis_size, "Z");
	if (anti_alias)
		glEnable(GL_MULTISAMPLE);

	glMatrixMode(GL_PROJECTION);
	glPopMatrix();

	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();

	//////////////////////////////////////////////////////////////////////////

	// The viewport and the scissor are restored.
	glScissor(scissor[0], scissor[1], scissor[2], scissor[3]);
	glViewport(viewport[0], viewport[1], viewport[2], viewport[3]);
}


dvec2 PaintCanvas::projectionOf(const dvec3& p) {    // point to screen
	qglviewer::Vec v = camera()->projectedCoordinatesOf(qglviewer::Vec(p.x, p.y, p.z));
	return dvec2(v.x, v.y);
}

dvec3 PaintCanvas::unProjectionOf(double winx, double winy, double winz) {  // screen to point
	qglviewer::Vec v = camera()->unprojectedCoordinatesOf(qglviewer::Vec(winx, winy, winz));
	return dvec3(v.x, v.y, v.z);
}


void PaintCanvas::setPointSet(PointCloud* pset, bool fit) {
	if (point_set_) {
		delete point_set_;
		render_->set_pointset(nil);
	}

	point_set_ = pset;
	point_set_->set_canvas(this);
	render_->set_pointset(point_set_);
	if (fit)
		fitScreen();
}


bool PaintCanvas::creatProject(const QString& file) {
	if (!project_)
		project_ = new Project;

	if (project_->create(file.toStdString())) {
		if (point_set_) {
			delete point_set_;
			point_set_ = nil;
		}
		cameras_.clear();
		update_graphics();
		return true;
	}

	return false;
}

bool PaintCanvas::loadProject(const QString& file) {
	if (!project_)
		project_ = new Project;

	if (project_ && project_->load(file.toStdString())) {
		if (point_set_) {
			delete point_set_;
			point_set_ = nil;
		}
		cameras_.clear();
		update_graphics();
		return true;
	} 

	return false;
}


bool PaintCanvas::saveProject(const QString& file) {
	if (project_)
		return project_->save(file.toStdString());

	return false;
}


void PaintCanvas::imageMatching() {
	if (!project_) {
		LOG(WARNING) << "empty project" << std::endl;
		return;
	}
	if (!project_->is_valid()) {
		LOG(WARNING) << "invalid project" << std::endl;
		return;
	}
	
	ImageMatching matching(project_);
	matching.apply();

	// some image could be ignored
	main_window_->listWidgetImages->updateImageList();
	main_window_->saveProject();
}


void PaintCanvas::sparseReconstruction() {
	if (!project_) {
		LOG(WARNING) << "empty project" << std::endl;
		return;
	}
	if (!project_->is_valid()) {
		LOG(WARNING) << "invalid project" << std::endl;
		return;
	}

	if (!point_set_) {
		PointCloud* pset = new PointCloud;
		setPointSet(pset, false);
	}

	SparseReconstruction sparse(project_);
	sparse.apply(point_set_);
	update_all();
}


void PaintCanvas::denseReconstruction() {
	if (!project_) {
		LOG(WARNING) << "empty project" << std::endl;
		return;
	}
	if (!project_->is_valid()) {
		LOG(WARNING) << "invalid project" << std::endl;
		return;
	}

	if (!point_set_) {
		PointCloud* pset = new PointCloud;
		setPointSet(pset, false);
	}

	DenseReconstruction dense(project_);
	dense.run_pmvs(point_set_);
	fitScreen();
	update_all();
}


void PaintCanvas::test() {


}

