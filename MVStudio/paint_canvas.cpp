#include "paint_canvas.h"
#include "main_window.h"

#include "../libs/opengl/opengl_info.h"

#include "../libs/pointset/point_set.h"
#include "../libs/pointset/point_set_io.h"
#include "../libs/pointset/point_set_render.h"
#include "../libs/algo/image_matching.h"
#include "../libs/algo/sparse_reconstruction.h"
#include "../libs/algo/dense_reconstruction.h"

#include "../3rd_party/QGLViewer/QGLViewer/manipulatedCameraFrame.h"

#include <QMouseEvent>
#include <QMessageBox>

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

	light_pos_ = vec3d(0.27f, 0.27f, 0.92f);

	//////////////////////////////////////////////////////////////////////////

	// Move camera according to viewer type (on X, Y or Z axis)
	// e.g. 'GL_FRONT' equals to 
	camera()->setViewDirection(qglviewer::Vec(0.0, 1.0, 0.0));

	camera()->lookAt(sceneCenter());
	camera()->setType(qglviewer::Camera::PERSPECTIVE);
	camera()->showEntireScene();
    ogf_check_gl;
}


PaintCanvas::~PaintCanvas() {
	delete render_;
}


void PaintCanvas::init()
{
    ogf_check_gl;
	Logger::out(title()) << "initializing..." << std::endl ;

 	//////////////////////////////////////////////////////////////////////////
 	GLenum err = glewInit();
 	if (GLEW_OK != err) {
 		// Problem: glewInit failed, something is seriously wrong.
 		Logger::err(title()) << glewGetErrorString(err) << std::endl ;
 	}
	Logger::out(title()) << "Using GLEW: " << GLInfo::glew_version() << std::endl ;
	Logger::out(title()) << "GL Vendor: " << GLInfo::gl_vendor() << std::endl ;
	Logger::out(title()) << "GL Renderer: " << GLInfo::gl_renderer() << std::endl ;
	Logger::out(title()) << "GL Version: " << GLInfo::gl_version() << std::endl ;

	Logger::out(title()) << "GLSL Version: " << GLInfo::glsl_version() << std::endl ;
	//Logger::out(title()) << "GL_Extensions: " << GLInfo::gl_extensions() << std::endl;

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
//	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);  //GL_ONE_MINUS_DST_ALPHA

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
    ogf_check_gl;
}

void PaintCanvas::setLightPosition(const vec3d& pos) { 
	light_pos_ = pos;

	makeCurrent();

	glPushMatrix();
	glLoadIdentity();
	// 0 = infinite light 
	// 1 = local light
	float p[] = {static_cast<float>(light_pos_.x), static_cast<float>(light_pos_.y), static_cast<float>(light_pos_.z), 0.0f };
	glLightfv(GL_LIGHT0, GL_POSITION, p);
	glPopMatrix();	

    doneCurrent();
	update_graphics();
    ogf_check_gl;
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
    ogf_check_gl;
	if (point_set_) {
		glDisable(GL_MULTISAMPLE);
		render_->draw();
	}

    if (show_coord_sys_)  {
        glEnable(GL_MULTISAMPLE);
        drawCornerAxis();
    }

    ogf_check_gl;
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
#if QT_VERSION >= QT_VERSION_CHECK(5, 14, 0)
	QStringList list = str.split(" ", Qt::SkipEmptyParts);
#else
    QStringList list = str.split(" ", QString::SkipEmptyParts);
#endif
	if(list.size() != 7) {
		Logger::err(title()) << "camera not available in clipboard" << std::endl;
		return;
	}

	float vec[3];
	for(int i=0; i<3; ++i){
		bool ok;
		vec[i] = list[i].toFloat(&ok);
		if(!ok) {
			Logger::err(title()) << "camera not available in clipboard" << std::endl;
			return;
		}
	}

	double orient[4];
	for(int i=0; i< 4; ++i) {
		bool ok;
		orient[i] = list[i+3].toDouble(&ok);
		if(!ok) {
			Logger::err(title()) << "camera not available in clipboard" << std::endl;
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
    ogf_check_gl;
}

void PaintCanvas::fit() {
	fitScreen();
    ogf_check_gl;
}

void PaintCanvas::update_graphics() {
	update();
    ogf_check_gl;
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
    ogf_check_gl;
}


void PaintCanvas::setProjectionMode(bool b) {
	if (b)
		camera()->setType(qglviewer::Camera::PERSPECTIVE);
	else
		camera()->setType(qglviewer::Camera::ORTHOGRAPHIC);

	update_graphics();
    ogf_check_gl;
}

void PaintCanvas::showCoordinateSystem(bool b) { 
	show_coord_sys_= b;
	update_graphics();
    ogf_check_gl;
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
    ogf_check_gl;
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

	glMatrixMode(GL_PROJECTION);
	glPopMatrix();

	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();

	//////////////////////////////////////////////////////////////////////////

	// The viewport and the scissor are restored.
	glScissor(scissor[0], scissor[1], scissor[2], scissor[3]);
	glViewport(viewport[0], viewport[1], viewport[2], viewport[3]);
    ogf_check_gl;
}


vec2d PaintCanvas::projectionOf(const vec3d& p) {    // point to screen
	qglviewer::Vec v = camera()->projectedCoordinatesOf(qglviewer::Vec(p.x, p.y, p.z));
	return vec2d(v.x, v.y);
}

vec3d PaintCanvas::unProjectionOf(double winx, double winy, double winz) {  // screen to point	
	qglviewer::Vec v = camera()->unprojectedCoordinatesOf(qglviewer::Vec(winx, winy, winz));
	return vec3d(v.x, v.y, v.z);
}


void PaintCanvas::setPointSet(PointSet* pset, bool fit) {
	if (point_set_) {
		delete point_set_;
		render_->set_pointset(nil);
	}

	point_set_ = pset;
	point_set_->set_canvas(this);
	render_->set_pointset(point_set_);
	if (fit)
		fitScreen();
    ogf_check_gl;
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
    ogf_check_gl;
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
        ogf_check_gl;
		return true;
	}

    ogf_check_gl;
	return false;
}


bool PaintCanvas::saveProject(const QString& file) {
	if (project_)
		return project_->save(file.toStdString());

	return false;
}


void PaintCanvas::imageMatching() {
	if (!project_) {
		Logger::warn(title()) << "empty project" << std::endl;
		return;
	}
	if (!project_->is_valid()) {
		Logger::warn(title()) << "invalid project" << std::endl;
		return;
	}

    makeCurrent();
	ImageMatching matching(project_);
	matching.apply();
    doneCurrent();
    ogf_check_gl;
	// some image could be ignored
	main_window_->listWidgetImages->updateImageList();
	main_window_->saveProject();
}


void PaintCanvas::sparseReconstruction() {
	if (!project_) {
		Logger::warn(title()) << "empty project" << std::endl;
		return;
	}
	if (!project_->is_valid()) {
		Logger::warn(title()) << "invalid project" << std::endl;
		return;
	}

	if (!point_set_) {
		PointSet* pset = new PointSet;
		setPointSet(pset, false);
	}

	SparseReconstruction sparse(project_);
	sparse.apply(point_set_);
	update_all();
    ogf_check_gl;
}


void PaintCanvas::denseReconstruction() {
	if (!project_) {
		Logger::warn(title()) << "empty project" << std::endl;
		return;
	}
	if (!project_->is_valid()) {
		Logger::warn(title()) << "invalid project" << std::endl;
		return;
	}

	if (!point_set_) {
		PointSet* pset = new PointSet;
		setPointSet(pset, false);
	}

	DenseReconstruction dense(project_);
	dense.run_pmvs(point_set_);
	fitScreen();
	update_all();
    ogf_check_gl;
}
