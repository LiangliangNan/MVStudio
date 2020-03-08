#ifndef PAINTCANVAS_H
#define PAINTCANVAS_H

#include "GL/glew.h"
#include "../3rd_party/QGLViewer/qglviewer.h"
#include "../math/math_types.h"
#include "../algo/sparse_reconstruction.h"
#include "../algo/project.h"
#include "../sfm/camera.h"
#include "../basic/canvas.h"

class PointSet;
class MainWindow;
class PointSetRender;


class PaintCanvas : public QGLViewer, public Canvas
{
	Q_OBJECT

public:
	PaintCanvas(QWidget *parent);
	~PaintCanvas();

public:
	std::string title() const { return "PaintCanvas"; }

	//////////////////////////////////////////////////////////////////////////
	void fit();
	void update_graphics();
	void update_all();

	//////////////////////////////////////////////////////////////////////////

	MainWindow* mainWindow() const { return main_window_; }
	void  setLightPosition(const vec3d& pos) ;
	const vec3d& lightPosition() const { return light_pos_; }

	vec2d projectionOf(const vec3d& p);          // point to screen
	vec3d unProjectionOf(double winx, double winy, double winz);  // screen to point

	//////////////////////////////////////////////////////////////////////////

	PointSet* pointSet() const { return point_set_; }
	void setPointSet(PointSet* pset, bool fit);

	bool creatProject(const QString& project_file);
	bool loadProject(const QString& file);
	bool saveProject(const QString& file);

	Project* project() const { return project_; }

protected:
	virtual void draw();
	virtual void init();

	// Mouse events functions
	virtual void mousePressEvent(QMouseEvent *e);
	virtual void mouseMoveEvent(QMouseEvent *e);
	virtual void mouseReleaseEvent(QMouseEvent *e);

	// Keyboard events functions
	virtual void keyPressEvent(QKeyEvent *e);

public Q_SLOTS:
	void fitScreen() ;

	void snapshotScreen();

	void copyCamera();
	void pasteCamera();

	void setLighting(bool b);

	void setProjectionMode(bool);
	void showCoordinateSystem(bool);

	void showPointAsSphere(bool b);
	void increasePointSize();
	void decreasePointSize();
	
	//////////////////////////////////////////////////////////////////////////

	void imageMatching();
	void sparseReconstruction();
	void denseReconstruction();

	void test();

private :
	void drawCornerAxis();

protected:
	MainWindow*	main_window_;
	vec3d		light_pos_;

	int		coord_system_region_size_;
	bool		show_coord_sys_;

	Project*	project_;

	PointSetRender*	render_;
	PointSet*		point_set_;
	std::vector<sfm::camera_params_t> cameras_;
};


#endif // PAINTCANVAS_H
