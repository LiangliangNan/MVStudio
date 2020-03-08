#ifndef MAIN_WINDOW_H
#define MAIN_WINDOW_H

#include <QMainWindow>

#include "ui_main_window.h"
#include "easy3d/core/types.h"

class QLabel;
class PaintCanvas;
class QProgressBar;

class MainWindow 
	: public QMainWindow
	, public Ui::MVStudioClass
{
	Q_OBJECT

public:
	MainWindow(QWidget *parent = 0, Qt::WindowFlags flags = 0);
	~MainWindow();

	static std::string title() { return "MainWindow"; }

	PaintCanvas* canvas() { return mainCanvas_; }

	void showCoordinateUnderMouse(const easy3d::dvec3 & p, bool found);

public Q_SLOTS:
	bool newProject();
	bool openProject();

	bool saveProject();
	bool saveProjectAs();

	bool openPointCloud();
	bool savePointCloud();

	void openRecentFile();
	void updateStatusBar();

	void setFullScreen(bool b);

	void about();
	void aboutQt();

	void setBackgroundColor();

private:
	void createMenus();
	void createActions(); 
	void createStatusBar();

	void createActionsForFileMenu();
	void createActionsForViewMenu();
	void createActionForReconstructionMenu();

	void readSettings();
	void writeSettings();

	void closeEvent(QCloseEvent *event);

	bool okToContinue();
	
	bool doOpenPointCloud(const QString &fileName);
	bool doSavePointCloud(const QString &fileName);

	void setCurrentFile(const QString &fileName);
	
	void updateRecentFileActions();
	QString strippedName(const QString &fullFileName);

protected:
	void dragEnterEvent(QDragEnterEvent *e);
	void dropEvent(QDropEvent *e);

private:
	PaintCanvas*	mainCanvas_;

	QStringList		recentFiles_;
	QString			curProjectFile_;
	QString			workingDirectory_;

	QLabel *statusLabel_,
		*coordinateUnderMouseLabel_,
		*numVerticesLabel_;

	enum { MaxRecentFiles = 5 };
	QAction *actionsRecentFile[MaxRecentFiles],
		*actionSeparator;
};

#endif // TESTQGLVIEWER_H
