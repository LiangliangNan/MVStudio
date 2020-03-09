#ifndef MAIN_WINDOW_H
#define MAIN_WINDOW_H

#include <QMainWindow>

#include "ui_main_window.h"
#include "../math/math_types.h"
#include "../basic/logger.h"
#include "../basic/progress.h"


class QLabel;
class PaintCanvas;
class QProgressBar;

class MainWindow 
	: public QMainWindow
	, public LoggerClient
	, public ProgressClient
	, public Ui::MVStudioClass
{
	Q_OBJECT

public:
	MainWindow(QWidget *parent = 0, Qt::WindowFlags flags = 0);
	~MainWindow();

	static std::string title() { return "MainWindow"; }

	PaintCanvas* canvas() { return mainCanvas_; }

	//////////////////////////////////////////////////////////////////////////

	virtual void out_message(const std::string& msg) ;
	virtual void warn_message(const std::string& msg) ;
	virtual void err_message(const std::string& msg) ;
	virtual void status_message(const std::string& msg, int timeout) ;
	virtual void notify_progress(std::size_t value, bool show_text);

	//////////////////////////////////////////////////////////////////////////

public Q_SLOTS:
	bool newProject();
	bool openProject();

	bool saveProject();
	bool saveProjectAs();

	bool openPointCloud();
	bool savePointCloud();

	void saveCameraToFile();
	void restoreCameraFromFile();

	void openRecentFile();
	void updateStatusBar();

	void setFullScreen(bool b);

	void cancelTask();

	void about();
	void aboutQt();

	void setBackgroundColor();

private:
	void createMenus();
	void createActions(); 
	void createStatusBar();

	void createActionsForFileMenu();
	void createActionsForViewMenu();
	void createActionsForCameraMenu();
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

	QProgressBar*	progress_bar_;

	QLabel *statusLabel_,
		*numVerticesLabel_;

	enum { MaxRecentFiles = 5 };
	QAction *actionsRecentFile[MaxRecentFiles],
		*actionSeparator;
};

#endif // TESTQGLVIEWER_H
