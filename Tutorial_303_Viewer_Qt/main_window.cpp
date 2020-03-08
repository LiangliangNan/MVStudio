#include <QMessageBox>
#include <QFileDialog>
#include <QLabel>
#include <QStatusBar>
#include <QSettings>
#include <QCloseEvent>
#include <QPlainTextEdit>
#include <QGroupBox>
#include <QColorDialog>
#include <QProgressBar>
#include <QMimeData>
#include <QPushButton>

#include "main_window.h"
#include "paint_canvas.h"

#include <easy3d/util/logging.h>
#include <easy3d/util/file_system.h>
#include <easy3d/core/point_cloud.h>
#include <easy3d/fileio/point_cloud_io.h>

using namespace easy3d;

MainWindow::MainWindow(QWidget *parent, Qt::WindowFlags flags)
: QMainWindow(parent, flags)
, workingDirectory_(".")
{	
	setupUi(this);
	listWidgetImages->setMainWindow(this);

	//////////////////////////////////////////////////////////////////////////

	mainCanvas_ = new PaintCanvas(this);
	mainCanvas_->setAttribute(Qt::WA_MouseTracking);
	mainCanvas_->setMouseTracking(true);
	setCentralWidget(mainCanvas_);

	//////////////////////////////////////////////////////////////////////////

	createActions();
	createMenus();
	createStatusBar();

	readSettings();
	setWindowTitle("MVStudio");
	setAcceptDrops(true);

	setContextMenuPolicy(Qt::CustomContextMenu);
	setWindowIcon(QIcon(":/Resources/MVStudio.png"));

	//setFixedSize(1200, 900);
	setFocusPolicy(Qt::StrongFocus);
	setWindowState(Qt::WindowMaximized);
}

MainWindow::~MainWindow()
{
}


void MainWindow::dragEnterEvent(QDragEnterEvent *e) {
	if (e->mimeData()->hasUrls()) {
		e->acceptProposedAction();
	}
}

void MainWindow::dropEvent(QDropEvent *e) {
    foreach (const QUrl &url, e->mimeData()->urls()) {
        const QString &fileName = url.toLocalFile();
	  doOpenPointCloud(fileName);
    }
}


void MainWindow::createMenus() {
	actionSeparator = menuFile->addSeparator();

	QList<QAction*> actions;
	for (int i = 0; i < MaxRecentFiles; ++i)
		actions.push_back(actionsRecentFile[i]);
	menuFile->insertActions(actionExit, actions);
}

void MainWindow::createActionsForFileMenu() {
	//////////////////////////////////////////////////////////////////////////
	for (int i = 0; i < MaxRecentFiles; ++i) {
		actionsRecentFile[i] = new QAction(this);
		actionsRecentFile[i]->setVisible(false);
		connect(actionsRecentFile[i], SIGNAL(triggered()), this, SLOT(openRecentFile()));
	}

	connect(actionNewProject, SIGNAL(triggered()), this, SLOT(newProject()));
	connect(actionOpenProject, SIGNAL(triggered()), this, SLOT(openProject()));
	connect(actionSaveProject, SIGNAL(triggered()), this, SLOT(saveProject()));
	connect(actionSaveProjectAs, SIGNAL(triggered()), this, SLOT(saveProjectAs()));
	
	connect(actionAddImages, SIGNAL(triggered()), listWidgetImages, SLOT(addImages()));

	connect(actionOpenPointCloud, SIGNAL(triggered()), this, SLOT(openPointCloud()));
	connect(actionSavePointCloud, SIGNAL(triggered()), this, SLOT(savePointCloud()));

	connect(actionExit, SIGNAL(triggered()), this, SLOT(close()));
	actionExit->setShortcut(QString("Ctrl+Q"));
}

void MainWindow::createActionsForViewMenu() {
	connect(actionPerspectiveOrthographic, SIGNAL(toggled(bool)), mainCanvas_, SLOT(setProjectionMode(bool)));
	mainCanvas_->setProjectionMode(true);

	menuView->insertSeparator(actionFullScreen);

	connect(actionFullScreen, SIGNAL(toggled(bool)), this, SLOT(setFullScreen(bool)));
	connect(actionFitWindow, SIGNAL(triggered()), mainCanvas_, SLOT(fitScreen()));

	connect(actionLighting, SIGNAL(toggled(bool)), mainCanvas_, SLOT(setLighting(bool)));

	connect(actionPointAsSphere, SIGNAL(toggled(bool)), mainCanvas_, SLOT(showPointAsSphere(bool)));
	connect(actionIncreasePointSize, SIGNAL(triggered()), mainCanvas_, SLOT(increasePointSize()));
	connect(actionDecreasePointSize, SIGNAL(triggered()), mainCanvas_, SLOT(decreasePointSize()));
}


void MainWindow::createActionForReconstructionMenu() {
	connect(actionImageMatching, SIGNAL(triggered()), mainCanvas_, SLOT(imageMatching()));
	connect(actionSparseReconstruction, SIGNAL(triggered()), mainCanvas_, SLOT(sparseReconstruction()));
	connect(actionDenseReconstruction, SIGNAL(triggered()), mainCanvas_, SLOT(denseReconstruction()));
}


void MainWindow::createActions() {
	// file menu
	createActionsForFileMenu();

	// view menu
	createActionsForViewMenu();

	// reconstruction menu
	createActionForReconstructionMenu();

	// about menu
	connect(actionAbout, SIGNAL(triggered()), this, SLOT(about()));
	connect(actionAboutQt, SIGNAL(triggered()), this, SLOT(aboutQt()));
}


void MainWindow::createStatusBar()
{	
	statusLabel_ = new QLabel("Ready");
	statusLabel_->setAlignment(Qt::AlignLeft | Qt::AlignVCenter);
	statusBar()->addWidget(statusLabel_, 1);

	int length = 220;

	coordinateUnderMouseLabel_ = new QLabel("XYZ = [-, -, -]");
	coordinateUnderMouseLabel_->setFixedWidth(length);
	coordinateUnderMouseLabel_->setAlignment(Qt::AlignLeft | Qt::AlignVCenter);
	statusBar()->addWidget(coordinateUnderMouseLabel_, 1);

	QLabel* space1 = new QLabel;
	statusBar()->addWidget(space1, 1);

	numVerticesLabel_ = new QLabel;
	numVerticesLabel_->setFixedWidth(length);
	numVerticesLabel_->setAlignment(Qt::AlignLeft | Qt::AlignVCenter);
	statusBar()->addPermanentWidget(numVerticesLabel_, 1);

	QLabel* space2 = new QLabel;
	statusBar()->addWidget(space2, 1);

	//////////////////////////////////////////////////////////////////////////

	updateStatusBar();
}


void MainWindow::updateStatusBar()
{
	QString vertices("");

	PointCloud* pset = mainCanvas_->pointCloud();
	if (pset) {
		int num_vertices = pset->n_vertices();
		vertices = QString("#vertices: %1").arg(num_vertices);
	} 

	numVerticesLabel_->setText( vertices );
}


void MainWindow::setBackgroundColor() {
	QColor color = QColorDialog::getColor(Qt::black, this);
	if (color.isValid() && mainCanvas_) {
		mainCanvas_->makeCurrent();
		mainCanvas_->setBackgroundColor(vec4(color.redF(), color.greenF(), color.blueF(), 1.0f));

		mainCanvas_->update();
	}
}

void MainWindow::showCoordinateUnderMouse(const dvec3& p, bool found) {
	QString coordString = "XYZ = [-, -, -]";
	if (found)
		coordString = QString("XYZ = [%1, %2, %3]").arg(p.x).arg(p.y).arg(p.z);
	coordinateUnderMouseLabel_->setText(coordString);
}

void MainWindow::about()
{
	QString title = QMessageBox::tr(
		"<p align=\"center\"><span style=\"font-style:italic;\">I'm a good software, though I have defects.</span></p>"
	);

#if defined (ENV_32_BIT)
	title += QMessageBox::tr("<h3>MVStudio (32-bit)</h3>");
#elif defined (ENV_64_BIT)
	title += QMessageBox::tr("<h3>MVStudio (64-bit)</h3>");
#else
	title += QMessageBox::tr("<h3>MVStudio</h3>");
#endif

#ifndef NDEBUG
	title += QMessageBox::tr(" (Debug Version)");
#endif

	QString text = QMessageBox::tr(
		"<p><h4> Build %1</h4></p>"
		"<p>MVStudio is a prototype program for SfM and MVS integration.</p>"
		"<p>Liangliang Nan<br>"
		"<a href=\"mailto:liangliang.nan@gmail.com\">liangliang.nan@gmail.com</a><br>"
		"<a href=\"https://3d.bk.tudelft.nl/liangliang/\">https://3d.bk.tudelft.nl/liangliang/</a></p>"
	).arg("20180326");

	QMessageBox::about(this, "About MVStudio", title + text);
}

void MainWindow::aboutQt()
{
	QMessageBox::aboutQt(this, tr("About Qt"));
}


void MainWindow::setFullScreen(bool b) {
// 	canvas()->toggleFullScreen();
 	if (b)
 		showFullScreen();
 	else
 		showMaximized();
}


void MainWindow::readSettings()
{
	QSettings settings("liangliang.nan@gmail.com", "MVStudio");

	recentFiles_ = settings.value("recentFiles").toStringList();
	updateRecentFileActions();

	workingDirectory_ = settings.value("workingDirectory").toString();
}

void MainWindow::writeSettings()
{
	QSettings settings("liangliang.nan@gmail.com", "MVStudio");
	settings.setValue("recentFiles", recentFiles_);
	settings.setValue("workingDirectory", workingDirectory_);
}

void MainWindow::closeEvent(QCloseEvent *event)
{
	if (okToContinue()) {
		writeSettings();
		event->accept();
	} else {
		event->ignore();
	}
}

bool MainWindow::okToContinue()
{
	if (isWindowModified()) {
		int r = QMessageBox::warning(this, tr("MVStudio"),
			tr("Do you want to save your project?"),
			QMessageBox::Yes | QMessageBox::Default,
			QMessageBox::No,
			QMessageBox::Cancel | QMessageBox::Escape);
		if (r == QMessageBox::Yes)
			return saveProject();
		else if (r == QMessageBox::Cancel)
			return false;
	}
	return true;
}


void MainWindow::setCurrentFile(const QString &fileName)
{
	curProjectFile_ = fileName;
	workingDirectory_ = fileName.left(fileName.lastIndexOf("/") + 1); // path includes "/"

	setWindowModified(false);

	QString shownName = "Untitled";
	if (!curProjectFile_.isEmpty()) {
		shownName = strippedName(curProjectFile_);
		recentFiles_.removeAll(curProjectFile_);
		recentFiles_.prepend(curProjectFile_);
		updateRecentFileActions();
	}

	setWindowTitle(tr("%1[*] - %2").arg(shownName).arg(tr("MVStudio")));
}

void MainWindow::openRecentFile()
{
	if (!okToContinue())
		return;

	QAction *action = qobject_cast<QAction *>(sender());
	if (action) {
		QString filename(action->data().toString());
		if (mainCanvas_->loadProject(filename)) {
			listWidgetImages->updateImageList();
			setCurrentFile(filename);
			LOG(INFO) << "project loaded";
		}
	}
}


bool MainWindow::openPointCloud()
{
	QString fileName = QFileDialog::getOpenFileName(this,
		tr("Open point cloud"), workingDirectory_,
		tr(
		"Point cloud format (*.ply *.pnc *.bpnc)\n"
		"All format (*.*)")
		);

	if (fileName.isEmpty())
		return false;

	return doOpenPointCloud(fileName);
}

bool MainWindow::savePointCloud()
{
	QString fileName = QFileDialog::getSaveFileName(this,
		tr("Save point cloud"), workingDirectory_, 
		tr(
		"Point cloud format (*.ply *.pnc *.bpnc)\n"
		"All format (*.*)")
		);

	if (fileName.isEmpty())
		return false;

	return doSavePointCloud(fileName);
}


bool MainWindow::newProject() {
	if (!okToContinue())
		return false;

	QString fileName = QFileDialog::getSaveFileName(this,
		tr("Please specify where to save your project file"), workingDirectory_,
		tr("Project file (*.mvsproj)")
		);

	if (fileName.isEmpty())
		return false;

	if (mainCanvas_->creatProject(fileName)) {
		listWidgetImages->updateImageList();
		setCurrentFile(fileName);
		return true;
	}
	else
		return false;
}


bool MainWindow::openProject() {
	if (!okToContinue())
		return false;

	QString fileName = QFileDialog::getOpenFileName(this,
		tr("Open MVS project"), workingDirectory_,
		tr(
		"Project file (*.mvsproj)\n"
		"All format (*.*)")
		);

	if (fileName.isEmpty())
		return false;

	if (fileName == curProjectFile_) {
		LOG(WARNING) << "project already opened";
		return false;
	}

	if (mainCanvas_->loadProject(fileName)) {
		listWidgetImages->updateImageList();
		setCurrentFile(fileName);
		LOG(INFO) << "project loaded" << std::endl;
		return true;
	}
	else
		return false;
}


bool MainWindow::saveProject() {
	if (!file_system::is_file(curProjectFile_.toStdString()))  {
		return saveProjectAs();
	}
	else {
		if (mainCanvas_->saveProject(curProjectFile_)) {
			LOG(INFO) << "project saved" << std::endl;
			setWindowModified(false);
			return true;
		}
		return false;
	}
}


bool MainWindow::saveProjectAs() {
	QString fileName = QFileDialog::getSaveFileName(this,
		tr("Save project"), curProjectFile_,
		tr(
		"Project file (*.mvsproj)\n"
		"All format (*.*)")
		);

	if (fileName.isEmpty())
		return false;

	if (mainCanvas_->saveProject(fileName)) {
		LOG(INFO) << "project saved" << std::endl;
		setWindowModified(false);
		setCurrentFile(fileName);
		return true;
	}
	else
		return false;
}

bool MainWindow::doOpenPointCloud(const QString &fileName)
{
	std::string name = fileName.toStdString();
	PointCloud* pset = PointCloudIO::load(name);

	if (pset) {
		mainCanvas_->setPointSet(pset, true);
		updateStatusBar();
		LOG(INFO) << "File loaded";

		return true;
	} 	
	else {
		LOG(WARNING) << "Open failed";
		return false;
	} 
}

bool MainWindow::doSavePointCloud(const QString &fileName)
{
	PointCloud* pset = mainCanvas_->pointCloud();
	if (pset == nullptr)
		return false;

	std::string name = fileName.toStdString();
	bool success = PointCloudIO::save(name, pset);

	if (success) {
		LOG(INFO) << "File saved";
		return true;
	} else {
		LOG(WARNING) << "Save failed";
		return false;
	}

	return true;
}


void MainWindow::updateRecentFileActions()
{
	QMutableStringListIterator i(recentFiles_);
	while (i.hasNext()) {
		if (!QFile::exists(i.next()))
			i.remove();
	}

	for (int j = 0; j < MaxRecentFiles; ++j) {
		if (j < recentFiles_.count()) {
			QString text = tr("&%1 %2").arg(j + 1).arg(strippedName(recentFiles_[j]));
			actionsRecentFile[j]->setText(text);
			actionsRecentFile[j]->setData(recentFiles_[j]);
			actionsRecentFile[j]->setVisible(true);
		} else {
			actionsRecentFile[j]->setVisible(false);
		}
	}

	actionSeparator->setVisible(!recentFiles_.isEmpty());
}

QString MainWindow::strippedName(const QString &fullFileName)
{
	return QFileInfo(fullFileName).fileName();
}
