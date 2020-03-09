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

#include "../basic/logger.h"
#include "../pointset/point_set.h"
#include "../pointset/point_set_io.h"


MainWindow::MainWindow(QWidget *parent, Qt::WindowFlags flags)
: QMainWindow(parent, flags)
, workingDirectory_(".")
{	
	setupUi(this);
	listWidgetImages->setMainWindow(this);

	//////////////////////////////////////////////////////////////////////////

    Logger::instance()->set_value(Logger::LOG_REGISTER_FEATURES, "*"); // log everything
	Logger::instance()->set_value(Logger::LOG_FILE_NAME, "MVStudio.log");
	Logger::instance()->register_client(this);
	Progress::instance()->set_client(this);

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

#ifndef NDEBUG
    setWindowTitle("MVStudio (Debug Version)");
#else
    setWindowTitle("MVStudio");
#endif

	setAcceptDrops(true);

	setContextMenuPolicy(Qt::CustomContextMenu);
	setWindowIcon(QIcon(":/Resources/MVStudio.png"));

	//setFixedSize(1200, 900);
	setFocusPolicy(Qt::StrongFocus);
	setWindowState(Qt::WindowMaximized);
}

MainWindow::~MainWindow()
{
	Progress::instance()->set_client(nil) ;
	Logger::instance()->unregister_client(this);
	Logger::terminate();
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

	connect(actionCombinePointCloudFiles, SIGNAL(triggered()), this, SLOT(combinePointCloudFiles()));

	connect(actionExit, SIGNAL(triggered()), this, SLOT(close()));
	actionExit->setShortcut(QString("Ctrl+Q"));
}

void MainWindow::createActionsForViewMenu() {
	connect(actionPerspectiveOrthographic, SIGNAL(toggled(bool)), mainCanvas_, SLOT(setProjectionMode(bool)));
	mainCanvas_->setProjectionMode(true);

	menuView->insertSeparator(actionFullScreen);

	connect(actionFullScreen, SIGNAL(toggled(bool)), this, SLOT(setFullScreen(bool)));
	connect(actionFitWindow, SIGNAL(triggered()), mainCanvas_, SLOT(fitScreen()));
	connect(actionShowCoordinateSystem, SIGNAL(toggled(bool)), mainCanvas_, SLOT(showCoordinateSystem(bool)));

	connect(actionLighting, SIGNAL(toggled(bool)), mainCanvas_, SLOT(setLighting(bool)));

	connect(actionSetBackgroundColor, SIGNAL(triggered()), this, SLOT(setBackgroundColor()));

	connect(actionIncreasePointSize, SIGNAL(triggered()), mainCanvas_, SLOT(increasePointSize()));
	connect(actionDecreasePointSize, SIGNAL(triggered()), mainCanvas_, SLOT(decreasePointSize()));
}


void MainWindow::createActionsForCameraMenu() {
	connect(actionCopyCamera, SIGNAL(triggered()), mainCanvas_, SLOT(copyCamera()));
	connect(actionPasteCamera, SIGNAL(triggered()), mainCanvas_, SLOT(pasteCamera()));

	connect(actionSaveCameraToFile, SIGNAL(triggered()), this, SLOT(saveCameraToFile()));
	connect(actionRestoreCameraFromFile, SIGNAL(triggered()), this, SLOT(restoreCameraFromFile()));

	connect(actionSnapshot, SIGNAL(triggered()), mainCanvas_, SLOT(snapshotScreen()));
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

	// camera menu
	createActionsForCameraMenu();

	// reconstruction menu
	createActionForReconstructionMenu();

	// about menu
	connect(actionHelp, SIGNAL(triggered()), mainCanvas_, SLOT(help()));
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

 	QPushButton* cancelButton = new QPushButton;
 	cancelButton->setIcon(QIcon(":/Resources/cancel.png"));
 	cancelButton->setFixedSize(20, 20);
 	statusBar()->addPermanentWidget(cancelButton, 1);
 	connect(cancelButton, SIGNAL(pressed()), this,  SLOT(cancelTask()));

// 	int s = dockWidgetOutput->sizeHint().width();
// 	std::cout << "dockWidgetOutput->sizeHint().width(): " << s << std::endl;
// 	s = dockWidgetImages->sizeHint().width();
// 	std::cout << "dockWidgetImages->sizeHint().width(): " << s << std::endl;

	progress_bar_ = new QProgressBar;
	progress_bar_->setAlignment(Qt::AlignRight | Qt::AlignVCenter);
	progress_bar_->setFixedWidth(dockWidgetOutput->sizeHint().width());
	statusBar()->addPermanentWidget(progress_bar_, 1);

	//////////////////////////////////////////////////////////////////////////

	updateStatusBar();
}


void MainWindow::updateStatusBar()
{
	QString vertices("");

	PointSet* pset = mainCanvas_->pointSet();
	if (pset) {
		int num_vertices = pset->num_points();
		vertices = QString("#vertices: %1").arg(num_vertices);
	} 

	numVerticesLabel_->setText( vertices );
}


void MainWindow::setBackgroundColor() {
	QColor color = QColorDialog::getColor(Qt::black, this);
	if (color.isValid() && mainCanvas_) {
		mainCanvas_->makeCurrent();
		mainCanvas_->setBackgroundColor(color);

		mainCanvas_->update_graphics();
	}
}

void MainWindow::showCoordinateUnderMouse(const vec3d& p, bool found) {
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


void MainWindow::out_message(const std::string& msg) {
	plainTextEditOutput->moveCursor(QTextCursor::End);
	plainTextEditOutput->insertPlainText(QString::fromStdString(msg));
	plainTextEditOutput->repaint();
	plainTextEditOutput->update();
}


void MainWindow::warn_message(const std::string& msg) {
	plainTextEditOutput->moveCursor(QTextCursor::End);
	plainTextEditOutput->insertPlainText(QString::fromStdString(msg));
	plainTextEditOutput->repaint();
	plainTextEditOutput->update();
}


void MainWindow::err_message(const std::string& msg) {
	plainTextEditOutput->moveCursor(QTextCursor::End);
	plainTextEditOutput->insertPlainText(QString::fromStdString(msg));
	plainTextEditOutput->repaint();
	plainTextEditOutput->update();
}



void MainWindow::status_message(const std::string& msg, int timeout) {
	statusBar()->showMessage(QString::fromStdString(msg), timeout);
}


void MainWindow::notify_progress(std::size_t value, bool show_text) {
	progress_bar_->setValue(int(value));
	progress_bar_->setTextVisible(show_text);
	mainCanvas_->update_all();
}



void MainWindow::cancelTask() {
	int value = progress_bar_->value() ;

	Progress::instance()->cancel();
	progress_bar_->reset();
	progress_bar_->setTextVisible(false);
	mainCanvas_->update_all();

	if (value != -1 && value != 0)
		Logger::warn("MVStudio") << "task canceled" << std::endl;
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

#ifndef NDEBUG
    setWindowTitle(tr("%1[*] - %2").arg(shownName).arg(tr("MVStudio (Debug Version)")));
#else
    setWindowTitle(tr("%1[*] - %2").arg(shownName).arg(tr("MVStudio")));
#endif
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
			Logger::out(title()) << "project loaded" << std::endl;
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
		Logger::warn(title()) << "project already opened" << std::endl;
		return false;
	}

	if (mainCanvas_->loadProject(fileName)) {
		listWidgetImages->updateImageList();
		setCurrentFile(fileName);
		Logger::out(title()) << "project loaded" << std::endl;
		return true;
	}
	else
		return false;
}


bool MainWindow::saveProject() {
	if (!FileUtils::is_file(curProjectFile_.toStdString()))  {
		return saveProjectAs();
	}
	else {
		if (mainCanvas_->saveProject(curProjectFile_)) {
			Logger::out(title()) << "project saved" << std::endl;
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
		Logger::out(title()) << "project saved" << std::endl;
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
	PointSet* pset = PointSetIO::read(name);

	if (pset) {
		mainCanvas_->setPointSet(pset, true);
		updateStatusBar();
		status_message("File loaded", 500);

		return true;
	} 	
	else {	
		status_message("Open failed", 500);
		return false;
	} 
}

bool MainWindow::doSavePointCloud(const QString &fileName)
{
	PointSet* pset = mainCanvas_->pointSet();
	if (pset == nil)
		return false;

	std::string name = fileName.toStdString();
	bool success = PointSetIO::save(name, pset);

	if (success) {
		status_message("File saved", 500);
		return true;
	} else {
		status_message("Saving failed", 500);
		return false;
	}

	return true;
}


void MainWindow::saveCameraToFile() {
	QString fileName = QFileDialog::getSaveFileName(this,
		tr("Save Camera State to File"), curProjectFile_ + "_camera.cs",
		tr("Camera state (*.cs)"));

	if (fileName.isEmpty())
		return ;

	if (mainCanvas_) {
		mainCanvas_->setStateFileName(fileName);
		mainCanvas_->saveStateToFile();

		mainCanvas_->setStateFileName("");
	} 
}


void MainWindow::restoreCameraFromFile() {
	QString fileName = QFileDialog::getOpenFileName(this,
		tr("Restore Camera State from File"), workingDirectory_,
		tr("Camera state (*.cs)"));

	if (fileName.isEmpty())
		return ;

	if (mainCanvas_) {
		mainCanvas_->setStateFileName(fileName);
		mainCanvas_->restoreStateFromFile();

		mainCanvas_->setStateFileName("");
		mainCanvas_->update_graphics();
	} 
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


void MainWindow::combinePointCloudFiles() {
	QStringList fileNames = QFileDialog::getOpenFileNames(this,
		tr("Please select files"), ".",
		tr("Points (*.ply *.pnc *.bpnc)")
		);

	if (fileNames.isEmpty())
		return;

	QFileInfo info(fileNames[0]);
	QString output_file = info.path() + "/combined_file.bpnc";
	// open file
	std::ofstream output(output_file.toStdString().c_str(), std::fstream::binary);
	if (output.fail()) {
		Logger::err("MainWindow") << "could not open file\'" << output_file.toStdString() << "\'" << std::endl;
		return;
	}

	size_t total_num = 0;
	ProgressLogger progress(fileNames.size());
	for (int i = 0; i < fileNames.size(); ++i) {
		if (progress.is_canceled())
			break;

		const std::string& in_file_name = fileNames[i].toStdString();
		const std::string& simple_name = FileUtils::simple_name(in_file_name);
		Logger::out("MainWindow") << "reading file " << i + 1 << "/" << fileNames.size() << ": " << simple_name << std::endl;
		PointSet* pset = PointSetIO::read(in_file_name);
		if (!pset) 
			continue;
		else {
			Logger::out("MainWindow") << "saving file " << i + 1 << "/" << fileNames.size() << std::endl;

			int num = pset->num_points();
			float* points = &(const_cast<PointSet*>(pset)->points()[0].x);
			float* normals = &(const_cast<PointSet*>(pset)->normals()[0].x);
			float* colors = &(const_cast<PointSet*>(pset)->colors()[0].x);

			for (int i = 0; i < num; ++i) {
				float* pt = points + (i * 3);	output.write((char*)pt, 12);
				float* nm = normals + (i * 3);	output.write((char*)nm, 12);
				float* cl = colors + (i * 3);	output.write((char*)cl, 12);
			}
		}

		if (pset)
			delete pset;

		progress.next();
	}

	Logger::out("MainWindow") << "done! Total number: " << total_num << std::endl;
}

