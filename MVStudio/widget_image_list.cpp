#include "widget_image_list.h"
#include "main_window.h"
#include "paint_canvas.h"
#include "image_item.h"

#include <easy3d/util/file_system.h>

#include <QMenu>
#include <QFileDialog>


WidgetImageList::WidgetImageList(QWidget * parent)
	: QListWidget(parent)
	, mainWindow_(nil)
{
 	actionIgnoreSelectedImages_ = new QAction("Ignore", this);
 	connect(actionIgnoreSelectedImages_, SIGNAL(triggered()), this, SLOT(ignoreSelectedImages()));
 	
 	actionUseSelectedImages_ = new QAction("Use", this);
 	connect(actionUseSelectedImages_, SIGNAL(triggered()), this, SLOT(useSelectedImages()));
	
	actionRemoveSelectedImages_ = new QAction("Remove", this);
	connect(actionRemoveSelectedImages_, SIGNAL(triggered()), this, SLOT(removeSelectedImages()));

	actionAddImages_ = new QAction("Add Images", this);
	connect(actionAddImages_, SIGNAL(triggered()), this, SLOT(addImages()));

	popupMenu_ = new QMenu(this);
	connect(this, SIGNAL(customContextMenuRequested(const QPoint&)), this, SLOT(showContextMenu(const QPoint&)));
}


WidgetImageList::~WidgetImageList()
{
}


void WidgetImageList::showContextMenu(const QPoint& p) {
	popupMenu_->clear();  // I want to customize the menu list
	
	if (selectedItems().size() > 0) {
		popupMenu_->addAction(actionIgnoreSelectedImages_);
		popupMenu_->addAction(actionUseSelectedImages_);
		popupMenu_->addSeparator();
		popupMenu_->addAction(actionRemoveSelectedImages_);
	}
	else if (count() == 0) {
		popupMenu_->addAction(actionAddImages_);
	}
	
	popupMenu_->popup(mapToGlobal(p));
}


void WidgetImageList::updateImageList() {
	const std::vector<ProjectImage>& images = mainWindow_->canvas()->project()->images;
	for (std::size_t i = 0; i < images.size(); ++i) {
		// make sure the 'item' widget exists
		ImageItem* item = dynamic_cast<ImageItem*>(this->item(i));
		if (!item) {
			item = new ImageItem;
			addItem(item);
		}

		item->set_image(images[i].file);
		const std::string& name = file_system::simple_name(images[i].file);
		if (images[i].ignored)
			item->setData(Qt::DisplayRole, QString::fromStdString(name + "	ignored"));
		else
			item->setData(Qt::DisplayRole, QString::fromStdString(name));
	}

	// remove redundant items
	int num_items = count();
	if (num_items > images.size()) {
		int delta = int(num_items - images.size());
		for (int i = 0; i < delta; ++i) {
			int idx = count() - 1;
			if (idx >= 0) {
				QListWidgetItem* item = takeItem(idx);
				delete item;
			}
			else
				ogf_assert_not_reached;
		}
	}

	assert(count() == images.size());

	QString title = QString("Images (Total: %1)").arg(num_items);
	mainWindow_->dockWidgetImages->setWindowTitle(title);
	update();
}

void WidgetImageList::addImages() {
	Project* project = mainWindow_->canvas()->project();
	if (!project) {
		Logger::warn("MainWindow") << "please create a project first" << std::endl;
		return;
	}

	QString imageDir = QFileDialog::getExistingDirectory(this,
		tr("Choose the folder containing the images"), QString::fromStdString(project->project_dir)
		);

	if (imageDir.isEmpty())
		return;

	project->add_images(imageDir.toStdString());
	updateImageList();
	mainWindow_->setWindowModified(true);
}

void WidgetImageList::removeSelectedImages() {
	if (!mainWindow_)
		return;

	std::vector<std::string> images;
	QList<QListWidgetItem*> itmes = selectedItems();
	for (int i = 0; i < itmes.size(); ++i) {
		std::string file = dynamic_cast<ImageItem*>(itmes[i])->image();
		images.push_back(file);
	}

	std::vector<ProjectImage>& images_files = mainWindow_->canvas()->project()->images;
	for (std::size_t i = 0; i < images.size(); ++i) {
		std::vector<ProjectImage>::iterator pos = images_files.begin();
		for (; pos != images_files.end(); ++pos) {
			if (pos->file == images[i]) {
				images_files.erase(pos);
				break;
			}
		}
	}

	mainWindow_->setWindowModified(true);
	updateImageList();
}



void WidgetImageList::ignoreSelectedImages() {
	if (!mainWindow_)
		return;

	std::vector<std::string> images;
	QList<QListWidgetItem*> itmes = selectedItems();
	for (int i = 0; i < itmes.size(); ++i) {
		std::string file = dynamic_cast<ImageItem*>(itmes[i])->image();
		images.push_back(file);
	}

	std::vector<ProjectImage>& images_files = mainWindow_->canvas()->project()->images;
	for (std::size_t i = 0; i < images.size(); ++i) {
		std::vector<ProjectImage>::iterator pos = images_files.begin();
		for (; pos != images_files.end(); ++pos) {
			if (pos->file == images[i]) {
				pos->ignored = true;
				break;
			}
		}
	}

	mainWindow_->setWindowModified(true);
	updateImageList();
}



void WidgetImageList::useSelectedImages() {
	if (!mainWindow_)
		return;

	std::vector<std::string> images;
	QList<QListWidgetItem*> itmes = selectedItems();
	for (int i = 0; i < itmes.size(); ++i) {
		std::string file = dynamic_cast<ImageItem*>(itmes[i])->image();
		images.push_back(file);
	}

	std::vector<ProjectImage>& images_files = mainWindow_->canvas()->project()->images;
	for (std::size_t i = 0; i < images.size(); ++i) {
		std::vector<ProjectImage>::iterator pos = images_files.begin();
		for (; pos != images_files.end(); ++pos) {
			if (pos->file == images[i]) {
				pos->ignored = false;
				break;
			}
		}
	}

	mainWindow_->setWindowModified(true);
	updateImageList();
}