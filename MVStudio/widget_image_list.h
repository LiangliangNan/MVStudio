#pragma once

#include <QListWidget>


class MainWindow;

class WidgetImageList : public QListWidget
{
	Q_OBJECT

public:
	WidgetImageList(QWidget * parent = 0);
	~WidgetImageList();
	
	void setMainWindow(MainWindow* w) { mainWindow_ = w; }

	void updateImageList();

private Q_SLOTS:
	void ignoreSelectedImages();
	void useSelectedImages();
	void removeSelectedImages();
	void addImages();

	void showContextMenu(const QPoint& p);

private:
	MainWindow*	mainWindow_;

	QMenu*		popupMenu_;

	QAction*	actionAddImages_;
	QAction*	actionRemoveSelectedImages_;

 	QAction*	actionIgnoreSelectedImages_;
 	QAction*	actionUseSelectedImages_;
};

