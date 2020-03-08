#ifndef IMAGE_ITEM_H
#define IMAGE_ITEM_H

#include <QListWidgetItem>

struct ProjectImage;
class ImageItem : public QListWidgetItem
{
public:
	ImageItem(QListWidget *parent = 0) {
		static QIcon iconImage("Resources/building.png");
		QListWidgetItem::setIcon(iconImage);
	}
	~ImageItem() {}

	void set_image(const std::string& img) { image_file_ = img; }
	const std::string& image() const { return image_file_; }

private:
	std::string	image_file_;
};

#endif // OBJECT_ITEM_H
