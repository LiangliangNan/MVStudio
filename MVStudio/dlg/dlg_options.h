#ifndef DLG_OPTIONS_H
#define DLG_OPTIONS_H

#include <QDialog>
#include "ui_dlg_options.h"



class PaintCanvas;

class DlgOptions : public QDialog, public Ui::DlgOptions
{
	Q_OBJECT

public:
	DlgOptions(PaintCanvas* cvs, QWidget* parent);
	~DlgOptions();

private Q_SLOTS:
	void showHint(bool);
	void resetParameters();
	void reconstruct();
	void trim();

private:
	PaintCanvas* canvas_;

	int			default_octree_depth_;
	QString		default_samples_per_node_;
	QString		default_trim_value_;
	QString		default_area_ratio_;

	std::string density_attr_name_;

	//////////////////////////////////////////////////////////////////////////

	std::string cur_model_name_;
};

#endif 
