#ifndef FIGURE3DWINDOW_H_
#define FIGURE3DWINDOW_H_

#include <vector>

#include <QWidget>
#include <QTimer>

#include "ui_Figure3DWindow.h"

#include "global.h"
#include "Matrix2.h"
#include "XZ.h"
#include "XZValue.h"



namespace Lab {

class Figure3DWindow : public QWidget {
	Q_OBJECT
public:
	Figure3DWindow(QWidget *parent = 0);
	~Figure3DWindow();

	void updateData(double valueScale, const Matrix2<XZValue<float> >* gridData = 0, const std::vector<XZ<float> >* pointList = 0);

	// Will be effective after the next data update.
	void setVisualization(Figure::Visualization visualization);
	// Will be effective after the next data update.
	void setColormap(Figure::Colormap colormap);
private slots:
	void on_colormapComboBox_currentIndexChanged(int index);
	void on_minDecibelsComboBox_currentIndexChanged(const QString& text);
	void on_visualizationComboBox_currentIndexChanged(int index);
	void on_showInfoCheckBox_stateChanged(int state);
	void on_rotationCheckBox_stateChanged(int state);
//	void updateScanWidget();
private:
	Figure::Visualization visualization_;
	int minDecibels_;
	Figure::Colormap colormap_;
	Ui::Figure3DWindowClass ui_;
	Matrix2<XZValue<float>> gridData_;
	std::vector<XZ<float>> pointList_;
	double valueScale_;
};

} // namespace Lab

#endif /* FIGURE3DWINDOW_H_ */
