/***************************************************************************
 *  Copyright 2014, 2017, 2018 Marcelo Y. Matuda                           *
 *                                                                         *
 *  This program is free software: you can redistribute it and/or modify   *
 *  it under the terms of the GNU General Public License as published by   *
 *  the Free Software Foundation, either version 3 of the License, or      *
 *  (at your option) any later version.                                    *
 *                                                                         *
 *  This program is distributed in the hope that it will be useful,        *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of         *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
 *  GNU General Public License for more details.                           *
 *                                                                         *
 *  You should have received a copy of the GNU General Public License      *
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.  *
 ***************************************************************************/

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
