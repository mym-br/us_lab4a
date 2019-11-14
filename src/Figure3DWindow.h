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

#include "ui_Figure3DWindow.h"

#include "Colormap.h"
#include "Matrix.h"
#include "Visualization.h"
#include "XYZ.h"
#include "XYZValue.h"



namespace Lab {

class Figure3DWindow : public QWidget {
	Q_OBJECT
public:
	Figure3DWindow(QWidget* parent=nullptr);
	~Figure3DWindow();

	// *gridData will contain old data after the call.
	void updateData(float valueScale,
			Matrix<XYZValue<float>>* gridData=nullptr,
			const std::vector<XYZ<float>>* pointList=nullptr);
	void setVisualization(Visualization::Value visualization);
	void setColormap(Colormap::Id colormap);
protected:
	virtual void closeEvent(QCloseEvent* event);
private slots:
	void on_colormapComboBox_currentIndexChanged(int index);
	void on_colormapComboBox_activated(int index); // caused by user interaction, called after on_colormapComboBox_currentIndexChanged
	void on_minDecibelsComboBox_currentIndexChanged(const QString& text);
	void on_visualizationComboBox_currentIndexChanged(int index);
	void on_visualizationComboBox_activated(int index); // caused by user interaction, called after on_visualizationComboBox_currentIndexChanged
	void on_showInfoCheckBox_stateChanged(int state);
	void on_rotationCheckBox_stateChanged(int state);
private:
	Ui::Figure3DWindowClass ui_;
};

} // namespace Lab

#endif /* FIGURE3DWINDOW_H_ */
