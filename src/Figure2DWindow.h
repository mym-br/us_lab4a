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

#ifndef FIGURE2DWINDOW_H_
#define FIGURE2DWINDOW_H_

#include <vector>

#include "ui_Figure2DWindow.h"



QT_FORWARD_DECLARE_CLASS(QMouseEvent)
class QCPAxis;

namespace Lab {

class Figure2DWindow : public QWidget {
	Q_OBJECT
public:
	Figure2DWindow(QWidget* parent=0);
	~Figure2DWindow();

	void updateData(const std::vector<double>& xList, const std::vector<double>& yList, bool markPoints);

private slots:
	void on_figure2DWidget_axisClick(QCPAxis* axis, QCPAxis::SelectablePart part, QMouseEvent* event);
	void on_figure2DWidget_selectionChangedByUser();
	void on_figure2DWidget_mouseDoubleClick(QMouseEvent* event);
private:
	Ui::Figure2DWindowClass ui_;
};

} // namespace Lab

#endif /* FIGURE2DWINDOW_H_ */
