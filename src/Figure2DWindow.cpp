/***************************************************************************
 *  Copyright 2014, 2017, 2018, 2019 Marcelo Y. Matuda                     *
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

#include "Figure2DWindow.h"

#include "Figure2DWidget.h"

namespace Lab {

Figure2DWindow::Figure2DWindow(QWidget* parent)
		: QWidget(parent)
{
	ui_.setupUi(this);
}

void
Figure2DWindow::updateData(const std::vector<double>& xList, const std::vector<double>& yList, bool markPoints)
{
	ui_.figure2DWidget->updateData(xList, yList);
	if (markPoints) {
		ui_.figure2DWidget->enablePointMarker(true);
	}
}

} // namespace Lab
