/***************************************************************************
 *  Copyright 2018 Marcelo Y. Matuda                                       *
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

#include "MultiLayer3DWindow.h"



namespace Lab {

MultiLayer3DWindow::MultiLayer3DWindow(QWidget* parent)
		: QWidget(parent)
{
	ui_.setupUi(this);
}

void
MultiLayer3DWindow::updateData(const std::vector<XYZValue<float>>& pointArray, const std::vector<unsigned int>& indexArray)
{
	ui_.oglMultiLayerWidget->updateData(pointArray, indexArray);
}

} // namespace Lab
