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

#include "Figure3DWindow.h"

#include <QCloseEvent>

namespace {

constexpr int MAX_MIN_DECIBELS = -1;
constexpr int MIN_MIN_DECIBELS = -100;
constexpr int MIN_DECIBELS_STEP = 1;

}

namespace Lab {

Figure3DWindow::Figure3DWindow(QWidget *parent)
		: QWidget(parent)
{
	ui_.setupUi(this);

	ui_.visualizationComboBox->blockSignals(true); // to avoid calling on_*ComboBox_currentIndexChanged
	for (const char** name = Visualization::valueNameList; *name; ++name) {
		ui_.visualizationComboBox->addItem(tr(*name));
	}
	ui_.visualizationComboBox->blockSignals(false);
	ui_.visualizationComboBox->setCurrentIndex(static_cast<int>(ui_.oglFigureWidget->visualization()));

	ui_.minDecibelsComboBox->blockSignals(true); // to avoid calling on_*ComboBox_currentIndexChanged
	for (int i = MAX_MIN_DECIBELS; i >= MIN_MIN_DECIBELS; i -= MIN_DECIBELS_STEP) {
		ui_.minDecibelsComboBox->addItem(QString::number(i));
	}
	ui_.minDecibelsComboBox->blockSignals(false);
	const int comboIdx = ui_.minDecibelsComboBox->findText(
				QString::number(static_cast<int>(ui_.oglFigureWidget->minDecibels())));
	if (comboIdx != -1) {
		ui_.minDecibelsComboBox->setCurrentIndex(comboIdx);
	}

	ui_.colormapComboBox->blockSignals(true); // to avoid calling on_*ComboBox_currentIndexChanged
	for (const char** name = Colormap::gradientNameList; *name; ++name) {
		ui_.colormapComboBox->addItem(tr(*name));
	}
	ui_.colormapComboBox->blockSignals(false);
	ui_.colormapComboBox->setCurrentIndex(static_cast<int>(ui_.oglFigureWidget->colormap()));
}

Figure3DWindow::~Figure3DWindow()
{
}

void
Figure3DWindow::updateData(float valueScale,
				Matrix<XYZValue<float>>* gridData,
				const std::vector<XYZ<float>>* pointList)
{
	ui_.oglFigureWidget->updateData(valueScale, gridData, pointList);
}

void
Figure3DWindow::on_visualizationComboBox_currentIndexChanged(int index)
{
	ui_.oglFigureWidget->setVisualization(static_cast<Visualization::Value>(index));
}

void Figure3DWindow::on_visualizationComboBox_activated(int /*index*/)
{
	ui_.oglFigureWidget->setUseManualSettings(true);
}

void
Figure3DWindow::on_showInfoCheckBox_stateChanged(int state)
{
	ui_.oglFigureWidget->setShowInfo(state == Qt::Checked);
}

void
Figure3DWindow::on_rotationCheckBox_stateChanged(int state)
{
	ui_.oglFigureWidget->setRotationMode(state == Qt::Checked);
}

void
Figure3DWindow::on_minDecibelsComboBox_currentIndexChanged(const QString& text)
{
	ui_.oglFigureWidget->setMinDecibels(static_cast<float>(text.toInt()));
}

void
Figure3DWindow::on_colormapComboBox_currentIndexChanged(int index)
{
	ui_.oglFigureWidget->setColormap(static_cast<Colormap::Gradient>(index));
}

void
Figure3DWindow::on_colormapComboBox_activated(int /*index*/)
{
	ui_.oglFigureWidget->setUseManualSettings(true);
}

void
Figure3DWindow::setVisualization(Visualization::Value visualization)
{
	if (ui_.oglFigureWidget->useManualSettings()) return;

	if (visualization != Visualization::VALUE_DEFAULT) {
		ui_.oglFigureWidget->setVisualization(visualization);
		ui_.visualizationComboBox->setCurrentIndex(static_cast<int>(visualization));
	}
}

void
Figure3DWindow::setColormap(Colormap::Gradient colormap)
{
	if (ui_.oglFigureWidget->useManualSettings()) return;

	if (colormap != Colormap::GRADIENT_DEFAULT) {
		ui_.oglFigureWidget->setColormap(colormap);
		ui_.colormapComboBox->setCurrentIndex(static_cast<int>(colormap));
	}
}

void
Figure3DWindow::closeEvent(QCloseEvent* event)
{
	ui_.oglFigureWidget->setUseManualSettings(false);
	event->accept();
}

} // namespace Lab
