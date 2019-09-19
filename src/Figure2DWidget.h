/***************************************************************************
 *  Copyright 2019 Marcelo Y. Matuda                                       *
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

#ifndef FIGURE_2D_WIDGET_H
#define FIGURE_2D_WIDGET_H

#include <algorithm>
#include <iostream>
#include <vector>

#include <QString>
#include <QWidget>

namespace Lab {

class Figure2DWidget : public QWidget {
	Q_OBJECT
public:
	explicit Figure2DWidget(QWidget* parent=0);

	template<typename T> void updateData(const std::vector<T>& x, const std::vector<T>& y);
	void setXLabel(QString label);
	void setYLabel(QString label);
	void enablePointMarker(bool enable) { drawPointMarker_ = enable; }
protected:
	virtual void paintEvent(QPaintEvent* event);
	virtual void mouseDoubleClickEvent(QMouseEvent* event);
	virtual void mousePressEvent(QMouseEvent* event);
	virtual void mouseMoveEvent(QMouseEvent* event);
	virtual void resizeEvent(QResizeEvent* event);
	virtual void wheelEvent(QWheelEvent* event);
	virtual void keyPressEvent(QKeyEvent* event);
private:
	void handleTransform();
	void autoSetAxesTicks(bool expand=false);
	bool resetFigure();

	static void autoSetAxisTicks(double minValue, double maxValue,
					std::vector<double>& ticks, double& offset, double& coef,
					bool expand);

	bool figureChanged_;
	bool drawPointMarker_;
	int leftMargin_;
	int rightMargin_;
	int topMargin_;
	int bottomMargin_;
	int xLabelWidth_;
	int yLabelWidth_;
	int textCapHeight_;
	int mainAreaWidth_;
	int mainAreaHeight_;
	double xScale_;
	double yScale_;
	double xBegin_;
	double xEnd_;
	double yBegin_;
	double yEnd_;
	double xTickOffset_;
	double yTickOffset_;
	double xTickCoef_;
	double yTickCoef_;
	double lastXBegin_;
	double lastXEnd_;
	double lastXScale_;
	std::vector<float> xList_;
	std::vector<float> yList_;
	std::vector<double> xTicks_;
	std::vector<double> yTicks_;
	std::vector<int> yTicksWidth_;
	QString xLabel_;
	QString yLabel_;
	QPointF lastMousePos_;
};

template<typename T>
void
Figure2DWidget::updateData(const std::vector<T>& x, const std::vector<T>& y)
{
	if (x.size() != y.size()) {
		std::cerr << "[Figure2DWidget::updateData] Arrays x and y with different sizes." << std::endl;
		return;
	}
	if (x.size() < 2U) {
		std::cerr << "[Figure2DWidget::updateData] Arrays x and y are too small." << std::endl;
		return;
	}

	xList_.resize(x.size());
	yList_.resize(y.size());
	for (unsigned int i = 0, end = x.size(); i < end; ++i) {
		xList_[i] = x[i];
		yList_[i] = y[i];
	}

	if (!resetFigure()) {
		xList_.clear();
		yList_.clear();
	}
}

} // namespace Lab

#endif // FIGURE_2D_WIDGET_H
