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

#include "Figure2DWidget.h"

#include <cmath>

#include <QMouseEvent>
#include <QPainter>
#include <QPointF>
#include <QRectF>

#define SPACING 10
#define TEXT_SPACING 10
#define TICK_SIZE 5
#define MIN_AXIS_DIV 8
#define WHEEL_ZOOM_FACTOR (0.01)
#define MOUSE_ZOOM_FACTOR (0.002)
#define DEFAULT_Y_DELTA (1.0)
#define MIN_VALUE_DELTA (1.0e-30)
#define MAX_VALUE_DELTA (1.0e30)
#define MARKER_SIZE (2.0)
#define EPS (1.0e-5)
#define MAX_TICK_POW (2.0)
#define MIN_TICK_POW (-2.0)

namespace Lab {

Figure2DWidget::Figure2DWidget(QWidget* parent)
		: QWidget(parent)
		, figureChanged_()
		, drawPointMarker_()
		, leftMargin_()
		, rightMargin_()
		, topMargin_()
		, bottomMargin_()
		, xLabelWidth_()
		, yLabelWidth_()
		, textCapHeight_()
		, mainAreaWidth_()
		, mainAreaHeight_()
		, xScale_()
		, yScale_()
		, xBegin_()
		, xEnd_()
		, yBegin_()
		, yEnd_()
		, xTickOffset_()
		, yTickOffset_()
		, xTickCoef_(1.0)
		, yTickCoef_(1.0)
		, xLabel_("x")
		, yLabel_("y")
{
	setAutoFillBackground(true);
	setBackgroundRole(QPalette::Base);

	setFocusPolicy(Qt::StrongFocus); // enable key events
}

// Note: with no antialiasing, the coordinates in QPointF are rounded to the nearest integer.
void
Figure2DWidget::paintEvent(QPaintEvent* /*event*/)
{
	if (xList_.empty() || xTicks_.empty()) return;

	QPainter painter(this);
	painter.setFont(QFont("monospace"));

	if (figureChanged_) {
		QFontMetrics fm = painter.fontMetrics();
		textCapHeight_ = fm.capHeight();

		yTicksWidth_.resize(yTicks_.size());
		int maxW = 0;
		for (unsigned int i = 0; i < yTicks_.size(); ++i) {
			yTicksWidth_[i] = fm.width(QString::number(yTicks_[i]));
			if (yTicksWidth_[i] > maxW) maxW = yTicksWidth_[i];
		}
		leftMargin_ = SPACING + textCapHeight_ + 2 * TEXT_SPACING + maxW + TICK_SIZE;

		rightMargin_ = fm.width(QString::number(xTicks_.back())) + SPACING;
		topMargin_ = textCapHeight_ + SPACING;
		bottomMargin_ = TICK_SIZE + 2 * TEXT_SPACING + 2 * textCapHeight_ + SPACING;
		xLabelWidth_ = fm.width(xLabel_);
		yLabelWidth_ = fm.width(yLabel_);

		handleTransform();

		figureChanged_ = false;
	}

	const double uBegin = leftMargin_;
	const double vBegin = topMargin_ + mainAreaHeight_;
	const double uEnd = leftMargin_ + mainAreaWidth_;
	const double vEnd = topMargin_;

	// X axis.
	const double vTick = vBegin + TICK_SIZE;
	const double vText = vTick + TEXT_SPACING + textCapHeight_;
	for (float tick : xTicks_) {
		const double u = uBegin + (tick * xTickCoef_ + xTickOffset_ - xBegin_) * xScale_;
		// Vertical line.
		painter.setPen(Qt::lightGray);
		painter.drawLine(QPointF(u, vBegin), QPointF(u, vEnd));
		painter.setPen(Qt::black);
		// Tick.
		painter.drawLine(QPointF(u, vBegin), QPointF(u, vTick));
		// Tich value.
		painter.drawText(QPointF(u, vText), QString::number(tick));
	}
	// X label.
	const double vLabelPos = vText + TEXT_SPACING + textCapHeight_;
	painter.drawText(QPointF(uBegin + mainAreaWidth_ / 2 - xLabelWidth_ / 2, vLabelPos), xLabel_);
	// X ticks transformation.
	QString xTickTransf;
	if (xTickCoef_   != 1.0) xTickTransf += QString("×1e%1 ").arg(std::log10(xTickCoef_));
	if (xTickOffset_ >  0.0) xTickTransf += "+";
	if (xTickOffset_ != 0.0) xTickTransf += QString("%1").arg(xTickOffset_);
	if (!xTickTransf.isEmpty()) {
		painter.drawText(QPointF(uBegin, vLabelPos), xTickTransf);
	}

	// Y axis.
	const double uTick = uBegin - TICK_SIZE;
	for (unsigned int i = 0; i < yTicks_.size(); ++i) {
		const double uText = uTick - TEXT_SPACING - yTicksWidth_[i];
		const double v = vBegin + (yTicks_[i] * yTickCoef_ + yTickOffset_ - yBegin_) * yScale_;
		// Horizontal line.
		painter.setPen(Qt::lightGray);
		painter.drawLine(QPointF(uBegin, v), QPointF(uEnd, v));
		painter.setPen(Qt::black);
		// Tick.
		painter.drawLine(QPointF(uBegin, v), QPointF(uTick, v));
		// Tick value.
		painter.drawText(QPointF(uText, v), QString::number(yTicks_[i]));
	}
	// Y label.
	painter.rotate(-90.0);
	painter.drawText(QPointF(-vBegin + mainAreaHeight_ / 2 - yLabelWidth_ / 2, SPACING + textCapHeight_), yLabel_);
	// Y ticks transformation.
	QString yTickTransf;
	if (yTickCoef_   != 1.0) yTickTransf += QString("×1e%1 ").arg(std::log10(yTickCoef_));
	if (yTickOffset_ >  0.0) yTickTransf += "+";
	if (yTickOffset_ != 0.0) yTickTransf += QString("%1").arg(yTickOffset_);
	if (!yTickTransf.isEmpty()) {
		painter.drawText(QPointF(-vBegin, SPACING + textCapHeight_), yTickTransf);
	}
	painter.resetTransform();

	// Frame.
	painter.drawRect(QRectF(QPointF(uBegin, vBegin), QPointF(uEnd, vEnd)));

	QPen pen2;
	pen2.setWidth(2);

	painter.setRenderHint(QPainter::Antialiasing);
	// When using antialiasing, 0.5 is added to (x, y) to match the case without antialiasing.

	// Draw curve.
	painter.setClipRect(uBegin, vEnd, mainAreaWidth_, mainAreaHeight_);
	painter.setPen(pen2);
	QPointF prevPoint(
			0.5 + uBegin + (xList_[0] - xBegin_) * xScale_,
			0.5 + vBegin + (yList_[0] - yBegin_) * yScale_);
	if (drawPointMarker_) {
		painter.drawEllipse(prevPoint, MARKER_SIZE, MARKER_SIZE);
	}
	for (unsigned int i = 1, size = xList_.size(); i < size; ++i) {
		QPointF currPoint(
			0.5 + uBegin + (xList_[i] - xBegin_) * xScale_,
			0.5 + vBegin + (yList_[i] - yBegin_) * yScale_);
		painter.drawLine(prevPoint, currPoint);
		if (drawPointMarker_) {
			painter.drawEllipse(currPoint, MARKER_SIZE, MARKER_SIZE);
		}
		prevPoint = currPoint;
	}
	painter.setClipping(false);
}

void
Figure2DWidget::mouseDoubleClickEvent(QMouseEvent* event)
{
	if (xList_.empty() || xTicks_.empty()) return;
	if (event->button() != Qt::LeftButton) return;

	resetFigure();
}

void
Figure2DWidget::mousePressEvent(QMouseEvent* event)
{
	if (xList_.empty() || xTicks_.empty()) return;

	lastMousePos_ = event->localPos();
}

void
Figure2DWidget::mouseMoveEvent(QMouseEvent* event)
{
	// Note: For mouse move events, event->button() == Qt::NoButton.

	if (xList_.empty() || xTicks_.empty()) return;

	QPointF delta = event->localPos() - lastMousePos_;
	if (event->buttons() & Qt::LeftButton) {
		const double dx = delta.x() / xScale_;
		const double dy = delta.y() / yScale_;
		xBegin_ -= dx;
		xEnd_   -= dx;
		yBegin_ -= dy;
		yEnd_   -= dy;
	} else if (event->buttons() & Qt::RightButton) {
		const double valueDelta = xEnd_ - xBegin_;
		const double factor = 1.0 + MOUSE_ZOOM_FACTOR * std::abs(delta.x());
		const double mean = 0.5 * (xBegin_ + xEnd_);
		if (delta.x() > 0.0) {
			if (valueDelta < MIN_VALUE_DELTA) {
				qDebug("Zoom limit.");
				return;
			}
			xBegin_ = mean + (xBegin_ - mean) / factor;
			xEnd_   = mean + (xEnd_   - mean) / factor;
		} else if (delta.x() < 0.0) {
			if (valueDelta > MAX_VALUE_DELTA) {
				qDebug("Zoom limit.");
				return;
			}
			xBegin_ = mean + (xBegin_ - mean) * factor;
			xEnd_   = mean + (xEnd_   - mean) * factor;
		}
	}

	lastMousePos_ = event->localPos();
	autoSetAxesTicks();
}

void
Figure2DWidget::wheelEvent(QWheelEvent* event)
{
	if (xList_.empty() || xTicks_.empty()) return;

	const double valueDelta = yEnd_ - yBegin_;
	QPoint angle = event->angleDelta() / 8; // degrees

	const double factor = 1.0 + WHEEL_ZOOM_FACTOR * std::abs(angle.y());
	const double mean = 0.5 * (yBegin_ + yEnd_);
	if (angle.y() > 0.0) {
		if (valueDelta < MIN_VALUE_DELTA) {
			qDebug("Zoom limit.");
			return;
		}
		yBegin_ = mean + (yBegin_ - mean) / factor;
		yEnd_   = mean + (yEnd_   - mean) / factor;
	} else if (angle.y() < 0.0) {
		if (valueDelta > MAX_VALUE_DELTA) {
			qDebug("Zoom limit.");
			return;
		}
		yBegin_ = mean + (yBegin_ - mean) * factor;
		yEnd_   = mean + (yEnd_   - mean) * factor;
	}

	autoSetAxesTicks();
}

void
Figure2DWidget::resizeEvent(QResizeEvent* /*event*/)
{
	if (xList_.empty() || xTicks_.empty()) return;

	handleTransform();
}

void
Figure2DWidget::keyPressEvent(QKeyEvent* event)
{
	if (xList_.empty() || xTicks_.empty()) {
		QWidget::keyPressEvent(event);
		return;
	}

	if (event->key() == Qt::Key_Space) {
		drawPointMarker_ = !drawPointMarker_;
	} else {
		QWidget::keyPressEvent(event);
		return;
	}

	update();
}

void
Figure2DWidget::handleTransform()
{
	mainAreaWidth_  =  width() - leftMargin_ -  rightMargin_;
	mainAreaHeight_ = height() -  topMargin_ - bottomMargin_;

	xScale_ =   mainAreaWidth_ / (xEnd_ - xBegin_);
	yScale_ = -mainAreaHeight_ / (yEnd_ - yBegin_);
}

void
Figure2DWidget::autoSetAxisTicks(double minValue, double maxValue,
					std::vector<double>& ticks, double& offset, double& coef,
					bool expand)
{
	ticks.clear();

	// Calculate step.
	const double range = maxValue - minValue;
	const double maxStep = range / MIN_AXIS_DIV;
	const double truncMaxStep = std::pow(10.0, std::floor(std::log10(maxStep)));
	double step = truncMaxStep;
	double aux;
	if ((aux = truncMaxStep * 5.0) < maxStep) {
		step = aux;
	} else if ((aux = truncMaxStep * 2.0) < maxStep) {
		step = aux;
	}
	//qDebug("maxStep: %f truncMaxStep: %f step: %f", maxStep, truncMaxStep, step);

	// Calculate tick values.
	double firstTick, lastTick;
	if (expand) {
		firstTick = std::floor(minValue / step) * step;
		lastTick  =  std::ceil(maxValue / step) * step;
	} else {
		firstTick =  std::ceil(minValue / step) * step;
		lastTick  = std::floor(maxValue / step) * step;
	}
	ticks.push_back(firstTick);
	for (unsigned int i = 0; ; ++i) {
		const double tick = firstTick + i * step;
		if (tick > lastTick + step * EPS) break;
		ticks.push_back(tick);
	}

	// Calculate offset.
	const double delta = ticks.back() - ticks.front();
	const double mean = 0.5 * (ticks.front() + ticks.back());
	if (mean > delta) {
		offset = ticks.front();
		for (double& t : ticks) t -= offset;
	} else if (mean < -delta) {
		offset = ticks.back();
		for (double& t : ticks) t -= offset;
	} else {
		offset = 0.0;
	}

	// Calculate coefficient.
	const double stepPow = std::log10(step);
	if (stepPow < MIN_TICK_POW) {
		coef = std::pow(10.0, std::floor(stepPow) - MIN_TICK_POW);
		for (double& t : ticks) t /= coef;
	} else if (stepPow > MAX_TICK_POW) {
		coef = std::pow(10.0, std::ceil(stepPow) - MAX_TICK_POW);
		for (double& t : ticks) t /= coef;
	} else {
		coef = 1.0;
	}
}

void
Figure2DWidget::autoSetAxesTicks(bool expand)
{
	autoSetAxisTicks(xBegin_, xEnd_, xTicks_, xTickOffset_, xTickCoef_, expand);
	autoSetAxisTicks(yBegin_, yEnd_, yTicks_, yTickOffset_, yTickCoef_, expand);

	figureChanged_ = true;
	update();
}

bool
Figure2DWidget::resetFigure()
{
	// The values may not be sorted.
	xEnd_   = *std::max_element(xList_.begin(), xList_.end());
	xBegin_ = *std::min_element(xList_.begin(), xList_.end());
	yEnd_   = *std::max_element(yList_.begin(), yList_.end());
	yBegin_ = *std::min_element(yList_.begin(), yList_.end());

	const double xValueDelta = xEnd_ - xBegin_;
	if (xValueDelta < MIN_VALUE_DELTA) {
		std::cerr << "[Figure2DWidget::resetFigure] Invalid x range (< " << MIN_VALUE_DELTA << ")." << std::endl;
		return false;
	}
	const double yValueDelta = yEnd_ - yBegin_;
	if (yValueDelta < MIN_VALUE_DELTA) {
		yBegin_ -= DEFAULT_Y_DELTA;
		yEnd_   += DEFAULT_Y_DELTA;
	}

	autoSetAxesTicks(true);

	xBegin_ = std::min(xBegin_, xTicks_.front() * xTickCoef_ + xTickOffset_);
	xEnd_   = std::max(xEnd_  , xTicks_.back()  * xTickCoef_ + xTickOffset_);
	yBegin_ = std::min(yBegin_, yTicks_.front() * yTickCoef_ + yTickOffset_);
	yEnd_   = std::max(yEnd_  , yTicks_.back()  * yTickCoef_ + yTickOffset_);

	figureChanged_ = true;
	update();
	return true;
}

void
Figure2DWidget::setXLabel(QString label)
{
	xLabel_ = label;
	figureChanged_ = true;
	update();
}

void
Figure2DWidget::setYLabel(QString label)
{
	yLabel_ = label;
	figureChanged_ = true;
	update();
}

} // namespace Lab
