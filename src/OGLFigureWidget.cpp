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

//
// Coordinates:
//
// Grid data:
//   rows: angles, x
//   cols: radiuses, z
//  x ^               x ^
//    |                 |
// |-||     z        |-||
// | |o----->   or   | ||     z
// |-|               |-|o----->
// array
//
// Model:
//          y ^
//            |
//            |     x
//          z o----->
//
// Widget:
//         x
//   o----->-----------|
//   |                 |
//   |                 |
// y V                 |
//   |                 |
//   |-----------------|
//
//==============================================================================
// Colormaps Viridis, Plasma, Inferno and Magma obtained from Matplotlib 2.02.
//

#include "OGLFigureWidget.h"

#include <algorithm> /* std::max, std::swap */
#include <cstddef> /* std::size_t */
#include <cmath>

#include <GL/glu.h>

#include <QColor>
#include <QKeyEvent>
#include <QMouseEvent>
#include <QPainter>
#include <QWheelEvent>

#include "Geometry.h"
#include "Log.h"
#include "ParallelHilbertEnvelope.h"
#include "Util.h"



namespace {

constexpr float DEFAULT_MIN_DECIBELS = -60.0f;
constexpr int MININUM_WIDTH = 400;
constexpr int MININUM_HEIGHT = 300;
constexpr float MIN_VALUE_SCALE = 0.1f;
constexpr float MAX_VALUE_SCALE = 10.0f;
constexpr float MIN_SCALE = 0.1f;
constexpr float MAX_SCALE = 100.0f;
constexpr float MARGIN = 0.005f;
constexpr float MAX_OFFSET = 1.0e5f;
constexpr float ZOOM_IN_FACTOR = 1.1f;
constexpr float ZOOM_OUT_FACTOR = 0.9f;
//constexpr float REFLECTOR_POINT_SIZE_SMALL = 1.0f;
constexpr float REFLECTOR_POINT_SIZE_BIG = 3.0f;
//constexpr unsigned int REFLECTOR_SIZE_THRESHOLD = 20;
constexpr float VALUE_EPS = 1e-30f;

} // namespace

namespace Lab {

OGLFigureWidget::OGLFigureWidget(QWidget* parent)
		: QOpenGLWidget(parent)
		, valuesInDecibel_()
		, editingDistanceMarkers_()
		, rotationMode_()
		, showPoints_(true)
		, showInfo_(true)
		, dataChanged_()
		, useManualSettings_()
		, colormap_(Colormap::GRADIENT_GRAY)
		, visualization_(Visualization::VALUE_RAW_LINEAR)
		, minDecibels_(DEFAULT_MIN_DECIBELS)
		, minValue_(Util::decibelsToLinear(minDecibels_))
		, scale_(1.0)
		, rotX_()
		, rotZ_()
		, minX_()
		, maxX_()
		, minY_()
		, maxY_()
		, minZ_()
		, maxZ_()
		, dataValueScale_()
		, valueScale_(1.0)
		, maxAbsLevel_()
		, maxAbsLevelDecibels_()
		, offsetX_()
		, offsetY_()
		, normal_{0.0, 0.0, 0.0}
#ifndef OGLFIGUREWIDGET_USE_VERTEX_ARRAY
		, oglDisplayList_()
#endif
{
	setMinimumWidth(MININUM_WIDTH);
	setMinimumHeight(MININUM_HEIGHT);
}

float
OGLFigureWidget::calcValueFactor(const Matrix<XYZValue<float>>& data)
{
	maxAbsLevel_ = Util::maxAbsoluteValueField<XYZValue<float>, float>(data);
	const float factor = (dataValueScale_ != 0.0) ? dataValueScale_ : 1.0f / std::max(maxAbsLevel_, VALUE_EPS);
	maxAbsLevel_ *= factor;
	maxAbsLevelDecibels_ = Util::linearToDecibels(maxAbsLevel_);
	return factor;
}

#ifdef OGLFIGUREWIDGET_USE_VERTEX_ARRAY
void
OGLFigureWidget::fillIndexArray(unsigned int iA, unsigned int iB, std::vector<GLuint>& indexArray)
{
	indexArray.clear();
	indexArray.push_back(iA++);
	indexArray.push_back(iB++);
	for (unsigned int j = 0, end = oglGridData_.n2() - 1; j < end; ++j, ++iA, ++iB) {
		if (j & 1U) {
			indexArray.push_back(iB);
			indexArray.push_back(iA);
			indexArray.push_back(iB); // create a degenerate triangle
		} else {
			indexArray.push_back(iA);
			indexArray.push_back(iB);
			indexArray.push_back(iA); // create a degenerate triangle
		}
	}
}
#endif

template<typename ColorScale>
void
OGLFigureWidget::fillOGLGridDataWithAbsValues(const Matrix<XYZValue<float>>& data, float valueFactor)
{
	auto srcIter = data.begin();
	auto srcEnd = data.end();
	auto destIter = oglGridData_.begin();
	for ( ; srcIter != srcEnd; ++srcIter, ++destIter) {
		const float value = std::abs(srcIter->value) * valueFactor;
		destIter->pos.x = srcIter->x + normal_.x * value;
		destIter->pos.y = srcIter->y + normal_.y * value;
		destIter->pos.z = srcIter->z + normal_.z * value;
		ColorScale::setColor(value, *destIter);
	}
}

template<typename ColorScale>
void
OGLFigureWidget::fillOGLGridDataWithLogAbsValues(const Matrix<XYZValue<float>>& data, float valueFactor)
{
	const float invMinusMinDecibels = -1.0f / minDecibels_;
	auto srcIter = data.begin();
	auto srcEnd = data.end();
	auto destIter = oglGridData_.begin();
	for ( ; srcIter != srcEnd; ++srcIter, ++destIter) {
		float value = std::abs(srcIter->value) * valueFactor;
		value = (value <= minValue_) ? 0.0f : (Util::linearToDecibels(value) - minDecibels_) * invMinusMinDecibels;
		destIter->pos.x = srcIter->x + normal_.x * value;
		destIter->pos.y = srcIter->y + normal_.y * value;
		destIter->pos.z = srcIter->z + normal_.z * value;
		ColorScale::setColor(value, *destIter);
	}
}

template<typename ColorScale>
void
OGLFigureWidget::fillOGLGridData()
{
	if (oglGridData_.n1() != gridData_.n1() || oglGridData_.n2() != gridData_.n2()) {
		oglGridData_.resize(gridData_.n1(), gridData_.n2());
	}

	switch (visualization_) {
	case Visualization::VALUE_DEFAULT: // falls through
	case Visualization::VALUE_RAW_LINEAR:
		{
			const float valueFactor = 0.5f * calcValueFactor(gridData_);
			auto srcIter = gridData_.begin();
			auto srcEnd = gridData_.end();
			auto destIter = oglGridData_.begin();
			for ( ; srcIter != srcEnd; ++srcIter, ++destIter) {
				const float value = srcIter->value * valueFactor + 0.5f;
				destIter->pos.x = srcIter->x + normal_.x * value;
				destIter->pos.y = srcIter->y + normal_.y * value;
				destIter->pos.z = srcIter->z + normal_.z * value;
				ColorScale::setColor(value, *destIter);
			}
		}
		break;
	case Visualization::VALUE_RECTIFIED_LINEAR:
		{
			const float valueFactor = calcValueFactor(gridData_);
			fillOGLGridDataWithAbsValues<ColorScale>(gridData_, valueFactor);
		}
		break;
	case Visualization::VALUE_RECTIFIED_LOG:
		{
			const float valueFactor = calcValueFactor(gridData_);
			fillOGLGridDataWithLogAbsValues<ColorScale>(gridData_, valueFactor);
		}
		break;
	case Visualization::VALUE_ENVELOPE_LINEAR:
		{
			Matrix<XYZValue<float>> envelope = gridData_;
			ParallelHilbertEnvelope<float>::calculateDim2Value(envelope);
			const float valueFactor = calcValueFactor(envelope);
			fillOGLGridDataWithAbsValues<ColorScale>(envelope, valueFactor);
		}
		break;
	case Visualization::VALUE_ENVELOPE_LOG:
		{
			Matrix<XYZValue<float>> envelope = gridData_;
			ParallelHilbertEnvelope<float>::calculateDim2Value(envelope);
			const float valueFactor = calcValueFactor(envelope);
			fillOGLGridDataWithLogAbsValues<ColorScale>(envelope, valueFactor);
		}
		break;
	}
}

void
OGLFigureWidget::updateGridData(float dataValueScale, Matrix<XYZValue<float>>& gridData)
{
	if (gridData.n1() < 2U || gridData.n2() < 2U) return;

	const auto& p0 = gridData(0, 0);
	const auto& p1 = gridData(gridData.n1() - 1U, 0);
	const auto& p2 = gridData(0, gridData.n2() - 1U);
	const float dx1 = p1.x - p0.x;
	const float dy1 = p1.y - p0.y;
	const float dz1 = p1.z - p0.z;
	const float dx2 = p2.x - p0.x;
	const float dy2 = p2.y - p0.y;
	const float dz2 = p2.z - p0.z;
	if (dx1 == 0.0 && dx2 == 0.0) {
		normal_ = {-1.0, 0.0, 0.0};
	} else if (dy1 == 0.0 && dy2 == 0.0) {
		normal_ = {0.0, -1.0, 0.0};
	} else if (dz1 == 0.0 && dz2 == 0.0) {
		normal_ = {0.0, 0.0, -1.0};
	} else {
		THROW_EXCEPTION(InvalidParameterException, "Invalid grid data.");
	}

	dataValueScale_ = dataValueScale;
	std::swap(gridData, gridData_);

	Matrix<XYZValue<float>>::ConstIterator iter = gridData_.begin(), end = gridData_.end();
	minX_ = maxX_ = iter->x;
	minY_ = maxY_ = iter->y;
	minZ_ = maxZ_ = iter->z;
	for (++iter; iter != end; ++iter) {
		const float x = iter->x;
		const float y = iter->y;
		const float z = iter->z;
		if (x < minX_) {
			minX_ = x;
		} else if (x > maxX_) {
			maxX_ = x;
		}
		if (y < minY_) {
			minY_ = y;
		} else if (y > maxY_) {
			maxY_ = y;
		}
		if (z < minZ_) {
			minZ_ = z;
		} else if (z > maxZ_) {
			maxZ_ = z;
		}
	}

	updateDataVisualization();

	update();

	dataChanged_ = true;
}

void
OGLFigureWidget::updatePointList(const std::vector<XYZ<float>>& pointList)
{
	if (pointList.empty()) return;

	pointList_.resize(pointList.size());

	std::vector<XYZ<float>>::const_iterator srcIter = pointList.begin(), srcIterEnd = pointList.end();
	std::vector<XYZ<float>>::iterator destIter = pointList_.begin();
	for ( ; srcIter != srcIterEnd; ++srcIter, ++destIter) {
		const float x = srcIter->x;
		const float y = srcIter->y;
		const float z = srcIter->z;
		destIter->x = x;
		destIter->y = y;
		destIter->z = z;
	}

	update();
}

void
OGLFigureWidget::setMinDecibels(float minDecibels)
{
	minDecibels_ = minDecibels;
	minValue_ = Util::decibelsToLinear(minDecibels_);

	if (visualization_ == Visualization::VALUE_RECTIFIED_LOG ||
			visualization_ == Visualization::VALUE_ENVELOPE_LOG) {
		updateDataVisualization();
		update();
	}
}

void
OGLFigureWidget::updateData(float dataValueScale,
				Matrix<XYZValue<float>>* gridData,
				const std::vector<XYZ<float>>* pointList)
{
	if (!gridData) return;

	updateGridData(dataValueScale, *gridData);
	if (pointList) updatePointList(*pointList);
}

void
OGLFigureWidget::updateDataVisualization()
{
	if (gridData_.empty()) return;

	switch (colormap_) {
	case Colormap::DEFAULT: // falls through
#define COLORMAP_ITEM(A, B, C) case Colormap::A: fillOGLGridData<C>(); break;
	COLORMAP_TABLE
#undef COLORMAP_ITEM
	}
}

void
OGLFigureWidget::setVisualization(Visualization::Value visualization)
{
	visualization_ = visualization;
	updateDataVisualization();
	update();
}

void
OGLFigureWidget::setColormap(Colormap::Id colormap)
{
	colormap_ = colormap;
	updateDataVisualization();
	update();
}

void
OGLFigureWidget::paintGL()
{
	QPainter painter(this);
	painter.setFont(QFont("monospace"));

	//==================================================
	painter.beginNativePainting();

#ifdef OGLFIGUREWIDGET_USE_VERTEX_ARRAY
	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_COLOR_ARRAY);
#endif

	// QPainter disables GL_DEPTH_TEST.
	glEnable(GL_DEPTH_TEST);

	// QPainter clears the depth buffer.
	//glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClear(GL_COLOR_BUFFER_BIT);

	const float dx = maxX_ - minX_;
	const float dy = maxY_ - minY_;
	const float dz = maxZ_ - minZ_;
	const float maxDelta = std::max(std::max(dx, dy), dz);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	const float aspect = static_cast<float>(width()) / height();
	const float yL = maxDelta * (0.5 + MARGIN);
	const float xL = yL * aspect;
	const float totalValueScale = 0.5f * yL * valueScale_;
	const float zL = scale_ * std::sqrt(xL * xL + yL * yL + totalValueScale * totalValueScale);
	glOrtho(-xL, xL, -yL, yL, -zL, zL); // l, r, b, t, near, far

	const float resX = (2.0 * xL) / width();
	const float resY = (2.0 * yL) / height();

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glTranslatef(resX * offsetX_, resY * offsetY_, 0.0f);
	glRotatef(rotX_, 1.0f, 0.0f, 0.0f);
	glRotatef(rotZ_, 0.0f, 0.0f, 1.0f);
	glScalef(scale_, scale_, scale_);
	glRotatef(-90.0f, 1.0f, 0.0f, 0.0f);

	glPushMatrix();

	glScalef(
		normal_.x != 0.0f ? totalValueScale : 1.0f,
		normal_.y != 0.0f ? totalValueScale : 1.0f,
		normal_.z != 0.0f ? totalValueScale : 1.0f);

	// Draw frame.
	glColor3f(0.5f, 0.5f, 0.5f);
	if (normal_.x != 0.0f) {
		glTranslatef(0.0f, -(minY_ + maxY_) * 0.5f, -(minZ_ + maxZ_) * 0.5f);
		glBegin(GL_LINE_LOOP);
			glVertex3f(normal_.x, minY_, minZ_);
			glVertex3f(normal_.x, minY_, maxZ_);
			glVertex3f(normal_.x, maxY_, maxZ_);
			glVertex3f(normal_.x, maxY_, minZ_);
		glEnd();
		glBegin(GL_LINE_LOOP);
			glVertex3f(0.0f, minY_, minZ_);
			glVertex3f(0.0f, minY_, maxZ_);
			glVertex3f(0.0f, maxY_, maxZ_);
			glVertex3f(0.0f, maxY_, minZ_);
		glEnd();
		glTranslatef(-(minX_ + maxX_) * 0.5f, 0.0f, 0.0f);
	} else if (normal_.y != 0.0f) {
		glTranslatef(-(minX_ + maxX_) * 0.5f, 0.0f, -(minZ_ + maxZ_) * 0.5f);
		glBegin(GL_LINE_LOOP);
			glVertex3f(minX_, normal_.y, minZ_);
			glVertex3f(minX_, normal_.y, maxZ_);
			glVertex3f(maxX_, normal_.y, maxZ_);
			glVertex3f(maxX_, normal_.y, minZ_);
		glEnd();
		glBegin(GL_LINE_LOOP);
			glVertex3f(minX_, 0.0f, minZ_);
			glVertex3f(minX_, 0.0f, maxZ_);
			glVertex3f(maxX_, 0.0f, maxZ_);
			glVertex3f(maxX_, 0.0f, minZ_);
		glEnd();
		glTranslatef(0.0f, -(minY_ + maxY_) * 0.5f, 0.0f);
	} else if (normal_.z != 0.0f) {
		glTranslatef(-(minX_ + maxX_) * 0.5f, -(minY_ + maxY_) * 0.5f, 0.0f);
		glBegin(GL_LINE_LOOP);
			glVertex3f(minX_, minY_, normal_.z);
			glVertex3f(minX_, maxY_, normal_.z);
			glVertex3f(maxX_, maxY_, normal_.z);
			glVertex3f(maxX_, minY_, normal_.z);
		glEnd();
		glBegin(GL_LINE_LOOP);
			glVertex3f(minX_, minY_, 0.0f);
			glVertex3f(minX_, maxY_, 0.0f);
			glVertex3f(maxX_, maxY_, 0.0f);
			glVertex3f(maxX_, minY_, 0.0f);
		glEnd();
		glTranslatef(0.0f, 0.0f, -(minZ_ + maxZ_) * 0.5f);
	}

	if (oglGridData_.n1() >= 2) {
#ifdef OGLFIGUREWIDGET_USE_VERTEX_ARRAY
		if (dataChanged_) {
			fillIndexArray(0, oglGridData_.n2(), evenIndexArray_);
			fillIndexArray(oglGridData_.n2(), 0, oddIndexArray_);
		}

		for (unsigned int i = 0; i < oglGridData_.n1() - 1; ++i) {
			glVertexPointer(3, GL_FLOAT, sizeof(OGLPoint3D), &oglGridData_(i, 0).pos);
			glColorPointer(3, GL_FLOAT, sizeof(OGLPoint3D), &oglGridData_(i, 0).color);
			if (i & 1U) {
				glDrawElements(GL_TRIANGLE_STRIP, oddIndexArray_.size(), GL_UNSIGNED_INT, oddIndexArray_.data());
			} else {
				glDrawElements(GL_TRIANGLE_STRIP, evenIndexArray_.size(), GL_UNSIGNED_INT, evenIndexArray_.data());
			}
		}
#else
		if (dataChanged_) {
			glDeleteLists(oglDisplayList_, 1);
			oglDisplayList_ = glGenLists(1);
			glNewList(oglDisplayList_, GL_COMPILE);

			for (unsigned int i = 0; i < oglGridData_.n1() - 1; ++i) {
				const bool odd = i & 1U;
				auto rangeA = oglGridData_.range2(odd ? i + 1 : i    );
				auto rangeB = oglGridData_.range2(odd ? i     : i + 1);

				glBegin(GL_TRIANGLE_STRIP);
				auto iterA = rangeA.begin();
				auto iterB = rangeB.begin();
				if (iterA != rangeA.end()) {
					const OGLPoint3D pA = *iterA;
					createVertex(pA);

					const OGLPoint3D pB = *iterB;
					createVertex(pB);

					++iterA; ++iterB;
					for (unsigned int j = 0; iterA != rangeA.end(); ++iterA, ++iterB, ++j) {
						const OGLPoint3D pA = *iterA;
						const OGLPoint3D pB = *iterB;
						if (j & 1U) {
							createVertex(pB);
							createVertex(pA);
							createVertex(pB); // creates a degenerate triangle
						} else {
							createVertex(pA);
							createVertex(pB);
							createVertex(pA); // creates a degenerate triangle
						}
					}
				}
				glEnd();
			}

			glEndList();
		}

		glCallList(oglDisplayList_);
#endif
	}

	glPopMatrix();

	glTranslatef(-(minX_ + maxX_) * 0.5f, -(minY_ + maxY_) * 0.5f, -(minZ_ + maxZ_) * 0.5f);

	glDisable(GL_DEPTH_TEST);

	if (showPoints_ && !pointList_.empty()) {
//		if (pointList_.size() > REFLECTOR_SIZE_THRESHOLD) {
//			glPointSize(REFLECTOR_POINT_SIZE_SMALL);
//		} else {
			glPointSize(REFLECTOR_POINT_SIZE_BIG);
//		}
		glBegin(GL_POINTS);
		glColor3f(0.0f, 0.5f, 0.0f);
		for (std::vector<XYZ<float>>::const_iterator iter = pointList_.begin(); iter != pointList_.end(); ++iter) {
			glVertex3f(iter->x, iter->y, iter->z);
		}
		glEnd();
	}

	const bool mustShowDistance = !distanceMarker1_.isNull() && !distanceMarker2_.isNull();
	float distance = 0.0;
	if (mustShowDistance && normal_.y != 0.0) {
		// Calculate the distance between distanceMarker1_ and distanceMarker2_.

		GLint viewport[4];
		GLdouble modelViewMatrix[16], projectionMatrix[16];
		glGetIntegerv(GL_VIEWPORT, viewport);
		glGetDoublev(GL_MODELVIEW_MATRIX, modelViewMatrix);
		glGetDoublev(GL_PROJECTION_MATRIX, projectionMatrix);

		const GLint maxGlY = viewport[3] - 1; // viewport[3] is the height of the window in pixels
		GLdouble p1X, p1Y, p1Z, p2X, p2Y, p2Z;
		gluUnProject(static_cast<GLfloat>(distanceMarker1_.x()), static_cast<GLfloat>(maxGlY - distanceMarker1_.y()), 0.0,
				modelViewMatrix, projectionMatrix, viewport, &p1X, &p1Y, &p1Z);
		gluUnProject(static_cast<GLfloat>(distanceMarker2_.x()), static_cast<GLfloat>(maxGlY - distanceMarker2_.y()), 0.0,
				modelViewMatrix, projectionMatrix, viewport, &p2X, &p2Y, &p2Z);
		distance = Geometry::distance2D(p1X, p1Z, p2X, p2Z);
	}

#ifdef OGLFIGUREWIDGET_USE_VERTEX_ARRAY
	glDisableClientState(GL_COLOR_ARRAY);
	glDisableClientState(GL_VERTEX_ARRAY);
#endif

	painter.endNativePainting();
	//==================================================

	painter.setRenderHint(QPainter::TextAntialiasing);
	painter.setPen(Qt::white);
	painter.setBackground(QBrush{Qt::black});
	painter.setBackgroundMode(Qt::OpaqueMode);

	if (mustShowDistance && normal_.y != 0.0) {
		painter.setPen(Qt::darkGray);
		painter.drawLine(distanceMarker1_, distanceMarker2_);
		painter.setPen(Qt::white);
		painter.drawText(10, 20, QString("Distance = ") + QString::number(distance));
	}
	if (showInfo_) {
		painter.drawText(10,  40, QString("off x ") + QString::number(offsetX_));
		painter.drawText(10,  60, QString("off y ") + QString::number(offsetY_));
		painter.drawText(10,  80, QString("rot x ") + QString::number(rotX_));
		painter.drawText(10, 100, QString("rot z ") + QString::number(rotZ_));
		painter.drawText(10, 120, QString("zoom  ") + QString::number(scale_));
		painter.drawText(10, 140, QString("scale ") + QString::number(valueScale_));
		painter.drawText(10, 160, QString("level ") + QString::number(maxAbsLevel_));
		painter.drawText(10, 180, QString("   dB ") + (std::abs(maxAbsLevelDecibels_) < 0.1f ?
								"0.0" :
								QString::number(maxAbsLevelDecibels_, 'f', 1)));
		painter.drawText(10, 200, QString("size  %1 x %2").arg(oglGridData_.n1()).arg(oglGridData_.n2()));
	}

	dataChanged_ = false;
}

void
OGLFigureWidget::mousePressEvent(QMouseEvent* e)
{
	mouseAnchor_ = e->pos();

	if (editingDistanceMarkers_) {
		distanceMarker1_ = QPoint(mouseAnchor_.x(), mouseAnchor_.y());
	}
}

void
OGLFigureWidget::mouseMoveEvent(QMouseEvent* e)
{
	QPoint diff = e->pos() - mouseAnchor_;
	if (e->buttons() & Qt::LeftButton) {
		mouseAnchor_ = e->pos();
		if (rotationMode_) {
			rotX_ += diff.y() * (180.0f / height());
			Util::clipAngleInDegrees180(rotX_);
			rotZ_ += diff.x() * (180.0f / width());
			Util::clipAngleInDegrees180(rotZ_);
		} else {
			if (editingDistanceMarkers_) {
				distanceMarker2_ = QPoint(mouseAnchor_.x(), mouseAnchor_.y());
			} else {
				offsetX_ += diff.x();
				Util::clip(offsetX_, -MAX_OFFSET, MAX_OFFSET);
				offsetY_ -= diff.y();
				Util::clip(offsetY_, -MAX_OFFSET, MAX_OFFSET);
			}
		}
	} else if (e->buttons() & Qt::RightButton) {
		if (rotationMode_) {
			valueScale_ -= diff.y() * ((MAX_VALUE_SCALE - MIN_VALUE_SCALE) / height());
			Util::clip(valueScale_, MIN_VALUE_SCALE, MAX_VALUE_SCALE);
			mouseAnchor_ = e->pos();
		} else {
			return;
		}
	}

	update();
}

void
OGLFigureWidget::wheelEvent(QWheelEvent* e)
{
	editingDistanceMarkers_ = false;
	distanceMarker1_ = QPoint(); // clear
	distanceMarker2_ = QPoint(); // clear

	const float oldScale = scale_;
	if (e->delta() > 0) {
		scale_ *= ZOOM_IN_FACTOR;
	} else {
		scale_ *= ZOOM_OUT_FACTOR;
	}
	Util::clip(scale_, MIN_SCALE, MAX_SCALE);

	const float factor = scale_ / oldScale;
	offsetX_ *= factor;
	Util::clip(offsetX_, -MAX_OFFSET, MAX_OFFSET);
	offsetY_ *= factor;
	Util::clip(offsetY_, -MAX_OFFSET, MAX_OFFSET);

	update();
}

void
OGLFigureWidget::keyPressEvent(QKeyEvent* e)
{
	switch (e->key()) {
	case Qt::Key_Space:
		scale_ = 1.0;
		rotX_ = 0.0;
		rotZ_ = 0.0;
		offsetX_ = 0.0;
		offsetY_ = 0.0;
		valueScale_ = 1.0;

		editingDistanceMarkers_ = false;
		distanceMarker1_ = QPoint(); // clear
		distanceMarker2_ = QPoint(); // clear

		break;
	case Qt::Key_Control:
		if (rotationMode_) return;
		rotX_ = 0.0;
		editingDistanceMarkers_ = true;
		distanceMarker1_ = QPoint(); // clear
		distanceMarker2_ = QPoint(); // clear
		break;
	case Qt::Key_P:
		showPoints_ = !showPoints_;
		break;
	case Qt::Key_X:
		rotX_ = -90.0;
		rotZ_ = 90.0;
		break;
	case Qt::Key_Z:
		rotX_ = -90.0;
		rotZ_ = 0.0;
		break;
	default:
		return;
	}

	update();
}

void
OGLFigureWidget::keyReleaseEvent(QKeyEvent* e)
{
	if (e->key() == Qt::Key_Control) {
		editingDistanceMarkers_ = false;
	}
}

void
OGLFigureWidget::initializeGL()
{
	//LOG_DEBUG << "initializeGL";

	//glShadeModel(GL_FLAT);

	//glClearColor(0.0, 0.0, 0.0, 1.0);
	glClearColor(1.0, 1.0, 1.0, 1.0);
}

void
OGLFigureWidget::resizeGL(int /*width*/, int /*height*/)
{
	// paintGL calls glViewport().
	//glViewport(0, 0, width, height);
}

} // namespace Lab
