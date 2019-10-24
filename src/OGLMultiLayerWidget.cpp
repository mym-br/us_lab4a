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

#include "OGLMultiLayerWidget.h"

#include <algorithm> /* copy, std::max */
#include <cmath>

#include <QKeyEvent>
#include <QMouseEvent>
#include <QPainter>
#include <QWheelEvent>

#include "Log.h"
#include "Util.h"



namespace {

constexpr int MININUM_WIDTH = 400;
constexpr int MININUM_HEIGHT = 300;
constexpr float MIN_SCALE = 0.1f;
constexpr float MAX_SCALE = 50.0f;
constexpr float MARGIN = 0.005f;
constexpr float ZOOM_IN_FACTOR = 1.1f;
constexpr float ZOOM_OUT_FACTOR = 0.9f;
constexpr float ARROWHEAD_SIZE_FACTOR = 0.015f;
constexpr float ARROW_SIZE_FACTOR = 0.1f;

}

namespace Lab {

OGLMultiLayerWidget::OGLMultiLayerWidget(QWidget* parent)
		: QOpenGLWidget(parent)
		, scale_(1.0)
		, rotX_()
		, rotY_()
		, minX_()
		, maxX_()
		, minY_()
		, maxY_()
		, minZ_()
		, maxZ_()
{
	setMinimumWidth(MININUM_WIDTH);
	setMinimumHeight(MININUM_HEIGHT);
}

void
OGLMultiLayerWidget::updateData(const std::vector<XYZValue<float>>& pointArray, const std::vector<unsigned int>& indexArray)
{
	if (pointArray.empty()) return;
	LOG_DEBUG << "pointArray.size(): " << pointArray.size();
	LOG_DEBUG << "indexArray.size(): " << indexArray.size();
	resetTransformation();

	minX_ = maxX_ = pointArray[0].x;
	minY_ = maxY_ = pointArray[0].z; // z -> y
	minZ_ = maxZ_ = pointArray[0].y; // y -> z

	pointArray_.resize(pointArray.size());
	for (unsigned int i = 0; i < pointArray.size(); ++i) {
		pointArray_[i].pos.x = pointArray[i].x;
		pointArray_[i].pos.y = pointArray[i].z; // z -> y
		pointArray_[i].pos.z = pointArray[i].y; // y -> z
		pointArray_[i].color.red   = 1.0f;
		pointArray_[i].color.green = 1.0f;
		pointArray_[i].color.blue  = 1.0f;
		pointArray_[i].color.alpha = qBound(0.0f, pointArray[i].value, 1.0f);

		if (pointArray_[i].pos.x < minX_) minX_ = pointArray_[i].pos.x;
		else if (pointArray_[i].pos.x > maxX_) maxX_ = pointArray_[i].pos.x;

		if (pointArray_[i].pos.y < minY_) minY_ = pointArray_[i].pos.y;
		else if (pointArray_[i].pos.y > maxY_) maxY_ = pointArray_[i].pos.y;

		if (pointArray_[i].pos.z < minZ_) minZ_ = pointArray_[i].pos.z;
		else if (pointArray_[i].pos.z > maxZ_) maxZ_ = pointArray_[i].pos.z;
	}

	indexArray_ = indexArray;

	revIndexArray_.resize(indexArray_.size());
	std::copy(indexArray_.rbegin(), indexArray_.rend(), revIndexArray_.begin());

	update();
}

void
OGLMultiLayerWidget::paintGL()
{
	QPainter painter(this);
	painter.setFont(QFont("monospace"));

	//==================================================
	painter.beginNativePainting();

	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_COLOR_ARRAY);

	// QPainter disables GL_DEPTH_TEST, GL_BLEND.
	//glEnable(GL_DEPTH_TEST); // not needed
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	// QPainter clears the depth buffer.
	//glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClear(GL_COLOR_BUFFER_BIT);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	const float aspect = static_cast<float>(width()) / height();
	const float yD = std::max((maxX_ - minX_) / aspect, maxY_ - minY_) * (0.5f + MARGIN);
	const float xD = yD * aspect;
	const float dz = (maxZ_ - minZ_) * (0.5f + MARGIN);
	const float zD = scale_ * std::sqrt(xD * xD + yD * yD + dz * dz);
	glOrtho(-xD, xD, -yD, yD, -zD, zD); // l, r, b, t, near, far

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glRotatef(rotX_, 1.0f, 0.0f, 0.0f);
	glRotatef(rotY_, 0.0f, 1.0f, 0.0f);
	glScalef(scale_, scale_, scale_);
	glTranslatef(
		-(minX_ + maxX_) * 0.5f,
		-(minY_ + maxY_) * 0.5f,
		-(minZ_ + maxZ_) * 0.5f);

	// Draw the image layers.
	if (pointArray_.size() >= 2) {
		glVertexPointer(3, GL_FLOAT, sizeof(OGLPoint3DA), &pointArray_[0].pos);
		glColorPointer(4, GL_FLOAT, sizeof(OGLPoint3DA), &pointArray_[0].color);
		// Use the right drawing order for the alpha blending.
		if (std::abs(rotY_) <= 90.0f) {
			glDrawElements(GL_TRIANGLES, indexArray_.size(), GL_UNSIGNED_INT, indexArray_.data());
		} else {
			glDrawElements(GL_TRIANGLES, revIndexArray_.size(), GL_UNSIGNED_INT, revIndexArray_.data());
		}
	}

	glDisable(GL_DEPTH_TEST);

	glColor3f(1.0f, 1.0f, 1.0f);

	// Bottom frame.
	glBegin(GL_LINE_LOOP);
		glVertex3f(minX_, minY_, minZ_);
		glVertex3f(minX_, minY_, maxZ_);
		glVertex3f(maxX_, minY_, maxZ_);
		glVertex3f(maxX_, minY_, minZ_);
	glEnd();
	// Top frame.
	glBegin(GL_LINE_LOOP);
		glVertex3f(minX_, maxY_, minZ_);
		glVertex3f(minX_, maxY_, maxZ_);
		glVertex3f(maxX_, maxY_, maxZ_);
		glVertex3f(maxX_, maxY_, minZ_);
	glEnd();

	const float arrowRef = std::max(std::max(maxX_ - minX_, maxY_ - minY_), maxZ_ - minZ_) * scale_;
	const float arrowSize = ARROW_SIZE_FACTOR * arrowRef;
	const float arrowheadSize = ARROWHEAD_SIZE_FACTOR * arrowRef;
	const float arrowheadSize2 = 0.5f * arrowheadSize;

	// X arrow.
	glBegin(GL_LINES);
		glVertex3f(minX_            , minY_, minZ_);
		glVertex3f(minX_ + arrowSize, minY_, minZ_);
	glEnd();
	glBegin(GL_LINE_LOOP);
		glVertex3f(minX_ + arrowSize                , minY_ + arrowheadSize2, minZ_);
		glVertex3f(minX_ + arrowSize + arrowheadSize, minY_                 , minZ_);
		glVertex3f(minX_ + arrowSize                , minY_ - arrowheadSize2, minZ_);
	glEnd();
	// X label.
	glBegin(GL_LINES);
		glVertex3f(minX_ + arrowSize                , minY_ -        arrowheadSize, minZ_);
		glVertex3f(minX_ + arrowSize + arrowheadSize, minY_ - 2.0f * arrowheadSize, minZ_);

		glVertex3f(minX_ + arrowSize                , minY_ - 2.0f * arrowheadSize, minZ_);
		glVertex3f(minX_ + arrowSize + arrowheadSize, minY_ -        arrowheadSize, minZ_);
	glEnd();

	// Y arrow.
	glBegin(GL_LINES);
		glVertex3f(minX_, minY_            , minZ_);
		glVertex3f(minX_, minY_ + arrowSize, minZ_);
	glEnd();
	glBegin(GL_LINE_LOOP);
		glVertex3f(minX_ + arrowheadSize2, minY_ + arrowSize                , minZ_);
		glVertex3f(minX_                 , minY_ + arrowSize + arrowheadSize, minZ_);
		glVertex3f(minX_ - arrowheadSize2, minY_ + arrowSize                , minZ_);
	glEnd();
	// Y label.
	glBegin(GL_LINE_STRIP);
		glVertex3f(minX_ - 2.0f * arrowheadSize, minY_ + arrowSize +        arrowheadSize, minZ_);
		glVertex3f(minX_ - 1.5f * arrowheadSize, minY_ + arrowSize + 0.5f * arrowheadSize, minZ_);
		glVertex3f(minX_ - 1.5f * arrowheadSize, minY_ + arrowSize                       , minZ_);
	glEnd();
	glBegin(GL_LINES);
		glVertex3f(minX_ - 1.5f * arrowheadSize, minY_ + arrowSize + 0.5f * arrowheadSize, minZ_);
		glVertex3f(minX_ -        arrowheadSize, minY_ + arrowSize +        arrowheadSize, minZ_);
	glEnd();

	// Z arrow.
	glBegin(GL_LINES);
		glVertex3f(minX_, minY_, minZ_);
		glVertex3f(minX_, minY_, minZ_ + arrowSize);
	glEnd();
	glBegin(GL_LINE_LOOP);
		glVertex3f(minX_, minY_ + arrowheadSize2, minZ_ + arrowSize);
		glVertex3f(minX_, minY_                 , minZ_ + arrowSize + arrowheadSize);
		glVertex3f(minX_, minY_ - arrowheadSize2, minZ_ + arrowSize);
	glEnd();
	// Z label.
	glBegin(GL_LINE_STRIP);
		if (std::abs(rotY_) <= 90.0f) {
			glVertex3f(minX_ - 2.0f * arrowheadSize, minY_ -        arrowheadSize, minZ_ + arrowSize);
			glVertex3f(minX_ -        arrowheadSize, minY_ -        arrowheadSize, minZ_ + arrowSize);
			glVertex3f(minX_ - 2.0f * arrowheadSize, minY_ - 2.0f * arrowheadSize, minZ_ + arrowSize);
			glVertex3f(minX_ -        arrowheadSize, minY_ - 2.0f * arrowheadSize, minZ_ + arrowSize);
		} else {
			glVertex3f(minX_ -        arrowheadSize, minY_ -        arrowheadSize, minZ_ + arrowSize);
			glVertex3f(minX_ - 2.0f * arrowheadSize, minY_ -        arrowheadSize, minZ_ + arrowSize);
			glVertex3f(minX_ -        arrowheadSize, minY_ - 2.0f * arrowheadSize, minZ_ + arrowSize);
			glVertex3f(minX_ - 2.0f * arrowheadSize, minY_ - 2.0f * arrowheadSize, minZ_ + arrowSize);
		}
	glEnd();

	glDisableClientState(GL_COLOR_ARRAY);
	glDisableClientState(GL_VERTEX_ARRAY);

	painter.endNativePainting();
	//==================================================

	painter.setRenderHint(QPainter::TextAntialiasing);
	painter.setPen(Qt::black);
	painter.setBackground(QBrush{Qt::white});
	painter.setBackgroundMode(Qt::OpaqueMode);

	painter.drawText(10,  20, QString("rot x ") + QString::number(rotX_));
	painter.drawText(10,  40, QString("rot y ") + QString::number(rotY_));
	painter.drawText(10,  60, QString("zoom  ") + QString::number(scale_));
}

void
OGLMultiLayerWidget::mousePressEvent(QMouseEvent* e)
{
	mouseAnchor_ = e->pos();
}

void
OGLMultiLayerWidget::mouseMoveEvent(QMouseEvent* e)
{
	QPoint diff = e->pos() - mouseAnchor_;
	if (e->buttons() & Qt::LeftButton) {
		mouseAnchor_ = e->pos();

		rotX_ += diff.y() * (180.0f / height());
		Util::clip(rotX_, -89.0f, 89.0f);
		rotY_ += diff.x() * (180.0f / width());
		Util::clipAngleInDegrees180(rotY_);

		update();
	}
}

void
OGLMultiLayerWidget::wheelEvent(QWheelEvent* e)
{
	if (e->delta() > 0) {
		scale_ *= ZOOM_IN_FACTOR;
	} else {
		scale_ *= ZOOM_OUT_FACTOR;
	}
	Util::clip(scale_, MIN_SCALE, MAX_SCALE);

	update();
}

void
OGLMultiLayerWidget::keyPressEvent(QKeyEvent* e)
{
	if (e->key() == Qt::Key_Space) {
		resetTransformation();
		update();
	}
}

void
OGLMultiLayerWidget::keyReleaseEvent(QKeyEvent* /*e*/)
{
}

void
OGLMultiLayerWidget::resetTransformation()
{
	scale_ = 1.0;
	rotX_ = 0.0;
	rotY_ = 0.0;
}

void
OGLMultiLayerWidget::initializeGL()
{
	//LOG_DEBUG << "initializeGL";

	//glShadeModel(GL_FLAT);

	glClearColor(0.0, 0.0, 0.0, 1.0);
}

void
OGLMultiLayerWidget::resizeGL(int /*width*/, int /*height*/)
{
	// paintGL calls glViewport().
	//glViewport(0, 0, width, height);
}

} // namespace Lab
