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

#define MININUM_WIDTH 400
#define MININUM_HEIGHT 300
#define MIN_SCALE 0.1f
#define MAX_SCALE 50.0f
#define MARGIN 0.005f
#define ZOOM_IN_FACTOR 1.1f
#define ZOOM_OUT_FACTOR 0.9f



namespace Lab {

OGLMultiLayerWidget::OGLMultiLayerWidget(QWidget* parent)
		: QOpenGLWidget{parent}
		, scale_{1.0}
		, rotX_{}
		, rotY_{}
		, minX_{}
		, maxX_{}
		, minY_{}
		, maxY_{}
		, minZ_{}
		, maxZ_{}
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
	const float zD = 2.0f * std::max(std::max(xD, yD), (maxZ_ - minZ_) * (0.5f + MARGIN)) * scale_;
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
	LOG_DEBUG << "initializeGL";

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
