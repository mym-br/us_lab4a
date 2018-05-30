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

#ifndef OGLMULTILAYERWIDGET_H
#define OGLMULTILAYERWIDGET_H

#include <vector>

#include <QOpenGLWidget>

#include "OGL.h"
#include "XYZValue.h"



namespace Lab {

class OGLMultiLayerWidget : public QOpenGLWidget {
	Q_OBJECT
public:
	OGLMultiLayerWidget(QWidget* parent = 0);
	virtual ~OGLMultiLayerWidget() {}

	void updateData(const std::vector<XYZValue<float>>& pointArray, const std::vector<unsigned int>& indexArray);
protected:
	virtual void initializeGL();
	virtual void paintGL();
	virtual void resizeGL(int width, int height);
	virtual void mousePressEvent(QMouseEvent* e);
	virtual void mouseMoveEvent(QMouseEvent* e);
	virtual void wheelEvent(QWheelEvent* e);
	virtual void keyPressEvent(QKeyEvent* e);
	virtual void keyReleaseEvent(QKeyEvent* e);
private:
	void resetTransformation();

	float scale_;
	float rotX_;
	float rotY_;
	float minX_;
	float maxX_;
	float minY_;
	float maxY_;
	float minZ_;
	float maxZ_;
	QPoint mouseAnchor_;
	std::vector<OGLPoint3DA> pointArray_;
	std::vector<GLuint> indexArray_;
};

} // namespace Lab

#endif // OGLMULTILAYERWIDGET_H
