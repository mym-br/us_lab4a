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

#ifndef OGLFIGUREWIDGET_H_
#define OGLFIGUREWIDGET_H_

#include <cmath>
#include <vector>

#include <GL/glu.h>

#include <QOpenGLWidget>

#include "global.h"
#include "Matrix2.h"
#include "XY.h"
#include "XZ.h"
#include "XZValue.h"

//QT_BEGIN_NAMESPACE
//class QPaintEvent;
//QT_END_NAMESPACE

#define OGLFIGUREWIDGET_USE_VERTEX_ARRAYS 1



namespace Lab {

struct OGLPos3D {
	float x;
	float y;
	float z;
	OGLPos3D(float ix, float iy, float iz) : x{ix}, y{iy}, z{iz} {}
};

struct OGLColor {
	float red;
	float green;
	float blue;
	OGLColor(float r, float g, float b) : red{r}, green{g}, blue{b} {}
};

struct OGLPoint3D {
	OGLPos3D pos;
	OGLColor color;
	OGLPoint3D() : pos{0.0, 0.0, 0.0}, color{0.0, 0.0, 0.0} {}
	OGLPoint3D(float x, float y, float z, float r, float g, float b)
		: pos{x, y, z}
		, color{r, g, b} {}
};

class OGLFigureWidget : public QOpenGLWidget {
	Q_OBJECT
public:
	OGLFigureWidget(QWidget* parent = 0);
	virtual ~OGLFigureWidget() {}

	void resetScale();
	void setVisualization(Figure::Visualization visualization) { visualization_ = visualization; }
	Figure::Visualization visualization() const { return visualization_; }
	void setMinDecibels(float minDecibels);
	float minDecibels() const { return minDecibels_; }
	void setColormap(Figure::Colormap colormap) { colormap_ = colormap; }
	Figure::Colormap colormap() const { return colormap_; }
	void updateGridData(double valueScale, const Matrix2<XZValue<float>>& gridData);
	void updatePointList(const std::vector<XZ<float>>& pointList);
	void setRotationMode(bool enabled) {
		editingDistanceMarkers_ = false;
		distanceMarker1_ = QPoint(); // clear
		distanceMarker2_ = QPoint(); // clear

		rotationMode_ = enabled;
		update();
	}
	void setShowInfo(bool enabled) {
		showInfo_ = enabled;
		update();
	}
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
	void saveGLState();
	void restoreGLState();

	float calcValueFactor(float valueScale, const Matrix2<XZValue<float>>& data);

	template<typename ColorScale>
	void fillOGLGridDataWithAbsValues(const Matrix2<XZValue<float>>& data, float valueFactor);

	template<typename ColorScale>
	void fillOGLGridDataWithLogAbsValues(const Matrix2<XZValue<float>>& data, float valueFactor);

	template<typename ColorScale> void fillOGLGridData(double valueScale, const Matrix2<XZValue<float>>& gridData);

	void createVertex(const OGLPoint3D& point, float valueScale);
#ifdef OGLFIGUREWIDGET_USE_VERTEX_ARRAYS
	void fillVertexArray(const OGLPoint3D& point, float valueScale);
#endif

	bool valuesInDecibel_;
	bool editingDistanceMarkers_;
	bool rotationMode_;
	bool showPoints_;
	bool showInfo_;
	Figure::Colormap colormap_;
	Figure::Visualization visualization_;
	float minDecibels_;
	float minValue_; // this variable is not used for linear scale
	float scale_;
	float rotX_;
	float rotZ_;
	float minX_;
	float maxX_;
	float minZ_;
	float maxZ_;
	float valueScale_;
	float maxAbsLevel_;
	float maxAbsLevelDecibels_;
	float offsetX_;
	float offsetY_;
	QPoint mouseAnchor_;
	QPoint distanceMarker1_;
	QPoint distanceMarker2_;
	Matrix2<OGLPoint3D> oglGridData_;
	std::vector<XY<float>> pointList_;
#ifdef OGLFIGUREWIDGET_USE_VERTEX_ARRAYS
	std::vector<OGLPos3D> vertexArray_;
	std::vector<OGLColor> colorArray_;
#endif
};

inline
void
OGLFigureWidget::createVertex(const OGLPoint3D& point, float valueScale)
{
	glColor3f(point.color.red, point.color.green, point.color.blue);
	glVertex3f(point.pos.x, point.pos.y, point.pos.z * valueScale);
}

#ifdef OGLFIGUREWIDGET_USE_VERTEX_ARRAYS
inline
void
OGLFigureWidget::fillVertexArray(const OGLPoint3D& point, float valueScale)
{
	vertexArray_.emplace_back(point.pos.x, point.pos.y, point.pos.z * valueScale);
	colorArray_.emplace_back(point.color.red, point.color.green, point.color.blue);
}
#endif

} // namespace Lab

#endif /* OGLFIGUREWIDGET_H_ */
