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

#include <vector>

#include <QOpenGLWidget>

#include "global.h"
#include "Matrix.h"
#include "OGL.h"
#include "XYZ.h"
#include "XYZValue.h"

#define OGLFIGUREWIDGET_USE_VERTEX_ARRAY 1



namespace Lab {

class OGLFigureWidget : public QOpenGLWidget {
	Q_OBJECT
public:
	OGLFigureWidget(QWidget* parent=nullptr);
	virtual ~OGLFigureWidget() {}

	void setVisualization(Figure::Visualization visualization);
	Figure::Visualization visualization() const { return visualization_; }
	void setMinDecibels(float minDecibels);
	float minDecibels() const { return minDecibels_; }
	void setColormap(Figure::Colormap colormap);
	Figure::Colormap colormap() const { return colormap_; }
	void updateData(float dataValueScale,
			const Matrix<XYZValue<float>>* gridData,
			const std::vector<XYZ<float>>* pointList);
	void updateDataVisualization();
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
	struct Vector {
		float x;
		float y;
		float z;
	};

	void updateGridData(float dataValueScale, const Matrix<XYZValue<float>>& gridData);
	void updatePointList(const std::vector<XYZ<float>>& pointList);

	float calcValueFactor(const Matrix<XYZValue<float>>& data);

	template<typename ColorScale>
	void fillOGLGridDataWithAbsValues(const Matrix<XYZValue<float>>& data, float valueFactor);

	template<typename ColorScale>
	void fillOGLGridDataWithLogAbsValues(const Matrix<XYZValue<float>>& data, float valueFactor);

	template<typename ColorScale> void fillOGLGridData();

#ifdef OGLFIGUREWIDGET_USE_VERTEX_ARRAY
	void fillIndexArray(unsigned int iA, unsigned int iB, std::vector<GLuint>& indexArray);
#else
	void createVertex(const OGLPoint3D& point);
#endif

	bool valuesInDecibel_;
	bool editingDistanceMarkers_;
	bool rotationMode_;
	bool showPoints_;
	bool showInfo_;
	bool dataChanged_;
	Figure::Colormap colormap_;
	Figure::Visualization visualization_;
	float minDecibels_;
	float minValue_; // this variable is not used for linear scale
	float scale_;
	float rotX_;
	float rotZ_;
	float minX_;
	float maxX_;
	float minY_;
	float maxY_;
	float minZ_;
	float maxZ_;
	float dataValueScale_;
	float valueScale_;
	float maxAbsLevel_;
	float maxAbsLevelDecibels_;
	float offsetX_;
	float offsetY_;
	Vector normal_;
	QPoint mouseAnchor_;
	QPoint distanceMarker1_;
	QPoint distanceMarker2_;
	Matrix<XYZValue<float>> gridData_;
	std::vector<XYZ<float>> pointList_;
	Matrix<OGLPoint3D> oglGridData_;
#ifdef OGLFIGUREWIDGET_USE_VERTEX_ARRAY
	std::vector<GLuint> evenIndexArray_;
	std::vector<GLuint> oddIndexArray_;
#else
	GLuint oglDisplayList_;
#endif
};

#ifndef OGLFIGUREWIDGET_USE_VERTEX_ARRAY
inline
void
OGLFigureWidget::createVertex(const OGLPoint3D& point)
{
	glColor3f(point.color.red, point.color.green, point.color.blue);
	glVertex3f(point.pos.x, point.pos.y, point.pos.z);
}
#endif

} // namespace Lab

#endif /* OGLFIGUREWIDGET_H_ */
