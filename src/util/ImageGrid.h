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

#ifndef IMAGEGRID_H
#define IMAGEGRID_H

#include <cmath>
#include <vector>

#include "Exception.h"
#include "Log.h"
#include "ParameterMap.h"
#include "Util.h"

#define IMAGE_GRID_MIN_LAMBDA (1.0e-6)



namespace Lab {

template<typename FloatType>
class ImageGrid {
public:
	template<typename T>
		static void get(ConstParameterMapPtr pm, FloatType lambda, T &grid);
private:
	ImageGrid();
	~ImageGrid();

	template<typename T>
		static void getRectangularGrid(ConstParameterMapPtr pm, FloatType lambda, T& grid);
	template<typename T>
		static void getSectorialGrid(ConstParameterMapPtr pm, FloatType lambda, T& grid);
};

template<typename FloatType>
template<typename T>
void
ImageGrid<FloatType>::getRectangularGrid(ConstParameterMapPtr pm, FloatType lambda, T& grid)
{
	const FloatType minX     = pm->value<FloatType>("rectangular_min_x"     ,     -500.0e-3, 500.0e-3);
	const FloatType maxX     = pm->value<FloatType>("rectangular_max_x"     , minX + 1.0e-4, 500.0e-3);
	const FloatType xStepDiv = pm->value<FloatType>("rectangular_x_step_div",        1.0e-2,   1000.0);
	const FloatType minZ     = pm->value<FloatType>("rectangular_min_z"     ,     -500.0e-3, 500.0e-3);
	const FloatType maxZ     = pm->value<FloatType>("rectangular_max_z"     , minZ + 1.0e-4, 500.0e-3);
	const FloatType zStepDiv = pm->value<FloatType>("rectangular_z_step_div",        1.0e-2,   1000.0);
	const FloatType originX  = pm->value<FloatType>("rectangular_origin_x"  ,     -200.0e-3, 200.0e-3);
	const FloatType originZ  = pm->value<FloatType>("rectangular_origin_z"  ,     -200.0e-3, 200.0e-3);

	std::vector<FloatType> xList;
	Util::fillSequenceFromStartToEndWithMaximumStep(xList, minX, maxX, lambda / xStepDiv);

	std::vector<FloatType> zList;
	Util::fillSequenceFromStartToEndWithMaximumStep(zList, minZ, maxZ, lambda / zStepDiv);

	LOG_DEBUG << "xList.size(): " << xList.size() << " zList.size(): " << zList.size();

	grid.resize(xList.size(), zList.size());
	for (std::size_t i = 0; i < xList.size(); ++i) {
		for (std::size_t j = 0; j < zList.size(); ++j) {
			grid(i, j).x = xList[i] + originX;
			grid(i, j).z = zList[j] + originZ;
		}
	}
}

template<typename FloatType>
template<typename T>
void
ImageGrid<FloatType>::getSectorialGrid(ConstParameterMapPtr pm, FloatType lambda, T& grid)
{
	const FloatType minRadius         = pm->value<FloatType>("sectorial_min_radius"     ,             1.0e-4,   100.0e-3);
	const FloatType maxRadius         = pm->value<FloatType>("sectorial_max_radius"     , minRadius + 1.0e-4, 10000.0e-3);
	const FloatType radiusStepDiv     = pm->value<FloatType>("sectorial_radius_step_div",             1.0e-2,      256.0);
	const FloatType minAngle          = pm->value<FloatType>("sectorial_min_angle"      ,              -80.0,       80.0);
	const FloatType maxAngle          = pm->value<FloatType>("sectorial_max_angle"      ,  minAngle + 1.0e-3,       80.0);
	const FloatType angleStep         = pm->value<FloatType>("sectorial_angle_step"     ,             1.0e-5,       10.0);
	const FloatType originX           = pm->value<FloatType>("sectorial_origin_x"       ,          -200.0e-3,   200.0e-3);
	const FloatType originZ           = pm->value<FloatType>("sectorial_origin_z"       ,          -200.0e-3,   200.0e-3);

	std::vector<FloatType> radiusList;
	Util::fillSequenceFromStartToEndWithMaximumStep(radiusList,
			minRadius, maxRadius, lambda / radiusStepDiv);

	std::vector<FloatType> angleList;
	Util::fillSequenceFromStartToEndWithMaximumStep(angleList,
			Util::degreeToRadian(minAngle), Util::degreeToRadian(maxAngle), Util::degreeToRadian(angleStep));

	LOG_DEBUG << "radiusList.size(): " << radiusList.size() << " angleList.size(): " << angleList.size();

	grid.resize(angleList.size(), radiusList.size());
	for (std::size_t i = 0; i < angleList.size(); ++i) {
		for (std::size_t j = 0; j < radiusList.size(); ++j) {
			const FloatType angle = angleList[i];
			const FloatType r = radiusList[j];
			grid(i, j).x = r * std::sin(angle) + originX;
			grid(i, j).z = r * std::cos(angle) + originZ;
		}
	}
}

template<typename FloatType>
template<typename T>
void
ImageGrid<FloatType>::get(ConstParameterMapPtr pm, FloatType lambda, T& grid)
{
	if (lambda < IMAGE_GRID_MIN_LAMBDA) {
		THROW_EXCEPTION(InvalidParameterException, "Lambda is too small: " << lambda << '.');
	}

	std::string gridType = pm->value<std::string>("image_grid_type");
	if (gridType == "rectangular") {
		getRectangularGrid(pm, lambda, grid);
	} else if (gridType == "sectorial") {
		getSectorialGrid(pm, lambda, grid);
	} else {
		THROW_EXCEPTION(InvalidParameterException, "Invalid grid type: " << gridType << '.');
	}
}

} // namespace Lab

#endif // IMAGEGRID_H
