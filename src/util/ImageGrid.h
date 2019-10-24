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



namespace Lab {

template<typename FloatType>
class ImageGrid {
public:
	template<typename T>
		static void get(ConstParameterMapPtr pm, FloatType lambda, T &grid);
private:
	static constexpr FloatType minLambda = 1.0e-6;

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
	const auto minX    = pm->value<FloatType>("rectangular_min_x"   , -10000.0, 10000.0);
	const auto maxX    = pm->value<FloatType>("rectangular_max_x"   ,     minX, 10000.0);
	const auto minY    = pm->value<FloatType>("rectangular_min_y"   , -10000.0, 10000.0);
	const auto maxY    = pm->value<FloatType>("rectangular_max_y"   ,     minY, 10000.0);
	const auto minZ    = pm->value<FloatType>("rectangular_min_z"   , -10000.0, 10000.0);
	const auto maxZ    = pm->value<FloatType>("rectangular_max_z"   ,     minZ, 10000.0);
	const auto stepDiv = pm->value<FloatType>("rectangular_step_div",   1.0e-2,  1000.0);
	const auto originX = pm->value<FloatType>("rectangular_origin_x", -10000.0, 10000.0);
	const auto originY = pm->value<FloatType>("rectangular_origin_y", -10000.0, 10000.0);
	const auto originZ = pm->value<FloatType>("rectangular_origin_z", -10000.0, 10000.0);

	const FloatType maxStep = lambda / stepDiv;
	const FloatType dx = maxX - minX;
	const FloatType dy = maxY - minY;
	const FloatType dz = maxZ - minZ;

	if (dx == 0.0) { // z-y
		std::vector<FloatType> zList;
		std::vector<FloatType> yList;
		Util::fillSequenceFromStartToEndWithMaximumStep(zList, minZ, maxZ, maxStep);
		Util::fillSequenceFromStartToEndWithMaximumStep(yList, minY, maxY, maxStep);
		grid.resize(zList.size(), yList.size());
		const FloatType x = minX + originX;
		for (std::size_t i = 0; i < zList.size(); ++i) {
			for (std::size_t j = 0; j < yList.size(); ++j) {
				grid(i, j).x = x;
				grid(i, j).y = yList[j] + originY;
				grid(i, j).z = zList[i] + originZ;
			}
		}
	} else if (dy == 0.0) { // x-z
		std::vector<FloatType> xList;
		std::vector<FloatType> zList;
		Util::fillSequenceFromStartToEndWithMaximumStep(xList, minX, maxX, maxStep);
		Util::fillSequenceFromStartToEndWithMaximumStep(zList, minZ, maxZ, maxStep);
		grid.resize(xList.size(), zList.size());
		const FloatType y = minY + originY;
		for (std::size_t i = 0; i < xList.size(); ++i) {
			for (std::size_t j = 0; j < zList.size(); ++j) {
				grid(i, j).x = xList[i] + originX;
				grid(i, j).y = y;
				grid(i, j).z = zList[j] + originZ;
			}
		}
	} else if (dz == 0.0) { // x-y
		std::vector<FloatType> xList;
		std::vector<FloatType> yList;
		Util::fillSequenceFromStartToEndWithMaximumStep(xList, minX, maxX, maxStep);
		Util::fillSequenceFromStartToEndWithMaximumStep(yList, minY, maxY, maxStep);
		grid.resize(xList.size(), yList.size());
		const FloatType z = minZ + originZ;
		for (std::size_t i = 0; i < xList.size(); ++i) {
			for (std::size_t j = 0; j < yList.size(); ++j) {
				grid(i, j).x = xList[i] + originX;
				grid(i, j).y = yList[j] + originY;
				grid(i, j).z = z;
			}
		}
	} else {
		THROW_EXCEPTION(InvalidParameterException, "Invalid grid configuration.");
	}
}

template<typename FloatType>
template<typename T>
void
ImageGrid<FloatType>::getSectorialGrid(ConstParameterMapPtr pm, FloatType lambda, T& grid)
{
	const auto minRadius     = pm->value<FloatType>("sectorial_min_radius"     ,             1.0e-4, 10000.0);
	const auto maxRadius     = pm->value<FloatType>("sectorial_max_radius"     , minRadius + 1.0e-6, 10000.0);
	const auto radiusStepDiv = pm->value<FloatType>("sectorial_radius_step_div",             1.0e-2,  1000.0);
	const auto minAngle      = pm->value<FloatType>("sectorial_min_angle"      ,              -88.0,    88.0);
	const auto maxAngle      = pm->value<FloatType>("sectorial_max_angle"      ,  minAngle + 1.0e-3,    88.0);
	const auto angleStep     = pm->value<FloatType>("sectorial_angle_step"     ,             1.0e-6,    10.0);
	const auto originX       = pm->value<FloatType>("sectorial_origin_x"       ,           -10000.0, 10000.0);
	const auto originZ       = pm->value<FloatType>("sectorial_origin_z"       ,           -10000.0, 10000.0);

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
			grid(i, j).y = 0.0;
			grid(i, j).z = r * std::cos(angle) + originZ;
		}
	}
}

template<typename FloatType>
template<typename T>
void
ImageGrid<FloatType>::get(ConstParameterMapPtr pm, FloatType lambda, T& grid)
{
	if (lambda < minLambda) {
		THROW_EXCEPTION(InvalidParameterException, "Lambda is too small: " << lambda << '.');
	}

	const auto gridType = pm->value<std::string>("image_grid_type");
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
