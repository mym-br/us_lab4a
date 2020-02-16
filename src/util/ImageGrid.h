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
#include "Matrix.h"
#include "ParameterMap.h"
#include "TemplateUtil.h"
#include "Util.h"



namespace Lab {

template<typename TFloat>
class ImageGrid {
public:
	template<typename T>
		static void get(const ParameterMap& pm, TFloat lambda, Matrix<T>& grid);
private:
	static constexpr TFloat minLambda = 1.0e-6;

	ImageGrid() = delete;

	template<typename T>
		static void getRectangularGrid(const ParameterMap& pm, TFloat lambda, Matrix<T>& grid);
	template<typename T>
		static void getSectorialGrid(const ParameterMap& pm, TFloat lambda, Matrix<T>& grid);
};

template<typename TFloat>
template<typename T>
void
ImageGrid<TFloat>::getRectangularGrid(const ParameterMap& pm, TFloat lambda, Matrix<T>& grid)
{
	if constexpr (has_y_member<T>::value) {

		const auto minX    = pm.value<TFloat>("rectangular_min_x"   , -10000.0, 10000.0);
		const auto maxX    = pm.value<TFloat>("rectangular_max_x"   ,     minX, 10000.0);
		const auto minY    = pm.value<TFloat>("rectangular_min_y"   , -10000.0, 10000.0);
		const auto maxY    = pm.value<TFloat>("rectangular_max_y"   ,     minY, 10000.0);
		const auto minZ    = pm.value<TFloat>("rectangular_min_z"   , -10000.0, 10000.0);
		const auto maxZ    = pm.value<TFloat>("rectangular_max_z"   ,     minZ, 10000.0);
		const auto stepDiv = pm.value<TFloat>("rectangular_step_div",   1.0e-2,  1000.0);
		const auto originX = pm.value<TFloat>("rectangular_origin_x", -10000.0, 10000.0);
		const auto originY = pm.value<TFloat>("rectangular_origin_y", -10000.0, 10000.0);
		const auto originZ = pm.value<TFloat>("rectangular_origin_z", -10000.0, 10000.0);

		const TFloat maxStep = lambda / stepDiv;
		const TFloat dx = maxX - minX;
		const TFloat dy = maxY - minY;
		const TFloat dz = maxZ - minZ;

		if (dx == 0.0) { // z-y
			std::vector<TFloat> zList;
			std::vector<TFloat> yList;
			Util::fillSequenceFromStartToEndWithMaximumStep(zList, minZ, maxZ, maxStep);
			Util::fillSequenceFromStartToEndWithMaximumStep(yList, minY, maxY, maxStep);
			grid.resize(zList.size(), yList.size());
			const TFloat x = minX + originX;
			for (std::size_t i = 0; i < zList.size(); ++i) {
				for (std::size_t j = 0; j < yList.size(); ++j) {
					grid(i, j).x = x;
					grid(i, j).y = yList[j] + originY;
					grid(i, j).z = zList[i] + originZ;
				}
			}
		} else if (dy == 0.0) { // x-z
			std::vector<TFloat> xList;
			std::vector<TFloat> zList;
			Util::fillSequenceFromStartToEndWithMaximumStep(xList, minX, maxX, maxStep);
			Util::fillSequenceFromStartToEndWithMaximumStep(zList, minZ, maxZ, maxStep);
			grid.resize(xList.size(), zList.size());
			const TFloat y = minY + originY;
			for (std::size_t i = 0; i < xList.size(); ++i) {
				for (std::size_t j = 0; j < zList.size(); ++j) {
					grid(i, j).x = xList[i] + originX;
					grid(i, j).y = y;
					grid(i, j).z = zList[j] + originZ;
				}
			}
		} else if (dz == 0.0) { // x-y
			std::vector<TFloat> xList;
			std::vector<TFloat> yList;
			Util::fillSequenceFromStartToEndWithMaximumStep(xList, minX, maxX, maxStep);
			Util::fillSequenceFromStartToEndWithMaximumStep(yList, minY, maxY, maxStep);
			grid.resize(xList.size(), yList.size());
			const TFloat z = minZ + originZ;
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

	} else {

		const auto minX    = pm.value<TFloat>("rectangular_min_x"   , -10000.0, 10000.0);
		const auto maxX    = pm.value<TFloat>("rectangular_max_x"   ,     minX, 10000.0);
		const auto minZ    = pm.value<TFloat>("rectangular_min_z"   , -10000.0, 10000.0);
		const auto maxZ    = pm.value<TFloat>("rectangular_max_z"   ,     minZ, 10000.0);
		const auto stepDiv = pm.value<TFloat>("rectangular_step_div",   1.0e-2,  1000.0);
		const auto originX = pm.value<TFloat>("rectangular_origin_x", -10000.0, 10000.0);
		const auto originZ = pm.value<TFloat>("rectangular_origin_z", -10000.0, 10000.0);

		const TFloat maxStep = lambda / stepDiv;

		// x-z
		std::vector<TFloat> xList;
		std::vector<TFloat> zList;
		Util::fillSequenceFromStartToEndWithMaximumStep(xList, minX, maxX, maxStep);
		Util::fillSequenceFromStartToEndWithMaximumStep(zList, minZ, maxZ, maxStep);
		grid.resize(xList.size(), zList.size());
		for (std::size_t i = 0; i < xList.size(); ++i) {
			for (std::size_t j = 0; j < zList.size(); ++j) {
				grid(i, j).x = xList[i] + originX;
				grid(i, j).z = zList[j] + originZ;
			}
		}

	}
}

template<typename TFloat>
template<typename T>
void
ImageGrid<TFloat>::getSectorialGrid(const ParameterMap& pm, TFloat lambda, Matrix<T>& grid)
{
	const auto minRadius     = pm.value<TFloat>("sectorial_min_radius"     ,             1.0e-4, 10000.0);
	const auto maxRadius     = pm.value<TFloat>("sectorial_max_radius"     , minRadius + 1.0e-6, 10000.0);
	const auto radiusStepDiv = pm.value<TFloat>("sectorial_radius_step_div",             1.0e-2,  1000.0);
	const auto minAngle      = pm.value<TFloat>("sectorial_min_angle"      ,              -88.0,    88.0);
	const auto maxAngle      = pm.value<TFloat>("sectorial_max_angle"      ,  minAngle + 1.0e-3,    88.0);
	const auto angleStep     = pm.value<TFloat>("sectorial_angle_step"     ,             1.0e-6,    10.0);
	const auto originX       = pm.value<TFloat>("sectorial_origin_x"       ,           -10000.0, 10000.0);
	const auto originZ       = pm.value<TFloat>("sectorial_origin_z"       ,           -10000.0, 10000.0);

	std::vector<TFloat> radiusList;
	Util::fillSequenceFromStartToEndWithMaximumStep(radiusList,
			minRadius, maxRadius, lambda / radiusStepDiv);

	std::vector<TFloat> angleList;
	Util::fillSequenceFromStartToEndWithMaximumStep(angleList,
			Util::degreeToRadian(minAngle), Util::degreeToRadian(maxAngle), Util::degreeToRadian(angleStep));

	LOG_DEBUG << "radiusList.size(): " << radiusList.size() << " angleList.size(): " << angleList.size();

	grid.resize(angleList.size(), radiusList.size());
	for (std::size_t i = 0; i < angleList.size(); ++i) {
		for (std::size_t j = 0; j < radiusList.size(); ++j) {
			const TFloat angle = angleList[i];
			const TFloat r = radiusList[j];
			grid(i, j).x = r * std::sin(angle) + originX;
			if constexpr (has_y_member<T>::value) {
				grid(i, j).y = 0.0;
			}
			grid(i, j).z = r * std::cos(angle) + originZ;
		}
	}
}

template<typename TFloat>
template<typename T>
void
ImageGrid<TFloat>::get(const ParameterMap& pm, TFloat lambda, Matrix<T>& grid)
{
	if (lambda < minLambda) {
		THROW_EXCEPTION(InvalidParameterException, "Lambda is too small: " << lambda << '.');
	}

	const auto gridType = pm.value<std::string>("image_grid_type");
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
