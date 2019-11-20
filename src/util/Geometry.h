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
#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <cmath>



namespace Lab {
namespace Geometry {

template<typename FloatType> FloatType distance2D(FloatType x0, FloatType y0,
							FloatType x1, FloatType y1);
template<typename FloatType> FloatType distance2DY0(FloatType x0,
							FloatType x1, FloatType y1);
template<typename FloatType> FloatType distance3D(FloatType x0, FloatType y0, FloatType z0,
							FloatType x1, FloatType y1, FloatType z1);
template<typename FloatType> FloatType distance3DZ0(FloatType x0, FloatType y0,
							FloatType x1, FloatType y1, FloatType z1);

//==============================================================================

template<typename FloatType>
FloatType
distance2D(FloatType x0, FloatType y0, FloatType x1, FloatType y1)
{
	const FloatType dx = x1 - x0;
	const FloatType dy = y1 - y0;
	return std::sqrt(dx * dx + dy * dy);
}

template<typename FloatType>
FloatType
distance2DY0(FloatType x0, FloatType x1, FloatType y1)
{
	const FloatType dx = x1 - x0;
	return std::sqrt(dx * dx + y1 * y1);
}

template<typename FloatType>
FloatType
distance3D(FloatType x0, FloatType y0, FloatType z0, FloatType x1, FloatType y1, FloatType z1)
{
	const FloatType dx = x1 - x0;
	const FloatType dy = y1 - y0;
	const FloatType dz = z1 - z0;
	return std::sqrt(dx * dx + dy * dy + dz * dz);
}

template<typename FloatType>
FloatType
distance3DZ0(FloatType x0, FloatType y0, FloatType x1, FloatType y1, FloatType z1)
{
	const FloatType dx = x1 - x0;
	const FloatType dy = y1 - y0;
	return std::sqrt(dx * dx + dy * dy + z1 * z1);
}

} // namespace Geometry
} // namespace Lab

#endif // GEOMETRY_H
