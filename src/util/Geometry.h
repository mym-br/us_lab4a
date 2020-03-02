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

#if USE_SIMD
# include "SIMD.h"
#endif



namespace Lab {
namespace Geometry {

template<typename TFloat> TFloat distance2D(TFloat x0, TFloat y0,
							TFloat x1, TFloat y1);
template<typename TFloat> TFloat distance2DY0(TFloat x0,
							TFloat x1, TFloat y1);
template<typename TFloat> TFloat distance3D(TFloat x0, TFloat y0, TFloat z0,
							TFloat x1, TFloat y1, TFloat z1);
template<typename TFloat> TFloat distance3DZ0(TFloat x0, TFloat y0,
							TFloat x1, TFloat y1, TFloat z1);

//==============================================================================

template<typename TFloat>
TFloat
distance2D(TFloat x0, TFloat y0, TFloat x1, TFloat y1)
{
#if USE_SIMD
	return SIMD::calcDistance(x0, y0, x1, y1);
#else
	const TFloat dx = x1 - x0;
	const TFloat dy = y1 - y0;
	return std::sqrt(dx * dx + dy * dy);
#endif
}

template<typename TFloat>
TFloat
distance2DY0(TFloat x0, TFloat x1, TFloat y1)
{
#if USE_SIMD
	return SIMD::calcDistance(x0, TFloat(0), x1, y1);
#else
	const TFloat dx = x1 - x0;
	return std::sqrt(dx * dx + y1 * y1);
#endif
}

template<typename TFloat>
TFloat
distance3D(TFloat x0, TFloat y0, TFloat z0, TFloat x1, TFloat y1, TFloat z1)
{
	const TFloat dx = x1 - x0;
	const TFloat dy = y1 - y0;
	const TFloat dz = z1 - z0;
	return std::sqrt(dx * dx + dy * dy + dz * dz);
}

template<typename TFloat>
TFloat
distance3DZ0(TFloat x0, TFloat y0, TFloat x1, TFloat y1, TFloat z1)
{
	const TFloat dx = x1 - x0;
	const TFloat dy = y1 - y0;
	return std::sqrt(dx * dx + dy * dy + z1 * z1);
}

} // namespace Geometry
} // namespace Lab

#endif // GEOMETRY_H
