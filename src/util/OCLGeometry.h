/***************************************************************************
 *  Copyright 2020 Marcelo Y. Matuda                                       *
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

#ifndef OCL_GEOMETRY_H
#define OCL_GEOMETRY_H

#include <string>

namespace Lab {
namespace OCLGeometry {

inline
std::string
code() {
	return R"CLC(

MFloat
distance2D(MFloat x0, MFloat y0, MFloat x1, MFloat y1)
{
	const MFloat dx = x1 - x0;
	const MFloat dy = y1 - y0;
	return sqrt(dx * dx + dy * dy);
}

MFloat
distance2DY0(MFloat x0, MFloat x1, MFloat y1)
{
	const MFloat dx = x1 - x0;
	return sqrt(dx * dx + y1 * y1);
}

)CLC";
}

} // namespace OCLGeometry
} // namespace Lab

#endif // OCL_GEOMETRY_H
