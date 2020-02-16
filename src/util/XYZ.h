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

#ifndef XYZ_H
#define XYZ_H

#include <iostream>

namespace Lab {

template<typename T>
struct XYZ {
	T x;
	T y;
	T z;

	XYZ() : x(), y(), z() {}
	XYZ(T x, T y, T z) : x(x), y(y), z(z) {}
};

template<typename T>
std::ostream&
operator<<(std::ostream& out, const XYZ<T>& data)
{
	out << data.x << ' ' << data.y << ' ' << data.z;
	return out;
}

} // namespace Lab

#endif // XYZ_H
