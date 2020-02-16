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

#ifndef XYZVALUE_H
#define XYZVALUE_H

namespace Lab {

template<typename T>
struct XYZValue {
	T x;
	T y;
	T z;
	T value;

	XYZValue() : x(), y(), z(), value() {}
	XYZValue(T x, T y, T z, T value) : x(x), y(y), z(z), value(value) {}

	bool operator==(const XYZValue<T>& o) const {
		return x == o.x && y == o.y && z == o.z && value == o.value;
	}
};

} // namespace Lab

#endif // XYZVALUE_H
