/***************************************************************************
 *  Copyright 2014, 2017, 2018, 2020 Marcelo Y. Matuda                     *
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

#ifndef XZ_H_
#define XZ_H_

#include <iostream>



namespace Lab {

template<typename T>
struct XZ {
	T x;
	T z;

	XZ() : x(), z() {}
	XZ(T x, T z) : x(x), z(z) {}
};

template<typename T>
std::ostream&
operator<<(std::ostream& out, const XZ<T>& data)
{
	out << data.x << ' ' << data.z;
	return out;
}

} // namespace Lab

#endif /* XZ_H_ */
