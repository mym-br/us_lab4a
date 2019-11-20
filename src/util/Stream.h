/***************************************************************************
 *  Copyright 2019 Marcelo Y. Matuda                                       *
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

#ifndef STREAM_H
#define STREAM_H

#include <ostream>
#include <vector>

namespace Lab {

template<typename T> std::ostream& operator<<(std::ostream& out, const std::vector<T>& v);

template<typename T>
std::ostream&
operator<<(std::ostream& out, const std::vector<T>& v)
{
	out << "{ ";
	if (v.size() > 0) {
		out << v[0];
	}
	for (std::size_t i = 1; i < v.size(); ++i) {
		out << ", " << v[i];
	}
	out << " }";
	return out;
}

} // namespace Lab

#endif // STREAM_H
