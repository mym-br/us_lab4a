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

#ifndef XZVALUEFACTOR_H_
#define XZVALUEFACTOR_H_

namespace Lab {

template<typename T>
struct XZValueFactor {
	T x;
	T z;
	T value;
	T factor;

	XZValueFactor() : x(), z(), value(), factor(1) {}
};

} // namespace Lab

#endif /* XZVALUEFACTOR_H_ */
