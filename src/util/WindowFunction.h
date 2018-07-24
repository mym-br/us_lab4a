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

#ifndef WINDOWFUNCTION_H_
#define WINDOWFUNCTION_H_

#include <vector>

#include "Util.h"



namespace Lab {
namespace WindowFunction {

template<typename FloatType>
void
hamming(int n, std::vector<FloatType>& w)
{
	w.resize(n);
	for (int i = 0; i < n; ++i) {
		w[i] = 0.54 - 0.46 * std::cos((2.0 * PI * i) / (n - 1));
	}
}

} // namespace WindowFunction
} // namespace Lab

#endif /* WINDOWFUNCTION_H_ */
