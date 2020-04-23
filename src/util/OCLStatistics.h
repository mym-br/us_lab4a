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

#ifndef OCL_STATISTICS_H
#define OCL_STATISTICS_H

#include <string>

namespace Lab {
namespace OCLStatistics {

inline
std::string
code() {
	return R"CLC(

MFloat
arithmeticMean(MFloat* data)
{
	MFloat sum = 0;
	for (unsigned int i = 0; i < STATISTICS_N; ++i) {
		sum += data[i];
	}
	return sum * ((MFloat) 1 / STATISTICS_N);
}

MFloat
standardDeviation(MFloat* data)
{
	MFloat sum = 0;
	MFloat mean = arithmeticMean(data);
	for (unsigned int i = 0; i < STATISTICS_N; ++i) {
		const MFloat e = data[i] - mean;
		sum += e * e;
	}
	return sqrt(sum * ((MFloat) 1 / STATISTICS_N));
}

)CLC";
}

} // namespace OCLStatistics
} // namespace Lab

#endif // OCL_STATISTICS_H
