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

#ifndef STATISTICS_H_
#define STATISTICS_H_

#include <cmath>
#include <cstddef> /* std::size_t */



namespace Lab {
namespace Statistics {

template<typename TFloat> TFloat standardDeviation(const TFloat* data, std::size_t size);
template<typename TFloat> TFloat arithmeticMean(const TFloat* data, std::size_t size);



template<typename TFloat>
TFloat
standardDeviation(const TFloat* data, std::size_t size)
{
	TFloat sum = 0.0;
	const TFloat* end = data + size;
	const TFloat mean = arithmeticMean(data, size);
	while (data != end) {
		const TFloat e = *data++ - mean;
		sum += e * e;
	}
	return std::sqrt(sum / size);
}

template<typename TFloat>
TFloat
arithmeticMean(const TFloat* data, std::size_t size)
{
	TFloat sum = 0.0;
	const TFloat* end = data + size;
	while (data != end) {
		sum += *data++;
	}
	return sum / size;
}

} // namespace Statistics
} // namespace Lab

#endif /* STATISTICS_H_ */
