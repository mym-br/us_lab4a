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

#ifndef CUDA_STATISTICS_CUH
#define CUDA_STATISTICS_CUH

namespace Lab {

__device__
inline
float
arithmeticMean(float* data, unsigned int size)
{
	float sum = 0.0;
	for (unsigned int i = 0; i < size; ++i) {
		sum += data[i];
	}
	return sum / size;
}

__device__
inline
float
standardDeviation(float* data, unsigned int size)
{
	float sum = 0.0;
	const float mean = arithmeticMean(data, size);
	for (unsigned int i = 0; i < size; ++i) {
		const float e = data[i] - mean;
		sum += e * e;
	}
	return sqrtf(sum / size);
}

} // namespace Lab

#endif // CUDA_STATISTICS_CUH
