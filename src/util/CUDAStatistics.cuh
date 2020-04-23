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

template<unsigned int N>
__device__
float
arithmeticMean(float* data)
{
	float sum = 0.0;
	for (unsigned int i = 0; i < N; ++i) {
		sum += data[i];
	}
	return sum * (static_cast<float>(1) / N);
}

template<unsigned int N>
__device__
float
standardDeviation(float* data)
{
	float sum = 0.0;
	const float mean = arithmeticMean<N>(data);
	for (unsigned int i = 0; i < N; ++i) {
		const float e = data[i] - mean;
		sum += e * e;
	}
	return sqrtf(sum * (static_cast<float>(1) / N));
}

} // namespace Lab

#endif // CUDA_STATISTICS_CUH
