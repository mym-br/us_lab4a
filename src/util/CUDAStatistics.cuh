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

template<typename TFloat, int N>
__device__
TFloat
arithmeticMean(TFloat* data)
{
	TFloat sum = 0;
	for (int i = 0; i < N; ++i) {
		sum += data[i];
	}
	return sum * (TFloat(1) / N);
}

template<typename TFloat, int N>
__device__
TFloat
standardDeviation(TFloat* data)
{
	TFloat sum = 0;
	const TFloat mean = arithmeticMean<TFloat, N>(data);
	for (int i = 0; i < N; ++i) {
		const TFloat e = data[i] - mean;
		sum += e * e;
	}
	return sqrt(sum * (TFloat(1) / N));
}

} // namespace Lab

#endif // CUDA_STATISTICS_CUH
