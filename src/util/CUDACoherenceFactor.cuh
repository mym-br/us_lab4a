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

#ifndef CUDA_COHERENCE_FACTOR_CUH
#define CUDA_COHERENCE_FACTOR_CUH

#include "CUDAStatistics.cuh"



namespace Lab {

template<typename TFloat, int N>
__device__
TFloat
calcPCF(TFloat* re, TFloat* im, TFloat factor)
{
	TFloat phi[N];
	TFloat phiAux[N];

	for (int i = 0; i < N; ++i) {
		phi[i] = atan2(im[i], re[i]);
	}
	for (int i = 0; i < N; ++i) {
		phiAux[i] = phi[i] - copysign(TFloat(M_PI), phi[i]);
	}

	const TFloat sf = fmin(standardDeviation<TFloat, N>(phi), standardDeviation<TFloat, N>(phiAux));
	return fmax(TFloat(0), TFloat(1) - factor * sf);
}

} // namespace Lab

#endif // CUDA_COHERENCE_FACTOR_CUH
