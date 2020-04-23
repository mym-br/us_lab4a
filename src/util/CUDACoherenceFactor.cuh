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

template<unsigned int N>
__device__
float
calcPCF(float* re, float* im, float factor)
{
	float phi[N];
	float phiAux[N];

	for (unsigned int i = 0; i < N; ++i) {
		phi[i] = atan2f(im[i], re[i]);
	}
	for (unsigned int i = 0; i < N; ++i) {
		phiAux[i] = phi[i] - copysignf(M_PI, phi[i]);
	}

	const float sf = fminf(standardDeviation<N>(phi), standardDeviation<N>(phiAux));
	return fmaxf(0.0, 1.0f - factor * sf);
}

} // namespace Lab

#endif // CUDA_COHERENCE_FACTOR_CUH
