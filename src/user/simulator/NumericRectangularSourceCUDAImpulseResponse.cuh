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

namespace Lab {

template<typename TFloat>
__global__
void
numericSourceIRKernel(
		unsigned int numSubElem,
		TFloat x,
		TFloat y,
		TFloat z,
		TFloat k1,
		TFloat k2,
		const TFloat* subElemX,
		const TFloat* subElemY,
		unsigned int* n0,
		TFloat* value)
{
	const unsigned int subElemIdx = blockIdx.x * blockDim.x + threadIdx.x;
	if (subElemIdx >= numSubElem) {
		// Get data from the first sub-element, to help min/max(n0).
		const TFloat dx = x - subElemX[0];
		const TFloat dy = y - subElemY[0];
		const TFloat r = sqrt(dx * dx + dy * dy + z * z);
		n0[subElemIdx] = rint(r * k1);
		return;
	}

	const TFloat dx = x - subElemX[subElemIdx];
	const TFloat dy = y - subElemY[subElemIdx];
	const TFloat r = sqrt(dx * dx + dy * dy + z * z);
	n0[subElemIdx] = rint(r * k1);
	value[subElemIdx] = k2 / r;
}

template<typename TFloat>
__global__
void
accumulateIRSamplesKernel(
		unsigned int numSubElem,
		unsigned int minN0,
		const unsigned int* n0,
		const TFloat* value,
		TFloat* h)
{
	const unsigned int subElemIdx = blockIdx.x * blockDim.x + threadIdx.x;
	if (subElemIdx >= numSubElem) {
		return;
	}

	// Different sub-elements may have the same value for n0.
	atomicAdd(h + n0[subElemIdx] - minN0, value[subElemIdx]);
}

} // namespace Lab
