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

#include "CUDAReduce.cuh"



namespace Lab {

__global__
void
blockReduceMinMaxKernel(const unsigned int* in, unsigned int* outMin, unsigned int* outMax, unsigned int n) {
	unsigned int minValue = 0xffffffff;
	unsigned int maxValue = 0;

	const unsigned int totalNumThreads = gridDim.x * blockDim.x;
	for (unsigned int i = blockIdx.x * blockDim.x + threadIdx.x; i < n; i += totalNumThreads) {
		minValue = min(minValue, in[i]);
		maxValue = max(maxValue, in[i]);
	}
	blockReduceMinMax(minValue, maxValue);
	if (threadIdx.x == 0) {
		outMin[blockIdx.x] = minValue; // one value for each block
		outMax[blockIdx.x] = maxValue;
	}
}

__global__
void
reduceMinMaxKernel(unsigned int* blockMinVal, unsigned int* blockMaxVal, unsigned int numBlocks) {
	unsigned int minValue = 0xffffffff;
	unsigned int maxValue = 0;

	const unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx != 0) return; // single thread

	for (int i = 0; i < numBlocks; ++i) {
		minValue = min(minValue, blockMinVal[i]);
		maxValue = max(maxValue, blockMaxVal[i]);
	}

	blockMinVal[0] = minValue;
	blockMaxVal[0] = maxValue;
}

//-----------------------------------------------------------------------------

__global__
void
blockReduceMinKernel(const unsigned int* in, unsigned int* out, unsigned int n) {
	unsigned int minValue = 0xffffffff;

	const unsigned int totalNumThreads = gridDim.x * blockDim.x;
	for (unsigned int i = blockIdx.x * blockDim.x + threadIdx.x; i < n; i += totalNumThreads) {
		minValue = min(minValue, in[i]);
	}
	minValue = blockReduceMin(minValue);
	if (threadIdx.x == 0) {
		out[blockIdx.x] = minValue; // one value for each block
	}
}

__global__
void
reduceMinKernel(unsigned int* blockVal, unsigned int numBlocks) {
	unsigned int minValue = 0xffffffff;

	const unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx != 0) return; // single thread

	for (int i = 0; i < numBlocks; ++i) {
		minValue = min(minValue, blockVal[i]);
	}

	blockVal[0] = minValue;
}

//-----------------------------------------------------------------------------

__global__
void
blockReduceMaxKernel(const unsigned int* in, unsigned int* out, unsigned int n) {
	unsigned int maxValue = 0;

	const unsigned int totalNumThreads = gridDim.x * blockDim.x;
	for (unsigned int i = blockIdx.x * blockDim.x + threadIdx.x; i < n; i += totalNumThreads) {
		maxValue = max(maxValue, in[i]);
	}
	maxValue = blockReduceMax(maxValue);
	if (threadIdx.x == 0) {
		out[blockIdx.x] = maxValue; // one value for each block
	}
}

__global__
void
reduceMaxKernel(unsigned int* blockVal, unsigned int numBlocks) {
	unsigned int maxValue = 0;

	const unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx != 0) return; // single thread

	for (int i = 0; i < numBlocks; ++i) {
		maxValue = max(maxValue, blockVal[i]);
	}

	blockVal[0] = maxValue;
}

} // namespace Lab
