/***************************************************************************
 *  Copyright 2014, 2017, 2018, 2020, 2025 Marcelo Y. Matuda               *
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

#include "cuda.h"

#include "CUDACoherenceFactor.cuh"



namespace Lab {

// Must have the same values as in the .cpp file.
#define TRANSP_BLOCK_SIZE 16
#define NUM_RX_ELEM 32

// NVIDIA sm_50 or newer:
//   - Shared memory has 32 banks of 32 bits.
template<typename TFloat>
__global__
void
transposeKernel(
		TFloat* rawData,
		TFloat* rawDataT,
		unsigned int oldSizeX,
		unsigned int oldSizeY)
{
	__shared__ TFloat temp[TRANSP_BLOCK_SIZE][TRANSP_BLOCK_SIZE + 1]; // +1 to avoid bank conflicts

	unsigned int iX = blockIdx.x * TRANSP_BLOCK_SIZE + threadIdx.x;
	unsigned int iY = blockIdx.y * TRANSP_BLOCK_SIZE + threadIdx.y;

	if (iX < oldSizeX && iY < oldSizeY) {
		temp[threadIdx.x][threadIdx.y] = rawData[iY * oldSizeX + iX];
	}

	__syncthreads();

	iX = blockIdx.y * TRANSP_BLOCK_SIZE + threadIdx.x;
	iY = blockIdx.x * TRANSP_BLOCK_SIZE + threadIdx.y;

	if (iX < oldSizeY && iY < oldSizeX) {
		rawDataT[iY * oldSizeY + iX] = temp[threadIdx.y][threadIdx.x];
	}
}

template<typename TFloat>
__global__
void
processImageKernel(
		TFloat* rawData,
		unsigned int numGridPoints,
		TFloat* gridValueRe,
		TFloat* gridValueIm,
		TFloat* rxApod)
{
	TFloat rxSignalListRe[NUM_RX_ELEM];
	TFloat rxSignalListIm[NUM_RX_ELEM];

	const unsigned int point = blockIdx.x * blockDim.x + threadIdx.x;
	if (point >= numGridPoints) return;

	unsigned int idx1 = point;
	unsigned int idx2 = point + numGridPoints;
	const unsigned int step = 2 * numGridPoints;
	for (int i = 0; i < NUM_RX_ELEM; ++i, idx1 += step, idx2 += step) {
		rxSignalListRe[i] = rawData[idx1];
		rxSignalListIm[i] = rawData[idx2];
	}

	TFloat sumRe = 0;
	TFloat sumIm = 0;
	for (int i = 0; i < NUM_RX_ELEM; ++i) {
		sumRe += rxSignalListRe[i] * rxApod[i];
		sumIm += rxSignalListIm[i] * rxApod[i];
	}

	gridValueRe[point] += sumRe;
	gridValueIm[point] += sumIm;
}

template<typename TFloat>
__global__
void
processImagePCFKernel(
		TFloat* rawData,
		unsigned int numGridPoints,
		TFloat* gridValueRe,
		TFloat* gridValueIm,
		TFloat* rxApod,
		TFloat pcfFactor)
{
	TFloat rxSignalListRe[NUM_RX_ELEM];
	TFloat rxSignalListIm[NUM_RX_ELEM];

	const unsigned int point = blockIdx.x * blockDim.x + threadIdx.x;
	if (point >= numGridPoints) return;

	unsigned int idx1 = point;
	unsigned int idx2 = point + numGridPoints;
	const unsigned int step = 2 * numGridPoints;
	for (int i = 0; i < NUM_RX_ELEM; ++i, idx1 += step, idx2 += step) {
		rxSignalListRe[i] = rawData[idx1];
		rxSignalListIm[i] = rawData[idx2];
	}

	const TFloat pcf = calcPCF<TFloat, NUM_RX_ELEM>(rxSignalListRe, rxSignalListIm, pcfFactor);

	TFloat sumRe = 0.0;
	TFloat sumIm = 0.0;
	for (int i = 0; i < NUM_RX_ELEM; ++i) {
		sumRe += rxSignalListRe[i] * rxApod[i];
		sumIm += rxSignalListIm[i] * rxApod[i];
	}

	gridValueRe[point] += sumRe * pcf;
	gridValueIm[point] += sumIm * pcf;
}

void
execTransposeKernel(const dim3& gridDim, const dim3& blockDim,
		float* rawData, float* rawDataT, unsigned int oldSizeX, unsigned int oldSizeY)
{
	transposeKernel<<<gridDim, blockDim>>>(rawData, rawDataT, oldSizeX, oldSizeY);
}

void
execProcessImageKernel(const dim3& gridDim, const dim3& blockDim,
		float* rawData, unsigned int numGridPoints, float* gridValueRe,
		float* gridValueIm, float* rxApod)
{
	processImageKernel<<<gridDim, blockDim>>>(rawData, numGridPoints, gridValueRe, gridValueIm, rxApod);
}

void
execProcessImagePCFKernel(const dim3& gridDim, const dim3& blockDim,
		float* rawData, unsigned int numGridPoints, float* gridValueRe,
		float* gridValueIm, float* rxApod, float pcfFactor)
{
	processImagePCFKernel<<<gridDim, blockDim>>>(rawData, numGridPoints, gridValueRe, gridValueIm, rxApod, pcfFactor);
}

} // namespace Lab
