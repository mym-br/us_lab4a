/***************************************************************************
 *  Copyright 2020, 2025 Marcelo Y. Matuda                                 *
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
#include "CUDAGeometry.cuh"



namespace Lab {

#define NUM_RX_ELEM 32

template<typename TFloat>
__global__
void
calculateDelaysSTAKernel(
		unsigned int numCols,
		unsigned int numRows,
		TFloat invCT,
		const TFloat* xArray,
		const TFloat (*gridXZ)[2],
		TFloat* delayTensor)
{
	const unsigned int row = blockIdx.x * blockDim.x + threadIdx.x;
	if (row >= numRows) return;

	const unsigned int col = blockIdx.y * blockDim.y + threadIdx.y;
	if (col >= numCols) return;

	const unsigned int elem = blockIdx.z * blockDim.z + threadIdx.z;
	if (elem >= NUM_RX_ELEM) return;

	const TFloat (*point)[2] = gridXZ + col * numRows + row;
	delayTensor[((elem * numCols) + col) * numRows + row] =
			distance2DY0<TFloat>(xArray[elem], (*point)[0], (*point)[1]) * invCT;
}

template<typename TFloat>
__global__
void
processRowColumnSTAKernel(
		unsigned int numCols,
		unsigned int numRows,
		TFloat signalOffset,
		const TFloat (*signalTensor)[2],
		unsigned int signalTensorN2,
		unsigned int signalTensorN3,
		unsigned int firstTxElem,
		unsigned int lastTxElem,
		const TFloat* rxApod,
		const TFloat* delayTensor,
		TFloat* gridValue)
{
	const unsigned int signalLength = signalTensorN3;
	const unsigned int maxPosition = signalLength - 2;

	const unsigned int row = blockIdx.x * blockDim.x + threadIdx.x;
	if (row >= numRows) return;

	const unsigned int col = blockIdx.y * blockDim.y + threadIdx.y;
	if (col >= numCols) return;

	TFloat sumRe = 0;
	TFloat sumIm = 0;
	for (int txElem = firstTxElem; txElem <= lastTxElem; ++txElem) {
		const TFloat txDelay = delayTensor[(txElem * numCols + col) * numRows + row];
		const TFloat txOffset = signalOffset + txDelay;
		for (int rxElem = 0; rxElem < NUM_RX_ELEM; ++rxElem) {
			const TFloat (*p)[2] = signalTensor + (((txElem - firstTxElem) * signalTensorN2) + rxElem) * signalLength;
			// Linear interpolation.
			const TFloat position = txOffset + delayTensor[(rxElem * numCols + col) * numRows + row];
			if (position >= 0.0f) {
				const unsigned int positionIdx = static_cast<unsigned int>(position);
				if (positionIdx <= maxPosition) {
					const TFloat k = position - positionIdx;
					TFloat v0[2]; // complex
					v0[0] = p[positionIdx][0];
					v0[1] = p[positionIdx][1];
					TFloat v1[2]; // complex
					v1[0] = p[positionIdx + 1][0];
					v1[1] = p[positionIdx + 1][1];
					TFloat v[2]; // complex
					v[0] = v0[0] + k * (v1[0] - v0[0]);
					v[1] = v0[1] + k * (v1[1] - v0[1]);

					sumRe += v[0] * rxApod[rxElem];
					sumIm += v[1] * rxApod[rxElem];
				}
			}
		}
	}
	const unsigned int point = col * numRows + row;
	gridValue[point] = sqrt(sumRe * sumRe + sumIm * sumIm);
}

template<typename TFloat>
__global__
void
processRowColumnSTAPCFKernel(
		unsigned int numCols,
		unsigned int numRows,
		TFloat signalOffset,
		const TFloat (*signalTensor)[2],
		unsigned int signalTensorN2,
		unsigned int signalTensorN3,
		unsigned int firstTxElem,
		unsigned int lastTxElem,
		const TFloat* rxApod,
		const TFloat* delayTensor,
		TFloat pcfFactor,
		TFloat* gridValue,
		TFloat* gridFactor)
{
	TFloat rxSignalListRe[NUM_RX_ELEM];
	TFloat rxSignalListIm[NUM_RX_ELEM];

	const unsigned int signalLength = signalTensorN3;
	const unsigned int maxPosition = signalLength - 2;

	const unsigned int row = blockIdx.x * blockDim.x + threadIdx.x;
	if (row >= numRows) return;

	const unsigned int col = blockIdx.y * blockDim.y + threadIdx.y;
	if (col >= numCols) return;

	for (int i = 0; i < NUM_RX_ELEM; ++i) {
		rxSignalListRe[i] = 0;
		rxSignalListIm[i] = 0;
	}
	for (int txElem = firstTxElem; txElem <= lastTxElem; ++txElem) {
		const TFloat txDelay = delayTensor[(txElem * numCols + col) * numRows + row];
		const TFloat txOffset = signalOffset + txDelay;
		for (int rxElem = 0; rxElem < NUM_RX_ELEM; ++rxElem) {
			const TFloat (*p)[2] = signalTensor + (((txElem - firstTxElem) * signalTensorN2) + rxElem) * signalLength;
			// Linear interpolation.
			const TFloat position = txOffset + delayTensor[(rxElem * numCols + col) * numRows + row];
			if (position >= 0.0f) {
				const unsigned int positionIdx = static_cast<unsigned int>(position);
				if (positionIdx <= maxPosition) {
					const TFloat k = position - positionIdx;
					TFloat v0[2]; // complex
					v0[0] = p[positionIdx][0];
					v0[1] = p[positionIdx][1];
					TFloat v1[2]; // complex
					v1[0] = p[positionIdx + 1][0];
					v1[1] = p[positionIdx + 1][1];
					TFloat v[2]; // complex
					v[0] = v0[0] + k * (v1[0] - v0[0]);
					v[1] = v0[1] + k * (v1[1] - v0[1]);

					rxSignalListRe[rxElem] += v[0] * rxApod[rxElem];
					rxSignalListIm[rxElem] += v[1] * rxApod[rxElem];
				}
			}
		}
	}

	const TFloat pcf = calcPCF<TFloat, NUM_RX_ELEM>(rxSignalListRe, rxSignalListIm, pcfFactor);
	TFloat sumRe = 0;
	TFloat sumIm = 0;
	for (int i = 0; i < NUM_RX_ELEM; ++i) {
		sumRe += rxSignalListRe[i];
		sumIm += rxSignalListIm[i];
	}

	const unsigned int point = col * numRows + row;
	gridValue[point] = sqrt(sumRe * sumRe + sumIm * sumIm);
	gridFactor[point] = pcf;
}

void execCalculateDelaysSTAKernel(const dim3& gridDim, const dim3& blockDim,
		unsigned int numCols, unsigned int numRows, float invCT,
		const float* xArray, const float (*gridXZ)[2], float* delayTensor)
{
	calculateDelaysSTAKernel<<<gridDim, blockDim>>>(numCols, numRows, invCT, xArray, gridXZ, delayTensor);
}

void execProcessRowColumnSTAKernel(const dim3& gridDim, const dim3& blockDim,
		unsigned int numCols, unsigned int numRows, float signalOffset,
		const float (*signalTensor)[2], unsigned int signalTensorN2, unsigned int signalTensorN3,
		unsigned int firstTxElem, unsigned int lastTxElem, const float* rxApod,
		const float* delayTensor, float* gridValue)
{
	processRowColumnSTAKernel<<<gridDim, blockDim>>>(numCols, numRows, signalOffset,
		signalTensor, signalTensorN2, signalTensorN3, firstTxElem, lastTxElem,
		rxApod, delayTensor, gridValue);
}

void execProcessRowColumnSTAPCFKernel(const dim3& gridDim, const dim3& blockDim,
		unsigned int numCols, unsigned int numRows, float signalOffset,
		const float (*signalTensor)[2], unsigned int signalTensorN2, unsigned int signalTensorN3,
		unsigned int firstTxElem, unsigned int lastTxElem, const float* rxApod, const float* delayTensor,
		float pcfFactor, float* gridValue, float* gridFactor)
{
	processRowColumnSTAPCFKernel<<<gridDim, blockDim>>>(numCols, numRows, signalOffset,
		signalTensor, signalTensorN2, signalTensorN3, firstTxElem, lastTxElem,
		rxApod, delayTensor, pcfFactor, gridValue, gridFactor);
}

} // namespace Lab
