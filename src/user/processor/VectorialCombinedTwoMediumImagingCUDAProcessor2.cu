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
#include "CUDAGeometry.cuh"



namespace Lab {

// NVIDIA sm_50 or newer:
//   - Shared memory has 32 banks of 32 bits.

// Must have the same value as in the .cpp file.
#define NUM_RX_ELEM 32

template<typename TFloat>
__device__
TFloat
calcTwoMediumTravelTime(TFloat x1, TFloat z1, TFloat xi, TFloat zi, TFloat x2, TFloat z2, TFloat invC1, TFloat invC2)
{
	const TFloat dx1 = xi - x1;
	const TFloat dz1 = zi - z1;
	const TFloat dx2 = x2 - xi;
	const TFloat dz2 = z2 - zi;
	return sqrt(dx1 * dx1 + dz1 * dz1) * invC1 + sqrt(dx2 * dx2 + dz2 * dz2) * invC2;
}

template<typename TFloat>
__device__
void
findMinTimeInTwoSteps(
		unsigned int blockSize,
		TFloat c1, TFloat c2,
		const TFloat (*interfacePointList)[2],
		unsigned int interfacePointListSize,
		TFloat x1, TFloat z1, TFloat x2, TFloat z2,
		TFloat& tMin, unsigned int& idxMin)
{
	const TFloat invC1 = 1 / c1;
	const TFloat invC2 = 1 / c2;

	// First step: step = blockSize
	{
		const TFloat (*point)[2] = interfacePointList;
		const TFloat t = calcTwoMediumTravelTime<TFloat>(x1, z1, (*point)[0], (*point)[1], x2, z2, invC1, invC2);
		tMin = t;
		idxMin = 0;
	}
	for (unsigned int i = 1; i < interfacePointListSize; i += blockSize) {
		const TFloat (*point)[2] = interfacePointList + i;
		const TFloat t = calcTwoMediumTravelTime<TFloat>(x1, z1, (*point)[0], (*point)[1], x2, z2, invC1, invC2);
		if (t < tMin) {
			tMin = t;
			idxMin = i;
		}
	}

	// Second step: step = 1
	const unsigned int iBegin = (idxMin > (blockSize - 1U)) ? idxMin - (blockSize - 1U) : 0;
	const unsigned int iEnd = umin(idxMin + blockSize, interfacePointListSize);
	for (unsigned int i = iBegin; i < iEnd; ++i) {
		const TFloat (*point)[2] = interfacePointList + i;
		const TFloat t = calcTwoMediumTravelTime<TFloat>(x1, z1, (*point)[0], (*point)[1], x2, z2, invC1, invC2);
		if (t < tMin) {
			tMin = t;
			idxMin = i;
		}
	}
}

// This function is not efficient in GPUs, but it avoids the transfer of the large delayTensor.
template<typename TFloat>
__global__
void
calculateDelaysTwoMediumKernel(
		unsigned int numCols,
		unsigned int numRows,
		unsigned int numElementsMux,
		TFloat fs,
		TFloat fsInvC2,
		TFloat c1,
		TFloat c2,
		unsigned int fermatBlockSize,
		const TFloat (*interfacePointList)[2],
		unsigned int interfacePointListSize,
		const TFloat* xArray,
		const unsigned int* minRowIdx,
		const TFloat* medium1DelayMatrix,
		const TFloat (*gridXZ)[2],
		TFloat* delayTensor)
{
	const unsigned int col = blockIdx.x * blockDim.x + threadIdx.x;
	if (col >= numCols) return;

	const unsigned int elem = blockIdx.y * blockDim.y + threadIdx.y;
	if (elem >= numElementsMux) return;

	if (minRowIdx[col] >= numRows) return;

	unsigned int lastInterfaceIdx = 0;

	// The first row above the interface.
	{
		const TFloat (*point)[2] = gridXZ + col * numRows + minRowIdx[col];

		// Fermat's principle. Find the fastest path.
		TFloat tMin;
		unsigned int idxMin;
		findMinTimeInTwoSteps<TFloat>(
				fermatBlockSize,
				c1, c2,
				interfacePointList,
				interfacePointListSize,
				xArray[elem], 0, (*point)[0], (*point)[1],
				tMin, idxMin);
		delayTensor[((elem * numCols) + col) * numRows + minRowIdx[col]] = tMin * fs;
		lastInterfaceIdx = idxMin;
	}

	const TFloat* medium1Delays = medium1DelayMatrix + elem * interfacePointListSize;

	for (unsigned int row = minRowIdx[col] + 1U; row < numRows; ++row) {
		const TFloat (*point)[2] = gridXZ + col * numRows + row;
		unsigned int idxMin = lastInterfaceIdx;
		TFloat tC2Min;
		{
			const TFloat (*ifPoint)[2] = interfacePointList + idxMin;
			tC2Min = medium1Delays[idxMin] + distance2D<TFloat>((*ifPoint)[0], (*ifPoint)[1], (*point)[0], (*point)[1]);
		}
		for (unsigned int idxSearch = idxMin + 1U; idxSearch < interfacePointListSize; ++idxSearch) {
			const TFloat (*ifPoint)[2] = interfacePointList + idxSearch;
			const TFloat tC2 = medium1Delays[idxSearch] + distance2D<TFloat>((*ifPoint)[0], (*ifPoint)[1], (*point)[0], (*point)[1]);
			if (tC2 >= tC2Min) {
				break;
			} else {
				tC2Min = tC2;
				idxMin = idxSearch;
			}
		}
		if (idxMin == lastInterfaceIdx) { // if the previous search was not successful
			for (int idxSearch = static_cast<int>(idxMin) - 1; idxSearch >= 0; --idxSearch) { // if idxMin = 0, idxSearch will start with -1
				const TFloat (*ifPoint)[2] = interfacePointList + idxSearch;
				const TFloat tC2 = medium1Delays[idxSearch] + distance2D<TFloat>((*ifPoint)[0], (*ifPoint)[1], (*point)[0], (*point)[1]);
				if (tC2 >= tC2Min) {
					break;
				} else {
					tC2Min = tC2;
					idxMin = idxSearch;
				}
			}
		}

		delayTensor[((elem * numCols) + col) * numRows + row] = tC2Min * fsInvC2;
		lastInterfaceIdx = idxMin;
	}
}

template<typename TFloat>
__global__
void
processRowColumnWithOneTxElemKernel(
		unsigned int numCols,
		unsigned int numRows,
		TFloat signalOffset,
		const TFloat (*signalTensor)[2],
		unsigned int signalTensorN2,
		unsigned int signalTensorN3,
		unsigned int baseElem,
		unsigned int baseElemIdx,
		unsigned int txElem,
		const unsigned int* minRowIdx,
		const TFloat* delayTensor,
		TFloat (*gridValue)[2],
		const TFloat* rxApod)
{
	const unsigned int signalLength = signalTensorN3;
	const unsigned int maxPosition = signalLength - 2;

	const unsigned int col = blockIdx.y * blockDim.y + threadIdx.y;
	if (col >= numCols) return;

	const unsigned int row = blockIdx.x * blockDim.x + threadIdx.x;
	if (row < minRowIdx[col] || row >= numRows) return;

	const TFloat txDelay = delayTensor[((baseElem + txElem) * numCols + col) * numRows + row];
	const TFloat txOffset = signalOffset + txDelay;

	TFloat rxSignalSumRe = 0;
	TFloat rxSignalSumIm = 0;
	for (unsigned int rxElem = 0; rxElem < NUM_RX_ELEM; ++rxElem) {
		const TFloat (*p)[2] = signalTensor + baseElemIdx * signalTensorN2 * signalTensorN3
					+ rxElem * signalLength;
		// Linear interpolation.
		const TFloat position = txOffset + delayTensor[((baseElem + rxElem) * numCols + col) * numRows + row];
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

				rxSignalSumRe += v[0] * rxApod[rxElem];
				rxSignalSumIm += v[1] * rxApod[rxElem];
			}
		}
	}
	const unsigned int point = col * numRows + row;
	gridValue[point][0] += rxSignalSumRe;
	gridValue[point][1] += rxSignalSumIm;
}

template<typename TFloat>
__global__
void
processRowColumnWithOneTxElemPCFKernel(
		unsigned int numCols,
		unsigned int numRows,
		TFloat signalOffset,
		const TFloat (*signalTensor)[2],
		unsigned int signalTensorN2,
		unsigned int signalTensorN3,
		unsigned int baseElem,
		unsigned int baseElemIdx,
		unsigned int txElem,
		const unsigned int* minRowIdx,
		const TFloat* delayTensor,
		TFloat (*gridValue)[2],
		const TFloat* rxApod,
		TFloat pcfFactor)
{
	TFloat rxSignalListRe[NUM_RX_ELEM];
	TFloat rxSignalListIm[NUM_RX_ELEM];

	const unsigned int signalLength = signalTensorN3;
	const unsigned int maxPosition = signalLength - 2;

	const unsigned int col = blockIdx.y * blockDim.y + threadIdx.y;
	if (col >= numCols) return;

	const unsigned int row = blockIdx.x * blockDim.x + threadIdx.x;
	if (row < minRowIdx[col] || row >= numRows) return;

	const TFloat txDelay = delayTensor[((baseElem + txElem) * numCols + col) * numRows + row];
	const TFloat txOffset = signalOffset + txDelay;

	for (unsigned int rxElem = 0; rxElem < NUM_RX_ELEM; ++rxElem) {
		const TFloat (*p)[2] = signalTensor + baseElemIdx * signalTensorN2 * signalTensorN3
					+ rxElem * signalLength;
		// Linear interpolation.
		const TFloat position = txOffset + delayTensor[((baseElem + rxElem) * numCols + col) * numRows + row];
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

				rxSignalListRe[rxElem] = v[0];
				rxSignalListIm[rxElem] = v[1];
			}
		}
	}

	const TFloat pcf = calcPCF<TFloat, NUM_RX_ELEM>(rxSignalListRe, rxSignalListIm, pcfFactor);
	TFloat sumRe = 0;
	TFloat sumIm = 0;
	for (int i = 0; i < NUM_RX_ELEM; ++i) {
		sumRe += rxSignalListRe[i] * rxApod[i];
		sumIm += rxSignalListIm[i] * rxApod[i];
	}

	const unsigned int point = col * numRows + row;
	gridValue[point][0] += sumRe * pcf;
	gridValue[point][1] += sumIm * pcf;
}

void
execCalculateDelaysTwoMediumKernel(const dim3& gridDim, const dim3& blockDim,
		unsigned int numCols, unsigned int numRows, unsigned int numElementsMux, float fs,
		float fsInvC2, float c1, float c2, unsigned int fermatBlockSize,
		const float (*interfacePointList)[2], unsigned int interfacePointListSize,
		const float* xArray, const unsigned int* minRowIdx,
		const float* medium1DelayMatrix, const float (*gridXZ)[2], float* delayTensor)
{
	calculateDelaysTwoMediumKernel<<<gridDim, blockDim>>>(numCols, numRows, numElementsMux,
		fs, fsInvC2, c1, c2, fermatBlockSize, interfacePointList, interfacePointListSize,
		xArray, minRowIdx, medium1DelayMatrix, gridXZ, delayTensor);
}

void
execProcessRowColumnWithOneTxElemKernel(const dim3& gridDim, const dim3& blockDim,
		unsigned int numCols, unsigned int numRows, float signalOffset, const float (*signalTensor)[2],
		unsigned int signalTensorN2, unsigned int signalTensorN3, unsigned int baseElem,
		unsigned int baseElemIdx, unsigned int txElem, const unsigned int* minRowIdx,
		const float* delayTensor, float (*gridValue)[2], const float* rxApod)
{
	processRowColumnWithOneTxElemKernel<<<gridDim, blockDim>>>(numCols, numRows, signalOffset, signalTensor,
		signalTensorN2, signalTensorN3, baseElem, baseElemIdx, txElem,
		minRowIdx, delayTensor, gridValue, rxApod);
}

void
execProcessRowColumnWithOneTxElemPCFKernel(const dim3& gridDim, const dim3& blockDim,
		unsigned int numCols, unsigned int numRows, float signalOffset, const float (*signalTensor)[2],
		unsigned int signalTensorN2, unsigned int signalTensorN3, unsigned int baseElem,
		unsigned int baseElemIdx, unsigned int txElem, const unsigned int* minRowIdx,
		const float* delayTensor, float (*gridValue)[2], const float* rxApod, float pcfFactor)
{
	processRowColumnWithOneTxElemPCFKernel<<<gridDim, blockDim>>>(numCols, numRows, signalOffset, signalTensor,
		signalTensorN2, signalTensorN3, baseElem, baseElemIdx, txElem, minRowIdx,
		delayTensor, gridValue, rxApod, pcfFactor);
}

} // namespace Lab
