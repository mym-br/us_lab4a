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

#include "NumericRectangularSourceCUDAImpulseResponse.h"

#include <cmath> /* ceil */
#include <limits>

#include "Log.h"
#include "Util.h" /* pi */

#include "CUDAUtil.h"

#ifndef MFloat
# define MFloat float
#endif

#define REDUCE_BLOCK_SIZE 1024
#define INITIAL_H_SIZE 10000000



namespace Lab {

__forceinline__
__device__
unsigned int
warpReduceMin(unsigned int value) {
	for (unsigned int offset = warpSize / 2; offset > 0; offset /= 2)  {
		// Get the value from thread n + offset, but only if it is in the same warp.
		// Otherwise return the current value.
		const unsigned int shiftedValue = __shfl_down_sync(0xffffffff, value, offset);
		value = min(value, shiftedValue);
	}
	// Here the first thread in the warp contains the minimum value.
	// The other threads in the warp contain garbage.
	return value;
}

// blockDim.x must be multiple of warpSize.
__forceinline__
__device__
unsigned int
blockReduceMin(unsigned int value) {
	__shared__ unsigned int shared[32];

	unsigned int warpThreadIdx = threadIdx.x % warpSize; // index of the thread in the warp
	unsigned int warpIdx = threadIdx.x / warpSize; // index of the warp in the block

	value = warpReduceMin(value);

	if (warpThreadIdx == 0) {
		shared[warpIdx] = value; // only the first thread in the warp has a valid value
	}

	__syncthreads();

	// The first warp gets the values in the shared memory.
	// Thread 0 gets the value from warp 0 in the block.
	// Thread 1 gets the value from warp 1 in the block.
	// ...
	const unsigned int numWarpsInBlock = blockDim.x / warpSize;
	value = (threadIdx.x < numWarpsInBlock) ? shared[threadIdx.x] : 0;

	if (warpIdx == 0) {
		value = warpReduceMin(value);
	}
	// Here the first thread in the block has the final value.
	// The other threads in the warp contain garbage.

	return value;
}

__global__
void
blockReduceMinKernel(const unsigned int* in, unsigned int* out, int n) {
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
reduceMinKernel(unsigned int* blockVal, int numBlocks) {
	unsigned int minValue = 0xffffffff;

	const unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx != 0) return; // single thread

	for (unsigned int i = 0; i < numBlocks; ++i) {
		minValue = min(minValue, blockVal[i]);
	}

	blockVal[0] = minValue;
}

//-----------------------------------------------------------------------------

__forceinline__
__device__
unsigned int
warpReduceMax(unsigned int value) {
	for (unsigned int offset = warpSize / 2; offset > 0; offset /= 2)  {
		// Get the value from thread n + offset, but only if it is in the same warp.
		// Otherwise return the current value.
		const unsigned int shiftedValue = __shfl_down_sync(0xffffffff, value, offset);
		value = max(value, shiftedValue);
	}
	// Here the first thread in the warp contains the maximum value.
	// The other threads in the warp contain garbage.
	return value;
}

// blockDim.x must be multiple of warpSize.
__forceinline__
__device__
unsigned int
blockReduceMax(unsigned int value) {
	__shared__ unsigned int shared[32];

	unsigned int warpThreadIdx = threadIdx.x % warpSize; // index of the thread in the warp
	unsigned int warpIdx = threadIdx.x / warpSize; // index of the warp in the block

	value = warpReduceMax(value);

	if (warpThreadIdx == 0) {
		shared[warpIdx] = value; // only the first thread in the warp has a valid value
	}

	__syncthreads();

	// The first warp gets the values in the shared memory.
	// Thread 0 gets the value from warp 0 in the block.
	// Thread 1 gets the value from warp 1 in the block.
	// ...
	const unsigned int numWarpsInBlock = blockDim.x / warpSize;
	value = (threadIdx.x < numWarpsInBlock) ? shared[threadIdx.x] : 0;

	if (warpIdx == 0) {
		value = warpReduceMax(value);
	}
	// Here the first thread in the block has the final value.
	// The other threads in the warp contain garbage.

	return value;
}

__global__
void
blockReduceMaxKernel(const unsigned int* in, unsigned int* out, int n) {
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
reduceMaxKernel(unsigned int* blockVal, int numBlocks) {
	unsigned int maxValue = 0;

	const unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx != 0) return; // single thread

	for (unsigned int i = 0; i < numBlocks; ++i) {
		maxValue = max(maxValue, blockVal[i]);
	}

	blockVal[0] = maxValue;
}

//-----------------------------------------------------------------------------

__global__
void
numericRectangularSourceIRKernel(
		unsigned int numSubElem,
		float x,
		float y,
		float z,
		float k1,
		float k2,
		float* subElemX,
		float* subElemY,
		unsigned int* n0,
		float* value)
{
	const unsigned int subElemIdx = blockIdx.x * blockDim.x + threadIdx.x;
	if (subElemIdx >= numSubElem) {
		n0[subElemIdx] = 0xffffffff;
		return;
	}

	const float dy = y - subElemY[subElemIdx];
	const float dx = x - subElemX[subElemIdx];
	const float r = sqrtf(dx * dx + dy * dy + z * z);
	n0[subElemIdx] = rintf(r * k1);
	value[subElemIdx] = k2 / r;
}

__global__
void
padKernel(
		unsigned int n1,
		unsigned int n2,
		unsigned int* data,
		unsigned int value)
{
	const unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx >= n1 && idx < n2) {
		data[idx] = value;
	}
}

__global__
void
accumulateIRSamplesKernel(
		unsigned int numSubElem,
		unsigned int minN0,
		const unsigned int* n0,
		const float* value,
		float* h)
{
	const unsigned int subElemIdx = blockIdx.x * blockDim.x + threadIdx.x;
	if (subElemIdx >= numSubElem) {
		return;
	}

	// Different sub-elements may have the same value for n0.
	atomicAdd(h + n0[subElemIdx] - minN0, value[subElemIdx]);
}

//=============================================================================

struct NumericRectangularSourceCUDAImpulseResponse::CUDAData {
	CUDAHostDevMem<MFloat> subElemX;
	CUDAHostDevMem<MFloat> subElemY;
	CUDADevMem<unsigned int> n0;
	CUDAHostDevMem<MFloat> value;
	CUDAHostDevMem<unsigned int> minN0;
	CUDAHostDevMem<unsigned int> maxN0;
	CUDAHostDevMem<MFloat> h;
};

NumericRectangularSourceCUDAImpulseResponse::NumericRectangularSourceCUDAImpulseResponse(
		MFloat samplingFreq,
		MFloat propagationSpeed,
		MFloat sourceWidth,
		MFloat sourceHeight,
		MFloat subElemSize)
			: samplingFreq_(samplingFreq)
			, propagationSpeed_(propagationSpeed)
			, subElemWidth_()
			, subElemHeight_()
{
	unsigned int numElemX, numElemY;
	if (subElemSize > sourceWidth) {
		numElemX = 1;
	} else {
		numElemX = static_cast<unsigned int>(std::ceil(sourceWidth / subElemSize));
	}
	if (subElemSize > sourceHeight) {
		numElemY = 1;
	} else {
		numElemY = static_cast<unsigned int>(std::ceil(sourceHeight / subElemSize));
	}
	const std::size_t n = numElemY * numElemX;
	if (n >= std::numeric_limits<unsigned int>::max()) {
		THROW_EXCEPTION(InvalidValueException, "Too many sub-elements in the rectangular area.");
	}
	numSubElem_ = n;
	subElemWidth_  = sourceWidth  / numElemX;
	subElemHeight_ = sourceHeight / numElemY;

	data_ = std::make_unique<CUDAData>();
	data_->subElemX = CUDAHostDevMem<MFloat>(numSubElem_);
	data_->subElemY = CUDAHostDevMem<MFloat>(numSubElem_);
	const unsigned int numReduceThreads = CUDAUtil::roundUpToMultipleOfBlockSize(numSubElem_, REDUCE_BLOCK_SIZE);
	data_->n0       = CUDADevMem<unsigned int>(numReduceThreads);
	data_->value    = CUDAHostDevMem<MFloat>(numSubElem_);
	const unsigned int numReduceBlocks = numReduceThreads / REDUCE_BLOCK_SIZE;
	data_->minN0    = CUDAHostDevMem<unsigned int>(numReduceBlocks);
	data_->maxN0    = CUDAHostDevMem<unsigned int>(numReduceBlocks);
	data_->h        = CUDAHostDevMem<MFloat>(INITIAL_H_SIZE);

	const MFloat halfW = 0.5 * (numElemX - 1);
	const MFloat halfH = 0.5 * (numElemY - 1);
	for (unsigned int iy = 0; iy < numElemY; ++iy) {
		for (unsigned int ix = 0; ix < numElemX; ++ix) {
			unsigned int idx = iy * numElemX + ix;
			data_->subElemX.hostPtr[idx] = (ix - halfW) * subElemWidth_;
			data_->subElemY.hostPtr[idx] = (iy - halfH) * subElemHeight_;
		}
	}
	exec(data_->subElemX.copyHostToDevice());
	exec(data_->subElemY.copyHostToDevice());

	LOG_DEBUG << "[NumericRectangularSourceCUDAImpulseResponse] numElemX=" << numElemX << " numElemY=" << numElemY;
}

// This destructor can't be inline.
NumericRectangularSourceCUDAImpulseResponse::~NumericRectangularSourceCUDAImpulseResponse()
{
}

void
NumericRectangularSourceCUDAImpulseResponse::getImpulseResponse(
							MFloat x,
							MFloat y,
							MFloat z,
							std::size_t& hOffset,
							std::vector<MFloat>& h)
{
	{
		const unsigned int blockSize = REDUCE_BLOCK_SIZE;
		const unsigned int numBlocks = CUDAUtil::numberOfBlocks(numSubElem_, blockSize);

		numericRectangularSourceIRKernel<<<numBlocks, blockSize>>>(
			 numSubElem_,
			 x, y, z,
			 samplingFreq_ / propagationSpeed_,
			 samplingFreq_ * subElemWidth_ * subElemHeight_ / (MFloat(2.0 * pi) * propagationSpeed_),
			 data_->subElemX.devPtr,
			 data_->subElemY.devPtr,
			 data_->n0.devPtr,
			 data_->value.devPtr);
		checkKernelLaunchError();
	}
	{
		const unsigned int blockSize = REDUCE_BLOCK_SIZE;
		const unsigned int numBlocks = CUDAUtil::numberOfBlocks(numSubElem_, blockSize);

		blockReduceMinKernel<<<numBlocks, blockSize>>>(
			data_->n0.devPtr,
			data_->minN0.devPtr,
			numSubElem_);
		checkKernelLaunchError();

		reduceMinKernel<<<1, 1>>>(
			data_->minN0.devPtr,
			numBlocks);
		checkKernelLaunchError();
	}
	exec(data_->minN0.copyDeviceToHost(sizeof(unsigned int)));
	const unsigned int minN0 = *(data_->minN0.hostPtr);

	{
		const unsigned int blockSize = REDUCE_BLOCK_SIZE;
		const unsigned int numBlocks = CUDAUtil::numberOfBlocks(numSubElem_, blockSize);

		// Prepare for reduceMaxKernel.
		padKernel<<<numBlocks, blockSize>>>(
			numSubElem_,
			numBlocks * blockSize,
			data_->n0.devPtr,
			0);
		checkKernelLaunchError();

		blockReduceMaxKernel<<<numBlocks, blockSize>>>(
			data_->n0.devPtr,
			data_->maxN0.devPtr,
			numSubElem_);
		checkKernelLaunchError();

		reduceMaxKernel<<<1, 1>>>(
			data_->maxN0.devPtr,
			numBlocks);
		checkKernelLaunchError();
	}
	exec(data_->maxN0.copyDeviceToHost(sizeof(unsigned int)));
	const unsigned int maxN0 = *(data_->maxN0.hostPtr);

	const unsigned int hSize = maxN0 - minN0 + 1U;
	if (hSize > data_->h.sizeInBytes / sizeof(MFloat)) {
		data_->h = CUDAHostDevMem<MFloat>(hSize * 2);
	}
	exec(cudaMemset(data_->h.devPtr, 0, hSize * sizeof(MFloat)));
	{
		const unsigned int blockSize = 64;
		const unsigned int numBlocks = CUDAUtil::numberOfBlocks(numSubElem_, blockSize);

		accumulateIRSamplesKernel<<<numBlocks, blockSize>>>(
			numSubElem_,
			minN0,
			data_->n0.devPtr,
			data_->value.devPtr,
			data_->h.devPtr);
	}
	exec(data_->h.copyDeviceToHost(hSize * sizeof(MFloat)));

	h.resize(hSize);
	for (unsigned int i = 0; i < hSize; ++i) {
		h[i] = data_->h.hostPtr[i];
	}
	hOffset = minN0;

	//LOG_DEBUG << "[NumericRectangularSourceCUDAImpulseResponse] minN0=" << minN0 << " maxN0=" << maxN0;
}

} // namespace Lab
