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

#ifndef CUDA_REDUCE_CUH
#define CUDA_REDUCE_CUH

namespace Lab {

// Call with <<<numBlocks, 1024>>>.
extern __global__ void blockReduceMinMaxKernel(const unsigned int* in, unsigned int* outMin, unsigned int* outMax, int n);
// Call with <<<1, 1>>>.
extern __global__ void reduceMinMaxKernel(unsigned int* blockMinVal, unsigned int* blockMaxVal, int numBlocks);

// Call with <<<numBlocks, 1024>>>.
extern __global__ void blockReduceMinKernel(const unsigned int* in, unsigned int* out, int n);
// Call with <<<1, 1>>>.
extern __global__ void reduceMinKernel(unsigned int* blockVal, int numBlocks);

// Call with <<<numBlocks, 1024>>>.
extern __global__ void blockReduceMaxKernel(const unsigned int* in, unsigned int* out, int n);
// Call with <<<1, 1>>>.
extern __global__ void reduceMaxKernel(unsigned int* blockVal, int numBlocks);



__forceinline__
__device__
void
warpReduceMinMax(unsigned int& minValue, unsigned int& maxValue) {
	for (unsigned int offset = warpSize / 2; offset > 0; offset /= 2)  {
		// Get the value from thread n + offset, but only if it is in the same warp.
		// Otherwise return the current value.
		const unsigned int shiftedMinValue = __shfl_down_sync(0xffffffff, minValue, offset);
		minValue = min(minValue, shiftedMinValue);

		const unsigned int shiftedMaxValue = __shfl_down_sync(0xffffffff, maxValue, offset);
		maxValue = max(maxValue, shiftedMaxValue);
	}
	// Here the first thread in the warp contains the minimum and maximum value.
	// The other threads in the warp contain garbage.
}

// blockDim.x must be multiple of warpSize.
__forceinline__
__device__
void
blockReduceMinMax(unsigned int& minValue, unsigned int& maxValue) {
	__shared__ unsigned int sharedMin[32];
	__shared__ unsigned int sharedMax[32];

	unsigned int warpThreadIdx = threadIdx.x % warpSize; // index of the thread in the warp
	unsigned int warpIdx = threadIdx.x / warpSize; // index of the warp in the block

	warpReduceMinMax(minValue, maxValue);

	if (warpThreadIdx == 0) {
		sharedMin[warpIdx] = minValue; // only the first thread in the warp has a valid value
		sharedMax[warpIdx] = maxValue;
	}

	__syncthreads();

	// The first warp gets the values in the shared memory.
	// Thread 0 gets the value from warp 0 in the block.
	// Thread 1 gets the value from warp 1 in the block.
	// ...
	const unsigned int numWarpsInBlock = blockDim.x / warpSize;
	minValue = (threadIdx.x < numWarpsInBlock) ? sharedMin[threadIdx.x] : 0;
	maxValue = (threadIdx.x < numWarpsInBlock) ? sharedMax[threadIdx.x] : 0;

	if (warpIdx == 0) {
		warpReduceMinMax(minValue, maxValue);
	}
	// Here the first thread in the block has the final value.
	// The other threads in the warp contain garbage.
}

//-----------------------------------------------------------------------------

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

} // namespace Lab

#endif // CUDA_REDUCE_CUH
