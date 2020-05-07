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

#include "NumericArrayOfRectangularSourcesCUDAImpulseResponse.h"

#include <cmath> /* ceil */
#include <limits>

#include "Exception.h"
#include "Log.h"
#include "Util.h" /* pi */

#include "CUDAReduce.cuh"
#include "CUDAUtil.h"
#include "NumericRectangularSourceCUDAImpulseResponse.cuh"

#ifndef MFloat
# define MFloat float
#endif

#define REDUCE_BLOCK_SIZE 1024
#define INITIAL_H_SIZE 10000



namespace Lab {

template<typename TFloat>
__global__
void
numericArraySourceIRKernel(
		unsigned int numElem,
		unsigned int numSubElem,
		TFloat x,
		TFloat y,
		TFloat z,
		TFloat k1,
		TFloat k2,
		const unsigned int* activeElem,
		const TFloat* elemDelay,
		const TFloat* elemPosX,
		const TFloat* elemPosY,
		const TFloat* subElemX,
		const TFloat* subElemY,
		unsigned int* n0,
		TFloat* value)
{
	const unsigned int subElemIdx = blockIdx.x * blockDim.x + threadIdx.x;
	if (subElemIdx >= numSubElem) {
		// Get data from the first sub-element, to help min/max(n0).
		const unsigned int activeElemIdx = activeElem[0];
		const TFloat dx = x - subElemX[0] - elemPosX[activeElemIdx];
		const TFloat dy = y - subElemY[0] - elemPosY[activeElemIdx];
		const TFloat r = sqrt(dx * dx + dy * dy + z * z);
		const unsigned int idx = (numElem - 1U) * numSubElem + subElemIdx;
		n0[idx] = rint(r * k1 + elemDelay[0]);
		return;
	}

	for (int i = 0; i < numElem; ++i) {
		const unsigned int activeElemIdx = activeElem[i];
		const TFloat dx = x - subElemX[subElemIdx] - elemPosX[activeElemIdx];
		const TFloat dy = y - subElemY[subElemIdx] - elemPosY[activeElemIdx];
		const TFloat r = sqrt(dx * dx + dy * dy + z * z);
		const unsigned int idx = i * numSubElem + subElemIdx;
		n0[idx] = rint(r * k1 + elemDelay[i]);
		value[idx] = k2 / r;
	}
}

//=============================================================================

struct NumericArrayOfRectangularSourcesCUDAImpulseResponse::CUDAData {
	CUDAContext c; // must be the first member
	CUDAHostDevMem<MFloat> subElemX;
	CUDAHostDevMem<MFloat> subElemY;
	CUDAHostDevMem<MFloat> elemDelay;
	CUDAHostDevMem<MFloat> elemPosX;
	CUDAHostDevMem<MFloat> elemPosY;
	CUDAHostDevMem<unsigned int> activeElem;
	CUDADevMem<unsigned int> n0;
	CUDADevMem<MFloat> value;
	CUDAHostDevMem<unsigned int> minN0;
	CUDAHostDevMem<unsigned int> maxN0;
	CUDAHostDevMem<MFloat> h;
};

NumericArrayOfRectangularSourcesCUDAImpulseResponse::NumericArrayOfRectangularSourcesCUDAImpulseResponse(
		MFloat samplingFreq,
		MFloat propagationSpeed,
		MFloat sourceWidth,
		MFloat sourceHeight,
		MFloat subElemSize,
		const std::vector<XY<MFloat>>& elemPos,
		const std::vector<MFloat>& focusDelay /* s */)
			: samplingFreq_(samplingFreq)
			, propagationSpeed_(propagationSpeed)
			, subElemWidth_()
			, subElemHeight_()
			, numElem_()
			, numSubElem_()
{
	if (elemPos.empty()) {
		THROW_EXCEPTION(InvalidParameterException, "Empty element position list.");
	}
	if (focusDelay.empty()) {
		THROW_EXCEPTION(InvalidParameterException, "Empty element delay list.");
	}
	if (focusDelay.size() > elemPos.size()) {
		THROW_EXCEPTION(InvalidParameterException, "Size of focusDelay is greater than size of elemPos.");
	}
	numElem_ = focusDelay.size(); // may be less than elemPos.size()

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

	if (Log::isDebugEnabled()) {
		int device;
		exec(cudaGetDevice(&device));

		cudaDeviceProp prop;
		exec(cudaGetDeviceProperties(&prop, device));
		LOG_DEBUG << "CUDA device: " << prop.name;
	}

	data_->subElemX   = CUDAHostDevMem<MFloat>(numSubElem_);
	data_->subElemY   = CUDAHostDevMem<MFloat>(numSubElem_);
	data_->elemDelay  = CUDAHostDevMem<MFloat>(numElem_);
	data_->elemPosX   = CUDAHostDevMem<MFloat>(elemPos.size());
	data_->elemPosY   = CUDAHostDevMem<MFloat>(elemPos.size());
	data_->activeElem = CUDAHostDevMem<unsigned int>(numElem_);
	const unsigned int numReduceThreads = CUDAUtil::roundUpToMultipleOfBlockSize(numElem_ * numSubElem_, REDUCE_BLOCK_SIZE);
	data_->n0         = CUDADevMem<unsigned int>(numReduceThreads);
	data_->value      = CUDADevMem<MFloat>(numElem_ * numSubElem_);
	const unsigned int numReduceBlocks = numReduceThreads / REDUCE_BLOCK_SIZE;
	data_->minN0      = CUDAHostDevMem<unsigned int>(numReduceBlocks);
	data_->maxN0      = CUDAHostDevMem<unsigned int>(numReduceBlocks);
	data_->h          = CUDAHostDevMem<MFloat>(INITIAL_H_SIZE);

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

	for (unsigned int i = 0; i < numElem_; ++i) {
		data_->elemDelay.hostPtr[i] = focusDelay[i] * samplingFreq_;
	}
	exec(data_->elemDelay.copyHostToDevice());

	for (unsigned int i = 0; i < elemPos.size(); ++i) {
		data_->elemPosX.hostPtr[i] = elemPos[i].x;
		data_->elemPosY.hostPtr[i] = elemPos[i].y;
	}
	exec(data_->elemPosX.copyHostToDevice());
	exec(data_->elemPosY.copyHostToDevice());

	LOG_DEBUG << "[NumericArrayOfRectangularSourcesCUDAImpulseResponse] numElemX=" << numElemX << " numElemY=" << numElemY;
}

// This destructor can't be inline.
NumericArrayOfRectangularSourcesCUDAImpulseResponse::~NumericArrayOfRectangularSourcesCUDAImpulseResponse()
{
}

void
NumericArrayOfRectangularSourcesCUDAImpulseResponse::getImpulseResponse(
							MFloat x,
							MFloat y,
							MFloat z,
							std::size_t& hOffset,
							std::vector<MFloat>& h,
							std::vector<unsigned int>* activeElemList)
{
	if (activeElemList) {
		if (activeElemList->empty()) {
			THROW_EXCEPTION(InvalidParameterException, "Empty active element list.");
		} else {
			if (activeElemList->size() != numElem_) {
				THROW_EXCEPTION(InvalidParameterException, "Active element size is not equal to number of delays.");
			}
			for (unsigned int i = 0; i < numElem_; ++i) {
				data_->activeElem.hostPtr[i] = (*activeElemList)[i];
			}
		}
	} else {
		for (unsigned int i = 0; i < numElem_; ++i) {
			data_->activeElem.hostPtr[i] = i;
		}
	}
	exec(data_->activeElem.copyHostToDevice());

	{
		const unsigned int blockSize = REDUCE_BLOCK_SIZE;
		const unsigned int numBlocks = CUDAUtil::numberOfBlocks(numSubElem_, blockSize);

		numericArraySourceIRKernel<MFloat><<<numBlocks, blockSize>>>(
			numElem_,
			numSubElem_,
			x, y, z,
			samplingFreq_ / propagationSpeed_,
			samplingFreq_ * subElemWidth_ * subElemHeight_ / (MFloat(2.0 * pi) * propagationSpeed_),
			data_->activeElem.devPtr,
			data_->elemDelay.devPtr,
			data_->elemPosX.devPtr,
			data_->elemPosY.devPtr,
			data_->subElemX.devPtr,
			data_->subElemY.devPtr,
			data_->n0.devPtr,
			data_->value.devPtr);
		checkKernelLaunchError();
	}
	{
		const unsigned int blockSize = REDUCE_BLOCK_SIZE;
		const unsigned int numBlocks = CUDAUtil::numberOfBlocks(numElem_ * numSubElem_, blockSize);

		blockReduceMinMaxKernel<<<numBlocks, blockSize>>>(
			data_->n0.devPtr,
			data_->minN0.devPtr,
			data_->maxN0.devPtr,
			numElem_ * numSubElem_);
		checkKernelLaunchError();

		reduceMinMaxKernel<<<1, 1>>>(
			data_->minN0.devPtr,
			data_->maxN0.devPtr,
			numBlocks);
		checkKernelLaunchError();
	}
	exec(data_->minN0.copyDeviceToHost(sizeof(unsigned int)));
	const unsigned int minN0 = *(data_->minN0.hostPtr);
	exec(data_->maxN0.copyDeviceToHost(sizeof(unsigned int)));
	const unsigned int maxN0 = *(data_->maxN0.hostPtr);

	const unsigned int hSize = maxN0 - minN0 + 1U;
	if (hSize > data_->h.sizeInBytes / sizeof(MFloat)) {
		data_->h = CUDAHostDevMem<MFloat>(hSize * 2);
	}
	exec(cudaMemset(data_->h.devPtr, 0, hSize * sizeof(MFloat)));
	{
		const unsigned int blockSize = 64;
		const unsigned int numBlocks = CUDAUtil::numberOfBlocks(numElem_ * numSubElem_, blockSize);

		accumulateIRSamplesKernel<MFloat><<<numBlocks, blockSize>>>(
			numElem_ * numSubElem_,
			minN0,
			data_->n0.devPtr,
			data_->value.devPtr,
			data_->h.devPtr);
		checkKernelLaunchError();
	}
	exec(data_->h.copyDeviceToHost(hSize * sizeof(MFloat)));

	h.resize(hSize);
	for (unsigned int i = 0; i < hSize; ++i) {
		h[i] = data_->h.hostPtr[i];
	}
	hOffset = minN0;

	//LOG_DEBUG << "[NumericArrayOfRectangularSourcesCUDAImpulseResponse] minN0=" << minN0 << " maxN0=" << maxN0;
}

} // namespace Lab
