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

struct NumericRectangularSourceCUDAImpulseResponse::CUDAData {
	CUDAContext c; // must be the first member
	CUDAHostDevMem<MFloat> subElemX;
	CUDAHostDevMem<MFloat> subElemY;
	CUDADevMem<unsigned int> n0;
	CUDADevMem<MFloat> value;
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
			, numSubElem_()
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

	if (Log::isDebugEnabled()) {
		int device;
		exec(cudaGetDevice(&device));

		cudaDeviceProp prop;
		exec(cudaGetDeviceProperties(&prop, device));
		LOG_DEBUG << "CUDA device: " << prop.name;
	}

	data_->subElemX = CUDAHostDevMem<MFloat>(numSubElem_);
	data_->subElemY = CUDAHostDevMem<MFloat>(numSubElem_);
	const unsigned int numReduceThreads = CUDAUtil::roundUpToMultipleOfBlockSize(numSubElem_, REDUCE_BLOCK_SIZE);
	data_->n0       = CUDADevMem<unsigned int>(numReduceThreads);
	data_->value    = CUDADevMem<MFloat>(numSubElem_);
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

		numericSourceIRKernel<MFloat><<<numBlocks, blockSize>>>(
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

		blockReduceMinMaxKernel<<<numBlocks, blockSize>>>(
			data_->n0.devPtr,
			data_->minN0.devPtr,
			data_->maxN0.devPtr,
			numSubElem_);
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
		const unsigned int numBlocks = CUDAUtil::numberOfBlocks(numSubElem_, blockSize);

		accumulateIRSamplesKernel<MFloat><<<numBlocks, blockSize>>>(
			numSubElem_,
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

	//LOG_DEBUG << "[NumericRectangularSourceCUDAImpulseResponse] minN0=" << minN0 << " maxN0=" << maxN0;
}

} // namespace Lab
