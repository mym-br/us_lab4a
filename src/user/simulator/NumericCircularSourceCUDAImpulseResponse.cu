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

#include "NumericCircularSourceCUDAImpulseResponse.h"

#include <cmath> /* abs, floor, sqrt */

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

#define NUMERIC_CIRCULAR_SOURCE_CUDA_IMPULSE_RESPONSE_USE_RANDOM 1

#ifdef NUMERIC_CIRCULAR_SOURCE_CUDA_IMPULSE_RESPONSE_USE_RANDOM
# include <random>
#endif



namespace Lab {

struct NumericCircularSourceCUDAImpulseResponse::CUDAData {
	CUDAHostDevMem<MFloat> subElemX;
	CUDAHostDevMem<MFloat> subElemY;
	CUDADevMem<unsigned int> n0;
	CUDADevMem<MFloat> value;
	CUDAHostDevMem<unsigned int> minN0;
	CUDAHostDevMem<unsigned int> maxN0;
	CUDAHostDevMem<MFloat> h;
};

NumericCircularSourceCUDAImpulseResponse::NumericCircularSourceCUDAImpulseResponse(
		MFloat samplingFreq,
		MFloat propagationSpeed,
		MFloat sourceRadius,
		MFloat numSubElemInRadius)
			: samplingFreq_(samplingFreq)
			, propagationSpeed_(propagationSpeed)
			, subElemArea_()
			, numSubElem_()
{
	std::vector<MFloat> subElemX;
	std::vector<MFloat> subElemY;

#ifdef NUMERIC_CIRCULAR_SOURCE_CUDA_IMPULSE_RESPONSE_USE_RANDOM
	const MFloat area = pi * (sourceRadius * sourceRadius);
	const MFloat subElemDensity = numSubElemInRadius * numSubElemInRadius / (sourceRadius * sourceRadius);
	numSubElem_ = static_cast<unsigned int>(subElemDensity * area);
	subElemArea_ = area / numSubElem_;

	std::mt19937 rndGen;
	std::random_device rd;
	rndGen.seed(rd());
	std::uniform_real_distribution<MFloat> dist(0.0, 1.0);

	// http://extremelearning.com.au/a-simple-method-to-construct-isotropic-quasirandom-blue-noise-point-sequences/
	// Infinite R2 jittered sequence.
	// i = 0, 1, 2, ...
	// u0, u1 = random numbers [0.0, 1.0)
	auto jitteredPoint2D = [&](int i, double u0, double u1, MFloat& x, MFloat& y) {
		constexpr double lambda = 0.5; // jitter parameter ( > 0)
		constexpr double phi = 1.324717957244746;
		constexpr double alpha0 = 1.0 / phi;
		constexpr double alpha1 = 1.0 / (phi * phi);
		constexpr double delta0 = 0.76;
		constexpr double i0 = 0.300;
		const double k = lambda * delta0 * std::sqrt(pi) / (4.0 * std::sqrt(i + i0));
		const double i1 = i + 1;
		const double x0 = alpha0 * i1 + k * u0; // x0 > 0
		const double y0 = alpha1 * i1 + k * u1; // y0 > 0
		x = x0 - std::floor(x0);
		y = y0 - std::floor(y0);
	};

	int i = 0;
	while (subElemX.size() < numSubElem_) {
		MFloat x, y;
		jitteredPoint2D(i, dist(rndGen), dist(rndGen), x, y);
		// [0.0, 1.0) --> [-1.0, 1.0)
		x = -1.0 + 2.0 * x;
		y = -1.0 + 2.0 * y;
		x *= sourceRadius;
		y *= sourceRadius;

		if (std::sqrt(x * x + y * y) < sourceRadius) {
			subElemX.push_back(x);
			subElemY.push_back(y);
		}
		++i;
	}
#else
	const MFloat d = sourceRadius / numSubElemInRadius; // sub-element side
	subElemArea_ = d * d * 2; // multiplied by 2 because only one half of the circle is used
	const MFloat d2 = d * 0.5;

	const unsigned int n = numSubElemInRadius;
	for (unsigned int iy = 0; iy < n; ++iy) {
		const MFloat yc = d2 + iy * d;
		const MFloat yt = yc + d2; // to test if the sub-element is inside the circle
		for (unsigned int ix = 0; ix < n; ++ix) {
			const MFloat xc = d2 + ix * d;
			const MFloat xt = xc + d2;
			if (std::sqrt(xt * xt + yt * yt) <= sourceRadius) {
				subElemX.push_back(xc);
				subElemY.push_back(yc);
				subElemX.push_back(-xc);
				subElemY.push_back(yc);
			} else {
				break;
			}
		}
	}
	numSubElem_ = subElemX.size();
#endif

	data_ = std::make_unique<CUDAData>();
	data_->subElemX = CUDAHostDevMem<MFloat>(numSubElem_);
	data_->subElemY = CUDAHostDevMem<MFloat>(numSubElem_);
	const unsigned int numReduceThreads = CUDAUtil::roundUpToMultipleOfBlockSize(numSubElem_, REDUCE_BLOCK_SIZE);
	data_->n0       = CUDADevMem<unsigned int>(numReduceThreads);
	data_->value    = CUDADevMem<MFloat>(numSubElem_);
	const unsigned int numReduceBlocks = numReduceThreads / REDUCE_BLOCK_SIZE;
	data_->minN0    = CUDAHostDevMem<unsigned int>(numReduceBlocks);
	data_->maxN0    = CUDAHostDevMem<unsigned int>(numReduceBlocks);
	data_->h        = CUDAHostDevMem<MFloat>(INITIAL_H_SIZE);

	for (unsigned int i = 0; i < numSubElem_; ++i) {
		data_->subElemX.hostPtr[i] = subElemX[i];
		data_->subElemY.hostPtr[i] = subElemY[i];
	}
	exec(data_->subElemX.copyHostToDevice());
	exec(data_->subElemY.copyHostToDevice());

	LOG_DEBUG << "[NumericCircularSourceCUDAImpulseResponse] numSubElem=" << numSubElem_;
}

// This destructor can't be inline.
NumericCircularSourceCUDAImpulseResponse::~NumericCircularSourceCUDAImpulseResponse()
{
}

void
NumericCircularSourceCUDAImpulseResponse::getImpulseResponse(
							MFloat x,
							MFloat y,
							MFloat z,
							std::size_t& hOffset,
							std::vector<MFloat>& h)
{
	// The field is symmetric.
	x = std::sqrt(x * x + y * y);
	y = 0;
	z = std::abs(z);

	{
		const unsigned int blockSize = REDUCE_BLOCK_SIZE;
		const unsigned int numBlocks = CUDAUtil::numberOfBlocks(numSubElem_, blockSize);

		numericSourceIRKernel<MFloat><<<numBlocks, blockSize>>>(
			 numSubElem_,
			 x, y, z,
			 samplingFreq_ / propagationSpeed_,
			 samplingFreq_ * subElemArea_ / (MFloat(2.0 * pi) * propagationSpeed_),
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

	//LOG_DEBUG << "[NumericCircularSourceCUDAImpulseResponse] minN0=" << minN0 << " maxN0=" << maxN0;
}

} // namespace Lab
