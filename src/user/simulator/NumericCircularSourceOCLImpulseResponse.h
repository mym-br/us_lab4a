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
#ifndef NUMERIC_CIRCULAR_SOURCE_OCL_IMPULSE_RESPONSE_H
#define NUMERIC_CIRCULAR_SOURCE_OCL_IMPULSE_RESPONSE_H

#include <cmath> /* abs, floor, sqrt */
#include <cstddef> /* std::size_t */
#include <memory>
#include <sstream>
#include <string>
#include <type_traits> /* is_same */
#include <vector>

#include "Exception.h"
#include "Log.h"
#include "Util.h" /* pi */

#include <CL/cl2.hpp>
#include "NumericRectangularSourceOCLImpulseResponse.h"
#include "OCLAtomic.h"
#include "OCLReduce.h"
#include "OCLUtil.h"

#define NUMERIC_CIRCULAR_SOURCE_OCL_IMPULSE_RESPONSE_ACCUMULATE_IN_GPU 1
#define NUMERIC_CIRCULAR_SOURCE_OCL_IMPULSE_RESPONSE_ACCUMULATE_WITH_LOCAL_MEM 1
#define NUMERIC_CIRCULAR_SOURCE_OCL_IMPULSE_RESPONSE_USE_RANDOM 1

#ifdef NUMERIC_CIRCULAR_SOURCE_OCL_IMPULSE_RESPONSE_USE_RANDOM
# include <random>
#endif



namespace Lab {

template<typename TFloat>
class NumericCircularSourceOCLImpulseResponse {
public:
	NumericCircularSourceOCLImpulseResponse(
					TFloat samplingFreq,
					TFloat propagationSpeed,
					TFloat sourceRadius,
					TFloat numSubElemInRadius);
	~NumericCircularSourceOCLImpulseResponse() = default;

	// Return h/c.
	void getImpulseResponse(TFloat x, TFloat y, TFloat z,
				std::size_t& hOffset /* samples */, std::vector<TFloat>& h);

private:
	enum {
		REDUCE_GROUP_SIZE = 1024,
		INITIAL_H_SIZE = 10000
	};

	NumericCircularSourceOCLImpulseResponse(const NumericCircularSourceOCLImpulseResponse&) = delete;
	NumericCircularSourceOCLImpulseResponse& operator=(const NumericCircularSourceOCLImpulseResponse&) = delete;
	NumericCircularSourceOCLImpulseResponse(NumericCircularSourceOCLImpulseResponse&&) = delete;
	NumericCircularSourceOCLImpulseResponse& operator=(NumericCircularSourceOCLImpulseResponse&&) = delete;

	TFloat samplingFreq_;
	TFloat propagationSpeed_;
	TFloat subElemArea_;
	unsigned int numSubElem_;

	// OpenCL.
	cl::Context clContext_;
	cl::Program clProgram_;
	cl::CommandQueue clCommandQueue_;
	cl::Buffer subElemXCLBuffer_;
	cl::Buffer subElemYCLBuffer_;
	cl::Buffer n0CLBuffer_;
	cl::Buffer valueCLBuffer_;
	std::unique_ptr<OCLPinnedHostMem<unsigned int>> minN0HostMem_;
	cl::Buffer minN0CLBuffer_;
	std::unique_ptr<OCLPinnedHostMem<unsigned int>> maxN0HostMem_;
	cl::Buffer maxN0CLBuffer_;
#ifdef NUMERIC_CIRCULAR_SOURCE_OCL_IMPULSE_RESPONSE_ACCUMULATE_IN_GPU
	std::unique_ptr<OCLPinnedHostMem<TFloat>> hHostMem_;
	cl::Buffer hCLBuffer_;
#else
	std::unique_ptr<OCLPinnedHostMem<unsigned int>> n0HostMem_;
	std::unique_ptr<OCLPinnedHostMem<TFloat>> valueHostMem_;
#endif
};



template<typename TFloat>
NumericCircularSourceOCLImpulseResponse<TFloat>::NumericCircularSourceOCLImpulseResponse(
		TFloat samplingFreq,
		TFloat propagationSpeed,
		TFloat sourceRadius,
		TFloat numSubElemInRadius)
			: samplingFreq_(samplingFreq)
			, propagationSpeed_(propagationSpeed)
			, subElemArea_()
			, numSubElem_()
{
	std::vector<TFloat> subElemX;
	std::vector<TFloat> subElemY;

#ifdef NUMERIC_CIRCULAR_SOURCE_OCL_IMPULSE_RESPONSE_USE_RANDOM
	const TFloat area = pi * (sourceRadius * sourceRadius);
	const TFloat subElemDensity = numSubElemInRadius * numSubElemInRadius / (sourceRadius * sourceRadius);
	numSubElem_ = static_cast<unsigned int>(subElemDensity * area);
	subElemArea_ = area / numSubElem_;

	std::mt19937 rndGen;
	std::random_device rd;
	rndGen.seed(rd());
	std::uniform_real_distribution<TFloat> dist(0.0, 1.0);

	// http://extremelearning.com.au/a-simple-method-to-construct-isotropic-quasirandom-blue-noise-point-sequences/
	// Infinite R2 jittered sequence.
	// i = 0, 1, 2, ...
	// u0, u1 = random numbers [0.0, 1.0)
	auto jitteredPoint2D = [&](int i, double u0, double u1, TFloat& x, TFloat& y) {
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
		TFloat x, y;
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
	const TFloat d = sourceRadius / numSubElemInRadius; // sub-element side
	subElemArea_ = d * d * 2; // multiplied by 2 because only one half of the circle is used
	const TFloat d2 = d * 0.5;

	const unsigned int n = numSubElemInRadius;
	for (unsigned int iy = 0; iy < n; ++iy) {
		const TFloat yc = d2 + iy * d;
		const TFloat yt = yc + d2; // to test if the sub-element is inside the circle
		for (unsigned int ix = 0; ix < n; ++ix) {
			const TFloat xc = d2 + ix * d;
			const TFloat xt = xc + d2;
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
	clContext_ = OCLUtil::initOpenCL();

	std::vector<std::string> kernelStrings = {
					OCLAtomic::code() +
					OCLReduce::code() +
					NumericRectangularSourceOCLImpulseResponse<TFloat>::getKernels()};
	clProgram_ = cl::Program(clContext_, kernelStrings);
	std::ostringstream progOpt;
	progOpt << OCL_PROGRAM_BUILD_OPTIONS
			<< " -DMFloat=" << TemplateUtil::typeName<TFloat>();
#ifdef NUMERIC_CIRCULAR_SOURCE_OCL_IMPULSE_RESPONSE_ACCUMULATE_WITH_LOCAL_MEM
	progOpt << " -DWITH_LOCAL_MEM=1";
#endif
	if (std::is_same<TFloat, float>::value) {
		progOpt << " -DSINGLE_PREC=1";
	}
	try {
		clProgram_.build(progOpt.str().c_str());
	} catch (...) {
		std::ostringstream msg;
		msg << "Error during OpenCL kernel compilation:\n";
		auto buildInfo = clProgram_.getBuildInfo<CL_PROGRAM_BUILD_LOG>();
		for (auto& pair : buildInfo) {
			msg << pair.second << "\n";
		}
		THROW_EXCEPTION(OCLException, msg.str());
	}

	clCommandQueue_ = cl::CommandQueue(clContext_);

	subElemXCLBuffer_ = cl::Buffer(clContext_, CL_MEM_READ_ONLY , sizeof(TFloat) * numSubElem_);
	subElemYCLBuffer_ = cl::Buffer(clContext_, CL_MEM_READ_ONLY , sizeof(TFloat) * numSubElem_);
	const unsigned int numReduceThreads = OCLUtil::roundUpToMultipleOfGroupSize(numSubElem_, REDUCE_GROUP_SIZE);
	n0CLBuffer_       = cl::Buffer(clContext_, CL_MEM_READ_WRITE, sizeof(unsigned int) * numReduceThreads);
	valueCLBuffer_    = cl::Buffer(clContext_, CL_MEM_READ_WRITE, sizeof(TFloat) * numSubElem_);
	const unsigned int numReduceGroups = numReduceThreads / REDUCE_GROUP_SIZE;
	minN0HostMem_     = std::make_unique<OCLPinnedHostMem<unsigned int>>(clContext_, clCommandQueue_,
									numReduceGroups,
									CL_MAP_READ);
	minN0CLBuffer_    = cl::Buffer(clContext_, CL_MEM_READ_WRITE, sizeof(unsigned int) * numReduceGroups);
	maxN0HostMem_     = std::make_unique<OCLPinnedHostMem<unsigned int>>(clContext_, clCommandQueue_,
									numReduceGroups,
									CL_MAP_READ);
	maxN0CLBuffer_    = cl::Buffer(clContext_, CL_MEM_READ_WRITE, sizeof(unsigned int) * numReduceGroups);
#ifdef NUMERIC_CIRCULAR_SOURCE_OCL_IMPULSE_RESPONSE_ACCUMULATE_IN_GPU
	hHostMem_         = std::make_unique<OCLPinnedHostMem<TFloat>>(clContext_, clCommandQueue_,
									INITIAL_H_SIZE,
									CL_MAP_READ);
	hCLBuffer_        = cl::Buffer(clContext_, CL_MEM_READ_WRITE, sizeof(TFloat) * INITIAL_H_SIZE);
#else
	n0HostMem_        = std::make_unique<OCLPinnedHostMem<unsigned int>>(clContext_, clCommandQueue_,
									numReduceThreads,
									CL_MAP_READ);
	valueHostMem_     = std::make_unique<OCLPinnedHostMem<TFloat>>(clContext_, clCommandQueue_,
									numSubElem_,
									CL_MAP_READ);
#endif
	clCommandQueue_.enqueueWriteBuffer(
		subElemXCLBuffer_, CL_BLOCKING, 0 /* offset */,
		subElemX.size() * sizeof(TFloat), subElemX.data());
	clCommandQueue_.enqueueWriteBuffer(
		subElemYCLBuffer_, CL_BLOCKING, 0 /* offset */,
		subElemY.size() * sizeof(TFloat), subElemY.data());

	LOG_DEBUG << "[NumericCircularSourceOCLImpulseResponse] numSubElem=" << numSubElem_;
}

template<typename TFloat>
void
NumericCircularSourceOCLImpulseResponse<TFloat>::getImpulseResponse(
							TFloat x,
							TFloat y,
							TFloat z,
							std::size_t& hOffset,
							std::vector<TFloat>& h)
{
	// The field is symmetric.
	x = std::sqrt(x * x + y * y);
	y = 0;
	z = std::abs(z);

	try {
		cl::Kernel kernel(clProgram_, "numericSourceIRKernel");
		kernel.setArg(0, numSubElem_);
		kernel.setArg(1, x);
		kernel.setArg(2, y);
		kernel.setArg(3, z);
		kernel.setArg(4, samplingFreq_ / propagationSpeed_);
		kernel.setArg(5, samplingFreq_ * subElemArea_ / (TFloat(2.0 * pi) * propagationSpeed_));
		kernel.setArg(6, subElemXCLBuffer_);
		kernel.setArg(7, subElemYCLBuffer_);
		kernel.setArg(8, n0CLBuffer_);
		kernel.setArg(9, valueCLBuffer_);

		const unsigned int groupSize = REDUCE_GROUP_SIZE;
		const unsigned int globalN0 = OCLUtil::roundUpToMultipleOfGroupSize(numSubElem_, groupSize);

		clCommandQueue_.enqueueNDRangeKernel(
			kernel,
			cl::NullRange, // offset
			cl::NDRange(globalN0), // global
			cl::NDRange(groupSize), // local
			nullptr /* previous events */, nullptr);
	} catch (cl::Error& e) {
		THROW_EXCEPTION(OCLException, "[numericSourceIRKernel] OpenCL error: " << e.what() << " (" << e.err() << ").");
	}

	const unsigned int reduceGlobalN0 = OCLUtil::roundUpToMultipleOfGroupSize(numSubElem_, REDUCE_GROUP_SIZE);
	try {
		cl::Kernel kernel(clProgram_, "groupReduceMinMaxKernel");
		kernel.setArg(0, n0CLBuffer_);
		kernel.setArg(1, minN0CLBuffer_);
		kernel.setArg(2, maxN0CLBuffer_);
		kernel.setArg(3, numSubElem_);
		kernel.setArg(4, cl::Local(REDUCE_GROUP_SIZE * sizeof(unsigned int)));
		kernel.setArg(5, cl::Local(REDUCE_GROUP_SIZE * sizeof(unsigned int)));

		clCommandQueue_.enqueueNDRangeKernel(
			kernel,
			cl::NullRange, // offset
			cl::NDRange(reduceGlobalN0), // global
			cl::NDRange(REDUCE_GROUP_SIZE), // local
			nullptr /* previous events */, nullptr);
	} catch (cl::Error& e) {
		THROW_EXCEPTION(OCLException, "[groupReduceMinMaxKernel] OpenCL error: " << e.what() << " (" << e.err() << ").");
	}

	try {
		cl::Kernel kernel(clProgram_, "reduceMinMaxKernel");
		kernel.setArg(0, minN0CLBuffer_);
		kernel.setArg(1, maxN0CLBuffer_);
		kernel.setArg(2, reduceGlobalN0 / REDUCE_GROUP_SIZE);

		clCommandQueue_.enqueueNDRangeKernel(
			kernel,
			cl::NullRange, // offset
			cl::NDRange(1), // global
			cl::NDRange(1), // local
			nullptr /* previous events */, nullptr);
	} catch (cl::Error& e) {
		THROW_EXCEPTION(OCLException, "[reduceMinMaxKernel] OpenCL error: " << e.what() << " (" << e.err() << ").");
	}

	try {
		clCommandQueue_.enqueueReadBuffer(
			minN0CLBuffer_, CL_BLOCKING, 0 /* offset */,
			sizeof(unsigned int), minN0HostMem_->hostPtr);
		clCommandQueue_.enqueueReadBuffer(
			maxN0CLBuffer_, CL_BLOCKING, 0 /* offset */,
			sizeof(unsigned int), maxN0HostMem_->hostPtr);
	} catch (cl::Error& e) {
		THROW_EXCEPTION(OCLException, "[Read minN0, maxN0] OpenCL error: " << e.what() << " (" << e.err() << ").");
	}

	const unsigned int minN0 = *(minN0HostMem_->hostPtr);
	const unsigned int maxN0 = *(maxN0HostMem_->hostPtr);
	const unsigned int hSize = maxN0 - minN0 + 1U;

#ifdef NUMERIC_CIRCULAR_SOURCE_OCL_IMPULSE_RESPONSE_ACCUMULATE_IN_GPU
	if (hSize > hHostMem_->sizeInBytes / sizeof(TFloat)) {
		hHostMem_  = std::make_unique<OCLPinnedHostMem<TFloat>>(clContext_, clCommandQueue_,
									hSize * 2,
									CL_MAP_READ);
		hCLBuffer_ = cl::Buffer(clContext_, CL_MEM_READ_WRITE, sizeof(TFloat) * hSize * 2);
	}

	clCommandQueue_.enqueueFillBuffer(
		hCLBuffer_, (TFloat) 0, 0 /* offset */, hSize * sizeof(TFloat));

	try {
		cl::Kernel kernel(clProgram_, "accumulateIRSamplesKernel");
		kernel.setArg(0, numSubElem_);
		kernel.setArg(1, minN0);
		kernel.setArg(2, n0CLBuffer_);
		kernel.setArg(3, valueCLBuffer_);
		kernel.setArg(4, hCLBuffer_);
#ifdef NUMERIC_CIRCULAR_SOURCE_OCL_IMPULSE_RESPONSE_ACCUMULATE_WITH_LOCAL_MEM
		kernel.setArg(5, hSize);
		kernel.setArg(6, cl::Local(hSize * sizeof(TFloat)));

		const unsigned int groupSize = 1024;
#else
		const unsigned int groupSize = 64;
#endif
		const unsigned int globalN0 = OCLUtil::roundUpToMultipleOfGroupSize(numSubElem_, groupSize);

		clCommandQueue_.enqueueNDRangeKernel(
			kernel,
			cl::NullRange, // offset
			cl::NDRange(globalN0), // global
			cl::NDRange(groupSize), // local
			nullptr /* previous events */, nullptr);
	} catch (cl::Error& e) {
		THROW_EXCEPTION(OCLException, "[accumulateIRSamplesKernel] OpenCL error: " << e.what() << " (" << e.err() << ").");
	}

	try {
		clCommandQueue_.enqueueReadBuffer(
			hCLBuffer_, CL_BLOCKING, 0 /* offset */,
			hSize * sizeof(TFloat), hHostMem_->hostPtr);
	} catch (cl::Error& e) {
		THROW_EXCEPTION(OCLException, "[Read h] OpenCL error: " << e.what() << " (" << e.err() << ").");
	}

	h.resize(hSize);
	for (unsigned int i = 0; i < hSize; ++i) {
		h[i] = hHostMem_->hostPtr[i];
	}
#else
	clCommandQueue_.enqueueReadBuffer(
		n0CLBuffer_, CL_BLOCKING, 0 /* offset */,
		numSubElem_ * sizeof(unsigned int), n0HostMem_->hostPtr);
	clCommandQueue_.enqueueReadBuffer(
		valueCLBuffer_, CL_BLOCKING, 0 /* offset */,
		numSubElem_ * sizeof(TFloat), valueHostMem_->hostPtr);

	h.assign(hSize, 0);
	for (unsigned int i = 0; i < numSubElem_; ++i) {
		h[n0HostMem_->hostPtr[i] - minN0] += valueHostMem_->hostPtr[i];
	}
#endif
	hOffset = minN0;

	//LOG_DEBUG << "[NumericCircularSourceOCLImpulseResponse] minN0=" << minN0 << " maxN0=" << maxN0;
}

} // namespace Lab

#endif // NUMERIC_CIRCULAR_SOURCE_OCL_IMPULSE_RESPONSE_H
