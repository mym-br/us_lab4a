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
#ifndef NUMERIC_RECTANGULAR_SOURCE_OCL_IMPULSE_RESPONSE_H
#define NUMERIC_RECTANGULAR_SOURCE_OCL_IMPULSE_RESPONSE_H

#include <cmath> /* ceil */
#include <cstddef> /* std::size_t */
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <type_traits> /* is_same */
#include <vector>

#include "Exception.h"
#include "Log.h"
#include "Util.h" /* pi */

#include <CL/cl2.hpp>
#include "OCLAtomic.h"
#include "OCLReduce.h"
#include "OCLUtil.h"

#define NUMERIC_RECTANGULAR_SOURCE_OCL_IMPULSE_RESPONSE_ACCUMULATE_IN_GPU 1
#define NUMERIC_RECTANGULAR_SOURCE_OCL_IMPULSE_RESPONSE_ACCUMULATE_WITH_LOCAL_MEM 1



namespace Lab {

// Calculate the acoustic field generated by a flat rectangular surface,
// using the numeric solution provided by:
//
// Piwakowski, B.
// Delannoy, B.
// Method for computing spatial pulse response: Time-domain approach
// J. Acoust. Soc. Am., vol. 86, no. 6, pp. 2422-2432, 1989.
// DOI: 10.1121/1.398449
//
// See also:
// Lasota, H.
// Salamon, R.
// Delannoy, B.
// Acoustic diffraction analysis by the impulse response method: A line impulse response approach},
// J. Acoust. Soc. Am., vol. 76, no. 1, pp. 280-290, 1984.
// DOI: 10.1121/1.391115
//
// Note:
// - The source is surrounded by a rigid baffle.
template<typename TFloat>
class NumericRectangularSourceOCLImpulseResponse {
public:
	NumericRectangularSourceOCLImpulseResponse(
					TFloat samplingFreq,
					TFloat propagationSpeed,
					TFloat sourceWidth,
					TFloat sourceHeight,
					TFloat subElemSize);
	~NumericRectangularSourceOCLImpulseResponse() = default;

	// Return h/c.
	void getImpulseResponse(TFloat x, TFloat y, TFloat z,
				std::size_t& hOffset /* samples */, std::vector<TFloat>& h);
private:
	enum {
		REDUCE_GROUP_SIZE = 1024,
		INITIAL_H_SIZE = 10000
	};

	NumericRectangularSourceOCLImpulseResponse(const NumericRectangularSourceOCLImpulseResponse&) = delete;
	NumericRectangularSourceOCLImpulseResponse& operator=(const NumericRectangularSourceOCLImpulseResponse&) = delete;
	NumericRectangularSourceOCLImpulseResponse(NumericRectangularSourceOCLImpulseResponse&&) = delete;
	NumericRectangularSourceOCLImpulseResponse& operator=(NumericRectangularSourceOCLImpulseResponse&&) = delete;

	std::string getKernels() const;

	TFloat samplingFreq_;
	TFloat propagationSpeed_;
	TFloat subElemWidth_;
	TFloat subElemHeight_;
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
#ifdef NUMERIC_RECTANGULAR_SOURCE_OCL_IMPULSE_RESPONSE_ACCUMULATE_IN_GPU
	std::unique_ptr<OCLPinnedHostMem<TFloat>> hHostMem_;
	cl::Buffer hCLBuffer_;
#else
	std::unique_ptr<OCLPinnedHostMem<unsigned int>> n0HostMem_;
	std::unique_ptr<OCLPinnedHostMem<TFloat>> valueHostMem_;
#endif
};



template<typename TFloat>
NumericRectangularSourceOCLImpulseResponse<TFloat>::NumericRectangularSourceOCLImpulseResponse(
		TFloat samplingFreq,
		TFloat propagationSpeed,
		TFloat sourceWidth,
		TFloat sourceHeight,
		TFloat subElemSize)
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

	clContext_ = OCLUtil::initOpenCL();

	std::vector<std::string> kernelStrings = {OCLAtomic::code() + OCLReduce::code() + getKernels()};
	clProgram_ = cl::Program(clContext_, kernelStrings);
	std::ostringstream progOpt;
	progOpt << OCL_PROGRAM_BUILD_OPTIONS
			<< " -DMFloat=" << TemplateUtil::typeName<TFloat>();
#ifdef NUMERIC_RECTANGULAR_SOURCE_OCL_IMPULSE_RESPONSE_ACCUMULATE_WITH_LOCAL_MEM
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

	subElemXCLBuffer_ = cl::Buffer(clContext_, CL_MEM_READ_ONLY, sizeof(TFloat) * numSubElem_);
	subElemYCLBuffer_ = cl::Buffer(clContext_, CL_MEM_READ_ONLY, sizeof(TFloat) * numSubElem_);
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
#ifdef NUMERIC_RECTANGULAR_SOURCE_OCL_IMPULSE_RESPONSE_ACCUMULATE_IN_GPU
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
	std::vector<TFloat> subElemX(numSubElem_);
	std::vector<TFloat> subElemY(numSubElem_);
	const TFloat halfW = 0.5 * (numElemX - 1);
	const TFloat halfH = 0.5 * (numElemY - 1);
	for (unsigned int iy = 0; iy < numElemY; ++iy) {
		for (unsigned int ix = 0; ix < numElemX; ++ix) {
			unsigned int idx = iy * numElemX + ix;
			subElemX[idx] = (ix - halfW) * subElemWidth_;
			subElemY[idx] = (iy - halfH) * subElemHeight_;
		}
	}
	clCommandQueue_.enqueueWriteBuffer(
		subElemXCLBuffer_, CL_TRUE /* blocking */, 0 /* offset */,
		subElemX.size() * sizeof(TFloat), subElemX.data());
	clCommandQueue_.enqueueWriteBuffer(
		subElemYCLBuffer_, CL_TRUE /* blocking */, 0 /* offset */,
		subElemY.size() * sizeof(TFloat), subElemY.data());

	LOG_DEBUG << "[NumericRectangularSourceOCLImpulseResponse] numElemX=" << numElemX << " numElemY=" << numElemY;
}

template<typename TFloat>
void
NumericRectangularSourceOCLImpulseResponse<TFloat>::getImpulseResponse(
							TFloat x,
							TFloat y,
							TFloat z,
							std::size_t& hOffset,
							std::vector<TFloat>& h)
{
	try {
		cl::Kernel kernel(clProgram_, "numericSourceIRKernel");
		kernel.setArg(0, numSubElem_);
		kernel.setArg(1, x);
		kernel.setArg(2, y);
		kernel.setArg(3, z);
		kernel.setArg(4, samplingFreq_ / propagationSpeed_);
		kernel.setArg(5, samplingFreq_ * subElemWidth_ * subElemHeight_ / (TFloat(2.0 * pi) * propagationSpeed_));
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

	clCommandQueue_.enqueueReadBuffer(
		minN0CLBuffer_, CL_TRUE /* blocking */, 0 /* offset */,
		sizeof(unsigned int), minN0HostMem_->hostPtr);
	const unsigned int minN0 = *(minN0HostMem_->hostPtr);
	clCommandQueue_.enqueueReadBuffer(
		maxN0CLBuffer_, CL_TRUE /* blocking */, 0 /* offset */,
		sizeof(unsigned int), maxN0HostMem_->hostPtr);
	const unsigned int maxN0 = *(maxN0HostMem_->hostPtr);

	const unsigned int hSize = maxN0 - minN0 + 1U;

#ifdef NUMERIC_RECTANGULAR_SOURCE_OCL_IMPULSE_RESPONSE_ACCUMULATE_IN_GPU
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
#ifdef NUMERIC_RECTANGULAR_SOURCE_OCL_IMPULSE_RESPONSE_ACCUMULATE_WITH_LOCAL_MEM
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
			hCLBuffer_, CL_TRUE /* blocking */, 0 /* offset */,
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
		n0CLBuffer_, CL_TRUE /* blocking */, 0 /* offset */,
		numSubElem_ * sizeof(unsigned int), n0HostMem_->hostPtr);
	clCommandQueue_.enqueueReadBuffer(
		valueCLBuffer_, CL_TRUE /* blocking */, 0 /* offset */,
		numSubElem_ * sizeof(TFloat), valueHostMem_->hostPtr);

	h.assign(hSize, 0);
	for (unsigned int i = 0; i < numSubElem_; ++i) {
		h[n0HostMem_->hostPtr[i] - minN0] += valueHostMem_->hostPtr[i];
	}
#endif
	hOffset = minN0;

	//LOG_DEBUG << "[NumericRectangularSourceOCLImpulseResponse] minN0=" << minN0 << " maxN0=" << maxN0;
}

template<typename TFloat>
std::string
NumericRectangularSourceOCLImpulseResponse<TFloat>::getKernels() const
{
	return R"CLC(

__kernel
void
numericSourceIRKernel(
		unsigned int numSubElem,
		MFloat x,
		MFloat y,
		MFloat z,
		MFloat k1,
		MFloat k2,
		__global MFloat* subElemX,
		__global MFloat* subElemY,
		__global unsigned int* n0,
		__global MFloat* value)
{
	const unsigned int subElemIdx = get_global_id(0);
	if (subElemIdx >= numSubElem) {
		// Get data from the first sub-element, to help min/max(n0).
		const MFloat dx = x - subElemX[0];
		const MFloat dy = y - subElemY[0];
		const MFloat r = sqrt(dx * dx + dy * dy + z * z);
		n0[subElemIdx] = rint(r * k1);
		return;
	}

	const MFloat dx = x - subElemX[subElemIdx];
	const MFloat dy = y - subElemY[subElemIdx];
	const MFloat r = sqrt(dx * dx + dy * dy + z * z);
	n0[subElemIdx] = rint(r * k1);
	value[subElemIdx] = k2 / r;
}

#ifdef WITH_LOCAL_MEM
__kernel
void
accumulateIRSamplesKernel(
		unsigned int numSubElem,
		unsigned int minN0,
		__global unsigned int* n0,
		__global MFloat* value,
		__global MFloat* h,
		unsigned int hSize,
		__local MFloat* localH)
{
	for (int hIdx = get_local_id(0); hIdx < hSize; hIdx += get_local_size(0)) {
		localH[hIdx] = 0;
	}

	barrier(CLK_LOCAL_MEM_FENCE);

	const unsigned int subElemIdx = get_global_id(0);
	if (subElemIdx < numSubElem) {
		// Different sub-elements may have the same value for n0.
# ifdef SINGLE_PREC
		atomicAddLocalFloat(localH + n0[subElemIdx] - minN0, value[subElemIdx]);
# else
		atomicAddLocalDouble(localH + n0[subElemIdx] - minN0, value[subElemIdx]);
# endif
	}

	barrier(CLK_LOCAL_MEM_FENCE);

	for (int hIdx = get_local_id(0); hIdx < hSize; hIdx += get_local_size(0)) {
		const MFloat sample = localH[hIdx];
		if (sample != 0) {
# ifdef SINGLE_PREC
			atomicAddGlobalFloat(h + hIdx, sample);
# else
			atomicAddGlobalDouble(h + hIdx, sample);
# endif
		}
	}
}
#else
__kernel
void
accumulateIRSamplesKernel(
		unsigned int numSubElem,
		unsigned int minN0,
		__global unsigned int* n0,
		__global MFloat* value,
		__global MFloat* h)
{
	const unsigned int subElemIdx = get_global_id(0);
	if (subElemIdx >= numSubElem) {
		return;
	}

	// Different sub-elements may have the same value for n0.
# ifdef SINGLE_PREC
	atomicAddGlobalFloat(h + n0[subElemIdx] - minN0, value[subElemIdx]);
# else
	atomicAddGlobalDouble(h + n0[subElemIdx] - minN0, value[subElemIdx]);
# endif
}
#endif

)CLC";
}

} // namespace Lab

#endif // NUMERIC_RECTANGULAR_SOURCE_OCL_IMPULSE_RESPONSE_H
