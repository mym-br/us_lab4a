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
#ifndef NUMERIC_ARRAY_OF_RECTANGULAR_SOURCES_OCL_IMPULSE_RESPONSE_H
#define NUMERIC_ARRAY_OF_RECTANGULAR_SOURCES_OCL_IMPULSE_RESPONSE_H

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
#include "XY.h"

#include <CL/cl2.hpp>
#include "NumericRectangularSourceOCLImpulseResponse.h"
#include "OCLAtomic.h"
#include "OCLReduce.h"
#include "OCLUtil.h"

#define NUMERIC_ARRAY_OF_RECTANGULAR_SOURCES_OCL_IMPULSE_RESPONSE_ACCUMULATE_IN_GPU 1
#define NUMERIC_ARRAY_OF_RECTANGULAR_SOURCES_OCL_IMPULSE_RESPONSE_ACCUMULATE_WITH_LOCAL_MEM 1



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
class NumericArrayOfRectangularSourcesOCLImpulseResponse {
public:
	NumericArrayOfRectangularSourcesOCLImpulseResponse(
					TFloat samplingFreq,
					TFloat propagationSpeed,
					TFloat sourceWidth,
					TFloat sourceHeight,
					TFloat subElemSize,
					const std::vector<XY<TFloat>>& elemPos,
					const std::vector<TFloat>& focusDelay /* s */);
	~NumericArrayOfRectangularSourcesOCLImpulseResponse() = default;

	// Return h/c.
	void getImpulseResponse(TFloat x, TFloat y, TFloat z,
				std::size_t& hOffset /* samples */, std::vector<TFloat>& h,
				std::vector<unsigned int>* activeElemList=nullptr);
private:
	enum {
		REDUCE_GROUP_SIZE = 1024,
		INITIAL_H_SIZE = 10000
	};

	NumericArrayOfRectangularSourcesOCLImpulseResponse(const NumericArrayOfRectangularSourcesOCLImpulseResponse&) = delete;
	NumericArrayOfRectangularSourcesOCLImpulseResponse& operator=(const NumericArrayOfRectangularSourcesOCLImpulseResponse&) = delete;
	NumericArrayOfRectangularSourcesOCLImpulseResponse(NumericArrayOfRectangularSourcesOCLImpulseResponse&&) = delete;
	NumericArrayOfRectangularSourcesOCLImpulseResponse& operator=(NumericArrayOfRectangularSourcesOCLImpulseResponse&&) = delete;

	static std::string getKernels();

	TFloat samplingFreq_;
	TFloat propagationSpeed_;
	TFloat subElemWidth_;
	TFloat subElemHeight_;
	unsigned int numElem_;
	unsigned int numSubElem_; // per element

	// OpenCL.
	cl::Context clContext_;
	cl::Program clProgram_;
	cl::CommandQueue clCommandQueue_;
	cl::Buffer subElemXCLBuffer_;
	cl::Buffer subElemYCLBuffer_;
	cl::Buffer elemDelayCLBuffer_;
	cl::Buffer elemPosXCLBuffer_;
	cl::Buffer elemPosYCLBuffer_;
	std::unique_ptr<OCLPinnedHostMem<unsigned int>> activeElemHostMem_;
	cl::Buffer activeElemCLBuffer_;
	cl::Buffer n0CLBuffer_;
	cl::Buffer valueCLBuffer_;
	std::unique_ptr<OCLPinnedHostMem<unsigned int>> minN0HostMem_;
	cl::Buffer minN0CLBuffer_;
	std::unique_ptr<OCLPinnedHostMem<unsigned int>> maxN0HostMem_;
	cl::Buffer maxN0CLBuffer_;
#ifdef NUMERIC_ARRAY_OF_RECTANGULAR_SOURCES_OCL_IMPULSE_RESPONSE_ACCUMULATE_IN_GPU
	std::unique_ptr<OCLPinnedHostMem<TFloat>> hHostMem_;
	cl::Buffer hCLBuffer_;
#else
	std::unique_ptr<OCLPinnedHostMem<unsigned int>> n0HostMem_;
	std::unique_ptr<OCLPinnedHostMem<TFloat>> valueHostMem_;
#endif
};



template<typename TFloat>
NumericArrayOfRectangularSourcesOCLImpulseResponse<TFloat>::NumericArrayOfRectangularSourcesOCLImpulseResponse(
		TFloat samplingFreq,
		TFloat propagationSpeed,
		TFloat sourceWidth,
		TFloat sourceHeight,
		TFloat subElemSize,
		const std::vector<XY<TFloat>>& elemPos,
		const std::vector<TFloat>& focusDelay /* s */)
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

	clContext_ = OCLUtil::initOpenCL();

	std::vector<std::string> kernelStrings = {
					OCLAtomic::code() +
					OCLReduce::code() +
					NumericRectangularSourceOCLImpulseResponse<TFloat>::getKernels() +
					getKernels()};
	clProgram_ = cl::Program(clContext_, kernelStrings);
	std::ostringstream progOpt;
	progOpt << OCL_PROGRAM_BUILD_OPTIONS
			<< " -DMFloat=" << TemplateUtil::typeName<TFloat>();
#ifdef NUMERIC_ARRAY_OF_RECTANGULAR_SOURCES_OCL_IMPULSE_RESPONSE_ACCUMULATE_WITH_LOCAL_MEM
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

	subElemXCLBuffer_   = cl::Buffer(clContext_, CL_MEM_READ_ONLY , sizeof(TFloat) * numSubElem_);
	subElemYCLBuffer_   = cl::Buffer(clContext_, CL_MEM_READ_ONLY , sizeof(TFloat) * numSubElem_);
	elemDelayCLBuffer_  = cl::Buffer(clContext_, CL_MEM_READ_ONLY , sizeof(TFloat) * numElem_);
	elemPosXCLBuffer_   = cl::Buffer(clContext_, CL_MEM_READ_ONLY , sizeof(TFloat) * elemPos.size());
	elemPosYCLBuffer_   = cl::Buffer(clContext_, CL_MEM_READ_ONLY , sizeof(TFloat) * elemPos.size());
	activeElemHostMem_  = std::make_unique<OCLPinnedHostMem<unsigned int>>(clContext_, clCommandQueue_,
									numElem_,
									CL_MAP_WRITE_INVALIDATE_REGION);
	activeElemCLBuffer_ = cl::Buffer(clContext_, CL_MEM_READ_ONLY , sizeof(unsigned int) * numElem_);
	const unsigned int numReduceThreads = OCLUtil::roundUpToMultipleOfGroupSize(numElem_ * numSubElem_, REDUCE_GROUP_SIZE);
	n0CLBuffer_         = cl::Buffer(clContext_, CL_MEM_READ_WRITE, sizeof(unsigned int) * numReduceThreads);
	valueCLBuffer_      = cl::Buffer(clContext_, CL_MEM_READ_WRITE, sizeof(TFloat) * numElem_ * numSubElem_);
	const unsigned int numReduceGroups = numReduceThreads / REDUCE_GROUP_SIZE;
	minN0HostMem_       = std::make_unique<OCLPinnedHostMem<unsigned int>>(clContext_, clCommandQueue_,
									numReduceGroups,
									CL_MAP_READ);
	minN0CLBuffer_      = cl::Buffer(clContext_, CL_MEM_READ_WRITE, sizeof(unsigned int) * numReduceGroups);
	maxN0HostMem_       = std::make_unique<OCLPinnedHostMem<unsigned int>>(clContext_, clCommandQueue_,
									numReduceGroups,
									CL_MAP_READ);
	maxN0CLBuffer_      = cl::Buffer(clContext_, CL_MEM_READ_WRITE, sizeof(unsigned int) * numReduceGroups);
#ifdef NUMERIC_ARRAY_OF_RECTANGULAR_SOURCES_OCL_IMPULSE_RESPONSE_ACCUMULATE_IN_GPU
	hHostMem_           = std::make_unique<OCLPinnedHostMem<TFloat>>(clContext_, clCommandQueue_,
									INITIAL_H_SIZE,
									CL_MAP_READ);
	hCLBuffer_          = cl::Buffer(clContext_, CL_MEM_READ_WRITE, sizeof(TFloat) * INITIAL_H_SIZE);
#else
	n0HostMem_          = std::make_unique<OCLPinnedHostMem<unsigned int>>(clContext_, clCommandQueue_,
									numReduceThreads,
									CL_MAP_READ);
	valueHostMem_       = std::make_unique<OCLPinnedHostMem<TFloat>>(clContext_, clCommandQueue_,
									numElem_ * numSubElem_,
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

	std::vector<TFloat> elemDelay(numElem_);
	for (unsigned int i = 0; i < elemDelay.size(); ++i) {
		elemDelay[i] = focusDelay[i] * samplingFreq_;
	}
	clCommandQueue_.enqueueWriteBuffer(
		elemDelayCLBuffer_, CL_TRUE /* blocking */, 0 /* offset */,
		elemDelay.size() * sizeof(TFloat), elemDelay.data());

	std::vector<TFloat> elemPosX(elemPos.size());
	std::vector<TFloat> elemPosY(elemPos.size());
	for (unsigned int i = 0; i < elemPos.size(); ++i) {
		elemPosX[i] = elemPos[i].x;
		elemPosY[i] = elemPos[i].y;
	}
	clCommandQueue_.enqueueWriteBuffer(
		elemPosXCLBuffer_, CL_TRUE /* blocking */, 0 /* offset */,
		elemPosX.size() * sizeof(TFloat), elemPosX.data());
	clCommandQueue_.enqueueWriteBuffer(
		elemPosYCLBuffer_, CL_TRUE /* blocking */, 0 /* offset */,
		elemPosY.size() * sizeof(TFloat), elemPosY.data());

	LOG_DEBUG << "[NumericArrayOfRectangularSourcesOCLImpulseResponse] numElemX=" << numElemX << " numElemY=" << numElemY;
}

template<typename TFloat>
void
NumericArrayOfRectangularSourcesOCLImpulseResponse<TFloat>::getImpulseResponse(
							TFloat x,
							TFloat y,
							TFloat z,
							std::size_t& hOffset,
							std::vector<TFloat>& h,
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
				activeElemHostMem_->hostPtr[i] = (*activeElemList)[i];
			}
		}
	} else {
		for (unsigned int i = 0; i < numElem_; ++i) {
			activeElemHostMem_->hostPtr[i] = i;
		}
	}
	clCommandQueue_.enqueueWriteBuffer(
		activeElemCLBuffer_, CL_TRUE /* blocking */, 0 /* offset */,
		activeElemHostMem_->sizeInBytes, activeElemHostMem_->hostPtr);

	try {
		cl::Kernel kernel(clProgram_, "numericArraySourceIRKernel");
		kernel.setArg( 0, numElem_);
		kernel.setArg( 1, numSubElem_);
		kernel.setArg( 2, x);
		kernel.setArg( 3, y);
		kernel.setArg( 4, z);
		kernel.setArg( 5, samplingFreq_ / propagationSpeed_);
		kernel.setArg( 6, samplingFreq_ * subElemWidth_ * subElemHeight_ / (TFloat(2.0 * pi) * propagationSpeed_));
		kernel.setArg( 7, activeElemCLBuffer_);
		kernel.setArg( 8, elemDelayCLBuffer_);
		kernel.setArg( 9, elemPosXCLBuffer_);
		kernel.setArg(10, elemPosYCLBuffer_);
		kernel.setArg(11, subElemXCLBuffer_);
		kernel.setArg(12, subElemYCLBuffer_);
		kernel.setArg(13, n0CLBuffer_);
		kernel.setArg(14, valueCLBuffer_);

		const unsigned int groupSize = REDUCE_GROUP_SIZE;
		const unsigned int globalN0 = OCLUtil::roundUpToMultipleOfGroupSize(numSubElem_, groupSize);

		clCommandQueue_.enqueueNDRangeKernel(
			kernel,
			cl::NullRange, // offset
			cl::NDRange(globalN0), // global
			cl::NDRange(groupSize), // local
			nullptr /* previous events */, nullptr);
	} catch (cl::Error& e) {
		THROW_EXCEPTION(OCLException, "[numericArraySourceIRKernel] OpenCL error: " << e.what() << " (" << e.err() << ").");
	}

	const unsigned int reduceGlobalN0 = OCLUtil::roundUpToMultipleOfGroupSize(numElem_ * numSubElem_, REDUCE_GROUP_SIZE);
	try {
		cl::Kernel kernel(clProgram_, "groupReduceMinMaxKernel");
		kernel.setArg(0, n0CLBuffer_);
		kernel.setArg(1, minN0CLBuffer_);
		kernel.setArg(2, maxN0CLBuffer_);
		kernel.setArg(3, numElem_ * numSubElem_);
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
			minN0CLBuffer_, CL_TRUE /* blocking */, 0 /* offset */,
			sizeof(unsigned int), minN0HostMem_->hostPtr);
		clCommandQueue_.enqueueReadBuffer(
			maxN0CLBuffer_, CL_TRUE /* blocking */, 0 /* offset */,
			sizeof(unsigned int), maxN0HostMem_->hostPtr);
	} catch (cl::Error& e) {
		THROW_EXCEPTION(OCLException, "[Read minN0, maxN0] OpenCL error: " << e.what() << " (" << e.err() << ").");
	}

	const unsigned int minN0 = *(minN0HostMem_->hostPtr);
	const unsigned int maxN0 = *(maxN0HostMem_->hostPtr);
	const unsigned int hSize = maxN0 - minN0 + 1U;

#ifdef NUMERIC_ARRAY_OF_RECTANGULAR_SOURCES_OCL_IMPULSE_RESPONSE_ACCUMULATE_IN_GPU
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
		kernel.setArg(0, numElem_ * numSubElem_);
		kernel.setArg(1, minN0);
		kernel.setArg(2, n0CLBuffer_);
		kernel.setArg(3, valueCLBuffer_);
		kernel.setArg(4, hCLBuffer_);
#ifdef NUMERIC_ARRAY_OF_RECTANGULAR_SOURCES_OCL_IMPULSE_RESPONSE_ACCUMULATE_WITH_LOCAL_MEM
		kernel.setArg(5, hSize);
		kernel.setArg(6, cl::Local(hSize * sizeof(TFloat)));

		const unsigned int groupSize = 1024;
#else
		const unsigned int groupSize = 64;
#endif
		const unsigned int globalN0 = OCLUtil::roundUpToMultipleOfGroupSize(numElem_ * numSubElem_, groupSize);

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
		numElem_ * numSubElem_ * sizeof(unsigned int), n0HostMem_->hostPtr);
	clCommandQueue_.enqueueReadBuffer(
		valueCLBuffer_, CL_TRUE /* blocking */, 0 /* offset */,
		numElem_ * numSubElem_ * sizeof(TFloat), valueHostMem_->hostPtr);

	h.assign(hSize, 0);
	for (unsigned int i = 0, iEnd = numElem_ * numSubElem_; i < iEnd; ++i) {
		h[n0HostMem_->hostPtr[i] - minN0] += valueHostMem_->hostPtr[i];
	}
#endif
	hOffset = minN0;

	//LOG_DEBUG << "[NumericArrayOfRectangularSourcesOCLImpulseResponse] minN0=" << minN0 << " maxN0=" << maxN0;
}

template<typename TFloat>
std::string
NumericArrayOfRectangularSourcesOCLImpulseResponse<TFloat>::getKernels()
{
	return R"CLC(

__kernel
void
numericArraySourceIRKernel(
		unsigned int numElem,
		unsigned int numSubElem,
		MFloat x,
		MFloat y,
		MFloat z,
		MFloat k1,
		MFloat k2,
		__constant unsigned int* activeElem,
		__constant MFloat* elemDelay,
		__constant MFloat* elemPosX,
		__constant MFloat* elemPosY,
		__global MFloat* subElemX,
		__global MFloat* subElemY,
		__global unsigned int* n0,
		__global MFloat* value)
{
	const unsigned int subElemIdx = get_global_id(0);
	if (subElemIdx >= numSubElem) {
		// Get data from the first sub-element, to help min/max(n0).
		const unsigned int activeElemIdx = activeElem[0];
		const MFloat dx = x - subElemX[0] - elemPosX[activeElemIdx];
		const MFloat dy = y - subElemY[0] - elemPosY[activeElemIdx];
		const MFloat r = sqrt(dx * dx + dy * dy + z * z);
		const unsigned int idx = (numElem - 1U) * numSubElem + subElemIdx;
		n0[idx] = rint(r * k1 + elemDelay[0]);
		return;
	}

	for (int i = 0; i < numElem; ++i) {
		const unsigned int activeElemIdx = activeElem[i];
		const MFloat dx = x - subElemX[subElemIdx] - elemPosX[activeElemIdx];
		const MFloat dy = y - subElemY[subElemIdx] - elemPosY[activeElemIdx];
		const MFloat r = sqrt(dx * dx + dy * dy + z * z);
		const unsigned int idx = i * numSubElem + subElemIdx;
		n0[idx] = rint(r * k1 + elemDelay[i]);
		value[idx] = k2 / r;
	}
}

)CLC";
}

} // namespace Lab

#endif // NUMERIC_ARRAY_OF_RECTANGULAR_SOURCES_OCL_IMPULSE_RESPONSE_H