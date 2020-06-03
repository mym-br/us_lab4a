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

#ifndef VECTORIAL_STA_OCL_PROCESSOR_H
#define VECTORIAL_STA_OCL_PROCESSOR_H

#include <algorithm> /* copy */
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <tbb/cache_aligned_allocator.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/partitioner.h>
#include <tbb/tbb.h>

#include "ArrayGeometry.h"
#include "ArrayProcessor.h"
#include "CoherenceFactor.h"
#include "Exception.h"
#include "ExecutionTimeMeasurement.h"
#include "HilbertEnvelope.h"
#include "Interpolator.h"
#include "Log.h"
#include "STAAcquisition.h"
#include "STAConfiguration.h"
#include "TemplateUtil.h"
#include "Timer.h"
#include "Util.h"
#include "XYZValueFactor.h"

#include <CL/cl2.hpp>
#include "OCLCoherenceFactor.h"
#include "OCLGeometry.h"
#include "OCLUtil.h"

// Depends on the signal.
// 1.0 --> pi radian / sample at the original sampling rate.
#define VECTORIAL_STA_OCL_PROCESSOR_UPSAMP_FILTER_HALF_TRANSITION_WIDTH (0.2)



namespace Lab {

// STA image formation, using analytic signals (each sample is a real-imag vector).
//
// Processing steps:
//   Signal preparation                                                    - CPU
//   Calculate delays                                                      - OpenCL
//   Apply delays, use apodization, [apply PCF] and accumulate the samples - OpenCL
//
// x = 0 is at the center of the element group.
// y = 0
// z = 0 is at the surface of the array.
template<typename TFloat>
class VectorialSTAOCLProcessor : public ArrayProcessor<XYZValueFactor<TFloat>> {
public:
	VectorialSTAOCLProcessor(
			const STAConfiguration<TFloat>& config,
			STAAcquisition<TFloat>& acquisition,
			unsigned int upsamplingFactor,
			AnalyticSignalCoherenceFactorProcessor<TFloat>& coherenceFactor,
			TFloat peakOffset,
			const std::vector<TFloat>& rxApod);
	virtual ~VectorialSTAOCLProcessor() = default;

	virtual void prepare(unsigned int baseElement);
	virtual void process(Matrix<XYZValueFactor<TFloat>>& gridData);

#ifdef LAB_ENABLE_EXECUTION_TIME_MEASUREMENT
	MeasurementList<double> tAcquisitionML;
	MeasurementList<double> tPrepareDataML;
	MeasurementList<double> tCalculateDelaysML;
	MeasurementList<double> tDelaySumML;
	void execTimeMeasReset(unsigned int n) {
		tAcquisitionML.reset(    EXECUTION_TIME_MEASUREMENT_ITERATIONS * n);
		tPrepareDataML.reset(    EXECUTION_TIME_MEASUREMENT_ITERATIONS * n);
		tCalculateDelaysML.reset(EXECUTION_TIME_MEASUREMENT_ITERATIONS);
		tDelaySumML.reset(       EXECUTION_TIME_MEASUREMENT_ITERATIONS);
	}
	void execTimeMeasShowResults(unsigned int n) {
		EXECUTION_TIME_MEASUREMENT_LOG_TIMES_X_N("tAcquisition:    ", tAcquisitionML, n);
		EXECUTION_TIME_MEASUREMENT_LOG_TIMES_X_N("tPrepareData:    ", tPrepareDataML, n);
		EXECUTION_TIME_MEASUREMENT_LOG_TIMES(    "tCalculateDelays:", tCalculateDelaysML);
		EXECUTION_TIME_MEASUREMENT_LOG_TIMES(    "tDelaySum:       ", tDelaySumML);
	}
#endif

private:
	struct PrepareDataThreadData {
		Interpolator<TFloat> interpolator;
		HilbertEnvelope<TFloat> envelope;
		std::vector<TFloat, tbb::cache_aligned_allocator<TFloat>> signal;
	};
	struct PrepareData;

	VectorialSTAOCLProcessor(const VectorialSTAOCLProcessor&) = delete;
	VectorialSTAOCLProcessor& operator=(const VectorialSTAOCLProcessor&) = delete;
	VectorialSTAOCLProcessor(VectorialSTAOCLProcessor&&) = delete;
	VectorialSTAOCLProcessor& operator=(VectorialSTAOCLProcessor&&) = delete;

	static std::string getKernels();

	const STAConfiguration<TFloat>& config_;
	unsigned int deadZoneSamplesUp_;
	STAAcquisition<TFloat>& acquisition_;
	unsigned int upsamplingFactor_;
	AnalyticSignalCoherenceFactorProcessor<TFloat>& coherenceFactor_;
	typename STAAcquisition<TFloat>::AcquisitionDataType acqData_;
	unsigned int signalLength_;
	TFloat signalOffset_;
	std::vector<TFloat, tbb::cache_aligned_allocator<TFloat>> xArray_;
	bool initialized_;
	const std::vector<TFloat> rxApod_;
	std::unique_ptr<tbb::enumerable_thread_specific<PrepareDataThreadData>> prepareDataTLS_;

	unsigned int gridN1_;
	unsigned int gridN2_;

	// OpenCL.
	bool clDataInitialized_;
	cl::Context clContext_;
	cl::Program clProgram_;
	cl::CommandQueue clCommandQueue_;
	std::unique_ptr<OCLPinnedHostMem<TFloat[2]>> gridXZHostMem_;
	cl::Buffer gridXZCLBuffer_;
	cl::Buffer rxApodCLBuffer_;
	cl::Buffer xArrayCLBuffer_;
	std::unique_ptr<OCLPinnedHostMem<TFloat[2]>> signalTensorHostMem_;
	cl::Buffer signalTensorCLBuffer_;
	cl::Buffer delayTensorCLBuffer_;
	std::unique_ptr<OCLPinnedHostMem<TFloat>> gridValueHostMem_;
	cl::Buffer gridValueCLBuffer_;
	std::unique_ptr<OCLPinnedHostMem<TFloat>> gridFactorHostMem_;
	cl::Buffer gridFactorCLBuffer_;
};



template<typename TFloat>
struct VectorialSTAOCLProcessor<TFloat>::PrepareData {
	void operator()(const tbb::blocked_range<unsigned int>& r) const {
		auto& local = prepareDataTLS.local();

		for (unsigned int rxElem = r.begin(); rxElem != r.end(); ++rxElem) {
			if (upsamplingFactor > 1) {
				local.interpolator.interpolate(&acqData(rxElem, 0), samplesPerChannelLow, local.signal.data());
			} else {
				auto range = acqData.range2(rxElem);
				std::copy(range.begin(), range.end(), local.signal.begin());
			}

			Util::removeDC(local.signal.data(), local.signal.size(), deadZoneSamplesUp);

			local.envelope.getAnalyticSignal(
					local.signal.data(),
					local.signal.size(),
					signalTensor + ((txElemIdx * signalTensorN2) + rxElem) * signalTensorN3);
		}
	}

	// Use only small types or references/pointers, because this object may be copied many times.
	const unsigned int samplesPerChannelLow;
	const Matrix<TFloat>& acqData;
	const unsigned int upsamplingFactor;
	tbb::enumerable_thread_specific<PrepareDataThreadData>& prepareDataTLS;
	TFloat (*signalTensor)[2];
	const unsigned int signalTensorN2;
	const unsigned int signalTensorN3;
	const unsigned int txElemIdx;
	const unsigned int deadZoneSamplesUp;
};

template<typename TFloat>
VectorialSTAOCLProcessor<TFloat>::VectorialSTAOCLProcessor(
			const STAConfiguration<TFloat>& config,
			STAAcquisition<TFloat>& acquisition,
			unsigned int upsamplingFactor,
			AnalyticSignalCoherenceFactorProcessor<TFloat>& coherenceFactor,
			TFloat peakOffset,
			const std::vector<TFloat>& rxApod)
		: config_(config)
		, deadZoneSamplesUp_((upsamplingFactor * config.samplingFrequency) * 2.0 * config.deadZoneM / config.propagationSpeed)
		, acquisition_(acquisition)
		, upsamplingFactor_(upsamplingFactor)
		, coherenceFactor_(coherenceFactor)
		, signalLength_()
		, signalOffset_()
		, initialized_()
		, rxApod_(rxApod)
		, gridN1_()
		, gridN2_()
		, clDataInitialized_()
{
	signalOffset_ = (config_.samplingFrequency * upsamplingFactor_) * (peakOffset / config_.centerFrequency);
	LOG_DEBUG << "signalOffset_: " << signalOffset_;

	clContext_ = OCLUtil::initOpenCL();

	std::vector<std::string> kernelStrings = {OCLCoherenceFactor::code() + OCLGeometry::code() + getKernels()};
	clProgram_ = cl::Program(clContext_, kernelStrings);
	std::ostringstream progOpt;
	progOpt << LAB_OPENCL_PROGRAM_BUILD_OPTIONS
			<< " -DNUM_RX_ELEM="        << config.numElements
			<< " -DSTATISTICS_N="       << config.numElements
			<< " -DCOHERENCE_FACTOR_N=" << config.numElements
			<< " -DMFloat="             << TemplateUtil::typeName<TFloat>();
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
}

template<typename TFloat>
void
VectorialSTAOCLProcessor<TFloat>::prepare(unsigned int baseElement)
{
	acquisition_.prepare(baseElement);
}

template<typename TFloat>
void
VectorialSTAOCLProcessor<TFloat>::process(Matrix<XYZValueFactor<TFloat>>& gridData)
{
	//LOG_DEBUG << "BEGIN ========== VectorialSTAOCLProcessor::process ==========";

	Util::resetValueFactor(gridData.begin(), gridData.end());

	const unsigned int numTxElem = config_.lastTxElem - config_.firstTxElem + 1U;
	const unsigned int numCols = gridData.n1();
	const unsigned int numRows = gridData.n2();
	if (numCols == 0) {
		THROW_EXCEPTION(InvalidValueException, "Zero columns in the grid.");
	}

	if (xArray_.empty()) {
		ArrayGeometry::getElementXCentered2D(config_.numElements, config_.pitch, xArray_);
	}

	// Prepare the signal tensor.
	for (unsigned int txElem = config_.firstTxElem; txElem <= config_.lastTxElem; ++txElem) {
		//LOG_INFO << "ACQ/PREP txElem: " << txElem << " <= " << config_.lastTxElem;

#ifdef LAB_ENABLE_EXECUTION_TIME_MEASUREMENT
		Timer acquisitionTimer;
#endif
		acquisition_.execute(txElem, acqData_);

#ifdef LAB_ENABLE_EXECUTION_TIME_MEASUREMENT
		tAcquisitionML.put(acquisitionTimer.getTime());
#endif
		if (!initialized_) {
			signalLength_ = acqData_.n2() * upsamplingFactor_;
			LOG_DEBUG << "signalLength: " << signalLength_;
			if (deadZoneSamplesUp_ >= signalLength_) {
				THROW_EXCEPTION(InvalidValueException,
						"Invalid dead zone: deadZoneSamplesUp_ (" << deadZoneSamplesUp_ <<
						") >= signalLength (" << signalLength_ << ").");
			}

			PrepareDataThreadData prepareDataThreadData;
			if (upsamplingFactor_ > 1) {
				prepareDataThreadData.interpolator.prepare(upsamplingFactor_, VECTORIAL_STA_OCL_PROCESSOR_UPSAMP_FILTER_HALF_TRANSITION_WIDTH);
			}
			prepareDataThreadData.signal.resize(signalLength_);
			prepareDataTLS_ = std::make_unique<tbb::enumerable_thread_specific<PrepareDataThreadData>>(prepareDataThreadData);

			initialized_ = true;
		}
		if (!clDataInitialized_) {
			gridN1_ = gridData.n1();
			gridN2_ = gridData.n2();

			gridXZHostMem_        = std::make_unique<OCLPinnedHostMem<TFloat[2]>>(clContext_, clCommandQueue_,
							gridN1_ * gridN2_,
							CL_MAP_WRITE_INVALIDATE_REGION);
			gridXZCLBuffer_       = cl::Buffer(clContext_, CL_MEM_READ_ONLY , sizeof(TFloat) * 2 * gridN1_ * gridN2_);
			rxApodCLBuffer_       = cl::Buffer(clContext_, CL_MEM_READ_ONLY , sizeof(TFloat) * rxApod_.size());
			xArrayCLBuffer_       = cl::Buffer(clContext_, CL_MEM_READ_ONLY , sizeof(TFloat) * xArray_.size());
			signalTensorHostMem_  = std::make_unique<OCLPinnedHostMem<TFloat[2]>>(clContext_, clCommandQueue_,
							numTxElem * config_.numElements * signalLength_,
							CL_MAP_WRITE_INVALIDATE_REGION);
			signalTensorCLBuffer_ = cl::Buffer(clContext_, CL_MEM_READ_ONLY , sizeof(TFloat) * 2 * numTxElem * config_.numElements * signalLength_);
			delayTensorCLBuffer_  = cl::Buffer(clContext_, CL_MEM_READ_WRITE, sizeof(TFloat) * config_.numElements * numCols * numRows);
			gridValueHostMem_     = std::make_unique<OCLPinnedHostMem<TFloat>>(clContext_, clCommandQueue_,
							numCols * numRows,
							CL_MAP_READ);
			gridValueCLBuffer_    = cl::Buffer(clContext_, CL_MEM_READ_WRITE, sizeof(TFloat) * numCols * numRows);
			if (coherenceFactor_.enabled()) {
				gridFactorHostMem_  = std::make_unique<OCLPinnedHostMem<TFloat>>(clContext_, clCommandQueue_,
								numCols * numRows,
								CL_MAP_READ);
				gridFactorCLBuffer_ = cl::Buffer(clContext_, CL_MEM_READ_WRITE, sizeof(TFloat) * numCols * numRows);
			}

			clDataInitialized_ = true;
		} else {
			if (gridN1_ != gridData.n1() ||
			    gridN2_ != gridData.n2()) {
				THROW_EXCEPTION(InvalidParameterException, "Data size mismatch.");
			}
		}

		const unsigned int samplesPerChannelLow = acqData_.n2();

#ifdef LAB_ENABLE_EXECUTION_TIME_MEASUREMENT
		Timer prepareDataTimer;
#endif
		PrepareData prepareDataOp = {
			samplesPerChannelLow,
			acqData_,
			upsamplingFactor_,
			*prepareDataTLS_,
			signalTensorHostMem_->hostPtr,
			config_.numElements,
			signalLength_,
			txElem - config_.firstTxElem,
			deadZoneSamplesUp_
		};
		//tbb::parallel_for(tbb::blocked_range<unsigned int>(0, config_.numElements), prepareDataOp);
		tbb::parallel_for(tbb::blocked_range<unsigned int>(0, config_.numElements, 1 /* grain size */), prepareDataOp, tbb::simple_partitioner());
		//prepareDataOp(tbb::blocked_range<unsigned int>(0, config_.numElements)); // single-thread

#ifdef LAB_ENABLE_EXECUTION_TIME_MEASUREMENT
		tPrepareDataML.put(prepareDataTimer.getTime());
#endif
	}

	if (!initialized_) THROW_EXCEPTION(InvalidStateException, "Not initialized.");

	// Prepare OpenCL buffers.
	Timer prepareBuffersTimer;
	for (std::size_t i = 0; i < gridData.size(); ++i) {
		const XYZValueFactor<TFloat>& point = *(gridData.data() + i);
		TFloat (*dest)[2] = gridXZHostMem_->hostPtr + i;
		(*dest)[0] = point.x;
		(*dest)[1] = point.z;
	}
	clCommandQueue_.enqueueWriteBuffer(
		gridXZCLBuffer_, CL_BLOCKING, 0 /* offset */,
		gridXZHostMem_->sizeInBytes, gridXZHostMem_->hostPtr);
	clCommandQueue_.enqueueWriteBuffer(
		rxApodCLBuffer_, CL_BLOCKING, 0 /* offset */,
		rxApod_.size() * sizeof(TFloat), rxApod_.data());
	clCommandQueue_.enqueueWriteBuffer(
		xArrayCLBuffer_, CL_BLOCKING, 0 /* offset */,
		xArray_.size() * sizeof(TFloat), xArray_.data());
	clCommandQueue_.enqueueWriteBuffer(
		signalTensorCLBuffer_, CL_BLOCKING, 0 /* offset */,
		signalTensorHostMem_->sizeInBytes, signalTensorHostMem_->hostPtr);
	clCommandQueue_.enqueueFillBuffer(
		gridValueCLBuffer_, (TFloat) 0, 0 /* offset */, gridValueHostMem_->sizeInBytes);
	LOG_DEBUG << "PREPARE BUFFERS " << prepareBuffersTimer.getTime();

	try {
#ifdef LAB_ENABLE_EXECUTION_TIME_MEASUREMENT
		Timer calculateDelaysTimer;
#endif
		cl::Kernel kernel(clProgram_, "calculateDelaysSTAKernel");
		kernel.setArg(0, numCols);
		kernel.setArg(1, numRows);
		kernel.setArg(2, config_.samplingFrequency * upsamplingFactor_ / config_.propagationSpeed);
		kernel.setArg(3, xArrayCLBuffer_);
		kernel.setArg(4, gridXZCLBuffer_);
		kernel.setArg(5, delayTensorCLBuffer_);

		const std::size_t rowGroupSize = 32;
		const std::size_t colGroupSize = 1;
		const std::size_t elemGroupSize = 1;
		const std::size_t globalN0 = OCLUtil::roundUpToMultipleOfGroupSize(numRows, rowGroupSize);
		const std::size_t globalN1 = OCLUtil::roundUpToMultipleOfGroupSize(numCols, colGroupSize);
		const std::size_t globalN2 = OCLUtil::roundUpToMultipleOfGroupSize(config_.numElements, elemGroupSize);

		cl::Event kernelEvent;

		clCommandQueue_.enqueueNDRangeKernel(
			kernel,
			cl::NullRange, // offset
			cl::NDRange(globalN0, globalN1, globalN2), // global
			cl::NDRange(rowGroupSize, colGroupSize, elemGroupSize), // local
			nullptr /* previous events */, &kernelEvent);

#ifdef LAB_ENABLE_EXECUTION_TIME_MEASUREMENT
		kernelEvent.wait();
		tCalculateDelaysML.put(calculateDelaysTimer.getTime());
#endif
	} catch (cl::Error& e) {
		THROW_EXCEPTION(OCLException, "[calculateDelaysSTAKernel] " << e);
	}

	//==================================================
	// Delay and sum.
	//==================================================

#ifdef LAB_ENABLE_EXECUTION_TIME_MEASUREMENT
	Timer delaySumTimer;
#endif
	cl::Kernel procKernel;
	try {
		if (coherenceFactor_.enabled()) {
			procKernel = cl::Kernel(clProgram_, "processRowColumnSTAPCFKernel");
		} else {
			procKernel = cl::Kernel(clProgram_, "processRowColumnSTAKernel");
		}
		procKernel.setArg( 0, numCols);
		procKernel.setArg( 1, numRows);
		procKernel.setArg( 2, signalOffset_);
		procKernel.setArg( 3, signalTensorCLBuffer_);
		procKernel.setArg( 4, config_.numElements);
		procKernel.setArg( 5, signalLength_);
		procKernel.setArg( 6, config_.firstTxElem);
		procKernel.setArg( 7, config_.lastTxElem);
		procKernel.setArg( 8, rxApodCLBuffer_);
		procKernel.setArg( 9, delayTensorCLBuffer_);
		if (coherenceFactor_.enabled()) {
			std::vector<TFloat> cfConstants;
			coherenceFactor_.implementation().getConstants(cfConstants);
			procKernel.setArg(10, cfConstants[2] /* factor */);
			procKernel.setArg(11, gridValueCLBuffer_);
			procKernel.setArg(12, gridFactorCLBuffer_);
		} else {
			procKernel.setArg(10, gridValueCLBuffer_);
		}
	} catch (cl::Error& e) {
		THROW_EXCEPTION(OCLException, "[Kernel preparation] " << e);
	}

	cl::Event procKernelEvent;
	try {
		std::size_t rowGroupSize;
		std::size_t colGroupSize;
		// Adjusted for GTX-1660.
		if (coherenceFactor_.enabled()) {
			rowGroupSize = 4;
			colGroupSize = 16;
		} else {
			rowGroupSize = 16;
			colGroupSize = 16;
		}
		const std::size_t globalN0 = OCLUtil::roundUpToMultipleOfGroupSize(numRows, rowGroupSize);
		const std::size_t globalN1 = OCLUtil::roundUpToMultipleOfGroupSize(numCols, colGroupSize);

		clCommandQueue_.enqueueNDRangeKernel(
			procKernel,
			cl::NullRange, // offset
			cl::NDRange(globalN0, globalN1), // global
			cl::NDRange(rowGroupSize, colGroupSize), // local
			nullptr /* previous events */, &procKernelEvent);
	} catch (cl::Error& e) {
		THROW_EXCEPTION(OCLException, "[processRowColumnSTA*Kernel] " << e);
	}

#ifdef LAB_ENABLE_EXECUTION_TIME_MEASUREMENT
	procKernelEvent.wait();
	tDelaySumML.put(delaySumTimer.getTime());
#endif

	clCommandQueue_.enqueueReadBuffer(
		gridValueCLBuffer_, CL_BLOCKING, 0 /* offset */,
		gridValueHostMem_->sizeInBytes, gridValueHostMem_->hostPtr);
	if (coherenceFactor_.enabled()) {
		clCommandQueue_.enqueueReadBuffer(
			gridFactorCLBuffer_, CL_BLOCKING, 0 /* offset */,
			gridFactorHostMem_->sizeInBytes, gridFactorHostMem_->hostPtr);
	}

	const unsigned int numSignals = numTxElem * config_.numElements;
	const TFloat valueCoef = 1.0f / TFloat(numSignals);

	//==================================================
	// Read the formed image.
	//==================================================
	for (unsigned int col = 0; col < numCols; ++col) {
		unsigned int gridPointIdx = col * numRows;
		for (unsigned int row = 0; row < numRows; ++row, ++gridPointIdx) {
			gridData(col, row).value = gridValueHostMem_->hostPtr[gridPointIdx] * valueCoef;
		}
	}
	if (coherenceFactor_.enabled()) {
		for (unsigned int col = 0; col < numCols; ++col) {
			unsigned int gridPointIdx = col * numRows;
			for (unsigned int row = 0; row < numRows; ++row, ++gridPointIdx) {
				gridData(col, row).factor = gridFactorHostMem_->hostPtr[gridPointIdx];
			}
		}
	}

	//LOG_DEBUG << "END ========== VectorialSTAOCLProcessor::process ==========";
}

template<typename TFloat>
std::string
VectorialSTAOCLProcessor<TFloat>::getKernels()
{
	return R"CLC(

__kernel
void
calculateDelaysSTAKernel(
		unsigned int numCols,
		unsigned int numRows,
		MFloat invCT,
		__global const MFloat* xArray,
		__global const MFloat (*gridXZ)[2],
		__global MFloat* delayTensor)
{
	const unsigned int row = get_global_id(0);
	if (row >= numRows) return;

	const unsigned int col = get_global_id(1);
	if (col >= numCols) return;

	const unsigned int elem = get_global_id(2);
	if (elem >= NUM_RX_ELEM) return;

	__global const MFloat (*point)[2] = gridXZ + col * numRows + row;
	delayTensor[((elem * numCols) + col) * numRows + row] =
			distance2DY0(xArray[elem], (*point)[0], (*point)[1]) * invCT;
}

__kernel
void
processRowColumnSTAKernel(
		unsigned int numCols,
		unsigned int numRows,
		MFloat signalOffset,
		__global MFloat (*signalTensor)[2],
		unsigned int signalTensorN2,
		unsigned int signalTensorN3,
		unsigned int firstTxElem,
		unsigned int lastTxElem,
		__constant MFloat* rxApod,
		__global MFloat* delayTensor,
		__global MFloat* gridValue)
{
	const unsigned int signalLength = signalTensorN3;
	const unsigned int maxPosition = signalLength - 2;

	const unsigned int row = get_global_id(0);
	if (row >= numRows) return;

	const unsigned int col = get_global_id(1);
	if (col >= numCols) return;

	MFloat sumRe = 0;
	MFloat sumIm = 0;
	for (int txElem = firstTxElem; txElem <= lastTxElem; ++txElem) {
		const MFloat txDelay = delayTensor[(txElem * numCols + col) * numRows + row];
		const MFloat txOffset = signalOffset + txDelay;
		for (int rxElem = 0; rxElem < NUM_RX_ELEM; ++rxElem) {
			__global const MFloat (*p)[2] = signalTensor + (((txElem - firstTxElem) * signalTensorN2) + rxElem) * signalLength;
			// Linear interpolation.
			const MFloat position = txOffset + delayTensor[(rxElem * numCols + col) * numRows + row];
			if (position >= 0.0f) {
				const unsigned int positionIdx = (unsigned int) position;
				if (positionIdx <= maxPosition) {
					const MFloat k = position - positionIdx;
					MFloat v0[2]; // complex
					v0[0] = p[positionIdx][0];
					v0[1] = p[positionIdx][1];
					MFloat v1[2]; // complex
					v1[0] = p[positionIdx + 1][0];
					v1[1] = p[positionIdx + 1][1];
					MFloat v[2]; // complex
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

__kernel
void
processRowColumnSTAPCFKernel(
		unsigned int numCols,
		unsigned int numRows,
		MFloat signalOffset,
		__global MFloat (*signalTensor)[2],
		unsigned int signalTensorN2,
		unsigned int signalTensorN3,
		unsigned int firstTxElem,
		unsigned int lastTxElem,
		__constant MFloat* rxApod,
		__global MFloat* delayTensor,
		MFloat pcfFactor,
		__global MFloat* gridValue,
		__global MFloat* gridFactor)
{
	MFloat rxSignalListRe[NUM_RX_ELEM];
	MFloat rxSignalListIm[NUM_RX_ELEM];

	const unsigned int signalLength = signalTensorN3;
	const unsigned int maxPosition = signalLength - 2;

	const unsigned int row = get_global_id(0);
	if (row >= numRows) return;

	const unsigned int col = get_global_id(1);
	if (col >= numCols) return;

	for (int i = 0; i < NUM_RX_ELEM; ++i) {
		rxSignalListRe[i] = 0;
		rxSignalListIm[i] = 0;
	}
	for (int txElem = firstTxElem; txElem <= lastTxElem; ++txElem) {
		const MFloat txDelay = delayTensor[(txElem * numCols + col) * numRows + row];
		const MFloat txOffset = signalOffset + txDelay;
		for (int rxElem = 0; rxElem < NUM_RX_ELEM; ++rxElem) {
			__global const MFloat (*p)[2] = signalTensor + (((txElem - firstTxElem) * signalTensorN2) + rxElem) * signalLength;
			// Linear interpolation.
			const MFloat position = txOffset + delayTensor[(rxElem * numCols + col) * numRows + row];
			if (position >= 0.0f) {
				const unsigned int positionIdx = (unsigned int) position;
				if (positionIdx <= maxPosition) {
					const MFloat k = position - positionIdx;
					MFloat v0[2]; // complex
					v0[0] = p[positionIdx][0];
					v0[1] = p[positionIdx][1];
					MFloat v1[2]; // complex
					v1[0] = p[positionIdx + 1][0];
					v1[1] = p[positionIdx + 1][1];
					MFloat v[2]; // complex
					v[0] = v0[0] + k * (v1[0] - v0[0]);
					v[1] = v0[1] + k * (v1[1] - v0[1]);

					rxSignalListRe[rxElem] += v[0] * rxApod[rxElem];
					rxSignalListIm[rxElem] += v[1] * rxApod[rxElem];
				}
			}
		}
	}

	const MFloat pcf = calcPCF(rxSignalListRe, rxSignalListIm, pcfFactor);
	MFloat sumRe = 0;
	MFloat sumIm = 0;
	for (int i = 0; i < NUM_RX_ELEM; ++i) {
		sumRe += rxSignalListRe[i];
		sumIm += rxSignalListIm[i];
	}

	const unsigned int point = col * numRows + row;
	gridValue[point] = sqrt(sumRe * sumRe + sumIm * sumIm);
	gridFactor[point] = pcf;
}

)CLC";
}

} // namespace Lab

#endif // VECTORIAL_STA_OCL_PROCESSOR_H
