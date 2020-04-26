/***************************************************************************
 *  Copyright 2014, 2017, 2018, 2020 Marcelo Y. Matuda                     *
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

#ifndef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_OCL_PROCESSOR_2_H
#define VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_OCL_PROCESSOR_2_H

#include <algorithm> /* copy */
#include <cmath> /* ceil, sqrt */
#include <complex>
#include <cstddef> /* std::size_t */
#include <memory>
#include <string>
#include <vector>

#include <tbb/cache_aligned_allocator.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/partitioner.h>
#include <tbb/tbb.h>

#include "ArrayGeometry.h"
#include "CoherenceFactor.h"
#include "Exception.h"
#include "ExecutionTimeMeasurement.h"
#include "FermatPrinciple.h"
#include "Geometry.h"
#include "HilbertEnvelope.h"
#include "Interpolator.h"
#include "Log.h"
#include "Matrix.h"
#include "Timer.h"
#include "TwoMediumSTAConfiguration.h"
#include "Util.h"
#include "XZ.h"

#include <CL/cl2.hpp>
#include "OCLCoherenceFactor.h"
#include "OCLGeometry.h"
#include "OCLUtil.h"



// Depends on the signal.
// 1.0 --> pi radian / sample at the original sampling rate.
#define VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_OCL_PROCESSOR_2_UPSAMP_FILTER_HALF_TRANSITION_WIDTH (0.2)



namespace Lab {

// Two-medium image formation, using analytic signals (each sample is a real-imag vector).
// The final image is a combination of sub-images, using apodization.
//
// Processing steps:
//   Find row at the interface                                             - CPU
//   Calculate delays at the interface                                     - CPU
//   Calculate delays in the grid above the interface                      - OpenCL
//   Signal preparation                                                    - CPU
//   Apply delays, use apodization, [apply PCF] and accumulate the samples - OpenCL
//
// The grid must be rectangular.
//
template<typename TFloat>
class VectorialCombinedTwoMediumImagingOCLProcessor2 {
public:
	struct StepConfiguration {
		unsigned int baseElemIdx;
		unsigned int baseElem;
		unsigned int txElem;
	};

	VectorialCombinedTwoMediumImagingOCLProcessor2(
			const TwoMediumSTAConfiguration<TFloat>& config,
			std::vector<Matrix<TFloat>>& acqDataList,
			unsigned int upsamplingFactor,
			AnalyticSignalCoherenceFactorProcessor<TFloat>& coherenceFactor,
			TFloat maxFermatBlockSize,
			TFloat peakOffset,
			unsigned int signalStartOffset);
	~VectorialCombinedTwoMediumImagingOCLProcessor2() = default;

	void process(
		const std::vector<StepConfiguration>& stepConfigList,
		const std::vector<XZ<TFloat>>& interfacePointList,
		const std::vector<TFloat>& rxApod,
		const Matrix<XZ<TFloat>>& gridXZ,
		Matrix<std::complex<TFloat>>& gridValue);

#ifdef USE_EXECUTION_TIME_MEASUREMENT
	MeasurementList<double> tMinRowIdxML;
	MeasurementList<double> tMedium1DelayMatrixML;
	MeasurementList<double> tCalculateDelaysML;
	MeasurementList<double> tPrepareDataML;
	MeasurementList<double> tProcessColumnML;
	void execTimeMeasReset() {
		tMinRowIdxML.reset(         EXECUTION_TIME_MEASUREMENT_ITERATIONS);
		tMedium1DelayMatrixML.reset(EXECUTION_TIME_MEASUREMENT_ITERATIONS);
		tCalculateDelaysML.reset(   EXECUTION_TIME_MEASUREMENT_ITERATIONS);
		tPrepareDataML.reset(       EXECUTION_TIME_MEASUREMENT_ITERATIONS);
		tProcessColumnML.reset(     EXECUTION_TIME_MEASUREMENT_ITERATIONS);
	}
	void execTimeMeasShowResults() {
		EXECUTION_TIME_MEASUREMENT_LOG_TIMES("tMinRowIdx:         ", tMinRowIdxML);
		EXECUTION_TIME_MEASUREMENT_LOG_TIMES("tMedium1DelayMatrix:", tMedium1DelayMatrixML);
		EXECUTION_TIME_MEASUREMENT_LOG_TIMES("tCalculateDelays:   ", tCalculateDelaysML);
		EXECUTION_TIME_MEASUREMENT_LOG_TIMES("tPrepareData:       ", tPrepareDataML);
		EXECUTION_TIME_MEASUREMENT_LOG_TIMES("tProcessColumn:     ", tProcessColumnML);
	}
#endif

private:
	struct PrepareDataThreadData {
		Interpolator<TFloat> interpolator;
		HilbertEnvelope<TFloat> envelope;
		std::vector<TFloat, tbb::cache_aligned_allocator<TFloat>> signal;
	};

	struct PrepareDataWithOneTxElem;

	VectorialCombinedTwoMediumImagingOCLProcessor2(const VectorialCombinedTwoMediumImagingOCLProcessor2&) = delete;
	VectorialCombinedTwoMediumImagingOCLProcessor2& operator=(const VectorialCombinedTwoMediumImagingOCLProcessor2&) = delete;
	VectorialCombinedTwoMediumImagingOCLProcessor2(VectorialCombinedTwoMediumImagingOCLProcessor2&&) = delete;
	VectorialCombinedTwoMediumImagingOCLProcessor2& operator=(VectorialCombinedTwoMediumImagingOCLProcessor2&&) = delete;

	std::string getKernels() const;

	const TwoMediumSTAConfiguration<TFloat>& config_;
	std::vector<Matrix<TFloat>>& acqDataList_;
	unsigned int upsamplingFactor_;
	AnalyticSignalCoherenceFactorProcessor<TFloat>& coherenceFactor_;
	TFloat maxFermatBlockSize_;
	const TFloat lambda2_;
	unsigned int signalLength_;
	TFloat signalOffset_;
	std::vector<unsigned int, tbb::cache_aligned_allocator<unsigned int>> minRowIdx_; // for each column
	std::vector<TFloat, tbb::cache_aligned_allocator<TFloat>> xArray_;
	Matrix<TFloat, tbb::cache_aligned_allocator<TFloat>> medium1DelayMatrix_;
	std::unique_ptr<tbb::enumerable_thread_specific<PrepareDataThreadData>> prepareDataTLS_;

	unsigned int stepConfigListSize_;
	unsigned int interfacePointListSize_;
	unsigned int rxApodSize_;
	unsigned int gridXZN1_;
	unsigned int gridXZN2_;

	// OpenCL.
	bool clDataInitialized_;
	cl::Context clContext_;
	cl::Program clProgram_;
	cl::CommandQueue clCommandQueue_;
	cl::Buffer gridXZCLBuffer_;
	cl::Buffer rxApodCLBuffer_;
	cl::Buffer xArrayCLBuffer_;
	std::unique_ptr<OCLPinnedHostMem<TFloat[2]>> signalTensorHostMem_;
	cl::Buffer signalTensorCLBuffer_;
	cl::Buffer minRowIdxCLBuffer_;
	cl::Buffer interfacePointListCLBuffer_;
	cl::Buffer medium1DelayMatrixCLBuffer_;
	cl::Buffer delayTensorCLBuffer_;
	std::unique_ptr<OCLPinnedHostMem<TFloat[2]>> gridValueHostMem_;
	cl::Buffer gridValueCLBuffer_;
};



template<typename TFloat>
struct VectorialCombinedTwoMediumImagingOCLProcessor2<TFloat>::PrepareDataWithOneTxElem {
	void operator()(const tbb::blocked_range<unsigned int>& r) const {
		auto& local = prepareDataTLS.local();

		for (unsigned int rxElem = r.begin(); rxElem != r.end(); ++rxElem) {
			if (upsamplingFactor > 1) {
				// Interpolate the signal.
				local.interpolator.interpolate(&acqDataList[baseElementIdx](rxElem, 0), samplesPerChannelLow, local.signal.data());
			} else {
				auto range = acqDataList[baseElementIdx].range2(rxElem);
				std::copy(range.begin(), range.end(), local.signal.begin());
			}

			Util::removeDC(local.signal.data(), local.signal.size());

			// Obtain the analytic signal.
			local.envelope.getAnalyticSignal(
					local.signal.data(),
					local.signal.size(),
					signalTensor + ((stepIdx * signalTensorN2) + rxElem) * signalTensorN3);
		}
	}

	// Use only small types or references/pointers, because this object may be copied many times.
	const std::size_t samplesPerChannelLow;
	const std::vector<Matrix<TFloat>>& acqDataList;
	const unsigned int upsamplingFactor;
	const unsigned int stepIdx;
	const unsigned int baseElementIdx;
	tbb::enumerable_thread_specific<PrepareDataThreadData>& prepareDataTLS;
	TFloat (*signalTensor)[2];
	const unsigned int signalTensorN2;
	const unsigned int signalTensorN3;
};

template<typename TFloat>
VectorialCombinedTwoMediumImagingOCLProcessor2<TFloat>::VectorialCombinedTwoMediumImagingOCLProcessor2(
			const TwoMediumSTAConfiguration<TFloat>& config,
			std::vector<Matrix<TFloat>>& acqDataList,
			unsigned int upsamplingFactor,
			AnalyticSignalCoherenceFactorProcessor<TFloat>& coherenceFactor,
			TFloat maxFermatBlockSize,
			TFloat peakOffset,
			unsigned int signalStartOffset)
		: config_(config)
		, acqDataList_(acqDataList)
		, upsamplingFactor_(upsamplingFactor)
		, coherenceFactor_(coherenceFactor)
		, maxFermatBlockSize_(maxFermatBlockSize)
		, lambda2_(config_.propagationSpeed2 / config_.centerFrequency)
		, stepConfigListSize_()
		, interfacePointListSize_()
		, rxApodSize_()
		, gridXZN1_()
		, gridXZN2_()
		, clDataInitialized_()
{
	const std::size_t origSignalLength = acqDataList_[0].n2();

	signalOffset_ = (config_.samplingFrequency * upsamplingFactor_) * (peakOffset / config_.centerFrequency) - signalStartOffset * upsamplingFactor_;
	signalLength_ = origSignalLength * upsamplingFactor_;
	LOG_DEBUG << "signalOffset_: " << signalOffset_ << " signalLength_: " << signalLength_;

	PrepareDataThreadData prepareDataThreadData;
	if (upsamplingFactor_ > 1) {
		prepareDataThreadData.interpolator.prepare(upsamplingFactor_, VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_OCL_PROCESSOR_2_UPSAMP_FILTER_HALF_TRANSITION_WIDTH);
	}
	prepareDataThreadData.signal.resize(signalLength_);
	prepareDataTLS_ = std::make_unique<tbb::enumerable_thread_specific<PrepareDataThreadData>>(prepareDataThreadData);

	clContext_ = OCLUtil::initOpenCL();

	std::vector<std::string> kernelStrings = {OCLCoherenceFactor::code() + OCLGeometry::code() + getKernels()};
	clProgram_ = cl::Program(clContext_, kernelStrings);
	std::ostringstream progOpt;
	progOpt << OCL_PROGRAM_BUILD_OPTIONS
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
VectorialCombinedTwoMediumImagingOCLProcessor2<TFloat>::process(
						const std::vector<StepConfiguration>& stepConfigList,
						const std::vector<XZ<TFloat>>& interfacePointList,
						const std::vector<TFloat>& rxApod,
						const Matrix<XZ<TFloat>>& gridXZ,
						Matrix<std::complex<TFloat>>& gridValue)
{
	//LOG_DEBUG << "BEGIN ========== VectorialCombinedTwoMediumImagingOCLProcessor2::process ==========";

	if (stepConfigList.empty()) {
		THROW_EXCEPTION(InvalidParameterException, "The list of step configurations is empty.");
	}
	if (gridXZ.n1() != gridValue.n1() || gridXZ.n2() != gridValue.n2()) {
		THROW_EXCEPTION(InvalidParameterException, "gridXZ and gridValue have different sizes.");
	}

	const std::size_t samplesPerChannelLow = acqDataList_[0].n2();

	const unsigned int numCols = gridXZ.n1();
	const unsigned int numRows = gridXZ.n2();
	if (numCols == 0) {
		THROW_EXCEPTION(InvalidValueException, "Zero columns in the grid.");
	}

	minRowIdx_.resize(numCols);
	medium1DelayMatrix_.resize(config_.numElementsMux, interfacePointList.size());

	XZ<TFloat> p1 = interfacePointList[0];
	XZ<TFloat> p2 = interfacePointList[1];
	const TFloat dx = p2.x - p1.x;
	const TFloat dz = p2.z - p1.z;
	const TFloat r = std::sqrt(dx * dx + dz * dz);
	const unsigned int fermatBlockSize = FermatPrinciple::calcBlockSizeForTwoStepSearch(interfacePointList.size(), r, lambda2_, maxFermatBlockSize_);
	LOG_DEBUG << "fermatBlockSize: " << fermatBlockSize;

#ifdef USE_EXECUTION_TIME_MEASUREMENT
	Timer minRowIdxTimer;
#endif
	//==================================================
	// Find the z-coordinate of the interface in
	// each column (minimum row).
	//==================================================
	const TFloat zStepGrid = gridXZ(0, 1).z - gridXZ(0, 0).z;
	unsigned int gridPointIdx = 0;
	for (unsigned int col = 0; col < gridXZ.n1(); ++col) {
		auto& point = gridXZ(col, 0);

		// Find the z coordinate of the interface.
		TFloat zIdxMin;
		unsigned int idxMin;
		FermatPrinciple::findNearestPointInXInTwoSteps(
				fermatBlockSize,
				interfacePointList,
				point.x,
				zIdxMin, idxMin);

		if (zIdxMin <= gridXZ(col, 0).z) {
			minRowIdx_[col] = 1;
		} else if (zIdxMin >= gridXZ(col, gridXZ.n2() - 1).z) {
			minRowIdx_[col] = gridXZ.n2(); // after the last index
		} else {
			minRowIdx_[col] = static_cast<unsigned int>(std::ceil((zIdxMin - gridXZ(col, 0).z) / zStepGrid)) + 1;
		}
		if (minRowIdx_[col] >= gridXZ.n2()) {
			THROW_EXCEPTION(InvalidValueException, "No valid rows in column " << col << '.');
		}

		gridPointIdx += gridXZ.n2() - minRowIdx_[col];
	}
	LOG_DEBUG << "number of valid grid points: " << gridPointIdx;
#ifdef USE_EXECUTION_TIME_MEASUREMENT
	tMinRowIdxML.put(minRowIdxTimer.getTime());
#endif

	if (xArray_.empty()) {
		ArrayGeometry::getElementX2D(config_.numElementsMux, config_.pitch, TFloat(0) /* offset */, xArray_);
	}

	if (!clDataInitialized_) {
		stepConfigListSize_     = stepConfigList.size();
		interfacePointListSize_ = interfacePointList.size();
		rxApodSize_             = rxApod.size();
		gridXZN1_               = gridXZ.n1();
		gridXZN2_               = gridXZ.n2();

		gridXZCLBuffer_             = cl::Buffer(clContext_, CL_MEM_READ_ONLY , sizeof(TFloat[2]) * gridXZN1_ * gridXZN2_);
		rxApodCLBuffer_             = cl::Buffer(clContext_, CL_MEM_READ_ONLY , sizeof(TFloat) * rxApodSize_);
		xArrayCLBuffer_             = cl::Buffer(clContext_, CL_MEM_READ_ONLY , sizeof(TFloat) * xArray_.size());
		signalTensorHostMem_        = std::make_unique<OCLPinnedHostMem<TFloat[2]>>(clContext_, clCommandQueue_,
									stepConfigListSize_ * config_.numElements * signalLength_,
									CL_MAP_WRITE_INVALIDATE_REGION);
		signalTensorCLBuffer_       = cl::Buffer(clContext_, CL_MEM_READ_ONLY , sizeof(TFloat[2]) * stepConfigListSize_ * config_.numElements * signalLength_);
		minRowIdxCLBuffer_          = cl::Buffer(clContext_, CL_MEM_READ_ONLY , sizeof(unsigned int) * minRowIdx_.size());
		interfacePointListCLBuffer_ = cl::Buffer(clContext_, CL_MEM_READ_ONLY , sizeof(TFloat[2]) * interfacePointListSize_);
		medium1DelayMatrixCLBuffer_ = cl::Buffer(clContext_, CL_MEM_READ_ONLY , sizeof(TFloat) * medium1DelayMatrix_.n1() * medium1DelayMatrix_.n2());
		delayTensorCLBuffer_        = cl::Buffer(clContext_, CL_MEM_READ_WRITE, sizeof(TFloat) * config_.numElementsMux * numCols * numRows);
		gridValueHostMem_           = std::make_unique<OCLPinnedHostMem<TFloat[2]>>(clContext_, clCommandQueue_,
									numCols * numRows,
									CL_MAP_READ);
		gridValueCLBuffer_          = cl::Buffer(clContext_, CL_MEM_READ_WRITE, sizeof(TFloat[2]) * numCols * numRows);

		clDataInitialized_ = true;
	} else {
		if (stepConfigListSize_     != stepConfigList.size()     ||
		    interfacePointListSize_ != interfacePointList.size() ||
		    rxApodSize_             != rxApod.size()             ||
		    gridXZN1_               != gridXZ.n1()               ||
		    gridXZN2_               != gridXZ.n2()) {
			THROW_EXCEPTION(InvalidParameterException, "Data size mismatch.");
		}
	}

#ifdef USE_EXECUTION_TIME_MEASUREMENT
	Timer medium1DelayMatrixTimer;
#endif
	const TFloat c2ByC1 = config_.propagationSpeed2 / config_.propagationSpeed1;
	for (unsigned int elem = 0; elem < config_.numElementsMux; ++elem) {
		TFloat* delays = &medium1DelayMatrix_(elem, 0);
		for (unsigned int i = 0; i < interfacePointList.size(); ++i) {
			const XZ<TFloat>& ifPoint = interfacePointList[i];
			delays[i] = Geometry::distance2DY0(xArray_[elem], ifPoint.x, ifPoint.z) * c2ByC1;
		}
	}
#ifdef USE_EXECUTION_TIME_MEASUREMENT
	tMedium1DelayMatrixML.put(medium1DelayMatrixTimer.getTime());
#endif

	// Prepare buffers.
	Timer prepareBuffersTimer;
	clCommandQueue_.enqueueWriteBuffer(
		gridXZCLBuffer_, CL_TRUE /* blocking */, 0 /* offset */,
		gridXZ.size() * sizeof(TFloat[2]), gridXZ.data());
	clCommandQueue_.enqueueWriteBuffer(
		rxApodCLBuffer_, CL_TRUE /* blocking */, 0 /* offset */,
		rxApod.size() * sizeof(TFloat), rxApod.data());
	clCommandQueue_.enqueueWriteBuffer(
		xArrayCLBuffer_, CL_TRUE /* blocking */, 0 /* offset */,
		xArray_.size() * sizeof(TFloat), xArray_.data());
	clCommandQueue_.enqueueWriteBuffer(
		minRowIdxCLBuffer_, CL_TRUE /* blocking */, 0 /* offset */,
		minRowIdx_.size() * sizeof(unsigned int), minRowIdx_.data());
	clCommandQueue_.enqueueWriteBuffer(
		interfacePointListCLBuffer_, CL_TRUE /* blocking */, 0 /* offset */,
		interfacePointList.size() * 2 * sizeof(TFloat), interfacePointList.data());
	clCommandQueue_.enqueueWriteBuffer(
		medium1DelayMatrixCLBuffer_, CL_TRUE /* blocking */, 0 /* offset */,
		medium1DelayMatrix_.size() * sizeof(TFloat), medium1DelayMatrix_.data());
	clCommandQueue_.enqueueFillBuffer(
		gridValueCLBuffer_, (TFloat) 0, 0 /* offset */, gridValueHostMem_->sizeInBytes);
	LOG_DEBUG << "PREPARE BUFFERS " << prepareBuffersTimer.getTime();

	try {
#ifdef USE_EXECUTION_TIME_MEASUREMENT
		Timer calculateDelaysTimer;
#endif
		cl::Kernel kernel(clProgram_, "calculateDelaysTwoMediumKernel");
		kernel.setArg( 0, numCols);
		kernel.setArg( 1, numRows);
		kernel.setArg( 2, config_.numElementsMux);
		kernel.setArg( 3, config_.samplingFrequency * upsamplingFactor_);
		kernel.setArg( 4, config_.samplingFrequency * upsamplingFactor_ / config_.propagationSpeed2);
		kernel.setArg( 5, config_.propagationSpeed1);
		kernel.setArg( 6, config_.propagationSpeed2);
		kernel.setArg( 7, fermatBlockSize);
		kernel.setArg( 8, interfacePointListCLBuffer_);
		kernel.setArg( 9, static_cast<unsigned int>(interfacePointList.size()));
		kernel.setArg(10, xArrayCLBuffer_);
		kernel.setArg(11, minRowIdxCLBuffer_);
		kernel.setArg(12, medium1DelayMatrixCLBuffer_);
		kernel.setArg(13, gridXZCLBuffer_);
		kernel.setArg(14, delayTensorCLBuffer_);

		// Adjusted for GTX-1660.
		const std::size_t colGroupSize  = 2;
		const std::size_t elemGroupSize = 8;
		const std::size_t globalN0 = OCLUtil::roundUpToMultipleOfGroupSize(numCols, colGroupSize);
		const std::size_t globalN1 = OCLUtil::roundUpToMultipleOfGroupSize(config_.numElementsMux, elemGroupSize);

		cl::Event kernelEvent;

		clCommandQueue_.enqueueNDRangeKernel(
			kernel,
			cl::NullRange, // offset
			cl::NDRange(globalN0, globalN1), // global
			cl::NDRange(colGroupSize, elemGroupSize), // local
			nullptr /* previous events */, &kernelEvent);

#ifdef USE_EXECUTION_TIME_MEASUREMENT
		//kernelEvent.wait();
		tCalculateDelaysML.put(calculateDelaysTimer.getTime());
#endif
	} catch (cl::Error& e) {
		THROW_EXCEPTION(OCLException, "[calculateDelaysTwoMediumKernel] OpenCL error: " << e.what() << " (" << e.err() << ").");
	}

#ifdef USE_EXECUTION_TIME_MEASUREMENT
	Timer prepareDataTimer;
#endif
	// Only one transmit element.
	unsigned int stepIdx = 0;
	for (const auto& stepConfig : stepConfigList) {
		PrepareDataWithOneTxElem prepareDataOp = {
			samplesPerChannelLow,
			acqDataList_,
			upsamplingFactor_,
			stepIdx,
			stepConfig.baseElemIdx,
			*prepareDataTLS_,
			signalTensorHostMem_->hostPtr,
			config_.numElements,
			signalLength_
		};
		//tbb::parallel_for(tbb::blocked_range<unsigned int>(0, config_.numElements), prepareDataOp);
		tbb::parallel_for(tbb::blocked_range<unsigned int>(0, config_.numElements, 1 /* grain size */), prepareDataOp, tbb::simple_partitioner());
		//prepareDataOp(tbb::blocked_range<unsigned int>(0, config_.numElements)); // single-thread

		++stepIdx;
	}
#ifdef USE_EXECUTION_TIME_MEASUREMENT
	tPrepareDataML.put(prepareDataTimer.getTime());

	Timer processColumnTimer;
#endif

	//Timer signalTransferTimer;

	clCommandQueue_.enqueueWriteBuffer(
		signalTensorCLBuffer_, CL_TRUE /* blocking */, 0 /* offset */,
		signalTensorHostMem_->sizeInBytes, signalTensorHostMem_->hostPtr);

	//LOG_DEBUG << "SIGNAL TRANSFER " << signalTransferTimer.getTime();

	//==================================================
	// Step configuration loop.
	//==================================================
	cl::Event procKernelEvent;
	for (unsigned int i = 0; i < stepConfigList.size(); ++i) {
		const auto& stepConfig = stepConfigList[i];
		LOG_DEBUG << "stepConfig.baseElemIdx: " << stepConfig.baseElemIdx;

		if (procKernelEvent()) procKernelEvent.wait();

		//==================================================
		// Delay and sum.
		//==================================================

		//Timer delaySumTimer;

		cl::Kernel procKernel;
		try {
			if (coherenceFactor_.enabled()) {
				procKernel = cl::Kernel(clProgram_, "processRowColumnWithOneTxElemPCFKernel");
			} else {
				procKernel = cl::Kernel(clProgram_, "processRowColumnWithOneTxElemKernel");
			}
			procKernel.setArg( 0, numCols);
			procKernel.setArg( 1, numRows);
			procKernel.setArg( 2, signalOffset_);
			procKernel.setArg( 3, signalTensorCLBuffer_);
			procKernel.setArg( 4, config_.numElements);
			procKernel.setArg( 5, signalLength_);
			procKernel.setArg( 6, stepConfig.baseElem);
			procKernel.setArg( 7, stepConfig.baseElemIdx);
			procKernel.setArg( 8, stepConfig.txElem);
			procKernel.setArg( 9, minRowIdxCLBuffer_);
			procKernel.setArg(10, delayTensorCLBuffer_);
			procKernel.setArg(11, gridValueCLBuffer_);
			procKernel.setArg(12, rxApodCLBuffer_);
			if (coherenceFactor_.enabled()) {
				std::vector<TFloat> cfConstants;
				coherenceFactor_.implementation().getConstants(cfConstants);
				procKernel.setArg(13, cfConstants[2] /* factor */);
			}
		} catch (cl::Error& e) {
			THROW_EXCEPTION(OCLException, "[Kernel preparation] OpenCL error: " << e.what() << " (" << e.err() << ").");
		}

		try {
			std::size_t rowGroupSize;
			std::size_t colGroupSize;
			// Adjusted for GTX-1660.
			if (coherenceFactor_.enabled()) {
				rowGroupSize = 16;
				colGroupSize = 16;
			} else {
				rowGroupSize = 64;
				colGroupSize = 2;
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
			THROW_EXCEPTION(OCLException, "[processRowColumnWithOneTxElem*Kernel] OpenCL error: " << e.what() << " (" << e.err() << ").");
		}

		//procKernelEvent.wait();
		//LOG_DEBUG << "DELAY-SUM " << delaySumTimer.getTime();
	}

	//procKernelEvent.wait();
	//Timer gridValueTransferTimer;

	clCommandQueue_.enqueueReadBuffer(
		gridValueCLBuffer_, CL_TRUE /* blocking */, 0 /* offset */,
		gridValueHostMem_->sizeInBytes, gridValueHostMem_->hostPtr);

	//==================================================
	// Read the formed image.
	//==================================================
	for (unsigned int col = 0; col < numCols; ++col) {
		for (unsigned int row = 0; row < minRowIdx_[col]; ++row) {
			gridValue(col, row) = 0;
		}
		unsigned int gridPointIdx = col * numRows + minRowIdx_[col];
		for (unsigned int row = minRowIdx_[col]; row < numRows; ++row, ++gridPointIdx) {
			gridValue(col, row) = std::complex<TFloat>(
							gridValueHostMem_->hostPtr[gridPointIdx][0],
							gridValueHostMem_->hostPtr[gridPointIdx][1]);
		}
	}

	//LOG_DEBUG << "GRID VALUE TRANSFER " << gridValueTransferTimer.getTime();

#ifdef USE_EXECUTION_TIME_MEASUREMENT
	tProcessColumnML.put(processColumnTimer.getTime());
#endif
	//LOG_DEBUG << "END ========== VectorialCombinedTwoMediumImagingOCLProcessor2::process ==========";
}

template<typename TFloat>
std::string
VectorialCombinedTwoMediumImagingOCLProcessor2<TFloat>::getKernels() const
{
	return R"CLC(

MFloat
calcTwoMediumTravelTime(MFloat x1, MFloat z1, MFloat xi, MFloat zi, MFloat x2, MFloat z2, MFloat invC1, MFloat invC2)
{
	const MFloat dx1 = xi - x1;
	const MFloat dz1 = zi - z1;
	const MFloat dx2 = x2 - xi;
	const MFloat dz2 = z2 - zi;
	return sqrt(dx1 * dx1 + dz1 * dz1) * invC1 + sqrt(dx2 * dx2 + dz2 * dz2) * invC2;
}

void
findMinTimeInTwoSteps(
		unsigned int blockSize,
		MFloat c1, MFloat c2,
		__global MFloat (*interfacePointList)[2],
		unsigned int interfacePointListSize,
		MFloat x1, MFloat z1, MFloat x2, MFloat z2,
		MFloat* tMin, unsigned int* idxMin)
{
	const MFloat invC1 = 1 / c1;
	const MFloat invC2 = 1 / c2;

	// First step: step = blockSize
	{
		__global const MFloat (*point)[2] = interfacePointList;
		const MFloat t = calcTwoMediumTravelTime(x1, z1, (*point)[0], (*point)[1], x2, z2, invC1, invC2);
		*tMin = t;
		*idxMin = 0;
	}
	for (unsigned int i = 1; i < interfacePointListSize; i += blockSize) {
		__global const MFloat (*point)[2] = interfacePointList + i;
		const MFloat t = calcTwoMediumTravelTime(x1, z1, (*point)[0], (*point)[1], x2, z2, invC1, invC2);
		if (t < *tMin) {
			*tMin = t;
			*idxMin = i;
		}
	}

	// Second step: step = 1
	const unsigned int iBegin = (*idxMin > (blockSize - 1U)) ? *idxMin - (blockSize - 1U) : 0;
	const unsigned int iEnd = min(*idxMin + blockSize, interfacePointListSize);
	for (unsigned int i = iBegin; i < iEnd; ++i) {
		__global const MFloat (*point)[2] = interfacePointList + i;
		const MFloat t = calcTwoMediumTravelTime(x1, z1, (*point)[0], (*point)[1], x2, z2, invC1, invC2);
		if (t < *tMin) {
			*tMin = t;
			*idxMin = i;
		}
	}
}

// This function is not efficient in GPUs, but it avoids the transfer of the large delayTensor.
__kernel
void
calculateDelaysTwoMediumKernel(
		unsigned int numCols,
		unsigned int numRows,
		unsigned int numElementsMux,
		MFloat fs,
		MFloat fsInvC2,
		MFloat c1,
		MFloat c2,
		unsigned int fermatBlockSize,
		__global MFloat (*interfacePointList)[2],
		unsigned int interfacePointListSize,
		__constant MFloat* xArray,
		__global unsigned int* minRowIdx,
		__global MFloat* medium1DelayMatrix,
		__global MFloat (*gridXZ)[2],
		__global MFloat* delayTensor)
{
	const unsigned int col = get_global_id(0);
	if (col >= numCols) return;

	const unsigned int elem = get_global_id(1);
	if (elem >= numElementsMux) return;

	if (minRowIdx[col] >= numRows) return;

	unsigned int lastInterfaceIdx = 0;

	// The first row above the interface.
	{
		__global const MFloat (*point)[2] = gridXZ + col * numRows + minRowIdx[col];

		// Fermat's principle. Find the fastest path.
		MFloat tMin;
		unsigned int idxMin;
		findMinTimeInTwoSteps(
				fermatBlockSize,
				c1, c2,
				interfacePointList,
				interfacePointListSize,
				xArray[elem], 0, (*point)[0], (*point)[1],
				&tMin, &idxMin);
		delayTensor[((elem * numCols) + col) * numRows + minRowIdx[col]] = tMin * fs;
		lastInterfaceIdx = idxMin;
	}

	__global const MFloat* medium1Delays = medium1DelayMatrix + elem * interfacePointListSize;

	for (unsigned int row = minRowIdx[col] + 1U; row < numRows; ++row) {
		__global const MFloat (*point)[2] = gridXZ + col * numRows + row;
		unsigned int idxMin = lastInterfaceIdx;
		MFloat tC2Min;
		{
			__global const MFloat (*ifPoint)[2] = interfacePointList + idxMin;
			tC2Min = medium1Delays[idxMin] + distance2D((*ifPoint)[0], (*ifPoint)[1], (*point)[0], (*point)[1]);
		}
		for (unsigned int idxSearch = idxMin + 1U; idxSearch < interfacePointListSize; ++idxSearch) {
			__global const MFloat (*ifPoint)[2] = interfacePointList + idxSearch;
			const MFloat tC2 = medium1Delays[idxSearch] + distance2D((*ifPoint)[0], (*ifPoint)[1], (*point)[0], (*point)[1]);
			if (tC2 >= tC2Min) {
				break;
			} else {
				tC2Min = tC2;
				idxMin = idxSearch;
			}
		}
		if (idxMin == lastInterfaceIdx) { // if the previous search was not successful
			for (int idxSearch = (int) idxMin - 1; idxSearch >= 0; --idxSearch) { // if idxMin = 0, idxSearch will start with -1
				__global const MFloat (*ifPoint)[2] = interfacePointList + idxSearch;
				const MFloat tC2 = medium1Delays[idxSearch] + distance2D((*ifPoint)[0], (*ifPoint)[1], (*point)[0], (*point)[1]);
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

__kernel
void
processRowColumnWithOneTxElemKernel(
		unsigned int numCols,
		unsigned int numRows,
		MFloat signalOffset,
		__global MFloat (*signalTensor)[2],
		unsigned int signalTensorN2,
		unsigned int signalTensorN3,
		unsigned int baseElem,
		unsigned int baseElemIdx,
		unsigned int txElem,
		__global unsigned int* minRowIdx,
		__global MFloat* delayTensor,
		__global MFloat (*gridValue)[2],
		__constant MFloat* rxApod)
{
	const unsigned int signalLength = signalTensorN3;
	const unsigned int maxPosition = signalLength - 2;

	const unsigned int col = get_global_id(1);
	if (col >= numCols) return;

	const unsigned int row = get_global_id(0);
	if (row < minRowIdx[col] || row >= numRows) return;

	const MFloat txDelay = delayTensor[((baseElem + txElem) * numCols + col) * numRows + row];
	const MFloat txOffset = signalOffset + txDelay;

	MFloat rxSignalSumRe = 0;
	MFloat rxSignalSumIm = 0;
	for (unsigned int rxElem = 0; rxElem < NUM_RX_ELEM; ++rxElem) {
		__global const MFloat (*p)[2] = signalTensor + baseElemIdx * signalTensorN2 * signalTensorN3
						+ rxElem * signalLength;
		// Linear interpolation.
		const MFloat position = txOffset + delayTensor[((baseElem + rxElem) * numCols + col) * numRows + row];
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

				rxSignalSumRe += v[0] * rxApod[rxElem];
				rxSignalSumIm += v[1] * rxApod[rxElem];
			}
		}
	}
	const unsigned int point = col * numRows + row;
	gridValue[point][0] += rxSignalSumRe;
	gridValue[point][1] += rxSignalSumIm;
}

__kernel
void
processRowColumnWithOneTxElemPCFKernel(
		unsigned int numCols,
		unsigned int numRows,
		MFloat signalOffset,
		__global MFloat (*signalTensor)[2],
		unsigned int signalTensorN2,
		unsigned int signalTensorN3,
		unsigned int baseElem,
		unsigned int baseElemIdx,
		unsigned int txElem,
		__global unsigned int* minRowIdx,
		__global MFloat* delayTensor,
		__global MFloat (*gridValue)[2],
		__constant MFloat* rxApod,
		MFloat pcfFactor)
{
	MFloat rxSignalListRe[NUM_RX_ELEM];
	MFloat rxSignalListIm[NUM_RX_ELEM];

	const unsigned int signalLength = signalTensorN3;
	const unsigned int maxPosition = signalLength - 2;

	const unsigned int col = get_global_id(1);
	if (col >= numCols) return;

	const unsigned int row = get_global_id(0);
	if (row < minRowIdx[col] || row >= numRows) return;

	const MFloat txDelay = delayTensor[((baseElem + txElem) * numCols + col) * numRows + row];
	const MFloat txOffset = signalOffset + txDelay;

	for (unsigned int rxElem = 0; rxElem < NUM_RX_ELEM; ++rxElem) {
		__global const MFloat (*p)[2] = signalTensor + baseElemIdx * signalTensorN2 * signalTensorN3
						+ rxElem * signalLength;
		// Linear interpolation.
		const MFloat position = txOffset + delayTensor[((baseElem + rxElem) * numCols + col) * numRows + row];
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

				rxSignalListRe[rxElem] = v[0];
				rxSignalListIm[rxElem] = v[1];
			}
		}
	}

	const MFloat pcf = calcPCF(rxSignalListRe, rxSignalListIm, pcfFactor);
	MFloat sumRe = 0;
	MFloat sumIm = 0;
	for (int i = 0; i < NUM_RX_ELEM; ++i) {
		sumRe += rxSignalListRe[i] * rxApod[i];
		sumIm += rxSignalListIm[i] * rxApod[i];
	}

	const unsigned int point = col * numRows + row;
	gridValue[point][0] += sumRe * pcf;
	gridValue[point][1] += sumIm * pcf;
}

)CLC";
}

} // namespace Lab

#endif // VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_OCL_PROCESSOR_2_H
