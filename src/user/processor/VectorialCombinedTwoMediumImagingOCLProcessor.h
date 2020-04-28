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

#ifndef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_OCL_PROCESSOR_H
#define VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_OCL_PROCESSOR_H

#include <algorithm> /* copy, fill, min */
#include <cmath> /* abs, ceil, sqrt */
#include <complex>
#include <cstddef> /* std::size_t */
#include <cstring> /* memset */
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <type_traits> /* is_same */
#include <utility> /* make_pair */
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
#include "TemplateUtil.h"
#include "Tensor3.h"
#include "Timer.h"
#include "TwoMediumSTAConfiguration.h"
#include "Util.h"
#include "XZ.h"

#include <CL/cl2.hpp>
#include "OCLCoherenceFactor.h"
#include "OCLUtil.h"



// Depends on the signal.
// 1.0 --> pi radian / sample at the original sampling rate.
#define VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_OCL_PROCESSOR_UPSAMP_FILTER_HALF_TRANSITION_WIDTH (0.2)

#define VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_OCL_PROCESSOR_NUM_RAW_DATA_BUFFERS 2
#define VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_OCL_PROCESSOR_USE_TRANSPOSE 1



namespace Lab {

// Two-medium image formation, using analytic signals (each sample is a real-imag vector).
// The final image is a combination of sub-images, using apodization.
//
// Without PCF, this class is slower than VectorialCombinedTwoMediumImagingProcessor,
// using a Core i5-3470. With PCF, this class is faster.
//
// Processing steps:
//   Find row at the interface                               - CPU
//   Calculate delays at the interface                       - CPU
//   Calculate delays in the grid above the interface        - CPU
//   Signal preparation                                      - CPU
//   Apply delays and store the sample                       - CPU
//   Use apodization, [apply PCF] and accumulate the samples - OpenCL
//
// The grid must be rectangular.
//
// Tested only with numElements = 32.
//
template<typename TFloat>
class VectorialCombinedTwoMediumImagingOCLProcessor {
public:
	struct StepConfiguration {
		unsigned int baseElemIdx;
		unsigned int baseElem;
		unsigned int txElem;
	};

	VectorialCombinedTwoMediumImagingOCLProcessor(
			const TwoMediumSTAConfiguration<TFloat>& config,
			std::vector<Matrix<TFloat>>& acqDataList,
			unsigned int upsamplingFactor,
			AnalyticSignalCoherenceFactorProcessor<TFloat>& coherenceFactor,
			TFloat maxFermatBlockSize,
			TFloat peakOffset,
			unsigned int signalStartOffset);
	~VectorialCombinedTwoMediumImagingOCLProcessor() = default;

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
	enum {
		OCL_TRANSPOSE_GROUP_SIZE_DIM_0 = 16, // do not change
		OCL_WORK_ITEMS_PER_GROUP = 64
	};

	struct CalculateDelays;

	struct PrepareDataWithOneTxElem;
	struct PrepareDataThreadData {
		Interpolator<TFloat> interpolator;
		HilbertEnvelope<TFloat> envelope;
		std::vector<TFloat, tbb::cache_aligned_allocator<TFloat>> signal;
	};

	struct ProcessColumnWithOneTxElem;

	VectorialCombinedTwoMediumImagingOCLProcessor(const VectorialCombinedTwoMediumImagingOCLProcessor&) = delete;
	VectorialCombinedTwoMediumImagingOCLProcessor& operator=(const VectorialCombinedTwoMediumImagingOCLProcessor&) = delete;
	VectorialCombinedTwoMediumImagingOCLProcessor(VectorialCombinedTwoMediumImagingOCLProcessor&&) = delete;
	VectorialCombinedTwoMediumImagingOCLProcessor& operator=(VectorialCombinedTwoMediumImagingOCLProcessor&&) = delete;

	static std::string getKernels();

	const TwoMediumSTAConfiguration<TFloat>& config_;
	std::vector<Matrix<TFloat>>& acqDataList_;
	unsigned int upsamplingFactor_;
	AnalyticSignalCoherenceFactorProcessor<TFloat>& coherenceFactor_;
	TFloat maxFermatBlockSize_;
	const TFloat lambda2_;
	std::size_t signalLength_;
	Tensor3<std::complex<TFloat>, tbb::cache_aligned_allocator<std::complex<TFloat>>> signalTensor_;
	TFloat signalOffset_;
	std::vector<unsigned int, tbb::cache_aligned_allocator<unsigned int>> minRowIdx_; // for each column
	std::vector<unsigned int, tbb::cache_aligned_allocator<unsigned int>> firstGridPointIdx_; // for each column
	std::vector<TFloat, tbb::cache_aligned_allocator<TFloat>> xArray_;
	Matrix<TFloat, tbb::cache_aligned_allocator<TFloat>> medium1DelayMatrix_;
	Tensor3<TFloat, tbb::cache_aligned_allocator<TFloat>> delayTensor_;
	unsigned int rawDataN1_;
	unsigned int rawDataN2_;
	std::size_t rawDataSizeInBytes_;
	std::unique_ptr<tbb::enumerable_thread_specific<PrepareDataThreadData>> prepareDataTLS_;

	unsigned int rxApodSize_;
	unsigned int gridXZN1_;
	unsigned int gridXZN2_;

	// OpenCL.
	bool clDataInitialized_;
	cl::Context clContext_;
	cl::Program clProgram_;
	cl::CommandQueue clCommandQueue_;
	std::vector<OCLPinnedHostMem<TFloat>> rawDataHostMemList_;
	cl::Buffer rawDataCLBuffer_;
#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_OCL_PROCESSOR_USE_TRANSPOSE
	cl::Buffer rawDataTCLBuffer_;
#endif
	cl::Buffer gridValueReCLBuffer_;
	cl::Buffer gridValueImCLBuffer_;
	std::vector<OCLPinnedHostMem<TFloat>> gridValueHostMemList_;
	cl::Buffer rxApodCLBuffer_;
};



template<typename TFloat>
VectorialCombinedTwoMediumImagingOCLProcessor<TFloat>::VectorialCombinedTwoMediumImagingOCLProcessor(
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
		, rawDataN1_()
		, rawDataN2_()
		, rawDataSizeInBytes_()
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
		prepareDataThreadData.interpolator.prepare(upsamplingFactor_, VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_OCL_PROCESSOR_UPSAMP_FILTER_HALF_TRANSITION_WIDTH);
	}
	prepareDataThreadData.signal.resize(signalLength_);
	prepareDataTLS_ = std::make_unique<tbb::enumerable_thread_specific<PrepareDataThreadData>>(prepareDataThreadData);

	clContext_ = OCLUtil::initOpenCL();

	std::vector<std::string> kernelStrings = {OCLCoherenceFactor::code() + getKernels()};
	clProgram_ = cl::Program(clContext_, kernelStrings);
	std::ostringstream progOpt;
	progOpt << OCL_PROGRAM_BUILD_OPTIONS
			<< " -DNUM_RX_ELEM="        << config.numElements
			<< " -DSTATISTICS_N="       << config.numElements
			<< " -DCOHERENCE_FACTOR_N=" << config.numElements
			<< " -DMFloat="             << TemplateUtil::typeName<TFloat>()
			<< " -DGROUP_SIZE="         << OCL_TRANSPOSE_GROUP_SIZE_DIM_0;
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
VectorialCombinedTwoMediumImagingOCLProcessor<TFloat>::process(
							const std::vector<StepConfiguration>& stepConfigList,
							const std::vector<XZ<TFloat>>& interfacePointList,
							const std::vector<TFloat>& rxApod,
							const Matrix<XZ<TFloat>>& gridXZ,
							Matrix<std::complex<TFloat>>& gridValue)
{
	//LOG_DEBUG << "BEGIN ========== VectorialCombinedTwoMediumImagingOCLProcessor::process ==========";

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
	LOG_DEBUG << "numCols: " << numCols << " numRows: " << numRows;

	minRowIdx_.resize(numCols);
	firstGridPointIdx_.resize(numCols);
	delayTensor_.resize(numCols, numRows, config_.numElementsMux);
	signalTensor_.resize(stepConfigList.size(), config_.numElements, signalLength_);
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
	for (unsigned int col = 0; col < numCols; ++col) {
		auto& point = gridXZ(col, 0);
		firstGridPointIdx_[col] = gridPointIdx; // the points below the interface are not considered

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
	const std::size_t numGridPoints = gridPointIdx;
#ifdef USE_EXECUTION_TIME_MEASUREMENT
	tMinRowIdxML.put(minRowIdxTimer.getTime());
#endif

#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_OCL_PROCESSOR_USE_TRANSPOSE
	const std::size_t transpNumGridPoints = OCLUtil::roundUpToMultipleOfGroupSize(numGridPoints, OCL_TRANSPOSE_GROUP_SIZE_DIM_0);
	LOG_DEBUG << "transpNumGridPoints: " << transpNumGridPoints;
	rawDataN1_ = transpNumGridPoints;
	rawDataN2_ = 2 * config_.numElements /* real, imag */;
#else
	rawDataN1_ = 2 * config_.numElements /* real, imag */;
	rawDataN2_ = numGridPoints;
#endif
	rawDataSizeInBytes_ = rawDataN1_ * rawDataN2_ * sizeof(TFloat);

	if (!clDataInitialized_) {
		rxApodSize_ = rxApod.size();
		gridXZN1_   = gridXZ.n1();
		gridXZN2_   = gridXZ.n2();
#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_OCL_PROCESSOR_USE_TRANSPOSE
		const std::size_t reservedGridPoints = OCLUtil::roundUpToMultipleOfGroupSize(numCols * numRows, OCL_TRANSPOSE_GROUP_SIZE_DIM_0);
#else
		const std::size_t reservedGridPoints = numCols * numRows;
#endif
		const std::size_t reservedRawDataSize = reservedGridPoints * 2 * config_.numElements;
		const std::size_t reservedRawDataSizeInBytes = reservedRawDataSize * sizeof(TFloat);
		LOG_DEBUG << "reservedGridPoints: " << reservedGridPoints <<
				" reservedRawDataSize: " << reservedRawDataSize <<
				" reservedRawDataSizeInBytes: " << reservedRawDataSizeInBytes;

		for (unsigned int i = 0; i < VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_OCL_PROCESSOR_NUM_RAW_DATA_BUFFERS; ++i) {
			rawDataHostMemList_.emplace_back(clContext_, clCommandQueue_, reservedRawDataSize, CL_MAP_WRITE_INVALIDATE_REGION);
		}
		rawDataCLBuffer_     = cl::Buffer(clContext_, CL_MEM_READ_ONLY , reservedRawDataSizeInBytes);
#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_OCL_PROCESSOR_USE_TRANSPOSE
		rawDataTCLBuffer_    = cl::Buffer(clContext_, CL_MEM_READ_WRITE, reservedRawDataSizeInBytes);
#endif
		gridValueHostMemList_.emplace_back(clContext_, clCommandQueue_, numCols * numRows, CL_MAP_READ); // real
		gridValueHostMemList_.emplace_back(clContext_, clCommandQueue_, numCols * numRows, CL_MAP_READ); // imag
		gridValueReCLBuffer_ = cl::Buffer(clContext_, CL_MEM_READ_WRITE, numCols * numRows * sizeof(TFloat));
		gridValueImCLBuffer_ = cl::Buffer(clContext_, CL_MEM_READ_WRITE, numCols * numRows * sizeof(TFloat));
		rxApodCLBuffer_      = cl::Buffer(clContext_, CL_MEM_READ_ONLY , Util::sizeInBytes(rxApod));

		clDataInitialized_ = true;
	} else {
		if (rxApodSize_ != rxApod.size() ||
		    gridXZN1_   != gridXZ.n1()   ||
		    gridXZN2_   != gridXZ.n2()) {
			THROW_EXCEPTION(InvalidParameterException, "Data size mismatch.");
		}
	}

	clCommandQueue_.enqueueFillBuffer(
		gridValueReCLBuffer_, (TFloat) 0, 0 /* offset */, numGridPoints * sizeof(TFloat));
	clCommandQueue_.enqueueFillBuffer(
		gridValueImCLBuffer_, (TFloat) 0, 0 /* offset */, numGridPoints * sizeof(TFloat));
	clCommandQueue_.enqueueWriteBuffer(
		rxApodCLBuffer_, CL_BLOCKING, 0 /* offset */,
		Util::sizeInBytes(rxApod), rxApod.data());

	const TFloat c2ByC1 = config_.propagationSpeed2 / config_.propagationSpeed1;

	if (xArray_.empty()) {
		ArrayGeometry::getElementX2D(config_.numElementsMux, config_.pitch, TFloat(0) /* offset */, xArray_);
	}

#ifdef USE_EXECUTION_TIME_MEASUREMENT
	Timer medium1DelayMatrixTimer;
#endif
	for (unsigned int elem = 0; elem < config_.numElementsMux; ++elem) {
		TFloat* delays = &medium1DelayMatrix_(elem, 0);
		for (unsigned int i = 0; i < interfacePointList.size(); ++i) {
			const XZ<TFloat>& ifPoint = interfacePointList[i];
			delays[i] = Geometry::distance2DY0(xArray_[elem], ifPoint.x, ifPoint.z) * c2ByC1;
		}
	}
#ifdef USE_EXECUTION_TIME_MEASUREMENT
	tMedium1DelayMatrixML.put(medium1DelayMatrixTimer.getTime());

	Timer calculateDelaysTimer;
#endif
	CalculateDelays calculateDelaysOp = {
		numRows,
		config_,
		config_.samplingFrequency * upsamplingFactor_,
		config_.samplingFrequency * upsamplingFactor_ / config_.propagationSpeed2,
		fermatBlockSize,
		interfacePointList,
		xArray_,
		minRowIdx_,
		medium1DelayMatrix_,
		gridXZ,
		delayTensor_
	};
	//tbb::parallel_for(tbb::blocked_range<unsigned int>(0, numCols), calculateDelaysOp);
	tbb::parallel_for(tbb::blocked_range<unsigned int>(0, numCols, 1 /* grain size */), calculateDelaysOp, tbb::simple_partitioner());
	//calculateDelaysOp(tbb::blocked_range<unsigned int>(0, numCols)); // single-thread
#ifdef USE_EXECUTION_TIME_MEASUREMENT
	tCalculateDelaysML.put(calculateDelaysTimer.getTime());
#endif

	// Only one transmit element.
#ifdef USE_EXECUTION_TIME_MEASUREMENT
	Timer prepareDataTimer;
#endif
	unsigned int stepIdx = 0;
	for (const auto& stepConfig : stepConfigList) {
		PrepareDataWithOneTxElem prepareDataOp = {
			samplesPerChannelLow,
			acqDataList_,
			upsamplingFactor_,
			stepIdx,
			stepConfig.baseElemIdx,
			*prepareDataTLS_,
			signalTensor_
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

#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_OCL_PROCESSOR_USE_TRANSPOSE
	cl::Kernel transpKernel;
#endif
	cl::Kernel procImageKernel;

	//==================================================
	// [OpenCL] Kernel preparation.
	//==================================================
	try {
#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_OCL_PROCESSOR_USE_TRANSPOSE
		transpKernel = cl::Kernel(clProgram_, "transposeKernel");
		transpKernel.setArg(0, rawDataCLBuffer_);
		transpKernel.setArg(1, rawDataTCLBuffer_);
		transpKernel.setArg(2, rawDataN2_);
		transpKernel.setArg(3, rawDataN1_);
		transpKernel.setArg(4, cl::Local(OCL_TRANSPOSE_GROUP_SIZE_DIM_0 * (OCL_TRANSPOSE_GROUP_SIZE_DIM_0 + 1) * sizeof(TFloat)));
#endif
		procImageKernel = cl::Kernel(clProgram_, coherenceFactor_.enabled() ? "processImagePCFKernel" : "processImageKernel");
#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_OCL_PROCESSOR_USE_TRANSPOSE
		procImageKernel.setArg(0, rawDataTCLBuffer_);
		procImageKernel.setArg(1, rawDataN1_);
#else
		procImageKernel.setArg(0, rawDataCLBuffer_);
		procImageKernel.setArg(1, rawDataN2_);
#endif
		procImageKernel.setArg(2, gridValueReCLBuffer_);
		procImageKernel.setArg(3, gridValueImCLBuffer_);
		procImageKernel.setArg(4, rxApodCLBuffer_);
		if (coherenceFactor_.enabled()) {
			std::vector<TFloat> cfConstants;
			coherenceFactor_.implementation().getConstants(cfConstants);
			procImageKernel.setArg(5, cfConstants[2] /* factor */);
		}
	} catch (cl::Error& e) {
		THROW_EXCEPTION(OCLException, "[Kernel preparation] OpenCL error: " << e.what() << " (" << e.err() << ").");
	}

	std::vector<cl::Event> writeBufferEventList(VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_OCL_PROCESSOR_NUM_RAW_DATA_BUFFERS);
	std::size_t procImageKernelGlobalSize = OCLUtil::roundUpToMultipleOfGroupSize(numGridPoints, OCL_WORK_ITEMS_PER_GROUP);
	//LOG_DEBUG << numGridPoints << ':' << procImageKernelGlobalSize << ':' << OCL_WORK_ITEMS_PER_GROUP;

	//==================================================
	// Step configuration loop.
	//==================================================
	for (unsigned int i = 0; i < stepConfigList.size(); ++i) {
		const unsigned int rawBufferIdx = i % VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_OCL_PROCESSOR_NUM_RAW_DATA_BUFFERS;
		const auto& stepConfig = stepConfigList[i];
		LOG_DEBUG << "stepConfig.baseElemIdx: " << stepConfig.baseElemIdx << " rawBufferIdx: " << rawBufferIdx;

		if (writeBufferEventList[rawBufferIdx]()) writeBufferEventList[rawBufferIdx].wait();

		Timer delayStoreTimer;

		//==================================================
		// Delay and store.
		//==================================================
		ProcessColumnWithOneTxElem processColumnOp = {
			static_cast<unsigned int>(numRows),
			0,
			config_,
			signalOffset_,
			signalTensor_,
			stepConfig,
			minRowIdx_,
			firstGridPointIdx_,
			delayTensor_,
			rawDataHostMemList_[rawBufferIdx].hostPtr,
			rawDataN2_,
		};

		//tbb::parallel_for(tbb::blocked_range<unsigned int>(0, numCols), processColumnOp);
		tbb::parallel_for(tbb::blocked_range<unsigned int>(0, numCols, 1 /* grain size */), processColumnOp, tbb::simple_partitioner());
		//processColumnOp(tbb::blocked_range<unsigned int>(0, numCols)); // single-thread

		LOG_DEBUG << "TIME OCL DELAY-STORE " << delayStoreTimer.getTime();

		try {
			Timer transfTimer;

			//==================================================
			// [OpenCL] Memory transfer to device.
			//==================================================
			clCommandQueue_.enqueueWriteBuffer(
				rawDataCLBuffer_, CL_NON_BLOCKING, 0 /* offset */,
				rawDataSizeInBytes_, rawDataHostMemList_[rawBufferIdx].hostPtr,
				nullptr /* previous events */, &writeBufferEventList[rawBufferIdx]);

			LOG_DEBUG << "TIME OCL TRANSF " << transfTimer.getTime(); // useful only with CL_BLOCKING
		} catch (cl::Error& e) {
			THROW_EXCEPTION(OCLException, "[Memory transfer] OpenCL error: " << e.what() << " (" << e.err() << ").");
		}

		try {
			cl::Event kernelEvent;

#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_OCL_PROCESSOR_USE_TRANSPOSE
			Timer transpTimer;
			//==================================================
			// [OpenCL] Transpose kernel.
			//==================================================
			clCommandQueue_.enqueueNDRangeKernel(
				transpKernel,
				cl::NullRange, // offset
				cl::NDRange(rawDataN2_, rawDataN1_), // global
				cl::NDRange(OCL_TRANSPOSE_GROUP_SIZE_DIM_0, OCL_TRANSPOSE_GROUP_SIZE_DIM_0), // local
				nullptr /* previous events */, &kernelEvent);
			//kernelEvent.wait();
			LOG_DEBUG << "TIME OCL TRANSPOSE " << transpTimer.getTime(); // useful only with kernelEvent.wait()
#endif
			Timer procTimer;
			//==================================================
			// [OpenCL] Final processing kernel.
			//==================================================
			clCommandQueue_.enqueueNDRangeKernel(
				procImageKernel,
				cl::NullRange, // offset
				cl::NDRange(procImageKernelGlobalSize), // global
				cl::NDRange(OCL_WORK_ITEMS_PER_GROUP), // local
				nullptr /* previous events */, &kernelEvent);
			//kernelEvent.wait();
			LOG_DEBUG << "TIME OCL PROC " << procTimer.getTime(); // useful only with kernelEvent.wait()

		} catch (cl::Error& e) {
			THROW_EXCEPTION(OCLException, "[Imaging] OpenCL error: " << e.what() << " (" << e.err() << ").");
		}
	}

	//==================================================
	// [OpenCL] Read the formed image.
	//==================================================
	clCommandQueue_.enqueueReadBuffer(
		gridValueReCLBuffer_, CL_NON_BLOCKING, 0 /* offset */,
		numGridPoints * sizeof(TFloat), gridValueHostMemList_[0].hostPtr);
	clCommandQueue_.enqueueReadBuffer(
		gridValueImCLBuffer_, CL_BLOCKING, 0 /* offset */,
		numGridPoints * sizeof(TFloat), gridValueHostMemList_[1].hostPtr);
	for (unsigned int col = 0; col < numCols; ++col) {
		for (unsigned int row = 0; row < minRowIdx_[col]; ++row) {
			gridValue(col, row) = 0;
		}
		unsigned int gridPointIdx = firstGridPointIdx_[col];
		for (unsigned int row = minRowIdx_[col]; row < numRows; ++row, ++gridPointIdx) {
			gridValue(col, row) = std::complex<TFloat>(
						gridValueHostMemList_[0].hostPtr[gridPointIdx],
						gridValueHostMemList_[1].hostPtr[gridPointIdx]);
		}
	}

#ifdef USE_EXECUTION_TIME_MEASUREMENT
	tProcessColumnML.put(processColumnTimer.getTime());
#endif

	//LOG_DEBUG << "END ========== VectorialCombinedTwoMediumImagingOCLProcessor::process ==========";
}



template<typename TFloat>
struct VectorialCombinedTwoMediumImagingOCLProcessor<TFloat>::CalculateDelays {
	void operator()(const tbb::blocked_range<unsigned int>& r) const {
		//LOG_DEBUG << "col = " << r.begin() << " n = " << (r.end() - r.begin());

		for (unsigned int col = r.begin(); col != r.end(); ++col) {
			if (minRowIdx[col] >= numRows) continue;

			for (unsigned int elem = 0; elem < config.numElementsMux; ++elem) {
				unsigned int lastInterfaceIdx = 0;

				// The first row above the interface.
				{
					const auto& point = gridXZ(col, minRowIdx[col]);

					// Fermat's principle. Find the fastest path.
					TFloat tMin;
					unsigned int idxMin;
					FermatPrinciple::findMinTimeInTwoSteps(
							fermatBlockSize,
							config.propagationSpeed1, config.propagationSpeed2,
							interfacePointList,
							xArray[elem], TFloat(0), point.x, point.z,
							tMin, idxMin);
					delayTensor(col, minRowIdx[col], elem) = tMin * fs;
					lastInterfaceIdx = idxMin;
				}

				const TFloat* medium1Delays = &medium1DelayMatrix(elem, 0);

				for (unsigned int row = minRowIdx[col] + 1; row < numRows; ++row) {
					const auto& point = gridXZ(col, row);
					unsigned int idxMin = lastInterfaceIdx;
					TFloat tC2Min;
					{
						const XZ<TFloat>& ifPoint = interfacePointList[idxMin];
						tC2Min = medium1Delays[idxMin] + Geometry::distance2D(ifPoint.x, ifPoint.z, point.x, point.z);
					}
					for (unsigned int idxSearch = idxMin + 1, end = interfacePointList.size(); idxSearch < end; ++idxSearch) {
						const XZ<TFloat>& ifPoint = interfacePointList[idxSearch];
						const TFloat tC2 = medium1Delays[idxSearch] + Geometry::distance2D(ifPoint.x, ifPoint.z, point.x, point.z);
						if (tC2 >= tC2Min) {
							break;
						} else {
							tC2Min = tC2;
							idxMin = idxSearch;
						}
					}
					if (idxMin == lastInterfaceIdx) { // if the previous search was not successful
						for (int idxSearch = static_cast<int>(idxMin) - 1; idxSearch >= 0; --idxSearch) { // if idxMin = 0, idxSearch will start with -1
							const XZ<TFloat>& ifPoint = interfacePointList[idxSearch];
							const TFloat tC2 = medium1Delays[idxSearch] + Geometry::distance2D(ifPoint.x, ifPoint.z, point.x, point.z);
							if (tC2 >= tC2Min) {
								break;
							} else {
								tC2Min = tC2;
								idxMin = idxSearch;
							}
						}
					}

//					unsigned int diff = (idxMin > lastInterfaceIdx) ?
//								idxMin - lastInterfaceIdx :
//								lastInterfaceIdx - idxMin;
//					if (diff > 1) {
//						LOG_DEBUG << "########## DIFF " << diff << " idxMin: " << idxMin << " col: " << col << " row - minRowIdx[col]: " << row - minRowIdx[col];
//					}

					delayTensor(col, row, elem) = tC2Min * fsInvC2;
					lastInterfaceIdx = idxMin;
				}

			}
		}
	}

	// Use only small types or references/pointers, because this object may be copied many times.
	const std::size_t numRows;
	const TwoMediumSTAConfiguration<TFloat>& config;
	const TFloat fs;
	const TFloat fsInvC2;
	const unsigned int fermatBlockSize;
	const std::vector<XZ<TFloat>>& interfacePointList;
	const std::vector<TFloat, tbb::cache_aligned_allocator<TFloat>>& xArray;
	const std::vector<unsigned int, tbb::cache_aligned_allocator<unsigned int>>& minRowIdx;
	const Matrix<TFloat, tbb::cache_aligned_allocator<TFloat>>& medium1DelayMatrix;
	const Matrix<XZ<TFloat>>& gridXZ;
	Tensor3<TFloat, tbb::cache_aligned_allocator<TFloat>>& delayTensor;
};



template<typename TFloat>
struct VectorialCombinedTwoMediumImagingOCLProcessor<TFloat>::PrepareDataWithOneTxElem {
	void operator()(const tbb::blocked_range<unsigned int>& r) const {
		PrepareDataThreadData& local = prepareDataTLS.local();

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
					&signalTensor(stepIdx, rxElem, 0));
		}
	}

	// Use only small types or references/pointers, because this object may be copied many times.
	const std::size_t samplesPerChannelLow;
	const std::vector<Matrix<TFloat>>& acqDataList;
	const unsigned int upsamplingFactor;
	const unsigned int stepIdx;
	const unsigned int baseElementIdx;
	tbb::enumerable_thread_specific<PrepareDataThreadData>& prepareDataTLS;
	Tensor3<std::complex<TFloat>, tbb::cache_aligned_allocator<std::complex<TFloat>>>& signalTensor;
};



template<typename TFloat>
struct VectorialCombinedTwoMediumImagingOCLProcessor<TFloat>::ProcessColumnWithOneTxElem {
	void operator()(const tbb::blocked_range<unsigned int>& r) const {
		//LOG_DEBUG << "col = " << r.begin() << " n = " << (r.end() - r.begin());

		const unsigned int signalLength = signalTensor.n3();
		const unsigned int maxPosition = signalLength - 2;

		for (unsigned int col = r.begin(); col != r.end(); ++col) {
			unsigned int gridPointIdx = firstGridPointIdx[col] - firstGridPointIdx[firstCol];
			for (unsigned int row = minRowIdx[col]; row < numRows; ++row, ++gridPointIdx) {
				const TFloat* delays = &delayTensor(col, row, stepConfig.baseElem);
				const TFloat txDelay = delays[stepConfig.txElem];
				const TFloat txOffset = signalOffset + txDelay;
				const auto* p = &signalTensor(stepConfig.baseElemIdx, 0 /* rxElem */, 0);
				for (unsigned int rxElem = 0; rxElem < config.numElements; ++rxElem, p += signalLength) {
					const unsigned int rxIdx = rxElem * 2;

					// Linear interpolation.
					const TFloat position = txOffset + delays[rxElem];
					if (position >= 0.0f) {
						const unsigned int positionIdx = static_cast<unsigned int>(position);
						if (positionIdx <= maxPosition) {
							const TFloat k = position - positionIdx;
							const auto v0 = p[positionIdx];
							const auto v1 = p[positionIdx + 1];
							const std::complex<TFloat> v = v0 + k * (v1 - v0);
#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_OCL_PROCESSOR_USE_TRANSPOSE
							rawData[gridPointIdx * rawDataN2 + rxIdx    ] = v.real();
							rawData[gridPointIdx * rawDataN2 + rxIdx + 1] = v.imag();
#else
							rawData[ rxIdx      * rawDataN2 + gridPointIdx] = v.real();
							rawData[(rxIdx + 1) * rawDataN2 + gridPointIdx] = v.imag();
#endif
							continue;
						}
					}
#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_OCL_PROCESSOR_USE_TRANSPOSE
					rawData[gridPointIdx * rawDataN2 + rxIdx    ] = 0;
					rawData[gridPointIdx * rawDataN2 + rxIdx + 1] = 0;
#else
					rawData[ rxIdx      * rawDataN2 + gridPointIdx] = 0;
					rawData[(rxIdx + 1) * rawDataN2 + gridPointIdx] = 0;
#endif
				}
			}
		}
	}

	// Use only small types or references/pointers, because this object may be copied many times.
	const unsigned int numRows;
	const unsigned int firstCol;
	const TwoMediumSTAConfiguration<TFloat>& config;
	const TFloat signalOffset;
	const Tensor3<std::complex<TFloat>, tbb::cache_aligned_allocator<std::complex<TFloat>>>& signalTensor;
	const StepConfiguration stepConfig;
	const std::vector<unsigned int, tbb::cache_aligned_allocator<unsigned int>>& minRowIdx;
	const std::vector<unsigned int, tbb::cache_aligned_allocator<unsigned int>>& firstGridPointIdx;
	const Tensor3<TFloat, tbb::cache_aligned_allocator<TFloat>>& delayTensor;
	TFloat* rawData;
	unsigned int rawDataN2;
};



template<typename TFloat>
std::string
VectorialCombinedTwoMediumImagingOCLProcessor<TFloat>::getKernels()
{
	return R"CLC(

// NVIDIA sm_50 or newer:
//   - Local (shared) memory has 32 banks of 32 bits.
__kernel
void
transposeKernel(
		__global MFloat* rawData,
		__global MFloat* rawDataT,
		unsigned int oldSizeX,
		unsigned int oldSizeY,
		__local MFloat* temp) // GROUP_SIZE * (GROUP_SIZE + 1) -- +1 to avoid bank conflicts
{
	unsigned int iX = get_global_id(0);
	unsigned int iY = get_global_id(1);
	if (iX < oldSizeX && iY < oldSizeY) {
		temp[get_local_id(0) * (GROUP_SIZE + 1) + get_local_id(1)] = rawData[iY * oldSizeX + iX];
	}

	barrier(CLK_LOCAL_MEM_FENCE);

	iX = get_group_id(1) * GROUP_SIZE + get_local_id(0);
	iY = get_group_id(0) * GROUP_SIZE + get_local_id(1);
	if (iX < oldSizeY && iY < oldSizeX) {
		rawDataT[iY * oldSizeY + iX] = temp[get_local_id(1) * (GROUP_SIZE + 1) + get_local_id(0)];
	}
}

__kernel
void
processImageKernel(
		__global MFloat* rawData,
		unsigned int numGridPoints,
		__global MFloat* gridValueRe,
		__global MFloat* gridValueIm,
		__constant MFloat* rxApod)
{
	MFloat rxSignalListRe[NUM_RX_ELEM];
	MFloat rxSignalListIm[NUM_RX_ELEM];

	const unsigned int point = get_global_id(0);
	if (point >= numGridPoints) return;

	unsigned int idx1 = point;
	unsigned int idx2 = point + numGridPoints;
	const unsigned int step = 2 * numGridPoints;
	for (int i = 0; i < NUM_RX_ELEM; ++i, idx1 += step, idx2 += step) {
		rxSignalListRe[i] = rawData[idx1];
		rxSignalListIm[i] = rawData[idx2];
	}

	MFloat sumRe = 0;
	MFloat sumIm = 0;
	for (int i = 0; i < NUM_RX_ELEM; ++i) {
		sumRe += rxSignalListRe[i] * rxApod[i];
		sumIm += rxSignalListIm[i] * rxApod[i];
	}

	gridValueRe[point] += sumRe;
	gridValueIm[point] += sumIm;
}

__kernel
void
processImagePCFKernel(
		__global MFloat* rawData,
		unsigned int numGridPoints,
		__global MFloat* gridValueRe,
		__global MFloat* gridValueIm,
		__constant MFloat* rxApod,
		MFloat pcfFactor)
{
	MFloat rxSignalListRe[NUM_RX_ELEM];
	MFloat rxSignalListIm[NUM_RX_ELEM];

	const unsigned int point = get_global_id(0);
	if (point >= numGridPoints) return;

	unsigned int idx1 = point;
	unsigned int idx2 = point + numGridPoints;
	const unsigned int step = 2 * numGridPoints;
	for (int i = 0; i < NUM_RX_ELEM; ++i, idx1 += step, idx2 += step) {
		rxSignalListRe[i] = rawData[idx1];
		rxSignalListIm[i] = rawData[idx2];
	}

	const MFloat pcf = calcPCF(rxSignalListRe, rxSignalListIm, pcfFactor);

	MFloat sumRe = 0;
	MFloat sumIm = 0;
	for (int i = 0; i < NUM_RX_ELEM; ++i) {
		sumRe += rxSignalListRe[i] * rxApod[i];
		sumIm += rxSignalListIm[i] * rxApod[i];
	}

	gridValueRe[point] += sumRe * pcf;
	gridValueIm[point] += sumIm * pcf;
}

)CLC";
}

} // namespace Lab

#endif // VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_OCL_PROCESSOR_H
