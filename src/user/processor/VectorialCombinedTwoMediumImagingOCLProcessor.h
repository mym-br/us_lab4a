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
#include <memory>
#include <sstream>
#include <string>
#include <utility> /* make_pair */
#include <vector>

#include <tbb/cache_aligned_allocator.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/partitioner.h>
#include <tbb/tbb.h>

#include <CL/cl2.hpp>

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
#include "Tensor3.h"
#include "Timer.h"
#include "TwoMediumSTAConfiguration.h"
#include "Util.h"
#include "XZ.h"



// Depends on the signal.
// 1.0 --> pi radian / sample at the original sampling rate.
#define VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_OCL_PROCESSOR_UPSAMP_FILTER_HALF_TRANSITION_WIDTH (0.2)

// 1: single buffer (slowest)
// 2: double buffer (fastest)
// 3: triple buffer (middle)
#define VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_OCL_PROCESSOR_NUM_RAW_DATA_BUFFERS 2

// Faster.
#define VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_OCL_PROCESSOR_USE_TRANSPOSE 1

#define VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_OCL_PROCESSOR_PROGRAM_BUILD_OPTIONS "-cl-std=CL1.2"



namespace Lab {

// The grid must be rectangular.
// Requirements:
// - Single precision
//
// Tested only with numElements = 32.
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
	~VectorialCombinedTwoMediumImagingOCLProcessor();

	void process(
		const std::vector<StepConfiguration>& stepConfigList,
		const std::vector<XZ<TFloat>>& interfacePointList,
		const std::vector<TFloat>& rxApod,
		const Matrix<XZ<TFloat>>& gridXZ,
		Matrix<std::complex<TFloat>>& gridValue);

#ifdef USE_EXECUTION_TIME_MEASUREMENT
	MeasurementList<double> tMinRowIdx;
	MeasurementList<double> tMedium1DelayMatrix;
	MeasurementList<double> tCalculateDelays;
	MeasurementList<double> tPrepareData;
	MeasurementList<double> tProcessColumn;
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

	struct ProcessColumn2ThreadData {
		AnalyticSignalCoherenceFactorProcessor<TFloat> coherenceFactor;
		std::vector<std::complex<TFloat>, tbb::cache_aligned_allocator<std::complex<TFloat>>> rxSignalSumList;
	};

	VectorialCombinedTwoMediumImagingOCLProcessor(const VectorialCombinedTwoMediumImagingOCLProcessor&) = delete;
	VectorialCombinedTwoMediumImagingOCLProcessor& operator=(const VectorialCombinedTwoMediumImagingOCLProcessor&) = delete;
	VectorialCombinedTwoMediumImagingOCLProcessor(VectorialCombinedTwoMediumImagingOCLProcessor&&) = delete;
	VectorialCombinedTwoMediumImagingOCLProcessor& operator=(VectorialCombinedTwoMediumImagingOCLProcessor&&) = delete;

	static std::size_t roundUpToMultipleOfGroupSize(std::size_t n, std::size_t groupSize) {
		std::size_t numGroups = (n + (groupSize - 1)) / groupSize;
		return numGroups * groupSize;
	}

	std::string getKernel() const;

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
	Matrix<TFloat, tbb::cache_aligned_allocator<TFloat>> medium1DelayMatrix_; // (interface_idx, element)
	Tensor3<TFloat, tbb::cache_aligned_allocator<TFloat>> delayTensor_;
	unsigned int rawDataN1_;
	unsigned int rawDataN2_;
	std::size_t rawDataSizeInBytes_;
	std::unique_ptr<tbb::enumerable_thread_specific<PrepareDataThreadData>> prepareDataTLS_;
	std::unique_ptr<tbb::enumerable_thread_specific<ProcessColumn2ThreadData>> processColumn2TLS_;
	bool clDataInitialized_;
	std::vector<TFloat, tbb::cache_aligned_allocator<TFloat>> gridValueRe_;
	std::vector<TFloat, tbb::cache_aligned_allocator<TFloat>> gridValueIm_;
	cl::Context clContext_;
	cl::Program clProgram_;
	cl::CommandQueue clCommandQueue_;
	std::vector<cl::Buffer> pinnedRawDataCLBufferList_;
	std::vector<TFloat*> mappedRawDataPtrList_;
	cl::Buffer rawDataCLBuffer_;
#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_OCL_PROCESSOR_USE_TRANSPOSE
	cl::Buffer rawDataTCLBuffer_;
#endif
	cl::Buffer gridValueReCLBuffer_;
	cl::Buffer gridValueImCLBuffer_;
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
		, clDataInitialized_()
		, pinnedRawDataCLBufferList_(VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_OCL_PROCESSOR_NUM_RAW_DATA_BUFFERS)
		, mappedRawDataPtrList_(VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_OCL_PROCESSOR_NUM_RAW_DATA_BUFFERS)
{
	if (sizeof(TFloat) != sizeof(float)) {
		THROW_EXCEPTION(InvalidParameterException, "Only single precision is supported.");
	}

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

	ProcessColumn2ThreadData processColumn2ThreadData;
	processColumn2ThreadData.coherenceFactor = coherenceFactor_;
	processColumn2TLS_ = std::make_unique<tbb::enumerable_thread_specific<ProcessColumn2ThreadData>>(processColumn2ThreadData);

	std::vector<cl::Platform> platforms;
	cl::Platform::get(&platforms);
	if (platforms.empty()) {
		THROW_EXCEPTION(UnavailableResourceException, "No OpenCL platforms available.");
	}

	if (Log::isDebugEnabled()) {
		for (auto& plat : platforms) {
			std::string name = plat.getInfo<CL_PLATFORM_NAME>();
			LOG_DEBUG << "OpenCL platform: " << name;

			std::vector<cl::Device> devices;
			plat.getDevices(CL_DEVICE_TYPE_ALL, &devices);
			if (devices.empty()) {
				THROW_EXCEPTION(UnavailableResourceException, "No OpenCL devices available for platform " << name << '.');
			}

			for (auto& dev : devices) {
				cl_device_type deviceType = dev.getInfo<CL_DEVICE_TYPE>();
				std::string devName = dev.getInfo<CL_DEVICE_NAME>();
				LOG_DEBUG << "  device name: " << devName;
				switch (deviceType) {
				case CL_DEVICE_TYPE_CPU:
					LOG_DEBUG << "    type: CPU";
					break;
				case CL_DEVICE_TYPE_GPU:
					LOG_DEBUG << "    type: GPU";
					break;
				default:
					LOG_DEBUG << "    type: other (" << deviceType << ").";
				}
			}
		}
	}

	if (OPENCL_PLATFORM >= platforms.size()) {
		THROW_EXCEPTION(UnavailableResourceException, "Invalid OpenCL platform: " << OPENCL_PLATFORM << '.');
	}
	cl::Platform chosenPlatform = platforms[OPENCL_PLATFORM];
	cl_context_properties contextProp[] = {CL_CONTEXT_PLATFORM, (cl_context_properties) chosenPlatform(), 0 /* end of list */};
	std::vector<cl::Device> devices;
	chosenPlatform.getDevices(CL_DEVICE_TYPE_ALL, &devices);
	if (devices.empty()) {
		THROW_EXCEPTION(UnavailableResourceException, "No OpenCL devices available.");
	}
	if (OPENCL_DEVICE >= devices.size()) {
		THROW_EXCEPTION(UnavailableResourceException, "Invalid OpenCL device: " << OPENCL_DEVICE << '.');
	}
	cl::Device chosenDevice = devices[OPENCL_DEVICE];
	clContext_ = cl::Context(chosenDevice, contextProp);
	std::vector<std::string> kernelStrings = {getKernel()};
	clProgram_ = cl::Program(clContext_, kernelStrings);
	std::ostringstream progOpt;
	progOpt << VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_OCL_PROCESSOR_PROGRAM_BUILD_OPTIONS <<
			" -DNUM_RX_ELEM=" << config.numElements;
	try {
		clProgram_.build(progOpt.str().c_str());
	} catch (...) {
		std::ostringstream msg;
		msg << "Error during OpenCL kernel compilation:\n";
		auto buildInfo = clProgram_.getBuildInfo<CL_PROGRAM_BUILD_LOG>();
		for (auto& pair : buildInfo) {
			msg << pair.second << "\n";
		}
		THROW_EXCEPTION(Exception, msg.str());
	}

	clCommandQueue_ = cl::CommandQueue(clContext_);
}

template<typename TFloat>
VectorialCombinedTwoMediumImagingOCLProcessor<TFloat>::~VectorialCombinedTwoMediumImagingOCLProcessor()
{
	if (clCommandQueue_()) {
		LOG_DEBUG << "~VectorialCombinedTwoMediumImagingOCLProcessor: enqueueUnmapMemObject";
		try {
			for (unsigned int i = 0; i < pinnedRawDataCLBufferList_.size(); ++i) {
				if (pinnedRawDataCLBufferList_[i]() && mappedRawDataPtrList_[i]) {
					cl::Event event;
					clCommandQueue_.enqueueUnmapMemObject(
									pinnedRawDataCLBufferList_[i],
									mappedRawDataPtrList_[i],
									nullptr, &event);
					event.wait();
				}
			}
		} catch (cl::Error& e) {
			LOG_ERROR << "[~VectorialCombinedTwoMediumImagingOCLProcessor: Unmap mappedRawDataPtrList_] Error: " << e.what();
		} catch (...) {
			LOG_ERROR << "[~VectorialCombinedTwoMediumImagingOCLProcessor: Unmap mappedRawDataPtrList_] Caught an unknown exception.";
		}
	}
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

	minRowIdx_.resize(gridXZ.n1() /* number of columns */);
	firstGridPointIdx_.resize(gridXZ.n1() /* number of columns */);
	delayTensor_.resize(gridXZ.n1() /* number of columns */, gridXZ.n2() /* number of rows */, config_.numElementsMux);
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
	for (unsigned int col = 0; col < gridXZ.n1(); ++col) {
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
#ifdef USE_EXECUTION_TIME_MEASUREMENT
	tMinRowIdx.put(minRowIdxTimer.getTime());
#endif

	const unsigned int cols = gridXZ.n1();
	if (cols == 0) {
		THROW_EXCEPTION(InvalidValueException, "Zero columns in the grid.");
	}
	std::size_t pointSum = 0;
	for (unsigned int col = 0; col < cols; ++col) {
		pointSum += gridXZ.n2() - minRowIdx_[col];
	}
	const std::size_t numGridPoints = pointSum;
	LOG_DEBUG << "cols: " << cols << " numGridPoints: " << numGridPoints;

#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_OCL_PROCESSOR_USE_TRANSPOSE
	const std::size_t transpNumGridPoints = roundUpToMultipleOfGroupSize(numGridPoints, OCL_TRANSPOSE_GROUP_SIZE_DIM_0);
	LOG_DEBUG << "numGridPoints: " << numGridPoints << " transpNumGridPoints: " << transpNumGridPoints;
	rawDataN1_ = transpNumGridPoints;
	rawDataN2_ = 2 * config_.numElements /* real, imag */;
#else
	rawDataN1_ = 2 * config_.numElements /* real, imag */;
	rawDataN2_ = numGridPoints;
#endif
	rawDataSizeInBytes_ = rawDataN1_ * rawDataN2_ * sizeof(TFloat);

	gridValueRe_.resize(numGridPoints);
	gridValueIm_.resize(numGridPoints);

	if (!clDataInitialized_) {
		for (unsigned int i = 0; i < pinnedRawDataCLBufferList_.size(); ++i) {
			pinnedRawDataCLBufferList_[i] = cl::Buffer(clContext_, CL_MEM_ALLOC_HOST_PTR, rawDataSizeInBytes_);
			mappedRawDataPtrList_[i] = static_cast<TFloat*>(
							clCommandQueue_.enqueueMapBuffer(
								pinnedRawDataCLBufferList_[i], CL_TRUE /* blocking */, CL_MAP_WRITE,
								0 /* offset */, rawDataSizeInBytes_));
		}
		rawDataCLBuffer_     = cl::Buffer(clContext_, CL_MEM_READ_ONLY , rawDataSizeInBytes_);
#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_OCL_PROCESSOR_USE_TRANSPOSE
		rawDataTCLBuffer_    = cl::Buffer(clContext_, CL_MEM_READ_WRITE, rawDataSizeInBytes_);
#endif
		gridValueReCLBuffer_ = cl::Buffer(clContext_, CL_MEM_READ_WRITE, Util::sizeInBytes(gridValueRe_));
		gridValueImCLBuffer_ = cl::Buffer(clContext_, CL_MEM_READ_WRITE, Util::sizeInBytes(gridValueIm_));
		rxApodCLBuffer_ =      cl::Buffer(clContext_, CL_MEM_READ_ONLY , Util::sizeInBytes(rxApod));

		clDataInitialized_ = true;
	}

	for (unsigned int i = 0; i < pinnedRawDataCLBufferList_.size(); ++i) {
		std::memset(mappedRawDataPtrList_[i], 0, rawDataSizeInBytes_);
	}
	std::memset(gridValueRe_.data(), 0, Util::sizeInBytes(gridValueRe_));
	std::memset(gridValueIm_.data(), 0, Util::sizeInBytes(gridValueIm_));
	clCommandQueue_.enqueueWriteBuffer(
		gridValueReCLBuffer_, CL_TRUE /* blocking */, 0 /* offset */,
		Util::sizeInBytes(gridValueRe_), gridValueRe_.data());
	clCommandQueue_.enqueueWriteBuffer(
		gridValueImCLBuffer_, CL_TRUE /* blocking */, 0 /* offset */,
		Util::sizeInBytes(gridValueIm_), gridValueIm_.data());
	clCommandQueue_.enqueueWriteBuffer(
		rxApodCLBuffer_, CL_TRUE /* blocking */, 0 /* offset */,
		Util::sizeInBytes(rxApod), rxApod.data());

	const TFloat c2ByC1 = config_.propagationSpeed2 / config_.propagationSpeed1;

	ArrayGeometry::getElementX2D(config_.numElementsMux, config_.pitch, TFloat(0) /* offset */, xArray_);

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
	tMedium1DelayMatrix.put(medium1DelayMatrixTimer.getTime());

	Timer calculateDelaysTimer;
#endif
	CalculateDelays calculateDelaysOp = {
		gridXZ.n2(),
		config_,
		config_.samplingFrequency * upsamplingFactor_,
		config_.samplingFrequency * upsamplingFactor_ / config_.propagationSpeed2,
		1 / config_.propagationSpeed1,
		1 / config_.propagationSpeed2,
		fermatBlockSize,
		interfacePointList,
		xArray_,
		minRowIdx_,
		medium1DelayMatrix_,
		gridXZ,
		delayTensor_
	};
	//tbb::parallel_for(tbb::blocked_range<unsigned int>(0, gridXZ.n1()), calculateDelaysOp);
	tbb::parallel_for(tbb::blocked_range<unsigned int>(0, gridXZ.n1(), 1 /* grain size */), calculateDelaysOp, tbb::simple_partitioner());
	//calculateDelaysOp(tbb::blocked_range<unsigned int>(0, gridXZ.n1())); // single-thread
#ifdef USE_EXECUTION_TIME_MEASUREMENT
	tCalculateDelays.put(calculateDelaysTimer.getTime());
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
	tPrepareData.put(prepareDataTimer.getTime());

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
		transpKernel.setArg(4, cl::Local(OCL_TRANSPOSE_GROUP_SIZE_DIM_0 * OCL_TRANSPOSE_GROUP_SIZE_DIM_0 * sizeof(TFloat)));
#endif
		procImageKernel = cl::Kernel(clProgram_, coherenceFactor_.enabled() ? "processImagePCFKernel" : "processImageKernel");
#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_OCL_PROCESSOR_USE_TRANSPOSE
		procImageKernel.setArg(0, rawDataTCLBuffer_);
#else
		procImageKernel.setArg(0, rawDataCLBuffer_);
#endif
		procImageKernel.setArg(1, gridValueReCLBuffer_);
		procImageKernel.setArg(2, gridValueImCLBuffer_);
		procImageKernel.setArg(3, rxApodCLBuffer_);
#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_OCL_PROCESSOR_USE_TRANSPOSE
		procImageKernel.setArg(4, rawDataN1_);
#else
		procImageKernel.setArg(4, rawDataN2_);
#endif
		if (coherenceFactor_.enabled()) {
			std::vector<TFloat> cfConstants;
			coherenceFactor_.implementation().getConstants(cfConstants);
			procImageKernel.setArg(5, cfConstants[2] /* factor */);
		}
	} catch (cl::Error& e) {
		LOG_ERROR << "[Kernel preparation] OpenCL error: " << e.what() << " (" << e.err() << ").";
		throw;
	}

	std::vector<cl::Event> writeBufferEventList(pinnedRawDataCLBufferList_.size());
	std::size_t procImageKernelGlobalSize = roundUpToMultipleOfGroupSize(numGridPoints, OCL_WORK_ITEMS_PER_GROUP);
	//LOG_DEBUG << numGridPoints << ':' << procImageKernelGlobalSize << ':' << OCL_WORK_ITEMS_PER_GROUP;

	//==================================================
	// Step configuration loop.
	//==================================================
	for (unsigned int i = 0; i < stepConfigList.size(); ++i) {
		const unsigned int rawBufferIdx = i % pinnedRawDataCLBufferList_.size();
		const auto& stepConfig = stepConfigList[i];
		LOG_DEBUG << "stepConfig.baseElemIdx: " << stepConfig.baseElemIdx << " rawBufferIdx: " << rawBufferIdx;

		if (writeBufferEventList[rawBufferIdx]()) writeBufferEventList[rawBufferIdx].wait();

		Timer delayStoreTimer;

		//==================================================
		// Delay and store.
		//==================================================
		ProcessColumnWithOneTxElem processColumnOp = {
			static_cast<unsigned int>(gridXZ.n2()),
			0,
			config_,
			signalOffset_,
			signalTensor_,
			stepConfig,
			minRowIdx_,
			firstGridPointIdx_,
			delayTensor_,
			mappedRawDataPtrList_[rawBufferIdx],
			rawDataN2_,
		};

		//tbb::parallel_for(tbb::blocked_range<unsigned int>(0, cols), processColumnOp);
		tbb::parallel_for(tbb::blocked_range<unsigned int>(0, cols, 1 /* grain size */), processColumnOp, tbb::simple_partitioner());
		//processColumnOp(tbb::blocked_range<unsigned int>(0, cols)); // single-thread

		LOG_DEBUG << "OCL DELAY-STORE " << delayStoreTimer.getTime();

		try {
			Timer transfTimer;

			//==================================================
			// [OpenCL] Memory transfer to device.
			//==================================================
			clCommandQueue_.enqueueWriteBuffer(
				rawDataCLBuffer_, CL_FALSE /* blocking */, 0 /* offset */,
				rawDataSizeInBytes_, mappedRawDataPtrList_[rawBufferIdx],
				nullptr, &writeBufferEventList[rawBufferIdx]);

			LOG_DEBUG << "OCL TRANSF " << transfTimer.getTime(); // useful only if the command was run with blocking activated
		} catch (cl::Error& e) {
			LOG_ERROR << "[oclCommandQueue_.enqueueWriteBuffer()] OpenCL error: " << e.what() << " (" << e.err() << ").";
			throw;
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
				cl::NullRange, /* offset range / must be null */
				cl::NDRange(rawDataN2_, rawDataN1_), /* global range, defines the total number of work-items */
				cl::NDRange(OCL_TRANSPOSE_GROUP_SIZE_DIM_0, OCL_TRANSPOSE_GROUP_SIZE_DIM_0), /* local range, defines the number of work-items in a work-group */
				nullptr /* events */, &kernelEvent);
			//kernelEvent.wait();
			LOG_DEBUG << "OCL TRANSPOSE " << transpTimer.getTime(); // useful only with kernelEvent.wait()
#endif
			Timer procTimer;
			//==================================================
			// [OpenCL] Final processing kernel.
			//==================================================
			clCommandQueue_.enqueueNDRangeKernel(
				procImageKernel,
				cl::NullRange, /* offset range / must be null */
				cl::NDRange(procImageKernelGlobalSize), /* global range, defines the total number of work-items */
				cl::NDRange(OCL_WORK_ITEMS_PER_GROUP), /* local range, defines the number of work-items in a work-group */
				nullptr /* events */, &kernelEvent);
			//kernelEvent.wait();
			LOG_DEBUG << "OCL PROC " << procTimer.getTime(); // useful only with kernelEvent.wait()

		} catch (cl::Error& e) {
			LOG_ERROR << "OpenCL error: " << e.what() << " (" << e.err() << ").";
			throw;
		}
	}

	//==================================================
	// [OpenCL] Read the formed image.
	//==================================================
	clCommandQueue_.enqueueReadBuffer(
		gridValueReCLBuffer_, CL_FALSE /* blocking */, 0 /* offset */,
		Util::sizeInBytes(gridValueRe_), gridValueRe_.data());
	clCommandQueue_.enqueueReadBuffer(
		gridValueImCLBuffer_, CL_TRUE /* blocking */, 0 /* offset */,
		Util::sizeInBytes(gridValueIm_), gridValueIm_.data());
	for (unsigned int col = 0; col < cols; ++col) {
		for (unsigned int row = 0; row < minRowIdx_[col]; ++row) {
			gridValue(col, row) = 0;
		}
		unsigned int gridPointIdx = firstGridPointIdx_[col];
		for (unsigned int row = minRowIdx_[col]; row < gridXZ.n2(); ++row, ++gridPointIdx) {
			gridValue(col, row) = std::complex<TFloat>(gridValueRe_[gridPointIdx], gridValueIm_[gridPointIdx]);
		}
	}

#ifdef USE_EXECUTION_TIME_MEASUREMENT
	tProcessColumn.put(processColumnTimer.getTime());
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
	const TFloat invC1;
	const TFloat invC2;
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
						} else {
#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_OCL_PROCESSOR_USE_TRANSPOSE
							rawData[gridPointIdx * rawDataN2 + rxIdx    ] = 0;
							rawData[gridPointIdx * rawDataN2 + rxIdx + 1] = 0;
#else
							rawData[ rxIdx      * rawDataN2 + gridPointIdx] = 0;
							rawData[(rxIdx + 1) * rawDataN2 + gridPointIdx] = 0;
#endif
						}
					} else {
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
VectorialCombinedTwoMediumImagingOCLProcessor<TFloat>::getKernel() const
{
	return R"CLC(

#define GROUP_SIZE 16

// NVIDIA sm_12:
//   - Local (shared) memory has 16 banks.
__kernel
void
transposeKernel(
		__global float* rawData,
		__global float* rawDataT,
		unsigned int oldSizeX,
		unsigned int oldSizeY,
		__local float* temp) // GROUP_SIZE * GROUP_SIZE
{
	unsigned int iX = get_global_id(0);
	unsigned int iY = get_global_id(1);
	if (iX < oldSizeX && iY < oldSizeY) {
		temp[get_local_id(0) + GROUP_SIZE * get_local_id(1)] = rawData[iX + oldSizeX * iY];
	}

	barrier(CLK_LOCAL_MEM_FENCE);

	iX = get_group_id(1) * GROUP_SIZE + get_local_id(0);
	iY = get_group_id(0) * GROUP_SIZE + get_local_id(1);
	if (iX < oldSizeY && iY < oldSizeX) {
		rawDataT[iX + oldSizeY * iY] = temp[GROUP_SIZE * get_local_id(0) + get_local_id(1)];
	}
}

__kernel
void
processImageKernel(
		__global float* rawData,
		__global float* gridValueRe,
		__global float* gridValueIm,
		__constant float* rxApod,
		unsigned int numGridPoints)
{
	float rxSignalListRe[NUM_RX_ELEM];
	float rxSignalListIm[NUM_RX_ELEM];

	const unsigned int point = get_global_id(0);
	if (point >= numGridPoints) return;

	for (unsigned int i = 0; i < NUM_RX_ELEM; ++i) {
		rxSignalListRe[i] = rawData[ (i << 1)      * numGridPoints + point];
		rxSignalListIm[i] = rawData[((i << 1) + 1) * numGridPoints + point];
		//rxSignalListRe[i] = rawData[point * (NUM_RX_ELEM << 1) + (i << 1)];
		//rxSignalListIm[i] = rawData[point * (NUM_RX_ELEM << 1) + ((i << 1) + 1)];
	}

	float sumRe = 0.0f;
	float sumIm = 0.0f;
	for (unsigned int i = 0; i < NUM_RX_ELEM; ++i) {
		sumRe += rxSignalListRe[i] * rxApod[i];
		sumIm += rxSignalListIm[i] * rxApod[i];
	}

	gridValueRe[point] += sumRe;
	gridValueIm[point] += sumIm;
}

float
arithmeticMean(float* data)
{
	float sum = 0.0f;
	for (unsigned int i = 0; i < NUM_RX_ELEM; ++i) {
		sum += data[i];
	}
	return sum * (1.0f / NUM_RX_ELEM);
}

float
standardDeviation(float* data)
{
	float sum = 0.0f;
	float mean = arithmeticMean(data);
	for (unsigned int i = 0; i < NUM_RX_ELEM; ++i) {
		float e = data[i] - mean;
		sum += e * e;
	}
	return sqrt(sum * (1.0f / NUM_RX_ELEM));
}

float
calcPCF(float* re, float* im, float factor)
{
	float phi[NUM_RX_ELEM];
	float phiAux[NUM_RX_ELEM];

#pragma unroll
	for (unsigned int i = 0; i < NUM_RX_ELEM; ++i) {
		phi[i] = atan2(im[i], re[i]);
	}
	for (unsigned int i = 0; i < NUM_RX_ELEM; ++i) {
		phiAux[i] = phi[i] - copysign(M_PI_F, phi[i]);
	}

	float sf = fmin(standardDeviation(phi), standardDeviation(phiAux));
	return fmax(0.0f, 1.0f - factor * sf);
}

__kernel
void
processImagePCFKernel(
		__global float* rawData,
		__global float* gridValueRe,
		__global float* gridValueIm,
		__constant float* rxApod,
		unsigned int numGridPoints,
		float pcfFactor)
{
	float rxSignalListRe[NUM_RX_ELEM];
	float rxSignalListIm[NUM_RX_ELEM];

	const unsigned int point = get_global_id(0);
	if (point >= numGridPoints) return;

	for (unsigned int i = 0; i < NUM_RX_ELEM; ++i) {
		rxSignalListRe[i] = rawData[ (i << 1)      * numGridPoints + point];
		rxSignalListIm[i] = rawData[((i << 1) + 1) * numGridPoints + point];
	}

	float pcf = calcPCF(rxSignalListRe, rxSignalListIm, pcfFactor);

	float sumRe = 0.0f;
	float sumIm = 0.0f;
	for (unsigned int i = 0; i < NUM_RX_ELEM; ++i) {
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
