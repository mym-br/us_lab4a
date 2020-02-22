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

#ifndef VECTORIALCOMBINEDTWOMEDIUMIMAGINGPROCESSOR4_H
#define VECTORIALCOMBINEDTWOMEDIUMIMAGINGPROCESSOR4_H

#include <algorithm> /* copy, fill, min */
#include <cmath> /* abs, ceil, sqrt */
#include <complex>
#include <cstddef> /* std::size_t */
#include <cstring> /* memset */
#include <memory>
#include <string>
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
#include "HilbertEnvelope.h"
#include "Interpolator.h"
#include "Log.h"
#include "Matrix.h"
#include "PseudorandomNumberGenerator.h"
#include "Tensor3.h"
#include "Timer.h"
#include "TwoMediumSTAConfiguration.h"
#include "Util.h"
#include "XZ.h"



// Depends on the signal.
// 1.0 --> pi radian / sample at the original sampling rate.
#define VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_UPSAMP_FILTER_HALF_TRANSITION_WIDTH (0.2)

// Faster.
#define VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_USE_PINNED_MEMORY_FOR_GPU 1
#define VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_USE_DOUBLE_BUFFER_FOR_GPU 1

// Faster in GPU.
#define VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_USE_TRANSPOSE 1

//#define VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_USE_PRNG 1
//#define VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_USE_CPU 1
#define VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_CPU_COEF (0.0)
#define VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_GPU_PLATFORM 0
#define VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_CPU_PLATFORM 1

#ifdef USE_SIMD
# include "SIMD.h"
#endif

#ifdef USE_OPENCL
# include <CL/cl2.hpp>
# define VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_PROGRAM_BUILD_OPTIONS "-cl-std=CL1.2 -DNUM_RX_ELEM=32"
# ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_USE_DOUBLE_BUFFER_FOR_GPU
#  ifndef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_USE_PINNED_MEMORY_FOR_GPU
#   error Invalid configuration.
#  endif
# endif
#endif



namespace Lab {

// The grid must be rectangular.
// Requirements:
// - Single precision
// - numElements = 32
template<typename TFloat>
class VectorialCombinedTwoMediumImagingProcessor4 {
public:
	struct StepConfiguration {
		unsigned int baseElemIdx;
		unsigned int baseElem;
		unsigned int txElem;
	};

	VectorialCombinedTwoMediumImagingProcessor4(
			const TwoMediumSTAConfiguration<TFloat>& config,
			std::vector<Matrix<TFloat>>& acqDataList,
			unsigned int upsamplingFactor,
			AnalyticSignalCoherenceFactorProcessor<TFloat>& coherenceFactor,
			TFloat maxFermatBlockSize,
			TFloat peakOffset,
			unsigned int ascanStartOffset);
	~VectorialCombinedTwoMediumImagingProcessor4();

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

	struct ProcessColumn2WithOneTxElem;
	struct ProcessColumn2ThreadData {
		AnalyticSignalCoherenceFactorProcessor<TFloat> coherenceFactor;
		std::vector<std::complex<TFloat>, tbb::cache_aligned_allocator<std::complex<TFloat>>> rxSignalSumList;
	};

	VectorialCombinedTwoMediumImagingProcessor4(const VectorialCombinedTwoMediumImagingProcessor4&) = delete;
	VectorialCombinedTwoMediumImagingProcessor4& operator=(const VectorialCombinedTwoMediumImagingProcessor4&) = delete;
	VectorialCombinedTwoMediumImagingProcessor4(VectorialCombinedTwoMediumImagingProcessor4&&) = delete;
	VectorialCombinedTwoMediumImagingProcessor4& operator=(VectorialCombinedTwoMediumImagingProcessor4&&) = delete;

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
	Tensor3<std::complex<TFloat>, tbb::cache_aligned_allocator<std::complex<TFloat>>> signalMatrix_;
	TFloat signalOffset_;
	std::vector<unsigned int, tbb::cache_aligned_allocator<unsigned int>> minRowIdx_; // for each column
	std::vector<unsigned int, tbb::cache_aligned_allocator<unsigned int>> firstGridPointIdx_; // for each column
	std::vector<TFloat, tbb::cache_aligned_allocator<TFloat>> xArray_;
	Matrix<TFloat, tbb::cache_aligned_allocator<TFloat>> medium1DelayMatrix_; // (interface_idx, element)
	Tensor3<TFloat, tbb::cache_aligned_allocator<TFloat>> delayMatrix_;
	Matrix<TFloat, tbb::cache_aligned_allocator<TFloat>> rawDataMatrix_;
	unsigned int rawDataN1_;
	unsigned int rawDataN2_;
	std::unique_ptr<tbb::enumerable_thread_specific<PrepareDataThreadData>> prepareDataTLS_;
	std::unique_ptr<tbb::enumerable_thread_specific<ProcessColumn2ThreadData>> processColumn2TLS_;
#ifdef USE_OPENCL
	bool oclDataInitialized_;
	std::vector<TFloat, tbb::cache_aligned_allocator<TFloat>> gridValueRe_;
	std::vector<TFloat, tbb::cache_aligned_allocator<TFloat>> gridValueIm_;
	cl::Context oclContext_;
	cl::Program oclProgram_;
	cl::CommandQueue oclCommandQueue_;
# ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_USE_PINNED_MEMORY_FOR_GPU
	cl::Buffer oclPinnedRawData_;
	TFloat* mappedRawData_;
#  ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_USE_DOUBLE_BUFFER_FOR_GPU
	cl::Buffer oclPinnedRawData2_;
	TFloat* mappedRawData2_;
#  endif
# endif
	cl::Buffer oclRawData_;
# ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_USE_TRANSPOSE
	cl::Buffer oclRawDataT_;
# endif
	cl::Buffer oclGridValueRe_;
	cl::Buffer oclGridValueIm_;
	cl::Buffer oclRxApod_;
# ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_USE_PRNG
	std::vector<int, tbb::cache_aligned_allocator<int>> initialPrngState_;
	cl::Buffer oclPrngState_;
# endif
# ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_USE_CPU
	unsigned int cpuRawDataN1_;
	unsigned int cpuRawDataN2_;
	std::vector<TFloat, tbb::cache_aligned_allocator<TFloat>> cpuGridValueRe_;
	std::vector<TFloat, tbb::cache_aligned_allocator<TFloat>> cpuGridValueIm_;
	cl::Context oclCpuContext_;
	cl::Program oclCpuProgram_;
	cl::CommandQueue oclCpuCommandQueue_;
	TFloat* mappedCpuRawData_;
	cl::Buffer oclCpuRawData_;
#  ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_USE_TRANSPOSE
	cl::Buffer oclCpuRawDataT_;
#  endif
	cl::Buffer oclCpuGridValueRe_;
	cl::Buffer oclCpuGridValueIm_;
	cl::Buffer oclCpuRxApod_;
#  ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_USE_PRNG
	std::vector<int, tbb::cache_aligned_allocator<int>> initialCpuPrngState_;
	cl::Buffer oclCpuPrngState_;
#  endif
# endif
#endif
};



template<typename TFloat>
VectorialCombinedTwoMediumImagingProcessor4<TFloat>::VectorialCombinedTwoMediumImagingProcessor4(
			const TwoMediumSTAConfiguration<TFloat>& config,
			std::vector<Matrix<TFloat>>& acqDataList,
			unsigned int upsamplingFactor,
			AnalyticSignalCoherenceFactorProcessor<TFloat>& coherenceFactor,
			TFloat maxFermatBlockSize,
			TFloat peakOffset,
			unsigned int ascanStartOffset)
		: config_(config)
		, acqDataList_(acqDataList)
		, upsamplingFactor_(upsamplingFactor)
		, coherenceFactor_(coherenceFactor)
		, maxFermatBlockSize_(maxFermatBlockSize)
		, lambda2_(config_.propagationSpeed2 / config_.centerFrequency)
		, rawDataN1_()
		, rawDataN2_()
#ifdef USE_OPENCL
		, oclDataInitialized_()
# ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_USE_PINNED_MEMORY_FOR_GPU
		, mappedRawData_()
#  ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_USE_DOUBLE_BUFFER_FOR_GPU
		, mappedRawData2_()
#  endif
# endif
# ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_USE_CPU
		, cpuRawDataN1_()
		, cpuRawDataN2_()
		, mappedCpuRawData_()
# endif
#endif
{
	if (sizeof(TFloat) != sizeof(float)) {
		THROW_EXCEPTION(InvalidParameterException, "Only single precision is supported.");
	}

	const std::size_t ascanLength = acqDataList_[0].n2();

	signalOffset_ = (config_.samplingFrequency * upsamplingFactor_) * (peakOffset / config_.centerFrequency) - ascanStartOffset * upsamplingFactor_;
	signalLength_ = ascanLength * upsamplingFactor_;
	LOG_DEBUG << "signalOffset_: " << signalOffset_ << " signalLength_: " << signalLength_;

	PrepareDataThreadData prepareDataThreadData;
	if (upsamplingFactor_ > 1) {
		prepareDataThreadData.interpolator.prepare(upsamplingFactor_, VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_UPSAMP_FILTER_HALF_TRANSITION_WIDTH);
	}
	prepareDataThreadData.signal.resize(signalLength_);
	prepareDataTLS_ = std::make_unique<tbb::enumerable_thread_specific<PrepareDataThreadData>>(prepareDataThreadData);

	ProcessColumn2ThreadData processColumn2ThreadData;
	processColumn2ThreadData.coherenceFactor = coherenceFactor_;
	processColumn2TLS_ = std::make_unique<tbb::enumerable_thread_specific<ProcessColumn2ThreadData>>(processColumn2ThreadData);

#ifdef USE_OPENCL
	std::vector<cl::Platform> platforms;
	cl::Platform::get(&platforms);
	if (platforms.empty()) {
		THROW_EXCEPTION(UnavailableResourceException, "No OpenCL platforms available.");
	}
# ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_USE_CPU
	if (platforms.size() < 2) {
		THROW_EXCEPTION(UnavailableResourceException, "Not enough platforms (n < 2).");
	}
# endif

	if (Log::isDebugEnabled()) {
		try {
			for (std::size_t i = 0; i < platforms.size(); ++i) {
				std::string name = platforms[i].getInfo<CL_PLATFORM_NAME>();
				LOG_DEBUG << "OpenCL platform: " << name;

				std::vector<cl::Device> devices;
				platforms[i].getDevices(CL_DEVICE_TYPE_ALL, &devices);
				if (devices.empty()) {
					THROW_EXCEPTION(UnavailableResourceException, "No OpenCL devices available for platform " << name << '.');
				}

				for (std::size_t j = 0; j < devices.size(); ++j) {
					cl_device_type deviceType = devices[i].getInfo<CL_DEVICE_TYPE>();
					std::string devName = devices[i].getInfo<CL_DEVICE_NAME>();
					LOG_DEBUG << "  device name: " << devName;
					switch (deviceType) {
					case CL_DEVICE_TYPE_CPU:
						LOG_DEBUG << "         type: CPU";
						break;
					case CL_DEVICE_TYPE_GPU:
						LOG_DEBUG << "         type: GPU";
						break;
					default:
						LOG_DEBUG << "         type: other (" << deviceType << ").";
					}
				}
			}
		} catch (cl::Error& e) {
			LOG_ERROR << "Error occurred while obtaining information about OpenCL devices: " << e.what() << " (" << e.err() << ").";
		}
	}

	cl_context_properties properties[] = {
					CL_CONTEXT_PLATFORM,
					(cl_context_properties)(platforms[VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_GPU_PLATFORM])(),
					0};
	oclContext_ = cl::Context(CL_DEVICE_TYPE_GPU, properties);
	std::vector<std::string> kernelStrings;
	kernelStrings.push_back(getKernel());
	oclProgram_ = cl::Program(oclContext_, kernelStrings);
	std::vector<cl::Device> devices = oclContext_.getInfo<CL_CONTEXT_DEVICES>();
	try {
		oclProgram_.build(devices, VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_PROGRAM_BUILD_OPTIONS);
	} catch (cl::Error& e) {
		try {
			LOG_ERROR << "Error in cl::Program.build(): " << e.what();
			std::string msg = oclProgram_.getBuildInfo<CL_PROGRAM_BUILD_OPTIONS>(devices[0]);
			LOG_ERROR << "Build options:\n" << msg;
			msg = oclProgram_.getBuildInfo<CL_PROGRAM_BUILD_LOG>(devices[0]);
			LOG_ERROR << "Build log:\n" << msg;
		} catch (...) {} // ignore
		THROW_EXCEPTION(Exception, "[oclContext_] Error in cl::Program.build(): " << e.what());
	}

	oclCommandQueue_ = cl::CommandQueue(oclContext_, devices[0]);

# ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_USE_CPU
	cl_context_properties cpuProperties[] = {
					CL_CONTEXT_PLATFORM,
					(cl_context_properties)(platforms[VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_CPU_PLATFORM])(),
					0};
	oclCpuContext_ = cl::Context(CL_DEVICE_TYPE_CPU, cpuProperties);
	oclCpuProgram_ = cl::Program(oclCpuContext_, sources);
	std::vector<cl::Device> cpuDevices = oclCpuContext_.getInfo<CL_CONTEXT_DEVICES>();
	try {
		oclCpuProgram_.build(cpuDevices, VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_PROGRAM_BUILD_OPTIONS);
	} catch (cl::Error& e) {
		try {
			LOG_ERROR << "Error in cl::Program.build() (CPU): " << e.what();
			std::string msg = oclCpuProgram_.getBuildInfo<CL_PROGRAM_BUILD_OPTIONS>(cpuDevices[0]);
			LOG_ERROR << "Build options:\n" << msg;
			msg = oclCpuProgram_.getBuildInfo<CL_PROGRAM_BUILD_LOG>(cpuDevices[0]);
			LOG_ERROR << "Build log:\n" << msg;
		} catch (...) {} // ignore
		THROW_EXCEPTION(Exception, "[oclCpuContext_] Error in cl::Program.build(): " << e.what());
	}

	oclCpuCommandQueue_ = cl::CommandQueue(oclCpuContext_, cpuDevices[0]);
# endif
#endif
}

template<typename TFloat>
VectorialCombinedTwoMediumImagingProcessor4<TFloat>::~VectorialCombinedTwoMediumImagingProcessor4()
{
#ifdef USE_OPENCL
# ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_USE_PINNED_MEMORY_FOR_GPU
	if (oclPinnedRawData_() && mappedRawData_ != nullptr) {
		LOG_DEBUG << "~VectorialCombinedTwoMediumImagingProcessor4: enqueueUnmapMemObject (GPU)";
		try {
			if (oclCommandQueue_()) {
				cl::Event event;
				oclCommandQueue_.enqueueUnmapMemObject(oclPinnedRawData_, mappedRawData_, nullptr, &event);
				event.wait();
			}
		} catch (cl::Error& e) {
			LOG_ERROR << "[mappedRawData_] Error in cl::CommandQueue.enqueueUnmapMemObject(): " << e.what();
		} catch (...) {
			LOG_ERROR << "[mappedRawData_] Caught an unknown exception.";
		}
	}
#  ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_USE_DOUBLE_BUFFER_FOR_GPU
	if (oclPinnedRawData2_() && mappedRawData2_ != nullptr) {
		LOG_DEBUG << "~VectorialCombinedTwoMediumImagingProcessor4: enqueueUnmapMemObject (GPU buffer 2)";
		try {
			if (oclCommandQueue_()) {
				cl::Event event;
				oclCommandQueue_.enqueueUnmapMemObject(oclPinnedRawData2_, mappedRawData2_, nullptr, &event);
				event.wait();
			}
		} catch (cl::Error& e) {
			LOG_ERROR << "[mappedRawData2_] Error in cl::CommandQueue.enqueueUnmapMemObject(): " << e.what();
		} catch (...) {
			LOG_ERROR << "[mappedRawData2_] Caught an unknown exception.";
		}
	}
#  endif
# endif
# ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_USE_CPU
	if (oclCpuRawData_() && mappedCpuRawData_ != nullptr) {
		LOG_DEBUG << "~VectorialCombinedTwoMediumImagingProcessor4: enqueueUnmapMemObject (CPU)";
		try {
			if (oclCpuCommandQueue_()) {
				cl::Event event;
				oclCpuCommandQueue_.enqueueUnmapMemObject(oclCpuRawData_, mappedCpuRawData_, nullptr, &event);
				event.wait();
			}
		} catch (cl::Error& e) {
			LOG_ERROR << "[mappedCpuRawData_] Error in cl::CommandQueue.enqueueUnmapMemObject(): " << e.what();
		} catch (...) {
			LOG_ERROR << "[mappedCpuRawData_] Caught an unknown exception.";
		}
	}
# endif
#endif
}

template<typename TFloat>
void
VectorialCombinedTwoMediumImagingProcessor4<TFloat>::process(
							const std::vector<StepConfiguration>& stepConfigList,
							const std::vector<XZ<TFloat>>& interfacePointList,
							const std::vector<TFloat>& rxApod,
							const Matrix<XZ<TFloat>>& gridXZ,
							Matrix<std::complex<TFloat>>& gridValue)
{
	//LOG_DEBUG << "BEGIN ========== VectorialCombinedTwoMediumImagingProcessor4::process ==========";

	if (stepConfigList.empty()) {
		THROW_EXCEPTION(InvalidParameterException, "The list of step configurations is empty.");
	}
	if (gridXZ.n1() != gridValue.n1() || gridXZ.n2() != gridValue.n2()) {
		THROW_EXCEPTION(InvalidParameterException, "gridXZ and gridValue have different sizes.");
	}

	const std::size_t samplesPerChannelLow = acqDataList_[0].n2();

	minRowIdx_.resize(gridXZ.n1() /* number of columns */);
	firstGridPointIdx_.resize(gridXZ.n1() /* number of columns */);
	delayMatrix_.resize(gridXZ.n1() /* number of columns */, gridXZ.n2() /* number of rows */, config_.numElementsMux);
	signalMatrix_.resize(stepConfigList.size(), config_.numElements, signalLength_);
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

#if defined USE_OPENCL && defined VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_USE_CPU
	const unsigned int cpuCols = std::min(
					static_cast<unsigned int>(gridXZ.n1()),
					static_cast<unsigned int>(VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_CPU_COEF * gridXZ.n1() + 0.5));
	const unsigned int cols = gridXZ.n1() - cpuCols;
	std::size_t pointSum = 0;
	for (unsigned int col = 0; col < cols; ++col) {
		pointSum += gridXZ.n2() - minRowIdx_[col];
	}
	const std::size_t numGridPoints = pointSum;
	std::size_t cpuPointSum = 0;
	for (unsigned int col = 0; col < cpuCols; ++col) {
		cpuPointSum += gridXZ.n2() - minRowIdx_[cols + col];
	}
	const std::size_t numCpuGridPoints = cpuPointSum;
	LOG_DEBUG << "cols: " << cols << " cpuCols: " << cpuCols << " gridXZ.n1(): " << gridXZ.n1() <<
			" numGridPoints: " << numGridPoints << " numCpuGridPoints: " << numCpuGridPoints;
#else
	const unsigned int cols = gridXZ.n1();
	std::size_t pointSum = 0;
	for (unsigned int col = 0; col < cols; ++col) {
		pointSum += gridXZ.n2() - minRowIdx_[col];
	}
	const std::size_t numGridPoints = pointSum;
	LOG_DEBUG << "cols: " << cols << " numGridPoints: " << numGridPoints;
#endif

#ifdef USE_OPENCL
# ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_USE_TRANSPOSE
	if (cols > 0) {
		const std::size_t transpNumGridPoints = roundUpToMultipleOfGroupSize(numGridPoints, OCL_TRANSPOSE_GROUP_SIZE_DIM_0);
		LOG_DEBUG << "numGridPoints: " << numGridPoints << " transpNumGridPoints: " << transpNumGridPoints;
		rawDataN1_ = transpNumGridPoints;
		rawDataN2_ = 2 * config_.numElements /* real, imag */;
#  ifndef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_USE_PINNED_MEMORY_FOR_GPU
		rawDataMatrix_.resize(rawDataN1_, rawDataN2_);
#  endif
	}
#  ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_USE_CPU
	if (cpuCols > 0) {
		const std::size_t transpNumCpuGridPoints = roundUpToMultipleOfGroupSize(numCpuGridPoints, OCL_TRANSPOSE_GROUP_SIZE_DIM_0);
		LOG_DEBUG << "numCpuGridPoints: " << numCpuGridPoints << " transpNumCpuGridPoints: " << transpNumCpuGridPoints;
		cpuRawDataN1_ = transpNumCpuGridPoints;
		cpuRawDataN2_ = 2 * config_.numElements /* real, imag */;
	}
#  endif
# else
	if (cols > 0) {
		rawDataN1_ = 2 * config_.numElements /* real, imag */;
		rawDataN2_ = numGridPoints;
#  ifndef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_USE_PINNED_MEMORY_FOR_GPU
		rawDataMatrix_.resize(rawDataN1_, rawDataN2_);
#  endif
	}
#  ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_USE_CPU
	if (cpuCols > 0) {
		cpuRawDataN1_ = 2 * config_.numElements /* real, imag */;
		cpuRawDataN2_ = numCpuGridPoints;
	}
#  endif
# endif
#else
	rawDataN1_ = 2 * config_.numElements /* real, imag */;
	rawDataN2_ = numGridPoints;
	rawDataMatrix_.resize(rawDataN1_, rawDataN2_);
#endif

#ifdef USE_OPENCL
	if (cols > 0) {
		gridValueRe_.resize(numGridPoints);
		gridValueIm_.resize(numGridPoints);
	}
# ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_USE_CPU
	if (cpuCols > 0) {
		cpuGridValueRe_.resize(numCpuGridPoints);
		cpuGridValueIm_.resize(numCpuGridPoints);
	}
# endif
	if (!oclDataInitialized_) {
		if (cols > 0) {
# ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_USE_PINNED_MEMORY_FOR_GPU
			oclPinnedRawData_ = cl::Buffer(oclContext_, CL_MEM_ALLOC_HOST_PTR, rawDataN1_ * rawDataN2_ * sizeof(TFloat));
			mappedRawData_ = static_cast<TFloat*>(oclCommandQueue_.enqueueMapBuffer(
									oclPinnedRawData_, CL_TRUE /* blocking */, CL_MAP_WRITE,
									0 /* offset */, rawDataN1_ * rawDataN2_ * sizeof(TFloat)));
#  ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_USE_DOUBLE_BUFFER_FOR_GPU
			oclPinnedRawData2_ = cl::Buffer(oclContext_, CL_MEM_ALLOC_HOST_PTR, rawDataN1_ * rawDataN2_ * sizeof(TFloat));
			mappedRawData2_ = static_cast<TFloat*>(oclCommandQueue_.enqueueMapBuffer(
									oclPinnedRawData2_, CL_TRUE /* blocking */, CL_MAP_WRITE,
									0 /* offset */, rawDataN1_ * rawDataN2_ * sizeof(TFloat)));
#  endif
# endif
			oclRawData_     = cl::Buffer(oclContext_, CL_MEM_READ_ONLY , rawDataN1_ * rawDataN2_ * sizeof(TFloat));
# ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_USE_TRANSPOSE
			oclRawDataT_    = cl::Buffer(oclContext_, CL_MEM_READ_WRITE, rawDataN1_ * rawDataN2_ * sizeof(TFloat));
# endif
			oclGridValueRe_ = cl::Buffer(oclContext_, CL_MEM_READ_WRITE, numGridPoints           * sizeof(TFloat));
			oclGridValueIm_ = cl::Buffer(oclContext_, CL_MEM_READ_WRITE, numGridPoints           * sizeof(TFloat));
			oclRxApod_ =      cl::Buffer(oclContext_, CL_MEM_READ_ONLY , rxApod.size()           * sizeof(TFloat));
# ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_USE_PRNG
			if (coherenceFactor_.enabled()) {
				oclPrngState_ = cl::Buffer(oclContext_, CL_MEM_READ_WRITE, numGridPoints * sizeof(int));
			}
# endif
		}
# ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_USE_CPU
		if (cpuCols > 0) {
			oclCpuRawData_     = cl::Buffer(oclCpuContext_, CL_MEM_READ_ONLY | CL_MEM_ALLOC_HOST_PTR, cpuRawDataN1_ * cpuRawDataN2_ * sizeof(TFloat));
#  ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_USE_TRANSPOSE
			oclCpuRawDataT_    = cl::Buffer(oclCpuContext_, CL_MEM_READ_WRITE                       , cpuRawDataN1_ * cpuRawDataN2_ * sizeof(TFloat));
#  endif
			oclCpuGridValueRe_ = cl::Buffer(oclCpuContext_, CL_MEM_READ_WRITE, numCpuGridPoints * sizeof(TFloat));
			oclCpuGridValueIm_ = cl::Buffer(oclCpuContext_, CL_MEM_READ_WRITE, numCpuGridPoints * sizeof(TFloat));
			oclCpuRxApod_      = cl::Buffer(oclCpuContext_, CL_MEM_READ_ONLY , rxApod.size()    * sizeof(TFloat));
#  ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_USE_PRNG
			if (coherenceFactor_.enabled()) {
				oclCpuPrngState_ = cl::Buffer(oclCpuContext_, CL_MEM_READ_WRITE, numCpuGridPoints * sizeof(int));
			}
#  endif
		}
# endif
		oclDataInitialized_ = true;
	}

	if (cols > 0) {
# ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_USE_PINNED_MEMORY_FOR_GPU
		std::memset(mappedRawData_, 0, rawDataN1_ * rawDataN2_ * sizeof(TFloat));
#  ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_USE_DOUBLE_BUFFER_FOR_GPU
		std::memset(mappedRawData2_, 0, rawDataN1_ * rawDataN2_ * sizeof(TFloat));
#  endif
# else
		std::memset(&rawDataMatrix_(0, 0), 0, rawDataN1_ * rawDataN2_ * sizeof(TFloat));
# endif
		std::memset(&gridValueRe_[0], 0, gridValueRe_.size() * sizeof(TFloat));
		std::memset(&gridValueIm_[0], 0, gridValueIm_.size() * sizeof(TFloat));
		oclCommandQueue_.enqueueWriteBuffer(
			oclGridValueRe_, CL_TRUE /* blocking */, 0 /* offset */,
			gridValueRe_.size() * sizeof(TFloat), &gridValueRe_[0]);
		oclCommandQueue_.enqueueWriteBuffer(
			oclGridValueIm_, CL_TRUE /* blocking */, 0 /* offset */,
			gridValueIm_.size() * sizeof(TFloat), &gridValueIm_[0]);
		oclCommandQueue_.enqueueWriteBuffer(
			oclRxApod_, CL_TRUE /* blocking */, 0 /* offset */,
			rxApod.size() * sizeof(TFloat), &rxApod[0]);
# ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_USE_PRNG
		if (coherenceFactor_.enabled()) {
			MinstdPseudorandomNumberGenerator prng(1);
			initialPrngState_.resize(numGridPoints);
			for (auto& e : initialPrngState_) {
				union {
					float f;
					int i;
				} u;
				u.f = prng.get();
				int v = std::abs(u.i);
				e = (v % 2147483646) + 1;
			}
			oclCommandQueue_.enqueueWriteBuffer(
				oclPrngState_, CL_TRUE /* blocking */, 0 /* offset */,
				initialPrngState_.size() * sizeof(int), &initialPrngState_[0]);
		}
# endif
	}
# ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_USE_CPU
	if (cpuCols > 0) {

		mappedCpuRawData_ = static_cast<TFloat*>(oclCpuCommandQueue_.enqueueMapBuffer(
								oclCpuRawData_, CL_TRUE /* blocking */, CL_MAP_WRITE,
								0 /* offset */, cpuRawDataN1_ * cpuRawDataN2_ * sizeof(TFloat)));
		std::memset(mappedCpuRawData_, 0, cpuRawDataN1_ * cpuRawDataN2_ * sizeof(TFloat));
		cl::Event unmapEvent;
		oclCpuCommandQueue_.enqueueUnmapMemObject(oclCpuRawData_, mappedCpuRawData_, nullptr, &unmapEvent);
		unmapEvent.wait();
		mappedCpuRawData_ = nullptr;

		std::memset(&cpuGridValueRe_[0], 0, cpuGridValueRe_.size() * sizeof(TFloat));
		std::memset(&cpuGridValueIm_[0], 0, cpuGridValueIm_.size() * sizeof(TFloat));
		oclCpuCommandQueue_.enqueueWriteBuffer(
			oclCpuGridValueRe_, CL_TRUE /* blocking */, 0 /* offset */,
			cpuGridValueRe_.size() * sizeof(TFloat), &cpuGridValueRe_[0]);
		oclCpuCommandQueue_.enqueueWriteBuffer(
			oclCpuGridValueIm_, CL_TRUE /* blocking */, 0 /* offset */,
			cpuGridValueIm_.size() * sizeof(TFloat), &cpuGridValueIm_[0]);
		oclCpuCommandQueue_.enqueueWriteBuffer(
			oclCpuRxApod_, CL_TRUE /* blocking */, 0 /* offset */,
			rxApod.size() * sizeof(TFloat), &rxApod[0]);
#  ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_USE_PRNG
		if (coherenceFactor_.enabled()) {
			MinstdPseudorandomNumberGenerator prng(1);
			initialCpuPrngState_.resize(numCpuGridPoints);
			for (auto& e : initialPrngState_) {
				union {
					float f;
					int i;
				} u;
				u.f = prng.get();
				int v = std::abs(u.i);
				e = (v % 2147483646) + 1;
			}
			oclCpuCommandQueue_.enqueueWriteBuffer(
				oclCpuPrngState_, CL_TRUE /* blocking */, 0 /* offset */,
				initialCpuPrngState_.size() * sizeof(int), &initialCpuPrngState_[0]);
		}
#  endif
	}
# endif
#else
	rawDataMatrix_ = 0.0;
	gridValue = std::complex<TFloat>(0.0);
#endif // USE_OPENCL

	const TFloat c2ByC1 = config_.propagationSpeed2 / config_.propagationSpeed1;

	ArrayGeometry::getElementX2D(config_.numElementsMux, config_.pitch, TFloat(0) /* offset */, xArray_);

#ifdef USE_EXECUTION_TIME_MEASUREMENT
	Timer medium1DelayMatrixTimer;
#endif
	for (unsigned int elem = 0; elem < config_.numElementsMux; ++elem) {
		TFloat* delays = &medium1DelayMatrix_(elem, 0);
		for (unsigned int i = 0; i < interfacePointList.size(); ++i) {
			const XZ<TFloat>& ifPoint = interfacePointList[i];
#ifdef USE_SIMD
			delays[i] = SIMD::calcDistance(xArray_[elem], 0, ifPoint.x, ifPoint.z) * c2ByC1;
#else
			const TFloat dx = ifPoint.x - xArray_[elem];
			const TFloat dz = ifPoint.z;
			delays[i] = std::sqrt(dx * dx + dz * dz) * c2ByC1;
#endif
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
		delayMatrix_
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
			signalMatrix_
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

#ifdef USE_OPENCL
# ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_USE_TRANSPOSE
	cl::Kernel kernel0;
# endif
	cl::Kernel kernel1;
	if (cols > 0) {
		//==================================================
		// [OpenCL] Kernel preparation - GPU.
		//==================================================
		try {
# ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_USE_TRANSPOSE
			kernel0 = cl::Kernel(oclProgram_, "transposeKernel");
			kernel0.setArg(0, oclRawData_);
			kernel0.setArg(1, oclRawDataT_);
			kernel0.setArg(2, rawDataN2_);
			kernel0.setArg(3, rawDataN1_);
			kernel0.setArg(4, cl::Local(OCL_TRANSPOSE_GROUP_SIZE_DIM_0 * (OCL_TRANSPOSE_GROUP_SIZE_DIM_0 + 1) * sizeof(TFloat)));
# endif
			if (coherenceFactor_.enabled()) {
				std::vector<TFloat> cfConstants;
				coherenceFactor_.implementation().getConstants(cfConstants);

				cl::Kernel kernel1a(oclProgram_, "processImagePCFKernel");
# ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_USE_TRANSPOSE
				kernel1a.setArg(0, oclRawDataT_);
# else
				kernel1a.setArg(0, oclRawData_);
# endif
				kernel1a.setArg(1, oclGridValueRe_);
				kernel1a.setArg(2, oclGridValueIm_);
				kernel1a.setArg(3, oclRxApod_);
# ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_USE_TRANSPOSE
				kernel1a.setArg(4, rawDataN1_);
# else
				kernel1a.setArg(4, rawDataN2_);
# endif
				kernel1a.setArg(5, cfConstants[2] /* factor */);
# ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_USE_PRNG
				kernel1a.setArg(6, oclPrngState_);
# endif
				kernel1 = kernel1a;
			} else {
				cl::Kernel kernel1b(oclProgram_, "processImageKernel");
# ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_USE_TRANSPOSE
				kernel1b.setArg(0, oclRawDataT_);
# else
				kernel1b.setArg(0, oclRawData_);
# endif
				kernel1b.setArg(1, oclGridValueRe_);
				kernel1b.setArg(2, oclGridValueIm_);
				kernel1b.setArg(3, oclRxApod_);
				kernel1b.setArg(4, static_cast<unsigned int>(numGridPoints));
				kernel1 = kernel1b;
			}
		} catch (cl::Error& e) {
			LOG_ERROR << "[Kernel preparation - GPU] OpenCL error: " << e.what() << " (" << e.err() << ").";
			throw;
		}
	}

	cl::Event writeBufferEvent;
# ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_USE_DOUBLE_BUFFER_FOR_GPU
	cl::Event writeBufferEvent2;
# endif
	std::size_t kernel1GlobalSize = roundUpToMultipleOfGroupSize(numGridPoints, OCL_WORK_ITEMS_PER_GROUP);
	//LOG_DEBUG << numGridPoints << ':' << kernel1GlobalSize << ':' << OCL_WORK_ITEMS_PER_GROUP;

# ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_USE_CPU
	std::size_t cpuKernel1GlobalSize = roundUpToMultipleOfGroupSize(numCpuGridPoints, OCL_WORK_ITEMS_PER_GROUP);
	//LOG_DEBUG << numCpuGridPoints << ':' << cpuKernel1GlobalSize << ':' << OCL_WORK_ITEMS_PER_GROUP;

#  ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_USE_TRANSPOSE
	cl::Kernel cpuKernel0;
#  endif
	cl::Kernel cpuKernel1;
	if (cpuCols > 0) {
		//==================================================
		// [OpenCL] Kernel preparation - CPU.
		//==================================================
		try {
#  ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_USE_TRANSPOSE
			cpuKernel0 = cl::Kernel(oclCpuProgram_, "transposeKernel");
			cpuKernel0.setArg(0, oclCpuRawData_);
			cpuKernel0.setArg(1, oclCpuRawDataT_);
			cpuKernel0.setArg(2, cpuRawDataN2_);
			cpuKernel0.setArg(3, cpuRawDataN1_);
			cpuKernel0.setArg(4, cl::__local(OCL_TRANSPOSE_GROUP_SIZE_DIM_0 * (OCL_TRANSPOSE_GROUP_SIZE_DIM_0 + 1) * sizeof(TFloat))); // __local or Local
#  endif
			if (coherenceFactor_.enabled()) {
				std::vector<TFloat> cfConstants;
				coherenceFactor_.implementation().getConstants(cfConstants);

				cl::Kernel cpuKernel1a(oclCpuProgram_, "processImagePCFKernel");
#  ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_USE_TRANSPOSE
				cpuKernel1a.setArg(0, oclCpuRawDataT_);
#  else
				cpuKernel1a.setArg(0, oclCpuRawData_);
#  endif
				cpuKernel1a.setArg(1, oclCpuGridValueRe_);
				cpuKernel1a.setArg(2, oclCpuGridValueIm_);
				cpuKernel1a.setArg(3, oclCpuRxApod_);
#  ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_USE_TRANSPOSE
				cpuKernel1a.setArg(4, cpuRawDataN1_);
#  else
				cpuKernel1a.setArg(4, cpuRawDataN2_);
#  endif
				cpuKernel1a.setArg(5, cfConstants[2] /* factor */);
#  ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_USE_PRNG
				cpuKernel1a.setArg(6, oclCpuPrngState_);
#  endif
				cpuKernel1 = cpuKernel1a;
			} else {
				cl::Kernel cpuKernel1b(oclCpuProgram_, "processImageKernel");
#  ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_USE_TRANSPOSE
				cpuKernel1b.setArg(0, oclCpuRawDataT_);
#  else
				cpuKernel1b.setArg(0, oclCpuRawData_);
#  endif
				cpuKernel1b.setArg(1, oclCpuGridValueRe_);
				cpuKernel1b.setArg(2, oclCpuGridValueIm_);
				cpuKernel1b.setArg(3, oclCpuRxApod_);
				cpuKernel1b.setArg(4, static_cast<unsigned int>(numCpuGridPoints));
				cpuKernel1 = cpuKernel1b;
			}
		} catch (cl::Error& e) {
			LOG_ERROR << "[Kernel preparation - CPU] OpenCL error: " << e.what() << " (" << e.err() << ").";
			throw;
		}
	}
# endif
#endif // USE_OPENCL

	//==================================================
	// Step configuration loop.
	//==================================================
	bool evenIter = true;
#if defined USE_OPENCL && defined VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_USE_DOUBLE_BUFFER_FOR_GPU
	for (unsigned int i = 0; i < stepConfigList.size(); ++i, evenIter = !evenIter) {
#else
	for (unsigned int i = 0; i < stepConfigList.size(); ++i) {
#endif
		const auto& stepConfig = stepConfigList[i];
		LOG_DEBUG << "stepConfig.baseElemIdx: " << stepConfig.baseElemIdx << " evenIter: " << evenIter;

		if (cols > 0) {
#ifdef USE_OPENCL
# ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_USE_DOUBLE_BUFFER_FOR_GPU
			if (evenIter) {
				if (writeBufferEvent() != nullptr) writeBufferEvent.wait();
			} else {
				if (writeBufferEvent2() != nullptr) writeBufferEvent2.wait();
			}
# else
			if (writeBufferEvent() != nullptr) writeBufferEvent.wait();
# endif
#endif
			Timer delayStoreTimer;

			//==================================================
			// Delay and store - CPU/GPU.
			//==================================================
			ProcessColumnWithOneTxElem processColumnOp = {
				static_cast<unsigned int>(gridXZ.n2()),
				0,
				config_,
				signalOffset_,
				signalMatrix_,
				stepConfig,
				minRowIdx_,
				firstGridPointIdx_,
				delayMatrix_,
#if defined USE_OPENCL && defined VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_USE_PINNED_MEMORY_FOR_GPU
# ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_USE_DOUBLE_BUFFER_FOR_GPU
				evenIter ? mappedRawData_ : mappedRawData2_,
# else
				mappedRawData_,
# endif
#else
				&rawDataMatrix_(0, 0),
#endif
				rawDataN2_,
			};

			//tbb::parallel_for(tbb::blocked_range<unsigned int>(0, cols), processColumnOp);
			tbb::parallel_for(tbb::blocked_range<unsigned int>(0, cols, 1 /* grain size */), processColumnOp, tbb::simple_partitioner());
			//processColumnOp(tbb::blocked_range<unsigned int>(0, cols)); // single-thread

			LOG_DEBUG << "OCL DELAY-STORE " << delayStoreTimer.getTime();
		}

#ifdef USE_OPENCL
		if (cols > 0) {
			try {
				Timer transfTimer;

				//==================================================
				// [OpenCL] Memory transfer to GPU.
				//==================================================
# ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_USE_PINNED_MEMORY_FOR_GPU
				oclCommandQueue_.enqueueWriteBuffer(
					oclRawData_, CL_FALSE /* blocking */, 0 /* offset */,
#  ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_USE_DOUBLE_BUFFER_FOR_GPU
					rawDataN1_ * rawDataN2_ * sizeof(TFloat), evenIter ? mappedRawData_ : mappedRawData2_,
					nullptr, evenIter ? &writeBufferEvent : &writeBufferEvent2);
#  else
					rawDataN1_ * rawDataN2_ * sizeof(TFloat), mappedRawData_,
					nullptr, &writeBufferEvent);
#  endif
# else
				oclCommandQueue_.enqueueWriteBuffer(
					oclRawData_, CL_FALSE /* blocking */, 0 /* offset */,
					rawDataN1_ * rawDataN2_ * sizeof(TFloat), &rawDataMatrix_(0, 0),
					nullptr, &writeBufferEvent);
# endif

				LOG_DEBUG << "OCL TRANSF " << transfTimer.getTime(); // useful only if the command was run with blocking activated
			} catch (cl::Error& e) {
				LOG_ERROR << "[oclCommandQueue_.enqueueWriteBuffer()] OpenCL error: " << e.what() << " (" << e.err() << ").";
				throw;
			}
		}

# ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_USE_CPU
		if (cpuCols > 0) {
			mappedCpuRawData_ = static_cast<TFloat*>(oclCpuCommandQueue_.enqueueMapBuffer(
									oclCpuRawData_, CL_TRUE /* blocking */, CL_MAP_WRITE,
									0 /* offset */, cpuRawDataN1_ * cpuRawDataN2_ * sizeof(TFloat)));

			//==================================================
			// Delay and store - OpenCL CPU.
			//==================================================
			ProcessColumnWithOneTxElem cpuProcessColumnOp = {
				static_cast<unsigned int>(gridXZ.n2()),
				cols,
				config_,
				signalOffset_,
				signalMatrix_,
				stepConfig,
				minRowIdx_,
				firstGridPointIdx_,
				delayMatrix_,
				mappedCpuRawData_,
				cpuRawDataN2_,
			};

			//tbb::parallel_for(tbb::blocked_range<unsigned int>(cols, gridXZ.n1()), processColumnOp);
			tbb::parallel_for(tbb::blocked_range<unsigned int>(cols, gridXZ.n1(), 1 /* grain size */), cpuProcessColumnOp, tbb::simple_partitioner());
			//cpuProcessColumnOp(tbb::blocked_range<unsigned int>(cols, gridXZ.n1())); // single-thread

			cl::Event unmapEvent;
			oclCpuCommandQueue_.enqueueUnmapMemObject(oclCpuRawData_, mappedCpuRawData_, nullptr, &unmapEvent);
			unmapEvent.wait();
			mappedCpuRawData_ = nullptr;
		}
# endif
#endif

#ifdef USE_OPENCL
		try {
			cl::Event kernelEvent;

# ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_USE_TRANSPOSE
			Timer transpTimer;

			if (cols > 0) {
				//==================================================
				// [OpenCL] Transpose kernel - GPU.
				//==================================================
				oclCommandQueue_.enqueueNDRangeKernel(
					kernel0,
					cl::NullRange, /* offset range / must be null */
					cl::NDRange(rawDataN2_, rawDataN1_), /* global range, defines the total number of work-items */
					cl::NDRange(OCL_TRANSPOSE_GROUP_SIZE_DIM_0, OCL_TRANSPOSE_GROUP_SIZE_DIM_0), /* local range, defines the number of work-items in a work-group */
					0 /* events */, &kernelEvent);
				//kernelEvent.wait();
			}
#  ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_USE_CPU
			if (cpuCols > 0) {
				//==================================================
				// [OpenCL] Transpose kernel - CPU.
				//==================================================
				oclCpuCommandQueue_.enqueueNDRangeKernel(
					cpuKernel0,
					cl::NullRange, /* offset range / must be null */
					cl::NDRange(cpuRawDataN2_, cpuRawDataN1_), /* global range, defines the total number of work-items */
					cl::NDRange(OCL_TRANSPOSE_GROUP_SIZE_DIM_0, OCL_TRANSPOSE_GROUP_SIZE_DIM_0), /* local range, defines the number of work-items in a work-group */
					0 /* events */, &kernelEvent);
				//kernelEvent.wait();
			}
#  endif
			LOG_DEBUG << "OCL TRANSPOSE " << transpTimer.getTime(); // useful only with kernelEvent.wait()
# endif
			Timer procTimer;

			if (cols > 0) {
				//==================================================
				// [OpenCL] Final processing kernel - GPU.
				//==================================================
				oclCommandQueue_.enqueueNDRangeKernel(
					kernel1,
					cl::NullRange, /* offset range / must be null */
					cl::NDRange(kernel1GlobalSize), /* global range, defines the total number of work-items */
					cl::NDRange(OCL_WORK_ITEMS_PER_GROUP), /* local range, defines the number of work-items in a work-group */
					0 /* events */, &kernelEvent);
				//kernelEvent.wait();
			}

# ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_USE_CPU
			if (cpuCols > 0) {
				//==================================================
				// [OpenCL] Final processing kernel - CPU.
				//==================================================
				oclCpuCommandQueue_.enqueueNDRangeKernel(
					cpuKernel1,
					cl::NullRange, /* offset range / must be null */
					cl::NDRange(cpuKernel1GlobalSize), /* global range, defines the total number of work-items */
					cl::NDRange(OCL_WORK_ITEMS_PER_GROUP), /* local range, defines the number of work-items in a work-group */
					0 /* events */, &kernelEvent);
				//kernelEvent.wait();
			}
# endif
			LOG_DEBUG << "OCL PROC " << procTimer.getTime(); // useful only with kernelEvent.wait()

		} catch (cl::Error& e) {
			LOG_ERROR << "OpenCL error: " << e.what() << " (" << e.err() << ").";
			throw;
		}
#else
		Timer procTimer;

		//==================================================
		// Final processing - CPU.
		//==================================================
		ProcessColumn2WithOneTxElem processColumn2Op = {
			gridXZ.n2(),
			config_,
			rxApod,
			stepConfig,
			minRowIdx_,
			firstGridPointIdx_,
			rawDataMatrix_,
			*processColumn2TLS_,
			gridValue
		};

		tbb::parallel_for(tbb::blocked_range<unsigned int>(0, gridXZ.n1(), 1 /* grain size */), processColumn2Op, tbb::simple_partitioner());
		//processColumn2Op(tbb::blocked_range<unsigned int>(0, gridXZ.n1())); // single-thread

		LOG_DEBUG << "CPU PROC " << procTimer.getTime();
#endif
	}
#ifdef USE_OPENCL
	if (cols > 0) {
		//==================================================
		// [OpenCL] Read the formed image - GPU.
		//==================================================
		oclCommandQueue_.enqueueReadBuffer(
			oclGridValueRe_, CL_FALSE /* blocking */, 0 /* offset */,
			gridValueRe_.size() * sizeof(TFloat), &gridValueRe_[0]);
		oclCommandQueue_.enqueueReadBuffer(
			oclGridValueIm_, CL_TRUE /* blocking */, 0 /* offset */,
			gridValueIm_.size() * sizeof(TFloat), &gridValueIm_[0]);
		for (unsigned int col = 0; col < cols; ++col) {
			for (unsigned int row = 0; row < minRowIdx_[col]; ++row) {
				gridValue(col, row) = 0;
			}
			unsigned int gridPointIdx = firstGridPointIdx_[col];
			for (unsigned int row = minRowIdx_[col]; row < gridXZ.n2(); ++row, ++gridPointIdx) {
				gridValue(col, row) = std::complex<TFloat>(gridValueRe_[gridPointIdx], gridValueIm_[gridPointIdx]);
			}
		}
	}
# ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_USE_CPU
	if (cpuCols > 0) {
		//==================================================
		// [OpenCL] Read the formed image - CPU.
		//==================================================
		oclCpuCommandQueue_.enqueueReadBuffer(
			oclCpuGridValueRe_, CL_FALSE /* blocking */, 0 /* offset */,
			cpuGridValueRe_.size() * sizeof(TFloat), &cpuGridValueRe_[0]);
		oclCpuCommandQueue_.enqueueReadBuffer(
			oclCpuGridValueIm_, CL_TRUE /* blocking */, 0 /* offset */,
			cpuGridValueIm_.size() * sizeof(TFloat), &cpuGridValueIm_[0]);
		for (unsigned int col = cols; col < gridXZ.n1() /* cols + cpuCols */; ++col) {
			for (unsigned int row = 0; row < minRowIdx_[col]; ++row) {
				gridValue(col, row) = 0;
			}
			unsigned int gridPointIdx = firstGridPointIdx_[col] - firstGridPointIdx_[cols];
			for (unsigned int row = minRowIdx_[col]; row < gridXZ.n2(); ++row, ++gridPointIdx) {
				gridValue(col, row) = std::complex<TFloat>(cpuGridValueRe_[gridPointIdx], cpuGridValueIm_[gridPointIdx]);
			}
		}
	}
# endif
#endif

#ifdef USE_EXECUTION_TIME_MEASUREMENT
	tProcessColumn.put(processColumnTimer.getTime());
#endif

	//LOG_DEBUG << "END ========== VectorialCombinedTwoMediumImagingProcessor4::process ==========";
}



template<typename TFloat>
struct VectorialCombinedTwoMediumImagingProcessor4<TFloat>::CalculateDelays {
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
#ifdef USE_SIMD
					FermatPrinciple::findMinTimeInTwoSteps2(
							fermatBlockSize,
							invC1, invC2,
							interfacePointList,
							xArray[elem], TFloat(0), point.x, point.z,
							tMin, idxMin);
#else
					FermatPrinciple::findMinTimeInTwoSteps(
							fermatBlockSize,
							config.propagationSpeed1, config.propagationSpeed2,
							interfacePointList,
							xArray[elem], TFloat(0), point.x, point.z,
							tMin, idxMin);
#endif
					delayMatrix(col, minRowIdx[col], elem) = tMin * fs;
					lastInterfaceIdx = idxMin;
				}

				const TFloat* medium1Delays = &medium1DelayMatrix(elem, 0);

				for (unsigned int row = minRowIdx[col] + 1; row < numRows; ++row) {
					const auto& point = gridXZ(col, row);
					unsigned int idxMin = lastInterfaceIdx;
					TFloat tC2Min;
					{
						const XZ<TFloat>& ifPoint = interfacePointList[idxMin];
#ifdef USE_SIMD
						tC2Min = medium1Delays[idxMin] + SIMD::calcDistance(ifPoint.x, ifPoint.z, point.x, point.z);
#else
						const TFloat dx2 = point.x - ifPoint.x;
						const TFloat dz2 = point.z - ifPoint.z;
						tC2Min = medium1Delays[idxMin] + std::sqrt(dx2 * dx2 + dz2 * dz2);
#endif
					}
					for (unsigned int idxSearch = idxMin + 1, end = interfacePointList.size(); idxSearch < end; ++idxSearch) {
						const XZ<TFloat>& ifPoint = interfacePointList[idxSearch];
#ifdef USE_SIMD
						const TFloat tC2 = medium1Delays[idxSearch] + SIMD::calcDistance(ifPoint.x, ifPoint.z, point.x, point.z);
#else
						const TFloat dx2 = point.x - ifPoint.x;
						const TFloat dz2 = point.z - ifPoint.z;
						const TFloat tC2 = medium1Delays[idxSearch] + std::sqrt(dx2 * dx2 + dz2 * dz2);
#endif
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
#ifdef USE_SIMD
							const TFloat tC2 = medium1Delays[idxSearch] + SIMD::calcDistance(ifPoint.x, ifPoint.z, point.x, point.z);
#else
							const TFloat dx2 = point.x - ifPoint.x;
							const TFloat dz2 = point.z - ifPoint.z;
							const TFloat tC2 = medium1Delays[idxSearch] + std::sqrt(dx2 * dx2 + dz2 * dz2);
#endif
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

					delayMatrix(col, row, elem) = tC2Min * fsInvC2;
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
	Tensor3<TFloat, tbb::cache_aligned_allocator<TFloat>>& delayMatrix;
};



template<typename TFloat>
struct VectorialCombinedTwoMediumImagingProcessor4<TFloat>::PrepareDataWithOneTxElem {
	void operator()(const tbb::blocked_range<unsigned int>& r) const {
		PrepareDataThreadData& local = prepareDataTLS.local();

		for (unsigned int rxElem = r.begin(); rxElem != r.end(); ++rxElem) {
			if (upsamplingFactor > 1) {
				// Interpolate the signal.
				local.interpolator.interpolate(&acqDataList[baseElementIdx](rxElem, 0), samplesPerChannelLow, &local.signal[0]);
			} else {
				auto range = acqDataList[baseElementIdx].range2(rxElem);
				std::copy(range.begin(), range.end(), local.signal.begin());
			}

			Util::removeDC(&local.signal[0], local.signal.size());

			// Obtain the analytic signal.
			local.envelope.getAnalyticSignal(
					&local.signal[0],
					local.signal.size(),
					&signalMatrix(stepIdx, rxElem, 0));
		}
	}

	// Use only small types or references/pointers, because this object may be copied many times.
	const std::size_t samplesPerChannelLow;
	const std::vector<Matrix<TFloat>>& acqDataList;
	const unsigned int upsamplingFactor;
	const unsigned int stepIdx;
	const unsigned int baseElementIdx;
	tbb::enumerable_thread_specific<PrepareDataThreadData>& prepareDataTLS;
	Tensor3<std::complex<TFloat>, tbb::cache_aligned_allocator<std::complex<TFloat>>>& signalMatrix;
};



template<typename TFloat>
struct VectorialCombinedTwoMediumImagingProcessor4<TFloat>::ProcessColumnWithOneTxElem {
	void operator()(const tbb::blocked_range<unsigned int>& r) const {
		//LOG_DEBUG << "col = " << r.begin() << " n = " << (r.end() - r.begin());

		const unsigned int signalLength = signalMatrix.n3();
		const unsigned int maxPosition = signalLength - 2;

		for (unsigned int col = r.begin(); col != r.end(); ++col) {
			unsigned int gridPointIdx = firstGridPointIdx[col] - firstGridPointIdx[firstCol];
			for (unsigned int row = minRowIdx[col]; row < numRows; ++row, ++gridPointIdx) {
				const TFloat* delays = &delayMatrix(col, row, stepConfig.baseElem);
				const TFloat txDelay = delays[stepConfig.txElem];
				const TFloat txOffset = signalOffset + txDelay;
				const auto* p = &signalMatrix(stepConfig.baseElemIdx, 0 /* rxElem */, 0);
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
#if defined USE_OPENCL && defined VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_USE_TRANSPOSE
							rawData[gridPointIdx * rawDataN2 + rxIdx    ] = v.real();
							rawData[gridPointIdx * rawDataN2 + rxIdx + 1] = v.imag();
#else
							rawData[ rxIdx      * rawDataN2 + gridPointIdx] = v.real();
							rawData[(rxIdx + 1) * rawDataN2 + gridPointIdx] = v.imag();
#endif
						} else {
#if defined USE_OPENCL && defined VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_USE_TRANSPOSE
							rawData[gridPointIdx * rawDataN2 + rxIdx    ] = 0;
							rawData[gridPointIdx * rawDataN2 + rxIdx + 1] = 0;
#else
							rawData[ rxIdx      * rawDataN2 + gridPointIdx] = 0;
							rawData[(rxIdx + 1) * rawDataN2 + gridPointIdx] = 0;
#endif
						}
					} else {
#if defined USE_OPENCL && defined VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_4_OPENCL_USE_TRANSPOSE
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
	const Tensor3<std::complex<TFloat>, tbb::cache_aligned_allocator<std::complex<TFloat>>>& signalMatrix;
	const StepConfiguration stepConfig;
	const std::vector<unsigned int, tbb::cache_aligned_allocator<unsigned int>>& minRowIdx;
	const std::vector<unsigned int, tbb::cache_aligned_allocator<unsigned int>>& firstGridPointIdx;
	const Tensor3<TFloat, tbb::cache_aligned_allocator<TFloat>>& delayMatrix;
	TFloat* rawData;
	unsigned int rawDataN2;
};



template<typename TFloat>
struct VectorialCombinedTwoMediumImagingProcessor4<TFloat>::ProcessColumn2WithOneTxElem {
	void operator()(const tbb::blocked_range<unsigned int>& r) const {
		//LOG_DEBUG << "col = " << r.begin() << " n = " << (r.end() - r.begin());
		ProcessColumn2ThreadData& local = processColumn2TLS.local();

		local.rxSignalSumList.resize(config.numElements);

		for (unsigned int col = r.begin(); col != r.end(); ++col) {
			unsigned int gridPointIdx = firstGridPointIdx[col];
			for (unsigned int row = minRowIdx[col]; row < numRows; ++row, ++gridPointIdx) {
				for (unsigned int rxElem = 0; rxElem < config.numElements; ++rxElem) {
					const unsigned int rawDataRow = rxElem * 2;
					local.rxSignalSumList[rxElem] = std::complex<TFloat>(
										rawDataMatrix(rawDataRow    , gridPointIdx),
										rawDataMatrix(rawDataRow + 1, gridPointIdx));
				}

				std::complex<TFloat> sum = 0.0;
				for (unsigned int rxElem = 0; rxElem < config.numElements; ++rxElem) {
					sum += local.rxSignalSumList[rxElem] * rxApod[rxElem];
				}

				if (local.coherenceFactor.enabled()) {
					gridValue(col, row) += sum * local.coherenceFactor.calculate(&local.rxSignalSumList[0], local.rxSignalSumList.size());
				} else {
					gridValue(col, row) += sum;
				}
			}
		}
	}

	// Use only small types or references/pointers, because this object may be copied many times.
	const std::size_t numRows;
	const TwoMediumSTAConfiguration<TFloat>& config;
	const std::vector<TFloat>& rxApod;
	const StepConfiguration stepConfig;
	const std::vector<unsigned int, tbb::cache_aligned_allocator<unsigned int>>& minRowIdx;
	const std::vector<unsigned int, tbb::cache_aligned_allocator<unsigned int>>& firstGridPointIdx;
	const Matrix<TFloat, tbb::cache_aligned_allocator<TFloat>>& rawDataMatrix;
	tbb::enumerable_thread_specific<ProcessColumn2ThreadData>& processColumn2TLS;
	Matrix<std::complex<TFloat>>& gridValue;
};



template<typename TFloat>
std::string
VectorialCombinedTwoMediumImagingProcessor4<TFloat>::getKernel() const
{
	return R"CLC(

#define GROUP_SIZE 16

//#define USE_PRNG 1

#ifdef USE_PRNG
# define PRNG_M (2147483647.0f)
# define PRNG_INV_M (1.0f / PRNG_M)
#endif

// NVIDIA sm_12:
//   - Local (shared) memory has 16 banks.
__kernel
void
transposeKernel(
		__global float* rawData,
		__global float* rawDataT,
		unsigned int oldSizeX,
		unsigned int oldSizeY,
		__local float* temp) // (GROUP_SIZE + 1) * GROUP_SIZE
{
	unsigned int iX = get_global_id(0);
	unsigned int iY = get_global_id(1);
	if (iX < oldSizeX && iY < oldSizeY) {
		temp[get_local_id(0) + (GROUP_SIZE + 1) * get_local_id(1)] = rawData[iX + oldSizeX * iY];
	}

	barrier(CLK_LOCAL_MEM_FENCE);

	iX = get_group_id(1) * GROUP_SIZE + get_local_id(0);
	iY = get_group_id(0) * GROUP_SIZE + get_local_id(1);
	if (iX < oldSizeY && iY < oldSizeX) {
		rawDataT[iX + oldSizeY * iY] = temp[(GROUP_SIZE + 1) * get_local_id(0) + get_local_id(1)];
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

#ifdef USE_PRNG
float
minstdPRNG(int* x) // 1 <= x < m
{
	long a = 16807; // 7**5
	long m = 2147483647; // 2**31-1 (prime)
	*x = ((long) *x * a) % m;
	return *x * PRNG_INV_M;
}
#endif

float
calcPCF(float* re, float* im, float factor
#ifdef USE_PRNG
		, int* prngX
#endif
		)
{
	float phi[NUM_RX_ELEM];
	float phiAux[NUM_RX_ELEM];

#pragma unroll
	for (unsigned int i = 0; i < NUM_RX_ELEM; ++i) {
#ifdef USE_PRNG
		if (re[i] == 0.0f && im[i] == 0.0f) {
			float prn = minstdPRNG(prngX);
			phi[i] = (2.0f * prn - 1.0f) * M_PI_F;
		} else
#endif
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
		float pcfFactor
#ifdef USE_PRNG
		, __global int* prngState
#endif
		)
{
	float rxSignalListRe[NUM_RX_ELEM];
	float rxSignalListIm[NUM_RX_ELEM];

	const unsigned int point = get_global_id(0);
	if (point >= numGridPoints) return;

#ifdef USE_PRNG
	int prngX = prngState[point];
#endif
	for (unsigned int i = 0; i < NUM_RX_ELEM; ++i) {
		rxSignalListRe[i] = rawData[ (i << 1)      * numGridPoints + point];
		rxSignalListIm[i] = rawData[((i << 1) + 1) * numGridPoints + point];
	}

	float pcf = calcPCF(rxSignalListRe, rxSignalListIm, pcfFactor
#ifdef USE_PRNG
				, &prngX
#endif
				);

	float sumRe = 0.0f;
	float sumIm = 0.0f;
	for (unsigned int i = 0; i < NUM_RX_ELEM; ++i) {
		sumRe += rxSignalListRe[i] * rxApod[i];
		sumIm += rxSignalListIm[i] * rxApod[i];
	}

	gridValueRe[point] += sumRe * pcf;
	gridValueIm[point] += sumIm * pcf;

#ifdef USE_PRNG
	prngState[point] = prngX;
#endif
}

)CLC";
}

} // namespace Lab

#endif // VECTORIALCOMBINEDTWOMEDIUMIMAGINGPROCESSOR4_H
