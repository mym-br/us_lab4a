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

#include "VectorialCombinedTwoMediumImagingCUDAProcessor2.h"

#include <algorithm> /* copy */
#include <cmath> /* ceil, sqrt */
#include <cstring> /* memset */

#include <tbb/partitioner.h>
#include <tbb/tbb.h>

#include "cuda.h"

#include "ArrayGeometry.h"
#include "Exception.h"
#include "FermatPrinciple.h"
#include "Geometry.h"
#include "Log.h"
#include "Timer.h"
#include "Util.h"

#include "CUDAUtil.h"



// Depends on the signal.
// 1.0 --> pi radian / sample at the original sampling rate.
#define UPSAMP_FILTER_HALF_TRANSITION_WIDTH (0.2)

#define USE_TRANSPOSE 1

#define BLOCK_SIZE 64

#ifndef MFloat
# define MFloat float
#endif



namespace Lab {

// NVIDIA sm_50 or newer:
//   - Shared memory has 32 banks of 32 bits.

#define TRANSP_BLOCK_SIZE 16
#define NUM_RX_ELEM 32

// Defined in VectorialCombinedTwoMediumImagingCUDAProcessor.cu.
extern __global__ void transposeKernel(float* rawData, float* rawDataT, int oldSizeX, int oldSizeY);
extern __device__ float calcPCF(float* re, float* im, float factor);

__global__
void
processColumnWithOneTxElemKernel(
		unsigned int numCols,
		unsigned int numRows,
		unsigned int firstCol,
		unsigned int numElements,
		float signalOffset,
		float (*signalTensor)[2],
		unsigned int signalTensorN2,
		unsigned int signalTensorN3,
		unsigned int baseElem,
		unsigned int baseElemIdx,
		unsigned int txElem,
		const unsigned int* minRowIdx,
		const unsigned int* firstGridPointIdx,
		const float* delayTensor,
		unsigned int delayTensorN2,
		unsigned int delayTensorN3,
		float* rawData,
		unsigned int rawDataN2) {

	const unsigned int signalLength = signalTensorN3;
	const unsigned int maxPosition = signalLength - 2;

	const unsigned int col = blockIdx.y * blockDim.y + threadIdx.y;
	if (col >= numCols) return;

	const unsigned int rxElem = blockIdx.x * blockDim.x + threadIdx.x;
	if (rxElem >= numElements) return;

	unsigned int gridPointIdx = firstGridPointIdx[col] - firstGridPointIdx[firstCol];
	for (unsigned int row = minRowIdx[col]; row < numRows; ++row, ++gridPointIdx) {
		const float* delays = delayTensor + (col * delayTensorN2 + row) * delayTensorN3 + baseElem;

		const float txDelay = delays[txElem];
		const float txOffset = signalOffset + txDelay;
		const float (*p)[2] = signalTensor + baseElemIdx * signalTensorN2 * signalTensorN3
					+ rxElem * signalLength;
		const unsigned int rxIdx = rxElem * 2;

		// Linear interpolation.
		const float position = txOffset + delays[rxElem];
		if (position >= 0.0f) {
			const unsigned int positionIdx = static_cast<unsigned int>(position);
			if (positionIdx <= maxPosition) {
				const float k = position - positionIdx;
				float v0[2]; // complex
				v0[0] = p[positionIdx][0];
				v0[1] = p[positionIdx][1];
				float v1[2]; // complex
				v1[0] = p[positionIdx + 1][0];
				v1[1] = p[positionIdx + 1][1];
				float v[2]; // complex
				v[0] = v0[0] + k * (v1[0] - v0[0]);
				v[1] = v0[1] + k * (v1[1] - v0[1]);
#ifdef USE_TRANSPOSE
				rawData[gridPointIdx * rawDataN2 + rxIdx    ] = v[0];
				rawData[gridPointIdx * rawDataN2 + rxIdx + 1] = v[1];
#else
				rawData[ rxIdx      * rawDataN2 + gridPointIdx] = v[0];
				rawData[(rxIdx + 1) * rawDataN2 + gridPointIdx] = v[1];
#endif
				continue;
			}
		}
#ifdef USE_TRANSPOSE
		rawData[gridPointIdx * rawDataN2 + rxIdx    ] = 0;
		rawData[gridPointIdx * rawDataN2 + rxIdx + 1] = 0;
#else
		rawData[ rxIdx      * rawDataN2 + gridPointIdx] = 0;
		rawData[(rxIdx + 1) * rawDataN2 + gridPointIdx] = 0;
#endif
	}
}

__global__
void
processImageKernel2(
		float* rawData,
		int numGridPoints,
		float (*gridValue)[2],
		float* rxApod)
{
	float rxSignalListRe[NUM_RX_ELEM];
	float rxSignalListIm[NUM_RX_ELEM];

	const int point = blockIdx.x * blockDim.x + threadIdx.x;
	if (point >= numGridPoints) return;

	for (unsigned int i = 0; i < NUM_RX_ELEM; ++i) {
		rxSignalListRe[i] = rawData[ (i << 1)      * numGridPoints + point];
		rxSignalListIm[i] = rawData[((i << 1) + 1) * numGridPoints + point];
		//rxSignalListRe[i] = rawData[point * (NUM_RX_ELEM << 1) + (i << 1)];
		//rxSignalListIm[i] = rawData[point * (NUM_RX_ELEM << 1) + ((i << 1) + 1)];
	}

	float sumRe = 0.0;
	float sumIm = 0.0;
	for (unsigned int i = 0; i < NUM_RX_ELEM; ++i) {
		sumRe += rxSignalListRe[i] * rxApod[i];
		sumIm += rxSignalListIm[i] * rxApod[i];
	}

	gridValue[point][0] += sumRe;
	gridValue[point][1] += sumIm;
}

__global__
void
processImagePCFKernel2(
		float* rawData,
		int numGridPoints,
		float (*gridValue)[2],
		float* rxApod,
		float pcfFactor)
{
	float rxSignalListRe[NUM_RX_ELEM];
	float rxSignalListIm[NUM_RX_ELEM];

	const int point = blockIdx.x * blockDim.x + threadIdx.x;
	if (point >= numGridPoints) return;

	for (int i = 0; i < NUM_RX_ELEM; ++i) {
		rxSignalListRe[i] = rawData[ (i << 1)      * numGridPoints + point];
		rxSignalListIm[i] = rawData[((i << 1) + 1) * numGridPoints + point];
	}

	float pcf = calcPCF(rxSignalListRe, rxSignalListIm, pcfFactor);

	float sumRe = 0.0;
	float sumIm = 0.0;
	for (int i = 0; i < NUM_RX_ELEM; ++i) {
		sumRe += rxSignalListRe[i] * rxApod[i];
		sumIm += rxSignalListIm[i] * rxApod[i];
	}

	gridValue[point][0] += sumRe * pcf;
	gridValue[point][1] += sumIm * pcf;
}

//=============================================================================

struct VectorialCombinedTwoMediumImagingCUDAProcessor2Data {
	bool cudaDataInitialized;
	CUDADevMem<MFloat> rxApod;
	CUDADevMem<MFloat[2]> signalTensor;
	CUDADevMem<unsigned int> firstGridPointIdx;
	CUDADevMem<unsigned int> minRowIdx;
	CUDADevMem<MFloat> delayTensor;
	CUDADevMem<MFloat> rawData;
#ifdef USE_TRANSPOSE
	CUDADevMem<MFloat> rawDataT;
#endif
	CUDAHostDevMem<MFloat[2]> gridValue;

	VectorialCombinedTwoMediumImagingCUDAProcessor2Data()
		: cudaDataInitialized()
	{}
};

template<typename TFloat>
struct CalculateDelays2 {
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
struct PrepareDataWithOneTxElem2 {
	void operator()(const tbb::blocked_range<unsigned int>& r) const {
		typename VectorialCombinedTwoMediumImagingCUDAProcessor2::PrepareDataThreadData<TFloat>& local = prepareDataTLS.local();

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
	tbb::enumerable_thread_specific<typename VectorialCombinedTwoMediumImagingCUDAProcessor2::PrepareDataThreadData<TFloat>>& prepareDataTLS;
	Tensor3<std::complex<TFloat>, tbb::cache_aligned_allocator<std::complex<TFloat>>>& signalTensor;
};

VectorialCombinedTwoMediumImagingCUDAProcessor2::VectorialCombinedTwoMediumImagingCUDAProcessor2(
			const TwoMediumSTAConfiguration<MFloat>& config,
			std::vector<Matrix<MFloat>>& acqDataList,
			unsigned int upsamplingFactor,
			AnalyticSignalCoherenceFactorProcessor<MFloat>& coherenceFactor,
			MFloat maxFermatBlockSize,
			MFloat peakOffset,
			unsigned int signalStartOffset)
		: config_(config)
		, acqDataList_(acqDataList)
		, upsamplingFactor_(upsamplingFactor)
		, coherenceFactor_(coherenceFactor)
		, maxFermatBlockSize_(maxFermatBlockSize)
		, lambda2_(config_.propagationSpeed2 / config_.centerFrequency)
		, rawDataN1_()
		, rawDataN2_()
		, rawDataSize_()
{
	const std::size_t origSignalLength = acqDataList_[0].n2();

	signalOffset_ = (config_.samplingFrequency * upsamplingFactor_) * (peakOffset / config_.centerFrequency) - signalStartOffset * upsamplingFactor_;
	signalLength_ = origSignalLength * upsamplingFactor_;
	LOG_DEBUG << "signalOffset_: " << signalOffset_ << " signalLength_: " << signalLength_;

	PrepareDataThreadData<MFloat> prepareDataThreadData;
	if (upsamplingFactor_ > 1) {
		prepareDataThreadData.interpolator.prepare(upsamplingFactor_, UPSAMP_FILTER_HALF_TRANSITION_WIDTH);
	}
	prepareDataThreadData.signal.resize(signalLength_);
	prepareDataTLS_ = std::make_unique<tbb::enumerable_thread_specific<PrepareDataThreadData<MFloat>>>(prepareDataThreadData);

	data_ = std::make_unique<VectorialCombinedTwoMediumImagingCUDAProcessor2Data>();

	if (Log::isDebugEnabled()) {
		int device;
		exec(cudaGetDevice(&device));

		cudaDeviceProp prop;
		exec(cudaGetDeviceProperties(&prop, device));
		LOG_DEBUG << "CUDA device: " << prop.name;
	}
}

// This destructor can't be inline.
VectorialCombinedTwoMediumImagingCUDAProcessor2::~VectorialCombinedTwoMediumImagingCUDAProcessor2()
{
}

void
VectorialCombinedTwoMediumImagingCUDAProcessor2::process(
						const std::vector<StepConfiguration>& stepConfigList,
						const std::vector<XZ<MFloat>>& interfacePointList,
						const std::vector<MFloat>& rxApod,
						const Matrix<XZ<MFloat>>& gridXZ,
						Matrix<std::complex<MFloat>>& gridValue)
{
	//LOG_DEBUG << "BEGIN ========== VectorialCombinedTwoMediumImagingCUDAProcessor2::process ==========";

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

	XZ<MFloat> p1 = interfacePointList[0];
	XZ<MFloat> p2 = interfacePointList[1];
	const MFloat dx = p2.x - p1.x;
	const MFloat dz = p2.z - p1.z;
	const MFloat r = std::sqrt(dx * dx + dz * dz);
	const unsigned int fermatBlockSize = FermatPrinciple::calcBlockSizeForTwoStepSearch(interfacePointList.size(), r, lambda2_, maxFermatBlockSize_);
	LOG_DEBUG << "fermatBlockSize: " << fermatBlockSize;

#ifdef USE_EXECUTION_TIME_MEASUREMENT
	Timer minRowIdxTimer;
#endif
	//==================================================
	// Find the z-coordinate of the interface in
	// each column (minimum row).
	//==================================================
	const MFloat zStepGrid = gridXZ(0, 1).z - gridXZ(0, 0).z;
	unsigned int gridPointIdx = 0;
	for (unsigned int col = 0; col < gridXZ.n1(); ++col) {
		auto& point = gridXZ(col, 0);
		firstGridPointIdx_[col] = gridPointIdx; // the points below the interface are not considered

		// Find the z coordinate of the interface.
		MFloat zIdxMin;
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

#ifdef USE_TRANSPOSE
	const std::size_t transpNumGridPoints = CUDAUtil::roundUpToMultipleOfBlockSize(numGridPoints, TRANSP_BLOCK_SIZE);
	LOG_DEBUG << "numGridPoints: " << numGridPoints << " transpNumGridPoints: " << transpNumGridPoints;
	rawDataN1_ = transpNumGridPoints;
	rawDataN2_ = 2 * config_.numElements /* real, imag */;
#else
	rawDataN1_ = 2 * config_.numElements /* real, imag */;
	rawDataN2_ = numGridPoints;
#endif
	rawDataSize_ = rawDataN1_ * rawDataN2_;

	if (!data_->cudaDataInitialized) {
		data_->rxApod            = CUDADevMem<MFloat>(rxApod.size());
		data_->signalTensor      = CUDADevMem<MFloat[2]>(signalTensor_.n1() * signalTensor_.n2() * signalTensor_.n3());
		data_->firstGridPointIdx = CUDADevMem<unsigned int>(firstGridPointIdx_.size());
		data_->minRowIdx         = CUDADevMem<unsigned int>(minRowIdx_.size());
		data_->delayTensor       = CUDADevMem<MFloat>(delayTensor_.n1() * delayTensor_.n2() * delayTensor_.n3());
		data_->rawData           = CUDADevMem<MFloat>(rawDataSize_);
#ifdef USE_TRANSPOSE
		data_->rawDataT          = CUDADevMem<MFloat>(rawDataSize_);
#endif
		data_->gridValue         = CUDAHostDevMem<MFloat[2]>(numGridPoints);

		data_->cudaDataInitialized = true;
	}

	// Prepare buffers.
	exec(cudaMemset(data_->gridValue.devPtr, 0, data_->gridValue.sizeInBytes));
	exec(cudaMemcpy(data_->rxApod.devPtr, rxApod.data(), Util::sizeInBytes(rxApod), cudaMemcpyHostToDevice));

	const MFloat c2ByC1 = config_.propagationSpeed2 / config_.propagationSpeed1;

	ArrayGeometry::getElementX2D(config_.numElementsMux, config_.pitch, MFloat(0) /* offset */, xArray_);

#ifdef USE_EXECUTION_TIME_MEASUREMENT
	Timer medium1DelayMatrixTimer;
#endif
	for (unsigned int elem = 0; elem < config_.numElementsMux; ++elem) {
		MFloat* delays = &medium1DelayMatrix_(elem, 0);
		for (unsigned int i = 0; i < interfacePointList.size(); ++i) {
			const XZ<MFloat>& ifPoint = interfacePointList[i];
			delays[i] = Geometry::distance2DY0(xArray_[elem], ifPoint.x, ifPoint.z) * c2ByC1;
		}
	}
#ifdef USE_EXECUTION_TIME_MEASUREMENT
	tMedium1DelayMatrix.put(medium1DelayMatrixTimer.getTime());

	Timer calculateDelaysTimer;
#endif
	CalculateDelays2<MFloat> calculateDelaysOp = {
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
		PrepareDataWithOneTxElem2<MFloat> prepareDataOp = {
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
	std::size_t procImageKernelGlobalSize = CUDAUtil::roundUpToMultipleOfBlockSize(numGridPoints, BLOCK_SIZE);
	LOG_DEBUG << numGridPoints << ':' << procImageKernelGlobalSize << ':' << BLOCK_SIZE;

	exec(cudaMemcpy(data_->firstGridPointIdx.devPtr, firstGridPointIdx_.data(),
				Util::sizeInBytes(firstGridPointIdx_), cudaMemcpyHostToDevice));
	exec(cudaMemcpy(data_->minRowIdx.devPtr        ,         minRowIdx_.data(),
				Util::sizeInBytes(minRowIdx_)        , cudaMemcpyHostToDevice));
	exec(cudaMemcpy(data_->signalTensor.devPtr     ,      signalTensor_.data(),
				Util::sizeInBytes(signalTensor_)     , cudaMemcpyHostToDevice));
	exec(cudaMemcpy(data_->delayTensor.devPtr      ,       delayTensor_.data(),
				Util::sizeInBytes(delayTensor_)      , cudaMemcpyHostToDevice));

	//==================================================
	// Step configuration loop.
	//==================================================
	for (unsigned int i = 0; i < stepConfigList.size(); ++i) {
		const auto& stepConfig = stepConfigList[i];
		LOG_DEBUG << "stepConfig.baseElemIdx: " << stepConfig.baseElemIdx;

		exec(cudaDeviceSynchronize());

		//==================================================
		// Delay and store.
		//==================================================
		{
			Timer delayStoreTimer;

			// xBlockSize * yBlockSize must be <= maximum number of threads per block.
			const std::size_t xBlockSize = 32;
			const std::size_t yBlockSize = 32;
			const std::size_t xNumBlocks = CUDAUtil::numberOfBlocks(config_.numElements, xBlockSize);
			const std::size_t yNumBlocks = CUDAUtil::numberOfBlocks(gridXZ.n1(), yBlockSize);
			const dim3 gridDim(xNumBlocks, yNumBlocks);
			const dim3 blockDim(xBlockSize, yBlockSize);

			processColumnWithOneTxElemKernel<<<gridDim, blockDim>>>(
					gridXZ.n1(),
					gridXZ.n2(),
					0,
					config_.numElements,
					signalOffset_,
					data_->signalTensor.devPtr,
					signalTensor_.n2(),
					signalTensor_.n3(),
					stepConfig.baseElem,
					stepConfig.baseElemIdx,
					stepConfig.txElem,
					data_->minRowIdx.devPtr,
					data_->firstGridPointIdx.devPtr,
					data_->delayTensor.devPtr,
					delayTensor_.n2(),
					delayTensor_.n3(),
					data_->rawData.devPtr,
					rawDataN2_);
			checkKernelLaunchError();

			LOG_DEBUG << "DELAY-STORE " << delayStoreTimer.getTime();
		}

#ifdef USE_TRANSPOSE
		{
			const dim3 gridDim(rawDataN2_ / TRANSP_BLOCK_SIZE, rawDataN1_ / TRANSP_BLOCK_SIZE);
			const dim3 blockDim(TRANSP_BLOCK_SIZE, TRANSP_BLOCK_SIZE);

			transposeKernel<<<gridDim, blockDim>>>(
							data_->rawData.devPtr,
							data_->rawDataT.devPtr,
							rawDataN2_,
							rawDataN1_);
			checkKernelLaunchError();
		}
#endif
		if (coherenceFactor_.enabled()) {
			std::vector<MFloat> cfConstants;
			coherenceFactor_.implementation().getConstants(cfConstants);

			processImagePCFKernel2<<<procImageKernelGlobalSize / BLOCK_SIZE, BLOCK_SIZE>>>(
#ifdef USE_TRANSPOSE
							data_->rawDataT.devPtr,
							rawDataN1_,
#else
							data_->rawData.devPtr,
							rawDataN2_,
#endif
							data_->gridValue.devPtr,
							data_->rxApod.devPtr,
							cfConstants[2] /* factor */);
			checkKernelLaunchError();
		} else {
			processImageKernel2<<<procImageKernelGlobalSize / BLOCK_SIZE, BLOCK_SIZE>>>(
#ifdef USE_TRANSPOSE
							data_->rawDataT.devPtr,
							rawDataN1_,
#else
							data_->rawData.devPtr,
							rawDataN2_,
#endif
							data_->gridValue.devPtr,
							data_->rxApod.devPtr);
			checkKernelLaunchError();
		}
	}

	exec(data_->gridValue.copyDeviceToHost());

	//==================================================
	// Read the formed image.
	//==================================================
	for (unsigned int col = 0; col < cols; ++col) {
		for (unsigned int row = 0; row < minRowIdx_[col]; ++row) {
			gridValue(col, row) = 0;
		}
		unsigned int gridPointIdx = firstGridPointIdx_[col];
		for (unsigned int row = minRowIdx_[col]; row < gridXZ.n2(); ++row, ++gridPointIdx) {
			gridValue(col, row) = std::complex<MFloat>(
							data_->gridValue.hostPtr[gridPointIdx][0],
							data_->gridValue.hostPtr[gridPointIdx][1]);
		}
	}

#ifdef USE_EXECUTION_TIME_MEASUREMENT
	tProcessColumn.put(processColumnTimer.getTime());
#endif
	//LOG_DEBUG << "END ========== VectorialCombinedTwoMediumImagingCUDAProcessor2::process ==========";
}

} // namespace Lab
