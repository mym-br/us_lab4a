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

#include "VectorialCombinedTwoMediumImagingCUDAProcessor.h"

#include <algorithm> /* copy */
#include <cmath> /* ceil, sqrt */
#include <cstring> /* memset */

#include <tbb/blocked_range.h>
#include <tbb/cache_aligned_allocator.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>
#include <tbb/partitioner.h>

#include "cuda.h"

#include "ArrayGeometry.h"
#include "Exception.h"
#include "FermatPrinciple.h"
#include "Geometry.h"
#include "Log.h"
#include "Timer.h"
#include "Util.h"

#include "CUDACoherenceFactor.cuh"
#include "CUDAUtil.h"



// Depends on the signal.
// 1.0 --> pi radian / sample at the original sampling rate.
#define UPSAMP_FILTER_HALF_TRANSITION_WIDTH (0.2)

#define NUM_RAW_DATA_BUFFERS 2
#define USE_TRANSPOSE 1

#define BLOCK_SIZE 64

#ifndef MFloat
# define MFloat float
#endif



namespace Lab {

#define TRANSP_BLOCK_SIZE 16
#define NUM_RX_ELEM 32

// NVIDIA sm_50 or newer:
//   - Shared memory has 32 banks of 32 bits.
template<typename TFloat>
__global__
void
transposeKernel(
		TFloat* rawData,
		TFloat* rawDataT,
		unsigned int oldSizeX,
		unsigned int oldSizeY)
{
	__shared__ TFloat temp[TRANSP_BLOCK_SIZE][TRANSP_BLOCK_SIZE + 1]; // +1 to avoid bank conflicts

	unsigned int iX = blockIdx.x * TRANSP_BLOCK_SIZE + threadIdx.x;
	unsigned int iY = blockIdx.y * TRANSP_BLOCK_SIZE + threadIdx.y;

	if (iX < oldSizeX && iY < oldSizeY) {
		temp[threadIdx.x][threadIdx.y] = rawData[iY * oldSizeX + iX];
	}

	__syncthreads();

	iX = blockIdx.y * TRANSP_BLOCK_SIZE + threadIdx.x;
	iY = blockIdx.x * TRANSP_BLOCK_SIZE + threadIdx.y;

	if (iX < oldSizeY && iY < oldSizeX) {
		rawDataT[iY * oldSizeY + iX] = temp[threadIdx.y][threadIdx.x];
	}
}

template<typename TFloat>
__global__
void
processImageKernel(
		TFloat* rawData,
		unsigned int numGridPoints,
		TFloat* gridValueRe,
		TFloat* gridValueIm,
		TFloat* rxApod)
{
	TFloat rxSignalListRe[NUM_RX_ELEM];
	TFloat rxSignalListIm[NUM_RX_ELEM];

	const unsigned int point = blockIdx.x * blockDim.x + threadIdx.x;
	if (point >= numGridPoints) return;

	unsigned int idx1 = point;
	unsigned int idx2 = point + numGridPoints;
	const unsigned int step = 2 * numGridPoints;
	for (int i = 0; i < NUM_RX_ELEM; ++i, idx1 += step, idx2 += step) {
		rxSignalListRe[i] = rawData[idx1];
		rxSignalListIm[i] = rawData[idx2];
	}

	TFloat sumRe = 0;
	TFloat sumIm = 0;
	for (int i = 0; i < NUM_RX_ELEM; ++i) {
		sumRe += rxSignalListRe[i] * rxApod[i];
		sumIm += rxSignalListIm[i] * rxApod[i];
	}

	gridValueRe[point] += sumRe;
	gridValueIm[point] += sumIm;
}

template<typename TFloat>
__global__
void
processImagePCFKernel(
		TFloat* rawData,
		unsigned int numGridPoints,
		TFloat* gridValueRe,
		TFloat* gridValueIm,
		TFloat* rxApod,
		TFloat pcfFactor)
{
	TFloat rxSignalListRe[NUM_RX_ELEM];
	TFloat rxSignalListIm[NUM_RX_ELEM];

	const unsigned int point = blockIdx.x * blockDim.x + threadIdx.x;
	if (point >= numGridPoints) return;

	unsigned int idx1 = point;
	unsigned int idx2 = point + numGridPoints;
	const unsigned int step = 2 * numGridPoints;
	for (int i = 0; i < NUM_RX_ELEM; ++i, idx1 += step, idx2 += step) {
		rxSignalListRe[i] = rawData[idx1];
		rxSignalListIm[i] = rawData[idx2];
	}

	const TFloat pcf = calcPCF<TFloat, NUM_RX_ELEM>(rxSignalListRe, rxSignalListIm, pcfFactor);

	TFloat sumRe = 0.0;
	TFloat sumIm = 0.0;
	for (int i = 0; i < NUM_RX_ELEM; ++i) {
		sumRe += rxSignalListRe[i] * rxApod[i];
		sumIm += rxSignalListIm[i] * rxApod[i];
	}

	gridValueRe[point] += sumRe * pcf;
	gridValueIm[point] += sumIm * pcf;
}

//=============================================================================

struct VectorialCombinedTwoMediumImagingCUDAProcessor::CUDAData {
	CUDAContext c; // must be the first member
	bool cudaDataInitialized;
	CUDADevMem<MFloat> rxApod;
	std::vector<CUDAHostDevMem<MFloat>> rawDataList;
#ifdef USE_TRANSPOSE
	CUDADevMem<MFloat> rawDataT;
#endif
	CUDAHostDevMem<MFloat> gridValueRe;
	CUDAHostDevMem<MFloat> gridValueIm;

	CUDAData()
		: cudaDataInitialized()
		, rawDataList(NUM_RAW_DATA_BUFFERS)
	{}
};

template<typename TFloat>
struct VectorialCombinedTwoMediumImagingCUDAProcessor::CalculateDelays {
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
struct VectorialCombinedTwoMediumImagingCUDAProcessor::PrepareDataWithOneTxElem {
	void operator()(const tbb::blocked_range<unsigned int>& r) const {
		typename VectorialCombinedTwoMediumImagingCUDAProcessor::PrepareDataThreadData<TFloat>& local = prepareDataTLS.local();

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
	tbb::enumerable_thread_specific<typename VectorialCombinedTwoMediumImagingCUDAProcessor::PrepareDataThreadData<TFloat>>& prepareDataTLS;
	Tensor3<std::complex<TFloat>, tbb::cache_aligned_allocator<std::complex<TFloat>>>& signalTensor;
};

template<typename TFloat>
struct VectorialCombinedTwoMediumImagingCUDAProcessor::ProcessColumnWithOneTxElem {
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
#ifdef USE_TRANSPOSE
							rawData[gridPointIdx * rawDataN2 + rxIdx    ] = v.real();
							rawData[gridPointIdx * rawDataN2 + rxIdx + 1] = v.imag();
#else
							rawData[ rxIdx      * rawDataN2 + gridPointIdx] = v.real();
							rawData[(rxIdx + 1) * rawDataN2 + gridPointIdx] = v.imag();
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
		}
	}

	// Use only small types or references/pointers, because this object may be copied many times.
	const unsigned int numRows;
	const unsigned int firstCol;
	const TwoMediumSTAConfiguration<TFloat>& config;
	const TFloat signalOffset;
	const Tensor3<std::complex<TFloat>, tbb::cache_aligned_allocator<std::complex<TFloat>>>& signalTensor;
	const VectorialCombinedTwoMediumImagingCUDAProcessor::StepConfiguration stepConfig;
	const std::vector<unsigned int, tbb::cache_aligned_allocator<unsigned int>>& minRowIdx;
	const std::vector<unsigned int, tbb::cache_aligned_allocator<unsigned int>>& firstGridPointIdx;
	const Tensor3<TFloat, tbb::cache_aligned_allocator<TFloat>>& delayTensor;
	TFloat* rawData;
	unsigned int rawDataN2;
};

VectorialCombinedTwoMediumImagingCUDAProcessor::VectorialCombinedTwoMediumImagingCUDAProcessor(
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
		, rxApodSize_()
		, gridXZN1_()
		, gridXZN2_()
{
	if (config_.numElements != NUM_RX_ELEM) {
		THROW_EXCEPTION(InvalidParameterException, "The number of elements (" << config_.numElements <<
				") is not equal to " << NUM_RX_ELEM << '.');
	}
	if (NUM_RAW_DATA_BUFFERS < 1) {
		THROW_EXCEPTION(InvalidValueException, "The number of raw data buffers must be >= 1.");
	}

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

	data_ = std::make_unique<CUDAData>();

	if (Log::isDebugEnabled()) {
		int device;
		exec(cudaGetDevice(&device));

		cudaDeviceProp prop;
		exec(cudaGetDeviceProperties(&prop, device));
		LOG_DEBUG << "CUDA device: " << prop.name;
	}
}

// This destructor can't be inline.
VectorialCombinedTwoMediumImagingCUDAProcessor::~VectorialCombinedTwoMediumImagingCUDAProcessor()
{
}

void
VectorialCombinedTwoMediumImagingCUDAProcessor::process(
						const std::vector<StepConfiguration>& stepConfigList,
						const std::vector<XZ<MFloat>>& interfacePointList,
						const std::vector<MFloat>& rxApod,
						const Matrix<XZ<MFloat>>& gridXZ,
						Matrix<std::complex<MFloat>>& gridValue)
{
	//LOG_DEBUG << "BEGIN ========== VectorialCombinedTwoMediumImagingCUDAProcessor::process ==========";

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

	XZ<MFloat> p1 = interfacePointList[0];
	XZ<MFloat> p2 = interfacePointList[1];
	const MFloat dx = p2.x - p1.x;
	const MFloat dz = p2.z - p1.z;
	const MFloat r = std::sqrt(dx * dx + dz * dz);
	const unsigned int fermatBlockSize = FermatPrinciple::calcBlockSizeForTwoStepSearch(interfacePointList.size(), r, lambda2_, maxFermatBlockSize_);
	LOG_DEBUG << "fermatBlockSize: " << fermatBlockSize;

#ifdef LAB_ENABLE_EXECUTION_TIME_MEASUREMENT
	Timer minRowIdxTimer;
#endif
	//==================================================
	// Find the z-coordinate of the interface in
	// each column (minimum row).
	//==================================================
	const MFloat zStepGrid = gridXZ(0, 1).z - gridXZ(0, 0).z;
	unsigned int gridPointIdx = 0;
	for (unsigned int col = 0; col < numCols; ++col) {
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
	const std::size_t numGridPoints = gridPointIdx;
#ifdef LAB_ENABLE_EXECUTION_TIME_MEASUREMENT
	tMinRowIdxML.put(minRowIdxTimer.getTime());
#endif

#ifdef USE_TRANSPOSE
	const std::size_t transpNumGridPoints = CUDAUtil::roundUpToMultipleOfBlockSize(numGridPoints, TRANSP_BLOCK_SIZE);
	LOG_DEBUG << "transpNumGridPoints: " << transpNumGridPoints;
	rawDataN1_ = transpNumGridPoints;
	rawDataN2_ = 2 * config_.numElements /* real, imag */;
#else
	rawDataN1_ = 2 * config_.numElements /* real, imag */;
	rawDataN2_ = numGridPoints;
#endif
	rawDataSize_ = rawDataN1_ * rawDataN2_;

	if (!data_->cudaDataInitialized) {
		rxApodSize_ = rxApod.size();
		gridXZN1_   = gridXZ.n1();
		gridXZN2_   = gridXZ.n2();
#ifdef USE_TRANSPOSE
		const std::size_t reservedGridPoints = CUDAUtil::roundUpToMultipleOfBlockSize(numCols * numRows, TRANSP_BLOCK_SIZE);
#else
		const std::size_t reservedGridPoints = numCols * numRows;
#endif
		LOG_DEBUG << "reservedGridPoints: " << reservedGridPoints;

		data_->rxApod = CUDADevMem<MFloat>(rxApodSize_);
		for (auto& m : data_->rawDataList) {
			m = CUDAHostDevMem<MFloat>(reservedGridPoints * 2 * config_.numElements);
		}
#ifdef USE_TRANSPOSE
		data_->rawDataT = CUDADevMem<MFloat>(reservedGridPoints * 2 * config_.numElements);
#endif
		data_->gridValueRe = CUDAHostDevMem<MFloat>(numCols * numRows);
		data_->gridValueIm = CUDAHostDevMem<MFloat>(numCols * numRows);

		data_->cudaDataInitialized = true;
	} else {
		if (rxApodSize_ != rxApod.size() ||
		    gridXZN1_   != gridXZ.n1()   ||
		    gridXZN2_   != gridXZ.n2()) {
			THROW_EXCEPTION(InvalidParameterException, "Data size mismatch.");
		}
	}

	// Prepare buffers.
	exec(cudaMemset(data_->gridValueRe.devPtr, 0, data_->gridValueRe.sizeInBytes));
	exec(cudaMemset(data_->gridValueIm.devPtr, 0, data_->gridValueIm.sizeInBytes));
	exec(cudaMemcpy(data_->rxApod.devPtr, rxApod.data(), Util::sizeInBytes(rxApod), cudaMemcpyHostToDevice));

	const MFloat c2ByC1 = config_.propagationSpeed2 / config_.propagationSpeed1;

	if (xArray_.empty()) {
		ArrayGeometry::getElementX2D(config_.numElementsMux, config_.pitch, MFloat(0) /* offset */, xArray_);
	}

#ifdef LAB_ENABLE_EXECUTION_TIME_MEASUREMENT
	Timer medium1DelayMatrixTimer;
#endif
	for (unsigned int elem = 0; elem < config_.numElementsMux; ++elem) {
		MFloat* delays = &medium1DelayMatrix_(elem, 0);
		for (unsigned int i = 0; i < interfacePointList.size(); ++i) {
			const XZ<MFloat>& ifPoint = interfacePointList[i];
			delays[i] = Geometry::distance2DY0(xArray_[elem], ifPoint.x, ifPoint.z) * c2ByC1;
		}
	}
#ifdef LAB_ENABLE_EXECUTION_TIME_MEASUREMENT
	tMedium1DelayMatrixML.put(medium1DelayMatrixTimer.getTime());

	Timer calculateDelaysTimer;
#endif
	CalculateDelays<MFloat> calculateDelaysOp = {
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
#ifdef LAB_ENABLE_EXECUTION_TIME_MEASUREMENT
	tCalculateDelaysML.put(calculateDelaysTimer.getTime());
#endif
	// Only one transmit element.
#ifdef LAB_ENABLE_EXECUTION_TIME_MEASUREMENT
	Timer prepareDataTimer;
#endif
	unsigned int stepIdx = 0;
	for (const auto& stepConfig : stepConfigList) {
		PrepareDataWithOneTxElem<MFloat> prepareDataOp = {
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
#ifdef LAB_ENABLE_EXECUTION_TIME_MEASUREMENT
	tPrepareDataML.put(prepareDataTimer.getTime());

	Timer processColumnTimer;
#endif
	std::size_t procImageKernelGlobalSize = CUDAUtil::roundUpToMultipleOfBlockSize(numGridPoints, BLOCK_SIZE);
	LOG_DEBUG << "procImageKernelGlobalSize: " << procImageKernelGlobalSize << " BLOCK_SIZE: " << BLOCK_SIZE;

	//==================================================
	// Step configuration loop.
	//==================================================
	for (unsigned int i = 0; i < stepConfigList.size(); ++i) {
		const unsigned int rawBufferIdx = i % data_->rawDataList.size();
		const auto& stepConfig = stepConfigList[i];
		LOG_DEBUG << "stepConfig.baseElemIdx: " << stepConfig.baseElemIdx << " rawBufferIdx: " << rawBufferIdx;

		if (data_->rawDataList.size() == 1) {
			exec(cudaDeviceSynchronize());
		}

		Timer delayStoreTimer;

		//==================================================
		// Delay and store.
		//==================================================
		ProcessColumnWithOneTxElem<MFloat> processColumnOp = {
			static_cast<unsigned int>(numRows),
			0,
			config_,
			signalOffset_,
			signalTensor_,
			stepConfig,
			minRowIdx_,
			firstGridPointIdx_,
			delayTensor_,
			data_->rawDataList[rawBufferIdx].hostPtr,
			rawDataN2_,
		};

		//tbb::parallel_for(tbb::blocked_range<unsigned int>(0, numCols), processColumnOp);
		tbb::parallel_for(tbb::blocked_range<unsigned int>(0, numCols, 1 /* grain size */), processColumnOp, tbb::simple_partitioner());
		//processColumnOp(tbb::blocked_range<unsigned int>(0, numCols)); // single-thread
		LOG_DEBUG << "DELAY-STORE " << delayStoreTimer.getTime();

		if (data_->rawDataList.size() > 1) {
			exec(cudaDeviceSynchronize());
		}

		exec(data_->rawDataList[rawBufferIdx].copyHostToDeviceAsync(rawDataSize_ * sizeof(MFloat)));

#ifdef USE_TRANSPOSE
		{
			//Timer transposeTimer;

			dim3 gridDim(rawDataN2_ / TRANSP_BLOCK_SIZE, rawDataN1_ / TRANSP_BLOCK_SIZE);
			dim3 blockDim(TRANSP_BLOCK_SIZE, TRANSP_BLOCK_SIZE);

			transposeKernel<MFloat><<<gridDim, blockDim>>>(
							data_->rawDataList[rawBufferIdx].devPtr,
							data_->rawDataT.devPtr,
							rawDataN2_,
							rawDataN1_);
			checkKernelLaunchError();

			//exec(cudaDeviceSynchronize());
			//LOG_DEBUG << "TRANSPOSE " << transposeTimer.getTime();
		}
#endif
		if (coherenceFactor_.enabled()) {
			std::vector<MFloat> cfConstants;
			coherenceFactor_.implementation().getConstants(cfConstants);

			processImagePCFKernel<MFloat><<<procImageKernelGlobalSize / BLOCK_SIZE, BLOCK_SIZE>>>(
#ifdef USE_TRANSPOSE
							data_->rawDataT.devPtr,
							rawDataN1_,
#else
							data_->rawDataList[rawBufferIdx].devPtr,
							rawDataN2_,
#endif
							data_->gridValueRe.devPtr,
							data_->gridValueIm.devPtr,
							data_->rxApod.devPtr,
							cfConstants[2] /* factor */);
			checkKernelLaunchError();
		} else {
			processImageKernel<MFloat><<<procImageKernelGlobalSize / BLOCK_SIZE, BLOCK_SIZE>>>(
#ifdef USE_TRANSPOSE
							data_->rawDataT.devPtr,
							rawDataN1_,
#else
							data_->rawDataList[rawBufferIdx].devPtr,
							rawDataN2_,
#endif
							data_->gridValueRe.devPtr,
							data_->gridValueIm.devPtr,
							data_->rxApod.devPtr);
			checkKernelLaunchError();
		}
	}

	exec(data_->gridValueRe.copyDeviceToHost(numGridPoints * sizeof(MFloat)));
	exec(data_->gridValueIm.copyDeviceToHost(numGridPoints * sizeof(MFloat)));

	//==================================================
	// Read the formed image.
	//==================================================
	for (unsigned int col = 0; col < numCols; ++col) {
		for (unsigned int row = 0; row < minRowIdx_[col]; ++row) {
			gridValue(col, row) = 0;
		}
		unsigned int gridPointIdx = firstGridPointIdx_[col];
		for (unsigned int row = minRowIdx_[col]; row < numRows; ++row, ++gridPointIdx) {
			gridValue(col, row) = std::complex<MFloat>(
							data_->gridValueRe.hostPtr[gridPointIdx],
							data_->gridValueIm.hostPtr[gridPointIdx]);
		}
	}

#ifdef LAB_ENABLE_EXECUTION_TIME_MEASUREMENT
	tProcessColumnML.put(processColumnTimer.getTime());
#endif
	//LOG_DEBUG << "END ========== VectorialCombinedTwoMediumImagingCUDAProcessor::process ==========";
}

} // namespace Lab
