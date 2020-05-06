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

#include "CUDACoherenceFactor.cuh"
#include "CUDAGeometry.cuh"
#include "CUDAUtil.h"



// Depends on the signal.
// 1.0 --> pi radian / sample at the original sampling rate.
#define UPSAMP_FILTER_HALF_TRANSITION_WIDTH (0.2)

#ifndef MFloat
# define MFloat float
#endif



namespace Lab {

// NVIDIA sm_50 or newer:
//   - Shared memory has 32 banks of 32 bits.

#define NUM_RX_ELEM 32

template<typename TFloat>
__device__
TFloat
calcTwoMediumTravelTime(TFloat x1, TFloat z1, TFloat xi, TFloat zi, TFloat x2, TFloat z2, TFloat invC1, TFloat invC2)
{
	const TFloat dx1 = xi - x1;
	const TFloat dz1 = zi - z1;
	const TFloat dx2 = x2 - xi;
	const TFloat dz2 = z2 - zi;
	return sqrt(dx1 * dx1 + dz1 * dz1) * invC1 + sqrt(dx2 * dx2 + dz2 * dz2) * invC2;
}

template<typename TFloat>
__device__
void
findMinTimeInTwoSteps(
		unsigned int blockSize,
		TFloat c1, TFloat c2,
		const TFloat (*interfacePointList)[2],
		unsigned int interfacePointListSize,
		TFloat x1, TFloat z1, TFloat x2, TFloat z2,
		TFloat& tMin, unsigned int& idxMin)
{
	const TFloat invC1 = 1 / c1;
	const TFloat invC2 = 1 / c2;

	// First step: step = blockSize
	{
		const TFloat (*point)[2] = interfacePointList;
		const TFloat t = calcTwoMediumTravelTime<TFloat>(x1, z1, (*point)[0], (*point)[1], x2, z2, invC1, invC2);
		tMin = t;
		idxMin = 0;
	}
	for (unsigned int i = 1; i < interfacePointListSize; i += blockSize) {
		const TFloat (*point)[2] = interfacePointList + i;
		const TFloat t = calcTwoMediumTravelTime<TFloat>(x1, z1, (*point)[0], (*point)[1], x2, z2, invC1, invC2);
		if (t < tMin) {
			tMin = t;
			idxMin = i;
		}
	}

	// Second step: step = 1
	const unsigned int iBegin = (idxMin > (blockSize - 1U)) ? idxMin - (blockSize - 1U) : 0;
	const unsigned int iEnd = umin(idxMin + blockSize, interfacePointListSize);
	for (unsigned int i = iBegin; i < iEnd; ++i) {
		const TFloat (*point)[2] = interfacePointList + i;
		const TFloat t = calcTwoMediumTravelTime<TFloat>(x1, z1, (*point)[0], (*point)[1], x2, z2, invC1, invC2);
		if (t < tMin) {
			tMin = t;
			idxMin = i;
		}
	}
}

// This function is not efficient in GPUs, but it avoids the transfer of the large delayTensor.
template<typename TFloat>
__global__
void
calculateDelaysTwoMediumKernel(
		unsigned int numCols,
		unsigned int numRows,
		unsigned int numElementsMux,
		TFloat fs,
		TFloat fsInvC2,
		TFloat c1,
		TFloat c2,
		unsigned int fermatBlockSize,
		const TFloat (*interfacePointList)[2],
		unsigned int interfacePointListSize,
		const TFloat* xArray,
		const unsigned int* minRowIdx,
		const TFloat* medium1DelayMatrix,
		const TFloat (*gridXZ)[2],
		TFloat* delayTensor)
{
	const unsigned int col = blockIdx.x * blockDim.x + threadIdx.x;
	if (col >= numCols) return;

	const unsigned int elem = blockIdx.y * blockDim.y + threadIdx.y;
	if (elem >= numElementsMux) return;

	if (minRowIdx[col] >= numRows) return;

	unsigned int lastInterfaceIdx = 0;

	// The first row above the interface.
	{
		const TFloat (*point)[2] = gridXZ + col * numRows + minRowIdx[col];

		// Fermat's principle. Find the fastest path.
		TFloat tMin;
		unsigned int idxMin;
		findMinTimeInTwoSteps<TFloat>(
				fermatBlockSize,
				c1, c2,
				interfacePointList,
				interfacePointListSize,
				xArray[elem], 0, (*point)[0], (*point)[1],
				tMin, idxMin);
		delayTensor[((elem * numCols) + col) * numRows + minRowIdx[col]] = tMin * fs;
		lastInterfaceIdx = idxMin;
	}

	const TFloat* medium1Delays = medium1DelayMatrix + elem * interfacePointListSize;

	for (unsigned int row = minRowIdx[col] + 1U; row < numRows; ++row) {
		const TFloat (*point)[2] = gridXZ + col * numRows + row;
		unsigned int idxMin = lastInterfaceIdx;
		TFloat tC2Min;
		{
			const TFloat (*ifPoint)[2] = interfacePointList + idxMin;
			tC2Min = medium1Delays[idxMin] + distance2D<TFloat>((*ifPoint)[0], (*ifPoint)[1], (*point)[0], (*point)[1]);
		}
		for (unsigned int idxSearch = idxMin + 1U; idxSearch < interfacePointListSize; ++idxSearch) {
			const TFloat (*ifPoint)[2] = interfacePointList + idxSearch;
			const TFloat tC2 = medium1Delays[idxSearch] + distance2D<TFloat>((*ifPoint)[0], (*ifPoint)[1], (*point)[0], (*point)[1]);
			if (tC2 >= tC2Min) {
				break;
			} else {
				tC2Min = tC2;
				idxMin = idxSearch;
			}
		}
		if (idxMin == lastInterfaceIdx) { // if the previous search was not successful
			for (int idxSearch = static_cast<int>(idxMin) - 1; idxSearch >= 0; --idxSearch) { // if idxMin = 0, idxSearch will start with -1
				const TFloat (*ifPoint)[2] = interfacePointList + idxSearch;
				const TFloat tC2 = medium1Delays[idxSearch] + distance2D<TFloat>((*ifPoint)[0], (*ifPoint)[1], (*point)[0], (*point)[1]);
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

template<typename TFloat>
__global__
void
processRowColumnWithOneTxElemKernel(
		unsigned int numCols,
		unsigned int numRows,
		TFloat signalOffset,
		const TFloat (*signalTensor)[2],
		unsigned int signalTensorN2,
		unsigned int signalTensorN3,
		unsigned int baseElem,
		unsigned int baseElemIdx,
		unsigned int txElem,
		const unsigned int* minRowIdx,
		const TFloat* delayTensor,
		TFloat (*gridValue)[2],
		const TFloat* rxApod)
{
	const unsigned int signalLength = signalTensorN3;
	const unsigned int maxPosition = signalLength - 2;

	const unsigned int col = blockIdx.y * blockDim.y + threadIdx.y;
	if (col >= numCols) return;

	const unsigned int row = blockIdx.x * blockDim.x + threadIdx.x;
	if (row < minRowIdx[col] || row >= numRows) return;

	const TFloat txDelay = delayTensor[((baseElem + txElem) * numCols + col) * numRows + row];
	const TFloat txOffset = signalOffset + txDelay;

	TFloat rxSignalSumRe = 0;
	TFloat rxSignalSumIm = 0;
	for (unsigned int rxElem = 0; rxElem < NUM_RX_ELEM; ++rxElem) {
		const TFloat (*p)[2] = signalTensor + baseElemIdx * signalTensorN2 * signalTensorN3
					+ rxElem * signalLength;
		// Linear interpolation.
		const TFloat position = txOffset + delayTensor[((baseElem + rxElem) * numCols + col) * numRows + row];
		if (position >= 0.0f) {
			const unsigned int positionIdx = static_cast<unsigned int>(position);
			if (positionIdx <= maxPosition) {
				const TFloat k = position - positionIdx;
				TFloat v0[2]; // complex
				v0[0] = p[positionIdx][0];
				v0[1] = p[positionIdx][1];
				TFloat v1[2]; // complex
				v1[0] = p[positionIdx + 1][0];
				v1[1] = p[positionIdx + 1][1];
				TFloat v[2]; // complex
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

template<typename TFloat>
__global__
void
processRowColumnWithOneTxElemPCFKernel(
		unsigned int numCols,
		unsigned int numRows,
		TFloat signalOffset,
		const TFloat (*signalTensor)[2],
		unsigned int signalTensorN2,
		unsigned int signalTensorN3,
		unsigned int baseElem,
		unsigned int baseElemIdx,
		unsigned int txElem,
		const unsigned int* minRowIdx,
		const TFloat* delayTensor,
		TFloat (*gridValue)[2],
		const TFloat* rxApod,
		TFloat pcfFactor)
{
	TFloat rxSignalListRe[NUM_RX_ELEM];
	TFloat rxSignalListIm[NUM_RX_ELEM];

	const unsigned int signalLength = signalTensorN3;
	const unsigned int maxPosition = signalLength - 2;

	const unsigned int col = blockIdx.y * blockDim.y + threadIdx.y;
	if (col >= numCols) return;

	const unsigned int row = blockIdx.x * blockDim.x + threadIdx.x;
	if (row < minRowIdx[col] || row >= numRows) return;

	const TFloat txDelay = delayTensor[((baseElem + txElem) * numCols + col) * numRows + row];
	const TFloat txOffset = signalOffset + txDelay;

	for (unsigned int rxElem = 0; rxElem < NUM_RX_ELEM; ++rxElem) {
		const TFloat (*p)[2] = signalTensor + baseElemIdx * signalTensorN2 * signalTensorN3
					+ rxElem * signalLength;
		// Linear interpolation.
		const TFloat position = txOffset + delayTensor[((baseElem + rxElem) * numCols + col) * numRows + row];
		if (position >= 0.0f) {
			const unsigned int positionIdx = static_cast<unsigned int>(position);
			if (positionIdx <= maxPosition) {
				const TFloat k = position - positionIdx;
				TFloat v0[2]; // complex
				v0[0] = p[positionIdx][0];
				v0[1] = p[positionIdx][1];
				TFloat v1[2]; // complex
				v1[0] = p[positionIdx + 1][0];
				v1[1] = p[positionIdx + 1][1];
				TFloat v[2]; // complex
				v[0] = v0[0] + k * (v1[0] - v0[0]);
				v[1] = v0[1] + k * (v1[1] - v0[1]);

				rxSignalListRe[rxElem] = v[0];
				rxSignalListIm[rxElem] = v[1];
			}
		}
	}

	const TFloat pcf = calcPCF<TFloat, NUM_RX_ELEM>(rxSignalListRe, rxSignalListIm, pcfFactor);
	TFloat sumRe = 0;
	TFloat sumIm = 0;
	for (int i = 0; i < NUM_RX_ELEM; ++i) {
		sumRe += rxSignalListRe[i] * rxApod[i];
		sumIm += rxSignalListIm[i] * rxApod[i];
	}

	const unsigned int point = col * numRows + row;
	gridValue[point][0] += sumRe * pcf;
	gridValue[point][1] += sumIm * pcf;
}

//=============================================================================

struct VectorialCombinedTwoMediumImagingCUDAProcessor2::CUDAData {
	CUDAContext c; // must be the first member
	bool cudaDataInitialized;
	CUDADevMem<MFloat[2]> gridXZ;
	CUDADevMem<MFloat> rxApod;
	CUDADevMem<MFloat> xArray;
	CUDAHostDevMem<MFloat[2]> signalTensor;
	CUDADevMem<unsigned int> minRowIdx;
	CUDADevMem<MFloat[2]> interfacePointList;
	CUDADevMem<MFloat> medium1DelayMatrix;
	CUDADevMem<MFloat> delayTensor;
	CUDAHostDevMem<MFloat[2]> gridValue;

	CUDAData()
		: cudaDataInitialized()
	{}
};

template<typename TFloat>
struct VectorialCombinedTwoMediumImagingCUDAProcessor2::PrepareDataWithOneTxElem {
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
	tbb::enumerable_thread_specific<typename VectorialCombinedTwoMediumImagingCUDAProcessor2::PrepareDataThreadData<TFloat>>& prepareDataTLS;
	TFloat (*signalTensor)[2];
	const unsigned int signalTensorN2;
	const unsigned int signalTensorN3;
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
		, stepConfigListSize_()
		, interfacePointListSize_()
		, rxApodSize_()
		, gridXZN1_()
		, gridXZN2_()
{
	if (config_.numElements != NUM_RX_ELEM) {
		THROW_EXCEPTION(InvalidParameterException, "The number of elements (" << config_.numElements <<
				") is not equal to " << NUM_RX_ELEM << '.');
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

	const unsigned int numCols = gridXZ.n1();
	const unsigned int numRows = gridXZ.n2();
	if (numCols == 0) {
		THROW_EXCEPTION(InvalidValueException, "Zero columns in the grid.");
	}

	minRowIdx_.resize(numCols);
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
	tMinRowIdxML.put(minRowIdxTimer.getTime());
#endif

	if (xArray_.empty()) {
		ArrayGeometry::getElementX2D(config_.numElementsMux, config_.pitch, MFloat(0) /* offset */, xArray_);
	}

	if (!data_->cudaDataInitialized) {
		stepConfigListSize_     = stepConfigList.size();
		interfacePointListSize_ = interfacePointList.size();
		rxApodSize_             = rxApod.size();
		gridXZN1_               = gridXZ.n1();
		gridXZN2_               = gridXZ.n2();

		data_->gridXZ             = CUDADevMem<MFloat[2]>(gridXZN1_ * gridXZN2_);
		data_->rxApod             = CUDADevMem<MFloat>(rxApodSize_);
		data_->xArray             = CUDADevMem<MFloat>(xArray_.size());
		data_->signalTensor       = CUDAHostDevMem<MFloat[2]>(stepConfigListSize_ * config_.numElements * signalLength_);
		data_->minRowIdx          = CUDADevMem<unsigned int>(minRowIdx_.size());
		data_->interfacePointList = CUDADevMem<MFloat[2]>(interfacePointListSize_);
		data_->medium1DelayMatrix = CUDADevMem<MFloat>(medium1DelayMatrix_.n1() * medium1DelayMatrix_.n2());
		data_->delayTensor        = CUDADevMem<MFloat>(config_.numElementsMux * numCols * numRows);
		data_->gridValue          = CUDAHostDevMem<MFloat[2]>(numCols * numRows);

		data_->cudaDataInitialized = true;
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
	const MFloat c2ByC1 = config_.propagationSpeed2 / config_.propagationSpeed1;
	for (unsigned int elem = 0; elem < config_.numElementsMux; ++elem) {
		MFloat* delays = &medium1DelayMatrix_(elem, 0);
		for (unsigned int i = 0; i < interfacePointList.size(); ++i) {
			const XZ<MFloat>& ifPoint = interfacePointList[i];
			delays[i] = Geometry::distance2DY0(xArray_[elem], ifPoint.x, ifPoint.z) * c2ByC1;
		}
	}
#ifdef USE_EXECUTION_TIME_MEASUREMENT
	tMedium1DelayMatrixML.put(medium1DelayMatrixTimer.getTime());
#endif

	// Prepare buffers.
	Timer prepareBuffersTimer;
	exec(cudaMemcpy(data_->gridXZ.devPtr            ,              gridXZ.data(),
				data_->gridXZ.sizeInBytes            , cudaMemcpyHostToDevice));
	exec(cudaMemcpy(data_->rxApod.devPtr            ,              rxApod.data(),
				data_->rxApod.sizeInBytes            , cudaMemcpyHostToDevice));
	exec(cudaMemcpy(data_->xArray.devPtr            ,             xArray_.data(),
				data_->xArray.sizeInBytes            , cudaMemcpyHostToDevice));
	exec(cudaMemcpy(data_->minRowIdx.devPtr         ,          minRowIdx_.data(),
				data_->minRowIdx.sizeInBytes         , cudaMemcpyHostToDevice));
	exec(cudaMemcpy(data_->interfacePointList.devPtr,  interfacePointList.data(),
				data_->interfacePointList.sizeInBytes, cudaMemcpyHostToDevice));
	exec(cudaMemcpy(data_->medium1DelayMatrix.devPtr, medium1DelayMatrix_.data(),
				data_->medium1DelayMatrix.sizeInBytes, cudaMemcpyHostToDevice));
	exec(cudaMemset(data_->gridValue.devPtr, 0, data_->gridValue.sizeInBytes));
	LOG_DEBUG << "PREPARE BUFFERS " << prepareBuffersTimer.getTime();

	{
#ifdef USE_EXECUTION_TIME_MEASUREMENT
		Timer calculateDelaysTimer;
#endif
		// Adjusted for GTX-1660.
		// Note: Less than 32 threads!
		const std::size_t colBlockSize = 2;
		const std::size_t elemBlockSize = 8;
		const std::size_t colNumBlocks  = CUDAUtil::numberOfBlocks(numCols, colBlockSize);
		const std::size_t elemNumBlocks = CUDAUtil::numberOfBlocks(config_.numElementsMux, elemBlockSize);
		const dim3 gridDim(colNumBlocks, elemNumBlocks);
		const dim3 blockDim(colBlockSize, elemBlockSize);

		calculateDelaysTwoMediumKernel<<<gridDim, blockDim>>>(
			numCols,
			numRows,
			config_.numElementsMux,
			config_.samplingFrequency * upsamplingFactor_,
			config_.samplingFrequency * upsamplingFactor_ / config_.propagationSpeed2,
			config_.propagationSpeed1,
			config_.propagationSpeed2,
			fermatBlockSize,
			data_->interfacePointList.devPtr,
			interfacePointList.size(),
			data_->xArray.devPtr,
			data_->minRowIdx.devPtr,
			data_->medium1DelayMatrix.devPtr,
			data_->gridXZ.devPtr,
			data_->delayTensor.devPtr);
		checkKernelLaunchError();

#ifdef USE_EXECUTION_TIME_MEASUREMENT
		//exec(cudaDeviceSynchronize());
		tCalculateDelaysML.put(calculateDelaysTimer.getTime());
#endif
	}

#ifdef USE_EXECUTION_TIME_MEASUREMENT
	Timer prepareDataTimer;
#endif
	// Only one transmit element.
	unsigned int stepIdx = 0;
	for (const auto& stepConfig : stepConfigList) {
		PrepareDataWithOneTxElem<MFloat> prepareDataOp = {
			samplesPerChannelLow,
			acqDataList_,
			upsamplingFactor_,
			stepIdx,
			stepConfig.baseElemIdx,
			*prepareDataTLS_,
			data_->signalTensor.hostPtr,
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

	//exec(cudaDeviceSynchronize());
	//Timer signalTransferTimer;

	exec(data_->signalTensor.copyHostToDevice());

	//LOG_DEBUG << "SIGNAL TRANSFER " << signalTransferTimer.getTime();

	//==================================================
	// Step configuration loop.
	//==================================================
	for (unsigned int i = 0; i < stepConfigList.size(); ++i) {
		const auto& stepConfig = stepConfigList[i];
		LOG_DEBUG << "stepConfig.baseElemIdx: " << stepConfig.baseElemIdx;

		exec(cudaDeviceSynchronize());

		//==================================================
		// Delay and sum.
		//==================================================

		//Timer delaySumTimer;

		if (coherenceFactor_.enabled()) {
			std::vector<MFloat> cfConstants;
			coherenceFactor_.implementation().getConstants(cfConstants);

			// Adjusted for GTX-1660.
			const std::size_t rowBlockSize = 16;
			const std::size_t colBlockSize = 16;
			const std::size_t rowNumBlocks = CUDAUtil::numberOfBlocks(numRows, rowBlockSize);
			const std::size_t colNumBlocks = CUDAUtil::numberOfBlocks(numCols, colBlockSize);
			const dim3 gridDim(rowNumBlocks, colNumBlocks);
			const dim3 blockDim(rowBlockSize, colBlockSize);

			processRowColumnWithOneTxElemPCFKernel<<<gridDim, blockDim>>>(
					numCols,
					numRows,
					signalOffset_,
					data_->signalTensor.devPtr,
					config_.numElements,
					signalLength_,
					stepConfig.baseElem,
					stepConfig.baseElemIdx,
					stepConfig.txElem,
					data_->minRowIdx.devPtr,
					data_->delayTensor.devPtr,
					data_->gridValue.devPtr,
					data_->rxApod.devPtr,
					cfConstants[2] /* factor */);
			checkKernelLaunchError();
		} else {
			// Adjusted for GTX-1660.
			const std::size_t rowBlockSize = 64;
			const std::size_t colBlockSize = 2;
			const std::size_t rowNumBlocks = CUDAUtil::numberOfBlocks(numRows, rowBlockSize);
			const std::size_t colNumBlocks = CUDAUtil::numberOfBlocks(numCols, colBlockSize);
			const dim3 gridDim(rowNumBlocks, colNumBlocks);
			const dim3 blockDim(rowBlockSize, colBlockSize);

			processRowColumnWithOneTxElemKernel<<<gridDim, blockDim>>>(
					numCols,
					numRows,
					signalOffset_,
					data_->signalTensor.devPtr,
					config_.numElements,
					signalLength_,
					stepConfig.baseElem,
					stepConfig.baseElemIdx,
					stepConfig.txElem,
					data_->minRowIdx.devPtr,
					data_->delayTensor.devPtr,
					data_->gridValue.devPtr,
					data_->rxApod.devPtr);
			checkKernelLaunchError();
		}

		//exec(cudaDeviceSynchronize());
		//LOG_DEBUG << "DELAY-SUM " << delaySumTimer.getTime();
	}

	//exec(cudaDeviceSynchronize());
	//Timer gridValueTransferTimer;

	exec(data_->gridValue.copyDeviceToHost());

	//==================================================
	// Read the formed image.
	//==================================================
	for (unsigned int col = 0; col < numCols; ++col) {
		for (unsigned int row = 0; row < minRowIdx_[col]; ++row) {
			gridValue(col, row) = 0;
		}
		unsigned int gridPointIdx = col * numRows + minRowIdx_[col];
		for (unsigned int row = minRowIdx_[col]; row < numRows; ++row, ++gridPointIdx) {
			gridValue(col, row) = std::complex<MFloat>(
							data_->gridValue.hostPtr[gridPointIdx][0],
							data_->gridValue.hostPtr[gridPointIdx][1]);
		}
	}

	//LOG_DEBUG << "GRID VALUE TRANSFER " << gridValueTransferTimer.getTime();

#ifdef USE_EXECUTION_TIME_MEASUREMENT
	tProcessColumnML.put(processColumnTimer.getTime());
#endif
	//LOG_DEBUG << "END ========== VectorialCombinedTwoMediumImagingCUDAProcessor2::process ==========";
}

} // namespace Lab
