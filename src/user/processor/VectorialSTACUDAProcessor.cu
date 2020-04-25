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

#include "VectorialSTACUDAProcessor.h"

#include <algorithm> /* copy */

#include <tbb/partitioner.h>
#include <tbb/tbb.h>

#include "cuda.h"

#include "ArrayGeometry.h"
#include "Exception.h"
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

#define NUM_RX_ELEM 32

template<typename TFloat>
__global__
void
calculateDelaysSTAKernel(
		unsigned int numCols,
		unsigned int numRows,
		TFloat invCT,
		const TFloat* xArray,
		const TFloat (*gridXZ)[2],
		TFloat* delayTensor)
{
	const unsigned int row = blockIdx.x * blockDim.x + threadIdx.x;
	if (row >= numRows) return;

	const unsigned int col = blockIdx.y * blockDim.y + threadIdx.y;
	if (col >= numCols) return;

	const unsigned int elem = blockIdx.z * blockDim.z + threadIdx.z;
	if (elem >= NUM_RX_ELEM) return;

	const TFloat (*point)[2] = gridXZ + col * numRows + row;
	delayTensor[((elem * numCols) + col) * numRows + row] =
			distance2DY0<TFloat>(xArray[elem], (*point)[0], (*point)[1]) * invCT;
}

template<typename TFloat>
__global__
void
processRowColumnSTAKernel(
		unsigned int numCols,
		unsigned int numRows,
		TFloat signalOffset,
		const TFloat (*signalTensor)[2],
		unsigned int signalTensorN2,
		unsigned int signalTensorN3,
		unsigned int firstTxElem,
		unsigned int lastTxElem,
		const TFloat* rxApod,
		const TFloat* delayTensor,
		TFloat* gridValue)
{
	const unsigned int signalLength = signalTensorN3;
	const unsigned int maxPosition = signalLength - 2;

	const unsigned int row = blockIdx.x * blockDim.x + threadIdx.x;
	if (row >= numRows) return;

	const unsigned int col = blockIdx.y * blockDim.y + threadIdx.y;
	if (col >= numCols) return;

	TFloat sumRe = 0;
	TFloat sumIm = 0;
	for (int txElem = firstTxElem; txElem <= lastTxElem; ++txElem) {
		const TFloat txDelay = delayTensor[(txElem * numCols + col) * numRows + row];
		const TFloat txOffset = signalOffset + txDelay;
		for (int rxElem = 0; rxElem < NUM_RX_ELEM; ++rxElem) {
			const TFloat (*p)[2] = signalTensor + (((txElem - firstTxElem) * signalTensorN2) + rxElem) * signalLength;
			// Linear interpolation.
			const TFloat position = txOffset + delayTensor[(rxElem * numCols + col) * numRows + row];
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

					sumRe += v[0] * rxApod[rxElem];
					sumIm += v[1] * rxApod[rxElem];
				}
			}
		}
	}
	const unsigned int point = col * numRows + row;
	gridValue[point] = sqrt(sumRe * sumRe + sumIm * sumIm);
}

template<typename TFloat>
__global__
void
processRowColumnSTAPCFKernel(
		unsigned int numCols,
		unsigned int numRows,
		TFloat signalOffset,
		const TFloat (*signalTensor)[2],
		unsigned int signalTensorN2,
		unsigned int signalTensorN3,
		unsigned int firstTxElem,
		unsigned int lastTxElem,
		const TFloat* rxApod,
		const TFloat* delayTensor,
		TFloat pcfFactor,
		TFloat* gridValue,
		TFloat* gridFactor)
{
	TFloat rxSignalListRe[NUM_RX_ELEM];
	TFloat rxSignalListIm[NUM_RX_ELEM];

	const unsigned int signalLength = signalTensorN3;
	const unsigned int maxPosition = signalLength - 2;

	const unsigned int row = blockIdx.x * blockDim.x + threadIdx.x;
	if (row >= numRows) return;

	const unsigned int col = blockIdx.y * blockDim.y + threadIdx.y;
	if (col >= numCols) return;

	for (int i = 0; i < NUM_RX_ELEM; ++i) {
		rxSignalListRe[i] = 0;
		rxSignalListIm[i] = 0;
	}
	for (int txElem = firstTxElem; txElem <= lastTxElem; ++txElem) {
		const TFloat txDelay = delayTensor[(txElem * numCols + col) * numRows + row];
		const TFloat txOffset = signalOffset + txDelay;
		for (int rxElem = 0; rxElem < NUM_RX_ELEM; ++rxElem) {
			const TFloat (*p)[2] = signalTensor + (((txElem - firstTxElem) * signalTensorN2) + rxElem) * signalLength;
			// Linear interpolation.
			const TFloat position = txOffset + delayTensor[(rxElem * numCols + col) * numRows + row];
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

					rxSignalListRe[rxElem] += v[0] * rxApod[rxElem];
					rxSignalListIm[rxElem] += v[1] * rxApod[rxElem];
				}
			}
		}
	}

	const TFloat pcf = calcPCF<TFloat, NUM_RX_ELEM>(rxSignalListRe, rxSignalListIm, pcfFactor);
	TFloat sumRe = 0;
	TFloat sumIm = 0;
	for (int i = 0; i < NUM_RX_ELEM; ++i) {
		sumRe += rxSignalListRe[i];
		sumIm += rxSignalListIm[i];
	}

	const unsigned int point = col * numRows + row;
	gridValue[point] = sqrt(sumRe * sumRe + sumIm * sumIm);
	gridFactor[point] = pcf;
}

//=============================================================================

struct VectorialSTACUDAProcessor::CUDAData {
	bool cudaDataInitialized;
	CUDAHostDevMem<MFloat[2]> gridXZ;
	CUDADevMem<MFloat> rxApod;
	CUDADevMem<MFloat> xArray;
	CUDAHostDevMem<MFloat[2]> signalTensor;
	CUDADevMem<MFloat> delayTensor;
	CUDAHostDevMem<MFloat> gridValue;
	CUDAHostDevMem<MFloat> gridFactor;

	CUDAData()
		: cudaDataInitialized()
	{}
};

template<typename TFloat>
struct VectorialSTACUDAProcessor::PrepareData {
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
	tbb::enumerable_thread_specific<typename VectorialSTACUDAProcessor::PrepareDataThreadData>& prepareDataTLS;
	TFloat (*signalTensor)[2];
	const unsigned int signalTensorN2;
	const unsigned int signalTensorN3;
	const unsigned int txElemIdx;
	const unsigned int deadZoneSamplesUp;
};

VectorialSTACUDAProcessor::VectorialSTACUDAProcessor(
			const STAConfiguration<MFloat>& config,
			STAAcquisition<MFloat>& acquisition,
			unsigned int upsamplingFactor,
			AnalyticSignalCoherenceFactorProcessor<MFloat>& coherenceFactor,
			MFloat peakOffset,
			const std::vector<MFloat>& rxApod)
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
{
	if (config_.numElements != NUM_RX_ELEM) {
		THROW_EXCEPTION(InvalidParameterException, "The number of elements (" << config_.numElements <<
				") is not equal to " << NUM_RX_ELEM << '.');
	}

	signalOffset_ = (config_.samplingFrequency * upsamplingFactor_) * (peakOffset / config_.centerFrequency);
	LOG_DEBUG << "signalOffset_: " << signalOffset_;

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
VectorialSTACUDAProcessor::~VectorialSTACUDAProcessor()
{
}

void
VectorialSTACUDAProcessor::prepare(unsigned int baseElement)
{
	acquisition_.prepare(baseElement);
}

void
VectorialSTACUDAProcessor::process(Matrix<XYZValueFactor<MFloat>>& gridData)
{
	//LOG_DEBUG << "BEGIN ========== VectorialSTACUDAProcessor::process ==========";

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

#ifdef USE_EXECUTION_TIME_MEASUREMENT
		Timer acquisitionTimer;
#endif
		acquisition_.execute(txElem, acqData_);

#ifdef USE_EXECUTION_TIME_MEASUREMENT
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
				prepareDataThreadData.interpolator.prepare(upsamplingFactor_, UPSAMP_FILTER_HALF_TRANSITION_WIDTH);
			}
			prepareDataThreadData.signal.resize(signalLength_);
			prepareDataTLS_ = std::make_unique<tbb::enumerable_thread_specific<PrepareDataThreadData>>(prepareDataThreadData);

			initialized_ = true;
		}
		if (!data_->cudaDataInitialized) {
			gridN1_ = gridData.n1();
			gridN2_ = gridData.n2();

			data_->gridXZ       = CUDAHostDevMem<MFloat[2]>(gridN1_ * gridN2_);
			data_->rxApod       = CUDADevMem<MFloat>(rxApod_.size());
			data_->xArray       = CUDADevMem<MFloat>(xArray_.size());
			data_->signalTensor = CUDAHostDevMem<MFloat[2]>(numTxElem * config_.numElements * signalLength_);
			data_->delayTensor  = CUDADevMem<MFloat>(config_.numElements * numCols * numRows);
			data_->gridValue    = CUDAHostDevMem<MFloat>(numCols * numRows);
			if (coherenceFactor_.enabled()) {
				data_->gridFactor = CUDAHostDevMem<MFloat>(numCols * numRows);
			}

			data_->cudaDataInitialized = true;
		} else {
			if (gridN1_ != gridData.n1() ||
			    gridN2_ != gridData.n2()) {
				THROW_EXCEPTION(InvalidParameterException, "Data size mismatch.");
			}
		}

		const unsigned int samplesPerChannelLow = acqData_.n2();

#ifdef USE_EXECUTION_TIME_MEASUREMENT
		Timer prepareDataTimer;
#endif
		PrepareData<MFloat> prepareDataOp = {
			samplesPerChannelLow,
			acqData_,
			upsamplingFactor_,
			*prepareDataTLS_,
			data_->signalTensor.hostPtr,
			config_.numElements,
			signalLength_,
			txElem - config_.firstTxElem,
			deadZoneSamplesUp_
		};
		//tbb::parallel_for(tbb::blocked_range<unsigned int>(0, config_.numElements), prepareDataOp);
		tbb::parallel_for(tbb::blocked_range<unsigned int>(0, config_.numElements, 1 /* grain size */), prepareDataOp, tbb::simple_partitioner());
		//prepareDataOp(tbb::blocked_range<unsigned int>(0, config_.numElements)); // single-thread

#ifdef USE_EXECUTION_TIME_MEASUREMENT
		tPrepareDataML.put(prepareDataTimer.getTime());
#endif
	}

	if (!initialized_) THROW_EXCEPTION(InvalidStateException, "Not initialized.");

	// Prepare CUDA buffers.
	Timer prepareBuffersTimer;
	for (std::size_t i = 0; i < gridData.size(); ++i) {
		const XYZValueFactor<MFloat>& point = *(gridData.data() + i);
		data_->gridXZ.hostPtr[i][0] = point.x;
		data_->gridXZ.hostPtr[i][1] = point.z;
	}
	exec(data_->gridXZ.copyHostToDevice());
	exec(cudaMemcpy(data_->rxApod.devPtr, rxApod_.data(),
				data_->rxApod.sizeInBytes, cudaMemcpyHostToDevice));
	exec(cudaMemcpy(data_->xArray.devPtr, xArray_.data(),
				data_->xArray.sizeInBytes, cudaMemcpyHostToDevice));
	exec(data_->signalTensor.copyHostToDevice());
	exec(cudaMemset(data_->gridValue.devPtr, 0, data_->gridValue.sizeInBytes));
	LOG_DEBUG << "PREPARE BUFFERS " << prepareBuffersTimer.getTime();

	{
#ifdef USE_EXECUTION_TIME_MEASUREMENT
		Timer calculateDelaysTimer;
#endif
		const std::size_t rowBlockSize = 32;
		const std::size_t colBlockSize = 1;
		const std::size_t elemBlockSize = 1;
		const std::size_t rowNumBlocks  = CUDAUtil::numberOfBlocks(numRows, rowBlockSize);
		const std::size_t colNumBlocks  = CUDAUtil::numberOfBlocks(numCols, colBlockSize);
		const std::size_t elemNumBlocks = CUDAUtil::numberOfBlocks(NUM_RX_ELEM, elemBlockSize);
		const dim3 gridDim(rowNumBlocks, colNumBlocks, elemNumBlocks);
		const dim3 blockDim(rowBlockSize, colBlockSize, elemBlockSize);

		calculateDelaysSTAKernel<<<gridDim, blockDim>>>(
			numCols,
			numRows,
			config_.samplingFrequency * upsamplingFactor_ / config_.propagationSpeed,
			data_->xArray.devPtr,
			data_->gridXZ.devPtr,
			data_->delayTensor.devPtr);
		checkKernelLaunchError();

#ifdef USE_EXECUTION_TIME_MEASUREMENT
		exec(cudaDeviceSynchronize());
		tCalculateDelaysML.put(calculateDelaysTimer.getTime());
#endif
	}

	//==================================================
	// Delay and sum.
	//==================================================

#ifdef USE_EXECUTION_TIME_MEASUREMENT
	Timer delaySumTimer;
#endif
	if (coherenceFactor_.enabled()) {
		std::vector<MFloat> cfConstants;
		coherenceFactor_.implementation().getConstants(cfConstants);

		// Adjusted for GTX-1660.
		const std::size_t rowBlockSize = 4;
		const std::size_t colBlockSize = 16;
		const std::size_t rowNumBlocks = CUDAUtil::numberOfBlocks(numRows, rowBlockSize);
		const std::size_t colNumBlocks = CUDAUtil::numberOfBlocks(numCols, colBlockSize);
		const dim3 gridDim(rowNumBlocks, colNumBlocks);
		const dim3 blockDim(rowBlockSize, colBlockSize);

		processRowColumnSTAPCFKernel<<<gridDim, blockDim>>>(
				numCols,
				numRows,
				signalOffset_,
				data_->signalTensor.devPtr,
				config_.numElements,
				signalLength_,
				config_.firstTxElem,
				config_.lastTxElem,
				data_->rxApod.devPtr,
				data_->delayTensor.devPtr,
				cfConstants[2], /* factor */
				data_->gridValue.devPtr,
				data_->gridFactor.devPtr);
		checkKernelLaunchError();
	} else {
		// Adjusted for GTX-1660.
		const std::size_t rowBlockSize = 16;
		const std::size_t colBlockSize = 16;
		const std::size_t rowNumBlocks = CUDAUtil::numberOfBlocks(numRows, rowBlockSize);
		const std::size_t colNumBlocks = CUDAUtil::numberOfBlocks(numCols, colBlockSize);
		const dim3 gridDim(rowNumBlocks, colNumBlocks);
		const dim3 blockDim(rowBlockSize, colBlockSize);

		processRowColumnSTAKernel<<<gridDim, blockDim>>>(
				numCols,
				numRows,
				signalOffset_,
				data_->signalTensor.devPtr,
				config_.numElements,
				signalLength_,
				config_.firstTxElem,
				config_.lastTxElem,
				data_->rxApod.devPtr,
				data_->delayTensor.devPtr,
				data_->gridValue.devPtr);
		checkKernelLaunchError();
	}

#ifdef USE_EXECUTION_TIME_MEASUREMENT
	exec(cudaDeviceSynchronize());
	tDelaySumML.put(delaySumTimer.getTime());
#endif

	exec(data_->gridValue.copyDeviceToHost());
	if (coherenceFactor_.enabled()) {
		exec(data_->gridFactor.copyDeviceToHost());
	}

	const unsigned int numSignals = numTxElem * config_.numElements;
	const MFloat valueCoef = 1.0f / MFloat(numSignals);

	//==================================================
	// Read the formed image.
	//==================================================
	for (unsigned int col = 0; col < numCols; ++col) {
		unsigned int gridPointIdx = col * numRows;
		for (unsigned int row = 0; row < numRows; ++row, ++gridPointIdx) {
			gridData(col, row).value = data_->gridValue.hostPtr[gridPointIdx] * valueCoef;
		}
	}
	if (coherenceFactor_.enabled()) {
		for (unsigned int col = 0; col < numCols; ++col) {
			unsigned int gridPointIdx = col * numRows;
			for (unsigned int row = 0; row < numRows; ++row, ++gridPointIdx) {
				gridData(col, row).factor = data_->gridFactor.hostPtr[gridPointIdx];
			}
		}
	}

	//LOG_DEBUG << "END ========== VectorialSTACUDAProcessor::process ==========";
}

} // namespace Lab
