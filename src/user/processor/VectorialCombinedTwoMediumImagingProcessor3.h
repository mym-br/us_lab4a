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

#ifndef VECTORIALCOMBINEDTWOMEDIUMIMAGINGPROCESSOR3_H
#define VECTORIALCOMBINEDTWOMEDIUMIMAGINGPROCESSOR3_H

#include <algorithm> /* copy, fill */
#include <cmath> /* ceil, sqrt */
#include <cstddef> /* std::size_t */
#include <memory>
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
#include "Tensor3.h"
#include "Timer.h"
#include "TwoMediumSTAConfiguration.h"
#include "Util.h"
#include "XZ.h"



// Depends on the signal.
// 1.0 --> pi radian / sample at the original sampling rate.
#define VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_3_UPSAMP_FILTER_HALF_TRANSITION_WIDTH 0.2

#define VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_3_USE_ANALYTIC_SIGNAL 1
#define VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_3_USE_SIMD 1

#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_3_USE_ANALYTIC_SIGNAL
# include <complex>
#endif
#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_3_USE_SIMD
# include "SIMD.h"
#endif



namespace Lab {

// Faster version of VectorialCombinedTwoMediumImagingProcessor2.
// The grid must be rectangular.
template<typename TFloat>
class VectorialCombinedTwoMediumImagingProcessor3 {
public:
	struct StepConfiguration {
		unsigned int baseElemIdx;
		unsigned int baseElem;
		unsigned int firstTxElem;
		unsigned int lastTxElem;
	};

	VectorialCombinedTwoMediumImagingProcessor3(
			const TwoMediumSTAConfiguration<TFloat>& config,
			std::vector<Tensor3<TFloat>>& acqDataList,
			unsigned int upsamplingFactor,
#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_3_USE_ANALYTIC_SIGNAL
			AnalyticSignalCoherenceFactorProcessor<TFloat>& coherenceFactor,
#else
			CoherenceFactorProcessor<TFloat>& coherenceFactor,
#endif
			TFloat maxFermatBlockSize,
			TFloat peakOffset,
			unsigned int ascanStartOffset);
	~VectorialCombinedTwoMediumImagingProcessor3() = default;

	// In stepConfigList, only two cases are allowed.
	// One is with all the items using only one transmit element,
	// and the other is with all the items using three or more transmit elements.
	void process(
		const std::vector<StepConfiguration>& stepConfigList,
		const std::vector<XZ<TFloat>>& interfacePointList,
		const std::vector<TFloat>& rxApod,
		const Matrix<XZ<TFloat>>& gridXZ,
#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_3_USE_ANALYTIC_SIGNAL
		Matrix<std::complex<TFloat>>& gridValue
#else
		Matrix<TFloat>& gridValue
#endif
		);

#ifdef EXECUTION_TIME_MEASUREMENT_ACTIVE
	MeasurementList<double> tMinRowIdx;
	MeasurementList<double> tMedium1DelayMatrix;
	MeasurementList<double> tCalculateDelays;
	MeasurementList<double> tPrepareData;
	MeasurementList<double> tProcessColumn;
#endif

private:
	struct CalculateDelays;

	struct PrepareDataWithOneTxElem;
	struct PrepareData;
	struct PrepareDataThreadData {
		Interpolator<TFloat> interpolator;
		HilbertEnvelope<TFloat> envelope;
		std::vector<TFloat, tbb::cache_aligned_allocator<TFloat>> signal;
	};

	struct ProcessColumnWithOneTxElem;
	struct ProcessColumn;
	struct ProcessColumnThreadData {
#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_3_USE_ANALYTIC_SIGNAL
		AnalyticSignalCoherenceFactorProcessor<TFloat> coherenceFactor;
		std::vector<std::complex<TFloat>, tbb::cache_aligned_allocator<std::complex<TFloat>>> rxSignalSumList;
#else
		CoherenceFactorProcessor<TFloat> coherenceFactor;
		std::vector<TFloat, tbb::cache_aligned_allocator<TFloat>> rxSignalSumList;
#endif
	};

	VectorialCombinedTwoMediumImagingProcessor3(const VectorialCombinedTwoMediumImagingProcessor3&) = delete;
	VectorialCombinedTwoMediumImagingProcessor3& operator=(const VectorialCombinedTwoMediumImagingProcessor3&) = delete;
	VectorialCombinedTwoMediumImagingProcessor3(VectorialCombinedTwoMediumImagingProcessor3&&) = delete;
	VectorialCombinedTwoMediumImagingProcessor3& operator=(VectorialCombinedTwoMediumImagingProcessor3&&) = delete;

	const TwoMediumSTAConfiguration<TFloat>& config_;
	std::vector<Tensor3<TFloat>>& acqDataList_;
	unsigned int upsamplingFactor_;
#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_3_USE_ANALYTIC_SIGNAL
	AnalyticSignalCoherenceFactorProcessor<TFloat>& coherenceFactor_;
#else
	CoherenceFactorProcessor<TFloat>& coherenceFactor_;
#endif
	TFloat maxFermatBlockSize_;
	const TFloat lambda2_;
	std::size_t signalLength_;
#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_3_USE_ANALYTIC_SIGNAL
	Tensor3<std::complex<TFloat>, tbb::cache_aligned_allocator<std::complex<TFloat>>> signalMatrix_;
#else
	Tensor3<TFloat, tbb::cache_aligned_allocator<TFloat>> signalMatrix_;
#endif
	TFloat signalOffset_;
	std::vector<unsigned int, tbb::cache_aligned_allocator<unsigned int>> minRowIdx_; // for each column
	std::vector<TFloat, tbb::cache_aligned_allocator<TFloat>> xArray_;
	Matrix<TFloat, tbb::cache_aligned_allocator<TFloat>> medium1DelayMatrix_; // (interface_idx, element)
	Tensor3<TFloat, tbb::cache_aligned_allocator<TFloat>> delayMatrix_;
	std::unique_ptr<tbb::enumerable_thread_specific<PrepareDataThreadData>> prepareDataTLS_;
	std::unique_ptr<tbb::enumerable_thread_specific<ProcessColumnThreadData>> processColumnTLS_;
};



template<typename TFloat>
VectorialCombinedTwoMediumImagingProcessor3<TFloat>::VectorialCombinedTwoMediumImagingProcessor3(
			const TwoMediumSTAConfiguration<TFloat>& config,
			std::vector<Tensor3<TFloat>>& acqDataList,
			unsigned int upsamplingFactor,
#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_3_USE_ANALYTIC_SIGNAL
			AnalyticSignalCoherenceFactorProcessor<TFloat>& coherenceFactor,
#else
			CoherenceFactorProcessor<TFloat>& coherenceFactor,
#endif
			TFloat maxFermatBlockSize,
			TFloat peakOffset,
			unsigned int ascanStartOffset)
		: config_(config)
		, acqDataList_(acqDataList)
		, upsamplingFactor_(upsamplingFactor)
		, coherenceFactor_(coherenceFactor)
		, maxFermatBlockSize_(maxFermatBlockSize)
		, lambda2_(config_.propagationSpeed2 / config_.centerFrequency)
{
	const std::size_t ascanLength = acqDataList_[0].n3();

	signalOffset_ = (config_.samplingFrequency * upsamplingFactor_) * (peakOffset / config_.centerFrequency) - ascanStartOffset * upsamplingFactor_;
	signalLength_ = ascanLength * upsamplingFactor_;
	LOG_DEBUG << "signalOffset_: " << signalOffset_ << " signalLength_: " << signalLength_;

	PrepareDataThreadData prepareDataThreadData;
	if (upsamplingFactor_ > 1) {
		prepareDataThreadData.interpolator.prepare(upsamplingFactor_, VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_3_UPSAMP_FILTER_HALF_TRANSITION_WIDTH);
	}
	prepareDataThreadData.signal.resize(signalLength_);
	prepareDataTLS_ = std::make_unique<tbb::enumerable_thread_specific<PrepareDataThreadData>>(prepareDataThreadData);

	ProcessColumnThreadData processColumnThreadData;
	processColumnThreadData.coherenceFactor = coherenceFactor_;
	processColumnTLS_ = std::make_unique<tbb::enumerable_thread_specific<ProcessColumnThreadData>>(processColumnThreadData);
}

template<typename TFloat>
void
VectorialCombinedTwoMediumImagingProcessor3<TFloat>::process(
							const std::vector<StepConfiguration>& stepConfigList,
							const std::vector<XZ<TFloat>>& interfacePointList,
							const std::vector<TFloat>& rxApod,
							const Matrix<XZ<TFloat>>& gridXZ,
#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_3_USE_ANALYTIC_SIGNAL
							Matrix<std::complex<TFloat>>& gridValue
#else
							Matrix<TFloat>& gridValue
#endif
							)
{
	//LOG_DEBUG << "BEGIN ========== VectorialCombinedTwoMediumImagingProcessor3::process ==========";

	if (stepConfigList.empty()) {
		THROW_EXCEPTION(InvalidParameterException, "The list of step configurations is empty.");
	}
	if (gridXZ.n1() != gridValue.n1() || gridXZ.n2() != gridValue.n2()) {
		THROW_EXCEPTION(InvalidParameterException, "gridXZ and gridValue have different sizes.");
	}

	const std::size_t samplesPerChannelLow = acqDataList_[0].n3();

	// Clear the values in gridValue.
	for (auto& pointValue : gridValue) {
		pointValue = 0;
	}

	XZ<TFloat> p1 = interfacePointList[0];
	XZ<TFloat> p2 = interfacePointList[1];
	const TFloat dx = p2.x - p1.x;
	const TFloat dz = p2.z - p1.z;
	const TFloat r = std::sqrt(dx * dx + dz * dz);
	const unsigned int fermatBlockSize = FermatPrinciple::calcBlockSizeForTwoStepSearch(interfacePointList.size(), r, lambda2_, maxFermatBlockSize_);
	LOG_DEBUG << "fermatBlockSize: " << fermatBlockSize;

	minRowIdx_.resize(gridXZ.n1() /* number of columns */);
	delayMatrix_.resize(gridXZ.n1() /* number of columns */, gridXZ.n2() /* number of rows */, config_.numElementsMux);

#ifdef EXECUTION_TIME_MEASUREMENT_ACTIVE
	Timer minRowIdxTimer;
#endif
	//unsigned int numValidPixels = 0;
	const TFloat zStepGrid = gridXZ(0, 1).z - gridXZ(0, 0).z;
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
		//numValidPixels += gridXZ.n2() - minRowIdx_[col];
	}
	//LOG_DEBUG << "numValidPixels: " << numValidPixels;
#ifdef EXECUTION_TIME_MEASUREMENT_ACTIVE
	tMinRowIdx.put(minRowIdxTimer.getTime());
#endif
	const TFloat c2ByC1 = config_.propagationSpeed2 / config_.propagationSpeed1;

	ArrayGeometry::getElementX2D(config_.numElementsMux, config_.pitch, TFloat(0) /* offset */, xArray_);

#ifdef EXECUTION_TIME_MEASUREMENT_ACTIVE
	Timer medium1DelayMatrixTimer;
#endif
	medium1DelayMatrix_.resize(config_.numElementsMux, interfacePointList.size());
	for (unsigned int elem = 0; elem < config_.numElementsMux; ++elem) {
		TFloat* delays = &medium1DelayMatrix_(elem, 0);
		for (unsigned int i = 0; i < interfacePointList.size(); ++i) {
			const XZ<TFloat>& ifPoint = interfacePointList[i];
#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_3_USE_SIMD
			delays[i] = SIMD::calcDistance(xArray_[elem], 0, ifPoint.x, ifPoint.z) * c2ByC1;
#else
			const TFloat dx = ifPoint.x - xArray_[elem];
			const TFloat dz = ifPoint.z;
			delays[i] = std::sqrt(dx * dx + dz * dz) * c2ByC1;
#endif
		}
	}
#ifdef EXECUTION_TIME_MEASUREMENT_ACTIVE
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
#ifdef EXECUTION_TIME_MEASUREMENT_ACTIVE
	tCalculateDelays.put(calculateDelaysTimer.getTime());
#endif
	const auto& firstStepConfig = stepConfigList.front();
	if (firstStepConfig.firstTxElem == firstStepConfig.lastTxElem) {
		// Only one transmit element.
#ifdef EXECUTION_TIME_MEASUREMENT_ACTIVE
		Timer t0;
#endif
		signalMatrix_.resize(
				stepConfigList.size(),
				config_.numElements,
				signalLength_);

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
#ifdef EXECUTION_TIME_MEASUREMENT_ACTIVE
		tPrepareData.put(t0.getTime());

		Timer t1;
#endif
		ProcessColumnWithOneTxElem processColumnOp = {
			gridValue.n2(),
			config_,
			signalOffset_,
			signalMatrix_,
			rxApod,
			stepConfigList,
			minRowIdx_,
			delayMatrix_,
			*processColumnTLS_,
			gridValue
		};

		//tbb::parallel_for(tbb::blocked_range<unsigned int>(0, gridValue.n1()), processColumnOp);
		tbb::parallel_for(tbb::blocked_range<unsigned int>(0, gridValue.n1(), 1 /* grain size */), processColumnOp, tbb::simple_partitioner());
		//processColumnOp(tbb::blocked_range<unsigned int>(0, gridValue.n1())); // single-thread
#ifdef EXECUTION_TIME_MEASUREMENT_ACTIVE
		tProcessColumn.put(t1.getTime());
#endif
	} else {
		// Three or more transmit elements.

		for (const auto& stepConfig : stepConfigList) {
			Timer t0;

			signalMatrix_.resize(
					stepConfig.lastTxElem - stepConfig.firstTxElem + 1,
					config_.numElements,
					signalLength_);

			// Prepare the signal matrix.
			for (unsigned int txElem = stepConfig.firstTxElem; txElem <= stepConfig.lastTxElem; ++txElem) {
				PrepareData prepareDataOp = {
					samplesPerChannelLow,
					acqDataList_,
					upsamplingFactor_,
					stepConfig.baseElemIdx,
					txElem - stepConfig.firstTxElem,
					*prepareDataTLS_,
					signalMatrix_
				};
				tbb::parallel_for(tbb::blocked_range<unsigned int>(0, config_.numElements), prepareDataOp);
				//prepareDataOp(tbb::blocked_range<unsigned int>(0, config_.numElements)); // single-thread
			}
			LOG_DEBUG << "STEP PREP time = " << t0.getTime();

			Timer t1;
			ProcessColumn processColumnOp = {
				gridValue.n2(),
				config_,
				signalOffset_,
				signalMatrix_,
				rxApod,
				stepConfig,
				minRowIdx_,
				delayMatrix_,
				*processColumnTLS_,
				gridValue
			};
			tbb::parallel_for(tbb::blocked_range<unsigned int>(0, gridValue.n1()), processColumnOp);
			//processColumnOp(tbb::blocked_range<unsigned int>(0, gridValue.n1())); // single-thread
			LOG_DEBUG << "STEP PROC time = " << t1.getTime();
		}
	}

	//LOG_DEBUG << "END ========== VectorialCombinedTwoMediumImagingProcessor3::process ==========";
}



template<typename TFloat>
struct VectorialCombinedTwoMediumImagingProcessor3<TFloat>::CalculateDelays {
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
#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_3_USE_SIMD
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
#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_3_USE_SIMD
						tC2Min = medium1Delays[idxMin] + SIMD::calcDistance(ifPoint.x, ifPoint.z, point.x, point.z);
#else
						const TFloat dx2 = point.x - ifPoint.x;
						const TFloat dz2 = point.z - ifPoint.z;
						tC2Min = medium1Delays[idxMin] + std::sqrt(dx2 * dx2 + dz2 * dz2);
#endif
					}
					for (unsigned int idxSearch = idxMin + 1, end = interfacePointList.size(); idxSearch < end; ++idxSearch) {
						const XZ<TFloat>& ifPoint = interfacePointList[idxSearch];
#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_3_USE_SIMD
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
#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_3_USE_SIMD
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
struct VectorialCombinedTwoMediumImagingProcessor3<TFloat>::PrepareDataWithOneTxElem {
	void operator()(const tbb::blocked_range<unsigned int>& r) const {
		PrepareDataThreadData& local = prepareDataTLS.local();

		for (unsigned int rxElem = r.begin(); rxElem != r.end(); ++rxElem) {
			if (upsamplingFactor > 1) {
				// Interpolate the signal.
				local.interpolator.interpolate(&acqDataList[baseElementIdx](0 /* txElemIdx */, rxElem, 0), samplesPerChannelLow, &local.signal[0]);
			} else {
				auto range = acqDataList[baseElementIdx].range3(0 /* txElemIdx */, rxElem);
				std::copy(range.begin(), range.end(), local.signal.begin());
			}

			Util::removeDC(&local.signal[0], local.signal.size());

#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_3_USE_ANALYTIC_SIGNAL
			// Obtain the analytic signal.
			local.envelope.getAnalyticSignal(
					&local.signal[0],
					local.signal.size(),
					&signalMatrix(stepIdx, rxElem, 0));
#else
			std::copy(local.signal.begin(), local.signal.end(), &signalMatrix(stepIdx, rxElem, 0));
#endif
		}
	}

	// Use only small types or references/pointers, because this object may be copied many times.
	const std::size_t samplesPerChannelLow;
	const std::vector<Tensor3<TFloat>>& acqDataList;
	const unsigned int upsamplingFactor;
	const unsigned int stepIdx;
	const unsigned int baseElementIdx;
	tbb::enumerable_thread_specific<PrepareDataThreadData>& prepareDataTLS;
#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_3_USE_ANALYTIC_SIGNAL
	Tensor3<std::complex<TFloat>, tbb::cache_aligned_allocator<std::complex<TFloat>>>& signalMatrix;
#else
	Tensor3<TFloat, tbb::cache_aligned_allocator<TFloat>>& signalMatrix;
#endif
};



template<typename TFloat>
struct VectorialCombinedTwoMediumImagingProcessor3<TFloat>::PrepareData {
	void operator()(const tbb::blocked_range<unsigned int>& r) const {
		PrepareDataThreadData& local = prepareDataTLS.local();

		for (unsigned int rxElem = r.begin(); rxElem != r.end(); ++rxElem) {
			if (upsamplingFactor > 1) {
				// Interpolate the signal.
				local.interpolator.interpolate(&acqDataList[baseElementIdx](txElemIdx, rxElem, 0), samplesPerChannelLow, &local.signal[0]);
			} else {
				auto range = acqDataList[baseElementIdx].range3(txElemIdx, rxElem);
				std::copy(range.begin(), range.end(), local.signal.begin());
			}

			Util::removeDC(&local.signal[0], local.signal.size());

#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_3_USE_ANALYTIC_SIGNAL
			// Obtain the analytic signal.
			local.envelope.getAnalyticSignal(
					&local.signal[0],
					local.signal.size(),
					&signalMatrix(txElemIdx, rxElem, 0));
#else
			std::copy(local.signal.begin(), local.signal.end(), &signalMatrix(txElemIdx, rxElem, 0));
#endif
		}
	}

	// Use only small types or references/pointers, because this object may be copied many times.
	const std::size_t samplesPerChannelLow;
	const std::vector<Tensor3<TFloat>>& acqDataList;
	const unsigned int upsamplingFactor;
	const unsigned int baseElementIdx;
	const unsigned int txElemIdx;
	tbb::enumerable_thread_specific<PrepareDataThreadData>& prepareDataTLS;
#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_3_USE_ANALYTIC_SIGNAL
	Tensor3<std::complex<TFloat>, tbb::cache_aligned_allocator<std::complex<TFloat>>>& signalMatrix;
#else
	Tensor3<TFloat, tbb::cache_aligned_allocator<TFloat>>& signalMatrix;
#endif
};



template<typename TFloat>
struct VectorialCombinedTwoMediumImagingProcessor3<TFloat>::ProcessColumnWithOneTxElem {
	void operator()(const tbb::blocked_range<unsigned int>& r) const {
		//LOG_DEBUG << "col = " << r.begin() << " n = " << (r.end() - r.begin());
		ProcessColumnThreadData& local = processColumnTLS.local();

		local.rxSignalSumList.resize(config.numElements);
		const unsigned int signalLength = signalMatrix.n3();
		const unsigned int maxPosition = signalLength - 2;

		for (unsigned int col = r.begin(); col != r.end(); ++col) {
			for (unsigned int row = minRowIdx[col]; row < numRows; ++row) {
#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_3_USE_ANALYTIC_SIGNAL
				std::complex<TFloat> pointSum = 0.0;
#else
				TFloat pointSum = 0.0;
#endif
				for (const auto& stepConfig : stepConfigList) {
					const TFloat* delays = &delayMatrix(col, row, stepConfig.baseElem);
					const TFloat txDelay = delays[stepConfig.firstTxElem];
					const TFloat txOffset = signalOffset + txDelay;
					const auto* p = &signalMatrix(stepConfig.baseElemIdx, 0 /* rxElem */, 0);
					for (unsigned int rxElem = 0; rxElem < config.numElements; ++rxElem, p += signalLength) {
						// Linear interpolation.
						const TFloat position = txOffset + delays[rxElem];
						if (position >= 0.0f) {
							const unsigned int positionIdx = static_cast<unsigned int>(position);
							if (positionIdx <= maxPosition) {
								const TFloat k = position - positionIdx;
								const auto v0 = p[positionIdx];
								const auto v1 = p[positionIdx + 1];
								local.rxSignalSumList[rxElem] = v0 + k * (v1 - v0);
							} else {
								local.rxSignalSumList[rxElem] = 0;
							}
						} else {
							local.rxSignalSumList[rxElem] = 0;
						}
					}

#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_3_USE_ANALYTIC_SIGNAL
					std::complex<TFloat> sum = 0.0;
#else
					TFloat sum = 0.0;
#endif
					for (unsigned int rxElem = 0; rxElem < config.numElements; ++rxElem) {
						sum += local.rxSignalSumList[rxElem] * rxApod[rxElem];
					}

					if (local.coherenceFactor.enabled()) {
						pointSum += sum * local.coherenceFactor.calculate(&local.rxSignalSumList[0], local.rxSignalSumList.size());
					} else {
						pointSum += sum;
					}
				}

				gridValue(col, row) = pointSum;
			}
		}
	}

	// Use only small types or references/pointers, because this object may be copied many times.
	const std::size_t numRows;
	const TwoMediumSTAConfiguration<TFloat>& config;
	const TFloat signalOffset;
#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_3_USE_ANALYTIC_SIGNAL
	const Tensor3<std::complex<TFloat>, tbb::cache_aligned_allocator<std::complex<TFloat>>>& signalMatrix;
#else
	const Tensor3<TFloat, tbb::cache_aligned_allocator<TFloat>>& signalMatrix;
#endif
	const std::vector<TFloat>& rxApod;
	const std::vector<StepConfiguration>& stepConfigList;
	const std::vector<unsigned int, tbb::cache_aligned_allocator<unsigned int>>& minRowIdx;
	const Tensor3<TFloat, tbb::cache_aligned_allocator<TFloat>>& delayMatrix;
	tbb::enumerable_thread_specific<ProcessColumnThreadData>& processColumnTLS;
#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_3_USE_ANALYTIC_SIGNAL
	Matrix<std::complex<TFloat>>& gridValue;
#else
	Matrix<TFloat>& gridValue;
#endif
};



template<typename TFloat>
struct VectorialCombinedTwoMediumImagingProcessor3<TFloat>::ProcessColumn {
	void operator()(const tbb::blocked_range<unsigned int>& r) const {
		//LOG_DEBUG << "col = " << r.begin() << " n = " << (r.end() - r.begin());
		ProcessColumnThreadData& local = processColumnTLS.local();

		local.rxSignalSumList.resize(config.numElements);
		const unsigned int signalLength = signalMatrix.n3();
		const unsigned int maxPosition = signalLength - 2;

		for (unsigned int col = r.begin(); col != r.end(); ++col) {
			for (unsigned int row = minRowIdx[col]; row < numRows; ++row) {
#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_3_USE_ANALYTIC_SIGNAL
				std::fill(local.rxSignalSumList.begin(), local.rxSignalSumList.end(), std::complex<TFloat>(0));
#else
				std::fill(local.rxSignalSumList.begin(), local.rxSignalSumList.end(), TFloat(0));
#endif
				const TFloat* delays = &delayMatrix(col, row, stepConfig.baseElem);

				for (unsigned int txElem = stepConfig.firstTxElem; txElem <= stepConfig.lastTxElem; ++txElem) {
					const TFloat txDelay = delays[txElem];
					const TFloat txOffset = signalOffset + txDelay;
					const unsigned int txElemIdx = txElem - stepConfig.firstTxElem;
					const auto* p = &signalMatrix(txElemIdx, 0 /* rxElem */, 0);
					for (unsigned int rxElem = 0; rxElem < config.numElements; ++rxElem, p += signalLength) {
						// Linear interpolation.
						const TFloat position = txOffset + delays[rxElem];
						if (position >= 0.0f) {
							const unsigned int positionIdx = static_cast<unsigned int>(position);
							if (positionIdx <= maxPosition) {
								const TFloat k = position - positionIdx;
								const auto v0 = p[positionIdx];
								const auto v1 = p[positionIdx + 1];
								local.rxSignalSumList[rxElem] += v0 + k * (v1 - v0);
							}
						}
					}
				}

#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_3_USE_ANALYTIC_SIGNAL
				std::complex<TFloat> sum = 0.0;
#else
				TFloat sum = 0.0;
#endif
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
	const TFloat signalOffset;
#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_3_USE_ANALYTIC_SIGNAL
	const Tensor3<std::complex<TFloat>, tbb::cache_aligned_allocator<std::complex<TFloat>>>& signalMatrix;
#else
	const Tensor3<TFloat, tbb::cache_aligned_allocator<TFloat>>& signalMatrix;
#endif
	const std::vector<TFloat>& rxApod;
	const StepConfiguration& stepConfig;
	const std::vector<unsigned int, tbb::cache_aligned_allocator<unsigned int>>& minRowIdx;
	const Tensor3<TFloat, tbb::cache_aligned_allocator<TFloat>>& delayMatrix;
	tbb::enumerable_thread_specific<ProcessColumnThreadData>& processColumnTLS;
#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_3_USE_ANALYTIC_SIGNAL
	Matrix<std::complex<TFloat>>& gridValue;
#else
	Matrix<TFloat>& gridValue;
#endif
};

} // namespace Lab

#endif // VECTORIALCOMBINEDTWOMEDIUMIMAGINGPROCESSOR3_H
