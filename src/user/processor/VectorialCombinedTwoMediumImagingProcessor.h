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

#ifndef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_H
#define VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_H

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
#define VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_UPSAMP_FILTER_HALF_TRANSITION_WIDTH 0.2

#define VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_USE_ANALYTIC_SIGNAL 1

#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_USE_ANALYTIC_SIGNAL
# include <complex>
#endif



namespace Lab {

// Two-medium image formation, using analytic signals (each sample is a real-imag vector).
// The final image is a combination of sub-images, using apodization.
//
// The grid must be rectangular.
//
template<typename TFloat>
class VectorialCombinedTwoMediumImagingProcessor {
public:
	struct StepConfiguration {
		unsigned int baseElemIdx;
		unsigned int baseElem;
		unsigned int firstTxElem;
		unsigned int lastTxElem;
	};

	VectorialCombinedTwoMediumImagingProcessor(
			const TwoMediumSTAConfiguration<TFloat>& config,
			std::vector<Tensor3<TFloat>>& acqDataList,
			unsigned int upsamplingFactor,
#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_USE_ANALYTIC_SIGNAL
			AnalyticSignalCoherenceFactorProcessor<TFloat>& coherenceFactor,
#else
			CoherenceFactorProcessor<TFloat>& coherenceFactor,
#endif
			TFloat maxFermatBlockSize,
			TFloat peakOffset,
			unsigned int signalStartOffset);
	~VectorialCombinedTwoMediumImagingProcessor() = default;

	// In stepConfigList, only two cases are allowed.
	// One is with all the items using only one transmit element,
	// and the other is with all the items using three or more transmit elements.
	void process(
		const std::vector<StepConfiguration>& stepConfigList,
		const std::vector<XZ<TFloat>>& interfacePointList,
		const std::vector<TFloat>& rxApod,
		const Matrix<XZ<TFloat>>& gridXZ,
#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_USE_ANALYTIC_SIGNAL
		Matrix<std::complex<TFloat>>& gridValue
#else
		Matrix<TFloat>& gridValue
#endif
		);

#ifdef LAB_ENABLE_EXECUTION_TIME_MEASUREMENT
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
#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_USE_ANALYTIC_SIGNAL
		AnalyticSignalCoherenceFactorProcessor<TFloat> coherenceFactor;
		std::vector<std::complex<TFloat>, tbb::cache_aligned_allocator<std::complex<TFloat>>> rxSignalSumList;
#else
		CoherenceFactorProcessor<TFloat> coherenceFactor;
		std::vector<TFloat, tbb::cache_aligned_allocator<TFloat>> rxSignalSumList;
#endif
	};

	VectorialCombinedTwoMediumImagingProcessor(const VectorialCombinedTwoMediumImagingProcessor&) = delete;
	VectorialCombinedTwoMediumImagingProcessor& operator=(const VectorialCombinedTwoMediumImagingProcessor&) = delete;
	VectorialCombinedTwoMediumImagingProcessor(VectorialCombinedTwoMediumImagingProcessor&&) = delete;
	VectorialCombinedTwoMediumImagingProcessor& operator=(VectorialCombinedTwoMediumImagingProcessor&&) = delete;

	const TwoMediumSTAConfiguration<TFloat>& config_;
	std::vector<Tensor3<TFloat>>& acqDataList_;
	unsigned int upsamplingFactor_;
#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_USE_ANALYTIC_SIGNAL
	AnalyticSignalCoherenceFactorProcessor<TFloat>& coherenceFactor_;
#else
	CoherenceFactorProcessor<TFloat>& coherenceFactor_;
#endif
	TFloat maxFermatBlockSize_;
	const TFloat lambda2_;
	std::size_t signalLength_;
#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_USE_ANALYTIC_SIGNAL
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
VectorialCombinedTwoMediumImagingProcessor<TFloat>::VectorialCombinedTwoMediumImagingProcessor(
			const TwoMediumSTAConfiguration<TFloat>& config,
			std::vector<Tensor3<TFloat>>& acqDataList,
			unsigned int upsamplingFactor,
#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_USE_ANALYTIC_SIGNAL
			AnalyticSignalCoherenceFactorProcessor<TFloat>& coherenceFactor,
#else
			CoherenceFactorProcessor<TFloat>& coherenceFactor,
#endif
			TFloat maxFermatBlockSize,
			TFloat peakOffset,
			unsigned int signalStartOffset)
		: config_(config)
		, acqDataList_(acqDataList)
		, upsamplingFactor_(upsamplingFactor)
		, coherenceFactor_(coherenceFactor)
		, maxFermatBlockSize_(maxFermatBlockSize)
		, lambda2_(config_.propagationSpeed2 / config_.centerFrequency)
{
	const std::size_t origSignalLength = acqDataList_[0].n3();

	signalOffset_ = (config_.samplingFrequency * upsamplingFactor_) * (peakOffset / config_.centerFrequency) - signalStartOffset * upsamplingFactor_;
	signalLength_ = origSignalLength * upsamplingFactor_;
	LOG_DEBUG << "signalOffset_: " << signalOffset_ << " signalLength_: " << signalLength_;

	PrepareDataThreadData prepareDataThreadData;
	if (upsamplingFactor_ > 1) {
		prepareDataThreadData.interpolator.prepare(upsamplingFactor_, VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_UPSAMP_FILTER_HALF_TRANSITION_WIDTH);
	}
	prepareDataThreadData.signal.resize(signalLength_);
	prepareDataTLS_ = std::make_unique<tbb::enumerable_thread_specific<PrepareDataThreadData>>(prepareDataThreadData);

	ProcessColumnThreadData processColumnThreadData;
	processColumnThreadData.coherenceFactor = coherenceFactor_;
	processColumnTLS_ = std::make_unique<tbb::enumerable_thread_specific<ProcessColumnThreadData>>(processColumnThreadData);
}

template<typename TFloat>
void
VectorialCombinedTwoMediumImagingProcessor<TFloat>::process(
							const std::vector<StepConfiguration>& stepConfigList,
							const std::vector<XZ<TFloat>>& interfacePointList,
							const std::vector<TFloat>& rxApod,
							const Matrix<XZ<TFloat>>& gridXZ,
#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_USE_ANALYTIC_SIGNAL
							Matrix<std::complex<TFloat>>& gridValue
#else
							Matrix<TFloat>& gridValue
#endif
							)
{
	//LOG_DEBUG << "BEGIN ========== VectorialCombinedTwoMediumImagingProcessor::process ==========";

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

#ifdef LAB_ENABLE_EXECUTION_TIME_MEASUREMENT
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
#ifdef LAB_ENABLE_EXECUTION_TIME_MEASUREMENT
	tMinRowIdxML.put(minRowIdxTimer.getTime());
#endif
	const TFloat c2ByC1 = config_.propagationSpeed2 / config_.propagationSpeed1;

	ArrayGeometry::getElementX2D(config_.numElementsMux, config_.pitch, TFloat(0) /* offset */, xArray_);

#ifdef LAB_ENABLE_EXECUTION_TIME_MEASUREMENT
	Timer medium1DelayMatrixTimer;
#endif
	medium1DelayMatrix_.resize(config_.numElementsMux, interfacePointList.size());
	for (unsigned int elem = 0; elem < config_.numElementsMux; ++elem) {
		TFloat* delays = &medium1DelayMatrix_(elem, 0);
		for (unsigned int i = 0; i < interfacePointList.size(); ++i) {
			const XZ<TFloat>& ifPoint = interfacePointList[i];
			delays[i] = Geometry::distance2DY0(xArray_[elem], ifPoint.x, ifPoint.z) * c2ByC1;
		}
	}
#ifdef LAB_ENABLE_EXECUTION_TIME_MEASUREMENT
	tMedium1DelayMatrixML.put(medium1DelayMatrixTimer.getTime());

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
#ifdef LAB_ENABLE_EXECUTION_TIME_MEASUREMENT
	tCalculateDelaysML.put(calculateDelaysTimer.getTime());
#endif
	const auto& firstStepConfig = stepConfigList.front();
	if (firstStepConfig.firstTxElem == firstStepConfig.lastTxElem) {
		// Only one transmit element.
#ifdef LAB_ENABLE_EXECUTION_TIME_MEASUREMENT
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
#ifdef LAB_ENABLE_EXECUTION_TIME_MEASUREMENT
		tPrepareDataML.put(t0.getTime());

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
#ifdef LAB_ENABLE_EXECUTION_TIME_MEASUREMENT
		tProcessColumnML.put(t1.getTime());
#endif
	} else {
#ifdef LAB_ENABLE_EXECUTION_TIME_MEASUREMENT
		double t0Sum = 0.0;
		double t1Sum = 0.0;
#endif
		// Three or more transmit elements.
		for (const auto& stepConfig : stepConfigList) {
#ifdef LAB_ENABLE_EXECUTION_TIME_MEASUREMENT
			Timer t0;
#endif
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
#ifdef LAB_ENABLE_EXECUTION_TIME_MEASUREMENT
			t0Sum += t0.getTime();

			Timer t1;
#endif
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
#ifdef LAB_ENABLE_EXECUTION_TIME_MEASUREMENT
			t1Sum += t1.getTime();
#endif
		}
#ifdef LAB_ENABLE_EXECUTION_TIME_MEASUREMENT
		tPrepareDataML.put(t0Sum);
		tProcessColumnML.put(t1Sum);
#endif
	}

	//LOG_DEBUG << "END ========== VectorialCombinedTwoMediumImagingProcessor::process ==========";
}



template<typename TFloat>
struct VectorialCombinedTwoMediumImagingProcessor<TFloat>::CalculateDelays {
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
struct VectorialCombinedTwoMediumImagingProcessor<TFloat>::PrepareDataWithOneTxElem {
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

#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_USE_ANALYTIC_SIGNAL
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
#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_USE_ANALYTIC_SIGNAL
	Tensor3<std::complex<TFloat>, tbb::cache_aligned_allocator<std::complex<TFloat>>>& signalMatrix;
#else
	Tensor3<TFloat, tbb::cache_aligned_allocator<TFloat>>& signalMatrix;
#endif
};



template<typename TFloat>
struct VectorialCombinedTwoMediumImagingProcessor<TFloat>::PrepareData {
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

#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_USE_ANALYTIC_SIGNAL
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
#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_USE_ANALYTIC_SIGNAL
	Tensor3<std::complex<TFloat>, tbb::cache_aligned_allocator<std::complex<TFloat>>>& signalMatrix;
#else
	Tensor3<TFloat, tbb::cache_aligned_allocator<TFloat>>& signalMatrix;
#endif
};



template<typename TFloat>
struct VectorialCombinedTwoMediumImagingProcessor<TFloat>::ProcessColumnWithOneTxElem {
	void operator()(const tbb::blocked_range<unsigned int>& r) const {
		//LOG_DEBUG << "col = " << r.begin() << " n = " << (r.end() - r.begin());
		ProcessColumnThreadData& local = processColumnTLS.local();

		local.rxSignalSumList.resize(config.numElements);
		const unsigned int signalLength = signalMatrix.n3();
		const unsigned int maxPosition = signalLength - 2;

		for (unsigned int col = r.begin(); col != r.end(); ++col) {
			for (unsigned int row = minRowIdx[col]; row < numRows; ++row) {
#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_USE_ANALYTIC_SIGNAL
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

#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_USE_ANALYTIC_SIGNAL
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
#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_USE_ANALYTIC_SIGNAL
	const Tensor3<std::complex<TFloat>, tbb::cache_aligned_allocator<std::complex<TFloat>>>& signalMatrix;
#else
	const Tensor3<TFloat, tbb::cache_aligned_allocator<TFloat>>& signalMatrix;
#endif
	const std::vector<TFloat>& rxApod;
	const std::vector<StepConfiguration>& stepConfigList;
	const std::vector<unsigned int, tbb::cache_aligned_allocator<unsigned int>>& minRowIdx;
	const Tensor3<TFloat, tbb::cache_aligned_allocator<TFloat>>& delayMatrix;
	tbb::enumerable_thread_specific<ProcessColumnThreadData>& processColumnTLS;
#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_USE_ANALYTIC_SIGNAL
	Matrix<std::complex<TFloat>>& gridValue;
#else
	Matrix<TFloat>& gridValue;
#endif
};



template<typename TFloat>
struct VectorialCombinedTwoMediumImagingProcessor<TFloat>::ProcessColumn {
	void operator()(const tbb::blocked_range<unsigned int>& r) const {
		//LOG_DEBUG << "col = " << r.begin() << " n = " << (r.end() - r.begin());
		ProcessColumnThreadData& local = processColumnTLS.local();

		local.rxSignalSumList.resize(config.numElements);
		const unsigned int signalLength = signalMatrix.n3();
		const unsigned int maxPosition = signalLength - 2;

		for (unsigned int col = r.begin(); col != r.end(); ++col) {
			for (unsigned int row = minRowIdx[col]; row < numRows; ++row) {
#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_USE_ANALYTIC_SIGNAL
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

#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_USE_ANALYTIC_SIGNAL
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
#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_USE_ANALYTIC_SIGNAL
	const Tensor3<std::complex<TFloat>, tbb::cache_aligned_allocator<std::complex<TFloat>>>& signalMatrix;
#else
	const Tensor3<TFloat, tbb::cache_aligned_allocator<TFloat>>& signalMatrix;
#endif
	const std::vector<TFloat>& rxApod;
	const StepConfiguration& stepConfig;
	const std::vector<unsigned int, tbb::cache_aligned_allocator<unsigned int>>& minRowIdx;
	const Tensor3<TFloat, tbb::cache_aligned_allocator<TFloat>>& delayMatrix;
	tbb::enumerable_thread_specific<ProcessColumnThreadData>& processColumnTLS;
#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_USE_ANALYTIC_SIGNAL
	Matrix<std::complex<TFloat>>& gridValue;
#else
	Matrix<TFloat>& gridValue;
#endif
};

} // namespace Lab

#endif // VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_H
