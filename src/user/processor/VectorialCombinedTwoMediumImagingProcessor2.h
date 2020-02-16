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

#ifndef VECTORIALCOMBINEDTWOMEDIUMIMAGINGPROCESSOR2_H
#define VECTORIALCOMBINEDTWOMEDIUMIMAGINGPROCESSOR2_H

#include <algorithm> /* copy, fill */
#include <cmath> /* ceil, sqrt */
#include <cstddef> /* std::size_t */
#include <memory>
#include <vector>

#include <tbb/enumerable_thread_specific.h>
#include <tbb/tbb.h>

#include "ArrayGeometry.h"
#include "CoherenceFactor.h"
#include "Exception.h"
#include "FermatPrinciple.h"
#include "HilbertEnvelope.h"
#include "Interpolator.h"
#include "Log.h"
#include "Matrix.h"
#include "SIMD.h"
#include "Tensor3.h"
#include "TwoMediumSTAConfiguration.h"
#include "Util.h"
#include "XZ.h"

//TODO: remove test (timer)
#include "Timer.h"

// Depends on the signal.
// 1.0 --> pi radian / sample at the original sampling rate.
#define VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_2_UPSAMP_FILTER_HALF_TRANSITION_WIDTH 0.2

#define VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_2_USE_ANALYTIC_SIGNAL 1
#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_2_USE_ANALYTIC_SIGNAL
# include <complex>
# include "XZComplexValue.h"
#else
# include "XZValue.h"
#endif



namespace Lab {

// Faster version of VectorialCombinedTwoMediumImagingProcessor.
// The grid must be rectangular.
//
// Problem: the calculation for the Fermat's principle is executed more than once for some
// grid point-element combinations, due to the superposition of the groups.
// Sometimes it is executed three times for the combination.
// This is solved in VectorialCombinedTwoMediumImagingProcessor3.
template<typename TFloat>
class VectorialCombinedTwoMediumImagingProcessor2 {
public:
	struct StepConfiguration {
		unsigned int baseElemIdx;
		unsigned int baseElem;
		unsigned int firstTxElem;
		unsigned int lastTxElem;
	};

	VectorialCombinedTwoMediumImagingProcessor2(
			const TwoMediumSTAConfiguration<TFloat>& config,
			std::vector<Tensor3<TFloat>>& acqDataList,
			unsigned int upsamplingFactor,
#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_2_USE_ANALYTIC_SIGNAL
			AnalyticSignalCoherenceFactorProcessor<TFloat>& coherenceFactor,
#else
			CoherenceFactorProcessor<TFloat>& coherenceFactor,
#endif
			TFloat maxFermatBlockSize,
			TFloat peakOffset,
			unsigned int ascanStartOffset);
	~VectorialCombinedTwoMediumImagingProcessor2() = default;

	void process(
		const std::vector<StepConfiguration>& stepConfigList,
		const std::vector<XZ<TFloat>>& interfacePointList,
		const std::vector<TFloat>& rxApod,
#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_2_USE_ANALYTIC_SIGNAL
		Matrix<XZComplexValue<TFloat>>& gridData
#else
		Matrix<XZValue<TFloat>>& gridData
#endif
		);
private:
	struct PrepareData;
	struct PrepareDataThreadData {
		Interpolator<TFloat> interpolator;
		HilbertEnvelope<TFloat> envelope;
		std::vector<TFloat> signal;
	};

	struct ProcessColumn;
	struct ProcessColumnThreadData {
#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_2_USE_ANALYTIC_SIGNAL
		AnalyticSignalCoherenceFactorProcessor<TFloat> coherenceFactor;
		std::vector<std::complex<TFloat>> rxSignalSumList;
#else
		CoherenceFactorProcessor<TFloat> coherenceFactor;
		std::vector<TFloat> rxSignalSumList;
#endif
		std::vector<TFloat> delayList;
		std::vector<unsigned int> lastInterfaceIdxList;
	};

	VectorialCombinedTwoMediumImagingProcessor2(const VectorialCombinedTwoMediumImagingProcessor2&) = delete;
	VectorialCombinedTwoMediumImagingProcessor2& operator=(const VectorialCombinedTwoMediumImagingProcessor2&) = delete;
	VectorialCombinedTwoMediumImagingProcessor2(VectorialCombinedTwoMediumImagingProcessor2&&) = delete;
	VectorialCombinedTwoMediumImagingProcessor2& operator=(VectorialCombinedTwoMediumImagingProcessor2&&) = delete;

	const TwoMediumSTAConfiguration<TFloat>& config_;
	std::vector<Tensor3<TFloat>>& acqDataList_;
	unsigned int upsamplingFactor_;
#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_2_USE_ANALYTIC_SIGNAL
	AnalyticSignalCoherenceFactorProcessor<TFloat>& coherenceFactor_;
#else
	CoherenceFactorProcessor<TFloat>& coherenceFactor_;
#endif
	TFloat maxFermatBlockSize_;
	const TFloat lambda2_;
	std::size_t signalLength_;
#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_2_USE_ANALYTIC_SIGNAL
	Tensor3<std::complex<TFloat>> signalMatrix_;
#else
	Tensor3<TFloat> signalMatrix_;
#endif
	TFloat signalOffset_;
	std::vector<unsigned int> minRowIdx_; // for each column
	std::vector<TFloat> xArray_;
	Matrix<TFloat> medium1DelayMatrix_; // (interface_idx, element)
	std::unique_ptr<tbb::enumerable_thread_specific<PrepareDataThreadData>> prepareDataTLS_;
	std::unique_ptr<tbb::enumerable_thread_specific<ProcessColumnThreadData>> processColumnTLS_;
};



template<typename TFloat>
VectorialCombinedTwoMediumImagingProcessor2<TFloat>::VectorialCombinedTwoMediumImagingProcessor2(
			const TwoMediumSTAConfiguration<TFloat>& config,
			std::vector<Tensor3<TFloat>>& acqDataList,
			unsigned int upsamplingFactor,
#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_2_USE_ANALYTIC_SIGNAL
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
		prepareDataThreadData.interpolator.prepare(upsamplingFactor_, VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_2_UPSAMP_FILTER_HALF_TRANSITION_WIDTH);
	}
	prepareDataThreadData.signal.resize(signalLength_);
	prepareDataTLS_ = std::make_unique<tbb::enumerable_thread_specific<PrepareDataThreadData>>(prepareDataThreadData);

	ProcessColumnThreadData processColumnThreadData;
	processColumnThreadData.coherenceFactor = coherenceFactor_;
	processColumnTLS_ = std::make_unique<tbb::enumerable_thread_specific<ProcessColumnThreadData>>(processColumnThreadData);
}

template<typename TFloat>
void
VectorialCombinedTwoMediumImagingProcessor2<TFloat>::process(
							const std::vector<StepConfiguration>& stepConfigList,
							const std::vector<XZ<TFloat>>& interfacePointList,
							const std::vector<TFloat>& rxApod,
#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_2_USE_ANALYTIC_SIGNAL
							Matrix<XZComplexValue<TFloat>>& gridData
#else
							Matrix<XZValue<TFloat>>& gridData
#endif
							)
{
	//LOG_DEBUG << "BEGIN ========== VectorialCombinedTwoMediumImagingProcessor2::process ==========";

	const std::size_t samplesPerChannelLow = acqDataList_[0].n3();

	// Clear the values in gridData.
	for (auto& item : gridData) {
		item.value = 0;
	}

	XZ<TFloat> p1 = interfacePointList[0];
	XZ<TFloat> p2 = interfacePointList[1];
	const TFloat dx = p2.x - p1.x;
	const TFloat dz = p2.z - p1.z;
	const TFloat r = std::sqrt(dx * dx + dz * dz);
	const unsigned int fermatBlockSize = FermatPrinciple::calcBlockSizeForTwoStepSearch(interfacePointList.size(), r, lambda2_, maxFermatBlockSize_);
	LOG_DEBUG << "fermatBlockSize: " << fermatBlockSize;

	minRowIdx_.resize(gridData.n1() /* number of columns */);

	//unsigned int numValidPixels = 0;
	const TFloat zStepGrid = gridData(0, 1).z - gridData(0, 0).z;
	for (unsigned int col = 0; col < gridData.n1(); ++col) {
		auto& point = gridData(col, 0);

		// Find the z coordinate of the interface.
		TFloat zIdxMin;
		unsigned int idxMin;
		FermatPrinciple::findNearestPointInXInTwoSteps(
				fermatBlockSize,
				interfacePointList,
				point.x,
				zIdxMin, idxMin);

		if (zIdxMin <= gridData(col, 0).z) {
			minRowIdx_[col] = 1;
		} else if (zIdxMin >= gridData(col, gridData.n2() - 1).z) {
			minRowIdx_[col] = gridData.n2(); // after the last index
		} else {
			minRowIdx_[col] = static_cast<unsigned int>(std::ceil((zIdxMin - gridData(col, 0).z) / zStepGrid)) + 1;
		}
		//numValidPixels += gridData.n2() - minRowIdx_[col];
	}
	//LOG_DEBUG << "numValidPixels: " << numValidPixels;

	const TFloat c2ByC1 = config_.propagationSpeed2 / config_.propagationSpeed1;

	for (const auto& stepConfig : stepConfigList) {
		signalMatrix_.resize(
				stepConfig.lastTxElem - stepConfig.firstTxElem + 1,
				config_.numElements,
				signalLength_);

		// Prepare the signal matrix.
		Timer t0;
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
		LOG_DEBUG << "PREP time = " << t0.getTime();

		const TFloat xOffset = (stepConfig.baseElem + 0.5 * (config_.numElements - 1) - 0.5 * (config_.numElementsMux - 1)) * config_.pitch;
		ArrayGeometry::getElementX2D(config_.numElements, config_.pitch, xOffset, xArray_);

		Timer tMedium1DelayMatrix;
		medium1DelayMatrix_.resize(interfacePointList.size(), config_.numElements);
		//LOG_DEBUG << "interfacePointList.size(): " << interfacePointList.size();
		for (unsigned int i = 0; i < interfacePointList.size(); ++i) {
			for (unsigned int elem = 0; elem < config_.numElements; ++elem) {
				const XZ<TFloat>& ifPoint = interfacePointList[i];
//				const TFloat dx = ifPoint.x - xArray_[elem];
//				const TFloat dz = ifPoint.z;
//				medium1DelayMatrix_(i, elem) = SIMD::calcDistance(dx, dz) * c2ByC1;
//				medium1DelayMatrix_(i, elem) = std::sqrt(dx * dx + dz * dz) * c2ByC1;
				medium1DelayMatrix_(i, elem) = SIMD::calcDistance(xArray_[elem], 0, ifPoint.x, ifPoint.z) * c2ByC1;
			}
		}
		LOG_DEBUG << "PREP tMedium1DelayMatrix time = " << tMedium1DelayMatrix.getTime();

		Timer t1;
		ProcessColumn processColumnOp = {
			gridData.n2(),
			config_,
			signalOffset_,
			signalMatrix_,
			config_.samplingFrequency * upsamplingFactor_,
			config_.samplingFrequency * upsamplingFactor_ / config_.propagationSpeed2,
			1 / config_.propagationSpeed1,
			1 / config_.propagationSpeed2,
			fermatBlockSize,
			interfacePointList,
			xArray_,
			rxApod,
			stepConfig,
			minRowIdx_,
			medium1DelayMatrix_,
			*processColumnTLS_,
			gridData
		};
		tbb::parallel_for(tbb::blocked_range<unsigned int>(0, gridData.n1()), processColumnOp);
		//processColumnOp(tbb::blocked_range<unsigned int>(0, gridData.n1())); // single-thread
		LOG_DEBUG << "PARTIAL PROC time = " << t1.getTime();
	}

	//LOG_DEBUG << "END ========== VectorialCombinedTwoMediumImagingProcessor2::process ==========";
}



template<typename TFloat>
struct VectorialCombinedTwoMediumImagingProcessor2<TFloat>::PrepareData {
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

#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_2_USE_ANALYTIC_SIGNAL
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
#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_2_USE_ANALYTIC_SIGNAL
	Tensor3<std::complex<TFloat>>& signalMatrix;
#else
	Tensor3<TFloat>& signalMatrix;
#endif
};



template<typename TFloat>
struct VectorialCombinedTwoMediumImagingProcessor2<TFloat>::ProcessColumn {
	void operator()(const tbb::blocked_range<unsigned int>& r) const {
		//LOG_DEBUG << "col = " << r.begin() << " n = " << (r.end() - r.begin());
		ProcessColumnThreadData& local = processColumnTLS.local();

		local.rxSignalSumList.resize(config.numElements);
		local.delayList.resize(config.numElements);
		local.lastInterfaceIdxList.resize(config.numElements);

		for (unsigned int col = r.begin(); col != r.end(); ++col) {
			for (unsigned int row = minRowIdx[col]; row < numRows; ++row) {
#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_2_USE_ANALYTIC_SIGNAL
				std::fill(local.rxSignalSumList.begin(), local.rxSignalSumList.end(), std::complex<TFloat>(0));
#else
				std::fill(local.rxSignalSumList.begin(), local.rxSignalSumList.end(), TFloat(0));
#endif
				auto& point = gridData(col, row);

				// Calculate the delays.
				if (row == minRowIdx[col]) { // the first row
					for (unsigned int elem = 0; elem < config.numElements; ++elem) {
						// Fermat's principle. Find the fastest path.
						TFloat tMin;
						unsigned int idxMin;
//						FermatPrinciple::findMinTimeInTwoSteps(
//								fermatBlockSize,
//								config.propagationSpeed1, config.propagationSpeed2,
//								interfacePointList,
//								xArray[elem], TFloat(0), point.x, point.z,
//								tMin, idxMin);
						FermatPrinciple::findMinTimeInTwoSteps2(
								fermatBlockSize,
								invC1, invC2,
								interfacePointList,
								xArray[elem], TFloat(0), point.x, point.z,
								tMin, idxMin);
						local.delayList[elem] = tMin * fs;
						local.lastInterfaceIdxList[elem] = idxMin;
					}
				} else {
					for (unsigned int elem = 0; elem < config.numElements; ++elem) {
						unsigned int idxMin = local.lastInterfaceIdxList[elem];
						TFloat tC2Min;
						{
							const XZ<TFloat>& ifPoint = interfacePointList[idxMin];
//							const TFloat dx2 = point.x - ifPoint.x;
//							const TFloat dz2 = point.z - ifPoint.z;
//							tC2Min = medium1DelayMatrix(idxMin, elem) + SIMD::calcDistance(dx2, dz2);
//							tC2Min = medium1DelayMatrix(idxMin, elem) + std::sqrt(dx2 * dx2 + dz2 * dz2);
							tC2Min = medium1DelayMatrix(idxMin, elem) + SIMD::calcDistance(ifPoint.x, ifPoint.z, point.x, point.z);
						}
						for (int idxSearch = static_cast<int>(idxMin) + 1, end = interfacePointList.size(); idxSearch < end; ++idxSearch) {
							const XZ<TFloat>& ifPoint = interfacePointList[idxSearch];
//							const TFloat dx2 = point.x - ifPoint.x;
//							const TFloat dz2 = point.z - ifPoint.z;
//							const TFloat tC2 = medium1DelayMatrix(idxSearch, elem) + SIMD::calcDistance(dx2, dz2);
//							const TFloat tC2 = medium1DelayMatrix(idxSearch, elem) + std::sqrt(dx2 * dx2 + dz2 * dz2);
							const TFloat tC2 = medium1DelayMatrix(idxSearch, elem) + SIMD::calcDistance(ifPoint.x, ifPoint.z, point.x, point.z);
							if (tC2 >= tC2Min) {
								break;
							} else {
								tC2Min = tC2;
								idxMin = idxSearch;
							}
						}
						if (idxMin == local.lastInterfaceIdxList[elem]) { // if the previous search was not successful
							for (int idxSearch = static_cast<int>(idxMin) - 1; idxSearch >= 0; --idxSearch) {
								const XZ<TFloat>& ifPoint = interfacePointList[idxSearch];
//								const TFloat dx2 = point.x - ifPoint.x;
//								const TFloat dz2 = point.z - ifPoint.z;
//								const TFloat tC2 = medium1DelayMatrix(idxSearch, elem) + SIMD::calcDistance(dx2, dz2);
//								const TFloat tC2 = medium1DelayMatrix(idxSearch, elem) + std::sqrt(dx2 * dx2 + dz2 * dz2);
								const TFloat tC2 = medium1DelayMatrix(idxSearch, elem) + SIMD::calcDistance(ifPoint.x, ifPoint.z, point.x, point.z);
								if (tC2 >= tC2Min) {
									break;
								} else {
									tC2Min = tC2;
									idxMin = idxSearch;
								}
							}
						}

//						unsigned int diff = (idxMin > local.lastInterfaceIdxList[elem]) ?
//									idxMin - local.lastInterfaceIdxList[elem] :
//									local.lastInterfaceIdxList[elem] - idxMin;
//						if (diff > 1) {
//							LOG_DEBUG << "########## DIFF " << diff << " idxMin: " << idxMin << " col: " << col <<
//								" row - minRowIdx[col]: " << row - minRowIdx[col];
//						}

						local.delayList[elem] = tC2Min * fsInvC2;
						local.lastInterfaceIdxList[elem] = idxMin;
					}
				}

				for (unsigned int txElem = stepConfig.firstTxElem; txElem <= stepConfig.lastTxElem; ++txElem) {
					const TFloat txDelay = local.delayList[txElem];
					const unsigned int txElemIdx = txElem - stepConfig.firstTxElem;
					const TFloat txOffset = signalOffset + txDelay;
					const auto* p = &signalMatrix(txElemIdx, 0 /* rxElem */, 0);
					for (unsigned int rxElem = 0; rxElem < config.numElements; ++rxElem, p += signalMatrix.n3()) {
						// Linear interpolation.
						const TFloat position = txOffset + local.delayList[rxElem];
						if (position >= 0.0f) {
							const unsigned int positionIdx = static_cast<unsigned int>(position);
							if (positionIdx < signalMatrix.n3() - 1) {
								const TFloat k = position - positionIdx;
								const auto* pLocal = p + positionIdx;
								local.rxSignalSumList[rxElem] += (1 - k) * *pLocal + k * *(pLocal + 1);
							}
						}
					}
				}

#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_2_USE_ANALYTIC_SIGNAL
				std::complex<TFloat> sum = 0.0;
#else
				TFloat sum = 0.0;
#endif
				for (unsigned int rxElem = 0; rxElem < config.numElements; ++rxElem) {
					sum += local.rxSignalSumList[rxElem] * rxApod[rxElem];
				}
				if (local.coherenceFactor.enabled()) {
					point.value += sum * local.coherenceFactor.calculate(&local.rxSignalSumList[0], local.rxSignalSumList.size());
				} else {
					point.value += sum;
				}
			}
		}
	}

	// Use only small types or references/pointers, because this object may be copied many times.
	const std::size_t numRows;
	const TwoMediumSTAConfiguration<TFloat>& config;
	const TFloat signalOffset;
#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_2_USE_ANALYTIC_SIGNAL
	const Tensor3<std::complex<TFloat>>& signalMatrix;
#else
	const Tensor3<TFloat>& signalMatrix;
#endif
	const TFloat fs;
	const TFloat fsInvC2;
	const TFloat invC1;
	const TFloat invC2;
	const unsigned int fermatBlockSize;
	const std::vector<XZ<TFloat>>& interfacePointList;
	const std::vector<TFloat>& xArray;
	const std::vector<TFloat>& rxApod;
	const StepConfiguration& stepConfig;
	const std::vector<unsigned int>& minRowIdx;
	const Matrix<TFloat>& medium1DelayMatrix;
	tbb::enumerable_thread_specific<ProcessColumnThreadData>& processColumnTLS;
#ifdef VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_2_USE_ANALYTIC_SIGNAL
	Matrix<XZComplexValue<TFloat>>& gridData;
#else
	Matrix<XZValue<TFloat>>& gridData;
#endif
};

} // namespace Lab

#endif // VECTORIALCOMBINEDTWOMEDIUMIMAGINGPROCESSOR2_H
