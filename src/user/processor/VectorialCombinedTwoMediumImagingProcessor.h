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

#ifndef VECTORIALCOMBINEDTWOMEDIUMIMAGINGPROCESSOR_H
#define VECTORIALCOMBINEDTWOMEDIUMIMAGINGPROCESSOR_H

#include <algorithm> /* copy, fill */
#include <cmath> /* sqrt */
#include <complex>
#include <cstddef> /* std::size_t */
#include <numeric> /* accumulate */
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
#include "STAAcquisition.h"
#include "Tensor3.h"
#include "TwoMediumSTAConfiguration.h"
#include "Util.h"
#include "XZ.h"
#include "XZComplexValueFactor.h"

//TODO: remove test (timer)
#include "Timer.h"

// Depends on the signal.
// 1.0 --> pi radian / sample at the original sampling rate.
#define VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_UPSAMP_FILTER_HALF_TRANSITION_WIDTH 0.2



namespace Lab {

template<typename TFloat>
class VectorialCombinedTwoMediumImagingProcessor {
public:
	struct StepConfiguration {
		unsigned int baseElem;
		unsigned int firstTxElem;
		unsigned int lastTxElem;
		unsigned int firstRxElem;
		unsigned int lastRxElem;
	};

	VectorialCombinedTwoMediumImagingProcessor(
			const TwoMediumSTAConfiguration<TFloat>& config,
			STAAcquisition<TFloat>& acquisition,
			unsigned int upsamplingFactor,
			AnalyticSignalCoherenceFactorProcessor<TFloat>& coherenceFactor,
			TFloat maxFermatBlockSize, TFloat peakOffset);
	~VectorialCombinedTwoMediumImagingProcessor() = default;

	void process(const std::vector<StepConfiguration>& stepConfigList,
			const std::vector<XZ<TFloat>>& interfacePointList,
			const std::vector<TFloat>& rxApod,
			Matrix<XZComplexValueFactor<TFloat>>& gridData);

private:
	struct ThreadData {
		AnalyticSignalCoherenceFactorProcessor<TFloat> coherenceFactor;
		std::vector<std::complex<TFloat>> rxSignalSumList;
		std::vector<TFloat> delayList;
	};

	struct ProcessColumn;

	VectorialCombinedTwoMediumImagingProcessor(const VectorialCombinedTwoMediumImagingProcessor&) = delete;
	VectorialCombinedTwoMediumImagingProcessor& operator=(const VectorialCombinedTwoMediumImagingProcessor&) = delete;
	VectorialCombinedTwoMediumImagingProcessor(VectorialCombinedTwoMediumImagingProcessor&&) = delete;
	VectorialCombinedTwoMediumImagingProcessor& operator=(VectorialCombinedTwoMediumImagingProcessor&&) = delete;

	const TwoMediumSTAConfiguration<TFloat>& config_;
	unsigned int deadZoneSamplesUp_;
	STAAcquisition<TFloat>& acquisition_;
	unsigned int upsamplingFactor_;
	AnalyticSignalCoherenceFactorProcessor<TFloat>& coherenceFactor_;
	TFloat maxFermatBlockSize_;
	const TFloat lambda2_;
	std::size_t signalLength_;
	Tensor3<std::complex<TFloat>> analyticSignalMatrix_;
	typename STAAcquisition<TFloat>::AcquisitionDataType acqData_;
	std::vector<TFloat> signal_;
	TFloat signalOffset_;
	Interpolator<TFloat> interpolator_;
	HilbertEnvelope<TFloat> envelope_;
};



template<typename TFloat>
VectorialCombinedTwoMediumImagingProcessor<TFloat>::VectorialCombinedTwoMediumImagingProcessor(
			const TwoMediumSTAConfiguration<TFloat>& config,
			STAAcquisition<TFloat>& acquisition,
			unsigned int upsamplingFactor,
			AnalyticSignalCoherenceFactorProcessor<TFloat>& coherenceFactor,
			TFloat maxFermatBlockSize, TFloat peakOffset)
		: config_(config)
		, deadZoneSamplesUp_((upsamplingFactor * config.samplingFrequency) * 2.0 * config.deadZoneM / config.propagationSpeed1)
		, acquisition_(acquisition)
		, upsamplingFactor_(upsamplingFactor)
		, coherenceFactor_(coherenceFactor)
		, maxFermatBlockSize_(maxFermatBlockSize)
		, lambda2_(config_.propagationSpeed2 / config_.centerFrequency)
{
	signalOffset_ = (config_.samplingFrequency * upsamplingFactor_) * peakOffset / config_.centerFrequency;
	LOG_DEBUG << "signalOffset_=" << signalOffset_;

	// This acquisition is done to obtain the ascan length.
	acquisition_.prepare(0);
	acquisition_.execute(0, acqData_);
	signalLength_ = acqData_.n2() * upsamplingFactor_;

	signal_.resize(signalLength_);

	if (upsamplingFactor_ > 1) {
		interpolator_.prepare(upsamplingFactor_, VECTORIAL_COMBINED_TWO_MEDIUM_IMAGING_PROCESSOR_UPSAMP_FILTER_HALF_TRANSITION_WIDTH);
	}

	if (deadZoneSamplesUp_ >= signalLength_) {
		THROW_EXCEPTION(InvalidValueException, "Invalid dead zone: deadZoneSamplesUp_ (" << deadZoneSamplesUp_ <<
							") >= signalLength (" << signalLength_ << ").");
	}
}

template<typename TFloat>
void
VectorialCombinedTwoMediumImagingProcessor<TFloat>::process(const std::vector<StepConfiguration>& stepConfigList,
								const std::vector<XZ<TFloat>>& interfacePointList,
								const std::vector<TFloat>& rxApod,
								Matrix<XZComplexValueFactor<TFloat>>& gridData)
{
	LOG_DEBUG << "BEGIN ========== VectorialCombinedTwoMediumImagingProcessor::process ==========";

	// Clear the values in gridData.
	for (auto& item : gridData) {
		item.value = 0.0;
		item.factor = 1.0;
	}

	XZ<TFloat> p1 = interfacePointList[0];
	XZ<TFloat> p2 = interfacePointList[1];
	const TFloat dx = p2.x - p1.x;
	const TFloat dz = p2.z - p1.z;
	const TFloat r = std::sqrt(dx * dx + dz * dz);
	const unsigned int fermatBlockSize = FermatPrinciple::calcBlockSizeForTwoStepSearch(interfacePointList.size(), r, lambda2_, maxFermatBlockSize_);
	LOG_DEBUG << "fermatBlockSize: " << fermatBlockSize;

	std::vector<TFloat> xArray;
	ThreadData threadData;
	threadData.coherenceFactor = coherenceFactor_;
	tbb::enumerable_thread_specific<ThreadData> processColumnTLS(threadData);

	for (const auto& stepConfig : stepConfigList) {
		analyticSignalMatrix_.resize(
					stepConfig.lastTxElem - stepConfig.firstTxElem + 1,
					stepConfig.lastRxElem - stepConfig.firstRxElem + 1,
					signalLength_);

		//TODO: remove test (timer)
		Timer t0;

		acquisition_.prepare(stepConfig.baseElem);

		// Prepare the signal matrix.
		for (unsigned int txElem = stepConfig.firstTxElem; txElem <= stepConfig.lastTxElem; ++txElem) {
			acquisition_.execute(txElem, acqData_);
			const std::size_t samplesPerChannelLow = acqData_.n2();

			for (unsigned int rxElem = stepConfig.firstRxElem; rxElem <= stepConfig.lastRxElem; ++rxElem) {

				if (upsamplingFactor_ > 1) {
					// Interpolate the signal.
					interpolator_.interpolate(&acqData_(rxElem, 0), samplesPerChannelLow, &signal_[0]);
				} else {
					auto range = acqData_.range2(rxElem);
					std::copy(range.begin(), range.end(), signal_.begin());
				}

				Util::removeDC(&signal_[0], signal_.size(), deadZoneSamplesUp_);

				// Obtain the analytic signal.
				envelope_.getAnalyticSignal(&signal_[0], signal_.size(),
						&analyticSignalMatrix_(txElem - stepConfig.firstTxElem, rxElem - stepConfig.firstRxElem, 0));

			}

			//LOG_DEBUG << "PREP txElem = " << txElem;
		}

		//TODO: remove test (timer)
		LOG_DEBUG << "PREP time = " << t0.getTime();

		const TFloat xOffset = (stepConfig.baseElem + 0.5 * (config_.numElements - 1) - 0.5 * (config_.numElementsMux - 1)) * config_.pitch;
		ArrayGeometry::getElementX2D(config_.numElements, config_.pitch, xOffset, xArray);

		ProcessColumn op = {
			gridData.n2(),
			config_,
			signalOffset_,
			analyticSignalMatrix_,
			config_.samplingFrequency * upsamplingFactor_,
			1 / config_.propagationSpeed1,
			config_.samplingFrequency * upsamplingFactor_ / config_.propagationSpeed1,
			1 / config_.propagationSpeed2,
			fermatBlockSize,
			interfacePointList,
			xArray,
			rxApod,
			stepConfig,
			processColumnTLS,
			gridData
		};
		tbb::parallel_for(tbb::blocked_range<unsigned int>(0, gridData.n1()), op);
	}

	LOG_DEBUG << "END ========== VectorialCombinedTwoMediumImagingProcessor::process ==========";
}



template<typename TFloat>
struct VectorialCombinedTwoMediumImagingProcessor<TFloat>::ProcessColumn {
	void operator()(const tbb::blocked_range<unsigned int>& r) const {
		//LOG_DEBUG << "col = " << r.begin() << " n = " << (r.end() - r.begin());

		ThreadData& data = processColumnTLS.local();

		data.rxSignalSumList.resize(stepConfig.lastRxElem - stepConfig.firstRxElem + 1);
		data.delayList.resize(config.numElements);

		for (unsigned int i = r.begin(); i != r.end(); ++i) { // columns
			for (unsigned int j = 0; j < numRows; ++j) {

				std::fill(data.rxSignalSumList.begin(), data.rxSignalSumList.end(), std::complex<TFloat>(0));
				XZComplexValueFactor<TFloat>& point = gridData(i, j);

				// Find the z coordinate of the interface.
				TFloat zIdxMin;
				unsigned int idxMin;
				FermatPrinciple::findNearestPointInXInTwoSteps(
						fermatBlockSize,
						interfacePointList,
						point.x,
						zIdxMin, idxMin);
				const TFloat zInterface = zIdxMin;
				const TFloat dz = point.z - zInterface;

				if (dz > 0.0) { // the point is above the interface
					// Calculate the delays.
					for (unsigned int elem = 0; elem < config.numElements; ++elem) {
						// Fermat's principle. Find the fastest path.
						TFloat tMin;
						unsigned int idxMin;
						FermatPrinciple::findMinTimeInTwoSteps(
								fermatBlockSize,
								config.propagationSpeed1, config.propagationSpeed2,
								interfacePointList,
								xArray[elem],
								TFloat(0), point.x, point.z,
								tMin, idxMin);
						if (idxMin == 0 || idxMin == interfacePointList.size() - 1) {
							// Ignore the points for which the interface is too short.
							data.delayList[elem] = analyticSignalMatrix.n3(); // mark as invalid
						} else {
							data.delayList[elem] = tMin * fs;
						}
					}

					for (unsigned int txElem = stepConfig.firstTxElem; txElem <= stepConfig.lastTxElem; ++txElem) {
						const TFloat txDelay = data.delayList[txElem];
						for (unsigned int rxElem = stepConfig.firstRxElem; rxElem <= stepConfig.lastRxElem; ++rxElem) {
							// Linear interpolation.
							const TFloat delay = signalOffset + txDelay + data.delayList[rxElem];
							const std::size_t delayIdx = static_cast<std::size_t>(delay);
							const TFloat k = delay - delayIdx;
							if (delayIdx < analyticSignalMatrix.n3() - 1) {
								const std::complex<TFloat>* p = &analyticSignalMatrix(
													txElem - stepConfig.firstTxElem,
													rxElem - stepConfig.firstRxElem,
													delayIdx);
								data.rxSignalSumList[rxElem - stepConfig.firstRxElem] += rxApod[rxElem] * ((1 - k) * *p + k * *(p + 1));
							}
						}
					}

					//point.factor = 1;//TODO: remove factor?
					const std::complex<TFloat> sum = std::accumulate(data.rxSignalSumList.begin(), data.rxSignalSumList.end(), std::complex<TFloat>(0));
					if (data.coherenceFactor.enabled()) {
						point.value += sum * data.coherenceFactor.calculate(&data.rxSignalSumList[0], data.rxSignalSumList.size());
					} else {
						point.value += sum;
					}
				}
			}
		}
	}

	// Use only small types or references/pointers, because this object may be copied many times.
	const std::size_t numRows;
	const TwoMediumSTAConfiguration<TFloat>& config;
	const TFloat signalOffset;
	const Tensor3<std::complex<TFloat>>& analyticSignalMatrix;
	const TFloat fs;
	const TFloat invC1;
	const TFloat invC1T;
	const TFloat invC2;
	const unsigned int fermatBlockSize;
	const std::vector<XZ<TFloat>>& interfacePointList;
	const std::vector<TFloat>& xArray;
	const std::vector<TFloat>& rxApod;
	const StepConfiguration& stepConfig;
	tbb::enumerable_thread_specific<ThreadData>& processColumnTLS;
	Matrix<XZComplexValueFactor<TFloat>>& gridData;
};

} // namespace Lab

#endif // VECTORIALCOMBINEDTWOMEDIUMIMAGINGPROCESSOR_H
