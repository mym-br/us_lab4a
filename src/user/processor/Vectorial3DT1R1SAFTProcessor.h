/***************************************************************************
 *  Copyright 2014, 2017, 2018 Marcelo Y. Matuda                           *
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

#ifndef VECTORIAL3DT1R1SAFTPROCESSOR_H
#define VECTORIAL3DT1R1SAFTPROCESSOR_H

#include <algorithm> /* for_each */
#include <cassert>
#include <cmath>
#include <complex>
#include <cstddef> /* std::size_t */
#include <numeric> /* accumulate */
#include <vector>

#include <tbb/blocked_range.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>

#include "ArrayProcessor.h"
#include "CoherenceFactor.h"
#include "Exception.h"
#include "Geometry.h"
#include "HilbertEnvelope.h"
#include "Interpolator.h"
#include "IterationCounter.h"
#include "Log.h"
#include "Matrix.h"
#include "STAAcquisition.h"
#include "SA3DConfiguration.h"
#include "Util.h"
#include "XYZValueFactor.h"



namespace Lab {

template<typename TFloat>
class Vectorial3DT1R1SAFTProcessor : public ArrayProcessor<XYZValueFactor<TFloat>> {
public:
	Vectorial3DT1R1SAFTProcessor(
			const SA3DConfiguration<TFloat>& config,
			STAAcquisition<TFloat>& acquisition,
			unsigned int upsamplingFactor,
			AnalyticSignalCoherenceFactorProcessor<TFloat>& coherenceFactor,
			TFloat peakOffset);
	virtual ~Vectorial3DT1R1SAFTProcessor() = default;

	virtual void prepare(unsigned int baseElement);
	virtual void process(Matrix<XYZValueFactor<TFloat>>& gridData);
private:
	// Depends on the signal.
	// 1.0 --> pi radian / sample at the original sampling rate.
	static constexpr TFloat upsampFilterHalfTransitionWidth = 0.2;

	struct ThreadData {
		AnalyticSignalCoherenceFactorProcessor<TFloat> coherenceFactor;
		std::vector<std::complex<TFloat>> rxSignalSumList;
		std::vector<TFloat> txDelayList;
		std::vector<TFloat> rxDelayList;
	};

	Vectorial3DT1R1SAFTProcessor(const Vectorial3DT1R1SAFTProcessor&) = delete;
	Vectorial3DT1R1SAFTProcessor& operator=(const Vectorial3DT1R1SAFTProcessor&) = delete;
	Vectorial3DT1R1SAFTProcessor(Vectorial3DT1R1SAFTProcessor&&) = delete;
	Vectorial3DT1R1SAFTProcessor& operator=(Vectorial3DT1R1SAFTProcessor&&) = delete;

	const SA3DConfiguration<TFloat>& config_;
	unsigned int deadZoneSamplesUp_;
	STAAcquisition<TFloat>& acquisition_;
	unsigned int upsamplingFactor_;
	AnalyticSignalCoherenceFactorProcessor<TFloat>& coherenceFactor_;
	Matrix<std::complex<TFloat>> analyticSignalMatrix_;
	typename STAAcquisition<TFloat>::AcquisitionDataType acqData_;
	std::vector<TFloat> tempSignal_;
	TFloat signalOffset_;
	Interpolator<TFloat> interpolator_;
	HilbertEnvelope<TFloat> envelope_;
	bool initialized_;
	unsigned int baseElement_;
};



template<typename TFloat>
Vectorial3DT1R1SAFTProcessor<TFloat>::Vectorial3DT1R1SAFTProcessor(
			const SA3DConfiguration<TFloat>& config,
			STAAcquisition<TFloat>& acquisition,
			unsigned int upsamplingFactor,
			AnalyticSignalCoherenceFactorProcessor<TFloat>& coherenceFactor,
			TFloat peakOffset)
		: config_(config)
		, deadZoneSamplesUp_((upsamplingFactor * config.samplingFrequency) * 2.0 * config.deadZoneM / config.propagationSpeed)
		, acquisition_(acquisition)
		, upsamplingFactor_(upsamplingFactor)
		, coherenceFactor_(coherenceFactor)
		, initialized_()
		, baseElement_()
{
	if (upsamplingFactor_ > 1) {
		interpolator_.prepare(upsamplingFactor_, upsampFilterHalfTransitionWidth);
	}

	signalOffset_ = (config_.samplingFrequency * upsamplingFactor_) * peakOffset / config_.centerFrequency;
}

template<typename TFloat>
void
Vectorial3DT1R1SAFTProcessor<TFloat>::prepare(unsigned int baseElement)
{
	baseElement_ = baseElement;
	acquisition_.prepare(baseElement_);
}

template<typename TFloat>
void
Vectorial3DT1R1SAFTProcessor<TFloat>::process(Matrix<XYZValueFactor<TFloat>>& gridData)
{
	LOG_DEBUG << "BEGIN ========== Vectorial3DT1R1SAFTProcessor::process ==========";

	assert(config_.activeTxElem.size() == config_.activeRxElem.size());
	assert(config_.txElemPos.size() == config_.rxElemPos.size());

	Util::resetValueFactor(gridData.begin(), gridData.end());

	const std::size_t numSignals = config_.activeTxElem.size();

	// Prepare the signal matrix.
	for (unsigned int iTxElem = 0, txEnd = config_.activeTxElem.size(); iTxElem < txEnd; ++iTxElem) {
		LOG_INFO << "ACQ/PREP txElem: " << config_.activeTxElem[iTxElem];

		acquisition_.execute(config_.activeTxElem[iTxElem], acqData_);

		if (!initialized_) {
			const std::size_t signalLength = acqData_.n2() * upsamplingFactor_;
			tempSignal_.resize(signalLength);
			analyticSignalMatrix_.resize(config_.activeTxElem.size(), signalLength);
			LOG_DEBUG << "signalOffset_=" << signalOffset_ << " signalLength=" << signalLength;
			if (deadZoneSamplesUp_ >= signalLength) {
				THROW_EXCEPTION(InvalidValueException,
						"Invalid dead zone: deadZoneSamplesUp_ (" << deadZoneSamplesUp_ <<
						") >= signalLength (" << signalLength << ").");
			}
			initialized_ = true;
		}

		const std::size_t samplesPerChannelLow = acqData_.n2();

		if (upsamplingFactor_ > 1) {
			interpolator_.interpolate(&acqData_(0, 0), samplesPerChannelLow, &tempSignal_[0]);
		} else {
			auto range = acqData_.range2(0);
			std::copy(range.begin(), range.end(), tempSignal_.begin());
		}

		Util::removeDC(&tempSignal_[0], tempSignal_.size(), deadZoneSamplesUp_);

		envelope_.getAnalyticSignal(&tempSignal_[0], tempSignal_.size(),
						&analyticSignalMatrix_(iTxElem, 0));
	}

	ThreadData threadData;
	threadData.coherenceFactor = coherenceFactor_;
	tbb::enumerable_thread_specific<ThreadData> tls(threadData);

	const TFloat invCT = (config_.samplingFrequency * upsamplingFactor_) / config_.propagationSpeed;
	const std::size_t numRows = gridData.n2();

	IterationCounter::reset(gridData.n1());

	tbb::parallel_for(tbb::blocked_range<std::size_t>(0, gridData.n1()),
	[&, invCT, numRows](const tbb::blocked_range<std::size_t>& r) {
		auto& local = tls.local();

		local.rxSignalSumList.resize(config_.activeRxElem.size());
		local.txDelayList.resize(config_.activeTxElem.size());
		local.rxDelayList.resize(config_.activeRxElem.size());

		// For each column:
		for (std::size_t i = r.begin(); i != r.end(); ++i) {
			// For each row:
			for (std::size_t j = 0; j < numRows; ++j) {

				std::fill(local.rxSignalSumList.begin(), local.rxSignalSumList.end(), std::complex<TFloat>(0));
				XYZValueFactor<TFloat>& point = gridData(i, j);

				// Calculate the delays.
				for (unsigned int iTxElem = 0, end = config_.activeTxElem.size(); iTxElem < end; ++iTxElem) {
					const XY<TFloat>& elemPos = config_.txElemPos[baseElement_ + config_.activeTxElem[iTxElem]];
					local.txDelayList[iTxElem] = Geometry::distance3DZ0(elemPos.x, elemPos.y, point.x, point.y, point.z) * invCT;
				}
				for (unsigned int iRxElem = 0, end = config_.activeRxElem.size(); iRxElem < end; ++iRxElem) {
					const XY<TFloat>& elemPos = config_.rxElemPos[baseElement_ + config_.activeRxElem[iRxElem]];
					local.rxDelayList[iRxElem] = Geometry::distance3DZ0(elemPos.x, elemPos.y, point.x, point.y, point.z) * invCT;
				}

				for (unsigned int iTxElem = 0, txEnd = config_.activeTxElem.size(); iTxElem < txEnd; ++iTxElem) {
					const TFloat txDelay = local.txDelayList[iTxElem];
					// Linear interpolation.
					const TFloat delay = signalOffset_ + txDelay + local.rxDelayList[iTxElem];
					const std::size_t delayIdx = static_cast<std::size_t>(delay);
					const TFloat k = delay - delayIdx;
					if (delayIdx + 1U < analyticSignalMatrix_.n2()) {
						const std::complex<TFloat>* p = &analyticSignalMatrix_(iTxElem, delayIdx);
						local.rxSignalSumList[iTxElem] += (1 - k) * *p + k * *(p + 1);
					}
				}

				if (local.coherenceFactor.enabled()) {
					point.factor = local.coherenceFactor.calculate(&local.rxSignalSumList[0], local.rxSignalSumList.size());
				}
				point.value = std::abs(std::accumulate(local.rxSignalSumList.begin(), local.rxSignalSumList.end(), std::complex<TFloat>(0)));
			}
		}

		IterationCounter::add(r.end() - r.begin());
	});

	std::for_each(gridData.begin(), gridData.end(),
			Util::MultiplyValueBy<XYZValueFactor<TFloat>, TFloat>(TFloat(1) / numSignals));

	LOG_DEBUG << "END ========== Vectorial3DT1R1SAFTProcessor::process ==========";
}

} // namespace Lab

#endif // VECTORIAL3DT1R1SAFTPROCESSOR_H
