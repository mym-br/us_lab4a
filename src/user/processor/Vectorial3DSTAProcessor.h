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

#ifndef VECTORIAL3DSTAPROCESSOR_H
#define VECTORIAL3DSTAPROCESSOR_H

#include <algorithm> /* for_each */
#include <cmath>
#include <complex>
#include <cstddef> /* std::size_t */
#include <numeric> /* accumulate */
#include <vector>

#include <tbb/enumerable_thread_specific.h>
#include <tbb/tbb.h>

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
#include "Tensor3.h"
#include "Util.h"
#include "XYZValueFactor.h"

// Depends on the signal.
// 1.0 --> pi radian / sample at the original sampling rate.
#define VECTORIAL_3D_STA_PROCESSOR_UPSAMP_FILTER_HALF_TRANSITION_WIDTH 0.2



namespace Lab {

template<typename FloatType>
class Vectorial3DSTAProcessor : public ArrayProcessor<FloatType> {
public:
	Vectorial3DSTAProcessor(
			const SA3DConfiguration<FloatType>& config,
			STAAcquisition<FloatType>& acquisition,
			unsigned int upsamplingFactor,
			AnalyticSignalCoherenceFactorProcessor<FloatType>& coherenceFactor,
			FloatType peakOffset,
			const std::vector<FloatType>& txApod,
			const std::vector<FloatType>& rxApod);
	virtual ~Vectorial3DSTAProcessor() {}

	virtual void process(unsigned int baseElement, Matrix<XYZValueFactor<FloatType>>& gridData);

private:
	struct ThreadData {
		AnalyticSignalCoherenceFactorProcessor<FloatType> coherenceFactor;
		std::vector<std::complex<FloatType>> rxSignalSumList;
		std::vector<FloatType> txDelayList;
		std::vector<FloatType> rxDelayList;
	};

	Vectorial3DSTAProcessor(const Vectorial3DSTAProcessor&) = delete;
	Vectorial3DSTAProcessor& operator=(const Vectorial3DSTAProcessor&) = delete;

	const SA3DConfiguration<FloatType>& config_;
	unsigned int deadZoneSamplesUp_;
	STAAcquisition<FloatType>& acquisition_;
	unsigned int upsamplingFactor_;
	AnalyticSignalCoherenceFactorProcessor<FloatType>& coherenceFactor_;
	Tensor3<std::complex<FloatType>> analyticSignalTensor_;
	typename STAAcquisition<FloatType>::AcquisitionDataType acqData_;
	std::vector<FloatType> tempSignal_;
	FloatType signalOffset_;
	Interpolator<FloatType> interpolator_;
	HilbertEnvelope<FloatType> envelope_;
	bool initialized_;
	std::vector<FloatType> txApod_;
	std::vector<FloatType> rxApod_;
};



template<typename FloatType>
Vectorial3DSTAProcessor<FloatType>::Vectorial3DSTAProcessor(
			const SA3DConfiguration<FloatType>& config,
			STAAcquisition<FloatType>& acquisition,
			unsigned int upsamplingFactor,
			AnalyticSignalCoherenceFactorProcessor<FloatType>& coherenceFactor,
			FloatType peakOffset,
			const std::vector<FloatType>& txApod,
			const std::vector<FloatType>& rxApod)
		: config_(config)
		, deadZoneSamplesUp_((upsamplingFactor * config.samplingFrequency) * 2.0 * config.deadZoneM / config.propagationSpeed)
		, acquisition_(acquisition)
		, upsamplingFactor_(upsamplingFactor)
		, coherenceFactor_(coherenceFactor)
		, initialized_()
		, txApod_(txApod)
		, rxApod_(rxApod)
{
	if (upsamplingFactor_ > 1) {
		interpolator_.prepare(upsamplingFactor_, VECTORIAL_3D_STA_PROCESSOR_UPSAMP_FILTER_HALF_TRANSITION_WIDTH);
	}

	signalOffset_ = (config_.samplingFrequency * upsamplingFactor_) * peakOffset / config_.centerFrequency;

	if (txApod_.size() != config_.activeTxElem.size()) {
		THROW_EXCEPTION(InvalidValueException, "Wrong tx apodization size: " << txApod_.size()
				<< " (should be " << config_.activeTxElem.size() << ").");
	}
	if (rxApod_.size() != config_.activeRxElem.size()) {
		THROW_EXCEPTION(InvalidValueException, "Wrong rx apodization size: " << rxApod_.size()
				<< " (should be " << config_.activeRxElem.size() << ").");
	}
}

template<typename FloatType>
void
Vectorial3DSTAProcessor<FloatType>::process(unsigned int baseElement, Matrix<XYZValueFactor<FloatType>>& gridData)
{
	LOG_DEBUG << "BEGIN ========== Vectorial3DSTAProcessor::process ==========";

	Util::resetValueFactor(gridData.begin(), gridData.end());

	const std::size_t numSignals = config_.activeTxElem.size() * config_.activeRxElem.size();

	IterationCounter::reset(config_.activeTxElem.size());

	// Prepare the signal matrix.
	for (unsigned int iTxElem = 0, txEnd = config_.activeTxElem.size(); iTxElem < txEnd; ++iTxElem) {
		LOG_INFO << "ACQ/PREP txElem: " << config_.activeTxElem[iTxElem];

		acquisition_.execute(baseElement, config_.activeTxElem[iTxElem], acqData_);

		if (!initialized_) {
			const std::size_t signalLength = acqData_.n2() * upsamplingFactor_;
			tempSignal_.resize(signalLength);
			analyticSignalTensor_.resize(
						config_.activeTxElem.size(),
						config_.activeRxElem.size(),
						signalLength);
			LOG_DEBUG << "signalOffset_=" << signalOffset_ << " signalLength=" << signalLength;
			if (deadZoneSamplesUp_ >= signalLength) {
				THROW_EXCEPTION(InvalidValueException,
						"Invalid dead zone: deadZoneSamplesUp_ (" << deadZoneSamplesUp_ <<
						") >= signalLength (" << signalLength << ").");
			}
			initialized_ = true;
		}

		const std::size_t samplesPerChannelLow = acqData_.n2();

		for (unsigned int iRxElem = 0, rxEnd = config_.activeRxElem.size(); iRxElem < rxEnd; ++iRxElem) {
			if (upsamplingFactor_ > 1) {
				interpolator_.interpolate(&acqData_(iRxElem, 0), samplesPerChannelLow, &tempSignal_[0]);
			} else {
				typename Matrix<FloatType>::Dim2Interval interval = acqData_.dim2Interval(iRxElem);
				std::copy(interval.first, interval.second, tempSignal_.begin());
			}

			Util::removeDC(&tempSignal_[0], tempSignal_.size(), deadZoneSamplesUp_);

			envelope_.getAnalyticSignal(&tempSignal_[0], tempSignal_.size(),
							&analyticSignalTensor_(iTxElem, iRxElem, 0));
		}

		IterationCounter::add(1);
	}

	ThreadData threadData;
	threadData.coherenceFactor = coherenceFactor_;
	tbb::enumerable_thread_specific<ThreadData> tls(threadData);

	const FloatType invCT = (config_.samplingFrequency * upsamplingFactor_) / config_.propagationSpeed;
	const std::size_t numRows = gridData.n2();

	IterationCounter::reset(gridData.n1());

	tbb::parallel_for(tbb::blocked_range<std::size_t>(0, gridData.n1()),
	[&, invCT, numRows](const tbb::blocked_range<std::size_t>& r) {
		LOG_DEBUG << "IMG col range start = " << r.begin() << " n = " << (r.end() - r.begin());

		auto& local = tls.local();

		local.rxSignalSumList.resize(config_.activeRxElem.size());
		local.txDelayList.resize(config_.activeTxElem.size());
		local.rxDelayList.resize(config_.activeRxElem.size());

		// For each column:
		for (std::size_t i = r.begin(); i != r.end(); ++i) {
			// For each row:
			for (std::size_t j = 0; j < numRows; ++j) {

				std::fill(local.rxSignalSumList.begin(), local.rxSignalSumList.end(), std::complex<FloatType>(0));
				XYZValueFactor<FloatType>& point = gridData(i, j);

				// Calculate the delays.
				for (unsigned int iTxElem = 0, end = config_.activeTxElem.size(); iTxElem < end; ++iTxElem) {
					const XY<FloatType>& elemPos = config_.txElemPos[baseElement + config_.activeTxElem[iTxElem]];
					local.txDelayList[iTxElem] = Geometry::distance3DZ0(elemPos.x, elemPos.y, point.x, point.y, point.z) * invCT;
				}
				for (unsigned int iRxElem = 0, end = config_.activeRxElem.size(); iRxElem < end; ++iRxElem) {
					const XY<FloatType>& elemPos = config_.rxElemPos[baseElement + config_.activeRxElem[iRxElem]];
					local.rxDelayList[iRxElem] = Geometry::distance3DZ0(elemPos.x, elemPos.y, point.x, point.y, point.z) * invCT;
				}

				for (unsigned int iTxElem = 0, txEnd = config_.activeTxElem.size(); iTxElem < txEnd; ++iTxElem) {
					const FloatType txDelay = local.txDelayList[iTxElem];
					for (unsigned int iRxElem = 0, rxEnd = config_.activeRxElem.size(); iRxElem < rxEnd; ++iRxElem) {
						// Linear interpolation.
						const FloatType delay = signalOffset_ + txDelay + local.rxDelayList[iRxElem];
						const std::size_t delayIdx = static_cast<std::size_t>(delay);
						const FloatType k = delay - delayIdx;
						if (delayIdx + 1U < analyticSignalTensor_.n3()) {
							const std::complex<FloatType>* p = &analyticSignalTensor_(iTxElem, iRxElem, delayIdx);
							local.rxSignalSumList[iRxElem] +=
									txApod_[iTxElem] * rxApod_[iRxElem]
									* ((1 - k) * *p + k * *(p + 1));
						}
					}
				}

				if (local.coherenceFactor.enabled()) {
					point.factor = local.coherenceFactor.calculate(&local.rxSignalSumList[0], local.rxSignalSumList.size());
				}
				point.value = std::abs(std::accumulate(local.rxSignalSumList.begin(), local.rxSignalSumList.end(), std::complex<FloatType>(0)));
			}
		}

		IterationCounter::add(r.end() - r.begin());
	});

	std::for_each(gridData.begin(), gridData.end(),
			Util::MultiplyValueBy<XYZValueFactor<FloatType>, FloatType>(FloatType(1) / numSignals));

	LOG_DEBUG << "END ========== Vectorial3DSTAProcessor::process ==========";
}

} // namespace Lab

#endif // VECTORIAL3DSTAPROCESSOR_H
