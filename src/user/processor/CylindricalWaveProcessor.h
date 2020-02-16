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

#ifndef CYLINDRICALWAVEPROCESSOR_H
#define CYLINDRICALWAVEPROCESSOR_H

#include <algorithm> /* copy, fill, reverse, reverse_copy */
#include <cmath> /* abs, sqrt */
#include <complex>
#include <cstddef> /* std::size_t */
#include <numeric> /* accumulate */
#include <vector>

#include <tbb/enumerable_thread_specific.h>
#include <tbb/tbb.h>

#include "CoherenceFactor.h"
#include "Exception.h"
#include "FFTWFilter.h"
#include "HilbertEnvelope.h"
#include "Interpolator.h"
#include "Log.h"
#include "Matrix.h"
#include "STAAcquisition.h"
#include "STAConfiguration.h"
#include "Util.h"
#include "XZValueFactor.h"

// Depends on the signal.
// 1.0 --> pi radian / sample at the original sampling rate.
#define CYLINDRICAL_WAVE_PROCESSOR_UPSAMP_FILTER_HALF_TRANSITION_WIDTH 0.2

#define CYLINDRICAL_WAVE_PROCESSOR_USE_CROSS_CORRELATION 1



namespace Lab {
//###############################################################################################
//TODO: When using cross-correlation, don't calculate analytic signal / don't use complex signal.
//###############################################################################################

// x = 0 is at the center of the element group.
// y = 0
// z = 0 is at the surface of the array.
template<typename TFloat>
class CylindricalWaveProcessor {
public:
	CylindricalWaveProcessor(
			const STAConfiguration<TFloat>& config,
			STAAcquisition<TFloat>& acquisition,
			unsigned int upsamplingFactor,
			AnalyticSignalCoherenceFactorProcessor<TFloat>& coherenceFactor,
#ifdef CYLINDRICAL_WAVE_PROCESSOR_USE_CROSS_CORRELATION
			const std::vector<TFloat> &refPulse
#else
			TFloat peakOffset
#endif
			);
	~CylindricalWaveProcessor() = default;

	void process(
		unsigned int baseElement,
		unsigned int txElem,
		unsigned int firstRxElem,
		unsigned int lastRxElem,
		Matrix<XZValueFactor<TFloat>>& gridData);

private:
	class ProcessColumn;
	struct ThreadData {
		AnalyticSignalCoherenceFactorProcessor<TFloat> coherenceFactor;
		std::vector<std::complex<TFloat>> rxSignalSumList;
		std::vector<TFloat> delayList;
	};

	CylindricalWaveProcessor(const CylindricalWaveProcessor&) = delete;
	CylindricalWaveProcessor& operator=(const CylindricalWaveProcessor&) = delete;
	CylindricalWaveProcessor(CylindricalWaveProcessor&&) = delete;
	CylindricalWaveProcessor& operator=(CylindricalWaveProcessor&&) = delete;

	const STAConfiguration<TFloat>& config_;
	unsigned int deadZoneSamplesUp_;
	STAAcquisition<TFloat>& acquisition_;
	unsigned int upsamplingFactor_;
	AnalyticSignalCoherenceFactorProcessor<TFloat>& coherenceFactor_;
	unsigned int signalLength_;
	Matrix<std::complex<TFloat>> analyticSignalList_;
	typename STAAcquisition<TFloat>::AcquisitionDataType acqData_;
#ifdef CYLINDRICAL_WAVE_PROCESSOR_USE_CROSS_CORRELATION
	std::vector<TFloat> tempSignal1_;
#endif
	std::vector<TFloat> tempSignal2_;
	TFloat signalOffset_;
	Interpolator<TFloat> interpolator_;
	FFTWFilter<TFloat> revRefPulseFilter_;
	HilbertEnvelope<TFloat> envelope_;
};



template<typename TFloat>
CylindricalWaveProcessor<TFloat>::CylindricalWaveProcessor(
			const STAConfiguration<TFloat>& config,
			STAAcquisition<TFloat>& acquisition,
			unsigned int upsamplingFactor,
			AnalyticSignalCoherenceFactorProcessor<TFloat>& coherenceFactor,
#ifdef CYLINDRICAL_WAVE_PROCESSOR_USE_CROSS_CORRELATION
			const std::vector<TFloat>& refPulse
#else
			TFloat peakOffset
#endif
			)
		: config_(config)
		, deadZoneSamplesUp_((upsamplingFactor * config.samplingFrequency) * 2.0 * config.deadZoneM / config.propagationSpeed)
		, acquisition_(acquisition)
		, upsamplingFactor_(upsamplingFactor)
		, coherenceFactor_(coherenceFactor)
		, signalLength_()
{
	// This acquisition is done to obtain the ascan length.
	acquisition_.prepare(0);
	acquisition_.execute(0, acqData_);
	const std::size_t origSignalLength = acqData_.n2();

	if (upsamplingFactor_ > 1) {
		interpolator_.prepare(upsamplingFactor_, CYLINDRICAL_WAVE_PROCESSOR_UPSAMP_FILTER_HALF_TRANSITION_WIDTH);
	}

#ifdef CYLINDRICAL_WAVE_PROCESSOR_USE_CROSS_CORRELATION
	if (refPulse.size() == 0) {
		THROW_EXCEPTION(InvalidParameterException, "The vector refPulse is empty.");
	}

	// The original signal is upsampled and convolved with the reverse of the upsampled refPulse.
	signalLength_ = (origSignalLength + refPulse.size()) * upsamplingFactor_ - 1;

	tempSignal1_.resize(origSignalLength * upsamplingFactor_);

	std::vector<TFloat> revRefPulse(refPulse.size() * upsamplingFactor_);
	if (upsamplingFactor_ > 1) {
		interpolator_.interpolate(&refPulse[0], refPulse.size(), &revRefPulse[0]);
		std::reverse(revRefPulse.begin(), revRefPulse.end());
	} else {
		std::reverse_copy(refPulse.begin(), refPulse.end(), revRefPulse.begin());
	}
	revRefPulseFilter_.setCoefficients(revRefPulse);

	signalOffset_ = static_cast<TFloat>(revRefPulse.size() - 1); // cross-correlation using convolution (revRefPulseFilter_)
#else
	signalLength_ = origSignalLength * upsamplingFactor_;

	tempSignal2_.resize(signalLength_);

	signalOffset_ = (config_.samplingFrequency * upsamplingFactor_) * peakOffset / config_.centerFrequency;
#endif
	LOG_DEBUG << "signalOffset_: " << signalOffset_ << " origSignalLength: " << origSignalLength << " signalLength_: " << signalLength_;

	if (deadZoneSamplesUp_ >= signalLength_) {
		THROW_EXCEPTION(InvalidValueException, "Invalid dead zone: deadZoneSamplesUp_ (" << deadZoneSamplesUp_ <<
							") >= signalLength (" << signalLength_ << ").");
	}
}

template<typename TFloat>
void
CylindricalWaveProcessor<TFloat>::process(
		unsigned int baseElement,
		unsigned int txElem,
		unsigned int firstRxElem,
		unsigned int lastRxElem,
		Matrix<XZValueFactor<TFloat>>& gridData)
{
	LOG_DEBUG << "BEGIN ========== CylindricalWaveProcessor::process ==========";

	Util::resetValueFactor(gridData.begin(), gridData.end());

	unsigned int numActiveRxElements = lastRxElem - firstRxElem + 1;
	analyticSignalList_.resize(numActiveRxElements, signalLength_);

	// Prepare the signal list.
	acquisition_.prepare(baseElement);
	acquisition_.execute(txElem, acqData_);
	const std::size_t samplesPerChannelLow = acqData_.n2();
	for (unsigned int rxElem = firstRxElem; rxElem <= lastRxElem; ++rxElem) {
#ifdef CYLINDRICAL_WAVE_PROCESSOR_USE_CROSS_CORRELATION
		if (upsamplingFactor_ > 1) {
			// Interpolate the signal.
			interpolator_.interpolate(&acqData_(rxElem, 0), samplesPerChannelLow, &tempSignal1_[0]);
		} else {
			auto range = acqData_.range2(rxElem);
			std::copy(range.begin(), range.end(), tempSignal1_.begin());
		}

		Util::removeDC(&tempSignal1_[0], tempSignal1_.size(), deadZoneSamplesUp_);

		// Cross-correlation using convolution.
		revRefPulseFilter_.filter(tempSignal1_, tempSignal2_);
#else
		if (upsamplingFactor_ > 1) {
			// Interpolates the signal.
			interpolator_.interpolate(&acqData_(rxElem, 0), samplesPerChannelLow, &tempSignal2_[0]);
		} else {
			auto range = acqData_.range2(rxElem);
			std::copy(range.begin(), range.end(), tempSignal2_.begin());
		}

		Util::removeDC(&tempSignal2_[0], tempSignal2_.size(), deadZoneSamplesUp_);
#endif
		// Obtain the analytic signal.
		envelope_.getAnalyticSignal(&tempSignal2_[0], tempSignal2_.size(), &analyticSignalList_(rxElem - firstRxElem, 0));
	}
	LOG_DEBUG << "PRC PREP txElem = " << txElem;

	ThreadData threadData;
	threadData.coherenceFactor = coherenceFactor_;

	tbb::enumerable_thread_specific<ThreadData> processColumnTLS(threadData);
	tbb::parallel_for(
			tbb::blocked_range<std::size_t>(0, gridData.n1()),
			ProcessColumn(
				gridData.n2(),
				config_,
				txElem,
				firstRxElem,
				lastRxElem,
				upsamplingFactor_,
				processColumnTLS,
				signalOffset_,
				analyticSignalList_,
				gridData));

	LOG_DEBUG << "END ========== CylindricalWaveProcessor::process ==========";
}



template<typename TFloat>
class CylindricalWaveProcessor<TFloat>::ProcessColumn {
public:
	ProcessColumn(
			std::size_t numRows,
			const STAConfiguration<TFloat>& config,
			unsigned int txElem,
			unsigned int firstRxElem,
			unsigned int lastRxElem,
			unsigned int upsamplingFactor,
			tbb::enumerable_thread_specific<ThreadData>& processColumnTLS,
			TFloat signalOffset,
			const Matrix<std::complex<TFloat>>& analyticSignalList,
			Matrix<XZValueFactor<TFloat>>& gridData)
				: numRows_(numRows)
				, config_(config)
				, txElem_(txElem)
				, firstRxElem_(firstRxElem)
				, lastRxElem_(lastRxElem)
				, processColumnTLS_(processColumnTLS)
				, signalOffset_(signalOffset)
				, analyticSignalList_(analyticSignalList)
				, invCT_((config_.samplingFrequency * upsamplingFactor) / config_.propagationSpeed)
				, xArray_(config_.numElements)
				, gridData_(gridData) {

		const TFloat xArrayCenter = (config_.numElements - 1) * config_.pitch / 2.0;
		for (unsigned int elem = 0; elem < config_.numElements; ++elem) {
			xArray_[elem] = elem * config_.pitch - xArrayCenter;
		}
	}

	void operator()(const tbb::blocked_range<std::size_t>& r) const {
		//LOG_DEBUG << "PRC col = " << r.begin() << " n = " << (r.end() - r.begin());

		ThreadData& data = processColumnTLS_.local();

		data.rxSignalSumList.resize(config_.numElements);
		data.delayList.resize(config_.numElements);

		// For each point in x-z:
		for (std::size_t i = r.begin(); i != r.end(); ++i) {
			for (unsigned int j = 0; j < numRows_; ++j) {

				std::fill(data.rxSignalSumList.begin(), data.rxSignalSumList.end(), std::complex<TFloat>(0));
				XZValueFactor<TFloat>& point = gridData_(i, j);

				// Calculate the delays.
				for (unsigned int elem = 0; elem < config_.numElements; ++elem) {
					const TFloat dx = point.x - xArray_[elem];
					const TFloat dz = point.z /* - zArray*/; // zArray = 0
					data.delayList[elem] = std::sqrt(dx * dx + dz * dz) * invCT_;
				}

				const TFloat txDelay = data.delayList[txElem_];
				for (unsigned int rxElem = firstRxElem_; rxElem <= lastRxElem_; ++rxElem) {
					// Linear interpolation.
					const TFloat delay = signalOffset_ + txDelay + data.delayList[rxElem];
					const std::size_t delayIdx = static_cast<std::size_t>(delay);
					const TFloat k = delay - delayIdx;
					if (delayIdx < analyticSignalList_.n2() - 1) {
						const std::complex<TFloat>* p = &analyticSignalList_(rxElem - firstRxElem_, delayIdx);
						data.rxSignalSumList[rxElem] = (1 - k) * *p + k * *(p + 1);
					}
				}

				if (data.coherenceFactor.enabled()) {
					point.factor = data.coherenceFactor.calculate(&data.rxSignalSumList[0], data.rxSignalSumList.size());
				}
#ifdef CYLINDRICAL_WAVE_PROCESSOR_USE_CROSS_CORRELATION
				point.value = std::accumulate(data.rxSignalSumList.begin(), data.rxSignalSumList.end(), std::complex<TFloat>(0)).real();
#else
				point.value = std::abs(std::accumulate(data.rxSignalSumList.begin(), data.rxSignalSumList.end(), std::complex<TFloat>(0)));
#endif
			}
		}
	}
private:
	const std::size_t numRows_;
	const STAConfiguration<TFloat>& config_;
	unsigned int txElem_;
	unsigned int firstRxElem_;
	unsigned int lastRxElem_;
	tbb::enumerable_thread_specific<ThreadData>& processColumnTLS_;
	const TFloat signalOffset_;
	const Matrix<std::complex<TFloat>>& analyticSignalList_;
	const TFloat invCT_;
	std::vector<TFloat> xArray_;
	Matrix<XZValueFactor<TFloat>>& gridData_;
};

} // namespace Lab

#endif // CYLINDRICALWAVEPROCESSOR_H
