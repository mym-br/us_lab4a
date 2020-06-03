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

#ifndef CCBF_PITCH_CATCH_PROCESSOR_H
#define CCBF_PITCH_CATCH_PROCESSOR_H

#include <algorithm> /* copy, reverse, reverse_copy */
#include <cmath> /* sqrt */
#include <cstddef> /* std::size_t */
#include <memory>
#include <utility> /* pair */
#include <vector>

#include <tbb/enumerable_thread_specific.h>
#include <tbb/tbb.h>

#include "DirectFFTWFilter.h"
#include "Exception.h"
#include "ExecutionTimeMeasurement.h"
#include "Interpolator.h"
#include "Log.h"
#include "STAConfiguration.h"
#include "Tensor3.h"
#include "Util.h"



// Depends on the signal.
// 1.0 --> pi radian / sample at the original sampling rate.
#define CCBF_PITCH_CATCH_PROCESSOR_UPSAMP_FILTER_HALF_TRANSITION_WIDTH 0.2

#define CCBF_PITCH_CATCH_PROCESSOR_USE_CROSS_CORRELATION 1
//#define CCBF_PITCH_CATCH_PROCESSOR_USE_ENVELOPE_OF_CROSS_CORRELATION 1

#if !defined(CCBF_PITCH_CATCH_PROCESSOR_USE_CROSS_CORRELATION) || defined(CCBF_PITCH_CATCH_PROCESSOR_USE_ENVELOPE_OF_CROSS_CORRELATION)
# include "HilbertEnvelope.h"
#endif

//#define CCBF_PITCH_CATCH_PROCESSOR_GET_MAX_VALUE 1



namespace Lab {

// Camacho, J.
// Cruza, J. F.
// Brizuela, J.
// Fritsch, C.
// Automatic dynamic depth focusing for NDT.
// IEEE Transactions on Ultrasonics, Ferroelectrics, and Frequency Control,
// vol. 61, no. 4, pp. 673-684, 2014.
// DOI: 10.1109/TUFFC.2014.2955
//
// x = 0 is at the center of the element group.
// y = 0
// z = 0 is at the surface of the array.
template<typename TFloat>
class CCBFPitchCatchProcessor {
public:
	CCBFPitchCatchProcessor(
			const STAConfiguration<TFloat>& config,
			const Tensor3<TFloat>& acqDataList,
			unsigned int upsamplingFactor,
			unsigned int txElem,
			unsigned int firstRxElem,
			unsigned int lastRxElem,
#ifdef CCBF_PITCH_CATCH_PROCESSOR_USE_CROSS_CORRELATION
			const std::vector<TFloat>& refPulse,
#else
			TFloat peakOffset,
#endif
			unsigned int signalStartOffset);

	~CCBFPitchCatchProcessor() = default;

	void process(
		unsigned int baseElementIdx,
		unsigned int baseElement,
		std::vector<std::pair<TFloat, TFloat>>& pointPositionList);

#ifdef LAB_ENABLE_EXECUTION_TIME_MEASUREMENT
	MeasurementList<double> tPartialPrepareDataML;
	MeasurementList<double> tPartialProcessML;
	void execTimeMeasReset(unsigned int n) {
		tPartialPrepareDataML.reset(EXECUTION_TIME_MEASUREMENT_ITERATIONS * n);
		tPartialProcessML.reset(    EXECUTION_TIME_MEASUREMENT_ITERATIONS * n);
	}
	void execTimeMeasShowResults(unsigned int n) {
		EXECUTION_TIME_MEASUREMENT_LOG_TIMES_X_N("tPrepareData:", tPartialPrepareDataML, n);
		EXECUTION_TIME_MEASUREMENT_LOG_TIMES_X_N("tProcess:    ", tPartialProcessML, n);
	}
#endif

private:
	struct PrepareData;
	struct PrepareDataThreadData {
		Interpolator<TFloat> interpolator;
#ifdef CCBF_PITCH_CATCH_PROCESSOR_USE_CROSS_CORRELATION
//		FFTWFilter<TFloat> revRefPulseFilter;
		DirectFFTWFilter<TFloat> revRefPulseFilter;
		std::vector<TFloat> tempSignal1;
#endif
#if !defined(CCBF_PITCH_CATCH_PROCESSOR_USE_CROSS_CORRELATION) || defined(CCBF_PITCH_CATCH_PROCESSOR_USE_ENVELOPE_OF_CROSS_CORRELATION)
		HilbertEnvelope<TFloat> envelope;
#endif
		std::vector<TFloat> signal;
	};

	CCBFPitchCatchProcessor(const CCBFPitchCatchProcessor&) = delete;
	CCBFPitchCatchProcessor& operator=(const CCBFPitchCatchProcessor&) = delete;
	CCBFPitchCatchProcessor(CCBFPitchCatchProcessor&&) = delete;
	CCBFPitchCatchProcessor& operator=(CCBFPitchCatchProcessor&&) = delete;

	const STAConfiguration<TFloat>& config_;
	const Tensor3<TFloat>& acqDataList_;
	unsigned int upsamplingFactor_;
	unsigned int txElem_;
	unsigned int firstRxElem_;
	unsigned int lastRxElem_;
	TFloat signalOffset_;
	std::vector<TFloat> distRxElem_;
	std::unique_ptr<tbb::enumerable_thread_specific<PrepareDataThreadData>> prepareDataTLS_;
};



template<typename TFloat>
CCBFPitchCatchProcessor<TFloat>::CCBFPitchCatchProcessor(
			const STAConfiguration<TFloat>& config,
			const Tensor3<TFloat>& acqDataList,
			unsigned int upsamplingFactor,
			unsigned int txElem,
			unsigned int firstRxElem,
			unsigned int lastRxElem,
#ifdef CCBF_PITCH_CATCH_PROCESSOR_USE_CROSS_CORRELATION
			const std::vector<TFloat>& refPulse,
#else
			TFloat peakOffset,
#endif
			unsigned int signalStartOffset)
		: config_(config)
		, acqDataList_(acqDataList)
		, upsamplingFactor_(upsamplingFactor)
		, txElem_(txElem)
		, firstRxElem_(firstRxElem)
		, lastRxElem_(lastRxElem)
{
	const std::size_t origSignalLength = acqDataList_.n3();

	PrepareDataThreadData prepareDataThreadData;
	if (upsamplingFactor_ > 1) {
		prepareDataThreadData.interpolator.prepare(upsamplingFactor_, CCBF_PITCH_CATCH_PROCESSOR_UPSAMP_FILTER_HALF_TRANSITION_WIDTH);
	}

#ifdef CCBF_PITCH_CATCH_PROCESSOR_USE_CROSS_CORRELATION
	if (refPulse.size() == 0) {
		THROW_EXCEPTION(InvalidParameterException, "The vector refPulse is empty.");
	}

	// The original signal is upsampled and convolved with the reverse of the upsampled refPulse.
	const unsigned int signalLength = (origSignalLength + refPulse.size()) * upsamplingFactor_ - 1;

	std::vector<TFloat> revRefPulse(refPulse.size() * upsamplingFactor_);
	if (upsamplingFactor_ > 1) {
		prepareDataThreadData.interpolator.interpolate(&refPulse[0], refPulse.size(), &revRefPulse[0]);
		std::reverse(revRefPulse.begin(), revRefPulse.end());
	} else {
		std::reverse_copy(refPulse.begin(), refPulse.end(), revRefPulse.begin());
	}

	signalOffset_ = (revRefPulse.size() - 1.0) - static_cast<TFloat>(signalStartOffset) * upsamplingFactor_; // cross-correlation using convolution (revRefPulseFilter_)

	LOG_DEBUG << "signalOffset_: " << signalOffset_ << " origSignalLength: " << origSignalLength << " signalLength: " << signalLength;
#else
	signalOffset_ = ((config_.samplingFrequency * upsamplingFactor_) * peakOffset / config_.centerFrequency) - static_cast<TFloat>(signalStartOffset) * upsamplingFactor_;

	LOG_DEBUG << "signalOffset_: " << signalOffset_ << " origSignalLength: " << origSignalLength;
#endif

	const unsigned int numActiveRxElements = lastRxElem_ - firstRxElem_ + 1;
	distRxElem_.resize(numActiveRxElements);

#ifdef CCBF_PITCH_CATCH_PROCESSOR_USE_CROSS_CORRELATION
	prepareDataThreadData.revRefPulseFilter.setCoefficients(revRefPulse);
	prepareDataThreadData.tempSignal1.resize(origSignalLength * upsamplingFactor_);
#else
	prepareDataThreadData.signal.resize(origSignalLength * upsamplingFactor_);
#endif
	prepareDataTLS_ = std::make_unique<tbb::enumerable_thread_specific<PrepareDataThreadData>>(prepareDataThreadData);
}

template<typename TFloat>
void
CCBFPitchCatchProcessor<TFloat>::process(
		unsigned int baseElementIdx,
		unsigned int baseElement,
		std::vector<std::pair<TFloat, TFloat>>& pointPositionList)
{
	//LOG_DEBUG << "BEGIN ========== CCBFPitchCatchProcessor::process ==========";

	const std::size_t samplesPerChannelLow = acqDataList_.n3();

	// Obtain the traveled distances.
#ifdef LAB_ENABLE_EXECUTION_TIME_MEASUREMENT
	Timer prepareDataTimer;
#endif
	PrepareData prepareDataOp = {
		config_,
		samplesPerChannelLow,
		acqDataList_,
		upsamplingFactor_,
		baseElementIdx,
		firstRxElem_,
		signalOffset_,
		*prepareDataTLS_,
		distRxElem_
	};
#ifdef CCBF_PITCH_CATCH_PROCESSOR_GET_MAX_VALUE
	prepareDataOp(tbb::blocked_range<unsigned int>(firstRxElem_, lastRxElem_ + 1)); // single-thread
#else
	tbb::parallel_for(tbb::blocked_range<unsigned int>(firstRxElem_, lastRxElem_ + 1), prepareDataOp);
#endif
#ifdef LAB_ENABLE_EXECUTION_TIME_MEASUREMENT
	tPartialPrepareDataML.put(prepareDataTimer.getTime());

	Timer processTimer;
#endif

	const TFloat rTx = distRxElem_[txElem_ - firstRxElem_] * 0.5f;
	const TFloat xTxElem = config_.pitch * (baseElement + txElem_);

	for (unsigned int rxElem = firstRxElem_; rxElem <= lastRxElem_; ++rxElem) {
		if (rxElem == txElem_) continue;

		const TFloat r = distRxElem_[rxElem - firstRxElem_];
		const TFloat d = ((rxElem > txElem_) ? rxElem - txElem_ : txElem_ - rxElem) * config_.pitch;
		const TFloat d2 = d * d;

		const TFloat rTx2 = rTx * rTx;
		const TFloat r1 = (4.0f * rTx2 * r) / (4.0f * rTx2 + r * r - d2);
		const TFloat r2 = r - r1;

		const TFloat r1_2 = r1 * r1;
		const TFloat g = (r1_2 - r2 * r2 + d2) / (2.0f * d);
		const TFloat h = std::sqrt(r1_2 - g * g);

		if (rxElem > txElem_) {
			pointPositionList.emplace_back(xTxElem + g, h);

		} else {
			pointPositionList.emplace_back(xTxElem - g, h);
		}
	}

#ifdef LAB_ENABLE_EXECUTION_TIME_MEASUREMENT
	tPartialProcessML.put(processTimer.getTime());
#endif
	//LOG_DEBUG << "END ========== CCBFPitchCatchProcessor::process ==========";
}



template<typename TFloat>
struct CCBFPitchCatchProcessor<TFloat>::PrepareData {
	void operator()(const tbb::blocked_range<unsigned int>& r) const {
		PrepareDataThreadData& local = prepareDataTLS.local();

#ifdef CCBF_PITCH_CATCH_PROCESSOR_GET_MAX_VALUE
		TFloat minMaxValue = 1e9;
#endif

		for (unsigned int rxElem = r.begin(); rxElem != r.end(); ++rxElem) {
#ifdef CCBF_PITCH_CATCH_PROCESSOR_USE_CROSS_CORRELATION
			if (upsamplingFactor > 1) {
				// Interpolate the signal.
				local.interpolator.interpolate(&acqDataList(baseElementIdx, rxElem - firstRxElem, 0), samplesPerChannelLow, &local.tempSignal1[0]);
			} else {
				auto range = acqDataList.range3(baseElementIdx, rxElem - firstRxElem);
				std::copy(range.begin(), range.end(), local.tempSignal1.begin());
			}

			// Cross-correlation using convolution.
			local.revRefPulseFilter.filter(local.tempSignal1, local.signal);

# ifdef CCBF_PITCH_CATCH_PROCESSOR_USE_ENVELOPE_OF_CROSS_CORRELATION
			Util::removeDC(&local.signal[0], local.signal.size());
			local.envelope.calculate(&local.signal[0], local.signal.size());
# endif
#else
			if (upsamplingFactor > 1) {
				// Interpolates the signal.
				local.interpolator.interpolate(&acqDataList(baseElementIdx, rxElem - firstRxElem, 0), samplesPerChannelLow, &local.signal[0]);
			} else {
				auto range = acqDataList.range3(baseElementIdx, rxElem - firstRxElem);
				std::copy(range.begin(), range.end(), local.signal.begin());
			}

			Util::removeDC(&local.signal[0], local.signal.size());
			local.envelope.calculate(&local.signal[0], local.signal.size());
#endif
			// Find the peak.
			TFloat maxValue = 0;
			int idxMax = 0;
			for (unsigned int i = 0; i < local.signal.size(); ++i) {
				if (local.signal[i] > maxValue) {
					maxValue = local.signal[i];
					idxMax = i;
				}
			}
#ifdef CCBF_PITCH_CATCH_PROCESSOR_GET_MAX_VALUE
			if (maxValue < minMaxValue) {
				minMaxValue = maxValue;
			}
#endif
			const TFloat tPulseEcho = (idxMax - signalOffset) / (config.samplingFrequency * upsamplingFactor);
			if (tPulseEcho <= 0) {
				THROW_EXCEPTION(InvalidValueException, "tPulseEcho <= 0");
			}
			distRxElem[rxElem - firstRxElem] = config.propagationSpeed * tPulseEcho;
			//LOG_DEBUG << "tPulseEcho: " << tPulseEcho << " distRxElem: " << distRxElem[rxElem - firstRxElem];
		}
#ifdef CCBF_PITCH_CATCH_PROCESSOR_GET_MAX_VALUE
		LOG_DEBUG << "minMaxValue: " << minMaxValue;
#endif
	}

	// Use only small types or references/pointers, because this object may be copied many times.
	const STAConfiguration<TFloat>& config;
	const std::size_t samplesPerChannelLow;
	const Tensor3<TFloat>& acqDataList;
	const unsigned int upsamplingFactor;
	const unsigned int baseElementIdx;
	const unsigned int firstRxElem;
	const TFloat signalOffset;
	tbb::enumerable_thread_specific<PrepareDataThreadData>& prepareDataTLS;
	std::vector<TFloat>& distRxElem;
};

} // namespace Lab

#endif // CCBF_PITCH_CATCH_PROCESSOR_H
