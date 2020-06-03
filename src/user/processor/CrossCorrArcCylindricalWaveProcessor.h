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

#ifndef CROSSCORRARCCYLINDRICALWAVEPROCESSOR_H
#define CROSSCORRARCCYLINDRICALWAVEPROCESSOR_H

#include <algorithm> /* copy, fill, reverse, reverse_copy */
#include <cmath> /* abs, sqrt */
#include <complex>
#include <cstddef> /* std::size_t */
#include <memory>
#include <numeric> /* accumulate */
#include <vector>

#include <tbb/enumerable_thread_specific.h>
#include <tbb/tbb.h>

#include "ArrayGeometry.h"
#include "CoherenceFactor.h"
#include "DirectFFTWFilter.h"
#include "Exception.h"
#include "ExecutionTimeMeasurement.h"
#include "HilbertEnvelope.h"
#include "Interpolator.h"
#include "Log.h"
#include "Matrix.h"
#include "STAConfiguration.h"
#include "Tensor3.h"
#include "Timer.h"
#include "Util.h"
#include "XZValueFactor.h"



// Depends on the signal.
// 1.0 --> pi radian / sample at the original sampling rate.
#define CROSS_CORR_ARC_CYLINDRICAL_WAVE_PROCESSOR_UPSAMP_FILTER_HALF_TRANSITION_WIDTH 0.2

//#define CROSS_CORR_ARC_CYLINDRICAL_WAVE_PROCESSOR_USE_ANALYTIC_SIGNAL_FOR_ANGLE 1
//#define CROSS_CORR_ARC_CYLINDRICAL_WAVE_PROCESSOR_USE_ANALYTIC_SIGNAL_FOR_RANGE 1



// Uses the peak of the cross-correlation between the signal and a reference pulse to measure time-of-flight.
namespace Lab {

// x = 0 is at the center of the element group.
// y = 0
// z = 0 is at the surface of the array.
template<typename TFloat>
class CrossCorrArcCylindricalWaveProcessor {
public:
	CrossCorrArcCylindricalWaveProcessor(
			const STAConfiguration<TFloat>& config,
			const Tensor3<TFloat>& acqDataList,
			unsigned int upsamplingFactor,
			unsigned int highUpsamplingFactor,
			unsigned int txElem,
			unsigned int firstRxElem,
			unsigned int lastRxElem,
			TFloat arcStepDiv,
			TFloat minArcAngle,
			TFloat maxArcAngle,
#ifdef CROSS_CORR_ARC_CYLINDRICAL_WAVE_PROCESSOR_USE_ANALYTIC_SIGNAL_FOR_ANGLE
			AnalyticSignalCoherenceFactorProcessor<TFloat>& coherenceFactor,
#else
			CoherenceFactorProcessor<TFloat>& coherenceFactor,
#endif
			const std::vector<TFloat>& refPulse,
			unsigned int signalStartOffset);

	~CrossCorrArcCylindricalWaveProcessor() = default;

	void process(unsigned int baseElementIdx, std::vector<XZValueFactor<TFloat>>& gridData);

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
		DirectFFTWFilter<TFloat> revRefPulseFilter;
//		FFTWFilter<TFloat> revRefPulseFilter;
#ifdef CROSS_CORR_ARC_CYLINDRICAL_WAVE_PROCESSOR_USE_ANALYTIC_SIGNAL_FOR_ANGLE
		HilbertEnvelope<TFloat> envelope;
#endif
		std::vector<TFloat> tempSignal1;
		std::vector<TFloat> tempSignal2;
	};

	struct ProcessArcRange;
	struct ProcessArcRangeThreadData {
#ifdef CROSS_CORR_ARC_CYLINDRICAL_WAVE_PROCESSOR_USE_ANALYTIC_SIGNAL_FOR_ANGLE
		AnalyticSignalCoherenceFactorProcessor<TFloat> coherenceFactor;
		std::vector<std::complex<TFloat>> rxSignalSumList;
#else
		CoherenceFactorProcessor<TFloat> coherenceFactor;
		std::vector<TFloat> rxSignalSumList;
#endif
		std::vector<TFloat> delayList;
	};

	CrossCorrArcCylindricalWaveProcessor(const CrossCorrArcCylindricalWaveProcessor&) = delete;
	CrossCorrArcCylindricalWaveProcessor& operator=(const CrossCorrArcCylindricalWaveProcessor&) = delete;
	CrossCorrArcCylindricalWaveProcessor(CrossCorrArcCylindricalWaveProcessor&&) = delete;
	CrossCorrArcCylindricalWaveProcessor& operator=(CrossCorrArcCylindricalWaveProcessor&&) = delete;

	const STAConfiguration<TFloat>& config_;
	const Tensor3<TFloat>& acqDataList_;
	unsigned int upsamplingFactor_;
	unsigned int highUpsamplingFactor_;
	unsigned int txElem_;
	unsigned int firstRxElem_;
	unsigned int lastRxElem_;
	TFloat arcStepDiv_;
	TFloat minArcAngle_;
	TFloat maxArcAngle_;
#ifdef CROSS_CORR_ARC_CYLINDRICAL_WAVE_PROCESSOR_USE_ANALYTIC_SIGNAL_FOR_ANGLE
	AnalyticSignalCoherenceFactorProcessor<TFloat>& coherenceFactor_;
#else
	CoherenceFactorProcessor<TFloat>& coherenceFactor_;
#endif
	unsigned int signalLength_;
	unsigned int centerSignalLength_;
#ifdef CROSS_CORR_ARC_CYLINDRICAL_WAVE_PROCESSOR_USE_ANALYTIC_SIGNAL_FOR_ANGLE
	Matrix<std::complex<TFloat>> analyticSignalList_;
#else
	Matrix<TFloat> signalList_;
#endif
#ifdef CROSS_CORR_ARC_CYLINDRICAL_WAVE_PROCESSOR_USE_ANALYTIC_SIGNAL_FOR_RANGE
	std::vector<std::complex<TFloat>> centerAnalyticSignal_;
#else
	std::vector<TFloat> centerSignal_;
#endif
	std::vector<TFloat> centerTempSignal1_;
	std::vector<TFloat> centerTempSignal2_;
	TFloat signalOffset_;
	TFloat centerSignalOffset_;
	Interpolator<TFloat> centerInterpolator_;
	DirectFFTWFilter<TFloat> centerRevRefPulseFilter_;
//	FFTWFilter<TFloat> centerRevRefPulseFilter_;
	std::vector<TFloat> arcAngleList_;
#ifdef CROSS_CORR_ARC_CYLINDRICAL_WAVE_PROCESSOR_USE_ANALYTIC_SIGNAL_FOR_RANGE
	HilbertEnvelope<TFloat> envelope_;
#endif
	std::vector<TFloat> xArray_;
	std::unique_ptr<tbb::enumerable_thread_specific<PrepareDataThreadData>> prepareDataTLS_;
	std::unique_ptr<tbb::enumerable_thread_specific<ProcessArcRangeThreadData>> processArcRangeTLS_;
};



template<typename TFloat>
CrossCorrArcCylindricalWaveProcessor<TFloat>::CrossCorrArcCylindricalWaveProcessor(
			const STAConfiguration<TFloat>& config,
			const Tensor3<TFloat>& acqDataList,
			unsigned int upsamplingFactor,
			unsigned int highUpsamplingFactor,
			unsigned int txElem,
			unsigned int firstRxElem,
			unsigned int lastRxElem,
			TFloat arcStepDiv,
			TFloat minArcAngle,
			TFloat maxArcAngle,
#ifdef CROSS_CORR_ARC_CYLINDRICAL_WAVE_PROCESSOR_USE_ANALYTIC_SIGNAL_FOR_ANGLE
			AnalyticSignalCoherenceFactorProcessor<TFloat>& coherenceFactor,
#else
			CoherenceFactorProcessor<TFloat>& coherenceFactor,
#endif
			const std::vector<TFloat>& refPulse,
			unsigned int signalStartOffset)
		: config_(config)
		, acqDataList_(acqDataList)
		, upsamplingFactor_(upsamplingFactor)
		, highUpsamplingFactor_(highUpsamplingFactor)
		, txElem_(txElem)
		, firstRxElem_(firstRxElem)
		, lastRxElem_(lastRxElem)
		, arcStepDiv_(arcStepDiv)
		, minArcAngle_(minArcAngle)
		, maxArcAngle_(maxArcAngle)
		, coherenceFactor_(coherenceFactor)
		, signalLength_()
		, centerSignalLength_()
		, signalOffset_()
		, centerSignalOffset_()
{
	const std::size_t origSignalLength = acqDataList_.n3();

	if (highUpsamplingFactor_ > 1) {
		centerInterpolator_.prepare(highUpsamplingFactor_, CROSS_CORR_ARC_CYLINDRICAL_WAVE_PROCESSOR_UPSAMP_FILTER_HALF_TRANSITION_WIDTH);
	}

	if (refPulse.size() == 0) {
		THROW_EXCEPTION(InvalidParameterException, "The vector refPulse is empty.");
	}

	// The original signal is upsampled and convolved with the reverse of the upsampled refPulse.
	signalLength_ = (origSignalLength + refPulse.size()) * upsamplingFactor_ - 1;
	centerSignalLength_ = (origSignalLength + refPulse.size()) * highUpsamplingFactor_ - 1;

	PrepareDataThreadData prepareDataThreadData;
	if (upsamplingFactor_ > 1) {
		prepareDataThreadData.interpolator.prepare(upsamplingFactor_, CROSS_CORR_ARC_CYLINDRICAL_WAVE_PROCESSOR_UPSAMP_FILTER_HALF_TRANSITION_WIDTH);
	}

	prepareDataThreadData.tempSignal1.resize(origSignalLength * upsamplingFactor_);
	centerTempSignal1_.resize(origSignalLength * highUpsamplingFactor_);

	std::vector<TFloat> revRefPulse(refPulse.size() * upsamplingFactor_);
	if (upsamplingFactor_ > 1) {
		prepareDataThreadData.interpolator.interpolate(&refPulse[0], refPulse.size(), &revRefPulse[0]);
		std::reverse(revRefPulse.begin(), revRefPulse.end());
	} else {
		std::reverse_copy(refPulse.begin(), refPulse.end(), revRefPulse.begin());
	}

	std::vector<TFloat> centerRevRefPulse(refPulse.size() * highUpsamplingFactor_);
	if (highUpsamplingFactor_ > 1) {
		centerInterpolator_.interpolate(&refPulse[0], refPulse.size(), &centerRevRefPulse[0]);
		std::reverse(centerRevRefPulse.begin(), centerRevRefPulse.end());
	} else {
		std::reverse_copy(refPulse.begin(), refPulse.end(), centerRevRefPulse.begin());
	}
	centerRevRefPulseFilter_.setCoefficients(centerRevRefPulse);

	signalOffset_ = (revRefPulse.size() - 1.0) - static_cast<TFloat>(signalStartOffset) * upsamplingFactor_; // cross-correlation using convolution (revRefPulseFilter_)
	centerSignalOffset_ = (centerRevRefPulse.size() - 1.0) - static_cast<TFloat>(signalStartOffset) * highUpsamplingFactor_;

	const unsigned int numActiveRxElements = lastRxElem_ - firstRxElem_ + 1;
#ifdef CROSS_CORR_ARC_CYLINDRICAL_WAVE_PROCESSOR_USE_ANALYTIC_SIGNAL_FOR_ANGLE
	analyticSignalList_.resize(numActiveRxElements, signalLength_);
#else
	signalList_.resize(numActiveRxElements, signalLength_);
#endif
#ifdef CROSS_CORR_ARC_CYLINDRICAL_WAVE_PROCESSOR_USE_ANALYTIC_SIGNAL_FOR_RANGE
	centerAnalyticSignal_.resize(centerSignalLength_);
#else
	centerSignal_.resize(centerSignalLength_);
#endif

	LOG_DEBUG << "signalOffset_: " << signalOffset_ << " centerSignalOffset_: " << centerSignalOffset_ << " origSignalLength: " << origSignalLength << '\n' <<
			" signalLength_: " << signalLength_ << " centerSignalLength_: " << centerSignalLength_;

	ArrayGeometry::getElementXCentered2D(config_.numElements, config_.pitch, xArray_);
	LOG_DEBUG << "xArray_[txElem_] = " << xArray_[txElem_];

	prepareDataThreadData.revRefPulseFilter.setCoefficients(revRefPulse);
	prepareDataThreadData.tempSignal1.resize(origSignalLength * upsamplingFactor_);
	prepareDataTLS_ = std::make_unique<tbb::enumerable_thread_specific<PrepareDataThreadData>>(prepareDataThreadData);

	ProcessArcRangeThreadData processArcRangeThreadData;
	processArcRangeThreadData.coherenceFactor = coherenceFactor_;
	processArcRangeTLS_ = std::make_unique<tbb::enumerable_thread_specific<ProcessArcRangeThreadData>>(processArcRangeThreadData);

	//TODO: arcAngleList.reserve / arcData.reserve
}

template<typename TFloat>
void
CrossCorrArcCylindricalWaveProcessor<TFloat>::process(unsigned int baseElementIdx, std::vector<XZValueFactor<TFloat>>& arcData)
{
	//LOG_DEBUG << "BEGIN ========== CrossCorrArcCylindricalWaveProcessor::process ==========";

	Util::resetValueFactor(arcData.begin(), arcData.end());
	const std::size_t samplesPerChannelLow = acqDataList_.n3();

	// Prepare the signal list.
#ifdef LAB_ENABLE_EXECUTION_TIME_MEASUREMENT
	Timer prepareDataTimer;
#endif
	PrepareData prepareDataOp = {
		samplesPerChannelLow,
		acqDataList_,
		upsamplingFactor_,
		baseElementIdx,
		firstRxElem_,
		*prepareDataTLS_,
#ifdef CROSS_CORR_ARC_CYLINDRICAL_WAVE_PROCESSOR_USE_ANALYTIC_SIGNAL_FOR_ANGLE
		analyticSignalList_
#else
		signalList_
#endif
	};
	tbb::parallel_for(tbb::blocked_range<unsigned int>(firstRxElem_, lastRxElem_ + 1), prepareDataOp);
	//prepareDataOp(tbb::blocked_range<unsigned int>(firstRxElem_, lastRxElem_ + 1)); // single-thread
#ifdef LAB_ENABLE_EXECUTION_TIME_MEASUREMENT
	tPartialPrepareDataML.put(prepareDataTimer.getTime());
#endif

#ifdef CROSS_CORR_ARC_CYLINDRICAL_WAVE_PROCESSOR_USE_ANALYTIC_SIGNAL_FOR_RANGE
# ifdef CROSS_CORR_ARC_CYLINDRICAL_WAVE_PROCESSOR_USE_ANALYTIC_SIGNAL_FOR_ANGLE
	if (highUpsamplingFactor_ != upsamplingFactor_) {
# else
	{
# endif
#else
# ifndef CROSS_CORR_ARC_CYLINDRICAL_WAVE_PROCESSOR_USE_ANALYTIC_SIGNAL_FOR_ANGLE
	if (highUpsamplingFactor_ != upsamplingFactor_) {
# else
	{
# endif
#endif
		// Prepare the center signal.

		Timer t2;

		const unsigned int rxElem = txElem_;
		if (highUpsamplingFactor_ > 1) {
			// Interpolate the signal.
			centerInterpolator_.interpolate(&acqDataList_(baseElementIdx, rxElem - firstRxElem_, 0), samplesPerChannelLow, &centerTempSignal1_[0]);
		} else {
			auto range = acqDataList_.range3(baseElementIdx, rxElem - firstRxElem_);
			std::copy(range.begin(), range.end(), centerTempSignal1_.begin());
		}

		// Cross-correlation using convolution.
		centerRevRefPulseFilter_.filter(centerTempSignal1_, centerTempSignal2_);

#ifdef CROSS_CORR_ARC_CYLINDRICAL_WAVE_PROCESSOR_USE_ANALYTIC_SIGNAL_FOR_RANGE
		Util::removeDC(&centerTempSignal2_[0], centerTempSignal2_.size());
		envelope_.getAnalyticSignal(&centerTempSignal2_[0], centerTempSignal2_.size(), &centerAnalyticSignal_[0]);
#else
		std::copy(centerTempSignal2_.begin(), centerTempSignal2_.end(), centerSignal_.begin());
#endif
		LOG_DEBUG << "PREP 2 time = " << t2.getTime();

#ifdef CROSS_CORR_ARC_CYLINDRICAL_WAVE_PROCESSOR_USE_ANALYTIC_SIGNAL_FOR_RANGE
# ifdef CROSS_CORR_ARC_CYLINDRICAL_WAVE_PROCESSOR_USE_ANALYTIC_SIGNAL_FOR_ANGLE
	} else {
		auto srcRange = analyticSignalList_.range2(txElem_ - firstRxElem_);
		std::copy(srcRange.begin(), srcRange.end(), centerAnalyticSignal_.begin());
	}
# else
	}
# endif
#else
# ifndef CROSS_CORR_ARC_CYLINDRICAL_WAVE_PROCESSOR_USE_ANALYTIC_SIGNAL_FOR_ANGLE
	} else {
		auto srcRange = signalList_.range2(txElem_ - firstRxElem_);
		std::copy(srcRange.begin(), srcRange.end(), centerSignal_.begin());
	}
# else
	}
# endif
#endif

#ifdef LAB_ENABLE_EXECUTION_TIME_MEASUREMENT
	Timer processTimer;
#endif
	// Find the radius of the center element arc using the center signal.
	TFloat maxValue = 0;
	int idxMax = 0;
#ifdef CROSS_CORR_ARC_CYLINDRICAL_WAVE_PROCESSOR_USE_ANALYTIC_SIGNAL_FOR_RANGE
	for (unsigned int i = 0; i < centerAnalyticSignal_.size(); ++i) {
		const TFloat absValue = std::abs(centerAnalyticSignal_[i]);
		if (absValue > maxValue) {
			maxValue = absValue;
			idxMax = i;
		}
	}
#else
	for (unsigned int i = 0; i < centerSignal_.size(); ++i) {
		if (centerSignal_[i] > maxValue) {
			maxValue = centerSignal_[i];
			idxMax = i;
		}
	}
#endif

	const TFloat tPulseEcho = (idxMax - centerSignalOffset_) / (config_.samplingFrequency * highUpsamplingFactor_);
	if (tPulseEcho <= 0) {
		THROW_EXCEPTION(InvalidValueException, "tPulseEcho <= 0");
	}
	const TFloat radius = (config_.propagationSpeed * tPulseEcho) / 2;
	LOG_DEBUG << "[RADIUS] baseElementIdx: " << baseElementIdx << " tPulseEcho: " << tPulseEcho << " radius: " << radius;

	const TFloat lambda = config_.propagationSpeed / config_.centerFrequency;
	Util::fillSequenceFromStartToEndWithMaximumStep(
				arcAngleList_,
				Util::degreeToRadian(minArcAngle_),
				Util::degreeToRadian(maxArcAngle_),
				(lambda / arcStepDiv_) / radius);

	// Fill the point coordinates in the arc.
	arcData.resize(arcAngleList_.size());
	const TFloat xArcCenter = (-(config_.numElements - 1.0) / 2.0 + txElem_) * config_.pitch;
	//LOG_DEBUG << "xArcCenter = " << xArcCenter;
	for (unsigned int i = 0; i < arcAngleList_.size(); ++i) {
		const TFloat angle = arcAngleList_[i];
		arcData[i].x = radius * std::cos(angle) + xArcCenter;
		arcData[i].z = radius * std::sin(angle);
		//LOG_DEBUG << "angle = " << angle << " x = " << arcData[i].x << " z = " << arcData[i].z;
	}

	ProcessArcRange processArcRangeOp = {
		config_,
		txElem_,
		firstRxElem_,
		lastRxElem_,
		signalOffset_,
#ifdef CROSS_CORR_ARC_CYLINDRICAL_WAVE_PROCESSOR_USE_ANALYTIC_SIGNAL_FOR_ANGLE
		analyticSignalList_,
#else
		signalList_,
#endif
		(config_.samplingFrequency * upsamplingFactor_) / config_.propagationSpeed,
		xArray_,
		*processArcRangeTLS_,
		arcData
	};
	tbb::parallel_for(tbb::blocked_range<std::size_t>(0, arcData.size()), processArcRangeOp);

#ifdef LAB_ENABLE_EXECUTION_TIME_MEASUREMENT
	tPartialProcessML.put(processTimer.getTime());
#endif

	//LOG_DEBUG << "END ========== CrossCorrArcCylindricalWaveProcessor::process ==========";
}



template<typename TFloat>
struct CrossCorrArcCylindricalWaveProcessor<TFloat>::PrepareData {
	void operator()(const tbb::blocked_range<unsigned int>& r) const {
		PrepareDataThreadData& local = prepareDataTLS.local();

		for (unsigned int rxElem = r.begin(); rxElem != r.end(); ++rxElem) {
			if (upsamplingFactor > 1) {
				// Interpolate the signal.
				local.interpolator.interpolate(&acqDataList(baseElementIdx, rxElem - firstRxElem, 0), samplesPerChannelLow, &local.tempSignal1[0]);
			} else {
				auto range = acqDataList.range3(baseElementIdx, rxElem - firstRxElem);
				std::copy(range.begin(), range.end(), local.tempSignal1.begin());
			}

			// Cross-correlation using convolution.
			local.revRefPulseFilter.filter(local.tempSignal1, local.tempSignal2);

			Util::removeDC(&local.tempSignal2[0], local.tempSignal2.size());

#ifdef CROSS_CORR_ARC_CYLINDRICAL_WAVE_PROCESSOR_USE_ANALYTIC_SIGNAL_FOR_ANGLE
			local.envelope.getAnalyticSignal(&local.tempSignal2[0], local.tempSignal2.size(), &analyticSignalList(rxElem - firstRxElem, 0));
#else
			std::copy(local.tempSignal2.begin(), local.tempSignal2.end(), signalList.range2(rxElem - firstRxElem).begin());
#endif
		}
	}

	// Use only small types or references/pointers, because this object may be copied many times.
	const std::size_t samplesPerChannelLow;
	const Tensor3<TFloat>& acqDataList;
	const unsigned int upsamplingFactor;
	const unsigned int baseElementIdx;
	const unsigned int firstRxElem;
	tbb::enumerable_thread_specific<PrepareDataThreadData>& prepareDataTLS;
#ifdef CROSS_CORR_ARC_CYLINDRICAL_WAVE_PROCESSOR_USE_ANALYTIC_SIGNAL_FOR_ANGLE
	Matrix<std::complex<TFloat>>& analyticSignalList;
#else
	Matrix<TFloat>& signalList;
#endif
};



template<typename TFloat>
struct CrossCorrArcCylindricalWaveProcessor<TFloat>::ProcessArcRange {
	void operator()(const tbb::blocked_range<std::size_t>& r) const {
		//LOG_DEBUG << "PRC col = " << r.begin() << " n = " << (r.end() - r.begin());

		ProcessArcRangeThreadData& data = processArcRangeTLS.local();

		data.rxSignalSumList.resize(config.numElements);
		const unsigned int numRxElem = lastRxElem - firstRxElem + 1;
		data.delayList.resize(numRxElem);

		// For each point in the arc:
		for (std::size_t i = r.begin(); i != r.end(); ++i) {
#ifdef CROSS_CORR_ARC_CYLINDRICAL_WAVE_PROCESSOR_USE_ANALYTIC_SIGNAL_FOR_ANGLE
			std::fill(data.rxSignalSumList.begin(), data.rxSignalSumList.end(), std::complex<TFloat>(0));
#else
			std::fill(data.rxSignalSumList.begin(), data.rxSignalSumList.end(), TFloat(0));
#endif
			XZValueFactor<TFloat>& point = arcData[i];

			// Calculate the delays.
			for (unsigned int elem = firstRxElem; elem <= lastRxElem; ++elem) {
				const TFloat dx = point.x - xArray[elem];
				const TFloat dz = point.z /* - zArray*/; // zArray = 0
				data.delayList[elem - firstRxElem] = std::sqrt(dx * dx + dz * dz) * invCT;
			}

			const TFloat txDelay = data.delayList[txElem - firstRxElem];
			for (unsigned int rxElem = firstRxElem; rxElem <= lastRxElem; ++rxElem) {
				// Linear interpolation.
				const TFloat delay = signalOffset + (txDelay + data.delayList[rxElem - firstRxElem]);
				//LOG_DEBUG << "delay = " << delay;
				if (delay >= 0) {
					const unsigned int delayIdx = static_cast<unsigned int>(delay);
#ifdef CROSS_CORR_ARC_CYLINDRICAL_WAVE_PROCESSOR_USE_ANALYTIC_SIGNAL_FOR_ANGLE
					if (delayIdx < analyticSignalList.n2() - 1) {
						const TFloat k = delay - delayIdx;
						const std::complex<TFloat>* p = &analyticSignalList(rxElem - firstRxElem, delayIdx);
#else
					if (delayIdx < signalList.n2() - 1) {
						const TFloat k = delay - delayIdx;
						const TFloat* p = &signalList(rxElem - firstRxElem, delayIdx);
#endif
						data.rxSignalSumList[rxElem] = (1 - k) * *p + k * *(p + 1);
					} else {
						data.rxSignalSumList[rxElem] = 0;
					}
				}
			}

			if (data.coherenceFactor.enabled()) {
				point.factor = data.coherenceFactor.calculate(&data.rxSignalSumList[0], data.rxSignalSumList.size());
			}
#ifdef CROSS_CORR_ARC_CYLINDRICAL_WAVE_PROCESSOR_USE_ANALYTIC_SIGNAL_FOR_ANGLE
			point.value = std::abs(std::accumulate(data.rxSignalSumList.begin(), data.rxSignalSumList.end(), std::complex<TFloat>(0)));
#else
			point.value = std::accumulate(data.rxSignalSumList.begin(), data.rxSignalSumList.end(), TFloat(0));
#endif
		}
	}

	// Use only small types or references/pointers, because this object may be copied many times.
	const STAConfiguration<TFloat>& config;
	const unsigned int txElem;
	const unsigned int firstRxElem;
	const unsigned int lastRxElem;
	const TFloat signalOffset;
#ifdef CROSS_CORR_ARC_CYLINDRICAL_WAVE_PROCESSOR_USE_ANALYTIC_SIGNAL_FOR_ANGLE
	const Matrix<std::complex<TFloat>>& analyticSignalList;
#else
	const Matrix<TFloat>& signalList;
#endif
	const TFloat invCT;
	const std::vector<TFloat>& xArray;
	tbb::enumerable_thread_specific<ProcessArcRangeThreadData>& processArcRangeTLS;
	std::vector<XZValueFactor<TFloat>>& arcData;
};

} // namespace Lab

#endif // CROSSCORRARCCYLINDRICALWAVEPROCESSOR_H
