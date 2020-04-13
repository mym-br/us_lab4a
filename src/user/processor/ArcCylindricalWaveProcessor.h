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

#ifndef ARC_CYLINDRICAL_WAVE_PROCESSOR_H
#define ARC_CYLINDRICAL_WAVE_PROCESSOR_H

#include <algorithm> /* copy, fill */
#include <cmath> /* sqrt */
#include <complex>
#include <cstddef> /* std::size_t */
#include <memory>
#include <numeric> /* accumulate */
#include <vector>

#include <tbb/enumerable_thread_specific.h>
#include <tbb/tbb.h>

#include "ArrayGeometry.h"
#include "CoherenceFactor.h"
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
#define ARC_CYLINDRICAL_WAVE_PROCESSOR_UPSAMP_FILTER_HALF_TRANSITION_WIDTH 0.2

//#define ARC_CYLINDRICAL_WAVE_PROCESSOR_USE_ANALYTIC_SIGNAL_FOR_ANGLE 1



// Uses the peak of the signal envelope to measure time-of-flight.
namespace Lab {

// x = 0 is at the center of the element group.
// y = 0
// z = 0 is at the surface of the array.
template<typename TFloat>
class ArcCylindricalWaveProcessor {
public:
	ArcCylindricalWaveProcessor(
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
#ifdef ARC_CYLINDRICAL_WAVE_PROCESSOR_USE_ANALYTIC_SIGNAL_FOR_ANGLE
			AnalyticSignalCoherenceFactorProcessor<TFloat>& coherenceFactor,
#else
			CoherenceFactorProcessor<TFloat>& coherenceFactor,
#endif
			TFloat peakOffset,
			unsigned int signalStartOffset);

	~ArcCylindricalWaveProcessor() = default;

	void process(unsigned int baseElementIdx, std::vector<XZValueFactor<TFloat>>& gridData);

#ifdef USE_EXECUTION_TIME_MEASUREMENT
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
#ifdef ARC_CYLINDRICAL_WAVE_PROCESSOR_USE_ANALYTIC_SIGNAL_FOR_ANGLE
		HilbertEnvelope<TFloat> envelope;
#endif
		std::vector<TFloat> tempSignal1;
	};

	struct ProcessArcRange;
	struct ProcessArcRangeThreadData {
#ifdef ARC_CYLINDRICAL_WAVE_PROCESSOR_USE_ANALYTIC_SIGNAL_FOR_ANGLE
		AnalyticSignalCoherenceFactorProcessor<TFloat> coherenceFactor;
		std::vector<std::complex<TFloat>> rxSignalSumList;
#else
		CoherenceFactorProcessor<TFloat> coherenceFactor;
		std::vector<TFloat> rxSignalSumList;
#endif
		std::vector<TFloat> delayList;
	};

	ArcCylindricalWaveProcessor(const ArcCylindricalWaveProcessor&) = delete;
	ArcCylindricalWaveProcessor& operator=(const ArcCylindricalWaveProcessor&) = delete;
	ArcCylindricalWaveProcessor(ArcCylindricalWaveProcessor&&) = delete;
	ArcCylindricalWaveProcessor& operator=(ArcCylindricalWaveProcessor&&) = delete;

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
#ifdef ARC_CYLINDRICAL_WAVE_PROCESSOR_USE_ANALYTIC_SIGNAL_FOR_ANGLE
	AnalyticSignalCoherenceFactorProcessor<TFloat>& coherenceFactor_;
#else
	CoherenceFactorProcessor<TFloat>& coherenceFactor_;
#endif
	unsigned int signalLength_;
	unsigned int centerSignalLength_;
#ifdef ARC_CYLINDRICAL_WAVE_PROCESSOR_USE_ANALYTIC_SIGNAL_FOR_ANGLE
	Matrix<std::complex<TFloat>> analyticSignalList_;
#else
	Matrix<TFloat> signalList_;
#endif
	std::vector<TFloat> centerSignal_;
	TFloat signalOffset_;
	TFloat centerSignalOffset_;
	Interpolator<TFloat> centerInterpolator_;
	std::vector<TFloat> arcAngleList_;
	HilbertEnvelope<TFloat> envelope_;
	std::vector<TFloat> xArray_;
	std::unique_ptr<tbb::enumerable_thread_specific<PrepareDataThreadData>> prepareDataTLS_;
	std::unique_ptr<tbb::enumerable_thread_specific<ProcessArcRangeThreadData>> processArcRangeTLS_;
};



template<typename TFloat>
ArcCylindricalWaveProcessor<TFloat>::ArcCylindricalWaveProcessor(
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
#ifdef ARC_CYLINDRICAL_WAVE_PROCESSOR_USE_ANALYTIC_SIGNAL_FOR_ANGLE
			AnalyticSignalCoherenceFactorProcessor<TFloat>& coherenceFactor,
#else
			CoherenceFactorProcessor<TFloat>& coherenceFactor,
#endif
			TFloat peakOffset,
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
		centerInterpolator_.prepare(highUpsamplingFactor_, ARC_CYLINDRICAL_WAVE_PROCESSOR_UPSAMP_FILTER_HALF_TRANSITION_WIDTH);
	}

	signalLength_ = origSignalLength * upsamplingFactor_;
	centerSignalLength_ = origSignalLength * highUpsamplingFactor_;

	PrepareDataThreadData prepareDataThreadData;
	if (upsamplingFactor_ > 1) {
		prepareDataThreadData.interpolator.prepare(upsamplingFactor_, ARC_CYLINDRICAL_WAVE_PROCESSOR_UPSAMP_FILTER_HALF_TRANSITION_WIDTH);
	}

	prepareDataThreadData.tempSignal1.resize(origSignalLength * upsamplingFactor_);

	signalOffset_ = ((config.samplingFrequency * upsamplingFactor_) * peakOffset / config.centerFrequency) - static_cast<TFloat>(signalStartOffset) * upsamplingFactor_;
	centerSignalOffset_ = ((config.samplingFrequency * highUpsamplingFactor_) * peakOffset / config.centerFrequency) - static_cast<TFloat>(signalStartOffset) * highUpsamplingFactor_;

	const unsigned int numActiveRxElements = lastRxElem_ - firstRxElem_ + 1;
#ifdef ARC_CYLINDRICAL_WAVE_PROCESSOR_USE_ANALYTIC_SIGNAL_FOR_ANGLE
	analyticSignalList_.resize(numActiveRxElements, signalLength_);
#else
	signalList_.resize(numActiveRxElements, signalLength_);
#endif
	centerSignal_.resize(centerSignalLength_);

	LOG_DEBUG << "signalOffset_: " << signalOffset_ << " centerSignalOffset_: " << centerSignalOffset_ << " origSignalLength: " << origSignalLength << '\n' <<
			" signalLength_: " << signalLength_ << " centerSignalLength_: " << centerSignalLength_;

	ArrayGeometry::getElementXCentered2D(config_.numElements, config_.pitch, xArray_);
	LOG_DEBUG << "xArray_[txElem_] = " << xArray_[txElem_];

	prepareDataThreadData.tempSignal1.resize(origSignalLength * upsamplingFactor_);
	prepareDataTLS_ = std::make_unique<tbb::enumerable_thread_specific<PrepareDataThreadData>>(prepareDataThreadData);

	ProcessArcRangeThreadData processArcRangeThreadData;
	processArcRangeThreadData.coherenceFactor = coherenceFactor_;
	processArcRangeTLS_ = std::make_unique<tbb::enumerable_thread_specific<ProcessArcRangeThreadData>>(processArcRangeThreadData);

	//TODO: arcAngleList.reserve / arcData.reserve
}

template<typename TFloat>
void
ArcCylindricalWaveProcessor<TFloat>::process(unsigned int baseElementIdx, std::vector<XZValueFactor<TFloat>>& arcData)
{
	//LOG_DEBUG << "BEGIN ========== ArcCylindricalWaveProcessor::process ==========";

	Util::resetValueFactor(arcData.begin(), arcData.end());
	const std::size_t samplesPerChannelLow = acqDataList_.n3();

	// Prepare the signal list.
#ifdef USE_EXECUTION_TIME_MEASUREMENT
	Timer prepareDataTimer;
#endif
	PrepareData prepareDataOp = {
		samplesPerChannelLow,
		acqDataList_,
		upsamplingFactor_,
		baseElementIdx,
		firstRxElem_,
		*prepareDataTLS_,
#ifdef ARC_CYLINDRICAL_WAVE_PROCESSOR_USE_ANALYTIC_SIGNAL_FOR_ANGLE
		analyticSignalList_
#else
		signalList_
#endif
	};
	tbb::parallel_for(tbb::blocked_range<unsigned int>(firstRxElem_, lastRxElem_ + 1), prepareDataOp);
	//prepareDataOp(tbb::blocked_range<unsigned int>(firstRxElem_, lastRxElem_ + 1)); // single-thread
#ifdef USE_EXECUTION_TIME_MEASUREMENT
	tPartialPrepareDataML.put(prepareDataTimer.getTime());
#endif

	// Prepare the center signal.
	Timer t2;
	const unsigned int rxElem = txElem_;
	if (highUpsamplingFactor_ > 1) {
		// Interpolate the signal.
		centerInterpolator_.interpolate(&acqDataList_(baseElementIdx, rxElem - firstRxElem_, 0), samplesPerChannelLow, &centerSignal_[0]);
	} else {
		auto range = acqDataList_.range3(baseElementIdx, rxElem - firstRxElem_);
		std::copy(range.begin(), range.end(), centerSignal_.begin());
	}
	Util::removeDC(&centerSignal_[0], centerSignal_.size());
	envelope_.calculate(&centerSignal_[0], centerSignal_.size());
	LOG_DEBUG << "PREP 2 time = " << t2.getTime();

#ifdef USE_EXECUTION_TIME_MEASUREMENT
	Timer processTimer;
#endif
	// Find the radius of the center element arc using the center signal.
	TFloat maxValue = 0;
	int idxMax = 0;
	for (unsigned int i = 0; i < centerSignal_.size(); ++i) {
		if (centerSignal_[i] > maxValue) {
			maxValue = centerSignal_[i];
			idxMax = i;
		}
	}

	const TFloat tPulseEcho = (idxMax - centerSignalOffset_) / (config_.samplingFrequency * highUpsamplingFactor_);
	if (tPulseEcho <= 0) {
		THROW_EXCEPTION(InvalidValueException, "tPulseEcho <= 0");
	}
	const TFloat radius = (config_.propagationSpeed * tPulseEcho) / 2;
	LOG_DEBUG << "[RADIUS] baseElementIdx: " << baseElementIdx << " tPulseEcho: " << tPulseEcho << " radius: " << radius;

	const TFloat lambda = config_.propagationSpeed / config_.centerFrequency;
	Util::fillSequence(arcAngleList_, Util::degreeToRadian(minArcAngle_), Util::degreeToRadian(maxArcAngle_), (lambda / arcStepDiv_) / radius);

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
#ifdef ARC_CYLINDRICAL_WAVE_PROCESSOR_USE_ANALYTIC_SIGNAL_FOR_ANGLE
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

#ifdef USE_EXECUTION_TIME_MEASUREMENT
	tPartialProcessML.put(processTimer.getTime());
#endif

	//LOG_DEBUG << "END ========== ArcCylindricalWaveProcessor::process ==========";
}



template<typename TFloat>
struct ArcCylindricalWaveProcessor<TFloat>::PrepareData {
	void operator()(const tbb::blocked_range<unsigned int>& r) const {
		PrepareDataThreadData& local = prepareDataTLS.local();

		for (unsigned int rxElem = r.begin(); rxElem != r.end(); ++rxElem) {
			if (upsamplingFactor > 1) {
				// Interpolates the signal.
				local.interpolator.interpolate(&acqDataList(baseElementIdx, rxElem - firstRxElem, 0), samplesPerChannelLow, &local.tempSignal1[0]);
			} else {
				auto range = acqDataList.range3(baseElementIdx, rxElem - firstRxElem);
				std::copy(range.begin(), range.end(), local.tempSignal1.begin());
			}

			Util::removeDC(&local.tempSignal1[0], local.tempSignal1.size());

#ifdef ARC_CYLINDRICAL_WAVE_PROCESSOR_USE_ANALYTIC_SIGNAL_FOR_ANGLE
			local.envelope.getAnalyticSignal(&local.tempSignal1[0], local.tempSignal1.size(), &analyticSignalList(rxElem - firstRxElem, 0));
#else
			std::copy(local.tempSignal1.begin(), local.tempSignal1.end(), signalList.range2(rxElem - firstRxElem).begin());
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
#ifdef ARC_CYLINDRICAL_WAVE_PROCESSOR_USE_ANALYTIC_SIGNAL_FOR_ANGLE
	Matrix<std::complex<TFloat>>& analyticSignalList;
#else
	Matrix<TFloat>& signalList;
#endif
};



template<typename TFloat>
struct ArcCylindricalWaveProcessor<TFloat>::ProcessArcRange {
	void operator()(const tbb::blocked_range<std::size_t>& r) const {
		//LOG_DEBUG << "PRC col = " << r.begin() << " n = " << (r.end() - r.begin());

		ProcessArcRangeThreadData& data = processArcRangeTLS.local();

		data.rxSignalSumList.resize(config.numElements);
		const unsigned int numRxElem = lastRxElem - firstRxElem + 1;
		data.delayList.resize(numRxElem);

		// For each point in the arc:
		for (std::size_t i = r.begin(); i != r.end(); ++i) {
#ifdef ARC_CYLINDRICAL_WAVE_PROCESSOR_USE_ANALYTIC_SIGNAL_FOR_ANGLE
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
#ifdef ARC_CYLINDRICAL_WAVE_PROCESSOR_USE_ANALYTIC_SIGNAL_FOR_ANGLE
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
#ifdef ARC_CYLINDRICAL_WAVE_PROCESSOR_USE_ANALYTIC_SIGNAL_FOR_ANGLE
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
#ifdef ARC_CYLINDRICAL_WAVE_PROCESSOR_USE_ANALYTIC_SIGNAL_FOR_ANGLE
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

#endif // ARC_CYLINDRICAL_WAVE_PROCESSOR_H
