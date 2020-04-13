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

#ifndef TANGENTCURVEGEOMETRYPROCESSOR_H
#define TANGENTCURVEGEOMETRYPROCESSOR_H

#include <algorithm> /* copy, reverse, reverse_copy */
#include <cmath> /* abs, acos, cos, sin */
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
#define TANGENT_CURVE_GEOMETRY_PROCESSOR_UPSAMP_FILTER_HALF_TRANSITION_WIDTH 0.2

#define TANGENT_CURVE_GEOMETRY_PROCESSOR_USE_CROSS_CORRELATION 1
//#define TANGENT_CURVE_GEOMETRY_PROCESSOR_USE_ENVELOPE_OF_CROSS_CORRELATION 1

#if !defined(TANGENT_CURVE_GEOMETRY_PROCESSOR_USE_CROSS_CORRELATION) || defined(TANGENT_CURVE_GEOMETRY_PROCESSOR_USE_ENVELOPE_OF_CROSS_CORRELATION)
# include "HilbertEnvelope.h"
#endif

//#define TANGENT_CURVE_GEOMETRY_PROCESSOR_GET_MAX_VALUE 1



namespace Lab {

// x = 0 is at the center of the element group.
// y = 0
// z = 0 is at the surface of the array.
template<typename TFloat>
class TangentCurveGeometryProcessor {
public:
	TangentCurveGeometryProcessor(
			const STAConfiguration<TFloat>& config,
			const Tensor3<TFloat>& acqDataList,
			unsigned int upsamplingFactor,
			unsigned int txElem,
			unsigned int firstRxElem,
			unsigned int lastRxElem,
			TFloat arcStep,
#ifdef TANGENT_CURVE_GEOMETRY_PROCESSOR_USE_CROSS_CORRELATION
			const std::vector<TFloat>& refPulse,
#else
			TFloat peakOffset,
#endif
			unsigned int signalStartOffset);

	~TangentCurveGeometryProcessor() = default;

	void process(
		unsigned int baseElementIdx,
		unsigned int baseElement,
		std::vector<std::pair<TFloat, TFloat>>& pointPositionList);

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
#ifdef TANGENT_CURVE_GEOMETRY_PROCESSOR_USE_CROSS_CORRELATION
//		FFTWFilter<TFloat> revRefPulseFilter;
		DirectFFTWFilter<TFloat> revRefPulseFilter;
		std::vector<TFloat> tempSignal1;
#endif
#if !defined(TANGENT_CURVE_GEOMETRY_PROCESSOR_USE_CROSS_CORRELATION) || defined(TANGENT_CURVE_GEOMETRY_PROCESSOR_USE_ENVELOPE_OF_CROSS_CORRELATION)
		HilbertEnvelope<TFloat> envelope;
#endif
		std::vector<TFloat> signal;
	};

	TangentCurveGeometryProcessor(const TangentCurveGeometryProcessor&) = delete;
	TangentCurveGeometryProcessor& operator=(const TangentCurveGeometryProcessor&) = delete;
	TangentCurveGeometryProcessor(TangentCurveGeometryProcessor&&) = delete;
	TangentCurveGeometryProcessor& operator=(TangentCurveGeometryProcessor&&) = delete;

	const STAConfiguration<TFloat>& config_;
	const Tensor3<TFloat>& acqDataList_;
	unsigned int upsamplingFactor_;
	unsigned int txElem_;
	unsigned int firstRxElem_;
	unsigned int lastRxElem_;
	TFloat arcStep_;
	TFloat signalOffset_;
	std::vector<TFloat> distRxElem_;
	std::vector<TFloat> intersecAngleList_;
	std::vector<TFloat> arcAngleList_; //TODO: reserve?
	std::vector<TFloat> arcRadiusList_; //TODO: reserve?
	std::unique_ptr<tbb::enumerable_thread_specific<PrepareDataThreadData>> prepareDataTLS_;
};



template<typename TFloat>
TangentCurveGeometryProcessor<TFloat>::TangentCurveGeometryProcessor(
			const STAConfiguration<TFloat>& config,
			const Tensor3<TFloat>& acqDataList,
			unsigned int upsamplingFactor,
			unsigned int txElem,
			unsigned int firstRxElem,
			unsigned int lastRxElem,
			TFloat arcStep,
#ifdef TANGENT_CURVE_GEOMETRY_PROCESSOR_USE_CROSS_CORRELATION
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
		, arcStep_(arcStep)
{
	const std::size_t origSignalLength = acqDataList_.n3();

	PrepareDataThreadData prepareDataThreadData;
	if (upsamplingFactor_ > 1) {
		prepareDataThreadData.interpolator.prepare(upsamplingFactor_, TANGENT_CURVE_GEOMETRY_PROCESSOR_UPSAMP_FILTER_HALF_TRANSITION_WIDTH);
	}

#ifdef TANGENT_CURVE_GEOMETRY_PROCESSOR_USE_CROSS_CORRELATION
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
	intersecAngleList_.assign(numActiveRxElements, 0); // the position (txElem - firstRxElem) will always contain zero

#ifdef TANGENT_CURVE_GEOMETRY_PROCESSOR_USE_CROSS_CORRELATION
	prepareDataThreadData.revRefPulseFilter.setCoefficients(revRefPulse);
	prepareDataThreadData.tempSignal1.resize(origSignalLength * upsamplingFactor_);
#else
	prepareDataThreadData.signal.resize(origSignalLength * upsamplingFactor_);
#endif
	prepareDataTLS_ = std::make_unique<tbb::enumerable_thread_specific<PrepareDataThreadData>>(prepareDataThreadData);
}

template<typename TFloat>
void
TangentCurveGeometryProcessor<TFloat>::process(
		unsigned int baseElementIdx,
		unsigned int baseElement,
		std::vector<std::pair<TFloat, TFloat>>& pointPositionList)
{
	//LOG_DEBUG << "BEGIN ========== TangentCurveGeometryProcessor::process ==========";

	const std::size_t samplesPerChannelLow = acqDataList_.n3();

	// Obtain the traveled distances.
#ifdef USE_EXECUTION_TIME_MEASUREMENT
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
#ifdef TANGENT_CURVE_GEOMETRY_PROCESSOR_GET_MAX_VALUE
	prepareDataOp(tbb::blocked_range<unsigned int>(firstRxElem_, lastRxElem_ + 1)); // single-thread
#else
	tbb::parallel_for(tbb::blocked_range<unsigned int>(firstRxElem_, lastRxElem_ + 1), prepareDataOp);
#endif
#ifdef USE_EXECUTION_TIME_MEASUREMENT
	tPartialPrepareDataML.put(prepareDataTimer.getTime());

	Timer processTimer;
#endif

	std::vector<TFloat> c(config_.numElements - 1, 0);
	for (unsigned int rxElem = firstRxElem_; rxElem <= lastRxElem_; ++rxElem) {
		if (txElem_ >= rxElem) {
			c[txElem_ - rxElem] = ((txElem_ - rxElem) * config_.pitch) / 2;
		} else {
			c[rxElem - txElem_] = ((rxElem - txElem_) * config_.pitch) / 2;
		}
	}

	const TFloat txRadius = distRxElem_[txElem_ - firstRxElem_] / 2; // radius of the circle centered on the tx element
	TFloat minIntersecAngle = pi;
	TFloat maxIntersecAngle = 0;

	//TODO: More than one intersection? #########
	// Find the intersections between the tx elem. circle and the ellipses.
	for (unsigned int rxElem = firstRxElem_; rxElem < txElem_; ++rxElem) {
		const unsigned rxElemIdx = rxElem - firstRxElem_;
		const unsigned int dist = txElem_ - rxElem;
		const TFloat a = distRxElem_[rxElemIdx] / 2;
		const TFloat e = c[dist] / a;
		const TFloat cosTheta = (a * (1 - e * e) - txRadius) / (txRadius * e);
//		if (std::abs(cosTheta) > 1) { // would generate NaN
//			intersecAngleList_[rxElemIdx] = PI;
//		} else {
			const TFloat angle = std::acos(cosTheta);
			if (angle > maxIntersecAngle) maxIntersecAngle = angle;
			intersecAngleList_[rxElemIdx] = angle;
//		}
	}
	for (unsigned int rxElem = txElem_ + 1; rxElem <= lastRxElem_; ++rxElem) {
		const unsigned rxElemIdx = rxElem - firstRxElem_;
		const unsigned int dist = rxElem - txElem_;
		const TFloat a = distRxElem_[rxElemIdx] / 2;
		const TFloat e = c[dist] / a;
		const TFloat cosTheta = -(a * (1 - e * e) - txRadius) / (txRadius * e);
//		if (std::abs(cosTheta) > 1) { // would generate NaN
//			intersecAngleList_[rxElemIdx] = 0;
//		} else {
			const TFloat angle = std::acos(cosTheta);
			if (angle < minIntersecAngle) minIntersecAngle = angle;
			intersecAngleList_[rxElemIdx] = angle;
//		}
	}
	//LOG_DEBUG << "intersecAngleList = " << intersecAngleList;

	Util::fillSequenceFromStartToEndWithMaximumStep(arcAngleList_,
							minIntersecAngle, maxIntersecAngle,
							arcStep_ / txRadius);
	arcRadiusList_.assign(arcAngleList_.size(), txRadius);

	// Search in the set of ellipses for the points that are farthest from the tx element.
	for (unsigned int rxElem = firstRxElem_; rxElem < txElem_; ++rxElem) {
		const unsigned rxElemIdx = rxElem - firstRxElem_;
		const unsigned int dist = txElem_ - rxElem;
		const TFloat a = distRxElem_[rxElemIdx] / 2;
		const TFloat e = c[dist] / a;
		const TFloat kAngle = a * (1 - e * e);
		for (unsigned int i = 0; i < arcAngleList_.size(); ++i) {
			if (arcAngleList_[i] > intersecAngleList_[rxElemIdx]) {
				const TFloat r = kAngle / (1 + e * std::cos(arcAngleList_[i]));
				if (r > arcRadiusList_[i]) arcRadiusList_[i] = r;
			}
		}
	}
	for (unsigned int rxElem = txElem_ + 1; rxElem <= lastRxElem_; ++rxElem) {
		const unsigned rxElemIdx = rxElem - firstRxElem_;
		const unsigned int dist = rxElem - txElem_;
		const TFloat a = distRxElem_[rxElemIdx] / 2;
		const TFloat e = c[dist] / a;
		const TFloat kAngle = a * (1 - e * e);
		for (unsigned int i = 0; i < arcAngleList_.size(); ++i) {
			if (arcAngleList_[i] < intersecAngleList_[rxElemIdx]) {
				const TFloat r = kAngle / (1 - e * std::cos(arcAngleList_[i]));
				if (r > arcRadiusList_[i]) arcRadiusList_[i] = r;
			}
		}
	}
	// At this point, arcRadiusList[i] may be = txRadius even if it is not the tangent point,
	// if too many data are invalid (low level echoes / noise).

	for (unsigned int i = 0; i < arcAngleList_.size(); ++i) {
		const TFloat x = (baseElement + txElem_) * config_.pitch + arcRadiusList_[i] * std::cos(arcAngleList_[i]);
		const TFloat z = arcRadiusList_[i] * std::sin(arcAngleList_[i]);
		pointPositionList.push_back(std::make_pair(x, z));
	}

#ifdef USE_EXECUTION_TIME_MEASUREMENT
	tPartialProcessML.put(processTimer.getTime());
#endif

	//LOG_DEBUG << "END ========== TangentCurveGeometryProcessor::process ==========";
}



template<typename TFloat>
struct TangentCurveGeometryProcessor<TFloat>::PrepareData {
	void operator()(const tbb::blocked_range<unsigned int>& r) const {
		PrepareDataThreadData& local = prepareDataTLS.local();

#ifdef TANGENT_CURVE_GEOMETRY_PROCESSOR_GET_MAX_VALUE
		TFloat minMaxValue = 1e9;
#endif

		for (unsigned int rxElem = r.begin(); rxElem != r.end(); ++rxElem) {
#ifdef TANGENT_CURVE_GEOMETRY_PROCESSOR_USE_CROSS_CORRELATION
			if (upsamplingFactor > 1) {
				// Interpolate the signal.
				local.interpolator.interpolate(&acqDataList(baseElementIdx, rxElem - firstRxElem, 0), samplesPerChannelLow, &local.tempSignal1[0]);
			} else {
				auto range = acqDataList.range3(baseElementIdx, rxElem - firstRxElem);
				std::copy(range.begin(), range.end(), local.tempSignal1.begin());
			}

			// Cross-correlation using convolution.
			local.revRefPulseFilter.filter(local.tempSignal1, local.signal);

# ifdef TANGENT_CURVE_GEOMETRY_PROCESSOR_USE_ENVELOPE_OF_CROSS_CORRELATION
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
#ifdef TANGENT_CURVE_GEOMETRY_PROCESSOR_GET_MAX_VALUE
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
#ifdef TANGENT_CURVE_GEOMETRY_PROCESSOR_GET_MAX_VALUE
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

#endif // TANGENTCURVEGEOMETRYPROCESSOR_H
