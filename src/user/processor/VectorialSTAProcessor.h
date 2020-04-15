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

#ifndef VECTORIALSTAPROCESSOR_H_
#define VECTORIALSTAPROCESSOR_H_

#include <algorithm> /* for_each */
#include <cmath>
#include <complex>
#include <cstddef> /* std::size_t */
#include <numeric> /* accumulate */
#include <vector>

#include <tbb/enumerable_thread_specific.h>
#include <tbb/tbb.h>

#include "ArrayGeometry.h"
#include "ArrayProcessor.h"
#include "CoherenceFactor.h"
#include "Exception.h"
#include "ExecutionTimeMeasurement.h"
#include "Geometry.h"
#include "HilbertEnvelope.h"
#include "Interpolator.h"
#include "IterationCounter.h"
#include "Log.h"
#include "Matrix.h"
#include "STAAcquisition.h"
#include "STAConfiguration.h"
#include "Tensor3.h"
#include "Util.h"



namespace Lab {

// STA image formation, using analytic signals (each sample is a real-imag vector).
//
// x = 0 is at the center of the element group.
// y = 0
// z = 0 is at the surface of the array.
template<typename TFloat, typename TPoint>
class VectorialSTAProcessor : public ArrayProcessor<TPoint> {
public:
	VectorialSTAProcessor(
			const STAConfiguration<TFloat>& config,
			STAAcquisition<TFloat>& acquisition,
			unsigned int upsamplingFactor,
			AnalyticSignalCoherenceFactorProcessor<TFloat>& coherenceFactor,
			TFloat peakOffset,
			bool calculateEnvelope,
			const std::vector<TFloat>& txApod,
			const std::vector<TFloat>& rxApod);
	virtual ~VectorialSTAProcessor() = default;

	virtual void prepare(unsigned int baseElement);
	virtual void process(Matrix<TPoint>& gridData);

#ifdef USE_EXECUTION_TIME_MEASUREMENT
	MeasurementList<double> tAcquisitionML;
	MeasurementList<double> tPrepareDataML;
	MeasurementList<double> tDelaySumML;
	void execTimeMeasReset(unsigned int n) {
		tAcquisitionML.reset(    EXECUTION_TIME_MEASUREMENT_ITERATIONS * n);
		tPrepareDataML.reset(    EXECUTION_TIME_MEASUREMENT_ITERATIONS * n);
		tDelaySumML.reset(       EXECUTION_TIME_MEASUREMENT_ITERATIONS);
	}
	void execTimeMeasShowResults(unsigned int n) {
		EXECUTION_TIME_MEASUREMENT_LOG_TIMES_X_N("tAcquisition:    ", tAcquisitionML, n);
		EXECUTION_TIME_MEASUREMENT_LOG_TIMES_X_N("tPrepareData:    ", tPrepareDataML, n);
		EXECUTION_TIME_MEASUREMENT_LOG_TIMES(    "tDelaySum:       ", tDelaySumML);
	}
#endif

private:
	// Depends on the signal.
	// 1.0 --> pi radian / sample at the original sampling rate.
	static constexpr TFloat upsampFilterHalfTransitionWidth = 0.2;

	struct ThreadData {
		AnalyticSignalCoherenceFactorProcessor<TFloat> coherenceFactor;
		std::vector<std::complex<TFloat>> rxSignalSumList;
		std::vector<TFloat> delayList;
	};

	VectorialSTAProcessor(const VectorialSTAProcessor&) = delete;
	VectorialSTAProcessor& operator=(const VectorialSTAProcessor&) = delete;
	VectorialSTAProcessor(VectorialSTAProcessor&&) = delete;
	VectorialSTAProcessor& operator=(VectorialSTAProcessor&&) = delete;

	const STAConfiguration<TFloat>& config_;
	unsigned int deadZoneSamplesUp_;
	STAAcquisition<TFloat>& acquisition_;
	unsigned int upsamplingFactor_;
	AnalyticSignalCoherenceFactorProcessor<TFloat>& coherenceFactor_;
	Tensor3<std::complex<TFloat>> analyticSignalTensor_;
	typename STAAcquisition<TFloat>::AcquisitionDataType acqData_;
	std::vector<TFloat> tempSignal_;
	TFloat signalOffset_;
	Interpolator<TFloat> interpolator_;
	HilbertEnvelope<TFloat> envelope_;
	bool calculateEnvelope_;
	bool initialized_;
	std::vector<TFloat> txApod_;
	std::vector<TFloat> rxApod_;
	std::vector<TFloat, tbb::cache_aligned_allocator<TFloat>> xArray_;
};



template<typename TFloat, typename TPoint>
VectorialSTAProcessor<TFloat, TPoint>::VectorialSTAProcessor(
			const STAConfiguration<TFloat>& config,
			STAAcquisition<TFloat>& acquisition,
			unsigned int upsamplingFactor,
			AnalyticSignalCoherenceFactorProcessor<TFloat>& coherenceFactor,
			TFloat peakOffset,
			bool calculateEnvelope,
			const std::vector<TFloat>& txApod,
			const std::vector<TFloat>& rxApod)
		: config_(config)
		, deadZoneSamplesUp_((upsamplingFactor * config.samplingFrequency) * 2.0 * config.deadZoneM / config.propagationSpeed)
		, acquisition_(acquisition)
		, upsamplingFactor_(upsamplingFactor)
		, coherenceFactor_(coherenceFactor)
		, calculateEnvelope_(calculateEnvelope)
		, initialized_()
		, txApod_(txApod)
		, rxApod_(rxApod)
{
	if (upsamplingFactor_ > 1) {
		interpolator_.prepare(upsamplingFactor_, upsampFilterHalfTransitionWidth);
	}

	signalOffset_ = (config_.samplingFrequency * upsamplingFactor_) * peakOffset / config_.centerFrequency;
}

template<typename TFloat, typename TPoint>
void
VectorialSTAProcessor<TFloat, TPoint>::prepare(unsigned int baseElement)
{
	acquisition_.prepare(baseElement);
}

template<typename TFloat, typename TPoint>
void
VectorialSTAProcessor<TFloat, TPoint>::process(Matrix<TPoint>& gridData)
{
	LOG_DEBUG << "BEGIN ========== VectorialSTAProcessor::process ==========";

	Util::resetValueFactor(gridData.begin(), gridData.end());

	const std::size_t numSignals = (config_.lastTxElem - config_.firstTxElem + 1U) * config_.numElements;

	// Prepare the signal matrix.
	for (unsigned int txElem = config_.firstTxElem; txElem <= config_.lastTxElem; ++txElem) {
		//LOG_INFO << "ACQ/PREP txElem: " << txElem << " <= " << config_.lastTxElem;

#ifdef USE_EXECUTION_TIME_MEASUREMENT
		Timer acquisitionTimer;
#endif
		acquisition_.execute(txElem, acqData_);

#ifdef USE_EXECUTION_TIME_MEASUREMENT
		tAcquisitionML.put(acquisitionTimer.getTime());
#endif
		if (!initialized_) {
			const std::size_t signalLength = acqData_.n2() * upsamplingFactor_;
			tempSignal_.resize(signalLength);
			analyticSignalTensor_.resize(
						config_.lastTxElem - config_.firstTxElem + 1,
						config_.numElements, signalLength);
			LOG_DEBUG << "signalOffset_=" << signalOffset_ << " signalLength=" << signalLength;
			if (deadZoneSamplesUp_ >= signalLength) {
				THROW_EXCEPTION(InvalidValueException,
						"Invalid dead zone: deadZoneSamplesUp_ (" << deadZoneSamplesUp_ <<
						") >= signalLength (" << signalLength << ").");
			}
			initialized_ = true;
		}

		const std::size_t samplesPerChannelLow = acqData_.n2();

#ifdef USE_EXECUTION_TIME_MEASUREMENT
		Timer prepareDataTimer;
#endif
		const unsigned int localTxElem = txElem - config_.firstTxElem;
		for (unsigned int rxElem = 0; rxElem < config_.numElements; ++rxElem) {
			if (upsamplingFactor_ > 1) {
				interpolator_.interpolate(&acqData_(rxElem, 0), samplesPerChannelLow, &tempSignal_[0]);
			} else {
				auto range = acqData_.range2(rxElem);
				std::copy(range.begin(), range.end(), tempSignal_.begin());
			}

			Util::removeDC(&tempSignal_[0], tempSignal_.size(), deadZoneSamplesUp_);

			envelope_.getAnalyticSignal(&tempSignal_[0], tempSignal_.size(), &analyticSignalTensor_(localTxElem, rxElem, 0));
		}

#ifdef USE_EXECUTION_TIME_MEASUREMENT
		tPrepareDataML.put(prepareDataTimer.getTime());
#endif
	}

	if (xArray_.empty()) {
		ArrayGeometry::getElementXCentered2D(config_.numElements, config_.pitch, xArray_);
	}

	ThreadData threadData;
	threadData.coherenceFactor = coherenceFactor_;
	tbb::enumerable_thread_specific<ThreadData> tls(threadData);

	const TFloat invCT = (config_.samplingFrequency * upsamplingFactor_) / config_.propagationSpeed;
	const std::size_t numRows = gridData.n2();

	IterationCounter::reset(gridData.n1());

#ifdef USE_EXECUTION_TIME_MEASUREMENT
	Timer delaySumTimer;
#endif
	tbb::parallel_for(tbb::blocked_range<std::size_t>(0, gridData.n1()),
	[&, invCT, numRows](const tbb::blocked_range<std::size_t>& r) {
		auto& local = tls.local();

		local.rxSignalSumList.resize(config_.numElements);
		local.delayList.resize(config_.numElements);

		// For each column:
		for (std::size_t i = r.begin(); i != r.end(); ++i) {
			// For each row:
			for (std::size_t j = 0; j < numRows; ++j) {

				std::fill(local.rxSignalSumList.begin(), local.rxSignalSumList.end(), std::complex<TFloat>(0));
				TPoint& point = gridData(i, j);

				// Calculate the delays.
				for (unsigned int elem = 0; elem < config_.numElements; ++elem) {
					local.delayList[elem] = Geometry::distance2DY0(xArray_[elem], point.x, point.z) * invCT;
				}

				for (unsigned int txElem = config_.firstTxElem; txElem <= config_.lastTxElem; ++txElem) {
					const TFloat txDelay = local.delayList[txElem];
					const unsigned int localTxElem = txElem - config_.firstTxElem;
					for (unsigned int rxElem = 0; rxElem < config_.numElements; ++rxElem) {
						// Linear interpolation.
						const TFloat delay = signalOffset_ + txDelay + local.delayList[rxElem];
						const std::size_t delayIdx = static_cast<std::size_t>(delay);
						const TFloat k = delay - delayIdx;
						if (delayIdx + 1U < analyticSignalTensor_.n3()) {
							const std::complex<TFloat>* p = &analyticSignalTensor_(localTxElem, rxElem, delayIdx);
							local.rxSignalSumList[rxElem] +=
									txApod_[txElem] * rxApod_[rxElem]
									* ((1 - k) * *p + k * *(p + 1));
						}
					}
				}

				if (local.coherenceFactor.enabled()) {
					point.factor = local.coherenceFactor.calculate(&local.rxSignalSumList[0], local.rxSignalSumList.size());
				}
				if (calculateEnvelope_) {
					point.value = std::abs(std::accumulate(local.rxSignalSumList.begin(), local.rxSignalSumList.end(), std::complex<TFloat>(0)));
				} else {
					point.value = std::accumulate(local.rxSignalSumList.begin(), local.rxSignalSumList.end(), std::complex<TFloat>(0)).real();
				}
			}
		}

		IterationCounter::add(r.end() - r.begin());
	});
#ifdef USE_EXECUTION_TIME_MEASUREMENT
	tDelaySumML.put(delaySumTimer.getTime());
#endif

	std::for_each(gridData.begin(), gridData.end(),
			Util::MultiplyValueBy<TPoint, TFloat>(TFloat(1) / numSignals));

	LOG_DEBUG << "END ========== VectorialSTAProcessor::process ==========";
}

} // namespace Lab

#endif /* VECTORIALSTAPROCESSOR_H_ */
