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

	struct PrepareDataThreadData {
		Interpolator<TFloat> interpolator;
		HilbertEnvelope<TFloat> envelope;
		std::vector<TFloat, tbb::cache_aligned_allocator<TFloat>> signal;
	};

	struct DelaySumThreadData {
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
	unsigned int signalLength_;
	TFloat signalOffset_;
	bool initialized_;
	const std::vector<TFloat> txApod_;
	const std::vector<TFloat> rxApod_;
	std::vector<TFloat, tbb::cache_aligned_allocator<TFloat>> xArray_;
	std::unique_ptr<tbb::enumerable_thread_specific<PrepareDataThreadData>> prepareDataTLS_;
	std::unique_ptr<tbb::enumerable_thread_specific<DelaySumThreadData>> delaySumTLS_;
};



template<typename TFloat, typename TPoint>
VectorialSTAProcessor<TFloat, TPoint>::VectorialSTAProcessor(
			const STAConfiguration<TFloat>& config,
			STAAcquisition<TFloat>& acquisition,
			unsigned int upsamplingFactor,
			AnalyticSignalCoherenceFactorProcessor<TFloat>& coherenceFactor,
			TFloat peakOffset,
			const std::vector<TFloat>& txApod,
			const std::vector<TFloat>& rxApod)
		: config_(config)
		, deadZoneSamplesUp_((upsamplingFactor * config.samplingFrequency) * 2.0 * config.deadZoneM / config.propagationSpeed)
		, acquisition_(acquisition)
		, upsamplingFactor_(upsamplingFactor)
		, coherenceFactor_(coherenceFactor)
		, signalLength_()
		, signalOffset_()
		, initialized_()
		, txApod_(txApod)
		, rxApod_(rxApod)
{
	signalOffset_ = (config_.samplingFrequency * upsamplingFactor_) * peakOffset / config_.centerFrequency;
	LOG_DEBUG << "signalOffset: " << signalOffset_;
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
	//LOG_DEBUG << "BEGIN ========== VectorialSTAProcessor::process ==========";

	Util::resetValueFactor(gridData.begin(), gridData.end());

	const std::size_t numSignals = (config_.lastTxElem - config_.firstTxElem + 1U) * config_.numElements;

	// Prepare the signal tensor.
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
			signalLength_ = acqData_.n2() * upsamplingFactor_;
			LOG_DEBUG << "signalLength: " << signalLength_;
			if (deadZoneSamplesUp_ >= signalLength_) {
				THROW_EXCEPTION(InvalidValueException,
						"Invalid dead zone: deadZoneSamplesUp_ (" << deadZoneSamplesUp_ <<
						") >= signalLength (" << signalLength_ << ").");
			}

			PrepareDataThreadData prepareDataThreadData;
			if (upsamplingFactor_ > 1) {
				prepareDataThreadData.interpolator.prepare(upsamplingFactor_, upsampFilterHalfTransitionWidth);
			}
			prepareDataThreadData.signal.resize(signalLength_);
			prepareDataTLS_ = std::make_unique<tbb::enumerable_thread_specific<PrepareDataThreadData>>(prepareDataThreadData);

			DelaySumThreadData delaySumThreadData;
			delaySumThreadData.coherenceFactor = coherenceFactor_;
			delaySumTLS_ = std::make_unique<tbb::enumerable_thread_specific<DelaySumThreadData>>(delaySumThreadData);

			analyticSignalTensor_.resize(
						config_.lastTxElem - config_.firstTxElem + 1,
						config_.numElements, signalLength_);

			initialized_ = true;
		}

		const unsigned int samplesPerChannelLow = acqData_.n2();

#ifdef USE_EXECUTION_TIME_MEASUREMENT
		Timer prepareDataTimer;
#endif
		tbb::parallel_for(tbb::blocked_range<unsigned int>(0, config_.numElements),
		[&, samplesPerChannelLow, txElem](const tbb::blocked_range<unsigned int>& r) {
			auto& local = prepareDataTLS_->local();

			for (unsigned int rxElem = r.begin(); rxElem != r.end(); ++rxElem) {
				if (upsamplingFactor_ > 1) {
					local.interpolator.interpolate(&acqData_(rxElem, 0), samplesPerChannelLow, local.signal.data());
				} else {
					auto range = acqData_.range2(rxElem);
					std::copy(range.begin(), range.end(), local.signal.begin());
				}

				Util::removeDC(local.signal.data(), local.signal.size(), deadZoneSamplesUp_);

				local.envelope.getAnalyticSignal(
						local.signal.data(),
						local.signal.size(),
						&analyticSignalTensor_(txElem - config_.firstTxElem, rxElem, 0));
			}
		});
#ifdef USE_EXECUTION_TIME_MEASUREMENT
		tPrepareDataML.put(prepareDataTimer.getTime());
#endif
	}

	if (xArray_.empty()) {
		ArrayGeometry::getElementXCentered2D(config_.numElements, config_.pitch, xArray_);
	}

	const TFloat invCT = (config_.samplingFrequency * upsamplingFactor_) / config_.propagationSpeed;
	const unsigned int numCols = gridData.n1();
	const unsigned int numRows = gridData.n2();

	IterationCounter::reset(numCols);

#ifdef USE_EXECUTION_TIME_MEASUREMENT
	Timer delaySumTimer;
#endif
	if (coherenceFactor_.enabled()) {
		tbb::parallel_for(tbb::blocked_range<unsigned int>(0, numCols),
		[&, invCT, numRows](const tbb::blocked_range<unsigned int>& r) {
			auto& local = delaySumTLS_->local();

			local.rxSignalSumList.resize(config_.numElements);
			local.delayList.resize(config_.numElements);

			for (unsigned int col = r.begin(); col != r.end(); ++col) {
				for (unsigned int row = 0; row < numRows; ++row) {

					std::fill(local.rxSignalSumList.begin(), local.rxSignalSumList.end(), std::complex<TFloat>(0));
					TPoint& point = gridData(col, row);

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
										((1 - k) * *p + k * *(p + 1))
										* txApod_[txElem] * rxApod_[rxElem];
							}
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
	} else {
		tbb::parallel_for(tbb::blocked_range<unsigned int>(0, numCols),
		[&, invCT, numRows](const tbb::blocked_range<unsigned int>& r) {
			auto& local = delaySumTLS_->local();

			local.delayList.resize(config_.numElements);

			for (unsigned int col = r.begin(); col != r.end(); ++col) {
				for (unsigned int row = 0; row < numRows; ++row) {
					TPoint& point = gridData(col, row);

					// Calculate the delays.
					for (unsigned int elem = 0; elem < config_.numElements; ++elem) {
						local.delayList[elem] = Geometry::distance2DY0(xArray_[elem], point.x, point.z) * invCT;
					}

					std::complex<TFloat> signalSum;
					for (unsigned int txElem = config_.firstTxElem; txElem <= config_.lastTxElem; ++txElem) {
						const TFloat txDelay = local.delayList[txElem];
						const unsigned int localTxElem = txElem - config_.firstTxElem;
						std::complex<TFloat> rxSignalSum;
						for (unsigned int rxElem = 0; rxElem < config_.numElements; ++rxElem) {
							// Linear interpolation.
							const TFloat delay = signalOffset_ + txDelay + local.delayList[rxElem];
							const std::size_t delayIdx = static_cast<std::size_t>(delay);
							const TFloat k = delay - delayIdx;
							if (delayIdx + 1U < analyticSignalTensor_.n3()) {
								const std::complex<TFloat>* p = &analyticSignalTensor_(localTxElem, rxElem, delayIdx);
								rxSignalSum += ((1 - k) * *p + k * *(p + 1)) * rxApod_[rxElem];
							}
						}
						signalSum += rxSignalSum * txApod_[txElem];
					}
					point.value = std::abs(signalSum);
				}
			}

			IterationCounter::add(r.end() - r.begin());
		});
	}
#ifdef USE_EXECUTION_TIME_MEASUREMENT
	tDelaySumML.put(delaySumTimer.getTime());
#endif

	std::for_each(gridData.begin(), gridData.end(),
			Util::MultiplyValueBy<TPoint, TFloat>(TFloat(1) / numSignals));

	//LOG_DEBUG << "END ========== VectorialSTAProcessor::process ==========";
}

} // namespace Lab

#endif /* VECTORIALSTAPROCESSOR_H_ */
