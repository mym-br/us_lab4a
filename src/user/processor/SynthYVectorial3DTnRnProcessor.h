/***************************************************************************
 *  Copyright 2014, 2017, 2018, 2019 Marcelo Y. Matuda                     *
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

#ifndef SYNTH_Y_SYNTH_Y_VECTORIAL_3D_TN_RN_PROCESSOR_H
#define SYNTH_Y_SYNTH_Y_VECTORIAL_3D_TN_RN_PROCESSOR_H

#include <algorithm> /* for_each */
#include <cmath>
#include <complex>
#include <cstddef> /* std::size_t */
#include <numeric> /* accumulate */
#include <vector>

#include <tbb/enumerable_thread_specific.h>
#include <tbb/tbb.h>

#include "CoherenceFactor.h"
#include "Exception.h"
#include "Geometry.h"
#include "HilbertEnvelope.h"
#include "Interpolator.h"
#include "IterationCounter.h"
#include "Log.h"
#include "Matrix.h"
#include "TnRnAcquisition.h"
#include "TnRnConfiguration.h"
#include "Util.h"
#include "XY.h"
#include "XYZValueFactor.h"



namespace Lab {

// abs(focusZ_) must be much smaller than the distance between the array and the reflector(s).
template<typename FloatType>
class SynthYVectorial3DTnRnProcessor {
public:
	SynthYVectorial3DTnRnProcessor(
			const TnRnConfiguration<FloatType>& config,
			TnRnAcquisition<FloatType>& acquisition,
			unsigned int upsamplingFactor,
			AnalyticSignalCoherenceFactorProcessor<FloatType>& coherenceFactor,
			FloatType peakOffset,
			const std::vector<FloatType>& rxApod);
	virtual ~SynthYVectorial3DTnRnProcessor() = default;

	void setTxDelays(FloatType focusX, FloatType focusY, FloatType focusZ /* can be negative */,
				const std::vector<FloatType>& txDelays);
	void prepare(unsigned int baseElement);
	void getAcqData(FloatType y);
	void process(unsigned int firstAcq, unsigned int lastAcq,
			Matrix<XYZValueFactor<FloatType>>& gridData);
private:
	// Depends on the signal.
	// 1.0 --> pi radian / sample at the original sampling rate.
	static constexpr FloatType upsampFilterHalfTransitionWidth = 0.2;

	struct ThreadData {
		AnalyticSignalCoherenceFactorProcessor<FloatType> coherenceFactor;
		std::vector<std::complex<FloatType>> rxSignalSumList;
		std::vector<FloatType> rxDelayList;
	};
	struct AcqData {
		Matrix<std::complex<FloatType>> analyticSignalMatrix;
		FloatType y;
	};

	SynthYVectorial3DTnRnProcessor(const SynthYVectorial3DTnRnProcessor&) = delete;
	SynthYVectorial3DTnRnProcessor& operator=(const SynthYVectorial3DTnRnProcessor&) = delete;
	SynthYVectorial3DTnRnProcessor(SynthYVectorial3DTnRnProcessor&&) = delete;
	SynthYVectorial3DTnRnProcessor& operator=(SynthYVectorial3DTnRnProcessor&&) = delete;

	const TnRnConfiguration<FloatType>& config_;
	unsigned int deadZoneSamplesUp_;
	TnRnAcquisition<FloatType>& acquisition_;
	unsigned int upsamplingFactor_;
	AnalyticSignalCoherenceFactorProcessor<FloatType>& coherenceFactor_;
	std::vector<AcqData> acqData_; // this data may require lots of memory
	typename TnRnAcquisition<FloatType>::AcquisitionDataType singleAcqData_;
	std::vector<FloatType> tempSignal_;
	FloatType signalOffset_;
	Interpolator<FloatType> interpolator_;
	HilbertEnvelope<FloatType> envelope_;
	bool initialized_;
	std::vector<FloatType> rxApod_;
	FloatType focusX_;
	FloatType focusY_;
	FloatType focusZ_;
	std::vector<FloatType> txDelays_; // focalization / divergence delays
	unsigned int baseElement_;
};



template<typename FloatType>
SynthYVectorial3DTnRnProcessor<FloatType>::SynthYVectorial3DTnRnProcessor(
			const TnRnConfiguration<FloatType>& config,
			TnRnAcquisition<FloatType>& acquisition,
			unsigned int upsamplingFactor,
			AnalyticSignalCoherenceFactorProcessor<FloatType>& coherenceFactor,
			FloatType peakOffset,
			const std::vector<FloatType>& rxApod)
		: config_(config)
		, deadZoneSamplesUp_((upsamplingFactor * config.samplingFrequency) * 2.0 * config.deadZoneM / config.propagationSpeed)
		, acquisition_(acquisition)
		, upsamplingFactor_(upsamplingFactor)
		, coherenceFactor_(coherenceFactor)
		, initialized_()
		, rxApod_(rxApod)
		, focusX_()
		, focusY_()
		, focusZ_()
		, baseElement_()
{
	if (upsamplingFactor_ > 1) {
		interpolator_.prepare(upsamplingFactor_, upsampFilterHalfTransitionWidth);
	}

	signalOffset_ = (config_.samplingFrequency * upsamplingFactor_) * peakOffset / config_.centerFrequency;

	if (rxApod_.size() != config_.numElements) {
		THROW_EXCEPTION(InvalidValueException, "Wrong rx apodization size: " << rxApod_.size()
				<< " (should be " << config_.numElements << ").");
	}
}

template<typename FloatType>
void
SynthYVectorial3DTnRnProcessor<FloatType>::setTxDelays(FloatType focusX, FloatType focusY, FloatType focusZ,
							const std::vector<FloatType>& txDelays)
{
	focusX_ = focusX;
	focusY_ = focusY;
	focusZ_ = focusZ;
	txDelays_ = txDelays;
}

template<typename FloatType>
void
SynthYVectorial3DTnRnProcessor<FloatType>::prepare(unsigned int baseElement)
{
	if (baseElement + config_.numElements > config_.numElementsMux) {
		THROW_EXCEPTION(InvalidParameterException, "Invalid base element: " << baseElement << '.');
	}
	if (txDelays_.empty()) {
		THROW_EXCEPTION(InvalidStateException, "Empty tx delay list.");
	}

	baseElement_ = baseElement;
	acquisition_.prepare(baseElement_, txDelays_);
}

template<typename FloatType>
void
SynthYVectorial3DTnRnProcessor<FloatType>::getAcqData(FloatType y)
{
	AcqData data;
	data.y = y;

	acquisition_.execute(singleAcqData_);

	const std::size_t signalLength = singleAcqData_.n2() * upsamplingFactor_;
	tempSignal_.resize(signalLength);
	data.analyticSignalMatrix.resize(config_.numElements, signalLength);
	LOG_DEBUG << "signalOffset_=" << signalOffset_ << " signalLength=" << signalLength;
	if (deadZoneSamplesUp_ >= signalLength) {
		THROW_EXCEPTION(InvalidValueException,
				"Invalid dead zone: deadZoneSamplesUp_ (" << deadZoneSamplesUp_ <<
				") >= signalLength (" << signalLength << ").");
	}

	const std::size_t samplesPerChannelLow = singleAcqData_.n2();

	for (unsigned int iRxElem = 0; iRxElem < config_.numElements; ++iRxElem) {
		if (upsamplingFactor_ > 1) {
			interpolator_.interpolate(&singleAcqData_(iRxElem, 0), samplesPerChannelLow, &tempSignal_[0]);
		} else {
			auto range = singleAcqData_.range2(iRxElem);
			std::copy(range.begin(), range.end(), tempSignal_.begin());
		}

		Util::removeDC(&tempSignal_[0], tempSignal_.size(), deadZoneSamplesUp_);

		envelope_.getAnalyticSignal(&tempSignal_[0], tempSignal_.size(),
						&data.analyticSignalMatrix(iRxElem, 0));
	}

	acqData_.push_back(std::move(data));
}

template<typename FloatType>
void
SynthYVectorial3DTnRnProcessor<FloatType>::process(unsigned int firstAcq, unsigned int lastAcq,
							Matrix<XYZValueFactor<FloatType>>& gridData)
{
	LOG_DEBUG << "BEGIN ========== SynthYVectorial3DTnRnProcessor::process ==========";

	if (txDelays_.empty()) {
		THROW_EXCEPTION(InvalidStateException, "Empty tx delay list.");
	}

	Util::resetValueFactor(gridData.begin(), gridData.end());

	const std::size_t numSignals = config_.numElements;

	ThreadData threadData;
	threadData.coherenceFactor = coherenceFactor_;
	tbb::enumerable_thread_specific<ThreadData> tls(threadData);

	const FloatType fsUp = config_.samplingFrequency * upsamplingFactor_;
	const FloatType invCT = fsUp / config_.propagationSpeed;
	const std::size_t numRows = gridData.n2();
	const unsigned int centerAcq = (firstAcq + lastAcq) / 2U;
	const FloatType centerY = acqData_[centerAcq].y;
	const unsigned int numAcq = lastAcq - firstAcq + 1U;
	const FloatType valueFactor = 1.0 / numAcq;

	IterationCounter::reset(gridData.n1());

	// Delay of the first active element.
	const FloatType d0 = txDelays_[0] * fsUp;

	const XY<FloatType>& firstElem = config_.txElemPos[baseElement_];
	// Travel time between the first active element and the focus.
	const FloatType t0 = Geometry::distance3DZ0(firstElem.x, firstElem.y, focusX_, focusY_, focusZ_) * invCT;

	tbb::parallel_for(tbb::blocked_range<std::size_t>(0, gridData.n1()),
	[&, invCT, numRows](const tbb::blocked_range<std::size_t>& r) {
		auto& local = tls.local();

		local.rxSignalSumList.resize(config_.numElements * numAcq);
		local.rxDelayList.resize(config_.numElements);
		FloatType txDelay;

		// For each column:
		for (std::size_t i = r.begin(); i != r.end(); ++i) {
			// For each row:
			for (std::size_t j = 0; j < numRows; ++j) {
				XYZValueFactor<FloatType>& point = gridData(i, j);
				std::fill(local.rxSignalSumList.begin(), local.rxSignalSumList.end(), std::complex<FloatType>(0));
				for (unsigned int acq = firstAcq; acq <= lastAcq; ++acq) { // synthetic array in y
					const FloatType y = acqData_[acq].y;
					const FloatType yOffset = y - centerY;
					const unsigned int acqOffset = (acq - firstAcq) * config_.numElements;

					// Calculate the delays.
					{
						// Travel time between the focus and the point.
						const FloatType t1 = Geometry::distance3D(focusX_, focusY_, focusZ_,
												point.x, point.y - yOffset, point.z) * invCT;
						if (focusZ_ > 0) {
							txDelay = t1 + t0 + d0;
						} else {
							txDelay = t1 - (t0 - d0);
						}
					}
					for (unsigned int iRxElem = 0, end = config_.numElements; iRxElem < end; ++iRxElem) {
						const XY<FloatType>& elemPos = config_.rxElemPos[baseElement_ + iRxElem];
						local.rxDelayList[iRxElem] = Geometry::distance3DZ0(elemPos.x, elemPos.y,
												point.x, point.y - yOffset, point.z) * invCT;
					}

					for (unsigned int iRxElem = 0, rxEnd = config_.numElements; iRxElem < rxEnd; ++iRxElem) {
						// Linear interpolation.
						const FloatType delay = signalOffset_ + txDelay + local.rxDelayList[iRxElem];
						const std::size_t delayIdx = static_cast<std::size_t>(delay);
						const FloatType k = delay - delayIdx;
						if (delayIdx + 1U < acqData_[acq].analyticSignalMatrix.n2()) {
							const std::complex<FloatType>* p = &acqData_[acq].analyticSignalMatrix(iRxElem, delayIdx);
							local.rxSignalSumList[acqOffset + iRxElem] = rxApod_[iRxElem] * ((1 - k) * *p + k * *(p + 1));
						}
					}
				}
				if (local.coherenceFactor.enabled()) {
					point.factor = local.coherenceFactor.calculate(&local.rxSignalSumList[0], local.rxSignalSumList.size());
				}
				point.value = valueFactor * std::abs(std::accumulate(local.rxSignalSumList.begin(), local.rxSignalSumList.end(), std::complex<FloatType>(0)));
				point.y = centerY;
			}
		}

		IterationCounter::add(r.end() - r.begin());
	});

	std::for_each(gridData.begin(), gridData.end(),
			Util::MultiplyValueBy<XYZValueFactor<FloatType>, FloatType>(FloatType(1) / numSignals));

	LOG_DEBUG << "END ========== SynthYVectorial3DTnRnProcessor::process ==========";
}

} // namespace Lab

#endif // SYNTH_Y_SYNTH_Y_VECTORIAL_3D_TN_RN_PROCESSOR_H
