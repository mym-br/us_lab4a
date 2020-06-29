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

#ifndef VECTORIALTWOMEDIUMSTAPROCESSOR_H_
#define VECTORIALTWOMEDIUMSTAPROCESSOR_H_

#include <algorithm> /* fill, min */
#include <cmath> /* abs, sqrt */
#include <complex>
#include <cstddef> /* std::size_t */
#include <numeric> /* accumulate */
#include <vector>

#include <tbb/blocked_range.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>

#include "CoherenceFactor.h"
#include "Exception.h"
#include "FermatPrinciple.h"
#include "HilbertEnvelope.h"
#include "Interpolator.h"
#include "Log.h"
#include "Matrix.h"
#include "STAAcquisition.h"
#include "Tensor3.h"
#include "TwoMediumSTAConfiguration.h"
#include "Util.h"
#include "XZ.h"
#include "XZValueFactor.h"

//TODO: remove test (timer)
#include "Timer.h"

// Depends on the signal.
// 1.0 --> pi radian / sample at the original sampling rate.
#define VECTORIAL_TWO_MEDIUM_STA_PROCESSOR_UPSAMP_FILTER_HALF_TRANSITION_WIDTH 0.2



namespace Lab {

template<typename TFloat>
class VectorialTwoMediumSTAProcessor {
public:
	VectorialTwoMediumSTAProcessor(
			const TwoMediumSTAConfiguration<TFloat>& config,
			STAAcquisition<TFloat>& acquisition,
			unsigned int upsamplingFactor,
			AnalyticSignalCoherenceFactorProcessor<TFloat>& coherenceFactor,
			TFloat maxFermatBlockSize, TFloat peakOffset);
	~VectorialTwoMediumSTAProcessor() = default;

	void prepare(unsigned int baseElement);
	void process(const std::vector<XZ<TFloat>>& interfacePointList, Matrix<XZValueFactor<TFloat>>& gridData);

private:
	struct ThreadData {
		AnalyticSignalCoherenceFactorProcessor<TFloat> coherenceFactor;
		std::vector<std::complex<TFloat>> rxSignalSumList;
		std::vector<TFloat> delayList;
	};

	class ProcessColumn;

	VectorialTwoMediumSTAProcessor(const VectorialTwoMediumSTAProcessor&) = delete;
	VectorialTwoMediumSTAProcessor& operator=(const VectorialTwoMediumSTAProcessor&) = delete;
	VectorialTwoMediumSTAProcessor(VectorialTwoMediumSTAProcessor&&) = delete;
	VectorialTwoMediumSTAProcessor& operator=(VectorialTwoMediumSTAProcessor&&) = delete;

	const TwoMediumSTAConfiguration<TFloat>& config_;
	unsigned int deadZoneSamplesUp_;
	STAAcquisition<TFloat>& acquisition_;
	unsigned int upsamplingFactor_;
	AnalyticSignalCoherenceFactorProcessor<TFloat>& coherenceFactor_;
	TFloat maxFermatBlockSize_;
	const TFloat minLambda_;
	Tensor3<std::complex<TFloat>> analyticSignalMatrix_;
	typename STAAcquisition<TFloat>::AcquisitionDataType acqData_;
	std::vector<TFloat> signal_;
	TFloat signalOffset_;
	Interpolator<TFloat> interpolator_;
	HilbertEnvelope<TFloat> envelope_;
};



template<typename TFloat>
VectorialTwoMediumSTAProcessor<TFloat>::VectorialTwoMediumSTAProcessor(
			const TwoMediumSTAConfiguration<TFloat>& config,
			STAAcquisition<TFloat>& acquisition,
			unsigned int upsamplingFactor,
			AnalyticSignalCoherenceFactorProcessor<TFloat>& coherenceFactor,
			TFloat maxFermatBlockSize, TFloat peakOffset)
		: config_(config)
		, deadZoneSamplesUp_((upsamplingFactor * config.samplingFrequency) * 2.0 * config.deadZoneM / config.propagationSpeed1)
		, acquisition_(acquisition)
		, upsamplingFactor_(upsamplingFactor)
		, coherenceFactor_(coherenceFactor)
		, maxFermatBlockSize_(maxFermatBlockSize)
		, minLambda_(std::min(config_.propagationSpeed1, config_.propagationSpeed2) / config_.centerFrequency)
{
	signalOffset_ = (config_.samplingFrequency * upsamplingFactor_) * peakOffset / config_.centerFrequency;
	LOG_DEBUG << "signalOffset_=" << signalOffset_;

	// This acquisition is done to obtain the signal length.
	acquisition_.prepare(0);
	acquisition_.execute(0, acqData_);
	const std::size_t signalLength = acqData_.n2() * upsamplingFactor_;

	signal_.resize(signalLength);
	analyticSignalMatrix_.resize(config_.numElements, config_.numElements, signalLength);

	if (upsamplingFactor_ > 1) {
		interpolator_.prepare(upsamplingFactor_, VECTORIAL_TWO_MEDIUM_STA_PROCESSOR_UPSAMP_FILTER_HALF_TRANSITION_WIDTH);
	}

	if (deadZoneSamplesUp_ >= signalLength) {
		THROW_EXCEPTION(InvalidValueException, "Invalid dead zone: deadZoneSamplesUp_ (" << deadZoneSamplesUp_ <<
							") >= signalLength (" << signalLength << ").");
	}
}

template<typename TFloat>
void
VectorialTwoMediumSTAProcessor<TFloat>::prepare(unsigned int baseElement)
{
	acquisition_.prepare(baseElement);
}

template<typename TFloat>
void
VectorialTwoMediumSTAProcessor<TFloat>::process(const std::vector<XZ<TFloat>>& interfacePointList,
							Matrix<XZValueFactor<TFloat>>& gridData)
{
	LOG_DEBUG << "BEGIN ========== VectorialTwoMediumSTAProcessor::process ==========";

	// Clear the values in gridData.
	for (typename Matrix<XZValueFactor<TFloat>>::Iterator iter = gridData.begin(); iter != gridData.end(); ++iter) {
		iter->value = 0.0;
		iter->factor = 1.0;
	}

	XZ<TFloat> p1 = interfacePointList[0];
	XZ<TFloat> p2 = interfacePointList[1];
	const TFloat dx = p2.x - p1.x;
	const TFloat dz = p2.z - p1.z;
	const TFloat r = std::sqrt(dx * dx + dz * dz);
	const unsigned int fermatBlockSize = FermatPrinciple::calcBlockSizeForTwoStepSearch(interfacePointList.size(), r, minLambda_, maxFermatBlockSize_);

	//TODO: remove test (timer)
	Timer t0;

	// Prepare the signal matrix.
	for (std::size_t txElem = 0; txElem < config_.numElements; ++txElem) {
		acquisition_.execute(txElem, acqData_);

		const std::size_t samplesPerChannelLow = acqData_.n2();
		//const std::size_t samplesPerChannel = samplesPerChannelLow * STA_PROCESSOR_UPSAMPLING_FACTOR;

		for (std::size_t rxElem = 0; rxElem < config_.numElements; ++rxElem) {

			// Interpolate the signal.
			interpolator_.interpolate(&acqData_(rxElem, 0), samplesPerChannelLow, &signal_[0]);

			Util::removeDC(&signal_[0], signal_.size(), deadZoneSamplesUp_);

			// Obtain the analytic signal.
			envelope_.getAnalyticSignal(&signal_[0], signal_.size(), &analyticSignalMatrix_(txElem, rxElem, 0));
		}

		//LOG_DEBUG << "PREP txElem = " << txElem;
	}

	//TODO: remove test (timer)
	LOG_DEBUG << "PREP time = " << t0.getTime();

	ThreadData threadData;
	threadData.coherenceFactor = coherenceFactor_;

	tbb::enumerable_thread_specific<ThreadData> processColumnTLS(threadData);
	tbb::parallel_for(
		tbb::blocked_range<std::size_t>(0, gridData.n1()),
		ProcessColumn(
			*this,
			gridData.n2(),
			config_,
			upsamplingFactor_,
			processColumnTLS,
			signalOffset_,
			analyticSignalMatrix_,
			fermatBlockSize,
			interfacePointList,
			gridData));

	LOG_DEBUG << "END ========== VectorialTwoMediumSTAProcessor::process ==========";
}



template<typename TFloat>
class VectorialTwoMediumSTAProcessor<TFloat>::ProcessColumn {
public:
	ProcessColumn(
			VectorialTwoMediumSTAProcessor<TFloat>& p,
			std::size_t numRows,
			const TwoMediumSTAConfiguration<TFloat>& config,
			unsigned int upsamplingFactor,
			tbb::enumerable_thread_specific<ThreadData>& processColumnTLS,
			TFloat signalOffset,
			const Tensor3<std::complex<TFloat>>& analyticSignalMatrix,
			unsigned int fermatBlockSize,
			const std::vector<XZ<TFloat>>& interfacePointList,
			Matrix<XZValueFactor<TFloat>>& gridData)
				: p_(p)
				, numRows_(numRows)
				, config_(config)
				, processColumnTLS_(processColumnTLS)
				, signalOffset_(signalOffset)
				, analyticSignalMatrix_(analyticSignalMatrix)
				, fs_(config_.samplingFrequency * upsamplingFactor)
				, invC1_(1.0 / config_.propagationSpeed1)
				, invC1T_(fs_ * invC1_)
				, invC2_(1.0 / config_.propagationSpeed2)
				, fermatBlockSize_(fermatBlockSize)
				, interfacePointList_(interfacePointList)
				, xArray_(config_.numElements)
				, gridData_(gridData) {

		const TFloat xArrayCenter = (config_.numElements - 1) * config_.pitch / 2.0;
		for (unsigned int elem = 0; elem < config_.numElements; ++elem) {
			xArray_[elem] = elem * config_.pitch - xArrayCenter;
		}
	}

	void operator()(const tbb::blocked_range<std::size_t>& r) const {
		LOG_DEBUG << "col = " << r.begin() << " n = " << (r.end() - r.begin());

		ThreadData& data = processColumnTLS_.local();

		data.rxSignalSumList.resize(config_.numElements);
		data.delayList.resize(config_.numElements);

		for (std::size_t i = r.begin(); i != r.end(); ++i) {
			for (unsigned int j = 0; j < numRows_; ++j) {

				std::fill(data.rxSignalSumList.begin(), data.rxSignalSumList.end(), std::complex<TFloat>(0));
				XZValueFactor<TFloat>& point = gridData_(i, j);

				// Find the z coordinate of the interface.
				TFloat zIdxMin;
				unsigned int idxMin;
				FermatPrinciple::findNearestPointInXInTwoSteps(
						fermatBlockSize_,
						interfacePointList_,
						point.x,
						zIdxMin, idxMin);
				const TFloat zInterface = zIdxMin;
				const TFloat dz = point.z - zInterface;

				if (dz > 0.0) { // the point is above the interface

					// Calculate the delays.
					for (unsigned int elem = 0; elem < config_.numElements; ++elem) {
						// Fermat's principle. Find the fastest path.
						TFloat tMin;
						unsigned int idxMin;
						FermatPrinciple::findMinTimeInTwoSteps(
								fermatBlockSize_,
								p_.config_.propagationSpeed1, p_.config_.propagationSpeed2,
								interfacePointList_,
								xArray_[elem], TFloat(0), point.x, point.z,
								tMin, idxMin);

						if (idxMin == 0 || idxMin == interfacePointList_.size() - 1) {
							// Ignore the points for which the interface is too short.
							data.delayList[elem] = analyticSignalMatrix_.n3(); // mark as invalid
						} else {
							data.delayList[elem] = tMin * fs_;
						}
					}

					for (unsigned int txElem = 0; txElem < config_.numElements; ++txElem) {
						const TFloat txDelay = data.delayList[txElem];
						for (unsigned int rxElem = 0; rxElem < config_.numElements; ++rxElem) {
							// Linear interpolation.
							const TFloat delay = signalOffset_ + txDelay + data.delayList[rxElem];
							const std::size_t delayIdx = static_cast<std::size_t>(delay);
							const TFloat k = delay - delayIdx;
							if (delayIdx < analyticSignalMatrix_.n3() - 1) {
								const std::complex<TFloat>* p = &analyticSignalMatrix_(txElem, rxElem, delayIdx);
								data.rxSignalSumList[rxElem] += (1 - k) * *p + k * *(p + 1);
							}
						}
					}

					if (data.coherenceFactor.enabled()) {
						point.factor = data.coherenceFactor.calculate(&data.rxSignalSumList[0], data.rxSignalSumList.size());
					}
					point.value = std::abs(std::accumulate(data.rxSignalSumList.begin(), data.rxSignalSumList.end(), std::complex<TFloat>(0)));
				}
			}
		}
	}
private:
	VectorialTwoMediumSTAProcessor<TFloat>& p_;
	const std::size_t numRows_;
	const TwoMediumSTAConfiguration<TFloat>& config_;
	tbb::enumerable_thread_specific<ThreadData>& processColumnTLS_;
	const TFloat signalOffset_;
	const Tensor3<std::complex<TFloat>>& analyticSignalMatrix_;
	const TFloat fs_;
	const TFloat invC1_;
	const TFloat invC1T_;
	const TFloat invC2_;
	const unsigned int fermatBlockSize_;
	const std::vector<XZ<TFloat>>& interfacePointList_;
	std::vector<TFloat> xArray_;
	Matrix<XZValueFactor<TFloat>>& gridData_;
};

} // namespace Lab

#endif /* VECTORIALTWOMEDIUMSTAPROCESSOR_H_ */
