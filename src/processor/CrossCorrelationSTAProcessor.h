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

#ifndef CROSSCORRELATIONSTAPROCESSOR_H_
#define CROSSCORRELATIONSTAPROCESSOR_H_

#include <algorithm> /* copy, fill, reverse_copy */
#include <cmath> /* sqrt */
#include <cstddef> /* std::size_t */
#include <numeric> /* accumulate */
#include <vector>

#include <tbb/enumerable_thread_specific.h>
#include <tbb/tbb.h>

#include "ArrayGeometry.h"
#include "CoherenceFactor.h"
#include "Exception.h"
#include "FFTWFilter.h"
#include "Interpolator.h"
#include "Log.h"
#include "Matrix2.h"
#include "Matrix3.h"
#include "STAAcquisition.h"
#include "STAConfiguration.h"
#include "STAProcessor.h"
#include "Util.h"
#include "XZValueFactor.h"

// Depends on the signal.
// 1.0 --> pi radian / sample at the original sampling rate.
#define CROSS_CORRELATION_STA_PROCESSOR_UPSAMP_FILTER_HALF_TRANSITION_WIDTH 0.2



namespace Lab {

// x = 0 is at the center of the element group.
// y = 0
// z = 0 is at the surface of the array.
template<typename FloatType>
class CrossCorrelationSTAProcessor : public STAProcessor<FloatType> {
public:
	CrossCorrelationSTAProcessor(
			const STAConfiguration<FloatType>& config,
			STAAcquisition<FloatType>& acquisition,
			unsigned int upsamplingFactor,
			CoherenceFactorProcessor<FloatType>& coherenceFactor,
			const std::vector<FloatType> &refPulse);
	virtual ~CrossCorrelationSTAProcessor() {}

	virtual void process(unsigned int baseElement, Matrix2<XZValueFactor<FloatType> >& gridData);

private:
	struct ProcessColumn;
	struct ThreadData {
		CoherenceFactorProcessor<FloatType> coherenceFactor;
		std::vector<FloatType> rxSignalSumList;
		std::vector<FloatType> delayList;
	};

	CrossCorrelationSTAProcessor(const CrossCorrelationSTAProcessor&);
	CrossCorrelationSTAProcessor& operator=(const CrossCorrelationSTAProcessor&);

	const STAConfiguration<FloatType>& config_;
	unsigned int deadZoneSamplesUp_;
	STAAcquisition<FloatType>& acquisition_;
	unsigned int upsamplingFactor_;
	CoherenceFactorProcessor<FloatType>& coherenceFactor_;
	bool getEnvelope_;
	Matrix3<FloatType> signalMatrix_;
	typename STAAcquisition<FloatType>::AcquisitionDataType acqData_;
	std::vector<FloatType> tempSignal1_;
	std::vector<FloatType> tempSignal2_;
	unsigned int signalOffset_;
	Interpolator<FloatType> interpolator_;
	FFTWFilter<FloatType> revRefPulseFilter_;
};



template<typename FloatType>
CrossCorrelationSTAProcessor<FloatType>::CrossCorrelationSTAProcessor(
			const STAConfiguration<FloatType>& config,
			STAAcquisition<FloatType>& acquisition,
			unsigned int upsamplingFactor,
			CoherenceFactorProcessor<FloatType>& coherenceFactor,
			const std::vector<FloatType>& refPulse)
		: config_(config)
		, deadZoneSamplesUp_((upsamplingFactor * config.samplingFrequency) * 2.0 * config.deadZoneM / config.propagationSpeed)
		, acquisition_(acquisition)
		, upsamplingFactor_(upsamplingFactor)
		, coherenceFactor_(coherenceFactor)
{
	if (refPulse.size() == 0) {
		THROW_EXCEPTION(InvalidParameterException, "The vector refPulse is empty.");
	}

	// This acquisition is done to obtain the ascan length.
	acquisition_.execute(0, 0, acqData_);
	const std::size_t origSignalLength = acqData_.n2();

	tempSignal1_.resize(origSignalLength * upsamplingFactor_);

	// The original signal is upsampled and convolved with the reverse of the upsampled refPulse.
	const std::size_t signalLength = (origSignalLength + refPulse.size()) * upsamplingFactor_ - 1;
	signalMatrix_.resize(config_.numElements, config_.numElements, signalLength);

	std::vector<FloatType> revRefPulse(refPulse.size() * upsamplingFactor_);
	if (upsamplingFactor_ > 1) {
		interpolator_.prepare(upsamplingFactor_, CROSS_CORRELATION_STA_PROCESSOR_UPSAMP_FILTER_HALF_TRANSITION_WIDTH);

		interpolator_.interpolate(&refPulse[0], refPulse.size(), &revRefPulse[0]);
		std::reverse(revRefPulse.begin(), revRefPulse.end());
		revRefPulseFilter_.setCoefficients(revRefPulse);
	} else {
		std::reverse_copy(refPulse.begin(), refPulse.end(), revRefPulse.begin());
		revRefPulseFilter_.setCoefficients(revRefPulse);
	}

	signalOffset_ = revRefPulse.size() - 1; // cross-correlation using convolution (revRefPulseFilter_)
	LOG_DEBUG << "signalOffset_: " << signalOffset_ << " origSignalLength: " << origSignalLength << " signalLength: " << signalLength;

	if (deadZoneSamplesUp_ >= signalLength) {
		THROW_EXCEPTION(InvalidValueException, "Invalid dead zone: deadZoneSamplesUp_ (" << deadZoneSamplesUp_ <<
							") >= signalLength (" << signalLength << ").");
	}
}

template<typename FloatType>
void
CrossCorrelationSTAProcessor<FloatType>::process(unsigned int baseElement, Matrix2<XZValueFactor<FloatType> >& gridData)
{
	LOG_DEBUG << "BEGIN ========== CrossCorrelationSTAProcessor::process ==========";

	Util::resetValueFactor(gridData.begin(), gridData.end());

	// Prepares the signal matrix.
	for (unsigned int txElem = 0; txElem < config_.numElements; ++txElem) {
		acquisition_.execute(baseElement, txElem, acqData_);

		const std::size_t samplesPerChannelLow = acqData_.n2();

		for (unsigned int rxElem = 0; rxElem < config_.numElements; ++rxElem) {
			if (upsamplingFactor_ > 1) {
				// Interpolates the signal.
				interpolator_.interpolate(&acqData_(rxElem, 0), samplesPerChannelLow, &tempSignal1_[0]);
			} else {
				typename Matrix2<FloatType>::Dim2Interval interval = acqData_.dim2Interval(rxElem);
				std::copy(interval.first, interval.second, tempSignal1_.begin());
			}

			Util::removeDC(&tempSignal1_[0], tempSignal1_.size(), deadZoneSamplesUp_);

			// Cross-correlation using convolution.
			revRefPulseFilter_.filter(tempSignal1_, tempSignal2_);

			std::copy(tempSignal2_.begin(), tempSignal2_.end(), &signalMatrix_(txElem, rxElem, 0));
		}

		LOG_DEBUG << "PRC PREP txElem = " << txElem;
	}

	std::vector<FloatType> xArray;
	ArrayGeometry::getElementXCentered2D(config_.numElements, config_.pitch, xArray);

	ThreadData threadData;
	threadData.coherenceFactor = coherenceFactor_;
	tbb::enumerable_thread_specific<ThreadData> processColumnTLS(threadData);

	ProcessColumn op = {
		gridData.n2(),
		config_,
		signalOffset_,
		signalMatrix_,
		(config_.samplingFrequency * upsamplingFactor_) / config_.propagationSpeed,
		xArray,
		processColumnTLS,
		gridData
	};
	tbb::parallel_for(tbb::blocked_range<std::size_t>(0, gridData.n1()), op);

	LOG_DEBUG << "END ========== CrossCorrelationSTAProcessor::process ==========";
}

template<typename FloatType>
struct CrossCorrelationSTAProcessor<FloatType>::ProcessColumn {
	void operator()(const tbb::blocked_range<std::size_t>& r) const {
		LOG_DEBUG << "PRC col = " << r.begin() << " n = " << (r.end() - r.begin());

		ThreadData& data = processColumnTLS.local();

		data.rxSignalSumList.resize(config.numElements);
		data.delayList.resize(config.numElements);

		// For each point in x-z:
		for (std::size_t i = r.begin(); i != r.end(); ++i) {
			for (unsigned int j = 0; j < numRows; ++j) {

				std::fill(data.rxSignalSumList.begin(), data.rxSignalSumList.end(), FloatType(0));
				XZValueFactor<FloatType>& point = gridData(i, j);

				// Calculates the delays.
				for (unsigned int elem = 0; elem < config.numElements; ++elem) {
					const FloatType dx = point.x - xArray[elem];
					const FloatType dz = point.z /* - zArray*/; // zArray = 0
					data.delayList[elem] = std::sqrt(dx * dx + dz * dz) * invCT;
				}

				for (unsigned int txElem = 0; txElem < config.numElements; ++txElem) {
					const FloatType txDelay = data.delayList[txElem];
					for (unsigned int rxElem = 0; rxElem < config.numElements; ++rxElem) {
						// Linear interpolation.
						const FloatType delay = signalOffset + (txDelay + data.delayList[rxElem]);
						const std::size_t delayIdx = static_cast<std::size_t>(delay);
						const FloatType k = delay - delayIdx;
						if (delayIdx < signalMatrix.n3() - 1) {
							const FloatType* p = &signalMatrix(txElem, rxElem, delayIdx);
							data.rxSignalSumList[rxElem] += (1 - k) * *p + k * *(p + 1);
						}
					}
				}

				if (data.coherenceFactor.enabled()) {
					point.factor = data.coherenceFactor.calculate(&data.rxSignalSumList[0], data.rxSignalSumList.size());
				}
				point.value = std::accumulate(data.rxSignalSumList.begin(), data.rxSignalSumList.end(), FloatType(0));
			}
		}
	}

	// Use only small types or references/pointers, because this object may be copied many times.
	const std::size_t numRows;
	const STAConfiguration<FloatType>& config;
	const unsigned int signalOffset;
	const Matrix3<FloatType>& signalMatrix;
	const FloatType invCT;
	const std::vector<FloatType>& xArray;
	tbb::enumerable_thread_specific<ThreadData>& processColumnTLS;
	Matrix2<XZValueFactor<FloatType>>& gridData;
};

} // namespace Lab

#endif // CROSSCORRELATIONSTAPROCESSOR_H_
