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

#include <cmath>
#include <complex>
#include <cstddef> /* std::size_t */
#include <numeric> /* accumulate */
#include <vector>

#include <tbb/enumerable_thread_specific.h>
#include <tbb/tbb.h>

#include "ArrayGeometry.h"
#include "CoherenceFactor.h"
#include "Exception.h"
#include "HilbertEnvelope.h"
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
#define VECTORIAL_STA_PROCESSOR_UPSAMP_FILTER_HALF_TRANSITION_WIDTH 0.2



namespace Lab {

// x = 0 is at the center of the element group.
// y = 0
// z = 0 is at the surface of the array.
template<typename FloatType>
class VectorialSTAProcessor : public STAProcessor<FloatType> {
public:
	struct ProcessColumn;
	struct ThreadData {
		AnalyticSignalCoherenceFactorProcessor<FloatType> coherenceFactor;
		std::vector<std::complex<FloatType>> rxSignalSumList;
		std::vector<FloatType> delayList;
	};

	VectorialSTAProcessor(
			const STAConfiguration<FloatType>& config,
			STAAcquisition<FloatType>& acquisition,
			unsigned int upsamplingFactor,
			AnalyticSignalCoherenceFactorProcessor<FloatType>& coherenceFactor,
			FloatType peakOffset,
			bool calculateEnvelope);
	virtual ~VectorialSTAProcessor() {}

	virtual void process(unsigned int baseElement, Matrix2<XZValueFactor<FloatType>>& gridData);

private:
	VectorialSTAProcessor(const VectorialSTAProcessor&);
	VectorialSTAProcessor& operator=(const VectorialSTAProcessor&);

	const STAConfiguration<FloatType>& config_;
	unsigned int deadZoneSamplesUp_;
	STAAcquisition<FloatType>& acquisition_;
	unsigned int upsamplingFactor_;
	AnalyticSignalCoherenceFactorProcessor<FloatType>& coherenceFactor_;
	Matrix3<std::complex<FloatType>> analyticSignalMatrix_;
	typename STAAcquisition<FloatType>::AcquisitionDataType acqData_;
	std::vector<FloatType> tempSignal_;
	FloatType signalOffset_;
	Interpolator<FloatType> interpolator_;
	HilbertEnvelope<FloatType> envelope_;
	bool calculateEnvelope_;
};



template<typename FloatType>
VectorialSTAProcessor<FloatType>::VectorialSTAProcessor(
			const STAConfiguration<FloatType>& config,
			STAAcquisition<FloatType>& acquisition,
			unsigned int upsamplingFactor,
			AnalyticSignalCoherenceFactorProcessor<FloatType>& coherenceFactor,
			FloatType peakOffset,
			bool calculateEnvelope)
		: config_(config)
		, deadZoneSamplesUp_((upsamplingFactor * config.samplingFrequency) * 2.0 * config.deadZoneM / config.propagationSpeed)
		, acquisition_(acquisition)
		, upsamplingFactor_(upsamplingFactor)
		, coherenceFactor_(coherenceFactor)
		, calculateEnvelope_(calculateEnvelope)
{
	// This acquisition is done to obtain the ascan length.
	acquisition_.execute(0, 0, acqData_);
	const std::size_t signalLength = acqData_.n2() * upsamplingFactor_;

	tempSignal_.resize(signalLength);
	analyticSignalMatrix_.resize(config_.numElements, config_.numElements, signalLength);

	if (upsamplingFactor_ > 1) {
		interpolator_.prepare(upsamplingFactor_, VECTORIAL_STA_PROCESSOR_UPSAMP_FILTER_HALF_TRANSITION_WIDTH);
	}

	signalOffset_ = (config_.samplingFrequency * upsamplingFactor_) * peakOffset / config_.centerFrequency;
	LOG_DEBUG << "signalOffset_=" << signalOffset_ << " signalLength: " << signalLength;

	if (deadZoneSamplesUp_ >= signalLength) {
		THROW_EXCEPTION(InvalidValueException, "Invalid dead zone: deadZoneSamplesUp_ (" << deadZoneSamplesUp_ <<
							") >= signalLength (" << signalLength << ").");
	}
}

template<typename FloatType>
void
VectorialSTAProcessor<FloatType>::process(unsigned int baseElement, Matrix2<XZValueFactor<FloatType> >& gridData)
{
	LOG_DEBUG << "BEGIN ========== VectorialSTAProcessor::process ==========";

	Util::resetValueFactor(gridData.begin(), gridData.end());

	// Prepares the signal matrix.
	for (unsigned int txElem = 0; txElem < config_.numElements; ++txElem) {
		acquisition_.execute(baseElement, txElem, acqData_);

		const std::size_t samplesPerChannelLow = acqData_.n2();

		for (unsigned int rxElem = 0; rxElem < config_.numElements; ++rxElem) {
			if (upsamplingFactor_ > 1) {
				// Interpolates the signal.
				interpolator_.interpolate(&acqData_(rxElem, 0), samplesPerChannelLow, &tempSignal_[0]);
			} else {
				typename Matrix2<FloatType>::Dim2Interval interval = acqData_.dim2Interval(rxElem);
				std::copy(interval.first, interval.second, tempSignal_.begin());
			}

			Util::removeDC(&tempSignal_[0], tempSignal_.size(), deadZoneSamplesUp_);

			// Obtains the analytic signal.
			envelope_.getAnalyticSignal(&tempSignal_[0], tempSignal_.size(), &analyticSignalMatrix_(txElem, rxElem, 0));
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
		analyticSignalMatrix_,
		(config_.samplingFrequency * upsamplingFactor_) / config_.propagationSpeed,
		xArray,
		calculateEnvelope_,
		processColumnTLS,
		gridData
	};
	tbb::parallel_for(tbb::blocked_range<std::size_t>(0, gridData.n1()), op);

	LOG_DEBUG << "END ========== VectorialSTAProcessor::process ==========";
}



template<typename FloatType>
struct VectorialSTAProcessor<FloatType>::ProcessColumn {
	void operator()(const tbb::blocked_range<std::size_t>& r) const {
		LOG_DEBUG << "PRC col = " << r.begin() << " n = " << (r.end() - r.begin());

		ThreadData& data = processColumnTLS.local();

		data.rxSignalSumList.resize(config.numElements);
		data.delayList.resize(config.numElements);

		// For each point in x-z:
		for (std::size_t i = r.begin(); i != r.end(); ++i) {
			for (unsigned int j = 0; j < numRows; ++j) {

				std::fill(data.rxSignalSumList.begin(), data.rxSignalSumList.end(), std::complex<FloatType>(0));
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
						const FloatType delay = signalOffset + txDelay + data.delayList[rxElem];
						const std::size_t delayIdx = static_cast<std::size_t>(delay);
						const FloatType k = delay - delayIdx;
						if (delayIdx < analyticSignalMatrix.n3() - 1) {
							const std::complex<FloatType>* p = &analyticSignalMatrix(txElem, rxElem, delayIdx);
							data.rxSignalSumList[rxElem] += (1 - k) * *p + k * *(p + 1);
						}
					}
				}

				if (data.coherenceFactor.enabled()) {
					point.factor = data.coherenceFactor.calculate(&data.rxSignalSumList[0], data.rxSignalSumList.size());
				}
				if (calculateEnvelope) {
					point.value = std::abs(std::accumulate(data.rxSignalSumList.begin(), data.rxSignalSumList.end(), std::complex<FloatType>(0)));
				} else {
					point.value = std::accumulate(data.rxSignalSumList.begin(), data.rxSignalSumList.end(), std::complex<FloatType>(0)).real();
				}
			}
		}
	}

	// Use only small types or references/pointers, because this object may be copied many times.
	const std::size_t numRows;
	const STAConfiguration<FloatType>& config;
	const FloatType signalOffset;
	const Matrix3<std::complex<FloatType>>& analyticSignalMatrix;
	const FloatType invCT;
	const std::vector<FloatType>& xArray;
	const bool calculateEnvelope;
	tbb::enumerable_thread_specific<ThreadData>& processColumnTLS;
	Matrix2<XZValueFactor<FloatType>>& gridData;
};

} // namespace Lab

#endif /* VECTORIALSTAPROCESSOR_H_ */
