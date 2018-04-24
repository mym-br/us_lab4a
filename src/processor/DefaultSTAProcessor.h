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

#ifndef DEFAULTSTAPROCESSOR_H_
#define DEFAULTSTAPROCESSOR_H_

#include <cmath>
#include <cstddef>
#include <numeric> /* accumulate */
#include <vector>

#include <tbb/enumerable_thread_specific.h>
#include <tbb/tbb.h>

#include "CoherenceFactor.h"
#include "Exception.h"
#include "Interpolator4X.h"
#include "Log.h"
#include "Matrix2.h"
#include "Matrix3.h"
#include "STAAcquisition.h"
#include "STAConfiguration.h"
#include "STAProcessor.h"
#include "Util.h"
#include "XZValueFactor.h"

// Do not change.
#define DEFAULT_STA_PROCESSOR_UPSAMPLING_FACTOR 4



namespace Lab {

// x = 0 is at the center of the element group.
// y = 0
// z = 0 is at the surface of the array.
template<typename FloatType>
class DefaultSTAProcessor : public STAProcessor<FloatType> {
public:
	DefaultSTAProcessor(
			const STAConfiguration<FloatType>& config,
			STAAcquisition<FloatType>& acquisition,
			CoherenceFactorProcessor<FloatType>& coherenceFactor,
			FloatType peakOffset);
	virtual ~DefaultSTAProcessor() {}

	virtual void process(unsigned int baseElement, Matrix2<XZValueFactor<FloatType>>& gridData);

private:
	struct ThreadData {
		CoherenceFactorProcessor<FloatType> coherenceFactor;
		std::vector<FloatType> rxSignalSumList;
		std::vector<FloatType> delayList;
	};

	class ProcessColumn;

	DefaultSTAProcessor(const DefaultSTAProcessor&);
	DefaultSTAProcessor& operator=(const DefaultSTAProcessor&);

	const STAConfiguration<FloatType>& config_;
	unsigned int deadZoneSamplesUp_;
	STAAcquisition<FloatType>& acquisition_;
	CoherenceFactorProcessor<FloatType>& coherenceFactor_;
	Matrix3<FloatType> signalMatrix_;
	typename STAAcquisition<FloatType>::AcquisitionDataType acqData_;
	std::vector<FloatType> signal_;
	FloatType signalOffset_;
	Interpolator4X<FloatType> interpolator_;
	//LinearInterpolator interpolator_;
};



template<typename FloatType>
DefaultSTAProcessor<FloatType>::DefaultSTAProcessor(
			const STAConfiguration<FloatType>& config,
			STAAcquisition<FloatType>& acquisition,
			CoherenceFactorProcessor<FloatType>& coherenceFactor,
			FloatType peakOffset)
		: config_(config)
		, deadZoneSamplesUp_((DEFAULT_STA_PROCESSOR_UPSAMPLING_FACTOR * config.samplingFrequency) * 2.0 * config.deadZoneM / config.propagationSpeed)
		, acquisition_(acquisition)
		, coherenceFactor_(coherenceFactor)
		//, interpolator_(DEFAULT_STA_PROCESSOR_UPSAMPLING_FACTOR)
{
	signalOffset_ = (config_.samplingFrequency * DEFAULT_STA_PROCESSOR_UPSAMPLING_FACTOR) * peakOffset / config_.centerFrequency;
	LOG_DEBUG << "signalOffset_=" << signalOffset_;

	// This acquisition is done to obtain the ascan length.
	acquisition_.execute(0, 0, acqData_);
	const std::size_t signalLength = acqData_.n2() * DEFAULT_STA_PROCESSOR_UPSAMPLING_FACTOR;

	signal_.resize(signalLength);
	signalMatrix_.resize(config_.numElements, config_.numElements, signalLength);
	//signalMatrix_ = 0.0;

	if (deadZoneSamplesUp_ >= signalLength) {
		THROW_EXCEPTION(InvalidValueException, "Invalid dead zone: deadZoneSamplesUp_ (" << deadZoneSamplesUp_ <<
							") >= signalLength (" << signalLength << ").");
	}
}

template<typename FloatType>
void
DefaultSTAProcessor<FloatType>::process(unsigned int baseElement, Matrix2<XZValueFactor<FloatType>>& gridData)
{
	LOG_DEBUG << "BEGIN ========== DefaultSTAProcessor::process ==========";

	// Clear the values in gridData.
	for (auto iter = gridData.begin(); iter != gridData.end(); ++iter) {
		iter->value = 0.0;
		iter->factor = 1.0;
	}

	// Prepare the signal matrix.
	for (std::size_t txElem = 0; txElem < config_.numElements; ++txElem) {
		LOG_INFO << "ACQ txElem: " << txElem << " < " << config_.numElements;

		acquisition_.execute(baseElement, txElem, acqData_);

		const std::size_t samplesPerChannelLow = acqData_.n2();
		//const std::size_t samplesPerChannel = samplesPerChannelLow * STA_PROCESSOR_UPSAMPLING_FACTOR;

		for (std::size_t rxElem = 0; rxElem < config_.numElements; ++rxElem) {

			// Interpolate the signal.
			interpolator_.interpolate(&acqData_(rxElem, 0), samplesPerChannelLow, &signal_[0]);

			// Copy the signal to the signal matrix.
			typename Matrix3<FloatType>::Dim3Interval interval = signalMatrix_.dim3Interval(txElem, rxElem);
			std::copy(signal_.begin(), signal_.end(), interval.first);

			Util::removeDC(&signalMatrix_(txElem, rxElem, 0), signalMatrix_.n3(), deadZoneSamplesUp_);
		}
	}

	ThreadData threadData;
	threadData.coherenceFactor = coherenceFactor_;
	tbb::enumerable_thread_specific<ThreadData> tls{threadData};

	std::vector<FloatType> xArray(config_.numElements);
	const FloatType xArrayCenter = (config_.numElementsMux - 1) * config_.pitch / 2.0;
	for (unsigned int elem = 0; elem < config_.numElements; ++elem) {
		xArray[elem] = (baseElement + elem) * config_.pitch - xArrayCenter;
	}

	const FloatType invCT = (config_.samplingFrequency * DEFAULT_STA_PROCESSOR_UPSAMPLING_FACTOR) / config_.propagationSpeed;
	const std::size_t numRows = gridData.n2();

	tbb::parallel_for(tbb::blocked_range<std::size_t>(0, gridData.n1()),
	[&, invCT, numRows](const tbb::blocked_range<std::size_t>& r) {
		LOG_DEBUG << "IMG col range start = " << r.begin() << " n = " << (r.end() - r.begin());

		auto& local = tls.local();

		local.rxSignalSumList.resize(config_.numElements);
		local.delayList.resize(config_.numElements);

		// For each column:
		for (std::size_t i = r.begin(); i != r.end(); ++i) {
			// For each row:
			for (std::size_t j = 0; j < numRows; ++j) {

				std::fill(local.rxSignalSumList.begin(), local.rxSignalSumList.end(), FloatType{0});
				XZValueFactor<FloatType>& point = gridData(i, j);

				// Calculate the delays.
				for (unsigned int elem = 0; elem < config_.numElements; ++elem) {
					const FloatType dx = point.x - xArray[elem];
					const FloatType dz = point.z /* - zArray*/; // zArray = 0
					local.delayList[elem] = std::sqrt(dx * dx + dz * dz) * invCT;
				}

				for (std::size_t txElem = 0; txElem < config_.numElements; ++txElem) {
					const FloatType txDelay = local.delayList[txElem];
					for (std::size_t rxElem = 0; rxElem < config_.numElements; ++rxElem) {
#if 0
						// Nearest neighbor.
						const std::size_t delayIdx = static_cast<std::size_t>(FloatType{0.5} + signalOffset_ + txDelay + local.delayList[rxElem]);
						if (delayIdx < signalMatrix_.n3()) {
							local.rxSignalSumList[rxElem] += signalMatrix_(txElem, rxElem, delayIdx);
						}
#else
						// Linear interpolation.
						const FloatType delay = signalOffset_ + txDelay + local.delayList[rxElem];
						const std::size_t delayIdx = static_cast<std::size_t>(delay);
						const FloatType k = delay - delayIdx;
						if (delayIdx < signalMatrix_.n3() - 1) {
							const FloatType* p = &signalMatrix_(txElem, rxElem, delayIdx);
							local.rxSignalSumList[rxElem] += (1 - k) * *p + k * *(p + 1);
						}
#endif
					}
				}

				if (local.coherenceFactor.enabled()) {
					point.factor = local.coherenceFactor.calculate(&local.rxSignalSumList[0], local.rxSignalSumList.size());
				}
				point.value = std::accumulate(local.rxSignalSumList.begin(), local.rxSignalSumList.end(), FloatType{0});
			}
		}
	});

	LOG_DEBUG << "END ========== DefaultSTAProcessor::process ==========";
}

} // namespace Lab

#endif /* DEFAULTSTAPROCESSOR_H_ */
