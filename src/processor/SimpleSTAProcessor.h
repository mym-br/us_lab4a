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

#ifndef SIMPLESTAPROCESSOR_H_
#define SIMPLESTAPROCESSOR_H_

#include <algorithm> /* for_each */
#include <cmath>
#include <cstddef> /* std::size_t */
#include <vector>

#include <tbb/enumerable_thread_specific.h>
#include <tbb/tbb.h>

#include "ArrayGeometry.h"
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
#define SIMPLE_STA_PROCESSOR_UPSAMPLING_FACTOR 4



namespace Lab {

// This processor does not support coherence factors.
// x = 0 is at the center of the element group.
// y = 0
// z = 0 is at the surface of the array.
template<typename FloatType>
class SimpleSTAProcessor : public STAProcessor<FloatType> {
public:
	SimpleSTAProcessor(
			const STAConfiguration<FloatType>& config,
			STAAcquisition<FloatType>& acquisition,
			FloatType peakOffset);
	virtual ~SimpleSTAProcessor() {}

	virtual void process(unsigned int baseElement, Matrix2<XZValueFactor<FloatType>>& gridData);

private:
	struct ThreadData {
		std::vector<FloatType> delayList;
	};

	SimpleSTAProcessor(const SimpleSTAProcessor&) = delete;
	SimpleSTAProcessor& operator=(const SimpleSTAProcessor&) = delete;

	const STAConfiguration<FloatType>& config_;
	STAAcquisition<FloatType>& acquisition_;
	Matrix3<FloatType> signalMatrix_;
	typename STAAcquisition<FloatType>::AcquisitionDataType acqData_;
	std::vector<FloatType> tempSignal_;
	FloatType signalOffset_;
	Interpolator4X<FloatType> interpolator_;
};



template<typename FloatType>
SimpleSTAProcessor<FloatType>::SimpleSTAProcessor(
			const STAConfiguration<FloatType>& config,
			STAAcquisition<FloatType>& acquisition,
			FloatType peakOffset)
		: config_(config)
		, acquisition_(acquisition)
{
	// This acquisition is done to obtain the signal length.
	acquisition_.execute(0, 0, acqData_);
	const std::size_t signalLength = acqData_.n2() * SIMPLE_STA_PROCESSOR_UPSAMPLING_FACTOR;

	tempSignal_.resize(signalLength);
	signalMatrix_.resize(config_.lastTxElem - config_.firstTxElem + 1, config_.numElements, signalLength);

	signalOffset_ = (config_.samplingFrequency * SIMPLE_STA_PROCESSOR_UPSAMPLING_FACTOR) * peakOffset / config_.centerFrequency;
	LOG_DEBUG << "signalOffset_=" << signalOffset_ << " signalLength=" << signalLength;
}

template<typename FloatType>
void
SimpleSTAProcessor<FloatType>::process(unsigned int baseElement, Matrix2<XZValueFactor<FloatType>>& gridData)
{
	LOG_DEBUG << "BEGIN ========== SimpleSTAProcessor::process ==========";

	Util::resetValueFactor(gridData.begin(), gridData.end());

	const std::size_t numSignals = (config_.lastTxElem - config_.firstTxElem + 1U) * config_.numElements;

	// Prepare the signal matrix.
	for (unsigned int txElem = config_.firstTxElem; txElem <= config_.lastTxElem; ++txElem) {
		LOG_INFO << "ACQ/PREP txElem: " << txElem << " <= " << config_.lastTxElem;

		acquisition_.execute(baseElement, txElem, acqData_);

		const std::size_t samplesPerChannelLow = acqData_.n2();

		const unsigned int localTxElem = txElem - config_.firstTxElem;
		for (unsigned int rxElem = 0; rxElem < config_.numElements; ++rxElem) {

			interpolator_.interpolate(&acqData_(rxElem, 0), samplesPerChannelLow, &tempSignal_[0]);

			// Copy the signal to the signal matrix.
			typename Matrix3<FloatType>::Dim3Interval interval = signalMatrix_.dim3Interval(localTxElem, rxElem);
			std::copy(tempSignal_.begin(), tempSignal_.end(), interval.first);

			Util::removeDC(&signalMatrix_(localTxElem, rxElem, 0), signalMatrix_.n3());
		}
	}

	std::vector<FloatType> xArray;
	ArrayGeometry::getElementXCentered2D(config_.numElements, config_.pitch, xArray);

	ThreadData threadData;
	tbb::enumerable_thread_specific<ThreadData> tls{threadData};

	const FloatType invCT = (config_.samplingFrequency * SIMPLE_STA_PROCESSOR_UPSAMPLING_FACTOR) / config_.propagationSpeed;
	const std::size_t numRows = gridData.n2();

	tbb::parallel_for(tbb::blocked_range<std::size_t>(0, gridData.n1()),
	[&, invCT, numRows](const tbb::blocked_range<std::size_t>& r) {
		auto& local = tls.local();

		local.delayList.resize(config_.numElements);

		// For each column:
		for (std::size_t i = r.begin(); i != r.end(); ++i) {
			// For each row:
			for (std::size_t j = 0; j < numRows; ++j) {

				FloatType pointValue = 0.0;
				XZValueFactor<FloatType>& point = gridData(i, j);

				// Calculate the delays.
				for (unsigned int elem = 0; elem < config_.numElements; ++elem) {
					const FloatType dx = point.x - xArray[elem];
					const FloatType dz = point.z /* - zArray*/; // zArray = 0
					local.delayList[elem] = std::sqrt(dx * dx + dz * dz) * invCT;
				}

				for (unsigned int txElem = config_.firstTxElem; txElem <= config_.lastTxElem; ++txElem) {
					const FloatType txDelay = local.delayList[txElem];
					const unsigned int localTxElem = txElem - config_.firstTxElem;
					for (unsigned int rxElem = 0; rxElem < config_.numElements; ++rxElem) {
#if 0
						// Nearest neighbor.
						const std::size_t delayIdx = static_cast<std::size_t>(FloatType{0.5} + signalOffset_ + txDelay + local.delayList[rxElem]);
						if (delayIdx < signalMatrix_.n3()) {
							pointValue += signalMatrix_(localTxElem, rxElem, delayIdx);
						}
#else
						// Linear interpolation.
						const FloatType delay = signalOffset_ + txDelay + local.delayList[rxElem];
						const std::size_t delayIdx = static_cast<std::size_t>(delay);
						const FloatType k = delay - delayIdx;
						if (delayIdx < signalMatrix_.n3() - 1) {
							const FloatType* p = &signalMatrix_(localTxElem, rxElem, delayIdx);
							pointValue += (1 - k) * *p + k * *(p + 1);
						}
#endif
					}
				}
				point.value = pointValue;
			}
		}
	});

	std::for_each(gridData.begin(), gridData.end(),
			Util::MultiplyValueBy<XZValueFactor<FloatType>, FloatType>{FloatType{1} / numSignals});

	LOG_DEBUG << "END ========== SimpleSTAProcessor::process ==========";
}

} // namespace Lab

#endif /* SIMPLESTAPROCESSOR_H_ */
