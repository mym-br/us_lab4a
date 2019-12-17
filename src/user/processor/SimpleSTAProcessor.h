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
#include <cstddef> /* std::size_t */
#include <vector>

#include <tbb/enumerable_thread_specific.h>
#include <tbb/tbb.h>

#include "ArrayGeometry.h"
#include "ArrayProcessor.h"
#include "Exception.h"
#include "Geometry.h"
#include "Interpolator4X.h"
#include "IterationCounter.h"
#include "Log.h"
#include "Matrix.h"
#include "STAAcquisition.h"
#include "STAConfiguration.h"
#include "Tensor3.h"
#include "Util.h"
#include "XYZValueFactor.h"



namespace Lab {

// This processor does not support coherence factors.
// x = 0 is at the center of the element group.
// y = 0
// z = 0 is at the surface of the array.
template<typename TFloat>
class SimpleSTAProcessor : public ArrayProcessor<TFloat> {
public:
	SimpleSTAProcessor(
			const STAConfiguration<TFloat>& config,
			STAAcquisition<TFloat>& acquisition,
			TFloat peakOffset);
	virtual ~SimpleSTAProcessor() = default;

	virtual void prepare(unsigned int baseElement);
	virtual void process(Matrix<XYZValueFactor<TFloat>>& gridData);
private:
	// Do not change.
	static constexpr unsigned int upsamplingFactor = 4;

	struct ThreadData {
		std::vector<TFloat> delayList;
	};

	SimpleSTAProcessor(const SimpleSTAProcessor&) = delete;
	SimpleSTAProcessor& operator=(const SimpleSTAProcessor&) = delete;
	SimpleSTAProcessor(SimpleSTAProcessor&&) = delete;
	SimpleSTAProcessor& operator=(SimpleSTAProcessor&&) = delete;

	const STAConfiguration<TFloat>& config_;
	bool initialized_;
	STAAcquisition<TFloat>& acquisition_;
	Tensor3<TFloat> signalTensor_;
	typename STAAcquisition<TFloat>::AcquisitionDataType acqData_;
	std::vector<TFloat> tempSignal_;
	TFloat signalOffset_;
	Interpolator4X<TFloat> interpolator_;
};



template<typename TFloat>
SimpleSTAProcessor<TFloat>::SimpleSTAProcessor(
			const STAConfiguration<TFloat>& config,
			STAAcquisition<TFloat>& acquisition,
			TFloat peakOffset)
		: config_(config)
		, initialized_()
		, acquisition_(acquisition)
{
	signalOffset_ = (config_.samplingFrequency * upsamplingFactor) * peakOffset / config_.centerFrequency;
}

template<typename TFloat>
void
SimpleSTAProcessor<TFloat>::prepare(unsigned int baseElement)
{
	acquisition_.prepare(baseElement);
}

template<typename TFloat>
void
SimpleSTAProcessor<TFloat>::process(Matrix<XYZValueFactor<TFloat>>& gridData)
{
	LOG_DEBUG << "BEGIN ========== SimpleSTAProcessor::process ==========";

	Util::resetValueFactor(gridData.begin(), gridData.end());

	const std::size_t numSignals = (config_.lastTxElem - config_.firstTxElem + 1U) * config_.numElements;

	// Prepare the signal matrix.
	for (unsigned int txElem = config_.firstTxElem; txElem <= config_.lastTxElem; ++txElem) {
		LOG_INFO << "ACQ/PREP txElem: " << txElem << " <= " << config_.lastTxElem;

		acquisition_.execute(txElem, acqData_);

		if (!initialized_) {
			const std::size_t signalLength = acqData_.n2() * upsamplingFactor;
			tempSignal_.resize(signalLength);
			signalTensor_.resize(config_.lastTxElem - config_.firstTxElem + 1, config_.numElements, signalLength);
			LOG_DEBUG << "signalOffset_=" << signalOffset_ << " signalLength=" << signalLength;
			initialized_ = true;
		}

		const std::size_t samplesPerChannelLow = acqData_.n2();

		const unsigned int localTxElem = txElem - config_.firstTxElem;
		for (unsigned int rxElem = 0; rxElem < config_.numElements; ++rxElem) {

			interpolator_.interpolate(&acqData_(rxElem, 0), samplesPerChannelLow, &tempSignal_[0]);

			// Copy the signal to the signal matrix.
			auto range = signalTensor_.range3(localTxElem, rxElem);
			std::copy(tempSignal_.begin(), tempSignal_.end(), range.begin());

			Util::removeDC(&signalTensor_(localTxElem, rxElem, 0), signalTensor_.n3());
		}
	}

	std::vector<TFloat> xArray;
	ArrayGeometry::getElementXCentered2D(config_.numElements, config_.pitch, xArray);

	ThreadData threadData;
	tbb::enumerable_thread_specific<ThreadData> tls(threadData);

	const TFloat invCT = (config_.samplingFrequency * upsamplingFactor) / config_.propagationSpeed;
	const std::size_t numRows = gridData.n2();

	IterationCounter::reset(gridData.n1());

	tbb::parallel_for(tbb::blocked_range<std::size_t>(0, gridData.n1()),
	[&, invCT, numRows](const tbb::blocked_range<std::size_t>& r) {
		auto& local = tls.local();

		local.delayList.resize(config_.numElements);

		// For each column:
		for (std::size_t i = r.begin(); i != r.end(); ++i) {
			// For each row:
			for (std::size_t j = 0; j < numRows; ++j) {

				TFloat pointValue = 0.0;
				XYZValueFactor<TFloat>& point = gridData(i, j);

				// Calculate the delays.
				for (unsigned int elem = 0; elem < config_.numElements; ++elem) {
					local.delayList[elem] = Geometry::distance2DY0(xArray[elem], point.x, point.z) * invCT;
				}

				for (unsigned int txElem = config_.firstTxElem; txElem <= config_.lastTxElem; ++txElem) {
					const TFloat txDelay = local.delayList[txElem];
					const unsigned int localTxElem = txElem - config_.firstTxElem;
					for (unsigned int rxElem = 0; rxElem < config_.numElements; ++rxElem) {
#if 0
						// Nearest neighbor.
						const std::size_t delayIdx = static_cast<std::size_t>(TFloat(0.5) + signalOffset_ + txDelay + local.delayList[rxElem]);
						if (delayIdx < signalTensor_.n3()) {
							pointValue += signalTensor_(localTxElem, rxElem, delayIdx);
						}
#else
						// Linear interpolation.
						const TFloat delay = signalOffset_ + txDelay + local.delayList[rxElem];
						const std::size_t delayIdx = static_cast<std::size_t>(delay);
						const TFloat k = delay - delayIdx;
						if (delayIdx + 1U < signalTensor_.n3()) {
							const TFloat* p = &signalTensor_(localTxElem, rxElem, delayIdx);
							pointValue += (1 - k) * *p + k * *(p + 1);
						}
#endif
					}
				}
				point.value = pointValue;
			}
		}

		IterationCounter::add(r.end() - r.begin());
	});

	std::for_each(gridData.begin(), gridData.end(),
			Util::MultiplyValueBy<XYZValueFactor<TFloat>, TFloat>(TFloat(1) / numSignals));

	LOG_DEBUG << "END ========== SimpleSTAProcessor::process ==========";
}

} // namespace Lab

#endif /* SIMPLESTAPROCESSOR_H_ */
