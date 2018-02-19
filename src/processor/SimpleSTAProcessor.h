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

#include <cmath>
#include <cstddef> /* std::size_t */

#include <tbb/enumerable_thread_specific.h>
#include <tbb/tbb.h>

#include "Exception.h"
#include "Interpolator4X.h"
#include "Log.h"
#include "Matrix2.h"
#include "STAAcquisition.h"
#include "STAConfiguration.h"
#include "STAProcessor.h"
#include "Util.h"
#include "XZValueFactor.h"

// Do not change.
#define SIMPLE_STA_PROCESSOR_UPSAMPLING_FACTOR 4



namespace Lab {

// This processor does not support SCF.
template<typename FloatType>
class SimpleSTAProcessor : public STAProcessor<FloatType> {
public:
	SimpleSTAProcessor(
			const STAConfiguration<FloatType>& config,
			STAAcquisition<FloatType>& acquisition,
			FloatType peakOffset);
	virtual ~SimpleSTAProcessor() {}

	virtual void process(unsigned int baseElement, Matrix2<XZValueFactor<FloatType> >& gridData);

private:
	struct ThreadData {
		std::vector<FloatType> delayList;
	};

	class CalculatePoints;

	SimpleSTAProcessor(const SimpleSTAProcessor&);
	SimpleSTAProcessor& operator=(const SimpleSTAProcessor&);

	const STAConfiguration<FloatType>& config_;
	STAAcquisition<FloatType>& acquisition_;
	Matrix2<FloatType> signalList_;
	typename STAAcquisition<FloatType>::AcquisitionDataType acqData_;
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
	signalOffset_ = (config_.samplingFrequency * SIMPLE_STA_PROCESSOR_UPSAMPLING_FACTOR) * peakOffset / config_.centerFrequency;
	LOG_DEBUG << "signalOffset_=" << signalOffset_;
}

template<typename FloatType>
void
SimpleSTAProcessor<FloatType>::process(unsigned int baseElement, Matrix2<XZValueFactor<FloatType> >& gridData)
{
	LOG_DEBUG << "BEGIN ========== SimpleSTAProcessor::process ==========";

	// Clears the values in gridData.
	for (typename Matrix2<XZValueFactor<FloatType> >::Iterator iter = gridData.begin(); iter != gridData.end(); ++iter) {
		iter->value = 0.0;
		iter->factor = 1.0;
	}

	//const std::size_t numDelayElem = delays_.n3();

	ThreadData threadData;

	// For each transmit element:
	for (std::size_t txElem = 0; txElem < config_.numElements; ++txElem) {

		acquisition_.execute(baseElement, txElem, acqData_);

		const std::size_t samplesPerChannelLow = acqData_.n2();
		const std::size_t samplesPerChannel = samplesPerChannelLow * SIMPLE_STA_PROCESSOR_UPSAMPLING_FACTOR;
		signalList_.resize(config_.numElements, samplesPerChannel);

		// For each receive element:
		for (std::size_t rxElem = 0; rxElem < config_.numElements; ++rxElem) {
			// Interpolates the signal.
			interpolator_.interpolate(&acqData_(rxElem, 0), samplesPerChannelLow, &signalList_(rxElem, 0));

			Util::removeDC(&signalList_(rxElem, 0), signalList_.n2());
		}

		tbb::enumerable_thread_specific<ThreadData> calculatePointsTLS(threadData);
		tbb::parallel_for(
			tbb::blocked_range<std::size_t>(0, gridData.n1()),
			CalculatePoints(
				gridData.n2(),
				config_,
				calculatePointsTLS,
				signalOffset_,
				txElem,
				signalList_,
				gridData));

		LOG_DEBUG << "PRC txElem = " << txElem;
	}

	LOG_DEBUG << "END ========== SimpleSTAProcessor::process ==========";
}



template<typename FloatType>
class SimpleSTAProcessor<FloatType>::CalculatePoints {
public:
	CalculatePoints(
			std::size_t numRows,
			const STAConfiguration<FloatType>& config,
			tbb::enumerable_thread_specific<ThreadData>& calculatePointsTLS,
			FloatType signalOffset,
			int txElem,
			const Matrix2<FloatType>& signalList,
			Matrix2<XZValueFactor<FloatType> >& gridData)
				: numRows_(numRows)
				, config_(config)
				, calculatePointsTLS_(calculatePointsTLS)
				, signalOffset_(signalOffset)
				, txElem_(txElem)
				, signalList_(signalList)
				, invCT_((config_.samplingFrequency * SIMPLE_STA_PROCESSOR_UPSAMPLING_FACTOR) / config_.propagationSpeed)
				, xArray_(config_.numElements)
				, gridData_(gridData) {

		const FloatType xArrayCenter = (config_.numElements - 1) * config_.pitch / 2.0;
		for (unsigned int elem = 0; elem < config_.numElements; ++elem) {
			xArray_[elem] = elem * config_.pitch - xArrayCenter;
		}
	}

	void operator()(const tbb::blocked_range<std::size_t>& r) const {
		ThreadData& data = calculatePointsTLS_.local();

		data.delayList.resize(config_.numElements);

		// For each point in x-z:
		for (std::size_t i = r.begin(); i != r.end(); ++i) {
			for (std::size_t j = 0; j < numRows_; ++j) {

				XZValueFactor<FloatType>& point = gridData_(i, j);

				// Calculates the delays.
				for (unsigned int elem = 0; elem < config_.numElements; ++elem) {
					const FloatType dx = point.x - xArray_[elem];
					const FloatType dz = point.z /* - zArray*/; // zArray = 0
					data.delayList[elem] = std::sqrt(dx * dx + dz * dz) * invCT_;
				}

				const FloatType txDelay = data.delayList[txElem_];
				FloatType pointValue = 0.0;
				for (std::size_t rxElem = 0; rxElem < config_.numElements; ++rxElem) {
#if 0
					// Nearest neighbor.
					const std::size_t delayIdx = static_cast<std::size_t>(0.5 + signalOffset_ + txDelay + data.delayList[rxElem]);
					if (delayIdx < signalList_.n2()) {
						pointValue += signalList_(rxElem, delayIdx);
					}
#else
					// Linear interpolation.
					const FloatType delay = signalOffset_ + txDelay + data.delayList[rxElem];
					const std::size_t delayIdx = static_cast<std::size_t>(delay);
					const FloatType k = delay - delayIdx;
					if (delayIdx < signalList_.n2() - 1) {
						const FloatType* p = &signalList_(rxElem, delayIdx);
						pointValue += (1.0 - k) * *p + k * *(p + 1);
					}
#endif
				}
				point.value += pointValue;
			}
		}
	}
private:
	const std::size_t numRows_;
	const STAConfiguration<FloatType>& config_;
	tbb::enumerable_thread_specific<ThreadData>& calculatePointsTLS_;
	const FloatType signalOffset_;
	const int txElem_;
	const Matrix2<FloatType>& signalList_;
	const FloatType invCT_;
	std::vector<FloatType> xArray_;
	Matrix2<XZValueFactor<FloatType> >& gridData_;
};

} // namespace Lab

#endif /* SIMPLESTAPROCESSOR_H_ */
