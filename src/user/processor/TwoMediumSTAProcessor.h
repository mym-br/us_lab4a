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

#ifndef TWOMEDIUMSTAPROCESSOR_H_
#define TWOMEDIUMSTAPROCESSOR_H_

#include <algorithm> /* fill */
#include <cmath> /* abs, sqrt */
#include <cstddef> /* std::size_t */
#include <limits>
#include <numeric> /* accumulate */
#include <vector>

#include <tbb/enumerable_thread_specific.h>
#include <tbb/tbb.h>

#include "CoherenceFactor.h"
#include "Exception.h"
#include "Interpolator4X.h"
#include "Log.h"
#include "Matrix.h"
#include "STAAcquisition.h"
#include "Tensor3.h"
#include "TwoMediumSTAConfiguration.h"
#include "Util.h"
#include "XZ.h"
#include "XZValueFactor.h"



namespace Lab {

template<typename TFloat>
class TwoMediumSTAProcessor {
public:
	TwoMediumSTAProcessor(const TwoMediumSTAConfiguration<TFloat>& config,
			STAAcquisition<TFloat>& acquisition,
			CoherenceFactorProcessor<TFloat>& coherenceFactor,
			TFloat peakOffset);
	~TwoMediumSTAProcessor() = default;

	void prepare(unsigned int baseElement);
	void process(const std::vector<XZ<TFloat>>& interfacePointList, Matrix<XZValueFactor<TFloat>>& gridData);

private:
	enum {
		UPSAMPLING_FACTOR = 4 // do not change
	};

	struct ThreadData {
		CoherenceFactorProcessor<TFloat> coherenceFactor;
		std::vector<TFloat> rxSignalSumList;
		std::vector<TFloat> delayList;
	};

	class ProcessColumn;

	TwoMediumSTAProcessor(const TwoMediumSTAProcessor&) = delete;
	TwoMediumSTAProcessor& operator=(const TwoMediumSTAProcessor&) = delete;
	TwoMediumSTAProcessor(TwoMediumSTAProcessor&&) = delete;
	TwoMediumSTAProcessor& operator=(TwoMediumSTAProcessor&&) = delete;

	const TwoMediumSTAConfiguration<TFloat>& config_;
	unsigned int deadZoneSamplesUp_;
	STAAcquisition<TFloat>& acquisition_;
	CoherenceFactorProcessor<TFloat>& coherenceFactor_;
	Tensor3<TFloat> signalMatrix_;
	typename STAAcquisition<TFloat>::AcquisitionDataType acqData_;
	TFloat signalOffset_;
	Interpolator4X<TFloat> interpolator_;
};



template<typename TFloat>
TwoMediumSTAProcessor<TFloat>::TwoMediumSTAProcessor(
			const TwoMediumSTAConfiguration<TFloat>& config,
			STAAcquisition<TFloat>& acquisition,
			CoherenceFactorProcessor<TFloat>& coherenceFactor,
			TFloat peakOffset)
		: config_(config)
		, deadZoneSamplesUp_((UPSAMPLING_FACTOR * config.samplingFrequency) * 2.0 * config.deadZoneM / config.propagationSpeed1)
		, acquisition_(acquisition)
		, coherenceFactor_(coherenceFactor)
{
	signalOffset_ = (config_.samplingFrequency * UPSAMPLING_FACTOR) * peakOffset / config_.centerFrequency;
	LOG_DEBUG << "signalOffset_=" << signalOffset_;

	// This acquisition is done to obtain the ascan length.
	acquisition_.prepare(0);
	acquisition_.execute(0, acqData_);
	const std::size_t signalLength = acqData_.n2() * UPSAMPLING_FACTOR;

	signalMatrix_.resize(config_.numElements, config_.numElements, signalLength);
	//signalMatrix_ = 0.0;

	if (deadZoneSamplesUp_ >= signalLength) {
		THROW_EXCEPTION(InvalidValueException, "Invalid dead zone: deadZoneSamplesUp_ (" << deadZoneSamplesUp_ <<
							") >= signalLength (" << signalLength << ").");
	}
}

template<typename TFloat>
void
TwoMediumSTAProcessor<TFloat>::prepare(unsigned int baseElement)
{
	acquisition_.prepare(baseElement);
}

template<typename TFloat>
void
TwoMediumSTAProcessor<TFloat>::process(const std::vector<XZ<TFloat>>& interfacePointList,
						Matrix<XZValueFactor<TFloat>>& gridData)
{
	LOG_DEBUG << "BEGIN ========== TwoMediumSTAProcessor::process ==========";

	// Clear the values in gridData.
	for (typename Matrix<XZValueFactor<TFloat>>::Iterator iter = gridData.begin(); iter != gridData.end(); ++iter) {
		iter->value = 0.0;
		iter->factor = 1.0;
	}

	// Prepare the signal matrix.
	for (std::size_t txElem = 0; txElem < config_.numElements; ++txElem) {
		acquisition_.execute(txElem, acqData_);

		const std::size_t samplesPerChannelLow = acqData_.n2();
		//const std::size_t samplesPerChannel = samplesPerChannelLow * STA_PROCESSOR_UPSAMPLING_FACTOR;

		for (std::size_t rxElem = 0; rxElem < config_.numElements; ++rxElem) {

			auto range = signalMatrix_.range3(txElem, rxElem);

			// Interpolate the signal.
			interpolator_.interpolate(&acqData_(rxElem, 0), samplesPerChannelLow, &(*range.begin()));

			Util::removeDC(&signalMatrix_(txElem, rxElem, 0), signalMatrix_.n3(), deadZoneSamplesUp_);
		}

		LOG_DEBUG << "PREP txElem = " << txElem;
	}

	ThreadData threadData;
	threadData.coherenceFactor = coherenceFactor_;

	tbb::enumerable_thread_specific<ThreadData> processColumnTLS(threadData);

	tbb::parallel_for(
		tbb::blocked_range<std::size_t>(0, gridData.n1()),
		ProcessColumn(
			gridData.n2(),
			config_,
			processColumnTLS,
			signalOffset_,
			signalMatrix_,
			interfacePointList,
			gridData));

	LOG_DEBUG << "END ========== TwoMediumSTAProcessor::process ==========";
}



template<typename TFloat>
class TwoMediumSTAProcessor<TFloat>::ProcessColumn {
public:
	ProcessColumn(
			std::size_t numRows,
			const TwoMediumSTAConfiguration<TFloat>& config,
			tbb::enumerable_thread_specific<ThreadData>& processColumnTLS,
			TFloat signalOffset,
			const Tensor3<TFloat>& signalMatrix,
			const std::vector<XZ<TFloat>>& interfacePointList,
			Matrix<XZValueFactor<TFloat>>& gridData)
				: numRows_(numRows)
				, config_(config)
				, processColumnTLS_(processColumnTLS)
				, signalOffset_(signalOffset)
				, signalMatrix_(signalMatrix)
				, fs_(config_.samplingFrequency * UPSAMPLING_FACTOR)
				, invC1_(1.0 / config_.propagationSpeed1)
				, invC1T_(fs_ * invC1_)
				, invC2_(1.0 / config_.propagationSpeed2)
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
			for (std::size_t j = 0; j < numRows_; ++j) {

				std::fill(data.rxSignalSumList.begin(), data.rxSignalSumList.end(), TFloat(0));
				XZValueFactor<TFloat>& point = gridData_(i, j);

				// Calculate the delays.
				for (unsigned int elem = 0; elem < config_.numElements; ++elem) {
					// Find the nearest point in the interface to the focus point, in the x direction.
					//TODO: use binary search???
					TFloat dxMin = std::numeric_limits<TFloat>::max();
					TFloat zDxMin = std::numeric_limits<TFloat>::max();
					for (typename std::vector<XZ<TFloat>>::const_iterator iter = interfacePointList_.begin(); iter != interfacePointList_.end(); ++iter) {
						const TFloat dx = std::abs(iter->x - point.x);
						if (dx < dxMin) {
							dxMin = dx;
							zDxMin = iter->z;
						}
					}

					const TFloat dz = point.z - zDxMin;
					if (dz <= 0.0) { // the point is below or at the interface
						const TFloat dx = point.x - xArray_[elem];
						const TFloat r = std::sqrt(dx * dx + point.z * point.z);
						data.delayList[elem] = r * invC1T_;
					} else {
						// Fermat's principle. Find the fastest path.
						//TODO: use binary search???
						TFloat tMin = std::numeric_limits<TFloat>::max();
						std::size_t idxMin = std::numeric_limits<std::size_t>::max();
						for (std::size_t pointIdx = 0, size = interfacePointList_.size(); pointIdx < size; ++pointIdx) {
							const XZ<TFloat>& interfacePoint = interfacePointList_[pointIdx];
							const TFloat dx1 = interfacePoint.x - xArray_[elem];
							const TFloat dz1 = interfacePoint.z;
							const TFloat dx2 = point.x - interfacePoint.x;
							const TFloat dz2 = point.z - interfacePoint.z;
							const TFloat t = std::sqrt(dx1 * dx1 + dz1 * dz1) * invC1_ + std::sqrt(dx2 * dx2 + dz2 * dz2) * invC2_;
							if (t < tMin) {
								tMin = t;
								idxMin = pointIdx;
							}
						}
						if (idxMin == 0 || idxMin == interfacePointList_.size() - 1) {
							// Ignore the points for which the interface is too short.
							data.delayList[elem] = signalMatrix_.n3(); // mark as invalid
						} else {
							data.delayList[elem] = tMin * fs_;
						}
					}
				}

				for (std::size_t txElem = 0; txElem < config_.numElements; ++txElem) {
					const TFloat txDelay = data.delayList[txElem];
					for (std::size_t rxElem = 0; rxElem < config_.numElements; ++rxElem) {
						// Linear interpolation.
						const TFloat delay = signalOffset_ + txDelay + data.delayList[rxElem];
						const std::size_t delayIdx = static_cast<std::size_t>(delay);
						const TFloat k = delay - delayIdx;
						if (delayIdx < signalMatrix_.n3() - 1) {
							const TFloat* p = &signalMatrix_(txElem, rxElem, delayIdx);
							data.rxSignalSumList[rxElem] += (1.0 - k) * *p + k * *(p + 1);
						}
					}
				}

				if (data.coherenceFactor.enabled()) {
					point.factor = data.coherenceFactor.calculate(&data.rxSignalSumList[0], data.rxSignalSumList.size());
				}
				point.value = std::accumulate(data.rxSignalSumList.begin(), data.rxSignalSumList.end(), TFloat(0));
			}
		}
	}
private:
	const std::size_t numRows_;
	const TwoMediumSTAConfiguration<TFloat>& config_;
	tbb::enumerable_thread_specific<ThreadData>& processColumnTLS_;
	const TFloat signalOffset_;
	const Tensor3<TFloat>& signalMatrix_;
	const TFloat fs_;
	const TFloat invC1_;
	const TFloat invC1T_;
	const TFloat invC2_;
	const std::vector<XZ<TFloat>>& interfacePointList_;
	std::vector<TFloat> xArray_;
	Matrix<XZValueFactor<TFloat>>& gridData_;
};

} // namespace Lab

#endif /* TWOMEDIUMSTAPROCESSOR_H_ */
