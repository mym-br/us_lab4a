/***************************************************************************
 *  Copyright 2018, 2019 Marcelo Y. Matuda                                 *
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
#ifndef SIMTRANSIENTPROPAGATION_H
#define SIMTRANSIENTPROPAGATION_H

#include <complex>
#include <vector>

#include <tbb/enumerable_thread_specific.h>
#include <tbb/tbb.h>

#include "ArrayOfRectangularFlatSourcesImpulseResponse.h"
#include "FFTWFilter2.h"
#include "IterationCounter.h"
#include "Log.h"
#include "Matrix.h"
#include "XYZValueArray.h"

namespace Lab {

template<typename FloatType, typename ImpulseResponse>
class SimTransientPropagation {
public:
	SimTransientPropagation();
	~SimTransientPropagation() {}

	struct CircularSourceThreadData {
		CircularSourceThreadData(
			FloatType samplingFreq,
			FloatType propagationSpeed,
			FloatType sourceRadius,
			FloatType discretization,
			const std::vector<FloatType>& dvdt)
				: ir(samplingFreq, propagationSpeed, sourceRadius, discretization)
		{
			filter.setCoefficients(dvdt, filterFreqCoeff);
		}
		ImpulseResponse ir;
		std::vector<std::complex<FloatType>> filterFreqCoeff;
		std::vector<FloatType> h;
		std::vector<FloatType> signal;
		FFTWFilter2<FloatType> filter;
	};

	struct RectangularSourceThreadData {
		RectangularSourceThreadData(
			FloatType samplingFreq,
			FloatType propagationSpeed,
			FloatType sourceWidth,
			FloatType sourceHeight,
			FloatType discretization,
			const std::vector<FloatType>& dvdt)
				: ir(samplingFreq, propagationSpeed, sourceWidth, sourceHeight, discretization)
		{
			filter.setCoefficients(dvdt, filterFreqCoeff);
		}
		ImpulseResponse ir;
		std::vector<std::complex<FloatType>> filterFreqCoeff;
		std::vector<FloatType> h;
		std::vector<FloatType> signal;
		FFTWFilter2<FloatType> filter;
	};

	struct ArrayOfRectangularSourcesThreadData {
		ArrayOfRectangularSourcesThreadData(
			FloatType samplingFreq,
			FloatType propagationSpeed,
			FloatType sourceWidth,
			FloatType sourceHeight,
			FloatType discretization,
			const std::vector<XY<FloatType>>& elemPos,
			const std::vector<FloatType>& focusDelay,
			const std::vector<FloatType>& dvdt)
				: ir(samplingFreq, propagationSpeed, sourceWidth, sourceHeight, discretization,
					elemPos, focusDelay)
		{
			filter.setCoefficients(dvdt, filterFreqCoeff);
		}
		ArrayOfRectangularFlatSourcesImpulseResponse<FloatType, ImpulseResponse> ir;
		std::vector<std::complex<FloatType>> filterFreqCoeff;
		std::vector<FloatType> h;
		std::vector<FloatType> signal;
		FFTWFilter2<FloatType> filter;
	};

	void getCircularSourcePropagation(
			FloatType samplingFreq,
			FloatType propagationSpeed,
			FloatType sourceRadius,
			FloatType discretization,
			const std::vector<FloatType>& dvdt,
			const std::vector<unsigned int>& propagIndexList,
			Matrix<XYZValueArray<FloatType>>& gridData);

	void getRectangularFlatSourcePropagation(
			FloatType samplingFreq,
			FloatType propagationSpeed,
			FloatType sourceWidth,
			FloatType sourceHeight,
			FloatType discretization,
			const std::vector<FloatType>& dvdt,
			const std::vector<unsigned int>& propagIndexList,
			Matrix<XYZValueArray<FloatType>>& gridData);

	void getArrayOfRectangularFlatSourcesPropagation(
			FloatType samplingFreq,
			FloatType propagationSpeed,
			FloatType sourceWidth,
			FloatType sourceHeight,
			FloatType discretization,
			const std::vector<FloatType>& dvdt,
			const std::vector<XY<FloatType>>& elemPos,
			const std::vector<FloatType>& focusDelay /* s */,
			const std::vector<unsigned int>& propagIndexList,
			Matrix<XYZValueArray<FloatType>>& gridData);
private:
	SimTransientPropagation(const SimTransientPropagation&) = delete;
	SimTransientPropagation& operator=(const SimTransientPropagation&) = delete;

	FFTWFilter2<FloatType> filter_;
};



template<typename FloatType, typename ImpulseResponse>
SimTransientPropagation<FloatType, ImpulseResponse>::SimTransientPropagation()
{
}

template<typename FloatType, typename ImpulseResponse>
void
SimTransientPropagation<FloatType, ImpulseResponse>::getCircularSourcePropagation(
					FloatType samplingFreq,
					FloatType propagationSpeed,
					FloatType sourceRadius,
					FloatType discretization,
					const std::vector<FloatType>& dvdt,
					const std::vector<unsigned int>& propagIndexList,
					Matrix<XYZValueArray<FloatType>>& gridData)
{
	CircularSourceThreadData threadData{
		samplingFreq,
		propagationSpeed,
		sourceRadius,
		discretization,
		dvdt
	};
	tbb::enumerable_thread_specific<CircularSourceThreadData> tls(threadData);

	IterationCounter::reset(gridData.n1());

	for (std::size_t i = 0, iEnd = gridData.n1(); i < iEnd; ++i) {
		tbb::parallel_for(tbb::blocked_range<std::size_t>(0, gridData.n2()),
			[&, i](const tbb::blocked_range<std::size_t>& r) {
				auto& local = tls.local();

				std::size_t hOffset;

				for (std::size_t j = r.begin(); j != r.end(); ++j) {
					XYZValueArray<FloatType>& point = gridData(i, j);
					local.ir.getImpulseResponse(point.x, point.y, point.z, hOffset, local.h);

					local.filter.filter(local.filterFreqCoeff, local.h, local.signal);

					point.values.resize(propagIndexList.size());
					for (unsigned int i = 0, end = propagIndexList.size(); i < end; ++i) {
						const unsigned int index = propagIndexList[i];
						if (index < hOffset) {
							point.values[i] = 0;
						} else {
							const unsigned int localIndex = index - hOffset;
							if (localIndex < local.signal.size()) {
								point.values[i] = local.signal[localIndex];
							} else {
								point.values[i] = 0;
							}
						}
					}
				}
		});

		IterationCounter::add(1);
	}
}

template<typename FloatType, typename ImpulseResponse>
void
SimTransientPropagation<FloatType, ImpulseResponse>::getRectangularFlatSourcePropagation(
					FloatType samplingFreq,
					FloatType propagationSpeed,
					FloatType sourceWidth,
					FloatType sourceHeight,
					FloatType discretization,
					const std::vector<FloatType>& dvdt,
					const std::vector<unsigned int>& propagIndexList,
					Matrix<XYZValueArray<FloatType>>& gridData)
{
	RectangularSourceThreadData threadData{
		samplingFreq,
		propagationSpeed,
		sourceWidth,
		sourceHeight,
		discretization,
		dvdt
	};
	tbb::enumerable_thread_specific<RectangularSourceThreadData> tls(threadData);

	IterationCounter::reset(gridData.n1());

	for (std::size_t i = 0, iEnd = gridData.n1(); i < iEnd; ++i) {
		tbb::parallel_for(tbb::blocked_range<std::size_t>(0, gridData.n2()),
			[&, i](const tbb::blocked_range<std::size_t>& r) {
				auto& local = tls.local();

				std::size_t hOffset;

				for (std::size_t j = r.begin(); j != r.end(); ++j) {
					XYZValueArray<FloatType>& point = gridData(i, j);
					local.ir.getImpulseResponse(point.x, point.y, point.z, hOffset, local.h);

					local.filter.filter(local.filterFreqCoeff, local.h, local.signal);

					point.values.resize(propagIndexList.size());
					for (unsigned int i = 0, end = propagIndexList.size(); i < end; ++i) {
						const unsigned int index = propagIndexList[i];
						if (index < hOffset) {
							point.values[i] = 0;
						} else {
							const unsigned int localIndex = index - hOffset;
							if (localIndex < local.signal.size()) {
								point.values[i] = local.signal[localIndex];
							} else {
								point.values[i] = 0;
							}
						}
					}
				}
		});

		IterationCounter::add(1);
	}
}

template<typename FloatType, typename ImpulseResponse>
void
SimTransientPropagation<FloatType, ImpulseResponse>::getArrayOfRectangularFlatSourcesPropagation(
					FloatType samplingFreq,
					FloatType propagationSpeed,
					FloatType sourceWidth,
					FloatType sourceHeight,
					FloatType discretization,
					const std::vector<FloatType>& dvdt,
					const std::vector<XY<FloatType>>& elemPos,
					const std::vector<FloatType>& focusDelay,
					const std::vector<unsigned int>& propagIndexList,
					Matrix<XYZValueArray<FloatType>>& gridData)
{
	ArrayOfRectangularSourcesThreadData threadData{
		samplingFreq,
		propagationSpeed,
		sourceWidth,
		sourceHeight,
		discretization,
		elemPos,
		focusDelay,
		dvdt
	};
	tbb::enumerable_thread_specific<ArrayOfRectangularSourcesThreadData> tls(threadData);

	IterationCounter::reset(gridData.n1());

	for (std::size_t i = 0, iEnd = gridData.n1(); i < iEnd; ++i) {
		tbb::parallel_for(tbb::blocked_range<std::size_t>(0, gridData.n2()),
			[&, i](const tbb::blocked_range<std::size_t>& r) {
				auto& local = tls.local();

				std::size_t hOffset;

				for (std::size_t j = r.begin(); j != r.end(); ++j) {
					XYZValueArray<FloatType>& point = gridData(i, j);
					local.ir.getImpulseResponse(point.x, point.y, point.z, hOffset, local.h);

					local.filter.filter(local.filterFreqCoeff, local.h, local.signal);

					point.values.resize(propagIndexList.size());
					for (unsigned int i = 0, end = propagIndexList.size(); i < end; ++i) {
						const unsigned int index = propagIndexList[i];
						if (index < hOffset) {
							point.values[i] = 0;
						} else {
							const unsigned int localIndex = index - hOffset;
							if (localIndex < local.signal.size()) {
								point.values[i] = local.signal[localIndex];
							} else {
								point.values[i] = 0;
							}
						}
					}
				}
		});

		IterationCounter::add(1);
	}
}

} // namespace Lab

#endif // SIMTRANSIENTPROPAGATION_H
