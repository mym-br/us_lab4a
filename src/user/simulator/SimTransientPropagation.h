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

#include "ArrayOfRectangularSourcesImpulseResponse.h"
#include "FFTWFilter2.h"
#include "IterationCounter.h"
#include "Log.h"
#include "Matrix.h"
#include "XYZValueArray.h"

namespace Lab {

template<typename TFloat, typename ImpulseResponse>
class SimTransientPropagation {
public:
	struct CircularSourceThreadData {
		CircularSourceThreadData(
			TFloat samplingFreq,
			TFloat propagationSpeed,
			TFloat sourceRadius,
			TFloat discretization,
			const std::vector<TFloat>& dvdt)
				: ir(samplingFreq, propagationSpeed, sourceRadius, discretization)
		{
			filter.setCoefficients(dvdt, filterFreqCoeff);
		}
		ImpulseResponse ir;
		std::vector<std::complex<TFloat>> filterFreqCoeff;
		std::vector<TFloat> h;
		std::vector<TFloat> signal;
		FFTWFilter2<TFloat> filter;
	};

	struct RectangularSourceThreadData {
		RectangularSourceThreadData(
			TFloat samplingFreq,
			TFloat propagationSpeed,
			TFloat sourceWidth,
			TFloat sourceHeight,
			TFloat discretization,
			const std::vector<TFloat>& dvdt)
				: ir(samplingFreq, propagationSpeed, sourceWidth, sourceHeight, discretization)
		{
			filter.setCoefficients(dvdt, filterFreqCoeff);
		}
		ImpulseResponse ir;
		std::vector<std::complex<TFloat>> filterFreqCoeff;
		std::vector<TFloat> h;
		std::vector<TFloat> signal;
		FFTWFilter2<TFloat> filter;
	};

	struct ArrayOfRectangularSourcesThreadData {
		ArrayOfRectangularSourcesThreadData(
			TFloat samplingFreq,
			TFloat propagationSpeed,
			TFloat sourceWidth,
			TFloat sourceHeight,
			TFloat discretization,
			const std::vector<XY<TFloat>>& elemPos,
			const std::vector<TFloat>& focusDelay,
			const std::vector<TFloat>& dvdt)
				: ir(samplingFreq, propagationSpeed, sourceWidth, sourceHeight, discretization,
					elemPos, focusDelay)
		{
			filter.setCoefficients(dvdt, filterFreqCoeff);
		}
		ArrayOfRectangularSourcesImpulseResponse<TFloat, ImpulseResponse> ir;
		std::vector<std::complex<TFloat>> filterFreqCoeff;
		std::vector<TFloat> h;
		std::vector<TFloat> signal;
		FFTWFilter2<TFloat> filter;
	};

	static void getCircularSourcePropagation(
			TFloat samplingFreq,
			TFloat propagationSpeed,
			TFloat sourceRadius,
			TFloat discretization,
			const std::vector<TFloat>& dvdt,
			const std::vector<unsigned int>& propagIndexList,
			Matrix<XYZValueArray<TFloat>>& gridData);

	static void getCircularSourcePropagationSingleThread(
			TFloat samplingFreq,
			TFloat propagationSpeed,
			TFloat sourceRadius,
			TFloat discretization,
			const std::vector<TFloat>& dvdt,
			const std::vector<unsigned int>& propagIndexList,
			Matrix<XYZValueArray<TFloat>>& gridData);

	static void getRectangularSourcePropagation(
			TFloat samplingFreq,
			TFloat propagationSpeed,
			TFloat sourceWidth,
			TFloat sourceHeight,
			TFloat discretization,
			const std::vector<TFloat>& dvdt,
			const std::vector<unsigned int>& propagIndexList,
			Matrix<XYZValueArray<TFloat>>& gridData);

	static void getRectangularSourcePropagationSingleThread(
			TFloat samplingFreq,
			TFloat propagationSpeed,
			TFloat sourceWidth,
			TFloat sourceHeight,
			TFloat discretization,
			const std::vector<TFloat>& dvdt,
			const std::vector<unsigned int>& propagIndexList,
			Matrix<XYZValueArray<TFloat>>& gridData);

	static void getArrayOfRectangularSourcesPropagation(
			TFloat samplingFreq,
			TFloat propagationSpeed,
			TFloat sourceWidth,
			TFloat sourceHeight,
			TFloat discretization,
			const std::vector<TFloat>& dvdt,
			const std::vector<XY<TFloat>>& elemPos,
			const std::vector<TFloat>& focusDelay /* s */,
			const std::vector<unsigned int>& propagIndexList,
			Matrix<XYZValueArray<TFloat>>& gridData);

	static void getArrayOfRectangularSourcesPropagationSingleThread(
			TFloat samplingFreq,
			TFloat propagationSpeed,
			TFloat sourceWidth,
			TFloat sourceHeight,
			TFloat discretization,
			const std::vector<TFloat>& dvdt,
			const std::vector<XY<TFloat>>& elemPos,
			const std::vector<TFloat>& focusDelay /* s */,
			const std::vector<unsigned int>& propagIndexList,
			Matrix<XYZValueArray<TFloat>>& gridData);
private:
	SimTransientPropagation() = delete;
};



template<typename TFloat, typename ImpulseResponse>
void
SimTransientPropagation<TFloat, ImpulseResponse>::getCircularSourcePropagation(
					TFloat samplingFreq,
					TFloat propagationSpeed,
					TFloat sourceRadius,
					TFloat discretization,
					const std::vector<TFloat>& dvdt,
					const std::vector<unsigned int>& propagIndexList,
					Matrix<XYZValueArray<TFloat>>& gridData)
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
					XYZValueArray<TFloat>& point = gridData(i, j);
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

template<typename TFloat, typename ImpulseResponse>
void
SimTransientPropagation<TFloat, ImpulseResponse>::getCircularSourcePropagationSingleThread(
					TFloat samplingFreq,
					TFloat propagationSpeed,
					TFloat sourceRadius,
					TFloat discretization,
					const std::vector<TFloat>& dvdt,
					const std::vector<unsigned int>& propagIndexList,
					Matrix<XYZValueArray<TFloat>>& gridData)
{
	CircularSourceThreadData threadData{
		samplingFreq,
		propagationSpeed,
		sourceRadius,
		discretization,
		dvdt
	};

	IterationCounter::reset(gridData.n1());

	for (std::size_t i = 0, iEnd = gridData.n1(); i < iEnd; ++i) {
		std::size_t hOffset;
		for (std::size_t j = 0, jEnd = gridData.n2(); j < jEnd; ++j) {
			XYZValueArray<TFloat>& point = gridData(i, j);
			threadData.ir.getImpulseResponse(point.x, point.y, point.z, hOffset, threadData.h);

			threadData.filter.filter(threadData.filterFreqCoeff, threadData.h, threadData.signal);

			point.values.resize(propagIndexList.size());
			for (unsigned int i = 0, end = propagIndexList.size(); i < end; ++i) {
				const unsigned int index = propagIndexList[i];
				if (index < hOffset) {
					point.values[i] = 0;
				} else {
					const unsigned int localIndex = index - hOffset;
					if (localIndex < threadData.signal.size()) {
						point.values[i] = threadData.signal[localIndex];
					} else {
						point.values[i] = 0;
					}
				}
			}
		}

		IterationCounter::add(1);
	}
}

template<typename TFloat, typename ImpulseResponse>
void
SimTransientPropagation<TFloat, ImpulseResponse>::getRectangularSourcePropagation(
					TFloat samplingFreq,
					TFloat propagationSpeed,
					TFloat sourceWidth,
					TFloat sourceHeight,
					TFloat discretization,
					const std::vector<TFloat>& dvdt,
					const std::vector<unsigned int>& propagIndexList,
					Matrix<XYZValueArray<TFloat>>& gridData)
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
					XYZValueArray<TFloat>& point = gridData(i, j);
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

template<typename TFloat, typename ImpulseResponse>
void
SimTransientPropagation<TFloat, ImpulseResponse>::getRectangularSourcePropagationSingleThread(
					TFloat samplingFreq,
					TFloat propagationSpeed,
					TFloat sourceWidth,
					TFloat sourceHeight,
					TFloat discretization,
					const std::vector<TFloat>& dvdt,
					const std::vector<unsigned int>& propagIndexList,
					Matrix<XYZValueArray<TFloat>>& gridData)
{
	RectangularSourceThreadData threadData{
		samplingFreq,
		propagationSpeed,
		sourceWidth,
		sourceHeight,
		discretization,
		dvdt
	};

	IterationCounter::reset(gridData.n1());

	for (std::size_t i = 0, iEnd = gridData.n1(); i < iEnd; ++i) {
		std::size_t hOffset;
		for (std::size_t j = 0, jEnd = gridData.n2(); j < jEnd; ++j) {
			XYZValueArray<TFloat>& point = gridData(i, j);
			threadData.ir.getImpulseResponse(point.x, point.y, point.z, hOffset, threadData.h);

			threadData.filter.filter(threadData.filterFreqCoeff, threadData.h, threadData.signal);

			point.values.resize(propagIndexList.size());
			for (unsigned int i = 0, end = propagIndexList.size(); i < end; ++i) {
				const unsigned int index = propagIndexList[i];
				if (index < hOffset) {
					point.values[i] = 0;
				} else {
					const unsigned int localIndex = index - hOffset;
					if (localIndex < threadData.signal.size()) {
						point.values[i] = threadData.signal[localIndex];
					} else {
						point.values[i] = 0;
					}
				}
			}
		}

		IterationCounter::add(1);
	}
}

template<typename TFloat, typename ImpulseResponse>
void
SimTransientPropagation<TFloat, ImpulseResponse>::getArrayOfRectangularSourcesPropagation(
					TFloat samplingFreq,
					TFloat propagationSpeed,
					TFloat sourceWidth,
					TFloat sourceHeight,
					TFloat discretization,
					const std::vector<TFloat>& dvdt,
					const std::vector<XY<TFloat>>& elemPos,
					const std::vector<TFloat>& focusDelay,
					const std::vector<unsigned int>& propagIndexList,
					Matrix<XYZValueArray<TFloat>>& gridData)
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
					XYZValueArray<TFloat>& point = gridData(i, j);
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

template<typename TFloat, typename ImpulseResponse>
void
SimTransientPropagation<TFloat, ImpulseResponse>::getArrayOfRectangularSourcesPropagationSingleThread(
					TFloat samplingFreq,
					TFloat propagationSpeed,
					TFloat sourceWidth,
					TFloat sourceHeight,
					TFloat discretization,
					const std::vector<TFloat>& dvdt,
					const std::vector<XY<TFloat>>& elemPos,
					const std::vector<TFloat>& focusDelay,
					const std::vector<unsigned int>& propagIndexList,
					Matrix<XYZValueArray<TFloat>>& gridData)
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

	IterationCounter::reset(gridData.n1());

	for (std::size_t i = 0, iEnd = gridData.n1(); i < iEnd; ++i) {
		std::size_t hOffset;
		for (std::size_t j = 0, jEnd = gridData.n2(); j < jEnd; ++j) {
			XYZValueArray<TFloat>& point = gridData(i, j);
			threadData.ir.getImpulseResponse(point.x, point.y, point.z, hOffset, threadData.h);

			threadData.filter.filter(threadData.filterFreqCoeff, threadData.h, threadData.signal);

			point.values.resize(propagIndexList.size());
			for (unsigned int i = 0, end = propagIndexList.size(); i < end; ++i) {
				const unsigned int index = propagIndexList[i];
				if (index < hOffset) {
					point.values[i] = 0;
				} else {
					const unsigned int localIndex = index - hOffset;
					if (localIndex < threadData.signal.size()) {
						point.values[i] = threadData.signal[localIndex];
					} else {
						point.values[i] = 0;
					}
				}
			}
		}

		IterationCounter::add(1);
	}
}

} // namespace Lab

#endif // SIMTRANSIENTPROPAGATION_H
