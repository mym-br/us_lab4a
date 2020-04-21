/***************************************************************************
 *  Copyright 2018 Marcelo Y. Matuda                                       *
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
#ifndef SIMTRANSIENTRADIATIONPATTERN_H
#define SIMTRANSIENTRADIATIONPATTERN_H

#include <complex>
#include <vector>

#include <tbb/enumerable_thread_specific.h>
#include <tbb/tbb.h>

#include "ArrayOfRectangularSourcesImpulseResponse.h"
#include "FFTWFilter2.h"
#include "IterationCounter.h"
#include "Log.h"
#include "Matrix.h"
#include "Util.h"
#include "XYZ.h"
#include "XYZValue.h"

#define SIM_TRANSIENT_RADIATION_PATTERN_USE_MULTITHREADING 1



namespace Lab {

template<typename TFloat, typename ImpulseResponse>
class SimTransientRadiationPattern {
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

	struct DirectArrayOfRectangularSourcesThreadData {
		DirectArrayOfRectangularSourcesThreadData(
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
		ImpulseResponse ir;
		std::vector<std::complex<TFloat>> filterFreqCoeff;
		std::vector<TFloat> h;
		std::vector<TFloat> signal;
		FFTWFilter2<TFloat> filter;
	};

	static void getCircularSourceRadiationPattern(
			TFloat samplingFreq,
			TFloat propagationSpeed,
			TFloat sourceRadius,
			TFloat discretization,
			const std::vector<TFloat>& dvdt,
			const std::vector<XYZ<TFloat>>& inputData,
			std::vector<TFloat>& radData);

	static void getCircularSourceRadiationPatternSingleThread(
			TFloat samplingFreq,
			TFloat propagationSpeed,
			TFloat sourceRadius,
			TFloat discretization,
			const std::vector<TFloat>& dvdt,
			const std::vector<XYZ<TFloat>>& inputData,
			std::vector<TFloat>& radData);

	static void getRectangularSourceRadiationPattern(
			TFloat samplingFreq,
			TFloat propagationSpeed,
			TFloat sourceWidth,
			TFloat sourceHeight,
			TFloat discretization,
			const std::vector<TFloat>& dvdt,
			const Matrix<XYZ<TFloat>>& inputData,
			Matrix<XYZValue<TFloat>>& gridData);

	static void getRectangularSourceRadiationPatternSingleThread(
			TFloat samplingFreq,
			TFloat propagationSpeed,
			TFloat sourceWidth,
			TFloat sourceHeight,
			TFloat discretization,
			const std::vector<TFloat>& dvdt,
			const Matrix<XYZ<TFloat>>& inputData,
			Matrix<XYZValue<TFloat>>& gridData);

	static void getArrayOfRectangularSourcesRadiationPattern(
			TFloat samplingFreq,
			TFloat propagationSpeed,
			TFloat sourceWidth,
			TFloat sourceHeight,
			TFloat discretization,
			const std::vector<TFloat>& dvdt,
			const std::vector<XY<TFloat>>& elemPos,
			const std::vector<TFloat>& focusDelay /* s */,
			const Matrix<XYZ<TFloat>>& inputData,
			Matrix<XYZValue<TFloat>>& gridData);

	static void getArrayOfRectangularSourcesRadiationPatternSingleThread(
			TFloat samplingFreq,
			TFloat propagationSpeed,
			TFloat sourceWidth,
			TFloat sourceHeight,
			TFloat discretization,
			const std::vector<TFloat>& dvdt,
			const std::vector<XY<TFloat>>& elemPos,
			const std::vector<TFloat>& focusDelay /* s */,
			const Matrix<XYZ<TFloat>>& inputData,
			Matrix<XYZValue<TFloat>>& gridData);

	static void getArrayOfRectangularSourcesRadiationPatternDirectSingleThread(
			TFloat samplingFreq,
			TFloat propagationSpeed,
			TFloat sourceWidth,
			TFloat sourceHeight,
			TFloat discretization,
			const std::vector<TFloat>& dvdt,
			const std::vector<XY<TFloat>>& elemPos,
			const std::vector<TFloat>& focusDelay /* s */,
			const Matrix<XYZ<TFloat>>& inputData,
			Matrix<XYZValue<TFloat>>& gridData);
private:
	SimTransientRadiationPattern() = delete;
};



template<typename TFloat, typename ImpulseResponse>
void
SimTransientRadiationPattern<TFloat, ImpulseResponse>::getCircularSourceRadiationPattern(
					TFloat samplingFreq,
					TFloat propagationSpeed,
					TFloat sourceRadius,
					TFloat discretization,
					const std::vector<TFloat>& dvdt,
					const std::vector<XYZ<TFloat>>& inputData,
					std::vector<TFloat>& radData)
{
	CircularSourceThreadData threadData{
		samplingFreq,
		propagationSpeed,
		sourceRadius,
		discretization,
		dvdt
	};
	tbb::enumerable_thread_specific<CircularSourceThreadData> tls(threadData);

	tbb::parallel_for(tbb::blocked_range<std::size_t>(0, inputData.size()),
		[&](const tbb::blocked_range<std::size_t>& r) {
			auto& local = tls.local();

			std::size_t hOffset;

			for (std::size_t i = r.begin(); i != r.end(); ++i) {
				const XYZ<TFloat>& id = inputData[i];
				local.ir.getImpulseResponse(id.x, id.y, id.z, hOffset, local.h);

				local.filter.filter(local.filterFreqCoeff, local.h, local.signal);

				TFloat minValue, maxValue;
				Util::minMax(local.signal, minValue, maxValue);
				radData[i] = maxValue - minValue;
			}
	});
}

template<typename TFloat, typename ImpulseResponse>
void
SimTransientRadiationPattern<TFloat, ImpulseResponse>::getCircularSourceRadiationPatternSingleThread(
					TFloat samplingFreq,
					TFloat propagationSpeed,
					TFloat sourceRadius,
					TFloat discretization,
					const std::vector<TFloat>& dvdt,
					const std::vector<XYZ<TFloat>>& inputData,
					std::vector<TFloat>& radData)
{
	CircularSourceThreadData threadData{
		samplingFreq,
		propagationSpeed,
		sourceRadius,
		discretization,
		dvdt
	};

	std::size_t hOffset;
	for (std::size_t i = 0, iEnd = inputData.size(); i < iEnd; ++i) {
		const XYZ<TFloat>& id = inputData[i];
		threadData.ir.getImpulseResponse(id.x, id.y, id.z, hOffset, threadData.h);

		threadData.filter.filter(threadData.filterFreqCoeff, threadData.h, threadData.signal);

		TFloat minValue, maxValue;
		Util::minMax(threadData.signal, minValue, maxValue);
		radData[i] = maxValue - minValue;
	}
}

template<typename TFloat, typename ImpulseResponse>
void
SimTransientRadiationPattern<TFloat, ImpulseResponse>::getRectangularSourceRadiationPattern(
					TFloat samplingFreq,
					TFloat propagationSpeed,
					TFloat sourceWidth,
					TFloat sourceHeight,
					TFloat discretization,
					const std::vector<TFloat>& dvdt,
					const Matrix<XYZ<TFloat>>& inputData,
					Matrix<XYZValue<TFloat>>& gridData)
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

	IterationCounter::reset(inputData.n1());

	for (std::size_t i = 0, iEnd = inputData.n1(); i < iEnd; ++i) {
		tbb::parallel_for(tbb::blocked_range<std::size_t>(0, inputData.n2()),
			[&, i](const tbb::blocked_range<std::size_t>& r) {
				auto& local = tls.local();

				std::size_t hOffset;

				for (std::size_t j = r.begin(); j != r.end(); ++j) {
					const XYZ<TFloat>& id = inputData(i, j);
					local.ir.getImpulseResponse(id.x, id.y, id.z, hOffset, local.h);

					local.filter.filter(local.filterFreqCoeff, local.h, local.signal);

					TFloat minValue, maxValue;
					Util::minMax(local.signal, minValue, maxValue);
					gridData(i, j).value = maxValue - minValue;
				}
		});

		IterationCounter::add(1);
	}
}

template<typename TFloat, typename ImpulseResponse>
void
SimTransientRadiationPattern<TFloat, ImpulseResponse>::getRectangularSourceRadiationPatternSingleThread(
					TFloat samplingFreq,
					TFloat propagationSpeed,
					TFloat sourceWidth,
					TFloat sourceHeight,
					TFloat discretization,
					const std::vector<TFloat>& dvdt,
					const Matrix<XYZ<TFloat>>& inputData,
					Matrix<XYZValue<TFloat>>& gridData)
{
	RectangularSourceThreadData threadData{
		samplingFreq,
		propagationSpeed,
		sourceWidth,
		sourceHeight,
		discretization,
		dvdt
	};

	IterationCounter::reset(inputData.n1());

	for (std::size_t i = 0, iEnd = inputData.n1(); i < iEnd; ++i) {
		std::size_t hOffset;
		for (std::size_t j = 0, jEnd = inputData.n2(); j < jEnd; ++j) {
			const XYZ<TFloat>& id = inputData(i, j);
			threadData.ir.getImpulseResponse(id.x, id.y, id.z, hOffset, threadData.h);

			threadData.filter.filter(threadData.filterFreqCoeff, threadData.h, threadData.signal);

			TFloat minValue, maxValue;
			Util::minMax(threadData.signal, minValue, maxValue);
			gridData(i, j).value = maxValue - minValue;
		}

		IterationCounter::add(1);
	}
}

template<typename TFloat, typename ImpulseResponse>
void
SimTransientRadiationPattern<TFloat, ImpulseResponse>::getArrayOfRectangularSourcesRadiationPattern(
					TFloat samplingFreq,
					TFloat propagationSpeed,
					TFloat sourceWidth,
					TFloat sourceHeight,
					TFloat discretization,
					const std::vector<TFloat>& dvdt,
					const std::vector<XY<TFloat>>& elemPos,
					const std::vector<TFloat>& focusDelay,
					const Matrix<XYZ<TFloat>>& inputData,
					Matrix<XYZValue<TFloat>>& gridData)
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

	IterationCounter::reset(inputData.n1());

	for (std::size_t i = 0, iEnd = inputData.n1(); i < iEnd; ++i) {
		tbb::parallel_for(tbb::blocked_range<std::size_t>(0, inputData.n2()),
			[&, i](const tbb::blocked_range<std::size_t>& r) {
				auto& local = tls.local();

				std::size_t hOffset;

				for (std::size_t j = r.begin(); j != r.end(); ++j) {
					const XYZ<TFloat>& id = inputData(i, j);
					local.ir.getImpulseResponse(id.x, id.y, id.z, hOffset, local.h);

					local.filter.filter(local.filterFreqCoeff, local.h, local.signal);

					TFloat minValue, maxValue;
					Util::minMax(local.signal, minValue, maxValue);
					gridData(i, j).value = maxValue - minValue;
				}
		});

		IterationCounter::add(1);
	}
}

template<typename TFloat, typename ImpulseResponse>
void
SimTransientRadiationPattern<TFloat, ImpulseResponse>::getArrayOfRectangularSourcesRadiationPatternSingleThread(
					TFloat samplingFreq,
					TFloat propagationSpeed,
					TFloat sourceWidth,
					TFloat sourceHeight,
					TFloat discretization,
					const std::vector<TFloat>& dvdt,
					const std::vector<XY<TFloat>>& elemPos,
					const std::vector<TFloat>& focusDelay,
					const Matrix<XYZ<TFloat>>& inputData,
					Matrix<XYZValue<TFloat>>& gridData)
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

	IterationCounter::reset(inputData.n1());

	for (std::size_t i = 0, iEnd = inputData.n1(); i < iEnd; ++i) {
		std::size_t hOffset;
		for (std::size_t j = 0, jEnd = inputData.n2(); j < jEnd; ++j) {
			const XYZ<TFloat>& id = inputData(i, j);
			threadData.ir.getImpulseResponse(id.x, id.y, id.z, hOffset, threadData.h);

			threadData.filter.filter(threadData.filterFreqCoeff, threadData.h, threadData.signal);

			TFloat minValue, maxValue;
			Util::minMax(threadData.signal, minValue, maxValue);
			gridData(i, j).value = maxValue - minValue;
		}

		IterationCounter::add(1);
	}
}

template<typename TFloat, typename ImpulseResponse>
void
SimTransientRadiationPattern<TFloat, ImpulseResponse>::getArrayOfRectangularSourcesRadiationPatternDirectSingleThread(
					TFloat samplingFreq,
					TFloat propagationSpeed,
					TFloat sourceWidth,
					TFloat sourceHeight,
					TFloat discretization,
					const std::vector<TFloat>& dvdt,
					const std::vector<XY<TFloat>>& elemPos,
					const std::vector<TFloat>& focusDelay,
					const Matrix<XYZ<TFloat>>& inputData,
					Matrix<XYZValue<TFloat>>& gridData)
{
	DirectArrayOfRectangularSourcesThreadData threadData{
		samplingFreq,
		propagationSpeed,
		sourceWidth,
		sourceHeight,
		discretization,
		elemPos,
		focusDelay,
		dvdt
	};

	IterationCounter::reset(inputData.n1());

	for (std::size_t i = 0, iEnd = inputData.n1(); i < iEnd; ++i) {
		std::size_t hOffset;
		for (std::size_t j = 0, jEnd = inputData.n2(); j < jEnd; ++j) {
			const XYZ<TFloat>& id = inputData(i, j);
			threadData.ir.getImpulseResponse(id.x, id.y, id.z, hOffset, threadData.h);

			threadData.filter.filter(threadData.filterFreqCoeff, threadData.h, threadData.signal);

			TFloat minValue, maxValue;
			Util::minMax(threadData.signal, minValue, maxValue);
			gridData(i, j).value = maxValue - minValue;
		}

		IterationCounter::add(1);
	}
}

} // namespace Lab

#endif // SIMTRANSIENTRADIATIONPATTERN_H
