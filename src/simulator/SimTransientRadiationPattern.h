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
#include <memory>

#include <tbb/enumerable_thread_specific.h>
#include <tbb/tbb.h>

#include "ArrayOfRectangularFlatSourcesImpulseResponse.h"
#include "FFTWFilter2.h"
#include "Log.h"
#include "Matrix.h"
#include "Util.h"
#include "XYZ.h"
#include "XYZValue.h"

#define SIM_TRANSIENT_RADIATION_PATTERN_USE_MULTITHREADING 1



namespace Lab {

template<typename FloatType, typename ImpulseResponse>
class SimTransientRadiationPattern {
public:
	SimTransientRadiationPattern();
	~SimTransientRadiationPattern() {}

#ifdef SIM_TRANSIENT_RADIATION_PATTERN_USE_MULTITHREADING
	struct ThreadData {
		ThreadData(
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
#endif
	struct ArrayThreadData {
		ArrayThreadData(
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

	void getRectangularFlatSourceRadiationPattern(
			FloatType samplingFreq,
			FloatType propagationSpeed,
			FloatType sourceWidth,
			FloatType sourceHeight,
			FloatType discretization,
			const std::vector<FloatType>& dvdt,
			const Matrix<XYZ<FloatType>>& inputData,
			Matrix<XYZValue<FloatType>>& gridData);

	void getArrayOfRectangularFlatSourcesRadiationPattern(
			FloatType samplingFreq,
			FloatType propagationSpeed,
			FloatType sourceWidth,
			FloatType sourceHeight,
			FloatType discretization,
			const std::vector<FloatType>& dvdt,
			const std::vector<XY<FloatType>>& elemPos,
			const std::vector<FloatType>& focusDelay /* s */,
			const Matrix<XYZ<FloatType>>& inputData,
			Matrix<XYZValue<FloatType>>& gridData);
private:
	SimTransientRadiationPattern(const SimTransientRadiationPattern&) = delete;
	SimTransientRadiationPattern& operator=(const SimTransientRadiationPattern&) = delete;

	FFTWFilter2<FloatType> filter_;
};



template<typename FloatType, typename ImpulseResponse>
SimTransientRadiationPattern<FloatType, ImpulseResponse>::SimTransientRadiationPattern()
{
}

template<typename FloatType, typename ImpulseResponse>
void
SimTransientRadiationPattern<FloatType, ImpulseResponse>::getRectangularFlatSourceRadiationPattern(
					FloatType samplingFreq,
					FloatType propagationSpeed,
					FloatType sourceWidth,
					FloatType sourceHeight,
					FloatType discretization,
					const std::vector<FloatType>& dvdt,
					const Matrix<XYZ<FloatType>>& inputData,
					Matrix<XYZValue<FloatType>>& gridData)
{
#ifdef SIM_TRANSIENT_RADIATION_PATTERN_USE_MULTITHREADING
	ThreadData threadData{
		samplingFreq,
		propagationSpeed,
		sourceWidth,
		sourceHeight,
		discretization,
		dvdt
	};
	tbb::enumerable_thread_specific<ThreadData> tls(threadData);

	for (std::size_t i = 0, iEnd = inputData.n1(); i < iEnd; ++i) {
		LOG_INFO << "i: " << i << " < " << iEnd;

		tbb::parallel_for(tbb::blocked_range<std::size_t>(0, inputData.n2()),
			[&, i](const tbb::blocked_range<std::size_t>& r) {
				auto& local = tls.local();

				std::size_t hOffset;

				for (std::size_t j = r.begin(); j != r.end(); ++j) {
					const XYZ<FloatType>& id = inputData(i, j);
					local.ir.getImpulseResponse(id.x, id.y, id.z, hOffset, local.h);

					local.filter.filter(local.filterFreqCoeff, local.h, local.signal);

					FloatType minValue, maxValue;
					Util::minMax(local.signal, minValue, maxValue);
					gridData(i, j).value = maxValue - minValue;
				}
		});
	}
#else
	std::size_t hOffset;
	std::vector<FloatType> h;
	auto impResp = std::make_unique<NumericRectangularFlatSourceImpulseResponse<FloatType>>(
									samplingFreq,
									propagationSpeed,
									sourceWidth,
									sourceHeight,
									discretization);

	std::vector<std::complex<FloatType>> filterFreqCoeff;
	filter_.setCoefficients(dvdt, filterFreqCoeff);
	std::vector<FloatType> signal;

	for (std::size_t i = 0, iEnd = inputData.n1(); i < iEnd; ++i) {
		for (std::size_t j = 0, jEnd = inputData.n2(); j < jEnd; ++j) {
			const XYZ<FloatType>& id = inputData(i, j);
			impResp->getImpulseResponse(id.x, id.y, id.z, hOffset, h);

			filter_.filter(filterFreqCoeff, h, signal);

			FloatType minValue, maxValue;
			Util::minMax(signal, minValue, maxValue);
			gridData(i, j).value = maxValue - minValue;
		}
	}
#endif
}

template<typename FloatType, typename ImpulseResponse>
void
SimTransientRadiationPattern<FloatType, ImpulseResponse>::getArrayOfRectangularFlatSourcesRadiationPattern(
					FloatType samplingFreq,
					FloatType propagationSpeed,
					FloatType sourceWidth,
					FloatType sourceHeight,
					FloatType discretization,
					const std::vector<FloatType>& dvdt,
					const std::vector<XY<FloatType>>& elemPos,
					const std::vector<FloatType>& focusDelay,
					const Matrix<XYZ<FloatType>>& inputData,
					Matrix<XYZValue<FloatType>>& gridData)
{
	ArrayThreadData threadData{
		samplingFreq,
		propagationSpeed,
		sourceWidth,
		sourceHeight,
		discretization,
		elemPos,
		focusDelay,
		dvdt
	};
	tbb::enumerable_thread_specific<ArrayThreadData> tls(threadData);

	for (std::size_t i = 0, iEnd = inputData.n1(); i < iEnd; ++i) {
		LOG_INFO << "i: " << i << " < " << iEnd;

		tbb::parallel_for(tbb::blocked_range<std::size_t>(0, inputData.n2()),
			[&, i](const tbb::blocked_range<std::size_t>& r) {
				auto& local = tls.local();

				std::size_t hOffset;

				for (std::size_t j = r.begin(); j != r.end(); ++j) {
					const XYZ<FloatType>& id = inputData(i, j);
					local.ir.getImpulseResponse(id.x, id.y, id.z, hOffset, local.h);

					local.filter.filter(local.filterFreqCoeff, local.h, local.signal);

					FloatType minValue, maxValue;
					Util::minMax(local.signal, minValue, maxValue);
					gridData(i, j).value = maxValue - minValue;
				}
		});
	}
}

} // namespace Lab

#endif // SIMTRANSIENTRADIATIONPATTERN_H
