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
#ifndef SIMTRANSIENTACOUSTICFIELD_H
#define SIMTRANSIENTACOUSTICFIELD_H

#include <complex>
#include <vector>

#include <tbb/enumerable_thread_specific.h>
#include <tbb/tbb.h>

#include "ArrayOfRectangularFlatSourcesImpulseResponse.h"
#include "FFTWFilter2.h"
#include "IterationCounter.h"
#include "Log.h"
#include "Matrix.h"
#include "Util.h"
#include "XYZValue.h"

namespace Lab {

template<typename FloatType, typename ImpulseResponse>
class SimTransientAcousticField {
public:
	SimTransientAcousticField();
	~SimTransientAcousticField() {}

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

	void getRectangularFlatSourceAcousticField(
			FloatType samplingFreq,
			FloatType propagationSpeed,
			FloatType sourceWidth,
			FloatType sourceHeight,
			FloatType discretization,
			const std::vector<FloatType>& dvdt,
			Matrix<XYZValue<FloatType>>& gridData);

	void getArrayOfRectangularFlatSourcesAcousticField(
			FloatType samplingFreq,
			FloatType propagationSpeed,
			FloatType sourceWidth,
			FloatType sourceHeight,
			FloatType discretization,
			const std::vector<FloatType>& dvdt,
			const std::vector<XY<FloatType>>& elemPos,
			const std::vector<FloatType>& focusDelay /* s */,
			Matrix<XYZValue<FloatType>>& gridData);
private:
	SimTransientAcousticField(const SimTransientAcousticField&) = delete;
	SimTransientAcousticField& operator=(const SimTransientAcousticField&) = delete;

	FFTWFilter2<FloatType> filter_;
};



template<typename FloatType, typename ImpulseResponse>
SimTransientAcousticField<FloatType, ImpulseResponse>::SimTransientAcousticField()
{
}

template<typename FloatType, typename ImpulseResponse>
void
SimTransientAcousticField<FloatType, ImpulseResponse>::getRectangularFlatSourceAcousticField(
					FloatType samplingFreq,
					FloatType propagationSpeed,
					FloatType sourceWidth,
					FloatType sourceHeight,
					FloatType discretization,
					const std::vector<FloatType>& dvdt,
					Matrix<XYZValue<FloatType>>& gridData)
{
	ThreadData threadData{
		samplingFreq,
		propagationSpeed,
		sourceWidth,
		sourceHeight,
		discretization,
		dvdt
	};
	tbb::enumerable_thread_specific<ThreadData> tls(threadData);

	IterationCounter::reset(gridData.n1());

	for (std::size_t i = 0, iEnd = gridData.n1(); i < iEnd; ++i) {
		tbb::parallel_for(tbb::blocked_range<std::size_t>(0, gridData.n2()),
			[&, i](const tbb::blocked_range<std::size_t>& r) {
				auto& local = tls.local();

				std::size_t hOffset;

				for (std::size_t j = r.begin(); j != r.end(); ++j) {
					XYZValue<FloatType>& point = gridData(i, j);
					local.ir.getImpulseResponse(point.x, point.y, point.z, hOffset, local.h);

					local.filter.filter(local.filterFreqCoeff, local.h, local.signal);

					FloatType minValue, maxValue;
					Util::minMax(local.signal, minValue, maxValue);
					point.value = maxValue - minValue;
				}
		});

		IterationCounter::add(1);
	}
}

template<typename FloatType, typename ImpulseResponse>
void
SimTransientAcousticField<FloatType, ImpulseResponse>::getArrayOfRectangularFlatSourcesAcousticField(
					FloatType samplingFreq,
					FloatType propagationSpeed,
					FloatType sourceWidth,
					FloatType sourceHeight,
					FloatType discretization,
					const std::vector<FloatType>& dvdt,
					const std::vector<XY<FloatType>>& elemPos,
					const std::vector<FloatType>& focusDelay,
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

	IterationCounter::reset(gridData.n1());

	for (std::size_t i = 0, iEnd = gridData.n1(); i < iEnd; ++i) {
		tbb::parallel_for(tbb::blocked_range<std::size_t>(0, gridData.n2()),
			[&, i](const tbb::blocked_range<std::size_t>& r) {
				auto& local = tls.local();

				std::size_t hOffset;

				for (std::size_t j = r.begin(); j != r.end(); ++j) {
					XYZValue<FloatType>& point = gridData(i, j);
					local.ir.getImpulseResponse(point.x, point.y, point.z, hOffset, local.h);

					local.filter.filter(local.filterFreqCoeff, local.h, local.signal);

					FloatType minValue, maxValue;
					Util::minMax(local.signal, minValue, maxValue);
					point.value = maxValue - minValue;
				}
		});

		IterationCounter::add(1);
	}
}

} // namespace Lab

#endif // SIMTRANSIENTACOUSTICFIELD_H
