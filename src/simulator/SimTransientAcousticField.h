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
#include <memory>

#include <tbb/enumerable_thread_specific.h>
#include <tbb/tbb.h>

#include "FFTWFilter2.h"
#include "Log.h"
#include "Matrix2.h"
#include "NumericRectangularFlatSourceImpulseResponse.h"
#include "Util.h"
#include "XZValue.h"



namespace Lab {

template<typename FloatType>
class SimTransientAcousticField {
public:
	SimTransientAcousticField();
	~SimTransientAcousticField() {}

	struct ThreadData {
		ThreadData(
			FloatType sourceWidth,
			FloatType sourceHeight,
			FloatType samplingFreq,
			FloatType propagationSpeed,
			FloatType subElemSize,
			const std::vector<FloatType>& dvdt)
				: ir(sourceWidth, sourceHeight, samplingFreq, propagationSpeed, subElemSize)
		{
			filter.setCoefficients(dvdt, filterFreqCoeff);
		}
		NumericRectangularFlatSourceImpulseResponse<FloatType> ir;
		std::vector<std::complex<FloatType>> filterFreqCoeff;
		std::vector<FloatType> h;
		std::vector<FloatType> signal;
		FFTWFilter2<FloatType> filter;
	};

	void getRectangularFlatSourceAcousticField(
			FloatType sourceWidth,
			FloatType sourceHeight,
			FloatType samplingFreq,
			FloatType propagationSpeed,
			FloatType subElemSize,
			std::vector<FloatType>& dvdt,
			FloatType y,
			Matrix2<XZValue<FloatType>>& gridData);
private:
	SimTransientAcousticField(const SimTransientAcousticField&) = delete;
	SimTransientAcousticField& operator=(const SimTransientAcousticField&) = delete;

	FFTWFilter2<FloatType> filter_;
};



template<typename FloatType>
SimTransientAcousticField<FloatType>::SimTransientAcousticField()
{
}

template<typename FloatType>
void
SimTransientAcousticField<FloatType>::getRectangularFlatSourceAcousticField(
					FloatType sourceWidth,
					FloatType sourceHeight,
					FloatType samplingFreq,
					FloatType propagationSpeed,
					FloatType subElemSize,
					std::vector<FloatType>& dvdt,
					FloatType y,
					Matrix2<XZValue<FloatType>>& gridData)
{
	ThreadData threadData{
		sourceWidth,
		sourceHeight,
		samplingFreq,
		propagationSpeed,
		subElemSize,
		dvdt
	};
	tbb::enumerable_thread_specific<ThreadData> tls{threadData};

	for (std::size_t i = 0, iEnd = gridData.n1(); i < iEnd; ++i) {
		LOG_DEBUG << "i = " << i << " / " << iEnd - 1;

		tbb::parallel_for(tbb::blocked_range<std::size_t>(0, gridData.n2()),
			[&, i, y](const tbb::blocked_range<std::size_t>& r) {
				auto& local = tls.local();

				std::size_t hOffset;

				for (std::size_t j = r.begin(); j != r.end(); ++j) {
					XZValue<FloatType>& point = gridData(i, j);
					local.ir.getImpulseResponse(point.x, y, point.z, hOffset, local.h);

					local.filter.filter(local.filterFreqCoeff, local.h, local.signal);

					FloatType minValue, maxValue;
					Util::minMax(local.signal, minValue, maxValue);
					point.value = maxValue - minValue;
				}
		});
	}
}

} // namespace Lab

#endif // SIMTRANSIENTACOUSTICFIELD_H
