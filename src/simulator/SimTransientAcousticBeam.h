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
#ifndef SIMTRANSIENTACOUSTICBEAM_H
#define SIMTRANSIENTACOUSTICBEAM_H

#include <complex>
#include <memory>

#include <tbb/enumerable_thread_specific.h>
#include <tbb/tbb.h>

#include "FFTWFilter2.h"
#include "Log.h"
#include "Matrix2.h"
#include "NumericArrayOfRectangularFlatSourcesImpulseResponse.h"
#include "NumericRectangularFlatSourceImpulseResponse.h"
#include "Util.h"
#include "XYZ.h"
#include "XZValue.h"

#define SIM_TRANSIENT_ACOUSTIC_BEAM_USE_MULTITHREADING 1



namespace Lab {

template<typename FloatType>
class SimTransientAcousticBeam {
public:
	SimTransientAcousticBeam();
	~SimTransientAcousticBeam() {}

#ifdef SIM_TRANSIENT_ACOUSTIC_BEAM_USE_MULTITHREADING
	struct ThreadData {
		ThreadData(
			FloatType samplingFreq,
			FloatType propagationSpeed,
			FloatType sourceWidth,
			FloatType sourceHeight,
			FloatType subElemSize,
			const std::vector<FloatType>& dvdt)
				: ir(samplingFreq, propagationSpeed, sourceWidth, sourceHeight, subElemSize)
		{
			filter.setCoefficients(dvdt, filterFreqCoeff);
		}
		NumericRectangularFlatSourceImpulseResponse<FloatType> ir;
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
			FloatType subElemSize,
			const std::vector<XY<FloatType>>& elemPos,
			const std::vector<FloatType>& focusDelay,
			const std::vector<FloatType>& dvdt)
				: ir(samplingFreq, propagationSpeed, sourceWidth, sourceHeight, subElemSize,
					elemPos, focusDelay)
		{
			filter.setCoefficients(dvdt, filterFreqCoeff);
		}
		NumericArrayOfRectangularFlatSourcesImpulseResponse<FloatType> ir;
		std::vector<std::complex<FloatType>> filterFreqCoeff;
		std::vector<FloatType> h;
		std::vector<FloatType> signal;
		FFTWFilter2<FloatType> filter;
	};

	void getRectangularFlatSourceAcousticBeam(
			FloatType samplingFreq,
			FloatType propagationSpeed,
			FloatType sourceWidth,
			FloatType sourceHeight,
			FloatType subElemSize,
			const std::vector<FloatType>& dvdt,
			const Matrix2<XYZ<FloatType>>& inputData,
			Matrix2<XZValue<FloatType>>& gridData);

	void getArrayOfRectangularFlatSourcesAcousticBeam(
			FloatType samplingFreq,
			FloatType propagationSpeed,
			FloatType sourceWidth,
			FloatType sourceHeight,
			FloatType subElemSize,
			const std::vector<FloatType>& dvdt,
			const std::vector<XY<FloatType>>& elemPos,
			const std::vector<FloatType>& focusDelay /* s */,
			const Matrix2<XYZ<FloatType>>& inputData,
			Matrix2<XZValue<FloatType>>& gridData);
private:
	SimTransientAcousticBeam(const SimTransientAcousticBeam&) = delete;
	SimTransientAcousticBeam& operator=(const SimTransientAcousticBeam&) = delete;

	FFTWFilter2<FloatType> filter_;
};



template<typename FloatType>
SimTransientAcousticBeam<FloatType>::SimTransientAcousticBeam()
{
}

template<typename FloatType>
void
SimTransientAcousticBeam<FloatType>::getRectangularFlatSourceAcousticBeam(
					FloatType samplingFreq,
					FloatType propagationSpeed,
					FloatType sourceWidth,
					FloatType sourceHeight,
					FloatType subElemSize,
					const std::vector<FloatType>& dvdt,
					const Matrix2<XYZ<FloatType>>& inputData,
					Matrix2<XZValue<FloatType>>& gridData)
{
#ifdef SIM_TRANSIENT_ACOUSTIC_BEAM_USE_MULTITHREADING
	ThreadData threadData{
		samplingFreq,
		propagationSpeed,
		sourceWidth,
		sourceHeight,
		subElemSize,
		dvdt
	};
	tbb::enumerable_thread_specific<ThreadData> tls{threadData};

	for (std::size_t i = 0, iEnd = inputData.n1(); i < iEnd; ++i) {
		LOG_DEBUG << "i = " << i << " / " << iEnd - 1;

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
									subElemSize);

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

template<typename FloatType>
void
SimTransientAcousticBeam<FloatType>::getArrayOfRectangularFlatSourcesAcousticBeam(
					FloatType samplingFreq,
					FloatType propagationSpeed,
					FloatType sourceWidth,
					FloatType sourceHeight,
					FloatType subElemSize,
					const std::vector<FloatType>& dvdt,
					const std::vector<XY<FloatType>>& elemPos,
					const std::vector<FloatType>& focusDelay,
					const Matrix2<XYZ<FloatType>>& inputData,
					Matrix2<XZValue<FloatType>>& gridData)
{
	ArrayThreadData threadData{
		samplingFreq,
		propagationSpeed,
		sourceWidth,
		sourceHeight,
		subElemSize,
		elemPos,
		focusDelay,
		dvdt
	};
	tbb::enumerable_thread_specific<ArrayThreadData> tls{threadData};

	for (std::size_t i = 0, iEnd = inputData.n1(); i < iEnd; ++i) {
		LOG_DEBUG << "i = " << i << " / " << iEnd - 1;

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

#endif // SIMTRANSIENTACOUSTICBEAM_H
