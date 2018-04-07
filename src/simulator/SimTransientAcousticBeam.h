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

#include "FFTWFilter2.h"
#include "Matrix2.h"
#include "NumericRectangularFlatSourceImpulseResponse.h"
#include "Util.h"
#include "XYZ.h"
#include "XZValue.h"



namespace Lab {

template<typename FloatType>
class SimTransientAcousticBeam {
public:
	SimTransientAcousticBeam();
	~SimTransientAcousticBeam() {}

	void getRectangularFlatSourceAcousticBeam(FloatType propagationSpeed,
			FloatType samplingFreq,
			FloatType sourceWidth,
			FloatType sourceHeight,
			FloatType subElemSize,
			Matrix2<XYZ<FloatType>>& inputData,
			std::vector<FloatType>& dvdt,
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
					FloatType propagationSpeed,
					FloatType samplingFreq,
					FloatType sourceWidth,
					FloatType sourceHeight,
					FloatType subElemSize,
					Matrix2<XYZ<FloatType>>& inputData,
					std::vector<FloatType>& dvdt,
					Matrix2<XZValue<FloatType>>& gridData)
{
	std::size_t hOffset;
	std::vector<FloatType> h;
	auto impResp = std::make_unique<NumericRectangularFlatSourceImpulseResponse<FloatType>>(
									sourceWidth,
									sourceHeight,
									samplingFreq,
									propagationSpeed,
									subElemSize);

	std::vector<std::complex<FloatType>> filterFreqCoeff;
	filter_.setCoefficients(dvdt, filterFreqCoeff);
	std::vector<FloatType> signal;

	for (std::size_t i = 0, iEnd = inputData.n1(); i < iEnd; ++i) {
		for (std::size_t j = 0, jEnd = inputData.n2(); j < jEnd; ++j) {
			XYZ<FloatType>& id = inputData(i, j);
			impResp->getImpulseResponse(id.x, id.y, id.z, hOffset, h);

			filter_.filter(filterFreqCoeff, h, signal);

			FloatType minValue, maxValue;
			Util::minMax(signal, minValue, maxValue);
			gridData(i, j).value = maxValue - minValue;
		}
	}
}

} // namespace Lab

#endif // SIMTRANSIENTACOUSTICBEAM_H
