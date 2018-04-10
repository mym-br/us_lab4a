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

#ifndef WAVEFORM_H_
#define WAVEFORM_H_

#include <cmath>
#include <vector>

#include "Util.h"



namespace Lab {
namespace Waveform {

template<typename FloatType> void getType1(FloatType centerFreq, FloatType samplingFreq,
						std::vector<FloatType>& v, FloatType numPeriods);
template<typename FloatType> void getType2(FloatType centerFreq, FloatType samplingFreq, FloatType numPeriods,
						FloatType k, std::vector<FloatType>& v);
template<typename FloatType> void getType2a(FloatType centerFreq, FloatType samplingFreq,
						std::vector<FloatType>& v, FloatType numPeriods);
template<typename FloatType> void getType2b(FloatType centerFreq, FloatType samplingFreq,
						std::vector<FloatType>& v, FloatType numPeriods);
template<typename FloatType> void getType2c(FloatType centerFreq, FloatType samplingFreq,
						std::vector<FloatType>& v, FloatType numPeriods);



template<typename FloatType>
void getType1(FloatType centerFreq, FloatType samplingFreq, std::vector<FloatType>& v, FloatType numPeriods)
{
	if (numPeriods <= 0.0) numPeriods = 3.0;

	const FloatType end = numPeriods / centerFreq;
	const FloatType dt = 1.0 / samplingFreq;
	const FloatType w = 2.0 * PI * centerFreq;
	const FloatType k = 1.0 / numPeriods;
	const unsigned int numSamples = static_cast<unsigned int>(end / dt) + 1U;
	v.resize(numSamples);
	for (unsigned int i = 0; i < numSamples; ++i) {
		const FloatType t = dt * i;
		v[i] = std::sin(w * t) * (FloatType{0.5} - FloatType{0.5} * std::cos(w * k * t));
	}
}

// Emeterio, J. L. S.
// Ullate, L. G.
// Diffraction impulse response of rectangular transducers.
// J. Acoust. Soc. Am., vol. 92, no. 2, pp. 651-662, 1992.
// DOI: 10.1121/1.403990
template<typename FloatType>
void getType2(FloatType centerFreq, FloatType samplingFreq, FloatType numPeriods, FloatType k, std::vector<FloatType>& v)
{
	const FloatType end = numPeriods / centerFreq;
	const FloatType dt = 1.0 / samplingFreq;
	const FloatType w = 2.0 * PI * centerFreq;
	const unsigned int numSamples = static_cast<unsigned int>(end / dt) + 1U;
	v.resize(numSamples);
	FloatType maxVal = 0.0;
	for (unsigned int i = 0; i < numSamples; ++i) {
		const FloatType t = dt * i;
		v[i] = t * t * t * std::exp(-k * centerFreq * t) * std::cos(w * t);
		const FloatType absVal = std::abs(v[i]);
		if (absVal > maxVal) maxVal = absVal;
	}
	for (auto& elem : v) {
		elem /= maxVal;
	}
}
template<typename FloatType>
void getType2a(FloatType centerFreq, FloatType samplingFreq, std::vector<FloatType>& v, FloatType numPeriods)
{
	if (numPeriods <= 0.0) numPeriods = 3.25;

	getType2(centerFreq, samplingFreq, numPeriods, FloatType{3.833}, v);
}
template<typename FloatType>
void getType2b(FloatType centerFreq, FloatType samplingFreq, std::vector<FloatType>& v, FloatType numPeriods)
{
	if (numPeriods <= 0.0) numPeriods = 8.25;

	getType2(centerFreq, samplingFreq, numPeriods, FloatType{1.437}, v);
}

template<typename FloatType>
void getType2c(FloatType centerFreq, FloatType samplingFreq, std::vector<FloatType>& v, FloatType numPeriods)
{
	if (numPeriods <= 0.0) numPeriods = 4.25;

	getType2(centerFreq, samplingFreq, numPeriods, FloatType{3.5}, v);
}

} // namespace Waveform
} // namespace Lab

#endif /* WAVEFORM_H_ */
