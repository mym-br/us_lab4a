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

#include "Exception.h"
#include "Util.h"



namespace Lab {
namespace Waveform {

// For types 2a, 2b and 2c, numPeriods doesn't significantly change the shape,
// if increased it only gives more time for the fade-out
// (reducing the discontinuity in the derivative at the end).
template<typename TFloat> void get(const std::string& type,
					TFloat centerFreq, TFloat samplingFreq, TFloat numPeriods,
					std::vector<TFloat>& v);

template<typename TFloat> void getType1(TFloat centerFreq, TFloat samplingFreq, TFloat numPeriods,
						std::vector<TFloat>& v);
template<typename TFloat> void getType2(TFloat centerFreq, TFloat samplingFreq, TFloat numPeriods, TFloat k,
						std::vector<TFloat>& v);
template<typename TFloat> void getType2a(TFloat centerFreq, TFloat samplingFreq, TFloat numPeriods,
						std::vector<TFloat>& v);
template<typename TFloat> void getType2b(TFloat centerFreq, TFloat samplingFreq, TFloat numPeriods,
						std::vector<TFloat>& v);
template<typename TFloat> void getType2c(TFloat centerFreq, TFloat samplingFreq, TFloat numPeriods,
						std::vector<TFloat>& v);



template<typename TFloat>
void
get(const std::string& type,
	TFloat centerFreq, TFloat samplingFreq, TFloat numPeriods,
	std::vector<TFloat>& v)
{
	if (type == "1") {
		getType1(centerFreq, samplingFreq, numPeriods, v);
	} else if (type == "2a") {
		getType2a(centerFreq, samplingFreq, numPeriods, v);
	} else if (type == "2b") {
		getType2b(centerFreq, samplingFreq, numPeriods, v);
	} else if (type == "2c") {
		getType2c(centerFreq, samplingFreq, numPeriods, v);
	} else {
		THROW_EXCEPTION(InvalidParameterException, "Invalid waveform type: " << type << '.');
	}
}

template<typename TFloat>
void
getType1(TFloat centerFreq, TFloat samplingFreq, TFloat numPeriods, std::vector<TFloat>& v)
{
	if (numPeriods <= 0.0) numPeriods = 3.0;

	const TFloat end = numPeriods / centerFreq;
	const TFloat dt = 1.0 / samplingFreq;
	const TFloat w = 2.0 * pi * centerFreq;
	const TFloat k = 1.0 / numPeriods;
	const unsigned int numSamples = static_cast<unsigned int>(end / dt) + 1U;
	v.resize(numSamples);
	for (unsigned int i = 0; i < numSamples; ++i) {
		const TFloat t = dt * i;
		v[i] = std::sin(w * t) * (TFloat(0.5) - TFloat(0.5) * std::cos(w * k * t));
	}
}

// Emeterio, J. L. S.
// Ullate, L. G.
// Diffraction impulse response of rectangular transducers.
// J. Acoust. Soc. Am., vol. 92, no. 2, pp. 651-662, 1992.
// DOI: 10.1121/1.403990
template<typename TFloat>
void
getType2(TFloat centerFreq, TFloat samplingFreq, TFloat numPeriods, TFloat k, std::vector<TFloat>& v)
{
	const TFloat end = numPeriods / centerFreq;
	const TFloat dt = 1.0 / samplingFreq;
	const TFloat w = 2.0 * pi * centerFreq;
	const unsigned int numSamples = static_cast<unsigned int>(end / dt) + 1U;
	v.resize(numSamples);
	TFloat maxVal = 0.0;
	for (unsigned int i = 0; i < numSamples; ++i) {
		const TFloat t = dt * i;
		v[i] = t * t * t * std::exp(-k * centerFreq * t) * std::cos(w * t);
		const TFloat absVal = std::abs(v[i]);
		if (absVal > maxVal) maxVal = absVal;
	}
	for (auto& elem : v) {
		elem /= maxVal;
	}
}
template<typename TFloat>
void
getType2a(TFloat centerFreq, TFloat samplingFreq, TFloat numPeriods, std::vector<TFloat>& v)
{
	if (numPeriods <= 0.0) numPeriods = 3.25;

	getType2(centerFreq, samplingFreq, numPeriods, TFloat(3.833), v);
}
template<typename TFloat>
void
getType2b(TFloat centerFreq, TFloat samplingFreq, TFloat numPeriods, std::vector<TFloat>& v)
{
	if (numPeriods <= 0.0) numPeriods = 8.25;

	getType2(centerFreq, samplingFreq, numPeriods, TFloat(1.437), v);
}

template<typename TFloat>
void
getType2c(TFloat centerFreq, TFloat samplingFreq, TFloat numPeriods, std::vector<TFloat>& v)
{
	if (numPeriods <= 0.0) numPeriods = 4.25;

	getType2(centerFreq, samplingFreq, numPeriods, TFloat(3.5), v);
}

} // namespace Waveform
} // namespace Lab

#endif /* WAVEFORM_H_ */
