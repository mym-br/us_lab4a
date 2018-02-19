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

#include "Util.h"



namespace Lab {
namespace Waveform {

template<typename FloatType> void createPulseA(FloatType centerFrequency /* Hz */, FloatType samplingFreq /* Hz */, std::vector<FloatType>& v);



template<typename FloatType>
void
createPulseA(FloatType centerFrequency, FloatType samplingFreq, std::vector<FloatType>& v)
{
	const FloatType end = 4.25 / centerFrequency;
	const FloatType period = 1.0 / samplingFreq;
	const FloatType twoPi = 2.0 * PI;
	const unsigned int numPoints = static_cast<unsigned int>(end / period) + 1;
	v.resize(numPoints);
	for (unsigned int i = 0; i < numPoints; ++i) {
		const FloatType t = period * i;
		v[i] = t * t * t * std::exp(-3.5 * centerFrequency * t) * std::cos(twoPi * centerFrequency * t);
	}
}

} // namespace Waveform
} // namespace Lab

#endif /* WAVEFORM_H_ */
