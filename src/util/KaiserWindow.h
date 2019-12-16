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

#ifndef KAISERWINDOW_H_
#define KAISERWINDOW_H_

#include <cmath> /* ceil, cyl_bessel_i, pow */
#include <climits>
#include <vector>

#include "Exception.h"
#include "Util.h"



namespace Lab {
namespace KaiserWindow {

// tolerance_dB > 0.0
template<typename FloatType> FloatType getBeta(FloatType tolerance_dB);

// tolerance_dB > 0.0
// transition_width: 1.0 --> fs/2
template<typename FloatType> unsigned int getSize(FloatType tolerance_dB, FloatType transitionWidth);

template<typename FloatType> void getWindow(unsigned int size, FloatType beta, std::vector<FloatType>& window);



template<typename FloatType>
FloatType
getBeta(FloatType tolerance_dB)
{
	FloatType beta;
	if (tolerance_dB > 50) {
		beta = 0.1102 * (tolerance_dB - 8.7);
	} else if (tolerance_dB > 21) {
		const double k = tolerance_dB - 21.0;
		beta = 0.5842 * std::pow(k, 0.4) + 0.07886 * k;
	} else {
		beta = 0.0;
	}

	return beta;
}

template<typename FloatType>
unsigned int
getSize(FloatType tolerance_dB, FloatType transitionWidth)
{
	// Scipy and Matlab use 7.95.
	// Oppenheim and Octave use 8.0.
	double size = std::ceil((tolerance_dB - 7.95) / (2.285 * transitionWidth * pi) + 1.0);
	if (size < 2.0 || size >= static_cast<double>(std::numeric_limits<unsigned int>::max())) {
		THROW_EXCEPTION(InvalidValueException, "Invalid Kaiser window size: " << size << '.');
	}
	return static_cast<unsigned int>(size);
}

template<typename FloatType>
void
getWindow(unsigned int size, FloatType beta, std::vector<FloatType>& window)
{
	if (size < 2) {
		THROW_EXCEPTION(InvalidValueException, "Invalid Kaiser window size: " << size << '.');
	}

	window.resize(size);

	const FloatType c1 = 2 / static_cast<FloatType>(size - 1U);
	const FloatType c2 = 1 / std::cyl_bessel_i(0, beta);
	for (unsigned int n = 0; n < size; ++n) {
		const FloatType k = c1 * n - 1;
		FloatType x = 1 - k * k;
		if (x < 0.0) x = 0.0;
		window[n] = c2 * std::cyl_bessel_i(0, beta * std::sqrt(x));
	}
}

} /* namespace KaiserWindow */
} /* namespace Lab */

#endif // KAISERWINDOW_H_
