#ifndef KAISERWINDOW_H_
#define KAISERWINDOW_H_

#include <cmath> /* ceil, pow */
#include <climits>
#include <vector>

#include <boost/math/special_functions/bessel.hpp>

#include "Exception.h"
#include "Util.h"



namespace Lab {
namespace KaiserWindow {

// tolerance_dB > 0.0
template<typename FloatType> FloatType getBeta(FloatType tolerance_dB);

// tolerance_dB > 0.0
// transition_width: 1.0 --> fs/2
template<typename FloatType> unsigned int getSize(FloatType tolerance_dB, FloatType transitionWidth);



template<typename FloatType>
FloatType getBeta(FloatType tolerance_dB)
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
unsigned int getSize(FloatType tolerance_dB, FloatType transitionWidth)
{
	// Scipy and Matlab use 7.95.
	// Oppenheim and Octave use 8.0.
	double size = std::ceil((tolerance_dB - 7.95) / (2.285 * transitionWidth * PI) + 1.0);
	if (size < 2.0 || size >= static_cast<double>(std::numeric_limits<unsigned int>::max())) {
		THROW_EXCEPTION(InvalidValueException, "Invalid Kaiser window size: " << size << '.');
	}
	return static_cast<unsigned int>(size);
}

template<typename FloatType>
void getWindow(unsigned int size, FloatType beta, std::vector<FloatType>& window)
{
	if (size < 2) {
		THROW_EXCEPTION(InvalidValueException, "Invalid Kaiser window size: " << size << '.');
	}

	window.resize(size);

	const FloatType c1 = 2 / static_cast<FloatType>(size - 1U);
	const FloatType c2 = 1 / boost::math::cyl_bessel_i(0, beta);
	for (unsigned int n = 0; n < size; ++n) {
		const FloatType k = c1 * n - 1;
		window[n] = c2 * boost::math::cyl_bessel_i(0, beta * std::sqrt(1 - k * k));
	}
}

} /* namespace KaiserWindow */
} /* namespace Lab */

#endif // KAISERWINDOW_H_
