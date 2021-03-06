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

#ifndef INTERPOLATOR_H_
#define INTERPOLATOR_H_

#include <cmath> /* ceil, log10 */
#include <vector>

#include <boost/math/special_functions/sinc.hpp>

#include "Exception.h"
#include "FFTWFilter.h"
#include "KaiserWindow.h"
#include "Util.h"



namespace Lab {

// This class is copy constructible and assignable.
template<typename TFloat>
class Interpolator {
public:
	Interpolator();

	// lpFilterTransitionWidth:
	//     half the total transition width
	//     1.0 -> pi radian / sample at the original sampling rate
	void prepare(unsigned int upsamplingFactor, TFloat lpFilterHalfTransitionWidth);

	void interpolate(const TFloat* input, std::size_t inputLength, TFloat* output);
private:
	static constexpr double kaiserTolerance = 1.0e-4;
	static constexpr unsigned int maxUpFactor = 64;

	unsigned int upsamplingFactor_;
	std::vector<TFloat> lowPassFIRFilter_;
	std::vector<TFloat> inputVector_;
	std::vector<TFloat> outputVector_;
	FFTWFilter<TFloat> filter_;
};

template<typename TFloat>
Interpolator<TFloat>::Interpolator()
		: upsamplingFactor_()
{
}

template<typename TFloat>
void
Interpolator<TFloat>::prepare(unsigned int upsamplingFactor, TFloat lpFilterHalfTransitionWidth)
{
	if (upsamplingFactor_ != 0) {
		THROW_EXCEPTION(InvalidStateException, "The object is already configured.");
	}
	if (upsamplingFactor < 2) {
		THROW_EXCEPTION(InvalidValueException, "Invalid upsampling factor: " << upsamplingFactor << ". Must be >= 2.");
	}
	if (upsamplingFactor > maxUpFactor) {
		THROW_EXCEPTION(InvalidValueException, "Invalid upsampling factor: " << upsamplingFactor
				<< ". Must be <= " << maxUpFactor << '.');
	}
	if (lpFilterHalfTransitionWidth > 1) {
		THROW_EXCEPTION(InvalidValueException, "Invalid transition width: " << lpFilterHalfTransitionWidth
				<< ". Must be <= 1.0");
	}

	const TFloat tol_dB = -20.0 * std::log10(kaiserTolerance);
	const TFloat kaiserBeta = KaiserWindow::getBeta(tol_dB);

	const TFloat finalTransitionWidth = (lpFilterHalfTransitionWidth * 2) / upsamplingFactor;
	const unsigned int windowSize = KaiserWindow::getSize(tol_dB, finalTransitionWidth);

	const TFloat twoUpFactor = 2 * upsamplingFactor;
	const unsigned int finalWindowSize = static_cast<unsigned int>(
				std::ceil((windowSize - 1) / twoUpFactor) * twoUpFactor) + 1U;

	// Kaiser window.
	std::vector<TFloat> window;
	KaiserWindow::getWindow(finalWindowSize, kaiserBeta, window);

	// Sinc.
	const TFloat numPeriods = (finalWindowSize - 1U) / twoUpFactor;
	std::vector<TFloat> x;
	Util::fillSequenceFromStartToEndWithSize(x, -numPeriods, numPeriods, finalWindowSize);
	lowPassFIRFilter_.resize(x.size());
	for (unsigned int i = 0; i < lowPassFIRFilter_.size(); ++i) {
		lowPassFIRFilter_[i] = boost::math::sinc_pi(pi * x[i]);
	}

	Util::multiply(window, lowPassFIRFilter_);

	filter_.setCoefficients(lowPassFIRFilter_);

	upsamplingFactor_ = upsamplingFactor;
}

// The argument "ouput" must point to an array of size inputLength * upsamplingFactor.
template<typename TFloat>
void
Interpolator<TFloat>::interpolate(const TFloat* input, std::size_t inputLength, TFloat* output)
{
	if (upsamplingFactor_ == 0) THROW_EXCEPTION(InvalidStateException, "The interpolator has not been initialized.");

	// Upsample.
	inputVector_.resize(inputLength * upsamplingFactor_);
	for (std::size_t i = 0, j = 0; i < inputLength; ++i, j += upsamplingFactor_) {
		inputVector_[j] = input[i];
		for (std::size_t k = j + 1, end = j + upsamplingFactor_; k < end; ++k) {
			inputVector_[k] = 0;
		}
	}

	// Apply the anti-aliasing filter.
	filter_.filter(inputVector_, outputVector_);

	// Copy to the output, compensating for the FIR filter delay. The signal is truncated at both ends.
	const std::size_t offset = (lowPassFIRFilter_.size() - 1) / 2;
//	for (std::size_t i = 0; i < inputLength * UPSAMPLING_FACTOR; ++i) {
//		output[i] = outputVector_[i + offset];
//	}
	memcpy(output, &outputVector_[offset], inputLength * upsamplingFactor_ * sizeof(TFloat));
}

} // namespace Lab

#endif // INTERPOLATOR_H_
