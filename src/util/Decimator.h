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
#ifndef DECIMATOR_H
#define DECIMATOR_H

#include <cstddef> /* std::size_t */
#include <cmath> /* ceil, log10 */
#include <vector>

#include <boost/math/special_functions/sinc.hpp>

#include "Exception.h"
#include "FFTWFilter.h"
#include "KaiserWindow.h"
#include "Util.h"

#define DECIMATOR_KAISER_TOLERANCE (1.0e-4)



namespace Lab {

// This class is copy constructible and assignable.
template<typename FloatType>
class Decimator {
public:
	Decimator();
	~Decimator();

	// lpFilterTransitionWidth:
	//     half the total transition width
	//     1.0 -> pi radian / sample at the destination sampling rate
	void prepare(unsigned int downsamplingFactor, FloatType lpFilterHalfTransitionWidth);

	void decimate(std::size_t inputOffset, const std::vector<FloatType>& input,
			std::size_t& outputOffset, std::vector<FloatType>& output);
private:
	unsigned int downsamplingFactor_;
	std::vector<FloatType> lowPassFIRFilter_;
	std::vector<FloatType> filteredSignal_;
	FFTWFilter<FloatType> filter_;
};

template<typename FloatType>
Decimator<FloatType>::Decimator()
		: downsamplingFactor_{}
{
}

template<typename FloatType>
Decimator<FloatType>::~Decimator()
{
}

template<typename FloatType>
void
Decimator<FloatType>::prepare(unsigned int downsamplingFactor, FloatType lpFilterHalfTransitionWidth)
{
	if (downsamplingFactor_ != 0) {
		THROW_EXCEPTION(InvalidStateException, "The object is already configured.");
	}
	if (downsamplingFactor < 2) {
		THROW_EXCEPTION(InvalidValueException, "Invalid downsampling factor: " << downsamplingFactor << ". Must be >= 2.");
	}

	//TODO: limit downsampling factor

	const FloatType tol_dB = -20.0 * std::log10(DECIMATOR_KAISER_TOLERANCE);
	const FloatType kaiserBeta = KaiserWindow::getBeta(tol_dB);

	//TODO:(???) Check finalTransitionWidth <= 1.0/downsamplingFactor or lpFilterTransitionWidth <= 0.5
	const FloatType finalTransitionWidth = (lpFilterHalfTransitionWidth * 2) / downsamplingFactor;
	const unsigned int windowSize = KaiserWindow::getSize(tol_dB, finalTransitionWidth);

	const FloatType twoDownFactor = 2 * downsamplingFactor;
	const unsigned int finalWindowSize = static_cast<unsigned int>(
				std::ceil((windowSize - 1) / twoDownFactor) * twoDownFactor) + 1U;

	// Kaiser window.
	std::vector<FloatType> window;
	KaiserWindow::getWindow(finalWindowSize, kaiserBeta, window);

	// Sinc.
	const FloatType numPeriods = (finalWindowSize - 1U) / twoDownFactor;
	std::vector<FloatType> x;
	Util::fillSequenceFromStartToEndWithSize(x, -numPeriods, numPeriods, finalWindowSize);
	lowPassFIRFilter_.resize(x.size());
	for (unsigned int i = 0; i < lowPassFIRFilter_.size(); ++i) {
		lowPassFIRFilter_[i] = boost::math::sinc_pi(PI * x[i]);
	}

	Util::multiply(window, lowPassFIRFilter_);

	filter_.setCoefficients(lowPassFIRFilter_);

	downsamplingFactor_ = downsamplingFactor;
}

template<typename FloatType>
void
Decimator<FloatType>::decimate(std::size_t inputOffset, const std::vector<FloatType>& input,
				std::size_t& outputOffset, std::vector<FloatType>& output)
{
	if (downsamplingFactor_ == 0) THROW_EXCEPTION(InvalidStateException, "The decimator has not been initialized.");

	// Apply the low-pass filter.
	filter_.filter(input, filteredSignal_);

	const unsigned int m = inputOffset % downsamplingFactor_;
	const unsigned int firstIndex = (m == 0) ? 0 : downsamplingFactor_ - m;

	// Copy to the output, compensating for the FIR filter delay. The signal is truncated at both ends.
	output.resize(1U + (input.size() - firstIndex - 1U) / downsamplingFactor_);
	const std::size_t filterOffset = (lowPassFIRFilter_.size() - 1) / 2;
	const FloatType outputFactor = 1.0 / downsamplingFactor_;
	for (std::size_t i = filterOffset + firstIndex, iEnd = filteredSignal_.size() - filterOffset, j = 0;
			i < iEnd;
			i += downsamplingFactor_, ++j) {
		output[j] = filteredSignal_[i] * outputFactor;
	}

	outputOffset = (inputOffset + firstIndex) / downsamplingFactor_;
}

} // namespace Lab

#endif // DECIMATOR_H
