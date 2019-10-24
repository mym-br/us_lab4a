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

	// Reentrant.
	void downsample(std::size_t filterOffset,
			std::size_t inputOffset, const std::vector<FloatType>& input,
			std::size_t& outputOffset, std::vector<FloatType>& output);

	void decimate(std::size_t inputOffset, const std::vector<FloatType>& input,
			std::size_t& outputOffset, std::vector<FloatType>& output);

	const std::vector<FloatType>& lowPassFIRFilter() const { return lowPassFIRFilter_; }
	unsigned int downsamplingFactor() const { return downsamplingFactor_; }
private:
	static constexpr unsigned int maxDownFactor = 10;
	static constexpr FloatType kaiserTolerance = 1.0e-4;

	unsigned int downsamplingFactor_;
	std::vector<FloatType> lowPassFIRFilter_;
	std::vector<FloatType> filteredSignal_;
	FFTWFilter<FloatType> filter_;
};

template<typename FloatType>
Decimator<FloatType>::Decimator()
		: downsamplingFactor_()
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
	if (downsamplingFactor > maxDownFactor) {
		THROW_EXCEPTION(InvalidValueException, "Invalid downsampling factor: " << downsamplingFactor
				<< ". Must be <= " << maxDownFactor << '.');
	}
	if (lpFilterHalfTransitionWidth > 1) {
		THROW_EXCEPTION(InvalidValueException, "Invalid transition width: " << lpFilterHalfTransitionWidth
				<< ". Must be <= 1.0");
	}

	const FloatType tol_dB = -20.0 * std::log10(kaiserTolerance);
	const FloatType kaiserBeta = KaiserWindow::getBeta(tol_dB);

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
		lowPassFIRFilter_[i] = boost::math::sinc_pi(pi * x[i]) / downsamplingFactor;
	}

	Util::multiply(window, lowPassFIRFilter_);

	filter_.setCoefficients(lowPassFIRFilter_);

	downsamplingFactor_ = downsamplingFactor;
}

template<typename FloatType>
void
Decimator<FloatType>::downsample(std::size_t filterOffset,
					std::size_t inputOffset, const std::vector<FloatType>& input,
					std::size_t& outputOffset, std::vector<FloatType>& output)
{
	const unsigned int m = inputOffset % downsamplingFactor_;
	// The first input relative index that is multiple of downsamplingFactor_, considering the absolute value.
	const unsigned int firstIndex = (m == 0) ? 0 : downsamplingFactor_ - m;

	if (input.size() < 2 * filterOffset + firstIndex + 1) {
		THROW_EXCEPTION(InvalidValueException, "The input is too short.");
	}

	// Copy to the output, compensating for the FIR filter delay. The signal is truncated at both ends.
	output.resize(1 + (input.size() - 2 * filterOffset - firstIndex - 1) / downsamplingFactor_);
	for (std::size_t i = filterOffset + firstIndex, iEnd = input.size() - filterOffset, j = 0;
			i < iEnd;
			i += downsamplingFactor_, ++j) {
		output[j] = input[i];
	}

	outputOffset = (inputOffset + firstIndex) / downsamplingFactor_;
}

template<typename FloatType>
void
Decimator<FloatType>::decimate(std::size_t inputOffset, const std::vector<FloatType>& input,
				std::size_t& outputOffset, std::vector<FloatType>& output)
{
	if (downsamplingFactor_ == 0) THROW_EXCEPTION(InvalidStateException, "The decimator has not been initialized.");

	// Apply the low-pass filter.
	filter_.filter(input, filteredSignal_);

	const std::size_t filterOffset = (lowPassFIRFilter_.size() - 1) / 2;

	downsample(filterOffset, inputOffset, filteredSignal_, outputOffset, output);
}

} // namespace Lab

#endif // DECIMATOR_H
