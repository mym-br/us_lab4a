#ifndef INTERPOLATOR_H_
#define INTERPOLATOR_H_

#include <cmath> /* ceil, log10 */
#include <vector>

#include <boost/math/special_functions/sinc.hpp>

#include "Exception.h"
#include "FFTWFilter.h"
#include "KaiserWindow.h"
#include "Util.h"

#define INTERPOLATOR_KAISER_TOLERANCE (1.0e-4)



namespace Lab {

// This class is copy constructible and assignable.
template<typename FloatType>
class Interpolator {
public:
	Interpolator();
	~Interpolator();

	// lpFilterTransitionWidth:
	//     half the total transition width
	//     1.0 -> pi radian / sample at the original sampling rate
	void prepare(unsigned int upsamplingFactor, FloatType lpFilterHalfTransitionWidth);

	void interpolate(const FloatType* input, std::size_t inputLength, FloatType* output);
private:
	unsigned int upsamplingFactor_;
	std::vector<FloatType> lowPassFIRFilter_;
	std::vector<FloatType> inputVector_;
	std::vector<FloatType> outputVector_;
	FFTWFilter<FloatType> filter_;
};

template<typename FloatType>
Interpolator<FloatType>::Interpolator()
		: upsamplingFactor_(0)
{
}

template<typename FloatType>
Interpolator<FloatType>::~Interpolator()
{
}

template<typename FloatType>
void
Interpolator<FloatType>::prepare(unsigned int upsamplingFactor, FloatType lpFilterHalfTransitionWidth)
{
	if (upsamplingFactor_ != 0) {
		THROW_EXCEPTION(InvalidStateException, "The object is already configured.");
	}
	if (upsamplingFactor < 2) {
		THROW_EXCEPTION(InvalidValueException, "Invalid upsampling factor: " << upsamplingFactor << ". Must be >= 2.");
	}

	//TODO: limit upsampling factor

	const FloatType tol_dB = -20.0 * std::log10(INTERPOLATOR_KAISER_TOLERANCE);

	const FloatType kaiserBeta = KaiserWindow::getBeta(tol_dB);

	//TODO:(???) Check finalTransitionWidth <= 1.0/upsamplingFactor or lpFilterTransitionWidth <= 0.5
	const FloatType finalTransitionWidth = (lpFilterHalfTransitionWidth * 2) / static_cast<FloatType>(upsamplingFactor);
	const unsigned int windowSize = KaiserWindow::getSize(tol_dB, finalTransitionWidth);

	const FloatType twoUpFactor = 2 * upsamplingFactor;
	const unsigned int finalWindowSize = static_cast<unsigned int>(std::ceil((windowSize - 1) / twoUpFactor) * twoUpFactor) + 1U;

	// Kaiser window.
	std::vector<FloatType> window;
	KaiserWindow::getWindow(finalWindowSize, kaiserBeta, window);

	// Sinc.
	const FloatType numPeriods = (finalWindowSize - 1U) / twoUpFactor;
	std::vector<FloatType> x;
	Util::fillSequenceWithSize(x, -numPeriods, numPeriods, finalWindowSize);
	lowPassFIRFilter_.resize(x.size());
	for (unsigned int i = 0; i < lowPassFIRFilter_.size(); ++i) {
		lowPassFIRFilter_[i] = boost::math::sinc_pi(PI * x[i]);
	}

	Util::multiply(window, lowPassFIRFilter_);

	filter_.setCoefficients(lowPassFIRFilter_);

	upsamplingFactor_ = upsamplingFactor;
}

// The argument "ouput" must point to an array of size inputLength * upsamplingFactor.
template<typename FloatType>
void
Interpolator<FloatType>::interpolate(const FloatType* input, std::size_t inputLength, FloatType* output)
{
	if (upsamplingFactor_ == 0) THROW_EXCEPTION(InvalidStateException, "The interpolator has not been initialized.");

	// Upsamples.
	inputVector_.resize(inputLength * upsamplingFactor_);
	for (std::size_t i = 0, j = 0; i < inputLength; ++i, j += upsamplingFactor_) {
		inputVector_[j] = input[i];
		for (std::size_t k = j + 1, end = j + upsamplingFactor_; k < end; ++k) {
			inputVector_[k] = 0;
		}
	}

	// Applies the anti-aliasing filter.
	filter_.filter(inputVector_, outputVector_);

	// Copies to the output, compensating for the FIR filter delay. The signal is truncated at both ends.
	const std::size_t offset = (lowPassFIRFilter_.size() - 1) / 2;
//	for (std::size_t i = 0; i < inputLength * UPSAMPLING_FACTOR; ++i) {
//		output[i] = outputVector_[i + offset];
//	}
	memcpy(output, &outputVector_[offset], inputLength * upsamplingFactor_ * sizeof(FloatType));
}

} // namespace Lab

#endif // INTERPOLATOR_H_
