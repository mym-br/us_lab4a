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

#ifndef HILBERTENVELOPE_H_
#define HILBERTENVELOPE_H_

#include <cmath> /* max */

#include "Exception.h"
#include "FFTUtil.h"
#include "FFTW.h"
//#include "Log.h"
#include "Util.h"
#include "Value.h"



namespace Lab {

// Calculates the envelope using the Hilbert transform.
template<typename FloatType>
class HilbertEnvelope {
public:
	HilbertEnvelope();
	HilbertEnvelope(const HilbertEnvelope& o);
	~HilbertEnvelope();

	HilbertEnvelope& operator=(const HilbertEnvelope&);

	template<typename ExternalElementType>
		void calculate(ExternalElementType* data, unsigned int size);

	// destData must point to an allocated sequence of complex numbers, with `size` elements.
	template<typename InputElementType, typename OutputElementType>
		void getAnalyticSignal(InputElementType* origData, unsigned int size, OutputElementType* destData);
private:
	// Must be a power of two.
	static constexpr unsigned int minPadding = 1024; /* reduces the aliasing to around -60 dB */
	static constexpr unsigned int maxInputSize = 1 << 20; // arbitrary

	void clean();
	void prepare(unsigned int numInputSamples);
	void process();

	bool initialized_;
	unsigned int numInputSamples_;
	unsigned int fftSize_;
	unsigned int freqDataSize_;
	FloatType* inputTimeData_;
	FloatType (*frequencyData_)[2];  // pointer to array of size two
	FloatType (*outputTimeData_)[2]; // pointer to array of size two
	FFTWPlan fftPlan_;
	FFTWPlan ifftPlan_;
};



template<typename FloatType>
HilbertEnvelope<FloatType>::HilbertEnvelope()
		: initialized_()
		, numInputSamples_()
		, fftSize_()
		, freqDataSize_()
		, inputTimeData_()
		, frequencyData_()
		, outputTimeData_()
		, fftPlan_()
		, ifftPlan_()
{
}

template<typename FloatType>
HilbertEnvelope<FloatType>::HilbertEnvelope(const HilbertEnvelope& o)
		: initialized_()
		, numInputSamples_()
		, fftSize_()
		, freqDataSize_()
		, inputTimeData_()
		, frequencyData_()
		, outputTimeData_()
		, fftPlan_()
		, ifftPlan_()
{
	*this = o;
}

template<typename FloatType>
HilbertEnvelope<FloatType>::~HilbertEnvelope()
{
	clean();
}

template<typename FloatType>
HilbertEnvelope<FloatType>&
HilbertEnvelope<FloatType>::operator=(const HilbertEnvelope& o)
{
	if (&o != this) {
		if (o.initialized_) {
			prepare(o.numInputSamples_);
		} else {
			clean();
		}
	}
	return *this;
}

template<typename FloatType>
void
HilbertEnvelope<FloatType>::clean()
{
	{
		FFTW fftw;
		fftw.destroy_plan(ifftPlan_);
		fftw.destroy_plan(fftPlan_);
	}
	FFTW::free(outputTimeData_); outputTimeData_ = nullptr;
	FFTW::free(frequencyData_);  frequencyData_  = nullptr;
	FFTW::free(inputTimeData_);  inputTimeData_  = nullptr;
	freqDataSize_ = 0;
	fftSize_ = 0;
	numInputSamples_ = 0;
	initialized_ = false;
}

template<typename FloatType>
void
HilbertEnvelope<FloatType>::prepare(unsigned int numInputSamples)
{
	clean();

	// The Hilbert impulse response decays with a factor of 1/n.
	// Padding is added to reduce the aliasing.
	// Warning: some memory is wasted here. The code could consider powers of 3/5/7 too.
	//fftSize_ = std::max<unsigned int>(Util::nextPowerOf2(numInputSamples + minPadding), 2 * minPadding);
	fftSize_ = std::max<unsigned int>(FFTUtil::nextFastEvenSize(numInputSamples + minPadding), 2 * minPadding);
	freqDataSize_ = fftSize_ / 2 + 1;
//	LOG_DEBUG << "[HilbertEnvelope::prepare] numInputSamples: "<< numInputSamples <<
//		", fftSize_: " << fftSize_ <<
//		", freqDataSize_: " << freqDataSize_;

	try {
		inputTimeData_  = FFTW::alloc_real<FloatType>(fftSize_);
		frequencyData_  = FFTW::alloc_complex<FloatType>(fftSize_);
		outputTimeData_ = FFTW::alloc_complex<FloatType>(fftSize_);
		FFTW fftw;
		fftPlan_ = fftw.plan_dft_r2c_1d(fftSize_, inputTimeData_, frequencyData_);
		ifftPlan_ = fftw.plan_idft_1d(fftSize_, frequencyData_, outputTimeData_);
	} catch (std::exception& e) {
		clean();
		THROW_EXCEPTION(Exception, "Error in FFTW preparation: " << e.what() << '.');
	}

	numInputSamples_ = numInputSamples;
	initialized_ = true;
}

template<typename FloatType>
void
HilbertEnvelope<FloatType>::process()
{
	// Time --> frequency.
	FFTW::execute(fftPlan_); // fftw_plan_dft_r2c_1d only returns floor(N/2)+1 complex values

	frequencyData_[0][FFTW::REAL] *= 0.5;
	// Here we are using 0.5 to avoid many multiplications by 2.
	// The results will need to be multiplied by 2 later.
	if ((fftSize_ & 1) == 0) { // when N is even
		frequencyData_[freqDataSize_ - 1][FFTW::REAL] *= 0.5;
	}
	// For "negative" frequencies:
	FloatType (*p)[2] = frequencyData_ + freqDataSize_;
	FloatType (*pEnd)[2] = frequencyData_ + fftSize_;
	while (p != pEnd) {
		(*p)[FFTW::REAL] = 0.0;
		(*p++)[FFTW::IMAG] = 0.0;
	}

	// Frequency --> time.
	FFTW::execute(ifftPlan_);
}

template<typename FloatType>
template<typename ExternalElementType>
void
HilbertEnvelope<FloatType>::calculate(ExternalElementType* data, unsigned int size)
{
	if (size == 0) {
		THROW_EXCEPTION(InvalidParameterException, "The number of input samples (" << size << ") is = 0.");
	}
	if (size > maxInputSize) {
		THROW_EXCEPTION(InvalidParameterException, "The number of input samples (" << size <<
				") is greater than " << maxInputSize << '.');
	}
	if (!initialized_ || size != numInputSamples_) {
		prepare(size);
	}

	// Fills the input vector.
	Value::copySequenceWithPadding(data, data + numInputSamples_, inputTimeData_, fftSize_ - numInputSamples_);

	process();

	// Gets the output.
	const FloatType coef = 2 / static_cast<FloatType>(fftSize_); // The values will be divided by fftSize because FFTW produces unnormalized results
	Value::transformSequence(outputTimeData_, outputTimeData_ + numInputSamples_, data, Value::ComplexToScaledAbsoluteOp<FloatType>(coef));
}



template<typename FloatType>
template<typename InputElementType, typename OutputElementType>
void
HilbertEnvelope<FloatType>::getAnalyticSignal(InputElementType* origData, unsigned int size, OutputElementType* destData)
{
	if (size == 0) {
		THROW_EXCEPTION(InvalidParameterException, "The number of input samples (" << size << ") is = 0.");
	}
	if (size > maxInputSize) {
		THROW_EXCEPTION(InvalidParameterException, "The number of input samples (" << size <<
				") is greater than " << maxInputSize << '.');
	}
	if (!initialized_ || size != numInputSamples_) {
		prepare(size);
	}

	// Fills the input vector.
	Value::copySequenceWithPadding(origData, origData + numInputSamples_, inputTimeData_, fftSize_ - numInputSamples_);

	process();

	// Gets the output.
	const FloatType coef = 2 / static_cast<FloatType>(fftSize_); // The values will be divided by fftSize because FFTW produces unnormalized results
	Value::transformSequence(outputTimeData_, outputTimeData_ + numInputSamples_, destData, Value::ScaleComplexOp<FloatType>(coef));
}

} // namespace Lab

#endif /* HILBERTENVELOPE_H_ */
