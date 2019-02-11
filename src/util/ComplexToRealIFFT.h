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

#ifndef COMPLEXTOREALIFFT_H_
#define COMPLEXTOREALIFFT_H_

#include <cmath> /* min */
#include <complex>
#include <cstddef> /* std::size_t */

#include "Exception.h"
#include "FFTW.h"
#include "Log.h"
#include "Value.h"



namespace Lab {

template<typename FloatType>
class ComplexToRealIFFT {
public:
	ComplexToRealIFFT();
	ComplexToRealIFFT(const ComplexToRealIFFT& o);
	~ComplexToRealIFFT();

	ComplexToRealIFFT& operator=(const ComplexToRealIFFT&);

	// destData must point to an allocated sequence, with `size` elements.
	template<typename InputElementType, typename OutputElementType>
		void calculate(InputElementType* origData, unsigned int size, OutputElementType* destData);
private:
	enum {
		// Must be a power of two.
		MIN_PADDING = 1024, /* reduces the aliasing to around -60 dB */

		MAX_INPUT_SIZE = 1 << 20 // arbitrary
	};

	void clean();
	void prepare(unsigned int numInputSamples);

	bool initialized_;
	unsigned int fftSize_;
	unsigned int freqDataSize_;
	FFTWPlan ifftPlan_;
	FloatType (*frequencyData_)[2];  // pointer to array of size two
	FloatType* timeData_;
};



template<typename FloatType>
ComplexToRealIFFT<FloatType>::ComplexToRealIFFT()
		: initialized_()
		, fftSize_()
		, freqDataSize_()
		, ifftPlan_()
		, frequencyData_()
		, timeData_()
{
}

template<typename FloatType>
ComplexToRealIFFT<FloatType>::ComplexToRealIFFT(const ComplexToRealIFFT& o)
{
	*this = o;
}

template<typename FloatType>
ComplexToRealIFFT<FloatType>::~ComplexToRealIFFT()
{
	clean();
}

template<typename FloatType>
ComplexToRealIFFT<FloatType>&
ComplexToRealIFFT<FloatType>::operator=(const ComplexToRealIFFT& o)
{
	if (&o != this) {
		if (o.initialized_) {
			prepare(o.fftSize_);
		} else {
			clean();
		}
	}
	return *this;
}

template<typename FloatType>
void
ComplexToRealIFFT<FloatType>::clean()
{
	{
		FFTW fftw;
		fftw.destroy_plan(ifftPlan_);
	}
	FFTW::free(timeData_);      timeData_      = nullptr;
	FFTW::free(frequencyData_); frequencyData_ = nullptr;
	freqDataSize_ = 0;
	fftSize_ = 0;
	initialized_ = false;
}

template<typename FloatType>
void
ComplexToRealIFFT<FloatType>::prepare(unsigned int numInputSamples)
{
	//LOG_DEBUG << "[ComplexToRealIFFT::prepare] numInputSamples: " << numInputSamples;

	clean();

	fftSize_ = numInputSamples;
	freqDataSize_ = fftSize_ / 2 + 1;

	try {
		frequencyData_ = FFTW::alloc_complex<FloatType>(freqDataSize_);
		timeData_      = FFTW::alloc_real<FloatType>(fftSize_);
		FFTW fftw;
		ifftPlan_ = fftw.plan_idft_c2r_1d(fftSize_, frequencyData_, timeData_);
	} catch (std::exception& e) {
		clean();
		THROW_EXCEPTION(Exception, "Error in FFTW preparation: " << e.what() << '.');
	}

	initialized_ = true;
}

template<typename FloatType>
template<typename InputElementType, typename OutputElementType>
void
ComplexToRealIFFT<FloatType>::calculate(InputElementType* origData, unsigned int size, OutputElementType* destData)
{
	if (size == 0) {
		THROW_EXCEPTION(InvalidParameterException, "The number of input samples (" << size << ") is = 0.");
	}
	if (size > MAX_INPUT_SIZE) {
		THROW_EXCEPTION(InvalidParameterException, "The number of input samples (" << size <<
				") is greater than " << MAX_INPUT_SIZE << '.');
	}
	if (!initialized_ || size != fftSize_) {
		prepare(size);
	}

	// Fills the input vector.
	// fftw_plan_idft_c2r_1d uses only floor(N/2)+1 complex values.
	Value::copySequence(origData, origData + freqDataSize_, frequencyData_);

	// Frequency --> time.
	FFTW::execute(ifftPlan_);

	// Gets the output.
	const FloatType coef = 1 / static_cast<FloatType>(fftSize_); // The values will be divided by fftSize because FFTW produces unnormalized results
	Value::transformSequence(timeData_, timeData_ + fftSize_, destData, Value::ScaleOp<FloatType>(coef));
}

} // namespace Lab

#endif // COMPLEXTOREALIFFT_H_
