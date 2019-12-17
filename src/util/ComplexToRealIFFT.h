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

template<typename TFloat>
class ComplexToRealIFFT {
public:
	ComplexToRealIFFT();
	ComplexToRealIFFT(const ComplexToRealIFFT& o);
	~ComplexToRealIFFT();

	ComplexToRealIFFT& operator=(const ComplexToRealIFFT& o);

	// destData must point to an allocated sequence, with `size` elements.
	template<typename InputElementType, typename OutputElementType>
		void calculate(InputElementType* origData, unsigned int size, OutputElementType* destData);
private:
	static constexpr unsigned int maxInputSize = 1 << 20; // arbitrary

	ComplexToRealIFFT(ComplexToRealIFFT&&) = delete;
	ComplexToRealIFFT& operator=(ComplexToRealIFFT&&) = delete;

	void clean();
	void prepare(unsigned int numInputSamples);

	bool initialized_;
	unsigned int fftSize_;
	unsigned int freqDataSize_;
	FFTWPlan ifftPlan_;
	TFloat (*frequencyData_)[2];  // pointer to array of size two
	TFloat* timeData_;
};



template<typename TFloat>
ComplexToRealIFFT<TFloat>::ComplexToRealIFFT()
		: initialized_()
		, fftSize_()
		, freqDataSize_()
		, ifftPlan_()
		, frequencyData_()
		, timeData_()
{
}

template<typename TFloat>
ComplexToRealIFFT<TFloat>::ComplexToRealIFFT(const ComplexToRealIFFT& o)
{
	*this = o;
}

template<typename TFloat>
ComplexToRealIFFT<TFloat>::~ComplexToRealIFFT()
{
	clean();
}

template<typename TFloat>
ComplexToRealIFFT<TFloat>&
ComplexToRealIFFT<TFloat>::operator=(const ComplexToRealIFFT& o)
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

template<typename TFloat>
void
ComplexToRealIFFT<TFloat>::clean()
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

template<typename TFloat>
void
ComplexToRealIFFT<TFloat>::prepare(unsigned int numInputSamples)
{
	//LOG_DEBUG << "[ComplexToRealIFFT::prepare] numInputSamples: " << numInputSamples;

	clean();

	fftSize_ = numInputSamples;
	freqDataSize_ = fftSize_ / 2 + 1;

	try {
		frequencyData_ = FFTW::alloc_complex<TFloat>(freqDataSize_);
		timeData_      = FFTW::alloc_real<TFloat>(fftSize_);
		FFTW fftw;
		ifftPlan_ = fftw.plan_idft_c2r_1d(fftSize_, frequencyData_, timeData_);
	} catch (std::exception& e) {
		clean();
		THROW_EXCEPTION(Exception, "Error in FFTW preparation: " << e.what() << '.');
	}

	initialized_ = true;
}

template<typename TFloat>
template<typename InputElementType, typename OutputElementType>
void
ComplexToRealIFFT<TFloat>::calculate(InputElementType* origData, unsigned int size, OutputElementType* destData)
{
	if (size == 0) {
		THROW_EXCEPTION(InvalidParameterException, "The number of input samples (" << size << ") is = 0.");
	}
	if (size > maxInputSize) {
		THROW_EXCEPTION(InvalidParameterException, "The number of input samples (" << size <<
				") is greater than " << maxInputSize << '.');
	}
	if (!initialized_ || size != fftSize_) {
		prepare(size);
	}

	// Fill the input vector.
	// fftw_plan_idft_c2r_1d uses only floor(N/2)+1 complex values.
	Value::copySequence(origData, origData + freqDataSize_, frequencyData_);

	// Frequency --> time.
	FFTW::execute(ifftPlan_);

	// Get the output.
	const TFloat coef = 1 / static_cast<TFloat>(fftSize_); // The values will be divided by fftSize because FFTW produces unnormalized results
	Value::transformSequence(timeData_, timeData_ + fftSize_, destData, Value::ScaleOp<TFloat>(coef));
}

} // namespace Lab

#endif // COMPLEXTOREALIFFT_H_
