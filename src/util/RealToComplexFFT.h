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

#ifndef REALTOCOMPLEXFFT_H_
#define REALTOCOMPLEXFFT_H_

#include <cmath> /* min */
#include <complex>
#include <cstddef> /* std::size_t */

#include "Exception.h"
#include "FFTW.h"
#include "Log.h"
#include "Value.h"



namespace Lab {

template<typename TFloat>
class RealToComplexFFT {
public:
	RealToComplexFFT();
	RealToComplexFFT(const RealToComplexFFT& o);
	~RealToComplexFFT();

	RealToComplexFFT& operator=(const RealToComplexFFT& o);

	// destData must point to an allocated sequence of complex numbers, with `size` elements.
	template<typename InputElementType, typename OutputElementType>
		void calculate(InputElementType* origData, unsigned int size, OutputElementType* destData);
private:
	static constexpr unsigned int maxInputSize = 1 << 20; // arbitrary

	RealToComplexFFT(RealToComplexFFT&&) = delete;
	RealToComplexFFT& operator=(RealToComplexFFT&&) = delete;

	void clean();
	void prepare(unsigned int numInputSamples);
	void process();

	bool initialized_;
	unsigned int fftSize_;
	unsigned int freqDataSize_;
	FFTWPlan fftPlan_;
	TFloat* timeData_;
	TFloat (*frequencyData_)[2];  // pointer to array of size two
};



template<typename TFloat>
RealToComplexFFT<TFloat>::RealToComplexFFT()
		: initialized_()
		, fftSize_()
		, freqDataSize_()
		, fftPlan_()
		, timeData_()
		, frequencyData_()
{
}

template<typename TFloat>
RealToComplexFFT<TFloat>::RealToComplexFFT(const RealToComplexFFT& o)
{
	*this = o;
}

template<typename TFloat>
RealToComplexFFT<TFloat>::~RealToComplexFFT()
{
	clean();
}

template<typename TFloat>
RealToComplexFFT<TFloat>&
RealToComplexFFT<TFloat>::operator=(const RealToComplexFFT& o)
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
RealToComplexFFT<TFloat>::clean()
{
	{
		FFTW fftw;
		fftw.destroy_plan(fftPlan_);
	}
	FFTW::free(frequencyData_); frequencyData_ = nullptr;
	FFTW::free(timeData_);      timeData_      = nullptr;
	freqDataSize_ = 0;
	fftSize_ = 0;
	initialized_ = false;
}

template<typename TFloat>
void
RealToComplexFFT<TFloat>::prepare(unsigned int numInputSamples)
{
	//LOG_DEBUG << "[RealToComplexFFT::prepare] numInputSamples: " << numInputSamples;

	clean();

	fftSize_ = numInputSamples;
	freqDataSize_ = fftSize_ / 2 + 1;

	try {
		timeData_      = FFTW::alloc_real<TFloat>(fftSize_);
		frequencyData_ = FFTW::alloc_complex<TFloat>(fftSize_);
		FFTW fftw;
		fftPlan_ = fftw.plan_dft_r2c_1d(fftSize_, timeData_, frequencyData_);
	} catch (std::exception& e) {
		clean();
		THROW_EXCEPTION(Exception, "Error in FFTW preparation: " << e.what() << '.');
	}

	initialized_ = true;
}

template<typename TFloat>
void
RealToComplexFFT<TFloat>::process()
{
	// Time --> frequency.
	FFTW::execute(fftPlan_); // fftw_plan_dft_r2c_1d only returns floor(N/2)+1 complex values

	bool nIsEven = ((fftSize_ & 1) == 0);

	// For "negative" frequencies:
	TFloat (*orig)[2] = frequencyData_ + freqDataSize_ - (nIsEven ? 2 : 1);
	TFloat (*dest)[2] = frequencyData_ + freqDataSize_;
	TFloat (*destEnd)[2] = frequencyData_ + fftSize_;
	while (dest != destEnd) {
		(*dest  )[FFTW::REAL] =  (*orig  )[FFTW::REAL];
		(*dest++)[FFTW::IMAG] = -(*orig--)[FFTW::IMAG];
	}
}

template<typename TFloat>
template<typename InputElementType, typename OutputElementType>
void
RealToComplexFFT<TFloat>::calculate(InputElementType* origData, unsigned int size, OutputElementType* destData)
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
	Value::copySequence(origData, origData + fftSize_, timeData_);

	process();

	// Get the output.
	Value::copySequence(frequencyData_, frequencyData_ + fftSize_, destData);
}

} // namespace Lab

#endif // REALTOCOMPLEXFFT_H_
