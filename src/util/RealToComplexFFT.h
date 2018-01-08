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

template<typename FloatType>
class RealToComplexFFT {
public:
	RealToComplexFFT();
	RealToComplexFFT(const RealToComplexFFT& o);
	~RealToComplexFFT();

	RealToComplexFFT& operator=(const RealToComplexFFT&);

	// destData must point to an allocated sequence of complex numbers, with `size` elements.
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
	void process();

	bool initialized_;
	unsigned int fftSize_;
	unsigned int freqDataSize_;
	FFTWPlan fftPlan_;
	FloatType* timeData_;
	FloatType (*frequencyData_)[2];  // pointer to array of size two
};



template<typename FloatType>
RealToComplexFFT<FloatType>::RealToComplexFFT()
		: initialized_(false)
		, fftSize_(0)
		, freqDataSize_(0)
		, fftPlan_()
		, timeData_(nullptr)
		, frequencyData_(nullptr)
{
}

template<typename FloatType>
RealToComplexFFT<FloatType>::RealToComplexFFT(const RealToComplexFFT& o)
{
	*this = o;
}

template<typename FloatType>
RealToComplexFFT<FloatType>::~RealToComplexFFT()
{
	clean();
}

template<typename FloatType>
RealToComplexFFT<FloatType>&
RealToComplexFFT<FloatType>::operator=(const RealToComplexFFT& o)
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
RealToComplexFFT<FloatType>::clean()
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

template<typename FloatType>
void
RealToComplexFFT<FloatType>::prepare(unsigned int numInputSamples)
{
	//LOG_DEBUG << "[RealToComplexFFT::prepare] numInputSamples: " << numInputSamples;

	clean();

	fftSize_ = numInputSamples;
	freqDataSize_ = fftSize_ / 2 + 1;

	try {
		timeData_      = FFTW::alloc_real<FloatType>(fftSize_);
		frequencyData_ = FFTW::alloc_complex<FloatType>(fftSize_);
		FFTW fftw;
		fftPlan_ = fftw.plan_dft_r2c_1d(fftSize_, timeData_, frequencyData_);
	} catch (std::exception& e) {
		clean();
		THROW_EXCEPTION(Exception, "Error in FFTW preparation: " << e.what() << '.');
	}

	initialized_ = true;
}

template<typename FloatType>
void
RealToComplexFFT<FloatType>::process()
{
	// Time --> frequency.
	FFTW::execute(fftPlan_); // fftw_plan_dft_r2c_1d only returns floor(N/2)+1 complex values

	bool nIsEven = ((fftSize_ & 1) == 0);

	// For "negative" frequencies:
	FloatType (*orig)[2] = frequencyData_ + freqDataSize_ - (nIsEven ? 2 : 1);
	FloatType (*dest)[2] = frequencyData_ + freqDataSize_;
	FloatType (*destEnd)[2] = frequencyData_ + fftSize_;
	while (dest != destEnd) {
		(*dest  )[FFTW::REAL] =  (*orig  )[FFTW::REAL];
		(*dest++)[FFTW::IMAG] = -(*orig--)[FFTW::IMAG];
	}
}

template<typename FloatType>
template<typename InputElementType, typename OutputElementType>
void
RealToComplexFFT<FloatType>::calculate(InputElementType* origData, unsigned int size, OutputElementType* destData)
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
	Value::copySequence(origData, origData + fftSize_, timeData_);

	process();

	// Gets the output.
	Value::copySequence(frequencyData_, frequencyData_ + fftSize_, destData);
}

} // namespace Lab

#endif // REALTOCOMPLEXFFT_H_
