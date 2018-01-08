#ifndef FFTWFILTER_H_
#define FFTWFILTER_H_

#include <algorithm> /* max, min */
#include <cassert>
#include <cstddef> /* std::size_t */
#include <cstring> /* memcpy, memset */
#include <vector>

#include "Exception.h"
#include "FFTUtil.h"
#include "FFTW.h"
//#include "Log.h"
#include "Util.h"



namespace Lab {

// Filter (convolver) that uses the overlap-save method.
//
// This class is copy constructible and assignable.
template<typename FloatType>
class FFTWFilter {
public:
	FFTWFilter();
	FFTWFilter(const FFTWFilter& o);
	~FFTWFilter();

	FFTWFilter& operator=(const FFTWFilter&);

	void setCoefficients(const std::vector<FloatType>& filterCoeff);

	// y.size() will be x.size() + filterCoeff.size() - 1.
	void filter(const std::vector<FloatType>& x, std::vector<FloatType>& y);
private:
	enum {
		MIN_FFT_SIZE = 512
	};

	void clean();
	void copyInput(const std::vector<FloatType>& v, long offset);
	void prepare();

	bool initialized_;
	unsigned int fftSize_;
	unsigned int freqDataSize_;
	unsigned int filterLength_;
	FloatType* fftIn1_;
	FloatType (*fftOut1_)[2]; // pointer to array of size two
	FloatType* fftIn2_;
	FloatType (*fftOut2_)[2]; // pointer to array of size two
	FloatType* ifftOut_;
	FFTWPlan fftPlan1_;
	FFTWPlan fftPlan2_;
	FFTWPlan ifftPlan_;
};



template<typename FloatType>
FFTWFilter<FloatType>::FFTWFilter()
		: initialized_(false)
		, fftSize_(0)
		, freqDataSize_(0)
		, filterLength_(0)
		, fftIn1_(nullptr)
		, fftOut1_(nullptr)
		, fftIn2_(nullptr)
		, fftOut2_(nullptr)
		, ifftOut_(nullptr)
		, fftPlan1_()
		, fftPlan2_()
		, ifftPlan_()
{
}

template<typename FloatType>
FFTWFilter<FloatType>::FFTWFilter(const FFTWFilter& o)
		: initialized_(false)
		, fftSize_(0)
		, freqDataSize_(0)
		, filterLength_(0)
		, fftIn1_(nullptr)
		, fftOut1_(nullptr)
		, fftIn2_(nullptr)
		, fftOut2_(nullptr)
		, ifftOut_(nullptr)
		, fftPlan1_()
		, fftPlan2_()
		, ifftPlan_()
{
	*this = o;
}

template<typename FloatType>
FFTWFilter<FloatType>::~FFTWFilter()
{
	clean();
}

template<typename FloatType>
FFTWFilter<FloatType>&
FFTWFilter<FloatType>::operator=(const FFTWFilter& o)
{
	if (&o != this) {
		if (o.initialized_) {
			this->fftSize_ = o.fftSize_;
			this->freqDataSize_ = o.freqDataSize_;
			this->filterLength_ = o.filterLength_;
			prepare();
			memcpy(this->fftOut1_, o.fftOut1_, sizeof(*this->fftOut1_) * this->freqDataSize_);
			this->initialized_ = true;
		} else {
			clean();
		}
	}
	return *this;
}

template<typename FloatType>
void
FFTWFilter<FloatType>::clean()
{
//	LOG_DEBUG << "FFTWFilter::clean()";
	{
		FFTW fftw;
		fftw.destroy_plan(ifftPlan_);
		fftw.destroy_plan(fftPlan2_);
		fftw.destroy_plan(fftPlan1_);
	}
	FFTW::free(ifftOut_); ifftOut_ = nullptr;
	FFTW::free(fftOut2_); fftOut2_ = nullptr;
	FFTW::free(fftIn2_);  fftIn2_  = nullptr;
	FFTW::free(fftOut1_); fftOut1_ = nullptr;
	FFTW::free(fftIn1_);  fftIn1_  = nullptr;
	filterLength_ = 0;
	freqDataSize_ = 0;
	fftSize_ = 0;
	initialized_ = false;
}

template<typename FloatType>
void
FFTWFilter<FloatType>::prepare()
{
	assert(!initialized_);
	assert(fftSize_ > 0);
	assert(freqDataSize_ > 0);
	assert(filterLength_ > 0);
	//LOG_DEBUG << "FFTWFilter::prepare()";

	try {
		fftIn1_   = FFTW::alloc_real<FloatType>(fftSize_);
		fftOut1_  = FFTW::alloc_complex<FloatType>(freqDataSize_);
		fftIn2_   = FFTW::alloc_real<FloatType>(fftSize_);
		fftOut2_  = FFTW::alloc_complex<FloatType>(freqDataSize_);
		ifftOut_  = FFTW::alloc_real<FloatType>(fftSize_);
		FFTW fftw;
		fftPlan1_ = fftw.plan_dft_r2c_1d(fftSize_, fftIn1_, fftOut1_);
		fftPlan2_ = fftw.plan_dft_r2c_1d(fftSize_, fftIn2_, fftOut2_);
		ifftPlan_ = fftw.plan_idft_c2r_1d(fftSize_, fftOut2_, ifftOut_);
	} catch (std::exception& e) {
		clean();
		THROW_EXCEPTION(Exception, "Error in FFTW preparation: " << e.what() << '.');
	}
}

template<typename FloatType>
void
FFTWFilter<FloatType>::setCoefficients(const std::vector<FloatType>& filterCoeff)
{
	if (initialized_) THROW_EXCEPTION(InvalidStateException, "The filter is already initialized.");

	const unsigned int filterSize = filterCoeff.size();
	if (filterSize == 0) {
		THROW_EXCEPTION(InvalidValueException, "The filter coefficients vector is empty.");
	}

	//fftSize_ = Util::nextPowerOf2(std::max(filterSize * 2, static_cast<unsigned int>(MIN_FFT_SIZE)));
	fftSize_ = FFTUtil::nextFastEvenSize(std::max(filterSize * 2, static_cast<unsigned int>(MIN_FFT_SIZE)));
	//LOG_DEBUG << "[FFTWFilter::setCoefficients] fftSize: " << fftSize_ << " filterSize: " << filterSize;
	freqDataSize_ = fftSize_ / 2 + 1;
	filterLength_ = filterSize;

	prepare();

	memcpy(fftIn1_, &filterCoeff[0], sizeof(FloatType) * filterSize);
	memset(fftIn1_ + filterSize, 0, sizeof(FloatType) * (fftSize_ - filterSize));
	FFTW::execute(fftPlan1_);

	initialized_ = true;
}

template<typename FloatType>
void
FFTWFilter<FloatType>::copyInput(const std::vector<FloatType>& v, long offset)
{
	if (offset < 0) {
		const std::size_t dataSize = v.size();
		memset(fftIn2_, 0, sizeof(FloatType) * (fftSize_ / 2));
		if (dataSize < fftSize_ / 2) {
			memcpy(fftIn2_ + fftSize_ / 2           , &v[0], sizeof(FloatType) * dataSize);
			memset(fftIn2_ + fftSize_ / 2 + dataSize,     0, sizeof(FloatType) * (fftSize_ / 2 - dataSize));
		} else {
			memcpy(fftIn2_ + fftSize_ / 2           , &v[0], sizeof(FloatType) * (fftSize_ / 2));
		}
	} else {
		const unsigned long base = offset + 1;
		const unsigned long dataSize = v.size() - base;
		fftIn2_[0] = 0.0;
		if (dataSize < fftSize_ - 1) {
			memcpy(fftIn2_ + 1           , &v[base], sizeof(FloatType) * dataSize);
			memset(fftIn2_ + 1 + dataSize,        0, sizeof(FloatType) * (fftSize_ - 1 - dataSize));
		} else {
			memcpy(fftIn2_ + 1           , &v[base], sizeof(FloatType) * (fftSize_ - 1));
		}
	}
}

template<typename FloatType>
void
FFTWFilter<FloatType>::filter(const std::vector<FloatType>& x, std::vector<FloatType>& y)
{
	if (!initialized_) THROW_EXCEPTION(InvalidStateException, "The filter has not been initialized.");

	const std::size_t xSize = x.size();
	if (xSize == 0) {
		THROW_EXCEPTION(InvalidValueException, "x is empty.");
	}
	const std::size_t ySize = filterLength_ + xSize - 1;
	y.assign(ySize, 0);
	const int halfFFTSize = fftSize_ / 2;

	// For each x block:
	for (long i = -halfFFTSize, end = xSize - 1; i < end; i += halfFFTSize) {
		copyInput(x, i);

		FFTW::execute(fftPlan2_);

		{
			const FloatType (*p1)[2] = fftOut1_;
			const FloatType (*p1End)[2] = fftOut1_ + freqDataSize_;
			FloatType (*p2)[2] = fftOut2_;
			while (p1 != p1End) {
				// Multiplication of two complex numbers.
				const FloatType re = (*p1)[0] * (*p2)[0] - (*p1)[1] * (*p2)[1];
				const FloatType im = (*p1)[0] * (*p2)[1] + (*p1)[1] * (*p2)[0];
				(*p2)[0] = re;
				(*p2)[1] = im;
				++p1; ++p2;
			}
		}

		FFTW::execute(ifftPlan_);

		const std::size_t offset = i + halfFFTSize;
		const long numOutElem = std::min<long>(halfFFTSize, static_cast<long>(ySize) - static_cast<long>(offset));
		if (numOutElem > 0) {
			const FloatType* orig = ifftOut_ + halfFFTSize;
			const FloatType* origEnd = orig + numOutElem;
			FloatType* dest = &y[offset];
			while (orig != origEnd) {
				*dest++ = *orig++;
			}
		}
	}

	// Normalization.
	Util::multiply<FloatType>(y, 1 / static_cast<FloatType>(fftSize_));
}

} // namespace Lab

#endif /* FFTWFILTER_H_ */
