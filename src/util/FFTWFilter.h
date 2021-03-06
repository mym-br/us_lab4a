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
// The FFT size is greater than the size of filterCoeff.
//
// This class is copy constructible and assignable.
template<typename TFloat>
class FFTWFilter {
public:
	FFTWFilter();
	FFTWFilter(const FFTWFilter& o);
	~FFTWFilter();

	FFTWFilter& operator=(const FFTWFilter& o);

	void setCoefficients(const std::vector<TFloat>& filterCoeff);

	// y.size() will be x.size() + filterCoeff.size() - 1.
	void filter(const std::vector<TFloat>& x, std::vector<TFloat>& y);
private:
	static constexpr unsigned int minFFTSize = 512;

	FFTWFilter(FFTWFilter&&) = delete;
	FFTWFilter& operator=(FFTWFilter&&) = delete;

	void clean();
	void copyInput(const std::vector<TFloat>& v, long offset);
	void prepare();

	bool initialized_;
	unsigned int fftSize_;
	unsigned int freqDataSize_;
	unsigned int filterLength_;
	TFloat* fftIn1_;
	TFloat (*fftOut1_)[2]; // pointer to array of size two
	TFloat* fftIn2_;
	TFloat (*fftOut2_)[2]; // pointer to array of size two
	TFloat* ifftOut_;
	FFTWPlan fftPlan1_;
	FFTWPlan fftPlan2_;
	FFTWPlan ifftPlan_;
};



template<typename TFloat>
FFTWFilter<TFloat>::FFTWFilter()
		: initialized_()
		, fftSize_()
		, freqDataSize_()
		, filterLength_()
		, fftIn1_()
		, fftOut1_()
		, fftIn2_()
		, fftOut2_()
		, ifftOut_()
		, fftPlan1_()
		, fftPlan2_()
		, ifftPlan_()
{
}

template<typename TFloat>
FFTWFilter<TFloat>::FFTWFilter(const FFTWFilter& o)
		: initialized_()
		, fftSize_()
		, freqDataSize_()
		, filterLength_()
		, fftIn1_()
		, fftOut1_()
		, fftIn2_()
		, fftOut2_()
		, ifftOut_()
		, fftPlan1_()
		, fftPlan2_()
		, ifftPlan_()
{
	*this = o;
}

template<typename TFloat>
FFTWFilter<TFloat>::~FFTWFilter()
{
	clean();
}

template<typename TFloat>
FFTWFilter<TFloat>&
FFTWFilter<TFloat>::operator=(const FFTWFilter& o)
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

template<typename TFloat>
void
FFTWFilter<TFloat>::clean()
{
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

template<typename TFloat>
void
FFTWFilter<TFloat>::prepare()
{
	assert(!initialized_);
	assert(fftSize_ > 0);
	assert(freqDataSize_ > 0);
	assert(filterLength_ > 0);
	//LOG_DEBUG << "FFTWFilter::prepare()";

	try {
		fftIn1_   = FFTW::alloc_real<TFloat>(fftSize_);
		fftOut1_  = FFTW::alloc_complex<TFloat>(freqDataSize_);
		fftIn2_   = FFTW::alloc_real<TFloat>(fftSize_);
		fftOut2_  = FFTW::alloc_complex<TFloat>(freqDataSize_);
		ifftOut_  = FFTW::alloc_real<TFloat>(fftSize_);
		FFTW fftw;
		fftPlan1_ = fftw.plan_dft_r2c_1d(fftSize_, fftIn1_, fftOut1_);
		fftPlan2_ = fftw.plan_dft_r2c_1d(fftSize_, fftIn2_, fftOut2_);
		ifftPlan_ = fftw.plan_idft_c2r_1d(fftSize_, fftOut2_, ifftOut_);
	} catch (std::exception& e) {
		clean();
		THROW_EXCEPTION(Exception, "Error in FFTW preparation: " << e.what() << '.');
	}
}

template<typename TFloat>
void
FFTWFilter<TFloat>::setCoefficients(const std::vector<TFloat>& filterCoeff)
{
	if (initialized_) THROW_EXCEPTION(InvalidStateException, "The filter is already initialized.");

	const unsigned int filterSize = filterCoeff.size();
	if (filterSize == 0) {
		THROW_EXCEPTION(InvalidValueException, "The filter coefficients vector is empty.");
	}

	//fftSize_ = Util::nextPowerOf2(std::max(filterSize * 2, minFFTSize));
	fftSize_ = FFTUtil::nextFastEvenSize(std::max(filterSize * 2, minFFTSize));
	//LOG_DEBUG << "[FFTWFilter::setCoefficients] fftSize: " << fftSize_ << " filterSize: " << filterSize;
	freqDataSize_ = fftSize_ / 2 + 1;
	filterLength_ = filterSize;

	prepare();

	memcpy(fftIn1_, &filterCoeff[0], sizeof(TFloat) * filterSize);
	memset(fftIn1_ + filterSize, 0, sizeof(TFloat) * (fftSize_ - filterSize));
	FFTW::execute(fftPlan1_);

	initialized_ = true;
}

template<typename TFloat>
void
FFTWFilter<TFloat>::copyInput(const std::vector<TFloat>& v, long offset)
{
	if (offset < 0) {
		const std::size_t dataSize = v.size();
		memset(fftIn2_, 0, sizeof(TFloat) * (fftSize_ / 2));
		if (dataSize < fftSize_ / 2) {
			memcpy(fftIn2_ + fftSize_ / 2           , &v[0], sizeof(TFloat) * dataSize);
			memset(fftIn2_ + fftSize_ / 2 + dataSize,     0, sizeof(TFloat) * (fftSize_ / 2 - dataSize));
		} else {
			memcpy(fftIn2_ + fftSize_ / 2           , &v[0], sizeof(TFloat) * (fftSize_ / 2));
		}
	} else {
		const unsigned long base = offset + 1;
		const unsigned long dataSize = v.size() - base;
		fftIn2_[0] = 0.0;
		if (dataSize < fftSize_ - 1) {
			memcpy(fftIn2_ + 1           , &v[base], sizeof(TFloat) * dataSize);
			memset(fftIn2_ + 1 + dataSize,        0, sizeof(TFloat) * (fftSize_ - 1 - dataSize));
		} else {
			memcpy(fftIn2_ + 1           , &v[base], sizeof(TFloat) * (fftSize_ - 1));
		}
	}
}

template<typename TFloat>
void
FFTWFilter<TFloat>::filter(const std::vector<TFloat>& x, std::vector<TFloat>& y)
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
			const TFloat (*p1)[2] = fftOut1_;
			const TFloat (*p1End)[2] = fftOut1_ + freqDataSize_;
			TFloat (*p2)[2] = fftOut2_;
			while (p1 != p1End) {
				// Multiplication of two complex numbers.
				const TFloat re = (*p1)[0] * (*p2)[0] - (*p1)[1] * (*p2)[1];
				const TFloat im = (*p1)[0] * (*p2)[1] + (*p1)[1] * (*p2)[0];
				(*p2)[0] = re;
				(*p2)[1] = im;
				++p1; ++p2;
			}
		}

		FFTW::execute(ifftPlan_);

		const std::size_t offset = i + halfFFTSize;
		const long numOutElem = std::min<long>(halfFFTSize, static_cast<long>(ySize) - static_cast<long>(offset));
		if (numOutElem > 0) {
			const TFloat* orig = ifftOut_ + halfFFTSize;
			const TFloat* origEnd = orig + numOutElem;
			TFloat* dest = &y[offset];
			while (orig != origEnd) {
				*dest++ = *orig++;
			}
		}
	}

	// Normalization.
	Util::multiply<TFloat>(y, 1 / static_cast<TFloat>(fftSize_));
}

} // namespace Lab

#endif /* FFTWFILTER_H_ */
