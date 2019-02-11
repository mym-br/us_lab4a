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

#ifndef DIRECTFFTWFILTER_H
#define DIRECTFFTWFILTER_H

#include <cassert>
#include <complex>
#include <cstddef> /* std::size_t */
#include <cstring> /* memcpy, memset */
#include <vector>

#include "Exception.h"
#include "FFTUtil.h"
#include "FFTW.h"
//#include "Log.h"
#include "Util.h"



namespace Lab {

// Filter (convolver) that executes DFT -> IDFT only once.
//
// This class is copy constructible and assignable.
template<typename FloatType>
class DirectFFTWFilter {
public:
	DirectFFTWFilter();
	DirectFFTWFilter(const DirectFFTWFilter& o);
	~DirectFFTWFilter();

	DirectFFTWFilter& operator=(const DirectFFTWFilter&);

	void setCoefficients(const std::vector<FloatType>& filterCoeff);

	// y.size() will be x.size() + filterCoeff_.size() - 1.
	void filter(const std::vector<FloatType>& x, std::vector<FloatType>& y);
private:
	enum {
		MAX_DATA_SIZE = 1000000 // arbitrary
	};

	void clean();
	void prepare();
	void prepareFFTW();

	bool initialized_;
	unsigned int dataSize_;
	unsigned int fftSize_;
	unsigned int freqDataSize_;
	FloatType* fftIn1_;
	FloatType (*fftOut1_)[2]; // pointer to array of size two
	FloatType* fftIn2_;
	FloatType (*fftOut2_)[2]; // pointer to array of size two
	FloatType* ifftOut_;
	FFTWPlan fftPlan1_;
	FFTWPlan fftPlan2_;
	FFTWPlan ifftPlan_;
	std::vector<FloatType> filterCoeff_;
};



template<typename FloatType>
DirectFFTWFilter<FloatType>::DirectFFTWFilter()
		: initialized_()
		, dataSize_()
		, fftSize_()
		, freqDataSize_()
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

template<typename FloatType>
DirectFFTWFilter<FloatType>::DirectFFTWFilter(const DirectFFTWFilter& o)
		: initialized_()
		, dataSize_()
		, fftSize_()
		, freqDataSize_()
		, fftIn1_()
		, fftOut1_()
		, fftIn2_()
		, fftOut2_()
		, ifftOut_()
		, fftPlan1_()
		, fftPlan2_()
		, ifftPlan_()
{
	//LOG_DEBUG << "DirectFFTWFilter: copy constructor";

	*this = o;
}

template<typename FloatType>
DirectFFTWFilter<FloatType>::~DirectFFTWFilter()
{
	clean();
}

template<typename FloatType>
DirectFFTWFilter<FloatType>&
DirectFFTWFilter<FloatType>::operator=(const DirectFFTWFilter& o)
{
	//LOG_DEBUG << "DirectFFTWFilter: assignment";

	if (&o != this) {
		if (o.initialized_) {
			this->dataSize_ = o.dataSize_;
			this->fftSize_ = o.fftSize_;
			this->freqDataSize_ = o.freqDataSize_;
			this->filterCoeff_ = o.filterCoeff_;
			prepareFFTW();
			memcpy(this->fftOut1_, o.fftOut1_, sizeof(*this->fftOut1_) * this->freqDataSize_);
			this->initialized_ = true;
		} else if (!o.filterCoeff_.empty()) {
			clean();
			this->filterCoeff_ = o.filterCoeff_;
		} else {
			clean();
		}
	}
	return *this;
}

template<typename FloatType>
void
DirectFFTWFilter<FloatType>::clean()
{
	//LOG_DEBUG << "DirectFFTWFilter::clean()";
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
	freqDataSize_ = 0;
	fftSize_ = 0;
	dataSize_ = 0;
	initialized_ = false;
}

template<typename FloatType>
void
DirectFFTWFilter<FloatType>::prepare()
{
	//LOG_DEBUG << "DirectFFTWFilter::prepare()";
	assert(!filterCoeff_.empty());
	assert(!initialized_);
	assert(dataSize_ > 0);

	const unsigned int filterSize = filterCoeff_.size();

	//fftSize_ = Util::nextPowerOf2(filterSize + dataSize_ - 1);
	fftSize_ = FFTUtil::nextFastEvenSize(filterSize + dataSize_ - 1);
	//LOG_DEBUG << "[DirectFFTWFilter::prepare] fftSize: " << fftSize_ << " filterSize: " << filterSize << " dataSize_: " << dataSize_;
	freqDataSize_ = fftSize_ / 2 + 1;

	prepareFFTW();

	memcpy(fftIn1_, &filterCoeff_[0], sizeof(FloatType) * filterSize);
	if (fftSize_ > filterSize) {
		memset(fftIn1_ + filterSize, 0, sizeof(FloatType) * (fftSize_ - filterSize));
	}
	FFTW::execute(fftPlan1_);

	initialized_ = true;
}

template<typename FloatType>
void
DirectFFTWFilter<FloatType>::prepareFFTW()
{
	//LOG_DEBUG << "DirectFFTWFilter::prepareFFTW()";
	assert(!initialized_);
	assert(fftSize_ > 0);
	assert(freqDataSize_ > 0);

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
DirectFFTWFilter<FloatType>::setCoefficients(const std::vector<FloatType>& filterCoeff)
{
	if (filterCoeff.empty()) THROW_EXCEPTION(InvalidParameterException, "The list of filter coefficients is empty.");
	if (filterCoeff.size() > MAX_DATA_SIZE) THROW_EXCEPTION(InvalidValueException, "The filter length is greater than " << MAX_DATA_SIZE << '.');

	filterCoeff_ = filterCoeff;
	initialized_ = false;
}


template<typename FloatType>
void
DirectFFTWFilter<FloatType>::filter(const std::vector<FloatType>& x, std::vector<FloatType>& y)
{
	if (filterCoeff_.empty()) THROW_EXCEPTION(InvalidStateException, "The filter coefficients have not been set.");

	const std::size_t xSize = x.size();
	if (xSize == 0) THROW_EXCEPTION(InvalidValueException, "x is empty.");
	if (xSize > MAX_DATA_SIZE) THROW_EXCEPTION(InvalidValueException, "The size of x is greater than " << MAX_DATA_SIZE << '.');
	if (!initialized_ || dataSize_ != xSize) {
		dataSize_ = xSize;
		prepare();
	}

	memcpy(fftIn2_, &x[0], sizeof(FloatType) * dataSize_);
	if (fftSize_ > dataSize_) {
		memset(fftIn2_ + dataSize_, 0, sizeof(FloatType) * (fftSize_ - dataSize_));
	}

	FFTW::execute(fftPlan2_);

	const std::complex<FloatType>* p1    = reinterpret_cast<std::complex<FloatType>*>(fftOut1_);
	const std::complex<FloatType>* p1End = reinterpret_cast<std::complex<FloatType>*>(fftOut1_ + freqDataSize_);
	std::complex<FloatType>*       p2    = reinterpret_cast<std::complex<FloatType>*>(fftOut2_);
	for ( ; p1 != p1End; ++p1, ++p2) {
		*p2 = (*p1) * (*p2);
	}

	FFTW::execute(ifftPlan_);

	const std::size_t ySize = filterCoeff_.size() + dataSize_ - 1;
	y.resize(ySize);
	const FloatType k = 1 / static_cast<FloatType>(fftSize_); // normalization factor
	const FloatType* orig = ifftOut_;
	const FloatType* origEnd = ifftOut_ + ySize;
	FloatType* dest = &y[0];
	while (orig != origEnd) {
		*dest++ = *orig++ * k;
	}
}

} // namespace Lab

#endif // DIRECTFFTWFILTER_H
