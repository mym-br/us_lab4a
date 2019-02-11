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

#ifndef FFTWFILTER2_H_
#define FFTWFILTER2_H_

#include <algorithm> /* min */
#include <cassert>
#include <complex>
#include <cstddef> /* std::size_t */
#include <cstring> /* memcpy, memset */
#include <vector>

#include "Exception.h"
#include "FFTW.h"
#include "Util.h"



namespace Lab {

// Convolver that uses the overlap-save method.
//
// The DFT of the filter is stored outside the object.
// The FFT size is fixed.
//
// This class is copy constructible and assignable.
template<typename FloatType>
class FFTWFilter2 {
public:
	FFTWFilter2();
	FFTWFilter2(const FFTWFilter2& o);
	~FFTWFilter2();

	FFTWFilter2& operator=(const FFTWFilter2&);

	void setCoefficients(const std::vector<FloatType>& filterCoeff, std::vector<std::complex<FloatType>>& filterFreqCoeff);

	// y.size() will be x.size() + filterCoeff.size() - 1.
	void filter(const std::vector<std::complex<FloatType>>& filterFreqCoeff, const std::vector<FloatType>& x, std::vector<FloatType>& y);
private:
	enum {
		FFT_SIZE = 512 // must be a power of two
	};

	void clean();
	void copyInput(const std::vector<FloatType>& v, long offset);
	void prepare();

	bool initialized_;
	unsigned int freqDataSize_;
	unsigned int filterLength_;
	unsigned int numFilterBlocks_;
	unsigned int filterFreqCoeffSize_;
	FloatType* fftIn_;
	FloatType (*fftOut_)[2]; // pointer to array of size two
	FloatType* ifftOut_;
	FFTWPlan fftPlan_;
	FFTWPlan ifftPlan_;
};



template<typename FloatType>
FFTWFilter2<FloatType>::FFTWFilter2()
		: initialized_()
		, freqDataSize_()
		, filterLength_()
		, numFilterBlocks_()
		, filterFreqCoeffSize_()
		, fftIn_()
		, fftOut_()
		, ifftOut_()
		, fftPlan_()
		, ifftPlan_()
{
}

template<typename FloatType>
FFTWFilter2<FloatType>::FFTWFilter2(const FFTWFilter2& o)
		: initialized_()
		, freqDataSize_()
		, filterLength_()
		, numFilterBlocks_()
		, filterFreqCoeffSize_()
		, fftIn_()
		, fftOut_()
		, ifftOut_()
		, fftPlan_()
		, ifftPlan_()
{
	*this = o;
}

template<typename FloatType>
FFTWFilter2<FloatType>::~FFTWFilter2()
{
	clean();
}

template<typename FloatType>
FFTWFilter2<FloatType>&
FFTWFilter2<FloatType>::operator=(const FFTWFilter2& o)
{
	if (&o != this) {
		if (o.initialized_) {
			this->freqDataSize_ = o.freqDataSize_;
			this->filterLength_ = o.filterLength_;
			prepare();
			this->numFilterBlocks_ = o.numFilterBlocks_;
			this->filterFreqCoeffSize_ = o.filterFreqCoeffSize_;
			this->initialized_ = true;
		} else {
			clean();
		}
	}
	return *this;
}

template<typename FloatType>
void
FFTWFilter2<FloatType>::clean()
{
	{
		FFTW fftw;
		fftw.destroy_plan(ifftPlan_);
		fftw.destroy_plan(fftPlan_);
	}
	FFTW::free(ifftOut_); ifftOut_ = nullptr;
	FFTW::free(fftOut_);  fftOut_  = nullptr;
	FFTW::free(fftIn_);   fftIn_   = nullptr;
	filterFreqCoeffSize_ = 0;
	numFilterBlocks_ = 0;
	filterLength_ = 0;
	freqDataSize_ = 0;
	initialized_ = false;
}

template<typename FloatType>
void
FFTWFilter2<FloatType>::prepare()
{
	assert(!initialized_);
	assert(freqDataSize_ > 0);
	//LOG_DEBUG << "FFTW2Filter::prepare()";

	try {
		fftIn_    = FFTW::alloc_real<FloatType>(FFT_SIZE);
		fftOut_   = FFTW::alloc_complex<FloatType>(freqDataSize_);
		ifftOut_  = FFTW::alloc_real<FloatType>(FFT_SIZE);
		FFTW fftw;
		fftPlan_  = fftw.plan_dft_r2c_1d(FFT_SIZE, fftIn_, fftOut_);
		ifftPlan_ = fftw.plan_idft_c2r_1d(FFT_SIZE, fftOut_, ifftOut_);
	} catch (std::exception& e) {
		clean();
		THROW_EXCEPTION(Exception, "Error in FFTW preparation: " << e.what() << '.');
	}
}

template<typename FloatType>
void
FFTWFilter2<FloatType>::setCoefficients(const std::vector<FloatType>& filterCoeff, std::vector<std::complex<FloatType>>& filterFreqCoeff)
{
	if (initialized_) THROW_EXCEPTION(InvalidStateException, "The filter is already initialized.");

	const unsigned int filterSize = filterCoeff.size();
	if (filterSize == 0) {
		THROW_EXCEPTION(InvalidValueException, "The filter coefficients vector is empty.");
	}

	freqDataSize_ = FFT_SIZE / 2 + 1;
	filterLength_ = filterSize;

	prepare();

	const unsigned int blockSize = FFT_SIZE / 2;
	numFilterBlocks_ = (filterSize % blockSize != 0) ? (filterSize / blockSize + 1) : (filterSize / blockSize);
	filterFreqCoeffSize_ = numFilterBlocks_ * freqDataSize_;
	filterFreqCoeff.resize(filterFreqCoeffSize_);

	// For each filter block:
	for (unsigned int offset = 0, block = 0; offset < filterSize; offset += blockSize, ++block) {

		const unsigned int remainingDataSize = filterSize - offset;
		if (remainingDataSize < blockSize) {
			std::copy(&filterCoeff[offset], &filterCoeff[offset] + remainingDataSize, fftIn_);
			std::fill(fftIn_ + remainingDataSize, fftIn_ + FFT_SIZE, FloatType(0));
		} else {
			std::copy(&filterCoeff[offset], &filterCoeff[offset] + blockSize, fftIn_);
			std::fill(fftIn_ + blockSize, fftIn_ + FFT_SIZE, FloatType(0));
		}

		FFTW::execute(fftPlan_);

		memcpy(&filterFreqCoeff[block * freqDataSize_], fftOut_, sizeof(FloatType) * 2 * freqDataSize_);
	}

	initialized_ = true;
}

template<typename FloatType>
void
FFTWFilter2<FloatType>::copyInput(const std::vector<FloatType>& v, long offset)
{
	if (offset < 0) {
		const std::size_t dataSize = v.size();
		memset(fftIn_, 0, sizeof(FloatType) * (FFT_SIZE / 2));
		if (dataSize < FFT_SIZE / 2) {
			memcpy(fftIn_ + FFT_SIZE / 2           , &v[0], sizeof(FloatType) * dataSize);
			memset(fftIn_ + FFT_SIZE / 2 + dataSize,     0, sizeof(FloatType) * (FFT_SIZE / 2 - dataSize));
		} else {
			memcpy(fftIn_ + FFT_SIZE / 2           , &v[0], sizeof(FloatType) * (FFT_SIZE / 2));
		}
	} else {
		const std::size_t base = offset + 1;
		const std::size_t dataSize = v.size() - base;
		fftIn_[0] = 0.0;
		if (dataSize < FFT_SIZE - 1) {
			memcpy(fftIn_ + 1           , &v[base], sizeof(FloatType) * dataSize);
			memset(fftIn_ + 1 + dataSize,        0, sizeof(FloatType) * (FFT_SIZE - 1 - dataSize));
		} else {
			memcpy(fftIn_ + 1           , &v[base], sizeof(FloatType) * (FFT_SIZE - 1));
		}
	}
}

template<typename FloatType>
void
FFTWFilter2<FloatType>::filter(const std::vector<std::complex<FloatType>>& filterFreqCoeff, const std::vector<FloatType>& x, std::vector<FloatType>& y)
{
	if (!initialized_) THROW_EXCEPTION(InvalidStateException, "The filter has not been initialized.");

	if (filterFreqCoeff.size() != filterFreqCoeffSize_) {
		THROW_EXCEPTION(InvalidParameterException, "The filter coefficients are invalid (wrong size).");
	}

	const int blockSize = FFT_SIZE / 2;
	const std::size_t xSize = x.size();
	const std::size_t ySize = filterLength_ + xSize - 1;
	y.assign(ySize, 0.0);

	for (unsigned int block = 0; block < numFilterBlocks_; ++block) {

		const std::size_t filterOffset = block * static_cast<std::size_t>(blockSize);
		const std::size_t filterFreqOffset = block * static_cast<std::size_t>(freqDataSize_);

		// For each x block:
		for (long j = -blockSize, end = xSize - 1; j < end; j += blockSize) {
			copyInput(x, j);

			FFTW::execute(fftPlan_);

			{
				const std::complex<FloatType>* p1 = &filterFreqCoeff[filterFreqOffset];
				const std::complex<FloatType>* p1End = p1 + freqDataSize_;
				FloatType (*p2)[2] = fftOut_;
				while (p1 != p1End) {
					// Multiplication of two complex numbers.
					const FloatType re = p1->real() * (*p2)[0] - p1->imag() * (*p2)[1];
					const FloatType im = p1->real() * (*p2)[1] + p1->imag() * (*p2)[0];
					(*p2)[0] = re;
					(*p2)[1] = im;
					++p1; ++p2;
				}
			}

			FFTW::execute(ifftPlan_);

			const std::size_t offset = filterOffset + (j + blockSize);
			const long numOutElem = std::min<long>(blockSize, static_cast<long>(ySize) - static_cast<long>(offset));
			if (numOutElem > 0) {
				const FloatType* orig = ifftOut_ + blockSize;
				const FloatType* origEnd = orig + numOutElem;
				FloatType* dest = &y[offset];
				while (orig != origEnd) {
					*dest++ += *orig++;
				}
			}
		}
	}

	// Normalization.
	Util::multiply<FloatType>(y, 1 / FloatType(FFT_SIZE));
}

} // namespace Lab

#endif /* FFTWFILTER2_H_ */
