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

#include "KissFFTFilter.h"

#include <complex>
#include <cstring>

#include "_kiss_fft_guts.h"

#include "Exception.h"
#include "Log.h"



//TODO: avoid copies???

namespace Lab {

KissFFTFilter::KissFFTFilter()
		: initialized_(false)
		, filterLength_(0)
		, fftIn1_(FFT_SIZE)
		, fftOut1_(FFT_SIZE / 2 + 1)
		, fftIn2_(FFT_SIZE)
		, fftOut2_(FFT_SIZE / 2 + 1)
		, ifftIn_(FFT_SIZE / 2 + 1)
		, ifftOut_(FFT_SIZE)
		, fftCfg1_(nullptr)
		, fftCfg2_(nullptr)
		, ifftCfg_(nullptr)
{
}

void
KissFFTFilter::prepare()
{
	if (initialized_) return;
	//LOG_DEBUG << "KissFFTFilter::prepare()";

	//TODO: enum for not inverse / inverse
	fftCfg1_ = kiss_fftr_alloc(FFT_SIZE, 0 /* not inverse */, nullptr, nullptr);
	if (fftCfg1_ == nullptr) {
		THROW_EXCEPTION(UnavailableResourceException, "Error in kiss_fftr_alloc.");
	}

	fftCfg2_ = kiss_fftr_alloc(FFT_SIZE, 0 /* not inverse */, nullptr, nullptr);
	if (fftCfg2_ == nullptr) {
		THROW_EXCEPTION(UnavailableResourceException, "Error in kiss_fftr_alloc.");
	}

	ifftCfg_ = kiss_fftr_alloc(FFT_SIZE, 1 /* inverse */, nullptr, nullptr);
	if (ifftCfg_ == nullptr) {
		THROW_EXCEPTION(UnavailableResourceException, "Error in kiss_fftr_alloc.");
	}

	initialized_ = true;
}

void
KissFFTFilter::clean()
{
	free(ifftCfg_);
	free(fftCfg2_);
	free(fftCfg1_);
}

KissFFTFilter::~KissFFTFilter()
{
	clean();
}

inline
void
KissFFTFilter::copyInput(const std::vector<double>& v, int offset)
{
	if (offset < 0) {
		const std::size_t dataSize = v.size();
		memset(&fftIn2_[0], 0, sizeof(double) * (FFT_SIZE / 2));
		if (dataSize < FFT_SIZE / 2) {
			memcpy(&fftIn2_[0] + FFT_SIZE / 2, &v[0], sizeof(double) * dataSize);
			memset(&fftIn2_[0] + FFT_SIZE / 2 + dataSize, 0, sizeof(double) * (FFT_SIZE / 2 - dataSize));
		} else {
			memcpy(&fftIn2_[0] + (FFT_SIZE / 2), &v[0], sizeof(double) * (FFT_SIZE / 2));
		}
	} else {
		const int base = offset + 1;
		const int dataSize = v.size() - base;
		fftIn2_[0] = 0.0f;
		if (dataSize < FFT_SIZE - 1) {
			memcpy(&fftIn2_[0] + 1, &v[base], sizeof(double) * dataSize);
			memset(&fftIn2_[0] + 1 + dataSize, 0, sizeof(double) * (FFT_SIZE - 1 - dataSize));
		} else {
			memcpy(&fftIn2_[0] + 1, &v[base], sizeof(double) * (FFT_SIZE - 1));
		}
	}
}

void
KissFFTFilter::setCoefficients(const std::vector<double>& filterCoeff)
{
	filterLength_ = filterCoeff.size();
	if (filterLength_ == 0) {
		THROW_EXCEPTION(InvalidParameterException, "filterCoeff is empty.");
	}
	if (filterLength_ > FFT_SIZE / 2) {
		THROW_EXCEPTION(InvalidParameterException, "The filter length is > " << (FFT_SIZE / 2) << '.');
	}

	memcpy(&fftIn1_[0], &filterCoeff[0], sizeof(double) * filterLength_);
	memset(&fftIn1_[0] + filterLength_, 0, sizeof(double) * (FFT_SIZE - filterLength_));
	kiss_fftr(fftCfg1_, &fftIn1_[0], &fftOut1_[0]);
}

void
KissFFTFilter::filter(const std::vector<double>& x, std::vector<double>& y)
{
	if (!initialized_) THROW_EXCEPTION(InvalidStateException, "The filter has not been initialized.");
	if (filterLength_ == 0) {
		THROW_EXCEPTION(InvalidStateException, "The filter coefficients have not been set.");
	}

	const std::size_t xLen = x.size();
	if (xLen == 0) {
		THROW_EXCEPTION(InvalidParameterException, "x is empty.");
	}
	const int yLen = filterLength_ + xLen - 1;
	y.assign(yLen, 0.0f);

	// For each x block:
	for (int i = -FFT_SIZE / 2, end = xLen - 1; i < end; i += FFT_SIZE / 2) {
		copyInput(x, i);
		kiss_fftr(fftCfg2_, &fftIn2_[0], &fftOut2_[0]);

		for (int j = 0; j < FFT_SIZE / 2 + 1; ++j) {
			C_MUL(ifftIn_[j], fftOut1_[j], fftOut2_[j]);
		}

		kiss_fftri(ifftCfg_, &ifftIn_[0], &ifftOut_[0]);

		const int offset = i + FFT_SIZE / 2;
		const int jMax = std::min(FFT_SIZE / 2, yLen - offset);
		for (int j = 0; j < jMax; ++j) {
			y[offset + j] += ifftOut_[FFT_SIZE / 2 + j];
		}
	}

	for (std::size_t i = 0, end = y.size(); i < end; ++i) {
		y[i] /= FFT_SIZE;
	}
}

} // namespace Lab
