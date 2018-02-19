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

#ifndef KISSFFTFILTER_H_
#define KISSFFTFILTER_H_

#include <cstddef> /* std::size_t */
#include <vector>

#include "kiss_fftr.h"

#include "Exception.h"



namespace Lab {

// Filter (convolver) that uses the overlap-save method.
//
class KissFFTFilter {
public:
	KissFFTFilter();
	~KissFFTFilter();

	KissFFTFilter(const KissFFTFilter& o) {
		*this = o;
	}

	KissFFTFilter& operator=(const KissFFTFilter& o) {
		if (&o != this) {
			if (o.initialized_) {
				THROW_EXCEPTION(InvalidStateException, "An initialized KissFFTFilter can't be copied.");
			}
			this->initialized_ = false;
			this->filterLength_ = o.filterLength_;
			this->fftIn1_ = o.fftIn1_;
			this->fftOut1_ = o.fftOut1_;
			this->fftIn2_ = o.fftIn2_;
			this->fftOut2_ = o.fftOut2_;
			this->ifftIn_ = o.ifftIn_;
			this->ifftOut_ = o.ifftOut_;
			this->fftCfg1_ = o.fftCfg1_;
			this->fftCfg2_ = o.fftCfg2_;
			this->ifftCfg_ = o.ifftCfg_;
		}
		return *this;
	}

	void prepare();

	// filterCoeff.size() must be <= FFT_SIZE / 2.
	void setCoefficients(const std::vector<double>& filterCoeff);

	// y.size() will be x.size() + filterCoeff.size() - 1.
	void filter(const std::vector<double>& x, std::vector<double>& y);
private:
	enum {
		FFT_SIZE = 512 // must be a power of two
	};

	void clean();
	void copyInput(const std::vector<double>& v, int offset);

	bool initialized_;
	std::size_t filterLength_;
	std::vector<kiss_fft_scalar> fftIn1_;
	std::vector<kiss_fft_cpx> fftOut1_;
	std::vector<kiss_fft_scalar> fftIn2_;
	std::vector<kiss_fft_cpx> fftOut2_;
	std::vector<kiss_fft_cpx> ifftIn_;
	std::vector<kiss_fft_scalar> ifftOut_;
	kiss_fftr_cfg fftCfg1_;
	kiss_fftr_cfg fftCfg2_;
	kiss_fftr_cfg ifftCfg_;
};

} // namespace Lab

#endif /* KISSFFTFILTER_H_ */
