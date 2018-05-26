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

#include "FFTW.h"



namespace Lab {

std::mutex FFTW::mutex_;

FFTW::FFTW()
{
	mutex_.lock();
}

FFTW::~FFTW()
{
	mutex_.unlock();
}

template<>
float*
FFTW::alloc_real<float>(size_t n) {
	float* p = fftwf_alloc_real(n);
	if (p == nullptr) THROW_EXCEPTION(Exception, "Error in fftw_alloc_real.");
	return p;
}
template<>
double*
FFTW::alloc_real<double>(size_t n) {
	double* p = fftw_alloc_real(n);
	if (p == nullptr) THROW_EXCEPTION(Exception, "Error in fftw_alloc_real.");
	return p;
}

template<>
fftwf_complex*
FFTW::alloc_complex<float>(size_t n) {
	fftwf_complex* p = fftwf_alloc_complex(n);
	if (p == nullptr) THROW_EXCEPTION(Exception, "Error in fftw_alloc_complex.");
	return p;
}
template<>
fftw_complex*
FFTW::alloc_complex<double>(size_t n) {
	fftw_complex* p = fftw_alloc_complex(n);
	if (p == nullptr) THROW_EXCEPTION(Exception, "Error in fftw_alloc_complex.");
	return p;
}

} // namespace Lab
