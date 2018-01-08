#include "FFTW.h"



namespace Lab {

boost::mutex FFTW::mutex_;

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
