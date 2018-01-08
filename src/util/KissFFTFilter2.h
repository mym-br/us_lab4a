#ifndef KISSFFTFILTER2_H_
#define KISSFFTFILTER2_H_

#include <cstddef> /* std::size_t */
#include <vector>

#include "kiss_fftr.h"



namespace Lab {

// Convolver that uses the overlap-save method.
//
class KissFFTFilter2 {
public:
	KissFFTFilter2();
	~KissFFTFilter2();

	void prepare();
	void setCoefficients(const std::vector<double>& filterCoeff, std::vector<kiss_fft_cpx>& filterFreqCoeff);

	// y.size() will be x.size() + filterCoeff.size() - 1.
	void filter(const std::vector<kiss_fft_cpx>& filterFreqCoeff, std::size_t filterSize, const std::vector<double>& x, std::vector<double>& y);//TODO: remove filterSize???
private:
	enum {
		FFT_SIZE = 512 // must be a power of two
	};

	void clean();
	void copyInput(const std::vector<double>& v, int offset);

	bool initialized_;
//	std::size_t filterSize_;
//	std::size_t numFilterBlocks_;
	std::vector<kiss_fft_scalar> fftIn1_;
//	std::vector<kiss_fft_cpx> fftOut1_;
	std::vector<kiss_fft_scalar> fftIn2_;
	std::vector<kiss_fft_cpx> fftOut2_;
	std::vector<kiss_fft_cpx> ifftIn_;
	std::vector<kiss_fft_scalar> ifftOut_;
//	std::vector<kiss_fft_cpx> filterFreqCoeff_;
	kiss_fftr_cfg fftCfg1_;
	kiss_fftr_cfg fftCfg2_;
	kiss_fftr_cfg ifftCfg_;
};

} // namespace Lab

#endif /* KISSFFTFILTER2_H_ */
