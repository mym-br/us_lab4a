#include "KissFFTFilter2.h"

#include <complex>
#include <cstring>

#include "_kiss_fft_guts.h"

#include "Log.h"
#include "Util.h"

//TODO: avoid copies???

namespace Lab {

KissFFTFilter2::KissFFTFilter2()
		: initialized_(false)
//		, filterSize_(0)
//		, numFilterBlocks_(0)
		, fftIn1_(FFT_SIZE)
//		, fftOut1_(FFT_SIZE / 2 + 1)
		, fftIn2_(FFT_SIZE)
		, fftOut2_(FFT_SIZE / 2 + 1)
		, ifftIn_(FFT_SIZE / 2 + 1)
		, ifftOut_(FFT_SIZE)
//		, filterFreqCoeff_()
		, fftCfg1_(nullptr)
		, fftCfg2_(nullptr)
		, ifftCfg_(nullptr)
{
	//LOG_DEBUG << "KissFFTFilter2()";

	memset(&fftIn1_[0] + (FFT_SIZE / 2), 0, sizeof(double) * (FFT_SIZE / 2));
}

void
KissFFTFilter2::prepare()
{
	if (initialized_) return;
	//LOG_DEBUG << "KissFFTFilter2::prepare()";

	//TODO: enum for not inverse / inverse
	fftCfg1_ = kiss_fftr_alloc(FFT_SIZE, 0 /* not inverse */, nullptr, nullptr);
	if (fftCfg1_ == nullptr) {
		clean();
		THROW_EXCEPTION(UnavailableResourceException, "Error in kiss_fftr_alloc.");
	}

	fftCfg2_ = kiss_fftr_alloc(FFT_SIZE, 0 /* not inverse */, nullptr, nullptr);
	if (fftCfg2_ == nullptr) {
		clean();
		THROW_EXCEPTION(UnavailableResourceException, "Error in kiss_fftr_alloc.");
	}

	ifftCfg_ = kiss_fftr_alloc(FFT_SIZE, 1 /* inverse */, nullptr, nullptr);
	if (ifftCfg_ == nullptr) {
		clean();
		THROW_EXCEPTION(UnavailableResourceException, "Error in kiss_fftr_alloc.");
	}

	initialized_ = true;
}

void
KissFFTFilter2::clean()
{
	free(ifftCfg_);
	free(fftCfg2_);
	free(fftCfg1_);
}

KissFFTFilter2::~KissFFTFilter2()
{
	clean();
}

inline
void
KissFFTFilter2::copyInput(const std::vector<double>& v, int offset)
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
		const std::size_t dataSize = v.size() - base;
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
KissFFTFilter2::setCoefficients(const std::vector<double>& filterCoeff, std::vector<kiss_fft_cpx>& filterFreqCoeff)
{
	const std::size_t blockSize = FFT_SIZE / 2;
	const std::size_t freqBlockSize = FFT_SIZE / 2 + 1;

	const std::size_t filterSize = filterCoeff.size();
	const std::size_t numFilterBlocks = (filterSize % blockSize != 0) ? (filterSize / blockSize + 1UL) : (filterSize / blockSize);
	filterFreqCoeff.resize(numFilterBlocks * freqBlockSize);

	// For each filter block:
	for (std::size_t offset = 0, block = 0; offset < filterSize; offset += FFT_SIZE / 2, ++block) {

		const std::size_t dataSize = filterCoeff.size() - offset;
		if (dataSize < blockSize) {
			memcpy(&fftIn1_[0], &filterCoeff[offset], sizeof(double) * dataSize);
			memset(&fftIn1_[0] + dataSize, 0, sizeof(double) * (blockSize - dataSize));
		} else {
			memcpy(&fftIn1_[0], &filterCoeff[offset], sizeof(double) * blockSize);
		}

		kiss_fftr(fftCfg1_, &fftIn1_[0], &filterFreqCoeff[block * freqBlockSize]);
	}
}

void
KissFFTFilter2::filter(const std::vector<kiss_fft_cpx>& filterFreqCoeff, std::size_t filterSize, const std::vector<double>& x, std::vector<double>& y)
{
	if (filterFreqCoeff.empty()) {
		THROW_EXCEPTION(InvalidParameterException, "filterFreqCoeff is empty.");
	}
	if (filterSize == 0) {
		THROW_EXCEPTION(InvalidParameterException, "filterSize = 0.");
	}

	const std::size_t blockSize = FFT_SIZE / 2;
	const std::size_t freqBlockSize = FFT_SIZE / 2 + 1;
	const std::size_t numFilterBlocks = (filterSize % blockSize != 0) ? (filterSize / blockSize + 1UL) : (filterSize / blockSize);

	const std::size_t xSize = x.size();
	if (xSize == 0) {
		THROW_EXCEPTION(InvalidParameterException, "x is empty.");
	}

	const std::size_t ySize = filterSize + xSize - 1UL;
	y.assign(ySize, 0.0f);

	for (std::size_t block = 0; block < numFilterBlocks; ++block) {

		const int filterOffset = block * blockSize;
		const std::size_t filterFreqOffset = block * freqBlockSize;

		// For each v2 block:
		for (int j = -FFT_SIZE / 2; j < static_cast<int>(xSize) - 1; j += FFT_SIZE / 2) {
			copyInput(x, j);
			kiss_fftr(fftCfg2_, &fftIn2_[0], &fftOut2_[0]);

			for (int k = 0; k < FFT_SIZE / 2 + 1; ++k) {
				C_MUL(ifftIn_[k], filterFreqCoeff[filterFreqOffset + k], fftOut2_[k]);
			}

			kiss_fftri(ifftCfg_, &ifftIn_[0], &ifftOut_[0]);

			const int offset = filterOffset + (j + FFT_SIZE / 2);
			const int kMax = std::min(FFT_SIZE / 2, static_cast<int>(ySize) - offset);
			for (int k = 0; k < kMax; ++k) {
				y[offset + k] += ifftOut_[FFT_SIZE / 2 + k];
			}
		}
	}

	for (std::size_t i = 0; i < ySize; ++i) {
		y[i] /= FFT_SIZE;
	}
}

} // namespace Lab
