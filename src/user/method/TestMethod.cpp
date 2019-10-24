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

#include "TestMethod.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef> /* std::size_t */
#include <fstream>
#include <vector>

#include "CoherenceFactor.h"
#include "ComplexToRealIFFT.h"
#include "ContainerDumper.h"
#include "Decimator.h"
#include "DirectFFTWFilter.h"
#include "FFTWFilter.h"
#include "FFTWFilter2.h"
#include "HilbertEnvelope.h"
#include "Interpolator.h"
#include "Interpolator4X.h"
#include "KaiserWindow.h"
#include "LinearInterpolator.h"
#include "Log.h"
#include "Matrix.h"
#include "Project.h"
#include "RealToComplexFFT.h"
#include "Statistics.h"
#include "Util.h"

//#define TEST_INTERPOLATOR_SHOW_FIGURES 1
//#define TEST_HILBERT_TRANSFORM_SHOW_FIGURES 1
//#define TEST_FFTW_FILTER_SHOW_FIGURES 1
//#define TEST_DIRECT_FFTW_FILTER_SHOW_FIGURES 1

#define CENTRAL_DIFF_MAX_RELATIVE_ABS_ERROR (1.0e-15)
#define BESSEL_I_MAX_RELATIVE_ABS_ERROR (1.0e-14)
#define FILL_SEQUENCE_MAX_ABS_ERROR (1.0e-15)
#define FFT_MAX_RELATIVE_ABS_ERROR (1.0e-12)
#define FILTER_MAX_RELATIVE_ABS_ERROR (1.0e-14)
#define HILBERT_MAX_RELATIVE_ABS_ERROR (1.0e-15)
#define IFFT_MAX_RELATIVE_ABS_ERROR (1.0e-10)
#define INTERPOLATOR_MAX_RELATIVE_ABS_ERROR (1.0e-14)
#define KAISER_MAX_RELATIVE_ABS_ERROR (1.0e-15)

// Fraction of the original bandwidth.
#define INTERPOLATOR_LP_FILTER_TRANSITION_WIDTH (0.45)

// Fraction of the destination bandwidth.
#define DECIMATOR_LP_FILTER_TRANSITION_WIDTH (0.3)

#define DECIMATOR_OFFSET 32
#define DECIMATOR_MAX_RELATIVE_ABS_ERROR (1.0e-14)



namespace {

unsigned int figureNumber = 0;

void
testCentralDiff(Lab::Project& p)
{
	std::vector<double> x, y, yRef;
	p.loadHDF5("central_diff_x", "v", x);
	p.loadHDF5("central_diff_y", "v", yRef);

	Lab::Util::centralDiff(x, 1.0 / 40e6, y);

	if (y.size() != yRef.size()) {
		THROW_EXCEPTION(Lab::TestException, "Wrong y size.");
	}

	double maxAbs = Lab::Util::maxAbsolute(yRef);
	for (unsigned int i = 0; i < y.size(); ++i) {
		const double error = std::abs(y[i] - yRef[i]) / maxAbs;
		if (error > CENTRAL_DIFF_MAX_RELATIVE_ABS_ERROR) {
			THROW_EXCEPTION(Lab::TestException, "Wrong value for index = " <<
					i << ": " << y[i] << " (error: " << error << ").");
		}
	}

	LOG_INFO << "Test ok: Central diff.";
}

void
testFFT(Lab::Project& p)
{
	std::vector<double> x1, x2, x3;
	std::vector<std::complex<double>> y1, y2, y3;
	std::vector<double> yr1Ref, yi1Ref, yr2Ref, yi2Ref, yr3Ref, yi3Ref;

	p.loadHDF5("fft_x_4000", "v", x1);
	p.loadHDF5("fft_x_2048", "v", x2);
	p.loadHDF5("fft_x_3571", "v", x3);
	p.loadHDF5("fft_yr_4000", "v", yr1Ref);
	p.loadHDF5("fft_yi_4000", "v", yi1Ref);
	p.loadHDF5("fft_yr_2048", "v", yr2Ref);
	p.loadHDF5("fft_yi_2048", "v", yi2Ref);
	p.loadHDF5("fft_yr_3571", "v", yr3Ref);
	p.loadHDF5("fft_yi_3571", "v", yi3Ref);
	y1.resize(x1.size());
	y2.resize(x2.size());
	y3.resize(x3.size());

	Lab::RealToComplexFFT<double> fft;

	fft.calculate(&x1[0], x1.size(), &y1[0]);
	double maxAbs = Lab::Util::maxAbsolute(y1);
	//LOG_DEBUG << "maxAbs: " << maxAbs;
	for (unsigned int i = 0; i < x1.size(); ++i) {
		const double error = std::max(std::abs(y1[i].real() - yr1Ref[i]), std::abs(y1[i].imag() - yi1Ref[i])) / maxAbs;
		if (error > FFT_MAX_RELATIVE_ABS_ERROR) {
			THROW_EXCEPTION(Lab::TestException, "Wrong value for index = " <<
					i << ": (" << y1[i].real() << ", " << y1[i].imag() << ") (error: " << error << ").");
		}
	}

	for (int t = 0; t < 3; ++t) {
		fft.calculate(&x2[0], x2.size(), &y2[0]);
		maxAbs = Lab::Util::maxAbsolute(y2);
		//LOG_DEBUG << "maxAbs: " << maxAbs;
		for (unsigned int i = 0; i < x2.size(); ++i) {
			const double error = std::max(std::abs(y2[i].real() - yr2Ref[i]), std::abs(y2[i].imag() - yi2Ref[i])) / maxAbs;
			if (error > FFT_MAX_RELATIVE_ABS_ERROR) {
				THROW_EXCEPTION(Lab::TestException, "Wrong value for index = " <<
						i << ": (" << y2[i].real() << ", " << y2[i].imag() << ") (error: " << error << ").");
			}
		}
	}

	Lab::RealToComplexFFT<double> fft2;
	fft2.calculate(&x3[0], x3.size(), &y3[0]);
	maxAbs = Lab::Util::maxAbsolute(y3);
	//LOG_DEBUG << "maxAbs: " << maxAbs;
	for (unsigned int i = 0; i < x3.size(); ++i) {
		const double error = std::max(std::abs(y3[i].real() - yr3Ref[i]), std::abs(y3[i].imag() - yi3Ref[i])) / maxAbs;
		if (error > FFT_MAX_RELATIVE_ABS_ERROR) {
			THROW_EXCEPTION(Lab::TestException, "Wrong value for index = " <<
					i << ": (" << y3[i].real() << ", " << y3[i].imag() << ") (error: " << error << ").");
		}
	}

	std::vector<std::complex<double>> yc1(yr1Ref.size());
	std::vector<std::complex<double>> yc2(yr2Ref.size());
	std::vector<std::complex<double>> yc3(yr3Ref.size());
	Lab::Value::copyRealImagToComplexSequence(yr1Ref.begin(), yr1Ref.end(), yi1Ref.begin(), yc1);
	Lab::Value::copyRealImagToComplexSequence(yr2Ref.begin(), yr2Ref.end(), yi2Ref.begin(), yc2);
	Lab::Value::copyRealImagToComplexSequence(yr3Ref.begin(), yr3Ref.end(), yi3Ref.begin(), yc3);
	std::vector<double> x1i(yr1Ref.size());
	std::vector<double> x2i(yr2Ref.size());
	std::vector<double> x3i(yr3Ref.size());

	Lab::ComplexToRealIFFT<double> ifft;

	for (int t = 0; t < 3; ++t) {
		ifft.calculate(&yc1[0], yc1.size(), &x1i[0]);
		maxAbs = Lab::Util::maxAbsolute(x1);
		//LOG_DEBUG << "maxAbs: " << maxAbs;
		for (unsigned int i = 0; i < x1i.size(); ++i) {
			const double error = std::abs(x1[i] - x1i[i]) / maxAbs;
			if (error > IFFT_MAX_RELATIVE_ABS_ERROR) {
				THROW_EXCEPTION(Lab::TestException, "Wrong value for index = " <<
						i << ": " << x1i[i] << " (error: " << error << ").");
			}
		}
	}

	//std::vector<double> idx;
	//Lab::Util::fillSequenceWithSize(idx, 0.0, x1i.size() - 1.0, x1i.size());
	//p.showFigure2D(0, "After IFFT / x1i", idx, x1i);

	Lab::ComplexToRealIFFT<double> ifft2;

	ifft2.calculate(&yc2[0], yc2.size(), &x2i[0]);
	maxAbs = Lab::Util::maxAbsolute(x2);
	//LOG_DEBUG << "maxAbs: "<< maxAbs;
	for (unsigned int i = 0; i < x2i.size(); ++i) {
		const double error = std::abs(x2[i] - x2i[i]) / maxAbs;
		if (error > IFFT_MAX_RELATIVE_ABS_ERROR) {
			THROW_EXCEPTION(Lab::TestException, "Wrong value for index = " <<
					i << ": " << x2i[i] << " (error: " << error << ").");
		}
	}

	//Lab::Util::fillSequenceWithSize(idx, 0.0, x2i.size() - 1.0, x2i.size());
	//p.showFigure2D(1, "After IFFT / x2i", idx, x2i);

	ifft2.calculate(&yc3[0], yc3.size(), &x3i[0]);
	maxAbs = Lab::Util::maxAbsolute(x3);
	//LOG_DEBUG << "maxAbs: " << maxAbs;
	for (unsigned int i = 0; i < x3i.size(); ++i) {
		const double error = std::abs(x3[i] - x3i[i]) / maxAbs;
		if (error > IFFT_MAX_RELATIVE_ABS_ERROR) {
			THROW_EXCEPTION(Lab::TestException, "Wrong value for index = " <<
					i << ": " << x3i[i] << " (error: " << error << ").");
		}
	}

	//Lab::Util::fillSequenceWithSize(idx, 0.0, x3i.size() - 1.0, x3i.size());
	//p.showFigure2D(2, "After IFFT / x3i", idx, x3i);

	LOG_INFO << "Test ok: FFT / IFFT";
}

void
testInterpolator(Lab::Project& p)
{
	const std::vector<unsigned int> upsampFactorList = {2, 8, 64};
	std::vector<double> x, y, y2, yRef;
	p.loadHDF5("interp_source", "v", x);

	for (unsigned int i = 0; i < upsampFactorList.size(); ++i) {
		Lab::Interpolator<double> interp;

		interp.prepare(upsampFactorList[i], INTERPOLATOR_LP_FILTER_TRANSITION_WIDTH);
		y.resize(x.size() * upsampFactorList[i]);
		interp.interpolate(&x[0], x.size(), &y[0]);

		std::ostringstream fileName;
		fileName << "interp_" << upsampFactorList[i] << 'x';
		p.loadHDF5(fileName.str().c_str(), "v", yRef);
#ifdef TEST_INTERPOLATOR_SHOW_FIGURES
		std::vector<double> t;
		Lab::Util::fillSequenceWithSize(t, 0.0, static_cast<double>(y.size() - 1), y.size());
		p.showFigure2D(figureNumber++, "testInterpolator: y", t, y);
#endif
		if (y.size() != yRef.size()) {
			THROW_EXCEPTION(Lab::TestException, "Wrong size: " << y.size() <<
					" (upsampling factor: " << upsampFactorList[i] << ").");
		}
		const double maxAbs = Lab::Util::maxAbsolute(yRef);
		for (unsigned int j = 0; j < y.size(); ++j) {
			const double error = std::abs(y[j] - yRef[j]) / maxAbs;
			if (error > INTERPOLATOR_MAX_RELATIVE_ABS_ERROR) {
				THROW_EXCEPTION(Lab::TestException, "Wrong value for index = " <<
						j << ": " << y[j] << " (error: " << error <<
						", upsampling factor: " << upsampFactorList[i] << ").");
			}
		}

		Lab::Interpolator<double> interp2(interp);
		y2.resize(x.size() * upsampFactorList[i]);
		interp2.interpolate(&x[0], x.size(), &y2[0]);
#ifdef TEST_INTERPOLATOR_SHOW_FIGURES
		p.showFigure2D(figureNumber++, "testInterpolator: y2", t, y2);
#endif
		for (unsigned int j = 0; j < y2.size(); ++j) {
			if (y2[j] != y[j]) {
				THROW_EXCEPTION(Lab::TestException, "y2 != y");
			}
		}
	}

	LOG_INFO << "Test ok: Interpolator";
}

void
testDecimator(Lab::Project& p)
{
	const std::vector<unsigned int> downsampFactorList = {2, 5, 10};
	std::vector<double> x, y, y2, yRef;
	p.loadHDF5("decimation_source", "v", x);

	for (unsigned int i = 0; i < downsampFactorList.size(); ++i) {
		Lab::Decimator<double> decimator;
		std::size_t yOffset;

		decimator.prepare(downsampFactorList[i], DECIMATOR_LP_FILTER_TRANSITION_WIDTH);
		decimator.decimate(DECIMATOR_OFFSET, x, yOffset, y);

		std::ostringstream fileName;
		fileName << "decimated_" << downsampFactorList[i] << 'x';
		p.loadHDF5(fileName.str().c_str(), "v", yRef);

		if (y.size() != yRef.size()) {
			THROW_EXCEPTION(Lab::TestException, "Wrong size: " << y.size() <<
					" (downsampling factor: " << downsampFactorList[i] << ").");
		}
		const double maxAbs = Lab::Util::maxAbsolute(yRef);
		for (unsigned int j = 0; j < y.size(); ++j) {
			if (std::isnan(y[j])) {
				THROW_EXCEPTION(Lab::TestException, "NaN value for y[" <<
						j << "] (downsampling factor: " << downsampFactorList[i] << ").");
			}
			const double error = std::abs(y[j] - yRef[j]) / maxAbs;
			if (error > DECIMATOR_MAX_RELATIVE_ABS_ERROR) {
				THROW_EXCEPTION(Lab::TestException, "Wrong value for index = " <<
						j << ": " << y[j] << " (error: " << error <<
						", downsampling factor: " << downsampFactorList[i] << ").");
			}
		}

		std::ostringstream offsetFileName;
		offsetFileName << p.directory() << "/decimated_" << downsampFactorList[i] << "x-offset.txt";
		std::ifstream in(offsetFileName.str().c_str(), std::ios_base::binary);
		if (!in) THROW_EXCEPTION(Lab::TestException, "Could not open the file: " << offsetFileName.str() << '.');
		unsigned int offsetRef;
		in >> offsetRef;
		if (!in) THROW_EXCEPTION(Lab::TestException, "Could not read the offset from file: " << offsetFileName.str() << '.');
		if (offsetRef != yOffset) {
			THROW_EXCEPTION(Lab::TestException, "Error: yOffset != offsetRef.");
		}

		Lab::Decimator<double> decimator2(decimator);
		std::size_t yOffset2;
		decimator2.decimate(DECIMATOR_OFFSET, x, yOffset2, y2);
		if (y2.size() != y.size()) {
			THROW_EXCEPTION(Lab::TestException, "Error: y2.size() != y.size().");
		}
		for (unsigned int j = 0; j < y2.size(); ++j) {
			if (y2[j] != y[j]) {
				THROW_EXCEPTION(Lab::TestException, "Error: y2 != y (y2[" << j << "] = " << y2[j] <<
						" y[" << j << "] = " << y[j] <<
						" downsampling factor = " << downsampFactorList[i] << ").");
			}
		}
		if (yOffset != yOffset2) {
			THROW_EXCEPTION(Lab::TestException, "Error: yOffset2 != yOffset.");
		}
	}

	LOG_INFO << "Test ok: Decimator";
}

void
testKaiserWindow(Lab::Project& p)
{
	std::vector<double> tol_dB, betaRef;
	p.loadHDF5("kaiser_tol_db", "v", tol_dB);
	p.loadHDF5("kaiser_beta"  , "v", betaRef);
	//p.showFigure2D(0, "Kaiser beta", tol_dB, betaRef);
	for (unsigned int i = 0; i < tol_dB.size(); ++i) {
		const double beta = Lab::KaiserWindow::getBeta(tol_dB[i]);
		const double error = std::abs(beta - betaRef[i]) / betaRef[i];
		if (error > KAISER_MAX_RELATIVE_ABS_ERROR) {
			THROW_EXCEPTION(Lab::TestException, "Wrong beta value for tol_dB = "
							<< tol_dB[i] << ": " << beta << " (error: " << error << ").");
		}
	}

	const std::vector<double> transWidthList = {0.05, 0.1, 0.5};
	std::vector<double> sizeRef;
	for (unsigned int i = 0; i < transWidthList.size(); ++i) {
		const double transWidth = transWidthList[i];
		std::ostringstream fileName;
		fileName << "kaiser_size-trans_width_" << transWidth;
		p.loadHDF5(fileName.str().c_str(), "v", sizeRef);

		for (unsigned int j = 0; j < sizeRef.size(); ++j) {
			const unsigned int size = Lab::KaiserWindow::getSize(tol_dB[i], transWidth);
			const double error = std::abs(size - sizeRef[i]) / sizeRef[i];
			if (error > 0.0) {
				THROW_EXCEPTION(Lab::TestException, "Wrong size value for tol_dB = "
								<< tol_dB[i] << ": " << size << " (error: " << error << ").");
			}
		}
	}

	LOG_INFO << "Test ok: Kaiser window";
}

void
testBessel(Lab::Project& p)
{
	std::vector<double> xBesselI0, besselI0;
	p.loadHDF5("bessel_i0_x", "v", xBesselI0);
	p.loadHDF5("bessel_i0"  , "v", besselI0);

	for (unsigned int i = 0; i < xBesselI0.size(); ++i) {
		const double x = xBesselI0[i];
		const double f = std::cyl_bessel_i(0, x);
		const double error = std::abs(f - besselI0[i]) / besselI0[i];
		if (error > BESSEL_I_MAX_RELATIVE_ABS_ERROR) {
			THROW_EXCEPTION(Lab::TestException, "Wrong value for x = " << x << ": " << f << " (error: " << error << ").");
		}
	}

	LOG_INFO << "Test ok: Bessel";
}

void
testHilbertTransform(Lab::Project& p)
{
	std::vector<double> x, yaRef, yrRef, yiRef;
	p.loadHDF5("hilbert_x", "v", x);
	p.loadHDF5("hilbert_ya", "v", yaRef);
	p.loadHDF5("hilbert_yr", "v", yrRef);
	p.loadHDF5("hilbert_yi", "v", yiRef);
	std::vector<double> ya = x;
	std::vector<std::complex<double>> yc(x.size());

	Lab::HilbertEnvelope<double> e;
	e.calculate(&ya[0], ya.size());

#ifdef TEST_HILBERT_TRANSFORM_SHOW_FIGURES
	std::vector<double> idx;
	Lab::Util::fillSequenceWithSize(idx, 0.0, x.size() - 1.0, x.size());
	p.showFigure2D(figureNumber++, "testHilbertTransform: ya", idx, ya);

	p.showFigure2D(figureNumber++, "testHilbertTransform: yaRef", idx, yaRef);
#endif

	if (ya.size() != yaRef.size()) {
		THROW_EXCEPTION(Lab::TestException, "ya.size() != yaRef.size() (ya.size() = " << ya.size() << ").");
	}
	double maxAbs = Lab::Util::maxAbsolute(yaRef);
	for (unsigned int i = 0; i < ya.size(); ++i) {
		const double error = std::abs(ya[i] - yaRef[i]) / maxAbs;
		if (error > HILBERT_MAX_RELATIVE_ABS_ERROR) {
			THROW_EXCEPTION(Lab::TestException, "Wrong value for i = " << i << ": " << ya[i] << " (error: " << error << ").");
		}
	}

	e.getAnalyticSignal(&x[0], x.size(), &yc[0]);
	if (yc.size() != yrRef.size()) {
		THROW_EXCEPTION(Lab::TestException, "yc.size() != yrRef.size() (yc.size() = " << yc.size() << ").");
	}
	maxAbs = Lab::Util::maxAbsolute(yc);
	for (unsigned int i = 0; i < yc.size(); ++i) {
		const double error = std::max(std::abs(yc[i].real() - yrRef[i]), std::abs(yc[i].imag() - yiRef[i])) / maxAbs;
		if (error > HILBERT_MAX_RELATIVE_ABS_ERROR) {
			THROW_EXCEPTION(Lab::TestException, "Wrong value for i = " << i << ": " << yc[i] << " (error: " << error << ").");
		}
	}

	Lab::HilbertEnvelope<double> e2(e);
	std::vector<double> ya2 = x;
	e2.calculate(&ya2[0], ya2.size());
	for (std::size_t i = 0; i < ya2.size(); ++i) {
		if (ya2[i] != ya[i]) {
			THROW_EXCEPTION(Lab::TestException, "ya2 != ya");
		}
	}

	std::vector<std::complex<double>> yc2(x.size());
	e2.getAnalyticSignal(&x[0], x.size(), &yc2[0]);
	for (std::size_t i = 0; i < yc2.size(); ++i) {
		if (yc2[i] != yc[i]) {
			THROW_EXCEPTION(Lab::TestException, "yc2 != yc");
		}
	}

	LOG_DEBUG << "Test ok: Hilbert transform";
}

void
testStatistics()
{
	const std::vector<double> a = {-1.0, 2.0, 4.0, 7.0};

	const double stdDev = Lab::Statistics::standardDeviation(a.data(), a.size());
	if (std::abs(stdDev - 2.9154759474) > 1e-10) {
		THROW_EXCEPTION(Lab::TestException, "Wrong stdDev: " << stdDev << '.');
	}

	const double mean = Lab::Statistics::arithmeticMean(a.data(), a.size());
	if (mean != 3.0) {
		THROW_EXCEPTION(Lab::TestException, "Wrong mean: " << mean << '.');
	}
	LOG_INFO << "Test ok: Statistics";
}

void
testLinearInterpolator()
{
	const unsigned int upsamplingFactor = 4;
	Lab::LinearInterpolator<double, upsamplingFactor> interp;
	std::vector<double> x;
	x.push_back(1);
	x.push_back(-1);
	x.push_back(0);
	x.push_back(0.5);
	std::vector<double> y(x.size() * upsamplingFactor);
	interp.interpolate(&x[0], 4, &y[0]);

	if (
			y[0]  !=  1.0 || y[1]  !=  0.5   || y[2]  !=  0.0  || y[3]  != -0.5   ||
			y[4]  != -1.0 || y[5]  != -0.75  || y[6]  != -0.5  || y[7]  != -0.25  ||
			y[8]  !=  0.0 || y[9]  !=  0.125 || y[10] !=  0.25 || y[11] !=  0.375 ||
			y[12] !=  0.5 || y[13] !=  0.375 || y[14] !=  0.25 || y[15] !=  0.125) {
		//THROW_EXCEPTION(Lab::TestException, "Wrong values: " << y << '.');
		THROW_EXCEPTION(Lab::TestException, "[testLinearInterpolator] Wrong values.");
	}

	LOG_INFO << "Test ok: LinearInterpolator";
}

void
testFFTWFilter(Lab::Project& p)
{
	std::vector<double> x;
	p.loadHDF5("filter_x", "x", x);
	std::vector<double> b;
	p.loadHDF5("filter_b", "b", b);
	std::vector<double> yRef;
	p.loadHDF5("filter_y", "y", yRef);

	std::vector<double> y;
	Lab::FFTWFilter<double> filter;
	filter.setCoefficients(b);
	filter.filter(x, y);
#ifdef TEST_FFTW_FILTER_SHOW_FIGURES
	std::vector<double> t;
	Lab::Util::fillSequenceWithSize(t, 0.0, static_cast<double>(y.size() - 1), y.size());
	p.showFigure2D(figureNumber++, "testFFTWFilter: y", t, y);
#endif
	if (y.size() != x.size() + b.size() - 1) {
		THROW_EXCEPTION(Lab::TestException, "y.size() != x.size() + b.size() - 1 [y.size()=" << y.size() <<
			" x.size()=" << x.size() <<
			" b.size()=" << b.size() << "].");
	}
	double maxError = 0.0, maxY = 0.0;
	for (std::vector<double>::iterator iter1 = yRef.begin(), iter2 = y.begin(); iter1 != yRef.end(); ++iter1, ++iter2) {
		const double error = std::abs(*iter1 - *iter2);
		if (error > maxError) maxError = error;
		if (std::abs(*iter1) > maxY) maxY = std::abs(*iter1);
	}
	const double maxRelError = maxError / maxY;
	//LOG_DEBUG << "maxError = " << maxError << " maxY = " << maxY << " maxRelError = " << maxRelError;
	if (maxRelError > FILTER_MAX_RELATIVE_ABS_ERROR) {
		THROW_EXCEPTION(Lab::TestException, "The maximum relative error is > " << FILTER_MAX_RELATIVE_ABS_ERROR << '.');
	}

	Lab::FFTWFilter<double> filter2(filter);
//	Lab::FFTWFilter<double> filter2;
//	filter2 = filter;
	std::vector<double> y2;
	filter2.filter(x, y2);
	if (y2.size() != y.size()) {
		THROW_EXCEPTION(Lab::TestException, "y2.size() != y.size()");
	}
#ifdef TEST_FFTW_FILTER_SHOW_FIGURES
	p.showFigure2D(figureNumber++, "testFFTWFilter: y2", t, y2);
#endif
	for (std::size_t i = 0; i < y2.size(); ++i) {
		if (y2[i] != y[i]) {
			THROW_EXCEPTION(Lab::TestException, "y2 != y" << y2[i] - y[i]);
		}
	}

	LOG_INFO << "Test ok: FFTWFilter";
}

void
testDirectFFTWFilter(Lab::Project& p)
{
	std::vector<double> x;
	p.loadHDF5("filter_x", "x", x);
	std::vector<double> b;
	p.loadHDF5("filter_b", "b", b);
	std::vector<double> yRef;
	p.loadHDF5("filter_y", "y", yRef);

	std::vector<double> y;
	Lab::DirectFFTWFilter<double> filter;
	filter.setCoefficients(b);
	filter.filter(x, y);
#ifdef TEST_DIRECT_FFTW_FILTER_SHOW_FIGURES
	std::vector<double> t;
	Lab::Util::fillSequenceWithSize(t, 0.0, static_cast<double>(y.size() - 1), y.size());
	p.showFigure2D(figureNumber++, "testDirectFFTWFilter: y", t, y);
#endif
	if (y.size() != x.size() + b.size() - 1) {
		THROW_EXCEPTION(Lab::TestException, "y.size() != x.size() + b.size() - 1 [y.size()=" << y.size() <<
			" x.size()=" << x.size() <<
			" b.size()=" << b.size() << "].");
	}
	double maxError = 0.0, maxY = 0.0;
	for (std::vector<double>::iterator iter1 = yRef.begin(), iter2 = y.begin(); iter1 != yRef.end(); ++iter1, ++iter2) {
		const double error = std::abs(*iter1 - *iter2);
		if (error > maxError) maxError = error;
		if (std::abs(*iter1) > maxY) maxY = std::abs(*iter1);
	}
	const double maxRelError = maxError / maxY;
	//LOG_DEBUG << "maxError = " << maxError << " maxY = " << maxY << " maxRelError = " << maxRelError;
	if (maxRelError > FILTER_MAX_RELATIVE_ABS_ERROR) {
		THROW_EXCEPTION(Lab::TestException, "The maximum relative error is > " << FILTER_MAX_RELATIVE_ABS_ERROR << '.');
	}

	Lab::DirectFFTWFilter<double> filter2(filter);
//	Lab::DirectFFTWFilter<double> filter2;
//	filter2 = filter;
	std::vector<double> y2;
	filter2.filter(x, y2);
	if (y2.size() != y.size()) {
		THROW_EXCEPTION(Lab::TestException, "y2.size() != y.size()");
	}
#ifdef TEST_DIRECT_FFTW_FILTER_SHOW_FIGURES
	p.showFigure2D(figureNumber++, "testDirectFFTWFilter: y2", t, y2);
#endif
	for (std::size_t i = 0; i < y2.size(); ++i) {
		if (y2[i] != y[i]) {
			THROW_EXCEPTION(Lab::TestException, "y2 != y" << y2[i] - y[i]);
		}
	}

	LOG_INFO << "Test ok: DirectFFTWFilter";
}

void
testFFTWFilter2(Lab::Project& p)
{
	std::vector<double> x;
	p.loadHDF5("filter_x", "x", x);
	std::vector<double> b;
	p.loadHDF5("filter_b", "b", b);
	std::vector<double> yRef;
	p.loadHDF5("filter_y", "y", yRef);
	std::vector<double> y2Ref;
	p.loadHDF5("filter_y2", "y", y2Ref);

	std::vector<std::complex<double>> filterFreqCoeff;
	std::vector<double> y;
	Lab::FFTWFilter2<double> filter;
	filter.setCoefficients(b, filterFreqCoeff);
	filter.filter(filterFreqCoeff, x, y);
	{
		if (y.size() != x.size() + b.size() - 1) {
			THROW_EXCEPTION(Lab::TestException, "y.size() != x.size() + b.size() - 1 [y.size()=" << y.size() <<
				" x.size()=" << x.size() <<
				" b.size()=" << b.size() << "].");
		}
		double maxError = 0.0, maxY = 0.0;
		for (std::vector<double>::iterator iter1 = yRef.begin(), iter2 = y.begin(); iter1 != yRef.end(); ++iter1, ++iter2) {
			const double error = std::abs(*iter1 - *iter2);
			if (error > maxError) maxError = error;
			if (std::abs(*iter1) > maxY) maxY = std::abs(*iter1);
		}
		const double maxRelError = maxError / maxY;
		//LOG_DEBUG << "maxError = " << maxError << " maxY = " << maxY << " maxRelError = " << maxRelError;
		if (maxRelError > FILTER_MAX_RELATIVE_ABS_ERROR) {
			THROW_EXCEPTION(Lab::TestException, "The maximum relative error is > " << FILTER_MAX_RELATIVE_ABS_ERROR << '.');
		}
	}

	Lab::FFTWFilter2<double> filter2;
	filter2.setCoefficients(yRef, filterFreqCoeff);
	Lab::FFTWFilter2<double> filter3;
	std::vector<std::complex<double>> dummyFilterFreqCoeff;
	filter3.setCoefficients(yRef, dummyFilterFreqCoeff);
	std::vector<double> y2;
	filter3.filter(filterFreqCoeff, x, y2);
	{
		if (y2.size() != x.size() + y.size() - 1) {
			THROW_EXCEPTION(Lab::TestException, "y2.size() != x.size() + y.size() - 1 [y2.size()=" << y2.size() <<
				" x.size()=" << x.size() <<
				" y.size()=" << y.size() << "].");
		}
		double maxError = 0.0, maxY = 0.0;
		for (std::vector<double>::iterator iter1 = y2Ref.begin(), iter2 = y2.begin(); iter1 != y2Ref.end(); ++iter1, ++iter2) {
			const double error = std::abs(*iter1 - *iter2);
			if (error > maxError) maxError = error;
			if (std::abs(*iter1) > maxY) maxY = std::abs(*iter1);
		}
		const double maxRelError = maxError / maxY;
		//LOG_DEBUG << "maxError = " << maxError << " maxY = " << maxY << " maxRelError = " << maxRelError;
		if (maxRelError > FILTER_MAX_RELATIVE_ABS_ERROR) {
			THROW_EXCEPTION(Lab::TestException, "The maximum relative error is > " << FILTER_MAX_RELATIVE_ABS_ERROR << '.');
		}
	}

	LOG_INFO << "Test ok: FFTWFilter2";
}

void
testInterpolator4X(Lab::Project& p)
{
	std::vector<double> x;
	p.loadHDF5("interp_x", "x", x);

	std::vector<double> y(x.size() * 4);

	Lab::Interpolator4X<double> interpolator;
	interpolator.interpolate(&x[0], x.size(), &y[0]);
	p.saveHDF5(y, "interp_y", "y");

	LOG_INFO << "Test ok: Interpolator4X";
}

void
testMultiplyBy()
{
	Lab::Matrix<double> m(3, 4);
	m(1, 0) = 1.0;
	m(1, 1) = 1.5;
	m(1, 2) = 2.0;
	m(1, 3) = 8.0;
	Lab::Matrix<double>::Dim2Interval interval = m.dim2Interval(1);
	std::for_each(interval.first, interval.second, Lab::Util::MultiplyBy<double>(-0.5));
	if (m(1, 0) != -0.5 ||
			m(1, 1) != -0.75 ||
			m(1, 2) != -1.0 ||
			m(1, 3) != -4.0) {
		THROW_EXCEPTION(Lab::TestException, "Wrong results.");
	}
	LOG_INFO << "Test ok: MultiplyBy";
}

void
testAdd()
{
	Lab::Matrix<double> m(3, 4);
	m(1, 0) = 1.0;
	m(1, 1) = 1.5;
	m(1, 2) = 2.0;
	m(1, 3) = 18.0;
	Lab::Matrix<double>::Dim2Interval interval = m.dim2Interval(1);
	std::for_each(interval.first, interval.second, Lab::Util::Add<double>(-10.0));
	if (m(1, 0) != -9.0 ||
			m(1, 1) != -8.5 ||
			m(1, 2) != -8 ||
			m(1, 3) != 8.0) {
		THROW_EXCEPTION(Lab::TestException, "Wrong results.");
	}
	LOG_INFO << "Test ok: Add";
}

void
testAddElements()
{
	std::vector<double> v1(5);
	v1[0] = 10.0;
	v1[1] = 11.0;
	v1[2] = 12.0;
	v1[3] = 13.0;
	v1[4] = 14.0;
	std::vector<double> v2(4);
	Lab::Util::addElements(v1.begin(), v2.begin(), v2.end());
	if (v2[0] != 10.0 ||
			v2[1] != 11.0 ||
			v2[2] != 12.0 ||
			v2[3] != 13.0) {
		THROW_EXCEPTION(Lab::TestException, "Wrong results.");
	}
	LOG_INFO << "Test ok: addElements(InputIterator first1, InputOutputIterator first2, InputOutputIterator last2)";
}

void
testAddElements2()
{
	std::vector<double> v1(4);
	v1[0] = 100.0;
	v1[1] = 101.0;
	v1[2] = 102.0;
	v1[3] = 103.0;
	std::vector<double> v2(4);
	v2[0] = 20.0;
	v2[1] = 30.0;
	v2[2] = 40.0;
	v2[3] = 50.0;
	std::vector<double> v3(3);
	Lab::Util::addElements(v1.begin(), v2.begin(), v3.begin(), v3.end());
	if (v3[0] != 120.0 ||
			v3[1] != 131.0 ||
			v3[2] != 142.0) {
		THROW_EXCEPTION(Lab::TestException, "Wrong results.");
	}
	LOG_INFO << "Test ok: addElements(InputIterator1 first1, InputIterator2 first2, OutputIterator first3, OutputIterator last3)";
}

void
testSCF()
{
	std::vector<double> v1;
	v1.push_back(1.0);
	v1.push_back(0.0);
	v1.push_back(-1.0);
	v1.push_back(1.0);
	v1.push_back(-1.0);
	v1.push_back(1.0);

	auto scf = std::make_unique<Lab::SignCoherenceFactor<double>>(2.0);

	const double scfValue = scf->calculate(&v1[0], v1.size());
	if (scfValue != 3.270805724762158e-3) {
		THROW_EXCEPTION(Lab::TestException, "Wrong result.");
	}
	LOG_INFO << "Test ok: calculateSCF";
}

void
testFillSequence()
{
	std::vector<double> v;
	Lab::Util::fillSequenceFromStartToEndWithMaximumStep(v, 1.0, 5.0, 0.99);
//	for (std::size_t i = 0; i < v.size(); ++i) {
//		LOG_DEBUG << "v " << v[i];
//	}
	if (v.size() != 6) {
		THROW_EXCEPTION(Lab::TestException, "Wrong size.");
	}
	if (std::abs(v[0] - 1.0) > FILL_SEQUENCE_MAX_ABS_ERROR ||
			std::abs(v[1] - 1.8) > FILL_SEQUENCE_MAX_ABS_ERROR ||
			std::abs(v[2] - 2.6) > FILL_SEQUENCE_MAX_ABS_ERROR ||
			std::abs(v[3] - 3.4) > FILL_SEQUENCE_MAX_ABS_ERROR ||
			std::abs(v[4] - 4.2) > FILL_SEQUENCE_MAX_ABS_ERROR ||
			std::abs(v[5] - 5.0) > FILL_SEQUENCE_MAX_ABS_ERROR) {
		THROW_EXCEPTION(Lab::TestException, "Wrong results.");
	}
	LOG_INFO << "Test ok: fillSequence(std::vector<T>& v, T startValue, T endValue, T step)";
}

void
testFillSequence2()
{
	std::vector<double> v;
	Lab::Util::fillSequenceFromStartToEndWithSize(v, 1.0, 5.0, 6);
	//for (std::size_t i = 0; i < v.size(); ++i) {
	//	LOG_DEBUG << "v " << v[i];
	//}
	if (v.size() != 6) {
		THROW_EXCEPTION(Lab::TestException, "Wrong size.");
	}
	if (std::abs(v[0] - 1.0) > FILL_SEQUENCE_MAX_ABS_ERROR ||
			std::abs(v[1] - 1.8) > FILL_SEQUENCE_MAX_ABS_ERROR ||
			std::abs(v[2] - 2.6) > FILL_SEQUENCE_MAX_ABS_ERROR ||
			std::abs(v[3] - 3.4) > FILL_SEQUENCE_MAX_ABS_ERROR ||
			std::abs(v[4] - 4.2) > FILL_SEQUENCE_MAX_ABS_ERROR ||
			std::abs(v[5] - 5.0) > FILL_SEQUENCE_MAX_ABS_ERROR) {
		THROW_EXCEPTION(Lab::TestException, "Wrong results.");
	}
	LOG_INFO << "Test ok: fillSequence(std::vector<T>& v, T startValue, T endValue, unsigned int size)";
}

void
testFillSequence3()
{
	{
		std::vector<int> v;

		std::vector<int> vRef{2, 3, 4, 5, 6};
		Lab::Util::fillSequence(v, 2, 6, 1);
		if (v != vRef) {
			THROW_EXCEPTION(Lab::TestException, "[1] v != vRef.");
		}
	}
	{

		std::vector<int> v;
		std::vector<int> vRef{-2, -3, -4, -5};
		Lab::Util::fillSequence(v, -2, -5, -1);
		if (v != vRef) {
			THROW_EXCEPTION(Lab::TestException, "[2] v != vRef.");
		}
	}

	LOG_INFO << "Test ok: fillSequence3(std::vector<int>& v, int startValue, int endValue, int step)";
}

void
testFillSequence4()
{
	{
		std::vector<unsigned int> v;

		std::vector<unsigned int> vRef{2, 3, 4, 5, 6, 7};
		Lab::Util::fillSequence(v, 2U, 7U, 1);
		if (v != vRef) {
			THROW_EXCEPTION(Lab::TestException, "[1] v != vRef.");
		}
	}
	{
		std::vector<unsigned int> v;

		std::vector<unsigned int> vRef{10, 9 , 8, 7, 6};
		Lab::Util::fillSequence(v, 10U, 6U, -1);
		if (v != vRef) {
			THROW_EXCEPTION(Lab::TestException, "[2] v != vRef.");
		}
	}

	LOG_INFO << "Test ok: fillSequence4(std::vector<unsigned int>& v, unsigned int startValue, unsigned int endValue, int step)";
}

} // namespace

namespace Lab {

TestMethod::TestMethod(Project& project)
		: project_(project)
{
}

TestMethod::~TestMethod()
{
}

void
TestMethod::execute()
{
	figureNumber = 0;

	testCentralDiff(project_);
	testFFT(project_);
	testInterpolator(project_);
	testDecimator(project_);
	testKaiserWindow(project_);
	testBessel(project_);
	testHilbertTransform(project_);
	testStatistics();
	testLinearInterpolator();
	testFFTWFilter(project_);
	testDirectFFTWFilter(project_);
	testFFTWFilter2(project_);
	testInterpolator4X(project_);
	testMultiplyBy();
	testAdd();
	testAddElements();
	testAddElements2();
	testSCF();
	testFillSequence();
	testFillSequence2();
	testFillSequence3();
	testFillSequence4();
}

} // namespace Lab
