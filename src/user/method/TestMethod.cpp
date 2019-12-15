/***************************************************************************
 *  Copyright 2014, 2017, 2018, 2019 Marcelo Y. Matuda                     *
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
#include <functional>
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

//#define TEST_DIRECT_FFTW_FILTER_SHOW_FIGURES 1
//#define TEST_FFTW_FILTER_SHOW_FIGURES 1
//#define TEST_HILBERT_TRANSFORM_SHOW_FIGURES 1
//#define TEST_INTERPOLATOR_SHOW_FIGURES 1
//#define TEST_KAISER_WINDOW_SHOW_FIGURES 1
//#define TEST_DEBUG 1

#define EXEC(A) LOG_INFO << #A; call(&TestMethod::A);

namespace {

constexpr double BESSEL_I_MAX_RELATIVE_ABS_ERROR = 1.0e-14;
constexpr double CENTRAL_DIFF_MAX_RELATIVE_ABS_ERROR = 1.0e-15;
constexpr double DECIMATOR_LP_FILTER_TRANSITION_WIDTH = 0.3; // fraction of the destination bandwidth
constexpr double DECIMATOR_MAX_RELATIVE_ABS_ERROR = 1.0e-14;
constexpr unsigned int DECIMATOR_OFFSET = 32;
constexpr double FFT_MAX_RELATIVE_ABS_ERROR = 1.0e-12;
constexpr double FILL_SEQUENCE_MAX_ABS_ERROR = 1.0e-15;
constexpr double FILTER_MAX_RELATIVE_ABS_ERROR = 1.0e-14;
constexpr double HILBERT_MAX_RELATIVE_ABS_ERROR = 1.0e-15;
constexpr double IFFT_MAX_RELATIVE_ABS_ERROR = 1.0e-10;
constexpr double INTERPOLATOR_LP_FILTER_TRANSITION_WIDTH = 0.45; // fraction of the original bandwidth
constexpr double INTERPOLATOR_MAX_RELATIVE_ABS_ERROR = 1.0e-14;
constexpr double KAISER_MAX_RELATIVE_ABS_ERROR = 1.0e-15;

} // namespace

namespace Lab {

TestMethod::TestMethod(Project& project)
		: project_(project)
		, errorCount_()
		, figureNumber_()
{
}

void
TestMethod::execute()
{
	errorCount_ = 0;
	figureNumber_ = 0;

	EXEC(testAdd)
	EXEC(testAddElements)
	EXEC(testAddElements2)
	EXEC(testBessel)
	EXEC(testCentralDiff)
	EXEC(testDecimator)
	EXEC(testDirectFFTWFilter)
	EXEC(testFFT)
	EXEC(testFFTWFilter)
	EXEC(testFFTWFilter2)
	EXEC(testFillSequence)
	EXEC(testFillSequence2)
	EXEC(testFillSequence3)
	EXEC(testFillSequence4)
	EXEC(testHilbertTransform)
	EXEC(testInterpolator)
	EXEC(testInterpolator4X)
	EXEC(testKaiserWindow)
	EXEC(testLinearInterpolator)
	EXEC(testMultiplyBy)
	EXEC(testSCF)
	EXEC(testStatistics)

	if (errorCount_ > 0) {
		LOG_ERROR << "Number of errors: " << errorCount_;
	} else {
		LOG_INFO << "No errors.";
	}
}

template<typename T>
void
TestMethod::call(T pmf)
{
	try {
		auto f = std::mem_fn(pmf);
		f(this);
		LOG_INFO << "  [OK]";
	} catch (std::exception& e) {
		++errorCount_;
		LOG_ERROR << e.what();
	}
}

void
TestMethod::testAdd()
{
	Matrix<double> m(3, 4);
	m(1, 0) = 1.0;
	m(1, 1) = 1.5;
	m(1, 2) = 2.0;
	m(1, 3) = 18.0;
	auto interval = m.dim2Interval(1);
	std::for_each(interval.first, interval.second, Util::Add<double>(-10.0));
	if (m(1, 0) != -9.0 ||
			m(1, 1) != -8.5 ||
			m(1, 2) != -8 ||
			m(1, 3) != 8.0) {
		THROW_EXCEPTION(TestException, "Wrong results.");
	}
}

void
TestMethod::testAddElements()
{
	std::vector<double> v1 = { 10.0, 11.0, 12.0, 13.0, 14.0 };
	std::vector<double> v2(4);
	Util::addElements(v1.begin(), v2.begin(), v2.end());
	if (v2[0] != 10.0 ||
			v2[1] != 11.0 ||
			v2[2] != 12.0 ||
			v2[3] != 13.0) {
		THROW_EXCEPTION(TestException, "Wrong results.");
	}
}

void
TestMethod::testAddElements2()
{
	std::vector<double> v1 = { 100.0, 101.0, 102.0, 103.0 };
	std::vector<double> v2 = {  20.0,  30.0,  40.0,  50.0 };
	std::vector<double> v3(3);
	Util::addElements(v1.begin(), v2.begin(), v3.begin(), v3.end());
	if (v3[0] != 120.0 ||
			v3[1] != 131.0 ||
			v3[2] != 142.0) {
		THROW_EXCEPTION(TestException, "Wrong results.");
	}
}

void
TestMethod::testBessel()
{
	std::vector<double> xBesselI0, besselI0;
	project_.loadHDF5("bessel_i0_x", "v", xBesselI0);
	project_.loadHDF5("bessel_i0"  , "v", besselI0);

	for (unsigned int i = 0; i < xBesselI0.size(); ++i) {
		const double x = xBesselI0[i];
		const double f = std::cyl_bessel_i(0, x);
		const double error = std::abs(f - besselI0[i]) / besselI0[i];
		if (error > BESSEL_I_MAX_RELATIVE_ABS_ERROR) {
			THROW_EXCEPTION(TestException, "Wrong value for x = " << x << ": " << f << " (error: " << error << ").");
		}
	}
}

void
TestMethod::testCentralDiff()
{
	std::vector<double> x, y, yRef;
	project_.loadHDF5("central_diff_x", "v", x);
	project_.loadHDF5("central_diff_y", "v", yRef);

	Util::centralDiff(x, 1.0 / 40e6, y);

	if (y.size() != yRef.size()) {
		THROW_EXCEPTION(TestException, "Wrong y size.");
	}

	double maxAbs = Util::maxAbsolute(yRef);
	for (unsigned int i = 0; i < y.size(); ++i) {
		const double error = std::abs(y[i] - yRef[i]) / maxAbs;
		if (error > CENTRAL_DIFF_MAX_RELATIVE_ABS_ERROR) {
			THROW_EXCEPTION(TestException, "Wrong value for index = " <<
					i << ": " << y[i] << " (error: " << error << ").");
		}
	}
}

void
TestMethod::testDecimator()
{
	const std::vector<unsigned int> downsampFactorList = { 2, 5, 10 };
	std::vector<double> x, y, y2, yRef;
	project_.loadHDF5("decimation_source", "v", x);

	for (unsigned int i = 0; i < downsampFactorList.size(); ++i) {
		Decimator<double> decimator;
		std::size_t yOffset;

		decimator.prepare(downsampFactorList[i], DECIMATOR_LP_FILTER_TRANSITION_WIDTH);
		decimator.decimate(DECIMATOR_OFFSET, x, yOffset, y);

		std::ostringstream fileName;
		fileName << "decimated_" << downsampFactorList[i] << 'x';
		project_.loadHDF5(fileName.str().c_str(), "v", yRef);

		if (y.size() != yRef.size()) {
			THROW_EXCEPTION(TestException, "Wrong size: " << y.size() <<
					" (downsampling factor: " << downsampFactorList[i] << ").");
		}
		const double maxAbs = Util::maxAbsolute(yRef);
		for (unsigned int j = 0; j < y.size(); ++j) {
			if (std::isnan(y[j])) {
				THROW_EXCEPTION(TestException, "NaN value for y[" <<
						j << "] (downsampling factor: " << downsampFactorList[i] << ").");
			}
			const double error = std::abs(y[j] - yRef[j]) / maxAbs;
			if (error > DECIMATOR_MAX_RELATIVE_ABS_ERROR) {
				THROW_EXCEPTION(TestException, "Wrong value for index = " <<
						j << ": " << y[j] << " (error: " << error <<
						", downsampling factor: " << downsampFactorList[i] << ").");
			}
		}

		std::ostringstream offsetFileName;
		offsetFileName << project_.directory() << "/decimated_" << downsampFactorList[i] << "x-offset.txt";
		std::ifstream in(offsetFileName.str().c_str(), std::ios_base::binary);
		if (!in) THROW_EXCEPTION(TestException, "Could not open the file: " << offsetFileName.str() << '.');
		unsigned int offsetRef;
		in >> offsetRef;
		if (!in) THROW_EXCEPTION(TestException, "Could not read the offset from file: " << offsetFileName.str() << '.');
		if (offsetRef != yOffset) {
			THROW_EXCEPTION(TestException, "Error: yOffset != offsetRef.");
		}

		Decimator<double> decimator2(decimator);
		std::size_t yOffset2;
		decimator2.decimate(DECIMATOR_OFFSET, x, yOffset2, y2);
		if (y2.size() != y.size()) {
			THROW_EXCEPTION(TestException, "Error: y2.size() != y.size().");
		}
		for (unsigned int j = 0; j < y2.size(); ++j) {
			if (y2[j] != y[j]) {
				THROW_EXCEPTION(TestException, "Error: y2 != y (y2[" << j << "] = " << y2[j] <<
						" y[" << j << "] = " << y[j] <<
						" downsampling factor = " << downsampFactorList[i] << ").");
			}
		}
		if (yOffset != yOffset2) {
			THROW_EXCEPTION(TestException, "Error: yOffset2 != yOffset.");
		}
	}
}

void
TestMethod::testDirectFFTWFilter()
{
	std::vector<double> x;
	project_.loadHDF5("filter_x", "x", x);
	std::vector<double> b;
	project_.loadHDF5("filter_b", "b", b);
	std::vector<double> yRef;
	project_.loadHDF5("filter_y", "y", yRef);

	std::vector<double> y;
	DirectFFTWFilter<double> filter;
	filter.setCoefficients(b);
	filter.filter(x, y);
#ifdef TEST_DIRECT_FFTW_FILTER_SHOW_FIGURES
	std::vector<double> t;
	Util::fillSequenceFromStartToEndWithSize(t, 0.0, static_cast<double>(y.size() - 1), y.size());
	project_.showFigure2D(figureNumber_++, "testDirectFFTWFilter: y", t, y);
#endif
	if (y.size() != x.size() + b.size() - 1) {
		THROW_EXCEPTION(TestException, "y.size() != x.size() + b.size() - 1 [y.size()=" << y.size() <<
			" x.size()=" << x.size() <<
			" b.size()=" << b.size() << "].");
	}
	double maxError = 0.0, maxY = 0.0;
	for (auto iter1 = yRef.cbegin(), iter2 = y.cbegin(); iter1 != yRef.cend(); ++iter1, ++iter2) {
		const double error = std::abs(*iter1 - *iter2);
		if (error > maxError) maxError = error;
		if (std::abs(*iter1) > maxY) maxY = std::abs(*iter1);
	}
	const double maxRelError = maxError / maxY;
#ifdef TEST_DEBUG
	LOG_DEBUG << "maxError = " << maxError << " maxY = " << maxY << " maxRelError = " << maxRelError;
#endif
	if (maxRelError > FILTER_MAX_RELATIVE_ABS_ERROR) {
		THROW_EXCEPTION(TestException, "The maximum relative error is > " << FILTER_MAX_RELATIVE_ABS_ERROR << '.');
	}
#if 1
	DirectFFTWFilter<double> filter2(filter);
#else
	DirectFFTWFilter<double> filter2;
	filter2 = filter;
#endif
	std::vector<double> y2;
	filter2.filter(x, y2);
	if (y2.size() != y.size()) {
		THROW_EXCEPTION(TestException, "y2.size() != y.size()");
	}
#ifdef TEST_DIRECT_FFTW_FILTER_SHOW_FIGURES
	project_.showFigure2D(figureNumber_++, "testDirectFFTWFilter: y2", t, y2);
#endif
	for (std::size_t i = 0; i < y2.size(); ++i) {
		if (y2[i] != y[i]) {
			THROW_EXCEPTION(TestException, "y2 != y" << y2[i] - y[i]);
		}
	}
}

void
TestMethod::testFFT()
{
	std::vector<double> x1, x2, x3;
	std::vector<std::complex<double>> y1, y2, y3;
	std::vector<double> yr1Ref, yi1Ref, yr2Ref, yi2Ref, yr3Ref, yi3Ref;

	project_.loadHDF5("fft_x_4000", "v", x1);
	project_.loadHDF5("fft_x_2048", "v", x2);
	project_.loadHDF5("fft_x_3571", "v", x3);
	project_.loadHDF5("fft_yr_4000", "v", yr1Ref);
	project_.loadHDF5("fft_yi_4000", "v", yi1Ref);
	project_.loadHDF5("fft_yr_2048", "v", yr2Ref);
	project_.loadHDF5("fft_yi_2048", "v", yi2Ref);
	project_.loadHDF5("fft_yr_3571", "v", yr3Ref);
	project_.loadHDF5("fft_yi_3571", "v", yi3Ref);
	y1.resize(x1.size());
	y2.resize(x2.size());
	y3.resize(x3.size());

	RealToComplexFFT<double> fft;

	fft.calculate(&x1[0], x1.size(), &y1[0]);
	double maxAbs = Util::maxAbsolute(y1);
#ifdef TEST_DEBUG
	LOG_DEBUG << "maxAbs: " << maxAbs;
#endif
	for (unsigned int i = 0; i < x1.size(); ++i) {
		const double error = std::max(
					std::abs(y1[i].real() - yr1Ref[i]),
					std::abs(y1[i].imag() - yi1Ref[i])) / maxAbs;
		if (error > FFT_MAX_RELATIVE_ABS_ERROR) {
			THROW_EXCEPTION(TestException, "Wrong value for index = " <<
					i << ": (" << y1[i].real() << ", " << y1[i].imag() << ") (error: " << error << ").");
		}
	}

	for (int t = 0; t < 3; ++t) {
		fft.calculate(&x2[0], x2.size(), &y2[0]);
		maxAbs = Util::maxAbsolute(y2);
#ifdef TEST_DEBUG
		LOG_DEBUG << "maxAbs: " << maxAbs;
#endif
		for (unsigned int i = 0; i < x2.size(); ++i) {
			const double error = std::max(
						std::abs(y2[i].real() - yr2Ref[i]),
						std::abs(y2[i].imag() - yi2Ref[i])) / maxAbs;
			if (error > FFT_MAX_RELATIVE_ABS_ERROR) {
				THROW_EXCEPTION(TestException, "Wrong value for index = " <<
						i << ": (" << y2[i].real() << ", " << y2[i].imag() << ") (error: " << error << ").");
			}
		}
	}

	RealToComplexFFT<double> fft2;
	fft2.calculate(&x3[0], x3.size(), &y3[0]);
	maxAbs = Util::maxAbsolute(y3);
#ifdef TEST_DEBUG
	LOG_DEBUG << "maxAbs: " << maxAbs;
#endif
	for (unsigned int i = 0; i < x3.size(); ++i) {
		const double error = std::max(
					std::abs(y3[i].real() - yr3Ref[i]),
					std::abs(y3[i].imag() - yi3Ref[i])) / maxAbs;
		if (error > FFT_MAX_RELATIVE_ABS_ERROR) {
			THROW_EXCEPTION(TestException, "Wrong value for index = " <<
					i << ": (" << y3[i].real() << ", " << y3[i].imag() << ") (error: " << error << ").");
		}
	}

	std::vector<std::complex<double>> yc1(yr1Ref.size());
	std::vector<std::complex<double>> yc2(yr2Ref.size());
	std::vector<std::complex<double>> yc3(yr3Ref.size());
	Value::copyRealImagToComplexSequence(yr1Ref.begin(), yr1Ref.end(), yi1Ref.begin(), yc1);
	Value::copyRealImagToComplexSequence(yr2Ref.begin(), yr2Ref.end(), yi2Ref.begin(), yc2);
	Value::copyRealImagToComplexSequence(yr3Ref.begin(), yr3Ref.end(), yi3Ref.begin(), yc3);
	std::vector<double> x1i(yr1Ref.size());
	std::vector<double> x2i(yr2Ref.size());
	std::vector<double> x3i(yr3Ref.size());

	ComplexToRealIFFT<double> ifft;

	for (int t = 0; t < 3; ++t) {
		ifft.calculate(&yc1[0], yc1.size(), &x1i[0]);
		maxAbs = Util::maxAbsolute(x1);
#ifdef TEST_DEBUG
		LOG_DEBUG << "maxAbs: " << maxAbs;
#endif
		for (unsigned int i = 0; i < x1i.size(); ++i) {
			const double error = std::abs(x1[i] - x1i[i]) / maxAbs;
			if (error > IFFT_MAX_RELATIVE_ABS_ERROR) {
				THROW_EXCEPTION(TestException, "Wrong value for index = " <<
						i << ": " << x1i[i] << " (error: " << error << ").");
			}
		}
	}

	ComplexToRealIFFT<double> ifft2;

	ifft2.calculate(&yc2[0], yc2.size(), &x2i[0]);
	maxAbs = Util::maxAbsolute(x2);
#ifdef TEST_DEBUG
	LOG_DEBUG << "maxAbs: "<< maxAbs;
#endif
	for (unsigned int i = 0; i < x2i.size(); ++i) {
		const double error = std::abs(x2[i] - x2i[i]) / maxAbs;
		if (error > IFFT_MAX_RELATIVE_ABS_ERROR) {
			THROW_EXCEPTION(TestException, "Wrong value for index = " <<
					i << ": " << x2i[i] << " (error: " << error << ").");
		}
	}

	ifft2.calculate(&yc3[0], yc3.size(), &x3i[0]);
	maxAbs = Util::maxAbsolute(x3);
#ifdef TEST_DEBUG
	LOG_DEBUG << "maxAbs: " << maxAbs;
#endif
	for (unsigned int i = 0; i < x3i.size(); ++i) {
		const double error = std::abs(x3[i] - x3i[i]) / maxAbs;
		if (error > IFFT_MAX_RELATIVE_ABS_ERROR) {
			THROW_EXCEPTION(TestException, "Wrong value for index = " <<
					i << ": " << x3i[i] << " (error: " << error << ").");
		}
	}
}

void
TestMethod::testFFTWFilter()
{
	std::vector<double> x;
	project_.loadHDF5("filter_x", "x", x);
	std::vector<double> b;
	project_.loadHDF5("filter_b", "b", b);
	std::vector<double> yRef;
	project_.loadHDF5("filter_y", "y", yRef);

	std::vector<double> y;
	FFTWFilter<double> filter;
	filter.setCoefficients(b);
	filter.filter(x, y);
#ifdef TEST_FFTW_FILTER_SHOW_FIGURES
	std::vector<double> t;
	Util::fillSequenceFromStartToEndWithSize(t, 0.0, static_cast<double>(y.size() - 1), y.size());
	project_.showFigure2D(figureNumber_++, "testFFTWFilter: y", t, y);
#endif
	if (y.size() != x.size() + b.size() - 1) {
		THROW_EXCEPTION(TestException, "y.size() != x.size() + b.size() - 1 [y.size()=" << y.size() <<
			" x.size()=" << x.size() <<
			" b.size()=" << b.size() << "].");
	}
	double maxError = 0.0, maxY = 0.0;
	for (auto iter1 = yRef.cbegin(), iter2 = y.cbegin(); iter1 != yRef.cend(); ++iter1, ++iter2) {
		const double error = std::abs(*iter1 - *iter2);
		if (error > maxError) maxError = error;
		if (std::abs(*iter1) > maxY) maxY = std::abs(*iter1);
	}
	const double maxRelError = maxError / maxY;
#ifdef TEST_DEBUG
	LOG_DEBUG << "maxError = " << maxError << " maxY = " << maxY << " maxRelError = " << maxRelError;
#endif
	if (maxRelError > FILTER_MAX_RELATIVE_ABS_ERROR) {
		THROW_EXCEPTION(TestException, "The maximum relative error is > " << FILTER_MAX_RELATIVE_ABS_ERROR << '.');
	}
#if 1
	FFTWFilter<double> filter2(filter);
#else
	FFTWFilter<double> filter2;
	filter2 = filter;
#endif
	std::vector<double> y2;
	filter2.filter(x, y2);
	if (y2.size() != y.size()) {
		THROW_EXCEPTION(TestException, "y2.size() != y.size()");
	}
#ifdef TEST_FFTW_FILTER_SHOW_FIGURES
	project_.showFigure2D(figureNumber_++, "testFFTWFilter: y2", t, y2);
#endif
	for (std::size_t i = 0; i < y2.size(); ++i) {
		if (y2[i] != y[i]) {
			THROW_EXCEPTION(TestException, "y2 != y" << y2[i] - y[i]);
		}
	}
}

void
TestMethod::testFFTWFilter2()
{
	std::vector<double> x;
	project_.loadHDF5("filter_x", "x", x);
	std::vector<double> b;
	project_.loadHDF5("filter_b", "b", b);
	std::vector<double> yRef;
	project_.loadHDF5("filter_y", "y", yRef);
	std::vector<double> y2Ref;
	project_.loadHDF5("filter_y2", "y", y2Ref);

	std::vector<std::complex<double>> filterFreqCoeff;
	std::vector<double> y;
	FFTWFilter2<double> filter;
	filter.setCoefficients(b, filterFreqCoeff);
	filter.filter(filterFreqCoeff, x, y);
	{
		if (y.size() != x.size() + b.size() - 1) {
			THROW_EXCEPTION(TestException, "y.size() != x.size() + b.size() - 1 [y.size()=" << y.size() <<
				" x.size()=" << x.size() <<
				" b.size()=" << b.size() << "].");
		}
		double maxError = 0.0, maxY = 0.0;
		for (auto iter1 = yRef.cbegin(), iter2 = y.cbegin(); iter1 != yRef.cend(); ++iter1, ++iter2) {
			const double error = std::abs(*iter1 - *iter2);
			if (error > maxError) maxError = error;
			if (std::abs(*iter1) > maxY) maxY = std::abs(*iter1);
		}
		const double maxRelError = maxError / maxY;
#ifdef TEST_DEBUG
		LOG_DEBUG << "maxError = " << maxError << " maxY = " << maxY << " maxRelError = " << maxRelError;
#endif
		if (maxRelError > FILTER_MAX_RELATIVE_ABS_ERROR) {
			THROW_EXCEPTION(TestException, "The maximum relative error is > " << FILTER_MAX_RELATIVE_ABS_ERROR << '.');
		}
	}

	FFTWFilter2<double> filter2;
	filter2.setCoefficients(yRef, filterFreqCoeff);
	FFTWFilter2<double> filter3;
	std::vector<std::complex<double>> dummyFilterFreqCoeff;
	filter3.setCoefficients(yRef, dummyFilterFreqCoeff);
	std::vector<double> y2;
	filter3.filter(filterFreqCoeff, x, y2);
	{
		if (y2.size() != x.size() + y.size() - 1) {
			THROW_EXCEPTION(TestException, "y2.size() != x.size() + y.size() - 1 [y2.size()=" << y2.size() <<
				" x.size()=" << x.size() <<
				" y.size()=" << y.size() << "].");
		}
		double maxError = 0.0, maxY = 0.0;
		for (auto iter1 = y2Ref.cbegin(), iter2 = y2.cbegin(); iter1 != y2Ref.cend(); ++iter1, ++iter2) {
			const double error = std::abs(*iter1 - *iter2);
			if (error > maxError) maxError = error;
			if (std::abs(*iter1) > maxY) maxY = std::abs(*iter1);
		}
		const double maxRelError = maxError / maxY;
#ifdef TEST_DEBUG
		LOG_DEBUG << "maxError = " << maxError << " maxY = " << maxY << " maxRelError = " << maxRelError;
#endif
		if (maxRelError > FILTER_MAX_RELATIVE_ABS_ERROR) {
			THROW_EXCEPTION(TestException, "The maximum relative error is > " << FILTER_MAX_RELATIVE_ABS_ERROR << '.');
		}
	}
}

void
TestMethod::testFillSequence()
{
	std::vector<double> v;
	Util::fillSequenceFromStartToEndWithMaximumStep(v, 1.0, 5.0, 0.99);
#ifdef TEST_DEBUG
	LOG_DEBUG << "v = " << v;
#endif
	if (v.size() != 6) {
		THROW_EXCEPTION(TestException, "Wrong size.");
	}
	if (std::abs(v[0] - 1.0) > FILL_SEQUENCE_MAX_ABS_ERROR ||
			std::abs(v[1] - 1.8) > FILL_SEQUENCE_MAX_ABS_ERROR ||
			std::abs(v[2] - 2.6) > FILL_SEQUENCE_MAX_ABS_ERROR ||
			std::abs(v[3] - 3.4) > FILL_SEQUENCE_MAX_ABS_ERROR ||
			std::abs(v[4] - 4.2) > FILL_SEQUENCE_MAX_ABS_ERROR ||
			std::abs(v[5] - 5.0) > FILL_SEQUENCE_MAX_ABS_ERROR) {
		THROW_EXCEPTION(TestException, "Wrong results.");
	}
}

void
TestMethod::testFillSequence2()
{
	std::vector<double> v;
	Util::fillSequenceFromStartToEndWithSize(v, 1.0, 5.0, 6);
#ifdef TEST_DEBUG
	LOG_DEBUG << "v = " << v;
#endif
	if (v.size() != 6) {
		THROW_EXCEPTION(TestException, "Wrong size.");
	}
	if (std::abs(v[0] - 1.0) > FILL_SEQUENCE_MAX_ABS_ERROR ||
			std::abs(v[1] - 1.8) > FILL_SEQUENCE_MAX_ABS_ERROR ||
			std::abs(v[2] - 2.6) > FILL_SEQUENCE_MAX_ABS_ERROR ||
			std::abs(v[3] - 3.4) > FILL_SEQUENCE_MAX_ABS_ERROR ||
			std::abs(v[4] - 4.2) > FILL_SEQUENCE_MAX_ABS_ERROR ||
			std::abs(v[5] - 5.0) > FILL_SEQUENCE_MAX_ABS_ERROR) {
		THROW_EXCEPTION(TestException, "Wrong results.");
	}
}

void
TestMethod::testFillSequence3()
{
	{
		std::vector<int> v;
		std::vector<int> vRef = { 2, 3, 4, 5, 6 };
		Util::fillSequence(v, 2, 6, 1);
		if (v != vRef) {
			THROW_EXCEPTION(TestException, "[1] v != vRef.");
		}
	}
	{

		std::vector<int> v;
		std::vector<int> vRef = { -2, -3, -4, -5 };
		Util::fillSequence(v, -2, -5, -1);
		if (v != vRef) {
			THROW_EXCEPTION(TestException, "[2] v != vRef.");
		}
	}
}

void
TestMethod::testFillSequence4()
{
	{
		std::vector<unsigned int> v;
		std::vector<unsigned int> vRef = { 2, 3, 4, 5, 6, 7 };
		Util::fillSequence(v, 2U, 7U, 1);
		if (v != vRef) {
			THROW_EXCEPTION(TestException, "[1] v != vRef.");
		}
	}
	{
		std::vector<unsigned int> v;
		std::vector<unsigned int> vRef = { 10, 9 , 8, 7, 6 };
		Util::fillSequence(v, 10U, 6U, -1);
		if (v != vRef) {
			THROW_EXCEPTION(TestException, "[2] v != vRef.");
		}
	}
}

void
TestMethod::testHilbertTransform()
{
	std::vector<double> x, yaRef, yrRef, yiRef;
	project_.loadHDF5("hilbert_x", "v", x);
	project_.loadHDF5("hilbert_ya", "v", yaRef);
	project_.loadHDF5("hilbert_yr", "v", yrRef);
	project_.loadHDF5("hilbert_yi", "v", yiRef);
	std::vector<double> ya = x;
	std::vector<std::complex<double>> yc(x.size());

	HilbertEnvelope<double> e;
	e.calculate(&ya[0], ya.size());

#ifdef TEST_HILBERT_TRANSFORM_SHOW_FIGURES
	std::vector<double> idx;
	Util::fillSequenceFromStartToEndWithSize(idx, 0.0, x.size() - 1.0, x.size());
	project_.showFigure2D(figureNumber_++, "testHilbertTransform: ya", idx, ya);

	project_.showFigure2D(figureNumber_++, "testHilbertTransform: yaRef", idx, yaRef);
#endif

	if (ya.size() != yaRef.size()) {
		THROW_EXCEPTION(TestException, "ya.size() != yaRef.size() (ya.size() = " << ya.size() << ").");
	}
	double maxAbs = Util::maxAbsolute(yaRef);
	for (unsigned int i = 0; i < ya.size(); ++i) {
		const double error = std::abs(ya[i] - yaRef[i]) / maxAbs;
		if (error > HILBERT_MAX_RELATIVE_ABS_ERROR) {
			THROW_EXCEPTION(TestException, "Wrong value for i = " << i << ": " << ya[i] << " (error: " << error << ").");
		}
	}

	e.getAnalyticSignal(&x[0], x.size(), &yc[0]);
	if (yc.size() != yrRef.size()) {
		THROW_EXCEPTION(TestException, "yc.size() != yrRef.size() (yc.size() = " << yc.size() << ").");
	}
	maxAbs = Util::maxAbsolute(yc);
	for (unsigned int i = 0; i < yc.size(); ++i) {
		const double error = std::max(
					std::abs(yc[i].real() - yrRef[i]),
					std::abs(yc[i].imag() - yiRef[i])) / maxAbs;
		if (error > HILBERT_MAX_RELATIVE_ABS_ERROR) {
			THROW_EXCEPTION(TestException, "Wrong value for i = " << i << ": " << yc[i] << " (error: " << error << ").");
		}
	}

	HilbertEnvelope<double> e2(e);
	std::vector<double> ya2 = x;
	e2.calculate(&ya2[0], ya2.size());
	for (std::size_t i = 0; i < ya2.size(); ++i) {
		if (ya2[i] != ya[i]) {
			THROW_EXCEPTION(TestException, "ya2 != ya");
		}
	}

	std::vector<std::complex<double>> yc2(x.size());
	e2.getAnalyticSignal(&x[0], x.size(), &yc2[0]);
	for (std::size_t i = 0; i < yc2.size(); ++i) {
		if (yc2[i] != yc[i]) {
			THROW_EXCEPTION(TestException, "yc2 != yc");
		}
	}
}

void
TestMethod::testInterpolator()
{
	const std::vector<unsigned int> upsampFactorList = { 2, 8, 64 };
	std::vector<double> x, y, y2, yRef;
	project_.loadHDF5("interp_source", "v", x);

	for (unsigned int i = 0; i < upsampFactorList.size(); ++i) {
		Interpolator<double> interp;

		interp.prepare(upsampFactorList[i], INTERPOLATOR_LP_FILTER_TRANSITION_WIDTH);
		y.resize(x.size() * upsampFactorList[i]);
		interp.interpolate(&x[0], x.size(), &y[0]);

		std::ostringstream fileName;
		fileName << "interp_" << upsampFactorList[i] << 'x';
		project_.loadHDF5(fileName.str().c_str(), "v", yRef);
#ifdef TEST_INTERPOLATOR_SHOW_FIGURES
		std::vector<double> t;
		Util::fillSequenceFromStartToEndWithSize(t, 0.0, static_cast<double>(y.size() - 1), y.size());
		project_.showFigure2D(figureNumber_++, "testInterpolator: y", t, y);
#endif
		if (y.size() != yRef.size()) {
			THROW_EXCEPTION(TestException, "Wrong size: " << y.size() <<
					" (upsampling factor: " << upsampFactorList[i] << ").");
		}
		const double maxAbs = Util::maxAbsolute(yRef);
		for (unsigned int j = 0; j < y.size(); ++j) {
			const double error = std::abs(y[j] - yRef[j]) / maxAbs;
			if (error > INTERPOLATOR_MAX_RELATIVE_ABS_ERROR) {
				THROW_EXCEPTION(TestException, "Wrong value for index = " <<
						j << ": " << y[j] << " (error: " << error <<
						", upsampling factor: " << upsampFactorList[i] << ").");
			}
		}

		Interpolator<double> interp2(interp);
		y2.resize(x.size() * upsampFactorList[i]);
		interp2.interpolate(&x[0], x.size(), &y2[0]);
#ifdef TEST_INTERPOLATOR_SHOW_FIGURES
		project_.showFigure2D(figureNumber_++, "testInterpolator: y2", t, y2);
#endif
		for (unsigned int j = 0; j < y2.size(); ++j) {
			if (y2[j] != y[j]) {
				THROW_EXCEPTION(TestException, "y2 != y");
			}
		}
	}
}

void
TestMethod::testInterpolator4X()
{
	std::vector<double> x;
	project_.loadHDF5("interp_x", "x", x);

	std::vector<double> y(x.size() * 4);

	Interpolator4X<double> interpolator;
	interpolator.interpolate(&x[0], x.size(), &y[0]);
	project_.saveHDF5(y, "interp_y", "y");
}

void
TestMethod::testKaiserWindow()
{
	std::vector<double> tol_dB, betaRef;
	project_.loadHDF5("kaiser_tol_db", "v", tol_dB);
	project_.loadHDF5("kaiser_beta"  , "v", betaRef);
#ifdef TEST_KAISER_WINDOW_SHOW_FIGURES
	project_.showFigure2D(figureNumber_++, "Kaiser beta", tol_dB, betaRef);
#endif
	for (unsigned int i = 0; i < tol_dB.size(); ++i) {
		const double beta = KaiserWindow::getBeta(tol_dB[i]);
		const double error = std::abs(beta - betaRef[i]) / betaRef[i];
		if (error > KAISER_MAX_RELATIVE_ABS_ERROR) {
			THROW_EXCEPTION(TestException, "Wrong beta value for tol_dB = "
							<< tol_dB[i] << ": " << beta << " (error: " << error << ").");
		}
	}

	const std::vector<double> transWidthList = { 0.05, 0.1, 0.5 };
	std::vector<double> sizeRef;
	for (unsigned int i = 0; i < transWidthList.size(); ++i) {
		const double transWidth = transWidthList[i];
		std::ostringstream fileName;
		fileName << "kaiser_size-trans_width_" << transWidth;
		project_.loadHDF5(fileName.str().c_str(), "v", sizeRef);

		for (unsigned int j = 0; j < sizeRef.size(); ++j) {
			const unsigned int size = KaiserWindow::getSize(tol_dB[i], transWidth);
			const double error = std::abs(size - sizeRef[i]) / sizeRef[i];
			if (error > 0.0) {
				THROW_EXCEPTION(TestException, "Wrong size value for tol_dB = "
								<< tol_dB[i] << ": " << size << " (error: " << error << ").");
			}
		}
	}
}

void
TestMethod::testLinearInterpolator()
{
	const unsigned int upsamplingFactor = 4;
	std::vector<double> x = { 1.0, -1.0, 0.0, 0.5 };
	std::vector<double> y(x.size() * upsamplingFactor);
	LinearInterpolator<double, upsamplingFactor>::interpolate(&x[0], 4, &y[0]);

	if (
			y[0]  !=  1.0 || y[1]  !=  0.5   || y[2]  !=  0.0  || y[3]  != -0.5   ||
			y[4]  != -1.0 || y[5]  != -0.75  || y[6]  != -0.5  || y[7]  != -0.25  ||
			y[8]  !=  0.0 || y[9]  !=  0.125 || y[10] !=  0.25 || y[11] !=  0.375 ||
			y[12] !=  0.5 || y[13] !=  0.375 || y[14] !=  0.25 || y[15] !=  0.125) {
		THROW_EXCEPTION(TestException, "Wrong values: " << y << '.');
	}
}

void
TestMethod::testMultiplyBy()
{
	Matrix<double> m(3, 4);
	m(1, 0) = 1.0;
	m(1, 1) = 1.5;
	m(1, 2) = 2.0;
	m(1, 3) = 8.0;
	auto interval = m.dim2Interval(1);
	std::for_each(interval.first, interval.second, Util::MultiplyBy<double>(-0.5));
	if (m(1, 0) != -0.5 ||
			m(1, 1) != -0.75 ||
			m(1, 2) != -1.0 ||
			m(1, 3) != -4.0) {
		THROW_EXCEPTION(TestException, "Wrong results.");
	}
}

void
TestMethod::testSCF()
{
	std::vector<double> v1 = { 1.0, 0.0, -1.0, 1.0, -1.0, 1.0 };

	auto scf = std::make_unique<SignCoherenceFactor<double>>(2.0);

	const double scfValue = scf->calculate(&v1[0], v1.size());
	if (scfValue != 3.270805724762158e-3) {
		THROW_EXCEPTION(TestException, "Wrong result.");
	}
}

void
TestMethod::testStatistics()
{
	const std::vector<double> a = { -1.0, 2.0, 4.0, 7.0 };

	const double stdDev = Statistics::standardDeviation(a.data(), a.size());
	if (std::abs(stdDev - 2.9154759474) > 1e-10) {
		THROW_EXCEPTION(TestException, "Wrong stdDev: " << stdDev << '.');
	}

	const double mean = Statistics::arithmeticMean(a.data(), a.size());
	if (mean != 3.0) {
		THROW_EXCEPTION(TestException, "Wrong mean: " << mean << '.');
	}
}

} // namespace Lab
