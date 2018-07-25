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

#ifndef WINDOWFUNCTION_H_
#define WINDOWFUNCTION_H_

#include <cmath> /* ceil */
#include <sstream>
#include <string>
#include <vector>

#include "Exception.h"



namespace Lab {
namespace WindowFunction {

// n can be a non-integer.
template<typename FloatType> void fadeIn(double n, FloatType* data);
// The first zeroed element is not included in the data.
// n can be a non-integer.
template<typename FloatType> void fadeIn2(double n, FloatType* data);
// n can be a non-integer.
template<typename FloatType> void fadeOut(double n, FloatType* data);
// The last zeroed element is not included in the data.
// n can be a non-integer.
template<typename FloatType> void fadeOut2(double n, FloatType* data);
template<typename FloatType> void rectangular(unsigned int n, std::vector<FloatType>& w);
template<typename FloatType> void trapezoidal(unsigned int nFadeIn, unsigned int nTop, unsigned int nFadeOut, std::vector<FloatType>& w);
// The first and the last zeroed elements are not included in the array.
template<typename FloatType> void trapezoidal2(unsigned int nFadeIn, unsigned int nTop, unsigned int nFadeOut, std::vector<FloatType>& w);
template<typename FloatType> void triangular(unsigned int n, std::vector<FloatType>& w);
// The first and the last zeroed elements are not included in the array.
template<typename FloatType> void triangular2(unsigned int n, std::vector<FloatType>& w);
template<typename FloatType> void hanning(unsigned int n, std::vector<FloatType>& w);
// The first and the last zeroed elements are not included in the array.
template<typename FloatType> void hanning2(unsigned int n, std::vector<FloatType>& w);
template<typename FloatType> void hamming(unsigned int n, std::vector<FloatType>& w);
template<typename FloatType> void blackman(unsigned int n, std::vector<FloatType>& w);
// The first and the last zeroed elements are not included in the array.
template<typename FloatType> void blackman2(unsigned int n, std::vector<FloatType>& w);



template<typename FloatType>
void
fadeIn(double n, FloatType* data)
{
	const double c1 = 1.0 / n;
	for (unsigned int i = 0; i < n; ++i) {
		data[i] = i * c1;
	}
}

template<typename FloatType>
void
fadeIn2(double n, FloatType* data)
{
	const double c1 = 1.0 / (n + 1.0);
	for (unsigned int i = 0; i < n; ++i) {
		data[i] = (i + 1.0) * c1;
	}
}

template<typename FloatType>
void
fadeOut(double n, FloatType* data)
{
	const double c1 = 1.0 / n;
	for (unsigned int i = 0; i < n; ++i) {
		data[static_cast<unsigned int>(std::ceil(n)) - i - 1U] = i * c1;
	}
}

template<typename FloatType>
void
fadeOut2(double n, FloatType* data)
{
	const double c1 = 1.0 / (n + 1.0);
	for (unsigned int i = 0; i < n; ++i) {
		data[static_cast<unsigned int>(std::ceil(n)) - i - 1U] = (i + 1.0) * c1;
	}
}

template<typename FloatType>
void
rectangular(unsigned int n, std::vector<FloatType>& w)
{
	if (n < 2U) THROW_EXCEPTION(InvalidParameterException, "The window function size must be >= 2.");
	w.assign(n, 1.0);
}

template<typename FloatType>
void
trapezoidal(unsigned int nFadeIn, unsigned int nTop, unsigned int nFadeOut, std::vector<FloatType>& w)
{
	const unsigned int n = nFadeIn + nTop + nFadeOut;
	if (n < 2U) THROW_EXCEPTION(InvalidParameterException, "The window function size must be >= 2.");
	w.resize(n);
	fadeIn(nFadeIn, &w[0]);
	for (unsigned int i = nFadeIn; i < nFadeIn + nTop; ++i) {
		w[i] = 1.0;
	}
	fadeOut(nFadeOut, &w[nFadeIn + nTop]);
}

template<typename FloatType>
void
trapezoidal2(unsigned int nFadeIn, unsigned int nTop, unsigned int nFadeOut, std::vector<FloatType>& w)
{
	const unsigned int n = nFadeIn + nTop + nFadeOut;
	if (n < 2U) THROW_EXCEPTION(InvalidParameterException, "The window function size must be >= 2.");
	w.resize(n);
	fadeIn2(nFadeIn, &w[0]);
	for (unsigned int i = nFadeIn; i < nFadeIn + nTop; ++i) {
		w[i] = 1.0;
	}
	fadeOut2(nFadeOut, &w[nFadeIn + nTop]);
}

template<typename FloatType>
void
triangular(unsigned int n, std::vector<FloatType>& w)
{
	if (n < 2U) THROW_EXCEPTION(InvalidParameterException, "The window function size must be >= 2.");
	w.resize(n);
	const double halfWidth = 0.5 * (n - 1.0);
	if (n % 2 == 0) {
		fadeIn(halfWidth, &w[0]);
		fadeOut(halfWidth, &w[n / 2]);
	} else {
		fadeIn(halfWidth, &w[0]);
		const unsigned int n2 = (n - 1) / 2;
		w[n2] = 1.0;
		fadeOut(halfWidth, &w[n2 + 1]);
	}
}

template<typename FloatType>
void
triangular2(unsigned int n, std::vector<FloatType>& w)
{
	if (n < 2U) THROW_EXCEPTION(InvalidParameterException, "The window function size must be >= 2.");
	w.resize(n);
	const double halfWidth = 0.5 * (n - 1.0);
	if (n % 2 == 0) {
		fadeIn2(halfWidth, &w[0]);
		fadeOut2(halfWidth, &w[n / 2]);
	} else {
		fadeIn2(halfWidth, &w[0]);
		const unsigned int n2 = (n - 1) / 2;
		w[n2] = 1.0;
		fadeOut2(halfWidth, &w[n2 + 1]);
	}
}

template<typename FloatType>
void
hanning(unsigned int n, std::vector<FloatType>& w)
{
	if (n < 2U) THROW_EXCEPTION(InvalidParameterException, "The window function size must be >= 2.");
	w.resize(n);
	const double c1 = 1.0 / (n - 1.0);
	for (unsigned int i = 0; i < n; ++i) {
		const double c2 = i * c1;
		w[i] = 0.5 - 0.5 * std::cos(2.0 * PI * c2);
	}
}

template<typename FloatType>
void
hanning2(unsigned int n, std::vector<FloatType>& w)
{
	if (n < 2U) THROW_EXCEPTION(InvalidParameterException, "The window function size must be >= 2.");
	w.resize(n);
	const double c1 = 1.0 / (n + 1.0);
	for (unsigned int i = 0; i < n; ++i) {
		const double c2 = (i + 1.0) * c1;
		w[i] = 0.5 - 0.5 * std::cos(2.0 * PI * c2);
	}
}

template<typename FloatType>
void
hamming(unsigned int n, std::vector<FloatType>& w)
{
	if (n < 2U) THROW_EXCEPTION(InvalidParameterException, "The window function size must be >= 2.");
	w.resize(n);
	const double c1 = 1.0 / (n - 1.0);
	for (unsigned int i = 0; i < n; ++i) {
		const double c2 = i * c1;
		w[i] = 0.54 - 0.46 * std::cos(2.0 * PI * c2);
	}
}

template<typename FloatType>
void
blackman(unsigned int n, std::vector<FloatType>& w)
{
	if (n < 2U) THROW_EXCEPTION(InvalidParameterException, "The window function size must be >= 2.");
	w.resize(n);
	const double c1 = 1.0 / (n - 1.0);
	for (unsigned int i = 0; i < n; ++i) {
		const double c2 = i * c1;
		w[i] = 0.42 - 0.5 * std::cos(2.0 * PI * c2) + 0.08 * std::cos(4.0 * PI * c2);
	}
}

template<typename FloatType>
void
blackman2(unsigned int n, std::vector<FloatType>& w)
{
	if (n < 2U) THROW_EXCEPTION(InvalidParameterException, "The window function size must be >= 2.");
	w.resize(n);
	const double c1 = 1.0 / (n + 1.0);
	for (unsigned int i = 0; i < n; ++i) {
		const double c2 = (i + 1.0) * c1;
		w[i] = 0.42 - 0.5 * std::cos(2.0 * PI * c2) + 0.08 * std::cos(4.0 * PI * c2);
	}
}

template<typename FloatType>
void
get(const std::string& description, unsigned int size, std::vector<FloatType>& window)
{
	std::istringstream in{description};
	std::string name;
	in >> name;
	if (!in) THROW_EXCEPTION(MissingValueException, "Missing window function name.");

	if (name == "rectangular") {
		rectangular(size, window);
	} else if (name == "triangular") {
		triangular(size, window);
	} else if (name == "triangular2") {
		triangular2(size, window);
	} else if (name == "hanning") {
		hanning(size, window);
	} else if (name == "hanning2") {
		hanning2(size, window);
	} else if (name == "hamming") {
		hamming(size, window);
	} else if (name == "blackman") {
		blackman(size, window);
	} else if (name == "blackman2") {
		blackman2(size, window);
	} else if (name == "trapezoidal") {
		unsigned int nFadeIn, nTop, nFadeOut;
		in >> nFadeIn;
		if (!in) THROW_EXCEPTION(InvalidParameterException, "Missing fade in size for the trapezoidal window function.");
		in >> nTop;
		if (!in) THROW_EXCEPTION(InvalidParameterException, "Missing top size for the trapezoidal window function.");
		in >> nFadeOut;
		if (!in) THROW_EXCEPTION(InvalidParameterException, "Missing fade out size for the trapezoidal window function.");
		if (size != nFadeIn + nTop + nFadeOut) {
			THROW_EXCEPTION(InvalidParameterException, "Invalid total size for the trapezoidal window function.");
		}
		trapezoidal(nFadeIn, nTop, nFadeOut, window);
	} else if (name == "trapezoidal2") {
		unsigned int nFadeIn, nTop, nFadeOut;
		in >> nFadeIn;
		if (!in) THROW_EXCEPTION(InvalidParameterException, "Missing fade in size for the trapezoidal window function.");
		in >> nTop;
		if (!in) THROW_EXCEPTION(InvalidParameterException, "Missing top size for the trapezoidal window function.");
		in >> nFadeOut;
		if (!in) THROW_EXCEPTION(InvalidParameterException, "Missing fade out size for the trapezoidal window function.");
		if (size != nFadeIn + nTop + nFadeOut) {
			THROW_EXCEPTION(InvalidParameterException, "Invalid total size for the trapezoidal window function.");
		}
		trapezoidal2(nFadeIn, nTop, nFadeOut, window);
	} else {
		THROW_EXCEPTION(InvalidParameterException, "Invalid window function name: " << name << '.');
	}
}

} // namespace WindowFunction
} // namespace Lab

#endif /* WINDOWFUNCTION_H_ */
