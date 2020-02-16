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

#ifndef VALUE_H_
#define VALUE_H_

#include <cmath> /* abs, sqrt */
#include <complex>
#include <vector>

#include "XYZValue.h"
#include "XYZValueFactor.h"
#include "XZValue.h"
#include "XZValueFactor.h"



namespace Lab {
namespace Value {

template<typename TFloat> void copy(const XYZValue<TFloat>& orig, TFloat& dest);
template<typename TFloat> void copy(const TFloat& orig, XYZValue<TFloat>& dest);
template<typename TFloat> void copy(const XZValue<TFloat>& orig, TFloat& dest);
template<typename TFloat> void copy(const TFloat& orig, XZValue<TFloat>& dest);

template<typename TFloat> void copy(const XYZValueFactor<TFloat>& orig, TFloat& dest);
template<typename TFloat> void copy(const TFloat& orig, XYZValueFactor<TFloat>& dest);
template<typename TFloat> void copy(const XZValueFactor<TFloat>& orig, TFloat& dest);
template<typename TFloat> void copy(const TFloat& orig, XZValueFactor<TFloat>& dest);

// array of size two --> complex
template<typename TFloat> void copy(const TFloat (&orig)[2], std::complex<TFloat>& dest);
// complex --> array of size two
template<typename TFloat> void copy(const std::complex<TFloat>& orig, TFloat (&dest)[2]);

template<typename TFloat1, typename TFloat2> void copy(const TFloat1& orig, TFloat2& dest);
template<typename TFloat1, typename TFloat2> void copy(const std::complex<TFloat1>& orig, TFloat2 &dest);

template<typename InputIterator, typename OutputIterator>
	void copySequence(InputIterator input, InputIterator inputEnd, OutputIterator output);
template<typename InputIterator, typename OutputIterator>
	void copySequenceWithPadding(InputIterator input, InputIterator inputEnd, OutputIterator output, unsigned int padding);
template<typename InputIterator, typename OutputIterator, typename T>
	void transformSequence(InputIterator orig, InputIterator origEnd, OutputIterator dest, T transform);
template<typename Iterator, typename TFloat>
	void copyRealImagToComplexSequence(Iterator re, Iterator reEnd, Iterator im, std::vector<std::complex<TFloat>>& cpx);
template<typename InputIterator, typename OutputIterator>
	void copyXYZValueSequence(InputIterator orig, InputIterator origEnd, OutputIterator dest);
template<typename InputIterator, typename OutputIterator>
	void copyXYZSequence(InputIterator orig, InputIterator origEnd, OutputIterator dest);
template<typename InputIterator, typename OutputIterator>
	void copyXZSequence(InputIterator orig, InputIterator origEnd, OutputIterator dest);
template<typename InputIterator, typename OutputIterator>
	void copyXZValueSequence(InputIterator orig, InputIterator origEnd, OutputIterator dest);



template<typename TFloat>
class ComplexToScaledAbsoluteOp {
public:
	ComplexToScaledAbsoluteOp(TFloat scale) : scale_(scale) { }

	template<typename DestType>
	void operator()(const std::complex<TFloat>& orig, DestType& dest) {
		copy(std::sqrt(orig.real() * orig.real() + orig.imag() * orig.imag()) * scale_, dest);
	}

	template<typename DestType>
	void operator()(const TFloat (&orig)[2], DestType& dest) {
		copy(std::sqrt(orig[0] * orig[0] + orig[1] * orig[1]) * scale_, dest);
	}
private:
	TFloat scale_;
};

template<typename T, typename U>
struct ScalarToValueFieldOp {
	void operator()(T& orig, U& dest) {
		dest.value = orig;
	}
};

template<typename TFloat>
class ScaleComplexOp {
public:
	ScaleComplexOp(TFloat scale) : scale_(scale) { }

	void operator()(const std::complex<TFloat>& orig, std::complex<TFloat>& dest) {
		dest = orig * scale_;
	}
	void operator()(const TFloat (&orig)[2], std::complex<TFloat>& dest) {
		dest = std::complex<TFloat>(orig[0], orig[1]) * scale_;
	}
	void operator()(const std::complex<TFloat>& orig, TFloat (&dest)[2]) {
		dest[0] = orig.real() * scale_;
		dest[0] = orig.imag() * scale_;
	}
	void operator()(const TFloat (&orig)[2], TFloat (&dest)[2]) {
		dest[0] = orig[0] * scale_;
		dest[0] = orig[1] * scale_;
	}
private:
	TFloat scale_;
};

template<typename T>
class ScaleOp {
public:
	ScaleOp(T scale) : scale_(scale) { }

	void operator()(const T& orig, T& dest) {
		dest = orig * scale_;
	}
private:
	T scale_;
};



template<typename TFloat>
void
copy(const XYZValue<TFloat>& orig, TFloat& dest)
{
	dest = orig.value;
}

template<typename TFloat>
void
copy(const TFloat& orig, XYZValue<TFloat>& dest)
{
	dest.value = orig;
}

template<typename TFloat>
void
copy(const XZValue<TFloat>& orig, TFloat& dest)
{
	dest = orig.value;
}

template<typename TFloat>
void
copy(const TFloat& orig, XZValue<TFloat>& dest)
{
	dest.value = orig;
}

template<typename TFloat>
void
copy(const XYZValueFactor<TFloat>& orig, TFloat& dest)
{
	dest = orig.value;
}

template<typename TFloat>
void
copy(const TFloat& orig, XYZValueFactor<TFloat>& dest)
{
	dest.value = orig;
}

template<typename TFloat>
void
copy(const XZValueFactor<TFloat>& orig, TFloat& dest)
{
	dest = orig.value;
}

template<typename TFloat>
void
copy(const TFloat& orig, XZValueFactor<TFloat>& dest)
{
	dest.value = orig;
}

template<typename TFloat>
void
copy(const TFloat (&orig)[2], std::complex<TFloat>& dest)
{
	dest = std::complex<TFloat>(orig[0], orig[1]);
}

template<typename TFloat>
void
copy(const std::complex<TFloat>& orig, TFloat (&dest)[2])
{
	dest[0] = orig.real();
	dest[1] = orig.imag();
}

template<typename TFloat1, typename TFloat2>
void
copy(const TFloat1& orig, TFloat2& dest)
{
	dest = orig;
}

template<typename TFloat1, typename TFloat2>
void
copy(const std::complex<TFloat1>& orig, TFloat2 &dest)
{
	dest = std::abs(orig);
}

template<typename InputIterator, typename OutputIterator>
void
copySequence(InputIterator orig, InputIterator origEnd, OutputIterator dest)
{
	while (orig != origEnd) {
		copy(*orig++, *dest++);
	}
}

template<typename InputIterator, typename OutputIterator>
void
copySequenceWithPadding(InputIterator orig, InputIterator origEnd, OutputIterator dest, unsigned int padding)
{
	while (orig != origEnd) {
		copy(*orig++, *dest++);
	}

	// Padding.
	OutputIterator destEnd = dest + padding;
	while (dest != destEnd) {
		*dest++ = 0;
	}
}

template<typename InputIterator, typename OutputIterator, typename T>
void
transformSequence(InputIterator orig, InputIterator origEnd, OutputIterator dest, T transform)
{
	while (orig != origEnd) {
		transform(*orig++, *dest++);
	}
}

template<typename Iterator, typename TFloat>
void
copyRealImagToComplexSequence(Iterator re, Iterator reEnd, Iterator im, std::vector<std::complex<TFloat>>& cpx)
{
	typename std::vector<std::complex<TFloat>>::iterator cpxIter = cpx.begin();
	while (re != reEnd) {
		*cpxIter++ = std::complex<TFloat>(*re++, *im++);
	}
}

template<typename InputIterator, typename OutputIterator>
void
copyXYZValueSequence(InputIterator orig, InputIterator origEnd, OutputIterator dest)
{
	for ( ; orig != origEnd; ++orig, ++dest) {
		dest->x = orig->x;
		dest->y = orig->y;
		dest->z = orig->z;
		copy(orig->value, dest->value);
	}
}

template<typename InputIterator, typename OutputIterator>
void
copyXYZSequence(InputIterator orig, InputIterator origEnd, OutputIterator dest)
{
	for ( ; orig != origEnd; ++orig, ++dest) {
		dest->x = orig->x;
		dest->y = orig->y;
		dest->z = orig->z;
	}
}

template<typename InputIterator, typename OutputIterator>
void
copyXZSequence(InputIterator orig, InputIterator origEnd, OutputIterator dest)
{
	for ( ; orig != origEnd; ++orig, ++dest) {
		dest->x = orig->x;
		dest->z = orig->z;
	}
}

template<typename InputIterator, typename OutputIterator>
void
copyXZValueSequence(InputIterator orig, InputIterator origEnd, OutputIterator dest)
{
	for ( ; orig != origEnd; ++orig, ++dest) {
		dest->x = orig->x;
		dest->z = orig->z;
		copy(orig->value, dest->value);
	}
}

} // namespace Value
} // namespace Lab

#endif // VALUE_H_
