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

#ifndef UTIL_H_
#define UTIL_H_

#include <algorithm> /* for_each, max_element */
#include <cstddef> /* std::size_t */
#include <vector>
#ifdef _WIN32
# define WIN32_LEAN_AND_MEAN
# include <windows.h> /* Sleep */
#else
# include <ctime> /* nanosleep */
#endif

#include <cmath>
#include <complex>
#include <cstring>
#include <limits>
#include <numeric> /* accumulate */

#include "Exception.h"
#include "Matrix2.h"
#include "XZ.h"
#include "XZValue.h"
#include "XZValueFactor.h"

#define PI 3.14159265358979323846



namespace Lab {
namespace Util {

void fillSequence(std::vector<unsigned int>& v, unsigned int startValue, unsigned int endValue, int step = 1);
void fillSequence(std::vector<int>& v, int startValue, int endValue, int step = 1);

// This function may not give good results with integer values.
// The sequence will contain _startValue_ and _endValue_. The step may be different.
template<typename T> void fillSequenceFromStartToEndWithMaximumStep(std::vector<T>& v, T startValue, T endValue, T step = 1);

// This function may not give good results with integer values.
// The sequence will contain _startValue_ and _endValue_, with exactly _size_ items.
template<typename T> void fillSequenceFromStartToEndWithSize(std::vector<T>& v, T startValue, T endValue, unsigned int size);

// The sequence will contain _startValue_, with exactly _size_ items.
template<typename T> void fillSequenceFromStartWithStepAndSize(std::vector<T>& v, T startValue, T step, unsigned int size);

void sleepMs(unsigned long milliseconds);
template<typename T> T meterToMillimeter(T value);
template<typename T> T millimeterToMeter(T value);
template<typename T> void clipAngleInDegrees180(T& angle);
template<typename T> void clip(T& value, T minValue, T maxValue);
template<typename T> T max(const std::vector<T>& list);
template<typename T> T maxAbsolute(const std::vector<T>& list);
template<typename T> T maxAbsolute(const std::vector<std::complex<T> >& list);
template<typename T> T maxAbsolute(const Matrix2<T>& data);
template<typename T, typename U> U maxValueField(const Matrix2<T>& data);
template<typename T, typename U> U maxAbsoluteValueField(const Matrix2<T>& data);
template<typename T> void minMax(const std::vector<T>& list, T& min, T& max);
template<typename T, typename U> void minMaxValueField(const Matrix2<T>& data, U& min, U& max);
template<typename T> void multiply(std::vector<T>& list, T coefficient);
template<typename T, typename U> void multiply(const std::vector<T>& factorList, std::vector<U>& data);
template<typename T> void centralDiff(const std::vector<T>& inputList, T period, std::vector<T>& outputList);
//template<typename FloatType, typename IntType, typename InputIterator, typename OutputIterator>
//	void copyFloatToInt(InputIterator input, InputIterator inputEnd, OutputIterator output);
template<typename T> void deleteObjects(T& container);
template<typename T> T decibelsToLinear(T decibels);
template<typename T> T linearToDecibels(T linear);
template<typename T> T degreeToRadian(T d);
template<typename T> T radianToDegree(T r);

// Adds the values pointed by first1 to the elements pointed by first2.
template<typename InputIterator, typename InputOutputIterator>
void addElements(InputIterator first1, InputOutputIterator first2, InputOutputIterator last2);

// Adds the values pointed by first1 to the values pointed by first2
// and sends the results to the elements pointed by first3.
template<typename InputIterator1, typename InputIterator2, typename OutputIterator>
void addElements(InputIterator1 first1, InputIterator2 first2, OutputIterator first3, OutputIterator last3);

template<typename T, typename U> void copyXZToSimpleMatrices(const Matrix2<T>& m, Matrix2<U>& mX, Matrix2<U>& mZ);
template<typename T, typename U> void copyXZFromSimpleMatrices(const Matrix2<T>& mX, const Matrix2<T>& mZ, Matrix2<U>& m);

template<typename T, typename U> void copyValueToSimpleMatrix(const Matrix2<T>& m, Matrix2<U>& mV);
template<typename T, typename U> void copyValueFromSimpleMatrix(const Matrix2<T>& mV, Matrix2<U>& m);

template<typename T, typename U> void copyFactorToSimpleMatrix(const Matrix2<T>& m, Matrix2<U>& mF);

template<typename T, typename U> void copyXZToSimpleVectors(const std::vector<T>& v, std::vector<U>& vX, std::vector<U>& vZ);
template<typename T, typename U> void copyXZFromSimpleVectors(const std::vector<T>& vX, const std::vector<T>& vZ, std::vector<U>& v);

template<typename T, typename U> void copyXZValue(const Matrix2<T>& orig, Matrix2<U>& dest);
template<typename T, typename U> void copyXZFactor(const Matrix2<T>& orig, Matrix2<U>& dest);
template<typename T, typename U> void copyXZ(const Matrix2<T>& orig, Matrix2<U>& dest);
template<typename T, typename U> void copyXZ(const std::vector<T>& orig, std::vector<U>& dest);

template<typename T> T minValue();
template<typename T> T maxValue();

template<typename T> T nextPowerOf2(T v);

template<typename T> T multiplyElements(const std::vector<T>& v);

template<typename InputIterator, typename OutputIterator, typename T> void copyUsingOperator(InputIterator orig, InputIterator origEnd, OutputIterator dest, T copyOp);

// arc step = arcStepFactor * lambda
template<typename FloatType>
	void generateArc2D(FloatType centerX, FloatType centerZ,
				FloatType radius, FloatType angle1, FloatType angle2,
				FloatType arcStepFactor, FloatType lambda,
				std::vector<XZ<FloatType> >& pointList);

template<typename FloatType>
	void calculateDelays2D(FloatType propagationSpeed, XZ<FloatType> focus,
				const std::vector<XZ<FloatType> >& arrayElemPosList, unsigned int numElements,
				std::vector<FloatType>& delays);

template<typename Iterator> void resetValueFactor(Iterator iter, Iterator iterEnd);

template<typename T> void removeDC(T* data, std::size_t size);
template<typename T> void removeDC(T* data, std::size_t size, std::size_t beginOffset);

//#############################################################################

template<typename T>
struct MultiplyBy {
	explicit MultiplyBy(T factor) : factor(factor) {}
	void operator()(T& value) { value *= factor; }
	T factor;
};

template<typename T>
struct Add {
	explicit Add(T term) : term(term) {}
	void operator()(T& value) { value += term; }
	T term;
};

template<typename T>
struct ValueFieldAbsOp {
	explicit ValueFieldAbsOp() {}
	void operator()(T& o) { o.value = std::abs(o.value); }
};



template<typename T>
void
deleteObjects(T& container)
{
	for (typename T::size_type i = 0, size = container.size(); i < size; ++i) {
		delete container[i];
	}
}

inline
void
fillSequence(std::vector<unsigned int>& v, unsigned int startValue, unsigned int endValue, int step)
{
	if (step == 0) {
		THROW_EXCEPTION(InvalidParameterException, "The step must not be zero.");
	}
	if (startValue < endValue && step < 0) {
		THROW_EXCEPTION(InvalidParameterException, "The step must be positive.");
	}
	if (startValue > endValue && step > 0) {
		THROW_EXCEPTION(InvalidParameterException, "The step must be negative.");
	}

	const std::size_t n = (step > 0) ? (endValue - startValue) / step + 1 : (startValue - endValue) / -step + 1;
	v.resize(n);
	v[0] = startValue;
	if (step > 0) {
		for (std::size_t i = 1; i < n; ++i) {
			v[i] = v[i - 1] + step;
		}
	} else {
		for (std::size_t i = 1; i < n; ++i) {
			v[i] = v[i - 1] - static_cast<unsigned int>(-step);
		}
	}
}

inline
void
fillSequence(std::vector<int>& v, int startValue, int endValue, int step)
{
	if (step == 0) {
		THROW_EXCEPTION(InvalidParameterException, "The step must not be zero.");
	}
	if (startValue < endValue && step < 0) {
		THROW_EXCEPTION(InvalidParameterException, "The step must be positive.");
	}
	if (startValue > endValue && step > 0) {
		THROW_EXCEPTION(InvalidParameterException, "The step must be negative.");
	}

	const std::size_t n = (endValue - startValue) / step + 1;
	v.resize(n);
	v[0] = startValue;
	for (std::size_t i = 1; i < n; ++i) {
		v[i] = v[i - 1] + step;
	}
}

template<typename T>
void
fillSequenceFromStartToEndWithMaximumStep(std::vector<T>& v, T startValue, T endValue, T step)
{
	if (step == 0) {
		THROW_EXCEPTION(InvalidParameterException, "The step must not be zero.");
	}
	if (startValue < endValue && step < 0) {
		THROW_EXCEPTION(InvalidParameterException, "The step must be positive.");
	}
	if (startValue > endValue && step > 0) {
		THROW_EXCEPTION(InvalidParameterException, "The step must be negative.");
	}

	const double numSteps = std::ceil((static_cast<double>(endValue) - startValue) / step);
	if (numSteps > static_cast<double>(std::numeric_limits<int>::max() - 2)) {
		THROW_EXCEPTION(InvalidParameterException, "The step is too small.");
	}

	const int n = 1 + static_cast<int>(numSteps);
	const double newStep = (endValue - startValue) / (n - 1);
	v.resize(n);
	for (int i = 0; i < n; ++i) {
		v[i] = startValue + i * newStep;
	}
}

template<typename T>
void
fillSequenceFromStartToEndWithSize(std::vector<T>& v, T startValue, T endValue, unsigned int size)
{
	if (size < 1) {
		THROW_EXCEPTION(InvalidParameterException, "The size must be >= 1.");
	}

	if (size == 1) {
		v.resize(1);
		v[0] = startValue;
	} else {
		const double step = static_cast<double>(endValue - startValue) / (size - 1U);
		v.resize(size);
		v[0] = startValue;
		for (unsigned int i = 1; i < size; ++i) {
			v[i] = startValue + i * step;
		}
	}
}

template<typename T>
void
fillSequenceFromStartWithStepAndSize(std::vector<T>& v, T startValue, T step, unsigned int size)
{
	if (size < 1) {
		THROW_EXCEPTION(InvalidParameterException, "The size must be >= 1.");
	}

	if (size == 1) {
		v.resize(1);
		v[0] = startValue;
	} else {
		v.resize(size);
		v[0] = startValue;
		for (unsigned int i = 1; i < size; ++i) {
			v[i] = startValue + i * step;
		}
	}
}

inline
void
sleepMs(unsigned long milliseconds)
{
#ifdef _WIN32
	Sleep(milliseconds);
#else
	struct timespec tspec;
	if (milliseconds >= 1000) {
		tspec.tv_sec = milliseconds / 1000;
		tspec.tv_nsec = (milliseconds % 1000) * 1000000;
	} else {
		tspec.tv_sec = 0;
		tspec.tv_nsec = milliseconds * 1000000;
	}
	nanosleep(&tspec, 0);
#endif
}

template<typename T>
T
meterToMillimeter(T value)
{
	return value * T(1.0e3);
}

template<typename T>
T
millimeterToMeter(T value)
{
	return value * T(1.0e-3);
}

template<typename T>
void
clipAngleInDegrees180(T& angle)
{
	if (angle > 180.0) {
		do {
			angle -= 360.0;
		} while (angle > 180.0);
	} else if (angle < -180.0) {
		do {
			angle += 360.0;
		} while (angle < -180.0);
	}
}

template<typename T>
void
clip(T& value, T minValue, T maxValue)
{
	if (value > maxValue) {
		value = maxValue;
	} else if (value < minValue) {
		value = minValue;
	}
}

template<typename T>
T
max(const std::vector<T>& list)
{
	T max = minValue<T>();

	for (typename std::vector<T>::const_iterator iter = list.begin(); iter != list.end(); ++iter) {
		if (max < *iter) {
			max = *iter;
		}
	}

	return max;
}

template<typename T>
T
maxAbsolute(const std::vector<T>& list)
{
	T max = 0;
	for (typename std::vector<T>::const_iterator iter = list.begin(); iter != list.end(); ++iter) {
		const T a = std::abs(*iter);
		if (max < a) max = a;
	}
	return max;
}

template<typename T>
T
maxAbsolute(const std::vector<std::complex<T> >& list)
{
	T max = 0;
	for (typename std::vector<std::complex<T> >::const_iterator iter = list.begin(); iter != list.end(); ++iter) {
		const T a = std::abs(*iter);
		if (max < a) max = a;
	}
	return max;
}

template<typename T>
T
maxAbsolute(const Matrix2<T>& data)
{
	T max = minValue<T>();
	for (typename Matrix2<T>::ConstIterator iter = data.begin(); iter != data.end(); ++iter) {
		const T v = std::abs(*iter);
		if (v > max) max = v;
	}
	return max;
}

template<typename T, typename U>
U
maxValueField(const Matrix2<T>& data)
{
	U max = minValue<U>();
	for (typename Matrix2<T>::ConstIterator iter = data.begin(); iter != data.end(); ++iter) {
		const U v = iter->value;
		if (v > max) max = v;
	}
	return max;
}

template<typename T, typename U>
U
maxAbsoluteValueField(const Matrix2<T>& data)
{
	U max = 0;
	for (typename Matrix2<T>::ConstIterator iter = data.begin(); iter != data.end(); ++iter) {
		const U v = std::abs(iter->value);
		if (v > max) max = v;
	}
	return max;
}

template<typename T>
void
minMax(const std::vector<T>& list, T& min, T& max)
{
	min = std::numeric_limits<T>::max();
	max = minValue<T>();
	for (typename std::vector<T>::const_iterator iter = list.begin(); iter != list.end(); ++iter) {
		const T a = *iter;
		if (max < a) max = a;
		if (min > a) min = a;
	}
}

template<typename T, typename U>
void
minMaxValueField(const Matrix2<T>& data, U& min, U& max)
{
	min = std::numeric_limits<U>::max();
	max = minValue<U>();
	for (typename Matrix2<T>::ConstIterator iter = data.begin(); iter != data.end(); ++iter) {
		const U v = iter->value;
		if (v > max) max = v;
		if (v < min) min = v;
	}
}

template<typename T>
void
multiply(std::vector<T>& list, T coefficient)
{
	for (typename std::vector<T>::iterator iter = list.begin(); iter != list.end(); ++iter) {
		*iter *= coefficient;
	}
}

template<typename T, typename U>
void
multiply(const std::vector<T>& factorList, std::vector<U>& data)
{
	if (factorList.size() != data.size()) {
		THROW_EXCEPTION(InvalidParameterException, "The input vectors have different sizes.");
	}
	typename std::vector<T>::const_iterator factorIter = factorList.begin();
	typename std::vector<T>::const_iterator factorEnd = factorList.end();
	typename std::vector<U>::iterator dataIter = data.begin();
	while (factorIter != factorEnd) {
		*dataIter++ *= *factorIter++;
	}
}

template<typename T>
void
centralDiff(const std::vector<T>& inputList, T period, std::vector<T>& outputList)
{
	const std::size_t inputSize = inputList.size();
	if (inputSize == 0) {
		THROW_EXCEPTION(InvalidParameterException, "inputList is empty.");
	}
	if (inputSize > (std::numeric_limits<unsigned int>::max() / sizeof(T)) / 2) { // arbitrary limit
		THROW_EXCEPTION(InvalidParameterException, "inputList is too long.");
	}
	if (period == 0) {
		THROW_EXCEPTION(InvalidParameterException, "period = 0.");
	}
	outputList.resize(inputSize + 2);

	const T coeff = T(0.5) / period;
	outputList[0] = coeff * inputList[0];
	outputList[1] = coeff * inputList[1];
	for (std::size_t i = 2, end = outputList.size() - 2; i < end; ++i) {
		outputList[i] = coeff * (inputList[i] - inputList[i - 2]);
	}
	outputList[outputList.size() - 2] = -coeff * inputList[inputSize - 2];
	outputList[outputList.size() - 1] = -coeff * inputList[inputSize - 1];
}

//template<typename FloatType, typename IntType, typename InputIterator, typename OutputIterator>
//void
//copyFloatToInt(InputIterator input, InputIterator inputEnd, OutputIterator output)
//{
//	for (; input != inputEnd; ++input, ++output) {
//		if (*input >= 0.0) {
//			*output = static_cast<IntType>(FloatType(0.5) + *input);
//		} else {
//			*output = static_cast<IntType>(FloatType(-0.5) + *input);
//		}
//	}
//}

template<typename T>
T
decibelsToLinear(T decibels) {
	return std::pow(T(10), decibels * T(0.05));
}

template<typename T>
T
linearToDecibels(T linear) {
	return T(20) * std::log10(linear);
}

template<typename T>
T
degreeToRadian(T d)
{
	return T(PI) * d / T(180);
}

template<typename T>
T
radianToDegree(T r)
{
	return T(180) * r / T(PI);
}

template<typename InputIterator, typename InputOutputIterator>
void
addElements(InputIterator first1, InputOutputIterator first2, InputOutputIterator last2)
{
	while (first2 != last2) {
		*first2++ += *first1++;
	}
}

template<typename InputIterator1, typename InputIterator2, typename OutputIterator>
void
addElements(InputIterator1 first1,
	    InputIterator2 first2,
	    OutputIterator first3, OutputIterator last3)
{
	while (first3 != last3) {
		*first3++ = *first1++ + *first2++;
	}
}

template<typename T, typename U>
void
copyXZToSimpleMatrices(const Matrix2<T>& m, Matrix2<U>& mX, Matrix2<U>& mZ)
{
	mX.resize(m.n1(), m.n2());
	mZ.resize(m.n1(), m.n2());
	typename Matrix2<T>::ConstIterator orig = m.begin();
	typename Matrix2<U>::Iterator destX = mX.begin();
	typename Matrix2<U>::Iterator destZ = mZ.begin();
	while (orig != m.end()) {
		*destX++ = orig->x;
		*destZ++ = orig->z;
		++orig;
	}
}

template<typename T, typename U>
void
copyXZFromSimpleMatrices(const Matrix2<T>& mX, const Matrix2<T>& mZ, Matrix2<U>& m)
{
	if (mX.n1() != mZ.n1() || mX.n2() != mZ.n2()
			|| mX.n1() != m.n1() || mX.n2() != m.n2()) {
		THROW_EXCEPTION(InvalidParameterException, "The input and output matrices must have the same sizes.");
	}
	typename Matrix2<T>::ConstIterator origX = mX.begin();
	typename Matrix2<T>::ConstIterator origZ = mZ.begin();
	typename Matrix2<U>::Iterator dest = m.begin();
	for ( ; origX != mX.end(); ++origX, ++origZ, ++dest) {
		dest->x = *origX;
		dest->z = *origZ;
	}
}

template<typename T, typename U>
void
copyValueToSimpleMatrix(const Matrix2<T>& m, Matrix2<U>& mV)
{
	mV.resize(m.n1(), m.n2());
	typename Matrix2<T>::ConstIterator orig = m.begin();
	typename Matrix2<U>::Iterator dest = mV.begin();
	while (orig != m.end()) {
		*dest++ = (orig++)->value;
	}
}

template<typename T, typename U>
void
copyValueFromSimpleMatrix(const Matrix2<T>& mV, Matrix2<U>& m)
{
	if (mV.n1() != m.n1() || mV.n2() != m.n2()) {
		THROW_EXCEPTION(InvalidParameterException, "The input and output matrices must have the same sizes.");
	}
	typename Matrix2<T>::ConstIterator orig = mV.begin();
	typename Matrix2<U>::Iterator dest = m.begin();
	while (orig != mV.end()) {
		(dest++)->value = *orig++;
	}
}

template<typename T, typename U>
void
copyFactorToSimpleMatrix(const Matrix2<T>& m, Matrix2<U>& mF)
{
	mF.resize(m.n1(), m.n2());
	typename Matrix2<T>::ConstIterator orig = m.begin();
	typename Matrix2<U>::Iterator dest = mF.begin();
	while (orig != m.end()) {
		*dest++ = (orig++)->factor;
	}
}

template<typename T, typename U>
void
copyXZToSimpleVectors(const std::vector<T>& v, std::vector<U>& vX, std::vector<U>& vZ)
{
	vX.resize(v.size());
	vZ.resize(v.size());
	typename std::vector<T>::const_iterator orig = v.begin();
	typename std::vector<U>::iterator destX = vX.begin();
	typename std::vector<U>::iterator destZ = vZ.begin();
	for ( ; orig != v.end(); ++orig) {
		*destX++ = orig->x;
		*destZ++ = orig->z;
	}
}

template<typename T, typename U>
void
copyXZFromSimpleVectors(const std::vector<T>& vX, const std::vector<T>& vZ, std::vector<U>& v)
{
	if (vX.size() != vZ.size() || vX.size() != v.size()) {
		THROW_EXCEPTION(InvalidParameterException, "The input and output vectors must have the same sizes.");
	}
	typename std::vector<T>::const_iterator origX = vX.begin();
	typename std::vector<T>::const_iterator origZ = vZ.begin();
	typename std::vector<U>::iterator dest = v.begin();
	for ( ; origX != vX.end(); ++origX, ++origZ, ++dest) {
		dest->x = *origX;
		dest->z = *origZ;
	}
}

template<typename T, typename U>
void
copyXZValue(const Matrix2<T>& orig, Matrix2<U>& dest)
{
	dest.resize(orig.n1(), orig.n2());
	typename Matrix2<T>::ConstIterator origIter = orig.begin();
	typename Matrix2<U>::Iterator destIter = dest.begin();
	while (origIter != orig.end()) {
		destIter->x     = origIter->x;
		destIter->z     = origIter->z;
		destIter->value = origIter->value;
		++origIter;
		++destIter;
	}
}

template<typename T, typename U>
void
copyXZFactor(const Matrix2<T>& orig, Matrix2<U>& dest)
{
	dest.resize(orig.n1(), orig.n2());
	typename Matrix2<T>::ConstIterator origIter = orig.begin();
	typename Matrix2<U>::Iterator destIter = dest.begin();
	while (origIter != orig.end()) {
		destIter->x     = origIter->x;
		destIter->z     = origIter->z;
		destIter->value = origIter->factor;
		++origIter;
		++destIter;
	}
}

template<typename T, typename U>
void
copyXZ(const Matrix2<T>& orig, Matrix2<U>& dest)
{
	dest.resize(orig.n1(), orig.n2());
	typename Matrix2<T>::ConstIterator origIter = orig.begin();
	typename Matrix2<U>::Iterator destIter = dest.begin();
	while (origIter != orig.end()) {
		destIter->x     = origIter->x;
		destIter->z     = origIter->z;
		++origIter;
		++destIter;
	}
}

template<typename T, typename U>
void
copyXZ(const std::vector<T>& orig, std::vector<U>& dest)
{
	dest.resize(orig.size());
	typename std::vector<T>::const_iterator origIter = orig.begin();
	typename std::vector<U>::iterator destIter = dest.begin();
	while (origIter != orig.end()) {
		destIter->x = origIter->x;
		destIter->z = origIter->z;
		++origIter;
		++destIter;
	}
}

// This template has specializations for float types.
template<typename T>
T
minValue()
{
	return std::numeric_limits<T>::min();
}

template<typename T>
T
maxValue()
{
	return std::numeric_limits<T>::max();
}

template<typename T>
T
nextPowerOf2(T v)
{
	v = std::abs(v);
	const T limit = std::numeric_limits<T>::max() / 2;
	T n = 1;
	while (n < v) {
		if (n > limit) {
			THROW_EXCEPTION(InvalidParameterException, "The input value is too big: " << v << '.');
		}
		n *= 2;
	}
	return n;
}

template<typename T>
T
multiplyElements(const std::vector<T>& v)
{
	T res = 1;
	for (typename std::vector<T>::const_iterator iter = v.begin(); iter != v.end(); ++iter) {
		res *= *iter;
	}
	return res;
}

struct CopyAbsValueOp {
	template<typename T, typename U>
	void operator()(const T& orig, U& dest) {
		dest = std::abs(orig.value);
	}
};

struct CopyValueOp {
	template<typename T, typename U>
	void operator()(const T& orig, U& dest) {
		dest = orig.value;
	}
};

struct CopyToValueOp {
	template<typename T, typename U>
	void operator()(const T& orig, U& dest) {
		dest.value = orig;
	}
};

struct CopyFactorOp {
	template<typename T, typename U>
	void operator()(const T& orig, U& dest) {
		dest = orig.factor;
	}
};

struct CopyToFactorOp {
	template<typename T, typename U>
	void operator()(const T& orig, U& dest) {
		dest.factor = orig;
	}
};

struct CopyXOp {
	template<typename T, typename U>
	void operator()(const T& orig, U& dest) {
		dest = orig.x;
	}
};

struct CopyToXOp {
	template<typename T, typename U>
	void operator()(const T& orig, U& dest) {
		dest.x = orig;
	}
};

struct CopyZOp {
	template<typename T, typename U>
	void operator()(const T& orig, U& dest) {
		dest = orig.z;
	}
};

struct CopyToZOp {
	template<typename T, typename U>
	void operator()(const T& orig, U& dest) {
		dest.z = orig;
	}
};

template<typename InputIterator, typename OutputIterator, typename T>
void
copyUsingOperator(InputIterator orig, InputIterator origEnd, OutputIterator dest, T copyOp)
{
	while (orig != origEnd) {
		copyOp(*orig++, *dest++);
	}
}

template<typename FloatType>
void
generateArc2D(FloatType centerX, FloatType centerZ,
		FloatType radius, FloatType angle1, FloatType angle2,
		FloatType arcStepFactor, FloatType lambda,
		std::vector<XZ<FloatType> >& pointList)
{
	if (radius == 0) THROW_EXCEPTION(InvalidParameterException, "radius = 0.");

	const FloatType angleStep = ((angle2 < angle1) ? -arcStepFactor : arcStepFactor) * lambda / radius;
	std::vector<FloatType> angleList;
	fillSequence(angleList, angle1, angle2, angleStep);
	pointList.resize(angleList.size());
	for (std::size_t i = 0; i < pointList.size(); ++i) {
		const FloatType ang = angleList[i];
		pointList[i] = XZ<FloatType>(centerX + radius * std::cos(ang), centerZ + radius * std::sin(ang));
	}
}

template<typename FloatType>
void
calculateDelays2D(FloatType propagationSpeed, XZ<FloatType> focus,
			const std::vector<XZ<FloatType> >& arrayElemPosList, unsigned int numElements,
			std::vector<FloatType>& delays)
{
	if (propagationSpeed == 0) THROW_EXCEPTION(InvalidParameterException, "propagationSpeed = 0.");
	if (arrayElemPosList.empty()) THROW_EXCEPTION(InvalidParameterException, "arrayElemPosList is empty.");

	if (delays.size() != numElements) delays.resize(numElements);

	for (unsigned int elem = 0; elem < numElements; ++elem) {
		const FloatType dx = focus.x - arrayElemPosList[elem].x;
		const FloatType dz = focus.z - arrayElemPosList[elem].z;
		const FloatType r = std::sqrt(dx * dx + dz * dz);
		const FloatType dt = r / propagationSpeed;
		delays[elem] = dt;
	}
	const FloatType maxDelay = *std::max_element(delays.begin(), delays.end());
	for (unsigned int elem = 0; elem < numElements; ++elem) {
		delays[elem] = maxDelay - delays[elem];
	}
}

template<typename Iterator>
void
resetValueFactor(Iterator iter, Iterator iterEnd)
{
	while (iter != iterEnd) {
		iter->value = 0;
		iter->factor = 1;
		++iter;
	}
}

template<typename T>
void
removeDC(T* data, std::size_t size)
{
	const T sum = std::accumulate(data, data + size, T(0));
	std::for_each(data, data + size, Add<T>(-sum / size));
}

template<typename T>
void
removeDC(T* data, std::size_t size, std::size_t beginOffset)
{
	if (beginOffset >= size) {
		THROW_EXCEPTION(InvalidParameterException, "Invalid beginOffset (" << beginOffset << " >= " << size << ").");
	}
	const T sum = std::accumulate(data + beginOffset, data + size, T(0)); // doesn't use the entire data
	std::for_each(data, data + size, Add<T>(-sum / (size - beginOffset))); // applies to the entire data
}

} // namespace Util
} // namespace Lab

#endif /* UTIL_H_ */
