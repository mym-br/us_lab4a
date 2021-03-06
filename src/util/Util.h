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

#ifndef UTIL_H_
#define UTIL_H_

#include <algorithm> /* for_each, max_element */
#include <chrono>
#include <cmath>
#include <complex>
#include <cstddef> /* std::size_t */
#include <cstring>
#include <limits>
#include <numeric> /* accumulate */
#include <ratio>
#include <thread>
#include <vector>

#include "Exception.h"
#include "Matrix.h"
#include "Tensor3.h"
#include "XYZValue.h"
#include "XYZValueArray.h"
#include "XZValue.h"

constexpr double pi = M_PI;



namespace Lab {
namespace Util {

template<typename T> T wavelength(T speed, T frequency);
// The minimum sampling rate without aliasing.
template<typename T> T nyquistRate(T maxFrequency);
// Wavelength at Nyquist rate.
template<typename T> T nyquistLambda(T speed, T maxFrequency);

void fillSequence(std::vector<unsigned int>& v, unsigned int startValue, unsigned int endValue, int step = 1);
void fillSequence(std::vector<int>& v, int startValue, int endValue, int step = 1);

// This function may not give good results with integer values.
// The sequence will contain _startValue_ and _endValue_. The effective step may be less than the requested.
template<typename T> void fillSequenceFromStartToEndWithMaximumStep(std::vector<T>& v, T startValue, T endValue, T step = 1);

// The sequence will contain _startValue_ and may not contain _endValue_.
template<typename T> void fillSequenceFromStartWithStep(std::vector<T>& v, T startValue, T endValue, T step = 1);

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
template<typename T> T maxAbsolute(const std::vector<std::complex<T>>& list);
template<typename T> T maxAbsolute(const Matrix<T>& data);
template<typename T, typename U> U maxValueField(const Matrix<T>& data);
template<typename T, typename U> U maxAbsoluteValueField(const Matrix<T>& data);
template<typename T> T maxAbsoluteValueField(const Matrix<XYZValueArray<T>>& data);
template<typename T> void minMax(const std::vector<T>& list, T& minVal, T& maxVal);
template<typename T, typename U> void minMaxValueField(const Matrix<T>& data, U& minVal, U& maxVal);
template<typename T> void multiply(std::vector<T>& list, T coefficient);
template<typename T> void multiply(Matrix<T>& data, T coefficient);
template<typename T, typename U> void multiply(const std::vector<T>& factorList, std::vector<U>& data);
template<typename T> void centralDiff(const std::vector<T>& inputList, T period, std::vector<T>& outputList);
template<typename T> void deleteObjects(T& container);
template<typename T> T decibelsToLinear(T decibels);
template<typename T> T linearToDecibels(T linear);
template<typename T> void linearToDecibels(std::vector<T>& data, T minDecibels);
template<typename T> T degreeToRadian(T d);
template<typename T> T radianToDegree(T r);

// Add the values pointed by first1 to the elements pointed by first2.
template<typename InputIterator, typename InputOutputIterator>
void addElements(InputIterator first1, InputOutputIterator first2, InputOutputIterator last2);

// Add the values pointed by first1 to the values pointed by first2
// and send the results to the elements pointed by first3.
template<typename InputIterator1, typename InputIterator2, typename OutputIterator>
void addElements(InputIterator1 first1, InputIterator2 first2, OutputIterator first3, OutputIterator last3);

template<typename T, typename U> void copyXYZToSimpleMatrices(const Matrix<T>& m, Matrix<U>& mX, Matrix<U>& mZ);
template<typename T, typename U> void copyXYZFromSimpleMatrices(const Matrix<T>& mX, const Matrix<T>& mZ, Matrix<U>& m);

template<typename T, typename U> void copyValueToSimpleMatrix(const Matrix<T>& m, Matrix<U>& mV);
template<typename T, typename U> void copyValueFromSimpleMatrix(const Matrix<T>& mV, Matrix<U>& m);

template<typename T, typename U> void copyFactorToSimpleMatrix(const Matrix<T>& m, Matrix<U>& mF);

template<typename T, typename U> void copyXYZToSimpleVectors(const std::vector<T>& v,
								std::vector<U>& vX,
								std::vector<U>& vY,
								std::vector<U>& vZ);
template<typename T, typename U> void copyXYZFromSimpleVectors(const std::vector<T>& vX,
								const std::vector<T>& vY,
								const std::vector<T>& vZ,
								std::vector<U>& v);

template<typename T, typename U> void copyXYZ(const Matrix<T>& orig, Matrix<U>& dest);
template<typename T, typename U> void copyXYZValue(const Matrix<T>& orig, Matrix<U>& dest);
template<typename T, typename U> void copyXYZFactor(const Matrix<T>& orig, Matrix<U>& dest);
template<typename T, typename U> void copyXZ(const Matrix<T>& orig, Matrix<U>& dest);
template<typename T, typename U> void copyXZ(const std::vector<T>& orig, std::vector<U>& dest);
template<typename T, typename U> void copy(const Matrix<XZValue<T>>& orig, Matrix<XYZValue<U>>& dest);

template<typename T> T minValue();
template<typename T> T maxValue();

template<typename T> T nextPowerOf2(T v);

template<typename T> T multiplyElements(const std::vector<T>& v);

template<typename InputIterator, typename OutputIterator, typename T> void copyUsingOperator(InputIterator orig, InputIterator origEnd, OutputIterator dest, T copyOp);

template<typename Iterator> void resetValueFactor(Iterator iter, Iterator iterEnd);

template<typename T> void removeDC(T* data, std::size_t size);
template<typename T> void removeDC(T* data, std::size_t size, std::size_t beginOffset);

template<typename T> T sign(T value);

template<typename T> void normalize(T& data);
template<typename T> void normalizeBySumOfAbs(std::vector<T>& data);

unsigned int numberOfDigits(unsigned int value);

template<typename Iterator> void applyFactorToValue(Iterator iter, Iterator iterEnd);

template<typename T, typename U> std::size_t sizeInBytes(const std::vector<T, U>& v);
template<typename T, typename U> std::size_t sizeInBytes(const Matrix<T, U>& m);
template<typename T, typename U> std::size_t sizeInBytes(const Tensor3<T, U>& t);

//#############################################################################

template<typename T>
struct MultiplyBy {
	explicit MultiplyBy(T factor) : factor(factor) {}
	void operator()(T& value) { value *= factor; }
	T factor;
};

template<typename T, typename U>
struct MultiplyValueBy {
	explicit MultiplyValueBy(U factor) : factor(factor) {}
	void operator()(T& item) { item.value *= factor; }
	U factor;
};

template<typename T>
struct Add {
	explicit Add(T term) : term(term) {}
	void operator()(T& value) { value += term; }
	T term;
};

template<typename T>
struct ValueFieldAbsOp {
	ValueFieldAbsOp() = default;
	void operator()(T& o) { o.value = std::abs(o.value); }
};

template<typename T>
T
wavelength(T speed, T frequency)
{
	return speed / frequency;
}

template<typename T>
T
nyquistRate(T maxFrequency)
{
	return 2 * maxFrequency;
}

template<typename T>
T
nyquistLambda(T speed, T maxFrequency)
{
	return wavelength(speed, nyquistRate(maxFrequency));
}

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
	if (n == 1) {
		THROW_EXCEPTION(InvalidParameterException, "Sequence with only one element.");
	}
	const double newStep = (endValue - startValue) / (n - 1);
	v.resize(n);
	for (int i = 0; i < n; ++i) {
		v[i] = startValue + i * newStep;
	}
}

template<typename T>
void
fillSequenceFromStartWithStep(std::vector<T>& v, T startValue, T endValue, T step)
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

	v.clear();
	T value = startValue;
	std::size_t i = 0;
	if (step > 0) {
		while (value <= endValue) {
			v.push_back(value);
			++i;
			value = startValue + i * step;
		}
	} else {
		while (value >= endValue) {
			v.push_back(value);
			++i;
			value = startValue + i * step;
		}
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
	std::this_thread::sleep_for(std::chrono::duration<unsigned long, std::milli>(milliseconds));
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
	T maxVal = minValue<T>();

	for (auto iter = list.cbegin(); iter != list.cend(); ++iter) {
		if (maxVal < *iter) {
			maxVal = *iter;
		}
	}

	return maxVal;
}

template<typename T>
T
maxAbsolute(const std::vector<T>& list)
{
	T maxVal = 0;
	for (auto iter = list.cbegin(); iter != list.cend(); ++iter) {
		const T a = std::abs(*iter);
		if (maxVal < a) maxVal = a;
	}
	return maxVal;
}

template<typename T>
T
maxAbsolute(const std::vector<std::complex<T>>& list)
{
	T maxVal = 0;
	for (auto iter = list.cbegin(); iter != list.cend(); ++iter) {
		const T a = std::abs(*iter);
		if (maxVal < a) maxVal = a;
	}
	return maxVal;
}

template<typename T>
T
maxAbsolute(const Matrix<T>& data)
{
	T maxVal = 0;
	for (auto iter = data.cbegin(); iter != data.cend(); ++iter) {
		const T a = std::abs(*iter);
		if (maxVal < a) maxVal = a;
	}
	return maxVal;
}

template<typename T, typename U>
U
maxValueField(const Matrix<T>& data)
{
	U maxVal = minValue<U>();
	for (auto iter = data.cbegin(); iter != data.cend(); ++iter) {
		const U v = iter->value;
		if (v > maxVal) maxVal = v;
	}
	return maxVal;
}

template<typename T, typename U>
U
maxAbsoluteValueField(const Matrix<T>& data)
{
	U maxVal = 0;
	for (auto iter = data.cbegin(); iter != data.cend(); ++iter) {
		const U a = std::abs(iter->value);
		if (maxVal < a) maxVal = a;
	}
	return maxVal;
}

template<typename T>
T
maxAbsoluteValueField(const Matrix<XYZValueArray<T>>& data)
{
	T maxVal = 0;
	for (auto iter = data.cbegin(); iter != data.cend(); ++iter) {
		for (auto v : iter->values) {
			const T a = std::abs(v);
			if (maxVal < a) maxVal = a;
		}
	}
	return maxVal;
}

template<typename T>
void
minMax(const std::vector<T>& list, T& minVal, T& maxVal)
{
	minVal = std::numeric_limits<T>::max();
	maxVal = minValue<T>();
	for (auto iter = list.cbegin(); iter != list.cend(); ++iter) {
		const T a = *iter;
		if (maxVal < a) maxVal = a;
		if (minVal > a) minVal = a;
	}
}

template<typename T, typename U>
void
minMaxValueField(const Matrix<T>& data, U& minVal, U& maxVal)
{
	minVal = std::numeric_limits<U>::max();
	maxVal = minValue<U>();
	for (auto iter = data.cbegin(); iter != data.cend(); ++iter) {
		const U v = iter->value;
		if (v > maxVal) maxVal = v;
		if (v < minVal) minVal = v;
	}
}

template<typename T>
void
multiply(std::vector<T>& list, T coefficient)
{
	for (auto iter = list.begin(); iter != list.end(); ++iter) {
		*iter *= coefficient;
	}
}

template<typename T>
void
multiply(Matrix<T>& data, T coefficient)
{
	for (auto iter = data.begin(); iter != data.end(); ++iter) {
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
	auto factorIter = factorList.cbegin();
	auto factorEnd = factorList.cend();
	auto dataIter = data.begin();
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
void
linearToDecibels(std::vector<T>& data, T minDecibels)
{
	normalize(data);
	auto minVal = decibelsToLinear(minDecibels);
	for (auto& value : data) {
		if (value <= minVal) {
			value = minDecibels;
		} else {
			value = linearToDecibels(value);
		}
	}
}

template<typename T>
T
degreeToRadian(T d)
{
	return T(pi) * d / T(180);
}

template<typename T>
T
radianToDegree(T r)
{
	return T(180) * r / T(pi);
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
copyXYZToSimpleMatrices(const Matrix<T>& m, Matrix<U>& mX, Matrix<U>& mY, Matrix<U>& mZ)
{
	mX.resize(m.n1(), m.n2());
	mY.resize(m.n1(), m.n2());
	mZ.resize(m.n1(), m.n2());
	auto orig = m.cbegin();
	auto destX = mX.begin();
	auto destY = mY.begin();
	auto destZ = mZ.begin();
	while (orig != m.end()) {
		*destX++ = orig->x;
		*destY++ = orig->y;
		*destZ++ = orig->z;
		++orig;
	}
}

template<typename T, typename U>
void
copyXYZFromSimpleMatrices(const Matrix<T>& mX, const Matrix<T>& mY, const Matrix<T>& mZ, Matrix<U>& m)
{
	if (mX.n1() != mY.n1() || mX.n2() != mY.n2() || mX.n1() != mZ.n1() || mX.n2() != mZ.n2()
			|| mX.n1() != m.n1() || mX.n2() != m.n2()) {
		THROW_EXCEPTION(InvalidParameterException, "The input and output matrices must have the same sizes.");
	}
	auto origX = mX.cbegin();
	auto origY = mY.cbegin();
	auto origZ = mZ.cbegin();
	auto dest = m.begin();
	for ( ; origX != mX.end(); ++origX, ++origY, ++origZ, ++dest) {
		dest->x = *origX;
		dest->y = *origY;
		dest->z = *origZ;
	}
}

template<typename T, typename U>
void
copyValueToSimpleMatrix(const Matrix<T>& m, Matrix<U>& mV)
{
	mV.resize(m.n1(), m.n2());
	auto orig = m.cbegin();
	auto dest = mV.begin();
	while (orig != m.end()) {
		*dest++ = (orig++)->value;
	}
}

template<typename T, typename U>
void
copyValueFromSimpleMatrix(const Matrix<T>& mV, Matrix<U>& m)
{
	if (mV.n1() != m.n1() || mV.n2() != m.n2()) {
		THROW_EXCEPTION(InvalidParameterException, "The input and output matrices must have the same sizes.");
	}
	auto orig = mV.cbegin();
	auto dest = m.begin();
	while (orig != mV.end()) {
		(dest++)->value = *orig++;
	}
}

template<typename T, typename U>
void
copyFactorToSimpleMatrix(const Matrix<T>& m, Matrix<U>& mF)
{
	mF.resize(m.n1(), m.n2());
	auto orig = m.cbegin();
	auto dest = mF.begin();
	while (orig != m.end()) {
		*dest++ = (orig++)->factor;
	}
}

template<typename T, typename U>
void
copyXYZToSimpleVectors(const std::vector<T>& v,
			std::vector<U>& vX,
			std::vector<U>& vY,
			std::vector<U>& vZ)
{
	vX.resize(v.size());
	vY.resize(v.size());
	vZ.resize(v.size());
	auto orig = v.cbegin();
	auto destX = vX.begin();
	auto destY = vY.begin();
	auto destZ = vZ.begin();
	for ( ; orig != v.end(); ++orig) {
		*destX++ = orig->x;
		*destY++ = orig->y;
		*destZ++ = orig->z;
	}
}

template<typename T, typename U>
void
copyXYZFromSimpleVectors(const std::vector<T>& vX,
				const std::vector<T>& vY,
				const std::vector<T>& vZ,
				std::vector<U>& v)
{
	if (vX.size() != vY.size() || vX.size() != vZ.size() || vX.size() != v.size()) {
		THROW_EXCEPTION(InvalidParameterException, "The input and output vectors must have the same sizes.");
	}
	auto origX = vX.cbegin();
	auto origY = vY.cbegin();
	auto origZ = vZ.cbegin();
	auto dest = v.begin();
	for ( ; origX != vX.end(); ++origX, ++origY, ++origZ, ++dest) {
		dest->x = *origX;
		dest->y = *origY;
		dest->z = *origZ;
	}
}

template<typename T, typename U>
void
copyXYZ(const Matrix<T>& orig, Matrix<U>& dest)
{
	dest.resize(orig.n1(), orig.n2());
	auto origIter = orig.cbegin();
	auto destIter = dest.begin();
	while (origIter != orig.end()) {
		destIter->x     = origIter->x;
		destIter->y     = origIter->y;
		destIter->z     = origIter->z;
		++origIter;
		++destIter;
	}
}

template<typename T, typename U>
void
copyXYZValue(const Matrix<T>& orig, Matrix<U>& dest)
{
	dest.resize(orig.n1(), orig.n2());
	auto origIter = orig.cbegin();
	auto destIter = dest.begin();
	while (origIter != orig.end()) {
		destIter->x     = origIter->x;
		destIter->y     = origIter->y;
		destIter->z     = origIter->z;
		destIter->value = origIter->value;
		++origIter;
		++destIter;
	}
}

template<typename T, typename U>
void
copyXYZFactor(const Matrix<T>& orig, Matrix<U>& dest)
{
	dest.resize(orig.n1(), orig.n2());
	auto origIter = orig.cbegin();
	auto destIter = dest.begin();
	while (origIter != orig.end()) {
		destIter->x     = origIter->x;
		destIter->y     = origIter->y;
		destIter->z     = origIter->z;
		destIter->value = origIter->factor;
		++origIter;
		++destIter;
	}
}

template<typename T, typename U>
void
copyXZ(const Matrix<T>& orig, Matrix<U>& dest)
{
	dest.resize(orig.n1(), orig.n2());
	auto origIter = orig.cbegin();
	auto destIter = dest.begin();
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
	auto origIter = orig.cbegin();
	auto destIter = dest.begin();
	while (origIter != orig.end()) {
		destIter->x = origIter->x;
		destIter->z = origIter->z;
		++origIter;
		++destIter;
	}
}

template<typename T, typename U>
void
copy(const Matrix<XZValue<T>>& orig, Matrix<XYZValue<U>>& dest)
{
	dest.resize(orig.n1(), orig.n2());
	auto origIter = orig.cbegin();
	auto destIter = dest.begin();
	while (origIter != orig.end()) {
		destIter->x = origIter->x;
		destIter->y = 0.0;
		destIter->z = origIter->z;
		++origIter;
		++destIter;
	}
}

template<typename T>
T
minValue()
{
	return std::numeric_limits<T>::lowest();
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
	for (auto iter = v.cbegin(); iter != v.cend(); ++iter) {
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

struct CopyYOp {
	template<typename T, typename U>
	void operator()(const T& orig, U& dest) {
		dest = orig.y;
	}
};

struct CopyToYOp {
	template<typename T, typename U>
	void operator()(const T& orig, U& dest) {
		dest.y = orig;
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

template<typename T>
T
sign(T value)
{
	if (value > 0) return 1;
	if (value < 0) return -1;
	return 0;
}

template<typename T>
void
normalize(T& data)
{
	const auto maxAbs = maxAbsolute(data);
	if (maxAbs == 0) return;
	const auto coeff = 1 / maxAbs;
	multiply(data, coeff);
}

template<typename T>
void
normalizeBySumOfAbs(std::vector<T>& data)
{
	T sumAbs = 0;
	for (const auto item : data) {
		sumAbs += std::abs(item);
	}
	if (sumAbs == 0) return;
	const auto coeff = 1 / sumAbs;
	multiply(data, coeff);
}

template<typename Iterator>
void
applyFactorToValue(Iterator iter, Iterator iterEnd)
{
	for ( ; iter != iterEnd; ++iter) {
		iter->value *= iter->factor;
		iter->factor = 1;
	}
}

template<typename T, typename U>
std::size_t
sizeInBytes(const std::vector<T, U>& v)
{
	return v.size() * sizeof(typename std::vector<T, U>::value_type);
}

template<typename T, typename U>
std::size_t
sizeInBytes(const Matrix<T, U>& m)
{
	return m.size() * sizeof(typename Matrix<T, U>::ValueType);
}

template<typename T, typename U>
std::size_t
sizeInBytes(const Tensor3<T, U>& t)
{
	return t.size() * sizeof(typename Tensor3<T, U>::ValueType);
}

} // namespace Util
} // namespace Lab

#endif /* UTIL_H_ */
