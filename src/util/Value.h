#ifndef VALUE_H_
#define VALUE_H_

#include <cmath> /* abs, sqrt */
#include <complex>

#include "XZValue.h"
#include "XZValueFactor.h"



namespace Lab {
namespace Value {

template<typename FloatType> void copy(const XZValue<FloatType>& orig, FloatType& dest);
template<typename FloatType> void copy(const FloatType& orig, XZValue<FloatType>& dest);

template<typename FloatType> void copy(const XZValueFactor<FloatType>& orig, FloatType& dest);
template<typename FloatType> void copy(const FloatType& orig, XZValueFactor<FloatType>& dest);

// array of size two --> complex
template<typename FloatType> void copy(const FloatType (&orig)[2], std::complex<FloatType>& dest);
// complex --> array of size two
template<typename FloatType> void copy(const std::complex<FloatType>& orig, FloatType (&dest)[2]);

template<typename FloatType1, typename FloatType2> void copy(const FloatType1& orig, FloatType2& dest);
template<typename FloatType1, typename FloatType2> void copy(const std::complex<FloatType1>& orig, FloatType2 &dest);

template<typename InputIterator, typename OutputIterator>
	void copySequence(InputIterator input, InputIterator inputEnd, OutputIterator output);
template<typename InputIterator, typename OutputIterator>
	void copySequenceWithPadding(InputIterator input, InputIterator inputEnd, OutputIterator output, unsigned int padding);
template<typename InputIterator, typename OutputIterator, typename T>
	void transformSequence(InputIterator orig, InputIterator origEnd, OutputIterator dest, T transform);
template<typename Iterator, typename FloatType>
	void copyRealImagToComplexSequence(Iterator re, Iterator reEnd, Iterator im, std::vector<std::complex<FloatType> >& cpx);
template<typename InputIterator, typename OutputIterator>
	void copyXZValueSequence(InputIterator orig, InputIterator origEnd, OutputIterator dest);
template<typename InputIterator, typename OutputIterator>
	void copyXZSequence(InputIterator orig, InputIterator origEnd, OutputIterator dest);



template<typename FloatType>
class ComplexToScaledAbsoluteOp {
public:
	ComplexToScaledAbsoluteOp(FloatType scale) : scale_(scale) { }

	template<typename DestType>
	void operator()(const std::complex<FloatType>& orig, DestType& dest) {
		copy(std::sqrt(orig.real() * orig.real() + orig.imag() * orig.imag()) * scale_, dest);
	}

	template<typename DestType>
	void operator()(const FloatType (&orig)[2], DestType& dest) {
		copy(std::sqrt(orig[0] * orig[0] + orig[1] * orig[1]) * scale_, dest);
	}
private:
	FloatType scale_;
};

template<typename T, typename U>
struct ScalarToValueFieldOp {
	void operator()(T& orig, U& dest) {
		dest.value = orig;
	}
};

template<typename FloatType>
class ScaleComplexOp {
public:
	ScaleComplexOp(FloatType scale) : scale_(scale) { }

	void operator()(const std::complex<FloatType>& orig, std::complex<FloatType>& dest) {
		dest = orig * scale_;
	}
	void operator()(const FloatType (&orig)[2], std::complex<FloatType>& dest) {
		dest = std::complex<FloatType>(orig[0], orig[1]) * scale_;
	}
	void operator()(const std::complex<FloatType>& orig, FloatType (&dest)[2]) {
		dest[0] = orig.real() * scale_;
		dest[0] = orig.imag() * scale_;
	}
	void operator()(const FloatType (&orig)[2], FloatType (&dest)[2]) {
		dest[0] = orig[0] * scale_;
		dest[0] = orig[1] * scale_;
	}
private:
	FloatType scale_;
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



template<typename FloatType>
void
copy(const XZValue<FloatType>& orig, FloatType& dest)
{
	dest = orig.value;
}

template<typename FloatType>
void
copy(const FloatType& orig, XZValue<FloatType>& dest)
{
	dest.value = orig;
}

template<typename FloatType>
void
copy(const XZValueFactor<FloatType>& orig, FloatType& dest)
{
	dest = orig.value;
}

template<typename FloatType>
void
copy(const FloatType& orig, XZValueFactor<FloatType>& dest)
{
	dest.value = orig;
}

template<typename FloatType>
void
copy(const FloatType (&orig)[2], std::complex<FloatType>& dest)
{
	dest = std::complex<FloatType>(orig[0], orig[1]);
}

template<typename FloatType>
void
copy(const std::complex<FloatType>& orig, FloatType (&dest)[2])
{
	dest[0] = orig.real();
	dest[1] = orig.imag();
}


template<typename FloatType1, typename FloatType2>
void
copy(const FloatType1& orig, FloatType2& dest)
{
	dest = orig;
}

template<typename FloatType1, typename FloatType2>
void
copy(const std::complex<FloatType1>& orig, FloatType2 &dest)
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

template<typename Iterator, typename FloatType>
void
copyRealImagToComplexSequence(Iterator re, Iterator reEnd, Iterator im, std::vector<std::complex<FloatType> >& cpx)
{
	typename std::vector<std::complex<FloatType> >::iterator cpxIter = cpx.begin();
	while (re != reEnd) {
		*cpxIter++ = std::complex<FloatType>(*re++, *im++);
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

template<typename InputIterator, typename OutputIterator>
void
copyXZSequence(InputIterator orig, InputIterator origEnd, OutputIterator dest)
{
	for ( ; orig != origEnd; ++orig, ++dest) {
		dest->x = orig->x;
		dest->z = orig->z;
	}
}

} // namespace Value
} // namespace Lab

#endif // VALUE_H_
