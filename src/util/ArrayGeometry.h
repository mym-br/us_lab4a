#ifndef ARRAYGEOMETRY_H
#define ARRAYGEOMETRY_H

#include <vector>

#include "Exception.h"



namespace Lab {
namespace ArrayGeometry {

template<typename FloatType> void getElementXCentered2D(unsigned int numElements, FloatType pitch, std::vector<FloatType>& x);
template<typename FloatType> void getElementX2D(unsigned int numElements, FloatType pitch, FloatType offset, std::vector<FloatType>& x);
template<typename FloatType, typename VectorType> void getElementX2D(unsigned int numElements, FloatType pitch, FloatType offset, VectorType& x);



template<typename FloatType>
void
getElementXCentered2D(unsigned int numElements, FloatType pitch, std::vector<FloatType>& x)
{
	if (numElements == 0) {
		THROW_EXCEPTION(InvalidParameterException, "The number of elements is zero.");
	}
	x.resize(numElements);
	const FloatType xCenter = ((numElements - 1) * pitch) / 2;
	for (unsigned int elem = 0; elem < numElements; ++elem) {
		x[elem] = elem * pitch - xCenter;
	}
}

template<typename FloatType>
void
getElementX2D(unsigned int numElements, FloatType pitch, FloatType offset, std::vector<FloatType>& x)
{
	if (numElements == 0) {
		THROW_EXCEPTION(InvalidParameterException, "The number of elements is zero.");
	}
	x.resize(numElements);
	const FloatType xCenter = ((numElements - 1) * pitch) / 2 - offset;
	for (unsigned int elem = 0; elem < numElements; ++elem) {
		x[elem] = elem * pitch - xCenter;
	}
}

template<typename FloatType, typename VectorType>
void
getElementX2D(unsigned int numElements, FloatType pitch, FloatType offset, VectorType& x)
{
	if (numElements == 0) {
		THROW_EXCEPTION(InvalidParameterException, "The number of elements is zero.");
	}
	x.resize(numElements);
	const FloatType xCenter = ((numElements - 1) * pitch) / 2 - offset;
	for (unsigned int elem = 0; elem < numElements; ++elem) {
		x[elem] = elem * pitch - xCenter;
	}
}

} // namespace ArrayGeometry
} // namespace Lab

#endif // ARRAYGEOMETRY_H
