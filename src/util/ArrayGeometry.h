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

#ifndef ARRAYGEOMETRY_H
#define ARRAYGEOMETRY_H

#include <vector>

#include "Exception.h"



namespace Lab {
namespace ArrayGeometry {

template<typename TFloat> void getElementXCentered2D(unsigned int numElements, TFloat pitch, std::vector<TFloat>& x);
template<typename TFloat> void getElementX2D(unsigned int numElements, TFloat pitch, TFloat offset, std::vector<TFloat>& x);
template<typename TFloat, typename VectorType> void getElementX2D(unsigned int numElements, TFloat pitch, TFloat offset, VectorType& x);



template<typename TFloat>
void
getElementXCentered2D(unsigned int numElements, TFloat pitch, std::vector<TFloat>& x)
{
	if (numElements == 0) {
		THROW_EXCEPTION(InvalidParameterException, "The number of elements is zero.");
	}
	x.resize(numElements);
	const TFloat xCenter = ((numElements - 1) * pitch) / 2;
	for (unsigned int elem = 0; elem < numElements; ++elem) {
		x[elem] = elem * pitch - xCenter;
	}
}

template<typename TFloat>
void
getElementX2D(unsigned int numElements, TFloat pitch, TFloat offset, std::vector<TFloat>& x)
{
	if (numElements == 0) {
		THROW_EXCEPTION(InvalidParameterException, "The number of elements is zero.");
	}
	x.resize(numElements);
	const TFloat xCenter = ((numElements - 1) * pitch) / 2 - offset;
	for (unsigned int elem = 0; elem < numElements; ++elem) {
		x[elem] = elem * pitch - xCenter;
	}
}

template<typename TFloat, typename VectorType>
void
getElementX2D(unsigned int numElements, TFloat pitch, TFloat offset, VectorType& x)
{
	if (numElements == 0) {
		THROW_EXCEPTION(InvalidParameterException, "The number of elements is zero.");
	}
	x.resize(numElements);
	const TFloat xCenter = ((numElements - 1) * pitch) / 2 - offset;
	for (unsigned int elem = 0; elem < numElements; ++elem) {
		x[elem] = elem * pitch - xCenter;
	}
}

} // namespace ArrayGeometry
} // namespace Lab

#endif // ARRAYGEOMETRY_H
