/***************************************************************************
 *  Copyright 2018 Marcelo Y. Matuda                                       *
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
#ifndef NUMERICARRAYOFRECTANGULARFLATSOURCESIMPULSERESPONSE_H
#define NUMERICARRAYOFRECTANGULARFLATSOURCESIMPULSERESPONSE_H

#include <limits>
#include <vector>

#include "NumericRectangularFlatSourceImpulseResponse.h"
#include "Util.h"
#include "XY.h"



namespace Lab {

template<typename FloatType>
class NumericArrayOfRectangularFlatSourcesImpulseResponse {
public:
	NumericArrayOfRectangularFlatSourcesImpulseResponse(
					FloatType sourceWidth,
					FloatType sourceHeight,
					FloatType samplingFreq,
					FloatType propagationSpeed,
					FloatType subElemSize,
					const std::vector<XY<FloatType>>& elemPos,
					const std::vector<FloatType>& focusDelay);
	~NumericArrayOfRectangularFlatSourcesImpulseResponse() {}

	void getImpulseResponse(FloatType x, FloatType y, FloatType z, std::size_t& hOffset, std::vector<FloatType>& h);
private:
	NumericRectangularFlatSourceImpulseResponse<FloatType> ir_;
	const std::vector<XY<FloatType>>& elemPos_;
	const std::vector<FloatType>& focusDelay_;
	std::vector<std::size_t> offsetList_;
	std::vector<std::vector<FloatType>> hList_;
};



template<typename FloatType>
NumericArrayOfRectangularFlatSourcesImpulseResponse<FloatType>::NumericArrayOfRectangularFlatSourcesImpulseResponse(
				FloatType sourceWidth,
				FloatType sourceHeight,
				FloatType samplingFreq,
				FloatType propagationSpeed,
				FloatType subElemSize,
				const std::vector<XY<FloatType>>& elemPos,
				const std::vector<FloatType>& focusDelay)
		: ir_{sourceWidth, sourceHeight, samplingFreq, propagationSpeed, subElemSize}
		, elemPos_{elemPos}
		, focusDelay_{focusDelay}
{
}

template<typename FloatType>
void
NumericArrayOfRectangularFlatSourcesImpulseResponse<FloatType>::getImpulseResponse(
		FloatType x, FloatType y, FloatType z, std::size_t& hOffset, std::vector<FloatType>& h)
{
	offsetList_.resize(elemPos_.size());
	hList_.resize(elemPos_.size());

	std::size_t iMin = std::numeric_limits<std::size_t>::max();
	std::size_t iMax = 0; // the index after the last

	// Obtain the impulse responses.
	for (std::size_t i = 0, iEnd = elemPos_.size(); i < iEnd; ++i) {
		const XY<FloatType>& pos = elemPos_[i];

		ir_.getImpulseResponse(x - pos.x, y - pos.y, z, offsetList_[i], hList_[i]);

		offsetList_[i] += static_cast<std::size_t>(std::nearbyint(focusDelay_[i]));
		const std::size_t start = offsetList_[i];
		if (start < iMin) iMin = start;
		const std::size_t end = offsetList_[i] + hList_[i].size();
		if (end > iMax) iMax = end;
	}

	// Accumulate the impulse responses.
	h.resize(iMax - iMin);
	for (std::size_t i = 0, iEnd = elemPos_.size(); i < iEnd; ++i) {
		const std::size_t hBegin = offsetList_[i];
		const std::size_t hEnd = hBegin + hList_[i].size();
		Util::addElements(hList_[i].begin(), h.begin() + hBegin - iMin, h.begin() + hEnd - iMin);
	}
	hOffset = iMin;
}

} // namespace Lab

#endif // NUMERICARRAYOFRECTANGULARFLATSOURCESIMPULSERESPONSE_H
