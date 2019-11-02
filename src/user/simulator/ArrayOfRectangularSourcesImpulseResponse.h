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
#ifndef ARRAYOFRECTANGULARSOURCESIMPULSERESPONSE_H
#define ARRAYOFRECTANGULARSOURCESIMPULSERESPONSE_H

#include <limits>
#include <vector>

#include "Exception.h"
#include "Util.h"
#include "XY.h"



namespace Lab {

template<typename FloatType, typename ImpulseResponse>
class ArrayOfRectangularSourcesImpulseResponse {
public:
	ArrayOfRectangularSourcesImpulseResponse(
					FloatType samplingFreq,
					FloatType propagationSpeed,
					FloatType sourceWidth,
					FloatType sourceHeight,
					FloatType discretization,
					const std::vector<XY<FloatType>>& elemPos,
					const std::vector<FloatType>& focusDelay);
	~ArrayOfRectangularSourcesImpulseResponse() {}

	void getImpulseResponse(FloatType x, FloatType y, FloatType z, std::size_t& hOffset, std::vector<FloatType>& h,
				std::vector<unsigned int>* activeElemList=nullptr);
private:
	const FloatType samplingFreq_;
	ImpulseResponse ir_;
	const std::vector<XY<FloatType>>& elemPos_;
	const std::vector<FloatType>& focusDelay_;
	std::vector<std::size_t> offsetList_;
	std::vector<std::vector<FloatType>> hList_;
};



template<typename FloatType, typename ImpulseResponse>
ArrayOfRectangularSourcesImpulseResponse<FloatType, ImpulseResponse>::ArrayOfRectangularSourcesImpulseResponse(
				FloatType samplingFreq,
				FloatType propagationSpeed,
				FloatType sourceWidth,
				FloatType sourceHeight,
				FloatType discretization,
				const std::vector<XY<FloatType>>& elemPos,
				const std::vector<FloatType>& focusDelay /* s */)
		: samplingFreq_(samplingFreq)
		, ir_(samplingFreq, propagationSpeed, sourceWidth, sourceHeight, discretization)
		, elemPos_(elemPos)
		, focusDelay_(focusDelay)
		, offsetList_(elemPos_.size())
		, hList_(elemPos_.size())
{
	if (elemPos_.empty()) {
		THROW_EXCEPTION(InvalidParameterException, "Empty element position list.");
	}
	if (focusDelay_.empty()) {
		THROW_EXCEPTION(InvalidParameterException, "Empty element delay list.");
	}
}

template<typename FloatType, typename ImpulseResponse>
void
ArrayOfRectangularSourcesImpulseResponse<FloatType, ImpulseResponse>::getImpulseResponse(
		FloatType x, FloatType y, FloatType z, std::size_t& hOffset, std::vector<FloatType>& h,
		std::vector<unsigned int>* activeElemList)
{
	if (activeElemList) {
		if (activeElemList->empty()) {
			THROW_EXCEPTION(InvalidParameterException, "Empty active element list.");
		} else {
			if (activeElemList->size() > elemPos_.size()) {
				THROW_EXCEPTION(InvalidParameterException, "Active element list size > number of elements.");
			}
			offsetList_.resize(activeElemList->size());
			hList_.resize(activeElemList->size());
		}
	} else {
		offsetList_.resize(elemPos_.size());
		hList_.resize(elemPos_.size());
	}
	if (focusDelay_.size() != offsetList_.size()) {
		THROW_EXCEPTION(InvalidParameterException, "Invalid focus delay list size: " << focusDelay_.size()
				<< " (should be " << offsetList_.size() << ").");
	}

	std::size_t iMin = std::numeric_limits<std::size_t>::max();
	std::size_t iMax = 0; // the index after the last

	// Obtain the impulse responses.
	for (std::size_t i = 0, iEnd = elemPos_.size(), j = 0; i < iEnd; ++i) {
		if (activeElemList) {
			if (j >= activeElemList->size()) break;
			if ((*activeElemList)[j] != i) {
				continue;
			}
		}
		const XY<FloatType>& pos = elemPos_[i];

		ir_.getImpulseResponse(x - pos.x, y - pos.y, z, offsetList_[j], hList_[j]);

		offsetList_[j] += static_cast<std::size_t>(std::nearbyint(focusDelay_[j] * samplingFreq_));
		const std::size_t start = offsetList_[j];
		if (start < iMin) iMin = start;
		const std::size_t end = offsetList_[j] + hList_[j].size();
		if (end > iMax) iMax = end;
		++j;
	}

	// Accumulate the impulse responses.
	h.assign(iMax - iMin, 0.0);
	for (std::size_t i = 0, iEnd = elemPos_.size(), j = 0; i < iEnd; ++i) {
		if (activeElemList) {
			if (j >= activeElemList->size()) break;
			if ((*activeElemList)[j] != i) {
				continue;
			}
		}
		const std::size_t hBegin = offsetList_[j];
		const std::size_t hEnd = hBegin + hList_[j].size();
		Util::addElements(hList_[j].begin(), h.begin() + hBegin - iMin, h.begin() + hEnd - iMin);
		++j;
	}
	hOffset = iMin;
}

} // namespace Lab

#endif // ARRAYOFRECTANGULARSOURCESIMPULSERESPONSE_H
