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
#ifndef NUMERIC_ARRAY_OF_RECTANGULAR_SOURCES_IMPULSE_RESPONSE_H
#define NUMERIC_ARRAY_OF_RECTANGULAR_SOURCES_IMPULSE_RESPONSE_H

#include <cmath> /* ceil, nearbyint, sqrt */
#include <cstddef> /* std::size_t */
#include <limits>
#include <vector>

#include "Exception.h"
#include "Log.h"
#include "Matrix.h"
#include "Tensor3.h"
#include "Util.h" /* pi */
#include "XY.h"



namespace Lab {

// Calculate the acoustic field generated by a flat rectangular surface,
// using the numeric solution provided by:
//
// Piwakowski, B.
// Delannoy, B.
// Method for computing spatial pulse response: Time-domain approach
// J. Acoust. Soc. Am., vol. 86, no. 6, pp. 2422-2432, 1989.
// DOI: 10.1121/1.398449
//
// See also:
// Lasota, H.
// Salamon, R.
// Delannoy, B.
// Acoustic diffraction analysis by the impulse response method: A line impulse response approach},
// J. Acoust. Soc. Am., vol. 76, no. 1, pp. 280-290, 1984.
// DOI: 10.1121/1.391115
//
// Note:
// - The source is surrounded by a rigid baffle.
template<typename TFloat>
class NumericArrayOfRectangularSourcesImpulseResponse {
public:
	NumericArrayOfRectangularSourcesImpulseResponse(
					TFloat samplingFreq,
					TFloat propagationSpeed,
					TFloat sourceWidth,
					TFloat sourceHeight,
					TFloat subElemSize,
					const std::vector<XY<TFloat>>& elemPos,
					const std::vector<TFloat>& focusDelay /* s */);

	// Return h/c.
	void getImpulseResponse(TFloat x, TFloat y, TFloat z,
				std::size_t& hOffset /* samples */, std::vector<TFloat>& h,
				std::vector<unsigned int>* activeElemList=nullptr);
private:
	struct SubElemData {
		std::size_t n0;
		TFloat value;
	};

	const std::vector<XY<TFloat>>& elemPos_;
	TFloat samplingFreq_;
	TFloat propagationSpeed_;
	TFloat subElemWidth_;
	TFloat subElemHeight_;
	unsigned int numElem_;
	unsigned int numSubElemX_; // per element
	unsigned int numSubElemY_; // per element
	Matrix<XY<TFloat>> subElemPos_;
	Tensor3<SubElemData> subElemData_;
	std::vector<unsigned int> activeElem_;
	std::vector<TFloat> elemDelay_;
};



template<typename TFloat>
NumericArrayOfRectangularSourcesImpulseResponse<TFloat>::NumericArrayOfRectangularSourcesImpulseResponse(
		TFloat samplingFreq,
		TFloat propagationSpeed,
		TFloat sourceWidth,
		TFloat sourceHeight,
		TFloat subElemSize,
		const std::vector<XY<TFloat>>& elemPos,
		const std::vector<TFloat>& focusDelay /* s */)
			: elemPos_(elemPos)
			, samplingFreq_(samplingFreq)
			, propagationSpeed_(propagationSpeed)
			, subElemWidth_()
			, subElemHeight_()
			, numElem_()
			, numSubElemX_()
			, numSubElemY_()
{
	if (elemPos.empty()) {
		THROW_EXCEPTION(InvalidParameterException, "Empty element position list.");
	}
	if (focusDelay.empty()) {
		THROW_EXCEPTION(InvalidParameterException, "Empty element delay list.");
	}
	if (focusDelay.size() > elemPos.size()) {
		THROW_EXCEPTION(InvalidParameterException, "Size of focusDelay is greater than size of elemPos.");
	}
	numElem_ = focusDelay.size(); // may be less than elemPos.size()
	activeElem_.resize(numElem_);

	if (subElemSize > sourceWidth) {
		numSubElemX_ = 1;
	} else {
		numSubElemX_ = static_cast<unsigned int>(std::ceil(sourceWidth / subElemSize));
	}
	if (subElemSize > sourceHeight) {
		numSubElemY_ = 1;
	} else {
		numSubElemY_ = static_cast<unsigned int>(std::ceil(sourceHeight / subElemSize));
	}
	subElemWidth_  = sourceWidth  / numSubElemX_;
	subElemHeight_ = sourceHeight / numSubElemY_;

	const TFloat halfW = 0.5 * (numSubElemX_ - 1);
	const TFloat halfH = 0.5 * (numSubElemY_ - 1);
	subElemPos_.resize(numSubElemY_, numSubElemX_);
	for (unsigned int iy = 0; iy < numSubElemY_; ++iy) {
		for (unsigned int ix = 0; ix < numSubElemX_; ++ix) {
			subElemPos_(iy, ix).x = (ix - halfW) * subElemWidth_;
			subElemPos_(iy, ix).y = (iy - halfH) * subElemHeight_;
		}
	}

	subElemData_.resize(numElem_, numSubElemY_, numSubElemX_);

	elemDelay_.resize(numElem_);
	for (unsigned int i = 0; i < numElem_; ++i) {
		elemDelay_[i] = focusDelay[i] * samplingFreq_;
	}

	LOG_DEBUG << "[NumericArrayOfRectangularSourcesImpulseResponse] numSubElemX=" << numSubElemX_ << " numSubElemY=" << numSubElemY_;
}

template<typename TFloat>
void
NumericArrayOfRectangularSourcesImpulseResponse<TFloat>::getImpulseResponse(
								TFloat x,
								TFloat y,
								TFloat z,
								std::size_t& hOffset,
								std::vector<TFloat>& h,
								std::vector<unsigned int>* activeElemList)
{
	if (activeElemList) {
		if (activeElemList->empty()) {
			THROW_EXCEPTION(InvalidParameterException, "Empty active element list.");
		} else {
			if (activeElemList->size() != numElem_) {
				THROW_EXCEPTION(InvalidParameterException, "Active element size is not equal to number of delays.");
			}
			for (unsigned int i = 0; i < numElem_; ++i) {
				activeElem_[i] = (*activeElemList)[i];
			}
		}
	} else {
		for (unsigned int i = 0; i < numElem_; ++i) {
			activeElem_[i] = i;
		}
	}

	std::size_t minN0 = std::numeric_limits<std::size_t>::max();
	std::size_t maxN0 = 0;
	const TFloat k1 = samplingFreq_ / propagationSpeed_;
	const TFloat k2 = samplingFreq_ * subElemWidth_ * subElemHeight_ / (TFloat(2.0 * pi) * propagationSpeed_);
	const TFloat z2 = z * z;
	for (unsigned int elemIdx = 0; elemIdx < numElem_; ++elemIdx) {
		const unsigned int activeElemIdx = activeElem_[elemIdx];
		for (unsigned int iy = 0; iy < numSubElemY_; ++iy) {
			const TFloat dy = y - subElemPos_(iy, 0).y - elemPos_[activeElemIdx].y;
			const TFloat dy2_z2 = dy * dy + z2;
			for (unsigned int ix = 0; ix < numSubElemX_; ++ix) {
				const TFloat dx = x - subElemPos_(iy, ix).x - elemPos_[activeElemIdx].x;
				const TFloat r = std::sqrt(dx * dx + dy2_z2);
				SubElemData& data = subElemData_(elemIdx, iy, ix);
				data.n0 = static_cast<std::size_t>(std::nearbyint(r * k1 + elemDelay_[elemIdx]));
				data.value = k2 / r;
				if (data.n0 < minN0) minN0 = data.n0;
				if (data.n0 > maxN0) maxN0 = data.n0;
			}
		}
	}

	h.assign(maxN0 - minN0 + 1U, 0);
	for (auto& data : subElemData_) {
		h[data.n0 - minN0] += data.value;
	}

	hOffset = minN0;

	//LOG_DEBUG << "[NumericArrayOfRectangularSourcesImpulseResponse] minN0=" << minN0 << " maxN0=" << maxN0;
}

} // namespace Lab

#endif // NUMERIC_ARRAY_OF_RECTANGULAR_SOURCES_IMPULSE_RESPONSE_H
