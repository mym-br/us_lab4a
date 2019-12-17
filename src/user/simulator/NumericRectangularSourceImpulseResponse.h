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
#ifndef NUMERICRECTANGULARSOURCEIMPULSERESPONSE_H
#define NUMERICRECTANGULARSOURCEIMPULSERESPONSE_H

#include <cmath>
#include <limits>
#include <vector>

#include "Log.h"
#include "Matrix.h"
#include "Util.h"



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
class NumericRectangularSourceImpulseResponse {
public:
	NumericRectangularSourceImpulseResponse(
					TFloat samplingFreq,
					TFloat propagationSpeed,
					TFloat sourceWidth,
					TFloat sourceHeight,
					TFloat subElemSize);

	// Return h/c.
	void getImpulseResponse(TFloat x, TFloat y, TFloat z,
				std::size_t& hOffset /* samples */, std::vector<TFloat>& h);
private:
	struct SubElem {
		TFloat x;
		TFloat y;
		std::size_t n0;
		TFloat value;
	};

	TFloat samplingFreq_;
	TFloat propagationSpeed_;
	TFloat subElemWidth_;
	TFloat subElemHeight_;
	unsigned int numElemX_;
	unsigned int numElemY_;
	Matrix<SubElem> subElem_;
};



template<typename TFloat>
NumericRectangularSourceImpulseResponse<TFloat>::NumericRectangularSourceImpulseResponse(
		TFloat samplingFreq,
		TFloat propagationSpeed,
		TFloat sourceWidth,
		TFloat sourceHeight,
		TFloat subElemSize)
			: samplingFreq_(samplingFreq)
			, propagationSpeed_(propagationSpeed)
			, subElemWidth_()
			, subElemHeight_()
			, numElemX_()
			, numElemY_()
{
	if (subElemSize > sourceWidth) {
		numElemX_ = 1;
	} else {
		numElemX_ = static_cast<unsigned int>(std::ceil(sourceWidth / subElemSize));
	}
	if (subElemSize > sourceHeight) {
		numElemY_ = 1;
	} else {
		numElemY_ = static_cast<unsigned int>(std::ceil(sourceHeight / subElemSize));
	}
	subElemWidth_  = sourceWidth  / numElemX_;
	subElemHeight_ = sourceHeight / numElemY_;

	const TFloat halfW = 0.5 * (numElemX_ - 1);
	const TFloat halfH = 0.5 * (numElemY_ - 1);
	subElem_.resize(numElemY_, numElemX_);
	for (unsigned int iy = 0; iy < numElemY_; ++iy) {
		for (unsigned int ix = 0; ix < numElemX_; ++ix) {
			subElem_(iy, ix).x = (ix - halfW) * subElemWidth_;
			subElem_(iy, ix).y = (iy - halfH) * subElemHeight_;
		}
	}

	LOG_DEBUG << "[NumericRectangularSourceImpulseResponse] numElemX=" << numElemX_ << " numElemY=" << numElemY_;
}

template<typename TFloat>
void
NumericRectangularSourceImpulseResponse<TFloat>::getImpulseResponse(
								TFloat x,
								TFloat y,
								TFloat z,
								std::size_t& hOffset,
								std::vector<TFloat>& h)
{
	std::size_t minN0 = std::numeric_limits<std::size_t>::max();
	std::size_t maxN0 = 0;
	const TFloat k1 = samplingFreq_ / propagationSpeed_;
	const TFloat k2 = samplingFreq_ * subElemWidth_ * subElemHeight_ / (TFloat(2.0 * pi) * propagationSpeed_);
	const TFloat z2 = z * z;
	for (unsigned int iy = 0; iy < numElemY_; ++iy) {
		const TFloat dy = y - subElem_(iy, 0).y;
		const TFloat dy2_z2 = dy * dy + z2;
		for (unsigned int ix = 0; ix < numElemX_; ++ix) {
			SubElem& se = subElem_(iy, ix);
			const TFloat dx = x - se.x;
			const TFloat r = std::sqrt(dx * dx + dy2_z2);
			se.n0 = static_cast<std::size_t>(std::nearbyint(r * k1));
			se.value = k2 / r;
			if (se.n0 < minN0) minN0 = se.n0;
			if (se.n0 > maxN0) maxN0 = se.n0;
		}
	}

	h.assign(maxN0 - minN0 + 1, 0);
	for (auto it = subElem_.begin(); it != subElem_.end(); ++it) {
		h[it->n0 - minN0] += it->value;
	}

	hOffset = minN0;
}

} // namespace Lab

#endif // NUMERICRECTANGULARSOURCEIMPULSERESPONSE_H
