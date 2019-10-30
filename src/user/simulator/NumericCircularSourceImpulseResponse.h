/***************************************************************************
 *  Copyright 2019 Marcelo Y. Matuda                                       *
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
#ifndef NUMERICCIRCULARSOURCEIMPULSERESPONSE_H
#define NUMERICCIRCULARSOURCEIMPULSERESPONSE_H

#include <cmath>
#include <limits>
#include <random>
#include <vector>

#include "Log.h"
#include "Util.h"



namespace Lab {

template<typename FloatType>
class NumericCircularSourceImpulseResponse {
public:
	NumericCircularSourceImpulseResponse(
					FloatType samplingFreq,
					FloatType propagationSpeed,
					FloatType sourceRadius,
					unsigned int numSubElem);
	~NumericCircularSourceImpulseResponse() {}

	// Return h/c.
	void getImpulseResponse(FloatType x, FloatType y, FloatType z,
				std::size_t& hOffset /* samples */, std::vector<FloatType>& h);
private:
	struct SubElem {
		SubElem(FloatType x, FloatType y) : x(x), y(y), n0(), value() {}
		FloatType x;
		FloatType y;
		std::size_t n0;
		FloatType value;
	};

	FloatType samplingFreq_;
	FloatType propagationSpeed_;
	FloatType subElemArea_;
	std::vector<SubElem> subElem_;
};



template<typename FloatType>
NumericCircularSourceImpulseResponse<FloatType>::NumericCircularSourceImpulseResponse(
		FloatType samplingFreq,
		FloatType propagationSpeed,
		FloatType sourceRadius,
		unsigned int numSubElem)
			: samplingFreq_(samplingFreq)
			, propagationSpeed_(propagationSpeed)
			, subElemArea_(pi * (sourceRadius * sourceRadius) / numSubElem)
{
	std::mt19937 rndGen;
	std::random_device rd;
	rndGen.seed(rd());
	std::uniform_real_distribution<FloatType> xDist(-1.0, 1.0);
	std::uniform_real_distribution<FloatType> yDist(0.0, 1.0);

	while (subElem_.size() < numSubElem) {
		const FloatType x = xDist(rndGen) * sourceRadius;
		const FloatType y = yDist(rndGen) * sourceRadius;
		if (std::sqrt(x * x + y * y) < sourceRadius) {
			subElem_.emplace_back(x, y);
		}
	}

	LOG_DEBUG << "[NumericCircularSourceImpulseResponse] numSubElem=" << numSubElem;
	//TODO: Export sub elem positions.
}

template<typename FloatType>
void
NumericCircularSourceImpulseResponse<FloatType>::getImpulseResponse(
								FloatType x,
								FloatType y,
								FloatType z,
								std::size_t& hOffset,
								std::vector<FloatType>& h)
{
	// The field is symmetric.
	x = std::sqrt(x * x + y * y);
	y = 0;
	z = std::abs(z);

	std::size_t minN0 = std::numeric_limits<std::size_t>::max();
	std::size_t maxN0 = 0;
	const FloatType k1 = samplingFreq_ / propagationSpeed_;
	const FloatType k2 = samplingFreq_ * subElemArea_ / (FloatType(2.0 * pi) * propagationSpeed_);
	const FloatType z2 = z * z;
	for (unsigned int i = 0, size = subElem_.size(); i < size; ++i) {
		SubElem& se = subElem_[i];
		const FloatType dx = x - se.x;
		const FloatType dy = y - se.y;
		const FloatType r = std::sqrt(dx * dx + dy * dy + z2);
		se.n0 = static_cast<std::size_t>(std::nearbyint(r * k1));
		se.value = k2 / r;
		if (se.n0 < minN0) minN0 = se.n0;
		if (se.n0 > maxN0) maxN0 = se.n0;
	}

	h.assign(maxN0 - minN0 + 1, 0);
	for (auto it = subElem_.begin(); it != subElem_.end(); ++it) {
		h[it->n0 - minN0] += it->value;
	}

	hOffset = minN0;
}

} // namespace Lab

#endif // NUMERICCIRCULARSOURCEIMPULSERESPONSE_H
