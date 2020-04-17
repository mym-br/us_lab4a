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

#include <cmath> /* abs, floor, nearbyint, sqrt */
#include <cstddef> /* std::size_t */
#include <limits>
#include <vector>

#include "Log.h"
#include "Util.h" /* pi */

#define NUMERIC_CIRCULAR_SOURCE_IMPULSE_RESPONSE_USE_RANDOM 1

#ifdef NUMERIC_CIRCULAR_SOURCE_IMPULSE_RESPONSE_USE_RANDOM
# include <random>
#endif

namespace Lab {

template<typename TFloat>
class NumericCircularSourceImpulseResponse {
public:
	NumericCircularSourceImpulseResponse(
					TFloat samplingFreq,
					TFloat propagationSpeed,
					TFloat sourceRadius,
					TFloat numSubElemInRadius);

	// Return h/c.
	void getImpulseResponse(TFloat x, TFloat y, TFloat z,
				std::size_t& hOffset /* samples */, std::vector<TFloat>& h);
private:
	struct SubElem {
		SubElem(TFloat x, TFloat y) : x(x), y(y), n0(), value() {}
		TFloat x;
		TFloat y;
		std::size_t n0;
		TFloat value;
	};

	TFloat samplingFreq_;
	TFloat propagationSpeed_;
	TFloat subElemArea_;
	std::vector<SubElem> subElem_;
};



template<typename TFloat>
NumericCircularSourceImpulseResponse<TFloat>::NumericCircularSourceImpulseResponse(
		TFloat samplingFreq,
		TFloat propagationSpeed,
		TFloat sourceRadius,
		TFloat numSubElemInRadius)
			: samplingFreq_(samplingFreq)
			, propagationSpeed_(propagationSpeed)
			, subElemArea_()
{
#ifdef NUMERIC_CIRCULAR_SOURCE_IMPULSE_RESPONSE_USE_RANDOM
	const TFloat area = pi * (sourceRadius * sourceRadius);
	const TFloat subElemDensity = numSubElemInRadius * numSubElemInRadius / (sourceRadius * sourceRadius);
	const unsigned int numSubElem = static_cast<unsigned int>(subElemDensity * area);
	subElemArea_ = area / numSubElem;

	std::mt19937 rndGen;
	std::random_device rd;
	rndGen.seed(rd());
	std::uniform_real_distribution<TFloat> dist(0.0, 1.0);

	// http://extremelearning.com.au/a-simple-method-to-construct-isotropic-quasirandom-blue-noise-point-sequences/
	// Infinite R2 jittered sequence.
	// i = 0, 1, 2, ...
	// u0, u1 = random numbers [0.0, 1.0)
	auto jitteredPoint2D = [&](int i, double u0, double u1, TFloat& x, TFloat& y) {
		constexpr double lambda = 0.5; // jitter parameter ( > 0)
		constexpr double phi = 1.324717957244746;
		constexpr double alpha0 = 1.0 / phi;
		constexpr double alpha1 = 1.0 / (phi * phi);
		constexpr double delta0 = 0.76;
		constexpr double i0 = 0.300;
		const double k = lambda * delta0 * std::sqrt(pi) / (4.0 * std::sqrt(i + i0));
		const double i1 = i + 1;
		const double x0 = alpha0 * i1 + k * u0; // x0 > 0
		const double y0 = alpha1 * i1 + k * u1; // y0 > 0
		x = x0 - std::floor(x0);
		y = y0 - std::floor(y0);
	};

	int i = 0;
	while (subElem_.size() < numSubElem) {
		TFloat x, y;
		jitteredPoint2D(i, dist(rndGen), dist(rndGen), x, y);
		// [0.0, 1.0) --> [-1.0, 1.0)
		x = -1.0 + 2.0 * x;
		y = -1.0 + 2.0 * y;
		x *= sourceRadius;
		y *= sourceRadius;

		if (std::sqrt(x * x + y * y) < sourceRadius) {
			subElem_.emplace_back(x, y);
		}
		++i;
	}
#else
	const TFloat d = sourceRadius / numSubElemInRadius; // sub-element side
	subElemArea_ = d * d * 2; // multiplied by 2 because only one half of the circle is used
	const TFloat d2 = d * 0.5;

	const unsigned int n = numSubElemInRadius;
	for (unsigned int iy = 0; iy < n; ++iy) {
		const TFloat yc = d2 + iy * d;
		const TFloat yt = yc + d2; // to test if the sub-element is inside the circle
		for (unsigned int ix = 0; ix < n; ++ix) {
			const TFloat xc = d2 + ix * d;
			const TFloat xt = xc + d2;
			if (std::sqrt(xt * xt + yt * yt) <= sourceRadius) {
				subElem_.emplace_back(xc, yc);
				subElem_.emplace_back(-xc, yc);
			} else {
				break;
			}
		}
	}
#endif

	LOG_DEBUG << "[NumericCircularSourceImpulseResponse] subElem_.size()=" << subElem_.size();
}

template<typename TFloat>
void
NumericCircularSourceImpulseResponse<TFloat>::getImpulseResponse(
								TFloat x,
								TFloat y,
								TFloat z,
								std::size_t& hOffset,
								std::vector<TFloat>& h)
{
	// The field is symmetric.
	x = std::sqrt(x * x + y * y);
	y = 0;
	z = std::abs(z);

	std::size_t minN0 = std::numeric_limits<std::size_t>::max();
	std::size_t maxN0 = 0;
	const TFloat k1 = samplingFreq_ / propagationSpeed_;
	const TFloat k2 = samplingFreq_ * subElemArea_ / (TFloat(2.0 * pi) * propagationSpeed_);
	const TFloat z2 = z * z;
	for (unsigned int i = 0, size = subElem_.size(); i < size; ++i) {
		SubElem& se = subElem_[i];
		const TFloat dx = x - se.x;
		const TFloat dy = y - se.y;
		const TFloat r = std::sqrt(dx * dx + dy * dy + z2);
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
