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
#ifndef ANALYTICRECTANGULARSOURCEIMPULSERESPONSE_H
#define ANALYTICRECTANGULARSOURCEIMPULSERESPONSE_H

#include <algorithm> /* max, min */
#include <cmath>
#include <cstddef> /* std::size_t */
#include <limits>
#include <vector>
#include <utility> /* swap */

#include "Exception.h"
#include "Util.h"



namespace Lab {

// Calculate the acoustic field generated by a rectangular surface,
// using the analytic solution provided by:
// Emeterio, J. L. S.
// Ullate, L. G.
// Diffraction impulse response of rectangular transducers.
// J. Acoust. Soc. Am., vol. 92, no. 2, pp. 651-662, 1992.
// DOI: 10.1121/1.403990
//
// Note:
// - The source is surrounded by a rigid baffle.
template<typename FloatType>
class AnalyticRectangularSourceImpulseResponse {
public:
	AnalyticRectangularSourceImpulseResponse(
					FloatType samplingFreq,
					FloatType propagationSpeed,
					FloatType sourceWidth,
					FloatType sourceHeight,
					FloatType minEdgeDivisor);
	~AnalyticRectangularSourceImpulseResponse() {}

	// Return h/c.
	void getImpulseResponse(FloatType x, FloatType y, FloatType z,
				std::size_t& hOffset /* samples */, std::vector<FloatType>& h);
private:
	static FloatType alphaP(FloatType v);
	static FloatType alphaN(FloatType v);

	FloatType samplingFreq_;
	FloatType propagationSpeed_;
	FloatType a_;
	FloatType b_;
	FloatType minADivisor_;
	bool swapXY_;
};



template<typename FloatType>
AnalyticRectangularSourceImpulseResponse<FloatType>::AnalyticRectangularSourceImpulseResponse(
		FloatType samplingFreq,
		FloatType propagationSpeed,
		FloatType sourceWidth,
		FloatType sourceHeight,
		FloatType minEdgeDivisor)
			: samplingFreq_(samplingFreq)
			, propagationSpeed_(propagationSpeed)
			, a_(sourceWidth / 2)
			, b_(sourceHeight / 2)
			, minADivisor_(minEdgeDivisor / 2)
			, swapXY_()
{
	if (a_ > b_) {
		std::swap(a_, b_);
		swapXY_ = true;
	}
}

template<typename FloatType>
FloatType
AnalyticRectangularSourceImpulseResponse<FloatType>::alphaP(FloatType v)
{
	return std::asin(std::min(v, FloatType(1)));
}

template<typename FloatType>
FloatType
AnalyticRectangularSourceImpulseResponse<FloatType>::alphaN(FloatType v)
{
	return -std::asin(std::min(-v, FloatType(1)));
}

template<typename FloatType>
void
AnalyticRectangularSourceImpulseResponse<FloatType>::getImpulseResponse(
								FloatType x,
								FloatType y,
								FloatType z,
								std::size_t& hOffset,
								std::vector<FloatType>& h)
{
	// The field is symmetric.
	x = std::abs(x);
	y = std::abs(y);
	z = std::abs(z);

	if (swapXY_) std::swap(x, y);

	const FloatType z2 = z * z;
	const FloatType c2 = propagationSpeed_ * propagationSpeed_;
	constexpr FloatType halfPi = pi / 2.0;
	constexpr FloatType halfPi3 = 3.0 * pi / 2.0;

	// Figure 2.
	const FloatType d1 = x - a_;
	const FloatType d2 = y - b_;
	const FloatType d3 = x + a_; // d3 > 0
	const FloatType d4 = y + b_; // d4 > 0

	const FloatType d1_2 = d1 * d1;
	const FloatType d2_2 = d2 * d2;
	const FloatType d3_2 = d3 * d3;
	const FloatType d4_2 = d4 * d4;

	const FloatType invC = 1 / propagationSpeed_;

	// (7)
	const FloatType ta = std::sqrt(d1_2 + d2_2 + z2) * invC;
	const FloatType tb = std::sqrt(d2_2 + d3_2 + z2) * invC;
	const FloatType tc = std::sqrt(d1_2 + d4_2 + z2) * invC;
	const FloatType td = std::sqrt(d3_2 + d4_2 + z2) * invC;
	// (8)
	const FloatType ts1 = std::sqrt(d1_2 + z2) * invC;
	const FloatType ts2 = std::sqrt(d2_2 + z2) * invC;
	// (9)
	const FloatType t0 = z * invC;

	// Determine the region and the start time of the impulse response.
	unsigned int region;
	FloatType tMin;
	if (x >= a_) {
		if (y >= b_) {
			region = 1;
			tMin = ta;
		} else {
			region = 3;
			tMin = ts1;
		}
	} else {
		if (y >= b_) {
			region = 2;
			tMin = ts2;
		} else {
			region = 4;
			tMin = t0;
		}
	}

	const FloatType dt = 1 / samplingFreq_;
	if (minADivisor_ > 0.0) {
		const FloatType deltaA = a_ / minADivisor_;
		const FloatType sigma1 = std::sqrt(std::max(c2 * tMin * tMin - z2, FloatType(0)));
		const FloatType sigma2 = sigma1 + deltaA;
		const FloatType maxDt = std::sqrt(sigma2 * sigma2 + z2) * invC - tMin;
		if (dt > maxDt) {
			THROW_EXCEPTION(InvalidValueException,
				"The sampling rate is too low (fs=" << samplingFreq_ << " min_fs=" << 1.0 / maxDt << ").");
		}
	}

	// (13)
	const FloatType tm1 = std::min(tb, tc);
	const FloatType tm2 = std::max(tb, tc);

	const std::size_t minAbsoluteIndex = static_cast<std::size_t>(std::ceil(tMin * samplingFreq_));
	const std::size_t maxAbsoluteIndex = static_cast<std::size_t>(std::ceil(td * samplingFreq_));
	h.assign(maxAbsoluteIndex - minAbsoluteIndex + 1, 0);
	const FloatType tOffset = minAbsoluteIndex * dt;

	const FloatType sigmaEps = a_ * std::numeric_limits<FloatType>::epsilon();

	switch (region) {
	case 1:
	{
		const std::size_t i1 = 0;
		const std::size_t i2 = static_cast<std::size_t>(std::ceil(tm1 * samplingFreq_)) - minAbsoluteIndex;
		const std::size_t i3 = static_cast<std::size_t>(std::ceil(tm2 * samplingFreq_)) - minAbsoluteIndex;
		const std::size_t i4 = maxAbsoluteIndex - minAbsoluteIndex;

		for (std::size_t i = i1; i < i2; ++i) {
			const FloatType t = tOffset + i * dt;
			const FloatType sigma = std::sqrt(std::max(c2 * t * t - z2, FloatType(0)));
			if (sigma <= sigmaEps) continue;
			const FloatType invSigma = 1 / sigma;
			const FloatType alpha1 = alphaP(d1 * invSigma);
			const FloatType alpha2 = alphaP(d2 * invSigma);
			h[i] = halfPi - alpha1 - alpha2;
		}
		if (tb <= tc) {
			for (std::size_t i = i2; i < i3; ++i) {
				const FloatType t = tOffset + i * dt;
				const FloatType invSigma = 1 / std::sqrt(c2 * t * t - z2);
				const FloatType alpha1 = alphaP(d1 * invSigma);
				const FloatType alpha3 = alphaP(d3 * invSigma);
				h[i] = alpha3 - alpha1;
			}
		} else {
			for (std::size_t i = i2; i < i3; ++i) {
				const FloatType t = tOffset + i * dt;
				const FloatType invSigma = 1 / std::sqrt(c2 * t * t - z2);
				const FloatType alpha2 = alphaP(d2 * invSigma);
				const FloatType alpha4 = alphaP(d4 * invSigma);
				h[i] = alpha4 - alpha2;
			}
		}
		for (std::size_t i = i3; i < i4; ++i) {
			const FloatType t = tOffset + i * dt;
			const FloatType invSigma = 1 / std::sqrt(c2 * t * t - z2);
			const FloatType alpha3 = alphaP(d3 * invSigma);
			const FloatType alpha4 = alphaP(d4 * invSigma);
			h[i] = alpha3 + alpha4 - halfPi;
		}
		break;
	}
	case 2:
	{
		const std::size_t i0 = 0;
		const std::size_t i1 = static_cast<std::size_t>(std::ceil(ta  * samplingFreq_)) - minAbsoluteIndex;
		const std::size_t i2 = static_cast<std::size_t>(std::ceil(tm1 * samplingFreq_)) - minAbsoluteIndex;
		const std::size_t i3 = static_cast<std::size_t>(std::ceil(tm2 * samplingFreq_)) - minAbsoluteIndex;
		const std::size_t i4 = maxAbsoluteIndex - minAbsoluteIndex;

		for (std::size_t i = i0; i < i1; ++i) {
			const FloatType t = tOffset + i * dt;
			const FloatType sigma = std::sqrt(std::max(c2 * t * t - z2, FloatType(0)));
			if (sigma <= sigmaEps) continue;
			const FloatType alpha2 = alphaP(d2 / sigma);
			h[i] = FloatType(pi) - 2 * alpha2;
		}
		for (std::size_t i = i1; i < i2; ++i) {
			const FloatType t = tOffset + i * dt;
			const FloatType sigma = std::sqrt(std::max(c2 * t * t - z2, FloatType(0)));
			if (sigma <= sigmaEps) continue;
			const FloatType invSigma = 1 / sigma;
			const FloatType alpha1 = alphaN(d1 * invSigma);
			const FloatType alpha2 = alphaP(d2 * invSigma);
			h[i] = halfPi - alpha1 - alpha2;
		}
		for (std::size_t i = i2; i < i3; ++i) {
			const FloatType t = tOffset + i * dt;
			const FloatType invSigma = 1 / std::sqrt(c2 * t * t - z2);
			const FloatType alpha1 = alphaN(d1 * invSigma);
			const FloatType alpha3 = alphaP(d3 * invSigma);
			const FloatType alpha4 = alphaP(d4 * invSigma);
			h[i] = alpha3 - alpha1 - FloatType(pi) + 2 * alpha4;
		}
		for (std::size_t i = i3; i < i4; ++i) {
			const FloatType t = tOffset + i * dt;
			const FloatType invSigma = 1 / std::sqrt(c2 * t * t - z2);
			const FloatType alpha3 = alphaP(d3 * invSigma);
			const FloatType alpha4 = alphaP(d4 * invSigma);
			h[i] = alpha3 + alpha4 - halfPi;
		}
		break;
	}
	case 3:
	{
		const std::size_t i0 = 0;
		const std::size_t i1 = static_cast<std::size_t>(std::ceil(ta  * samplingFreq_)) - minAbsoluteIndex;
		const std::size_t i2 = static_cast<std::size_t>(std::ceil(tm1 * samplingFreq_)) - minAbsoluteIndex;
		const std::size_t i3 = static_cast<std::size_t>(std::ceil(tm2 * samplingFreq_)) - minAbsoluteIndex;
		const std::size_t i4 = maxAbsoluteIndex - minAbsoluteIndex;

		for (std::size_t i = i0; i < i1; ++i) {
			const FloatType t = tOffset + i * dt;
			const FloatType sigma = std::sqrt(std::max(c2 * t * t - z2, FloatType(0)));
			if (sigma <= sigmaEps) continue;
			const FloatType invSigma = 1 / sigma;
			const FloatType alpha1 = alphaP(d1 * invSigma);
			const FloatType alpha3 = alphaP(d3 * invSigma);
			h[i] = 2 * (alpha3 - alpha1);
		}
		for (std::size_t i = i1; i < i2; ++i) {
			const FloatType t = tOffset + i * dt;
			const FloatType sigma = std::sqrt(std::max(c2 * t * t - z2, FloatType(0)));
			if (sigma <= sigmaEps) continue;
			const FloatType invSigma = 1 / sigma;
			const FloatType alpha1 = alphaP(d1 * invSigma);
			const FloatType alpha2 = alphaN(d2 * invSigma);
			const FloatType alpha3 = alphaP(d3 * invSigma);
			h[i] = 2 * alpha3 - halfPi - alpha1 - alpha2;
		}
		if (tb <= tc ) {
			for (std::size_t i = i2; i < i3; ++i) {
				const FloatType t = tOffset + i * dt;
				const FloatType invSigma = 1 / std::sqrt(c2 * t * t - z2);
				const FloatType alpha1 = alphaP(d1 * invSigma);
				const FloatType alpha3 = alphaP(d3 * invSigma);
				h[i] = alpha3 - alpha1;
			}
		} else {
			for (std::size_t i = i2; i < i3; ++i) {
				const FloatType t = tOffset + i * dt;
				const FloatType invSigma = 1 / std::sqrt(c2 * t * t - z2);
				const FloatType alpha2 = alphaN(d2 * invSigma);
				const FloatType alpha3 = alphaP(d3 * invSigma);
				const FloatType alpha4 = alphaP(d4 * invSigma);
				h[i] = alpha4 - alpha2 - FloatType(pi) + 2 * alpha3;
			}
		}
		for (std::size_t i = i3; i < i4; ++i) {
			const FloatType t = tOffset + i * dt;
			const FloatType invSigma = 1 / std::sqrt(c2 * t * t - z2);
			const FloatType alpha3 = alphaP(d3 * invSigma);
			const FloatType alpha4 = alphaP(d4 * invSigma);
			h[i] = alpha3 + alpha4 - halfPi;
		}
		break;
	}
	case 4:
	{
		const std::size_t i0 = 0;
		const std::size_t i1 = static_cast<std::size_t>(std::ceil(ta  * samplingFreq_)) - minAbsoluteIndex;
		const std::size_t i2 = static_cast<std::size_t>(std::ceil(tm1 * samplingFreq_)) - minAbsoluteIndex;
		const std::size_t i3 = static_cast<std::size_t>(std::ceil(tm2 * samplingFreq_)) - minAbsoluteIndex;
		const std::size_t i4 = maxAbsoluteIndex - minAbsoluteIndex;

		for (std::size_t i = i0; i < i1; ++i) {
			const FloatType t = tOffset + i * dt;
			const FloatType sigma = std::sqrt(std::max(c2 * t * t - z2, FloatType(0)));
			if (sigma <= sigmaEps) continue;
			const FloatType invSigma = 1 / sigma;
			const FloatType alpha1 = alphaN(d1 * invSigma);
			const FloatType alpha2 = alphaN(d2 * invSigma);
			const FloatType alpha3 = alphaP(d3 * invSigma);
			const FloatType alpha4 = alphaP(d4 * invSigma);
			h[i] = 2 * (alpha3 + alpha4 - alpha1 - alpha2 - FloatType(pi));
		}
		for (std::size_t i = i1; i < i2; ++i) {
			const FloatType t = tOffset + i * dt;
			const FloatType sigma = std::sqrt(std::max(c2 * t * t - z2, FloatType(0)));
			if (sigma <= sigmaEps) continue;
			const FloatType invSigma = 1 / sigma;
			const FloatType alpha1 = alphaN(d1 * invSigma);
			const FloatType alpha2 = alphaN(d2 * invSigma);
			const FloatType alpha3 = alphaP(d3 * invSigma);
			const FloatType alpha4 = alphaP(d4 * invSigma);
			h[i] = 2 * (alpha3 + alpha4) - alpha1 - alpha2 - halfPi3;
		}
		if (tb <= tc) {
			for (std::size_t i = i2; i < i3; ++i) {
				const FloatType t = tOffset + i * dt;
				const FloatType invSigma = 1 / std::sqrt(c2 * t * t - z2);
				const FloatType alpha1 = alphaN(d1 * invSigma);
				const FloatType alpha3 = alphaP(d3 * invSigma);
				const FloatType alpha4 = alphaP(d4 * invSigma);
				h[i] = alpha3 - alpha1 + 2 * alpha4 - FloatType(pi);
			}
		} else {
			for (std::size_t i = i2; i < i3; ++i) {
				const FloatType t = tOffset + i * dt;
				const FloatType invSigma = 1 / std::sqrt(c2 * t * t - z2);
				const FloatType alpha2 = alphaN(d2 * invSigma);
				const FloatType alpha3 = alphaP(d3 * invSigma);
				const FloatType alpha4 = alphaP(d4 * invSigma);
				h[i] = alpha4 - alpha2 + 2 * alpha3 - FloatType(pi);
			}
		}
		for (std::size_t i = i3; i < i4; ++i) {
			const FloatType t = tOffset + i * dt;
			const FloatType invSigma = 1 / std::sqrt(c2 * t * t - z2);
			const FloatType alpha3 = alphaP(d3 * invSigma);
			const FloatType alpha4 = alphaP(d4 * invSigma);
			h[i] = alpha3 + alpha4 - halfPi;
		}
		break;
	}
	}

	//Util::multiply(h, propagationSpeed_ / FloatType(2.0 * pi));
	Util::multiply(h, 1 / FloatType(2.0 * pi));
	hOffset = minAbsoluteIndex;
}

} // namespace Lab

#endif // ANALYTICRECTANGULARSOURCEIMPULSERESPONSE_H
