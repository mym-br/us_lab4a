#ifndef FITTING_H_
#define FITTING_H_

#include <cmath> /* abs, isnan, sqrt */
#include <cstddef> /* std::size_t */
#include <iostream>
#include <utility> /* pair */
#include <vector>

#include "Exception.h"



namespace Lab {
namespace Fitting {

template<typename TFloat> void circleFittingByPratt(const std::vector<std::pair<TFloat, TFloat>>& pointList /* x, y */,
							TFloat& xCenter, TFloat& yCenter, TFloat& radius, unsigned int* numPointsWithNaN=nullptr);



//-----------------------------------------------------------------------------
// Function converted from CircleFitByPratt.m.
//
// V. Pratt, "Direct least-squares fitting of algebraic surfaces",
// Computer Graphics, Vol. 21, pages 145-152 (1987)
//
// 2011-05-15 Converted by Marcelo Y. Matuda.
//-----------------------------------------------------------------------------
// Copyright (c) 2009, Nikolai Chernov
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//    * Redistributions of source code must retain the above copyright
//      notice, this list of conditions and the following disclaimer.
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in
//      the documentation and/or other materials provided with the distribution
//    * Neither the name of the University of Alabama at Birmingham nor the names
//      of its contributors may be used to endorse or promote products derived
//      from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
template<typename TFloat>
void
circleFittingByPratt(const std::vector<std::pair<TFloat, TFloat>>& pointList /* x, y */, TFloat& xCenter, TFloat& yCenter, TFloat& radius, unsigned int* numPointsWithNaN)
{
	if (pointList.size() < 3) {
		THROW_EXCEPTION(InvalidParameterException, "The number of points (" << pointList.size() << ") is not >= 3.");
	}

	double sumX = 0.0, sumY = 0.0;
	std::size_t n = 0; // number of valid data points
	for (typename std::vector<std::pair<TFloat, TFloat>>::const_iterator iter = pointList.begin(); iter != pointList.end(); ++iter) {
		if (std::isnan(iter->first) || std::isnan(iter->second)) continue; // ignore points with NaN

		++n;
		sumX += iter->first;
		sumY += iter->second;
	}
	if (numPointsWithNaN != nullptr) {
		if (n != pointList.size()) {
			std::cerr << "[Fitting::circleFittingByPratt] #################### Number of points with NaN: " << (pointList.size() - n) << std::endl;
			*numPointsWithNaN = pointList.size() - n;
		} else {
			*numPointsWithNaN = 0;
		}
	}

	const double centroidX = sumX / n; // the centroid of the data set
	const double centroidY = sumY / n;

	// Computing moments (note: all moments will be normed, i.e. divided by n).
	double mxx = 0.0, myy = 0.0, mxy = 0.0, mxz = 0.0, myz = 0.0, mzz = 0.0;
	for (std::size_t i = 0, size = pointList.size(); i < size; ++i) {
		if (pointList[i].first != pointList[i].first || pointList[i].second != pointList[i].second) continue; // ignore points with NaN

		const double xi = pointList[i].first - centroidX; // centering data
		const double yi = pointList[i].second - centroidY; // centering data
		const double zi = xi * xi + yi * yi;
		mxx += xi * xi;
		myy += yi * yi;
		mxy += xi * yi;
		mxz += xi * zi;
		myz += yi * zi;
		mzz += zi * zi;
	}
	mxx /= n;
	myy /= n;
	mxy /= n;
	mxz /= n;
	myz /= n;
	mzz /= n;

	// Computing the coefficients of the characteristic polynomial.
	const double mz = mxx + myy;
	const double cov_xy = mxx * myy - mxy * mxy;
	const double mxz2 = mxz * mxz;
	const double myz2 = myz * myz;
	const double mz2 = mz * mz;
	const double a2 = 4 * cov_xy - 3 * mz2 - mzz;
	const double a1 = (mzz + 4 * cov_xy - mz2) * mz - mxz2 - myz2 ;
	const double a0 = mxz2 * myy + myz2 * mxx + (mz2 - mzz) * cov_xy - 2 * mxz * myz * mxy;
	const double a22 = a2 + a2;

	const double epsilon = 1.0e-12;
	double xOld, xNew = 0.0;
	double yOld, yNew = 1.0e20;
	const unsigned int iterMax = 20;

	// Newton's method starting at x=0.
	for (unsigned int iter = 0; iter < iterMax; ++iter) {
		yOld = yNew;
		yNew = a0 + xNew * (a1 + xNew * (a2 + 4.0 * xNew * xNew));
		if (std::abs(yNew) > std::abs(yOld)) {
			THROW_EXCEPTION(InvalidStateException, "Newton-Pratt goes in the wrong direction: |yNew| > |yOld|.");
		}
		const double dy = a1 + xNew * (a22 + 16.0 * xNew * xNew);
		xOld = xNew;
		xNew = xOld - yNew / dy;
		if (std::abs((xNew - xOld) / xNew) < epsilon) {
			break;
		}
		if (iter == iterMax - 1) {
			THROW_EXCEPTION(InvalidStateException, "Newton-Pratt will not converge.");
		}
		if (xNew < 0.0) {
			THROW_EXCEPTION(InvalidStateException, "Newton-Pratt: negative root: " << xNew << '.');
		}
	}

	// Computing the circle parameters.
	const double det = xNew * (xNew - mz) + cov_xy;
	const double xCenterRel = (mxz * (myy - xNew) - myz * mxy) / det / 2.0;
	const double yCenterRel = (myz * (mxx - xNew) - mxz * mxy) / det / 2.0;
	xCenter = xCenterRel + centroidX;
	yCenter = yCenterRel + centroidY;
	radius = std::sqrt(xCenterRel * xCenterRel + yCenterRel * yCenterRel + mz + 2.0 * xNew);
}

} // namespace Fitting
} // namespace Lab

#endif /* FITTING_H_ */
