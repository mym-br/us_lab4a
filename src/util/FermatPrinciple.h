/***************************************************************************
 *  Copyright 2014, 2017, 2018, 2020 Marcelo Y. Matuda                     *
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

#ifndef FERMATPRINCIPLE_H_
#define FERMATPRINCIPLE_H_

#include <algorithm> /* min, max */
#include <cmath> /* abs, sqrt */
#include <limits>
#include <vector>

#include "SIMD.h"
#include "XZ.h"



namespace Lab {
namespace FermatPrinciple {

// Brute-force.
template<typename TFloat> void findMinTimeForParallelPlanarInterface(
						TFloat c1, TFloat c2,
						const std::vector<TFloat>& xInterfaceList, TFloat zInterface,
						TFloat x1, TFloat z1, TFloat x2, TFloat z2,
						TFloat& tMin, unsigned int& idxMin);

template<typename TFloat> unsigned int calcBlockSizeForTwoStepSearch(
						unsigned int interfaceDataSize, TFloat interfaceStep,
						TFloat wavelength, unsigned int maxBlockSizeInWavelengths);

// In the first step, tests only positions separated by blockSize.
// In the second step, tests positions (blockSize - 1) before and after the position k obtained in the first step.
//
// It is assumed that the minimum will be in that region around the position k.
//
// The number of tests is N/M + 2*M - 2, where:
// N = interfacePointList.size()
// M = blockSize
//
// The optimum M = sqrt(N/2). For this M, the number of tests is 2*sqrt(2*N).
//
template<typename TFloat> void findMinTimeInTwoSteps(
						unsigned int blockSize,
						TFloat c1, TFloat c2,
						const std::vector<XZ<TFloat>>& interfacePointList,
						TFloat x1, TFloat z1, TFloat x2, TFloat z2,
						TFloat& tMin, unsigned int& idxMin);
// Uses SIMD.
template<typename TFloat> void findMinTimeInTwoSteps2(
						unsigned int blockSize,
						TFloat invC1, TFloat invC2,
						const std::vector<XZ<TFloat>>& interfacePointList,
						TFloat x1, TFloat z1, TFloat x2, TFloat z2,
						TFloat& tMin, unsigned int& idxMin);
template<typename TFloat> void findMinTimeInTwoSteps(
						unsigned int blockSize,
						TFloat c1, TFloat c2X, TFloat c2Z,
						const std::vector<XZ<TFloat>>& interfacePointList,
						TFloat x1, TFloat z1, TFloat x2, TFloat z2,
						TFloat& tMin, unsigned int& idxMin);

// Auxiliary.
//
// Executed in two steps like findMinTimeInTwoSteps().
template<typename TFloat> void findNearestPointInXInTwoSteps(
						unsigned int blockSize,
						const std::vector<XZ<TFloat>>& interfacePointList,
						TFloat x,
						TFloat& zIdxMin, unsigned int& idxMin);



template<typename TFloat>
void
findMinTimeForParallelPlanarInterface(TFloat c1, TFloat c2,
		const std::vector<TFloat>& xInterfaceList, TFloat zInterface,
		TFloat x1, TFloat z1, TFloat x2, TFloat z2,
		TFloat& tMin, unsigned int &idxMin)
{
	const TFloat invC1 = 1 / c1;
	const TFloat invC2 = 1 / c2;
	const TFloat dz1 = zInterface - z1;
	const TFloat dz2 = z2 - zInterface;
	const TFloat dz12 = dz1 * dz1;
	const TFloat dz22 = dz2 * dz2;
	tMin = std::numeric_limits<TFloat>::max();
	for (unsigned int i = 0, size = xInterfaceList.size(); i < size; ++i) {
		const TFloat xInterface = xInterfaceList[i];
		const TFloat dx1 = xInterface - x1;
		const TFloat dx2 = x2 - xInterface;
		const TFloat t = std::sqrt(dx1 * dx1 + dz12) * invC1 + std::sqrt(dx2 * dx2 + dz22) * invC2;
		if (t < tMin) {
			tMin = t;
			idxMin = i;
		}
	}
}

template<typename TFloat>
unsigned int
calcBlockSizeForTwoStepSearch(unsigned int interfaceDataSize, TFloat interfaceStep, TFloat wavelength, unsigned int maxBlockSizeInWavelengths)
{
	const TFloat idealBlockSize = std::sqrt(interfaceDataSize / 2.0) + 0.5;
	const TFloat maxBlockSize = (wavelength * maxBlockSizeInWavelengths) / interfaceStep;
	const TFloat blockSize = std::min(idealBlockSize, maxBlockSize);
	return std::max(2U, static_cast<unsigned int>(blockSize));
}

template<typename TFloat>
void
findMinTimeInTwoSteps(
		unsigned int blockSize,
		TFloat c1, TFloat c2,
		const std::vector<XZ<TFloat>>& interfacePointList,
		TFloat x1, TFloat z1, TFloat x2, TFloat z2,
		TFloat& tMin, unsigned int& idxMin)
{
	const TFloat invC1 = 1 / c1;
	const TFloat invC2 = 1 / c2;

	// First step: step = blockSize
	tMin = std::numeric_limits<TFloat>::max();
	idxMin = std::numeric_limits<unsigned int>::max();
	for (unsigned int i = 0, end = interfacePointList.size(); i < end; i += blockSize) {
		const XZ<TFloat>& point = interfacePointList[i];
		const TFloat dx1 = point.x - x1;
		const TFloat dz1 = point.z - z1;
		const TFloat dx2 = x2 - point.x;
		const TFloat dz2 = z2 - point.z;
		const TFloat t = std::sqrt(dx1 * dx1 + dz1 * dz1) * invC1 + std::sqrt(dx2 * dx2 + dz2 * dz2) * invC2;
		if (t < tMin) {
			tMin = t;
			idxMin = i;
		}
	}

	// Second step: step = 1
	const unsigned int iBegin = (idxMin > (blockSize - 1)) ? idxMin - (blockSize - 1) : 0;
	const unsigned int iEnd = std::min<unsigned int>(idxMin + blockSize, interfacePointList.size());
	for (unsigned int i = iBegin; i < iEnd; ++i) {
		const XZ<TFloat>& point = interfacePointList[i];
		const TFloat dx1 = point.x - x1;
		const TFloat dz1 = point.z - z1;
		const TFloat dx2 = x2 - point.x;
		const TFloat dz2 = z2 - point.z;
		const TFloat t = std::sqrt(dx1 * dx1 + dz1 * dz1) * invC1 + std::sqrt(dx2 * dx2 + dz2 * dz2) * invC2;
		if (t < tMin) {
			tMin = t;
			idxMin = i;
		}
	}
}

template<typename TFloat>
void
findMinTimeInTwoSteps2(
		unsigned int blockSize,
		TFloat invC1, TFloat invC2,
		const std::vector<XZ<TFloat>>& interfacePointList,
		TFloat x1, TFloat z1, TFloat x2, TFloat z2,
		TFloat& tMin, unsigned int& idxMin)
{
	// First step: step = blockSize
	tMin = std::numeric_limits<TFloat>::max();
	idxMin = std::numeric_limits<unsigned int>::max();
	for (unsigned int i = 0, end = interfacePointList.size(); i < end; i += blockSize) {
		const XZ<TFloat>& point = interfacePointList[i];
		const TFloat t = SIMD::calcTwoMediumTravelTime(x1, z1, point.x, point.z, x2, z2, invC1, invC2);
		if (t < tMin) {
			tMin = t;
			idxMin = i;
		}
	}

	// Second step: step = 1
	const unsigned int iBegin = (idxMin > (blockSize - 1)) ? idxMin - (blockSize - 1) : 0;
	const unsigned int iEnd = std::min<unsigned int>(idxMin + blockSize, interfacePointList.size());
	for (unsigned int i = iBegin; i < iEnd; ++i) {
		const XZ<TFloat>& point = interfacePointList[i];
		const TFloat t = SIMD::calcTwoMediumTravelTime(x1, z1, point.x, point.z, x2, z2, invC1, invC2);
		if (t < tMin) {
			tMin = t;
			idxMin = i;
		}
	}
}

template<typename TFloat>
void
findMinTimeInTwoSteps(
		unsigned int blockSize,
		TFloat c1, TFloat c2X, TFloat c2Z,
		const std::vector<XZ<TFloat>>& interfacePointList,
		TFloat x1, TFloat z1, TFloat x2, TFloat z2,
		TFloat& tMin, unsigned int& idxMin)
{
	const TFloat invC1 = 1 / c1;

	// First step: step = blockSize
	tMin = std::numeric_limits<TFloat>::max();
	idxMin = std::numeric_limits<unsigned int>::max();
	for (unsigned int i = 0, end = interfacePointList.size(); i < end; i += blockSize) {
		const XZ<TFloat>& point = interfacePointList[i];
		const TFloat dx1 = point.x - x1;
		const TFloat dz1 = point.z - z1;
		const TFloat dx2 = x2 - point.x;
		const TFloat dz2 = z2 - point.z;
		const TFloat r2 = std::sqrt(dx2 * dx2 + dz2 * dz2);
		const TFloat cX = (dx2 / r2) * c2X;
		const TFloat cZ = (dz2 / r2) * c2Z;
		const TFloat c2 = std::sqrt(cX * cX + cZ * cZ);
		const TFloat t = std::sqrt(dx1 * dx1 + dz1 * dz1) * invC1 + r2 / c2;
		if (t < tMin) {
			tMin = t;
			idxMin = i;
		}
	}

	// Second step: step = 1
	const unsigned int iBegin = (idxMin > (blockSize - 1)) ? idxMin - (blockSize - 1) : 0;
	const unsigned int iEnd = std::min<unsigned int>(idxMin + blockSize, interfacePointList.size());
	for (unsigned int i = iBegin; i < iEnd; ++i) {
		const XZ<TFloat>& point = interfacePointList[i];
		const TFloat dx1 = point.x - x1;
		const TFloat dz1 = point.z - z1;
		const TFloat dx2 = x2 - point.x;
		const TFloat dz2 = z2 - point.z;
		const TFloat r2 = std::sqrt(dx2 * dx2 + dz2 * dz2);
		const TFloat cX = (dx2 / r2) * c2X;
		const TFloat cZ = (dz2 / r2) * c2Z;
		const TFloat c2 = std::sqrt(cX * cX + cZ * cZ);
		const TFloat t = std::sqrt(dx1 * dx1 + dz1 * dz1) * invC1 + r2 / c2;
		if (t < tMin) {
			tMin = t;
			idxMin = i;
		}
	}
}

template<typename TFloat>
void
findNearestPointInXInTwoSteps(
		unsigned int blockSize,
		const std::vector<XZ<TFloat>>& interfacePointList,
		TFloat x,
		TFloat& zIdxMin, unsigned int& idxMin)
{
	// First step: step = blockSize
	TFloat dxMin = std::numeric_limits<TFloat>::max();
	idxMin = 0;
	for (unsigned int i = 0; i < interfacePointList.size(); i += blockSize) {
		const TFloat dx = std::abs(interfacePointList[i].x - x);
		if (dx < dxMin) {
			dxMin = dx;
			idxMin = i;
		}
	}

	// Second step: step = 1
	const unsigned int iBegin = (idxMin > (blockSize - 1)) ? idxMin - (blockSize - 1) : 0;
	const unsigned int iEnd = std::min<unsigned int>(idxMin + blockSize, interfacePointList.size());
	for (unsigned int i = iBegin; i < iEnd; ++i) {
		const TFloat dx = std::abs(interfacePointList[i].x - x);
		if (dx < dxMin) {
			dxMin = dx;
			idxMin = i;
		}
	}

	zIdxMin = interfacePointList[idxMin].z;
}

} // namespace FermatPrinciple
} // namespace Lab

#endif /* FERMATPRINCIPLE_H_ */
