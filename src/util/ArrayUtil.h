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
#ifndef ARRAYUTIL_H
#define ARRAYUTIL_H

#include <limits>
#include <vector>

#include "Exception.h"
#include "Geometry.h"
#include "ParameterMap.h"
#include "XY.h"



namespace Lab {
namespace ArrayUtil {

template<typename TFloat> void calculateElementPositions(
					TFloat pitchX, unsigned int numElemX, TFloat offsetX,
					TFloat pitchY, unsigned int numElemY, TFloat offsetY,
					std::vector<XY<TFloat>>& elemPos);

template<typename TFloat> void calculateTxElementPositions(const ParameterMap& pm, std::vector<XY<TFloat>>& elemPos);
template<typename TFloat> void calculateRxElementPositions(const ParameterMap& pm, std::vector<XY<TFloat>>& elemPos);

template<typename TFloat> void calculateTx3DFocusDelay(
					TFloat focusX, TFloat focusY, TFloat focusZ, TFloat propagationSpeed,
					const std::vector<XY<TFloat>>& elemPos, std::vector<TFloat>& focusDelay /* s */);
template<typename TFloat> void calculateTx3DFocusDelay(
					TFloat focusX, TFloat focusY, TFloat focusZ, TFloat propagationSpeed,
					const std::vector<XY<TFloat>>& elemPos, unsigned int baseElement,
					unsigned int numGroupElements, std::vector<TFloat>& focusDelay /* s */);

//==============================================================================

template<typename TFloat>
void
calculateElementPositions(
		TFloat pitchX, unsigned int numElemX, TFloat offsetX,
		TFloat pitchY, unsigned int numElemY, TFloat offsetY,
		std::vector<XY<TFloat>>& elemPos)
{
	// Calculate the center of each element.
	const TFloat halfW = (numElemX - 1) * 0.5 * pitchX;
	const TFloat halfH = (numElemY - 1) * 0.5 * pitchY;
	elemPos.resize(numElemX * numElemY);
	for (unsigned int iy = 0; iy < numElemY; ++iy) {
		for (unsigned int ix = 0; ix < numElemX; ++ix) {
			XY<TFloat>& pos = elemPos[iy * numElemX + ix];
			pos.x = ix * pitchX - halfW + offsetX;
			pos.y = iy * pitchY - halfH + offsetY;
		}
	}
}

template<typename TFloat>
void
calculateTxElementPositions(const ParameterMap& pm, std::vector<XY<TFloat>>& elemPos)
{
	const auto pitchX   = pm.value<TFloat>("tx_pitch_x"   ,  1.0e-6, 1000.0);
	const auto pitchY   = pm.value<TFloat>("tx_pitch_y"   ,  1.0e-6, 1000.0);
	const auto numElemX = pm.value<TFloat>("tx_num_elem_x",       1,   1024);
	const auto numElemY = pm.value<TFloat>("tx_num_elem_y",       1,   1024);
	const auto offsetX  = pm.value<TFloat>("tx_offset_x"  , -1000.0, 1000.0);
	const auto offsetY  = pm.value<TFloat>("tx_offset_y"  , -1000.0, 1000.0);

	calculateElementPositions(pitchX, numElemX, offsetX, pitchY, numElemY, offsetY, elemPos);
}

template<typename TFloat>
void
calculateRxElementPositions(const ParameterMap& pm, std::vector<XY<TFloat>>& elemPos)
{
	const auto pitchX   = pm.value<TFloat>("rx_pitch_x"   ,  1.0e-6, 1000.0);
	const auto pitchY   = pm.value<TFloat>("rx_pitch_y"   ,  1.0e-6, 1000.0);
	const auto numElemX = pm.value<TFloat>("rx_num_elem_x",       1,   1024);
	const auto numElemY = pm.value<TFloat>("rx_num_elem_y",       1,   1024);
	const auto offsetX  = pm.value<TFloat>("rx_offset_x"  , -1000.0, 1000.0);
	const auto offsetY  = pm.value<TFloat>("rx_offset_y"  , -1000.0, 1000.0);

	calculateElementPositions(pitchX, numElemX, offsetX, pitchY, numElemY, offsetY, elemPos);
}

template<typename TFloat>
void
calculateTx3DFocusDelay(TFloat focusX, TFloat focusY, TFloat focusZ, TFloat propagationSpeed,
			const std::vector<XY<TFloat>>& elemPos, std::vector<TFloat>& focusDelay /* s */)
{
	calculateTx3DFocusDelay(focusX, focusY, focusZ, propagationSpeed,
				elemPos, 0, elemPos.size(), focusDelay);
}

template<typename TFloat>
void
calculateTx3DFocusDelay(TFloat focusX, TFloat focusY, TFloat focusZ, TFloat propagationSpeed,
			const std::vector<XY<TFloat>>& elemPos, unsigned int baseElement,
			unsigned int numGroupElements, std::vector<TFloat>& focusDelay /* s */)
{
	if (baseElement + numGroupElements > elemPos.size()) {
		THROW_EXCEPTION(InvalidParameterException, "Error: baseElement + numGroupElements > elemPos.size() ("
				<< " baseElement=" << baseElement << " numGroupElements=" << numGroupElements
				<< " elemPos.size()=" << elemPos.size() << ").");
	}
	focusDelay.assign(numGroupElements, 0.0);

	const TFloat invC = 1 / propagationSpeed;
	if (focusZ > 0.0) {
		TFloat maxDt = 0.0;
		for (unsigned int i = 0, iEnd = focusDelay.size(); i < iEnd; ++i) {
			const XY<TFloat>& pos = elemPos[baseElement + i];
			const TFloat dt = Geometry::distance3DZ0(pos.x, pos.y, focusX, focusY, focusZ) * invC;
			if (dt > maxDt) maxDt = dt;
			focusDelay[i] = dt;
		}
		for (unsigned int i = 0, iEnd = focusDelay.size(); i < iEnd; ++i) {
			focusDelay[i] = maxDt - focusDelay[i];
		}
	} else {
		TFloat minDt = std::numeric_limits<TFloat>::max();
		for (unsigned int i = 0, iEnd = focusDelay.size(); i < iEnd; ++i) {
			const XY<TFloat>& pos = elemPos[baseElement + i];
			const TFloat dt = Geometry::distance3DZ0(pos.x, pos.y, focusX, focusY, focusZ) * invC;
			if (dt < minDt) minDt = dt;
			focusDelay[i] = dt;
		}
		for (unsigned int i = 0, iEnd = focusDelay.size(); i < iEnd; ++i) {
			focusDelay[i] -= minDt;
		}
	}
}

} // namespace ArrayUtil
} // namespace Lab

#endif // ARRAYUTIL_H
