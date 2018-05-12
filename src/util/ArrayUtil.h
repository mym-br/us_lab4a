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

#include <vector>

#include "Geometry.h"
#include "ParameterMap.h"
#include "XY.h"



namespace Lab {
namespace ArrayUtil {

template<typename FloatType> void calculateElementPositions(
					FloatType arrayPitchX, unsigned int numArrayElemX,
					FloatType arrayPitchY, unsigned int numArrayElemY,
					std::vector<XY<FloatType>>& elemPos);

template<typename FloatType> void calculateTxElementPositions(const ParameterMap& pm, std::vector<XY<FloatType>>& elemPos);
template<typename FloatType> void calculateRxElementPositions(const ParameterMap& pm, std::vector<XY<FloatType>>& elemPos);

template<typename FloatType> void calculateTx3DFocusDelay(
					const ParameterMap& pm, FloatType propagationSpeed,
					const std::vector<XY<FloatType>>& elemPos, std::vector<FloatType>& focusDelay /* s */);

//==============================================================================

template<typename FloatType>
void
calculateElementPositions(
		FloatType arrayPitchX, unsigned int numArrayElemX,
		FloatType arrayPitchY, unsigned int numArrayElemY,
		std::vector<XY<FloatType>>& elemPos)
{
	// Calculate the center of each element.
	const FloatType halfW = (numArrayElemX - 1) * 0.5 * arrayPitchX;
	const FloatType halfH = (numArrayElemY - 1) * 0.5 * arrayPitchY;
	elemPos.resize(numArrayElemX * numArrayElemY);
	for (unsigned int iy = 0; iy < numArrayElemY; ++iy) {
		for (unsigned int ix = 0; ix < numArrayElemX; ++ix) {
			XY<FloatType>& pos = elemPos[iy * numArrayElemX + ix];
			pos.x = ix * arrayPitchX - halfW;
			pos.y = iy * arrayPitchY - halfH;
		}
	}
}

template<typename FloatType>
void
calculateTxElementPositions(const ParameterMap& pm, std::vector<XY<FloatType>>& elemPos)
{
	const FloatType arrayPitchX      = pm.value<FloatType>("tx_array_pitch_x", 0.0, 10.0);
	const FloatType arrayPitchY      = pm.value<FloatType>("tx_array_pitch_y", 0.0, 10.0);
	const unsigned int numArrayElemX = pm.value<FloatType>("tx_num_array_elem_x", 1, 1024);
	const unsigned int numArrayElemY = pm.value<FloatType>("tx_num_array_elem_y", 1, 1024);

	calculateElementPositions(arrayPitchX, numArrayElemX, arrayPitchY, numArrayElemY, elemPos);
}

template<typename FloatType>
void
calculateRxElementPositions(const ParameterMap& pm, std::vector<XY<FloatType>>& elemPos)
{
	const FloatType arrayPitchX      = pm.value<FloatType>("rx_array_pitch_x", 0.0, 10.0);
	const FloatType arrayPitchY      = pm.value<FloatType>("rx_array_pitch_y", 0.0, 10.0);
	const unsigned int numArrayElemX = pm.value<FloatType>("rx_num_array_elem_x", 1, 1024);
	const unsigned int numArrayElemY = pm.value<FloatType>("rx_num_array_elem_y", 1, 1024);

	calculateElementPositions(arrayPitchX, numArrayElemX, arrayPitchY, numArrayElemY, elemPos);
}

template<typename FloatType>
void
calculateTx3DFocusDelay(
		const ParameterMap& pm,
		FloatType propagationSpeed,
		const std::vector<XY<FloatType>>& elemPos,
		std::vector<FloatType>& focusDelay)
{
	focusDelay.assign(elemPos.size(), 0.0);

	const bool useFocus = pm.value<bool>("use_tx_focus");
	if (useFocus) {
		const FloatType focusX = pm.value<FloatType>("tx_focus_x", -10000.0, 10000.0);
		const FloatType focusY = pm.value<FloatType>("tx_focus_y", -10000.0, 10000.0);
		const FloatType focusZ = pm.value<FloatType>("tx_focus_z", -10000.0, 10000.0);

		const FloatType invC = 1 / propagationSpeed;
		FloatType maxDt = 0.0;
		for (unsigned int i = 0, iEnd = focusDelay.size(); i < iEnd; ++i) {
			const XY<FloatType>& pos = elemPos[i];
			const FloatType dt = Geometry::distance3DZ0(pos.x, pos.y, focusX, focusY, focusZ) * invC;
			if (dt > maxDt) maxDt = dt;
			focusDelay[i] = dt;
		}
		for (unsigned int i = 0, iEnd = focusDelay.size(); i < iEnd; ++i) {
			focusDelay[i] = maxDt - focusDelay[i];
		}
	}
}

} // namespace ArrayUtil
} // namespace Lab

#endif // ARRAYUTIL_H
