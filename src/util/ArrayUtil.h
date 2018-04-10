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

#include "ParameterMap.h"
#include "XY.h"



namespace Lab {
namespace ArrayUtil {

template<typename FloatType> void calculateXYArrayParameters(const ParameterMap& taskPM, FloatType propagationSpeed,
						FloatType samplingFreq, std::vector<XY<FloatType>>& elemPos,
						std::vector<FloatType>& focusDelay);



template<typename FloatType>
void
calculateXYArrayParameters(
		const ParameterMap& taskPM,
		FloatType propagationSpeed,
		FloatType samplingFreq,
		std::vector<XY<FloatType>>& elemPos,
		std::vector<FloatType>& focusDelay)
{
	const bool useFocus              = taskPM.value<bool>("use_focus");
	const FloatType focusX           = taskPM.value<FloatType>("focus_x", 0.0, 10000.0);
	const FloatType focusY           = taskPM.value<FloatType>("focus_y", 0.0, 10000.0);
	const FloatType focusZ           = taskPM.value<FloatType>("focus_z", 0.0, 10000.0);
	const FloatType arrayPitchX      = taskPM.value<FloatType>("array_pitch_x", 0.0, 10.0);
	const unsigned int numArrayElemX = taskPM.value<FloatType>("num_array_elem_x", 1, 1024);
	const FloatType arrayPitchY      = taskPM.value<FloatType>("array_pitch_y", 0.0, 10.0);
	const unsigned int numArrayElemY = taskPM.value<FloatType>("num_array_elem_y", 1, 1024);

	const unsigned int numElem = numArrayElemX * numArrayElemY;

	// Calculate the center of each element.
	const FloatType halfW = (numArrayElemX - 1) * 0.5 * arrayPitchX;
	const FloatType halfH = (numArrayElemY - 1) * 0.5 * arrayPitchY;
	elemPos.resize(numElem);
	for (unsigned int iy = 0; iy < numArrayElemY; ++iy) {
		for (unsigned int ix = 0; ix < numArrayElemX; ++ix) {
			XY<FloatType>& pos = elemPos[iy * numArrayElemX + ix];
			pos.x = ix * arrayPitchX - halfW;
			pos.y = iy * arrayPitchY - halfH;
		}
	}

	if (useFocus) {
		focusDelay.resize(numElem);
		FloatType maxDt = 0.0;
		for (unsigned int iy = 0; iy < numArrayElemY; ++iy) {
			for (unsigned int ix = 0; ix < numArrayElemX; ++ix) {
				const unsigned int index = iy * numArrayElemX + ix;
				XY<FloatType>& pos = elemPos[index];
				const FloatType dx = focusX - pos.x;
				const FloatType dy = focusY - pos.y;
				const FloatType focusDt = std::sqrt(dx * dx + dy * dy + focusZ * focusZ) / propagationSpeed;
				focusDelay[index] = focusDt;
				if (focusDt > maxDt) maxDt = focusDt;
			}
		}
		for (unsigned int iy = 0; iy < numArrayElemY; ++iy) {
			for (unsigned int ix = 0; ix < numArrayElemX; ++ix) {
				const unsigned int index = iy * numArrayElemX + ix;
				focusDelay[index] = (maxDt - focusDelay[index]) * samplingFreq;
			}
		}
	} else {
		focusDelay.assign(numElem, 0.0);
	}
}

} // namespace ArrayUtil
} // namespace Lab

#endif // ARRAYUTIL_H
