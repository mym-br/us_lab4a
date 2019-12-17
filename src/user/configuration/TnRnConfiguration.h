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

#ifndef TNRNCONFIGURATION_H
#define TNRNCONFIGURATION_H

#include <vector>

#include "ArrayUtil.h"
#include "Exception.h"
#include "ParameterMap.h"
#include "XY.h"



namespace Lab {

// Configuration for imaging using all the active group elements in transmission and reception - one medium.
template<typename TFloat>
struct TnRnConfiguration {
	TnRnConfiguration() = default;
	TnRnConfiguration(const ParameterMap& imgPM, const ParameterMap& arrayPM) { load(imgPM, arrayPM); }

	unsigned int numElements;
	unsigned int numElementsMux;
	unsigned int numPulses;
	TFloat centerFrequency; // Hz
	TFloat maxFrequency; // Hz
	TFloat acquisitionTime; // s
	TFloat minGain; // dB
	TFloat maxGain; // dB
	TFloat propagationSpeed; // m/s
	TFloat acquisitionDelay; // s
	TFloat samplingFrequency; // Hz
	TFloat deadZoneM; // m
	TFloat valueScale;

	std::vector<XY<TFloat>> txElemPos; // m
	std::vector<XY<TFloat>> rxElemPos; // m

	void load(const ParameterMap& imgPM, const ParameterMap& arrayPM);
};

template<typename TFloat>
void
TnRnConfiguration<TFloat>::load(const ParameterMap& imgPM, const ParameterMap& arrayPM)
{
	imgPM.getValue(numElementsMux   , "num_elements_mux"   ,       1,    1024);
	imgPM.getValue(numElements      , "num_elements"       ,       1, numElementsMux);
	imgPM.getValue(numPulses        , "num_pulses"         ,       1,     100);
	imgPM.getValue(centerFrequency  , "center_frequency"   ,   100.0, 100.0e6);
	imgPM.getValue(maxFrequency     , "max_frequency"      ,   100.0, 100.0e6);
	imgPM.getValue(acquisitionTime  , "acquisition_time"   ,  1.0e-6,     1.0);
	imgPM.getValue(minGain          , "min_gain"           , -2000.0,  2000.0);
	imgPM.getValue(maxGain          , "max_gain"           , minGain,  2000.0);
	imgPM.getValue(propagationSpeed , "propagation_speed_1",   100.0, 10000.0);
	imgPM.getValue(acquisitionDelay , "acquisition_delay"  ,     0.0,     1.0);
	imgPM.getValue(samplingFrequency, "sampling_frequency" ,   100.0, 200.0e6);
	imgPM.getValue(deadZoneM        , "dead_zone_m"        ,     0.0, 50.0e-3);
	imgPM.getValue(valueScale       , "value_scale"        ,     0.0,  1.0e30);

	ArrayUtil::calculateTxElementPositions(arrayPM, txElemPos);
	ArrayUtil::calculateRxElementPositions(arrayPM, rxElemPos);
	if (txElemPos.size() > numElementsMux || rxElemPos.size() > numElementsMux) {
		THROW_EXCEPTION(InvalidParameterException, "Error: numElementsMux is less than the number of array elements.");
	}
}

} // namespace Lab

#endif // TNRNCONFIGURATION_H
