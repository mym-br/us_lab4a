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
template<typename FloatType>
struct TnRnConfiguration {
	TnRnConfiguration() { }
	~TnRnConfiguration() { }

	TnRnConfiguration(ConstParameterMapPtr imgPM, ConstParameterMapPtr arrayPM) { load(imgPM, arrayPM); }

	unsigned int numElements;
	unsigned int numElementsMux;
	unsigned int numPulses;
	FloatType centerFrequency; // Hz
	FloatType maxFrequency; // Hz
	FloatType acquisitionTime; // s
	FloatType minGain; // dB
	FloatType maxGain; // dB
	FloatType propagationSpeed; // m/s
	FloatType acquisitionDelay; // s
	FloatType samplingFrequency; // Hz
	FloatType deadZoneM; // m
	FloatType valueScale;

	std::vector<XY<FloatType>> txElemPos; // m
	std::vector<XY<FloatType>> rxElemPos; // m

	void load(ConstParameterMapPtr imgPM, ConstParameterMapPtr arrayPM);
};

template<typename FloatType>
void
TnRnConfiguration<FloatType>::load(ConstParameterMapPtr imgPM, ConstParameterMapPtr arrayPM)
{
	numElementsMux    = imgPM->value<unsigned int>("num_elements_mux"   ,       1,    1024);
	numElements       = imgPM->value<unsigned int>("num_elements"       ,       1, numElementsMux);
	numPulses         = imgPM->value<unsigned int>("num_pulses"         ,       1,     100);
	centerFrequency   = imgPM->value<FloatType>(   "center_frequency"   ,   100.0, 100.0e6);
	maxFrequency      = imgPM->value<FloatType>(   "max_frequency"      ,   100.0, 100.0e6);
	acquisitionTime   = imgPM->value<FloatType>(   "acquisition_time"   ,  1.0e-6,     1.0);
	minGain           = imgPM->value<FloatType>(   "min_gain"           , -2000.0,  2000.0);
	maxGain           = imgPM->value<FloatType>(   "max_gain"           , minGain,  2000.0);
	propagationSpeed  = imgPM->value<FloatType>(   "propagation_speed_1",   100.0, 10000.0);
	acquisitionDelay  = imgPM->value<FloatType>(   "acquisition_delay"  ,     0.0,     1.0);
	samplingFrequency = imgPM->value<FloatType>(   "sampling_frequency" ,   100.0, 200.0e6);
	deadZoneM         = imgPM->value<FloatType>(   "dead_zone_m"        ,     0.0, 50.0e-3);
	valueScale        = imgPM->value<FloatType>(   "value_scale"        ,     0.0,  1.0e30);

	ArrayUtil::calculateTxElementPositions(*arrayPM, txElemPos);
	ArrayUtil::calculateRxElementPositions(*arrayPM, rxElemPos);
	if (txElemPos.size() > numElementsMux || rxElemPos.size() > numElementsMux) {
		THROW_EXCEPTION(InvalidParameterException, "Error: numElementsMux is less than the number of array elements.");
	}
}

} // namespace Lab

#endif // TNRNCONFIGURATION_H
