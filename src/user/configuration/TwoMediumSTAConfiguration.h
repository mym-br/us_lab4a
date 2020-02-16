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

#ifndef TWOMEDIUMSTACONFIGURATION_H_
#define TWOMEDIUMSTACONFIGURATION_H_

#include "ParameterMap.h"



namespace Lab {

// Configuration for Synthetic Transmit Aperture - two-medium.
template<typename TFloat>
struct TwoMediumSTAConfiguration {
	TwoMediumSTAConfiguration() = default;
	TwoMediumSTAConfiguration(const ParameterMap& pm) { load(pm); }

	unsigned int numElements;
	unsigned int numElementsMux;
	TFloat pitch; // m
	TFloat centerFrequency; // Hz
	TFloat acquisitionTime; // s
	TFloat minGain; // dB
	TFloat maxGain; // dB
	TFloat propagationSpeed1; // m/s
	TFloat propagationSpeed2; // m/s
	TFloat acquisitionDelay; // s
	TFloat samplingFrequency; // Hz
	TFloat deadZoneM; // m

	void load(const ParameterMap& pm);
};



template<typename TFloat>
void
TwoMediumSTAConfiguration<TFloat>::load(const ParameterMap& pm)
{
	pm.getValue(numElementsMux   , "num_elements_mux"   ,       1,    1024);
	pm.getValue(numElements      , "num_elements"       ,       1, numElementsMux);
	pm.getValue(pitch            , "pitch"              , 0.01e-3, 10.0e-3);
	pm.getValue(centerFrequency  , "center_frequency"   ,   100.0, 100.0e6);
	pm.getValue(acquisitionTime  , "acquisition_time"   ,  1.0e-6,     1.0);
	pm.getValue(minGain          , "min_gain"           ,     0.0,    48.0);
	pm.getValue(maxGain          , "max_gain"           ,     0.0,    48.0);
	pm.getValue(propagationSpeed1, "propagation_speed_1",   100.0, 10000.0);
	pm.getValue(propagationSpeed2, "propagation_speed_2",   100.0, 10000.0);
	pm.getValue(acquisitionDelay , "acquisition_delay"  ,     0.0,     1.0);
	pm.getValue(samplingFrequency, "sampling_frequency" ,   100.0, 200.0e6);
	pm.getValue(deadZoneM        , "dead_zone_m"        ,     0.0, 50.0e-3);
}

} // namespace Lab

#endif /* TWOMEDIUMSTACONFIGURATION_H_ */
