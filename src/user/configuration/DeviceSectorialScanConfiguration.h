/***************************************************************************
 *  Copyright 2014, 2017, 2018 Marcelo Y. Matuda                           *
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

#ifndef DEVICESECTORIALSCANCONFIGURATION_H
#define DEVICESECTORIALSCANCONFIGURATION_H

#include "ParameterMap.h"



namespace Lab {

template<typename TFloat>
struct DeviceSectorialScanConfiguration {
	unsigned int numElements;
	unsigned int numElementsMux;
	unsigned int baseElement;      // first element: 0
	unsigned int signalMode;       // 0: Raw / 1: Envelope
	TFloat pitch;                  // m
	TFloat centerFrequency;        // Hz
	TFloat gain;                   // dB
	TFloat propagationSpeed;       // m/s
	TFloat acquisitionDelay;       // s
	TFloat samplingFrequency;      // Hz
	TFloat focalEmissionDistance;  // m
	TFloat focalReceptionDistance; // m
	TFloat valueScale;

	// Sectorial scan grid.
	TFloat rangeStart; // m
	TFloat rangeEnd;   // m
	TFloat startAngle; // degree
	TFloat endAngle;   // degree
	TFloat angleStep;  // degree

	bool enableFocusing;

	void load(const ParameterMap& pm);
};



template<typename TFloat>
void
DeviceSectorialScanConfiguration<TFloat>::load(const ParameterMap& pm)
{
	pm.getValue(numElementsMux        , "num_elements_mux"        ,       8,      1024);
	pm.getValue(numElements           , "num_elements"            ,       8, numElementsMux);
	pm.getValue(baseElement           , "base_element"            ,       0, numElementsMux - numElements);
	pm.getValue(signalMode            , "signal_mode"             ,       0,         1);
	pm.getValue(pitch                 , "pitch"                   , 0.01e-3,   10.0e-3);
	pm.getValue(centerFrequency       , "center_frequency"        ,   100.0,   100.0e6);
	pm.getValue(gain                  , "gain"                    ,     0.0,     100.0);
	pm.getValue(propagationSpeed      , "propagation_speed"       ,   100.0,   10000.0);
	pm.getValue(acquisitionDelay      , "acquisition_delay"       ,     0.0,       1.0);
	pm.getValue(samplingFrequency     , "sampling_frequency"      ,   100.0,   200.0e6);
	pm.getValue(focalEmissionDistance , "focus_emission_distance" ,  1.0e-3, 1000.0e-3);
	pm.getValue(focalReceptionDistance, "focus_reception_distance",  1.0e-3, 1000.0e-3);
	pm.getValue(valueScale            , "value_scale"             ,     0.0,    1.0e10);

	// Sectorial scan grid.
	pm.getValue(rangeStart            , "range_start"             ,              1.0e-3,  99.0);
	pm.getValue(rangeEnd              , "range_end"               , rangeStart + 1.0e-2, 100.0);
	pm.getValue(startAngle            , "start_angle"             ,               -90.0,  89.0);
	pm.getValue(endAngle              , "end_angle"               ,    startAngle + 1.0,  90.0);
	pm.getValue(angleStep             , "angle_step"              ,              1.0e-3,  10.0);

	pm.getValue(enableFocusing        , "enable_focusing");

	if (numElementsMux & 0x01) {
		THROW_EXCEPTION(InvalidParameterException, "The value of num_elements_mux is not even.");
	}
	if (numElements & 0x01) {
		THROW_EXCEPTION(InvalidParameterException, "The value of num_elements is not even.");
	}
}

} // namespace Lab

#endif // DEVICESECTORIALSCANCONFIGURATION_H
