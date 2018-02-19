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

#ifndef STACONFIGURATION_H_
#define STACONFIGURATION_H_

#include "ParameterMap.h"



namespace Lab {

// Configuration for Synthetic Transmit Aperture - one medium.
template<typename FloatType>
struct STAConfiguration {
	STAConfiguration() { }
	~STAConfiguration() { }
	STAConfiguration(ConstParameterMapPtr pm) { load(pm); }

	unsigned int numElements;
	unsigned int numElementsMux;
	FloatType pitch; // m
	FloatType centerFrequency; // Hz
	FloatType acquisitionTime; // s
	FloatType minGain; // (0 - 48) dB
	FloatType maxGain; // (0 - 48) dB
	FloatType propagationSpeed; // m/s
	FloatType acquisitionDelay; // s
	FloatType samplingFrequency; // Hz
	FloatType deadZoneM; // m

	void load(ConstParameterMapPtr pm);
};



template<typename FloatType>
void
STAConfiguration<FloatType>::load(ConstParameterMapPtr pm)
{
	numElementsMux    = pm->value<unsigned int>("num_elements_mux"   ,       1,    1024);
	numElements       = pm->value<unsigned int>("num_elements"       ,       1, numElementsMux);
	pitch             = pm->value<FloatType>(   "pitch"              , 0.01e-3, 10.0e-3);
	centerFrequency   = pm->value<FloatType>(   "center_frequency"   ,   100.0, 100.0e6);
	acquisitionTime   = pm->value<FloatType>(   "acquisition_time"   ,  1.0e-6,     1.0);
	minGain           = pm->value<FloatType>(   "min_gain"           ,     0.0,    48.0);
	maxGain           = pm->value<FloatType>(   "max_gain"           ,     0.0,    48.0);
	propagationSpeed  = pm->value<FloatType>(   "propagation_speed_1",   100.0, 10000.0);
	acquisitionDelay  = pm->value<FloatType>(   "acquisition_delay"  ,     0.0,     1.0);
	samplingFrequency = pm->value<FloatType>(   "sampling_frequency" ,   100.0, 200.0e6);
	deadZoneM         = pm->value<FloatType>(   "dead_zone_m"        ,     0.0, 50.0e-3);
}

} // namespace Lab

#endif /* STACONFIGURATION_H_ */
