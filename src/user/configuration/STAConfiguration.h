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
	STAConfiguration(ParamMapPtr pm) { load(pm); }

	unsigned int numElements;
	unsigned int numElementsMux;
	unsigned int firstTxElem;
	unsigned int lastTxElem;
	unsigned int numPulses;
	FloatType pitch; // m
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

	void load(ParamMapPtr pm);
};



template<typename FloatType>
void
STAConfiguration<FloatType>::load(ParamMapPtr pm)
{
	pm->getValue(numElementsMux   , "num_elements_mux"   ,       1,    1024);
	pm->getValue(numElements      , "num_elements"       ,       1, numElementsMux);
	pm->getValue(firstTxElem      , "first_tx_elem"      ,       0, numElements - 1);
	pm->getValue(lastTxElem       , "last_tx_elem"       , firstTxElem, numElements - 1);
	pm->getValue(numPulses        , "num_pulses"         ,       1,     100);
	pm->getValue(pitch            , "pitch"              ,  1.0e-6,  1000.0);
	pm->getValue(centerFrequency  , "center_frequency"   ,   100.0, 100.0e6);
	pm->getValue(maxFrequency     , "max_frequency"      ,   100.0, 100.0e6);
	pm->getValue(acquisitionTime  , "acquisition_time"   ,  1.0e-6,     1.0);
	pm->getValue(minGain          , "min_gain"           , -2000.0,  2000.0);
	pm->getValue(maxGain          , "max_gain"           , minGain,  2000.0);
	pm->getValue(propagationSpeed , "propagation_speed_1",   100.0, 10000.0);
	pm->getValue(acquisitionDelay , "acquisition_delay"  ,     0.0,     1.0);
	pm->getValue(samplingFrequency, "sampling_frequency" ,   100.0, 200.0e6);
	pm->getValue(deadZoneM        , "dead_zone_m"        ,     0.0, 50.0e-3);
	pm->getValue(valueScale       , "value_scale"        ,     0.0,  1.0e30);
}

} // namespace Lab

#endif /* STACONFIGURATION_H_ */
