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

#ifndef SA3DCONFIGURATION_H
#define SA3DCONFIGURATION_H

#include <sstream>
#include <string>
#include <vector>

#include "ArrayUtil.h"
#include "Exception.h"
#include "Log.h"
#include "ParameterMap.h"
#include "XY.h"



namespace Lab {

// Configuration for 3D Synthetic Aperture - one medium.
template<typename FloatType>
struct SA3DConfiguration {
	SA3DConfiguration() = default;
	SA3DConfiguration(ParamMapPtr saPM, ParamMapPtr arrayPM) { load(saPM, arrayPM); }

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

	std::vector<unsigned int> activeTxElem; // relative to the base element
	std::vector<unsigned int> activeRxElem; // relative to the base element
	std::vector<XY<FloatType>> txElemPos; // m
	std::vector<XY<FloatType>> rxElemPos; // m

	void load(ParamMapPtr saPM, ParamMapPtr arrayPM);

private:
	void fillActiveElem(const std::string& listStr, unsigned int maxElem, std::vector<unsigned int>& activeElem);
};

template<typename FloatType>
void
SA3DConfiguration<FloatType>::fillActiveElem(const std::string& listStr, unsigned int maxElem,
						std::vector<unsigned int>& activeElem)
{
	std::istringstream in(listStr);
	std::string item;
	while (in >> item) {
		auto pos = item.find(':');
		if (pos == std::string::npos) {
			unsigned long elem;
			try {
				elem = std::stoul(item);
			} catch (...) {
				THROW_EXCEPTION(InvalidValueException, "Missing element in the list of active elements.");
			}
			if (elem > maxElem) {
				THROW_EXCEPTION(InvalidValueException, "Invalid element in the list of active elements.");
			}
			LOG_DEBUG << "elem: " << elem;
			activeElem.push_back(static_cast<unsigned int>(elem));
		} else {
			unsigned long elem1, elem2;
			try {
				elem1 = std::stoul(item);
				elem2 = std::stoul(item.substr(pos + 1));
			} catch (...) {
				THROW_EXCEPTION(InvalidValueException, "Missing element in the list of active elements.");
			}
			if (elem1 > maxElem || elem2 > maxElem) {
				THROW_EXCEPTION(InvalidValueException, "Invalid element in the list of active elements.");
			}
			for (unsigned int elem = elem1; elem <= elem2; ++elem) {
				LOG_DEBUG << "elem: " << elem;
				activeElem.push_back(static_cast<unsigned int>(elem));
			}
		}
	}
	if (activeElem.empty()) {
		THROW_EXCEPTION(InvalidValueException, "Empty list of active elements.");
	}
}

template<typename FloatType>
void
SA3DConfiguration<FloatType>::load(ParamMapPtr saPM, ParamMapPtr arrayPM)
{
	saPM->getValue(numElementsMux   , "num_elements_mux"   ,       1,    1024);
	saPM->getValue(numPulses        , "num_pulses"         ,       1,     100);
	saPM->getValue(centerFrequency  , "center_frequency"   ,   100.0, 100.0e6);
	saPM->getValue(maxFrequency     , "max_frequency"      ,   100.0, 100.0e6);
	saPM->getValue(acquisitionTime  , "acquisition_time"   ,  1.0e-6,     1.0);
	saPM->getValue(minGain          , "min_gain"           , -2000.0,  2000.0);
	saPM->getValue(maxGain          , "max_gain"           , minGain,  2000.0);
	saPM->getValue(propagationSpeed , "propagation_speed_1",   100.0, 10000.0);
	saPM->getValue(acquisitionDelay , "acquisition_delay"  ,     0.0,     1.0);
	saPM->getValue(samplingFrequency, "sampling_frequency" ,   100.0, 200.0e6);
	saPM->getValue(deadZoneM        , "dead_zone_m"        ,     0.0, 50.0e-3);
	saPM->getValue(valueScale       , "value_scale"        ,     0.0,  1.0e30);

	ArrayUtil::calculateTxElementPositions(*arrayPM, txElemPos);
	ArrayUtil::calculateRxElementPositions(*arrayPM, rxElemPos);
	if (txElemPos.size() > numElementsMux || rxElemPos.size() > numElementsMux) {
		THROW_EXCEPTION(InvalidParameterException, "Error: numElementsMux is less than the number of array elements.");
	}

	LOG_DEBUG << "Active tx elements:";
	fillActiveElem(saPM->value<std::string>("active_tx_elem"), txElemPos.size() - 1U, activeTxElem);
	LOG_DEBUG << "Active rx elements:";
	fillActiveElem(saPM->value<std::string>("active_rx_elem"), rxElemPos.size() - 1U, activeRxElem);
}

} // namespace Lab

#endif // SA3DCONFIGURATION_H
