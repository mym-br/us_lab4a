/***************************************************************************
 *  Copyright 2019, 2020 Marcelo Y. Matuda                                 *
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

#ifndef SIMCIRCULARSOURCEMETHOD_H
#define SIMCIRCULARSOURCEMETHOD_H

#include <string>
#include <vector>

#include "Method.h"



namespace Lab {

class ParameterMap;
class Project;

template<typename TFloat>
class SimCircularSourceMethod : public Method {
public:
	SimCircularSourceMethod(Project& project);
	virtual ~SimCircularSourceMethod() = default;

	virtual void execute();
private:
	struct MainData {
		TFloat propagationSpeed;
		TFloat centerFreq;
		TFloat maxFreq;
		TFloat nyquistRate;
		std::string outputDir;
	};
	struct SimulationData {
		TFloat samplingFreq;
		TFloat excNumPeriods;
		TFloat discretFactor;
		std::string irMethod;
		std::string excitationType;
		std::vector<TFloat> exc;
	};
	struct SourceData {
		TFloat sourceRadius;
	};

	SimCircularSourceMethod(const SimCircularSourceMethod&) = delete;
	SimCircularSourceMethod& operator=(const SimCircularSourceMethod&) = delete;
	SimCircularSourceMethod(SimCircularSourceMethod&&) = delete;
	SimCircularSourceMethod& operator=(SimCircularSourceMethod&&) = delete;

	void loadData(const ParameterMap& taskPM, MainData& data, SimulationData& simData);
	void loadSourceData(SourceData& srcData);
	void loadSimulationData(const MainData& data, const std::string& irMethod, SimulationData& simData);
	void prepareExcitation(TFloat dt, const SimulationData& simData, std::vector<TFloat>& tExc,
				std::vector<TFloat>& dvdt, std::vector<TFloat>& tDvdt);

	void execImpulseResponse(); // calculate p/(c*density)
	void execTransientAcousticField();
	void execTransientPropagation();
	void execTransientRadiationPattern();

	Project& project_;
};

} // namespace Lab

#endif // SIMCIRCULARSOURCEMETHOD_H
