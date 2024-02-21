/***************************************************************************
 *  Copyright 2018, 2019 Marcelo Y. Matuda                                 *
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

#ifndef SIMRECTANGULARSOURCEMETHOD_H
#define SIMRECTANGULARSOURCEMETHOD_H

#include <string>
#include <vector>

#include "Method.h"
#include "XY.h"



namespace Lab {

class ParameterMap;
class Project;

template<typename TFloat>
class SimRectangularSourceMethod : public Method {
public:
	SimRectangularSourceMethod(Project& project);
	virtual ~SimRectangularSourceMethod() = default;

	virtual void execute();
private:
	struct MainData {
		TFloat propagationSpeed;
		TFloat centerFreq;
		TFloat maxFreq;
		TFloat nyquistRate;
		std::string outputDir;
		unsigned int numThreadsNumeric;
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
		TFloat sourceWidth;
		TFloat sourceHeight;

		// For arrays.
		TFloat focusX;
		TFloat focusY;
		TFloat focusZ;
		std::vector<XY<TFloat>> elemPos;
		std::vector<TFloat> focusDelay;
		bool useFocus;
	};

	SimRectangularSourceMethod(const SimRectangularSourceMethod&) = delete;
	SimRectangularSourceMethod& operator=(const SimRectangularSourceMethod&) = delete;
	SimRectangularSourceMethod(SimRectangularSourceMethod&&) = delete;
	SimRectangularSourceMethod& operator=(SimRectangularSourceMethod&&) = delete;

	void loadData(const ParameterMap& taskPM, MainData& data, SimulationData& simData);
	void loadSourceData(MainData& data, bool sourceIsArray, SourceData& srcData);
	void loadSimulationData(const MainData& data, const std::string& irMethod, SimulationData& simData);
	void prepareExcitation(TFloat dt, const SimulationData& simData, std::vector<TFloat>& tExc,
				std::vector<TFloat>& dvdt, std::vector<TFloat>& tDvdt);

	void execImpulseResponse(bool sourceIsArray); // calculate p/(c*density)
	void execTransientRadiationPattern(bool sourceIsArray);
	void execTransientAcousticField(bool sourceIsArray);
	void execTransientPropagation(bool sourceIsArray);

	Project& project_;
};

} // namespace Lab

#endif // SIMRECTANGULARSOURCEMETHOD_H
