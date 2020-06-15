/***************************************************************************
 *  Copyright 2014, 2017, 2018, 2019 Marcelo Y. Matuda                     *
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

#ifndef SINGLEVIRTUALSOURCEMETHOD_H
#define SINGLEVIRTUALSOURCEMETHOD_H

#include <string>
#include <vector>

#include "Matrix.h"
#include "Method.h"
#include "XYZ.h"
#include "XYZValueFactor.h"



namespace Lab {

template<typename TPoint> class ArrayProcessor;
template<typename TFloat> class TnRnAcquisition;
class ParameterMap;
class Project;

template<typename TFloat>
class SingleVirtualSourceMethod : public Method {
public:
	SingleVirtualSourceMethod(Project& project);
	virtual ~SingleVirtualSourceMethod() = default;

	virtual void execute();
private:
	static constexpr const char* timeFile    = "/time";
	static constexpr const char* timeDataset = "time";

	SingleVirtualSourceMethod(const SingleVirtualSourceMethod&) = delete;
	SingleVirtualSourceMethod& operator=(const SingleVirtualSourceMethod&) = delete;
	SingleVirtualSourceMethod(SingleVirtualSourceMethod&&) = delete;
	SingleVirtualSourceMethod& operator=(SingleVirtualSourceMethod&&) = delete;

	void process(TFloat valueScale, ArrayProcessor<XYZValueFactor<TFloat>>& processor, unsigned int baseElement,
			bool saveCoordinates, const std::string& outputDir);
	void useCoherenceFactor(TFloat valueScale, const std::string& outputDir);
	void execContinuousNetworkImaging(TFloat valueScale, ArrayProcessor<XYZValueFactor<TFloat>>& processor,
						unsigned int baseElement, bool coherenceFactorEnabled);
	void saveSignalSequence(const ParameterMap& taskPM, unsigned int baseElement,
				const std::vector<TFloat>& txDelays,
				TnRnAcquisition<TFloat>& acquisition);
	void createImagesFromSavedSignalSequence(const ParameterMap& taskPM,
							unsigned int baseElement, TFloat valueScale, bool coherenceFactorEnabled,
							TnRnAcquisition<TFloat>& acq, ArrayProcessor<XYZValueFactor<TFloat>>& processor);

	Project& project_;
	Matrix<XYZValueFactor<TFloat>> gridData_;
	std::vector<XYZ<float>> pointList_;
};

} // namespace Lab

#endif // SINGLEVIRTUALSOURCEMETHOD_H
