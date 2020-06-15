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

#ifndef STAMETHOD_H_
#define STAMETHOD_H_

#include <string>
#include <vector>

#include "Matrix.h"
#include "Method.h"
#include "Visualization.h"
#include "XYZ.h"
#include "XYZValueFactor.h"



namespace Lab {

template<typename TPoint> class ArrayProcessor;
class Project;
template<typename TFloat> struct STAConfiguration;

template<typename TFloat>
class STAMethod : public Method {
public:
	STAMethod(Project& project);
	virtual ~STAMethod() = default;

	virtual void execute();
private:
	STAMethod(const STAMethod&) = delete;
	STAMethod& operator=(const STAMethod&) = delete;
	STAMethod(STAMethod&&) = delete;
	STAMethod& operator=(STAMethod&&) = delete;

	void process(TFloat valueScale, ArrayProcessor<XYZValueFactor<TFloat>>& processor, unsigned int baseElement,
			const std::string& outputDir);
	template<typename P>
		void process(TFloat valueScale, P& processor, unsigned int baseElement,
				const STAConfiguration<TFloat>& config, const std::string& outputDir);
	void useCoherenceFactor(TFloat valueScale, bool calculateEnvelope, const std::string& outputDir);

	Project& project_;
	Matrix<XYZValueFactor<TFloat>> gridData_;
	std::vector<XYZ<float>> pointList_;
	Visualization::Value visual_;
};

} // namespace Lab

#endif /* STAMETHOD_H_ */
