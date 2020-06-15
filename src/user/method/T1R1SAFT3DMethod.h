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

#ifndef T1R1SAFT3DMETHOD_H
#define T1R1SAFT3DMETHOD_H

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

template<typename TFloat>
class T1R1SAFT3DMethod : public Method {
public:
	T1R1SAFT3DMethod(Project& project);
	virtual ~T1R1SAFT3DMethod() = default;

	virtual void execute();
private:
	T1R1SAFT3DMethod(const T1R1SAFT3DMethod&) = delete;
	T1R1SAFT3DMethod& operator=(const T1R1SAFT3DMethod&) = delete;
	T1R1SAFT3DMethod(T1R1SAFT3DMethod&&) = delete;
	T1R1SAFT3DMethod& operator=(T1R1SAFT3DMethod&&) = delete;

	void process(TFloat valueScale, ArrayProcessor<XYZValueFactor<TFloat>>& processor, unsigned int baseElement, const std::string& outputDir);
	void useCoherenceFactor(TFloat valueScale, const std::string& outputDir);

	Project& project_;
	Matrix<XYZValueFactor<TFloat>> gridData_;
	std::vector<XYZ<float>> pointList_;
	Visualization::Value visual_;
};

} // namespace Lab

#endif // T1R1SAFT3DMETHOD_H
