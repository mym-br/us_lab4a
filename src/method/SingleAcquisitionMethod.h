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

#ifndef SINGLEACQUISITIONMETHOD_H_
#define SINGLEACQUISITIONMETHOD_H_

#include <string>

#include <boost/scoped_ptr.hpp>

#include "Method.h"



namespace Lab {

class ArrayAcqClient;
class Project;

struct SingleAcquisitionMethodConfiguration {
	double samplingFrequency;
	double centerFrequency;
	double acquisitionTime;
	double minGain;
	//double maxGain;
	unsigned int numElementsMux;
	unsigned int numElements;
	unsigned int baseElement;
	unsigned int txGroupElement;
	unsigned int rxGroupElement;
	unsigned int averageN;
	std::string savedAcqDir;
};

class SingleAcquisitionMethod : public Method {
public:
	SingleAcquisitionMethod(Project& project);
	virtual ~SingleAcquisitionMethod();

	virtual void execute();

private:
	SingleAcquisitionMethod(const SingleAcquisitionMethod&);
	SingleAcquisitionMethod& operator=(const SingleAcquisitionMethod&);

	void fillConfiguration();

	Project& project_;
	SingleAcquisitionMethodConfiguration config_;
	boost::scoped_ptr<ArrayAcqClient> acq_;
};

} // namespace Lab

#endif /* SINGLEACQUISITIONMETHOD_H_ */
