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

#ifndef SAVEDSTAACQUISITION_H_
#define SAVEDSTAACQUISITION_H_

#include <iomanip>
#include <sstream>
#include <string>

#include "Exception.h"
#include "ParameterMap.h"
#include "Project.h"
#include "STAAcquisition.h"



namespace Lab {

template<typename FloatType>
class SavedSTAAcquisition : public STAAcquisition<FloatType> {
public:
	SavedSTAAcquisition(
		const Project& project,
		unsigned int numElements);
	virtual ~SavedSTAAcquisition();

	virtual void execute(unsigned int baseElement, unsigned int txElement,
				typename STAAcquisition<FloatType>::AcquisitionDataType& acqData);
private:
	SavedSTAAcquisition(const SavedSTAAcquisition&);
	SavedSTAAcquisition& operator=(const SavedSTAAcquisition&);

	const Project& project_;
	const unsigned int numElements_;
	std::string dataDir_;
};



template<typename FloatType>
SavedSTAAcquisition<FloatType>::SavedSTAAcquisition(
			const Project& project,
			unsigned int numElements)
		: project_(project)
		, numElements_(numElements)
{
	ConstParameterMapPtr pm = project_.loadParameterMap("config-saved_sta_acquisition.txt");
	dataDir_ = pm->value<std::string>("data_dir");
}

template<typename FloatType>
SavedSTAAcquisition<FloatType>::~SavedSTAAcquisition()
{
}

template<typename FloatType>
void
SavedSTAAcquisition<FloatType>::execute(unsigned int baseElement, unsigned int txElement,
						typename STAAcquisition<FloatType>::AcquisitionDataType& acqData)
{
	if (txElement >= numElements_) {
		THROW_EXCEPTION(InvalidParameterException, "Invalid tx element (" << txElement <<"), should be < " << numElements_ << '.');
	}

	std::ostringstream filePath;
	filePath << dataDir_ << std::setfill('0') << "/signal-base" << std::setw(4) << baseElement << "-tx" << std::setw(4) << txElement;
	project_.loadHDF5(filePath.str(), "signal", acqData);

	if (numElements_ != acqData.n1()) {
		THROW_EXCEPTION(InvalidFileException, "Invalid number of rx elements (" << acqData.n1() <<
			") in file " << filePath.str() << ", should be " << numElements_ << '.');
	}
}

} // namespace Lab

#endif /* SAVEDSTAACQUISITION_H_ */
