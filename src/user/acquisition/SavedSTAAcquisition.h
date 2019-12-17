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

#include <string>

#include "Exception.h"
#include "FileUtil.h"
#include "Project.h"
#include "STAAcquisition.h"



namespace Lab {

template<typename TFloat>
class SavedSTAAcquisition : public STAAcquisition<TFloat> {
public:
	SavedSTAAcquisition(const Project& project,
				unsigned int numRxElements,
				const std::string& dataDir);
	virtual ~SavedSTAAcquisition() = default;

	void setDataDir(const std::string& dataDir) { dataDir_ = dataDir; }

	virtual void prepare(unsigned int baseElement);
	virtual void execute(unsigned int txElement,
				typename STAAcquisition<TFloat>::AcquisitionDataType& acqData);
private:
	SavedSTAAcquisition(const SavedSTAAcquisition&) = delete;
	SavedSTAAcquisition& operator=(const SavedSTAAcquisition&) = delete;
	SavedSTAAcquisition(SavedSTAAcquisition&&) = delete;
	SavedSTAAcquisition& operator=(SavedSTAAcquisition&&) = delete;

	const Project& project_;
	const unsigned int numRxElements_;
	unsigned int baseElement_;
	std::string dataDir_;
};



template<typename TFloat>
SavedSTAAcquisition<TFloat>::SavedSTAAcquisition(
			const Project& project,
			unsigned int numRxElements,
			const std::string& dataDir)
		: project_(project)
		, numRxElements_(numRxElements)
		, baseElement_()
		, dataDir_(dataDir)
{
}

template<typename TFloat>
void
SavedSTAAcquisition<TFloat>::prepare(unsigned int baseElement)
{
	baseElement_ = baseElement;
}

template<typename TFloat>
void
SavedSTAAcquisition<TFloat>::execute(unsigned int txElement,
						typename STAAcquisition<TFloat>::AcquisitionDataType& acqData)
{
	std::string filePath = FileUtil::txElemSignalsPath(dataDir_, baseElement_, txElement);
	project_.loadHDF5(filePath, "signal", acqData);

	if (numRxElements_ != acqData.n1()) {
		THROW_EXCEPTION(InvalidFileException, "Invalid number of rx elements (" << acqData.n1() <<
			") in file " << filePath << ", should be " << numRxElements_ << '.');
	}
}

} // namespace Lab

#endif /* SAVEDSTAACQUISITION_H_ */
