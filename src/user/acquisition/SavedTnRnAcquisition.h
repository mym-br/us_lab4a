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

#ifndef SAVEDTNRNACQUISITION_H
#define SAVEDTNRNACQUISITION_H

#include <string>
#include <vector>

#include "Exception.h"
#include "FileUtil.h"
#include "Project.h"
#include "TnRnAcquisition.h"



namespace Lab {

template<typename FloatType>
class SavedTnRnAcquisition : public TnRnAcquisition<FloatType> {
public:
	SavedTnRnAcquisition(const Project& project,
				unsigned int numRxElements,
				const std::string& dataDir);
	virtual ~SavedTnRnAcquisition();

	void setDataDir(const std::string& dataDir) { dataDir_ = dataDir; }
	virtual void execute(unsigned int baseElement, const std::vector<FloatType>& txDelays,
				typename TnRnAcquisition<FloatType>::AcquisitionDataType& acqData);
private:
	SavedTnRnAcquisition(const SavedTnRnAcquisition&) = delete;
	SavedTnRnAcquisition& operator=(const SavedTnRnAcquisition&) = delete;

	const Project& project_;
	const unsigned int numRxElements_;
	std::string dataDir_;
};



template<typename FloatType>
SavedTnRnAcquisition<FloatType>::SavedTnRnAcquisition(
			const Project& project,
			unsigned int numRxElements,
			const std::string& dataDir)
		: project_(project)
		, numRxElements_(numRxElements)
		, dataDir_(dataDir)
{
}

template<typename FloatType>
SavedTnRnAcquisition<FloatType>::~SavedTnRnAcquisition()
{
}

template<typename FloatType>
void
SavedTnRnAcquisition<FloatType>::execute(unsigned int baseElement, const std::vector<FloatType>& /*txDelays*/,
						typename TnRnAcquisition<FloatType>::AcquisitionDataType& acqData)
{
	std::string filePath = FileUtil::signalsPath(dataDir_, baseElement);
	project_.loadHDF5(filePath, "signal", acqData);

	if (numRxElements_ != acqData.n1()) {
		THROW_EXCEPTION(InvalidFileException, "Invalid number of rx elements (" << acqData.n1() <<
			") in file " << filePath << ", should be " << numRxElements_ << '.');
	}
}

} // namespace Lab

#endif // SAVEDTNRNACQUISITION_H
