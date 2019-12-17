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

template<typename TFloat>
class SavedTnRnAcquisition : public TnRnAcquisition<TFloat> {
public:
	SavedTnRnAcquisition(const Project& project,
				unsigned int numRxElements,
				const std::string& dataDir);
	virtual ~SavedTnRnAcquisition() = default;

	void setDataDir(const std::string& dataDir) { dataDir_ = dataDir; }

	virtual void prepare(unsigned int baseElement, const std::vector<TFloat>& txDelays);
	virtual void execute(typename TnRnAcquisition<TFloat>::AcquisitionDataType& acqData);
private:
	SavedTnRnAcquisition(const SavedTnRnAcquisition&) = delete;
	SavedTnRnAcquisition& operator=(const SavedTnRnAcquisition&) = delete;
	SavedTnRnAcquisition(SavedTnRnAcquisition&&) = delete;
	SavedTnRnAcquisition& operator=(SavedTnRnAcquisition&&) = delete;

	const Project& project_;
	const unsigned int numRxElements_;
	unsigned int baseElement_;
	std::string dataDir_;
};



template<typename TFloat>
SavedTnRnAcquisition<TFloat>::SavedTnRnAcquisition(
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
SavedTnRnAcquisition<TFloat>::prepare(unsigned int baseElement, const std::vector<TFloat>& /*txDelays*/)
{
	baseElement_ = baseElement;
}

template<typename TFloat>
void
SavedTnRnAcquisition<TFloat>::execute(typename TnRnAcquisition<TFloat>::AcquisitionDataType& acqData)
{
	std::string filePath = FileUtil::signalsPath(dataDir_, baseElement_);
	project_.loadHDF5(filePath, "signal", acqData);

	if (numRxElements_ != acqData.n1()) {
		THROW_EXCEPTION(InvalidFileException, "Invalid number of rx elements (" << acqData.n1() <<
			") in file " << filePath << ", should be " << numRxElements_ << '.');
	}
}

} // namespace Lab

#endif // SAVEDTNRNACQUISITION_H
