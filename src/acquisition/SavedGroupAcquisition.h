#ifndef SAVEDGROUPACQUISITION_H_
#define SAVEDGROUPACQUISITION_H_

#include <sstream>
#include <string>

#include "Exception.h"
#include "Matrix2.h"
#include "Project.h"



namespace Lab {

template<typename FloatType>
class SavedGroupAcquisition {
public:
	typedef Matrix2<FloatType> AcquisitionDataType;

	SavedGroupAcquisition(
		const Project& project,
		unsigned int numGroupElements,
		const std::string& savedDataDir);
	virtual ~SavedGroupAcquisition();

	void execute(unsigned int acquisitionNumber, AcquisitionDataType& acqData);
private:
	SavedGroupAcquisition(const SavedGroupAcquisition&);
	SavedGroupAcquisition& operator=(const SavedGroupAcquisition&);

	const Project& project_;
	const unsigned int numGroupElements_;
	std::string savedDataDir_;
};



template<typename FloatType>
SavedGroupAcquisition<FloatType>::SavedGroupAcquisition(
			const Project& project,
			unsigned int numGroupElements,
			const std::string& savedDataDir)
		: project_(project)
		, numGroupElements_(numGroupElements)
		, savedDataDir_(savedDataDir)
{
}

template<typename FloatType>
SavedGroupAcquisition<FloatType>::~SavedGroupAcquisition()
{
}

template<typename FloatType>
void
SavedGroupAcquisition<FloatType>::execute(unsigned int acquisitionNumber, AcquisitionDataType& acqData)
{
	std::ostringstream filePath;
	filePath << savedDataDir_ << "/signal" << acquisitionNumber;
	project_.loadHDF5(filePath.str(), "signal", acqData);

	if (numGroupElements_ != acqData.n1()) {
		THROW_EXCEPTION(InvalidFileException, "Invalid number of rx elements (" << acqData.n1() <<
			") in file " << filePath.str() << ", should be " << numGroupElements_ << '.');
	}
}

} // namespace Lab

#endif /* SAVEDGROUPACQUISITION_H_ */
