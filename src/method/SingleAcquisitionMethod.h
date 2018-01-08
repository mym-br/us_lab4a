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
