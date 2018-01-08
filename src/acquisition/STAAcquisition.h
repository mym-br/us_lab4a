#ifndef STAACQUISITION_H_
#define STAACQUISITION_H_

#include "Matrix2.h"



namespace Lab {

template<typename FloatType>
class STAAcquisition {
public:
	typedef Matrix2<FloatType> AcquisitionDataType;

	STAAcquisition() {}
	virtual ~STAAcquisition() {}

	virtual void execute(unsigned int baseElement, unsigned int txElement /* relative to baseElement */, AcquisitionDataType& acqData) = 0;
};

} // namespace Lab

#endif /* ACQUISITION_H_ */
