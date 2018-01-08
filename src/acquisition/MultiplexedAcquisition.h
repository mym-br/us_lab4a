#ifndef MULTIPLEXEDACQUISITION_H_
#define MULTIPLEXEDACQUISITION_H_

#include <vector>

#include "Matrix2.h"



namespace Lab {

template<typename FloatType>
class MultiplexedAcquisition {
public:
	typedef Matrix2<FloatType> AcquisitionDataType;

	MultiplexedAcquisition() {}
	virtual ~MultiplexedAcquisition() {}

	virtual void execute(
			unsigned int baseElement,
			const std::vector<FloatType>& txDelays,
			AcquisitionDataType& acqData) = 0;
};

} // namespace Lab

#endif /* MULTIPLEXEDACQUISITION_H_ */
