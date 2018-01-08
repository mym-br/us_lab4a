#ifndef SECTORIALSCANACQUISITION_H
#define SECTORIALSCANACQUISITION_H

#include "Matrix2.h"
#include "XZValue.h"



namespace Lab {

template<typename FloatType>
class SectorialScanAcquisition {
public:
	typedef Matrix2<XZValue<FloatType>> AcquisitionDataType;

	SectorialScanAcquisition() {}
	virtual ~SectorialScanAcquisition() {}

	virtual void execute(AcquisitionDataType& acqData) = 0;
};

} // namespace Lab

#endif // SECTORIALSCANACQUISITION_H
