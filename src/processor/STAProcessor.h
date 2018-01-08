#ifndef STAPROCESSOR_H_
#define STAPROCESSOR_H_

#include "Matrix2.h"
#include "XZValueFactor.h"



namespace Lab {

template<typename FloatType>
class STAProcessor {
public:
	STAProcessor() {}
	virtual ~STAProcessor() {}

	virtual void process(unsigned int baseElement, Matrix2<XZValueFactor<FloatType> >& gridData) = 0;
};

} // namespace Lab

#endif /* STAPROCESSOR_H_ */
