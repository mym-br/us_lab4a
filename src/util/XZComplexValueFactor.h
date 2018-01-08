#ifndef XZCOMPLEXVALUEFACTOR_H
#define XZCOMPLEXVALUEFACTOR_H

#include <complex>



namespace Lab {

template<typename FloatType>
struct XZComplexValueFactor {
	typedef std::complex<FloatType> ValueType;

	FloatType x;
	FloatType z;
	std::complex<FloatType> value;
	FloatType factor;
};

} // namespace Lab

#endif // XZCOMPLEXVALUEFACTOR_H
