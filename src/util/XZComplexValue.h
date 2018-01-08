#ifndef XZCOMPLEXVALUE_H
#define XZCOMPLEXVALUE_H

#include <complex>



namespace Lab {

template<typename FloatType>
struct XZComplexValue {
	typedef std::complex<FloatType> ValueType;

	FloatType x;
	FloatType z;
	std::complex<FloatType> value;
};

} // namespace Lab

#endif // XZCOMPLEXVALUE_H
