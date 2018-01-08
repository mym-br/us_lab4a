#ifndef XZVALUEFACTOR_H_
#define XZVALUEFACTOR_H_

namespace Lab {

template<typename FloatType>
struct XZValueFactor {
	typedef FloatType ValueType;

	FloatType x;
	FloatType z;
	FloatType value;
	FloatType factor;
};

} // namespace Lab

#endif /* XZVALUEFACTOR_H_ */
