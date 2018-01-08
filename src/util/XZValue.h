#ifndef XZVALUE_H_
#define XZVALUE_H_

namespace Lab {

template<typename FloatType>
struct XZValue {
	typedef FloatType ValueType;

	XZValue() {}
	XZValue(FloatType x, FloatType z, FloatType value) : x(x), z(z), value(value) {}

	FloatType x;
	FloatType z;
	FloatType value;

	template<typename FloatType2> XZValue<FloatType>& operator=(const XZValue<FloatType2>& o) {
		if (&o != this) {
			this->x = o.x;
			this->z = o.z;
			this->value = o.value;
		}
		return *this;
	}
};

} // namespace Lab

#endif /* XZVALUE_H_ */
