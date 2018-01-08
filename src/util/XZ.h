#ifndef XZ_H_
#define XZ_H_

#include <iostream>



namespace Lab {

template<typename FloatType>
struct XZ {
	typedef FloatType ValueType;

	XZ() {}
	XZ(FloatType x, FloatType z) : x(x), z(z) {}

	template<typename FloatType2> XZ<FloatType>& operator=(const XZ<FloatType2>& o) {
		if (&o != this) {
			this->x = o.x;
			this->z = o.z;
		}
		return *this;
	}

	FloatType x;
	FloatType z;
};

template<typename T>
std::ostream&
operator<<(std::ostream& out, const XZ<T>& data)
{
	out << data.x << ' ' << data.z;
	return out;
}

} // namespace Lab

#endif /* XZ_H_ */
