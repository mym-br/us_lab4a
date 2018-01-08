#ifndef XY_H_
#define XY_H_

#include <iostream>



namespace Lab {

template<typename T>
struct XY {
	T x;
	T y;
};

template<typename T>
std::ostream&
operator<<(std::ostream& out, const XY<T>& data)
{
	out << data.x << ' ' << data.y;
	return out;
}

} // namespace Lab

#endif /* XY_H_ */
