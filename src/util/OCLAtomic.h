// This file is in the public domain.

#ifndef OCL_ATOMIC_H
#define OCL_ATOMIC_H

#include <string>

namespace Lab {
namespace OCLAtomic {

inline
std::string
code() {
	return R"CLC(

// https://streamhpc.com/blog/2016-02-09/atomic-operations-for-floats-in-opencl-improved/
// Can be very slow!
inline
float
atomicAddGlobalFloat(volatile __global float* p, float value)
{
	union {
		unsigned int ui;
		float f;
	} old, expected, next;

	old.f = *p;
	do {
		expected.f = old.f;
		next.f = expected.f + value;
		old.ui = atomic_cmpxchg(
				(volatile __global unsigned int*) p,
				expected.ui,
				next.ui);
	} while (old.ui != expected.ui);
}

)CLC";
}

} // namespace OCLAtomic
} // namespace Lab

#endif // OCL_ATOMIC_H
