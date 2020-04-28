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

inline
void
atomicAddGlobalFloat(volatile __global float* p, float value)
{
	union {
		unsigned int ui;
		float f;
	} old, /*expected,*/ next;

	// https://streamhpc.com/blog/2016-02-09/atomic-operations-for-floats-in-opencl-improved/
	// Can be very slow!
//	old.f = *p;
//	do {
//		expected.f = old.f;
//		next.f = expected.f + value;
//		old.ui = atomic_cmpxchg(
//				(volatile __global unsigned int*) p,
//				expected.ui,
//				next.ui);
//	} while (old.ui != expected.ui);

	// https://community.khronos.org/t/roadmap-for-atomic-floating-point-support/7619
	// Slower than atomicAdd in CUDA.
	old.f = value;
	do {
		next.ui = atomic_xchg((volatile __global unsigned int*) p, 0U);
		next.f += old.f;
		old.ui = atomic_xchg((volatile __global unsigned int*) p, next.ui);
	} while (old.ui != 0U);
}

inline
void
atomicAddLocalFloat(volatile __local float* p, float value)
{
	union {
		unsigned int ui;
		float f;
	} old, next;

	old.f = value;
	do {
		next.ui = atomic_xchg((volatile __local unsigned int*) p, 0U);
		next.f += old.f;
		old.ui = atomic_xchg((volatile __local unsigned int*) p, next.ui);
	} while (old.ui != 0U);
}

#pragma OPENCL EXTENSION cl_khr_int64_base_atomics: enable

inline
void
atomicAddGlobalDouble(volatile __global double* p, double value)
{
	union {
		unsigned long ul;
		double d;
	} old, next;

	old.d = value;
	do {
		next.ul = atomic_xchg((volatile __global unsigned long*) p, 0UL);
		next.d += old.d;
		old.ul = atomic_xchg((volatile __global unsigned long*) p, next.ul);
	} while (old.ul != 0UL);
}

inline
void
atomicAddLocalDouble(volatile __local double* p, double value)
{
	union {
		unsigned long ul;
		double d;
	} old, next;

	old.d = value;
	do {
		next.ul = atomic_xchg((volatile __local unsigned long*) p, 0UL);
		next.d += old.d;
		old.ul = atomic_xchg((volatile __local unsigned long*) p, next.ul);
	} while (old.ul != 0UL);
}

)CLC";
}

} // namespace OCLAtomic
} // namespace Lab

#endif // OCL_ATOMIC_H
