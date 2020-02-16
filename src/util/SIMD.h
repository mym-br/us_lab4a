/***************************************************************************
 *  Copyright 2014, 2017, 2018, 2020 Marcelo Y. Matuda                     *
 *                                                                         *
 *  This program is free software: you can redistribute it and/or modify   *
 *  it under the terms of the GNU General Public License as published by   *
 *  the Free Software Foundation, either version 3 of the License, or      *
 *  (at your option) any later version.                                    *
 *                                                                         *
 *  This program is distributed in the hope that it will be useful,        *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of         *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
 *  GNU General Public License for more details.                           *
 *                                                                         *
 *  You should have received a copy of the GNU General Public License      *
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.  *
 ***************************************************************************/

#ifndef LAB_SIMD_H
#define LAB_SIMD_H

#include <x86intrin.h>



namespace Lab {
namespace SIMD {

float calcDistance(float dx, float dz);
double calcDistance(double dx, double dz);
float calcDistance(float x1, float z1, float x2, float z2);
double calcDistance(double x1, double z1, double x2, double z2);
float calcTwoMediumTravelTime(float x1, float z1, float xi, float zi, float x2, float z2, float invC1, float invC2);
double calcTwoMediumTravelTime(double x1, double z1, double xi, double zi, double x2, double z2, double invC1, double invC2);



// std::sqrt(dx * dx + dz * dz)
inline
float
calcDistance(float dx, float dz)
{
	__m128 v1 = _mm_set_ps(0, 0, dz, dx);
	v1 = _mm_mul_ps(v1, v1);

	__m128 dummy;
	v1 = _mm_hadd_ps(v1, dummy); // v1[0] = v1[0] + v1[1]; v1[1] = v1[2] + v1[3]; v1[2] = dummy[0] + dummy[1]; v1[3] = dummy[2] + dummy[3]
	v1 = _mm_sqrt_ss(v1);

	float f;
	_mm_store_ss(&f, v1);
	return f;
}

// std::sqrt(dx * dx + dz * dz)
inline
double
calcDistance(double dx, double dz)
{
	// dx and dz are received in xmm registers. When using _mm_load_pd,
	// instead of using them directly, gcc stores them on the stack.
	// When using the built-in method, or _mm_set_pd, the values are used directly.
	// ##### VERY SLOW
	//const double d1[2] __attribute__ ((aligned (16))) = {dx, dz};
	//__m128d v1 = _mm_load_pd(d1);

	__m128d v1 = _mm_set_pd(dz, dx);
	v1 = _mm_mul_pd(v1, v1);

//	double d2[2] __attribute__ ((aligned (16)));
//	_mm_store_pd(d2, v1);
//	__m128d v2 = _mm_load_sd(&d2[1]);
//	v1 = _mm_add_pd(v1, v2);

	__m128d dummy;
	v1 = _mm_hadd_pd(v1, dummy); // v1[0] = v1[0] + v1[1]; v1[1] = dummy[0] + dummy[1]
	//v1 = _mm_sqrt_pd(v1);
	v1 = _mm_sqrt_sd(dummy, v1);

	double d;
	_mm_store_sd(&d, v1);
	return d;
}

// std::sqrt((x2 - x1)^2 + (z2 - z1)^2)
inline
float
calcDistance(float x1, float z1, float x2, float z2)
{
	__m128 v1 = _mm_set_ps(0, 0, z2, x2);
	__m128 v2 = _mm_set_ps(0, 0, z1, x1);
	v1 = _mm_sub_ps(v1, v2); // v1 - v2
	v1 = _mm_mul_ps(v1, v1);

	__m128 dummy;
	v1 = _mm_hadd_ps(v1, dummy); // v1[0] = v1[0] + v1[1]; v1[1] = v1[2] + v1[3]; v1[2] = dummy[0] + dummy[1]; v1[3] = dummy[2] + dummy[3]
	v1 = _mm_sqrt_ss(v1);

	float f;
	_mm_store_ss(&f, v1);
	return f;
}

// std::sqrt((x2 - x1)^2 + (z2 - z1)^2)
inline
double
calcDistance(double x1, double z1, double x2, double z2)
{
	__m128d v1 = _mm_set_pd(z2, x2);
	__m128d v2 = _mm_set_pd(z1, x1);
	v1 = _mm_sub_pd(v1, v2); // v1 - v2
	v1 = _mm_mul_pd(v1, v1);

	__m128d dummy;
	v1 = _mm_hadd_pd(v1, dummy); // v1[0] = v1[0] + v1[1]; v1[1] = dummy[0] + dummy[1];
	v1 = _mm_sqrt_sd(dummy, v1);

	double d;
	_mm_store_sd(&d, v1);
	return d;
}

// sqrt((xi - x1)^2 + (zi - z1)^2) * invC1 + sqrt((x2 - xi)^2 + (z2 - zi)^2) * invC2
inline
float
calcTwoMediumTravelTime(float x1, float z1, float xi, float zi, float x2, float z2, float invC1, float invC2)
{
	__m128 v1 = _mm_set_ps(zi, xi, z1, x1);
	__m128 v2 = _mm_set_ps(z2, x2, zi, xi);
	v1 = _mm_sub_ps(v2, v1); // v2 - v1
	v1 = _mm_mul_ps(v1, v1);

	__m128 dummy;
	v1 = _mm_hadd_ps(v1, dummy); // v1[0] = v1[0] + v1[1]; v1[1] = v1[2] + v1[3]; v1[2] = dummy[0] + dummy[1]; v1[3] = dummy[2] + dummy[3]
	v1 = _mm_sqrt_ps(v1);

	v2 = _mm_set_ps(0, 0, invC2, invC1);
	v1 = _mm_mul_ps(v1, v2);
	v1 = _mm_hadd_ps(v1, dummy); // v1[0] = v1[0] + v1[1]; v1[1] = v1[2] + v1[3]; v1[2] = dummy[0] + dummy[1]; v1[3] = dummy[2] + dummy[3]

	float f;
	_mm_store_ss(&f, v1);
	return f;
}

// sqrt((xi - x1)^2 + (zi - z1)^2) * invC1 + sqrt((x2 - xi)^2 + (z2 - zi)^2) * invC2
inline
double
calcTwoMediumTravelTime(double x1, double z1, double xi, double zi, double x2, double z2, double invC1, double invC2)
{
	__m128d v1 = _mm_set_pd(z1, x1);
	__m128d v2 = _mm_set_pd(zi, xi);
	v1 = _mm_sub_pd(v2, v1); // v2 - v1
	v1 = _mm_mul_pd(v1, v1);

	__m128d v3 = _mm_set_pd(z2, x2);
	v3 = _mm_sub_pd(v3, v2); // v3 - v2
	v3 = _mm_mul_pd(v3, v3);

	v1 = _mm_hadd_pd(v1, v3); // v1[0] = v1[0] + v1[1]; v1[1] = v3[0] + v3[1];
	v1 = _mm_sqrt_pd(v1);

	v2 = _mm_set_pd(invC2, invC1);
	v1 = _mm_mul_pd(v1, v2);
	__m128d dummy;
	v1 = _mm_hadd_pd(v1, dummy); // v1[0] = v1[0] + v1[1]; v1[1] = dummy[0] + dummy[1];

	double d;
	_mm_store_sd(&d, v1);
	return d;
}

} // namespace SIMD
} // namespace Lab

#endif // LAB_SIMD_H
