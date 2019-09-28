/***************************************************************************
 *  Copyright 2014, 2017, 2018, 2019 Marcelo Y. Matuda                     *
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

#ifndef INTERPOLATOR4X_H_
#define INTERPOLATOR4X_H_

#include <cassert>
#include <cstddef> /* std::size_t */
#include <cstring> /* memcpy */
#include <vector>

#include "FFTWFilter.h"



namespace Lab {

// This class is copy constructible and assignable.
template<typename FloatType>
class Interpolator4X {
public:
	Interpolator4X();
	~Interpolator4X();

	// The argument "ouput" must point to an array of size inputLength * UPSAMPLING_FACTOR.
	void interpolate(const FloatType* input, std::size_t inputLength, FloatType* output);
private:
	enum {
		UPSAMPLING_FACTOR = 4 /* do not change */
	};

	void prepare();

	bool initialized_;
	std::vector<FloatType> lowPassFIRFilter_;
	std::vector<FloatType> inputVector_;
	std::vector<FloatType> outputVector_;
	FFTWFilter<FloatType> filter_;
};



template<typename FloatType>
Interpolator4X<FloatType>::Interpolator4X()
		: initialized_()
{
	// Calculated in Octave using:
	//-----------------------------------------------------------------------------
	// upFactor = 4;
	// nPeriods = 10;
	// beta = 5;
	// n = nPeriods * (2 * upFactor);
	// wd = kaiser(n + 1, beta);
	// w = linspace(-nPeriods, nPeriods, n + 1);
	// s = sinc(w);
	// h = s' .* wd;
	//-----------------------------------------------------------------------------
	// UPSAMPLING_FACTOR must be = 4.
	lowPassFIRFilter_ = {
		-1.43105366198583e-18,
		-1.12987833807950e-03,
		-2.10182991017381e-03,
		-1.90045816956871e-03,
		 3.70668154995274e-18,
		 2.92465237736067e-03,
		 5.01533833417820e-03,
		 4.24961391099716e-03,
		-6.98467610339543e-18,
		-5.92827810890104e-03,
		-9.78460543296953e-03,
		-8.02101956023436e-03,
		 1.12071029784912e-17,
		 1.05983435413046e-02,
		 1.71042670496974e-02,
		 1.37452093502169e-02,
		-1.61736567970827e-17,
		-1.75679820870747e-02,
		-2.79654751657672e-02,
		-2.22057501623689e-02,
		 2.15511119619965e-17,
		 2.78491885383373e-02,
		 4.40234603019062e-02,
		 3.47731193965138e-02,
		-2.69054320758213e-17,
		-4.33954624057215e-02,
		-6.86448380782130e-02,
		-5.43908677412120e-02,
		 3.17527702409652e-17,
		 6.89284613460359e-02,
		 1.10519377295669e-01,
		 8.92265500252782e-02,
		-3.56219808482348e-17,
		-1.20057925470521e-01,
		-2.01755992002028e-01,
		-1.73867065267645e-01,
		 3.81188682834355e-17,
		 2.96354170885177e-01,
		 6.33073094648964e-01,
		 8.99060260635281e-01,
		 1.00000000000000e+00,
		 8.99060260635281e-01,
		 6.33073094648964e-01,
		 2.96354170885177e-01,
		 3.81188682834355e-17,
		-1.73867065267645e-01,
		-2.01755992002028e-01,
		-1.20057925470521e-01,
		-3.56219808482348e-17,
		 8.92265500252782e-02,
		 1.10519377295669e-01,
		 6.89284613460359e-02,
		 3.17527702409652e-17,
		-5.43908677412120e-02,
		-6.86448380782130e-02,
		-4.33954624057215e-02,
		-2.69054320758213e-17,
		 3.47731193965138e-02,
		 4.40234603019062e-02,
		 2.78491885383373e-02,
		 2.15511119619965e-17,
		-2.22057501623689e-02,
		-2.79654751657672e-02,
		-1.75679820870747e-02,
		-1.61736567970827e-17,
		 1.37452093502169e-02,
		 1.71042670496974e-02,
		 1.05983435413046e-02,
		 1.12071029784912e-17,
		-8.02101956023436e-03,
		-9.78460543296953e-03,
		-5.92827810890104e-03,
		-6.98467610339543e-18,
		 4.24961391099716e-03,
		 5.01533833417820e-03,
		 2.92465237736067e-03,
		 3.70668154995274e-18,
		-1.90045816956871e-03,
		-2.10182991017381e-03,
		-1.12987833807950e-03,
		-1.43105366198583e-18
	};
}

template<typename FloatType>
Interpolator4X<FloatType>::~Interpolator4X()
{
}

template<typename FloatType>
void
Interpolator4X<FloatType>::prepare()
{
	assert(!initialized_);

	filter_.setCoefficients(lowPassFIRFilter_);

	initialized_ = true;
}

template<typename FloatType>
void
Interpolator4X<FloatType>::interpolate(const FloatType* input, std::size_t inputLength, FloatType* output)
{
	if (!initialized_) prepare();

	// Upsample.
	inputVector_.resize(inputLength * UPSAMPLING_FACTOR);
	for (std::size_t i = 0, j = 0; i < inputLength; ++i, j += UPSAMPLING_FACTOR) {
		inputVector_[j] = input[i];
		for (std::size_t k = j + 1, end = j + UPSAMPLING_FACTOR; k < end; ++k) {
			inputVector_[k] = 0;
		}
	}

	// Apply the anti-aliasing filter.
	filter_.filter(inputVector_, outputVector_);

	// Copy to the output, compensating for the FIR filter delay. The signal is truncated at both ends.
	const std::size_t offset = (lowPassFIRFilter_.size() - 1) / 2;
//	for (std::size_t i = 0; i < inputLength * UPSAMPLING_FACTOR; ++i) {
//		output[i] = outputVector_[i + offset];
//	}
	memcpy(output, &outputVector_[offset], inputLength * UPSAMPLING_FACTOR * sizeof(FloatType));
}

} // namespace Lab

#endif /* INTERPOLATOR4X_H_ */
