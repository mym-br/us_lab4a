/***************************************************************************
 *  Copyright 2014, 2017, 2018 Marcelo Y. Matuda                           *
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

#ifndef LINEARINTERPOLATOR_H_
#define LINEARINTERPOLATOR_H_

#include <cstddef> /* std::size_t */



namespace Lab {

template<typename FloatType, unsigned int UpsamplingFactor>
class LinearInterpolator {
public:
	LinearInterpolator();
	~LinearInterpolator();

	void interpolate(const FloatType* input, std::size_t inputLength, FloatType* output);
private:
	LinearInterpolator(const LinearInterpolator&);
	LinearInterpolator& operator=(const LinearInterpolator&);
};

template<typename FloatType, unsigned int UpsamplingFactor>
LinearInterpolator<FloatType, UpsamplingFactor>::LinearInterpolator()
{
}

template<typename FloatType, unsigned int UpsamplingFactor>
LinearInterpolator<FloatType, UpsamplingFactor>::~LinearInterpolator()
{
}

// The argument "ouput" must point to an array of size inputLength * upsamplingFactor_.
template<typename FloatType, unsigned int UpsamplingFactor>
void
LinearInterpolator<FloatType, UpsamplingFactor>::interpolate(const FloatType* input, std::size_t inputLength, FloatType* output)
{
	for (std::size_t i = 0; i < inputLength - 1; ++i) {
		const FloatType step = (input[i + 1] - input[i]) / UpsamplingFactor;
		FloatType* destPtr = output + i * UpsamplingFactor;
		FloatType value = input[i];
		for (unsigned int j = 0; j < UpsamplingFactor; ++j) {
			*destPtr++ = value;
			value += step;
		}
	}

	// The last interval.
	const std::size_t i = inputLength - 1;
	const FloatType step = -input[i] / UpsamplingFactor;
	FloatType* destPtr = output + i * UpsamplingFactor;
	FloatType value = input[i];
	for (unsigned int j = 0; j < UpsamplingFactor; ++j) {
		*destPtr++ = value;
		value += step;
	}
}

} // namespace Lab

#endif /* LINEARINTERPOLATOR_H_ */
