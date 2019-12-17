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

template<typename TFloat, unsigned int UpsamplingFactor>
class LinearInterpolator {
public:
	static void interpolate(const TFloat* input, std::size_t inputLength, TFloat* output);
private:
	LinearInterpolator() = delete;
};

// The argument "ouput" must point to an array of size inputLength * upsamplingFactor_.
template<typename TFloat, unsigned int UpsamplingFactor>
void
LinearInterpolator<TFloat, UpsamplingFactor>::interpolate(const TFloat* input, std::size_t inputLength, TFloat* output)
{
	for (std::size_t i = 0; i < inputLength - 1; ++i) {
		const TFloat step = (input[i + 1] - input[i]) / UpsamplingFactor;
		TFloat* destPtr = output + i * UpsamplingFactor;
		TFloat value = input[i];
		for (unsigned int j = 0; j < UpsamplingFactor; ++j) {
			*destPtr++ = value;
			value += step;
		}
	}

	// The last interval.
	const std::size_t i = inputLength - 1;
	const TFloat step = -input[i] / UpsamplingFactor;
	TFloat* destPtr = output + i * UpsamplingFactor;
	TFloat value = input[i];
	for (unsigned int j = 0; j < UpsamplingFactor; ++j) {
		*destPtr++ = value;
		value += step;
	}
}

} // namespace Lab

#endif /* LINEARINTERPOLATOR_H_ */
