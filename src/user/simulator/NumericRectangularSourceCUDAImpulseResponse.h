/***************************************************************************
 *  Copyright 2020 Marcelo Y. Matuda                                       *
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
#ifndef NUMERIC_RECTANGULAR_SOURCE_CUDA_IMPULSE_RESPONSE_H
#define NUMERIC_RECTANGULAR_SOURCE_CUDA_IMPULSE_RESPONSE_H

#include <cstddef> /* std::size_t */
#include <memory>
#include <vector>



namespace Lab {

// Calculate the acoustic field generated by a flat rectangular surface,
// using the numeric solution provided by:
//
// Piwakowski, B.
// Delannoy, B.
// Method for computing spatial pulse response: Time-domain approach
// J. Acoust. Soc. Am., vol. 86, no. 6, pp. 2422-2432, 1989.
// DOI: 10.1121/1.398449
//
// See also:
// Lasota, H.
// Salamon, R.
// Delannoy, B.
// Acoustic diffraction analysis by the impulse response method: A line impulse response approach},
// J. Acoust. Soc. Am., vol. 76, no. 1, pp. 280-290, 1984.
// DOI: 10.1121/1.391115
//
// Note:
// - The source is surrounded by a rigid baffle.
class NumericRectangularSourceCUDAImpulseResponse {
public:
	NumericRectangularSourceCUDAImpulseResponse(
					float samplingFreq,
					float propagationSpeed,
					float sourceWidth,
					float sourceHeight,
					float subElemSize);
	~NumericRectangularSourceCUDAImpulseResponse();

	// Return h/c.
	void getImpulseResponse(float x, float y, float z,
				std::size_t& hOffset /* samples */, std::vector<float>& h);
private:
	struct CUDAData;

	float samplingFreq_;
	float propagationSpeed_;
	float subElemWidth_;
	float subElemHeight_;
	unsigned int numSubElem_;
	std::unique_ptr<CUDAData> data_;
};

} // namespace Lab

#endif // NUMERIC_RECTANGULAR_SOURCE_CUDA_IMPULSE_RESPONSE_H
