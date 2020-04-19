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
#ifndef NUMERIC_CIRCULAR_SOURCE_CUDA_IMPULSE_RESPONSE_H
#define NUMERIC_CIRCULAR_SOURCE_CUDA_IMPULSE_RESPONSE_H

#include <cstddef> /* std::size_t */
#include <memory>
#include <vector>



namespace Lab {

class NumericCircularSourceCUDAImpulseResponse {
public:
	NumericCircularSourceCUDAImpulseResponse(
					float samplingFreq,
					float propagationSpeed,
					float sourceRadius,
					float numSubElemInRadius);
	~NumericCircularSourceCUDAImpulseResponse();

	// Return h/c.
	void getImpulseResponse(float x, float y, float z,
				std::size_t& hOffset /* samples */, std::vector<float>& h);
private:
	struct CUDAData;

	NumericCircularSourceCUDAImpulseResponse(const NumericCircularSourceCUDAImpulseResponse&) = delete;
	NumericCircularSourceCUDAImpulseResponse& operator=(const NumericCircularSourceCUDAImpulseResponse&) = delete;
	NumericCircularSourceCUDAImpulseResponse(NumericCircularSourceCUDAImpulseResponse&&) = delete;
	NumericCircularSourceCUDAImpulseResponse& operator=(NumericCircularSourceCUDAImpulseResponse&&) = delete;

	float samplingFreq_;
	float propagationSpeed_;
	float subElemArea_;
	unsigned int numSubElem_;
	std::unique_ptr<CUDAData> data_;
};

} // namespace Lab

#endif // NUMERIC_CIRCULAR_SOURCE_CUDA_IMPULSE_RESPONSE_H
