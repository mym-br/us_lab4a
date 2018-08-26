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

#include "HDF5Util.h"

#include <algorithm> /* max_element, min */

#include "Util.h"

#define HDF5UTIL_MAX_CHUNK_SIZE_2D 64
#define HDF5UTIL_MAX_CHUNK_SIZE_1D (64 * 64 * 2)



namespace Lab {
namespace HDF5Util {

template<>
PredType
hdf5MemoryType<std::vector<float>>()
{
	return PredType::NATIVE_FLOAT;
}

template<>
PredType
hdf5MemoryType<std::vector<std::pair<float, float>>>()
{
	return PredType::NATIVE_FLOAT;
}

template<>
PredType
hdf5MemoryType<Matrix<float>>()
{
	return PredType::NATIVE_FLOAT;
}

template<>
PredType
hdf5MemoryType<std::vector<double>>()
{
	return PredType::NATIVE_DOUBLE;
}

template<>
PredType
hdf5MemoryType<std::vector<std::pair<double, double>>>()
{
	return PredType::NATIVE_DOUBLE;
}

template<>
PredType
hdf5MemoryType<Matrix<double>>()
{
	return PredType::NATIVE_DOUBLE;
}

template<>
PredType
hdf5MemoryType<std::vector<unsigned int>>()
{
	return PredType::NATIVE_UINT;
}

void
calcChunkDims(const std::vector<hsize_t>& dataDims, std::vector<hsize_t>& chunkDims)
{
	chunkDims.assign(dataDims.size(), 0);
	unsigned int rank = dataDims.size();
	switch (rank) {
	case 1:
		chunkDims[0] = std::min<hsize_t>(dataDims[0], HDF5UTIL_MAX_CHUNK_SIZE_1D);
		return;
	case 2:
		if (dataDims[0] == 1) {
			chunkDims[0] = 1;
			chunkDims[1] = std::min<hsize_t>(dataDims[1], HDF5UTIL_MAX_CHUNK_SIZE_1D);
			return;
		} else if (dataDims[1] == 1) {
			chunkDims[0] = std::min<hsize_t>(dataDims[0], HDF5UTIL_MAX_CHUNK_SIZE_1D);
			chunkDims[1] = 1;
			return;
		} else if (dataDims[0] >= HDF5UTIL_MAX_CHUNK_SIZE_2D && dataDims[1] >= HDF5UTIL_MAX_CHUNK_SIZE_2D) {
			chunkDims[0] = std::min<hsize_t>(dataDims[0], HDF5UTIL_MAX_CHUNK_SIZE_2D);
			chunkDims[1] = std::min<hsize_t>(dataDims[1], HDF5UTIL_MAX_CHUNK_SIZE_2D);
			return;
		} else {
			chunkDims = dataDims;
			hsize_t size = Util::multiplyElements(chunkDims);
			while (size > HDF5UTIL_MAX_CHUNK_SIZE_1D) {
				std::vector<hsize_t>::iterator maxElement = std::max_element(chunkDims.begin(), chunkDims.end());
				*maxElement /= 2;
				size = Util::multiplyElements(chunkDims);
			}
			return;
		}
	default:
		THROW_EXCEPTION(InvalidValueException, "Invalid rank: " << rank << '.');
	}
}

} // namespace HDF5Util
} // namespace Lab
