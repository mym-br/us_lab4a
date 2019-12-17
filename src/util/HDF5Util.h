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

#ifndef HDF5UTIL_H_
#define HDF5UTIL_H_

#include <cstddef> /* std::size_t */
#include <string>
#include <utility>
#include <vector>

#include <H5Cpp.h>

#include "Exception.h"
#include "lzf_filter.h"
#include "Matrix.h"

#define HDF5UTIL_USE_LZF 1

#ifdef HDF5UTIL_USE_LZF
# define HDF5UTIL_FILTER H5PY_FILTER_LZF
#else
# define HDF5UTIL_FILTER H5Z_FILTER_DEFLATE
#endif

namespace Lab {
namespace HDF5Util {

constexpr const char* fileSuffix = ".h5";
constexpr int deflateLevel = 1;

#ifndef H5_NO_NAMESPACE
using namespace H5;
#endif

template<typename TFloat> void resize(std::vector<TFloat>& container, hsize_t n1, hsize_t n2);
template<typename TFloat> void resize(std::vector<std::pair<TFloat, TFloat>>& container, hsize_t n1, hsize_t n2);
template<typename TFloat> void resize(Matrix<TFloat>& container, hsize_t n1, hsize_t n2);

template<typename TFloat> TFloat* getBeginPtr(std::vector<TFloat>& container);
template<typename TFloat> const TFloat* getBeginPtr(const std::vector<TFloat>& container);

template<typename TFloat> TFloat* getBeginPtr(std::vector<std::pair<TFloat, TFloat>>& container);
template<typename TFloat> const TFloat* getBeginPtr(const std::vector<std::pair<TFloat, TFloat>>& container);

template<typename TFloat> TFloat* getBeginPtr(Matrix<TFloat>& container);
template<typename TFloat> const TFloat* getBeginPtr(const Matrix<TFloat>& container);

template<typename TFloat> void getSize(const std::vector<TFloat>& container, hsize_t& n1, hsize_t& n2);
template<typename TFloat> void getSize(const std::vector<std::pair<TFloat, TFloat>>& container, hsize_t& n1, hsize_t& n2);
template<typename TFloat> void getSize(const Matrix<TFloat>& container, hsize_t& n1, hsize_t& n2);

template<typename T> void load(const std::string& filePath, const std::string& dataSetName, T& container);
template<typename T> void save(const T& container, const std::string& filePath, const std::string& datasetName);

template<typename T> PredType hdf5MemoryType();

void calcChunkDims(const std::vector<hsize_t>& dataDims, std::vector<hsize_t>& chunkDims);



template<typename TFloat>
void
resize(std::vector<TFloat>& container, hsize_t n1, hsize_t n2)
{
	if (n1 != 1) {
		THROW_EXCEPTION(InvalidParameterException, "The first dimension (" << n1 << ") is not equal to 1.");
	}
	if (n2 == 0) {
		THROW_EXCEPTION(InvalidParameterException, "The second dimension (" << n2 << ") is equal to 0.");
	}
	container.resize(n2);
}

template<typename TFloat>
void
resize(std::vector<std::pair<TFloat, TFloat>>& container, hsize_t n1, hsize_t n2)
{
	if (n1 == 0) {
		THROW_EXCEPTION(InvalidParameterException, "The first dimension (" << n1 << ") is equal to 0.");
	}
	if (n2 != 2) {
		THROW_EXCEPTION(InvalidParameterException, "The second dimension (" << n2 << ") is not equal to 2.");
	}
	container.resize(n1);
}

template<typename TFloat>
void
resize(Matrix<TFloat>& container, hsize_t n1, hsize_t n2)
{
	if (n1 == 0) {
		THROW_EXCEPTION(InvalidParameterException, "The first dimension (" << n1 << ") is equal to 0.");
	}
	if (n2 == 0) {
		THROW_EXCEPTION(InvalidParameterException, "The second dimension (" << n2 << ") is equal to 0.");
	}
	container.resize(n1, n2);
}

template<typename TFloat>
TFloat*
getBeginPtr(std::vector<TFloat>& v)
{
	return &v[0];
}

template<typename TFloat>
const TFloat*
getBeginPtr(const std::vector<TFloat>& v)
{
	return &v[0];
}

template<typename TFloat>
TFloat*
getBeginPtr(std::vector<std::pair<TFloat, TFloat>>& v)
{
	return &(v[0].first);
}

template<typename TFloat>
const TFloat*
getBeginPtr(const std::vector<std::pair<TFloat, TFloat>>& v)
{
	return &(v[0].first);
}

template<typename TFloat>
TFloat*
getBeginPtr(Matrix<TFloat>& m)
{
	return &m(0, 0);
}

template<typename TFloat>
const TFloat*
getBeginPtr(const Matrix<TFloat>& m)
{
	return &m(0, 0);
}

template<typename TFloat>
void
getSize(const std::vector<TFloat>& v, hsize_t& n1, hsize_t& n2)
{
	n1 = 1;
	n2 = v.size();
}

template<typename TFloat>
void
getSize(const std::vector<std::pair<TFloat, TFloat>>& v, hsize_t& n1, hsize_t& n2)
{
	n1 = v.size();
	n2 = 2;
}

template<typename TFloat>
void
getSize(const Matrix<TFloat>& m, hsize_t& n1, hsize_t& n2)
{
	n1 = m.n1();
	n2 = m.n2();
}

template<typename T>
void
load(const std::string& filePath, const std::string& dataSetName, T& container)
{
	try {
		H5::Exception::dontPrint();

		H5File file(filePath, H5F_ACC_RDONLY);

		DataSet dataSet = file.openDataSet(dataSetName);
		H5T_class_t typeClass = dataSet.getTypeClass();
		if (typeClass != H5T_FLOAT) {
			THROW_EXCEPTION(InvalidFileException, "The data type class in the file " << filePath << " is not floating point.");
		}

		H5::FloatType type = dataSet.getFloatType();
		//H5std_string orderString;
		//H5T_order_t order = type.getOrder(orderString);
		//std::cout << orderString << std::endl;
		std::size_t typeSize = type.getSize();
		if (typeSize != sizeof(double)) {
			THROW_EXCEPTION(InvalidFileException, "The data type in the file " << filePath << " is not double.");
		}

		DataSpace dataSetSpace = dataSet.getSpace();
		if (!dataSetSpace.isSimple()) {
			THROW_EXCEPTION(InvalidFileException, "The dataset in the file " << filePath << " is not simple.");
		}

		int rank = dataSetSpace.getSimpleExtentNdims();
		switch (rank) {
		case 1:
			{
				hsize_t dim;
				dataSetSpace.getSimpleExtentDims(&dim);
				resize(container, 1, dim);
			}
			break;
		case 2:
			{
				hsize_t dims[2];
				dataSetSpace.getSimpleExtentDims(dims);
				resize(container, dims[0], dims[1]);
			}
			break;
		default:
			THROW_EXCEPTION(InvalidFileException, "The rank in the file " << filePath << " is not 1 or 2.");
		}

		dataSet.read(getBeginPtr(container), hdf5MemoryType<T>());

	} catch (const H5::Exception& e) {
		THROW_EXCEPTION(IOException, "An error ocurred in HDF5 library with file " << filePath << ": " << e.getCDetailMsg());
	}
}

template<typename T>
void
save(const T& container, const std::string& filePath, const std::string& datasetName)
{
	const int rank = 2;
	hsize_t n1, n2;
	getSize(container, n1, n2);
	if (n1 == 0) {
		THROW_EXCEPTION(InvalidParameterException, "The container's first dimension (" << n1 << ") is equal to 0.");
	}
	if (n2 == 0) {
		THROW_EXCEPTION(InvalidParameterException, "The container's second dimension (" << n2 << ") is equal to 0.");
	}

	try {
		H5::Exception::dontPrint();

		H5File file(filePath, H5F_ACC_TRUNC);
		DSetCreatPropList dcpl;

		std::vector<hsize_t> dataDims;
		dataDims.push_back(n1);
		dataDims.push_back(n2);

		// Check if the filter is available.
		htri_t avail = H5Zfilter_avail(HDF5UTIL_FILTER);
		if (avail == 1) {
			unsigned int filterConfigFlags;
			herr_t retVal = H5Zget_filter_info(HDF5UTIL_FILTER, &filterConfigFlags);
			if (retVal >= 0 &&
					(filterConfigFlags & H5Z_FILTER_CONFIG_ENCODE_ENABLED) &&
					(filterConfigFlags & H5Z_FILTER_CONFIG_DECODE_ENABLED)) {
				dcpl.setShuffle();
#ifdef HDF5UTIL_USE_LZF
				dcpl.setFilter(HDF5UTIL_FILTER, H5Z_FLAG_OPTIONAL, 0, nullptr);
#else
				dcpl.setDeflate(deflateLevel); // compress with DEFLATE (zlib)
#endif
				std::vector<hsize_t> chunkDims;
				calcChunkDims(dataDims, chunkDims);

				dcpl.setChunk(rank, &chunkDims[0]); // needed when using compression
			}
		}

		DataSpace datasetSpace(rank, &dataDims[0]);
		DataSet dataset = file.createDataSet(datasetName, PredType::IEEE_F64LE /* dataset (in file) type */, datasetSpace, dcpl);
		dataset.write(getBeginPtr(container), hdf5MemoryType<T>());
	} catch (const H5::Exception& e) {
		THROW_EXCEPTION(IOException, "An error ocurred in HDF5 library with file " << filePath << ": " << e.getCDetailMsg());
	}
}

} // namespace HDF5Util
} // namespace Lab

#endif /* HDF5UTIL_H_ */
