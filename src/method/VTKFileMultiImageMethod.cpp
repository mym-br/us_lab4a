/***************************************************************************
 *  Copyright 2018 Marcelo Y. Matuda                                       *
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

#include "VTKFileMultiImageMethod.h"

#include <algorithm> /* swap */
#include <cmath> /* abs */
#include <cstddef> /* std::size_t */
#include <fstream>
#include <limits>
#include <string>
#include <vector>

#include "Exception.h"
#include "Log.h"
#include "Matrix2.h"
#include "ParameterMap.h"
#include "Project.h"
#include "Util.h"

#define VTK_FILE_NAME "multi_image.vtk"



namespace {

struct Cell {
	std::size_t i0;
	std::size_t i1;
	std::size_t i2;
	std::size_t i3;
	std::size_t i4;
	std::size_t i5;
	std::size_t i6;
	std::size_t i7;
	Cell(std::size_t i0, std::size_t i1, std::size_t i2, std::size_t i3,
			std::size_t i4, std::size_t i5, std::size_t i6, std::size_t i7)
		: i0(i0), i1(i1), i2(i2), i3(i3)
		, i4(i4), i5(i5), i6(i6), i7(i7) {}
};

} // namespace

namespace Lab {

VTKFileMultiImageMethod::VTKFileMultiImageMethod(Project& project)
		: project_{project}
{
}

VTKFileMultiImageMethod::~VTKFileMultiImageMethod()
{
}

void
VTKFileMultiImageMethod::execute()
{
	ConstParameterMapPtr taskPM = project_.taskParameterMap();
	const std::string imageBaseDir = taskPM->value<std::string>("image_dir");
	const std::string xFile        = taskPM->value<std::string>("x_file");
	const std::string xDataset     = taskPM->value<std::string>("x_dataset");
	const std::string zFile        = taskPM->value<std::string>("z_file");
	const std::string zDataset     = taskPM->value<std::string>("z_dataset");
	const std::string imageFile    = taskPM->value<std::string>("image_file");
	const std::string imageDataset = taskPM->value<std::string>("image_dataset");
	const float minY               = taskPM->value<float>(      "min_y"       , -10000.0, 10000.0);
	const float yStep              = taskPM->value<float>(      "y_step"      ,      0.0,   100.0);
	const float minDecibels        = taskPM->value<float>(      "min_decibels",   -100.0,    -1.0);
	const bool logScale            = taskPM->value<bool>(       "log_scale");
	const bool invertZ             = taskPM->value<bool>(       "invert_z");

	Project::GridDataType prevGridData;
	Matrix2<std::size_t> prevPointIndex;
	Project::GridDataType gridData;
	Matrix2<std::size_t> pointIndex;
	std::vector<XYZValue<float>> pointArray;
	std::vector<Cell> cellArray;
	const float minValue = Util::decibelsToLinear(minDecibels);
	const float valueCoeff = logScale ? (-1.0f / minDecibels) : (1.0f / (1.0f - minValue));
	float y = minY;

	auto calcValue = [&](float value) -> float {
		value = std::abs(value);
		Util::clip(value, minValue, 1.0f);
		if (logScale) {
			return (Util::linearToDecibels(value) - minDecibels) * valueCoeff;
		} else {
			return (value - minValue) * valueCoeff;
		}
	};
	auto storePoint = [&](const XZValue<float>& point, std::size_t& index) -> std::size_t {
		if (index != std::numeric_limits<std::size_t>::max()) {
			return index;
		}

		pointArray.emplace_back(point.x,
					y,
					invertZ ? -point.z : point.z,
					calcValue(point.value));
		index = pointArray.size() - 1U;
		return index;
	};
	auto addHexahedron = [&](std::size_t i, std::size_t j) {
		if (
				prevGridData(i    , j    ).value >= minValue ||
				prevGridData(i    , j + 1).value >= minValue ||
				prevGridData(i + 1, j    ).value >= minValue ||
				prevGridData(i + 1, j + 1).value >= minValue ||
				    gridData(i    , j    ).value >= minValue ||
				    gridData(i    , j + 1).value >= minValue ||
				    gridData(i + 1, j    ).value >= minValue ||
				    gridData(i + 1, j + 1).value >= minValue) {

			cellArray.emplace_back(
				storePoint(prevGridData(i    , j    ), prevPointIndex(i    , j    )),
				storePoint(prevGridData(i + 1, j    ), prevPointIndex(i + 1, j    )),
				storePoint(prevGridData(i + 1, j + 1), prevPointIndex(i + 1, j + 1)),
				storePoint(prevGridData(i    , j + 1), prevPointIndex(i    , j + 1)),
				storePoint(    gridData(i    , j    ),     pointIndex(i    , j    )),
				storePoint(    gridData(i + 1, j    ),     pointIndex(i + 1, j    )),
				storePoint(    gridData(i + 1, j + 1),     pointIndex(i + 1, j + 1)),
				storePoint(    gridData(i    , j + 1),     pointIndex(i    , j + 1)));
		}
	};

	for (unsigned int acqNumber = 0; ; ++acqNumber, y += yStep) {
		std::string imageDir = FileUtil::path(imageBaseDir, "/", acqNumber);
		if (!project_.directoryExists(imageDir)) {
			break;
		}
		LOG_DEBUG << "##### acqNumber: " << acqNumber;

		LOG_DEBUG << "Loading the image...";
		project_.loadHDF5(imageDir + '/' + imageFile, imageDataset, gridData, Util::CopyToValueOp());
		LOG_DEBUG << "Loading the X coordinates...";
		project_.loadHDF5(imageDir + '/' + xFile, xDataset, gridData, Util::CopyToXOp());
		LOG_DEBUG << "Loading the Z coordinates...";
		project_.loadHDF5(imageDir + '/' + zFile, zDataset, gridData, Util::CopyToZOp());
		if (gridData.n1() < 2 || gridData.n2() < 2) {
			THROW_EXCEPTION(InvalidValueException, "The image is too small.");
		}
		pointIndex.resize(gridData.n1(), gridData.n2());
		pointIndex = std::numeric_limits<std::size_t>::max();
		if (acqNumber == 0) {
			std::swap(prevGridData, gridData);
			std::swap(prevPointIndex, pointIndex);
			continue;
		}

		for (std::size_t i = 0; i < gridData.n1() - 1; ++i) {
			for (std::size_t j = 0; j < gridData.n2() - 1; ++j) {
				addHexahedron(i, j);
			}
		}

		std::swap(prevGridData, gridData);
		std::swap(prevPointIndex, pointIndex);
	}

	std::ofstream out{project_.directory() + '/' + VTK_FILE_NAME};
	out <<
		"# vtk DataFile Version 3.0\n"
		"Multi image\n"
		"ASCII\n"
		"DATASET UNSTRUCTURED_GRID\n"
		"POINTS " << pointArray.size() << " float\n";
	for (XYZValue<float>& p : pointArray) {
		out << p.x << ' ' << p.y << ' ' << p.z << '\n';
	}

	out << "CELLS " << cellArray.size() << ' '
		<< cellArray.size() * (1U /* the first column indicates the number of points per cell */ *
					8U /* indexes of the cell points */);
	for (Cell& c : cellArray) {
		out << "8 " << c.i0;
		out << ' '  << c.i1;
		out << ' '  << c.i2;
		out << ' '  << c.i3;
		out << ' '  << c.i4;
		out << ' '  << c.i5;
		out << ' '  << c.i6;
		out << ' '  << c.i7 << '\n';
	}

	out << "CELL_TYPES " << cellArray.size() << '\n';
	for (unsigned int i = 0; i < cellArray.size(); ++i) {
		out << "12\n"; // VTK_HEXAHEDRON
	}

	out << "POINT_DATA " << pointArray.size() << '\n';
	out << "SCALARS pointValues float\n";
	out << "LOOKUP_TABLE pointValuesTable\n";
	for (XYZValue<float>& p : pointArray) {
		out << p.value << '\n';
	}
}

} // namespace Lab
