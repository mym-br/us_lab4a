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

#include <cmath> /* abs */
#include <cstddef> /* std::size_t */
#include <fstream>
#include <limits>
#include <string>
#include <vector>

#include "Exception.h"
#include "Log.h"
#include "Matrix.h"
#include "ParameterMap.h"
#include "Project.h"
#include "Util.h"
#include "XYZValue.h"



namespace {

constexpr const char* VTK_FILE_NAME = "multi_image.vtk";
constexpr int NUM_CELL_POINTS = 8;

struct Cell {
	std::size_t i[NUM_CELL_POINTS];
};

} // namespace

namespace Lab {

VTKFileMultiImageMethod::VTKFileMultiImageMethod(Project& project)
		: project_(project)
{
}

void
VTKFileMultiImageMethod::execute()
{
	const auto& taskPM = project_.taskParameterMap();
	const auto imageBaseDir = taskPM.value<std::string>("image_dir");
	const auto imagPM = project_.loadChildParameterMap(taskPM, "imag_config_file");
	const auto xFile        = imagPM->value<std::string>("x_file");
	const auto xDataset     = imagPM->value<std::string>("x_dataset");
	const auto yFile        = imagPM->value<std::string>("y_file");
	const auto yDataset     = imagPM->value<std::string>("y_dataset");
	const auto zFile        = imagPM->value<std::string>("z_file");
	const auto zDataset     = imagPM->value<std::string>("z_dataset");
	const auto imageFile    = imagPM->value<std::string>("image_file");
	const auto imageDataset = imagPM->value<std::string>("image_dataset");
	const auto minDecibels  = imagPM->value<float>(      "min_decibels", -100.0, -1.0);
	const auto logScale     = imagPM->value<bool>(       "log_scale");
	const auto invertZ      = imagPM->value<bool>(       "invert_z");

	Project::GridDataType prevGridData;
	Matrix<std::size_t> prevPointIndex;
	Project::GridDataType gridData;
	Matrix<std::size_t> pointIndex;
	std::vector<XYZValue<float>> pointArray;
	std::vector<Cell> cellArray;
	const float minValue = Util::decibelsToLinear(minDecibels);
	const float valueCoeff = logScale ? (-1.0f / minDecibels) : (1.0f / (1.0f - minValue));

	auto calcValue = [&](float value) -> float {
		value = std::abs(value);
		Util::clip(value, minValue, 1.0f);
		if (logScale) {
			return (Util::linearToDecibels(value) - minDecibels) * valueCoeff;
		} else {
			return (value - minValue) * valueCoeff;
		}
	};
	auto storePoint = [&](const XYZValue<float>& point, std::size_t& index) -> std::size_t {
		if (index != std::numeric_limits<std::size_t>::max()) {
			return index;
		}

		pointArray.emplace_back(point.x,
					invertZ ? -point.z : point.z, // z -> y
					point.y,                      // y -> z
					calcValue(point.value));
		index = pointArray.size() - 1U;
		return index;
	};
	auto addHexahedron = [&](std::size_t i, std::size_t j) {
		if (
				prevGridData(i    , j    ).value >= minValue ||
				prevGridData(i + 1, j    ).value >= minValue ||
				prevGridData(i + 1, j + 1).value >= minValue ||
				prevGridData(i    , j + 1).value >= minValue ||
				    gridData(i    , j    ).value >= minValue ||
				    gridData(i + 1, j    ).value >= minValue ||
				    gridData(i + 1, j + 1).value >= minValue ||
				    gridData(i    , j + 1).value >= minValue) {

			Cell cell;
			cell.i[0] = storePoint(prevGridData(i    , j    ), prevPointIndex(i    , j    ));
			cell.i[1] = storePoint(prevGridData(i + 1, j    ), prevPointIndex(i + 1, j    ));
			cell.i[2] = storePoint(prevGridData(i + 1, j + 1), prevPointIndex(i + 1, j + 1));
			cell.i[3] = storePoint(prevGridData(i    , j + 1), prevPointIndex(i    , j + 1));
			cell.i[4] = storePoint(    gridData(i    , j    ),     pointIndex(i    , j    ));
			cell.i[5] = storePoint(    gridData(i + 1, j    ),     pointIndex(i + 1, j    ));
			cell.i[6] = storePoint(    gridData(i + 1, j + 1),     pointIndex(i + 1, j + 1));
			cell.i[7] = storePoint(    gridData(i    , j + 1),     pointIndex(i    , j + 1));
			cellArray.push_back(cell);
		}
	};

	for (unsigned int acqNumber = 0; ; ++acqNumber) {
		std::string imageDir = FileUtil::path(imageBaseDir, "/", acqNumber);
		if (!project_.directoryExists(imageDir)) {
			if (acqNumber == 0) {
				THROW_EXCEPTION(UnavailableResourceException, "No file found at " << imageDir << '.');
			}
			break;
		}
		LOG_INFO << "ACQ " << acqNumber;

		LOG_DEBUG << "Loading the image...";
		project_.loadHDF5(imageDir + '/' + imageFile, imageDataset, gridData, Util::CopyToValueOp());
		LOG_DEBUG << "Loading the X coordinates...";
		project_.loadHDF5(imageDir + '/' + xFile, xDataset, gridData, Util::CopyToXOp());
		LOG_DEBUG << "Loading the Y coordinates...";
		project_.loadHDF5(imageDir + '/' + yFile, yDataset, gridData, Util::CopyToYOp());
		LOG_DEBUG << "Loading the Z coordinates...";
		project_.loadHDF5(imageDir + '/' + zFile, zDataset, gridData, Util::CopyToZOp());
		if (gridData.n1() < 2 || gridData.n2() < 2) {
			THROW_EXCEPTION(InvalidValueException, "The image is too small.");
		}
		pointIndex.resize(gridData.n1(), gridData.n2());
		pointIndex = std::numeric_limits<std::size_t>::max();
		if (acqNumber == 0) {
			prevGridData.swap(gridData);
			prevPointIndex.swap(pointIndex);
			continue;
		}

		for (std::size_t i = 0; i < gridData.n1() - 1; ++i) {
			for (std::size_t j = 0; j < gridData.n2() - 1; ++j) {
				addHexahedron(i, j);
			}
		}

		prevGridData.swap(gridData);
		prevPointIndex.swap(pointIndex);
	}

	std::ofstream out(project_.expDirectory() + '/' + VTK_FILE_NAME);
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
		<< cellArray.size() * (1U /* the first column indicates the number of points per cell */ +
					8U /* indexes of the cell points */) << '\n';
	for (Cell& c : cellArray) {
		out << "8"; // number of points per cell
		for (int i = 0; i < NUM_CELL_POINTS; ++i) {
			out << ' '  << c.i[i];
		}
		out << '\n';
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
