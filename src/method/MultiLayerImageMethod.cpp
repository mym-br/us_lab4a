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

#include "MultiLayerImageMethod.h"

#include <cmath>
#include <string>
#include <vector>

#include "Log.h"
#include "ParameterMap.h"
#include "Project.h"
#include "Util.h"



namespace Lab {

MultiLayerImageMethod::MultiLayerImageMethod(Project& project)
		: project_(project)
{
}

MultiLayerImageMethod::~MultiLayerImageMethod()
{
}

void
MultiLayerImageMethod::execute()
{
	ConstParameterMapPtr taskPM = project_.taskParameterMap();
	const std::string imageBaseDir = taskPM->value<std::string>("image_dir");
	const std::string xFile        = taskPM->value<std::string>("x_file");
	const std::string xDataset     = taskPM->value<std::string>("x_dataset");
	const std::string zFile        = taskPM->value<std::string>("z_file");
	const std::string zDataset     = taskPM->value<std::string>("z_dataset");
	const std::string imageFile    = taskPM->value<std::string>("image_file");
	const std::string imageDataset = taskPM->value<std::string>("image_dataset");
	const float minDecibels        = taskPM->value<float>(      "min_decibels",   -100.0,    -1.0);
	const bool logScale            = taskPM->value<bool>(       "log_scale");
	const bool invertZ             = taskPM->value<bool>(       "invert_z");

	ConstParameterMapPtr scanPM = project_.loadChildParameterMap(taskPM, "scan_config_file");
	const float minY  = scanPM->value<float>("min_y" , -10000.0, 10000.0);
	const float yStep = scanPM->value<float>("y_step",      0.0,   100.0);

	Project::GridDataType projGridData;
	std::vector<XYZValue<float>> pointArray;
	std::vector<unsigned int> indexArray;
	const float minValue = Util::decibelsToLinear(minDecibels);
	const float valueCoeff = logScale ? (-1.0f / minDecibels) : (1.0f / (1.0f - minValue));
	unsigned int j1, j2, j3;
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
	auto storePoint = [&](unsigned int endIndex, const XYZValue<float>& point) {
		for (unsigned int i = 1; i <= 3; ++i) {
			if (endIndex < i) break;
			if (pointArray[endIndex - i] == point) {
				indexArray.push_back(endIndex - i);
				return;
			}
		}
		pointArray.push_back(point);
		indexArray.push_back(pointArray.size() - 1U);
	};
	auto addTriangle = [&](unsigned int i) {
		if (projGridData(i, j1).value >= minValue ||
				projGridData(i, j2).value >= minValue ||
				projGridData(i, j3).value >= minValue) {

			const unsigned int endIndex = pointArray.size();
			storePoint(endIndex, XYZValue<float>{
						projGridData(i, j1).x,
						y,
						invertZ ? -projGridData(i, j1).z : projGridData(i, j1).z,
						calcValue(projGridData(i, j1).value)});
			storePoint(endIndex, XYZValue<float>{
						projGridData(i, j2).x,
						y,
						invertZ ? -projGridData(i, j2).z : projGridData(i, j2).z,
						calcValue(projGridData(i, j2).value)});
			storePoint(endIndex, XYZValue<float>{
						projGridData(i, j3).x,
						y,
						invertZ ? -projGridData(i, j3).z : projGridData(i, j3).z,
						calcValue(projGridData(i, j3).value)});
		}
	};

	for (unsigned int acqNumber = 0; ; ++acqNumber, y += yStep) {
		std::string imageDir = FileUtil::path(imageBaseDir, "/", acqNumber);
		if (!project_.directoryExists(imageDir)) {
			break;
		}
		LOG_DEBUG << "##### acqNumber: " << acqNumber;

		LOG_DEBUG << "Loading the image...";
		project_.loadHDF5(imageDir + '/' + imageFile, imageDataset, projGridData, Util::CopyToValueOp());
		LOG_DEBUG << "Loading the X coordinates...";
		project_.loadHDF5(imageDir + '/' + xFile, xDataset, projGridData, Util::CopyToXOp());
		LOG_DEBUG << "Loading the Z coordinates...";
		project_.loadHDF5(imageDir + '/' + zFile, zDataset, projGridData, Util::CopyToZOp());

		if (projGridData.n1() < 2 || projGridData.n2() < 2) continue;

		for (unsigned int i = 0; i < projGridData.n1() - 1; ++i) {
			unsigned int jA, jB;
			if (i & 1U) {
				jA = projGridData.n2(); // note: exceeds the maximum column
				jB = 0;
			} else {
				jA = 0;
				jB = projGridData.n2(); // note: exceeds the maximum column
			}
			for (unsigned int j = 0; j < projGridData.n2() - 1; ++j, ++jA, ++jB) {
				if (j & 1U) {
					j1 = jA;
					j2 = jB;
					j3 = jA + 1;
					addTriangle(i);

					j1 = jA + 1;
					j2 = jB;
					j3 = jB + 1;
					addTriangle(i);
				} else {
					j1 = jB;
					j2 = jB + 1;
					j3 = jA;
					addTriangle(i);

					j1 = jA;
					j2 = jB + 1;
					j3 = jA + 1;
					addTriangle(i);
				}
			}
		}
	}

	// Add points to indicate the original limits.
	if (!pointArray.empty() && projGridData.n1() >= 2 && projGridData.n2() >= 2) {
		const float maxY = y - yStep;
		const auto& firstPoint = projGridData(0, 0);
		float minX = firstPoint.x;
		float maxX = firstPoint.x;
		float z = invertZ ? -firstPoint.z : firstPoint.z;
		float minZ = z;
		float maxZ = z;
		for (auto& elem : projGridData) {
			if (elem.x < minX) minX = elem.x;
			else if (elem.x > maxX) maxX = elem.x;

			z = invertZ ? -elem.z : elem.z;
			if (elem.z < minZ) minZ = z;
			else if (elem.z > maxZ) maxZ = z;
		}

		pointArray.push_back(XYZValue<float>{minX, minY, minZ, 0.0});
		pointArray.push_back(XYZValue<float>{maxX, maxY, maxZ, 0.0});
	}

	project_.showMultiLayer3D(1, "Figure", pointArray, indexArray);
}

} // namespace Lab
