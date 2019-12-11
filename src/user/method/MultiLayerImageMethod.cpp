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
#include <limits>
#include <string>
#include <vector>

#include "Exception.h"
#include "Log.h"
#include "ParameterMap.h"
#include "Project.h"
#include "Util.h"



namespace Lab {

MultiLayerImageMethod::MultiLayerImageMethod(Project& project)
		: project_(project)
{
}

void
MultiLayerImageMethod::execute()
{
	const ParameterMap& taskPM = project_.taskParamMap();
	const auto imageBaseDir = taskPM.value<std::string>("image_dir");
	const ParamMapPtr imagPM = project_.getSubParamMap("imag_config_file");
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

	Project::GridDataType projGridData;
	std::vector<XYZValue<float>> pointArray;
	std::vector<unsigned int> indexArray;
	const float minValue = Util::decibelsToLinear(minDecibels);
	const float valueCoeff = logScale ? (-1.0f / minDecibels) : (1.0f / (1.0f - minValue));
	unsigned int j1, j2, j3;

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
						projGridData(i, j1).y,
						invertZ ? -projGridData(i, j1).z : projGridData(i, j1).z,
						calcValue(projGridData(i, j1).value)});
			storePoint(endIndex, XYZValue<float>{
						projGridData(i, j2).x,
						projGridData(i, j2).y,
						invertZ ? -projGridData(i, j2).z : projGridData(i, j2).z,
						calcValue(projGridData(i, j2).value)});
			storePoint(endIndex, XYZValue<float>{
						projGridData(i, j3).x,
						projGridData(i, j3).y,
						invertZ ? -projGridData(i, j3).z : projGridData(i, j3).z,
						calcValue(projGridData(i, j3).value)});
		}
	};

	float minY = std::numeric_limits<float>::max();
	float maxY = std::numeric_limits<float>::lowest();

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
		project_.loadHDF5(imageDir + '/' + imageFile, imageDataset, projGridData, Util::CopyToValueOp());
		LOG_DEBUG << "Loading the X coordinates...";
		project_.loadHDF5(imageDir + '/' + xFile, xDataset, projGridData, Util::CopyToXOp());
		LOG_DEBUG << "Loading the Y coordinates...";
		project_.loadHDF5(imageDir + '/' + yFile, yDataset, projGridData, Util::CopyToYOp());
		LOG_DEBUG << "Loading the Z coordinates...";
		project_.loadHDF5(imageDir + '/' + zFile, zDataset, projGridData, Util::CopyToZOp());

		if (projGridData.n1() < 2 || projGridData.n2() < 2) continue;

		if (projGridData(0, 0).y < minY) minY = projGridData(0, 0).y;
		if (projGridData(0, 0).y > maxY) maxY = projGridData(0, 0).y;

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
