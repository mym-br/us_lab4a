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

#include <string>
#include <vector>

#include "Log.h"
#include "ParameterMap.h"
#include "Project.h"
#include "Util.h"



namespace Lab {

MultiLayerImageMethod::MultiLayerImageMethod(Project& project)
		: project_{project}
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
	const float minY               = taskPM->value<float>(      "min_y"       , -10000.0, 10000.0);
	const float yStep              = taskPM->value<float>(      "y_step"      ,      0.0,   100.0);
	const float minDecibels        = taskPM->value<float>(      "min_decibels",   -100.0,    -1.0);
	const bool logScale            = taskPM->value<bool>(       "log_scale");

	Project::GridDataType projGridData;
	std::vector<XYZValue<float>> pointArray;
	std::vector<unsigned int> indexArray;
	const float minValue = Util::decibelsToLinear(minDecibels);
	unsigned int j1, j2, j3;

	float y = minY;
	for (unsigned int acqNumber = 0; ; ++acqNumber, y += yStep) {
		std::string imageDir = FileUtil::path(imageBaseDir, "/", acqNumber);
		if (!project_.directoryExists(imageDir)) {
			break;
		}

		LOG_DEBUG << "Loading the image...";
		project_.loadHDF5(imageDir + '/' + imageFile, imageDataset, projGridData, Util::CopyToValueOp());
		LOG_DEBUG << "Loading the X coordinates...";
		project_.loadHDF5(imageDir + '/' + xFile, xDataset, projGridData, Util::CopyToXOp());
		LOG_DEBUG << "Loading the Z coordinates...";
		project_.loadHDF5(imageDir + '/' + zFile, zDataset, projGridData, Util::CopyToZOp());

		for (unsigned int i = 0; i < projGridData.n1() - 1; ++i) {
			unsigned int jA, jB;
			if (i & 1U) {
				jA = projGridData.n2();
				jB = 0;
			} else {
				jA = 0;
				jB = projGridData.n2();
			}
//TODO: Don't repeat vertices.
			for (unsigned int j = 0; j < projGridData.n2() - 1; ++j, ++jA, ++jB) {
				if (j & 1U) {
					j1 = jA;
					j2 = jB;
					j3 = jA + 1;
					if (projGridData(i, j1).value >= minValue &&
							projGridData(i, j2).value >= minValue &&
							projGridData(i, j3).value >= minValue) {
						const float coeff = 1.0f / (1.0f - minValue);
						// The vertices of a triangle.
						pointArray.emplace_back(
								projGridData(i, j1).x,
								y,
								projGridData(i, j1).z,
								(projGridData(i, j1).value - minValue) * coeff); // linear
						pointArray.emplace_back(
								projGridData(i, j2).x,
								y,
								projGridData(i, j2).z,
								(projGridData(i, j2).value - minValue) * coeff); // linear
						pointArray.emplace_back(
								projGridData(i, j3).x,
								y,
								projGridData(i, j3).z,
								(projGridData(i, j3).value - minValue) * coeff); // linear
						// The indexes of the triangle.
						indexArray.push_back(pointArray.size() - 3);
						indexArray.push_back(pointArray.size() - 2);
						indexArray.push_back(pointArray.size() - 1);
					}

					j1 = jA + 1;
					j2 = jB;
					j3 = jB + 1;
					if (projGridData(i, j1).value >= minValue &&
							projGridData(i, j2).value >= minValue &&
							projGridData(i, j3).value >= minValue) {
						const float coeff = 1.0f / (1.0f - minValue);
						// The vertices of a triangle.
						pointArray.emplace_back(
								projGridData(i, j1).x,
								y,
								projGridData(i, j1).z,
								(projGridData(i, j1).value - minValue) * coeff); // linear
						pointArray.emplace_back(
								projGridData(i, j2).x,
								y,
								projGridData(i, j2).z,
								(projGridData(i, j2).value - minValue) * coeff); // linear
						pointArray.emplace_back(
								projGridData(i, j3).x,
								y,
								projGridData(i, j3).z,
								(projGridData(i, j3).value - minValue) * coeff); // linear
						// The indexes of the triangle.
						indexArray.push_back(pointArray.size() - 3);
						indexArray.push_back(pointArray.size() - 2);
						indexArray.push_back(pointArray.size() - 1);
					}
				} else {
					j1 = jB;
					j2 = jB + 1;
					j3 = jA;
					if (projGridData(i, j1).value >= minValue &&
							projGridData(i, j2).value >= minValue &&
							projGridData(i, j3).value >= minValue) {
						const float coeff = 1.0f / (1.0f - minValue);
						// The vertices of a triangle.
						pointArray.emplace_back(
								projGridData(i, j1).x,
								y,
								projGridData(i, j1).z,
								(projGridData(i, j1).value - minValue) * coeff); // linear
						pointArray.emplace_back(
								projGridData(i, j2).x,
								y,
								projGridData(i, j2).z,
								(projGridData(i, j2).value - minValue) * coeff); // linear
						pointArray.emplace_back(
								projGridData(i, j3).x,
								y,
								projGridData(i, j3).z,
								(projGridData(i, j3).value - minValue) * coeff); // linear
						// The indexes of the triangle.
						indexArray.push_back(pointArray.size() - 3);
						indexArray.push_back(pointArray.size() - 2);
						indexArray.push_back(pointArray.size() - 1);
					}

					j1 = jA;
					j2 = jB + 1;
					j3 = jA + 1;
					if (projGridData(i, j1).value >= minValue &&
							projGridData(i, j2).value >= minValue &&
							projGridData(i, j3).value >= minValue) {
						const float coeff = 1.0f / (1.0f - minValue);
						// The vertices of a triangle.
						pointArray.emplace_back(
								projGridData(i, j1).x,
								y,
								projGridData(i, j1).z,
								(projGridData(i, j1).value - minValue) * coeff); // linear
						pointArray.emplace_back(
								projGridData(i, j2).x,
								y,
								projGridData(i, j2).z,
								(projGridData(i, j2).value - minValue) * coeff); // linear
						pointArray.emplace_back(
								projGridData(i, j3).x,
								y,
								projGridData(i, j3).z,
								(projGridData(i, j3).value - minValue) * coeff); // linear
						// The indexes of the triangle.
						indexArray.push_back(pointArray.size() - 3);
						indexArray.push_back(pointArray.size() - 2);
						indexArray.push_back(pointArray.size() - 1);
					}
				}
			}
		}
	}

	project_.showMultiLayer3D(1, "Figure", pointArray, indexArray);
}

} // namespace Lab
