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

#include "ShowImageMethod.h"

#include <string>
#include <vector>

#include "Log.h"
#include "Matrix.h"
#include "ParameterMap.h"
#include "Project.h"
#include "Util.h"



namespace Lab {

ShowImageMethod::ShowImageMethod(Project& project)
		: project_(project)
{
}

ShowImageMethod::~ShowImageMethod()
{
}

void
ShowImageMethod::execute()
{
	ParamMapPtr taskPM = project_.taskParameterMap();
	const auto dataDir      = taskPM->value<std::string>("data_dir");
	const auto xFile        = taskPM->value<std::string>("x_file");
	const auto xDataset     = taskPM->value<std::string>("x_dataset");
	const auto yFile        = taskPM->value<std::string>("y_file");
	const auto yDataset     = taskPM->value<std::string>("y_dataset");
	const auto zFile        = taskPM->value<std::string>("z_file");
	const auto zDataset     = taskPM->value<std::string>("z_dataset");
	const auto imageFile    = taskPM->value<std::string>("image_file");
	const auto imageDataset = taskPM->value<std::string>("image_dataset");

	Project::GridDataType projGridData;

	project_.loadImageFromHDF5(dataDir, imageFile, imageDataset, projGridData);
	project_.loadXYZFromHDF5(dataDir,
					xFile, xDataset,
					yFile, yDataset,
					zFile, zDataset,
					projGridData);

	if (taskPM->contains("points_x_file")) {
		std::vector<Project::PointType> projPointList;
		std::vector<float> pointsX, pointsY, pointsZ;
		const auto pointsXFile = taskPM->value<std::string>("points_x_file");
		project_.loadHDF5(dataDir + '/' + pointsXFile, "x", pointsX);
		const auto pointsYFile = taskPM->value<std::string>("points_y_file");
		project_.loadHDF5(dataDir + '/' + pointsYFile, "y", pointsY);
		const auto pointsZFile = taskPM->value<std::string>("points_z_file");
		project_.loadHDF5(dataDir + '/' + pointsZFile, "z", pointsZ);
		projPointList.resize(pointsX.size());
		Util::copyXYZFromSimpleVectors(pointsX, pointsY, pointsZ, projPointList);
		project_.showFigure3D(1, "Image", &projGridData, &projPointList,
					true, Figure::VISUALIZATION_RECTIFIED_LOG, Figure::COLORMAP_VIRIDIS);
	} else {
		project_.showFigure3D(1, "Image", &projGridData, Project::emptyPointList,
					true, Figure::VISUALIZATION_RECTIFIED_LOG, Figure::COLORMAP_VIRIDIS);
	}
}

} // namespace Lab
