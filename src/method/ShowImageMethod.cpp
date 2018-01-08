#include "ShowImageMethod.h"

#include <string>
#include <vector>

#include "Log.h"
#include "Matrix2.h"
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
	ConstParameterMapPtr taskPM = project_.taskParameterMap();
	std::string dataDir   = taskPM->value<std::string>("data_dir");
	std::string xFile     = taskPM->value<std::string>("x_file");
	std::string zFile     = taskPM->value<std::string>("z_file");
	std::string imageFile = taskPM->value<std::string>("image_file");

	Project::GridDataType projGridData;
	{
		Matrix2<float> values;
		LOG_DEBUG << "Loading the image...";
		project_.loadHDF5(dataDir + '/' + imageFile, "image", values);
		projGridData.resize(values.n1(), values.n2());
		Util::copyValueFromSimpleMatrix(values, projGridData);
	}
	{
		Matrix2<float> x, z;
		LOG_DEBUG << "Loading the X coordinates...";
		project_.loadHDF5(dataDir + '/' + xFile, "x", x);
		LOG_DEBUG << "Loading the Z coordinates...";
		project_.loadHDF5(dataDir + '/' + zFile, "z", z);
		Util::copyXZFromSimpleMatrices(x, z, projGridData);
	}

	if (taskPM->contains("points_x_file")) {
		std::vector<Project::PointType> projPointList;
		std::vector<float> pointsX, pointsZ;
		std::string pointsXFile = taskPM->value<std::string>("points_x_file");
		project_.loadHDF5(dataDir + '/' + pointsXFile, "x", pointsX);
		std::string pointsZFile = taskPM->value<std::string>("points_z_file");
		project_.loadHDF5(dataDir + '/' + pointsZFile, "z", pointsZ);
		projPointList.resize(pointsX.size());
		Util::copyXZFromSimpleVectors(pointsX, pointsZ, projPointList);
		project_.showFigure3D(1, "Image", &projGridData, &projPointList,
					true, Figure::VISUALIZATION_ENVELOPE_LOG, Figure::COLORMAP_VIRIDIS);
	} else {
		project_.showFigure3D(1, "Image", &projGridData, Project::emptyPointList,
					true, Figure::VISUALIZATION_ENVELOPE_LOG, Figure::COLORMAP_VIRIDIS);
	}
}

} // namespace Lab
