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

#ifndef DEVICESECTORIALSCANMETHOD_H
#define DEVICESECTORIALSCANMETHOD_H

#include <cstddef> /* std::size_t */
#include <memory>
#include <sstream>
#include <string>

#include "Colormap.h"
#include "DeviceSectorialScanConfiguration.h"
#include "Log.h"
#include "Matrix.h"
#include "Method.h"
#include "NetworkDeviceSectorialScanAcquisition.h"
#include "ParameterMap.h"
#include "Project.h"
#include "Timer.h"
#include "Util.h"
#include "XYZ.h"
#include "XYZValue.h"
#include "XZValue.h"



namespace Lab {

template<typename FloatType>
class DeviceSectorialScanMethod : public Method {
public:
	DeviceSectorialScanMethod(Project& project);
	virtual ~DeviceSectorialScanMethod();

	virtual void execute();

private:
	DeviceSectorialScanMethod(const DeviceSectorialScanMethod&);
	DeviceSectorialScanMethod& operator=(const DeviceSectorialScanMethod&);

	void getSingleImageFromNetwork();
	void showSavedImage();
	void execContinuousNetworkImaging();
	void execTriggeredNetworkImaging();

	Project& project_;
	DeviceSectorialScanConfiguration<FloatType> config_;
};



template<typename FloatType>
DeviceSectorialScanMethod<FloatType>::DeviceSectorialScanMethod(Project& project)
		: project_(project)
{
	config_.load(project_.loadParameterMap("config-device_sectorial_scan_method.txt"));
}

template<typename FloatType>
DeviceSectorialScanMethod<FloatType>::~DeviceSectorialScanMethod()
{
}

template<typename FloatType>
void
DeviceSectorialScanMethod<FloatType>::execute()
{
	switch (project_.method()) {
	case MethodEnum::device_sectorial_scan_sp_network:
		getSingleImageFromNetwork();
		break;
	case MethodEnum::device_sectorial_scan_sp_saved:
		showSavedImage();
		break;
	case MethodEnum::device_sectorial_scan_sp_network_continuous:
		execContinuousNetworkImaging();
		break;
	case MethodEnum::device_sectorial_scan_sp_network_trigger:
		execTriggeredNetworkImaging();
		break;
	default:
		THROW_EXCEPTION(InvalidParameterException, "Invalid method: " << static_cast<int>(project_.method()) << '.');
	}
}

template<typename FloatType>
void
DeviceSectorialScanMethod<FloatType>::getSingleImageFromNetwork()
{
	ParamMapPtr pm = project_.taskParameterMap();
	const auto outputDir = pm->value<std::string>("output_dir");

	project_.createDirectory(outputDir, false);

	auto acquisition = std::make_unique<NetworkDeviceSectorialScanAcquisition<FloatType>>(project_, config_);
	Matrix<XZValue<FloatType>> acqImageData;
	acquisition->execute(acqImageData);

	Matrix<XYZValue<FloatType>> imageData;
	Util::copy(acqImageData, imageData);

	std::vector<XYZ<float>> pointList = {{((config_.numElements - 1U) / 2.0f) * config_.pitch, 0.0, 0.0}};

	project_.showFigure3D(1, "Image", &imageData, &pointList,
				true, Figure::VISUALIZATION_RECTIFIED_LOG, Colormap::GRADIENT_VIRIDIS, config_.valueScale);

	project_.saveImageToHDF5(imageData, outputDir);
	project_.saveXYZToHDF5(imageData, outputDir);
}

template<typename FloatType>
void
DeviceSectorialScanMethod<FloatType>::showSavedImage()
{
	ParamMapPtr pm = project_.taskParameterMap();
	const auto outputDir = pm->value<std::string>("output_dir");

	Matrix<XYZValue<FloatType>> imageData;

	project_.loadImageFromHDF5(outputDir, "image_value", "value", imageData);
	project_.loadXYZFromHDF5(outputDir,
					"image_x", "x",
					"image_y", "y",
					"image_z", "z",
					imageData);

	std::vector<XYZ<float>> pointList = {{((config_.numElements - 1U) / 2.0f) * config_.pitch, 0.0, 0.0}};

	project_.showFigure3D(1, "Image", &imageData, &pointList,
				true, Figure::VISUALIZATION_RECTIFIED_LOG, Colormap::GRADIENT_VIRIDIS, config_.valueScale);
}

template<typename FloatType>
void
DeviceSectorialScanMethod<FloatType>::execContinuousNetworkImaging()
{
	auto acquisition = std::make_unique<NetworkDeviceSectorialScanAcquisition<FloatType>>(project_, config_);
	Matrix<XZValue<FloatType>> acqImageData;
	Matrix<XYZValue<FloatType>> imageData;

	std::vector<XYZ<float>> pointList = {{((config_.numElements - 1U) / 2.0f) * config_.pitch, 0.0, 0.0}};

	int n = 0;
	Timer t;
	do {
		acquisition->execute(acqImageData);

		Util::copy(acqImageData, imageData);

		project_.showFigure3D(1, "Image", &imageData, &pointList,
					true, Figure::VISUALIZATION_RECTIFIED_LOG, Colormap::GRADIENT_VIRIDIS, config_.valueScale);

		if (++n == 10) {
			LOG_INFO << 10.0 / t.getTime() << " image/s";
			n = 0;
			t.reset();
		}
	} while (!project_.processingCancellationRequested());
}

template<typename FloatType>
void
DeviceSectorialScanMethod<FloatType>::execTriggeredNetworkImaging()
{
	project_.resetTrigger();

	ParamMapPtr pm = project_.taskParameterMap();
	const auto outputDirPrefix = pm->value<std::string>("output_dir_prefix");

	auto acquisition = std::make_unique<NetworkDeviceSectorialScanAcquisition<FloatType>>(project_, config_);
	Matrix<XZValue<FloatType>> acqImageData;
	Matrix<XYZValue<FloatType>> imageData;

	std::vector<XYZ<float>> pointList = {{((config_.numElements - 1U) / 2.0f) * config_.pitch, 0.0, 0.0}};

	std::size_t triggerCount = 0;
	while (project_.waitForTrigger(&triggerCount)) {

		acquisition->execute(acqImageData);

		Util::copy(acqImageData, imageData);

		project_.showFigure3D(1, "Image", &imageData, &pointList,
					true, Figure::VISUALIZATION_RECTIFIED_LOG, Colormap::GRADIENT_VIRIDIS, config_.valueScale);

		std::ostringstream out;
		out << outputDirPrefix << triggerCount;
		std::string outputDir = out.str();

		project_.createDirectory(outputDir, true);
		project_.saveImageToHDF5(imageData, outputDir);
		project_.saveXYZToHDF5(imageData, outputDir);
	}
}

} // namespace Lab

#endif // DEVICESECTORIALSCANMETHOD_H
