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
#include "Visualization.h"
#include "XYZ.h"
#include "XYZValue.h"
#include "XZValue.h"



namespace Lab {

template<typename TFloat>
class DeviceSectorialScanMethod : public Method {
public:
	DeviceSectorialScanMethod(Project& project);
	virtual ~DeviceSectorialScanMethod() = default;

	virtual void execute();
private:
	DeviceSectorialScanMethod(const DeviceSectorialScanMethod&) = delete;
	DeviceSectorialScanMethod& operator=(const DeviceSectorialScanMethod&) = delete;
	DeviceSectorialScanMethod(DeviceSectorialScanMethod&&) = delete;
	DeviceSectorialScanMethod& operator=(DeviceSectorialScanMethod&&) = delete;

	void getSingleImageFromNetwork();
	void showSavedImage();
	void execContinuousNetworkImaging();
	void execTriggeredNetworkImaging();

	Project& project_;
	DeviceSectorialScanConfiguration<TFloat> config_;
};



template<typename TFloat>
DeviceSectorialScanMethod<TFloat>::DeviceSectorialScanMethod(Project& project)
		: project_(project)
{
	config_.load(*project_.getParamMap("config-device_sectorial_scan_method.txt"));
}

template<typename TFloat>
void
DeviceSectorialScanMethod<TFloat>::execute()
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

template<typename TFloat>
void
DeviceSectorialScanMethod<TFloat>::getSingleImageFromNetwork()
{
	const ParameterMap& pm = project_.taskParamMap();
	const auto outputDir = pm.value<std::string>("output_dir");

	project_.createDirectory(outputDir, false);

	auto acquisition = std::make_unique<NetworkDeviceSectorialScanAcquisition<TFloat>>(project_, config_);
	Matrix<XZValue<TFloat>> acqImageData;
	acquisition->execute(acqImageData);

	Matrix<XYZValue<TFloat>> imageData;
	Util::copy(acqImageData, imageData);

	std::vector<XYZ<float>> pointList = {{((config_.numElements - 1U) / 2.0f) * config_.pitch, 0.0, 0.0}};

	project_.showFigure3D(1, "Image", &imageData, &pointList,
				true, Visualization::Value::RECTIFIED_LOG, Colormap::Id::GRADIENT_VIRIDIS, config_.valueScale);

	project_.saveImageToHDF5(imageData, outputDir);
	project_.saveXYZToHDF5(imageData, outputDir);
}

template<typename TFloat>
void
DeviceSectorialScanMethod<TFloat>::showSavedImage()
{
	const ParameterMap& pm = project_.taskParamMap();
	const auto outputDir = pm.value<std::string>("output_dir");

	Matrix<XYZValue<TFloat>> imageData;

	project_.loadImageFromHDF5(outputDir, "image_value", "value", imageData);
	project_.loadXYZFromHDF5(outputDir,
					"image_x", "x",
					"image_y", "y",
					"image_z", "z",
					imageData);

	std::vector<XYZ<float>> pointList = {{((config_.numElements - 1U) / 2.0f) * config_.pitch, 0.0, 0.0}};

	project_.showFigure3D(1, "Image", &imageData, &pointList,
				true, Visualization::Value::RECTIFIED_LOG, Colormap::Id::GRADIENT_VIRIDIS, config_.valueScale);
}

template<typename TFloat>
void
DeviceSectorialScanMethod<TFloat>::execContinuousNetworkImaging()
{
	auto acquisition = std::make_unique<NetworkDeviceSectorialScanAcquisition<TFloat>>(project_, config_);
	Matrix<XZValue<TFloat>> acqImageData;
	Matrix<XYZValue<TFloat>> imageData;

	std::vector<XYZ<float>> pointList = {{((config_.numElements - 1U) / 2.0f) * config_.pitch, 0.0, 0.0}};

	int n = 0;
	Timer t;
	do {
		acquisition->execute(acqImageData);

		Util::copy(acqImageData, imageData);

		project_.showFigure3D(1, "Image", &imageData, &pointList,
					true, Visualization::Value::RECTIFIED_LOG, Colormap::Id::GRADIENT_VIRIDIS, config_.valueScale);

		if (++n == 10) {
			LOG_INFO << 10.0 / t.getTime() << " image/s";
			n = 0;
			t.reset();
		}
	} while (!project_.processingCancellationRequested());
}

template<typename TFloat>
void
DeviceSectorialScanMethod<TFloat>::execTriggeredNetworkImaging()
{
	project_.resetTrigger();

	const ParameterMap& pm = project_.taskParamMap();
	const auto outputDirPrefix = pm.value<std::string>("output_dir_prefix");

	auto acquisition = std::make_unique<NetworkDeviceSectorialScanAcquisition<TFloat>>(project_, config_);
	Matrix<XZValue<TFloat>> acqImageData;
	Matrix<XYZValue<TFloat>> imageData;

	std::vector<XYZ<float>> pointList = {{((config_.numElements - 1U) / 2.0f) * config_.pitch, 0.0, 0.0}};

	std::size_t triggerCount = 0;
	while (project_.waitForTrigger(&triggerCount)) {

		acquisition->execute(acqImageData);

		Util::copy(acqImageData, imageData);

		project_.showFigure3D(1, "Image", &imageData, &pointList,
					true, Visualization::Value::RECTIFIED_LOG, Colormap::Id::GRADIENT_VIRIDIS, config_.valueScale);

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
