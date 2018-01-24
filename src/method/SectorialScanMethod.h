#ifndef SECTORIALSCANMETHOD_H
#define SECTORIALSCANMETHOD_H

#include <cstddef> /* std::size_t */
#include <memory>
#include <sstream>
#include <string>

#include "Log.h"
#include "Matrix2.h"
#include "Method.h"
#include "NetworkSectorialScanAcquisition.h"
#include "ParameterMap.h"
#include "Project.h"
#include "SectorialScanConfiguration.h"
#include "Timer.h"
#include "Util.h"
#include "XZValue.h"



namespace Lab {

template<typename FloatType>
class SectorialScanMethod : public Method {
public:
	SectorialScanMethod(Project& project);
	virtual ~SectorialScanMethod();

	virtual void execute();

private:
	SectorialScanMethod(const SectorialScanMethod&);
	SectorialScanMethod& operator=(const SectorialScanMethod&);

	void getSingleImageFromNetwork();
	void showSavedImage();
	void execContinuousNetworkImaging();
	void execTriggeredNetworkImaging();

	Project& project_;
	SectorialScanConfiguration<FloatType> config_;
};



template<typename FloatType>
SectorialScanMethod<FloatType>::SectorialScanMethod(Project& project)
		: project_(project)
{
	config_.load(project_.loadParameterMap("config-sectorial_scan_method.txt"));
}

template<typename FloatType>
SectorialScanMethod<FloatType>::~SectorialScanMethod()
{
}

template<typename FloatType>
void
SectorialScanMethod<FloatType>::execute()
{
	switch (project_.method()) {
	case MethodType::sectorial_scan_sp_network:
		getSingleImageFromNetwork();
		break;
	case MethodType::sectorial_scan_sp_saved:
		showSavedImage();
		break;
	case MethodType::sectorial_scan_sp_network_continuous:
		execContinuousNetworkImaging();
		break;
	case MethodType::sectorial_scan_sp_network_trigger:
		execTriggeredNetworkImaging();
		break;
	default:
		THROW_EXCEPTION(InvalidParameterException, "Invalid method: " << static_cast<int>(project_.method()) << '.');
	}
}

template<typename FloatType>
void
SectorialScanMethod<FloatType>::getSingleImageFromNetwork()
{
	const std::string outputDir = project_.taskParameterMap()->value<std::string>("output_dir");

	auto acquisition = std::make_unique<NetworkSectorialScanAcquisition<FloatType>>(project_, config_);
	Matrix2<XZValue<FloatType>> imageData;
	acquisition->execute(imageData);

	std::vector<XZ<float>> pointList = {{((config_.numElements - 1U) / 2.0f) * config_.pitch, 0.0}};

	project_.showFigure3D(1, "Image", &imageData, &pointList,
				true, Figure::VISUALIZATION_RECTIFIED_LOG, Figure::COLORMAP_VIRIDIS, config_.valueScale);

	project_.saveImageToHDF5(imageData, outputDir);
}

template<typename FloatType>
void
SectorialScanMethod<FloatType>::showSavedImage()
{
	const std::string outputDir = project_.taskParameterMap()->value<std::string>("output_dir");

	Matrix2<XZValue<FloatType>> imageData;

	project_.loadImageFromHDF5(outputDir, imageData);

	std::vector<XZ<float>> pointList = {{((config_.numElements - 1U) / 2.0f) * config_.pitch, 0.0}};

	project_.showFigure3D(1, "Image", &imageData, &pointList,
				true, Figure::VISUALIZATION_RECTIFIED_LOG, Figure::COLORMAP_VIRIDIS, config_.valueScale);
}

template<typename FloatType>
void
SectorialScanMethod<FloatType>::execContinuousNetworkImaging()
{
	auto acquisition = std::make_unique<NetworkSectorialScanAcquisition<FloatType>>(project_, config_);
	Matrix2<XZValue<FloatType>> imageData;

	std::vector<XZ<float>> pointList = {{((config_.numElements - 1U) / 2.0f) * config_.pitch, 0.0}};

	int n = 0;
	Timer t;
	do {
		acquisition->execute(imageData);

		project_.showFigure3D(1, "Image", &imageData, &pointList,
					true, Figure::VISUALIZATION_RECTIFIED_LOG, Figure::COLORMAP_VIRIDIS, config_.valueScale);

		if (++n == 10) {
			LOG_INFO << 10.0 / t.getTime() << " image/s";
			n = 0;
			t.reset();
		}
	} while (!project_.processingCancellationRequested());
}

template<typename FloatType>
void
SectorialScanMethod<FloatType>::execTriggeredNetworkImaging()
{
	project_.resetTrigger();

	const std::string outputDirPrefix = project_.taskParameterMap()->value<std::string>("output_dir_prefix");

	auto acquisition = std::make_unique<NetworkSectorialScanAcquisition<FloatType>>(project_, config_);
	Matrix2<XZValue<FloatType>> imageData;

	std::vector<XZ<float>> pointList = {{((config_.numElements - 1U) / 2.0f) * config_.pitch, 0.0}};

	std::size_t triggerCount = 0;
	while (project_.waitForTrigger(&triggerCount)) {

		acquisition->execute(imageData);

		project_.showFigure3D(1, "Image", &imageData, &pointList,
					true, Figure::VISUALIZATION_RECTIFIED_LOG, Figure::COLORMAP_VIRIDIS, config_.valueScale);

		std::ostringstream out;
		out << outputDirPrefix << triggerCount;
		std::string outputDir = out.str();

		project_.createDirectory(outputDir, true);
		project_.saveImageToHDF5(imageData, outputDir);
	}
}

} // namespace Lab

#endif // SECTORIALSCANMETHOD_H
