#ifndef SECTORIALSCANMETHOD_H
#define SECTORIALSCANMETHOD_H

#include <memory>
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

	void fillConfiguration();
	void getSingleImageFromNetwork();
	void showSavedImage();
	void execContinuousNetworkImaging();

	Project& project_;
	SectorialScanConfiguration<FloatType> config_;
};



template<typename FloatType>
void
SectorialScanMethod<FloatType>::fillConfiguration()
{
	ConstParameterMapPtr pm = project_.loadParameterMap("config-sectorial_scan_method.txt");
	config_.numElementsMux         = pm->value<unsigned int>("num_elements_mux"        ,       8,      1024);
	config_.numElements            = pm->value<unsigned int>("num_elements"            ,       8, config_.numElementsMux);
	config_.baseElement            = pm->value<unsigned int>("base_element"            ,       0, config_.numElementsMux - config_.numElements);
	config_.signalMode             = pm->value<unsigned int>("signal_mode"             ,       0,         1);
	config_.pitch                  = pm->value<FloatType>(   "pitch"                   , 0.01e-3,   10.0e-3);
	config_.centerFrequency        = pm->value<FloatType>(   "center_frequency"        ,   100.0,   100.0e6);
	config_.gain                   = pm->value<FloatType>(   "gain"                    ,     0.0,     100.0);
	config_.propagationSpeed       = pm->value<FloatType>(   "propagation_speed"       ,   100.0,   10000.0);
	config_.acquisitionDelay       = pm->value<FloatType>(   "acquisition_delay"       ,     0.0,       1.0);
	config_.samplingFrequency      = pm->value<FloatType>(   "sampling_frequency"      ,   100.0,   200.0e6);
	config_.focalEmissionDistance  = pm->value<FloatType>(   "focus_emission_distance" ,  1.0e-3, 1000.0e-3);
	config_.focalReceptionDistance = pm->value<FloatType>(   "focus_reception_distance",  1.0e-3, 1000.0e-3);

	// Sectorial scan grid.
	config_.rangeStart             = pm->value<FloatType>(   "range_start"             ,  1.0e-3, 1000.0e-3);
	config_.rangeEnd               = pm->value<FloatType>(   "range_end"               ,  1.0e-3, 1000.0e-3);
	config_.startAngle             = pm->value<FloatType>(   "start_angle"             ,   -90.0,      89.0);
	config_.endAngle               = pm->value<FloatType>(   "end_angle"               , config_.startAngle + 1.0, 90.0);
	config_.angleStep              = pm->value<FloatType>(   "angle_step"              ,  1.0e-3,      10.0);

	config_.enableFocusing         = pm->value<bool>(        "enable_focusing");

	if (config_.numElementsMux & 0x01) {
		THROW_EXCEPTION(InvalidParameterException, "The value of num_elements_mux is not even.");
	}
	if (config_.numElements & 0x01) {
		THROW_EXCEPTION(InvalidParameterException, "The value of num_elements is not even.");
	}
}

template<typename FloatType>
SectorialScanMethod<FloatType>::SectorialScanMethod(Project& project)
		: project_(project)
{
	fillConfiguration();
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
	case Method::SECTORIAL_SCAN_SP_NETWORK:
		getSingleImageFromNetwork();
		break;
	case Method::SECTORIAL_SCAN_SP_SAVED:
		showSavedImage();
		break;
	case Method::SECTORIAL_SCAN_SP_NETWORK_CONTINUOUS:
		execContinuousNetworkImaging();
		break;
	default:
		THROW_EXCEPTION(InvalidParameterException, "Invalid method: " << project_.method() << '.');
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
				true, Figure::VISUALIZATION_RECTIFIED_LOG, Figure::COLORMAP_VIRIDIS);

	{
		Matrix2<FloatType> aux;
		Util::copyValueToSimpleMatrix(imageData, aux);
		LOG_DEBUG << "Saving the image...";
		project_.saveHDF5(aux, outputDir + "/image_value", "value");

		Matrix2<FloatType> gridX, gridZ;
		Util::copyXZToSimpleMatrices(imageData, gridX, gridZ);
		LOG_DEBUG << "Saving the X coordinates...";
		project_.saveHDF5(gridX, outputDir + "/image_x", "x");
		LOG_DEBUG << "Saving the Z coordinates...";
		project_.saveHDF5(gridZ, outputDir + "/image_z", "z");
	}
}

template<typename FloatType>
void
SectorialScanMethod<FloatType>::showSavedImage()
{
	const std::string outputDir = project_.taskParameterMap()->value<std::string>("output_dir");

	Matrix2<XZValue<FloatType>> imageData;

	{
		Matrix2<FloatType> aux;
		LOG_DEBUG << "Loading the image...";
		project_.loadHDF5(outputDir + "/image_value", "value", aux);
		imageData.resize(aux.n1(), aux.n2());
		Util::copyValueFromSimpleMatrix(aux, imageData);

		Matrix2<FloatType> gridX, gridZ;
		LOG_DEBUG << "Loading the X coordinates...";
		project_.loadHDF5(outputDir + "/image_x", "x", gridX);
		LOG_DEBUG << "Loading the Z coordinates...";
		project_.loadHDF5(outputDir + "/image_z", "z", gridZ);
		Util::copyXZFromSimpleMatrices(gridX, gridZ, imageData);
	}

	std::vector<XZ<float>> pointList = {{((config_.numElements - 1U) / 2.0f) * config_.pitch, 0.0}};

	project_.showFigure3D(1, "Image", &imageData, &pointList,
				true, Figure::VISUALIZATION_RECTIFIED_LOG, Figure::COLORMAP_VIRIDIS);
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
					true, Figure::VISUALIZATION_RECTIFIED_LOG, Figure::COLORMAP_VIRIDIS);

		if (++n == 10) {
			LOG_INFO << 10.0 / t.getTime() << " image/s";
			n = 0;
			t.reset();
		}
	} while (!project_.processingCancellationRequested());
}

} // namespace Lab

#endif // SECTORIALSCANMETHOD_H
