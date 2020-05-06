/***************************************************************************
 *  Copyright 2019, 2020 Marcelo Y. Matuda                                 *
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

#ifndef SIMCIRCULARSOURCEMETHOD_H
#define SIMCIRCULARSOURCEMETHOD_H

#include <cmath>
#include <iomanip>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "AnalyticCircularSourceImpulseResponse.h"
#include "Colormap.h"
#include "Exception.h"
#include "FFTWFilter2.h"
#include "ImageGrid.h"
#include "Log.h"
#include "Matrix.h"
#include "Method.h"
#include "NumericCircularSourceImpulseResponse.h"
#include "ParameterMap.h"
#include "Project.h"
#include "SimTransientAcousticField.h"
#include "SimTransientPropagation.h"
#include "SimTransientRadiationPattern.h"
#include "Timer.h"
#include "Util.h"
#include "Visualization.h"
#include "Waveform.h"
#include "XYZ.h"
#include "XYZValue.h"
#include "XYZValueArray.h"
#ifdef USE_CUDA
# include "NumericCircularSourceCUDAImpulseResponse.h"
#endif
#ifdef USE_OPENCL
# include "NumericCircularSourceOCLImpulseResponse.h"
#endif



namespace Lab {

template<typename TFloat>
class SimCircularSourceMethod : public Method {
public:
	SimCircularSourceMethod(Project& project);
	virtual ~SimCircularSourceMethod() = default;

	virtual void execute();
private:
	struct MainData {
		TFloat propagationSpeed;
		TFloat centerFreq;
		TFloat maxFreq;
		TFloat nyquistRate;
		std::string outputDir;
	};
	struct SimulationData {
		TFloat samplingFreq;
		TFloat excNumPeriods;
		TFloat discretFactor;
		std::string irMethod;
		std::string excitationType;
		std::vector<TFloat> exc;
	};
	struct SourceData {
		TFloat sourceRadius;
	};

	SimCircularSourceMethod(const SimCircularSourceMethod&) = delete;
	SimCircularSourceMethod& operator=(const SimCircularSourceMethod&) = delete;
	SimCircularSourceMethod(SimCircularSourceMethod&&) = delete;
	SimCircularSourceMethod& operator=(SimCircularSourceMethod&&) = delete;

	void loadData(const ParameterMap& taskPM, MainData& data, SimulationData& simData);
	void loadSourceData(SourceData& srcData);
	void loadSimulationData(const MainData& data, const std::string& irMethod, SimulationData& simData);
	void prepareExcitation(TFloat dt, const SimulationData& simData, std::vector<TFloat>& tExc,
				std::vector<TFloat>& dvdt, std::vector<TFloat>& tDvdt);

	void execImpulseResponse(); // calculate p/(c*density)
	void execTransientAcousticField();
	void execTransientPropagation();
	void execTransientRadiationPattern();

	Project& project_;
};

template<typename TFloat>
SimCircularSourceMethod<TFloat>::SimCircularSourceMethod(Project& project)
		: project_(project)
{
}

template<typename TFloat>
void
SimCircularSourceMethod<TFloat>::loadData(const ParameterMap& taskPM, MainData& data, SimulationData& simData)
{
	const ParamMapPtr mainPM = project_.getSubParamMap("main_config_file");
	//mainPM->getValue(data.density         , "density"          , 0.0, 100000.0);
	mainPM->getValue(data.propagationSpeed, "propagation_speed", 0.0, 100000.0);
	mainPM->getValue(data.centerFreq      , "center_frequency" , 0.0,  100.0e6);
	mainPM->getValue(data.maxFreq         , "max_frequency"    , 0.0,  200.0e6);
	data.nyquistRate = Util::nyquistRate(data.maxFreq);
	taskPM.getValue(data.outputDir, "output_dir");

	loadSimulationData(data, taskPM.value<std::string>("impulse_response_method"), simData);
}

template<typename TFloat>
void
SimCircularSourceMethod<TFloat>::loadSourceData(SourceData& srcData)
{
	const ParamMapPtr singlePM = project_.getSubParamMap("source_config_file");
	singlePM->getValue(srcData.sourceRadius, "source_radius" , 0.0, 10.0);
}

template<typename TFloat>
void
SimCircularSourceMethod<TFloat>::loadSimulationData(const MainData& data, const std::string& irMethod, SimulationData& simData)
{
	const ParamMapPtr simPM = project_.getSubParamMap("simulation_config_file");
	simData.samplingFreq = simPM->value<TFloat>("sampling_frequency_factor", 0.0, 10000.0) * data.nyquistRate;
	simPM->getValue(simData.excitationType, "excitation_type");
	simPM->getValue(simData.excNumPeriods , "excitation_num_periods", 0.0, 100.0);
	simData.irMethod = irMethod;

	Waveform::get(simData.excitationType, data.centerFreq, simData.samplingFreq, simData.excNumPeriods, simData.exc);

	if (simData.irMethod == "numeric" ||
			simData.irMethod == "numeric_cuda" ||
			simData.irMethod == "numeric_ocl") {
		simPM->getValue(simData.discretFactor , "num_sub_elem_per_lambda", 0.01, 100.0);
	} else if (simData.irMethod == "analytic") {
		// Empty.
	} else {
		THROW_EXCEPTION(InvalidParameterException, "Invalid impulse response method: " << simData.irMethod << '.');
	}
}

template<typename TFloat>
void
SimCircularSourceMethod<TFloat>::prepareExcitation(TFloat dt,
						const SimulationData& simData, std::vector<TFloat>& tExc,
						std::vector<TFloat>& dvdt, std::vector<TFloat>& tDvdt)
{
	Util::fillSequenceFromStartWithStepAndSize(tExc, TFloat(0.0), dt, simData.exc.size());
	project_.showFigure2D(1, "v", tExc, simData.exc);

	Util::centralDiff(simData.exc, dt, dvdt);
	Util::normalizeBySumOfAbs(dvdt);
	Util::fillSequenceFromStartWithStepAndSize(tDvdt, TFloat(0.0), dt, dvdt.size());
	project_.showFigure2D(2, "dv/dt", tDvdt, dvdt);
}

template<typename TFloat>
void
SimCircularSourceMethod<TFloat>::execTransientRadiationPattern()
{
	const ParameterMap& taskPM = project_.taskParamMap();
	MainData mainData;
	SimulationData simData;
	loadData(taskPM, mainData, simData);
	SourceData srcData;
	loadSourceData(srcData);

	project_.createDirectory(mainData.outputDir, false);

	const ParamMapPtr radPM = project_.getSubParamMap("rad_config_file");
	const auto distance   = radPM->value<TFloat>("distance", 0.0, 100.0);
	const auto thetaYStep = radPM->value<TFloat>("theta_y_step", 0.0, 10.0);
	const auto thetaYMin  = TFloat(0);
	const auto thetaYMax  = radPM->value<TFloat>("theta_y_max" , thetaYMin + 0.1, 90.0);

	const TFloat dt = 1.0 / simData.samplingFreq;
	std::vector<TFloat> tExc, dvdt, tDvdt;
	prepareExcitation(dt, simData, tExc, dvdt, tDvdt);

	std::vector<TFloat> thetaYList;
	Util::fillSequenceFromStartToEndWithSize(thetaYList, thetaYMin, thetaYMax, std::ceil((thetaYMax - thetaYMin) / thetaYStep) + 1);

	std::vector<TFloat> radData(thetaYList.size());
	std::vector<XYZ<TFloat>> inputData(thetaYList.size());

	for (unsigned int i = 0, size = thetaYList.size(); i < size; ++i) {
		const TFloat tY = Util::degreeToRadian(thetaYList[i]);
		const TFloat x = distance * std::sin(tY);
		const TFloat y = 0.0;
		const TFloat z = distance * std::cos(tY);
		XYZ<TFloat>& id = inputData[i];
		id.x = x;
		id.y = y;
		id.z = z;
	}

	if (simData.irMethod == "numeric") {
		const TFloat nyquistLambda = Util::wavelength(mainData.propagationSpeed, mainData.nyquistRate);
		const TFloat numSubElemPerLambda = simData.discretFactor;
		const TFloat numSubElemInRadius = srcData.sourceRadius * (numSubElemPerLambda / nyquistLambda);
		SimTransientRadiationPattern<
			TFloat,
			NumericCircularSourceImpulseResponse<TFloat>>::getCircularSourceRadiationPattern(
					simData.samplingFreq, mainData.propagationSpeed, srcData.sourceRadius,
					numSubElemInRadius,
					dvdt, inputData, radData);
#ifdef USE_CUDA
	} else if (simData.irMethod == "numeric_cuda") {
		if constexpr (std::is_same<TFloat, float>::value) {
			const TFloat nyquistLambda = Util::wavelength(mainData.propagationSpeed, mainData.nyquistRate);
			const TFloat numSubElemPerLambda = simData.discretFactor;
			const TFloat numSubElemInRadius = srcData.sourceRadius * (numSubElemPerLambda / nyquistLambda);
			SimTransientRadiationPattern<
				TFloat,
				NumericCircularSourceCUDAImpulseResponse>::getCircularSourceRadiationPatternSingleThread(
						simData.samplingFreq, mainData.propagationSpeed, srcData.sourceRadius,
						numSubElemInRadius,
						dvdt, inputData, radData);
		} else {
			THROW_EXCEPTION(InvalidValueException, "Invalid float type.");
		}
#endif
#ifdef USE_OPENCL
	} else if (simData.irMethod == "numeric_ocl") {
		const TFloat nyquistLambda = Util::wavelength(mainData.propagationSpeed, mainData.nyquistRate);
		const TFloat numSubElemPerLambda = simData.discretFactor;
		const TFloat numSubElemInRadius = srcData.sourceRadius * (numSubElemPerLambda / nyquistLambda);
		SimTransientRadiationPattern<
			TFloat,
			NumericCircularSourceOCLImpulseResponse<TFloat>>::getCircularSourceRadiationPatternSingleThread(
					simData.samplingFreq, mainData.propagationSpeed, srcData.sourceRadius,
					numSubElemInRadius,
					dvdt, inputData, radData);
#endif
	} else if (simData.irMethod == "analytic") {
		SimTransientRadiationPattern<
			TFloat,
			AnalyticCircularSourceImpulseResponse<TFloat>>::getCircularSourceRadiationPattern(
					simData.samplingFreq, mainData.propagationSpeed, srcData.sourceRadius,
					0.0,
					dvdt, inputData, radData);
	} else {
		THROW_EXCEPTION(InvalidParameterException, "Invalid impulse response method: " << simData.irMethod << '.');
	}

	const TFloat maxAbsValue = Util::maxAbsolute(radData);
	const TFloat k = 1.0 / maxAbsValue;
	Util::multiply(radData, k);

	Util::linearToDecibels(radData, TFloat(-100.0));
	project_.showFigure2D(3, "Pattern", thetaYList, radData);

	project_.saveHDF5(simData.exc, mainData.outputDir + "/v"              , "value");
	project_.saveHDF5(tExc       , mainData.outputDir + "/v_time"         , "value");
	project_.saveHDF5(dvdt       , mainData.outputDir + "/dvdt"           , "value");
	project_.saveHDF5(tDvdt      , mainData.outputDir + "/dvdt_time"      , "value");
	project_.saveHDF5(thetaYList , mainData.outputDir + "/theta_y"        , "value");
	project_.saveHDF5(radData    , mainData.outputDir + "/pattern_theta_y", "value");
}

template<typename TFloat>
void
SimCircularSourceMethod<TFloat>::execTransientAcousticField()
{
	const ParameterMap& taskPM = project_.taskParamMap();
	MainData mainData;
	SimulationData simData;
	loadData(taskPM, mainData, simData);
	SourceData srcData;
	loadSourceData(srcData);

	project_.createDirectory(mainData.outputDir, false);

	const TFloat dt = 1.0 / simData.samplingFreq;
	std::vector<TFloat> tExc, dvdt, tDvdt;
	prepareExcitation(dt, simData, tExc, dvdt, tDvdt);

	Matrix<XYZValue<TFloat>> gridData;

	const TFloat nyquistLambda = Util::wavelength(mainData.propagationSpeed, mainData.nyquistRate);
	ImageGrid<TFloat>::get(*project_.getSubParamMap("grid_config_file"), nyquistLambda, gridData);

	if (simData.irMethod == "numeric") {
		const TFloat numSubElemPerLambda = simData.discretFactor;
		const TFloat numSubElemInRadius = srcData.sourceRadius * (numSubElemPerLambda / nyquistLambda);
		SimTransientAcousticField<TFloat>::template getCircularSourceAcousticField<NumericCircularSourceImpulseResponse<TFloat>>(
					simData.samplingFreq, mainData.propagationSpeed, srcData.sourceRadius,
					numSubElemInRadius,
					dvdt, gridData);
#ifdef USE_CUDA
	} else if (simData.irMethod == "numeric_cuda") {
		if constexpr (std::is_same<TFloat, float>::value) {
			const TFloat numSubElemPerLambda = simData.discretFactor;
			const TFloat numSubElemInRadius = srcData.sourceRadius * (numSubElemPerLambda / nyquistLambda);
			SimTransientAcousticField<TFloat>::template getCircularSourceAcousticFieldSingleThread<NumericCircularSourceCUDAImpulseResponse>(
						simData.samplingFreq, mainData.propagationSpeed, srcData.sourceRadius,
						numSubElemInRadius,
						dvdt, gridData);
		} else {
			THROW_EXCEPTION(InvalidValueException, "Invalid float type.");
		}
#endif
#ifdef USE_OPENCL
	} else if (simData.irMethod == "numeric_ocl") {
		const TFloat numSubElemPerLambda = simData.discretFactor;
		const TFloat numSubElemInRadius = srcData.sourceRadius * (numSubElemPerLambda / nyquistLambda);
		SimTransientAcousticField<TFloat>::template getCircularSourceAcousticFieldSingleThread<NumericCircularSourceOCLImpulseResponse<TFloat>>(
					simData.samplingFreq, mainData.propagationSpeed, srcData.sourceRadius,
					numSubElemInRadius,
					dvdt, gridData);
#endif
	} else if (simData.irMethod == "analytic") {
		SimTransientAcousticField<TFloat>::template getCircularSourceAcousticField<AnalyticCircularSourceImpulseResponse<TFloat>>(
					simData.samplingFreq, mainData.propagationSpeed, srcData.sourceRadius,
					0.0,
					dvdt, gridData);
	} else {
		THROW_EXCEPTION(InvalidParameterException, "Invalid impulse response method: " << simData.irMethod << '.');
	}

	const TFloat maxAbsValue = Util::maxAbsoluteValueField<XYZValue<TFloat>, TFloat>(gridData);
	const TFloat k = 1.0 / maxAbsValue;
	for (auto it = gridData.begin(); it != gridData.end(); ++it) {
		it->value *= k;
	}

	std::vector<XYZ<float>> pointList = {{0.0, 0.0, 0.0}};

	project_.showFigure3D(1, "Acoustic field", &gridData, &pointList,
					true, Visualization::VALUE_RECTIFIED_LINEAR, Colormap::GRADIENT_VIRIDIS);

	project_.saveHDF5(simData.exc, mainData.outputDir + "/v"        , "value");
	project_.saveHDF5(tExc       , mainData.outputDir + "/v_time"   , "value");
	project_.saveHDF5(dvdt       , mainData.outputDir + "/dvdt"     , "value");
	project_.saveHDF5(tDvdt      , mainData.outputDir + "/dvdt_time", "value");
	project_.saveImageToHDF5(gridData, mainData.outputDir);
	project_.saveXYZToHDF5(gridData, mainData.outputDir);
}

template<typename TFloat>
void
SimCircularSourceMethod<TFloat>::execTransientPropagation()
{
	Timer tProc;

	const ParameterMap& taskPM = project_.taskParamMap();
	MainData mainData;
	SimulationData simData;
	loadData(taskPM, mainData, simData);
	SourceData srcData;
	loadSourceData(srcData);

	project_.createDirectory(mainData.outputDir, false);

	const ParamMapPtr propagPM = project_.getSubParamMap("propag_config_file");
	const auto propagDistanceStep = propagPM->value<TFloat>("propagation_distance_step", 1.0e-6, 100.0);
	const auto propagMinDistance  = propagPM->value<TFloat>("propagation_min_distance", 0.0, 100.0);
	const auto propagMaxDistance  = propagPM->value<TFloat>("propagation_max_distance", propagMinDistance + propagDistanceStep, 100.0);
	const auto propagRepet        = propagPM->value<unsigned int>("propagation_repetitions", 1, 100);
	const auto propagPause        = propagPM->value<unsigned int>("propagation_pause", 0, 10000);

	std::vector<TFloat> propagDist;
	Util::fillSequenceFromStartWithStep(propagDist, propagMinDistance, propagMaxDistance, propagDistanceStep);
	std::vector<unsigned int> propagIndexList(propagDist.size());
	const TFloat coef = simData.samplingFreq / mainData.propagationSpeed;
	for (unsigned int i = 0, end = propagIndexList.size(); i < end; ++i) {
		propagIndexList[i] = std::rint(propagDist[i] * coef);
	}

	const TFloat dt = 1.0 / simData.samplingFreq;
	std::vector<TFloat> tExc, dvdt, tDvdt;
	prepareExcitation(dt, simData, tExc, dvdt, tDvdt);

	Matrix<XYZValueArray<TFloat>> gridData;

	const TFloat nyquistLambda = Util::wavelength(mainData.propagationSpeed, mainData.nyquistRate);
	ImageGrid<TFloat>::get(*project_.getSubParamMap("grid_config_file"), nyquistLambda, gridData);

	if (simData.irMethod == "numeric") {
		const TFloat numSubElemPerLambda = simData.discretFactor;
		const TFloat numSubElemInRadius = srcData.sourceRadius * (numSubElemPerLambda / nyquistLambda);
		SimTransientPropagation<
			TFloat,
			NumericCircularSourceImpulseResponse<TFloat>>::getCircularSourcePropagation(
					simData.samplingFreq, mainData.propagationSpeed, srcData.sourceRadius,
					numSubElemInRadius,
					dvdt, propagIndexList, gridData);
#ifdef USE_CUDA
	} else if (simData.irMethod == "numeric_cuda") {
		if constexpr (std::is_same<TFloat, float>::value) {
			const TFloat numSubElemPerLambda = simData.discretFactor;
			const TFloat numSubElemInRadius = srcData.sourceRadius * (numSubElemPerLambda / nyquistLambda);
			SimTransientPropagation<
				TFloat,
				NumericCircularSourceCUDAImpulseResponse>::getCircularSourcePropagationSingleThread(
						simData.samplingFreq, mainData.propagationSpeed, srcData.sourceRadius,
						numSubElemInRadius,
						dvdt, propagIndexList, gridData);
		} else {
			THROW_EXCEPTION(InvalidValueException, "Invalid float type.");
		}
#endif
#ifdef USE_OPENCL
	} else if (simData.irMethod == "numeric_ocl") {
		const TFloat numSubElemPerLambda = simData.discretFactor;
		const TFloat numSubElemInRadius = srcData.sourceRadius * (numSubElemPerLambda / nyquistLambda);
		SimTransientPropagation<
			TFloat,
			NumericCircularSourceOCLImpulseResponse<TFloat>>::getCircularSourcePropagationSingleThread(
					simData.samplingFreq, mainData.propagationSpeed, srcData.sourceRadius,
					numSubElemInRadius,
					dvdt, propagIndexList, gridData);
#endif
	} else if (simData.irMethod == "analytic") {
		SimTransientPropagation<
			TFloat,
			AnalyticCircularSourceImpulseResponse<TFloat>>::getCircularSourcePropagation(
					simData.samplingFreq, mainData.propagationSpeed, srcData.sourceRadius,
					0.0,
					dvdt, propagIndexList, gridData);
	} else {
		THROW_EXCEPTION(InvalidParameterException, "Invalid impulse response method: " << simData.irMethod << '.');
	}

	// Normalize.
	const TFloat maxAbsValue = Util::maxAbsoluteValueField(gridData);
	const TFloat k = 1.0 / maxAbsValue;
	for (auto it = gridData.begin(); it != gridData.end(); ++it) {
		for (auto& v : it->values) {
			v *= k;
		}
	}
	LOG_INFO << "gridData: maxAbsValue = " << maxAbsValue;

	LOG_DEBUG << ">>> Processing time: " << tProc.getTime();

	std::vector<XYZ<float>> pointList = {{0.0, 0.0, 0.0}};
	Project::GridDataType projGridData;

	// Show the propagation.
	Util::copyXYZ(gridData, projGridData);
	for (unsigned int n = 0; n < propagRepet; ++n) {
		for (unsigned int i = 0, end = propagIndexList.size(); i < end; ++i) {
			auto origIter = gridData.begin();
			auto destIter = projGridData.begin();
			while (origIter != gridData.end()) {
				destIter->value = origIter->values[i]; // copy just one value
				++origIter; ++destIter;
			}

			project_.showFigure3D(1, "Propagation", &projGridData, &pointList,
						true, Visualization::VALUE_RAW_LINEAR, Colormap::GRADIENT_GRAY);

			Util::sleepMs(propagPause);
		}
	}

	project_.saveHDF5(simData.exc, mainData.outputDir + "/v"        , "value");
	project_.saveHDF5(tExc       , mainData.outputDir + "/v_time"   , "value");
	project_.saveHDF5(dvdt       , mainData.outputDir + "/dvdt"     , "value");
	project_.saveHDF5(tDvdt      , mainData.outputDir + "/dvdt_time", "value");

	// Save images.
	project_.saveXYZToHDF5(gridData, mainData.outputDir);
	std::string imagePrefix = mainData.outputDir + "/image_value-propag-";
	const unsigned int numDigits = Util::numberOfDigits(propagIndexList.size() - 1U);
	for (unsigned int i = 0, end = propagIndexList.size(); i < end; ++i) {
		LOG_DEBUG << "Saving image " << i << "...";
		std::ostringstream imageFileName;
		imageFileName << imagePrefix << std::setw(numDigits) << std::setfill('0') << i;
		project_.saveHDF5(gridData, imageFileName.str(),
					"value",
					[&i](const XYZValueArray<TFloat>& orig, double& dest) {
						dest = orig.values[i];
					});
	}
}

template<typename TFloat>
void
SimCircularSourceMethod<TFloat>::execImpulseResponse()
{
	const ParameterMap& taskPM = project_.taskParamMap();
	MainData mainData;
	SimulationData simData;
	loadData(taskPM, mainData, simData);
	SourceData srcData;
	loadSourceData(srcData);

	project_.createDirectory(mainData.outputDir, false);

	const ParamMapPtr irPM = project_.getSubParamMap("ir_config_file");
	const auto pointX = irPM->value<TFloat>("point_x", 0.0, 10000.0);
	const auto pointZ = irPM->value<TFloat>("point_z", 0.0, 10000.0);

	const TFloat dt = 1.0 / simData.samplingFreq;
	std::vector<TFloat> tExc, dvdt, tDvdt;
	prepareExcitation(dt, simData, tExc, dvdt, tDvdt);

	std::size_t hOffset;
	std::vector<TFloat> h;
	if (simData.irMethod == "numeric") {
		const TFloat nyquistLambda = Util::wavelength(mainData.propagationSpeed, mainData.nyquistRate);
		const TFloat numSubElemPerLambda = simData.discretFactor;
		const TFloat numSubElemInRadius = srcData.sourceRadius * (numSubElemPerLambda / nyquistLambda);
		auto impResp = std::make_unique<NumericCircularSourceImpulseResponse<TFloat>>(
					simData.samplingFreq, mainData.propagationSpeed, srcData.sourceRadius,
					numSubElemInRadius);
		impResp->getImpulseResponse(pointX, 0.0, pointZ, hOffset, h);
#ifdef USE_CUDA
	} else if (simData.irMethod == "numeric_cuda") {
		if constexpr (std::is_same<TFloat, float>::value) {
			const TFloat nyquistLambda = Util::wavelength(mainData.propagationSpeed, mainData.nyquistRate);
			const TFloat numSubElemPerLambda = simData.discretFactor;
			const TFloat numSubElemInRadius = srcData.sourceRadius * (numSubElemPerLambda / nyquistLambda);
			auto impResp = std::make_unique<NumericCircularSourceCUDAImpulseResponse>(
						simData.samplingFreq, mainData.propagationSpeed, srcData.sourceRadius,
						numSubElemInRadius);
			impResp->getImpulseResponse(pointX, 0.0, pointZ, hOffset, h);
		} else {
			THROW_EXCEPTION(InvalidValueException, "Invalid float type.");
		}
#endif
#ifdef USE_OPENCL
	} else if (simData.irMethod == "numeric_ocl") {
		const TFloat nyquistLambda = Util::wavelength(mainData.propagationSpeed, mainData.nyquistRate);
		const TFloat numSubElemPerLambda = simData.discretFactor;
		const TFloat numSubElemInRadius = srcData.sourceRadius * (numSubElemPerLambda / nyquistLambda);
		auto impResp = std::make_unique<NumericCircularSourceOCLImpulseResponse<TFloat>>(
					simData.samplingFreq, mainData.propagationSpeed, srcData.sourceRadius,
					numSubElemInRadius);
		impResp->getImpulseResponse(pointX, 0.0, pointZ, hOffset, h);
#endif
	} else if (simData.irMethod == "analytic") {
		auto impResp = std::make_unique<AnalyticCircularSourceImpulseResponse<TFloat>>(
					simData.samplingFreq, mainData.propagationSpeed, srcData.sourceRadius);
		impResp->getImpulseResponse(pointX, 0.0, pointZ, hOffset, h);
	} else {
		THROW_EXCEPTION(InvalidParameterException, "Invalid impulse response method: " << simData.irMethod << '.');
	}

	std::vector<TFloat> tH;
	Util::fillSequenceFromStartWithStepAndSize(tH, hOffset / simData.samplingFreq, dt, h.size());
	project_.showFigure2D(3, "Impulse response", tH, h);

	std::vector<std::complex<TFloat>> filterFreqCoeff;
	auto filter = std::make_unique<FFTWFilter2<TFloat>>();
	filter->setCoefficients(dvdt, filterFreqCoeff);
	std::vector<TFloat> signal;

	filter->filter(filterFreqCoeff, h, signal);

	std::vector<TFloat> tSignal;
	Util::fillSequenceFromStartWithStepAndSize(tSignal, hOffset / simData.samplingFreq, dt, signal.size());
	//Util::multiply(signal, density);
	project_.showFigure2D(4, "Pressure", tSignal, signal);

	const TFloat maxAbsValue = Util::maxAbsolute(signal);
	LOG_INFO << "signal: maxAbsValue = " << maxAbsValue;

	project_.saveHDF5(simData.exc, mainData.outputDir + "/v"                    , "value");
	project_.saveHDF5(tExc       , mainData.outputDir + "/v_time"               , "value");
	project_.saveHDF5(dvdt       , mainData.outputDir + "/dvdt"                 , "value");
	project_.saveHDF5(tDvdt      , mainData.outputDir + "/dvdt_time"            , "value");
	project_.saveHDF5(h          , mainData.outputDir + "/impulse_response"     , "value");
	project_.saveHDF5(tH         , mainData.outputDir + "/impulse_response_time", "value");
	project_.saveHDF5(signal     , mainData.outputDir + "/pressure"             , "value");
	project_.saveHDF5(tSignal    , mainData.outputDir + "/pressure_time"        , "value");
}

template<typename TFloat>
void
SimCircularSourceMethod<TFloat>::execute()
{
	switch (project_.method()) {
	case MethodEnum::sim_acoustic_field_circular_source_transient:
	case MethodEnum::sim_acoustic_field_circular_source_transient_sp:
		execTransientAcousticField();
		break;
	case MethodEnum::sim_impulse_response_circular_source:
	case MethodEnum::sim_impulse_response_circular_source_sp:
		execImpulseResponse();
		break;
	case MethodEnum::sim_propagation_circular_source_transient:
	case MethodEnum::sim_propagation_circular_source_transient_sp:
		execTransientPropagation();
		break;
	case MethodEnum::sim_radiation_pattern_circular_source_transient:
	case MethodEnum::sim_radiation_pattern_circular_source_transient_sp:
		execTransientRadiationPattern();
		break;
	default:
		THROW_EXCEPTION(InvalidParameterException, "Invalid method: " << static_cast<int>(project_.method()) << '.');
	}
}

} // namespace Lab

#endif // SIMCIRCULARSOURCEMETHOD_H
