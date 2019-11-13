/***************************************************************************
 *  Copyright 2019 Marcelo Y. Matuda                                       *
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

namespace Lab {

template<typename FloatType>
class SimCircularSourceMethod : public Method {
public:
	SimCircularSourceMethod(Project& project);
	virtual ~SimCircularSourceMethod();

	virtual void execute();
private:
	struct MainData {
		FloatType propagationSpeed;
		FloatType centerFreq;
		FloatType maxFreq;
		FloatType nyquistRate;
		std::string outputDir;
	};
	struct SimulationData {
		FloatType samplingFreq;
		FloatType excNumPeriods;
		FloatType discretFactor;
		std::string irMethod;
		std::string excitationType;
		std::vector<FloatType> exc;
	};
	struct SourceData {
		FloatType sourceRadius;
	};

	SimCircularSourceMethod(const SimCircularSourceMethod&) = delete;
	SimCircularSourceMethod& operator=(const SimCircularSourceMethod&) = delete;

	void loadData(ParamMapPtr taskPM, MainData& data, SimulationData& simData);
	void loadSourceData(ParamMapPtr taskPM, SourceData& srcData);
	void loadSimulationData(ParamMapPtr taskPM, const MainData& data, SimulationData& simData);
	void prepareExcitation(FloatType dt, const SimulationData& simData, std::vector<FloatType>& tExc,
				std::vector<FloatType>& dvdt, std::vector<FloatType>& tDvdt);

	// Returns p/(c*density).
	void execImpulseResponse();

	void execTransientAcousticField();
	void execTransientPropagation();
	void execTransientRadiationPattern();

	Project& project_;
};

template<typename FloatType>
SimCircularSourceMethod<FloatType>::SimCircularSourceMethod(Project& project)
		: project_(project)
{
}

template<typename FloatType>
SimCircularSourceMethod<FloatType>::~SimCircularSourceMethod()
{
}

template<typename FloatType>
void
SimCircularSourceMethod<FloatType>::loadData(ParamMapPtr taskPM, MainData& data, SimulationData& simData)
{
	ParamMapPtr mainPM = project_.loadChildParameterMap(taskPM, "main_config_file");
	//mainPM->getValue(data.density         , "density"          , 0.0, 100000.0);
	mainPM->getValue(data.propagationSpeed, "propagation_speed", 0.0, 100000.0);
	mainPM->getValue(data.centerFreq      , "center_frequency" , 0.0,  100.0e6);
	mainPM->getValue(data.maxFreq         , "max_frequency"    , 0.0,  200.0e6);
	data.nyquistRate = 2.0 * data.maxFreq;
	taskPM->getValue(data.outputDir, "output_dir");

	loadSimulationData(taskPM, data, simData);
}

template<typename FloatType>
void
SimCircularSourceMethod<FloatType>::loadSourceData(ParamMapPtr taskPM, SourceData& srcData)
{
	ParamMapPtr singlePM = project_.loadChildParameterMap(taskPM, "source_config_file");
	singlePM->getValue(srcData.sourceRadius, "source_radius" , 0.0, 10.0);
}

template<typename FloatType>
void
SimCircularSourceMethod<FloatType>::loadSimulationData(ParamMapPtr taskPM, const MainData& data, SimulationData& simData)
{
	ParamMapPtr simPM = project_.loadChildParameterMap(taskPM, "simulation_config_file");
	simData.samplingFreq = simPM->value<FloatType>("sampling_frequency_factor", 0.0, 10000.0) * data.nyquistRate;
	simPM->getValue(simData.irMethod      , "impulse_response_method");
	simPM->getValue(simData.excitationType, "excitation_type");
	simPM->getValue(simData.excNumPeriods , "excitation_num_periods", 0.0, 100.0);

	Waveform::get(simData.excitationType, data.centerFreq, simData.samplingFreq, simData.excNumPeriods, simData.exc);

	if (simData.irMethod == "numeric") {
		simPM->getValue(simData.discretFactor , "num_sub_elem_per_lambda", 0.01, 100.0);
	} else if (simData.irMethod == "analytic") {
		// Empty.
	} else {
		THROW_EXCEPTION(InvalidParameterException, "Invalid impulse response method: " << simData.irMethod << '.');
	}
}

template<typename FloatType>
void
SimCircularSourceMethod<FloatType>::prepareExcitation(FloatType dt,
						const SimulationData& simData, std::vector<FloatType>& tExc,
						std::vector<FloatType>& dvdt, std::vector<FloatType>& tDvdt)
{
	Util::fillSequenceFromStartWithStepAndSize(tExc, 0.0, dt, simData.exc.size());
	project_.showFigure2D(1, "v", tExc, simData.exc);

	Util::centralDiff(simData.exc, dt, dvdt);
	Util::normalizeBySumOfAbs(dvdt);
	Util::fillSequenceFromStartWithStepAndSize(tDvdt, 0.0, dt, dvdt.size());
	project_.showFigure2D(2, "dv/dt", tDvdt, dvdt);
}

template<typename FloatType>
void
SimCircularSourceMethod<FloatType>::execTransientRadiationPattern()
{
	ParamMapPtr taskPM = project_.taskParameterMap();
	MainData mainData;
	SimulationData simData;
	loadData(taskPM, mainData, simData);
	SourceData srcData;
	loadSourceData(taskPM, srcData);

	project_.createDirectory(mainData.outputDir, false);

	ParamMapPtr radPM = project_.loadChildParameterMap(taskPM, "rad_config_file");
	const auto distance   = radPM->value<FloatType>("distance", 0.0, 100.0);
	const auto thetaYStep = radPM->value<FloatType>("theta_y_step", 0.0, 10.0);
	const auto thetaYMin  = FloatType(0);
	const auto thetaYMax  = radPM->value<FloatType>("theta_y_max" , thetaYMin + 0.1, 90.0);

	const FloatType dt = 1.0 / simData.samplingFreq;
	std::vector<FloatType> tExc, dvdt, tDvdt;
	prepareExcitation(dt, simData, tExc, dvdt, tDvdt);

	std::vector<FloatType> thetaYList;
	Util::fillSequenceFromStartToEndWithSize(thetaYList, thetaYMin, thetaYMax, std::ceil((thetaYMax - thetaYMin) / thetaYStep) + 1);

	std::vector<FloatType> radData(thetaYList.size());
	std::vector<XYZ<FloatType>> inputData(thetaYList.size());

	for (unsigned int i = 0, size = thetaYList.size(); i < size; ++i) {
		const FloatType tY = Util::degreeToRadian(thetaYList[i]);
		const FloatType x = distance * std::sin(tY);
		const FloatType y = 0.0;
		const FloatType z = distance * std::cos(tY);
		XYZ<FloatType>& id = inputData[i];
		id.x = x;
		id.y = y;
		id.z = z;
	}

	if (simData.irMethod == "numeric") {
		const FloatType nyquistLambda = mainData.propagationSpeed / mainData.nyquistRate;
		const FloatType numSubElemPerLambda = simData.discretFactor;
		const FloatType numSubElemInRadius = srcData.sourceRadius * (numSubElemPerLambda / nyquistLambda);
		auto radPat = std::make_unique<SimTransientRadiationPattern<FloatType, NumericCircularSourceImpulseResponse<FloatType>>>();
		radPat->getCircularSourceRadiationPattern(
					simData.samplingFreq, mainData.propagationSpeed, srcData.sourceRadius,
					numSubElemInRadius,
					dvdt, inputData, radData);
	} else if (simData.irMethod == "analytic") {
		auto radPat = std::make_unique<SimTransientRadiationPattern<FloatType, AnalyticCircularSourceImpulseResponse<FloatType>>>();
		radPat->getCircularSourceRadiationPattern(
					simData.samplingFreq, mainData.propagationSpeed, srcData.sourceRadius,
					0.0,
					dvdt, inputData, radData);
	} else {
		THROW_EXCEPTION(InvalidParameterException, "Invalid impulse response method: " << simData.irMethod << '.');
	}

	const FloatType maxAbsValue = Util::maxAbsolute(radData);
	const FloatType k = 1.0 / maxAbsValue;
	Util::multiply(radData, k);

	Util::linearToDecibels(radData, -100.0);
	project_.showFigure2D(3, "Pattern", thetaYList, radData);

	project_.saveHDF5(simData.exc, mainData.outputDir + "/v"              , "value");
	project_.saveHDF5(tExc       , mainData.outputDir + "/v_time"         , "value");
	project_.saveHDF5(dvdt       , mainData.outputDir + "/dvdt"           , "value");
	project_.saveHDF5(tDvdt      , mainData.outputDir + "/dvdt_time"      , "value");
	project_.saveHDF5(thetaYList , mainData.outputDir + "/theta_y"        , "value");
	project_.saveHDF5(radData    , mainData.outputDir + "/pattern_theta_y", "value");
}

template<typename FloatType>
void
SimCircularSourceMethod<FloatType>::execTransientAcousticField()
{
	ParamMapPtr taskPM = project_.taskParameterMap();
	MainData mainData;
	SimulationData simData;
	loadData(taskPM, mainData, simData);
	SourceData srcData;
	loadSourceData(taskPM, srcData);

	project_.createDirectory(mainData.outputDir, false);

	const FloatType dt = 1.0 / simData.samplingFreq;
	std::vector<FloatType> tExc, dvdt, tDvdt;
	prepareExcitation(dt, simData, tExc, dvdt, tDvdt);

	Matrix<XYZValue<FloatType>> gridData;

	const FloatType nyquistLambda = mainData.propagationSpeed / mainData.nyquistRate;
	ImageGrid<FloatType>::get(project_.loadChildParameterMap(taskPM, "grid_config_file"), nyquistLambda, gridData);

	if (simData.irMethod == "numeric") {
		const FloatType numSubElemPerLambda = simData.discretFactor;
		const FloatType numSubElemInRadius = srcData.sourceRadius * (numSubElemPerLambda / nyquistLambda);
		auto acField = std::make_unique<SimTransientAcousticField<FloatType, NumericCircularSourceImpulseResponse<FloatType>>>();
		acField->getCircularSourceAcousticField(
					simData.samplingFreq, mainData.propagationSpeed, srcData.sourceRadius,
					numSubElemInRadius,
					dvdt, gridData);
	} else if (simData.irMethod == "analytic") {
		auto acField = std::make_unique<SimTransientAcousticField<FloatType, AnalyticCircularSourceImpulseResponse<FloatType>>>();
		acField->getCircularSourceAcousticField(
					simData.samplingFreq, mainData.propagationSpeed, srcData.sourceRadius,
					0.0,
					dvdt, gridData);
	} else {
		THROW_EXCEPTION(InvalidParameterException, "Invalid impulse response method: " << simData.irMethod << '.');
	}

	const FloatType maxAbsValue = Util::maxAbsoluteValueField<XYZValue<FloatType>, FloatType>(gridData);
	const FloatType k = 1.0 / maxAbsValue;
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

template<typename FloatType>
void
SimCircularSourceMethod<FloatType>::execTransientPropagation()
{
	ParamMapPtr taskPM = project_.taskParameterMap();
	MainData mainData;
	SimulationData simData;
	loadData(taskPM, mainData, simData);
	SourceData srcData;
	loadSourceData(taskPM, srcData);

	project_.createDirectory(mainData.outputDir, false);

	ParamMapPtr propagPM = project_.loadChildParameterMap(taskPM, "propag_config_file");
	const auto propagDistanceStep = propagPM->value<FloatType>("propagation_distance_step", 1.0e-6, 100.0);
	const auto propagMinDistance  = propagPM->value<FloatType>("propagation_min_distance", 0.0, 100.0);
	const auto propagMaxDistance  = propagPM->value<FloatType>("propagation_max_distance", propagMinDistance + propagDistanceStep, 100.0);
	const auto propagRepet        = propagPM->value<unsigned int>("propagation_repetitions", 1, 100);
	const auto propagPause        = propagPM->value<unsigned int>("propagation_pause", 0, 10000);

	std::vector<FloatType> propagDist;
	Util::fillSequenceFromStartWithStep(propagDist, propagMinDistance, propagMaxDistance, propagDistanceStep);
	std::vector<unsigned int> propagIndexList(propagDist.size());
	const FloatType coef = simData.samplingFreq / mainData.propagationSpeed;
	for (unsigned int i = 0, end = propagIndexList.size(); i < end; ++i) {
		propagIndexList[i] = std::rint(propagDist[i] * coef);
	}

	const FloatType dt = 1.0 / simData.samplingFreq;
	std::vector<FloatType> tExc, dvdt, tDvdt;
	prepareExcitation(dt, simData, tExc, dvdt, tDvdt);

	Matrix<XYZValueArray<FloatType>> gridData;

	const FloatType nyquistLambda = mainData.propagationSpeed / mainData.nyquistRate;
	ImageGrid<FloatType>::get(project_.loadChildParameterMap(taskPM, "grid_config_file"), nyquistLambda, gridData);

	if (simData.irMethod == "numeric") {
		const FloatType numSubElemPerLambda = simData.discretFactor;
		const FloatType numSubElemInRadius = srcData.sourceRadius * (numSubElemPerLambda / nyquistLambda);
		auto propag = std::make_unique<SimTransientPropagation<FloatType, NumericCircularSourceImpulseResponse<FloatType>>>();
		propag->getCircularSourcePropagation(
					simData.samplingFreq, mainData.propagationSpeed, srcData.sourceRadius,
					numSubElemInRadius,
					dvdt, propagIndexList, gridData);
	} else if (simData.irMethod == "analytic") {
		auto propag = std::make_unique<SimTransientPropagation<FloatType, AnalyticCircularSourceImpulseResponse<FloatType>>>();
		propag->getCircularSourcePropagation(
					simData.samplingFreq, mainData.propagationSpeed, srcData.sourceRadius,
					0.0,
					dvdt, propagIndexList, gridData);
	} else {
		THROW_EXCEPTION(InvalidParameterException, "Invalid impulse response method: " << simData.irMethod << '.');
	}

	// Normalize.
	const FloatType maxAbsValue = Util::maxAbsoluteValueField(gridData);
	const FloatType k = 1.0 / maxAbsValue;
	for (auto it = gridData.begin(); it != gridData.end(); ++it) {
		for (auto& v : it->values) {
			v *= k;
		}
	}
	LOG_INFO << "gridData: maxAbsValue = " << maxAbsValue;

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
					[&i](const XYZValueArray<FloatType>& orig, double& dest) {
						dest = orig.values[i];
					});
	}
}

template<typename FloatType>
void
SimCircularSourceMethod<FloatType>::execImpulseResponse()
{
	ParamMapPtr taskPM = project_.taskParameterMap();
	MainData mainData;
	SimulationData simData;
	loadData(taskPM, mainData, simData);
	SourceData srcData;
	loadSourceData(taskPM, srcData);

	project_.createDirectory(mainData.outputDir, false);

	ParamMapPtr irPM = project_.loadChildParameterMap(taskPM, "ir_config_file");
	const auto pointX = irPM->value<FloatType>("point_x", 0.0, 10000.0);
	const auto pointZ = irPM->value<FloatType>("point_z", 0.0, 10000.0);

	const FloatType dt = 1.0 / simData.samplingFreq;
	std::vector<FloatType> tExc, dvdt, tDvdt;
	prepareExcitation(dt, simData, tExc, dvdt, tDvdt);

	std::size_t hOffset;
	std::vector<FloatType> h;
	if (simData.irMethod == "numeric") {
		const FloatType nyquistLambda = mainData.propagationSpeed / mainData.nyquistRate;
		const FloatType numSubElemPerLambda = simData.discretFactor;
		const FloatType numSubElemInRadius = srcData.sourceRadius * (numSubElemPerLambda / nyquistLambda);
		auto impResp = std::make_unique<NumericCircularSourceImpulseResponse<FloatType>>(
					simData.samplingFreq, mainData.propagationSpeed, srcData.sourceRadius,
					numSubElemInRadius);
		impResp->getImpulseResponse(pointX, 0.0, pointZ, hOffset, h);
	} else if (simData.irMethod == "analytic") {
		auto impResp = std::make_unique<AnalyticCircularSourceImpulseResponse<FloatType>>(
					simData.samplingFreq, mainData.propagationSpeed, srcData.sourceRadius);
		impResp->getImpulseResponse(pointX, 0.0, pointZ, hOffset, h);
	} else {
		THROW_EXCEPTION(InvalidParameterException, "Invalid impulse response method: " << simData.irMethod << '.');
	}

	std::vector<FloatType> tH;
	Util::fillSequenceFromStartWithStepAndSize(tH, hOffset / simData.samplingFreq, dt, h.size());
	project_.showFigure2D(3, "Impulse response", tH, h);

	std::vector<std::complex<FloatType>> filterFreqCoeff;
	auto filter = std::make_unique<FFTWFilter2<FloatType>>();
	filter->setCoefficients(dvdt, filterFreqCoeff);
	std::vector<FloatType> signal;

	filter->filter(filterFreqCoeff, h, signal);

	std::vector<FloatType> tSignal;
	Util::fillSequenceFromStartWithStepAndSize(tSignal, hOffset / simData.samplingFreq, dt, signal.size());
	//Util::multiply(signal, density);
	project_.showFigure2D(4, "Pressure", tSignal, signal);

	const FloatType maxAbsValue = Util::maxAbsolute(signal);
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

template<typename FloatType>
void
SimCircularSourceMethod<FloatType>::execute()
{
	Timer tProc;

	switch (project_.method()) {
	case MethodEnum::sim_acoustic_field_circular_source_transient:
		execTransientAcousticField();
		break;
	case MethodEnum::sim_impulse_response_circular_source:
		execImpulseResponse();
		break;
	case MethodEnum::sim_propagation_circular_source_transient:
		execTransientPropagation();
		break;
	case MethodEnum::sim_radiation_pattern_circular_source_transient:
		execTransientRadiationPattern();
		break;
	default:
		THROW_EXCEPTION(InvalidParameterException, "Invalid method: " << static_cast<int>(project_.method()) << '.');
	}

	LOG_DEBUG << ">>> Processing time: " << tProc.getTime();
}

} // namespace Lab

#endif // SIMCIRCULARSOURCEMETHOD_H
