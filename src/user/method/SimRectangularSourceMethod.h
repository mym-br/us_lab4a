/***************************************************************************
 *  Copyright 2018, 2019 Marcelo Y. Matuda                                 *
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

#ifndef SIMRECTANGULARSOURCEMETHOD_H
#define SIMRECTANGULARSOURCEMETHOD_H

#include <cmath>
#include <iomanip>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "AnalyticRectangularSourceImpulseResponse.h"
#include "ArrayOfRectangularSourcesImpulseResponse.h"
#include "ArrayUtil.h"
#include "Colormap.h"
#include "Exception.h"
#include "FFTWFilter2.h"
#include "ImageGrid.h"
#include "Log.h"
#include "Matrix.h"
#include "Method.h"
#include "NumericRectangularSourceImpulseResponse.h"
#include "ParameterMap.h"
#include "Project.h"
#include "SimTransientAcousticField.h"
#include "SimTransientPropagation.h"
#include "SimTransientRadiationPattern.h"
#include "Timer.h"
#include "Util.h"
#include "Visualization.h"
#include "Waveform.h"
#include "WavefrontObjFileWriter.h"
#include "XY.h"
#include "XYZ.h"
#include "XYZValue.h"
#include "XYZValueArray.h"

namespace Lab {

template<typename FloatType>
class SimRectangularSourceMethod : public Method {
public:
	SimRectangularSourceMethod(Project& project);
	virtual ~SimRectangularSourceMethod() = default;

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
		FloatType sourceWidth;
		FloatType sourceHeight;

		// For arrays.
		FloatType focusX;
		FloatType focusY;
		FloatType focusZ;
		std::vector<XY<FloatType>> elemPos;
		std::vector<FloatType> focusDelay;
		bool useFocus;
	};

	SimRectangularSourceMethod(const SimRectangularSourceMethod&) = delete;
	SimRectangularSourceMethod& operator=(const SimRectangularSourceMethod&) = delete;
	SimRectangularSourceMethod(SimRectangularSourceMethod&&) = delete;
	SimRectangularSourceMethod& operator=(SimRectangularSourceMethod&&) = delete;

	void loadData(const ParameterMap& taskPM, MainData& data, SimulationData& simData);
	void loadSourceData(MainData& data, bool sourceIsArray, SourceData& srcData);
	void loadSimulationData(const MainData& data, SimulationData& simData);
	void prepareExcitation(FloatType dt, const SimulationData& simData, std::vector<FloatType>& tExc,
				std::vector<FloatType>& dvdt, std::vector<FloatType>& tDvdt);

	void execImpulseResponse(bool sourceIsArray); // calculate p/(c*density)
	void execTransientRadiationPattern(bool sourceIsArray);
	void execTransientAcousticField(bool sourceIsArray);
	void execTransientPropagation(bool sourceIsArray);

	Project& project_;
};

template<typename FloatType>
SimRectangularSourceMethod<FloatType>::SimRectangularSourceMethod(Project& project)
		: project_(project)
{
}

template<typename FloatType>
void
SimRectangularSourceMethod<FloatType>::loadData(const ParameterMap& taskPM, MainData& data, SimulationData& simData)
{
	const ParamMapPtr mainPM = project_.getSubParamMap("main_config_file");
	//mainPM->getValue(data.density         , "density"          , 0.0, 100000.0);
	mainPM->getValue(data.propagationSpeed, "propagation_speed", 0.0, 100000.0);
	mainPM->getValue(data.centerFreq      , "center_frequency" , 0.0,  100.0e6);
	mainPM->getValue(data.maxFreq         , "max_frequency"    , 0.0,  200.0e6);
	data.nyquistRate = Util::nyquistRate(data.maxFreq);
	taskPM.getValue(data.outputDir, "output_dir");

	loadSimulationData(data, simData);
}

template<typename FloatType>
void
SimRectangularSourceMethod<FloatType>::loadSourceData(MainData& data, bool sourceIsArray, SourceData& srcData)
{
	WavefrontObjFileWriter<FloatType> fw((project_.expDirectory() + "/source_geometry.obj").c_str());

	if (sourceIsArray) {
		const ParamMapPtr arrayPM = project_.getSubParamMap("array_config_file");
		ArrayUtil::calculateTxElementPositions(*arrayPM, srcData.elemPos);
		arrayPM->getValue(srcData.sourceWidth , "element_width" , 0.0, 10.0);
		arrayPM->getValue(srcData.sourceHeight, "element_height", 0.0, 10.0);

		arrayPM->getValue(srcData.useFocus, "use_tx_focus");
		if (srcData.useFocus) {
			arrayPM->getValue(srcData.focusX, "tx_focus_x", -10000.0, 10000.0);
			arrayPM->getValue(srcData.focusY, "tx_focus_y", -10000.0, 10000.0);
			arrayPM->getValue(srcData.focusZ, "tx_focus_z", -10000.0, 10000.0);
			ArrayUtil::calculateTx3DFocusDelay(srcData.focusX, srcData.focusY, srcData.focusZ,
								data.propagationSpeed, srcData.elemPos, srcData.focusDelay);
		} else {
			srcData.focusX = srcData.focusY = srcData.focusZ = 0;
			srcData.focusDelay.assign(srcData.elemPos.size(), 0.0);
		}

		const FloatType hw = 0.5 * srcData.sourceWidth;
		const FloatType hh = 0.5 * srcData.sourceHeight;
		for (const auto& pos : srcData.elemPos) {
			fw.addPoint(pos.x - hw, pos.y + hh, 0.0);
			fw.addPoint(pos.x + hw, pos.y + hh, 0.0);
			fw.addPoint(pos.x + hw, pos.y - hh, 0.0);
			fw.addPoint(pos.x - hw, pos.y - hh, 0.0);
			fw.addQuad(-4, -3, -2, -1);
		}
		fw.write();
	} else {
		const ParamMapPtr singlePM = project_.getSubParamMap("single_source_config_file");
		singlePM->getValue(srcData.sourceWidth , "source_width" , 0.0, 10.0);
		singlePM->getValue(srcData.sourceHeight, "source_height", 0.0, 10.0);

		srcData.useFocus = false;
		srcData.focusX = srcData.focusY = srcData.focusZ = 0;

		const FloatType hw = 0.5 * srcData.sourceWidth;
		const FloatType hh = 0.5 * srcData.sourceHeight;
		fw.addPoint(-hw,  hh, 0.0);
		fw.addPoint( hw,  hh, 0.0);
		fw.addPoint( hw, -hh, 0.0);
		fw.addPoint(-hw, -hh, 0.0);
		fw.addQuad(-4, -3, -2, -1);
		fw.write();
	}
}

template<typename FloatType>
void
SimRectangularSourceMethod<FloatType>::loadSimulationData(const MainData& data, SimulationData& simData)
{
	const ParamMapPtr simPM = project_.getSubParamMap("simulation_config_file");
	simData.samplingFreq = simPM->value<FloatType>("sampling_frequency_factor", 0.0, 10000.0) * data.nyquistRate;
	simPM->getValue(simData.irMethod      , "impulse_response_method");
	simPM->getValue(simData.excitationType, "excitation_type");
	simPM->getValue(simData.excNumPeriods , "excitation_num_periods", 0.0, 100.0);

	Waveform::get(simData.excitationType, data.centerFreq, simData.samplingFreq, simData.excNumPeriods, simData.exc);

	if (simData.irMethod == "numeric") {
		simPM->getValue(simData.discretFactor, "sub_elem_size_factor", 0.0, 1.0e3);
	} else if (simData.irMethod == "analytic") {
		simPM->getValue(simData.discretFactor, "min_edge_divisor"    , 0.0, 1.0e6);
	} else {
		THROW_EXCEPTION(InvalidParameterException, "Invalid impulse response method: " << simData.irMethod << '.');
	}
}

template<typename FloatType>
void
SimRectangularSourceMethod<FloatType>::prepareExcitation(FloatType dt,
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
SimRectangularSourceMethod<FloatType>::execTransientRadiationPattern(bool sourceIsArray)
{
	const ParameterMap& taskPM = project_.taskParamMap();
	MainData mainData;
	SimulationData simData;
	loadData(taskPM, mainData, simData);
	SourceData srcData;
	loadSourceData(mainData, sourceIsArray, srcData);

	project_.createDirectory(mainData.outputDir, false);

	const ParamMapPtr radPM = project_.getSubParamMap("rad_config_file");
	const auto distance   = radPM->value<FloatType>("distance", 0.0, 100.0);
	const auto thetaYStep = radPM->value<FloatType>("theta_y_step", 0.0, 10.0);
	const auto thetaYMin  = sourceIsArray ? radPM->value<FloatType>("theta_y_min" , -90.0, 90.0) : FloatType(0);
	const auto thetaYMax  = radPM->value<FloatType>("theta_y_max" , thetaYMin + 0.1, 90.0);
	const auto thetaXStep = radPM->value<FloatType>("theta_x_step", 0.0, 10.0);
	const auto thetaXMin  = sourceIsArray ? radPM->value<FloatType>("theta_x_min" , -90.0, 90.0) : FloatType(0);
	const auto thetaXMax  = radPM->value<FloatType>("theta_x_max" , thetaXMin + 0.1, 90.0);

	const FloatType dt = 1.0 / simData.samplingFreq;
	std::vector<FloatType> tExc, dvdt, tDvdt;
	prepareExcitation(dt, simData, tExc, dvdt, tDvdt);

	std::vector<FloatType> thetaXList;
	Util::fillSequenceFromStartToEndWithSize(thetaXList, thetaXMin, thetaXMax, std::ceil((thetaXMax - thetaXMin) / thetaXStep) + 1);
	std::vector<FloatType> thetaYList;
	Util::fillSequenceFromStartToEndWithSize(thetaYList, thetaYMin, thetaYMax, std::ceil((thetaYMax - thetaYMin) / thetaYStep) + 1);

	Matrix<XYZValue<FloatType>> gridData{thetaXList.size(), thetaYList.size()};
	Matrix<XYZ<FloatType>> inputData{thetaXList.size(), thetaYList.size()};

	for (unsigned int ix = 0, xSize = thetaXList.size(); ix < xSize; ++ix) {
		const FloatType tX = Util::degreeToRadian(thetaXList[ix]);
		const FloatType sinTX = std::sin(tX);
		const FloatType cosTX = std::cos(tX);
		for (unsigned int iy = 0, ySize = thetaYList.size(); iy < ySize; ++iy) {
			const FloatType tY = Util::degreeToRadian(thetaYList[iy]);
			const FloatType x = distance * std::sin(tY);
			const FloatType rx = distance * std::cos(tY);
			const FloatType y = rx * sinTX;
			const FloatType z = rx * cosTX;
			XYZValue<FloatType>& gd = gridData(ix, iy);
			gd.x = thetaYList[iy];
			gd.y = 0.0;
			gd.z = thetaXList[ix];
			gd.value = 0.0;
			XYZ<FloatType>& id = inputData(ix, iy);
			id.x = x;
			id.y = y;
			id.z = z;
		}
	}

	if (simData.irMethod == "numeric") {
		const FloatType subElemSize = mainData.propagationSpeed / (mainData.nyquistRate * simData.discretFactor);
		auto radPat = std::make_unique<SimTransientRadiationPattern<FloatType, NumericRectangularSourceImpulseResponse<FloatType>>>();
		if (sourceIsArray) {
			radPat->getArrayOfRectangularSourcesRadiationPattern(
						simData.samplingFreq, mainData.propagationSpeed, srcData.sourceWidth, srcData.sourceHeight,
						subElemSize,
						dvdt, srcData.elemPos, srcData.focusDelay, inputData, gridData);
		} else {
			radPat->getRectangularSourceRadiationPattern(
						simData.samplingFreq, mainData.propagationSpeed, srcData.sourceWidth, srcData.sourceHeight,
						subElemSize,
						dvdt, inputData, gridData);
		}
	} else if (simData.irMethod == "analytic") {
		const FloatType minEdgeDivisor = simData.discretFactor;
		auto radPat = std::make_unique<SimTransientRadiationPattern<FloatType, AnalyticRectangularSourceImpulseResponse<FloatType>>>();
		if (sourceIsArray) {
			radPat->getArrayOfRectangularSourcesRadiationPattern(
						simData.samplingFreq, mainData.propagationSpeed, srcData.sourceWidth, srcData.sourceHeight,
						minEdgeDivisor,
						dvdt, srcData.elemPos, srcData.focusDelay, inputData, gridData);
		} else {
			radPat->getRectangularSourceRadiationPattern(
						simData.samplingFreq, mainData.propagationSpeed, srcData.sourceWidth, srcData.sourceHeight,
						minEdgeDivisor,
						dvdt, inputData, gridData);
		}
	} else {
		THROW_EXCEPTION(InvalidParameterException, "Invalid impulse response method: " << simData.irMethod << '.');
	}

	const FloatType maxAbsValue = Util::maxAbsoluteValueField<XYZValue<FloatType>, FloatType>(gridData);
	const FloatType k = 1.0 / maxAbsValue;
	for (auto it = gridData.begin(); it != gridData.end(); ++it) {
		it->value *= k;
	}

	std::vector<XYZ<float>> pointList = {{0.0, 0.0, 0.0}};

	project_.showFigure3D(1, "Pattern", &gridData, &pointList,
					true, Visualization::VALUE_RECTIFIED_LOG, Colormap::GRADIENT_VIRIDIS);

	std::size_t sectionTYIndex = 0;
	std::size_t sectionTXIndex = 0;
	if (sourceIsArray) {
		FloatType sectionTY = 0.0;
		FloatType sectionTX = 0.0;
		if (srcData.useFocus) {
			sectionTY = Util::radianToDegree(std::atan2(srcData.focusX, srcData.focusZ));
			sectionTX = Util::radianToDegree(std::atan2(srcData.focusY, srcData.focusZ));
		}
		sectionTYIndex = std::rint(((sectionTY - thetaYMin) / (thetaYMax - thetaYMin)) * (thetaYList.size() - 1));
		sectionTXIndex = std::rint(((sectionTX - thetaXMin) / (thetaXMax - thetaXMin)) * (thetaXList.size() - 1));
		LOG_DEBUG << "Section theta-x: " << thetaXList[sectionTXIndex] << " theta-y: " << thetaYList[sectionTYIndex];
	}

	std::vector<FloatType> patternTY(gridData.n2());
	auto tyInterval = gridData.dim2Interval(sectionTXIndex);
	Util::copyUsingOperator(tyInterval.first, tyInterval.second, patternTY.begin(), Util::CopyValueOp{});
	Util::linearToDecibels(patternTY, -100.0);
	project_.showFigure2D(3, "Pattern theta-y", thetaYList, patternTY);

	std::vector<FloatType> patternTX(gridData.n1());
	auto txInterval = gridData.dim1Interval(sectionTYIndex);
	Util::copyUsingOperator(txInterval.first, txInterval.second, patternTX.begin(), Util::CopyValueOp{});
	Util::linearToDecibels(patternTX, -100.0);
	project_.showFigure2D(4, "Pattern theta-x", thetaXList, patternTX);

	project_.saveHDF5(simData.exc, mainData.outputDir + "/v"              , "value");
	project_.saveHDF5(tExc       , mainData.outputDir + "/v_time"         , "value");
	project_.saveHDF5(dvdt       , mainData.outputDir + "/dvdt"           , "value");
	project_.saveHDF5(tDvdt      , mainData.outputDir + "/dvdt_time"      , "value");
	project_.saveImageToHDF5(gridData, mainData.outputDir);
	project_.saveXYZToHDF5(gridData, mainData.outputDir);
	project_.saveHDF5(thetaYList , mainData.outputDir + "/theta_y"        , "value");
	project_.saveHDF5(patternTY  , mainData.outputDir + "/pattern_theta_y", "value");
	project_.saveHDF5(thetaXList , mainData.outputDir + "/theta_x"        , "value");
	project_.saveHDF5(patternTX  , mainData.outputDir + "/pattern_theta_x", "value");
}

template<typename FloatType>
void
SimRectangularSourceMethod<FloatType>::execTransientAcousticField(bool sourceIsArray)
{
	const ParameterMap& taskPM = project_.taskParamMap();
	MainData mainData;
	SimulationData simData;
	loadData(taskPM, mainData, simData);
	SourceData srcData;
	loadSourceData(mainData, sourceIsArray, srcData);

	project_.createDirectory(mainData.outputDir, false);

	const FloatType dt = 1.0 / simData.samplingFreq;
	std::vector<FloatType> tExc, dvdt, tDvdt;
	prepareExcitation(dt, simData, tExc, dvdt, tDvdt);

	Matrix<XYZValue<FloatType>> gridData;

	const FloatType nyquistLambda = Util::wavelength(mainData.propagationSpeed, mainData.nyquistRate);
	ImageGrid<FloatType>::get(*project_.getSubParamMap("grid_config_file"), nyquistLambda, gridData);

	if (simData.irMethod == "numeric") {
		const FloatType subElemSize = mainData.propagationSpeed / (mainData.nyquistRate * simData.discretFactor);
		auto acField = std::make_unique<SimTransientAcousticField<FloatType, NumericRectangularSourceImpulseResponse<FloatType>>>();
		if (sourceIsArray) {
			acField->getArrayOfRectangularSourcesAcousticField(
						simData.samplingFreq, mainData.propagationSpeed, srcData.sourceWidth, srcData.sourceHeight,
						subElemSize,
						dvdt, srcData.elemPos, srcData.focusDelay, gridData);
		} else {
			acField->getRectangularSourceAcousticField(
						simData.samplingFreq, mainData.propagationSpeed, srcData.sourceWidth, srcData.sourceHeight,
						subElemSize,
						dvdt, gridData);
		}
	} else if (simData.irMethod == "analytic") {
		const FloatType minEdgeDivisor = simData.discretFactor;
		auto acField = std::make_unique<SimTransientAcousticField<FloatType, AnalyticRectangularSourceImpulseResponse<FloatType>>>();
		if (sourceIsArray) {
			acField->getArrayOfRectangularSourcesAcousticField(
						simData.samplingFreq, mainData.propagationSpeed, srcData.sourceWidth, srcData.sourceHeight,
						minEdgeDivisor,
						dvdt, srcData.elemPos, srcData.focusDelay, gridData);
		} else {
			acField->getRectangularSourceAcousticField(
						simData.samplingFreq, mainData.propagationSpeed, srcData.sourceWidth, srcData.sourceHeight,
						minEdgeDivisor,
						dvdt, gridData);
		}
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
SimRectangularSourceMethod<FloatType>::execTransientPropagation(bool sourceIsArray)
{
	const ParameterMap& taskPM = project_.taskParamMap();
	MainData mainData;
	SimulationData simData;
	loadData(taskPM, mainData, simData);
	SourceData srcData;
	loadSourceData(mainData, sourceIsArray, srcData);

	project_.createDirectory(mainData.outputDir, false);

	const ParamMapPtr propagPM = project_.getSubParamMap("propag_config_file");
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

	const FloatType nyquistLambda = Util::wavelength(mainData.propagationSpeed, mainData.nyquistRate);
	ImageGrid<FloatType>::get(*project_.getSubParamMap("grid_config_file"), nyquistLambda, gridData);

	if (simData.irMethod == "numeric") {
		const FloatType subElemSize = mainData.propagationSpeed / (mainData.nyquistRate * simData.discretFactor);
		auto propag = std::make_unique<SimTransientPropagation<FloatType, NumericRectangularSourceImpulseResponse<FloatType>>>();
		if (sourceIsArray) {
			propag->getArrayOfRectangularSourcesPropagation(
						simData.samplingFreq, mainData.propagationSpeed, srcData.sourceWidth, srcData.sourceHeight,
						subElemSize,
						dvdt, srcData.elemPos, srcData.focusDelay, propagIndexList, gridData);
		} else {
			propag->getRectangularSourcePropagation(
						simData.samplingFreq, mainData.propagationSpeed, srcData.sourceWidth, srcData.sourceHeight,
						subElemSize,
						dvdt, propagIndexList, gridData);
		}
	} else if (simData.irMethod == "analytic") {
		const FloatType minEdgeDivisor = simData.discretFactor;
		auto propag = std::make_unique<SimTransientPropagation<FloatType, AnalyticRectangularSourceImpulseResponse<FloatType>>>();
		if (sourceIsArray) {
			propag->getArrayOfRectangularSourcesPropagation(
						simData.samplingFreq, mainData.propagationSpeed, srcData.sourceWidth, srcData.sourceHeight,
						minEdgeDivisor,
						dvdt, srcData.elemPos, srcData.focusDelay, propagIndexList, gridData);
		} else {
			propag->getRectangularSourcePropagation(
						simData.samplingFreq, mainData.propagationSpeed, srcData.sourceWidth, srcData.sourceHeight,
						minEdgeDivisor,
						dvdt, propagIndexList, gridData);
		}
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
SimRectangularSourceMethod<FloatType>::execImpulseResponse(bool sourceIsArray)
{
	const ParameterMap& taskPM = project_.taskParamMap();
	MainData mainData;
	SimulationData simData;
	loadData(taskPM, mainData, simData);
	SourceData srcData;
	loadSourceData(mainData, sourceIsArray, srcData);

	project_.createDirectory(mainData.outputDir, false);

	const ParamMapPtr irPM = project_.getSubParamMap("ir_config_file");
	const auto pointX = irPM->value<FloatType>("point_x", sourceIsArray ? -10000.0 : 0.0, 10000.0);
	const auto pointY = irPM->value<FloatType>("point_y", sourceIsArray ? -10000.0 : 0.0, 10000.0);
	const auto pointZ = irPM->value<FloatType>("point_z", 0.0, 10000.0);

	const FloatType dt = 1.0 / simData.samplingFreq;
	std::vector<FloatType> tExc, dvdt, tDvdt;
	prepareExcitation(dt, simData, tExc, dvdt, tDvdt);

	std::size_t hOffset;
	std::vector<FloatType> h;
	if (simData.irMethod == "numeric") {
		const FloatType subElemSize = mainData.propagationSpeed / (mainData.nyquistRate * simData.discretFactor);
		if (sourceIsArray) {
			auto impResp = std::make_unique<ArrayOfRectangularSourcesImpulseResponse<FloatType, NumericRectangularSourceImpulseResponse<FloatType>>>(
						simData.samplingFreq, mainData.propagationSpeed, srcData.sourceWidth, srcData.sourceHeight,
						subElemSize,
						srcData.elemPos, srcData.focusDelay);
			impResp->getImpulseResponse(pointX, pointY, pointZ, hOffset, h);
		} else {
			auto impResp = std::make_unique<NumericRectangularSourceImpulseResponse<FloatType>>(
						simData.samplingFreq, mainData.propagationSpeed, srcData.sourceWidth, srcData.sourceHeight,
						subElemSize);
			impResp->getImpulseResponse(pointX, pointY, pointZ, hOffset, h);
		}
	} else if (simData.irMethod == "analytic") {
		const FloatType minEdgeDivisor = simData.discretFactor;
		if (sourceIsArray) {
			auto impResp = std::make_unique<ArrayOfRectangularSourcesImpulseResponse<FloatType, AnalyticRectangularSourceImpulseResponse<FloatType>>>(
						simData.samplingFreq, mainData.propagationSpeed, srcData.sourceWidth, srcData.sourceHeight,
						minEdgeDivisor,
						srcData.elemPos, srcData.focusDelay);
			impResp->getImpulseResponse(pointX, pointY, pointZ, hOffset, h);
		} else {
			auto impResp = std::make_unique<AnalyticRectangularSourceImpulseResponse<FloatType>>(
						simData.samplingFreq, mainData.propagationSpeed, srcData.sourceWidth, srcData.sourceHeight,
						minEdgeDivisor);
			impResp->getImpulseResponse(pointX, pointY, pointZ, hOffset, h);
		}
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
SimRectangularSourceMethod<FloatType>::execute()
{
	Timer tProc;

	switch (project_.method()) {
	case MethodEnum::sim_acoustic_field_rectangular_source_transient:
		execTransientAcousticField(false);
		break;
	case MethodEnum::sim_acoustic_field_array_of_rectangular_sources_transient:
		execTransientAcousticField(true);
		break;
	case MethodEnum::sim_impulse_response_rectangular_source:
		execImpulseResponse(false);
		break;
	case MethodEnum::sim_impulse_response_array_of_rectangular_sources:
		execImpulseResponse(true);
		break;
	case MethodEnum::sim_propagation_rectangular_source_transient:
		execTransientPropagation(false);
		break;
	case MethodEnum::sim_propagation_array_of_rectangular_sources_transient:
		execTransientPropagation(true);
		break;
	case MethodEnum::sim_radiation_pattern_rectangular_source_transient:
		execTransientRadiationPattern(false);
		break;
	case MethodEnum::sim_radiation_pattern_array_of_rectangular_sources_transient:
		execTransientRadiationPattern(true);
		break;
	default:
		THROW_EXCEPTION(InvalidParameterException, "Invalid method: " << static_cast<int>(project_.method()) << '.');
	}

	LOG_DEBUG << ">>> Processing time: " << tProc.getTime();
}

} // namespace Lab

#endif // SIMRECTANGULARSOURCEMETHOD_H
