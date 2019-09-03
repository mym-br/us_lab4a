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

#ifndef SIMRECTANGULARFLATSOURCEMETHOD_H
#define SIMRECTANGULARFLATSOURCEMETHOD_H

#include <cmath>
#include <iomanip>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "AnalyticRectangularFlatSourceImpulseResponse.h"
#include "ArrayOfRectangularFlatSourcesImpulseResponse.h"
#include "ArrayUtil.h"
#include "Exception.h"
#include "FFTWFilter2.h"
#include "ImageGrid.h"
#include "Log.h"
#include "Matrix.h"
#include "Method.h"
#include "NumericRectangularFlatSourceImpulseResponse.h"
#include "ParameterMap.h"
#include "Project.h"
#include "SimTransientAcousticField.h"
#include "SimTransientPropagation.h"
#include "SimTransientRadiationPattern.h"
#include "Timer.h"
#include "Util.h"
#include "Waveform.h"
#include "WavefrontObjFileWriter.h"
#include "XY.h"
#include "XYZ.h"
#include "XYZValue.h"
#include "XYZValueArray.h"

namespace Lab {

template<typename FloatType>
class SimRectangularFlatSourceMethod : public Method {
public:
	SimRectangularFlatSourceMethod(Project& project);
	virtual ~SimRectangularFlatSourceMethod();

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

	SimRectangularFlatSourceMethod(const SimRectangularFlatSourceMethod&) = delete;
	SimRectangularFlatSourceMethod& operator=(const SimRectangularFlatSourceMethod&) = delete;

	void loadData(ConstParameterMapPtr taskPM, MainData& data, SimulationData& simData);
	void loadSourceData(ConstParameterMapPtr taskPM, MainData& data, bool sourceIsArray, SourceData& srcData);
	void loadSimulationData(ConstParameterMapPtr mainPM, const MainData& data, SimulationData& simData);
	void prepareExcitation(FloatType dt, const SimulationData& simData, std::vector<FloatType>& tExc,
				std::vector<FloatType>& dvdt, std::vector<FloatType>& tDvdt);

	void execTransientRadiationPattern(bool sourceIsArray);
	void execTransientAcousticField(bool sourceIsArray);
	void execTransientPropagation(bool sourceIsArray);
	// Returns p/(c*density).
	void execImpulseResponse(bool sourceIsArray);

	Project& project_;
};

template<typename FloatType>
SimRectangularFlatSourceMethod<FloatType>::SimRectangularFlatSourceMethod(Project& project)
		: project_(project)
{
}

template<typename FloatType>
SimRectangularFlatSourceMethod<FloatType>::~SimRectangularFlatSourceMethod()
{
}

template<typename FloatType>
void
SimRectangularFlatSourceMethod<FloatType>::loadData(ConstParameterMapPtr taskPM, MainData& data, SimulationData& simData)
{
	ConstParameterMapPtr mainPM = project_.loadChildParameterMap(taskPM, "main_config_file");
	//data.density        = mainPM->value<FloatType>("density", 0.0, 100000.0);
	data.propagationSpeed = mainPM->value<FloatType>("propagation_speed", 0.0, 100000.0);
	data.centerFreq       = mainPM->value<FloatType>("center_frequency", 0.0, 100.0e6);
	data.maxFreq          = mainPM->value<FloatType>("max_frequency", 0.0, 200.0e6);
	data.nyquistRate = 2.0 * data.maxFreq;
	data.outputDir        = mainPM->value<std::string>("output_dir");

	loadSimulationData(mainPM, data, simData);
}

template<typename FloatType>
void
SimRectangularFlatSourceMethod<FloatType>::loadSourceData(ConstParameterMapPtr taskPM, MainData& data, bool sourceIsArray, SourceData& srcData)
{
	WavefrontObjFileWriter<FloatType> fw((project_.directory() + "/source_geometry.obj").c_str());

	if (sourceIsArray) {
		ConstParameterMapPtr arrayPM = project_.loadChildParameterMap(taskPM, "array_config_file");
		ArrayUtil::calculateTxElementPositions(*arrayPM, srcData.elemPos);
		srcData.sourceWidth  = arrayPM->value<FloatType>("element_width", 0.0, 10.0);
		srcData.sourceHeight = arrayPM->value<FloatType>("element_height", 0.0, 10.0);

		srcData.useFocus = arrayPM->value<bool>("use_tx_focus");
		if (srcData.useFocus) {
			srcData.focusX = arrayPM->value<FloatType>("tx_focus_x", -10000.0, 10000.0);
			srcData.focusY = arrayPM->value<FloatType>("tx_focus_y", -10000.0, 10000.0);
			srcData.focusZ = arrayPM->value<FloatType>("tx_focus_z", -10000.0, 10000.0);
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
		ConstParameterMapPtr singlePM = project_.loadChildParameterMap(taskPM, "single_source_config_file");
		srcData.sourceWidth  = singlePM->value<FloatType>("source_width", 0.0, 10.0);
		srcData.sourceHeight = singlePM->value<FloatType>("source_height", 0.0, 10.0);

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
SimRectangularFlatSourceMethod<FloatType>::loadSimulationData(ConstParameterMapPtr mainPM, const MainData& data, SimulationData& simData)
{
	ConstParameterMapPtr  simPM = project_.loadChildParameterMap(mainPM, "simulation_config_file");
	simData.samplingFreq   = simPM->value<FloatType>("sampling_frequency_factor", 0.0, 10000.0) * data.nyquistRate;
	simData.irMethod       = simPM->value<std::string>("impulse_response_method");
	simData.excitationType = simPM->value<std::string>("excitation_type");
	simData.excNumPeriods  = simPM->value<FloatType>("excitation_num_periods", 0.0, 100.0);

	if (simData.excitationType == "1") {
		Waveform::getType1(data.centerFreq, simData.samplingFreq, simData.excNumPeriods, simData.exc);
	} else if (simData.excitationType == "2a") {
		Waveform::getType2a(data.centerFreq, simData.samplingFreq, simData.excNumPeriods, simData.exc);
	} else if (simData.excitationType == "2b") {
		Waveform::getType2b(data.centerFreq, simData.samplingFreq, simData.excNumPeriods, simData.exc);
	} else {
		THROW_EXCEPTION(InvalidParameterException, "Invalid excitation type: " << simData.excitationType << '.');
	}

	if (simData.irMethod == "numeric") {
		simData.discretFactor = simPM->value<FloatType>("sub_elem_size_factor", 0.0, 1.0e3);
	} else if (simData.irMethod == "analytic") {
		simData.discretFactor = simPM->value<FloatType>("min_edge_divisor", 0.0, 1.0e6);
	} else {
		THROW_EXCEPTION(InvalidParameterException, "Invalid impulse response method: " << simData.irMethod << '.');
	}
}

template<typename FloatType>
void
SimRectangularFlatSourceMethod<FloatType>::prepareExcitation(FloatType dt,
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
SimRectangularFlatSourceMethod<FloatType>::execTransientRadiationPattern(bool sourceIsArray)
{
	ConstParameterMapPtr taskPM = project_.taskParameterMap();
	MainData mainData;
	SimulationData simData;
	loadData(taskPM, mainData, simData);
	SourceData srcData;
	loadSourceData(taskPM, mainData, sourceIsArray, srcData);

	project_.createDirectory(mainData.outputDir, false);

	const FloatType distance   = taskPM->value<FloatType>("distance", 0.0, 100.0);
	const FloatType thetaYStep = taskPM->value<FloatType>("theta_y_step", 0.0, 10.0);
	const FloatType thetaYMin  = sourceIsArray ? taskPM->value<FloatType>("theta_y_min" , -90.0, 90.0) : FloatType(0);
	const FloatType thetaYMax  = taskPM->value<FloatType>("theta_y_max" , thetaYMin + 0.1, 90.0);
	const FloatType thetaXStep = taskPM->value<FloatType>("theta_x_step", 0.0, 10.0);
	const FloatType thetaXMin  = sourceIsArray ? taskPM->value<FloatType>("theta_x_min" , -90.0, 90.0) : FloatType(0);
	const FloatType thetaXMax  = taskPM->value<FloatType>("theta_x_max" , thetaXMin + 0.1, 90.0);

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
		auto radPat = std::make_unique<SimTransientRadiationPattern<FloatType, NumericRectangularFlatSourceImpulseResponse<FloatType>>>();
		if (sourceIsArray) {
			radPat->getArrayOfRectangularFlatSourcesRadiationPattern(
						simData.samplingFreq, mainData.propagationSpeed, srcData.sourceWidth, srcData.sourceHeight,
						subElemSize,
						dvdt, srcData.elemPos, srcData.focusDelay, inputData, gridData);
		} else {
			radPat->getRectangularFlatSourceRadiationPattern(
						simData.samplingFreq, mainData.propagationSpeed, srcData.sourceWidth, srcData.sourceHeight,
						subElemSize,
						dvdt, inputData, gridData);
		}
	} else if (simData.irMethod == "analytic") {
		const FloatType minEdgeDivisor = simData.discretFactor;
		auto radPat = std::make_unique<SimTransientRadiationPattern<FloatType, AnalyticRectangularFlatSourceImpulseResponse<FloatType>>>();
		if (sourceIsArray) {
			radPat->getArrayOfRectangularFlatSourcesRadiationPattern(
						simData.samplingFreq, mainData.propagationSpeed, srcData.sourceWidth, srcData.sourceHeight,
						minEdgeDivisor,
						dvdt, srcData.elemPos, srcData.focusDelay, inputData, gridData);
		} else {
			radPat->getRectangularFlatSourceRadiationPattern(
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
					true, Figure::VISUALIZATION_RECTIFIED_LOG, Figure::COLORMAP_VIRIDIS);

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
	project_.saveHDF5(thetaYList , mainData.outputDir + "/theta_y"        , "value");
	project_.saveHDF5(patternTY  , mainData.outputDir + "/pattern_theta_y", "value");
	project_.saveHDF5(thetaXList , mainData.outputDir + "/theta_x"        , "value");
	project_.saveHDF5(patternTX  , mainData.outputDir + "/pattern_theta_x", "value");
}

template<typename FloatType>
void
SimRectangularFlatSourceMethod<FloatType>::execTransientAcousticField(bool sourceIsArray)
{
	ConstParameterMapPtr taskPM = project_.taskParameterMap();
	MainData mainData;
	SimulationData simData;
	loadData(taskPM, mainData, simData);
	SourceData srcData;
	loadSourceData(taskPM, mainData, sourceIsArray, srcData);

	project_.createDirectory(mainData.outputDir, false);

	const FloatType dt = 1.0 / simData.samplingFreq;
	std::vector<FloatType> tExc, dvdt, tDvdt;
	prepareExcitation(dt, simData, tExc, dvdt, tDvdt);

	Matrix<XYZValue<FloatType>> gridData;

	const FloatType nyquistLambda = mainData.propagationSpeed / mainData.nyquistRate;
	ImageGrid<FloatType>::get(project_.loadChildParameterMap(taskPM, "grid_config_file"), nyquistLambda, gridData);

	if (simData.irMethod == "numeric") {
		const FloatType subElemSize = mainData.propagationSpeed / (mainData.nyquistRate * simData.discretFactor);
		auto acField = std::make_unique<SimTransientAcousticField<FloatType, NumericRectangularFlatSourceImpulseResponse<FloatType>>>();
		if (sourceIsArray) {
			acField->getArrayOfRectangularFlatSourcesAcousticField(
						simData.samplingFreq, mainData.propagationSpeed, srcData.sourceWidth, srcData.sourceHeight,
						subElemSize,
						dvdt, srcData.elemPos, srcData.focusDelay, gridData);
		} else {
			acField->getRectangularFlatSourceAcousticField(
						simData.samplingFreq, mainData.propagationSpeed, srcData.sourceWidth, srcData.sourceHeight,
						subElemSize,
						dvdt, gridData);
		}
	} else if (simData.irMethod == "analytic") {
		const FloatType minEdgeDivisor = simData.discretFactor;
		auto acField = std::make_unique<SimTransientAcousticField<FloatType, AnalyticRectangularFlatSourceImpulseResponse<FloatType>>>();
		if (sourceIsArray) {
			acField->getArrayOfRectangularFlatSourcesAcousticField(
						simData.samplingFreq, mainData.propagationSpeed, srcData.sourceWidth, srcData.sourceHeight,
						minEdgeDivisor,
						dvdt, srcData.elemPos, srcData.focusDelay, gridData);
		} else {
			acField->getRectangularFlatSourceAcousticField(
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
					true, Figure::VISUALIZATION_RECTIFIED_LINEAR, Figure::COLORMAP_VIRIDIS);

	project_.saveHDF5(simData.exc, mainData.outputDir + "/v"        , "value");
	project_.saveHDF5(tExc       , mainData.outputDir + "/v_time"   , "value");
	project_.saveHDF5(dvdt       , mainData.outputDir + "/dvdt"     , "value");
	project_.saveHDF5(tDvdt      , mainData.outputDir + "/dvdt_time", "value");
	project_.saveImageToHDF5(gridData, mainData.outputDir);
}

template<typename FloatType>
void
SimRectangularFlatSourceMethod<FloatType>::execTransientPropagation(bool sourceIsArray)
{
	ConstParameterMapPtr taskPM = project_.taskParameterMap();
	MainData mainData;
	SimulationData simData;
	loadData(taskPM, mainData, simData);
	SourceData srcData;
	loadSourceData(taskPM, mainData, sourceIsArray, srcData);

	project_.createDirectory(mainData.outputDir, false);

	const FloatType propagDistanceStep = taskPM->value<FloatType>("propagation_distance_step", 1.0e-6, 100.0);
	const FloatType propagMinDistance  = taskPM->value<FloatType>("propagation_min_distance", 0.0, 100.0);
	const FloatType propagMaxDistance  = taskPM->value<FloatType>("propagation_max_distance", propagMinDistance + propagDistanceStep, 100.0);
	const unsigned int propagRepet     = taskPM->value<unsigned int>("propagation_repetitions", 1, 100);
	const unsigned int propagPause     = taskPM->value<unsigned int>("propagation_pause", 0, 10000);

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
		const FloatType subElemSize = mainData.propagationSpeed / (mainData.nyquistRate * simData.discretFactor);
		auto propag = std::make_unique<SimTransientPropagation<FloatType, NumericRectangularFlatSourceImpulseResponse<FloatType>>>();
		if (sourceIsArray) {
			propag->getArrayOfRectangularFlatSourcesPropagation(
						simData.samplingFreq, mainData.propagationSpeed, srcData.sourceWidth, srcData.sourceHeight,
						subElemSize,
						dvdt, srcData.elemPos, srcData.focusDelay, propagIndexList, gridData);
		} else {
			propag->getRectangularFlatSourcePropagation(
						simData.samplingFreq, mainData.propagationSpeed, srcData.sourceWidth, srcData.sourceHeight,
						subElemSize,
						dvdt, propagIndexList, gridData);
		}
	} else if (simData.irMethod == "analytic") {
		const FloatType minEdgeDivisor = simData.discretFactor;
		auto propag = std::make_unique<SimTransientPropagation<FloatType, AnalyticRectangularFlatSourceImpulseResponse<FloatType>>>();
		if (sourceIsArray) {
			propag->getArrayOfRectangularFlatSourcesPropagation(
						simData.samplingFreq, mainData.propagationSpeed, srcData.sourceWidth, srcData.sourceHeight,
						minEdgeDivisor,
						dvdt, srcData.elemPos, srcData.focusDelay, propagIndexList, gridData);
		} else {
			propag->getRectangularFlatSourcePropagation(
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
						true, Figure::VISUALIZATION_RAW_LINEAR, Figure::COLORMAP_GRAY);

			Util::sleepMs(propagPause);
		}
	}

	project_.saveHDF5(simData.exc, mainData.outputDir + "/v"        , "value");
	project_.saveHDF5(tExc       , mainData.outputDir + "/v_time"   , "value");
	project_.saveHDF5(dvdt       , mainData.outputDir + "/dvdt"     , "value");
	project_.saveHDF5(tDvdt      , mainData.outputDir + "/dvdt_time", "value");

	// Save images.
	LOG_DEBUG << "Saving the X coordinates...";
	project_.saveHDF5(gridData, mainData.outputDir + "/image_x", "x", Util::CopyXOp());
	LOG_DEBUG << "Saving the Y coordinates...";
	project_.saveHDF5(gridData, mainData.outputDir + "/image_y", "y", Util::CopyYOp());
	LOG_DEBUG << "Saving the Z coordinates...";
	project_.saveHDF5(gridData, mainData.outputDir + "/image_z", "z", Util::CopyZOp());
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
SimRectangularFlatSourceMethod<FloatType>::execImpulseResponse(bool sourceIsArray)
{
	ConstParameterMapPtr taskPM = project_.taskParameterMap();
	MainData mainData;
	SimulationData simData;
	loadData(taskPM, mainData, simData);
	SourceData srcData;
	loadSourceData(taskPM, mainData, sourceIsArray, srcData);

	project_.createDirectory(mainData.outputDir, false);

	const FloatType pointX = taskPM->value<FloatType>("point_x", sourceIsArray ? -10000.0 : 0.0, 10000.0);
	const FloatType pointY = taskPM->value<FloatType>("point_y", sourceIsArray ? -10000.0 : 0.0, 10000.0);
	const FloatType pointZ = taskPM->value<FloatType>("point_z", 0.0, 10000.0);

	const FloatType dt = 1.0 / simData.samplingFreq;
	std::vector<FloatType> tExc, dvdt, tDvdt;
	prepareExcitation(dt, simData, tExc, dvdt, tDvdt);

	std::size_t hOffset;
	std::vector<FloatType> h;
	if (simData.irMethod == "numeric") {
		const FloatType subElemSize = mainData.propagationSpeed / (mainData.nyquistRate * simData.discretFactor);
		if (sourceIsArray) {
			auto impResp = std::make_unique<ArrayOfRectangularFlatSourcesImpulseResponse<FloatType, NumericRectangularFlatSourceImpulseResponse<FloatType>>>(
						simData.samplingFreq, mainData.propagationSpeed, srcData.sourceWidth, srcData.sourceHeight,
						subElemSize,
						srcData.elemPos, srcData.focusDelay);
			impResp->getImpulseResponse(pointX, pointY, pointZ, hOffset, h);
		} else {
			auto impResp = std::make_unique<NumericRectangularFlatSourceImpulseResponse<FloatType>>(
						simData.samplingFreq, mainData.propagationSpeed, srcData.sourceWidth, srcData.sourceHeight,
						subElemSize);
			impResp->getImpulseResponse(pointX, pointY, pointZ, hOffset, h);
		}
	} else if (simData.irMethod == "analytic") {
		const FloatType minEdgeDivisor = simData.discretFactor;
		if (sourceIsArray) {
			auto impResp = std::make_unique<ArrayOfRectangularFlatSourcesImpulseResponse<FloatType, AnalyticRectangularFlatSourceImpulseResponse<FloatType>>>(
						simData.samplingFreq, mainData.propagationSpeed, srcData.sourceWidth, srcData.sourceHeight,
						minEdgeDivisor,
						srcData.elemPos, srcData.focusDelay);
			impResp->getImpulseResponse(pointX, pointY, pointZ, hOffset, h);
		} else {
			auto impResp = std::make_unique<AnalyticRectangularFlatSourceImpulseResponse<FloatType>>(
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
SimRectangularFlatSourceMethod<FloatType>::execute()
{
	Timer tProc;

	switch (project_.method()) {
	case MethodEnum::sim_acoustic_field_rectangular_flat_source_transient:
		execTransientAcousticField(false);
		break;
	case MethodEnum::sim_acoustic_field_array_of_rectangular_flat_sources_transient:
		execTransientAcousticField(true);
		break;
	case MethodEnum::sim_impulse_response_rectangular_flat_source:
		execImpulseResponse(false);
		break;
	case MethodEnum::sim_impulse_response_array_of_rectangular_flat_sources:
		execImpulseResponse(true);
		break;
	case MethodEnum::sim_propagation_rectangular_flat_source_transient:
		execTransientPropagation(false);
		break;
	case MethodEnum::sim_propagation_array_of_rectangular_flat_sources_transient:
		execTransientPropagation(true);
		break;
	case MethodEnum::sim_radiation_pattern_rectangular_flat_source_transient:
		execTransientRadiationPattern(false);
		break;
	case MethodEnum::sim_radiation_pattern_array_of_rectangular_flat_sources_transient:
		execTransientRadiationPattern(true);
		break;
	default:
		THROW_EXCEPTION(InvalidParameterException, "Invalid method: " << static_cast<int>(project_.method()) << '.');
	}

	LOG_DEBUG << ">>> Processing time: " << tProc.getTime();
}

} // namespace Lab

#endif // SIMRECTANGULARFLATSOURCEMETHOD_H
