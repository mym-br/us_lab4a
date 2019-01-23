/***************************************************************************
 *  Copyright 2018 Marcelo Y. Matuda                                       *
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
#include <memory>
#include <string>

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
	SimRectangularFlatSourceMethod(const SimRectangularFlatSourceMethod&) = delete;
	SimRectangularFlatSourceMethod& operator=(const SimRectangularFlatSourceMethod&) = delete;

	void execTransientRadiationPattern(bool sourceIsArray);
	void execTransientAcousticField();
	void execTransientArrayAcousticField();
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
SimRectangularFlatSourceMethod<FloatType>::execTransientRadiationPattern(bool sourceIsArray)
{
	ConstParameterMapPtr taskPM = project_.taskParameterMap();

	const std::string outputDir = taskPM->value<std::string>("output_dir");
	project_.createDirectory(outputDir, false);

	const std::string irMethod       = taskPM->value<std::string>("impulse_response_method");
	const FloatType distance         = taskPM->value<FloatType>("distance", 0.0, 100.0);
	const FloatType thetaYStep       = taskPM->value<FloatType>("theta_y_step", 0.0, 10.0);
	const FloatType thetaYMin        = sourceIsArray ? taskPM->value<FloatType>("theta_y_min" , -90.0, 90.0) : FloatType(0);
	const FloatType thetaYMax        = taskPM->value<FloatType>("theta_y_max" , thetaYMin + 0.1, 90.0);
	const FloatType thetaXStep       = taskPM->value<FloatType>("theta_x_step", 0.0, 10.0);
	const FloatType thetaXMin        = sourceIsArray ? taskPM->value<FloatType>("theta_x_min" , -90.0, 90.0) : FloatType(0);
	const FloatType thetaXMax        = taskPM->value<FloatType>("theta_x_max" , thetaXMin + 0.1, 90.0);
	const FloatType sourceWidth      = taskPM->value<FloatType>("source_width", 0.0, 10.0);
	const FloatType sourceHeight     = taskPM->value<FloatType>("source_height", 0.0, 10.0);
	const FloatType propagationSpeed = taskPM->value<FloatType>("propagation_speed", 0.0, 100000.0);
	const FloatType centerFreq       = taskPM->value<FloatType>("center_frequency", 0.0, 100.0e6);
	const FloatType maxFreq          = taskPM->value<FloatType>("max_frequency", 0.0, 200.0e6);
	const FloatType nyquistRate = 2.0 * maxFreq;
	const FloatType samplingFreq     = taskPM->value<FloatType>("sampling_frequency_factor", 0.0, 10000.0) * nyquistRate;
	const std::string excitationType = taskPM->value<std::string>("excitation_type");
	const FloatType excNumPeriods    = taskPM->value<FloatType>("excitation_num_periods", 0.0, 100.0);

	std::vector<XY<FloatType>> elemPos;
	std::vector<FloatType> focusDelay;
	if (sourceIsArray) {
		ConstParameterMapPtr arrayPM = project_.loadChildParameterMap(taskPM, "array_config_file");
		ArrayUtil::calculateTxElementPositions(*arrayPM, elemPos);
		ArrayUtil::calculateTx3DFocusDelay(*taskPM, propagationSpeed, elemPos, focusDelay);
	}

	std::vector<FloatType> exc;
	if (excitationType == "1") {
		Waveform::getType1(centerFreq, samplingFreq, exc, excNumPeriods);
	} else if (excitationType == "2a") {
		Waveform::getType2a(centerFreq, samplingFreq, exc, excNumPeriods);
	} else if (excitationType == "2b") {
		Waveform::getType2b(centerFreq, samplingFreq, exc, excNumPeriods);
	} else {
		THROW_EXCEPTION(InvalidParameterException, "Invalid excitation type: " << excitationType << '.');
	}

	const FloatType dt = 1.0 / samplingFreq;
	std::vector<FloatType> tExc(exc.size());
	for (unsigned int i = 0; i < tExc.size(); ++i) {
		tExc[i] = dt * i;
	}
	project_.showFigure2D(1, "Excitation", tExc, exc);

	std::vector<FloatType> dvdt;
	Util::centralDiff(exc, dt, dvdt);

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

	if (irMethod == "numeric") {
		const FloatType subElemSize = propagationSpeed /
				(nyquistRate * taskPM->value<FloatType>("sub_elem_size_factor", 0.0, 1.0e3));
		auto radPat = std::make_unique<SimTransientRadiationPattern<FloatType, NumericRectangularFlatSourceImpulseResponse<FloatType>>>();
		if (sourceIsArray) {
			radPat->getArrayOfRectangularFlatSourcesRadiationPattern(
						samplingFreq, propagationSpeed, sourceWidth, sourceHeight,
						subElemSize,
						dvdt, elemPos, focusDelay, inputData, gridData);
		} else {
			radPat->getRectangularFlatSourceRadiationPattern(
						samplingFreq, propagationSpeed, sourceWidth, sourceHeight,
						subElemSize,
						dvdt, inputData, gridData);
		}
	} else if (irMethod == "analytic") {
		const FloatType minEdgeDivisor = taskPM->value<FloatType>("min_edge_divisor", 0.0, 1.0e6);
		auto radPat = std::make_unique<SimTransientRadiationPattern<FloatType, AnalyticRectangularFlatSourceImpulseResponse<FloatType>>>();
		if (sourceIsArray) {
			radPat->getArrayOfRectangularFlatSourcesRadiationPattern(
						samplingFreq, propagationSpeed, sourceWidth, sourceHeight,
						minEdgeDivisor,
						dvdt, elemPos, focusDelay, inputData, gridData);
		} else {
			radPat->getRectangularFlatSourceRadiationPattern(
						samplingFreq, propagationSpeed, sourceWidth, sourceHeight,
						minEdgeDivisor,
						dvdt, inputData, gridData);
		}
	} else {
		THROW_EXCEPTION(InvalidParameterException, "Invalid impulse response method: " << irMethod << '.');
	}

	const FloatType maxAbsValue = Util::maxAbsoluteValueField<XYZValue<FloatType>, FloatType>(gridData);
	const FloatType k = 1.0 / maxAbsValue;
	for (auto it = gridData.begin(); it != gridData.end(); ++it) {
		it->value *= k;
	}

	std::vector<XYZ<float>> pointList = {{0.0, 0.0, 0.0}};

	Project::GridDataType projGridData;
	Util::copyXYZValue(gridData, projGridData);
	project_.showFigure3D(1, "Pattern", &projGridData, &pointList,
					true, Figure::VISUALIZATION_RECTIFIED_LOG, Figure::COLORMAP_VIRIDIS);

	std::size_t sectionTYIndex = 0;
	std::size_t sectionTXIndex = 0;
	if (sourceIsArray) {
		FloatType sectionTY = 0.0;
		FloatType sectionTX = 0.0;
		const bool useFocus = taskPM->value<bool>("use_tx_focus");
		if (useFocus) {
			const FloatType focusX = taskPM->value<FloatType>("tx_focus_x", -10000.0, 10000.0);
			const FloatType focusY = taskPM->value<FloatType>("tx_focus_y", -10000.0, 10000.0);
			const FloatType focusZ = taskPM->value<FloatType>("tx_focus_z", -10000.0, 10000.0);
			sectionTY = Util::radianToDegree(std::atan2(focusX, focusZ));
			sectionTX = Util::radianToDegree(std::atan2(focusY, focusZ));
		}
		sectionTYIndex = std::rint(((sectionTY - thetaYMin) / (thetaYMax - thetaYMin)) * (thetaYList.size() - 1));
		sectionTXIndex = std::rint(((sectionTX - thetaXMin) / (thetaXMax - thetaXMin)) * (thetaXList.size() - 1));
		LOG_DEBUG << "Section theta-x: " << thetaXList[sectionTXIndex] << " theta-y: " << thetaYList[sectionTYIndex];
	}

	std::vector<FloatType> patternTY(gridData.n2());
	auto tyInterval = gridData.dim2Interval(sectionTXIndex);
	Util::copyUsingOperator(tyInterval.first, tyInterval.second, patternTY.begin(), Util::CopyValueOp{});
	Util::linearToDecibels(patternTY, -100.0);
	project_.showFigure2D(1, "Pattern theta-y", thetaYList, patternTY);

	std::vector<FloatType> patternTX(gridData.n1());
	auto txInterval = gridData.dim1Interval(sectionTYIndex);
	Util::copyUsingOperator(txInterval.first, txInterval.second, patternTX.begin(), Util::CopyValueOp{});
	Util::linearToDecibels(patternTX, -100.0);
	project_.showFigure2D(2, "Pattern theta-x", thetaXList, patternTX);

	project_.saveHDF5(exc , outputDir + "/excitation"     , "value");
	project_.saveHDF5(tExc, outputDir + "/excitation_time", "value");
	project_.saveImageToHDF5(gridData, outputDir);
	project_.saveHDF5(thetaYList, outputDir + "/theta_y"        , "value");
	project_.saveHDF5(patternTY , outputDir + "/pattern_theta_y", "value");
	project_.saveHDF5(thetaXList, outputDir + "/theta_x"        , "value");
	project_.saveHDF5(patternTX , outputDir + "/pattern_theta_x", "value");
}

template<typename FloatType>
void
SimRectangularFlatSourceMethod<FloatType>::execTransientAcousticField()
{
	ConstParameterMapPtr taskPM = project_.taskParameterMap();

	const std::string outputDir = taskPM->value<std::string>("output_dir");
	project_.createDirectory(outputDir, false);

	const std::string irMethod       = taskPM->value<std::string>("impulse_response_method");
	const FloatType sourceWidth      = taskPM->value<FloatType>("source_width", 0.0, 10.0);
	const FloatType sourceHeight     = taskPM->value<FloatType>("source_height", 0.0, 10.0);
	const FloatType propagationSpeed = taskPM->value<FloatType>("propagation_speed", 0.0, 100000.0);
	const FloatType centerFreq       = taskPM->value<FloatType>("center_frequency", 0.0, 100.0e6);
	const FloatType maxFreq          = taskPM->value<FloatType>("max_frequency", 0.0, 200.0e6);
	const FloatType nyquistRate = 2.0 * maxFreq;
	const FloatType samplingFreq     = taskPM->value<FloatType>("sampling_frequency_factor", 0.0, 10000.0) * nyquistRate;
	const std::string excitationType = taskPM->value<std::string>("excitation_type");
	const FloatType excNumPeriods    = taskPM->value<FloatType>("excitation_num_periods", 0.0, 100.0);

	std::vector<FloatType> exc;
	if (excitationType == "1") {
		Waveform::getType1(centerFreq, samplingFreq, exc, excNumPeriods);
	} else if (excitationType == "2a") {
		Waveform::getType2a(centerFreq, samplingFreq, exc, excNumPeriods);
	} else if (excitationType == "2b") {
		Waveform::getType2b(centerFreq, samplingFreq, exc, excNumPeriods);
	} else {
		THROW_EXCEPTION(InvalidParameterException, "Invalid excitation type: " << excitationType << '.');
	}

	const FloatType dt = 1.0 / samplingFreq;
	std::vector<FloatType> tExc(exc.size());
	for (unsigned int i = 0; i < tExc.size(); ++i) {
		tExc[i] = dt * i;
	}
	project_.showFigure2D(1, "Excitation", tExc, exc);

	std::vector<FloatType> dvdt;
	Util::centralDiff(exc, dt, dvdt);

	Matrix<XYZValue<FloatType>> gridData;

	const FloatType nyquistLambda = propagationSpeed / nyquistRate;
	ImageGrid<FloatType>::get(project_.loadChildParameterMap(taskPM, "grid_config_file"), nyquistLambda, gridData);

	if (irMethod == "numeric") {
		const FloatType subElemSize = propagationSpeed /
				(nyquistRate * taskPM->value<FloatType>("sub_elem_size_factor", 0.0, 1.0e3));
		auto acField = std::make_unique<SimTransientAcousticField<FloatType, NumericRectangularFlatSourceImpulseResponse<FloatType>>>();
		acField->getRectangularFlatSourceAcousticField(
					samplingFreq, propagationSpeed, sourceWidth, sourceHeight,
					subElemSize,
					dvdt, gridData);
	} else if (irMethod == "analytic") {
		const FloatType minEdgeDivisor = taskPM->value<FloatType>("min_edge_divisor", 0.0, 1.0e6);
		auto acField = std::make_unique<SimTransientAcousticField<FloatType, AnalyticRectangularFlatSourceImpulseResponse<FloatType>>>();
		acField->getRectangularFlatSourceAcousticField(
					samplingFreq, propagationSpeed, sourceWidth, sourceHeight,
					minEdgeDivisor,
					dvdt, gridData);
	} else {
		THROW_EXCEPTION(InvalidParameterException, "Invalid impulse response method: " << irMethod << '.');
	}

	const FloatType maxAbsValue = Util::maxAbsoluteValueField<XYZValue<FloatType>, FloatType>(gridData);
	const FloatType k = 1.0 / maxAbsValue;
	for (auto it = gridData.begin(); it != gridData.end(); ++it) {
		it->value *= k;
	}

	std::vector<XYZ<float>> pointList = {{0.0, 0.0, 0.0}};

	Project::GridDataType projGridData;
	Util::copyXYZValue(gridData, projGridData);
	project_.showFigure3D(1, "Acoustic field", &projGridData, &pointList,
					true, Figure::VISUALIZATION_RECTIFIED_LINEAR, Figure::COLORMAP_VIRIDIS);

	project_.saveHDF5(exc , outputDir + "/excitation"     , "value");
	project_.saveHDF5(tExc, outputDir + "/excitation_time", "value");
	project_.saveImageToHDF5(gridData, outputDir);
}

template<typename FloatType>
void
SimRectangularFlatSourceMethod<FloatType>::execTransientArrayAcousticField()
{
	ConstParameterMapPtr taskPM = project_.taskParameterMap();

	const std::string outputDir = taskPM->value<std::string>("output_dir");
	project_.createDirectory(outputDir, false);

	const std::string irMethod       = taskPM->value<std::string>("impulse_response_method");
	const FloatType sourceWidth      = taskPM->value<FloatType>("source_width", 0.0, 10.0);
	const FloatType sourceHeight     = taskPM->value<FloatType>("source_height", 0.0, 10.0);
	const FloatType propagationSpeed = taskPM->value<FloatType>("propagation_speed", 0.0, 100000.0);
	const FloatType centerFreq       = taskPM->value<FloatType>("center_frequency", 0.0, 100.0e6);
	const FloatType maxFreq          = taskPM->value<FloatType>("max_frequency", 0.0, 200.0e6);
	const FloatType nyquistRate = 2.0 * maxFreq;
	const FloatType samplingFreq     = taskPM->value<FloatType>("sampling_frequency_factor", 0.0, 10000.0) * nyquistRate;
	const std::string excitationType = taskPM->value<std::string>("excitation_type");
	const FloatType excNumPeriods    = taskPM->value<FloatType>("excitation_num_periods", 0.0, 100.0);

	std::vector<XY<FloatType>> elemPos;
	std::vector<FloatType> focusDelay;
	ConstParameterMapPtr arrayPM = project_.loadChildParameterMap(taskPM, "array_config_file");
	ArrayUtil::calculateTxElementPositions(*arrayPM, elemPos);
	ArrayUtil::calculateTx3DFocusDelay(*taskPM, propagationSpeed, elemPos, focusDelay);

	std::vector<FloatType> exc;
	if (excitationType == "1") {
		Waveform::getType1(centerFreq, samplingFreq, exc, excNumPeriods);
	} else if (excitationType == "2a") {
		Waveform::getType2a(centerFreq, samplingFreq, exc, excNumPeriods);
	} else if (excitationType == "2b") {
		Waveform::getType2b(centerFreq, samplingFreq, exc, excNumPeriods);
	} else {
		THROW_EXCEPTION(InvalidParameterException, "Invalid excitation type: " << excitationType << '.');
	}

	const FloatType dt = 1.0 / samplingFreq;
	std::vector<FloatType> tExc(exc.size());
	for (unsigned int i = 0; i < tExc.size(); ++i) {
		tExc[i] = dt * i;
	}
	project_.showFigure2D(1, "Excitation", tExc, exc);

	std::vector<FloatType> dvdt;
	Util::centralDiff(exc, dt, dvdt);

	Matrix<XYZValue<FloatType>> gridData;

	const FloatType nyquistLambda = propagationSpeed / nyquistRate;
	ImageGrid<FloatType>::get(project_.loadChildParameterMap(taskPM, "grid_config_file"), nyquistLambda, gridData);

	if (irMethod == "numeric") {
		const FloatType subElemSize = propagationSpeed /
				(nyquistRate * taskPM->value<FloatType>("sub_elem_size_factor", 0.0, 1.0e3));
		auto acField = std::make_unique<SimTransientAcousticField<FloatType, NumericRectangularFlatSourceImpulseResponse<FloatType>>>();
		acField->getArrayOfRectangularFlatSourcesAcousticField(
					samplingFreq, propagationSpeed, sourceWidth, sourceHeight,
					subElemSize,
					dvdt, elemPos, focusDelay, gridData);
	} else if (irMethod == "analytic") {
		const FloatType minEdgeDivisor = taskPM->value<FloatType>("min_edge_divisor", 0.0, 1.0e6);
		auto acField = std::make_unique<SimTransientAcousticField<FloatType, AnalyticRectangularFlatSourceImpulseResponse<FloatType>>>();
		acField->getArrayOfRectangularFlatSourcesAcousticField(
					samplingFreq, propagationSpeed, sourceWidth, sourceHeight,
					minEdgeDivisor,
					dvdt, elemPos, focusDelay, gridData);
	} else {
		THROW_EXCEPTION(InvalidParameterException, "Invalid impulse response method: " << irMethod << '.');
	}

	const FloatType maxAbsValue = Util::maxAbsoluteValueField<XYZValue<FloatType>, FloatType>(gridData);
	const FloatType k = 1.0 / maxAbsValue;
	for (auto it = gridData.begin(); it != gridData.end(); ++it) {
		it->value *= k;
	}

	std::vector<XYZ<float>> pointList = {{0.0, 0.0, 0.0}};

	Project::GridDataType projGridData;
	Util::copyXYZValue(gridData, projGridData);
	project_.showFigure3D(1, "Acoustic field", &projGridData, &pointList,
					true, Figure::VISUALIZATION_RECTIFIED_LINEAR, Figure::COLORMAP_VIRIDIS);

	project_.saveHDF5(exc , outputDir + "/excitation"     , "value");
	project_.saveHDF5(tExc, outputDir + "/excitation_time", "value");
	project_.saveImageToHDF5(gridData, outputDir);
}

template<typename FloatType>
void
SimRectangularFlatSourceMethod<FloatType>::execTransientPropagation(bool sourceIsArray)
{
	ConstParameterMapPtr taskPM = project_.taskParameterMap();

	const std::string outputDir = taskPM->value<std::string>("output_dir");
	project_.createDirectory(outputDir, false);

	const std::string irMethod         = taskPM->value<std::string>("impulse_response_method");
	const FloatType sourceWidth        = taskPM->value<FloatType>("source_width", 0.0, 10.0);
	const FloatType sourceHeight       = taskPM->value<FloatType>("source_height", 0.0, 10.0);
	const FloatType propagationSpeed   = taskPM->value<FloatType>("propagation_speed", 0.0, 100000.0);
	const FloatType centerFreq         = taskPM->value<FloatType>("center_frequency", 0.0, 100.0e6);
	const FloatType maxFreq            = taskPM->value<FloatType>("max_frequency", 0.0, 200.0e6);
	const FloatType nyquistRate = 2.0 * maxFreq;
	const FloatType samplingFreq       = taskPM->value<FloatType>("sampling_frequency_factor", 0.0, 10000.0) * nyquistRate;
	const std::string excitationType   = taskPM->value<std::string>("excitation_type");
	const FloatType excNumPeriods      = taskPM->value<FloatType>("excitation_num_periods", 0.0, 100.0);
	const FloatType propagDistanceStep = taskPM->value<FloatType>("propagation_distance_step", 1.0e-6, 100.0);
	const FloatType propagMinDistance  = taskPM->value<FloatType>("propagation_min_distance", 0.0, 100.0);
	const FloatType propagMaxDistance  = taskPM->value<FloatType>("propagation_max_distance", propagMinDistance + propagDistanceStep, 100.0);
	const unsigned int propagRepet     = taskPM->value<unsigned int>("propagation_repetitions", 1, 100);
	const unsigned int propagPause     = taskPM->value<unsigned int>("propagation_pause", 0, 10000);

	std::vector<FloatType> propagDist;
	Util::fillSequenceFromStartWithStep(propagDist, propagMinDistance, propagMaxDistance, propagDistanceStep);
	std::vector<unsigned int> propagIndexList(propagDist.size());
	const FloatType coef = samplingFreq / propagationSpeed;
	for (unsigned int i = 0, end = propagIndexList.size(); i < end; ++i) {
		propagIndexList[i] = std::rint(propagDist[i] * coef);
	}

	std::vector<XY<FloatType>> elemPos;
	std::vector<FloatType> focusDelay;
	if (sourceIsArray) {
		ConstParameterMapPtr arrayPM = project_.loadChildParameterMap(taskPM, "array_config_file");
		ArrayUtil::calculateTxElementPositions(*arrayPM, elemPos);
		ArrayUtil::calculateTx3DFocusDelay(*taskPM, propagationSpeed, elemPos, focusDelay);
	}

	std::vector<FloatType> exc;
	if (excitationType == "1") {
		Waveform::getType1(centerFreq, samplingFreq, exc, excNumPeriods);
	} else if (excitationType == "2a") {
		Waveform::getType2a(centerFreq, samplingFreq, exc, excNumPeriods);
	} else if (excitationType == "2b") {
		Waveform::getType2b(centerFreq, samplingFreq, exc, excNumPeriods);
	} else {
		THROW_EXCEPTION(InvalidParameterException, "Invalid excitation type: " << excitationType << '.');
	}

	const FloatType dt = 1.0 / samplingFreq;
	std::vector<FloatType> tExc(exc.size());
	for (unsigned int i = 0; i < tExc.size(); ++i) {
		tExc[i] = dt * i;
	}
	project_.showFigure2D(1, "Excitation", tExc, exc);

	std::vector<FloatType> dvdt;
	Util::centralDiff(exc, dt, dvdt);

	Matrix<XYZValueArray<FloatType>> gridData;

	const FloatType nyquistLambda = propagationSpeed / nyquistRate;
	ImageGrid<FloatType>::get(project_.loadChildParameterMap(taskPM, "grid_config_file"), nyquistLambda, gridData);

	if (irMethod == "numeric") {
		const FloatType subElemSize = propagationSpeed /
				(nyquistRate * taskPM->value<FloatType>("sub_elem_size_factor", 0.0, 1.0e3));
		auto propag = std::make_unique<SimTransientPropagation<FloatType, NumericRectangularFlatSourceImpulseResponse<FloatType>>>();
		if (sourceIsArray) {
			propag->getArrayOfRectangularFlatSourcesPropagation(
						samplingFreq, propagationSpeed, sourceWidth, sourceHeight,
						subElemSize,
						dvdt, elemPos, focusDelay, propagIndexList, gridData);
		} else {
			propag->getRectangularFlatSourcePropagation(
						samplingFreq, propagationSpeed, sourceWidth, sourceHeight,
						subElemSize,
						dvdt, propagIndexList, gridData);
		}
	} else if (irMethod == "analytic") {
		const FloatType minEdgeDivisor = taskPM->value<FloatType>("min_edge_divisor", 0.0, 1.0e6);
		auto propag = std::make_unique<SimTransientPropagation<FloatType, AnalyticRectangularFlatSourceImpulseResponse<FloatType>>>();
		if (sourceIsArray) {
			propag->getArrayOfRectangularFlatSourcesPropagation(
						samplingFreq, propagationSpeed, sourceWidth, sourceHeight,
						minEdgeDivisor,
						dvdt, elemPos, focusDelay, propagIndexList, gridData);
		} else {
			propag->getRectangularFlatSourcePropagation(
						samplingFreq, propagationSpeed, sourceWidth, sourceHeight,
						minEdgeDivisor,
						dvdt, propagIndexList, gridData);
		}
	} else {
		THROW_EXCEPTION(InvalidParameterException, "Invalid impulse response method: " << irMethod << '.');
	}

	// Normalize.
	const FloatType maxAbsValue = Util::maxAbsoluteValueField(gridData);
	const FloatType k = 1.0 / maxAbsValue;
	for (auto it = gridData.begin(); it != gridData.end(); ++it) {
		for (auto& v : it->values) {
			v *= k;
		}
	}

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

	project_.saveHDF5(exc , outputDir + "/excitation"     , "value");
	project_.saveHDF5(tExc, outputDir + "/excitation_time", "value");

	// Save images.
	LOG_DEBUG << "Saving the X coordinates...";
	project_.saveHDF5(gridData, outputDir + "/image_x", "x", Util::CopyXOp());
	LOG_DEBUG << "Saving the Y coordinates...";
	project_.saveHDF5(gridData, outputDir + "/image_y", "y", Util::CopyYOp());
	LOG_DEBUG << "Saving the Z coordinates...";
	project_.saveHDF5(gridData, outputDir + "/image_z", "z", Util::CopyZOp());
	std::string imagePrefix = outputDir + "/image_value-";
	for (unsigned int i = 0, end = propagIndexList.size(); i < end; ++i) {
		LOG_DEBUG << "Saving image " << i << "...";
		project_.saveHDF5(gridData, imagePrefix + std::to_string(i), "value",
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

	const std::string outputDir = taskPM->value<std::string>("output_dir");
	project_.createDirectory(outputDir, false);

	const std::string irMethod       = taskPM->value<std::string>("impulse_response_method");
	const FloatType sourceWidth      = taskPM->value<FloatType>("source_width", 0.0, 10.0);
	const FloatType sourceHeight     = taskPM->value<FloatType>("source_height", 0.0, 10.0);
	const FloatType propagationSpeed = taskPM->value<FloatType>("propagation_speed", 0.0, 100000.0);
	//const FloatType density          = taskPM->value<FloatType>("density", 0.0, 100000.0);
	const FloatType centerFreq       = taskPM->value<FloatType>("center_frequency", 0.0, 100.0e6);
	const FloatType maxFreq          = taskPM->value<FloatType>("max_frequency", 0.0, 200.0e6);
	const FloatType nyquistRate = 2.0 * maxFreq;
	const FloatType samplingFreq     = taskPM->value<FloatType>("sampling_frequency_factor", 0.0, 10000.0) * nyquistRate;
	const std::string excitationType = taskPM->value<std::string>("excitation_type");
	const FloatType excNumPeriods    = taskPM->value<FloatType>("excitation_num_periods", 0.0, 100.0);
	const FloatType pointX           = taskPM->value<FloatType>("point_x", 0.0, 10000.0);
	const FloatType pointY           = taskPM->value<FloatType>("point_y", 0.0, 10000.0);
	const FloatType pointZ           = taskPM->value<FloatType>("point_z", 0.0, 10000.0);

	std::vector<XY<FloatType>> elemPos;
	std::vector<FloatType> focusDelay;
	if (sourceIsArray) {
		ConstParameterMapPtr arrayPM = project_.loadChildParameterMap(taskPM, "array_config_file");
		ArrayUtil::calculateTxElementPositions(*arrayPM, elemPos);
		ArrayUtil::calculateTx3DFocusDelay(*taskPM, propagationSpeed, elemPos, focusDelay);
	}

	std::vector<FloatType> exc;
	if (excitationType == "1") {
		Waveform::getType1(centerFreq, samplingFreq, exc, excNumPeriods);
	} else if (excitationType == "2a") {
		Waveform::getType2a(centerFreq, samplingFreq, exc, excNumPeriods);
	} else if (excitationType == "2b") {
		Waveform::getType2b(centerFreq, samplingFreq, exc, excNumPeriods);
	} else {
		THROW_EXCEPTION(InvalidParameterException, "Invalid excitation type: " << excitationType << '.');
	}

	const FloatType dt = 1.0 / samplingFreq;
	std::vector<FloatType> tExc;
	Util::fillSequenceFromStartWithStepAndSize(tExc, 0.0, dt, exc.size());
	project_.showFigure2D(1, "v", tExc, exc);

	std::vector<FloatType> dvdt;
	Util::centralDiff(exc, dt, dvdt);
	Util::normalizeBySumOfAbs(dvdt);
	std::vector<FloatType> tDvdt;
	Util::fillSequenceFromStartWithStepAndSize(tDvdt, 0.0, dt, dvdt.size());
	project_.showFigure2D(2, "dv/dt", tDvdt, dvdt);

	std::size_t hOffset;
	std::vector<FloatType> h;
	if (irMethod == "numeric") {
		const FloatType subElemSize = propagationSpeed /
				(nyquistRate * taskPM->value<FloatType>("sub_elem_size_factor", 0.0, 1.0e3));
		if (sourceIsArray) {
			auto impResp = std::make_unique<ArrayOfRectangularFlatSourcesImpulseResponse<FloatType, NumericRectangularFlatSourceImpulseResponse<FloatType>>>(
						samplingFreq, propagationSpeed, sourceWidth, sourceHeight,
						subElemSize,
						elemPos, focusDelay);
			impResp->getImpulseResponse(pointX, pointY, pointZ, hOffset, h);
		} else {
			auto impResp = std::make_unique<NumericRectangularFlatSourceImpulseResponse<FloatType>>(
						samplingFreq, propagationSpeed, sourceWidth, sourceHeight,
						subElemSize);
			impResp->getImpulseResponse(pointX, pointY, pointZ, hOffset, h);
		}
	} else if (irMethod == "analytic") {
		const FloatType minEdgeDivisor = taskPM->value<FloatType>("min_edge_divisor", 0.0, 1.0e6);
		if (sourceIsArray) {
			auto impResp = std::make_unique<ArrayOfRectangularFlatSourcesImpulseResponse<FloatType, AnalyticRectangularFlatSourceImpulseResponse<FloatType>>>(
						samplingFreq, propagationSpeed, sourceWidth, sourceHeight,
						minEdgeDivisor,
						elemPos, focusDelay);
			impResp->getImpulseResponse(pointX, pointY, pointZ, hOffset, h);
		} else {
			auto impResp = std::make_unique<AnalyticRectangularFlatSourceImpulseResponse<FloatType>>(
						samplingFreq, propagationSpeed, sourceWidth, sourceHeight,
						minEdgeDivisor);
			impResp->getImpulseResponse(pointX, pointY, pointZ, hOffset, h);
		}
	} else {
		THROW_EXCEPTION(InvalidParameterException, "Invalid impulse response method: " << irMethod << '.');
	}

	std::vector<FloatType> tH;
	Util::fillSequenceFromStartWithStepAndSize(tH, hOffset / samplingFreq, dt, h.size());
	project_.showFigure2D(3, "Impulse response", tH, h);

	std::vector<std::complex<FloatType>> filterFreqCoeff;
	auto filter = std::make_unique<FFTWFilter2<FloatType>>();
	filter->setCoefficients(dvdt, filterFreqCoeff);
	std::vector<FloatType> signal;

	filter->filter(filterFreqCoeff, h, signal);

	std::vector<FloatType> tSignal;
	Util::fillSequenceFromStartWithStepAndSize(tSignal, hOffset / samplingFreq, dt, signal.size());
	//Util::multiply(signal, density);
	project_.showFigure2D(4, "Pressure", tSignal, signal);

	project_.saveHDF5(exc    , outputDir + "/excitation"           , "value");
	project_.saveHDF5(tExc   , outputDir + "/excitation_time"      , "value");
	project_.saveHDF5(h      , outputDir + "/impulse_response"     , "value");
	project_.saveHDF5(tH     , outputDir + "/impulse_response_time", "value");
	project_.saveHDF5(signal , outputDir + "/pressure"             , "value");
	project_.saveHDF5(tSignal, outputDir + "/pressure_time"        , "value");
}

template<typename FloatType>
void
SimRectangularFlatSourceMethod<FloatType>::execute()
{
	Timer tProc;

	switch (project_.method()) {
	case MethodType::sim_acoustic_field_rectangular_flat_source_transient:
		execTransientAcousticField();
		break;
	case MethodType::sim_acoustic_field_array_of_rectangular_flat_sources_transient:
		execTransientArrayAcousticField();
		break;
	case MethodType::sim_impulse_response_rectangular_flat_source:
		execImpulseResponse(false);
		break;
	case MethodType::sim_impulse_response_array_of_rectangular_flat_sources:
		execImpulseResponse(true);
		break;
	case MethodType::sim_propagation_rectangular_flat_source_transient:
		execTransientPropagation(false);
		break;
	case MethodType::sim_propagation_array_of_rectangular_flat_sources_transient:
		execTransientPropagation(true);
		break;
	case MethodType::sim_radiation_pattern_rectangular_flat_source_transient:
		execTransientRadiationPattern(false);
		break;
	case MethodType::sim_radiation_pattern_array_of_rectangular_flat_sources_transient:
		execTransientRadiationPattern(true);
		break;
	default:
		THROW_EXCEPTION(InvalidParameterException, "Invalid method: " << static_cast<int>(project_.method()) << '.');
	}

	LOG_DEBUG << ">>> Processing time: " << tProc.getTime();
}

} // namespace Lab

#endif // SIMRECTANGULARFLATSOURCEMETHOD_H
