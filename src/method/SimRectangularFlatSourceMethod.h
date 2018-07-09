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

#include "AnalyticRectangularFlatSourceImpulseResponse.h"
#include "ArrayOfRectangularFlatSourcesImpulseResponse.h"
#include "ArrayUtil.h"
#include "Exception.h"
#include "FFTWFilter2.h"
#include "ImageGrid.h"
#include "Log.h"
#include "Matrix2.h"
#include "Method.h"
#include "NumericRectangularFlatSourceImpulseResponse.h"
#include "ParameterMap.h"
#include "Project.h"
#include "SimTransientAcousticBeam.h"
#include "SimTransientAcousticField.h"
#include "Timer.h"
#include "Util.h"
#include "Waveform.h"
#include "XY.h"
#include "XYZ.h"
#include "XYZValue.h"



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

	void execTransientAcousticBeam();
	void execTransientArrayAcousticBeam();
	void execTransientAcousticField();
	void execTransientArrayAcousticField();
	// Return p/(c*density).
	void execImpulseResponse();
	// Return p/(c*density).
	void execArrayImpulseResponse();

	Project& project_;
};



template<typename FloatType>
SimRectangularFlatSourceMethod<FloatType>::SimRectangularFlatSourceMethod(Project& project)
		: project_{project}
{
}

template<typename FloatType>
SimRectangularFlatSourceMethod<FloatType>::~SimRectangularFlatSourceMethod()
{
}

template<typename FloatType>
void
SimRectangularFlatSourceMethod<FloatType>::execTransientAcousticBeam()
{
	ConstParameterMapPtr taskPM = project_.taskParameterMap();

	const std::string outputDir = taskPM->value<std::string>("output_dir");
	project_.createDirectory(outputDir, false);

	const std::string irMethod       = taskPM->value<std::string>("impulse_response_method");
	const FloatType beamDistance     = taskPM->value<FloatType>("beam_distance", 0.0, 100.0);
	const FloatType beamThetaYStep   = taskPM->value<FloatType>("beam_theta_y_step", 0.0, 10.0);
	const FloatType beamThetaYMax    = taskPM->value<FloatType>("beam_theta_y_max" , 0.0, 90.0);
	const FloatType beamThetaXStep   = taskPM->value<FloatType>("beam_theta_x_step", 0.0, 10.0);
	const FloatType beamThetaXMax    = taskPM->value<FloatType>("beam_theta_x_max" , 0.0, 90.0);
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

	std::vector<FloatType> thetaXList;
	Util::fillSequenceFromStartToEndWithSize(thetaXList, 0.0, beamThetaXMax, std::ceil(beamThetaXMax / beamThetaXStep) + 1);
	std::vector<FloatType> thetaYList;
	Util::fillSequenceFromStartToEndWithSize(thetaYList, 0.0, beamThetaYMax, std::ceil(beamThetaYMax / beamThetaYStep) + 1);

	Matrix2<XYZValue<FloatType>> gridData{thetaXList.size(), thetaYList.size()};
	Matrix2<XYZ<FloatType>> inputData{thetaXList.size(), thetaYList.size()};

	for (unsigned int ix = 0, xSize = thetaXList.size(); ix < xSize; ++ix) {
		const FloatType tX = Util::degreeToRadian(thetaXList[ix]);
		const FloatType sinTX = std::sin(tX);
		const FloatType cosTX = std::cos(tX);
		for (unsigned int iy = 0, ySize = thetaYList.size(); iy < ySize; ++iy) {
			const FloatType tY = Util::degreeToRadian(thetaYList[iy]);
			const FloatType x = beamDistance * std::sin(tY);
			const FloatType rx = beamDistance * std::cos(tY);
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
		auto acBeam = std::make_unique<SimTransientAcousticBeam<FloatType, NumericRectangularFlatSourceImpulseResponse<FloatType>>>();
		acBeam->getRectangularFlatSourceAcousticBeam(
					samplingFreq, propagationSpeed, sourceWidth, sourceHeight,
					subElemSize,
					dvdt, inputData, gridData);
	} else if (irMethod == "analytic") {
		const FloatType minEdgeDivisor = taskPM->value<FloatType>("min_edge_divisor", 0.0, 1.0e6);
		auto acBeam = std::make_unique<SimTransientAcousticBeam<FloatType, AnalyticRectangularFlatSourceImpulseResponse<FloatType>>>();
		acBeam->getRectangularFlatSourceAcousticBeam(
					samplingFreq, propagationSpeed, sourceWidth, sourceHeight,
					minEdgeDivisor,
					dvdt, inputData, gridData);
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
	project_.showFigure3D(1, "Beam", &projGridData, &pointList,
					true, Figure::VISUALIZATION_RECTIFIED_LOG, Figure::COLORMAP_VIRIDIS);

	std::vector<FloatType> beamTY(gridData.n2());
	auto tyInterval = gridData.dim2Interval(0);
	Util::copyUsingOperator(tyInterval.first, tyInterval.second, beamTY.begin(), Util::CopyValueOp{});
	Util::linearToDecibels(beamTY, -100.0);
	project_.showFigure2D(1, "Beam theta-y", thetaYList, beamTY);

	std::vector<FloatType> beamTX(gridData.n1());
	auto txInterval = gridData.dim1Interval(0);
	Util::copyUsingOperator(txInterval.first, txInterval.second, beamTX.begin(), Util::CopyValueOp{});
	Util::linearToDecibels(beamTX, -100.0);
	project_.showFigure2D(2, "Beam theta-x", thetaXList, beamTX);

	project_.saveHDF5(exc , outputDir + "/excitation"     , "value");
	project_.saveHDF5(tExc, outputDir + "/excitation_time", "value");
	project_.saveImageToHDF5(gridData, outputDir);
	project_.saveHDF5(thetaYList, outputDir + "/theta_y"     , "value");
	project_.saveHDF5(beamTY    , outputDir + "/beam_theta_y", "value");
	project_.saveHDF5(thetaXList, outputDir + "/theta_x"     , "value");
	project_.saveHDF5(beamTX    , outputDir + "/beam_theta_x", "value");
}

template<typename FloatType>
void
SimRectangularFlatSourceMethod<FloatType>::execTransientArrayAcousticBeam()
{
	ConstParameterMapPtr taskPM = project_.taskParameterMap();

	const std::string outputDir = taskPM->value<std::string>("output_dir");
	project_.createDirectory(outputDir, false);

	const std::string irMethod       = taskPM->value<std::string>("impulse_response_method");
	const FloatType beamDistance     = taskPM->value<FloatType>("beam_distance", 0.0, 100.0);
	const FloatType beamThetaYStep   = taskPM->value<FloatType>("beam_theta_y_step", 0.0, 10.0);
	const FloatType beamThetaYMin    = taskPM->value<FloatType>("beam_theta_y_min" , -90.0, 90.0);
	const FloatType beamThetaYMax    = taskPM->value<FloatType>("beam_theta_y_max" , beamThetaYMin + 0.1, 90.0);
	const FloatType beamThetaXStep   = taskPM->value<FloatType>("beam_theta_x_step", 0.0, 10.0);
	const FloatType beamThetaXMin    = taskPM->value<FloatType>("beam_theta_x_min" , -90.0, 90.0);
	const FloatType beamThetaXMax    = taskPM->value<FloatType>("beam_theta_x_max" , beamThetaXMin + 0.1, 90.0);
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

	std::vector<FloatType> thetaXList;
	Util::fillSequenceFromStartToEndWithSize(thetaXList, beamThetaXMin, beamThetaXMax, std::ceil((beamThetaXMax - beamThetaXMin) / beamThetaXStep) + 1);
	std::vector<FloatType> thetaYList;
	Util::fillSequenceFromStartToEndWithSize(thetaYList, beamThetaYMin, beamThetaYMax, std::ceil((beamThetaYMax - beamThetaYMin) / beamThetaYStep) + 1);

	Matrix2<XYZValue<FloatType>> gridData{thetaXList.size(), thetaYList.size()};
	Matrix2<XYZ<FloatType>> inputData{thetaXList.size(), thetaYList.size()};

	for (unsigned int ix = 0, xSize = thetaXList.size(); ix < xSize; ++ix) {
		const FloatType tX = Util::degreeToRadian(thetaXList[ix]);
		const FloatType sinTX = std::sin(tX);
		const FloatType cosTX = std::cos(tX);
		for (unsigned int iy = 0, ySize = thetaYList.size(); iy < ySize; ++iy) {
			const FloatType tY = Util::degreeToRadian(thetaYList[iy]);
			const FloatType x = beamDistance * std::sin(tY);
			const FloatType rx = beamDistance * std::cos(tY);
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
		auto acBeam = std::make_unique<SimTransientAcousticBeam<FloatType, NumericRectangularFlatSourceImpulseResponse<FloatType>>>();
		acBeam->getArrayOfRectangularFlatSourcesAcousticBeam(
					samplingFreq, propagationSpeed, sourceWidth, sourceHeight,
					subElemSize,
					dvdt, elemPos, focusDelay, inputData, gridData);
	} else if (irMethod == "analytic") {
		const FloatType minEdgeDivisor = taskPM->value<FloatType>("min_edge_divisor", 0.0, 1.0e6);
		auto acBeam = std::make_unique<SimTransientAcousticBeam<FloatType, AnalyticRectangularFlatSourceImpulseResponse<FloatType>>>();
		acBeam->getArrayOfRectangularFlatSourcesAcousticBeam(
					samplingFreq, propagationSpeed, sourceWidth, sourceHeight,
					minEdgeDivisor,
					dvdt, elemPos, focusDelay, inputData, gridData);
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
	project_.showFigure3D(1, "Beam", &projGridData, &pointList,
					true, Figure::VISUALIZATION_RECTIFIED_LOG, Figure::COLORMAP_VIRIDIS);

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
	std::size_t sectionTYIndex = std::rint(((sectionTY - beamThetaYMin) / (beamThetaYMax - beamThetaYMin)) * (thetaYList.size() - 1));
	std::size_t sectionTXIndex = std::rint(((sectionTX - beamThetaXMin) / (beamThetaXMax - beamThetaXMin)) * (thetaXList.size() - 1));
	LOG_DEBUG << "Section theta-x: " << thetaXList[sectionTXIndex] << " theta-y: " << thetaYList[sectionTYIndex];

	std::vector<FloatType> beamTY(gridData.n2());
	auto tyInterval = gridData.dim2Interval(sectionTXIndex);
	Util::copyUsingOperator(tyInterval.first, tyInterval.second, beamTY.begin(), Util::CopyValueOp{});
	Util::linearToDecibels(beamTY, -100.0);
	project_.showFigure2D(1, "Beam theta-y", thetaYList, beamTY);

	std::vector<FloatType> beamTX(gridData.n1());
	auto txInterval = gridData.dim1Interval(sectionTYIndex);
	Util::copyUsingOperator(txInterval.first, txInterval.second, beamTX.begin(), Util::CopyValueOp{});
	Util::linearToDecibels(beamTX, -100.0);
	project_.showFigure2D(2, "Beam theta-x", thetaXList, beamTX);

	project_.saveHDF5(exc , outputDir + "/excitation"     , "value");
	project_.saveHDF5(tExc, outputDir + "/excitation_time", "value");
	project_.saveImageToHDF5(gridData, outputDir);
	project_.saveHDF5(thetaYList, outputDir + "/theta_y"     , "value");
	project_.saveHDF5(beamTY    , outputDir + "/beam_theta_y", "value");
	project_.saveHDF5(thetaXList, outputDir + "/theta_x"     , "value");
	project_.saveHDF5(beamTX    , outputDir + "/beam_theta_x", "value");
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

	Matrix2<XYZValue<FloatType>> gridData;

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

	Matrix2<XYZValue<FloatType>> gridData;

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
SimRectangularFlatSourceMethod<FloatType>::execImpulseResponse()
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
		auto impResp = std::make_unique<NumericRectangularFlatSourceImpulseResponse<FloatType>>(
					samplingFreq, propagationSpeed, sourceWidth, sourceHeight,
					subElemSize);
		impResp->getImpulseResponse(pointX, pointY, pointZ, hOffset, h);
	} else if (irMethod == "analytic") {
		const FloatType minEdgeDivisor = taskPM->value<FloatType>("min_edge_divisor", 0.0, 1.0e6);
		auto impResp = std::make_unique<AnalyticRectangularFlatSourceImpulseResponse<FloatType>>(
					samplingFreq, propagationSpeed, sourceWidth, sourceHeight,
					minEdgeDivisor);
		impResp->getImpulseResponse(pointX, pointY, pointZ, hOffset, h);
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
SimRectangularFlatSourceMethod<FloatType>::execArrayImpulseResponse()
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
		auto impResp = std::make_unique<ArrayOfRectangularFlatSourcesImpulseResponse<FloatType, NumericRectangularFlatSourceImpulseResponse<FloatType>>>(
					samplingFreq, propagationSpeed, sourceWidth, sourceHeight,
					subElemSize,
					elemPos, focusDelay);
		impResp->getImpulseResponse(pointX, pointY, pointZ, hOffset, h);
	} else if (irMethod == "analytic") {
		const FloatType minEdgeDivisor = taskPM->value<FloatType>("min_edge_divisor", 0.0, 1.0e6);
		auto impResp = std::make_unique<ArrayOfRectangularFlatSourcesImpulseResponse<FloatType, AnalyticRectangularFlatSourceImpulseResponse<FloatType>>>(
					samplingFreq, propagationSpeed, sourceWidth, sourceHeight,
					minEdgeDivisor,
					elemPos, focusDelay);
		impResp->getImpulseResponse(pointX, pointY, pointZ, hOffset, h);
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
	case MethodType::sim_acoustic_beam_rectangular_flat_source_transient:
		execTransientAcousticBeam();
		break;
	case MethodType::sim_acoustic_beam_array_of_rectangular_flat_sources_transient:
		execTransientArrayAcousticBeam();
		break;
	case MethodType::sim_acoustic_field_rectangular_flat_source_transient:
		execTransientAcousticField();
		break;
	case MethodType::sim_acoustic_field_array_of_rectangular_flat_sources_transient:
		execTransientArrayAcousticField();
		break;
	case MethodType::sim_impulse_response_rectangular_flat_source:
		execImpulseResponse();
		break;
	case MethodType::sim_impulse_response_array_of_rectangular_flat_sources:
		execArrayImpulseResponse();
		break;
	default:
		THROW_EXCEPTION(InvalidParameterException, "Invalid method: " << static_cast<int>(project_.method()) << '.');
	}

	LOG_DEBUG << ">>> Processing time: " << tProc.getTime();
}

} // namespace Lab

#endif // SIMRECTANGULARFLATSOURCEMETHOD_H
