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
#include <limits>
#include <memory>

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
#include "XZ.h"
#include "XZValue.h"



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

	void calculateArrayParameters(
			const ParameterMap& taskPM,
			FloatType propagationSpeed,
			FloatType samplingFreq,
			std::vector<XY<FloatType>>& elemPos,
			std::vector<FloatType>& focusDelay);

	void execTransientAcousticBeam();
	void execTransientArrayAcousticBeam();
	void execTransientAcousticField();
	void execTransientArrayAcousticField();
	void execImpulseResponse();
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
SimRectangularFlatSourceMethod<FloatType>::calculateArrayParameters(
		const ParameterMap& taskPM,
		FloatType propagationSpeed,
		FloatType samplingFreq,
		std::vector<XY<FloatType>>& elemPos,
		std::vector<FloatType>& focusDelay)
{
	const bool useFocus              = taskPM.value<bool>("use_focus");
	const FloatType focusX           = taskPM.value<FloatType>("focus_x", 0.0, 10000.0);
	const FloatType focusY           = taskPM.value<FloatType>("focus_y", 0.0, 10000.0);
	const FloatType focusZ           = taskPM.value<FloatType>("focus_z", 0.0, 10000.0);
	const FloatType arrayPitchX      = taskPM.value<FloatType>("array_pitch_x", 0.0, 10.0);
	const unsigned int numArrayElemX = taskPM.value<FloatType>("num_array_elem_x", 1, 1024);
	const FloatType arrayPitchY      = taskPM.value<FloatType>("array_pitch_y", 0.0, 10.0);
	const unsigned int numArrayElemY = taskPM.value<FloatType>("num_array_elem_y", 1, 1024);

	const unsigned int numElem = numArrayElemX * numArrayElemY;

	// Calculate the center of each element.
	const FloatType halfW = (numArrayElemX - 1) * 0.5 * arrayPitchX;
	const FloatType halfH = (numArrayElemY - 1) * 0.5 * arrayPitchY;
	elemPos.resize(numElem);
	for (unsigned int iy = 0; iy < numArrayElemY; ++iy) {
		for (unsigned int ix = 0; ix < numArrayElemX; ++ix) {
			XY<FloatType>& pos = elemPos[iy * numArrayElemX + ix];
			pos.x = ix * arrayPitchX - halfW;
			pos.y = iy * arrayPitchY - halfH;
		}
	}

	if (useFocus) {
		focusDelay.resize(numElem);
		FloatType maxDt = 0.0;
		for (unsigned int iy = 0; iy < numArrayElemY; ++iy) {
			for (unsigned int ix = 0; ix < numArrayElemX; ++ix) {
				const unsigned int index = iy * numArrayElemX + ix;
				XY<FloatType>& pos = elemPos[index];
				const FloatType dx = focusX - pos.x;
				const FloatType dy = focusY - pos.y;
				const FloatType focusDt = std::sqrt(dx * dx + dy * dy + focusZ * focusZ) / propagationSpeed;
				focusDelay[index] = focusDt;
				if (focusDt > maxDt) maxDt = focusDt;
			}
		}
		for (unsigned int iy = 0; iy < numArrayElemY; ++iy) {
			for (unsigned int ix = 0; ix < numArrayElemX; ++ix) {
				const unsigned int index = iy * numArrayElemX + ix;
				focusDelay[index] = (maxDt - focusDelay[index]) * samplingFreq;
			}
		}
	} else {
		focusDelay.assign(numElem, 0.0);
	}
}

template<typename FloatType>
void
SimRectangularFlatSourceMethod<FloatType>::execTransientAcousticBeam()
{
	ConstParameterMapPtr taskPM = project_.taskParameterMap();

	const std::string outputDir = taskPM->value<std::string>("output_dir");

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
	const FloatType subElemSize      = propagationSpeed / (nyquistRate * taskPM->value<FloatType>("sub_elem_size_factor", 0.0, 1000.0));
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

	Matrix2<XZValue<FloatType>> gridData{thetaXList.size(), thetaYList.size()};
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
			XZValue<FloatType>& gd = gridData(ix, iy);
			gd.x = thetaYList[iy];
			gd.z = thetaXList[ix];
			gd.value = 0.0;
			XYZ<FloatType>& id = inputData(ix, iy);
			id.x = x;
			id.y = y;
			id.z = z;
		}
	}

	auto acBeam = std::make_unique<SimTransientAcousticBeam<FloatType>>();
	acBeam->getRectangularFlatSourceAcousticBeam(
				sourceWidth,
				sourceHeight,
				samplingFreq,
				propagationSpeed,
				subElemSize,
				dvdt,
				inputData,
				gridData);

	const FloatType maxAbsValue = Util::maxAbsoluteValueField<XZValue<FloatType>, FloatType>(gridData);
	const FloatType k = 1.0 / maxAbsValue;
	for (auto it = gridData.begin(); it != gridData.end(); ++it) {
		it->value *= k;
	}

	std::vector<XZ<float>> pointList = {{0.0, 0.0}};

	Project::GridDataType projGridData;
	Util::copyXZValue(gridData, projGridData);
	project_.showFigure3D(1, "Beam", &projGridData, &pointList,
					true, Figure::VISUALIZATION_RECTIFIED_LOG, Figure::COLORMAP_VIRIDIS);
}

template<typename FloatType>
void
SimRectangularFlatSourceMethod<FloatType>::execTransientArrayAcousticBeam()
{

}

template<typename FloatType>
void
SimRectangularFlatSourceMethod<FloatType>::execTransientAcousticField()
{
	ConstParameterMapPtr taskPM = project_.taskParameterMap();

	const std::string outputDir = taskPM->value<std::string>("output_dir");

	const FloatType sourceWidth      = taskPM->value<FloatType>("source_width", 0.0, 10.0);
	const FloatType sourceHeight     = taskPM->value<FloatType>("source_height", 0.0, 10.0);
	const FloatType propagationSpeed = taskPM->value<FloatType>("propagation_speed", 0.0, 100000.0);
	const FloatType centerFreq       = taskPM->value<FloatType>("center_frequency", 0.0, 100.0e6);
	const FloatType maxFreq          = taskPM->value<FloatType>("max_frequency", 0.0, 200.0e6);
	const FloatType nyquistRate = 2.0 * maxFreq;
	const FloatType samplingFreq     = taskPM->value<FloatType>("sampling_frequency_factor", 0.0, 10000.0) * nyquistRate;
	const FloatType subElemSize      = propagationSpeed / (nyquistRate * taskPM->value<FloatType>("sub_elem_size_factor", 0.0, 1000.0));
	const std::string excitationType = taskPM->value<std::string>("excitation_type");
	const FloatType excNumPeriods    = taskPM->value<FloatType>("excitation_num_periods", 0.0, 100.0);
	const FloatType y                = taskPM->value<FloatType>("y", -1000.0, 1000.0);

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

	Matrix2<XZValue<FloatType>> gridData;

	const FloatType nyquistLambda = propagationSpeed / nyquistRate;
	ImageGrid<FloatType>::get(project_.loadChildParameterMap(taskPM, "grid_config_file"), nyquistLambda, gridData);

	auto acField = std::make_unique<SimTransientAcousticField<FloatType>>();
	acField->getRectangularFlatSourceAcousticField(
				sourceWidth,
				sourceHeight,
				samplingFreq,
				propagationSpeed,
				subElemSize,
				dvdt,
				y,
				gridData);

	const FloatType maxAbsValue = Util::maxAbsoluteValueField<XZValue<FloatType>, FloatType>(gridData);
	const FloatType k = 1.0 / maxAbsValue;
	for (auto it = gridData.begin(); it != gridData.end(); ++it) {
		it->value *= k;
	}

	std::vector<XZ<float>> pointList = {{0.0, 0.0}};

	Project::GridDataType projGridData;
	Util::copyXZValue(gridData, projGridData);
	project_.showFigure3D(1, "Acoustic field", &projGridData, &pointList,
					true, Figure::VISUALIZATION_RECTIFIED_LINEAR, Figure::COLORMAP_VIRIDIS);
}

template<typename FloatType>
void
SimRectangularFlatSourceMethod<FloatType>::execTransientArrayAcousticField()
{

}

template<typename FloatType>
void
SimRectangularFlatSourceMethod<FloatType>::execImpulseResponse()
{
	ConstParameterMapPtr taskPM = project_.taskParameterMap();

	const std::string outputDir = taskPM->value<std::string>("output_dir");

	const FloatType sourceWidth      = taskPM->value<FloatType>("source_width", 0.0, 10.0);
	const FloatType sourceHeight     = taskPM->value<FloatType>("source_height", 0.0, 10.0);
	const FloatType propagationSpeed = taskPM->value<FloatType>("propagation_speed", 0.0, 100000.0);
	const FloatType density          = taskPM->value<FloatType>("density", 0.0, 100000.0);
	const FloatType centerFreq       = taskPM->value<FloatType>("center_frequency", 0.0, 100.0e6);
	const FloatType maxFreq          = taskPM->value<FloatType>("max_frequency", 0.0, 200.0e6);
	const FloatType nyquistRate = 2.0 * maxFreq;
	const FloatType samplingFreq     = taskPM->value<FloatType>("sampling_frequency_factor", 0.0, 10000.0) * nyquistRate;
	const FloatType subElemSize      = propagationSpeed / (nyquistRate * taskPM->value<FloatType>("sub_elem_size_factor", 0.0, 1000.0));
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
	project_.showFigure2D(1, "Excitation", tExc, exc);

	std::vector<FloatType> dvdt;
	Util::centralDiff(exc, dt, dvdt);

	std::size_t hOffset;
	std::vector<FloatType> h;
	auto impResp = std::make_unique<NumericRectangularFlatSourceImpulseResponse<FloatType>>(
									sourceWidth,
									sourceHeight,
									samplingFreq,
									propagationSpeed,
									subElemSize);
	impResp->getImpulseResponse(pointX, pointY, pointZ, hOffset, h);

	std::vector<FloatType> tH;
	Util::fillSequenceFromStartWithStepAndSize(tH, hOffset / samplingFreq, dt, h.size());
	project_.showFigure2D(2, "Impulse response", tH, h);

	std::vector<std::complex<FloatType>> filterFreqCoeff;
	auto filter = std::make_unique<FFTWFilter2<FloatType>>();
	filter->setCoefficients(dvdt, filterFreqCoeff);
	std::vector<FloatType> signal;

	filter->filter(filterFreqCoeff, h, signal);

	std::vector<FloatType> tSignal;
	Util::fillSequenceFromStartWithStepAndSize(tSignal, hOffset / samplingFreq, dt, signal.size());
	Util::multiply(signal, density);
	project_.showFigure2D(3, "Pressure", tSignal, signal);
}

template<typename FloatType>
void
SimRectangularFlatSourceMethod<FloatType>::execArrayImpulseResponse()
{
	ConstParameterMapPtr taskPM = project_.taskParameterMap();

	const std::string outputDir = taskPM->value<std::string>("output_dir");

	const FloatType sourceWidth      = taskPM->value<FloatType>("source_width", 0.0, 10.0);
	const FloatType sourceHeight     = taskPM->value<FloatType>("source_height", 0.0, 10.0);
	const FloatType propagationSpeed = taskPM->value<FloatType>("propagation_speed", 0.0, 100000.0);
	const FloatType density          = taskPM->value<FloatType>("density", 0.0, 100000.0);
	const FloatType centerFreq       = taskPM->value<FloatType>("center_frequency", 0.0, 100.0e6);
	const FloatType maxFreq          = taskPM->value<FloatType>("max_frequency", 0.0, 200.0e6);
	const FloatType nyquistRate = 2.0 * maxFreq;
	const FloatType samplingFreq     = taskPM->value<FloatType>("sampling_frequency_factor", 0.0, 10000.0) * nyquistRate;
	const FloatType subElemSize      = propagationSpeed / (nyquistRate * taskPM->value<FloatType>("sub_elem_size_factor", 0.0, 1000.0));
	const std::string excitationType = taskPM->value<std::string>("excitation_type");
	const FloatType excNumPeriods    = taskPM->value<FloatType>("excitation_num_periods", 0.0, 100.0);
	const FloatType pointX           = taskPM->value<FloatType>("point_x", 0.0, 10000.0);
	const FloatType pointY           = taskPM->value<FloatType>("point_y", 0.0, 10000.0);
	const FloatType pointZ           = taskPM->value<FloatType>("point_z", 0.0, 10000.0);

	std::vector<XY<FloatType>> elemPos;
	std::vector<FloatType> focusDelay;
	calculateArrayParameters(*taskPM, propagationSpeed, samplingFreq, elemPos, focusDelay);

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
	project_.showFigure2D(1, "Excitation", tExc, exc);

	std::vector<FloatType> dvdt;
	Util::centralDiff(exc, dt, dvdt);

	auto impResp = std::make_unique<NumericRectangularFlatSourceImpulseResponse<FloatType>>(
									sourceWidth,
									sourceHeight,
									samplingFreq,
									propagationSpeed,
									subElemSize);

	std::vector<std::size_t> offsetList(elemPos.size());
	std::vector<std::vector<FloatType>> hList(elemPos.size());
	std::size_t iMin = std::numeric_limits<std::size_t>::max();
	std::size_t iMax = 0; // the index after the last

	// Obtain the impulse responses.
	for (std::size_t i = 0, iEnd = elemPos.size(); i < iEnd; ++i) {
		XY<FloatType>& pos = elemPos[i];

		impResp->getImpulseResponse(pointX - pos.x, pointY - pos.y, pointZ, offsetList[i], hList[i]);

		offsetList[i] += static_cast<std::size_t>(std::nearbyint(focusDelay[i]));
		const std::size_t start = offsetList[i];
		if (start < iMin) iMin = start;
		const std::size_t end = offsetList[i] + hList[i].size();
		if (end > iMax) iMax = end;
	}

	// Accumulate the impulse responses.
	std::vector<FloatType> hSum(iMax - iMin);
	for (std::size_t i = 0, iEnd = elemPos.size(); i < iEnd; ++i) {
		const std::size_t hBegin = offsetList[i];
		const std::size_t hEnd = hBegin + hList[i].size();
		Util::addElements(hList[i].begin(), hSum.begin() + hBegin - iMin, hSum.begin() + hEnd - iMin);
	}
	const std::size_t hOffset = iMin;

	std::vector<FloatType> tH;
	Util::fillSequenceFromStartWithStepAndSize(tH, hOffset / samplingFreq, dt, hSum.size());
	project_.showFigure2D(2, "Impulse response", tH, hSum);

	std::vector<std::complex<FloatType>> filterFreqCoeff;
	auto filter = std::make_unique<FFTWFilter2<FloatType>>();
	filter->setCoefficients(dvdt, filterFreqCoeff);
	std::vector<FloatType> signal;

	filter->filter(filterFreqCoeff, hSum, signal);

	std::vector<FloatType> tSignal;
	Util::fillSequenceFromStartWithStepAndSize(tSignal, hOffset / samplingFreq, dt, signal.size());
	Util::multiply(signal, density);
	project_.showFigure2D(3, "Pressure", tSignal, signal);
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
