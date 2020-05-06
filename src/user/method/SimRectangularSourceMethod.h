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
#include "NumericArrayOfRectangularSourcesImpulseResponse.h"
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
#ifdef USE_CUDA
# include "NumericArrayOfRectangularSourcesCUDAImpulseResponse.h"
# include "NumericRectangularSourceCUDAImpulseResponse.h"
#endif
#ifdef USE_OPENCL
# include "NumericArrayOfRectangularSourcesOCLImpulseResponse.h"
# include "NumericRectangularSourceOCLImpulseResponse.h"
#endif



namespace Lab {

template<typename TFloat>
class SimRectangularSourceMethod : public Method {
public:
	SimRectangularSourceMethod(Project& project);
	virtual ~SimRectangularSourceMethod() = default;

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
		TFloat sourceWidth;
		TFloat sourceHeight;

		// For arrays.
		TFloat focusX;
		TFloat focusY;
		TFloat focusZ;
		std::vector<XY<TFloat>> elemPos;
		std::vector<TFloat> focusDelay;
		bool useFocus;
	};

	SimRectangularSourceMethod(const SimRectangularSourceMethod&) = delete;
	SimRectangularSourceMethod& operator=(const SimRectangularSourceMethod&) = delete;
	SimRectangularSourceMethod(SimRectangularSourceMethod&&) = delete;
	SimRectangularSourceMethod& operator=(SimRectangularSourceMethod&&) = delete;

	void loadData(const ParameterMap& taskPM, MainData& data, SimulationData& simData);
	void loadSourceData(MainData& data, bool sourceIsArray, SourceData& srcData);
	void loadSimulationData(const MainData& data, const std::string& irMethod, SimulationData& simData);
	void prepareExcitation(TFloat dt, const SimulationData& simData, std::vector<TFloat>& tExc,
				std::vector<TFloat>& dvdt, std::vector<TFloat>& tDvdt);

	void execImpulseResponse(bool sourceIsArray); // calculate p/(c*density)
	void execTransientRadiationPattern(bool sourceIsArray);
	void execTransientAcousticField(bool sourceIsArray);
	void execTransientPropagation(bool sourceIsArray);

	Project& project_;
};

template<typename TFloat>
SimRectangularSourceMethod<TFloat>::SimRectangularSourceMethod(Project& project)
		: project_(project)
{
}

template<typename TFloat>
void
SimRectangularSourceMethod<TFloat>::loadData(const ParameterMap& taskPM, MainData& data, SimulationData& simData)
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
SimRectangularSourceMethod<TFloat>::loadSourceData(MainData& data, bool sourceIsArray, SourceData& srcData)
{
	WavefrontObjFileWriter<TFloat> fw((project_.expDirectory() + "/source_geometry.obj").c_str());

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

		const TFloat hw = 0.5 * srcData.sourceWidth;
		const TFloat hh = 0.5 * srcData.sourceHeight;
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

		const TFloat hw = 0.5 * srcData.sourceWidth;
		const TFloat hh = 0.5 * srcData.sourceHeight;
		fw.addPoint(-hw,  hh, 0.0);
		fw.addPoint( hw,  hh, 0.0);
		fw.addPoint( hw, -hh, 0.0);
		fw.addPoint(-hw, -hh, 0.0);
		fw.addQuad(-4, -3, -2, -1);
		fw.write();
	}
}

template<typename TFloat>
void
SimRectangularSourceMethod<TFloat>::loadSimulationData(const MainData& data, const std::string& irMethod, SimulationData& simData)
{
	const ParamMapPtr simPM = project_.getSubParamMap("simulation_config_file");
	simData.samplingFreq = simPM->value<TFloat>("sampling_frequency_factor", 0.0, 10000.0) * data.nyquistRate;
	simPM->getValue(simData.excitationType, "excitation_type");
	simPM->getValue(simData.excNumPeriods , "excitation_num_periods", 0.0, 100.0);
	simData.irMethod = irMethod;

	Waveform::get(simData.excitationType, data.centerFreq, simData.samplingFreq, simData.excNumPeriods, simData.exc);

	if (simData.irMethod == "numeric" ||
			simData.irMethod == "numeric_cuda"     ||
			simData.irMethod == "numeric_cuda_cpu" ||
			simData.irMethod == "numeric_ocl"      ||
			simData.irMethod == "numeric_ocl_cpu") {
		simPM->getValue(simData.discretFactor, "sub_elem_size_factor", 0.0, 1.0e3);
	} else if (simData.irMethod == "analytic") {
		simPM->getValue(simData.discretFactor, "min_edge_divisor"    , 0.0, 1.0e6);
	} else {
		THROW_EXCEPTION(InvalidParameterException, "Invalid impulse response method: " << simData.irMethod << '.');
	}
}

template<typename TFloat>
void
SimRectangularSourceMethod<TFloat>::prepareExcitation(TFloat dt,
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
SimRectangularSourceMethod<TFloat>::execTransientRadiationPattern(bool sourceIsArray)
{
	const ParameterMap& taskPM = project_.taskParamMap();
	MainData mainData;
	SimulationData simData;
	loadData(taskPM, mainData, simData);
	SourceData srcData;
	loadSourceData(mainData, sourceIsArray, srcData);

	project_.createDirectory(mainData.outputDir, false);

	const ParamMapPtr radPM = project_.getSubParamMap("rad_config_file");
	const auto distance   = radPM->value<TFloat>("distance", 0.0, 100.0);
	const auto thetaYStep = radPM->value<TFloat>("theta_y_step", 0.0, 10.0);
	const auto thetaYMin  = sourceIsArray ? radPM->value<TFloat>("theta_y_min" , -90.0, 90.0) : TFloat(0);
	const auto thetaYMax  = radPM->value<TFloat>("theta_y_max" , thetaYMin + 0.1, 90.0);
	const auto thetaXStep = radPM->value<TFloat>("theta_x_step", 0.0, 10.0);
	const auto thetaXMin  = sourceIsArray ? radPM->value<TFloat>("theta_x_min" , -90.0, 90.0) : TFloat(0);
	const auto thetaXMax  = radPM->value<TFloat>("theta_x_max" , thetaXMin + 0.1, 90.0);

	const TFloat dt = 1.0 / simData.samplingFreq;
	std::vector<TFloat> tExc, dvdt, tDvdt;
	prepareExcitation(dt, simData, tExc, dvdt, tDvdt);

	std::vector<TFloat> thetaXList;
	Util::fillSequenceFromStartToEndWithSize(thetaXList, thetaXMin, thetaXMax, std::ceil((thetaXMax - thetaXMin) / thetaXStep) + 1);
	std::vector<TFloat> thetaYList;
	Util::fillSequenceFromStartToEndWithSize(thetaYList, thetaYMin, thetaYMax, std::ceil((thetaYMax - thetaYMin) / thetaYStep) + 1);

	Matrix<XYZValue<TFloat>> gridData{thetaXList.size(), thetaYList.size()};
	Matrix<XYZ<TFloat>> inputData{thetaXList.size(), thetaYList.size()};

	for (unsigned int ix = 0, xSize = thetaXList.size(); ix < xSize; ++ix) {
		const TFloat tX = Util::degreeToRadian(thetaXList[ix]);
		const TFloat sinTX = std::sin(tX);
		const TFloat cosTX = std::cos(tX);
		for (unsigned int iy = 0, ySize = thetaYList.size(); iy < ySize; ++iy) {
			const TFloat tY = Util::degreeToRadian(thetaYList[iy]);
			const TFloat x = distance * std::sin(tY);
			const TFloat rx = distance * std::cos(tY);
			const TFloat y = rx * sinTX;
			const TFloat z = rx * cosTX;
			XYZValue<TFloat>& gd = gridData(ix, iy);
			gd.x = thetaYList[iy];
			gd.y = 0.0;
			gd.z = thetaXList[ix];
			gd.value = 0.0;
			XYZ<TFloat>& id = inputData(ix, iy);
			id.x = x;
			id.y = y;
			id.z = z;
		}
	}

	if (simData.irMethod == "numeric") {
		const TFloat subElemSize = mainData.propagationSpeed / (mainData.nyquistRate * simData.discretFactor);
		if (sourceIsArray) {
			SimTransientRadiationPattern<
				TFloat,
				NumericArrayOfRectangularSourcesImpulseResponse<TFloat>>::getArrayOfRectangularSourcesRadiationPatternDirect(
						simData.samplingFreq, mainData.propagationSpeed,
						srcData.sourceWidth, srcData.sourceHeight,
						subElemSize,
						dvdt, srcData.elemPos, srcData.focusDelay, inputData, gridData);
		} else {
			SimTransientRadiationPattern<
				TFloat,
				NumericRectangularSourceImpulseResponse<TFloat>>::getRectangularSourceRadiationPattern(
						simData.samplingFreq, mainData.propagationSpeed,
						srcData.sourceWidth, srcData.sourceHeight,
						subElemSize,
						dvdt, inputData, gridData);
		}
#ifdef USE_CUDA
	} else if (simData.irMethod == "numeric_cuda") {
		if constexpr (std::is_same<TFloat, float>::value) {
			const TFloat subElemSize = mainData.propagationSpeed / (mainData.nyquistRate * simData.discretFactor);
			if (sourceIsArray) {
				SimTransientRadiationPattern<
					TFloat,
					NumericArrayOfRectangularSourcesCUDAImpulseResponse>::getArrayOfRectangularSourcesRadiationPatternDirectSingleThread(
							simData.samplingFreq, mainData.propagationSpeed,
							srcData.sourceWidth, srcData.sourceHeight,
							subElemSize,
							dvdt, srcData.elemPos, srcData.focusDelay, inputData, gridData);
			} else {
				SimTransientRadiationPattern<
					TFloat,
					NumericRectangularSourceCUDAImpulseResponse>::getRectangularSourceRadiationPatternSingleThread(
							simData.samplingFreq, mainData.propagationSpeed,
							srcData.sourceWidth, srcData.sourceHeight,
							subElemSize,
							dvdt, inputData, gridData);
			}
		} else {
			THROW_EXCEPTION(InvalidValueException, "Invalid float type.");
		}
#endif
#ifdef USE_OPENCL
	} else if (simData.irMethod == "numeric_ocl") {
		const TFloat subElemSize = mainData.propagationSpeed / (mainData.nyquistRate * simData.discretFactor);
		if (sourceIsArray) {
			SimTransientRadiationPattern<
				TFloat,
				NumericArrayOfRectangularSourcesOCLImpulseResponse<TFloat>>::getArrayOfRectangularSourcesRadiationPatternDirectSingleThread(
						simData.samplingFreq, mainData.propagationSpeed,
						srcData.sourceWidth, srcData.sourceHeight,
						subElemSize,
						dvdt, srcData.elemPos, srcData.focusDelay, inputData, gridData);
		} else {
			SimTransientRadiationPattern<
				TFloat,
				NumericRectangularSourceOCLImpulseResponse<TFloat>>::getRectangularSourceRadiationPatternSingleThread(
						simData.samplingFreq, mainData.propagationSpeed,
						srcData.sourceWidth, srcData.sourceHeight,
						subElemSize,
						dvdt, inputData, gridData);
		}
#endif
	} else if (simData.irMethod == "analytic") {
		const TFloat minEdgeDivisor = simData.discretFactor;
		if (sourceIsArray) {
			SimTransientRadiationPattern<
				TFloat,
				AnalyticRectangularSourceImpulseResponse<TFloat>>::getArrayOfRectangularSourcesRadiationPattern(
						simData.samplingFreq, mainData.propagationSpeed,
						srcData.sourceWidth, srcData.sourceHeight,
						minEdgeDivisor,
						dvdt, srcData.elemPos, srcData.focusDelay, inputData, gridData);
		} else {
			SimTransientRadiationPattern<
				TFloat,
				AnalyticRectangularSourceImpulseResponse<TFloat>>::getRectangularSourceRadiationPattern(
						simData.samplingFreq, mainData.propagationSpeed,
						srcData.sourceWidth, srcData.sourceHeight,
						minEdgeDivisor,
						dvdt, inputData, gridData);
		}
	} else {
		THROW_EXCEPTION(InvalidParameterException, "Invalid impulse response method: " << simData.irMethod << '.');
	}

	const TFloat maxAbsValue = Util::maxAbsoluteValueField<XYZValue<TFloat>, TFloat>(gridData);
	const TFloat k = 1.0 / maxAbsValue;
	for (auto it = gridData.begin(); it != gridData.end(); ++it) {
		it->value *= k;
	}

	std::vector<XYZ<float>> pointList = {{0.0, 0.0, 0.0}};

	project_.showFigure3D(1, "Pattern", &gridData, &pointList,
					true, Visualization::VALUE_RECTIFIED_LOG, Colormap::GRADIENT_VIRIDIS);

	std::size_t sectionTYIndex = 0;
	std::size_t sectionTXIndex = 0;
	if (sourceIsArray) {
		TFloat sectionTY = 0.0;
		TFloat sectionTX = 0.0;
		if (srcData.useFocus) {
			sectionTY = Util::radianToDegree(std::atan2(srcData.focusX, srcData.focusZ));
			sectionTX = Util::radianToDegree(std::atan2(srcData.focusY, srcData.focusZ));
		}
		sectionTYIndex = std::rint(((sectionTY - thetaYMin) / (thetaYMax - thetaYMin)) * (thetaYList.size() - 1));
		sectionTXIndex = std::rint(((sectionTX - thetaXMin) / (thetaXMax - thetaXMin)) * (thetaXList.size() - 1));
		LOG_DEBUG << "Section theta-x: " << thetaXList[sectionTXIndex] << " theta-y: " << thetaYList[sectionTYIndex];
	}

	std::vector<TFloat> patternTY(gridData.n2());
	auto tyRange = gridData.range2(sectionTXIndex);
	Util::copyUsingOperator(tyRange.begin(), tyRange.end(), patternTY.begin(), Util::CopyValueOp{});
	Util::linearToDecibels(patternTY, TFloat(-100.0));
	project_.showFigure2D(3, "Pattern theta-y", thetaYList, patternTY);

	std::vector<TFloat> patternTX(gridData.n1());
	auto txRange = gridData.range1(sectionTYIndex);
	Util::copyUsingOperator(txRange.begin(), txRange.end(), patternTX.begin(), Util::CopyValueOp{});
	Util::linearToDecibels(patternTX, TFloat(-100.0));
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

template<typename TFloat>
void
SimRectangularSourceMethod<TFloat>::execTransientAcousticField(bool sourceIsArray)
{
	const ParameterMap& taskPM = project_.taskParamMap();
	MainData mainData;
	SimulationData simData;
	loadData(taskPM, mainData, simData);
	SourceData srcData;
	loadSourceData(mainData, sourceIsArray, srcData);

	project_.createDirectory(mainData.outputDir, false);

	const TFloat dt = 1.0 / simData.samplingFreq;
	std::vector<TFloat> tExc, dvdt, tDvdt;
	prepareExcitation(dt, simData, tExc, dvdt, tDvdt);

	Matrix<XYZValue<TFloat>> gridData;

	const TFloat nyquistLambda = Util::wavelength(mainData.propagationSpeed, mainData.nyquistRate);
	ImageGrid<TFloat>::get(*project_.getSubParamMap("grid_config_file"), nyquistLambda, gridData);

	if (simData.irMethod == "numeric") {
		const TFloat subElemSize = mainData.propagationSpeed / (mainData.nyquistRate * simData.discretFactor);
		if (sourceIsArray) {
			SimTransientAcousticField<TFloat>::template getArrayOfRectangularSourcesAcousticFieldDirect<NumericArrayOfRectangularSourcesImpulseResponse<TFloat>>(
						simData.samplingFreq, mainData.propagationSpeed,
						srcData.sourceWidth, srcData.sourceHeight,
						subElemSize,
						dvdt, srcData.elemPos, srcData.focusDelay, gridData);
		} else {
			SimTransientAcousticField<TFloat>::template getRectangularSourceAcousticField<NumericRectangularSourceImpulseResponse<TFloat>>(
						simData.samplingFreq, mainData.propagationSpeed,
						srcData.sourceWidth, srcData.sourceHeight,
						subElemSize,
						dvdt, gridData);
		}
#ifdef USE_CUDA
	} else if (simData.irMethod == "numeric_cuda") {
		if constexpr (std::is_same<TFloat, float>::value) {
			const TFloat subElemSize = mainData.propagationSpeed / (mainData.nyquistRate * simData.discretFactor);
			if (sourceIsArray) {
				SimTransientAcousticField<TFloat>::template getArrayOfRectangularSourcesAcousticFieldDirectSingleThread<NumericArrayOfRectangularSourcesCUDAImpulseResponse>(
							simData.samplingFreq, mainData.propagationSpeed,
							srcData.sourceWidth, srcData.sourceHeight,
							subElemSize,
							dvdt, srcData.elemPos, srcData.focusDelay, gridData);
			} else {
				SimTransientAcousticField<TFloat>::template getRectangularSourceAcousticFieldSingleThread<NumericRectangularSourceCUDAImpulseResponse>(
							simData.samplingFreq, mainData.propagationSpeed,
							srcData.sourceWidth, srcData.sourceHeight,
							subElemSize,
							dvdt, gridData);
			}
		} else {
			THROW_EXCEPTION(InvalidValueException, "Invalid float type.");
		}
	} else if (simData.irMethod == "numeric_cuda_cpu") {
		if constexpr (std::is_same<TFloat, float>::value) {
			const TFloat subElemSize = mainData.propagationSpeed / (mainData.nyquistRate * simData.discretFactor);
			if (sourceIsArray) {
				SimTransientAcousticField<TFloat>::template getArrayOfRectangularSourcesAcousticFieldDirectSTMT<
										NumericArrayOfRectangularSourcesCUDAImpulseResponse,
										NumericArrayOfRectangularSourcesImpulseResponse<TFloat>>(
							simData.samplingFreq, mainData.propagationSpeed,
							srcData.sourceWidth, srcData.sourceHeight,
							subElemSize,
							dvdt, srcData.elemPos, srcData.focusDelay, gridData);
			} else {
				SimTransientAcousticField<TFloat>::template getRectangularSourceAcousticFieldSTMT<
										NumericRectangularSourceCUDAImpulseResponse,
										NumericRectangularSourceImpulseResponse<TFloat>>(
							simData.samplingFreq, mainData.propagationSpeed,
							srcData.sourceWidth, srcData.sourceHeight,
							subElemSize,
							dvdt, gridData);
			}
		} else {
			THROW_EXCEPTION(InvalidValueException, "Invalid float type.");
		}
#endif
#ifdef USE_OPENCL
	} else if (simData.irMethod == "numeric_ocl") {
		const TFloat subElemSize = mainData.propagationSpeed / (mainData.nyquistRate * simData.discretFactor);
		if (sourceIsArray) {
			SimTransientAcousticField<TFloat>::template getArrayOfRectangularSourcesAcousticFieldDirectSingleThread<NumericArrayOfRectangularSourcesOCLImpulseResponse<TFloat>>(
						simData.samplingFreq, mainData.propagationSpeed,
						srcData.sourceWidth, srcData.sourceHeight,
						subElemSize,
						dvdt, srcData.elemPos, srcData.focusDelay, gridData);
		} else {
			SimTransientAcousticField<TFloat>::template getRectangularSourceAcousticFieldSingleThread<NumericRectangularSourceOCLImpulseResponse<TFloat>>(
						simData.samplingFreq, mainData.propagationSpeed,
						srcData.sourceWidth, srcData.sourceHeight,
						subElemSize,
						dvdt, gridData);
		}
	} else if (simData.irMethod == "numeric_ocl_cpu") {
		const TFloat subElemSize = mainData.propagationSpeed / (mainData.nyquistRate * simData.discretFactor);
		if (sourceIsArray) {
			SimTransientAcousticField<TFloat>::template getArrayOfRectangularSourcesAcousticFieldDirectSTMT<
									NumericArrayOfRectangularSourcesOCLImpulseResponse<TFloat>,
									NumericArrayOfRectangularSourcesImpulseResponse<TFloat>>(
						simData.samplingFreq, mainData.propagationSpeed,
						srcData.sourceWidth, srcData.sourceHeight,
						subElemSize,
						dvdt, srcData.elemPos, srcData.focusDelay, gridData);
		} else {
			SimTransientAcousticField<TFloat>::template getRectangularSourceAcousticFieldSTMT<
									NumericRectangularSourceOCLImpulseResponse<TFloat>,
									NumericRectangularSourceImpulseResponse<TFloat>>(
						simData.samplingFreq, mainData.propagationSpeed,
						srcData.sourceWidth, srcData.sourceHeight,
						subElemSize,
						dvdt, gridData);
		}
#endif
	} else if (simData.irMethod == "analytic") {
		const TFloat minEdgeDivisor = simData.discretFactor;
		if (sourceIsArray) {
			SimTransientAcousticField<TFloat>::template getArrayOfRectangularSourcesAcousticField<AnalyticRectangularSourceImpulseResponse<TFloat>>(
						simData.samplingFreq, mainData.propagationSpeed,
						srcData.sourceWidth, srcData.sourceHeight,
						minEdgeDivisor,
						dvdt, srcData.elemPos, srcData.focusDelay, gridData);
		} else {
			SimTransientAcousticField<TFloat>::template getRectangularSourceAcousticField<AnalyticRectangularSourceImpulseResponse<TFloat>>(
						simData.samplingFreq, mainData.propagationSpeed,
						srcData.sourceWidth, srcData.sourceHeight,
						minEdgeDivisor,
						dvdt, gridData);
		}
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
SimRectangularSourceMethod<TFloat>::execTransientPropagation(bool sourceIsArray)
{
	Timer tProc;

	const ParameterMap& taskPM = project_.taskParamMap();
	MainData mainData;
	SimulationData simData;
	loadData(taskPM, mainData, simData);
	SourceData srcData;
	loadSourceData(mainData, sourceIsArray, srcData);

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
		const TFloat subElemSize = mainData.propagationSpeed / (mainData.nyquistRate * simData.discretFactor);
		if (sourceIsArray) {
			SimTransientPropagation<
				TFloat,
				NumericArrayOfRectangularSourcesImpulseResponse<TFloat>>::getArrayOfRectangularSourcesPropagationDirect(
						simData.samplingFreq, mainData.propagationSpeed,
						srcData.sourceWidth, srcData.sourceHeight,
						subElemSize,
						dvdt, srcData.elemPos, srcData.focusDelay, propagIndexList, gridData);
		} else {
			SimTransientPropagation<
				TFloat,
				NumericRectangularSourceImpulseResponse<TFloat>>::getRectangularSourcePropagation(
						simData.samplingFreq, mainData.propagationSpeed,
						srcData.sourceWidth, srcData.sourceHeight,
						subElemSize,
						dvdt, propagIndexList, gridData);
		}
#ifdef USE_CUDA
	} else if (simData.irMethod == "numeric_cuda") {
		if constexpr (std::is_same<TFloat, float>::value) {
			const TFloat subElemSize = mainData.propagationSpeed / (mainData.nyquistRate * simData.discretFactor);
			if (sourceIsArray) {
				SimTransientPropagation<
					TFloat,
					NumericArrayOfRectangularSourcesCUDAImpulseResponse>::getArrayOfRectangularSourcesPropagationDirectSingleThread(
							simData.samplingFreq, mainData.propagationSpeed,
							srcData.sourceWidth, srcData.sourceHeight,
							subElemSize,
							dvdt, srcData.elemPos, srcData.focusDelay, propagIndexList, gridData);
			} else {
				SimTransientPropagation<
					TFloat,
					NumericRectangularSourceCUDAImpulseResponse>::getRectangularSourcePropagationSingleThread(
							simData.samplingFreq, mainData.propagationSpeed,
							srcData.sourceWidth, srcData.sourceHeight,
							subElemSize,
							dvdt, propagIndexList, gridData);
			}
		} else {
			THROW_EXCEPTION(InvalidValueException, "Invalid float type.");
		}
#endif
#ifdef USE_OPENCL
	} else if (simData.irMethod == "numeric_ocl") {
		const TFloat subElemSize = mainData.propagationSpeed / (mainData.nyquistRate * simData.discretFactor);
		if (sourceIsArray) {
			SimTransientPropagation<
				TFloat,
				NumericArrayOfRectangularSourcesOCLImpulseResponse<TFloat>>::getArrayOfRectangularSourcesPropagationDirectSingleThread(
						simData.samplingFreq, mainData.propagationSpeed,
						srcData.sourceWidth, srcData.sourceHeight,
						subElemSize,
						dvdt, srcData.elemPos, srcData.focusDelay, propagIndexList, gridData);
		} else {
			SimTransientPropagation<
				TFloat,
				NumericRectangularSourceOCLImpulseResponse<TFloat>>::getRectangularSourcePropagationSingleThread(
						simData.samplingFreq, mainData.propagationSpeed,
						srcData.sourceWidth, srcData.sourceHeight,
						subElemSize,
						dvdt, propagIndexList, gridData);
		}
#endif
	} else if (simData.irMethod == "analytic") {
		const TFloat minEdgeDivisor = simData.discretFactor;
		if (sourceIsArray) {
			SimTransientPropagation<
				TFloat,
				AnalyticRectangularSourceImpulseResponse<TFloat>>::getArrayOfRectangularSourcesPropagation(
						simData.samplingFreq, mainData.propagationSpeed,
						srcData.sourceWidth, srcData.sourceHeight,
						minEdgeDivisor,
						dvdt, srcData.elemPos, srcData.focusDelay, propagIndexList, gridData);
		} else {
			SimTransientPropagation<
				TFloat,
				AnalyticRectangularSourceImpulseResponse<TFloat>>::getRectangularSourcePropagation(
						simData.samplingFreq, mainData.propagationSpeed,
						srcData.sourceWidth, srcData.sourceHeight,
						minEdgeDivisor,
						dvdt, propagIndexList, gridData);
		}
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
SimRectangularSourceMethod<TFloat>::execImpulseResponse(bool sourceIsArray)
{
	const ParameterMap& taskPM = project_.taskParamMap();
	MainData mainData;
	SimulationData simData;
	loadData(taskPM, mainData, simData);
	SourceData srcData;
	loadSourceData(mainData, sourceIsArray, srcData);

	project_.createDirectory(mainData.outputDir, false);

	const ParamMapPtr irPM = project_.getSubParamMap("ir_config_file");
	const auto pointX = irPM->value<TFloat>("point_x", sourceIsArray ? -10000.0 : 0.0, 10000.0);
	const auto pointY = irPM->value<TFloat>("point_y", sourceIsArray ? -10000.0 : 0.0, 10000.0);
	const auto pointZ = irPM->value<TFloat>("point_z", 0.0, 10000.0);

	const TFloat dt = 1.0 / simData.samplingFreq;
	std::vector<TFloat> tExc, dvdt, tDvdt;
	prepareExcitation(dt, simData, tExc, dvdt, tDvdt);

	std::size_t hOffset;
	std::vector<TFloat> h;
	if (simData.irMethod == "numeric") {
		const TFloat subElemSize = mainData.propagationSpeed / (mainData.nyquistRate * simData.discretFactor);
		if (sourceIsArray) {
			auto impResp = std::make_unique<NumericArrayOfRectangularSourcesImpulseResponse<TFloat>>(
						simData.samplingFreq, mainData.propagationSpeed, srcData.sourceWidth, srcData.sourceHeight,
						subElemSize,
						srcData.elemPos, srcData.focusDelay);
			impResp->getImpulseResponse(pointX, pointY, pointZ, hOffset, h);
		} else {
			auto impResp = std::make_unique<NumericRectangularSourceImpulseResponse<TFloat>>(
						simData.samplingFreq, mainData.propagationSpeed, srcData.sourceWidth, srcData.sourceHeight,
						subElemSize);
			impResp->getImpulseResponse(pointX, pointY, pointZ, hOffset, h);
		}
#ifdef USE_CUDA
	} else if (simData.irMethod == "numeric_cuda") {
		if constexpr (std::is_same<TFloat, float>::value) {
			const TFloat subElemSize = mainData.propagationSpeed / (mainData.nyquistRate * simData.discretFactor);
			if (sourceIsArray) {
				auto impResp = std::make_unique<NumericArrayOfRectangularSourcesCUDAImpulseResponse>(
							simData.samplingFreq, mainData.propagationSpeed, srcData.sourceWidth, srcData.sourceHeight,
							subElemSize,
							srcData.elemPos, srcData.focusDelay);
				impResp->getImpulseResponse(pointX, pointY, pointZ, hOffset, h);
			} else {
				auto impResp = std::make_unique<NumericRectangularSourceCUDAImpulseResponse>(
							simData.samplingFreq, mainData.propagationSpeed, srcData.sourceWidth, srcData.sourceHeight,
							subElemSize);
				impResp->getImpulseResponse(pointX, pointY, pointZ, hOffset, h);
			}
		} else {
			THROW_EXCEPTION(InvalidValueException, "Invalid float type.");
		}
#endif
#ifdef USE_OPENCL
	} else if (simData.irMethod == "numeric_ocl") {
		const TFloat subElemSize = mainData.propagationSpeed / (mainData.nyquistRate * simData.discretFactor);
		if (sourceIsArray) {
			auto impResp = std::make_unique<NumericArrayOfRectangularSourcesOCLImpulseResponse<TFloat>>(
						simData.samplingFreq, mainData.propagationSpeed, srcData.sourceWidth, srcData.sourceHeight,
						subElemSize,
						srcData.elemPos, srcData.focusDelay);
			impResp->getImpulseResponse(pointX, pointY, pointZ, hOffset, h);
		} else {
			auto impResp = std::make_unique<NumericRectangularSourceOCLImpulseResponse<TFloat>>(
						simData.samplingFreq, mainData.propagationSpeed, srcData.sourceWidth, srcData.sourceHeight,
						subElemSize);
			impResp->getImpulseResponse(pointX, pointY, pointZ, hOffset, h);
		}
#endif
	} else if (simData.irMethod == "analytic") {
		const TFloat minEdgeDivisor = simData.discretFactor;
		if (sourceIsArray) {
			auto impResp = std::make_unique<ArrayOfRectangularSourcesImpulseResponse<TFloat, AnalyticRectangularSourceImpulseResponse<TFloat>>>(
						simData.samplingFreq, mainData.propagationSpeed, srcData.sourceWidth, srcData.sourceHeight,
						minEdgeDivisor,
						srcData.elemPos, srcData.focusDelay);
			impResp->getImpulseResponse(pointX, pointY, pointZ, hOffset, h);
		} else {
			auto impResp = std::make_unique<AnalyticRectangularSourceImpulseResponse<TFloat>>(
						simData.samplingFreq, mainData.propagationSpeed, srcData.sourceWidth, srcData.sourceHeight,
						minEdgeDivisor);
			impResp->getImpulseResponse(pointX, pointY, pointZ, hOffset, h);
		}
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
SimRectangularSourceMethod<TFloat>::execute()
{
	switch (project_.method()) {
	case MethodEnum::sim_acoustic_field_rectangular_source_transient:
	case MethodEnum::sim_acoustic_field_rectangular_source_transient_sp:
		execTransientAcousticField(false);
		break;
	case MethodEnum::sim_acoustic_field_array_of_rectangular_sources_transient:
	case MethodEnum::sim_acoustic_field_array_of_rectangular_sources_transient_sp:
		execTransientAcousticField(true);
		break;
	case MethodEnum::sim_impulse_response_rectangular_source:
	case MethodEnum::sim_impulse_response_rectangular_source_sp:
		execImpulseResponse(false);
		break;
	case MethodEnum::sim_impulse_response_array_of_rectangular_sources:
	case MethodEnum::sim_impulse_response_array_of_rectangular_sources_sp:
		execImpulseResponse(true);
		break;
	case MethodEnum::sim_propagation_rectangular_source_transient:
	case MethodEnum::sim_propagation_rectangular_source_transient_sp:
		execTransientPropagation(false);
		break;
	case MethodEnum::sim_propagation_array_of_rectangular_sources_transient:
	case MethodEnum::sim_propagation_array_of_rectangular_sources_transient_sp:
		execTransientPropagation(true);
		break;
	case MethodEnum::sim_radiation_pattern_rectangular_source_transient:
	case MethodEnum::sim_radiation_pattern_rectangular_source_transient_sp:
		execTransientRadiationPattern(false);
		break;
	case MethodEnum::sim_radiation_pattern_array_of_rectangular_sources_transient:
	case MethodEnum::sim_radiation_pattern_array_of_rectangular_sources_transient_sp:
		execTransientRadiationPattern(true);
		break;
	default:
		THROW_EXCEPTION(InvalidParameterException, "Invalid method: " << static_cast<int>(project_.method()) << '.');
	}
}

} // namespace Lab

#endif // SIMRECTANGULARSOURCEMETHOD_H
