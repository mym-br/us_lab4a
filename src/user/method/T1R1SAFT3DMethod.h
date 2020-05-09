/***************************************************************************
 *  Copyright 2014, 2017, 2018 Marcelo Y. Matuda                           *
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

#ifndef T1R1SAFT3DMETHOD_H
#define T1R1SAFT3DMETHOD_H

#include <algorithm> /* max */
#include <cstddef> /* std::size_t */
#include <memory>
#include <string>
#include <vector>

#include "CoherenceFactor.h"
#include "Colormap.h"
#include "Exception.h"
#include "ImageGrid.h"
#include "Log.h"
#include "Matrix.h"
#include "Method.h"
#include "Project.h"
#include "Simulated3DT1R1SAFTAcquisition.h"
#include "STAAcquisition.h"
#include "SA3DConfiguration.h"
#include "Timer.h"
#include "Util.h"
#include "Vectorial3DT1R1SAFTProcessor.h"
#include "Visualization.h"
#include "XYZ.h"
#include "XYZValueFactor.h"



namespace Lab {

template<typename TFloat>
class T1R1SAFT3DMethod : public Method {
public:
	T1R1SAFT3DMethod(Project& project);
	virtual ~T1R1SAFT3DMethod() = default;

	virtual void execute();
private:
	T1R1SAFT3DMethod(const T1R1SAFT3DMethod&) = delete;
	T1R1SAFT3DMethod& operator=(const T1R1SAFT3DMethod&) = delete;
	T1R1SAFT3DMethod(T1R1SAFT3DMethod&&) = delete;
	T1R1SAFT3DMethod& operator=(T1R1SAFT3DMethod&&) = delete;

	void process(TFloat valueScale, ArrayProcessor<XYZValueFactor<TFloat>>& processor, unsigned int baseElement, const std::string& outputDir);
	void useCoherenceFactor(TFloat valueScale, const std::string& outputDir);

	Project& project_;
	Matrix<XYZValueFactor<TFloat>> gridData_;
	std::vector<XYZ<float>> pointList_;
	Visualization::Value visual_;
};



template<typename TFloat>
T1R1SAFT3DMethod<TFloat>::T1R1SAFT3DMethod(Project& project)
		: project_(project)
		, pointList_{{0.0, 0.0, 0.0}}
		, visual_(Visualization::Value::RECTIFIED_LOG)
{
}

template<typename TFloat>
void
T1R1SAFT3DMethod<TFloat>::useCoherenceFactor(TFloat valueScale, const std::string& outputDir)
{
	project_.saveFactorToHDF5(gridData_, outputDir, "image_factor", "factor");

	Util::applyFactorToValue(gridData_.begin(), gridData_.end());

	project_.saveImageToHDF5(gridData_, outputDir, "image_cf", "cf");

	project_.showFigure3D(2, "Coherence factor image", &gridData_, &pointList_,
				true, Visualization::Value::RECTIFIED_LOG, Colormap::Id::GRADIENT_VIRIDIS, valueScale);
}

template<typename TFloat>
void
T1R1SAFT3DMethod<TFloat>::process(TFloat valueScale, ArrayProcessor<XYZValueFactor<TFloat>>& processor, unsigned int baseElement, const std::string& outputDir)
{
	Timer tProc;

	processor.prepare(baseElement);
	processor.process(gridData_);

	project_.saveImageToHDF5(gridData_, outputDir);
	project_.saveXYZToHDF5(gridData_, outputDir);

	project_.showFigure3D(1, "Raw image", &gridData_, &pointList_,
				true, visual_, Colormap::Id::GRADIENT_VIRIDIS, valueScale);

	LOG_DEBUG << ">>> Acquisition + processing time: " << tProc.getTime();
}

template<typename TFloat>
void
T1R1SAFT3DMethod<TFloat>::execute()
{
	const ParameterMap& taskPM = project_.taskParamMap();
	const ParamMapPtr saPM    = project_.getSubParamMap("sa_config_file");
	const ParamMapPtr arrayPM = project_.getSubParamMap("array_config_file");
	const SA3DConfiguration<TFloat> config(*saPM, *arrayPM);
	if (config.txElemPos.size() != config.rxElemPos.size()) {
		THROW_EXCEPTION(InvalidParameterException, "The number of receive elements is not equal to the number of transmit elements.");
	}
	if (config.activeTxElem.size() != config.activeRxElem.size()) {
		THROW_EXCEPTION(InvalidParameterException, "The number of active receive elements is not equal to the number of active transmit elements.");
	}

	const auto baseElement = saPM->value<unsigned int>("base_element", 0, config.numElementsMux - 1U);

	std::unique_ptr<STAAcquisition<TFloat>> acquisition;

	switch (project_.method()) {
	case MethodEnum::t1r1saft_3d_simulated_save_signals:
	case MethodEnum::t1r1saft_3d_simulated_seq_y_save_signals:
	case MethodEnum::t1r1saft_3d_vectorial_simulated:
		acquisition = std::make_unique<Simulated3DT1R1SAFTAcquisition<TFloat>>(project_, config);
		break;
	default:
		THROW_EXCEPTION(InvalidParameterException, "Invalid method: " << static_cast<int>(project_.method()) << '.');
	}

	if (project_.method() == MethodEnum::t1r1saft_3d_simulated_save_signals) {
		const auto dataDir = taskPM.value<std::string>("data_dir");
		typename STAAcquisition<TFloat>::AcquisitionDataType acqData;
		acquisition->prepare(baseElement);
		for (unsigned int txElem : config.activeTxElem) {
			acquisition->execute(txElem, acqData);
			project_.saveTxElemSignalsToHDF5(acqData, dataDir, 0, baseElement, txElem);
		}
		return;
	} else if (project_.method() == MethodEnum::t1r1saft_3d_simulated_seq_y_save_signals) {
		const auto dataDir = taskPM.value<std::string>("data_dir");
		const ParamMapPtr seqYPM = project_.getSubParamMap("seq_y_config_file");
		const auto yStep = seqYPM->value<TFloat>("y_step",          0.0,   100.0);
		const auto minY  = seqYPM->value<TFloat>("min_y" ,     -10000.0, 10000.0);
		const auto maxY  = seqYPM->value<TFloat>("max_y" , minY + yStep, 10000.0);
		project_.createDirectory(dataDir, true);
		typename STAAcquisition<TFloat>::AcquisitionDataType acqData;
		std::vector<TFloat> yList;
		Util::fillSequenceFromStartWithStep(yList, minY, maxY, yStep);
		auto& simAcq = dynamic_cast<Simulated3DT1R1SAFTAcquisition<TFloat>&>(*acquisition);
		acquisition->prepare(baseElement);
		for (std::size_t i = 0, end = yList.size(); i < end; ++i) {
			const TFloat y = yList[i];
			simAcq.modifyReflectorsOffset(0.0, -y);
			for (unsigned int txElem : config.activeTxElem) {
				acquisition->execute(txElem, acqData);
				project_.saveTxElemSignalsToHDF5(acqData, dataDir, i, baseElement, txElem);
			}
		}
		return;
	}

	const auto outputDir = taskPM.value<std::string>("output_dir");
	project_.createDirectory(outputDir, false);

	const TFloat nyquistLambda = Util::nyquistLambda(config.propagationSpeed, config.maxFrequency);
	ImageGrid<TFloat>::get(*project_.getSubParamMap("grid_config_file"), nyquistLambda, gridData_);

	if (project_.method() == MethodEnum::t1r1saft_3d_vectorial_simulated) {
		const ParamMapPtr imagPM = project_.getSubParamMap("imag_config_file");
		const auto peakOffset       = imagPM->value<TFloat>(      "peak_offset", 0.0, 50.0);
		const auto upsamplingFactor = imagPM->value<unsigned int>("upsampling_factor", 1, 128);
		AnalyticSignalCoherenceFactorProcessor<TFloat> coherenceFactor(*project_.getSubParamMap("coherence_factor_config_file"));
		auto processor = std::make_unique<Vectorial3DT1R1SAFTProcessor<TFloat>>(
							config, *acquisition, upsamplingFactor,
							coherenceFactor, peakOffset);
		process(config.valueScale, *processor, baseElement, outputDir);
		if (coherenceFactor.enabled()) {
			useCoherenceFactor(config.valueScale, outputDir);
		}
	} else {
		THROW_EXCEPTION(InvalidParameterException, "Invalid method: " << static_cast<int>(project_.method()) << '.');
	}
}

} // namespace Lab

#endif // T1R1SAFT3DMETHOD_H
