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
#include "Exception.h"
#include "global.h"
#include "ImageGrid.h"
#include "Log.h"
#include "Matrix.h"
#include "Method.h"
#include "Project.h"
#include "Simulated3DT1R1SAFTAcquisition.h"
#include "STAAcquisition.h"
#include "SA3DConfiguration.h"
#include "Timer.h"
#include "Vectorial3DT1R1SAFTProcessor.h"
#include "Util.h"
#include "XYZ.h"
#include "XYZValueFactor.h"



namespace Lab {

template<typename FloatType>
class T1R1SAFT3DMethod : public Method {
public:
	T1R1SAFT3DMethod(Project& project);
	virtual ~T1R1SAFT3DMethod();

	virtual void execute();

private:
	T1R1SAFT3DMethod(const T1R1SAFT3DMethod&) = delete;
	T1R1SAFT3DMethod& operator=(const T1R1SAFT3DMethod&) = delete;

	void process(FloatType valueScale, ArrayProcessor<FloatType>& processor, unsigned int baseElement, const std::string& outputDir);
	void useCoherenceFactor(FloatType valueScale, const std::string& outputDir);

	Project& project_;
	Matrix<XYZValueFactor<FloatType>> gridData_;
	Project::GridDataType projGridData_;
	std::vector<XYZ<float>> pointList_;
	Figure::Visualization visual_;
};



template<typename FloatType>
T1R1SAFT3DMethod<FloatType>::T1R1SAFT3DMethod(Project& project)
		: project_(project)
		, pointList_{{0.0, 0.0, 0.0}}
		, visual_(Figure::VISUALIZATION_RECTIFIED_LOG)
{
}

template<typename FloatType>
T1R1SAFT3DMethod<FloatType>::~T1R1SAFT3DMethod()
{
}

template<typename FloatType>
void
T1R1SAFT3DMethod<FloatType>::useCoherenceFactor(FloatType valueScale, const std::string& outputDir)
{
	LOG_DEBUG << "Saving the image factors...";
	project_.saveHDF5(gridData_, outputDir + "/image_factor", "factor", Util::CopyFactorOp());

	// Applies the coherence factor method.
	for (auto iter = gridData_.begin(); iter != gridData_.end(); ++iter) {
		iter->value *= iter->factor;
		iter->factor = 1.0;
	}

	LOG_DEBUG << "Saving the CF image...";
	project_.saveHDF5(gridData_, outputDir + "/image_cf", "cf", Util::CopyValueOp());

	Util::copyXYZValue(gridData_, projGridData_);
	project_.showFigure3D(2, "Coherence factor image", &projGridData_, &pointList_,
				true, Figure::VISUALIZATION_RECTIFIED_LOG, Figure::COLORMAP_VIRIDIS, valueScale);
}

template<typename FloatType>
void
T1R1SAFT3DMethod<FloatType>::process(FloatType valueScale, ArrayProcessor<FloatType>& processor, unsigned int baseElement, const std::string& outputDir)
{
	Timer tProc;

	processor.process(baseElement, gridData_);

	project_.saveImageToHDF5(gridData_, outputDir);

	Util::copyXYZValue(gridData_, projGridData_);
	project_.showFigure3D(1, "Raw image", &projGridData_, &pointList_,
				true, visual_, Figure::COLORMAP_VIRIDIS, valueScale);

	LOG_DEBUG << ">>> Acquisition + processing time: " << tProc.getTime();
}

template<typename FloatType>
void
T1R1SAFT3DMethod<FloatType>::execute()
{
	ConstParameterMapPtr taskPM = project_.taskParameterMap();
	ConstParameterMapPtr saPM    = project_.loadChildParameterMap(taskPM, "sa_config_file");
	ConstParameterMapPtr arrayPM = project_.loadChildParameterMap(taskPM, "array_config_file");
	const SA3DConfiguration<FloatType> config(saPM, arrayPM);
	if (config.txElemPos.size() != config.rxElemPos.size()) {
		THROW_EXCEPTION(InvalidParameterException, "The number of receive elements is not equal to the number of transmit elements.");
	}
	if (config.activeTxElem.size() != config.activeRxElem.size()) {
		THROW_EXCEPTION(InvalidParameterException, "The number of active receive elements is not equal to the number of active transmit elements.");
	}

	const unsigned int baseElement = taskPM->value<unsigned int>("base_element", 0, config.numElementsMux - 1U);

	std::unique_ptr<STAAcquisition<FloatType>> acquisition;

	switch (project_.method()) {
	case MethodType::t1r1saft_3d_simulated_save_signals:       // falls through
	case MethodType::t1r1saft_3d_simulated_seq_y_save_signals: // falls through
	case MethodType::t1r1saft_3d_vectorial_simulated:
		acquisition = std::make_unique<Simulated3DT1R1SAFTAcquisition<FloatType>>(project_, config);
		break;
	default:
		THROW_EXCEPTION(InvalidParameterException, "Invalid method: " << static_cast<int>(project_.method()) << '.');
	}

	if (project_.method() == MethodType::t1r1saft_3d_simulated_save_signals) {
		const std::string dataDir = taskPM->value<std::string>("data_dir");
		typename STAAcquisition<FloatType>::AcquisitionDataType acqData;
		for (unsigned int txElem : config.activeTxElem) {
			acquisition->execute(baseElement, txElem, acqData);
			project_.saveTxElemSignalsToHDF5(acqData, dataDir, 0, baseElement, txElem);
		}
		return;
	} else if (project_.method() == MethodType::t1r1saft_3d_simulated_seq_y_save_signals) {
		const std::string dataDir = taskPM->value<std::string>("data_dir");
		const FloatType yStep     = taskPM->value<FloatType>(  "y_step",          0.0,   100.0);
		const FloatType minY      = taskPM->value<FloatType>(  "min_y" ,     -10000.0, 10000.0);
		const FloatType maxY      = taskPM->value<FloatType>(  "max_y" , minY + yStep, 10000.0);
		project_.createDirectory(dataDir, true);
		typename STAAcquisition<FloatType>::AcquisitionDataType acqData;
		std::vector<FloatType> yList;
		Util::fillSequenceFromStartWithStep(yList, minY, maxY, yStep);
		auto& simAcq = dynamic_cast<Simulated3DT1R1SAFTAcquisition<FloatType>&>(*acquisition);
		for (std::size_t i = 0, end = yList.size(); i < end; ++i) {
			const FloatType y = yList[i];
			simAcq.modifyReflectorsOffset(0.0, -y);
			for (unsigned int txElem : config.activeTxElem) {
				acquisition->execute(baseElement, txElem, acqData);
				project_.saveTxElemSignalsToHDF5(acqData, dataDir, i, baseElement, txElem);
			}
		}
		return;
	}

	const FloatType peakOffset  = taskPM->value<FloatType>(  "peak_offset" , 0.0, 50.0);
	const std::string outputDir = taskPM->value<std::string>("output_dir");
	project_.createDirectory(outputDir, false);

	const FloatType nyquistRate = 2.0 * config.maxFrequency;
	const FloatType nyquistLambda = config.propagationSpeed / nyquistRate;
	ImageGrid<FloatType>::get(project_.loadChildParameterMap(taskPM, "grid_config_file"), nyquistLambda, gridData_);

	if (project_.method() == MethodType::t1r1saft_3d_vectorial_simulated) {
		const unsigned int upsamplingFactor = taskPM->value<unsigned int>("upsampling_factor", 1, 128);
		AnalyticSignalCoherenceFactorProcessor<FloatType> coherenceFactor(project_.loadChildParameterMap(taskPM, "coherence_factor_config_file"));
		auto processor = std::make_unique<Vectorial3DT1R1SAFTProcessor<FloatType>>(
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
