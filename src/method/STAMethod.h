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

#ifndef STAMETHOD_H_
#define STAMETHOD_H_

#include <algorithm>
#include <iomanip>
#include <memory>

#include "CoherenceFactor.h"
#include "DefaultSTAProcessor.h"
#include "HilbertEnvelope.h"
#include "ImageGrid.h"
#include "Log.h"
#include "Method.h"
#include "NetworkSTAAcquisition.h"
#include "ParallelHilbertEnvelope.h"
#include "Project.h"
#include "SavedSTAAcquisition.h"
#include "SimpleSTAProcessor.h"
#include "Simulated3DSTAAcquisition.h"
#include "SimulatedSTAAcquisition.h"
#include "STAAcquisition.h"
#include "STAConfiguration.h"
#include "STAProcessor.h"
#include "Timer.h"
#include "VectorialSTAProcessor.h"
#include "Util.h"
#include "XZValue.h"
#include "XZValueFactor.h"



namespace Lab {

template<typename FloatType>
class STAMethod : public Method {
public:
	STAMethod(Project& project);
	virtual ~STAMethod();

	virtual void execute();

private:
	STAMethod(const STAMethod&);
	STAMethod& operator=(const STAMethod&);

	void saveSignals(const STAConfiguration<FloatType>& config, STAAcquisition<FloatType>& acq, unsigned int baseElement, const std::string& dataDir);

	Project& project_;
};



template<typename FloatType>
STAMethod<FloatType>::STAMethod(Project& project)
		: project_(project)
{
}

template<typename FloatType>
STAMethod<FloatType>::~STAMethod()
{
}

template<typename FloatType>
void
STAMethod<FloatType>::execute()
{
	ConstParameterMapPtr taskPM = project_.taskParameterMap();

	const STAConfiguration<FloatType> config(project_.loadChildParameterMap(taskPM, "sta_config_file"));
	const unsigned int baseElement = taskPM->value<unsigned int>("base_element", 0, config.numElementsMux - config.numElements);

	std::unique_ptr<STAAcquisition<FloatType>> acquisition;

	switch (project_.method()) {
	case MethodType::sta_simple_simulated: // falls through
	case MethodType::sta_simulated:
		acquisition = std::make_unique<SimulatedSTAAcquisition<FloatType>>(project_, config);
		break;
	case MethodType::sta_dp_network:   // falls through
	case MethodType::sta_save_signals:
		acquisition = std::make_unique<NetworkSTAAcquisition<FloatType>>(project_, config);
		break;
	case MethodType::sta_simple_saved:       // falls through
	case MethodType::sta_dp_saved:           // falls through
	case MethodType::sta_vectorial_dp_saved: // falls through
	case MethodType::sta_sp_saved:           // falls through
	case MethodType::sta_vectorial_sp_saved:
		acquisition = std::make_unique<SavedSTAAcquisition<FloatType>>(project_, config.numElements);
		break;
	case MethodType::sta_simulated_3d:           // falls through
	case MethodType::sta_vectorial_simulated_3d:
		acquisition = std::make_unique<Simulated3DSTAAcquisition<FloatType>>(project_, config);
		break;
	default:
		THROW_EXCEPTION(InvalidParameterException, "Invalid method: " << static_cast<int>(project_.method()) << '.');
	}

	if (project_.method() == MethodType::sta_save_signals) {
		std::string dataDir = taskPM->value<std::string>("data_dir");
		saveSignals(config, *acquisition, baseElement, dataDir);
		return;
	}

	Matrix2<XZValueFactor<FloatType>> gridData;

	const FloatType lambda = config.propagationSpeed / config.centerFrequency;
	ImageGrid<FloatType>::get(project_.loadChildParameterMap(taskPM, "grid_config_file"), lambda, gridData);

	bool coherenceFactorEnabled = false;
	bool vectorialProcessingWithEnvelope = false;

	const FloatType peakOffset  = taskPM->value<FloatType>(  "peak_offset", 0.0, 50.0);
	const std::string outputDir = taskPM->value<std::string>("output_dir");

	switch (project_.method()) {
	case MethodType::sta_simple_simulated: // falls through
	case MethodType::sta_simple_saved:
		{
			auto processor = std::make_unique<SimpleSTAProcessor<FloatType>>(config, *acquisition, peakOffset);

			Timer tProc;
			processor->process(baseElement, gridData);
			LOG_DEBUG << ">>> Acquisition + processing time: " << tProc.getTime();
		}
		break;
	case MethodType::sta_vectorial_dp_saved:     // falls through
	case MethodType::sta_vectorial_sp_saved:     // falls through
	case MethodType::sta_vectorial_simulated_3d:
		{
			vectorialProcessingWithEnvelope     = taskPM->value<bool>(        "calculate_envelope_in_processing");
			const unsigned int upsamplingFactor = taskPM->value<unsigned int>("upsampling_factor", 1, 128);
			AnalyticSignalCoherenceFactorProcessor<FloatType> coherenceFactor(project_.loadChildParameterMap(taskPM, "coherence_factor_config_file"));
			coherenceFactorEnabled = coherenceFactor.enabled();
			auto processor = std::make_unique<VectorialSTAProcessor<FloatType>>(config, *acquisition, upsamplingFactor, coherenceFactor, peakOffset, vectorialProcessingWithEnvelope);

			Timer tProc;
			processor->process(baseElement, gridData);
			LOG_DEBUG << ">>> Acquisition + processing time: " << tProc.getTime();
		}
		break;
	default:
		{
			CoherenceFactorProcessor<FloatType> coherenceFactor(project_.loadChildParameterMap(taskPM, "coherence_factor_config_file"));
			coherenceFactorEnabled = coherenceFactor.enabled();
			auto processor = std::make_unique<DefaultSTAProcessor<FloatType>>(config, *acquisition, coherenceFactor, peakOffset);

			Timer tProc;
			processor->process(baseElement, gridData);
			LOG_DEBUG << ">>> Acquisition + processing time: " << tProc.getTime();
		}
	}

	project_.saveImageToHDF5(gridData, outputDir);

	Figure::Visualization visual;
	if (vectorialProcessingWithEnvelope) {
		visual = Figure::VISUALIZATION_RECTIFIED_LOG;
	} else {
		visual = Figure::VISUALIZATION_ENVELOPE_LOG;
	}

	std::vector<XZ<float>> pointList = {{0.0, 0.0}};

	Project::GridDataType projGridData;
	Util::copyXZValue(gridData, projGridData);
	project_.showFigure3D(1, "Raw image", &projGridData, &pointList,
				true, visual, Figure::COLORMAP_VIRIDIS);

	if (coherenceFactorEnabled) {
		LOG_DEBUG << "Saving the image factors...";
		project_.saveHDF5(gridData, outputDir + "/image_factor", "image", Util::CopyFactorOp());

		if (!vectorialProcessingWithEnvelope) {
			ParallelHilbertEnvelope<FloatType>::calculateDim2(gridData);
		}

		// Applies the coherence factor method.
		for (auto iter = gridData.begin(); iter != gridData.end(); ++iter) {
			iter->value *= iter->factor;
			iter->factor = 1.0;
		}

		LOG_DEBUG << "Saving the CF image...";
		project_.saveHDF5(gridData, outputDir + "/image_cf", "image", Util::CopyValueOp());

		Util::copyXZValue(gridData, projGridData);
		project_.showFigure3D(2, "Coherence factor image", &projGridData, &pointList,
					true, Figure::VISUALIZATION_RECTIFIED_LOG, Figure::COLORMAP_VIRIDIS);
	}
}

template<typename FloatType>
void
STAMethod<FloatType>::saveSignals(const STAConfiguration<FloatType>& config, STAAcquisition<FloatType>& acq, unsigned int baseElement, const std::string& dataDir)
{
	typename STAAcquisition<FloatType>::AcquisitionDataType acqData;

	for (unsigned int txElem = 0; txElem < config.numElements; ++txElem) {
		acq.execute(baseElement, txElem, acqData);

		std::ostringstream filePath;
		filePath << dataDir << std::setfill('0') << "/signal-base" << std::setw(4) << baseElement << "-tx" << std::setw(4) << txElem;
		std::string filePathStr = filePath.str();
		LOG_DEBUG << "Saving " << filePathStr << "...";
		project_.saveHDF5(acqData, filePathStr, "signal");
	}
}

} // namespace Lab

#endif /* STAMETHOD_H_ */
