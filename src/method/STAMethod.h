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
#include <sstream>
#include <string>

#include "CoherenceFactor.h"
#include "DefaultSTAProcessor.h"
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
#include "Timer.h"
#include "VectorialSTAProcessor.h"
#include "Util.h"
#include "XZ.h"
#include "XZValueFactor.h"



namespace Lab {

template<typename FloatType>
class STAMethod : public Method {
public:
	STAMethod(Project& project);
	virtual ~STAMethod();

	virtual void execute();

private:
	STAMethod(const STAMethod&) = delete;
	STAMethod& operator=(const STAMethod&) = delete;

	void process(STAProcessor<FloatType>& processor, unsigned int baseElement, const std::string& outputDir);
	void useCoherenceFactor(bool calculateEnvelope, const std::string& outputDir);

	Project& project_;
	Matrix2<XZValueFactor<FloatType>> gridData_;
	Project::GridDataType projGridData_;
	std::vector<XZ<float>> pointList_;
	Figure::Visualization visual_;
	std::string outputDir_;
};



template<typename FloatType>
STAMethod<FloatType>::STAMethod(Project& project)
		: project_{project}
		, pointList_{{0.0, 0.0}}
		, visual_{Figure::VISUALIZATION_ENVELOPE_LOG}
{
}

template<typename FloatType>
STAMethod<FloatType>::~STAMethod()
{
}

template<typename FloatType>
void
STAMethod<FloatType>::useCoherenceFactor(bool calculateEnvelope, const std::string& outputDir)
{
	LOG_DEBUG << "Saving the image factors...";
	project_.saveHDF5(gridData_, outputDir + "/image_factor", "image", Util::CopyFactorOp());

	if (calculateEnvelope) {
		ParallelHilbertEnvelope<FloatType>::calculateDim2(gridData_);
	}

	// Applies the coherence factor method.
	for (auto iter = gridData_.begin(); iter != gridData_.end(); ++iter) {
		iter->value *= iter->factor;
		iter->factor = 1.0;
	}

	LOG_DEBUG << "Saving the CF image...";
	project_.saveHDF5(gridData_, outputDir_ + "/image_cf", "image", Util::CopyValueOp());

	Util::copyXZValue(gridData_, projGridData_);
	project_.showFigure3D(2, "Coherence factor image", &projGridData_, &pointList_,
				true, Figure::VISUALIZATION_RECTIFIED_LOG, Figure::COLORMAP_VIRIDIS);
}

template<typename FloatType>
void
STAMethod<FloatType>::process(STAProcessor<FloatType>& processor, unsigned int baseElement, const std::string& outputDir)
{
	Timer tProc;

	processor.process(baseElement, gridData_);

	project_.saveImageToHDF5(gridData_, outputDir);

	Util::copyXZValue(gridData_, projGridData_);
	project_.showFigure3D(1, "Raw image", &projGridData_, &pointList_,
				true, visual_, Figure::COLORMAP_VIRIDIS);

	LOG_DEBUG << ">>> Acquisition + processing time: " << tProc.getTime();
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
		acquisition = std::make_unique<SavedSTAAcquisition<FloatType>>(
					project_, config.numElements,
					taskPM->value<std::string>("data_dir"));
		break;
	case MethodType::sta_simulated_3d:                 // falls through
	case MethodType::sta_simulated_3d_save_signals:    // falls through
	case MethodType::sta_vectorial_simulated_3d:       // falls through
	case MethodType::sta_vectorial_simulated_3d_seq_y:
		acquisition = std::make_unique<Simulated3DSTAAcquisition<FloatType>>(project_, config);
		break;
	default:
		THROW_EXCEPTION(InvalidParameterException, "Invalid method: " << static_cast<int>(project_.method()) << '.');
	}

	if (project_.method() == MethodType::sta_save_signals ||
			project_.method() == MethodType::sta_simulated_3d_save_signals) {
		std::string dataDir = taskPM->value<std::string>("data_dir");
		typename STAAcquisition<FloatType>::AcquisitionDataType acqData;
		for (unsigned int txElem = config.firstTxElem; txElem <= config.lastTxElem; ++txElem) {
			acquisition->execute(baseElement, txElem, acqData);
			project_.saveSTASignalsToHDF5(acqData, dataDir, 0, baseElement, txElem);
		}
		return;
	}

	const FloatType peakOffset  = taskPM->value<FloatType>(  "peak_offset" , 0.0, 50.0);
	const std::string outputDir = taskPM->value<std::string>("output_dir");
	project_.createDirectory(outputDir, false);

	const FloatType nyquistRate = 2.0 * config.maxFrequency;
	const FloatType nyquistLambda = config.propagationSpeed / nyquistRate;
	ImageGrid<FloatType>::get(project_.loadChildParameterMap(taskPM, "grid_config_file"), nyquistLambda, gridData_);

	visual_ = Figure::VISUALIZATION_ENVELOPE_LOG;

	switch (project_.method()) {
	case MethodType::sta_simple_simulated: // falls through
	case MethodType::sta_simple_saved:
		{
			auto processor = std::make_unique<SimpleSTAProcessor<FloatType>>(config, *acquisition, peakOffset);
			process(*processor, baseElement, outputDir);
		}
		break;
	case MethodType::sta_vectorial_dp_saved:     // falls through
	case MethodType::sta_vectorial_sp_saved:     // falls through
	case MethodType::sta_vectorial_simulated_3d:
		{
			bool vectorialProcessingWithEnvelope = taskPM->value<bool>(        "calculate_envelope_in_processing");
			const unsigned int upsamplingFactor  = taskPM->value<unsigned int>("upsampling_factor", 1, 128);
			if (vectorialProcessingWithEnvelope) {
				visual_ = Figure::VISUALIZATION_RECTIFIED_LOG;
			}
			AnalyticSignalCoherenceFactorProcessor<FloatType> coherenceFactor(project_.loadChildParameterMap(taskPM, "coherence_factor_config_file"));
			auto processor = std::make_unique<VectorialSTAProcessor<FloatType>>(
							config, *acquisition, upsamplingFactor,
							coherenceFactor, peakOffset, vectorialProcessingWithEnvelope);
			process(*processor, baseElement, outputDir);
			if (coherenceFactor.enabled()) {
				useCoherenceFactor(!vectorialProcessingWithEnvelope, outputDir);
			}
		}
		break;
	default:
		{
			CoherenceFactorProcessor<FloatType> coherenceFactor(project_.loadChildParameterMap(taskPM, "coherence_factor_config_file"));
			auto processor = std::make_unique<DefaultSTAProcessor<FloatType>>(config, *acquisition,
												coherenceFactor, peakOffset);
			process(*processor, baseElement, outputDir);
			if (coherenceFactor.enabled()) {
				useCoherenceFactor(true, outputDir);
			}
		}
	}
}

} // namespace Lab

#endif /* STAMETHOD_H_ */
