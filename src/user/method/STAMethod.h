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

#include <cstddef> /* std::size_t */
#include <memory>
#include <string>
#include <vector>

#include "CoherenceFactor.h"
#include "DefaultSTAProcessor.h"
#include "Exception.h"
#include "FileUtil.h"
#include "global.h"
#include "ImageGrid.h"
#include "Log.h"
#include "Matrix.h"
#include "Method.h"
#include "NetworkSTAAcquisition.h"
#include "ParallelHilbertEnvelope.h"
#include "Project.h"
#include "SavedSTAAcquisition.h"
#include "SimpleSTAProcessor.h"
#include "SimulatedSTAAcquisition.h"
#include "STAAcquisition.h"
#include "STAConfiguration.h"
#include "Timer.h"
#include "VectorialSTAProcessor.h"
#include "Util.h"
#include "XYZ.h"
#include "XYZValueFactor.h"



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

	void process(FloatType valueScale, ArrayProcessor<FloatType>& processor, unsigned int baseElement, const std::string& outputDir);
	void useCoherenceFactor(FloatType valueScale, bool calculateEnvelope, const std::string& outputDir);

	Project& project_;
	Matrix<XYZValueFactor<FloatType>> gridData_;
	std::vector<XYZ<float>> pointList_;
	Figure::Visualization visual_;
};



template<typename FloatType>
STAMethod<FloatType>::STAMethod(Project& project)
		: project_(project)
		, pointList_{{0.0, 0.0, 0.0}}
		, visual_(Figure::VISUALIZATION_ENVELOPE_LOG)
{
}

template<typename FloatType>
STAMethod<FloatType>::~STAMethod()
{
}

template<typename FloatType>
void
STAMethod<FloatType>::useCoherenceFactor(FloatType valueScale, bool calculateEnvelope, const std::string& outputDir)
{
	LOG_DEBUG << "Saving the image factors...";
	project_.saveHDF5(gridData_, outputDir + "/image_factor", "factor", Util::CopyFactorOp());

	if (calculateEnvelope) {
		ParallelHilbertEnvelope<FloatType>::calculateDim2(gridData_);
	}

	// Applies the coherence factor method.
	for (auto iter = gridData_.begin(); iter != gridData_.end(); ++iter) {
		iter->value *= iter->factor;
		iter->factor = 1.0;
	}

	LOG_DEBUG << "Saving the CF image...";
	project_.saveHDF5(gridData_, outputDir + "/image_cf", "cf", Util::CopyValueOp());

	project_.showFigure3D(2, "Coherence factor image", &gridData_, &pointList_,
				true, Figure::VISUALIZATION_RECTIFIED_LOG, Figure::COLORMAP_VIRIDIS, valueScale);
}

template<typename FloatType>
void
STAMethod<FloatType>::process(FloatType valueScale, ArrayProcessor<FloatType>& processor, unsigned int baseElement, const std::string& outputDir)
{
	Timer tProc;

	processor.prepare(baseElement);
	processor.process(gridData_);

	project_.saveImageToHDF5(gridData_, outputDir);

	project_.showFigure3D(1, "Raw image", &gridData_, &pointList_,
				true, visual_, Figure::COLORMAP_VIRIDIS, valueScale);

	LOG_DEBUG << ">>> Acquisition + processing time: " << tProc.getTime();
}

template<typename FloatType>
void
STAMethod<FloatType>::execute()
{
	ParamMapPtr taskPM = project_.taskParameterMap();
	ParamMapPtr staPM = project_.loadChildParameterMap(taskPM, "sta_config_file");
	const STAConfiguration<FloatType> config(staPM);
	const auto baseElement = staPM->value<unsigned int>("base_element", 0, config.numElementsMux - config.numElements);

	std::unique_ptr<STAAcquisition<FloatType>> acquisition;

	switch (project_.method()) {
	case MethodEnum::sta_simple_simulated: // falls through
	case MethodEnum::sta_simulated:
		acquisition = std::make_unique<SimulatedSTAAcquisition<FloatType>>(project_, config);
		break;
	case MethodEnum::sta_dp_network:           // falls through
	case MethodEnum::sta_vectorial_dp_network: // falls through
	case MethodEnum::sta_save_signals:
		acquisition = std::make_unique<NetworkSTAAcquisition<FloatType>>(project_, config);
		break;
	case MethodEnum::sta_simple_saved:       // falls through
	case MethodEnum::sta_dp_saved:           // falls through
	case MethodEnum::sta_vectorial_dp_saved: // falls through
	case MethodEnum::sta_sp_saved:           // falls through
	case MethodEnum::sta_vectorial_sp_saved:
		acquisition = std::make_unique<SavedSTAAcquisition<FloatType>>(
					project_, config.numElements,
					FileUtil::path(taskPM->value<std::string>("data_dir"), "/", 0));
		break;
	default:
		THROW_EXCEPTION(InvalidParameterException, "Invalid method: " << static_cast<int>(project_.method()) << '.');
	}

	if (project_.method() == MethodEnum::sta_save_signals) {
		const auto dataDir = taskPM->value<std::string>("data_dir");
		typename STAAcquisition<FloatType>::AcquisitionDataType acqData;
		acquisition->prepare(baseElement);
		for (unsigned int txElem = config.firstTxElem; txElem <= config.lastTxElem; ++txElem) {
			acquisition->execute(txElem, acqData);
			project_.saveTxElemSignalsToHDF5(acqData, dataDir, 0, baseElement, txElem);
		}
		return;
	}

	const auto outputDir = taskPM->value<std::string>("output_dir");
	project_.createDirectory(outputDir, false);

	ParamMapPtr imagPM = project_.loadChildParameterMap(taskPM, "imag_config_file");
	const auto peakOffset = imagPM->value<FloatType>("peak_offset", 0.0, 50.0);

	const FloatType nyquistRate = 2.0 * config.maxFrequency;
	const FloatType nyquistLambda = config.propagationSpeed / nyquistRate;
	ImageGrid<FloatType>::get(project_.loadChildParameterMap(taskPM, "grid_config_file"), nyquistLambda, gridData_);

	visual_ = Figure::VISUALIZATION_ENVELOPE_LOG;

	switch (project_.method()) {
	case MethodEnum::sta_simple_simulated: // falls through
	case MethodEnum::sta_simple_saved:
		{
			auto processor = std::make_unique<SimpleSTAProcessor<FloatType>>(config, *acquisition, peakOffset);
			process(config.valueScale, *processor, baseElement, outputDir);
		}
		break;
	case MethodEnum::sta_vectorial_dp_network: // falls through
	case MethodEnum::sta_vectorial_dp_saved:   // falls through
	case MethodEnum::sta_vectorial_sp_saved:
		{
			const auto processingWithEnvelope = imagPM->value<bool>(        "calculate_envelope_in_processing");
			const auto upsamplingFactor       = imagPM->value<unsigned int>("upsampling_factor", 1, 128);
			if (processingWithEnvelope) {
				visual_ = Figure::VISUALIZATION_RECTIFIED_LOG;
			}
			AnalyticSignalCoherenceFactorProcessor<FloatType> coherenceFactor(project_.loadChildParameterMap(taskPM, "coherence_factor_config_file"));
			std::vector<FloatType> txApod(config.numElements, 1.0);
			std::vector<FloatType> rxApod(config.numElements, 1.0);
			auto processor = std::make_unique<VectorialSTAProcessor<FloatType>>(
							config, *acquisition, upsamplingFactor,
							coherenceFactor, peakOffset, processingWithEnvelope, txApod, rxApod);
			process(config.valueScale, *processor, baseElement, outputDir);
			if (coherenceFactor.enabled()) {
				useCoherenceFactor(config.valueScale, !processingWithEnvelope, outputDir);
			}
		}
		break;
	default:
		{
			CoherenceFactorProcessor<FloatType> coherenceFactor(project_.loadChildParameterMap(taskPM, "coherence_factor_config_file"));
			auto processor = std::make_unique<DefaultSTAProcessor<FloatType>>(config, *acquisition,
												coherenceFactor, peakOffset);
			process(config.valueScale, *processor, baseElement, outputDir);
			if (coherenceFactor.enabled()) {
				useCoherenceFactor(config.valueScale, true, outputDir);
			}
		}
	}
}

} // namespace Lab

#endif /* STAMETHOD_H_ */
