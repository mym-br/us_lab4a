/***************************************************************************
 *  Copyright 2014, 2017, 2018, 2020 Marcelo Y. Matuda                     *
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
#include "Colormap.h"
#include "DefaultSTAProcessor.h"
#include "Exception.h"
#include "ExecutionTimeMeasurement.h"
#include "FileUtil.h"
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
#include "Visualization.h"
#include "Util.h"
#include "XYZ.h"
#include "XYZValueFactor.h"
#ifdef USE_CUDA
# include "VectorialSTACUDAProcessor.h"
#endif



namespace Lab {

template<typename TFloat>
class STAMethod : public Method {
public:
	STAMethod(Project& project);
	virtual ~STAMethod() = default;

	virtual void execute();
private:
	STAMethod(const STAMethod&) = delete;
	STAMethod& operator=(const STAMethod&) = delete;
	STAMethod(STAMethod&&) = delete;
	STAMethod& operator=(STAMethod&&) = delete;

	void process(TFloat valueScale, ArrayProcessor<XYZValueFactor<TFloat>>& processor, unsigned int baseElement,
			const std::string& outputDir);
	template<typename P>
		void process(TFloat valueScale, P& processor, unsigned int baseElement,
				const STAConfiguration<TFloat>& config, const std::string& outputDir);
	void useCoherenceFactor(TFloat valueScale, bool calculateEnvelope, const std::string& outputDir);

	Project& project_;
	Matrix<XYZValueFactor<TFloat>> gridData_;
	std::vector<XYZ<float>> pointList_;
	Visualization::Value visual_;
};



template<typename TFloat>
STAMethod<TFloat>::STAMethod(Project& project)
		: project_(project)
		, pointList_{{0.0, 0.0, 0.0}}
		, visual_(Visualization::VALUE_ENVELOPE_LOG)
{
}

template<typename TFloat>
void
STAMethod<TFloat>::useCoherenceFactor(TFloat valueScale, bool calculateEnvelope, const std::string& outputDir)
{
	project_.saveFactorToHDF5(gridData_, outputDir, "image_factor", "factor");

	if (calculateEnvelope) {
		ParallelHilbertEnvelope<TFloat>::calculateDim2(gridData_);
	}

	Util::applyFactorToValue(gridData_.begin(), gridData_.end());

	project_.saveImageToHDF5(gridData_, outputDir, "image_cf", "cf");

	project_.showFigure3D(2, "Coherence factor image", &gridData_, &pointList_,
				true, Visualization::VALUE_RECTIFIED_LOG, Colormap::GRADIENT_VIRIDIS, valueScale);
}

template<typename TFloat>
void
STAMethod<TFloat>::process(TFloat valueScale, ArrayProcessor<XYZValueFactor<TFloat>>& processor, unsigned int baseElement,
				const std::string& outputDir)
{
	processor.prepare(baseElement);

	processor.process(gridData_);

	project_.saveImageToHDF5(gridData_, outputDir);
	project_.saveXYZToHDF5(gridData_, outputDir);

	project_.showFigure3D(1, "Raw image", &gridData_, &pointList_,
				true, visual_, Colormap::GRADIENT_VIRIDIS, valueScale);
}

template<typename TFloat>
template<typename P>
void
STAMethod<TFloat>::process(TFloat valueScale, P& processor, unsigned int baseElement,
				const STAConfiguration<TFloat>& config, const std::string& outputDir)
{
	processor.prepare(baseElement);

	const unsigned int numTxElem = config.lastTxElem - config.firstTxElem + 1U;
	LOG_DEBUG << "numTxElem: " << numTxElem;

	BEGIN_EXECUTION_TIME_MEASUREMENT_WITH_PARTIAL_X_N(&processor, numTxElem)

	processor.process(gridData_);

	END_EXECUTION_TIME_MEASUREMENT_WITH_PARTIAL_X_N(&processor, numTxElem)

	project_.saveImageToHDF5(gridData_, outputDir);
	project_.saveXYZToHDF5(gridData_, outputDir);

	project_.showFigure3D(1, "Raw image", &gridData_, &pointList_,
				true, visual_, Colormap::GRADIENT_VIRIDIS, valueScale);
}

template<typename TFloat>
void
STAMethod<TFloat>::execute()
{
	const ParameterMap& taskPM = project_.taskParamMap();
	const ParamMapPtr staPM = project_.getSubParamMap("sta_config_file");
	const STAConfiguration<TFloat> config(*staPM);
	const auto baseElement = staPM->value<unsigned int>("base_element", 0, config.numElementsMux - config.numElements);

	std::unique_ptr<STAAcquisition<TFloat>> acquisition;

	switch (project_.method()) {
	case MethodEnum::sta_simple_simulated:
	case MethodEnum::sta_simulated:
		acquisition = std::make_unique<SimulatedSTAAcquisition<TFloat>>(project_, config);
		break;
	case MethodEnum::sta_dp_network:
	case MethodEnum::sta_vectorial_dp_network:
	case MethodEnum::sta_save_signals:
		acquisition = std::make_unique<NetworkSTAAcquisition<TFloat>>(project_, config);
		break;
	case MethodEnum::sta_simple_saved:
	case MethodEnum::sta_dp_saved:
	case MethodEnum::sta_vectorial_dp_saved:
	case MethodEnum::sta_sp_saved:
	case MethodEnum::sta_vectorial_sp_saved:
#ifdef USE_CUDA
	case MethodEnum::sta_vectorial_cuda_sp_saved:
#endif
		acquisition = std::make_unique<SavedSTAAcquisition<TFloat>>(
					project_, config.numElements,
					FileUtil::path(taskPM.value<std::string>("data_dir"), "/", 0));
		break;
	default:
		THROW_EXCEPTION(InvalidParameterException, "Invalid method: " << static_cast<int>(project_.method()) << '.');
	}

	if (project_.method() == MethodEnum::sta_save_signals) {
		const auto dataDir = taskPM.value<std::string>("data_dir");
		typename STAAcquisition<TFloat>::AcquisitionDataType acqData;
		acquisition->prepare(baseElement);
		for (unsigned int txElem = config.firstTxElem; txElem <= config.lastTxElem; ++txElem) {
			acquisition->execute(txElem, acqData);
			project_.saveTxElemSignalsToHDF5(acqData, dataDir, 0, baseElement, txElem);
		}
		return;
	}

	const auto outputDir = taskPM.value<std::string>("output_dir");
	project_.createDirectory(outputDir, false);

	const ParamMapPtr imagPM = project_.getSubParamMap("imag_config_file");
	const auto peakOffset = imagPM->value<TFloat>("peak_offset", 0.0, 50.0);

	const TFloat nyquistLambda = Util::nyquistLambda(config.propagationSpeed, config.maxFrequency);
	ImageGrid<TFloat>::get(*project_.getSubParamMap("grid_config_file"), nyquistLambda, gridData_);

	visual_ = Visualization::VALUE_ENVELOPE_LOG;

	switch (project_.method()) {
	case MethodEnum::sta_simple_simulated:
	case MethodEnum::sta_simple_saved:
		{
			auto processor = std::make_unique<SimpleSTAProcessor<TFloat>>(config, *acquisition, peakOffset);
			process(config.valueScale, *processor, baseElement, outputDir);
		}
		break;
	case MethodEnum::sta_vectorial_dp_network:
	case MethodEnum::sta_vectorial_dp_saved:
	case MethodEnum::sta_vectorial_sp_saved:
		{
			const auto upsamplingFactor = imagPM->value<unsigned int>("upsampling_factor", 1, 128);
			visual_ = Visualization::VALUE_RECTIFIED_LOG;
			AnalyticSignalCoherenceFactorProcessor<TFloat> coherenceFactor(*project_.getSubParamMap("coherence_factor_config_file"));
			std::vector<TFloat> txApod(config.numElements, 1.0); //TODO: Fill txApod
			std::vector<TFloat> rxApod(config.numElements, 1.0); //TODO: Fill rxApod
			auto processor = std::make_unique<VectorialSTAProcessor<TFloat, XYZValueFactor<TFloat>>>(
							config, *acquisition, upsamplingFactor,
							coherenceFactor, peakOffset, txApod, rxApod);
			process(config.valueScale, *processor, baseElement, config, outputDir);
			if (coherenceFactor.enabled()) {
				useCoherenceFactor(config.valueScale, false, outputDir);
			}
		}
		break;
#ifdef USE_CUDA
	case MethodEnum::sta_vectorial_cuda_sp_saved:
		if constexpr (std::is_same<TFloat, float>::value) {
			const auto upsamplingFactor = imagPM->value<unsigned int>("upsampling_factor", 1, 128);
			visual_ = Visualization::VALUE_RECTIFIED_LOG;
			AnalyticSignalCoherenceFactorProcessor<TFloat> coherenceFactor(*project_.getSubParamMap("coherence_factor_config_file"));
			std::vector<TFloat> txApod(config.numElements, 1.0); //TODO: Fill txApod
			std::vector<TFloat> rxApod(config.numElements, 1.0); //TODO: Fill rxApod
			auto processor = std::make_unique<VectorialSTACUDAProcessor>(
							config, *acquisition, upsamplingFactor,
							coherenceFactor, peakOffset, txApod, rxApod);
			process(config.valueScale, *processor, baseElement, config, outputDir);
			if (coherenceFactor.enabled()) {
				useCoherenceFactor(config.valueScale, false, outputDir);
			}
		} else {
			THROW_EXCEPTION(InvalidValueException, "Invalid float type.");
		}
		break;
#endif
	default:
		{
			CoherenceFactorProcessor<TFloat> coherenceFactor(*project_.getSubParamMap("coherence_factor_config_file"));
			auto processor = std::make_unique<DefaultSTAProcessor<TFloat, XYZValueFactor<TFloat>>>(
							config, *acquisition,
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
