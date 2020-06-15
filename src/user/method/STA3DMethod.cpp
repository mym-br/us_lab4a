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

#include "STA3DMethod.h"

#include <cstddef> /* std::size_t */
#include <memory>

#include "CoherenceFactor.h"
#include "Colormap.h"
#include "Exception.h"
#include "ImageGrid.h"
#include "IterationCounter.h"
#include "Log.h"
#include "Project.h"
#include "Simulated3DSTAAcquisition.h"
#include "STAAcquisition.h"
#include "SA3DConfiguration.h"
#include "Timer.h"
#include "Util.h"
#include "Vectorial3DSTAProcessor.h"
#include "Visualization.h"



namespace Lab {

template<typename TFloat>
STA3DMethod<TFloat>::STA3DMethod(Project& project)
		: project_(project)
		, pointList_{{0.0, 0.0, 0.0}}
{
}

template<typename TFloat>
void
STA3DMethod<TFloat>::useCoherenceFactor(TFloat valueScale, const std::string& outputDir)
{
	project_.saveFactorToHDF5(gridData_, outputDir, "image_factor", "factor");

	Util::applyFactorToValue(gridData_.begin(), gridData_.end());

	project_.saveImageToHDF5(gridData_, outputDir, "image_cf", "cf");

	project_.showFigure3D(2, "Coherence factor image", &gridData_, &pointList_,
				true, Visualization::Value::RECTIFIED_LOG, Colormap::Id::GRADIENT_VIRIDIS, valueScale);
}

template<typename TFloat>
void
STA3DMethod<TFloat>::process(TFloat valueScale, ArrayProcessor<XYZValueFactor<TFloat>>& processor, unsigned int baseElement, const std::string& outputDir)
{
	Timer tProc;

	processor.prepare(baseElement);
	processor.process(gridData_);

	project_.saveImageToHDF5(gridData_, outputDir);
	project_.saveXYZToHDF5(gridData_, outputDir);

	project_.showFigure3D(1, "Raw image", &gridData_, &pointList_,
				true, Visualization::Value::RECTIFIED_LOG, Colormap::Id::GRADIENT_VIRIDIS, valueScale);

	LOG_DEBUG << ">>> Acquisition + processing time: " << tProc.getTime();
}

template<typename TFloat>
void
STA3DMethod<TFloat>::execute()
{
	const ParameterMap& taskPM = project_.taskParamMap();
	const ParamMapPtr saPM    = project_.getSubParamMap("sa_config_file");
	const ParamMapPtr arrayPM = project_.getSubParamMap("array_config_file");
	const SA3DConfiguration<TFloat> config(*saPM, *arrayPM);

	const auto baseElement = saPM->value<unsigned int>("base_element", 0, config.numElementsMux - 1U);

	std::unique_ptr<STAAcquisition<TFloat>> acquisition;

	switch (project_.method()) {
	case MethodEnum::sta_3d_simulated_save_signals:
	case MethodEnum::sta_3d_simulated_seq_y_save_signals:
	case MethodEnum::sta_3d_vectorial_simulated:
		acquisition = std::make_unique<Simulated3DSTAAcquisition<TFloat>>(project_, config);
		break;
	default:
		THROW_EXCEPTION(InvalidParameterException, "Invalid method: " << static_cast<int>(project_.method()) << '.');
	}

	if (project_.method() == MethodEnum::sta_3d_simulated_save_signals) {
		const auto dataDir = taskPM.value<std::string>("data_dir");
		typename STAAcquisition<TFloat>::AcquisitionDataType acqData;
		IterationCounter::reset(config.activeTxElem.size());
		acquisition->prepare(baseElement);
		for (unsigned int txElem : config.activeTxElem) {
			acquisition->execute(txElem, acqData);
			project_.saveTxElemSignalsToHDF5(acqData, dataDir, 0, baseElement, txElem);
			IterationCounter::add(1);
		}
		return;
	} else if (project_.method() == MethodEnum::sta_3d_simulated_seq_y_save_signals) {
		const auto dataDir = taskPM.value<std::string>("data_dir");
		const ParamMapPtr seqYCylPM = project_.getSubParamMap("seq_y_cyl_config_file");
		const auto yStep = seqYCylPM->value<TFloat>("y_step",          0.0,   100.0);
		const auto minY  = seqYCylPM->value<TFloat>("min_y" ,     -10000.0, 10000.0);
		const auto maxY  = seqYCylPM->value<TFloat>("max_y" , minY + yStep, 10000.0);
		project_.createDirectory(dataDir, true);
		typename STAAcquisition<TFloat>::AcquisitionDataType acqData;
		std::vector<TFloat> yList;
		Util::fillSequenceFromStartWithStep(yList, minY, maxY, yStep);
		auto& simAcq = dynamic_cast<Simulated3DSTAAcquisition<TFloat>&>(*acquisition);
		IterationCounter::reset(yList.size());
		acquisition->prepare(baseElement);
		for (std::size_t i = 0, end = yList.size(); i < end; ++i) {
			const TFloat y = yList[i];
			simAcq.modifyReflectorsOffset(0.0, -y);
			for (unsigned int txElem : config.activeTxElem) {
				acquisition->execute(txElem, acqData);
				project_.saveTxElemSignalsToHDF5(acqData, dataDir, i, baseElement, txElem);
			}
			IterationCounter::add(1);
		}
		return;
	}

	const auto outputDir = taskPM.value<std::string>("output_dir");
	project_.createDirectory(outputDir, false);

	const TFloat nyquistLambda = Util::nyquistLambda(config.propagationSpeed, config.maxFrequency);
	ImageGrid<TFloat>::get(*project_.getSubParamMap("grid_config_file"), nyquistLambda, gridData_);

	if (project_.method() == MethodEnum::sta_3d_vectorial_simulated) {
		const ParamMapPtr imagPM = project_.getSubParamMap("imag_config_file");
		const auto upsamplingFactor = imagPM->value<unsigned int>("upsampling_factor", 1, 128);
		const auto peakOffset       = imagPM->value<TFloat>(      "peak_offset", 0.0, 50.0);

		std::vector<TFloat> rxApod;
		const auto rxApodFile = imagPM->value<std::string>("rx_apodization_file");
		project_.loadHDF5(rxApodFile, "apod", rxApod);

		AnalyticSignalCoherenceFactorProcessor<TFloat> coherenceFactor(*project_.getSubParamMap("coherence_factor_config_file"));
		auto processor = std::make_unique<Vectorial3DSTAProcessor<TFloat>>(
							config, *acquisition, upsamplingFactor,
							coherenceFactor, peakOffset,
							rxApod);
		process(config.valueScale, *processor, baseElement, outputDir);
		if (coherenceFactor.enabled()) {
			useCoherenceFactor(config.valueScale, outputDir);
		}
	} else {
		THROW_EXCEPTION(InvalidParameterException, "Invalid method: " << static_cast<int>(project_.method()) << '.');
	}
}

// Explicit instantiation.
template class STA3DMethod<float>;
template class STA3DMethod<double>;

} // namespace Lab
