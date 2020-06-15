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

#include "NetworkSyncSTAMethod.h"

#include <algorithm> /* for_each */
#include <memory>
#include <vector>

#include "CoherenceFactor.h"
#include "Colormap.h"
#include "FileUtil.h"
#include "ImageGrid.h"
#include "Log.h"
#include "Matrix.h"
#include "NetworkSTAAcquisition.h"
#include "Project.h"
#include "SavedSTAAcquisition.h"
#include "STAAcquisition.h"
#include "STAConfiguration.h"
#include "SyncServer.h"
#include "Timer.h"
#include "Util.h"
#include "VectorialSTAProcessor.h"
#include "Visualization.h"
#include "WindowFunction.h"
#include "XYZ.h"
#include "XYZValueFactor.h"



namespace Lab {

template<typename TFloat>
NetworkSyncSTAMethod<TFloat>::NetworkSyncSTAMethod(Project& project)
		: project_(project)
{
}

template<typename TFloat>
void
NetworkSyncSTAMethod<TFloat>::execute()
{
	const ParameterMap& taskPM = project_.taskParamMap();
	const ParamMapPtr staPM = project_.getSubParamMap("sta_config_file");
	const STAConfiguration<TFloat> config(*staPM);
	const auto baseElement = staPM->value<unsigned int>("base_element", 0, config.numElementsMux - config.numElements);
	const auto dataDir = taskPM.value<std::string>("data_dir");

	if (project_.method() == MethodEnum::sta_network_sync_save_signals) {
		project_.createDirectory(dataDir, true);
		auto acquisition = std::make_unique<NetworkSTAAcquisition<TFloat>>(project_, config);
		saveSignals(config, *acquisition, baseElement, dataDir);
		return;
	}

	const auto outputDir = taskPM.value<std::string>("output_dir");
	const ParamMapPtr imagPM = project_.getSubParamMap("imag_config_file");
	const auto peakOffset       = imagPM->value<TFloat>(      "peak_offset"      , 0.0, 50.0);
	const auto upsamplingFactor = imagPM->value<unsigned int>("upsampling_factor",   1,  128);
	const auto rxApodDesc       = imagPM->value<std::string>( "rx_apodization");

	std::vector<TFloat> rxApod(config.numElements);
	WindowFunction::get(rxApodDesc, config.numElements, rxApod);
	std::vector<TFloat> elemIndex;
	Util::fillSequenceFromStartToEndWithSize(elemIndex,
		TFloat(0), static_cast<TFloat>(config.numElements - 1U), static_cast<TFloat>(config.numElements));
	project_.showFigure2D(1, "RX Apodization", elemIndex, rxApod, true, true);

	project_.createDirectory(outputDir, true);

	Matrix<XYZValueFactor<TFloat>> gridData;

	const TFloat nyquistLambda = Util::nyquistLambda(config.propagationSpeed, config.maxFrequency);
	ImageGrid<TFloat>::get(*project_.getSubParamMap("grid_config_file"), nyquistLambda, gridData);

	AnalyticSignalCoherenceFactorProcessor<TFloat> coherenceFactor(*project_.getSubParamMap("coherence_factor_config_file"));
	bool coherenceFactorEnabled = coherenceFactor.enabled();
	auto acquisition = std::make_unique<SavedSTAAcquisition<TFloat>>(project_, config.numElements, "");
	auto processor = std::make_unique<VectorialSTAProcessor<TFloat, XYZValueFactor<TFloat>>>(
					config, *acquisition,
					upsamplingFactor, coherenceFactor, peakOffset,
					rxApod);
	std::vector<XYZ<float>> pointList = {{0.0, 0.0, 0.0}};
	TFloat valueLevel = 0.0;

	// Load y.
	std::vector<double> yList;
	const std::string yFileName = dataDir + yFile;
	project_.loadHDF5(yFileName, yDataset, yList);

	Timer timer;

	processor->prepare(baseElement);
	for (unsigned int acqNumber = 0; acqNumber < yList.size(); ++acqNumber) {
		std::string acqDataDir = FileUtil::path(dataDir, "/", acqNumber);
		if (!project_.directoryExists(acqDataDir)) {
			break;
		}

		LOG_INFO << "ACQ " << acqNumber;

		for (auto iter = gridData.begin(); iter != gridData.end(); ++iter) {
			iter->y = yList[acqNumber];
		}

		acquisition->setDataDir(acqDataDir);
		processor->process(gridData);

		std::string acqOutputDir = FileUtil::path(outputDir, "/", acqNumber);
		project_.createDirectory(acqOutputDir, true);

		if (config.valueScale != 0.0) {
			std::for_each(gridData.begin(), gridData.end(),
				Util::MultiplyValueBy<XYZValueFactor<TFloat>, TFloat>(config.valueScale));
		}

		project_.saveImageToHDF5(gridData, acqOutputDir);
		project_.saveXYZToHDF5(gridData, acqOutputDir);
		const TFloat maxAbsValue = Util::maxAbsoluteValueField<XYZValueFactor<TFloat>, TFloat>(gridData);
		if (maxAbsValue > valueLevel) valueLevel = maxAbsValue;

		project_.showFigure3D(1, "Raw image", &gridData, &pointList,
					true, Visualization::Value::RECTIFIED_LOG, Colormap::Id::GRADIENT_VIRIDIS,
					config.valueScale != 0.0 ? 1.0 : 0.0);

		if (coherenceFactorEnabled) {
			project_.saveFactorToHDF5(gridData, acqOutputDir, "image_factor", "factor");

			Util::applyFactorToValue(gridData.begin(), gridData.end());

			project_.saveImageToHDF5(gridData, acqOutputDir, "image_cf", "cf");

			project_.showFigure3D(2, "Coherence factor image", &gridData, &pointList,
						true, Visualization::Value::RECTIFIED_LOG, Colormap::Id::GRADIENT_VIRIDIS,
						config.valueScale != 0.0 ? 1.0 : 0.0);
		}
	}

	LOG_INFO << "=== value level: " << valueLevel << " (" << Util::linearToDecibels(valueLevel) << " dB)";

	LOG_DEBUG << ">>> Acquisition + processing + saving data time: " << timer.getTime();
}

template<typename TFloat>
void
NetworkSyncSTAMethod<TFloat>::saveSignals(const STAConfiguration<TFloat>& config, STAAcquisition<TFloat>& acq,
						unsigned int baseElement, const std::string& dataDir)
{
	const ParamMapPtr scanPM = project_.getSubParamMap("scan_config_file");
	const auto serverPort = scanPM->value<unsigned int>("sync_server_port",     1024,   65535);
	const auto minY       = scanPM->value<TFloat>(      "min_y"           , -10000.0, 10000.0);
	const auto yStep      = scanPM->value<TFloat>(      "y_step"          ,   1.0e-6,  1000.0);

	SyncServer server(serverPort);
	typename STAAcquisition<TFloat>::AcquisitionDataType acqData;
	std::vector<double> timeList;
	timeList.reserve(1000);
	std::vector<double> yList;
	yList.reserve(1000);

	// Capture and save signals.
	Timer timer;
	unsigned int acqNumber = 0;
	const double t0 = timer.getTime();
	acq.prepare(baseElement);
	while (true) {
		LOG_DEBUG << "Waiting for trigger...";
		if (!server.waitForTrigger()) break; // trigger abort
		if (acqNumber == maxSteps) break;
		LOG_INFO << "ACQ " << acqNumber;

		for (unsigned int txElem = config.firstTxElem; txElem <= config.lastTxElem; ++txElem) {
			timeList.push_back(timer.getTime() - t0);
			yList.push_back(minY + acqNumber * yStep); // fixed step
			acq.execute(txElem, acqData);
			project_.saveTxElemSignalsToHDF5(acqData, dataDir, acqNumber, baseElement, txElem);
		}

		server.freeTrigger(); // let the array move to the next position

		++acqNumber;
	}

	// Save times.
	const std::string timeFileName = dataDir + timeFile;
	project_.saveHDF5(timeList, timeFileName, timeDataset);
	// Save y.
	const std::string yFileName = dataDir + yFile;
	project_.saveHDF5(yList, yFileName, yDataset);
}

// Explicit instantiation.
template class NetworkSyncSTAMethod<float>;
template class NetworkSyncSTAMethod<double>;

} // namespace Lab
