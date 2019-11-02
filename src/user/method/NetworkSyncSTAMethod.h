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

#ifndef NETWORKSYNCSTAMETHOD_H
#define NETWORKSYNCSTAMETHOD_H

#include <algorithm> /* for_each */
#include <memory>
#include <string>
#include <vector>

#include "CoherenceFactor.h"
#include "FileUtil.h"
#include "ImageGrid.h"
#include "Log.h"
#include "Matrix.h"
#include "Method.h"
#include "NetworkSTAAcquisition.h"
#include "ParallelHilbertEnvelope.h"
#include "Project.h"
#include "SavedSTAAcquisition.h"
#include "STAAcquisition.h"
#include "STAConfiguration.h"
#include "SyncServer.h"
#include "Timer.h"
#include "VectorialSTAProcessor.h"
#include "Util.h"
#include "WindowFunction.h"
#include "XYZ.h"
#include "XYZValueFactor.h"



namespace Lab {

template<typename FloatType>
class NetworkSyncSTAMethod : public Method {
public:
	NetworkSyncSTAMethod(Project& project);
	virtual ~NetworkSyncSTAMethod();

	virtual void execute();

private:
	static constexpr unsigned int maxSteps = 10000;
	static constexpr const char* timeFile    = "/time";
	static constexpr const char* timeDataset = "time";
	static constexpr const char* yFile    = "/y";
	static constexpr const char* yDataset = "y";

	NetworkSyncSTAMethod(const NetworkSyncSTAMethod&) = delete;
	NetworkSyncSTAMethod& operator=(const NetworkSyncSTAMethod&) = delete;

	void saveSignals(ParamMapPtr taskPM, const STAConfiguration<FloatType>& config, STAAcquisition<FloatType>& acq,
				unsigned int baseElement, const std::string& dataDir);

	Project& project_;
};



template<typename FloatType>
NetworkSyncSTAMethod<FloatType>::NetworkSyncSTAMethod(Project& project)
		: project_(project)
{
}

template<typename FloatType>
NetworkSyncSTAMethod<FloatType>::~NetworkSyncSTAMethod()
{
}

template<typename FloatType>
void
NetworkSyncSTAMethod<FloatType>::execute()
{
	ParamMapPtr taskPM = project_.taskParameterMap();
	ParamMapPtr staPM = project_.loadChildParameterMap(taskPM, "sta_config_file");
	const STAConfiguration<FloatType> config(staPM);
	const auto baseElement = staPM->value<unsigned int>("base_element", 0, config.numElementsMux - config.numElements);
	const auto dataDir = taskPM->value<std::string>("data_dir");

	if (project_.method() == MethodEnum::sta_network_sync_save_signals) {
		project_.createDirectory(dataDir, true);
		auto acquisition = std::make_unique<NetworkSTAAcquisition<FloatType>>(project_, config);
		saveSignals(taskPM, config, *acquisition, baseElement, dataDir);
		return;
	}

	const auto outputDir = taskPM->value<std::string>("output_dir");
	ParamMapPtr imagPM = project_.loadChildParameterMap(taskPM, "imag_config_file");
	const auto peakOffset                      = imagPM->value<FloatType>(   "peak_offset"      , 0.0, 50.0);
	const auto vectorialProcessingWithEnvelope = imagPM->value<bool>(        "calculate_envelope_in_processing");
	const auto upsamplingFactor                = imagPM->value<unsigned int>("upsampling_factor",   1,  128);
	const auto txApodDesc                      = imagPM->value<std::string>( "tx_apodization");
	const auto rxApodDesc                      = imagPM->value<std::string>( "rx_apodization");

	std::vector<FloatType> txApod(config.numElements);
	WindowFunction::get(txApodDesc, config.numElements, txApod);
	std::vector<FloatType> rxApod(config.numElements);
	WindowFunction::get(rxApodDesc, config.numElements, rxApod);
	std::vector<FloatType> elemIndex;
	Util::fillSequenceFromStartToEndWithSize(elemIndex,
		FloatType(0), static_cast<FloatType>(config.numElements - 1U), static_cast<FloatType>(config.numElements));
	project_.showFigure2D(1, "TX Apodization", elemIndex, txApod, true, true);
	project_.showFigure2D(2, "RX Apodization", elemIndex, rxApod, true, true);

	project_.createDirectory(outputDir, true);

	Matrix<XYZValueFactor<FloatType>> gridData;

	const FloatType nyquistRate = 2.0 * config.maxFrequency;
	const FloatType nyquistLambda = config.propagationSpeed / nyquistRate;
	ImageGrid<FloatType>::get(project_.loadChildParameterMap(taskPM, "grid_config_file"), nyquistLambda, gridData);

	AnalyticSignalCoherenceFactorProcessor<FloatType> coherenceFactor(project_.loadChildParameterMap(taskPM, "coherence_factor_config_file"));
	bool coherenceFactorEnabled = coherenceFactor.enabled();
	auto acquisition = std::make_unique<SavedSTAAcquisition<FloatType>>(project_, config.numElements, "");
	auto processor = std::make_unique<VectorialSTAProcessor<FloatType>>(config, *acquisition,
					upsamplingFactor, coherenceFactor, peakOffset,
					vectorialProcessingWithEnvelope, txApod, rxApod);
	Figure::Visualization visual;
	if (vectorialProcessingWithEnvelope) {
		visual = Figure::VISUALIZATION_RECTIFIED_LOG;
	} else {
		visual = Figure::VISUALIZATION_ENVELOPE_LOG;
	}
	std::vector<XYZ<float>> pointList = {{0.0, 0.0, 0.0}};
	FloatType valueLevel = 0.0;

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
				Util::MultiplyValueBy<XYZValueFactor<FloatType>, FloatType>(config.valueScale));
		}

		project_.saveImageToHDF5(gridData, acqOutputDir);
		project_.saveXYZToHDF5(gridData, acqOutputDir);
		const FloatType maxAbsValue = Util::maxAbsoluteValueField<XYZValueFactor<FloatType>, FloatType>(gridData);
		if (maxAbsValue > valueLevel) valueLevel = maxAbsValue;

		project_.showFigure3D(1, "Raw image", &gridData, &pointList,
					true, visual, Figure::COLORMAP_VIRIDIS,
					config.valueScale != 0.0 ? 1.0 : 0.0);

		if (coherenceFactorEnabled) {
			project_.saveFactorToHDF5(gridData, acqOutputDir, "image_factor", "factor");

			if (!vectorialProcessingWithEnvelope) {
				ParallelHilbertEnvelope<FloatType>::calculateDim2(gridData);
			}

			Util::applyFactorToValue(gridData.begin(), gridData.end());

			project_.saveImageToHDF5(gridData, acqOutputDir, "image_cf", "cf");

			project_.showFigure3D(2, "Coherence factor image", &gridData, &pointList,
						true, Figure::VISUALIZATION_RECTIFIED_LOG, Figure::COLORMAP_VIRIDIS,
						config.valueScale != 0.0 ? 1.0 : 0.0);
		}
	}

	LOG_INFO << "=== value level: " << valueLevel << " (" << Util::linearToDecibels(valueLevel) << " dB)";

	LOG_DEBUG << ">>> Acquisition + processing + saving data time: " << timer.getTime();
}

template<typename FloatType>
void
NetworkSyncSTAMethod<FloatType>::saveSignals(ParamMapPtr taskPM, const STAConfiguration<FloatType>& config, STAAcquisition<FloatType>& acq,
						unsigned int baseElement, const std::string& dataDir)
{
	ParamMapPtr scanPM = project_.loadChildParameterMap(taskPM, "scan_config_file");
	const auto serverPort = scanPM->value<unsigned int>("sync_server_port",     1024,   65535);
	const auto minY       = scanPM->value<FloatType>(   "min_y"           , -10000.0, 10000.0);
	const auto yStep      = scanPM->value<FloatType>(   "y_step"          ,   1.0e-6,  1000.0);

	SyncServer server(serverPort);
	typename STAAcquisition<FloatType>::AcquisitionDataType acqData;
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

} // namespace Lab

#endif // NETWORKSYNCSTAMETHOD_H
