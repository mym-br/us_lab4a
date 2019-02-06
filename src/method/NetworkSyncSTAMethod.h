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

#include <memory>
#include <string>
#include <vector>

#include "CoherenceFactor.h"
#include "FileUtil.h"
#include "ImageGrid.h"
#include "Log.h"
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
	NetworkSyncSTAMethod(const NetworkSyncSTAMethod&) = delete;
	NetworkSyncSTAMethod& operator=(const NetworkSyncSTAMethod&) = delete;

	void saveSignals(const STAConfiguration<FloatType>& config, STAAcquisition<FloatType>& acq,
				unsigned int syncServerPort, unsigned int baseElement, const std::string& dataDir);

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
	ConstParameterMapPtr taskPM = project_.taskParameterMap();
	ConstParameterMapPtr staPM = project_.loadChildParameterMap(taskPM, "sta_config_file");
	const STAConfiguration<FloatType> config(staPM);
	const unsigned int baseElement = staPM->value<unsigned int>("base_element", 0, config.numElementsMux - config.numElements);
	const std::string dataDir      = taskPM->value<std::string>( "data_dir");

	if (project_.method() == MethodType::sta_network_sync_save_signals) {
		auto acquisition = std::make_unique<NetworkSTAAcquisition<FloatType>>(project_, config);
		const unsigned int serverPort = taskPM->value<unsigned int>("sync_server_port" , 1024, 65535);
		project_.createDirectory(dataDir, true);
		saveSignals(config, *acquisition, serverPort, baseElement, dataDir);
		return;
	}

	const FloatType peakOffset           = taskPM->value<FloatType>(   "peak_offset"      ,      0.0,    50.0);
	bool vectorialProcessingWithEnvelope = taskPM->value<bool>(        "calculate_envelope_in_processing");
	const unsigned int upsamplingFactor  = taskPM->value<unsigned int>("upsampling_factor",        1,     128);
	const std::string outputDir          = taskPM->value<std::string>( "output_dir");
	const std::string txApodDesc         = taskPM->value<std::string>( "tx_apodization");
	const std::string rxApodDesc         = taskPM->value<std::string>( "rx_apodization");

	ConstParameterMapPtr scanPM = project_.loadChildParameterMap(taskPM, "scan_config_file");
	const FloatType minY  = scanPM->value<FloatType>("min_y" , -10000.0, 10000.0);
	const FloatType yStep = scanPM->value<FloatType>("y_step",   1.0e-6,  1000.0);

	std::vector<FloatType> txApod(config.numElements);
	WindowFunction::get(txApodDesc, config.numElements, txApod);
	std::vector<FloatType> rxApod(config.numElements);
	WindowFunction::get(rxApodDesc, config.numElements, rxApod);
	std::vector<FloatType> elemIndex;
	Util::fillSequenceFromStartToEndWithSize(elemIndex,
		FloatType{0}, static_cast<FloatType>(config.numElements - 1U), static_cast<FloatType>(config.numElements));
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
	Project::GridDataType projGridData;
	FloatType valueLevel = 0.0;

	Timer timer;

	FloatType y = minY;
	for (unsigned int acqNumber = 0; ; ++acqNumber, y += yStep) {
		std::string acqDataDir = FileUtil::path(dataDir, "/", acqNumber);
		if (!project_.directoryExists(acqDataDir)) {
			break;
		}

		LOG_DEBUG << "ACQ number: " << acqNumber;

		for (auto iter = gridData.begin(); iter != gridData.end(); ++iter) {
			iter->y = y;
		}
		acquisition->setDataDir(acqDataDir);
		processor->process(baseElement, gridData);

		std::string acqOutputDir = FileUtil::path(outputDir, "/", acqNumber);
		project_.createDirectory(acqOutputDir, true);

		project_.saveImageToHDF5(gridData, acqOutputDir);
		const FloatType maxAbsValue = Util::maxAbsoluteValueField<XYZValueFactor<FloatType>, FloatType>(gridData);
		if (maxAbsValue > valueLevel) valueLevel = maxAbsValue;

		Util::copyXYZValue(gridData, projGridData);
		project_.showFigure3D(1, "Raw image", &projGridData, &pointList,
					true, visual, Figure::COLORMAP_VIRIDIS, config.valueScale);

		if (coherenceFactorEnabled) {
			LOG_DEBUG << "Saving the image factors...";
			project_.saveHDF5(gridData, acqOutputDir + "/image_factor", "factor", Util::CopyFactorOp());

			if (!vectorialProcessingWithEnvelope) {
				ParallelHilbertEnvelope<FloatType>::calculateDim2(gridData);
			}

			// Apply the coherence factor method.
			for (auto iter = gridData.begin(); iter != gridData.end(); ++iter) {
				iter->value *= iter->factor;
				iter->factor = 1.0;
			}

			LOG_DEBUG << "Saving the CF image...";
			project_.saveHDF5(gridData, acqOutputDir + "/image_cf", "cf", Util::CopyValueOp());

			Util::copyXYZValue(gridData, projGridData);
			project_.showFigure3D(2, "Coherence factor image", &projGridData, &pointList,
						true, Figure::VISUALIZATION_RECTIFIED_LOG, Figure::COLORMAP_VIRIDIS, config.valueScale);
		}
	}

	LOG_INFO << "=== value level: " << valueLevel << " (" << Util::linearToDecibels(valueLevel) << " dB)";

	LOG_DEBUG << ">>> Acquisition + processing + saving data time: " << timer.getTime();
}

template<typename FloatType>
void
NetworkSyncSTAMethod<FloatType>::saveSignals(const STAConfiguration<FloatType>& config, STAAcquisition<FloatType>& acq,
						unsigned int syncServerPort, unsigned int baseElement, const std::string& dataDir)
{
	SyncServer server(syncServerPort);
	typename STAAcquisition<FloatType>::AcquisitionDataType acqData;

	unsigned int acqNumber = 0;
	while (true) {
		if (!server.waitForTrigger()) break;

		for (unsigned int txElem = config.firstTxElem; txElem <= config.lastTxElem; ++txElem) {
			acq.execute(baseElement, txElem, acqData);
			project_.saveTxElemSignalsToHDF5(acqData, dataDir, acqNumber, baseElement, txElem);
		}

		server.freeTrigger();

		++acqNumber;
	}
}

} // namespace Lab

#endif // NETWORKSYNCSTAMETHOD_H
