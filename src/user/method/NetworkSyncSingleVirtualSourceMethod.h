/***************************************************************************
 *  Copyright 2018, 2019 Marcelo Y. Matuda                                 *
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

#ifndef NETWORK_SYNC_SINGLE_VIRTUAL_SOURCE_METHOD_H
#define NETWORK_SYNC_SINGLE_VIRTUAL_SOURCE_METHOD_H

#include <algorithm> /* for_each */
#include <memory>
#include <string>
#include <vector>

#include "ArrayUtil.h"
#include "CoherenceFactor.h"
#include "FileUtil.h"
#include "ImageGrid.h"
#include "Log.h"
#include "Matrix.h"
#include "Method.h"
#include "NetworkTnRnAcquisition.h"
#include "Project.h"
#include "SavedTnRnAcquisition.h"
#include "SyncServer.h"
#include "Timer.h"
#include "TnRnAcquisition.h"
#include "TnRnConfiguration.h"
#include "Vectorial3DTnRnProcessor.h"
#include "Util.h"
#include "XYZ.h"
#include "XYZValueFactor.h"

#define NETWORK_SYNC_SINGLE_VIRTUAL_SOURCE_METHOD_MAX_STEPS 10000
#define NETWORK_SYNC_SINGLE_VIRTUAL_SOURCE_METHOD_TIME_FILE "/time"
#define NETWORK_SYNC_SINGLE_VIRTUAL_SOURCE_METHOD_TIME_DATASET "time"
#define NETWORK_SYNC_SINGLE_VIRTUAL_SOURCE_METHOD_Y_FILE "/y"
#define NETWORK_SYNC_SINGLE_VIRTUAL_SOURCE_METHOD_Y_DATASET "y"



namespace Lab {

template<typename FloatType>
class NetworkSyncSingleVirtualSourceMethod : public Method {
public:
	NetworkSyncSingleVirtualSourceMethod(Project& project);
	virtual ~NetworkSyncSingleVirtualSourceMethod();

	virtual void execute();

private:
	NetworkSyncSingleVirtualSourceMethod(const NetworkSyncSingleVirtualSourceMethod&) = delete;
	NetworkSyncSingleVirtualSourceMethod& operator=(const NetworkSyncSingleVirtualSourceMethod&) = delete;

	void saveSignals(ConstParameterMapPtr taskPM, TnRnAcquisition<FloatType>& acq,
				unsigned int baseElement, const std::vector<FloatType>& txDelays,
				const std::string& dataDir);

	Project& project_;
};



template<typename FloatType>
NetworkSyncSingleVirtualSourceMethod<FloatType>::NetworkSyncSingleVirtualSourceMethod(Project& project)
		: project_(project)
{
}

template<typename FloatType>
NetworkSyncSingleVirtualSourceMethod<FloatType>::~NetworkSyncSingleVirtualSourceMethod()
{
}

template<typename FloatType>
void
NetworkSyncSingleVirtualSourceMethod<FloatType>::execute()
{
	ConstParameterMapPtr taskPM = project_.taskParameterMap();
	ConstParameterMapPtr saPM    = project_.loadChildParameterMap(taskPM, "sa_config_file");
	ConstParameterMapPtr arrayPM = project_.loadChildParameterMap(taskPM, "array_config_file");
	const TnRnConfiguration<FloatType> config(saPM, arrayPM);
	const unsigned int baseElement = saPM->value<unsigned int>("base_element", 0, config.numElementsMux - config.numElements);
	const FloatType focusZ         = saPM->value<FloatType>(   "tx_focus_z", -10000.0, 10000.0);
	const std::string dataDir      = taskPM->value<std::string>("data_dir");

	FloatType focusX = 0, focusY = 0;
	// Set the focus at the mean x, y.
	for (unsigned int i = baseElement, end = baseElement + config.numElements; i < end; ++i) {
		if (i >= config.txElemPos.size()) {
			THROW_EXCEPTION(InvalidValueException, "Invalid tx element number: " << i << '.');
		}
		const XY<FloatType>& pos = config.txElemPos[i];
		focusX += pos.x;
		focusY += pos.y;
	}
	focusX /= config.numElements;
	focusY /= config.numElements;
	std::vector<FloatType> txDelays;
	ArrayUtil::calculateTx3DFocusDelay(focusX, focusY, focusZ, config.propagationSpeed,
						config.txElemPos, baseElement, config.numElements, txDelays);

	if (project_.method() == MethodEnum::single_virtual_source_network_sync_save_signals) {
		project_.createDirectory(dataDir, true);
		auto acquisition = std::make_unique<NetworkTnRnAcquisition<FloatType>>(project_, config);
		saveSignals(taskPM, *acquisition, baseElement, txDelays, dataDir);
		return;
	}

	const std::string outputDir = taskPM->value<std::string>("output_dir");
	ConstParameterMapPtr imagPM = project_.loadChildParameterMap(taskPM, "imag_config_file");
	const FloatType peakOffset          = imagPM->value<FloatType>(   "peak_offset"      , 0.0, 50.0);
	const unsigned int upsamplingFactor = imagPM->value<unsigned int>("upsampling_factor",   1,  128);

	const std::string rxApodFile = imagPM->value<std::string>("rx_apodization_file");
	std::vector<FloatType> rxApod;
	project_.loadHDF5(rxApodFile, "apod", rxApod);

	project_.createDirectory(outputDir, true);

	Matrix<XYZValueFactor<FloatType>> gridData;

	const FloatType nyquistRate = 2.0 * config.maxFrequency;
	const FloatType nyquistLambda = config.propagationSpeed / nyquistRate;
	ImageGrid<FloatType>::get(project_.loadChildParameterMap(taskPM, "grid_config_file"), nyquistLambda, gridData);

	AnalyticSignalCoherenceFactorProcessor<FloatType> coherenceFactor(project_.loadChildParameterMap(taskPM, "coherence_factor_config_file"));
	bool coherenceFactorEnabled = coherenceFactor.enabled();
	auto acquisition = std::make_unique<SavedTnRnAcquisition<FloatType>>(project_, config.numElements, "");
	auto processor = std::make_unique<Vectorial3DTnRnProcessor<FloatType>>(config, *acquisition,
					upsamplingFactor, coherenceFactor, peakOffset, rxApod);
	processor->setTxDelays(focusX, focusY, focusZ, txDelays);

	std::vector<XYZ<float>> pointList = {{0.0, 0.0, 0.0}};
	FloatType valueLevel = 0.0;
	FloatType cfValueLevel = 0.0;

	// Load y (acquisition position).
	std::vector<double> yList;
	const std::string yFileName = dataDir + NETWORK_SYNC_SINGLE_VIRTUAL_SOURCE_METHOD_Y_FILE;
	project_.loadHDF5(yFileName, NETWORK_SYNC_SINGLE_VIRTUAL_SOURCE_METHOD_Y_DATASET, yList);

	Timer timer;

	for (unsigned int acqNumber = 0; acqNumber < yList.size(); ++acqNumber) {
		std::string acqDataDir = FileUtil::path(dataDir, "/", acqNumber);
		if (!project_.directoryExists(acqDataDir)) {
			break;
		}

		LOG_INFO << "ACQ " << acqNumber;

		for (auto iter = gridData.begin(); iter != gridData.end(); ++iter) {
			iter->y = focusY; // this is the y value for processing
		}

		acquisition->setDataDir(acqDataDir);
		processor->process(baseElement, gridData);

		for (auto iter = gridData.begin(); iter != gridData.end(); ++iter) {
			iter->y = yList[acqNumber]; // this is the y value that will be saved to file
		}

		std::string acqOutputDir = FileUtil::path(outputDir, "/", acqNumber);
		project_.createDirectory(acqOutputDir, true);

		if (config.valueScale != 0.0) {
			std::for_each(gridData.begin(), gridData.end(),
				Util::MultiplyValueBy<XYZValueFactor<FloatType>, FloatType>(config.valueScale));
		}

		project_.saveImageToHDF5(gridData, acqOutputDir);
		const FloatType maxAbsValue = Util::maxAbsoluteValueField<XYZValueFactor<FloatType>, FloatType>(gridData);
		if (maxAbsValue > valueLevel) valueLevel = maxAbsValue;

		project_.showFigure3D(1, "Raw image", &gridData, &pointList,
					true, Figure::VISUALIZATION_RECTIFIED_LOG, Figure::COLORMAP_VIRIDIS,
					config.valueScale != 0.0 ? 1.0 : 0.0);

		if (coherenceFactorEnabled) {
			LOG_DEBUG << "Saving the image factors...";
			project_.saveHDF5(gridData, acqOutputDir + "/image_factor", "factor", Util::CopyFactorOp());

			// Apply the coherence factor method.
			for (auto iter = gridData.begin(); iter != gridData.end(); ++iter) {
				iter->value *= iter->factor;
				iter->factor = 1.0;
			}
			const FloatType maxAbsCFValue = Util::maxAbsoluteValueField<XYZValueFactor<FloatType>, FloatType>(gridData);
			if (maxAbsCFValue > cfValueLevel) cfValueLevel = maxAbsCFValue;

			LOG_DEBUG << "Saving the CF image...";
			project_.saveHDF5(gridData, acqOutputDir + "/image_cf", "cf", Util::CopyValueOp());

			project_.showFigure3D(2, "Coherence factor image", &gridData, &pointList,
						true, Figure::VISUALIZATION_RECTIFIED_LOG, Figure::COLORMAP_VIRIDIS,
						config.valueScale != 0.0 ? 1.0 : 0.0);
		}
	}

	LOG_INFO << "=== value level: " << valueLevel << " (" << Util::linearToDecibels(valueLevel) << " dB)";
	if (coherenceFactorEnabled) {
		LOG_INFO << "=== CF value level: " << cfValueLevel << " (" << Util::linearToDecibels(cfValueLevel) << " dB)";
	}

	LOG_DEBUG << ">>> Acquisition + processing + saving data time: " << timer.getTime();
}

template<typename FloatType>
void
NetworkSyncSingleVirtualSourceMethod<FloatType>::saveSignals(ConstParameterMapPtr taskPM, TnRnAcquisition<FloatType>& acq,
								unsigned int baseElement, const std::vector<FloatType>& txDelays,
								const std::string& dataDir)
{
	ConstParameterMapPtr scanPM = project_.loadChildParameterMap(taskPM, "scan_config_file");
	const unsigned int serverPort = scanPM->value<unsigned int>("sync_server_port", 1024, 65535);
	const bool asyncAcq           = scanPM->value<bool>(        "async_acquisition");
	const FloatType minY          = scanPM->value<FloatType>(   "min_y", -10000.0, 10000.0);

	SyncServer server(serverPort);
	std::vector<typename TnRnAcquisition<FloatType>::AcquisitionDataType> acqDataList;
	acqDataList.reserve(1000);
	std::vector<double> timeList;
	timeList.reserve(1000);
	std::vector<double> yList;
	yList.reserve(1000);

	// Capture signals.
	Timer timer;
	unsigned int acqNumber = 0;
	if (asyncAcq) {
		const FloatType ySpeed = scanPM->value<FloatType>("y_speed", 1.0, 100000.0) * (1.0e-3 / 60.0);
		const FloatType maxY   = scanPM->value<FloatType>("max_y", minY, 10000.0);
		LOG_DEBUG << "Waiting for trigger...";
		if (server.waitForTrigger()) {
			server.freeTrigger(); // start the array movement
			const double t0 = timer.getTime();
			FloatType y = minY;
			while (y <= maxY) {
				if (acqNumber == NETWORK_SYNC_SINGLE_VIRTUAL_SOURCE_METHOD_MAX_STEPS) break;
				LOG_INFO << "ACQ " << acqNumber;

				acqDataList.emplace_back();
				timeList.push_back(timer.getTime() - t0);
				y = minY + timeList.back() * ySpeed; // variable step
				yList.push_back(y);
				acq.execute(baseElement, txDelays, acqDataList.back());

				++acqNumber;
			}
		}
		LOG_DEBUG << "Average y step: " << (yList.back() - yList.front()) / (yList.size() - 1);
	} else {
		const FloatType yStep = scanPM->value<FloatType>("y_step", 1.0e-6, 1000.0);
		const double t0 = timer.getTime();
		while (true) {
			LOG_DEBUG << "Waiting for trigger...";
			if (!server.waitForTrigger()) break; // trigger abort
			if (acqNumber == NETWORK_SYNC_SINGLE_VIRTUAL_SOURCE_METHOD_MAX_STEPS) break;
			LOG_INFO << "ACQ " << acqNumber;

			acqDataList.emplace_back();
			timeList.push_back(timer.getTime() - t0);
			yList.push_back(minY + acqNumber * yStep); // fixed step
			acq.execute(baseElement, txDelays, acqDataList.back());

			server.freeTrigger(); // let the array move to the next position

			++acqNumber;
		}
	}

	// Save signals.
	for (unsigned int i = 0, iEnd = acqDataList.size(); i < iEnd; ++i) {
		project_.saveSignalsToHDF5(acqDataList[i], dataDir, i, baseElement);
	}
	// Save times.
	const std::string timeFileName = dataDir + NETWORK_SYNC_SINGLE_VIRTUAL_SOURCE_METHOD_TIME_FILE;
	project_.saveHDF5(timeList, timeFileName, NETWORK_SYNC_SINGLE_VIRTUAL_SOURCE_METHOD_TIME_DATASET);
	// Save y.
	const std::string yFileName = dataDir + NETWORK_SYNC_SINGLE_VIRTUAL_SOURCE_METHOD_Y_FILE;
	project_.saveHDF5(yList, yFileName, NETWORK_SYNC_SINGLE_VIRTUAL_SOURCE_METHOD_Y_DATASET);
}

} // namespace Lab

#endif // NETWORK_SYNC_SINGLE_VIRTUAL_SOURCE_METHOD_H
