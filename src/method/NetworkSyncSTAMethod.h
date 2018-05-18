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

#include <iomanip>
#include <memory>
#include <sstream>
#include <string>

#include <chrono>
#include <thread>//TODO: Remove test
using namespace std::chrono_literals;

#include "CoherenceFactor.h"
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
#include "XZ.h"
#include "XZValueFactor.h"



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

	void saveSignals(const STAConfiguration<FloatType>& config, STAAcquisition<FloatType>& acq, unsigned int baseElement, const std::string& dataDir);

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

	const STAConfiguration<FloatType> config(project_.loadChildParameterMap(taskPM, "sta_config_file"));
	const unsigned int baseElement = taskPM->value<unsigned int>("base_element", 0, config.numElementsMux - config.numElements);

	std::unique_ptr<STAAcquisition<FloatType>> acquisition;

	switch (project_.method()) {
	case MethodType::sta_network_sync_save_signals:
		//acquisition = std::make_unique<NetworkSTAAcquisition<FloatType>>(project_, config);
		break;
	case MethodType::sta_network_sync:
		acquisition = std::make_unique<SavedSTAAcquisition<FloatType>>(project_, config.numElements);
		break;
	default:
		THROW_EXCEPTION(InvalidParameterException, "Invalid method: " << static_cast<int>(project_.method()) << '.');
	}

	if (project_.method() == MethodType::sta_network_sync_save_signals) {
		std::string dataDir = taskPM->value<std::string>("data_dir");
		project_.createDirectory(dataDir, true);
		saveSignals(config, *acquisition, baseElement, dataDir);
		return;
	}

	const FloatType peakOffset  = taskPM->value<FloatType>(  "peak_offset" , 0.0, 50.0);
	const std::string outputDir = taskPM->value<std::string>("output_dir");
	project_.createDirectory(outputDir, false);

	Matrix2<XZValueFactor<FloatType>> gridData;

	const FloatType nyquistRate = 2.0 * config.maxFrequency;
	const FloatType nyquistLambda = config.propagationSpeed / nyquistRate;
	ImageGrid<FloatType>::get(project_.loadChildParameterMap(taskPM, "grid_config_file"), nyquistLambda, gridData);

	bool coherenceFactorEnabled = false;
	bool vectorialProcessingWithEnvelope = false;

	vectorialProcessingWithEnvelope     = taskPM->value<bool>(        "calculate_envelope_in_processing");
	const unsigned int upsamplingFactor = taskPM->value<unsigned int>("upsampling_factor", 1, 128);
	AnalyticSignalCoherenceFactorProcessor<FloatType> coherenceFactor(project_.loadChildParameterMap(taskPM, "coherence_factor_config_file"));
	coherenceFactorEnabled = coherenceFactor.enabled();
	auto processor = std::make_unique<VectorialSTAProcessor<FloatType>>(config, *acquisition, upsamplingFactor, coherenceFactor, peakOffset, vectorialProcessingWithEnvelope);

	Timer tProc;
	processor->process(baseElement, gridData);
	LOG_DEBUG << ">>> Acquisition + processing time: " << tProc.getTime();

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
NetworkSyncSTAMethod<FloatType>::saveSignals(const STAConfiguration<FloatType>& config, STAAcquisition<FloatType>& acq, unsigned int baseElement, const std::string& dataDir)
{
	SyncServer server(5544);

	while (true) {
		if (!server.waitForTrigger()) break;
		std::this_thread::sleep_for(3s);
		server.freeTrigger();
	}

#if 0
	typename STAAcquisition<FloatType>::AcquisitionDataType acqData;

	for (unsigned int txElem = 0; txElem < config.numElements; ++txElem) {
		acq.execute(baseElement, txElem, acqData);

		std::ostringstream filePath;
		filePath << dataDir << std::setfill('0') << "/signal-base" << std::setw(4) << baseElement << "-tx" << std::setw(4) << txElem;
		std::string filePathStr = filePath.str();
		LOG_DEBUG << "Saving " << filePathStr << "...";
		project_.saveHDF5(acqData, filePathStr, "signal");
	}
#endif
}

} // namespace Lab

#endif // NETWORKSYNCSTAMETHOD_H
