/***************************************************************************
 *  Copyright 2014, 2017, 2018, 2019 Marcelo Y. Matuda                     *
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

#ifndef SINGLEVIRTUALSOURCEMETHOD_H
#define SINGLEVIRTUALSOURCEMETHOD_H

#include <memory>
#include <string>
#include <vector>

#include "ArrayUtil.h"
#include "CoherenceFactor.h"
#include "Colormap.h"
#include "ContainerDumper.h"
#include "Exception.h"
#include "FileUtil.h"
#include "ImageGrid.h"
#include "Log.h"
#include "Matrix.h"
#include "Method.h"
#include "NetworkTnRnAcquisition.h"
#include "Project.h"
#include "SavedTnRnAcquisition.h"
#include "Simulated3DTnRnAcquisition.h"
#include "Timer.h"
#include "TnRnAcquisition.h"
#include "TnRnConfiguration.h"
#include "Util.h"
#include "Vectorial3DTnRnProcessor.h"
#include "Visualization.h"
#include "XY.h"
#include "XYZ.h"
#include "XYZValueFactor.h"



namespace Lab {

template<typename TFloat>
class SingleVirtualSourceMethod : public Method {
public:
	SingleVirtualSourceMethod(Project& project);
	virtual ~SingleVirtualSourceMethod() = default;

	virtual void execute();
private:
	static constexpr const char* timeFile    = "/time";
	static constexpr const char* timeDataset = "time";

	SingleVirtualSourceMethod(const SingleVirtualSourceMethod&) = delete;
	SingleVirtualSourceMethod& operator=(const SingleVirtualSourceMethod&) = delete;
	SingleVirtualSourceMethod(SingleVirtualSourceMethod&&) = delete;
	SingleVirtualSourceMethod& operator=(SingleVirtualSourceMethod&&) = delete;

	void process(TFloat valueScale, ArrayProcessor<XYZValueFactor<TFloat>>& processor, unsigned int baseElement,
			bool saveCoordinates, const std::string& outputDir);
	void useCoherenceFactor(TFloat valueScale, const std::string& outputDir);
	void execContinuousNetworkImaging(TFloat valueScale, ArrayProcessor<XYZValueFactor<TFloat>>& processor,
						unsigned int baseElement, bool coherenceFactorEnabled);
	void saveSignalSequence(const ParameterMap& taskPM, unsigned int baseElement,
				const std::vector<TFloat>& txDelays,
				TnRnAcquisition<TFloat>& acquisition);
	void createImagesFromSavedSignalSequence(const ParameterMap& taskPM,
							unsigned int baseElement, TFloat valueScale, bool coherenceFactorEnabled,
							TnRnAcquisition<TFloat>& acq, ArrayProcessor<XYZValueFactor<TFloat>>& processor);

	Project& project_;
	Matrix<XYZValueFactor<TFloat>> gridData_;
	std::vector<XYZ<float>> pointList_;
};



template<typename TFloat>
SingleVirtualSourceMethod<TFloat>::SingleVirtualSourceMethod(Project& project)
		: project_(project)
		, pointList_{{0.0, 0.0, 0.0}}
{
}

template<typename TFloat>
void
SingleVirtualSourceMethod<TFloat>::useCoherenceFactor(TFloat valueScale, const std::string& outputDir)
{
	project_.saveFactorToHDF5(gridData_, outputDir, "image_factor", "factor");

	Util::applyFactorToValue(gridData_.begin(), gridData_.end());

	project_.saveImageToHDF5(gridData_, outputDir, "image_cf", "cf");

	project_.showFigure3D(2, "Coherence factor image", &gridData_, &pointList_,
				true, Visualization::Value::RECTIFIED_LOG, Colormap::Id::GRADIENT_VIRIDIS, valueScale);
}

template<typename TFloat>
void
SingleVirtualSourceMethod<TFloat>::process(TFloat valueScale, ArrayProcessor<XYZValueFactor<TFloat>>& processor,
						unsigned int baseElement, bool saveCoordinates, const std::string& outputDir)
{
	Timer tProc;

	processor.prepare(baseElement);
	processor.process(gridData_);

	project_.saveImageToHDF5(gridData_, outputDir);
	if (saveCoordinates) {
		project_.saveXYZToHDF5(gridData_, outputDir);
	};

	project_.showFigure3D(1, "Raw image", &gridData_, &pointList_,
				true, Visualization::Value::RECTIFIED_LOG, Colormap::Id::GRADIENT_VIRIDIS, valueScale);

	LOG_DEBUG << ">>> Acquisition + processing time: " << tProc.getTime();
}

template<typename TFloat>
void
SingleVirtualSourceMethod<TFloat>::execContinuousNetworkImaging(TFloat valueScale, ArrayProcessor<XYZValueFactor<TFloat>>& processor,
									unsigned int baseElement, bool coherenceFactorEnabled)
{
	int n = 0;
	Timer t;
	processor.prepare(baseElement);
	do {
		processor.process(gridData_);

		if (coherenceFactorEnabled) {
			Util::applyFactorToValue(gridData_.begin(), gridData_.end());
		}

		project_.showFigure3D(1, "Image", &gridData_, &pointList_,
					true, Visualization::Value::RECTIFIED_LOG, Colormap::Id::GRADIENT_VIRIDIS, valueScale);

		if (++n == 10) {
			LOG_INFO << 10.0 / t.getTime() << " image/s";
			n = 0;
			t.reset();
		}
	} while (!project_.processingCancellationRequested());
}

template<typename TFloat>
void
SingleVirtualSourceMethod<TFloat>::saveSignalSequence(const ParameterMap& taskPM, unsigned int baseElement,
								const std::vector<TFloat>& txDelays,
								TnRnAcquisition<TFloat>& acquisition)
{
	const auto dataDir = taskPM.value<std::string>("data_dir");
	const ParamMapPtr seqPM = project_.getSubParamMap("seq_config_file");
	const auto acqTime = seqPM->value<double>("acquisition_time", 1.0, 60.0);

	std::vector<typename TnRnAcquisition<TFloat>::AcquisitionDataType> acqDataList;
	std::vector<double> timeList;

	// Capture signals.
	unsigned int acqNumber = 0;
	acquisition.prepare(baseElement, txDelays);
	Timer timer;
	const double t0 = timer.getTime();
	double t = 0.0;
	do {
		LOG_INFO << "ACQ " << acqNumber;

		timeList.push_back(t);
		acqDataList.emplace_back();
		acquisition.execute(acqDataList.back());

		t = timer.getTime() - t0;
		++acqNumber;
	} while (t <= acqTime);

	// Save signals.
	for (unsigned int i = 0, iEnd = acqDataList.size(); i < iEnd; ++i) {
		project_.saveSignalsToHDF5(acqDataList[i], dataDir, i, baseElement);
	}
	// Save times.
	const std::string fileName = dataDir + timeFile;
	project_.saveHDF5(timeList, fileName, timeDataset);
}

template<typename TFloat>
void
SingleVirtualSourceMethod<TFloat>::createImagesFromSavedSignalSequence(const ParameterMap& taskPM,
							unsigned int baseElement, TFloat valueScale, bool coherenceFactorEnabled,
							TnRnAcquisition<TFloat>& acq, ArrayProcessor<XYZValueFactor<TFloat>>& processor)
{
	const auto dataDir   = taskPM.value<std::string>("data_dir");
	const auto outputDir = taskPM.value<std::string>("output_dir");
	auto& savedAcq = dynamic_cast<SavedTnRnAcquisition<TFloat>&>(acq);

	// Load times.
	std::vector<double> timeList;
	project_.loadHDF5(dataDir + timeFile, timeDataset, timeList);

	for (unsigned int acqNumber = 0, end = timeList.size(); acqNumber < end; ++acqNumber) {
		LOG_INFO << "ACQ " << acqNumber;

		const std::string imgDir = FileUtil::path(outputDir, "/", acqNumber);
		project_.createDirectory(imgDir, false);

		savedAcq.setDataDir(FileUtil::path(dataDir, "/", acqNumber));

		// Process and save images.
		process(valueScale, processor, baseElement, acqNumber == 0, imgDir);
		if (coherenceFactorEnabled) {
			useCoherenceFactor(valueScale, imgDir);
		}
	}

	// Save times.
	project_.saveHDF5(timeList, outputDir + timeFile, timeDataset);
}

template<typename TFloat>
void
SingleVirtualSourceMethod<TFloat>::execute()
{
	const ParameterMap& taskPM = project_.taskParamMap();
	const ParamMapPtr svsPM   = project_.getSubParamMap("svs_config_file");
	const ParamMapPtr arrayPM = project_.getSubParamMap("array_config_file");
	const TnRnConfiguration<TFloat> config(*svsPM, *arrayPM);
	const auto baseElement = svsPM->value<unsigned int>("base_element", 0, config.numElementsMux - 1U);
	const auto focusZ      = svsPM->value<TFloat>(      "tx_focus_z", -10000.0, 10000.0);

	TFloat focusX = 0, focusY = 0;
	// Set the focus at the mean x, y.
	for (unsigned int i = baseElement, end = baseElement + config.numElements; i < end; ++i) {
		if (i >= config.txElemPos.size()) {
			THROW_EXCEPTION(InvalidValueException, "Invalid tx element number: " << i << '.');
		}
		const XY<TFloat>& pos = config.txElemPos[i];
		focusX += pos.x;
		focusY += pos.y;
	}
	focusX /= config.numElements;
	focusY /= config.numElements;
	std::vector<TFloat> txDelays;
	ArrayUtil::calculateTx3DFocusDelay(focusX, focusY, focusZ, config.propagationSpeed,
						config.txElemPos, baseElement, config.numElements, txDelays);
	ContainerDumper::save((project_.expDirectory() + "/tx_delays.txt").c_str(), txDelays.begin(), txDelays.end());

	std::unique_ptr<TnRnAcquisition<TFloat>> acquisition;
	std::string savedDataDir;

	switch (project_.method()) {
	case MethodEnum::single_virtual_source_3d_simulated_save_signals:
	case MethodEnum::single_virtual_source_3d_vectorial_simulated:
		acquisition = std::make_unique<Simulated3DTnRnAcquisition<TFloat>>(project_, config);
		break;
	case MethodEnum::single_virtual_source_3d_network_save_signals:
	case MethodEnum::single_virtual_source_3d_network_save_signal_sequence:
	case MethodEnum::single_virtual_source_3d_vectorial_dp_network:
	case MethodEnum::single_virtual_source_3d_vectorial_sp_network_continuous:
		acquisition = std::make_unique<NetworkTnRnAcquisition<TFloat>>(project_, config);
		break;
	case MethodEnum::single_virtual_source_3d_vectorial_dp_saved:
		savedDataDir = FileUtil::path(taskPM.value<std::string>("data_dir"), "/", 0);
		[[fallthrough]];
	case MethodEnum::single_virtual_source_3d_vectorial_dp_saved_sequence:
		acquisition = std::make_unique<SavedTnRnAcquisition<TFloat>>(project_, config.numElements, savedDataDir);
		break;
	default:
		THROW_EXCEPTION(InvalidParameterException, "Invalid method: " << static_cast<int>(project_.method()) << '.');
	}

	if (project_.method() == MethodEnum::single_virtual_source_3d_simulated_save_signals ||
			project_.method() == MethodEnum::single_virtual_source_3d_network_save_signals) {
		const auto dataDir = taskPM.value<std::string>("data_dir");
		typename TnRnAcquisition<TFloat>::AcquisitionDataType acqData;
		acquisition->prepare(baseElement, txDelays);
		acquisition->execute(acqData);
		project_.saveSignalsToHDF5(acqData, dataDir, 0, baseElement);
		return;
	} else if (project_.method() == MethodEnum::single_virtual_source_3d_network_save_signal_sequence) {
		saveSignalSequence(taskPM, baseElement, txDelays, *acquisition);
		return;
	}

	const TFloat nyquistLambda = Util::nyquistLambda(config.propagationSpeed, config.maxFrequency);
	ImageGrid<TFloat>::get(*project_.getSubParamMap("grid_config_file"), nyquistLambda, gridData_);

	if (project_.method() == MethodEnum::single_virtual_source_3d_vectorial_simulated ||
			project_.method() == MethodEnum::single_virtual_source_3d_vectorial_dp_network ||
			project_.method() == MethodEnum::single_virtual_source_3d_vectorial_dp_saved ||
			project_.method() == MethodEnum::single_virtual_source_3d_vectorial_dp_saved_sequence ||
			project_.method() == MethodEnum::single_virtual_source_3d_vectorial_sp_network_continuous) {

		const ParamMapPtr imagPM = project_.getSubParamMap("imag_config_file");
		const auto peakOffset       = imagPM->value<TFloat>(      "peak_offset", 0.0, 50.0);
		const auto upsamplingFactor = imagPM->value<unsigned int>("upsampling_factor", 1, 128);
		const auto rxApodFile       = imagPM->value<std::string>( "rx_apodization_file");
		std::vector<TFloat> rxApod;
		project_.loadHDF5(rxApodFile, "apod", rxApod);

		AnalyticSignalCoherenceFactorProcessor<TFloat> coherenceFactor(*project_.getSubParamMap("coherence_factor_config_file"));
		auto processor = std::make_unique<Vectorial3DTnRnProcessor<TFloat>>(
							config, *acquisition, upsamplingFactor,
							coherenceFactor, peakOffset,
							rxApod);
		processor->setTxDelays(focusX, focusY, focusZ, txDelays);
		if (project_.method() == MethodEnum::single_virtual_source_3d_vectorial_sp_network_continuous) {
			execContinuousNetworkImaging(config.valueScale, *processor, baseElement, coherenceFactor.enabled());
		} else if (project_.method() == MethodEnum::single_virtual_source_3d_vectorial_dp_saved_sequence) {
			createImagesFromSavedSignalSequence(taskPM, baseElement, config.valueScale, coherenceFactor.enabled(),
								*acquisition, *processor);
		} else {
			const auto outputDir = taskPM.value<std::string>("output_dir");
			project_.createDirectory(outputDir, false);
			process(config.valueScale, *processor, baseElement, true, outputDir);
			if (coherenceFactor.enabled()) {
				useCoherenceFactor(config.valueScale, outputDir);
			}
		}
	} else {
		THROW_EXCEPTION(InvalidParameterException, "Invalid method: " << static_cast<int>(project_.method()) << '.');
	}
}

} // namespace Lab

#endif // SINGLEVIRTUALSOURCEMETHOD_H
