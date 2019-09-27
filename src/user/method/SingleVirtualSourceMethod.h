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
#include "Exception.h"
#include "FileUtil.h"
#include "global.h"
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
#include "XY.h"
#include "XYZ.h"
#include "XYZValueFactor.h"

#define SINGLE_VIRTUAL_SOURCE_METHOD_TIME_FILE "/time"
#define SINGLE_VIRTUAL_SOURCE_METHOD_TIME_DATASET "time"



namespace Lab {

template<typename FloatType>
class SingleVirtualSourceMethod : public Method {
public:
	SingleVirtualSourceMethod(Project& project);
	virtual ~SingleVirtualSourceMethod();

	virtual void execute();

private:
	SingleVirtualSourceMethod(const SingleVirtualSourceMethod&) = delete;
	SingleVirtualSourceMethod& operator=(const SingleVirtualSourceMethod&) = delete;

	void process(FloatType valueScale, ArrayProcessor<FloatType>& processor, unsigned int baseElement,
			bool saveCoordinates, const std::string& outputDir);
	void applyCoherenceFactor();
	void useCoherenceFactor(FloatType valueScale, const std::string& outputDir);
	void execContinuousNetworkImaging(FloatType valueScale, ArrayProcessor<FloatType>& processor,
						unsigned int baseElement, bool coherenceFactorEnabled);
	void saveSignalSequence(ConstParameterMapPtr taskPM, unsigned int baseElement,
				const std::vector<FloatType>& txDelays,
				TnRnAcquisition<FloatType>& acquisition);
	void createImagesFromSavedSignalSequence(ConstParameterMapPtr taskPM,
							unsigned int baseElement, FloatType valueScale, bool coherenceFactorEnabled,
							TnRnAcquisition<FloatType>& acq, ArrayProcessor<FloatType>& processor);

	Project& project_;
	Matrix<XYZValueFactor<FloatType>> gridData_;
	std::vector<XYZ<float>> pointList_;
};



template<typename FloatType>
SingleVirtualSourceMethod<FloatType>::SingleVirtualSourceMethod(Project& project)
		: project_(project)
		, pointList_{{0.0, 0.0, 0.0}}
{
}

template<typename FloatType>
SingleVirtualSourceMethod<FloatType>::~SingleVirtualSourceMethod()
{
}

template<typename FloatType>
void
SingleVirtualSourceMethod<FloatType>::applyCoherenceFactor()
{
	for (auto iter = gridData_.begin(); iter != gridData_.end(); ++iter) {
		iter->value *= iter->factor;
		//iter->factor = 1.0;
	}
}

template<typename FloatType>
void
SingleVirtualSourceMethod<FloatType>::useCoherenceFactor(FloatType valueScale, const std::string& outputDir)
{
	LOG_DEBUG << "Saving the image factors...";
	project_.saveHDF5(gridData_, outputDir + "/image_factor", "factor", Util::CopyFactorOp());

	applyCoherenceFactor();

	LOG_DEBUG << "Saving the CF image...";
	project_.saveHDF5(gridData_, outputDir + "/image_cf", "cf", Util::CopyValueOp());

	project_.showFigure3D(2, "Coherence factor image", &gridData_, &pointList_,
				true, Figure::VISUALIZATION_RECTIFIED_LOG, Figure::COLORMAP_VIRIDIS, valueScale);
}

template<typename FloatType>
void
SingleVirtualSourceMethod<FloatType>::process(FloatType valueScale, ArrayProcessor<FloatType>& processor,
						unsigned int baseElement, bool saveCoordinates, const std::string& outputDir)
{
	Timer tProc;

	processor.prepare(baseElement);
	processor.process(gridData_);

	project_.saveImageToHDF5(gridData_, outputDir, saveCoordinates);

	project_.showFigure3D(1, "Raw image", &gridData_, &pointList_,
				true, Figure::VISUALIZATION_RECTIFIED_LOG, Figure::COLORMAP_VIRIDIS, valueScale);

	LOG_DEBUG << ">>> Acquisition + processing time: " << tProc.getTime();
}

template<typename FloatType>
void
SingleVirtualSourceMethod<FloatType>::execContinuousNetworkImaging(FloatType valueScale, ArrayProcessor<FloatType>& processor,
									unsigned int baseElement, bool coherenceFactorEnabled)
{
	int n = 0;
	Timer t;
	processor.prepare(baseElement);
	do {
		processor.process(gridData_);

		if (coherenceFactorEnabled) {
			applyCoherenceFactor();
		}

		project_.showFigure3D(1, "Image", &gridData_, &pointList_,
					true, Figure::VISUALIZATION_RECTIFIED_LOG, Figure::COLORMAP_VIRIDIS, valueScale);

		if (++n == 10) {
			LOG_INFO << 10.0 / t.getTime() << " image/s";
			n = 0;
			t.reset();
		}
	} while (!project_.processingCancellationRequested());
}

template<typename FloatType>
void
SingleVirtualSourceMethod<FloatType>::saveSignalSequence(ConstParameterMapPtr taskPM, unsigned int baseElement,
								const std::vector<FloatType>& txDelays,
								TnRnAcquisition<FloatType>& acquisition)
{
	const std::string dataDir = taskPM->value<std::string>("data_dir");
	ConstParameterMapPtr seqPM = project_.loadChildParameterMap(taskPM, "seq_config_file");
	const double acqTime = seqPM->value<double>("acquisition_time", 1.0, 60.0);

	std::vector<typename TnRnAcquisition<FloatType>::AcquisitionDataType> acqDataList;
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
	const std::string fileName = dataDir + SINGLE_VIRTUAL_SOURCE_METHOD_TIME_FILE;
	project_.saveHDF5(timeList, fileName, SINGLE_VIRTUAL_SOURCE_METHOD_TIME_DATASET);
}

template<typename FloatType>
void
SingleVirtualSourceMethod<FloatType>::createImagesFromSavedSignalSequence(ConstParameterMapPtr taskPM,
							unsigned int baseElement, FloatType valueScale, bool coherenceFactorEnabled,
							TnRnAcquisition<FloatType>& acq, ArrayProcessor<FloatType>& processor)
{
	const std::string dataDir   = taskPM->value<std::string>("data_dir");
	const std::string outputDir = taskPM->value<std::string>("output_dir");
	auto& savedAcq = dynamic_cast<SavedTnRnAcquisition<FloatType>&>(acq);

	// Load times.
	std::vector<double> timeList;
	project_.loadHDF5(dataDir + SINGLE_VIRTUAL_SOURCE_METHOD_TIME_FILE,
				SINGLE_VIRTUAL_SOURCE_METHOD_TIME_DATASET, timeList);

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
	project_.saveHDF5(timeList, outputDir + SINGLE_VIRTUAL_SOURCE_METHOD_TIME_FILE,
				SINGLE_VIRTUAL_SOURCE_METHOD_TIME_DATASET);
}

template<typename FloatType>
void
SingleVirtualSourceMethod<FloatType>::execute()
{
	ConstParameterMapPtr taskPM = project_.taskParameterMap();
	ConstParameterMapPtr svsPM   = project_.loadChildParameterMap(taskPM, "svs_config_file");
	ConstParameterMapPtr arrayPM = project_.loadChildParameterMap(taskPM, "array_config_file");
	const TnRnConfiguration<FloatType> config(svsPM, arrayPM);
	const unsigned int baseElement = svsPM->value<unsigned int>("base_element", 0, config.numElementsMux - 1U);
	const FloatType focusZ         = svsPM->value<FloatType>(   "tx_focus_z", -10000.0, 10000.0);

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

	std::unique_ptr<TnRnAcquisition<FloatType>> acquisition;
	std::string savedDataDir;

	switch (project_.method()) {
	case MethodEnum::single_virtual_source_3d_simulated_save_signals: // falls through
	case MethodEnum::single_virtual_source_3d_vectorial_simulated:
		acquisition = std::make_unique<Simulated3DTnRnAcquisition<FloatType>>(project_, config);
		break;
	case MethodEnum::single_virtual_source_3d_network_save_signals:            // falls through
	case MethodEnum::single_virtual_source_3d_network_save_signal_sequence:    // falls through
	case MethodEnum::single_virtual_source_3d_vectorial_dp_network:            // falls through
	case MethodEnum::single_virtual_source_3d_vectorial_sp_network_continuous:
		acquisition = std::make_unique<NetworkTnRnAcquisition<FloatType>>(project_, config);
		break;
	case MethodEnum::single_virtual_source_3d_vectorial_dp_saved:
		savedDataDir = FileUtil::path(taskPM->value<std::string>("data_dir"), "/", 0); // falls through
	case MethodEnum::single_virtual_source_3d_vectorial_dp_saved_sequence:
		acquisition = std::make_unique<SavedTnRnAcquisition<FloatType>>(project_, config.numElements, savedDataDir);
		break;
	default:
		THROW_EXCEPTION(InvalidParameterException, "Invalid method: " << static_cast<int>(project_.method()) << '.');
	}

	if (project_.method() == MethodEnum::single_virtual_source_3d_simulated_save_signals ||
			project_.method() == MethodEnum::single_virtual_source_3d_network_save_signals) {
		const std::string dataDir = taskPM->value<std::string>("data_dir");
		typename TnRnAcquisition<FloatType>::AcquisitionDataType acqData;
		acquisition->prepare(baseElement, txDelays);
		acquisition->execute(acqData);
		project_.saveSignalsToHDF5(acqData, dataDir, 0, baseElement);
		return;
	} else if (project_.method() == MethodEnum::single_virtual_source_3d_network_save_signal_sequence) {
		saveSignalSequence(taskPM, baseElement, txDelays, *acquisition);
		return;
	}

	const FloatType nyquistRate = 2.0 * config.maxFrequency;
	const FloatType nyquistLambda = config.propagationSpeed / nyquistRate;
	ImageGrid<FloatType>::get(project_.loadChildParameterMap(taskPM, "grid_config_file"), nyquistLambda, gridData_);

	if (project_.method() == MethodEnum::single_virtual_source_3d_vectorial_simulated ||
			project_.method() == MethodEnum::single_virtual_source_3d_vectorial_dp_network ||
			project_.method() == MethodEnum::single_virtual_source_3d_vectorial_dp_saved ||
			project_.method() == MethodEnum::single_virtual_source_3d_vectorial_dp_saved_sequence ||
			project_.method() == MethodEnum::single_virtual_source_3d_vectorial_sp_network_continuous) {

		ConstParameterMapPtr imagPM = project_.loadChildParameterMap(taskPM, "imag_config_file");
		const FloatType peakOffset          = imagPM->value<FloatType>(   "peak_offset", 0.0, 50.0);
		const unsigned int upsamplingFactor = imagPM->value<unsigned int>("upsampling_factor", 1, 128);
		const std::string rxApodFile        = imagPM->value<std::string>( "rx_apodization_file");
		std::vector<FloatType> rxApod;
		project_.loadHDF5(rxApodFile, "apod", rxApod);

		AnalyticSignalCoherenceFactorProcessor<FloatType> coherenceFactor(project_.loadChildParameterMap(taskPM, "coherence_factor_config_file"));
		auto processor = std::make_unique<Vectorial3DTnRnProcessor<FloatType>>(
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
			const std::string outputDir = taskPM->value<std::string>("output_dir");
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
