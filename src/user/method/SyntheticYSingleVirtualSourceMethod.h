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

#ifndef SYNTHETIC_Y_SINGLE_VIRTUAL_SOURCE_METHOD_H
#define SYNTHETIC_Y_SINGLE_VIRTUAL_SOURCE_METHOD_H

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
#include "Project.h"
#include "SavedTnRnAcquisition.h"
#include "Timer.h"
#include "TnRnConfiguration.h"
#include "SynthYVectorial3DTnRnProcessor.h"
#include "Util.h"
#include "XYZ.h"
#include "XYZValueFactor.h"



namespace Lab {

template<typename FloatType>
class SyntheticYSingleVirtualSourceMethod : public Method {
public:
	SyntheticYSingleVirtualSourceMethod(Project& project);
	virtual ~SyntheticYSingleVirtualSourceMethod();

	virtual void execute();

private:
	static constexpr const char* yFile    = "/y";
	static constexpr const char* yDataset = "y";

	SyntheticYSingleVirtualSourceMethod(const SyntheticYSingleVirtualSourceMethod&) = delete;
	SyntheticYSingleVirtualSourceMethod& operator=(const SyntheticYSingleVirtualSourceMethod&) = delete;

	Project& project_;
};



template<typename FloatType>
SyntheticYSingleVirtualSourceMethod<FloatType>::SyntheticYSingleVirtualSourceMethod(Project& project)
		: project_(project)
{
}

template<typename FloatType>
SyntheticYSingleVirtualSourceMethod<FloatType>::~SyntheticYSingleVirtualSourceMethod()
{
}

template<typename FloatType>
void
SyntheticYSingleVirtualSourceMethod<FloatType>::execute()
{
	ConstParameterMapPtr taskPM = project_.taskParameterMap();
	ConstParameterMapPtr svsPM   = project_.loadChildParameterMap(taskPM, "svs_config_file");
	ConstParameterMapPtr arrayPM = project_.loadChildParameterMap(taskPM, "array_config_file");
	const TnRnConfiguration<FloatType> config(svsPM, arrayPM);
	const unsigned int baseElement = svsPM->value<unsigned int>("base_element", 0, config.numElementsMux - config.numElements);
	const FloatType focusZ         = svsPM->value<FloatType>(   "tx_focus_z", -10000.0, 10000.0);
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
	auto processor = std::make_unique<SynthYVectorial3DTnRnProcessor<FloatType>>(config, *acquisition,
					upsamplingFactor, coherenceFactor, peakOffset, rxApod);
	processor->setTxDelays(focusX, focusY, focusZ, txDelays);
	processor->prepare(baseElement);

	std::vector<XYZ<float>> pointList = {{0.0, 0.0, 0.0}};
	FloatType valueLevel = 0.0;
	FloatType cfValueLevel = 0.0;

	// Load y (acquisition position).
	std::vector<double> yList;
	const std::string yFileName = dataDir + yFile;
	project_.loadHDF5(yFileName, yDataset, yList);

	const unsigned int synthYSize = imagPM->value<unsigned int>("synthetic_y_array_size", 1, yList.size());

	Timer timer;

	for (unsigned int acqNumber = 0; acqNumber < yList.size(); ++acqNumber) {
		std::string acqDataDir = FileUtil::path(dataDir, "/", acqNumber);
		if (!project_.directoryExists(acqDataDir)) {
			break;
		}

		LOG_INFO << "ACQ " << acqNumber;

		acquisition->setDataDir(acqDataDir);
		processor->getAcqData(yList[acqNumber]);
	}

	for (unsigned int acqNumber = 0; acqNumber <= yList.size() - synthYSize; ++acqNumber) {
		LOG_INFO << "ACQ " << acqNumber << '/' << yList.size() - synthYSize;

		for (auto iter = gridData.begin(); iter != gridData.end(); ++iter) {
			iter->y = focusY; // this is the y value for processing
		}

		processor->process(acqNumber, acqNumber + synthYSize - 1U, gridData);

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

} // namespace Lab

#endif // SYNTHETIC_Y_SINGLE_VIRTUAL_SOURCE_METHOD_H
