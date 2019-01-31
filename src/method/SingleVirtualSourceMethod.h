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
#include "global.h"
#include "ImageGrid.h"
#include "Log.h"
#include "Matrix.h"
#include "Method.h"
#include "Project.h"
#include "Simulated3DTnRnAcquisition.h"
#include "Timer.h"
#include "TnRnAcquisition.h"
#include "TnRnConfiguration.h"
#include "Util.h"
#include "Vectorial3DTnRnProcessor.h"
#include "XY.h"
#include "XYZ.h"
#include "XYZValueFactor.h"



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

	void process(FloatType valueScale, ArrayProcessor<FloatType>& processor, unsigned int baseElement, const std::string& outputDir);
	void useCoherenceFactor(FloatType valueScale, const std::string& outputDir);

	Project& project_;
	Matrix<XYZValueFactor<FloatType>> gridData_;
	Project::GridDataType projGridData_;
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
SingleVirtualSourceMethod<FloatType>::useCoherenceFactor(FloatType valueScale, const std::string& outputDir)
{
	LOG_DEBUG << "Saving the image factors...";
	project_.saveHDF5(gridData_, outputDir + "/image_factor", "factor", Util::CopyFactorOp());

	// Applies the coherence factor method.
	for (auto iter = gridData_.begin(); iter != gridData_.end(); ++iter) {
		iter->value *= iter->factor;
		iter->factor = 1.0;
	}

	LOG_DEBUG << "Saving the CF image...";
	project_.saveHDF5(gridData_, outputDir + "/image_cf", "cf", Util::CopyValueOp());

	Util::copyXYZValue(gridData_, projGridData_);
	project_.showFigure3D(2, "Coherence factor image", &projGridData_, &pointList_,
				true, Figure::VISUALIZATION_RECTIFIED_LOG, Figure::COLORMAP_VIRIDIS, valueScale);
}

template<typename FloatType>
void
SingleVirtualSourceMethod<FloatType>::process(FloatType valueScale, ArrayProcessor<FloatType>& processor, unsigned int baseElement, const std::string& outputDir)
{
	Timer tProc;

	processor.process(baseElement, gridData_);

	project_.saveImageToHDF5(gridData_, outputDir);

	Util::copyXYZValue(gridData_, projGridData_);
	project_.showFigure3D(1, "Raw image", &projGridData_, &pointList_,
				true, Figure::VISUALIZATION_RECTIFIED_LOG, Figure::COLORMAP_VIRIDIS, valueScale);

	LOG_DEBUG << ">>> Acquisition + processing time: " << tProc.getTime();
}

template<typename FloatType>
void
SingleVirtualSourceMethod<FloatType>::execute()
{
	ConstParameterMapPtr taskPM = project_.taskParameterMap();
	ConstParameterMapPtr imgPM   = project_.loadChildParameterMap(taskPM, "img_config_file");
	ConstParameterMapPtr arrayPM = project_.loadChildParameterMap(taskPM, "array_config_file");
	const TnRnConfiguration<FloatType> config(imgPM, arrayPM);
	const unsigned int baseElement = taskPM->value<unsigned int>("base_element", 0, config.numElementsMux - 1U);
	const FloatType focusZ         = taskPM->value<FloatType>("tx_focus_z", -10000.0, 10000.0);

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
	ArrayUtil::calculateTx3DFocusDelay(focusX, focusY, focusZ, config.propagationSpeed, config.txElemPos, txDelays);

	std::unique_ptr<TnRnAcquisition<FloatType>> acquisition;

	switch (project_.method()) {
	case MethodType::single_virtual_source_3d_simulated_save_signals: // falls through
	case MethodType::single_virtual_source_3d_vectorial_simulated:
		acquisition = std::make_unique<Simulated3DTnRnAcquisition<FloatType>>(project_, config);
		break;
	default:
		THROW_EXCEPTION(InvalidParameterException, "Invalid method: " << static_cast<int>(project_.method()) << '.');
	}

	if (project_.method() == MethodType::single_virtual_source_3d_simulated_save_signals) {
		const std::string dataDir = taskPM->value<std::string>("data_dir");
		typename TnRnAcquisition<FloatType>::AcquisitionDataType acqData;
		acquisition->execute(baseElement, txDelays, acqData);
		project_.saveSignalsToHDF5(acqData, dataDir, 0, baseElement);
		return;
	}

	const FloatType peakOffset  = taskPM->value<FloatType>(  "peak_offset" , 0.0, 50.0);
	const std::string outputDir = taskPM->value<std::string>("output_dir");
	project_.createDirectory(outputDir, false);

	const std::string rxApodFile = taskPM->value<std::string>("rx_apodization_file");
	std::vector<FloatType> rxApod;
	project_.loadHDF5(rxApodFile, "apod", rxApod);

	const FloatType nyquistRate = 2.0 * config.maxFrequency;
	const FloatType nyquistLambda = config.propagationSpeed / nyquistRate;
	ImageGrid<FloatType>::get(project_.loadChildParameterMap(taskPM, "grid_config_file"), nyquistLambda, gridData_);

	if (project_.method() == MethodType::single_virtual_source_3d_vectorial_simulated) {
		const unsigned int upsamplingFactor = taskPM->value<unsigned int>("upsampling_factor", 1, 128);
		AnalyticSignalCoherenceFactorProcessor<FloatType> coherenceFactor(project_.loadChildParameterMap(taskPM, "coherence_factor_config_file"));
		auto processor = std::make_unique<Vectorial3DTnRnProcessor<FloatType>>(
							config, *acquisition, upsamplingFactor,
							coherenceFactor, peakOffset,
							rxApod);
		processor->setTxDelays(focusX, focusY, focusZ, txDelays);
		process(config.valueScale, *processor, baseElement, outputDir);
		if (coherenceFactor.enabled()) {
			useCoherenceFactor(config.valueScale, outputDir);
		}
	} else {
		THROW_EXCEPTION(InvalidParameterException, "Invalid method: " << static_cast<int>(project_.method()) << '.');
	}
}

} // namespace Lab

#endif // SINGLEVIRTUALSOURCEMETHOD_H
