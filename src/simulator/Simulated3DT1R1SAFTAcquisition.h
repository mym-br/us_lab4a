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
#ifndef SIMULATED3DT1R1SAFTACQUISITION_H
#define SIMULATED3DT1R1SAFTACQUISITION_H

#include <memory>
#include <string>
#include <vector>

#include "Exception.h"
#include "Log.h"
#include "Matrix.h"
#include "ParameterMap.h"
#include "Project.h"
#include "Simulated3DAcquisitionDevice.h"
#include "SA3DConfiguration.h"
#include "STAAcquisition.h"
#include "Util.h"
#include "XYZValue.h"



namespace Lab {

template<typename FloatType>
class Simulated3DT1R1SAFTAcquisition : public STAAcquisition<FloatType> {
public:
	Simulated3DT1R1SAFTAcquisition(Project& project, const SA3DConfiguration<FloatType>& config);
	virtual ~Simulated3DT1R1SAFTAcquisition();

	virtual void execute(unsigned int baseElement, unsigned int txElement,
				typename STAAcquisition<FloatType>::AcquisitionDataType& acqData);

	void modifyReflectorsOffset(FloatType offsetX, FloatType offsetY);
private:
	Simulated3DT1R1SAFTAcquisition(const Simulated3DT1R1SAFTAcquisition&) = delete;
	Simulated3DT1R1SAFTAcquisition& operator=(const Simulated3DT1R1SAFTAcquisition&) = delete;

	Project& project_;
	const SA3DConfiguration<FloatType>& config_;
	FloatType maxAbsValue_; // auxiliar
	std::unique_ptr<Simulated3DAcquisitionDevice<FloatType>> acqDevice_;
	std::vector<XYZValue<FloatType>> reflectorList_;
	FloatType reflectorsOffsetX_;
	FloatType reflectorsOffsetY_;
};



template<typename FloatType>
Simulated3DT1R1SAFTAcquisition<FloatType>::Simulated3DT1R1SAFTAcquisition(Project& project, const SA3DConfiguration<FloatType>& config)
		: project_{project}
		, config_{config}
		, maxAbsValue_{}
{
	//TODO: check numChannels/numChannelsMux and other params

	ConstParameterMapPtr taskPM  = project_.taskParameterMap();
	ConstParameterMapPtr pm      = project_.loadChildParameterMap(taskPM, "simulated_3d_t1r1saft_acquisition_config_file");
	ConstParameterMapPtr arrayPM = project_.loadChildParameterMap(taskPM, "array_config_file");

	const std::string reflectorsFileName = pm->value<std::string>("reflectors_file");
	reflectorsOffsetX_                   = pm->value<FloatType>(  "reflectors_offset_x", -10000.0, 10000.0);
	reflectorsOffsetY_                   = pm->value<FloatType>(  "reflectors_offset_y", -10000.0, 10000.0);

	Matrix<FloatType> inputReflectorList;
	project_.loadHDF5(reflectorsFileName, "reflectors", inputReflectorList);
	if (inputReflectorList.n2() != 4) {
		THROW_EXCEPTION(InvalidValueException, "Wrong number of columns (" << inputReflectorList.n2() <<
				") in the file \"" << reflectorsFileName << ".h5\" (should be 4). ");
	}

	reflectorList_.resize(inputReflectorList.n1());
	for (std::size_t i = 0, end = inputReflectorList.n1(); i < end; ++i) {
		XYZValue<FloatType>& data = reflectorList_[i];
		data.x     = inputReflectorList(i, 0);
		data.y     = inputReflectorList(i, 1);
		data.z     = inputReflectorList(i, 2);
		data.value = inputReflectorList(i, 3);
	}

	acqDevice_ = std::make_unique<Simulated3DAcquisitionDevice<FloatType>>(
								*pm,
								*arrayPM,
								config_.samplingFrequency,
								config_.propagationSpeed,
								config_.maxFrequency);
	acqDevice_->setAcquisitionTime(config_.acquisitionTime);
	acqDevice_->setExcitationWaveform(config_.centerFrequency);
	acqDevice_->setReflectorList(reflectorList_);
	acqDevice_->setReflectorOffset(reflectorsOffsetX_, reflectorsOffsetY_);
	acqDevice_->setGain(config_.minGain);
}

template<typename FloatType>
Simulated3DT1R1SAFTAcquisition<FloatType>::~Simulated3DT1R1SAFTAcquisition()
{
}

template<typename FloatType>
void
Simulated3DT1R1SAFTAcquisition<FloatType>::execute(unsigned int baseElement, unsigned int txElement,
						typename STAAcquisition<FloatType>::AcquisitionDataType& acqData)
{
	LOG_DEBUG << "ACQ baseElement=" << baseElement << " txElement=" << txElement;

	const std::size_t signalLength = acqDevice_->signalLength();
	if (signalLength == 0) {
		THROW_EXCEPTION(InvalidValueException, "signalLength = 0.");
	}
	if (acqData.n1() != config_.activeRxElem.size() || acqData.n2() != signalLength) {
		acqData.resize(config_.activeRxElem.size(), signalLength);
		LOG_DEBUG << "RESIZE acqData(channels, signalLength): channels=" << config_.activeRxElem.size()
				<< " signalLength=" << signalLength;
	}

	std::vector<bool> txMask(config_.txElemPos.size());
	bool found = false;
	unsigned int iTxElem = 0;
	for (unsigned int localElem : config_.activeTxElem) {
		if (txElement == localElem) {
			found = true;
			const unsigned int elem = baseElement + localElem;
			if (elem >= txMask.size()) {
				THROW_EXCEPTION(InvalidValueException, "Invalid active tx element: " << elem << '.');
			}
			break;
		}
		++iTxElem;
	}
	if (!found) THROW_EXCEPTION(InvalidValueException, "Invalid tx element: " << baseElement + txElement << '.');
	txMask[baseElement + txElement] = true;
	acqDevice_->setActiveTxElements(txMask);

	std::vector<bool> rxMask(config_.rxElemPos.size());
	if (iTxElem >= config_.activeRxElem.size()) {
		THROW_EXCEPTION(InvalidValueException, "Receive element not found.");
	}
	const unsigned int rxElem = baseElement + config_.activeRxElem[iTxElem];
	if (rxElem >= rxMask.size()) {
		THROW_EXCEPTION(InvalidValueException, "Invalid active rx element: " << rxElem << '.');
	}
	rxMask[rxElem] = true;
	acqDevice_->setActiveRxElements(rxMask);

	const std::vector<FloatType>& signalList = acqDevice_->getSignalList();

	const FloatType mx = Util::maxAbsolute(signalList);
	if (mx > maxAbsValue_) maxAbsValue_ = mx;
	LOG_DEBUG << "########## max(abs(signalList)): " << mx << " global: " << maxAbsValue_;

	std::copy(signalList.begin(), signalList.end(), acqData.begin()); //TODO: memcpy? how to avoid this copy?
}

template<typename FloatType>
void
Simulated3DT1R1SAFTAcquisition<FloatType>::modifyReflectorsOffset(FloatType offsetX, FloatType offsetY)
{
	acqDevice_->setReflectorOffset(reflectorsOffsetX_ + offsetX, reflectorsOffsetY_ + offsetY);
}

} // namespace Lab

#endif // SIMULATED3DT1R1SAFTACQUISITION_H
