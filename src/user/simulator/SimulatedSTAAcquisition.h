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

#ifndef SIMULATEDSTAACQUISITION_H_
#define SIMULATEDSTAACQUISITION_H_

#include <memory>
#include <string>
#include <vector>

#include "Exception.h"
#include "Log.h"
#include "Matrix.h"
#include "ParameterMap.h"
#include "Project.h"
#include "Simulated3DAcquisitionDevice.h"
#include "STAAcquisition.h"
#include "STAConfiguration.h"
#include "Util.h"
#include "XYZValue.h"



namespace Lab {

template<typename TFloat>
class SimulatedSTAAcquisition : public STAAcquisition<TFloat> {
public:
	SimulatedSTAAcquisition(Project& project, const STAConfiguration<TFloat>& config);
	virtual ~SimulatedSTAAcquisition() = default;

	virtual void prepare(unsigned int baseElement);
	virtual void execute(unsigned int txElement,
				typename STAAcquisition<TFloat>::AcquisitionDataType& acqData);
private:
	SimulatedSTAAcquisition(const SimulatedSTAAcquisition&) = delete;
	SimulatedSTAAcquisition& operator=(const SimulatedSTAAcquisition&) = delete;
	SimulatedSTAAcquisition(SimulatedSTAAcquisition&&) = delete;
	SimulatedSTAAcquisition& operator=(SimulatedSTAAcquisition&&) = delete;

	Project& project_;
	const STAConfiguration<TFloat>& config_;
	TFloat maxAbsValue_; // auxiliar
	std::unique_ptr<Simulated3DAcquisitionDevice<TFloat>> acqDevice_;
	std::vector<XYZValue<TFloat>> reflectorList_;
	unsigned int baseElement_;
};



template<typename TFloat>
SimulatedSTAAcquisition<TFloat>::SimulatedSTAAcquisition(Project& project, const STAConfiguration<TFloat>& config)
		: project_(project)
		, config_(config)
		, maxAbsValue_()
		, baseElement_()
{
	const ParamMapPtr pm      = project_.getSubParamMap("simulated_3d_acquisition_config_file");
	const ParamMapPtr arrayPM = project_.getSubParamMap("array_config_file");

	const auto reflectorsFileName = pm->value<std::string>("reflectors_file");
	const auto reflectorsOffsetX  = pm->value<TFloat>(     "reflectors_offset_x", -10000.0, 10000.0);

	Matrix<TFloat> inputReflectorList;
	project_.loadHDF5(reflectorsFileName, "reflectors", inputReflectorList);
	if (inputReflectorList.n2() != 2) {
		THROW_EXCEPTION(InvalidValueException, "Wrong number of columns (" << inputReflectorList.n2() <<
				") in the file \"" << reflectorsFileName << ".h5\" (should be 2). ");
	}

	reflectorList_.resize(inputReflectorList.n1());
	for (std::size_t i = 0, end = inputReflectorList.n1(); i < end; ++i) {
		XYZValue<TFloat>& data = reflectorList_[i];
		data.x     = inputReflectorList(i, 0);
		data.y     = 0.0;
		data.z     = inputReflectorList(i, 1);
		data.value = 1.0;
	}

	acqDevice_ = std::make_unique<Simulated3DAcquisitionDevice<TFloat>>(
								*pm,
								*arrayPM,
								config_.samplingFrequency,
								config_.propagationSpeed,
								config_.maxFrequency,
								project_.expDirectory());
	acqDevice_->setAcquisitionTime(config_.acquisitionTime);
	acqDevice_->setExcitationWaveform(config_.centerFrequency);
	acqDevice_->setReflectorList(reflectorList_);
	acqDevice_->setReflectorOffset(reflectorsOffsetX, 0.0);
	acqDevice_->setGain(config_.minGain);
}

template<typename TFloat>
void
SimulatedSTAAcquisition<TFloat>::prepare(unsigned int baseElement)
{
	if (baseElement + config_.numElements > config_.numElementsMux) {
		THROW_EXCEPTION(InvalidParameterException, "Invalid base element:" << baseElement << '.');
	}

	std::vector<bool> rxMask(config_.numElementsMux);
	for (unsigned int i = baseElement, end = baseElement + config_.numElements; i < end; ++i) {
		rxMask[i] = true;
	}
	acqDevice_->setActiveRxElements(rxMask);

	baseElement_ = baseElement;
}

template<typename TFloat>
void
SimulatedSTAAcquisition<TFloat>::execute(unsigned int txElement,
						typename STAAcquisition<TFloat>::AcquisitionDataType& acqData)
{
	LOG_DEBUG << "ACQ baseElement=" << baseElement_ << " txElement=" << txElement;

	const std::size_t signalLength = acqDevice_->signalLength();
	if (signalLength == 0) {
		THROW_EXCEPTION(InvalidValueException, "signalLength = 0.");
	}
	if (baseElement_ + txElement >= config_.numElementsMux) {
		THROW_EXCEPTION(InvalidValueException, "Invalid tx element: " << txElement << '.');
	}
	if (acqData.n1() != config_.numElements || acqData.n2() != signalLength) {
		acqData.resize(config_.numElements, signalLength);
		LOG_DEBUG << "RESIZE acqData(channels, signalLength): channels=" << config_.numElements
				<< " signalLength=" << signalLength;
	}

	std::vector<bool> txMask(config_.numElementsMux);
	txMask[baseElement_ + txElement] = true;
	acqDevice_->setActiveTxElements(txMask);

	const std::vector<TFloat>& signalList = acqDevice_->getSignalList();

	const TFloat mx = Util::maxAbsolute(signalList);
	if (mx > maxAbsValue_) maxAbsValue_ = mx;
	LOG_DEBUG << "########## max(abs(signalList)): " << mx << " global: " << maxAbsValue_;

	std::copy(signalList.begin(), signalList.end(), acqData.begin());
}

} // namespace Lab

#endif /* SIMULATEDSTAACQUISITION_H_ */
