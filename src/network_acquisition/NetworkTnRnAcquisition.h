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

#ifndef NETWORKTNRNACQUISITION_H
#define NETWORKTNRNACQUISITION_H

#include <algorithm> /* copy */
#include <cstddef> /* std::size_t */
#include <memory>
#include <string>
#include <vector>

#include "ArrayAcqClient.h"
#include "Exception.h"
#include "global.h"
#include "ParameterMap.h"
#include "Project.h"
#include "TnRnAcquisition.h"
#include "TnRnConfiguration.h"
#include "Util.h"



namespace Lab {

template<typename FloatType>
class NetworkTnRnAcquisition : public TnRnAcquisition<FloatType> {
public:
	NetworkTnRnAcquisition(
		const Project& project,
		const TnRnConfiguration<FloatType>& config);
	virtual ~NetworkTnRnAcquisition();

	virtual void execute(unsigned int baseElement, const std::vector<FloatType>& txDelays,
				typename TnRnAcquisition<FloatType>::AcquisitionDataType& acqData);
private:
	NetworkTnRnAcquisition(const NetworkTnRnAcquisition&) = delete;
	NetworkTnRnAcquisition& operator=(const NetworkTnRnAcquisition&) = delete;

	const Project& project_;
	const TnRnConfiguration<FloatType>& config_;
	std::unique_ptr<ArrayAcqClient> acq_;
	std::vector<float> signalBuffer_;
	double valueFactor_;
	unsigned int averageN_;
	std::vector<float> txDelays_;
};



template<typename FloatType>
NetworkTnRnAcquisition<FloatType>::NetworkTnRnAcquisition(const Project& project, const TnRnConfiguration<FloatType>& config)
		: project_(project)
		, config_(config)
		, acq_()
{
	ConstParameterMapPtr pm = project_.loadParameterMap(NETWORK_AQUISITION_CONFIG_FILE);
	const std::string serverIpAddress = pm->value<std::string>(   "server_ip_address");
	const unsigned short portNumber   = pm->value<unsigned short>("server_port_number",   49152,  65535);
	valueFactor_                = 1.0 / pm->value<double>(        "value_scale"       , 1.0e-30, 1.0e30);
	averageN_                         = pm->value<unsigned int>(  "average_n"         ,       1,    256);

	acq_ = std::make_unique<ArrayAcqClient>(serverIpAddress.c_str(), portNumber);

	acq_->execPreConfiguration();

	acq_->setSamplingFrequency(config_.samplingFrequency);
	acq_->setCenterFrequency(config_.centerFrequency, config_.numPulses);

	std::string txMask(config_.numElements, '1');
	acq_->setActiveTransmitElements(txMask);

	std::string rxMask(config_.numElements, '1');
	acq_->setActiveReceiveElements(rxMask);

	acq_->setGain(config_.minGain);
	acq_->setAcquisitionTime(config_.acquisitionTime);

	acq_->execPostConfiguration();
}

template<typename FloatType>
NetworkTnRnAcquisition<FloatType>::~NetworkTnRnAcquisition()
{
}

template<typename FloatType>
void
NetworkTnRnAcquisition<FloatType>::execute(unsigned int baseElement, const std::vector<FloatType>& txDelays,
						typename TnRnAcquisition<FloatType>::AcquisitionDataType& acqData)
{
	LOG_DEBUG << "ACQ baseElement=" << baseElement;
	if (baseElement + config_.numElements > config_.numElementsMux) {
		THROW_EXCEPTION(InvalidParameterException, "Invalid base element:" << baseElement << '.');
	}

	const std::size_t signalLength = acq_->getSignalLength();
	if (signalLength == 0) {
		THROW_EXCEPTION(InvalidValueException, "signalLength = 0.");
	}
	if (acqData.n1() != config_.numElements || acqData.n2() != signalLength) {
		acqData.resize(config_.numElements, signalLength);
		LOG_DEBUG << "RESIZE acqData(channels,signalLength): channels=" << config_.numElements << " signalLength=" << signalLength;
	}

	acq_->setBaseElement(baseElement);

	txDelays_.resize(txDelays.size());
	std::copy(txDelays.begin(), txDelays.end(), txDelays_.begin());
	acq_->setTransmitDelays(txDelays_);

	if (averageN_ > 1U) {
		acqData = 0.0;
		for (unsigned int i = 0; i < averageN_; ++i) {
			LOG_DEBUG << "ACQ " << i;
			acq_->getSignal(signalBuffer_);
			Util::addElements(signalBuffer_.begin(), acqData.begin(), acqData.end());
		}
		const FloatType factor = valueFactor_ / static_cast<FloatType>(averageN_);
		Util::multiply(acqData, factor);
	} else {
		acq_->getSignal(signalBuffer_);
		std::copy(signalBuffer_.begin(), signalBuffer_.end(), acqData.begin());
		Util::multiply(acqData, static_cast<FloatType>(valueFactor_));
	}
}

} // namespace Lab

#endif // NETWORKTNRNACQUISITION_H
