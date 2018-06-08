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

#ifndef NETWORKSTAACQUISITION_H_
#define NETWORKSTAACQUISITION_H_

#include <cstddef> /* std::size_t */
#include <memory>
#include <string>

#include "ArrayAcqClient.h"
#include "Exception.h"
#include "global.h"
#include "Matrix2.h"
#include "ParameterMap.h"
#include "Project.h"
#include "STAAcquisition.h"
#include "STAConfiguration.h"
#include "Util.h"



namespace Lab {

template<typename FloatType>
class NetworkSTAAcquisition : public STAAcquisition<FloatType> {
public:
	NetworkSTAAcquisition(
		const Project& project,
		const STAConfiguration<FloatType>& config);
	virtual ~NetworkSTAAcquisition();

	virtual void execute(unsigned int baseElement, unsigned int txElement,
				typename STAAcquisition<FloatType>::AcquisitionDataType& acqData);
private:
	NetworkSTAAcquisition(const NetworkSTAAcquisition&);
	NetworkSTAAcquisition& operator=(const NetworkSTAAcquisition&);

	const Project& project_;
	const STAConfiguration<FloatType>& config_;
	std::unique_ptr<ArrayAcqClient> acq_;
	std::vector<float> signalBuffer_;
	float valueFactor_;
};



template<typename FloatType>
NetworkSTAAcquisition<FloatType>::NetworkSTAAcquisition(const Project& project, const STAConfiguration<FloatType>& config)
		: project_{project}
		, config_{config}
		, acq_{}
{
	ConstParameterMapPtr pm = project_.loadParameterMap(NETWORK_AQUISITION_CONFIG_FILE);
	const std::string serverIpAddress = pm->value<std::string>(   "server_ip_address");
	const unsigned short portNumber   = pm->value<unsigned short>("server_port_number",   49152,  65535);
	valueFactor_               = 1.0f / pm->value<float>(         "value_scale"       , 1.0e-30, 1.0e30);

	acq_ = std::make_unique<ArrayAcqClient>(serverIpAddress.c_str(), portNumber);

	acq_->execPreConfiguration();

	acq_->setSamplingFrequency(config_.samplingFrequency);
	acq_->setCenterFrequency(config_.centerFrequency);

	std::string rxMask(config_.numElements, '1');
	acq_->setActiveReceiveElements(rxMask);

	acq_->setGain(config_.minGain);
	acq_->setAcquisitionTime(config_.acquisitionTime);

	acq_->execPostConfiguration();
}

template<typename FloatType>
NetworkSTAAcquisition<FloatType>::~NetworkSTAAcquisition()
{
}

template<typename FloatType>
void
NetworkSTAAcquisition<FloatType>::execute(unsigned int baseElement, unsigned int txElement,
						typename STAAcquisition<FloatType>::AcquisitionDataType& acqData)
{
	LOG_DEBUG << "ACQ baseElement=" << baseElement << " txElement=" << txElement;

	const std::size_t signalLength = acq_->getSignalLength();
	if (signalLength == 0) {
		THROW_EXCEPTION(InvalidValueException, "signalLength = 0.");
	}
	if (acqData.n1() != config_.numElements || acqData.n2() != signalLength) {
		acqData.resize(config_.numElements, signalLength);
		LOG_DEBUG << "RESIZE acqData(channels,signalLength): channels=" << config_.numElements << " signalLength=" << signalLength;
	}

	acq_->setBaseElement(baseElement);

	std::string txMask(config_.numElements, '0');
	txMask[txElement] = '1';
	acq_->setActiveTransmitElements(txMask);

	acq_->getSignal(signalBuffer_);
	Util::multiply(signalBuffer_, valueFactor_);

	std::copy(signalBuffer_.begin(), signalBuffer_.end(), acqData.begin());
}

} // namespace Lab

#endif /* NETWORKSTAACQUISITION_H_ */
