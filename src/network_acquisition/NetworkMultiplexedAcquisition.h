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

#ifndef NETWORKMULTIPLEXEDACQUISITION_H_
#define NETWORKMULTIPLEXEDACQUISITION_H_

#include <algorithm> /* std::copy */
#include <cstddef> /* std::size_t */
#include <string>
#include <vector>

#include <boost/scoped_ptr.hpp>

#include "ArrayAcqClient.h"
#include "Exception.h"
#include "ParameterMap.h"
#include "Project.h"
#include "MultiplexedAcquisition.h"
#include "global.h"
#include "Log.h"



namespace Lab {

// All the elements in the group emit and receive in each acquisition.
// Focalization only in emission.
template<typename FloatType, typename Config>
class NetworkMultiplexedAcquisition : public MultiplexedAcquisition<FloatType> {
public:
	NetworkMultiplexedAcquisition(
		const Project& project,
		const Config& config);
	virtual ~NetworkMultiplexedAcquisition();

	virtual void execute(
			unsigned int baseElement,
			const std::vector<FloatType>& txDelays,
			typename MultiplexedAcquisition<FloatType>::AcquisitionDataType& acqData);
private:
	NetworkMultiplexedAcquisition(const NetworkMultiplexedAcquisition&);
	NetworkMultiplexedAcquisition& operator=(const NetworkMultiplexedAcquisition&);

	void setTxDelays(const std::vector<float>& txDelays);
	void setTxDelays(const std::vector<double>& txDelays);
	void getSignal(MultiplexedAcquisition<float>::AcquisitionDataType& acqData);
	void getSignal(MultiplexedAcquisition<double>::AcquisitionDataType& acqData);

	const Project& project_;
	const Config& config_;
	boost::scoped_ptr<ArrayAcqClient> acq_;
	std::vector<float> txDelays_;
	std::vector<float> signalBuffer_;
};



template<typename FloatType, typename Config>
NetworkMultiplexedAcquisition<FloatType, Config>::NetworkMultiplexedAcquisition(const Project& project, const Config& config)
		: project_(project)
		, config_(config)
		, txDelays_(config_.numElements)
{
	ConstParameterMapPtr pm = project_.loadParameterMap(NETWORK_AQUISITION_CONFIG_FILE);
	std::string serverIpAddress = pm->value<std::string>(   "server_ip_address");
	unsigned short portNumber   = pm->value<unsigned short>("server_port_number", 49152, 65535);

	acq_.reset(new ArrayAcqClient(serverIpAddress.c_str(), portNumber));

	acq_->execPreConfiguration();

	acq_->setSamplingFrequency(config_.samplingFrequency);
	acq_->setCenterFrequency(config_.centerFrequency);

	std::string allEnabledMask(config_.numElements, '1');
	acq_->setActiveReceiveElements(allEnabledMask);
	acq_->setActiveTransmitElements(allEnabledMask);

	acq_->setGain(config_.minGain);
	acq_->setAcquisitionTime(config_.acquisitionTime);

	// Doesn't use receive delays.
	std::vector<float> rxDelays(config_.numElements, 0.0);
	acq_->setReceiveDelays(rxDelays);

	acq_->execPostConfiguration();
}

template<typename FloatType, typename Config>
NetworkMultiplexedAcquisition<FloatType, Config>::~NetworkMultiplexedAcquisition()
{
}

template<typename FloatType, typename Config>
void
NetworkMultiplexedAcquisition<FloatType, Config>::setTxDelays(const std::vector<float>& txDelays)
{
	acq_->setTransmitDelays(txDelays);
}

template<typename FloatType, typename Config>
void
NetworkMultiplexedAcquisition<FloatType, Config>::setTxDelays(const std::vector<double>& txDelays)
{
	std::copy(txDelays.begin(), txDelays.end(), txDelays_.begin()); // double --> float
	acq_->setTransmitDelays(txDelays_);
}

template<typename FloatType, typename Config>
void
NetworkMultiplexedAcquisition<FloatType, Config>::getSignal(MultiplexedAcquisition<float>::AcquisitionDataType& acqData)
{
	acq_->getSignal(&acqData(0, 0), acqData.size());
}

template<typename FloatType, typename Config>
void
NetworkMultiplexedAcquisition<FloatType, Config>::getSignal(MultiplexedAcquisition<double>::AcquisitionDataType& acqData)
{
	acq_->getSignal(signalBuffer_);
	std::copy(signalBuffer_.begin(), signalBuffer_.end(), acqData.begin()); // float --> double
}

template<typename FloatType, typename Config>
void
NetworkMultiplexedAcquisition<FloatType, Config>::execute(
		unsigned int baseElement,
		const std::vector<FloatType>& txDelays,
		typename MultiplexedAcquisition<FloatType>::AcquisitionDataType& acqData)
{
	LOG_DEBUG << "ACQ baseElement=" << baseElement;

	//TODO: check
//	if (txDelays.size() != config_.numElements) {
//		THROW_EXCEPTION(InvalidParameterException, "Invalid txDelays size: " << txDelays.size() << '.');
//	}

	const std::size_t signalLength = acq_->getSignalLength();
	if (signalLength == 0) {
		THROW_EXCEPTION(InvalidParameterException, "signalLength = 0.");
	}
	if (acqData.n1() != config_.numElements || acqData.n2() != signalLength) {
		acqData.resize(config_.numElements, signalLength);
		LOG_DEBUG << "RESIZE acqData(channels=" << config_.numElements << ", signalLength=" << signalLength << ')';
	}

	acq_->setBaseElement(baseElement);

	setTxDelays(txDelays);

	getSignal(acqData);
}

} // namespace Lab

#endif /* NETWORKMULTIPLEXEDACQUISITION_H_ */
