#ifndef NETWORKSTAACQUISITION_H_
#define NETWORKSTAACQUISITION_H_

#include <cstddef> /* std::size_t */
#include <string>

#include <boost/scoped_ptr.hpp>

#include "ArrayAcqClient.h"
#include "Exception.h"
#include "global.h"
#include "Matrix2.h"
#include "ParameterMap.h"
#include "Project.h"
#include "STAAcquisition.h"
#include "STAConfiguration.h"



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
	boost::scoped_ptr<ArrayAcqClient> acq_;
	std::vector<float> ascanBuffer_;
};



template<typename FloatType>
NetworkSTAAcquisition<FloatType>::NetworkSTAAcquisition(const Project& project, const STAConfiguration<FloatType>& config)
		: project_(project)
		, config_(config)
		, acq_(0)
{
	ConstParameterMapPtr pm = project_.loadParameterMap(NETWORK_AQUISITION_CONFIG_FILE);
	std::string serverIpAddress = pm->value<std::string>(   "server_ip_address");
	unsigned short portNumber   = pm->value<unsigned short>("server_port_number", 49152, 65535);

	acq_.reset(new ArrayAcqClient(serverIpAddress.c_str(), portNumber));

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

	const std::size_t ascanLength = acq_->getAscanLength();
	if (ascanLength == 0) {
		THROW_EXCEPTION(InvalidValueException, "ascanLength = 0.");
	}
	if (acqData.n1() != config_.numElements || acqData.n2() != ascanLength) {
		acqData.resize(config_.numElements, ascanLength);
		LOG_DEBUG << "RESIZE acqData(channels,ascanLength): channels=" << config_.numElements << " ascanLength=" << ascanLength;
	}

	acq_->setBaseElement(baseElement);

	std::string txMask(config_.numElements, '0');
	txMask[txElement] = '1';
	acq_->setActiveTransmitElements(txMask);

	acq_->getAscan(ascanBuffer_);

	std::copy(ascanBuffer_.begin(), ascanBuffer_.end(), acqData.begin());
}

} // namespace Lab

#endif /* NETWORKSTAACQUISITION_H_ */
