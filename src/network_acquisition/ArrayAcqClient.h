#ifndef ARRAYACQCLIENT_H_
#define ARRAYACQCLIENT_H_

#include <cstddef> /* std::size_t */
#include <vector>

#include <boost/asio/ip/tcp.hpp>
#include <boost/cstdint.hpp>

#include "ArrayAcqProtocol.h"
#include "Exception.h"



namespace Lab {

struct ArrayNetworkConnectionException : public virtual Exception {};

class ArrayAcqClient : private ArrayAcqProtocol {
public:
	ArrayAcqClient(const char* serverAddress, unsigned short portNumber);
	~ArrayAcqClient();

	void getAscan(std::vector<float>& buffer);
	void getAscan(float* buffer, std::size_t size);
	boost::uint32_t getAscanLength();
	float getMaxSampleValue();
	float getMinSampleValue();
	float getSamplingFrequency();

	void setAcquisitionTime(float acqTime /* s */);
	void setActiveReceiveElements(const std::string& mask /* '0', '1'*/);
	void setActiveTransmitElements(const std::string& mask /* '0', '1'*/);
	void setBaseElement(boost::uint32_t baseElement /* 0 ... (NUM_CHANNELS_MUX - NUM_CHANNELS)*/);
	void setCenterFrequency(float fc /* MHz */);
	void setGain(float gain /* device-dependent */);
	void setReceiveDelays(const std::vector<float>& delays /* s */);
	void setSamplingFrequency(float fs /* MHz */);
	void setTransmitDelays(const std::vector<float>& delays /* s */);

	void execPreConfiguration();
	void execPostConfiguration();
	void execPreLoopConfiguration();
	void execPostLoopConfiguration();
private:
	ArrayAcqClient(const ArrayAcqClient&);
	ArrayAcqClient& operator=(const ArrayAcqClient&);

	void connect();
	void disconnect();
	void handleError(boost::uint32_t messageType);
	void handleOkOrErrorResponse();

	boost::asio::io_service ioService_;
	boost::asio::ip::tcp::socket socket_;
};

} // namespace Lab

#endif /* ARRAYACQCLIENT_H_ */
