#ifndef PHASEDARRAYACQCLIENT_H_
#define PHASEDARRAYACQCLIENT_H_

#include <vector>

#include <boost/asio/ip/tcp.hpp>
#include <boost/cstdint.hpp>

#include "PhasedArrayAcqProtocol.h"
#include "Exception.h"



namespace Lab {

struct PhasedArrayNetworkConnectionException : public virtual Exception {};

class PhasedArrayAcqClient : private PhasedArrayAcqProtocol {
public:
	PhasedArrayAcqClient(const char* serverAddress, unsigned short portNumber);
	~PhasedArrayAcqClient();

	void getImageBuffer(std::vector<float>& buffer);
	boost::uint32_t getImageNumRows();
	boost::uint32_t getImageNumCols();
	void setGain(float value /* dB */);
	void setAcquisitionDelay(float value /* 0 - 1600 µs, 25 ns step */);
	// 0: Raw / 1: Envelope
	void setSignalMode(int mode);
	void setSamplingFrequency(float fs /* Hz */);
	float getSamplingFrequency();
	void setRange(float start /* m */, float end /* m */);
	void setSectorialScan(unsigned short baseElement /* first element: 0 */, unsigned short numActiveElements,
				float startAngle /* degree */, float endAngle /* degree */, float angleStep /* degree */);
	void setFocalPoint(float emissionDistance /* m */, float receptionDistance /* m */);
	void setMaterialVelocity(float value /* m/s */);
	void setPhasedArrayConfiguration(unsigned short numElements, float pitch /* m */, float centerFreq /* Hz */);
	void setApodization(const std::vector<float>& coeff);
	void resetApodization();
	void setTimeGainCompensation(const std::vector<float>& gain /* dB */, const std::vector<float>& time /* us */);
	void getImageLineGeometry(std::vector<float>& startX /* mm */, std::vector<float>& startZ /* mm */, std::vector<float>& angle /* radian */);

	void execPreConfiguration();
	void execPostConfiguration();
	void execPreLoopConfiguration();
	void execPostLoopConfiguration();

private:
	PhasedArrayAcqClient(const PhasedArrayAcqClient&);
	PhasedArrayAcqClient& operator=(const PhasedArrayAcqClient&);

	void connect();
	void disconnect();
	void handleError(boost::uint32_t messageType);
	void handleOkOrErrorResponse();

	boost::asio::io_service ioService_;
	boost::asio::ip::tcp::socket socket_;
};

} // namespace Lab

#endif /* PHASEDARRAYACQCLIENT_H_ */
