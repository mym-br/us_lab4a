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

	void getSignal(std::vector<float>& buffer);
	void getSignal(float* buffer, std::size_t size);
	boost::uint32_t getSignalLength();
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
