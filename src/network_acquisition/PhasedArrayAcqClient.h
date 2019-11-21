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

#ifndef PHASEDARRAYACQCLIENT_H_
#define PHASEDARRAYACQCLIENT_H_

#include <vector>

#include <boost/asio/io_service.hpp>
#include <boost/asio/ip/tcp.hpp>
#include <boost/cstdint.hpp>

#include "PhasedArrayAcqProtocol.h"
#include "Exception.h"



namespace Lab {

struct PhasedArrayNetworkConnectionException : std::runtime_error {
	using std::runtime_error::runtime_error;
};

class PhasedArrayAcqClient : private PhasedArrayAcqProtocol {
public:
	PhasedArrayAcqClient(const char* serverAddress, unsigned short portNumber);
	~PhasedArrayAcqClient();

	void getImageBuffer(std::vector<float>& buffer);
	boost::uint32_t getImageNumRows();
	boost::uint32_t getImageNumCols();
	void setGain(float value /* dB */);
	void setAcquisitionDelay(float value /* 0 - 1600 us, 25 ns step */);
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
	PhasedArrayAcqClient(const PhasedArrayAcqClient&) = delete;
	PhasedArrayAcqClient& operator=(const PhasedArrayAcqClient&) = delete;
	PhasedArrayAcqClient(PhasedArrayAcqClient&&) = delete;
	PhasedArrayAcqClient& operator=(PhasedArrayAcqClient&&) = delete;

	void connect();
	void disconnect();
	void handleError(boost::uint32_t messageType);
	void handleOkOrErrorResponse();

	boost::asio::io_service ioService_;
	boost::asio::ip::tcp::socket socket_;
};

} // namespace Lab

#endif /* PHASEDARRAYACQCLIENT_H_ */
