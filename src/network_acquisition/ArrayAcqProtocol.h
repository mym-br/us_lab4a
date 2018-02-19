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

#ifndef ARRAYACQPROTOCOL_H_
#define ARRAYACQPROTOCOL_H_

#include "RawBuffer.h"

#include <boost/asio/read.hpp>
#include <boost/asio/write.hpp>
#include <boost/cstdint.hpp>

#include "Exception.h"



namespace Lab {

struct InvalidArrayNetworkMessageTypeException : public virtual Exception {};
struct InvalidArrayNetworkRequestException : public virtual Exception {};
struct ArrayNetworkServerException : public virtual Exception {};

/*******************************************************************************
 *
 */
class ArrayAcqProtocol {
protected:
	enum {
		HEADER_RAW_BUFFER_SIZE = 8,
		PROTOCOL_VERSION = 1006
	};
	enum MessageType {
		CONNECT_REQUEST = 2001,
		DISCONNECT_REQUEST,

		OK_RESPONSE,
		ERROR_RESPONSE,

		GET_ASCAN_LENGTH_REQUEST,
		GET_ASCAN_LENGTH_RESPONSE,
		GET_ASCAN_REQUEST,
		GET_ASCAN_RESPONSE,
		GET_MAX_SAMPLE_VALUE_REQUEST,
		GET_MAX_SAMPLE_VALUE_RESPONSE,
		GET_MIN_SAMPLE_VALUE_REQUEST,
		GET_MIN_SAMPLE_VALUE_RESPONSE,
		GET_SAMPLING_FREQUENCY_REQUEST,
		GET_SAMPLING_FREQUENCY_RESPONSE,

		SET_ACQUISITION_TIME_REQUEST,
		SET_ACTIVE_RECEIVE_ELEMENTS_REQUEST,
		SET_ACTIVE_TRANSMIT_ELEMENTS_REQUEST,
		SET_BASE_ELEMENT_REQUEST,
		SET_CENTER_FREQUENCY_REQUEST,
		SET_GAIN_REQUEST,
		SET_RECEIVE_DELAYS_REQUEST,
		SET_SAMPLING_FREQUENCY_REQUEST,
		SET_TRANSMIT_DELAYS_REQUEST,

		EXEC_PRE_CONFIGURATION_REQUEST,
		EXEC_POST_CONFIGURATION_REQUEST,
		EXEC_PRE_LOOP_CONFIGURATION_REQUEST,
		EXEC_POST_LOOP_CONFIGURATION_REQUEST
	};

	ArrayAcqProtocol() {}
	~ArrayAcqProtocol() {}

	void prepareMessage(MessageType type);
	void sendMessage(boost::asio::ip::tcp::socket& socket);
	boost::uint32_t receiveMessage(boost::asio::ip::tcp::socket& socket);

	RawBuffer headerRawBuffer_;
	RawBuffer dataRawBuffer_;
private:
	ArrayAcqProtocol(const ArrayAcqProtocol&);
	ArrayAcqProtocol& operator=(const ArrayAcqProtocol&);
};

/*******************************************************************************
 *
 */
inline
void
ArrayAcqProtocol::prepareMessage(MessageType type)
{
	headerRawBuffer_.reset();
	headerRawBuffer_.putUInt32(type);

	dataRawBuffer_.reset();
}

/*******************************************************************************
 *
 */
inline
void
ArrayAcqProtocol::sendMessage(boost::asio::ip::tcp::socket& socket)
{
	headerRawBuffer_.putUInt32(dataRawBuffer_.size());

	std::size_t n = boost::asio::write(socket, boost::asio::buffer(&(headerRawBuffer_.front()), headerRawBuffer_.size()));
	if (n != headerRawBuffer_.size()) {
		THROW_EXCEPTION(IOException, "Wrong number of header bytes sent: " << n << " (expected: " << headerRawBuffer_.size() << ").");
	}
	if (dataRawBuffer_.size() > 0) {
		n = boost::asio::write(socket, boost::asio::buffer(&dataRawBuffer_.front(), dataRawBuffer_.size()));
		if (n != dataRawBuffer_.size()) {
			THROW_EXCEPTION(IOException, "Wrong number of data bytes sent: " << n << " (expected: " << dataRawBuffer_.size() << ").");
		}
	}
}

/*******************************************************************************
 *
 */
inline
boost::uint32_t
ArrayAcqProtocol::receiveMessage(boost::asio::ip::tcp::socket& socket)
{
	headerRawBuffer_.reserve(HEADER_RAW_BUFFER_SIZE);
	std::size_t n = boost::asio::read(socket, boost::asio::buffer(&headerRawBuffer_.front(), headerRawBuffer_.size()));
	if (n != headerRawBuffer_.size()) {
		THROW_EXCEPTION(IOException, "Wrong number of header bytes received: " << n << " (expected: " << headerRawBuffer_.size() << ").");
	}

	boost::uint32_t messageType = headerRawBuffer_.getUInt32();
	const boost::uint32_t dataSize = headerRawBuffer_.getUInt32();
	if (dataSize > 0) {
		dataRawBuffer_.reserve(dataSize);
		n = boost::asio::read(socket, boost::asio::buffer(&dataRawBuffer_.front(), dataRawBuffer_.size()));
		if (n != dataRawBuffer_.size()) {
			THROW_EXCEPTION(IOException, "Wrong number of data bytes received: " << n << " (expected: " << dataRawBuffer_.size() << ").");
		}
	}

	return messageType;
}

} // namespace Lab

#endif /* ARRAYACQPROTOCOL_H_ */
