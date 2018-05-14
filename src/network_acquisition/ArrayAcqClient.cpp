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

#include "ArrayAcqClient.h"

#include <iostream>
#include <string>



namespace Lab {

/*******************************************************************************
 * Constructor.
 */
ArrayAcqClient::ArrayAcqClient(const char* serverAddress, unsigned short portNumber)
		: ioService_()
		, socket_(ioService_)
{
	boost::asio::ip::tcp::endpoint endPoint(boost::asio::ip::address::from_string(serverAddress), portNumber);
	try {
		socket_.connect(endPoint);

		boost::asio::ip::tcp::no_delay option(true);
		socket_.set_option(option);

		connect();
	} catch (std::exception& e) {
		THROW_EXCEPTION(ArrayNetworkConnectionException, "Could not connect to the server. Cause: " << e);
	}
}

/*******************************************************************************
 * Destructor.
 */
ArrayAcqClient::~ArrayAcqClient()
{
	try {
		disconnect();
	} catch (std::exception& e) {
		std::cerr << "Error in ~ArrayAcqClient(): " << e.what() << std::endl;
	}
}

void
ArrayAcqClient::connect()
{
	prepareMessage(CONNECT_REQUEST);
	dataRawBuffer_.putUInt32(PROTOCOL_VERSION);
	sendMessage(socket_);

	handleOkOrErrorResponse();
}

void
ArrayAcqClient::disconnect()
{
	prepareMessage(DISCONNECT_REQUEST);
	sendMessage(socket_);
}

void
ArrayAcqClient::handleError(boost::uint32_t messageType)
{
	if (messageType == ERROR_RESPONSE) {
		std::string msg;
		dataRawBuffer_.getString(msg);
		THROW_EXCEPTION(ArrayNetworkServerException, "Error message from the server: " << msg);
	} else {
		THROW_EXCEPTION(InvalidArrayNetworkMessageTypeException, "Invalid message type: " << messageType);
	}
}

void
ArrayAcqClient::handleOkOrErrorResponse()
{
	boost::uint32_t messageType = receiveMessage(socket_);
	if (messageType != OK_RESPONSE) {
		handleError(messageType);
	}
}

boost::uint32_t
ArrayAcqClient::getSignalLength()
{
	prepareMessage(GET_SIGNAL_LENGTH_REQUEST);
	sendMessage(socket_);

	boost::uint32_t messageType = receiveMessage(socket_);
	if (messageType != GET_SIGNAL_LENGTH_RESPONSE) {
		handleError(messageType);
	}
	return dataRawBuffer_.getUInt32();
}

void
ArrayAcqClient::getSignal(std::vector<float>& buffer)
{
	prepareMessage(GET_SIGNAL_REQUEST);
	sendMessage(socket_);

	boost::uint32_t messageType = receiveMessage(socket_);
	if (messageType != GET_SIGNAL_RESPONSE) {
		handleError(messageType);
	}
	//dataRawBuffer_.getFloatArray(buffer);
	dataRawBuffer_.getInt16Array(buffer);
}

void
ArrayAcqClient::getSignal(float* buffer, std::size_t size)
{
	prepareMessage(GET_SIGNAL_REQUEST);
	sendMessage(socket_);

	boost::uint32_t messageType = receiveMessage(socket_);
	if (messageType != GET_SIGNAL_RESPONSE) {
		handleError(messageType);
	}
	//dataRawBuffer_.getFloatArray(buffer);
	dataRawBuffer_.getInt16Array(buffer, size);
}

float
ArrayAcqClient::getMaxSampleValue()
{
	prepareMessage(GET_MAX_SAMPLE_VALUE_REQUEST);
	sendMessage(socket_);

	boost::uint32_t messageType = receiveMessage(socket_);
	if (messageType != GET_MAX_SAMPLE_VALUE_RESPONSE) {
		handleError(messageType);
	}
	return static_cast<float>(dataRawBuffer_.getInt16());
}

float
ArrayAcqClient::getMinSampleValue()
{
	prepareMessage(GET_MIN_SAMPLE_VALUE_REQUEST);
	sendMessage(socket_);

	boost::uint32_t messageType = receiveMessage(socket_);
	if (messageType != GET_MIN_SAMPLE_VALUE_RESPONSE) {
		handleError(messageType);
	}
	return static_cast<float>(dataRawBuffer_.getInt16());
}

float
ArrayAcqClient::getSamplingFrequency()
{
	prepareMessage(GET_SAMPLING_FREQUENCY_REQUEST);
	sendMessage(socket_);

	boost::uint32_t messageType = receiveMessage(socket_);
	if (messageType != GET_SAMPLING_FREQUENCY_RESPONSE) {
		handleError(messageType);
	}
	return dataRawBuffer_.getFloat();
}

void
ArrayAcqClient::setAcquisitionTime(float acqTime)
{
	prepareMessage(SET_ACQUISITION_TIME_REQUEST);
	dataRawBuffer_.putFloat(acqTime);
	sendMessage(socket_);

	handleOkOrErrorResponse();
}

void
ArrayAcqClient::setActiveReceiveElements(const std::string& mask)
{
	prepareMessage(SET_ACTIVE_RECEIVE_ELEMENTS_REQUEST);
	dataRawBuffer_.putString(mask);
	sendMessage(socket_);

	handleOkOrErrorResponse();
}

void
ArrayAcqClient::setActiveTransmitElements(const std::string& mask)
{
	prepareMessage(SET_ACTIVE_TRANSMIT_ELEMENTS_REQUEST);
	dataRawBuffer_.putString(mask);
	sendMessage(socket_);

	handleOkOrErrorResponse();
}

void
ArrayAcqClient::setBaseElement(boost::uint32_t baseElement)
{
	prepareMessage(SET_BASE_ELEMENT_REQUEST);
	dataRawBuffer_.putUInt32(baseElement);
	sendMessage(socket_);

	handleOkOrErrorResponse();
}

void
ArrayAcqClient::setCenterFrequency(float fc)
{
	prepareMessage(SET_CENTER_FREQUENCY_REQUEST);
	dataRawBuffer_.putFloat(fc);
	sendMessage(socket_);

	handleOkOrErrorResponse();
}

void
ArrayAcqClient::setGain(float gain)
{
	prepareMessage(SET_GAIN_REQUEST);
	dataRawBuffer_.putFloat(gain);
	sendMessage(socket_);

	handleOkOrErrorResponse();
}

void
ArrayAcqClient::setReceiveDelays(const std::vector<float>& delays)
{
	prepareMessage(SET_RECEIVE_DELAYS_REQUEST);
	dataRawBuffer_.putFloatArray(delays);
	sendMessage(socket_);

	handleOkOrErrorResponse();
}

void
ArrayAcqClient::setSamplingFrequency(float fs)
{
	prepareMessage(SET_SAMPLING_FREQUENCY_REQUEST);
	dataRawBuffer_.putFloat(fs);
	sendMessage(socket_);

	handleOkOrErrorResponse();
}

void
ArrayAcqClient::setTransmitDelays(const std::vector<float>& delays)
{
	prepareMessage(SET_TRANSMIT_DELAYS_REQUEST);
	dataRawBuffer_.putFloatArray(delays);
	sendMessage(socket_);

	handleOkOrErrorResponse();
}

void
ArrayAcqClient::execPreConfiguration()
{
	prepareMessage(EXEC_PRE_CONFIGURATION_REQUEST);
	sendMessage(socket_);

	handleOkOrErrorResponse();
}

void
ArrayAcqClient::execPostConfiguration()
{
	prepareMessage(EXEC_POST_CONFIGURATION_REQUEST);
	sendMessage(socket_);

	handleOkOrErrorResponse();
}

void
ArrayAcqClient::execPreLoopConfiguration()
{
	prepareMessage(EXEC_PRE_LOOP_CONFIGURATION_REQUEST);
	sendMessage(socket_);

	handleOkOrErrorResponse();
}

void
ArrayAcqClient::execPostLoopConfiguration()
{
	prepareMessage(EXEC_POST_LOOP_CONFIGURATION_REQUEST);
	sendMessage(socket_);

	handleOkOrErrorResponse();
}

} // namespace Lab
