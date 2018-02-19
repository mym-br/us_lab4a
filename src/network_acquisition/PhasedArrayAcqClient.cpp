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

#include "PhasedArrayAcqClient.h"

#include <iostream>
#include <string>



namespace Lab {

/*******************************************************************************
 * Constructor.
 */
PhasedArrayAcqClient::PhasedArrayAcqClient(const char* serverAddress, unsigned short portNumber)
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
		THROW_EXCEPTION(PhasedArrayNetworkConnectionException, "Could not connect to the server. Cause: " << e);
	}
}

/*******************************************************************************
 * Destructor.
 */
PhasedArrayAcqClient::~PhasedArrayAcqClient()
{
	try {
		disconnect();
	} catch (std::exception& e) {
		std::cerr << "Error in ~PhasedArrayAcqClient(): " << e.what() << std::endl;
	}
}

void
PhasedArrayAcqClient::connect()
{
	prepareMessage(CONNECT_REQUEST);
	dataRawBuffer_.putUInt32(PROTOCOL_VERSION);
	sendMessage(socket_);

	handleOkOrErrorResponse();
}

void
PhasedArrayAcqClient::disconnect()
{
	prepareMessage(DISCONNECT_REQUEST);
	sendMessage(socket_);
}

void
PhasedArrayAcqClient::handleError(boost::uint32_t messageType)
{
	if (messageType == ERROR_RESPONSE) {
		std::string msg;
		dataRawBuffer_.getString(msg);
		THROW_EXCEPTION(PhasedArrayNetworkServerException, "Error message from the server: " << msg);
	} else {
		THROW_EXCEPTION(InvalidPhasedArrayNetworkMessageTypeException, "Invalid message type: " << messageType);
	}
}

void
PhasedArrayAcqClient::handleOkOrErrorResponse()
{
	boost::uint32_t messageType = receiveMessage(socket_);
	if (messageType != OK_RESPONSE) {
		handleError(messageType);
	}
}

void
PhasedArrayAcqClient::getImageBuffer(std::vector<float>& buffer)
{
	prepareMessage(GET_IMAGE_BUFFER_REQUEST);
	sendMessage(socket_);

	boost::uint32_t messageType = receiveMessage(socket_);
	if (messageType != GET_IMAGE_BUFFER_RESPONSE) {
		handleError(messageType);
	}
	dataRawBuffer_.getInt16Array(buffer);
}

boost::uint32_t
PhasedArrayAcqClient::getImageNumRows()
{
	prepareMessage(GET_IMAGE_NUM_ROWS_REQUEST);
	sendMessage(socket_);

	boost::uint32_t messageType = receiveMessage(socket_);
	if (messageType != GET_IMAGE_NUM_ROWS_RESPONSE) {
		handleError(messageType);
	}
	return dataRawBuffer_.getUInt32();
}

boost::uint32_t
PhasedArrayAcqClient::getImageNumCols()
{
	prepareMessage(GET_IMAGE_NUM_COLS_REQUEST);
	sendMessage(socket_);

	boost::uint32_t messageType = receiveMessage(socket_);
	if (messageType != GET_IMAGE_NUM_COLS_RESPONSE) {
		handleError(messageType);
	}
	return dataRawBuffer_.getUInt32();
}

void
PhasedArrayAcqClient::setGain(float value)
{
	prepareMessage(SET_GAIN_REQUEST);
	dataRawBuffer_.putFloat(value);
	sendMessage(socket_);

	handleOkOrErrorResponse();
}

void
PhasedArrayAcqClient::setAcquisitionDelay(float value)
{
	prepareMessage(SET_ACQUISITION_DELAY_REQUEST);
	dataRawBuffer_.putFloat(value);
	sendMessage(socket_);

	handleOkOrErrorResponse();
}

void
PhasedArrayAcqClient::setSignalMode(int mode)
{
	prepareMessage(SET_SIGNAL_MODE_REQUEST);
	dataRawBuffer_.putInt32(mode);
	sendMessage(socket_);

	handleOkOrErrorResponse();
}

void
PhasedArrayAcqClient::setSamplingFrequency(float fs)
{
	prepareMessage(SET_SAMPLING_FREQUENCY_REQUEST);
	dataRawBuffer_.putFloat(fs);
	sendMessage(socket_);

	handleOkOrErrorResponse();
}

float
PhasedArrayAcqClient::getSamplingFrequency()
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
PhasedArrayAcqClient::setRange(float start, float end)
{
	prepareMessage(SET_RANGE_REQUEST);
	dataRawBuffer_.putFloat(start);
	dataRawBuffer_.putFloat(end);
	sendMessage(socket_);

	handleOkOrErrorResponse();
}

void
PhasedArrayAcqClient::setSectorialScan(unsigned short baseElement, unsigned short numActiveElements, float startAngle, float endAngle, float angleStep)
{
	prepareMessage(SET_SECTORIAL_SCAN_REQUEST);
	dataRawBuffer_.putUInt16(baseElement);
	dataRawBuffer_.putUInt16(numActiveElements);
	dataRawBuffer_.putFloat(startAngle);
	dataRawBuffer_.putFloat(endAngle);
	dataRawBuffer_.putFloat(angleStep);
	sendMessage(socket_);

	handleOkOrErrorResponse();
}

void
PhasedArrayAcqClient::setFocalPoint(float emissionDistance, float receptionDistance)
{
	prepareMessage(SET_FOCAL_POINT_REQUEST);
	dataRawBuffer_.putFloat(emissionDistance);
	dataRawBuffer_.putFloat(receptionDistance);
	sendMessage(socket_);

	handleOkOrErrorResponse();
}

void
PhasedArrayAcqClient::setMaterialVelocity(float value)
{
	prepareMessage(SET_MATERIAL_VELOCITY_REQUEST);
	dataRawBuffer_.putFloat(value);
	sendMessage(socket_);

	handleOkOrErrorResponse();
}

void
PhasedArrayAcqClient::setPhasedArrayConfiguration(unsigned short numElements, float pitch, float centerFreq)
{
	prepareMessage(SET_PHASED_ARRAY_CONFIGURATION_REQUEST);
	dataRawBuffer_.putUInt16(numElements);
	dataRawBuffer_.putFloat(pitch);
	dataRawBuffer_.putFloat(centerFreq);
	sendMessage(socket_);

	handleOkOrErrorResponse();
}

void
PhasedArrayAcqClient::setApodization(const std::vector<float>& coeff)
{
	prepareMessage(SET_APODIZATION_REQUEST);
	dataRawBuffer_.putFloatArray(coeff);
	sendMessage(socket_);

	handleOkOrErrorResponse();
}

void
PhasedArrayAcqClient::resetApodization()
{
	prepareMessage(SET_APODIZATION_REQUEST);
	sendMessage(socket_);

	handleOkOrErrorResponse();
}

void
PhasedArrayAcqClient::setTimeGainCompensation(const std::vector<float>& gain, const std::vector<float>& time)
{
	prepareMessage(SET_TIME_GAIN_COMPENSATION_REQUEST);
	dataRawBuffer_.putFloatArray(gain);
	dataRawBuffer_.putFloatArray(time);
	sendMessage(socket_);

	handleOkOrErrorResponse();
}

void
PhasedArrayAcqClient::getImageLineGeometry(std::vector<float>& startX, std::vector<float>& startZ, std::vector<float>& angle)
{
	prepareMessage(GET_IMAGE_LINE_GEOMETRY_REQUEST);
	sendMessage(socket_);

	boost::uint32_t messageType = receiveMessage(socket_);
	if (messageType != GET_IMAGE_LINE_GEOMETRY_RESPONSE) {
		handleError(messageType);
	}
	dataRawBuffer_.getFloatArray(startX);
	dataRawBuffer_.getFloatArray(startZ);
	dataRawBuffer_.getFloatArray(angle);
}

void
PhasedArrayAcqClient::execPreConfiguration()
{
	prepareMessage(EXEC_PRE_CONFIGURATION_REQUEST);
	sendMessage(socket_);

	handleOkOrErrorResponse();
}

void
PhasedArrayAcqClient::execPostConfiguration()
{
	prepareMessage(EXEC_POST_CONFIGURATION_REQUEST);
	sendMessage(socket_);

	handleOkOrErrorResponse();
}

void
PhasedArrayAcqClient::execPreLoopConfiguration()
{
	prepareMessage(EXEC_PRE_LOOP_CONFIGURATION_REQUEST);
	sendMessage(socket_);

	handleOkOrErrorResponse();
}

void
PhasedArrayAcqClient::execPostLoopConfiguration()
{
	prepareMessage(EXEC_POST_LOOP_CONFIGURATION_REQUEST);
	sendMessage(socket_);

	handleOkOrErrorResponse();
}

} // namespace Lab
