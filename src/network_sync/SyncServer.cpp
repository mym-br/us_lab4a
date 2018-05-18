/***************************************************************************
 *  Copyright 2018 Marcelo Y. Matuda                                       *
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

#include "SyncServer.h"

#include <iostream>
#include <string>

#include <boost/asio/buffer.hpp>
#include <boost/asio/read_until.hpp>
#include <boost/asio/streambuf.hpp>
#include <boost/asio/write.hpp>
#include <boost/system/error_code.hpp>

#include "Log.h"



namespace Lab {

/*******************************************************************************
 * Constructor.
 */
SyncServer::SyncServer(unsigned short portNumber)
		: ioService_()
		, acceptor_(ioService_)
		, socket_(ioService_)
{
	boost::asio::ip::tcp::endpoint endPoint(boost::asio::ip::tcp::v4(), portNumber);
	acceptor_.open(endPoint.protocol());

	acceptor_.set_option(boost::asio::ip::tcp::no_delay(true));

	boost::system::error_code ec;
	acceptor_.bind(endPoint, ec);
	if (ec) THROW_EXCEPTION(SyncBindException, "Bind error: " << ec.message() << '.');

	acceptor_.listen();
}

/*******************************************************************************
 * Destructor.
 */
SyncServer::~SyncServer()
{
}

/*******************************************************************************
 *
 */
bool
SyncServer::waitForTrigger()
{
	LOG_DEBUG << "[SyncServer] waitForTrigger()";

	if (socket_.is_open()) {
		socket_.shutdown(boost::asio::ip::tcp::socket::shutdown_both);
		socket_.close();
	}

	boost::system::error_code ec;

	// Wait for a connection.
	acceptor_.accept(socket_, ec);
	if (ec) {
		THROW_EXCEPTION(IOException, "Error in acceptor_.accept(): " << ec.message() << '.');
	}

	boost::asio::streambuf buf;
	/*size_t n =*/ boost::asio::read_until(socket_, buf, '\n', ec);
	if (ec) {
		THROW_EXCEPTION(IOException, "Error occurred while reading from client: " << ec.message() << '.');
	}

	std::istream in(&buf);
	std::string line;
	std::getline(in, line); // line will not contain the delimiter '\n'
	if (line == "TRIGGER ABORT\r") {
		LOG_DEBUG << "[SyncServer] Trigger aborted.";
		return false;
	}
	if (line != "TRIGGER\r") {
		THROW_EXCEPTION(IOException, "Invalid message: " << line << '.');
	}

	LOG_DEBUG << "[SyncServer] Trigger received.";
	return true;
}

/*******************************************************************************
 *
 */
void
SyncServer::freeTrigger()
{
	LOG_DEBUG << "[SyncServer] freeTrigger()";

	boost::system::error_code ec;

	std::string msg = "TRIGGER OK\r\n";
	boost::asio::write(socket_, boost::asio::buffer(msg), ec);
	if (ec) {
		THROW_EXCEPTION(IOException, "Error occurred while writing to the client: " << ec.message() << '.');
	}

	LOG_DEBUG << "[SyncServer] Waiting for the client to disconnect...";

	boost::asio::streambuf buf;
	size_t n = boost::asio::read_until(socket_, buf, '\n', ec);
	if (ec && ec != boost::asio::error::eof) {
		THROW_EXCEPTION(IOException, "Error occurred while reading from client: " << ec.message() << '.');
	}
	if (n != 0) {
		THROW_EXCEPTION(IOException, "Unexpected data from client.");
	}

	socket_.shutdown(boost::asio::ip::tcp::socket::shutdown_both);
	socket_.close();

	LOG_DEBUG << "[SyncServer] Connection closed by the client.";
}

} // namespace Lab
