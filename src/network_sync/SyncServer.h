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

#ifndef SYNCSERVER_H_
#define SYNCSERVER_H_

#include <boost/asio/io_service.hpp>
#include <boost/asio/ip/tcp.hpp>

#include "Exception.h"



namespace Lab {

struct SyncBindException : std::runtime_error {
	using std::runtime_error::runtime_error;
};

/*******************************************************************************
 *
 */
class SyncServer {
public:
	SyncServer(unsigned short portNumber);
	~SyncServer();

	bool waitForTrigger();
	void freeTrigger();
private:
	SyncServer(const SyncServer&) = delete;
	SyncServer& operator=(const SyncServer&) = delete;

	boost::asio::io_service ioService_;
	boost::asio::ip::tcp::acceptor acceptor_;
	boost::asio::ip::tcp::socket socket_;
};

} // namespace Lab

#endif /* SYNCSERVER_H_ */
