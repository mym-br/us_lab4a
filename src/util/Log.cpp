/***************************************************************************
 *  Copyright 2014, 2017, 2018, 2019 Marcelo Y. Matuda                     *
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

#include "Log.h"

#include <iostream>

#include <QMutexLocker>



namespace Lab {

/*******************************************************************************
 * Destructor.
 */
ErrorLogMessage::~ErrorLogMessage()
{
	try {
		buffer_ << '\n' << ERROR_LOG_SUFFIX;
		Log::add(buffer_);
	} catch (...) {
		std::cerr << "Error in ~ErrorLogMessage()." << std::endl;
	}
}

/*******************************************************************************
 * Destructor.
 */
NormalLogMessage::~NormalLogMessage()
{
	try {
		Log::add(buffer_);
	} catch (...) {
		std::cerr << "Error in ~NormalLogMessage()." << std::endl;
	}
}

/*******************************************************************************
 * Static members.
 */
std::atomic<Log::Level> Log::level_;
std::ostringstream Log::buffer_;
QMutex Log::logMutex_;

/*******************************************************************************
 *
 */
void
Log::add(const std::ostringstream& inputBuffer)
{
	QMutexLocker locker(&logMutex_);
	std::streampos pos = buffer_.tellp();
	if (pos > 0) {
		buffer_ << '\n';
	}
	buffer_ << inputBuffer.str();
}

/*******************************************************************************
 *
 */
void
Log::transferTo(std::string& out)
{
	QMutexLocker locker(&logMutex_);
	out = buffer_.str();
	buffer_.str(""); // clear
}

} // namespace Lab
