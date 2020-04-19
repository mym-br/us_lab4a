/***************************************************************************
 *  Copyright 2014, 2017, 2018, 2019, 2020 Marcelo Y. Matuda               *
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

#ifndef LOG_H_
#define LOG_H_

#include <atomic>
#include <cstddef> /* std::size_t */
#include <mutex>
#include <sstream>
#include <vector>

#include "Stream.h"

#define LOG_ERROR Lab::ErrorLogMessage()
#define LOG_WARNING if(Lab::Log::isWarningEnabled())Lab::NormalLogMessage()
#define LOG_INFO if(Lab::Log::isInfoEnabled())Lab::NormalLogMessage()
#define LOG_DEBUG if(Lab::Log::isDebugEnabled())Lab::NormalLogMessage()

constexpr const char* ERROR_LOG_PREFIX = "ERROR >>>";
constexpr const char* ERROR_LOG_SUFFIX = "<<<";

namespace Lab {

class ErrorLogMessage {
public:
	ErrorLogMessage() { buffer_ << ERROR_LOG_PREFIX << '\n'; }
	~ErrorLogMessage();
	template<typename T> ErrorLogMessage& operator<<(const T& item) {
		buffer_ << item;
		return *this;
	}
private:
	std::ostringstream buffer_;
};

class NormalLogMessage {
public:
	~NormalLogMessage();
	template<typename T> NormalLogMessage& operator<<(const T& item) {
		buffer_ << item;
		return *this;
	}
private:
	std::ostringstream buffer_;
};

class Log {
public:
	enum Level {
		LEVEL_ERROR,
		LEVEL_WARNING,
		LEVEL_INFO,
		LEVEL_DEBUG
	};

	static bool isWarningEnabled() { return level_ >= LEVEL_WARNING; }
	static bool isInfoEnabled() { return level_ >= LEVEL_INFO; }
	static bool isDebugEnabled() { return level_ >= LEVEL_DEBUG; }
	static void setLevel(Level level) { level_ = level; }
	static void add(const std::ostringstream& inputBuffer);
	static void transferTo(std::string& out);
private:
	Log() = delete;

	static std::atomic<Level> level_;
	static std::ostringstream buffer_;
	static std::mutex logMutex_;
};

} // namespace Lab

#endif /* LOG_H_ */
