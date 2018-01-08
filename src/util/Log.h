#ifndef LOG_H_
#define LOG_H_

#include <cstddef> /* std::size_t */
#include <sstream>
#include <vector>

#include <QMutex>

#include <tbb/atomic.h>



#define LOG_ERROR Lab::ErrorLog()
#define LOG_WARNING if(Lab::Log::isWarningEnabled())Lab::WarningLog()
#define LOG_INFO if(Lab::Log::isInfoEnabled())Lab::InfoLog()
#define LOG_DEBUG if(Lab::Log::isDebugEnabled())Lab::DebugLog()

#define ERROR_LOG_PREFIX "ERROR >>>"
#define ERROR_LOG_SUFFIX "<<<"



namespace Lab {

/*******************************************************************************
 *
 */
class ErrorLog {
public:
	ErrorLog() { buffer_ << ERROR_LOG_PREFIX << '\n'; }
	~ErrorLog();

	template<typename T> ErrorLog& operator<<(const T& item);
private:
	std::ostringstream buffer_;
};

/*******************************************************************************
 *
 */
template<typename T>
ErrorLog&
ErrorLog::operator<<(const T& item)
{
	try {
		buffer_ << item;
	} catch (...) {
		// Ignore.
	}
	return *this;
}

/*******************************************************************************
 *
 */
class WarningLog {
public:
	WarningLog() {}
	~WarningLog();

	template<typename T> WarningLog& operator<<(const T& item);
private:
	std::ostringstream buffer_;
};

/*******************************************************************************
 *
 */
template<typename T>
WarningLog&
WarningLog::operator<<(const T& item)
{
	try {
		buffer_ << item;
	} catch (...) {
		// Ignore.
	}
	return *this;
}

/*******************************************************************************
 *
 */
class InfoLog {
public:
	InfoLog() {}
	~InfoLog();

	template<typename T> InfoLog& operator<<(const T& item);
private:
	std::ostringstream buffer_;
};

/*******************************************************************************
 *
 */
template<typename T>
InfoLog&
InfoLog::operator<<(const T& item)
{
	try {
		buffer_ << item;
	} catch (...) {
		// Ignore.
	}
	return *this;
}

/*******************************************************************************
 *
 */
class DebugLog {
public:
	DebugLog() {}
	~DebugLog();

	template<typename T> DebugLog& operator<<(const T& item);
private:
	std::ostringstream buffer_;
};

/*******************************************************************************
 *
 */
template<typename T>
DebugLog&
DebugLog::operator<<(const T& item)
{
	try {
		buffer_ << item;
	} catch (...) {
		// Ignore.
	}
	return *this;
}

/*******************************************************************************
 *
 */
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
	static tbb::atomic<Level> level_;
	static std::ostringstream buffer_;
	static QMutex logMutex_;

	Log() {}
	Log(const Log&);
	Log& operator=(const Log&);
};

template<typename T>
std::ostringstream&
operator<<(std::ostringstream& out, const std::vector<T>& v)
{
	out << '(';
	if (v.size() > 0) {
		out << v[0];
	}
	for (std::size_t i = 1; i < v.size(); ++i) {
		out << ", " << v[i];
	}
	out << ')';
	return out;
}

} // namespace Lab

#endif /* LOG_H_ */
