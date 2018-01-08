#include "Log.h"

#include <QMutexLocker>



namespace Lab {

/*******************************************************************************
 * Static members.
 */
tbb::atomic<Log::Level> Log::level_;
std::ostringstream Log::buffer_;
QMutex Log::logMutex_;



/*******************************************************************************
 * Destructor.
 */
ErrorLog::~ErrorLog()
{
	try {
		buffer_ << '\n' << ERROR_LOG_SUFFIX;
		Log::add(buffer_);
	} catch (...) {
		// Ignore.
	}
}

/*******************************************************************************
 * Destructor.
 */
WarningLog::~WarningLog()
{
	try {
		Log::add(buffer_);
	} catch (...) {
		// Ignore.
	}
}

/*******************************************************************************
 * Destructor.
 */
InfoLog::~InfoLog()
{
	try {
		Log::add(buffer_);
	} catch (...) {
		// Ignore.
	}
}

/*******************************************************************************
 * Destructor.
 */
DebugLog::~DebugLog()
{
	try {
		Log::add(buffer_);
	} catch (...) {
		// Ignore.
	}
}

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
