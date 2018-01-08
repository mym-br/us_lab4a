#include "Controller.h"

#include <QDate>
#include <QTime>

#include "Log.h"



namespace Lab {

Controller::Controller(Project& project)
		: state_(STATE_PROCESSING_STOPPED)
		, processingThread_(*this, project)
{
	processingThread_.start();

	// Needed for queued connections and types that are not known to Qt's meta-object system.
	//qRegisterMetaType<ProcessingNode::Function>("ProcessingNode::Function");
}

Controller::~Controller()
{
	exit();
}

// This function must be called from the processing thread.
void
Controller::connectProcessingNode(const ProcessingNode* node)
{
	connect(this, SIGNAL(processingRequested()),
	        node, SLOT(process()));
	connect(node, SIGNAL(processingComplete()),
	        this, SLOT(execAfterProcessing()));
	connect(node, SIGNAL(error()),
	        this, SLOT(execAfterError()));
}

void
Controller::exit()
{
	if (processingThread_.isRunning()) {
		processingThread_.quit();
		processingThread_.wait();
	}
}

void
Controller::enableProcessing()
{
	if (state_ == STATE_PROCESSING_STOPPED) {
		state_ = STATE_PROCESSING_ENABLED;
		LOG_DEBUG << "STATE: PROCESSING_ENABLED";
		showTimestamp();
		emit processingRequested();
	}
}

void
Controller::execAfterProcessing()
{
	state_ = STATE_PROCESSING_STOPPED;
	showTimestamp();
	LOG_DEBUG << "STATE: PROCESSING_STOPPED";
	emit processingFinished();
}

void
Controller::execAfterError()
{
	state_ = STATE_PROCESSING_STOPPED;
	showTimestamp();
	LOG_DEBUG << "STATE: PROCESSING_STOPPED (after error)";
	emit processingFinished();
}

void
Controller::showTimestamp()
{
	LOG_INFO << QDate::currentDate().toString("yyyy-MM-dd").toStdString() << ' ' << QTime::currentTime().toString("hh:mm:ss").toStdString();
}

} // namespace Lab
