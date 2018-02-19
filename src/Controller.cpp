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
