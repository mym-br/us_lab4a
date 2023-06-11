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

#ifndef CONTROLLER_H_
#define CONTROLLER_H_

#include <QObject>

#include "ProcessingNode.h"
#include "ProcessingThread.h"



namespace Lab {

class Project;

class Controller : public QObject {
	Q_OBJECT
public:
	explicit Controller(Project& project);
	virtual ~Controller();

	void connectProcessingNode(const ProcessingNode* node);
	void exit();
	void enableProcessing();
	bool processingEnabled() const { return (state_ == State::PROCESSING_ENABLED); }
Q_SIGNALS:
	void processingFinished();
	void processingRequested();
private Q_SLOTS:
	void execAfterProcessing();
	void execAfterError();
private:
	enum class State {
		PROCESSING_STOPPED,
		PROCESSING_ENABLED
	};

	static void showTimestamp();

	State state_; // this state can only be read/changed by the main thread
	ProcessingThread processingThread_;
};

} // namespace Lab

#endif /* CONTROLLER_H_ */
